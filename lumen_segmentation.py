import numpy as np
import pandas as pd
import os
import glob
import sys
import matplotlib.pyplot as plt
import skimage
from typing import Union
import alphashape
import shapely
from pathlib import Path
from mpl_toolkits.axes_grid1 import ImageGrid
from zipfile import ZipFile
from scipy import spatial

USE_GPU = False
# These are imported in this way so that some operations can be forced to CPU even when GPU is used for compute
if USE_GPU:
    from cucim import skimage as gpu_skimage
    import cupy as cp

    sk: gpu_skimage = gpu_skimage
    xp: cp = cp
else:
    sk: skimage = skimage
    xp: np = np

# Unzip Dataset Function
def unzip_dataset(path="dataset.zip"):
    """
    Unzip dataset into project directory if no dataset exists
    :param path: str, path to dataset directory
    """
    print("[+] Unzipping dataset")

    unzipped_dataset_path = "dataset"

    if not os.path.exists(unzipped_dataset_path):
        os.mkdir(unzipped_dataset_path)

    if len(os.listdir(unzipped_dataset_path)) == 0:
        with ZipFile(path, "r") as dataset:
            dataset.extractall(path=unzipped_dataset_path)

def prepare_df_for_annotations(csv_spatial_data_path, arbitrary_scale=2):
    """
    Reads csv contents into a data frame, scales centroids, and casts values to
    int.
    :param csv_spatial_data_path: str, path to csv file
    :param arbitrary_scale: int, scale factor for centroids
    :returns df: DataFrame, data frame with scaled int centroids
    """
    df = pd.read_csv(csv_spatial_data_path)

    df["x_scaled_centroid"] = (df["x_centroid"].values * arbitrary_scale).astype(int)
    df["y_scaled_centroid"] = (df["y_centroid"].values * arbitrary_scale).astype(int)

    return df

def prepare_df_for_loading_labels(csv_spatial_data_path, arbitrary_scale=2):
    """
    Reads csv contents into a data frame, scales mRNA locations, and casts values to
    int.
    :param csv_spatial_data_path: str, path to csv file
    :param arbitrary_scale: int, scale factor for centroids
    :returns df: DataFrame, data frame with scaled int centroids
    """
    df = pd.read_csv(csv_spatial_data_path)

    df["x_loc_int"] = (df["x_location"].values * arbitrary_scale).astype(int)
    df["y_loc_int"] = (df["y_location"].values * arbitrary_scale).astype(int)

    return df


# Transcript To Image Functions
def load_spatial_data_clust_label(filename: str, label_exclude: list[int, ...] = None, cache: bool = True, recreate: bool = False, genes_include: list[str, ...] = None) -> np.ndarray:
    """
    Load a spatial transcriptome file and return a 2D array
    :param filename: str, input file to be rad with pd.read_csv
    :param label_exclude: list, which clusters/labels to exclude
    :param cache: bool, store a .npy of the image for faster loading later
    :param recreate: bool, overwrite and recreate the cached file
    :return: xp.xdarray (2D), labelmap of transcripts where value signifies cell
    cluster number
    """
    cache_fname = os.path.splitext(filename)[0] + ".npy"

    if (cache and not os.path.exists(cache_fname)) or recreate:
        df = pd.read_csv(filename)

        arbitrary_scale: int = 2  # Adjust this to scale the image to get additional 'accuracy' on localization

        df["x_loc_int"] = (df["x_location"].values * arbitrary_scale).astype(int)
        df["y_loc_int"] = (df["y_location"].values * arbitrary_scale).astype(int)

        # Create a new empty image that is the right size to hold the data
        img_arr = np.zeros((np.max(df["x_loc_int"]) + 1, np.max(df["y_loc_int"]) + 1)).astype(int)

        if genes_include is not None:
            print(f"Only keeping specified genes")
            df = df[df['feature_name'].isin(genes_include)]

        for label in np.unique(df.loc[:, "gmm10_xenium_fullpanel_30ktrained_3NB"]):
            if label_exclude is not None and label in label_exclude:
                print(f"[+] Skipping {label}")
                continue

            # Vectorized image load operation, specify X and Y locations in the slice operation
            df_match_label = df.loc[df.loc[:, "gmm10_xenium_fullpanel_30ktrained_3NB"] == label, ]
            img_arr[df_match_label.loc[:,"x_loc_int"], df_match_label.loc[:,"y_loc_int"]] = label

        np.save(cache_fname, img_arr)
    else:
        print("Loading cached file")
        img_arr = np.load(cache_fname)


    return(img_arr)

def load_spatial_data_clust_layer(filename: str, label_exclude: list[int, ...] = None, cache: bool = True, recreate: bool = False) -> np.ndarray:
    """
    Load a spatial transcriptome file and return a 3D array, first dimension is
    the transcriptional cluster
    :param filename: str, input file to be rad with pd.read_csv
    :param label_exclude: list, which clusters/labels to exclude
    :param cache: bool, store a .npy of the image for faster loading later
    :param recreate: bool, overwrite and recreate the cached file
    :return: xp.xdarray (3D), dimensions are the transcriptional cluster, the x
    location, and the y location respectively
    """
    cache_fname = os.path.splitext(filename)[0] + "_layers.npy"

    if (cache and not os.path.exists(cache_fname)) or recreate:
        df = pd.read_csv(filename)

        arbitrary_scale: int = 2  # Adjust this to scale the image to get additional 'accuracy' on localization

        df["x_loc_int"] = (df["x_location"].values * arbitrary_scale).astype(int)
        df["y_loc_int"] = (df["y_location"].values * arbitrary_scale).astype(int)

        # Create a new empty image that is the right size to hold the data
        img_arr = np.zeros((np.max(df["gmm10_xenium_fullpanel_30ktrained_3NB"]), np.max(df["x_loc_int"]) + 1, np.max(df["y_loc_int"]) + 1)).astype(int)

        for label in np.unique(df.loc[:, "gmm10_xenium_fullpanel_30ktrained_3NB"]):
            if label_exclude is not None and label in label_exclude:
                print(f"Skipping {label}")
                continue

            # Vectorized image load operation, specify X and Y locations in the slice operation
            df_match_label = df.loc[df.loc[:, "gmm10_xenium_fullpanel_30ktrained_3NB"] == label, ]
            img_arr[label, df_match_label.loc[:, "x_loc_int"], df_match_label.loc[:,"y_loc_int"]] = 1

        np.save(cache_fname, img_arr)
    else:
        print("[+] Loading cached file")

        img_arr = np.load(cache_fname)

    return(img_arr)


# Segmentation Helper Functions
def mask_with_alpha_shapes(img: np.ndarray) -> np.ndarray:
    """
    This function efficiently identifies foreground from background using
    alphashapes from a binary image where 1 = Foreground. Operation downscales
    image 10x for processing, and upscales 10x at the end
    :param img: xp.xdarray (2D), binary labelmap
    :return: xp.xdarray (2D), binary labelmap with alphashape mask applied
    """
    print("[+] Masking points with alpha shapes")

    pad_size = 500
    # Pad so that binary_closing doesn't cause weird things on the edges
    img_pad = xp.pad(img, ((pad_size, pad_size),(pad_size, pad_size)))
    img_filt = sk.morphology.remove_small_objects(img_pad, 5000)
    img_filt_closed = sk.morphology.binary_closing(img_filt, sk.morphology.disk(100))
    img_filt_closed_cleaned = sk.morphology.remove_small_objects(img_filt_closed, 50000)
    img_labeled = sk.measure.label(img_filt_closed_cleaned)
    img_downscale_10x = sk.transform.resize(img_labeled, [i // 10 for i in img_labeled.shape],
                                            order=0, preserve_range=True, anti_aliasing=False).astype(int)

    del img_pad
    del img_filt
    del img_filt_closed
    del img_filt_closed_cleaned
    # Send to CPU for alphashape
    if USE_GPU: img_downscale_10x = cp.asnumpy(img_downscale_10x)

    # Split each connected, labeled, object into a different matrix
    label_list = []
    img_list = []
    for label in range(1, np.max(img_downscale_10x) + 1):
        individual_obj = img_downscale_10x.copy()
        individual_obj[individual_obj != label] = 0
        img_list.append(individual_obj)
        # mask is a 2D boolean np.array containing the segmentation result
        contour_img = np.logical_xor(individual_obj, skimage.morphology.binary_erosion(individual_obj))
        contour_y, contour_x = np.where(contour_img)
        label_list.append((contour_y, contour_x))

    img_out = np.zeros(img_downscale_10x.shape, dtype=bool)

    for label_idx, labelmask in enumerate(label_list):
        points = np.array([labelmask[0], labelmask[1]]).T

        # Use an alphashape and decrease its sensitivity until a single polygon is found
        # At worst (alpha = 0), this uses the convex hull of the connected object
        # At best, this is a close approximation of the shape
        iter_step = 0.01
        alpha = 0.05
        alpha_shape = alphashape.alphashape(points, alpha)
        in_img = img_list[label_idx] > 0

        while type(alpha_shape) != shapely.geometry.polygon.Polygon:
            alpha = alpha - iter_step

            if alpha == 0: print("No ideal solution found, using convex hull")

            alpha_shape = alphashape.alphashape(points, alpha)

        # Create a mask from the polygon, not implemented for GPU
        # Create a mask from the polygon, not implemented for GPU
        raster_mask = skimage.draw.polygon2mask(img_downscale_10x.shape, np.asarray(list(zip(*alpha_shape.exterior.xy))))
        raster_mask = np.bitwise_or(raster_mask, in_img)  # Make sure we don't drop anything
        img_out = np.bitwise_or(img_out, raster_mask) # Assemble all objects into the final mask

    # Resize output
    if USE_GPU: img_out = cp.asarray(img_out) # Return to GPU for downstream if needed
    scaleup = sk.transform.resize(img_out, img_labeled.shape,
                                  order=0, preserve_range=True, anti_aliasing=False).astype(bool)

    scaleup = scaleup[pad_size:-pad_size, pad_size:-pad_size]

    return scaleup

def remove_too_large(labeled, size_threshold=2000000):
    """
    Removes airspaces that are larger than a certain number of pixels
    :param labeled: xp.xdarray (2D), segmented image labelmap
    :return: xp.xdarray (2D), labelmap with large airspace removed
    """
    print("[+] Removing large blobs")

    component_sizes = xp.bincount(labeled.ravel())
    too_big = component_sizes > size_threshold
    too_small_mask = too_big[labeled]
    labeled[too_small_mask] = 0

    return labeled


def segment_airways(spatial_data):
    """
    Segments airways seperately
    :param spatial_data: xp.xdarray (2D), spatial data for a singular image
    :return: xp.xdarray (2D), mask of airway tissue and lumen
    """
    print("[+] Segmenting airways")

    if USE_GPU: spatial_data = cp.asarray(spatial_data) # Send data to GPU

    only_airway = xp.where(spatial_data == 1, 1, 0)
    dialated = sk.morphology.binary_dilation(only_airway, sk.morphology.disk(10))
    small_objects_removed = sk.morphology.remove_small_objects(dialated, 60)
    masked_airway = mask_with_alpha_shapes(small_objects_removed)
    dilated_again = sk.morphology.binary_dilation(masked_airway, sk.morphology.disk(10))

    inverted = xp.invert(dilated_again)
    labeled = sk.measure.label(inverted)

    labeled = xp.where(labeled > 1, 1, 0)
    masked_airway = xp.bitwise_or(labeled, masked_airway)

    # Erosion is 15 not 20 to compensate for the erosion step in the main pipline, will change to variable
    masked_airway = sk.morphology.binary_erosion(masked_airway, sk.morphology.disk(15))

    inverted_only_airway = xp.invert(dialated)

    lumen = xp.bitwise_and(inverted_only_airway, masked_airway)

    return lumen


def expand_fibroblasts(spatial_data, dialation):
    """
    Expand signals from fibroblasts (clusters 5 and 6)
    :param spatial_data: xp.xdarray (2D), spatial data for a singular image
    :param dialation: sk.morphology footprint for dialation shape and size
    :return: xp.xdarray (2D), image with expanded fibroblasts alone
    """
    print("[+] Expanding fibroblasts")

    if USE_GPU: spatial_data = cp.asarray(spatial_data) # Send data to GPU

    only_fibroblasts = xp.where((spatial_data == 5) | (spatial_data == 6), 1, 0)

    expanded_fibroblasts = sk.morphology.binary_dilation(only_fibroblasts, dialation)

    return expanded_fibroblasts

# Segmention Pipeline
def run_segmentation(spatial_data: np.ndarray) -> list[np.ndarray, ...]:
    """
    Run segmentation on spatial data
    :param spatial_data: xp.xdarray, spatial data for a singular image
    :return: list, list of labelmaps (xp.xdarray) from each step in segmentation
    pipeline
    """
    print("[+] Running segmentation")

    initial = spatial_data.astype(int)

    expanded_fibroblasts = expand_fibroblasts(initial, sk.morphology.disk(10))
    binarized = initial > 0

    if USE_GPU: binarized = cp.asarray(binarized)

    expanded_fibroblasts_binary = xp.bitwise_or(expanded_fibroblasts, binarized)
    closed_one = sk.morphology.binary_closing(expanded_fibroblasts_binary, sk.morphology.disk(4))
    dialated_one = sk.morphology.binary_dilation(closed_one, sk.morphology.square(5))
    small_objects_removed = sk.morphology.remove_small_objects(dialated_one, 60)
    dilated_two = sk.morphology.binary_dilation(small_objects_removed, sk.morphology.square(5))
    closed_two = sk.morphology.binary_closing(dilated_two, sk.morphology.disk(5))

    airway_lumens = segment_airways(initial)
    inverse_airway_lumens = xp.invert(airway_lumens)
    closed_with_segmented_airways = xp.bitwise_and(closed_two, inverse_airway_lumens)

    foreground_cleaned = sk.morphology.remove_small_objects(closed_with_segmented_airways, 200)
    background_cleaned = sk.morphology.remove_small_holes(foreground_cleaned, 200)
    eroded = sk.morphology.binary_erosion(background_cleaned, sk.morphology.disk(5))
    small_spaces_removed = sk.morphology.remove_small_holes(background_cleaned, 400)  # Set lower bound on size of airspace

    mask = mask_with_alpha_shapes(small_spaces_removed)
    inverse_mask = xp.invert(mask)
    dialated_inverse_mask = sk.morphology.binary_dilation(inverse_mask, sk.morphology.disk(15))
    background_masked = xp.bitwise_or(small_spaces_removed, dialated_inverse_mask)

    inverted = xp.invert(background_masked)
    labeled = sk.measure.label(inverted)

    too_large_removed = labeled

    if USE_GPU: too_large_removed = cp.asnumpy(too_large_removed)

    colored = skimage.color.label2rgb(too_large_removed, colors = colors)  # Don't run this on GPU, needs too much memory


    if USE_GPU:
        return [initial, background_masked.get(), inverted.get(),
                labeled.get(), colored]
    else:
        return [initial, binarized, expanded_fibroblasts_binary,
                closed_one, dialated_one, small_objects_removed,
                dilated_two, closed_two, closed_with_segmented_airways, foreground_cleaned,
                background_cleaned, eroded, background_masked, inverted,
                labeled, colored]

#Display Settings
# Create a list of more colors to apply to labels
cm = plt.get_cmap("tab20b")
colors = [cm(i) for i in range(cm.N)]
cm2 = plt.get_cmap("tab20c")
colors2 = [cm2(i) for i in range(cm2.N)]
colors.extend(colors2)

plt.rcParams["figure.figsize"] = [8, 8]
# Display Functions
def display_all(subplots, filename = None):
    """
    Display all steps of processing in a table (designed to assist parameter tuning)
    :param subplots: list, list of each output of each step in the segmentation process
    :param filename: str, name of file to save image to
    """
    print("[+] Displaying all steps in segmentation")

    figure = plt.figure(figsize=(12, 12), dpi=300)
    grid = ImageGrid(figure, 111, nrows_ncols=(4, 4), axes_pad=0.3)

    names = ["Initial", "Binarized", "Expanded Fibroblasts", "Closed (1)", "Dilated (1)", "Remove Small Objects",
             "Dilated (2)", "Closed (2)", "Segmented Airways", "Foreground Cleaned",
             "Background Cleaned", "Eroded", "Background Masked", "Inverted",
             "Labeled", "Colored"]

    for axes, image, name in zip(grid, subplots, names):
        axes.set_title(name)
        axes.imshow(image)

    if filename is not None:
        print(f"[+] Saving to: {filename}")
        plt.savefig(filename)

    # plt.show()


def display(image, name="", size=(12, 12), filename=None):
    """
    Display image on screen
    :param image: xp.xdarray, spatial data for a singular image
    :param name: str, title of image to display
    :param size: tuple, width and height of image in inches
    :param filename: str, name of file to save image to
    """
    print(f"[+] Displaying {name}")

    plt.figure(figsize=size, dpi = 300)
    plt.title(name)
    plt.imshow(image)

    if filename is not None:
        print(f"[+] Saving to: {filename}")
        plt.savefig(filename)

    # plt.show()

def visualize_cells_assigned_array(cells_assigned, name="Cells Assigned", size=(12, 12), dilation=10, filename=None, dpi=300):
    """
    Display an image of cells labeled by airspace affiliation
    :param cells_assigned: xp.xdarray, cells labeled by airspace affiliation
    :param name: str, title of image to display
    :param dilation: int, amount to grow cells in image for better visualization
    :param filename: str, name of file to save image to
    """
    cells_assigned = sk.morphology.dilation(cells_assigned, sk.morphology.disk(dilation))

    if USE_GPU: cells_assigned = cp.asnumpy(cells_assigned)

    colored = skimage.color.label2rgb(cells_assigned, colors=colors)
    del cells_assigned
    display(colored, name, size, filename)


def visualize_cells_assigned(spatial_data, shape, name="Cells Assigned", size=(12, 12), dilation=10, filename=None, dpi = 300):
    """
    Prepares a cells_assigned xp.xdarray from csv data or a df to be passed into
    visualize_cells_assigned_array function for displaying.
    :param spatial_data: str, path to csv file or a df
    :param shape: tuple, shape of labelmap
    :param name: str, title of image to display
    :param dilation: int, cells are dilated according to dilation so they become
    visible when displayed
    :param filename: name of file to save image to
    """
    if type(spatial_data) is str:
        df = prepare_df(spatial_data).copy()
    else:
        df = spatial_data.copy()
        arbitrary_scale = 2
        df["x_scaled_centroid"] = (df["x_centroid"].values * arbitrary_scale).astype(int)
        df["y_scaled_centroid"] = (df["y_centroid"].values * arbitrary_scale).astype(int)

    cells_assigned = xp.zeros(shape, dtype=df['lumen_id'].dtype)

    cells_assigned[df['x_scaled_centroid'], df['y_scaled_centroid']] = df['lumen_id']

    visualize_cells_assigned_array(cells_assigned, name, size, dilation, filename)

# Cell Annotation Functions
def calculate_distances(points1, points2):
    """
    Calculates all disances between two sets of points
    :param points1: xp.xdarray (2D), list of points of cells
    :param points2: xp.xdarray (2D), list of points in airspace
    :returns distances: xp.xdarray (2D), matrix of distances shape (N1, N2);
    each element represents the distance between a point in points1 and a
    point in points2
    """
    diff = points1[:, None, :] - points2[None, :, :]
    distances = xp.linalg.norm(diff, axis=-1)

    return distances

def process_batch(x_scaled_centroid, y_scaled_centroid, segmented, distance_threshold):
    """
    Assigns cells in provided batch with label of closest airspace within or at
    30px.
    :param x_scaled_centroid: list, list of all centroid x locations
    :param y_scaled_centroid: list, list of all centroud y locations
    :param segmented: xp.xdarray (2D), segmented labelmap
    :param distance_threshold: int, max distance from airspace allowed for the
    transcript to be assigned that airspace id
    :returns assignments: xp.xdarray (1D), list where each index corsponds to a cell
    in the batch, and each value is an airspace label
    """
    segment_labels = xp.unique(segmented)[1:]
    label_distances = xp.empty((0, len(x_scaled_centroid)))

    for label in segment_labels:
        points = xp.argwhere(segmented == label)
        distances = calculate_distances(xp.stack([x_scaled_centroid, y_scaled_centroid], axis=-1), points)
        min_distances = xp.min(distances, axis=1)
        label_distances = xp.vstack((label_distances, min_distances))

    closest_label = xp.argmin(label_distances, axis=0) + 1
    closest_distance = xp.min(label_distances, axis=0)
    distance_mask = closest_distance <= distance_threshold
    assignments = closest_label * distance_mask

    return assignments

def annotate_cells(csv_spatial_data_path, segmented, distance_threshold=30, percent_of_df=100):
    """
    Takes dataframe with spatial data, calculates each transcripts distance to
    the closest airpace, assigns that transcript an airspace label if that
    distance is within or at 30px. Transcripts are assigned within the dataframe
    in an added column titled "airspace_id."
    :param csv_spatial_data_path: str, path to csv file
    :param segmented: xp.xparray (2D), segmented labelmap
    :param batch_size: int, batch size for GPU processing
    :param distance_threshold: int, max distance a cell can be from an airspace
    to be assigned that airspace id
    :returns df: DataFrame, data frame with cells assigned in "airspace_id" column
    """
    print("[+] Running cell annotation")

    df = prepare_df_for_annotations(csv_spatial_data_path)
    defined_df = df.dropna(subset=['x_centroid', 'y_centroid'])
    unique_df = defined_df.drop_duplicates(subset='cell_id').copy()

    unique_df = downsize_df(unique_df, 100)

    # Pad image to account for rounding in centroid location?
    pad_x = (np.max(unique_df["x_scaled_centroid"]) + 1) - segmented.shape[0]
    pad_y = (np.max(unique_df["y_scaled_centroid"]) + 1) - segmented.shape[1]
    if pad_x < 1: pad_x = 0
    if pad_y < 1: pad_y = 0
    segmented = np.pad(segmented, ((0, pad_x), (0, pad_y)))

    unique_df['lumen_id'] = 0
    if USE_GPU: segmented = cp.asarray(segmented).astype(int)

    if USE_GPU:
        np_segmented = segmented.get()
    else:
        np_segmented = segmented
    max_dist = distance_threshold
    arbitrary_scale = 2
    np_segmented_realzise = skimage.transform.resize(np_segmented, (i // arbitrary_scale for i in segmented.shape),
                                                     order=0, preserve_range=True, anti_aliasing=False)
    # Grabbing metrics for each lumen
    props_table = skimage.measure.regionprops_table(
        np_segmented_realzise,
        properties=('area',
                    'area_filled',
                    'axis_major_length',
                    'axis_minor_length',
                    'eccentricity',
                    'feret_diameter_max',
                    'orientation',
                    'perimeter')
    )
    props_table_df = pd.DataFrame.from_dict(props_table)

    lumen_id_list = []
    lumen_dist = []
    for i in range(unique_df.shape[0]):
        first_cell = unique_df.iloc[i,:]
        cell_x_pos = first_cell['x_scaled_centroid']
        cell_y_pos = first_cell['y_scaled_centroid']

        x_min = cell_x_pos - max_dist
        y_min = cell_y_pos - max_dist
        if x_min < 0: x_min = 0
        if y_min < 0: y_min = 0

        x_max = cell_x_pos + max_dist
        y_max = cell_y_pos + max_dist
        if x_max > np_segmented.shape[1]: x_max = np_segmented.shape[0]
        if y_max > np_segmented.shape[0]: y_max = np_segmented.shape[1]

        search_area = np_segmented[x_min:x_max, y_min:y_max]
        has_marker = np.sum(search_area)

        if has_marker != 0:
            y_pxls, x_pxls = np.where(search_area > 0)
            coords = list(zip(y_pxls, x_pxls))

            x = [(max_dist, max_dist)]
            distance, index = spatial.KDTree(coords).query(x)
            if distance <= max_dist:
                px_color = search_area[coords[index[0]]]
                lumen_id_list.append(px_color)
                lumen_dist.append(distance[0] / arbitrary_scale)
            else:
                lumen_id_list.append(0)
                lumen_dist.append(0)
        else:
            lumen_id_list.append(0)
            lumen_dist.append(0)

    unique_df['lumen_id'] = lumen_id_list
    unique_df['lumen_dist'] = lumen_dist
    df_merge = df.merge(unique_df[['cell_id', 'lumen_id']], on='cell_id', how='left')
    df_merge2 = df_merge.merge(props_table_df, right_index=True, left_on = 'lumen_id', how='left')
    columns_to_remove = ['x_scaled_centroid', 'y_scaled_centroid', 'x_loc_int', 'y_loc_int']
    return df_merge2.loc[:, ~df_merge2.columns.isin(columns_to_remove)]


def prepare_df(csv_spatial_data_path, arbitrary_scale=2):
    """
    Reads csv contents into a data frame, scales centroids, and casts values to
    int.
    :param csv_spatial_data_path: str, path to csv file
    :param arbitrary_scale: int, scale factor for centroids
    :returns df: DataFrame, data frame with scaled int centroids
    """
    df = pd.read_csv(csv_spatial_data_path)

    df["x_scaled_centroid"] = (df["x_centroid"].values * arbitrary_scale).astype(int)
    df["y_scaled_centroid"] = (df["y_centroid"].values * arbitrary_scale).astype(int)

    return df


def downsize_df(df, percentage):
    """
    Downsizes data frame (for testing code more quickly)
    :param df: DataFrame, spatial data frame
    :param percentage: int, percentage of cells to use
    :returns df: DataFrame, downsized data frame
    """
    if percentage == 100:
        return df
    num_rows = int(len(df) * ((100 - percentage) / 100))
    rows = df.sample(num_rows)
    df = df.drop(rows.index)

    return df


base_path = sys.argv[1]
all_images = glob.glob(f"{base_path}/*.csv")

suffix_out = "out"

def run_image(image_path):
    print(f"Running {image_path}")
    folder = str(Path(image_path).parent)
    fname = str(Path(image_path).name)
    if not os.path.exists(f"{folder}/{suffix_out}/"):
        os.makedirs(f"{folder}/{suffix_out}/")
    if os.path.exists(f"{folder}/{suffix_out}/{fname[:-8]}_annotated_{suffix_out}.csv.gz"):
        print(f"Skipping:: {folder}/{suffix_out}/{fname[:-8]}_annotated_{suffix_out}.csv.gz")
        return
    image = load_spatial_data_clust_label(image_path, label_exclude=[9, 2],
                                          cache=False, recreate = True)
    segmented = run_segmentation(image)
    display(segmented[-1], filename = f"{folder}/{suffix_out}/{fname[:-8]}_segment_{suffix_out}.png")
    df = annotate_cells(image_path, segmented[-2], 30, 100)
    df.to_csv(f"{folder}/{suffix_out}/{fname[:-8]}_annotated_{suffix_out}.csv.gz", index=False)
    visualize_cells_assigned(df, segmented[-2].shape, "Cells Assigned", (12, 12), 10, f"{image_path[:-8]}_cells_{suffix_out}.png")

for image_path in all_images:
    try:
        run_image(image_path)
    except Exception as e:
        print(e)
