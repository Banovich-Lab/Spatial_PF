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
from functools import reduce

# import ray
# ray.init()

USE_GPU = False
# These are imported in this way so that some operations can be forced to CPU even when GPU is used for compute
if USE_GPU:
    # from cucim import skimage as gpu_skimage
    import cupy as cp
    from cupyx.scipy import spatial, ndimage

    # sk: gpu_skimage = gpu_skimage
    sk: skimage = skimage
    xp: cp = cp
else:
    from scipy import spatial, ndimage
    sk: skimage = skimage
    xp: np = np


#Display Settings
# Create a list of more colors to apply to labels
cm = plt.get_cmap("tab20b")
colors = [cm(i) for i in range(cm.N)]
cm2 = plt.get_cmap("tab20c")
colors2 = [cm2(i) for i in range(cm2.N)]
colors.extend(colors2)

plt.rcParams["figure.figsize"] = [8, 8]
# Display Functions


def prepare_df_for_annotations(csv_spatial_data_path):
    """
    Reads csv contents into a data frame, scales centroids, and casts values to
    int.
    :param csv_spatial_data_path: str, path to csv file
    :param arbitrary_scale: int, scale factor for centroids
    :returns df: DataFrame, data frame with scaled int centroids
    """
    if type(csv_spatial_data_path) is str:
        df = pd.read_csv(csv_spatial_data_path)
    else:
        df = csv_spatial_data_path.copy()
    x_min = np.min(df["x_location"].values.astype(float))
    y_min = np.min(df["y_location"].values.astype(float))
    df = df.dropna(subset=['x_centroid', 'y_centroid'])
    df["x_loc_int_norelative"] = df["x_location"].values.astype(float) - x_min
    df["y_loc_int_norelative"] = df["y_location"].values.astype(float) - y_min
    df["x_loc_int"] = (df["x_loc_int_norelative"]).astype(int)
    df["y_loc_int"] = (df["y_loc_int_norelative"]).astype(int)

    df["x_scaled_centroid_norelative"] = df["x_centroid"].values.astype(float) - x_min
    df["y_scaled_centroid_norelative"] = df["y_centroid"].values.astype(float) - y_min
    df["x_scaled_centroid"] = (df["x_scaled_centroid_norelative"]).astype(int)
    df["y_scaled_centroid"] = (df["y_scaled_centroid_norelative"]).astype(int)

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


def calc_hist_fullsize(filtered_data, output_shape: tuple[int, int] = ()):
    xmin_out = 0
    xmax_out = output_shape[0]
    ymin_out = 0
    ymax_out = output_shape[1]

    filt_xmin = np.min(filtered_data['x_loc_int'])
    filt_xmin = filt_xmin if filt_xmin > 0 else 0
    filt_xmax = np.max(filtered_data['x_loc_int'])
    filt_ymin = np.min(filtered_data['y_loc_int'])
    filt_ymin = filt_ymin if filt_ymin > 0 else 0
    filt_ymax = np.max(filtered_data['y_loc_int'])

    pad_x_l = filt_xmin
    pad_x_r = xmax_out - filt_xmax
    pad_y_l = filt_ymin
    pad_y_r = ymax_out - filt_ymax

    density_filt = sk.exposure.rescale_intensity(
        np.sqrt(np.histogram2d(filtered_data['x_location'], filtered_data['y_location'],
                               bins=[filt_xmax - filt_xmin, filt_ymax - filt_ymin], density=True)[0])
    )
    out = np.pad(density_filt, ((pad_x_l, pad_x_r), (pad_y_l, pad_y_r)))
    return out


def load_spatial_data_density(filename: str, label_exclude: list[int, ...] = None, cache: bool = True, recreate: bool = False) -> np.ndarray:
    """
    Load a spatial transcriptome file and return a 3D array, first dimension is
    the transcriptional cluster
    :param filename: str, input file to be rad with pd.read_csv
    :param label_exclude: list, which clusters/labels to exclude
    :param cache: bool, store a .npy of the image for faster loading later
    :param recreate: bool, overwrite and recreate the cached file
    :return: xp.xdarray (3D), dimensions are the transcriptional cluster, last slot is the sum, other slots match
    location, and the y location respectively
    """
    cache_fname = os.path.splitext(filename)[0] + "_density.npy"

    if (cache and not os.path.exists(cache_fname)) or recreate:
        # if DEBUG:
        #     df = prepare_df_for_annotations(downsize_df(pd.read_csv(filename), 10))
        # else:
        df = prepare_df_for_annotations(filename)
        xmin = np.min(df['x_loc_int'])
        xmin = xmin if xmin > 0 else 0
        xmax = np.max(df['x_loc_int'])
        ymin = np.min(df['y_loc_int'])
        ymin = ymin if ymin > 0 else 0
        ymax = np.max(df['y_loc_int'])

        # Create a new empty image that is the right size to hold the data
        img_arr = np.zeros((np.max(df["gmm12"]) + 1, np.max(df["x_loc_int"]), np.max(df["y_loc_int"]))).astype(float)

        for label in np.unique(df.loc[:, "gmm12"]):
            if label_exclude is not None and label in label_exclude:
                print(f"Skipping {label}")
                continue
            # Vectorized image load operation, specify X and Y locations in the slice operation
            df_match_label = df.loc[df.loc[:, "gmm12"] == label, ]
            density = calc_hist_fullsize(df_match_label, (xmax, ymax))
            img_arr[label, :] = density

        img_arr[-1, :] = calc_hist_fullsize(df, (xmax, ymax))

        np.save(cache_fname, img_arr)
    else:
        print("[+] Loading cached file")

        img_arr = np.load(cache_fname)

    return img_arr


def fastplot(data):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 600
    plt.figsize=(10.0, 10.0)
    fig, ax = plt.subplots()
    # ax.imshow(np.rot90(data), cmap=plt.cm.gist_earth_r,)
    ax.imshow(data, cmap=plt.cm.gist_earth_r,)
    # ax.plot(df['x_loc_int'], df['y_loc_int'], 'k.', markersize=0.01)
    ax.set_xlim([0, data.shape[-1]])
    ax.set_ylim([data.shape[-2], 0])
    plt.show()


def get_disk(size: int = 10):
    if USE_GPU:
        return cp.asarray(sk.morphology.disk(size))  # Send data to GPU
    else:
        return sk.morphology.disk(size)


def mask_alphashape_binary(img, alpha: float = 0.1, downscale: int = 2):
    '''
    input is a (mostly) closed binary array
    :param img:
    :param alpha:
    :param downscale:
    :return: binary ndarray of masked shapes
    '''
    print("[+] Masking with alphashapes")
    # img = connect_close_spaces.copy()
    img = img.astype(bool)

    pad_size = 50
    # Pad so that edges are handled correctly
    img_pad = xp.pad(img, ((pad_size, pad_size),(pad_size, pad_size)))
    downscale_factor = downscale  # Increase for speed, decrease for accuracy
    img_downscale = sk.transform.resize(img_pad, [i // downscale_factor for i in img_pad.shape],
                                        order=0, preserve_range=True, anti_aliasing=False).astype(int)
    img_downscale = ndimage.binary_closing(img_downscale, get_disk(5))
    img_downscale = sk.morphology.remove_small_objects(img_downscale, 100)
    img_labeled = sk.measure.label(img_downscale)
    # fastplot(img_labeled)
    if USE_GPU: img_downscale = cp.asnumpy(img_downscale)

    # Split each connected, labeled, object into a different matrix
    label_list = []
    img_list = []
    for label in range(1, np.max(img_labeled) + 1):
        individual_obj = img_labeled.copy()
        individual_obj[individual_obj != label] = 0
        img_list.append(individual_obj)
        # mask is a 2D boolean np.array containing the segmentation result
        contour_img = np.logical_xor(individual_obj, skimage.morphology.binary_erosion(individual_obj))
        contour_y, contour_x = np.where(contour_img)
        label_list.append((contour_y, contour_x))

    img_out = np.zeros(img_downscale.shape, dtype=bool)

    for label_idx, labelmask in enumerate(label_list):
        points = np.array([labelmask[0], labelmask[1]]).T

        # Use an alphashape and decrease its sensitivity until a single polygon is found
        # At worst (alpha = 0), this uses the convex hull of the connected object
        # At best, this is a close approximation of the shape
        iter_step = 0.01
        alpha = alpha  # Starting at 0.1, airways are bigger than whole tissues
        alpha_shape = alphashape.alphashape(points, alpha)
        in_img = img_list[label_idx] > 0

        while type(alpha_shape) != shapely.geometry.polygon.Polygon:
            alpha = alpha - iter_step

            if alpha == 0: print("No ideal solution found, using convex hull")

            alpha_shape = alphashape.alphashape(points, alpha)

        # Create a mask from the polygon, not implemented for GPU
        # Create a mask from the polygon, not implemented for GPU
        raster_mask = skimage.draw.polygon2mask(img_downscale.shape, np.asarray(list(zip(*alpha_shape.exterior.xy))))
        raster_mask = np.bitwise_or(raster_mask, in_img)  # Make sure we don't drop anything
        img_out = np.bitwise_or(img_out, raster_mask)  # Assemble all objects into the final mask

    # Resize output
    if USE_GPU: img_out = cp.asarray(img_out) # Return to GPU for downstream if needed
    scaleup = sk.transform.resize(img_out, img_pad.shape,
                                  order=0, preserve_range=True, anti_aliasing=False).astype(bool)

    scaleup = scaleup[pad_size:-pad_size, pad_size:-pad_size]

    return scaleup


def segment_airways(spatial_data, airway_cluster = 6):
    """
    Draws alphashape around airways so that they are always enclosed
    :param spatial_data: xp.xdarray (2D), spatial data for a singular image
    :return: xp.xdarray (2D), mask of airway tissue and lumen
    """

    if USE_GPU:
        only_airway = cp.asarray(spatial_data[airway_cluster, :]) # Send data to GPU
    else:
        only_airway = spatial_data[airway_cluster, :]

    blur = sk.exposure.rescale_intensity(sk.filters.gaussian(only_airway, sigma = 5))

    closed = ndimage.grey_closing(blur, footprint = get_disk(10))
    binarize = closed > 0.15
    small_objects_removed = sk.morphology.remove_small_objects(binarize, 2500)

    masked_airway = mask_alphashape_binary(small_objects_removed, alpha=1, downscale=2)
    eroded = ndimage.binary_erosion(masked_airway, get_disk(15))
    masked = np.bitwise_or(eroded, small_objects_removed)
    presumed_lumen = np.bitwise_and(masked, np.invert(ndimage.binary_erosion(masked, get_disk(50))))

    return presumed_lumen, masked, masked_airway


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


def expand_fibroblasts(spatial_data, dialation = 10, fibroblast_clusts: list[int, int] | None = None) -> list[np.ndarray, ...]:
    """
    Expand signals from fibroblasts (clusters 5 and 6)
    :param spatial_data: xp.xdarray (2D), spatial data for a singular image
    :param dialation: sk.morphology footprint for dialation shape and size
    :param fibroblast_clusts: list of clusters
    :return: list of xp.array with n corresponding in n clusters
    """
    if fibroblast_clusts is None:
        return spatial_data

    out = []
    for clust in fibroblast_clusts:
        if USE_GPU:
            only_selected = cp.asarray(spatial_data[clust, :]) # Send data to GPU
        else:
            only_selected = spatial_data[clust, :]

        blur = sk.exposure.rescale_intensity(sk.filters.gaussian(only_selected, sigma = 15))
        closed = ndimage.grey_dilation(blur, footprint = get_disk(dialation))
        binarize = closed > 0.15
        out.append(ndimage.binary_closing(binarize, get_disk(dialation)))

    return out


# Segmention Pipeline
def run_segmentation(spatial_data: np.ndarray) -> tuple[np.ndarray, tuple[np.ndarray]]:
    """
    Run segmentation on spatial data
    :param spatial_data: xp.xdarray, spatial data for a singular image
    :return: list, list of labelmaps (xp.xdarray) from each step in segmentation
    pipeline
    """
    # spatial_data = image_array.copy()
    print("[+] Running segmentation")
    airway_cluster = 6
    fb_clusts = [8, 5]
    immune_clust = 7

    print("[+] Segmenting airways")
    airways = segment_airways(spatial_data, airway_cluster = airway_cluster)

    print("[+] Expanding fibroblasts")
    bigger_fibroblasts = expand_fibroblasts(spatial_data, dialation=20, fibroblast_clusts = fb_clusts)

    all_others_noimmune = np.delete(spatial_data, [airway_cluster, *fb_clusts, immune_clust], axis = 0)
    merged = np.sum(all_others_noimmune, axis=0)
    blur = sk.exposure.rescale_intensity(sk.filters.gaussian(merged, sigma = 10))

    all_dialation = 5  # The bigger this is, the more 'merged' things will be, but smaller airspaces get split
    closed = ndimage.grey_dilation(blur, footprint = get_disk(all_dialation))
    binarize_all = closed > 0.15

    print("[+] Combining masks")
    binarize_all_together = reduce(lambda x, y: np.bitwise_or(x, y), [binarize_all, airways[0], *bigger_fibroblasts])

    foreground_cleaned = sk.morphology.remove_small_objects(binarize_all_together, 1000)
    background_cleaned = sk.morphology.remove_small_holes(foreground_cleaned, 250)

    connect_close_spaces = ndimage.binary_erosion(background_cleaned, get_disk(5))
    # fastplot(connect_close_spaces)
    # connect_close_spaces = background_cleaned
    mask = mask_alphashape_binary(connect_close_spaces, downscale=10, alpha=0.03)
    mask_erode = ndimage.binary_erosion(mask, get_disk(25))
    inverse_mask = xp.invert(mask_erode)

    background_masked = xp.bitwise_or(connect_close_spaces, inverse_mask)

    inverted = xp.invert(background_masked)
    labeled = sk.measure.label(inverted)
    print("[+] Labeling mask")
    labeled_dialated = ndimage.maximum_filter(labeled, 15)
    labeled_dialated[labeled != 0] = labeled[labeled != 0]

    labeled_airways = sk.morphology.label(np.logical_or(airways[1], airways[2]))
    labeled_updated_airway = labeled_dialated.copy()
    if len(np.unique(labeled_airways)) > 1:
        print("[+] Updating airway masks")
        for airway in np.unique(labeled_airways)[1:]:
            select_airway_mask = labeled_airways.copy()
            select_airway_mask[np.where(labeled_airways != airway)] = 0
            select_airway_mask = ndimage.binary_erosion(select_airway_mask.astype(bool), get_disk(15))
            labs_in_airway = labeled_dialated.copy()
            labs_in_airway[np.invert(select_airway_mask.astype(bool))] = 0

            labs = labs_in_airway[labs_in_airway != 0]
            counts = np.bincount(labs)
            if len(counts) > 0:
                most_common_label = np.argmax(counts)
                # Binary dialate, then apply same label to all
                dialated_mask = ndimage.binary_dilation(labs_in_airway.astype(bool), get_disk(10))
                mask_label_merge = np.bitwise_or(select_airway_mask, dialated_mask)
                mask_label_merge = ndimage.binary_erosion(mask_label_merge, get_disk(5))
                labeled_updated_airway[mask_label_merge] = most_common_label
        final_out = labeled_updated_airway
        # fastplot(final_out)
    else:
        final_out = labeled_dialated
    print("[+] Creating color image")
    colored = skimage.color.label2rgb(final_out, colors=colors)  # Don't run this on GPU, needs too much memory

    return (final_out, (airways, bigger_fibroblasts, all_others_noimmune,
                        merged, blur, binarize_all_together,
                        foreground_cleaned, background_cleaned,
                        background_masked, labeled, labeled_dialated, labeled_updated_airway,
                        colored))



def display_all(subplots, filename = None):
    """
    Display all steps of processing in a table (designed to assist parameter tuning)
    :param subplots: list, list of each output of each step in the segmentation process
    :param filename: str, name of file to save image to
    """
    print("[+] Displaying all steps in segmentation")

    figure = plt.figure(figsize=(12, 12), dpi=300)
    grid = ImageGrid(figure, 111, nrows_ncols=(4, 4), axes_pad=0.3)

    # (airways, bigger_fibroblasts, all_others_noimmune,
    #  merged, blur, binarize_all_together,
    #  foreground_cleaned, background_cleaned, connect_close_spaces,
    #  background_masked, labeled, labeled_dialated, labeled_updated_airway,
    #  colored)
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
    if USE_GPU: cells_assigned = cp.asarray(cells_assigned)
    cells_assigned = sk.morphology.dilation(cells_assigned, sk.morphology.disk(dilation))

    if USE_GPU: cells_assigned = cp.asnumpy(cells_assigned)

    colored = skimage.color.label2rgb(cells_assigned, colors=colors)
    del cells_assigned
    display(colored, name, size, filename)


def visualize_cells_assigned(spatial_data, name="Cells Assigned", size=(12, 12), dilation=10, filename=None, dpi = 300):
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

    df = prepare_df_for_annotations(spatial_data).copy()
    df = df.dropna(subset=['x_centroid', 'y_centroid'])
    df = df.drop_duplicates(subset='cell_id')
    array_shape = (np.max(df['x_scaled_centroid'])  + 1, np.max(df['y_scaled_centroid']) + 1)
    cells_assigned = xp.zeros(array_shape, dtype=df['lumen_id'].dtype)

    cells_assigned[df['x_scaled_centroid'], df['y_scaled_centroid']] = df['lumen_id']

    visualize_cells_assigned_array(cells_assigned, name, size, dilation, filename)


def annotate_cells(csv_spatial_data_path, segmented, distance_threshold=20):
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
    # segmented = segmented[-2].copy()
    # image_path, segmented[-2], 30, 100
    # csv_spatial_data_path = image_path
    # distance_threshold=30
    # percent_of_df=10
    df = prepare_df_for_annotations(csv_spatial_data_path)
    np.max(df["x_scaled_centroid"])
    defined_df = df.dropna(subset=['x_centroid', 'y_centroid'])
    unique_df = defined_df.drop_duplicates(subset='cell_id').copy()

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

    # Grabbing metrics for each lumen
    props_table = skimage.measure.regionprops_table(
        np_segmented,
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

    print("[+] Finding label of neighboring lumens")
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

        if x_min < 0: x_min = 0
        if y_min < 0: y_min = 0

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
                lumen_dist.append(distance[0])
            else:
                lumen_id_list.append(0)
                lumen_dist.append(0)
        else:
            lumen_id_list.append(0)
            lumen_dist.append(0)

    print("[+] Merging IDs into dataframe")
    unique_df['lumen_id'] = lumen_id_list
    unique_df['lumen_dist'] = lumen_dist
    df_orig = pd.read_csv(csv_spatial_data_path)
    df_merge = df_orig.merge(unique_df[['cell_id', 'lumen_id']], on='cell_id', how='left')
    df_merge2 = df_merge.merge(props_table_df, right_index=True, left_on = 'lumen_id', how='left')
    columns_to_remove = ['x_scaled_centroid', 'y_scaled_centroid', 'x_loc_int', 'y_loc_int',
                         'x_loc_int_norelative', 'y_loc_int_norelative',
                         'x_scaled_centroid_norelative', 'y_scaled_centroid_norelative']
    return df_merge2.loc[:, ~df_merge2.columns.isin(columns_to_remove)]



suffix_out = "out"
# @ray.remote
def run_image(image_path):
    print(f"Running {image_path}")
    try:
        folder = str(Path(image_path).parent)
        fname = str(Path(image_path).name)
        if not os.path.exists(f"{folder}/{suffix_out}/"):
            os.makedirs(f"{folder}/{suffix_out}/")
        if os.path.exists(f"{folder}/{suffix_out}/{fname[:-4]}_cells_{suffix_out}.png"):
            print(f"Skipping:: {folder}/{suffix_out}/{fname[:-4]}_annotated_{suffix_out}.csv.gz")
            return
        image_array = load_spatial_data_density(image_path, cache=True)
        visualize_cells_assigned_array(image_array[-1] * 255, dilation = 1, filename = f"{folder}/{suffix_out}/{fname[:-4]}_all_{suffix_out}.png", name = "Transcripts")
                                       # name=f"Transcripts. Sum = {np.sum(image_array[-1])}, Mean = {np.mean(image_array[-1])}, stdev = {np.std(image_array[-1])}")

        segmented = run_segmentation(image_array)
        # Save the whole image output
        display(segmented[1][-1], filename = f"{folder}/{suffix_out}/{fname[:-4]}_segment_{suffix_out}.png")
        df = annotate_cells(image_path, segmented[0], 20)
        # out_fname = f"{folder}/{suffix_out}/{fname[:-8]}_annotated_{suffix_out}.csv"
        df.to_csv(f"{folder}/{suffix_out}/{fname[:-4]}_annotated_{suffix_out}.csv.gz", index=False)
        visualize_cells_assigned(df, "Cells Assigned", (12, 12), 5, f"{folder}/{suffix_out}/{fname[:-4]}_cells_{suffix_out}.png")
    except Exception as e:
        print(f"Error: {e}")



base_path = "/data/aspera/new dataset jun 2024/upload_for_lumen_seg_June2024/"
all_images = glob.glob(f"{base_path}/*.csv")

for path in all_images:
    run_image(path)
