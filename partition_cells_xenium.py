#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 14:52:07 2023

@author: emee
"""

#This script is intended to be used with standard 10x Xenium output.

#This script partitions cells based on whether their centroid falls within
#an annotated region. It appends a column for each region to the metadata csv 
#containing a bool for each cell. If true, the cell centroid falls within the 
#annotated region. If false, the cell centroid does not fall within the 
#annotated region.

#This script assumes that:
#Annotations are in geojson format.
#Each geojson file corresponds to one annotation.
#All geojson files are in the same directory.

import pandas as pd
import geopandas as gpd
import shapely
import os

#Path to directory containing geojson files (one file per annotated region)
annotations_path = ""

#Path to cells file from Xenium output
cells_path = ""

#Origin of coordinate system
x_origin = 0
y_origin = 0
 

x_origin_transform = x_origin * 0.2125
y_origin_transform = y_origin * 0.2125
#Load cell metadata csv (assumed to be in same directory as geojsons)
cells = pd.read_csv(cells_path)
#Create list of geojson files found in directory
files = []
for file in os.listdir(annotations_path):
    if file.endswith(".geojson"):
        files.append(file)
        
#Loops through each geojson in list
for file in files:
    #Load and scale geojson (Xenium is always .2125 microns/pixel)
    geojson = gpd.read_file(annotations_path + "/" + file)
    scaled_geojson = geojson.scale(origin=(0, 0), xfact = .2125, yfact = .2125)
    #Create shapely points of cell centroids
    points = shapely.points(cells["x_centroid"] - x_origin_transform, cells["y_centroid"] - y_origin_transform)
    
    #Check if point is within geometry
    point_list = []
    for point in points:
        point_list.append(shapely.contains(scaled_geojson, point))
    
    #Cast to bool and append column to cell metadata 
    #Column name is geojson file name.   
    point_list_float = list(map(float, point_list))
    point_list_bool = list(map(bool, point_list_float))
    cells[file] = point_list_bool

#Write new metadata csv to same directory as geojsons
cells.to_csv(annotations_path + "/partitioned_cells.csv")
