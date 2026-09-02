"""
map_of_Mexico.py

@author: Andrea Gregorio
@date: Agosto 2026

This script is for review the results of algorithm.py (background rate) u_xy and show them over a map of Mexico, 
because the data were from Guerrero
"""
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import griddata

# for contourf
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np


def prepare_grid(gdf, file_earthquakes, file_uxy):
    xmin, ymin, xmax, ymax = gdf.total_bounds
    X, Y = np.meshgrid(np.linspace(xmin, xmax, 256), np.linspace(ymin, ymax, 256))
    Z = -X**2 - Y**2  # this function is the result from evaluate u_xy over the grid
    return X, Y, Z

def show_grid_results(gdf, X, Y, Z):
    plt.style.use('_mpl-gallery-nogrid')
    levels = np.linspace(Z.min(), Z.max())
    fig, ax = plt.subplots()
    ax.contourf(X, Y, Z, levels=levels)
    gdf.boundary.plot(ax=ax, color='black')
    plt.show()


if __name__ == '__main__':
    earthquakes = "Earthquakes.csv"
    uxy = "U_xy.csv"
    shp_mexico = "Mapa base a nivel estatal y mapa general. Formato Raster/mbtifgw.shp"

    gdf = gpd.read_file(shp_mexico)
    x, y, z = prepare_grid(gdf, earthquakes, uxy)
    show_grid_results(gdf, x, y, z)
