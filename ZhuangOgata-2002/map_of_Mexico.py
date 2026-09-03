"""
map_of_Mexico.py

@author: Andrea Gregorio
@date: Agosto 2026

This script is for review the results of algorithm.py (background rate) u_xy and show them over a map of Mexico, 
because the data were from Guerrero
"""
# for contourf
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from algorithm_1 import ETAS_Declustering


def prepare_grid(gdf, file_earthquakes):
    xmin, ymin, xmax, ymax = gdf.total_bounds
    X, Y = np.meshgrid(np.linspace(xmin, xmax, 256), np.linspace(ymin, ymax, 256))

    obj = ETAS_Declustering(file_earthquakes)
    Z = obj.evaluate_u_over_grid(X, Y)

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
    x, y, z = prepare_grid(gdf, earthquakes)
    show_grid_results(gdf, x, y, z)
