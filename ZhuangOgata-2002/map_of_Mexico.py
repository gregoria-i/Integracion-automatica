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
from shapely.geometry import box

from algorithm_1 import ETAS_Declustering
from datetime import datetime

def prepare_grid(gdf, file_earthquakes):
    xmin = -106
    xmax = -96
    ymin = 15
    ymax = 20

    area = box(xmin, ymin, xmax, ymax)  # reduce the area to create a smaller grid
    gdf = gdf.clip(area)

    X, Y = np.meshgrid(np.linspace(xmin, xmax, 256), np.linspace(ymin, ymax, 256))
    obj = ETAS_Declustering(file_earthquakes)
    Z = obj.evaluate_u_over_grid(X, Y)
    return X, Y, Z, gdf

def show_grid_results(gdf, X, Y, Z):
    # plt.style.use('_mpl-gallery-nogrid')
    levels = np.linspace(Z.min(), Z.max())
    fig, ax = plt.subplots()
    plt.contour(X, Y, Z, levels=levels, cmap='inferno')
    plt.colorbar()
    plt.grid()
    plt.title("Background intensity")
    gdf.boundary.plot(ax=ax, color='black')
    plt.show()


if __name__ == '__main__':
    inicio = datetime.now()
    earthquakes = "Earthquakes.csv"
    shp_mexico = "Mapa base a nivel estatal y mapa general. Formato Raster/mbtifgw.shp"

    gdf = gpd.read_file(shp_mexico)
    x, y, z, gdf = prepare_grid(gdf, earthquakes)
    np.savetxt("U_grid.csv", z, delimiter=",",
               comments="", fmt="%.15g")
    fin = datetime.now()
    print(f"Tiempo utilizado: {str(fin-inicio)}")
    show_grid_results(gdf, x, y, z)
