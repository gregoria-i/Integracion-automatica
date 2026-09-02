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


def prepare_data(file_earthquakes, file_uxy, M0):
    # Earthquakes df and u_xy df joined
    df = pd.read_csv(file_earthquakes)
    u_xy_series = pd.read_csv(file_uxy)

    df = df[df['Magnitude']>= M0].copy()
    df.reset_index(drop=True, inplace=True)

    df['u_xy'] = u_xy_series['u_xy']

    return df

def show_results(df: pd.DataFrame):
    # Points in a color scale
    #fig = px.density_map(df, lat=df["Latitude"], lon=df["Longitude"], z=df["u_xy"], 
    fig = px.scatter_map(df, lat=df["Latitude"], lon=df["Longitude"], color=df["u_xy"],
    #fig = px.scatter(df, y=df["Latitude"], x=df["Longitude"], color=df["u_xy"],
                         hover_name=df["u_xy"],  # text that appears over 
                         hover_data=df[["Magnitude", "Year", "Month", "Day"]],  # text that appears after latitude and longitude
                         color_continuous_scale='plasma'
                         )

    fig.update_layout(map_style="open-street-map")
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})
    fig.update_layout(map_bounds={"west": -180, "east": -80, "south": 10, "north": 90})
    fig.show()

def show_grid_results(shp_file):
    plt.style.use('_mpl-gallery-nogrid')

    gdf = gpd.read_file(shp_file)
    #  print(gdf.crs)  # type of coordenates
    xmin, ymin, xmax, ymax = gdf.total_bounds
    X, Y = np.meshgrid(np.linspace(xmin, xmax, 256), np.linspace(ymin, ymax, 256))
    Z = -X**2 - Y**2  # this cannot be the same value for all the grid, because levels needs increasing values

    levels = np.linspace(Z.min(), Z.max())

    fig, ax = plt.subplots()

    ax.contourf(X, Y, Z, levels=levels)
    
    gdf.boundary.plot(ax=ax, color='black')
    ax.set_aspect('equal')
    plt.show()


if __name__ == '__main__':
    earthquakes = "Earthquakes.csv"
    uxy = "U_xy.csv"
    shp_mexico = "Mapa base a nivel estatal y mapa general. Formato Raster/mbtifgw.shp"

    show_grid_results(shp_mexico)
