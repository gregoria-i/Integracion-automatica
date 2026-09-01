"""
map_of_Mexico.py

@author: Andrea Gregorio
@date: Agosto 2026

This script is for review the results of algorithm.py (background rate) u_xy and show them over a map of Mexico, 
because the data were from Guerrero
"""
import plotly.express as px
import pandas as pd


def prepare_data(file_earthquakes, file_uxy, M0):
    # Earthquakes df and u_xy df joined
    df = pd.read_csv(file_earthquakes)
    u_xy_series = pd.read_csv(file_uxy)

    df = df[df['Magnitude']>= M0].copy()
    df.reset_index(drop=True, inplace=True)
    print(df.head())

    df['u_xy'] = u_xy_series['u_xy']

    return df

def show_results(df: pd.DataFrame):
    fig = px.scatter_map(df, lat=df["Latitude"], lon=df["Longitude"], color=df["u_xy"],
                         hover_name=df["u_xy"],  # text that appears over 
                         hover_data=df[["Magnitude", "Year", "Month", "Day"]],  # text that appears after latitude and longitude
                         )
    fig.update_layout(map_style="open-street-map")
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})
    fig.update_layout(map_bounds={"west": -180, "east": -80, "south": 10, "north": 90})
    fig.show()


if __name__ == '__main__':
    earthquakes = "Earthquakes.csv"
    uxy = "U_xy.csv"

    df = prepare_data(earthquakes, uxy, M0=4.3)
    show_results(df)
