"""
map_of_Mexico.py

@author: Andrea Gregorio
@date: Agosto 2026

This script is for review the results of algorithm.py (background rate) u_xy and show them over a map of Mexico, 
because the data were from Guerrero
"""
import plotly.express as px
import pandas as pd

if __name__ == '__main__':
    cities = {"City":["1", "2", "3"], "State": ["Guerrero", "Ciudad de México", "Puebla"], "lat":[16.84942, 19.42847, 19.04778], "lon":[-99.90891, -99.12766, -98.20723]}
    df = pd.DataFrame(cities)
    fig = px.scatter_map(
        df,
        lat=df["lat"],
        lon=df["lon"],
        hover_name=df["City"],  # text that appears over 
        hover_data=df[["State"]],  # text that appears after latitude and longitude
        color_discrete_sequence=["red"],
        zoom=3,
            height=300,
    )
    fig.update_layout(map_style="open-street-map")
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})
    fig.update_layout(map_bounds={"west": -180, "east": -80, "south": 10, "north": 90})
    fig.show()
