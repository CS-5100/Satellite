import matplotlib.pyplot as plt
import geopandas as gpd

def plot_initial_final_satellites(existing_satellites: gpd.GeoDataFrame,
                                  initial_positions: gpd.GeoDataFrame,
                                  final_positions: gpd.GeoDataFrame,
                                  land: gpd.GeoDataFrame,
                                  ocean: gpd.GeoDataFrame,
                                  buffer=None):
    """Generates two plots of the earth with satellites on it, the top plot mapping initial satellite placements
    and the bottom plot mapping where the satellites ended up after processing

    Args:
        existing_satellites (gpd.GeoDataFrame): GeoDataFrame defining the satellites already in the satellite constellation
        initial_positions (gpd.GeoDataFrame): GeoDataFrame defining the initial positions of the proposed candidate satellites
        final_positions (gpd.GeoDataFrame): GeoDataFrame defining the final positions of the proposed candidate satellites
        land (gpd.GeoDataFrame): GeoDataFrame defining the land masses
        ocean (gpd.GeoDataFrame): GeoDataFrame defining the ocean
        buffer (_type_, optional): The radius of the circles to plot around each satellite in meters.
            If none, the points shown do not correlate to actual coverage areas. Defaults to None.
    """
    initial_fig, initial_ax = plt.subplots(2, 1)
    before_ax = initial_ax[0]
    after_ax = initial_ax[1]

    # add map data to both plots
    land.plot(ax=before_ax, color="#228B22")
    ocean.plot(ax=before_ax, color="#246BCE")
    land.plot(ax=after_ax, color="#228B22")
    ocean.plot(ax=after_ax, color="#246BCE")

    # if a plot of buffered circles is requested
    if buffer is not None:

        # copy current GeoDataFrames to preserve their initial state
        existing_satellites_buffered = existing_satellites.copy()
        initial_satellites_buffered = initial_positions.copy()
        final_satellites_buffered = final_positions.copy()

        # buffer the Point objects they create by the pre-specified buffer
        existing_satellites_buffered["geometry"] = existing_satellites_buffered[
            "geometry"
        ].buffer(buffer)
        initial_satellites_buffered["geometry"] = initial_satellites_buffered[
            "geometry"
        ].buffer(buffer)
        final_satellites_buffered["geometry"] = final_satellites_buffered[
            "geometry"
        ].buffer(buffer)

        # add the initial satellites to the initial map
        existing_satellites_buffered.plot(
            ax=before_ax, color="#0E0E10"
        )  # Jet Black
        initial_satellites_buffered.plot(ax=before_ax, color="#FF4F00")

        # add the final satellites to the final map
        existing_satellites_buffered.plot(ax=after_ax, color="#0E0E10")  # Jet Black
        final_satellites_buffered.plot(ax=after_ax, color="#FF4F00")

    # otherwise only add points and set their markersize to be 1 for easy visibility
    else:
        existing_satellites.plot(
            ax=before_ax, color="#0E0E10", markersize=1
        )  # Jet Black
        # the below color is apparently known as International Orange (Aerospace) and used in the aerospace industry
        initial_positions.plot(ax=before_ax, color="#FF4F00", markersize=1)

        existing_satellites.plot(
            ax=after_ax, color="#0E0E10", markersize=1
        )  # Jet Black
        # the below color is apparently known as International Orange (Aerospace) and used in the aerospace industry
        final_positions.plot(ax=after_ax, color="#FF4F00", markersize=1)

    # show the plot
    plt.show()