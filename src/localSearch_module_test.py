from random import randint
import pandas as pd

import geodataframe_processing as gdfp
from custom_plotting import plot_initial_final_satellites

# Constants
EARTH_SURFACE_AREA_SQ_KM = 509600000
EARTH_LAND_AREA_SQ_KM = 148326000
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Global Parameters
BUFFER_RADIUS = 12065
PERTURB_DISTANCE_KM = (
    500  # we may want to make this something that decays exponentially
)

def local_search_optimization(
    satellite_gdf,
    map,
    new_satellite_column,
    buffer_radius,
    num_iterations=10,
):
    best_gdf = satellite_gdf.copy()
    best_land_coverage = gdfp.calculate_land_coverage(best_gdf, map, buffer_radius)

    for iteration in range(num_iterations):
        print(f"Iteration {iteration + 1}/{num_iterations}")

        # Perturb positions of new satellites only
        new_gdf = best_gdf.copy()
        gdfp.perturb_positions(
            new_gdf,
            new_satellite_column,
            eq_dist_epsg=EQUAL_DISTANCE_EPSG,
            eq_area_epsg=EQUAL_AREA_EPSG,
            max_shift_km=PERTURB_DISTANCE_KM,
            random_state=randint(0, 10000),
        )

        # Calculate new land coverage area
        new_land_coverage = gdfp.calculate_land_coverage(
            new_gdf, land_map, buffer_radius
        )

        if new_land_coverage > best_land_coverage:
            print(
                f"Improvement found: {new_land_coverage:.2f} km² vs {best_land_coverage:.2f} km²"
            )
            best_gdf, best_land_coverage = new_gdf, new_land_coverage
        else:
            print("No improvement in land coverage.")

    print(f"Optimized land coverage area: {best_land_coverage:.2f} km²")
    return best_gdf, best_land_coverage


# Load existing satellite data from TLEs and flag new satellites
existing_satellites_gdf = gdfp.load_existing_satellites(EPSG=EQUAL_AREA_EPSG)
existing_satellites_gdf["new_satellite"] = False  # Mark existing satellites

# Load map data
land_map, ocean_map = gdfp.load_land_ocean_data(EPSG=EQUAL_AREA_EPSG)

# Generate new satellites, create a copy of them for plotting, and combine with existing ones
new_satellites_gdf = gdfp.generate_new_satellites(
    EPSG=EQUAL_AREA_EPSG,
    num_satellites=30,
    input_map=land_map,
    true_random=False,
    sample_separate=False,
)
initial_satellites_gdf = new_satellites_gdf.copy()
satellite_gdf = pd.concat([existing_satellites_gdf, new_satellites_gdf])

# Run local search with perturbation targeted to new satellites
optimized_gdf, optimized_land_coverage = local_search_optimization(
    satellite_gdf,
    land_map,
    "new_satellite",
    num_iterations=30,
    buffer_radius=BUFFER_RADIUS,
)

# extracting the new satellites from the GeoDataFrame for plotting purposes
final_satellites_gdf = optimized_gdf[optimized_gdf["new_satellite"]].copy()

# re-project data onto a geodetic reference system for printing
# and updating coordinates
optimized_gdf_geodetic = optimized_gdf.to_crs(epsg=4326)
optimized_gdf_geodetic_coordinates = optimized_gdf_geodetic.get_coordinates()
optimized_gdf_geodetic["Latitude"] = optimized_gdf_geodetic_coordinates.y
optimized_gdf_geodetic["Longitude"] = optimized_gdf_geodetic_coordinates.x

# Print out details of the newly generated satellites
print("\nNew Satellites' Final Positions After Optimization:")
print(
    optimized_gdf_geodetic[optimized_gdf_geodetic["new_satellite"]][
        ["Satellite Name", "Latitude", "Longitude", "geometry"]
    ]
)

plot_initial_final_satellites(existing_satellites=existing_satellites_gdf,
                              initial_positions=initial_satellites_gdf,
                              final_positions=final_satellites_gdf,
                              land=land_map,
                              ocean=ocean_map)