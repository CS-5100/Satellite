import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.stats import median_abs_deviation
from geodataframe_processing import (
    load_satellites_from_file,
    load_land_ocean_data,
    generate_new_satellites,
    calculate_land_coverage,
)
from search_functions import hill_climbing, random_restart_simulated_annealing

# EPSG Codes
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Time the TLE data was collected
TLE_TIME = datetime(2024, 11, 14, 16, 12, 32, 202983)

# Other Global Parameters
PERTURB_DISTANCE_KM = 500


# Importing a static set of TLEs for sensitivity analysis
base_gdf = load_satellites_from_file(
    EPSG=EQUAL_AREA_EPSG,
    input_filename="parameter_sensitivity_analysis_tles.txt",
    time=TLE_TIME,
)

# Loading land coverage area
land_map, ocean_map = load_land_ocean_data(EPSG=EQUAL_AREA_EPSG)


# defining methods for looping across many episodes and gathering performance data
def buffer_variation_hc(
    existing_satellite_gdf,
    equal_area_epsg,
    equal_distance_epsg,
    input_map,
    perturb_dist,
    ls_iter=10,
    initial_buffer=12065,
    multipliers=[0.5, 1, 2, 10, 20],
    iterations_per_multiplier=3,
    print_progress=True,
):

    earth_land_area_square_km = 148326000

    # initializing lists for the recording of parameters

    # location parameters (medians)
    average_durations = []
    average_total_coverage = []
    average_added_coverage_above_existing = []
    average_added_coverage_above_initial = []
    average_percent_coverage = []
    average_percent_above_existing = []
    average_percent_above_initial = []

    # scale parameters (MADs)
    median_absolute_deviation_duration = []
    median_absolute_deviation_coverage = []
    median_absolute_deviation_above_existing = []
    median_absolute_deviation_above_initial = []
    median_absolute_deviation_percent_coverage = []
    median_absolute_deviation_percent_above_existing = []
    median_absolute_deviation_percent_above_initial = []

    # loop over a given buffer radius
    for multiplier in multipliers:

        br = initial_buffer * multiplier

        multiplier_duration = []
        multiplier_coverage_total = []
        multiplier_coverage_above_existing = []
        multiplier_coverage_above_initial = []
        multiplier_percent_coverage = []
        multiplier_percent_above_existing = []
        multiplier_percent_above_initial = []

        # loop a pre-specified number of times for the given multiplier
        for _ in range(iterations_per_multiplier):

            # generating a random initial set of candidate satellites
            initial_candidate_gdf = generate_new_satellites(
                EPSG=equal_area_epsg,
                num_satellites=30,
                input_map=input_map,
                true_random=False,
                sample_separate=False,
            )

            # concatenate the GeoDataFrames and copy it
            sat_gdf = pd.concat([existing_satellite_gdf, initial_candidate_gdf])
            # initial_sat_gdf = sat_gdf.copy()

            # get the initial land coverage
            initial_land_coverage = calculate_land_coverage(
                gdf=sat_gdf, map=input_map, buffer_radius=br
            )

            # Using the regular hill climbing
            best_gdf, best_land_coverage, iteration_time_list, _ = hill_climbing(
                satellite_gdf=sat_gdf,
                map=input_map,
                new_satellite_column="new_satellite",
                buffer_radius=br,
                equal_distance_epsg=equal_distance_epsg,
                equal_area_epsg=equal_area_epsg,
                perturbation_distance=perturb_dist,
                print_progress=print_progress,
                num_iterations=ls_iter,
            )
            # get the time for the episode and append it to the multiplier list
            episode_time = np.sum(iteration_time_list)
            multiplier_duration.append(episode_time)

            # append the total coverage to the multiplier list
            multiplier_coverage_total.append(best_land_coverage)
            percent_total_land_coverage = 100 * (
                best_land_coverage / earth_land_area_square_km
            )
            multiplier_percent_coverage.append(percent_total_land_coverage)

            # get the total added land coverage above existing for the episode
            new_satellites = best_gdf[best_gdf["new_satellite"]]
            added_land_coverage = calculate_land_coverage(
                gdf=new_satellites, map=input_map, buffer_radius=br
            )
            percent_land_coverage_above_existing = 100 * (
                added_land_coverage / earth_land_area_square_km
            )
            multiplier_coverage_above_existing.append(added_land_coverage)
            multiplier_percent_above_existing.append(
                percent_land_coverage_above_existing
            )

            # get the total added land coverage above the initial configuration
            multiplier_coverage_above_initial.append(
                best_land_coverage - initial_land_coverage
            )
            percent_land_coverage_above_initial = 100 * (
                (best_land_coverage - initial_land_coverage) / earth_land_area_square_km
            )
            multiplier_percent_above_initial.append(percent_land_coverage_above_initial)

        # once the iterations for the multiplier have completed,
        # get medians and MADs
        average_durations.append(np.median(multiplier_duration))
        average_total_coverage.append(np.median(multiplier_coverage_total))
        average_added_coverage_above_existing.append(
            np.median(multiplier_coverage_above_existing)
        )
        average_added_coverage_above_initial.append(
            np.median(multiplier_coverage_above_initial)
        )
        average_percent_coverage.append(np.median(multiplier_percent_coverage))
        average_percent_above_existing.append(
            np.median(multiplier_percent_above_existing)
        )
        average_percent_above_initial.append(
            np.median(multiplier_percent_above_initial)
        )

        median_absolute_deviation_duration.append(
            median_abs_deviation(multiplier_duration)
        )
        median_absolute_deviation_coverage.append(
            median_abs_deviation(multiplier_coverage_total)
        )
        median_absolute_deviation_above_existing.append(
            median_abs_deviation(multiplier_coverage_above_existing)
        )
        median_absolute_deviation_above_initial.append(
            median_abs_deviation(multiplier_coverage_above_initial)
        )
        median_absolute_deviation_percent_coverage.append(
            median_abs_deviation(multiplier_percent_coverage)
        )
        median_absolute_deviation_percent_above_existing.append(
            median_abs_deviation(multiplier_percent_above_existing)
        )
        median_absolute_deviation_percent_above_initial.append(
            median_abs_deviation(multiplier_percent_above_initial)
        )

    return (
        average_durations,
        median_absolute_deviation_duration,
        average_total_coverage,
        median_absolute_deviation_coverage,
        average_added_coverage_above_existing,
        median_absolute_deviation_above_existing,
        average_added_coverage_above_initial,
        median_absolute_deviation_above_initial,
        average_percent_coverage,
        median_absolute_deviation_percent_coverage,
        average_percent_above_existing,
        median_absolute_deviation_percent_above_existing,
        average_percent_above_initial,
        median_absolute_deviation_percent_above_initial,
    )


def buffer_variation_rrsa(
    existing_satellite_gdf,
    equal_area_epsg,
    equal_distance_epsg,
    input_map,
    perturb_dist,
    ls_iter=10,
    init_temp=1000,
    cr=0.95,
    restart_thresh=10,
    initial_buffer=12065,
    multipliers=[0.5, 1, 2, 10, 20],
    iterations_per_multiplier=3,
):

    earth_land_area_square_km = 148326000

    # initializing lists for the recording of parameters

    # location parameters (medians)
    average_durations = []
    average_total_coverage = []
    average_added_coverage_above_existing = []
    average_added_coverage_above_initial = []
    average_percent_coverage = []
    average_percent_above_existing = []
    average_percent_above_initial = []

    # scale parameters (MADs)
    median_absolute_deviation_duration = []
    median_absolute_deviation_coverage = []
    median_absolute_deviation_above_existing = []
    median_absolute_deviation_above_initial = []
    median_absolute_deviation_percent_coverage = []
    median_absolute_deviation_percent_above_existing = []
    median_absolute_deviation_percent_above_initial = []

    # loop over a given buffer radius
    for multiplier in multipliers:

        br = initial_buffer * multiplier

        multiplier_duration = []
        multiplier_coverage_total = []
        multiplier_coverage_above_existing = []
        multiplier_coverage_above_initial = []
        multiplier_percent_coverage = []
        multiplier_percent_above_existing = []
        multiplier_percent_above_initial = []

        # loop a pre-specified number of times for the given multiplier
        for _ in range(iterations_per_multiplier):

            # generating a random initial set of candidate satellites
            initial_candidate_gdf = generate_new_satellites(
                EPSG=equal_area_epsg,
                num_satellites=30,
                input_map=input_map,
                true_random=False,
                sample_separate=False,
            )

            # concatenate the GeoDataFrames and copy it
            sat_gdf = pd.concat([existing_satellite_gdf, initial_candidate_gdf])
            # initial_sat_gdf = sat_gdf.copy()

            # get the initial land coverage
            initial_land_coverage = calculate_land_coverage(
                gdf=sat_gdf, map=input_map, buffer_radius=br
            )

            # Using the regular hill climbing
            best_gdf, best_land_coverage, iteration_time_list, _ = (
                random_restart_simulated_annealing(
                    satellite_gdf=sat_gdf,
                    map=input_map,
                    new_satellite_column="new_satellite",
                    buffer_radius=br,
                    equal_distance_epsg=equal_distance_epsg,
                    equal_area_epsg=equal_area_epsg,
                    perturbation_distance=perturb_dist,
                    num_iterations=ls_iter,
                    initial_temp=init_temp,
                    cooling_rate=cr,
                    restart_threshold=restart_thresh,
                )
            )
            # get the time for the episode and append it to the multiplier list
            episode_time = np.sum(iteration_time_list)
            multiplier_duration.append(episode_time)

            # append the total coverage to the multiplier list
            multiplier_coverage_total.append(best_land_coverage)
            percent_total_land_coverage = 100 * (
                best_land_coverage / earth_land_area_square_km
            )
            multiplier_percent_coverage.append(percent_total_land_coverage)

            # get the total added land coverage above existing for the episode
            new_satellites = best_gdf[best_gdf["new_satellite"]]
            added_land_coverage = calculate_land_coverage(
                gdf=new_satellites, map=input_map, buffer_radius=br
            )
            percent_land_coverage_above_existing = 100 * (
                added_land_coverage / earth_land_area_square_km
            )
            multiplier_coverage_above_existing.append(added_land_coverage)
            multiplier_percent_above_existing.append(
                percent_land_coverage_above_existing
            )

            # get the total added land coverage above the initial configuration
            multiplier_coverage_above_initial.append(
                best_land_coverage - initial_land_coverage
            )
            percent_land_coverage_above_initial = 100 * (
                (best_land_coverage - initial_land_coverage) / earth_land_area_square_km
            )
            multiplier_percent_above_initial.append(percent_land_coverage_above_initial)

        # once the iterations for the multiplier have completed,
        # get medians and MADs
        average_durations.append(np.median(multiplier_duration))
        average_total_coverage.append(np.median(multiplier_coverage_total))
        average_added_coverage_above_existing.append(
            np.median(multiplier_coverage_above_existing)
        )
        average_added_coverage_above_initial.append(
            np.median(multiplier_coverage_above_initial)
        )
        average_percent_coverage.append(np.median(multiplier_percent_coverage))
        average_percent_above_existing.append(
            np.median(multiplier_percent_above_existing)
        )
        average_percent_above_initial.append(
            np.median(multiplier_percent_above_initial)
        )

        median_absolute_deviation_duration.append(
            median_abs_deviation(multiplier_duration)
        )
        median_absolute_deviation_coverage.append(
            median_abs_deviation(multiplier_coverage_total)
        )
        median_absolute_deviation_above_existing.append(
            median_abs_deviation(multiplier_coverage_above_existing)
        )
        median_absolute_deviation_above_initial.append(
            median_abs_deviation(multiplier_coverage_above_initial)
        )
        median_absolute_deviation_percent_coverage.append(
            median_abs_deviation(multiplier_percent_coverage)
        )
        median_absolute_deviation_percent_above_existing.append(
            median_abs_deviation(multiplier_percent_above_existing)
        )
        median_absolute_deviation_percent_above_initial.append(
            median_abs_deviation(multiplier_percent_above_initial)
        )

    return (
        average_durations,
        median_absolute_deviation_duration,
        average_total_coverage,
        median_absolute_deviation_coverage,
        average_added_coverage_above_existing,
        median_absolute_deviation_above_existing,
        average_added_coverage_above_initial,
        median_absolute_deviation_above_initial,
        average_percent_coverage,
        median_absolute_deviation_percent_coverage,
        average_percent_above_existing,
        median_absolute_deviation_percent_above_existing,
        average_percent_above_initial,
        median_absolute_deviation_percent_above_initial,
    )


# running the loops for each algorithm
(
    hc_average_durations,
    hc_median_absolute_deviation_duration,
    hc_average_total_coverage,
    hc_median_absolute_deviation_coverage,
    hc_average_added_coverage_above_existing,
    hc_median_absolute_deviation_above_existing,
    hc_average_added_coverage_above_initial,
    hc_median_absolute_deviation_above_initial,
    hc_average_percent_coverage,
    hc_median_absolute_deviation_percent_coverage,
    hc_average_percent_above_existing,
    hc_median_absolute_deviation_percent_above_existing,
    hc_average_percent_above_initial,
    hc_median_absolute_deviation_percent_above_initial,
) = buffer_variation_hc(
    existing_satellite_gdf=base_gdf,
    equal_area_epsg=EQUAL_AREA_EPSG,
    equal_distance_epsg=EQUAL_DISTANCE_EPSG,
    input_map=land_map,
    perturb_dist=PERTURB_DISTANCE_KM,
)

(
    rrsa_average_durations,
    rrsa_median_absolute_deviation_duration,
    rrsa_average_total_coverage,
    rrsa_median_absolute_deviation_coverage,
    rrsa_average_added_coverage_above_existing,
    rrsa_median_absolute_deviation_above_existing,
    rrsa_average_added_coverage_above_initial,
    rrsa_median_absolute_deviation_above_initial,
    rrsa_average_percent_coverage,
    rrsa_median_absolute_deviation_percent_coverage,
    rrsa_average_percent_above_existing,
    rrsa_median_absolute_deviation_percent_above_existing,
    rrsa_average_percent_above_initial,
    rrsa_median_absolute_deviation_percent_above_initial,
) = buffer_variation_rrsa(
    existing_satellite_gdf=base_gdf,
    equal_area_epsg=EQUAL_AREA_EPSG,
    equal_distance_epsg=EQUAL_DISTANCE_EPSG,
    input_map=land_map,
    perturb_dist=PERTURB_DISTANCE_KM,
)


# plotting function for this parameter sensitivity analysis
def plot_for_sensitivity_analysis(
    multipliers,
    initial_buffer,
    hc_parameter,
    hc_error,
    rrsa_parameter,
    rrsa_error,
    y_label,
    title,
):

    # create the figure
    fig, ax = plt.subplots()

    multiplier_np = np.asarray(multipliers)
    radii = initial_buffer * multiplier_np

    hc_line = ax.errorbar(
        x=radii,
        y=hc_parameter,
        yerr=hc_error,
        color="orange",
        label="First-Choice Hill Climbing",
    )
    rrsa_line = ax.errorbar(
        x=radii,
        y=rrsa_parameter,
        yerr=rrsa_error,
        color="blue",
        label="Random Restart Hill Climbing +\nSimulated Annealing",
    )

    ax.legend(handles=[hc_line, rrsa_line], fontsize="x-small")
    ax.set_xlabel("Coverage Radius, Meters")
    ax.set_ylabel(y_label)
    plt.title(title)

    return fig


default_multipliers = [0.5, 1, 2, 10, 20]
default_buffer = 12065

durations_versus_buffer = plot_for_sensitivity_analysis(
    multipliers=default_multipliers,
    initial_buffer=default_buffer,
    hc_parameter=hc_average_durations,
    hc_error=hc_median_absolute_deviation_duration,
    rrsa_parameter=rrsa_average_durations,
    rrsa_error=rrsa_median_absolute_deviation_duration,
    y_label="Duration,\nSeconds",
    title="Runtime Duration as a function of Coverage Radius",
)
total_coverage_versus_buffer = plot_for_sensitivity_analysis(
    multipliers=default_multipliers,
    initial_buffer=default_buffer,
    hc_parameter=hc_average_total_coverage,
    hc_error=hc_median_absolute_deviation_coverage,
    rrsa_parameter=rrsa_average_total_coverage,
    rrsa_error=rrsa_median_absolute_deviation_coverage,
    y_label="Total Coverage Area,\nSquare Km",
    title="Total Satellite Coverage Area\nas a Function of Coverage Radius",
)
coverage_above_existing_versus_buffer = plot_for_sensitivity_analysis(
    multipliers=default_multipliers,
    initial_buffer=default_buffer,
    hc_parameter=hc_average_added_coverage_above_existing,
    hc_error=hc_median_absolute_deviation_above_existing,
    rrsa_parameter=rrsa_average_added_coverage_above_existing,
    rrsa_error=rrsa_median_absolute_deviation_above_existing,
    y_label="Added Coverage Area\nSquare Km",
    title="Coverage Area Added Above Existing Satellite Constellation\nas a Function of Coverage Radius",
)
coverage_above_initial_versus_buffer = plot_for_sensitivity_analysis(
    multipliers=default_multipliers,
    initial_buffer=default_buffer,
    hc_parameter=hc_average_added_coverage_above_initial,
    hc_error=hc_median_absolute_deviation_above_initial,
    rrsa_parameter=rrsa_average_added_coverage_above_initial,
    rrsa_error=rrsa_median_absolute_deviation_above_initial,
    y_label="Added Coverage Area\nSquare Km",
    title="Coverage Area Added Above Initial Random Placement\nas a Function of Coverage Radius",
)
percent_coverage_versus_buffer = plot_for_sensitivity_analysis(
    multipliers=default_multipliers,
    initial_buffer=default_buffer,
    hc_parameter=hc_average_percent_coverage,
    hc_error=hc_median_absolute_deviation_percent_coverage,
    rrsa_parameter=rrsa_average_percent_coverage,
    rrsa_error=rrsa_median_absolute_deviation_percent_coverage,
    y_label="Percent Land Coverage",
    title="Total Percent Land Area Covered\nas a Function of Coverage Radius",
)
percent_above_existing_versus_buffer = plot_for_sensitivity_analysis(
    multipliers=default_multipliers,
    initial_buffer=default_buffer,
    hc_parameter=hc_average_percent_above_existing,
    hc_error=hc_median_absolute_deviation_percent_above_existing,
    rrsa_parameter=rrsa_average_percent_above_existing,
    rrsa_error=rrsa_median_absolute_deviation_percent_above_existing,
    y_label="Percent Added Land Coverage",
    title="Percent Land Area Added Above Existing Satellites\nas a Function of Coverage Radius",
)
percent_above_initial_versus_buffer = plot_for_sensitivity_analysis(
    multipliers=default_multipliers,
    initial_buffer=default_buffer,
    hc_parameter=hc_average_percent_above_initial,
    hc_error=hc_median_absolute_deviation_percent_above_initial,
    rrsa_parameter=rrsa_average_percent_above_initial,
    rrsa_error=rrsa_median_absolute_deviation_percent_above_initial,
    y_label="Percent Added Land Coverage",
    title="Percent Land Area Added Above Initial Configuration\nas a Function of Coverage Radius",
)
