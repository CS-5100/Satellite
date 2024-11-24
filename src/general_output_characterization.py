import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import time
from scipy.stats import mannwhitneyu, median_abs_deviation
from datetime import datetime
from geodataframe_processing import (
    load_satellites_from_file,
    load_land_ocean_data,
    generate_new_satellites,
    calculate_land_coverage,
)
from search_functions import hill_climbing, random_restart_simulated_annealing
procedure_start = time.time()
# EPSG Codes
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Time the TLE data was collected
TLE_TIME = datetime(2024, 11, 14, 16, 12, 32, 202983)

# Other Global Parameters
PERTURB_DISTANCE_KM = 100
LOCAL_SEARCH_ITERATIONS = 200
NUMBER_EPISODES = 50


# Importing a static set of TLEs for sensitivity analysis
base_gdf = load_satellites_from_file(
    EPSG=EQUAL_AREA_EPSG,
    input_filename="parameter_sensitivity_analysis_tles.txt",
    time=TLE_TIME,
)

# Loading land coverage area
land_map, ocean_map = load_land_ocean_data(EPSG=EQUAL_AREA_EPSG)


# defining methods for looping across many episodes and gathering performance data
def generate_episodes_hc(
    existing_satellite_gdf,
    equal_area_epsg,
    equal_distance_epsg,
    input_map,
    perturb_dist,
    num_episodes,
    ls_iter=10,
    buffer=12065,
    print_progress=True,
):

    earth_land_area_square_km = 148326000

    # initializing lists for the recording of parameters

    # location parameters (medians)
    durations = []
    total_coverage = []
    added_coverage_above_existing = []
    added_coverage_above_initial = []
    percent_coverage = []
    percent_above_existing = []
    percent_above_initial = []
    coverage_lists = []

    # loop over a set number of episodes
    for episode in range(num_episodes):

        print(
            "First-Choice Hill Climbing Iteration:",
            (episode + 1),
            "out of",
            num_episodes,
        )

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

        # get the initial land coverage
        initial_land_coverage = calculate_land_coverage(
            gdf=sat_gdf, map=input_map, buffer_radius=buffer
        )

        # Using the regular hill climbing
        best_gdf, best_land_coverage, iteration_time_list, iteration_coverage_list = (
            hill_climbing(
                satellite_gdf=sat_gdf,
                map=input_map,
                new_satellite_column="new_satellite",
                buffer_radius=buffer,
                equal_distance_epsg=equal_distance_epsg,
                equal_area_epsg=equal_area_epsg,
                perturbation_distance=perturb_dist,
                print_progress=print_progress,
                num_iterations=ls_iter,
            )
        )

        # get the time for the episode in minutes and append it to the multiplier list
        episode_time = np.sum(iteration_time_list) / 60.0
        durations.append(episode_time)

        # append the total coverage to the multiplier list
        total_coverage.append(best_land_coverage)
        percent_total_land_coverage = 100 * (
            best_land_coverage / earth_land_area_square_km
        )
        percent_coverage.append(percent_total_land_coverage)

        # get the total added land coverage above existing for the episode
        new_satellites = best_gdf[best_gdf["new_satellite"]]
        added_land_coverage = calculate_land_coverage(
            gdf=new_satellites, map=input_map, buffer_radius=buffer
        )

        percent_land_coverage_above_existing = 100 * (
            added_land_coverage / earth_land_area_square_km
        )

        added_coverage_above_existing.append(added_land_coverage)
        percent_above_existing.append(percent_land_coverage_above_existing)

        # get the total added land coverage above the initial configuration
        added_coverage_above_initial.append(best_land_coverage - initial_land_coverage)
        percent_land_coverage_above_initial = 100 * (
            (best_land_coverage - initial_land_coverage) / earth_land_area_square_km
        )
        percent_above_initial.append(percent_land_coverage_above_initial)

        coverage_lists.append(iteration_coverage_list)

    return (
        durations,
        total_coverage,
        added_coverage_above_existing,
        added_coverage_above_initial,
        percent_coverage,
        percent_above_existing,
        percent_above_initial,
        coverage_lists,
    )


def generate_episodes_rrsa(
    existing_satellite_gdf,
    equal_area_epsg,
    equal_distance_epsg,
    input_map,
    perturb_dist,
    num_episodes,
    ls_iter=10,
    init_temp=1000,
    cr=0.95,
    restart_thresh=10,
    buffer=12065,
):

    earth_land_area_square_km = 148326000

    # initializing lists for the recording of parameters

    # initializing lists for the recording of parameters

    # location parameters (medians)
    durations = []
    total_coverage = []
    added_coverage_above_existing = []
    added_coverage_above_initial = []
    percent_coverage = []
    percent_above_existing = []
    percent_above_initial = []
    coverage_lists = []

    # loop over a set number of episodes
    for episode in range(num_episodes):

        print(
            "Random Restart Hill Climbing With Simulated Annealing:",
            (episode + 1),
            "out of",
            num_episodes,
        )

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

        # get the initial land coverage
        initial_land_coverage = calculate_land_coverage(
            gdf=sat_gdf, map=input_map, buffer_radius=buffer
        )

        # Using the regular hill climbing
        best_gdf, best_land_coverage, iteration_time_list, iteration_coverage_list = (
            random_restart_simulated_annealing(
                satellite_gdf=sat_gdf,
                map=input_map,
                new_satellite_column="new_satellite",
                buffer_radius=buffer,
                equal_distance_epsg=equal_distance_epsg,
                equal_area_epsg=equal_area_epsg,
                perturbation_distance=perturb_dist,
                num_iterations=ls_iter,
                initial_temp=init_temp,
                cooling_rate=cr,
                restart_threshold=restart_thresh,
            )
        )

        # get the time for the episode in minutes and append it to the multiplier list
        episode_time = np.sum(iteration_time_list) / 60.0
        durations.append(episode_time)

        # append the total coverage to the multiplier list
        total_coverage.append(best_land_coverage)
        percent_total_land_coverage = 100 * (
            best_land_coverage / earth_land_area_square_km
        )
        percent_coverage.append(percent_total_land_coverage)

        # get the total added land coverage above existing for the episode
        new_satellites = best_gdf[best_gdf["new_satellite"]]
        added_land_coverage = calculate_land_coverage(
            gdf=new_satellites, map=input_map, buffer_radius=buffer
        )

        percent_land_coverage_above_existing = 100 * (
            added_land_coverage / earth_land_area_square_km
        )

        added_coverage_above_existing.append(added_land_coverage)
        percent_above_existing.append(percent_land_coverage_above_existing)

        # get the total added land coverage above the initial configuration
        added_coverage_above_initial.append(best_land_coverage - initial_land_coverage)
        percent_land_coverage_above_initial = 100 * (
            (best_land_coverage - initial_land_coverage) / earth_land_area_square_km
        )
        percent_above_initial.append(percent_land_coverage_above_initial)

        coverage_lists.append(iteration_coverage_list)

    return (
        durations,
        total_coverage,
        added_coverage_above_existing,
        added_coverage_above_initial,
        percent_coverage,
        percent_above_existing,
        percent_above_initial,
        coverage_lists,
    )


# running the loops for each algorithm
(
    hc_durations,
    hc_total_coverage,
    hc_added_coverage_above_existing,
    hc_added_coverage_above_initial,
    hc_percent_coverage,
    hc_percent_above_existing,
    hc_percent_above_initial,
    hc_coverage_lists,
) = generate_episodes_hc(
    existing_satellite_gdf=base_gdf,
    equal_area_epsg=EQUAL_AREA_EPSG,
    equal_distance_epsg=EQUAL_DISTANCE_EPSG,
    input_map=land_map,
    perturb_dist=PERTURB_DISTANCE_KM,
    num_episodes=NUMBER_EPISODES,
    ls_iter=LOCAL_SEARCH_ITERATIONS,
)

(
    rrsa_durations,
    rrsa_total_coverage,
    rrsa_added_coverage_above_existing,
    rrsa_added_coverage_above_initial,
    rrsa_percent_coverage,
    rrsa_percent_above_existing,
    rrsa_percent_above_initial,
    rrsa_coverage_lists,
) = generate_episodes_rrsa(
    existing_satellite_gdf=base_gdf,
    equal_area_epsg=EQUAL_AREA_EPSG,
    equal_distance_epsg=EQUAL_DISTANCE_EPSG,
    input_map=land_map,
    perturb_dist=PERTURB_DISTANCE_KM,
    num_episodes=NUMBER_EPISODES,
    ls_iter=LOCAL_SEARCH_ITERATIONS,
)


hc_dict = {
    "Duration": hc_durations,
    "Total_Coverage": hc_total_coverage,
    "Coverage_Above_Existing": hc_added_coverage_above_existing,
    "Coverage_Above_Initial": hc_added_coverage_above_initial,
    "Percent_Coverage": hc_percent_coverage,
    "Percent_Coverage_Above_Existing": hc_percent_above_existing,
    "Percent_Coverage_Above_Initial": hc_percent_above_initial,
}

rrsa_dict = {
    "Duration": rrsa_durations,
    "Total_Coverage": rrsa_total_coverage,
    "Coverage_Above_Existing": rrsa_added_coverage_above_existing,
    "Coverage_Above_Initial": rrsa_added_coverage_above_initial,
    "Percent_Coverage": rrsa_percent_coverage,
    "Percent_Coverage_Above_Existing": rrsa_percent_above_existing,
    "Percent_Coverage_Above_Initial": rrsa_percent_above_initial,
}

hc_df = pd.DataFrame(hc_dict)
hc_df["Algorithm"] = "First-Choice Hill Climbing"

rrsa_df = pd.DataFrame(rrsa_dict)
rrsa_df["Algorithm"] = "Random Restart Hill Climbing with Simulated Annealing"

hc_episode_df = pd.DataFrame(data=hc_coverage_lists)
rrsa_episode_df = pd.DataFrame(data=rrsa_coverage_lists)

df = pd.concat([hc_df, rrsa_df])
summary_filename = (
    "Distribution_"
    + str(NUMBER_EPISODES)
    + "episodes_"
    + str(LOCAL_SEARCH_ITERATIONS)
    + "iterations.csv"
)

df.to_csv(summary_filename)
hc_episode_df.to_csv("hc_episodes.csv")
rrsa_episode_df.to_csv("rrsa_episodes.csv")

# df_2 = pd.read_csv("Distribution_50episodes_200iterations.csv")
# hc_df_2 = df_2[df_2["Algorithm"] == "First-Choice Hill Climbing"]
# rrsa_df_2 = df_2[
#     df_2["Algorithm"] == "Random Restart Hill Climbing with Simulated Annealing"
# ]

# hc_durations = hc_df_2.iloc[:, 1]
# hc_total_coverage = hc_df_2.iloc[:, 2]
# hc_added_coverage_above_existing = hc_df_2.iloc[:, 3]
# hc_added_coverage_above_initial = hc_df_2.iloc[:, 4]
# hc_percent_coverage = hc_df_2.iloc[:, 5]
# hc_percent_above_existing = hc_df_2.iloc[:, 6]
# hc_percent_above_initial = hc_df_2.iloc[:, 7]

# rrsa_durations = rrsa_df_2.iloc[:, 1]
# rrsa_total_coverage = rrsa_df_2.iloc[:, 2]
# rrsa_added_coverage_above_existing = rrsa_df_2.iloc[:, 3]
# rrsa_added_coverage_above_initial = rrsa_df_2.iloc[:, 4]
# rrsa_percent_coverage = rrsa_df_2.iloc[:, 5]
# rrsa_percent_above_existing = rrsa_df_2.iloc[:, 6]
# rrsa_percent_above_initial = rrsa_df_2.iloc[:, 7]

# print(
#     "Median, First-Choice Hill Climbing: {:.2e} +- {:.2e}".format(
#         np.median(hc_added_coverage_above_initial), median_abs_deviation(hc_added_coverage_above_initial)
#     )
# )
# print(
#     "Median, Random Restart Hill Climbing with Simulated Annealing: {:.2e} +- {:.2e}".format(
#         np.median(rrsa_added_coverage_above_initial), median_abs_deviation(rrsa_added_coverage_above_initial)
#     )
# )

hc_episodes_df = pd.read_csv("hc_episodes.csv", index_col=0)
rrsa_episodes_df = pd.read_csv("rrsa_episodes.csv", index_col=0)

hc_episodes_col_average = hc_episodes_df.apply(np.mean, axis=0)
hc_episodes_col_mad = hc_episodes_df.apply(median_abs_deviation, axis=0)

rrsa_episodes_col_average = rrsa_episodes_df.apply(np.mean, axis=0)
rrsa_episodes_col_mad = rrsa_episodes_df.apply(median_abs_deviation, axis=0)

hc_coverage_lists = hc_episodes_df.values.tolist()
rrsa_coverage_lists = rrsa_episodes_df.values.tolist()

procedure_end = time.time()
procedure_duration_seconds = procedure_end - procedure_start
procedure_duration_minutes = procedure_duration_seconds / 60.0
procedure_duration_hours = procedure_duration_minutes / 60.0

print("The procedure took", procedure_duration_hours, "hours")

def plot_episode_lists(
    hc_lists,
    rrsa_lists,
    iterations,
    title,
):

    # create the figure
    fig, ax = plt.subplots()

    x_values = [i for i in range(1, iterations + 1)]

    for y_list in hc_lists:
        ax.plot(x_values, y_list, color="orange", alpha = 0.5)

    for y_list in rrsa_lists:
        ax.plot(x_values, y_list, color="blue", alpha = 0.5)

    orange_line_handle = mlines.Line2D(
        [], [], color="orange", label="First-Choice Hill Climbing"
    )
    blue_line_handle = mlines.Line2D(
        [],
        [],
        color="blue",
        label="Random Restart Hill Climbing +\nSimulated Annealing",
    )
    ax.legend(
        handles=[orange_line_handle, blue_line_handle],
        fontsize="x-small",
        loc="lower right",
    )
    ax.set_xlabel("Local Search Iterations")
    ax.set_ylabel("Coverage Area, Square Kilometers")
    plt.title(title)

    return fig


def plot_episode_average(
    hc_median,
    hc_mad,
    rrsa_median,
    rrsa_mad,
    iterations
):

    # create the figure
    fig, ax = plt.subplots()

    x_values = [i for i in range(1, iterations + 1)]

    hc_line = ax.errorbar(
        x=x_values,
        y=hc_median,
        yerr=hc_mad,
        color="orange",
        label="First-Choice Hill Climbing",
        fmt='o-',
        alpha = 0.5
    )

    rrsa_line = ax.errorbar(
        x=x_values,
        y=rrsa_median,
        yerr=rrsa_mad,
        color="blue",
        label="Random Restart Hill Climbing +\nSimulated Annealing",
        fmt='o-',
        alpha = 0.5
    )

    ax.legend(handles=[hc_line, rrsa_line], fontsize="x-small")
    ax.set_ylabel("Coverage Radius, Kilometers")
    ax.set_xlabel("Local Search Iterations")
    plt.title("Episode Coverage Area Medians")

    return fig


def violin_plot_compare(
    hc_parameter,
    rrsa_parameter,
    y_label,
    title,
):

    # create the figure
    fig, ax = plt.subplots()

    colors = ["orange", "blue"]

    _, p = mannwhitneyu(hc_parameter, rrsa_parameter)
    p_value = "Two Sided p Value: {:.2e}".format(p)

    plot_data = [hc_parameter, rrsa_parameter]

    hc_x_var = [1 for _ in range(len(hc_parameter))]
    rrsa_x_var = [2 for _ in range(len(rrsa_parameter))]
    x_var = hc_x_var + rrsa_x_var

    quartile1, _, quartile3 = np.percentile(plot_data, [25, 50, 75], axis=1)

    ax.set_xticks(
        [y + 1 for y in range(len(plot_data))],
        labels=[
            "First-Choice Hill Climbing",
            "Random Restart Hill Climbing +\nSimulated Annealing",
        ],
    )
    
    ax.set_ylabel(y_label)
    ax.text(
        0.5,
        0.05,
        p_value,
        transform=ax.transAxes,
        horizontalalignment="center",
        verticalalignment="bottom",
    )

    ax.set_title(title)
    ax.scatter(x_var, plot_data, marker="o", color="black", zorder=2.5)
    ax.vlines(np.array([1, 2]), quartile1, quartile3, color=colors, linestyle="-", lw=5)

    violin_plot_parts = ax.violinplot(plot_data, showmedians=True, showextrema=False)

    for i, pc in enumerate(violin_plot_parts["bodies"]):
        pc.set_color(colors[i])

    return fig


episode_iterations = plot_episode_lists(
    hc_lists=hc_coverage_lists,
    rrsa_lists=rrsa_coverage_lists,
    iterations=LOCAL_SEARCH_ITERATIONS,
    title="Coverage Area Versus Iteration",
)

episode_averages = plot_episode_average(
    hc_median=hc_episodes_col_average,
    hc_mad=hc_episodes_col_mad,
    rrsa_median=rrsa_episodes_col_average,
    rrsa_mad=rrsa_episodes_col_mad,
    iterations=LOCAL_SEARCH_ITERATIONS
)

duration_distributions = violin_plot_compare(
    hc_parameter=hc_durations,
    rrsa_parameter=rrsa_durations,
    y_label="Duration, Minutes",
    title="Empirical Distribution of Durations",
)
coverage_distributions = violin_plot_compare(
    hc_parameter=hc_total_coverage,
    rrsa_parameter=rrsa_total_coverage,
    y_label="Total Coverage Area\nSquare Kilometers",
    title="Empirical Distribution of\nTotal Coverage Areas",
)
coverage_above_existing_distributions = violin_plot_compare(
    hc_parameter=hc_added_coverage_above_existing,
    rrsa_parameter=rrsa_added_coverage_above_existing,
    y_label="Added Coverage Area\nSquare Kilometers",
    title="Empirical Distribution of\nAdded Coverage Areas Above Existing",
)
coverage_above_initial_distributions = violin_plot_compare(
    hc_parameter=hc_added_coverage_above_initial,
    rrsa_parameter=rrsa_added_coverage_above_initial,
    y_label="Added Coverage Area\nSquare Kilometers",
    title="Empirical Distribution of\nAdded Coverage Areas Above Initial Configuration",
)
percent_coverage_distributions = violin_plot_compare(
    hc_parameter=hc_percent_coverage,
    rrsa_parameter=rrsa_percent_coverage,
    y_label="Percent Earth Land Area, Total",
    title="Empirical Distribution of\nTotal Earth Coverage Percent",
)
percent_coverage_above_existing_distributions = violin_plot_compare(
    hc_parameter=hc_percent_above_existing,
    rrsa_parameter=rrsa_percent_above_existing,
    y_label="Percent Earth Land Area, Added",
    title="Empirical Distribution of\nAdded Earth Coverage Percent Above Existing",
)
percent_coverage_above_initial_distributions = violin_plot_compare(
    hc_parameter=hc_percent_above_initial,
    rrsa_parameter=rrsa_percent_above_initial,
    y_label="Percent Earth Land Area, Added",
    title="Empirical Distribution of\nAdded Earth Coverage Percent Above Initial Configuration",
)


episode_iterations.savefig("episode_iterations.png")
episode_averages.savefig("episode_medians.png")
duration_distributions.savefig("duration_distributions.png")
coverage_distributions.savefig("coverage_distributions.png")
coverage_above_existing_distributions.savefig(
    "coverage_above_existing_distributions.png"
)
coverage_above_initial_distributions.savefig("coverage_above_initial_distributions.png")
percent_coverage_distributions.savefig("percent_coverage_distributions.png")
percent_coverage_above_existing_distributions.savefig(
    "percent_coverage_above_existing_distributions.png"
)
percent_coverage_above_initial_distributions.savefig(
    "percent_coverage_above_initial_distributions.png"
)

plt.show()
