import geodataframe_processing as gdfp
import pandas as pd
import numpy as np
import time
import math
from geopandas import GeoDataFrame


def hill_climbing(
    satellite_gdf: GeoDataFrame,
    input_map: GeoDataFrame,
    new_satellite_column: str,
    buffer_radius: float,
    equal_distance_epsg: (str | int),
    equal_area_epsg: (str | int),
    perturbation_distance: float,
    num_iterations=10,
    print_progress=False,
    random_seed=None,
):
    """Given a geopandas GeoDataFrame of existing satellites and new candidate satellite positions, and a
    GeoDataFrame defining a land area to cover, use First-Choice Hill Climbing to improve upon the random
    candidate satellite positions with respect to land area coverage. Current perturbation logic in modules may
    introduce candidate satellite positions that have empty or NaN coordinates, so beware of that.

    Args:
        satellite_gdf (GeoDataFrame): A set of satellite coordinates in GeoDataFrame format, including both a set of
            existing satellites within a satellite constellation (generated with either gdfp.load_existing_satellites or
            gdfp.load_existing_satellites) concatenated with a set of new candidate satellites generated with
            gdfp.generate_new_satellites.
        input_map (GeoDataFrame): The vector map data defining the geometric shape of the land mass to cover.
        new_satellite_column (str): The string denoting which GeoDataFrame column contains the boolean variable marking new satellites.
        buffer_radius (float): the assumed coverage radius of each satellite in meters
        equal_distance_epsg (str  |  int): an equal distance coordinate reference system specified by an EPSG code, used for accurate positional perturbation calculations
        equal_area_epsg (str  |  int): an equal area coordinate reference system specified by an EPSG code, used for accurate coverage area intersection calculations
        perturbation_distance (float): The maximum shift in latitude and longitude a satellite position can be perturbed by, in kilometers
        num_iterations (int, optional): The number of local search iterations to execute. Defaults to 10.
        print_progress (bool, optional): Whether or not to print the program's progress towards a more optimal solution to the standard output. Defaults to False.
        random_seed (int, optional): Whether or not to use a random seed for the random number generator. Not guaranteed to be used internally at this point of development.
            Defaults to None.

    Returns:
        (GeoDataFrame, float, list[float], list[float]): A tuple of values
            - best_gdf (GeoDataFrame): a GeoDataFrame containing the existing satellites within a satellite constellation in addition to new candidate satellite positions that
                that have a land coverage area equal to or greater than the initial random satellite positions in the original GeoDataFrame
            - best_land_coverage (float): the total land area covered by all the satellites in the output GeoDataFrame
            - iteration_time_list (list[float]): the time that each local search iteration took in seconds as a list. Mainly for debugging/trace purposes.
            - best_coverage_list (float): the land area covered by the best set of existing candidate satellite positions so far at each
                local search iteration as a list. Mainly for debugging/trace purposes.
            
    """
    # create random number generator
    rng = np.random.default_rng(seed=random_seed)

    best_gdf = satellite_gdf.copy()
    best_land_coverage = gdfp.calculate_land_coverage(
        gdf=best_gdf, map=input_map, buffer_radius=buffer_radius
    )
    
    iteration_time_list = []
    best_coverage_list = []

    for iteration in range(num_iterations):
        iteration_start = time.time()
        if print_progress:
            print(f"Iteration {iteration + 1}/{num_iterations}")

        # Perturb positions of new satellites only
        new_gdf = best_gdf.copy()
        gdfp.perturb_positions(
            gdf=new_gdf,
            new_satellite_column=new_satellite_column,
            eq_dist_epsg=equal_distance_epsg,
            eq_area_epsg=equal_area_epsg,
            max_shift_km=perturbation_distance,
            random_state=rng.integers(low=0, high=10001),
        )

        # Calculate new land coverage area
        new_land_coverage = gdfp.calculate_land_coverage(
            gdf=new_gdf, map=input_map, buffer_radius=buffer_radius
        )

        if new_land_coverage > best_land_coverage:
            if print_progress:
                print(
                    f"Improvement found: {new_land_coverage:.2f} km² vs {best_land_coverage:.2f} km²"
                )
            best_gdf, best_land_coverage = new_gdf, new_land_coverage
        else:
            if print_progress:
                print("No improvement in land coverage.")
        iteration_end = time.time()
        
        iteration_duration = iteration_end - iteration_start
        
        best_coverage_list.append(best_land_coverage)
        iteration_time_list.append(iteration_duration)
                
    if print_progress:
        print(f"Optimized land coverage area: {best_land_coverage:.2f} km²")
    return best_gdf, best_land_coverage, iteration_time_list, best_coverage_list


def random_restart_simulated_annealing(
    satellite_gdf: GeoDataFrame,
    input_map: GeoDataFrame,
    new_satellite_column: str,
    buffer_radius: float,
    equal_distance_epsg: (str | int),
    equal_area_epsg: (str | int),
    perturbation_distance: float,
    num_iterations=100,
    initial_temp=1000,
    cooling_rate=0.95,
    restart_threshold=10,
    random_seed=None,
):
    """Given a geopandas GeoDataFrame of existing satellites and new candidate satellite positions, and a
    GeoDataFrame defining a land area to cover, use First-Choice Hill Climbing in conjunction with random restarting and
    simulated annealing to improve upon the random candidate satellite positions with respect to land area coverage.
    Current perturbation logic in modules may introduce candidate satellite positions that have empty or NaN coordinates,
    so beware of that.

    Args:
        satellite_gdf (GeoDataFrame): A set of satellite coordinates in GeoDataFrame format, including both a set of
            existing satellites within a satellite constellation (generated with either gdfp.load_existing_satellites or
            gdfp.load_existing_satellites) concatenated with a set of new candidate satellites generated with
            gdfp.generate_new_satellites.
        input_map (GeoDataFrame): The vector map data defining the geometric shape of the land mass to cover.
        new_satellite_column (str): The string denoting which GeoDataFrame column contains the boolean variable marking new satellites.
        buffer_radius (float): the assumed coverage radius of each satellite in meters
        equal_distance_epsg (str  |  int): an equal distance coordinate reference system specified by an EPSG code, used for accurate positional perturbation calculations
        equal_area_epsg (str  |  int): an equal area coordinate reference system specified by an EPSG code, used for accurate coverage area intersection calculations
        perturbation_distance (float): The maximum shift in latitude and longitude a satellite position can be perturbed by, in kilometers
        num_iterations (int, optional): The number of local search iterations to execute. Defaults to 100.
        initial_temp (int, optional): The initial temperature to use for the simulated annealing acceptance probability calculations. Defaults to 1000.
        cooling_rate (float, optional): The float to multiply the temperature by after each local search iteration. Must be between zero and one for proper functionality.
            Defaults to 0.95.
        restart_threshold (int, optional): The threshold for the number of local search iterations to execute without land coverage area improvement
            before randomly restarting with new candidate satellite positions and resetting the temperature parameter. Defaults to 10.
        random_seed (int, optional): Whether or not to use a random seed for the random number generator. Not guaranteed to be used internally at this point of development.
            Defaults to None.

    Returns:
        (GeoDataFrame, float, list[float], list[float]): A tuple of values
            - best_gdf (GeoDataFrame): a GeoDataFrame containing the existing satellites within a satellite constellation in addition to new candidate satellite positions that
                that have a land coverage area equal to or greater than the initial random satellite positions in the original GeoDataFrame
            - best_land_coverage (float): the total land area covered by all the satellites in the output GeoDataFrame
            - iteration_times (list[float]): the time that each local search iteration took in seconds as a list. Mainly for debugging/trace purposes.
            - best_land_coverage_list (float): the land area covered by the best set of existing candidate satellite positions so far at each
                local search iteration as a list. Mainly for debugging/trace purposes.
    """
    # create random number generator
    rng = np.random.default_rng(seed=random_seed)

    # copy initial satellite geodataframe
    # into two separate geodataframes
    best_gdf = satellite_gdf.copy()
    current_gdf = best_gdf.copy()

    # initialize the best land coverage section
    # and current land coverage section
    best_land_coverage = gdfp.calculate_land_coverage(
        gdf=best_gdf, map=input_map, buffer_radius=buffer_radius
    )
    current_land_coverage = best_land_coverage

    # initialize temperature variables
    temperature = initial_temp

    # Tracking iteration times and land coverage for plotting
    iteration_times = []
    best_land_coverage_list = []

    # Initialize restart counter
    stagnation_counter = 0

    for iteration in range(num_iterations):
        start_time = time.time()  # Start timing the iteration
        print(
            f"Iteration {iteration + 1}/{num_iterations}, Temperature: {temperature:.2f}"
        )

        # Perturb positions of new satellites only
        new_gdf = current_gdf.copy()
        gdfp.perturb_positions(
            new_gdf,
            new_satellite_column,
            eq_dist_epsg=equal_distance_epsg,
            eq_area_epsg=equal_area_epsg,
            max_shift_km=perturbation_distance,
            random_state=rng.integers(low=0, high=10001),
        )

        # Calculate new land coverage area
        new_land_coverage = gdfp.calculate_land_coverage(
            gdf=new_gdf, map=input_map, buffer_radius=buffer_radius
        )

        # Calculate the improvement or deterioration
        delta_coverage = new_land_coverage - current_land_coverage

        # Acceptance probability for worse solutions
        acceptance_probability = (
            math.exp(delta_coverage / temperature) if delta_coverage < 0 else 1
        )

        # if we find a configuration that improves the land coverage
        # or the temperature allows us to accept a worse solution
        # re-set the current geodataframe and coverage
        if delta_coverage > 0 or rng.random() < acceptance_probability:
            print(
                f"Accepted new configuration with land coverage: {new_land_coverage:.2f} km²"
            )
            current_gdf, current_land_coverage = new_gdf, new_land_coverage
            stagnation_counter = 0  # Reset stagnation counter

            # Update best solution if this is the best encountered so far
            if new_land_coverage > best_land_coverage:
                best_gdf, best_land_coverage = new_gdf, new_land_coverage
        else:
            print("Retained current configuration.")

        # Cooling step
        temperature *= cooling_rate

        # Increment stagnation counter if no improvement
        stagnation_counter += 1

        # Check if we should restart due to stagnation
        if stagnation_counter >= restart_threshold:
            print("Random restart due to stagnation.")

            # get the current existing satellites from the input dataframe
            current_existing_satellite_gdf = satellite_gdf[
                satellite_gdf[new_satellite_column] == False
            ].copy()

            # generate a new set of candidate satellites
            new_candidate_satellite_gdf = gdfp.generate_new_satellites(
                EPSG=equal_area_epsg,
                num_satellites=30,
                input_map=input_map,
                true_random=False,
                sample_separate=False,
            )

            # concatenate both geodataframes to replace the current geodataframe
            current_gdf = pd.concat(
                [current_existing_satellite_gdf, new_candidate_satellite_gdf]
            )

            # get the land coverage
            current_land_coverage = gdfp.calculate_land_coverage(
                current_gdf, map, buffer_radius
            )

            # reset the stagnation counter annd temperature
            stagnation_counter = 0  # Reset the stagnation counter
            temperature = initial_temp  # Reset temperature for new start

        # Track iteration times and best land coverage
        elapsed_time = time.time() - start_time
        iteration_times.append(elapsed_time)
        best_land_coverage_list.append(best_land_coverage)

    print(
        f"Optimized land coverage area after simulated annealing: {best_land_coverage:.2f} km²"
    )
    return best_gdf, best_land_coverage, iteration_times, best_land_coverage_list
