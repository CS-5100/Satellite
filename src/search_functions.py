import geodataframe_processing as gdfp
import pandas as pd
import numpy as np
import time
import math


def hill_climbing(
    satellite_gdf,
    map,
    new_satellite_column,
    buffer_radius,
    equal_distance_epsg,
    equal_area_epsg,
    perturbation_distance,
    num_iterations=10,
    print_progress=False,
    random_seed=None,
):
    # create random number generator
    rng = np.random.default_rng(seed=random_seed)

    best_gdf = satellite_gdf.copy()
    best_land_coverage = gdfp.calculate_land_coverage(
        gdf=best_gdf, map=map, buffer_radius=buffer_radius
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
            gdf=new_gdf, map=map, buffer_radius=buffer_radius
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
    satellite_gdf,
    map,
    new_satellite_column,
    buffer_radius,
    equal_distance_epsg,
    equal_area_epsg,
    perturbation_distance,
    num_iterations=100,
    initial_temp=1000,
    cooling_rate=0.95,
    restart_threshold=10,
    random_seed=None,
):
    # create random number generator
    rng = np.random.default_rng(seed=random_seed)

    # copy initial satellite geodataframe
    # into two separate geodataframes
    best_gdf = satellite_gdf.copy()
    current_gdf = best_gdf.copy()

    # initialize the best land coverage section
    # and current land coverage section
    best_land_coverage = gdfp.calculate_land_coverage(
        gdf=best_gdf, map=map, buffer_radius=buffer_radius
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
            gdf=new_gdf, map=map, buffer_radius=buffer_radius
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
                input_map=map,
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
