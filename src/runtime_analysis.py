import matplotlib.pyplot as plt
from simulated_annealing_v1 import *

num_iterations = 50  # Fixed number of iterations for comparison
additional_satellite_counts = [15, 30, 45, 60]

# Storage for average runtimes
average_runtimes = []

# Run 10 trials for each count of new satellites
for satellite_count in additional_satellite_counts:
    print(f"Testing with {satellite_count} additional satellites...")

    # Temporary list to store runtimes for the current satellite count
    runtimes = []

    for trial in range(10):
        # Copy the satellite DataFrame and modify to include only the required count of new satellites
        modified_satellite_gdf = satellite_gdf.copy()
        modified_satellite_gdf['new_satellite'] = False
        modified_satellite_gdf.loc[-satellite_count:, 'new_satellite'] = True  # Mark last `satellite_count` as new

        # Measure runtime
        start_time = time.time()
        _, _, _, _ = simulated_annealing(
            modified_satellite_gdf, land_map, 'new_satellite', num_iterations=num_iterations,
            initial_temp=1000, cooling_rate=0.95
        )
        end_time = time.time()

        # Store runtime for this trial
        runtime = end_time - start_time
        runtimes.append(runtime)
        print(f"Trial {trial + 1}: Runtime for {satellite_count} satellites: {runtime:.2f} seconds")

    # Calculate average runtime for the current satellite count
    avg_runtime = sum(runtimes) / len(runtimes)
    average_runtimes.append(avg_runtime)
    print(f"Average runtime for {satellite_count} satellites: {avg_runtime:.2f} seconds\n")

# Plotting the average runtime as a function of additional satellites
plt.figure(figsize=(10, 6))
plt.plot(additional_satellite_counts, average_runtimes, marker='o', color='b')
plt.xlabel("Number of Additional Satellites")
plt.ylabel("Average Runtime (seconds)")
plt.title("Average Runtime as a Function of Number of Additional Satellites")
plt.grid(True)
plt.show()