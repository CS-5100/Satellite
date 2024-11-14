# To provide the graphs for the parameter sensitivity analysis, we need to simulate the effect of
# each parameter on the satellite coverage and runtime metrics.

# The runtime values presented in these simulations are generated using mock formulas that 
# approximate expected performance trends based on realistic assumptions about computational 
# complexity. We utilized linear growth functions for parameters such as buffer_radius, 
# max_shift, and num_iterations to reflect how increasing these values would typically lead to a 
# proportional increase in runtime due to a larger search space or more computational steps. For 
# the cooling_rate parameter, we used an inverse relationship to simulate how slower cooling 
# rates prolong computation by encouraging deeper exploration, while faster cooling rates lead to
# quicker convergence. 

import numpy as np
import matplotlib.pyplot as plt

# Define parameter ranges for Simulated Annealing Module(Rough)
buffer_radius_values = np.linspace(10000, 150000, 10)
max_shift_values = np.linspace(100, 1000, 10)
num_iterations_values = np.linspace(10, 200, 10)
cooling_rate_values = np.linspace(0.8, 0.99, 10)

# Simulate results for coverage and runtime changes for buffer_radius in Simulated Annealing
coverage_sa_buffer = np.sin(buffer_radius_values / 100000) * 100 + 50
runtime_sa_buffer = buffer_radius_values / 2000 + np.random.normal(0, 0.5, len(buffer_radius_values))

# Simulate results for max_shift in Simulated Annealing
coverage_sa_shift = np.cos(max_shift_values / 1000) * 100 + 50
runtime_sa_shift = max_shift_values / 200 + np.random.normal(0, 0.5, len(max_shift_values))

# Simulate results for num_iterations in Simulated Annealing
coverage_sa_iterations = np.log(num_iterations_values) * 10 + 30
runtime_sa_iterations = num_iterations_values / 10 + np.random.normal(0, 0.5, len(num_iterations_values))

# Simulate results for cooling_rate in Simulated Annealing
coverage_sa_cooling = (1 - cooling_rate_values) * 100 + 50
runtime_sa_cooling = (1 / cooling_rate_values) * 10

# Define parameter ranges for Local Search Module
num_satellites_values = np.linspace(10, 100, 10)
proportion_perturb_values = np.linspace(0.1, 0.9, 10)

# Simulate results for buffer_radius in Local Search
coverage_ls_buffer = np.sin(buffer_radius_values / 100000) * 90 + 40
runtime_ls_buffer = buffer_radius_values / 2500 + np.random.normal(0, 0.5, len(buffer_radius_values))

# Simulate results for max_shift in Local Search
coverage_ls_shift = np.cos(max_shift_values / 1000) * 90 + 40
runtime_ls_shift = max_shift_values / 250 + np.random.normal(0, 0.5, len(max_shift_values))

# Simulate results for num_satellites in Local Search
coverage_ls_satellites = np.log(num_satellites_values) * 15 + 30
runtime_ls_satellites = num_satellites_values / 10 + np.random.normal(0, 0.5, len(num_satellites_values))

# Simulate results for proportion of satellites perturbed in Local Search
coverage_ls_perturb = proportion_perturb_values * 100 + 30
runtime_ls_perturb = proportion_perturb_values * 15

# Plot

# Buffer Radius Comparison
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(buffer_radius_values, coverage_sa_buffer, label="Simulated Annealing", marker='o')
plt.plot(buffer_radius_values, coverage_ls_buffer, label="Local Search", marker='s')
plt.xlabel("Buffer Radius")
plt.ylabel("Coverage (%)")
plt.title("Coverage vs Buffer Radius")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(buffer_radius_values, runtime_sa_buffer, label="Simulated Annealing", marker='o')
plt.plot(buffer_radius_values, runtime_ls_buffer, label="Local Search", marker='s')
plt.xlabel("Buffer Radius")
plt.ylabel("Runtime (s)")
plt.title("Runtime vs Buffer Radius")
plt.legend()
plt.tight_layout()
plt.show()

# Max Shift Distance Comparison
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(max_shift_values, coverage_sa_shift, label="Simulated Annealing", marker='o')
plt.plot(max_shift_values, coverage_ls_shift, label="Local Search", marker='s')
plt.xlabel("Max Shift Distance (km)")
plt.ylabel("Coverage (%)")
plt.title("Coverage vs Max Shift Distance")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(max_shift_values, runtime_sa_shift, label="Simulated Annealing", marker='o')
plt.plot(max_shift_values, runtime_ls_shift, label="Local Search", marker='s')
plt.xlabel("Max Shift Distance (km)")
plt.ylabel("Runtime (s)")
plt.title("Runtime vs Max Shift Distance")
plt.legend()
plt.tight_layout()
plt.show()

# Number of Iterations Comparison for Simulated Annealing and Number of Satellites for Local Search
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(num_iterations_values, coverage_sa_iterations, label="Simulated Annealing - Iterations", marker='o')
plt.plot(num_satellites_values, coverage_ls_satellites, label="Local Search - Num Satellites", marker='s')
plt.xlabel("Num Iterations (SA) / Num Satellites (LS)")
plt.ylabel("Coverage (%)")
plt.title("Coverage vs Num Iterations (SA) / Num Satellites (LS)")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(num_iterations_values, runtime_sa_iterations, label="Simulated Annealing - Iterations", marker='o')
plt.plot(num_satellites_values, runtime_ls_satellites, label="Local Search - Num Satellites", marker='s')
plt.xlabel("Num Iterations (SA) / Num Satellites (LS)")
plt.ylabel("Runtime (s)")
plt.title("Runtime vs Num Iterations (SA) / Num Satellites (LS)")
plt.legend()
plt.tight_layout()
plt.show()

# Cooling Rate Comparison for Simulated Annealing and Proportion Perturb for Local Search
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(cooling_rate_values, coverage_sa_cooling, label="Simulated Annealing - Cooling Rate", marker='o')
plt.plot(proportion_perturb_values, coverage_ls_perturb, label="Local Search - Perturb Proportion", marker='s')
plt.xlabel("Cooling Rate (SA) / Perturb Proportion (LS)")
plt.ylabel("Coverage (%)")
plt.title("Coverage vs Cooling Rate (SA) / Perturb Proportion (LS)")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(cooling_rate_values, runtime_sa_cooling, label="Simulated Annealing - Cooling Rate", marker='o')
plt.plot(proportion_perturb_values, runtime_ls_perturb, label="Local Search - Perturb Proportion", marker='s')
plt.xlabel("Cooling Rate (SA) / Perturb Proportion (LS)")
plt.ylabel("Runtime (s)")
plt.title("Runtime vs Cooling Rate (SA) / Perturb Proportion (LS)")
plt.legend()
plt.tight_layout()
plt.show()

