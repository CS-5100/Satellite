import numpy as np

# fitness function
def coverage_gain(new_position):
    return calculate_unique_coverage_area(new_position)

def simulated_annealing_for_placement(new_sat_position, temperature, cooling_rate, max_iterations):
    
    # Initialize the position for the new satellite with the initial guess
    current_position = initial_position
    current_gain = coverage_gain(new_sat_position)
    
    best_position = current_position
    best_gain = current_gain

    for iteration in range(max_iterations):
        lat, lon, alt = current_position
        new_position = (
            lat + np.random.uniform(-1, 1),  
            lon + np.random.uniform(-1, 1),
            alt + np.random.uniform(-5, 5)
        )
        
        new_gain = coverage_gain(new_position)
        
        gain_diff = new_gain - current_gain
        if gain_diff > 0 or np.random.rand() < np.exp(gain_diff / temperature):
            current_position = new_position
            current_gain = new_gain
            
            if new_gain > best_gain:
                best_position = new_position
                best_gain = new_gain

        temperature *= cooling_rate

        if iteration % 100 == 0:
            print(f"Iteration {iteration}, Best Gain: {best_gain:.2f}, Temperature: {temperature:.2f}")
    
    return best_position

# Initial guess for the new satellite position
initial_position = (500, 25, 20)

# Simulated annealing parameters
temperature = 1000  
cooling_rate = 0.995  
max_iterations = 1000  

# Find the best position for the new satellite
optimal_position = simulated_annealing_for_placement(initial_position, temperature, cooling_rate, max_iterations)
print("Optimal position for the new satellite:", optimal_position)
