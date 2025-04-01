import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform

def generate_points_on_sphere(num_points, dim):
    """Generate random points uniformly distributed on the unit sphere."""
    points = np.random.randn(num_points, dim)  
    points /= np.linalg.norm(points, axis=1)[:, np.newaxis]  
    return points

def pairwise_vecs_dists(points):
    """Compute pairwise difference vectors and distances between all points."""
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    distances = np.linalg.norm(diff, axis=-1)
    normalized_diff = diff / distances[..., np.newaxis]  
    return normalized_diff, distances

def potential_energy(distances):
    """Compute total Coulomb potential energy (sum of 1/d for all pairs)."""
    return np.sum(1 / distances)

#def gradient_force(normalized_diff, distances, points):
#    """Compute net Coulomb force for each point, constrained to sphere."""
#    np.fill_diagonal(distances, np.nan)  
#    force_magnitudes = 1 / distances**2  
#    forces = normalized_diff * force_magnitudes[..., np.newaxis]  
#    net_forces = np.nansum(forces, axis=1)  
#
#    radial_unit_vectors = points  
#    radial_component = (np.sum(net_forces * radial_unit_vectors, axis=1)
#                        [:, np.newaxis] * radial_unit_vectors)
#
#    tangent_forces = net_forces - radial_component  
#    return tangent_forces

def objective_function(points, n, dim):
    """Objective function for optimization based on Coulomb potential energy."""
    points_mat = np.reshape(points, (n, dim))  
    pairwise_distances = pdist(points_mat)  
    return potential_energy(pairwise_distances)

def gradient_function(points, n, dim):
    """Compute the gradient (forces) for optimization."""
    points_mat = np.reshape(points, (n, dim))  
    normalized_vectors, pairwise_distances = pairwise_vecs_dists(points_mat)
    forces = gradient_force(normalized_vectors, pairwise_distances, points_mat)
    return forces.flatten()  

def constraint_unit_sphere(points, n, dim):#, lamb):
    """Ensure all points remain on the unit sphere."""
    points_mat = np.reshape(points, (n, dim))  
    return np.linalg.norm(points_mat, axis=1) - 1  

def nearest_neighbor_distance(positions):
    """Calculate the distance to the nearest neighbor for each point."""
    n = positions.shape[0]
    distances = np.zeros(n)

    for i in range(n):
        # Calculate distances from point i to all other points
        dists = np.linalg.norm(positions[i] - positions, axis=1)
        # Set distance to self to a large value to ignore it
        dists[i] = np.inf
        # Find the minimum distance to any other point
        distances[i] = np.min(dists)

    return distances

def plot_histogram(initial_points, optimized_points):
    """Plot histograms of pairwise distances before and after optimization."""
    initial_dists = nearest_neighbor_distance(initial_points)
    optimized_dists = nearest_neighbor_distance(optimized_points)
    plt.hist(optimized_dists, bins=30, alpha=0.5, label="Optimized", color="red")
    plt.hist(initial_dists, bins=30, alpha=0.5, label="Initial", color="blue")

    plt.xlabel("Pairwise Distance")
    plt.ylabel("Frequency")
    plt.title("Histogram of Pairwise Distances")
    plt.legend()
    plt.show()

def get_best_existing_cost(n, dim):
    """
    Check existing solution files and find the best (lowest) cost value.

    Parameters:
    n (int): Number of points.
    dim (int): Dimensionality.

    Returns:
    float: Lowest cost found in files, or infinity if no valid files exist.
    """
    best_cost = float('inf')
    pattern = re.compile(rf"optimized_n{n}_d{dim}_cost_([\d\.]+)\.txt")

    for filename in os.listdir():
        match = pattern.match(filename)
        if match:
            try:
                cost = float(match.group(1))
                if cost < best_cost:
                    best_cost = cost
            except ValueError:
                continue  

    return best_cost

def save_solution_if_best(optimized_points, cost, n, dim):
    """
    Save the solution only if it's the best found so far.

    Parameters:
    optimized_points (numpy array): Optimized point positions.
    cost (float): Optimization cost value.
    n (int): Number of points.
    dim (int): Dimensionality.
    """
    best_existing_cost = get_best_existing_cost(n, dim)

    if cost < best_existing_cost:
        filename = f"optimized_n{n}_d{dim}_cost_{cost:.6f}.txt"
        np.savetxt(filename, optimized_points, fmt="%.6f")
        print(f"New best solution saved: {filename}")
    else:
        print("No improvement, solution not saved.")
        
# Loop for n from 2 to 7, with 250 points each time
n = 25
for dim in range(2, 8):
    points = generate_points_on_sphere(n, dim)

    initial_guess = points.flatten()
    constraints = [{'type': 'eq', 'fun': constraint_unit_sphere, 'args': (n, dim)}]

    result = minimize(
        fun=objective_function,        
        x0=initial_guess,              
        args=(n, dim),                 
        constraints=constraints,       
        method='SLSQP',                
        options={'disp': True, 'maxiter': 1000}  
    )

    optimized_points = result.x.reshape((n, dim))
    save_solution_if_best(optimized_points, result.fun, n, dim)

    print("Distance to nearest neighbor")
    print(nearest_neighbor_distance(optimized_points))
    plot_histogram(points, optimized_points)


