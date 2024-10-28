import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def potential_energy(flat_positions, n, k):
    """Calculate the potential energy of the given positions."""
    positions = flat_positions.reshape(n, k)
    energy = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            distance = np.linalg.norm(positions[i] - positions[j])
            if distance > 0:
                energy += 1 / distance
    return energy

def jac(flat_postitions, n, k):
    positions = flat_positions.reshape(n, k)
    return gradient_vector 

def constraint(flat_positions, n, k):
    """Ensure points lie on the surface of a k-dimensional sphere."""
    positions = flat_positions.reshape(n, k)
    return np.sum(positions**2, axis=1) - 1  # Each point should be on the sphere

def thompson_problem(n, k, initial_positions):
    """Solve the Thompson problem for n points on a k-dimensional sphere."""

    # Constraints for optimization (points should lie on the sphere)
    con = {'type': 'eq', 'fun': lambda x: constraint(x, n, k)}

    # Optimize the positions using SLSQP method
    result = minimize(potential_energy, 
                      initial_positions.flatten(), 
                      args=(n, k), 
                      method='SLSQP', 
                      constraints=con,
                      options={'disp':True})

    # Reshape the result back into a 2D array of positions
    optimized_positions = result.x.reshape(n, k)
    # Normalize to ensure they lie on the sphere
    optimized_positions /= np.linalg.norm(optimized_positions, axis=1)[:, np.newaxis]
    
    return optimized_positions, result.fun


def plot_results(positions1, positions2):
    """Plot two sets of positions based on the number of dimensions."""
    k1 = positions1.shape[1]
    k2 = positions2.shape[1]

    # Create a figure with two subplots
    fig = plt.figure(figsize=(12, 6))

    # Plot for the first set of positions
    ax1 = fig.add_subplot(1, 2, 1)
    if k1 == 2:
        ax1.scatter(positions1[:, 0], positions1[:, 1], c='blue', marker='o')
        ax1.set_title("Set 1: Initial Points on 2D Sphere")
        ax1.set_xlim(-1, 1)
        ax1.set_ylim(-1, 1)
    elif k1 == 3:
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax1.scatter(positions1[:, 0], positions1[:, 1], positions1[:, 2], c='blue', marker='o')
        ax1.set_title("Set 1: Initial Points on 3D Sphere")
        ax1.set_xlim([-1, 1])
        ax1.set_ylim([-1, 1])
        ax1.set_zlim([-1, 1])

    ax1.grid()

    # Plot for the second set of positions
    ax2 = fig.add_subplot(1, 2, 2)
    if k2 == 2:
        ax2.scatter(positions2[:, 0], positions2[:, 1], c='orange', marker='o')
        ax2.set_title("Set 2: Optimized Points on 2D Sphere")
        ax2.set_xlim(-1, 1)
        ax2.set_ylim(-1, 1)
    elif k2 == 3:
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        ax2.scatter(positions2[:, 0], positions2[:, 1], positions2[:, 2], c='orange', marker='o')
        ax2.set_title("Set 2: Optimized Points on 3D Sphere")
        ax2.set_xlim([-1, 1])
        ax2.set_ylim([-1, 1])
        ax2.set_zlim([-1, 1])

    ax2.grid()

    plt.tight_layout()
    plt.show()

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

# Example usage
n = 10  # Number of points
k = 5   # Number of dimensions

# Sample initial positions
initial_positions = np.random.randn(n, k)
initial_positions /= np.linalg.norm(initial_positions, axis=1)[:, np.newaxis]
initial_nn_distances = nearest_neighbor_distance(initial_positions)
print(potential_energy(initial_positions, n, k))

positions, energy = thompson_problem(n, k, initial_positions)

print("Optimized Positions:\n", positions)
print("Minimum Potential Energy:", energy)

# Check the nearest neighbor distances
nn_distances = nearest_neighbor_distance(positions)
print("Nearest Neighbor Distances:", nn_distances)

# Create the histogram
fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_subplot(1, 2, 1)
ax1.hist(initial_nn_distances, bins=30, alpha=0.5, label='Vector 1', color='blue')
ax2 = fig.add_subplot(1, 2, 2)
ax2.hist(nn_distances, bins=30, alpha=0.5, label='Vector 2', color='orange')
plt.show()

# Plot the results if k <= 3
if k <= 3:
    plot_results(initial_positions, positions)

