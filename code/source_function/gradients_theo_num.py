import jax.numpy as jnp

def pairwise_distances(points):
    """Compute pairwise Euclidean distances between points."""
    diff = points[:, None, :] - points[None, :, :]
    distances = jnp.linalg.norm(diff, axis=-1)
    return distances

def energy(points):
    """Calculate the energy of the system as the sum of inverse distances."""
    distances = pairwise_distances(points)
    distances = jnp.maximum(distances, 1e-8)  # Prevent division by zero
    return jnp.sum(1.0 / distances[jnp.triu_indices_from(distances, k=1)])  # Sum upper triangle without diagonal

def analytical_gradient(points):
    """Compute the gradient of the energy function analytically."""
    n, k = points.shape
    distances = pairwise_distances(points)
    
    # Prevent division by zero by clamping distances
    distances = jnp.maximum(distances, 1e-8)
    
    # Initialize gradients
    grads = jnp.zeros_like(points)
    
    # Vectorized gradient calculation
    for i in range(n):
        for j in range(n):
            if i != j:
                r_ij = distances[i, j]
                grad_ij = - (1 / r_ij**3) * (points[i] - points[j])
                grads = grads.at[i].add(grad_ij)  # Accumulate the gradient for point i
    
    return grads

def numerical_gradient(points, epsilon=1e-6):
    """Compute the gradient of the energy function using the finite difference method (vectorized)."""
    n, k = points.shape
    
    # Ensure the points and epsilon are in float64 precision
    points = points.astype(jnp.float64)
    epsilon = jnp.float64(epsilon)
    
    grads = jnp.zeros_like(points, dtype=jnp.float64)
    
    # Perturb each point and calculate energy
    perturbed_points_plus = points + epsilon
    perturbed_points_minus = points - epsilon
    
    # Vectorized energy computation
    energy_plus = energy(perturbed_points_plus)
    energy_minus = energy(perturbed_points_minus)
    
    # Compute numerical gradient using central difference
    grads = (energy_plus - energy_minus) / (2 * epsilon)
    
    return grads

def compare_gradients(points):
    """Compare the numerical and analytical gradients."""
    # Ensure points are in float64 precision
    points = points.astype(jnp.float64)
    
    analytical_grads = analytical_gradient(points)
    numerical_grads = numerical_gradient(points)
    
    print("Analytical Gradients:\n", analytical_grads)
    print("\nNumerical Gradients:\n", numerical_grads)
    
    # Compute the difference between the gradients
    difference = jnp.linalg.norm(analytical_grads - numerical_grads)
    print("\nGradient Difference (L2 norm):", difference)

# Example usage:
n = 5  # Number of points
k = 2  # Dimension of the sphere (2D sphere)
points = jnp.array([[1.0, 0.0], [0.5, 0.5], [-0.5, 0.5], [-1.0, 0.0], [0.0, -1.0]])  # Example points

# Compare the gradients
compare_gradients(points)

