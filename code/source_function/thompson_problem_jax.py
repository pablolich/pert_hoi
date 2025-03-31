import jax.numpy as jnp
import jax
import optax


def initialize_points(k, n, key):
    """Generate n points uniformly distributed on the unit k-sphere."""
    key, subkey = jax.random.split(key)
    points = jax.random.normal(subkey, (n, k + 1))
    points /= jnp.linalg.norm(points, axis=1, keepdims=True)
    return points

def energy(points):
    """Calculate the energy of the system as the sum of inverse distances."""
    distances = pairwise_distances(points)
    distances = distances[jnp.triu_indices_from(distances, k=1)]
    return jnp.sum(1.0 / distances)

def pairwise_distances(points):
    """Compute pairwise Euclidean distances between points."""
    diff = points[:, None, :] - points[None, :, :]
    distances = jnp.linalg.norm(diff, axis=-1)
    # Extract upper triangle without diagonal
    return distances 

def analytical_gradient(points):
    """Compute the gradient of the energy function with respect to each point."""
    n, k = points.shape
    
    # Compute pairwise distances
    distances = pairwise_distances(points)
    
    # Ensure no division by zero by clamping distances to a minimum value
    
    # Initialize gradients (same shape as points)
    grads = jnp.zeros_like(points)
    
    # Calculate the gradient of the energy function with respect to each point
    for i in range(n):
        for j in range(n):
            if i > j:
                r_ij = distances[i, j]
                grad_ij = - (1 / r_ij**3) * (points[i] - points[j])
                grads = grads.at[i].add(grad_ij)  # Accumulate the gradient for point i
    
    return grads

def numerical_gradient(points, epsilon=1e-6):
    """Compute the gradient of the energy function using the finite difference method."""
    n, k = points.shape
    points = points.astype(jnp.float64)
    grads = jnp.zeros_like(points, dtype = jnp.float64)

    for i in range(n):
        for d in range(k):  # Iterate over the dimensions of the points
            # Perturb the i-th point in the d-th direction
            perturbed_points = points.at[i, d].add(epsilon)
            energy_plus = energy(perturbed_points)
            
            perturbed_points = points.at[i, d].subtract(epsilon)
            energy_minus = energy(perturbed_points)
            
            # Numerical gradient using the central difference formula
            grads = grads.at[i, d].set((energy_plus - energy_minus) / (2 * epsilon))
    
    return grads

# Example usage:
n = 2  # Number of points
k = 2  # Dimension of the sphere (2D sphere)
points = jnp.array([[0.0, 1], [0.0, -1.0]])

import ipdb; ipdb.set_trace(context = 20)

grads = analytical_gradient(points)
num_grads = numerical_gradient(points)
print("Analytical gradients:\n", grads)
print("Numerical gradients:\n", num_grads)

import ipdb; ipdb.set_trace(context = 20)




def optimize_points(k, n, steps=1000, lr=0.01):
    """Optimize the positions of n points on the k-sphere to minimize energy."""
    key = jax.random.PRNGKey(42)
    points = initialize_points(k, n, key)
    
    # Define optimizer
    optimizer = optax.adam(lr)
    opt_state = optimizer.init(points)
    
    @jax.jit
    def step(points, opt_state):
        loss, grads = jax.value_and_grad(energy)(points)
        
        # Project gradients onto the tangent space of the sphere
        grads_tangent = grads - (jnp.sum(grads * points, axis=1, keepdims=True) * points)
        jax.debug.print("grads: {}", grads)
        updates, opt_state = optimizer.update(grads_tangent, opt_state)
        points = optax.apply_updates(points, updates)
        points /= jnp.linalg.norm(points, axis=1, keepdims=True)  # Project back to sphere
        
        return points, opt_state, loss
    
    for i in range(steps):
        points, opt_state, loss = step(points, opt_state)
        import ipdb; ipdb.set_trace(context = 20)
        if i % 10 == 0:
            print(f"Step {i}: Energy = {loss:.6f}")
    
    return points, loss


# Example usage
k = 2  # 2D sphere (3D space)
n = 10  # Number of points
optimized_points, final_energy = optimize_points(k, n)
print("Optimized points:", optimized_points)
print("Final energy:", final_energy)

