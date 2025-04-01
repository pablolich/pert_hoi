import numpy as np
from scipy.special import sph_harm
from scipy.linalg import lstsq

# Generate sample data
np.random.seed(42)  # Reproducibility
num_samples = 100  # Number of points on S^2
l_max = 5  # Maximum degree of spherical harmonics

# Random points on S^2
theta = np.arccos(2 * np.random.rand(num_samples) - 1)  # Polar angle
phi = 2 * np.pi * np.random.rand(num_samples)  # Azimuthal angle

# Define a function on the sphere (e.g., f = Y_3^2 + noise)
f_true = sph_harm(2, 3, phi, theta).real + 0.1 * np.random.randn(num_samples)

# Construct spherical harmonic basis up to l_max
harmonics = []
for l in range(l_max + 1):
    for m in range(-l, l + 1):
        harmonics.append(sph_harm(m, l, phi, theta).real)

H = np.column_stack(harmonics)  # Design matrix

# Solve least squares problem H * c = f
coeffs, _, _, _ = lstsq(H, f_true)

# Define reconstruction function
def f_approx(theta, phi, coeffs, l_max):
    approx = 0
    idx = 0
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            approx += coeffs[idx] * sph_harm(m, l, phi, theta).real
            idx += 1
    return approx

# Compute reconstructed values for the same training points

f_reconstructed = f_approx(theta, phi, coeffs, l_max)

# Visualization (Optional)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.plot([-1, 1], [-1, 1], 'r--', label="Ideal Fit (y=x)") 
plt.scatter(f_true, f_reconstructed)
plt.show()
