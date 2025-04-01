import numpy as np
from scipy.special import eval_gegenbauer
from scipy.linalg import lstsq

def sample_sphere(s, num_samples):
    """
    Generate `num_samples` points uniformly distributed on the unit s-sphere S^s.
    Uses the Marsaglia method.
    """
    X = np.random.randn(num_samples, s + 1)  # Gaussian samples in R^(s+1)
    X /= np.linalg.norm(X, axis=1, keepdims=True)  # Normalize to project onto S^s
    return X

def hyperspherical_harmonics(s, L_max, X):
    """
    Compute hyperspherical harmonics up to degree L_max on S^s.
    
    Parameters:
    - s: Dimension of the sphere S^s.
    - L_max: Maximum degree of harmonics.
    - X: Points on S^s (shape: [num_samples, s+1])
    
    Returns:
    - Harmonic matrix H of shape [num_samples, num_harmonics]
    """
    num_samples = X.shape[0]
    harmonics = []

    # Compute harmonic basis up to L_max
    for L in range(L_max + 1):
        for m in range(-L, L + 1):  # Fourier components
            # Compute Gegenbauer polynomials in higher dimensions
            C_L = eval_gegenbauer(L, s / 2 - 1, X[:, 0])  # Generalized Legendre poly
            
            # Azimuthal dependence: Use sin/cos for m ≠ 0
            if m == 0:
                harmonics.append(C_L)
            else:
                harmonics.append(C_L * np.cos(m * np.arctan2(X[:, 2], X[:, 1])))
                harmonics.append(C_L * np.sin(m * np.arctan2(X[:, 2], X[:, 1])))

    return np.column_stack(harmonics)  # Design matrix

# Parameters
s = 3  # Dimension of the sphere S^s
num_samples = 200  # Number of observations
L_max = 20  # Max harmonic degree

# Generate data on S^s
X = sample_sphere(s, num_samples)

# Define a function on S^s (for example, a simple harmonic-like function)
f_true = np.sin(2 * np.pi * X[:, 1]) + 0.1 * np.random.randn(num_samples)

# Construct the hyperspherical harmonic basis
H = hyperspherical_harmonics(s, L_max, X)

# Solve least squares problem H * c = f
coeffs, _, _, _ = lstsq(H, f_true)

# Function to approximate f at new points
def f_approx(X_new, coeffs, s, L_max):
    H_new = hyperspherical_harmonics(s, L_max, X_new)
    return H_new @ coeffs

# Example: Approximate at new points
X_test = sample_sphere(s, 50)
f_reconstructed = f_approx(X_test, coeffs, s, L_max)

# Print first few results
print("True f values (first 5):", f_true[:5])
print("Reconstructed f values (first 5):", f_reconstructed[:5])

