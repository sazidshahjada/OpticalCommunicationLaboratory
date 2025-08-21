"""
Exp-6. For Pure Silica Fibers:
a. Observe the graphical representation of refractive index(n) vs wavelength(λ).
b. Observe the graphical representation of material dispersion vs wavelength(λ).
"""

import numpy as np
import matplotlib.pyplot as plt

# Sellmeier coefficients for pure silica
A1, A2, A3 = 0.6961663, 0.4079426, 0.8974794
B1, B2, B3 = 0.0684043, 0.1162414, 9.896161

# Wavelength range (μm) based on the figure
lambda_um = np.linspace(0.6, 2.0, 1000)

# Refractive index using Sellmeier equation
def refractive_index(lambda_um):
    lambda2 = lambda_um**2
    n2 = 1 + (A1 * lambda2 / (lambda2 - B1**2)) + (A2 * lambda2 / (lambda2 - B2**2)) + (A3 * lambda2 / (lambda2 - B3**2))
    return np.sqrt(n2)

# Material dispersion: D_m = -lambda/c * d^2n/d(lambda)^2 (ps/(nm·km))
def material_dispersion(lambda_um, n):
    c = 3e8  # Speed of light (m/s)
    # Convert lambda to meters for derivatives
    lambda_m = lambda_um * 1e-6
    # Numerical second derivative of n w.r.t. lambda
    d_lambda = lambda_m[1] - lambda_m[0]
    dn_dlambda = np.gradient(n, d_lambda)
    d2n_dlambda2 = np.gradient(dn_dlambda, d_lambda)
    # Material dispersion in ps/(nm·km)
    D_m = - (lambda_m / c) * d2n_dlambda2 * 1e12 / 1e-9  # Convert to ps/(nm·km)
    return D_m[5:-5], lambda_um[5:-5]

# Compute refractive index
n = refractive_index(lambda_um)

# Compute material dispersion
D_m, lambda_um_disp = material_dispersion(lambda_um, n)

# Create 1x2 subplot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Subplot 1: Refractive Index vs. Wavelength
ax1.plot(lambda_um, n, color='blue')
ax1.set_xlabel('Wavelength (μm)')
ax1.set_ylabel('Refractive Index (n)')
ax1.set_title('Refractive Index vs. Wavelength for Pure Silica')
ax1.grid(True)

# Subplot 2: Material Dispersion vs. Wavelength
ax2.plot(lambda_um_disp, D_m, color='red')
ax2.set_xlabel('Wavelength (μm)')
ax2.set_ylabel('Material Dispersion (ps/(nm·km))')
ax2.set_title('Material Dispersion vs. Wavelength for Pure Silica')
ax2.axhline(0, color='black', linestyle='--', linewidth=0.5)  # Zero dispersion line
ax2.grid(True)

# Adjust layout
plt.tight_layout()
plt.show()