"""
Exp-3. For Step Index Fibers:
a) Calculate the Acceptance angle(θa).
b) Calculate the waveguide dispersion at given wavelength and plot of
waveguide dispersion vs wavelength(λ).
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Parameters
n1 = 1.46  # Core refractive index
n2 = 1.45  # Cladding refractive index
a = 5e-6   # Core radius (m)
c = 3e8    # Speed of light (m/s)

# Function to calculate Numerical Aperture (NA)
def calculate_na(n1, n2):
    return np.sqrt(n1**2 - n2**2)

# Function to calculate Acceptance angle (θa) in degrees
def calculate_acceptance_angle(na):
    return np.degrees(np.arcsin(na))

# Function to calculate V number
def calculate_v_number(a, na, wavelength):
    return (2 * np.pi * a * na) / wavelength

# Function to calculate normalized propagation constant (b)
def calculate_b(v):
    return 1 - (1.1428 / v)**2

# Function to calculate effective refractive index
def calculate_n_eff(n2, b, n1):
    return np.sqrt(n2**2 + b * (n1**2 - n2**2))

# Function to calculate propagation constant (β)
def calculate_beta(k0, n_eff):
    return k0 * n_eff

# Function to calculate waveguide dispersion (D_w) in ps/(nm·km)
def calculate_waveguide_dispersion(lambda_range, beta_range, c):
    d_lambda = lambda_range[1] - lambda_range[0]
    d_beta = np.gradient(beta_range, d_lambda)  # First derivative
    d2_beta = np.gradient(d_beta, d_lambda)    # Second derivative
    lambda_range_disp = lambda_range[5:-5]     # Adjust for length reduction
    d2_beta = d2_beta[5:-5]                    # Trim to match
    return - (lambda_range_disp / c) * d2_beta * 1e6, lambda_range_disp

# Calculate NA
NA = calculate_na(n1, n2)

# Calculate Acceptance angle
theta_a_deg = calculate_acceptance_angle(NA)

# Wavelength range (500 nm to 2000 nm)
lambda_range = np.linspace(500e-9, 2000e-9, 100)

# Calculate V number, b, n_eff, and beta for the range
V_range = calculate_v_number(a, NA, lambda_range)
b_range = calculate_b(V_range)
n_eff_range = calculate_n_eff(n2, b_range, n1)
k0_range = 2 * np.pi / lambda_range
beta_range = calculate_beta(k0_range, n_eff_range)

# Calculate waveguide dispersion
D_w_range, lambda_range_disp = calculate_waveguide_dispersion(lambda_range, beta_range, c)

# Calculate V number at a reference wavelength (e.g., 650 nm for consistency)
V_0 = calculate_v_number(a, NA, 650e-9)

# Print results
print(f"Numerical Aperture (NA): {NA:.4f}")
print(f"Acceptance angle (θa): {theta_a_deg:.4f} degrees")
print(f"V number at 650 nm: {V_0:.4f}")
# Note: Waveguide dispersion is plotted over the range, no single value printed

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
sns.lineplot(x=lambda_range_disp*1e6, y=D_w_range, ax=ax, color='blue')
ax.set_xlabel("Wavelength (μm)")
ax.set_ylabel("Waveguide Dispersion D_w (ps/nm/km)")
ax.set_title("Waveguide Dispersion vs Wavelength (Step-Index Fiber)")
ax.grid(True)
plt.tight_layout()
plt.show()