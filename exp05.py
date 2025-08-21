"""
Exp-5. For Graded Index Fibers:
a. Observe the graphical representation of acceptance angle(θa) w.r.t
radius of core for different profile parameters.
b. Observe the graphical representation of core refractive index(n1) w.r.t. radius
of core for different profile parameters.
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters for graded-index fiber
n1 = 1.5  # Core refractive index at center
n2 = 1.48  # Cladding refractive index
delta = (n1**2 - n2**2) / (2 * n1**2)  # Relative refractive index difference
alpha_values = [1, 2, 4]  # Profile parameters
r_a = np.linspace(0, 1, 100)  # Normalized radius (r/a) from 0 to 1

# Function to calculate refractive index n(r) for a given alpha
def refractive_index(r_a, alpha, n1, delta):
    return n1 * np.sqrt(1 - 2 * delta * (r_a ** alpha))

# Function to calculate acceptance angle theta_a (in degrees)
def acceptance_angle(r_a, alpha, n1, n2, delta):
    n_r = refractive_index(r_a, alpha, n1, delta)
    NA = np.sqrt(n_r**2 - n2**2)  # Numerical aperture
    theta_a = np.arcsin(NA) * 180 / np.pi  # Convert to degrees
    return np.where(n_r > n2, theta_a, 0)  # Set to 0 if n(r) <= n2

# Create a 1x2 subplot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

# Subplot 1: Acceptance Angle vs. Normalized Core Radius
for alpha in alpha_values:
    theta_a = acceptance_angle(r_a, alpha, n1, n2, delta)
    ax1.plot(r_a, theta_a, label=f'α = {alpha}')
ax1.set_xlabel('Normalized Core Radius (r/a)')
ax1.set_ylabel('Acceptance Angle (θa, degrees)')
ax1.set_title('Acceptance Angle vs. Core Radius')
ax1.legend()
ax1.grid(True)

# Subplot 2: Refractive Index vs. Normalized Core Radius
for alpha in alpha_values:
    n_r = refractive_index(r_a, alpha, n1, delta)
    ax2.plot(r_a, n_r, label=f'α = {alpha}')
ax2.set_xlabel('Normalized Core Radius (r/a)')
ax2.set_ylabel('Refractive Index (n)')
ax2.set_title('Refractive Index vs. Core Radius')
ax2.set_ylim(1.47, 1.51)  # Set y-axis limits for clarity
ax2.legend()
ax2.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()

