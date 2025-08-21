"""
Exp-1. For Step Index Fibers:
a) Calculate the Numerical aperture(NA), Propagation constant(β), V number
and Cutoff wavelength(λc).
b) Observe the graphical representation of Propagation constant vs V number.
c) Observe the graphical representation of Wavelength (λ) vs V number.
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Define the parameters
n1 = 1.48  # Refractive index of core
n2 = 1.45  # Refractive index of cladding
lambda_0 = 0.65e-6  # Wavelength in meters
a = 5e-6  # Core radius in meters
k = 2 * np.pi / lambda_0  # Wave number

# Calculate Numerical Aperture (NA)
NA = np.sqrt(n1**2 - n2**2)

# Calculate V number
V = (2 * np.pi * a * NA) / lambda_0

# Calculate propagation constant (β)
b = 1 - (1.1428 / V)**2
n_eff = np.sqrt(n2**2 + b * (n1**2 - n2**2))
beta = k * n_eff

# Calculate cutoff wavelength (λc) for single mode (V=2.405)
lambda_c = (2 * np.pi * a * NA) / 2.405

print(f"Numerical Aperture (NA): {NA:.4f}")
print(f"V number: {V:.4f}")
print(f"Propagation constant (β): {beta:.4e} m^-1")
print(f"Cutoff wavelength (λc): {lambda_c*1e6:.4f} μm")



# Sweep wavelength to get a range of V numbers and propagation constants
lambdas = np.linspace(500e-9, 2000e-9, 100)
Vs = (2 * np.pi * a * NA) / lambdas
bs = 1 - (1.1428 / Vs)**2
n_effs = np.sqrt(n2**2 + bs * (n1**2 - n2**2))
betas = k * n_effs

fig, axs = plt.subplots(1, 2, figsize=(14, 5))

# Subplot 1: Propagation Constant vs V number
sns.lineplot(x=Vs, y=betas, ax=axs[0], color='blue')
axs[0].set_xlabel("V Number")
axs[0].set_ylabel("Propagation Constant β (m$^{-1}$)")
axs[0].set_title("Propagation Constant vs V Number")
axs[0].grid(True)

# Subplot 2: Wavelength (λ) vs V number
sns.lineplot(x=lambdas * 1e6, y=Vs, ax=axs[1], color='red')
axs[1].set_xlabel("Wavelength (μm)")
axs[1].set_ylabel("V Number")
axs[1].set_title("Wavelength vs V Number")
axs[1].grid(True)

plt.tight_layout()
plt.show()