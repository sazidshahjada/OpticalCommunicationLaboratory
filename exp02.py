"""
Exp-2. For Graded Index Fibers:
a. Calculate the Numerical aperture(NA), Propagation constant(β), V number
and Cutoff wavelength(λc).
b. Observe the graphical representation of Propagation constant vs V number.
c. Observe the graphical representation of Wavelength (λ) vs V number.
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Parameters for graded-index fiber
n1 = 1.48  # Core refractive index
n2 = 1.45  # Cladding refractive index
a = 5e-6   # Core radius in meters
lambda_0 = 0.65e-6  # Wavelength in meters
delta = (n1**2 - n2**2) / (2 * n1**2) # Relative index difference

# Calculate Numerical Aperture (NA)
NA = np.sqrt(n1**2 - n2**2)

# Calculate V number
V = (2 * np.pi * a * NA) / lambda_0

# Calculate propagation constant (β) for LP_01 mode
n_eff = n1 * np.sqrt(1 - 2 * delta / np.sqrt(V))
k = 2 * np.pi / lambda_0
beta = k * n_eff

# Calculate cutoff wavelength (λc) for single mode (V ≈ 3.4 for parabolic profile)
V_c = 3.4
lambda_c = (2 * np.pi * a * NA) / V_c

print(f"Numerical Aperture (NA): {NA:.4f}")
print(f"V number: {V:.4f}")
print(f"Propagation constant (β): {beta:.4e} m^-1")
print(f"Cutoff wavelength (λc): {lambda_c*1e6:.4f} μm")



# Sweep wavelength for plots
lambdas = np.linspace(500e-9, 2000e-9, 100)
Vs = (2 * np.pi * a * NA) / lambdas
n_effs = n1 * np.sqrt(1 - 2 * delta / np.sqrt(Vs))
ks = 2 * np.pi / lambdas  # Use wavelength-dependent k
betas = ks * n_effs

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