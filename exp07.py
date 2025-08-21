"""
Exp-7. For Silica Fibers, calculate the total loss and plot a graph for loss vs
wavelength(λ).
"""

import numpy as np
import matplotlib.pyplot as plt

# Wavelength range (0.8 to 1.8 μm)
wavelength = np.linspace(0.8, 1.8, 1000)

# Rayleigh scattering loss (A_R / λ^4)
A_R = 0.9
rayleigh_loss = A_R / (wavelength ** 4)

# UV absorption loss (A_UV * exp(-B_UV * λ))
A_UV = 0.1
B_UV = 5.0
uv_loss = A_UV * np.exp(-B_UV * wavelength)

# IR absorption loss (A_IR * exp(B_IR * (λ - λ_0)))
A_IR = 0.01
B_IR = 7.0
lambda_0 = 1.27
ir_loss = A_IR * np.exp(B_IR * (wavelength - lambda_0))

# Total loss
total_loss = rayleigh_loss + uv_loss + ir_loss

# Find minimum loss point
min_loss_idx = np.argmin(total_loss)
min_wavelength = wavelength[min_loss_idx]
min_loss = total_loss[min_loss_idx]

print(f"Minimum Loss: {min_loss:.2f} dB/km at Wavelength: {min_wavelength*1000:.3f} nm")

# Create the plot (normal scale)
plt.figure(figsize=(12, 8))

# Plot individual loss components (linear scale)
plt.plot(wavelength, rayleigh_loss, 'b--', label='Rayleigh Scattering', linewidth=2)
plt.plot(wavelength, uv_loss, 'g--', label='UV Absorption', linewidth=2)
plt.plot(wavelength, ir_loss, 'm--', label='IR Absorption', linewidth=2)
plt.plot(wavelength, total_loss, 'r-', label='Total Loss', linewidth=3)

# Mark minimum loss point
plt.plot(min_wavelength, min_loss, 'ro', markersize=8,
         label=f'Min Loss: {min_loss:.2f} dB/km @ {min_wavelength*1000:.3f} nm')

# Vertical line at minimum loss
plt.axvline(x=min_wavelength, color='gray', linestyle='--', linewidth=2)

# Formatting
plt.xlabel('Wavelength (μm)', fontsize=12)
plt.ylabel('Loss (dB/km)', fontsize=12)
plt.title('Silica Optical Fiber Loss vs Wavelength (Linear Scale)', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10)
plt.xlim(0.8, 1.8)

# Scale adjusted to fit the loss values better
plt.ylim(0, max(total_loss) * 1.1)

plt.tight_layout()
plt.show()