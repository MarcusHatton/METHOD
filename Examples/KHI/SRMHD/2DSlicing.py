import h5py
import numpy as np
import matplotlib.pyplot as plt

# Load HDF5 file
filename = "3d/KHI/dp_80x80x80_9.hdf5"
with h5py.File(filename, 'r') as f:
    print("Available keys:", list(f.keys()))
    # Try to find the density field (common names: 'rho', 'density', 'prims', etc.)
    rho = f['Primitive/rho'][:]

# Check shape
print("rho shape:", rho.shape)
# Pick center slices
nx, ny, nz = rho.shape
x_slice = rho[nx // 2, :, :]
y_slice = rho[:, ny // 2, :]
z_slice = rho[:, :, nz // 2]
# Plot slices
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for ax, data, title in zip(axs, [x_slice, y_slice, z_slice], ['X slice', 'Y slice', 'Z slice']):
    im = ax.imshow(data.T, origin='lower', cmap='viridis')
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
plt.suptitle("Density Slices (rho)")
plt.tight_layout()
plt.savefig('plots/3D_Density_Slices')
plt.show()
