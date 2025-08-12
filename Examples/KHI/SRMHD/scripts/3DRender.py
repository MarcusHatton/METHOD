import h5py
import numpy as np
import pyvista as pv

filename = "3d/KHI/dp_40x40x40_9.hdf5" 
with h5py.File(filename, 'r') as f:
    rho = f['Primitive/rho'][:]

# Create a uniform grid
nx, ny, nz = rho.shape
grid = pv.UniformGrid()
grid.dimensions = np.array(rho.shape) + 1
grid.origin = (0, 0, 0)
grid.spacing = (1.0/nx, 1.0/ny, 1.0/nz)
grid.cell_data["rho"] = rho.flatten(order="F") # Fortran order

# Plot with volume rendering
plotter = pv.Plotter()
plotter.add_volume(grid, scalars="rho", opacity="linear", cmap="viridis")
plotter.add_axes()
plotter.show()
