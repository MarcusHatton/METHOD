import h5py
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np

fs = []

for n in range(0,6,1):
    fs.append(h5py.File(f'3d/KHI/dp_80x80x80_{n}.hdf5', 'r'))

Nx, Ny, Nz = 80, 80, 80
X, Y, Z = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(Nz))

values = fs[-1]['Primitive/rho'][:]

fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=values.flatten(),
    isomin=0.1,
    isomax=0.8,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
    ))
#fig.show()
fig.write_image("80x80x80_3D_go.png")

plt.savefig('80x80x80_3D_go.pdf',bbox_inches='tight')
