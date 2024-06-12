import h5py
import matplotlib.pyplot as plt
import numpy as np

paper_width=(3.65,2.3)

fss = []

for n in range(0,6,1):
    fss.append(h5py.File(f'3d/KHRandom/dp_40x40x40_{n}.hdf5', 'r'))
    #fss.append(h5py.File(f'2d/KHRandom/dp_200x200x0_{n}.hdf5', 'r'))

fig, axes = plt.subplots(2,3,figsize=paper_width)
axes = axes.flatten()

#for fs, ax in zip(fss, axes):
#    ax.imshow(fs['Primitive/n'][:])

print(fss[-1]['Primitive/n'].shape)
print(fss[-1]['Primitive/n'])

z_slices = np.arange(0,20,4)
print(z_slices)
for ax, z_slice in zip(axes, z_slices):
    ax.imshow(fss[-1]['Primitive/n'][:,:,z_slice])
    #print(fss[-1]['Primitive/n'][:,:,z_slice])


#fig.savefig('MISCE_TI_KHI_Ideal_t30_800x800_n.pdf',bbox_inches='tight',dpi=1200)
plt.savefig('40x40x40_Ideal_t10_slices_z.pdf',bbox_inches='tight')
