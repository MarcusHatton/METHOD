import h5py
import matplotlib.pyplot as plt

paper_width=(3.65,2.3)

fss = []

for n in range(0,6,1):
    fss.append(h5py.File(f'2d/KHRandom/dp_200x200x0_{n}.hdf5', 'r'))

fig, axes = plt.subplots(2,3,figsize=paper_width)
axes = axes.flatten()

for fs, ax in zip(fss, axes):
    ax.imshow(fs['Conserved/D'][:])

#fig.savefig('MISCE_TI_KHI_Ideal_t30_800x800_n.pdf',bbox_inches='tight',dpi=1200)
plt.savefig('200x200_Shear_t20.pdf',bbox_inches='tight')
