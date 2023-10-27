from ISCE_plotting import plot_khi
import h5py

paper_width=(3.65,2.3)

fss = []

for n in range(6):
    fss.append(h5py.File(f'2d/Shear/dp_800x1600x0_{n}.hdf5', 'r'))

print(fss[0]['Optional'].attrs['eta'][0])

fig = plot_khi(fss, group_str = 'Primitive', quant_str='n', quant_legend_str = r"$n$",
             model_str='MISCE', coeff_str='eta', coeff_legend_str = r"$\eta$", figsize=paper_width)

#fig.savefig('MISCE_KHI_Ideal_t30_800x1600_n.pdf',bbox_inches='tight',dpi=1200)
fig.savefig('testing.pdf')
