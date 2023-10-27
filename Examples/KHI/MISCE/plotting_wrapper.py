from ISCE_plotting import plot_khi
import h5py

paper_width=(3.65,2.3)

fss = []

for n in range(6):
    fss.append(h5py.File(f'2d/1em5/dp_200x400x0_{n}.hdf5', 'r'))


fig = plot_khi(fss, group_str = 'Primitive', quant_str='n', quant_legend_str = r"$n$",
             model_str='MISCE', coeff_str='eta', coeff_legend_str = r"$\eta$", figsize=paper_width)

fig.savefig('MISCE_KHI_Shear_1em5_t6.25_200x400_n.pdf',bbox_inches='tight',dpi=1200)
