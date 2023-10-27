from ISCE_plotting import plot_khi
import h5py

paper_width=(3.65,2.3)

fss = []

for n in range(0,11,2):
    fss.append(h5py.File(f'2d/Ideal/TI/dp_200x200x0_{n}.hdf5', 'r'))
    #print(n)

#print(fss)

fig = plot_khi(fss, group_str = 'Primitive', quant_str='n', quant_legend_str = r"$n$",
             model_str='MISCE', coeff_str='eta', coeff_legend_str = r"$\eta$", figsize=paper_width)

#fig.savefig('MISCE_TI_KHI_Ideal_t30_800x800_n.pdf',bbox_inches='tight',dpi=1200)
fig.savefig('testing.pdf',bbox_inches='tight')
