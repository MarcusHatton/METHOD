from MIS_plotting import plot_khi
import h5py

paper_width=(3.65,2.3)

fss = []

for n in range(0,6,1):
    fss.append(h5py.File(f'/scratch/mjh1n20/KHI/MIS/Shear/t30/dp_800x1600x0_{n}.hdf5', 'r'))
    #fss.append(h5py.File(f'2d/Shear/1em3_5em3/dp_400x800x0_{n}.hdf5', 'r'))
    #print(n)

#print(fss)

fig = plot_khi(fss, group_str = 'Primitive', quant_str='n', quant_legend_str = r"$n$",
             model_str='MISCE', coeff_str='eta', coeff_legend_str = r"$\eta$", figsize=paper_width)

#fig.savefig('MIS_KHI_Ideal_400x800_n.pdf',bbox_inches='tight',dpi=1200)
fig.savefig('MIS_KHI_Shear_1em3_5em3_800x1600_t30_n.pdf',bbox_inches='tight',dpi=1200)
#fig.savefig('testing.pdf',bbox_inches='tight')
