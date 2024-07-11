from ISCE_plotting import plot_khi
import h5py

paper_width=(3.65,2.3)

fss = []

for n in range(0,11,2):
    fss.append(h5py.File(f'SubGrid/CoefficientTesting/2d/KHRandom/dp_200x200x0_{n}.hdf5', 'r'))
    #print(n)

#print(fss)

#fig = plot_khi(fss, group_str = 'Primitive', quant_str='n', quant_legend_str = r"$n$",
#             model_str='MISCE', coeff_str='eta', coeff_legend_str = r"$n$")#, figsize=paper_width)
#fig.savefig('testing.pdf',bbox_inches='tight')

#fig = plot_khi(fss, group_str = 'Auxiliary', quant_str='kappa', quant_legend_str = r"$\kappa$",
#             model_str='NS', coeff_str='kappa', coeff_legend_str = r"$\kappa$")#, figsize=paper_width)
#fig.savefig('testing_kappa.pdf',bbox_inches='tight')

fig = plot_khi(fss, group_str = 'Auxiliary', quant_str='zeta', quant_legend_str = r"$\zeta$",
             model_str='NS', coeff_str='zeta', coeff_legend_str = r"$\zeta$")#, figsize=paper_width)
fig.savefig('testing_zeta.pdf',bbox_inches='tight')

fig = plot_khi(fss, group_str = 'Auxiliary', quant_str='eta', quant_legend_str = r"$\eta$",
             model_str='NS', coeff_str='eta', coeff_legend_str = r"$\eta$")#, figsize=paper_width)
fig.savefig('testing_eta.pdf',bbox_inches='tight')
