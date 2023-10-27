from ISCE_plotting import plot_shock_compare_wrapper
from ISCE_plotting import plot_shock_compare
import h5py

paper_width=(3.65,2.3)

fss = []

# load for time-series
#for n in [1,5]:
#    fss.append(h5py.File(f'1d/Bulk/5em2_1em1/ds_800_{n}.hdf5', 'r'))

# load for different dissipation values
n = 5
fss.append(h5py.File(f'1d/Ideal/ds_800_{n}.hdf5', 'r'))
#fss.append(h5py.File(f'1d/Bulk/5em2_1em1/ds_800_{n}.hdf5', 'r'))
fss.append(h5py.File(f'1d/Bulk/5em2_1em1/ds_800_{n}.hdf5', 'r'))
#fss.append(h5py.File(f'1d/Bulk/1em2_5em3/ds_800_{n}.hdf5', 'r'))
#fss.append(h5py.File(f'1d/Bulk/NLO_Testing/ds_800_{n}.hdf5', 'r'))

#plot_shock_compare_wrapper(fss,
#                       quant_strs=('n',),
#                       model_strs=('Ideal', 'MISCE',), coeff_strs=(None,'zeta'),
#                       figsize=paper_width, lss=('-', '--'))

fig=plot_shock_compare(fss, group_strs=('Primitive', 'Primitive', 'Auxiliary', 'Auxiliary'), quant_strs=('v1', 'rho', 'PiNS', 'T'), \
		#model_strs=('Bulk', 'Bulk'), quant_legend_strs=(r"$v_x$", r"$\rho$", r"$\Pi_{NS}$", r"$T$"))
		model_strs=('Ideal', 'Bulk'), quant_legend_strs=(r"$v_x$", r"$\rho$", r"$\Pi_{NS}$", r"$T$"))

fig.savefig('MISCE_ST_IdvsBulk_5em2_1em1_800_compare.pdf',bbox_inches='tight',dpi=1200)
#fig.savefig('testing.pdf',bbox_inches='tight')
