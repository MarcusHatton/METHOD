import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import os

fss = []

for n in range(4):
    #fss.append(h5py.File(f'2d/Shear/5em4_1em3/t_6.25/dp_400x800x0_{n}.hdf5', 'r'))
    #fss.append(h5py.File(f'2d/Shear/1em3_5em3/t_6.25/dp_400x800x0_{n}.hdf5', 'r'))
    #fss.append(h5py.File(f'2d/Shear/1em3_5em3/dp_800x1600x0_{n}.hdf5', 'r'))
    fss.append(h5py.File(f'2d/Shear/5em4_1em3/dp_800x1600x0_{n}.hdf5', 'r'))
    #fss.append(h5py.File(f'2d/Shear/1em3_5em3/dp_400x800x0_{n}.hdf5', 'r'))
    #fss.append(h5py.File(f'2d/Ideal/dp_800x1600x0_{n}.hdf5', 'r'))

#Options
#params = {'text.usetex' : True, 'font.size' : 20, 'font.family' : 'lmodern'              }
#plt.rcParams.update(params) 
#plt.rcParams.update(plt.rcParamsDefault)

quant_str= 'n'

fig, axes = plt.subplots(3,2,figsize=(20,16),sharex=True, sharey=True)#, constrained_layout=True)
for i, fs, axis in zip(range(6), fss, axes.flatten()):
    axes = axes.flatten()
    extent = [-1.0,1.0,-0.5,0.5]
    axes[i].imshow(fs['Primitive/'+quant_str][:], extent=extent)
    axes[i].set_title(r'$t = ~ $'+str(fs.attrs['t'][0]))
    if i == 0:
        #axes[i].set_title('MISCE, '+quant_str+' ('+r'$t = ~$'+str(fs.attrs['t'][0])+')')
        axes[i].set_xticks([-1.0,-0.5,0.0,0.5,1.0])
        axes[i].set_yticks([-0.5,-0.25,0.0,0.25,0.5])
        axes[i].tick_params(axis='x', which='both', bottom=False)
        axes[i].set_ylabel(r'$x$')
    if i == 1:
        axes[i].tick_params(axis='x', which='both', bottom=False)
        axes[i].tick_params(axis='y', which='both', left=False)
        #axes[i].set_title(r'$t = ~$'+str(fs.attrs['t'][0])+r'$, ~ \eta=1 \times 10^{-15}$')
    if i == 2:
        axes[i].tick_params(axis='x', which='both', bottom=False)
        axes[i].set_ylabel(r'$x$')
    if i == 3:
        axes[i].tick_params(axis='x', which='both', bottom=False)
        axes[i].tick_params(axis='y', which='both', left=False)
        #pass
    if i == 4:
        axes[i].set_ylabel(r'$x$')
        axes[i].set_xlabel(r'$y$')
    if i == 5:
        axes[i].set_xlabel(r'$y$')
        axes[i].tick_params(axis='y', which='both', left=False)

#axes[0,0].set_title(r'$Viscous, 400x800$') 
#axes[0,0].set_title(r'$Viscous, 400x800$') 
#axes[0,0].set_title(r'$Ideal, 800x1600$') 
#axes[0,0].set_title(r'$Ideal, 200x400$') 
#axes[0,0].set_title(r'$Ideal, 80x160$') 
#axes[-1].set_title(r'$t = 6.25$')

fig.tight_layout()
#plt.savefig('IS_KHI_Shear_1em3_5em3_400x800_n.pdf', dpi=800)
#plt.savefig('IS_KHI_Shear_5em4_1em3_400x800_n.pdf', dpi=1200)
plt.savefig('testing.pdf')
