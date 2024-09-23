# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 21:38:23 2024

@author: Marcus
"""

from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import h5py
import glob
import os
import pickle
#from os import system
import sys as system
import seaborn as sns

if __name__ == "__main__":


    fss = []
    for n in range(0,11,1):
        fss.append(h5py.File(f'SubGrid/2d/KHRandom/dp_200x200x0_{n}.hdf5', 'r'))

    data_to_plot = []
    all_data = np.array([])

    for fs in fss[::]:
        data_to_plot.append(np.array(fs['Auxiliary/zeta'][:]))

    #t_labels = [r'$t=2.0$',r'$t=4.0$',r'$t=6.0$',r'$t=8.0$',r'$t=12.0$',r'$t=16.0$']
    #t_labels = [r'$t=0.0$', r'$t=2.0$',r'$t=4.0$',r'$t=6.0$',r'$t=8.0$',r'$t=10.0$']
    t_labels = [r'$t=0.0$', r'$t=4.0$',r'$t=8.0$',r'$t=12.0$',r'$t=16.0$',r'$t=20.0$']

    fig, axes = plt.subplots(3,4,figsize=(12,8), sharex=True, sharey=True)#, constrained_layout=True)
    for ax, data, i in zip(axes.flatten(), data_to_plot, range(len(data_to_plot))):
        data = data[~np.isnan(data)]
        data = data[data != 0.0]
        data = data[data != np.inf]
        all_data = np.concatenate((all_data,data))
        #print(data)
        #print(np.min(data), np.max(data))
        #ax.hist(data.flatten())

    #plt.savefig('DissCoeffStats_Testing.pdf', bbox_inches='tight',dpi=500)
    #plt.close()
    #exit()

    plt.figure()
    sns.displot(all_data, log_scale=True)
    #plt.xlim(3e-6, 3e-4)
    plt.savefig('SG_DisPlot_Zeta_200x200_t20.pdf', bbox_inches='tight',dpi=1200)
    
























