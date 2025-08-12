from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import h5py
import glob
import os
import pickle

if __name__ == "__main__":

    pickle_dir = 'pickles/'
    pickle_files = ['KESpec_400x400_t=20.pickle', 'KESpec_400x400_t=20_CG=2,2.pickle', 'KESpec_200x200_t=20.pickle', \
    #'KESpec_400x400_t=20_CG=4,4.pickle', 'KESpec_400x400_t=20_CG=8,8.pickle', 'KESpec_400x400_t=20_CG=16,16.pickle',\
    'SG_KESpec_200x200_t=20.pickle']

    fig = plt.figure()

    for pickle_file in pickle_files:
        with open(pickle_dir+pickle_file, 'rb') as filehandle:
            KE_Spec = pickle.load(filehandle)
            Nx = len(KE_Spec)
  
        plt.loglog(np.arange(1, Nx+1), np.arange(1, Nx+1)*KE_Spec,label=pickle_file)
        #plt.loglog([5.0, 60.0], [5.0*10**-3, (5.0*10**-3)*(12.0**(-5/3))], 'k--')
        plt.loglog(np.arange(1, Nx+1), (np.arange(1, Nx+1)**(8/3))*np.arange(1, Nx+1)*KE_Spec, 'k--')
    plt.annotate(r'$k^{-5/3}$', xy=(20, 0.01), fontsize=15)
    plt.legend()
    plt.ylabel(r"$k|P_{T}(k)|^2$", {'fontsize':'large'})
    plt.xlabel(r'$k$')
    plt.savefig('./plots/KESpec_Compare3.pdf')
    plt.close()

