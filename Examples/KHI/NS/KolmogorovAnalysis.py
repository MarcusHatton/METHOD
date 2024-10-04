# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 21:38:23 2024

@author: marcu
"""

from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import h5py
import glob
import os
import pickle
from skimage.util import view_as_blocks
#from scipy.ndimage import uniform_filter

def SimpleCoarseGrainer(data, CGfactors=(1,1)):
    block_view = view_as_blocks(data,CGfactors)
    coarse_shape = np.shape(block_view)[:2] 
    # Don't want sub-shape within blocks to be averaged
    coarse_data = np.zeros(coarse_shape)
    for i in range(coarse_shape[0]):
        for j in range(coarse_shape[1]):
            coarse_data[i][j] = np.average(block_view[i][j])
    return coarse_data

def getFourierTrans(u):

    nx, ny = np.shape(u)
    NN = nx // 2

    uhat = np.zeros((NN, ny), dtype=np.complex_)
    for k in range(NN):
        for y in range(ny):
            # Sum over all x adding to uhat
            for i in range(nx):
                uhat[k, y] += u[i, y] * np.exp(-(2*np.pi*1j*k*i)/nx)

    return uhat / nx

def getFourierTrans(intPlot, u, CG_factors=(1,1)):
    """
    Returns the 1D discrete fourier transform of the variable u along the x-direction
    ready for the power spectrum method.
    Parameters
    ----------
    intPlot : object
        interactivePlot object containing all the simulation data, normally the final instance
    u : ndarray
        Two dimensional array of the variable we want the power spectrum of
    Returns
    -------
    uhat : array (N,)
        Fourier transform of u
    """
    nx, ny = intPlot['Domain'].attrs['nx'][0] // CG_factors[0], intPlot['Domain'].attrs['ny'][0] // CG_factors[1]
    NN = nx // 2
    print(NN, nx, ny)
    print(np.shape(u))
    uhat = np.zeros((NN, ny), dtype=np.complex_)

    for k in range(NN):
        for y in range(ny):
            # Sum over all x adding to uhat
            for i in range(nx):
                uhat[k, y] += u[i, y] * np.exp(-(2*np.pi*1j*k*i)/nx)

    return uhat / nx

def getPowerSpectrumSq(intPlot, u, CG_factors=(1,1)):
    """
    Returns the integrated power spectrum of the variable u, up to the Nyquist frequency = nx/2
    Parameters
    ----------
    intPlot : object
        interactivePlot object containing all the simulation data, normally the final instance
    u : ndarray
        Two dimensional array of the variable we want the power spectrum of
    """
    NN = intPlot['Domain'].attrs['nx'][0] // (2*CG_factors[0])
    dy = intPlot['Domain'].attrs['dy'][0] * CG_factors[1]
    uhat = getFourierTrans(intPlot, u, CG_factors)
    P = np.zeros(NN)

    for k in range(NN):
        for j in range(intPlot['Domain'].attrs['ny'][0] // CG_factors[1]):
            P[k] += (np.absolute(uhat[k, j])**2) * dy

    P = P / np.sum(P)
    return P

#def getPowerSpectrumSq(u):

#    nx, ny = np.shape(u)
#    NN = nx // 2

#def GetKESF(rho, W):
#    return rho * W * (W-1)

def GetKESF(frame, CG_factors=(1,1)):
    """
    Retrieves and computes the kinetic energy density for each frame in a single fluid animation.
    Parameters
    ----------
    anim : object
        animation class containing all user def variables
    frame : Array
        Frame from the animation class containing all user def variables at the time we want
    """
    #Nx = frame['Domain'].attrs['nx'][0] // 2

    # Why is this needed?
    #vx = SimpleCoarseGrainer(frame['Primitive/v1'][:], CG_factors)
    #vy = frame['Primitive/v2'][:]
    #vsq = vx**2 + vy**2
    #W = 1 / np.sqrt(1 - vsq)

    rho = SimpleCoarseGrainer(frame['Primitive/rho'][:], CG_factors)
    W = SimpleCoarseGrainer(frame['Auxiliary/W'][:], CG_factors)

    KE = rho * W * (W-1) ## this expression ??

    #TE = gamma*m_0*c**2 = rho
    #KE = gamma*m_0*c**2 - m_0*c**2 - e*n*c**2  = (gamma-1)*m_0*c**2 - e*n*c**2 = (W-1)*n - e*n

    return KE

if __name__ == "__main__":


    fs = []
    n_files = 11
    for n in range(n_files):
        fs.append(h5py.File(f'./2d/KHRandom/Ideal/dp_400x400x0_{n}.hdf5', 'r'))
    # print(fs)
    
    KE_Spec = 0
    f = fs[5]
    CG_factors = (1,1)
    KE_Spec = getPowerSpectrumSq(f, GetKESF(f, CG_factors), CG_factors)

    with open('./pickles/KESpec_400x400_t=10.pickle', 'wb') as filehandle:
        pickle.dump(KE_Spec, filehandle)

    Nx = f['Domain'].attrs['nx'][0] // (2*CG_factors[0])
    plt.loglog(np.arange(1, Nx+1), np.arange(1, Nx+1)*KE_Spec)
    #plt.loglog([5.0, 60.0], [5.0*10**-3, (5.0*10**-3)*(12.0**(-5/3))], 'k--')
    plt.loglog(np.arange(1, Nx+1), (np.arange(1, Nx+1)**(8/3))*np.arange(1, Nx+1)*KE_Spec, 'k--')
    plt.annotate(r'$k^{-5/3}$', xy=(20, 0.01), fontsize=15)
    plt.ylabel(r"$k|P_{T}(k)|^2$", {'fontsize':'large'})
    plt.xlabel(r'$k$')
    plt.savefig('./plots/KESpec_400x400_t=10.pdf')
    plt.close()
































