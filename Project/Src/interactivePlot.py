"""
    Script gathers the state vectors stored in the Data directory and offers 
    functionality to plot various elements.
"""


import numpy as np
import pylab as plt
from matplotlib import cm
import warnings

warnings.filterwarnings('ignore', "No labelled objects found. ")

# Change this to the relative path to the data you want to plot
# File names must start with e.g. `primitive`, anything between this
# and `.dat` should be stored in appendix
DataDirectory = '../Data/'
appendix = ''

def gatherData():
    """
    Collects and stores all the data required for plotting the final state of
    the system.
    
    Returns
    -------
        cons : array of float
            (Ncons, Nx, Ny, Nz) Array containing the conserved vector
        consLabels : array of string
            (Ncons,) The labels of the conserved elements
        prims : array of float
            (Nprims, Nx, Ny, Nz) Array containing the primitive vector
        primLabels : array of string
            (Nprims,) The labels of the primitive elements
        aux : array of float
            (Naux, Nx, Ny, Nz) Array containing the auxilliary vector
        auxLabels : array of string
            (Naux,) The labels of the auxilliary elements
        c : dictionary
            Dictionary containing all constant data saved in simData. Access
            elements by typing as an argument the constant you want as a string.
            E.g. to get zmax, enter  -->    c['zmax']
            All links are the same as the constant name in the SimData class.
    
    """
    
    # Dictionary to hold constants
    c = {}
    # Get constants first
    print("Fetching constants...")
    with open(DataDirectory + 'constants' + appendix + '.dat', 'r') as f:
        for i, line in enumerate(f):
            if not i==0:
                line=line.split()
                c['nx'] = int(line[0])
                c['ny'] = int(line[1])
                c['nz'] = int(line[2])
                c['Nx'] = int(line[3])
                c['Ny'] = int(line[4])
                c['Nz'] = int(line[5])
                c['xmin'] = float(line[6])
                c['xmax'] = float(line[7])
                c['ymin'] = float(line[8])
                c['ymax'] = float(line[9])
                c['zmin'] = float(line[10])
                c['zmax'] = float(line[11])
                c['endTime'] = float(line[12])
                c['cfl'] = float(line[13])
                c['Ng'] = int(line[14])
                c['gamma'] = float(line[15])
                c['sigma'] = float(line[16])
                c['Ncons'] = int(line[17])
                c['Nprims'] = int(line[18])
                c['Naux'] = int(line[19])
                c['cp'] = float(line[20])
                c['dt'] = float(line[21])
                c['t'] = float(line[22])
                c['dx'] = float(line[23])
                c['dy'] = float(line[24])
                c['dz'] = float(line[25])

    print("{} conserved variables".format(c['Ncons']))
    print("{} primitive variables".format(c['Nprims']))
    print("{} auxilliary variables".format(c['Naux']))

    # Now get primitive variables and store the data in array...
    prims = np.zeros([c['Nprims'], c['Nx'], c['Ny'], c['Nz']])
    print("Fetching primitive variables...")
    with open(DataDirectory + 'primitive' + appendix + '.dat', 'r') as f:
        for i, line in enumerate(f):
            # Get primitive var labels
            if i==0:
                primLabels = line.split()[2:]
            # Get primitive var data
            else:
                temp = line.split()
                for k in range(c['Nz']):
                    prims[getVarFromLine(i, c['Nx'], c['Ny'])][getXIndexFromLine(i, c['Nx'], c['Ny'])][getYIndexFromLine(i, c['Nx'], c['Ny'])] = float(temp[k])

    # Clean up labels (remove the commas)
    cleanPrimLabels = []
    for i in range(len(primLabels)-1):
        cleanPrimLabels.append(primLabels[i][:-1])
    cleanPrimLabels.append(primLabels[-1])
    
    
    
    # Now gather conserved data
    cons = np.zeros([c['Ncons'], c['Nx'], c['Ny'], c['Nz']])
    print("Fetching conserved variables...")
    with open(DataDirectory + 'conserved' + appendix + '.dat', 'r') as f:
        for i, line in enumerate(f):
            # Get cons var labels
            if i==0:
                consLabels = line.split()[2:]
            # Get cons var data
            else:
                temp = line.split()
                for k in range(c['Nz']):
                    cons[getVarFromLine(i, c['Nx'], c['Ny'])][getXIndexFromLine(i, c['Nx'], c['Ny'])][getYIndexFromLine(i, c['Nx'], c['Ny'])] = float(temp[k])

    # Clean up labels (remove the commas)
    cleanConsLabels = []
    for i in range(len(consLabels)-1):
        cleanConsLabels.append(consLabels[i][:-1])
    cleanConsLabels.append(consLabels[-1])
    
    
    
    # And finally the aux vars
    aux = np.zeros([c['Naux'], c['Nx'], c['Ny'], c['Nz']])
    print("Fetching auxilliary variables...")
    with open(DataDirectory + 'auxilliary' + appendix +'.dat', 'r') as f:
        for i, line in enumerate(f):
            # Get cons var labels
            if i==0:
                auxLabels = line.split()[2:]
            # Get cons var data
            else:
                temp = line.split()
                for k in range(c['Nz']):
                    aux[getVarFromLine(i, c['Nx'], c['Ny'])][getXIndexFromLine(i, c['Nx'], c['Ny'])][getYIndexFromLine(i, c['Nx'], c['Ny'])] = float(temp[k])

    # Clean up labels (remove the commas)
    cleanAuxLabels = []
    for i in range(len(auxLabels)-1):
        cleanAuxLabels.append(auxLabels[i][:-1])
    cleanAuxLabels.append(auxLabels[-1])
    
    return cons, cleanConsLabels, prims, cleanPrimLabels, aux, cleanAuxLabels, c
    
    
    

def getVarFromLine(line, Nx, Ny):
    """
    Given the line number that the iterator is on, and the size of the x-domain,
    returns the index of the primitive variable this data belongs to.

    Parameters
    ----------
        line: int
            The line number the file pointer is pointing to. We want to know which
            primitive variable this line's data corresponds to.
        Nx: int
            The total number (incl ghost cells) of domain cells in the x-direction.
        Ny: int
            The total number (incl ghost cells) of domain cells in the y-direction.

    Returns
    -------
        var:
            The primitive variable index of this line's data.

    Other
    -----
        Function will throw a ValueError if trying to get the primitive index
        of the first (zero'th) line.
    """
    if line == 0:
        raise ValueError('Line zero does not contain any data')
    else:
        return ((line-1)//Ny)//Nx


def getXIndexFromLine(line, Nx, Ny):
    """
    Given the line number that the iterator is on, and the size of the x-domain,
    returns the x-index of this line's data.

    Parameters
    ----------
        line: int
            The line number the file pointer is pointing to. We want to know which
            primitive variable this line's data corresponds to.
        Nx: int
            The total number (incl ghost cells) of domain cells in the x-direction.
        Ny: int
            The total number (incl ghost cells) of domain cells in the y-direction.

    Returns
    -------
        index:
            The x-index of the current line's data.
    """
    return ((line-1)//Ny)%Nx

def getYIndexFromLine(line, Nx, Ny):
    """
    Given the line number that the iterator is on, and the size of the y-domain,
    returns the y-index of this line's data.

    Parameters
    ----------
        line: int
            The line number the file pointer is pointing to. We want to know which
            primitive variable this line's data corresponds to.
        Nx: int
            The total number (incl ghost cells) of domain cells in the x-direction.
        Ny: int
            The total number (incl ghost cells) of domain cells in the y-direction.

    Returns
    -------
        index:
            The y-index of the current line's data.
    """
    return (line-1)%Ny

def plotHeatMaps(data, dataLabels, c, color=None, axis=2):
    """
    Plots the 2D heatmap of the given data. The axes to be plotted can be
    selected via the axis parameter---this corresponds to the axis you want
    to ignore.
    
    Parameters
    ----------
        data: array of float
            (Nvars, Nx, Ny, Nz) The data to be plotted (e.g. prims)
        dataLabels: array of string
            (Nvars,) The corresponding labels for the elements of data
        c: dictionary
            The simulation constant data
        color: matplotlib color map
            The colour theme to be plotting in. This can take string arguments
            but best to stick to variants of cm.somecolourscheme
            E.g. cm.magma
        axis: int
            The axis the user wants to ignore.
            (0, 1, 2) = (x, y, z)
    """
    for i in range(data.shape[0]):
        fig = plt.figure()
        if (axis == 0):
            plotVars = data[i, c['Nx']//2, c['Ng']:-c['Ng'], c['Ng']:-c['Ng']]
            axisLabel1 = r'$y$'
            axisLabel2 = r'$z$'
        if (axis == 1):
            plotVars = data[i, c['Ng']:-c['Ng'], c['Ny']//2, c['Ng']:-c['Ng']]
            axisLabel1 = r'$x$'
            axisLabel2 = r'$z$'
        if (axis == 2):
            plotVars = data[i, c['Ng']:-c['Ng'], c['Ng']:-c['Ng'], c['Nz']//2]
            axisLabel1 = r'$x$'
            axisLabel2 = r'$y$'
            
        if color==None:
            color = cm.plasma
        surf = plt.imshow(plotVars, cmap=color, interpolation='nearest')
        plt.title(r'Time Evolution for {}: $t = {}$'.format(dataLabels[i], c['t']))
        plt.xlabel(axisLabel2)
        plt.ylabel(axisLabel1)
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.legend()
        plt.show()
        
        
def plotSlice(data, dataLabels, c, axis=0):
    """
    Plots the variation of data in the `axis` direction. 
    
    Parameters
    ----------
        data: array of float
            (Nvars, Nx, Ny, Nz) The data to be plotted (e.g. prims)
        dataLabels: array of string
            (Nvars,) The corresponding labels for the elements of data
        c: dictionary
            The simulation constant data
        color: matplotlib color map
            The colour theme to be plotting in. This can take string arguments
            but best to stick to variants of cm.somecolourscheme
            E.g. cm.magma
        axis: int, optional
            The axis the user wants to plot in.
            (0, 1, 2) = (x, y, z)
            Defaults to axis=0, x-direction.
    """
    Nx, Ny, Nz, Ng= c['Nx'], c['Ny'], c['Nz'], c['Ng']
    
    for i in range(len(data)):
        plt.figure()
        if (axis == 0):
            plotVars = data[i, Ng:-Ng, Ny//2, Nz//2]
            axisLabel = r'$x$'
            step = c['dx']
            n = c['nx']
            left, right = c['xmin'], c['xmax']
        if (axis == 1):
            plotVars = data[i, Nx//2, Ng:-Ng, Nz//2]
            axisLabel = r'$y$'
            step = c['dy']
            n = c['ny']
            left, right = c['ymin'], c['ymax']
        if (axis == 2):
            plotVars = data[i, Nx//2, Ny//2, Ng:-Ng]
            axisLabel = r'$z$'
            step = c['dz']
            n = c['nz']
            left, right = c['zmin'], c['zmax']
            
        ymin = np.min(plotVars)
        ymax = np.max(plotVars)
        rangeY = ymax - ymin
        ylower = ymin - 0.025 * rangeY
        yupper = ymax + 0.025 * rangeY
        xs = np.linspace(left + step/2, right - step/2, n)
        plt.plot(xs, plotVars)
        plt.title(r'Time Evolution for {}: $t = {}$'.format(dataLabels[i], c['t']))
        plt.xlabel(axisLabel)
        plt.ylabel(r'$q_{}(x)$'.format(i+1))
        plt.ylim((ylower, yupper))
        plt.legend(loc='lower center', fontsize=10)
        plt.show()


def plotTwoFluidAverage(data, dataLabels, c, axis=0):
    """
    Plots the variation of average data in the `axis` direction of the two fluids. 
    
    Parameters
    ----------
        data: array of float
            (Nvars, Nx, Ny, Nz) The data to be plotted (e.g. prims)
        dataLabels: array of string
            (Nvars,) The corresponding labels for the elements of data
        c: dictionary
            The simulation constant data
        color: matplotlib color map
            The colour theme to be plotting in. This can take string arguments
            but best to stick to variants of cm.somecolourscheme
            E.g. cm.magma
        axis: int, optional
            The axis the user wants to plot in.
            (0, 1, 2) = (x, y, z)
            Defaults to axis=0, x-direction.
    """
    Nx, Ny, Nz, Ng= c['Nx'], c['Ny'], c['Nz'], c['Ng']
    Nvars = len(data)
    if (Nvars == 18):
        print("Cannot take average of conserved variables. Conserved variables are summed values...")
        return None
    elif (Nvars == 16):
        print("Determining average values of primitive variables...")
        singles = 5
    elif (Nvars == 35):
        print("Determining average values of auxilliary variables...")
        singles = 10
    else:
        print("Cannot interpret the vector you have given as prims or aux.")
        print("Exiting...")
        return None
    
    avgPlotVars = (data[:singles] + data[singles:2*singles]) / 2
    
    for i in range(singles):
        plt.figure()
        if (axis == 0):
            plotVars = avgPlotVars[i, Ng:-Ng, Ny//2, Nz//2]
            axisLabel = r'$x$'
            step = c['dx']
            n = c['nx']
            left, right = c['xmin'], c['xmax']
        if (axis == 1):
            plotVars = avgPlotVars[i, Nx//2, Ng:-Ng, Nz//2]
            axisLabel = r'$y$'
            step = c['dy']
            n = c['ny']
            left, right = c['ymin'], c['ymax']
        if (axis == 2):
            plotVars = avgPlotVars[i, Nx//2, Ny//2, Ng:-Ng]
            axisLabel = r'$z$'
            step = c['dz']
            n = c['nz']
            left, right = c['zmin'], c['zmax']
            
        ymin = np.min(plotVars)
        ymax = np.max(plotVars)
        rangeY = ymax - ymin
        ylower = ymin - 0.025 * rangeY
        yupper = ymax + 0.025 * rangeY
        xs = np.linspace(left + step/2, right - step/2, n)
        plt.plot(xs, plotVars)
        plt.title(r'Time Evolution for {}: $t = {}$'.format(dataLabels[i], c['t']))
        plt.xlabel(axisLabel)
        plt.ylabel(r'$q_{}(x)$'.format(i+1))
        plt.ylim((ylower, yupper))
        plt.legend(loc='lower center', fontsize=10)
        plt.show()


# Function declarations over, access data and plot!

if __name__ == '__main__':

    cons, consLabels, prims, primLabels, aux, auxLabels, c = gatherData()
    

