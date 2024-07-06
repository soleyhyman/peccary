"""
Generating and plotting on HC plane.

The class ``HCplots`` is used for plotting the HC plane and
for creating the bases of other plot formats that can be useful for
diagnostics, interpretation, and analysis with PECCARY.
"""

import numpy as np
from math import factorial
import matplotlib.pylab as plt

### USE THIS LINE FOR DISTRIBUTION ###
# import peccary.utils

### USING THIS LOCAL VERSION FOR NOW ###
import utils

__all__ = ["plotBoundsHC","HCplane"]

# def plotBoundsHC(ax, n=5, nsteps=1000, lw=1., ls='--', color='k', alpha=0.5):
def plotBoundsHC(ax, n=5, nsteps=1000, **kwargs):
    """
    Plots region boundary lines on HC plane

    Parameters
    ----------
    ax : Matplotlib axis
        Axis on which to plot empty HC curves
    n : int, optional
        Embedding dimension/pattern length, by default 5
    nsteps : int, optional
        Number of steps to use for generating HC bounding curves,
        by default 1000
    **kwargs : dict, optional
        Style arguments passed to ``matplotlib.pyplot.plot``,
        if none given uses PECCARY defaults
    """
    minH, minC, maxH, maxC = utils.calcHCplane(n=n, nsteps=nsteps) # H and C values for HC min/max curves
    minHVal = utils.HminPer(n=n) # Min. H value for periodic fnct
    maxHVal = utils.HmaxPer(n=n) # Max. H value for periodic fnct
    argMinH_CminCurve = np.argmin(np.abs(minH-minHVal)) # Min. H value argument for H values for min HC curve
    argMinH_CmaxCurve = np.argmin(np.abs(maxH-minHVal)) # Min. H value argument for H values for max HC curve
    argMaxH_CminCurve = np.argmin(np.abs(minH-maxHVal)) # Max. H value argument for H values for min HC curve
    argMaxH_CmaxCurve = np.argmin(np.abs(maxH-maxHVal)) # Max. H value argument for H values for max HC curve

    if 'lw' not in kwargs.keys():
        kwargs['lw'] = 1.
    if 'ls' not in kwargs.keys():
        kwargs['ls'] = '--'
    if 'color' not in kwargs.keys():
        kwargs['color'] = 'k'
    if 'alpha' not in kwargs.keys():
        kwargs['alpha'] = 0.5

    ax.plot(np.array([minHVal,minHVal]), np.array([minC[argMinH_CminCurve], maxC[argMinH_CmaxCurve]]), **kwargs)
    ax.plot(np.array([maxHVal,maxHVal]), np.array([minC[argMaxH_CminCurve], maxC[argMaxH_CmaxCurve]]), **kwargs)
    ax.plot(np.array([minHVal,maxHVal]), np.array([maxC[argMinH_CmaxCurve], maxC[argMaxH_CmaxCurve]]), **kwargs)

def HCplane(H=None, C=None, ax=None, n=5, nsteps=1000, fontsize=12, showAxLabels=True, showBoundaryLines=True, savePlot=False, savePath='', **kwargs):
    """
    Creates a blank HC plane with maximum and minimum curves for the given embedding dimension

    Parameters
    ----------
    ax : Matplotlib axis
        Axis on which to plot empty HC curves
    n : int, optional
        Embedding dimension/pattern length, by default 5
    nsteps : int, optional
        Number of steps to use for generating HC bounding curves,
        by default 1000
    fontsize : integer or float, optional
        Fontsize of axis labels, by default 12
    showAxLabels : bool, optional
        Show pre-defined axis labels, by default True
    showBoundaryLines : bool, optional
        Show HC plane boundary lines, by default True
    savePlot : bool, optional
        Saves HC plot if set to True, by default False
    savePath : str, optional
        Path to save plot if savePlot set to True, by default ''
        Note: Use only forward slashes in savePath
    **kwargs : dict, optional
        Style arguments for region boundary lines passed to 
        ``matplotlib.pyplot.plot``, if none given uses PECCARY defaults
    """
    # Calculate HC plane envelope
    Cminx, Cminy, Cmaxx, Cmaxy = utils.calcHCplane(n=n, nsteps=nsteps)

    # Check if an existing axis has been inputted
    # Otherwise, create figure and subplot
    if ax is None:
        fig, ax = plt.subplots(1,1)
        returnAx = True
    else:
        returnAx = False

    # Plot HC envelope
    ax.plot(Cminx,Cminy,'k-',Cmaxx,Cmaxy,'k-', zorder=0)

    # Check whether or not to plot HC plane boundaries
    if showBoundaryLines:
        plotBoundsHC(ax, **kwargs)
    else:
        pass

    # Check whether H and C values have been inputted
    # If so, plot them
    if (H is None) or (C is None):
        pass
    else:
        ax.scatter(H, C)

    # Check whether to include axis labels
    if showAxLabels:
        ax.set_xlabel(r"Normalized Permutation Entropy, $H$", fontsize=fontsize)
        ax.set_ylabel(r"Statistical Complexity, $C$", fontsize=fontsize)
        ax.tick_params(axis='both', labelsize=fontsize-2)
    else:
        pass
    
    # Check whether to save plot
    if savePlot:
        ax.set(xlim=(0,1.0), ylim=(0,0.45))
        ax.set_xticks(np.arange(0,1.1,0.1))
        ax.set_yticks(np.arange(0,0.45,0.05))
        if savePath.endswith('/'):
            savefile=savePath + 'HCn5_blank'
        else:
            savefile=savePath + '/HCn5_blank'
        plt.savefig(savefile+'.png')
    elif returnAx:
        return fig, ax