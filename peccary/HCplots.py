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
import peccary.utils

### USING THIS LOCAL VERSION FOR NOW ###
import utils

__all__ = ["HCplots"]

def getHC_boundary_lines(ax, n=5, nsteps=1000, lw=1., ls='--', color='k', alpha=0.5):
    """
    Creates boundary lines on HC plane

    Parameters
    ----------
    ax : Matplotlib axis
        Axis on which to plot empty HC curves
    n : int, optional
        Embedding dimension/pattern length, by default 5
    nsteps : int, optional
        Number of steps to use for generating HC bounding curves,
        by default 1000
    lw : float, optional
        Line width of boundary lines (must be float >0), by default 1.
    ls : string, optional
        Line style (Matplotlib linestyle string) of boundary lines, by default '--'
    color : string, optional
        Color of boundary lines (Matplotlib color string), by default 'k'
    alpha : float, optional
        Opacity of boundary lines (float between 0 and 1), by default 0.5
    """
    minH, minC, maxH, maxC = utils.calcBoundsHC(n=n, nsteps=nsteps) # H and C values for HC min/max curves
    minHVal = utils.HminPer(n=n) # Min. H value for periodic fnct
    maxHVal = utils.HmaxPer(n=n) # Max. H value for periodic fnct
    argMinH_CminCurve = np.argmin(np.abs(minH-minHVal)) # Min. H value argument for H values for min HC curve
    argMinH_CmaxCurve = np.argmin(np.abs(maxH-minHVal)) # Min. H value argument for H values for max HC curve
    argMaxH_CminCurve = np.argmin(np.abs(minH-maxHVal)) # Max. H value argument for H values for min HC curve
    argMaxH_CmaxCurve = np.argmin(np.abs(maxH-maxHVal)) # Max. H value argument for H values for max HC curve
    ax.plot(np.array([minHVal,minHVal]), np.array([minC[argMinH_CminCurve], maxC[argMinH_CmaxCurve]]), lw=lw, ls=ls, color=color, alpha=alpha)
    ax.plot(np.array([maxHVal,maxHVal]), np.array([minC[argMaxH_CminCurve], maxC[argMaxH_CmaxCurve]]), lw=lw, ls=ls, color=color, alpha=alpha)
    ax.plot(np.array([minHVal,maxHVal]), np.array([maxC[argMinH_CmaxCurve], maxC[argMaxH_CmaxCurve]]), lw=lw, ls=ls, color=color, alpha=alpha)  

def generateCurves(ax, n=5, nsteps=1000, lw=1., ls='--', color='k', alpha=0.5, fontsize=12, showAxLabels=True, showBoundaryLines=True, savePlot=False, savePath=''):
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
    lw : float, optional
        Line width of boundary lines (must be float >0), by default 1.
    ls : string, optional
        Line style (Matplotlib linestyle string) of boundary lines, by default '--'
    color : string, optional
        Color of boundary lines (Matplotlib color string), by default 'k'
    alpha : float, optional
        Opacity of boundary lines (float between 0 and 1), by default 0.5
    fonstize : integer or float, optional
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
    """
    Cminx, Cminy, Cmaxx, Cmaxy = utils.calcBoundsHC(n=n, nsteps=nsteps)

    ax.plot(Cminx,Cminy,'k-',Cmaxx,Cmaxy,'k-')

    if showBoundaryLines == True:
        getHC_boundary_lines(ax, lw=lw, ls=ls, color=color, alpha=alpha)
    else:
        pass

    if showAxLabels == True:
        ax.set_xlabel(r"Normalized Permutation Entropy, $H$", fontsize=fontsize)
        ax.set_ylabel(r"Statistical Complexity, $C$", fontsize=fontsize)
    else:
        pass
    
    if savePlot==True:
        ax.set(xlim=(0,1.0), ylim=(0,0.45))
        ax.set_xticks(np.arange(0,1.1,0.1))
        ax.set_yticks(np.arange(0,0.45,0.05))
        if savePath.endswith('/'):
            savefile=savePath + 'HCn5_blank'
        else:
            savefile=savePath + '/HCn5_blank'
        plt.savefig(savefile+'.png')