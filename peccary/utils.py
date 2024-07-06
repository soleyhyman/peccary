"""
Utility functions for PECCARY
"""

import numpy as np
from math import factorial
import matplotlib.pylab as plt

__all__ = ["ell2tpat", "tpat2ell"]

def ell2tpat(ell,n,dt):
    """
    Convert sampling interval to pattern timescale

    Parameters
    ----------
    ell : float
        Sampling interval
    n : float
        Sampling size
    dt : float
        Timestep or timeseries resolution

    Returns
    -------
    float
        Equivalent pattern timescale
    """
    return ell*dt*(n-1.)

def tpat2ell(tpat,n,dt):
    """
    Convert pattern timescale to sampling interval

    Parameters
    ----------
    tpat : float
        Pattern timescale
    n : float
        Sampling size
    dt : float
        Timestep or timeseries resolution

    Returns
    -------
    float
        Equivalent sampling interval
    """
    return tpat/(dt*(n-1.))

def HmaxPer(n=5): 
    """
    Calculate maximum Permutation Entropy value (H) for a
    periodic function given the specified sampling size n.
    
    Equation: [insert equation here]

    Parameters
    ----------
    n : int
        Sampling size, by default 5

    Returns
    -------
    float
        Maximum Permuation Entropy for periodic function
    """
    Nper = 2.*(2.*(n-2.)+1.)
    return np.log(Nper)/np.log(factorial(n))

def HminPer(n=5):
    """
    Calculate minimum Permutation Entropy value (H) for a
    periodic function given the specified sampling size n.
    
    Equation: [insert equation here]

    Parameters
    ----------
    n : int
        Sampling size, by default 5

    Returns
    -------
    float
        Minimum Permuation Entropy for periodic function
    """
    return np.log(2)/np.log(factorial(n))