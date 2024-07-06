"""
Utility functions for PECCARY
"""

import numpy as np
from math import factorial
import matplotlib.pylab as plt

__all__ = ["ell2tpat", "tpat2ell", "HmaxPer", "HminPer", "calcHCplane", "getMaxC"]

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

def calcHCplane(n=5, nsteps=1000):	
    """
    Get maximum and minimum C(H) curves based on pattern length
    for plotting on the HC plane

    Parameters
    ----------
    n : int, optional
        Embedding dimension/pattern length, by default 5
    nsteps : int, optional
        Number of steps to use for generating HC bounding curves,
        by default 1000

    Returns
    -------
    4-tuple
        4-tuple of ndarrays consisting (respectively) of:
        - H values for minimum HC curve
        - C values for minimum HC curve
        - H values for maximum HC curve
        - C values for maximum HC curve
    """
    nsteps = nsteps
    n = n
    N = factorial(n)
    invN = 1./N
    log2_N = np.log2(N)
    log2_Np1 = np.log2(N+1.)     

    # Set up blank arrays for x- and y-values of max/min curves
    Cmaxx = np.zeros((N-1)*nsteps)
    Cmaxy = np.zeros((N-1)*nsteps)
    Cminx = np.zeros(nsteps)
    Cminy = np.zeros(nsteps)
    
    # Calculate H and C values for minimum HC curve
    for i in np.arange(nsteps):
        pk = invN + i*(1.-(invN))/nsteps
        pj = (1. - pk)/(N - 1.)
        S = -pk * np.log2(pk) - (N - 1.) * pj * np.log2(pj)
        qk = pk/2. + 1./(2.*N)
        qj = pj/2. + 1./(2.*N)
        Scom = -qk * np.log2(qk) - (N - 1.) * qj * np.log2(qj)
        Cminx[i] = S / log2_N
        Cminy[i] = -2. * (S/log2_N) * (Scom - 0.5*S - 0.5*log2_N)\
        /((1 + invN)*log2_Np1 - 2*np.log2(2.*N) + log2_N)	
        
    # Calculate H and C values for maximum HC curve
    for i in np.arange(1,N):
        for l in np.arange(nsteps):
            pk = l*(1./(N-i+1.))/nsteps
            pj = (1. - pk)/(N - i)
            if pk ==0.:
                S = -(N - i) * pj * np.log2(pj)
            else:
                S = -pk * np.log2(pk) - (N - i) * pj * np.log2(pj)
            qk = pk/2. + 1./(2.*N)
            qj = pj/2. + 1./(2.*N)
            Scom = -qk * np.log2(qk) - (N - i) * qj * np.log2(qj) - \
            (i-1)*(1./(2.*N))*np.log2(1./(2.*N))
            Cmaxx[(i-1)*nsteps+l] = S / log2_N
            Cmaxy[(i-1)*nsteps+l] = -2.*(S/log2_N)*(Scom - 0.5*S - 0.5*log2_N) \
            /((1. + invN)*log2_Np1 - 2.*np.log2(2.*N) + log2_N)
            
    return Cminx, Cminy, Cmaxx, Cmaxy

def getMaxC(n=5, nsteps=1000):
    """
    Get maximum possible Statistical Complexity value
    from HC bounding curves for sampling size n

    Parameters
    ----------
    n : int, optional
        Embedding dimension/pattern length, by default 5
    nsteps : int, optional
        Number of steps to use for generating HC bounding curves,
        by default 1000

    Returns
    -------
    float
        Maximum possible C value for given n
    """
    Cminx, Cminy, Cmaxx, Cmaxy = calcHCplane(n=n, nsteps=nsteps)
    return np.max(Cmaxy)