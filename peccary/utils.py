"""
Utility functions for PECCARY
"""

import numpy as np
from math import factorial
from scipy.signal import argrelmax, argrelmin

from .timeseries import Timeseries

__all__ = ["ell2tpat", "tpat2ell", "HmaxPer", "HminPer", "calcHCplane", "getMaxC"]

def ell2tpat(ell,dt,n=5):
    """
    Convert sampling interval to pattern timescale

    Parameters
    ----------
    ell : float
        Sampling interval
    dt : float
        Timestep or timeseries resolution
    n : int, optional
        Sampling size, by default 5

    Returns
    -------
    float
        Equivalent pattern timescale
    """
    return ell*dt*(n-1.)

def tpat2ell(tpat,dt,n=5,returnInt=True):
    """
    Convert pattern timescale to sampling interval

    Parameters
    ----------
    tpat : float or ndarray
        Pattern timescale
    dt : float
        Timestep or timeseries resolution
    n : int, optional
        Sampling size, by default 5
    returnInt : boolean, optional
        Return integer or 

    Returns
    -------
    int, list of int, or ndarray
        Equivalent sampling interval(s); returns int if
        tpat is a single value, list of ints if tpat is
        array/list, and ndarray if tpat is array/list and
        returnInt is set to False
    """    
    if type(tpat)==int or type(tpat)==float:
        return int(np.round(tpat/(dt*(n-1.))))
    elif returnInt:
        return [int(np.round(tp/(dt*(n-1.)))) for tp in tpat]
    else:
        return np.array(tpat)/(dt*(n-1.))

def HmaxPer(n=5): 
    """
    Calculate maximum Permutation Entropy value (H) for a
    periodic function given the specified sampling size n.

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
        Sampling size, by default 5
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
        Sampling, by default 5
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

def tNatApprox(t, data, n=5, method='ncross', attr=None, dt=None, ptcl=None, crossVal=None):
    """
    Approximate the natural timescale for an inputted timeseries.

    Parameters
    ----------
    t : 1-D array
        Timesteps data
    data : 1-D array
        Timeseries data
    n : int, optional
        Sampling size, by default 5, by default 5
    method : str, optional
        Method to use for natural timescale approximation, choose
        from 'minavg', 'maxavg', 'ncross', by default 'ncross'
    attr : str, optional
        Timeseries attribute; needed if extracting a coordinate
        or data attribute from a Timeseries object, by default None
    dt : float, optional
        Timestep resolution; needed if using built-in, 
        by default None
    ptcl: int, optional
        Particle index, by default None
    crossVal: float, optional
        Central value used for 'ncross' method, by default None
    """
    # Check if data is Timeseries object
    if isinstance(data, Timeseries):
        tser = getattr(data, attr) # data timeseries
    else:
        # Setting up data
        tser=np.array(data) # data timeseries

    if ptcl is not None:
        if type(ptcl) is int:
            tser = tser[ptcl]
        else:
            raise TypeError('Particle index (ptcl) must be integer')
    
    if len(tser.shape)>1:
        raise TypeError('Data must be a 1-D array')
    
    if method == 'minavg':
        return np.mean(np.diff(t[argrelmin(tser)]))
    elif method == 'maxavg':
        return np.mean(np.diff(t[argrelmax(tser)]))
    elif method == 'ncross':
        if crossVal is None:
            crossVal = np.mean(tser)
        else:
            pass
        crossLocs = np.where(np.sign(tser[1:]-crossVal) != np.sign(tser[:-1]-crossVal))[0]
        return t[-1]/(len(crossLocs)/2.)
    else:
        raise KeyError("Method {} does not exist. Please use either 'ncross', 'minavg', or 'maxavg'.".format(method))