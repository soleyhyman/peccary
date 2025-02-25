"""
_summary_
"""

import numpy as np
import warnings

from .core import peccary
from . import utils 
from .exceptions import peccaryUserWarning

def doPeccary(timeseries, attr, tIndStart=0, tIndEnd=None, tNat=None, tNatMethod='maxavg', verbose=False):
    """
    Run PECCARY routine and automatically calculate ideal sampling timescales
    and check that they meet the requirements for reliable calculations.

    NOTE: as of 2025-02-25, this routine works for timeseries that have some
    period or dominant repeating patter (e.g., an orbit). Work is underway 
    to make this work for random timeseries.

    Parameters
    ----------
    timeseries : Timeseries
        Timeseries data
    attr : str, optional
        Timeseries attribute; needed if extracting a coordinate
        or data attribute from a Timeseries object, by default None
    tIndStart : int, optional
        Index of corresponding start time of PECCARY analysis window 
        for timeseries, by default 0
    tIndEnd : _type_, optional
        Index of corresponding end time of PECCARY analysis window 
        for timeseries, by default None
    tNat : float, optional
        Natural timescale of timeseries attribute; if not specified,
        will be calculated using the specified or defualt tNatMethod, 
        by default None
    tNatMethod : str, optional
        Method for approximating natural timescale; choose between
        'maxavg', 'minavg', and 'ncross', by default 'maxavg'
    verbose : bool, optional
        Print out timescales and timescale ratios for debugging various
        checks, by default False

    Returns
    -------
    Hvals : ndarray
            Normalized Shannon Perumation Entropy as function of sampling interval, :math:`H(\\ell)`
    Cvals : ndarray
        Normalized Jensen-Shannon complexity as function of sampling interval, :math:`C(\\ell)`
    sampInts : ndarray
        Sampling interval (:math:`\\ell`) values

    Raises
    ------
    RuntimeError
        Duration timescale is less than natural timescale (tDur < tNat). This will result in inaccurate [H,C] values.
    peccaryUserWarning
        Duration timescale < 1.5 * natural timescale. This could result in unreliable [H,C] values.
    """
    # approximate natural timescale
    times = timeseries.t[tIndStart:tIndEnd]
    pecc = peccary(timeseries, attr=str(attr), indStart=tIndStart, indEnd=tIndEnd)
    if tNat==None:
        tNat = utils.tNatApprox(timeseries.t, pecc.tser, method=tNatMethod)
    else:
        pass

    # raise error if tDur/tNat>1 or raise warning if tDur/tNat<1.5
    tDur = times[-1]-times[0]
    if tDur/tNat<1.:
        # warningMsg = 'WARNING: Duration timescale is less than natural timescale. This will result in inaccurate [H,C] values.'
        # warnings.warn(warningMsg, peccaryUserWarning)
        raise RuntimeError('Duration timescale is less than natural timescale (tDur < tNat). This will result in inaccurate [H,C] values.')
    elif tDur/tNat<1.5:
        warningMsg = 'WARNING: Duration timescale < 1.5 * natural timescale. This could result in unreliable [H,C] values.'
        warnings.warn(warningMsg, peccaryUserWarning)
    else:
        pass
    
    # debugging printout if verbose is True
    if verbose==True:
        # print as a check for initial tests
        print('tNat =', tNat)
        print('tDur =', tDur)
        print('tDur/tNat =', tDur/tNat)
    else:
        pass

    # get ideal sampling intervals
    sampRange = utils.calcSampInt(tnat=tNat, dt=pecc.dt)
    sampRangeArr = np.arange(sampRange[0], sampRange[1])
    condition = (tDur/sampRangeArr)>0.3
    sampRangeArr = np.extract(condition, sampRangeArr)

    # debugging printout if verbose is True
    if verbose==True:
        # print as a check for initial tests
        tPats = pecc.tPat(sampRangeArr)
        print('tPat =', tPats)
        print('tDur/tPat =', tDur/tPats)
    else:
        pass

    if len(sampRangeArr)==0:
        raise RuntimeError("The duration timescale of the provided timeseries is too small for any of the appropriate pattern timescales/sampling intervals (tDur < 0.3 tPat).")
    else:
        pass
    
    # calc H, C
    Hvals,Cvals,sampInts = pecc.calcHCcurves(sampIntArray=sampRangeArr)

    return Hvals,Cvals,sampInts