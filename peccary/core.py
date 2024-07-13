"""
PECCARY method.

The class ``peccary`` is the core of the package and does
all of the Permutation Entropy and Statistical Complexity calculations.
"""

import numpy as np
from math import factorial

from .timeseries import Timeseries
from . import utils

__all__ = ["peccary"]

class peccary:
    """
    The ``peccary`` class is the core of the PECCARY method, which is based
    on work by Bandt & Pompe (2002), Rosso et al. (2007), and Weck et al. (2015).
    The class extracts ordinal patterns of a specified sampling size :math:`n` 
    (default :math:`n=5`) and calculates a probability distribution of all possible
    ordinal pattern permutations. After the pattern distribution has been created,
    the values for Permutation Entropy and Jensen-Shanon Statistical Complexity can 
    be calculated. 

    The ``peccary`` class can be initialized with either a 1-D array or a
    ``Timeseries`` instance. In the case of a timeseries object, the desired attribute
    to use for analysis should be indicated with the ``attr`` parameter. If there
    are multiple particles for the data, the desired particle should be indexed with
    the ``ptcl`` parameter. Although the ``dt`` parameter is optional, it is recommended
    to specifiy it if ``data`` is a standard 1-D array, since it is used in the convenience
    methods of ``tPat`` and ``ell_from_tPat``.
    
    See :ref:`Getting started <start>` for a demostration of how to initialize ``peccary`` with a 
    ``Timeseries`` instance. 

    References
    ----------
    [1] `Bandt, C., & Pompe, B. 2002, Phys Rev Lett, 88 (American Physical Society), 174102 <https://link.aps.org/doi/10.1103/PhysRevLett.88.174102>`__

    [2] `Rosso, O. A., Larrondo, H. A., Martin, M. T., Plastino, A., & Fuentes, M. A. 2007, Phys Rev Lett, 99 (American Physical Society), 154102 <https://link.aps.org/doi/10.1103/PhysRevLett.99.154102>`__
    
    [3] `Weck, P. J., Schaffner, D. A., Brown, M. R., & Wicks, R. T. 2015, Phys Rev E, 91 (American Physical Society), 023101 <https://link.aps.org/doi/10.1103/PhysRevE.91.023101>`__
    """
    def __init__(self, data, n=5, attr=None, dt=None, ptcl=None):
        """
        Initialize PECCARY class

        Parameters
        ----------
        data : 1-D array
            Timeseries data for PECCARY 
        n : int, optional
            Sampling size, by default 5, by default 5
        attr : str, optional
            Timeseries attribute; needed if extracting a coordinate
            or data attribute from a Timeseries object, by default None
        dt : float, optional
            Timestep resolution; needed if using built-in, 
            by default None
        ptcl: int, optional
            Particle index, by default None

        Attributes
        ----------
        tser : 1-D array
            Data timeseries
        tlen : int
            Length (number of timesteps of timeseries)
        n : int
            Sampling size
        N : int
            Total number of possible ordinal pattern permutations (:math:`n!`)
        invN : float
            Probability of all ordinal patterns for unifor distribution (:math:`1/n!`)
        log2_N : float
            Calculated value of :math:`\log_2 (n!)`
        log2_Np1 : float
            Calculated value of :math:`\log_2 (n! + 1)`

        Raises
        ------
        TypeError
            Data must be a 1-D array
        """
        # Check if data is Timeseries object
        if isinstance(data, Timeseries):
            # # Check if there's more than one data attribute stored
            # attrCheck = np.array([x is not None for x in [data.x, data.y, data.z, data.data]])
            # if np.sum(attrCheck) > 1:
            #     self.tser = np.where(attrCheck)
            self.tser = getattr(data, attr) # data timeseries
            self.dt = data.dt
        else:
            # Setting up data
            self.tser=np.array(data) # data timeseries

        if ptcl is not None:
            if type(ptcl) is int:
                self.tser = self.tser[ptcl]
            else:
                raise TypeError('Particle index (ptcl) must be integer')
        
        if len(self.tser.shape)>1:
            raise TypeError('Data must be a 1-D array')

        # Precalculate constant values for later use
        self.tlen = len(self.tser) # length of timeseries
        self.n = n # sampling size
        self.N = factorial(n) # n!
        self.invN = 1./self.N # 1/n!
        self.log2_N = np.log2(self.N) # log_2(n!)
        self.log2_Np1 = np.log2(self.N+1.) # log_2(n! + 1)

    def tPat(self, sampInt=1):
        """
        Calculates pattern time of in PECCARY routine, based on 
        sampling size and sampling interval

        Parameters
        ----------
        dt : float
            Timestep interval
        sampInt : int or ndarray, optional
            Sampling interval, by default 1

        Returns
        -------
        float or ndarray
            Pattern time(s) corresponding to inputted sampling interval(s)
        """
        return utils.ell2tpat(ell=sampInt, n=self.n, dt=self.dt)
    
    def ell_from_tPat(self, tPat):
        """
        Calculates the sampling interval in PECCARY routine based 
        on the pattern time and sampling size

        Parameters
        ----------
        dt : float
            Timestep interval
        tPat : float or ndarray
            Pattern time

        Returns
        -------
        float or ndarray
            Sampling interval(s) corresponding to inputted pattern time(s)
        """
        return utils.tpat2ell(tpat=tPat, n=self.n, dt=self.dt)
    
    def getPtot(self,sampInt=1):
        """
        Get total number of permutations based on sampling interval

        Parameters
        ----------
        sampInt : int, optional
            Sampling interval, by default 1

        Returns
        -------
        int
            Total number of permutations
        """
        return self.tlen - sampInt*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator

    def constructPatternCount(self,sampInt=1):
        """
        Count occurance of patterns and total number of permutations

        Parameters
        ----------
        sampInt : int, optional
            Sampling interval, by default 1

        Returns
        -------
        count : ndarray
            A Count occurance of patterns
        Ptot : int 
            Total number of permutations
        """
        Ptot = self.tlen - sampInt*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator
        
        patterns = np.array([self.tser[i:i+(self.n-1)*sampInt+1:sampInt].argsort(kind='stable') for i in range(Ptot)]) # Get patterns based on sampling size and interval
        count = np.unique(patterns, axis=0, return_counts=True) # Count the number of unique ordinal patterns
        
        return count, Ptot

    def getPatternsCounts(self,sampInt=1):
        """
        Return patterns, respective counts, and total number of permutations

        Parameters
        ----------
        sampInt : int, optional
            Sampling interval, by default 1

        Returns
        -------
        patterns : ndarray
            Patterns
        count : ndarray
            Count occurances of patterns
        Ptot : int 
            Total number of permutations
        """
        Ptot = self.tlen - sampInt*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator
        
        patterns = np.array([self.tser[i:i+(self.n-1)*sampInt+1:sampInt].argsort(kind='stable') for i in range(Ptot)]) # Get patterns based on sampling size and interval
        count = np.unique(patterns, axis=0, return_counts=True) # Count the number of unique ordinal patterns
        
        return count[0], count[1], Ptot
    
    def calcS(self, sampInt=1):		
        """
        Calculate Shannon permutation entropy (S) and disqeulibrium (S_e)

        Parameters
        ----------
        sampInt : int, optional
            Sampling interval, by default 1

        Returns
        -------
        Sp : float
            Shannon permutation entropy
        Se : float
            Disequilibrium
        """
        count,Ptot = self.constructPatternCount(sampInt) # Get pattern counts and total number of permutations
        invPtot=1./Ptot # Inverse of total number of permutations

        # Calculate S from the count
        probabilities = count[1]*invPtot # Probability of each pattern occurence; Weck+2015 Eq.(1)
        S = -1.*np.sum([p*np.log2(p) for p in probabilities]) # Shannon entropy formula for pattern probabilities
        Se = -1.*np.sum([P_Pe_sum_2*np.log2(P_Pe_sum_2) for P_Pe_sum_2 in 0.5*(probabilities+self.invN)]) + 0.5*(self.N-len(probabilities))*self.invN*(1+self.log2_N) # Disequilibrium between distribution and uniform distribution; Schaffner & Daniel (in prep), Eq.(8)
        return S, Se

    def calcH(self,sampInt=1):		
        """
        Calculate normalized Permutation Entropy

        Parameters
        ----------
        sampInt : int, optional
            Sampling interval, by default 1

        Returns
        -------
        Sp : float
            Normalized Shannon permutation entropy
        Se : float
            Normalized Shannon + Uniform permutation entropy
        """
        S, Se = self.calcS(sampInt) # Shannon entropy and disequilibrium
        return S/self.log2_N,Se/self.log2_N 

    def calcHC(self,sampInt=1):
        """
        Calculate normalized Jensen-Shannon statistical complexity 
        and normalized permutation entropy

        Parameters
        ----------
        sampInt : int, optional
            Sampling interval, by default 1

        Returns
        -------
        H : float
            Normalized Shannon permutation entropy
        C : float
            Normalized Jensen-Shannon statistical complexity
        """
        S, Se  = self.calcS(sampInt) # Shannon entropy and disequilibrium
        H = S/self.log2_N # Normalized permutation entropy; Schaffner & Daniel (in prep) Eq.(3)
        C = -2.*((Se - 0.5*S - 0.5*self.log2_N)
                /((1. + self.invN)*self.log2_Np1 - 2.*np.log2(2.*self.N) 
                + self.log2_N)*(H)) # Jensen-Shannon statistical complexity; Weck+15 Eq.(4)

        return H, C

    def calcC_fromSSe(self, S, Se):
        """
        Calculate normalized Jensen-Shannon statistical complexity
        from pre-calculated Shannon permutation entropy and 
        disequilibrium values

        Parameters
        ----------
        S : float
            Shannon permutation entropy
        Se : float
            Disequilibrium

        Returns
        -------
        float
            Normalized Jensen-Shannon statistical complexity
        """
        return -2.*((Se - 0.5*S - 0.5*self.log2_N)/((1. + self.invN)*self.log2_Np1 - 2.*np.log2(2.*self.N) + self.log2_N)*(S/self.log2_N)) # Jensen-Shannon statistical complexity; Weck+15 Eq.(4)

    def calcS_fromPatternCount(self, count, Ptot):
        """
        Calculate Shannon permutation entropy and disequilibrium 
        from pattern count and total number of permutations

        Parameters
        ----------
        count : ndarray
            Count occurance result (e.g., from constructPatternCount)
        Ptot : int
            Total number of permutations (e.g., from constructPatternCount)

        Returns
        -------
        Sp : float
            Normalized Shannon permutation entropy
        Se : float
            Normalized Shannon + Uniform permutation entropy
        """
        invPtot=1./Ptot # Inverse of total number of permutations

        probabilities = count[1]*invPtot # Probability of each pattern occurence; Weck+2015 Eq.(1)
        S = -1.*np.sum([p*np.log2(p) for p in probabilities]) # Shannon entropy formula for pattern probabilities
        Se = -1.*np.sum([P_Pe_sum_2*np.log2(P_Pe_sum_2) for P_Pe_sum_2 in 0.5*(probabilities+self.invN)]) + 0.5*(self.N-len(probabilities))*self.invN*(1+self.log2_N) # Disequilibrium between distribution and uniform distribution; Schaffner & Daniel (in prep), Eq.(8)
        return S,Se

    def calcHCcurves(self, min_sampInt=1, max_sampInt=100, step_sampInt=1, sampIntArray=None):
        """
        Returns Permutation Entropy and Statistical Complexity values for multiple
        sampling interval values, i.e., :math:`H(\ell)` and :math:`C(\ell)`

        Parameters
        ----------
        data : array
            Timeseries data for PECCARY
        n : int, optional
            Sampling size, by default 5
        min_sampInt : int, optional
            Smallest sampling interval to loop through, by default 1
        max_sampInt : int, optional
            Largest sampling interval to loop through, by default 100
        step_sampInt : int, optional
            How many sampling interval values to skip over
        sampIntArray: ndarray or list
            Custom array of sampling intervals, supersedes min_sampInt,
            max_sampInt, and step_sampInt, by default None

        Returns
        -------
        Hvals : ndarray
            Normalized Shannon Perumation Entropy as function of sampling interval, :math:`H(\ell)`
        Cvals : ndarray
            Normalized Jensen-Shannon complexity as function of sampling interval, :math:`C(\ell)`
        sampInts : ndarray
            Sampling interval (:math:`\ell`) values
        """
        if type(sampIntArray) != type(None):
            sampInts = sampIntArray # Array of sampling intervals from sampIntArray
        else:
            sampInts = np.arange(min_sampInt,max_sampInt,step_sampInt) # Make array of sampling intervals from min_sampInt to max_sampInt, increasing by 1
        Hvals, Cvals = np.array([self.calcHC(sampInt=int(sampInts[i])) for i in range(len(sampInts))]).T # Get normalized permutation entropy and statisical complexity for each of the sampling intervals
        return Hvals, Cvals, sampInts