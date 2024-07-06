"""
PECCARY method.

The class ``peccary`` is the core of the packagke and does
all of the Permutation Entropy and Statistical Complexity calculations.
"""

import numpy as np
from math import factorial
import matplotlib.pylab as plt

__all__ = ["peccary"]

class peccary:
    def __init__(self, data, n=5):
        """
        Initialize PECCARY class

        Parameters
        ----------
        data : 1-D array
            Timeseries data for PECCARY 
        n : int, optional
            Sampling size, by default 5

        Attributes
        ----------
        T : 1-D array
            Data timeseries
        t : int
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
        # Setting up data
        self.T=np.array(data) # data timeseries
        
        if len(self.T.shape)>1:
            raise TypeError( 'Data must be a 1-D array')

        # Precalculate constant values for later use
        self.t = len(self.T) # length of timeseries
        self.n = n # sampling size
        self.N = factorial(n) # n!
        self.invN = 1./self.N # 1/n!
        self.log2_N = np.log2(self.N) # log_2(n!)
        self.log2_Np1 = np.log2(self.N+1.) # log_2(n! + 1)

    def tPat(self, dt, sampInt=1):
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
        return sampInt * dt * (self.n - 1.)
    
    def ell_from_tPat(self, dt, tPat):
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
        return tPat/(dt * (self.n - 1.))
    
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
        return self.t - sampInt*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator

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
        Ptot = self.t - sampInt*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator
        
        patterns = np.array([self.T[i:i+(self.n-1)*sampInt+1:sampInt].argsort(kind='stable') for i in range(Ptot)]) # Get patterns based on sampling size and interval
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
        Ptot = self.t - sampInt*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator
        
        patterns = np.array([self.T[i:i+(self.n-1)*sampInt+1:sampInt].argsort(kind='stable') for i in range(Ptot)]) # Get patterns based on sampling size and interval
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

        #Calculate S from the count
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

    def calcCofH(self,sampInt=1):
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

    def calcPESCcurves(self, min_sampInt=1, max_sampInt=100, step_sampInt=1, sampIntArray=None):
        """
        Returns Permutation Entropy and Statistical Complexity values for multiple
        sampling interval values, i.e., :math:`H(\ell)` and :math:`C(\ell) `

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
        :math:`H(\ell)` : array
            Normalized Shannon Perumation Entropy as function of sampling interval :math:`\ell`
        :math:`C(\ell)` : array
            Normalized Jensen-Shannon complexity as function of sampling interval :math:`\ell`
        """
        if type(sampIntArray) != type(None):
            sampInts = sampIntArray # Array of sampling intervals from sampIntArray
        else:
            sampInts = np.arange(min_sampInt,max_sampInt,step_sampInt) # Make array of sampling intervals from min_sampInt to max_sampInt, increasing by 1
        Hvals, Cvals = np.array([self.calcCofH(sampInt=int(sampInts[i])) for i in range(len(sampInts))]).T # Get normalized permutation entropy and statisical complexity for each of the sampling intervals
        return Hvals, Cvals