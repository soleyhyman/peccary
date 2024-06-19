"""
PECCARY method and plotting.

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
            Embedding dimension, by default 5

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
        self.N=factorial(n) # n!
        self.invN = 1./self.N # 1/n!
        self.log2_N = np.log2(self.N) # log_2(n!)
        self.log2_Np1 = np.log2(self.N+1.) # log_2(n! + 1)

    def tPattern(self, dt, delay=1):
        """
        Calculates pattern time of in PECCARY routine, based on 
        embedding dimension and embed delay

        Parameters
        ----------
        dt : float
            Timestep interval
        delay : int or ndarray, optional
            Embed delay/pattern length, by default 1

        Returns
        -------
        float or ndarray
            Pattern time(s) corresponding to inputted embed delay(s)
        """
        return delay * dt * (self.n - 1.)
    
    def tau_from_tPattern(self, dt, tPattern):
        """
        Calculates the embed delay in PECCARY routine based 
        on the pattern time and embedding dimension

        Parameters
        ----------
        dt : float
            Timestep interval
        tPattern : float or ndarray
            Pattern time

        Returns
        -------
        float or ndarray
            Embed delay(s) corresponding to inputted pattern time(s)
        """
        return tPattern/(dt * (self.n - 1.))
    
    def getPtot(self,delay=1):
        """
        Get total number of permutations based on embed delay

        Parameters
        ----------
        delay : int, optional
            Embed delay/pattern length, by default 1

        Returns
        -------
        int
            Total number of permutations
        """
        return self.t - delay*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator

    def constructPatternCount(self,delay=1):
        """
        Count occurance of patterns and total number of permutations

        Parameters
        ----------
        delay : int, optional
            Embed delay, by default 1

        Returns
        -------
        count : ndarray
            A Count occurance of patterns
        Ptot : int 
            Total number of permutations
        """
        Ptot = self.t - delay*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator
        
        patterns = np.array([self.T[i:i+(self.n-1)*delay+1:delay].argsort(kind='stable') for i in range(Ptot)]) # Get patterns based on n and delay (aka tau)
        count = np.unique(patterns, axis=0, return_counts=True) # Count the number of unique ordinal patterns
        
        return count, Ptot

    def getPatternsCounts(self,delay=1):
        """
        Return patterns, respective counts, and total number of permutations

        Parameters
        ----------
        delay : int, optional
            Embed delay, by default 1

        Returns
        -------
        patterns : ndarray
            Patterns
        count : ndarray
            Count occurances of patterns
        Ptot : int 
            Total number of permutations
        """
        Ptot = self.t - delay*(self.n - 1) # Total number of order n permutations in T; Weck+2015 Eq.(1) denominator
        
        patterns = np.array([self.T[i:i+(self.n-1)*delay+1:delay].argsort(kind='stable') for i in range(Ptot)]) # Get patterns based on n and delay (aka tau)
        count = np.unique(patterns, axis=0, return_counts=True) # Count the number of unique ordinal patterns
        
        return count[0],count[1], Ptot
    
    def calcS(self, delay=1):		
        """
        Calculate Shannon permutation entropy (S) and disqeulibrium (S_e)

        Parameters
        ----------
        delay : int, optional
            Integer delay, by default 1

        Returns
        -------
        Sp : float
            Shannon permutation entropy
        Se : float
            Disequilibrium
        """
        count,Ptot = self.constructPatternCount(delay) # Get pattern counts and total number of permutations
        invPtot=1./Ptot # Inverse of total number of permutations

        #Calculate S from the count
        probabilities = count[1]*invPtot # Probability of each pattern occurence; Weck+2015 Eq.(1)
        S = -1.*np.sum([p*np.log2(p) for p in probabilities]) # Shannon entropy formula for pattern probabilities
        Se = -1.*np.sum([P_Pe_sum_2*np.log2(P_Pe_sum_2) for P_Pe_sum_2 in 0.5*(probabilities+self.invN)]) + 0.5*(self.N-len(probabilities))*self.invN*(1+self.log2_N) # Disequilibrium between distribution and uniform distribution; Schaffner & Daniel (in prep), Eq.(8)
        return S, Se

    def calcH(self,delay=1):		
        """
        Calculate normalized Permutation Entropy

        Parameters
        ----------
        delay : int, optional
            Integer delay, by default 1

        Returns
        -------
        Sp : float
            Normalized Shannon permutation entropy
        Se : float
            Normalized Shannon + Uniform permutation entropy
        """
        S, Se = self.calcS(delay) # Shannon entropy and disequilibrium
        return S/self.log2_N,Se/self.log2_N 

    def calcCofH(self,delay=1):
        """
        Calculate normalized Jensen-Shannon statistical complexity 
        and normalized permutation entropy

        Parameters
        ----------
        delay : int, optional
            Integer delay, by default 1

        Returns
        -------
        H : float
            Normalized Shannon permutation entropy
        C : float
            Normalized Jensen-Shannon statistical complexity
        """
        S, Se  = self.calcS(delay) # Shannon entropy and disequilibrium
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

    def calcPESCcurves(self, min_delay=1, max_delay=100, delayInt=1):
        """
        Returns PE(tau) and SC(tau) for specified range of tau values

        Parameters
        ----------
        data : array
            Timeseries data for PECCARY
        n : int, optional
            Embedding dimension, by default 5
        min_delay : int, optional
            Smallest embedding delay to loop through, by default 1
        max_delay : int, optional
            Largest embedding delay to loop through, by default 100
        delayInt : int, optional
            How many sampling interval values to skip over

        Returns
        -------
        H : array
            Normalized Shannon Perumation Entropy as function of embedding delay tau
        C(tau) : array
            Normalized Jensen-Shannon complexity as function of embedding delay tau
        """
        delays = np.arange(min_delay,max_delay,delayInt) # Make array of embed delay integers from min_delay to max_delay, increasing by 1
        Hvals, Cvals = np.array([self.calcCofH(delay=delays[i]) for i in range(len(delays))]).T # Get normalized permutation entropy and statisical complexity for each of the embed delays
        return Hvals, Cvals