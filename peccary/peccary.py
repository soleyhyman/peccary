"""
Within the ``peccary`` module, there are two main classes.
The first class, ``peccary``, is the core of the packagke and does
all of the Permutation Entropy and Statistical Complexity calculations.
The second class, ``HCplots``, is used for plotting the HC plane and
for creating the bases of other plot formats that can be useful for
diagnostics, interpretation, and analysis with PECCARY.
"""

import numpy as np
from math import factorial
import matplotlib.pylab as plt

class peccary:
    def __init__(self, data, n=5):
        """
        Initialize PESCy class

        Parameters
        ----------
        data : 1-D array
            Timeseries data for PESCy 
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
        Calculates pattern time of in PESCy routine, based on 
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
        Calculates the embed delay in PESCy routine based 
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
            Timeseries data for PESCy
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
    
class HCplots:
    def __init__(self, nsteps, n=5):
        """
        Initialize HCplots class

        Parameters
        ----------
        nsteps : int
            Number of timesteps
        n : int, optional
            Embedding dimension/pattern length, by default 5
        """
        # Precalculate constant values for later use
        self.nsteps = nsteps
        self.n = n
        self.N=factorial(n)
        self.invN = 1./self.N
        self.log2_N = np.log2(self.N)
        self.log2_Np1 = np.log2(self.N+1.)

    def HmaxPer(self): # eq. 6 in PESCy paper (updated)
        Nper = 2.*(2.*(self.n-2.)+1.)
        return np.log(Nper)/np.log(self.N)

    def HminPer(self): # eq. 6 in PESCy paper (updated)
        return np.log(2)/np.log(self.N)

    def Cmaxmin(self):	
        """
        Get maximum and minimum C(H) curves based on pattern length
        for plotting on the HC plane

        Returns
        -------
        Cminx
            H values for minimum HC curve
        Cminy
            C values for minimum HC curve
        Cmaxx
            H values for maximum HC curve
        Cmaxy
            C values for maximum HC curve
        """
        # Set up blank arrays for x- and y-values of max/min curves
        Cmaxx = np.zeros((self.N-1)*self.nsteps)
        Cmaxy = np.zeros((self.N-1)*self.nsteps)
        Cminx = np.zeros(self.nsteps)
        Cminy = np.zeros(self.nsteps)
        
        # Calculate H and C values for minimum HC curve
        for i in np.arange(self.nsteps):
            pk = self.invN + i*(1.-(self.invN))/self.nsteps
            pj = (1. - pk)/(self.N - 1.)
            S = -pk * np.log2(pk) - (self.N - 1.) * pj * np.log2(pj)
            qk = pk/2. + 1./(2.*self.N)
            qj = pj/2. + 1./(2.*self.N)
            Scom = -qk * np.log2(qk) - (self.N - 1.) * qj * np.log2(qj)
            Cminx[i] = S / self.log2_N
            Cminy[i] = -2. * (S/self.log2_N) * (Scom - 0.5*S - 0.5*self.log2_N)\
            /((1 + self.invN)*self.log2_Np1 - 2*np.log2(2.*self.N) + self.log2_N)	
            
        # Calculate H and C values for maximum HC curve
        for i in np.arange(1,self.N):
            for l in np.arange(self.nsteps):
                pk = l*(1./(self.N-i+1.))/self.nsteps
                pj = (1. - pk)/(self.N - i)
                if pk ==0.:
                    S = -(self.N - i) * pj * np.log2(pj)
                else:
                    S = -pk * np.log2(pk) - (self.N - i) * pj * np.log2(pj)
                qk = pk/2. + 1./(2.*self.N)
                qj = pj/2. + 1./(2.*self.N)
                Scom = -qk * np.log2(qk) - (self.N - i) * qj * np.log2(qj) - \
                (i-1)*(1./(2.*self.N))*np.log2(1./(2.*self.N))
                Cmaxx[(i-1)*self.nsteps+l] = S / self.log2_N
                Cmaxy[(i-1)*self.nsteps+l] = -2.*(S/self.log2_N)*(Scom - 0.5*S - 0.5*self.log2_N) \
                /((1. + self.invN)*self.log2_Np1 - 2.*np.log2(2.*self.N) + self.log2_N)
                
        return Cminx, Cminy, Cmaxx, Cmaxy     

    def getHC_boundary_lines(self, ax, lw=1., ls='--', color='k', alpha=0.5):
        """
        Creates boundary lines on HC plane

        Parameters
        ----------
        ax : Matplotlib axis
            Axis on which to plot empty HC curves
        lw : float, optional
            Line width of boundary lines (must be float >0), by default 1.
        ls : string, optional
            Line style (Matplotlib linestyle string) of boundary lines, by default '--'
        color : string, optional
            Color of boundary lines (Matplotlib color string), by default 'k'
        alpha : float, optional
            Opacity of boundary lines (float between 0 and 1), by default 0.5
        """
        minH, minC, maxH, maxC = self.Cmaxmin() # H and C values for HC min/max curves
        minHVal = self.HminPer() # Min. H value for periodic fnct
        maxHVal = self.HmaxPer() # Max. H value for periodic fnct
        argMinH_CminCurve = np.argmin(np.abs(minH-minHVal)) # Min. H value argument for H values for min HC curve
        argMinH_CmaxCurve = np.argmin(np.abs(maxH-minHVal)) # Min. H value argument for H values for max HC curve
        argMaxH_CminCurve = np.argmin(np.abs(minH-maxHVal)) # Max. H value argument for H values for min HC curve
        argMaxH_CmaxCurve = np.argmin(np.abs(maxH-maxHVal)) # Max. H value argument for H values for max HC curve
        ax.plot(np.array([minHVal,minHVal]), np.array([minC[argMinH_CminCurve], maxC[argMinH_CmaxCurve]]), lw=lw, ls=ls, color=color, alpha=alpha)
        ax.plot(np.array([maxHVal,maxHVal]), np.array([minC[argMaxH_CminCurve], maxC[argMaxH_CmaxCurve]]), lw=lw, ls=ls, color=color, alpha=alpha)
        ax.plot(np.array([minHVal,maxHVal]), np.array([maxC[argMinH_CmaxCurve], maxC[argMaxH_CmaxCurve]]), lw=lw, ls=ls, color=color, alpha=alpha)  

    def generateCurves(self, ax, lw=1., ls='--', color='k', alpha=0.5, fontsize=12, showAxLabels=True, savePlot=False, savePath=''):
        """
        Creates a blank HC plane with maximum and minimum curves for the given embedding dimension

        Parameters
        ----------
        ax : Matplotlib axis
            Axis on which to plot empty HC curves
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
        savePlot : bool, optional
            Saves HC plot if set to True, by default False
        savePath : str, optional
            Path to save plot if savePlot set to True, by default ''
            Note: Use only forward slashes in savePath
        """
        Cminx, Cminy, Cmaxx, Cmaxy = self.Cmaxmin()

        ax.plot(Cminx,Cminy,'k-',Cmaxx,Cmaxy,'k-')
        self.getHC_boundary_lines(ax, lw=lw, ls=ls, color=color, alpha=alpha)
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