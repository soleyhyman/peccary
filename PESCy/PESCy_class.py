
"""Defines the PESCy functions."""
__all__ = [
    "calcS",
    "calcH",
    "calcCofH",
    "Cmaxmin",
    "generateCurves",
    "constructPatternCount",
    "calcS_fromPatternCount",
    "calcPESCcurves" 
]

import numpy as np
from math import factorial
from collections import Counter
import matplotlib.pylab as plt

class PESCy:
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
        self.T=np.array(data)
        if len(self.T.shape)>1:
            raise TypeError( 'Data must be a 1-D array')
        self.t = len(self.T)

        self.n = n
        self.N=factorial(n)
        self.invN = 1./self.N
        self.log2_N = np.log2(self.N)
        self.log2_Np1 = np.log2(self.N+1.)

    def constructPatternCount(self,delay=1):
        """
        Count occurance of patterns and total number of permutations

        Parameters
        ----------
        delay : int, optional
            Embed delay, by default 1

        Returns
        -------
        count : Counter object
            A Count occurance of patterns
        Ptot : int 
            Total number of permutations
        """
        Ptot = self.t - delay*(self.n - 1)    #Total number of order n permutations in T
        #print 'Number of permutations = ', Ptot
        A = []			 #Array to store each permutation
        
        for i in range(Ptot):	#Will run through all possible n segments of T
            A.append(''.join(self.T[i:i+(self.n-1)*delay+1:delay].argsort().astype(str)))
        
        #Count occurance of patterns
        count=Counter(A)
        
        return count,Ptot
    
    def calcS(self, delay=1):		
        """
        function Cjs - Returns the Shannon Permutation Energy

        Parameters
        ----------
        delay : int, optional
            Integer delay, by default 1

        Returns
        -------
        Sp : float
            Shannon permutation entropy
        Se : float
            Shannon + Uniform permutation entropy
        """
        count,Ptot = self.constructPatternCount(delay)
        print('Number of permutations = ', Ptot)
        invPtot=1./Ptot     #Inverse for later calcuations

        #Calculate S from the count
        S = 0.
        Se = 0.
        for q in iter(count.values()):
            q*=invPtot #convert to probability
            S += -q * np.log2(q)
            q+=self.invN
            q/=2.
            Se += -q * np.log2(q)
        for i in range(len(count),self.N):
            q=1./2./self.N
            Se += -q * np.log2(q)
        return S,Se

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
        S, Se = self.calcS(delay)
        return S/self.log2_N,Se/self.log2_N

    def calcCofH(self,delay=1):
        '''
        function Cjs - Returns the normalized Jensen-Shannon statistical complexity
        Input:
            data  - array
            n     - embedding dimension (default=5)
            delay - embedding delay (default=1)
        Output:
            C - Normalized Jensen-Shannon complexity
            H - Normalized Shannon Perumation Entropy
        '''		
        S, Se  = self.calcS(delay)
        H = S/self.log2_N
        C = -2.*((Se - 0.5*S - 0.5*self.log2_N)
                /((1. + self.invN)*self.log2_Np1 - 2.*np.log2(2.*self.N) 
                + self.log2_N)*(H))

        return H, C

    def Cmaxmin(self, nsteps):	
        """
        _summary_

        Parameters
        ----------
        nsteps : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        Cmaxx = np.zeros((self.N-1)*nsteps)
        Cmaxy = np.zeros((self.N-1)*nsteps)
        Cminx = np.zeros(nsteps)
        Cminy = np.zeros(nsteps)
        
        for i in np.arange(nsteps):
            pk = self.invN + i*(1.-(self.invN))/nsteps
            pj = (1. - pk)/(self.N - 1.)
            S = -pk * np.log2(pk) - (self.N - 1.) * pj * np.log2(pj)
            qk = pk/2. + 1./(2.*self.N)
            qj = pj/2. + 1./(2.*self.N)
            Scom = -qk * np.log2(qk) - (self.N - 1.) * qj * np.log2(qj)
            Cminx[i] = S / self.log2_N
            Cminy[i] = -2. * (S/self.log2_N) * (Scom - 0.5*S - 0.5*self.log2_N) \
            /((1 + self.invN)*self.log2_Np1 - 2*np.log2(2.*self.N) + self.log2_N)	
            
        for i in np.arange(1,self.N):
            for l in np.arange(nsteps):
                pk = l*(1./(self.N-i+1.))/nsteps
                pj = (1. - pk)/(self.N - i)
                if pk ==0.:
                    S = -(self.N - i) * pj * np.log2(pj)
                else:
                    S = -pk * np.log2(pk) - (self.N - i) * pj * np.log2(pj)
                qk = pk/2. + 1./(2.*self.N)
                qj = pj/2. + 1./(2.*self.N)
                Scom = -qk * np.log2(qk) - (self.N - i) * qj * np.log2(qj) - \
                (i-1)*(1./(2.*self.N))*np.log2(1./(2.*self.N))
                #print (i-1.)*nsteps+l
                Cmaxx[(i-1)*nsteps+l] = S / self.log2_N
                Cmaxy[(i-1)*nsteps+l] = -2.*(S/self.log2_N)*(Scom - 0.5*S - 0.5*self.log2_N) \
                /((1. + self.invN)*self.log2_Np1 - 2.*np.log2(2.*self.N) + self.log2_N)
                
        return Cminx, Cminy, Cmaxx, Cmaxy       

    def calcS_fromPatternCount(self, count, tot_perms):
        """
        Calculate S from pattern count and total number of permutations

        Parameters
        ----------
        count : Counter object
            Count occurance result from constructPatternCount()
        tot_perms : int
            total number of permutations from constructPatternCount()

        Returns
        -------
        Sp : float
            Normalized Shannon permutation entropy
        Se : float
            Normalized Shannon + Uniform permutation entropy
        """
        Ptot=tot_perms
        invPtot=1./Ptot     #Inverse for later calcuations
        S = 0.
        Se = 0.
        for q in iter(count.values()):
            q*=invPtot #convert to probability
            S += -q * np.log2(q)
            q+=self.invN
            q/=2.
            Se += -q * np.log2(q)
        for i in range(len(count),self.N):
            q=1./2./self.N
            Se += -q * np.log2(q)
        return S,Se

    def calcPESCcurves(self, min_delay=1, max_delay=100):
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

        Returns
        -------
        C(tau) : array
            Normalized Jensen-Shannon complexity as function of embedding delay tau
        H : array
            Normalized Shannon Perumation Entropy as function of embedding delay tau
        """
        delay_array = np.arange(min_delay,max_delay)
        num_delays=len(delay_array)
        PEs=np.zeros([num_delays])
        SCs=np.zeros([num_delays])
        for loop_delay in delay_array-1:		
            if (loop_delay%100)==0: print('On Delay ',delay_array[loop_delay])
            permstore_counter = []
            permstore_counter = Counter(permstore_counter)
            tot_perms = 0
            arr,nperms = self.constructPatternCount(delay=delay_array[loop_delay])
            permstore_counter = permstore_counter+arr
            tot_perms = tot_perms+nperms
            PE_tot,PE_tot_Se = self.calcS_fromPatternCount(permstore_counter,tot_perms)
            C =  -2.*((PE_tot_Se - 0.5*PE_tot - 0.5*self.log2_N)
                        /((1 + self.invN)*self.log2_Np1 - 2.*np.log2(2.*self.N) 
                        + self.log2_N)*(PE_tot/self.log2_N))
            PEs[loop_delay]=PE_tot/self.log2_N
            SCs[loop_delay]=C
        return PEs,SCs
    
    def generateCurves(self, nsteps=1000, savePlot=False, savePath=''):
        """
        Creates a blank CH plane with maximum and minimum curves for the given embedding dimension

        Parameters
        ----------
        nsteps : int, optional
            Integer number of steps, by default 1000
        savePlot : bool, optional
            Saves CH plot if set to True, by default False
        savePath : str, optional
            Path to save plot if savePlot set to True, by default ''
            Note: Use only forward slashes in savePath

        Returns
        -------
        fig, ax
            Matplotlib figure and axis objects
        """
        Cminx, Cminy, Cmaxx, Cmaxy = self.Cmaxmin(nsteps)

        fig,ax = plt.subplot()
        ax.plot(Cminx,Cminy,'k-',Cmaxx,Cmaxy,'k-')
        ax.set_xlabel(r"Normalized Permutation Entropy, $H$", fontsize=12)
        ax.set_ylabel(r"Statistical Complexity, $C$", fontsize=12)
        
        if savePlot==True:
            ax.set(xlim=(0,1.0), ylim=(0,0.45))
            ax.set_xticks(np.arange(0,1.1,0.1))
            ax.set_yticks(np.arange(0,0.45,0.05))
            if savePath.endswith('/'):
                savefile=savePath + 'CHn5_blank'
            else:
                savefile=savePath + '/CHn5_blank'
            plt.savefig(savefile+'.png')

        return fig, ax 