"""
Generating and plotting on HC plane.

The class ``HCplots`` is used for plotting the HC plane and
for creating the bases of other plot formats that can be useful for
diagnostics, interpretation, and analysis with PECCARY.
"""

import numpy as np
from math import factorial
import matplotlib.pylab as plt

__all__ = ["HCplots"]

class HCplots:
    def __init__(self, nsteps=1000, n=5):
        """
        Initialize HCplots class

        Parameters
        ----------
        nsteps : int
            Number of steps to use for generating HC bounding curves,
            by default 1000
        n : int, optional
            Embedding dimension/pattern length, by default 5

        Attributes
        ----------
        nsteps : int
            Number of steps used in generating HC bounding curves
        n : int
            Sampling size
        N : int
            Total number of possible ordinal pattern permutations (:math:`n!`)
        invN : float
            Probability of all ordinal patterns for uniform distribution (:math:`1/n!`)
        log2_N : float
            Calculated value of :math:`\log_2 (n!)`
        log2_Np1 : float
            Calculated value of :math:`\log_2 (n! + 1)`

        """
        # Precalculate constant values for later use
        self.nsteps = nsteps
        self.n = n
        self.N = factorial(n)
        self.invN = 1./self.N
        self.log2_N = np.log2(self.N)
        self.log2_Np1 = np.log2(self.N+1.)

    def HmaxPer(self): 
        """
        Calculate maximum Permutation Entropy value (H) for a
        periodic function given the specified sampling size n.
        
        Equation: [insert equation here]

        Returns
        -------
        float
            Maximum Permuation Entropy for periodic function
        """
        Nper = 2.*(2.*(self.n-2.)+1.)
        return np.log(Nper)/np.log(self.N)

    def HminPer(self):
        """
        Calculate minimum Permutation Entropy value (H) for a
        periodic function given the specified sampling size n.
        
        Equation: [insert equation here]

        Returns
        -------
        float
            Minimum Permuation Entropy for periodic function
        """
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

    def getMaxC(self):
        """
        Get maximum possible Statistical Complexity value
        from HC bounding curves for sampling size n

        Returns
        -------
        float
            Maximum possible C value for given n
        """
        Cminx, Cminy, Cmaxx, Cmaxy = self.Cmaxmin()
        return np.max(Cmaxy)

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

    def generateCurves(self, ax, lw=1., ls='--', color='k', alpha=0.5, fontsize=12, showAxLabels=True, showBoundaryLines=True, savePlot=False, savePath=''):
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
        Cminx, Cminy, Cmaxx, Cmaxy = self.Cmaxmin()

        ax.plot(Cminx,Cminy,'k-',Cmaxx,Cmaxy,'k-')

        if showBoundaryLines == True:
            self.getHC_boundary_lines(ax, lw=lw, ls=ls, color=color, alpha=alpha)
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