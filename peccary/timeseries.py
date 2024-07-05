"""
Test functions for PECCARY method.

The functions and classes in the ``timeseries`` module consist of
different functions and physical systems used for testing PECCARY
"""

import numpy as np
from scipy.integrate import solve_ivp
import scipy.fft as sfft
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

__all__ = ["generateHenon", "generateTent", "generateAsymmTent", "generateLogisticMap", "lorenz", "doublePendulum", "noiseColors"]

def generateHenon(n, a=1.4, b=0.3):
    """
    Generate timeseries from Hénon map. 
    
    Default parameters are for the classical Hénon map, 
    where values are chaotic.

    Parameters
    ----------
    n : int
        Number of timesteps to generate
    a : float, optional
        Parameter for Hénon map, by default 1.4
    b : float, optional
        Parameter for Hénon map, by default 0.3

    Returns
    -------
    ndarray
        1D array (length n) of timeseries for Hénon map
    """
    X = np.zeros((2,n))
    X[0,0] = 1.
    X[1,0] = 1.
    for i in range(1,n):
        X[0,i] = 1. - a * X[0,i-1] ** 2. + X[1,i-1]
        X[1,i] = b * X[0,i-1]
    return X[0,:]

def generateTent(n, mu=2.):
    """
    Generate timeseries from tent map with parameter mu

    Parameters
    ----------
    n : int
        Number of timesteps to generate
    mu : float, optional
        Parameter for changing the map, by default 2.

    Returns
    -------
    ndarray
        1D array (length n) of timeseries for tent map
    """
    X = np.zeros([n])
    X[0] = 0.1

    for i in range(1,n):
        if X[i-1] < mu:
            X[i] = mu*X[i-1]
        else:
            X[i] = mu*(1 - X[i-1])
    return X

def generateAsymmTent(n, a=0.1847):
    """
    Generate timeseries from asymmetric tent map with parameter a

    Parameters
    ----------
    n : int
        Number of timesteps to generate
    a : float, optional
        Parameter for changing the map, by default 0.1847

    Returns
    -------
    ndarray
        1D array (length n) of timeseries for asymmetric tent map
    """
    X = np.zeros([n])
    X[0] = 0.1

    for i in range(1,n):
        if X[i-1] < a:
            X[i] = X[i-1]/a
        else:
            X[i] = (1 - X[i-1])/(1 - a)
    return X

def generateLogisticMap(n, r=4.):
    """
    Generate timeseries from logistic map with growth rate parameter r

    Parameters
    ----------
    n : int
        Number of timesteps to generate
    r : float, optional
        Growth rate parameter, by default 4.

    Returns
    -------
    ndarray
        1D array (length n) of timeseries for logisitic map with parameter r
    """
    X = np.zeros([n])
    X[0] = 0.1
    for i in range(1,n):
        X[i] = r * X[i-1] * (1 - X[i-1])
    return X

class lorenz:
    def __init__(self, dt=0.01, nsteps=10000, tDur=None, t0=0., s=10, r=20, b=2.667, initialVals=(0.,1.,1.05)):
        """
        Initialize Lorenz Strange Attractor class.

        Parameters
        ----------
        dt : float, optional
            Timestep resolution, by default 0.01
        nsteps : int, optional
            Number of timesteps to integrate, by default 10000
        tDur : float, optional
            Duration of time to integrate over, supersedes nsteps,
            by default None
        t0 : float, optional
            Initial time value, only used if tDur is not None,
            by default 0.
        s : int, optional
            Sigma parameter, by default 10
        r : int, optional
            Rho parameter, by default 20
        b : float, optional
            Beta parameter, by default 2.667
        initialVals: 3-tuple or similar three-element array
            Initial values for system, by default (0.,1.,1.05)

        Attributes
        ----------
        dt : float
            Timestep resolution
        nsteps : int
            Number of timesteps used in integration
        tDur : None or float
            Duration of timseries, None unless specified
        t : None or ndarray
            Array of timesteps corresponding to timeseries,
            None unless tDur has been specified
        t0 : float
            Initial time for integration
        s : int
            Sigma parameter
        r : int
            Rho parameter
        b : float
            Beta parameter
        initialVals: 3-tuple
            Initial values for system
        x : ndarray
            x-coordinates of system after integration
        y : ndarray
            y-coordinates of system after integration
        z : ndarray
            z-coordinates of system after integration
        xyzSim : 3-tuple of ndarrays
            x-, y-, and z- coordinates of system after integration

        Notes
        -----
        Modified from `Matplotlib tutorial <https://matplotlib.org/stable/gallery/mplot3d/lorenz_attractor.html>`_

        """
        self.dt = dt
        self.nsteps = nsteps 
        self.tDur = tDur
        self.t0 = t0
        self.t = None
        self.s = s
        self.r = r
        self.b = b
        self.initialVals = tuple(initialVals)

        self.xyzSim = self.generate()
        self.x, self.y, self.z = self.xyzSim

    def getPartials(self, xyz):
        """
        Get partial derivatives for initialized Lorenz system

        Parameters
        ----------
        xyz : ndarray,
            3D array containing points of interest in three-dimensional space

        Returns
        -------
        ndarray
            3D array containing partial derivatives at xyz
        """
        x, y, z = xyz
        x_dot = self.s*(y - x)
        y_dot = self.r*x - y - x*z
        z_dot = x*y - self.b*z
        return np.array([x_dot, y_dot, z_dot])
    
    def generate(self):
        """
        Generate x/y/z timeseries for Lorenz strange attractor

        Returns
        -------
        lorenzX: ndarray
            x-coordinate timeseries for Lorenz system
        lorenzY: ndarray
            y-coordinate timeseries for Lorenz system
        lorenzZ: ndarray
            z-coordinate timeseries for Lorenz system
        """
        if type(self.tDur) != type(None):
            self.t = np.arange(self.t0,self.tDur+self.dt,self.dt)
            self.nsteps = len(self.t[:-1])
        else:
            pass
        xyzs = np.empty((self.nsteps + 1, 3))  # Need one more for the initial values
        xyzs[0] = self.initialVals  # Set initial values
        # Step through "time", calculating the partial derivatives at the current point
        # and using them to estimate the next point
        for i in range(self.nsteps):
            xyzs[i + 1] = xyzs[i] + self.getPartials(xyzs[i]) * self.dt

        return xyzs.T[0], xyzs.T[1], xyzs.T[2]
    
class doublePendulum:
    def __init__(self, tf=2.5, dt=np.power(2.,-6.), L1=1.0, L2=1.0, M1=1.0, M2=1.0, th1=120.0, w1=0.0, th2=-10.0, w2=0.0):
        """
        Initialize double pendulum function

        Notes
        -----
        Code based on a `Matplotlib tutorial <https://matplotlib.org/stable/gallery/animation/double_pendulum.html>`_
        and translated into a Python class by Sóley Hyman.
        The formulae in that tutorial turn were translated from the C code by `Michael S. Wheatland <http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c>`_

        Parameters
        ----------
        tf : float, optional
            How many seconds to simulate, by default 2.5
        dt : float, optional
            Time resolution in seconds, by default 0.015625
        L1 : float, optional
            Length of pendulum 1 in m, by default 1.0
        L2 : float, optional
            Length of pendulum 2 in m, by default 1.0
        M1 : float, optional
            Mass of pendulum 1 in kg, by default 1.0
        M2 : float, optional
            Mass of pendulum 2 in kg, by default 1.0
        th1 : float, optional
            Initial angle of mass 1 in degrees, by default 120.0*u.deg
        w1 : float, optional
            Initial angular velocity of mass 1 in degrees/second, by default 0.0
        th2 : float, optional
            Initial angle of mass 2 in degrees, by default -10.0*u.deg
        w2 : float, optional
            Initial angular velocity of mass 2 in degrees/second, by default 0.0

        Attributes
        ----------
        tf : float
            Final timestep in seconds
        dt : float
            Time resolution in seconds
        L1 : float
            Length of pendulum 1 in m
        L2 : float
            Length of pendulum 2 in m
        L : float
            Combined length of both pendulums in m
        M1 : float
            Mass of pendulum 1 in kg
        M2 : float
            Mass of pendulum 2 in kg
        th1 : float
            Initial angle of mass 1 in degrees
        w1 : float
            Initial angular velocity of mass 1 in degrees/second
        th2 : float
            Initial angle of mass 2 in degrees
        w2 : float
            Initial angular velocity of mass 2 in degrees/second
        g : 9.80665
            Hardcoded value of Earth's gravitational constant in m/s:math:`^2`
        t : ndarray
            Array of timesteps used in integration
        x1 : ndarray
            x-coordinates of mass 1 after integration
        y1 : ndarray
            y-coordinates of mass 1 after integration
        x2 : ndarray
            x-coordinates of mass 1 after integration
        y2 : ndarray
            y-coordinates of mass 1 after integration
        simPos : 4-tuple
            Timeseries of x- and y- coordinates for masses 1 and 2, respectively
        line : Matplotlib Line2D object
            Lines representing the pendulum rods, 
            *only exists when plotAnimate is used*
        trace : Matplotlib Line2D object
            Points representing history,
            *only exists when plotAnimate is used* 
        textCurrentTime : Matplotlib Text object
            Text object of current frame timestep,
            *only exists when plotAnimate is used*
        time_template : str
            String format for frame timestep labels,
            *only exists when plotAnimate is used*
        """
        self.tf = tf
        self.dt = dt
        self.L1 = L1
        self.L2 = L2
        self.L = L1 + L2
        self.M1 = M1
        self.M2 = M2
        self.th1 = th1
        self.w1 = w1
        self.th2 = th2
        self.w2 = w2
        self.g = 9.80665 # m/s^2

        # create a time array from 0 to tf, sampled at intervals of dt
        self.t = np.arange(0, self.tf, self.dt)
        self.simPos = self.integrate()
        self.x1, self.y1, self.x2, self.y2 = self.simPos

    def derivs(self, t, state):
        """
        Get partial derivatives for double pendulum system

        Parameters
        ----------
        t : ndarray
            Time array
        state : _type_
            Initial conditions

        Returns
        -------
        ndarray
            Partial derivatives
        """
        dydx = np.zeros_like(state)

        dydx[0] = state[1]

        delta = state[2] - state[0]
        den1 = (self.M1+self.M2) * self.L1 - self.M2 * self.L1 * np.cos(delta) * np.cos(delta)
        dydx[1] = ((self.M2 * self.L1 * state[1] * state[1] * np.sin(delta) * np.cos(delta)
                    + self.M2 * self.g * np.sin(state[2]) * np.cos(delta)
                    + self.M2 * self.L2 * state[3] * state[3] * np.sin(delta)
                    - (self.M1+self.M2) * self.g * np.sin(state[0]))
                / den1)

        dydx[2] = state[3]

        den2 = (self.L2/self.L1) * den1
        dydx[3] = ((- self.M2 * self.L2 * state[3] * state[3] * np.sin(delta) * np.cos(delta)
                    + (self.M1+self.M2) * self.g * np.sin(state[0]) * np.cos(delta)
                    - (self.M1+self.M2) * self.L1 * state[1] * state[1] * np.sin(delta)
                    - (self.M1+self.M2) * self.g * np.sin(state[2]))
                / den2)

        return dydx
    
    def integrate(self):
        """
        Integrate double pendulum system

        Returns
        -------
        ndarray
            Timeseries of x-coordinates for mass 1
        ndarray
            Timeseries of y-coordinates for mass 1
        ndarray
            Timeseries of x-coordinates for mass 2
        ndarray
            Timeseries of y-coordinates for mass 2
        """
        # initial state
        state = np.radians([self.th1, self.w1, self.th2, self.w2])

        y = solve_ivp(self.derivs, self.t[[0, -1]], state, t_eval=self.t).y.T

        x1 = self.L1*np.sin(y[:, 0])
        y1 = -self.L1*np.cos(y[:, 0])

        x2 = self.L2*np.sin(y[:, 2]) + x1
        y2 = -self.L2*np.cos(y[:, 2]) + y1

        return x1, y1, x2, y2

    def plotStatic(self):
        """
        Quick-plot y(x), x(t), and y(t) for the second mass (lower mass)
        of the double pendulum system.
        """
        fig,axs = plt.subplots(1,3, figsize=(15, 4))
        axs[0].set_aspect('equal')
        axs[0].plot(self.x2,self.y2)
        axs[0].set_title('XY plane of pendulum movment')

        axs[1].plot(self.t,self.x2)
        axs[1].set_xlabel('Time (seconds)')
        axs[1].set_title('X coordinate')

        axs[2].plot(self.t,self.y2)
        axs[2].set_xlabel('Time (seconds)')
        axs[2].set_title('Y coordinate')

        return fig,axs
    
    def animateFrame(self, i):
        """
        Animate a frame of the double pendulum system

        Notes
        -----
        This will very rarely be used alone; please use 
        the plotAnimate function instead, which is a wrapper
        for this function        

        Parameters
        ----------
        i : int
            Index of the frame to animate

        Returns
        -------
        Matplotlib Line2D object
            Lines representing the pendulum rods
        Matplotlib Line2D object
            Points representing history 
        Matplotlib Text object
            Text object of current frame timestep
        """
        presentX = [0, self.x1[i], self.x2[i]]
        presentY = [0, self.y1[i], self.y2[i]]

        pastX = self.x2[:i]
        pastY = self.y2[:i]

        self.line.set_data(presentX, presentY)
        self.trace.set_data(pastX, pastY)
        self.textCurrentTime.set_text(self.time_template % (i*self.dt))
        return self.line, self.trace, self.textCurrentTime

    def plotAnimate(self):
        """
        Animate the integrated double pendulum system
        """
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(autoscale_on=False, xlim=(-self.L, self.L), ylim=(-self.L, 1.))
        ax.set_aspect('equal')
        ax.grid()

        self.line, = ax.plot([], [], 'o-', lw=2)
        self.trace, = ax.plot([], [], '.-', lw=1, ms=2)
        self.time_template = 'time = %.1fs'
        self.textCurrentTime = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        
        ani = FuncAnimation(fig, self.animateFrame, len(self.y2), interval=self.dt*1000, blit=True)
        plt.show()

class noiseColors:
    def __init__(self, N):
        """
        Initialize noiseColors class to generate noisy timeseries
        with different power spectra.

        Parameters
        ----------
        N : int
            Number of data points to generate

        Attributes
        ----------
        N : int
            Number of data points used in generating timeseries
        whiteNoiseDat : ndarray
            Default white noise generated with specified size
        whiteFTT : complex ndarray
            1-D discrete Fourier transform of whiteNoiseDat
        whiteFreq : ndarray
            Frequencies corresponding to whiteFFT
        nonzeroFreq : ndarray
            whiteFreq with all frequencies equal to 0 replaced with inf
            (needed for brownianNoise and pinkNoise)
        """
        self.N = int(N)
        self.whiteNoiseDat = np.random.random(self.N)
        self.whiteFFT = sfft.rfft(self.whiteNoiseDat)
        self.whiteFreq = sfft.rfftfreq(self.N)
        self.nonzeroFreq = np.where(self.whiteFreq == 0, np.inf, self.whiteFreq)
        
    def whiteNoise(self):
        """
        Generate white noise (flat frequency power spectrum)

        Returns
        -------
        ndarray
            1D array of length N containing white noise
        """
        freqs = np.power(self.whiteFreq,0.)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return sfft.irfft(noiseSpec)
      
    def blueNoise(self):
        """
        Generate blue noise (density :math:`\propto f`)

        Returns
        -------
        ndarray
            1D array of length N containing blue noise
        """
        freqs = np.power(self.whiteFreq,0.5)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return sfft.irfft(noiseSpec)
    
    def violetNoise(self):
        """
        Generate violet noise (density :math:`\propto f^2`)

        Returns
        -------
        ndarray
            1D array of length N containing violet noise
        """
        freqs = np.power(self.whiteFreq,1.)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return sfft.irfft(noiseSpec)
    
    def brownianNoise(self):
        """
        Generate Brownian noise (density :math:`\propto 1/f^2`)

        Returns
        -------
        ndarray
            1D array of length N containing Brownian/red noise
        """
        freqs = np.power(self.nonzeroFreq,-1.)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return sfft.irfft(noiseSpec)
    
    def redNoise(self):
        """
        Generate red noise (density :math:`\propto 1/f^2`).
        This is an alias of the Brownian noise function

        Returns
        -------
        ndarray
            1D array of length N containing red/Brownian noise
        """
        return self.brownianNoise()
    
    def pinkNoise(self):
        """
        Generate pink noise (density :math:`\propto 1/f`)

        Returns
        -------
        ndarray
            1D array of length N containing pink noise
        """
        freqs = np.power(self.nonzeroFreq,-0.5)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return sfft.irfft(noiseSpec)
    
    def getNoiseTypes(self):
        """
        Generate all noise types

        Returns
        -------
        list
            List containing generated white, blue, violet, 
            Brownian/red, and pink noise timeseries
        """
        return [self.whiteNoise(), self.blueNoise(), self.violetNoise(), self.brownianNoise(), self.pinkNoise()]