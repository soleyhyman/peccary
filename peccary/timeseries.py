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
    def __init__(self, s=10, r=20, b=2.667):
        """
        Initialize Lorenz Strange Attractor class.

        Parameters
        ----------
        s : int, optional
            Sigma parameter, by default 10
        r : int, optional
            Rho parameter, by default 20
        b : float, optional
            Beta parameter, by default 2.667

        Notes
        -----
        Modified from `Matplotlib tutorial <https://matplotlib.org/stable/gallery/mplot3d/lorenz_attractor.html>`_

        """
        self.s = s
        self.r = r
        self.b = b

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
    
    def generate(self, dt=0.01, nsteps=10000):
        """
        Generate x/y/z timeseries for Lorenz strange attractor

        Parameters
        ----------
        dt : float, optional
            Timestep resolution, by default 0.01
        nsteps : int, optional
            Number of timesteps to integrate, by default 10000

        Returns
        -------
        lorenzX: ndarray
            x-coordinate timeseries for Lorenz system
        lorenzY: ndarray
            y-coordinate timeseries for Lorenz system
        lorenzZ: ndarray
            z-coordinate timeseries for Lorenz system
        """
        xyzs = np.empty((nsteps + 1, 3))  # Need one more for the initial values
        xyzs[0] = (0., 1., 1.05)  # Set initial values
        # Step through "time", calculating the partial derivatives at the current point
        # and using them to estimate the next point
        for i in range(nsteps):
            xyzs[i + 1] = xyzs[i] + self.getPartials(xyzs[i]) * dt

        return xyzs.T[0], xyzs.T[1], xyzs.T[2]
    
class doublePendulum:
    def __init__(self, tf=2.5, dt=np.power(2.,-6.), L1=1.0, L2=1.0, M1=1.0, M2=1.0):
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
        """
        self.tf = tf
        self.dt = dt
        self.L1 = L1
        self.L2 = L2
        self.L = L1 + L2
        self.M1 = M1
        self.M2 = M2
        self.g = 9.80665 # m/s^2

        # create a time array from 0 to tf, sampled at intervals of dt
        self.t = np.arange(0, self.tf, self.dt)

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
    
    def integrate(self, th1=120.0, w1=0.0, th2=-10.0, w2=0.0):
        """
        Integrate double pendulum system

        Parameters
        ----------
        th1 : float, optional
            Initial angle of mass 1 in degrees, by default 120.0*u.deg
        w1 : float, optional
            Initial angular velocity of mass 1 in degrees/second, by default 0.0
        th2 : float, optional
            Initial angle of mass 2 in degrees, by default -10.0*u.deg
        w2 : float, optional
            Initial angular velocity of mass 2 in degrees/second, by default 0.0
        """
        # initial state
        state = np.radians([th1, w1, th2, w2])

        y = solve_ivp(self.derivs, self.t[[0, -1]], state, t_eval=self.t).y.T

        self.x1 = self.L1*np.sin(y[:, 0])
        self.y1 = -self.L1*np.cos(y[:, 0])

        self.x2 = self.L2*np.sin(y[:, 2]) + self.x1
        self.y2 = -self.L2*np.cos(y[:, 2]) + self.y1

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
        """
        self.N = int(N)
        self.whiteNoiseDat = np.random.random(self.N)
        self.whiteFTT = sfft.rfft(self.whiteNoiseDat)
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
        noiseSpec = self.whiteFTT*freqs
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
        noiseSpec = self.whiteFTT*freqs
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
        noiseSpec = self.whiteFTT*freqs
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
        noiseSpec = self.whiteFTT*freqs
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
        noiseSpec = self.whiteFTT*freqs
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