"""
Test functions for PECCARY method.

The functions and classes in the ``examples`` module consist of
different functions and physical systems used for testing PECCARY
"""

import numpy as np
from scipy.integrate import solve_ivp
import scipy.fft as sfft
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from .timeseries import Timeseries

__all__ = ["henonMap", "tentMap", "asymmTentMap", "logisticMap", "lorenz", "doublePendulum", "noiseColors"]

def henonMap(n, a=1.4, b=0.3):
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
    peccary.timeseries.Timeseries
        Timeseries for Hénon map, stored in ``x`` attribute
        of ``Timeseries`` class

    References
    ----------
    [1] For more information on the Hénon map, see Wolfram Mathworld's `entry <https://mathworld.wolfram.com/HenonMap.html>`__.

    Examples
    --------
    To create a chaotic timeseries of 10 points using the Hénon map, e.g.,
    
    >>> import peccary.examples as ex
    >>> henon = ex.henonMap(10)

    The data can be called as follows:

    >>> henon.x
    array([ 1.        ,  0.6       ,  0.796     ,  0.2929376 ,  1.11866259,
       -0.6640871 ,  0.71818243,  0.07867346,  1.20678941, -1.01527491])

    The default parameters are :math:`a=1.4` and :math:b=0.3`, which give a
    chaotic timeseries (i.e., the classical Hénon map). Other values for the
    parameters can result in chaotic, intermittent chaotic and periodic, or
    converging periodic behavior.
    """
    X = np.zeros((2,n))
    X[0,0] = 1.
    X[1,0] = 1.
    for i in range(1,n):
        X[0,i] = 1. - a * X[0,i-1] ** 2. + X[1,i-1]
        X[1,i] = b * X[0,i-1]

    return Timeseries(t=np.arange(n), x=X[0,:], dt=1.)

def tentMap(n, mu=2.):
    """
    Generate timeseries from tent map with parameter :math:`\mu`

    Parameters
    ----------
    n : int
        Number of timesteps to generate
    mu : float, optional
        Parameter for changing the map, by default 2.

    Returns
    -------
    peccary.timeseries.Timeseries
        Timeseries for tent map, stored in ``x`` attribute
        of ``Timeseries`` class

    References
    ----------
    [1] For more information on the tent map, see Wolfram Mathworld's `entry <https://mathworld.wolfram.com/TentMap.html>`__.


    Examples
    --------
    To create a timeseries of 10 points using the tent map, e.g.,
    
    >>> import peccary.examples as ex
    >>> tent = ex.tentMap(10)

    The data can be called as follows:

    >>> tent.x
    array([  0.1,   0.2,   0.4,   0.8,   1.6,   3.2,  -4.4,  -8.8, -17.6, -35.2])

    The :math:`\mu` parameter can be changed to produce different types of behavior.
    """
    X = np.zeros([n])
    X[0] = 0.1

    for i in range(1,n):
        if X[i-1] < 0.5:
            X[i] = mu*X[i-1]
        else:
            X[i] = mu*(1 - X[i-1])
    return Timeseries(t=np.arange(n), x=X, dt=1.)

def asymmTentMap(n, a=0.1847):
    """
    Generate timeseries from asymmetric tent map with parameter :math:`a`

    Parameters
    ----------
    n : int
        Number of timesteps to generate
    a : float, optional
        Parameter for changing the map, by default 0.1847

    Returns
    -------
    peccary.timeseries.Timeseries
        Timeseries for asymmetric tent map, stored in ``x`` attribute
        of ``Timeseries`` class

    Examples
    --------
    To create a timeseries of 10 points using the asymmetric tent map, e.g.,
    
    >>> import peccary.examples as ex
    >>> atent = ex.asymmTentMap(10)

    The data can be called as follows:

    >>> atent.x
    array([0.1       , 0.54141852, 0.56246962, 0.53664955, 0.56831896,
           0.52947508, 0.57711875, 0.51868178, 0.5903572 , 0.50244425])

    The :math:`a` parameter can be changed to produce different types of behavior.
    """
    X = np.zeros([n])
    X[0] = 0.1

    for i in range(1,n):
        if X[i-1] < a:
            X[i] = X[i-1]/a
        else:
            X[i] = (1 - X[i-1])/(1 - a)
    return Timeseries(t=np.arange(n), x=X, dt=1.)

def logisticMap(n, r=4.):
    """
    Generate timeseries from logistic map with growth rate parameter :math:`r`

    Parameters
    ----------
    n : int
        Number of timesteps to generate
    r : float, optional
        Growth rate parameter, by default 4.

    Returns
    -------
    peccary.timeseries.Timeseries
        Timeseries for logisitic map with parameter r, 
        stored in ``x`` attribute of ``Timeseries`` class

    References
    ----------
    [1] For more information on the logistic map, see Wolfram Mathworld's `entry <https://mathworld.wolfram.com/LogisticMap.html>`__.

    Examples
    --------
    To create a timeseries of 10 points using the logistic map, e.g.,
    
    >>> import peccary.examples as ex
    >>> logis = ex.logisticMap(10)

    The data can be called as follows:

    >>> logis.x
    array([0.1       , 0.36      , 0.9216    , 0.28901376, 0.82193923,
           0.58542054, 0.97081333, 0.11333925, 0.40197385, 0.9615635 ])

    The parameter :math:`r=4` is used to produce chaotic timeseries, but
    other values may be used for different behavior. 
    """
    X = np.zeros([n])
    X[0] = 0.1
    for i in range(1,n):
        X[i] = r * X[i-1] * (1 - X[i-1])
    return Timeseries(t=np.arange(n), x=X, dt=1.)

class lorenz:
    """
    The ``lorenz`` class can be used to generate the x-, y-, and z-coordinates
    of the Lorenz strange attractor, which is chaotic for certain values of the
    input parameters :math:`\sigma`, :math:`\\rho`, and :math:`\\beta`.

    References
    ----------
    [1] Code modified and expanded from a `Matplotlib tutorial <https://matplotlib.org/stable/gallery/mplot3d/lorenz_attractor.html>`__

    Examples
    --------
    To initialize the class using the default input parameters (:math:`\sigma = 10`, :math:`\\rho = 20`, :math:`\\beta = 2.667`)
    and initial values (:math:`x_0 = 0`, :math:`y_0 = 1`, :math:`z = 1.05`), run the following code:

    >>> import peccary.examples as ex
    >>> lorenzSys = ex.lorenz()

    To integrate the system, use the ``integrate`` method. The timestep 
    resolution can be controlled with the ``dt`` parameter (default 0.01) 
    and either the number of timesteps ``nsteps`` (default 10000) or the
    duration parameter ``tDur`` for the number of "seconds". The ``tDur`` 
    parameter is by default not used, but when specified, it supersedes 
    ``nsteps``. To integrate the initialized system for 250.0 seconds, run:

    >>> lor = lorenzSys.integrate(tDur=250.0)

    The 3D data is stored in a ``peccary.timeseries.Timeseries`` class with
    the attributes ``x``, ``y``, and ``z``, with the timesteps stored in 
    the ``t`` attribute. The data can be extracted and plotted, e.g.,

    >>> import matplotlib.pyplot as plt
    >>> ax = plt.figure().add_subplot(projection='3d')
    >>> ax.plot(lor.x, lor.y, lor.z, lw=0.5)
    >>> plt.show()
    """
    def __init__(self, s=10, r=20, b=2.667, initialVals=(0.,1.,1.05)):
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
        initialVals: 3-tuple or similar three-element array
            Initial values for system, by default (0.,1.,1.05)

        Attributes
        ----------
        s : int
            Sigma parameter
        r : int
            Rho parameter
        b : float
            Beta parameter
        initialVals: 3-tuple
            Initial values for system
        """
        self.s = s
        self.r = r
        self.b = b
        self.initialVals = tuple(initialVals)

    def getPartials(self, xyz):
        """
        Get partial derivatives for initialized Lorenz system

        Parameters
        ----------
        xyz : ndarray
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
    
    def integrate(self, dt=0.01, nsteps=10000, tDur=None):
        """
        Integrate x/y/z timeseries for Lorenz strange attractor

        Parameters
        ----------
        dt : float, optional
            Timestep resolution, by default 0.01
        nsteps : int, optional
            Number of timesteps to integrate, by default 10000
        tDur : float, optional
            Duration of time to integrate over, supersedes nsteps,
            by default None

        Attributes
        ----------
        dt : float
            Timestep resolution
        nsteps : int
            Number of timesteps used in integration
        tDur : None or float
            Duration of timseries, None unless specified
        t : None or ndarray
            Array of timesteps corresponding to timeseries
        t0 : float
            Initial time for integration

        Returns
        -------
        peccary.timeseries.Timeseries
            Timeseries for x-, y-, and z- coordinates of Lorenz strange
            attractor, stored in ``x``, ``y``, and ``z`` attributes of ``Timeseries`` class
        """
        self.dt = dt
        self.nsteps = nsteps 
        self.tDur = tDur
        if self.tDur is not None:
            # check if time duration has been set, if so create times array from that
            self.t = np.arange(0.,self.tDur+self.dt,self.dt)
            self.nsteps = len(self.t[:-1])
        else:
            # if duration has not be set, create "times" array from nsteps and dt
            self.t = np.arange(0.,(self.nsteps+1)*self.dt,self.dt)
        xyzs = np.empty((self.nsteps + 1, 3))  # Need one more for the initial values
        xyzs[0] = self.initialVals  # Set initial values
        
        # Integrate for each timestep
        for i in range(self.nsteps):
            xyzs[i + 1] = xyzs[i] + self.getPartials(xyzs[i]) * self.dt

        x,y,z = xyzs.T[0], xyzs.T[1], xyzs.T[2]

        return Timeseries(t=self.t, x=x, y=y, z=z, dt=self.dt)
    
class doublePendulum:
    """
    The ``doublePendulum`` class can be used to generate the x- and y-timeseries
    of a specified double pendulum system. 
    
    References
    ----------
    [1] Code based on a `Matplotlib tutorial <https://matplotlib.org/stable/gallery/animation/double_pendulum.html>`__
    
    [2] The formulae in that tutorial turn were translated from the C code by `Michael S. Wheatland <http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c>`__
    
    Examples
    --------
    By default, the lengths of the pendulum rods are 1 meter and the two masses 
    are 1 kilogram each, but these can be changed using the parameters ``L1``, 
    ``L2``, ``M1``, and ``M2``. To initialize the class with the defaults, use: 

    >>> import peccary.examples as ex
    >>> pendSys = ex.doublePendulum()

    The ``integrate`` method allows for control of the timestep resolution
    (default :math:`2^{-6}`), the time to integrate over ``tDur``, and the initial
    conditions (:math:`\\theta_1=120^\circ`, :math:`\omega_1=0^\circ ~\\text{s}^{-1}`, :math:`\\theta_2=-10^\circ`, :math:`\omega_1=0^\circ ~\\text{s}^{-1}`).
    To integrate the systems for 10.0 seconds with the default resolution and
    initial conditions, use:

    >>> pend = pendSys.integrate(tDur=10.0)

    The x- and y- coordinates are stored in a ``peccary.timeseries.Timeseries``
    class in the ``x`` and ``y`` attributes in (2,n)-shape arrays. For example,
    to index the x-coordinates of mass 2, use:

    >>> pend.x[1]

    There are also two plotting methods built in to the ``doublePendulum`` class:
    ``plotStatic`` and ``plotAnimate``. The method ``plotStatic`` plots the xy-plane
    of the pendulum system, as well as the x- and y-timeseires of the specified 
    pendulum mass. By default, it plots the lower mass, mass 2, but mass 1 can be 
    selected setting the ``mass`` argument to ``mass=1``. The method returns
    ``matplotlib.figure.Figure`` and ``matplotlib.axes.Axes`` objects, which can
    be plotted or saved to a file. To plot, call the method for the initialized
    system and use the integrated ``Timeseries`` for the parameter ``tser``, i.e.,

    >>> fig,ax = pendSys.plotStatic(pend)
    >>> plt.show() # to show plot
    >>> ### OR ###
    >>> plt.savefig('pend.png') # to save plot
    
    The ``plotAnimate`` method works the same way, except that it animates the
    xy-plane only. It can be used as follows:

    >>> pendSys.plotAnimate(pend)
    """
    def __init__(self, L1=1.0, L2=1.0, M1=1.0, M2=1.0):
        """
        Initialize double pendulum function

        Parameters
        ----------
        L1 : float, optional
            Length of pendulum 1 in m, by default 1.0
        L2 : float, optional
            Length of pendulum 2 in m, by default 1.0
        M1 : float, optional
            Mass of pendulum 1 in kg, by default 1.0
        M2 : float, optional
            Mass of pendulum 2 in kg, by default 1.0

        Attributes
        ----------
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
        g : 9.80665
            Hardcoded value of Earth's gravitational constant in m/s:math:`^2`
        line : ``matplotlib.lines.Line2D`` object
            Lines representing the pendulum rods, 
            *only exists when plotAnimate is used*
        trace : ``matplotlib.lines.Line2D`` object
            Points representing history,
            *only exists when plotAnimate is used* 
        textCurrentTime : Matplotlib Text object
            Text object of current frame timestep,
            *only exists when plotAnimate is used*
        time_template : str
            String format for frame timestep labels,
            *only exists when plotAnimate is used*
        """
        self.L1 = L1
        self.L2 = L2
        self.L = L1 + L2
        self.M1 = M1
        self.M2 = M2
        self.g = 9.80665 # m/s^2
    
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
    
    def integrate(self, tDur=2.5, dt=np.power(2.,-6.), th1=120.0, w1=0.0, th2=-10.0, w2=0.0):
        """
        Integrate double pendulum system

        Parameters
        ----------
        tDur : float, optional
            How many seconds to simulate, by default 2.5
        dt : float, optional
            Time resolution in seconds, by default 0.015625
        th1 : float, optional
            Initial angle of mass 1 in degrees, by default 120.0*u.deg
        w1 : float, optional
            Initial angular velocity of mass 1 in degrees/second, by default 0.0
        th2 : float, optional
            Initial angle of mass 2 in degrees, by default -10.0*u.deg
        w2 : float, optional
            Initial angular velocity of mass 2 in degrees/second, by default 0.0

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
        # create a time array from 0 to tDur, sampled at intervals of dt
        t = np.arange(0, tDur, dt)

        # initial state
        state = np.radians([th1, w1, th2, w2])

        y = solve_ivp(self.derivs, t[[0, -1]], state, t_eval=t).y.T

        x1 = self.L1*np.sin(y[:, 0])
        y1 = -self.L1*np.cos(y[:, 0])

        x2 = self.L2*np.sin(y[:, 2]) + x1
        y2 = -self.L2*np.cos(y[:, 2]) + y1

        return Timeseries(t=t, x=np.vstack((x1,x2)), y=np.vstack((y1,y2)), dt=dt)

    def plotStatic(self, tser, mass=2):
        """
        Quick-plot y(x), x(t), and y(t) for one of the masses
        of a double pendulum system.

        Parameters
        ----------
        tser : peccary.timeseries.Timeseries
            ``Timeseries`` class created by doublePendulum.integrate()
        mass : int, optional
            Pendulum mass to plot, by default 2 (lower mass)

        Returns
        -------
        fig : ``matplotlib.figure.Figure``
            ``matplotlib`` figure created for plot
        ax : array of ``matplotlib.axes.Axes``
            Array of ``matplotlib`` Axes created for plot
        """
        m = int(mass-1)

        fig,axs = plt.subplots(1,3, figsize=(15, 4))
        plt.subplots_adjust(wspace=0.3)
        axs[0].axis('equal')
        axs[0].plot(tser.x[m],tser.y[m])
        axs[0].set_title('Double pendulum Mass {}'.format(int(mass)))
        axs[0].set_xlabel('x (meters)')
        axs[0].set_ylabel('y (meters)')

        axs[1].plot(tser.t,tser.x[m])
        axs[1].set_xlabel('Time (seconds)')
        axs[1].set_ylabel('x (meters)')
        axs[1].set_title('Mass {} x-coordinate'.format(int(mass)))

        axs[2].plot(tser.t,tser.y[m])
        axs[2].set_xlabel('Time (seconds)')
        axs[2].set_ylabel('y (meters)')
        axs[2].set_title('Mass {} y-coordinate'.format(int(mass)))

        return fig,axs
    
    def animateFrame(self, i, tser):
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
        tser : peccary.timeseries.Timeseries
            ``Timeseries`` class created by doublePendulum.integrate()

        Returns
        -------
        ``matplotlib.lines.Line2D`` object
            Lines representing the pendulum rods
        ``matplotlib.lines.Line2D`` object
            Points representing history 
        ``matplotlib.text.Text`` instance
            Text object of current frame timestep
        """
        presentX = [0, tser.x[0][i], tser.x[1][i]]
        presentY = [0, tser.y[0][i], tser.y[1][i]]

        pastX = tser.x[1][:i]
        pastY = tser.y[1][:i]

        self.line.set_data(presentX, presentY)
        self.trace.set_data(pastX, pastY)
        self.textCurrentTime.set_text(self.time_template % (i*tser.dt))
        return self.line, self.trace, self.textCurrentTime

    def plotAnimate(self, tser):
        """
        Animate a integrated double pendulum system

        Parameters
        ----------
        tser : peccary.timeseries.Timeseries
            ``Timeseries`` class created by doublePendulum.integrate()
        """
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(autoscale_on=False, xlim=(-self.L, self.L), ylim=(-self.L, 1.))
        ax.set_aspect('equal')
        ax.grid()

        self.line, = ax.plot([], [], 'o-', lw=2)
        self.trace, = ax.plot([], [], '.-', lw=1, ms=2)
        self.time_template = 'time = %.1fs'
        self.textCurrentTime = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        
        ani = FuncAnimation(fig, self.animateFrame, len(tser.y[1]), 
                            fargs=(tser,), interval=tser.dt*1000, blit=True)
        plt.show()

class noiseColors:
    """
    The ``noiseColors`` class generates noisy timeseries of length ``n``
    with different power spectra. Various methods generate timesereis for
    white noise (flat power density spectrum), blue noise (density :math:`\propto \\nu`), 
    violet noise (density :math:`\propto \\nu^2`), Brownian/red noise (density :math:`\propto \\nu^{-2}`), 
    and pink noise (density :math:`\propto \\nu^{-1}`).

    Examples
    --------
    To initialize the ``noiseColors`` class for 1000 points, use:

    >>> import peccary.examples as ex
    >>> noise = ex.noiseColors(1000)

    To generate the different noise power spectra, use:

    >>> white = noise.white()
    >>> blue = noise.blue()
    >>> violet = noise.violet()
    >>> red = noise.red() # this is equivalent to using noise.brownian()
    >>> pink = noise.pink()

    The timeseries are stored in the ``x`` attribute with all of the methods.
    
    Alternatively, the method ``getNoiseTypes`` can be used to generate a list
    of all available noise colors, e.g.,

    >>> noises = ex.noiseColors(1000).getNoiseTypes()
    """
    def __init__(self, n):
        """
        Initialize noiseColors class to generate noisy timeseries
        with different power spectra.

        Parameters
        ----------
        n : int
            Number of data points to generate

        Attributes
        ----------
        n : int
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
        self.n = int(n)
        self.whiteNoiseDat = np.random.random(self.n)
        self.whiteFFT = sfft.rfft(self.whiteNoiseDat)
        self.whiteFreq = sfft.rfftfreq(self.n)
        self.nonzeroFreq = np.where(self.whiteFreq == 0, np.inf, self.whiteFreq)
        self.tArt = np.arange(self.n) # create artificial "timesteps" for generated timeseries
        
    def white(self):
        """
        Generate white noise (flat frequency power spectrum)

        Returns
        -------
        peccary.timeseries.Timeseries
            Timeseries for white noise, stored in ``x`` attribute
            of ``Timeseries`` class
        """
        freqs = np.power(self.whiteFreq,0.)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return Timeseries(t=self.tArt, x=sfft.irfft(noiseSpec), dt=1.)
      
    def blue(self):
        """
        Generate blue noise (density :math:`\propto \\nu`)

        Returns
        -------
        peccary.timeseries.Timeseries
            Timeseries for blue noise, stored in ``x`` attribute
            of ``Timeseries`` class
        """
        freqs = np.power(self.whiteFreq,0.5)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return Timeseries(t=self.tArt, x=sfft.irfft(noiseSpec), dt=1.)
    
    def violet(self):
        """
        Generate violet noise (density :math:`\propto \\nu^2`)

        Returns
        -------
        peccary.timeseries.Timeseries
            Timeseries for violet noise, stored in ``x`` attribute
            of ``Timeseries`` class
        """
        freqs = np.power(self.whiteFreq,1.)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return Timeseries(t=self.tArt, x=sfft.irfft(noiseSpec), dt=1.)
    
    def brownian(self):
        """
        Generate Brownian noise (density :math:`\propto \\nu^{-2}`)

        Returns
        -------
        peccary.timeseries.Timeseries
            Timeseries for Brownian/red noise noise, stored in 
            ``x`` attribute of ``Timeseries`` class
        """
        freqs = np.power(self.nonzeroFreq,-1.)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return Timeseries(t=self.tArt, x=sfft.irfft(noiseSpec), dt=1.)
    
    def red(self):
        """
        Generate red noise (density :math:`\propto \\nu^{-2}`).
        This is an alias of the Brownian noise function

        Returns
        -------
        peccary.timeseries.Timeseries
            Timeseries for red/Brownian noise noise, stored in 
            ``x`` attribute of ``Timeseries`` class
        """
        return self.brownian()
    
    def pink(self):
        """
        Generate pink noise (density :math:`\propto \\nu^{-1}`)

        Returns
        -------
        peccary.timeseries.Timeseries
            Timeseries for pink noise, stored in ``x`` attribute
            of ``Timeseries`` class
        """
        freqs = np.power(self.nonzeroFreq,-0.5)
        freqs = freqs/np.sqrt(np.mean(freqs**2))
        noiseSpec = self.whiteFFT*freqs
        return Timeseries(t=self.tArt, x=sfft.irfft(noiseSpec), dt=1.)
    
    def getNoiseTypes(self):
        """
        Generate all noise types

        Returns
        -------
        list
            List containing generated white, blue, violet, 
            Brownian/red, and pink noise timeseries
        """
        return [self.white(), self.blue(), self.violet(), self.brownian(), self.pink()]