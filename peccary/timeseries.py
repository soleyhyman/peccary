
# """Timeseries functions for testing PECCARY"""
# __all__ = [
#     "lorenz",
#     "lorenz_mod",
#     "generateLorenz"
# ]

import numpy as np

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
        Modified from `Matplotlib tutorial<https://matplotlib.org/stable/gallery/mplot3d/lorenz_attractor.html>`_

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