"""
Timeseries class wrapper
"""

import numpy as np
from math import factorial
import matplotlib.pylab as plt

__all__ = ["Timeseries"]

class Timeseries:
    def __init__(self, t, x=None, y=None, z=None, data=None, dt=None):
        """
        Create Timeseries object to store data and coordinates

        Parameters
        ----------
        t : ndarray
            Timesteps corresponding to inputted timeseries
        x : ndarray, optional
            x-coordinates of the system, can be defined for multiple
            particles as [timestep, particle], by default None
        y : ndarray, optional
            y-coordinates of the system, can be defined for multiple
            particles as [timestep, particle], by default None
        z : ndarray, optional
            z-coordinates of the system, can be defined for multiple
            particles as [timestep, particle], by default None
        data : 1D-array, optional
            Generic array of data, by default None
        dt : float, optional
            Timestep resolution, by default None

        Attributes
        ----------
        t : 1D-array
            Timesteps corresponding to inputted timeseries
        x : ndarray, optional
            x-coordinate timeseries
        y : ndarray, optional
            y-coordinate timeseries
        z : ndarray, optional
            z-coordinate timeseries
        data : ndarray, optional
            Generic 1D-array timeseries
        dt : float, optional
            Timestep resolution
        """
        if not isinstance(t, np.ndarray) or len(t.shape) != 1:
            raise ValueError("'t' attribute must be 1D-array")
        else:
            self.t = t
        self.x = x
        self.y = y
        self.z = z
        if data is not None:
            if not isinstance(data, np.ndarray) or len(data.shape) != 1:
                raise ValueError("'data' attribute must be 1D-array")
            else:
                self.data = data
        else:
            pass
        if dt is None:
            self.dt = np.mean(np.diff(t))
        else:
            self.dt = dt