"""
Timeseries class wrapper
"""

import numpy as np

__all__ = ["Timeseries"]

class Timeseries:
    """
    The ``Timeseries`` class is a flexible class used to store data
    and timesteps for a simulation. 
    
    The only required parameter is the timesteps parameter, ``t``. 
    If no ``dt`` is inputted in the class initialization, the timestep
    resolution ``dt`` will be calculated from the inputted timesteps
    array.

    The other parameter are ``x``, ``y``, ``z``, and ``data``. The
    first three are intended to be used for systems with one, two, or
    three dimensions and can consist of multidimensional arrays for
    scenarios where the inputted arrays contain data on multiple particles
    (e.g., ``examples.doublePendulum``). The shapes of these arrays should
    be (particles, timesteps).
    
    The ``data`` parameter requires a 1D array and is intended to be the
    catch-all input for datasets that do not fall under the other parameters.

    Examples
    --------
    Suppose we have a random, 1D array of white noise that we'd like to store.
    Since an array of random noise does not have a particular timescale, we
    can approximate it by using the index numbers as the "timesteps". For example:

    >>> from peccary.timeseries import Timeseries
    >>> import numpy as np
    >>> data = np.random.random(1000)
    >>> times = np.arange(len(data))
    >>> tser = Timeseries(t=times, x=data)
    >>> ### OR ###
    >>> tser = Timeseries(t=times, x=data)

    Let's now try making a 
    """
    def __init__(self, t, x=None, y=None, z=None, data=None, dt=None):
        """
        Create Timeseries object to store data and coordinates

        Parameters
        ----------
        t : ndarray
            Timesteps corresponding to inputted timeseries
        x : ndarray, optional
            x-coordinates of the system, can be defined for multiple
            particles as [particle, timestep], by default None
        y : ndarray, optional
            y-coordinates of the system, can be defined for multiple
            particles as [particle, timestep], by default None
        z : ndarray, optional
            z-coordinates of the system, can be defined for multiple
            particles as [particle, timestep], by default None
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

        Raises
        ------
        ValueError:
            If inputted t or data are not a 1D-array
        """
        # Check that timesteps array are 1D array
        if not isinstance(t, np.ndarray) or len(t.shape) != 1:
            raise ValueError("'t' attribute must be 1D-array")
        else:
            self.t = t

        # Set x/y/z attributes
        self.x = x
        self.y = y
        self.z = z

        # Check if data parameter is being used
        if data is not None:
            # Check that data is a 1D array
            if not isinstance(data, np.ndarray) or len(data.shape) != 1:
                raise ValueError("'data' attribute must be 1D-array")
            else:
                self.data = data
        else:
            pass

        # Check if dt has been inputted, if not estimate it from t
        if dt is None:
            self.dt = np.mean(np.diff(t))
        else:
            self.dt = dt