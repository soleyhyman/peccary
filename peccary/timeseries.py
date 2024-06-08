
"""Timeseries functions for testing PECCARY"""
__all__ = [
    "lorenz",
    "lorenz_mod",
    "generateLorenz"
]

import numpy as np

def lorenz(x, y, z, s=10, r=20, b=2.667):
    """
    Given:
       x, y, z: a point of interest in three dimensional space
       s, r, b: parameters defining the lorenz attractor
    Returns:
       x_dot, y_dot, z_dot: values of the lorenz attractor's partial
           derivatives at the point x, y, z
    """
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot

def lorenz_mod(x, y, z, s=10, r=28, b=2.667):
    """
    Given:
       x, y, z: a point of interest in three dimensional space
       s, r, b: parameters defining the lorenz attractor
    Returns:
       x_dot, y_dot, z_dot: values of the lorenz attractor's partial
           derivatives at the point x, y, z
    """
    x_dot = s*(y - x)
    y_dot = 2.3*r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot

def generateLorenz(dt=0.01,num_steps=10000,s=10,r=20,b=2.667):
    # Need one more for the initial values
    xs = np.empty(num_steps + 1)
    ys = np.empty(num_steps + 1)
    zs = np.empty(num_steps + 1)
    
    # Set initial values
    xs[0], ys[0], zs[0] = (0., 1., 1.05)
    
    # Step through "time", calculating the partial derivatives at the current point
    # and using them to estimate the next point
    for i in range(num_steps):
        x_dot, y_dot, z_dot = lorenz(xs[i], ys[i], zs[i])
        xs[i + 1] = xs[i] + (x_dot * dt)
        ys[i + 1] = ys[i] + (y_dot * dt)
        zs[i + 1] = zs[i] + (z_dot * dt)
        
    return xs,ys,zs

def generateHenon(N):

    X = np.zeros((2,N))
    X[0,0] = 1.
    X[1,0] = 1.
    a = 1.4
    b = 0.3
    for i in range(1,N):
        X[0,i] = 1. - a * X[0,i-1] ** 2. + X[1,i-1]
        X[1,i] = b * X[0,i-1]
    return X[0,:]

def generateTent(N):
    w=0.1847
    X = np.zeros([N])
    X[0] = 0.1

    for i in range(1,N):
        if X[i-1] < w:
            X[i] = X[i-1]/w
        else:
            X[i] = (1 - X[i-1])/(1 - w)
    return X
	#address = '/Users/peterweck/Documents/extimeseries/tent'+str(N)+'.npz'
	#np.savez_compressed(address, x=X)

def generateLogisticMap(N,r=4.):
	
    X = np.zeros([N])
    X[0] = 0.1
    for i in range(1,N):
        X[i] = r * X[i-1] * (1 - X[i-1])
    return X
	#address = '/Users/peterweck/Documents/extimeseries/logis'+str(N)+'.npz'
	#np.savez_compressed(address, x=X)