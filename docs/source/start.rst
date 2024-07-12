.. _start:

Getting started with PECCARY
============================

The following example will walk you through the basic functions of the
``peccary`` package. Additional descriptions and examples for all
functions and class can be found in the :ref:`User Guide <user-guide>`.

First, import ``peccary`` (the core analysis of the package), 
``examples`` (which can generate a variety of different
timeseries), and ``HCplots`` (used for plotting on the :math:`HC`-plane):

>>> from peccary import peccary
>>> import peccary.examples as ex
>>> from peccary import HCplots

To start with PECCARY, you need some sort of timeseries to analyze. 
Let's set up a double pendulum system where each mass is 1 kg and 
both rods are 1 meter long. To initialize it, we use:

>>> pendSys =  ex.doublePendulum(L1=1.0, L2=1.0, M1=1.0, M2=1.0)

Now let's integrate the system for 10 seconds at a time resolution
of :math:`2^{-6}:` seconds (approx 0.016 seconds). Right now, we'll 
just use the default initial angle and angular velocity conditions
of :math:`\theta_1 = 120^{\circ}`, :math:`\theta_2 = -10^{\circ}`, 
and :math:`\omega_1 = \omega_2 = 0^{\circ} \textrm{s}^{-1}`, but 
these values  can be changed with the ``th1``, ``w1``, ``th2``, 
and ``w2`` arguments of the function. We can do this with:

>>> pend = pendSys.integrate(tDur=10., dt=np.power(2.,-6.))

The resulting object of the integration is a ``Timeseries`` object,
which wraps the timesteps, time resolution, and x- and y-coordinates
into one object. You can call the time resolution with the attribute
and print it out, e.g.,:

>>> print(pend.dt)  
0.015625

We can check that everything looks right by plotting using the 
``plotStatic`` method of the ``doublePendulum`` class. To do this, 
we apply this method to the ``pendSys`` class we initialized, and 
then input the ``Timeseries`` class returned by the ``integrate``
method, i.e.,:

>>> fig, ax = pendSys.plotStatic(pend)
>>> plt.show()

This gives us the plot:

.. raw:: html

   <img src="_static/pendEx.png" width="100%"
    style="margin-bottom: 32px;"/>

That looks good! Now, let's run PECCARY on the x-coordinates of the second mass (mass 2)
of our double pendulum system. Since the ``doublePendulum`` class stores x- and y-data
for both masses in a ``Timeseries`` class, we need to specify that we want to analyze
the x-coordinates for the the second particle (or mass) in the system:

>>> pecc = peccary(pend, attr='x', ptcl=1)

To run PECCARY for a range of sampling intervals, we use the ``calcHCcurves`` method.
By default, it will run the analysis for sampling intervals from :math:`\ell = 1` to 
:math:`\ell = 100` at increments of 1. This can be changed via the ``min_sampInt``,
``max_sampInt``, and ``step_sampInt`` parameters or with the ``sampIntArray`` parameter.

Right now, we'll just use the default parameters, to get:

>>> H, C, ells = pecc.calcHCcurves()

To plot :math:`H` and :math:`C` values on the :math:`HC`-plane with the default settings,
we simply run:

>>> HCplots.HCplane(H,C)
>>> plt.show()

This gives us the plot:

.. raw:: html

   <img src="_static/pendEx-HC.png" width="100%"
    style="margin-bottom: 32px;"/>

The documentation for ``HCplots.HCplane`` lists additional arguments for modifying the styles
of the various elements of the plot.

Now supposed we'd like to calculate the :math:`[H, C]` coordinate for the maximum Statistical
Complexity of the x-coordinates of mass 2 of this double pendulum system. One way to do it is
to plot the :math:`H(\ell)` and :math:`C(\ell)` curves to identify what timescale corresponds
to maximum complexity, e.g.,

>>> HCplots.HCcurves(H=H, C=C, sampInts=ells)
>>> plt.show()

This gives us the plot:

.. raw:: html

   <img src="_static/pendEx-HC.png" width="100%"
    style="margin-bottom: 32px;"/>
[continue writeup]