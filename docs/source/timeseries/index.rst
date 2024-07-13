.. module:: peccary.timeseries

.. _peccary-timeseries:

*******************************************
Storing timeseries (``peccary.timeseries``)
*******************************************

Introduction
============

The ``timeseries`` module defines the ``Timeseries`` class wrapper,
which can be used in concert with ``peccary.peccary``. It can store
2D or 3D Cartesian coordinates in the ``x``, ``y``, and ``z`` attributes,
or generic 1D arrays in the ``data`` attribute. The timesteps corresponding
to the timeseries data are stored in the ``t`` attribute, and the timesteps
resolution is stored in the ``dt`` attribute.

All members of the ``peccary.examples`` module store their timeseries with
the ``Timeseries`` class.

.. _peccary-timeseries-api:

API
===

.. automodapi:: peccary.timeseries
    :no-inheritance-diagram:
    :automodsumm_included_members: __init__