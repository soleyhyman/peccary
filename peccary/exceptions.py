"""
Warning classes for PECCARY package

This module contains the base parent warning, peccaryWarning,
which is inherited by other warning subclasses.
"""

import warnings

__all__ = ["peccary"]

class peccaryWarning(Warning):
    """
    The base warning class from which all PECCARY warnings inherit.
    """

class peccaryUserWarning(UserWarning, peccaryWarning):
    """
    The main warning class for PECCARY.
    """