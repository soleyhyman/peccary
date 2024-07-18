|logo|

*******
PECCARY
*******
|PyPI| |License|

PECCARY (Permutation Entropy and statistiCal Complexity Analysis for astRophYsics) 
is a pure-python package for distinguishing between regular, complex, and stochastic
behavior in timeseries. It is based on the work by 
`Bandt & Pompe (2002) <https://ui.adsabs.harvard.edu/#abs/2002PhRvL..88q4102B/abstract>`__ , 
`Rosso et al. (2007) <https://ui.adsabs.harvard.edu/#abs/2007PhRvL..99o4102R/abstract>`__ , 
and `Weck et al. (2015) <https://ui.adsabs.harvard.edu/#abs/2015PhRvE..91b3101W/abstract>`__.
This code is also based on work by collaborator David Schaffner, who wrote the initial 
version of some of the method, called `PESCy <https://github.com/dschaffner/PESCy>`__.

In addition to calculating the Permutation Entropy ($H$) and Statistical Complexity
($C$) values, this package also has plotting tools for the $HC$-plane and visualizing the 
resulting $[H,C]$ values for various timeseries.

A detailed summary of the PECCARY method can be found in Hyman, Daniel, & Schaffner (`arXiv:2407.11970 <https://arxiv.org/abs/2407.11970>`__). 
If you make use of PECCARY, please include a citation to Hyman, Daniel, & Schaffner (`arXiv:2407.11970 <https://arxiv.org/abs/2407.11970>`__)
in any publications.

Documentation
-------------
|Documentation Status|

The documentation for ``peccary`` is hosted on `Read the Docs <http://peccary.readthedocs.io>`__.

Installation and Dependencies
-----------------------------

The recommended way to install the latest stable version of ``peccary`` 
is with ``pip`` via the terminal with the command:

>>> pip install peccary

You can also use the command:

>>> python -m pip install peccary

See the `installation instructions <https://peccary.readthedocs.io/en/latest/installation.html>`__
in the `documentation <https://peccary.readthedocs.io>`__ for more instructions.

.. |PyPI| image:: https://badge.fury.io/py/peccary.svg
   :target: https://pypi.org/project/peccary/
.. |Documentation Status| image:: https://readthedocs.org/projects/peccary/badge/?version=latest
   :target: http://peccary.readthedocs.io/en/latest/?badge=latest
.. |logo| image:: https://peccary.readthedocs.io/en/latest/_static/peccary-logo-banner.png
   :target: https://github.com/soleyhyman/peccary
   :width: 400
.. |License| image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/soleyhyman/peccary/blob/main/LICENSE