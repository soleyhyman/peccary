.. PECCARY documentation master file, created by
   sphinx-quickstart on Fri Jun  7 16:52:07 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

   <img src="_static/peccary-logo-banner.png" width="70%"
    style="margin-bottom: 32px;"/>


*******
PECCARY
*******

``peccary`` (Permutation Entropy and statistiCal Complexity Analysis for astRophYsics) 
is a Python package for distinguishing between regular, complex, and stochastic
behavior in timeseries. It is based on the work by 
`Bandt & Pompe <https://ui.adsabs.harvard.edu/#abs/2002PhRvL..88q4102B/abstract>`_ (2002), 
`Rosso et al. <https://ui.adsabs.harvard.edu/#abs/2007PhRvL..99o4102R/abstract>`_ (2007), 
and 
`Weck et al. <https://ui.adsabs.harvard.edu/#abs/2015PhRvE..91b3101W/abstract>`_ (2015).
This code is also based on work by collaborator David Schaffner, who wrote the initial 
version of some of the method, called `PESCy <https://github.com/dschaffner/PESCy>`__.

In addition to calculating the Permutation Entropy (:math:`H`) and Statistical Complexity
(:math:`C`) values, this package also has plotting tools for the :math:`HC`-plane and 
visualizing the resulting :math:`[H,C]` values for various timeseries, examples timeseries, 
and utility functions.

A detailed summary of the PECCARY method can be found in Hyman, Daniel, & Schaffner (`arXiv:2407.11970 <https://arxiv.org/abs/2407.11970>`__). 
If you make use of PECCARY, please include a :ref:`citation <cite>` to 
Hyman, Daniel, & Schaffner (`arXiv:2407.11970 <https://arxiv.org/abs/2407.11970>`__) in any publications.

.. note::
   This project is under active development.


Contents
--------

.. toctree::
   :maxdepth: 1

   About <about.rst>
   Installation instructions <installation.rst>
   Getting started <start.rst>
   user-guide.rst
   

.. _cite:

Citation and Attribution
------------------------
|Zenodo|

If you make use of this code, please cite the paper:

.. code-block::

   @article{peccaryPaper,
      author        = {{Hyman}, S{\'o}ley {\'O}. and {Daniel}, Kathryne J. and {Schaffner}, David A.},
      title         = "{PECCARY: A Novel Approach for Characterizing Orbital Complexity, Stochasticity, and Regularity}",
      journal       = {\apj},
      keywords      = {Theoretical techniques, Galaxy dynamics, Orbits, Orbit determination, Time series analysis, Astronomical methods, Astronomy software, 2093, 591, 1184, 1175, 1916, 1043, 1855, Instrumentation and Methods for Astrophysics, Astrophysics of Galaxies},
      year          = 2025,
      month         = jul,
      volume        = {987},
      number        = {2},
      eid           = {195},
      pages         = {195},
      doi           = {10.3847/1538-4357/adda3e},
      archivePrefix = {arXiv},
      eprint        = {2407.11970},
      primaryClass  = {astro-ph.IM},
      adsurl        = {https://ui.adsabs.harvard.edu/abs/2025ApJ...987..195H},
      adsnote       = {Provided by the SAO/NASA Astrophysics Data System}
   }

Please also cite the PECCARY version you used as a software citation using the Zenodo DOI |Zenodo|:

.. code-block::

   @software{peccaryZenodo,
      author       = {Hyman, Sóley and
                        Schaffner, David},
      title        = {soleyhyman/peccary: v0.1.1},
      month        = jul,
      year         = 2024,
      publisher    = {Zenodo},
      version      = {v0.1.1},
      doi          = {10.5281/zenodo.13168299},
      url          = {https://doi.org/10.5281/zenodo.13168299},
      swhid        = {swh:1:dir:2268b12e64f82631437df8024f5cd6904a85c396
                        ;origin=https://doi.org/10.5281/zenodo.13168298;vi
                        sit=swh:1:snp:d7900e21d7d057b5b231c93c616cdb6b8dcd
                        b9ae;anchor=swh:1:rel:73c52e8ebc846b6f6ba39fcad830
                        8c60fcc9481c;path=/
                        },
      }

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.13168299.svg
   :target: https://doi.org/10.5281/zenodo.13168299