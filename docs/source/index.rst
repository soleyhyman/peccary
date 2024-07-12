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
is a pure-python package for distinguishing between regular, complex, and stochastic
behavior in timeseries. It is based on the work by 
`Bandt & Pompe <https://ui.adsabs.harvard.edu/#abs/2002PhRvL..88q4102B/abstract>`_ (2002), 
`Rosso et al. <https://ui.adsabs.harvard.edu/#abs/2007PhRvL..99o4102R/abstract>`_ (2007). 
and 
`Weck et al. <https://ui.adsabs.harvard.edu/#abs/2015PhRvE..91b3101W/abstract>`_ (2015),
In addition to calculating the Permutation Entropy (:math:`H`) and Statistical Complexity
(:math:`C`) values, this package also has plotting tools for the :math:`HC`-plane and visualizing the 
resulting :math:`[H,C]` values for various timeseries.

A detailed summary of the PECCARY method can be found in Hyman, Daniel, & Schaffner (in prep). 
If you make use of PECCARY, please include a citation to Hyman, Daniel, & Schaffner (in prep) 
in any publications.

.. note::
   This project is under active development.


Contents
--------

.. toctree::
   :maxdepth: 1

   Introduction <intro.rst>
   Installation instructions <installation.rst>
   Getting started <start.rst>
   user-guide.rst
   

Citation and Attribution
------------------------

If you make use of this code, please cite the paper:

.. .. code-block:: bibtex

..     @article{gala,
..       doi = {10.21105/joss.00388},
..       url = {https://doi.org/10.21105%2Fjoss.00388},
..       year = 2017,
..       month = {oct},
..       publisher = {The Open Journal},
..       volume = {2},
..       number = {18},
..       author = {Adrian M. Price-Whelan},
..       title = {Gala: A Python package for galactic dynamics},
..       journal = {The Journal of Open Source Software}}

.. Please also cite the Zenodo DOI |DOI| of the version you used as a software
.. citation:

.. .. include:: ZENODO.rst

.. .. |JOSS| image:: http://joss.theoj.org/papers/10.21105/joss.00388/status.svg
..    :target: http://joss.theoj.org/papers/10.21105/joss.00388
.. .. |DOI| image:: https://zenodo.org/badge/17577779.svg
..    :target: https://zenodo.org/badge/latestdoi/17577779


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`