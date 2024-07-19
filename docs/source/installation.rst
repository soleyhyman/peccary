.. _installation:

Installation Instructions
=========================


Dependencies
------------

``peccary`` requires the use of 
`numpy <https://numpy.org/>`__ (version 1.14 or greater),
`scipy <https://scipy.org/>`__ (version 0.19 or greater), 
and `matplotlib <https://matplotlib.org/>`__.

If for some reason these packages are not automatically installed 
with the ``pip`` installation described :ref:`below <pipInstall>`,
you can install them via:
    
>>> pip install numpy>=1.14
>>> pip install scipy>=0.19
>>> pip install matplotlib

Installation
------------

.. _pipInstall:

Using ``pip`` (recommended)
:::::::::::::::::::::::::::

The recommended way to install the latest stable version of ``peccary`` 
is with ``pip`` via the terminal with the command:

>>> pip install peccary

You can also use the command:

>>> python -m pip install peccary


.. _gitInstall:

From source: 
::::::::::::
You can also clone the latest development version of ``peccary`` from 
`GitHub <https://github.com/>`_ using ``git``:

>>> git clone git://github.com/soleyhyman/peccary.git

To build and install the project (from the root of the source tree, e.g., inside
the cloned ``git`` directory):

>>> python -m pip install .