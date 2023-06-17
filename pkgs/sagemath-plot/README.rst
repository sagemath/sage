===============================================================================================
 Sage: Open Source Mathematics Software: Plotting and graphics with Matplotlib, Three.JS, etc.
===============================================================================================

About SageMath
--------------

   "Creating a Viable Open Source Alternative to
    Magma, Maple, Mathematica, and MATLAB"

   Copyright (C) 2005-2023 The Sage Development Team

   https://www.sagemath.org

SageMath fully supports all major Linux distributions, recent versions of macOS, and Windows (using Cygwin or Windows Subsystem for Linux).

The traditional and recommended way to install SageMath is from source via Sage-the-distribution (https://www.sagemath.org/download-source.html).  Sage-the-distribution first builds a large number of open source packages from source (unless it finds suitable versions installed in the system) and then installs the Sage Library (sagelib, implemented in Python and Cython).


About this experimental pip-installable source distribution
-----------------------------------------------------------

This pip-installable source distribution `sagemath-plot` is an experimental distribution of a part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`).


What is included
----------------

* `2D Graphics <https://doc.sagemath.org/html/en/reference/plotting/index.html>`_

* `3D Graphics <https://doc.sagemath.org/html/en/reference/plot3d/index.html>`_

* Interfaces: `Gnuplot <https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/gnuplot.html>`_, `Jmol <https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/jmoldata.html>`_, `POV-Ray <https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/povray.html>`_, `Tachyon <https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/tachyon.html>`_
