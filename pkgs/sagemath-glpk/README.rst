=========================================================================================================
 Sage: Open Source Mathematics Software: Linear and mixed integer linear optimization backend using GLPK
=========================================================================================================

About SageMath
--------------

   "Creating a Viable Open Source Alternative to
    Magma, Maple, Mathematica, and MATLAB"

   Copyright (C) 2005-2022 The Sage Development Team

   https://www.sagemath.org

SageMath fully supports all major Linux distributions, recent versions of
macOS, and Windows (using Cygwin or Windows Subsystem for Linux).

The traditional and recommended way to install SageMath is from source via
Sage-the-distribution (https://www.sagemath.org/download-source.html).
Sage-the-distribution first builds a large number of open source packages from
source (unless it finds suitable versions installed in the system) and then
installs the Sage Library (sagelib, implemented in Python and Cython).


About this pip-installable source distribution
----------------------------------------------

This pip-installable source distribution ``sagemath-glpk`` provides
a backend for linear and mixed integer linear optimization backend using GLPK.

It can be installed as an extra of the distribution
`sagemath-polyhedra <https://pypi.org/project/sagemath-polyhedra>`_::

  $ pip install "sagemath-polyhedra[glpk]"


What is included
----------------

* `GLPK backends <https://doc.sagemath.org/html/en/reference/numerical/index.html#linear-optimization-lp-and-mixed-integer-linear-optimization-mip-solver-backends>`_ for LP, MILP, and graphs
