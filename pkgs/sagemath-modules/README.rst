===========================================================================================================================================================================================================
 Sage: Open Source Mathematics Software: Vectors, matrices, tensors, vector spaces, affine spaces, modules and algebras, additive groups, quadratic forms, root systems, homology, coding theory, matroids
===========================================================================================================================================================================================================

About SageMath
--------------

   "Creating a Viable Open Source Alternative to
    Magma, Maple, Mathematica, and MATLAB"

   Copyright (C) 2005-2020 The Sage Development Team

   https://www.sagemath.org

SageMath fully supports all major Linux distributions, recent versions of macOS, and Windows (using Cygwin or Windows Subsystem for Linux).

The traditional and recommended way to install SageMath is from source via Sage-the-distribution (https://www.sagemath.org/download-source.html).  Sage-the-distribution first builds a large number of open source packages from source (unless it finds suitable versions installed in the system) and then installs the Sage Library (sagelib, implemented in Python and Cython).


About this experimental pip-installable source distribution
-----------------------------------------------------------

This pip-installable source distribution `sagemath-polyhedra` is an experimental distribution of a part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`), sufficient for computations with convex polyhedra in arbitrary dimension, in exact rational arithmetic.

By default, `sagemath-polyhedra` makes use of the Parma Polyhedra Library (PPL) as a backend for the computations, via https://pypi.org/project/pplpy/, which is a dependency of this distribution.

Using `sagemath-polyhedra`'s unified interface, other optional backends can also be used:

* Normaliz, via https://pypi.org/project/PyNormaliz/,
* Polymake, via https://pypi.org/project/JuPyMake/;
* `cddlib <https://github.com/cddlib/cddlib>`_

`sagemath-polyhedra` also provides integration with other packages for additional functionality:

* `LattE integrale <https://www.math.ucdavis.edu/~latte/software.php>`_
* `lrslib <http://cgm.cs.mcgill.ca/~avis/C/lrs.html>`_

Documentation
-------------

* `Combinatorial and Discrete Geometry <https://doc.sagemath.org/html/en/reference/discrete_geometry/index.html>`_
