==========================================================================
 Sage: Open Source Mathematics Software: Combinatorial matrix recognition
==========================================================================

About SageMath
--------------

   "Creating a Viable Open Source Alternative to
    Magma, Maple, Mathematica, and MATLAB"

   Copyright (C) 2005-2023 The Sage Development Team

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

This pip-installable source distribution ``sagemath-cmr`` is a small
optional distribution for use with ``sagemath-standard``.

It provides a Cython interface to the CMR library (https://github.com/discopt/cmr),
providing recognition and decomposition algorithms for:

- Totally Unimodular Matrices
- Network Matrices
- Complement Totally Unimodular Matrices
- (Strongly) k-Modular and Unimodular Matrices
- Regular Matroids
- Graphic / Cographic / Planar Matrices
- Series-Parallel Matroids
