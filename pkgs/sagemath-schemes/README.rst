=============================================================================================================================================
 Sage: Open Source Mathematics Software: Schemes, varieties, elliptic curves, algebraic Riemann surfaces, modular forms, arithmetic dynamics
=============================================================================================================================================

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

This pip-installable source distribution `sagemath-schemes` is an experimental distribution of a part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`).


Status
------

The wheel builds. Some Cython modules that depend on FLINT or NTL are excluded.

`sage.all__sagemath_schemes` cannot be imported yet; tests disabled in `tox.ini`.
