========================================================================================
 Sage: Open Source Mathematics Software: Matrices over small finite fields with meataxe
========================================================================================

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

This pip-installable source distribution ``sagemath-meataxe`` is a small
optional distribution for use with ``sagemath-standard``.

This distribution provides the SageMath modules ``sage.libs.meataxe``
and ``sage.matrix.matrix_gfpn_dense``.

It provides a specialized implementation of matrices over the finite field F_q, where
q <= 255, using the `SharedMeatAxe <http://users.minet.uni-jena.de/~king/SharedMeatAxe/>`
library.
