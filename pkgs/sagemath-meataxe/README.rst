========================================================================================
 Sage: Open Source Mathematics Software: Matrices over small finite fields with meataxe
========================================================================================

About SageMath
--------------

   "Creating a Viable Open Source Alternative to
    Magma, Maple, Mathematica, and MATLAB"

   Copyright (C) 2005-2024 The Sage Development Team

   https://www.sagemath.org

SageMath fully supports all major Linux distributions, recent versions of
macOS, and Windows (Windows Subsystem for Linux).

See https://doc.sagemath.org/html/en/installation/index.html
for general installation instructions.


About this pip-installable source distribution
----------------------------------------------

This pip-installable source distribution ``sagemath-meataxe`` is a small
optional distribution for use with ``sagemath-standard``.

This distribution provides the SageMath modules ``sage.libs.meataxe``
and ``sage.matrix.matrix_gfpn_dense``.

It provides a specialized implementation of matrices over the finite field F_q, where
q <= 255, using the `SharedMeatAxe <http://users.minet.uni-jena.de/~king/SharedMeatAxe/>`
library.
