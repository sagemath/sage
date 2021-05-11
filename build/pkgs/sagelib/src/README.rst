=================================================================
 Sage: Open Source Mathematics Software: Standard Python Library
=================================================================

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

This pip-installable source distribution `sagemath-standard` is an experimental distribution of the Sage Library.  Use at your own risk.

Building `sagemath-standard` has a large number of system packages as prerequisites. See https://doc.sagemath.org/html/en/installation/source.html#linux-recommended-installation
for partial lists for various systems.

A modularization effort is in progress with the goal of making it possible to install parts of the Sage Library with fewer prerequisites. https://trac.sagemath.org/ticket/29705
