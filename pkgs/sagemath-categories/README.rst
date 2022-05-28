=========================================================================
 Sage: Open Source Mathematics Software: Sage categories and basic rings
=========================================================================

About SageMath
--------------

   "Creating a Viable Open Source Alternative to
    Magma, Maple, Mathematica, and MATLAB"

   Copyright (C) 2005-2022 The Sage Development Team

   https://www.sagemath.org

SageMath fully supports all major Linux distributions, recent versions of macOS, and Windows (using Cygwin or Windows Subsystem for Linux).

The traditional and recommended way to install SageMath is from source via Sage-the-distribution (https://www.sagemath.org/download-source.html).  Sage-the-distribution first builds a large number of open source packages from source (unless it finds suitable versions installed in the system) and then installs the Sage Library (sagelib, implemented in Python and Cython).


About this experimental pip-installable source distribution
-----------------------------------------------------------

This pip-installable source distribution `sagemath-categories` is an experimental distribution of a small part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`).  It is a superset of the `sagemath-objects` (providing Sage objects, the element/parent framework, categories, the coercion system and the related metaclasses), making various additional categories available without introducing dependencies on additional mathematical libraries.


Dependencies
------------

When building from source, development packages of `gmp`, `mpfr`, and `mpc` are needed.


Documentation
-------------

* `Categories <https://doc.sagemath.org/html/en/reference/categories/index.html>`_

* `Structure <https://doc.sagemath.org/html/en/reference/structure/index.html>`_

* `Coercion <https://doc.sagemath.org/html/en/reference/coercion/index.html>`_

* `Classes, Metaclasses <https://doc.sagemath.org/html/en/reference/misc/index.html#special-base-classes-decorators-etc>`_
