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


What is included
----------------

* `Ideals and Varieties <https://doc.sagemath.org/html/en/reference/polynomial_rings/sage/rings/polynomial/multi_polynomial_ideal.html>`_

* `Schemes <https://doc.sagemath.org/html/en/reference/schemes/index.html>`_

* `Plane and Space Curves <https://doc.sagemath.org/html/en/reference/curves/index.html>`_

* `Elliptic and Hyperelliptic Curves <https://doc.sagemath.org/html/en/reference/arithmetic_curves/index.html>`_

* `Modular Forms <https://doc.sagemath.org/html/en/reference/modfrm/index.html>`_

* `Modular Symbols <https://doc.sagemath.org/html/en/reference/modsym/index.html>`_

* `Modular Abelian Varieties <https://doc.sagemath.org/html/en/reference/modabvar/index.html>`_

* `Arithmetic Dynamical Systems <https://doc.sagemath.org/html/en/reference/dynamics/index.html#arithmetic-dynamical-systems>`_


Status
------

The wheel builds. Some Cython modules that depend on FLINT or NTL are excluded.

`sage.all__sagemath_schemes` can be imported.

Many tests fail; see ``known-test-failures.json``.
