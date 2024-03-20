===========================================================
 Sage: Open Source Mathematics Software: Symbolic calculus
===========================================================

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

This pip-installable source distribution `sagemath-symbolics` is an experimental distribution of a part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`).


What is included
----------------

* `Symbolic Calculus <https://doc.sagemath.org/html/en/reference/calculus/index.html>`_

* `Pynac <http://pynac.org/>`_ (fork of GiNaC)

* `Maxima <https://doc.sagemath.org/html/en/reference/spkg/maxima.html>`_, running on `ECL <https://doc.sagemath.org/html/en/reference/spkg/ecl.html>`_

* Arithmetic Functions, `Elementary and Special Functions <https://doc.sagemath.org/html/en/reference/functions/index.html>`_
  (via `sagemath-categories <https://doc.sagemath.org/html/en/reference/spkg/sagemath_categories.html>`_)

* `Asymptotic Expansions <https://doc.sagemath.org/html/en/reference/asymptotic/index.html>`_

* `Topological, Differentiable, Pseudo-Riemannian, Poisson Manifolds <https://doc.sagemath.org/html/en/reference/manifolds/index.html>`_

* `Hyperbolic Geometry <https://doc.sagemath.org/html/en/reference/hyperbolic_geometry/index.html>`_



Available as extras, from other distributions
---------------------------------------------

`pip install "sagemath-symbolics[giac]"`
 Computer algebra system `GIAC <https://doc.sagemath.org/html/en/reference/spkg/giac.html>`_, via `sagemath-giac <https://doc.sagemath.org/html/en/reference/spkg/sagemath_giac.html>`_

`pip install "sagemath-symbolics[primecount]"`
 `Prime counting function <https://doc.sagemath.org/html/en/reference/functions/sage/functions/prime_pi.html>`_
 implementation `primecount <https://doc.sagemath.org/html/en/reference/spkg/primecount.html>`_, via `primecountpy <https://doc.sagemath.org/html/en/reference/spkg/primecountpy.html>`_

`pip install "sagemath-symbolics[sympy]"`
 Python library for symbolic mathematics / computer algebra system `SymPy <https://doc.sagemath.org/html/en/reference/spkg/sympy.html>`_
