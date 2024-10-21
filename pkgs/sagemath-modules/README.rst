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

This pip-installable source distribution `sagemath-modules` is an experimental distribution of a part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`).


What is included
----------------

* `Vectors, Vector Spaces, Modules <https://doc.sagemath.org/html/en/reference/modules/index.html>`_

* `Matrices and Spaces of Matrices <https://doc.sagemath.org/html/en/reference/matrices/index.html>`_

* Fields of real and complex numbers in arbitrary precision floating point arithmetic (using MPFR, GSL, mpmath, MPC)

* `Free Modules with Combinatorial Bases <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/free_module.html>`_

* `Tensor Modules <https://doc.sagemath.org/html/en/reference/tensor_free_modules/index.html>`_

* `Additive Abelian Groups <https://doc.sagemath.org/html/en/reference/groups/sage/groups/additive_abelian/additive_abelian_group.html>`_

* `Matrix and Affine Groups <https://doc.sagemath.org/html/en/reference/groups/index.html#matrix-and-affine-groups>`_

* `Root Systems <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/all.html#sage-combinat-root-system-all>`_

* `Quadratic Forms <https://doc.sagemath.org/html/en/reference/quadratic_forms/index.html>`_

* `Ring Extensions <https://doc.sagemath.org/html/en/reference/rings/sage/rings/ring_extension.html>`_ and `Derivations <https://doc.sagemath.org/html/en/reference/rings/sage/rings/derivation.html>`_

* `Clifford, Exterior <https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/clifford_algebra.html>`_, and  `Weyl Algebras <https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/weyl_algebra.html>`_

* `Chain Complexes, Homology <https://doc.sagemath.org/html/en/reference/homology/index.html>`_, `Free Resolutions <https://doc.sagemath.org/html/en/reference/resolutions/index.html>`_

* `Matroid Theory <https://doc.sagemath.org/html/en/reference/matroids/index.html>`_

* `Coding Theory <https://doc.sagemath.org/html/en/reference/coding/index.html>`_

* `Cryptography <https://doc.sagemath.org/html/en/reference/cryptography/index.html>`_

* `Probability Spaces and Distributions <https://doc.sagemath.org/html/en/reference/probability/index.html>`_, `Statistics <https://doc.sagemath.org/html/en/reference/stats/index.html>`_


Available as extras, from other distributions
---------------------------------------------

`pip install "sagemath-modules[RDF,CDF]"`
 Linear algebra over fields of real and complex numbers using NumPy

`pip install "sagemath-modules[RBF,CBF]"`
 Linear algebra over fields of real and complex numbers with ball arithmetic using FLINT/arb

`pip install "sagemath-modules[GF,GF2,GF2e,GFpn]"`
 Linear algebra over finite fields (various implementations)

`pip install "sagemath-modules[QQbar,NumberField,CyclotomicField]"`
 Linear algebra over the algebraic numbers or number fields

`pip install "sagemath-modules[padics]"`
 Linear algebra over p-adic rings and fields

`pip install "sagemath-modules[combinat]"`
 Modules and algebras with combinatorial bases; algebraic combinatorics

`pip install "sagemath-modules[invariant]"`
 Submodules invariant under group actions

`pip install "sagemath-modules[standard]"`
 All related features as in a standard installation of SageMath
