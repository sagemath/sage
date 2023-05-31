==============================================================================================
 Sage: Open Source Mathematics Software: Sage categories, basic rings, polynomials, functions
==============================================================================================

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


About this experimental pip-installable source distribution
-----------------------------------------------------------

This pip-installable source distribution `sagemath-categories` is an
experimental distribution of a small part of the Sage Library.  Use at your own
risk.  It provides a small subset of the modules of the Sage library
("sagelib", `sagemath-standard`).  It is a superset of the `sagemath-objects`
(providing Sage objects, the element/parent framework, categories, the coercion
system and the related metaclasses), making various additional categories
available without introducing dependencies on additional mathematical
libraries.


What is included
----------------

* `Structure <https://doc.sagemath.org/html/en/reference/structure/index.html>`_, `Coercion framework <https://doc.sagemath.org/html/en/reference/coercion/index.html>`_, `Base Classes, Metaclasses <https://doc.sagemath.org/html/en/reference/misc/index.html#special-base-classes-decorators-etc>`_

* `Categories and functorial constructions <https://doc.sagemath.org/html/en/reference/categories/index.html>`_

* `Sets <https://doc.sagemath.org/html/en/reference/sets/index.html>`_ (except `RealSet`)

* Basic Combinatorial and Data Structures: `Binary trees <https://doc.sagemath.org/html/en/reference/data_structures/sage/misc/binary_tree.html>`_, `Bitsets <https://doc.sagemath.org/html/en/reference/data_structures/sage/data_structures/bitset.html>`_, `Permutations <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/permutation.html>_`,

* Basic Rings and Fields: `Integers, Rationals <https://doc.sagemath.org/html/en/reference/rings_standard/index.html>`_

* `Commutative Polynomials <https://doc.sagemath.org/html/en/reference/polynomial_rings/index.html>`_, `Power Series and Laurent Series <https://doc.sagemath.org/html/en/reference/power_series/index.html>`_,

* Arithmetic Functions, `Elementary and Special Functions <https://doc.sagemath.org/html/en/reference/functions/index.html>`_ as generic entry points


Available in other distribution packages
----------------------------------------

* `sagemath-combinat <https://pypi.org/project/sagemath-combinat>`_:
  Algebraic combinatorics, combinatorial representation theory

* `sagemath-graphs <https://pypi.org/project/sagemath-graphs>`_:
  Graphs, posets, hypergraphs, designs, abstract complexes, combinatorial polyhedra, abelian sandpiles

* `sagemath-groups <https://pypi.org/project/sagemath-groups>`_:
  Groups, semigroups, invariant theory

* `sagemath-modules <https://pypi.org/project/sagemath-modules>`_:
  Vectors, matrices, tensors, vector spaces, affine spaces,
  modules and algebras, additive groups, quadratic forms, homology, coding theory, matroids

* `sagemath-plot <https://pypi.org/project/sagemath-plot>`_:
  Plotting and graphics with Matplotlib, Three.JS, etc.

* `sagemath-polyhedra <https://pypi.org/project/sagemath-polyhedra>`_:
  Convex polyhedra in arbitrary dimension, triangulations, polyhedral fans, lattice points, geometric complexes, hyperplane arrangements

* `sagemath-repl <https://pypi.org/project/sagemath-repl>`_:
  IPython REPL, the interactive language of SageMath (preparser), interacts, development tools

* `sagemath-symbolics <https://pypi.org/project/sagemath-symbolics>`_:
  Symbolic expressions, calculus, differentiable manifolds, schemes, varieties, Groebner bases, asymptotics



Dependencies
------------

When building from source, development packages of `gmp`, `mpfr`, and `mpc` are needed.
