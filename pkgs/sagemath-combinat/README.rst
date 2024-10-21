======================================================================================================
 Sage: Open Source Mathematics Software: Algebraic combinatorics, combinatorial representation theory
======================================================================================================

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

This pip-installable source distribution `sagemath-combinat` is an experimental distribution of a part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`).


What is included
----------------

* `Enumerative Combinatorics <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/enumerated_sets.html#sage-combinat-enumerated-sets>`_: `Partitions, Tableaux <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/catalog_partitions.html>`_

* `Combinatorics on Words <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/words/all.html#sage-combinat-words-all>`_, `Free Monoids <https://doc.sagemath.org/html/en/reference/monoids/index.html>`_, `Automatic Semigroups <https://doc.sagemath.org/html/en/reference/monoids/sage/monoids/automatic_semigroup.html>`_

* `Symmetric Functions <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/sf/all.html#sage-combinat-sf-all>`_, other `Algebras with combinatorial bases <https://doc.sagemath.org/html/en/reference/algebras/index.html>`_


Available as extras, from other distribution packages
-----------------------------------------------------

* `sagemath-graphs <https://pypi.org/project/sagemath-graphs>`_:
  Graphs, posets, finite state machines, combinatorial designs, incidence structures, quivers

* `sagemath-modules <https://pypi.org/project/sagemath-modules>`_:
  Modules and algebras, root systems, coding theory

* `sagemath-polyhedra <https://pypi.org/project/sagemath-polyhedra>`_:
  Polyhedra, lattice points, hyperplane arrangements
