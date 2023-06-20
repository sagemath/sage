====================================================================================================================
 Sage: Open Source Mathematics Software: Convex polyhedra in arbitrary dimension, mixed integer linear optimization
====================================================================================================================

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

This pip-installable source distribution `sagemath-polyhedra` is an experimental distribution of a part of the Sage Library.  Use at your own risk.  It provides a small subset of the modules of the Sage library ("sagelib", `sagemath-standard`), sufficient for computations with convex polyhedra in arbitrary dimension (in exact rational arithmetic), and linear and mixed integer linear optimization (in floating point arithmetic).


What is included
----------------

* `Combinatorial and Discrete Geometry <https://doc.sagemath.org/html/en/reference/discrete_geometry/index.html>`_: Polyhedra, lattice polyhedra, lattice points in polyhedra, triangulations, fans, polyhedral complexes, hyperplane arrrangements

* `Parma Polyhedra Library (PPL) backends for rational polyhedra <https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/backend_ppl.html>`_, `lattice polygons <https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/ppl_lattice_polygon.html>`_, `lattice polytopes <https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/ppl_lattice_polytope.html>`_; via `pplpy <https://doc.sagemath.org/html/en/reference/spkg/pplpy.html#spkg-pplpy>`_

* `Python backend for polyhedra over general ordered fields <https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/backend_field.html>`_

* `Linear, Mixed Integer Linear, and Semidefinite Optimization frontends <https://doc.sagemath.org/html/en/reference/numerical/index.html#numerical-optimization>`_

* `GNU Linear Programming Kit (GLPK) backend for large-scale linear and mixed integer linear optimization (floating point arithmetic) <https://doc.sagemath.org/html/en/reference/numerical/sage/numerical/backends/glpk_backend.html>`_

* `Interactive Simplex Method <https://doc.sagemath.org/html/en/reference/numerical/sage/numerical/interactive_simplex_method.html>`_


Available as extras, from other distributions
---------------------------------------------

Additional features:

`pip install "sagemath-polyhedra[graphs]"`
 Face lattices, combinatorial polyhedra, graph-theoretic constructions

`pip install "sagemath-polyhedra[groups]"`
 Constructing symmetric polyhedra, computing automorphisms, lattice point counting modulo group actions

`pip install "sagemath-polyhedra[toric]"`
 `Toric Varieties <https://doc.sagemath.org/html/en/reference/schemes/index.html#toric-varieties>`_

Other backends for polyhedral computations can be installed:

`pip install "sagemath-polyhedra[normaliz]"`
 `Normaliz <https://doc.sagemath.org/html/en/reference/spkg/normaliz.html#spkg-normaliz>`_, via `PyNormaliz <https://doc.sagemath.org/html/en/reference/spkg/pynormaliz.html#spkg-pynormaliz>`_

`pip install "sagemath-polyhedra[polymake]"`
 `Polymake <https://doc.sagemath.org/html/en/reference/spkg/polymake.html#spkg-polymake>`_, via `JuPyMake <https://pypi.org/project/JuPyMake/>`_

`sagemath-polyhedra` also provides integration with other packages for additional functionality:

* `LattE integrale <https://doc.sagemath.org/html/en/reference/spkg/latte_int.html#spkg-latte-int>`_
* `lrslib <https://doc.sagemath.org/html/en/reference/spkg/lrslib.html#spkg-lrslib>`_

Optional backends for optimization:

`pip install "sagemath-polyhedra[cbc]"`
 `COIN/OR CBC <https://doc.sagemath.org/html/en/reference/spkg/cbc.html#spkg-cbc>`_ Mixed Integer Linear Optimization solver,
 via `sage_numerical_backends_coin <https://doc.sagemath.org/html/en/reference/spkg/sage_numerical_backends_coin.html#spkg-sage-numerical-backends-coin>`_

`pip install "sagemath-polyhedra[cplex]"`
 CPLEX Mixed Integer Optimization solver (proprietary; requires licensed installation),
 via `sage_numerical_backends_cplex <https://doc.sagemath.org/html/en/reference/spkg/sage_numerical_backends_cplex.html#spkg-sage-numerical-backends-cplex>`_

`pip install "sagemath-polyhedra[cvxpy]"`
 `CVXPy <https://doc.sagemath.org/html/en/reference/spkg/cvxpy.html#spkg-cvxpy>`_ as middle-end for `various backends <https://www.cvxpy.org/install/>`_

`pip install "sagemath-polyhedra[gurobi]"`
 Gurobi Mixed Integer Optimization solver (proprietary; requires licensed installation), via `sage_numerical_backends_gurobi <https://doc.sagemath.org/html/en/reference/spkg/sage_numerical_backends_gurobi.html#spkg-sage-numerical-backends-gurobi>`_

`pip install "sagemath-polyhedra[scip]"`
 `SCIP <https://doc.sagemath.org/html/en/reference/spkg/scip.html#spkg-scip>`_ Mixed Integer Optimization solver,
 via `PySCIPOpt <https://doc.sagemath.org/html/en/reference/spkg/pyscipopt.html#spkg-pyscipopt>`_
