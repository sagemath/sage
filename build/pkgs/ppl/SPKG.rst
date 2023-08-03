ppl: Parma Polyhedra Library
============================

Description
-----------

The Parma Polyhedra Library (PPL) provides numerical abstractions
especially targeted at applications in the field of analysis and
verification of complex systems. These abstractions include convex
polyhedra, defined as the intersection of a finite number of (open or
closed) halfspaces, each described by a linear inequality (strict or
non-strict) with rational coefficients; some special classes of
polyhedra shapes that offer interesting complexity/precision tradeoffs;
and grids which represent regularly spaced points that satisfy a set of
linear congruence relations. The library also supports finite powersets
and products of (any kind of) polyhedra and grids, a mixed integer
linear programming problem solver using an exact-arithmetic version of
the simplex algorithm, a parametric integer programming solver, and
primitives for the termination analysis via the automatic synthesis of
linear ranking functions.

It is written in C++, but comes with interfaces to C, Java, OCaml, and
Prolog. PPL is one of the fastest implementations of polyhedral
computations.

Benchmarks are included in this paper: https://arxiv.org/abs/cs/0612085

License
-------

GPL v3+


Upstream Contact
----------------

- https://www.bugseng.com/ppl

Core Development Team

- Roberto Bagnara (University of Parma)
- Patricia M. Hill (University of Parma)
- Enea Zaffanella (University of Parma)
