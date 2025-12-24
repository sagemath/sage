highs: Linear optimization solver
=================================

Description
-----------

HiGHS is a high performance serial and parallel solver for large scale sparse
linear optimization problems of the form:

    min c'x subject to L <= Ax <= U; l <= x <= u

where c, x, L, U, l, u are vectors and A is a matrix.

HiGHS has implementations of the dual revised simplex method, primal and dual
revised simplex solvers, an interior point solver, and a MIP solver.

This package builds the HiGHS C/C++ library with C API headers for use in
Sage's numerical backends.

License
-------

- MIT License

Upstream Contact
----------------

- https://github.com/ERGO-Code/HiGHS
- https://www.highs.dev
- Email: hello@highs.dev

Dependencies
------------

- CMake (build dependency)

Special Update/Build Instructions
----------------------------------

The library is built using CMake and provides the C API interface
(highs_c_api.h) for integration with Sage's backend system.

