highs: Linear optimization solver with Python bindings
======================================================

Description
-----------

HiGHS is a high performance serial and parallel solver for large scale sparse
linear optimization problems of the form:

    min c'x subject to L <= Ax <= U; l <= x <= u

where c, x, L, U, l, u are vectors and A is a matrix.

HiGHS has implementations of the dual revised simplex method, primal and dual
revised simplex solvers, an interior point solver, and a MIP solver.

This package includes the highspy Python bindings built using pybind11.

License
-------

- MIT License

Upstream Contact
----------------

- https://github.com/ERGO-Code/HiGHS
- https://www.highs.dev
- Email: highsopt@gmail.com

Dependencies
------------

- CMake (build dependency)
- scikit-build-core (build dependency)
- pybind11 (build dependency)
- numpy

Special Update/Build Instructions
----------------------------------

The Python bindings (highspy) are built automatically via pip using the
pyproject.toml included in the HiGHS source distribution.
