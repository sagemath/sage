# cython: wraparound=False, boundscheck=False
r"""
Cython helper methods to compute integral points in polyhedra

Note that while the URL of this documentation page ends with
``integral_points_generic_dense``, this is merely to allow Sphinx to generate
the documentation automatically. Imports should be from
:mod:`sage.geometry.integral_points`, as can be seen in the examples below.
Furthermore, not all functions are exported to the public interface.
"""

from sage.modules.vector_integer_dense cimport Vector_integer_dense as VectorClass
from sage.matrix.matrix_dense cimport Matrix_dense as MatrixClass

include "integral_points.pxi"
