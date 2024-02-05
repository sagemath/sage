Matrices and Spaces of Matrices
===============================

Sage provides native support for working with matrices over any
commutative or noncommutative ring. The parent object for a matrix
is a matrix space ``MatrixSpace(R, n, m)`` of all
:math:`n\times
m` matrices over a ring :math:`R`.

To create a matrix, either use the ``matrix(...)``
function or create a matrix space using the
``MatrixSpace`` command and coerce an object into it.

Matrices also act on row vectors, which you create using the
``vector(...)`` command or by making a
``VectorSpace`` and coercing lists into it. The natural
action of matrices on row vectors is from the right. Sage currently
does not have a column vector class (on which matrices would act
from the left), but this is planned.

In addition to native Sage matrices, Sage also includes the
following additional ways to compute with matrices:


-  Several math software systems included with Sage have their own
   native matrix support, which can be used from Sage. E.g., PARI,
   GAP, Maxima, and Singular all have a notion of matrices.

-  The GSL C-library is included with Sage, and can be used via
   Cython.

-  The :mod:`scipy:scipy` module provides support for
   *sparse* numerical linear algebra, among many other things.

-  The ``numpy`` module, which you load by typing
   ``import numpy`` is included standard with Sage. It
   contains a very sophisticated and well developed array class, plus
   optimized support for *numerical linear algebra*.  Sage's matrices
   over RDF and CDF (native floating-point real and complex numbers)
   use numpy.

Finally, this module contains some data-structures for matrix-like
objects like operation tables (e.g. the multiplication table of a group).

.. toctree::
   :maxdepth: 1


   sage/matrix/matrix_space
   sage/matrix/constructor
   sage/matrix/special
   sage/matrix/args

   sage/matrix/docs
   sage/matrix/matrix0
   sage/matrix/matrix1
   sage/matrix/matrix2

   sage/matrix/strassen
   sage/matrix/berlekamp_massey

   sage/matrix/matrix_dense
   sage/matrix/matrix_sparse

   sage/matrix/matrix_generic_dense
   sage/matrix/matrix_generic_sparse

   sage/matrix/matrix_integer_dense
   sage/matrix/matrix_integer_sparse
   sage/matrix/matrix_integer_dense_hnf
   sage/matrix/matrix_integer_dense_saturation
   sage/matrix/matrix_rational_dense
   sage/matrix/matrix_rational_sparse
   sage/matrix/matrix_double_dense
   sage/matrix/matrix_real_double_dense
   sage/matrix/matrix_mod2_dense
   sage/matrix/matrix_gf2e_dense
   sage/matrix/matrix_modn_dense_double
   sage/matrix/matrix_modn_dense_float
   sage/matrix/matrix_modn_sparse
   sage/matrix/matrix_symbolic_dense
   sage/matrix/matrix_symbolic_sparse
   sage/matrix/matrix_complex_double_dense
   sage/matrix/matrix_complex_ball_dense
   sage/matrix/matrix_polynomial_dense
   sage/matrix/matrix_mpolynomial_dense
   sage/matrix/matrix_cyclo_dense
   sage/matrix/matrix_cmr_sparse

   sage/matrix/operation_table

   sage/matrix/action
   sage/matrix/change_ring
   sage/matrix/echelon_matrix
   sage/matrix/matrix_misc
   sage/matrix/matrix_window
   sage/matrix/misc
   sage/matrix/misc_mpfr
   sage/matrix/misc_flint
   sage/matrix/symplectic_basis
   sage/matrix/compute_J_ideal
   sage/matrix/seymour_decomposition

   sage/matrix/benchmark

.. include:: ../footer.txt
