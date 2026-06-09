r"""
Utility functions for matrices

Inline utility functions for matrices. This module is mainly to reduce code
duplication, as well as ensure error messages are kept consistent.

AUTHORS:

- Vincent Macri (2026): refactored old code into `check_matrix_multiplication_sizes`
"""

# ****************************************************************************
#       Copyright (C) 2026 Vincent Macri <vincent.macri@ucalgary.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element cimport Matrix

cdef inline void check_matrix_multiplication_sizes(Matrix left, Matrix right) except *:
    if left._ncols != right._nrows:
        raise ArithmeticError("number of columns of left must equal number of rows of right")
