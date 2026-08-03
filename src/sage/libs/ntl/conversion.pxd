r"""
Inline conversions between NTL and Sage

Each NTL type has a corresponding Sage types and we use the following
conventions for conversion functions

- ``new_ntl_XXX`` : create a new ntl object
- ``new_sage_XXX``   : create a new Sage object
- ``set_ntl_XXX`` : set the entries of the ntl object
- ``set_sage_XXX``   : set the entries of the Sage object
"""
# ***************************************************************************
#       Copyright (C) 2007 Martin Albrecht
#                     2008 Clement Pernet
#                     2018 Vincent Delecroix
#                     2018 Alex J. Best
#                     2026 Vincent Macri <vincent.macri@ucalgary.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ***************************************************************************

from sage.libs.ntl.types cimport mat_ZZ_p_c
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p

from sage.matrix.matrix_modn_dense_float cimport Matrix_modn_dense_float
from sage.matrix.matrix_modn_dense_double cimport Matrix_modn_dense_double
from sage.matrix.matrix_modn_dense_flint cimport Matrix_modn_dense_flint
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense


ctypedef fused Matrix_modn_dense_sage:
    Matrix_modn_dense_float
    Matrix_modn_dense_double
    Matrix_modn_dense_flint
    Matrix_generic_dense

cdef inline void set_ntl_matrix_modn_dense_sage(mat_ZZ_p_c& A, ntl_ZZ_pContext_class c, Matrix_modn_dense_sage m) noexcept:
    r"""
    Set the entries of a NTL matrix from a Sage matrix.

    INPUT:

    - ``A`` -- NTL matrix
    - ``m`` -- Sage matrix
    """
    cdef size_t i, j
    cdef ntl_ZZ_p tmp
    A.SetDims(m._nrows, m._ncols)
    for i in range(m._nrows):
        for j in range(m._ncols):
            tmp = ntl_ZZ_p(m[i, j], c)
            A.put(i, j, tmp.x)

cdef inline void set_ntl_matrix_modn_dense(mat_ZZ_p_c& A, ntl_ZZ_pContext_class c, m) noexcept:
    r"""
    set the entries of a NTL matrix from a Sage matrix.

    INPUT:

    - ``A`` -- NTL matrix
    - ``m`` -- Sage matrix
    """

    if isinstance(m, Matrix_modn_dense_float):
        set_ntl_matrix_modn_dense_sage(A, c, <Matrix_modn_dense_float> m)
    elif isinstance(m, Matrix_modn_dense_double):
        set_ntl_matrix_modn_dense_sage(A, c, <Matrix_modn_dense_double> m)
    elif isinstance(m, Matrix_modn_dense_flint):
        set_ntl_matrix_modn_dense_sage(A, c, <Matrix_modn_dense_flint> m)
    elif isinstance(m, Matrix_generic_dense):
        set_ntl_matrix_modn_dense_sage(A, c, <Matrix_generic_dense> m)
    else:
        raise NotImplementedError("matrix type not yet implemented")
