r"""
Inline conversions between LinBox and Sage

Each LinBox type has a corresponding Sage types and we use the following
conventions for conversion functions

- ``new_linbox_XXX`` : create a new linbox object
- ``new_sage_XXX``   : create a new Sage object
- ``set_linbox_XXX`` : set the entries of the linbox object
- ``set_sage_XXX``   : set the entries of the Sage object

For matrices that uses a flint datastructure, see the lower level conversions
in the module ``linbox_flint_interface``.
"""
#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht
#       Copyright (C) 2008 Clement Pernet
#       Copyright (C) 2018 Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from libcpp.vector cimport vector as cppvector

from sage.libs.gmp.mpz cimport mpz_set

from sage.libs.linbox.givaro cimport Modular_uint64, ZRing, Integer
from sage.libs.linbox.linbox cimport SparseMatrix_Modular_uint64, SparseMatrix_integer, DenseVector_integer

from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse
from sage.matrix.matrix_integer_sparse cimport Matrix_integer_sparse

from sage.modules.vector_modn_sparse cimport c_vector_modint
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.modules.vector_integer_sparse cimport mpz_vector,  mpz_vector_get_entry, mpz_vector_set_entry

########################################
# algorithm for solving linear systems #
########################################

ctypedef enum linbox_specifier:
    METHOD_DEFAULT              # no specification
    METHOD_DENSE_ELIMINATION     # DenseElimination
    METHOD_SPARSE_ELIMINATION   # SparseElimination
    METHOD_BLACKBOX             # Blackbox
    METHOD_WIEDEMANN            # Wiedeman
    ERROR

cdef inline linbox_specifier get_method(str algo) except ERROR:
    if algo is None or algo == "default":
        return METHOD_DEFAULT
    elif algo == "dense_elimination" or \
         algo == "linbox_dense_elimination" or \
         algo == "LinBox::DenseElimination":
        return METHOD_DENSE_ELIMINATION
    elif algo == "sparse_elimination" or \
         algo == "linbox_sparse_elimination" or \
         algo == "LinBox::SparseElimination":
        return METHOD_SPARSE_ELIMINATION
    elif algo == "blackbox" or \
         algo == "linbox_blackbox" or \
         algo == "LinBox::Blackbox":
        return METHOD_BLACKBOX
    elif algo == 'wiedemann' or \
         algo == "linbox_wiedemann" or \
         algo == "LinBox::Wiedeman":
        return METHOD_WIEDEMANN
    else:
        raise ValueError("unknown algorithm")

################################################
# matrix_modn_sparse (sparse matrix over Z/nZ) #
################################################

cdef inline void set_linbox_matrix_modn_sparse(SparseMatrix_Modular_uint64& A, Matrix_modn_sparse m) noexcept:
    r"""
    Set the entries of a LinBox matrix from a Sage matrix.

    INPUT:

    - ``A`` -- LinBox matrix
    - ``m`` -- Sage matrix
    """
    cdef c_vector_modint * row
    cdef size_t i, j
    for i in range(<size_t> m._nrows):
        row = m.rows + i
        for j in range(<size_t> row.num_nonzero):
            A.setEntry(i, row.positions[j], row.entries[j])

cdef inline SparseMatrix_Modular_uint64 * new_linbox_matrix_modn_sparse(Modular_uint64 &F, Matrix_modn_sparse m) noexcept:
    r"""
    Return a new LinBox matrix from a Sage matrix.

    Such matrix has to be deallocated with a "del" statement.

    INPUT:

    - ``F`` -- LinBox field
    - ``m`` -- Sage matrix
    """
    cdef SparseMatrix_Modular_uint64 * A = new SparseMatrix_Modular_uint64(F, <size_t> m._nrows, <size_t> m._ncols)
    set_linbox_matrix_modn_sparse(A[0], m)
    return A

#########################
# matrix_integer_sparse #
#########################

cdef inline void set_linbox_matrix_integer_sparse(SparseMatrix_integer& A, Matrix_integer_sparse m) noexcept:
    r"""
    Set the entries of a LinBox matrix from a Sage matrix.

    INPUT:

    - ``A`` -- LinBox matrix
    - ``m`` -- Sage matrix
    """
    cdef size_t i, j, k
    cdef mpz_vector * v
    cdef Integer t
    for i in range(<size_t> m._nrows):
        v = m._matrix + i
        for k in range(<size_t> v.num_nonzero):
            j = v.positions[k]
            mpz_set(t.get_mpz(), v.entries[k])
            A.setEntry(i, j, t)

cdef inline SparseMatrix_integer * new_linbox_matrix_integer_sparse(ZRing &ZZ, Matrix_integer_sparse m) noexcept:
    r"""
    Return a new LinBox matrix from a Sage matrix.

    Suc matrix has to be deallocated with a "del" statement.

    INPUT:

    - ``m`` -- Sage matrix
    """
    cdef SparseMatrix_integer * A = new SparseMatrix_integer(ZZ, <size_t> m._nrows, <size_t> m._ncols)
    set_linbox_matrix_integer_sparse(A[0], m)
    return A

########################
# vector integer dense #
########################

cdef inline DenseVector_integer * new_linbox_vector_integer_dense(ZRing &ZZ, Vector_integer_dense v) noexcept:
    r"""
    Return a new linbox vector from a sage one.

    INPUT:

    - ``v`` -- a Sage dense integer vector
    """
    cdef cppvector[Integer] * vec = new cppvector[Integer](<size_t> v._degree)
    cdef size_t i
    for i in range(<size_t> v._degree):
        mpz_set(vec[0][i].get_mpz(), v._entries[i])

    cdef DenseVector_integer * V = new DenseVector_integer(ZZ, vec[0])
    del vec

    return V

cdef inline Vector_integer_dense new_sage_vector_integer_dense(P, DenseVector_integer &v):
    r"""
    Return a new Sage vector from a LinBox one.

    INPUT:

    - ``P`` -- parent for the Sage vector
    - ``v`` -- linbox vector
    """
    cdef Vector_integer_dense res = P()
    cdef size_t i
    for i in range(<size_t> res._degree):
        mpz_set(res._entries[i], v.getEntry(i).get_mpz_const())

    return res
