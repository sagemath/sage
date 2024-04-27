"""
Generic matrices
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport sage.structure.element
cimport sage.structure.mutability

cdef class Matrix(sage.structure.element.Matrix):
    # Properties of any matrix  (plus _parent, inherited from base class)
    cdef public object _subdivisions
    cdef public object _base_ring
    cdef bint _is_immutable

    cpdef _add_(self, other)
    cpdef _sub_(self, other)

    cdef bint _will_use_strassen(self, Matrix right) except -2
    cdef bint _will_use_strassen_echelon(self) except -2
    cdef int _strassen_default_cutoff(self, Matrix right) except -2
    cdef int _strassen_default_echelon_cutoff(self) except -2

    # Implementation of hash function
    cdef long _hash_(self) except -1
    cdef void get_hash_constants(self, long C[5]) noexcept

    # Cache
    cdef public object _cache
    cdef long hash  # cached hash value
    cdef void clear_cache(self) noexcept
    cdef fetch(self, key)
    cdef cache(self, key, x)

    # Mutability and bounds checking
    cdef check_bounds(self, Py_ssize_t i, Py_ssize_t j)
    cdef check_mutability(self)
    cdef check_bounds_and_mutability(self, Py_ssize_t i, Py_ssize_t j)

    # Unsafe entry access
    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x)
    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j)
    cdef _coerce_element(self, x)
    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1

    # Row and column operations
    cdef check_row_bounds(self, Py_ssize_t r1, Py_ssize_t r2)
    cdef check_column_bounds(self, Py_ssize_t c1, Py_ssize_t c2)
    cdef check_row_bounds_and_mutability(self, Py_ssize_t r1, Py_ssize_t r2)
    cdef check_column_bounds_and_mutability(self, Py_ssize_t c1, Py_ssize_t c2)
    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2)
    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2)
    cdef add_multiple_of_row_c(self, Py_ssize_t i, Py_ssize_t j,    s, Py_ssize_t col_start)
    cdef add_multiple_of_column_c(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t row_start)
    cdef rescale_row_c(self, Py_ssize_t i, s, Py_ssize_t start_col)
    cdef rescale_col_c(self, Py_ssize_t i, s, Py_ssize_t start_row)

    # Helper function for inverse of sparse matrices
    cdef build_inverse_from_augmented_sparse(self, A)
