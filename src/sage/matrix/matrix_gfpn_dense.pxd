# sage_setup: distribution = sagemath-meataxe

# ***************************************************************************
#       Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.matrix.matrix_dense cimport Matrix_dense
from sage.structure.element cimport Matrix
from sage.libs.meataxe cimport *


cdef class FieldConverter_class:
    cdef field  # A function converting an int to a field element
    cdef FEL zero_FEL  # the FEL representation of zero
    cpdef fel_to_field(self, FEL x) noexcept
    cpdef FEL field_to_fel(self, x) except 255

cdef FieldConverter_class FieldConverter(field) noexcept

cdef class Matrix_gfpn_dense(Matrix_dense):
    cdef Matrix_t *Data
    cdef readonly FieldConverter_class _converter

    cdef set_slice_unsafe(self, Py_ssize_t i, Matrix_gfpn_dense S) noexcept
    cdef inline int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j) noexcept
    cpdef Matrix_gfpn_dense get_slice(self, Py_ssize_t i, Py_ssize_t j) noexcept
    cpdef list _rowlist_(self, i, j=*) noexcept
    cpdef Matrix_gfpn_dense _multiply_classical(Matrix_gfpn_dense self, Matrix_gfpn_dense right) noexcept
    cpdef Matrix_gfpn_dense _multiply_strassen(Matrix_gfpn_dense self, Matrix_gfpn_dense right, cutoff=*) noexcept

cdef Matrix_gfpn_dense new_mtx(Matrix_t* mat, Matrix_gfpn_dense template) noexcept
