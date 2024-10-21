# sage_setup: distribution = sagemath-modules
"""
Matrix windows
"""

#*****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.tuple cimport *

#########################################################################
# Generic matrix windows, which are used for block echelon and strassen #
# algorithms.                                                           #
#########################################################################
cdef class MatrixWindow:
    ############################
    # creation and initialization
    ############################

    cpdef MatrixWindow new_matrix_window(MatrixWindow self, Matrix matrix,
                                         Py_ssize_t row, Py_ssize_t col,
                                        Py_ssize_t n_rows, Py_ssize_t n_cols):
        """
        This method is here only to provide a fast cdef way of
        constructing new matrix windows. The only implicit assumption
        is that self._matrix and matrix are over the same base ring
        (so share the zero).
        """
        cdef type t = type(self)
        cdef MatrixWindow M = <MatrixWindow>t.__new__(t)
        M._matrix = matrix
        M._row = row
        M._col = col
        M._nrows = n_rows
        M._ncols = n_cols
        M._cached_zero = self._cached_zero
        return M

    def __init__(MatrixWindow self, Matrix matrix,
                 Py_ssize_t row, Py_ssize_t col, Py_ssize_t nrows, Py_ssize_t ncols):
        self._matrix = matrix
        self._row = row
        self._col = col
        self._nrows = nrows
        self._ncols = ncols

    cdef object _zero(self):
        if self._cached_zero is None:
            self._cached_zero = self._matrix.base_ring()(0)  # expensive
        return self._cached_zero

    cpdef MatrixWindow matrix_window(MatrixWindow self, Py_ssize_t row, Py_ssize_t col,
                                    Py_ssize_t n_rows, Py_ssize_t n_cols):
        """
        Return a matrix window relative to this window of the
        underlying matrix.
        """
        if row == 0 and col == 0 and n_rows == self._nrows and n_cols == self._ncols:
            return self
        return self.new_matrix_window(self._matrix, self._row + row, self._col + col, n_rows, n_cols)

    cpdef new_empty_window(MatrixWindow self, Py_ssize_t nrows, Py_ssize_t ncols):
        a = self._matrix.new_matrix(nrows, ncols)
        return self.new_matrix_window(a, 0, 0, nrows, ncols)

    def __repr__(self):
        return "Matrix window of size %s x %s at (%s,%s):\n%s"%(
            self._nrows, self._ncols, self._row, self._col, self._matrix)

    ############################
    # Getting and setting entries
    ############################

    def set(MatrixWindow self, MatrixWindow src):
        if self._matrix._parent.base_ring() is not src._matrix._parent.base_ring():
            raise TypeError("Parents must be equal.")
        self.set_to(src)

    cpdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        self._matrix.set_unsafe(i + self._row, j + self._col, x)

    cpdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._matrix.get_unsafe(i + self._row, j + self._col)

    def __setitem__(self, ij, x):
        cdef Py_ssize_t i, j
        if isinstance(ij, tuple):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError("index must be an integer or pair of integers")
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
                raise IndexError("matrix index out of range")
            self.set_unsafe(i, j, x)
        else:
            # If ij is not a tuple, coerce to an integer and set the row.
            i = ij
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, x)

    def __getitem__(self, ij):
        cdef Py_ssize_t i, j

        if isinstance(ij, tuple):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError("index must be an integer or pair of integers")
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
                raise IndexError("matrix index out of range")
            return self.get_unsafe(i, j)
        else:
            # If ij is not a tuple, coerce to an integer and get the row.
            i = ij
            return self.row(i)

    cpdef matrix(MatrixWindow self):
        """
        Return the underlying matrix that this window is a view of.
        """
        return self._matrix

    cpdef to_matrix(MatrixWindow self):
        """
        Return an actual matrix object representing this view.
        """
        cdef MatrixWindow w
        a = self._matrix.new_matrix(self._nrows, self._ncols)  # zero matrix
        w = self.new_matrix_window(a, 0, 0, self._nrows, self._ncols)
        w.set_to(self)
        return a

    def nrows(MatrixWindow self):
        return self._nrows

    def ncols(MatrixWindow self):
        return self._ncols

    cpdef set_to(MatrixWindow self, MatrixWindow A):
        """
        Change self, making it equal A.
        """
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError("incompatible dimensions")
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j))
        return 0

    cpdef set_to_zero(MatrixWindow self):
        cdef Py_ssize_t i, j
        z = self._zero()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, z)

    cpdef add(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError("incompatible dimensions")
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) + A.get_unsafe(i, j))

    cpdef subtract(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError("incompatible dimensions")
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) - A.get_unsafe(i, j))

    cpdef set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError("incompatible dimensions")
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError("incompatible dimensions")
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) + B.get_unsafe(i, j))

    cpdef set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) - B.get_unsafe(i, j))

    cpdef set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError("incompatible dimensions")
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self._zero()
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    cpdef add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError("incompatible dimensions")
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self.get_unsafe(i, j)
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    cpdef subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError("incompatible dimensions")
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self.get_unsafe(i, j)
                for k from 0 <= k < A._ncols:
                    s = s - A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    cpdef swap_rows(MatrixWindow self, Py_ssize_t a, Py_ssize_t b):
        self._matrix.swap_rows_c(self._row + a, self._row + b)

    def echelon_in_place(MatrixWindow self):
        """
        Calculate the echelon form of this matrix, returning the list of pivot columns
        """
        echelon = self.to_matrix()
        echelon.echelonize() # TODO: read only, only need to copy pointers
        self.set_to(echelon.matrix_window())
        return echelon.pivots()

    cpdef bint element_is_zero(MatrixWindow self, Py_ssize_t i, Py_ssize_t j) noexcept:
        return self._matrix.get_unsafe(i+self._row, j+self._col) == self._zero()
