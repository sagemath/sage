"""
Matrices over the $\GF{2}$ via NTL

This class is only provided to have a complete NTL interface and for
comparison purposes. Sage's native matrices over $F_2$ are much faster
for many problems like matrix multiplication and Gaussian elimination.

AUTHORS:
 - Martin Albrecht <malb@informatik.uni-bremen.de>
   2008-09: initial version
"""

#*****************************************************************************
#         Copyright (C) 2005 William Stein <wstein@gmail.com>
#    Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/signals.pxi"
include 'misc.pxi'
include 'decl.pxi'

from cpython.object cimport Py_EQ, Py_NE
from ntl_GF2 cimport ntl_GF2
from sage.rings.integer cimport Integer
from sage.libs.ntl.ntl_ZZ import unpickle_class_args

cdef class ntl_mat_GF2(object):
    r"""
    The \class{mat_GF2} class implements arithmetic with matrices over $F_2$.
    """
    def __init__(self, nrows=0, ncols=0, v=None):
        """
        Constructs a matrix over ntl.GF2.

        INPUT:
            nrows -- number of rows
            ncols -- number of columns
            v     -- either a list or a matrix over GF(2^x)

        EXAMPLES::

            sage: A = ntl.mat_GF2(4,4); A
            [[0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            ]

            sage: A = random_matrix(GF(2),4,4); A
            [0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]

            sage: B = ntl.mat_GF2(A); B
            [[0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]
            ]

            sage: B = ntl.mat_GF2(4, 4, A.list()); B
            [[0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]
            ]
        """
        cdef Py_ssize_t _nrows, _ncols
        cdef Py_ssize_t i, j
        cdef GF2_c _elem

        from sage.matrix.matrix import is_Matrix

        if is_Matrix(nrows):
            _nrows = nrows.nrows()
            _ncols = nrows.ncols()
            v = nrows
            self.x.SetDims(_nrows, _ncols)
            sig_on()
            for i from 0 <= i < _nrows:
                for j from 0 <= j < _ncols:
                    GF2_conv_long(_elem, int(v[i,j])%2)
                    mat_GF2_setitem(&self.x, i, j, &_elem)
            sig_off()
            return

        _nrows = nrows
        _ncols = ncols
        self.x.SetDims(_nrows, _ncols)

        if v is not None:
            sig_on()
            for i from 0 <= i < _nrows:
                for j from 0 <= j < _ncols:
                    elem = v[i*_ncols+j]
                    if not isinstance(elem, ntl_GF2):
                        elem = ntl_GF2(elem)
                    mat_GF2_setitem(&self.x, i, j, &(<ntl_GF2>elem).x)
            sig_off()

    cdef ntl_GF2 _new_element(self):
        cdef ntl_GF2 r
        r = ntl_GF2.__new__(ntl_GF2)
        return r

    cdef ntl_mat_GF2 _new(self):
        cdef ntl_mat_GF2 r
        r = ntl_mat_GF2.__new__(ntl_mat_GF2)
        r.x.SetDims(self.x.NumRows(),self.x.NumCols())
        return r

    def __reduce__(self):
        """
            sage: A = random_matrix(GF(2),4,4)
            sage: B = ntl.mat_GF2(A)
            sage: loads(dumps(B)) == B # indirect doctest
            True
        """
        return unpickle_class_args, (ntl_mat_GF2, (self.x.NumRows(), self.x.NumCols(), self.list()))

    def __repr__(self):
        """
        Return the string representation of this matrix.

        EXAMPLE::

            sage: A = random_matrix(GF(2),4,4)
            sage: B = ntl.mat_GF2(A); B # indirect doctest
            [[0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]
            ]
        """
        return mat_GF2_to_PyString(&self.x)

    def __mul__(ntl_mat_GF2 self, other):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: B = random_matrix(GF(2),4,4)
            sage: ntl.mat_GF2(A)*ntl.mat_GF2(B)
            [[0 0 1 0]
            [1 1 0 1]
            [0 0 0 1]
            [1 1 0 0]
            ]

            sage: A*B
            [0 0 1 0]
            [1 1 0 1]
            [0 0 0 1]
            [1 1 0 0]
        """
        cdef ntl_mat_GF2 r = self._new()
        if not isinstance(other, ntl_mat_GF2):
            other = ntl_mat_GF2(other)
        sig_on()
        mat_GF2_mul(r.x, self.x, (<ntl_mat_GF2>other).x)
        sig_off()
        return r

    def __sub__(ntl_mat_GF2 self, other):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: B = random_matrix(GF(2),4,4)
            sage: ntl.mat_GF2(A) - ntl.mat_GF2(B)
            [[0 1 0 0]
            [0 1 0 0]
            [1 1 1 0]
            [0 1 1 1]
            ]

            sage: A - B
            [0 1 0 0]
            [0 1 0 0]
            [1 1 1 0]
            [0 1 1 1]
        """
        cdef ntl_mat_GF2 r = self._new()
        if not isinstance(other, ntl_mat_GF2):
            other = ntl_mat_GF2(other,)
        sig_on()
        mat_GF2_sub(r.x, self.x, (<ntl_mat_GF2>other).x)
        sig_off()
        return r

    def __add__(ntl_mat_GF2 self, other):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: B = random_matrix(GF(2),4,4)
            sage: ntl.mat_GF2(A) + ntl.mat_GF2(B)
            [[0 1 0 0]
            [0 1 0 0]
            [1 1 1 0]
            [0 1 1 1]
            ]

            sage: A + B
            [0 1 0 0]
            [0 1 0 0]
            [1 1 1 0]
            [0 1 1 1]

        """
        cdef ntl_mat_GF2 r = self._new()
        if not isinstance(other, ntl_mat_GF2):
            other = ntl_mat_GF2(other)
        sig_on()
        mat_GF2_add(r.x, self.x, (<ntl_mat_GF2>other).x)
        sig_off()
        return r

    def __neg__(ntl_mat_GF2 self):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: -ntl.mat_GF2(A)
            [[0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]
            ]

            sage: -A
            [0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]
        """
        cdef ntl_mat_GF2 r = self._new()
        sig_on()
        mat_GF2_negate(r.x, self.x)
        sig_off()
        return r

    def __pow__(ntl_mat_GF2 self, long e, ignored):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: ntl.mat_GF2(A)^0
            [[1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            ]

            sage: A^0
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: ntl.mat_GF2(A)^3
            [[0 1 1 0]
            [0 0 0 0]
            [0 1 1 0]
            [0 1 1 0]
            ]

            sage: A^3
            [0 1 1 0]
            [0 0 0 0]
            [0 1 1 0]
            [0 1 1 0]
        """
        cdef ntl_mat_GF2 r = self._new()
        sig_on()
        mat_GF2_power(r.x, self.x, e)
        sig_off()
        return r

    def __richcmp__(ntl_mat_GF2 self, other, int op):
        """
        Compare self to other.

        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: A1 = ntl.mat_GF2(A)
            sage: A2 = ntl.mat_GF2(A)
            sage: A1 == A2
            True
            sage: A1[0,0] += 1
            sage: A1 == A2
            False
            sage: A1 == "x"
            False
        """
        if op != Py_EQ and op != Py_NE:
            raise TypeError("matrices over GF(2) are not ordered")

        cdef ntl_mat_GF2 b
        try:
            b = <ntl_mat_GF2?>other
        except TypeError:
            return NotImplemented

        return (op == Py_EQ) == (self.x == b.x)

    def NumRows(self):
        """
        Return the number of rows of this matrix.

        EXAMPLES::

            sage: A = ntl.mat_GF2(10,10)
            sage: A.NumRows()
            10
        """
        return int(self.x.NumRows())

    def NumCols(self):
        """
        Return the number of columns of this matrix.

        EXAMPLES::

            sage: A = ntl.mat_GF2(10,10)
            sage: A.NumCols()
            10
        """
        return int(self.x.NumCols())

    def __setitem__(self, ij, x):
        """
        EXAMPLES::

            sage: A = ntl.mat_GF2(5,5)
            sage: A[0,0] = 1
            sage: A[0,2] = 1
            sage: A
            [[1 0 1 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            ]
        """
        cdef int i, j
        if not isinstance(x, ntl_GF2):
            x = ntl_GF2(x)

        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.x.NumCols()==1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = ij
            j = 0
        elif self.x.NumRows()==1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = 0
            j = ij
        else:
            raise TypeError('ij is not a matrix index')

        if i < 0 or i >= self.x.NumRows() or j < 0 or j >= self.x.NumCols():
            raise IndexError("array index out of range")

        mat_GF2_setitem(&self.x, i, j, &(<ntl_GF2>x).x)

    def __getitem__(self, ij):
        """
        EXAMPLE::

            sage: A = ntl.mat_GF2(3,3,range(9))
            sage: A[0,0]
            0
            sage: A[1,2]
            1
        """
        cdef int i, j
        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.x.NumCols() == 1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = ij
            j = 0
        elif self.x.NumRows() == 1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = 0
            j = ij
        else:
            raise TypeError('ij is not a matrix index')

        if i < 0 or i >= self.x.NumRows() or j < 0 or j >= self.x.NumCols():
            raise IndexError("array index out of range")

        cdef ntl_GF2 e = self._new_element()
        e.x = self.x.get( i+1, j+1 )
        return e

    def determinant(self):
        """
        Returns the determinant.

        EXAMPLES::

            sage: A = ntl.mat_GF2(3,3,range(9))
            sage: A.determinant()
            0
            sage: A = ntl.mat_GF2(3,3,[1,0,0, 0,1,0, 0,0,1])
            sage: A.determinant()
            1
        """
        cdef ntl_GF2 r = self._new_element()
        sig_on()
        r.x = mat_GF2_determinant(self.x)
        sig_off()
        return r

    def gauss(self,ncols=-1):
        """
        Performs unitary row operations so as to bring this matrix
        into row echelon form (not reduced!).  If the optional
        argument \code{ncols} is supplied, stops when first ncols
        columns are in echelon form.  The return value is the rank (or
        the rank of the first ncols columns).

        INPUT:
           ncols -- number of columns to process (default: all)

        EXAMPLES::
            sage: A = random_matrix(GF(2), 10, 10)
            sage: Abar = ntl.mat_GF2(A)
            sage: A.echelon_form()
            [1 0 0 0 0 0 1 0 1 0]
            [0 1 0 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 1 0 1 0]
            [0 0 0 1 0 0 1 0 1 0]
            [0 0 0 0 1 0 1 0 0 0]
            [0 0 0 0 0 1 1 0 0 0]
            [0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            sage: A.rank()
            8

            sage: Abar.gauss()
            8

            sage: Abar
            [[1 1 1 1 0 1 0 1 1 0]
            [0 1 1 1 0 1 1 0 0 1]
            [0 0 1 1 1 1 0 0 0 0]
            [0 0 0 1 0 0 1 1 1 1]
            [0 0 0 0 1 1 0 1 0 0]
            [0 0 0 0 0 1 1 1 0 1]
            [0 0 0 0 0 0 0 1 0 1]
            [0 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0]
            ]
        """
        if ncols == -1:
            ncols = self.x.NumCols()
        return int(mat_GF2_gauss(self.x, int(ncols)))

    def list(self):
        """
        Returns a list of the entries in this matrix

        EXAMPLES::
            sage: A = random_matrix(GF(2), 4, 4)
            sage: Abar = ntl.mat_GF2(A)
            sage: A.list()
            [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0]

            sage: Abar.list()
            [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0]

        """
        cdef Py_ssize_t i, j
        return [self[i,j] for i in range(self.NumRows()) for j in range(self.x.NumCols())]

    def IsZero(self):
        r"""
        Return \code{True} if this matrix contains only zeroes, and \code{False} otherwise.

        EXAMPLES::
            sage: A = random_matrix(GF(2), 10, 10)
            sage: Abar = ntl.mat_GF2(A)
            sage: Abar.IsZero()
            False
            sage: Abar = ntl.mat_GF2(10,10)
            sage: Abar.IsZero()
            True
        """
        cdef long isZero
        sig_on()
        isZero = mat_GF2_IsZero(self.x)
        sig_off()
        return bool(isZero)

    def _sage_(ntl_mat_GF2 self):
        r"""
        Returns a \class{Matrix} over GF(2).

        EXAMPLES::

            sage: A = random_matrix(GF(2), 6, 6); A
            [0 1 0 1 1 0]
            [0 1 1 1 0 1]
            [0 0 0 1 0 1]
            [0 1 1 0 0 1]
            [0 0 0 1 1 1]
            [0 0 1 1 1 1]

            sage: Abar = ntl.mat_GF2(A); Abar
            [[0 1 0 1 1 0]
            [0 1 1 1 0 1]
            [0 0 0 1 0 1]
            [0 1 1 0 0 1]
            [0 0 0 1 1 1]
            [0 0 1 1 1 1]
            ]

            sage: Abar._sage_()
            [0 1 0 1 1 0]
            [0 1 1 1 0 1]
            [0 0 0 1 0 1]
            [0 1 1 0 0 1]
            [0 0 0 1 1 1]
            [0 0 1 1 1 1]
        """
        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        from sage.matrix.constructor import Matrix
        m =  Matrix(FiniteField(2),self.x.NumRows(),self.x.NumCols())

        cdef Py_ssize_t i, j

        for i from 0 <= i < self.x.NumRows():
            for j from 0 <= j < self.x.NumCols():
                m[i,j] = GF2_conv_to_long(self.x.get( i+1, j+1))
        return m

    def transpose(ntl_mat_GF2 self):
        """
        Returns the transposed matrix of this matrix.

        EXAMPLES::
            sage: A = random_matrix(GF(2), 10, 10)
            sage: Abar = ntl.mat_GF2(A); Abar
            [[0 1 0 1 1 0 0 0 1 1]
            [0 1 1 1 0 1 1 0 0 1]
            [0 0 0 1 0 1 0 0 1 0]
            [0 1 1 0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1 0 1 1]
            [0 0 1 1 1 1 0 0 0 0]
            [1 1 1 1 0 1 0 1 1 0]
            [0 0 0 1 1 0 0 0 1 1]
            [1 0 0 0 1 1 1 0 1 1]
            [1 0 0 1 1 0 1 0 0 0]
            ]

            sage: Abar.transpose()
            [[0 0 0 0 0 0 1 0 1 1]
            [1 1 0 1 0 0 1 0 0 0]
            [0 1 0 1 0 1 1 0 0 0]
            [1 1 1 0 1 1 1 1 0 1]
            [1 0 0 0 1 1 0 1 1 1]
            [0 1 1 1 1 1 1 0 1 0]
            [0 1 0 0 1 0 0 0 1 1]
            [0 0 0 1 0 0 1 0 0 0]
            [1 0 1 1 1 0 1 1 1 0]
            [1 1 0 0 1 0 0 1 1 0]
            ]
        """
        cdef ntl_mat_GF2 r = self._new()
        sig_on()
        mat_GF2_transpose(r.x, self.x)
        sig_off()
        return r

    def __invert__(self):
        """
        Return $X = A^{-1}$; an error is raised if A is singular.

        EXAMPLES::
            sage: l = [0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, \
                       0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, \
                       1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, \
                       0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0]
            sage: A = ntl.mat_GF2(8,8,l)
            sage: E = ~A*A
            sage: E.IsIdent()
            True
        """
        cdef ntl_mat_GF2 r = self._new()
        sig_on()
        mat_GF2_inv(r.x, self.x)
        sig_off()
        return r

    def IsIdent(self, n = -1):
        """
        test if this matrix is the n x n identity matrix.

        EXAMPLES::
            sage: A = ntl.mat_GF2(4,4)
            sage: A[0,0] = 1
            sage: A[1,1] = 1
            sage: A[2,2] = 1
            sage: A.IsIdent()
            False
            sage: A[3,3] = 1
            sage: A.IsIdent()
            True
        """
        if n < 0:
            n = self.NumRows()
        return bool(mat_GF2_IsIdent(self.x, n))

    def IsDiag(self, long n, ntl_GF2 d):
        """
        test if X is an  n x n diagonal matrix with d on diagonal.

        EXAMPLES::
            sage: A = ntl.mat_GF2(4,4)
            sage: A[0,0] = 1
            sage: A[1,1] = 1
            sage: A[2,2] = 1
            sage: A.IsDiag(3, ntl.GF2(1))
            False
            sage: A[3,3] = 1
            sage: A.IsDiag(4, ntl.GF2(1))
            True
        """
        return bool(mat_GF2_IsDiag(self.x, n, d.x))

    def image(self):
        """
        If A is this matrix and X the matrix returned by this function
        then, the rows of X are computed as basis of A's row space.  X
        is in row echelon form.

        EXAMPLE::
            sage: A = random_matrix(GF(2),10,10)
            sage: Abar = ntl.mat_GF2(A)
            sage: A.image()
            Vector space of degree 10 and dimension 8 over Finite Field of size 2
            Basis matrix:
            [1 0 0 0 0 0 1 0 1 0]
            [0 1 0 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 1 0 1 0]
            [0 0 0 1 0 0 1 0 1 0]
            [0 0 0 0 1 0 1 0 0 0]
            [0 0 0 0 0 1 1 0 0 0]
            [0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 1]


            sage: Abar.image()
            [[1 1 1 1 0 1 0 1 1 0]
            [0 1 1 1 0 1 1 0 0 1]
            [0 0 1 1 1 1 0 0 0 0]
            [0 0 0 1 0 0 1 1 1 1]
            [0 0 0 0 1 1 0 1 0 0]
            [0 0 0 0 0 1 1 1 0 1]
            [0 0 0 0 0 0 0 1 0 1]
            [0 0 0 0 0 0 0 0 0 1]
            ]
        """
        cdef ntl_mat_GF2 X = self._new()
        sig_on()
        mat_GF2_image(X.x, self.x)
        sig_off()
        return X

    def kernel(self):
        """
        Computes a basis for the kernel of the map x -> x*A. where x
        is a row vector.

        EXAMPLE::

            sage: A = random_matrix(GF(2),10,10)
            sage: Abar = ntl.mat_GF2(A)
            sage: A.kernel()
            Vector space of degree 10 and dimension 2 over Finite Field of size 2
            Basis matrix:
            [1 1 1 0 1 1 0 1 0 0]
            [0 0 0 1 1 0 1 0 1 0]
            sage: Abar.kernel()
            [[0 0 0 1 1 0 1 0 1 0]
            [1 1 1 0 1 1 0 1 0 0]
            ]
        """
        cdef ntl_mat_GF2 X = self._new()
        sig_on()
        mat_GF2_kernel(X.x, self.x)
        sig_off()
        return X
