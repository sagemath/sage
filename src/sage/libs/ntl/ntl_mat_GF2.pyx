# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

r"""
Matrices over the `\GF{2}` via NTL

This class is only provided to have a complete NTL interface and for
comparison purposes. Sage's native matrices over `F_2` are much faster
for many problems like matrix multiplication and Gaussian elimination.

AUTHORS:

- Martin Albrecht <malb@informatik.uni-bremen.de> 2008-09: initial version
"""

# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from sage.ext.cplusplus cimport ccrepr

include 'misc.pxi'
include 'decl.pxi'

from cpython.object cimport Py_EQ, Py_NE
from sage.libs.ntl.ntl_GF2 cimport ntl_GF2
from sage.rings.integer cimport Integer
from sage.libs.ntl.ntl_ZZ import unpickle_class_args


cdef class ntl_mat_GF2():
    r"""
    The \class{mat_GF2} class implements arithmetic with matrices over `F_2`.
    """
    def __init__(self, nrows=0, ncols=0, v=None):
        """
        Construct a matrix over ntl.GF2.

        INPUT:

        - ``nrows`` -- number of rows
        - ``ncols`` -- number of columns
        - ``v`` -- either a list or a matrix over GF(2^x)

        EXAMPLES::

            sage: A = ntl.mat_GF2(4,4); A
            [[0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            ]

            sage: A = random_matrix(GF(2),4,4); A  # random
            [0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]

            sage: B = ntl.mat_GF2(A); B  # random
            [[0 1 0 1]
            [0 1 1 1]
            [0 0 0 1]
            [0 1 1 0]
            ]

            sage: B = ntl.mat_GF2(4, 4, A.list())
            sage: B == A
            True
        """
        cdef Py_ssize_t _nrows, _ncols
        cdef Py_ssize_t i, j
        cdef GF2_c _elem

        from sage.structure.element import Matrix

        if isinstance(nrows, Matrix):
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
        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: B = ntl.mat_GF2(A)
            sage: loads(dumps(B)) == B # indirect doctest
            True
        """
        return unpickle_class_args, (ntl_mat_GF2, (self.x.NumRows(), self.x.NumCols(), self.list()))

    def __repr__(self):
        """
        Return the string representation of this matrix.

        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: B = ntl.mat_GF2(A)
            sage: B.__repr__()[1:-2] == A.__repr__()
            True
        """
        return ccrepr(self.x)

    def __mul__(ntl_mat_GF2 self, other):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),4,4)
            sage: B = random_matrix(GF(2),4,4)
            sage: c = ntl.mat_GF2(A)*ntl.mat_GF2(B)
            sage: c._sage_() == A*B
            True
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
            sage: c = ntl.mat_GF2(A) - ntl.mat_GF2(B)
            sage: c._sage_() == A - B
            True
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
            sage: c = ntl.mat_GF2(A) + ntl.mat_GF2(B)
            sage: c._sage_() == A + B
            True
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
            sage: (-ntl.mat_GF2(A))._sage_() == -A
            True
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
            sage: Abar = ntl.mat_GF2(A)
            sage: (Abar^0)._sage_() == A^0
            True
            sage: (Abar^1)._sage_() == A^1
            True
            sage: (Abar^2)._sage_() == A^2
            True
            sage: (Abar^3)._sage_() == A^3
            True
        """
        cdef ntl_mat_GF2 r = self._new()
        sig_on()
        mat_GF2_power(r.x, self.x, e)
        sig_off()
        return r

    def __richcmp__(ntl_mat_GF2 self, other, int op):
        """
        Compare ``self`` to ``other``.

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
        EXAMPLES::

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
        Return the determinant.

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

    def gauss(self, ncols=-1):
        r"""
        Perform unitary row operations so as to bring this matrix
        into row echelon form (not reduced!).

        If the optional argument ``ncols`` is supplied, stops when
        first ``ncols`` columns are in echelon form.  The return value is
        the rank (or the rank of the first ``ncols`` columns).

        INPUT:

        - ``ncols`` -- number of columns to process (default: all)

        EXAMPLES::

            sage: A = random_matrix(GF(2), 10, 10)
            sage: Abar = ntl.mat_GF2(A)
            sage: A.rank() == Abar.gauss()
            True
            sage: Abar  # random
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

        ``Abar`` is in row echelon form now::

            sage: first_nonzero_indices = [Abar._sage_().row(i).nonzero_positions()[0] for i in range(A.rank())]
            sage: all(first_nonzero_indices[i] < first_nonzero_indices[i+1] for i in range(A.rank()-1))
            True

        ``Abar`` is not reduced::

            sage: all(Abar._sage_().row(i).nonzero_positions() == [] for i in range(A.rank(), Abar.NumRows()))
            True
        """
        if ncols == -1:
            ncols = self.x.NumCols()
        return int(mat_GF2_gauss(self.x, int(ncols)))

    def list(self):
        """
        Return a list of the entries in this matrix.

        EXAMPLES::

            sage: A = random_matrix(GF(2), 4, 4)
            sage: Abar = ntl.mat_GF2(A)
            sage: A.list() == Abar.list()
            True
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
        Return a \class{Matrix} over GF(2).

        EXAMPLES::

            sage: A = random_matrix(GF(2), 6, 6)
            sage: Abar = ntl.mat_GF2(A)
            sage: Abar._sage_() == A
            True
        """
        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        from sage.matrix.constructor import matrix
        m = matrix(FiniteField(2),self.x.NumRows(),self.x.NumCols())

        cdef Py_ssize_t i, j

        for i from 0 <= i < self.x.NumRows():
            for j from 0 <= j < self.x.NumCols():
                m[i,j] = GF2_conv_to_long(self.x.get( i+1, j+1))
        return m

    def transpose(ntl_mat_GF2 self):
        """
        Return the transposed matrix of this matrix.

        EXAMPLES::

            sage: A = random_matrix(GF(2), 10, 10)
            sage: Abar = ntl.mat_GF2(A)
            sage: Abar_t = Abar.transpose()
            sage: A_t = A.transpose()
            sage: A_t == Abar_t._sage_()
            True
        """
        cdef ntl_mat_GF2 r = self._new()
        sig_on()
        mat_GF2_transpose(r.x, self.x)
        sig_off()
        return r

    def __invert__(self):
        """
        Return `X = A^{-1}`; an error is raised if A is singular.

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

    def IsIdent(self, n=-1):
        """
        Test if this matrix is the n x n identity matrix.

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
        Test if X is an  n x n diagonal matrix with d on diagonal.

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
        then, the rows of X are computed as basis of A's row space.
        X is in row echelon form.

        EXAMPLES::

            sage: A = random_matrix(GF(2),10,10)
            sage: Abar = ntl.mat_GF2(A)
            sage: A_image = A.image().matrix()
            sage: Abar_image =  Abar.image()._sage_()
            sage: A_image.row_space() == Abar_image.row_space()
            True

        X is in row echelon form::

            sage: first_nonzero_indices = [row.nonzero_positions()[0] for row in Abar_image.rows()]
            sage: all(first_nonzero_indices[i] < first_nonzero_indices[i+1] for i in range(Abar_image.nrows() - 1))
            True
        """
        cdef ntl_mat_GF2 X = self._new()
        sig_on()
        mat_GF2_image(X.x, self.x)
        sig_off()
        return X

    def kernel(self):
        """
        Compute a basis for the kernel of the map x -> x*A. where x
        is a row vector.

        EXAMPLES::

            sage: A = random_matrix(GF(2),10,10)
            sage: Abar = ntl.mat_GF2(A)
            sage: K_abar = Abar.kernel(); K_abar  # random
            [[0 0 0 1 1 0 1 0 1 0]
            [1 1 1 0 1 1 0 1 0 0]
            ]
            sage: (K_abar*Abar).IsZero()
            True
            sage: K_a = A.kernel().matrix()
            sage: K_a.row_space() == K_abar._sage_().row_space()
            True
        """
        cdef ntl_mat_GF2 X = self._new()
        sig_on()
        mat_GF2_kernel(X.x, self.x)
        sig_off()
        return X
