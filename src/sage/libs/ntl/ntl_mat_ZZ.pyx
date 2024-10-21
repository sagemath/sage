# sage_setup: distribution = sagemath-ntl
# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

from cysignals.signals cimport sig_on, sig_off
from sage.ext.cplusplus cimport ccrepr

include 'misc.pxi'
include 'decl.pxi'

from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from cpython.object cimport PyObject_RichCompare

from sage.libs.ntl.ntl_ZZ import unpickle_class_args

cdef inline ntl_ZZ make_ZZ(ZZ_c* x):
    cdef ntl_ZZ y
    y = ntl_ZZ()
    y.x = x[0]
    del x
    return y

# You must do sig_on() before calling this function
cdef inline ntl_ZZ make_ZZ_sig_off(ZZ_c* x):
    cdef ntl_ZZ y = make_ZZ(x)
    sig_off()
    return y

cdef inline ntl_mat_ZZ make_mat_ZZ(mat_ZZ_c* x):
    cdef ntl_mat_ZZ y
    y = ntl_mat_ZZ.__new__(ntl_mat_ZZ)
    y.x = x[0]
    del x
    y.__nrows = y.x.NumRows()
    y.__ncols = y.x.NumCols()
    return y

# You must do sig_on() before calling this function
cdef inline ntl_mat_ZZ make_mat_ZZ_sig_off(mat_ZZ_c* x):
    cdef ntl_mat_ZZ y = make_mat_ZZ(x)
    sig_off()
    return y


##############################################################################
#
# ntl_mat_ZZ: Matrices over the integers via NTL
#
##############################################################################

cdef class ntl_mat_ZZ():
    # see ntl_mat_ZZ.pxd for data members
    r"""
    The \class{mat_ZZ} class implements arithmetic with matrices over `\Z`.
    """
    def __init__(self, nrows=0,  ncols=0, v=None):
        r"""
        The \class{mat_ZZ} class implements arithmetic with matrices over `\Z`.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(3,3) ; M
            [
            [0 0 0]
            [0 0 0]
            [0 0 0]
            ]
            sage: ntl.mat_ZZ(3,3,[1..9])
            [
            [1 2 3]
            [4 5 6]
            [7 8 9]
            ]
        """
        cdef unsigned long i, j
        cdef ntl_ZZ tmp
        if nrows == 0 and ncols == 0:
            return
        nrows = int(nrows)
        ncols = int(ncols)
        self.x.SetDims(nrows, ncols)
        self.__nrows = nrows
        self.__ncols = ncols
        if v is not None:
            for i from 0 <= i < nrows:
                for j from 0 <= j < ncols:
                    tmp = ntl_ZZ(v[i*ncols+j])
                    mat_ZZ_setitem(&self.x, i, j, &tmp.x)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: m = ntl.mat_ZZ(3, 2, range(6)); m
            [
            [0 1]
            [2 3]
            [4 5]
            ]
            sage: loads(dumps(m))
            [
            [0 1]
            [2 3]
            [4 5]
            ]
            sage: loads(dumps(m)) == m
            True
        """
        return unpickle_class_args, (ntl_mat_ZZ, (self.__nrows, self.__ncols, self.list()))

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,3,[5..10]) ; M.__repr__()
            '[\n[5 6 7]\n[8 9 10]\n]'
        """
        return ccrepr(self.x).replace('[[','[\n[',1)

    def __mul__(ntl_mat_ZZ self, other):
        """
        Multiply two matrices.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,2,[8..11]) ; N = ntl.mat_ZZ(2,2,[-1..2])
            sage: M*N
            [
            [1 18]
            [1 22]
            ]
        """
        cdef ntl_mat_ZZ r = ntl_mat_ZZ.__new__(ntl_mat_ZZ)
        if not isinstance(self, ntl_mat_ZZ):
            self = ntl_mat_ZZ(self)
        if not isinstance(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        sig_on()
        mat_ZZ_mul(r.x, (<ntl_mat_ZZ>self).x, (<ntl_mat_ZZ>other).x)
        sig_off()
        return r

    def __sub__(ntl_mat_ZZ self, other):
        """
        Return ``self - other``.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,2,[8..11]) ; N = ntl.mat_ZZ(2,2,[-1..2])
            sage: M-N
            [
            [9 9]
            [9 9]
            ]
        """
        cdef ntl_mat_ZZ r = ntl_mat_ZZ.__new__(ntl_mat_ZZ)
        if not isinstance(self, ntl_mat_ZZ):
            self = ntl_mat_ZZ(self)
        if not isinstance(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        sig_on()
        mat_ZZ_sub(r.x, (<ntl_mat_ZZ>self).x, (<ntl_mat_ZZ>other).x)
        sig_off()
        return r

    def __add__(ntl_mat_ZZ self, other):
        """
        Return ``self + other``.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,2,[8..11]) ; N = ntl.mat_ZZ(2,2,[-1..2])
            sage: M+N
            [
            [7 9]
            [11 13]
            ]
        """
        cdef ntl_mat_ZZ r = ntl_mat_ZZ.__new__(ntl_mat_ZZ)
        if not isinstance(self, ntl_mat_ZZ):
            self = ntl_mat_ZZ(self)
        if not isinstance(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        sig_on()
        mat_ZZ_add(r.x, (<ntl_mat_ZZ>self).x, (<ntl_mat_ZZ>other).x)
        sig_off()
        return r

    def __richcmp__(ntl_mat_ZZ self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,2,[3..6])
            sage: N = M^2
            sage: M == M
            True
            sage: M == N
            False
            sage: M == 0
            False
            sage: M < ntl.mat_ZZ(2,2,[4,4,4,4])
            True
            sage: M != ntl.mat_ZZ(1,1,[3])
            True
            sage: M != 0
            True
        """
        cdef ntl_mat_ZZ b
        try:
            b = <ntl_mat_ZZ?>other
        except TypeError:
            return NotImplemented
        return PyObject_RichCompare(self.list(), other.list(), op)

    def __pow__(ntl_mat_ZZ self, long e, ignored):
        """
        Return ``self`` to the e power.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,2,[8..11])
            sage: M**4
            [
            [56206 62415]
            [69350 77011]
            ]
            sage: M**0
            [
            [1 0]
            [0 1]
            ]
            sage: M**(-1)
            Traceback (most recent call last):
            ...
            ValueError: cannot take negative powers of matrices.
        """
        if self.__nrows != self.__ncols:
            raise TypeError("cannot take powers of non-square matrices.")
        if e < 0:
            raise ValueError("cannot take negative powers of matrices.")
        cdef ntl_mat_ZZ r = ntl_mat_ZZ.__new__(ntl_mat_ZZ)
        sig_on()
        mat_ZZ_power(r.x, (<ntl_mat_ZZ>self).x, e)
        sig_off()
        return r

    def nrows(self):
        """
        Return the number of rows in ``self``.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(5,5,range(25))
            sage: M.nrows()
            5
        """
        return self.__nrows

    def ncols(self):
        """
        Return the number of columns in ``self``.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(5,8,range(40))
            sage: M.ncols()
            8
        """
        return self.__ncols

    def __setitem__(self, ij, x):
        """
        Given a tuple (i, j), return self[i,j].

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,9,[3..20])
            sage: M[1,7] ## indirect doctest
            19
        """
        cdef ntl_ZZ y
        cdef int i, j
        if not isinstance(x, ntl_ZZ):
            y = ntl_ZZ(x)
        else:
            y = x
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError('ij must be a 2-tuple')
        i, j = int(ij[0]),int(ij[1])
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError("array index out of range")
        sig_on()
        mat_ZZ_setitem(&self.x, i, j, &y.x)
        sig_off()

    def __getitem__(self, ij):
        """
        EXAMPLES::

            sage: m = ntl.mat_ZZ(3, 2, range(6))
            sage: m[0,0] ## indirect doctest
            0
            sage: m[2,1]
            5
            sage: m[3,2] #  oops, 0 based
            Traceback (most recent call last):
            ...
            IndexError: array index out of range
        """
        cdef int i, j
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError('ij must be a 2-tuple')
        i, j = ij
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError("array index out of range")
        sig_on()
        return make_ZZ_sig_off(mat_ZZ_getitem(&self.x, i+1, j+1))

    def list(self):
        """
        EXAMPLES::

            sage: m = ntl.mat_ZZ(3, 4, range(12)); m
            [
            [0 1 2 3]
            [4 5 6 7]
            [8 9 10 11]
            ]
            sage: m.list()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        """
        cdef int i, j
        sig_on()
        L = [make_ZZ(mat_ZZ_getitem(&self.x, i+1, j+1))
                    for i from 0 <= i < self.__nrows
                        for j from 0 <= j < self.__ncols]
        sig_off()
        return L

    def determinant(self, deterministic=True):
        """
        Return the determinant of ``self``.

        EXAMPLES::

            sage: ntl.mat_ZZ(8,8,[3..66]).determinant()
            0
            sage: ntl.mat_ZZ(2,5,range(10)).determinant()
            Traceback (most recent call last):
            ...
            TypeError: cannot take determinant of non-square matrix.
            sage: ntl.mat_ZZ(4,4,[next_prime(2**i) for i in range(16)]).determinant()
            -10248
            sage: ntl.mat_ZZ(4,4,[ ZZ.random_element() for _ in range(16) ]).determinant()  # random
            678
        """
        if self.__nrows != self.__ncols:
            raise TypeError("cannot take determinant of non-square matrix.")
        sig_on()
        return make_ZZ_sig_off(mat_ZZ_determinant(&self.x, deterministic))

    def HNF(self, D=None):
        r"""
        The input matrix A=self is an n x m matrix of rank m (so n >=
        m), and D is a multiple of the determinant of the lattice L
        spanned by the rows of A.  The output W is the Hermite Normal
        Form of A; that is, W is the unique m x m matrix whose rows
        span L, such that

        - W is lower triangular,
        - the diagonal entries are positive,
        - any entry below the diagonal is a nonnegative number
          strictly less than the diagonal entry in its column.

        This is implemented using the algorithm of [P. Domich,
        R. Kannan and L. Trotter, Math. Oper. Research 12:50-59,
        1987].

        TIMINGS:

        NTL is not very good compared to MAGMA, unfortunately::

            sage: a = MatrixSpace(ZZ,200).random_element(x=-2, y=2)    # -2 to 2
            sage: A = ntl.mat_ZZ(200,200)
            sage: for i in range(a.nrows()):
            ....:     for j in range(a.ncols()):
            ....:         A[i,j] = a[i,j]
            sage: t = cputime(); d = A.determinant()
            sage: cputime(t)          # random
            0.33201999999999998
            sage: t = cputime(); B = A.HNF(d)  # long time (5s on sage.math, 2011)
            sage: cputime(t)          # random
            6.4924050000000006

        In comparison, MAGMA does this much more quickly:
        \begin{verbatim}
            > A := MatrixAlgebra(IntegerRing(),200)![Random(-2,2) : i in [1..200^2]];
            > time d := Determinant(A);
            Time: 0.140
            > time H := HermiteForm(A);
            Time: 0.290
        \end{verbatim}

        Also, PARI is also faster than NTL if one uses the flag 1 to
        the mathnf routine.  The above takes 16 seconds in PARI.

        TESTS::

            sage: ntl.mat_ZZ(2,2,[1..4]).HNF()
            [
            [1 0]
            [0 2]
            ]
        """
        cdef ntl_ZZ _D
        if D is None:
            _D = self.determinant()
        else:
            _D = ntl_ZZ(D)
        sig_on()
        return make_mat_ZZ_sig_off(mat_ZZ_HNF(&self.x, &_D.x))

    def charpoly(self):
        """
        Find the characteristic polynomial of self, and return it
        as an NTL ZZX.

        EXAMPLES::

            sage: M = ntl.mat_ZZ(2,2,[1,2,3,4])
            sage: M.charpoly()
            [-2 -5 1]
            sage: type(_)
            <class 'sage.libs.ntl.ntl_ZZX.ntl_ZZX'>
            sage: M.determinant()
            -2
        """
        cdef ntl_ZZX r = ntl_ZZX()
        sig_on()
        mat_ZZ_CharPoly(r.x, self.x)
        sig_off()
        return r

    def BKZ_FP(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.BKZ_FP(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.BKZ_FP(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_BKZ_FP(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_BKZ_FP_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def BKZ_QP(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.BKZ_QP(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.BKZ_QP(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_BKZ_QP(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_BKZ_QP_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def BKZ_QP1(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.BKZ_QP1(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.BKZ_QP1(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_BKZ_QP1(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_BKZ_QP1_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def BKZ_XD(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.BKZ_XD(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.BKZ_XD(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_BKZ_XD(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_BKZ_XD_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def BKZ_RR(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.BKZ_RR(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.BKZ_RR(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_BKZ_RR(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_BKZ_RR_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def G_BKZ_FP(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.G_BKZ_FP(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.G_BKZ_FP(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_G_BKZ_FP(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_G_BKZ_FP_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def G_BKZ_QP(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.G_BKZ_QP(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.G_BKZ_QP(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_G_BKZ_QP(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_G_BKZ_QP_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def G_BKZ_QP1(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.G_BKZ_QP1(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.G_BKZ_QP1(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_G_BKZ_QP1(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_G_BKZ_QP1_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def G_BKZ_XD(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.G_BKZ_XD(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.G_BKZ_XD(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_G_BKZ_XD(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_G_BKZ_XD_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def G_BKZ_RR(self, U=None, delta=0.99, BlockSize=10, prune=0, verbose=False):
        r"""
        All BKZ methods are equivalent to the LLL routines,
        except that Block Korkin-Zolotarev reduction is applied. We
        describe here only the differences in the calling syntax.

        * The optional parameter "BlockSize" specifies the size of the
        blocks in the reduction. High values yield shorter vectors,
        but the running time increases exponentially with BlockSize.
        BlockSize should be between 2 and the number of rows of B.

        * The optional parameter "prune" can be set to any positive
        number to invoke the Volume Heuristic from [Schnorr and
        Horner, Eurocrypt '95].  This can significantly reduce the
        running time, and hence allow much bigger block size, but the
        quality of the reduction is of course not as good in general.
        Higher values of prune mean better quality, and slower running
        time.  When prune == 0, pruning is disabled.  Recommended
        usage: for BlockSize >= 30, set 10 <= prune <= 15.

        * The QP1 variant uses quad_float precision to compute
        Gram-Schmidt, but uses double precision in the search phase
        of the block reduction algorithm.  This seems adequate for
        most purposes, and is faster than QP, which uses quad_float
        precision uniformly throughout.

        INPUT:

        - ``U`` -- permutation matrix (see LLL, default: ``None``)
        - ``delta`` -- reduction parameter (default: 0.99)
        - ``BlockSize`` -- see above (default: 10)
        - ``prune`` -- see above (default: 0)
        - ``verbose`` -- print verbose output (default: ``False``)

        EXAMPLES::

            sage: A = Matrix(ZZ,5,5,range(25))
            sage: a = A._ntl_()
            sage: a.G_BKZ_RR(); a
            2
            [
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 1 2 3 4]
            [5 3 1 -1 -3]
            ]

            sage: U = ntl.mat_ZZ(2,2) # note that the dimension doesn't matter
            sage: r = a.G_BKZ_RR(U=U); U
            [
            [0 1 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ]
        """
        if U is None:
            sig_on()
            rank = mat_ZZ_G_BKZ_RR(self.x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        elif isinstance(U, ntl_mat_ZZ):
            sig_on()
            rank = mat_ZZ_G_BKZ_RR_U(self.x, (<ntl_mat_ZZ>U).x, float(delta), int(BlockSize), int(prune), 0, int(verbose))
            sig_off()
            return rank
        else:
            raise TypeError("parameter U has wrong type.")

    def LLL(self, a=3, b=4, return_U=False, verbose=False):
        r"""
        Perform LLL reduction of ``self`` (puts \code{self} in an LLL form).

        \code{self} is an `m x n` matrix, viewed as `m` rows of
        `n`-vectors.  `m` may be less than, equal to, or greater than `n`,
        and the rows need not be linearly independent. ``self`` is
        transformed into an LLL-reduced basis, and the return value is
        the rank r of ``self`` so as det2 (see below).  The first `m-r` rows
        of ``self`` are zero.

        More specifically, elementary row transformations are
        performed on \code{self} so that the nonzero rows of
        new-\code{self} form an LLL-reduced basis for the lattice
        spanned by the rows of old-\code{self}.  The default reduction
        parameter is `\delta=3/4`, which means that the squared length
        of the first nonzero basis vector is no more than `2^{r-1}`
        times that of the shortest vector in the lattice.

        det2 is calculated as the \emph{square} of the determinant of
        the lattice---note that sqrt(det2) is in general an integer
        only when r = n.

        If return_U is True, a value U is returned which is the
        transformation matrix, so that U is a unimodular m x m matrix
        with U * old-\code{self} = new-\code{self}. Note that the
        first m-r rows of U form a basis (as a lattice) for the kernel
        of old-B.

        The parameters a and b allow an arbitrary reduction parameter
        `\delta=a/b`, where `1/4 < a/b \leq 1`, where a and b are positive
        integers.  For a basis reduced with parameter delta, the
        squared length of the first nonzero basis vector is no more
        than `1/(delta-1/4)^{r-1}` times that of the shortest vector in
        the lattice.

        The algorithm employed here is essentially the one in Cohen's
        book: [H. Cohen, A Course in Computational Algebraic Number
        Theory, Springer, 1993]

        INPUT:

        - ``a`` -- parameter a as described above (default: 3)
        - ``b`` -- parameter b as described above (default: 4)
        - ``return_U`` -- return U as described above
        - ``verbose`` -- if ``True`` NTL will produce some verbatim messages on
          what's going on internally (default: ``False``)

        OUTPUT:

        (rank,det2,[U]) where rank,det2, and U are as described
        above and U is an optional return value if return_U is ``True``.

        EXAMPLES::

            sage: M=ntl.mat_ZZ(3,3,[1,2,3,4,5,6,7,8,9])
            sage: M.LLL()
            (2, 54)
            sage: M
            [
            [0 0 0]
            [2 1 0]
            [-1 1 3]
            ]
            sage: M=ntl.mat_ZZ(4,4,[-6,9,-15,-18,4,-6,10,12,10,-16,18,35,-24,36,-46,-82]); M
            [
            [-6 9 -15 -18]
            [4 -6 10 12]
            [10 -16 18 35]
            [-24 36 -46 -82]
            ]
            sage: M.LLL()
            (3, 19140)
            sage: M
            [
            [0 0 0 0]
            [0 -2 0 0]
            [-2 1 -5 -6]
            [0 -1 -7 5]
            ]

        WARNING: This method modifies \code{self}. So after applying
        this method your matrix will be a vector of vectors.
        """
        cdef ZZ_c *det2
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_LLL_U(&det2, &self.x, &U.x, int(a), int(b), int(verbose)))

            return rank, make_ZZ_sig_off(det2), U
        else:
            sig_on()
            rank = int(mat_ZZ_LLL(&det2,&self.x,int(a),int(b),int(verbose)))
            return rank, make_ZZ_sig_off(det2)

    def LLL_FP(self, delta=0.75, return_U=False, verbose=False):
        r"""
        Perform approximate LLL reduction of \code{self} (puts
        \code{self} in an LLL form) subject to the following
        conditions:

        The precision is double.

        The return value is the rank of B.

        Classical Gram-Schmidt Orthogonalization is used:

        This choice uses classical methods for computing the
        Gram-Schmidt orthogonalization.  It is fast but prone to
        stability problems.  This strategy was first proposed by
        Schnorr and Euchner [C. P. Schnorr and M. Euchner,
        Proc. Fundamentals of Computation Theory, LNCS 529, pp. 68-85,
        1991].  The version implemented here is substantially
        different, improving both stability and performance.

        If return_U is True, then also U is returned which is
        the transition matrix: `U * self_{old} = self_{new}`

        The optional argument 'delta' is the reduction parameter, and
        may be set so that 0.50 <= delta < 1.  Setting it close to 1
        yields shorter vectors, and also improves the stability, but
        increases the running time.  Recommended value: delta =
        0.99.

        The optional parameter 'verbose' can be set to see all kinds
        of fun things printed while the routine is executing.  A
        status report is also printed every once in a while.

        INPUT:

        - ``delta`` -- as described above (0.5 <= delta < 1.0) (default: 0.75)
        - ``return_U`` -- return U as described above
        - ``verbose`` -- if ``True`` NTL will produce some verbatim messages on
          what's going on internally (default: ``False``)

        OUTPUT:

        (rank,[U]) where rank and U are as described above and U
        is an optional return value if ``return_U`` is ``True``.

        EXAMPLES::

            sage: M=ntl.mat_ZZ(3,3,[1,2,3,4,5,6,7,8,9])
            sage: M.LLL_FP()
            2
            sage: M
            [
            [0 0 0]
            [2 1 0]
            [-1 1 3]
            ]
            sage: M=ntl.mat_ZZ(4,4,[-6,9,-15,-18,4,-6,10,12,10,-16,18,35,-24,36,-46,-82]); M
            [
            [-6 9 -15 -18]
            [4 -6 10 12]
            [10 -16 18 35]
            [-24 36 -46 -82]
            ]
            sage: M.LLL_FP()
            3
            sage: M
            [
            [0 0 0 0]
            [0 -2 0 0]
            [-2 1 -5 -6]
            [0 -1 -7 5]
            ]

        WARNING: This method modifies \code{self}. So after applying this
        method your matrix will be a vector of vectors.
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_LLL_FP_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_LLL_FP(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank

    def LLL_QP(self, delta, return_U=False, verbose=False):
        r"""
        Perform the same reduction as \code{self.LLL_FP} using the
        same calling conventions but with quad float precision.

        EXAMPLES::

            sage: M=ntl.mat_ZZ(3,3,[1,2,3,4,5,6,7,8,9])
            sage: M.LLL_QP(delta=0.75)
            2
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_LLL_QP_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_LLL_QP(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank

    def LLL_XD(self, delta, return_U=False, verbose=False):
        r"""
        Perform the same reduction as \code{self.LLL_FP} using the
        same calling conventions but with extended exponent double
        precision.

        EXAMPLES::

            sage: M=ntl.mat_ZZ(3,3,[1,2,3,4,5,6,7,8,9])
            sage: M.LLL_XD(delta=0.75)
            2
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_LLL_XD_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_LLL_XD(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank

    def LLL_RR(self, delta, return_U=False, verbose=False):
        r"""
        Perform the same reduction as \code{self.LLL_FP} using the
        same calling conventions but with arbitrary precision floating
        point numbers.

        EXAMPLES::

            sage: M=ntl.mat_ZZ(3,3,[1,2,3,4,5,6,7,8,9])
            sage: M.LLL_RR(delta=0.75)
            2
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ.__new__(ntl_mat_ZZ)
            sig_on()
            rank = int(mat_ZZ_LLL_RR_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_LLL_RR(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank

    # Givens Orthogonalization.  This is a bit slower, but generally
    # much more stable, and is really the preferred orthogonalization
    # strategy.  For a nice description of this, see Chapter 5 of
    # [G. Golub and C. van Loan, Matrix Computations, 3rd edition,
    # Johns Hopkins Univ. Press, 1996].

    def G_LLL_FP(self, delta, return_U=False, verbose=False):
        r"""
        Perform the same reduction as self.LLL_FP using the same
        calling conventions but uses the Givens Orthogonalization.

        Givens Orthogonalization.  This is a bit slower, but generally
        much more stable, and is really the preferred
        orthogonalization strategy.  For a nice description of this,
        see Chapter 5 of [G. Golub and C. van Loan, Matrix
        Computations, 3rd edition, Johns Hopkins Univ. Press, 1996].
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_G_LLL_FP_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_G_LLL_FP(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank

    def G_LLL_QP(self, delta, return_U=False, verbose=False):
        r"""
        Perform the same reduction as self.G_LLL_FP using the same
        calling conventions but with quad float precision.
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_G_LLL_QP_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_G_LLL_QP(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank

    def G_LLL_XD(self, delta, return_U=False, verbose=False):
        r"""
        Perform the same reduction as self.G_LLL_FP using the same
        calling conventions but with extended exponent double
        precision.
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_G_LLL_XD_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_G_LLL_XD(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank

    def G_LLL_RR(self, delta, return_U=False, verbose=False):
        r"""
        Perform the same reduction as self.G_LLL_FP using the same
        calling conventions but with arbitrary precision floating
        point numbers.
        """
        cdef ntl_mat_ZZ U
        if return_U:
            U = ntl_mat_ZZ(self.__nrows, self.__nrows)
            sig_on()
            rank = int(mat_ZZ_G_LLL_RR_U(self.x, U.x, float(delta), 0, 0, int(verbose)))
            sig_off()
            return rank, U
        else:
            sig_on()
            rank = int(mat_ZZ_G_LLL_RR(self.x,float(delta),0,0,int(verbose)))
            sig_off()
            return rank
