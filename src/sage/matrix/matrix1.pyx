"""
Base class for matrices, part 1

For design documentation see :mod:`sage.matrix.docs`.

TESTS::

    sage: A = Matrix(GF(5),3,3,srange(9))
    sage: TestSuite(A).run()
"""

################################################################################
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL).
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
################################################################################

include "sage/ext/stdsage.pxi"
include "sage/ext/python.pxi"

import sage.modules.free_module

cdef class Matrix(matrix0.Matrix):
    ###################################################
    # Coercion to Various Systems
    ###################################################

    def _pari_init_(self):
        """
        Return a string defining a GP representation of self.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]
            sage: b = gp(a); b   # indirect doctest
            [x + 1, 2/3; 1/2*x^2, x^3 + 1]
            sage: a.determinant()
            x^4 + x^3 - 1/3*x^2 + x + 1
            sage: b.matdet()
            x^4 + x^3 - 1/3*x^2 + x + 1
        """
        w = self.list()
        cdef Py_ssize_t nr, nc, i, j
        nr = self._nrows
        nc = self._ncols
        v = []
        for i from 0 <= i < nr:
            tmp = []
            for j from 0 <= j < nc:
                tmp.append(w[i*nc + j]._pari_init_())
            v.append( ','.join(tmp))
        return 'Mat([%s])'%(';'.join(v))

    def _pari_(self):
        """
        Return the Pari matrix corresponding to self.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]
            sage: b = pari(a); b  # indirect doctest
            [x + 1, 2/3; 1/2*x^2, x^3 + 1]
            sage: a.determinant()
            x^4 + x^3 - 1/3*x^2 + x + 1
            sage: b.matdet()
            x^4 + x^3 - 1/3*x^2 + x + 1

        This function preserves precision for entries of inexact type (e.g.
        reals)::

            sage: R = RealField(4)       # 4 bits of precision
            sage: a = matrix(R, 2, [1, 2, 3, 1]); a
            [1.0 2.0]
            [3.0 1.0]
            sage: b = pari(a); b
            [1.000000000, 2.000000000; 3.000000000, 1.000000000] # 32-bit
            [1.00000000000000, 2.00000000000000; 3.00000000000000, 1.00000000000000] # 64-bit
            sage: b[0][0].precision()    # in words
            3
        """
        from sage.libs.pari.all import pari
        return pari.matrix(self._nrows, self._ncols, self._list())

    def _gap_init_(self):
        """
        Returns a string defining a gap representation of self.

        EXAMPLES::

            sage: A = MatrixSpace(QQ,3,3)([0,1,2,3,4,5,6,7,8])
            sage: g = gap(A) # indirect doctest
            sage: g
            [ [ 0, 1, 2 ], [ 3, 4, 5 ], [ 6, 7, 8 ] ]
            sage: g.CharacteristicPolynomial()
            x_1^3-12*x_1^2-18*x_1
            sage: A.characteristic_polynomial()
            x^3 - 12*x^2 - 18*x
            sage: matrix(g,QQ) == A
            True

        Particularly difficult is the case of matrices over
        cyclotomic fields and general number fields. See
        tickets #5618 and #8909.
        ::

            sage: K.<zeta> = CyclotomicField(8)
            sage: A = MatrixSpace(K,2,2)([0,1+zeta,2*zeta,3])
            sage: g = gap(A); g
            [ [ 0, 1+E(8) ], [ 2*E(8), 3 ] ]
            sage: matrix(g,K) == A
            True
            sage: g.IsMatrix()
            true

            sage: L.<tau> = NumberField(x^3-2)
            sage: A = MatrixSpace(L,2,2)([0,1+tau,2*tau,3])
            sage: g = gap(A); g
            [ [ !0, tau+1 ], [ 2*tau, !3 ] ]
            sage: matrix(g,L) == A
            True

        """
        cdef Py_ssize_t i, j
        v = []
        for i from 0 <= i < self._nrows:
            tmp = []
            for j from 0 <= j < self._ncols:
                tmp.append(self.get_unsafe(i,j)._gap_init_())
            v.append( '[%s]'%(','.join(tmp)) )
        # It is needed to multiply with 'One(...)', because
        # otherwise the result would not be a gap matrix
        return '[%s]*One(%s)'%(','.join(v),sage.interfaces.gap.gap(self.base_ring()).name())

    def _libgap_(self):
        """
        Construct a LibGAP matrix.

        INPUT:

        - ``M`` -- a matrix.

        OUTPUT:

        A GAP matrix, that is, a list of lists with entries over a
        common ring.

        EXAMPLES::

            sage: libgap(identity_matrix(ZZ,2))
            [ [ 1, 0 ], [ 0, 1 ] ]
            sage: libgap(matrix(GF(3),2,2,[4,5,6,7]))
            [ [ Z(3)^0, Z(3) ], [ 0*Z(3), Z(3)^0 ] ]
        """
        from sage.libs.gap.libgap import libgap
        return libgap._construct_matrix(self)

    def _giac_init_(self):
        """
        Return a Giac string representation of this matrix.

        EXAMPLES::

            sage: M = matrix(ZZ,2,range(4))             #optional
            sage: giac(M)                              #optional (indirect doctest)
            [[0,1],[2,3]]

        ::

            sage: M = matrix(QQ,3,[1,2,3,4/3,5/3,6/4,7,8,9])    #optional
            sage: giac(M)                                      #optional
            [[1,2,3],[4/3,5/3,3/2],[7,8,9]]

        ::

            sage: P.<x> = ZZ[]                          #optional
            sage: M = matrix(P, 2, [-9*x^2-2*x+2, x-1, x^2+8*x, -3*x^2+5]) #optional
            sage: giac(M)                             #optional
            [[-9*x^2-2*x+2,x-1],[x^2+8*x,-3*x^2+5]]
        """
        s = str(self.rows()).replace('(','[').replace(')',']')
        return "(%s)"%(s)

    def _maxima_init_(self):
        """
        Return a string representation of this matrix in Maxima.

        EXAMPLES::

            sage: m = matrix(3,range(9)); m
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: m._maxima_init_()
            'matrix([0,1,2],[3,4,5],[6,7,8])'
            sage: a = maxima(m); a
            matrix([0,1,2],[3,4,5],[6,7,8])
            sage: a.charpoly('x').expand()
            -x^3+12*x^2+18*x
            sage: m.charpoly()
            x^3 - 12*x^2 - 18*x
        """
        cdef Py_ssize_t i, j
        v = []
        for i from 0 <= i < self._nrows:
            tmp = []
            for j from 0 <= j < self._ncols:
                tmp.append(self.get_unsafe(i,j)._maxima_init_())
            v.append( '[%s]'%(','.join(tmp)) )
        return 'matrix(%s)'%(','.join(v))

    def _mathematica_init_(self):
        """
        Return Mathematica string representation of this matrix.

        EXAMPLES::

            sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
            sage: g = mathematica(A); g                  # optional
            {{1, 2, 3}, {4/3, 5/3, 3/2}, {7, 8, 9}}
            sage: A._mathematica_init_()
            '{{1/1, 2/1, 3/1}, {4/3, 5/3, 3/2}, {7/1, 8/1, 9/1}}'

        ::

            sage: A = matrix([[1,2],[3,4]])
            sage: g = mathematica(A); g                  # optional
            {{1, 2}, {3, 4}}

        ::

            sage: a = matrix([[pi, sin(x)], [cos(x), 1/e]]); a
            [    pi sin(x)]
            [cos(x) e^(-1)]
            sage: a._mathematica_init_()
            '{{Pi, Sin[x]}, {Cos[x], Exp[-1]}}'
        """
        return '{' + ', '.join([v._mathematica_init_() for v in self.rows()]) + '}'

    def _magma_init_(self, magma):
        r"""
        Return a string that evaluates in the given Magma session to this
        matrix.

        EXAMPLES:

        We first coerce a square matrix.

        ::

            sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
            sage: B = magma(A); B                       # (indirect doctest) optional - magma
            [  1   2   3]
            [4/3 5/3 3/2]
            [  7   8   9]
            sage: B.Type()                              # optional - magma
            AlgMatElt
            sage: B.Parent()                            # optional - magma
            Full Matrix Algebra of degree 3 over Rational Field

        We coerce a non-square matrix over
        `\ZZ/8\ZZ`.

        ::

            sage: A = MatrixSpace(Integers(8),2,3)([-1,2,3,4,4,-2])
            sage: B = magma(A); B                       # optional - magma
            [7 2 3]
            [4 4 6]
            sage: B.Type()                              # optional - magma
            ModMatRngElt
            sage: B.Parent()                            # optional - magma
            Full RMatrixSpace of 2 by 3 matrices over IntegerRing(8)

        ::

            sage: R.<x,y> = QQ[]
            sage: A = MatrixSpace(R,2,2)([x+y,x-1,y+5,x*y])
            sage: B = magma(A); B                       # optional - magma
            [x + y x - 1]
            [y + 5   x*y]

        ::

            sage: R.<x,y> = ZZ[]
            sage: A = MatrixSpace(R,2,2)([x+y,x-1,y+5,x*y])
            sage: B = magma(A); B                       # optional - magma
            [x + y x - 1]
            [y + 5   x*y]

        We coerce a matrix over a cyclotomic field, where the generator
        must be named during the coercion.

        ::

            sage: K = CyclotomicField(9) ; z = K.0
            sage: M = matrix(K,3,3,[0,1,3,z,z**4,z-1,z**17,1,0])
            sage: M
            [                 0                  1                  3]
            [             zeta9            zeta9^4          zeta9 - 1]
            [-zeta9^5 - zeta9^2                  1                  0]
            sage: magma(M)                             # optional - magma
            [                 0                  1                  3]
            [             zeta9            zeta9^4          zeta9 - 1]
            [-zeta9^5 - zeta9^2                  1                  0]
            sage: magma(M**2) == magma(M)**2           # optional - magma
            True
        """
        P = magma(self.parent())
        v = [x._magma_init_(magma) for x in self.list()]
        return '%s![%s]'%(P.name(), ','.join(v))

    def _maple_init_(self):
        """
        Return a Maple string representation of this matrix.

        EXAMPLES::

            sage: M = matrix(ZZ,2,range(4))             #optional
            sage: maple(M)                              #optional (indirect doctest)
            Matrix(2, 2, [[0,1],[2,3]])

        ::

            sage: M = matrix(QQ,3,[1,2,3,4/3,5/3,6/4,7,8,9])    #optional
            sage: maple(M)                                      #optional
            Matrix(3, 3, [[1,2,3],[4/3,5/3,3/2],[7,8,9]])

        ::

            sage: P.<x> = ZZ[]                          #optional
            sage: M = matrix(P, 2, [-9*x^2-2*x+2, x-1, x^2+8*x, -3*x^2+5]) #optional
            sage: maple(M)                             #optional
            Matrix(2, 2, [[-9*x^2-2*x+2,x-1],[x^2+8*x,-3*x^2+5]])
        """
        s = str(self.rows()).replace('(','[').replace(')',']')
        return "Matrix(%s,%s,%s)"%(self.nrows(), self.ncols(), s)

    def _singular_(self, singular=None):
        """
        Tries to coerce this matrix to a singular matrix.
        """
        if singular is None:
            from sage.interfaces.all import singular as singular_default
            singular = singular_default
        try:
            self.base_ring()._singular_(singular)
        except (NotImplementedError, AttributeError):
            raise TypeError, "Cannot coerce to Singular"

        return singular.matrix(self.nrows(),self.ncols(),singular(self.list()))

    def _macaulay2_(self, macaulay2=None):
        """
        EXAMPLES::

            sage: m = matrix(ZZ, [[1,2],[3,4]])
            sage: macaulay2(m)                  #optional (indirect doctest)
            | 1 2 |
            | 3 4 |

        ::

            sage: R.<x,y> = QQ[]
            sage: m = matrix([[x,y],[1+x,1+y]])
            sage: macaulay2(m)                  #optional
            | x   y   |
            | x+1 y+1 |
        """
        base_ring = macaulay2(self.base_ring())
        entries = map(list, self)
        return macaulay2(entries).matrix()


    def _scilab_init_(self):
        """
        Returns a string defining a Scilab representation of self.

        EXAMPLES:

            sage: a = matrix([[1,2,3],[4,5,6],[7,8,9]]); a  # optional - scilab
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: a._scilab_init_()         # optional - scilab
            '[1,2,3;4,5,6;7,8,9]'

        AUTHORS:

        - Ronan Paixao (2008-12-12)
        """
        w = self.list()
        cdef Py_ssize_t nr, nc, i, j
        nr = self._nrows
        nc = self._ncols
        v = []
        for i from 0 <= i < nr:
            tmp = []
            for j from 0 <= j < nc:
                tmp.append(w[i*nc + j]._pari_init_())
            v.append( ','.join(tmp))
        return '[%s]'%(';'.join(v))

    def _scilab_(self, scilab=None):
        """
        Creates a ScilabElement object based on self and returns it.

        EXAMPLES:

            sage: a = matrix([[1,2,3],[4,5,6],[7,8,9]]); a  # optional - scilab
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: b = scilab(a); b      # optional - scilab (indirect doctest)
                1.    2.    3.
                4.    5.    6.
                7.    8.    9.

        AUTHORS:

        - Ronan Paixao (2008-12-12)
        """
        return scilab(self._scilab_init_())


    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(matrix(QQ, 3, 3, [5..13])/7, verify=True)
            # Verified
            matrix(QQ, [[5/7, 6/7, 1], [8/7, 9/7, 10/7], [11/7, 12/7, 13/7]])
            sage: sage_input(MatrixSpace(GF(5), 50, 50, sparse=True).random_element(density=0.002), verify=True)
            # Verified
            matrix(GF(5), 50, 50, {(4,44):2, (5,25):1, (26,9):3, (43,24):3, (44,38):4})
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: matrix(RDF, [[3, 1], [4, 1]])._sage_input_(SageInputBuilder(), False)
            {call: {atomic:matrix}({atomic:RDF}, {list: ({list: ({atomic:3}, {atomic:1})}, {list: ({atomic:4}, {atomic:1})})})}
            sage: matrix(ZZ, 50, 50, {(9,17):1})._sage_input_(SageInputBuilder(), False)
            {call: {atomic:matrix}({atomic:ZZ}, {atomic:50}, {atomic:50}, {dict: {{atomic:(9,17)}:{atomic:1}}})}

        TESTS::

            sage: sage_input(matrix(RR, 0, 3, []), verify=True)
            # Verified
            matrix(RR, 0, 3)
            sage: sage_input(matrix(RR, 3, 0, []), verify=True)
            # Verified
            matrix(RR, 3, 0)
            sage: sage_input(matrix(RR, 0, 0, []), verify=True)
            # Verified
            matrix(RR, 0, 0)
        """
        if self.is_sparse():
            entries = list(self.dict().items())
            entries.sort()
            # We hand-format the keys to get rid of the space that would
            # normally follow the comma
            entries = [(sib.name('(%d,%d)'%k), sib(v, 2)) for k,v in entries]
            return sib.name('matrix')(self.base_ring(),
                                      sib.int(self.nrows()),
                                      sib.int(self.ncols()),
                                      sib.dict(entries))
        elif self.nrows() == 0 or self.ncols() == 0:
            return sib.name('matrix')(self.base_ring(),
                                      sib.int(self.nrows()),
                                      sib.int(self.ncols()))
        else:
            entries = [[sib(v, 2) for v in row] for row in self.rows()]
            return sib.name('matrix')(self.base_ring(), entries)


    def numpy(self, dtype=None):
        """
        Return the Numpy matrix associated to this matrix.

        INPUT:

        - ``dtype`` - The desired data-type for the array. If not given,
          then the type will be determined as the minimum type required
          to hold the objects in the sequence.

        EXAMPLES::

            sage: a = matrix(3,range(12))
            sage: a.numpy()
            array([[ 0,  1,  2,  3],
                   [ 4,  5,  6,  7],
                   [ 8,  9, 10, 11]])
            sage: a.numpy('f')
            array([[  0.,   1.,   2.,   3.],
                   [  4.,   5.,   6.,   7.],
                   [  8.,   9.,  10.,  11.]], dtype=float32)
            sage: a.numpy('d')
            array([[  0.,   1.,   2.,   3.],
                   [  4.,   5.,   6.,   7.],
                   [  8.,   9.,  10.,  11.]])
            sage: a.numpy('B')
            array([[ 0,  1,  2,  3],
                   [ 4,  5,  6,  7],
                   [ 8,  9, 10, 11]], dtype=uint8)

        Type ``numpy.typecodes`` for a list of the possible
        typecodes::

            sage: import numpy
            sage: sorted(numpy.typecodes.items())
            [('All', '?bhilqpBHILQPefdgFDGSUVOMm'), ('AllFloat', 'efdgFDG'), ('AllInteger', 'bBhHiIlLqQpP'), ('Character', 'c'), ('Complex', 'FDG'), ('Datetime', 'Mm'), ('Float', 'efdg'), ('Integer', 'bhilqp'), ('UnsignedInteger', 'BHILQP')]

        Alternatively, numpy automatically calls this function (via
        the magic :meth:`__array__` method) to convert Sage matrices
        to numpy arrays::

            sage: import numpy
            sage: b=numpy.array(a); b
            array([[ 0,  1,  2,  3],
                   [ 4,  5,  6,  7],
                   [ 8,  9, 10, 11]])
            sage: b.dtype
            dtype('int32')  # 32-bit
            dtype('int64')  # 64-bit
            sage: b.shape
            (3, 4)
        """
        import numpy
        A = numpy.matrix(self.list(), dtype=dtype)
        return numpy.resize(A,(self.nrows(), self.ncols()))

    # Define the magic "__array__" function so that numpy.array(m) can convert
    # a matrix m to a numpy array.
    # See http://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html#converting-an-arbitrary-sequence-object
    __array__=numpy

    ###################################################
    # Construction functions
    ###################################################

    def matrix_over_field(self):
        """
        Return copy of this matrix, but with entries viewed as elements of
        the fraction field of the base ring (assuming it is defined).

        EXAMPLES::

            sage: A = MatrixSpace(IntegerRing(),2)([1,2,3,4])
            sage: B = A.matrix_over_field()
            sage: B
            [1 2]
            [3 4]
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        """
        return self.change_ring(self.base_ring().fraction_field())


    def lift(self):
        """
        Return lift of self to the covering ring of the base ring R,
        which is by definition the ring returned by calling
        cover_ring() on R, or just R itself if the cover_ring method
        is not defined.

        EXAMPLES::

            sage: M = Matrix(Integers(7), 2, 2, [5, 9, 13, 15]) ; M
            [5 2]
            [6 1]
            sage: M.lift()
            [5 2]
            [6 1]
            sage: parent(M.lift())
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

        The field QQ doesn't have a cover_ring method::

            sage: hasattr(QQ, 'cover_ring')
            False

        So lifting a matrix over QQ gives back the same exact matrix.

        ::

            sage: B = matrix(QQ, 2, [1..4])
            sage: B.lift()
            [1 2]
            [3 4]
            sage: B.lift() is B
            True
        """
        if hasattr(self._base_ring, 'cover_ring'):
            S = self._base_ring.cover_ring()
            if S is not self._base_ring:
                return self.change_ring(S)
        return self

    #############################################################################################
    # rows, columns, sparse_rows, sparse_columns, dense_rows, dense_columns, row, column
    #############################################################################################

    def columns(self, copy=True):
        r"""
        Return a list of the columns of self.

        INPUT:

        - ``copy`` - (default: True) if True, return a copy of the list
          of columns which is safe to change.

        If ``self`` is a sparse matrix, columns are returned as sparse vectors,
        otherwise returned vectors are dense.

        EXAMPLES::

            sage: matrix(3, [1..9]).columns()
            [(1, 4, 7), (2, 5, 8), (3, 6, 9)]
            sage: matrix(RR, 2, [sqrt(2), pi, exp(1), 0]).columns()
            [(1.41421356237310, 2.71828182845905), (3.14159265358979, 0.000000000000000)]
            sage: matrix(RR, 0, 2, []).columns()
            [(), ()]
            sage: matrix(RR, 2, 0, []).columns()
            []
            sage: m = matrix(RR, 3, 3, {(1,2): pi, (2, 2): -1, (0,1): sqrt(2)})
            sage: parent(m.columns()[0])
            Sparse vector space of dimension 3 over Real Field with 53 bits of precision

        Sparse matrices produce sparse columns.  ::

            sage: A = matrix(QQ, 2, range(4), sparse=True)
            sage: v = A.columns()[0]
            sage: v.is_sparse()
            True

        TESTS::

            sage: A = matrix(QQ, 4, range(16))
            sage: A.columns('junk')
            Traceback (most recent call last):
            ...
            ValueError: 'copy' must be True or False, not junk
        """
        if not copy in [True, False]:
            msg = "'copy' must be True or False, not {0}"
            raise ValueError(msg.format(copy))
        x = self.fetch('columns')
        if not x is None:
            if copy: return list(x)
            return x
        if self.is_sparse():
            columns = self.sparse_columns(copy=copy)
        else:
            columns = self.dense_columns(copy=copy)
        self.cache('columns', columns)
        if copy: return list(columns)
        return columns

    def rows(self, copy=True):
        r"""
        Return a list of the rows of self.

        INPUT:

        - ``copy`` - (default: True) if True, return a copy of the list
          of rows which is safe to change.

        If ``self`` is a sparse matrix, rows are returned as sparse vectors,
        otherwise returned vectors are dense.

        EXAMPLES::

            sage: matrix(3, [1..9]).rows()
            [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
            sage: matrix(RR, 2, [sqrt(2), pi, exp(1), 0]).rows()
            [(1.41421356237310, 3.14159265358979), (2.71828182845905, 0.000000000000000)]
            sage: matrix(RR, 0, 2, []).rows()
            []
            sage: matrix(RR, 2, 0, []).rows()
            [(), ()]
            sage: m = matrix(RR, 3, 3, {(1,2): pi, (2, 2): -1, (0,1): sqrt(2)})
            sage: parent(m.rows()[0])
            Sparse vector space of dimension 3 over Real Field with 53 bits of precision

        Sparse matrices produce sparse rows.  ::

            sage: A = matrix(QQ, 2, range(4), sparse=True)
            sage: v = A.rows()[0]
            sage: v.is_sparse()
            True

        TESTS::

            sage: A = matrix(QQ, 4, range(16))
            sage: A.rows('junk')
            Traceback (most recent call last):
            ...
            ValueError: 'copy' must be True or False, not junk
        """
        if not copy in [True, False]:
            msg = "'copy' must be True or False, not {0}"
            raise ValueError(msg.format(copy))
        x = self.fetch('rows')
        if not x is None:
            if copy: return list(x)
            return x
        if self.is_sparse():
            rows = self.sparse_rows(copy=copy)
        else:
            rows = self.dense_rows(copy=copy)
        self.cache('rows', rows)
        if copy: return list(rows)
        return rows

    def dense_columns(self, copy=True):
        """
        Return list of the dense columns of self.

        INPUT:

        - ``copy`` - (default: True) if True, return a copy so you can
          modify it safely

        EXAMPLES:

        An example over the integers::

            sage: a = matrix(3,3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.dense_columns()
            [(0, 3, 6), (1, 4, 7), (2, 5, 8)]

        We do an example over a polynomial ring::

            sage: R.<x> = QQ[ ]
            sage: a = matrix(R, 2, [x,x^2, 2/3*x,1+x^5]); a
            [      x     x^2]
            [  2/3*x x^5 + 1]
            sage: a.dense_columns()
            [(x, 2/3*x), (x^2, x^5 + 1)]
            sage: a = matrix(R, 2, [x,x^2, 2/3*x,1+x^5], sparse=True)
            sage: c = a.dense_columns(); c
            [(x, 2/3*x), (x^2, x^5 + 1)]
            sage: parent(c[1])
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

        TESTS:

        Check that the returned rows are immutable as per :trac:`14874`::

            sage: m = Mat(ZZ,3,3)(range(9))
            sage: v = m.dense_columns()
            sage: map(lambda x: x.is_mutable(), v)
            [False, False, False]
        """
        x = self.fetch('dense_columns')
        if not x is None:
            if copy: return list(x)
            return x
        cdef Py_ssize_t i
        A = self if self.is_dense() else self.dense_matrix()
        C = [A.column(i) for i in range(self._ncols)]

        # Make the vectors immutable since we are caching them
        map(lambda x: x.set_immutable(), C)

        # cache result
        self.cache('dense_columns', C)
        if copy:
            return list(C)
        else:
            return C

    def dense_rows(self, copy=True):
        """
        Return list of the dense rows of self.

        INPUT:

        - ``copy`` - (default: True) if True, return a copy so you can
          modify it safely (note that the individual vectors in the copy
          should not be modified since they are mutable!)

        EXAMPLES::

            sage: m = matrix(3, range(9)); m
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: v = m.dense_rows(); v
            [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
            sage: v is m.dense_rows()
            False
            sage: m.dense_rows(copy=False) is m.dense_rows(copy=False)
            True
            sage: m[0,0] = 10
            sage: m.dense_rows()
            [(10, 1, 2), (3, 4, 5), (6, 7, 8)]

        TESTS:

        Check that the returned rows are immutable as per :trac:`14874`::

            sage: m = Mat(ZZ,3,3)(range(9))
            sage: v = m.dense_rows()
            sage: map(lambda x: x.is_mutable(), v)
            [False, False, False]
        """
        x = self.fetch('dense_rows')
        if not x is None:
            if copy: return list(x)
            return x

        cdef Py_ssize_t i
        A = self if self.is_dense() else self.dense_matrix()
        R = [A.row(i) for i in range(self._nrows)]

        # Make the vectors immutable since we are caching them
        map(lambda x: x.set_immutable(), R)

        # cache result
        self.cache('dense_rows', R)
        if copy:
            return list(R)
        else:
            return R


    def sparse_columns(self, copy=True):
        r"""
        Return a list of the columns of ``self`` as sparse vectors (or free module elements).

        INPUT:

        - ``copy`` - (default: True) if True, return a copy so you can
           modify it safely

        EXAMPLES::

            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: v = a.sparse_columns(); v
            [(0, 3), (1, 4), (2, 5)]
            sage: v[1].is_sparse()
            True

        TESTS:

        Columns of sparse matrices having no columns were fixed on :trac:`10714`::

            sage: m = matrix(10, 0, sparse=True)
            sage: m.ncols()
            0
            sage: m.columns()
            []

        Check that the returned columns are immutable as per :trac:`14874`::

            sage: m = Mat(ZZ,3,3,sparse=True)(range(9))
            sage: v = m.sparse_columns()
            sage: map(lambda x: x.is_mutable(), v)
            [False, False, False]
        """
        x = self.fetch('sparse_columns')
        if not x is None:
            if copy: return list(x)
            return x

        cdef Py_ssize_t i, j
        C = []
        if self._ncols > 0:
            F = sage.modules.free_module.FreeModule(self._base_ring, self._nrows, sparse=True)

            k = 0
            entries = {}
            for i, j in self.nonzero_positions(copy=False, column_order=True):
                if j > k:
                    # new column -- emit vector
                    while len(C) < k:
                        C.append(F(0))
                    C.append(F(entries, coerce=False, copy=False, check=False))
                    entries = {}
                    k = j
                entries[i] = self.get_unsafe(i, j)

            # finish up
            while len(C) < k:
                C.append(F(0))
            C.append(F(entries, coerce=False, copy=False, check=False))
            while len(C) < self._ncols:
                C.append(F(0))

        # Make the vectors immutable since we are caching them
        map(lambda x: x.set_immutable(), C)

        # cache and return result
        self.cache('sparse_columns', C)
        if copy:
            return list(C)
        else:
            return C

    def sparse_rows(self, copy=True):
        r"""
        Return a list of the rows of ``self`` as sparse vectors (or free module elements).

        INPUT:

        - ``copy`` - (default: True) if True, return a copy so you can
           modify it safely

        EXAMPLES::

            sage: m = Mat(ZZ,3,3,sparse=True)(range(9)); m
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: v = m.sparse_rows(); v
            [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
            sage: m.sparse_rows(copy=False) is m.sparse_rows(copy=False)
            True
            sage: v[1].is_sparse()
            True
            sage: m[0,0] = 10
            sage: m.sparse_rows()
            [(10, 1, 2), (3, 4, 5), (6, 7, 8)]

        TESTS:

        Rows of sparse matrices having no rows were fixed on :trac:`10714`::

            sage: m = matrix(0, 10, sparse=True)
            sage: m.nrows()
            0
            sage: m.rows()
            []

        Check that the returned rows are immutable as per :trac:`14874`::

            sage: m = Mat(ZZ,3,3,sparse=True)(range(9))
            sage: v = m.sparse_rows()
            sage: map(lambda x: x.is_mutable(), v)
            [False, False, False]
        """
        x = self.fetch('sparse_rows')
        if not x is None:
            if copy: return list(x)
            return x

        cdef Py_ssize_t i, j
        R = []
        if self._nrows > 0:
            F = sage.modules.free_module.FreeModule(self._base_ring, self._ncols, sparse=True)

            k = 0
            entries = {}
            for i, j in self.nonzero_positions(copy=False):
                if i > k:
                    # new row -- emit vector
                    while len(R) < k:
                        R.append(F(0))
                    R.append(F(entries, coerce=False, copy=False, check=False))
                    entries = {}
                    k = i
                entries[j] = self.get_unsafe(i, j)

            # finish up
            while len(R) < k:
                R.append(F(0))
            R.append(F(entries, coerce=False, copy=False, check=False))
            while len(R) < self._nrows:
                R.append(F(0))

        # Make the vectors immutable since we are caching them
        map(lambda x: x.set_immutable(), R)

        # cache and return result
        self.cache('sparse_rows', R)
        if copy:
            return list(R)
        else:
            return R

    def column(self, Py_ssize_t i, from_list=False):
        """
        Return the ``i``'th column of this matrix as a vector.

        This column is a dense vector if and only if the matrix is a dense
        matrix.

        INPUT:

        - ``i`` - integer

        - ``from_list`` - bool (default: False); if true, returns the
          ``i``'th element of ``self.columns()`` (see :func:`columns()`),
          which may be faster, but requires building a list of all
          columns the first time it is called after an entry of the
          matrix is changed.


        EXAMPLES::

            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.column(1)
            (1, 4)

        If the column is negative, it wraps around, just like with list
        indexing, e.g., -1 gives the right-most column::

            sage: a.column(-1)
            (2, 5)

        TESTS::

            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.column(3)
            Traceback (most recent call last):
            ...
            IndexError: column index out of range
            sage: a.column(-4)
            Traceback (most recent call last):
            ...
            IndexError: column index out of range
        """
        if self._ncols == 0:
            raise IndexError, "matrix has no columns"
        if i >= self._ncols or i < -self._ncols:
            raise IndexError, "column index out of range"
        i = i % self._ncols
        if i < 0:
            i = i + self._ncols
        if from_list:
            return self.columns(copy=False)[i]
        cdef Py_ssize_t j
        V = sage.modules.free_module.FreeModule(self._base_ring,
                                     self._nrows, sparse=self.is_sparse())
        tmp = [self.get_unsafe(j,i) for j in range(self._nrows)]
        return V(tmp, coerce=False, copy=False, check=False)

    def row(self, Py_ssize_t i, from_list=False):
        """
        Return the ``i``'th row of this matrix as a vector.

        This row is a dense vector if and only if the matrix is a dense
        matrix.

        INPUT:

        - ``i`` - integer

        - ``from_list`` - bool (default: False); if true, returns the
          ``i``'th element of ``self.rows()`` (see :func:`rows`), which
          may be faster, but requires building a list of all rows the
          first time it is called after an entry of the matrix is
          changed.

        EXAMPLES::

            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.row(0)
            (0, 1, 2)
            sage: a.row(1)
            (3, 4, 5)
            sage: a.row(-1)  # last row
            (3, 4, 5)

        TESTS::

            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.row(2)
            Traceback (most recent call last):
            ...
            IndexError: row index out of range
            sage: a.row(-3)
            Traceback (most recent call last):
            ...
            IndexError: row index out of range
        """
        if self._nrows == 0:
            raise IndexError, "matrix has no rows"
        if i >= self._nrows or i < -self._nrows:
            raise IndexError, "row index out of range"
        i = i % self._nrows
        if i < 0:
            i = i + self._nrows
        if from_list:
            return self.rows(copy=False)[i]
        cdef Py_ssize_t j
        V = sage.modules.free_module.FreeModule(self._base_ring,
                                      self._ncols, sparse=self.is_sparse())
        tmp = [self.get_unsafe(i,j) for j in range(self._ncols)]
        return V(tmp, coerce=False, copy=False, check=False)


    ###########################################################################
    # Building matrices out of other matrices, rows, or columns
    ###########################################################################
    def stack(self, bottom, subdivide=False):
        r"""
        Returns a new matrix formed by appending the matrix
        (or vector) ``bottom`` beneath ``self``.

        INPUT:

        - ``bottom`` - a matrix, vector or free module element, whose
          dimensions are compatible with ``self``.

        - ``subdivide`` - default: ``False`` - request the resulting
          matrix to have a new subdivision, separating ``self`` from ``bottom``.

        OUTPUT:

        A new matrix formed by appending ``bottom`` beneath ``self``.
        If ``bottom`` is a vector (or free module element) then in this context
        it is appropriate to consider it as a row vector.  (The code first
        converts a vector to a 1-row matrix.)

        If ``subdivide`` is ``True`` then any row subdivisions for
        the two matrices are preserved, and a new subdivision is added
        between ``self`` and ``bottom``.  If the column divisions are
        identical, then they are preserved, otherwise they are discarded.
        When ``subdivide`` is ``False`` there is no subdivision information
        in the result.

        .. warning::
            If ``subdivide`` is ``True`` then unequal column subdivisions
            will be discarded, since it would be ambiguous how to interpret
            them.  If the subdivision behavior is not what you need,
            you can manage subdivisions yourself with methods like
            :meth:`~sage.matrix.matrix2.Matrix.subdivisions`
            and
            :meth:`~sage.matrix.matrix2.Matrix.subdivide`.
            You might also find :func:`~sage.matrix.constructor.block_matrix`
            or
            :func:`~sage.matrix.constructor.block_diagonal_matrix`
            useful and simpler in some instances.


        EXAMPLES:

        Stacking with a matrix. ::

            sage: A = matrix(QQ, 4, 3, range(12))
            sage: B = matrix(QQ, 3, 3, range(9))
            sage: A.stack(B)
            [ 0  1  2]
            [ 3  4  5]
            [ 6  7  8]
            [ 9 10 11]
            [ 0  1  2]
            [ 3  4  5]
            [ 6  7  8]

        Stacking with a vector. ::

            sage: A = matrix(QQ, 3, 2, [0, 2, 4, 6, 8, 10])
            sage: v = vector(QQ, 2, [100, 200])
            sage: A.stack(v)
            [  0   2]
            [  4   6]
            [  8  10]
            [100 200]

        Errors are raised if the sizes are incompatible. ::

            sage: A = matrix(RR, [[1, 2],[3, 4]])
            sage: B = matrix(RR, [[10, 20, 30], [40, 50, 60]])
            sage: A.stack(B)
            Traceback (most recent call last):
            ...
            TypeError: number of columns must be the same, 2 != 3

            sage: v = vector(RR, [100, 200, 300])
            sage: A.stack(v)
            Traceback (most recent call last):
            ...
            TypeError: number of columns must be the same, 2 != 3

        Setting ``subdivide`` to ``True`` will, in its simplest form,
        add a subdivision between ``self`` and ``bottom``. ::

            sage: A = matrix(QQ, 2, 5, range(10))
            sage: B = matrix(QQ, 3, 5, range(15))
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [--------------]
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [10 11 12 13 14]

        Row subdivisions are preserved by stacking, and enriched,
        if subdivisions are requested.  (So multiple stackings can
        be recorded.) ::

            sage: A = matrix(QQ, 2, 4, range(8))
            sage: A.subdivide([1], None)
            sage: B = matrix(QQ, 3, 4, range(12))
            sage: B.subdivide([2], None)
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3]
            [-----------]
            [ 4  5  6  7]
            [-----------]
            [ 0  1  2  3]
            [ 4  5  6  7]
            [-----------]
            [ 8  9 10 11]

        Column subdivisions can be preserved, but only if they are identical.
        Otherwise, this information is discarded and must be managed
        separately. ::

            sage: A = matrix(QQ, 2, 5, range(10))
            sage: A.subdivide(None, [2,4])
            sage: B = matrix(QQ, 3, 5, range(15))
            sage: B.subdivide(None, [2,4])
            sage: A.stack(B, subdivide=True)
            [ 0  1| 2  3| 4]
            [ 5  6| 7  8| 9]
            [-----+-----+--]
            [ 0  1| 2  3| 4]
            [ 5  6| 7  8| 9]
            [10 11|12 13|14]

            sage: A.subdivide(None, [1,2])
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [--------------]
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [10 11 12 13 14]

        The result retains the base ring of ``self`` by coercing
        the elements of ``bottom`` into the base ring of ``self``. ::

            sage: A = matrix(QQ, 1, 2, [1,2])
            sage: B = matrix(RR, 1, 2, [sin(1.1), sin(2.2)])
            sage: C = A.stack(B); C
            [                  1                   2]
            [183017397/205358938 106580492/131825561]
            sage: C.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: D = B.stack(A); D
            [0.891207360061435 0.808496403819590]
            [ 1.00000000000000  2.00000000000000]
            sage: D.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Real Field with 53 bits of precision

        Sometimes it is not possible to coerce into the base ring
        of ``self``.  A solution is to change the base ring of ``self``
        to a more expansive ring.  Here we mix the rationals with
        a ring of polynomials with rational coefficients.  ::

            sage: R = PolynomialRing(QQ, 'y')
            sage: A = matrix(QQ, 1, 2, [1,2])
            sage: B = matrix(R, 1, 2, ['y', 'y^2'])

            sage: C = B.stack(A); C
            [  y y^2]
            [  1   2]
            sage: C.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field

            sage: D = A.stack(B)
            Traceback (most recent call last):
            ...
            TypeError: not a constant polynomial

            sage: E = A.change_ring(R)
            sage: F = E.stack(B); F
            [  1   2]
            [  y y^2]
            sage: F.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in y over Rational Field

        TESTS:

        A legacy test from the original implementation.  ::

            sage: M = Matrix(QQ, 2, 3, range(6))
            sage: N = Matrix(QQ, 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]

        AUTHOR:

        - Rob Beezer (2011-03-19) - rewritten to mirror code for :meth:`augment`
        """
        from sage.matrix.constructor import matrix

        if hasattr(bottom, '_vector_'):
            bottom = bottom.row()
        if not isinstance(bottom, sage.matrix.matrix1.Matrix):
            raise TypeError('a matrix must be stacked with another matrix, or '
                'a vector')

        cdef Matrix other
        other = bottom

        if self._ncols != other._ncols:
            raise TypeError('number of columns must be the same, '
                '{0} != {1}'.format(self._ncols, other._ncols))
        if not (self._base_ring is other.base_ring()):
            other = other.change_ring(self._base_ring)

        cdef Matrix Z
        Z = self.new_matrix(nrows = self._nrows + other._nrows)

        cdef Py_ssize_t r, c
        for r from 0 <= r < self._nrows:
            for c from 0 <= c < self._ncols:
                Z.set_unsafe(r,c, self.get_unsafe(r,c))
        nr = self.nrows()

        for r from 0 <= r < other._nrows:
            for c from 0 <= c < other._ncols:
                Z.set_unsafe(r+nr, c, other.get_unsafe(r,c))

        if subdivide:
            Z._subdivide_on_stack(self, other)

        return Z

    def augment(self, right, subdivide=False):
        r"""
        Returns a new matrix formed by appending the matrix (or vector)
        ``right`` on the right side of ``self``.

        INPUT:

        - ``right`` - a matrix, vector or free module element, whose
          dimensions are compatible with ``self``.

        - ``subdivide`` - default: ``False`` - request the resulting
          matrix to have a new subdivision, separating ``self`` from
          ``right``.

        OUTPUT:

        A new matrix formed by appending ``right`` onto the right side
        of ``self``.  If ``right`` is a vector (or free module element)
        then in this context it is appropriate to consider it as a
        column vector.  (The code first converts a vector to a 1-column
        matrix.)

        If ``subdivide`` is ``True`` then any column subdivisions for
        the two matrices are preserved, and a new subdivision is added
        between ``self`` and ``right``.  If the row divisions are
        identical, then they are preserved, otherwise they are
        discarded.  When ``subdivide`` is ``False`` there is no
        subdivision information in the result.

        .. warning::
            If ``subdivide`` is ``True`` then unequal row subdivisions
            will be discarded, since it would be ambiguous how to
            interpret them.  If the subdivision behavior is not what you
            need, you can manage subdivisions yourself with methods like
            :meth:`~sage.matrix.matrix2.Matrix.get_subdivisions` and
            :meth:`~sage.matrix.matrix2.Matrix.subdivide`.  You might
            also find :func:`~sage.matrix.constructor.block_matrix` or
            :func:`~sage.matrix.constructor.block_diagonal_matrix`
            useful and simpler in some instances.


        EXAMPLES:

        Augmenting with a matrix. ::

            sage: A = matrix(QQ, 3, range(12))
            sage: B = matrix(QQ, 3, range(9))
            sage: A.augment(B)
            [ 0  1  2  3  0  1  2]
            [ 4  5  6  7  3  4  5]
            [ 8  9 10 11  6  7  8]

        Augmenting with a vector. ::

            sage: A = matrix(QQ, 2, [0, 2, 4, 6, 8, 10])
            sage: v = vector(QQ, 2, [100, 200])
            sage: A.augment(v)
            [  0   2   4 100]
            [  6   8  10 200]

        Errors are raised if the sizes are incompatible. ::

            sage: A = matrix(RR, [[1, 2],[3, 4]])
            sage: B = matrix(RR, [[10, 20], [30, 40], [50, 60]])
            sage: A.augment(B)
            Traceback (most recent call last):
            ...
            TypeError: number of rows must be the same, 2 != 3

            sage: v = vector(RR, [100, 200, 300])
            sage: A.augment(v)
            Traceback (most recent call last):
            ...
            TypeError: number of rows must be the same, 2 != 3

        Setting ``subdivide`` to ``True`` will, in its simplest form,
        add a subdivision between ``self`` and ``right``. ::

            sage: A = matrix(QQ, 3, range(12))
            sage: B = matrix(QQ, 3, range(15))
            sage: A.augment(B, subdivide=True)
            [ 0  1  2  3| 0  1  2  3  4]
            [ 4  5  6  7| 5  6  7  8  9]
            [ 8  9 10 11|10 11 12 13 14]

        Column subdivisions are preserved by augmentation, and enriched,
        if subdivisions are requested.  (So multiple augmentations can
        be recorded.) ::

            sage: A = matrix(QQ, 3, range(6))
            sage: A.subdivide(None, [1])
            sage: B = matrix(QQ, 3, range(9))
            sage: B.subdivide(None, [2])
            sage: A.augment(B, subdivide=True)
            [0|1|0 1|2]
            [2|3|3 4|5]
            [4|5|6 7|8]

        Row subdivisions can be preserved, but only if they are
        identical.  Otherwise, this information is discarded and must be
        managed separately. ::

            sage: A = matrix(QQ, 3, range(6))
            sage: A.subdivide([1,3], None)
            sage: B = matrix(QQ, 3, range(9))
            sage: B.subdivide([1,3], None)
            sage: A.augment(B, subdivide=True)
            [0 1|0 1 2]
            [---+-----]
            [2 3|3 4 5]
            [4 5|6 7 8]
            [---+-----]

            sage: A.subdivide([1,2], None)
            sage: A.augment(B, subdivide=True)
            [0 1|0 1 2]
            [2 3|3 4 5]
            [4 5|6 7 8]

        The result retains the base ring of ``self`` by coercing the
        elements of ``right`` into the base ring of ``self``. ::

            sage: A = matrix(QQ, 2, [1,2])
            sage: B = matrix(RR, 2, [sin(1.1), sin(2.2)])
            sage: C = A.augment(B); C
            [                  1 183017397/205358938]
            [                  2 106580492/131825561]
            sage: C.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: D = B.augment(A); D
            [0.89120736006...  1.00000000000000]
            [0.80849640381...  2.00000000000000]
            sage: D.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Real Field with 53 bits of precision

        Sometimes it is not possible to coerce into the base ring of
        ``self``.  A solution is to change the base ring of ``self`` to
        a more expansive ring.  Here we mix the rationals with a ring of
        polynomials with rational coefficients.  ::

            sage: R = PolynomialRing(QQ, 'y')
            sage: A = matrix(QQ, 1, [1,2])
            sage: B = matrix(R, 1, ['y', 'y^2'])

            sage: C = B.augment(A); C
            [  y y^2   1   2]
            sage: C.parent()
            Full MatrixSpace of 1 by 4 dense matrices over Univariate Polynomial Ring in y over Rational Field

            sage: D = A.augment(B)
            Traceback (most recent call last):
            ...
            TypeError: not a constant polynomial

            sage: E = A.change_ring(R)
            sage: F = E.augment(B); F
            [  1   2   y y^2]
            sage: F.parent()
            Full MatrixSpace of 1 by 4 dense matrices over Univariate Polynomial Ring in y over Rational Field

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples
        - Rob Beezer (2010-12-07): vector argument, docstring, subdivisions
        """
        from sage.matrix.constructor import matrix

        if hasattr(right, '_vector_'):
            right = right.column()
        if not isinstance(right, sage.matrix.matrix1.Matrix):
            raise TypeError("a matrix must be augmented with another matrix, "
                "or a vector")

        cdef Matrix other
        other = right

        if self._nrows != other._nrows:
            raise TypeError('number of rows must be the same, '
                '{0} != {1}'.format(self._nrows, other._nrows))
        if not (self._base_ring is other.base_ring()):
            other = other.change_ring(self._base_ring)

        cdef Matrix Z
        Z = self.new_matrix(ncols = self._ncols + other._ncols)

        cdef Py_ssize_t r, c
        for r from 0 <= r < self._nrows:
            for c from 0 <= c < self._ncols:
                Z.set_unsafe(r,c, self.get_unsafe(r,c))
        nc = self.ncols()

        for r from 0 <= r < other._nrows:
            for c from 0 <= c < other._ncols:
                Z.set_unsafe(r, c+nc, other.get_unsafe(r,c))

        if subdivide:
            Z._subdivide_on_augment(self, other)

        return Z

    def matrix_from_columns(self, columns):
        """
        Return the matrix constructed from self using columns with indices
        in the columns list.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_columns([2,1])
            [2 1]
            [5 4]
            [0 7]
        """
        if not (PY_TYPE_CHECK(columns, list) or PY_TYPE_CHECK(columns, tuple)):
            raise TypeError, "columns (=%s) must be a list of integers"%columns
        cdef Matrix A
        cdef Py_ssize_t ncols,k,r

        ncols = PyList_GET_SIZE(columns)
        A = self.new_matrix(ncols = ncols)
        k = 0
        for i from 0 <= i < ncols:
            if columns[i] < 0 or columns[i] >= self._ncols:
                raise IndexError, "column %s out of range"%columns[i]
            for r from 0 <= r < self._nrows:
                A.set_unsafe(r,k, self.get_unsafe(r,columns[i]))
            k = k + 1
        return A

    def delete_columns(self, dcols, check=True):
        """
        Return the matrix constructed from deleting the columns with indices in the ``dcols`` list.

        INPUT:

        * ``dcols`` - list of indices of columns to be deleted from self.
        * ``check`` - checks whether any index in ``dcols`` is out of range. Defaults to ``True``.

        SEE ALSO:
            The methods :meth:`delete_rows` and :meth:`matrix_from_columns` are related.

        EXAMPLES::

            sage: A = Matrix(3,4,range(12)); A
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            sage: A.delete_columns([0,2])
            [ 1  3]
            [ 5  7]
            [ 9 11]

        ``dcols`` can be a tuple. But only the underlying set of indices matters. ::

            sage: A.delete_columns((2,0,2))
            [ 1  3]
            [ 5  7]
            [ 9 11]

        The default is to check whether any index in ``dcols`` is out of range. ::

            sage: A.delete_columns([-1,2,4])
            Traceback (most recent call last):
            ...
            IndexError: [4, -1] contains invalid indices.
            sage: A.delete_columns([-1,2,4], check=False)
            [ 0  1  3]
            [ 4  5  7]
            [ 8  9 11]

        TESTS:

        The list of indices is checked.  ::

            sage: A.delete_columns('junk')
            Traceback (most recent call last):
            ...
            TypeError: The argument must be a list or a tuple, not junk

        AUTHORS:
            - Wai Yan Pong (2012-03-05)
        """
        if not (PY_TYPE_CHECK(dcols, list) or PY_TYPE_CHECK(dcols, tuple)):
            raise TypeError("The argument must be a list or a tuple, not {l}".format(l=dcols))
        cdef list cols, diff_cols

        if check:
            diff_cols = list(set(dcols).difference(set(range(self._ncols))))
            if not (diff_cols == []):
                raise IndexError("{d} contains invalid indices.".format(d=diff_cols))
        cols = [k for k in range(self._ncols) if not k in dcols]
        return self.matrix_from_columns(cols)

    def matrix_from_rows(self, rows):
        """
        Return the matrix constructed from self using rows with indices in
        the rows list.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows([2,1])
            [6 7 0]
            [3 4 5]
        """
        if not (PY_TYPE_CHECK(rows, list) or PY_TYPE_CHECK(rows, tuple)):
            raise TypeError, "rows must be a list of integers"
        cdef Matrix A
        cdef Py_ssize_t nrows,k,c

        nrows = PyList_GET_SIZE(rows)
        A = self.new_matrix(nrows = nrows)
        k = 0
        for i from 0 <= i < nrows:
            if rows[i] < 0 or rows[i] >= self._nrows:
                raise IndexError, "row %s out of range"%rows[i]
            for c from 0 <= c < self._ncols:
                A.set_unsafe(k,c, self.get_unsafe(rows[i],c))
            k += 1
        return A


    def delete_rows(self, drows, check=True):
        """
        Return the matrix constructed from deleting the rows with indices in the ``drows`` list.

        INPUT:

        * ``drows`` - list of indices of rows to be deleted from self.
        * ``check`` - checks whether any index in ``drows`` is out of range. Defaults to ``True``.

        SEE ALSO:
            The methods :meth:`delete_columns` and :meth:`matrix_from_rows` are related.

        EXAMPLES::

            sage: A = Matrix(4,3,range(12)); A
            [ 0  1  2]
            [ 3  4  5]
            [ 6  7  8]
            [ 9 10 11]
            sage: A.delete_rows([0,2])
            [ 3  4  5]
            [ 9 10 11]

        ``drows`` can be a tuple. But only the underlying set of indices matters. ::

            sage: A.delete_rows((2,0,2))
            [ 3  4  5]
            [ 9 10 11]

        The default is to check whether the any index in ``drows`` is out of range. ::

            sage: A.delete_rows([-1,2,4])
            Traceback (most recent call last):
            ...
            IndexError: [4, -1] contains invalid indices.
            sage: A.delete_rows([-1,2,4], check=False)
            [ 0  1  2]
            [ 3  4  5]
            [ 9 10 11]

        TESTS:

        The list of indices is checked.  ::

            sage: A.delete_rows('junk')
            Traceback (most recent call last):
            ...
            TypeError: The argument must be a list or a tuple, not junk

        AUTHORS:
            - Wai Yan Pong (2012-03-05)
        """
        if not (PY_TYPE_CHECK(drows, list) or PY_TYPE_CHECK(drows, tuple)):
            raise TypeError("The argument must be a list or a tuple, not {l}".format(l=drows))
        cdef list rows, diff_rows

        if check:
            diff_rows = list(set(drows).difference(set(range(self._nrows))))
            if not (diff_rows == []):
                raise IndexError("{d} contains invalid indices.".format(d=diff_rows))
        rows = [k for k in range(self._nrows) if not k in drows]
        return self.matrix_from_rows(rows)

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from self from the given rows and
        columns.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows_and_columns([1], [0,2])
            [3 5]
            sage: A.matrix_from_rows_and_columns([1,2], [1,2])
            [4 5]
            [7 0]

        Note that row and column indices can be reordered or repeated::

            sage: A.matrix_from_rows_and_columns([2,1], [2,1])
            [0 7]
            [5 4]

        For example here we take from row 1 columns 2 then 0 twice, and do
        this 3 times.

        ::

            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0,0])
            [5 3 3]
            [5 3 3]
            [5 3 3]

        AUTHORS:

        - Jaap Spies (2006-02-18)

        - Didier Deshommes: some Pyrex speedups implemented
        """
        if not PY_TYPE_CHECK(rows, list):
            raise TypeError, "rows must be a list of integers"
        if not PY_TYPE_CHECK(columns, list):
            raise TypeError, "columns must be a list of integers"

        cdef Matrix A
        cdef Py_ssize_t nrows, ncols,k,r,i,j

        r = 0
        ncols = PyList_GET_SIZE(columns)
        nrows = PyList_GET_SIZE(rows)
        A = self.new_matrix(nrows = nrows, ncols = ncols)

        tmp = [el for el in columns if el >= 0 and el < self._ncols]
        columns = tmp
        if ncols != PyList_GET_SIZE(columns):
            raise IndexError, "column index out of range"

        for i from 0 <= i < nrows:
            if rows[i] < 0 or rows[i] >= self._nrows:
                raise IndexError, "row %s out of range"%i
            k = 0
            for j from 0 <= j < ncols:
                A.set_unsafe(r,k, self.get_unsafe(rows[i],columns[j]))
                k += 1
            r += 1
        return A

    def submatrix(self, Py_ssize_t row=0, Py_ssize_t col=0,
                        Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        """
        Return the matrix constructed from self using the specified
        range of rows and columns.

        INPUT:

        - ``row``, ``col`` - index of the starting row and column.
          Indices start at zero.

        - ``nrows``, ``ncols`` - (optional) number of rows and columns to
          take. If not provided, take all rows below and all columns to
          the right of the starting entry.

        SEE ALSO:

        The functions :func:`matrix_from_rows`,
        :func:`matrix_from_columns`, and
        :func:`matrix_from_rows_and_columns` allow one to select
        arbitrary subsets of rows and/or columns.

        EXAMPLES:

        Take the `3 \\times 3` submatrix starting from entry (1,1) in a
        `4 \\times 4` matrix::

            sage: m = matrix(4, [1..16])
            sage: m.submatrix(1, 1)
            [ 6  7  8]
            [10 11 12]
            [14 15 16]

        Same thing, except take only two rows::

            sage: m.submatrix(1, 1, 2)
            [ 6  7  8]
            [10 11 12]

        And now take only one column::

            sage: m.submatrix(1, 1, 2, 1)
            [ 6]
            [10]

        You can take zero rows or columns if you want::

            sage: m.submatrix(1, 1, 0)
            []
            sage: parent(m.submatrix(1, 1, 0))
            Full MatrixSpace of 0 by 3 dense matrices over Integer Ring
        """
        if nrows == -1:
            nrows = self._nrows - row
        if ncols == -1:
            ncols = self._ncols - col
        return self.matrix_from_rows_and_columns(range(row, row+nrows), range(col, col+ncols))



    def set_row(self, row, v):
        r"""
        Sets the entries of row ``row`` to the entries of ``v``.

        INPUT:

        - ``row`` - index of row to be set.

        - ``v`` - a list or vector of the new entries.

        OUTPUT:

        Changes the matrix in-place, so there is no output.

        EXAMPLES:

        New entries may be contained in a vector.::

            sage: A = matrix(QQ, 5, range(25))
            sage: u = vector(QQ, [0, -1, -2, -3, -4])
            sage: A.set_row(2, u)
            sage: A
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [ 0 -1 -2 -3 -4]
            [15 16 17 18 19]
            [20 21 22 23 24]

        New entries may be in any sort of list.::

            sage: A = matrix([[1, 2], [3, 4]]); A
            [1 2]
            [3 4]
            sage: A.set_row(0, [0, 0]); A
            [0 0]
            [3 4]
            sage: A.set_row(1, (0, 0)); A
            [0 0]
            [0 0]

        TESTS::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: A.set_row(2, [0, 0]); A
            Traceback (most recent call last):
            ...
            ValueError: row number must be between 0 and 1 (inclusive), not 2

            sage: A.set_row(0, [0, 0, 0])
            Traceback (most recent call last):
            ...
            ValueError: list of new entries must be of length 2 (not 3)

            sage: A = matrix(2, [1, 2, 3, 4])
            sage: A.set_row(0, [1/3, 1]); A
            Traceback (most recent call last):
            ...
            TypeError: Cannot set row with Rational Field elements over Integer Ring, use change_ring first.
        """
        if len(v) != self._ncols:
            msg = "list of new entries must be of length {0} (not {1})"
            raise ValueError(msg.format(self._ncols, len(v)))
        if (row < 0) or (row >= self._nrows):
            msg = "row number must be between 0 and {0} (inclusive), not {1}"
            raise ValueError(msg.format(self._nrows-1, row))

        try:
            for j in range(self._ncols):
                self[row, j] = v[j]
        except TypeError:
            msg = "Cannot set row with {0} elements over {1}, use change_ring first."
            raise TypeError(msg.format(v[j].parent(), self.base_ring()))

    def set_column(self, col, v):
        r"""
        Sets the entries of column ``col`` to the entries of ``v``.

        INPUT:

        - ``col`` - index of column to be set.

        - ``v`` - a list or vector of the new entries.

        OUTPUT:

        Changes the matrix in-place, so there is no output.

        EXAMPLES:

        New entries may be contained in a vector.::

            sage: A = matrix(QQ, 5, range(25))
            sage: u = vector(QQ, [0, -1, -2, -3, -4])
            sage: A.set_column(2, u)
            sage: A
            [ 0  1  0  3  4]
            [ 5  6 -1  8  9]
            [10 11 -2 13 14]
            [15 16 -3 18 19]
            [20 21 -4 23 24]

        New entries may be in any sort of list.::

            sage: A = matrix([[1, 2], [3, 4]]); A
            [1 2]
            [3 4]
            sage: A.set_column(0, [0, 0]); A
            [0 2]
            [0 4]
            sage: A.set_column(1, (0, 0)); A
            [0 0]
            [0 0]

        TESTS::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: A.set_column(2, [0, 0]); A
            Traceback (most recent call last):
            ...
            ValueError: column number must be between 0 and 1 (inclusive), not 2

            sage: A.set_column(0, [0, 0, 0])
            Traceback (most recent call last):
            ...
            ValueError: list of new entries must be of length 2 (not 3)

            sage: A = matrix(2, [1, 2, 3, 4])
            sage: A.set_column(0, [1/4, 1]); A
            Traceback (most recent call last):
            ...
            TypeError: Cannot set column with Rational Field elements over Integer Ring, use change_ring first.
        """
        if len(v) != self._nrows:
            msg = "list of new entries must be of length {0} (not {1})"
            raise ValueError(msg.format(self._nrows, len(v)))
        if (col < 0) or (col >= self._ncols):
            msg = "column number must be between 0 and {0} (inclusive), not {1}"
            raise ValueError(msg.format(self._ncols-1, col))

        try:
            for i in range(self._nrows):
                self[i, col] = v[i]
        except TypeError:
            msg = "Cannot set column with {0} elements over {1}, use change_ring first."
            raise TypeError(msg.format(v[i].parent(), self.base_ring()))


    ####################################################################################
    # Change of representation between dense and sparse.
    ####################################################################################

    def dense_matrix(self):
        """
        If this matrix is sparse, return a dense matrix with the same
        entries. If this matrix is dense, return this matrix (not a copy).

        .. note::

           The definition of "dense" and "sparse" in Sage have nothing to
           do with the number of nonzero entries. Sparse and dense are
           properties of the underlying representation of the matrix.

        EXAMPLES::

            sage: A = MatrixSpace(QQ,2, sparse=True)([1,2,0,1])
            sage: A.is_sparse()
            True
            sage: B = A.dense_matrix()
            sage: B.is_sparse()
            False
            sage: A*B
            [1 4]
            [0 1]
            sage: A.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        In Sage, the product of a sparse and a dense matrix is always
        dense::

            sage: (A*B).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: (B*A).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        TESTS:

        Make sure that subdivisions are preserved when switching
        between dense and sparse matrices::

            sage: a = matrix(ZZ, 3, range(9))
            sage: a.subdivide([1,2],2)
            sage: a.subdivisions()
            ([1, 2], [2])
            sage: b = a.sparse_matrix().dense_matrix()
            sage: b.subdivisions()
            ([1, 2], [2])
        """
        if self.is_dense():
            return self
        cdef Matrix A
        A = self.new_matrix(self._nrows, self._ncols, 0, coerce=False,
                               copy = False, sparse=False)
        for i,j in self.nonzero_positions():
            A.set_unsafe(i,j,self.get_unsafe(i,j))
        A.subdivide(self.subdivisions())
        return A

    def sparse_matrix(self):
        """
        If this matrix is dense, return a sparse matrix with the same
        entries. If this matrix is sparse, return this matrix (not a
        copy).

        .. note::

           The definition of "dense" and "sparse" in Sage have nothing
           to do with the number of nonzero entries. Sparse and dense are
           properties of the underlying representation of the matrix.

        EXAMPLES::

            sage: A = MatrixSpace(QQ,2, sparse=False)([1,2,0,1])
            sage: A.is_sparse()
            False
            sage: B = A.sparse_matrix()
            sage: B.is_sparse()
            True
            sage: A
            [1 2]
            [0 1]
            sage: B
            [1 2]
            [0 1]
            sage: A*B
            [1 4]
            [0 1]
            sage: A.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: B.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: (A*B).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: (B*A).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        """
        if self.is_sparse():
            return self
        A = self.new_matrix(self._nrows, self._ncols, entries = self.dict(), coerce=False,
                               copy = False, sparse=True)
        A.subdivide(self.subdivisions())
        return A

    def matrix_space(self, nrows=None, ncols=None, sparse=None):
        """
        Return the ambient matrix space of self.

        INPUT:

        - ``nrows``, ``ncols`` - (optional) number of rows and columns in
          returned matrix space.
        - ``sparse`` - whether the returned matrix space uses sparse or
          dense matrices.

        EXAMPLES::

            sage: m = matrix(3, [1..9])
            sage: m.matrix_space()
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
            sage: m.matrix_space(ncols=2)
            Full MatrixSpace of 3 by 2 dense matrices over Integer Ring
            sage: m.matrix_space(1)
            Full MatrixSpace of 1 by 3 dense matrices over Integer Ring
            sage: m.matrix_space(1, 2, True)
            Full MatrixSpace of 1 by 2 sparse matrices over Integer Ring
        """
        from sage.matrix.matrix_space import MatrixSpace
        if nrows is None:
            nrows = self._nrows
        if ncols is None:
            ncols = self._ncols
        if sparse is None:
            sparse = self.is_sparse()
        base_ring = self._base_ring
        return MatrixSpace(base_ring, nrows, ncols, sparse)

    def new_matrix(self, nrows=None, ncols=None, entries=None,
                   coerce=True, copy=True, sparse=None):
        """
        Create a matrix in the parent of this matrix with the given number
        of rows, columns, etc. The default parameters are the same as for
        self.

        INPUT:

        These three variables get sent to :func:`matrix_space`:

        - ``nrows``, ``ncols`` - number of rows and columns in returned
          matrix. If not specified, defaults to ``None`` and will give a
          matrix of the same size as self.
        - ``sparse`` - whether returned matrix is sparse or not. Defaults
          to same value as self.

        The remaining three variables (``coerce``, ``entries``, and
        ``copy``) are used by
        :func:`sage.matrix.matrix_space.MatrixSpace` to construct the
        new matrix.

        .. warning::

           This function called with no arguments returns the zero
           matrix of the same dimension and sparseness of self.

        EXAMPLES::

            sage: A = matrix(ZZ,2,2,[1,2,3,4]); A
            [1 2]
            [3 4]
            sage: A.new_matrix()
            [0 0]
            [0 0]
            sage: A.new_matrix(1,1)
            [0]
            sage: A.new_matrix(3,3).parent()
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

        ::

            sage: A = matrix(RR,2,3,[1.1,2.2,3.3,4.4,5.5,6.6]); A
            [1.10000000000000 2.20000000000000 3.30000000000000]
            [4.40000000000000 5.50000000000000 6.60000000000000]
            sage: A.new_matrix()
            [0.000000000000000 0.000000000000000 0.000000000000000]
            [0.000000000000000 0.000000000000000 0.000000000000000]
            sage: A.new_matrix().parent()
            Full MatrixSpace of 2 by 3 dense matrices over Real Field with 53 bits of precision

        """
        if self._nrows == nrows and self._ncols == ncols and (sparse is None or self.is_sparse() == sparse):
            return self._parent(entries=entries, coerce=coerce, copy=copy)
        return self.matrix_space(nrows, ncols, sparse=sparse)(entries=entries,
                                             coerce=coerce, copy=copy)
    def block_sum(self, Matrix other):
        """
        Return the block matrix that has self and other on the diagonal::

            [ self     0 ]
            [    0 other ]

        EXAMPLES::

            sage: A = matrix(QQ[['t']], 2, range(1, 5))
            sage: A.block_sum(100*A)
            [  1   2   0   0]
            [  3   4   0   0]
            [  0   0 100 200]
            [  0   0 300 400]
        """
        if not isinstance(other, Matrix):
            raise TypeError, "other must be a Matrix"
        top = self.augment(self.new_matrix(ncols=other._ncols))
        bottom = other.new_matrix(ncols=self._ncols).augment(other)
        return top.stack(bottom)

