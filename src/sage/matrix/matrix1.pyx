"""
Base class for matrices, part 1

For design documentation see :mod:`sage.matrix.docs`.

TESTS::

    sage: A = Matrix(GF(5), 3, 3, srange(9))
    sage: TestSuite(A).run()
"""

#*****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.sequence cimport PySequence_Fast

import sage.modules.free_module
from sage.structure.coerce cimport coercion_model


_MISSING = object()


cdef class Matrix(Matrix0):
    ###################################################
    # Coercion to Various Systems
    ###################################################

    def _pari_init_(self):
        """
        Return a string defining a GP representation of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]
            sage: b = gp(a); b   # indirect doctest                                     # needs sage.libs.pari
            [x + 1, 2/3; 1/2*x^2, x^3 + 1]
            sage: a.determinant()
            x^4 + x^3 - 1/3*x^2 + x + 1
            sage: b.matdet()                                                            # needs sage.libs.pari
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

    def __pari__(self):
        """
        Return the Pari matrix corresponding to ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]
            sage: b = pari(a); b  # indirect doctest                                    # needs sage.libs.pari
            [x + 1, 2/3; 1/2*x^2, x^3 + 1]
            sage: a.determinant()
            x^4 + x^3 - 1/3*x^2 + x + 1
            sage: b.matdet()                                                            # needs sage.libs.pari
            x^4 + x^3 - 1/3*x^2 + x + 1

        This function preserves precision for entries of inexact type (e.g.
        reals)::

            sage: R = RealField(4)       # 4 bits of precision
            sage: a = matrix(R, 2, [1, 2, 3, 1]); a
            [1.0 2.0]
            [3.0 1.0]
            sage: b = pari(a); b                                                        # needs sage.libs.pari
            [1.000000000, 2.000000000; 3.000000000, 1.000000000] # 32-bit
            [1.00000000000000, 2.00000000000000; 3.00000000000000, 1.00000000000000] # 64-bit
        """
        from sage.libs.pari import pari
        return pari.matrix(self._nrows, self._ncols, self._list())

    def _gap_init_(self):
        """
        Return a string defining a gap representation of ``self``.

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: A = MatrixSpace(QQ,3,3)([0,1,2,3,4,5,6,7,8])
            sage: g = gap(A)  # indirect doctest
            sage: g
            [ [ 0, 1, 2 ], [ 3, 4, 5 ], [ 6, 7, 8 ] ]
            sage: g.CharacteristicPolynomial()
            x_1^3-12*x_1^2-18*x_1
            sage: A.characteristic_polynomial()                                         # needs sage.libs.pari
            x^3 - 12*x^2 - 18*x
            sage: matrix(QQ, g) == A
            True

        Particularly difficult is the case of matrices over cyclotomic
        fields and general number fields. See :issue:`5618` and :issue:`8909`::

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: K.<zeta> = CyclotomicField(8)
            sage: A = MatrixSpace(K, 2, 2)([0, 1+zeta, 2*zeta, 3])
            sage: g = gap(A); g
            [ [ 0, 1+E(8) ], [ 2*E(8), 3 ] ]
            sage: matrix(K, g) == A
            True
            sage: g.IsMatrix()
            true

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: L.<tau> = NumberField(x^3 - 2)
            sage: A = MatrixSpace(L, 2, 2)([0, 1+tau, 2*tau, 3])
            sage: g = gap(A); g
            [ [ !0, tau+1 ], [ 2*tau, !3 ] ]
            sage: matrix(L, g) == A
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

        - ``M`` -- a matrix

        OUTPUT:

        A GAP matrix, that is, a list of lists with entries over a
        common ring.

        EXAMPLES::

            sage: libgap(identity_matrix(ZZ, 2))                                        # needs sage.libs.gap
            [ [ 1, 0 ], [ 0, 1 ] ]
            sage: libgap(matrix(GF(3), 2, 2, [4,5,6,7]))                                # needs sage.libs.gap
            [ [ Z(3)^0, Z(3) ], [ 0*Z(3), Z(3)^0 ] ]
        """
        from sage.libs.gap.libgap import libgap
        return libgap._construct_matrix(self)

    def _fricas_init_(self):
        """
        Return a FriCAS string representation of this matrix.

        EXAMPLES::

            sage: M = matrix(ZZ,2,range(4))
            sage: fricas(M)                                                     # optional - fricas
            +0  1+
            |    |
            +2  3+

        ::

            sage: M = matrix(QQ,3,[1,2,3,4/3,5/3,6/4,7,8,9])
            sage: fricas(M).sage().parent()                                     # optional - fricas
            Full MatrixSpace of 3 by 3 dense matrices over Rational Field

        ::

            sage: P.<x> = ZZ[]
            sage: M = matrix(P, 2, [-9*x^2-2*x+2, x-1, x^2+8*x, -3*x^2+5])
            sage: fricas(M)                                                     # optional - fricas
            +     2                      +
            |- 9 x  - 2 x + 2    x - 1   |
            |                            |
            |     2                 2    |
            +    x  + 8 x      - 3 x  + 5+

        ::

            sage: y = var('y')                                                          # needs sage.symbolic
            sage: M = matrix(SR, 2, [y+sin(y), y - 4, 1/y, dilog(y)])                   # needs sage.symbolic
            sage: M == fricas(M).sage()         # optional - fricas                     # needs sage.symbolic
            True
        """
        s = ','.join('[' + ','.join(cf._fricas_init_() for cf in row) + ']'
                     for row in self.rows())
        R = self.base_ring()
        try:
            R._fricas_()
        except TypeError:
            return "matrix([%s])" % s
        return "matrix([%s])$Matrix(%s)" % (s, R._fricas_init_())

    def _giac_init_(self):
        """
        Return a Giac string representation of this matrix.

        EXAMPLES::

            sage: M = matrix(ZZ, 2, range(4))
            sage: giac(M)                                                               # needs giac
            [[0,1],[2,3]]

            sage: M = matrix(QQ, 3, [1,2,3, 4/3,5/3,6/4, 7,8,9])
            sage: giac(M)                                                               # needs giac
            [[1,2,3],[4/3,5/3,3/2],[7,8,9]]

            sage: P.<x> = ZZ[]
            sage: M = matrix(P, 2, [-9*x^2-2*x+2, x-1, x^2+8*x, -3*x^2+5])
            sage: giac(M)                                                               # needs giac
            [[-9*sageVARx^2-2*sageVARx+2,sageVARx-1],[sageVARx^2+8*sageVARx,-3*sageVARx^2+5]]

            sage: y = var('y')                                                          # needs sage.symbolic
            sage: M = matrix(SR, 2, [y+sin(y), y - 4, 1/y, dilog(y)])                   # needs sage.symbolic
            sage: giac(M).det().sage()                                                  # needs giac sage.symbolic
            (y^2*dilog(y) + y*dilog(y)*sin(y) - y + 4)/y
        """
        s = ','.join('[' + ','.join(cf._giac_init_() for cf in row) + ']'
                     for row in self.rows())
        return "([%s])" % s

    def _maxima_init_(self):
        """
        Return a string representation of this matrix in Maxima.

        EXAMPLES::

            sage: m = matrix(3,range(9)); m
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: m._maxima_init_()                                                     # needs sage.symbolic
            'matrix([0,1,2],[3,4,5],[6,7,8])'
            sage: a = maxima(m); a                                                      # needs sage.symbolic
            matrix([0,1,2],[3,4,5],[6,7,8])
            sage: a.charpoly('x').expand()                                              # needs sage.symbolic
            ...-x^3...+12*x^2+18*x
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
            sage: g = mathematica(A); g                  # optional - mathematica
            {{1, 2, 3}, {4/3, 5/3, 3/2}, {7, 8, 9}}
            sage: A._mathematica_init_()                                                # needs sage.symbolic
            '{{1/1, 2/1, 3/1}, {4/3, 5/3, 3/2}, {7/1, 8/1, 9/1}}'

        ::

            sage: A = matrix([[1,2],[3,4]])
            sage: g = mathematica(A); g                  # optional - mathematica
            {{1, 2}, {3, 4}}

        ::

            sage: a = matrix([[pi, sin(x)], [cos(x), 1/e]]); a                          # needs sage.symbolic
            [    pi sin(x)]
            [cos(x) e^(-1)]
            sage: a._mathematica_init_()                                                # needs sage.symbolic
            '{{Pi, Sin[x]}, {Cos[x], Exp[-1]}}'
        """
        return '{' + ', '.join(v._mathematica_init_() for v in self.rows()) + '}'

    def _magma_init_(self, magma) -> str:
        r"""
        Return a string that evaluates in the given Magma session to this
        matrix.

        EXAMPLES:

        We first coerce a square matrix. ::

            sage: # optional - magma
            sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
            sage: B = magma(A); B  # indirect doctest
            [  1   2   3]
            [4/3 5/3 3/2]
            [  7   8   9]
            sage: B.Type()
            AlgMatElt
            sage: B.Parent()
            Full Matrix Algebra of degree 3 over Rational Field

        We coerce a non-square matrix over
        `\ZZ/8\ZZ`. ::

            sage: # optional - magma
            sage: A = MatrixSpace(Integers(8),2,3)([-1,2,3,4,4,-2])
            sage: B = magma(A); B
            [7 2 3]
            [4 4 6]
            sage: B.Type()
            ModMatRngElt
            sage: B.Parent()
            Full RMatrixSpace of 2 by 3 matrices over IntegerRing(8)

            sage: R.<x,y> = QQ[]
            sage: A = MatrixSpace(R,2,2)([x+y,x-1,y+5,x*y])
            sage: B = magma(A); B                       # optional - magma
            [x + y x - 1]
            [y + 5   x*y]

            sage: R.<x,y> = ZZ[]
            sage: A = MatrixSpace(R,2,2)([x+y,x-1,y+5,x*y])
            sage: B = magma(A); B                       # optional - magma
            [x + y x - 1]
            [y + 5   x*y]

        We coerce a matrix over a cyclotomic field, where the generator
        must be named during the coercion. ::

            sage: # optional - magma, needs sage.rings.number_field
            sage: K = CyclotomicField(9); z = K.0
            sage: M = matrix(K, 3, 3, [0,1,3,z,z**4,z-1,z**17,1,0]); M
            [                 0                  1                  3]
            [             zeta9            zeta9^4          zeta9 - 1]
            [-zeta9^5 - zeta9^2                  1                  0]
            sage: magma(M)
            [                 0                  1                  3]
            [             zeta9            zeta9^4          zeta9 - 1]
            [-zeta9^5 - zeta9^2                  1                  0]
            sage: magma(M**2) == magma(M)**2
            True

        One sparse matrix::

            sage: M = matrix(QQ,2,2,[4,6,55,0],sparse=True)
            sage: T = magma(M); T    # optional - magma
            Sparse matrix with 2 rows and 2 columns over Rational Field
            sage: T.Determinant()    # optional - magma
            -330
        """
        if self.is_sparse():
            R = magma(self.base_ring())
            s = "SparseMatrix({}, {}, {}, ["
            s = s.format(R.name(), self.nrows(), self.ncols())
            entries = ("<{}, {}, {}>".format(ij[0] + 1, ij[1] + 1,
                                             mij._magma_init_(magma))
                       for ij, mij in self.dict().items())
            return s + ', '.join(entries) + "])"
        else:
            P = magma(self.parent())
            v = [x._magma_init_(magma) for x in self.list()]
            return '%s![%s]' % (P.name(), ','.join(v))

    def _maple_init_(self):
        """
        Return a Maple string representation of this matrix.

        EXAMPLES::

            sage: M = matrix(ZZ,2,range(4))
            sage: maple(M)  # optional - maple
            Matrix(2, 2, [[0,1],[2,3]])

            sage: M = matrix(QQ,3,[1,2,3,4/3,5/3,6/4,7,8,9])
            sage: maple(M)  # optional - maple
            Matrix(3, 3, [[1,2,3],[4/3,5/3,3/2],[7,8,9]])

            sage: P.<x> = ZZ[]
            sage: M = matrix(P, 2, [-9*x^2-2*x+2, x-1, x^2+8*x, -3*x^2+5])
            sage: maple(M)  # optional - maple
            Matrix(2, 2, [[-9*x^2-2*x+2,x-1],[x^2+8*x,-3*x^2+5]])

            sage: y = var('y')                                                          # needs sage.symbolic
            sage: M = matrix(SR, 2, [y+sin(y), y - 4, 1/y, dilog(y)])                   # needs sage.symbolic
            sage: M == maple(M).sage()          # optional - maple                      # needs sage.symbolic
            True
        """
        s = ','.join('[' + ','.join(cf._maple_init_() for cf in row) + ']'
                     for row in self.rows())
        return "Matrix(%s,%s,[%s])" % (self.nrows(), self.ncols(), s)

    def _polymake_(self, polymake=None):
        """
        Try to coerce this matrix to a polymake matrix.

        EXAMPLES::

            sage: M = matrix(ZZ,2,range(4))
            sage: polymake(M)                           # optional - jupymake
            0 1
            2 3
            sage: K.<sqrt5> = QuadraticField(5)                                         # needs sage.rings.number_field
            sage: M = matrix(K, [[1, 2], [sqrt5, 3]])                                   # needs sage.rings.number_field
            sage: polymake(M)                           # optional - jupymake           # needs sage.rings.number_field
            1 2
            0+1r5 3
        """
        P = polymake(self.parent())
        return polymake.new_object(P, [ list(r) for r in self.rows(copy=False) ])

    def _singular_(self, singular=None):
        """
        Try to coerce this matrix to a singular matrix.
        """
        if singular is None:
            from sage.interfaces.singular import singular as singular_default
            singular = singular_default
        try:
            self.base_ring()._singular_(singular)
        except (NotImplementedError, AttributeError):
            raise TypeError("Cannot coerce to Singular")

        return singular.matrix(self.nrows(),self.ncols(),singular(self.list()))

    def _macaulay2_(self, macaulay2=None):
        """
        EXAMPLES::

            sage: m = matrix(ZZ, [[1,2],[3,4]])
            sage: macaulay2(m)  # indirect doctest              # optional - macaulay2
            | 1 2 |
            | 3 4 |

        ::

            sage: R.<x,y> = QQ[]
            sage: m = matrix([[x,y],[1+x,1+y]])
            sage: macaulay2(m)                                  # optional - macaulay2
            | x   y   |
            | x+1 y+1 |

        TESTS:

        Entries of the matrix get promoted to the base ring (:issue:`28566`)::

            sage: R.<x,y> = QQ[]
            sage: m = macaulay2(matrix(R, [[1, 2], [3, 4]]))    # optional - macaulay2
            sage: m.ring()._operator('===', R).sage()           # optional - macaulay2
            True

        Check that degenerate matrix dimensions are handled correctly
        (:issue:`28591`)::

            sage: macaulay2(matrix(QQ, 2, 0)).numrows()         # optional - macaulay2
            2
            sage: macaulay2(matrix(QQ, 0, 2)).numcols()         # optional - macaulay2
            2
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default
        if not self.nrows():
            return (macaulay2(self.base_ring())
                    .matrix('toList(%s:{})' % self.ncols()).transpose())
        entries = [list(row) for row in self]
        return macaulay2(self.base_ring()).matrix(entries)

    def _scilab_init_(self):
        """
        Return a string defining a Scilab representation of ``self``.

        EXAMPLES::

            sage: a = matrix([[1,2,3],[4,5,6],[7,8,9]]); a
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: a._scilab_init_()                                                     # needs sage.libs.pari
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
        Create a ScilabElement object based on ``self`` and returns it.

        EXAMPLES::

            sage: a = matrix([[1,2,3],[4,5,6],[7,8,9]]); a
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: b = scilab(a); b  # indirect doctest      # optional - scilab
                1.    2.    3.
                4.    5.    6.
                7.    8.    9.

        AUTHORS:

        - Ronan Paixao (2008-12-12)
        """
        return scilab(self._scilab_init_())

    def _sympy_(self):
        r"""
        Return a SymPy matrix corresponding to ``self``.

        OUTPUT:

        - An instance of either an ``ImmutableMatrix`` or ``ImmutableSparseMatrix``,
          regardless of whether ``self`` is mutable or not.

        EXAMPLES::

            sage: A = matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]); A
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: sA = A._sympy_(); sA                                                  # needs sympy
            Matrix([
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]])
            sage: type(sA)                                                              # needs sympy
            <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>

            sage: I = MatrixSpace(QQ, 5, 5, sparse=True).identity_matrix()
            sage: sI = I._sympy_(); sI                                                  # needs sympy
            Matrix([
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]])
            sage: type(sI)                                                              # needs sympy
            <class 'sympy.matrices.immutable.ImmutableSparseMatrix'>

        If ``self`` was immutable, then converting the result to Sage gives
        back ``self``::

            sage: immA = matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]], immutable=True)
            sage: immA._sympy_()._sage_() is immA                                       # needs sympy
            True

        If ``self`` was mutable, then converting back to Sage creates a new matrix::

            sage: sA._sage_() is A                                                      # needs sympy sage.symbolic
            False
            sage: sA._sage_() == A                                                      # needs sympy sage.symbolic
            True

        Symbolic matrices are supported::

            sage: # needs sympy sage.symbolic
            sage: M = matrix([[sin(x), cos(x)], [-cos(x), sin(x)]]); M
            [ sin(x)  cos(x)]
            [-cos(x)  sin(x)]
            sage: sM = M._sympy_(); sM
            Matrix([
            [ sin(x), cos(x)],
            [-cos(x), sin(x)]])
            sage: sM.subs(x, pi/4)
            Matrix([
            [ sqrt(2)/2, sqrt(2)/2],
            [-sqrt(2)/2, sqrt(2)/2]])

        TESTS:

        Dense 0-column/0-row matrices::

            sage: ZeroCol = matrix(QQ, 3, 0, sparse=False); ZeroCol
            []
            sage: sZeroCol = ZeroCol._sympy_(); sZeroCol                                # needs sympy
            Matrix(3, 0, [])

            sage: ZeroRow = matrix(QQ, 0, 2, sparse=False); ZeroRow
            []
            sage: sZeroRow = ZeroRow._sympy_(); sZeroRow                                # needs sympy
            Matrix(0, 2, [])
        """
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        from sympy.matrices import ImmutableMatrix, ImmutableSparseMatrix
        if self.is_sparse():
            matrix = ImmutableSparseMatrix(self.nrows(), self.ncols(),
                                           self.dict(copy=False))
        else:
            if not self.nrows() or not self.ncols():
                matrix = ImmutableMatrix(self.nrows(), self.ncols(), ())
            else:
                matrix = ImmutableMatrix(self.rows())
        if self.is_immutable():
            matrix._sage_object = self
        return matrix

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(matrix(QQ, 3, 3, [5..13])/7, verify=True)
            # Verified
            matrix(QQ, [[5/7, 6/7, 1], [8/7, 9/7, 10/7], [11/7, 12/7, 13/7]])
            sage: M = MatrixSpace(GF(5), 50, 50, sparse=True).random_element(density=0.002)
            sage: input = sage_input(M, verify=True)
            sage: sage_eval(input) == M
            True
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

    def numpy(self, dtype=None, copy=_MISSING):
        """
        Return the Numpy matrix associated to this matrix.

        INPUT:

        - ``dtype`` -- the desired data-type for the array. If not given,
          then the type will be determined as the minimum type required
          to hold the objects in the sequence.

        EXAMPLES::

            sage: # needs numpy
            sage: a = matrix(3, range(12))
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

            sage: import numpy                                                          # needs numpy
            sage: numpy.typecodes.items()                                               # needs numpy  # random
            [('All', '?bhilqpBHILQPefdgFDGSUVOMm'), ('AllFloat', 'efdgFDG'),
            ...

        For instance, you can see possibilities for real floating point numbers::

            sage: numpy.typecodes['Float']                                              # needs numpy
            'efdg'

        Alternatively, numpy automatically calls this function (via
        the magic :meth:`__array__` method) to convert Sage matrices
        to numpy arrays::

            sage: # needs numpy
            sage: import numpy
            sage: b = numpy.array(a); b
            array([[ 0,  1,  2,  3],
                   [ 4,  5,  6,  7],
                   [ 8,  9, 10, 11]])
            sage: b.dtype
            dtype('int32')  # 32-bit
            dtype('int64')  # 64-bit
            sage: b.shape
            (3, 4)

        TESTS::

            sage: # needs numpy
            sage: matrix(3, range(12)).numpy(copy=False)
            doctest:warning...
            DeprecationWarning: passing copy argument to numpy() is deprecated
            See https://github.com/sagemath/sage/issues/39152 for details.
            array([[ 0,  1,  2,  3],
                   [ 4,  5,  6,  7],
                   [ 8,  9, 10, 11]])
        """
        if copy is not _MISSING:
            from sage.misc.superseded import deprecation
            deprecation(39152, "passing copy argument to numpy() is deprecated")
        import numpy
        return numpy.asarray(self.list(), dtype=dtype).reshape(self.nrows(), self.ncols())

    def __array__(self, dtype=None, copy=None):
        """
        Define the magic ``__array__`` function so that ``numpy.array(m)`` can convert
        a matrix ``m`` to a numpy array. See
        `Interoperability with NumPy <https://numpy.org/doc/1.26/user/basics.interoperability.html>`_.

        Note that subclasses should override :meth:`numpy`, but usually not this method.

        INPUT:

        - ``dtype`` -- the desired data-type for the array. If not given,
          then the type will be determined automatically.

        - ``copy`` -- required for numpy 2.0 compatibility.
          See <https://numpy.org/devdocs/numpy_2_0_migration_guide.html#adapting-to-changes-in-the-copy-keyword>`_.
          Note that ``copy=False`` is not supported.

        TESTS::

            sage: # needs numpy
            sage: import numpy as np
            sage: a = matrix(3, range(12))
            sage: if np.lib.NumpyVersion(np.__version__) >= '2.0.0':
            ....:     try:
            ....:         np.array(a, copy=False)  # in numpy 2.0, this raises an error
            ....:     except ValueError:
            ....:         pass
            ....:     else:
            ....:         assert False
            ....: else:
            ....:     b = np.array(a, copy=False)  # in numpy 1.26, this means "avoid copy if possible"
            ....:     # https://numpy.org/doc/1.26/reference/generated/numpy.array.html#numpy.array
            ....:     # but no-copy is not supported so it will copy anyway
            ....:     a[0,0] = 1
            ....:     assert b[0,0] == 0
            ....:     b = np.asarray(a)
            ....:     a[0,0] = 2
            ....:     assert b[0,0] == 1
        """
        import numpy as np
        if np.lib.NumpyVersion(np.__version__) >= '2.0.0':
            if copy is False:
                raise ValueError("Sage matrix cannot be converted to numpy array without copying")
        else:
            assert copy is None  # numpy versions before 2.0 should not pass copy argument
        return self.numpy(dtype)

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
        Return lift of ``self`` to the covering ring of the base ring R,
        which is by definition the ring returned by calling
        cover_ring() on R, or just R itself if the cover_ring method
        is not defined.

        EXAMPLES::

            sage: M = Matrix(Integers(7), 2, 2, [5, 9, 13, 15]); M
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

    def lift_centered(self):
        """
        Apply the lift_centered method to every entry of ``self``.

        OUTPUT:

        If ``self`` is a matrix over the Integers mod `n`, this method returns the
        unique matrix `m` such that `m` is congruent to ``self`` mod `n` and for
        every entry `m[i,j]` we have `-n/2 < m[i,j] \\leq n/2`. If the
        coefficient ring does not have a cover_ring method, return ``self``.

        EXAMPLES::

            sage: M = Matrix(Integers(8), 2, 4, range(8)); M
            [0 1 2 3]
            [4 5 6 7]
            sage: L = M.lift_centered(); L
            [ 0  1  2  3]
            [ 4 -3 -2 -1]
            sage: parent(L)
            Full MatrixSpace of 2 by 4 dense matrices over Integer Ring

        The returned matrix is congruent to M modulo 8.::

            sage: L.mod(8)
            [0 1 2 3]
            [4 5 6 7]

        The field QQ doesn't have a cover_ring method::

            sage: hasattr(QQ, 'cover_ring')
            False

        So lifting a matrix over QQ gives back the same exact matrix.

        ::

            sage: B = matrix(QQ, 2, [1..4])
            sage: B.lift_centered()
            [1 2]
            [3 4]
            sage: B.lift_centered() is B
            True
        """
        try:
            S = self._base_ring.cover_ring()
            if S is not self._base_ring:
                return self.parent().change_ring(S)([v.lift_centered() for v in self])
        except AttributeError:
            pass
        return self

    #############################################################################################
    # rows, columns, sparse_rows, sparse_columns, dense_rows, dense_columns, row, column
    #############################################################################################
    cpdef row_ambient_module(self, base_ring=None, sparse=None):
        r"""
        Return the free module that contains the rows of the matrix.

        EXAMPLES::

            sage: M = matrix(Zmod(5), 2, 3)
            sage: M.row_ambient_module()
            Vector space of dimension 3 over Ring of integers modulo 5
            sage: M.row(1).parent() == M.row_ambient_module()
            True

            sage: M = Matrix(ZZ, 3, 4)
            sage: M.row_ambient_module()
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
            sage: M.row_ambient_module(QQ)
            Vector space of dimension 4 over Rational Field

            sage: M = Matrix(QQ, 4, 5)
            sage: M.row_ambient_module()
            Vector space of dimension 5 over Rational Field
            sage: M.row_ambient_module(ZZ)
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
        """
        if base_ring is not None or sparse is not None:
            if base_ring is None:
                base_ring = self._base_ring
            if sparse is None:
                sparse = self.is_sparse_c()
            if base_ring is self._base_ring and sparse == self.is_sparse_c():
                return self.row_ambient_module()
            return sage.modules.free_module.FreeModule(base_ring, self._ncols, sparse=sparse)

        x = self.fetch('row_ambient_module')
        if x is not None:
            return x
        x = sage.modules.free_module.FreeModule(self._base_ring, self._ncols,
                                                sparse=self.is_sparse_c())
        self.cache('row_ambient_module', x)
        return x

    def _row_ambient_module(self, base_ring=None):
        r"""
        TESTS::

            sage: M = matrix(Zmod(5), 2, 3)
            sage: M._row_ambient_module()
            doctest:warning
            ...
            DeprecationWarning: the method _row_ambient_module is deprecated use row_ambient_module (without underscore) instead
            See https://github.com/sagemath/sage/issues/32984 for details.
            Vector space of dimension 3 over Ring of integers modulo 5
        """
        from sage.misc.superseded import deprecation
        deprecation(32984, 'the method _row_ambient_module is deprecated use row_ambient_module (without underscore) instead')
        return self.row_ambient_module(base_ring)

    cpdef column_ambient_module(self, base_ring=None, sparse=None):
        r"""
        Return the free module that contains the columns of the matrix.

        EXAMPLES::

            sage: M = matrix(Zmod(5), 2, 3)
            sage: M.column_ambient_module()
            Vector space of dimension 2 over Ring of integers modulo 5
            sage: M.column(1).parent() == M.column_ambient_module()
            True

            sage: M = Matrix(ZZ, 3, 4)
            sage: M.column_ambient_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: M.column_ambient_module(QQ)
            Vector space of dimension 3 over Rational Field

            sage: M = Matrix(QQ, 4, 5)
            sage: M.column_ambient_module()
            Vector space of dimension 4 over Rational Field
            sage: M.column_ambient_module(ZZ)
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
        """
        if base_ring is not None or sparse is not None:
            if base_ring is None:
                base_ring = self._base_ring
            if sparse is None:
                sparse = self.is_sparse_c()
            if base_ring is self._base_ring and sparse == self.is_sparse_c():
                return self.column_ambient_module()
            return sage.modules.free_module.FreeModule(base_ring, self._nrows, sparse=sparse)

        x = self.fetch('column_ambient_module')
        if x is not None:
            return x
        x = sage.modules.free_module.FreeModule(self._base_ring, self._nrows,
                                                sparse=self.is_sparse_c())
        self.cache('column_ambient_module', x)
        return x

    def _column_ambient_module(self):
        r"""
        TESTS::

            sage: M = matrix(Zmod(5), 2, 3)
            sage: M._column_ambient_module()
            doctest:warning
            ...
            DeprecationWarning: the method _column_ambient_module is deprecated use column_ambient_module (without underscore) instead
            See https://github.com/sagemath/sage/issues/32984 for details.
            Vector space of dimension 2 over Ring of integers modulo 5
        """
        from sage.misc.superseded import deprecation
        deprecation(32984, 'the method _column_ambient_module is deprecated use column_ambient_module (without underscore) instead')
        return self.column_ambient_module()

    def columns(self, copy=True):
        r"""
        Return a list of the columns of ``self``.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``True``, return a copy of the list
          of columns which is safe to change

        If ``self`` is a sparse matrix, columns are returned as sparse vectors,
        otherwise returned vectors are dense.

        EXAMPLES::

            sage: matrix(3, [1..9]).columns()
            [(1, 4, 7), (2, 5, 8), (3, 6, 9)]
            sage: matrix(RR, 2, [sqrt(2), pi, exp(1), 0]).columns()                     # needs sage.symbolic
            [(1.41421356237310, 2.71828182845905), (3.14159265358979, 0.000000000000000)]
            sage: matrix(RR, 0, 2, []).columns()
            [(), ()]
            sage: matrix(RR, 2, 0, []).columns()
            []
            sage: m = matrix(RR, 3, 3, {(1,2): pi, (2, 2): -1, (0,1): sqrt(2)})         # needs sage.symbolic
            sage: parent(m.columns()[0])                                                # needs sage.symbolic
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
            ValueError: 'copy' must be ``True`` or False, not junk
        """
        if copy not in [True, False]:
            msg = "'copy' must be ``True`` or False, not {0}"
            raise ValueError(msg.format(copy))
        x = self.fetch('columns')
        if x is not None:
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
        Return a list of the rows of ``self``.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``True``, return a copy of the list
          of rows which is safe to change

        If ``self`` is a sparse matrix, rows are returned as sparse vectors,
        otherwise returned vectors are dense.

        EXAMPLES::

            sage: matrix(3, [1..9]).rows()
            [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
            sage: matrix(RR, 2, [sqrt(2), pi, exp(1), 0]).rows()                        # needs sage.symbolic
            [(1.41421356237310, 3.14159265358979), (2.71828182845905, 0.000000000000000)]
            sage: matrix(RR, 0, 2, []).rows()
            []
            sage: matrix(RR, 2, 0, []).rows()
            [(), ()]
            sage: m = matrix(RR, 3, 3, {(1,2): pi, (2, 2): -1, (0,1): sqrt(2)})         # needs sage.symbolic
            sage: parent(m.rows()[0])                                                   # needs sage.symbolic
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
            ValueError: 'copy' must be ``True`` or False, not junk
        """
        if copy not in [True, False]:
            msg = "'copy' must be ``True`` or False, not {0}"
            raise ValueError(msg.format(copy))
        x = self.fetch('rows')
        if x is not None:
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
        Return list of the dense columns of ``self``.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``True``, return a copy so you can
          modify it safely

        EXAMPLES:

        An example over the integers::

            sage: a = matrix(3, 3, range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.dense_columns()
            [(0, 3, 6), (1, 4, 7), (2, 5, 8)]

        We do an example over a polynomial ring::

            sage: R.<x> = QQ[]
            sage: a = matrix(R, 2, [x,x^2, 2/3*x,1+x^5]); a
            [      x     x^2]
            [  2/3*x x^5 + 1]
            sage: a.dense_columns()
            [(x, 2/3*x), (x^2, x^5 + 1)]
            sage: a = matrix(R, 2, [x,x^2, 2/3*x,1+x^5], sparse=True)
            sage: c = a.dense_columns(); c
            [(x, 2/3*x), (x^2, x^5 + 1)]
            sage: parent(c[1])
            Ambient free module of rank 2 over the principal ideal domain
             Univariate Polynomial Ring in x over Rational Field

        TESTS:

        Check that the returned rows are immutable as per :issue:`14874`::

            sage: m = Mat(ZZ, 3, 3)(range(9))
            sage: v = m.dense_columns()
            sage: [x.is_mutable() for x in v]
            [False, False, False]
        """
        x = self.fetch('dense_columns')
        if x is not None:
            if copy: return list(x)
            return x
        cdef Py_ssize_t i
        A = self if self.is_dense() else self.dense_matrix()
        C = [A.column(i) for i in range(self._ncols)]

        # Make the vectors immutable since we are caching them
        for x in C:
            x.set_immutable()

        # cache result
        self.cache('dense_columns', C)
        if copy:
            return list(C)
        else:
            return C

    def dense_rows(self, copy=True):
        """
        Return list of the dense rows of ``self``.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``True``, return a copy so you can
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

        Check that the returned rows are immutable as per :issue:`14874`::

            sage: m = Mat(ZZ, 3, 3)(range(9))
            sage: v = m.dense_rows()
            sage: [x.is_mutable() for x in v]
            [False, False, False]
        """
        x = self.fetch('dense_rows')
        if x is not None:
            if copy: return list(x)
            return x

        cdef Py_ssize_t i
        A = self if self.is_dense() else self.dense_matrix()
        R = [A.row(i) for i in range(self._nrows)]

        # Make the vectors immutable since we are caching them
        for x in R:
            x.set_immutable()

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

        - ``copy`` -- boolean (default: ``True``); if ``True``, return a copy so you can
          modify it safely

        EXAMPLES::

            sage: a = matrix(2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: v = a.sparse_columns(); v
            [(0, 3), (1, 4), (2, 5)]
            sage: v[1].is_sparse()
            True

        TESTS:

        Columns of sparse matrices having no columns were fixed on :issue:`10714`::

            sage: m = matrix(10, 0, sparse=True)
            sage: m.ncols()
            0
            sage: m.columns()
            []

        Check that the returned columns are immutable as per :issue:`14874`::

            sage: m = Mat(ZZ, 3, 3, sparse=True)(range(9))
            sage: v = m.sparse_columns()
            sage: [x.is_mutable() for x in v]
            [False, False, False]
        """
        x = self.fetch('sparse_columns')
        if x is not None:
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
        for x in C:
            x.set_immutable()

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

        - ``copy`` -- boolean (default: ``True``); if ``True``, return a copy so you can
          modify it safely

        EXAMPLES::

            sage: m = Mat(ZZ, 3, 3, sparse=True)(range(9)); m
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

        Rows of sparse matrices having no rows were fixed on :issue:`10714`::

            sage: m = matrix(0, 10, sparse=True)
            sage: m.nrows()
            0
            sage: m.rows()
            []

        Check that the returned rows are immutable as per :issue:`14874`::

            sage: m = Mat(ZZ, 3, 3, sparse=True)(range(9))
            sage: v = m.sparse_rows()
            sage: [x.is_mutable() for x in v]
            [False, False, False]
        """
        x = self.fetch('sparse_rows')
        if x is not None:
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
        for x in R:
            x.set_immutable()

        # cache and return result
        self.cache('sparse_rows', R)
        if copy:
            return list(R)
        else:
            return R

    def column(self, Py_ssize_t i, from_list=False):
        """
        Return the ``i``-th column of this matrix as a vector.

        This column is a dense vector if and only if the matrix is a dense
        matrix.

        INPUT:

        - ``i`` -- integer

        - ``from_list`` -- boolean (default: ``False``); if ``True``, returns the
          ``i``-th element of ``self.columns()`` (see :func:`columns()`),
          which may be faster, but requires building a list of all
          columns the first time it is called after an entry of the
          matrix is changed.

        EXAMPLES::

            sage: a = matrix(2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.column(1)
            (1, 4)

        If the column is negative, it wraps around, just like with list
        indexing, e.g., -1 gives the right-most column::

            sage: a.column(-1)
            (2, 5)

        TESTS::

            sage: a = matrix(2, 3, range(6)); a
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
            raise IndexError("matrix has no columns")
        if i >= self._ncols or i < -self._ncols:
            raise IndexError("column index out of range")
        i = i % self._ncols
        if i < 0:
            i = i + self._ncols
        if from_list:
            return self.columns(copy=False)[i]
        cdef Py_ssize_t j
        V = self.column_ambient_module()
        tmp = [self.get_unsafe(j, i) for j in range(self._nrows)]
        return V(tmp, coerce=False, copy=False, check=False)

    def row(self, Py_ssize_t i, from_list=False):
        """
        Return the ``i``-th row of this matrix as a vector.

        This row is a dense vector if and only if the matrix is a dense
        matrix.

        INPUT:

        - ``i`` -- integer

        - ``from_list`` -- boolean (default: ``False``); if ``True``, returns the
          ``i``-th element of ``self.rows()`` (see :func:`rows`), which
          may be faster, but requires building a list of all rows the
          first time it is called after an entry of the matrix is
          changed.

        EXAMPLES::

            sage: a = matrix(2, 3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.row(0)
            (0, 1, 2)
            sage: a.row(1)
            (3, 4, 5)
            sage: a.row(-1)  # last row
            (3, 4, 5)

        TESTS::

            sage: a = matrix(2, 3, range(6)); a
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
            raise IndexError("matrix has no rows")
        if i >= self._nrows or i < -self._nrows:
            raise IndexError("row index out of range")
        i = i % self._nrows
        if i < 0:
            i = i + self._nrows
        if from_list:
            return self.rows(copy=False)[i]
        cdef Py_ssize_t j
        V = self.row_ambient_module()
        tmp = [self.get_unsafe(i,j) for j in range(self._ncols)]
        return V(tmp, coerce=False, copy=False, check=False)

    ###########################################################################
    # Building matrices out of other matrices, rows, or columns
    ###########################################################################

    def stack(self, bottom, subdivide=False):
        r"""
        Return a new matrix formed by appending the matrix (or vector)
        ``bottom`` below ``self``::

            [  self  ]
            [ bottom ]

        INPUT:

        - ``bottom`` -- a matrix, vector or free module element, whose
          dimensions are compatible with ``self``

        - ``subdivide`` -- (default: ``False``) request the resulting
          matrix to have a new subdivision, separating ``self`` from ``bottom``

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
            TypeError: number of columns must be the same, not 2 and 3

            sage: v = vector(RR, [100, 200, 300])
            sage: A.stack(v)
            Traceback (most recent call last):
            ...
            TypeError: number of columns must be the same, not 2 and 3

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

        The base ring of the result is the common parent for the base
        rings of ``self`` and ``bottom``. In particular, the parent for
        ``A.stack(B)`` and ``B.stack(A)`` should be equal::

            sage: A = matrix(QQ, 1, 2, [1,2])
            sage: B = matrix(RR, 1, 2, [sin(1.1), sin(2.2)])
            sage: C = A.stack(B); C
            [ 1.00000000000000  2.00000000000000]
            [0.891207360061435 0.808496403819590]
            sage: C.parent()
            Full MatrixSpace of 2 by 2 dense matrices
             over Real Field with 53 bits of precision

            sage: D = B.stack(A); D
            [0.891207360061435 0.808496403819590]
            [ 1.00000000000000  2.00000000000000]
            sage: D.parent()
            Full MatrixSpace of 2 by 2 dense matrices
             over Real Field with 53 bits of precision

        ::

            sage: R.<y> = PolynomialRing(ZZ)
            sage: A = matrix(QQ, 1, 2, [1, 2/3])
            sage: B = matrix(R, 1, 2, [y, y^2])

            sage: C = A.stack(B); C
            [  1 2/3]
            [  y y^2]
            sage: C.parent()
            Full MatrixSpace of 2 by 2 dense matrices over
             Univariate Polynomial Ring in y over Rational Field

        Stacking a dense matrix atop a sparse one returns a sparse
        matrix::

            sage: M = Matrix(ZZ, 2, 3, range(6), sparse=False)
            sage: N = diagonal_matrix([10,11,12], sparse=True)
            sage: P = M.stack(N); P
            [ 0  1  2]
            [ 3  4  5]
            [10  0  0]
            [ 0 11  0]
            [ 0  0 12]
            sage: P.is_sparse()
            True
            sage: P = N.stack(M); P
            [10  0  0]
            [ 0 11  0]
            [ 0  0 12]
            [ 0  1  2]
            [ 3  4  5]
            sage: P.is_sparse()
            True

        One can stack matrices over different rings (:issue:`16399`). ::

            sage: M = Matrix(ZZ, 2, 3, range(6))
            sage: N = Matrix(QQ, 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]
            sage: N.stack(M)
            [10 11 12]
            [ 0  1  2]
            [ 3  4  5]

        TESTS:

        A legacy test from the original implementation.  ::

            sage: M = Matrix(QQ, 2, 3, range(6))
            sage: N = Matrix(QQ, 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]

        Non-matrices fail gracefully::

            sage: M.stack(polygen(QQ))
            Traceback (most recent call last):
            ...
            TypeError: a matrix must be stacked with another matrix or a vector

        AUTHORS:

        - Rob Beezer (2011-03-19): rewritten to mirror code for :meth:`augment`

        - Jeroen Demeyer (2015-01-06): refactor, see :issue:`16399`.
          Put all boilerplate in one place (here) and put the actual
          type-dependent implementation in ``_stack_impl``.
        """
        cdef Matrix other
        if isinstance(bottom, Matrix):
            other = <Matrix>bottom
        else:
            if hasattr(bottom, '_vector_'):
                bottom = bottom.row()
            else:
                raise TypeError('a matrix must be stacked with '
                        'another matrix or a vector')
            other = <Matrix?>bottom

        if self._ncols != other._ncols:
            raise TypeError("number of columns must be the same, not %s and %s" %
                    (self.ncols(), bottom.ncols()) )

        top_ring = self._base_ring
        bottom_ring = other._base_ring
        if top_ring is not bottom_ring:
            R = coercion_model.common_parent(top_ring, bottom_ring)
            if top_ring is not R:
                self = self.change_ring(R)
            if bottom_ring is not R:
                other = other.change_ring(R)

        if type(self) is not type(other):
            # If one of the matrices is sparse, return a sparse matrix
            if self.is_sparse_c() and not other.is_sparse_c():
                other = other.sparse_matrix()
            elif other.is_sparse_c() and not self.is_sparse_c():
                self = self.sparse_matrix()

        Z = self._stack_impl(other)
        if subdivide:
            Z._subdivide_on_stack(self, other)
        return Z

    cdef _stack_impl(self, bottom):
        """
        Implementation of :meth:`stack`.

        Assume that ``self`` and ``other`` are compatible in the sense
        that they have the same base ring and that both are either
        dense or sparse.
        """
        cdef Matrix other = <Matrix>bottom
        cdef Matrix Z
        Z = self.new_matrix(nrows=self._nrows + other._nrows, ncols=self._ncols)

        cdef Py_ssize_t r, c
        cdef Py_ssize_t nr = self._nrows
        for r in range(self._nrows):
            for c in range(self._ncols):
                Z.set_unsafe(r, c, self.get_unsafe(r,c))
        for r in range(other._nrows):
            for c in range(other._ncols):
                Z.set_unsafe(r+nr, c, other.get_unsafe(r,c))

        return Z

    def augment(self, right, subdivide=False):
        r"""
        Return a new matrix formed by appending the matrix (or vector)
        ``right`` on the right side of ``self``.

        INPUT:

        - ``right`` -- a matrix, vector or free module element, whose
          dimensions are compatible with ``self``

        - ``subdivide`` -- (default: ``False``) request the resulting
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
            sage: C = A.augment(B); C                                                   # needs sage.symbolic
            [                  1 183017397/205358938]
            [                  2 106580492/131825561]
            sage: C.parent()                                                            # needs sage.symbolic
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: D = B.augment(A); D
            [0.89120736006...  1.00000000000000]
            [0.80849640381...  2.00000000000000]
            sage: D.parent()
            Full MatrixSpace of 2 by 2 dense matrices
             over Real Field with 53 bits of precision

        Sometimes it is not possible to coerce into the base ring of
        ``self``.  A solution is to change the base ring of ``self`` to
        a more expansive ring.  Here we mix the rationals with a ring of
        polynomials with rational coefficients.  ::

            sage: R.<y> = PolynomialRing(QQ)
            sage: A = matrix(QQ, 1, [1,2])
            sage: B = matrix(R, 1, [y, y^2])

            sage: C = B.augment(A); C
            [  y y^2   1   2]
            sage: C.parent()
            Full MatrixSpace of 1 by 4 dense matrices over
             Univariate Polynomial Ring in y over Rational Field

            sage: D = A.augment(B)
            Traceback (most recent call last):
            ...
            TypeError: y is not a constant polynomial

            sage: E = A.change_ring(R)
            sage: F = E.augment(B); F
            [  1   2   y y^2]
            sage: F.parent()
            Full MatrixSpace of 1 by 4 dense matrices over
             Univariate Polynomial Ring in y over Rational Field

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples
        - Rob Beezer (2010-12-07): vector argument, docstring, subdivisions
        """
        from sage.matrix.constructor import matrix

        if not isinstance(right, sage.matrix.matrix1.Matrix):
            if hasattr(right, '_vector_'):
                right = right.column()
            else:
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
        Return the matrix constructed from ``self`` using columns with indices
        in the columns list.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8), 3, 3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_columns([2,1])
            [2 1]
            [5 4]
            [0 7]
        """
        columns = PySequence_Fast(columns, "columns is not iterable")
        cdef Py_ssize_t ncols = len(columns)

        # Construct new matrix
        cdef Matrix A = self.new_matrix(ncols=ncols)
        cdef Py_ssize_t i, j, col
        for j, col in enumerate(columns):
            if col < 0 or col >= self._ncols:
                raise IndexError("column index out of range")
            for i in range(self._nrows):
                A.set_unsafe(i, j, self.get_unsafe(i, col))
        return A

    def delete_columns(self, dcols, check=True):
        """
        Return the matrix constructed from deleting the columns with indices in the ``dcols`` list.

        INPUT:

        - ``dcols`` -- list of indices of columns to be deleted from ``self``
        - ``check`` -- boolean (default: ``True``); check whether any index in
          ``dcols`` is out of range

        .. SEEALSO::

            The methods :meth:`delete_rows` and :meth:`matrix_from_columns`
            are related.

        EXAMPLES::

            sage: A = Matrix(3, 4, range(12)); A
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
            IndexError: [-1, 4] contains invalid indices
            sage: A.delete_columns([-1,2,4], check=False)
            [ 0  1  3]
            [ 4  5  7]
            [ 8  9 11]

        TESTS:

        The list of indices is checked.  ::

            sage: A.delete_columns("junk")
            Traceback (most recent call last):
            ...
            IndexError: ['j', 'k', 'n', 'u'] contains invalid indices

        AUTHORS:

        - Wai Yan Pong (2012-03-05)
        """
        if not isinstance(dcols, (list, tuple)):
            dcols = list(dcols)

        cdef list cols, diff_cols

        if check:
            diff_cols = sorted(set(dcols).difference(set(range(self._ncols))))
            if diff_cols:
                raise IndexError("{d} contains invalid indices".format(d=diff_cols))
        cols = [k for k in range(self._ncols) if k not in dcols]
        return self.matrix_from_columns(cols)

    def matrix_from_rows(self, rows):
        """
        Return the matrix constructed from ``self`` using rows with indices in
        the rows list.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8), 3, 3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows([2,1])
            [6 7 0]
            [3 4 5]
        """
        rows = PySequence_Fast(rows, "rows is not iterable")
        cdef Py_ssize_t nrows = len(rows)

        # Construct new matrix
        cdef Matrix A = self.new_matrix(nrows=nrows)
        cdef Py_ssize_t i, j, row
        for i, row in enumerate(rows):
            if row < 0 or row >= self._nrows:
                raise IndexError("row index out of range")
            for j in range(self._ncols):
                A.set_unsafe(i, j, self.get_unsafe(row, j))
        return A

    def delete_rows(self, drows, check=True):
        """
        Return the matrix constructed from deleting the rows with indices in the ``drows`` list.

        INPUT:

        - ``drows`` -- list of indices of rows to be deleted from ``self``
        - ``check`` -- boolean (default: ``True``); whether to check if any
          index in ``drows`` is out of range

        .. SEEALSO::

            The methods :meth:`delete_columns` and :meth:`matrix_from_rows`
            are related.

        EXAMPLES::

            sage: A = Matrix(4, 3, range(12)); A
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
            IndexError: [-1, 4] contains invalid indices
            sage: A.delete_rows([-1,2,4], check=False)
            [ 0  1  2]
            [ 3  4  5]
            [ 9 10 11]

        TESTS:

        The list of indices is checked.  ::

            sage: A.delete_rows("junk")
            Traceback (most recent call last):
            ...
            IndexError: ['j', 'k', 'n', 'u'] contains invalid indices

        AUTHORS

        - Wai Yan Pong (2012-03-05)
        """
        if not isinstance(drows, (list, tuple)):
            drows = list(drows)

        cdef list rows, diff_rows

        if check:
            diff_rows = sorted(set(drows).difference(set(range(self._nrows))))
            if diff_rows:
                raise IndexError("{d} contains invalid indices".format(d=diff_rows))
        rows = [k for k in range(self._nrows) if k not in drows]
        return self.matrix_from_rows(rows)

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from ``self`` from the given rows and
        columns.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8), 3, 3)
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
        this 3 times::

            sage: A.matrix_from_rows_and_columns([1,1,1], [2,0,0])
            [5 3 3]
            [5 3 3]
            [5 3 3]

        AUTHORS:

        - Jaap Spies (2006-02-18)

        - Didier Deshommes: some Pyrex speedups implemented
        """
        rows = PySequence_Fast(rows, "rows is not iterable")
        columns = PySequence_Fast(columns, "columns is not iterable")

        cdef Py_ssize_t ncols = len(columns)
        cdef Py_ssize_t nrows = len(rows)

        # Check whether column indices are valid
        cdef Py_ssize_t i, j, row, col
        for col in columns:
            if col < 0 or col >= self._ncols:
                raise IndexError("column index out of range")

        # Construct new matrix
        cdef Matrix A = self.new_matrix(nrows=nrows, ncols=ncols)
        for i, row in enumerate(rows):
            if row < 0 or row >= self._nrows:
                raise IndexError("row index out of range")
            for j, col in enumerate(columns):
                A.set_unsafe(i, j, self.get_unsafe(row, col))
        return A

    def submatrix(self, Py_ssize_t row=0, Py_ssize_t col=0,
                        Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        """
        Return the matrix constructed from ``self`` using the specified
        range of rows and columns.

        INPUT:

        - ``row``, ``col`` -- index of the starting row and column (indices
          start at zero)

        - ``nrows``, ``ncols`` -- (optional) number of rows and columns to
          take. If not provided, take all rows below and all columns to
          the right of the starting entry.

        .. SEEALSO::

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
        Set the entries of row ``row`` to the entries of ``v``.

        INPUT:

        - ``row`` -- index of row to be set

        - ``v`` -- list or vector of the new entries

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
        Set the entries of column ``col`` to the entries of ``v``.

        INPUT:

        - ``col`` -- index of column to be set

        - ``v`` -- list or vector of the new entries

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
            TypeError: Cannot set column with Rational Field elements
            over Integer Ring, use change_ring first.
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

    def zero_pattern_matrix(self, ring=None):
        """
        Return a matrix that contains one for corresponding zero entries.

        All other entries are zero.

        INPUT:

        - ``ring`` -- (optional); base ring of the output; default is ``ZZ``

        OUTPUT:

        A new dense matrix with same dimensions as ``self``
        and with base ring ``ring``.

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, [1,2,-2,0])
            sage: M.zero_pattern_matrix()
            [0 0]
            [0 1]

            sage: M = Matrix(QQ, 2, [1,2/3,-2,0])
            sage: M.zero_pattern_matrix()
            [0 0]
            [0 1]

        Default base ring for the output is ``ZZ``::

            sage: M.zero_pattern_matrix().base_ring()
            Integer Ring

        Specify a different base ring for the output::

            sage: M.zero_pattern_matrix(GF(2)).base_ring()
            Finite Field of size 2

        Examples for different base rings for ``self``::

            sage: M = Matrix(Zmod(8), 3, 2, [2, 3, 9, 8, 1, 0]); M
            [2 3]
            [1 0]
            [1 0]
            sage: M.zero_pattern_matrix()
            [0 0]
            [0 1]
            [0 1]

        ::

            sage: W.<a> = CyclotomicField(100)                                          # needs sage.rings.number_field
            sage: M = Matrix(2, 3, [a, a/2, 0, a^2, a^100-1, a^2 - a]); M               # needs sage.rings.number_field
            [      a   1/2*a       0]
            [    a^2       0 a^2 - a]
            sage: M.zero_pattern_matrix()                                               # needs sage.rings.number_field
            [0 0 1]
            [0 1 0]

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(2^4)
            sage: l = [a^2 + 1, a^3 + 1, 0, 0, a, a^3 + a + 1, a + 1,
            ....:      a + 1, a^2, a^3 + a + 1, a^3 + a, a^3 + a]
            sage: M = Matrix(K, 3, 4, l); M
            [    a^2 + 1     a^3 + 1           0           0]
            [          a a^3 + a + 1       a + 1       a + 1]
            [        a^2 a^3 + a + 1     a^3 + a     a^3 + a]
            sage: M.zero_pattern_matrix()
            [0 0 1 1]
            [0 0 0 0]
            [0 0 0 0]

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(25)
            sage: M = Matrix(K, 2, 3, [0, 2, 3, 5, a, a^2])
            sage: M
            [    0     2     3]
            [    0     a a + 3]
            sage: M.zero_pattern_matrix()
            [1 0 0]
            [1 0 0]

        .. NOTE::

            This method can be optimized by improving
            :meth:`get_is_zero_unsafe` for derived matrix classes.
        """
        if ring is None:
            from sage.rings.integer_ring import ZZ
            ring = ZZ

        cdef object one = ring.one()
        cdef Py_ssize_t i, j

        from sage.matrix.matrix_space import MatrixSpace
        MZ = MatrixSpace(ring, self._nrows, self._ncols, sparse=False)
        cdef Matrix M = MZ(ring.zero())

        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if self.get_is_zero_unsafe(i, j):
                    M.set_unsafe(i, j, one)
        return M

    ######################################################################
    # Change of representation between dense and sparse.
    ######################################################################

    def dense_matrix(self):
        """
        If this matrix is sparse, return a dense matrix with the same
        entries. If this matrix is dense, return this matrix (not a copy).

        .. NOTE::

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
            sage: A == B
            True
            sage: B.dense_matrix() is B
            True
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

        Ensure we can compute the correct dense matrix even if the
        dict items are ETuples (see :issue:`17658`)::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: matrix(GF(5^2, "z"), {ETuple((1, 1)): 2}).dense_matrix()              # needs sage.rings.finite_rings
            [0 0]
            [0 2]
        """
        if self.is_dense():
            return self
        cdef Matrix A
        A = self.new_matrix(self._nrows, self._ncols, self,
                coerce=False, sparse=False)
        if self._subdivisions is not None:
            A.subdivide(self.subdivisions())
        return A

    def sparse_matrix(self):
        """
        If this matrix is dense, return a sparse matrix with the same
        entries. If this matrix is sparse, return this matrix (not a
        copy).

        .. NOTE::

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
            sage: A == B
            True
            sage: B.sparse_matrix() is B
            True
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
        A = self.new_matrix(self._nrows, self._ncols, self,
                coerce=False, sparse=True)
        if self._subdivisions is not None:
            A.subdivide(self.subdivisions())
        return A

    def matrix_space(self, nrows=None, ncols=None, sparse=None):
        """
        Return the ambient matrix space of ``self``.

        INPUT:

        - ``nrows``, ``ncols`` -- (optional) number of rows and columns in
          returned matrix space

        - ``sparse`` -- whether the returned matrix space uses sparse or
          dense matrices

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

            sage: M = MatrixSpace(QQ, 3, implementation='generic')
            sage: m = M.an_element()
            sage: m.matrix_space()
            Full MatrixSpace of 3 by 3 dense matrices over Rational Field
             (using Matrix_generic_dense)
            sage: m.matrix_space(nrows=2, ncols=12)
            Full MatrixSpace of 2 by 12 dense matrices over Rational Field
             (using Matrix_generic_dense)
            sage: m.matrix_space(nrows=2, sparse=True)
            Full MatrixSpace of 2 by 3 sparse matrices over Rational Field
        """
        if nrows is None:
            nrows = self._nrows
        if ncols is None:
            ncols = self._ncols
        if sparse is None:
            sparse = self.is_sparse()

        base_ring = self._base_ring

        if nrows == self._nrows and ncols == self._ncols and sparse == self.is_sparse():
            return self._parent
        else:
            if sparse == self.is_sparse():
                implementation = self.__class__
            else:
                implementation = None
            from sage.matrix.matrix_space import MatrixSpace
            return MatrixSpace(base_ring, nrows, ncols, sparse, implementation)

    def new_matrix(self, nrows=None, ncols=None, entries=None,
                   coerce=True, copy=True, sparse=None):
        """
        Create a matrix in the parent of this matrix with the given number
        of rows, columns, etc. The default parameters are the same as for
        ``self``.

        INPUT:

        These three variables get sent to :func:`matrix_space`:

        - ``nrows``, ``ncols`` -- number of rows and columns in returned
          matrix. If not specified, defaults to ``None`` and will give a
          matrix of the same size as ``self``.
        - ``sparse`` -- whether returned matrix is sparse or not. Defaults
          to same value as self

        The remaining three variables (``coerce``, ``entries``, and
        ``copy``) are used by
        :func:`sage.matrix.matrix_space.MatrixSpace` to construct the
        new matrix.

        .. warning::

           This function called with no arguments returns the zero
           matrix of the same dimension and sparseness of ``self``.

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
            Full MatrixSpace of 2 by 3 dense matrices
             over Real Field with 53 bits of precision

        ::

            sage: M = MatrixSpace(ZZ, 2, 3, implementation='generic')
            sage: m = M.an_element()
            sage: m.new_matrix().parent()
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
             (using Matrix_generic_dense)
            sage: m.new_matrix(3,3).parent()
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
             (using Matrix_generic_dense)
            sage: m.new_matrix(3,3, sparse=True).parent()
            Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
        """
        if (sparse is None or self.is_sparse() == sparse):
            if self._nrows == nrows and self._ncols == ncols:
                return self._parent(entries, coerce=coerce, copy=copy)
            elif self._nrows == ncols and self._ncols == nrows:
                return self._parent.transposed(entries, coerce=coerce, copy=copy)
        return self.matrix_space(nrows, ncols, sparse)(entries, coerce=coerce, copy=copy)

    def block_sum(self, Matrix other):
        """
        Return the block matrix that has ``self`` and ``other`` on the diagonal::

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
            raise TypeError("other must be a Matrix")
        top = self.augment(self.new_matrix(ncols=other._ncols))
        bottom = other.new_matrix(ncols=self._ncols).augment(other)
        return top.stack(bottom)
