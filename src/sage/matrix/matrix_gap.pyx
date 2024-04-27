r"""
Wrappers on GAP matrices
"""
# ****************************************************************************
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.fields import Fields
from sage.libs.gap.libgap import libgap
from sage.structure.element cimport Matrix
from sage.matrix.args cimport MatrixArgs_init


cdef class Matrix_gap(Matrix_dense):
    r"""
    A Sage matrix wrapper over a GAP matrix.

    EXAMPLES::

        sage: M = MatrixSpace(ZZ, 2, implementation='gap')
        sage: m1 = M([1, 0, 2, -3])
        sage: m2 = M([2, 2, 5, -1])
        sage: type(m1)
        <class 'sage.matrix.matrix_gap.Matrix_gap'>

        sage: m1 * m2
        [  2   2]
        [-11   7]
        sage: type(m1 * m2)
        <class 'sage.matrix.matrix_gap.Matrix_gap'>

        sage: M = MatrixSpace(QQ, 5, 3, implementation='gap')
        sage: m = M(range(15))
        sage: m.left_kernel()
        Vector space of degree 5 and dimension 3 over Rational Field
        Basis matrix:
        [ 1  0  0 -4  3]
        [ 0  1  0 -3  2]
        [ 0  0  1 -2  1]

        sage: M = MatrixSpace(ZZ, 10, implementation='gap')
        sage: m = M(range(100))
        sage: m.transpose().parent() is M
        True

        sage: # needs sage.rings.number_field
        sage: UCF = UniversalCyclotomicField()
        sage: M = MatrixSpace(UCF, 3, implementation='gap')
        sage: m = M([UCF.zeta(i) for i in range(1,10)])
        sage: m
        [               1               -1             E(3)]
        [            E(4)             E(5)          -E(3)^2]
        [            E(7)             E(8) -E(9)^4 - E(9)^7]
        sage: (m^2)[1,2]
        E(180)^32 - E(180)^33 + E(180)^68 - E(180)^69 + E(180)^104 - E(180)^141 - E(180)^156 + E(180)^176 - E(180)^177

    TESTS::

        sage: rings = [ZZ, QQ, UniversalCyclotomicField(), GF(2), GF(3)]
        sage: rings += [UniversalCyclotomicField()]                                     # needs sage.rings.number_field
        sage: for ring in rings:
        ....:     M = MatrixSpace(ring, 2, implementation='gap')
        ....:     TestSuite(M).run(skip=['_test_construction'])
        ....:     M = MatrixSpace(ring, 2, 3, implementation='gap')
        ....:     TestSuite(M).run(skip=['_test_construction'])
    """
    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``True`` (the default), convert elements to the
          base ring before passing them to GAP. If ``False``, pass the
          elements to GAP as given.

        TESTS::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: M(0)
            [0 0]
            [0 0]
            sage: M(1)
            [1 0]
            [0 1]
            sage: M(2)
            [2 0]
            [0 2]
            sage: type(M(0))
            <class 'sage.matrix.matrix_gap.Matrix_gap'>
            sage: type(M(1))
            <class 'sage.matrix.matrix_gap.Matrix_gap'>
            sage: type(M(2))
            <class 'sage.matrix.matrix_gap.Matrix_gap'>

            sage: M = MatrixSpace(QQ, 2, 3, implementation='gap')
            sage: M(0)
            [0 0 0]
            [0 0 0]
            sage: M(1)
            Traceback (most recent call last):
            ...
            TypeError: nonzero scalar matrix must be square
            sage: MatrixSpace(QQ, 1, 2, implementation='gap')(0)
            [0 0]
            sage: MatrixSpace(QQ, 2, 1, implementation='gap')(0)
            [0]
            [0]
        """
        ma = MatrixArgs_init(parent, entries)
        it = ma.iter(coerce)
        cdef list mat = []
        cdef long i, j
        for i in range(ma.nrows):
            row = [next(it) for j in range(ma.ncols)]
            mat.append(row)
        self._libgap = libgap(mat)

    cdef Matrix_gap _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        if nrows == self._nrows and ncols == self._ncols:
            P = self._parent
        else:
            P = self.matrix_space(nrows, ncols)

        return Matrix_gap.__new__(Matrix_gap, P)

    def __copy__(self):
        r"""
        TESTS::

            sage: M = MatrixSpace(QQ, 2, implementation='gap')
            sage: m1 = M([1,2,0,3])
            sage: m2 = m1.__copy__()
            sage: m2
            [1 2]
            [0 3]
            sage: m1[0,1] = -2
            sage: m1
            [ 1 -2]
            [ 0  3]
            sage: m2
            [1 2]
            [0 3]
        """
        cdef Matrix_gap M = self._new(self._nrows, self._ncols)
        M._libgap = self._libgap.deepcopy(1)
        return M

    def __reduce__(self):
        r"""
        TESTS::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: m = M([1,2,1,2])
            sage: loads(dumps(m)) == m
            True
        """
        return self._parent, (self.list(),)

    cpdef GapElement gap(self):
        r"""
        Return the underlying gap object.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: m = M([1,2,2,1]).gap()
            sage: m
            [ [ 1, 2 ], [ 2, 1 ] ]
            sage: type(m)
            <class 'sage.libs.gap.element.GapElement_List'>

            sage: m.MatrixAutomorphisms()
            Group([ (1,2) ])
        """
        return self._libgap

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._base_ring(self._libgap[i,j])

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        r"""
        TESTS::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: m = M(0)
            sage: m[0,1] = 13
            sage: m
            [ 0 13]
            [ 0  0]
            sage: m[1,0] = -1/2
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
        """
        self._libgap[i,j] = x

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare ``self`` and ``right``.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: m1 = M([1,2,3,4])
            sage: m2 = M([1,2,3,4])
            sage: m3 = M([1,2,0,4])
            sage: m1 == m2
            True
            sage: m1 != m2
            False
            sage: m1 == m3
            False
            sage: m1 != m3
            True

            sage: M = MatrixSpace(QQ, 2, implementation='gap')
            sage: m1 = M([1/2, 1/3, 2, -5])
            sage: m2 = M([1/2, 1/3, 2, -5])
            sage: m3 = M([1/2, 0, 2, -5])
            sage: m1 == m2
            True
            sage: m1 != m2
            False
            sage: m1 == m3
            False
            sage: m1 != m3
            True

            sage: # needs sage.rings.number_field
            sage: UCF = UniversalCyclotomicField()
            sage: M = MatrixSpace(UCF, 2, implementation='gap')
            sage: m1 = M([E(2), E(3), 0, E(4)])
            sage: m2 = M([E(2), E(3), 0, E(4)])
            sage: m3 = M([E(2), E(3), 0, E(5)])
            sage: m1 == m2
            True
            sage: m1 != m2
            False
            sage: m1 == m3
            False
            sage: m1 != m3
            True
        """
        return (<Matrix_gap> self)._libgap._richcmp_((<Matrix_gap> other)._libgap, op)

    def __neg__(self):
        r"""
        TESTS::

            sage: M = MatrixSpace(ZZ, 2, 3, implementation='gap')
            sage: m = M([1, -1, 3, 2, -5, 1])
            sage: -m
            [-1  1 -3]
            [-2  5 -1]
        """
        cdef Matrix_gap M = self._new(self._nrows, self._ncols)
        M._libgap = self._libgap.AdditiveInverse()
        return M

    def __invert__(self):
        r"""
        TESTS::

            sage: M = MatrixSpace(QQ, 2, implementation='gap')
            sage: ~M([4,2,2,2])
            [ 1/2 -1/2]
            [-1/2    1]
        """
        cdef Matrix_gap M
        if self._base_ring in Fields():
            M = self._new(self._nrows, self._ncols)
            M._libgap = self._libgap.Inverse()
            return M
        return Matrix_dense.__invert__(self)

    cpdef _add_(left, right):
        r"""
        TESTS::

            sage: M = MatrixSpace(ZZ, 2, 3, implementation='gap')
            sage: M([1,2,3,4,3,2]) + M([1,1,1,1,1,1]) == M([2,3,4,5,4,3])
            True
        """
        cdef Matrix_gap cleft = <Matrix_gap> left
        cdef Matrix_gap ans = cleft._new(cleft._nrows, cleft._ncols)
        ans._libgap = left._libgap + (<Matrix_gap> right)._libgap
        return ans

    cpdef _sub_(left, right):
        r"""
        TESTS::

            sage: M = MatrixSpace(ZZ, 2, 3, implementation='gap')
            sage: M([1,2,3,4,3,2]) - M([1,1,1,1,1,1]) == M([0,1,2,3,2,1])
            True
        """
        cdef Matrix_gap cleft = <Matrix_gap> left
        cdef Matrix_gap ans = cleft._new(cleft._nrows, cleft._ncols)
        ans._libgap = left._libgap - (<Matrix_gap> right)._libgap
        return ans

    cdef Matrix _matrix_times_matrix_(left, Matrix right):
        r"""
        TESTS::

            sage: M = MatrixSpace(QQ, 2, implementation='gap')
            sage: m1 = M([1,2,-4,3])
            sage: m2 = M([-1,1,1,-1])
            sage: m1 * m2
            [ 1 -1]
            [ 7 -7]
        """
        if left._ncols != right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right.")
        cdef Matrix_gap M = left._new(left._nrows, right._ncols)
        M._libgap = <Matrix_gap> ((<Matrix_gap> left)._libgap * (<Matrix_gap> right)._libgap)
        return M

    def transpose(self):
        r"""
        Return the transpose of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(QQ, 2, implementation='gap')
            sage: M([4,2,23,52]).transpose()
            [ 4 23]
            [ 2 52]

            sage: M = MatrixSpace(QQ, 1, 3, implementation='gap')
            sage: M([4,2,52]).transpose()
            [ 4]
            [ 2]
            [52]
        """
        cdef Matrix_gap M
        M = self._new(self._ncols, self._nrows)
        M._libgap = self._libgap.TransposedMat()
        return M

    def determinant(self):
        r"""
        Return the determinant of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: M([2, 1, 1, 1]).determinant()
            1
            sage: M([2, 1, 3, 3]).determinant()
            3

        TESTS::

            sage: M = MatrixSpace(ZZ, 1, implementation='gap')
            sage: parent(M(1).determinant())
            Integer Ring

            sage: M = MatrixSpace(QQ, 1, implementation='gap')
            sage: parent(M(1).determinant())
            Rational Field

            sage: # needs sage.rings.number_field
            sage: M = MatrixSpace(UniversalCyclotomicField(), 1, implementation='gap')
            sage: parent(M(1).determinant())
            Universal Cyclotomic Field
        """
        return self._base_ring(self._libgap.DeterminantMat())

    def trace(self):
        r"""
        Return the trace of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: M([2, 1, 1, 1]).trace()
            3
            sage: M([2, 1, 3, 3]).trace()
            5

        TESTS::

            sage: M = MatrixSpace(ZZ, 1, implementation='gap')
            sage: parent(M(1).trace())
            Integer Ring

            sage: M = MatrixSpace(QQ, 1, implementation='gap')
            sage: parent(M(1).trace())
            Rational Field

            sage: # needs sage.rings.number_field
            sage: M = MatrixSpace(UniversalCyclotomicField(), 1, implementation='gap')
            sage: parent(M(1).trace())
            Universal Cyclotomic Field
        """
        return self._base_ring(self._libgap.TraceMat())

    def rank(self):
        r"""
        Return the rank of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: M([2, 1, 1, 1]).rank()
            2
            sage: M([2, 1, 4, 2]).rank()
            1
        """
        return int(self._libgap.RankMat())

    def minpoly(self, var='x', **kwds):
        """
        Compute the minimal polynomial.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, implementation='gap')
            sage: M([0, 1, -1, -1]).minpoly()
            x^2 + x + 1
        """
        po = self._libgap.MinimalPolynomial().sage()
        return po.change_variable_name(var)

    minimal_polynomial = minpoly

    def elementary_divisors(self):
        """
        Return the list of elementary divisors of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 5, implementation='gap')
            sage: T = M([(i+1)**(j+1) for i in range(5) for j in range(5)])
            sage: T.elementary_divisors()
            [1, 2, 6, 24, 120]
        """
        return [self._base_ring(u)
                for u in self._libgap.ElementaryDivisorsMat()]
