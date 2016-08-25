"""
Actions used by the coercion model for matrix and vector multiplications

.. WARNING::

    The class :class:`MatrixMulAction` and its descendants extends the class
    :class:`Action`. As a consequence objects from these classes only keep weak
    references to the underlying sets which are acted upon. This decision was
    made in :trac:`715` in order to allow garbage collection within the coercion
    framework, where actions are mainly used, and avoid memory leaks.

    To ensure that the underlying set of such an object does not get garbage
    collected, it is sufficient to explicitly create a strong reference to it
    before creating the action.

    ::

        sage: MSQ = MatrixSpace(QQ, 2)
        sage: MSZ = MatrixSpace(ZZ['x'], 2)
        sage: A = MSQ.get_action(MSZ)
        sage: A
        Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field on Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring
        sage: import gc
        sage: _ = gc.collect()
        sage: A
        Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field on Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring

.. NOTE::

    The :func:`MatrixSpace` function caches the objects it creates. Therefore,
    the underlying set ``MSZ`` in the above example will not be garbage
    collected, even if it is not strongly ref'ed. Nonetheless, there is no
    guarantee that the set that is acted upon will always be cached in such a
    way, so that following the above example is good practice.

EXAMPLES:

An action requires a common parent for the base rings, so the following
doesn't work (see :trac:`17859`)::

    sage: vector(QQ, [1]) * matrix(Zmod(2), [[1]])
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Vector space of
    dimension 1 over Rational Field' and 'Full MatrixSpace of 1 by 1
    dense matrices over Ring of integers modulo 2'

AUTHOR:

- Robert Bradshaw (2007-09): Initial version.
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import operator

from matrix_space import MatrixSpace, is_MatrixSpace
from sage.modules.free_module import FreeModule, is_FreeModule
from sage.structure.element cimport coercion_model


cdef class MatrixMulAction(Action):
    def __init__(self, G, S, is_left):
        if not is_MatrixSpace(G):
            raise TypeError("Not a matrix space: %s" % G)
        if G.base_ring() is not S.base_ring():
            base = coercion_model.common_parent(G.base_ring(), S.base_ring())
        else:
            base = G.base_ring()
        Action.__init__(self, G, S, is_left, operator.mul)
        self._codomain = self._create_codomain(base)
        self.fix_sparseness = G.is_sparse() != S.is_sparse()

    def codomain(self):
        return self._codomain

    def domain(self):
        """
        EXAMPLES:

        By :trac:`715`, there only is a weak reference on the underlying set,
        so that it can be garbage collected if only the action itself is
        explicitly referred to. Hence, we first assign the involved matrix
        spaces to a variable::

            sage: MSQ = MatrixSpace(QQ, 2)
            sage: MSZ = MatrixSpace(ZZ['x'], 2)
            sage: A = MSQ.get_action(MSZ); A
            Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field on Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring
            sage: A.actor()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: A.domain()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring
            sage: A.codomain()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

        .. NOTE::

            The :func:`MatrixSpace` function caches the object it creates.
            Therefore, the underlying set ``MSZ`` in the above example will not
            be garbage collected, even if it is not strongly ref'ed.
            Nonetheless, there is no guarantee that the set that is acted upon
            will always be cached in such a way, so that following the above
            example is good practice.

        """
        return self.underlying_set()


cdef class MatrixMatrixAction(MatrixMulAction):
    def __init__(self, G, S):
        """
        EXAMPLES:

        By :trac:`715`, there only is a weak reference on the underlying set,
        so that it can be garbage collected if only the action itself is
        explicitly referred to. Hence, we first assign the involved matrix
        spaces to a variable::

            sage: R.<x> = ZZ[]
            sage: MSR = MatrixSpace(R, 3, 3)
            sage: MSQ = MatrixSpace(QQ, 3, 2)
            sage: from sage.matrix.action import MatrixMatrixAction
            sage: A = MatrixMatrixAction(MSR, MSQ); A
            Left action by Full MatrixSpace of 3 by 3 dense matrices over Univariate Polynomial Ring in x over Integer Ring on Full MatrixSpace of 3 by 2 dense matrices over Rational Field
            sage: A.codomain()
            Full MatrixSpace of 3 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field
            sage: A(matrix(R, 3, 3, x), matrix(QQ, 3, 2, range(6)))
            [  0   x]
            [2*x 3*x]
            [4*x 5*x]

        .. NOTE::

            The :func:`MatrixSpace` function caches the object it creates.
            Therefore, the underlying set ``MSZ`` in the above example will not
            be garbage collected, even if it is not strongly ref'ed.
            Nonetheless, there is no guarantee that the set that is acted upon
            will always be cached in such a way, so that following the above
            example is good practice.

        """
        if not is_MatrixSpace(S):
            raise TypeError("Not a matrix space: %s" % S)
        MatrixMulAction.__init__(self, G, S, True)

    def _create_codomain(self, base):
        """
        EXAMPLES:

        By :trac:`715`, there only is a weak reference on the underlying set,
        so that it can be garbage collected if only the action itself is
        explicitly referred to. Hence, we first assign the involved matrix
        spaces to a variable::

            sage: from sage.matrix.action import MatrixMatrixAction
            sage: R.<x> = ZZ[]
            sage: MSR = MatrixSpace(R, 3, 3)
            sage: MSQ = MatrixSpace(QQ, 3, 2)
            sage: A = MatrixMatrixAction(MSR, MSQ); A
            Left action by Full MatrixSpace of 3 by 3 dense matrices over Univariate Polynomial Ring in x over Integer Ring on Full MatrixSpace of 3 by 2 dense matrices over Rational Field
            sage: A.codomain()
            Full MatrixSpace of 3 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

        .. NOTE::

            The :func:`MatrixSpace` function caches the object it creates.
            Therefore, the underlying set ``MSZ`` in the above example will not
            be garbage collected, even if it is not strongly ref'ed.
            Nonetheless, there is no guarantee that the set that is acted upon
            will always be cached in such a way, so that following the above
            example is good practice.

        """
        if self.G.ncols() != self.underlying_set().nrows():
            raise TypeError("incompatible dimensions %s, %s" %
                    (self.G.ncols(),  self.underlying_set().nrows()))
        return MatrixSpace(base, self.G.nrows(), self.underlying_set().ncols(),
                           sparse = self.G.is_sparse() and self.underlying_set().is_sparse())

    cpdef _call_(self, g, s):
        """
        EXAMPLES:

        Respects compatible subdivisions::

            sage: M = matrix(5, 5, prime_range(100))
            sage: M.subdivide(2,3); M
            [ 2  3  5| 7 11]
            [13 17 19|23 29]
            [--------+-----]
            [31 37 41|43 47]
            [53 59 61|67 71]
            [73 79 83|89 97]
            sage: N = matrix(5,2,[n^2 for n in range(10)])
            sage: N.subdivide(3,1); N
            [ 0| 1]
            [ 4| 9]
            [16|25]
            [--+--]
            [36|49]
            [64|81]
            sage: M*N
            [ 1048| 1388]
            [ 3056| 4117]
            [-----+-----]
            [ 5360| 7303]
            [ 8168|11143]
            [11056|15077]

        Note that this is just like block matrix multiplication::

            sage: M.subdivision(0,0) * N.subdivision(0,0) + M.subdivision(0,1) * N.subdivision(1,0)
            [1048]
            [3056]

        If the subdivisions aren't compatible, ignore them.
        ::

            sage: N.subdivide(1,1); N
            [ 0| 1]
            [--+--]
            [ 4| 9]
            [16|25]
            [36|49]
            [64|81]
            sage: M*N
            [ 1048  1388]
            [ 3056  4117]
            [ 5360  7303]
            [ 8168 11143]
            [11056 15077]

        """
        cdef Matrix A = g #<Matrix>g
        cdef Matrix B = s #<Matrix>s
        if A._parent._base is not self._codomain._base:
            A = A.change_ring(self._codomain._base)
        if B._parent._base is not self._codomain._base:
            B = B.change_ring(self._codomain._base)
        if self.fix_sparseness:
            if B.is_sparse_c():
                B = B.dense_matrix()
            else:
                A = A.dense_matrix()
        prod = A._matrix_times_matrix_(B)
        if A._subdivisions is not None or B._subdivisions is not None:
            Asubs = A.subdivisions()
            Bsubs = B.subdivisions()
            if Asubs[1] == Bsubs[0]:
                prod.subdivide(Asubs[0], Bsubs[1])
        return prod


cdef class MatrixVectorAction(MatrixMulAction):
    def __init__(self, G, S):
        """
        EXAMPLES::

            sage: from sage.matrix.action import MatrixVectorAction
            sage: A = MatrixVectorAction(MatrixSpace(QQ, 3, 3), VectorSpace(CDF, 4)); A
            Traceback (most recent call last):
            ...
            TypeError: incompatible dimensions 3, 4
            """
        if not is_FreeModule(S):
            raise TypeError("Not a free module: %s" % S)
        MatrixMulAction.__init__(self, G, S, True)

    def _create_codomain(self, base):
        """
        EXAMPLES::

            sage: from sage.matrix.action import MatrixVectorAction
            sage: M = MatrixSpace(QQ, 5, 3)
            sage: V = VectorSpace(CDF, 3)    # strong reference prevents garbage collection
            sage: A = MatrixVectorAction(M, V); A
            Left action by Full MatrixSpace of 5 by 3 dense matrices over Rational Field on Vector space of dimension 3 over Complex Double Field
            sage: A.codomain()
            Vector space of dimension 5 over Complex Double Field
        """
        if self.G.ncols() != self.underlying_set().degree():
            raise TypeError("incompatible dimensions %s, %s" % (self.G.ncols(),
                                                                 self.underlying_set().degree()))
        return FreeModule(base, self.G.nrows(), sparse = self.G.is_sparse())

    cpdef _call_(self, g, s):
        cdef Matrix A = g #<Matrix>g
        cdef Vector v = s #<Vector>s
        if A._parent._base is not self._codomain._base:
            A = A.change_ring(self._codomain._base)
        if v._parent._base is not self._codomain._base:
            v = v.change_ring(self._codomain._base)
        if self.fix_sparseness:
            if A.is_sparse_c():
                v = v.sparse_vector()
            else:
                v = v.dense_vector()
        return A._matrix_times_vector_(v)


cdef class VectorMatrixAction(MatrixMulAction):
    def __init__(self, G, S):
        """
        EXAMPLES::

            sage: from sage.matrix.action import VectorMatrixAction
            sage: A = VectorMatrixAction(MatrixSpace(QQ, 5, 3), VectorSpace(CDF, 3)); A
            Traceback (most recent call last):
            ...
            TypeError: incompatible dimensions 5, 3
        """
        if not is_FreeModule(S):
            raise TypeError("Not a free module: %s" % S)
        MatrixMulAction.__init__(self, G, S, False)

    def _create_codomain(self, base):
        """
        EXAMPLES::

            sage: from sage.matrix.action import VectorMatrixAction
            sage: A = VectorMatrixAction(MatrixSpace(QQ, 3, 5), VectorSpace(CDF, 3)); A
            Right action by Full MatrixSpace of 3 by 5 dense matrices over Rational Field on Vector space of dimension 3 over Complex Double Field
            sage: A.codomain()
            Vector space of dimension 5 over Complex Double Field
        """
        if self.G.nrows() != self.underlying_set().degree():
            raise TypeError("incompatible dimensions %s, %s" % (self.G.nrows(),
                                                                 self.underlying_set().degree()))
        return FreeModule(base, self.G.ncols(), sparse = self.G.is_sparse())

    cpdef _call_(self, s, g):
        cdef Matrix A = g #<Matrix>g
        cdef Vector v = s #<Vector>s
        if A._parent._base is not self._codomain._base:
            A = A.change_ring(self._codomain._base)
        if v._parent._base is not self._codomain._base:
            v = v.change_ring(self._codomain._base)
        if self.fix_sparseness:
            if A.is_sparse_c():
                v = v.sparse_vector()
            else:
                v = v.dense_vector()
        return (<Matrix>A)._vector_times_matrix_(v) # v * A

