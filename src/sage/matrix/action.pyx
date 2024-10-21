# sage_setup: distribution = sagemath-modules
"""
Actions used by the coercion model for matrix and vector multiplications

.. WARNING::

    The class :class:`MatrixMulAction` and its descendants extends the class
    :class:`Action`. As a consequence objects from these classes only keep weak
    references to the underlying sets which are acted upon. This decision was
    made in :issue:`715` in order to allow garbage collection within the coercion
    framework, where actions are mainly used, and avoid memory leaks.

    To ensure that the underlying set of such an object does not get garbage
    collected, it is sufficient to explicitly create a strong reference to it
    before creating the action.

    ::

        sage: MSQ = MatrixSpace(QQ, 2)
        sage: MSZ = MatrixSpace(ZZ['x'], 2)
        sage: A = MSQ.get_action(MSZ)
        sage: A
        Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field
         on Full MatrixSpace of 2 by 2 dense matrices
          over Univariate Polynomial Ring in x over Integer Ring
        sage: import gc
        sage: _ = gc.collect()
        sage: A
        Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field
         on Full MatrixSpace of 2 by 2 dense matrices
          over Univariate Polynomial Ring in x over Integer Ring

.. NOTE::

    The :func:`MatrixSpace` function caches the objects it creates. Therefore,
    the underlying set ``MSZ`` in the above example will not be garbage
    collected, even if it is not strongly ref'ed. Nonetheless, there is no
    guarantee that the set that is acted upon will always be cached in such a
    way, so that following the above example is good practice.

EXAMPLES:

An action requires a common parent for the base rings, so the following
doesn't work (see :issue:`17859`)::

    sage: vector(QQ, [1]) * matrix(Zmod(2), [[1]])
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for *: 'Vector space of
    dimension 1 over Rational Field' and 'Full MatrixSpace of 1 by 1
    dense matrices over Ring of integers modulo 2'

AUTHOR:

- Robert Bradshaw (2007-09): Initial version.
"""

# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator

from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import FreeModule, FreeModule_generic
from sage.structure.coerce cimport coercion_model
from sage.categories.homset import Hom, End


try:
    from sage.schemes.generic.homset import SchemeHomset_generic, SchemeHomset_points
except ImportError:
    SchemeHomset_generic = SchemeHomset_points = ()


cdef class MatrixMulAction(Action):
    """
    Abstract base class for a matrix space acting on something.

    EXAMPLES::

        sage: MSQ = MatrixSpace(QQ, 2)
        sage: MSZ = MatrixSpace(ZZ['x'], 2)
        sage: A = MSQ.get_action(MSZ); A
        Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field
         on Full MatrixSpace of 2 by 2 dense matrices
            over Univariate Polynomial Ring in x over Integer Ring
        sage: A.actor()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: A.domain()
        Full MatrixSpace of 2 by 2 dense matrices
         over Univariate Polynomial Ring in x over Integer Ring
        sage: A.codomain()
        Full MatrixSpace of 2 by 2 dense matrices
         over Univariate Polynomial Ring in x over Rational Field
    """
    def __init__(self, G, S, is_left):
        if not isinstance(G, MatrixSpace):
            raise TypeError("Not a matrix space: %s" % G)
        if isinstance(S, SchemeHomset_generic):
            if G.base_ring() is not S.domain().base_ring():
                base = coercion_model.common_parent(G.base_ring(), S.domain().base_ring())
            else:
                base = G.base_ring()
        else:
            if G.base_ring() is not S.base_ring():
                base = coercion_model.common_parent(G.base_ring(), S.base_ring())
            else:
                base = G.base_ring()
            self.fix_sparseness = G.is_sparse() != S.is_sparse()
        Action.__init__(self, G, S, is_left, operator.mul)
        self._codomain = self._create_codomain(base)

    def codomain(self):
        return self._codomain


cdef class MatrixMatrixAction(MatrixMulAction):
    """
    Action of a matrix on another matrix.

    This is always implemented as a left action.

    EXAMPLES:

    By :issue:`715`, there only is a weak reference on the underlying set,
    so that it can be garbage collected if only the action itself is
    explicitly referred to. Hence, we first assign the involved matrix
    spaces to a variable::

        sage: R.<x> = ZZ[]
        sage: MSR = MatrixSpace(R, 3, 3)
        sage: MSQ = MatrixSpace(QQ, 3, 2)
        sage: from sage.matrix.action import MatrixMatrixAction
        sage: A = MatrixMatrixAction(MSR, MSQ); A
        Left action
         by Full MatrixSpace of 3 by 3 dense matrices
            over Univariate Polynomial Ring in x over Integer Ring
         on Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        sage: A.codomain()
        Full MatrixSpace of 3 by 2 dense matrices
         over Univariate Polynomial Ring in x over Rational Field
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
    def __init__(self, G, S):
        if not isinstance(S, MatrixSpace):
            raise TypeError("Not a matrix space: %s" % S)

        MatrixMulAction.__init__(self, G, S, True)

        # disallow multiplication on different backends (same size and rings)
        if (G.base_ring() is S.base_ring() and
           G.is_sparse() == S.is_sparse() and
           G.Element is not S.Element):
            raise TypeError("no matrix multiplication between different implementations")

        # disallow multiplication (sparse) x (dense) when the densification is not the default
        # implementation
        if self.fix_sparseness:
            if G.is_sparse():
                if not S._has_default_implementation():
                    raise TypeError("matrix multiplication not allowed")
            else:
                if not G._has_default_implementation():
                    raise TypeError("matrix multiplication not allowed")

    def _create_codomain(self, base):
        """
        EXAMPLES:

        By :issue:`715`, there only is a weak reference on the underlying set,
        so that it can be garbage collected if only the action itself is
        explicitly referred to. Hence, we first assign the involved matrix
        spaces to a variable::

            sage: from sage.matrix.action import MatrixMatrixAction
            sage: R.<x> = ZZ[]
            sage: MSR = MatrixSpace(R, 3, 3)
            sage: MSQ = MatrixSpace(QQ, 3, 2)
            sage: A = MatrixMatrixAction(MSR, MSQ); A
            Left action by Full MatrixSpace of 3 by 3 dense matrices over Univariate Polynomial Ring in x over Integer Ring
            on Full MatrixSpace of 3 by 2 dense matrices over Rational Field
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
                    (self.G.ncols(), self.underlying_set().nrows()))
        return MatrixSpace(base, self.G.nrows(), self.underlying_set().ncols(),
                           sparse = self.G.is_sparse() and self.underlying_set().is_sparse())

    cpdef _act_(self, g, s):
        """
        EXAMPLES:

        Respects compatible subdivisions::

            sage: M = matrix(5, 5, prime_range(100))                                    # needs sage.libs.pari
            sage: M.subdivide(2, 3); M                                                  # needs sage.libs.pari
            [ 2  3  5| 7 11]
            [13 17 19|23 29]
            [--------+-----]
            [31 37 41|43 47]
            [53 59 61|67 71]
            [73 79 83|89 97]
            sage: N = matrix(5, 2, [n^2 for n in range(10)])
            sage: N.subdivide(3, 1); N
            [ 0| 1]
            [ 4| 9]
            [16|25]
            [--+--]
            [36|49]
            [64|81]
            sage: M*N                                                                   # needs sage.libs.pari
            [ 1048| 1388]
            [ 3056| 4117]
            [-----+-----]
            [ 5360| 7303]
            [ 8168|11143]
            [11056|15077]

        Note that this is just like block matrix multiplication::

            sage: M.subdivision(0,0) * N.subdivision(0,0) + M.subdivision(0,1) * N.subdivision(1,0)                     # needs sage.libs.pari
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
            sage: M*N                                                                   # needs sage.libs.pari
            [ 1048  1388]
            [ 3056  4117]
            [ 5360  7303]
            [ 8168 11143]
            [11056 15077]
        """
        cdef Matrix A = <Matrix>g
        cdef Matrix B = <Matrix>s
        if A._parent._base is not self._codomain._base:
            A = A.change_ring(self._codomain._base)
        if B._parent._base is not self._codomain._base:
            B = B.change_ring(self._codomain._base)
        if self.fix_sparseness:
            if B.is_sparse_c():
                B = B.dense_matrix()
            else:
                A = A.dense_matrix()
        assert type(A) is type(B), (type(A), type(B))
        prod = A._matrix_times_matrix_(B)
        if A._subdivisions is not None or B._subdivisions is not None:
            Asubs = A.subdivisions()
            Bsubs = B.subdivisions()
            if Asubs[1] == Bsubs[0]:
                prod.subdivide(Asubs[0], Bsubs[1])
        return prod


cdef class MatrixVectorAction(MatrixMulAction):
    """
    Left action of a matrix on a vector
    """
    def __init__(self, G, S):
        """
        EXAMPLES::

            sage: from sage.matrix.action import MatrixVectorAction
            sage: A = MatrixVectorAction(MatrixSpace(QQ, 3, 3), VectorSpace(CDF, 4)); A
            Traceback (most recent call last):
            ...
            TypeError: incompatible dimensions 3, 4
            """
        if not isinstance(S, FreeModule_generic):
            raise TypeError("Not a free module: %s" % S)
        MatrixMulAction.__init__(self, G, S, True)

    def _create_codomain(self, base):
        """
        EXAMPLES::

            sage: from sage.matrix.action import MatrixVectorAction
            sage: M = MatrixSpace(QQ, 5, 3)
            sage: V = VectorSpace(CDF, 3)    # strong reference prevents garbage collection
            sage: A = MatrixVectorAction(M, V); A
            Left action by Full MatrixSpace of 5 by 3 dense matrices over Rational Field
             on Vector space of dimension 3 over Complex Double Field
            sage: A.codomain()
            Vector space of dimension 5 over Complex Double Field
        """
        if self.G.ncols() != self.underlying_set().degree():
            raise TypeError("incompatible dimensions %s, %s" % (self.G.ncols(),
                                                                 self.underlying_set().degree()))
        return FreeModule(base, self.G.nrows(), sparse = self.G.is_sparse())

    cpdef _act_(self, g, s):
        cdef Matrix A = <Matrix>g
        cdef Vector v = <Vector>s
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
    """
    Right action of a matrix on a vector
    """
    def __init__(self, G, S):
        """
        EXAMPLES::

            sage: from sage.matrix.action import VectorMatrixAction
            sage: A = VectorMatrixAction(MatrixSpace(QQ, 5, 3), VectorSpace(CDF, 3)); A
            Traceback (most recent call last):
            ...
            TypeError: incompatible dimensions 5, 3
        """
        if not isinstance(S, FreeModule_generic):
            raise TypeError("Not a free module: %s" % S)
        MatrixMulAction.__init__(self, G, S, False)

    def _create_codomain(self, base):
        """
        EXAMPLES::

            sage: from sage.matrix.action import VectorMatrixAction
            sage: M = MatrixSpace(QQ, 3, 5)
            sage: V = VectorSpace(CDF, 3)
            sage: A = VectorMatrixAction(M, V)
            sage: A
            Right action by Full MatrixSpace of 3 by 5 dense matrices over Rational Field
             on Vector space of dimension 3 over Complex Double Field
            sage: A.codomain()
            Vector space of dimension 5 over Complex Double Field
        """
        if self.G.nrows() != self.underlying_set().degree():
            raise TypeError("incompatible dimensions %s, %s" % (self.G.nrows(),
                                                                 self.underlying_set().degree()))
        return FreeModule(base, self.G.ncols(), sparse = self.G.is_sparse())

    cpdef _act_(self, g, s):
        cdef Matrix A = <Matrix>g
        cdef Vector v = <Vector>s
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

cdef class MatrixPolymapAction(MatrixMulAction):
    """
    Left action of a matrix on a scheme polynomial morphism
    """
    def __init__(self, G, S):
        """
        Initialize the action.

        EXAMPLES::

            sage: from sage.matrix.action import MatrixPolymapAction
            sage: M = MatrixSpace(QQ,2,2)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: A = MatrixPolymapAction(M,H)
            sage: A
            Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational
            Field on Set of morphisms
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 1 over Rational Field
        """
        if not isinstance(S, SchemeHomset_generic):
            raise TypeError("not a scheme polynomial morphism: %s"% S)
        MatrixMulAction.__init__(self, G, S, True)

    def _create_codomain(self, base):
        """
        EXAMPLES::

            sage: from sage.matrix.action import MatrixPolymapAction
            sage: M = MatrixSpace(QQ,2,2)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: A = MatrixPolymapAction(M,H)
            sage: A.codomain()
            Set of morphisms
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 1 over Rational Field
            sage: A.codomain().is_endomorphism_set()
            True
        """
        if self.underlying_set().is_endomorphism_set():
            return End(self.underlying_set().domain().change_ring(base))
        return Hom(self.underlying_set().domain().change_ring(base), self.underlying_set().codomain().change_ring(base))

    cpdef _act_(self, mat, f):
        """
        Call the action.

        INPUT:

        - ``mat`` -- a matrix

        - ``f`` -- a scheme homomorphism

        EXAMPLES::

            sage: from sage.matrix.action import MatrixPolymapAction
            sage: M = MatrixSpace(QQ, 2, 2)
            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2])
            sage: A = MatrixPolymapAction(M, H)
            sage: m = matrix([[1,1], [0,1]])
            sage: A._act_(m, f)
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 2*y^2 : y^2)
        """
        return f._matrix_times_polymap_(mat, self._codomain)

cdef class PolymapMatrixAction(MatrixMulAction):
    """
    Right action of a matrix on a scheme polynomial morphism
    """
    def __init__(self, G, S):
        """
        Initialize the action.

        EXAMPLES::

            sage: from sage.matrix.action import PolymapMatrixAction
            sage: M = MatrixSpace(QQ,2,2)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: A = PolymapMatrixAction(M,H)
            sage: A
            Right action by Full MatrixSpace of 2 by 2 dense matrices over Rational
            Field on Set of morphisms
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 1 over Rational Field
        """
        if not isinstance(S, SchemeHomset_generic):
            raise TypeError("not a scheme polynomial morphism: %s"% S)
        MatrixMulAction.__init__(self, G, S, False  )

    def _create_codomain(self, base):
        """
        Create the codomain.

        EXAMPLES::

            sage: from sage.matrix.action import PolymapMatrixAction
            sage: M = MatrixSpace(QQ,2,2)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: A = PolymapMatrixAction(M,H)
            sage: A.codomain()
            Set of morphisms
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 1 over Rational Field
            sage: A.codomain().is_endomorphism_set()
            True
        """
        if self.underlying_set().is_endomorphism_set():
            return End(self.underlying_set().domain().change_ring(base))
        return Hom(self.underlying_set().domain().change_ring(base), self.underlying_set().codomain().change_ring(base))

    cpdef _act_(self, mat, f):
        """
        Call the action.

        INPUT:

        - ``mat`` -- a matrix

        - ``f`` -- a scheme homomorphism

        EXAMPLES::

            sage: from sage.matrix.action import PolymapMatrixAction
            sage: M = MatrixSpace(QQ, 2, 2)
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2])
            sage: A = PolymapMatrixAction(M, H)
            sage: m = matrix([[1,1], [0,1]])
            sage: A._act_(m, f)
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 2*x*y + 2*y^2 : y^2)
        """
        return f._polymap_times_matrix_(mat, self._codomain)


cdef class MatrixSchemePointAction(MatrixMulAction):
    r"""
    Action class for left multiplication of schemes points by matrices.
    """
    def __init__(self, G, S):
        """
        Initialization of action class.

        EXAMPLES::

            sage: from sage.matrix.action import MatrixSchemePointAction
            sage: M = MatrixSpace(QQ, 2, 2)
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: A = MatrixSchemePointAction(M, P(QQ))
            sage: A
            Left action by Full MatrixSpace of 2 by 2 dense matrices over
            Rational Field on Set of rational points of Projective Space
            of dimension 1 over Rational Field
        """
        if not isinstance(S, SchemeHomset_points):
            raise TypeError("not a homset of scheme points: %s"% S)
        MatrixMulAction.__init__(self, G, S, True)

    def _create_codomain(self, base):
        """
        Create the point homset for the resulting point.

        EXAMPLES::

            sage: from sage.matrix.action import MatrixSchemePointAction
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: M = MatrixSpace(QQ, 2, 2)
            sage: A = MatrixSchemePointAction(M, P(QQ))
            sage: A.codomain()
            Set of rational points of Projective Space of dimension 1 over Rational Field
        """
        #need to extend the base of the ambient space
        #and return the set of point over the base
        amb = self.underlying_set().codomain()
        return amb.change_ring(base)(base)

    cpdef _act_(self, mat, P):
        """
        Action of matrices on scheme points.

        EXAMPLES::

            sage: from sage.matrix.action import MatrixSchemePointAction
            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: Q = P(1,1)
            sage: M = MatrixSpace(QQ, 2, 2)
            sage: A = MatrixSchemePointAction(M, Q.parent())
            sage: m = matrix([[1,1], [0,1]])
            sage: A._act_(m, Q)
            (2 : 1)
        """
        return P._matrix_times_point_(mat, self._codomain)
