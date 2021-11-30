"""
Matrix Group Elements

EXAMPLES::

    sage: F = GF(3); MS = MatrixSpace(F,2,2)
    sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
    sage: G = MatrixGroup(gens); G
    Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )
    sage: g = G([[1,1],[0,1]])
    sage: h = G([[1,2],[0,1]])
    sage: g*h
    [1 0]
    [0 1]

You cannot add two matrices, since this is not a group operation.
You can coerce matrices back to the matrix space and add them
there::

    sage: g + h
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +:
    'Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )' and
    'Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )'

    sage: g.matrix() + h.matrix()
    [2 0]
    [0 2]

Similarly, you cannot multiply group elements by scalars but you can
do it with the underlying matrices::

    sage: 2*g
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for *: 'Integer Ring' and 'Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )'

AUTHORS:

- David Joyner (2006-05): initial version David Joyner

- David Joyner (2006-05): various modifications to address William
  Stein's TODO's.

- William Stein (2006-12-09): many revisions.

- Volker Braun (2013-1) port to new Parent, libGAP.

- Travis Scrimshaw (2016-01): reworks class hierarchy in order
  to cythonize
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element cimport MultiplicativeGroupElement, Element, MonoidElement, Matrix
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport richcmp
from sage.libs.gap.element cimport GapElement, GapElement_List
from sage.groups.libgap_wrapper cimport ElementLibGAP

from sage.structure.element import is_Matrix
from sage.structure.factorization import Factorization
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ


cpdef is_MatrixGroupElement(x):
    """
    Test whether ``x`` is a matrix group element

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.group_element import is_MatrixGroupElement
        sage: is_MatrixGroupElement('helloooo')
        False

        sage: G = GL(2,3)
        sage: is_MatrixGroupElement(G.an_element())
        True
    """
    return isinstance(x, (MatrixGroupElement_generic, MatrixGroupElement_gap))

###################################################################
#
# Matrix group elements implemented in Sage
#
###################################################################

cdef class MatrixGroupElement_generic(MultiplicativeGroupElement):
    """
    Element of a matrix group over a generic ring.

    The group elements are implemented as Sage matrices.

    INPUT:

    - ``M`` -- a matrix

    - ``parent`` -- the parent

    - ``check`` -- bool (default: ``True``); if ``True``, then
       does some type checking

    - ``convert`` -- bool (default: ``True``); if ``True``, then
      convert ``M`` to the right matrix space

    EXAMPLES::

        sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
        sage: g = W.an_element()
        sage: g
        [ 0  0 -1]
        [ 1  0 -1]
        [ 0  1 -1]
    """
    def __init__(self, parent, M, check=True, convert=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: TestSuite(g).run()
        """
        if convert:
            M = parent.matrix_space()(M)
        if check:
            if not is_Matrix(M):
                raise TypeError('M must be a matrix')
            if M.parent() is not parent.matrix_space():
                raise TypeError('M must be a in the matrix space of the group')
            parent._check_matrix(M)
        super(MatrixGroupElement_generic, self).__init__(parent)
        if M.is_immutable():
            self._matrix = M
        else:
            self._matrix = M.__copy__()
            self._matrix.set_immutable()

    def __hash__(self):
        r"""
        TESTS::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: hash(g)
            660522311176098153  # 64-bit
            -606138007          # 32-bit
        """
        return hash(self._matrix)

    def __reduce__(self):
        """
        Implement pickling.

        TESTS::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: loads(g.dumps()) == g
            True
        """
        return (_unpickle_generic_element, (self.parent(), self._matrix,))

    def _repr_(self):
        """
        Return string representation of this matrix.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: W.an_element()
            [ 0  0 -1]
            [ 1  0 -1]
            [ 0  1 -1]
        """
        return str(self._matrix)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: latex(g)
            \left(\begin{array}{rrr}
            0 & 0 & -1 \\
            1 & 0 & -1 \\
            0 & 1 & -1
            \end{array}\right)
        """
        return self._matrix._latex_()

    cpdef _act_on_(self, x, bint self_on_left):
        """
        EXAMPLES::

            sage: W = CoxeterGroup(['A',4], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: g * vector([1,1,1,1])
            (0, 1, 1, 1)
            sage: v = vector([3,2,1,-1])
            sage: g = W.gen(1)
            sage: v * g == v * g.matrix()   # indirect doctest
            True
        """
        if not is_MatrixGroupElement(x) and x not in self.parent().base_ring():
            try:
                if self_on_left:
                    return self._matrix * x
                else:
                    return x * self._matrix
            except TypeError:
                return None

    cpdef _richcmp_(self, other, int op):
        """
        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: TestSuite(g).run()
            sage: h = W.gen(0) * W.gen(1) * W.gen(2)
            sage: g == h
            True
            sage: a = W.gen(0)
            sage: a == g
            False
            sage: a != g
            True
        """
        cdef MatrixGroupElement_generic x = <MatrixGroupElement_generic>self
        cdef MatrixGroupElement_generic y = <MatrixGroupElement_generic>other
        return richcmp(x._matrix, y._matrix, op)

    cpdef list list(self):
        """
        Return list representation of this matrix.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: g
            [-1  1  0]
            [ 0  1  0]
            [ 0  0  1]
            sage: g.list()
            [[-1, 1, 0], [0, 1, 0], [0, 0, 1]]
        """
        return [r.list() for r in self._matrix.rows()]

    def matrix(self):
        """
        Obtain the usual matrix (as an element of a matrix space)
        associated to this matrix group element.

        One reason to compute the associated matrix is that matrices
        support a huge range of functionality.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: g.matrix()
            [-1  1  0]
            [ 0  1  0]
            [ 0  0  1]
            sage: parent(g.matrix())
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

        Matrices have extra functionality that matrix group elements
        do not have::

            sage: g.matrix().charpoly('t')
            t^3 - t^2 - t + 1
        """
        return self._matrix

    def _matrix_(self, base=None):
        """
        Method used by the :func:`matrix` constructor.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: matrix(RDF, g)
            [-1.0  1.0  0.0]
            [ 0.0  1.0  0.0]
            [ 0.0  0.0  1.0]
        """
        return self.matrix()

    cpdef _mul_(self, other):
        """
        Return the product of ``self`` and`` other``, which must
        have identical parents.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: h = W.an_element()
            sage: g * h
            [ 1  0  0]
            [ 1  0 -1]
            [ 0  1 -1]
        """
        cdef Parent parent = self.parent()
        cdef MatrixGroupElement_generic y = <MatrixGroupElement_generic>other
        cdef Matrix M = self._matrix * y._matrix
        # Make it immutable so the constructor doesn't make a copy
        M.set_immutable()
        return parent.element_class(parent, M, check=False, convert=False)

    def is_one(self):
        """
        Return whether ``self`` is the identity of the group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3])
            sage: g = W.gen(0)
            sage: g.is_one()
            False

            sage: W.an_element().is_one()
            False
            sage: W.one().is_one()
            True
        """
        return self._matrix.is_one()

    def __invert__(self):
        """
        Return the inverse group element

        OUTPUT:

        A matrix group element.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: ~g
            [-1  1  0]
            [-1  0  1]
            [-1  0  0]
            sage: g * ~g == W.one()
            True
            sage: ~g * g == W.one()
            True

            sage: W = CoxeterGroup(['B',3])
            sage: W.base_ring()
            Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?
            sage: g = W.an_element()
            sage: ~g
            [-1  1  0]
            [-1  0  a]
            [-a  0  1]
        """
        cdef Parent parent = self.parent()
        cdef Matrix M = self._matrix
        # We have a special method for dense matrices over ZZ
        if M.base_ring() is ZZ and M.is_dense():
            M = M.inverse_of_unit()
        else:
            M = ~M
            if M.base_ring() is not parent.base_ring():
                M = M.change_ring(parent.base_ring())
        # Make it immutable so the constructor doesn't make a copy
        M.set_immutable()
        return parent.element_class(parent, M, check=False, convert=False)

    inverse = __invert__

###################################################################
#
# Matrix group elements implemented in GAP
#
###################################################################

cdef class MatrixGroupElement_gap(ElementLibGAP):
    """
    Element of a matrix group over a generic ring.

    The group elements are implemented as wrappers around libGAP matrices.

    INPUT:

    - ``M`` -- a matrix

    - ``parent`` -- the parent

    - ``check`` -- bool (default: ``True``); if ``True`` does some
      type checking

    - ``convert`` -- bool (default: ``True``); if ``True`` convert
      ``M`` to the right matrix space
    """
    def __init__(self, parent, M, check=True, convert=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: MS = MatrixSpace(GF(3),2,2)
            sage: G = MatrixGroup(MS([[1,0],[0,1]]), MS([[1,1],[0,1]]))
            sage: G.gen(0)
            [1 0]
            [0 1]
            sage: g = G.random_element()
            sage: TestSuite(g).run()
        """
        if isinstance(M, GapElement):
            ElementLibGAP.__init__(self, parent, M)
            return
        if convert:
            M = parent.matrix_space()(M)
        from sage.libs.gap.libgap import libgap
        M_gap = libgap(M)
        if check:
            if not is_Matrix(M):
                raise TypeError('M must be a matrix')
            if M.parent() is not parent.matrix_space():
                raise TypeError('M must be a in the matrix space of the group')
            parent._check_matrix(M, M_gap)
        ElementLibGAP.__init__(self, parent, M_gap)

    def __reduce__(self):
        """
        Implement pickling.

        TESTS::

            sage: MS = MatrixSpace(GF(3), 2, 2)
            sage: G = MatrixGroup(MS([[1,0],[0,1]]), MS([[1,1],[0,1]]))
            sage: loads(G.gen(0).dumps())
            [1 0]
            [0 1]
        """
        return (self.parent(), (self.matrix(),))

    def __hash__(self):
        r"""
        TESTS::

            sage: MS = MatrixSpace(GF(3), 2)
            sage: G = MatrixGroup([MS([1,1,0,1]), MS([1,0,1,1])])
            sage: g = G.an_element()
            sage: hash(g)
            -5306160029685893860  # 64-bit
            -181258980            # 32-bit
        """
        return hash(self.matrix())

    def _repr_(self):
        r"""
        Return string representation of this matrix.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: g  # indirect doctest
            [1 1]
            [0 1]
            sage: g._repr_()
            '[1 1]\n[0 1]'
        """
        return str(self.matrix())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: print(g._latex_())
            \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right)

        Type ``view(g._latex_())`` to see the object in an
        xdvi window (assuming you have latex and xdvi installed).
        """
        return self.matrix()._latex_()

    cpdef _act_on_(self, x, bint self_on_left):
        """
        EXAMPLES::

            sage: G = GL(4,7)
            sage: G.0 * vector([1,2,3,4])
            (3, 2, 3, 4)
            sage: v = vector(GF(7), [3,2,1,-1])
            sage: g = G.1
            sage: v * g == v * g.matrix()   # indirect doctest
            True
        """
        if not is_MatrixGroupElement(x) and x not in self.parent().base_ring():
            try:
                if self_on_left:
                    return self.matrix() * x
                else:
                    return x * self.matrix()
            except TypeError:
                return None

    cpdef _richcmp_(self, other, int op):
        """
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = G([1,1, 0,1])
            sage: g == h
            True
            sage: g == G.one()
            False
        """
        return richcmp(self.matrix(), other.matrix(), op)

    @cached_method
    def matrix(self):
        """
        Obtain the usual matrix (as an element of a matrix space)
        associated to this matrix group element.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: m = G.gen(0).matrix(); m
            [1 0]
            [0 1]
            sage: m.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 3

            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])])
            sage: g = G.0
            sage: g.matrix()
            [1 1]
            [0 1]
            sage: parent(g.matrix())
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7

        Matrices have extra functionality that matrix group elements
        do not have::

            sage: g.matrix().charpoly('t')
            t^2 + 5*t + 1
        """
        # We do a slightly specialized version of sage.libs.gap.element.GapElement.matrix()
        #   in order to use our current matrix space directly and avoid
        #   some overhead safety checks.
        entries = self.gap().Flat()
        MS = self.parent().matrix_space()
        ring = MS.base_ring()
        m = MS([x.sage(ring=ring) for x in entries])
        m.set_immutable()
        return m

    def _matrix_(self, base=None):
        """
        Method used by the :func:`matrix` constructor.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: g = G.gen(0)
            sage: M = matrix(GF(9), g); M; parent(M)
            [1 1]
            [0 1]
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field in z2 of size 3^2
        """
        return self.matrix()

    cpdef list list(self):
        """
        Return list representation of this matrix.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G.0
            sage: g
            [1 0]
            [0 1]
            sage: g.list()
            [[1, 0], [0, 1]]
        """
        return [r.list() for r in self.matrix().rows()]

    @cached_method
    def multiplicative_order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer `n` such that `g^n = 1`, or
        +Infinity if no such integer exists.

        EXAMPLES::

            sage: k = GF(7)
            sage: G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])]); G
            Matrix group over Finite Field of size 7 with 2 generators (
            [1 1]  [1 0]
            [0 1], [0 2]
            )
            sage: G.order()
            21
            sage: G.gen(0).multiplicative_order(), G.gen(1).multiplicative_order()
            (7, 3)

        ``order`` is just an alias for ``multiplicative_order``::

            sage: G.gen(0).order(), G.gen(1).order()
            (7, 3)

            sage: k = QQ
            sage: G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])]); G
            Matrix group over Rational Field with 2 generators (
            [1 1]  [1 0]
            [0 1], [0 2]
            )
            sage: G.order()
            +Infinity
            sage: G.gen(0).order(), G.gen(1).order()
            (+Infinity, +Infinity)

            sage: gl = GL(2, ZZ);  gl
            General Linear Group of degree 2 over Integer Ring
            sage: g = gl.gen(2);  g
            [1 1]
            [0 1]
            sage: g.order()
            +Infinity
        """
        order = self.gap().Order()
        if order.IsInt():
            return order.sage()
        else:
            assert order.IsInfinity()
            from sage.rings.all import Infinity
            return Infinity

    def word_problem(self, gens=None):
        r"""
        Solve the word problem.

        This method writes the group element as a product of the
        elements of the list ``gens``, or the standard generators of
        the parent of self if ``gens`` is None.

        INPUT:

        - ``gens`` -- a list/tuple/iterable of elements (or objects
          that can be converted to group elements), or ``None``
          (default). By default, the generators of the parent group
          are used.

        OUTPUT:

        A factorization object that contains information about the
        order of factors and the exponents. A ``ValueError`` is raised
        if the group element cannot be written as a word in ``gens``.

        ALGORITHM:

        Use GAP, which has optimized algorithms for solving the word
        problem (the GAP functions ``EpimorphismFromFreeGroup`` and
        ``PreImagesRepresentative``).

        EXAMPLES::

            sage: G = GL(2,5); G
            General Linear Group of degree 2 over Finite Field of size 5
            sage: G.gens()
            (
            [2 0]  [4 1]
            [0 1], [4 0]
            )
            sage: G(1).word_problem([G.gen(0)])
            1
            sage: type(_)
            <class 'sage.structure.factorization.Factorization'>

            sage: g = G([0,4,1,4])
            sage: g.word_problem()
            ([4 1]
             [4 0])^-1

        Next we construct a more complicated element of the group from the
        generators::

            sage: s,t = G.0, G.1
            sage: a = (s * t * s); b = a.word_problem(); b
            ([2 0]
             [0 1]) *
            ([4 1]
             [4 0]) *
            ([2 0]
             [0 1])
            sage: flatten(b)
            [
            [2 0]     [4 1]     [2 0]
            [0 1], 1, [4 0], 1, [0 1], 1
            ]
            sage: b.prod() == a
            True

        We solve the word problem using some different generators::

            sage: s = G([2,0,0,1]); t = G([1,1,0,1]); u = G([0,-1,1,0])
            sage: a.word_problem([s,t,u])
            ([2 0]
             [0 1])^-1 *
            ([1 1]
             [0 1])^-1 *
            ([0 4]
             [1 0]) *
            ([2 0]
             [0 1])^-1

        We try some elements that don't actually generate the group::

            sage: a.word_problem([t,u])
            Traceback (most recent call last):
            ...
            ValueError: word problem has no solution

        AUTHORS:

        - David Joyner and William Stein
        - David Loeffler (2010): fixed some bugs
        - Volker Braun (2013): LibGAP
        """
        from sage.libs.gap.libgap import libgap
        G = self.parent()
        if gens:
            gen = lambda i:gens[i]
            H = libgap.Group([G(x).gap() for x in gens])
        else:
            gen = G.gen
            H = G.gap()
        hom = H.EpimorphismFromFreeGroup()
        preimg = hom.PreImagesRepresentative(self.gap())

        if preimg.is_bool():
            assert preimg == libgap.eval('fail')
            raise ValueError('word problem has no solution')

        result = []
        n = preimg.NumberSyllables().sage()
        exponent_syllable  = libgap.eval('ExponentSyllable')
        generator_syllable = libgap.eval('GeneratorSyllable')
        for i in range(n):
            exponent  = exponent_syllable(preimg, i+1).sage()
            generator = gen(generator_syllable(preimg, i+1).sage() - 1)
            result.append( (generator, exponent) )
        result = Factorization(result)
        result._set_cr(True)
        return result

def _unpickle_generic_element(G, mat):
    """
    Unpickle the element in ``G`` given by ``mat``.

    EXAMPLES::

        sage: m1 = matrix(SR, [[1,2],[3,4]])
        sage: m2 = matrix(SR, [[1,3],[-1,0]])
        sage: G = MatrixGroup(m1, m2)
        sage: m = G.an_element()
        sage: from sage.groups.matrix_gps.group_element import _unpickle_generic_element
        sage: _unpickle_generic_element(G, m.matrix()) == m
        True
    """
    return G.element_class(G, mat, False, False)

