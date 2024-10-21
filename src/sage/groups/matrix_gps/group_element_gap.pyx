# sage_setup: distribution = sagemath-gap

r"""
Matrix group elements implemented in GAP
"""

#*****************************************************************************
#       Copyright (C) 2006      David Joyner and William Stein <wstein@gmail.com>
#                     2013      Volker Braun <vbraun.name@gmail.com>
#                     2015-2017 Vincent Delecroix
#                     2016      Travis Scrimshaw <tscrimsh at umn.edu>
#                     2018      Jeroen Demeyer
#                     2023      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.matrix_gps.group_element cimport is_MatrixGroupElement
from sage.libs.gap.element cimport GapElement
from sage.misc.cachefunc import cached_method
from sage.structure.element import Matrix
from sage.structure.factorization import Factorization
from sage.structure.richcmp cimport richcmp


cdef class MatrixGroupElement_gap(ElementLibGAP):
    """
    Element of a matrix group over a generic ring.

    The group elements are implemented as wrappers around libGAP matrices.

    INPUT:

    - ``M`` -- a matrix

    - ``parent`` -- the parent

    - ``check`` -- boolean (default: ``True``); if ``True``, do some
      type checking

    - ``convert`` -- boolean (default: ``True``); if ``True``, convert
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
            if not isinstance(M, Matrix):
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
        m = MS([ring(x) for x in entries])
        m.set_immutable()
        return m

    def _matrix_(self, base=None):
        """
        Method used by the :func:`matrix` constructor.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: g = G.gen(0)
            sage: M = matrix(GF(9), g); M; parent(M)                                    # needs sage.rings.finite_rings
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
            from sage.rings.infinity import Infinity
            return Infinity

    def word_problem(self, gens=None):
        r"""
        Solve the word problem.

        This method writes the group element as a product of the
        elements of the list ``gens``, or the standard generators of
        the parent of ``self`` if ``gens`` is ``None``.

        INPUT:

        - ``gens`` -- list/tuple/iterable of elements (or objects
          that can be converted to group elements), or ``None``
          (default); by default, the generators of the parent group
          are used

        OUTPUT:

        A factorization object that contains information about the
        order of factors and the exponents. A :exc:`ValueError` is raised
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
