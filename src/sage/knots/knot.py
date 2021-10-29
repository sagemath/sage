# -*- coding: utf-8 -*-
r"""
Knots

AUTHORS:

- Miguel Angel Marco Buzunariz
- Amit Jamadagni
"""

# ****************************************************************************
#       Copyright (C) 2014   Travis Scrimshaw <tscrim at ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.knots.link import Link
from sage.knots.knot_table import small_knots_table
from sage.knots.gauss_code import (recover_orientations, dowker_to_gauss,
                                   rectangular_diagram)

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.misc.fast_methods import Singleton
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.categories.monoids import Monoids

# We need Link to be first in the MRO in order to use its equality, hash, etc.
class Knot(Link, Element, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A knot.

    A knot is defined as embedding of the circle `\mathbb{S}^1` in the
    3-dimensional sphere `\mathbb{S}^3`, considered up to ambient isotopy.
    They represent the physical idea of a knotted rope, but with the
    particularity that the rope is closed. That is, the ends of the rope
    are joined.

    .. SEEALSO::

        :class:`Link`

    INPUT:

    - ``data`` -- see :class:`Link` for the allowable inputs
    - ``check`` -- optional, default ``True``. If ``True``, make sure
      that the data define a knot, not a link

    EXAMPLES:

    We construct the knot `8_{14}` and compute some invariants::

        sage: B = BraidGroup(4)
        sage: K = Knot(B([1,1,1,2,-1,2,-3,2,-3]))

    .. PLOT::
        :width: 300 px

        B = BraidGroup(4)
        K = Knot(B([1,1,1,2,-1,2,-3,2,-3]))
        sphinx_plot(K.plot())

    ::

        sage: K.alexander_polynomial()
        -2*t^-2 + 8*t^-1 - 11 + 8*t - 2*t^2
        sage: K.jones_polynomial()
        t^7 - 3*t^6 + 4*t^5 - 5*t^4 + 6*t^3 - 5*t^2 + 4*t + 1/t - 2
        sage: K.determinant()
        31
        sage: K.signature()
        -2
        sage: K.colored_jones_polynomial(2)  # long time
        q^-1 - 2 + 4*q - 5*q^2 + 6*q^3 - 5*q^4 + 4*q^5 - 3*q^6 + q^7

    REFERENCES:

    - :wikipedia:`Knot_(mathematics)`
    """
    @staticmethod
    def __classcall_private__(self, data, check=True):
        """
        Make sure this is an instance of the element class
        of :class:`Knots`.

        EXAMPLES::

            sage: B = BraidGroup(8)
            sage: K = Knot(B([-1, -1, -1, 2, 1, -2, 3, -2, 3]))
            sage: type(K)
            <class 'sage.knots.knot.Knots_with_category.element_class'>
        """
        return Knots().element_class(data, check=check)

    def __init__(self, data, check=True):
        """
        Initialize ``self``.

        TESTS::

            sage: B = BraidGroup(8)
            sage: K = Knot(B([1, -2, 1, -2]))
            sage: TestSuite(K).run()
            sage: K = Knot([[1, 1, 2, 2]])
            sage: TestSuite(K).run()

        The following is not a knot: it has two components. ::

            sage: Knot([[[1, 2], [-2, -1]], [1, -1]])
            Traceback (most recent call last):
            ...
            ValueError: the input has more than 1 connected component

            sage: Knot([[[1, 2], [-2, -1]], [1, -1]], check=False)
            Knot represented by 2 crossings
        """
        Element.__init__(self, Knots())
        Link.__init__(self, data)
        if check:
            if self.number_of_components() != 1:
                raise ValueError("the input has more than 1 connected "
                                 "component")

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: B = BraidGroup(8)
            sage: K = Knot(B([1, 2, 1, 2]))
            sage: K
            Knot represented by 4 crossings
            sage: K = Knot([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10],
            ....:           [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8],
            ....:           [12, 9, 13, 10]])
            sage: K
            Knot represented by 7 crossings
        """
        pd_len = len(self.pd_code())
        return 'Knot represented by {} crossings'.format(pd_len)

    def _unicode_art_(self):
        """
        Return unicode art for the knot.

        INPUT:

        - a knot

        OUTPUT:

        - unicode art for the knot

        EXAMPLES::

            sage: W = Knots()
            sage: K = W.from_dowker_code([-4,-6,-2])
            sage: unicode_art(K)
             ╭─╮
            ╭──│╮
            │╰╮││
            │ ╰─╯
            ╰──╯

            sage: G = [-1, 2, -3, 4, 6, -7, 8, 1, -2, 3, -4, -5, 7, -8, 5, -6]
            sage: K = Knots().from_gauss_code(G)
            sage: unicode_art(K)
               ╭─────╮
               │╭─────╮
              ╭─│╮   ││
             ╭─╯││   ││
             │╰─╯│   ││
             │   ╰──╮││
             │     ╭│╯│
             │    ╭│╯ │
            ╭─────│╯  │
            │╰────╯   │
            ╰─────────╯

        TESTS::

            sage: W = Knots()
            sage: unicode_art(W.one())
            ╭╮
            ╰╯
        """
        style = 2  # among 0, 1, 2 (how to display crossings, see below)
        gauss = self.gauss_code()
        if not gauss:
            gauss = []
        else:
            gauss = gauss[0]

        graphe, (hori, vert) = rectangular_diagram(gauss)
        maxx, maxy = 0, 0
        for a, b in graphe:
            maxx = max(a, maxx)
            maxy = max(b, maxy)
        M = [[" " for a in range(maxy + 1)] for b in range(maxx + 1)]
        for a, b in graphe:
            (x, y), (xx, yy) = graphe.neighbors((a, b))
            if x != a:
                x, y, xx, yy = xx, yy, x, y
            if y < b:
                if xx < a:
                    M[a][b] = u"╯"
                else:
                    M[a][b] = u"╮"
            else:
                if xx < a:
                    M[a][b] = u"╰"
                else:
                    M[a][b] = u"╭"

        for ab, cd in graphe.edge_iterator(labels=False):
            a, b = ab
            c, d = cd
            if a == c:
                b, d = sorted((b, d))
                for i in range(b + 1, d):
                    M[a][i] = u"─"
            else:
                a, c = sorted((a, c))
                for i in range(a + 1, c):
                    M[i][b] = u"│"

        if style == 0:
            H = u"┿"
            V = u"╂"
        elif style == 1:
            H = u"━"
            V = u"┃"
        elif style == 2:
            H = u"─"
            V = u"│"

        for x, y in hori:
            M[x][y] = H
        for x, y in vert:
            M[x][y] = V

        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt([''.join(ligne) for ligne in M])

    def dt_code(self):
        """
        Return the DT code of ``self``.

        ALGORITHM:

        The DT code is generated by the following way:

        Start moving along the knot, as we encounter the crossings we
        start numbering them, so every crossing has two numbers assigned to
        it once we have traced the entire knot. Now we take the even number
        associated with every crossing.

        The following sign convention is to be followed:

        Take the even number with a negative sign if it is an overcrossing
        that we are encountering.

        OUTPUT: DT code representation of the knot

        EXAMPLES::

            sage: K = Knot([[1,5,2,4],[5,3,6,2],[3,1,4,6]])
            sage: K.dt_code()
            [4, 6, 2]
            sage: B = BraidGroup(4)
            sage: K = Knot(B([1, 2, 1, 2]))
            sage: K.dt_code()
            [4, -6, 8, -2]
            sage: K = Knot([[[1, -2, 3, -4, 5, -1, 2, -3, 4, -5]],
            ....:          [1, 1, 1, 1, 1]])
            sage: K.dt_code()
            [6, 8, 10, 2, 4]
        """
        b = self.braid().Tietze()
        N = len(b)
        label = [0 for i in range(2 * N)]
        string = 1
        next_label = 1
        type1 = 0
        crossing = 0
        while next_label <= 2 * N:
            string_found = False
            for i in range(crossing, N):
                if abs(b[i]) == string or abs(b[i]) == string - 1:
                    string_found = True
                    crossing = i
                    break
            if not string_found:
                for i in range(0, crossing):
                    if abs(b[i]) == string or abs(b[i]) == string - 1:
                        string_found = True
                        crossing = i
                        break
            assert label[2 * crossing + next_label % 2] != 1, "invalid knot"

            label[2 * crossing + next_label % 2] = next_label
            next_label += 1
            if type1 == 0:
                if b[crossing] < 0:
                    type1 = 1
                else:
                    type1 = -1
            else:
                type1 = -type1
                if (abs(b[crossing]) == string) ^ (b[crossing] * type1 < 0):
                    if next_label % 2:
                        label[2 * crossing] = label[2 * crossing] * -1
            if abs(b[crossing]) == string:
                string += 1
            else:
                string -= 1
            crossing += 1
        code = [0 for _ in range(N)]
        for i in range(N):
            for j in range(N):
                if label[2 * j + 1] == 2 * i + 1:
                    code[i] = label[2 * j]
                    break
        return code

    def arf_invariant(self):
        """
        Return the Arf invariant.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: K = Knot(B([-1, 2, 1, 2]))
            sage: K.arf_invariant()
            0
            sage: B = BraidGroup(8)
            sage: K = Knot(B([-2, 3, 1, 2, 1, 4]))
            sage: K.arf_invariant()
            0
            sage: K = Knot(B([1, 2, 1, 2]))
            sage: K.arf_invariant()
            1
        """
        a = self.alexander_polynomial()(-1)
        if (a % 8) == 1 or (a % 8) == 7:
            return 0
        return 1

    def colored_jones_polynomial(self, N, variab=None, try_inverse=True):
        r"""
        Return the colored Jones polynomial of the trace closure of the braid.

        INPUT:

        - ``N`` -- integer; the number of colors
        - ``variab`` -- (default: `q`) the variable in the resulting
          Laurent polynomial
        - ``try_inverse`` -- boolean (default: ``True``); if ``True``,
          attempt a faster calculation by using the inverse of the braid

        ALGORITHM:

        The algorithm used is described in [HL2018]_ for the corresponding
        braid representation. We follow their notation, but work in a
        suitable free algebra over a Laurent polynomial ring in one
        variable to simplify bookkeeping.

        EXAMPLES::

            sage: W = Knots()
            sage: K = W.from_dowker_code([-4,-6,-2])
            sage: K.colored_jones_polynomial(2)
            -q^-4 + q^-3 + q^-1
            sage: K.colored_jones_polynomial(2, 't')
            -t^-4 + t^-3 + t^-1
            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: K.colored_jones_polynomial(2, -t)
            -t^-4 - t^-3 - t^-1

            sage: R.<t> = ZZ[]
            sage: K.colored_jones_polynomial(2, t+1)
            (t^3 + 3*t^2 + 4*t + 1)/(t^4 + 4*t^3 + 6*t^2 + 4*t + 1)
        """
        return self.braid().colored_jones_polynomial(N=N, variab=variab,
                                                     try_inverse=try_inverse)

    def connected_sum(self, other):
        r"""
        Return the oriented connected sum of ``self`` and ``other``.

        .. NOTE::

            We give the knots an orientation based upon the braid
            representation.

        INPUT:

        - ``other`` -- a knot

        OUTPUT:

        A knot equivalent to the connected sum of ``self`` and ``other``.

        EXAMPLES::

            sage: B = BraidGroup(2)
            sage: trefoil = Knot(B([1,1,1]))
            sage: K = trefoil.connected_sum(trefoil); K
            Knot represented by 6 crossings
            sage: K.braid()
            s0^3*s1^-1*s0^3*s1

        .. PLOT::
            :width: 300 px

            B = BraidGroup(2)
            trefoil = Knot(B([1,1,1]))
            K = trefoil.connected_sum(trefoil)
            sphinx_plot(K.plot())

        ::

            sage: rev_trefoil = Knot(B([-1,-1,-1]))
            sage: K2 = trefoil.connected_sum(rev_trefoil); K2
            Knot represented by 6 crossings
            sage: K2.braid()
            s0^3*s1^-1*s0^-3*s1

        .. PLOT::
            :width: 300 px

            B = BraidGroup(2)
            t = Knot(B([1,1,1]))
            tr = Knot(B([-1,-1,-1]))
            K2 = t.connected_sum(tr)
            sphinx_plot(K2.plot())

        Observe that both knots have according ``dowker_notation`` (showing that
        the constructing from DT-code may not be unique for non prime knots, see
        :meth:`from_dowker_code`)::

            sage: K.dowker_notation()
            [(4, 1), (2, 5), (6, 3), (10, 7), (8, 11), (12, 9)]
            sage: K2.dowker_notation()
            [(4, 1), (2, 5), (6, 3), (7, 10), (11, 8), (9, 12)]
            sage: K.homfly_polynomial() == K2.homfly_polynomial()
            False

        TESTS::

            sage: B = BraidGroup(2)
            sage: trivial = Knots().one()
            sage: trivial * trivial
            Knot represented by 0 crossings
            sage: trefoil = Knot(B([1,1,1]))
            sage: trefoil * trivial
            Knot represented by 3 crossings
            sage: trefoil * trefoil
            Knot represented by 6 crossings

        REFERENCES:

        - :wikipedia:`Connected_sum`
        """
        from sage.functions.generalized import sign
        ogc1 = self.oriented_gauss_code()
        ogc2 = other.oriented_gauss_code()
        if not ogc1[0]:
            return other
        if not ogc2[0]:
            return self
        # how much we have to "displace" the numbering of the
        # crossings of other
        m1 = max(abs(i) for i in ogc1[0][0])
        m2 = min(abs(i) for i in ogc2[0][0])
        n = m1 - m2 + 1
        # construct the oriented Gauss code of the result
        ogc2_0_0 = [a + int(sign(a)) * n for a in ogc2[0][0]]
        nogc = [[ogc1[0][0] + ogc2_0_0], ogc1[1] + ogc2[1]]
        return type(self)(nogc)

    _mul_ = connected_sum


class Knots(Singleton, Parent):
    """
    The set for all knots, as a monoid for the connected sum.
    """
    def __init__(self):
        """
        TESTS::

            sage: S = Knots()
            sage: S.cardinality()
            +Infinity
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=Monoids().Infinite())

    def _repr_(self):
        r"""
        TESTS::

            sage: Knots()
            Knots
        """
        return "Knots"

    def one(self):
        """
        Return the unit of the monoid.

        This is the trivial knot.

        EXAMPLES::

            sage: Knots().one()
            Knot represented by 0 crossings
        """
        return self.element_class([])

    def an_element(self):
        """
        Return the trefoil knot.

        EXAMPLES::

            sage: Knots().an_element()
            Knot represented by 3 crossings
        """
        return self.element_class([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])

    def from_gauss_code(self, gauss):
        """
        Build a knot from a signed Gauss code.

        This makes some arbitrary choice of orientation.

        INPUT:

        - a signed Gauss code

        OUTPUT:

        - a knot

        EXAMPLES::

            sage: W = Knots()
            sage: K1 = W.from_gauss_code([2, -1, 3, -2, 1, -3])
            sage: K1.alexander_polynomial()
            t^-1 - 1 + t
        """
        orientations = recover_orientations(gauss)[3]
        return Knot([[gauss], orientations])

    def from_dowker_code(self, code):
        """
        Build a knot from a Dowker-Thistlethwaite code.

        The Dowker-Thistlethwaite code of a knot diagram is defined as follows.

        Start following the knot diagram at some regular point. Label the
        crossings by a number (starting from number 1) in the order in
        which they are met. At the end, every crossing gets numbered
        twice, once by an even number and once by an odd number. When
        meeting an over-crossing with even number, use instead the
        negative of this even number as label.

        Then the set of crossings gives a set of pairs (odd,
        even). Sort this set according to the odd component, and then
        keep only the even components in the same order. This is the
        Dowker-Thistlethwaite code.

        INPUT:

        a list of signed even numbers, the Dowker-Thistlethwaite code of a knot

        OUTPUT:

        a knot

        .. WARNING::

            In general the Dowker-Thistlethwaite code does not describe a knot
            uniquely. It is not only insensitive on mirror images, but may also
            mix up non prime knots. For example ``[4, 6, 2, 10, 12, 8]`` describes
            the connected sum of two trefoil knots, as well as the connected sum
            of a trefoil with its mirror (see the corresponding example in the
            documentation of :meth:`connected_sum`).

        EXAMPLES::

            sage: W = Knots()
            sage: K1 = W.from_dowker_code([8,10,2,12,4,6])
            sage: K1.dowker_notation()
            [(5, 2), (9, 4), (11, 6), (1, 8), (3, 10), (7, 12)]

            sage: W.from_dowker_code([6,10,14,12,16,2,18,4,8])
            Knot represented by 9 crossings

            sage: W.from_dowker_code([4,8,10,-14,2,-16,-18,-6,-12])
            Knot represented by 9 crossings

            sage: K3 = W.from_dowker_code([6,-12,2,8,-4,-10]); K3
            Knot represented by 6 crossings
            sage: K3.dowker_notation()
            [(5, 2), (4, 9), (1, 6), (7, 8), (10, 11), (12, 3)]

        .. SEEALSO:: :meth:`~sage.knots.knot.Knot.dowker_notation`

        REFERENCES:

        - :wikipedia:`Dowker_notation`

        - http://katlas.org/wiki/DT_(Dowker-Thistlethwaite)_Codes
        """
        gauss = dowker_to_gauss(code)
        orientations = recover_orientations(gauss)[3]
        return Knot([[gauss], orientations])

    def from_table(self, n, k):
        """
        Return a knot from its index in the Rolfsen table.

        INPUT:

        - ``n`` -- the crossing number
        - ``k`` -- a positive integer

        OUTPUT:

        the knot `K_{n,k}` in the Rolfsen table

        EXAMPLES::

            sage: K1 = Knots().from_table(6,3); K1
            Knot represented by 6 crossings
            sage: K1.alexander_polynomial()
            t^-2 - 3*t^-1 + 5 - 3*t + t^2

            sage: K2 = Knots().from_table(8,4); K2
            Knot represented by 9 crossings
            sage: K2.determinant()
            19
            sage: K2.signature()
            2

            sage: K3 = Knots().from_table(10,56); K3
            Knot represented by 11 crossings
            sage: K3.jones_polynomial()
            t^10 - 3*t^9 + 6*t^8 - 9*t^7 + 10*t^6 - 11*t^5 + 10*t^4 - 7*t^3
            + 5*t^2 - 2*t + 1

            sage: K4 = Knots().from_table(10,100)
            sage: K4.genus()
            4

        TESTS::

            sage: Knots().from_table(6,6)
            Traceback (most recent call last):
            ...
            ValueError: not found in the knot table

            sage: Knots().from_table(12,6)
            Traceback (most recent call last):
            ...
            ValueError: more than 10 crossings, not in the knot table

        REFERENCES:

        - KnotAtlas, http://katlas.math.toronto.edu/wiki/The_Rolfsen_Knot_Table
        """
        if n > 10:
            raise ValueError('more than 10 crossings, not in the knot table')
        from sage.groups.braid import BraidGroup
        if (n, k) in small_knots_table:
            m, word = small_knots_table[(n, k)]
            G = BraidGroup(m)
            return Knot(G(word))
        else:
            raise ValueError('not found in the knot table')

    Element = Knot
