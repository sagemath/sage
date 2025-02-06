# -*- coding: utf-8 -*-
r"""
Methods for calculations with elements of a ChowRing

Provide methods for calculations with elements of a ChowRing.
For example, suppose we work with the ChowRing of the projective plane blown
up in a point `\widehat{\mathbb{P}^2}` and consider the special element
`x = (2*H - E)^2 + E`::

    sage: A.<H,E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
    sage: x = (2 * H - E)^2 + E; x
    -3*E^2 + E

Compute the decomposition of x in the graded ring A and return its part in
degree equal to the dimension::

    sage: x.by_degrees()
    [0, E, -3*E^2]

The codimension is the top degree of an element::

    sage: x.codimension()
    2
    sage: (E + H).codimension()
    1

Compute the integral of x with respect to the point class `H^2`::

    sage: A.set_point_class(H^2)        # Need to set the point class first.
    sage: x.integral()
    3

The other methods are related to the elements defined as the Chern character
of a vector bundle. Recall that if `X` is a smooth algebraic variety, then the
Chern character gives an isomorphism

.. math::

    ch: K_{0}(X)\otimes\QQ\longrightarrow A^{*}(X)\otimes\QQ

In our example if `T` is the tangent bundle on `\mathbb{P}^2` then `T` is of
rank `2` with Chern classes `[1, 3*H - E, 4*H^2]`  and its Chern character is
`t = 2 - E + 3 * H`::

    sage: A.set_dimension(2)    # Need to set the dimension first.
    sage: t = 2 - E + 3 * H
    sage: t._expp()
    -4*E^2 + 3*H - E + 1
    sage: t.wedge(2)._expp()
    3*H - E + 1
    sage: t.symm(2)._expp()
    -32*E^2 + 9*H - 3*E + 1


AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""

# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.functions.other import factorial
from sage.symbolic.ring import SR
from sage.functions.log import exp
from sage.rings.rational_field import QQ
from sage.libs.singular.function import singular_function
from sage.structure.richcmp import richcmp


class ChowRingElement(QuotientRingElement):
    """
    Class representing the elements of a ChowRing.
    """
    def __init__(self, parent, rep):
        r"""
        Construct a :class:`ChowRingElement`.


        INPUT:

        - ``parent`` -- a ChowRing;

        - ``rep'' -- a representing of an element of the ChowRing.

        OUTPUT:

        - :class:`ChowRing <ChowRing_generic>`.

        EXAMPLES::

            sage: A = ChowRing('h', 1, 'h^6')
            sage: h = A('h')

        TESTS::

            sage: A = ChowRing()
            sage: A.set_point_class(1)
            sage: A(5).integral()
            5

            sage: A.<h> = ChowRing('h', 1, 'h^6')
            sage: x = 2 + 3*h - 5*h^3
            sage: TestSuite(x).run()
        """
        self._parentchowring = parent
        # Reduce is broken (and not needed) for parents with no generators.
        red = False if parent.ngens() == 0 else True
        QuotientRingElement.__init__(self, parent, rep, reduce=red)

    def _richcmp_(self, other, op):
        r"""
        Compares two elements.

        EXAMPLES::

            sage: A.<H,E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: H^2 == -E^2
            True
            sage: 2 + H^2 == E - E^2
            False

        TESTS::

            sage: A = ChowRing()
            sage: A(5) == A(0)
            False
        """
        A = self._parentchowring
        if A.ngens() == 0:
            return richcmp(QQ(str(self)), QQ(str(other)), op)
        return QuotientRingElement._richcmp_(self, other, op)

    @cached_method
    def by_degrees(self):
        r"""
        Returns this Chow ring element as a list `[e_0,\dots,e_{d}]` where `e_k`
        has degree `k`.

        EXAMPLES::

            sage: A.<H,E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: (-2 - E + 3*H - 4*H^2).by_degrees()
            [-2, 3*H - E, 4*E^2]
            sage: (-E + 3*H - 4*H^2).by_degrees()
            [0, 3*H - E, 4*E^2]
            sage: (-2 - E + 3*H).by_degrees()
            [-2, 3*H - E]

        TESTS::

            sage: A = ChowRing()
            sage: A(5).by_degrees()
            [5]
        """
        A = self._parentchowring
        # The computation is done in the cover ring in order to be able
        # to use the methods monomials and coefficients, then the list
        # is returned as a list of elements of the parent chow ring.
        monoms, coeffs = self.lift().monomials(), self.lift().coefficients()
        l = len(monoms)
        d = max([monoms[i].degree() for i in range(l)]) if l else 0
        p = [0] * (d + 1)
        for i in range(len(monoms)):
            k = monoms[i].degree()
            if k <= d:
                p[k] += coeffs[i] * monoms[i]
        return [A(v) for v in p]

    def truncate(self, a=0, b=0):
        r"""
        Return this Chow ring element truncated in degrees below a and above b
        (strictly).

        INPUT:

        - ``a`` -- an optional integer (defaults to 0)

        - ``b`` -- an optional integer (defaults to 0)

        OUTPUT:

        - a Chow ring element.

        EXAMPLES::

            sage: A.<H,E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: x = -2 - E + 3*H - 4*H^2; x.by_degrees()
            [-2, 3*H - E, 4*E^2]
            sage: x.truncate()
            -2
            sage: x.truncate(0, 1)
            3*H - E - 2
            sage: x.truncate(0, 2)
            4*E^2 + 3*H - E - 2
            sage: x.truncate(1, 2)
            4*E^2 + 3*H - E
        """
        return sum([x for x in self.by_degrees()[a:b + 1]])

    @cached_method
    def integral(self):
        r"""
        Return the integral of this Chow ring element with respect to the
        parents point_class.

        EXAMPLES::

            sage: A.<H,E> = ChowRing(['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'])
            sage: A.set_point_class('H^2')
            sage: x = (2*H-E)^2; x.integral()
            3

        TESTS::

            sage: A = ChowRing()
            sage: A.set_point_class(2)
            sage: A(6).integral()
            3
        """
        A = self._parentchowring
        if A.point_class() is None:
            raise ValueError("Need a point_class before calling integral.")
        # We compute in the covering (in order to be able to call 'division')
        R, p, s = A.cover_ring(), A.point_class().lift(), self.lift()
        division = singular_function('division')
        return s / p if R.ngens() == 0 else QQ(division(s, p, ring=R)[0][0, 0])

    @cached_method
    def codimension(self):
        """
        Return the codimension of this Chow ring element.

        EXAMPLES::

            sage: A.<H> = ChowRing('H', 1, 'H^5')
            sage: H.codimension()
            1
            sage: (H + H^2).codimension()
            2

        TESTS::

            sage: A = ChowRing(0)
            sage: A(5).codimension()
            0
        """
        return self.lift().degree()

    @cached_method
    def _expp(self):
        r"""
        Internal method (see also the Schubert package by Katz and Strømme).

        Chern character --> Total Chern Class

        INPUT:

        - ``dim`` -- an integer

        OUTPUT:

        - A list of elements of the parent Chow ring.

        EXAMPLES::

            sage: A.<h> = ChowRing('h', 1, 'h^3')  # P2
            sage: A.set_dimension(2)
            sage: t = 3/2*h^2 + 3*h + 2  # ch(TP^2)
            sage: t._expp()
            3*h^2 + 3*h + 1

        TESTS::

            sage: A = ChowRing()
            sage: A.set_dimension(0)
            sage: A(5)._expp()
            1
        """
        A, sbd = self._parentchowring, self.by_degrees()
        d = A.dimension()
        if d is None:
            err = "Need a dimension before calling chern_classes."
            raise ValueError(err)
        # We calculate in the cover ring (no relations).
        R = A.cover_ring()
        p = [x.lift() for x in sbd] + [R(0)] * (d - len(sbd) + 1)
        for k in range(d + 1):
            p[k] *= (-1) ** k * factorial(k)
        cc = [R(1)] + [R(0)] * d
        for k in range(1, d + 1):
            cc[k] = - sum(cc[k - i] * p[i] for i in range(1, k + 1)) / k
        return sum([A(c) for c in cc])

    @cached_method
    def _logg(self):
        r"""
        Internal method (see also the Schubert package by Katz and Strømme).

        Total Chern Class --> Chern Character

        EXAMPLES::

            sage: A.<h> = ChowRing('h', 1, 'h^3')  # P2
            sage: A.set_dimension(2)
            sage: t = A(1) + 3*h + 3*h^2
            sage: t._logg()
            3/2*h^2 + 3*h

        TESTS:

            sage: t = A(1)
            sage: t._logg()
            0

        """
        A, sbd = self._parentchowring, self.by_degrees()
        d = A.dimension()
        if d is None:
            err = "Need a dimension before calling _chern_character."
            raise ValueError(err)
        # We calculate in the cover ring (no relations).
        R = A.cover_ring()
        cc = [R(c.lift()) for c in sbd] + [R(0)] * (d + 1 - len(sbd))
        p = [0] * (d + 1)
        for k in range(1, d + 1):
            p[k] = - k * cc[k] - sum(p[i] * cc[k - i] for i in range(1, k))
        result = sum((-1) ** i * p[i] / factorial(i) for i in range(1, d + 1))
        return A(result)

    @cached_method
    def todd(self):
        r"""
        Return the todd class of this Chow ring element.

        OUTPUT:

        - An element of the parent Chow ring.

        EXAMPLE::

            sage: A.<c1, c2, c3> = ChowRing(['c1', 'c2', 'c3'], [1, 2, 3])
            sage: A.set_dimension(3)
            sage: x = 1 + c1 + 1/2*c1^2 + 1/6*c1^3

            sage: x._expp()  # Check that x=ch(L) with L l.b. c1(L)=c1
            c1 + 1
            sage: x.todd()
            1/12*c1^2 + 1/2*c1 + 1

        TEST::

            sage: A = ChowRing()
            sage: A.set_dimension(0)
            sage: A(7).todd()
            1
        """
        A, sbd = self._parentchowring, self.by_degrees()
        if A.dimension() is None:
            err = "Need a dimension before calling todd."
            raise ValueError(err)
        dim = A.dimension()
        ch = [x.lift() for x in sbd] + [A(0)] * (dim - len(sbd) + 1)
        # Get the Taylor development
        x = SR('x')
        f = x / (SR(1) - exp(-x))
        g = f.taylor(x, 0, dim + 2)
        # Apply char
        cc = [g.coefficient(x, i) for i in range(dim + 1)]
        p, cc = [0] * (dim + 1), cc + [0] * (dim + 1 - len(cc))
        for k in range(1, dim + 1):
            p[k] = - k * cc[k] - sum(p[i] * cc[k - i] for i in range(1, k))
        res = sum((-1) ** i * p[i] * SR(ch[i]) for i in range(1, dim + 1))
        # Apply classes
        return A(str(res))._expp()

    @cached_method
    def adams(self, k):
        r"""
        Return the result of the k-th Adams operator applied to this Chow ring
        element.

        INPUT :

        - ``k``-- an integer

        EXAMPLE::

            sage: A.<c1, c2, c3> = ChowRing(['c1', 'c2', 'c3'], [1, 2, 3])
            sage: (1 + c1 + c2).adams(2)
            4*c2 + 2*c1 + 1
        """
        ch = self.by_degrees()
        return sum([k ** i * ch[i] for i in range(len(ch))])

    @cached_method
    def wedge(self, p):
        r"""
        Return the p-th exterior power of this Chow ring element.
        Truncate above degree d.

        INPUT:

        - ``p``-- an integer

        By convention, `0` is returned for negative values of p.

        EXAMPLES::

            sage: A.<h> = ChowRing('h', 1, 'h^3')  # P2
            sage: A.set_dimension(2)
            sage: t = 3/2*h^2 + 3*h + 2  # ch(TP^2)
            sage: t._expp()
            3*h^2 + 3*h + 1
            sage: t.wedge(0)._expp()
            1
            sage: t.wedge(1)._expp()
            3*h^2 + 3*h + 1
            sage: t.wedge(2)._expp()
            3*h + 1
            sage: t.wedge(3)._expp()
            Traceback (most recent call last):
            ...
            ValueError: Wedge for powers outside 0..rank not allowed.

        TESTS::

            sage A = ChowRing(0)
            sage A(5).wedge(0)
            1
        """
        A, rank = self._parentchowring, int(self.by_degrees()[0])
        if p < 0 or p > rank:
            raise ValueError('Wedge for powers outside 0..rank not allowed.')
        if p == 0:
            return A(1)
        if p == 1:
            return self
        if A.dimension() is None:
            err = "Need a dimension before calling wedge."
            raise ValueError(err)
        d = A.dimension()
        return A(sum((-1) ** (p - j + 1) *
                     (self.wedge(j) * self.adams(p - j)).truncate(0, d)
                     for j in range(p)) / p)

    @cached_method
    def symm(self, p):
        r"""
        Return the p-th symmetric power of this Chow ring element.
        Truncate above degree d.

        INPUT:

        - ``p``-- an integer

        - ``d``-- an integer

        OUTPUT:

        - A Chow ring element.

        EXAMPLES::

            sage: A.<h> = ChowRing('h', 1, 'h^3')  # P2
            sage: A.set_dimension(2)
            sage: t = 1/2*h^2 + h + 1  # ch(O(1))
            sage: t._expp()
            h + 1
            sage: t.symm(0)._expp()
            1
            sage: t.symm(1)._expp()
            h + 1
            sage: t.symm(2)._expp()
            2*h + 1
            sage: t.symm(5) == t^5  # t corresponds to a line bundle.
            True
        """
        A, rank = self._parentchowring, int(self.by_degrees()[0])
        if p < 0:
            raise ValueError('Symm for negative powers not allowed.')
        if p == 0:
            return A(1)
        if p == 1:
            return self
        r = min(p, rank)
        if A.dimension() is None:
            err = "Need a dimension before calling wedge."
            raise ValueError(err)
        d = A.dimension()
        return A(sum((-1) ** (j + 1) *
                     (self.wedge(j) * self.symm(p - j)).truncate(0, d)
                     for j in range(1, r + 1)))
