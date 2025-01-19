r"""
This module implements `\Q/n\Z` for `n \in \Q`.

When `n \in \Z`, you can construct these groups as follows::

    sage: G = QQ/ZZ; G
    Q/Z
    sage: QQ/(2*ZZ)
    Q/2Z

You can create random elements::

    sage: all(G.random_element().parent() is G for _ in range(4))
    True

There is an iterator over the (infinitely many) elements::

    sage: import itertools
    sage: list(itertools.islice(G, 10r))
    [0, 1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5]
"""

# ****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from .qmodnz_element import QmodnZ_Element


class QmodnZ(Parent, UniqueRepresentation):
    r"""
    The ``QmodnZ`` class represents the abelian group `\Q/n\Z`.

    INPUT:

    The constructor may be called in any of the following ways.

    #. ``QmodnZ(n)``, where

        - ``n`` -- a rational number (including 0 or negative rational numbers)

    #. ``QQ/(n*ZZ)``, where

        - ``n`` -- integer (including 0 or negative integers)

    OUTPUT: the abelian group `\Q/n\Z`

    EXAMPLES::

        sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
        sage: QQ/(19*ZZ)
        Q/19Z

        sage: QmodnZ(19)
        Q/19Z

        sage: QmodnZ(2/3)
        Q/(2/3)Z
    """
    Element = QmodnZ_Element

    def __init__(self, n=1):
        r"""
        Initialization.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(2)
            sage: G
            Q/2Z

        TESTS::

            sage: G = QQ/(19*ZZ)
            sage: TestSuite(G).run()
        """
        self.n = QQ(n).abs()
        category = CommutativeAdditiveGroups().Infinite()
        Parent.__init__(self, base=ZZ, category=category)
        self._populate_coercion_lists_(coerce_list=[QQ])

    def _repr_(self):
        r"""
        Display the group.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(1); G
            Q/Z

            sage: G = QQ/(3*ZZ); G
            Q/3Z

            sage: G = QmodnZ(1/5); G
            Q/(1/5)Z
        """
        if self.n == 1:
            return "Q/Z"
        elif self.n in ZZ:
            return "Q/%sZ" % self.n
        else:
            return "Q/(%s)Z" % self.n

    def _coerce_map_from_(self, S):
        r"""
        Coercion from a parent ``S``.

        There is a coercion from ``S`` if ``S`` has a coerce map to `\Q`
        or if `S = \Q/m\Z` for `m` a multiple of `n`.

        TESTS::

            sage: G2 = QQ/(2*ZZ)
            sage: G3 = QQ/(3*ZZ)
            sage: G4 = QQ/(4*ZZ)
            sage: G2.has_coerce_map_from(QQ)
            True
            sage: G2.has_coerce_map_from(ZZ)
            True
            sage: G2.has_coerce_map_from(ZZ['x'])
            False
            sage: G2.has_coerce_map_from(G3)
            False
            sage: G2.has_coerce_map_from(G4)
            True
            sage: G4.has_coerce_map_from(G2)
            False
        """
        if QQ.has_coerce_map_from(S):
            return True
        if isinstance(S, QmodnZ) and (S.n / self.n in ZZ):
            return True

    def _element_constructor_(self, x):
        r"""
        Construct an element in `\Q/n\Z`.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(2/3)
            sage: G(5/6)
            1/6
        """
        return self.element_class(self, QQ(x))

    def an_element(self):
        """
        Return an element, for use in coercion system.

        TESTS::

           sage: (QQ/ZZ).an_element()
           0
        """
        return self(0)

    def some_elements(self):
        """
        Return some elements, for use in testing.

        TESTS::

            sage: L = (QQ/ZZ).some_elements()
            sage: len(L)
            92
        """
        return list({self(x) for x in QQ.some_elements()})

    def random_element(self):
        r"""
        Return a random element of `\Q/n\Z`.

        The denominator is selected
        using the ``1/n`` distribution on integers, modified to return
        a positive value.  The numerator is then selected uniformly.

        EXAMPLES::

            sage: G = QQ/(6*ZZ)
            sage: G.random_element().parent() is G
            True
        """
        if self.n == 0:
            return self(QQ.random_element())
        d = ZZ.random_element()
        if d >= 0:
            d = 2 * d + 1
        else:
            d = -2 * d
        n = ZZ.random_element((self.n * d).ceil())
        return self(n / d)

    def __iter__(self):
        r"""
        Create an iterator that generates the elements of `\Q/n\Z` without
        repetition, organized by increasing denominator.

        For a fixed denominator, elements are listed by increasing numerator.

        EXAMPLES:

        The first 19 elements of `\Q/5\Z`::

            sage: import itertools
            sage: list(itertools.islice(QQ/(5*ZZ), 19r))
            [0, 1, 2, 3, 4, 1/2, 3/2, 5/2, 7/2, 9/2, 1/3, 2/3, 4/3, 5/3,
             7/3, 8/3, 10/3, 11/3, 13/3]
        """
        if self.n == 0:
            for x in QQ:
                yield self(x)
        else:
            d = ZZ(0)
            while True:
                for a in d.coprime_integers((d * self.n).floor()):
                    yield self(a / d)
                d += 1
