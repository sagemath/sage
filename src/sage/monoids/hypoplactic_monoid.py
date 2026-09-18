r"""
Hypoplactic monoid

This file implements the hypoplactic monoid on the alphabet
`\{1, 2, \ldots, n\}`. Elements are represented by words, with equality
determined by comparing the quasi-ribbon tableaux obtained from
Krob--Thibon insertion. Multiplication is induced by concatenation of
words, followed by replacing the product with its quasi-ribbon reading word
representative. For references, see [KT1997]_ and [Nov2000]_.


AUTHORS:

- Daniel Chen, Lisa Johnston, Junbok Lee, Evuilynn Nguyen, Heather Ross,
  Anne Schilling, Chenchen Zhao (2026): initial version
"""
# ****************************************************************************
#       Copyright (C) 2026 Daniel Chen, Lisa Johnston, Junbok Lee, Evuilynn Nguyen, Heather Ross, Anne Schilling, Chenchen Zhao
#
#  Distributed under the terms of the GNU General Public License (GPL)
#              https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau, QuasiRibbonTableaux
from sage.monoids.plactic_monoid import WordMonoid, WordMonoidElement

from sage.structure.parent import Parent
from sage.categories.sets_with_grading import SetsWithGrading
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.combinat.family import Family
from sage.misc.cachefunc import cached_method


class HypoplacticMonoid(WordMonoid):
    r"""
    The hypoplactic monoid on the alphabet `\{1, 2, \ldots, n\}`.

    INPUT:

    - ``n`` -- a positive integer; the size of the alphabet

    OUTPUT:

    The hypoplactic monoid of rank ``n``.

    The hypoplactic monoid is a quotient of the free monoid on the alphabet
    `\{1, 2, \ldots, n\}`. In this implementation, elements are represented
    by words in the alphabet `\{1, 2, \ldots, n\}`. Equality is determined by
    comparing the quasi-ribbon tableaux obtained from Krob--Thibon insertion.

    The identity element is the empty word. Multiplication is induced by
    concatenation of words. The product is stored using the quasi-ribbon
    reading word representative obtained from hypoplactic insertion.

    EXAMPLES::

        sage: H = HypoplacticMonoid(4)
        sage: H
        Hypoplactic monoid of rank 4
        sage: H.rank()
        4

    Elements are constructed from tuples::

        sage: x = H([3, 2, 2, 1])
        sage: x
        3221
        sage: x.to_tableau()
        [[1], [2, 2], [None, 3]]
        sage: x.to_word()
        2132

    Two words represent the same hypoplactic element when they have the same
    quasi-ribbon insertion tableau::

        sage: H([3, 2, 2, 1]) == H([2, 3, 1, 2])
        True
        sage: H([3, 2, 2, 1]) == H([1, 2, 2, 3])
        False

    Multiplication is induced by concatenation, followed by replacing the
    result with its quasi-ribbon reading word representative::

        sage: H([3]) * H([4])
        34
        sage: H([4]) * H([3])
        43

    TESTS::

        sage: H = HypoplacticMonoid(4)
        sage: H([]) == H.one()
        True
        sage: len(H.one())
        0
        sage: H.rank()
        4
        sage: H.one() * H([3, 2, 2, 1]) == H([3, 2, 2, 1])
        True
        sage: H([3, 2, 2, 1]) * H.one() == H([3, 2, 2, 1])
        True

        sage: HypoplacticMonoid(0)
        Traceback (most recent call last):
        ...
        ValueError: the rank must be a positive integer

        sage: H([1, 2, 5])
        Traceback (most recent call last):
        ...
        ValueError: letters must be integers from 1 to 4

        sage: H([0, 1])
        Traceback (most recent call last):
        ...
        ValueError: letters must be integers from 1 to 4

        sage: H([-1, 2])
        Traceback (most recent call last):
        ...
        ValueError: letters must be integers from 1 to 4

        sage: H([1, 'a'])
        Traceback (most recent call last):
        ...
        ValueError: letters must be integers from 1 to 4

        sage: H([1, 2.5])
        Traceback (most recent call last):
        ...
        ValueError: letters must be integers from 1 to 4

        sage: TestSuite(H).run() # long time
    """

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: HypoplacticMonoid(4)
            Hypoplactic monoid of rank 4
        """
        return f"Hypoplactic monoid of rank {self._n}"

    def subset(self, k):
        r"""
        Return the hypoplactic monoid elements represented by words of length ``k``.

        Since the hypoplactic monoid is infinite, this returns the finite set of
        elements of a fixed size, using their canonical reading word representatives.

        EXAMPLES::

            sage: H = HypoplacticMonoid(2)
            sage: H.subset(1)
            Lazy family (to_word(i))_{i in Quasi-ribbon tableaux of size 1 with entries at most 2}
            sage: list(H.subset(1))
            [1, 2]
            sage: H.subset(2)
            Lazy family (to_word(i))_{i in Quasi-ribbon tableaux of size 2 with entries at most 2}
            sage: list(H.subset(2))
            [11, 12, 22, 21]
        """
        if not isinstance(k, (int, Integer)) or k < 0:
            raise ValueError("the size must be a nonnegative integer")

        quasiribbontableaux = QuasiRibbonTableaux(size=k, max_entry=self.rank())

        def to_word(t):
            return self(t.to_word_by_column())
        return Family(quasiribbontableaux, to_word, lazy=True)

    class Element(WordMonoidElement):
        r"""
        An element of a hypoplactic monoid.

        Elements are represented by words in the alphabet
        `\{1, 2, \ldots, n\}`.

        EXAMPLES::

            sage: H = HypoplacticMonoid(4)
            sage: x = H([3, 2, 2, 1])
            sage: x
            3221
            sage: parent(x)
            Hypoplactic monoid of rank 4
        """
        @cached_method
        def to_tableau(self):
            """
            Return the quasi-ribbon insertion tableau corresponding to ``self``.

            The tableau is computed using Krob--Thibon insertion.

            OUTPUT:

            The quasi-ribbon tableau obtained by inserting the word
            representing ``self``.

            EXAMPLES::

                sage: H = HypoplacticMonoid(4)
                sage: H([3, 2, 2, 1]).to_tableau()
                [[1], [2, 2], [None, 3]]
                sage: H([3, 4, 3, 2, 1, 2]).to_tableau()
                [[1], [2, 2], [None, 3, 3], [None, None, 4]]

            TESTS::

                sage: H = HypoplacticMonoid(4)
                sage: H([]).to_tableau()
                []
            """
            Q = QuasiRibbonTableaux()
            return Q.insert_word(self.value)

        def to_word(self):
            """
            Return the quasi-ribbon reading word representative of ``self``.

            The reading word is obtained from the quasi-ribbon insertion
            tableau by reading columns from left to right, and from bottom to
            top within each column.

            OUTPUT:

            A tuple containing the quasi-ribbon reading word of ``self``.

            EXAMPLES::

                sage: H = HypoplacticMonoid(4)
                sage: H([3, 2, 2, 1]).to_word()
                2132
                sage: H([3, 4, 3, 2, 1, 2]).to_word()
                213243

            TESTS::

                sage: H = HypoplacticMonoid(4)
                sage: H([]).to_word()
            """
            parent = self.parent()
            return parent(list(self.to_tableau().to_word_by_column()))
