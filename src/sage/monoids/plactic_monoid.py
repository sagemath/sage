r"""
Plactic monoid

AUTHORS:

- Daniel Chen, Lisa Johnston, Junbok Lee, Evuilynn Nguyen, Heather Ross, Chenchen Zhao (2026): initial version

This file implements the plactic monoid on the alphabet
`\{1, 2, \ldots, n\}`. Elements are represented by words, with equality
determined by their RSK insertion tableaux. Multiplication is given by
concatenation of words, and the identity element is the empty word.

This file consists of the following major classes:

Parent classes:

* :class:`PlacticMonoid`

Element classes:

* :class:`PlacticMonoid.Element`

The main functionality includes constructing plactic monoid elements,
computing their RSK insertion tableaux, converting elements to their row
reading word representatives, computing shapes and equivalence classes,
testing canonical representatives, and listing all elements of a fixed word
length.
"""

# ****************************************************************************
#       Copyright (C) 2026 Daniel Chen, Lisa Johnston, Junbok Lee, Evuilynn Nguyen, Heather Ross, Chenchen Zhao
#
#  Distributed under the terms of the GNU General Public License (GPL)
#              https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.sets_with_grading import SetsWithGrading

from sage.misc.cachefunc import cached_method
from itertools import permutations, chain
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.sets.family import Family
from sage.combinat.rsk import RSK, RSK_inverse
from sage.combinat.tableau import StandardTableaux, SemistandardTableaux
from sage.combinat.partition import Partitions
from sage.combinat.permutation import Permutations


class WordMonoid(UniqueRepresentation, Parent):
    r"""
    This class is an ancestor class for the plactic and hypoplactic monoid.

    INPUT:

    - ``n`` -- a positive integer; the size of the alphabet

    Elements are represented by words in `\{1, 2, \ldots, n\}`.
    It is assumed that the methods `to_tableau`, `to_word` and `equivalence_class` are
    implemented. Equality is determined by methods `to_tableau`.
    """

    def __init__(self, n):
        """
        Initialize ``self``.

        INPUT:

        - ``n`` -- a positive integer; the size of the alphabet

        EXAMPLES::

            sage: P = PlacticMonoid(4)
            sage: P.rank()
            4
            sage: TestSuite(PlacticMonoid(2)).run()

        TESTS::

            sage: PlacticMonoid(4) is PlacticMonoid(ZZ(4))
            True
            sage: PlacticMonoid(-1)
            Traceback (most recent call last):
            ...
            ValueError: the rank must be a positive integer
        """
        from sage.categories.monoids import Monoids
        if not isinstance(n, (int, Integer)):
            raise ValueError("the rank must be a positive integer")
        n = ZZ(n)
        if n <= 0:
            raise ValueError("the rank must be a positive integer")
        self._n = n
        Parent.__init__(self, category=(Monoids().FinitelyGenerated().Infinite(),
                                        SetsWithGrading().Infinite()))

    def rank(self):
        """
        Return the rank of ``self``.

        EXAMPLES::

            sage: PlacticMonoid(4).rank()
            4
            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: HypoplacticMonoid(4).rank()
            4
        """
        return self._n

    @cached_method
    def monoid_generators(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: M = PlacticMonoid(4)
            sage: G = M.monoid_generators()
            sage: G[1], G[2], G[3], G[4]
            (1, 2, 3, 4)
            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: G = H.monoid_generators()
            sage: G
            Finite family {1: 1, 2: 2, 3: 3, 4: 4}
            sage: G[1], G[2], G[3], G[4]
            (1, 2, 3, 4)
        """
        from sage.sets.family import Family
        return Family({i: self.element_class(self, (i,))
                       for i in range(1, self._n + 1)})

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: M = PlacticMonoid(3)
            sage: M.one() == M([])
            True
            sage: len(M.one())
            0
            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(3)
            sage: H.one() == H([])
            True
            sage: len(H.one())
            0
        """
        return self.element_class(self, ())

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: M = PlacticMonoid(3)
            sage: M.an_element()
            1
            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(3)
            sage: H.an_element()
            1
        """
        return self.monoid_generators()[1]


class WordMonoidElement(ElementWrapper):
    r"""
    An element of a word monoid.

    Elements are represented by words in the alphabet
    `\{1, 2, \ldots, n\}`.

    EXAMPLES::

        sage: M = PlacticMonoid(4)
        sage: M([2, 1, 3])
        213

        sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
        sage: H = HypoplacticMonoid(4)
        sage: x = H([3, 2, 2, 1])
        sage: x
        3221
        sage: parent(x)
        Hypoplactic monoid of rank 4
    """
    def __init__(self, parent, value):
        """
        Initialize ``self``.

        INPUT:

        - ``parent`` -- word monoid
        - ``value`` -- word given as a list or tuple of letters in the
          alphabet of ``parent``

        TESTS::

            sage: M = PlacticMonoid(4)
            sage: M([1, 2, 4])
            124
            sage: M([1, 2.5])
            Traceback (most recent call last):
            ...
            ValueError: letters must be integers from 1 to 4

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: x = H([3, 2, 2, 1]); x
            3221
            sage: H([1, 2.5])
            Traceback (most recent call last):
            ...
            ValueError: letters must be integers from 1 to 4
        """
        r = parent.rank()
        try:
            value = tuple(map(ZZ, value))
        except TypeError:
            raise ValueError("letters must be integers from 1 to %s" % r)
        if not all(1 <= i <= r for i in value):
            raise ValueError("letters must be integers from 1 to %s" % r)
        ElementWrapper.__init__(self, parent, value)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = PlacticMonoid(4)
            sage: M([2, 1, 3])
            213
            sage: M([2, 3, 1])
            231
            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: H([2, 1, 3])
            213
        """
        if not self.value:
            return ''
        return ''.join(str(x) for x in self.value)

    def __len__(self):
        """
        Return the length of ``self`` as a word.

        This is also the grade of ``self``.

        EXAMPLES::

            sage: M = PlacticMonoid(4)
            sage: len(M([3, 1, 2]))
            3
            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: len(H([2, 1, 3]))
            3
        """
        return len(self.value)

    grade = __len__

    def __hash__(self):
        """
        TESTS::

            sage: M = PlacticMonoid(4)
            sage: x = M([3, 1, 2])
            sage: hash(x) == hash(M([3,1,2]))
            True
            sage: y = M([1, 3, 2])
            sage: hash(x) == hash(y)
            True
            sage: y = M([1, 2, 3])
            sage: hash(x)==hash(y)
            False

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: x = H([2, 3, 1, 2])
            sage: hash(x) == hash(H([2, 3, 1, 2]))
            True
        """
        return hash(self.to_tableau())

    def __iter__(self):
        """
        Iterate over the letters of ``self``.

        EXAMPLES::

            sage: M = PlacticMonoid(4)
            sage: list(M([3, 1, 2]))
            [3, 1, 2]

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: list(H([3, 2, 2, 1]))
            [3, 2, 2, 1]
        """
        return iter(self.value)

    def _mul_(self, other):
        """
        Multiply ``self`` by ``other``.

        Multiplication is induced by concatenation of words. The
        concatenated word is inserted using the method `to_tableau`,
        and the product is stored using the resulting representative
        using the word method.

        INPUT:

        - ``other`` -- an element of the hypoplactic monoid

        OUTPUT:

        The product of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: a = H([3])
            sage: b = H([4])
            sage: a * b
            34

            sage: c = H([4])
            sage: d = H([3])
            sage: c * d
            43

            sage: M = PlacticMonoid(4)
            sage: a = M([2, 1]); b = M([3, 2])
            sage: a * b
            2132
            sage: (a * b).to_word()
            2312

        TESTS::

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: H([]) * H([3, 2, 2, 1]) == H([3, 2, 2, 1])
            True
            sage: H([3, 2, 2, 1]) * H([]) == H([3, 2, 2, 1])
            True
        """
        parent = self.parent()
        word = self.value + other.value
        return self.__class__(parent, word)

    def __eq__(self, other):
        """
        Return whether ``self`` and ``other`` are equal.

        EXAMPLES::

            sage: M = PlacticMonoid(4)
            sage: M([2, 1, 3]) == M([2, 3, 1])
            True
            sage: M([2, 1, 3]) == M([3, 2, 1])
            False
            sage: M3 = PlacticMonoid(3)
            sage: M([2, 1, 3]) == M3([2, 1, 3])
            False

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: H([3, 2, 2, 1]) == H([2, 3, 1, 2])
            True
            sage: H([3, 2, 2, 1]) == H([1, 2, 2, 3])
            False
            sage: H3 = HypoplacticMonoid(3)
            sage: H([2, 1, 3]) == H3([2, 1, 3])
            False
        """
        return (isinstance(other, self.parent().Element)
                and self.parent() == other.parent()
                and self.to_tableau() == other.to_tableau())

    def shape(self):
        """
        Return the shape of the insertion tableau of ``self``.

        EXAMPLES::

            sage: M = PlacticMonoid(4)
            sage: M([2, 1, 3]).shape()
            [2, 1]
            sage: M([]).shape()
            []
        """
        return self.to_tableau().shape()

    def is_canonical(self):
        """
        Return whether ``self`` is its row reading word representative.

        EXAMPLES::

            sage: M = PlacticMonoid(3)
            sage: M([3, 2, 1]).is_canonical()
            True
            sage: M([1, 3, 2]).is_canonical()
            False

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: H([3, 2, 2, 1]).is_canonical()
            False
            sage: H([2,1,3,2]).is_canonical()
            True
        """
        return self.value == self.to_word().value

    def equivalence_class(self):
        r"""
        Return the equivalence class of ``self``.

        This is the list of all words with the same insertion tableau as ``self``.

        EXAMPLES::

            sage: from sage.monoids.hypoplactic_monoid import HypoplacticMonoid
            sage: H = HypoplacticMonoid(3)
            sage: H([2, 1, 3]).equivalence_class()
            [213, 231]
            sage: H = HypoplacticMonoid(4)
            sage: H([3, 1, 4, 2]).equivalence_class()
            [3142, 3124, 3412, 1342, 1324]
            sage: H([3, 1, 1, 2]).equivalence_class()
            [3112, 1312, 1132]
        """
        parent = self.parent()
        tab = self.to_tableau()
        return [m for w in Permutations(self.value)
                if (m := parent(w)).to_tableau() == tab]


class PlacticMonoid(WordMonoid):
    r"""
    The plactic monoid on the alphabet `\{1, 2, \ldots, n\}`.

    INPUT:

    - ``n`` -- a positive integer; the size of the alphabet

    Elements are represented by words in `\{1, 2, \ldots, n\}`. Equality is
    determined by comparing RSK insertion tableaux. Multiplication is induced
    by concatenation of words, and the identity is the empty word.

    EXAMPLES::

        sage: M = PlacticMonoid(4)
        sage: M
        Plactic monoid of rank 4
        sage: M.rank()
        4
        sage: M([2, 1, 3]).to_tableau()
        [[1, 3], [2]]
        sage: M([2, 1, 3]) == M([2, 3, 1])
        True
        sage: M([2, 1]) * M([3, 2])
        2132
        sage: (M([2, 1]) * M([3, 2])).to_word()
        2312

    TESTS::

        sage: M = PlacticMonoid(4)
        sage: M([]) == M.one()
        True
        sage: PlacticMonoid(0)
        Traceback (most recent call last):
        ...
        ValueError: the rank must be a positive integer
        sage: PlacticMonoid(4.1)
        Traceback (most recent call last):
        ...
        ValueError: the rank must be a positive integer
        sage: M([1, 2, 5])
        Traceback (most recent call last):
        ...
        ValueError: letters must be integers from 1 to 4
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlacticMonoid(4)
            Plactic monoid of rank 4
        """
        return f"Plactic monoid of rank {self._n}"

    def subset(self, k):
        r"""
        Return the plactic monoid elements represented by words of length ``k``.

        Since the plactic monoid is infinite, this returns the finite set of
        elements of a fixed size, using their row reading word representatives.

        EXAMPLES::

            sage: M = PlacticMonoid(2)
            sage: M.subset(1)
            Lazy family (to_word(i))_{i in Semistandard tableaux of size 1 and maximum entry 2}
            sage: list(M.subset(1))
            [1, 2]
            sage: M.subset(2)
            Lazy family (to_word(i))_{i in Semistandard tableaux of size 2 and maximum entry 2}
            sage: list(M.subset(2))
            [11, 12, 22, 21]
        """
        if not isinstance(k, (int, Integer)) or k < 0:
            raise ValueError("the size must be a nonnegative integer")

        # Plactic monoid elements correspond to semistandard tableaux.
        # For each partition shape of k, Sage generates all tableaux of that
        # shape with entries bounded by the rank.
        def to_word(t):
            return self(t.to_word())
        tableaux = SemistandardTableaux(k, max_entry=self.rank())
        return Family(tableaux, to_word, lazy=True)

    class Element(WordMonoidElement):
        r"""
        An element of a plactic monoid, represented by a word.

        EXAMPLES::

            sage: M = PlacticMonoid(4)
            sage: M([2, 1, 3])
            213
        """
        def to_word(self):
            """
            Return the row reading word representative of ``self``.

            EXAMPLES::

                sage: M = PlacticMonoid(4)
                sage: M([2, 3, 1]).to_word()
                213
            """
            tableau_flipped = list(reversed(self.to_tableau()))
            row_list = list(chain.from_iterable(tableau_flipped))
            parent = self.parent()
            return parent(row_list)

        def to_tableau(self):
            """
            Return the RSK insertion tableau corresponding to ``self``.

            EXAMPLES::

                sage: M = PlacticMonoid(4)
                sage: M([1, 3, 2]).to_tableau()
                [[1, 2], [3]]
                sage: M([]).to_tableau()
                []
            """
            return RSK(self.value)[0]

        def equivalence_class(self):
            r"""
            Return the plactic equivalence class of ``self``.

            This is the list of all words with the same RSK insertion tableau
            as ``self``.

            EXAMPLES::

                sage: M = PlacticMonoid(3)
                sage: M([2, 1, 3]).equivalence_class()
                [213, 231]
            """
            parent = self.parent()
            P = self.to_tableau()
            shape = P.shape()
            return [parent(RSK_inverse(P, Q)[1]) for Q in StandardTableaux(shape)]
