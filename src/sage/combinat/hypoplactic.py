r"""
Hypoplactic monoid

AUTHORS:

- Pranita Urlam (2026): initial version

The *hypoplactic monoid* `\mathrm{Hypo}(A)` on a totally ordered finite
alphabet `A = \{1 < 2 < \cdots < n\}` is the quotient of the free monoid
`A^*` by the *hypoplactic congruence*: two words `u, v \in A^*` are
*hypoplactically congruent* if and only if they produce the same
quasi-ribbon insertion tableau under the Krob-Thibon (KT) insertion
algorithm.

The KT insertion algorithm is analogous to the RSK correspondence for the
plactic monoid, but uses quasi-ribbon tableaux instead of Young tableaux.

REFERENCES:

- [KT1997]_ D. Krob and J.-Y. Thibon, *Noncommutative symmetric functions IV:
  Quantum linear groups and Hecke algebras at q = 0*, J. Algebraic Combin. **6**
  (1997), 339--376.

- [Nov2000]_ J.-C. Novelli, *On the hypoplactic monoid*, Discrete Math. **217**
  (2000), 315--336.

EXAMPLES::

    sage: from sage.combinat.hypoplactic import HypoplacticMonoid, krob_thibon_insertion
    sage: H = HypoplacticMonoid(4)
    sage: H
    Hypoplactic monoid on {1, 2, 3, 4}
    sage: a = H([2, 1, 3]); b = H([3, 1]); a * b
    [2, 1, 3, 3, 1]
    sage: a == H([2, 3, 1])   # 213 ~ 231 by a Knuth relation (yxz ~ yzx, y=2, x=1, z=3)
    True
    sage: P, Q = krob_thibon_insertion([3, 1, 3, 2])
    sage: P
    [[1, 2], [3, 3]]
    sage: Q
    [[2, 4], [1, 3]]
"""
# ****************************************************************************
#       Copyright (C) 2026 Pranita Urlam
#
#  Distributed under the terms of the GNU General Public License (GPL)
#              https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau, _kt_insert
from sage.structure.parent import Parent
from sage.structure.element import MonoidElement
from sage.categories.monoids import Monoids
from sage.structure.unique_representation import UniqueRepresentation


# ---------------------------------------------------------------------------
# KT insertion
# ---------------------------------------------------------------------------

def krob_thibon_insertion(word):
    r"""
    Perform the Krob-Thibon (KT) insertion on ``word``.

    Given a word `w = w_1 w_2 \cdots w_n` over a totally ordered alphabet,
    this algorithm returns a pair `(P, Q)` where:

    - `P` is the *quasi-ribbon insertion tableau* (a
      :class:`~sage.combinat.quasi_ribbon_tableau.QuasiRibbonTableau`).
    - `Q` is the *recording tableau*, stored as a list of lists with the same
      shape as `P`. Entry `i` of `Q` is placed in the cell added when `w_i`
      was inserted.

    Two words are hypoplactically congruent if and only if they yield the
    same insertion tableau `P`.

    INPUT:

    - ``word`` -- a list (or iterable) of positive integers

    OUTPUT:

    A pair ``(P, Q)`` where ``P`` is a
    :class:`~sage.combinat.quasi_ribbon_tableau.QuasiRibbonTableau` and
    ``Q`` is a list of lists of integers.

    EXAMPLES::

        sage: from sage.combinat.hypoplactic import krob_thibon_insertion
        sage: krob_thibon_insertion([])
        ([], [])
        sage: krob_thibon_insertion([3])
        ([[3]], [[1]])
        sage: krob_thibon_insertion([3, 1, 3, 2])
        ([[1, 2], [3, 3]], [[2, 4], [1, 3]])

    The insertion tableau is the same for hypoplactically equivalent words
    (e.g. ``[2,1,3]`` and ``[2,3,1]`` are related by the Knuth move
    `yxz \sim yzx` with `x=1 < y=2 \leq z=3`)::

        sage: P1, _ = krob_thibon_insertion([2, 1, 3])
        sage: P2, _ = krob_thibon_insertion([2, 3, 1])
        sage: P1 == P2
        True

    An increasing word produces a single row::

        sage: krob_thibon_insertion([1, 2, 3, 4])[0]
        [[1, 2, 3, 4]]

    A decreasing word produces a single column::

        sage: krob_thibon_insertion([4, 3, 2, 1])[0]
        [[1], [2], [3], [4]]

    TESTS::

        sage: from sage.combinat.hypoplactic import krob_thibon_insertion
        sage: P, Q = krob_thibon_insertion([1, 3, 2, 4])
        sage: P
        [[1, 2], [3, 4]]
        sage: Q
        [[1, 3], [2, 4]]
        sage: sum(len(r) for r in Q) == 4
        True
    """
    word = list(word)
    if not word:
        return QuasiRibbonTableau([]), []

    # Build P and Q simultaneously.  Q tracks which cell is added at each step.
    p_rows = []   # rows of P, bottom to top
    q_rows = []   # rows of Q, bottom to top (same shape as P)

    for step, x in enumerate(word, start=1):
        old_shape = [len(r) for r in p_rows]
        p_rows = _kt_insert(p_rows, x)
        new_shape = [len(r) for r in p_rows]

        # Find which cell was added: the unique position where new_shape > old_shape.
        # The new rows list may have one more row than the old one (new bottom row)
        # or an existing row may have grown by one cell.
        if len(new_shape) > len(old_shape):
            # A new row was added at the bottom.
            # Shift Q rows up by one and prepend a new row [step].
            q_rows = [[step]] + q_rows
        else:
            # An existing row grew by one cell.
            for i, (old_len, new_len) in enumerate(zip(old_shape, new_shape)):
                if new_len > old_len:
                    q_rows[i] = q_rows[i] + [step]
                    break

    return QuasiRibbonTableau(p_rows), q_rows


# ---------------------------------------------------------------------------
# HypoplacticMonoid
# ---------------------------------------------------------------------------

class HypoplacticMonoidElement(MonoidElement):
    r"""
    An element of a :class:`HypoplacticMonoid`.

    Elements are represented by words over `\{1, \ldots, n\}`.  Two
    elements are equal if and only if their KT insertion tableaux agree.

    EXAMPLES::

        sage: from sage.combinat.hypoplactic import HypoplacticMonoid
        sage: H = HypoplacticMonoid(3)
        sage: a = H([1, 2, 3]); a
        [1, 2, 3]
        sage: b = H([2, 1, 3]); b
        [2, 1, 3]
        sage: a == b  # different insertion tableaux
        False
        sage: H([2, 1, 3]) == H([2, 3, 1])  # same insertion tableau (213 ~ 231)
        True
    """

    def __init__(self, parent, word):
        r"""
        Initialise ``self``.

        TESTS::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: H = HypoplacticMonoid(3)
            sage: H([1, 2])
            [1, 2]
        """
        MonoidElement.__init__(self, parent)
        self._word = tuple(word)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(3)([1, 2, 3])
            [1, 2, 3]
        """
        return repr(list(self._word))

    def word(self):
        r"""
        Return the underlying word as a tuple.

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(4)([2, 1, 3]).word()
            (2, 1, 3)
        """
        return self._word

    def insertion_tableau(self):
        r"""
        Return the KT insertion tableau of ``self``.

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(4)([2, 1, 3]).insertion_tableau()
            [[1], [2, 3]]
            sage: HypoplacticMonoid(4)([2, 3, 1]).insertion_tableau()
            [[1], [2, 3]]
        """
        P, _ = krob_thibon_insertion(self._word)
        return P

    def _mul_(self, other):
        r"""
        Return the product of ``self`` and ``other`` (concatenation of words).

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: H([1, 2]) * H([2, 1])
            [1, 2, 2, 1]
            sage: H([1]) * H.one()
            [1]
        """
        return self.__class__(self.parent(), self._word + other._word)

    def __eq__(self, other):
        r"""
        Two elements are equal iff their insertion tableaux agree.

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: H([2, 1, 3]) == H([2, 3, 1])   # 213 ~ 231 by Knuth
            True
            sage: H([1, 2, 3]) == H([3, 2, 1])
            False
        """
        if not isinstance(other, HypoplacticMonoidElement):
            return False
        return self.insertion_tableau() == other.insertion_tableau()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        r"""
        Hash via the insertion tableau's reading word.

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: hash(H([2, 1, 3])) == hash(H([2, 3, 1]))
            True
        """
        P = self.insertion_tableau()
        return hash(tuple(P.reading_word()))


class HypoplacticMonoid(UniqueRepresentation, Parent):
    r"""
    The hypoplactic monoid on a totally ordered alphabet `\{1, \ldots, n\}`.

    The hypoplactic monoid `\mathrm{Hypo}(n)` is the quotient of the free
    monoid `\{1, \ldots, n\}^*` by the hypoplactic congruence: two words are
    congruent iff they have the same quasi-ribbon insertion tableau under the
    Krob-Thibon (KT) insertion algorithm.

    Multiplication is by concatenation of words.

    INPUT:

    - ``n`` -- positive integer; the size of the alphabet

    EXAMPLES::

        sage: from sage.combinat.hypoplactic import HypoplacticMonoid
        sage: H = HypoplacticMonoid(4)
        sage: H
        Hypoplactic monoid on {1, 2, 3, 4}
        sage: H.alphabet()
        (1, 2, 3, 4)
        sage: a = H([1, 2]); b = H([3, 2]); a * b
        [1, 2, 3, 2]
        sage: H.one()
        []

    Knuth-equivalent words are equal (``[2,1,3]`` and ``[2,3,1]`` are related
    by `yxz \sim yzx` with `x=1 < y=2 \leq z=3`)::

        sage: H([2, 1, 3]) == H([2, 3, 1])
        True

    The identity element is the empty word::

        sage: H([1, 2, 3]) * H([]) == H([1, 2, 3])
        True
        sage: H([]) * H([1, 2, 3]) == H([1, 2, 3])
        True

    TESTS::

        sage: from sage.combinat.hypoplactic import HypoplacticMonoid
        sage: H = HypoplacticMonoid(3)
        sage: TestSuite(H).run(skip=['_test_elements'])
    """

    Element = HypoplacticMonoidElement

    def __init__(self, n):
        r"""
        Initialise ``self``.

        TESTS::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(3)
            Hypoplactic monoid on {1, 2, 3}
        """
        if n < 1:
            raise ValueError("n must be a positive integer")
        self._n = n
        Parent.__init__(self, category=Monoids())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(4)
            Hypoplactic monoid on {1, 2, 3, 4}
        """
        return f"Hypoplactic monoid on {{{', '.join(str(i) for i in range(1, self._n + 1))}}}"

    def alphabet(self):
        r"""
        Return the alphabet as a tuple `(1, 2, \ldots, n)`.

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(4).alphabet()
            (1, 2, 3, 4)
        """
        return tuple(range(1, self._n + 1))

    def _element_constructor_(self, word):
        r"""
        Construct an element from a word (list of letters).

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: H = HypoplacticMonoid(4)
            sage: H([1, 2, 3])
            [1, 2, 3]
            sage: H([])
            []
        """
        word = list(word)
        if any(not (1 <= x <= self._n) for x in word):
            raise ValueError(f"all letters must be in {{1, ..., {self._n}}}")
        return self.element_class(self, word)

    def one(self):
        r"""
        Return the identity element (the empty word).

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(3).one()
            []
        """
        return self.element_class(self, [])

    def an_element(self):
        r"""
        Return a sample element.

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(3).an_element()
            [1]
        """
        return self.element_class(self, [1])

    def generators(self):
        r"""
        Return the generators of ``self``, one for each letter.

        EXAMPLES::

            sage: from sage.combinat.hypoplactic import HypoplacticMonoid
            sage: HypoplacticMonoid(3).generators()
            ([1], [2], [3])
        """
        return tuple(self.element_class(self, [i]) for i in range(1, self._n + 1))
