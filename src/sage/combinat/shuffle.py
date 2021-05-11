# -*- coding: utf-8 -*-
r"""
Shuffle product of iterables

The shuffle product of two sequences of lengths `m` and `n` is a
sum over the `\binom{m+n}{n}` ways of interleaving the two sequences.

That could be defined inductively by:

.. MATH::

    (a_n)_{n \geqslant 0} \Cup (b_m)_{m \geqslant 0} =
    a_0 \cdot \left((a_n)_{n \geqslant 1} \Cup (b_m)_{m \geqslant 0}\right)
    + b_0 \cdot \left((a_n)_{n \geqslant 0} \Cup (b_m)_{m \geqslant 1}\right)

with `(a_n)` and `(b_m)` two non-empty sequences and if one of them is empty
then the product is equals to the other.

The shuffle product has been introduced by S. Eilenberg and S. Mac Lane in
1953 [EilLan53]_.

EXAMPLES::

    sage: from sage.combinat.shuffle import ShuffleProduct
    sage: list(ShuffleProduct([1,2], ["a", "b", "c"]))
    [[1, 2, 'a', 'b', 'c'],
     ['a', 1, 2, 'b', 'c'],
     [1, 'a', 2, 'b', 'c'],
     ['a', 'b', 1, 2, 'c'],
     ['a', 1, 'b', 2, 'c'],
     [1, 'a', 'b', 2, 'c'],
     ['a', 'b', 'c', 1, 2],
     ['a', 'b', 1, 'c', 2],
     ['a', 1, 'b', 'c', 2],
     [1, 'a', 'b', 'c', 2]]

References:

.. [EilLan53] On the groups `H(\pi, n)`, I,
    Samuel Eilenberg and
    Saunders Mac Lane,
    1953.

Author:

- Jean-Baptiste Priez
"""
# ****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Iterable
import itertools
import operator

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.words.abstract_word import Word_class
from sage.rings.integer import Integer
from sage.structure.parent import Parent


# TODO: Think about Parent/Element for this and the category
# sage.categories.finite_enumerated_sets.FiniteEnumeratedSets

class ShuffleProduct_abstract(Parent):
    """
    Abstract base class for shuffle products.
    """
    def __init__(self, l1, l2, element_constructor=None):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: SP = ShuffleProduct([1,2],[4,5,7,8,9])
            sage: TestSuite(SP).run(skip="_test_an_element")
        """
        self._l1 = l1
        self._l2 = l2
        if element_constructor is None:
            self._element_constructor_ = self._constructor_
        else:
            self._element_constructor_ = element_constructor
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _constructor_(self, x):
        """
        This is used by default to produce elements when iterating.

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: S = ShuffleProduct([1,2,3],[4,5])
            sage: S._constructor_((1,2,3))
            [1, 2, 3]

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: S = ShuffleProduct_overlapping(w,u)
            sage: S._constructor_((1,2,3))
            [1, 2, 3]

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: S = SetShuffleProduct([[1,2],[3,4]], [[1,4]])
            sage: S._constructor_((1,2,3))
            [1, 2, 3]

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping_r
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: S = ShuffleProduct_overlapping_r(w,u,1)
            sage: S._constructor_((1,2,3))
            [1, 2, 3]
        """
        return list(x)

    def __eq__(self, other):
        """
        Test for equality.

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: SP = ShuffleProduct([1,2],[4,5,7,8,9])
            sage: loads(dumps(SP)) == SP
            True
            sage: SP == ShuffleProduct([1,2],[4,5,7])
            False
        """
        if not isinstance(other, type(self)):
            return False
        return (self._l1 == other._l1 and self._l2 == other._l2)

    def __ne__(self, other):
        """
        Test for unequality

        EXAMPLES::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping_r
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: A = ShuffleProduct_overlapping_r(w,u,1)
            sage: B = ShuffleProduct_overlapping_r(u,w,2)
            sage: A != A
            False
            sage: A != B
            True
        """
        return not (self == other)

    def __contains__(self, x):
        """
        Check containment.

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: SP = ShuffleProduct([1,2],[4,5,7])
            sage: [1,4,5,2,7] in SP
            True
            sage: (1,4,5,2,7) in SP
            False
            sage: [2,1,4,5,7] in SP
            False
        """
        return x in self.list()


class SetShuffleProduct(ShuffleProduct_abstract):
    """
    The union of all possible shuffle products of two sets of iterables.

    EXAMPLES::

        sage: from sage.combinat.shuffle import SetShuffleProduct
        sage: sorted(SetShuffleProduct({(1,), (2,3)}, {(4,5), (6,)}))
        [[1, 4, 5],
         [1, 6],
         [2, 3, 4, 5],
         [2, 3, 6],
         [2, 4, 3, 5],
         [2, 4, 5, 3],
         [2, 6, 3],
         [4, 1, 5],
         [4, 2, 3, 5],
         [4, 2, 5, 3],
         [4, 5, 1],
         [4, 5, 2, 3],
         [6, 1],
         [6, 2, 3]]
    """
    def __init__(self, l1, l2, element_constructor=None):
        """
        Construct the set of all possible shuffle products of two sets of iterables.

        INPUT:

        - ``l1``, ``l2`` -- iterable: the sets to shuffle

        - ``element_constructor`` --  constructor for the returned elements

        TESTS::

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: X = SetShuffleProduct({(1,2,3), (2,3,4)}, {(5,)})
            sage: X   # random
            Shuffle set product of: [(2, 3, 4), (1, 2, 3)] and [(5,)]
            sage: TestSuite(X).run(skip="_test_an_element")

            sage: list(SetShuffleProduct({(1,2,3), (2,3,4)}, {(5,)}))   # random
            [[2, 3, 4, 5], [2, 5, 3, 4], [5, 2, 3, 4], [2, 3, 5, 4],
             [1, 2, 3, 5], [1, 5, 2, 3], [5, 1, 2, 3], [1, 2, 5, 3]]
        """
        assert(isinstance(l1, Iterable) and
               isinstance(l2, Iterable))
        assert(all(isinstance(elem, Iterable) for elem in l1))
        assert(all(isinstance(elem, Iterable) for elem in l2))

        if element_constructor is None:
            try:
                e = next(iter(l1))
                try:
                    element_constructor = e.parent()._element_constructor_
                except AttributeError:
                    pass
            except StopIteration:
                pass
        ShuffleProduct_abstract.__init__(self, list(l1), list(l2), element_constructor)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: SetShuffleProduct([[1,2],[3,4]], [[1,4]])
            Shuffle set product of: [[1, 2], [3, 4]] and [[1, 4]]
            sage: SetShuffleProduct([()], [[1,4]])
            Shuffle set product of: [()] and [[1, 4]]

        """
        return "Shuffle set product of: %s and %s" % (self._element_constructor_(self._l1),
                                                      self._element_constructor_(self._l2))

    def _ascii_art_(self):
        r"""
        TESTS::

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: ascii_art(SetShuffleProduct([[BinaryTree()], [BinaryTree([]), BinaryTree([[],[]])]],
            ....: [[1,4]]))
            Set shuffle product of:
            [       [ o,   o   ] ]
            [       [     / \  ] ]
            [ [  ], [    o   o ] ] and [ [ 1, 4 ] ]

        """
        from sage.typeset.ascii_art import ascii_art
        return (ascii_art("Set shuffle product of:") *
                (ascii_art(self._l1) + ascii_art(" and ") +
                 ascii_art(self._l2)))

    def __iter__(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: list(SetShuffleProduct([[],[]], [[]]))
            [[], []]
            sage: list(SetShuffleProduct([[1,2],[3]], [[4]]))
            [[1, 2, 4], [4, 1, 2], [1, 4, 2], [3, 4], [4, 3]]
            sage: list(SetShuffleProduct([[1,2],[3,4]], [[1,4]], element_constructor=set))
            [{1, 2, 4},
             {1, 2, 4},
             {1, 2, 4},
             {1, 2, 4},
             {1, 2, 4},
             {1, 2, 4},
             {1, 3, 4},
             {1, 3, 4},
             {1, 3, 4},
             {1, 3, 4},
             {1, 3, 4},
             {1, 3, 4}]
        """
        return itertools.chain.from_iterable(
            ShuffleProduct(*pair,
                           element_constructor=self._element_constructor_)
            for pair in itertools.product(self._l1, self._l2))

    def cardinality(self):
        """
        The cardinality is defined by the sum of the cardinality of all shuffles.
        That means by a sum of binomials.

        TESTS::

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: SetShuffleProduct([[1,2],[3,4]], [[1,4]], element_constructor=set).cardinality()
            12
        """
        def comp_binom(el1, el2):
            ll1 = Integer(len(el1))
            ll2 = Integer(len(el2))
            return (ll1 + ll2).binomial(ll2)

        return sum(comp_binom(el1, el2)
                for (el1, el2) in itertools.product(self._l1, self._l2))


class ShuffleProduct(ShuffleProduct_abstract):
    """
    Shuffle product of two iterables.

    EXAMPLES::

        sage: from sage.combinat.shuffle import ShuffleProduct
        sage: list(ShuffleProduct("abc", "de", element_constructor="".join))
        ['abcde',
         'adbce',
         'dabce',
         'abdce',
         'adebc',
         'daebc',
         'deabc',
         'adbec',
         'dabec',
         'abdec']
        sage: list(ShuffleProduct("", "de", element_constructor="".join))
        ['de']

    """
    def __init__(self, l1, l2, element_constructor=None):
        """
        Construct the shuffle product of two iterable.

        INPUT:

        - ``l1``, ``l2`` -- iterable: iterables to shuffle

        - ``element_constructor``:  constructor for the returned elements

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: SP = ShuffleProduct([1,2,3],[4,5])
            sage: SP
            Shuffle product of: [1, 2, 3] and [4, 5]
            sage: TestSuite(SP).run(skip="_test_an_element")

            sage: list(ShuffleProduct(Word("aa"), Word("bbb"), Word))
            [word: aabbb, word: baabb, word: ababb, word: bbaab, word: babab, word: abbab,
             word: bbbaa, word: bbaba, word: babba, word: abbba]
        """
        assert(isinstance(l1, Iterable) and isinstance(l2, Iterable))

        if element_constructor is None:
            try:
                element_constructor = l1.parent()._element_constructor_
            except AttributeError:
                pass
        ShuffleProduct_abstract.__init__(self, list(l1), list(l2), element_constructor)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: ShuffleProduct([1,2,3],[4,5])
            Shuffle product of: [1, 2, 3] and [4, 5]
            sage: B = BinaryTree
            sage: ShuffleProduct([B(), B([[],[]])], [])
            Shuffle product of: [., [[., .], [., .]]] and []
        """
        return "Shuffle product of: %s and %s" % (self._l1, self._l2)

    def _ascii_art_(self):
        r"""
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: ascii_art(ShuffleProduct([1,2,3],[4,5]))
            Shuffle product of:
            [ 1, 2, 3 ] and [ 4, 5 ]
            sage: B = BinaryTree
            sage: ascii_art(ShuffleProduct([B([]), B([[],[]])],
            ....:   [B([[[],[]],[[],None]])]))
            Shuffle product of:
                             [     __o__   ]
                             [    /     \  ]
            [ o,   o   ]     [   o       o ]
            [     / \  ]     [  / \     /  ]
            [    o   o ] and [ o   o   o   ]
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art("Shuffle product of:") * \
            (ascii_art(self._l1) + ascii_art(" and ") +
             ascii_art(self._l2))

    def __iter__(self):
        r"""
        Efficient iteration from a gray code on binary words in `B(n,k)`.

        (with `B(n,k)` the number of binary words of size `n` with `k` one.

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: list(ShuffleProduct([1,2,3],[4,5]))
            [[1, 2, 3, 4, 5], [1, 4, 2, 3, 5], [4, 1, 2, 3, 5], [1, 2, 4, 3, 5], [1, 4, 5, 2, 3],
             [4, 1, 5, 2, 3], [4, 5, 1, 2, 3], [1, 4, 2, 5, 3], [4, 1, 2, 5, 3], [1, 2, 4, 5, 3]]
            sage: B = BinaryTree
            sage: ascii_art(list(ShuffleProduct([B([]), B([[],[]])],
            ....:   [B([[[],[]],[[],None]])])))
            [ [ o,   o  ,     __o__   ]  [     __o__  , o,   o   ]
            [ [     / \      /     \  ]  [    /     \       / \  ]
            [ [    o   o    o       o ]  [   o       o     o   o ]
            [ [            / \     /  ]  [  / \     /            ]
            [ [           o   o   o   ], [ o   o   o             ],
            <BLANKLINE>
             [ o,     __o__  ,   o   ] ]
             [       /     \    / \  ] ]
             [      o       o  o   o ] ]
             [     / \     /         ] ]
             [    o   o   o          ] ]
        """

        # -------------- Gray code --------------
        def swap(i, j):
            l[i - 1], l[j - 1] = l[j - 1], l[i - 1]

        def gen(n, k):
            if 0 < k < n:
                for _ in gen(n - 1, k):
                    yield

                if k == 1:
                    swap(n, n - 1)
                else:
                    swap(n, k - 1)
                yield

                for _ in neg(n - 1, k - 1):
                    yield

        def neg(n, k):
            if 0 < k < n:
                for _ in gen(n - 1, k - 1):
                    yield

                if k == 1:
                    swap(n, n - 1)
                else:
                    swap(n, k - 1)
                yield

                for _ in neg(n - 1, k):
                    yield

        ####################################

        m = len(self._l1)
        n = len(self._l2)
        mn = m + n
        l = [0] * m + [1] * n  # [0, 0 ... m times, 1, 1, 1 ... n times]

        EC = self._element_constructor_
        yield EC(self._l1 + self._l2)

        for _ in gen(mn, m):
            l1 = iter(self._l1)
            l2 = iter(self._l2)
            yield EC([next(l2) if l[k] else next(l1) for k in range(mn)])

    def __contains__(self, iterable):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: sh = ShuffleProduct([1,2,3],[4,5])
            sage: list(range(1,6)) in sh
            True
            sage: list(range(1,7)) in sh
            False
            sage: [3,4,5,1,2] in sh
            False
            sage: [1,4,2,5,3] in sh
            True
            sage: [1,4,2,2,5,3] in sh
            False
        """
        if not isinstance(iterable, type(self._element_constructor_([]))):
            return False

        l1 = self._l1
        l2 = self._l2
        len_l1 = len(l1)
        len_l2 = len(l2)
        i_l1 = i_l2 = 0
        iterable = list(iterable)

        for i, el in enumerate(iterable):
            if l1[i_l1] == el:
                i_l1 += 1
            elif l2[i_l2] == el:
                i_l2 += 1
            else:
                return False
            if i_l1 == len_l1:
                return iterable[i + 1:] == l2[i_l2:]
            if i_l2 == len_l2:
                return iterable[i + 1:] == l1[i_l1:]
        return (i_l1 + 1 == len_l1) and (i_l2 + 1 == len_l2)

    def cardinality(self):
        r"""
        Return the number of shuffles of `l_1` and `l_2`, respectively of
        lengths `m` and `n`, which is `\binom{m+n}{n}`.

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: ShuffleProduct([3,1,2], [4,2,1,3]).cardinality()
            35
            sage: ShuffleProduct([3,1,2,5,6,4], [4,2,1,3]).cardinality() == binomial(10,4)
            True
        """
        ll1 = Integer(len(self._l1))
        ll2 = Integer(len(self._l2))
        return (ll1 + ll2).binomial(ll1)


class ShuffleProduct_overlapping_r(ShuffleProduct_abstract):
    """
    The overlapping shuffle product of the two words ``w1`` and ``w2``
    with precisely ``r`` overlaps.

    See :class:`ShuffleProduct_overlapping` for a definition.

    EXAMPLES::

        sage: from sage.combinat.shuffle import ShuffleProduct_overlapping_r
        sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
        sage: S = ShuffleProduct_overlapping_r(w,u,1)
        sage: list(S)
        [word: 11,9,1,
         word: 2,18,1,
         word: 11,1,9,
         word: 2,9,10,
         word: 939,
         word: 9,2,10]
    """
    def __init__(self, w1, w2, r, element_constructor=None, add=operator.add):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping_r
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: S = ShuffleProduct_overlapping_r(w,u,1)
            sage: TestSuite(S).run(skip="_test_an_element")
        """
        self.r = r
        self.add = add

        if element_constructor is None:
            # Special case for words since their parent has a __call__ but
            #   not an _element_constructor_
            if isinstance(w1, Word_class):
                element_constructor = w1.parent()
            else:
                try:
                    element_constructor = w1.parent()._element_constructor_
                except AttributeError:
                    pass
        ShuffleProduct_abstract.__init__(self, w1, w2, element_constructor)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping_r
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: ShuffleProduct_overlapping_r(w,u,1).__repr__()
            'Overlapping shuffle product of word: 29 and word: 91 with 1 overlaps'
        """
        return "Overlapping shuffle product of %s and %s with %s overlaps" % (repr(self._l1), repr(self._l2), self.r)

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping_r
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: A = ShuffleProduct_overlapping_r(w,u,1)
            sage: B = ShuffleProduct_overlapping_r(u,w,2)
            sage: A == A
            True
            sage: A == B
            False
        """
        return ShuffleProduct_abstract.__eq__(self, other) and self.r == other.r

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping_r
            sage: w, u = Word([1,2]), Word([3,4])
            sage: ShuffleProduct_overlapping_r(w,u,1).list()
            [word: 424, word: 154, word: 442, word: 136, word: 352, word: 316]
            sage: w, u = map(Words(range(1,7)), [[1,2], [3,4]])
            sage: W = Words(range(1,7))
            sage: w, u = W([1,2]), W([3,4])
            sage: ShuffleProduct_overlapping_r(w, u, 1).list() #indirect doctest
            [word: 424, word: 154, word: 442, word: 136, word: 352, word: 316]

            sage: I, J = Composition([2, 2]), Composition([1, 1])
            sage: S = ShuffleProduct_overlapping_r(I, J, 1)
            sage: S.list()
            [[3, 2, 1], [2, 3, 1], [3, 1, 2], [2, 1, 3], [1, 3, 2], [1, 2, 3]]

        TESTS:

        We need to be explicitly more generic about the resulting parent
        when shuffling two compositions `I` and `J` (:trac:`15131`)::

            sage: I, J = Compositions(4)([2, 2]), Composition([1, 1])
            sage: S = ShuffleProduct_overlapping_r(I, J, 1, Compositions())
            sage: S.list()
            [[3, 2, 1], [2, 3, 1], [3, 1, 2], [2, 1, 3], [1, 3, 2], [1, 2, 3]]
        """
        EC = self._element_constructor_

        m = len(self._l1)
        n = len(self._l2)
        r = self.r
        add = self.add

        wc1, wc2 = self._l1, self._l2

        blank = [None] * (m + n - r)
        for iv in IntegerVectors(m, m + n - r, max_part=1):
            w = blank[:]
            filled_places = []
            unfilled_places = []
            # Fill in w1 into the iv slots
            i = 0
            for j in range(len(iv)):
                if iv[j] == 1:
                    w[j] = wc1[i]
                    i += 1
                    filled_places.append(j)
                else:
                    unfilled_places.append(j)

            # Choose r of these filled places
            for subset in itertools.combinations(filled_places, r):
                places_to_fill = sorted(unfilled_places + list(subset))

                # Fill in w2 into the places
                i = 0
                res = w[:]
                for j in places_to_fill:
                    if res[j] is not None:
                        res[j] = add(res[j], wc2[i])
                    else:
                        res[j] = wc2[i]
                    i += 1

                yield EC(res)


class ShuffleProduct_overlapping(ShuffleProduct_abstract):
    r"""
    The overlapping shuffle product of the two words ``w1`` and
    ``w2``.

    If `u` and `v` are two words whose letters belong to an
    additive monoid or to another kind of alphabet on which addition
    is well-defined, then the *overlapping shuffle product* of
    `u` and `v` is a certain multiset of words defined as follows:
    Let `a` and `b` be the lengths of `u` and `v`, respectively.
    Let `A` be the set `\{(0, 1), (0, 2), \cdots, (0, a)\}`, and
    let `B` be the set `\{(1, 1), (1, 2), \cdots, (1, b)\}`.
    Notice that the sets `A` and `B` are disjoint. We can make
    `A` and `B` into posets by setting `(k, i) \leq (k, j)` for
    all `k \in \{0, 1\}` and `i \leq j`. Then, `A \cup B` becomes
    a poset by disjoint union (we don't set `(0, i) \leq (1, i)`).
    Let `p` be the map from `A \cup B` to the set of all letters
    which sends every `(0, i)` to the `i`-th letter of `u`, and
    every `(1, j)` to the `j`-th letter of `v`. For every
    nonnegative integer `c` and every surjective map
    `f : A \cup B \to \{ 1, 2, \cdots, c \}` for which both
    restrictions `f \mid_A` and `f \mid_B` are strictly increasing,
    let `w(f)` be the length-`c` word such that for every
    `1 \leq k \leq c`, the `k`-th letter of `w(f)` equals
    `\sum_{j \in f^{-1}(k)} p(j)` (this sum always has either
    one or two addends). The overlapping shuffle product of `u`
    and `v` is then the multiset of all `w(f)` with `c` ranging
    over all nonnegative integers and `f` ranging
    over the surjective maps
    `f : A \cup B \to \{ 1, 2, \cdots, c \}` for which both
    restrictions `f \mid_A` and `f \mid_B` are strictly increasing.

    If one restricts `c` to a particular fixed nonnegative
    integer, then the multiset is instead called the *overlapping
    shuffle product with precisely `a + b - c` overlaps*. This is
    nonempty only if `\max \{a, b\} \leq c \leq a + b`.

    If `c = a + b`, then the overlapping shuffle product with
    precisely `a + b - c` overlaps is plainly the shuffle product
    (:class:`ShuffleProduct_w1w2`).

    INPUT:

    - ``w1``, ``w2`` -- iterables
    - ``element_constructor`` -- (default: the parent of ``w1``)
      the function used to construct the output
    - ``add`` -- (default: ``+``) the addition function

    EXAMPLES::

        sage: from sage.combinat.shuffle import ShuffleProduct_overlapping
        sage: w, u = [[2, 9], [9, 1]]
        sage: S = ShuffleProduct_overlapping(w, u)
        sage: sorted(S)
        [[2, 9, 1, 9],
         [2, 9, 9, 1],
         [2, 9, 9, 1],
         [2, 9, 10],
         [2, 18, 1],
         [9, 1, 2, 9],
         [9, 2, 1, 9],
         [9, 2, 9, 1],
         [9, 2, 10],
         [9, 3, 9],
         [11, 1, 9],
         [11, 9, 1],
         [11, 10]]
        sage: A = [{1,2}, {3,4}]
        sage: B = [{2,3}, {4,5,6}]
        sage: S = ShuffleProduct_overlapping(A, B, add=lambda X,Y: X.union(Y))
        sage: list(S)
        [[{1, 2}, {3, 4}, {2, 3}, {4, 5, 6}],
         [{1, 2}, {2, 3}, {3, 4}, {4, 5, 6}],
         [{1, 2}, {2, 3}, {4, 5, 6}, {3, 4}],
         [{2, 3}, {1, 2}, {3, 4}, {4, 5, 6}],
         [{2, 3}, {1, 2}, {4, 5, 6}, {3, 4}],
         [{2, 3}, {4, 5, 6}, {1, 2}, {3, 4}],
         [{1, 2, 3}, {3, 4}, {4, 5, 6}],
         [{1, 2}, {2, 3, 4}, {4, 5, 6}],
         [{1, 2, 3}, {4, 5, 6}, {3, 4}],
         [{1, 2}, {2, 3}, {3, 4, 5, 6}],
         [{2, 3}, {1, 2, 4, 5, 6}, {3, 4}],
         [{2, 3}, {1, 2}, {3, 4, 5, 6}],
         [{1, 2, 3}, {3, 4, 5, 6}]]
    """
    def __init__(self, w1, w2, element_constructor=None, add=operator.add):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: S = ShuffleProduct_overlapping(w,u)
            sage: TestSuite(S).run(skip="_test_an_element")
        """
        self._add = add

        if element_constructor is None:
            # Special case for words since their parent has a __call__ but
            #   not an _element_constructor_
            if isinstance(w1, Word_class):
                element_constructor = w1.parent()
            else:
                try:
                    element_constructor = self._l1.parent()._element_constructor_
                except AttributeError:
                    pass
        ShuffleProduct_abstract.__init__(self, w1, w2, element_constructor)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping
            sage: w, u = map(Words(range(20)), [[2, 9], [9, 1]])
            sage: ShuffleProduct_overlapping(w,u).__repr__()
            'Overlapping shuffle product of word: 29 and word: 91'
        """
        return "Overlapping shuffle product of %s and %s" % (repr(self._l1),
                                                             repr(self._l2))

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.shuffle import ShuffleProduct_overlapping
            sage: w, u = map(Words(range(10)), [[0,1],[2,3]])
            sage: S = ShuffleProduct_overlapping(w,u)
            sage: S.list()
            [word: 0123, word: 0213, word: 0231, word: 2013, word: 2031,
             word: 2301, word: 213, word: 033, word: 231, word: 024,
             word: 231, word: 204, word: 24]

            sage: w, u = map(Words(range(1,10)), [[1,2],[3,4]])
            sage: S = ShuffleProduct_overlapping(w,u)
            sage: S.list()
            [word: 1234, word: 1324, word: 1342, word: 3124, word: 3142,
             word: 3412, word: 424, word: 154, word: 442, word: 136,
             word: 352, word: 316, word: 46]
        """
        m = len(self._l1)
        n = len(self._l2)
        for r in range(min(m, n) + 1):
            yield from ShuffleProduct_overlapping_r(self._l1, self._l2, r,
                                                    self._element_constructor_,
                                                    add=self._add)
