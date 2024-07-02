r"""
Shuffle product of words

.. SEEALSO::

    The module :mod:`sage.combinat.shuffle` contains a more general
    implementation of shuffle product.
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.combinat.words.word import Word_class, Word
from sage.arith.misc import binomial
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.composition import Composition
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class ShuffleProduct_w1w2(Parent, UniqueRepresentation):
    def __init__(self, w1, w2, check=True):
        r"""
        The shuffle product of the two words ``w1`` and ``w2``.

        If `u` and `v` are two words, then the *shuffle product* of
        `u` and `v` is a certain multiset of words defined as follows:
        Let `a` and `b` be the lengths of `u` and `v`, respectively.
        For every `a`-element subset `I` of `\{1, 2, \cdots, a+b\}`,
        let `w(I)` be the length-`a+b` word such that:

        - for every `1 \leq k \leq a`, the `i_k`-th letter of `w(I)`
          is the `k`-th letter of `u`, where `i_k` is the
          `k`-th smallest element of `I`;

        - for every `1 \leq l \leq b`, the `j_l`-th letter of `w(I)`
          is the `l`-th letter of `v`, where `j_l` is the
          `l`-th smallest element of
          `\{1, 2, \cdots, a+b\} \setminus I`.

        The shuffle product of `u` and `v` is then the multiset of
        all `w(I)` with `I` ranging over the `a`-element subsets of
        `\{1, 2, \cdots, a+b\}`.

        INPUT:

        - ``check`` -- boolean (default: ``True``); whether to check that
          all words in the shuffle product belong to the correct parent

        EXAMPLES::

            sage: from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            sage: W = Words([1,2,3,4])
            sage: s = ShuffleProduct_w1w2(W([1,2]),W([3,4]))
            sage: sorted(s)
            [word: 1234, word: 1324, word: 1342, word: 3124,
             word: 3142, word: 3412]
            sage: s == loads(dumps(s))
            True
            sage: TestSuite(s).run()

            sage: s = ShuffleProduct_w1w2(W([1,4,3]),W([2]))
            sage: sorted(s)
            [word: 1243, word: 1423, word: 1432, word: 2143]

            sage: s = ShuffleProduct_w1w2(W([1,4,3]),W([]))
            sage: sorted(s)
            [word: 143]

        TESTS::

            sage: W = Words([1,2,3,4])
            sage: s = ShuffleProduct_w1w2(W([1,2]), W([3,4]), check=False)
            sage: len(list(s))
            6
        """
        self._w1 = w1
        self._w2 = w2
        self._check = bool(check)
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            sage: W = Words("abcd")
            sage: repr(ShuffleProduct_w1w2(W("ab"),W("cd")))
            'Shuffle product of word: ab and word: cd'
        """
        return "Shuffle product of %s and %s" % (repr(self._w1), repr(self._w2))

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            sage: W = Words("abcd")
            sage: w = W("ab")
            sage: u = W("cd")
            sage: S = ShuffleProduct_w1w2(w,u)
            sage: w*u in S
            True
            sage: all(w.is_subword_of(x) for x in S)
            True
            sage: w in S
            False

        We check that :issue:`14121` is solved::

            sage: w = W('ab')
            sage: x = W('ac')
            sage: x*w in w.shuffle(x)
            True
        """
        if not isinstance(x, Word_class):
            return False
        if x.length() != self._w1.length() + self._w2.length():
            return False
        w1 = list(self._w1)
        w2 = list(self._w2)
        wx = list(x)
        for _ in range(len(wx)):
            try:
                letter = wx.pop(0)
            except IndexError:
                return False
            if w1 and w2 and letter == w1[0] == w2[0]:
                return (Word(wx) in self._w1[1:].shuffle(self._w2) or
                        Word(wx) in self._w1.shuffle(self._w2[1:]))
            if w1 and letter == w1[0]:
                w1.pop(0)
            elif w2 and letter == w2[0]:
                w2.pop(0)
            else:
                return False
        return not wx

    def cardinality(self):
        r"""
        Return the number of words in the shuffle product
        of ``w1`` and ``w2``.

        This is understood as a multiset cardinality, not as a
        set cardinality; it does not count the distinct words only.

        It is given by `\binom{l_1+l_2}{l_1}`, where `l_1` is the
        length of ``w1`` and where `l_2` is the length of ``w2``.

        EXAMPLES::

            sage: from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            sage: w, u = map(Words("abcd"), ["ab", "cd"])
            sage: S = ShuffleProduct_w1w2(w,u)
            sage: S.cardinality()
            6

            sage: w, u = map(Words("ab"), ["ab", "ab"])
            sage: S = ShuffleProduct_w1w2(w,u)
            sage: S.cardinality()
            6
        """
        len_w1 = self._w1.length()
        len_w2 = self._w2.length()
        return binomial(len_w1 + len_w2, len_w1)

    def __iter__(self):
        """
        Return an iterator for the words in the
        shuffle product of ``w1`` and ``w2``.

        EXAMPLES::

            sage: from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            sage: w, u = map(Words("abcd"), ["ab", "cd"])
            sage: S = ShuffleProduct_w1w2(w,u)
            sage: S.list() #indirect test
            [word: abcd, word: acbd, word: acdb, word: cabd,
             word: cadb, word: cdab]

            sage: I = Composition([1, 1])
            sage: J = Composition([2])
            sage: S = ShuffleProduct_w1w2(I, J)
            sage: next(iter(S))
            [1, 1, 2]

        TESTS:

        Sage is no longer confused by a too-restrictive parent of `I`
        when shuffling compositions `I` and `J` (cf. :issue:`15131`)::

            sage: I = Compositions(2)([1, 1])
            sage: J = Composition([2])
            sage: S = ShuffleProduct_w1w2(I, J)
            sage: S.list()
            [[1, 1, 2], [1, 2, 1], [2, 1, 1]]
        """
        n1 = len(self._w1)
        n2 = len(self._w2)
        w1_parent = self._w1.parent()
        use_w1_parent = True
        try:
            w1_parent(list(self._w1) + list(self._w2), check=self._check)
        except (ValueError, TypeError):
            use_w1_parent = False
            if isinstance(self._w1, Composition):
                large_parent = Composition
            elif isinstance(self._w1, Word_class):
                large_parent = Word
        for iv in IntegerVectors(n1, n1 + n2, max_part=1):
            it1 = iter(self._w1)
            it2 = iter(self._w2)
            w = [next(it1) if v else next(it2) for v in iv]
            if use_w1_parent:
                yield w1_parent(w, check=self._check)
            else:
                yield large_parent(w)


class ShuffleProduct_shifted(ShuffleProduct_w1w2):
    def __init__(self, w1, w2, check=True):
        """
        Shifted shuffle product of ``w1`` with ``w2``.

        This is the shuffle product of ``w1`` with the word
        obtained by adding the length of ``w1`` to every letter
        of ``w2``.

        Note that this class is meant to be used for words; it
        misbehaves when ``w1`` is a permutation or composition.

        INPUT:

        - ``check`` -- boolean (default: ``True``); whether to check that
          all words in the shuffle product belong to the correct parent

        EXAMPLES::

            sage: from sage.combinat.words.shuffle_product import ShuffleProduct_shifted
            sage: w, u = Word([1,2]), Word([3,4])
            sage: S = ShuffleProduct_shifted(w,u)
            sage: S == loads(dumps(S))
            True

        TESTS::

            sage: w, u = Word([1,2]), Word([3,4])
            sage: S = ShuffleProduct_shifted(w, u, check=False)
            sage: len(list(S))
            6
        """
        shift = w1.length()
        shifted_w2 = w1.parent()([x + shift for x in w2], check=check)
        ShuffleProduct_w1w2.__init__(self, w1, shifted_w2, check)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.shuffle_product import ShuffleProduct_shifted
            sage: w, u = Word([0,1]), Word([2,3])
            sage: ShuffleProduct_shifted(w,u).__repr__()
            'Shuffle product of word: 01 and word: 45'
        """
        return "Shuffle product of %s and %s" % (repr(self._w1), repr(self._w2))
