r"""
Tuples
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.libs.gap.libgap import libgap
from sage.rings.integer_ring import ZZ
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets


class Tuples(Parent, UniqueRepresentation):
    """
    Return the enumerated set of ordered tuples of S of length k.

    An ordered tuple of length k of set is an ordered selection with
    repetition and is represented by a list of length k containing
    elements of set.

    EXAMPLES::

        sage: S = [1,2]
        sage: Tuples(S,3).list()
        [[1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1], [1, 1, 2],
         [2, 1, 2], [1, 2, 2], [2, 2, 2]]
        sage: mset = ["s","t","e","i","n"]
        sage: Tuples(mset,2).list()
        [['s', 's'], ['t', 's'], ['e', 's'], ['i', 's'], ['n', 's'],
         ['s', 't'], ['t', 't'], ['e', 't'], ['i', 't'], ['n', 't'],
         ['s', 'e'], ['t', 'e'], ['e', 'e'], ['i', 'e'], ['n', 'e'],
         ['s', 'i'], ['t', 'i'], ['e', 'i'], ['i', 'i'], ['n', 'i'],
         ['s', 'n'], ['t', 'n'], ['e', 'n'], ['i', 'n'], ['n', 'n']]

    ::

        sage: K.<a> = GF(4, 'a')
        sage: mset = [x for x in K if x != 0]
        sage: Tuples(mset,2).list()
        [[a, a], [a + 1, a], [1, a], [a, a + 1], [a + 1, a + 1], [1, a + 1],
         [a, 1], [a + 1, 1], [1, 1]]
    """
    @staticmethod
    def __classcall_private__(cls, S, k):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: T = Tuples(['l','i','t'],2); T
            Tuples of ('l', 'i', 't') of length 2
        """
        return super(Tuples, cls).__classcall__(cls, tuple(S), k)

    def __init__(self, S, k):
        """
        TESTS::

            sage: T = Tuples([1,2,3],2)
            sage: T == loads(dumps(T))
            True
        """
        self.S = S
        self.k = k
        self._index_list = [S.index(s) for s in S]
        category = FiniteEnumeratedSets()
        Parent.__init__(self, category=category)

    def __repr__(self):
        """
        TESTS::

            sage: repr(Tuples([1,2,3],2))
            'Tuples of (1, 2, 3) of length 2'
        """
        return "Tuples of %s of length %s" % (self.S, self.k)

    def __iter__(self):
        """
        EXAMPLES::

            sage: S = [1,2]
            sage: Tuples(S,3).list()
            [[1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1], [1, 1, 2],
             [2, 1, 2], [1, 2, 2], [2, 2, 2]]
            sage: mset = ["s","t","e","i","n"]
            sage: Tuples(mset,2).list()
            [['s', 's'], ['t', 's'], ['e', 's'], ['i', 's'], ['n', 's'],
             ['s', 't'], ['t', 't'], ['e', 't'], ['i', 't'],
             ['n', 't'], ['s', 'e'], ['t', 'e'], ['e', 'e'], ['i', 'e'],
             ['n', 'e'], ['s', 'i'], ['t', 'i'], ['e', 'i'],
             ['i', 'i'], ['n', 'i'], ['s', 'n'], ['t', 'n'], ['e', 'n'],
             ['i', 'n'], ['n', 'n']]
        """
        S = self.S
        k = self.k
        import copy
        if k <= 0:
            yield []
            return
        if k == 1:
            for x in S:
                yield [x]
            return

        for s in S:
            for x in Tuples(S, k - 1):
                y = copy.copy(x)
                y.append(s)
                yield y

    def cardinality(self):
        """
        EXAMPLES::

            sage: S = [1,2,3,4,5]
            sage: Tuples(S,2).cardinality()
            25
            sage: S = [1,1,2,3,4,5]
            sage: Tuples(S,2).cardinality()
            25
        """
        return ZZ(libgap.NrTuples(self._index_list, ZZ(self.k)))


Tuples_sk = Tuples


class UnorderedTuples(Parent, UniqueRepresentation):
    """
    Return the enumerated set of unordered tuples of S of length k.

    An unordered tuple of length k of set is a unordered selection with
    repetitions of set and is represented by a sorted list of length k
    containing elements from set.

    EXAMPLES::

        sage: S = [1,2]
        sage: UnorderedTuples(S,3).list()
        [[1, 1, 1], [1, 1, 2], [1, 2, 2], [2, 2, 2]]
        sage: UnorderedTuples(["a","b","c"],2).list()
        [['a', 'a'], ['a', 'b'], ['a', 'c'], ['b', 'b'], ['b', 'c'],
         ['c', 'c']]
    """
    @staticmethod
    def __classcall_private__(cls, S, k):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: T = UnorderedTuples(['l','i','t'],2); T
            Unordered tuples of ('l', 'i', 't') of length 2
        """
        return super(UnorderedTuples, cls).__classcall__(cls, tuple(S), k)

    def __init__(self, S, k):
        """
        TESTS::

            sage: T = Tuples([1,2,3],2)
            sage: T == loads(dumps(T))
            True
        """
        self.S = S
        self.k = k
        self._index_list = [S.index(s) for s in S]
        category = FiniteEnumeratedSets()
        Parent.__init__(self, category=category)

    def __repr__(self):
        """
        TESTS::

            sage: repr(UnorderedTuples([1,2,3],2))
            'Unordered tuples of (1, 2, 3) of length 2'
        """
        return "Unordered tuples of %s of length %s" % (self.S, self.k)

    def list(self):
        """
        EXAMPLES::

            sage: S = [1,2]
            sage: UnorderedTuples(S,3).list()
            [[1, 1, 1], [1, 1, 2], [1, 2, 2], [2, 2, 2]]
            sage: UnorderedTuples(["a","b","c"],2).list()
            [['a', 'a'], ['a', 'b'], ['a', 'c'], ['b', 'b'], ['b', 'c'],
             ['c', 'c']]
        """
        ans = libgap.UnorderedTuples(self._index_list, ZZ(self.k))
        return [[self.S[i] for i in l] for l in ans]

    def cardinality(self):
        """
        EXAMPLES::

            sage: S = [1,2,3,4,5]
            sage: UnorderedTuples(S,2).cardinality()
            15
        """
        return ZZ(libgap.NrUnorderedTuples(self._index_list, ZZ(self.k)))


UnorderedTuples_sk = UnorderedTuples
