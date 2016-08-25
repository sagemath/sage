"""
Weighted Integer Vectors

AUTHORS:

 - Mike Hansen (2007): initial version, ported from MuPAD-Combinat
 - Nicolas M. Thiery (2010-10-30): WeightedIntegerVectors(weights) + cleanup

.. WARNING::

    The list(self) function in this file used the :class:`Permutation` class improperly, returning
    a list of, generally speaking, invalid permutations (repeated entries, including 0).
"""
#*****************************************************************************
#  Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from six.moves.builtins import list as builtinlist
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.combinat.words.word import Word
from .permutation import Permutation

def WeightedIntegerVectors(n = None, weight = None):
    """
    Returns the combinatorial class of integer vectors of ``n``
    weighted by ``weight``, that is, the nonnegative integer vectors
    `(v_1,\\dots,v_{length(weight)})` satisfying `\\sum_i v_i
    weight[i]==n`.

    INPUT:

     - ``n`` -- a non negative integer (optional)

     - ``weight`` -- a tuple (or list or iterable) of positive integers

    EXAMPLES::

        sage: WeightedIntegerVectors(8, [1,1,2])
        Integer vectors of 8 weighted by [1, 1, 2]
        sage: WeightedIntegerVectors(8, [1,1,2]).first()
        [0, 0, 4]
        sage: WeightedIntegerVectors(8, [1,1,2]).last()
        [8, 0, 0]
        sage: WeightedIntegerVectors(8, [1,1,2]).cardinality()
        25
        sage: WeightedIntegerVectors(8, [1,1,2]).random_element()
        [1, 1, 3]

        sage: WeightedIntegerVectors([1,1,2])
        Integer vectors weighted by [1, 1, 2]
        sage: WeightedIntegerVectors([1,1,2]).cardinality()
        +Infinity
        sage: WeightedIntegerVectors([1,1,2]).first()
        [0, 0, 0]

    TESTS::

        sage: WeightedIntegerVectors(None,None)
        Traceback (most recent call last):
        ...
        ValueError: weights should be specified

    .. TODO::

        should the order of the arguments ``n`` and ``weight`` be
        exchanged to simplify the logic ?
    """
    if weight is None and n is not None:
        weight = n
        n = None
    if weight is None:
        raise ValueError("weights should be specified")
    weight = tuple(weight)
    if n is None:
        return WeightedIntegerVectors_all(weight)
    else:
        return WeightedIntegerVectors_nweight(n, weight)

class WeightedIntegerVectors_all(DisjointUnionEnumeratedSets):
    r"""
    Set of weighted integer vectors.

    EXAMPLES::

        sage: W = WeightedIntegerVectors([3,1,1,2,1]); W
        Integer vectors weighted by [3, 1, 1, 2, 1]
        sage: W.cardinality()
        +Infinity

        sage: W12 = W.graded_component(12)
        sage: W12.an_element()
        [4, 0, 0, 0, 0]
        sage: W12.last()
        [0, 12, 0, 0, 0]
        sage: W12.cardinality()
        441
        sage: for w in W12: print(w)
        [4, 0, 0, 0, 0]
        [3, 0, 0, 1, 1]
        [3, 0, 1, 1, 0]
        ...
        [0, 11, 1, 0, 0]
        [0, 12, 0, 0, 0]
    """
    def __init__(self, weights):
        """
        TESTS::

            sage: C = WeightedIntegerVectors([2,1,3])
            sage: C.__class__
            <class 'sage.combinat.integer_vector_weighted.WeightedIntegerVectors_all_with_category'>
            sage: C.category()
            Join of Category of sets with grading and Category of infinite enumerated sets
            sage: TestSuite(C).run()
        """
        self._weights = weights
        from sage.sets.all import Family, NonNegativeIntegers
        # Use "partial" to make the basis function (with the weights
        # argument specified) pickleable.  Otherwise, it seems to
        # cause problems...
        from functools import partial
        F = Family(NonNegativeIntegers(), partial(WeightedIntegerVectors, weight = weights))
        DisjointUnionEnumeratedSets.__init__(self, F, facade=True, keepkey=False,
                                             category = (SetsWithGrading(), InfiniteEnumeratedSets()))

    def _repr_(self):
        """
        EXAMPLES::

            sage: WeightedIntegerVectors([2,1,3])
            Integer vectors weighted by [2, 1, 3]
        """
        return "Integer vectors weighted by %s"%list(self._weights)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [] in WeightedIntegerVectors([])
            True
            sage: [3,0,0] in WeightedIntegerVectors([2,1,1])
            True
            sage: [3,0] in WeightedIntegerVectors([2,1,1])
            False
            sage: [3,-1,0] in WeightedIntegerVectors([2,1,1])
            False
        """
        return (isinstance(x, (builtinlist, Permutation))
                and len(x) == len(self._weights)
                and all(isinstance(i, (int, Integer)) and i >= 0 for i in x))

    def subset(self, size = None):
        """
        EXAMPLES::

            sage: C = WeightedIntegerVectors([2,1,3])
            sage: C.subset(4)
            Integer vectors of 4 weighted by [2, 1, 3]
        """
        if size is None:
            return self
        return self._family[size]

    def grading(self, x): # or degree / grading
        """
        EXAMPLES::

            sage: C = WeightedIntegerVectors([2,1,3])
            sage: C.grading((2,1,1))
            8
        """
        return sum([exp*deg for exp,deg in zip(x, self._weights)])

class WeightedIntegerVectors_nweight(UniqueRepresentation, Parent):
    def __init__(self, n, weight):
        """
        TESTS::

            sage: C = WeightedIntegerVectors(8, [1,1,2])
            sage: C.__class__
            <class 'sage.combinat.integer_vector_weighted.WeightedIntegerVectors_nweight_with_category'>
            sage: TestSuite(C).run()
        """
        Parent.__init__(self, category = FiniteEnumeratedSets())
        self._n = n
        self._weights = weight

    def _repr_(self):
        """
        TESTS::

            sage: repr(WeightedIntegerVectors(8, [1,1,2]))
            'Integer vectors of 8 weighted by [1, 1, 2]'
        """
        return "Integer vectors of %s weighted by %s"%(self._n, list(self._weights))

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [] in WeightedIntegerVectors(0, [])
            True
            sage: [] in WeightedIntegerVectors(1, [])
            False
            sage: [3,0,0] in WeightedIntegerVectors(6, [2,1,1])
            True
            sage: [1] in WeightedIntegerVectors(1, [1])
            True
            sage: [1] in WeightedIntegerVectors(2, [2])
            True
            sage: [2] in WeightedIntegerVectors(4, [2])
            True
            sage: [2, 0] in WeightedIntegerVectors(4, [2, 2])
            True
            sage: [2, 1] in WeightedIntegerVectors(4, [2, 2])
            False
            sage: [2, 1] in WeightedIntegerVectors(6, [2, 2])
            True
            sage: [2, 1, 0] in WeightedIntegerVectors(6, [2, 2])
            False
            sage: [0] in WeightedIntegerVectors(0, [])
            False
        """
        if not isinstance(x, (builtinlist, Permutation)):
            return False
        if len(self._weights) != len(x):
            return False
        s = 0
        for i in range(len(x)):
            if not isinstance(x[i], (int, Integer)):
                return False
            s += x[i]*self._weights[i]
        if s != self._n:
            return False

        return True

    def __iter__(self):
        """
        TESTS::

            sage: WeightedIntegerVectors(7, [2,2]).list()
            []
            sage: WeightedIntegerVectors(3, [2,1,1]).list()
            [[1, 0, 1], [1, 1, 0], [0, 0, 3], [0, 1, 2], [0, 2, 1], [0, 3, 0]]

        ::

            sage: ivw = [ WeightedIntegerVectors(k, [1,1,1]) for k in range(11) ]
            sage: iv  = [ IntegerVectors(k, 3) for k in range(11) ]
            sage: all( [ sorted(iv[k].list()) == sorted(ivw[k].list()) for k in range(11) ] )
            True

        ::

            sage: ivw = [ WeightedIntegerVectors(k, [2,3,7]) for k in range(11) ]
            sage: all( [ i.cardinality() == len(i.list()) for i in ivw] )
            True
        """
        if not self._weights:
            if self._n == 0:
                yield []
            return

        perm = Word(self._weights).standard_permutation()
        perm = [len(self._weights)-i for i in perm]
        l = [x for x in sorted(self._weights, reverse=True)]
        for x in iterator_fast(self._n, l):
            yield [x[i] for i in perm]
            #.action(x)
            #_left_to_right_multiply_on_right(Permutation(x))

def iterator_fast(n, l):
    """
    Iterate over all ``l`` weighted integer vectors with total weight ``n``.

    INPUT:

    - ``n`` -- an integer
    - ``l`` -- the weights in weakly decreasing order

    EXAMPLES::

        sage: from sage.combinat.integer_vector_weighted import iterator_fast
        sage: list(iterator_fast(3, [2,1,1]))
        [[1, 1, 0], [1, 0, 1], [0, 3, 0], [0, 2, 1], [0, 1, 2], [0, 0, 3]]
        sage: list(iterator_fast(2, [2]))
        [[1]]

    Test that :trac:`20491` is fixed::

        sage: type(list(iterator_fast(2, [2]))[0][0])
        <type 'sage.rings.integer.Integer'>
    """
    if n < 0:
        return

    zero = ZZ.zero()
    one = ZZ.one()

    if not l:
        if n == 0:
            yield []
        return
    if len(l) == 1:
        if n % l[0] == 0:
            yield [n // l[0]]
        return

    k = 0
    cur = [n // l[k] + one]
    rem = n - cur[-1] * l[k] # Amount remaining
    while cur:
        cur[-1] -= one
        rem += l[k]
        if rem == zero:
            yield cur + [zero] * (len(l) - len(cur))
        elif cur[-1] < zero or rem < zero:
            rem += cur.pop() * l[k]
            k -= 1
        elif len(l) == len(cur) + 1:
            if rem % l[-1] == zero:
                yield cur + [rem // l[-1]]
        else:
            k += 1
            cur.append(rem // l[k] + one)
            rem -= cur[-1] * l[k]

