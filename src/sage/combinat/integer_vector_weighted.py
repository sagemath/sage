"""
Weighted Integer Vectors

AUTHORS:

- Mike Hansen (2007): initial version, ported from MuPAD-Combinat
- Nicolas M. Thiery (2010-10-30): WeightedIntegerVectors(weights) + cleanup
"""
# ****************************************************************************
#  Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.misc.lazy_import import lazy_import
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.combinat.integer_vector import IntegerVector
from sage.combinat.permutation import Permutation

lazy_import('sage.combinat.words.word', 'Word')


class WeightedIntegerVectors(Parent, UniqueRepresentation):
    r"""
    The class of integer vectors of `n` weighted by ``weight``, that is, the
    nonnegative integer vectors `(v_1, \ldots, v_{\ell})`
    satisfying `\sum_{i=1}^{\ell} v_i w_i = n` where `\ell` is
    ``length(weight)`` and `w_i` is ``weight[i]``.

    INPUT:

    - ``n`` -- nonnegative integer (optional)

    - ``weight`` -- tuple (or list or iterable) of positive integers

    EXAMPLES::

        sage: WeightedIntegerVectors(8, [1,1,2])
        Integer vectors of 8 weighted by [1, 1, 2]
        sage: WeightedIntegerVectors(8, [1,1,2]).first()
        [0, 0, 4]
        sage: WeightedIntegerVectors(8, [1,1,2]).last()
        [8, 0, 0]
        sage: WeightedIntegerVectors(8, [1,1,2]).cardinality()
        25
        sage: w = WeightedIntegerVectors(8, [1,1,2]).random_element()
        sage: w.parent() is WeightedIntegerVectors(8, [1,1,2])
        True

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
        ValueError: the weights must be specified

    .. TODO::

        Should the order of the arguments ``n`` and ``weight`` be
        exchanged to simplify the logic?
    """
    @staticmethod
    def __classcall_private__(cls, n=None, weight=None):
        """
        Normalize inputs to ensure a unique representation.

        TESTS::

            sage: W = WeightedIntegerVectors(8, [1,1,2])
            sage: W2 = WeightedIntegerVectors(int(8), (1,1,2))
            sage: W is W2
            True
        """
        if weight is None:
            if n is None:
                raise ValueError("the weights must be specified")
            if n in ZZ:
                weight = (n,)
            else:
                weight = tuple(n)
            n = None

        weight = tuple(weight)
        if n is None:
            return WeightedIntegerVectors_all(weight)

        return super().__classcall__(cls, n, weight)

    def __init__(self, n, weight):
        """
        TESTS::

            sage: WIV = WeightedIntegerVectors(8, [1,1,2])
            sage: TestSuite(WIV).run()
        """
        self._n = n
        self._weights = weight
        Parent.__init__(self, category=FiniteEnumeratedSets())

    Element = IntegerVector

    def _element_constructor_(self, lst):
        """
        Construct an element of ``self`` from ``lst``.

        EXAMPLES::

            sage: WIV = WeightedIntegerVectors(3, [2,1,1])
            sage: elt = WIV([1, 1, 0]); elt
            [1, 1, 0]
            sage: elt.parent() is WIV
            True
            sage: WIV([1, 1, 0])
            [1, 1, 0]
            sage: WIV([1, 2, 0])
            Traceback (most recent call last):
            ...
            ValueError: cannot convert [1, 2, 0] into Integer vectors of 3
             weighted by [2, 1, 1]
        """
        if isinstance(lst, IntegerVector):
            if lst.parent() is self:
                return lst
        if lst not in self:
            raise ValueError("cannot convert %s into %s" % (lst, self))
        return self.element_class(self, lst)

    def _repr_(self):
        """
        TESTS::

            sage: WeightedIntegerVectors(8, [1,1,2])
            Integer vectors of 8 weighted by [1, 1, 2]
        """
        return "Integer vectors of %s weighted by %s" % (self._n, list(self._weights))

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
        if not isinstance(x, (list, IntegerVector, Permutation)):
            return False
        if len(self._weights) != len(x):
            return False
        s = 0
        for i, xi in enumerate(x):
            if not isinstance(xi, (int, Integer)) and xi not in ZZ:
                return False
            s += xi * self._weights[i]
        return s == self._n

    def _recfun(self, n, l):
        """
        EXAMPLES::

            sage: w = WeightedIntegerVectors(3, [2,1,1])
            sage: list(w._recfun(3, [1,1,2]))
            [[0, 1, 1], [1, 0, 1], [0, 3, 0], [1, 2, 0], [2, 1, 0], [3, 0, 0]]
        """
        w = l[-1]
        l = l[:-1]
        if l == []:
            d = int(n) // int(w)
            if n % w == 0:
                yield [d]
                # Otherwise: bad branch
            return

        for d in range(int(n) // int(w), -1, -1):
            for x in self._recfun(n - d * w, l):
                yield x + [d]

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
            sage: all(sorted(map(list, iv[k])) == sorted(map(list, ivw[k])) for k in range(11))
            True

        ::

            sage: ivw = [ WeightedIntegerVectors(k, [2,3,7]) for k in range(11) ]
            sage: all(i.cardinality() == len(i.list()) for i in ivw)
            True
        """
        if not self._weights:
            if self._n == 0:
                yield self.element_class(self, [])
            return

        perm = Word(self._weights).standard_permutation()
        perm = [len(self._weights) - i for i in perm]
        l = sorted(self._weights, reverse=True)
        for x in iterator_fast(self._n, l):
            yield self.element_class(self, [x[i] for i in perm])
            # .action(x)
            # _left_to_right_multiply_on_right(Permutation(x))


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

    def __init__(self, weight):
        """
        TESTS::

            sage: C = WeightedIntegerVectors([2,1,3])
            sage: C.category()
            Category of facade infinite enumerated sets with grading
            sage: TestSuite(C).run()
        """
        self._weights = weight
        from sage.sets.family import Family
        from sage.sets.non_negative_integers import NonNegativeIntegers
        # Use "partial" to make the basis function (with the weights
        # argument specified) pickleable.  Otherwise, it seems to
        # cause problems...
        from functools import partial
        F = Family(NonNegativeIntegers(), partial(WeightedIntegerVectors, weight=weight))
        cat = (SetsWithGrading(), InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self, F, facade=True, keepkey=False,
                                             category=cat)

    def _repr_(self):
        """
        EXAMPLES::

            sage: WeightedIntegerVectors([2,1,3])
            Integer vectors weighted by [2, 1, 3]
        """
        return "Integer vectors weighted by %s" % list(self._weights)

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
        return (isinstance(x, (list, IntegerVector, Permutation))
                and len(x) == len(self._weights)
                and all(i in ZZ and i >= 0 for i in x))

    def subset(self, size=None):
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
        return sum(exp * deg for exp, deg in zip(x, self._weights))


def iterator_fast(n, l):
    """
    Iterate over all ``l`` weighted integer vectors with total weight ``n``.

    INPUT:

    - ``n`` -- integer
    - ``l`` -- the weights in weakly decreasing order

    EXAMPLES::

        sage: from sage.combinat.integer_vector_weighted import iterator_fast
        sage: list(iterator_fast(3, [2,1,1]))
        [[1, 1, 0], [1, 0, 1], [0, 3, 0], [0, 2, 1], [0, 1, 2], [0, 0, 3]]
        sage: list(iterator_fast(2, [2]))
        [[1]]

    Test that :issue:`20491` is fixed::

        sage: type(list(iterator_fast(2, [2]))[0][0])
        <class 'sage.rings.integer.Integer'>
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


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.integer_vector_weighted', 'WeightedIntegerVectors_nweight', WeightedIntegerVectors)
