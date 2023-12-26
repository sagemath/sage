# sage_setup: distribution = sagemath-categories
from sage.sets.family import Family
from sage.combinat.integer_lists import IntegerListsLex
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets


def IntegerListsNN(**kwds):
    """
    Lists of nonnegative integers with constraints.

    This function returns the union of ``IntegerListsLex(n, **kwds)``
    where `n` ranges over all nonnegative integers.

    EXAMPLES::

        sage: from sage.combinat.integer_lists.nn import IntegerListsNN
        sage: L = IntegerListsNN(max_length=3, max_slope=-1)
        sage: L
        Disjoint union of Lazy family (<lambda>(i))_{i in Non negative integer semiring}
        sage: it = iter(L)
        sage: for _ in range(20):
        ....:     print(next(it))
        []
        [1]
        [2]
        [3]
        [2, 1]
        [4]
        [3, 1]
        [5]
        [4, 1]
        [3, 2]
        [6]
        [5, 1]
        [4, 2]
        [3, 2, 1]
        [7]
        [6, 1]
        [5, 2]
        [4, 3]
        [4, 2, 1]
        [8]
    """
    return DisjointUnionEnumeratedSets(Family(NN, lambda i: IntegerListsLex(i, **kwds)))
