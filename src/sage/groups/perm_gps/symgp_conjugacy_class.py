r"""
Conjugacy Classes Of The Symmetric Group

AUTHORS:

- Vincent Delecroix, Travis Scrimshaw (2014-11-23)
"""

from sage.groups.conjugacy_classes import ConjugacyClass, ConjugacyClassGAP
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.combinat.partition import Partitions_n, _Partitions
from sage.combinat.set_partition import SetPartitions
from sage.combinat.permutation import Permutation, from_cycles
from sage.sets.set import Set
import itertools


class SymmetricGroupConjugacyClassMixin(object):
    r"""
    Mixin class which contains methods for conjugacy classes of
    the symmetric group.
    """
    def __init__(self, domain, part):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,2), (3,4,5)])
            sage: C = G.conjugacy_class(Partition([3,2]))
            sage: type(C._part)
            <class 'sage.combinat.partition.Partitions_n_with_category.element_class'>
        """
        P = Partitions_n(len(domain))
        self._part = P(part)
        self._domain = domain
        self._set = None

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: G.conjugacy_class(Partition([4]))
            Conjugacy class of cycle type [4] in
             Symmetric group of order 4! as a permutation group
        """
        return "Conjugacy class of cycle type %s in %s"%(self._part, self._parent)

    def __eq__(self, other):
        r"""
        Comparison of conjugacy classes is done by comparing the
        defining cycle types.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,2), (3,4,5)])
            sage: C = G.conjugacy_class(Partition([3,2]))
            sage: Cp = G.conjugacy_class(g)
            sage: C == Cp
            True
        """
        if not isinstance(other, SymmetricGroupConjugacyClassMixin):
            return False
        return self._part == other._part

    def __ne__(self, other):
        """
        Test for unequality.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,3), (2,4,5)])
            sage: C = G.conjugacy_class(Partition([3,2]))
            sage: Cp = G.conjugacy_class(g)
            sage: C != Cp
            False
        """
        return not (self == other)

    def partition(self):
        """
        Return the partition of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,2), (3,4,5)])
            sage: C = G.conjugacy_class(g)
        """
        return self._part


class SymmetricGroupConjugacyClass(SymmetricGroupConjugacyClassMixin, ConjugacyClassGAP):
    """
    A conjugacy class of the symmetric group.

    INPUT:

    - ``group`` -- the symmetric group
    - ``part`` -- a partition or an element of ``group``
    """
    def __init__(self, group, part):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,2), (3,4,5)])
            sage: C = G.conjugacy_class(g)
            sage: TestSuite(C).run()
            sage: C = G.conjugacy_class(Partition([3,2]))
            sage: TestSuite(C).run()
        """
        if isinstance(part, PermutationGroupElement) and part.parent() is group:
            elt = part
            part = sorted([len(x) for x in part.cycle_tuples(True)], reverse=True)
        else:
            elt = default_representative(part, group)
        SymmetricGroupConjugacyClassMixin.__init__(self, group.domain(), part)
        ConjugacyClassGAP.__init__(self, group, elt)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: C = G.conjugacy_class(Partition([3,1]))
            sage: for x in C: x
            (2,3,4)
            (2,4,3)
            (1,2,3)
            (1,3,2)
            (1,2,4)
            (1,4,2)
            (1,3,4)
            (1,4,3)
        """
        if self._set:
            for x in self._set:
                yield x
            return

        for x in conjugacy_class_iterator(self._part, self._parent.domain()):
            yield self._parent.element_class(x, self._parent, check=False)

    def set(self):
        r"""
        The set of all elements in the conjugacy class ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: g = G((1,2))
            sage: C = G.conjugacy_class(g)
            sage: S = [(2,3), (1,2), (1,3)]
            sage: C.set() == Set(G(x) for x in S)
            True
        """
        if not self._set:
            self._set = Set(self._parent.element_class(x, self._parent, check=False)
                            for x in conjugacy_class_iterator(self._part, self._domain) )
        return self._set

class PermutationsConjugacyClass(SymmetricGroupConjugacyClassMixin, ConjugacyClass):
    """
    A conjugacy class of the permutations of `n`.

    INPUT:

    - ``P`` -- the permutations of `n`
    - ``part`` -- a partition or an element of ``P``
    """
    def __init__(self, P, part):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = Permutations(5)
            sage: g = G([2, 1, 4, 5, 3])
            sage: C = G.conjugacy_class(g)
            sage: TestSuite(C).run()
            sage: C = G.conjugacy_class(Partition([3,2]))
            sage: TestSuite(C).run()
        """
        if isinstance(part, Permutation) and part.parent() is P:
            elt = part
            part = elt.cycle_type()
        else:
            elt = P.element_in_conjugacy_classes(part)
        SymmetricGroupConjugacyClassMixin.__init__(self, range(1, P.n+1), part)
        ConjugacyClass.__init__(self, P, elt)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: G = Permutations(4)
            sage: C = G.conjugacy_class(Partition([3,1]))
            sage: for x in C: x
            [1, 3, 4, 2]
            [1, 4, 2, 3]
            [2, 3, 1, 4]
            [3, 1, 2, 4]
            [2, 4, 3, 1]
            [4, 1, 3, 2]
            [3, 2, 4, 1]
            [4, 2, 1, 3]
        """
        if self._set:
            for x in self._set:
                yield x
            return

        for x in conjugacy_class_iterator(self._part, self._domain):
            yield from_cycles(self._parent.n, x, self._parent)

    def set(self):
        r"""
        The set of all elements in the conjugacy class ``self``.

        EXAMPLES::

            sage: G = Permutations(3)
            sage: g = G([2, 1, 3])
            sage: C = G.conjugacy_class(g)
            sage: S = [[1, 3, 2], [2, 1, 3], [3, 2, 1]]
            sage: C.set() == Set(G(x) for x in S)
            True
        """
        if not self._set:
            self._set = Set(from_cycles(self._parent.n, x, self._parent)
                            for x in conjugacy_class_iterator(self._part, self._domain) )
        return self._set

#####################################################################
## Helper functions

def default_representative(part, G):
    r"""
    Construct the default representative for the conjugacy class of
    cycle type ``part`` of a symmetric group ``G``.

    Let `\lambda` be a partition of `n`. We pick a representative by

    .. MATH::

        (1, 2, \ldots, \lambda_1)
        (\lambda_1 + 1, \ldots, \lambda_1 + \lambda_2)
        (\lambda_1 + \lambda_2 + \cdots + \lambda_{\ell-1}, \ldots, n),

    where `\ell` is the length (or number of parts) of `\lambda`.

    INPUT:

    - ``part`` -- partition

    - ``G`` -- a symmetric group

    EXAMPLES::

        sage: from sage.groups.perm_gps.symgp_conjugacy_class import default_representative
        sage: S = SymmetricGroup(4)
        sage: for p in Partitions(4):
        ....:     print(default_representative(p, S))
        (1,2,3,4)
        (1,2,3)
        (1,2)(3,4)
        (1,2)
        ()
    """
    D = G.domain()
    total = 0
    cycles = []
    for p in part:
        cycles.append(tuple(D[total:total+p]))
        total += p
    return G.element_class(cycles, G, check=False)


def conjugacy_class_iterator(part, S=None):
    r"""
    Return an iterator over the conjugacy class associated to
    the partition ``part``.

    The elements are given as a list of tuples, each tuple being a cycle.

    INPUT:

    - ``part`` -- partition

    - ``S`` -- (optional, default: `\{ 1, 2, \ldots, n \}`, where `n`
      is the size of ``part``) a set

    OUTPUT:

    An iterator over the conjugacy class consisting of all
    permutations of the set ``S`` whose cycle type is ``part``.

    EXAMPLES::

        sage: from sage.groups.perm_gps.symgp_conjugacy_class import conjugacy_class_iterator
        sage: for p in conjugacy_class_iterator([2,2]): print(p)
        [(1, 2), (3, 4)]
        [(1, 4), (2, 3)]
        [(1, 3), (2, 4)]

    In order to get permutations, one just has to wrap::

        sage: S = SymmetricGroup(5)
        sage: for p in conjugacy_class_iterator([3,2]): print(S(p))
        (1,3)(2,4,5)
        (1,3)(2,5,4)
        (1,2)(3,4,5)
        (1,2)(3,5,4)
        ...
        (1,4)(2,3,5)
        (1,4)(2,5,3)

    Check that the number of elements is the number of elements in
    the conjugacy class::

        sage: s = lambda p: sum(1 for _ in conjugacy_class_iterator(p))
        sage: all(s(p) == p.conjugacy_class_size() for p in Partitions(5))
        True

    It is also possible to specify any underlying set::

        sage: it = conjugacy_class_iterator([2,2,2], 'abcdef')
        sage: sorted(flatten(next(it)))
        ['a', 'b', 'c', 'd', 'e', 'f']
        sage: all(len(x) == 2 for x in next(it))
        True
    """
    n = sum(part)
    if part not in _Partitions:
        raise ValueError("invalid partition")
    if S is None:
        S = range(1, n + 1)
    else:
        S = list(S)
        if n != len(S):
            raise ValueError("the sum of the partition %s does not match the size of %s" % (part, S))

    m = len(part)
    for s in SetPartitions(S, part):
        its = [iter(t) for t in s]
        firsts = [next(t) for t in its]
        rests = [list(t) for t in its]
        iterator = tuple(itertools.permutations(r) for r in rests)
        for r in itertools.product(*iterator):
            yield [(firsts[i],) + r[i] for i in range(m)]

