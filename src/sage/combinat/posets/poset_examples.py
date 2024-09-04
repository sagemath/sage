# sage.doctest: needs sage.combinat sage.modules
r"""
Catalog of posets and lattices

Some common posets can be accessed through the ``posets.<tab>`` object::

    sage: posets.PentagonPoset()
    Finite lattice containing 5 elements

Moreover, the set of all posets of order `n` is represented by ``Posets(n)``::

    sage: Posets(5)
    Posets containing 5 elements

The infinite set of all posets can be used to find minimal examples::

    sage: for P in Posets():
    ....:     if not P.is_series_parallel():
    ....:         break
    sage: P
    Finite poset containing 4 elements

**Catalog of common posets:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~Posets.AntichainPoset` | Return an antichain on `n` elements.
    :meth:`~Posets.BooleanLattice` | Return the Boolean lattice on `2^n` elements.
    :meth:`~Posets.ChainPoset` | Return a chain on `n` elements.
    :meth:`~Posets.Crown` | Return the crown poset on `2n` elements.
    :meth:`~Posets.DexterSemilattice` | Return the Dexter semilattice.
    :meth:`~Posets.DiamondPoset` | Return the lattice of rank two on `n` elements.
    :meth:`~Posets.DivisorLattice` | Return the divisor lattice of an integer.
    :meth:`~Posets.DoubleTailedDiamond` | Return the double tailed diamond poset on `2n + 2` elements.
    :meth:`~Posets.IntegerCompositions` | Return the poset of integer compositions of `n`.
    :meth:`~Posets.IntegerPartitions` | Return the poset of integer partitions of ``n``.
    :meth:`~Posets.IntegerPartitionsDominanceOrder` | Return the lattice of integer partitions on the integer `n` ordered by dominance.
    :meth:`~Posets.MobilePoset` | Return the mobile poset formed by the `ribbon` with `hangers` below and an `anchor` above.
    :meth:`~Posets.NoncrossingPartitions` | Return the poset of noncrossing partitions of a finite Coxeter group ``W``.
    :meth:`~Posets.PentagonPoset` | Return the Pentagon poset.
    :meth:`~Posets.PermutationPattern` | Return the Permutation pattern poset.
    :meth:`~Posets.PermutationPatternInterval` | Return an interval in the Permutation pattern poset.
    :meth:`~Posets.PermutationPatternOccurrenceInterval` | Return the occurrence poset for a pair of comparable elements in the Permutation pattern poset.
    :meth:`~Posets.PowerPoset` | Return a power poset.
    :meth:`~Posets.ProductOfChains` | Return a product of chain posets.
    :meth:`~Posets.RandomLattice` | Return a random lattice on `n` elements.
    :meth:`~Posets.RandomPoset` | Return a random poset on `n` elements.
    :meth:`~Posets.RibbonPoset` | Return a ribbon on `n` elements with descents at `descents`.
    :meth:`~Posets.RestrictedIntegerPartitions` | Return the poset of integer partitions of `n`, ordered by restricted refinement.
    :meth:`~Posets.SetPartitions` | Return the poset of set partitions of the set `\{1,\dots,n\}`.
    :meth:`~Posets.ShardPoset` | Return the shard intersection order.
    :meth:`~Posets.SSTPoset` | Return the poset on semistandard tableaux of shape `s` and largest entry `f` that is ordered by componentwise comparison.
    :meth:`~Posets.StandardExample` | Return the standard example of a poset with dimension `n`.
    :meth:`~Posets.SymmetricGroupAbsoluteOrderPoset` | The poset of permutations with respect to absolute order.
    :meth:`~Posets.SymmetricGroupBruhatIntervalPoset` | The poset of permutations with respect to Bruhat order.
    :meth:`~Posets.SymmetricGroupBruhatOrderPoset` | The poset of permutations with respect to Bruhat order.
    :meth:`~Posets.SymmetricGroupWeakOrderPoset` | The poset of permutations of `\{ 1, 2, \ldots, n \}` with respect to the weak order.
    :meth:`~Posets.TamariLattice` | Return the Tamari lattice.
    :meth:`~Posets.TetrahedralPoset` | Return the Tetrahedral poset with `n-1` layers based on the input colors.
    :meth:`~Posets.UpDownPoset` | Return the up-down poset on `n` elements.
    :meth:`~Posets.YoungDiagramPoset` | Return the poset of cells in the Young diagram of a partition.
    :meth:`~Posets.YoungsLattice` | Return Young's Lattice up to rank `n`.
    :meth:`~Posets.YoungsLatticePrincipalOrderIdeal` | Return the principal order ideal of the partition `lam` in Young's Lattice.
    :meth:`~Posets.YoungFibonacci` | Return the Young-Fibonacci lattice up to rank `n`.

**Other available posets:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.face_lattice` | Return the face lattice of a polyhedron.
    :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron.face_lattice` | Return the face lattice of a combinatorial polyhedron.

Constructions
-------------
"""
# ****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
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

from sage.misc.classcall_metaclass import ClasscallMetaclass
import sage.categories.posets
from sage.combinat.permutation import Permutations, Permutation, to_standard
from sage.combinat.posets.posets import Poset, FinitePoset, FinitePosets_n
from sage.combinat.posets.d_complete import DCompletePoset
from sage.combinat.posets.mobile import MobilePoset as Mobile
from sage.combinat.posets.lattices import (LatticePoset, MeetSemilattice,
                                           JoinSemilattice, FiniteLatticePoset)
from sage.categories.finite_posets import FinitePosets
from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.graphs.digraph import DiGraph
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers


def check_int(n, minimum=0):
    """
    Check that ``n`` is an integer at least equal to ``minimum``.

    This is a boilerplate function ensuring input safety.

    INPUT:

    - ``n`` -- anything

    - ``minimum`` -- an optional integer (default: 0)

    EXAMPLES::

        sage: from sage.combinat.posets.poset_examples import check_int
        sage: check_int(6, 3)
        6
        sage: check_int(6)
        6

        sage: check_int(-1)
        Traceback (most recent call last):
        ...
        ValueError: number of elements must be a nonnegative integer, not -1

        sage: check_int(1, 3)
        Traceback (most recent call last):
        ...
        ValueError: number of elements must be an integer at least 3, not 1

        sage: check_int('junk')
        Traceback (most recent call last):
        ...
        ValueError: number of elements must be a nonnegative integer, not junk
    """
    if minimum == 0:
        msg = "a nonnegative integer"
    else:
        msg = f"an integer at least {minimum}"
    if n not in NonNegativeIntegers() or n < minimum:
        raise ValueError("number of elements must be " + msg + f", not {n}")
    return Integer(n)


class Posets(metaclass=ClasscallMetaclass):
    r"""
    A collection of posets and lattices.

    EXAMPLES::

        sage: posets.BooleanLattice(3)
        Finite lattice containing 8 elements
        sage: posets.ChainPoset(3)
        Finite lattice containing 3 elements
        sage: posets.RandomPoset(17,.15)
        Finite poset containing 17 elements

    The category of all posets::

        sage: Posets()
        Category of posets

    The enumerated set of all posets on `3` elements, up to an
    isomorphism::

        sage: Posets(3)
        Posets containing 3 elements

    .. SEEALSO:: :class:`~sage.categories.posets.Posets`, :class:`FinitePosets`, :func:`Poset`

    TESTS::

        sage: P = Posets
        sage: TestSuite(P).run()
    """
    @staticmethod
    def __classcall__(cls, n=None):
        r"""
        Return either the category of all posets, or the finite
        enumerated set of all finite posets on ``n`` elements up to an
        isomorphism.

        EXAMPLES::

            sage: Posets()
            Category of posets
            sage: Posets(4)
            Posets containing 4 elements
        """
        if n is None:
            return sage.categories.posets.Posets()
        n = check_int(n)
        return FinitePosets_n(n)

    @staticmethod
    def BooleanLattice(n, facade=None, use_subsets=False):
        r"""
        Return the Boolean lattice containing `2^n` elements.

        - ``n`` -- integer; number of elements will be `2^n`
        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor
        - ``use_subsets`` -- boolean (default: ``False``); if ``True``,
          then label the elements by subsets of `\{1, 2, \ldots, n\}`;
          otherwise label the elements by `0, 1, 2, \ldots, 2^n-1`

        EXAMPLES::

            sage: posets.BooleanLattice(5)
            Finite lattice containing 32 elements

            sage: sorted(posets.BooleanLattice(2))
            [0, 1, 2, 3]
            sage: sorted(posets.BooleanLattice(2, use_subsets=True), key=list)
            [{}, {1}, {1, 2}, {2}]

        TESTS:

        Check isomorphism::

            sage: B5 = posets.BooleanLattice(5)
            sage: B5S = posets.BooleanLattice(5, use_subsets=True)
            sage: B5.is_isomorphic(B5S)
            True

        Check the corner cases::

            sage: list(posets.BooleanLattice(0, use_subsets=True))
            [{}]
            sage: list(posets.BooleanLattice(1, use_subsets=True))
            [{}, {1}]
        """
        n = check_int(n)
        if n == 0:
            if use_subsets:
                from sage.sets.set import Set
                return LatticePoset(([Set()], []), facade=facade)
            return LatticePoset(([0], []), facade=facade)
        if n == 1:
            if use_subsets:
                from sage.sets.set import Set
                V = [Set(), Set([1])]
                return LatticePoset((V, [V]), facade=facade)
            return LatticePoset(([0, 1], [[0, 1]]), facade=facade)

        if use_subsets:
            from sage.sets.set import Set
            cur_level = [frozenset(range(1, n + 1))]
            D = DiGraph()
            D.add_vertex(Set(cur_level[0]))
            while cur_level:
                next_level = set()
                for X in cur_level:
                    for i in X:
                        Y = X.difference([i])
                        D.add_edge(Set(Y), Set(X))
                        next_level.add(Y)
                cur_level = next_level
            return FiniteLatticePoset(D, category=FiniteLatticePosets(),
                                      facade=facade)

        D = DiGraph({v: [Integer(v | (1 << y))
                         for y in range(n) if v & (1 << y) == 0]
                     for v in range(2**n)})
        return FiniteLatticePoset(hasse_diagram=D,
                                  category=FiniteLatticePosets(),
                                  facade=facade)

    @staticmethod
    def ChainPoset(n, facade=None):
        r"""
        Return a chain (a totally ordered poset) containing `n` elements.

        - ``n`` -- integer; number of elements
        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        EXAMPLES::

            sage: C = posets.ChainPoset(6); C
            Finite lattice containing 6 elements
            sage: C.linear_extension()
            [0, 1, 2, 3, 4, 5]

        TESTS::

            sage: for i in range(5):
            ....:     for j in range(5):
            ....:         if C.covers(C(i),C(j)) and j != i+1:
            ....:             print("TEST FAILED")

        Check that :issue:`8422` is solved::

            sage: posets.ChainPoset(0)
            Finite lattice containing 0 elements
            sage: C = posets.ChainPoset(1); C
            Finite lattice containing 1 elements
            sage: C.cover_relations()
            []
            sage: C = posets.ChainPoset(2); C
            Finite lattice containing 2 elements
            sage: C.cover_relations()
            [[0, 1]]
        """
        n = check_int(n)
        D = DiGraph([range(n), [[x, x + 1] for x in range(n - 1)]],
                    format='vertices_and_edges')
        return FiniteLatticePoset(hasse_diagram=D,
                                  category=FiniteLatticePosets(),
                                  facade=facade)

    @staticmethod
    def AntichainPoset(n, facade=None):
        """
        Return an antichain (a poset with no comparable elements)
        containing `n` elements.

        INPUT:

        - ``n`` -- integer; number of elements
        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        EXAMPLES::

            sage: A = posets.AntichainPoset(6); A
            Finite poset containing 6 elements

        TESTS::

            sage: for i in range(5):
            ....:     for j in range(5):
            ....:         if A.covers(A(i),A(j)):
            ....:             print("TEST FAILED")

        TESTS:

        Check that :issue:`8422` is solved::

            sage: posets.AntichainPoset(0)
            Finite poset containing 0 elements
            sage: C = posets.AntichainPoset(1); C
            Finite poset containing 1 elements
            sage: C.cover_relations()
            []
            sage: C = posets.AntichainPoset(2); C
            Finite poset containing 2 elements
            sage: C.cover_relations()
            []
        """
        n = check_int(n)
        return Poset((range(n), []), facade=facade)

    @staticmethod
    def PentagonPoset(facade=None):
        """
        Return the Pentagon poset.

        INPUT:

        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        EXAMPLES::

            sage: P = posets.PentagonPoset(); P
            Finite lattice containing 5 elements
            sage: P.cover_relations()
            [[0, 1], [0, 2], [1, 4], [2, 3], [3, 4]]

        TESTS:

        This is smallest lattice that is not modular::

            sage: P.is_modular()
            False

        This poset and the :meth:`DiamondPoset` are the two smallest
        lattices which are not distributive::

            sage: P.is_distributive()
            False
            sage: posets.DiamondPoset(5).is_distributive()
            False
        """
        return LatticePoset([[1, 2], [4], [3], [4], []], facade=facade)

    @staticmethod
    def DiamondPoset(n, facade=None):
        """
        Return the lattice of rank two containing ``n`` elements.

        INPUT:

        - ``n`` -- number of elements, an integer at least 3

        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        EXAMPLES::

            sage: posets.DiamondPoset(7)
            Finite lattice containing 7 elements
        """
        n = check_int(n, 3)
        c = [[n - 1] for x in range(n)]
        c[0] = list(range(1, n - 1))
        c[n - 1] = []
        D = DiGraph({v: c[v] for v in range(n)}, format='dict_of_lists')
        return FiniteLatticePoset(hasse_diagram=D,
                                  category=FiniteLatticePosets(),
                                  facade=facade)

    @staticmethod
    def Crown(n, facade=None):
        r"""
        Return the crown poset of `2n` elements.

        In this poset every element `i` for `0 \leq i \leq n-1`
        is covered by elements `i+n` and `i+n+1`, except that
        `n-1` is covered by `n` and `n+1`.

        INPUT:

        - ``n`` -- number of elements, an integer at least 2

        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        EXAMPLES::

            sage: posets.Crown(3)
            Finite poset containing 6 elements
        """
        n = check_int(n, 2)
        D = {i: [i + n, i + n + 1] for i in range(n - 1)}
        D[n - 1] = [n, n + n - 1]
        return FinitePoset(hasse_diagram=DiGraph(D), category=FinitePosets(),
                           facade=facade)

    @staticmethod
    def DivisorLattice(n, facade=None):
        """
        Return the divisor lattice of an integer.

        Elements of the lattice are divisors of `n` and `x < y` in the
        lattice if `x` divides `y`.

        INPUT:

        - ``n`` -- integer
        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        EXAMPLES::

            sage: P = posets.DivisorLattice(12)
            sage: sorted(P.cover_relations())
            [[1, 2], [1, 3], [2, 4], [2, 6], [3, 6], [4, 12], [6, 12]]

            sage: P = posets.DivisorLattice(10, facade=False)
            sage: P(2) < P(5)
            False

        TESTS::

            sage: posets.DivisorLattice(1)
            Finite lattice containing 1 elements with distinguished linear extension
        """
        from sage.arith.misc import divisors, is_prime
        n = check_int(n, 1)
        Div_n = divisors(n)
        hasse = DiGraph([Div_n, lambda a, b: b % a == 0 and is_prime(b // a)])
        return FiniteLatticePoset(hasse, elements=Div_n, facade=facade,
                                  category=FiniteLatticePosets())

    @staticmethod
    def IntegerCompositions(n):
        """
        Return the poset of integer compositions of the integer ``n``.

        A composition of a positive integer `n` is a list of positive
        integers that sum to `n`. The order is reverse refinement:
        `[p_1,p_2,...,p_l] < [q_1,q_2,...,q_m]` if `q` consists
        of an integer composition of `p_1`, followed by an integer
        composition of `p_2`, and so on.

        EXAMPLES::

            sage: P = posets.IntegerCompositions(7); P
            Finite poset containing 64 elements
            sage: len(P.cover_relations())
            192
        """
        from sage.combinat.composition import Compositions
        C = Compositions(n)
        return Poset((C, [[c, d] for c in C for d in C if d.is_finer(c)]),
                     cover_relations=False)

    @staticmethod
    def IntegerPartitions(n):
        """
        Return the poset of integer partitions on the integer ``n``.

        A partition of a positive integer `n` is a non-increasing list
        of positive integers that sum to `n`. If `p` and `q` are
        integer partitions of `n`, then `p` covers `q` if and only
        if `q` is obtained from `p` by joining two parts of `p`
        (and sorting, if necessary).

        EXAMPLES::

            sage: P = posets.IntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            28
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers
            of elements in the poset of integer partitions.
            """
            lc = []
            for i in range(len(partition) - 1):
                for j in range(i + 1, len(partition)):
                    new_partition = partition[:]
                    del new_partition[j]
                    del new_partition[i]
                    new_partition.append(partition[i] + partition[j])
                    new_partition.sort(reverse=True)
                    tup = tuple(new_partition)
                    if tup not in lc:
                        lc.append(tup)
            return lc
        from sage.combinat.partition import Partitions
        H = DiGraph(dict([[tuple(p), lower_covers(p)] for p in Partitions(n)]))
        return Poset(H.reverse())

    @staticmethod
    def RestrictedIntegerPartitions(n):
        """
        Return the poset of integer partitions on the integer `n`
        ordered by restricted refinement.

        That is, if `p` and `q` are integer partitions of `n`, then
        `p` covers `q` if and only if `q` is obtained from `p` by
        joining two distinct parts of `p` (and sorting, if necessary).

        EXAMPLES::

            sage: P = posets.RestrictedIntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            17
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers of elements in the
            restricted poset of integer partitions.
            """
            lc = []
            for i in range(len(partition) - 1):
                for j in range(i + 1, len(partition)):
                    if partition[i] != partition[j]:
                        new_partition = partition[:]
                        del new_partition[j]
                        del new_partition[i]
                        new_partition.append(partition[i] + partition[j])
                        new_partition.sort(reverse=True)
                        tup = tuple(new_partition)
                        if tup not in lc:
                            lc.append(tup)
            return lc
        from sage.combinat.partition import Partitions
        H = DiGraph(dict([[tuple(p), lower_covers(p)] for p in Partitions(n)]))
        return Poset(H.reverse())

    @staticmethod
    def IntegerPartitionsDominanceOrder(n):
        r"""
        Return the lattice of integer partitions on the integer `n`
        ordered by dominance.

        That is, if `p=(p_1,\ldots,p_i)` and `q=(q_1,\ldots,q_j)` are
        integer partitions of `n`, then `p` is greater than `q` if and
        only if `p_1+\cdots+p_k > q_1+\cdots+q_k` for all `k`.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: P = posets.IntegerPartitionsDominanceOrder(6); P
            Finite lattice containing 11 elements
            sage: P.cover_relations()
            [[[1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1]],
             [[2, 1, 1, 1, 1], [2, 2, 1, 1]],
             [[2, 2, 1, 1], [2, 2, 2]],
             [[2, 2, 1, 1], [3, 1, 1, 1]],
             [[2, 2, 2], [3, 2, 1]],
             [[3, 1, 1, 1], [3, 2, 1]],
             [[3, 2, 1], [3, 3]],
             [[3, 2, 1], [4, 1, 1]],
             [[3, 3], [4, 2]],
             [[4, 1, 1], [4, 2]],
             [[4, 2], [5, 1]],
             [[5, 1], [6]]]
        """
        n = check_int(n)
        from sage.combinat.partition import Partitions, Partition
        return LatticePoset((Partitions(n), Partition.dominates)).dual()

    @staticmethod
    def PowerPoset(n):
        r"""
        Return the power poset on `n` element posets.

        Elements of the power poset are all posets on
        the set `\{0, 1, \ldots, n-1\}` ordered by extension.
        That is, the antichain of `n` elements is the bottom and
        `P_a \le P_b` in the power poset if `P_b` is an extension
        of `P_a`.

        These were studied in [Bru1994]_.

        EXAMPLES::

            sage: P3 = posets.PowerPoset(3); P3
            Finite meet-semilattice containing 19 elements
            sage: all(P.is_chain() for P in P3.maximal_elements())
            True

        TESTS::

            sage: P0 = posets.PowerPoset(0); P0
            Finite meet-semilattice containing 1 elements
            sage: P0[0]
            Finite poset containing 0 elements
            sage: P1 = posets.PowerPoset(1); P1
            Finite meet-semilattice containing 1 elements
            sage: P1[0]
            Finite poset containing 1 elements
            sage: P1[0][0]
            0
        """
        # Todo: Make this faster.
        n = check_int(n)
        all_pos_n = set()
        Pn = list(Posets(n))
        for P in Pn:
            for r in Permutations(P):
                all_pos_n.add(P.relabel(list(r)))

        return MeetSemilattice((all_pos_n,
                                lambda A, B: all(B.is_lequal(x, y)
                                                 for x, y in A.cover_relations_iterator())))

    @staticmethod
    def ProductOfChains(chain_lengths, facade=None):
        """
        Return a product of chains.

        - ``chain_lengths`` -- list of nonnegative integers; number of
          elements in each chain

        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        EXAMPLES::

            sage: P = posets.ProductOfChains([2, 2]); P
            Finite lattice containing 4 elements
            sage: P.linear_extension()
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: P.upper_covers((0,0))
            [(0, 1), (1, 0)]
            sage: P.lower_covers((1,1))
            [(0, 1), (1, 0)]

        TESTS::

            sage: P = posets.ProductOfChains([]); P
            Finite lattice containing 0 elements
            sage: P = posets.ProductOfChains([3, 0, 1]); P
            Finite lattice containing 0 elements
            sage: P = posets.ProductOfChains([1,1,1,1]); P
            Finite lattice containing 1 elements
        """
        try:
            l = [Integer(x) for x in chain_lengths]
        except TypeError:
            raise TypeError("parameter chain_lengths must be a list of integers, not {}".format(chain_lengths))
        if any(x < 0 for x in l):
            raise TypeError("parameter chain_lengths must be a list of nonnegative integers, not {}".format(l))

        # given the empty list, we expect the empty poset.
        if not chain_lengths:
            return LatticePoset(facade=facade)
        from sage.categories.cartesian_product import cartesian_product
        elements = cartesian_product([range(i) for i in l])

        def compare(a, b):
            return all(x <= y for x, y in zip(a, b))
        return LatticePoset([elements, compare], facade=facade)

    @staticmethod
    def RandomPoset(n, p):
        r"""
        Generate a random poset on ``n`` elements according to a
        probability ``p``.

        INPUT:

        - ``n`` -- number of elements, a nonnegative integer

        - ``p`` -- a probability, a real number between 0 and 1 (inclusive)

        OUTPUT:

        A poset on `n` elements. The probability `p` roughly measures
        width/height of the output: `p=0` always generates an antichain,
        `p=1` will return a chain. To create interesting examples,
        keep the probability small, perhaps on the order of `1/n`.

        EXAMPLES::

            sage: set_random_seed(0)  # Results are reproducible
            sage: P = posets.RandomPoset(5, 0.3)
            sage: P.cover_relations()
            [[5, 4], [4, 2], [1, 2]]

        .. SEEALSO:: :meth:`RandomLattice`

        TESTS::

            sage: posets.RandomPoset(6, 'garbage')
            Traceback (most recent call last):
            ...
            TypeError: probability must be a real number, not garbage

            sage: posets.RandomPoset(6, -0.5)
            Traceback (most recent call last):
            ...
            ValueError: probability must be between 0 and 1, not -0.5

            sage: posets.RandomPoset(0, 0.5)
            Finite poset containing 0 elements
        """
        from sage.misc.prandom import random
        n = check_int(n)
        try:
            p = float(p)
        except (TypeError, ValueError):
            raise TypeError(f"probability must be a real number, not {p}")
        if p < 0 or p > 1:
            raise ValueError(f"probability must be between 0 and 1, not {p}")

        D = DiGraph(loops=False, multiedges=False)
        D.add_vertices(range(n))
        for i in range(n):
            for j in range(i + 1, n):
                if random() < p:
                    D.add_edge(i, j)
        D.relabel(list(Permutations(n).random_element()))
        return Poset(D, cover_relations=False)

    @staticmethod
    def RandomLattice(n, p, properties=None):
        r"""
        Return a random lattice on ``n`` elements.

        INPUT:

        - ``n`` -- number of elements, a nonnegative integer

        - ``p`` -- a probability, a positive real number less than one

        - ``properties`` -- list of properties for the lattice. Currently
          implemented:

          * ``None``, no restrictions for lattices to create
          * ``'planar'``, the lattice has an upward planar drawing
          * ``'dismantlable'`` (implicated by ``'planar'``)
          * ``'distributive'`` (implicated by ``'stone'``)
          * ``'stone'``

        OUTPUT:

        A lattice on `n` elements. When ``properties`` is ``None``,
        the probability `p` roughly measures number of covering
        relations of the lattice. To create interesting examples, make
        the probability a little below one, for example `0.9`.

        Currently parameter ``p`` has no effect only when ``properties``
        is not ``None``.

        .. NOTE::

            Results are reproducible in same Sage version only. Underlying
            algorithm may change in future versions.

        EXAMPLES::

            sage: set_random_seed(0)  # Results are reproducible
            sage: L = posets.RandomLattice(8, 0.995); L
            Finite lattice containing 8 elements
            sage: L.cover_relations()
            [[7, 6], [7, 3], [7, 1], ..., [5, 4], [2, 4], [1, 4], [0, 4]]
            sage: L = posets.RandomLattice(10, 0, properties=['dismantlable'])
            sage: L.is_dismantlable()
            True

        .. SEEALSO:: :meth:`RandomPoset`

        TESTS::

            sage: posets.RandomLattice(6, 'garbage')
            Traceback (most recent call last):
            ...
            TypeError: probability must be a real number, not garbage

            sage: posets.RandomLattice(6, -0.5)
            Traceback (most recent call last):
            ...
            ValueError: probability must be a positive real number and below 1, not -0.5

            sage: posets.RandomLattice(10, 0.5, properties=['junk'])
            Traceback (most recent call last):
            ...
            ValueError: unknown value junk for 'properties'

            sage: posets.RandomLattice(0, 0.5)
            Finite lattice containing 0 elements
        """
        from copy import copy
        n = check_int(n)
        try:
            p = float(p)
        except Exception:
            raise TypeError(f"probability must be a real number, not {p}")
        if p < 0 or p >= 1:
            raise ValueError("probability must be a positive real number and below 1, not {}".format(p))

        if properties is None:
            # Basic case, no special properties for lattice asked.
            if n <= 3:
                return posets.ChainPoset(n)
            covers = _random_lattice(n, p)
            covers_dict = {i: covers[i] for i in range(n)}
            D = DiGraph(covers_dict)
            D.relabel([i - 1 for i in Permutations(n).random_element()])
            return LatticePoset(D, cover_relations=True)

        if isinstance(properties, str):
            properties = {properties}
        else:
            properties = set(properties)

        known_properties = {'planar', 'dismantlable', 'distributive', 'stone'}
        errors = properties.difference(known_properties)
        if errors:
            raise ValueError("unknown value %s for 'properties'" % errors.pop())

        if n <= 3:
            # Change this, if property='complemented' is added
            return posets.ChainPoset(n)

        # Handling properties: planar => dismantlable, stone => distributive
        if 'planar' in properties:
            properties.discard('dismantlable')
        if 'stone' in properties:
            properties.discard('distributive')

        # Test property combinations that are not implemented.
        if 'distributive' in properties and len(properties) > 1:
            raise NotImplementedError("combining 'distributive' with other properties is not implemented")
        if 'stone' in properties and len(properties) > 1:
            raise NotImplementedError("combining 'stone' with other properties is not implemented")

        if properties == {'planar'}:
            D = _random_planar_lattice(n)
            D.relabel([i - 1 for i in Permutations(n).random_element()])
            return LatticePoset(D)

        if properties == {'dismantlable'}:
            D = _random_dismantlable_lattice(n)
            D.relabel([i - 1 for i in Permutations(n).random_element()])
            return LatticePoset(D)

        if properties == {'stone'}:
            D = _random_stone_lattice(n)
            D.relabel([i - 1 for i in Permutations(n).random_element()])
            return LatticePoset(D)

        if properties == {'distributive'}:
            tmp = Poset(_random_distributive_lattice(n)).order_ideals_lattice(as_ideals=False)
            D = copy(tmp._hasse_diagram)
            D.relabel([i - 1 for i in Permutations(n).random_element()])
            return LatticePoset(D)

        raise AssertionError("bug in RandomLattice()")

    @staticmethod
    def SetPartitions(n):
        r"""
        Return the lattice of set partitions of the set `\{1,\ldots,n\}`
        ordered by refinement.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: posets.SetPartitions(4)
            Finite lattice containing 15 elements
        """
        from sage.combinat.set_partition import SetPartitions
        n = check_int(n)
        S = SetPartitions(n)

        def covers(x):
            for i, s in enumerate(x):
                for j in range(i + 1, len(x)):
                    L = list(x)
                    L[i] = s.union(x[j])
                    L.pop(j)
                    yield S(L)

        return LatticePoset({x: list(covers(x)) for x in S},
                            cover_relations=True)

    @staticmethod
    def SSTPoset(s, f=None):
        """
        The lattice poset on semistandard tableaux of shape ``s`` and largest
        entry ``f`` that is ordered by componentwise comparison of the
        entries.

        INPUT:

        - ``s`` -- shape of the tableaux

        - ``f`` -- integer (default: ``None``); the maximum fill number.
          By default (``None``), the method uses the number of cells in the shape.

        .. NOTE::

            This is a basic implementation and most certainly
            not the most efficient.

        EXAMPLES::

            sage: posets.SSTPoset([2,1])
            Finite lattice containing 8 elements

            sage: posets.SSTPoset([2,1],4)
            Finite lattice containing 20 elements

            sage: posets.SSTPoset([2,1],2).cover_relations()
            [[[[1, 1], [2]], [[1, 2], [2]]]]

            sage: posets.SSTPoset([3,2]).bottom()       # long time (6s on sage.math, 2012)
            [[1, 1, 1], [2, 2]]

            sage: posets.SSTPoset([3,2],4).maximal_elements()
            [[[3, 3, 4], [4, 4]]]
        """
        from sage.combinat.tableau import SemistandardTableaux

        def tableaux_is_less_than(a, b):
            return all(ix <= iy for x, y in zip(a, b) for ix, iy in zip(x, y))

        if f is None:
            f = sum(s)
        E = SemistandardTableaux(s, max_entry=f)
        return LatticePoset((E, tableaux_is_less_than))

    @staticmethod
    def StandardExample(n, facade=None):
        r"""
        Return the partially ordered set on `2n` elements with
        dimension `n`.

        Let `P` be the poset on `\{0, 1, 2, \ldots, 2n-1\}` whose defining
        relations are that `i < j` for every `0 \leq i < n \leq j < 2n`
        except when `i + n = j`. The poset `P` is the so-called
        *standard example* of a poset with dimension `n`.

        INPUT:

        - ``n`` -- integer `\ge 2`; dimension of the constructed poset

        - ``facade`` -- boolean; whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        OUTPUT: the standard example of a poset of dimension `n`

        EXAMPLES::

            sage: A = posets.StandardExample(3); A
            Finite poset containing 6 elements
            sage: A.dimension()                                                         # needs networkx
            3

        REFERENCES:

        - [Gar2015]_
        - [Ros1999]_

        TESTS::

            sage: A = posets.StandardExample(10); A
            Finite poset containing 20 elements
            sage: len(A.cover_relations())
            90

            sage: P = posets.StandardExample(5, facade=False)
            sage: P(4) < P(3), P(4) > P(3)
            (False, False)
        """
        n = check_int(n, 2)
        return Poset((range(2 * n), [[i, j + n] for i in range(n)
                                     for j in range(n) if i != j]),
                     facade=facade)

    @staticmethod
    def SymmetricGroupBruhatOrderPoset(n):
        """
        The poset of permutations with respect to Bruhat order.

        EXAMPLES::

            sage: posets.SymmetricGroupBruhatOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10:
            element_labels = {s: "".join(str(x) for x in s)
                              for s in Permutations(n)}
        return Poset({s: s.bruhat_succ() for s in Permutations(n)},
                     element_labels)

    @staticmethod
    def SymmetricGroupBruhatIntervalPoset(start, end):
        """
        The poset of permutations with respect to Bruhat order.

        INPUT:

        - ``start`` -- list permutation

        - ``end`` -- list permutation (same n, of course)

        .. NOTE::

           Must have ``start`` <= ``end``.

        EXAMPLES:

        Any interval is rank symmetric if and only if it avoids these
        permutations::

            sage: P1 = posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4], [3,4,1,2])
            sage: P2 = posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4], [4,2,3,1])
            sage: ranks1 = [P1.rank(v) for v in P1]
            sage: ranks2 = [P2.rank(v) for v in P2]
            sage: [ranks1.count(i) for i in sorted(set(ranks1))]
            [1, 3, 5, 4, 1]
            sage: [ranks2.count(i) for i in sorted(set(ranks2))]
            [1, 3, 5, 6, 4, 1]
        """
        start = Permutation(start)
        end = Permutation(end)
        if len(start) != len(end):
            raise TypeError(f"start ({start}) and end ({end}) must have same length")
        if not start.bruhat_lequal(end):
            raise TypeError(f"must have start ({start}) <= end ({end}) in Bruhat order")
        unseen = [start]
        nodes = {}
        while unseen:
            perm = unseen.pop(0)
            nodes[perm] = [succ_perm for succ_perm in perm.bruhat_succ()
                           if succ_perm.bruhat_lequal(end)]
            unseen.extend(succ_perm for succ_perm in nodes[perm]
                          if succ_perm not in nodes)
        return Poset(nodes)

    @staticmethod
    def SymmetricGroupWeakOrderPoset(n, labels='permutations', side='right'):
        r"""
        The poset of permutations of `\{ 1, 2, \ldots, n \}` with respect
        to the weak order (also known as the permutohedron order, cf.
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`).

        The optional variable ``labels`` (default: ``'permutations'``)
        determines the labelling of the elements if `n < 10`. The optional
        variable ``side`` (default: ``'right'``) determines whether the
        right or the left permutohedron order is to be used.

        EXAMPLES::

            sage: posets.SymmetricGroupWeakOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10 and labels == "permutations":
            element_labels = dict([[s, "".join(map(str, s))]
                                   for s in Permutations(n)])
        if n < 10 and labels == "reduced_words":
            element_labels = dict([[s, "".join(map(str, s.reduced_word_lexmin()))]
                                   for s in Permutations(n)])
        if side == "left":

            def weak_covers(s):
                r"""
                Nested function for computing the covers of elements in the
                poset of left weak order for the symmetric group.
                """
                return [v for v in s.bruhat_succ() if
                        s.length() + (s.inverse().right_action_product(v)).length() == v.length()]
        else:
            def weak_covers(s):
                r"""
                Nested function for computing the covers of elements in the
                poset of right weak order for the symmetric group.
                """
                return [v for v in s.bruhat_succ() if
                        s.length() + (s.inverse().left_action_product(v)).length() == v.length()]
        return Poset(dict([[s, weak_covers(s)] for s in Permutations(n)]),
                     element_labels)

    @staticmethod
    def TetrahedralPoset(n, *colors, **labels):
        r"""
        Return the tetrahedral poset based on the input colors.

        This method will return the tetrahedral poset with `n-1` layers and
        covering relations based on the input colors of 'green', 'red',
        'orange', 'silver', 'yellow' and 'blue' as defined in [Striker2011]_.
        For particular color choices, the order ideals of the resulting
        tetrahedral poset will be isomorphic to known combinatorial objects.

        For example, for the colors 'blue', 'yellow', 'orange', and 'green',
        the order ideals will be in bijection with alternating sign matrices.
        For the colors 'yellow', 'orange', and 'green', the order ideals will
        be in bijection with semistandard Young tableaux of staircase shape.
        For the colors 'red', 'orange', 'green', and optionally 'yellow', the
        order ideals will be in bijection with totally symmetric
        self-complementary plane partitions in a `2n \times 2n \times 2n` box.

        INPUT:

        - ``n`` -- defines the number (n-1) of layers in the poset

        - ``colors`` -- the colors that define the covering relations of the
          poset; colors used are 'green', 'red', 'yellow', 'orange', 'silver',
          and 'blue'

        - ``labels`` -- keyword variable used to determine whether the poset
          is labeled with integers or tuples.  To label with integers, the
          method should be called with ``labels='integers'``.  Otherwise, the
          labeling will default to tuples.

        EXAMPLES::

            sage: posets.TetrahedralPoset(4,'green','red','yellow','silver','blue','orange')
            Finite poset containing 10 elements

            sage: posets.TetrahedralPoset(4,'green','red','yellow','silver','blue','orange',
            ....:                         labels='integers')
            Finite poset containing 10 elements

            sage: A = AlternatingSignMatrices(3)
            sage: p = A.lattice()
            sage: ji = p.join_irreducibles_poset()
            sage: tet = posets.TetrahedralPoset(3, 'green','yellow','blue','orange')
            sage: ji.is_isomorphic(tet)
            True

        TESTS::

            sage: posets.TetrahedralPoset(4,'scarlet')
            Traceback (most recent call last):
            ...
            ValueError: color input must be among: 'green', 'red', 'yellow',
            'orange', 'silver', and 'blue'
        """
        n = check_int(n, 2)
        n = n - 1
        for c in colors:
            if c not in ('green', 'red', 'yellow', 'orange', 'silver', 'blue'):
                raise ValueError("color input must be among: 'green', 'red', 'yellow', 'orange', 'silver', and 'blue'")
        elem = [(i, j, k) for i in range(n)
                for j in range(n - i) for k in range(n - i - j)]
        rels = []
        elem_labels = {}
        if 'labels' in labels:
            if labels['labels'] == 'integers':
                labelcount = 0
                for (i, j, k) in elem:
                    elem_labels[(i, j, k)] = labelcount
                    labelcount += 1
        for c in colors:
            for (i, j, k) in elem:
                if i + j + k < n - 1:
                    if c == 'green':
                        rels.append([(i, j, k), (i + 1, j, k)])
                    if c == 'red':
                        rels.append([(i, j, k), (i, j, k + 1)])
                    if c == 'yellow':
                        rels.append([(i, j, k), (i, j + 1, k)])
                if j < n - 1 and k > 0:
                    if c == 'orange':
                        rels.append([(i, j, k), (i, j + 1, k - 1)])
                if i < n - 1 and j > 0:
                    if c == 'silver':
                        rels.append([(i, j, k), (i + 1, j - 1, k)])
                if i < n - 1 and k > 0:
                    if c == 'blue':
                        rels.append([(i, j, k), (i + 1, j, k - 1)])
        return Poset([elem, rels], elem_labels)

    # shard intersection order
    import sage.combinat.shard_order
    ShardPoset = staticmethod(sage.combinat.shard_order.shard_poset)

    # Tamari lattices
    import sage.combinat.tamari_lattices
    TamariLattice = staticmethod(sage.combinat.tamari_lattices.TamariLattice)
    DexterSemilattice = staticmethod(sage.combinat.tamari_lattices.DexterSemilattice)

    @staticmethod
    def CoxeterGroupAbsoluteOrderPoset(W, use_reduced_words=True):
        r"""
        Return the poset of elements of a Coxeter group with respect
        to absolute order.

        INPUT:

        - ``W`` -- a Coxeter group
        - ``use_reduced_words`` -- boolean (default: ``True``); if
          ``True``, then the elements are labeled by their lexicographically
          minimal reduced word

        EXAMPLES::

            sage: W = CoxeterGroup(['B', 3])                                            # needs sage.groups
            sage: posets.CoxeterGroupAbsoluteOrderPoset(W)                              # needs sage.groups
            Finite poset containing 48 elements

            sage: W = WeylGroup(['B', 2], prefix='s')                                   # needs sage.groups
            sage: posets.CoxeterGroupAbsoluteOrderPoset(W, False)                       # needs sage.groups
            Finite poset containing 8 elements
        """
        if use_reduced_words:
            element_labels = {s: tuple(s.reduced_word()) for s in W}
            return Poset({s: list(s.absolute_covers()) for s in W}, element_labels)
        return Poset({s: list(s.absolute_covers()) for s in W})

    @staticmethod
    def NoncrossingPartitions(W):
        """
        Return the lattice of noncrossing partitions.

        INPUT:

        - ``W`` -- a finite Coxeter group or a Weyl group

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3])                                            # needs sage.groups
            sage: posets.NoncrossingPartitions(W)                                       # needs sage.groups
            Finite lattice containing 14 elements

            sage: W = WeylGroup(['B', 2], prefix='s')                                   # needs sage.groups
            sage: posets.NoncrossingPartitions(W)                                       # needs sage.groups
            Finite lattice containing 6 elements
        """
        return W.noncrossing_partition_lattice()

    @staticmethod
    def SymmetricGroupAbsoluteOrderPoset(n, labels='permutations'):
        r"""
        Return the poset of permutations with respect to absolute order.

        INPUT:

        - ``n`` -- a positive integer

        - ``label`` -- (default: ``'permutations'``) a label for the elements
          of the poset returned by the function; the options are

          * ``'permutations'`` -- labels the elements by their
            one-line notation
          * ``'reduced_words'`` -- labels the elements by the
            lexicographically minimal reduced word
          * ``'cycles'`` -- labels the elements by their expression
            as a product of cycles

        EXAMPLES::

            sage: posets.SymmetricGroupAbsoluteOrderPoset(4)                            # needs sage.groups
            Finite poset containing 24 elements
            sage: posets.SymmetricGroupAbsoluteOrderPoset(3, labels='cycles')           # needs sage.groups
            Finite poset containing 6 elements
            sage: posets.SymmetricGroupAbsoluteOrderPoset(3, labels='reduced_words')    # needs sage.groups
            Finite poset containing 6 elements
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        W = SymmetricGroup(n)
        if labels == "permutations":
            element_labels = {s: s.tuple() for s in W}
        if labels == "reduced_words":
            element_labels = {s: tuple(s.reduced_word()) for s in W}
        if labels == "cycles":
            element_labels = {s: "".join(x for x in s.cycle_string() if x != ',')
                              for s in W}

        return Poset({s: list(s.absolute_covers()) for s in W}, element_labels)

    @staticmethod
    def UpDownPoset(n, m=1):
        r"""
        Return the up-down poset on `n` elements where every `(m+1)`
        step is down and the rest are up.

        The case where `m=1` is sometimes referred to as the zig-zag poset
        or the fence.

        INPUT:

        - ``n`` -- nonnegative integer; number of elements in the poset
        - ``m`` -- nonnegative integer (default: 1); how frequently down
          steps occur

        OUTPUT:

        The partially ordered set on `\{ 0, 1, \ldots, n-1 \}`
        where `i` covers `i+1` if `m` divides `i+1`, and `i+1` covers `i`
        otherwise.

        EXAMPLES::

            sage: P = posets.UpDownPoset(7, 2); P
            Finite poset containing 7 elements
            sage: sorted(P.cover_relations())
            [[0, 1], [1, 2], [3, 2], [3, 4], [4, 5], [6, 5]]

        Fibonacci numbers as the number of antichains of a poset::

            sage: [len(posets.UpDownPoset(n).antichains().list()) for n in range(6)]
            [1, 2, 3, 5, 8, 13]

        TESTS::

            sage: P = posets.UpDownPoset(0); P
            Finite poset containing 0 elements
        """
        n = check_int(n)
        try:
            m = Integer(m)
        except TypeError:
            raise TypeError(f"parameter m must be an integer, not {m}")
        if m < 1:
            raise ValueError(f"parameter m must be positive, not {m}")

        covers = [[i, i + 1] if (i + 1) % (m + 1) else [i + 1, i]
                  for i in range(n - 1)]
        return Poset((range(n), covers), cover_relations=True)

    @staticmethod
    def YoungDiagramPoset(lam, dual=False):
        """
        Return the poset of cells in the Young diagram of a partition.

        INPUT:

        - ``lam`` -- a partition
        - ``dual`` -- boolean (default: ``False``); determines the orientation
          of the poset. If ``True``, then it is a join semilattice,
          otherwise it is a meet semilattice.

        EXAMPLES::

            sage: P = posets.YoungDiagramPoset(Partition([2, 2])); P
            Finite meet-semilattice containing 4 elements

            sage: sorted(P.cover_relations())
            [[(0, 0), (0, 1)], [(0, 0), (1, 0)], [(0, 1), (1, 1)], [(1, 0), (1, 1)]]

            sage: posets.YoungDiagramPoset([3, 2], dual=True)
            Finite join-semilattice containing 5 elements
        """
        from sage.combinat.partition import Partition
        lam = Partition(lam)
        if dual:
            def cell_geq(a, b):
                """
                Nested function that returns ``True`` if the cell `a` is to the
                right or below the cell `b` in the (English) Young diagram.
                """
                return ((a[0] == b[0] + 1 and a[1] == b[1]) or
                        (a[1] == b[1] + 1 and a[0] == b[0]))
            return JoinSemilattice((lam.cells(), cell_geq), cover_relations=True)
        else:
            def cell_leq(a, b):
                """
                Nested function that returns ``True`` if the cell `a` is
                to the left or above
                the cell `b` in the (English) Young diagram.
                """
                return ((a[0] == b[0] - 1 and a[1] == b[1]) or
                        (a[1] == b[1] - 1 and a[0] == b[0]))
            return MeetSemilattice((lam.cells(), cell_leq), cover_relations=True)

    @staticmethod
    def YoungsLattice(n):
        """
        Return Young's Lattice up to rank `n`.

        In other words, the poset of partitions
        of size less than or equal to `n` ordered by inclusion.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: P = posets.YoungsLattice(3); P
            Finite meet-semilattice containing 7 elements
            sage: sorted(P.cover_relations())
            [[[], [1]],
             [[1], [1, 1]],
             [[1], [2]],
             [[1, 1], [1, 1, 1]],
             [[1, 1], [2, 1]],
             [[2], [2, 1]],
             [[2], [3]]]
        """
        from sage.combinat.partition import Partitions, Partition
        from sage.misc.flatten import flatten
        partitions = flatten([list(Partitions(i)) for i in range(n + 1)])
        return JoinSemilattice((partitions, Partition.contains)).dual()

    @staticmethod
    def YoungsLatticePrincipalOrderIdeal(lam):
        """
        Return the principal order ideal of the
        partition `lam` in Young's Lattice.

        INPUT:

        - ``lam`` -- a partition

        EXAMPLES::

            sage: P = posets.YoungsLatticePrincipalOrderIdeal(Partition([2,2]))
            sage: P
            Finite lattice containing 6 elements
            sage: P.cover_relations()
            [[[], [1]],
             [[1], [1, 1]],
             [[1], [2]],
             [[1, 1], [2, 1]],
             [[2], [2, 1]],
             [[2, 1], [2, 2]]]
        """
        ideal = {}
        level = [lam]
        while level:
            new_level = set()
            for mu in level:
                down = mu.down_list()
                ideal[mu] = down
                new_level.update(down)
            level = new_level

        H = DiGraph(ideal)
        return LatticePoset(H.reverse())

    @staticmethod
    def YoungFibonacci(n):
        """
        Return the Young-Fibonacci lattice up to rank `n`.

        Elements of the (infinite) lattice are words with letters '1'
        and '2'.  The covers of a word are the words with another '1'
        added somewhere not after the first occurrence of an existing
        '1' and, additionally, the words where the first '1' is replaced by a
        '2'. The lattice is truncated to have rank `n`.

        See :wikipedia:`Young-Fibonacci lattice`.

        EXAMPLES::

            sage: Y5 = posets.YoungFibonacci(5); Y5
            Finite meet-semilattice containing 20 elements
            sage: sorted(Y5.upper_covers(Word('211')))
            [word: 1211, word: 2111, word: 221]

        TESTS::

            sage: posets.YoungFibonacci(0)
            Finite meet-semilattice containing 1 elements
            sage: posets.YoungFibonacci(1)
            Finite meet-semilattice containing 2 elements
        """
        from sage.combinat.posets.lattices import FiniteMeetSemilattice
        from sage.categories.finite_posets import FinitePosets
        from sage.combinat.words.word import Word

        n = check_int(n)

        if n == 0:
            return MeetSemilattice({'': []})

        covers = []
        current_level = ['']
        for _ in range(1, n + 1):
            new_level = set()
            for low in current_level:
                ind = low.find('1')
                if ind != -1:  # = found a '1' -> change first '1' to '2'
                    up = low[:ind] + '2' + low[ind + 1:]
                    new_level.add(up)
                    covers.append((low, up))
                else:  # no '1' in low
                    ind = len(low)

                # add '1' to every position not after first existing '1'
                for j in range(ind + 1):
                    up = '2' * j + '1' + low[j:len(low)]
                    new_level.add(up)
                    covers.append((low, up))

            current_level = new_level

        D = DiGraph([[], covers], format='vertices_and_edges')
        D.relabel(Word, inplace=True)
        return FiniteMeetSemilattice(hasse_diagram=D, category=FinitePosets())

    @staticmethod
    def DoubleTailedDiamond(n):
        r"""
        Return a double-tailed diamond of `2n + 2` elements.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: P = posets.DoubleTailedDiamond(2); P
            Finite d-complete poset containing 6 elements
            sage: P.cover_relations()
            [[1, 2], [2, 3], [2, 4], [3, 5], [4, 5], [5, 6]]
        """
        n = check_int(n, 1)

        edges = [(i, i + 1) for i in range(1, n)]
        edges.extend([(n, n + 1), (n, n + 2), (n + 1, n + 3), (n + 2, n + 3)])
        edges.extend([(i, i + 1) for i in range(n + 3, 2 * n + 2)])
        p = DiGraph([list(range(1, 2 * n + 3)), edges])
        return DCompletePoset(p)

    @staticmethod
    def PermutationPattern(n):
        r"""
        Return the poset of permutations under pattern containment
        up to rank `n`.

        INPUT:

        - ``n`` -- positive integer

        A permutation `u = u_1 \cdots u_n` contains the pattern
        `v = v_1 \cdots v_m` if there is a (not necessarily consecutive)
        subsequence of `u`  of length `m` whose entries have the same
        relative order as `v`.

        See :wikipedia:`Permutation_pattern`.

        EXAMPLES::

            sage: P4 = posets.PermutationPattern(4); P4
            Finite poset containing 33 elements
            sage: sorted(P4.lower_covers(Permutation([2,4,1,3])))
            [[1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2]]

        .. SEEALSO::

            :meth:`~sage.combinat.permutation.Permutation.has_pattern`

        TESTS::

            sage: posets.PermutationPattern(1)
            Finite poset containing 1 elements
            sage: posets.PermutationPattern(2)
            Finite poset containing 3 elements
        """
        n = check_int(n, 1)
        elem = []
        for i in range(1, n + 1):
            elem += Permutations(i)
        return Poset((elem, lambda a, b: b.has_pattern(a)))

    @staticmethod
    def PermutationPatternInterval(bottom, top):
        r"""
        Return the poset consisting of an interval in the poset of permutations
        under pattern containment between ``bottom`` and ``top``.

        INPUT:

        - ``bottom``, ``top`` -- permutations where ``top`` contains
          ``bottom`` as a pattern

        A permutation `u = u_1 \cdots u_n` contains the pattern
        `v = v_1 \cdots v_m` if there is a (not necessarily consecutive)
        subsequence of `u`  of length `m` whose entries have the same
        relative order as `v`.

        See :wikipedia:`Permutation_pattern`.

        EXAMPLES::

            sage: t = Permutation([2,3,1])
            sage: b = Permutation([4,6,2,3,5,1])
            sage: R = posets.PermutationPatternInterval(t, b); R
            Finite poset containing 14 elements
            sage: R.moebius_function(R.bottom(),R.top())
            -4

        .. SEEALSO::

            :meth:`~sage.combinat.permutation.Permutation.has_pattern`,
            :meth:`PermutationPattern`

        TESTS::

            sage: p = Permutation([1])
            sage: posets.PermutationPatternInterval(p, p)
            Finite poset containing 1 elements
        """
        P = Permutations()
        top = P(top)
        bottom = P(bottom)
        if not top.has_pattern(bottom):
            raise ValueError(f"{top} doesn't contain {bottom} as a pattern")
        # Make a list of lists of elements in the interval divided by rank.
        # List will be flattened at the end
        elem = [[top]]
        level = 0    # Consider the top element to be level 0, and then go down from there.
        rel = []     # List of covering relations to be fed into poset constructor.
        while len(top) - len(bottom) >= level + 1:
            elem.append([])  # Add a new empty level
            for upper in elem[level]:
                # Run through all permutations on current level
                #   and find relations for which it is upper cover
                upper_perm = P(upper)
                for i in range(len(top) - level):
                    # Try and remove the ith element from the permutation
                    lower = list(upper)
                    j = lower.pop(i)
                    for k in range(len(top) - level - 1):  # Standardize result
                        if lower[k] > j:
                            lower[k] = lower[k] - 1
                    lower_perm = P(lower)
                    if lower_perm.has_pattern(bottom):  # Check to see if result is in interval
                        rel += [[lower_perm, upper_perm]]
                        if lower not in elem[level + 1]:
                            elem[level + 1].append(lower_perm)
            level += 1
        elem = [item for sublist in elem for item in sublist]
        return Poset((elem, rel))

    @staticmethod
    def PermutationPatternOccurrenceInterval(bottom, top, pos):
        r"""
        Return the poset consisting of an interval in the poset of
        permutations under pattern containment between ``bottom`` and
        ``top``, where a specified instance of ``bottom`` in ``top``
        must be maintained.

        INPUT:

        - ``bottom``, ``top`` -- permutations where ``top`` contains
           ``bottom`` as a pattern
        - ``pos`` -- list of indices indicating a distinguished copy of
           ``bottom`` inside ``top`` (indexed starting at 0)

        For further information (and picture illustrating included example),
        see [ST2010]_ .

        See :wikipedia:`Permutation_pattern`.

        EXAMPLES::

            sage: t = Permutation([3,2,1])
            sage: b = Permutation([6,3,4,5,2,1])
            sage: A = posets.PermutationPatternOccurrenceInterval(t, b, (0,2,4)); A
            Finite poset containing 8 elements

        .. SEEALSO::

            :meth:`~sage.combinat.permutation.Permutation.has_pattern`,
            :meth:`PermutationPattern`, :meth:`PermutationPatternInterval`
        """
        P = Permutations()
        top = P(top)
        bottom = P(bottom)
        if not to_standard([top[z] for z in pos]) == list(bottom):  # check input
            raise ValueError("cannot find 'bottom' in 'top' given by 'pos'")
        elem = [[(top, pos)]]
        level = 0
        rel = []
        while len(top) - len(bottom) >= level + 1:
            elem.append([])  # Add a new empty level
            for upper in elem[level]:
                for i in range(len(top) - level):
                    # Try and remove the ith element from the permutation
                    if i in upper[1]:
                        continue
                    lower_perm = list(upper[0])
                    j = lower_perm.pop(i)
                    for e in range(len(top) - level - 1):
                        if lower_perm[e] > j:
                            lower_perm[e] = lower_perm[e] - 1
                    lower_pos = list(upper[1])
                    for f in range(len(upper[1])):
                        if upper[1][f] > i:
                            lower_pos[f] = upper[1][f] - 1
                    rel += [[(P(lower_perm), tuple(lower_pos)),
                             (P(upper[0]), upper[1])]]
                    if (P(lower_perm), tuple(lower_pos)) not in elem[level + 1]:
                        elem[level + 1].append((P(lower_perm), tuple(lower_pos)))
            level += 1
        elem = [item for sublist in elem for item in sublist]
        return Poset([elem, rel])

    @staticmethod
    def RibbonPoset(n, descents):
        r"""
        Return a ribbon poset on ``n`` vertices with descents at ``descents``.

        INPUT:

        - ``n`` -- the number of vertices
        - ``descents`` -- an iterable; the indices on the ribbon where `y > x`

        EXAMPLES::

            sage: R = Posets.RibbonPoset(5, [1,2])
            sage: sorted(R.cover_relations())
            [[0, 1], [2, 1], [3, 2], [3, 4]]
        """
        n = check_int(n)
        return Mobile(DiGraph([list(range(n)),
                               [(i + 1, i) if i in descents else (i, i + 1)
                                for i in range(n - 1)]]))

    @staticmethod
    def MobilePoset(ribbon, hangers, anchor=None):
        r"""
        Return a mobile poset with the ribbon ``ribbon`` and
        with hanging d-complete posets specified in ``hangers``
        and a d-complete poset attached above, specified in ``anchor``.

        INPUT:

        - ``ribbon`` -- a finite poset that is a ribbon
        - ``hangers`` -- dictionary mapping an element on the ribbon
          to a list of d-complete posets that it covers
        - ``anchor`` -- (optional) a ``tuple`` (``ribbon_elmt``,
          ``anchor_elmt``, ``anchor_poset``), where ``anchor_elmt`` covers
          ``ribbon_elmt``, and ``anchor_elmt`` is an acyclic element of
          ``anchor_poset``

        EXAMPLES::

            sage: R = Posets.RibbonPoset(5, [1,2])
            sage: H = Poset([[5, 6, 7], [(5, 6), (6,7)]])
            sage: M = Posets.MobilePoset(R, {3: [H]})
            sage: len(M.cover_relations())
            7

            sage: P = posets.MobilePoset(posets.RibbonPoset(7, [1,3]),
            ....:         {1: [posets.YoungDiagramPoset([3, 2], dual=True)],
            ....:          3: [posets.DoubleTailedDiamond(6)]},
            ....:         anchor=(4, 2, posets.ChainPoset(6)))
            sage: len(P.cover_relations())
            33
        """
        elements = []
        cover_relations = []

        cover_relations.extend(ribbon.cover_relations())
        elements.extend(ribbon._elements)

        if anchor:
            cover_relations.extend(((anchor[0], cr[0]), (anchor[0], cr[1]))
                                   for cr in anchor[2].cover_relations())
            cover_relations.append((anchor[0], (anchor[0], anchor[1])))

            elements.extend((anchor[0], elmt) for elmt in anchor[2]._elements)

        for r, hangs in hangers.items():
            for i, h in enumerate(hangs):
                elements.extend((r, i, v) for v in h._elements)
                cover_relations.extend(((r, i, cr[0]), (r, i, cr[1]))
                                       for cr in h.cover_relations())
                cover_relations.append(((r, i, h.top()), r))

        return Mobile(DiGraph([elements, cover_relations]))


# RANDOM LATTICES

# Following are helper functions for random lattice generation.
# There is no parameter checking, 0, 1, ..., n may or may not be a
# linear extension, exact output type may vary, etc. Direct use is
# discouraged. Use by posets.RandomLattice(..., properties=[...]).


def _random_lattice(n, p):
    r"""
    Return a random lattice.

    INPUT:

    - ``n`` -- number of elements, a nonnegative integer
    - ``p`` -- a number at least zero and less than one; higher number
      means more covering relations

    OUTPUT:

    A list of lists. Interpreted as a list of lower covers
    for a poset, it is a lattice with ``0..n-1`` as a linear
    extension.

    EXAMPLES::

        sage: set_random_seed(42)  # Results are reproducible
        sage: sage.combinat.posets.poset_examples._random_lattice(7, 0.4)
        [[], [0], [0], [1, 2], [1], [0], [3, 4, 5]]

    ALGORITHM::

        We add elements one by one. We check that adding a maximal
        element `e` to a meet-semilattice `L` with maximal elements
        `M` will create a semilattice by checking that there is a
        meet for `e, m` for all `m \in M`. We do that by keeping
        track of meet matrix and list of maximal elements.
    """
    from sage.arith.misc import integer_floor as floor
    from sage.misc.functional import sqrt
    from sage.misc.prandom import random

    n = n - 1
    meets = [[None] * n for _ in range(n)]
    meets[0][0] = 0
    maxs = {0}
    lc_all = [[]]  # No lower covers for the bottom element.

    for i in range(1, n):

        # Look for an admissible lower cover for the next element i
        while True:
            # Generate a random antichain
            lc_list = [i - 1 - floor(i * sqrt(random()))]
            while random() < p and 0 not in lc_list:
                new = i - 1 - floor(i * sqrt(random()))
                if any(meets[new][lc] in [new, lc] for lc in lc_list):
                    continue
                lc_list.append(new)
            # Check whether it is admissible as a new lower cover
            if all(any(all(meets[m][meets[a][a1]] == meets[m][a1] for a1 in lc_list if a1 != a) for a in lc_list) for m in maxs):
                break

        # We've found a suitable lower cover for i
        maxs.difference_update(lc_list)

        # Now compute new row and column to meet matrix.
        meets[i][i] = i
        for lc in lc_list:
            meets[i][lc] = meets[lc][i] = lc
        for e in range(i):
            meets[i][e] = meets[e][i] = max(meets[e][lc] for lc in lc_list)

        maxs.add(i)
        lc_all.append(lc_list)

    lc_all.append(list(maxs))  # Add the top element.
    return lc_all


def _random_dismantlable_lattice(n):
    r"""
    Return a random dismantlable lattice on `n` elements.

    INPUT:

    - ``n`` -- number of elements, a nonnegative integer

    OUTPUT:

    A digraph that can be interpreted as the Hasse diagram of a random
    dismantlable lattice. It has `0` as the bottom element and `n-1` as
    the top element, but otherwise `0, \ldots, n-1` *is not* usually a
    linear extension of the lattice.

    EXAMPLES::

        sage: set_random_seed(78)  # Results are reproducible
        sage: D = sage.combinat.posets.poset_examples._random_dismantlable_lattice(10); D
        Digraph on 10 vertices
        sage: D.neighbors_in(8)
        [0]

    ALGORITHM::

        We add elements one by one by "de-dismantling", i.e. select
        a random pair of comparable elements and add a new element
        between them.
    """
    from sage.misc.prandom import randint

    D = DiGraph({0: [n - 1]})
    for i in range(1, n - 1):
        a = randint(0, i // 2)
        b_ = list(D.depth_first_search(a))
        b = b_[randint(1, len(b_) - 1)]
        D.add_vertex(i)
        D.add_edge(a, i)
        D.add_edge(i, b)
        D.delete_edge(a, b)
    return D


def _random_planar_lattice(n):
    r"""
    Return a random planar lattice on `n` elements.

    INPUT:

    - ``n`` -- number of elements, a nonnegative integer

    OUTPUT:

    A random planar lattice. It has `0` as the bottom
    element and `n-1` as the top element, but otherwise
    `0, \ldots, n-1` *is not* usually a linear extension of
    the lattice.

    EXAMPLES::

        sage: set_random_seed(78)  # Results are reproducible
        sage: D = sage.combinat.posets.poset_examples._random_planar_lattice(10); D
        Digraph on 10 vertices
        sage: D.neighbors_in(8)
        [1]

    ALGORITHM::

        Every planar lattice is dismantlable.

        We add elements one by one like when generating
        dismantlable lattices, and after every addition
        check that we still have a planar lattice.
    """
    from sage.misc.prandom import randint

    G = DiGraph({0: [n - 1]})
    while G.order() < n:
        i = G.order() - 1
        a = randint(0, i // 2)
        b_ = list(G.depth_first_search(a))
        b = b_[randint(1, len(b_) - 1)]
        G1 = G.copy()
        G.add_vertex(i)
        G.add_edge(a, i)
        G.add_edge(i, b)
        G.delete_edge(a, b)
        G2 = G.copy()
        G2.add_edge(n - 1, 0)
        if not G2.is_planar():
            G = G1.copy()
    return G


def _random_distributive_lattice(n):
    """
    Return a random poset that has `n` antichains.

    INPUT:

    - ``n`` -- number of elements, a nonnegative integer

    OUTPUT:

    A random poset (as DiGraph) that has `n` antichains; i.e. a poset
    that's order ideals lattice has `n` elements.

    EXAMPLES::

        sage: g = sage.combinat.posets.poset_examples._random_distributive_lattice(10)
        sage: Poset(g).order_ideals_lattice(as_ideals=False).cardinality()
        10

    ALGORITHM:

    Add elements until there are at least `n` antichains.
    Remove elements until there are at most `n` antichains.
    Repeat.
    """
    from sage.combinat.posets.hasse_diagram import HasseDiagram
    from copy import copy
    from sage.combinat.subset import Subsets
    from sage.graphs.digraph_generators import digraphs

    if n < 4:
        return digraphs.Path(n - 1)

    H = HasseDiagram({0: []})
    while sum(1 for _ in H.antichains_iterator()) < n:
        D = copy(H)
        newcover = Subsets(H).random_element()
        new_element = H.order()
        D.add_vertex(new_element)
        for e in newcover:
            D.add_edge(e, new_element)

        D = D.transitive_reduction()
        H = HasseDiagram(D)

        while sum(1 for _ in H.antichains_iterator()) > n:
            D = copy(H)
            to_delete = H.random_vertex()
            for a in D.neighbors_in(to_delete):
                for b in D.neighbors_out(to_delete):
                    D.add_edge(a, b)
            D.delete_vertex(to_delete)
            D.relabel({z: z - 1 for z in range(to_delete + 1, D.order() + 1)})
            H = HasseDiagram(D)
    return D


def _random_stone_lattice(n):
    """
    Return a random Stone lattice on `n` elements.

    INPUT:

    - ``n`` -- number of elements, a nonnegative integer

    OUTPUT:

    A random lattice (as a digraph) of `n` elements.

    EXAMPLES::

        sage: g = sage.combinat.posets.poset_examples._random_stone_lattice(10)
        sage: LatticePoset(g).is_stone()
        True

    ALGORITHM:

    Randomly split `n` to some factors. For every factor `p` generate
    a random distributive lattice on `p-1` elements and add a new bottom
    element to it. Compute the cartesian product of those lattices.
    """
    from sage.arith.misc import factor
    from sage.combinat.partition import Partitions
    from sage.misc.misc_c import prod
    from copy import copy

    factors = sum([[f[0]] * f[1] for f in factor(n)], [])
    sage.misc.prandom.shuffle(factors)

    part_lengths = list(Partitions(len(factors)).random_element())
    parts = []
    while part_lengths:
        x = part_lengths.pop()
        parts.append(prod(factors[:x]))
        factors = factors[x:]

    result = DiGraph(1)
    for p in parts:
        g = _random_distributive_lattice(p - 1)
        g = copy(Poset(g).order_ideals_lattice(as_ideals=False)._hasse_diagram)
        g.add_edge(-1, 0)
        result = result.cartesian_product(g)
        result.relabel()

    return result


posets = Posets
