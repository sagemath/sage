r"""
Shard intersection order

This file builds a combinatorial version of the shard intersection
order of type A (in the classification of finite Coxeter groups). This
is a lattice on the set of permutations, closely related to
noncrossing partitions and the weak order.

For technical reasons, the elements of the posets are not permutations,
but can be easily converted from and to permutations::

    sage: from sage.combinat.shard_order import ShardPosetElement
    sage: p0 = Permutation([1,3,4,2])
    sage: e0 = ShardPosetElement(p0); e0
    (1, 3, 4, 2)
    sage: Permutation(list(e0)) == p0
    True

.. SEEALSO::

    A general implementation for all finite Coxeter groups is available as
    :meth:`~sage.categories.finite_coxeter_groups.FiniteCoxeterGroups.ParentMethods.shard_poset`

REFERENCES:

.. [Banc2011] \E. E. Bancroft, *Shard Intersections and Cambrian Congruence
   Classes in Type A.*, Ph.D. Thesis, North Carolina State University. 2011.

.. [Pete2013] \T. Kyle Petersen, *On the shard intersection order of
   a Coxeter group*, SIAM J. Discrete Math. 27 (2013), no. 4, 1880-1912.

.. [Read2011] \N. Reading, *Noncrossing partitions and the shard intersection
   order*, J. Algebraic Combin., 33 (2011), 483-530.
"""
from sage.combinat.posets.posets import Poset
from sage.graphs.digraph import DiGraph
from sage.combinat.permutation import Permutations


class ShardPosetElement(tuple):
    r"""
    An element of the shard poset.

    This is basically a permutation with extra stored arguments:

    - ``p`` -- the permutation itself as a tuple
    - ``runs`` -- the decreasing runs as a tuple of tuples
    - ``run_indices`` -- list; ``integer -> index of the run``
    - ``dpg`` -- the transitive closure of the shard preorder graph
    - ``spg`` -- the transitive reduction of the shard preorder graph

    These elements can easily be converted from and to permutations::

        sage: from sage.combinat.shard_order import ShardPosetElement
        sage: p0 = Permutation([1,3,4,2])
        sage: e0 = ShardPosetElement(p0); e0
        (1, 3, 4, 2)
        sage: Permutation(list(e0)) == p0
        True
    """
    def __new__(cls, p):
        r"""
        Initialization of the underlying tuple.

        TESTS::

            sage: from sage.combinat.shard_order import ShardPosetElement
            sage: ShardPosetElement(Permutation([1,3,4,2]))
            (1, 3, 4, 2)
        """
        return tuple.__new__(cls, p)

    def __init__(self, p):
        r"""
        INPUT:

        - ``p`` -- a permutation

        EXAMPLES::

            sage: from sage.combinat.shard_order import ShardPosetElement
            sage: p0 = Permutation([1,3,4,2])
            sage: e0 = ShardPosetElement(p0); e0
            (1, 3, 4, 2)
            sage: e0.dpg
            Transitive closure of : Digraph on 3 vertices
            sage: e0.spg
            Digraph on 3 vertices
        """
        self.runs = p.decreasing_runs(as_tuple=True)
        self.run_indices = [None] * (len(p) + 1)
        for i, bloc in enumerate(self.runs):
            for j in bloc:
                self.run_indices[j] = i
        G = shard_preorder_graph(self.runs)
        self.dpg = G.transitive_closure()
        self.spg = G.transitive_reduction()

    def __le__(self, other):
        """
        Comparison between two elements of the poset.

        This is the core function in the implementation of the
        shard intersection order.

        One first compares the number of runs, then the set partitions,
        then the pre-orders.

        EXAMPLES::

            sage: from sage.combinat.shard_order import ShardPosetElement
            sage: p0 = Permutation([1,3,4,2])
            sage: p1 = Permutation([1,4,3,2])
            sage: e0 = ShardPosetElement(p0)
            sage: e1 = ShardPosetElement(p1)
            sage: e0 <= e1
            True
            sage: e1 <= e0
            False

            sage: p0 = Permutation([1,2,5,7,3,4,6,8])
            sage: p1 = Permutation([2,5,7,3,4,8,6,1])
            sage: e0 = ShardPosetElement(p0)
            sage: e1 = ShardPosetElement(p1)
            sage: e0 <= e1
            True
            sage: e1 <= e0
            False
        """
        if not isinstance(other, ShardPosetElement) or len(self) != len(other):
            raise TypeError("these are not comparable")
        if self.runs == other.runs:
            return True

        # r1 must have less runs than r0
        if len(other.runs) > len(self.runs):
            return False

        dico1 = other.run_indices

        # conversion: index of run in r0 -> index of run in r1
        dico0 = [None] * len(self.runs)
        for i, bloc in enumerate(self.runs):
            j0 = dico1[bloc[0]]
            for k in bloc:
                if dico1[k] != j0:
                    return False
            dico0[i] = j0

        # at this point, the set partitions given by tuples are comparable
        dg0 = self.spg
        dg1 = other.dpg

        for i, j in dg0.edge_iterator(labels=False):
            if dico0[i] != dico0[j] and not dg1.has_edge(dico0[i], dico0[j]):
                return False
        return True


def shard_preorder_graph(runs):
    """
    Return the preorder attached to a tuple of decreasing runs.

    This is a directed graph, whose vertices correspond to the runs.

    There is an edge from a run `R` to a run `S` if `R` is before `S`
    in the list of runs and the two intervals defined by the initial and
    final indices of `R` and `S` overlap.

    This only depends on the initial and final indices of the runs.
    For this reason, this input can also be given in that shorten way.

    INPUT:

    - ``runs`` -- either

      - a tuple of tuples, the runs of a permutation, or

      - a tuple of pairs `(i,j)`, each one standing for a run from `i` to `j`

    OUTPUT: a directed graph, with vertices labelled by integers

    EXAMPLES::

        sage: from sage.combinat.shard_order import shard_preorder_graph
        sage: s = Permutation([2,8,3,9,6,4,5,1,7])
        sage: def cut(lr):
        ....:     return tuple((r[0], r[-1]) for r in lr)
        sage: shard_preorder_graph(cut(s.decreasing_runs()))
        Digraph on 5 vertices
        sage: s = Permutation([9,4,3,2,8,6,5,1,7])
        sage: P = shard_preorder_graph(s.decreasing_runs())
        sage: P.is_isomorphic(digraphs.TransitiveTournament(3))
        True
    """
    N = len(runs)
    dg = DiGraph(N)
    dg.add_edges((i, j) for i in range(N - 1)
                 for j in range(i + 1, N)
                 if runs[i][-1] < runs[j][0] and runs[j][-1] < runs[i][0])
    return dg


def shard_poset(n):
    """
    Return the shard intersection order on permutations of size `n`.

    This is defined on the set of permutations. To every permutation,
    one can attach a pre-order, using the descending runs and their
    relative positions.

    The shard intersection order is given by the implication (or refinement)
    order on the set of pre-orders defined from all permutations.

    This can also be seen in a geometrical way. Every pre-order defines
    a cone in a vector space of dimension `n`. The shard poset is given by
    the inclusion of these cones.

    .. SEEALSO::

        :func:`~sage.combinat.shard_order.shard_preorder_graph`

    EXAMPLES::

        sage: P = posets.ShardPoset(4); P  # indirect doctest
        Finite poset containing 24 elements
        sage: P.chain_polynomial()
        34*q^4 + 90*q^3 + 79*q^2 + 24*q + 1
        sage: P.characteristic_polynomial()
        q^3 - 11*q^2 + 23*q - 13
        sage: P.zeta_polynomial()
        17/3*q^3 - 6*q^2 + 4/3*q
        sage: P.is_self_dual()
        False
    """
    import operator
    Sn = [ShardPosetElement(s) for s in Permutations(n)]
    return Poset([Sn, operator.le], cover_relations=False, facade=True)
