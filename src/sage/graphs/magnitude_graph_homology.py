"""
Graph Magnitude Homology Ranks Calculations

Original maple code by Simon Willerton

Translation into Python + SageMath by James Cranch

Further modifications by Simon Willerton

EXAMPLES::

    sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
    sage: h = graphs.CycleGraph(5)
    sage: Mh = MagnitudeHomology(h, 3)
    sage: Mh.total_rank()
    [ 5  0  0  0]
    [ 0 10  0  0]
    [ 0  0 10  0]
    [ 0  0 10 10]

The program just generates the chain groups and calculates the differentials
then uses the ChainComplex package to calculate the homology.

For a graph g the chain groups MC_{*,*}(g) break up in to subcomplexes
MC_{*,l}^{s,t}(g) where l is the length of the chain and s and t are the initial
and terminal vertices of the chain.

So here generators[s,t,k,l] is a list of the degree k generators of such a chain group.
Then differential[s,t,k,l] is a matrix giving the differential from generators[s,t,k,l] to generators[s,t,k-1,l].
The homology of each subcomplex is calculated then the ranks are added together to give the required output.
"""
from six.moves import range

from sage.graphs.distances_all_pairs import distances_all_pairs
from sage.homology.chain_complex import ChainComplex
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.structure.sage_object import SageObject


class MagnitudeHomology(SageObject):
    """
    Compute the rational magnitude homology of the graph.

    INPUT:

    - g -- a graph

    - lmax -- (default 6) an integer

    The coefficient ring is `\QQ`.

    EXAMPLES::

        sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
        sage: g = Graph([(1,2),(2,3),(3,4),(4,1),(3,1)])
        sage: Mg = MagnitudeHomology(g)
        sage: Mg.total_rank()
        [  4   0   0   0   0   0   0]
        [  0  10   0   0   0   0   0]
        [  0   0  24   0   0   0   0]
        [  0   0   0  58   0   0   0]
        [  0   0   0   0 140   0   0]
        [  0   0   0   0   0 338   0]
        [  0   0   0   0   0   0 816]

        sage: h = graphs.CycleGraph(5)
        sage: Mh = MagnitudeHomology(h, 7)
        sage: Mh.total_rank()
        [ 5  0  0  0  0  0  0  0]
        [ 0 10  0  0  0  0  0  0]
        [ 0  0 10  0  0  0  0  0]
        [ 0  0 10 10  0  0  0  0]
        [ 0  0  0 30 10  0  0  0]
        [ 0  0  0  0 50 10  0  0]
        [ 0  0  0  0 20 70 10  0]
        [ 0  0  0  0  0 80 90 10]

        sage: k = graphs.PetersenGraph()
        sage: Mk = MagnitudeHomology(k)
        sage: Mk.total_rank()
        [  10    0    0    0    0    0    0]
        [   0   30    0    0    0    0    0]
        [   0    0   30    0    0    0    0]
        [   0    0  120   30    0    0    0]
        [   0    0    0  480   30    0    0]
        [   0    0    0    0  840   30    0]
        [   0    0    0    0 1440 1200   30]

    REFERENCE: :arxiv:`1505.04125`
    """
    def __init__(self, g, lmax=6):
        """
        INPUT:

        - g -- a graph

        - lmax -- an integer

        EXAMPLES::

            sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
            sage: h = graphs.CycleGraph(5)
            sage: Mh = MagnitudeHomology(h, 3)
        """
        self._basering = QQ
        self._v = g.vertices()
        self._lm = lmax
        self._km = lmax + 1
        self._g = g
        self._d = distances_all_pairs(g)

    @cached_method
    def compute_generators(self):
        """
        Populate the generators recursively.

        EXAMPLES::

            sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
            sage: g = graphs.CycleGraph(7)
            sage: M = MagnitudeHomology(g)
            sage: sorted(M.compute_generators()[(6, 1, 6, 6)])
            [(6, 0, 1, 0, 1, 0, 1),
             (6, 0, 1, 0, 1, 2, 1),
             (6, 0, 1, 0, 6, 0, 1),
             (6, 0, 1, 2, 1, 0, 1),
             (6, 0, 1, 2, 1, 2, 1),
             (6, 0, 1, 2, 3, 2, 1),
             (6, 0, 6, 0, 1, 0, 1),
             (6, 0, 6, 0, 1, 2, 1),
             (6, 0, 6, 0, 6, 0, 1),
             (6, 0, 6, 5, 6, 0, 1),
             (6, 5, 4, 5, 6, 0, 1),
             (6, 5, 6, 0, 1, 0, 1),
             (6, 5, 6, 0, 1, 2, 1),
             (6, 5, 6, 0, 6, 0, 1),
             (6, 5, 6, 5, 6, 0, 1)]
        """
        generators = {(s, t, k, l): [] for s in self._v
                      for t in self._v
                      for k in range(self._km + 2) for l in range(self._lm + 1)}

        def add_generators(a, l, x):
            k = len(a) - 1
            if k <= self._km and l <= self._lm:
                generators[(a[0], a[len(a) - 1], k, l)].append(a)
                for y in self._v:
                    if x != y:
                        add_generators(a + [y], l + self._d[x][y], y)

        for x in self._v:
            add_generators([x], 0, x)

        # number the generators, so as to produce differentials rapidly
        for s in self._v:
            for t in self._v:
                for l in range(self._lm + 1):
                    for k in range(self._km + 1):
                        generators[(s, t, k, l)] = {tuple(a): i
                                                    for i, a in enumerate(generators[(s, t, k, l)])}
        return generators

    def differential(self, s, t, k, l):
        """
        Compute the differential ?
        """
        generators = self.compute_generators()
        m = {}
        h = generators[(s, t, k - 1, l)]
        for a, i in generators[(s, t, k, l)].items():
            for z in range(len(a) - 2):
                if self._d[a[z]][a[z + 1]] + self._d[a[z + 1]][a[z + 2]] == self._d[a[z]][a[z + 2]]:
                    j = h[a[:z + 1] + a[z + 2:]]
                    if z % 2:
                        m[(j, i)] = m.get((j, i), 0) + 1
                    else:
                        m[(j, i)] = m.get((j, i), 0) - 1
        return matrix(self._basering, len(h), len(generators[(s, t, k, l)]), m)

    def chains(self, s, t, l):
        """
        Return the chain complex from vertex s to vertex t with degree l.

        EXAMPLES::

            sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
            sage: g = graphs.CycleGraph(7)
            sage: M = MagnitudeHomology(g)
            sage: M.chains(3,4,1)
            Chain complex with at most 1 nonzero terms over Rational Field
        """
        generators = self.compute_generators()
        differentials = {k: self.differential(s, t, k, l)
                         for k in range(1, self._km + 1)
                         if generators[(s, t, k, l)] or
                         generators[(s, t, k - 1, l)]}
        return ChainComplex(differentials, base_ring=self._basering, degree=-1)

    def homology(self, s, t, l):
        """
        Return the homology from vertex s to vertex t with degree l.

        EXAMPLES::

            sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
            sage: g = graphs.CycleGraph(7)
            sage: M = MagnitudeHomology(g)
            sage: M.homology(3,4,1)
            {1: Vector space of dimension 1 over Rational Field}
            sage: M.homology(0,0,1)
            {}
        """
        return self.chains(s, t, l).homology(generators=False)

    @cached_method
    def total_rank(self):
        """
        Return the matrix of ranks of rational magnitude homology groups.

        The coefficients are the sums over pairs of vertices (s, t) of
        of the ranks of rational magnitude homology groups from s to t with
        fixed degree and homological degree.

        EXAMPLES::

            sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
            sage: g = graphs.CycleGraph(7)
            sage: M = MagnitudeHomology(g)
            sage: M.total_rank()
            [ 7  0  0  0  0  0  0]
            [ 0 14  0  0  0  0  0]
            [ 0  0 14  0  0  0  0]
            [ 0  0  0 14  0  0  0]
            [ 0  0 14  0 14  0  0]
            [ 0  0  0 42  0 14  0]
            [ 0  0  0  0 70  0 14]
        """
        t_rank = matrix(ZZ, self._lm + 1, self._lm + 1, 0)

        for s in self._v:
            for t in self._v:
                for l in range(self._lm + 1):
                    for degree, group in self.homology(s, t, l).items():
                        t_rank[l, degree] += group.rank()
        return t_rank

    @cached_method
    def magnitude_function(self):
        """
        Return the series expansion of the magnitude function.

        EXAMPLES::

            sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
            sage: g = graphs.CycleGraph(7)
            sage: M = MagnitudeHomology(g)
            sage: M.magnitude_function()
            7 - 14*t + 14*t^2 - 14*t^3 + 28*t^4 - 56*t^5 + 84*t^6
        """
        t = PowerSeriesRing(ZZ, 't').gen()
        M = self.total_rank()
        return sum(sum((-1) ** k * M[l, k] for k in range(self._km)) * t ** l
                   for l in range(self._lm + 1))

    @cached_method
    def full_homology(self):
        """
        Return the complete set of ranks of rational magnitude homology groups.

        EXAMPLES::

            sage: from sage.graphs.magnitude_graph_homology import MagnitudeHomology
            sage: g = graphs.CycleGraph(7)
            sage: M = MagnitudeHomology(g)
            sage: Z = M.full_homology()
            sage: Z[3,4,1]
            {1: Vector space of dimension 1 over Rational Field}
            sage: Z[0,0,1]
            {}
        """
        return {(s, t, l): self.homology(s, t, l) for s in self._v
                for t in self._v for l in range(self._lm + 1)}
