"""
Path Semigroups
"""

#*****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.semigroups import Semigroups
from sage.categories.monoids import Monoids
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.cachefunc import cached_method
from paths import QuiverPath
from representation import QuiverRep

class PathSemigroup(UniqueRepresentation, Parent):
    r"""
    The partial semigroup that is given by the directed paths of a quiver,
    subject to concatenation.

    See :mod:`~sage.quivers.representation` for a definition of this
    semigroup and of the notion of a path in a quiver.

    Note that a *partial semigroup* here is defined as a set `G` with a
    partial binary operation `G \times G \to G \cup \{\mbox{None}\}`,
    which is written infix as a `*` sign and satisfies associativity in
    the following sense: If `a`, `b` and `c` are three elements of `G`,
    and if one of the products `(a*b)*c` and `a*(b*c)` exists, then so
    does the other and the two products are equal. A partial semigroup
    is not required to have a neutral element (and this one usually has
    no such element).

    EXAMPLES::

        sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
        sage: S = Q.path_semigroup()
        sage: S
        Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices
        sage: S.variable_names()
        ('e_1', 'e_2', 'e_3', 'a', 'b', 'c', 'd')
        sage: S.gens()
        (e_1, e_2, e_3, a, b, c, d)
        sage: S.category()
        Join of Category of finite semigroups and Category of finite enumerated sets

    In the test suite, we skip the associativity test, as in this example the
    paths used for testing can't be concatenated::

        sage: TestSuite(S).run(skip=['_test_associativity'])

    If there is only a single vertex, the partial semigroup is a monoid. If
    the underlying quiver has cycles or loops, then the (partial) semigroup
    only is an infinite enumerated set. This time, there is no need to skip
    tests::

        sage: Q = DiGraph({1:{1:['a', 'b', 'c', 'd']}})
        sage: M = Q.path_semigroup()
        sage: M
        Monoid formed by the directed paths of Looped multi-digraph on 1 vertex
        sage: M.category()
        Join of Category of monoids and Category of infinite enumerated sets
        sage: TestSuite(M).run()
    """
    Element = QuiverPath

    @staticmethod
    def __classcall__(cls, Q):
        """
        Normalize the arguments passed to the constructor.

        The normalization consists of making an immutable copy of ``Q``
        that is made weighted.

        INPUT:

        - a :class:`~sage.graphs.digraph.DiGraph`.

        TESTS::

            sage: G1 = DiGraph({1:{2:['a']}})
            sage: G2 = DiGraph({1:{2:['b']}})
            sage: P1 = G1.path_semigroup()
            sage: P2 = G2.path_semigroup()
            sage: G1 == G2 # equality of unweighted graphs ignores edge labels
            True
            sage: P1.quiver() == P2.quiver() # edge labels no longer ignored
            False
            sage: P1 == P2
            False
        """
        # If self is immutable and weighted, then the copy is really cheap:
        # __copy__ just returns self.
        Q = Q.copy(immutable=True, weighted=True)
        return super(PathSemigroup, cls).__classcall__(cls, Q)

    def __init__(self, Q):
        """
        Initialize ``self``.

        INPUT:

        - a :class:`~sage.graphs.digraph.DiGraph`.

        EXAMPLES:

        Note that usually a path semigroup is created using
        :meth:`sage.graphs.digraph.DiGraph.path_semigroup`. Here, we
        demonstrate the direct construction::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}}, immutable=True)
            sage: from sage.quivers.path_semigroup import PathSemigroup
            sage: P = PathSemigroup(Q)
            sage: P is DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}}).path_semigroup() # indirect doctest
            True
            sage: P
            Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices

        While hardly of any use, it is possible to construct the path
        semigroup of an empty quiver (it is, of course, empty)::

            sage: D = DiGraph({})
            sage: A = D.path_semigroup(); A
            Partial semigroup formed by the directed paths of Digraph on 0 vertices
            sage: A.list()
            []
        """
        #########
        ## Verify that the graph labels are acceptable for this implementation ##
        # Check that edges are labelled with nonempty strings and don't begin
        # with 'e_' or contain '*'
        for x in Q.edge_labels():
            if not isinstance(x, str) or x == '':
                raise ValueError("edge labels of the digraph must be nonempty strings")
            if x[0:2] == 'e_':
                raise ValueError("edge labels of the digraph must not begin with 'e_'")
            if x.find('*') != -1:
                raise ValueError("edge labels of the digraph must not contain '*'")
        # Check validity of input: vertices have to be labelled 1,2,3,... and
        # edge labels must be unique
        for v in Q:
            if v not in ZZ:
                raise ValueError("vertices of the digraph must be labelled by integers")
        if len(set(Q.edge_labels())) != Q.num_edges():
            raise ValueError("edge labels of the digraph must be unique")
        nvert = len(Q.vertices())
        ## Determine the category which this (partial) semigroup belongs to
        if Q.is_directed_acyclic() and not Q.has_loops():
            cat = FiniteEnumeratedSets()
        else:
            cat = InfiniteEnumeratedSets()
        names = ['e_{0}'.format(v) for v in Q.vertices()] + [e[2] for e in Q.edges()]
        self._quiver = Q
        if nvert == 1:
            cat = cat.join([cat,Monoids()])
        else:
            # TODO: Make this the category of "partial semigroups"
            cat = cat.join([cat,Semigroups()])
        Parent.__init__(self, names=names, category=cat)

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q.path_semigroup()
            Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices
            sage: Q = DiGraph({1:{1:['a','b', 'c', 'd']}})
            sage: Q.path_semigroup()
            Monoid formed by the directed paths of Looped multi-digraph on 1 vertex
        """
        if self._quiver.num_verts() != 1:
            return "Partial semigroup formed by the directed paths of {}".format(self._quiver)
        return "Monoid formed by the directed paths of {}".format(self._quiver)

    def _coerce_map_from_(self, other):
        """
        A coercion from `A` to `B` exists if the underlying quiver
        of `A` is a sub-quiver of the underlying quiver of `B` (preserving
        names).

        EXAMPLES::

            sage: Q1 = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: Q2 = DiGraph({1:{2:['a'], 3:['c']}})
            sage: Q3 = DiGraph({1:{2:['a','x'], 3:['c']}, 3:{1:['d']}})
            sage: P1 = Q1.path_semigroup()
            sage: P2 = Q2.path_semigroup()
            sage: P3 = Q3.path_semigroup()
            sage: P1.has_coerce_map_from(P2)   # indirect doctest
            True
            sage: P1.has_coerce_map_from(P3)
            False
            sage: d = P1([(3,1,'d')]); d
            d
            sage: c = P2([(1,3,'c')]); c
            c
            sage: c.parent() is P1
            False
            sage: c in P1    # indirect doctest
            True
            sage: d*c        # indirect doctest
            d*c
            sage: (d*c).parent() is P1
            True
            sage: c3 = P3([(1,3,'c')]); c3
            c
            sage: c3 in P1   # indirect doctest
            False
            sage: d*c3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*':
             'Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices'
             and 'Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices'
        """
        if not isinstance(other, PathSemigroup):
            return
        # This is what we would like to do:
        #     return other.quiver().is_subgraph(self._quiver, induced=False)
        # However, this is deprecated for non-induced subgraphs
        # of directed multi-graphs.
        #
        # We ignore the deprecation and do what the deprecated method is doing
        # internally, directly using the backend to make things faster.
        sQ = self._quiver._backend
        oQ = other.quiver()._backend
        if sQ.num_verts() < oQ.num_verts():
            return False
        if any(not sQ.has_vertex(v) for v in oQ.iterator_verts(None)):
            return False
        return all(sQ.has_edge(*e) for e in oQ.iterator_out_edges(oQ.iterator_verts(None), True))

    @cached_method
    def arrows(self):
        """
        Return the elements corresponding to edges of the underlying quiver.

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).path_semigroup()
            sage: P.arrows()
            (a, b, c, d)
        """
        return tuple(self([e]) for e in self._quiver.edges())

    @cached_method
    def idempotents(self):
        """
        Return the idempotents corresponding to the vertices of the
        underlying quiver.

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).path_semigroup()
            sage: P.idempotents()
            (e_1, e_2, e_3)
        """
        return tuple(self((v,v)) for v in self._quiver.vertices())

    def ngens(self):
        """
        Return the number of generators (:meth:`arrows` and
        :meth:`idempotents`).

        EXAMPLES::

            sage: F = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).path_semigroup()
            sage: F.ngens()
            7
        """
        Q = self._quiver
        return Q.num_verts() + Q.num_edges()

    @cached_method
    def gen(self, i):
        """
        Return generator number `i`.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        An idempotent, if `i` is smaller than the number of vertices,
        or an arrow otherwise.

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).path_semigroup()
            sage: P.1         # indirect doctest
            e_2
            sage: P.idempotents()[1]
            e_2
            sage: P.5
            c
            sage: P.gens()[5]
            c
        """
        Q = self._quiver
        nv = Q.num_verts()
        if i < nv:
            v = Q.vertices()[i]
            return self((v,v))
        e = Q.edges()[i-nv]
        return self([e])

    def gens(self):
        """
        Return the tuple of generators.

        .. NOTE::

            This coincides with the sum of the output of
            :meth:`idempotents` and :meth:`arrows`.

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).path_semigroup()
            sage: P.gens()
            (e_1, e_2, e_3, a, b, c, d)
            sage: P.gens() == P.idempotents() + P.arrows()
            True
        """
        return self.idempotents() + self.arrows()

    def is_finite(self):
        """
        This partial semigroup is finite if and only if the underlying
        quiver is acyclic.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q.path_semigroup().is_finite()
            True
            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: Q.path_semigroup().is_finite()
            False
        """
        return self._quiver.is_directed_acyclic() and not self._quiver.has_loops()

    def __len__(self):
        """
        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F = Q.path_semigroup()
            sage: len(F)
            9
            sage: list(F)
            [e_1, e_2, e_3, a, b, c, d, a*d, b*d]
            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.path_semigroup()
            sage: len(F)
            Traceback (most recent call last):
            ...
            ValueError: the underlying quiver has cycles, thus, there may be an infinity of directed paths
        """
        return len(self.all_paths())

    def cardinality(self):
        """
        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F = Q.path_semigroup()
            sage: F.cardinality()
            9
            sage: A = F.algebra(QQ)
            sage: list(A.basis())
            [e_1, e_2, e_3, a, b, c, d, a*d, b*d]
            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.path_semigroup()
            sage: F.cardinality()
            +Infinity
            sage: A = F.algebra(QQ)
            sage: list(A.basis())
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list
        """
        from sage.all import ZZ
        if self._quiver.is_directed_acyclic() and not self._quiver.has_loops():
            return ZZ(len(self))
        from sage.rings.infinity import Infinity
        return Infinity

    def __iter__(self):
        """
        Iterate over the elements of ``self``, i.e., over quiver paths.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: P = Q.path_semigroup()
            sage: list(P)
            [e_1, e_2, e_3, a, b, c, d, a*d, b*d]

        The elements are sorted by length. Of course, the list of elements
        is infinite for quivers with cycles. ::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: P = Q.path_semigroup()
            sage: P.is_finite()
            False
            sage: list(P)
            Traceback (most recent call last):
            ...
            ValueError: the underlying quiver has cycles, thus, there may be an infinity of directed paths

         However, one can iterate::

            sage: counter = 0
            sage: for p in P:
            ....:     counter += 1
            ....:     print p
            ....:     if counter==20:
            ....:         break
            e_1
            e_2
            e_3
            a
            b
            d
            c
            a*d
            b*d
            d*c
            c*a
            c*b
            a*d*c
            b*d*c
            d*c*a
            d*c*b
            c*a*d
            c*b*d
            a*d*c*a
            a*d*c*b

        """
        d = 0
        length_d_available = True
        while length_d_available:
            length_d_available = False
            for v in self._quiver.vertices():
                for w in self.iter_paths_by_length_and_startpoint(d,v):
                    length_d_available = True
                    yield w
            d += 1

    def iter_paths_by_length_and_startpoint(self, d, v):
        """
        An iterator over quiver paths with a fixed length and start point.

        INPUT:

        - ``d`` -- an integer, the path length
        - ``v`` -- a vertex, start point of the paths

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: P = Q.path_semigroup()
            sage: P.is_finite()
            False
            sage: list(P.iter_paths_by_length_and_startpoint(4,1))
            [a*d*c*a, a*d*c*b, b*d*c*a, b*d*c*b]
            sage: list(P.iter_paths_by_length_and_startpoint(5,1))
            [a*d*c*a*d, a*d*c*b*d, b*d*c*a*d, b*d*c*b*d]
            sage: list(P.iter_paths_by_length_and_startpoint(5,2))
            [d*c*a*d*c, d*c*b*d*c]

        TEST::

             sage: Q = DiGraph({1:{1:['a','b', 'c', 'd']}})
             sage: P = Q.path_semigroup()
             sage: list(P.iter_paths_by_length_and_startpoint(2,1))
             [a*a,
              a*b,
              a*c,
              a*d,
              b*a,
              b*b,
              b*c,
              b*d,
              c*a,
              c*b,
              c*c,
              c*d,
              d*a,
              d*b,
              d*c,
              d*d]
             sage: len(list(P.iter_paths_by_length_and_startpoint(2,1)))
             16
        """
        # iterate over length d paths starting at vertex v
        if not d >= 0:
            raise ValueError("path length must be a non-negative integer")
        if v not in self._quiver:
            raise ValueError("the starting point {} is not a vertex of the underlying quiver".format(v))
        if not d:
            yield self([(v,v)], check=False)
        else:
            for w in self.iter_paths_by_length_and_startpoint(d-1, v):
                for a in self._quiver._backend.iterator_out_edges([w.terminal_vertex()], True):
                    yield self(list(w)+[a],check=False)

    def iter_paths_by_length_and_endpoint(self, d, v):
        """
        An iterator over quiver paths with a fixed length and end point.

        INPUT:

        - ``d`` -- an integer, the path length
        - ``v`` -- a vertex, end point of the paths

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.path_semigroup()
            sage: F.is_finite()
            False
            sage: list(F.iter_paths_by_length_and_endpoint(4,1))
            [c*a*d*c, c*b*d*c]
            sage: list(F.iter_paths_by_length_and_endpoint(5,1))
            [d*c*a*d*c, d*c*b*d*c]
            sage: list(F.iter_paths_by_length_and_endpoint(5,2))
            [c*a*d*c*a, c*b*d*c*a, c*a*d*c*b, c*b*d*c*b]
        """
        # iterate over length d paths ending at vertex v
        if not d >= 0:
            raise ValueError("path length must be a non-negative integer")
        if v not in self._quiver:
            raise ValueError("the starting point {} is not a vertex of the underlying quiver".format(v))
        if not d:
            yield self([(v,v)], check=False)
        else:
            for w in self.iter_paths_by_length_and_endpoint(d-1, v):
                for a in self._quiver._backend.iterator_in_edges([w.initial_vertex()],True):
                    yield self([a]+list(w), check=False)

    def quiver(self):
        """
        Return the underlying quiver (i.e., digraph) of this path semigroup.

        .. NOTE:

            The returned digraph always is an immutable copy of the originally
            given digraph that is made weighted.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}},
            ....:             weighted=False)
            sage: F = Q.path_semigroup()
            sage: F.quiver() == Q
            False
            sage: Q.weighted(True)
            sage: F.quiver() == Q
            True
        """
        return self._quiver

    def reverse(self):
        """
        The path semigroup of the reverse quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.path_semigroup()
            sage: F.reverse() is Q.reverse().path_semigroup()
            True
        """
        return self._quiver.reverse().path_semigroup()

    def algebra(self, k):
        """
        Return the path algebra of the underlying quiver.

        INPUT:

        - ``k`` -- a commutative ring

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: P = Q.path_semigroup()
            sage: P.algebra(GF(3))
            Path algebra of Multi-digraph on 3 vertices over Finite Field of size 3
        """
        from sage.quivers.algebra import PathAlgebra
        return PathAlgebra(k, self)

    ###########################################################################
    #                                                                         #
    # REPRESENTATION THEORETIC FUNCTIONS                                      #
    #    These functions involve the representation theory of quivers.        #
    #                                                                         #
    ###########################################################################

    def representation(self, k, *args, **kwds):
        """
        Return a representation of the quiver.

        For more information see the
        :class:`~sage.quivers.representation.QuiverRep` documentation.

        TESTS::

            sage: Q = DiGraph({1:{3:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: M
            Representation with dimension vector (2, 3, 2)
        """
        return QuiverRep(k, self, *args, **kwds)

    def S(self, k, vertex):
        """
        Return the simple module over `k` at the given vertex
        ``vertex``.

        This module is literally simple only when `k` is a field.

        INPUT:

        - `k` -- ring, the base ring of the representation

        - ``vertex`` -- integer, a vertex of the quiver

        OUTPUT:

        - :class:`~sage.quivers.representation.QuiverRep`, the simple module
          at ``vertex`` with base ring `k`

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a','b']}, 2:{3:['c','d']}}).path_semigroup()
            sage: S1 = P.S(GF(3), 1)
            sage: P.S(ZZ, 3).dimension_vector()
            (0, 0, 1)
            sage: P.S(ZZ, 1).dimension_vector()
            (1, 0, 0)

        The vertex given must be a vertex of the quiver::

            sage: P.S(QQ, 4)
            Traceback (most recent call last):
            ...
            ValueError: must specify a valid vertex of the quiver
        """
        if vertex not in self._quiver:
            raise ValueError("must specify a valid vertex of the quiver")

        # This is the module with k at the given vertex and zero elsewhere.  As
        # all maps are zero we only need to specify that the given vertex has
        # dimension 1 and the constructor will zero out everything else.
        return QuiverRep(k, self, {vertex: 1})

    simple = S

    def P(self, k, vertex):
        """
        Return the indecomposable projective module over `k` at the given
        vertex ``vertex``.

        This module is literally indecomposable only when `k` is a field.

        INPUT:

        - `k` -- ring, the base ring of the representation

        - ``vertex`` -- integer, a vertex of the quiver

        OUTPUT:

        - :class:`~sage.quivers.representation.QuiverRep`, the indecomposable
          projective module at ``vertex`` with base ring `k`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c','d']}}).path_semigroup()
            sage: P2 = Q.P(GF(3), 2)
            sage: Q.P(ZZ, 3).dimension_vector()
            (0, 0, 1)
            sage: Q.P(ZZ, 1).dimension_vector()
            (1, 2, 4)

        The vertex given must be a vertex of the quiver::

            sage: Q.P(QQ, 4)
            Traceback (most recent call last):
            ...
            ValueError: must specify a valid vertex of the quiver
        """
        if vertex not in self._quiver:
            raise ValueError("must specify a valid vertex of the quiver")
        return QuiverRep(k, self, [[(vertex, vertex)]], option='paths')

    projective = P

    def I(self, k, vertex):
        """
        Return the indecomposable injective module over `k` at the
        given vertex ``vertex``.

        This module is literally indecomposable only when `k` is a field.

        INPUT:

        - `k` -- ring, the base ring of the representation

        - ``vertex`` -- integer, a vertex of the quiver

        OUTPUT:

        - :class:`~sage.quivers.representation.QuiverRep`, the indecomposable
          injective module at vertex ``vertex`` with base ring `k`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c','d']}}).path_semigroup()
            sage: I2 = Q.I(GF(3), 2)
            sage: Q.I(ZZ, 3).dimension_vector()
            (4, 2, 1)
            sage: Q.I(ZZ, 1).dimension_vector()
            (1, 0, 0)

        The vertex given must be a vertex of the quiver::

            sage: Q.I(QQ, 4)
            Traceback (most recent call last):
            ...
            ValueError: must specify a valid vertex of the quiver
        """
        if vertex not in self._quiver:
            raise ValueError("must specify a valid vertex of the quiver")
        return QuiverRep(k, self, [[(vertex, vertex)]], option='dual paths')

    injective = I

    def free_module(self, k):
        """
        Return a free module of rank `1` over ``kP``, where `P` is
        ``self``. (In other words, the regular representation.)

        INPUT:

        - ``k`` -- ring, the base ring of the representation.

        OUTPUT:

        - :class:`~sage.quivers.representation.QuiverRep_with_path_basis`, the
          path algebra considered as a right module over itself.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}}).path_semigroup()
            sage: Q.free_module(GF(3)).dimension_vector()
            (1, 3, 6)
        """
        return QuiverRep(k, self, [[(v, v)] for v in self._quiver], option='paths')


    def all_paths(self, start=None, end=None):
        """
        List of all paths between a pair of vertices ``(start, end)``.

        INPUT:

        - ``start`` -- integer or ``None`` (default: ``None``); the initial
          vertex of the paths in the output; if ``None`` is given then
          the initial vertex is arbitrary.
        - ``end`` -- integer or ``None`` (default: ``None``); the terminal
          vertex of the paths in the output; if ``None`` is given then
          the terminal vertex is arbitrary

        OUTPUT:

        - list of paths, excluding the invalid path

        .. TODO::

            This currently does not work for quivers with cycles, even if
            there are only finitely many paths from ``start`` to ``end``.

        .. NOTE::

            If there are multiple edges between two vertices, the method
            :meth:`sage.graphs.digraph.all_paths` will not differentiate
            between them. But this method, which is not for digraphs but for
            their path semigroup associated with them, will.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F = Q.path_semigroup()
            sage: F.all_paths(1, 3)
            [a*d, b*d, c]

        If ``start=end`` then we expect only the trivial path at that vertex::

            sage: F.all_paths(1, 1)
            [e_1]

        The empty list is returned if there are no paths between the
        given vertices::

            sage: F.all_paths(3, 1)
            []

        If ``end=None`` then all edge paths beginning at ``start`` are
        returned, including the trivial path::

            sage: F.all_paths(2)
            [e_2, d]

        If ``start=None`` then all edge paths ending at ``end`` are
        returned, including the trivial path.  Note that the two edges
        from vertex 1 to vertex 2 count as two different edge paths::

            sage: F.all_paths(None, 2)
            [a, b, e_2]
            sage: F.all_paths(end=2)
            [a, b, e_2]

        If ``start=end=None`` then all edge paths are returned, including
        trivial paths::

            sage: F.all_paths()
            [e_1, a, b, a*d, b*d, c, e_2, d, e_3]

        The vertex given must be a vertex of the quiver::

            sage: F.all_paths(1, 4)
            Traceback (most recent call last):
            ...
            ValueError: the end vertex 4 is not a vertex of the quiver

        If the underlying quiver is cyclic, a ``ValueError``
        is raised::

            sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.path_semigroup()
            sage: F.all_paths()
            Traceback (most recent call last):
            ...
            ValueError: the underlying quiver has cycles, thus, there may be an infinity of directed paths

        TESTS:

        We check a single edge with a multi-character label::

            sage: Q = DiGraph([[1,2,'abc']])
            sage: PQ = Q.path_semigroup()
            sage: PQ.all_paths(1,2)
            [abc]

        An example with multiple edges::

            sage: Q = DiGraph([[1,2,'abc'], [1,2,'def']], multiedges=True)
            sage: PQ = Q.path_semigroup()
            sage: PQ.all_paths(1,2)
            [abc, def]
        """
        # Check that given arguments are vertices
        if start is not None and start not in self._quiver:
            raise ValueError("the start vertex {} is not a vertex of the quiver".format(start))
        if end is not None and end not in self._quiver:
            raise ValueError("the end vertex {} is not a vertex of the quiver".format(end))

        # Handle quivers with cycles
        Q = self._quiver
        if not (Q.is_directed_acyclic()):
            raise ValueError("the underlying quiver has cycles, thus, there may be an infinity of directed paths")

        # Handle start=None
        if start is None:
            results = []
            for v in Q:
                results += self.all_paths(v, end)
            return results

        # Handle end=None
        if end is None:
            results = []
            for v in Q:
                results += self.all_paths(start, v)
            return results

        # Handle the trivial case
        if start == end:
            return [self([(start, end)], check=False)]

        # This function will recursively convert a path given in terms of
        # vertices to a list of QuiverPaths.
        def _v_to_e(path):
            if len(path) == 1:
                return [self([(path[0], path[0])], check=False)]
            paths = []
            l = Q.edge_label(path[0], path[1])
            if isinstance(l, str):
                for b in _v_to_e(path[1:]):
                    paths.append(self([(path[0], path[1], l)] + list(b), check=False))
            else:
                for a in l:
                    for b in _v_to_e(path[1:]):
                        paths.append(self([(path[0], path[1], a)] + list(b), check=False))
            return paths

        # For each vertex path we append the resulting edge paths
        result = []
        for path in Q.all_paths(start, end):
            result += _v_to_e(path)

        # The result is all paths from start to end
        return result

