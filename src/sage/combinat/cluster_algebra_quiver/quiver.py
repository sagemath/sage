# sage.doctest: needs sage.graphs sage.modules
r"""
Quiver

A *quiver* is an oriented graph without loops, two-cycles, or multiple
edges. The edges are labelled by pairs `(i,-j)` (with `i` and `j` being
positive integers) such that the matrix `M = (m_{ab})` with
`m_{ab} = i, m_{ba} = -j` for an edge `(i,-j)` between vertices
`a` and `b` is skew-symmetrizable.

.. WARNING::

    This is not the standard definition of a quiver. Normally, in
    cluster algebra theory, a quiver is defined as an oriented graph
    without loops and two-cycles but with multiple edges allowed; the
    edges are unlabelled. This notion of quivers, however, can be seen
    as a particular case of our notion of quivers. Namely, if we have
    a quiver (in the regular sense of this word) with (precisely)
    `i` edges from `a` to `b`, then we represent it by a quiver
    (in our sense of this word) with an edge from `a` to `b` labelled
    by the pair `(i,-i)`.

For the compendium on the cluster algebra and quiver package see [MS2011]_

AUTHORS:

- Gregg Musiker
- Christian Stump

.. SEEALSO:: For mutation types of combinatorial quivers, see :meth:`~sage.combinat.cluster_algebra_quiver.quiver_mutation_type.QuiverMutationType`. Cluster seeds are closely related to :meth:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterSeed`.
"""
# ****************************************************************************
#       Copyright (C) 2011 Gregg Musiker <musiker@math.mit.edu>
#                          Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy
from itertools import product

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.graphs.views import EdgesView
from sage.arith.misc import gcd
from sage.misc.misc_c import prod
from sage.misc.lazy_import import lazy_import
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import QuiverMutationType, QuiverMutationType_Irreducible, QuiverMutationType_Reducible, _edge_list_to_matrix
from sage.combinat.cluster_algebra_quiver.mutation_class import _principal_part, _digraph_mutate, _matrix_to_digraph, _dg_canonical_form, _mutation_class_iter, _digraph_to_dig6, _dig6_to_matrix
from sage.combinat.cluster_algebra_quiver.mutation_type import _connected_mutation_type, _mutation_type_from_data, is_mutation_finite

from sage.combinat.cluster_algebra_quiver.interact import cluster_interact

lazy_import('sage.modules.free_module_element', 'vector')
lazy_import('sage.matrix.constructor', 'matrix')


class ClusterQuiver(SageObject):
    """
    The *quiver* associated to an *exchange matrix*.

    INPUT:

    - ``data`` -- can be any of the following::

      * :class:`QuiverMutationType`
      * :class:`str` -- string representing a :class:`QuiverMutationType` or a common quiver type (see Examples)
      * :class:`ClusterQuiver`
      * :class:`Matrix` -- a skew-symmetrizable matrix
      * :class:`DiGraph` -- must be the input data for a quiver
      * List of edges -- must be the edge list of a digraph for a quiver

    - ``frozen`` -- (default: ``None``) sets the list of frozen variables
      if the input type is a :class:`DiGraph`, it is ignored otherwise

    - ``user_labels`` -- (default: ``None``) sets the names of the labels for
      the vertices of the quiver if the input type is not a :class:`DiGraph`;
      otherwise it is ignored

    EXAMPLES:

    From a :class:`QuiverMutationType`::

        sage: Q = ClusterQuiver(['A',5]); Q
        Quiver on 5 vertices of type ['A', 5]

        sage: Q = ClusterQuiver(['B',2]); Q
        Quiver on 2 vertices of type ['B', 2]
        sage: Q2 = ClusterQuiver(['C',2]); Q2
        Quiver on 2 vertices of type ['B', 2]
        sage: MT = Q.mutation_type(); MT.standard_quiver() == Q
        True
        sage: MT = Q2.mutation_type(); MT.standard_quiver() == Q2
        False

        sage: Q = ClusterQuiver(['A',[2,5],1]); Q
        Quiver on 7 vertices of type ['A', [2, 5], 1]

        sage: Q = ClusterQuiver(['A', [5,0],1]); Q
        Quiver on 5 vertices of type ['D', 5]
        sage: Q.is_finite()
        True
        sage: Q.is_acyclic()
        False

        sage: Q = ClusterQuiver(['F', 4, [2,1]]); Q
        Quiver on 6 vertices of type ['F', 4, [1, 2]]
        sage: MT = Q.mutation_type(); MT.standard_quiver() == Q
        False
        sage: dg = Q.digraph(); Q.mutate([2,1,4,0,5,3])
        sage: dg2 = Q.digraph(); dg2.is_isomorphic(dg,edge_labels=True)
        False
        sage: dg2.is_isomorphic(MT.standard_quiver().digraph(),edge_labels=True)
        True

        sage: Q = ClusterQuiver(['G',2, (3,1)]); Q
        Quiver on 4 vertices of type ['G', 2, [1, 3]]
        sage: MT = Q.mutation_type(); MT.standard_quiver() == Q
        False

        sage: Q = ClusterQuiver(['GR',[3,6]]); Q
        Quiver on 4 vertices of type ['D', 4]
        sage: MT = Q.mutation_type(); MT.standard_quiver() == Q
        False

        sage: Q = ClusterQuiver(['GR',[3,7]]); Q
        Quiver on 6 vertices of type ['E', 6]

        sage: Q = ClusterQuiver(['TR',2]); Q
        Quiver on 3 vertices of type ['A', 3]
        sage: MT = Q.mutation_type(); MT.standard_quiver() == Q
        False
        sage: Q.mutate([1,0]); MT.standard_quiver() == Q
        True

        sage: Q = ClusterQuiver(['TR',3]); Q
        Quiver on 6 vertices of type ['D', 6]
        sage: MT = Q.mutation_type(); MT.standard_quiver() == Q
        False

    From a :class:`ClusterQuiver`::

        sage: Q = ClusterQuiver(['A',[2,5],1]); Q
        Quiver on 7 vertices of type ['A', [2, 5], 1]
        sage: T = ClusterQuiver( Q ); T
        Quiver on 7 vertices of type ['A', [2, 5], 1]

    From a Matrix::

        sage: Q = ClusterQuiver(['A',[2,5],1]); Q
        Quiver on 7 vertices of type ['A', [2, 5], 1]
        sage: T = ClusterQuiver( Q._M ); T
        Quiver on 7 vertices

        sage: Q = ClusterQuiver( matrix([[0,1,-1],[-1,0,1],[1,-1,0],[1,2,3]]) ); Q
        Quiver on 4 vertices with 1 frozen vertex

        sage: Q = ClusterQuiver( matrix([]) ); Q
        Quiver without vertices

    From a DiGraph::

        sage: Q = ClusterQuiver(['A',[2,5],1]); Q
        Quiver on 7 vertices of type ['A', [2, 5], 1]
        sage: T = ClusterQuiver( Q._digraph ); T
        Quiver on 7 vertices

        sage: Q = ClusterQuiver( DiGraph([[1,2],[2,3],[3,4],[4,1]]) ); Q
        Quiver on 4 vertices

        sage: Q = ClusterQuiver(DiGraph([['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'e']]),
        ....:          frozen=['c']); Q
        Quiver on 5 vertices with 1 frozen vertex
        sage: Q.mutation_type()
        [ ['A', 2], ['A', 2] ]
        sage: Q
        Quiver on 5 vertices of type [ ['A', 2], ['A', 2] ] with 1 frozen vertex

    From a List of edges::

        sage: Q = ClusterQuiver(['A',[2,5],1]); Q
        Quiver on 7 vertices of type ['A', [2, 5], 1]
        sage: T = ClusterQuiver( Q._digraph.edges(sort=True) ); T
        Quiver on 7 vertices

        sage: Q = ClusterQuiver([[1, 2], [2, 3], [3, 4], [4, 1]]); Q
        Quiver on 4 vertices

    TESTS::

        sage: Q = ClusterQuiver(DiGraph([[1,1]]))
        Traceback (most recent call last):
        ...
        ValueError: cannot add edge from 1 to 1 in graph without loops

        sage: Q = ClusterQuiver([[1,1]])
        Traceback (most recent call last):
        ...
        ValueError: cannot add edge from 1 to 1 in graph without loops

        sage: Q = ClusterQuiver(DiGraph([[1, 0],[0,1]]))
        Traceback (most recent call last):
        ...
        ValueError: the input DiGraph contains two-cycles

        sage: Q = ClusterQuiver('whatever')
        Traceback (most recent call last):
        ...
        ValueError: the input data was not recognized
    """

    def __init__(self, data, frozen=None, user_labels=None):
        """
        TESTS::

            sage: Q = ClusterQuiver(['A',4])
            sage: TestSuite(Q).run()
        """
        from sage.combinat.cluster_algebra_quiver.cluster_seed import ClusterSeed
        from sage.structure.element import Matrix

        if isinstance(user_labels, list):
            user_labels = [tuple(x) if isinstance(x, list) else x for x in user_labels]
        elif isinstance(user_labels, dict):
            values = [tuple(user_labels[x]) if isinstance(user_labels[x], list) else user_labels[x] for x in user_labels]
            user_labels = dict(zip(user_labels.keys(),
                                   values))

        # constructs a quiver from a mutation type
        if isinstance(data, (QuiverMutationType_Irreducible,
                             QuiverMutationType_Reducible)):
            if frozen is not None:
                print('The input specifies a mutation type, so the'
                      ' additional parameter frozen is ignored.'
                      ' Use set_frozen to freeze vertices.')

            mutation_type = data
            self.__init__( mutation_type.standard_quiver() )
            if user_labels:
                self.relabel(user_labels)
                self._nlist = list(user_labels)

        # constructs a quiver from string representing a mutation type
        #   or a common quiver type (see Examples)
        # NOTE: for now, any string representing a *reducible type* is
        #   coerced into the standard quiver, but there is now more flexibility
        #   in how to input a connected (irreducible) quiver.
        elif (isinstance(data, (list,tuple))
              and (isinstance(data[0], str) or all(isinstance(comp, (list,tuple))
                   and isinstance(comp[0], str) for comp in data))):
            if frozen is not None:
                print('The input specifies a mutation type, so the additional'
                      ' parameter frozen is ignored. Use set_frozen to freeze vertices.')
            mutation_type = QuiverMutationType(data)

            # The command QuiverMutationType_Irreducible (which is not imported
            #   globally) already creates the desired digraph as long as we
            #   bypass the mutation type checking of QuiverMutationType and
            #   format the input appropriately.  Thus we handle several
            #   special cases this way.
            if len(data) == 2 and isinstance(data[0], str):
                if data[0] == 'TR' or data[0] == 'GR' or (data[0] == 'C' and data[1] == 2):
                    if data[1] in ZZ:
                        quiv = ClusterQuiver( QuiverMutationType_Irreducible( data[0], data[1] )._digraph )
                        quiv._mutation_type = mutation_type
                        self.__init__( quiv )
                    elif isinstance(data[1], list):
                        quiv = ClusterQuiver( QuiverMutationType_Irreducible( data[0], tuple(data[1]) )._digraph )
                        quiv._mutation_type = mutation_type
                        self.__init__( quiv )
                else:
                    self.__init__( mutation_type.standard_quiver() )
            elif len(data) == 3 and isinstance(data[0], str):
                if (data[0] == 'F' and data[1] == 4 and data[2] == [2,1]) or (data[0] == 'G' and data[1] == 2 and data[2] == [3,1]):
                    quiv = ClusterQuiver( QuiverMutationType_Irreducible( data[0], data[1], tuple(data[2]) )._digraph )
                    quiv._mutation_type = mutation_type
                    self.__init__( quiv )
                elif (data[0] == 'F' and data[1] == 4 and data[2] == (2,1) ) or (data[0] == 'G' and data[1] == 2 and data[2] == (3,1) ):
                    quiv = ClusterQuiver( QuiverMutationType_Irreducible( data[0], data[1], data[2] )._digraph )
                    quiv._mutation_type = mutation_type
                    self.__init__( quiv )
                elif data[0] == 'A' and isinstance(data[1], list) and data[2] == 1:
                    if len(data[1]) == 2 and min(data[1]) == 0:
                        quiv = ClusterQuiver( QuiverMutationType_Irreducible( data[0], tuple(data[1]), data[2] )._digraph )
                        quiv._mutation_type = mutation_type
                        self.__init__( quiv )
                    else:
                        self.__init__( mutation_type.standard_quiver() )

                elif data[0] == 'A' and isinstance(data[1], tuple) and data[2] == 1:
                    if len(data[1]) == 2 and min(data[1]) == 0:
                        quiv = ClusterQuiver( QuiverMutationType_Irreducible( data[0], data[1], data[2] )._digraph )
                        quiv._mutation_type = mutation_type
                        self.__init__( quiv )
                    else:
                        self.__init__( mutation_type.standard_quiver() )

                else:
                    self.__init__( mutation_type.standard_quiver() )
            else:
                self.__init__( mutation_type.standard_quiver() )

            if user_labels:
                if isinstance(user_labels, dict):
                    self._nlist = list(user_labels.keys())
                else:
                    self._nlist = user_labels

                self.relabel(self._nlist)

        # constructs a quiver from a cluster seed
        elif isinstance(data, ClusterSeed):
            self.__init__( data.quiver() )

        # constructs a quiver from a quiver
        elif isinstance(data, ClusterQuiver):
            if frozen is not None:
                print('The input data is a quiver, therefore the additional'
                      ' parameter frozen is ignored. Use set_frozen to freeze vertices.')

            self._M = copy(data._M)
            self._M.set_immutable()
            self._n = data._n
            self._m = data._m
            self._mlist = list(data._mlist)
            self._nlist = list(data._nlist)
            self._digraph = copy( data._digraph )
            self._vertex_dictionary = data._vertex_dictionary
            self._mutation_type = data._mutation_type
            self._description = data._description

        # constructs a quiver from a matrix
        elif isinstance(data, Matrix):
            if not _principal_part(data).is_skew_symmetrizable( positive=True ):
                raise ValueError('The principal part of the matrix data must be skew-symmetrizable.')

            if frozen is not None:
                print('The input data is a matrix, therefore the additional parameter frozen is ignored. Frozen vertices read off accordingly if the matrix is not square.')

            self._M = copy(data).sparse_matrix()
            self._M.set_immutable()
            self._n = n = self._M.ncols()
            self._m = m = self._M.nrows() - self._n
            self._digraph = _matrix_to_digraph( self._M )
            self._vertex_dictionary = {}
            self._mutation_type = None

            if user_labels:
                if isinstance(user_labels, dict):
                    self._nlist = list(user_labels.keys())[0:n]
                    self._mlist = list(user_labels.keys())[n:n+m]
                elif isinstance(user_labels, list):
                    self._nlist = user_labels[0:n]
                    self._mlist = user_labels[n:n+m]
                self._digraph.relabel(self._nlist + self._mlist)
            else:
                self._mlist = list(range(n,n+m))
                self._nlist = list(range(n))
            if n+m == 0:
                self._description = 'Quiver without vertices'
            elif n+m == 1:
                self._description = 'Quiver on 1 vertex'
            else:
                self._description = 'Quiver on %d vertices' % (n+m)

        # constructs a quiver from a digraph
        elif isinstance(data, DiGraph):
            if frozen is None:
                frozen = []
            if not isinstance(frozen, (list, tuple)):
                raise ValueError("'frozen' must be a list of vertices")
            frozen = set(frozen)
            if not frozen.issubset(data.vertex_iterator()):
                raise ValueError("frozen elements must be vertices")

            mlist = self._mlist = list(frozen)
            m = self._m = len(mlist)

            try:
                nlist = sorted(x for x in data if x not in frozen)
            except TypeError:
                nlist = sorted([x for x in data if x not in frozen], key=str)
            self._nlist = nlist
            n = self._n = len(nlist)
            labelDict = {x: i for i, x in enumerate(nlist + mlist)}

            dg = copy(data)
            if data.has_loops():
                raise ValueError("the input DiGraph contains a loop")

            edges = set(data.edge_iterator(labels=False))
            if any((b, a) in edges for (a, b) in edges):
                raise ValueError("the input DiGraph contains two-cycles")

            dg_labelling = False
            if not set(dg.vertex_iterator()) == set(range(n + m)):
                # relabelling to integers
                # frozen vertices must be preserved
                dg_labelling = nlist + mlist
                dg.relabel(labelDict)

            multiple_edges = dg.multiple_edges()
            if multiple_edges:
                multi_edges = {}
                for v1, v2, label in multiple_edges:
                    if label not in ZZ:
                        raise ValueError("the input DiGraph contains multiple"
                                         " edges labeled by non-integers")
                    elif (v1, v2) in multi_edges:
                        multi_edges[(v1, v2)] += label
                    else:
                        multi_edges[(v1, v2)] = label
                    dg.delete_edge(v1, v2)
                dg.add_edges([(v1, v2, multi_edges[(v1,v2)])
                              for v1, v2 in multi_edges])

            for edge in dg.edge_iterator():
                if edge[0] >= n and edge[1] >= n:
                    raise ValueError("the input digraph contains edges"
                                     " within the frozen vertices")
                if edge[2] is None:
                    dg.set_edge_label( edge[0], edge[1], (1,-1) )
                    edge = (edge[0],edge[1],(1,-1))
                elif edge[2] in ZZ:
                    dg.set_edge_label( edge[0], edge[1], (edge[2],-edge[2]) )
                    edge = (edge[0],edge[1],(edge[2],-edge[2]))
                elif isinstance(edge[2], list) and len(edge[2]) != 2:
                    raise ValueError("the input digraph contains an edge with"
                                     " the wrong type of list as a label")
                elif isinstance(edge[2], list) and len(edge[2]) == 2:
                    dg.set_edge_label( edge[0], edge[1], (edge[2][0], edge[2][1]))
                    edge = (edge[0],edge[1],(edge[2][0],edge[2][1]))
                elif ( edge[0] >= n or edge[1] >= n ) and not edge[2][0] == - edge[2][1]:
                    raise ValueError("the input digraph contains an edge to or"
                                     " from a frozen vertex which is not skew-symmetric")
                if edge[2][0] < 0:
                    raise ValueError("the input digraph contains an edge of the"
                                     " form (a,-b) with negative a")

            M = _edge_list_to_matrix( dg.edge_iterator(), list(range(n)),
                                      list(range(n, n + m)) )
            if not _principal_part(M).is_skew_symmetrizable(positive=True):
                raise ValueError("the input digraph must be skew-symmetrizable")

            self._digraph = dg
            self._vertex_dictionary = {}
            if dg_labelling is not False:
                self.relabel(dg_labelling)  # relabelling back

            self._M = M
            self._M.set_immutable()
            if n + m == 0:
                self._description = 'Quiver without vertices'
            elif n + m == 1:
                self._description = 'Quiver on %d vertex' % (n+m)
            else:
                self._description = 'Quiver on %d vertices' % (n+m)
            self._mutation_type = None

        # if data is a list of edges, the appropriate digraph is constructed.

        elif isinstance(data, (list, EdgesView)) and all(isinstance(x, (list,tuple)) for x in data):
            dg = DiGraph(data)
            self.__init__(data=dg, frozen=frozen)

        # otherwise, an error is raised
        else:
            raise ValueError("the input data was not recognized")

        # stopgap for bugs arising from coefficients
        if self._m:
            from sage.misc.stopgap import stopgap
            stopgap("Having frozen nodes is known to produce wrong answers", 22381)

    def __eq__(self, other):
        """
        Return ``True`` if ``self`` and ``other`` represent the same quiver.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: T = Q.mutate( 2, inplace=False )
            sage: Q.__eq__( T )
            False
            sage: T.mutate( 2 )
            sage: Q.__eq__( T )
            True
        """
        return isinstance(other, ClusterQuiver) and self._M == other._M

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: hash(Q)  # indirect doctest
            7654921743699262111  # 64-bit
            -1264862561          # 32-bit
        """
        return hash(self._M)

    def _repr_(self):
        """
        Return the description of ``self``.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q._repr_()
            "Quiver on 5 vertices of type ['A', 5]"
        """
        name = self._description
        if self._mutation_type:
            if isinstance(self._mutation_type, str):
                name += ' of ' + self._mutation_type
            else:
                name += ' of type ' + str(self._mutation_type)
        if self._m == 1:
            name += ' with %s frozen vertex' % self._m
        elif self._m > 1:
            name += ' with %s frozen vertices' % self._m
        return name

    def plot(self, circular=True, center=(0, 0), directed=True, mark=None,
             save_pos=False, greens=[]):
        """
        Return the plot of the underlying digraph of ``self``.

        INPUT:

        - ``circular`` -- boolean (default: ``True``); if ``True``, the
          circular plot is chosen, otherwise >>spring<< is used
        - ``center`` -- (default: (0,0)) sets the center of the circular plot,
          otherwise it is ignored
        - ``directed`` -- boolean (default: ``True``); if ``True``, the
          directed version is shown, otherwise the undirected
        - ``mark`` -- (default: ``None``) if set to i, the vertex i is
          highlighted
        - ``save_pos`` -- boolean (default: ``False``); if ``True``, the
          positions of the vertices are saved
        - ``greens`` -- (default: ``[]``) if set to a list, will display the green
          vertices as green

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.plot()                                                              # needs sage.plot sage.symbolic
            Graphics object consisting of 15 graphics primitives
            sage: Q.plot(circular=True)                                                 # needs sage.plot sage.symbolic
            Graphics object consisting of 15 graphics primitives
            sage: Q.plot(circular=True, mark=1)                                         # needs sage.plot sage.symbolic
            Graphics object consisting of 15 graphics primitives
        """
        from sage.plot.colors import rainbow
        from sage.graphs.graph_generators import GraphGenerators
        from sage.symbolic.constants import e, pi
        from sage.rings.cc import CC
        from sage.rings.imaginary_unit import I
        graphs = GraphGenerators()
        # returns positions for graph vertices on two concentric cycles with radius 1 and 2

        def _graphs_concentric_circles(n, m):
            g1 = graphs.CycleGraph(n).get_pos()
            g2 = graphs.CycleGraph(m).get_pos()
            for i in g2:
                z = CC(g2[i])*e**(-pi*I/(2*m))
                g2[i] = (z.real_part(),z.imag_part())
            for i in range(m):
                g1[n+i] = [2*g2[i][0], 2*g2[i][1]]
            return g1

        n, m = self._n, self._m
        # So that we don't remove elements of these lists later
        nlist = copy(self._nlist)
        mlist = copy(self._mlist)
        colors = rainbow(11)
        color_dict = { colors[0]:[], colors[1]:[], colors[6]:[], colors[5]:[] }

        # Set up our graph. If it's directed we have a digraph, else just a normal graph
        if directed:
            dg = DiGraph( self._digraph )
        else:
            dg = Graph( self._digraph )

        # For each edge in our graph we assign a color
        for edge in dg.edges(sort=True):
            v1,v2,(a,b) = edge

            if v1 in nlist and v2 in nlist:
                if (a,b) == (1,-1):
                    color_dict[ colors[0] ].append((v1,v2))
                else:
                    color_dict[ colors[6] ].append((v1,v2))
            else:
                if (a,b) == (1,-1):
                    color_dict[ colors[1] ].append((v1,v2))
                else:
                    color_dict[ colors[5] ].append((v1,v2))
            if a == -b:
                if a == 1:
                    dg.set_edge_label(v1, v2, '')
                else:
                    dg.set_edge_label(v1, v2, a)

        # If a mark is given, then we set that mark apart from the rest
        # The mark is assumed to be a vertex
        if mark is not None:

            if mark in nlist:
                nlist.remove(mark)
                partition = (nlist, mlist, [mark])
            elif mark in mlist:
                mlist.remove(mark)
                partition = (nlist, mlist, [mark])
            else:
                raise ValueError("The given mark is not a vertex of self.")
        else:

            # Partition out the green vertices
            for i in greens:
                if i in nlist:
                    nlist.remove(i)
                else:
                    mlist.remove(i)
            partition = (nlist, mlist, greens)

        vertex_color_dict = {'tomato': partition[0],
                             'lightblue': partition[1],
                             'lightgreen': partition[2]}

        options = {
            'graph_border': True,
            'edge_colors': color_dict,
            'vertex_colors': vertex_color_dict,
            'edge_labels': True,
            'vertex_labels': True,
        }
        if circular:
            pp = _graphs_concentric_circles( n, m )
            options['pos'] = {}
            for v in pp:
                # If we're using vertex dictionary set that as key
                if v in self._vertex_dictionary:
                    vkey = self._vertex_dictionary[v]
                else:
                    vkey = v
                options['pos'][vkey] = (pp[v][0] + center[0], pp[v][1] + center[1])

        return dg.plot( **options )

    def show(self, fig_size=1, circular=False, directed=True, mark=None, save_pos=False, greens=[]):
        """
        Show the plot of the underlying digraph of ``self``.

        INPUT:

        - ``fig_size`` -- (default: 1) factor by which the size of the plot
          is multiplied
        - ``circular`` -- boolean (default: ``False``); if ``True``, the
          circular plot is chosen, otherwise >>spring<< is used
        - ``directed`` -- boolean (default: ``True``); if ``True``, the directed
          version is shown, otherwise the undirected
        - ``mark`` -- boolean (default: ``None``); if set to i, the vertex i is
          highlighted
        - ``save_pos`` -- boolean (default: ``False``); if ``True``, the
          positions of the vertices are saved
        - ``greens`` -- (default: ``[]``) if set to a list, will display the
          green vertices as green

        TESTS::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.show() # long time
        """
        n, m = self._n, self._m
        plot = self.plot( circular=circular, directed=directed, mark=mark, save_pos=save_pos, greens=greens)
        if circular:
            plot.show( figsize=[fig_size*3*(n+m)/4+1,fig_size*3*(n+m)/4+1] )
        else:
            plot.show( figsize=[fig_size*n+1,fig_size*n+1] )

    def interact(self, fig_size=1, circular=True):
        r"""
        Start an interactive window for cluster quiver mutations.

        Only in *Jupyter notebook mode*.

        INPUT:

        - ``fig_size`` -- (default: 1) factor by which the size of the
          plot is multiplied

        - ``circular`` -- boolean (default: ``True``); if ``True``, the
          circular plot is chosen, otherwise >>spring<< is used

        TESTS::

            sage: S = ClusterQuiver(['A',4])
            sage: S.interact()                                                          # needs sage.plot sage.symbolic
            ...VBox(children=...
        """
        return cluster_interact(self, fig_size, circular, kind='quiver')

    def save_image(self, filename, circular=False):
        """
        Save the plot of the underlying digraph of ``self``.

        INPUT:

        - ``filename`` -- the filename the image is saved to
        - ``circular`` -- boolean (default: ``False``); if ``True``, the
          circular plot is chosen, otherwise >>spring<< is used

        EXAMPLES::

            sage: Q = ClusterQuiver(['F',4,[1,2]])
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:                 # needs sage.plot sage.symbolic
            ....:     Q.save_image(f.name)
        """
        graph_plot = self.plot(circular=circular)
        graph_plot.save(filename=filename)

    def qmu_save(self, filename=None):
        """
        Save ``self`` in a ``.qmu`` file.

        This file can then be opened in Bernhard Keller's Quiver Applet.

        See https://webusers.imj-prg.fr/~bernhard.keller/quivermutation/

        INPUT:

        - ``filename`` -- the filename the image is saved to

        If a filename is not specified, the default name is
        ``from_sage.qmu`` in the current sage directory.

        EXAMPLES::

            sage: Q = ClusterQuiver(['F',4,[1,2]])
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.qmu') as f:                 # needs sage.plot sage.symbolic
            ....:     Q.qmu_save(f.name)

        Make sure we can save quivers with `m != n` frozen variables, see :issue:`14851`::

            sage: S = ClusterSeed(['A',3])
            sage: T1 = S.principal_extension()
            sage: Q = T1.quiver()
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.qmu') as f:                 # needs sage.plot sage.symbolic
            ....:     Q.qmu_save(f.name)
        """
        M = self.b_matrix()
        if self.m():
            from sage.matrix.constructor import matrix
            from sage.matrix.constructor import block_matrix
            M1 = M.matrix_from_rows(list(range(self.n())))
            M2 = M.matrix_from_rows(list(range(self.n(), self.n() + self.m())))
            M3 = matrix(self.m(), self.m())
            M = block_matrix([[M1, -M2.transpose()], [M2, M3]])
        dg = self.digraph()
        dg.plot(save_pos=True)
        PP = dg.get_pos()
        m = M.ncols()

        if filename is None:
            filename = 'from_sage.qmu'
        try:
            self._default_filename = filename
        except AttributeError:
            pass
        if filename[-4:] != '.qmu':
            filename += '.qmu'

        string = []
        string.append('//Number of points')
        string.append(str(m))
        string.append('//Vertex radius')
        string.append('9')
        string.append('//Labels shown')
        string.append('1')
        string.append('//Matrix')
        string.append(str(m) + ' ' + str(m))
        string.extend(' '.join(str(M[i, j]) for j in range(m))
                      for i in range(m))
        string.append('//Points')

        for i in range(m):
            x, y = PP[i]
            txt = '9 ' + str(100 * x) + ' ' + str(100 * y)
            if i >= self.n():
                txt += ' 1'
            string.append(txt)

        string.append('//Historycounter')
        string.append('-1')
        string.append('//History')
        string.append('')
        string.append('//Cluster is null')

        string = '\n'.join(string)

        with open(filename, 'w') as myfile:
            myfile.write(string)

    def b_matrix(self):
        """
        Return the b-matrix of ``self``.

        EXAMPLES::

            sage: ClusterQuiver(['A',4]).b_matrix()
            [ 0  1  0  0]
            [-1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -1  0]

            sage: ClusterQuiver(['B',4]).b_matrix()
            [ 0  1  0  0]
            [-1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -2  0]

            sage: ClusterQuiver(['D',4]).b_matrix()
            [ 0  1  0  0]
            [-1  0 -1 -1]
            [ 0  1  0  0]
            [ 0  1  0  0]

            sage: ClusterQuiver(QuiverMutationType([['A',2],['B',2]])).b_matrix()
            [ 0  1  0  0]
            [-1  0  0  0]
            [ 0  0  0  1]
            [ 0  0 -2  0]
        """
        return copy(self._M)

    def digraph(self):
        """
        Return the underlying digraph of ``self``.

        EXAMPLES::

            sage: ClusterQuiver(['A',1]).digraph()
            Digraph on 1 vertex
            sage: list(ClusterQuiver(['A',1]).digraph())
            [0]
            sage: ClusterQuiver(['A',1]).digraph().edges(sort=True)
            []

            sage: ClusterQuiver(['A',4]).digraph()
            Digraph on 4 vertices
            sage: ClusterQuiver(['A',4]).digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]

            sage: ClusterQuiver(['B',4]).digraph()
            Digraph on 4 vertices
            sage: ClusterQuiver(['A',4]).digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]

            sage: ClusterQuiver(QuiverMutationType([['A',2],['B',2]])).digraph()
            Digraph on 4 vertices

            sage: ClusterQuiver(QuiverMutationType([['A',2],['B',2]])).digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 3, (1, -2))]

            sage: ClusterQuiver(['C', 4], user_labels = ['x', 'y', 'z', 'w']).digraph().edges(sort=True)
            [('x', 'y', (1, -1)), ('z', 'w', (2, -1)), ('z', 'y', (1, -1))]
        """
        return copy(self._digraph)

    def mutation_type(self):
        """
        Return the mutation type of ``self``.

        Return the mutation_type of each connected component of ``self`` if it
        can be determined, otherwise, the mutation type of this component is
        set to be unknown.

        The mutation types of the components are ordered by vertex labels.

        If you do many type recognitions, you should consider to save
        exceptional mutation types using
        :meth:`~sage.combinat.cluster_algebra_quiver.quiver_mutation_type.save_quiver_data`

        WARNING:

        - All finite types can be detected,
        - All affine types can be detected, EXCEPT affine type D (the algorithm is not yet implemented)
        - All exceptional types can be detected.

        EXAMPLES::

            sage: ClusterQuiver(['A',4]).mutation_type()
            ['A', 4]
            sage: ClusterQuiver(['A',(3,1),1]).mutation_type()
            ['A', [1, 3], 1]
            sage: ClusterQuiver(['C',2]).mutation_type()
            ['B', 2]
            sage: ClusterQuiver(['B',4,1]).mutation_type()
            ['BD', 4, 1]

        finite types::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q._mutation_type = None
            sage: Q.mutation_type()
            ['A', 5]

            sage: Q = ClusterQuiver([(0,1),(1,2),(2,3),(3,4)])
            sage: Q.mutation_type()
            ['A', 5]

            sage: Q = ClusterQuiver(DiGraph([['a', 'b'], ['c', 'b'], ['c', 'd'], ['e', 'd']]),
            ....:                   frozen=['c'])
            sage: Q.mutation_type()
            [ ['A', 2], ['A', 2] ]

        affine types::

            sage: Q = ClusterQuiver(['E',8,[1,1]]); Q
            Quiver on 10 vertices of type ['E', 8, [1, 1]]
            sage: Q._mutation_type = None; Q
            Quiver on 10 vertices
            sage: Q.mutation_type() # long time
            ['E', 8, [1, 1]]

        the not yet working affine type D (unless user has saved small classical quiver data)::

            sage: Q = ClusterQuiver(['D',4,1])
            sage: Q._mutation_type = None
            sage: Q.mutation_type() # todo: not implemented
            ['D', 4, 1]

        the exceptional types::

            sage: Q = ClusterQuiver(['X',6])
            sage: Q._mutation_type = None
            sage: Q.mutation_type() # long time
            ['X', 6]

        examples from page 8 of [Ke2008]_::

            sage: dg = DiGraph(); dg.add_edges([(9,0),(9,4),(4,6),(6,7),(7,8),(8,3),(3,5),(5,6),(8,1),(2,3)])
            sage: ClusterQuiver( dg ).mutation_type() # long time
            ['E', 8, [1, 1]]

            sage: dg = DiGraph( { 0:[3], 1:[0,4], 2:[0,6], 3:[1,2,7], 4:[3,8], 5:[2], 6:[3,5], 7:[4,6], 8:[7] } )
            sage: ClusterQuiver( dg ).mutation_type() # long time
            ['E', 8, 1]

            sage: dg = DiGraph( { 0:[3,9], 1:[0,4], 2:[0,6], 3:[1,2,7], 4:[3,8], 5:[2], 6:[3,5], 7:[4,6], 8:[7], 9:[1] } )
            sage: ClusterQuiver( dg ).mutation_type() # long time
            ['E', 8, [1, 1]]

        infinite types::

            sage: Q = ClusterQuiver(['GR',[4,9]])
            sage: Q._mutation_type = None
            sage: Q.mutation_type()
            'undetermined infinite mutation type'

        reducible types::

            sage: Q = ClusterQuiver([['A', 3], ['B', 3]])
            sage: Q._mutation_type = None
            sage: Q.mutation_type()
            [ ['A', 3], ['B', 3] ]

            sage: Q = ClusterQuiver([['A', 3], ['T', [4,4,4]]])
            sage: Q._mutation_type = None
            sage: Q.mutation_type()
            [['A', 3], 'undetermined infinite mutation type']

            sage: Q = ClusterQuiver([['A', 3], ['B', 3], ['T', [4,4,4]]])
            sage: Q._mutation_type = None
            sage: Q.mutation_type()
            [['A', 3], ['B', 3], 'undetermined infinite mutation type']

            sage: Q = ClusterQuiver([[0,1,2],[1,2,2],[2,0,2],[3,4,1],[4,5,1]])
            sage: Q.mutation_type()
            ['undetermined finite mutation type', ['A', 3]]

        TESTS::

            sage: Q = ClusterQuiver(matrix([[0, 3], [-1, 0], [1, 0], [0, 1]]))
            sage: Q.mutation_type()
            ['G', 2]
            sage: Q = ClusterQuiver(matrix([[0, -1, -1, 1, 0], [1, 0, 1, 0, 1], [1, -1, 0, -1, 0], [-1, 0, 1, 0, 1], [0, -1, 0, -1, 0], [0, 1, 0, -1, -1], [0, 1, -1, 0, 0]]))
            sage: Q.mutation_type()
            'undetermined infinite mutation type'
        """
        # checking if the mutation type is known already
        if self._mutation_type is None:
            # checking mutation type only for the principal part
            if self._m > 0:
                dg = self._digraph.subgraph(self._nlist)
            else:
                dg = self._digraph

            # checking the type for each connected component
            mutation_type = []
            connected_components = sorted(dg.connected_components(sort=False))
            for component in connected_components:
                # constructing the digraph for this component
                dg_component = dg.subgraph( component )
                dg_component.relabel()
                # turning dg_component into a canonical form
                _dg_canonical_form(dg_component)
                # turning dg_component into a canonical form
                dig6 = _digraph_to_dig6(dg_component, hashable=True)
                # and getting the corresponding matrix
                M = _dig6_to_matrix(dig6)

                # checking if this quiver is mutation infinite
                is_finite, path = is_mutation_finite(M)
                if is_finite is False:
                    mut_type_part = 'undetermined infinite mutation type'
                else:
                    # checking if this quiver is in the database
                    mut_type_part = _mutation_type_from_data( dg_component.order(), dig6, compute_if_necessary=False )
                    # checking if the algorithm can determine the mutation type
                    if mut_type_part == 'unknown':
                        mut_type_part = _connected_mutation_type(dg_component)
                    # checking if this quiver is of exceptional type by computing the exceptional mutation classes
                    if mut_type_part == 'unknown':
                        mut_type_part = _mutation_type_from_data(dg_component.order(), dig6, compute_if_necessary=True)
                    if mut_type_part == 'unknown':
                        mut_type_part = 'undetermined finite mutation type'
                mutation_type.append( mut_type_part )

            # the empty quiver case
            if len(mutation_type) == 0:
                mutation_type = None
            # the connected quiver case
            elif len(mutation_type) == 1:
                mutation_type = mutation_type[0]
            # the reducible quiver case
            elif len(mutation_type) > 1:
                if any(isinstance(mut_type_part, str)
                       for mut_type_part in mutation_type):
                    pass
                else:
                    mutation_type = QuiverMutationType(mutation_type)
            self._mutation_type = mutation_type
        return self._mutation_type

    def n(self):
        """
        Return the number of free vertices of ``self``.

        EXAMPLES::

            sage: ClusterQuiver(['A',4]).n()
            4
            sage: ClusterQuiver(['A',(3,1),1]).n()
            4
            sage: ClusterQuiver(['B',4]).n()
            4
            sage: ClusterQuiver(['B',4,1]).n()
            5
        """
        return self._n

    def m(self):
        """
        Return the number of frozen vertices of ``self``.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',4])
            sage: Q.m()
            0

            sage: T = ClusterQuiver(Q.digraph().edges(sort=True), frozen=[3])
            sage: T.n()
            3
            sage: T.m()
            1
        """
        return self._m

    def free_vertices(self):
        """
        Return the list of free vertices of ``self``.

        EXAMPLES::

            sage: Q = ClusterQuiver(DiGraph([['a', 'b'], ['c', 'b'], ['c', 'd'], ['e', 'd']]),
            ....:                   frozen=['b', 'd'])
            sage: Q.free_vertices()
            ['a', 'c', 'e']
        """
        return self._nlist

    def frozen_vertices(self):
        """
        Return the list of frozen vertices of ``self``.

        EXAMPLES::

            sage: Q = ClusterQuiver(DiGraph([['a', 'b'], ['c', 'b'], ['c', 'd'], ['e', 'd']]),
            ....:                   frozen=['b', 'd'])
            sage: sorted(Q.frozen_vertices())
            ['b', 'd']
        """
        return self._mlist

    def canonical_label(self, certificate=False):
        """
        Return the canonical labelling of ``self``.

        See :meth:`sage.graphs.generic_graph.GenericGraph.canonical_label`.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); if ``True``, the
          dictionary from ``self.vertices()`` to the vertices of the returned
          quiver is returned as well

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',4]); Q.digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]

            sage: T = Q.canonical_label(); T.digraph().edges(sort=True)
            [(0, 3, (1, -1)), (1, 2, (1, -1)), (1, 3, (1, -1))]

            sage: T, iso = Q.canonical_label(certificate=True)
            sage: T.digraph().edges(sort=True); iso
            [(0, 3, (1, -1)), (1, 2, (1, -1)), (1, 3, (1, -1))]
            {0: 0, 1: 3, 2: 1, 3: 2}

            sage: Q = ClusterQuiver(QuiverMutationType([['B',2],['A',1]])); Q
            Quiver on 3 vertices of type [ ['B', 2], ['A', 1] ]

            sage: Q.canonical_label()
            Quiver on 3 vertices of type [ ['A', 1], ['B', 2] ]

            sage: Q.canonical_label(certificate=True)
            (Quiver on 3 vertices of type [ ['A', 1], ['B', 2] ], {0: 1, 1: 2, 2: 0})
        """
        # computing the canonical form respecting the frozen variables
        dg = copy(self._digraph)
        iso, _ = _dg_canonical_form(dg, frozen=self._mlist)
        frozen = [iso[i] for i in self._mlist]
        Q = ClusterQuiver(dg, frozen=frozen)
        # getting the new ordering for the mutation type if necessary
        if self._mutation_type:
            if dg.is_connected():
                Q._mutation_type = self._mutation_type
            else:
                CC = sorted(self._digraph.connected_components(sort=False))
                CC_new = sorted(zip([sorted(iso[i] for i in L) for L in CC],
                                    range(len(CC))))
                comp_iso = [L[1] for L in CC_new]
                Q._mutation_type = []
                for i in range( len( CC_new ) ):
                    Q._mutation_type.append( copy( self._mutation_type.irreducible_components()[ comp_iso[i] ] ) )
                Q._mutation_type = QuiverMutationType( Q._mutation_type )
        if certificate:
            return Q, iso
        else:
            return Q

    def is_acyclic(self):
        """
        Return true if ``self`` is acyclic.

        EXAMPLES::

            sage: ClusterQuiver(['A',4]).is_acyclic()
            True

            sage: ClusterQuiver(['A',[2,1],1]).is_acyclic()
            True

            sage: ClusterQuiver([[0,1],[1,2],[2,0]]).is_acyclic()
            False
        """
        return self._digraph.is_directed_acyclic()

    def is_bipartite(self, return_bipartition=False):
        """
        Return ``True`` if ``self`` is bipartite.

        EXAMPLES::

            sage: ClusterQuiver(['A',[3,3],1]).is_bipartite()
            True

            sage: ClusterQuiver(['A',[4,3],1]).is_bipartite()
            False
        """
        dg = copy(self._digraph)
        dg.delete_vertices(list(range(self._n, self._n + self._m)))
        if any(dg.in_degree(i) and dg.out_degree(i) for i in dg):
            return False
        if not return_bipartition:
            return True
        return dg.to_undirected().bipartite_sets()

    def exchangeable_part(self):
        """
        Return the restriction to the principal part (i.e. exchangeable part) of ``self``, the subquiver obtained by deleting the frozen vertices of ``self``.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',4])
            sage: T = ClusterQuiver(Q.digraph().edges(sort=True), frozen=[3])
            sage: T.digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]

            sage: T.exchangeable_part().digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 1, (1, -1))]

            sage: Q2 = Q.principal_extension()
            sage: Q3 = Q2.principal_extension()
            sage: Q2.exchangeable_part() == Q3.exchangeable_part()
            True
        """
        dg = DiGraph( self._digraph )
        dg.delete_vertices(list(range(self._n, self._n + self._m)))
        Q = ClusterQuiver( dg )
        Q._mutation_type = self._mutation_type
        return Q

    def principal_extension(self, inplace=False):
        """
        Return the principal extension of ``self``, adding n frozen vertices to any previously frozen vertices. I.e., the quiver obtained by adding an outgoing edge to every mutable vertex of ``self``.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',2]); Q
            Quiver on 2 vertices of type ['A', 2]
            sage: T = Q.principal_extension(); T
            Quiver on 4 vertices of type ['A', 2] with 2 frozen vertices
            sage: T2 = T.principal_extension(); T2
            Quiver on 6 vertices of type ['A', 2] with 4 frozen vertices
            sage: Q.digraph().edges(sort=True)
            [(0, 1, (1, -1))]
            sage: T.digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 0, (1, -1)), (3, 1, (1, -1))]
            sage: T2.digraph().edges(sort=True)
            [(0, 1, (1, -1)), (2, 0, (1, -1)), (3, 1, (1, -1)), (4, 0, (1, -1)), (5, 1, (1, -1))]
        """
        dg = DiGraph(self._digraph)
        dg.add_edges( [(self._n+self._m+i,i) for i in range(self._n)] )
        Q = ClusterQuiver( dg, frozen=list(range(self._n, self._n+self._m+self._n)) )
        Q._mutation_type = self._mutation_type
        if inplace:
            self.__init__(Q)
        else:
            return Q

    def first_sink(self):
        r"""
        Return the first vertex of ``self`` that is a sink.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.mutate([1,2,4,3,2])
            sage: Q.first_sink()
            0
        """
        sinks = self.digraph().sinks()

        if sinks:
            return sinks[0]
        return None

    def sinks(self):
        r"""
        Return all vertices of ``self`` that are sinks.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.mutate([1,2,4,3,2])
            sage: Q.sinks()
            [0, 2]

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.mutate([2,1,3,4,2])
            sage: Q.sinks()
            [3]
        """
        return self.digraph().sinks()

    def first_source(self):
        r"""
        Return the first vertex of ``self`` that is a source.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.mutate([2,1,3,4,2])
            sage: Q.first_source()
            1
        """
        sources = self.digraph().sources()

        if sources:
            return sources[0]
        return None

    def sources(self):
        r"""
        Return all vertices of ``self`` that are sources.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.mutate([1,2,4,3,2])
            sage: Q.sources()
            []

            sage: Q = ClusterQuiver(['A',5])
            sage: Q.mutate([2,1,3,4,2])
            sage: Q.sources()
            [1]
        """
        return self.digraph().sources()

    def mutate(self, data, inplace=True):
        """
        Mutates ``self`` at a sequence of vertices.

        INPUT:

        - ``sequence`` -- a vertex of ``self``, an iterator of vertices of
          ``self``, a function which takes in the ClusterQuiver and returns a
          vertex or an iterator of vertices, or a string of the parameter
          wanting to be called on ClusterQuiver that will return a vertex or
          an iterator of vertices
        - ``inplace`` -- boolean (default: ``True``); if ``False``, the result
          is returned, otherwise ``self`` is modified

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',4]); Q.b_matrix()
            [ 0  1  0  0]
            [-1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -1  0]

            sage: Q.mutate(0); Q.b_matrix()
            [ 0 -1  0  0]
            [ 1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -1  0]

            sage: T = Q.mutate(0, inplace=False); T
            Quiver on 4 vertices of type ['A', 4]

            sage: Q.mutate(0)
            sage: Q == T
            True

            sage: Q.mutate([0,1,0])
            sage: Q.b_matrix()
            [ 0 -1  1  0]
            [ 1  0  0  0]
            [-1  0  0  1]
            [ 0  0 -1  0]

            sage: Q = ClusterQuiver(QuiverMutationType([['A',1],['A',3]]))
            sage: Q.b_matrix()
            [ 0  0  0  0]
            [ 0  0  1  0]
            [ 0 -1  0 -1]
            [ 0  0  1  0]

            sage: T = Q.mutate(0,inplace=False)
            sage: Q == T
            True

            sage: Q = ClusterQuiver(['A',3]); Q.b_matrix()
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
            sage: Q.mutate('first_sink'); Q.b_matrix()
            [ 0 -1  0]
            [ 1  0  1]
            [ 0 -1  0]
            sage: Q.mutate('first_source'); Q.b_matrix()
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]

            sage: dg = DiGraph()
            sage: dg.add_vertices(['a','b','c','d','e'])
            sage: dg.add_edges([['a','b'], ['b','c'], ['c','d'], ['d','e']])
            sage: Q2 = ClusterQuiver(dg, frozen=['c']); Q2.b_matrix()
            [ 0  1  0  0]
            [-1  0  0  0]
            [ 0  0  0  1]
            [ 0  0 -1  0]
            [ 0 -1  1  0]
            sage: Q2.mutate('a'); Q2.b_matrix()
            [ 0 -1  0  0]
            [ 1  0  0  0]
            [ 0  0  0  1]
            [ 0  0 -1  0]
            [ 0 -1  1  0]

            sage: dg = DiGraph([['a', 'b'], ['b', 'c']], format='list_of_edges')
            sage: Q = ClusterQuiver(dg);Q
            Quiver on 3 vertices
            sage: Q.mutate(['a','b'],inplace=False).digraph().edges(sort=True)
            [('a', 'b', (1, -1)), ('c', 'b', (1, -1))]

        TESTS::

            sage: Q = ClusterQuiver(['A',4]); Q.mutate(0,1)
            Traceback (most recent call last):
            ...
            ValueError: The second parameter must be boolean.  To mutate at a sequence of length 2, input it as a list.

            sage: Q = ClusterQuiver(['A',4]); Q.mutate(0,0)
            Traceback (most recent call last):
            ...
            ValueError: The second parameter must be boolean.  To mutate at a sequence of length 2, input it as a list.
        """
        dg = self._digraph
        V = nlist = self._nlist
        mlist = self._mlist

        # If we get a string which is not a cluster variable, execute as a function
        if isinstance(data, str):
            if data not in V:
                data = getattr(self, data)()

        # If we get a function, execute it
        if callable(data):
            # function should return either integer or sequence
            data = data(self)

        if data is None:
            raise ValueError('Not mutating: No vertices given.')

        if data in V:
            seq = [data]
        else:
            seq = data
        if isinstance(seq, tuple):
            seq = list( seq )
        if not isinstance(seq, list):
            raise ValueError('The quiver can only be mutated at a vertex or at a sequence of vertices')
        if not isinstance(inplace, bool):
            raise ValueError('The second parameter must be boolean.  To mutate at a sequence of length 2, input it as a list.')
        if any(v not in V for v in seq):
            v = next(v for v in seq if v not in V)
            raise ValueError('The quiver cannot be mutated at the vertex %s' % v)

        for v in seq:
            dg = _digraph_mutate(dg, v, frozen=mlist)

        if inplace:
            self._M = _edge_list_to_matrix(dg.edge_iterator(), nlist, mlist)
            self._M.set_immutable()
            self._digraph = dg
        else:
            Q = ClusterQuiver(dg, frozen=mlist)
            Q._mutation_type = self._mutation_type
            return Q

    def mutation_sequence(self, sequence, show_sequence=False, fig_size=1.2 ):
        """
        Return a list containing the sequence of quivers obtained from ``self``
        by a sequence of mutations on vertices.

        INPUT:

        - ``sequence`` -- list or tuple of vertices of ``self``
        - ``show_sequence`` -- boolean (default: ``False``); if ``True``, a png
          containing the mutation sequence is shown
        - ``fig_size`` -- (default: 1.2) factor by which the size of the
          sequence is expanded

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',4])
            sage: seq = Q.mutation_sequence([0,1]); seq
            [Quiver on 4 vertices of type ['A', 4],
             Quiver on 4 vertices of type ['A', 4],
             Quiver on 4 vertices of type ['A', 4]]
            sage: [T.b_matrix() for T in seq]
            [
            [ 0  1  0  0]  [ 0 -1  0  0]  [ 0  1 -1  0]
            [-1  0 -1  0]  [ 1  0 -1  0]  [-1  0  1  0]
            [ 0  1  0  1]  [ 0  1  0  1]  [ 1 -1  0  1]
            [ 0  0 -1  0], [ 0  0 -1  0], [ 0  0 -1  0]
            ]
        """
        n = self._n
        m = self._m
        if m == 0:
            width_factor = 3
            fig_size = fig_size*2*n/3
        else:
            width_factor = 6
            fig_size = fig_size*4*n/3
        V = list(range(n))

        if isinstance(sequence, tuple):
            sequence = list( sequence )
        if not isinstance(sequence, list):
            raise ValueError('The quiver can only be mutated at a vertex or at a sequence of vertices')
        if any(v not in V for v in sequence):
            v = next(v for v in sequence if v not in V)
            raise ValueError('The quiver can only be mutated at the vertex %s' % v )

        quiver = copy( self )
        quiver_sequence = []
        quiver_sequence.append( copy( quiver ) )

        for v in sequence:
            quiver.mutate( v )
            quiver_sequence.append( copy( quiver ) )

        if show_sequence:
            from sage.plot.plot import Graphics
            from sage.plot.text import text

            def _plot_arrow( v, k, center=(0,0) ):
                return text(r"$\longleftrightarrow$",(center[0],center[1]), fontsize=25) + text(r"$\mu_"+str(v)+"$",(center[0],center[1]+0.15), fontsize=15) \
                    + text("$"+str(k)+"$",(center[0],center[1]-0.2), fontsize=15)

            plot_sequence = [ quiver_sequence[i].plot( circular=True, center=(i*width_factor,0) ) for i in range(len(quiver_sequence)) ]
            arrow_sequence = [ _plot_arrow( sequence[i],i+1,center=((i+0.5)*width_factor,0) ) for i in range(len(sequence)) ]
            sequence = []
            for i in range( len( plot_sequence ) ):
                if i < len( arrow_sequence ):
                    sequence.append( plot_sequence[i] + arrow_sequence[i] )
                else:
                    sequence.append( plot_sequence[i] )
            plot_obj = Graphics()
            for elem in sequence:
                plot_obj += elem
            plot_obj.show(axes=False, figsize=[fig_size*len(quiver_sequence),fig_size])
        return quiver_sequence

    def reorient(self, data):
        """
        Reorient ``self`` with respect to the given total order, or
        with respect to an iterator of edges in ``self`` to be
        reverted.

        .. WARNING::

            This operation might change the mutation type of ``self``.

        INPUT:

        - ``data`` -- an iterator defining a total order on
          ``self.vertices()``, or an iterator of edges in ``self`` to
          be reoriented.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',(2,3),1])
            sage: Q.mutation_type()
            ['A', [2, 3], 1]

            sage: Q.reorient([(0,1),(1,2),(2,3),(3,4)])
            sage: Q.mutation_type()
            ['D', 5]

            sage: Q.reorient([0,1,2,3,4])
            sage: Q.mutation_type()
            ['A', [1, 4], 1]

        TESTS::

            sage: Q = ClusterQuiver(['A',2])
            sage: Q.reorient([])
            Traceback (most recent call last):
            ...
            ValueError: empty input
            sage: Q.reorient([3,4])
            Traceback (most recent call last):
            ...
            ValueError: not a total order on the vertices of the quiver or
            a list of edges to be oriented
        """
        if not data:
            raise ValueError('empty input')
        first = data[0]

        if set(data) == set(range(self._n + self._m)):
            dg_new = DiGraph()
            for edge in self._digraph.edges(sort=True):
                if data.index(edge[0]) < data.index(edge[1]):
                    dg_new.add_edge(edge[0], edge[1], edge[2])
                else:
                    dg_new.add_edge(edge[1], edge[0], edge[2])
            self._digraph = dg_new
            self._M = _edge_list_to_matrix(dg_new.edges(sort=True),
                                           self._nlist, self._mlist)
            self._M.set_immutable()
            self._mutation_type = None
        elif isinstance(first, (list, tuple)) and len(first) == 2:
            edges = self._digraph.edges(sort=True, labels=False)
            for edge in data:
                if (edge[1], edge[0]) in edges:
                    label = self._digraph.edge_label(edge[1], edge[0])
                    self._digraph.delete_edge(edge[1], edge[0])
                    self._digraph.add_edge(edge[0], edge[1], label)
            self._M = _edge_list_to_matrix(self._digraph.edges(sort=True),
                                           self._nlist, self._mlist)
            self._M.set_immutable()
            self._mutation_type = None
        else:
            raise ValueError('not a total order on the vertices of the quiver or a list of edges to be oriented')

    def mutation_class_iter(self, depth=infinity, show_depth=False,
                            return_paths=False, data_type='quiver',
                            up_to_equivalence=True, sink_source=False):
        """
        Return an iterator for the mutation class of ``self`` together with certain constraints.

        INPUT:

        - ``depth`` -- integer (default: infinity); only quivers with distance
          at most depth from ``self`` are returned.
        - ``show_depth`` -- boolean (default: ``False``); if ``True``, the
          actual depth of the mutation is shown
        - ``return_paths`` -- boolean (default: ``False``); if ``True``, a
          shortest path of mutation sequences from ``self`` to the given quiver
          is returned as well
        - ``data_type`` -- (default: ``'quiver'``) can be one of the following::

            * "quiver"
            * "matrix"
            * "digraph"
            * "dig6"
            * "path"

        - ``up_to_equivalence`` -- boolean (default: ``True``); if ``True``,
          only one quiver for each graph-isomorphism class is recorded
        - ``sink_source`` -- boolean (default: ``False``); if ``True``, only
          mutations at sinks and sources are applied

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',3])
            sage: it = Q.mutation_class_iter()
            sage: for T in it: print(T)
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]

            sage: it = Q.mutation_class_iter(depth=1)
            sage: for T in it: print(T)
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]

            sage: it = Q.mutation_class_iter(show_depth=True)
            sage: for T in it: pass
            Depth: 0     found: 1          Time: ... s
            Depth: 1     found: 3          Time: ... s
            Depth: 2     found: 4          Time: ... s

            sage: it = Q.mutation_class_iter(return_paths=True)
            sage: for T in it: print(T)
            (Quiver on 3 vertices of type ['A', 3], [])
            (Quiver on 3 vertices of type ['A', 3], [1])
            (Quiver on 3 vertices of type ['A', 3], [0])
            (Quiver on 3 vertices of type ['A', 3], [0, 1])

            sage: it = Q.mutation_class_iter(up_to_equivalence=False)
            sage: for T in it: print(T)
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]

            sage: it = Q.mutation_class_iter(return_paths=True,up_to_equivalence=False)
            sage: mutation_class = list(it)
            sage: len(mutation_class)
            14
            sage: mutation_class[0]
            (Quiver on 3 vertices of type ['A', 3], [])

            sage: Q = ClusterQuiver(['A',3])
            sage: it = Q.mutation_class_iter(data_type='path')
            sage: for T in it: print(T)
            []
            [1]
            [0]
            [0, 1]

            sage: Q = ClusterQuiver(['A',3])
            sage: it = Q.mutation_class_iter(return_paths=True,data_type='matrix')
            sage: next(it)
            (
            [ 0  0  1]
            [ 0  0  1]
            [-1 -1  0], []
            )

            sage: dg = DiGraph([['a', 'b'], ['b', 'c']], format='list_of_edges')
            sage: S = ClusterQuiver(dg, frozen=['b'])
            sage: S.mutation_class()
            [Quiver on 3 vertices with 1 frozen vertex,
             Quiver on 3 vertices with 1 frozen vertex,
             Quiver on 3 vertices with 1 frozen vertex]
        """
        if data_type == 'path':
            return_paths = False
        if data_type == "dig6":
            return_dig6 = True
        else:
            return_dig6 = False

        # jump to the standard labelling convention
        dg = ClusterQuiver(self._M).digraph()
        frozen = list(range(self._n, self._n + self._m))

        MC_iter = _mutation_class_iter(dg, self._n, self._m, depth=depth,
                                       return_dig6=return_dig6,
                                       show_depth=show_depth,
                                       up_to_equivalence=up_to_equivalence,
                                       sink_source=sink_source)
        for data in MC_iter:
            if data_type == "quiver":
                next_element = ClusterQuiver(data[0], frozen=frozen)
                next_element._mutation_type = self._mutation_type
            elif data_type == "matrix":
                next_element = ClusterQuiver(data[0], frozen=frozen)._M
            elif data_type == "digraph":
                next_element = data[0]
            elif data_type == "dig6":
                next_element = data[0]
            elif data_type == "path":
                next_element = data[1]
            else:
                raise ValueError("the parameter for data_type was "
                                 "not recognized")
            if return_paths:
                yield (next_element, data[1])
            else:
                yield next_element

    def mutation_class(self, depth=infinity, show_depth=False, return_paths=False,
                       data_type='quiver', up_to_equivalence=True, sink_source=False):
        """
        Return the mutation class of ``self`` together with certain constraints.

        INPUT:

        - ``depth`` -- (default: ``infinity``) integer, only seeds with
          distance at most depth from ``self`` are returned
        - ``show_depth`` -- boolean (default: ``False``); if ``True``, the
          actual depth of the mutation is shown
        - ``return_paths`` -- boolean (default: ``False``); if ``True``, a
          shortest path of mutation sequences from ``self`` to the given quiver is
          returned as well
        - ``data_type`` -- (default: ``'quiver'``) can be one of
          the following:

          * ``'quiver'`` -- the quiver is returned
          * ``'dig6'`` -- the dig6-data is returned
          * ``'path'`` -- shortest paths of mutation sequences from
            ``self`` are returned

        - ``sink_source`` -- boolean (default: ``False``); if ``True``, only
          mutations at sinks and sources are applied

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',3])
            sage: Ts = Q.mutation_class()
            sage: for T in Ts: print(T)
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]

            sage: Ts = Q.mutation_class(depth=1)
            sage: for T in Ts: print(T)
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]

            sage: Ts = Q.mutation_class(show_depth=True)
            Depth: 0     found: 1          Time: ... s
            Depth: 1     found: 3          Time: ... s
            Depth: 2     found: 4          Time: ... s

            sage: Ts = Q.mutation_class(return_paths=True)
            sage: for T in Ts: print(T)
            (Quiver on 3 vertices of type ['A', 3], [])
            (Quiver on 3 vertices of type ['A', 3], [1])
            (Quiver on 3 vertices of type ['A', 3], [0])
            (Quiver on 3 vertices of type ['A', 3], [0, 1])

            sage: Ts = Q.mutation_class(up_to_equivalence=False)
            sage: for T in Ts: print(T)
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]
            Quiver on 3 vertices of type ['A', 3]

            sage: Ts = Q.mutation_class(return_paths=True,up_to_equivalence=False)
            sage: len(Ts)
            14
            sage: Ts[0]
            (Quiver on 3 vertices of type ['A', 3], [])

            sage: Ts = Q.mutation_class(show_depth=True)
            Depth: 0     found: 1          Time: ... s
            Depth: 1     found: 3          Time: ... s
            Depth: 2     found: 4          Time: ... s

            sage: Ts = Q.mutation_class(show_depth=True, up_to_equivalence=False)
            Depth: 0     found: 1          Time: ... s
            Depth: 1     found: 4          Time: ... s
            Depth: 2     found: 6          Time: ... s
            Depth: 3     found: 10        Time: ... s
            Depth: 4     found: 14        Time: ... s

        TESTS::

            sage: all(len(ClusterQuiver(['A',n]).mutation_class())
            ....:     == ClusterQuiver(['A',n]).mutation_type().class_size()
            ....:     for n in [2..6])
            True

            sage: all(len(ClusterQuiver(['B',n]).mutation_class())
            ....:     == ClusterQuiver(['B',n]).mutation_type().class_size()
            ....:     for n in [2..6])
            True
        """
        if depth is infinity and not self.is_mutation_finite():
            raise ValueError('the mutation class can - for infinite mutation'
                             ' types - only be computed up to a given depth')
        return list(self.mutation_class_iter(depth=depth, show_depth=show_depth,
                                             return_paths=return_paths,
                                             data_type=data_type,
                                             up_to_equivalence=up_to_equivalence,
                                             sink_source=sink_source))

    def is_finite(self):
        """
        Return ``True`` if ``self`` is of finite type.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',3])
            sage: Q.is_finite()
            True
            sage: Q = ClusterQuiver(['A',[2,2],1])
            sage: Q.is_finite()
            False
            sage: Q = ClusterQuiver([['A',3],['B',3]])
            sage: Q.is_finite()
            True
            sage: Q = ClusterQuiver(['T',[4,4,4]])
            sage: Q.is_finite()
            False
            sage: Q = ClusterQuiver([['A',3],['T',[4,4,4]]])
            sage: Q.is_finite()
            False
            sage: Q = ClusterQuiver([['A',3],['T',[2,2,3]]])
            sage: Q.is_finite()
            True
            sage: Q = ClusterQuiver([['A',3],['D',5]])
            sage: Q.is_finite()
            True
            sage: Q = ClusterQuiver([['A',3],['D',5,1]])
            sage: Q.is_finite()
            False

            sage: Q = ClusterQuiver([[0,1,2],[1,2,2],[2,0,2]])
            sage: Q.is_finite()
            False

            sage: Q = ClusterQuiver([[0,1,2],[1,2,2],[2,0,2],[3,4,1],[4,5,1]])
            sage: Q.is_finite()
            False
        """
        mt = self.mutation_type()
        return (type(mt) in [QuiverMutationType_Irreducible,
                             QuiverMutationType_Reducible] and mt.is_finite())

    def is_mutation_finite( self, nr_of_checks=None, return_path=False ):
        """
        Uses a non-deterministic method by random mutations in various
        directions. Can result in a wrong answer.

        INPUT:

        - ``nr_of_checks`` -- (default: ``None``) number of mutations applied;
          Standard is 500*(number of vertices of ``self``)
        - ``return_path`` -- (default: ``False``) if ``True``, in case of
          ``self`` not being mutation finite, a path from ``self`` to a quiver
          with an edge label `(a,-b)` and `a*b > 4` is returned

        ALGORITHM:

        A quiver is mutation infinite if and only if every edge label (a,-b) satisfy a*b > 4.
        Thus, we apply random mutations in random directions

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',10])
            sage: Q._mutation_type = None
            sage: Q.is_mutation_finite()
            True

            sage: Q = ClusterQuiver([(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(2,9)])
            sage: Q.is_mutation_finite()
            False
        """
        if self._n <= 2:
            is_finite = True
            path = None
        elif not return_path and self._mutation_type == 'undetermined infinite mutation type':
            is_finite = False
        elif (isinstance(self._mutation_type, (QuiverMutationType_Irreducible,
                                               QuiverMutationType_Reducible))
              and self._mutation_type.is_mutation_finite()):
            is_finite = True
            path = None
        elif (not return_path and isinstance(self._mutation_type,
                                             (QuiverMutationType_Irreducible,
                                              QuiverMutationType_Reducible))
              and not self._mutation_type.is_mutation_finite()):
            is_finite = False
        else:
            # turning dg_component into a canonical form
            dig6 = _digraph_to_dig6(self.digraph())
            # and getting the corresponding matrix
            M = _dig6_to_matrix(dig6)

            is_finite, path = is_mutation_finite(M, nr_of_checks=nr_of_checks)
        if return_path:
            return is_finite, path
        else:
            return is_finite

    def number_of_edges(self):
        r"""
        Return the total number of edges on the quiver.

        Note: This only works with non-valued quivers. If used on a
        non-valued quiver then the positive value is taken to be the number of edges added

        OUTPUT: integer of the number of edges

        EXAMPLES::

            sage: S = ClusterQuiver(['A',4]); S.number_of_edges()
            3

            sage: S = ClusterQuiver(['B',4]); S.number_of_edges()
            3
        """
        digraph_edges = self.digraph().edges(sort=True)

        total_edges = 0
        for edge in digraph_edges:
            total_edges += edge[2][0]

        return total_edges

    def relabel(self, relabelling, inplace=True):
        r"""
        Return the quiver after doing a relabelling.

        Will relabel the vertices of the quiver.

        INPUT:

        - ``relabelling`` -- dictionary of labels to move around
        - ``inplace`` -- boolean (default: ``True``); if ``True``, will return
          a duplicate of the quiver

        EXAMPLES::

            sage: S = ClusterQuiver(['A',4]).relabel({1:'5',2:'go'})
        """
        if inplace:
            quiver = self
        else:
            quiver = ClusterQuiver(self)

        # Instantiation variables
        old_vertices = list(quiver.digraph())
        digraph_labels = {}
        dict_labels = {}

        # Organize labels noting that for:
        #    _digraph: { old_vertex: new_vertex}
        #    _vertex_dictionary: {num: new_vertex}
        if isinstance(relabelling, list):
            digraph_labels = {old_vertices[i]: relabelling[i] for i in range(len(relabelling))}
            dict_labels = {range(len(relabelling))[i]: relabelling[i] for i in range(len(relabelling))}
        elif isinstance(relabelling, dict):
            # need to make sure we map correctly
            for key in relabelling:
                val = relabelling[key]

                if key in old_vertices:
                    # If the key is in the old vertices, use that mapping
                    digraph_labels[key] = val
                    # And place it in the right order for our dictionary
                    loc = [i for i, x in enumerate(old_vertices)
                           if x == key][0]
                    dict_labels[loc] = val
                elif isinstance(key, int) and len(old_vertices) > key:
                    # If the key is an integer, grab that particular vertex
                    digraph_labels[old_vertices[key]] = val
                    # And copy it over to our dictionary
                    dict_labels[key] = val

        quiver._digraph.relabel(digraph_labels)
        quiver._vertex_dictionary = dict_labels
        return quiver

    def poincare_semistable(self, theta, d):
        r"""
        Return the Poincaré polynomial of the moduli space of semi-stable
        representations of dimension vector `d`.

        INPUT:

        - ``theta`` -- stability weight, as list or vector of rationals
        - ``d`` -- dimension vector, as list or vector of coprime integers

        The semi-stability is taken with respect to the slope function

        .. MATH::

             \mu(d) = \theta(d) / \operatorname{dim}(d)

        where `d` is a dimension vector.

        This uses the matrix-inversion algorithm from [Rei2002]_.

        EXAMPLES::

            sage: Q = ClusterQuiver(['A',2])
            sage: Q.poincare_semistable([1,0],[1,0])
            1
            sage: Q.poincare_semistable([1,0],[1,1])
            1

            sage: K2 = ClusterQuiver(matrix([[0,2],[-2,0]]))
            sage: theta = (1, 0)
            sage: K2.poincare_semistable(theta, [1,0])
            1
            sage: K2.poincare_semistable(theta, [1,1])
            v^2 + 1
            sage: K2.poincare_semistable(theta, [1,2])
            1
            sage: K2.poincare_semistable(theta, [1,3])
            0

            sage: K3 = ClusterQuiver(matrix([[0,3],[-3,0]]))
            sage: theta = (1, 0)
            sage: K3.poincare_semistable(theta, (2,3))
            v^12 + v^10 + 3*v^8 + 3*v^6 + 3*v^4 + v^2 + 1
            sage: K3.poincare_semistable(theta, (3,4))(1)
            68

        TESTS::

            sage: Q = ClusterQuiver(['A',2])
            sage: Q.poincare_semistable([1,0],[2,2])
            Traceback (most recent call last):
            ...
            ValueError: dimension vector d is not coprime

            sage: Q = ClusterQuiver(['A',3])
            sage: Q.poincare_semistable([1,1,0],[2,3,4])
            0

        REFERENCES:

        .. [Rei2002] Markus Reineke, *The Harder-Narasimhan system in quantum
           groups and cohomology of quiver moduli*, :arxiv:`math/0204059`
        """
        if gcd([x for x in d if x]) != 1:
            raise ValueError("dimension vector d is not coprime")
        d = vector(ZZ, d)
        theta = vector(theta)

        n = self.n()
        b_mat = self.b_matrix()
        Eu = matrix(ZZ, n, n,
                    lambda i, j: -b_mat[i, j] if b_mat[i, j] > 0 else 0)
        Eu = 1 + Eu
        edges = list(self.digraph().edges(sort=True, labels=False))

        mu_d = theta.dot_product(d) / sum(d)

        Li = [0 * d]
        it = (vector(e) for e in product(*[range(d_i + 1)
                                           for d_i in d]))
        Li += [e for e in it if e.dot_product(theta) > mu_d * sum(e)]
        Li.append(d)
        N = len(Li) - 1

        q = polygen(QQ, 'v')  # q stands for v**2 until the last line

        def cardinal_RG(d):
            cardinal_G = prod(q**d_i - q**k for d_i in d for k in range(d_i))
            cardinal_R = prod(q**(b_mat[i, j] * d[i] * d[j])
                              for i, j in edges)
            return cardinal_R / cardinal_G

        Reineke_submat = matrix(q.parent().fraction_field(), N, N)

        for i, e in enumerate(Li[:-1]):
            for j, f in enumerate(Li[1:]):
                if e == f:
                    Reineke_submat[i, j] = 1
                    continue
                f_e = f - e
                if all(x >= 0 for x in f_e):
                    power = (-f_e) * Eu * e
                    Reineke_submat[i, j] = q**power * cardinal_RG(f_e)

        poly = (-1)**N * ((1 - q) * Reineke_submat.det()).numerator()
        return poly(q**2)  # replacing q by v**2

    def d_vector_fan(self):
        r"""
        Return the d-vector fan associated with the quiver.

        It is the fan whose maximal cones are generated by the
        d-matrices of the clusters.

        This is a complete simplicial fan (and even smooth when the
        initial quiver is acyclic). It only makes sense for quivers of
        finite type.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron sage.libs.singular
            sage: Fd = ClusterQuiver([[1,2]]).d_vector_fan(); Fd
            Rational polyhedral fan in 2-d lattice N
            sage: Fd.ngenerating_cones()
            5
            sage: Fd = ClusterQuiver([[1,2],[2,3]]).d_vector_fan(); Fd
            Rational polyhedral fan in 3-d lattice N
            sage: Fd.ngenerating_cones()
            14
            sage: Fd.is_smooth()
            True
            sage: Fd = ClusterQuiver([[1,2],[2,3],[3,1]]).d_vector_fan(); Fd
            Rational polyhedral fan in 3-d lattice N
            sage: Fd.ngenerating_cones()
            14
            sage: Fd.is_smooth()
            False

        TESTS::

            sage: ClusterQuiver(['A',[2,2],1]).d_vector_fan()
            Traceback (most recent call last):
            ...
            ValueError: only makes sense for quivers of finite type
        """
        if not self.is_finite():
            raise ValueError('only makes sense for quivers of finite type')

        from .cluster_seed import ClusterSeed
        from sage.geometry.fan import Fan
        from sage.geometry.cone import Cone

        seed = ClusterSeed(self)
        return Fan([Cone(s.d_matrix().columns())
                    for s in seed.mutation_class()])

    def g_vector_fan(self):
        r"""
        Return the g-vector fan associated with the quiver.

        It is the fan whose maximal cones are generated by the
        g-matrices of the clusters.

        This is a complete simplicial fan. It is only supported for
        quivers of finite type.

        EXAMPLES::

            sage: Fg = ClusterQuiver([[1,2]]).g_vector_fan(); Fg
            Rational polyhedral fan in 2-d lattice N
            sage: Fg.ngenerating_cones()
            5

            sage: Fg = ClusterQuiver([[1,2],[2,3]]).g_vector_fan(); Fg
            Rational polyhedral fan in 3-d lattice N
            sage: Fg.ngenerating_cones()
            14
            sage: Fg.is_smooth()
            True

            sage: Fg = ClusterQuiver([[1,2],[2,3],[3,1]]).g_vector_fan(); Fg
            Rational polyhedral fan in 3-d lattice N
            sage: Fg.ngenerating_cones()
            14
            sage: Fg.is_smooth()
            True

        TESTS::

            sage: ClusterQuiver(['A',[2,2],1]).g_vector_fan()
            Traceback (most recent call last):
            ...
            ValueError: only supported for quivers of finite type
        """
        from .cluster_seed import ClusterSeed
        from sage.geometry.fan import Fan
        from sage.geometry.cone import Cone

        if not (self.is_finite()):
            raise ValueError('only supported for quivers of finite type')
        seed = ClusterSeed(self).principal_extension()
        return Fan([Cone(s.g_matrix().columns())
                    for s in seed.mutation_class()])
