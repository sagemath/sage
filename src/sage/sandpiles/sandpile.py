r"""
Sandpiles

Functions and classes for mathematical sandpiles.

Version: 2.3

AUTHOR:

- Marshall Hampton (2010-1-10) modified for inclusion as a module within Sage
  library.

- David Perkinson (2010-12-14) added show3d(), fixed bug in resolution(),
  replaced elementary_divisors() with invariant_factors(), added show() for
  SandpileConfig and SandpileDivisor.

- David Perkinson (2010-9-18): removed is_undirected, added show(), added
  verbose arguments to several functions to display SandpileConfigs and
  divisors as lists of integers

- David Perkinson (2010-12-19): created separate SandpileConfig,
  SandpileDivisor, and Sandpile classes

- David Perkinson (2009-07-15): switched to using config_to_list instead of
  .values(), thus fixing a few bugs when not using integer labels for vertices.

- David Perkinson (2009): many undocumented improvements

- David Perkinson (2008-12-27): initial version

EXAMPLES:

A weighted directed graph given as a Python dictionary::

    sage: from sage.sandpiles import *
    sage: g = {0: {},                    \
               1: {0: 1, 2: 1, 3: 1},    \
               2: {1: 1, 3: 1, 4: 1},    \
               3: {1: 1, 2: 1, 4: 1},    \
               4: {2: 1, 3: 1}}

The associated sandpile with 0 chosen as the sink::

    sage: S = Sandpile(g,0)

A picture of the graph::

    sage: S.show()

The relevant Laplacian matrices::

    sage: S.laplacian()
    [ 0  0  0  0  0]
    [-1  3 -1 -1  0]
    [ 0 -1  3 -1 -1]
    [ 0 -1 -1  3 -1]
    [ 0  0 -1 -1  2]
    sage: S.reduced_laplacian()
    [ 3 -1 -1  0]
    [-1  3 -1 -1]
    [-1 -1  3 -1]
    [ 0 -1 -1  2]

The number of elements of the sandpile group for S::

    sage: S.group_order()
    8

The structure of the sandpile group::

    sage: S.invariant_factors()
    [1, 1, 1, 8]

The elements of the sandpile group for S::

    sage: S.recurrents()
    [{1: 2, 2: 2, 3: 2, 4: 1},
     {1: 2, 2: 2, 3: 2, 4: 0},
     {1: 2, 2: 1, 3: 2, 4: 0},
     {1: 2, 2: 2, 3: 0, 4: 1},
     {1: 2, 2: 0, 3: 2, 4: 1},
     {1: 2, 2: 2, 3: 1, 4: 0},
     {1: 2, 2: 1, 3: 2, 4: 1},
     {1: 2, 2: 2, 3: 1, 4: 1}]

The maximal stable element (2 grains of sand on vertices 1, 2, and 3, and 1
grain of sand on vertex 4::

    sage: S.max_stable()
    {1: 2, 2: 2, 3: 2, 4: 1}
    sage: S.max_stable().values()
    [2, 2, 2, 1]

The identity of the sandpile group for S::

    sage: S.identity()
    {1: 2, 2: 2, 3: 2, 4: 0}

Some group operations::

    sage: m = S.max_stable()
    sage: i = S.identity()
    sage: m.values()
    [2, 2, 2, 1]
    sage: i.values()
    [2, 2, 2, 0]
    sage: m+i    # coordinate-wise sum
    {1: 4, 2: 4, 3: 4, 4: 1}
    sage: m - i
    {1: 0, 2: 0, 3: 0, 4: 1}
    sage: m & i  # add, then stabilize
    {1: 2, 2: 2, 3: 2, 4: 1}
    sage: e = m + m
    sage: e
    {1: 4, 2: 4, 3: 4, 4: 2}
    sage: ~e   # stabilize
    {1: 2, 2: 2, 3: 2, 4: 0}
    sage: a = -m
    sage: a & m
    {1: 0, 2: 0, 3: 0, 4: 0}
    sage: a * m   # add, then find the equivalent recurrent
    {1: 2, 2: 2, 3: 2, 4: 0}
    sage: a^3  # a*a*a
    {1: 2, 2: 2, 3: 2, 4: 1}
    sage: a^(-1) == m
    True
    sage: a < m  # every coordinate of a is < that of m
    True

Firing an unstable vertex returns resulting configuration::

    sage: c = S.max_stable() + S.identity()
    sage: c.fire_vertex(1)
    {1: 1, 2: 5, 3: 5, 4: 1}
    sage: c
    {1: 4, 2: 4, 3: 4, 4: 1}

Fire all unstable vertices::

    sage: c.unstable()
    [1, 2, 3]
    sage: c.fire_unstable()
    {1: 3, 2: 3, 3: 3, 4: 3}

Stabilize c, returning the resulting configuration and the firing
vector::

    sage: c.stabilize(True)
    [{1: 2, 2: 2, 3: 2, 4: 1}, {1: 6, 2: 8, 3: 8, 4: 8}]
    sage: c
    {1: 4, 2: 4, 3: 4, 4: 1}
    sage: S.max_stable() & S.identity() == c.stabilize()
    True

The number of superstable configurations of each degree::

    sage: S.hilbert_function()
    [1, 4, 8]
    sage: S.postulation()
    2

the saturated, homogeneous sandpile ideal

    sage: S.ideal()
    Ideal (x1 - x0, x3*x2 - x0^2, x4^2 - x0^2, x2^3 - x4*x3*x0, x4*x2^2 - x3^2*x0, x3^3 - x4*x2*x0, x4*x3^2 - x2^2*x0) of Multivariate Polynomial Ring in x4, x3, x2, x1, x0 over Rational Field

its minimal free resolution

    sage: S.resolution()
    'R^1 <-- R^7 <-- R^15 <-- R^13 <-- R^4'

and its Betti numbers::

    sage: S.betti()
               0     1     2     3     4
    ------------------------------------
        0:     1     1     -     -     -
        1:     -     2     2     -     -
        2:     -     4    13    13     4
    ------------------------------------
    total:     1     7    15    13     4

Distribution of avalanche sizes::

    sage: S = grid_sandpile(10,10)
    sage: m = S.max_stable()
    sage: a = []
    sage: for i in range(1000):
    ...       m = m.add_random()
    ...       m, f = m.stabilize(True)
    ...       a.append(sum(f.values()))
    ...
    sage: p = list_plot([[log(i+1),log(a.count(i))] for i in [0..max(a)] if a.count(i)])
    sage: p.axes_labels(['log(N)','log(D(N))'])
    sage: t = text("Distribution of avalanche sizes", (2,2), rgbcolor=(1,0,0))
    sage: show(p+t,axes_labels=['log(N)','log(D(N))'])

To calculate linear systems associated with divisors, 4ti2 must be installed.
One way to do this is to run sage -i to install glpk, then 4ti2.  See
http://sagemath.org/download-packages.html to get the exact names of these
packages.  An alternative is to install 4ti2 separately, then point the
following variable to the correct path.
"""

#*****************************************************************************
#       Copyright (C) 2011 David Perkinson <davidp@reed.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from string import join
import os
from sage.symbolic.all import I, pi
from sage.functions.log import exp
from sage.graphs.all import DiGraph, Graph, graphs, digraphs
from copy import deepcopy
from sage.rings.all import PolynomialRing, QQ, ZZ, lcm
from sage.misc.all import prod, det, forall, tmp_filename, random, randint, exists, denominator, srange
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, identity_matrix
from sage.interfaces.singular import singular
from sage.combinat.combinat import CombinatorialClass
from sage.combinat.set_partition import SetPartitions
from sage.homology.simplicial_complex import SimplicialComplex
from sage.plot.colors import rainbow
from sage.env import SAGE_LOCAL

path_to_zsolve = os.path.join(SAGE_LOCAL,'bin','zsolve')

class Sandpile(DiGraph):
    """
    Class for Dhar's abelian sandpile model.
    """

    def __init__(self, g, sink):
        r"""
        Create a sandpile.

        A sandpile is always a weighted graph.

        INPUT:

         - ``g`` - dict for directed multigraph (see NOTES) edges weighted by
           nonnegative integers

         - ``sink`` - A sink vertex.  Any outgoing edges from the designated
           sink are ignored for the purposes of stabilization.  It is assumed
           that every vertex has a directed path into the sink.

        OUTPUT:

        Sandpile

        EXAMPLES:

        Here ``g`` represents a square with directed, multiple edges with three
        vertices, ``a``, ``b``, ``c``, and ``d``.  The vertex ``a`` has
        outgoing edges to itself (weight 2), to vertex ``b`` (weight 1), and
        vertex ``c`` (weight 3), for example.

        ::

            sage: g = {'a': {'a':2, 'b':1, 'c':3}, 'b': {'a':1, 'd':1},\
                       'c': {'a':1,'d': 1}, 'd': {'b':1, 'c':1}}
            sage: G = Sandpile(g,'d')

        Here is a square with unweighted edges.  In this example, the graph is
        also undirected.

        ::

            sage: g = {0:[1,2], 1:[0,3], 2:[0,3], 3:[1,2]}
            sage: G = Sandpile(g,3)

        NOTES::

        Loops are allowed.  There are two admissible formats, both of which are
        dictionaries whose keys are the vertex names.  In one format, the
        values are dictionaries with keys the names of vertices which are the
        tails of outgoing edges and with values the weights of the edges.  In
        the other format, the values are lists of names of vertices which are
        the tails of the outgoing edges.  All weights are then automatically
        assigned the value 1.

        TESTS::

            sage: S = complete_sandpile(4)
            sage: TestSuite(S).run()

        Make sure we cannot make an unweighted sandpile::

            sage: G = Sandpile({0:[]}, 0, weighted=False)
            Traceback (most recent call last):
            ...
            TypeError: __init__() got an unexpected keyword argument 'weighted'
        """
        # preprocess a graph, if necessary
        if isinstance(g, dict) and isinstance(g.values()[0], dict):
            pass # this is the default format
        elif isinstance(g, dict) and isinstance(g.values()[0], list):
            processed_g = {}
            for k in g.keys():
                temp = {}
                for vertex in g[k]:
                    temp[vertex] = 1
                processed_g[k] = temp
            g = processed_g
        elif isinstance(g, Graph):
            processed_g = {}
            for v in g.vertices():
                edges = {}
                for n in g.neighbors(v):
                    if isinstance(g.edge_label(v,n), type(1)) and g.edge_label(v,n) >=0:
                        edges[n] = g.edge_label(v,n)
                    else:
                        edges[n] = 1
                processed_g[v] = edges
            g = processed_g
        elif isinstance(g, DiGraph):
            processed_g = {}
            for v in g.vertices():
                edges = {}
                for n in g.neighbors_out(v):
                    if (isinstance(g.edge_label(v,n), type(1))
                        and g.edge_label(v,n)>=0):
                        edges[n] = g.edge_label(v,n)
                    else:
                        edges[n] = 1
                processed_g[v] = edges
            g = processed_g
        else:
            raise SyntaxError(g)

        # create digraph and initialize some variables
        DiGraph.__init__(self,g,weighted=True)
        self._dict = deepcopy(g)
        self._sink = sink  # key for sink
        self._sink_ind = self.vertices().index(sink)
        self._nonsink_vertices = deepcopy(self.vertices())
        del self._nonsink_vertices[self._sink_ind]
        # compute laplacians:
        self._laplacian = self.laplacian_matrix(indegree=False)
        temp = range(self.num_verts())
        del temp[self._sink_ind]
        self._reduced_laplacian = self._laplacian[temp,temp]

    def __copy__(self):
        """
        Make a copy of this sandpile

        OUTPUT:

        A new :class:`Sandpile` instance.

        EXAMPLES::

            sage: G = complete_sandpile(4)
            sage: G_copy = copy(G)
            sage: G_copy == G == G.__copy__()
            True
        """
        return self.__class__(self, self._sink)

    def __getattr__(self, name):
        """
        Set certain variables only when called.

        INPUT:

        ``name`` - name of an internal method

        OUTPUT:

        None.

        EXAMPLES::

            sage: S = complete_sandpile(5)
            sage: S.__getattr__('_max_stable')
            {1: 3, 2: 3, 3: 3, 4: 3}
        """
        if name not in self.__dict__:
            if name == '_max_stable':
                self._set_max_stable()
                return deepcopy(self.__dict__[name])
            if name == '_max_stable_div':
                self._set_max_stable_div()
                return deepcopy(self.__dict__[name])
            elif name == '_out_degrees':
                self._set_out_degrees()
                return deepcopy(self.__dict__[name])
            elif name == '_in_degrees':
                self._set_in_degrees()
                return deepcopy(self.__dict__[name])
            elif name == '_burning_config' or name == '_burning_script':
                self._set_burning_config()
                return deepcopy(self.__dict__[name])
            elif name == '_identity':
                self._set_identity()
                return deepcopy(self.__dict__[name])
            elif name == '_recurrents':
                self._set_recurrents()
                return deepcopy(self.__dict__[name])
            elif name == '_min_recurrents':
                self._set_min_recurrents()
                return deepcopy(self.__dict__[name])
            elif name == '_superstables':
                self._set_superstables()
                return deepcopy(self.__dict__[name])
            elif name == '_group_order':
                self.__dict__[name] = det(self._reduced_laplacian.dense_matrix())
                return self.__dict__[name]
            elif name == '_invariant_factors':
                self._set_invariant_factors()
                return deepcopy(self.__dict__[name])
            elif name == '_betti_complexes':
                self._set_betti_complexes()
                return deepcopy(self.__dict__[name])
            elif (name == '_postulation' or name == '_h_vector'
                   or name == '_hilbert_function'):
                self._set_hilbert_function()
                return deepcopy(self.__dict__[name])
            elif (name == '_ring' or name == '_unsaturated_ideal'):
                self._set_ring()
                return self.__dict__[name]
            elif name == '_ideal':
                self._set_ideal()
                return self.__dict__[name]
            elif (name == '_resolution' or name == '_betti' or name ==
            '_singular_resolution'):
                self._set_resolution()
                return self.__dict__[name]
            elif name == '_groebner':
                self._set_groebner()
                return self.__dict__[name]
            elif name == '_points':
                self._set_points()
                return self.__dict__[name]
            else:
                raise AttributeError(name)

    def version(self):
        r"""
        The version number of Sage Sandpiles

        INPUT:

        None

        OUTPUT:

        string


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.version()
            Sage Sandpiles Version 2.3
        """
        print 'Sage Sandpiles Version 2.3'

    def show(self,**kwds):
        r"""
        Draws the graph.

        INPUT:

        ``kwds`` - arguments passed to the show method for Graph or DiGraph

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.show()
            sage: S.show(graph_border=True, edge_labels=True)
        """

        if self.is_undirected():
            Graph(self).show(**kwds)
        else:
            DiGraph(self).show(**kwds)

    def show3d(self,**kwds):
        r"""
        Draws the graph.

        INPUT:

        ``kwds`` - arguments passed to the show method for Graph or DiGraph

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.show3d()
        """

        if self.is_undirected():
            Graph(self).show3d(**kwds)
        else:
            DiGraph(self).show3d(**kwds)

    def dict(self):
        r"""
        A dictionary of dictionaries representing a directed graph.

        INPUT:

        None

        OUTPUT:

        dict


        EXAMPLES::

            sage: G = sandlib('generic')
            sage: G.dict()
            {0: {},
             1: {0: 1, 3: 1, 4: 1},
             2: {0: 1, 3: 1, 5: 1},
             3: {2: 1, 5: 1},
             4: {1: 1, 3: 1},
             5: {2: 1, 3: 1}}
            sage: G.sink()
            0
        """
        return deepcopy(self._dict)

    def sink(self):
        r"""
        The identifier for the sink vertex.

        INPUT:

        None

        OUTPUT:

        Object (name for the sink vertex)

        EXAMPLES::

            sage: G = sandlib('generic')
            sage: G.sink()
            0
            sage: H = grid_sandpile(2,2)
            sage: H.sink()
            'sink'
            sage: type(H.sink())
            <type 'str'>
        """
        return self._sink

    def laplacian(self):
        r"""
        The Laplacian matrix of the graph.

        INPUT:

        None

        OUTPUT:

        matrix


        EXAMPLES::

            sage: G = sandlib('generic')
            sage: G.laplacian()
            [ 0  0  0  0  0  0]
            [-1  3  0 -1 -1  0]
            [-1  0  3 -1  0 -1]
            [ 0  0 -1  2  0 -1]
            [ 0 -1  0 -1  2  0]
            [ 0  0 -1 -1  0  2]

        NOTES:

        The function ``laplacian_matrix`` should be avoided.  It returns the
        indegree version of the laplacian.
        """
        return deepcopy(self._laplacian)

    def reduced_laplacian(self):
        r"""
        The reduced Laplacian matrix of the graph.

        INPUT:

        None

        OUTPUT:

        matrix


        EXAMPLES::

            sage: G = sandlib('generic')
            sage: G.laplacian()
            [ 0  0  0  0  0  0]
            [-1  3  0 -1 -1  0]
            [-1  0  3 -1  0 -1]
            [ 0  0 -1  2  0 -1]
            [ 0 -1  0 -1  2  0]
            [ 0  0 -1 -1  0  2]
            sage: G.reduced_laplacian()
            [ 3  0 -1 -1  0]
            [ 0  3 -1  0 -1]
            [ 0 -1  2  0 -1]
            [-1  0 -1  2  0]
            [ 0 -1 -1  0  2]

        NOTES:

        This is the Laplacian matrix with the row and column indexed by the
        sink vertex removed.
        """
        return deepcopy(self._reduced_laplacian)

    def group_order(self):
        r"""
        The size of the sandpile group.

        INPUT:

        None

        OUTPUT:

        int

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.group_order()
            15
        """
        return self._group_order

    def _set_max_stable(self):
        r"""
        Initialize the variable holding the maximal stable configuration.

        INPUT:

        None

        OUTPUT:

        NONE

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_max_stable()
        """
        m = {}
        for v in self._nonsink_vertices:
            m[v] = self.out_degree(v)-1
        self._max_stable = SandpileConfig(self,m)

    def max_stable(self):
        r"""
        The maximal stable configuration.

        INPUT:

        None

        OUTPUT:

        SandpileConfig (the maximal stable configuration)


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.max_stable()
            {1: 2, 2: 2, 3: 1, 4: 1, 5: 1}
        """
        return deepcopy(self._max_stable)

    def _set_max_stable_div(self):
        r"""
        Initialize the variable holding the maximal stable divisor.

        INPUT:

        None

        OUTPUT:

        NONE

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_max_stable_div()
        """
        m = {}
        for v in self.vertices():
            m[v] = self.out_degree(v)-1
        self._max_stable_div = SandpileDivisor(self,m)

    def max_stable_div(self):
        r"""
        The maximal stable divisor.

        INPUT:

        SandpileDivisor

        OUTPUT:

        SandpileDivisor (the maximal stable divisor)


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.max_stable_div()
            {0: -1, 1: 2, 2: 2, 3: 1, 4: 1, 5: 1}
            sage: S.out_degree()
            {0: 0, 1: 3, 2: 3, 3: 2, 4: 2, 5: 2}
        """
        return deepcopy(self._max_stable_div)

    def _set_out_degrees(self):
        r"""
        Initialize the variable holding the out-degrees.

        INPUT:

        None

        OUTPUT:

        NONE

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_out_degrees()
        """
        self._out_degrees = dict(self.zero_div())
        for v in self.vertices():
            self._out_degrees[v] = 0
            for e in self.edges_incident(v):
                self._out_degrees[v] += e[2]

    def out_degree(self, v=None):
        r"""
        The out-degree of a vertex or a list of all out-degrees.

        INPUT:

        ``v`` (optional) - vertex name

        OUTPUT:

        integer or dict

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.out_degree(2)
            3
            sage: S.out_degree()
            {0: 0, 1: 3, 2: 3, 3: 2, 4: 2, 5: 2}
        """
        if not v is None:
            return self._out_degrees[v]
        else:
            return self._out_degrees

    def _set_in_degrees(self):
        """
        Initialize the variable holding the in-degrees.

        INPUT:

        None

        OUTPUT:

        NONE

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_in_degrees()
        """
        self._in_degrees = dict(self.zero_div())
        for e in self.edges():
            self._in_degrees[e[1]] += e[2]

    def in_degree(self, v=None):
        r"""
        The in-degree of a vertex or a list of all in-degrees.

        INPUT:

        ``v`` - vertex name or None

        OUTPUT:

        integer or dict

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.in_degree(2)
            2
            sage: S.in_degree()
            {0: 2, 1: 1, 2: 2, 3: 4, 4: 1, 5: 2}
        """
        if not v is None:
            return self._in_degrees[v]
        else:
            return self._in_degrees

    def _set_burning_config(self):
        r"""
        Calculate the minimal burning configuration and its corresponding
        script.

        EXAMPLES::

            sage: g = {0:{},1:{0:1,3:1,4:1},2:{0:1,3:1,5:1}, \
                       3:{2:1,5:1},4:{1:1,3:1},5:{2:1,3:1}}
            sage: S = Sandpile(g,0)
            sage: S._set_burning_config()
        """
        # TODO: Cythonize!
        d = self._reduced_laplacian.nrows()
        burn = sum(self._reduced_laplacian)
        script=[1]*d  # d 1s
        done = False
        while not done:
            bad = -1
            for i in range(d):
                if burn[i] < 0:
                    bad = i
                    break
            if bad == -1:
                done = True
            else:
                burn += self._reduced_laplacian[bad]
                script[bad]+=1
        b = iter(burn)
        s = iter(script)
        bc = {} # burning config
        bs = {} # burning script
        for v in self._nonsink_vertices:
            bc[v] = next(b)
            bs[v] = next(s)
        self._burning_config = SandpileConfig(self,bc)
        self._burning_script = SandpileConfig(self,bs)

    def burning_config(self):
        r"""
        A minimal burning configuration.

        INPUT:

        None

        OUTPUT:

        dict (configuration)

        EXAMPLES::

            sage: g = {0:{},1:{0:1,3:1,4:1},2:{0:1,3:1,5:1}, \
                       3:{2:1,5:1},4:{1:1,3:1},5:{2:1,3:1}}
            sage: S = Sandpile(g,0)
            sage: S.burning_config()
            {1: 2, 2: 0, 3: 1, 4: 1, 5: 0}
            sage: S.burning_config().values()
            [2, 0, 1, 1, 0]
            sage: S.burning_script()
            {1: 1, 2: 3, 3: 5, 4: 1, 5: 4}
            sage: script = S.burning_script().values()
            sage: script
            [1, 3, 5, 1, 4]
            sage: matrix(script)*S.reduced_laplacian()
            [2 0 1 1 0]

        NOTES:

        The burning configuration and script are computed using a modified
        version of Speer's script algorithm.  This is a generalization to
        directed multigraphs of Dhar's burning algorithm.

        A *burning configuration* is a nonnegative integer-linear
        combination of the rows of the reduced Laplacian matrix having
        nonnegative entries and such that every vertex has a path from some
        vertex in its support.  The corresponding *burning script* gives
        the integer-linear combination needed to obtain the burning
        configuration.  So if `b` is the burning configuration, `\sigma` is its
        script, and `\tilde{L}` is the reduced Laplacian, then `\sigma\cdot
        \tilde{L} = b`.  The *minimal burning configuration* is the one
        with the minimal script (its components are no larger than the
        components of any other script
        for a burning configuration).

        The following are equivalent for a configuration `c` with burning
        configuration `b` having script `\sigma`:

         - `c` is recurrent;
         - `c+b` stabilizes to `c`;
         - the firing vector for the stabilization of `c+b` is `\sigma`.
        """
        return deepcopy(self._burning_config)

    def burning_script(self):
        r"""
        A script for the minimal burning configuration.

        INPUT:

        None

        OUTPUT:

        dict

        EXAMPLES::

            sage: g = {0:{},1:{0:1,3:1,4:1},2:{0:1,3:1,5:1},\
            3:{2:1,5:1},4:{1:1,3:1},5:{2:1,3:1}}
            sage: S = Sandpile(g,0)
            sage: S.burning_config()
            {1: 2, 2: 0, 3: 1, 4: 1, 5: 0}
            sage: S.burning_config().values()
            [2, 0, 1, 1, 0]
            sage: S.burning_script()
            {1: 1, 2: 3, 3: 5, 4: 1, 5: 4}
            sage: script = S.burning_script().values()
            sage: script
            [1, 3, 5, 1, 4]
            sage: matrix(script)*S.reduced_laplacian()
            [2 0 1 1 0]

        NOTES:

        The burning configuration and script are computed using a modified
        version of Speer's script algorithm.  This is a generalization to
        directed multigraphs of Dhar's burning algorithm.

        A *burning configuration* is a nonnegative integer-linear
        combination of the rows of the reduced Laplacian matrix having
        nonnegative entries and such that every vertex has a path from some
        vertex in its support.  The corresponding *burning script* gives the
        integer-linear combination needed to obtain the burning configuration.
        So if `b` is the burning configuration, `s` is its script, and
        `L_{\mathrm{red}}` is the reduced Laplacian, then `s\cdot
        L_{\mathrm{red}}= b`.  The *minimal burning configuration* is the one
        with the minimal script (its components are no larger than the
        components of any other script
        for a burning configuration).

        The following are equivalent for a configuration `c` with burning
        configuration `b` having script `s`:

         - `c` is recurrent;
         - `c+b` stabilizes to `c`;
         - the firing vector for the stabilization of `c+b` is `s`.
        """
        return deepcopy(self._burning_script)

    def nonsink_vertices(self):
        r"""
        The names of the nonsink vertices.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.nonsink_vertices()
            [1, 2, 3, 4, 5]
        """
        return self._nonsink_vertices

    def all_k_config(self,k):
        r"""
        The configuration with all values set to k.

        INPUT:

        ``k`` - integer

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.all_k_config(7)
            {1: 7, 2: 7, 3: 7, 4: 7, 5: 7}
        """
        return SandpileConfig(self,[k]*(self.num_verts()-1))

    def zero_config(self):
        r"""
        The all-zero configuration.

        INPUT:

        None

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.zero_config()
            {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
        """
        return self.all_k_config(0)

    # TODO: cythonize stabilization!
    # The following would presumably be moved to the SandpileConfig class
    #def new_stabilize(self, config):
    #    r"""
    #    Stabilize \code{config}, returning \code{[out_config, firing_vector]},
    #    where \code{out_config} is the modified configuration.
    #    """
    #    c, f = cython_stabilize(config, self.reduced_laplacian(),
    #        self.out_degree(), self.nonsink_vertices())
    #    self._config = c
    #    return [c, f]

    def _set_identity(self):
        r"""
        Computes ``_identity``, the variable holding the identity configuration
        of the sandpile group, when ``identity()`` is first called by a user.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_identity()
        """
        m = self._max_stable
        self._identity = (m&m).dualize()&m

    def identity(self):
        r"""
        The identity configuration.

        INPUT:

        None

        OUTPUT:

        dict (the identity configuration)

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: e = S.identity()
            sage: x = e & S.max_stable()  # stable addition
            sage: x
            {1: 2, 2: 2, 3: 1, 4: 1, 5: 1}
            sage: x == S.max_stable()
            True
        """
        return deepcopy(self._identity)

    def _set_recurrents(self):
        """
        Computes ``_recurrents``, the variable holding the list of recurrent
        configurations, when ``recurrents()`` is first called by a user.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_recurrents()
        """
        self._recurrents = []
        active = [self._max_stable]
        while active != []:
            c = active.pop()
            self._recurrents.append(c)
            for v in self._nonsink_vertices:
                cnext = deepcopy(c)
                cnext[v] += 1
                cnext = ~cnext
                if (cnext not in active) and (cnext not in self._recurrents):
                    active.append(cnext)

    def recurrents(self, verbose=True):
        r"""
        The list of recurrent configurations. If ``verbose`` is
        ``False``, the configurations are converted to lists of
        integers.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        list (of recurrent configurations)


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.recurrents()
            [{1: 2, 2: 2, 3: 1, 4: 1, 5: 1}, {1: 2, 2: 2, 3: 0, 4: 1, 5: 1}, {1: 0, 2: 2, 3: 1, 4: 1, 5: 0}, {1: 0, 2: 2, 3: 1, 4: 1, 5: 1}, {1: 1, 2: 2, 3: 1, 4: 1, 5: 1}, {1: 1, 2: 2, 3: 0, 4: 1, 5: 1}, {1: 2, 2: 2, 3: 1, 4: 0, 5: 1}, {1: 2, 2: 2, 3: 0, 4: 0, 5: 1}, {1: 2, 2: 2, 3: 1, 4: 0, 5: 0}, {1: 1, 2: 2, 3: 1, 4: 1, 5: 0}, {1: 1, 2: 2, 3: 1, 4: 0, 5: 0}, {1: 1, 2: 2, 3: 1, 4: 0, 5: 1}, {1: 0, 2: 2, 3: 0, 4: 1, 5: 1}, {1: 2, 2: 2, 3: 1, 4: 1, 5: 0}, {1: 1, 2: 2, 3: 0, 4: 0, 5: 1}]
            sage: S.recurrents(verbose=False)
            [[2, 2, 1, 1, 1], [2, 2, 0, 1, 1], [0, 2, 1, 1, 0], [0, 2, 1, 1, 1], [1, 2, 1, 1, 1], [1, 2, 0, 1, 1], [2, 2, 1, 0, 1], [2, 2, 0, 0, 1], [2, 2, 1, 0, 0], [1, 2, 1, 1, 0], [1, 2, 1, 0, 0], [1, 2, 1, 0, 1], [0, 2, 0, 1, 1], [2, 2, 1, 1, 0], [1, 2, 0, 0, 1]]
        """
        if verbose:
            return deepcopy(self._recurrents)
        else:
            return [r.values() for r in self._recurrents]

    def _set_superstables(self):
        r"""
        Computes ``_superstables``, the variable holding the list of superstable
        configurations, when ``superstables()`` is first called by a user.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_superstables()
        """
        self._superstables = [c.dualize() for c in self._recurrents]

    def superstables(self, verbose=True):
        r"""
        The list of superstable configurations as dictionaries if
        ``verbose`` is ``True``, otherwise as lists of integers.  The
        superstables are also known as `G`-parking functions.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        list (of superstable elements)


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.superstables()
            [{1: 0, 2: 0, 3: 0, 4: 0, 5: 0},
             {1: 0, 2: 0, 3: 1, 4: 0, 5: 0},
             {1: 2, 2: 0, 3: 0, 4: 0, 5: 1},
             {1: 2, 2: 0, 3: 0, 4: 0, 5: 0},
             {1: 1, 2: 0, 3: 0, 4: 0, 5: 0},
             {1: 1, 2: 0, 3: 1, 4: 0, 5: 0},
             {1: 0, 2: 0, 3: 0, 4: 1, 5: 0},
             {1: 0, 2: 0, 3: 1, 4: 1, 5: 0},
             {1: 0, 2: 0, 3: 0, 4: 1, 5: 1},
             {1: 1, 2: 0, 3: 0, 4: 0, 5: 1},
             {1: 1, 2: 0, 3: 0, 4: 1, 5: 1},
             {1: 1, 2: 0, 3: 0, 4: 1, 5: 0},
             {1: 2, 2: 0, 3: 1, 4: 0, 5: 0},
             {1: 0, 2: 0, 3: 0, 4: 0, 5: 1},
             {1: 1, 2: 0, 3: 1, 4: 1, 5: 0}]
            sage: S.superstables(False)
            [[0, 0, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [2, 0, 0, 0, 1],
             [2, 0, 0, 0, 0],
             [1, 0, 0, 0, 0],
             [1, 0, 1, 0, 0],
             [0, 0, 0, 1, 0],
             [0, 0, 1, 1, 0],
             [0, 0, 0, 1, 1],
             [1, 0, 0, 0, 1],
             [1, 0, 0, 1, 1],
             [1, 0, 0, 1, 0],
             [2, 0, 1, 0, 0],
             [0, 0, 0, 0, 1],
             [1, 0, 1, 1, 0]]
        """
        if verbose:
            return deepcopy(self._superstables)
        else:
            verts = self.nonsink_vertices()
            return [s.values() for s in self._superstables]

    def is_undirected(self):
        r"""
        ``True`` if ``(u,v)`` is and edge if and only if ``(v,u)`` is an
        edges, each edge with the same weight.

        INPUT:

        None

        OUTPUT:

        boolean

        EXAMPLES::

            sage: complete_sandpile(4).is_undirected()
            True
            sage: sandlib('gor').is_undirected()
            False
        """
        return self.laplacian().is_symmetric()

    def _set_min_recurrents(self):
        r"""
        Computes the minimal recurrent elements.  If the underlying graph is
        undirected, these are the recurrent elements of least degree.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: complete_sandpile(4)._set_min_recurrents()
        """
        if self.is_undirected():
            m = min([r.deg() for r in self.recurrents()])
            rec = [r for r in self.recurrents() if r.deg()==m]
        else:
            rec = self.recurrents()
            for r in self.recurrents():
                if exists(rec, lambda x: r>x)[0]:
                    rec.remove(r)
        self._min_recurrents = rec

    def min_recurrents(self, verbose=True):
        r"""
        The minimal recurrent elements.  If the underlying graph is
        undirected, these are the recurrent elements of least degree.  If
        ``verbose is ``False``, the configurations are converted to lists of
        integers.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        list of SandpileConfig

        EXAMPLES::

            sage: S=sandlib('riemann-roch2')
            sage: S.min_recurrents()
            [{1: 0, 2: 0, 3: 1}, {1: 1, 2: 1, 3: 0}]
            sage: S.min_recurrents(False)
            [[0, 0, 1], [1, 1, 0]]
            sage: S.recurrents(False)
            [[1, 1, 2],
             [0, 1, 1],
             [0, 1, 2],
             [1, 0, 1],
             [1, 0, 2],
             [0, 0, 2],
             [1, 1, 1],
             [0, 0, 1],
             [1, 1, 0]]
            sage: [i.deg() for i in S.recurrents()]
            [4, 2, 3, 2, 3, 2, 3, 1, 2]
        """
        if verbose:
            return deepcopy(self._min_recurrents)
        else:
            return [r.values() for r in self._min_recurrents]

    def max_superstables(self, verbose=True):
        r"""
        The maximal superstable configurations.  If the underlying graph is
        undirected, these are the superstables of highest degree.  If
        ``verbose`` is ``False``, the configurations are converted to lists of
        integers.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        list (of maximal superstables)

        EXAMPLES::

            sage: S=sandlib('riemann-roch2')
            sage: S.max_superstables()
            [{1: 1, 2: 1, 3: 1}, {1: 0, 2: 0, 3: 2}]
            sage: S.superstables(False)
            [[0, 0, 0],
             [1, 0, 1],
             [1, 0, 0],
             [0, 1, 1],
             [0, 1, 0],
             [1, 1, 0],
             [0, 0, 1],
             [1, 1, 1],
             [0, 0, 2]]
            sage: S.h_vector()
            [1, 3, 4, 1]
        """
        result = [r.dualize() for r in self.min_recurrents()]
        if verbose:
            return result
        else:
            return [r.values() for r in result]


    def nonspecial_divisors(self, verbose=True):
        r"""
        The nonspecial divisors: those divisors of degree ``g-1`` with
        empty linear system.  The term is only defined for undirected graphs.
        Here, ``g = |E| - |V| + 1`` is the genus of the graph.  If ``verbose``
        is ``False``, the divisors are converted to lists of integers.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        list (of divisors)

        EXAMPLES::

            sage: S = complete_sandpile(4)
            sage: ns = S.nonspecial_divisors() # optional - 4ti2
            sage: D = ns[0] # optional - 4ti2
            sage: D.values() # optional - 4ti2
            [-1, 1, 0, 2]
            sage: D.deg() # optional - 4ti2
            2
            sage: [i.effective_div() for i in ns] # optional - 4ti2
            [[], [], [], [], [], []]
        """
        if self.is_undirected():
            result = []
            for s in self.max_superstables():
                D = dict(s)
                D[self._sink] = -1
                D = SandpileDivisor(self, D)
                result.append(D)
            if verbose:
                return result
            else:
                return [r.values() for r in result]
        else:
            raise UserWarning("The underlying graph must be undirected.")

    def canonical_divisor(self):
        r"""
        The canonical divisor: the divisor ``deg(v)-2`` grains of sand
        on each vertex.  Only for undirected graphs.

        INPUT:

        None

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = complete_sandpile(4)
            sage: S.canonical_divisor()
            {0: 1, 1: 1, 2: 1, 3: 1}
        """
        if self.is_undirected():
            return SandpileDivisor(self,[self.out_degree(v)-2 for v in self.vertices()])
        else:
            raise UserWarning("Only for undirected graphs.")

    def _set_invariant_factors(self):
        r"""
        Computes the variable holding the elementary divisors of the sandpile
        group when ``invariant_factors()`` is first called by the user.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_invariant_factors()
        """
        # Sage seems to have issues with computing the Smith normal form and
        # elementary divisors of a sparse matrix, so we have to convert:
        A = self.reduced_laplacian().dense_matrix()
        self._invariant_factors = A.elementary_divisors()

    def invariant_factors(self):
        r"""
        The invariant factors of the sandpile group (a finite
        abelian group).

        INPUT:

        None

        OUTPUT:

        list of integers

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.invariant_factors()
            [1, 1, 1, 1, 15]
        """
        return deepcopy(self._invariant_factors)

    def _set_hilbert_function(self):
        """
        Computes the variables holding the Hilbert function of the homogeneous
        homogeneous sandpile ideal, the first differences of the Hilbert
        function, and the postulation number for the zero-set of the sandpile
        ideal when any one of these is called by the user.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_hilbert_function()
        """
        v = [i.deg() for i in self._superstables]
        self._postulation = max(v)
        self._h_vector = [v.count(i) for i in range(self._postulation+1)]
        self._hilbert_function = [1]
        for i in range(self._postulation):
            self._hilbert_function.append(self._hilbert_function[i]
                +self._h_vector[i+1])

    def h_vector(self):
        r"""
        The first differences of the Hilbert function of the homogeneous
        sandpile ideal.  It lists the number of superstable configurations in
        each degree.

        INPUT:

        None

        OUTPUT:

        list of nonnegative integers


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.hilbert_function()
            [1, 5, 11, 15]
            sage: S.h_vector()
            [1, 4, 6, 4]
        """
        return deepcopy(self._h_vector)

    def hilbert_function(self):
        r"""
        The Hilbert function of the homogeneous sandpile ideal.

        INPUT:

        None

        OUTPUT:

        list of nonnegative integers

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.hilbert_function()
            [1, 5, 11, 15]
        """
        return deepcopy(self._hilbert_function)

    def postulation(self):
        r"""
        The postulation number of the sandpile ideal.  This is the
        largest weight of a superstable configuration of the graph.

        INPUT:

        None

        OUTPUT:

        nonnegative integer

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.postulation()
            3
        """
        return self._postulation

    def reorder_vertices(self):
        r"""
        Create a copy of the sandpile but with the vertices ordered according
        to their distance from the sink, from greatest to least.

        INPUT:

        None

        OUTPUT:

        Sandpile

        EXAMPLES::
            sage: S = sandlib('kite')
            sage: S.dict()
            {0: {}, 1: {0: 1, 2: 1, 3: 1}, 2: {1: 1, 3: 1, 4: 1}, 3: {1: 1, 2: 1, 4: 1}, 4: {2: 1, 3: 1}}
            sage: T = S.reorder_vertices()
            sage: T.dict()
            {0: {1: 1, 2: 1}, 1: {0: 1, 2: 1, 3: 1}, 2: {0: 1, 1: 1, 3: 1}, 3: {1: 1, 2: 1, 4: 1}, 4: {}}
        """

        # first order the vertices according to their distance from the sink
        verts = self.vertices()
        verts = sorted(verts, self._compare_vertices)
        verts.reverse()
        perm = {}
        for i in range(len(verts)):
            perm[verts[i]]=i
        old = self.dict()
        new = {}
        for i in old:
            entry = {}
            for j in old[i]:
                entry[perm[j]]=old[i][j]
            new[perm[i]] = entry
        return Sandpile(new,len(verts)-1)


#################### Functions for divisors #####################

    def all_k_div(self,k):
        r"""
        The divisor with all values set to k.

        INPUT:

        ``k`` - integer

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.all_k_div(7)
            {0: 7, 1: 7, 2: 7, 3: 7, 4: 7, 5: 7}
        """
        return SandpileDivisor(self,[k]*self.num_verts())

    def zero_div(self):
        r"""
        The all-zero divisor.

        INPUT:

        None

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.zero_div()
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
        """
        return self.all_k_div(0)

    # FIX: save this in the __dict__
    def _set_betti_complexes(self):
        r"""
        Compute the value return by the ``betti_complexes`` method.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = Sandpile({0:{},1:{0: 1, 2: 1, 3: 4},2:{3: 5},3:{1: 1, 2: 1}},0)
            sage: S._set_betti_complexes() # optional - 4ti2
        """
        results = []
        verts = self.vertices()
        r = self.recurrents()
        for D in r:
            d = D.deg()
            # change D to a dict since SandpileConfig will not allow adding a key
            D = dict(D)
            D[self.sink()] = -d
            D = SandpileDivisor(self,D)
            test = True
            while test:
                D[self.sink()] += 1
                complex = D.Dcomplex()
                if sum(complex.betti().values()) > 1:  # change from 0 to 1
                    results.append([deepcopy(D), complex])
                if len(complex.maximal_faces()) == 1 and list(complex.maximal_faces()[0]) == verts:
                    test = False
        self._betti_complexes = results

    def betti_complexes(self):
        r"""
        A list of all the divisors with nonempty linear systems whose
        corresponding simplicial complexes have nonzero homology in some
        dimension. Each such divisors is returned with its corresponding
        simplicial complex.

        INPUT:

        None

        OUTPUT:

        list (of pairs [divisors, corresponding simplicial complex])


        EXAMPLES::

            sage: S = Sandpile({0:{},1:{0: 1, 2: 1, 3: 4},2:{3: 5},3:{1: 1, 2: 1}},0)
            sage: p = S.betti_complexes() # optional - 4ti2
            sage: p[0] # optional - 4ti2
            [{0: -8, 1: 5, 2: 4, 3: 1}, Simplicial complex with vertex set (1, 2, 3) and facets {(1, 2), (3,)}]
            sage: S.resolution() # optional - 4ti2
            'R^1 <-- R^5 <-- R^5 <-- R^1'
            sage: S.betti()
                       0     1     2     3
            ------------------------------
                0:     1     -     -     -
                1:     -     5     5     -
                2:     -     -     -     1
            ------------------------------
            total:     1     5     5     1
            sage: len(p) # optional - 4ti2
            11
            sage: p[0][1].homology() # optional - 4ti2
            {0: Z, 1: 0}
            sage: p[-1][1].homology() # optional - 4ti2
            {0: 0, 1: 0, 2: Z}
        """
        return deepcopy(self._betti_complexes)

#######################################
######### Algebraic Geometry ##########
#######################################

    def _compare_vertices(self, v, w):
        r"""
        Compare vertices based on their distance from the sink.

        INPUT:

        ``v``, ``w`` - vertices

        OUTPUT:

        integer

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.vertices()
            [0, 1, 2, 3, 4, 5]
            sage: S.distance(3,S.sink())
            2
            sage: S.distance(1,S.sink())
            1
            sage: S._compare_vertices(1,3)
            -1
        """
        return self.distance(v, self._sink) - self.distance(w, self._sink)

    def _set_ring(self):
        r"""
        Set up polynomial ring for the sandpile.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_ring()
        """
        # first order the vertices according to their distance from the sink
        verts = self.vertices()
        verts = sorted(verts, self._compare_vertices)
        verts.reverse()

        # variable i refers to the i-th vertex in self.vertices()
        names = [self.vertices().index(v) for v in verts]

        vars = ''
        for i in names:
            vars += 'x' + str(i) + ','
        vars = vars[:-1]
        # create the ring
        self._ring = PolynomialRing(QQ, vars)
        # create the ideal
        gens = []
        for i in self.nonsink_vertices():
            new_gen = 'x' + str(self.vertices().index(i))
            new_gen += '^' + str(self.out_degree(i))
            new_gen += '-'
            for j in self._dict[i]:
                new_gen += 'x' + str(self.vertices().index(j))
                new_gen += '^' + str(self._dict[i][j]) + '*'
            new_gen = new_gen[:-1]
            gens.append(new_gen)
        self._unsaturated_ideal = self._ring.ideal(gens)

    def _set_ideal(self):
        r"""
        Create the saturated lattice ideal for the sandpile.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_ideal()
        """
        R = self.ring()
        I = self._unsaturated_ideal._singular_()
        self._ideal = R.ideal(I.sat(prod(R.gens())._singular_())[1])

    def unsaturated_ideal(self):
        r"""
        The unsaturated, homogeneous sandpile ideal.

        INPUT:

        None

        OUTPUT:

        ideal

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.unsaturated_ideal().gens()
            [x1^3 - x4*x3*x0, x2^3 - x5*x3*x0, x3^2 - x5*x2, x4^2 - x3*x1, x5^2 - x3*x2]
            sage: S.ideal().gens()
            [x2 - x0, x3^2 - x5*x0, x5*x3 - x0^2, x4^2 - x3*x1, x5^2 - x3*x0, x1^3 - x4*x3*x0, x4*x1^2 - x5*x0^2]
        """
        return self._unsaturated_ideal

    def ideal(self, gens=False):
        r"""
        The saturated, homogeneous sandpile ideal (or its generators if
        ``gens=True``).

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        ideal or, optionally, the generators of an ideal

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.ideal()
            Ideal (x2 - x0, x3^2 - x5*x0, x5*x3 - x0^2, x4^2 - x3*x1, x5^2 - x3*x0, x1^3 - x4*x3*x0, x4*x1^2 - x5*x0^2) of Multivariate Polynomial Ring in x5, x4, x3, x2, x1, x0 over Rational Field
            sage: S.ideal(True)
            [x2 - x0, x3^2 - x5*x0, x5*x3 - x0^2, x4^2 - x3*x1, x5^2 - x3*x0, x1^3 - x4*x3*x0, x4*x1^2 - x5*x0^2]
            sage: S.ideal().gens()  # another way to get the generators
            [x2 - x0, x3^2 - x5*x0, x5*x3 - x0^2, x4^2 - x3*x1, x5^2 - x3*x0, x1^3 - x4*x3*x0, x4*x1^2 - x5*x0^2]
        """
        if gens:
            return self._ideal.gens()
        else:
            return self._ideal

    def ring(self):
        r"""
        The ring containing the homogeneous sandpile ideal.

        INPUT:

        None

        OUTPUT:

        ring

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.ring()
            Multivariate Polynomial Ring in x5, x4, x3, x2, x1, x0 over Rational Field
            sage: S.ring().gens()
            (x5, x4, x3, x2, x1, x0)

        NOTES:

        The indeterminate `xi` corresponds to the `i`-th vertex as listed my
        the method ``vertices``. The term-ordering is degrevlex with
        indeterminates ordered according to their distance from the sink (larger
        indeterminates are further from the sink).
        """
        return self._ring

    def _set_resolution(self):
        r"""
        Compute the free resolution of the homogeneous sandpile ideal.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_resolution()
        """
        # get the resolution in singular form
        res = self.ideal()._singular_().mres(0)
        # compute the betti numbers
        #self._betti = [1] + [len(res[i])
        #        for i in range(1,len(res)-2)]
        self._betti = [1] + [len(x) for x in res]
        # convert the resolution to a list of Sage poly matrices
        result = []
        zero = self._ring.gens()[0]*0
        for i in range(1,len(res)+1):
            syz_mat = []
            new = [res[i][j] for j in range(1,res[i].size()+1)]
            for j in range(self._betti[i]):
                row = new[j].transpose().sage_matrix(self._ring)
                row = [r for r in row[0]]
                if len(row)<self._betti[i-1]:
                    row += [zero]*(self._betti[i-1]-len(row))
                syz_mat.append(row)
            syz_mat = matrix(self._ring, syz_mat).transpose()
            result.append(syz_mat)
        self._resolution = result
        self._singular_resolution = res

    def resolution(self, verbose=False):
        r"""
        This function computes a minimal free resolution of the homogeneous
        sandpile ideal.  If ``verbose`` is ``True``, then all of the mappings
        are returned.  Otherwise, the resolution is summarized.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        free resolution of the sandpile ideal

        EXAMPLES::

            sage: S = sandlib('gor')
            sage: S.resolution()
            'R^1 <-- R^5 <-- R^5 <-- R^1'
            sage: S.resolution(True)
            [
            [ x1^2 - x3*x0 x3*x1 - x2*x0  x3^2 - x2*x1  x2*x3 - x0^2  x2^2 - x1*x0],
            <BLANKLINE>
            [ x3  x2   0  x0   0]  [ x2^2 - x1*x0]
            [-x1 -x3  x2   0 -x0]  [-x2*x3 + x0^2]
            [ x0  x1   0  x2   0]  [-x3^2 + x2*x1]
            [  0   0 -x1 -x3  x2]  [x3*x1 - x2*x0]
            [  0   0  x0  x1 -x3], [ x1^2 - x3*x0]
            ]
            sage: r = S.resolution(True)
            sage: r[0]*r[1]
            [0 0 0 0 0]
            sage: r[1]*r[2]
            [0]
            [0]
            [0]
            [0]
            [0]
        """
        if verbose:
            return self._resolution
        else:
            r = ['R^'+str(i) for i in self._betti]
            return join(r,' <-- ')

    def _set_groebner(self):
        r"""
        Computes a Groebner basis for the homogeneous sandpile ideal with
        respect to the standard sandpile ordering (see ``ring``).

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S._set_groebner()
        """
        self._groebner = self._ideal.groebner_basis()

    def groebner(self):
        r"""
        A Groebner basis for the homogeneous sandpile ideal with
        respect to the standard sandpile ordering (see ``ring``).

        INPUT:

        None

        OUTPUT:

        Groebner basis

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.groebner()
            [x4*x1^2 - x5*x0^2, x1^3 - x4*x3*x0, x5^2 - x3*x0, x4^2 - x3*x1, x5*x3 - x0^2, x3^2 - x5*x0, x2 - x0]
        """
        return self._groebner

    def betti(self, verbose=True):
        r"""
        Computes the Betti table for the homogeneous sandpile ideal.  If
        ``verbose`` is ``True``, it prints the standard Betti table, otherwise,
        it returns a less formated table.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        Betti numbers for the sandpile


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.betti()
                       0     1     2     3     4     5
            ------------------------------------------
                0:     1     1     -     -     -     -
                1:     -     4     6     2     -     -
                2:     -     2     7     7     2     -
                3:     -     -     6    16    14     4
            ------------------------------------------
            total:     1     7    19    25    16     4
            sage: S.betti(False)
            [1, 7, 19, 25, 16, 4]
        """
        if verbose:
            print singular.eval('print(betti(%s),"betti")'%self._singular_resolution.name())
        else:
            return self._betti

    def solve(self):
        r"""
        Approximations of the complex affine zeros of the sandpile
        ideal.

        INPUT:

        None

        OUTPUT:

        list of complex numbers

        EXAMPLES::

            sage: S = Sandpile({0: {}, 1: {2: 2}, 2: {0: 4, 1: 1}}, 0)
            sage: S.solve()
            [[-0.707107 + 0.707107*I, 0.707107 - 0.707107*I], [-0.707107 - 0.707107*I, 0.707107 + 0.707107*I], [-I, -I], [I, I], [0.707107 + 0.707107*I, -0.707107 - 0.707107*I], [0.707107 - 0.707107*I, -0.707107 + 0.707107*I], [1, 1], [-1, -1]]
            sage: len(_)
            8
            sage: S.group_order()
            8

        NOTES:

        The solutions form a multiplicative group isomorphic to the sandpile
        group.  Generators for this group are given exactly by ``points()``.
        """
        singular.setring(self._ring._singular_())
        v = [singular.var(i) for i in range(1,singular.nvars(self._ring))]
        vars = '('
        for i in v:
            vars += str(i)
            vars += ','
        vars = vars[:-1]  # to get rid of the final ,
        vars += ')'
        L = singular.subst(self._ideal,
                singular.var(singular.nvars(self._ring)),1)
        R = singular.ring(0,vars,'lp')
        K = singular.fetch(self._ring,L)
        K = singular.groebner(K)
        singular.LIB('solve.lib')
        M = K.solve(5,1)
        singular.setring(M)
        sol= singular('SOL').sage_structured_str_list()
        sol = sol[0][0]
        sol = [map(eval,[j.replace('i','I') for j in k]) for k in sol]
        return sol

    def _set_points(self):
        r"""
        Generators for the multiplicative group of zeros of the sandpile
        ideal.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::
            sage: S = sandlib('generic')
            sage: S._set_points()
        """
        L = self._reduced_laplacian.transpose().dense_matrix()
        n = self.num_verts()-1;
        D, U, V = L.smith_form()
        self._points = []
        one = [1]*n
        for k in range(n):
            x = [exp(2*pi*I*U[k,t]/D[k,k]) for t in range(n)]
            if x not in self._points and x != one:
                self._points.append(x)

    def points(self):
        r"""
        Generators for the multiplicative group of zeros of the sandpile
        ideal.

        INPUT:

        None

        OUTPUT:

        list of complex numbers

        EXAMPLES:

        The sandpile group in this example is cyclic, and hence there is a
        single generator for the group of solutions.

        ::

            sage: S = sandlib('generic')
            sage: S.points()
            [[e^(4/5*I*pi), 1, e^(2/3*I*pi), e^(-34/15*I*pi), e^(-2/3*I*pi)]]
        """
        return self._points

    # FIX: use the is_symmetric functions for configurations.
    def symmetric_recurrents(self, orbits):
        r"""
        The list of symmetric recurrent configurations.

        INPUT:

        ``orbits`` - list of lists partitioning the vertices

        OUTPUT:

        list of recurrent configurations

        EXAMPLES::

            sage: S = sandlib('kite')
            sage: S.dict()
            {0: {},
             1: {0: 1, 2: 1, 3: 1},
             2: {1: 1, 3: 1, 4: 1},
             3: {1: 1, 2: 1, 4: 1},
             4: {2: 1, 3: 1}}
            sage: S.symmetric_recurrents([[1],[2,3],[4]])
            [{1: 2, 2: 2, 3: 2, 4: 1}, {1: 2, 2: 2, 3: 2, 4: 0}]
            sage: S.recurrents()
            [{1: 2, 2: 2, 3: 2, 4: 1},
             {1: 2, 2: 2, 3: 2, 4: 0},
             {1: 2, 2: 1, 3: 2, 4: 0},
             {1: 2, 2: 2, 3: 0, 4: 1},
             {1: 2, 2: 0, 3: 2, 4: 1},
             {1: 2, 2: 2, 3: 1, 4: 0},
             {1: 2, 2: 1, 3: 2, 4: 1},
             {1: 2, 2: 2, 3: 1, 4: 1}]

        NOTES:

        The user is responsible for ensuring that the list of orbits comes from
        a group of symmetries of the underlying graph.
        """
        sym_recurrents = []
        active = [self._max_stable]
        while active != []:
            c = active.pop()
            sym_recurrents.append(c)
            for orb in orbits:
                cnext = deepcopy(c)
                for v in orb:
                    cnext[v] += 1
                cnext = cnext.stabilize()
                if (cnext not in active) and (cnext not in sym_recurrents):
                    active.append(cnext)
        return deepcopy(sym_recurrents)

#######################################
########### SandpileConfig Class ##############
#######################################
class SandpileConfig(dict):
    r"""
    Class for configurations on a sandpile.
    """

    def __init__(self, S, c):
        r"""
        Create a configuration on a Sandpile.

        INPUT:

        - ``S`` - Sandpile
        - ``c`` - dict or list representing a configuration

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: C = SandpileConfig(S,[1,2,3,3,2])
            sage: C.order()
            3
            sage: ~(C+C+C)==S.identity()
            True
        """
        if len(c)==S.num_verts()-1:
            if isinstance(c, dict) or isinstance(c, SandpileConfig):
                dict.__init__(self,c)
            elif isinstance(c, list):
                c.reverse()
                config = {}
                for v in S.vertices():
                    if v!=S.sink():
                        config[v] = c.pop()
                dict.__init__(self,config)
        else:
            raise SyntaxError(c)

        self._sandpile = S
        self._vertices = S.nonsink_vertices()

    def __deepcopy__(self, memo):
        r"""
        Overrides the deepcopy method for dict.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: c = SandpileConfig(S, [1,2,3,4,5])
            sage: d = deepcopy(c)
            sage: d[1] += 10
            sage: c
            {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
            sage: d
            {1: 11, 2: 2, 3: 3, 4: 4, 5: 5}
        """
        c = SandpileConfig(self._sandpile, dict(self))
        c.__dict__.update(self.__dict__)
        return c

    def __setitem__(self, key, item):
        r"""
        Overrides the setitem method for dict.

        INPUT:

        - ``key``, ``item`` - objects

        OUTPUT:

        None

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [4,1])
            sage: c.equivalent_recurrent()
            {1: 1, 2: 1}
            sage: c.__dict__
            {'_equivalent_recurrent': [{1: 1, 2: 1}, {1: 2, 2: 1}],
             '_sandpile': Digraph on 3 vertices,
             '_vertices': [1, 2]}

        NOTES:

        In the example, above, changing the value of ``c`` at some vertex makes
        a call to setitem, which resets some of the stored variables for ``c``.
        """
        if key in self.keys():
            dict.__setitem__(self,key,item)
            S = self._sandpile
            V = self._vertices
            self.__dict__ = {'_sandpile':S, '_vertices': V}
        else:
            pass

    def __getattr__(self, name):
        """
        Set certain variables only when called.

        INPUT:

        ``name`` - name of an internal method

        OUTPUT:

        None.

        EXAMPLES::

            sage: S = complete_sandpile(4)
            sage: C = SandpileConfig(S,[1,1,1])
            sage: C.__getattr__('_deg')
            3
        """
        if name not in self.__dict__:
            if name=='_deg':
                self._set_deg()
                return self.__dict__[name]
            if name=='_stabilize':
                self._set_stabilize()
                return self.__dict__[name]
            if name=='_equivalent_recurrent':
                self._set_equivalent_recurrent()
                return self.__dict__[name]
            if name=='_is_recurrent':
                self._set_is_recurrent()
                return self.__dict__[name]
            if name=='_equivalent_superstable':
                self._set_equivalent_superstable()
                return self.__dict__[name]
            if name=='_is_superstable':
                self._set_is_superstable()
                return self.__dict__[name]
            else:
                raise AttributeError(name)

    def _set_deg(self):
        r"""
        Compute and store the degree of the configuration.

        INPUT:

        None

        OUTPUT:

        integer

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: c._set_deg()
        """
        self._deg = sum(self.values())

    def deg(self):
        r"""
        The degree of the configuration.

        INPUT:

        None

        OUTPUT:

        integer

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: c.deg()
            3
        """
        return self._deg

    def __add__(self,other):
        r"""
        Addition of configurations.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [3,2])
            sage: c + d
            {1: 4, 2: 4}
        """
        result = deepcopy(self)
        for v in self:
            result[v] += other[v]
        return result

    def __sub__(self,other):
        r"""
        Subtraction of configurations.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [3,2])
            sage: c - d
            {1: -2, 2: 0}
        """
        sum = deepcopy(self)
        for v in self:
            sum[v] -= other[v]
        return sum

    def __rsub__(self,other):
        r"""
        Right-side subtraction of configurations.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = {1: 3, 2: 2}
            sage: d - c
            {1: 2, 2: 0}

        TESTS::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = {1: 3, 2: 2}
            sage: c.__rsub__(d)
            {1: 2, 2: 0}
        """
        sum = deepcopy(other)
        for v in self:
            sum[v] -= self[v]
        return sum

    def __neg__(self):
        r"""
        The additive inverse of the configuration.

        INPUT:

        None

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: -c
            {1: -1, 2: -2}
        """
        return SandpileConfig(self._sandpile, [-self[v] for v in self._vertices])

    # recurrent addition
    def __mul__(self, other):
        r"""
        The recurrent element equivalent to the sum.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(4), 0)
            sage: c = SandpileConfig(S, [1,0,0])
            sage: c + c  # ordinary addition
            {1: 2, 2: 0, 3: 0}
            sage: c & c  # add and stabilize
            {1: 0, 2: 1, 3: 0}
            sage: c*c  # add and find equivalent recurrent
            {1: 1, 2: 1, 3: 1}
            sage: (c*c).is_recurrent()
            True
            sage: c*(-c) == S.identity()
            True
        """
        return (self+other).equivalent_recurrent()

    def __le__(self, other):
        r"""
        ``True`` if every component of ``self`` is at most that of
        ``other``.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [2,3])
            sage: e = SandpileConfig(S, [2,0])
            sage: c <= c
            True
            sage: c <= d
            True
            sage: d <= c
            False
            sage: c <= e
            False
            sage: e <= c
            False
        """
        return forall(self._vertices, lambda v: self[v]<=other[v])[0]

    def __lt__(self,other):
        r"""
        ``True`` if every component of ``self`` is at most that
        of ``other`` and the two configurations are not equal.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [2,3])
            sage: c < c
            False
            sage: c < d
            True
            sage: d < c
            False
        """
        return self<=other and self!=other

    def __ge__(self, other):
        r"""
        ``True`` if every component of ``self`` is at least that of
        ``other``.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [2,3])
            sage: e = SandpileConfig(S, [2,0])
            sage: c >= c
            True
            sage: d >= c
            True
            sage: c >= d
            False
            sage: e >= c
            False
            sage: c >= e
            False
        """
        return forall(self._vertices, lambda v: self[v]>=other[v])[0]

    def __gt__(self,other):
        r"""
        ``True`` if every component of ``self`` is at least that
        of ``other`` and the two configurations are not equal.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [1,3])
            sage: c > c
            False
            sage: d > c
            True
            sage: c > d
            False
        """
        return self>=other and self!=other

    # recurrent power
    def __pow__(self, k):
        r"""
        The recurrent element equivalent to the sum of the
        configuration with itself ``k`` times.  If ``k`` is negative, do the
        same for the negation of the configuration.  If ``k`` is zero, return
        the identity of the sandpile group.

        INPUT:

        ``k`` - SandpileConfig

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(4), 0)
            sage: c = SandpileConfig(S, [1,0,0])
            sage: c^3
            {1: 1, 2: 1, 3: 0}
            sage: (c + c + c) == c^3
            False
            sage: (c + c + c).equivalent_recurrent() == c^3
            True
            sage: c^(-1)
            {1: 1, 2: 1, 3: 0}
            sage: c^0 == S.identity()
            True
        """
        result = self._sandpile.zero_config()
        if k == 0:
            return self._sandpile.identity()
        else:
            if k<0:
                k = -k
                for i in range(k):
                    result -= self
            else:
                for i in range(k):
                    result += self
            return result.equivalent_recurrent()

    # stable addition
    def __and__(self, other):
        r"""
        The stabilization of the sum.

        INPUT:

        ``other`` - SandpileConfig

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(4), 0)
            sage: c = SandpileConfig(S, [1,0,0])
            sage: c + c  # ordinary addition
            {1: 2, 2: 0, 3: 0}
            sage: c & c  # add and stabilize
            {1: 0, 2: 1, 3: 0}
            sage: c*c  # add and find equivalent recurrent
            {1: 1, 2: 1, 3: 1}
            sage: ~(c + c) == c & c
            True
        """
        return ~(self+other)

    def sandpile(self):
        r"""
        The configuration's underlying sandpile.

        INPUT:

        None

        OUTPUT:

        Sandpile

        EXAMPLES::

            sage: S = sandlib('genus2')
            sage: c = S.identity()
            sage: c.sandpile()
            Digraph on 4 vertices
            sage: c.sandpile() == S
            True
        """
        return self._sandpile

    def values(self):
        r"""
        The values of the configuration as a list, sorted in the order
        of the vertices.

        INPUT:

        None

        OUTPUT:

        list of integers

        boolean

        EXAMPLES::

            sage: S = Sandpile({'a':[1,'b'], 'b':[1,'a'], 1:['a']},'a')
            sage: c = SandpileConfig(S, {'b':1, 1:2})
            sage: c
            {1: 2, 'b': 1}
            sage: c.values()
            [2, 1]
            sage: S.nonsink_vertices()
            [1, 'b']
        """
        return [self[v] for v in self._vertices]

    def dualize(self):
        r"""
        The difference between the maximal stable configuration and the
        configuration.

        INPUT:

        None

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: S.max_stable()
            {1: 1, 2: 1}
            sage: c.dualize()
            {1: 0, 2: -1}
            sage: S.max_stable() - c == c.dualize()
            True
        """
        return self._sandpile.max_stable()-self

    def fire_vertex(self, v):
        r"""
        Fire the vertex ``v``.

        INPUT:

        ``v`` - vertex

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: c = SandpileConfig(S, [1,2])
            sage: c.fire_vertex(2)
            {1: 2, 2: 0}
        """
        c = dict(self)
        c[v] -= self._sandpile.out_degree(v)
        for e in self._sandpile.outgoing_edges(v):
            if e[1]!=self._sandpile.sink():
                c[e[1]]+=e[2]
        return SandpileConfig(self._sandpile,c)

    def fire_script(self, sigma):
        r"""
        Fire the script ``sigma``, i.e., fire each vertex the indicated number
        of times.

        INPUT:

        ``sigma`` - SandpileConfig or (list or dict representing a SandpileConfig)

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(4), 0)
            sage: c = SandpileConfig(S, [1,2,3])
            sage: c.unstable()
            [2, 3]
            sage: c.fire_script(SandpileConfig(S,[0,1,1]))
            {1: 2, 2: 1, 3: 2}
            sage: c.fire_script(SandpileConfig(S,[2,0,0])) == c.fire_vertex(1).fire_vertex(1)
            True
        """
        c = dict(self)
        if not isinstance(sigma, SandpileConfig):
            sigma = SandpileConfig(self._sandpile, sigma)
        sigma = sigma.values()
        for i in range(len(sigma)):
            v = self._vertices[i]
            c[v] -= sigma[i]*self._sandpile.out_degree(v)
            for e in self._sandpile.outgoing_edges(v):
                if e[1]!=self._sandpile.sink():
                    c[e[1]]+=sigma[i]*e[2]
        return SandpileConfig(self._sandpile, c)

    def unstable(self):
        r"""
        List of the unstable vertices.

        INPUT:

        None

        OUTPUT:

        list of vertices

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(4), 0)
            sage: c = SandpileConfig(S, [1,2,3])
            sage: c.unstable()
            [2, 3]
        """
        return [v for v in self._vertices if
                self[v]>=self._sandpile.out_degree(v)]

    def fire_unstable(self):
        r"""
        Fire all unstable vertices.

        INPUT:

        None

        OUTPUT:

        SandpileConfig

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(4), 0)
            sage: c = SandpileConfig(S, [1,2,3])
            sage: c.fire_unstable()
            {1: 2, 2: 1, 3: 2}
        """
        c = dict(self)
        for v in self.unstable():
            c[v] -= self._sandpile.out_degree(v)
            for e in self._sandpile.outgoing_edges(v):
                if e[1]!=self._sandpile.sink():
                    c[e[1]]+=e[2]
        return SandpileConfig(self._sandpile,c)

    # FIX: make this faster by not converting to SandpileConfig?
    def _set_stabilize(self):
        r"""
        Computes the stabilized configuration and its firing vector.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: c = S.max_stable() + S.identity()
            sage: c._set_stabilize()
        """
        c = deepcopy(self)
        firing_vector = self._sandpile.zero_config()
        unstable = c.unstable()
        while unstable:
            for v in unstable:
                firing_vector[v] += 1
            c = c.fire_unstable()
            unstable = c.unstable()
        self._stabilize = [c, firing_vector]

    def stabilize(self, with_firing_vector=False):
        r"""
        The stabilized configuration. Optionally returns the
        corresponding firing vector.

        INPUT: s

        ``with_firing_vector`` (optional) -  boolean

        OUTPUT:

        ``SandpileConfig`` or ``[SandpileConfig, firing_vector]``

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: c = S.max_stable() + S.identity()
            sage: c.stabilize(True)
            [{1: 2, 2: 2, 3: 1, 4: 1, 5: 1}, {1: 1, 2: 5, 3: 7, 4: 1, 5: 6}]
            sage: S.max_stable() & S.identity()
            {1: 2, 2: 2, 3: 1, 4: 1, 5: 1}
            sage: S.max_stable() & S.identity() == c.stabilize()
            True
            sage: ~c
            {1: 2, 2: 2, 3: 1, 4: 1, 5: 1}
        """
        if with_firing_vector:
            return self._stabilize
        else:
            return self._stabilize[0]

    def __invert__(self):
        r"""
        The stabilized configuration.

        INPUT:

        None

        OUTPUT:

        ``SandpileConfig``

        Returns the stabilized configuration.
        EXAMPLES::

            sage: S = sandlib('generic')
            sage: c = S.max_stable() + S.identity()
            sage: ~c
            {1: 2, 2: 2, 3: 1, 4: 1, 5: 1}
            sage: ~c == c.stabilize()
            True
        """
        return self._stabilize[0]

    def support(self):
        r"""
        The input is a dictionary of integers.  The output is a list of keys
        of nonzero values of the dictionary.

        INPUT:

        None

        OUTPUT:

        list - support of the config

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: c = S.identity()
            sage: c.values()
            [2, 2, 1, 1, 0]
            sage: c.support()
            [1, 2, 3, 4]
            sage: S.vertices()
            [0, 1, 2, 3, 4, 5]
        """
        return [i for i in self.keys() if self[i] !=0]


    def add_random(self):
        r"""
        Add one grain of sand to a random nonsink vertex.

        INPUT:

        None

        OUTPUT:

        SandpileConfig

        EXAMPLES:

        We compute the 'sizes' of the avalanches caused by adding random grains
        of sand to the maximal stable configuration on a grid graph.  The
        function ``stabilize()`` returns the firing vector of the
        stabilization, a dictionary whose values say how many times each vertex
        fires in the stabilization.

        ::

            sage: S = grid_sandpile(10,10)
            sage: m = S.max_stable()
            sage: a = []
            sage: for i in range(1000):
            ...       m = m.add_random()
            ...       m, f = m.stabilize(True)
            ...       a.append(sum(f.values()))
            ...
            sage: p = list_plot([[log(i+1),log(a.count(i))] for i in [0..max(a)] if a.count(i)])
            sage: p.axes_labels(['log(N)','log(D(N))'])
            sage: t = text("Distribution of avalanche sizes", (2,2), rgbcolor=(1,0,0))
            sage: show(p+t,axes_labels=['log(N)','log(D(N))'])
        """
        config = dict(self)
        C = CombinatorialClass()
        C.list = lambda: self.keys()
        config[C.random_element()] += 1
        return SandpileConfig(self._sandpile,config)

    def order(self):
        r"""
        The order of the recurrent element equivalent to ``config``.

        INPUT:

        ``config`` - configuration

        OUTPUT:

        integer

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: [r.order() for r in S.recurrents()]
            [3, 3, 5, 15, 15, 15, 5, 15, 15, 5, 15, 5, 15, 1, 15]
        """
        v = vector(self.values())
        w = v*self._sandpile.reduced_laplacian().dense_matrix()**(-1)
        return lcm([denominator(i) for i in w])

    def is_stable(self):
        r"""
        ``True`` if stable.

        INPUT:

        None

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.max_stable().is_stable()
            True
            sage: (S.max_stable() + S.max_stable()).is_stable()
            False
            sage: (S.max_stable() & S.max_stable()).is_stable()
            True
        """
        for v in self._vertices:
            if self[v] >= self._sandpile.out_degree(v):
                return False
        return True

    def _set_equivalent_recurrent(self):
        r"""
        Sets the equivalent recurrent configuration and the corresponding
        firing vector.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: a = -S.max_stable()
            sage: a._set_equivalent_recurrent()
        """
        old = self
        # revise next line with c._sandpile.zero_config()
        firing_vector = self._sandpile.zero_config()
        done = False
        bs = self._sandpile.burning_script()
        bc = self._sandpile.burning_config()
        while not done:
            firing_vector = firing_vector - bs
            new, new_fire = (old + bc).stabilize(True)
            firing_vector = firing_vector + new_fire
            if new == old:
                done = True
            else:
                old = new
        self._equivalent_recurrent = [new, firing_vector]

    def equivalent_recurrent(self, with_firing_vector=False):
        r"""
        The recurrent configuration equivalent to the given
        configuration and optionally returns the corresponding firing vector.

        INPUT:

        ``with_firing_vector`` (optional) -  boolean

        OUTPUT:

        ``SandpileConfig`` or ``[SandpileConfig, firing_vector]``


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: c = SandpileConfig(S, [0,0,0,0,0])
            sage: c.equivalent_recurrent() == S.identity()
            True
            sage: x = c.equivalent_recurrent(True)
            sage: r = vector([x[0][v] for v in S.nonsink_vertices()])
            sage: f = vector([x[1][v] for v in S.nonsink_vertices()])
            sage: cv = vector(c.values())
            sage: r == cv - f*S.reduced_laplacian()
            True

        NOTES:

        Let `L` be the reduced laplacian, `c` the initial configuration, `r` the
        returned configuration, and `f` the firing vector.  Then `r = c - f\cdot
        L`.
        """
        if with_firing_vector:
            return self._equivalent_recurrent
        else:
            return self._equivalent_recurrent[0]

    def _set_is_recurrent(self):
        r"""
        Computes and stores whether the configuration is recurrent.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.max_stable()._set_is_recurrent()
        """
        if '_recurrents' in self._sandpile.__dict__:
            self._is_recurrent = (self in self._sandpile._recurrents)
        elif '_equivalent_recurrent' in self.__dict__:
            self._is_recurrent = (self._equivalent_recurrent == self)
        else:
            # add the burning configuration to config
            # change once sandpile._burning_config is a SandpileConfig
            b = SandpileConfig(self._sandpile,self._sandpile._burning_config)
            c = ~(self + b)
            self._is_recurrent = (c == self)

    def is_recurrent(self):
        r"""
        ``True`` if the configuration is recurrent.

        INPUT:

        None

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.identity().is_recurrent()
            True
            sage: S.zero_config().is_recurrent()
            False
        """
        return self._is_recurrent

    def _set_equivalent_superstable(self):
        r"""
        Sets the superstable configuration equivalent to the given
        configuration and its corresponding firing vector.

        INPUT:

        None

        OUTPUT:

        ``[configuration, firing_vector]``


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: m = S.max_stable()
            sage: m._set_equivalent_superstable()
        """
        r, fv = self.dualize().equivalent_recurrent(with_firing_vector=True)
        self._equivalent_superstable = [r.dualize(), -fv]

    def equivalent_superstable(self, with_firing_vector=False):
        r"""
        The equivalent superstable configuration. Optionally
        returns the corresponding firing vector.

        INPUT:

        ``with_firing_vector`` (optional) - boolean

        OUTPUT:

        ``SandpileConfig`` or ``[SandpileConfig, firing_vector]``


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: m = S.max_stable()
            sage: m.equivalent_superstable().is_superstable()
            True
            sage: x = m.equivalent_superstable(True)
            sage: s = vector(x[0].values())
            sage: f = vector(x[1].values())
            sage: mv = vector(m.values())
            sage: s == mv - f*S.reduced_laplacian()
            True

        NOTES:

        Let `L` be the reduced laplacian, `c` the initial configuration, `s` the
        returned configuration, and `f` the firing vector.  Then `s = c - f\cdot
        L`.
        """
        if with_firing_vector:
            return self._equivalent_superstable
        else:
            return self._equivalent_superstable[0]

    def _set_is_superstable(self):
        r"""
        Computes and stores whether ``config`` is superstable.

        INPUT:

        None

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.zero_config()._set_is_superstable()
        """
        if '_superstables' in self._sandpile.__dict__:
            self._is_superstable = (self in self._sandpile._superstables)
        elif '_equivalent_superstable' in self.__dict__:
            self._is_superstable = (self._equivalent_superstable[0] == self)
        else:
            self._is_superstable = self.dualize().is_recurrent()

    def is_superstable(self):
        r"""
        ``True`` if ``config`` is superstable, i.e., whether its dual is
        recurrent.

        INPUT:

        None

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.zero_config().is_superstable()
            True
        """
        return self._is_superstable

    def is_symmetric(self, orbits):
        r"""
        This function checks if the values of the configuration are constant
        over the vertices in each sublist of ``orbits``.

        INPUT:

         ``orbits`` - list of lists of vertices

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = sandlib('kite')
            sage: S.dict()
            {0: {},
             1: {0: 1, 2: 1, 3: 1},
             2: {1: 1, 3: 1, 4: 1},
             3: {1: 1, 2: 1, 4: 1},
             4: {2: 1, 3: 1}}
            sage: c = SandpileConfig(S, [1, 2, 2, 3])
            sage: c.is_symmetric([[2,3]])
            True
        """
        for x in orbits:
            if len(set([self[v] for v in x])) > 1:
                return False
        return True

    def show(self,sink=True,colors=True,heights=False,directed=None,**kwds):
        r"""
        Show the configuration.

        INPUT:

        - ``sink`` - whether to show the sink
        - ``colors`` - whether to color-code the amount of sand on each vertex
        - ``heights`` - whether to label each vertex with the amount of sand
        - ``kwds`` - arguments passed to the show method for Graph
        - ``directed`` - whether to draw directed edges

        OUTPUT:

        None

        EXAMPLES::

            sage: S=sandlib('genus2')
            sage: c=S.identity()
            sage: S=sandlib('genus2')
            sage: c=S.identity()
            sage: c.show()
            sage: c.show(directed=False)
            sage: c.show(sink=False,colors=False,heights=True)
        """
        if directed==True:
            T = DiGraph(self.sandpile())
        elif directed==False:
            T = Graph(self.sandpile())
        elif self.sandpile().is_directed():
            T = DiGraph(self.sandpile())
        else:
            T = Graph(self.sandpile())

        max_height = max(self.sandpile().out_degree_sequence())
        if not sink:
            T.delete_vertex(self.sandpile().sink())
        if heights:
            a = {}
            for i in T.vertices():
                if i==self.sandpile().sink():
                    a[i] = str(i)
                else:
                    a[i] = str(i)+":"+str(self[i])
            T.relabel(a)
        if colors:
            vc = {}  # vertex colors
            r = rainbow(max_height) # colors
            for i in range(max_height):
                vc[r[i]] = []
            for i in self.sandpile().nonsink_vertices():
                if heights:
                    vc[r[self[i]]].append(a[i])
                else:
                    vc[r[self[i]]].append(i)
            T.show(vertex_colors=vc,**kwds)
        else:
            T.show(**kwds)

#######################################
########### SandpileDivisor Class #############
#######################################

class SandpileDivisor(dict):
    r"""
    Class for divisors on a sandpile.
    """

    def __init__(self, S, D):
        r"""
        Create a divisor on a Sandpile.

        INPUT:

        - ``S`` - Sandpile
        - ``D`` - dict or list representing a divisor

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S,[0,1,0,1,1,3])
            sage: D.support()
            [1, 3, 4, 5]

        """
        if len(D)==S.num_verts():
            if type(D) in [dict, SandpileDivisor, SandpileConfig]:
                dict.__init__(self,dict(D))
            elif isinstance(D, list):
                div = {}
                for i in range(S.num_verts()):
                    div[S.vertices()[i]] = D[i]
                    dict.__init__(self,div)
        else:
            raise SyntaxError(D)

        self._sandpile = S
        self._vertices = S.vertices()

    def __deepcopy__(self, memo):
        r"""
        Overrides the deepcopy method for dict.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S, [1,2,3,4,5,6])
            sage: E = deepcopy(D)
            sage: E[0] += 10
            sage: D
            {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6}
            sage: E
            {0: 11, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6}
        """
        D = SandpileDivisor(self._sandpile, dict(self))
        D.__dict__.update(self.__dict__)
        return D

    def __setitem__(self, key, item):
        r"""
        Overrides the setitem method for dict.

        INPUT:

        - ``key``, ``item`` - objects

        OUTPUT:

        None

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [0,1,1])
            sage: eff = D.effective_div() # optional - 4ti2
            sage: D.__dict__ # optional - 4ti2
            {'_effective_div': [{0: 2, 1: 0, 2: 0}, {0: 0, 1: 1, 2: 1}],
             '_linear_system': {'homog': [[1, 1, 1], [-1, -1, -1]],
              'inhomog': [[1, 0, 0], [0, 0, 0], [0, -1, -1]],
              'num_homog': 2,
              'num_inhomog': 3},
             '_sandpile': Digraph on 3 vertices,
             '_vertices': [0, 1, 2]}
            sage: D[0] += 1 # optional - 4ti2
            sage: D.__dict__ # optional - 4ti2
            {'_sandpile': Digraph on 3 vertices, '_vertices': [0, 1, 2]}

        NOTES:

        In the example, above, changing the value of ``D`` at some vertex makes
        a call to setitem, which resets some of the stored variables for ``D``.
        """
        if key in self.keys():
            dict.__setitem__(self,key,item)
            S = self._sandpile
            V = self._vertices
            self.__dict__ = {'_sandpile':S, '_vertices': V}
        else:
            pass

    def __getattr__(self, name):
        """
        Set certain variables only when called.

        INPUT:

        ``name`` - name of an internal method

        OUTPUT:

        None.

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S,[0,1,0,1,1,3])
            sage: D.__getattr__('_deg')
            6
        """
        if name not in self.__dict__:
            if name=='_deg':
                self._set_deg()
                return self.__dict__[name]
            if name=='_linear_system':
                self._set_linear_system()
                return self.__dict__[name]
            if name=='_effective_div':
                self._set_effective_div()
                return self.__dict__[name]
            if name=='_r_of_D':
                self._set_r_of_D()
                return self.__dict__[name]
            if name=='_Dcomplex':
                self._set_Dcomplex()
                return self.__dict__[name]
            if name=='_life':
                self._set_life()
                return self.__dict__[name]
            else:
                raise AttributeError(name)

    def _set_deg(self):
        r"""
        Compute and store the degree of the divisor.

        INPUT:

        None

        OUTPUT:

        integer

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D._set_deg()
        """
        self._deg = sum(self.values())

    def deg(self):
        r"""
        The degree of the divisor.

        INPUT:

        None

        OUTPUT:

        integer

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.deg()
            6
        """
        return self._deg

    def __add__(self, other):
        r"""
        Addition of divisors.

        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [3,2,1])
            sage: D + E
            {0: 4, 1: 4, 2: 4}
        """
        sum = deepcopy(self)
        for v in self:
            sum[v] += other[v]
        return sum

    def __radd__(self,other):
        r"""
        Right-side addition of divisors.

        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [3,2,1])
            sage: D.__radd__(E)
            {0: 4, 1: 4, 2: 4}
        """
        sum = deepcopy(other)
        for v in self:
            sum[v] += self[v]
        return sum

    def __sub__(self, other):
        r"""
        Subtraction of divisors.

        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        Difference of ``self`` and ``other``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [3,2,1])
            sage: D - E
            {0: -2, 1: 0, 2: 2}
        """
        sum = deepcopy(self)
        for v in self:
            sum[v] -= other[v]
        return sum

    def __rsub__(self, other):
        r"""
        Right-side subtraction of divisors.
        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        Difference of ``self`` and ``other``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = {0: 3, 1: 2, 2: 1}
            sage: D.__rsub__(E)
            {0: 2, 1: 0, 2: -2}
            sage: E - D
            {0: 2, 1: 0, 2: -2}
        """
        sum = deepcopy(other)
        for v in self:
            sum[v] -= self[v]
        return sum

    def __neg__(self):
        r"""
        The additive inverse of the divisor.

        INPUT:

        None

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: -D
            {0: -1, 1: -2, 2: -3}
        """
        return SandpileDivisor(self._sandpile, [-self[v] for v in self._vertices])

    def __le__(self, other):
        r"""
        ``True`` if every component of ``self`` is at most that of
        ``other``.

        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [2,3,4])
            sage: F = SandpileDivisor(S, [2,0,4])
            sage: D <= D
            True
            sage: D <= E
            True
            sage: E <= D
            False
            sage: D <= F
            False
            sage: F <= D
            False
        """
        return forall(self._vertices, lambda v: self[v]<=other[v])[0]

    def __lt__(self,other):
        r"""
        ``True`` if every component of ``self`` is at most that
        of ``other`` and the two divisors are not equal.

        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [2,3,4])
            sage: D < D
            False
            sage: D < E
            True
            sage: E < D
            False
        """
        return self<=other and self!=other

    def __ge__(self, other):
        r"""
        ``True`` if every component of ``self`` is at least that of
        ``other``.

        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [2,3,4])
            sage: F = SandpileDivisor(S, [2,0,4])
            sage: D >= D
            True
            sage: E >= D
            True
            sage: D >= E
            False
            sage: F >= D
            False
            sage: D >= F
            False
        """
        return forall(self._vertices, lambda v: self[v]>=other[v])[0]

    def __gt__(self,other):
        r"""
        ``True`` if every component of ``self`` is at least that
        of ``other`` and the two divisors are not equal.

        INPUT:

        ``other`` - SandpileDivisor

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [1,3,4])
            sage: D > D
            False
            sage: E > D
            True
            sage: D > E
            False
        """
        return self>=other and self!=other

    def sandpile(self):
        r"""
        The divisor's underlying sandpile.

        INPUT:

        None

        OUTPUT:

        Sandpile

        EXAMPLES::

            sage: S = sandlib('genus2')
            sage: D = SandpileDivisor(S,[1,-2,0,3])
            sage: D.sandpile()
            Digraph on 4 vertices
            sage: D.sandpile() == S
            True
        """
        return self._sandpile

    def values(self):
        r"""
        The values of the divisor as a list, sorted in the order of the
        vertices.

        INPUT:

        None

        OUTPUT:

        list of integers

        boolean

        EXAMPLES::

            sage: S = Sandpile({'a':[1,'b'], 'b':[1,'a'], 1:['a']},'a')
            sage: D = SandpileDivisor(S, {'a':0, 'b':1, 1:2})
            sage: D
            {'a': 0, 1: 2, 'b': 1}
            sage: D.values()
            [2, 0, 1]
            sage: S.vertices()
            [1, 'a', 'b']
        """
        return [self[v] for v in self._vertices]

    def dualize(self):
        r"""
        The difference between the maximal stable divisor and the
        divisor.

        INPUT:

        None

        OUTPUT:

        SandpileDivisor

        EXAMPLES::
            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.dualize()
            {0: 0, 1: -1, 2: -2}
            sage: S.max_stable_div() - D == D.dualize()
            True
        """
        return self._sandpile.max_stable_div() - self

    def fire_vertex(self,v):
        r"""
        Fire the vertex ``v``.

        INPUT:

        ``v`` - vertex

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.fire_vertex(1)
            {0: 2, 1: 0, 2: 4}
        """
        D = dict(self)
        D[v] -= self._sandpile.out_degree(v)
        for e in self._sandpile.outgoing_edges(v):
            D[e[1]]+=e[2]
        return SandpileDivisor(self._sandpile,D)

    def fire_script(self, sigma):
        r"""
        Fire the script ``sigma``, i.e., fire each vertex the indicated number
        of times.

        INPUT:

        ``sigma`` - SandpileDivisor or (list or dict representing a SandpileDivisor)

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.unstable()
            [1, 2]
            sage: D.fire_script([0,1,1])
            {0: 3, 1: 1, 2: 2}
            sage: D.fire_script(SandpileDivisor(S,[2,0,0])) == D.fire_vertex(0).fire_vertex(0)
            True
        """
        D = dict(self)
        if not isinstance(sigma, SandpileDivisor):
            sigma = SandpileDivisor(self._sandpile, sigma)
        sigma = sigma.values()
        for i in range(len(sigma)):
            v = self._vertices[i]
            D[v] -= sigma[i]*self._sandpile.out_degree(v)
            for e in self._sandpile.outgoing_edges(v):
                D[e[1]]+=sigma[i]*e[2]
        return SandpileDivisor(self._sandpile, D)

    def unstable(self):
        r"""
        List of the unstable vertices.

        INPUT:

        None

        OUTPUT:

        list of vertices

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.unstable()
            [1, 2]
        """
        return [v for v in self._vertices if
                self[v]>=self._sandpile.out_degree(v)]

    def fire_unstable(self):
        r"""
        Fire all unstable vertices.

        INPUT:

        None

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.fire_unstable()
            {0: 3, 1: 1, 2: 2}
        """
        D = dict(self)
        for v in self.unstable():
            D[v] -= self._sandpile.out_degree(v)
            for e in self._sandpile.outgoing_edges(v):
                D[e[1]]+=e[2]
        return SandpileDivisor(self._sandpile,D)

    def _set_linear_system(self):
        r"""
        Computes and stores the complete linear system of a divisor.

        INPUT: None

        OUTPUT:

        dict - ``{num_homog: int, homog:list, num_inhomog:int, inhomog:list}``

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [0,1,1])
            sage: D._set_linear_system() # optional - 4ti2

        WARNING:

        This method requires 4ti2.
        """
        # import os

        L = self._sandpile._laplacian.transpose()
        n = self._sandpile.num_verts()

        # temporary file names
        lin_sys = tmp_filename()
        lin_sys_mat = lin_sys + '.mat'
        lin_sys_rel = lin_sys + '.rel'
        lin_sys_rhs = lin_sys + '.rhs'
        lin_sys_sign= lin_sys + '.sign'
        lin_sys_zhom= lin_sys + '.zhom'
        lin_sys_zinhom= lin_sys + '.zinhom'
        lin_sys_log = lin_sys + '.log'

        mat_file = open(lin_sys_mat,'w')
        mat_file.write(str(n)+' ')
        mat_file.write(str(n)+'\n')
        for r in L:
            mat_file.write(join(map(str,r)))
            mat_file.write('\n')
        mat_file.close()
        # relations file
        rel_file = open(lin_sys_rel,'w')
        rel_file.write('1 ')
        rel_file.write(str(n)+'\n')
        rel_file.write(join(['>']*n))
        rel_file.write('\n')
        rel_file.close()
        # right-hand side file
        rhs_file = open(lin_sys_rhs,'w')
        rhs_file.write('1 ')
        rhs_file.write(str(n)+'\n')
        rhs_file.write(join([str(-i) for i in self.values()]))
        rhs_file.write('\n')
        rhs_file.close()
        # sign file
        sign_file = open(lin_sys_sign,'w')
        sign_file.write('1 ')
        sign_file.write(str(n)+'\n')
        """
        Conjecture: taking only 1s just below is OK, i.e., looking for solutions
        with nonnegative entries.  The Laplacian has kernel of dimension 1,
        generated by a nonnegative vector.  I would like to say that translating
        by this vector, we transform any solution into a nonnegative solution.
        What if the vector in the kernel does not have full support though?
        """
        sign_file.write(join(['2']*n))  # so maybe a 1 could go here
        sign_file.write('\n')
        sign_file.close()
        # compute
        try:
            os.system(path_to_zsolve+' -q ' + lin_sys + ' > ' + lin_sys_log)
            # process the results
            zhom_file = open(lin_sys_zhom,'r')
        except IOError:
            print """
                 **********************************
                 *** This method requires 4ti2. ***
                 **********************************

            Read the beginning of sandpile.sage or see the Sage Sandpiles
            documentation for installation instructions.
            """
            return
        ## first, the cone generators (the homogeneous points)
        a = zhom_file.read()
        zhom_file.close()
        a = a.split('\n')
        # a starts with two numbers. We are interested in the first one
        num_homog = int(a[0].split()[0])
        homog = [map(int,i.split()) for i in a[1:-1]]
        ## second, the inhomogeneous points
        zinhom_file = open(lin_sys_zinhom,'r')
        b = zinhom_file.read()
        zinhom_file.close()
        b = b.split('\n')
        num_inhomog = int(b[0].split()[0])
        inhomog = [map(int,i.split()) for i in b[1:-1]]
        self._linear_system = {'num_homog':num_homog, 'homog':homog,
                'num_inhomog':num_inhomog, 'inhomog':inhomog}

    def linear_system(self):
        r"""
        The complete linear system of a divisor.

        INPUT: None

        OUTPUT:

        dict - ``{num_homog: int, homog:list, num_inhomog:int, inhomog:list}``

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S, [0,0,0,0,0,2])
            sage: D.linear_system() # optional - 4ti2
            {'homog': [[1, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0, 0]],
             'inhomog': [[0, 0, 0, 0, 0, -1], [0, 0, -1, -1, 0, -2], [0, 0, 0, 0, 0, 0]],
             'num_homog': 2,
             'num_inhomog': 3}

        NOTES:

        If `L` is the Laplacian, an arbitrary `v` such that `v\cdot L\geq -D`
        has the form `v = w + t` where `w` is in ``inhomg`` and `t` is in the
        integer span of ``homog`` in the output of ``linear_system(D)``.

        WARNING:

        This method requires 4ti2.
        """
        return self._linear_system

    def _set_effective_div(self):
        r"""
        Computes all of the linearly equivalent effective divisors linearly.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S, [0,0,0,0,0,2]) # optional - 4ti2
            sage: D._set_effective_div() # optional - 4ti2
        """
        result = []
        r = self.linear_system()
        d = vector(self.values())
        for x in r['inhomog']:
            c = vector(x)*self._sandpile._laplacian + d
            c = SandpileDivisor(self._sandpile,list(c))
            if c not in result:
                result.append(c)
        self._effective_div = result

    def effective_div(self, verbose=True):
        r"""
        All linearly equivalent effective divisors.  If ``verbose``
        is ``False``, the divisors are converted to lists of integers.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        list (of divisors)

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S, [0,0,0,0,0,2]) # optional - 4ti2
            sage: D.effective_div() # optional - 4ti2
            [{0: 0, 1: 0, 2: 1, 3: 1, 4: 0, 5: 0}, {0: 1, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0}, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 2}]
            sage: D.effective_div(False) # optional - 4ti2
            [[0, 0, 1, 1, 0, 0], [1, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 2]]
        """
        if verbose:
            return deepcopy(self._effective_div)
        else:
            return [E.values() for E in self._effective_div]

    def _set_r_of_D(self, verbose=False):
        r"""
        Computes ``r(D)`` and an effective divisor ``F`` such
        that ``|D - F|`` is empty.

        INPUT:

        verbose (optional) - boolean

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S, [0,0,0,0,0,4]) # optional - 4ti2
            sage: D._set_r_of_D() # optional - 4ti2
        """
        eff = self.effective_div()
        n = self._sandpile.num_verts()
        r = -1
        if eff == []:
            self._r_of_D = (r, self)
            return
        else:
            d = vector(self.values())
            # standard basis vectors
            e = []
            for i in range(n):
                v = vector([0]*n)
                v[i] += 1
                e.append(v)
            level = [vector([0]*n)]
            while True:
                r += 1
                if verbose:
                    print r
                new_level = []
                for v in level:
                    for i in range(n):
                        w = v + e[i]
                        if w not in new_level:
                            new_level.append(w)
                            C = d - w
                            C = SandpileDivisor(self._sandpile,list(C))
                            eff = C.effective_div()
                            if eff == []:
                                self._r_of_D = (r, SandpileDivisor(self._sandpile,list(w)))
                                return
                level = new_level

    def r_of_D(self, verbose=False):
        r"""
        Returns ``r(D)`` and, if ``verbose`` is ``True``, an effective divisor
        ``F`` such that ``|D - F|`` is empty.

        INPUT:

        ``verbose`` (optional) - boolean

        OUTPUT:

        integer ``r(D)`` or tuple (integer ``r(D)``, divisor ``F``)

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S, [0,0,0,0,0,4]) # optional - 4ti2
            sage: E = D.r_of_D(True) # optional - 4ti2
            sage: E # optional - 4ti2
            (1, {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: 0})
            sage: F = E[1] # optional - 4ti2
            sage: (D - F).values() # optional - 4ti2
            [0, -1, 0, -1, 0, 4]
            sage: (D - F).effective_div() # optional - 4ti2
            []
            sage: SandpileDivisor(S, [0,0,0,0,0,-4]).r_of_D(True) # optional - 4ti2
            (-1, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: -4})
        """
        if verbose:
            return self._r_of_D
        else:
            return self._r_of_D[0]

    def support(self):
        r"""
        List of keys of the nonzero values of the divisor.

        INPUT:

        None

        OUTPUT:

        list - support of the divisor


        EXAMPLES::

            sage: S = sandlib('generic')
            sage: c = S.identity()
            sage: c.values()
            [2, 2, 1, 1, 0]
            sage: c.support()
            [1, 2, 3, 4]
            sage: S.vertices()
            [0, 1, 2, 3, 4, 5]
        """
        return [i for i in self.keys() if self[i] !=0]

    def _set_Dcomplex(self):
        r"""
        Computes the simplicial complex determined by the supports of the
        linearly equivalent effective divisors.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: D = SandpileDivisor(S, [0,1,2,0,0,1]) # optional - 4ti2
            sage: D._set_Dcomplex() # optional - 4ti2
        """
        simp = []
        for E in self.effective_div():
            supp_E = E.support()
            test = True
            for s in simp:
                if set(supp_E).issubset(set(s)):
                    test = False
                    break
            if test:
                simp.append(supp_E)
        result = []
        simp.reverse()
        while simp != []:
            supp = simp.pop()
            test = True
            for s in simp:
                if set(supp).issubset(set(s)):
                    test = False
                    break
            if test:
                result.append(supp)
        self._Dcomplex = SimplicialComplex(result)

    def Dcomplex(self):
        r"""
        The simplicial complex determined by the supports of the
        linearly equivalent effective divisors.

        INPUT:

        None

        OUTPUT:

        simplicial complex

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: p = SandpileDivisor(S, [0,1,2,0,0,1]).Dcomplex() # optional - 4ti2
            sage: p.homology() # optional - 4ti2
            {0: 0, 1: Z x Z, 2: 0, 3: 0}
            sage: p.f_vector() # optional - 4ti2
            [1, 6, 15, 9, 1]
            sage: p.betti() # optional - 4ti2
            {0: 1, 1: 2, 2: 0, 3: 0}
        """
        return self._Dcomplex

    def betti(self):
        r"""
        The Betti numbers for the simplicial complex associated with
        the divisor.

        INPUT:

        None

        OUTPUT:

        dictionary of integers

        EXAMPLES::

            sage: S = Sandpile(graphs.CycleGraph(3), 0)
            sage: D = SandpileDivisor(S, [2,0,1])
            sage: D.betti() # optional - 4ti2
            {0: 1, 1: 1}
        """
        return self.Dcomplex().betti()

    def add_random(self):
        r"""
        Add one grain of sand to a random vertex.

        INPUT:

        None

        OUTPUT:

        SandpileDivisor

        EXAMPLES::

            sage: S = sandlib('generic')
            sage: S.zero_div().add_random()  #random
            {0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0}
        """
        D = dict(self)
        C = CombinatorialClass()
        C.list = lambda: self.keys()
        D[C.random_element()] += 1
        return SandpileDivisor(self._sandpile,D)

    def is_symmetric(self, orbits):
        r"""
        This function checks if the values of the divisor are constant
        over the vertices in each sublist of ``orbits``.

        INPUT:

         - ``orbits`` - list of lists of vertices

        OUTPUT:

        boolean

        EXAMPLES::

            sage: S = sandlib('kite')
            sage: S.dict()
            {0: {},
             1: {0: 1, 2: 1, 3: 1},
             2: {1: 1, 3: 1, 4: 1},
             3: {1: 1, 2: 1, 4: 1},
             4: {2: 1, 3: 1}}
            sage: D = SandpileDivisor(S, [2,1, 2, 2, 3])
            sage: D.is_symmetric([[0,2,3]])
            True
        """
        for x in orbits:
            if len(set([self[v] for v in x])) > 1:
                return False
        return True

    def _set_life(self):
        r"""
        Will the sequence of divisors ``D_i`` where ``D_{i+1}`` is obtained from
        ``D_i`` by firing all unstable vertices of ``D_i`` stabilize?  If so,
        save the resulting cycle, otherwise save ``[]``.

        INPUT:

        None

        OUTPUT:

        None

        EXAMPLES::

            sage: S = complete_sandpile(4)
            sage: D = SandpileDivisor(S, {0: 4, 1: 3, 2: 3, 3: 2})
            sage: D._set_life()
        """
        oldD = deepcopy(self)
        result = [oldD]
        while True:
            if oldD.unstable()==[]:
                self._life = []
                return
            newD = oldD.fire_unstable()
            if newD not in result:
                result.append(newD)
                oldD = deepcopy(newD)
            else:
                self._life = result[result.index(newD):]
                return

    def is_alive(self, cycle=False):
        r"""
        Will the divisor stabilize under repeated firings of all unstable
        vertices?  Optionally returns the resulting cycle.

        INPUT:

        ``cycle`` (optional) - boolean

        OUTPUT:

        boolean or optionally, a list of SandpileDivisors

        EXAMPLES::

            sage: S = complete_sandpile(4)
            sage: D = SandpileDivisor(S, {0: 4, 1: 3, 2: 3, 3: 2})
            sage: D.is_alive()
            True
            sage: D.is_alive(True)
            [{0: 4, 1: 3, 2: 3, 3: 2}, {0: 3, 1: 2, 2: 2, 3: 5}, {0: 1, 1: 4, 2: 4, 3: 3}]
        """
        if cycle:
            return self._life
        else:
            return self._life != []

    def show(self,heights=True,directed=None,**kwds):
        r"""
        Show the divisor.

        INPUT:

        - ``heights`` - whether to label each vertex with the amount of sand
        - ``kwds`` - arguments passed to the show method for Graph
        - ``directed`` - whether to draw directed edges

        OUTPUT:

        None

        EXAMPLES::

            sage: S = sandlib('genus2')
            sage: D = SandpileDivisor(S,[1,-2,0,2])
            sage: D.show(graph_border=True,vertex_size=700,directed=False)
        """
        if directed==True:
            T = DiGraph(self.sandpile())
        elif directed==False:
            T = Graph(self.sandpile())
        elif self.sandpile().is_directed():
            T = DiGraph(self.sandpile())
        else:
            T = Graph(self.sandpile())

        max_height = max(self.sandpile().out_degree_sequence())
        if heights:
            a = {}
            for i in T.vertices():
                a[i] = str(i)+":"+str(T[i])
            T.relabel(a)
        T.show(**kwds)

#######################################
######### Some test graphs ############
#######################################

def sandlib(selector=None):
    r"""
    Returns the sandpile identified by ``selector``.  If no argument is
    given, a description of the sandpiles in the sandlib is printed.

    INPUT:

    ``selector`` - identifier or None

    OUTPUT:

    sandpile or description

    EXAMPLES::

            sage: sandlib()
            <BLANKLINE>
              Sandpiles in the sandlib:
                 kite : generic undirected graphs with 5 vertices
                 generic : generic digraph with 6 vertices
                 genus2 : Undirected graph of genus 2
                 ci1 : complete intersection, non-DAG but equivalent to a DAG
                 riemann-roch1 : directed graph with postulation 9 and 3 maximal weight superstables
                 riemann-roch2 : directed graph with a superstable not majorized by a maximal superstable
                 gor : Gorenstein but not a complete intersection


            sage: S = sandlib('gor')
            sage: S.resolution()
            'R^1 <-- R^5 <-- R^5 <-- R^1'
    """
    # The convention is for the sink to be zero.
    sandpiles = {
        'generic':{
                   'description':'generic digraph with 6 vertices',
                   'graph':{0:{},1:{0:1,3:1,4:1},2:{0:1,3:1,5:1},3:{2:1,5:1},4:{1:1,3:1},5:{2:1,3:1}}
                  },
        'kite':{
                'description':'generic undirected graphs with 5 vertices',
                'graph':{0:{}, 1:{0:1,2:1,3:1}, 2:{1:1,3:1,4:1}, 3:{1:1,2:1,4:1},
                         4:{2:1,3:1}}
               },
        'riemann-roch1':{
                         'description':'directed graph with postulation 9 and 3 maximal weight superstables',
                         'graph':{0: {1: 3, 3: 1},
                                  1: {0: 2, 2: 2, 3: 2},
                                  2: {0: 1, 1: 1},
                                  3: {0: 3, 1: 1, 2: 1}
                                 }
                        },
        'riemann-roch2':{
                          'description':'directed graph with a superstable not majorized by a maximal superstable',
                          'graph':{
                                   0: {},
                                   1: {0: 1, 2: 1},
                                   2: {0: 1, 3: 1},
                                   3: {0: 1, 1: 1, 2: 1}
                                  }
                        },
        'gor':{
               'description':'Gorenstein but not a complete intersection',
               'graph':{
                        0: {},
                        1: {0:1, 2: 1, 3: 4},
                        2: {3: 5},
                        3: {1: 1, 2: 1}
                       }
              },
        'ci1':{
               'description':'complete intersection, non-DAG but  equivalent to a DAG',
                   'graph':{0:{}, 1: {2: 2}, 2: {0: 4, 1: 1}}
              },
        'genus2':{
                  'description':'Undirected graph of genus 2',
                  'graph':{
                            0:[1,2],
                            1:[0,2,3],
                            2:[0,1,3],
                            3:[1,2]
                           }
                  },
    }
    if selector is None:
        print
        print '  Sandpiles in the sandlib:'
        for i in sandpiles:
            print '    ', i, ':', sandpiles[i]['description']
        print
    elif selector not in sandpiles.keys():
        print selector, 'is not in the sandlib.'
    else:
        return Sandpile(sandpiles[selector]['graph'], 0)

#################################################
########## Some useful functions ################
#################################################

def complete_sandpile(n):
    r"""
    The sandpile on the complete graph with n vertices.

    INPUT:

    ``n`` - positive integer

    OUTPUT:

    Sandpile

    EXAMPLES::

        sage: K = complete_sandpile(5)
        sage: K.betti(verbose=False)
        [1, 15, 50, 60, 24]
    """
    return Sandpile(graphs.CompleteGraph(n), 0)

def grid_sandpile(m,n):
    """
    The mxn grid sandpile.  Each nonsink vertex has degree 4.

    INPUT:

    ``m``, ``n`` - positive integers

    OUTPUT:

    Sandpile with sink named ``sink``.

    EXAMPLES::

        sage: G = grid_sandpile(3,4)
        sage: G.dict()
        {'sink': {},
         (1, 1): {'sink': 2, (1, 2): 1, (2, 1): 1},
         (1, 2): {'sink': 1, (1, 1): 1, (1, 3): 1, (2, 2): 1},
         (1, 3): {'sink': 1, (1, 2): 1, (1, 4): 1, (2, 3): 1},
         (1, 4): {'sink': 2, (1, 3): 1, (2, 4): 1},
         (2, 1): {'sink': 1, (1, 1): 1, (2, 2): 1, (3, 1): 1},
         (2, 2): {(1, 2): 1, (2, 1): 1, (2, 3): 1, (3, 2): 1},
         (2, 3): {(1, 3): 1, (2, 2): 1, (2, 4): 1, (3, 3): 1},
         (2, 4): {'sink': 1, (1, 4): 1, (2, 3): 1, (3, 4): 1},
         (3, 1): {'sink': 2, (2, 1): 1, (3, 2): 1},
         (3, 2): {'sink': 1, (2, 2): 1, (3, 1): 1, (3, 3): 1},
         (3, 3): {'sink': 1, (2, 3): 1, (3, 2): 1, (3, 4): 1},
         (3, 4): {'sink': 2, (2, 4): 1, (3, 3): 1}}
        sage: G.group_order()
        4140081
        sage: G.invariant_factors()
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1380027]
    """
    g = {}
    # corners first
    g[(1,1)] = {(1,2):1, (2,1):1, 'sink':2}
    g[(m,1)] = {(m-1,1):1, (m,2):1, 'sink':2}
    g[(1,n)] = {(1,n-1):1, (2,n):1, 'sink':2}
    g[(m,n)] = {(m-1,n):1, (m,n-1):1, 'sink':2}
    # top edge
    for col in range(2,n):
        g[(1,col)] = {(1,col-1):1, (1,col+1):1, (2,col):1, 'sink':1}
    # left edge
    for row in range (2,m):
        g[(row,1)] = {(row-1,1):1, (row+1,1):1, (row,2):1, 'sink':1}
    # right edge
    for row in range (2,m):
        g[(row,n)] = {(row-1,n):1, (row+1,n):1, (row,n-1):1, 'sink':1}
    # bottom edge
    for col in range(2,n):
        g[(m,col)] = {(m,col-1):1, (m,col+1):1, (m-1,col):1, 'sink':1}
    # inner vertices
    for row in range(2,m):
        for col in range(2,n):
            g[(row,col)] ={(row-1,col):1, (row+1,col):1, (row,col-1):1, (row,col+1):1}
    # the sink vertex
    g['sink'] = {}
    return Sandpile(g, 'sink')

def triangle_sandpile(n):
    r"""
    A triangular sandpile.  Each nonsink vertex has out-degree six.  The
    vertices on the boundary of the triangle are connected to the sink.

    INPUT:

    ``n`` - int

    OUTPUT:

    Sandpile

    EXAMPLES::

        sage: T = triangle_sandpile(5)
        sage: T.group_order()
        135418115000
    """
    T = {'sink':{}}
    for i in range(n):
        for j in range(n-i):
            T[(i,j)] = {}
            if i<n-j-1:
                T[(i,j)][(i+1,j)] = 1
                T[(i,j)][(i,j+1)] = 1
            if i>0:
                T[(i,j)][(i-1,j+1)] = 1
                T[(i,j)][(i-1,j)] = 1
            if j>0:
                T[(i,j)][(i,j-1)] = 1
                T[(i,j)][(i+1,j-1)] = 1
            d = len(T[(i,j)])
            if d<6:
                T[(i,j)]['sink'] = 6-d
    T = Sandpile(T,'sink')
    pos = {}
    for x in T.nonsink_vertices():
        coords = list(x)
        coords[0]+=QQ(1)/2*coords[1]
        pos[x] = coords
    pos['sink'] = (-1,-1)
    T.set_pos(pos)
    return T

def aztec_sandpile(n):
    r"""
    The aztec diamond graph.

    INPUT:

    ``n`` - integer

    OUTPUT:

    dictionary for the aztec diamond graph

    EXAMPLES::

        sage: aztec_sandpile(2)
        {'sink': {(-3/2, -1/2): 2,
          (-3/2, 1/2): 2,
          (-1/2, -3/2): 2,
          (-1/2, 3/2): 2,
          (1/2, -3/2): 2,
          (1/2, 3/2): 2,
          (3/2, -1/2): 2,
          (3/2, 1/2): 2},
         (-3/2, -1/2): {'sink': 2, (-3/2, 1/2): 1, (-1/2, -1/2): 1},
         (-3/2, 1/2): {'sink': 2, (-3/2, -1/2): 1, (-1/2, 1/2): 1},
         (-1/2, -3/2): {'sink': 2, (-1/2, -1/2): 1, (1/2, -3/2): 1},
         (-1/2, -1/2): {(-3/2, -1/2): 1,
          (-1/2, -3/2): 1,
          (-1/2, 1/2): 1,
          (1/2, -1/2): 1},
         (-1/2, 1/2): {(-3/2, 1/2): 1, (-1/2, -1/2): 1, (-1/2, 3/2): 1, (1/2, 1/2): 1},
         (-1/2, 3/2): {'sink': 2, (-1/2, 1/2): 1, (1/2, 3/2): 1},
         (1/2, -3/2): {'sink': 2, (-1/2, -3/2): 1, (1/2, -1/2): 1},
         (1/2, -1/2): {(-1/2, -1/2): 1, (1/2, -3/2): 1, (1/2, 1/2): 1, (3/2, -1/2): 1},
         (1/2, 1/2): {(-1/2, 1/2): 1, (1/2, -1/2): 1, (1/2, 3/2): 1, (3/2, 1/2): 1},
         (1/2, 3/2): {'sink': 2, (-1/2, 3/2): 1, (1/2, 1/2): 1},
         (3/2, -1/2): {'sink': 2, (1/2, -1/2): 1, (3/2, 1/2): 1},
         (3/2, 1/2): {'sink': 2, (1/2, 1/2): 1, (3/2, -1/2): 1}}
        sage: Sandpile(aztec_sandpile(2),'sink').group_order()
        4542720

    NOTES:

    This is the aztec diamond graph with a sink vertex added.  Boundary
    vertices have edges to the sink so that each vertex has degree 4.
    """
    aztec_sandpile = {}
    half = QQ(1)/2
    for i in srange(n):
        for j in srange(n-i):
            aztec_sandpile[(half+i,half+j)] = {}
            aztec_sandpile[(-half-i,half+j)] = {}
            aztec_sandpile[(half+i,-half-j)] = {}
            aztec_sandpile[(-half-i,-half-j)] = {}
    non_sinks = aztec_sandpile.keys()
    aztec_sandpile['sink'] = {}
    for vert in non_sinks:
        weight = abs(vert[0]) + abs(vert[1])
        x = vert[0]
        y = vert[1]
        if weight < n:
            aztec_sandpile[vert] = {(x+1,y):1, (x,y+1):1, (x-1,y):1, (x,y-1):1}
        else:
            if (x+1,y) in aztec_sandpile.keys():
                aztec_sandpile[vert][(x+1,y)] = 1
            if (x,y+1) in aztec_sandpile.keys():
                aztec_sandpile[vert][(x,y+1)] = 1
            if (x-1,y) in aztec_sandpile.keys():
                aztec_sandpile[vert][(x-1,y)] = 1
            if (x,y-1) in aztec_sandpile.keys():
                aztec_sandpile[vert][(x,y-1)] = 1
            if len(aztec_sandpile[vert]) < 4:
                out_degree = 4 - len(aztec_sandpile[vert])
                aztec_sandpile[vert]['sink'] = out_degree
                aztec_sandpile['sink'][vert] = out_degree
    return aztec_sandpile

def random_digraph(num_verts, p=0.5, directed=True, weight_max=1):
    """
    A random weighted digraph with a directed spanning tree rooted at `0`.  If
    ``directed = False``, the only difference is that if `(i,j,w)` is an edge with
    tail `i`, head `j`, and weight `w`, then `(j,i,w)` appears also.  The result
    is returned as a Sage digraph.

    INPUT:

     - ``num_verts`` - number of vertices
     - ``p`` - probability edges occur
     - ``directed`` - True if directed
     - ``weight_max`` - integer maximum for random weights

    OUTPUT:

    random graph

    EXAMPLES::

        sage: g = random_digraph(6,0.2,True,3)
        sage: S = Sandpile(g,0)
        sage: S.show(edge_labels = True)

    TESTS:

    Check that we can construct a random digraph with the
    default arguments (:trac:`12181`)::

        sage: random_digraph(5)
        Digraph on 5 vertices
    """

    a = digraphs.RandomDirectedGN(num_verts)
    b = graphs.RandomGNP(num_verts,p)
    a.add_edges(b.edges())
    if directed:
        c = graphs.RandomGNP(num_verts,p)
        # reverse the edges of c and add them in
        a.add_edges([(j,i,None) for i,j,k in c.edges()])
    else:
        a.add_edges([(j,i,None) for i,j,k in a.edges()])
        a.add_edges([(j,i,None) for i,j,k in b.edges()])
    # now handle the weights
    for i,j,k in a.edge_iterator():
        a.set_edge_label(i,j,ZZ.random_element(weight_max)+1)
    return a

def random_DAG(num_verts, p=0.5, weight_max=1):
    r"""
    A random directed acyclic graph with ``num_verts`` vertices.
    The method starts with the sink vertex and adds vertices one at a time.
    Each vertex is connected only to only previously defined vertices, and the
    probability of each possible connection is given by the argument ``p``.
    The weight of an edge is a random integer between ``1`` and
    ``weight_max``.

    INPUT:

     - ``num_verts`` - positive integer
     - ``p`` - real number such that `0 < p \leq 1`
     - ``weight_max`` - positive integer

    OUTPUT:

    a dictionary, encoding the edges of a directed acyclic graph with sink `0`

    EXAMPLES::

        sage: d = DiGraph(random_DAG(5, .5));d
        Digraph on 5 vertices

    TESTS:

    Check that we can construct a random DAG with the
    default arguments (:trac:`12181`)::

        sage: g = random_DAG(5);DiGraph(g)
        Digraph on 5 vertices

    Check that bad inputs are rejected::
        sage: g = random_DAG(5,1.1)
        Traceback (most recent call last):
        ...
        ValueError: The parameter p must satisfy 0 < p <= 1.
        sage: g = random_DAG(5,0.1,-1)
        Traceback (most recent call last):
        ...
        ValueError: The parameter weight_max must be positive.
    """
    if not(0 < p and p <= 1):
        raise ValueError("The parameter p must satisfy 0 < p <= 1.")
    weight_max=ZZ(weight_max)
    if not(0 < weight_max):
        raise ValueError("The parameter weight_max must be positive.")
    g = {0:{}}
    for i in range(1,num_verts):
        out_edges = {}
        while out_edges == {}:
            for j in range(i):
                if p > random():
                    out_edges[j] = randint(1,weight_max)
        g[i] = out_edges
    return g

def random_tree(n,d):
    r"""
    A random undirected tree with ``n`` nodes, no node having
    degree higher than ``d``.

    INPUT:

    ``n``, ``d`` - integers

    OUTPUT:

    Graph

    EXAMPLES::

        sage: T = random_tree(15,3)
        sage: T.show()
        sage: S = Sandpile(T,0)
        sage: U = S.reorder_vertices()
        sage: U.show()
    """
    g = Graph()
    # active vertices
    active = [0]
    g.add_vertex(0)
    next_vertex = 1
    while g.num_verts()<n:
        node = randint(0,g.num_verts()-1)
        if g.degree(node)>d:
            active.remove(node)
            break
        r = randint(0,d)
        if r>0:
            for i in range(r):
                g.add_vertex(next_vertex)
                g.add_edge((node,next_vertex))
                active.append(next_vertex)
                next_vertex+=1
    return g

def glue_graphs(g,h,glue_g,glue_h):
    r"""
    Glue two graphs together.

    INPUT:

     - ``g``, ``h`` - dictionaries for directed multigraphs
     - ``glue_h``, ``glue_g`` - dictionaries for a vertex

    OUTPUT:

    dictionary for a directed multigraph


    EXAMPLES::

        sage: x = {0: {}, 1: {0: 1}, 2: {0: 1, 1: 1}, 3: {0: 1, 1: 1, 2: 1}}
        sage: y = {0: {}, 1: {0: 2}, 2: {1: 2}, 3: {0: 1, 2: 1}}
        sage: glue_x = {1: 1, 3: 2}
        sage: glue_y = {0: 1, 1: 2, 3: 1}
        sage: z = glue_graphs(x,y,glue_x,glue_y)
        sage: z
        {0: {},
         'x0': {0: 1, 'x1': 1, 'x3': 2, 'y1': 2, 'y3': 1},
         'x1': {'x0': 1},
         'x2': {'x0': 1, 'x1': 1},
         'x3': {'x0': 1, 'x1': 1, 'x2': 1},
         'y1': {0: 2},
         'y2': {'y1': 2},
         'y3': {0: 1, 'y2': 1}}
        sage: S = Sandpile(z,0)
        sage: S.h_vector()
        [1, 6, 17, 31, 41, 41, 31, 17, 6, 1]
        sage: S.resolution()
        'R^1 <-- R^7 <-- R^21 <-- R^35 <-- R^35 <-- R^21 <-- R^7 <-- R^1'

    NOTES:

    This method makes a dictionary for a graph by combining those for
    ``g`` and ``h``.  The sink of ``g`` is replaced by a vertex that
    is connected to the vertices of ``g`` as specified by ``glue_g``
    the vertices of ``h`` as specified in ``glue_h``.  The sink of the glued
    graph is `0`.

    Both ``glue_g`` and ``glue_h`` are dictionaries with entries of the form
    ``v:w`` where ``v`` is the vertex to be connected to and ``w`` is the weight
    of the connecting edge.
    """
    # first find the sinks of g and h
    for i in g:
        if g[i] == {}:
            g_sink = i
            break
    for i in h:
        if h[i] == {}:
            h_sink = i
            break
    k = {0: {}}  # the new graph dictionary, starting with the sink
    for i in g:
        if i != g_sink:
            new_edges = {}
            for j in g[i]:
                new_edges['x'+str(j)] = g[i][j]
            k['x'+str(i)] = new_edges
    for i in h:
        if i != h_sink:
            new_edges = {}
            for j in h[i]:
                if j == h_sink:
                    new_edges[0] = h[i][j]
                else:
                    new_edges['y'+str(j)] = h[i][j]
            k['y'+str(i)] = new_edges
    # now handle the glue vertex (old g sink)
    new_edges = {}
    for i in glue_g:
        new_edges['x'+str(i)] = glue_g[i]
    for i in glue_h:
        if i == h_sink:
            new_edges[0] = glue_h[i]
        else:
            new_edges['y'+str(i)] = glue_h[i]
    k['x'+str(g_sink)] = new_edges
    return k

def firing_graph(S,eff):
    r"""
    Creates a digraph with divisors as vertices and edges between two
    divisors ``D`` and ``E`` if firing a single vertex in ``D`` gives
    ``E``.

    INPUT:

    ``S`` - sandpile
    ``eff`` - list of divisors

    OUTPUT:

    DiGraph

    EXAMPLES::

        sage: S = Sandpile(graphs.CycleGraph(6),0)
        sage: D = SandpileDivisor(S, [1,1,1,1,2,0])
        sage: eff = D.effective_div() # optional - 4ti2
        sage: firing_graph(S,eff).show3d(edge_size=.005,vertex_size=0.01) # optional - 4ti2
    """
    g = DiGraph()
    g.add_vertices(range(len(eff)))
    for i in g.vertices():
        for v in eff[i]:
            if eff[i][v]>=S.out_degree(v):
                new_div = deepcopy(eff[i])
                new_div[v] -= S.out_degree(v)
                for oe in S.outgoing_edges(v):
                    new_div[oe[1]]+=oe[2]
                if new_div in eff:
                    g.add_edge((i,eff.index(new_div)))
    return g

def parallel_firing_graph(S, eff):
    r"""
    Creates a digraph with divisors as vertices and edges between two
    divisors ``D`` and ``E`` if firing all unstable vertices in ``D`` gives
    ``E``.

    INPUT:

    ``S`` - Sandpile
    ``eff`` - list of divisors

    OUTPUT:

    DiGraph

    EXAMPLES::

        sage: S = Sandpile(graphs.CycleGraph(6),0)
        sage: D = SandpileDivisor(S, [1,1,1,1,2,0])
        sage: eff = D.effective_div() # optional - 4ti2
        sage: parallel_firing_graph(S,eff).show3d(edge_size=.005,vertex_size=0.01) # optional - 4ti2
    """
    g = DiGraph()
    g.add_vertices(range(len(eff)))
    for i in g.vertices():
        new_edge = False
        new_div = deepcopy(eff[i])
        for v in eff[i]:
            if eff[i][v]>=S.out_degree(v):
                new_edge = True
                new_div[v] -= S.out_degree(v)
                for oe in S.outgoing_edges(v):
                    new_div[oe[1]]+=oe[2]
        if new_edge and (new_div in eff):
            g.add_edge((i,eff.index(new_div)))
    return g

def admissible_partitions(S, k):
    r"""
    The partitions of the vertices of ``S`` into ``k`` parts,
    each of which is connected.

    INPUT:

    ``S`` - Sandpile
    ``k`` - integer

    OUTPUT:

    list of partitions

    EXAMPLES::

        sage: S = Sandpile(graphs.CycleGraph(4), 0)
        sage: P = [admissible_partitions(S, i) for i in [2,3,4]]
        sage: P
        [[{{0}, {1, 2, 3}},
          {{0, 2, 3}, {1}},
          {{0, 1, 3}, {2}},
          {{0, 1, 2}, {3}},
          {{0, 1}, {2, 3}},
          {{0, 3}, {1, 2}}],
         [{{0}, {1}, {2, 3}},
          {{0}, {1, 2}, {3}},
          {{0, 3}, {1}, {2}},
          {{0, 1}, {2}, {3}}],
         [{{0}, {1}, {2}, {3}}]]
        sage: for p in P:
        ...    sum([partition_sandpile(S, i).betti(verbose=False)[-1] for i in p])
        6
        8
        3
        sage: S.betti()
                   0     1     2     3
        ------------------------------
            0:     1     -     -     -
            1:     -     6     8     3
        ------------------------------
        total:     1     6     8     3
    """
    v = S.vertices()
    if S.is_directed():
        G = DiGraph(S)
    else:
        G = Graph(S)
    result = []
    for p in SetPartitions(v, k):
        if forall(p, lambda x : G.subgraph(list(x)).is_connected())[0]:
            result.append(p)
    return result

def partition_sandpile(S,p):
    r"""
    Each set of vertices in ``p`` is regarded as a single vertex, with and edge
    between ``A`` and ``B`` if some element of ``A`` is connected by an edge
    to  some element of ``B`` in ``S``.

    INPUT:

    ``S`` - Sandpile
    ``p`` - partition of the vertices of ``S``

    OUTPUT:

    Sandpile

    EXAMPLES::

        sage: S = Sandpile(graphs.CycleGraph(4), 0)
        sage: P = [admissible_partitions(S, i) for i in [2,3,4]]
        sage: for p in P:
        ...    sum([partition_sandpile(S, i).betti(verbose=False)[-1] for i in p])
        6
        8
        3
        sage: S.betti()
                   0     1     2     3
        ------------------------------
            0:     1     -     -     -
            1:     -     6     8     3
        ------------------------------
        total:     1     6     8     3
    """
    from sage.combinat.combination import Combinations
    g = Graph()
    g.add_vertices([tuple(i) for i in p])
    for u,v in Combinations(g.vertices(), 2):
        for i in u:
            for j in v:
                if (i,j,1) in S.edges():
                    g.add_edge((u, v))
                    break
    for i in g.vertices():
        if S.sink() in i:
            return Sandpile(g,i)

def firing_vector(S,D,E):
    r"""
    If ``D`` and ``E`` are linearly equivalent divisors, find the firing vector
    taking ``D`` to ``E``.

    INPUT:

    - ``S`` - Sandpile
    - ``D``, ``E`` - tuples (representing linearly equivalent divisors)

    OUTPUT:

    tuple (representing a firing vector from ``D`` to ``E``)

    EXAMPLES::

      sage: S = complete_sandpile(4)
      sage: D = SandpileDivisor(S, {0: 0, 1: 0, 2: 8, 3: 0})
      sage: E = SandpileDivisor(S, {0: 2, 1: 2, 2: 2, 3: 2})
      sage: v = firing_vector(S, D, E)
      sage: v
      (0, 0, 2, 0)

    The divisors must be linearly equivalent::

      sage: vector(D.values()) - S.laplacian()*vector(v) == vector(E.values())
      True
      sage: firing_vector(S, D, S.zero_div())
      Error. Are the divisors linearly equivalent?
  """
    try:
        v = vector(D.values())
        w = vector(E.values())
        return tuple(S.laplacian().solve_left(v-w))
    except ValueError:
        print "Error. Are the divisors linearly equivalent?"
        return

def min_cycles(G,v):
    r"""
    Minimal length cycles in the digraph ``G`` starting at vertex ``v``.

    INPUT:

    ``G`` - DiGraph
    ``v`` - vertex of ``G``

    OUTPUT:

    list of lists of vertices

    EXAMPLES::

        sage: T = sandlib('gor')
        sage: [min_cycles(T, i) for i in T.vertices()]
        [[], [[1, 3]], [[2, 3, 1], [2, 3]], [[3, 1], [3, 2]]]
    """
    pr = G.neighbors_in(v)
    sp = G.shortest_paths(v)
    return [sp[i] for i in pr if i in sp.keys()]

def wilmes_algorithm(M):
    r"""
    Computes an integer matrix ``L`` with the same integer row span as ``M``
    and such that ``L`` is the reduced laplacian of a directed multigraph.

    INPUT:

    ``M`` - square integer matrix of full rank

    OUTPUT:

    ``L`` - integer matrix

    EXAMPLES::

        sage: P = matrix([[2,3,-7,-3],[5,2,-5,5],[8,2,5,4],[-5,-9,6,6]])
        sage: wilmes_algorithm(P)
        [ 1642   -13 -1627    -1]
        [   -1  1980 -1582  -397]
        [    0    -1  1650 -1649]
        [    0     0 -1658  1658]

    NOTES:

    The algorithm is due to John Wilmes.
    """
    # find the gcd of the row-sums, and perform the corresponding row
    # operations on M
    if M.matrix_over_field().is_invertible():
        L = deepcopy(M)
        L = matrix(ZZ,L)
        U = matrix(ZZ,[sum(i) for i in L]).smith_form()[2].transpose()
        L = U*M
        for k in range(1,M.nrows()-1):
            smith = matrix(ZZ,[i[k-1] for i in L[k:]]).smith_form()[2].transpose()
            U = identity_matrix(ZZ,k).block_sum(smith)
            L = U*L
            L[k] = -L[k]
        if L[-1][-2]>0:
            L[-1] = -L[-1]
        for k in range(M.nrows()-2,-1,-1):
            for i in range(k+2,M.nrows()):
                while L[k][i-1]>0:
                    L[k] = L[k] + L[i]
            v = -L[k+1]
            for i in range(k+2,M.nrows()):
                v = abs(L[i,i-1])*v + v[i-1]*L[i]
            while L[k,k]<=0 or L[k,-1]>0:
                L[k] = L[k] + v
        return L
    else:
        raise UserWarning('matrix not of full rank')

######### Notes ################
"""
* pickling sandpiles?
* Does 'sat' return a minimal generating set
"""

