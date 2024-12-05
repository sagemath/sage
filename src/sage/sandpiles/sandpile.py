"""
Sandpiles

Functions and classes for mathematical sandpiles.

Version: 2.4

AUTHOR:

- David Perkinson (June 4, 2015) Upgraded from version 2.3 to 2.4.

MAJOR CHANGES

1. Eliminated dependence on 4ti2, substituting the use of Polyhedron methods. Thus, no optional packages are necessary.
2. Fixed bug in ``Sandpile.__init__`` so that now multigraphs are handled correctly.
3. Created ``sandpiles`` to handle examples of Sandpiles in analogy with ``graphs``, ``simplicial_complexes``, and ``polytopes``.  In the process, we implemented a much faster way of producing the sandpile grid graph.
4. Added support for open and closed sandpile Markov chains.
5. Added support for Weierstrass points.
6. Implemented the Cori-Le Borgne algorithm for computing ranks of divisors on complete graphs.

NEW METHODS

**Sandpile**: avalanche_polynomial, genus, group_gens, help, jacobian_representatives, markov_chain, picard_representatives, smith_form, stable_configs, stationary_density, tutte_polynomial.

**SandpileConfig**: burst_size, help.

**SandpileDivisor**:  help, is_linearly_equivalent, is_q_reduced, is_weierstrass_pt, polytope, polytope_integer_pts, q_reduced, rank, simulate_threshold, stabilize, weierstrass_div, weierstrass_gap_seq, weierstrass_pts, weierstrass_rank_seq.

MINOR CHANGES

* The ``sink`` argument to ``Sandpile.__init__`` now defaults to the first vertex.
* A SandpileConfig or SandpileDivisor may now be multiplied by an integer.
* Sped up ``__add__`` method for SandpileConfig and SandpileDivisor.
* Enhanced string representation of a Sandpile (via ``_repr_`` and the ``name`` methods).
* Recurrents for complete graphs and cycle graphs are computed more quickly.
* The stabilization code for SandpileConfig has been made more efficient.
* Added optional probability distribution arguments to ``add_random`` methods.

---------------------------------------

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

For general help, enter ``Sandpile.help()``, ``SandpileConfig.help()``, and
``SandpileDivisor.help()``.  Miscellaneous examples appear below.

A weighted directed graph given as a Python dictionary::

    sage: from sage.sandpiles import *
    sage: g = {0: {},
    ....:      1: {0: 1, 2: 1, 3: 1},
    ....:      2: {1: 1, 3: 1, 4: 1},
    ....:      3: {1: 1, 2: 1, 4: 1},
    ....:      4: {2: 1, 3: 1}}

The associated sandpile with 0 chosen as the sink::

    sage: S = Sandpile(g,0)

Or just::

    sage: S = Sandpile(g)

A picture of the graph::

    sage: S.show()                              # long time                             # needs sage.plot

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

An arbitrary sandpile configuration::

    sage: c = SandpileConfig(S,[1,0,4,-3])
    sage: c.equivalent_recurrent()
    {1: 2, 2: 2, 3: 2, 4: 0}

Some group operations::

    sage: m = S.max_stable()
    sage: i = S.identity()
    sage: m.values()
    [2, 2, 2, 1]
    sage: i.values()
    [2, 2, 2, 0]
    sage: m + i    # coordinate-wise sum
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

    sage: S.h_vector()
    [1, 3, 4]
    sage: S.postulation()
    2

the saturated homogeneous toppling ideal::

    sage: S.ideal()                                                                     # needs sage.libs.singular
    Ideal (x1 - x0, x3*x2 - x0^2, x4^2 - x0^2, x2^3 - x4*x3*x0,
           x4*x2^2 - x3^2*x0, x3^3 - x4*x2*x0, x4*x3^2 - x2^2*x0) of
     Multivariate Polynomial Ring in x4, x3, x2, x1, x0 over Rational Field

its minimal free resolution::

    sage: S.resolution()                                                                # needs sage.libs.singular
    'R^1 <-- R^7 <-- R^15 <-- R^13 <-- R^4'

and its Betti numbers::

    sage: S.betti()                                                                     # needs sage.libs.singular
               0     1     2     3     4
    ------------------------------------
        0:     1     1     -     -     -
        1:     -     2     2     -     -
        2:     -     4    13    13     4
    ------------------------------------
    total:     1     7    15    13     4

Some various ways of creating Sandpiles::

    sage: S = sandpiles.Complete(4) # for more options enter ``sandpile.TAB``
    sage: S = sandpiles.Wheel(6)

A multidigraph with loops (vertices 0, 1, 2; for example, there is a directed
edge from vertex 2 to vertex 1 of weight 3, which can be thought of as three
directed edges of the form (2,3).  There is also a single loop at vertex 2 and
an edge (2,0) of weight 2)::

    sage: S = Sandpile({0:[1,2], 1:[0,0,2], 2:[0,0,1,1,1,2], 3:[2]})

Using the graph library (vertex 1 is specified as the sink; omitting
this would make the sink vertex 0 by default)::

    sage: S = Sandpile(graphs.PetersenGraph(),1)

Distribution of avalanche sizes::

    sage: S = sandpiles.Grid(10,10)
    sage: m = S.max_stable()
    sage: a = []
    sage: for i in range(1000):
    ....:     m = m.add_random()
    ....:     m, f = m.stabilize(True)
    ....:     a.append(sum(f.values()))

    sage: # needs sage.plot
    sage: p = list_plot([[log(i + 1), log(a.count(i))]
    ....:                for i in [0..max(a)] if a.count(i)])
    sage: p.axes_labels(['log(N)', 'log(D(N))'])
    sage: t = text("Distribution of avalanche sizes", (2,2), rgbcolor=(1,0,0))
    sage: show(p + t, axes_labels=['log(N)', 'log(D(N))'])      # long time

Working with sandpile divisors::

    sage: S = sandpiles.Complete(4)
    sage: D = SandpileDivisor(S, [0,0,0,5])
    sage: E = D.stabilize(); E
    {0: 1, 1: 1, 2: 1, 3: 2}
    sage: D.is_linearly_equivalent(E)
    True
    sage: D.q_reduced()
    {0: 4, 1: 0, 2: 0, 3: 1}
    sage: S = sandpiles.Complete(4)
    sage: D = SandpileDivisor(S, [0,0,0,5])
    sage: E = D.stabilize(); E
    {0: 1, 1: 1, 2: 1, 3: 2}
    sage: D.is_linearly_equivalent(E)
    True
    sage: D.q_reduced()
    {0: 4, 1: 0, 2: 0, 3: 1}
    sage: D.rank()
    2

    sage: # needs sage.geometry.polyhedron
    sage: sorted(D.effective_div(), key=str)
    [{0: 0, 1: 0, 2: 0, 3: 5},
     {0: 0, 1: 0, 2: 4, 3: 1},
     {0: 0, 1: 4, 2: 0, 3: 1},
     {0: 1, 1: 1, 2: 1, 3: 2},
     {0: 4, 1: 0, 2: 0, 3: 1}]
    sage: sorted(D.effective_div(False))
    [[0, 0, 0, 5], [0, 0, 4, 1], [0, 4, 0, 1], [1, 1, 1, 2], [4, 0, 0, 1]]
    sage: D.rank()
    2
    sage: D.rank(True)
    (2, {0: 2, 1: 1, 2: 0, 3: 0})
    sage: E = D.rank(True)[1]  # E proves the rank is not 3
    sage: E.values()
    [2, 1, 0, 0]
    sage: E.deg()
    3
    sage: rank(D - E)
    -1
    sage: (D - E).effective_div()
    []

    sage: D.weierstrass_pts()
    (0, 1, 2, 3)
    sage: D.weierstrass_rank_seq(0)
    (2, 1, 0, 0, 0, -1)
    sage: D.weierstrass_pts()
    (0, 1, 2, 3)
    sage: D.weierstrass_rank_seq(0)
    (2, 1, 0, 0, 0, -1)
"""

# ****************************************************************************
#       Copyright (C) 2011 David Perkinson <davidp@reed.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from collections import Counter
from copy import deepcopy
from inspect import getdoc
from textwrap import dedent

from IPython.lib import pretty

from sage.calculus.functional import derivative
from sage.combinat.integer_vector import integer_vectors_nk_fast_iter
from sage.combinat.parking_functions import ParkingFunctions
from sage.combinat.set_partition import SetPartitions
from sage.combinat.vector_partition import IntegerVectorsIterator
from sage.functions.log import exp
from sage.arith.misc import binomial
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.probability.probability_distribution import GeneralDiscreteDistribution
from sage.topology.simplicial_complex import SimplicialComplex
from sage.interfaces.singular import singular
from sage.matrix.constructor import matrix, identity_matrix
from sage.misc.functional import det, denominator
from sage.misc.misc import exists
from sage.misc.misc_c import prod
from sage.misc.temporary_file import tmp_filename
from sage.arith.srange import xsrange
from sage.modules.free_module_element import vector
from sage.misc.lazy_import import lazy_import
from sage.arith.functions import lcm
from sage.arith.misc import falling_factorial
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.constants import I, pi
from sage.symbolic.ring import SR
from sage.features.four_ti_2 import FourTi2Executable
lazy_import("sage.plot.colors", "rainbow")


def _sandpile_help(cls, usage, verbose=True):
    """
    Print help text for classes in this module; see the ``help()`` methods on
    individual classes in this module for example usage.
    """

    # We collect the first sentence of each docstring.  The sentence is,
    # by definition, from the beginning of the string to the first
    # occurrence of a period or question mark.  If neither of these appear
    # in the string, take the sentence to be the empty string.  If the
    # latter occurs, something should be changed.
    from sage.misc.sagedoc import detex
    methods = []
    for attr in sorted(vars(cls)):
        if attr[0] != '_':
            doc = getdoc(getattr(cls, attr))
            period = doc.find('.')
            question = doc.find('?')
            if period == -1 and question == -1:
                doc = ''  # Neither appears!
            else:
                if period == -1:
                    period = len(doc) + 1
                if question == -1:
                    question = len(doc) + 1
                if period < question:
                    doc = doc.split('.')[0]
                    doc = detex(doc).strip() + '.'
                else:
                    doc = doc.split('?')[0]
                    doc = detex(doc).strip() + '?'
            methods.append((attr, doc))

    print(usage)
    print()

    mlen = max(len(attr) for attr, doc in methods)
    if verbose:
        for attr, doc in methods:
            print(attr.ljust(mlen), '--', doc)
    else:
        for attr, _ in methods:
            print(attr)


class Sandpile(DiGraph):
    """
    Class for Dhar's abelian sandpile model.
    """
    @staticmethod
    def version():
        r"""
        The version number of Sage Sandpiles.

        OUTPUT: string

        EXAMPLES::

            sage: Sandpile.version()
            Sage Sandpiles Version 2.4
            sage: S = sandpiles.Complete(3)
            sage: S.version()
            Sage Sandpiles Version 2.4
        """
        print('Sage Sandpiles Version 2.4')

    @staticmethod
    def help(verbose=True):
        r"""
        List of Sandpile-specific methods (not inherited from :class:`Graph`).
        If ``verbose``, include short descriptions.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: printed string

        EXAMPLES::

            sage: Sandpile.help() # long time
            For detailed help with any method FOO listed below,
            enter "Sandpile.FOO?" or enter "S.FOO?" for any Sandpile S.
            <BLANKLINE>
            all_k_config             -- The constant configuration with all values set to k.
            all_k_div                -- The divisor with all values set to k.
            avalanche_polynomial     -- The avalanche polynomial.
            betti                    -- The Betti table for the homogeneous toppling ideal.
            betti_complexes          -- The support-complexes with non-trivial homology.
            burning_config           -- The minimal burning configuration.
            burning_script           -- A script for the minimal burning configuration.
            canonical_divisor        -- The canonical divisor.
            dict                     -- A dictionary of dictionaries representing a directed graph.
            genus                    -- The genus: (# non-loop edges) - (# vertices) + 1.
            groebner                 -- A Groebner basis for the homogeneous toppling ideal.
            group_gens               -- A minimal list of generators for the sandpile group.
            group_order              -- The size of the sandpile group.
            h_vector                 -- The number of superstable configurations in each degree.
            help                     -- List of Sandpile-specific methods (not inherited from ...Graph...).
            hilbert_function         -- The Hilbert function of the homogeneous toppling ideal.
            ideal                    -- The saturated homogeneous toppling ideal.
            identity                 -- The identity configuration.
            in_degree                -- The in-degree of a vertex or a list of all in-degrees.
            invariant_factors        -- The invariant factors of the sandpile group.
            is_undirected            -- Is the underlying graph undirected?
            jacobian_representatives -- Representatives for the elements of the Jacobian group.
            laplacian                -- The Laplacian matrix of the graph.
            markov_chain             -- The sandpile Markov chain for configurations or divisors.
            max_stable               -- The maximal stable configuration.
            max_stable_div           -- The maximal stable divisor.
            max_superstables         -- The maximal superstable configurations.
            min_recurrents           -- The minimal recurrent elements.
            nonsink_vertices         -- The nonsink vertices.
            nonspecial_divisors      -- The nonspecial divisors.
            out_degree               -- The out-degree of a vertex or a list of all out-degrees.
            picard_representatives   -- Representatives of the divisor classes of degree d in the Picard group.
            points                   -- Generators for the multiplicative group of zeros of the sandpile ideal.
            postulation              -- The postulation number of the toppling ideal.
            recurrents               -- The recurrent configurations.
            reduced_laplacian        -- The reduced Laplacian matrix of the graph.
            reorder_vertices         -- A copy of the sandpile with vertex names permuted.
            resolution               -- A minimal free resolution of the homogeneous toppling ideal.
            ring                     -- The ring containing the homogeneous toppling ideal.
            show                     -- Draw the underlying graph.
            show3d                   -- Draw the underlying graph.
            sink                     -- The sink vertex.
            smith_form               -- The Smith normal form for the Laplacian.
            solve                    -- Approximations of the complex affine zeros of the sandpile ideal.
            stable_configs           -- Generator for all stable configurations.
            stationary_density       -- The stationary density of the sandpile.
            superstables             -- The superstable configurations.
            symmetric_recurrents     -- The symmetric recurrent configurations.
            tutte_polynomial         -- The Tutte polynomial of the underlying graph.
            unsaturated_ideal        -- The unsaturated, homogeneous toppling ideal.
            version                  -- The version number of Sage Sandpiles.
            zero_config              -- The all-zero configuration.
            zero_div                 -- The all-zero divisor.
        """
        _sandpile_help(Sandpile, dedent("""\
            For detailed help with any method FOO listed below,
            enter "Sandpile.FOO?" or enter "S.FOO?" for any Sandpile S."""),
                       verbose=verbose)

    def __init__(self, g, sink=None):
        r"""
        Create a sandpile.

        A sandpile is always a weighted graph.

        INPUT:

        - ``g`` -- dictionary for directed multigraph with edges weighted by
          nonnegative integers (see NOTE), a Graph or DiGraph

        - ``sink`` -- (optional) A sink vertex.  Any outgoing edges from the
          designated sink are ignored for the purposes of stabilization.  It is
          assumed that every vertex has a directed path into the sink.  If the
          ``sink`` argument is omitted, the first vertex in the list of the
          Sandpile's vertices is set as the sink.

        OUTPUT: Sandpile

        EXAMPLES:

        Below, ``g`` represents a square with directed, multiple edges with three
        vertices, ``a``, ``b``, ``c``, and ``d``.  The vertex ``a`` has
        outgoing edges to itself (weight 2), to vertex ``b`` (weight 1), and
        vertex ``c`` (weight 3), for example.

        ::

            sage: g = {'a': {'a':2, 'b':1, 'c':3}, 'b': {'a':1, 'd':1},
            ....:      'c': {'a':1,'d': 1}, 'd': {'b':1, 'c':1}}
            sage: G = Sandpile(g,'d')

        Here is a square with unweighted edges.  In this example, the graph is
        also undirected. ::

            sage: g = {0:[1,2], 1:[0,3], 2:[0,3], 3:[1,2]}
            sage: G = Sandpile(g,3)

        In the following example, multiple edges and loops in the dictionary
        become edge weights in the Sandpile. ::

            sage: s = Sandpile({0:[1,2,3], 1:[0,1,2,2,2], 2:[1,1,0,2,2,2,2]})
            sage: s.laplacian()
            [ 3 -1 -1 -1]
            [-1  4 -3  0]
            [-1 -2  3  0]
            [ 0  0  0  0]
            sage: s.dict()
            {0: {1: 1, 2: 1, 3: 1}, 1: {0: 1, 1: 1, 2: 3}, 2: {0: 1, 1: 2, 2: 4}}

        Sandpiles can be created from Graphs and DiGraphs. ::

            sage: g = DiGraph({0:{1:2,2:4}, 1:{1:3,2:1}, 2:{1:7}}, weighted=True)
            sage: s = Sandpile(g)
            sage: s.dict()
            {0: {1: 2, 2: 4}, 1: {0: 0, 1: 3, 2: 1}, 2: {0: 0, 1: 7}}
            sage: s.sink()
            0
            sage: s = sandpiles.Cycle(4)
            sage: s.laplacian()
            [ 2 -1  0 -1]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [-1  0 -1  2]

        .. NOTE::

            Loops are allowed.  There are four admissible input formats.  Two of
            these are dictionaries whose keys are the vertex names.  In one, the
            values are dictionaries with keys the names of vertices which are the
            heads of outgoing edges and with values the weights of the edges.  In
            the other format, the values are lists of names of vertices which are
            the heads of the outgoing edges, with weights determined by the number
            of times a name of a vertex appears in the list.  Both Graphs and
            DiGraphs can also be used as inputs.

        TESTS::

            sage: S = sandpiles.Complete(4)
            sage: TestSuite(S).run()

        Make sure we cannot make an unweighted sandpile::

            sage: G = Sandpile({0:[]}, 0, weighted=False)
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() got an unexpected keyword argument 'weighted'
        """
        # set graph name
        if isinstance(g, (Graph, DiGraph)):
            name = g.name()
            if name == '':
                name = 'sandpile graph'
            else:
                p = name.lower().find('graph')
                if p == -1:
                    name = name + ' sandpile graph'
                else:
                    name = name[:p] + 'sandpile graph' + name[p + 5:]
            self._name = name
        else:
            self._name = 'sandpile graph'
        # preprocess a graph, if necessary
        if isinstance(g, dict) and isinstance(next(iter(g.values())), dict):
            pass  # this is the default format
        elif isinstance(g, dict) and isinstance(next(iter(g.values())), list):
            processed_g = {i: dict(Counter(g[i])) for i in g}
            g = processed_g
        elif isinstance(g, (Graph, DiGraph)):
            if not g.weighted():
                h = g.to_dictionary(multiple_edges=True)
                g = {i: dict(Counter(h[i])) for i in h}
            else:
                vi = {v: i for i, v in enumerate(g.vertices(sort=True))}
                ad = g.weighted_adjacency_matrix()
                g = {v: {w: ad[vi[v], vi[w]] for w in g.neighbor_iterator(v)}
                     for v in g}
        else:
            raise SyntaxError(g)

        # create digraph and initialize some variables
        DiGraph.__init__(self, g, weighted=True)
        self._dict = deepcopy(g)
        vertices = self.vertices(sort=True)
        if sink is None:
            self._sink = vertices[0]
            self._sink_ind = 0
        else:
            self._sink = sink  # key for sink
            self._sink_ind = vertices.index(sink)
        self._nonsink_vertices = vertices
        del self._nonsink_vertices[self._sink_ind]
        # compute Laplacians:
        self._laplacian = self.laplacian_matrix(indegree=False)
        temp = list(range(self.num_verts()))
        del temp[self._sink_ind]
        self._reduced_laplacian = self._laplacian[temp, temp]

    def __copy__(self):
        """
        Make a copy of this sandpile.

        OUTPUT: a new :class:`Sandpile` instance

        EXAMPLES::

            sage: G = sandpiles.Complete(4)
            sage: G_copy = copy(G)
            sage: G_copy == G == G.__copy__()
            True
        """
        return self.__class__(self, self._sink)

    def __getattr__(self, name):
        """
        Set certain variables only when called.

        INPUT:

        - ``name`` -- name of an internal method

        EXAMPLES::

            sage: S = sandpiles.Complete(5)
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
            elif name == '_group_gens':
                self._set_group_gens()
                return deepcopy(self.__dict__[name])
            elif name == '_group_order':
                self.__dict__[name] = det(self._reduced_laplacian.dense_matrix())
                return self.__dict__[name]
            elif name == '_invariant_factors':
                self._set_invariant_factors()
                return deepcopy(self.__dict__[name])
            elif name == '_smith_form':
                self._set_smith_form()
                return deepcopy(self.__dict__[name])
            elif name == '_jacobian_representatives':
                self._set_jacobian_representatives()
                return deepcopy(self.__dict__[name])
            elif name == '_avalanche_polynomial':
                self._set_avalanche_polynomial()
                return deepcopy(self.__dict__[name])
            elif name == '_stationary_density':
                self._set_stationary_density()
                return self.__dict__[name]
            elif name == '_betti_complexes':
                self._set_betti_complexes()
                return deepcopy(self.__dict__[name])
            elif name in ['_postulation', '_h_vector', '_hilbert_function']:
                self._set_hilbert_function()
                return deepcopy(self.__dict__[name])
            elif (name == '_ring' or name == '_unsaturated_ideal'):
                self._set_ring()
                return self.__dict__[name]
            elif name == '_ideal':
                self._set_ideal()
                return self.__dict__[name]
            elif name in ['_resolution', '_betti', '_singular_resolution']:
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

    def __str__(self) -> str:
        r"""
        The name of the sandpile.

        OUTPUT: string

        EXAMPLES::

            sage: s = Sandpile(graphs.PetersenGraph(),2)
            sage: str(s)
            'Petersen sandpile graph'
            sage: str(sandpiles.House())
            'House sandpile graph'
            sage: str(Sandpile({0:[1,1], 1:[0]}))
            'sandpile graph'
        """
        return self.name()

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: s = Sandpile(graphs.PetersenGraph(),2)
            sage: repr(s)
            'Petersen sandpile graph: 10 vertices, sink = 2'
            sage: repr(sandpiles.Complete(4))
            'Complete sandpile graph: 4 vertices, sink = 0'
            sage: repr(Sandpile({0:[1,1], 1:[0]}))
            'sandpile graph: 2 vertices, sink = 0'
        """
        return self._name + ': ' + str(self.num_verts()) + ' vertices, sink = ' + str(self.sink())

    def show(self, **kwds):
        r"""
        Draw the underlying graph.

        INPUT:

        - ``kwds`` -- (optional) arguments passed to the show method for Graph or DiGraph

        EXAMPLES::

            sage: S = Sandpile({0:[], 1:[0,3,4], 2:[0,3,5], 3:[2,5], 4:[1,1], 5:[2,4]})
            sage: S.show()                                                              # needs sage.plot
            sage: S.show(graph_border=True, edge_labels=True)                           # needs sage.plot
        """
        if self.is_undirected():
            Graph(self).show(**kwds)
        else:
            DiGraph(self).show(**kwds)

    def show3d(self, **kwds):
        r"""
        Draw the underlying graph.

        INPUT:

        - ``kwds`` -- (optional) arguments passed to the show method for Graph or DiGraph

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: S.show3d()                    # long time                             # needs sage.plot
        """

        if self.is_undirected():
            Graph(self).show3d(**kwds)
        else:
            DiGraph(self).show3d(**kwds)

    def dict(self):
        r"""
        A dictionary of dictionaries representing a directed graph.

        OUTPUT: dictionary

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.dict()
            {0: {1: 1, 2: 1},
             1: {0: 1, 2: 1, 3: 1},
             2: {0: 1, 1: 1, 3: 1},
             3: {1: 1, 2: 1}}
            sage: S.sink()
            0
        """
        return deepcopy(self._dict)

    def sink(self):
        r"""
        The sink vertex.

        OUTPUT: sink vertex

        EXAMPLES::

            sage: G = sandpiles.House()
            sage: G.sink()
            0
            sage: H = sandpiles.Grid(2,2)
            sage: H.sink()
            (0, 0)
            sage: type(H.sink())
            <... 'tuple'>
        """
        return self._sink

    def laplacian(self):
        r"""
        The Laplacian matrix of the graph.  Its *rows* encode the vertex firing rules.

        OUTPUT: matrix

        EXAMPLES::

            sage: G = sandpiles.Diamond()
            sage: G.laplacian()
            [ 2 -1 -1  0]
            [-1  3 -1 -1]
            [-1 -1  3 -1]
            [ 0 -1 -1  2]

        .. WARNING::

            The function ``laplacian_matrix`` should be avoided.  It returns the
            indegree version of the Laplacian.
        """
        return deepcopy(self._laplacian)

    def reduced_laplacian(self):
        r"""
        The reduced Laplacian matrix of the graph.

        OUTPUT: matrix

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.laplacian()
            [ 2 -1 -1  0]
            [-1  3 -1 -1]
            [-1 -1  3 -1]
            [ 0 -1 -1  2]
            sage: S.reduced_laplacian()
            [ 3 -1 -1]
            [-1  3 -1]
            [-1 -1  2]

        .. NOTE::

            This is the Laplacian matrix with the row and column indexed by the
            sink vertex removed.
        """
        return deepcopy(self._reduced_laplacian)

    def group_order(self):
        r"""
        The size of the sandpile group.

        OUTPUT: integer

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: S.group_order()
            11
        """
        return self._group_order

    def _set_max_stable(self):
        r"""
        Initialize the variable holding the maximal stable configuration.

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: S._set_max_stable()
            sage: '_max_stable' in S.__dict__
            True
        """
        m = {v: self.out_degree(v) - 1 for v in self._nonsink_vertices}
        self._max_stable = SandpileConfig(self, m)

    def max_stable(self):
        r"""
        The maximal stable configuration.

        OUTPUT: SandpileConfig (the maximal stable configuration)

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: S.max_stable()
            {1: 1, 2: 2, 3: 2, 4: 1}
        """
        return deepcopy(self._max_stable)

    def _set_max_stable_div(self):
        r"""
        Initialize the variable holding the maximal stable divisor.

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S._set_max_stable_div()
            sage: '_max_stable_div' in S.__dict__
            True
        """
        m = {v: self.out_degree(v) - 1 for v in self.vertices(sort=False)}
        self._max_stable_div = SandpileDivisor(self, m)

    def max_stable_div(self):
        r"""
        The maximal stable divisor.

        OUTPUT:

        SandpileDivisor (the maximal stable divisor)

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s.max_stable_div()
            {0: 1, 1: 2, 2: 2, 3: 1}
            sage: s.out_degree()
            {0: 2, 1: 3, 2: 3, 3: 2}
        """
        return deepcopy(self._max_stable_div)

    def _set_out_degrees(self):
        r"""
        Initialize the variable holding the out-degrees.

        EXAMPLES::

            sage: s = sandpiles.House()
            sage: s._set_out_degrees()
            sage: '_out_degrees' in s.__dict__
            True
        """
        self._out_degrees = {v: 0 for v in self.vertices(sort=False)}
        for v in self.vertices(sort=False):
            for e in self.edges_incident(v):
                self._out_degrees[v] += e[2]

    def out_degree(self, v=None):
        r"""
        The out-degree of a vertex or a list of all out-degrees.

        INPUT:

        - ``v`` -- (optional) vertex name

        OUTPUT: integer or dict

        EXAMPLES::

            sage: s = sandpiles.House()
            sage: s.out_degree()
            {0: 2, 1: 2, 2: 3, 3: 3, 4: 2}
            sage: s.out_degree(2)
            3
        """
        if v is not None:
            return self._out_degrees[v]
        return self._out_degrees

    def _set_in_degrees(self):
        """
        Initialize the variable holding the in-degrees.

        EXAMPLES::

            sage: s = sandpiles.House()
            sage: s._set_in_degrees()
            sage: '_in_degrees' in s.__dict__
            True
        """
        self._in_degrees = {v: 0 for v in self}
        for e in self.edge_iterator():
            self._in_degrees[e[1]] += e[2]

    def in_degree(self, v=None):
        r"""
        The in-degree of a vertex or a list of all in-degrees.

        INPUT:

        - ``v`` -- (optional) vertex name

        OUTPUT: integer or dict

        EXAMPLES::

            sage: s = sandpiles.House()
            sage: s.in_degree()
            {0: 2, 1: 2, 2: 3, 3: 3, 4: 2}
            sage: s.in_degree(2)
            3
        """
        if v is not None:
            return self._in_degrees[v]
        return self._in_degrees

    def _set_burning_config(self):
        r"""
        Calculate the minimal burning configuration and its corresponding
        script.

        EXAMPLES::

            sage: g = {0:{},1:{0:1,3:1,4:1},2:{0:1,3:1,5:1},
            ....:      3:{2:1,5:1},4:{1:1,3:1},5:{2:1,3:1}}
            sage: S = Sandpile(g,0)
            sage: S._set_burning_config()
        """
        # TODO: Cythonize!
        d = self._reduced_laplacian.nrows()
        burn = sum(self._reduced_laplacian)
        script = [1] * d  # d 1s
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
                script[bad] += 1
        b = iter(burn)
        s = iter(script)
        bc = {}  # burning config
        bs = {}  # burning script
        for v in self._nonsink_vertices:
            bc[v] = next(b)
            bs[v] = next(s)
        self._burning_config = SandpileConfig(self, bc)
        self._burning_script = SandpileConfig(self, bs)

    def burning_config(self):
        r"""
        The minimal burning configuration.

        OUTPUT:

        dict (configuration)

        EXAMPLES::

            sage: g = {0:{},1:{0:1,3:1,4:1},2:{0:1,3:1,5:1},
            ....:      3:{2:1,5:1},4:{1:1,3:1},5:{2:1,3:1}}
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

        .. NOTE::

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

        OUTPUT: dictionary

        EXAMPLES::

            sage: g = {0:{},1:{0:1,3:1,4:1},2:{0:1,3:1,5:1},
            ....:      3:{2:1,5:1},4:{1:1,3:1},5:{2:1,3:1}}
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

        .. NOTE::

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
        The nonsink vertices.

        OUTPUT: list of vertices

        EXAMPLES::

            sage: s = sandpiles.Grid(2,3)
            sage: s.nonsink_vertices()
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3)]
        """
        return self._nonsink_vertices

    def all_k_config(self, k):
        r"""
        The constant configuration with all values set to `k`.

        INPUT:

        - ``k`` -- integer

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s.all_k_config(7)
            {1: 7, 2: 7, 3: 7}
        """
        return SandpileConfig(self, [k] * (self.num_verts() - 1))

    def zero_config(self):
        r"""
        The all-zero configuration.

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s.zero_config()
            {1: 0, 2: 0, 3: 0}
        """
        return self.all_k_config(0)

    # TODO: cythonize stabilization!
    # The following would presumably be moved to the SandpileConfig class
    # def new_stabilize(self, config):
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
        Compute ``_identity``, the variable holding the identity configuration
        of the sandpile group, when ``identity()`` is first called by a user.

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S._set_identity()
            sage: '_identity' in S.__dict__
            True
        """
        m = self._max_stable
        self._identity = (m & m).dualize() & m

    def identity(self, verbose=True):
        r"""
        The identity configuration.  If ``verbose`` is ``False``, the
        configuration is converted to a list of integers.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: SandpileConfig or a list of integers

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s.identity()
            {1: 2, 2: 2, 3: 0}
            sage: s.identity(False)
            [2, 2, 0]
            sage: s.identity() & s.max_stable() == s.max_stable()
            True
        """
        if verbose:
            return deepcopy(self._identity)
        else:
            return self._identity.values()

    def _set_recurrents(self):
        """
        Compute ``_recurrents``, the variable holding the list of recurrent
        configurations, when ``recurrents()`` is first called by a user.

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s._set_recurrents()
            sage: '_recurrents' in s.__dict__
            True
        """
        if self.name() == 'Complete sandpile graph':
            n = self.num_verts()
            self._recurrents = [SandpileConfig(self, [n - 1 - i for i in p])
                                for p in ParkingFunctions(n - 1)]
        elif self.name() == 'Cycle sandpile graph':
            n = self.num_verts()
            one = [1] * (n - 2)
            self._recurrents = [SandpileConfig(self, [1] * (n - 1))] + [SandpileConfig(self, one[:i] + [0] + one[i:]) for i in range(n - 1)]
        else:
            self._recurrents = []
            active = [self._max_stable]
            while active:
                c = active.pop()
                self._recurrents.append(c)
                for v in self._nonsink_vertices:
                    cnext = deepcopy(c)
                    cnext[v] += 1
                    cnext = ~cnext
                    if (cnext not in active) and (cnext not in self._recurrents):
                        active.append(cnext)
        self._recurrents = self._recurrents

    def recurrents(self, verbose=True):
        r"""
        The recurrent configurations. If ``verbose`` is ``False``, the
        configurations are converted to lists of integers.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: list of recurrent configurations

        EXAMPLES::

            sage: r = Sandpile(graphs.HouseXGraph(),0).recurrents()
            sage: r[:3]
            [{1: 2, 2: 3, 3: 3, 4: 1}, {1: 1, 2: 3, 3: 3, 4: 0}, {1: 1, 2: 3, 3: 3, 4: 1}]
            sage: sandpiles.Complete(4).recurrents(False)                               # needs sage.combinat
            [[2, 2, 2],
             [2, 2, 1],
             [2, 1, 2],
             [1, 2, 2],
             [2, 2, 0],
             [2, 0, 2],
             [0, 2, 2],
             [2, 1, 1],
             [1, 2, 1],
             [1, 1, 2],
             [2, 1, 0],
             [2, 0, 1],
             [1, 2, 0],
             [1, 0, 2],
             [0, 2, 1],
             [0, 1, 2]]
            sage: sandpiles.Cycle(4).recurrents(False)
            [[1, 1, 1], [0, 1, 1], [1, 0, 1], [1, 1, 0]]
        """
        if verbose:
            return deepcopy(self._recurrents)
        else:
            return [r.values() for r in self._recurrents]

    def _set_superstables(self):
        r"""
        Compute ``_superstables``, the variable holding the list of superstable
        configurations, when ``superstables()`` is first called by a user.

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s._set_superstables()
            sage: '_superstables' in s.__dict__
            True
        """
        self._superstables = [c.dualize() for c in self._recurrents]

    def superstables(self, verbose=True):
        r"""
        The superstable configurations.  If ``verbose`` is ``False``, the
        configurations are converted to lists of integers.  Superstables for
        undirected graphs are also known as ``G-parking functions``.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: list of SandpileConfig

        EXAMPLES::

            sage: sp = Sandpile(graphs.HouseXGraph(),0).superstables()
            sage: sp[:3]
            [{1: 0, 2: 0, 3: 0, 4: 0}, {1: 1, 2: 0, 3: 0, 4: 1}, {1: 1, 2: 0, 3: 0, 4: 0}]
            sage: sandpiles.Complete(4).superstables(False)                             # needs sage.combinat
            [[0, 0, 0],
             [0, 0, 1],
             [0, 1, 0],
             [1, 0, 0],
             [0, 0, 2],
             [0, 2, 0],
             [2, 0, 0],
             [0, 1, 1],
             [1, 0, 1],
             [1, 1, 0],
             [0, 1, 2],
             [0, 2, 1],
             [1, 0, 2],
             [1, 2, 0],
             [2, 0, 1],
             [2, 1, 0]]
            sage: sandpiles.Cycle(4).superstables(False)
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
        """
        if verbose:
            return deepcopy(self._superstables)
        else:
            return [s.values() for s in self._superstables]

    def _set_group_gens(self):
        r"""
        A minimal list of generators for the sandpile group.

        EXAMPLES::

            sage: s = sandpiles.Cycle(3)
            sage: s._set_group_gens()
            sage: '_group_gens' in s.__dict__
            True
        """
        D, U, _ = self.reduced_laplacian().transpose().smith_form()
        F = U.inverse()
        self._group_gens = [SandpileConfig(self, [Integer(j) for j in F.column(i)]).equivalent_recurrent()
                            for i in range(F.nrows()) if D[i][i] != 1]

    def group_gens(self, verbose=True):
        r"""
        A minimal list of generators for the sandpile group.  If ``verbose`` is ``False``
        then the generators are represented as lists of integers.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT:

        list of SandpileConfig (or of lists of integers if ``verbose`` is ``False``)

        EXAMPLES::

            sage: s = sandpiles.Cycle(5)
            sage: s.group_gens()
            [{1: 0, 2: 1, 3: 1, 4: 1}]
            sage: s.group_gens()[0].order()
            5
            sage: s = sandpiles.Complete(5)
            sage: s.group_gens(False)
            [[2, 3, 2, 2], [2, 2, 3, 2], [2, 2, 2, 3]]
            sage: [i.order() for i in s.group_gens()]
            [5, 5, 5]
            sage: s.invariant_factors()
            [1, 5, 5, 5]
        """
        if verbose:
            return deepcopy(self._group_gens)
        else:
            return [c.values() for c in self._group_gens]

    def genus(self):
        r"""
        The genus: (# non-loop edges) - (# vertices) + 1.

        This is only defined for undirected graphs.

        OUTPUT: integer

        EXAMPLES::

            sage: sandpiles.Complete(4).genus()
            3
            sage: sandpiles.Cycle(5).genus()
            1
        """
        if self.is_undirected():
            return self.laplacian().trace() / 2 - self.num_verts() + 1
        raise TypeError("the underlying graph must be undirected")

    def is_undirected(self):
        r"""
        Is the underlying graph undirected?  ``True`` if `(u,v)` is and edge if
        and only if `(v,u)` is an edge, each edge with the same weight.

        OUTPUT: boolean

        EXAMPLES::

            sage: sandpiles.Complete(4).is_undirected()
            True
            sage: s = Sandpile({0:[1,2], 1:[0,2], 2:[0]}, 0)
            sage: s.is_undirected()
            False
        """
        return self.laplacian().is_symmetric()

    def _set_min_recurrents(self):
        r"""
        Compute the minimal recurrent elements.  If the underlying graph is
        undirected, these are the recurrent elements of least degree.

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: s._set_min_recurrents()                                               # needs sage.combinat
            sage: '_min_recurrents' in s.__dict__                                       # needs sage.combinat
            True
        """
        if self.is_undirected():
            m = min([r.deg() for r in self.recurrents()])
            rec = [r for r in self.recurrents() if r.deg() == m]
        else:
            rec = list(self.recurrents())
            for r in self.recurrents():
                if exists(rec, lambda x: r > x)[0]:
                    rec.remove(r)
        self._min_recurrents = rec

    def min_recurrents(self, verbose=True):
        r"""
        The minimal recurrent elements.  If the underlying graph is
        undirected, these are the recurrent elements of least degree.
        If ``verbose`` is ``False``, the configurations are converted
        to lists of integers.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: list of SandpileConfig

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s.recurrents(False)
            [[2, 2, 1],
             [2, 2, 0],
             [1, 2, 0],
             [2, 0, 1],
             [0, 2, 1],
             [2, 1, 0],
             [1, 2, 1],
             [2, 1, 1]]
            sage: s.min_recurrents(False)
            [[1, 2, 0], [2, 0, 1], [0, 2, 1], [2, 1, 0]]
            sage: [i.deg() for i in s.recurrents()]
            [5, 4, 3, 3, 3, 3, 4, 4]
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

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: tuple of SandpileConfig

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s.superstables(False)
            [[0, 0, 0],
             [0, 0, 1],
             [1, 0, 1],
             [0, 2, 0],
             [2, 0, 0],
             [0, 1, 1],
             [1, 0, 0],
             [0, 1, 0]]
            sage: s.max_superstables(False)
            [[1, 0, 1], [0, 2, 0], [2, 0, 0], [0, 1, 1]]
            sage: s.h_vector()
            [1, 3, 4]
        """
        result = [r.dualize() for r in self.min_recurrents()]
        if verbose:
            return result
        else:
            return [r.values() for r in result]

    def tutte_polynomial(self):
        r"""
        The Tutte polynomial of the underlying graph.
        Only defined for undirected sandpile graphs.

        OUTPUT: polynomial

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: s.tutte_polynomial()
            x^3 + y^3 + 3*x^2 + 4*x*y + 3*y^2 + 2*x + 2*y
            sage: s.tutte_polynomial().subs(x=1)
            y^3 + 3*y^2 + 6*y + 6
            sage: s.tutte_polynomial().subs(x=1).coefficients() == s.h_vector()         # needs sage.combinat
            True
        """
        if self.is_undirected():
            return Graph(self).tutte_polynomial()
        raise TypeError("the underlying graph must be undirected")

    def _set_avalanche_polynomial(self):
        """
        Compute the avalanche polynomial.  See ``self.avalanche_polynomial`` for details.

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: s._set_avalanche_polynomial()                                         # needs sage.combinat
            sage: '_avalanche_polynomial' in s.__dict__                                 # needs sage.combinat
            True
        """
        n = self.num_verts() - 1
        R = PolynomialRing(QQ, "x", n)
        A = R(0)
        V = []
        for i in range(n):
            c = self.zero_config()
            c[self.nonsink_vertices()[i]] += 1
            V.append(c)
        for r in self.recurrents():
            for i in range(n):
                e = tuple((r + V[i]).stabilize(True)[1].values())
                A += R({e: 1})
        self._avalanche_polynomial = A

    def avalanche_polynomial(self, multivariable=True):
        r"""
        The avalanche polynomial.  See NOTE for details.

        INPUT:

        - ``multivariable`` -- boolean (default: ``True``)

        OUTPUT: polynomial

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: s.avalanche_polynomial()                                              # needs sage.combinat
            9*x0*x1*x2 + 2*x0*x1 + 2*x0*x2 + 2*x1*x2 + 3*x0 + 3*x1 + 3*x2 + 24
            sage: s.avalanche_polynomial(False)                                         # needs sage.combinat
            9*x0^3 + 6*x0^2 + 9*x0 + 24

        .. NOTE::

            For each nonsink vertex `v`, let `x_v` be an indeterminate.
            If `(r,v)` is a pair consisting of a recurrent `r` and nonsink
            vertex `v`, then for each nonsink vertex `w`, let `n_w` be the
            number of times vertex `w` fires in the stabilization of `r + v`.
            Let `M(r,v)` be the monomial `\prod_w x_w^{n_w}`, i.e., the exponent
            records the vector of `n_w` as `w` ranges over the nonsink vertices.
            The avalanche polynomial is then the sum of `M(r,v)` as `r` ranges
            over the recurrents and `v` ranges over the nonsink vertices.  If
            ``multivariable`` is ``False``, then set all the indeterminates equal
            to each other (and, thus, only count the number of vertex firings in the
            stabilizations, forgetting which particular vertices fired).
        """
        if multivariable:
            return deepcopy(self._avalanche_polynomial)
        X = self._avalanche_polynomial.parent().gens()
        return self._avalanche_polynomial.subs({X[i]: X[0]
            for i in range(1, self.num_verts() - 1)})

    def nonspecial_divisors(self, verbose=True):
        r"""
        The nonspecial divisors. Only for undirected graphs.  (See NOTE.)

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: list (of divisors)

        EXAMPLES::

            sage: # needs sage.combinat
            sage: S = sandpiles.Complete(4)
            sage: ns = S.nonspecial_divisors()
            sage: D = ns[0]
            sage: D.values()                                                            # needs sage.symbolic
            [-1, 0, 1, 2]
            sage: D.deg()                                                               # needs sage.symbolic
            2
            sage: [i.effective_div() for i in ns]
            [[], [], [], [], [], []]

        .. NOTE::

            The "nonspecial divisors" are those divisors of degree `g-1` with
            empty linear system.  The term is only defined for undirected graphs.
            Here, `g = |E| - |V| + 1` is the genus of the graph (not counting loops
            as part of `|E|`).  If ``verbose`` is ``False``, the divisors are converted
            to lists of integers.

        .. WARNING::

            The underlying graph must be undirected.
        """
        if self.is_undirected():
            result = []
            for s in self.max_superstables():
                D = dict(s)
                D[self._sink] = -1
                D = SandpileDivisor(self, D)
                result.append(D)
            return result if verbose else [r.values() for r in result]

        raise TypeError("the underlying graph must be undirected")

    def canonical_divisor(self):
        r"""
        The canonical divisor.  This is the divisor with `\deg(v)-2` grains of
        sand on each vertex (not counting loops).  Only for undirected graphs.

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Complete(4)
            sage: S.canonical_divisor()
            {0: 1, 1: 1, 2: 1, 3: 1}
            sage: s = Sandpile({0:[1,1],1:[0,0,1,1,1]},0)
            sage: s.canonical_divisor()  # loops are disregarded
            {0: 0, 1: 0}

        .. WARNING::

            The underlying graph must be undirected.
        """
        if self.is_undirected():
            return SandpileDivisor(self, [self.laplacian()[i][i] - 2 for i in range(self.num_verts())])
        raise TypeError("only for undirected graphs")

    def _set_invariant_factors(self):
        r"""
        Compute the variable holding the elementary divisors of the sandpile
        group when ``invariant_factors()`` is first called by the user.

        EXAMPLES::

            sage: s = sandpiles.Grid(2,2)
            sage: s._set_invariant_factors()
            sage: '_invariant_factors' in s.__dict__
            True
        """
        # Sage seems to have issues with computing the Smith normal form and
        # elementary divisors of a sparse matrix, so we have to convert:
        A = self.reduced_laplacian().dense_matrix()
        self._invariant_factors = A.elementary_divisors()

    def invariant_factors(self):
        r"""
        The invariant factors of the sandpile group.

        OUTPUT: list of integers

        EXAMPLES::

            sage: s = sandpiles.Grid(2,2)
            sage: s.invariant_factors()
            [1, 1, 8, 24]
        """
        return deepcopy(self._invariant_factors)

    def _set_hilbert_function(self):
        """
        Compute the variables holding the Hilbert function of the homogeneous
        homogeneous toppling ideal, the first differences of the Hilbert
        function, and the postulation number for the zero-set of the sandpile
        ideal when any one of these is called by the user.

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: s._set_hilbert_function()
            sage: '_hilbert_function' in s.__dict__
            True
        """
        v = [i.deg() for i in self._superstables]
        self._postulation = max(v)
        self._h_vector = [v.count(i) for i in range(self._postulation + 1)]
        self._hilbert_function = [1]
        for i in range(self._postulation):
            self._hilbert_function.append(self._hilbert_function[i]
                                          + self._h_vector[i + 1])

    def h_vector(self):
        r"""
        The number of superstable configurations in each degree.  Equivalently,
        this is the list of first differences of the Hilbert function of the
        (homogeneous) toppling ideal.

        OUTPUT: list of nonnegative integers

        EXAMPLES::

            sage: s = sandpiles.Grid(2,2)
            sage: s.hilbert_function()
            [1, 5, 15, 35, 66, 106, 146, 178, 192]
            sage: s.h_vector()
            [1, 4, 10, 20, 31, 40, 40, 32, 14]
        """
        return deepcopy(self._h_vector)

    def hilbert_function(self):
        r"""
        The Hilbert function of the homogeneous toppling ideal.

        OUTPUT: list of nonnegative integers

        EXAMPLES::

            sage: s = sandpiles.Wheel(5)
            sage: s.hilbert_function()
            [1, 5, 15, 31, 45]
            sage: s.h_vector()
            [1, 4, 10, 16, 14]
        """
        return deepcopy(self._hilbert_function)

    def postulation(self):
        r"""
        The postulation number of the toppling ideal.  This is the
        largest weight of a superstable configuration of the graph.

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: s.postulation()                                                       # needs sage.combinat
            3
        """
        return self._postulation

    def _set_smith_form(self):
        r"""
        Compute the Smith Normal Form for the transpose of the Laplacian.

        EXAMPLES::

            sage: s = sandpiles.Complete(3)
            sage: s._set_smith_form()
            sage: '_smith_form' in s.__dict__
            True
        """
        self._smith_form = self.laplacian().transpose().smith_form()

    def smith_form(self):
        r"""
        The Smith normal form for the Laplacian.  In detail: a list of integer
        matrices `D, U, V` such that `ULV = D` where `L` is the transpose of the
        Laplacian, `D` is diagonal, and  `U` and `V` are invertible over the
        integers.

        OUTPUT: list of integer matrices

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D,U,V = s.smith_form()
            sage: D
            [1 0 0 0]
            [0 4 0 0]
            [0 0 4 0]
            [0 0 0 0]
            sage: U*s.laplacian()*V == D  # Laplacian symmetric => transpose not necessary
            True
        """
        return deepcopy(self._smith_form)

    def reorder_vertices(self):
        r"""
        A copy of the sandpile with vertex names permuted.

        After reordering, vertex `u` comes before vertex `v` in the
        list of vertices if `u` is closer to the sink.

        OUTPUT: Sandpile

        EXAMPLES::

            sage: S = Sandpile({0:[1], 2:[0,1], 1:[2]})
            sage: S.dict()
            {0: {1: 1}, 1: {2: 1}, 2: {0: 1, 1: 1}}
            sage: T = S.reorder_vertices()

        The vertices 1 and 2 have been swapped::

            sage: T.dict()
            {0: {1: 1}, 1: {0: 1, 2: 1}, 2: {0: 1}}
        """
        # first order the vertices according to their distance to the sink
        distance_to_sink = self.reverse().shortest_path_lengths(self._sink)
        verts = sorted(self, key=lambda v: distance_to_sink[v], reverse=True)
        perm = {v: i for i, v in enumerate(verts)}
        old = self.dict()
        new = {}
        for i in old:
            entry = {}
            for j in old[i]:
                entry[perm[j]] = old[i][j]
            new[perm[i]] = entry
        return Sandpile(new, len(verts) - 1)

    def _set_jacobian_representatives(self):
        r"""
        Find representatives for the elements of the Jacobian group.

        EXAMPLES::

            sage: s = sandpiles.Complete(3)
            sage: s._set_jacobian_representatives()                                     # needs sage.combinat
            sage: '_jacobian_representatives' in s.__dict__                             # needs sage.combinat
            True
        """
        if self.is_undirected():
            easy = True
        else:
            ker = self.laplacian().left_kernel().basis()
            tau = abs(ker[self._sink_ind])
            easy = bool(tau == 1)
        if easy:
            result = []
            for r in self.superstables():
                D = {v: r[v] for v in self._nonsink_vertices}
                D[self._sink] = - r.deg()
                result.append(SandpileDivisor(self, D))
                self._jacobian_representatives = result
        else:
            result = []
            sr = self.superstables()
            order = self.group_order() / tau
            while len(result) < order:
                r = sr.pop()
                active = {v: r[v] for v in self._nonsink_vertices}
                active[self._sink] = -r.deg()
                active = SandpileDivisor(self, active)
                repeated = False
                for D in result:
                    if active.is_linearly_equivalent(D):
                        repeated = True
                        break  # active is repeated in new_result
                if not repeated:
                    result.append(active)
            self._jacobian_representatives = result

    def jacobian_representatives(self, verbose=True):
        r"""
        Representatives for the elements of the Jacobian group. If ``verbose``
        is ``False``, then lists representing the divisors are returned.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: list of SandpileDivisor (or of lists representing divisors)

        EXAMPLES:

        For an undirected graph, divisors of the form ``s - deg(s)*sink`` as
        ``s`` varies over the superstables forms a distinct set of
        representatives for the Jacobian group.::

            sage: s = sandpiles.Complete(3)
            sage: s.superstables(False)                                                 # needs sage.combinat
            [[0, 0], [0, 1], [1, 0]]
            sage: s.jacobian_representatives(False)                                     # needs sage.combinat
            [[0, 0, 0], [-1, 0, 1], [-1, 1, 0]]

        If the graph is directed, the representatives described above may by
        equivalent modulo the rowspan of the Laplacian matrix::

            sage: s = Sandpile({0: {1: 1, 2: 2}, 1: {0: 2, 2: 4}, 2: {0: 4, 1: 2}},0)
            sage: s.group_order()
            28
            sage: s.jacobian_representatives()                                          # needs sage.symbolic
            [{0: -5, 1: 3, 2: 2}, {0: -4, 1: 3, 2: 1}]

        Let `\tau` be the nonnegative generator of the kernel of the transpose of
        the Laplacian, and let `\tau_s` be its sink component, then the sandpile
        group is isomorphic to the direct sum of the cyclic group of order
        `\tau_s` and the Jacobian group.  In the example above, we have::

            sage: s.laplacian().left_kernel()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [14  5  8]

        .. NOTE::

            The Jacobian group is the set of all divisors of degree zero modulo the
            integer rowspan of the Laplacian matrix.
        """
        if verbose:
            return deepcopy(self._jacobian_representatives)
        else:
            return [D.values() for D in self._jacobian_representatives]

    def picard_representatives(self, d, verbose=True):
        r"""
        Representatives of the divisor classes of degree `d` in the Picard group.  (Also
        see the documentation for ``jacobian_representatives``.)

        INPUT:

        - ``d`` -- integer

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: slist of SandpileDivisors (or lists representing divisors)

        EXAMPLES::

            sage: s = sandpiles.Complete(3)
            sage: s.superstables(False)                                                 # needs sage.combinat
            [[0, 0], [0, 1], [1, 0]]
            sage: s.jacobian_representatives(False)                                     # needs sage.combinat
            [[0, 0, 0], [-1, 0, 1], [-1, 1, 0]]
            sage: s.picard_representatives(3,False)                                     # needs sage.combinat
            [[3, 0, 0], [2, 0, 1], [2, 1, 0]]
        """
        D = self.zero_div()
        D[self._sink] = d
        if verbose:
            return [E + D for E in self._jacobian_representatives]
        else:
            return [(E + D).values() for E in self._jacobian_representatives]

    def stable_configs(self, smax=None):
        r"""
        Generator for all stable configurations.  If ``smax`` is provided, then
        the generator gives all stable configurations less than or equal to
        ``smax``.  If ``smax`` does not represent a stable configuration, then each
        component of ``smax`` is replaced by the corresponding component of the
        maximal stable configuration.

        INPUT:

        - ``smax`` -- (optional) SandpileConfig or list representing a SandpileConfig

        OUTPUT: generator for all stable configurations

        EXAMPLES::

            sage: s = sandpiles.Complete(3)
            sage: a = s.stable_configs()
            sage: next(a)                                                               # needs sage.combinat
            {1: 0, 2: 0}
            sage: [i.values() for i in a]                                               # needs sage.combinat
            [[0, 1], [1, 0], [1, 1]]
            sage: b = s.stable_configs([1,0])
            sage: list(b)                                                               # needs sage.combinat
            [{1: 0, 2: 0}, {1: 1, 2: 0}]
        """
        if smax is None:
            smax = self.max_stable().values()
        else:
            c = SandpileConfig(self, smax)
            if not c <= self.max_stable():
                smax = [min(c[v], self.max_stable()[v]) for v in self.nonsink_vertices()]
            else:
                smax = c.values()
        for c in IntegerVectorsIterator(smax):
            yield SandpileConfig(self, c)

    def markov_chain(self, state, distrib=None):
        r"""
        The sandpile Markov chain for configurations or divisors.
        The chain starts at ``state``.  See NOTE for details.

        INPUT:

        - ``state`` -- SandpileConfig, SandpileDivisor, or list representing
          one of these

        - ``distrib`` -- (optional) list of nonnegative numbers summing to 1
          (representing a prob. dist.)

        OUTPUT: generator for Markov chain (see NOTE)

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: m = s.markov_chain([0,0,0])
            sage: next(m)          # random
            {1: 0, 2: 0, 3: 0}
            sage: next(m).values() # random
            [0, 0, 0]
            sage: next(m).values() # random
            [0, 0, 0]
            sage: next(m).values() # random
            [0, 0, 0]
            sage: next(m).values() # random
            [0, 1, 0]
            sage: next(m).values() # random
            [0, 2, 0]
            sage: next(m).values() # random
            [0, 2, 1]
            sage: next(m).values() # random
            [1, 2, 1]
            sage: next(m).values() # random
            [2, 2, 1]
            sage: m = s.markov_chain(s.zero_div(), [0.1,0.1,0.1,0.7])
            sage: next(m).values() # random
            [0, 0, 0, 1]
            sage: next(m).values() # random
            [0, 0, 1, 1]
            sage: next(m).values() # random
            [0, 0, 1, 2]
            sage: next(m).values() # random
            [1, 1, 2, 0]
            sage: next(m).values() # random
            [1, 1, 2, 1]
            sage: next(m).values() # random
            [1, 1, 2, 2]
            sage: next(m).values() # random
            [1, 1, 2, 3]
            sage: next(m).values() # random
            [1, 1, 2, 4]
            sage: next(m).values() # random
            [1, 1, 3, 4]

        .. NOTE::

            The ``closed sandpile Markov chain`` has state space consisting of the configurations
            on a sandpile.  It transitions from a state by choosing a vertex at random
            (according to the probability distribution ``distrib``), dropping a grain of sand at
            that vertex, and stabilizing.  If the chosen vertex is the sink, the chain stays
            at the current state.

            The ``open sandpile Markov chain`` has state space consisting of the recurrent elements,
            i.e., the state space is the sandpile group.  It transitions from the configuration `c`
            by choosing a vertex `v` at random according to ``distrib``.  The next state is the
            stabilization of `c+v`.  If `v` is the sink vertex, then the stabilization of `c+v`
            is defined to be `c`.

            Note that in either case, if ``distrib`` is specified, its length is equal to
            the total number of vertices (including the sink).

        REFERENCES:

        - [Lev2014]_
        """
        st = deepcopy(state)
        V = self.vertices(sort=True)
        n = len(V)
        if isinstance(st, list):
            if len(st) == self.num_verts() - 1:
                st = SandpileConfig(self, st)
            elif len(st) == self.num_verts():
                st = SandpileDivisor(self, st)
            else:
                raise SyntaxError(state)
        if distrib is None:  # default = uniform distribution
            distrib = [QQ.one() / n] * n
        X = GeneralDiscreteDistribution(distrib)
        if isinstance(st, SandpileConfig):
            while True:
                i = X.get_random_element()
                if V[i] != self.sink():
                    st[V[i]] += 1
                    st = st.stabilize()
                yield st
        elif isinstance(st, SandpileDivisor):
            alive = st.is_alive()
            while True:
                i = X.get_random_element()
                st[V[i]] += 1
                if alive:
                    yield st
                else:
                    if st.is_alive():
                        alive = True
                    else:
                        st = st.stabilize()
                    yield st
        else:
            raise SyntaxError(state)

    def _set_stationary_density(self):
        r"""
        Set the stationary density of the sandpile.

        EXAMPLES::

            sage: s = sandpiles.Complete(3)
            sage: s._set_stationary_density()
            sage: '_stationary_density' in s.__dict__
            True
        """
        if self.name() == 'Complete sandpile graph':
            n = Integer(self.num_verts())
            self._stationary_density = (n + QQ.one() / n + sum(falling_factorial(n, i) / n**i for i in range(1, n + 1)) - 3) / 2
        elif self.is_undirected() and '_h_vector' not in self.__dict__:
            t = Graph(self).tutte_polynomial().subs(x=1)
            myR = PolynomialRing(QQ, 'y')
            y = myR.gens()[0]
            t = myR(t)
            dt = derivative(t, y).subs(y=1)
            t = t.subs(y=1)
            self._stationary_density = (self.num_edges()/2 + dt/t)/self.num_verts()
        else:
            sink_deg = self.out_degree(self.sink())
            h = vector(ZZ, self.h_vector())
            m = self.max_stable().deg()
            d = vector(ZZ, range(m, m-len(h), -1))
            self._stationary_density = (h*d/self.group_order() + sink_deg)/self.num_verts()

    def stationary_density(self):
        r"""
        The stationary density of the sandpile.

        OUTPUT: rational number

        EXAMPLES::

            sage: s = sandpiles.Complete(3)
            sage: s.stationary_density()
            10/9

            sage: # needs sage.combinat
            sage: s = Sandpile(digraphs.DeBruijn(2,2),'00')
            sage: s.stationary_density()
            9/8

        .. NOTE::

            The stationary density of a sandpile is the sum `\sum_c (\deg(c) + \deg(s))`
            where `\deg(s)` is the degree of the sink and the sum is over all
            recurrent configurations.

        REFERENCES:

        - [Lev2014]_
        """
        return self._stationary_density

    # ---------------- Functions for divisors ----------------

    def all_k_div(self, k):
        r"""
        The divisor with all values set to `k`.

        INPUT:

        - ``k`` -- integer

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: S.all_k_div(7)
            {0: 7, 1: 7, 2: 7, 3: 7, 4: 7}
        """
        return SandpileDivisor(self, [k]*self.num_verts())

    def zero_div(self):
        r"""
        The all-zero divisor.

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: S.zero_div()
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
        """
        return self.all_k_div(0)

    def _set_betti_complexes(self):
        r"""
        Compute the value return by the ``betti_complexes`` method.

        EXAMPLES::

            sage: S = Sandpile({0:{},1:{0: 1, 2: 1, 3: 4},2:{3: 5},3:{1: 1, 2: 1}},0)
            sage: S._set_betti_complexes()                                              # needs sage.geometry.polyhedron
            sage: '_betti_complexes' in S.__dict__                                      # needs sage.geometry.polyhedron
            True
        """
        results = []
        verts = self.vertices(sort=True)
        r = self.recurrents()
        for D in r:
            d = D.deg()
            # change D to a dict since SandpileConfig will not allow adding a key
            D = dict(D)
            D[self.sink()] = -d
            D = SandpileDivisor(self, D)
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
        The support-complexes with non-trivial homology.  (See NOTE.)

        OUTPUT: list (of pairs [divisors, corresponding simplicial complex])

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: S = Sandpile({0:{},1:{0: 1, 2: 1, 3: 4},2:{3: 5},3:{1: 1, 2: 1}},0)
            sage: p = S.betti_complexes()
            sage: p[0]
            [{0: -8, 1: 5, 2: 4, 3: 1},
             Simplicial complex with vertex set (1, 2, 3) and facets {(3,), (1, 2)}]
            sage: S.resolution()                                                        # needs sage.libs.singular
            'R^1 <-- R^5 <-- R^5 <-- R^1'
            sage: S.betti()                                                             # needs sage.libs.singular
                       0     1     2     3
            ------------------------------
                0:     1     -     -     -
                1:     -     5     5     -
                2:     -     -     -     1
            ------------------------------
            total:     1     5     5     1
            sage: len(p)
            11
            sage: p[0][1].homology()
            {0: Z, 1: 0}
            sage: p[-1][1].homology()
            {0: 0, 1: 0, 2: Z}

        .. NOTE::

            A ``support-complex`` is the simplicial complex formed from the
            supports of the divisors in a linear system.
        """
        return deepcopy(self._betti_complexes)

    #######################################
    #         Algebraic Geometry          #
    #######################################

    def _set_ring(self):
        r"""
        Set up polynomial ring for the sandpile.

        EXAMPLES::

            sage: S = sandpiles.Complete(4)
            sage: S._set_ring()
            sage: '_ring' in S.__dict__
            True
        """
        # first order the vertices according to their distance to the sink
        distance_to_sink = self.reverse().shortest_path_lengths(self._sink)
        verts = sorted(self, key=lambda v: distance_to_sink[v])

        # variable i refers to the i-th vertex in self.vertices(sort=True)
        vertex_to_int = {v: i for i, v in enumerate(self.vertices(sort=True))}
        names = [vertex_to_int[v] for v in reversed(verts)]

        vars = ''
        for i in names:
            vars += 'x' + str(i) + ','
        vars = vars[:-1]
        # create the ring
        self._ring = PolynomialRing(QQ, vars)
        # create the ideal
        gens = []
        for i in self.nonsink_vertices():
            new_gen = 'x' + str(vertex_to_int[i])
            new_gen += '^' + str(self.out_degree(i))
            new_gen += '-'
            for j in self._dict[i]:
                new_gen += 'x' + str(vertex_to_int[j])
                new_gen += '^' + str(self._dict[i][j]) + '*'
            new_gen = new_gen[:-1]
            gens.append(new_gen)
        self._unsaturated_ideal = self._ring.ideal(gens)

    def _set_ideal(self):
        r"""
        Create the saturated lattice ideal for the sandpile.

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S._set_ideal()                                                        # needs sage.libs.singular
            sage: '_ideal' in S.__dict__                                                # needs sage.libs.singular
            True
        """
        from sage.libs.singular.function_factory import ff
        try:
            sat = ff.elim__lib.sat_with_exp
        except NameError:
            sat = ff.elim__lib.sat
        R = self.ring()
        I = self._unsaturated_ideal
        I_sat_gens = sat(I, prod(R.gens()))[0]
        self._ideal = R.ideal(I_sat_gens)

    def unsaturated_ideal(self):
        r"""
        The unsaturated, homogeneous toppling ideal.

        OUTPUT: ideal

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.unsaturated_ideal().gens()
            [x1^3 - x3*x2*x0, x2^3 - x3*x1*x0, x3^2 - x2*x1]
            sage: S.ideal().gens()                                                      # needs sage.libs.singular
            [x2*x1 - x0^2, x3^2 - x0^2, x1^3 - x3*x2*x0,
             x3*x1^2 - x2^2*x0, x2^3 - x3*x1*x0, x3*x2^2 - x1^2*x0]
        """
        return self._unsaturated_ideal

    def ideal(self, gens=False):
        r"""
        The saturated homogeneous toppling ideal.  If ``gens`` is ``True``, the
        generators for the ideal are returned instead.

        INPUT:

        - ``gens`` -- boolean (default: ``False``)

        OUTPUT: ideal or, optionally, the generators of an ideal

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.ideal()                                                             # needs sage.libs.singular
            Ideal (x2*x1 - x0^2, x3^2 - x0^2, x1^3 - x3*x2*x0,
                   x3*x1^2 - x2^2*x0, x2^3 - x3*x1*x0, x3*x2^2 - x1^2*x0)
             of Multivariate Polynomial Ring in x3, x2, x1, x0 over Rational Field
            sage: S.ideal(True)                                                         # needs sage.libs.singular
            [x2*x1 - x0^2, x3^2 - x0^2, x1^3 - x3*x2*x0,
             x3*x1^2 - x2^2*x0, x2^3 - x3*x1*x0, x3*x2^2 - x1^2*x0]
            sage: S.ideal().gens()  # another way to get the generators                 # needs sage.libs.singular
            [x2*x1 - x0^2, x3^2 - x0^2, x1^3 - x3*x2*x0,
             x3*x1^2 - x2^2*x0, x2^3 - x3*x1*x0, x3*x2^2 - x1^2*x0]
        """
        if gens:
            return self._ideal.gens()
        else:
            return self._ideal

    def ring(self):
        r"""
        The ring containing the homogeneous toppling ideal.

        OUTPUT: ring

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.ring()
            Multivariate Polynomial Ring in x3, x2, x1, x0 over Rational Field
            sage: S.ring().gens()
            (x3, x2, x1, x0)

        .. NOTE::

            The indeterminate ``xi`` corresponds to the `i`-th vertex as listed my
            the method ``vertices``. The term-ordering is degrevlex with
            indeterminates ordered according to their distance from the sink (larger
            indeterminates are further from the sink).
        """
        return self._ring

    def _set_resolution(self):
        r"""
        Compute the free resolution of the homogeneous toppling ideal.

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S._set_resolution()                                                   # needs sage.libs.singular
            sage: '_resolution' in S.__dict__                                           # needs sage.libs.singular
            True
        """
        # get the resolution in singular form
        res = self.ideal()._singular_().mres(0)
        # compute the betti numbers
        # self._betti = [1] + [len(res[i])
        #        for i in range(1,len(res)-2)]
        self._betti = [1] + [len(x) for x in res]
        # convert the resolution to a list of Sage poly matrices
        result = []
        zero = self._ring.gens()[0] * 0
        for i in range(1, len(res) + 1):
            syz_mat = []
            new = [res[i][j] for j in range(1, int(res[i].size()) + 1)]
            for j in range(self._betti[i]):
                row = new[j].transpose().sage_matrix(self._ring)
                row = list(row[0])
                if len(row) < self._betti[i-1]:
                    row += [zero]*(self._betti[i-1]-len(row))
                syz_mat.append(row)
            syz_mat = matrix(self._ring, syz_mat).transpose()
            result.append(syz_mat)
        self._resolution = result
        self._singular_resolution = res

    def resolution(self, verbose=False):
        r"""
        A minimal free resolution of the homogeneous toppling ideal.  If
        ``verbose`` is ``True``, then all of the mappings are returned.
        Otherwise, the resolution is summarized.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``)

        OUTPUT: free resolution of the toppling ideal

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: S = Sandpile({0: {}, 1: {0: 1, 2: 1, 3: 4}, 2: {3: 5}, 3: {1: 1, 2: 1}},0)
            sage: S.resolution()  # a Gorenstein sandpile graph
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
            return ' <-- '.join(r)

    def _set_groebner(self):
        r"""
        Compute a Groebner basis for the homogeneous toppling ideal with
        respect to the standard sandpile ordering (see ``ring``).

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S._set_groebner()                                                     # needs sage.libs.singular
            sage: '_groebner' in S.__dict__                                             # needs sage.libs.singular
            True
        """
        self._groebner = self._ideal.groebner_basis()

    def groebner(self):
        r"""
        A Groebner basis for the homogeneous toppling ideal.  It is computed
        with respect to the standard sandpile ordering (see ``ring``).

        OUTPUT: Groebner basis

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.groebner()                                                          # needs sage.libs.singular
            [x3*x2^2 - x1^2*x0, x2^3 - x3*x1*x0, x3*x1^2 - x2^2*x0,
             x1^3 - x3*x2*x0, x3^2 - x0^2, x2*x1 - x0^2]
        """
        return self._groebner

    def betti(self, verbose=True):
        r"""
        The Betti table for the homogeneous toppling ideal.  If
        ``verbose`` is ``True``, it prints the standard Betti table, otherwise,
        it returns a less formatted table.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: Betti numbers for the sandpile

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.betti()                                                             # needs sage.libs.singular
                       0     1     2     3
            ------------------------------
                0:     1     -     -     -
                1:     -     2     -     -
                2:     -     4     9     4
            ------------------------------
            total:     1     6     9     4
            sage: S.betti(False)                                                        # needs sage.libs.singular
            [1, 6, 9, 4]
        """
        if verbose:
            print(singular.eval('print(betti(%s), "betti")' %
                                self._singular_resolution.name()))
        else:
            return self._betti

    def solve(self):
        r"""
        Approximations of the complex affine zeros of the sandpile
        ideal.

        OUTPUT: list of complex numbers

        EXAMPLES::

            sage: S = Sandpile({0: {}, 1: {2: 2}, 2: {0: 4, 1: 1}}, 0)
            sage: Z = S.solve(); Z                                                      # needs sage.libs.singular
            [[-0.707107000000000 + 0.707107000000000*I,
              0.707107000000000 - 0.707107000000000*I],
             [-0.707107000000000 - 0.707107000000000*I,
              0.707107000000000 + 0.707107000000000*I],
             [-I, -I],
             [I, I],
             [0.707107000000000 + 0.707107000000000*I,
              -0.707107000000000 - 0.707107000000000*I],
             [0.707107000000000 - 0.707107000000000*I,
              -0.707107000000000 + 0.707107000000000*I],
             [1, 1],
             [-1, -1]]
            sage: len(Z)                                                                # needs sage.libs.singular
            8
            sage: S.group_order()
            8

        .. NOTE::

            The solutions form a multiplicative group isomorphic to the sandpile
            group.  Generators for this group are given exactly by ``points()``.
        """
        singular.setring(self._ring._singular_())
        v = [singular.var(i) for i in range(1, int(singular.nvars(self._ring)))]
        vars_ = '({})'.format(','.join(str(i) for i in v))

        L = singular.subst(self._ideal,
                           singular.var(singular.nvars(self._ring)), 1)
        _ = singular.ring(0, vars_, 'lp')
        K = singular.fetch(self._ring, L)
        K = singular.groebner(K)
        singular.LIB('solve.lib')
        M = K.solve(5, 1)
        singular.setring(M)
        sol = singular('SOL').sage_structured_str_list()
        return [[SR(j) for j in k] for k in sol[0][0]]

    def _set_points(self):
        r"""
        Generators for the multiplicative group of zeros of the sandpile
        ideal.

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S._set_points()                                                       # needs sage.symbolic
            sage: '_points' in S.__dict__                                               # needs sage.symbolic
            True
        """
        L = self._reduced_laplacian.transpose().dense_matrix()
        n = self.num_verts() - 1
        D, U, V = L.smith_form()
        self._points = []
        one = [1] * n
        twopii = 2 * pi * I
        for k in range(n):
            x = [exp(twopii * U[k, t] / D[k, k]) for t in range(n)]
            if x not in self._points and x != one:
                self._points.append(x)

    def points(self):
        r"""
        Generators for the multiplicative group of zeros of the sandpile
        ideal.

        OUTPUT: list of complex numbers

        EXAMPLES:

        The sandpile group in this example is cyclic, and hence there is a
        single generator for the group of solutions.

        ::

            sage: S = sandpiles.Complete(4)
            sage: S.points()                                                            # needs sage.symbolic
            [[-I, I, 1], [-I, 1, I]]
        """
        return self._points

    # FIX: use the is_symmetric functions for configurations.
    def symmetric_recurrents(self, orbits):
        r"""
        The symmetric recurrent configurations.

        INPUT:

        - ``orbits`` -- list of lists partitioning the vertices

        OUTPUT: list of recurrent configurations

        EXAMPLES::

            sage: S = Sandpile({0: {},
            ....:              1: {0: 1, 2: 1, 3: 1},
            ....:              2: {1: 1, 3: 1, 4: 1},
            ....:              3: {1: 1, 2: 1, 4: 1},
            ....:              4: {2: 1, 3: 1}})
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

        .. NOTE::

            The user is responsible for ensuring that the list of orbits comes from
            a group of symmetries of the underlying graph.
        """
        sym_recurrents = []
        active = [self._max_stable]
        while active:
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


##########################################
#           SandpileConfig Class         #
##########################################
class SandpileConfig(dict):
    r"""
    Class for configurations on a sandpile.
    """
    @staticmethod
    def help(verbose=True):
        r"""
        List of SandpileConfig methods.  If ``verbose``, include short descriptions.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: printed string

        EXAMPLES::

            sage: SandpileConfig.help()
            Shortcuts for SandpileConfig operations:
            ~c    -- stabilize
            c & d -- add and stabilize
            c * c -- add and find equivalent recurrent
            c^k   -- add k times and find equivalent recurrent
                     (taking inverse if k is negative)
            <BLANKLINE>
            For detailed help with any method FOO listed below,
            enter "SandpileConfig.FOO?" or enter "c.FOO?" for any SandpileConfig c.
            <BLANKLINE>
            add_random             -- Add one grain of sand to a random vertex.
            burst_size             -- The burst size of the configuration with respect to the given vertex.
            deg                    -- The degree of the configuration.
            dualize                -- The difference with the maximal stable configuration.
            equivalent_recurrent   -- The recurrent configuration equivalent to the given configuration.
            equivalent_superstable -- The equivalent superstable configuration.
            fire_script            -- Fire the given script.
            fire_unstable          -- Fire all unstable vertices.
            fire_vertex            -- Fire the given vertex.
            help                   -- List of SandpileConfig methods.
            is_recurrent           -- Return whether the configuration is recurrent.
            is_stable              -- Return whether the configuration is stable.
            is_superstable         -- Return whether the configuration is superstable.
            is_symmetric           -- Return whether the configuration is symmetric.
            order                  -- The order of the equivalent recurrent element.
            sandpile               -- The configuration's underlying sandpile.
            show                   -- Show the configuration.
            stabilize              -- The stabilized configuration.
            support                -- The vertices containing sand.
            unstable               -- The unstable vertices.
            values                 -- The values of the configuration as a list.
        """
        _sandpile_help(SandpileConfig, dedent("""\
            Shortcuts for SandpileConfig operations:
            ~c    -- stabilize
            c & d -- add and stabilize
            c * c -- add and find equivalent recurrent
            c^k   -- add k times and find equivalent recurrent
                     (taking inverse if k is negative)

            For detailed help with any method FOO listed below,
            enter "SandpileConfig.FOO?" or enter "c.FOO?" for any SandpileConfig c."""),
                       verbose=verbose)

    def __init__(self, S, c):
        r"""
        Create a configuration on a Sandpile.

        INPUT:

        - ``S`` -- Sandpile

        - ``c`` -- dictionary or list representing a configuration

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = SandpileConfig(S,[1,1,0])
            sage: 3*c
            {1: 3, 2: 3, 3: 0}
            sage: ~(3*c)  # stabilization
            {1: 2, 2: 2, 3: 0}
        """
        if len(c) == S.num_verts()-1:
            if isinstance(c, (dict, SandpileConfig)):
                dict.__init__(self, c)
            elif isinstance(c, list):
                c.reverse()
                config = {}
                for v in S.vertices(sort=True):
                    if v != S.sink():
                        config[v] = c.pop()
                dict.__init__(self, config)
        else:
            raise SyntaxError(c)

        self._sandpile = S
        self._vertices = S.nonsink_vertices()

    def __deepcopy__(self, memo):
        r"""
        Overrides the deepcopy method for dict.

        INPUT:

        - ``memo`` -- (optional) dict

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = SandpileConfig(S,[1,1,0])
            sage: d = deepcopy(c)
            sage: d[1] += 10
            sage: c
            {1: 1, 2: 1, 3: 0}
            sage: d
            {1: 11, 2: 1, 3: 0}
        """
        c = SandpileConfig(self._sandpile, dict(self))
        c.__dict__.update(self.__dict__)
        return c

    def __setitem__(self, key, item):
        r"""
        Overrides the setitem method for dict.

        INPUT:

        - ``key``, ``item`` -- objects

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [4,1])
            sage: c.equivalent_recurrent()
            {1: 1, 2: 1}
            sage: c.__dict__
            {'_equivalent_recurrent': [{1: 1, 2: 1}, {1: 2, 2: 1}],
             '_sandpile': Cycle sandpile graph: 3 vertices, sink = 0,
             '_vertices': [1, 2]}

        .. NOTE::

            In the example, above, changing the value of ``c`` at some vertex makes
            a call to setitem, which resets some of the stored variables for ``c``.
        """
        if key in self:
            dict.__setitem__(self, key, item)
            S = self._sandpile
            V = self._vertices
            self.__dict__ = {'_sandpile': S, '_vertices': V}

    def __getattr__(self, name):
        """
        Set certain variables only when called.

        INPUT:

        - ``name`` -- name of an internal method

        EXAMPLES::

            sage: S = sandpiles.Complete(4)
            sage: C = SandpileConfig(S,[1,1,1])
            sage: C.__getattr__('_deg')
            3
        """
        if name not in self.__dict__:
            if name == '_deg':
                self._set_deg()
                return self.__dict__[name]
            if name == '_stabilize':
                self._set_stabilize()
                return self.__dict__[name]
            if name == '_equivalent_recurrent':
                self._set_equivalent_recurrent()
                return self.__dict__[name]
            if name == '_is_recurrent':
                self._set_is_recurrent()
                return self.__dict__[name]
            if name == '_equivalent_superstable':
                self._set_equivalent_superstable()
                return self.__dict__[name]
            if name == '_is_superstable':
                self._set_is_superstable()
                return self.__dict__[name]
            else:
                raise AttributeError(name)

    def _set_deg(self):
        r"""
        Compute and store the degree of the configuration.

        OUTPUT: integer

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: c._set_deg()
            sage: '_deg' in c.__dict__
            True
        """
        self._deg = sum(self.values())

    def deg(self):
        r"""
        The degree of the configuration.

        OUTPUT: integer

        EXAMPLES::

            sage: S = sandpiles.Complete(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: c.deg()
            3
        """
        return self._deg

    def __add__(self, other):
        r"""
        Addition of configurations.

        INPUT:

        - ``other`` -- SandpileConfig

        OUTPUT: sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [3,2])
            sage: c + d
            {1: 4, 2: 4}
        """
        return SandpileConfig(self.sandpile(),
                              [i + j for i, j in zip(self.values(),
                                                     other.values())])

    def __sub__(self, other):
        r"""
        Subtraction of configurations.

        INPUT:

        - ``other`` -- SandpileConfig

        OUTPUT: sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [3,2])
            sage: c - d
            {1: -2, 2: 0}
        """
        sum = deepcopy(self)
        for v in self:
            sum[v] -= other[v]
        return sum

    def __rsub__(self, other):
        r"""
        Right-side subtraction of configurations.

        INPUT:

        - ``other`` -- SandpileConfig

        OUTPUT: sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = {1: 3, 2: 2}
            sage: d - c
            {1: 2, 2: 0}

        TESTS::

            sage: S = sandpiles.Cycle(3)
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

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: -c
            {1: -1, 2: -2}
        """
        return SandpileConfig(self._sandpile, [-self[v] for v in self._vertices])

    # recurrent addition or multiplication on the right by an integer
    def __mul__(self, other):
        r"""
        If ``other`` is a configuration, the recurrent element equivalent
        to the sum.  If ``other`` is an integer, the sum of configuration with
        itself ``other`` times.

        INPUT:

        - ``other`` -- SandpileConfig or Integer

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
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
            sage: c
            {1: 1, 2: 0, 3: 0}
            sage: c*3
            {1: 3, 2: 0, 3: 0}
        """
        if isinstance(other, SandpileConfig):
            return (self+other).equivalent_recurrent()
        if isinstance(other, Integer):
            return SandpileConfig(self.sandpile(), [other*i for i in self.values()])
        raise TypeError(other)

    def __rmul__(self, other):
        r"""
        The sum of configuration with itself ``other`` times.

        INPUT:

        - ``other`` -- integer

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
            sage: c = SandpileConfig(S,[1,2,3])
            sage: c
            {1: 1, 2: 2, 3: 3}
            sage: 3*c
            {1: 3, 2: 6, 3: 9}
            sage: 3*c == c*3
            True
        """
        return SandpileConfig(self.sandpile(), [other*i for i in self.values()])

    def __le__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at most that of
        ``other``.

        INPUT:

        - ``other`` -- SandpileConfig

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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
        return all(self[v] <= other[v] for v in self._vertices)

    def __lt__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at most that
        of ``other`` and the two configurations are not equal.

        INPUT:

        - ``other`` -- SandpileConfig

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [2,3])
            sage: c < c
            False
            sage: c < d
            True
            sage: d < c
            False
        """
        return self <= other and self != other

    def __ge__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at least that of
        ``other``.

        INPUT:

        - ``other`` -- SandpileConfig

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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
        return all(self[v] >= other[v] for v in self._vertices)

    def __gt__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at least that
        of ``other`` and the two configurations are not equal.

        INPUT:

        - ``other`` -- SandpileConfig

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: d = SandpileConfig(S, [1,3])
            sage: c > c
            False
            sage: d > c
            True
            sage: c > d
            False
        """
        return self >= other and self != other

    # recurrent power
    def __pow__(self, k):
        r"""
        The recurrent element equivalent to the sum of the
        configuration with itself `k` times.  If `k` is negative, do the
        same for the negation of the configuration.  If `k` is zero, return
        the identity of the sandpile group.

        INPUT:

        - ``k`` -- SandpileConfig

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
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
            if k < 0:
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

        - ``other`` -- SandpileConfig

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
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

        OUTPUT: Sandpile

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = S.identity()
            sage: c.sandpile()
            Diamond sandpile graph: 4 vertices, sink = 0
            sage: c.sandpile() == S
            True
        """
        return self._sandpile

    def values(self):
        r"""
        The values of the configuration as a list.

        The list is sorted in the order of the vertices.

        OUTPUT: list of integers

        boolean

        EXAMPLES::

            sage: S = Sandpile({'a':['c','b'], 'b':['c','a'], 'c':['a']},'a')
            sage: c = SandpileConfig(S, {'b':1, 'c':2})
            sage: c
            {'b': 1, 'c': 2}
            sage: c.values()
            [1, 2]
            sage: S.nonsink_vertices()
            ['b', 'c']
        """
        return [self[v] for v in self._vertices]

    def dualize(self):
        r"""
        The difference with the maximal stable configuration.

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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
        Fire the given vertex.

        INPUT:

        - ``v`` -- vertex

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: c = SandpileConfig(S, [1,2])
            sage: c.fire_vertex(2)
            {1: 2, 2: 0}
        """
        c = dict(self)
        c[v] -= self._sandpile.out_degree(v)
        for e in self._sandpile.outgoing_edge_iterator(v):
            if e[1] != self._sandpile.sink():
                c[e[1]] += e[2]
        return SandpileConfig(self._sandpile, c)

    def fire_script(self, sigma):
        r"""
        Fire the given script.  In other words,  fire each vertex the number of
        times indicated by ``sigma``.

        INPUT:

        - ``sigma`` -- SandpileConfig or (list or dict representing a SandpileConfig)

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
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
            for e in self._sandpile.outgoing_edge_iterator(v):
                if e[1] != self._sandpile.sink():
                    c[e[1]] += sigma[i]*e[2]
        return SandpileConfig(self._sandpile, c)

    def unstable(self):
        r"""
        The unstable vertices.

        OUTPUT: list of vertices

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
            sage: c = SandpileConfig(S, [1,2,3])
            sage: c.unstable()
            [2, 3]
        """
        return [v for v in self._vertices if
                self[v] >= self._sandpile.out_degree(v)]

    def fire_unstable(self):
        r"""
        Fire all unstable vertices.

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
            sage: c = SandpileConfig(S, [1,2,3])
            sage: c.fire_unstable()
            {1: 2, 2: 1, 3: 2}
        """
        c = dict(self)
        for v in self.unstable():
            c[v] -= self._sandpile.out_degree(v)
            for e in self._sandpile.outgoing_edge_iterator(v):
                if e[1] != self._sandpile.sink():
                    c[e[1]] += e[2]
        return SandpileConfig(self._sandpile, c)

    def _set_stabilize(self):
        r"""
        Compute the stabilized configuration and its firing vector.

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: c = 2*S.max_stable()
            sage: c._set_stabilize()
            sage: '_stabilize' in c.__dict__
            True
        """
        s = self._sandpile
        c = deepcopy(self)
        firing_vector = s.zero_config()
        unstable = c.unstable()
        while unstable:
            for v in unstable:
                dm = divmod(c[v], s.out_degree(v))
                c[v] = dm[1]
                firing_vector[v] += dm[0]
                for e in s.outgoing_edge_iterator(v):
                    if e[1] != s.sink():
                        c[e[1]] += dm[0] * e[2]
            unstable = c.unstable()
        self._stabilize = [c, firing_vector]

    def stabilize(self, with_firing_vector=False):
        r"""
        The stabilized configuration. Optionally returns the
        corresponding firing vector.

        INPUT:

        - ``with_firing_vector`` -- boolean (default: ``False``)

        OUTPUT: ``SandpileConfig`` or ``[SandpileConfig, firing_vector]``

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: c = 2*S.max_stable()
            sage: c._set_stabilize()
            sage: '_stabilize' in c.__dict__
            True
            sage: S = sandpiles.House()
            sage: c = S.max_stable() + S.identity()
            sage: c.stabilize(True)
            [{1: 1, 2: 2, 3: 2, 4: 1}, {1: 2, 2: 2, 3: 3, 4: 3}]
            sage: S.max_stable() & S.identity() == c.stabilize()
            True
            sage: ~c == c.stabilize()
            True
        """
        if with_firing_vector:
            return self._stabilize
        else:
            return self._stabilize[0]

    def __invert__(self):
        r"""
        The stabilized configuration.

        OUTPUT: ``SandpileConfig``

        Returns the stabilized configuration.

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: c = S.max_stable() + S.identity()
            sage: ~c == c.stabilize()
            True
        """
        return self._stabilize[0]

    def support(self):
        r"""
        The vertices containing sand.

        OUTPUT: list - support of the configuration

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = S.identity()
            sage: c
            {1: 2, 2: 2, 3: 0}
            sage: c.support()
            [1, 2]
        """
        return [i for i in self if self[i] != 0]

    def add_random(self, distrib=None):
        r"""
        Add one grain of sand to a random vertex.  Optionally, a probability
        distribution, ``distrib``, may be placed on the vertices or the nonsink vertices.
        See NOTE for details.

        INPUT:

        - ``distrib`` -- (optional) list of nonnegative numbers summing to 1 (representing a prob. dist.)

        OUTPUT: SandpileConfig

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: c = s.zero_config()
            sage: c.add_random() # random
            {1: 0, 2: 1, 3: 0}
            sage: c
            {1: 0, 2: 0, 3: 0}
            sage: c.add_random([0.1,0.1,0.8]) # random
            {1: 0, 2: 0, 3: 1}
            sage: c.add_random([0.7,0.1,0.1,0.1]) # random
            {1: 0, 2: 0, 3: 0}

        We compute the "sizes" of the avalanches caused by adding random grains
        of sand to the maximal stable configuration on a grid graph.  The
        function ``stabilize()`` returns the firing vector of the
        stabilization, a dictionary whose values say how many times each vertex
        fires in the stabilization.::

            sage: S = sandpiles.Grid(10,10)
            sage: m = S.max_stable()
            sage: a = []
            sage: for i in range(1000):
            ....:     m = m.add_random()
            ....:     m, f = m.stabilize(True)
            ....:     a.append(sum(f.values()))

            sage: # needs sage.plot
            sage: p = list_plot([[log(i + 1), log(a.count(i))]
            ....:                for i in [0..max(a)] if a.count(i)])
            sage: p.axes_labels(['log(N)', 'log(D(N))'])
            sage: t = text("Distribution of avalanche sizes", (2,2), rgbcolor=(1,0,0))
            sage: show(p + t, axes_labels=['log(N)', 'log(D(N))'])      # long time

        .. NOTE::

            If ``distrib`` is ``None``, then the probability is the uniform probability on the nonsink
            vertices.  Otherwise, there are two possibilities:

            (i) the length of ``distrib`` is equal to the number of vertices, and ``distrib`` represents
            a probability distribution on all of the vertices.  In that case, the sink may be chosen
            at random, in which case, the  configuration is unchanged.

            (ii) Otherwise, the length of ``distrib`` must be equal to the number of nonsink vertices,
            and ``distrib`` represents a probability distribution on the nonsink vertices.

        .. WARNING::

            If ``distrib != None``, the user is responsible for assuring the sum of its entries is
            1 and that its length is equal to the number of sink vertices or the number of nonsink vertices.
        """
        c = deepcopy(self)
        n = self._sandpile.num_verts()
        if distrib is None:  # default = uniform distribution on nonsink vertices
            distrib = [QQ.one() / (n - 1)] * (n - 1)
        if len(distrib) == n - 1:  # prob. dist. on nonsink vertices
            X = GeneralDiscreteDistribution(distrib)
            V = self._sandpile.nonsink_vertices()
            c[V[X.get_random_element()]] += 1
        else:  # prob. dist. on all the vertices
            X = GeneralDiscreteDistribution(distrib)
            V = self._sandpile.vertices(sort=True)
            i = X.get_random_element()
            if i != self._sandpile._sink_ind:  # not the sink
                c[V[i]] += 1
        return c

    def order(self):
        r"""
        The order of the equivalent recurrent element.

        OUTPUT: integer

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = SandpileConfig(S,[2,0,1])
            sage: c.order()
            4
            sage: ~(c + c + c + c) == S.identity()
            True
            sage: c = SandpileConfig(S,[1,1,0])
            sage: c.order()
            1
            sage: c.is_recurrent()
            False
            sage: c.equivalent_recurrent() == S.identity()
            True
        """
        v = vector(self.values())
        w = v*self._sandpile.reduced_laplacian().dense_matrix()**(-1)
        return lcm([denominator(i) for i in w])

    def is_stable(self):
        r"""
        Return whether the configuration is stable.

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.max_stable().is_stable()
            True
            sage: (2*S.max_stable()).is_stable()
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
        Set the equivalent recurrent configuration and the corresponding
        firing vector.

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: a = -S.max_stable()
            sage: a._set_equivalent_recurrent()
            sage: '_equivalent_recurrent' in a.__dict__
            True
        """
        old = self
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
        The recurrent configuration equivalent to the given configuration.
        Optionally, return the corresponding firing vector.

        INPUT:

        - ``with_firing_vector`` -- boolean (default: ``False``)

        OUTPUT: SandpileConfig or [SandpileConfig, firing_vector]

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = SandpileConfig(S, [0,0,0])
            sage: c.equivalent_recurrent() == S.identity()
            True
            sage: x = c.equivalent_recurrent(True)
            sage: r = vector([x[0][v] for v in S.nonsink_vertices()])
            sage: f = vector([x[1][v] for v in S.nonsink_vertices()])
            sage: cv = vector(c.values())
            sage: r == cv - f*S.reduced_laplacian()
            True

        .. NOTE::

            Let `L` be the reduced Laplacian, `c` the initial configuration, `r` the
            returned configuration, and `f` the firing vector.  Then `r = c - f\cdot
            L`.
        """
        if with_firing_vector:
            return self._equivalent_recurrent
        else:
            return self._equivalent_recurrent[0]

    def _set_is_recurrent(self):
        r"""
        Compute and store whether the configuration is recurrent.

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = S.max_stable()
            sage: c._set_is_recurrent()
            sage: '_is_recurrent' in c.__dict__
            True
        """
        if '_recurrents' in self._sandpile.__dict__:
            self._is_recurrent = (self in self._sandpile._recurrents)
        elif '_equivalent_recurrent' in self.__dict__:
            self._is_recurrent = (self._equivalent_recurrent == self)
        else:
            # add the burning configuration to config
            b = self._sandpile._burning_config
            c = ~(self + b)
            self._is_recurrent = (c == self)

    def is_recurrent(self):
        r"""
        Return whether the configuration is recurrent.

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.identity().is_recurrent()
            True
            sage: S.zero_config().is_recurrent()
            False
        """
        return self._is_recurrent

    def _set_equivalent_superstable(self):
        r"""
        Set the superstable configuration equivalent to the given
        configuration and its corresponding firing vector.

        OUTPUT: [configuration, firing_vector]

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: m = S.max_stable()
            sage: m._set_equivalent_superstable()
            sage: '_equivalent_superstable' in m.__dict__
            True
        """
        r, fv = self.dualize().equivalent_recurrent(with_firing_vector=True)
        self._equivalent_superstable = [r.dualize(), -fv]

    def equivalent_superstable(self, with_firing_vector=False):
        r"""
        The equivalent superstable configuration. Optionally, return the
        corresponding firing vector.

        INPUT:

        - ``with_firing_vector`` -- boolean (default: ``False``)

        OUTPUT: SandpileConfig or [SandpileConfig, firing_vector]

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: m = S.max_stable()
            sage: m.equivalent_superstable().is_superstable()
            True
            sage: x = m.equivalent_superstable(True)
            sage: s = vector(x[0].values())
            sage: f = vector(x[1].values())
            sage: mv = vector(m.values())
            sage: s == mv - f*S.reduced_laplacian()
            True

        .. NOTE::

            Let `L` be the reduced Laplacian, `c` the initial configuration, `s` the
            returned configuration, and `f` the firing vector.  Then `s = c - f\cdot
            L`.
        """
        if with_firing_vector:
            return self._equivalent_superstable
        else:
            return self._equivalent_superstable[0]

    def _set_is_superstable(self):
        r"""
        Compute and store whether ``config`` is superstable.

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: z = S.zero_config()
            sage: z._set_is_superstable()
            sage: '_is_superstable' in z.__dict__
            True
        """
        if '_superstables' in self._sandpile.__dict__:
            self._is_superstable = (self in self._sandpile._superstables)
        elif '_equivalent_superstable' in self.__dict__:
            self._is_superstable = (self._equivalent_superstable[0] == self)
        else:
            self._is_superstable = self.dualize().is_recurrent()

    def is_superstable(self):
        r"""
        Return whether the configuration is superstable.

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: S.zero_config().is_superstable()
            True
        """
        return self._is_superstable

    def is_symmetric(self, orbits) -> bool:
        r"""
        Return whether the configuration is symmetric.  Return ``True`` if the values of the
        configuration are constant over the vertices in each sublist of
        ``orbits``.

        INPUT:

        - ``orbits`` -- list of lists of vertices

        OUTPUT: boolean

        EXAMPLES::

            sage: S = Sandpile({0: {},
            ....:              1: {0: 1, 2: 1, 3: 1},
            ....:              2: {1: 1, 3: 1, 4: 1},
            ....:              3: {1: 1, 2: 1, 4: 1},
            ....:              4: {2: 1, 3: 1}})
            sage: c = SandpileConfig(S, [1, 2, 2, 3])
            sage: c.is_symmetric([[2,3]])
            True
        """
        return all(len({self[v] for v in x}) <= 1 for x in orbits)

    def burst_size(self, v):
        r"""
        The burst size of the configuration with respect to the given vertex.

        INPUT:

        - ``v`` -- vertex

        OUTPUT: integer

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: [i.burst_size(0) for i in s.recurrents()]
            [1, 1, 1, 1, 1, 1, 1, 1]
            sage: [i.burst_size(1) for i in s.recurrents()]
            [0, 0, 1, 2, 1, 2, 0, 2]

        .. NOTE::

            To define ``c.burst(v)``, if `v` is not the sink, let `c'` be the unique
            recurrent for which the stabilization of `c' + v` is `c`.  The
            burst size is then the amount of sand that goes into the sink during this
            stabilization.  If `v` is the sink, the burst size is defined to be 1.

        REFERENCES:

        - [Lev2014]_
        """
        if v == self.sandpile().sink():
            return 1
        else:
            w = deepcopy(self)
            w[v] -= 1
            w = w.equivalent_recurrent()
            return w.deg() - self.deg() + 1

    def show(self, sink=True, colors=True, heights=False, directed=None, **kwds):
        r"""
        Show the configuration.

        INPUT:

        - ``sink`` -- boolean (default: ``True``); whether to show the sink

        - ``colors`` -- boolean (default: ``True``); whether to color-code the
          amount of sand on each vertex

        - ``heights`` -- boolean (default: ``False``); whether to label each
          vertex with the amount of sand

        - ``directed`` -- (optional) whether to draw directed edges

        - ``kwds`` -- (optional) arguments passed to the show method for Graph

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: c = S.identity()
            sage: c.show()                                                              # needs sage.plot
            sage: c.show(directed=False)                                                # needs sage.plot
            sage: c.show(sink=False, colors=False, heights=True)                        # needs sage.plot
        """
        if directed:
            T = DiGraph(self.sandpile())
        elif directed is False:
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
            for i in T:
                if i == self.sandpile().sink():
                    a[i] = str(i)
                else:
                    a[i] = str(i)+":"+str(self[i])
            T.relabel(a)
        if colors:
            vc = {}  # vertex colors
            r = rainbow(max_height)  # colors # noqa: F821
            for i in range(max_height):
                vc[r[i]] = []
            for i in self.sandpile().nonsink_vertices():
                if heights:
                    vc[r[self[i]]].append(a[i])
                else:
                    vc[r[self[i]]].append(i)
            T.show(vertex_colors=vc, **kwds)
        else:
            T.show(**kwds)


# Note: There ought to be a better way to do this: sage.repl.display is
# intended to help extend pretty-printing capabilities but it still doesn't
# provide an interface to do something as simple as this (in this case we are
# informing IPython that SandpileConfig, being a dict subclass, should be
# pretty-printed in the same way a dict)
pretty.for_type(SandpileConfig, pretty.for_type(dict, None))


###############################################
#           SandpileDivisor Class             #
###############################################

class SandpileDivisor(dict):
    r"""
    Class for divisors on a sandpile.
    """

    @staticmethod
    def help(verbose=True):
        r"""
        List of SandpileDivisor methods.  If ``verbose``, include short descriptions.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: printed string

        EXAMPLES::

            sage: SandpileDivisor.help()
            For detailed help with any method FOO listed below,
            enter "SandpileDivisor.FOO?" or enter "D.FOO?" for any SandpileDivisor D.
            <BLANKLINE>
            Dcomplex               -- The support-complex.
            add_random             -- Add one grain of sand to a random vertex.
            betti                  -- The Betti numbers for the support-complex.
            deg                    -- The degree of the divisor.
            dualize                -- The difference with the maximal stable divisor.
            effective_div          -- All linearly equivalent effective divisors.
            fire_script            -- Fire the given script.
            fire_unstable          -- Fire all unstable vertices.
            fire_vertex            -- Fire the given vertex.
            help                   -- List of SandpileDivisor methods.
            is_alive               -- Return whether the divisor is stabilizable.
            is_linearly_equivalent -- Return whether the given divisor is linearly equivalent.
            is_q_reduced           -- Return whether the divisor is q-reduced.
            is_symmetric           -- Return whether the divisor is symmetric.
            is_weierstrass_pt      -- Return whether the given vertex is a Weierstrass point.
            polytope               -- The polytope determining the complete linear system.
            polytope_integer_pts   -- The integer points inside divisor's polytope.
            q_reduced              -- The linearly equivalent q-reduced divisor.
            rank                   -- The rank of the divisor.
            sandpile               -- The divisor's underlying sandpile.
            show                   -- Show the divisor.
            simulate_threshold     -- The first unstabilizable divisor in the closed Markov chain.
            stabilize              -- The stabilization of the divisor.
            support                -- List of vertices at which the divisor is nonzero.
            unstable               -- The unstable vertices.
            values                 -- The values of the divisor as a list.
            weierstrass_div        -- The Weierstrass divisor.
            weierstrass_gap_seq    -- The Weierstrass gap sequence at the given vertex.
            weierstrass_pts        -- The Weierstrass points (vertices).
            weierstrass_rank_seq   -- The Weierstrass rank sequence at the given vertex.
        """
        _sandpile_help(SandpileDivisor, dedent("""\
            For detailed help with any method FOO listed below,
            enter "SandpileDivisor.FOO?" or enter "D.FOO?" for any SandpileDivisor D."""),
                       verbose=verbose)

    def __init__(self, S, D):
        r"""
        Create a divisor on a Sandpile.

        INPUT:

        - ``S`` -- Sandpile

        - ``D`` -- dictionary or list representing a divisor

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(6)
            sage: D = SandpileDivisor(S,[0,1,0,1,1,3])
            sage: D.support()
            [1, 3, 4, 5]
        """
        if len(D) == S.num_verts():
            if type(D) in [dict, SandpileDivisor, SandpileConfig]:
                dict.__init__(self, dict(D))
            elif isinstance(D, list):
                div = {v: D[i] for i, v in enumerate(S.vertices(sort=True))}
                dict.__init__(self, div)
        else:
            raise SyntaxError(D)

        self._sandpile = S
        self._vertices = S.vertices(sort=True)
        self._weierstrass_rank_seq = {}

    def __deepcopy__(self, memo):
        r"""
        Overrides the deepcopy method for dict.

        INPUT:

        - ``memo`` -- (optional) dictionary

        EXAMPLES::

            sage: S = sandpiles.Cycle(6)
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

        - ``key``, ``item`` -- objects

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S,[0,1,1])
            sage: eff = D.effective_div()                                               # needs sage.geometry.polyhedron
            sage: D.__dict__                                                            # needs sage.geometry.polyhedron
            {'_effective_div': [{0: 0, 1: 1, 2: 1}, {0: 2, 1: 0, 2: 0}],
             '_polytope': A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
             '_polytope_integer_pts': ((0, 0), (1, 1)),
             '_sandpile': Cycle sandpile graph: 3 vertices, sink = 0,
             '_vertices': [0, 1, 2],
             '_weierstrass_rank_seq': {}}
            sage: D[0] += 1
            sage: D.__dict__
            {'_sandpile': Cycle sandpile graph: 3 vertices, sink = 0,
             '_vertices': [0, 1, 2]}

        .. NOTE::

            In the example, above, changing the value of `D` at some vertex makes
            a call to setitem, which resets some of the stored variables for `D`.
        """
        if key in self:
            dict.__setitem__(self, key, item)
            S = self._sandpile
            V = self._vertices
            self.__dict__ = {'_sandpile': S, '_vertices': V}

    def __getattr__(self, name):
        """
        Set certain variables only when called.

        INPUT:

        - ``name`` -- name of an internal method

        EXAMPLES::

            sage: S = sandpiles.Cycle(6)
            sage: D = SandpileDivisor(S,[0,1,0,1,1,3])
            sage: D.__getattr__('_deg')
            6
        """
        if name not in self.__dict__:
            if name == '_deg':
                self._set_deg()
                return self.__dict__[name]
            if name == '_q_reduced':
                self._set_q_reduced()
                return self.__dict__[name]
            if name == '_linear_system':
                self._set_linear_system()
                return self.__dict__[name]
            if name == '_effective_div':
                self._set_effective_div()
                return self.__dict__[name]
            if name == '_polytope':
                self._set_polytope()
                return self.__dict__[name]
            if name == '_polytope_integer_pts':
                self._set_polytope_integer_pts()
                return self.__dict__[name]
            if name == '_rank':
                self._set_rank()
                return self.__dict__[name]
            if name == '_rank_witness':
                self._set_rank(True)
                return self.__dict__[name]
            if name == '_r_of_D':
                self._set_r_of_D()
                return self.__dict__[name]
            if name == '_Dcomplex':
                self._set_Dcomplex()
                return self.__dict__[name]
            if name == '_life':
                self._set_life()
                return self.__dict__[name]
            if name == '_stabilize':
                self._set_stabilize()
                return self.__dict__[name]
            if name == '_weierstrass_pts':
                self._set_weierstrass_pts()
                return self.__dict__[name]
            else:
                raise AttributeError(name)

    def _set_deg(self):
        r"""
        Compute and store the degree of the divisor.

        OUTPUT: integer

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D._set_deg()
            sage: '_deg' in D.__dict__
            True
        """
        self._deg = sum(self.values())

    def deg(self):
        r"""
        The degree of the divisor.

        OUTPUT: integer

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.deg()
            6
        """
        return self._deg

    def __add__(self, other):
        r"""
        Addition of divisors.

        INPUT:

        - ``other`` -- SandpileDivisor

        OUTPUT: sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [3,2,1])
            sage: D + E
            {0: 4, 1: 4, 2: 4}
        """
        return SandpileDivisor(self.sandpile(),
                               [i + j for i, j in zip(self.values(),
                                                      other.values())])

    def __mul__(self, other):
        r"""
        Sum of the divisor with itself ``other`` times.

        INPUT:

        - ``other`` -- integer

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
            sage: D = SandpileDivisor(S,[1,2,3,4])
            sage: D
            {0: 1, 1: 2, 2: 3, 3: 4}
            sage: 3*D
            {0: 3, 1: 6, 2: 9, 3: 12}
            sage: 3*D == D*3
            True
        """
        return SandpileDivisor(self.sandpile(),[i*other for i in self.values()])

    def __rmul__(self, other):
        r"""
        The sum of divisor with itself ``other`` times.

        INPUT:

        - ``other`` -- integer

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
            sage: D = SandpileDivisor(S,[1,2,3,4])
            sage: D
            {0: 1, 1: 2, 2: 3, 3: 4}
            sage: 3*D
            {0: 3, 1: 6, 2: 9, 3: 12}
            sage: 3*D == D*3
            True
        """
        return SandpileDivisor(self.sandpile(),[other*i for i in self.values()])

    def __radd__(self, other):
        r"""
        Right-side addition of divisors.

        INPUT:

        - ``other`` -- SandpileDivisor

        OUTPUT: sum of ``self`` and ``other``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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

        - ``other`` -- SandpileDivisor

        OUTPUT: difference of ``self`` and ``other``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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

        - ``other`` -- SandpileDivisor

        OUTPUT: difference of ``self`` and ``other``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: -D
            {0: -1, 1: -2, 2: -3}
        """
        return SandpileDivisor(self._sandpile, [-self[v] for v in self._vertices])

    def __le__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at most that of
        ``other``.

        INPUT:

        - ``other`` -- SandpileDivisor

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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
        return all(self[v] <= other[v] for v in self._vertices)

    def __lt__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at most that
        of ``other`` and the two divisors are not equal.

        INPUT:

        - ``other`` -- SandpileDivisor

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [2,3,4])
            sage: D < D
            False
            sage: D < E
            True
            sage: E < D
            False
        """
        return self <= other and self != other

    def __ge__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at least that of
        ``other``.

        INPUT:

        - ``other`` -- SandpileDivisor

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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
        return all(self[v] >= other[v] for v in self._vertices)

    def __gt__(self, other):
        r"""
        Return ``True`` if every component of ``self`` is at least that
        of ``other`` and the two divisors are not equal.

        INPUT:

        - ``other`` -- SandpileDivisor

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: E = SandpileDivisor(S, [1,3,4])
            sage: D > D
            False
            sage: E > D
            True
            sage: D > E
            False
        """
        return self >= other and self != other

    def sandpile(self):
        r"""
        The divisor's underlying sandpile.

        OUTPUT: Sandpile

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: D = SandpileDivisor(S,[1,-2,0,3])
            sage: D.sandpile()
            Diamond sandpile graph: 4 vertices, sink = 0
            sage: D.sandpile() == S
            True
        """
        return self._sandpile

    def values(self):
        r"""
        The values of the divisor as a list.

        The list is sorted in the order of the vertices.

        OUTPUT: list of integers

        boolean

        EXAMPLES::

            sage: S = Sandpile({'a':['c','b'], 'b':['c','a'], 'c':['a']},'a')
            sage: D = SandpileDivisor(S, {'a':0, 'b':1, 'c':2})
            sage: D
            {'a': 0, 'b': 1, 'c': 2}
            sage: D.values()
            [0, 1, 2]
            sage: S.vertices(sort=True)
            ['a', 'b', 'c']
        """
        return [self[v] for v in self._vertices]

    def dualize(self):
        r"""
        The difference with the maximal stable divisor.

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.dualize()
            {0: 0, 1: -1, 2: -2}
            sage: S.max_stable_div() - D == D.dualize()
            True
        """
        return self._sandpile.max_stable_div() - self

    def fire_vertex(self, v):
        r"""
        Fire the given vertex.

        INPUT:

        - ``v`` -- vertex

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.fire_vertex(1)
            {0: 2, 1: 0, 2: 4}
        """
        D = dict(self)
        D[v] -= self._sandpile.out_degree(v)
        for e in self._sandpile.outgoing_edge_iterator(v):
            D[e[1]] += e[2]
        return SandpileDivisor(self._sandpile,D)

    def fire_script(self, sigma):
        r"""
        Fire the given script.  In other words, fire each vertex the number of
        times indicated by ``sigma``.

        INPUT:

        - ``sigma`` -- SandpileDivisor or (list or dict representing a SandpileDivisor)

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
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
            for e in self._sandpile.outgoing_edge_iterator(v):
                D[e[1]] += sigma[i] * e[2]
        return SandpileDivisor(self._sandpile, D)

    def unstable(self):
        r"""
        The unstable vertices.

        OUTPUT: list of vertices

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.unstable()
            [1, 2]
        """
        return [v for v in self._vertices if
                self[v] >= self._sandpile.out_degree(v)]

    def fire_unstable(self):
        r"""
        Fire all unstable vertices.

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [1,2,3])
            sage: D.fire_unstable()
            {0: 3, 1: 1, 2: 2}
        """
        D = dict(self)
        for v in self.unstable():
            D[v] -= self._sandpile.out_degree(v)
            for e in self._sandpile.outgoing_edge_iterator(v):
                D[e[1]] += e[2]
        return SandpileDivisor(self._sandpile,D)

    def _set_q_reduced(self):
        r"""
        The linearly equivalent `q`-reduced divisor.

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[2,-3,2,0])
            sage: D._set_q_reduced()
            sage: '_q_reduced' in D.__dict__
            True
        """
        S = self.sandpile()
        c = SandpileConfig(S,[self[i] for i in S.nonsink_vertices()])
        c = c.equivalent_superstable()
        D = {v:c[v] for v in S.nonsink_vertices()}
        D[S.sink()] = self.deg() - c.deg()
        self._q_reduced = SandpileDivisor(S,D)

    def q_reduced(self, verbose=True):
        r"""
        The linearly equivalent `q`-reduced divisor.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: SandpileDivisor or list representing SandpileDivisor

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[2,-3,2,0])
            sage: D.q_reduced()
            {0: -2, 1: 1, 2: 2, 3: 0}
            sage: D.q_reduced(False)
            [-2, 1, 2, 0]

        .. NOTE::

            The divisor `D` is `qreduced if `D = c + kq` where `c`
            is superstable, `k` is an integer, and `q` is the sink.
        """
        if verbose:
            return deepcopy(self._q_reduced)
        else:
            return self._q_reduced.values()

    def is_q_reduced(self):
        r"""
        Return whether the divisor is `q`-reduced.  This would mean that `self = c + kq` where
        `c` is superstable, `k` is an integer, and `q` is the sink vertex.

        OUTPUT: boolean

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[2,-3,2,0])
            sage: D.is_q_reduced()
            False
            sage: SandpileDivisor(s,[10,0,1,2]).is_q_reduced()
            True

        For undirected or, more generally, Eulerian graphs, `q`-reduced divisors are
        linearly equivalent if and only if they are equal.  The same does not hold for
        general directed graphs:

        ::

            sage: s = Sandpile({0:[1],1:[1,1]})
            sage: D = SandpileDivisor(s,[-1,1])
            sage: Z = s.zero_div()
            sage: D.is_q_reduced()
            True
            sage: Z.is_q_reduced()
            True
            sage: D == Z
            False
            sage: D.is_linearly_equivalent(Z)
            True
        """
        S = self.sandpile()
        c = SandpileConfig(S,[self[v] for v in S.nonsink_vertices()])
        return c.is_superstable()

    def is_linearly_equivalent(self, D, with_firing_vector=False):
        r"""
        Return whether the given divisor is linearly equivalent.  Optionally, returns the
        firing vector.  (See NOTE.)

        INPUT:

        - ``D`` -- SandpileDivisor or list, tuple, etc. representing a divisor

        - ``with_firing_vector`` -- boolean (default: ``False``)

        OUTPUT: boolean or integer vector

        EXAMPLES::

            sage: s = sandpiles.Complete(3)
            sage: D = SandpileDivisor(s,[2,0,0])
            sage: D.is_linearly_equivalent([0,1,1])
            True
            sage: D.is_linearly_equivalent([0,1,1],True)
            (0, -1, -1)
            sage: v = vector(D.is_linearly_equivalent([0,1,1],True))
            sage: vector(D.values()) - s.laplacian()*v
            (0, 1, 1)
            sage: D.is_linearly_equivalent([0,0,0])
            False
            sage: D.is_linearly_equivalent([0,0,0],True)
            ()

        .. NOTE::

            - If ``with_firing_vector`` is ``False``, returns either ``True`` or ``False``.

            - If ``with_firing_vector`` is ``True`` then: (i) if ``self`` is linearly
              equivalent to `D`, returns a vector `v` such that ``self - v*self.laplacian().transpose() = D``.
              Otherwise, (ii) if ``self`` is not linearly equivalent to `D`, the output is the empty vector, ``()``.
        """
        # First try to convert D into a vector.
        v = vector(self.values())
        if isinstance(D,SandpileDivisor):
            w = vector(D.values())
        else:
            w = vector(D)
        # Now test for linear equivalence and find firing vector
        D,U,V = self.sandpile()._smith_form
        b = v - w
        ub = U*b
        if ub[-1] != 0:
            if with_firing_vector:
                return vector([])
            else:
                return False
        else:
            try:
                x = vector(ZZ,[ub[i]/D[i][i] for i in range(D.nrows()-1)]+[0])
                if with_firing_vector:
                    return V*x
                else:
                    return True
            except Exception:
                if with_firing_vector:
                    return vector([])
                else:
                    return False

    def simulate_threshold(self, distrib=None):
        r"""
        The first unstabilizable divisor in the closed Markov chain.
        (See NOTE.)

        INPUT:

        - ``distrib`` -- (optional)  list of nonnegative numbers representing a probability distribution on the vertices

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = s.zero_div()
            sage: D.simulate_threshold()  # random
            {0: 2, 1: 3, 2: 1, 3: 2}
            sage: n(mean([D.simulate_threshold().deg() for _ in range(10)]))  # random
            7.10000000000000
            sage: n(s.stationary_density()*s.num_verts())
            6.93750000000000

        .. NOTE::

            Starting at ``self``, repeatedly choose a vertex and add a grain of
            sand to it.  Return the first unstabilizable divisor that is
            reached.  Also see the ``markov_chain`` method for the underlying
            sandpile.
        """
        E = deepcopy(self)
        S = E.sandpile()
        V = S.vertices(sort=True)
        n = S.num_verts()
        if distrib is None:  # default = uniform distribution
            distrib = [QQ.one() / n] * n
        X = GeneralDiscreteDistribution(distrib)
        while not E.is_alive():
            E = E.stabilize()
            i = X.get_random_element()
            E[V[i]] += 1
        return E

    def _set_linear_system(self):
        r"""
        Compute and store the complete linear system of a divisor.

        OUTPUT:

        dictionary; ``{num_homog: int, homog:list, num_inhomog:int, inhomog:list}``

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [0,1,1])
            sage: D._set_linear_system() # known bug (won't fix due to deprecation optional - 4ti2)

        .. WARNING::

            This method requires 4ti2.
        """
        L = self._sandpile._laplacian.transpose()
        n = self._sandpile.num_verts()

        # temporary file names
        lin_sys = tmp_filename()
        lin_sys_mat = lin_sys + '.mat'
        lin_sys_rel = lin_sys + '.rel'
        lin_sys_rhs = lin_sys + '.rhs'
        lin_sys_sign = lin_sys + '.sign'
        lin_sys_zhom = lin_sys + '.zhom'
        lin_sys_zinhom = lin_sys + '.zinhom'
        lin_sys_log = lin_sys + '.log'

        with open(lin_sys_mat, 'w') as mat_file:
            mat_file.write(str(n)+' ')
            mat_file.write(str(n)+'\n')
            for r in L:
                mat_file.write(''.join(map(str,r)))
                mat_file.write('\n')
        # relations file
        with open(lin_sys_rel, 'w') as rel_file:
            rel_file.write('1 ')
            rel_file.write(str(n)+'\n')
            rel_file.write('>'*n)
            rel_file.write('\n')
        # right-hand side file
        with open(lin_sys_rhs, 'w') as rhs_file:
            rhs_file.write('1 ')
            rhs_file.write(str(n)+'\n')
            rhs_file.write(''.join(str(-i) for i in self.values()))
            rhs_file.write('\n')
        # sign file
        with open(lin_sys_sign, 'w') as sign_file:
            sign_file.write('1 ')
            sign_file.write(str(n)+'\n')
            """
            Conjecture: taking only 1s just below is OK, i.e.,
            looking for solutions with nonnegative entries.  The
            Laplacian has kernel of dimension 1, generated by a
            nonnegative vector.  I would like to say that translating
            by this vector, we transform any solution into a
            nonnegative solution.  What if the vector in the kernel
            does not have full support though?
            """
            sign_file.write('2'*n)  # so maybe a 1 could go here
            sign_file.write('\n')
        # compute
        import os
        import shlex
        try:
            path_to_zsolve = FourTi2Executable('zsolve').absolute_filename()
            os.system(shlex.quote(path_to_zsolve) + ' -q ' + lin_sys + ' > ' + lin_sys_log)
            # process the results
            zhom_file = open(lin_sys_zhom)
        except OSError:
            print("""
                 **********************************
                 *** This method requires 4ti2. ***
                 **********************************
            """)
            return
        # first, the cone generators (the homogeneous points)
        a = zhom_file.read()
        zhom_file.close()
        a = a.split('\n')
        # a starts with two numbers. We are interested in the first one
        num_homog = int(a[0].split()[0])
        homog = [map(int, i.split()) for i in a[1:-1]]
        # second, the inhomogeneous points
        with open(lin_sys_zinhom) as zinhom_file:
            b = zinhom_file.read()
        b = b.split('\n')
        num_inhomog = int(b[0].split()[0])
        inhomog = [map(int, i.split()) for i in b[1:-1]]
        self._linear_system = {'num_homog': num_homog, 'homog': homog,
                               'num_inhomog': num_inhomog, 'inhomog': inhomog}

    def _set_polytope(self):
        r"""
        Compute the polyhedron determining the linear system for D.

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[4,2,0,0])
            sage: D._set_polytope()                                                     # needs sage.geometry.polyhedron
            sage: '_polytope' in D.__dict__                                             # needs sage.geometry.polyhedron
            True
        """
        S = self.sandpile()
        myL = S.laplacian().transpose().delete_columns([S._sink_ind])
        my_ieqs = [[self[v]] + list(-myL[i]) for i,v in enumerate(S.vertices(sort=True))]
        self._polytope = Polyhedron(ieqs=my_ieqs)

    def polytope(self):
        r"""
        The polytope determining the complete linear system.

        OUTPUT: polytope

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[4,2,0,0])
            sage: p = D.polytope()
            sage: p.inequalities()
            (An inequality (-3, 1, 1) x + 2 >= 0,
             An inequality (1, 1, 1) x + 4 >= 0,
             An inequality (1, -3, 1) x + 0 >= 0,
             An inequality (1, 1, -3) x + 0 >= 0)
            sage: D = SandpileDivisor(s,[-1,0,0,0])
            sage: D.polytope()
            The empty polyhedron in QQ^3

        .. NOTE::

            For a divisor `D`, this is the intersection of (i) the polyhedron
            determined by the system of inequalities `L^t x \leq D` where `L^t`
            is the transpose of the Laplacian with (ii) the hyperplane
            `x_{\mathrm{sink\_vertex}} = 0`. The polytope is thought of as sitting in
            `(n-1)`-dimensional Euclidean space where `n` is the number of
            vertices.
        """
        return deepcopy(self._polytope)

    def _set_polytope_integer_pts(self):
        r"""
        Record the integer lattice points inside the polytope determining the
        complete linear system (see the documentation for ``polytope``).

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[4,2,0,0])
            sage: D._set_polytope_integer_pts()                                         # needs sage.geometry.polyhedron
            sage: '_polytope_integer_pts' in D.__dict__                                 # needs sage.geometry.polyhedron
            True
        """
        self._polytope_integer_pts = self._polytope.integral_points()

    def polytope_integer_pts(self):
        r"""
        The integer points inside divisor's polytope.  The polytope referred to
        here is the one determining the divisor's complete linear system (see the
        documentation for ``polytope``).

        OUTPUT: tuple of integer vectors

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[4,2,0,0])
            sage: sorted(D.polytope_integer_pts())                                      # needs sage.geometry.polyhedron
            [(-2, -1, -1),
             (-1, -2, -1),
             (-1, -1, -2),
             (-1, -1, -1),
             (0, -1, -1),
             (0, 0, 0)]
            sage: D = SandpileDivisor(s,[-1,0,0,0])
            sage: D.polytope_integer_pts()                                              # needs sage.geometry.polyhedron
            ()
        """
        return deepcopy(self._polytope_integer_pts)

    def _set_effective_div(self):
        r"""
        Compute all of the linearly equivalent effective divisors linearly.

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[4,2,0,0])
            sage: D._set_effective_div()                                                # needs sage.geometry.polyhedron
            sage: '_effective_div' in D.__dict__                                        # needs sage.geometry.polyhedron
            True
        """
        S = self.sandpile()
        myL = S.laplacian().transpose().delete_columns([S._sink_ind])
        dv = vector(ZZ, self.values())
        self._effective_div = [SandpileDivisor(S,list(dv - myL*i))
                               for i in self._polytope_integer_pts]

    def effective_div(self, verbose=True, with_firing_vectors=False):
        r"""
        All linearly equivalent effective divisors.  If ``verbose``
        is ``False``, the divisors are converted to lists of integers.
        If ``with_firing_vectors`` is ``True`` then a list of firing vectors
        is also given, each of which prescribes the vertices to be fired
        in order to obtain an effective divisor.

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        - ``with_firing_vectors`` -- boolean (default: ``False``)

        OUTPUT: list (of divisors)

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[4,2,0,0])
            sage: sorted(D.effective_div(), key=str)
            [{0: 0, 1: 2, 2: 0, 3: 4},
             {0: 0, 1: 2, 2: 4, 3: 0},
             {0: 0, 1: 6, 2: 0, 3: 0},
             {0: 1, 1: 3, 2: 1, 3: 1},
             {0: 2, 1: 0, 2: 2, 3: 2},
             {0: 4, 1: 2, 2: 0, 3: 0}]
            sage: sorted(D.effective_div(False))
            [[0, 2, 0, 4],
             [0, 2, 4, 0],
             [0, 6, 0, 0],
             [1, 3, 1, 1],
             [2, 0, 2, 2],
             [4, 2, 0, 0]]
            sage: sorted(D.effective_div(with_firing_vectors=True), key=str)
            [({0: 0, 1: 2, 2: 0, 3: 4}, (0, -1, -1, -2)),
             ({0: 0, 1: 2, 2: 4, 3: 0}, (0, -1, -2, -1)),
             ({0: 0, 1: 6, 2: 0, 3: 0}, (0, -2, -1, -1)),
             ({0: 1, 1: 3, 2: 1, 3: 1}, (0, -1, -1, -1)),
             ({0: 2, 1: 0, 2: 2, 3: 2}, (0, 0, -1, -1)),
             ({0: 4, 1: 2, 2: 0, 3: 0}, (0, 0, 0, 0))]
            sage: a = _[2]
            sage: a[0].values()
            [0, 6, 0, 0]
            sage: vector(D.values()) - s.laplacian()*a[1]
            (0, 6, 0, 0)
            sage: sorted(D.effective_div(False, True))
            [([0, 2, 0, 4], (0, -1, -1, -2)),
             ([0, 2, 4, 0], (0, -1, -2, -1)),
             ([0, 6, 0, 0], (0, -2, -1, -1)),
             ([1, 3, 1, 1], (0, -1, -1, -1)),
             ([2, 0, 2, 2], (0, 0, -1, -1)),
             ([4, 2, 0, 0], (0, 0, 0, 0))]
            sage: D = SandpileDivisor(s,[-1,0,0,0])
            sage: D.effective_div(False,True)
            []
        """
        S = self.sandpile()
        eff = deepcopy(self._effective_div)
        if with_firing_vectors:
            fv = [vector(list(i)[:S._sink_ind] + [0] + list(i)[S._sink_ind:]) for i in self._polytope_integer_pts]
        if verbose and with_firing_vectors:
            return list(zip(eff, fv))
        elif verbose:   # verbose without firing vectors
            return eff
        elif with_firing_vectors:  # not verbose but with firing vectors
            return list(zip([i.values() for i in eff], fv))
        else:  # not verbose, no firing vectors
            return [i.values() for i in eff]

    def _set_rank(self, set_witness=False):
        r"""
        Find the rank of the divisor `D` and an effective divisor `E` such that
        `D - E` is unwinnable, i.e., has an empty complete linear system.  If
        Riemann-Roch applies, ``verbose`` is ``False``, and the degree of `D` is greater
        than `2g-2` (`g = ` genus), then the rank is `\deg(D) - g`.  In that case,
        the divisor `E` is not calculated.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[4,2,0,0])
            sage: D._set_rank()
            sage: '_rank' in D.__dict__
            True
            sage: '_rank_witness' in D.__dict__
            False
            sage: D._set_rank(True)                                                     # needs sage.geometry.polyhedron
            sage: '_rank_witness' in D.__dict__                                         # needs sage.geometry.polyhedron
            True
            sage: D = SandpileDivisor(s,[1,0,0,0])
            sage: D._set_rank()
            sage: '_rank' in D.__dict__
            True
            sage: '_rank_witness' in D.__dict__
            False
        """
        S = self.sandpile()
        # If undirected and D has high degree, use Riemann-Roch.
        if S.is_undirected() and not set_witness:
            # We've been careful about loops
            g = sum(S.laplacian().diagonal()) / 2 - S.num_verts() + 1
            if self.deg() > 2 * g - 2:
                self._rank = self.deg() - g
                return  # return early
        # If S is a complete sandpile graph and a witness is not needed, use
        # the Cori-Le Borgne algorithm
        if S.name() == 'Complete sandpile graph' and not set_witness:
            # Cori-LeBorgne algorithm
            n = S.num_verts()
            rk = -1
            E = self.q_reduced()
            k = E[S.sink()]
            c = [E[v] for v in S.nonsink_vertices()]
            c.sort()
            while k >= 0:
                rk += 1
                try:
                    d = next(i for i,j in enumerate(c) if i == j and i != 0)
                except Exception:
                    d = n - 1
                k = k - d
                if k >= 0:
                    c[0] = n - 1 - d
                    b1 = [c[i] + n - d for i in range(1,d)]
                    b2 = [c[i] - d for i in range(d,n-1)]
                    c = b2 + [c[0]] + b1
            self._rank = rk
        # All other cases.
        else:
            rk = -1
            while True:
                for e in integer_vectors_nk_fast_iter(rk + 1, S.num_verts()):
                    E = SandpileDivisor(S, e)
                    if not (self - E).effective_div():
                        self._rank = rk
                        self._rank_witness = E
                        return
                rk += 1

    def rank(self, with_witness=False):
        r"""
        The rank of the divisor.  Optionally returns an effective divisor `E` such
        that `D - E` is not winnable (has an empty complete linear system).

        INPUT:

        - ``with_witness`` -- boolean (default: ``False``)

        OUTPUT: integer or (integer, SandpileDivisor)

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: S = sandpiles.Complete(4)
            sage: D = SandpileDivisor(S,[4,2,0,0])
            sage: D.rank()
            3
            sage: D.rank(True)
            (3, {0: 3, 1: 0, 2: 1, 3: 0})
            sage: E = _[1]
            sage: (D - E).rank()                                                        # needs sage.rings.number_field
            -1

         Riemann-Roch theorem::

            sage: # needs sage.geometry.polyhedron
            sage: D.rank() - (S.canonical_divisor()-D).rank() == D.deg() + 1 - S.genus()
            True

         Riemann-Roch theorem::

            sage: # needs sage.geometry.polyhedron
            sage: D.rank() - (S.canonical_divisor()-D).rank() == D.deg() + 1 - S.genus()
            True
            sage: S = Sandpile({0:[1,1,1,2],1:[0,0,0,1,1,1,2,2],2:[2,2,1,1,0]},0)  # multigraph with loops
            sage: D = SandpileDivisor(S,[4,2,0])
            sage: D.rank(True)
            (2, {0: 1, 1: 1, 2: 1})
            sage: S = Sandpile({0:[1,2], 1:[0,2,2], 2: [0,1]},0) # directed graph
            sage: S.is_undirected()
            False
            sage: D = SandpileDivisor(S,[0,2,0])
            sage: D.effective_div()
            [{0: 0, 1: 2, 2: 0}, {0: 2, 1: 0, 2: 0}]
            sage: D.rank(True)
            (0, {0: 0, 1: 0, 2: 1})
            sage: E = D.rank(True)[1]
            sage: (D - E).effective_div()                                               # needs sage.rings.number_field
            []

        .. NOTE::

            The rank of a divisor `D` is -1 if `D` is not linearly equivalent to an effective divisor
            (i.e., the dollar game represented by `D` is unwinnable).  Otherwise, the rank of `D` is
            the largest integer `r` such that `D - E` is linearly equivalent to an effective divisor
            for all effective divisors `E` with `\deg(E) = r`.
        """
        if with_witness:
            return (self._rank, deepcopy(self._rank_witness))
        else:
            return self._rank

    def _set_r_of_D(self, verbose=False):
        r"""
        Compute `r(D)` and an effective divisor `F` such that `|D - F|` is
        empty.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: S = sandpiles.Cycle(6)
            sage: D = SandpileDivisor(S, [0,0,0,0,0,4])         # optional - 4ti2
            sage: D._set_r_of_D()                               # optional - 4ti2
        """
        eff = self.effective_div()
        n = self._sandpile.num_verts()
        r = -1
        if not eff:
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
                    print(r)
                new_level = []
                for v in level:
                    for i in range(n):
                        w = v + e[i]
                        if w not in new_level:
                            new_level.append(w)
                            C = d - w
                            C = SandpileDivisor(self._sandpile,list(C))
                            eff = C.effective_div()
                            if not eff:
                                self._r_of_D = (r, SandpileDivisor(self._sandpile,list(w)))
                                return
                level = new_level

    def weierstrass_rank_seq(self, v='sink'):
        r"""
        The Weierstrass rank sequence at the given vertex.  Computes the rank of
        the divisor `D - nv` starting with `n=0` and ending when the rank is
        `-1`.

        INPUT:

        - ``v`` -- (default: ``sink``) vertex

        OUTPUT: tuple of int

        EXAMPLES::

            sage: s = sandpiles.House()
            sage: K = s.canonical_divisor()
            sage: [K.weierstrass_rank_seq(v) for v in s.vertices(sort=True)]            # needs sage.geometry.polyhedron
            [(1, 0, -1), (1, 0, -1), (1, 0, -1), (1, 0, -1), (1, 0, 0, -1)]
        """
        s = self.sandpile()
        if v == 'sink':
            v = s.sink()
        try:
            seq = self._weierstrass_rank_seq[v]
        except Exception:
            D = deepcopy(self)
            verts = s.vertices(sort=True)
            Ei = s.zero_div()
            Ei[verts.index(v)] = 1
            Ei = SandpileDivisor(s,Ei)
            r = D.rank()
            seq = [r]
            while r != -1:
                D = D - Ei
                r = D.rank()
                seq.append(r)
            self._weierstrass_rank_seq[v] = seq
        return tuple(seq)

    def weierstrass_gap_seq(self, v='sink', weight=True):
        r"""
        The Weierstrass gap sequence at the given vertex.  If ``weight`` is
        ``True``, then also compute the weight of each gap value.

        INPUT:

        - ``v`` -- (default: ``sink``) vertex

        - ``weight`` -- boolean (default: ``True``)

        OUTPUT: list or (list of list) of integers

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: s = sandpiles.Cycle(4)
            sage: D = SandpileDivisor(s,[2,0,0,0])
            sage: [D.weierstrass_gap_seq(v,False) for v in s.vertices(sort=True)]
            [(1, 3), (1, 2), (1, 3), (1, 2)]
            sage: [D.weierstrass_gap_seq(v) for v in s.vertices(sort=True)]
            [((1, 3), 1), ((1, 2), 0), ((1, 3), 1), ((1, 2), 0)]
            sage: D.weierstrass_gap_seq()   # gap sequence at sink vertex, 0
            ((1, 3), 1)
            sage: D.weierstrass_rank_seq()  # rank sequence at the sink vertex
            (1, 0, 0, -1)

        .. NOTE::

            The integer `k` is a Weierstrass gap for the divisor `D` at vertex `v` if the rank
            of `D - (k-1)v` does not equal the rank of `D - kv`.  Let `r` be the rank of `D` and
            let `k_i` be the `i`-th gap at `v`.  The Weierstrass weight of `v` for `D` is the
            sum of `(k_i - i)` as `i` ranges from `1` to `r + 1`.  It measure the difference
            between the sequence `r, r - 1, ..., 0, -1, -1, ...` and the rank sequence
            `\mathrm{rank}(D), \mathrm{rank}(D - v), \mathrm{rank}(D - 2v), \dots`
        """
        L = self.weierstrass_rank_seq(v)
        gaps = [i for i in range(1, len(L)) if L[i] != L[i - 1]]
        gaps = tuple(gaps)
        if weight:
            return gaps, sum(gaps) - binomial(len(gaps) + 1, 2)
        else:
            return gaps

    def is_weierstrass_pt(self, v='sink'):
        r"""
        Return whether the given vertex is a Weierstrass point.

        INPUT:

        - ``v`` -- (default: ``sink``) vertex

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: s = sandpiles.House()
            sage: K = s.canonical_divisor()
            sage: K.weierstrass_rank_seq()  # sequence at the sink vertex, 0
            (1, 0, -1)
            sage: K.is_weierstrass_pt()
            False
            sage: K.weierstrass_rank_seq(4)
            (1, 0, 0, -1)
            sage: K.is_weierstrass_pt(4)
            True

        .. NOTE::

            The vertex `v` is a (generalized) Weierstrass point for divisor
            `D` if the sequence of ranks `r(D - nv)`
            for `n = 0, 1, 2, \dots` is not `r(D), r(D)-1, \dots, 0, -1, -1, \dots`
        """
        return self.weierstrass_gap_seq(v)[1] > 0

    def _set_weierstrass_pts(self):
        r"""
        Tuple of Weierstrass vertices.

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: D = SandpileDivisor(s, [2,1,0,0])
            sage: D._set_weierstrass_pts()                                              # needs sage.geometry.polyhedron
            sage: '_weierstrass_pts' in D.__dict__                                      # needs sage.geometry.polyhedron
            True
        """
        self._weierstrass_pts = tuple([v for v in self.sandpile().vertices(sort=True)
                                       if self.is_weierstrass_pt(v)])

    def weierstrass_pts(self, with_rank_seq=False):
        r"""
        The Weierstrass points (vertices). Optionally, return the corresponding
        rank sequences.

        INPUT:

        - ``with_rank_seq`` -- boolean (default: ``False``)

        OUTPUT: tuple of vertices or list of (vertex, rank sequence)

        EXAMPLES::

            sage: s = sandpiles.House()
            sage: K = s.canonical_divisor()
            sage: K.weierstrass_pts()                                                   # needs sage.geometry.polyhedron
            (4,)
            sage: K.weierstrass_pts(True)                                               # needs sage.geometry.polyhedron
            [(4, (1, 0, 0, -1))]

        .. NOTE::

            The vertex `v` is a (generalized) Weierstrass point for divisor
            `D` if the sequence of ranks `r(D - nv)`
            for `n = 0, 1, 2, \dots`` is not `r(D), r(D)-1, \dots, 0, -1, -1, \dots`
        """
        if with_rank_seq:
            rks = [self.weierstrass_rank_seq(v) for v in self._weierstrass_pts]
            return list(zip(self._weierstrass_pts, rks))
        return self._weierstrass_pts

    def weierstrass_div(self, verbose=True):
        r"""
        The Weierstrass divisor.  Its value at a vertex is the weight of that
        vertex as a Weierstrass point.  (See
        ``SandpileDivisor.weierstrass_gap_seq``.)

        INPUT:

        - ``verbose`` -- boolean (default: ``True``)

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: D = SandpileDivisor(s,[4,2,1,0])
            sage: [D.weierstrass_rank_seq(v) for v in s]                                # needs sage.geometry.polyhedron
            [(5, 4, 3, 2, 1, 0, 0, -1),
             (5, 4, 3, 2, 1, 0, -1),
             (5, 4, 3, 2, 1, 0, 0, 0, -1),
             (5, 4, 3, 2, 1, 0, 0, -1)]
            sage: D.weierstrass_div()                                                   # needs sage.geometry.polyhedron
            {0: 1, 1: 0, 2: 2, 3: 1}
            sage: k5 = sandpiles.Complete(5)
            sage: K = k5.canonical_divisor()
            sage: K.weierstrass_div()
            {0: 9, 1: 9, 2: 9, 3: 9, 4: 9}
        """
        s = self.sandpile()
        D = [self.weierstrass_gap_seq(v, True)[1] for v in s.vertices(sort=True)]
        D = SandpileDivisor(s, D)
        if verbose:
            return D
        else:
            return D.values()

    def support(self):
        r"""
        List of vertices at which the divisor is nonzero.

        OUTPUT: list representing the support of the divisor

        EXAMPLES::

            sage: S = sandpiles.Cycle(4)
            sage: D = SandpileDivisor(S, [0,0,1,1])
            sage: D.support()
            [2, 3]
            sage: S.vertices(sort=True)
            [0, 1, 2, 3]
        """
        return [i for i in self if self[i] != 0]

    def _set_Dcomplex(self):
        r"""
        Compute the simplicial complex determined by the supports of the
        linearly equivalent effective divisors.

        EXAMPLES::

            sage: S = sandpiles.Complete(4)
            sage: D = SandpileDivisor(S, [0,0,1,1])
            sage: D._set_Dcomplex()                                                     # needs sage.geometry.polyhedron
            sage: '_Dcomplex' in D.__dict__                                             # needs sage.geometry.polyhedron
            True
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
        while simp:
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
        The support-complex. (See NOTE.)

        OUTPUT: simplicial complex

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: S = sandpiles.House()
            sage: p = SandpileDivisor(S, [1,2,1,0,0]).Dcomplex()
            sage: p.homology()
            {0: 0, 1: Z x Z, 2: 0}
            sage: p.f_vector()
            [1, 5, 10, 4]
            sage: p.betti()
            {0: 1, 1: 2, 2: 0}

        .. NOTE::

            The "support-complex" is the simplicial complex determined by the
            supports of the linearly equivalent effective divisors.
        """
        return self._Dcomplex

    def betti(self):
        r"""
        The Betti numbers for the support-complex.  (See NOTE.)

        OUTPUT: dictionary of integers

        EXAMPLES::

            sage: S = sandpiles.Cycle(3)
            sage: D = SandpileDivisor(S, [2,0,1])
            sage: D.betti()                                                             # needs sage.geometry.polyhedron
            {0: 1, 1: 1}

        .. NOTE::

            The "support-complex" is the simplicial complex determined by the
            supports of the linearly equivalent effective divisors.
        """
        return self.Dcomplex().betti()

    def add_random(self, distrib=None):
        r"""
        Add one grain of sand to a random vertex.

        INPUT:

        - ``distrib`` -- (optional) list of nonnegative numbers representing a probability distribution on the vertices

        OUTPUT: SandpileDivisor

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = s.zero_div()
            sage: D.add_random() # random
            {0: 0, 1: 0, 2: 1, 3: 0}
            sage: D.add_random([0.1,0.1,0.1,0.7]) # random
            {0: 0, 1: 0, 2: 0, 3: 1}

        .. WARNING::

            If ``distrib`` is not ``None``, the user is responsible for assuring the sum of its entries is 1.
        """
        D = deepcopy(self)
        S = self.sandpile()
        V = S.vertices(sort=True)
        if distrib is None:  # default = uniform distribution
            n = S.num_verts()
            distrib = [QQ.one() / n] * n
        X = GeneralDiscreteDistribution(distrib)
        i = X.get_random_element()
        D[V[i]] += 1
        return D

    def is_symmetric(self, orbits):
        r"""
        Return whether the divisor is symmetric.  Return ``True`` if the values of the
        configuration are constant over the vertices in each sublist of
        ``orbits``.

        INPUT:

        - ``orbits`` -- list of lists of vertices

        OUTPUT: boolean

        EXAMPLES::

            sage: S = sandpiles.House()
            sage: S.dict()
            {0: {1: 1, 2: 1},
             1: {0: 1, 3: 1},
             2: {0: 1, 3: 1, 4: 1},
             3: {1: 1, 2: 1, 4: 1},
             4: {2: 1, 3: 1}}
            sage: D = SandpileDivisor(S, [0,0,1,1,3])
            sage: D.is_symmetric([[2,3], [4]])
            True
        """
        return all(len({self[v] for v in x}) <= 1 for x in orbits)

    def _set_life(self):
        r"""
        Will the sequence of divisors `D_i` where `D_{i+1}` is obtained from
        `D_i` by firing all unstable vertices of `D_i` stabilize?  If so,
        save the resulting cycle, otherwise save ``[]``.

        EXAMPLES::

            sage: S = sandpiles.Complete(4)
            sage: D = SandpileDivisor(S, {0: 4, 1: 3, 2: 3, 3: 2})
            sage: D._set_life()
            sage: '_life' in D.__dict__
            True
        """
        oldD = deepcopy(self)
        result = [oldD]
        while True:
            if oldD.unstable() == []:
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
        Return whether the divisor is stabilizable.  In other words, will the divisor stabilize
        under repeated firings of all unstable vertices?  Optionally returns the
        resulting cycle.

        INPUT:

        - ``cycle`` -- boolean (default: ``False``)

        OUTPUT: boolean or optionally, a list of SandpileDivisors

        EXAMPLES::

            sage: S = sandpiles.Complete(4)
            sage: D = SandpileDivisor(S, {0: 4, 1: 3, 2: 3, 3: 2})
            sage: D.is_alive()
            True
            sage: D.is_alive(True)
            [{0: 4, 1: 3, 2: 3, 3: 2}, {0: 3, 1: 2, 2: 2, 3: 5}, {0: 1, 1: 4, 2: 4, 3: 3}]
        """
        return self._life if cycle else bool(self._life)

    def _set_stabilize(self):
        r"""
        The stabilization of the divisor.  If not stabilizable, return an error.

        EXAMPLES::

            sage: s = sandpiles.Diamond()
            sage: D = SandpileDivisor(s, [2,1,0,0])
            sage: D._set_stabilize()
            sage: '_stabilize' in D.__dict__
            True
        """
        if self.is_alive():
            raise RuntimeError('divisor is not stabilizable')
        else:
            firing_vector = self._sandpile.zero_div()
            E = deepcopy(self)
            unstable = E.unstable()
            while unstable != []:
                E = E.fire_unstable()
                for v in unstable:
                    firing_vector[v] += 1
                    unstable = E.unstable()
            self._stabilize = [E, firing_vector]

    def stabilize(self, with_firing_vector=False):
        r"""
        The stabilization of the divisor.  If not stabilizable, return an error.

        INPUT:

        - ``with_firing_vector`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: s = sandpiles.Complete(4)
            sage: D = SandpileDivisor(s,[0,3,0,0])
            sage: D.stabilize()
            {0: 1, 1: 0, 2: 1, 3: 1}
            sage: D.stabilize(with_firing_vector=True)
            [{0: 1, 1: 0, 2: 1, 3: 1}, {0: 0, 1: 1, 2: 0, 3: 0}]
        """
        if with_firing_vector:
            return self._stabilize
        return self._stabilize[0]

    def show(self, heights=True, directed=None, **kwds):
        r"""
        Show the divisor.

        INPUT:

        - ``heights`` -- boolean (default: ``True``); whether to label each
          vertex with the amount of sand

        - ``directed`` -- (optional) whether to draw directed edges

        - ``kwds`` -- (optional) arguments passed to the show method for Graph

        EXAMPLES::

            sage: S = sandpiles.Diamond()
            sage: D = SandpileDivisor(S, [1,-2,0,2])
            sage: D.show(graph_border=True, vertex_size=700, directed=False)            # needs sage.plot
        """
        if directed:
            T = DiGraph(self.sandpile())
        elif directed is False:
            T = Graph(self.sandpile())
        elif self.sandpile().is_directed():
            T = DiGraph(self.sandpile())
        else:
            T = Graph(self.sandpile())

        if heights:
            a = {i: str(i) + ":" + str(T[i]) for i in T}
            T.relabel(a)
        T.show(**kwds)


# See note about this after the definition of SandpileConfig
pretty.for_type(SandpileDivisor, pretty.for_type(dict, None))


#######################################
#         Some test graphs            #
#######################################

def sandlib(selector=None):
    r"""
    Return the sandpile identified by ``selector``.  If no argument is
    given, a description of the sandpiles in the sandlib is printed.

    INPUT:

    - ``selector`` -- (optional) identifier or None

    OUTPUT: Sandpile or description

    EXAMPLES::

            sage: from sage.sandpiles.sandpile import sandlib
            sage: sandlib()
              Sandpiles in the sandlib:
                 ci1 : complete intersection, non-DAG but equivalent to a DAG
                 generic : generic digraph with 6 vertices
                 genus2 : Undirected graph of genus 2
                 gor : Gorenstein but not a complete intersection
                 kite : generic undirected graphs with 5 vertices
                 riemann-roch1 : directed graph with postulation 9 and 3 maximal weight superstables
                 riemann-roch2 : directed graph with a superstable not majorized by a maximal superstable
            sage: S = sandlib('gor')
            sage: S.resolution()                                                        # needs sage.libs.singular
            'R^1 <-- R^5 <-- R^5 <-- R^1'
    """
    # The convention is for the sink to be zero.
    sandpiles = {
        'generic': {'description': 'generic digraph with 6 vertices',
                    'graph': {0: {}, 1: {0: 1, 3: 1, 4: 1},
                              2: {0: 1, 3: 1, 5: 1},
                              3: {2: 1, 5: 1}, 4: {1: 1, 3: 1},
                              5: {2: 1, 3: 1}}},
        'kite': {'description': 'generic undirected graphs with 5 vertices',
                 'graph': {0: {}, 1: {0: 1, 2: 1, 3: 1},
                           2: {1: 1, 3: 1, 4: 1}, 3: {1: 1, 2: 1, 4: 1},
                           4: {2: 1, 3: 1}}},
        'riemann-roch1': {'description': 'directed graph with postulation 9 and 3 maximal weight superstables',
                          'graph': {0: {1: 3, 3: 1},
                                    1: {0: 2, 2: 2, 3: 2},
                                    2: {0: 1, 1: 1},
                                    3: {0: 3, 1: 1, 2: 1}}},
        'riemann-roch2': {'description': 'directed graph with a superstable not majorized by a maximal superstable',
                          'graph': {0: {},
                                    1: {0: 1, 2: 1},
                                    2: {0: 1, 3: 1},
                                    3: {0: 1, 1: 1, 2: 1}}},
        'gor': {'description': 'Gorenstein but not a complete intersection',
                'graph': {0: {},
                          1: {0: 1, 2: 1, 3: 4},
                          2: {3: 5},
                          3: {1: 1, 2: 1}}},
        'ci1': {'description': 'complete intersection, non-DAG but equivalent to a DAG',
                'graph': {0: {}, 1: {2: 2}, 2: {0: 4, 1: 1}}},
        'genus2': {'description': 'Undirected graph of genus 2',
                   'graph': {0: [1, 2],
                             1: [0, 2, 3],
                             2: [0, 1, 3],
                             3: [1, 2]}},
    }
    if selector is None:
        print('')
        print('  Sandpiles in the sandlib:')
        for i in sorted(sandpiles):
            print('    ', i, ':', sandpiles[i]['description'])
        print("")
    elif selector not in sandpiles:
        print(selector, 'is not in the sandlib.')
    else:
        return Sandpile(sandpiles[selector]['graph'], 0)


#################################################
#          Some useful functions                #
#################################################

def triangle_sandpile(n):
    r"""
    A triangular sandpile.  Each nonsink vertex has out-degree six.  The
    vertices on the boundary of the triangle are connected to the sink.

    INPUT:

    - ``n`` -- integer

    OUTPUT: Sandpile

    EXAMPLES::

        sage: from sage.sandpiles.sandpile import triangle_sandpile
        sage: T = triangle_sandpile(5)
        sage: T.group_order()
        135418115000
    """
    T = {(-1, -1):{}}
    for i in range(n):
        for j in range(n-i):
            T[(i,j)] = {}
            if i < n-j-1:
                T[(i,j)][(i+1,j)] = 1
                T[(i,j)][(i,j+1)] = 1
            if i > 0:
                T[(i,j)][(i-1,j+1)] = 1
                T[(i,j)][(i-1,j)] = 1
            if j > 0:
                T[(i,j)][(i,j-1)] = 1
                T[(i,j)][(i+1,j-1)] = 1
            d = len(T[(i,j)])
            if d < 6:
                T[(i,j)][(-1, -1)] = 6-d
    T = Sandpile(T, (-1, -1))
    pos = {}
    for x in T.nonsink_vertices():
        coords = list(x)
        coords[0] += QQ(1)/2*coords[1]
        pos[x] = coords
    pos[(-1, -1)] = (-1, -1)
    T.set_pos(pos)
    return T


def aztec_sandpile(n):
    r"""
    The aztec diamond graph.

    INPUT:

    - ``n`` -- integer

    OUTPUT: dictionary for the aztec diamond graph

    EXAMPLES::

        sage: from sage.sandpiles.sandpile import aztec_sandpile
        sage: T = aztec_sandpile(2)
        sage: sorted(len(v) for u, v in T.items())
        [3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 8]
        sage: Sandpile(T,(0, 0)).group_order()
        4542720

    .. NOTE::

        This is the aztec diamond graph with a sink vertex added.  Boundary
        vertices have edges to the sink so that each vertex has degree 4.
    """
    aztec_sandpile = {}
    half = QQ((1, 2))
    for i in xsrange(n):
        for j in xsrange(n-i):
            aztec_sandpile[(half+i, half+j)] = {}
            aztec_sandpile[(-half-i, half+j)] = {}
            aztec_sandpile[(half+i, -half-j)] = {}
            aztec_sandpile[(-half-i, -half-j)] = {}
    non_sinks = list(aztec_sandpile)
    aztec_sandpile[(0, 0)] = {}
    for vert in non_sinks:
        weight = abs(vert[0]) + abs(vert[1])
        x = vert[0]
        y = vert[1]
        if weight < n:
            aztec_sandpile[vert] = {(x+1, y): 1, (x, y+1): 1, (x-1, y): 1, (x, y-1): 1}
        else:
            if (x+1, y) in aztec_sandpile:
                aztec_sandpile[vert][(x+1, y)] = 1
            if (x, y+1) in aztec_sandpile:
                aztec_sandpile[vert][(x, y+1)] = 1
            if (x-1, y) in aztec_sandpile:
                aztec_sandpile[vert][(x-1, y)] = 1
            if (x, y-1) in aztec_sandpile:
                aztec_sandpile[vert][(x, y-1)] = 1
            if len(aztec_sandpile[vert]) < 4:
                out_degree = 4 - len(aztec_sandpile[vert])
                aztec_sandpile[vert][(0, 0)] = out_degree
                aztec_sandpile[(0, 0)][vert] = out_degree
    return aztec_sandpile


def glue_graphs(g, h, glue_g, glue_h):
    r"""
    Glue two graphs together.

    INPUT:

    - ``g``, ``h`` -- dictionaries for directed multigraphs

    - ``glue_h``, ``glue_g`` -- dictionaries for a vertex

    OUTPUT: dictionary for a directed multigraph

    EXAMPLES::

        sage: from sage.sandpiles.sandpile import glue_graphs
        sage: x = {0: {}, 1: {0: 1}, 2: {0: 1, 1: 1}, 3: {0: 1, 1: 1, 2: 1}}
        sage: y = {0: {}, 1: {0: 2}, 2: {1: 2}, 3: {0: 1, 2: 1}}
        sage: glue_x = {1: 1, 3: 2}
        sage: glue_y = {0: 1, 1: 2, 3: 1}
        sage: z = glue_graphs(x,y,glue_x,glue_y); z
        {'sink': {},
         'x0': {'sink': 1, 'x1': 1, 'x3': 2, 'y1': 2, 'y3': 1},
         'x1': {'x0': 1},
         'x2': {'x0': 1, 'x1': 1},
         'x3': {'x0': 1, 'x1': 1, 'x2': 1},
         'y1': {'sink': 2},
         'y2': {'y1': 2},
         'y3': {'sink': 1, 'y2': 1}}
        sage: S = Sandpile(z,'sink')
        sage: S.h_vector()
        [1, 6, 17, 31, 41, 41, 31, 17, 6, 1]
        sage: S.resolution()                                                            # needs sage.libs.singular
        'R^1 <-- R^7 <-- R^21 <-- R^35 <-- R^35 <-- R^21 <-- R^7 <-- R^1'

    .. NOTE::

        This method makes a dictionary for a graph by combining those for
        `g` and `h`.  The sink of `g` is replaced by a vertex that
        is connected to the vertices of `g` as specified by ``glue_g``
        the vertices of `h` as specified in ``glue_h``.  The sink of the glued
        graph is ``'sink'``.

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
    k = {'sink': {}}  # the new graph dictionary, starting with the sink
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
                    new_edges['sink'] = h[i][j]
                else:
                    new_edges['y'+str(j)] = h[i][j]
            k['y'+str(i)] = new_edges
    # now handle the glue vertex (old g sink)
    new_edges = {}
    for i in glue_g:
        new_edges['x'+str(i)] = glue_g[i]
    for i in glue_h:
        if i == h_sink:
            new_edges['sink'] = glue_h[i]
        else:
            new_edges['y'+str(i)] = glue_h[i]
    k['x'+str(g_sink)] = new_edges
    return k


def firing_graph(S, eff):
    r"""
    Create a digraph with divisors as vertices and edges between two divisors
    `D` and `E` if firing a single vertex in `D` gives `E`.

    INPUT:

    - ``S`` -- Sandpile

    - ``eff`` -- list of divisors

    OUTPUT: DiGraph

    EXAMPLES::

        sage: S = sandpiles.Cycle(6)
        sage: D = SandpileDivisor(S, [1,1,1,1,2,0])
        sage: eff = D.effective_div()                                                   # needs sage.geometry.polyhedron
        sage: firing_graph(S, eff).show3d(edge_size=.005,               # long time, needs sage.geometry.polyhedron sage.plot
        ....:                             vertex_size=0.01)
    """
    g = DiGraph()
    g.add_vertices(range(len(eff)))
    for i in g.vertices(sort=True):
        for v in eff[i]:
            if eff[i][v] >= S.out_degree(v):
                new_div = deepcopy(eff[i])
                new_div[v] -= S.out_degree(v)
                for oe in S.outgoing_edges(v):
                    new_div[oe[1]] += oe[2]
                if new_div in eff:
                    g.add_edge((i, eff.index(new_div)))
    return g


def parallel_firing_graph(S, eff):
    r"""
    Create a digraph with divisors as vertices and edges between two divisors
    `D` and `E` if firing all unstable vertices in `D` gives `E`.

    INPUT:

    - ``S`` -- Sandpile

    - ``eff`` -- list of divisors

    OUTPUT: DiGraph

    EXAMPLES::

        sage: S = sandpiles.Cycle(6)
        sage: D = SandpileDivisor(S, [1,1,1,1,2,0])
        sage: eff = D.effective_div()                                                   # needs sage.geometry.polyhedron
        sage: parallel_firing_graph(S, eff).show3d(edge_size=.005,      # long time, needs sage.geometry.polyhedron sage.plot
        ....:                                      vertex_size=0.01)
    """
    g = DiGraph()
    g.add_vertices(range(len(eff)))
    for i in g.vertices(sort=True):
        new_edge = False
        new_div = deepcopy(eff[i])
        for v in eff[i]:
            if eff[i][v] >= S.out_degree(v):
                new_edge = True
                new_div[v] -= S.out_degree(v)
                for oe in S.outgoing_edges(v):
                    new_div[oe[1]] += oe[2]
        if new_edge and (new_div in eff):
            g.add_edge((i, eff.index(new_div)))
    return g


def admissible_partitions(S, k):
    r"""
    The partitions of the vertices of `S` into `k` parts, each of which is
    connected.

    INPUT:

    - ``S`` -- Sandpile

    - ``k`` -- integer

    OUTPUT: partitions

    EXAMPLES::

        sage: from sage.sandpiles.sandpile import admissible_partitions
        sage: from sage.sandpiles.sandpile import partition_sandpile
        sage: S = sandpiles.Cycle(4)
        sage: P = [list(admissible_partitions(S, i)) for i in [2,3,4]]                  # needs sage.combinat
        sage: P                                                                         # needs sage.combinat
        [[{{0, 2, 3}, {1}},
          {{0, 3}, {1, 2}},
          {{0, 1, 3}, {2}},
          {{0}, {1, 2, 3}},
          {{0, 1}, {2, 3}},
          {{0, 1, 2}, {3}}],
         [{{0, 3}, {1}, {2}},
          {{0}, {1}, {2, 3}},
          {{0}, {1, 2}, {3}},
          {{0, 1}, {2}, {3}}],
         [{{0}, {1}, {2}, {3}}]]
        sage: for p in P:                                                               # needs sage.combinat
        ....:  sum([partition_sandpile(S, i).betti(verbose=False)[-1] for i in p])
        6
        8
        3
        sage: S.betti()                                                                 # needs sage.libs.singular
                   0     1     2     3
        ------------------------------
            0:     1     -     -     -
            1:     -     6     8     3
        ------------------------------
        total:     1     6     8     3
    """
    v = S.vertices(sort=True)
    G = DiGraph(S) if S.is_directed() else Graph(S)
    for p in SetPartitions(v, k):
        if all(G.subgraph(list(x)).is_connected() for x in p):
            yield p


def partition_sandpile(S, p):
    r"""
    Each set of vertices in `p` is regarded as a single vertex, with and edge
    between `A` and `B` if some element of `A` is connected by an edge to  some
    element of `B` in `S`.

    INPUT:

    - ``S`` -- Sandpile

    - ``p`` -- partition of the vertices of ``S``

    OUTPUT: Sandpile

    EXAMPLES::

        sage: from sage.sandpiles.sandpile import admissible_partitions, partition_sandpile
        sage: S = sandpiles.Cycle(4)
        sage: P = [list(admissible_partitions(S, i)) for i in [2,3,4]]                  # needs sage.combinat
        sage: for p in P:                                                               # needs sage.combinat
        ....:  sum([partition_sandpile(S, i).betti(verbose=False)[-1] for i in p])
        6
        8
        3
        sage: S.betti()                                                                 # needs sage.libs.singular
                   0     1     2     3
        ------------------------------
            0:     1     -     -     -
            1:     -     6     8     3
        ------------------------------
        total:     1     6     8     3
    """
    from itertools import combinations
    g = Graph()
    g.add_vertices(tuple(i) for i in p)
    for u, v in combinations(g, 2):
        for i in u:
            for j in v:
                if (i, j, 1) in S.edges(sort=False):
                    g.add_edge((u, v))
                    break
    for i in g.vertices(sort=False):
        if S.sink() in i:
            return Sandpile(g, i)


def min_cycles(G, v):
    r"""
    Minimal length cycles in the digraph `G` starting at vertex `v`.

    INPUT:

    - ``G`` -- DiGraph

    - ``v`` -- vertex of ``G``

    OUTPUT: list of lists of vertices

    EXAMPLES::

        sage: from sage.sandpiles.sandpile import min_cycles, sandlib
        sage: T = sandlib('gor')
        sage: [min_cycles(T, i) for i in T.vertices(sort=True)]
        [[], [[1, 3]], [[2, 3, 1], [2, 3]], [[3, 1], [3, 2]]]
    """
    pr = G.neighbors_in(v)
    sp = G.shortest_paths(v)
    return [sp[i] for i in pr if i in sp]


def wilmes_algorithm(M):
    r"""
    Compute an integer matrix `L` with the same integer row span as `M` and
    such that `L` is the reduced Laplacian of a directed multigraph.

    INPUT:

    - ``M`` -- square integer matrix of full rank

    OUTPUT: integer matrix (``L``)

    EXAMPLES::

        sage: P = matrix([[2,3,-7,-3],[5,2,-5,5],[8,2,5,4],[-5,-9,6,6]])
        sage: wilmes_algorithm(P)
        [ 3279   -79 -1599 -1600]
        [   -1  1539  -136 -1402]
        [    0    -1  1650 -1649]
        [    0     0 -1658  1658]

    REFERENCES:

    - [PPW2013]_
    """
    # find the gcd of the row-sums, and perform the corresponding row
    # operations on M
    if M.matrix_over_field().is_invertible():
        L = deepcopy(M)
        L = matrix(ZZ, L)
        U = matrix(ZZ, [sum(i) for i in L]).smith_form()[2].transpose()
        L = U*M
        for k in range(1, M.nrows()-1):
            smith = matrix(ZZ, [i[k-1] for i in L[k:]]).smith_form()[2].transpose()
            U = identity_matrix(ZZ, k).block_sum(smith)
            L = U*L
            L[k] = -L[k]
        if L[-1][-2] > 0:
            L[-1] = -L[-1]
        for k in range(M.nrows()-2, -1, -1):
            for i in range(k+2, M.nrows()):
                while L[k][i-1] > 0:
                    L[k] = L[k] + L[i]
            v = -L[k+1]
            for i in range(k+2, M.nrows()):
                v = abs(L[i, i-1])*v + v[i-1]*L[i]
            while L[k, k] <= 0 or L[k, -1] > 0:
                L[k] = L[k] + v
        return L

    raise ValueError('matrix not of full rank')
