r"""
Hypergraph generators

This module implements generators of hypergraphs. All hypergraphs can be built
through the ``hypergraphs`` object. For instance, to build a complete 3-uniform
hypergraph on 5 points, one can do::

    sage: H = hypergraphs.CompleteUniform(5, 3)

To enumerate hypergraphs with certain properties up to isomorphism, one can use
method :meth:`~nauty`, which calls Brendan McKay's Nauty
(`<http://cs.anu.edu.au/~bdm/nauty/>`_)::

    sage: list(hypergraphs.nauty(2, 2, connected=True))
    [((0,), (0, 1))]


**This module contains the following hypergraph generators**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~HypergraphGenerators.nauty` | Enumerate hypergraphs up to isomorphism using Nauty.
    :meth:`~HypergraphGenerators.CompleteUniform` | Return the complete `k`-uniform hypergraph on `n` points.
    :meth:`~HypergraphGenerators.UniformRandomUniform` | Return a uniformly sampled `k`-uniform hypergraph on `n` points with `m` hyperedges.


Functions and methods
---------------------
"""
from sage.env import SAGE_NAUTY_BINS_PREFIX as nautyprefix

class HypergraphGenerators():
    r"""
    A class consisting of constructors for common hypergraphs.
    """

    def nauty(self, number_of_sets, number_of_vertices,
              multiple_sets=False,
              vertex_min_degree=None, vertex_max_degree=None,
              set_max_size=None, set_min_size=None,
              regular=False, uniform=False,
              max_intersection=None,
              connected=False,
              debug=False, options=""):
        r"""
        Enumerate hypergraphs up to isomorphism using Nauty.

        INPUT:

        - ``number_of_sets`` -- integer, at most 64 minus ``number_of_vertices``

        - ``number_of_vertices`` -- integer, at most 30

        - ``multiple_sets`` -- boolean (default: ``False``); whether to allow
          several sets of the hypergraph to be equal.

        - ``vertex_min_degree``, ``vertex_max_degree`` -- integers (default:
          ``None``); define the maximum and minimum degree of an element from
          the ground set (i.e. the number of sets which contain it).

        - ``set_min_size``, ``set_max_size`` -- integers (default: ``None``);
          define the maximum and minimum size of a set.

        - ``regular`` -- integers (default: ``False``); if set to an integer
          value `k`, requires the hypergraphs to be `k`-regular. It is actually
          a shortcut for the corresponding min/max values.

        - ``uniform`` -- integers (default: ``False``); if set to an integer
          value `k`, requires the hypergraphs to be `k`-uniform. It is actually
          a shortcut for the corresponding min/max values.

        - ``max_intersection`` -- integers (default: ``None``); constraints the
          maximum cardinality of the intersection of two sets from the
          hypergraphs.

        - ``connected`` -- boolean (default: ``False``); whether to require the
          hypergraphs to be connected.

        - ``debug`` -- boolean (default: ``False``); if ``True`` the first line
          of genbgL's output to standard error is captured and the first call to
          the generator's ``next()`` function will return this line as a string.
          A line leading with ">A" indicates a successful initiation of the
          program with some information on the arguments, while a line beginning
          with ">E" indicates an error with the input.

        - ``options`` -- string (default: ``""``) -- anything else that should
          be forwarded as input to Nauty's genbgL. See its documentation for more
          information : `<http://cs.anu.edu.au/~bdm/nauty/>`_.

          .. NOTE::

              For genbgL the *first class* elements are vertices, and *second
              class* elements are the hypergraph's sets.

        OUTPUT:

        A tuple of tuples.

        EXAMPLES:

        Small hypergraphs::

            sage: list(hypergraphs.nauty(4, 2))
            [((), (0,), (1,), (0, 1))]

        Only connected ones::

            sage: list(hypergraphs.nauty(2, 2, connected=True))
            [((0,), (0, 1))]

        Non-empty sets only::

            sage: list(hypergraphs.nauty(3, 2, set_min_size=1))
            [((0,), (1,), (0, 1))]

        The Fano Plane, as the only 3-uniform hypergraph with 7 sets and 7
        vertices::

            sage: fano = next(hypergraphs.nauty(7, 7, uniform=3, max_intersection=1))
            sage: print(fano)
            ((0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (2, 4, 5), (2, 3, 6), (1, 4, 6))

        The Fano Plane, as the only 3-regular hypergraph with 7 sets and 7
        vertices::

            sage: fano = next(hypergraphs.nauty(7, 7, regular=3, max_intersection=1))
            sage: print(fano)
            ((0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (2, 4, 5), (2, 3, 6), (1, 4, 6))

        TESTS::

            sage: len(list(hypergraphs.nauty(20, 20, uniform=2, regular=2,max_intersection=1)))
            49
            sage: list(hypergraphs.nauty(40, 40, uniform=2, regular=2,max_intersection=1))
            Traceback (most recent call last):
            ...
            ValueError: cannot have more than 30 vertices
            sage: list(hypergraphs.nauty(40, 30, uniform=2, regular=2,max_intersection=1))
            Traceback (most recent call last):
            ...
            ValueError: cannot have more than 64 sets+vertices
        """
        if number_of_vertices > 30:
            raise ValueError("cannot have more than 30 vertices")
        if number_of_sets + number_of_vertices > 64:
            raise ValueError("cannot have more than 64 sets+vertices")

        import subprocess

        nauty_input = options

        if connected:
            nauty_input += " -c"

        if not multiple_sets:
            nauty_input += " -z"

        if max_intersection is not None:
            nauty_input += " -Z" + str(max_intersection)

        # degrees and sizes
        if regular is not False:
            vertex_max_degree = vertex_min_degree = regular
        if vertex_max_degree is None:
            vertex_max_degree = number_of_sets
        if vertex_min_degree is None:
            vertex_min_degree = 0

        if uniform is not False:
            set_max_size = set_min_size = uniform
        if set_max_size is None:
            set_max_size = number_of_vertices
        if set_min_size is None:
            set_min_size = 0

        nauty_input += " -d" + str(vertex_min_degree) + ":" + str(set_min_size)
        nauty_input += " -D" + str(vertex_max_degree) + ":" + str(set_max_size)

        nauty_input +=  " " + str(number_of_vertices) + " " + str(number_of_sets) + " "

        sp = subprocess.Popen(nautyprefix + "genbgL {0}".format(nauty_input), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        if debug:
            yield sp.stderr.readline()

        gen = sp.stdout
        total = number_of_sets + number_of_vertices
        from sage.graphs.graph import Graph
        while True:
            try:
                s = next(gen)
            except StopIteration:
                # Exhausted list of graphs from nauty geng
                return

            G = Graph(s[:-1], format='graph6')

            yield tuple(tuple(G.neighbor_iterator(v)) for v in range(number_of_vertices, total))

    def CompleteUniform(self, n, k):
        r"""
        Return the complete `k`-uniform hypergraph on `n` points.

        INPUT:

        - ``k,n`` -- nonnegative integers with `k\leq n`

        EXAMPLES::

            sage: h = hypergraphs.CompleteUniform(5, 2); h
            Incidence structure with 5 points and 10 blocks
            sage: len(h.packing())
            2
        """
        from sage.combinat.designs.incidence_structures import IncidenceStructure
        from itertools import combinations
        return IncidenceStructure(points=n, blocks=list(combinations(range(n), k)))

    def UniformRandomUniform(self, n, k, m):
        r"""
        Return a uniformly sampled `k`-uniform hypergraph on `n` points with
        `m` hyperedges.

        - ``n`` -- number of nodes of the graph

        - ``k`` -- uniformity

        - ``m`` -- number of edges

        EXAMPLES::

            sage: H = hypergraphs.UniformRandomUniform(52, 3, 17)
            sage: H
            Incidence structure with 52 points and 17 blocks
            sage: H.is_connected()
            False

        TESTS::

            sage: hypergraphs.UniformRandomUniform(-52, 3, 17)
            Traceback (most recent call last):
            ...
            ValueError: number of vertices should be non-negative
            sage: hypergraphs.UniformRandomUniform(52.9, 3, 17)
            Traceback (most recent call last):
            ...
            ValueError: number of vertices should be an integer
            sage: hypergraphs.UniformRandomUniform(52, -3, 17)
            Traceback (most recent call last):
            ...
            ValueError: the uniformity should be non-negative
            sage: hypergraphs.UniformRandomUniform(52, I, 17)
            Traceback (most recent call last):
            ...
            ValueError: the uniformity should be an integer
        """
        from sage.rings.integer import Integer
        from sage.combinat.subset import Subsets
        from sage.misc.prandom import sample

        # Construct the vertex set
        if n < 0:
            raise ValueError("number of vertices should be non-negative")
        try:
            nverts = Integer(n)
        except TypeError:
            raise ValueError("number of vertices should be an integer")
        vertices = list(range(nverts))

        # Construct the edge set
        if k < 0:
            raise ValueError("the uniformity should be non-negative")
        try:
            uniformity = Integer(k)
        except TypeError:
            raise ValueError("the uniformity should be an integer")
        all_edges = Subsets(vertices, uniformity)
        try:
            edges = [all_edges[t] for t in sample(range(len(all_edges)), m)]
        except OverflowError:
            raise OverflowError("binomial({}, {}) too large to be treated".format(n, k))
        except ValueError:
            raise ValueError("number of edges m must be between 0 and binomial({}, {})".format(n, k))

        from sage.combinat.designs.incidence_structures import IncidenceStructure
        return IncidenceStructure(points=vertices, blocks=edges)

    def BinomialRandomUniform(self, n, k, p):
        r"""
        Return a random `k`-uniform hypergraph on `n` points, in which each
        edge is inserted independently with probability `p`.

        - ``n`` -- number of nodes of the graph

        - ``k`` -- uniformity

        - ``p`` -- probability of an edge

        EXAMPLES::

            sage: hypergraphs.BinomialRandomUniform(50, 3, 1).num_blocks()
            19600
            sage: hypergraphs.BinomialRandomUniform(50, 3, 0).num_blocks()
            0

        TESTS::

            sage: hypergraphs.BinomialRandomUniform(50, 3, -0.1)
            Traceback (most recent call last):
            ...
            ValueError: edge probability should be in [0,1]
            sage: hypergraphs.BinomialRandomUniform(50, 3, 1.1)
            Traceback (most recent call last):
            ...
            ValueError: edge probability should be in [0,1]
            sage: hypergraphs.BinomialRandomUniform(-50, 3, 0.17)
            Traceback (most recent call last):
            ...
            ValueError: number of vertices should be non-negative
            sage: hypergraphs.BinomialRandomUniform(50.9, 3, 0.17)
            Traceback (most recent call last):
            ...
            ValueError: number of vertices should be an integer
            sage: hypergraphs.BinomialRandomUniform(50, -3, 0.17)
            Traceback (most recent call last):
            ...
            ValueError: the uniformity should be non-negative
            sage: hypergraphs.BinomialRandomUniform(50, I, 0.17)
            Traceback (most recent call last):
            ...
            ValueError: the uniformity should be an integer
        """
        from sage.rings.integer import Integer
        if n < 0:
            raise ValueError("number of vertices should be non-negative")
        try:
            nverts = Integer(n)
        except TypeError:
            raise ValueError("number of vertices should be an integer")
        if k < 0:
            raise ValueError("the uniformity should be non-negative")
        try:
            uniformity = Integer(k)
        except TypeError:
            raise ValueError("the uniformity should be an integer")
        if not 0 <= p <= 1:
            raise ValueError("edge probability should be in [0,1]")

        import numpy.random as nrn
        from sage.functions.other import binomial
        m = nrn.binomial(binomial(nverts, uniformity), p)
        return hypergraphs.UniformRandomUniform(n, k, m)


hypergraphs = HypergraphGenerators()
