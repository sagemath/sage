r"""
Common Digraphs

All digraphs in Sage can be built through the ``digraphs`` object. In order to
build a circuit on 15 elements, one can do::

    sage: g = digraphs.Circuit(15)

To get a circulant graph on 10 vertices in which a vertex `i` has `i+2` and
`i+3` as outneighbors::

    sage: p = digraphs.Circulant(10,[2,3])

More interestingly, one can get the list of all digraphs that Sage knows how to
build by typing ``digraphs.`` in Sage and then hitting tab.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraphGenerators.ButterflyGraph`      | Returns a n-dimensional butterfly graph.
    :meth:`~DiGraphGenerators.Circuit`             | Returns the circuit on `n` vertices.
    :meth:`~DiGraphGenerators.Circulant`           | Returns a circulant digraph on `n` vertices from a set of integers.
    :meth:`~DiGraphGenerators.Complete`            | Return a complete digraph on `n` vertices.
    :meth:`~DiGraphGenerators.DeBruijn`            | Returns the De Bruijn digraph with parameters `k,n`.
    :meth:`~DiGraphGenerators.GeneralizedDeBruijn` | Returns the generalized de Bruijn digraph of order `n` and degree `d`.
    :meth:`~DiGraphGenerators.ImaseItoh`           | Returns the digraph of Imase and Itoh of order `n` and degree `d`.
    :meth:`~DiGraphGenerators.Kautz`               | Returns the Kautz digraph of degree `d` and diameter `D`.
    :meth:`~DiGraphGenerators.Path`                | Returns a directed path on `n` vertices.
    :meth:`~DiGraphGenerators.RandomDirectedGNC`   | Returns a random GNC (growing network with copying) digraph with `n` vertices.
    :meth:`~DiGraphGenerators.RandomDirectedGNM`   | Returns a random labelled digraph on `n` nodes and `m` arcs.
    :meth:`~DiGraphGenerators.RandomDirectedGNP`   | Returns a random digraph on `n` nodes.
    :meth:`~DiGraphGenerators.RandomDirectedGN`    | Returns a random GN (growing network) digraph with `n` vertices.
    :meth:`~DiGraphGenerators.RandomDirectedGNR`   | Returns a random GNR (growing network with redirection) digraph.
    :meth:`~DiGraphGenerators.RandomSemiComplete`  | Return a random semi-complete digraph of order `n`.
    :meth:`~DiGraphGenerators.RandomTournament`    | Returns a random tournament on `n` vertices.
    :meth:`~DiGraphGenerators.TransitiveTournament`| Returns a transitive tournament on `n` vertices.
    :meth:`~DiGraphGenerators.tournaments_nauty`   | Returns all tournaments on `n` vertices using Nauty.


AUTHORS:

- Robert L. Miller (2006)
- Emily A. Kirkman (2006)
- Michael C. Yurko (2009)
- David Coudert    (2012)

Functions and methods
---------------------

"""

################################################################################
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################
from __future__ import print_function
from six.moves import range

from   math import sin, cos, pi
from sage.misc.randstate import current_randstate
from sage.graphs.digraph import DiGraph


class DiGraphGenerators():
    r"""
    A class consisting of constructors for several common digraphs,
    including orderly generation of isomorphism class representatives.

    A list of all graphs and graph structures in this database is
    available via tab completion. Type "digraphs." and then hit tab to
    see which graphs are available.

    The docstrings include educational information about each named
    digraph with the hopes that this class can be used as a reference.

    The constructors currently in this class include::

                Random Directed Graphs:
                    - RandomDirectedGN
                    - RandomDirectedGNC
                    - RandomDirectedGNP
                    - RandomDirectedGNM
                    - RandomDirectedGNR
                    - RandomTournament
                    - RandomSemiComplete

                Families of Graphs:
                    - Complete
                    - DeBruijn
                    - GeneralizedDeBruijn
                    - Kautz
                    - Path
                    - ImaseItoh
                    - RandomTournament
                    - TransitiveTournament
                    - tournaments_nauty



    ORDERLY GENERATION: digraphs(vertices, property=lambda x: True,
    augment='edges', size=None)

    Accesses the generator of isomorphism class representatives.
    Iterates over distinct, exhaustive representatives.

    INPUT:


    -  ``vertices`` - natural number

    -  ``property`` - any property to be tested on digraphs
       before generation.

    -  ``augment`` - choices:

    -  ``'vertices'`` - augments by adding a vertex, and
       edges incident to that vertex. In this case, all digraphs on up to
       n=vertices are generated. If for any digraph G satisfying the
       property, every subgraph, obtained from G by deleting one vertex
       and only edges incident to that vertex, satisfies the property,
       then this will generate all digraphs with that property. If this
       does not hold, then all the digraphs generated will satisfy the
       property, but there will be some missing.

    -  ``'edges'`` - augments a fixed number of vertices by
       adding one edge In this case, all digraphs on exactly n=vertices
       are generated. If for any graph G satisfying the property, every
       subgraph, obtained from G by deleting one edge but not the vertices
       incident to that edge, satisfies the property, then this will
       generate all digraphs with that property. If this does not hold,
       then all the digraphs generated will satisfy the property, but
       there will be some missing.

    -  ``implementation`` - which underlying implementation to use (see DiGraph?)

    -  ``sparse`` - ignored if implementation is not ``c_graph``

    EXAMPLES: Print digraphs on 2 or less vertices.

    ::

        sage: for D in digraphs(2, augment='vertices'):
        ....:     print(D)
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices

    Note that we can also get digraphs with underlying Cython implementation::

        sage: for D in digraphs(2, augment='vertices', implementation='c_graph'):
        ....:     print(D)
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices

    Print digraphs on 3 vertices.

    ::

        sage: for D in digraphs(3):
        ....:     print(D)
        Digraph on 3 vertices
        Digraph on 3 vertices
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices

    Generate all digraphs with 4 vertices and 3 edges.

    ::

        sage: L = digraphs(4, size=3)
        sage: len(list(L))
        13

    Generate all digraphs with 4 vertices and up to 3 edges.

    ::

        sage: L = list(digraphs(4, lambda G: G.size() <= 3))
        sage: len(L)
        20
        sage: graphs_list.show_graphs(L)  # long time

    Generate all digraphs with degree at most 2, up to 5 vertices.

    ::

        sage: property = lambda G: ( max([G.degree(v) for v in G] + [0]) <= 2 )
        sage: L = list(digraphs(5, property, augment='vertices'))
        sage: len(L)
        75

    Generate digraphs on the fly: (see http://oeis.org/classic/A000273)

    ::

        sage: for i in range(5):
        ....:     print(len(list(digraphs(i))))
        1
        1
        3
        16
        218

    REFERENCE:

    - Brendan D. McKay, Isomorph-Free Exhaustive generation.  Journal
      of Algorithms Volume 26, Issue 2, February 1998, pages 306-324.
    """

    def ButterflyGraph(self, n, vertices='strings'):
        """
        Returns a n-dimensional butterfly graph. The vertices consist of
        pairs (v,i), where v is an n-dimensional tuple (vector) with binary
        entries (or a string representation of such) and i is an integer in
        [0..n]. A directed edge goes from (v,i) to (w,i+1) if v and w are
        identical except for possibly v[i] != w[i].

        A butterfly graph has `(2^n)(n+1)` vertices and
        `n2^{n+1}` edges.

        INPUT:


        -  ``vertices`` - 'strings' (default) or 'vectors',
           specifying whether the vertices are zero-one strings or actually
           tuples over GF(2).


        EXAMPLES::

            sage: digraphs.ButterflyGraph(2).edges(labels=False)
            [(('00', 0), ('00', 1)),
            (('00', 0), ('10', 1)),
            (('00', 1), ('00', 2)),
            (('00', 1), ('01', 2)),
            (('01', 0), ('01', 1)),
            (('01', 0), ('11', 1)),
            (('01', 1), ('00', 2)),
            (('01', 1), ('01', 2)),
            (('10', 0), ('00', 1)),
            (('10', 0), ('10', 1)),
            (('10', 1), ('10', 2)),
            (('10', 1), ('11', 2)),
            (('11', 0), ('01', 1)),
            (('11', 0), ('11', 1)),
            (('11', 1), ('10', 2)),
            (('11', 1), ('11', 2))]
            sage: digraphs.ButterflyGraph(2,vertices='vectors').edges(labels=False)
            [(((0, 0), 0), ((0, 0), 1)),
            (((0, 0), 0), ((1, 0), 1)),
            (((0, 0), 1), ((0, 0), 2)),
            (((0, 0), 1), ((0, 1), 2)),
            (((0, 1), 0), ((0, 1), 1)),
            (((0, 1), 0), ((1, 1), 1)),
            (((0, 1), 1), ((0, 0), 2)),
            (((0, 1), 1), ((0, 1), 2)),
            (((1, 0), 0), ((0, 0), 1)),
            (((1, 0), 0), ((1, 0), 1)),
            (((1, 0), 1), ((1, 0), 2)),
            (((1, 0), 1), ((1, 1), 2)),
            (((1, 1), 0), ((0, 1), 1)),
            (((1, 1), 0), ((1, 1), 1)),
            (((1, 1), 1), ((1, 0), 2)),
            (((1, 1), 1), ((1, 1), 2))]
        """
        # We could switch to Sage integers to handle arbitrary n.
        if vertices=='strings':
            if n>=31:
                raise NotImplementedError("vertices='strings' is only valid for n<=30.")
            from sage.graphs.generic_graph_pyx import int_to_binary_string
            butterfly = {}
            for v in range(2 ** n):
                for i in range(n):
                    w = v
                    w ^= (1 << i)   # push 1 to the left by i and xor with w
                    bv = int_to_binary_string(v)
                    bw = int_to_binary_string(w)
                    # pad and reverse the strings
                    padded_bv = ('0'*(n-len(bv))+bv)[::-1]
                    padded_bw = ('0'*(n-len(bw))+bw)[::-1]
                    butterfly[(padded_bv,i)]=[(padded_bv,i+1), (padded_bw,i+1)]
        elif vertices=='vectors':
            from sage.modules.free_module import VectorSpace
            from sage.rings.finite_rings.finite_field_constructor import FiniteField
            from copy import copy
            butterfly = {}
            for v in VectorSpace(FiniteField(2),n):
                for i in range(n):
                    w=copy(v)
                    w[i] += 1 # Flip the ith bit
                    # We must call tuple since vectors are mutable.  To obtain
                    # a vector from the tuple t, just call vector(t).
                    butterfly[(tuple(v),i)]=[(tuple(v),i+1), (tuple(w),i+1)]
        else:
            raise NotImplementedError("vertices must be 'strings' or 'vectors'.")
        return DiGraph(butterfly)

    def Path(self,n):
        r"""
        Returns a directed path on `n` vertices.

        INPUT:

        - ``n`` (integer) -- number of vertices in the path.

        EXAMPLES::

            sage: g = digraphs.Path(5)
            sage: g.vertices()
            [0, 1, 2, 3, 4]
            sage: g.size()
            4
            sage: g.automorphism_group().cardinality()
            1
        """
        g = DiGraph(n)
        g.name("Path")

        if n:
            g.add_path(list(range(n)))

        g.set_pos({i:(i,0) for i in range(n)})
        return g

    def TransitiveTournament(self, n):
        r"""
        Returns a transitive tournament on `n` vertices.

        In this tournament there is an edge from `i` to `j` if `i<j`.

        See :wikipedia:`Tournament_(graph_theory)`

        INPUT:

        - ``n`` (integer) -- number of vertices in the tournament.

        EXAMPLES::

            sage: g = digraphs.TransitiveTournament(5)
            sage: g.vertices()
            [0, 1, 2, 3, 4]
            sage: g.size()
            10
            sage: g.automorphism_group().cardinality()
            1

        TESTS::

            sage: digraphs.TransitiveTournament(-1)
            Traceback (most recent call last):
            ...
            ValueError: The number of vertices cannot be strictly negative!
        """
        g = DiGraph(n)
        g.name("Transitive Tournament")

        for i in range(n-1):
            for j in range(i+1, n):
                g.add_edge(i, j)

        if n:
            from sage.graphs.graph_plot import _circle_embedding
            _circle_embedding(g, list(range(n)))

        return g

    def RandomTournament(self, n):
        r"""
        Returns a random tournament on `n` vertices.

        For every pair of vertices, the tournament has an edge from
        `i` to `j` with probability `1/2`, otherwise it has an edge
        from `j` to `i`.

        See :wikipedia:`Tournament_(graph_theory)`

        INPUT:

        - ``n`` (integer) -- number of vertices.

        .. SEEALSO::

            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.Complete`

            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.RandomSemiComplete`

        EXAMPLES::

            sage: T = digraphs.RandomTournament(10); T
            Random Tournament: Digraph on 10 vertices
            sage: T.size() == binomial(10, 2)
            True
            sage: digraphs.RandomTournament(-1)
            Traceback (most recent call last):
            ...
            ValueError: The number of vertices cannot be strictly negative!
        """
        from sage.misc.prandom import random
        g = DiGraph(n)
        g.name("Random Tournament")

        for i in range(n-1):
            for j in range(i+1, n):
                if random() <= .5:
                    g.add_edge(i, j)
                else:
                    g.add_edge(j, i)

        if n:
            from sage.graphs.graph_plot import _circle_embedding
            _circle_embedding(g, list(range(n)))

        return g

    def tournaments_nauty(self, n,
                          min_out_degree = None, max_out_degree = None,
                          strongly_connected = False, debug=False, options=""):
        r"""
        Returns all tournaments on `n` vertices using Nauty.

        INPUT:

        - ``n`` (integer) -- number of vertices.

        - ``min_out_degree``, ``max_out_degree`` (integers) -- if set to
          ``None`` (default), then the min/max out-degree is not constrained.

        - ``debug`` (boolean) -- if ``True`` the first line of genbg's output to
          standard error is captured and the first call to the generator's
          ``next()`` function will return this line as a string.  A line leading
          with ">A" indicates a successful initiation of the program with some
          information on the arguments, while a line beginning with ">E"
          indicates an error with the input.

        - ``options`` (string) -- anything else that should be forwarded as
          input to Nauty's genbg. See its documentation for more information :
          `<http://cs.anu.edu.au/~bdm/nauty/>`_.


        .. NOTE::

            To use this method you must first install the Nauty spkg.

        EXAMPLES::

            sage: for g in digraphs.tournaments_nauty(4):
            ....:     print(g.edges(labels = False))
            [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
            [(1, 0), (1, 3), (2, 0), (2, 1), (3, 0), (3, 2)]
            [(0, 2), (1, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
            [(0, 2), (0, 3), (1, 0), (2, 1), (3, 1), (3, 2)]
            sage: tournaments = digraphs.tournaments_nauty
            sage: [len(list(tournaments(x))) for x in range(1,8)]
            [1, 1, 2, 4, 12, 56, 456]
            sage: [len(list(tournaments(x, strongly_connected = True))) for x in range(1,9)]
            [1, 0, 1, 1, 6, 35, 353, 6008]
        """
        import subprocess

        nauty_input = options

        if min_out_degree is None:
            min_out_degree = 0
        if max_out_degree is None:
            max_out_degree = n-1

        nauty_input += " -d"+str(min_out_degree)
        nauty_input += " -D"+str(max_out_degree)

        if strongly_connected:
            nauty_input += " -c"

        nauty_input +=  " "+str(n) +" "

        sp = subprocess.Popen("gentourng {0}".format(nauty_input), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        if debug:
            yield sp.stderr.readline()

        gen = sp.stdout
        while True:
            try:
                s = next(gen)
            except StopIteration:
                raise StopIteration("Exhausted list of graphs from nauty geng")

            G = DiGraph(n)
            i = 0
            j = 1
            for b in s[:-1]:
                if b == '0':
                    G.add_edge(i,j)
                else:
                    G.add_edge(j,i)

                if j == n-1:
                    i += 1
                    j = i+1
                else:
                    j += 1

            yield G


    def Complete(self, n, loops=False):
        r"""
        Return the complete digraph on `n` vertices.

        INPUT:

        - ``n`` (integer) -- number of vertices.

        - ``loops`` (boolean) -- whether to add loops or not, i.e., edges from
          `u` to itself.

        .. SEEALSO::

            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.RandomSemiComplete`

            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.RandomTournament`

        EXAMPLES::

            sage: n = 10
            sage: G = digraphs.Complete(n); G
            Complete digraph: Digraph on 10 vertices
            sage: G.size() == n*(n-1)
            True
            sage: G = digraphs.Complete(n, loops=True); G
            Complete digraph with loops: Looped digraph on 10 vertices
            sage: G.size() == n*n
            True
            sage: digraphs.Complete(-1)
            Traceback (most recent call last):
            ...
            ValueError: The number of vertices cannot be strictly negative!
        """
        G = DiGraph(n, name="Complete digraph"+(" with loops" if loops else ''), loops=loops)

        if loops:
            G.add_edges((u,u) for u in range(n))

        G.add_edges((u,v) for u in range(n) for v in range(n) if u!=v)

        if n:
            from sage.graphs.graph_plot import _circle_embedding
            _circle_embedding(G, list(range(n)))

        return G

    def Circuit(self,n):
        r"""
        Returns the circuit on `n` vertices

        The circuit is an oriented ``CycleGraph``

        EXAMPLE:

        A circuit is the smallest strongly connected digraph::

            sage: circuit = digraphs.Circuit(15)
            sage: len(circuit.strongly_connected_components()) == 1
            True
        """
        g = DiGraph(n)
        g.name("Circuit")

        if n==0:
            return g
        elif n == 1:
            g.allow_loops(True)
            g.add_edge(0,0)
            return g
        else:
            g.add_edges([(i,i+1) for i in range(n-1)])
            g.add_edge(n-1,0)
            return g

    def Circulant(self,n,integers):
        r"""
        Returns a circulant digraph on `n` vertices from a set of integers.

        INPUT:

        - ``n`` (integer) -- number of vertices.

        - ``integers`` -- the list of integers such that there is an edge from
          `i` to `j` if and only if ``(j-i)%n in integers``.

        EXAMPLE::

            sage: digraphs.Circulant(13,[3,5,7])
            Circulant graph ([3, 5, 7]): Digraph on 13 vertices

        TESTS::

            sage: digraphs.Circulant(13,[3,5,7,"hey"])
            Traceback (most recent call last):
            ...
            ValueError: The list must contain only relative integers.
            sage: digraphs.Circulant(3,[3,5,7,3.4])
            Traceback (most recent call last):
            ...
            ValueError: The list must contain only relative integers.
        """
        from sage.graphs.graph_plot import _circle_embedding
        from sage.rings.integer_ring import ZZ

        # Bad input and loops
        loops = False
        for i in integers:
            if not i in ZZ:
                raise ValueError("The list must contain only relative integers.")
            if (i%n) == 0:
                loops = True

        G=DiGraph(n, name="Circulant graph ("+str(integers)+")", loops=loops)

        _circle_embedding(G, list(range(n)))
        for v in range(n):
            G.add_edges([(v,(v+j)%n) for j in integers])

        return G

    def DeBruijn(self, k, n, vertices = 'strings'):
        r"""
        Returns the De Bruijn digraph with parameters `k,n`.

        The De Bruijn digraph with parameters `k,n` is built upon a set of
        vertices equal to the set of words of length `n` from a dictionary of
        `k` letters.

        In this digraph, there is an arc `w_1w_2` if `w_2` can be obtained from
        `w_1` by removing the leftmost letter and adding a new letter at its
        right end.  For more information, see the
        :wikipedia:`Wikipedia article on De Bruijn graph <De_Bruijn_graph>`.

        INPUT:

        - ``k`` -- Two possibilities for this parameter :
              - An integer equal to the cardinality of the alphabet to use, that
                is the degree of the digraph to be produced.
              - An iterable object to be used as the set of letters. The degree
                of the resulting digraph is the cardinality of the set of
                letters.

        - ``n`` -- An integer equal to the length of words in the De Bruijn
          digraph when ``vertices == 'strings'``, and also to the diameter of
          the digraph.

        - ``vertices`` -- 'strings' (default) or 'integers', specifying whether
          the vertices are words build upon an alphabet or integers.

        EXAMPLES::

            sage: db=digraphs.DeBruijn(2,2); db
            De Bruijn digraph (k=2, n=2): Looped digraph on 4 vertices
            sage: db.order()
            4
            sage: db.size()
            8

        TESTS::

            sage: digraphs.DeBruijn(5,0)
            De Bruijn digraph (k=5, n=0): Looped multi-digraph on 1 vertex
            sage: digraphs.DeBruijn(0,0)
            De Bruijn digraph (k=0, n=0): Looped multi-digraph on 0 vertices
        """
        from sage.combinat.words.words import Words
        from sage.rings.integer import Integer

        W = Words(list(range(k)) if isinstance(k, Integer) else k, n)
        A = Words(list(range(k)) if isinstance(k, Integer) else k, 1)
        g = DiGraph(loops=True)

        if vertices == 'strings':
            if n == 0 :
                g.allow_multiple_edges(True)
                v = W[0]
                for a in A:
                    g.add_edge(v.string_rep(), v.string_rep(), a.string_rep())
            else:
                for w in W:
                    ww = w[1:]
                    for a in A:
                        g.add_edge(w.string_rep(), (ww*a).string_rep(), a.string_rep())
        else:
            d = W.size_of_alphabet()
            g = digraphs.GeneralizedDeBruijn(d**n, d)

        g.name( "De Bruijn digraph (k=%s, n=%s)"%(k,n) )
        return g

    def GeneralizedDeBruijn(self, n, d):
        r"""
        Returns the generalized de Bruijn digraph of order `n` and degree `d`.

        The generalized de Bruijn digraph has been defined in [RPK80]_
        [RPK83]_. It has vertex set `V=\{0, 1,..., n-1\}` and there is an arc
        from vertex `u \in V` to all vertices `v \in V` such that
        `v \equiv (u*d + a) \mod{n}` with `0 \leq a < d`.

        When `n = d^{D}`, the generalized de Bruijn digraph is isomorphic to the
        de Bruijn digraph of degree `d` and diameter `D`.

        INPUT:

        - ``n`` -- is the number of vertices of the digraph

        - ``d`` -- is the degree of the digraph

        .. SEEALSO::

            * :meth:`sage.graphs.generic_graph.GenericGraph.is_circulant` --
              checks whether a (di)graph is circulant, and/or returns all
              possible sets of parameters.

        EXAMPLE::

            sage: GB = digraphs.GeneralizedDeBruijn(8, 2)
            sage: GB.is_isomorphic(digraphs.DeBruijn(2, 3), certificate = True)
            (True, {0: '000', 1: '001', 2: '010', 3: '011', 4: '100', 5: '101', 6: '110', 7: '111'})

        TESTS:

        An exception is raised when the degree is less than one::

            sage: G = digraphs.GeneralizedDeBruijn(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: The generalized de Bruijn digraph is defined for degree at least one.

        An exception is raised when the order of the graph is less than one::

            sage: G = digraphs.GeneralizedDeBruijn(0, 2)
            Traceback (most recent call last):
            ...
            ValueError: The generalized de Bruijn digraph is defined for at least one vertex.


        REFERENCES:

        .. [RPK80] \S. M. Reddy, D. K. Pradhan, and J. Kuhl. Directed graphs with
          minimal diameter and maximal connectivity, School Eng., Oakland Univ.,
          Rochester MI, Tech. Rep., July 1980.

        .. [RPK83] \S. Reddy, P. Raghavan, and J. Kuhl. A Class of Graphs for
          Processor Interconnection. *IEEE International Conference on Parallel
          Processing*, pages 154-157, Los Alamitos, Ca., USA, August 1983.
        """
        if n < 1:
            raise ValueError("The generalized de Bruijn digraph is defined for at least one vertex.")
        if d < 1:
            raise ValueError("The generalized de Bruijn digraph is defined for degree at least one.")

        GB = DiGraph(loops = True)
        GB.allow_multiple_edges(True)
        for u in range(n):
            for a in range(u*d, u*d+d):
                GB.add_edge(u, a%n)

        GB.name( "Generalized de Bruijn digraph (n=%s, d=%s)"%(n,d) )
        return GB


    def ImaseItoh(self, n, d):
        r"""
        Returns the digraph of Imase and Itoh of order `n` and degree `d`.

        The digraph of Imase and Itoh has been defined in [II83]_. It has vertex
        set `V=\{0, 1,..., n-1\}` and there is an arc from vertex `u \in V` to
        all vertices `v \in V` such that `v \equiv (-u*d-a-1) \mod{n}` with
        `0 \leq a < d`.

        When `n = d^{D}`, the digraph of Imase and Itoh is isomorphic to the de
        Bruijn digraph of degree `d` and diameter `D`. When `n = d^{D-1}(d+1)`,
        the digraph of Imase and Itoh is isomorphic to the Kautz digraph
        [Kautz68]_ of degree `d` and diameter `D`.

        INPUT:

        - ``n`` -- is the number of vertices of the digraph

        - ``d`` -- is the degree of the digraph

        EXAMPLES::

            sage: II = digraphs.ImaseItoh(8, 2)
            sage: II.is_isomorphic(digraphs.DeBruijn(2, 3), certificate = True)
            (True, {0: '010', 1: '011', 2: '000', 3: '001', 4: '110', 5: '111', 6: '100', 7: '101'})

            sage: II = digraphs.ImaseItoh(12, 2)
            sage: II.is_isomorphic(digraphs.Kautz(2, 3), certificate = True)
            (True, {0: '010', 1: '012', 2: '021', 3: '020', 4: '202', 5: '201', 6: '210', 7: '212', 8: '121', 9: '120', 10: '102', 11: '101'})


        TESTS:

        An exception is raised when the degree is less than one::

            sage: G = digraphs.ImaseItoh(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: The digraph of Imase and Itoh is defined for degree at least one.

        An exception is raised when the order of the graph is less than two::

            sage: G = digraphs.ImaseItoh(1, 2)
            Traceback (most recent call last):
            ...
            ValueError: The digraph of Imase and Itoh is defined for at least two vertices.


        REFERENCE:

        .. [II83] \M. Imase and M. Itoh. A design for directed graphs with
          minimum diameter, *IEEE Trans. Comput.*, vol. C-32, pp. 782-784, 1983.
        """
        if n < 2:
            raise ValueError("The digraph of Imase and Itoh is defined for at least two vertices.")
        if d < 1:
            raise ValueError("The digraph of Imase and Itoh is defined for degree at least one.")

        II = DiGraph(loops = True)
        II.allow_multiple_edges(True)
        for u in range(n):
            for a in range(-u*d-d, -u*d):
                II.add_edge(u, a % n)

        II.name( "Imase and Itoh digraph (n=%s, d=%s)"%(n,d) )
        return II


    def Kautz(self, k, D, vertices = 'strings'):
        r"""
        Returns the Kautz digraph of degree `d` and diameter `D`.

        The Kautz digraph has been defined in [Kautz68]_. The Kautz digraph of
        degree `d` and diameter `D` has `d^{D-1}(d+1)` vertices. This digraph is
        build upon a set of vertices equal to the set of words of length `D`
        from an alphabet of `d+1` letters such that consecutive letters are
        differents. There is an arc from vertex `u` to vertex `v` if `v` can be
        obtained from `u` by removing the leftmost letter and adding a new
        letter, distinct from the rightmost letter of `u`, at the right end.

        The Kautz digraph of degree `d` and diameter `D` is isomorphic to the
        digraph of Imase and Itoh [II83]_ of degree `d` and order
        `d^{D-1}(d+1)`.

        See also the
        :wikipedia:`Wikipedia article on Kautz Graphs <Kautz_graph>`.

        INPUT:

        - ``k`` -- Two possibilities for this parameter :
            - An integer equal to the degree of the digraph to be produced, that
              is the cardinality minus one of the alphabet to use.
            - An iterable object to be used as the set of letters. The degree of
              the resulting digraph is the cardinality of the set of letters
              minus one.

        - ``D`` -- An integer equal to the diameter of the digraph, and also to
              the length of a vertex label when ``vertices == 'strings'``.

        - ``vertices`` -- 'strings' (default) or 'integers', specifying whether
                      the vertices are words build upon an alphabet or integers.


        EXAMPLES::

            sage: K = digraphs.Kautz(2, 3)
            sage: K.is_isomorphic(digraphs.ImaseItoh(12, 2), certificate = True)
            (True,
             {'010': 0,
              '012': 1,
              '020': 3,
              '021': 2,
              '101': 11,
              '102': 10,
              '120': 9,
              '121': 8,
              '201': 5,
              '202': 4,
              '210': 6,
              '212': 7})

            sage: K = digraphs.Kautz([1,'a','B'], 2)
            sage: K.edges()
            [('1B', 'B1', '1'), ('1B', 'Ba', 'a'), ('1a', 'a1', '1'), ('1a', 'aB', 'B'), ('B1', '1B', 'B'), ('B1', '1a', 'a'), ('Ba', 'a1', '1'), ('Ba', 'aB', 'B'), ('a1', '1B', 'B'), ('a1', '1a', 'a'), ('aB', 'B1', '1'), ('aB', 'Ba', 'a')]

            sage: K = digraphs.Kautz([1,'aA','BB'], 2)
            sage: K.edges()
            [('1,BB', 'BB,1', '1'), ('1,BB', 'BB,aA', 'aA'), ('1,aA', 'aA,1', '1'), ('1,aA', 'aA,BB', 'BB'), ('BB,1', '1,BB', 'BB'), ('BB,1', '1,aA', 'aA'), ('BB,aA', 'aA,1', '1'), ('BB,aA', 'aA,BB', 'BB'), ('aA,1', '1,BB', 'BB'), ('aA,1', '1,aA', 'aA'), ('aA,BB', 'BB,1', '1'), ('aA,BB', 'BB,aA', 'aA')]


        TESTS:

        An exception is raised when the degree is less than one::

            sage: G = digraphs.Kautz(0, 2)
            Traceback (most recent call last):
            ...
            ValueError: Kautz digraphs are defined for degree at least one.

            sage: G = digraphs.Kautz(['a'], 2)
            Traceback (most recent call last):
            ...
            ValueError: Kautz digraphs are defined for degree at least one.

        An exception is raised when the diameter of the graph is less than one::

            sage: G = digraphs.Kautz(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: Kautz digraphs are defined for diameter at least one.


        REFERENCE:

        .. [Kautz68] \W. H. Kautz. Bounds on directed (d, k) graphs. Theory of
          cellular logic networks and machines, AFCRL-68-0668, SRI Project 7258,
          Final Rep., pp. 20-28, 1968.
        """
        if D < 1:
            raise ValueError("Kautz digraphs are defined for diameter at least one.")

        from sage.combinat.words.words import Words
        from sage.rings.integer import Integer

        my_alphabet = Words([str(i) for i in range(k+1)] if isinstance(k, Integer) else k, 1)
        if my_alphabet.alphabet().cardinality() < 2:
            raise ValueError("Kautz digraphs are defined for degree at least one.")

        if vertices == 'strings':

            # We start building the set of vertices
            V = [i for i in my_alphabet]
            for i in range(D-1):
                VV = []
                for w in V:
                    VV += [w*a for a in my_alphabet if not w.has_suffix(a) ]
                V = VV

            # We now build the set of arcs
            G = DiGraph()
            for u in V:
                for a in my_alphabet:
                    if not u.has_suffix(a):
                        G.add_edge(u.string_rep(), (u[1:]*a).string_rep(), a.string_rep())

        else:
            d = my_alphabet.size_of_alphabet()-1
            G = digraphs.ImaseItoh( (d+1)*(d**(D-1)), d)

        G.name( "Kautz digraph (k=%s, D=%s)"%(k,D) )
        return G

    def RandomDirectedGN(self, n, kernel=lambda x:x, seed=None):
        """
        Returns a random GN (growing network) digraph with n vertices.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen with a
        preferential attachment model, i.e. probability is proportional to
        degree. The default attachment kernel is a linear function of
        degree. The digraph is always a tree, so in particular it is a
        directed acyclic graph.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``kernel`` - the attachment kernel

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGN(25)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (3, 1), (4, 0), (5, 0), (6, 1), (7, 0), (8, 3), (9, 0), (10, 8), (11, 3), (12, 9), (13, 8), (14, 0), (15, 11), (16, 11), (17, 5), (18, 11), (19, 6), (20, 5), (21, 14), (22, 5), (23, 18), (24, 11)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Organization of Growing
          Random Networks, Phys. Rev. E vol. 63 (2001), p. 066123.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return DiGraph(networkx.gn_graph(n, kernel, seed=seed))

    def RandomDirectedGNC(self, n, seed=None):
        """
        Returns a random GNC (growing network with copying) digraph with n
        vertices.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen with a
        preferential attachment model, i.e. probability is proportional to
        degree. The new vertex is also linked to all of the previously
        added vertex's successors.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNC(25)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (2, 1), (3, 0), (4, 0), (4, 1), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (7, 0), (7, 1), (7, 4), (8, 0), (9, 0), (9, 8), (10, 0), (10, 1), (10, 2), (10, 5), (11, 0), (11, 8), (11, 9), (12, 0), (12, 8), (12, 9), (13, 0), (13, 1), (14, 0), (14, 8), (14, 9), (14, 12), (15, 0), (15, 8), (15, 9), (15, 12), (16, 0), (16, 1), (16, 4), (16, 7), (17, 0), (17, 8), (17, 9), (17, 12), (18, 0), (18, 8), (19, 0), (19, 1), (19, 4), (19, 7), (20, 0), (20, 1), (20, 4), (20, 7), (20, 16), (21, 0), (21, 8), (22, 0), (22, 1), (22, 4), (22, 7), (22, 19), (23, 0), (23, 8), (23, 9), (23, 12), (23, 14), (24, 0), (24, 8), (24, 9), (24, 12), (24, 15)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Network Growth by
          Copying, Phys. Rev. E vol. 71 (2005), p. 036118.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return DiGraph(networkx.gnc_graph(n, seed=seed))

    def RandomDirectedGNP(self, n, p, loops = False, seed = None):
        r"""
        Returns a random digraph on `n` nodes. Each edge is inserted
        independently with probability `p`.

        INPUT:

        - ``n`` -- number of nodes of the digraph

        - ``p`` -- probability of an edge

        - ``loops`` -- is a boolean set to True if the random digraph may have
          loops, and False (default) otherwise.

        - ``seed`` -- integer seed for random number generator (default=None).

        REFERENCES:

        .. [1] \P. Erdos and A. Renyi, On Random Graphs, Publ.  Math. 6, 290
               (1959).

        .. [2] \E. N. Gilbert, Random Graphs, Ann. Math.  Stat., 30, 1141 (1959).


        PLOTTING: When plotting, this graph will use the default spring-layout
        algorithm, unless a position dictionary is specified.

        EXAMPLE::

            sage: set_random_seed(0)
            sage: D = digraphs.RandomDirectedGNP(10, .2)
            sage: D.num_verts()
            10
            sage: D.edges(labels=False)
            [(1, 0), (1, 2), (3, 6), (3, 7), (4, 5), (4, 7), (4, 8), (5, 2), (6, 0), (7, 2), (8, 1), (8, 9), (9, 4)]
        """
        from sage.graphs.graph_generators_pyx import RandomGNP
        if 0.0 > p or 1.0 < p:
            raise ValueError("The probability p must be in [0..1].")

        if seed is None:
            seed = current_randstate().long_seed()

        return RandomGNP(n, p, directed = True, loops = loops)

    def RandomDirectedGNM(self, n, m, loops = False):
        r"""
        Returns a random labelled digraph on `n` nodes and `m` arcs.

        INPUT:

        - ``n`` (integer) -- number of vertices.

        - ``m`` (integer) -- number of edges.

        - ``loops`` (boolean) -- whether to allow loops (set to ``False`` by
          default).

        PLOTTING: When plotting, this graph will use the default spring-layout
        algorithm, unless a position dictionary is specified.

        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNM(10, 5)
            sage: D.num_verts()
            10
            sage: D.edges(labels=False)
            [(0, 3), (1, 5), (5, 1), (7, 0), (8, 5)]

        With loops::

            sage: D = digraphs.RandomDirectedGNM(10, 100, loops = True)
            sage: D.num_verts()
            10
            sage: D.loops()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None), (4, 4, None), (5, 5, None), (6, 6, None), (7, 7, None), (8, 8, None), (9, 9, None)]

        TESTS::

            sage: digraphs.RandomDirectedGNM(10,-3)
            Traceback (most recent call last):
            ...
            ValueError: The number of edges must satisfy 0<= m <= n(n-1) when no loops are allowed, and 0<= m <= n^2 otherwise.

            sage: digraphs.RandomDirectedGNM(10,100)
            Traceback (most recent call last):
            ...
            ValueError: The number of edges must satisfy 0<= m <= n(n-1) when no loops are allowed, and 0<= m <= n^2 otherwise.
        """
        n, m = int(n), int(m)

        # The random graph is built by drawing randomly and uniformly two
        # integers u,v, and adding the corresponding edge if it does not exist,
        # as many times as necessary.

        # When the graph is dense, we actually compute its complement. This will
        # prevent us from drawing the same pair u,v too many times.

        from sage.misc.prandom import _pyrand
        rand = _pyrand()
        D = DiGraph(n, loops = loops)

        # Ensuring the parameters n,m make sense.
        #
        # If the graph is dense, we actually want to build its complement. We
        # update m accordingly.

        good_input = True
        is_dense = False

        if m < 0:
            good_input = False

        if loops:
            if m > n*n:
                good_input = False
            elif m > n*n/2:
                is_dense = True
                m = n*n - m

        else:
            if m > n*(n-1):
                good_input = False
            elif m > n*(n-1)/2:
                is_dense = True
                m = n*(n-1) - m

        if not good_input:
            raise ValueError("The number of edges must satisfy 0<= m <= n(n-1) when no loops are allowed, and 0<= m <= n^2 otherwise.")

        # When the given number of edges defines a density larger than 1/2, it
        # should be faster to compute the complement of the graph (less edges to
        # generate), then to return its complement. This being said, the
        # .complement() method for sparse graphs is very slow at the moment.

        # Similarly, it is faster to test whether a pair belongs to a dictionary
        # than to test the adjacency of two vertices in a graph. For these
        # reasons, the following code mainly works on dictionaries.

        adj = dict( (i, dict()) for i in range(n) )

        # We fill the dictionary structure, but add the corresponding edge in
        # the graph only if is_dense is False. If it is true, we will add the
        # edges in a second phase.


        while m > 0:

            # It is better to obtain random numbers this way than by calling the
            # randint or randrange method. This, because they are very expensive
            # when trying to compute MANY random integers, and because the
            # following lines is precisely what they do anyway, after checking
            # their parameters are correct.

            u=int(rand.random()*n)
            v=int(rand.random()*n)

            if (u != v or loops) and (not v in adj[u]):
                adj[u][v] = 1
                m -= 1
                if not is_dense:
                    D.add_edge(u,v)

        # If is_dense is True, it means the graph has not been built. We fill D
        # with the complement of the edges stored in the adj dictionary

        if is_dense:
            for u in range(n):
                for v in range(n):
                    if ((u != v) or loops) and (not (v in adj[u])):
                        D.add_edge(u,v)

        return D

    def RandomDirectedGNR(self, n, p, seed=None):
        """
        Returns a random GNR (growing network with redirection) digraph
        with n vertices and redirection probability p.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen uniformly.
        With probability p, the arc is instead redirected to the successor
        vertex. The digraph is always a tree.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``p`` - redirection probability

        -  ``seed`` - for the random number generator.


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNR(25, .2)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (2, 1), (3, 0), (4, 0), (4, 1), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (7, 0), (7, 1), (7, 4), (8, 0), (9, 0), (9, 8), (10, 0), (10, 1), (10, 2), (10, 5), (11, 0), (11, 8), (11, 9), (12, 0), (12, 8), (12, 9), (13, 0), (13, 1), (14, 0), (14, 8), (14, 9), (14, 12), (15, 0), (15, 8), (15, 9), (15, 12), (16, 0), (16, 1), (16, 4), (16, 7), (17, 0), (17, 8), (17, 9), (17, 12), (18, 0), (18, 8), (19, 0), (19, 1), (19, 4), (19, 7), (20, 0), (20, 1), (20, 4), (20, 7), (20, 16), (21, 0), (21, 8), (22, 0), (22, 1), (22, 4), (22, 7), (22, 19), (23, 0), (23, 8), (23, 9), (23, 12), (23, 14), (24, 0), (24, 8), (24, 9), (24, 12), (24, 15)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Organization of Growing
          Random Networks, Phys. Rev. E vol. 63 (2001), p. 066123.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return DiGraph(networkx.gnc_graph(n, seed=seed))

    def RandomSemiComplete(self, n):
        r"""
        Return a random semi-complete digraph on `n` vertices.

        A directed graph `G=(V,E)` is *semi-complete* if for any pair of
        vertices `u` and `v`, there is *at least* one arc between them.

        To generate randomly a semi-complete digraph, we have to ensure, for any
        pair of distinct vertices `u` and `v`, that with probability `1/3` we
        have only arc `uv`, with probability `1/3` we have only arc `vu`, and
        with probability `1/3` we have both arc `uv` and arc `vu`. We do so by
        selecting a random integer `coin` in `[1,3]`. When `coin==1` we select
        only arc `uv`, when `coin==3` we select only arc `vu`, and when
        `coin==2` we select both arcs. In other words, we select arc `uv` when
        `coin\leq 2` and arc `vu` when `coin\geq 2`.

        INPUT:

        - ``n`` (integer) -- the number of nodes

        .. SEEALSO::

            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.Complete`

            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.RandomTournament`

        EXAMPLES::

            sage: SC = digraphs.RandomSemiComplete(10); SC
            Random Semi-Complete digraph: Digraph on 10 vertices
            sage: SC.size() >= binomial(10, 2)
            True
            sage: digraphs.RandomSemiComplete(-1)
            Traceback (most recent call last):
            ...
            ValueError: The number of vertices cannot be strictly negative!
        """
        G = DiGraph(n, name="Random Semi-Complete digraph")

        # For each pair u,v we choose a randon number ``coin`` in [1,3].
        # We select edge `(u,v)` if `coin==1` or `coin==2`.
        # We select edge `(v,u)` if `coin==2` or `coin==3`.
        import itertools
        from sage.misc.prandom import randint
        for u,v in itertools.combinations(range(n), 2):
            coin = randint(1,3)
            if coin<=2:
                G.add_edge(u,v)
            if coin>=2:
                G.add_edge(v,u)

        if n:
            from sage.graphs.graph_plot import _circle_embedding
            _circle_embedding(G, list(range(n)))

        return G

################################################################################
#   DiGraph Iterators
################################################################################

    def __call__(self, vertices, property=lambda x: True, augment='edges', size=None, implementation='c_graph', sparse=True):
        """
        Accesses the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives.

        INPUT:


        -  ``vertices`` - natural number

        -  ``property`` - any property to be tested on digraphs
           before generation.

        -  ``augment`` - choices:

        -  ``'vertices'`` - augments by adding a vertex, and
           edges incident to that vertex. In this case, all digraphs on up to
           n=vertices are generated. If for any digraph G satisfying the
           property, every subgraph, obtained from G by deleting one vertex
           and only edges incident to that vertex, satisfies the property,
           then this will generate all digraphs with that property. If this
           does not hold, then all the digraphs generated will satisfy the
           property, but there will be some missing.

        -  ``'edges'`` - augments a fixed number of vertices by
           adding one edge In this case, all digraphs on exactly n=vertices
           are generated. If for any graph G satisfying the property, every
           subgraph, obtained from G by deleting one edge but not the vertices
           incident to that edge, satisfies the property, then this will
           generate all digraphs with that property. If this does not hold,
           then all the digraphs generated will satisfy the property, but
           there will be some missing.

        -  ``implementation`` - which underlying implementation to use (see DiGraph?)

        -  ``sparse`` - ignored if implementation is not ``c_graph``

        EXAMPLES: Print digraphs on 2 or less vertices.

        ::

            sage: for D in digraphs(2, augment='vertices'):
            ....:     print(D)
            Digraph on 0 vertices
            Digraph on 1 vertex
            Digraph on 2 vertices
            Digraph on 2 vertices
            Digraph on 2 vertices

        Print digraphs on 3 vertices.

        ::

            sage: for D in digraphs(3):
            ....:     print(D)
            Digraph on 3 vertices
            Digraph on 3 vertices
            ...
            Digraph on 3 vertices
            Digraph on 3 vertices

        For more examples, see the class level documentation, or type

        ::

            sage: digraphs? # not tested

        REFERENCE:

        - Brendan D. McKay, Isomorph-Free Exhaustive generation.
          Journal of Algorithms Volume 26, Issue 2, February 1998,
          pages 306-324.
        """
        if size is not None:
            extra_property = lambda x: x.size() == size
        else:
            extra_property = lambda x: True
        if augment == 'vertices':
            from sage.graphs.graph_generators import canaug_traverse_vert
            g = DiGraph(implementation=implementation, sparse=sparse)
            for gg in canaug_traverse_vert(g, [], vertices, property, dig=True, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        elif augment == 'edges':
            from sage.graphs.graph_generators import canaug_traverse_edge
            g = DiGraph(vertices, implementation=implementation, sparse=sparse)
            gens = []
            for i in range(vertices-1):
                gen = list(range(i))
                gen.append(i+1); gen.append(i)
                gen += list(range(i + 2, vertices))
                gens.append(gen)
            for gg in canaug_traverse_edge(g, gens, property, dig=True, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        else:
            raise NotImplementedError()


# Easy access to the graph generators from the command line:
digraphs = DiGraphGenerators()




