# -*- coding: utf-8 -*-
r"""
Common Graphs (Graph Generators)

All graphs in Sage can be built through the ``graphs`` object. In order to
build a complete graph on 15 elements, one can do::

    sage: g = graphs.CompleteGraph(15)

To get a path with 4 vertices, and the house graph::

    sage: p = graphs.PathGraph(4)
    sage: h = graphs.HouseGraph()

More interestingly, one can get the list of all graphs that Sage knows how to
build by typing ``graphs.`` in Sage and then hitting tab.
"""

# This method appends a list of methods to the doc as a 3xN table.

# Here's the point :
#
# we just have to insert the method's name in this file to add it to
# the tab, and in exchange the doc contains a table of width 3 with
# all methods listed, so that the reading order is Column1, then
# Column2, then Column3. Doing this by hand is hell with Sphinx when
# you need to insert a new method inside of the list !

def __append_to_doc(methods):
    global __doc__
    __doc__ += ("\n.. csv-table::\n"
    "    :class: contentstable\n"
    "    :widths: 33, 33, 33\n"
    "    :delim: |\n\n")

    h = (len(methods)+2)//3
    # Reorders the list of methods for horizontal reading, the only one Sphinx understands
    reordered_methods = [0]*3*h
    for i, m in enumerate(methods):
        reordered_methods[3*(i%h)+(i//h)] = m
    methods = reordered_methods

    # Adding the list to the __doc__ string
    wrap_name = lambda x : ":meth:`"+str(x)+" <GraphGenerators."+str(x)+">`" if x else ""
    while methods:
        a = methods.pop(0)
        b = methods.pop(0)
        c = methods.pop(0)
        __doc__ += "    "+wrap_name(a)+" | "+wrap_name(b)+" | "+wrap_name(c)+"\n"

__doc__ += """
**Basic structures**
"""

__append_to_doc(
    ["BullGraph",
     "ButterflyGraph",
     "CircularLadderGraph",
     "ClawGraph",
     "CycleGraph",
     "CompleteBipartiteGraph",
     "CompleteGraph",
     "CompleteMultipartiteGraph",
     "DiamondGraph",
     "EmptyGraph",
     "Grid2dGraph",
     "GridGraph",
     "HouseGraph",
     "HouseXGraph",
     "LadderGraph",
     "LollipopGraph",
     "PathGraph",
     "StarGraph",
     "ToroidalGrid2dGraph",
     "Toroidal6RegularGrid2dGraph"]
    )

__doc__ += """
**Small Graphs**

A small graph is just a single graph and has no parameter influencing
the number of edges or vertices.
"""

__append_to_doc(
    ["Balaban10Cage",
     "Balaban11Cage",
     "BidiakisCube",
     "BiggsSmithGraph",
     "BlanusaFirstSnarkGraph",
     "BlanusaSecondSnarkGraph",
     "BrinkmannGraph",
     "BrouwerHaemersGraph",
     "BuckyBall",
     "CameronGraph",
     "ChvatalGraph",
     "ClebschGraph",
     "CoxeterGraph",
     "DesarguesGraph",
     "DoubleStarSnark",
     "DurerGraph",
     "DyckGraph",
     "EllinghamHorton54Graph",
     "EllinghamHorton78Graph",
     "ErreraGraph",
     "FlowerSnark",
     "FolkmanGraph",
     "FosterGraph",
     "FranklinGraph",
     "FruchtGraph",
     "GoldnerHararyGraph",
     "GrayGraph",
     "GrotzschGraph",
     "HallJankoGraph",
     "HarriesGraph",
     "HarriesWongGraph",
     "HeawoodGraph",
     "HerschelGraph",
     "HigmanSimsGraph",
     "HoffmanGraph",
     "HoffmanSingletonGraph",
     "HoltGraph",
     "HortonGraph",
     "KittellGraph",
     "KrackhardtKiteGraph",
     "LjubljanaGraph",
     "M22Graph",
     "MarkstroemGraph",
     "McGeeGraph",
     "McLaughlinGraph",
     "MeredithGraph",
     "MoebiusKantorGraph",
     "MoserSpindle",
     "NauruGraph",
     "PappusGraph",
     "PoussinGraph",
     "PetersenGraph",
     "RobertsonGraph",
     "SchlaefliGraph",
     "ShrikhandeGraph",
     "SimsGewirtzGraph",
     "SousselierGraph",
     "SylvesterGraph",
     "SzekeresSnarkGraph",
     "ThomsenGraph",
     "TietzeGraph",
     "Tutte12Cage",
     "TutteCoxeterGraph",
     "TutteGraph",
     "WagnerGraph",
     "WatkinsSnarkGraph",
     "WellsGraph",
     "WienerArayaGraph"])

__doc__ += """
*Platonic solids* (ordered ascending by number of vertices)
"""

__append_to_doc(
    ["TetrahedralGraph",
     "OctahedralGraph",
     "HexahedralGraph",
     "IcosahedralGraph",
     "DodecahedralGraph"])

__doc__ += """
**Families of graphs**

A family of graph is an infinite set of graphs which can be indexed by fixed
number of parameters, e.g. two integer parameters. (A method whose name starts
with a small letter does not return a single graph object but a graph iterator
or a list of graphs or ...)
"""

__append_to_doc(
    ["BalancedTree",
     "BarbellGraph",
     "BubbleSortGraph",
     "CirculantGraph",
     "cospectral_graphs",
     "CubeGraph",
     "DorogovtsevGoltsevMendesGraph",
     "FibonacciTree",
     "FoldedCubeGraph",
     "FriendshipGraph",
     "fullerenes",
     "fusenes",
     "FuzzyBallGraph",
     "GeneralizedPetersenGraph",
     "HanoiTowerGraph",
     "HararyGraph",
     "HyperStarGraph",
     "JohnsonGraph",
     "KneserGraph",
     "LCFGraph",
     "line_graph_forbidden_subgraphs",
     "MycielskiGraph",
     "MycielskiStep",
     "NKStarGraph",
     "NStarGraph",
     "OddGraph",
     "PaleyGraph",
     "RingedTree",
     "SymplecticGraph",
     "trees",
     "WheelGraph"])

__doc__ += """
*Chessboard Graphs*
"""

__append_to_doc(
    ["BishopGraph",
     "KingGraph",
     "KnightGraph",
     "QueenGraph",
     "RookGraph"])

__doc__ += """
**Intersection graphs**

These graphs are generated by geometric representations. The objects of
the representation correspond to the graph vertices and the intersections
of objects yield the graph edges.
"""

__append_to_doc(
    ["IntervalGraph",
     "PermutationGraph",
     "ToleranceGraph"])

__doc__ += """
**Random graphs**
"""

__append_to_doc(
    ["RandomBarabasiAlbert",
     "RandomBipartite",
     "RandomBoundedToleranceGraph",
     "RandomGNM",
     "RandomGNP",
     "RandomHolmeKim",
     "RandomIntervalGraph",
     "RandomLobster",
     "RandomNewmanWattsStrogatz",
     "RandomRegular",
     "RandomShell",
     "RandomToleranceGraph",
     "RandomTree",
     "RandomTreePowerlaw"])

__doc__ += """
**Graphs with a given degree sequence**
"""

__append_to_doc(
    ["DegreeSequence",
     "DegreeSequenceBipartite",
     "DegreeSequenceConfigurationModel",
     "DegreeSequenceExpected",
     "DegreeSequenceTree"])

__doc__ += """
**Miscellaneous**
"""

__append_to_doc(
    ["WorldMap"]
    )

__doc__ += """

AUTHORS:

- Robert Miller (2006-11-05): initial version, empty, random, petersen

- Emily Kirkman (2006-11-12): basic structures, node positioning for
  all constructors

- Emily Kirkman (2006-11-19): docstrings, examples

- William Stein (2006-12-05): Editing.

- Robert Miller (2007-01-16): Cube generation and plotting

- Emily Kirkman (2007-01-16): more basic structures, docstrings

- Emily Kirkman (2007-02-14): added more named graphs

- Robert Miller (2007-06-08-11): Platonic solids, random graphs,
  graphs with a given degree sequence, random directed graphs

- Robert Miller (2007-10-24): Isomorph free exhaustive generation

- Nathann Cohen (2009-08-12): WorldMap

- Michael Yurko (2009-9-01): added hyperstar, (n,k)-star, n-star, and
  bubblesort graphs

- Anders Jonsson (2009-10-15): added generalized Petersen graphs

- Harald Schilly and Yann Laigle-Chapuy (2010-03-24): added Fibonacci Tree

- Jason Grout (2010-06-04): cospectral_graphs

- Edward Scheinerman (2010-08-11): RandomTree

- Ed Scheinerman (2010-08-21): added Grotzsch graph and Mycielski graphs

- Minh Van Nguyen (2010-11-26): added more named graphs

- Keshav Kini (2011-02-16): added Shrikhande and Dyck graphs

- David Coudert (2012-02-10): new RandomGNP generator

- David Coudert (2012-08-02): added chessboard graphs: Queen, King,
  Knight, Bishop, and Rook graphs

- Nico Van Cleemput (2013-05-26): added fullerenes

- Nico Van Cleemput (2013-07-01): added benzenoids

- Birk Eisermann (2013-07-29): new section 'intersection graphs',
  added (random, bounded) tolerance graphs


Functions and methods
---------------------
"""

###########################################################################

#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Python standard library

# import from Sage library
import graph

class GraphGenerators():
    r"""
    A class consisting of constructors for several common graphs, as well as
    orderly generation of isomorphism class representatives. See the
    :mod:`module's help <sage.graphs.graph_generators>` for a list of supported
    constructors.

    A list of all graphs and graph structures (other than isomorphism class
    representatives) in this database is available via tab completion. Type
    "graphs." and then hit the tab key to see which graphs are available.

    The docstrings include educational information about each named
    graph with the hopes that this class can be used as a reference.

    For all the constructors in this class (except the octahedral,
    dodecahedral, random and empty graphs), the position dictionary is
    filled to override the spring-layout algorithm.


    ORDERLY GENERATION::

        graphs(vertices, property=lambda x: True, augment='edges', size=None)

    This syntax accesses the generator of isomorphism class
    representatives. Iterates over distinct, exhaustive
    representatives.

    Also: see the use of the optional nauty package for generating graphs
    at the :meth:`nauty_geng` method.

    INPUT:

    - ``vertices`` -- natural number.

    - ``property`` -- (default: ``lambda x: True``) any property to be
      tested on graphs before generation, but note that in general the
      graphs produced are not the same as those produced by using the
      property function to filter a list of graphs produced by using
      the ``lambda x: True`` default. The generation process assumes
      the property has certain characteristics set by the ``augment``
      argument, and only in the case of inherited properties such that
      all subgraphs of the relevant kind (for ``augment='edges'`` or
      ``augment='vertices'``) of a graph with the property also
      possess the property will there be no missing graphs.  (The
      ``property`` argument is ignored if ``degree_sequence`` is
      specified.)

    - ``augment`` -- (default: ``'edges'``) possible values:

      - ``'edges'`` -- augments a fixed number of vertices by
        adding one edge. In this case, all graphs on exactly ``n=vertices`` are
        generated. If for any graph G satisfying the property, every
        subgraph, obtained from G by deleting one edge but not the vertices
        incident to that edge, satisfies the property, then this will
        generate all graphs with that property. If this does not hold, then
        all the graphs generated will satisfy the property, but there will
        be some missing.

      - ``'vertices'`` -- augments by adding a vertex and
        edges incident to that vertex. In this case, all graphs up to
        ``n=vertices`` are generated. If for any graph G satisfying the
        property, every subgraph, obtained from G by deleting one vertex
        and only edges incident to that vertex, satisfies the property,
        then this will generate all graphs with that property. If this does
        not hold, then all the graphs generated will satisfy the property,
        but there will be some missing.

    - ``size`` -- (default: ``None``) the size of the graph to be generated.

    - ``degree_sequence`` -- (default: ``None``) a sequence of non-negative integers,
      or ``None``. If specified, the generated graphs will have these
      integers for degrees. In this case, property and size are both
      ignored.

    - ``loops`` -- (default: ``False``) whether to allow loops in the graph
      or not.

    - ``implementation`` -- (default: ``'c_graph'``) which underlying
      implementation to use (see ``Graph?``).

    - ``sparse`` -- (default: ``True``) ignored if implementation is not
      ``'c_graph'``.

    - ``copy`` (boolean) -- If set to ``True`` (default)
      this method makes copies of the graphs before returning
      them. If set to ``False`` the method returns the graph it
      is working on. The second alternative is faster, but modifying
      any of the graph instances returned by the method may break
      the function's behaviour, as it is using these graphs to
      compute the next ones : only use ``copy_graph = False`` when
      you stick to *reading* the graphs returned.

    EXAMPLES:

    Print graphs on 3 or less vertices::

        sage: for G in graphs(3, augment='vertices'):
        ...    print G
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Note that we can also get graphs with underlying Cython implementation::

        sage: for G in graphs(3, augment='vertices', implementation='c_graph'):
        ...    print G
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Print graphs on 3 vertices.

    ::

        sage: for G in graphs(3):
        ...    print G
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices

    Generate all graphs with 5 vertices and 4 edges.

    ::

        sage: L = graphs(5, size=4)
        sage: len(list(L))
        6

    Generate all graphs with 5 vertices and up to 4 edges.

    ::

        sage: L = list(graphs(5, lambda G: G.size() <= 4))
        sage: len(L)
        14
        sage: graphs_list.show_graphs(L) # long time

    Generate all graphs with up to 5 vertices and up to 4 edges.

    ::

        sage: L = list(graphs(5, lambda G: G.size() <= 4, augment='vertices'))
        sage: len(L)
        31
        sage: graphs_list.show_graphs(L)              # long time

    Generate all graphs with degree at most 2, up to 6 vertices.

    ::

        sage: property = lambda G: ( max([G.degree(v) for v in G] + [0]) <= 2 )
        sage: L = list(graphs(6, property, augment='vertices'))
        sage: len(L)
        45

    Generate all bipartite graphs on up to 7 vertices: (see
    http://oeis.org/classic/A033995)

    ::

        sage: L = list( graphs(7, lambda G: G.is_bipartite(), augment='vertices') )
        sage: [len([g for g in L if g.order() == i]) for i in [1..7]]
        [1, 2, 3, 7, 13, 35, 88]

    Generate all bipartite graphs on exactly 7 vertices::

        sage: L = list( graphs(7, lambda G: G.is_bipartite()) )
        sage: len(L)
        88

    Generate all bipartite graphs on exactly 8 vertices::

        sage: L = list( graphs(8, lambda G: G.is_bipartite()) ) # long time
        sage: len(L)                                            # long time
        303

    Remember that the property argument does not behave as a filter,
    except for appropriately inheritable properties::

        sage: property = lambda G: G.is_vertex_transitive()
        sage: len(list(graphs(4, property)))
        1
        sage: len(filter(property, graphs(4)))
        4
        sage: property = lambda G: G.is_bipartite()
        sage: len(list(graphs(4, property)))
        7
        sage: len(filter(property, graphs(4)))
        7

    Generate graphs on the fly: (see http://oeis.org/classic/A000088)

    ::

        sage: for i in range(0, 7):
        ...    print len(list(graphs(i)))
        1
        1
        2
        4
        11
        34
        156

    Generate all simple graphs, allowing loops: (see
    http://oeis.org/classic/A000666)

    ::

        sage: L = list(graphs(5,augment='vertices',loops=True))               # long time
        sage: for i in [0..5]: print i, len([g for g in L if g.order() == i]) # long time
        0 1
        1 2
        2 6
        3 20
        4 90
        5 544

    Generate all graphs with a specified degree sequence (see
    http://oeis.org/classic/A002851)::

        sage: for i in [4,6,8]:  # long time (4s on sage.math, 2012)
        ...       print i, len([g for g in graphs(i, degree_sequence=[3]*i) if g.is_connected()])
        4 1
        6 2
        8 5
        sage: for i in [4,6,8]:  # long time (7s on sage.math, 2012)
        ...       print i, len([g for g in graphs(i, augment='vertices', degree_sequence=[3]*i) if g.is_connected()])
        4 1
        6 2
        8 5

    ::

        sage: print 10, len([g for g in graphs(10,degree_sequence=[3]*10) if g.is_connected()]) # not tested
        10 19

    Make sure that the graphs are really independent and the generator
    survives repeated vertex removal (trac 8458)::

        sage: for G in graphs(3):
        ...       G.delete_vertex(0)
        ...       print(G.order())
        2
        2
        2
        2

    REFERENCE:

    - Brendan D. McKay, Isomorph-Free Exhaustive generation.  *Journal
      of Algorithms*, Volume 26, Issue 2, February 1998, pages 306-324.
    """

###########################################################################
#   Graph Iterators
###########################################################################

    def __call__(self, vertices=None, property=lambda x: True, augment='edges',
        size=None, deg_seq=None, degree_sequence=None, loops=False, implementation='c_graph',
        sparse=True, copy = True):
        """
        Accesses the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives. See the docstring
        of this class for full documentation.

        EXAMPLES:

        Print graphs on 3 or less vertices::

            sage: for G in graphs(3, augment='vertices'):
            ...    print G
            Graph on 0 vertices
            Graph on 1 vertex
            Graph on 2 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 2 vertices
            Graph on 3 vertices

        ::

            sage: for g in graphs():
            ...    if g.num_verts() > 3: break
            ...    print g
            Graph on 0 vertices
            Graph on 1 vertex
            Graph on 2 vertices
            Graph on 2 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 3 vertices
            Graph on 3 vertices

        For more examples, see the class level documentation, or type::

            sage: graphs? # not tested

        REFERENCE:

        - Brendan D. McKay, Isomorph-Free Exhaustive generation.
          Journal of Algorithms Volume 26, Issue 2, February 1998,
          pages 306-324.
        """
        from sage.graphs.all import Graph
        from sage.misc.superseded import deprecation
        from copy import copy as copyfun

        if deg_seq is not None:
            deprecation(11927, "The argument name deg_seq is deprecated. It will be "
                        "removed in a future release of Sage. So, please use "
                        "degree_sequence instead.")
        if degree_sequence is None:
            degree_sequence=deg_seq
        if degree_sequence is not None:
            if vertices is None:
                raise NotImplementedError
            if len(degree_sequence) != vertices or sum(degree_sequence)%2 or sum(degree_sequence) > vertices*(vertices-1):
                raise ValueError("Invalid degree sequence.")
            degree_sequence = sorted(degree_sequence)
            if augment == 'edges':
                property = lambda x: all([degree_sequence[i] >= d for i,d in enumerate(sorted(x.degree()))])
                extra_property = lambda x: degree_sequence == sorted(x.degree())
            else:
                property = lambda x: all([degree_sequence[i] >= d for i,d in enumerate(sorted(x.degree() + [0]*(vertices-x.num_verts()) ))])
                extra_property = lambda x: x.num_verts() == vertices and degree_sequence == sorted(x.degree())
        elif size is not None:
            extra_property = lambda x: x.size() == size
        else:
            extra_property = lambda x: True

        if augment == 'vertices':
            if vertices is None:
                raise NotImplementedError
            g = Graph(loops=loops, implementation=implementation, sparse=sparse)
            for gg in canaug_traverse_vert(g, [], vertices, property, loops=loops, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield copyfun(gg) if copy else gg
        elif augment == 'edges':
            if vertices is None:
                from sage.rings.all import Integer
                vertices = Integer(0)
                while True:
                    for g in self(vertices, loops=loops, implementation=implementation, sparse=sparse):
                        yield copyfun(g) if copy else g
                    vertices += 1
            g = Graph(vertices, loops=loops, implementation=implementation, sparse=sparse)
            gens = []
            for i in range(vertices-1):
                gen = range(i)
                gen.append(i+1); gen.append(i)
                gen += range(i+2, vertices)
                gens.append(gen)
            for gg in canaug_traverse_edge(g, gens, property, loops=loops, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield copyfun(gg) if copy else gg
        else:
            raise NotImplementedError


    def nauty_geng(self, options="", debug=False):
        r"""
        Returns a generator which creates graphs from nauty's geng program.

        .. note::

            Due to license restrictions, the nauty package is distributed
            as a Sage optional package.  At a system command line, execute
            ``sage -i nauty`` to see the nauty license and install the
            package.

        INPUT:

        - ``options`` - a string passed to  geng  as if it was run at
          a system command line. At a minimum, you *must* pass the
          number of vertices you desire.  Sage expects the graphs to be
          in nauty's "graph6" format, do not set an option to change
          this default or results will be unpredictable.

        - ``debug`` - default: ``False`` - if ``True`` the first line of
          geng's output to standard error is captured and the first call
          to the generator's ``next()`` function will return this line
          as a string.  A line leading with ">A" indicates a successful
          initiation of the program with some information on the arguments,
          while a line beginning with ">E" indicates an error with the input.

        The possible options, obtained as output of ``geng --help``::

                 n    : the number of vertices
            mine:maxe : a range for the number of edges
                        #:0 means '# or more' except in the case 0:0
              res/mod : only generate subset res out of subsets 0..mod-1

                -c    : only write connected graphs
                -C    : only write biconnected graphs
                -t    : only generate triangle-free graphs
                -f    : only generate 4-cycle-free graphs
                -b    : only generate bipartite graphs
                            (-t, -f and -b can be used in any combination)
                -m    : save memory at the expense of time (only makes a
                            difference in the absence of -b, -t, -f and n <= 28).
                -d#   : a lower bound for the minimum degree
                -D#   : a upper bound for the maximum degree
                -v    : display counts by number of edges
                -l    : canonically label output graphs

                -q    : suppress auxiliary output (except from -v)

        Options which cause geng to use an output format different
        than the graph6 format are not listed above (-u, -g, -s, -y, -h)
        as they will confuse the creation of a Sage graph.  The res/mod
        option can be useful when using the output in a routine run
        several times in parallel.

        OUTPUT:

        A generator which will produce the graphs as Sage graphs.
        These will be simple graphs: no loops, no multiple edges, no
        directed edges.

        .. SEEALSO::

            :meth:`Graph.is_strongly_regular` -- tests whether a graph is
            strongly regular and/or returns its parameters.

        EXAMPLES:

        The generator can be used to construct graphs for testing,
        one at a time (usually inside a loop).  Or it can be used to
        create an entire list all at once if there is sufficient memory
        to contain it.  ::

            sage: gen = graphs.nauty_geng("2") # optional nauty
            sage: gen.next() # optional nauty
            Graph on 2 vertices
            sage: gen.next() # optional nauty
            Graph on 2 vertices
            sage: gen.next() # optional nauty
            Traceback (most recent call last):
            ...
            StopIteration: Exhausted list of graphs from nauty geng

        A list of all graphs on 7 vertices.  This agrees with
        Sloane's OEIS sequence A000088.  ::

            sage: gen = graphs.nauty_geng("7") # optional nauty
            sage: len(list(gen))  # optional nauty
            1044

        A list of just the connected graphs on 7 vertices.  This agrees with
        Sloane's OEIS sequence A001349.  ::

            sage: gen = graphs.nauty_geng("7 -c") # optional nauty
            sage: len(list(gen))  # optional nauty
            853

        The ``debug`` switch can be used to examine geng's reaction
        to the input in the ``options`` string.  We illustrate success.
        (A failure will be a string beginning with ">E".)  Passing the
        "-q" switch to geng will supress the indicator of a
        successful initiation.  ::

            sage: gen = graphs.nauty_geng("4", debug=True) # optional nauty
            sage: print gen.next() # optional nauty
            >A nauty-geng -d0D3 n=4 e=0-6
        """
        import subprocess
        from sage.misc.package import is_package_installed
        if not is_package_installed("nauty"):
            raise TypeError("the optional nauty package is not installed")
        sp = subprocess.Popen("nauty-geng {0}".format(options), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)
        if debug:
            yield sp.stderr.readline()
        gen = sp.stdout
        while True:
            try:
                s = gen.next()
            except StopIteration:
                raise StopIteration("Exhausted list of graphs from nauty geng")
            G = graph.Graph(s[:-1], format='graph6')
            yield G


    def cospectral_graphs(self, vertices, matrix_function=lambda g: g.adjacency_matrix(), graphs=None):
        r"""
        Find all sets of graphs on ``vertices`` vertices (with
        possible restrictions) which are cospectral with respect to a
        constructed matrix.

        INPUT:

        - ``vertices`` - The number of vertices in the graphs to be tested

        - ``matrix_function`` - A function taking a graph and giving back
          a matrix.  This defaults to the adjacency matrix.  The spectra
          examined are the spectra of these matrices.

        - ``graphs`` - One of three things:

           - ``None`` (default) - test all graphs having ``vertices``
             vertices

           - a function taking a graph and returning ``True`` or ``False``
             - test only the graphs on ``vertices`` vertices for which
             the function returns ``True``

           - a list of graphs (or other iterable object) - these graphs
             are tested for cospectral sets.  In this case,
             ``vertices`` is ignored.

        OUTPUT:

           A list of lists of graphs.  Each sublist will be a list of
           cospectral graphs (lists of cadinality 1 being omitted).


        .. SEEALSO::

            :meth:`Graph.is_strongly_regular` -- tests whether a graph is
            strongly regular and/or returns its parameters.

        EXAMPLES::

            sage: g=graphs.cospectral_graphs(5)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Dr?', 'Ds_']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True

        There are two sets of cospectral graphs on six vertices with no isolated vertices::

            sage: g=graphs.cospectral_graphs(6, graphs=lambda x: min(x.degree())>0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Ep__', 'Er?G'], ['ExGg', 'ExoG']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True
            sage: g[1][1].am().charpoly()==g[1][1].am().charpoly()
            True

        There is one pair of cospectral trees on eight vertices::

            sage: g=graphs.cospectral_graphs(6, graphs=graphs.trees(8))
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['GiPC?C', 'GiQCC?']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True

        There are two sets of cospectral graphs (with respect to the
        Laplacian matrix) on six vertices::

            sage: g=graphs.cospectral_graphs(6, matrix_function=lambda g: g.laplacian_matrix())
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Edq_', 'ErcG'], ['Exoo', 'EzcG']]
            sage: g[0][1].laplacian_matrix().charpoly()==g[0][1].laplacian_matrix().charpoly()
            True
            sage: g[1][1].laplacian_matrix().charpoly()==g[1][1].laplacian_matrix().charpoly()
            True

        To find cospectral graphs with respect to the normalized
        Laplacian, assuming the graphs do not have an isolated vertex, it
        is enough to check the spectrum of the matrix `D^{-1}A`, where `D`
        is the diagonal matrix of vertex degrees, and A is the adjacency
        matrix.  We find two such cospectral graphs (for the normalized
        Laplacian) on five vertices::

            sage: def DinverseA(g):
            ...     A=g.adjacency_matrix().change_ring(QQ)
            ...     for i in range(g.order()):
            ...         A.rescale_row(i, 1/len(A.nonzero_positions_in_row(i)))
            ...     return A
            sage: g=graphs.cospectral_graphs(5, matrix_function=DinverseA, graphs=lambda g: min(g.degree())>0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Dlg', 'Ds_']]
            sage: g[0][1].laplacian_matrix(normalized=True).charpoly()==g[0][1].laplacian_matrix(normalized=True).charpoly()
            True
        """
        from sage.graphs.all import graphs as graph_gen
        if graphs is None:
            graph_list=graph_gen(vertices)
        elif callable(graphs):
            graph_list=iter(g for g in graph_gen(vertices) if graphs(g))
        else:
            graph_list=iter(graphs)

        from collections import defaultdict
        charpolys=defaultdict(list)
        for g in graph_list:
            cp=matrix_function(g).charpoly()
            charpolys[cp].append(g)

        cospectral_graphs=[]
        for cp,g_list in charpolys.items():
            if len(g_list)>1:
                cospectral_graphs.append(g_list)

        return cospectral_graphs

    def fullerenes(self, order, ipr=False):
        r"""
        Returns a generator which creates fullerene graphs using
        the buckygen generator (see [buckygen]_).

        INPUT:

        - ``order`` - a positive even integer smaller than or equal to 254.
          This specifies the number of vertices in the generated fullerenes.

        - ``ipr`` - default: ``False`` - if ``True`` only fullerenes that
          satisfy the Isolated Pentagon Rule are generated. This means that
          no pentagonal faces share an edge.

        OUTPUT:

        A generator which will produce the fullerene graphs as Sage graphs
        with an embedding set. These will be simple graphs: no loops, no
        multiple edges, no directed edges.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        There are 1812 isomers of `\textrm{C}_{60}`, i.e., 1812 fullerene graphs
        on 60 vertices:  ::

            sage: gen = graphs.fullerenes(60)  # optional buckygen
            sage: len(list(gen))  # optional buckygen
            1812

        However, there is only one IPR fullerene graph on 60 vertices: the famous
        Buckminster Fullerene:  ::

            sage: gen = graphs.fullerenes(60, ipr=True)  # optional buckygen
            sage: gen.next()  # optional buckygen
            Graph on 60 vertices
            sage: gen.next()  # optional buckygen
            Traceback (most recent call last):
            ...
            StopIteration

        The unique fullerene graph on 20 vertices is isomorphic to the dodecahedron
        graph. ::

            sage: gen = graphs.fullerenes(20)  # optional buckygen
            sage: g = gen.next()  # optional buckygen
            sage: g.is_isomorphic(graphs.DodecahedralGraph()) # optional buckygen
            True
            sage: g.get_embedding()  # optional buckygen
            {1: [2, 3, 4],
             2: [1, 5, 6],
             3: [1, 7, 8],
             4: [1, 9, 10],
             5: [2, 10, 11],
             6: [2, 12, 7],
             7: [3, 6, 13],
             8: [3, 14, 9],
             9: [4, 8, 15],
             10: [4, 16, 5],
             11: [5, 17, 12],
             12: [6, 11, 18],
             13: [7, 18, 14],
             14: [8, 13, 19],
             15: [9, 19, 16],
             16: [10, 15, 17],
             17: [11, 16, 20],
             18: [12, 20, 13],
             19: [14, 20, 15],
             20: [17, 19, 18]}
            sage: g.plot3d(layout='spring')  # optional buckygen

        REFERENCE:

        .. [buckygen] G. Brinkmann, J. Goedgebeur and B.D. McKay, Generation of Fullerenes,
          Journal of Chemical Information and Modeling, 52(11):2910-2918, 2012.
        """
        from sage.misc.package import is_package_installed
        if not is_package_installed("buckygen"):
            raise TypeError("the optional buckygen package is not installed")

        # number of vertices should be positive
        if order < 0:
            raise ValueError("Number of vertices should be positive.")

        # buckygen can only output fullerenes on up to 254 vertices
        if order > 254:
            raise ValueError("Number of vertices should be at most 254.")

        # fullerenes only exist for an even number of vertices, larger than 20
        # and different from 22
        if order % 2 == 1 or order < 20 or order == 22:
            return

        command = 'buckygen -'+('I' if ipr else '')+'d {0}d'.format(order)

        import subprocess
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)
        out = sp.stdout

        #start of code to read planar code

        header = out.read(15)
        assert header == '>>planar_code<<', 'Not a valid planar code header'

        #read graph per graph
        while True:
            c = out.read(1)
            if len(c)==0:
                return

            # Each graph is stored in the following way :
            #
            # The first character is the number of vertices, followed by
            # n11,...,n1k,null character,n21,...,n2k',null character, ...
            #
            # where the n1* are all neighbors of n1 and all n2* are the
            # neighbors of n2, ...
            #
            # Besides, these neighbors are enumerated in clockwise order.
            order = ord(c)

            zeroCount = 0

            g = [[] for i in range(order)]

            while zeroCount < order:
                c = out.read(1)
                if ord(c)==0:
                    zeroCount += 1
                else:
                    g[zeroCount].append(ord(c))

            #construct graph based on g
            g = {i+1:di for i,di in enumerate(g)}
            G = graph.Graph(g)
            G.set_embedding(g)
            yield(G)

    def fusenes(self, hexagon_count, benzenoids=False):
        r"""
        Returns a generator which creates fusenes and benzenoids using
        the benzene generator (see [benzene]_). Fusenes are planar
        polycyclic hydrocarbons with all bounded faces hexagons. Benzenoids
        are fusenes that are subgraphs of the hexagonal lattice.

        INPUT:

        - ``hexagon_count`` - a positive integer smaller than or equal to 30.
          This specifies the number of hexagons in the generated benzenoids.

        - ``benzenoids`` - default: ``False`` - if ``True`` only benzenoids are
          generated.

        OUTPUT:

        A generator which will produce the fusenes as Sage graphs
        with an embedding set. These will be simple graphs: no loops, no
        multiple edges, no directed edges.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        There is a unique fusene with 2 hexagons:  ::

            sage: gen = graphs.fusenes(2)  # optional benzene
            sage: len(list(gen))  # optional benzene
            1

        This fusene is naphtalene (`\textrm{C}_{10}\textrm{H}_{8}`).
        In the fusene graph the H-atoms are not stored, so this is
        a graph on just 10 vertices:  ::

            sage: gen = graphs.fusenes(2)  # optional benzene
            sage: gen.next()  # optional benzene
            Graph on 10 vertices
            sage: gen.next()  # optional benzene
            Traceback (most recent call last):
            ...
            StopIteration

        There are 6505 benzenoids with 9 hexagons:  ::

            sage: gen = graphs.fusenes(9, benzenoids=True)  # optional benzene
            sage: len(list(gen))  # optional benzene
            6505

        REFERENCE:

        .. [benzene] G. Brinkmann, G. Caporossi and P. Hansen, A Constructive Enumeration of Fusenes and Benzenoids,
          Journal of Algorithms, 45:155-166, 2002.
        """
        from sage.misc.package import is_package_installed
        if not is_package_installed("benzene"):
            raise TypeError("the optional benzene package is not installed")

        # number of hexagons should be positive
        if hexagon_count < 0:
            raise ValueError("Number of hexagons should be positive.")

        # benzene is only built for fusenes with up to 30 hexagons
        if hexagon_count > 30:
            raise ValueError("Number of hexagons should be at most 30.")

        # there are no fusenes with 0 hexagons
        if hexagon_count == 0:
            return

        # there is only one unique fusene with 1 hexagon (and benzene doesn't generate it)
        if hexagon_count == 1:
            g = {1:[6, 2], 2:[1, 3], 3:[2, 4], 4:[3, 5], 5:[4, 6], 6:[5, 1]}
            G = graph.Graph(g)
            G.set_embedding(g)
            yield(G)
            return

        command = 'benzene '+('b' if benzenoids else '')+' {0} p'.format(hexagon_count)

        import subprocess
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)
        out = sp.stdout

        #start of code to read planar code

        header = out.read(15)
        assert header == '>>planar_code<<', 'Not a valid planar code header'

        #read graph per graph
        while True:
            c = out.read(1)
            if len(c)==0:
                return

            # Each graph is stored in the following way :
            #
            # The first character is the number of vertices, followed by
            # n11,...,n1k,null character,n21,...,n2k',null character, ...
            #
            # where the n1* are all neighbors of n1 and all n2* are the
            # neighbors of n2, ...
            #
            # Besides, these neighbors are enumerated in clockwise order.
            order = ord(c)

            zeroCount = 0

            g = [[] for i in range(order)]

            while zeroCount < order:
                c = out.read(1)
                if ord(c)==0:
                    zeroCount += 1
                else:
                    g[zeroCount].append(ord(c))

            #construct graph based on g
            g = {i+1:di for i,di in enumerate(g)}
            G = graph.Graph(g)
            G.set_embedding(g)
            yield(G)

###########################################################################
# Basic Graphs
###########################################################################
    import sage.graphs.generators.basic
    BullGraph                = staticmethod(sage.graphs.generators.basic.BullGraph)
    ButterflyGraph           = staticmethod(sage.graphs.generators.basic.ButterflyGraph)
    CircularLadderGraph      = staticmethod(sage.graphs.generators.basic.CircularLadderGraph)
    ClawGraph                = staticmethod(sage.graphs.generators.basic.ClawGraph)
    CycleGraph               = staticmethod(sage.graphs.generators.basic.CycleGraph)
    CompleteGraph            = staticmethod(sage.graphs.generators.basic.CompleteGraph)
    CompleteBipartiteGraph   = staticmethod(sage.graphs.generators.basic.CompleteBipartiteGraph)
    CompleteMultipartiteGraph= staticmethod(sage.graphs.generators.basic.CompleteMultipartiteGraph)
    DiamondGraph             = staticmethod(sage.graphs.generators.basic.DiamondGraph)
    EmptyGraph               = staticmethod(sage.graphs.generators.basic.EmptyGraph)
    Grid2dGraph              = staticmethod(sage.graphs.generators.basic.Grid2dGraph)
    GridGraph                = staticmethod(sage.graphs.generators.basic.GridGraph)
    HouseGraph               = staticmethod(sage.graphs.generators.basic.HouseGraph)
    HouseXGraph              = staticmethod(sage.graphs.generators.basic.HouseXGraph)
    LadderGraph              = staticmethod(sage.graphs.generators.basic.LadderGraph)
    LollipopGraph            = staticmethod(sage.graphs.generators.basic.LollipopGraph)
    PathGraph                = staticmethod(sage.graphs.generators.basic.PathGraph)
    StarGraph                = staticmethod(sage.graphs.generators.basic.StarGraph)
    Toroidal6RegularGrid2dGraph = staticmethod(sage.graphs.generators.basic.Toroidal6RegularGrid2dGraph)
    ToroidalGrid2dGraph      = staticmethod(sage.graphs.generators.basic.ToroidalGrid2dGraph)

###########################################################################
# Small Graphs
###########################################################################
    import sage.graphs.generators.smallgraphs
    Balaban10Cage            = staticmethod(sage.graphs.generators.smallgraphs.Balaban10Cage)
    Balaban11Cage            = staticmethod(sage.graphs.generators.smallgraphs.Balaban11Cage)
    BidiakisCube             = staticmethod(sage.graphs.generators.smallgraphs.BidiakisCube)
    BiggsSmithGraph          = staticmethod(sage.graphs.generators.smallgraphs.BiggsSmithGraph)
    BlanusaFirstSnarkGraph   = staticmethod(sage.graphs.generators.smallgraphs.BlanusaFirstSnarkGraph)
    BlanusaSecondSnarkGraph  = staticmethod(sage.graphs.generators.smallgraphs.BlanusaSecondSnarkGraph)
    BrinkmannGraph           = staticmethod(sage.graphs.generators.smallgraphs.BrinkmannGraph)
    BrouwerHaemersGraph      = staticmethod(sage.graphs.generators.smallgraphs.BrouwerHaemersGraph)
    BuckyBall                = staticmethod(sage.graphs.generators.smallgraphs.BuckyBall)
    CameronGraph             = staticmethod(sage.graphs.generators.smallgraphs.CameronGraph)
    ChvatalGraph             = staticmethod(sage.graphs.generators.smallgraphs.ChvatalGraph)
    ClebschGraph             = staticmethod(sage.graphs.generators.smallgraphs.ClebschGraph)
    CoxeterGraph             = staticmethod(sage.graphs.generators.smallgraphs.CoxeterGraph)
    DesarguesGraph           = staticmethod(sage.graphs.generators.smallgraphs.DesarguesGraph)
    DoubleStarSnark          = staticmethod(sage.graphs.generators.smallgraphs.DoubleStarSnark)
    DurerGraph               = staticmethod(sage.graphs.generators.smallgraphs.DurerGraph)
    DyckGraph                = staticmethod(sage.graphs.generators.smallgraphs.DyckGraph)
    EllinghamHorton54Graph   = staticmethod(sage.graphs.generators.smallgraphs.EllinghamHorton54Graph)
    EllinghamHorton78Graph   = staticmethod(sage.graphs.generators.smallgraphs.EllinghamHorton78Graph)
    ErreraGraph              = staticmethod(sage.graphs.generators.smallgraphs.ErreraGraph)
    FlowerSnark              = staticmethod(sage.graphs.generators.smallgraphs.FlowerSnark)
    FolkmanGraph             = staticmethod(sage.graphs.generators.smallgraphs.FolkmanGraph)
    FosterGraph              = staticmethod(sage.graphs.generators.smallgraphs.FosterGraph)
    FranklinGraph            = staticmethod(sage.graphs.generators.smallgraphs.FranklinGraph)
    FruchtGraph              = staticmethod(sage.graphs.generators.smallgraphs.FruchtGraph)
    GoldnerHararyGraph       = staticmethod(sage.graphs.generators.smallgraphs.GoldnerHararyGraph)
    GrayGraph                = staticmethod(sage.graphs.generators.smallgraphs.GrayGraph)
    GrotzschGraph            = staticmethod(sage.graphs.generators.smallgraphs.GrotzschGraph)
    HallJankoGraph           = staticmethod(sage.graphs.generators.smallgraphs.HallJankoGraph)
    WellsGraph               = staticmethod(sage.graphs.generators.smallgraphs.WellsGraph)
    HarriesGraph             = staticmethod(sage.graphs.generators.smallgraphs.HarriesGraph)
    HarriesWongGraph         = staticmethod(sage.graphs.generators.smallgraphs.HarriesWongGraph)
    HeawoodGraph             = staticmethod(sage.graphs.generators.smallgraphs.HeawoodGraph)
    HerschelGraph            = staticmethod(sage.graphs.generators.smallgraphs.HerschelGraph)
    HigmanSimsGraph          = staticmethod(sage.graphs.generators.smallgraphs.HigmanSimsGraph)
    HoffmanGraph             = staticmethod(sage.graphs.generators.smallgraphs.HoffmanGraph)
    HoffmanSingletonGraph    = staticmethod(sage.graphs.generators.smallgraphs.HoffmanSingletonGraph)
    HoltGraph                = staticmethod(sage.graphs.generators.smallgraphs.HoltGraph)
    HortonGraph              = staticmethod(sage.graphs.generators.smallgraphs.HortonGraph)
    KittellGraph             = staticmethod(sage.graphs.generators.smallgraphs.KittellGraph)
    KrackhardtKiteGraph      = staticmethod(sage.graphs.generators.smallgraphs.KrackhardtKiteGraph)
    LjubljanaGraph           = staticmethod(sage.graphs.generators.smallgraphs.LjubljanaGraph)
    M22Graph                 = staticmethod(sage.graphs.generators.smallgraphs.M22Graph)
    MarkstroemGraph          = staticmethod(sage.graphs.generators.smallgraphs.MarkstroemGraph)
    McGeeGraph               = staticmethod(sage.graphs.generators.smallgraphs.McGeeGraph)
    McLaughlinGraph          = staticmethod(sage.graphs.generators.smallgraphs.McLaughlinGraph)
    MeredithGraph            = staticmethod(sage.graphs.generators.smallgraphs.MeredithGraph)
    MoebiusKantorGraph       = staticmethod(sage.graphs.generators.smallgraphs.MoebiusKantorGraph)
    MoserSpindle             = staticmethod(sage.graphs.generators.smallgraphs.MoserSpindle)
    NauruGraph               = staticmethod(sage.graphs.generators.smallgraphs.NauruGraph)
    PappusGraph              = staticmethod(sage.graphs.generators.smallgraphs.PappusGraph)
    PoussinGraph             = staticmethod(sage.graphs.generators.smallgraphs.PoussinGraph)
    PetersenGraph            = staticmethod(sage.graphs.generators.smallgraphs.PetersenGraph)
    RobertsonGraph           = staticmethod(sage.graphs.generators.smallgraphs.RobertsonGraph)
    SchlaefliGraph           = staticmethod(sage.graphs.generators.smallgraphs.SchlaefliGraph)
    ShrikhandeGraph          = staticmethod(sage.graphs.generators.smallgraphs.ShrikhandeGraph)
    SimsGewirtzGraph         = staticmethod(sage.graphs.generators.smallgraphs.SimsGewirtzGraph)
    SousselierGraph          = staticmethod(sage.graphs.generators.smallgraphs.SousselierGraph)
    SylvesterGraph           = staticmethod(sage.graphs.generators.smallgraphs.SylvesterGraph)
    SzekeresSnarkGraph       = staticmethod(sage.graphs.generators.smallgraphs.SzekeresSnarkGraph)
    ThomsenGraph             = staticmethod(sage.graphs.generators.smallgraphs.ThomsenGraph)
    TietzeGraph              = staticmethod(sage.graphs.generators.smallgraphs.TietzeGraph)
    Tutte12Cage              = staticmethod(sage.graphs.generators.smallgraphs.Tutte12Cage)
    TutteCoxeterGraph        = staticmethod(sage.graphs.generators.smallgraphs.TutteCoxeterGraph)
    TutteGraph               = staticmethod(sage.graphs.generators.smallgraphs.TutteGraph)
    WagnerGraph              = staticmethod(sage.graphs.generators.smallgraphs.WagnerGraph)
    WatkinsSnarkGraph        = staticmethod(sage.graphs.generators.smallgraphs.WatkinsSnarkGraph)
    WienerArayaGraph         = staticmethod(sage.graphs.generators.smallgraphs.WienerArayaGraph)

###########################################################################
# Platonic Solids
###########################################################################
    import sage.graphs.generators.platonic_solids
    DodecahedralGraph        = staticmethod(sage.graphs.generators.platonic_solids.DodecahedralGraph)
    HexahedralGraph          = staticmethod(sage.graphs.generators.platonic_solids.HexahedralGraph)
    IcosahedralGraph         = staticmethod(sage.graphs.generators.platonic_solids.IcosahedralGraph)
    OctahedralGraph          = staticmethod(sage.graphs.generators.platonic_solids.OctahedralGraph)
    TetrahedralGraph         = staticmethod(sage.graphs.generators.platonic_solids.TetrahedralGraph)

###########################################################################
# Families
###########################################################################
    import sage.graphs.generators.families
    BalancedTree           = staticmethod(sage.graphs.generators.families.BalancedTree)
    BarbellGraph           = staticmethod(sage.graphs.generators.families.BarbellGraph)
    BubbleSortGraph        = staticmethod(sage.graphs.generators.families.BubbleSortGraph)
    CirculantGraph         = staticmethod(sage.graphs.generators.families.CirculantGraph)
    CubeGraph              = staticmethod(sage.graphs.generators.families.CubeGraph)
    DorogovtsevGoltsevMendesGraph = staticmethod(sage.graphs.generators.families.DorogovtsevGoltsevMendesGraph)
    FibonacciTree          = staticmethod(sage.graphs.generators.families.FibonacciTree)
    FoldedCubeGraph        = staticmethod(sage.graphs.generators.families.FoldedCubeGraph)
    FriendshipGraph        = staticmethod(sage.graphs.generators.families.FriendshipGraph)
    FuzzyBallGraph         = staticmethod(sage.graphs.generators.families.FuzzyBallGraph)
    GeneralizedPetersenGraph = staticmethod(sage.graphs.generators.families.GeneralizedPetersenGraph)
    HanoiTowerGraph        = staticmethod(sage.graphs.generators.families.HanoiTowerGraph)
    HararyGraph            = staticmethod(sage.graphs.generators.families.HararyGraph)
    HyperStarGraph         = staticmethod(sage.graphs.generators.families.HyperStarGraph)
    JohnsonGraph           = staticmethod(sage.graphs.generators.families.JohnsonGraph)
    KneserGraph            = staticmethod(sage.graphs.generators.families.KneserGraph)
    LCFGraph               = staticmethod(sage.graphs.generators.families.LCFGraph)
    line_graph_forbidden_subgraphs = staticmethod(sage.graphs.generators.families.line_graph_forbidden_subgraphs)
    MycielskiGraph         = staticmethod(sage.graphs.generators.families.MycielskiGraph)
    MycielskiStep          = staticmethod(sage.graphs.generators.families.MycielskiStep)
    NKStarGraph            = staticmethod(sage.graphs.generators.families.NKStarGraph)
    NStarGraph             = staticmethod(sage.graphs.generators.families.NStarGraph)
    OddGraph               = staticmethod(sage.graphs.generators.families.OddGraph)
    PaleyGraph             = staticmethod(sage.graphs.generators.families.PaleyGraph)
    RingedTree             = staticmethod(sage.graphs.generators.families.RingedTree)
    SymplecticGraph        = staticmethod(sage.graphs.generators.families.SymplecticGraph)
    trees                  = staticmethod(sage.graphs.generators.families.trees)
    WheelGraph             = staticmethod(sage.graphs.generators.families.WheelGraph)

###########################################################################
# Chessboard Graphs
###########################################################################
    import sage.graphs.generators.chessboard
    ChessboardGraphGenerator = staticmethod(sage.graphs.generators.chessboard.ChessboardGraphGenerator)
    BishopGraph              = staticmethod(sage.graphs.generators.chessboard.BishopGraph)
    KingGraph                = staticmethod(sage.graphs.generators.chessboard.KingGraph)
    KnightGraph              = staticmethod(sage.graphs.generators.chessboard.KnightGraph)
    QueenGraph               = staticmethod(sage.graphs.generators.chessboard.QueenGraph)
    RookGraph                = staticmethod(sage.graphs.generators.chessboard.RookGraph)

###########################################################################
# Intersection graphs
###########################################################################
    import sage.graphs.generators.intersection
    IntervalGraph            = staticmethod(sage.graphs.generators.intersection.IntervalGraph)
    PermutationGraph         = staticmethod(sage.graphs.generators.intersection.PermutationGraph)
    ToleranceGraph           = staticmethod(sage.graphs.generators.intersection.ToleranceGraph)

###########################################################################
# Random Graphs
###########################################################################
    import sage.graphs.generators.random
    RandomBarabasiAlbert     = staticmethod(sage.graphs.generators.random.RandomBarabasiAlbert)
    RandomBipartite          = staticmethod(sage.graphs.generators.random.RandomBipartite)
    RandomBoundedToleranceGraph = staticmethod(sage.graphs.generators.random.RandomBoundedToleranceGraph)
    RandomGNM                = staticmethod(sage.graphs.generators.random.RandomGNM)
    RandomGNP                = staticmethod(sage.graphs.generators.random.RandomGNP)
    RandomHolmeKim           = staticmethod(sage.graphs.generators.random.RandomHolmeKim)
    RandomInterval           = staticmethod(sage.graphs.generators.random.RandomInterval) # deprecated
    RandomIntervalGraph      = staticmethod(sage.graphs.generators.random.RandomIntervalGraph)
    RandomLobster            = staticmethod(sage.graphs.generators.random.RandomLobster)
    RandomNewmanWattsStrogatz = staticmethod(sage.graphs.generators.random.RandomNewmanWattsStrogatz)
    RandomRegular            = staticmethod(sage.graphs.generators.random.RandomRegular)
    RandomShell              = staticmethod(sage.graphs.generators.random.RandomShell)
    RandomToleranceGraph     = staticmethod(sage.graphs.generators.random.RandomToleranceGraph)
    RandomTreePowerlaw       = staticmethod(sage.graphs.generators.random.RandomTreePowerlaw)
    RandomTree               = staticmethod(sage.graphs.generators.random.RandomTree)

###########################################################################
# World Map
###########################################################################
    import sage.graphs.generators.world_map
    WorldMap = staticmethod(sage.graphs.generators.world_map.WorldMap)

###########################################################################
# Degree Sequence
###########################################################################
    import sage.graphs.generators.degree_sequence
    DegreeSequence           = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequence)
    DegreeSequenceBipartite  = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceBipartite)
    DegreeSequenceConfigurationModel = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceConfigurationModel)
    DegreeSequenceTree       = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceTree)
    DegreeSequenceExpected   = staticmethod(sage.graphs.generators.degree_sequence.DegreeSequenceExpected)

def canaug_traverse_vert(g, aut_gens, max_verts, property, dig=False, loops=False, implementation='c_graph', sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free (di)graphs satisfying a
    given property.

    INPUT:


    -  ``g`` - current position on the tree.

    -  ``aut_gens`` - list of generators of Aut(g), in
       list notation.

    -  ``max_verts`` - when to retreat.

    -  ``property`` - check before traversing below g.

    -  ``degree_sequence`` - specify a degree sequence to try to
       obtain.


    EXAMPLES::

        sage: from sage.graphs.graph_generators import canaug_traverse_vert
        sage: list(canaug_traverse_vert(Graph(), [], 3, lambda x: True))
        [Graph on 0 vertices, ... Graph on 3 vertices]

    The best way to access this function is through the graphs()
    iterator:

    Print graphs on 3 or less vertices.

    ::

        sage: for G in graphs(3, augment='vertices'):
        ...    print G
        ...
        Graph on 0 vertices
        Graph on 1 vertex
        Graph on 2 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 2 vertices
        Graph on 3 vertices

    Print digraphs on 2 or less vertices.

    ::

        sage: for D in digraphs(2, augment='vertices'):
        ...    print D
        ...
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices
    """
    from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree

    if not property(g):
        return
    yield g

    n = g.order()
    if n < max_verts:

        # build a list representing C(g) - the vertex to be added
        # is at the end, so only specify which edges...
        # in the case of graphs, there are n possibilities,
        # and in the case of digraphs, there are 2*n.
        if dig:
            possibilities = 2*n
        else:
            possibilities = n
        num_roots = 2**possibilities
        children = [-1]*num_roots

        # union-find C(g) under Aut(g)
        for gen in aut_gens:
            for i in xrange(len(children)):
                k = 0
                for j in xrange(possibilities):
                    if (1 << j)&i:
                        if dig and j >= n:
                            k += (1 << (gen[j-n]+n))
                        else:
                            k += (1 << gen[j])
                while children[k] != -1:
                    k = children[k]
                while children[i] != -1:
                    i = children[i]
                if i != k:
                    # union i & k
                    smaller, larger = sorted([i,k])
                    children[larger] = smaller
                    num_roots -= 1

        # find representatives of orbits of C(g)
        roots = []
        found_roots = 0
        i = 0
        while found_roots < num_roots:
            if children[i] == -1:
                found_roots += 1
                roots.append(i)
            i += 1
        for i in roots:
            # construct a z for each number in roots...
            z = g.copy(implementation=implementation, sparse=sparse)
            z.add_vertex(n)
            edges = []
            if dig:
                index = 0
                while index < possibilities/2:
                    if (1 << index)&i:
                        edges.append((index,n))
                    index += 1
                while index < possibilities:
                    if (1 << index)&i:
                        edges.append((n,index-n))
                    index += 1
            else:
                index = 0
                while (1 << index) <= i:
                    if (1 << index)&i:
                        edges.append((index,n))
                    index += 1
            z.add_edges(edges)
            z_s = []
            if property(z):
                z_s.append(z)
            if loops:
                z = z.copy(implementation=implementation, sparse=sparse)
                z.add_edge((n,n))
                if property(z):
                    z_s.append(z)
            for z in z_s:
                z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices()], certify=True, dig=(dig or loops))
                cut_vert = 0
                while canonical_relabeling[cut_vert] != n:
                    cut_vert += 1
                sub_verts = [v for v in z if v != cut_vert]
                m_z = z.subgraph(sub_verts)

                if m_z == g:
                    for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                        yield a
                else:
                    for possibility in check_aut(z_aut_gens, cut_vert, n):
                        if m_z.relabel(possibility, inplace=False) == g:
                            for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                                yield a
                            break

def check_aut(aut_gens, cut_vert, n):
    """
    Helper function for exhaustive generation.

    At the start, check_aut is given a set of generators for the
    automorphism group, aut_gens. We already know we are looking for
    an element of the auto- morphism group that sends cut_vert to n,
    and check_aut generates these for the canaug_traverse function.

    EXAMPLE: Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3.

    ::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = range(n+1)
    seen_perms = [perm]
    unchecked_perms = [perm]
    while len(unchecked_perms) != 0:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for i in xrange(len(perm)):
                new_perm[i] = gen[perm[i]]
            if new_perm not in seen_perms:
                seen_perms.append(new_perm)
                unchecked_perms.append(new_perm)
                if new_perm[cut_vert] == n:
                    yield new_perm

def canaug_traverse_edge(g, aut_gens, property, dig=False, loops=False, implementation='c_graph', sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free graphs satisfying a
    given property.

    INPUT:


    -  ``g`` - current position on the tree.

    -  ``aut_gens`` - list of generators of Aut(g), in
       list notation.

    -  ``property`` - check before traversing below g.


    EXAMPLES::

        sage: from sage.graphs.graph_generators import canaug_traverse_edge
        sage: G = Graph(3)
        sage: list(canaug_traverse_edge(G, [], lambda x: True))
        [Graph on 3 vertices, ... Graph on 3 vertices]

    The best way to access this function is through the graphs()
    iterator:

    Print graphs on 3 or less vertices.

    ::

        sage: for G in graphs(3):
        ...    print G
        ...
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices

    Print digraphs on 3 or less vertices.

    ::

        sage: for G in digraphs(3):
        ...    print G
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices
    """
    from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree

    if not property(g):
        return
    yield g
    n = g.order()
    if dig:
        max_size = n*(n-1)
    else:
        max_size = (n*(n-1))>>1 # >> 1 is just / 2 (this is n choose 2)
    if loops: max_size += n
    if g.size() < max_size:
        # build a list representing C(g) - the edge to be added
        # is one of max_size choices
        if dig:
            children = [[(j,i) for i in xrange(n)] for j in xrange(n)]
        else:
            children = [[(j,i) for i in xrange(j)] for j in xrange(n)]
        # union-find C(g) under Aut(g)
        orbits = range(n)
        for gen in aut_gens:
            for iii in xrange(n):
                if orbits[gen[iii]] != orbits[iii]:
                    temp = orbits[gen[iii]]
                    for jjj in xrange(n):
                        if orbits[jjj] == temp:
                            orbits[jjj] = orbits[iii]
                if dig:
                    jjj_range = range(iii) + range(iii+1, n)
                else:
                    jjj_range = xrange(iii) # iii > jjj
                for jjj in jjj_range:
                    i, j = iii, jjj
                    if dig:
                        x, y = gen[i], gen[j]
                    else:
                        y, x = sorted([gen[i], gen[j]])
                    if children[i][j] != children[x][y]:
                        x_val, y_val = x, y
                        i_val, j_val = i, j
                        if dig:
                            while (x_val, y_val) != children[x_val][y_val]:
                                x_val, y_val = children[x_val][y_val]
                            while (i_val, j_val) != children[i_val][j_val]:
                                i_val, j_val = children[i_val][j_val]
                        else:
                            while (x_val, y_val) != children[x_val][y_val]:
                                y_val, x_val = sorted(children[x_val][y_val])
                            while (i_val, j_val) != children[i_val][j_val]:
                                j_val, i_val = sorted(children[i_val][j_val])
                        while (x, y) != (x_val, y_val):
                            xx, yy = x, y
                            x, y = children[x][y]
                            children[xx][yy] = (x_val, y_val)
                        while (i, j) != (i_val, j_val):
                            ii, jj = i, j
                            i, j = children[i][j]
                            children[ii][jj] = (i_val, j_val)
                        if x < i:
                            children[i][j] = (x, y)
                        elif x > i:
                            children[x][y] = (i, j)
                        elif y < j:
                            children[i][j] = (x, y)
                        elif y > j:
                            children[x][y] = (i, j)
                        else:
                            continue
        # find representatives of orbits of C(g)
        roots = []
        for i in range(n):
            if dig:
                j_range = range(i) + range(i+1, n)
            else:
                j_range = range(i)
            for j in j_range:
                if children[i][j] == (i, j):
                    roots.append((i,j))
        if loops:
            seen = []
            for i in xrange(n):
                if orbits[i] not in seen:
                    roots.append((i,i))
                    seen.append(orbits[i])
        for i, j in roots:
            if g.has_edge(i, j):
                continue
            # construct a z for each edge in roots...
            z = g.copy(implementation=implementation, sparse=sparse)
            z.add_edge(i, j)
            if not property(z):
                continue
            z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices()], certify=True, dig=(dig or loops))
            relabel_inverse = [0]*n
            for ii in xrange(n):
                relabel_inverse[canonical_relabeling[ii]] = ii
            z_can = z.relabel(canonical_relabeling, inplace=False)
            cut_edge_can = z_can.edges(labels=False, sort=True)[-1]
            cut_edge = [relabel_inverse[cut_edge_can[0]], relabel_inverse[cut_edge_can[1]]]
            if dig:
                cut_edge = tuple(cut_edge)
            else:
                cut_edge = tuple(sorted(cut_edge))

            from copy import copy
            m_z = copy(z)
            m_z.delete_edge(cut_edge)
            if m_z == g:
                for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                    yield a
            else:
                for possibility in check_aut_edge(z_aut_gens, cut_edge, i, j, n, dig=dig):
                    if m_z.relabel(possibility, inplace=False) == g:
                        for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, implementation=implementation, sparse=sparse):
                            yield a
                        break

def check_aut_edge(aut_gens, cut_edge, i, j, n, dig=False):
    """
    Helper function for exhaustive generation.

    At the start, check_aut_edge is given a set of generators for the
    automorphism group, aut_gens. We already know we are looking for
    an element of the auto- morphism group that sends cut_edge to {i,
    j}, and check_aut generates these for the canaug_traverse
    function.

    EXAMPLE: Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3.

    ::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = range(n)
    seen_perms = [perm]
    unchecked_perms = [perm]
    while len(unchecked_perms) != 0:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for ii in xrange(n):
                new_perm[ii] = gen[perm[ii]]
            if new_perm not in seen_perms:
                seen_perms.append(new_perm)
                unchecked_perms.append(new_perm)
                if new_perm[cut_edge[0]] == i and new_perm[cut_edge[1]] == j:
                    yield new_perm
                if not dig and new_perm[cut_edge[0]] == j and new_perm[cut_edge[1]] == i:
                    yield new_perm


# Easy access to the graph generators from the command line:
graphs = GraphGenerators()
