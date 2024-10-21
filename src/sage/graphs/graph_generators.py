# sage_setup: distribution = sagemath-graphs
r"""
Common graphs

All graphs in Sage can be built through the ``graphs`` object. In order to
build a complete graph on 15 elements, one can do::

    sage: g = graphs.CompleteGraph(15)

To get a path with 4 vertices, and the house graph::

    sage: p = graphs.PathGraph(4)
    sage: h = graphs.HouseGraph()

More interestingly, one can get the list of all graphs that Sage knows how to
build by typing ``graphs.`` in Sage and then hitting :kbd:`Tab`.
"""

import subprocess


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
        reordered_methods[3*(i % h) + (i//h)] = m
    methods = reordered_methods

    # Adding the list to the __doc__ string
    def wrap_name(x):
        if x:
            return ":meth:`" + str(x) + " <GraphGenerators." + str(x) + ">`"
        return ""

    while methods:
        a = methods.pop(0)
        b = methods.pop(0)
        c = methods.pop(0)
        __doc__ += "    " + wrap_name(a) + " | " + wrap_name(b) + " | " + wrap_name(c) + "\n"


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
     "CorrelationGraph",
     "DiamondGraph",
     "GemGraph",
     "DartGraph",
     "ForkGraph",
     "DipoleGraph",
     "EmptyGraph",
     "Grid2dGraph",
     "GridGraph",
     "HouseGraph",
     "HouseXGraph",
     "LadderGraph",
     "LollipopGraph",
     "MoebiusLadderGraph",
     "PathGraph",
     "StarGraph",
     "TadpoleGraph",
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
     "Cell600",
     "Cell120",
     "ChvatalGraph",
     "ClebschGraph",
     "cocliques_HoffmannSingleton",
     "ConwaySmith_for_3S7",
     "CoxeterGraph",
     "CubeplexGraph",
     "DesarguesGraph",
     "DejterGraph",
     "distance_3_doubly_truncated_Golay_code_graph",
     "DoubleStarSnark",
     "DoublyTruncatedWittGraph",
     "DurerGraph",
     "DyckGraph",
     "EllinghamHorton54Graph",
     "EllinghamHorton78Graph",
     "ErreraGraph",
     "F26AGraph",
     "FlowerSnark",
     "FolkmanGraph",
     "FosterGraph",
     "FosterGraph3S6",
     "FranklinGraph",
     "FruchtGraph",
     "GoldnerHararyGraph",
     "GolombGraph",
     "GossetGraph",
     "graph_3O73",
     "GrayGraph",
     "GritsenkoGraph",
     "GrotzschGraph",
     "HallJankoGraph",
     "HarborthGraph",
     "HarriesGraph",
     "HarriesWongGraph",
     "HeawoodGraph",
     "HerschelGraph",
     "HigmanSimsGraph",
     "HoffmanGraph",
     "HoffmanSingletonGraph",
     "HoltGraph",
     "HortonGraph",
     "IoninKharaghani765Graph",
     "IvanovIvanovFaradjevGraph",
     "J2Graph",
     "JankoKharaghaniGraph",
     "JankoKharaghaniTonchevGraph",
     "KittellGraph",
     "KrackhardtKiteGraph",
     "Klein3RegularGraph",
     "Klein7RegularGraph",
     "LargeWittGraph",
     "LeonardGraph",
     "LjubljanaGraph",
     "vanLintSchrijverGraph",
     "LivingstoneGraph",
     "locally_GQ42_distance_transitive_graph",
     "LocalMcLaughlinGraph",
     "M22Graph",
     "MarkstroemGraph",
     "MathonStronglyRegularGraph",
     "McGeeGraph",
     "McLaughlinGraph",
     "MeredithGraph",
     "MoebiusKantorGraph",
     "MoserSpindle",
     "MurtyGraph",
     "NauruGraph",
     "PappusGraph",
     "PoussinGraph",
     "PerkelGraph",
     "PetersenGraph",
     "RobertsonGraph",
     "SchlaefliGraph",
     "shortened_00_11_binary_Golay_code_graph",
     "shortened_000_111_extended_binary_Golay_code_graph",
     "ShrikhandeGraph",
     "SimsGewirtzGraph",
     "SousselierGraph",
     "SylvesterGraph",
     "SzekeresSnarkGraph",
     "ThomsenGraph",
     "TietzeGraph",
     "TricornGraph",
     "TruncatedIcosidodecahedralGraph",
     "TruncatedTetrahedralGraph",
     "TruncatedWittGraph",
     "Tutte12Cage",
     "TutteCoxeterGraph",
     "TutteGraph",
     "TwinplexGraph",
     "U42Graph216",
     "U42Graph540",
     "WagnerGraph",
     "WatkinsSnarkGraph",
     "WellsGraph",
     "WienerArayaGraph",
     "SuzukiGraph"])

__doc__ += """
**Platonic solids** (ordered ascending by number of vertices)
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
    ["AlternatingFormsGraph",
     "AztecDiamondGraph",
     "BalancedTree",
     "BarbellGraph",
     "BilinearFormsGraph",
     "BiwheelGraph",
     "BubbleSortGraph",
     "CaiFurerImmermanGraph",
     "chang_graphs",
     "CirculantGraph",
     "cographs",
     "cospectral_graphs",
     "CubeGraph",
     "CubeConnectedCycle",
     "distance_regular_graph",
     "DorogovtsevGoltsevMendesGraph",
     "DoubleGrassmannGraph",
     "DoubleOddGraph",
     "EgawaGraph",
     "FibonacciTree",
     "FoldedCubeGraph",
     "FriendshipGraph",
     "fullerenes",
     "FurerGadget",
     "fusenes",
     "FuzzyBallGraph",
     "GeneralisedDodecagonGraph",
     "GeneralisedHexagonGraph",
     "GeneralisedOctagonGraph",
     "GeneralizedPetersenGraph",
     "GeneralizedSierpinskiGraph",
     "GoethalsSeidelGraph",
     "GrassmannGraph",
     "HalfCube",
     "HammingGraph",
     "HanoiTowerGraph",
     "HararyGraph",
     "HermitianFormsGraph",
     "HyperStarGraph",
     "JohnsonGraph",
     "KneserGraph",
     "LCFGraph",
     "line_graph_forbidden_subgraphs",
     "MathonPseudocyclicMergingGraph",
     "MathonPseudocyclicStronglyRegularGraph",
     "MuzychukS6Graph",
     "MycielskiGraph",
     "MycielskiStep",
     "nauty_geng",
     "nauty_genbg",
     "NKStarGraph",
     "NStarGraph",
     "OddGraph",
     "PaleyGraph",
     "PasechnikGraph",
     "petersen_family",
     "planar_graphs",
     "plantri_gen",
     "quadrangulations",
     "RingedTree",
     "SierpinskiGasketGraph",
     "SquaredSkewHadamardMatrixGraph",
     "SwitchedSquaredSkewHadamardMatrixGraph",
     "StaircaseGraph",
     "strongly_regular_graph",
     "trees",
     "TruncatedBiwheelGraph",
     "nauty_gentreeg",
     "triangulations",
     "TuranGraph",
     "UstimenkoGraph",
     "WheelGraph",
     "WindmillGraph"])


__doc__ += """
**Graphs from classical geometries over finite fields**

A number of classes of graphs related to geometries over finite fields and
quadrics and Hermitean varieties there.
"""

__append_to_doc(
    ["AffineOrthogonalPolarGraph",
     "AhrensSzekeresGeneralizedQuadrangleGraph",
     "NonisotropicOrthogonalPolarGraph",
     "NonisotropicUnitaryPolarGraph",
     "OrthogonalDualPolarGraph",
     "OrthogonalPolarGraph",
     "SymplecticDualPolarGraph",
     "SymplecticPolarGraph",
     "TaylorTwographDescendantSRG",
     "TaylorTwographSRG",
     "T2starGeneralizedQuadrangleGraph",
     "Nowhere0WordsTwoWeightCodeGraph",
     "HaemersGraph",
     "CossidentePenttilaGraph",
     "UnitaryDualPolarGraph",
     "UnitaryPolarGraph"])

__doc__ += """
**Chessboard Graphs**
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
    ["IntersectionGraph",
     "IntervalGraph",
     "OrthogonalArrayBlockGraph",
     "PermutationGraph",
     "ToleranceGraph"])

__doc__ += """
**Random graphs**
"""

__append_to_doc(
    ["RandomBarabasiAlbert",
     "RandomBicubicPlanar",
     "RandomBipartite",
     "RandomRegularBipartite",
     "RandomBlockGraph",
     "RandomBoundedToleranceGraph",
     "RandomGNM",
     "RandomGNP",
     "RandomHolmeKim",
     "RandomChordalGraph",
     "RandomIntervalGraph",
     "RandomKTree",
     "RandomPartialKTree",
     "RandomLobster",
     "RandomNewmanWattsStrogatz",
     "RandomProperIntervalGraph",
     "RandomRegular",
     "RandomShell",
     "RandomToleranceGraph",
     "RandomTree",
     "RandomTreePowerlaw",
     "RandomTriangulation",
     "RandomUnitDiskGraph"])

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
    ["WorldMap",
     "EuropeMap",
     "AfricaMap",
     "USAMap"]
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

- Ed Scheinerman (2010-11-15): added RandomTriangulation

- Minh Van Nguyen (2010-11-26): added more named graphs

- Keshav Kini (2011-02-16): added Shrikhande and Dyck graphs

- David Coudert (2012-02-10): new RandomGNP generator

- David Coudert (2012-08-02): added chessboard graphs: Queen, King,
  Knight, Bishop, and Rook graphs

- Nico Van Cleemput (2013-05-26): added fullerenes

- Nico Van Cleemput (2013-07-01): added benzenoids

- Birk Eisermann (2013-07-29): new section 'intersection graphs',
  added (random, bounded) tolerance graphs

- Marco Cognetta (2016-03-03): added TuranGraph

- Janmenjaya Panda (2024-05-26): added MoebiusLadderGraph

- Janmenjaya Panda (2024-06-09): added StaircaseGraph, BiwheelGraph and
  TruncatedBiwheelGraph


Functions and methods
---------------------
"""

# ****************************************************************************
#       Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                          Emily A. Kirkman
#                     2009 Michael C. Yurko <myurko@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# import from Python standard library

# import from Sage library
from . import graph


class GraphGenerators:
    r"""
    A class consisting of constructors for several common graphs, as well as
    orderly generation of isomorphism class representatives. See the
    :mod:`module's help <sage.graphs.graph_generators>` for a list of supported
    constructors.

    A list of all graphs and graph structures (other than isomorphism class
    representatives) in this database is available via tab completion. Type
    "graphs." and then hit the :kbd:`Tab` key to see which graphs are available.

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

    Also: see the use of the nauty package for generating graphs
    at the :meth:`nauty_geng` method.

    INPUT:

    - ``vertices`` -- a natural number or ``None`` to infinitely generate
      bigger and bigger graphs

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
        adding one edge. In this case, all graphs on *exactly* ``n=vertices`` are
        generated. If for any graph G satisfying the property, every
        subgraph, obtained from G by deleting one edge but not the vertices
        incident to that edge, satisfies the property, then this will
        generate all graphs with that property. If this does not hold, then
        all the graphs generated will satisfy the property, but there will
        be some missing.

      - ``'vertices'`` -- augments by adding a vertex and
        edges incident to that vertex. In this case, all graphs *up to*
        ``n=vertices`` are generated. If for any graph G satisfying the
        property, every subgraph, obtained from G by deleting one vertex
        and only edges incident to that vertex, satisfies the property,
        then this will generate all graphs with that property. If this does
        not hold, then all the graphs generated will satisfy the property,
        but there will be some missing.

    - ``size`` -- (default: ``None``) the size of the graph to be generated

    - ``degree_sequence`` -- (default: ``None``) a sequence of nonnegative integers,
      or ``None``. If specified, the generated graphs will have these
      integers for degrees. In this case, property and size are both
      ignored.

    - ``loops`` -- boolean (default: ``False``); whether to allow loops in the graph
      or not

    - ``sparse`` -- (default: ``True``) whether to use a sparse or dense data
      structure. See the documentation of :class:`~sage.graphs.graph.Graph`.

    - ``copy`` -- boolean; if set to ``True`` (default)
      this method makes copies of the graphs before returning
      them. If set to ``False`` the method returns the graph it
      is working on. The second alternative is faster, but modifying
      any of the graph instances returned by the method may break
      the function's behaviour, as it is using these graphs to
      compute the next ones: only use ``copy = False`` when
      you stick to *reading* the graphs returned.

    EXAMPLES:

    Print graphs on 3 or less vertices::

        sage: for G in graphs(3, augment='vertices'):
        ....:     print(G)
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
        ....:    print(G)
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
        sage: graphs_list.show_graphs(L)        # long time                             # needs sage.plot

    Generate all graphs with up to 5 vertices and up to 4 edges.

    ::

        sage: L = list(graphs(5, lambda G: G.size() <= 4, augment='vertices'))
        sage: len(L)
        31
        sage: graphs_list.show_graphs(L)        # long time                             # needs sage.plot

    Generate all graphs with degree at most 2, up to 6 vertices.

    ::

        sage: property = lambda G: ( max([G.degree(v) for v in G] + [0]) <= 2 )
        sage: L = list(graphs(6, property, augment='vertices'))
        sage: len(L)
        45

    Generate all bipartite graphs on up to 7 vertices: (see
    :oeis:`A033995`)

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
        sage: len(list(graphs(4, property)))                                            # needs sage.groups
        1
        sage: sum(1 for g in graphs(4) if property(g))                                  # needs sage.groups
        4

        sage: property = lambda G: G.is_bipartite()
        sage: len(list(graphs(4, property)))
        7
        sage: sum(1 for g in graphs(4) if property(g))
        7

    Generate graphs on the fly: (see :oeis:`A000088`)

    ::

        sage: for i in range(7):
        ....:     print(len(list(graphs(i))))
        1
        1
        2
        4
        11
        34
        156

    Generate all simple graphs, allowing loops: (see :oeis:`A000666`)

    ::

        sage: L = list(graphs(5,augment='vertices',loops=True))               # long time
        sage: for i in [0..5]:  # long time
        ....:     print((i, len([g for g in L if g.order() == i])))
        (0, 1)
        (1, 2)
        (2, 6)
        (3, 20)
        (4, 90)
        (5, 544)

    Generate all graphs with a specified degree sequence (see :oeis:`A002851`)::

        sage: for i in [4,6,8]:  # long time (4s on sage.math, 2012)
        ....:     print((i, len([g for g in graphs(i, degree_sequence=[3]*i) if g.is_connected()])))
        (4, 1)
        (6, 2)
        (8, 5)
        sage: for i in [4,6,8]:  # long time (7s on sage.math, 2012)
        ....:     print((i, len([g for g in graphs(i, augment='vertices', degree_sequence=[3]*i) if g.is_connected()])))
        (4, 1)
        (6, 2)
        (8, 5)

    ::

        sage: print((10, len([g for g in graphs(10,degree_sequence=[3]*10) if g.is_connected()]))) # not tested
        (10, 19)

    Make sure that the graphs are really independent and the generator
    survives repeated vertex removal (:issue:`8458`)::

        sage: for G in graphs(3):
        ....:     G.delete_vertex(0)
        ....:     print(G.order())
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

    def __call__(self, vertices=None, property=None, augment='edges', size=None,
                 degree_sequence=None, loops=False, sparse=True, copy=True):
        """
        Access the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives. See the docstring
        of this class for full documentation.

        EXAMPLES:

        Print graphs on 3 or less vertices::

            sage: for G in graphs(3, augment='vertices'):
            ....:    print(G)
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
            ....:    if g.num_verts() > 3: break
            ....:    print(g)
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
        # Use nauty for the basic case, as it is much faster.
        if (vertices and property is None and size is None and
                degree_sequence is None and not loops and augment == 'edges' and
                sparse and copy):
            for g in graphs.nauty_geng(vertices):
                yield g
            return

        if property is None:
            def property(x):
                return True

        from sage.graphs.graph import Graph
        from copy import copy as copyfun

        if degree_sequence is not None:
            if vertices is None:
                raise NotImplementedError
            if (len(degree_sequence) != vertices or sum(degree_sequence) % 2
                    or sum(degree_sequence) > vertices*(vertices - 1)):
                raise ValueError("Invalid degree sequence.")
            degree_sequence = sorted(degree_sequence)
            if augment == 'edges':
                def property(x):
                    D = sorted(x.degree())
                    return all(degree_sequence[i] >= d for i, d in enumerate(D))

                def extra_property(x):
                    return degree_sequence == sorted(x.degree())
            else:
                def property(x):
                    D = sorted(x.degree() + [0] * (vertices - x.num_verts()))
                    return all(degree_sequence[i] >= d for i, d in enumerate(D))

                def extra_property(x):
                    if x.num_verts() != vertices:
                        return False
                    return degree_sequence == sorted(x.degree())
        elif size is not None:
            def extra_property(x):
                return x.size() == size
        else:
            def extra_property(x):
                return True

        if augment == 'vertices':
            if vertices is None:
                raise NotImplementedError
            g = Graph(loops=loops, sparse=sparse)
            for gg in canaug_traverse_vert(g, [], vertices, property, loops=loops, sparse=sparse):
                if extra_property(gg):
                    yield copyfun(gg) if copy else gg
        elif augment == 'edges':
            if vertices is None:
                from sage.rings.integer import Integer
                vertices = Integer(0)
                while True:
                    for g in self(vertices, loops=loops, sparse=sparse):
                        yield copyfun(g) if copy else g
                    vertices += 1
            g = Graph(vertices, loops=loops, sparse=sparse)
            gens = []
            for i in range(vertices - 1):
                gen = list(range(i))
                gen.append(i + 1)
                gen.append(i)
                gen += list(range(i + 2, vertices))
                gens.append(gen)
            for gg in canaug_traverse_edge(g, gens, property, loops=loops, sparse=sparse):
                if extra_property(gg):
                    yield copyfun(gg) if copy else gg
        else:
            raise NotImplementedError

    def nauty_geng(self, options='', debug=False):
        r"""
        Return a generator which creates graphs from nauty's geng program.

        INPUT:

        - ``options`` -- string (default: ``''``); a string passed to ``geng``
          as if it was run at a system command line. At a minimum, you *must*
          pass the number of vertices you desire.  Sage expects the graphs to be
          in nauty's "graph6" format, do not set an option to change this
          default or results will be unpredictable.

        - ``debug`` -- boolean (default: ``False``); if ``True`` the first line
          of ``geng``'s output to standard error is captured and the first call
          to the generator's ``next()`` function will return this line as a
          string.  A line leading with ">A" indicates a successful initiation of
          the program with some information on the arguments, while a line
          beginning with ">E" indicates an error with the input.

        The possible options, obtained as output of ``geng --help``::

                 n       : the number of vertices
            mine:maxe    : <int>:<int> a range for the number of edges
                            <int>:0 means '<int> or more' except in the case 0:0
              res/mod : only generate subset res out of subsets 0..mod-1

                -c       : only write connected graphs
                -C       : only write biconnected graphs
                -t       : only generate triangle-free graphs
                -f       : only generate 4-cycle-free graphs
                -b       : only generate bipartite graphs
                              (-t, -f and -b can be used in any combination)
                -m       : save memory at the expense of time (only makes a
                              difference in the absence of -b, -t, -f and n <= 28).
                -d<int>  : a lower bound for the minimum degree
                -D<int>  : a upper bound for the maximum degree
                -v       : display counts by number of edges
                -l       : canonically label output graphs

                -q       : suppress auxiliary output (except from -v)

        Options which cause ``geng`` to use an output format different than the
        graph6 format are not listed above (-u, -g, -s, -y, -h) as they will
        confuse the creation of a Sage graph.  The res/mod option can be useful
        when using the output in a routine run several times in parallel.

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

            sage: gen = graphs.nauty_geng("2")
            sage: next(gen)
            Graph on 2 vertices
            sage: next(gen)
            Graph on 2 vertices
            sage: next(gen)
            Traceback (most recent call last):
            ...
            StopIteration

        A list of all graphs on 7 vertices.  This agrees with
        :oeis:`A000088`.  ::

            sage: gen = graphs.nauty_geng("7")
            sage: len(list(gen))
            1044

        A list of just the connected graphs on 7 vertices.  This agrees with
        :oeis:`A001349`.  ::

            sage: gen = graphs.nauty_geng("7 -c")
            sage: len(list(gen))
            853

        A list of connected degree exactly 2 graphs on 5 vertices. ::

            sage: gen = graphs.nauty_geng("5 -c -d2 -D2")
            sage: len(list(gen))
            1

        The ``debug`` switch can be used to examine ``geng``'s reaction to the
        input in the ``options`` string.  We illustrate success.  (A failure
        will be a string beginning with ">E".)  Passing the "-q" switch to
        ``geng`` will suppress the indicator of a successful initiation, and so
        the first returned value might be an empty string if ``debug`` is
        ``True``::

            sage: gen = graphs.nauty_geng("4", debug=True)
            sage: print(next(gen))
            >A ...geng -d0D3 n=4 e=0-6
            sage: gen = graphs.nauty_geng("4 -q", debug=True)
            sage: next(gen)
            ''

        TESTS:

        Wrong input, ``"-c3"`` instead of ``"-c 3"`` (:issue:`14068`)::

            sage: list(graphs.nauty_geng("-c3", debug=False))
            Traceback (most recent call last):
            ...
            ValueError: wrong format of parameter option
            sage: list(graphs.nauty_geng("-c3", debug=True))
            ['>E Usage: ...geng ...\n']
            sage: list(graphs.nauty_geng("-c 3", debug=True))
            ['>A ...geng -cd1D2 n=3 e=2-3\n', Graph on 3 vertices, Graph on 3 vertices]
        """
        import shlex
        from sage.features.nauty import NautyExecutable
        geng_path = NautyExecutable("geng").absolute_filename()
        sp = subprocess.Popen(shlex.quote(geng_path) + " {0}".format(options), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True,
                              encoding='latin-1')
        msg = sp.stderr.readline()
        if debug:
            yield msg
        elif msg.startswith('>E'):
            raise ValueError('wrong format of parameter option')
        gen = sp.stdout
        while True:
            try:
                s = next(gen)
            except StopIteration:
                # Exhausted list of graphs from nauty geng
                return
            G = graph.Graph(s[:-1], format='graph6')
            yield G

    def nauty_genbg(self, options='', debug=False):
        r"""
        Return a generator which creates bipartite graphs from nauty's ``genbgL``
        program.

        INPUT:

        - ``options`` -- string (default: ``""``); a string passed to ``genbgL``
          as if it was run at a system command line. At a minimum, you *must*
          pass the number of vertices you desire in each side. Sage expects the
          bipartite graphs to be in nauty's "graph6" format, do not set an
          option to change this default or results will be unpredictable.

        - ``debug`` -- boolean (default: ``False``); if ``True`` the first line
          of ``geng``'s output to standard error is captured and the first call
          to the generator's ``next()`` function will return this line as a
          string. A line leading with ">A" indicates a successful initiation of
          the program with some information on the arguments, while a line
          beginning with ">E" indicates an error with the input.

        The possible options, obtained as output of ``genbgL --help``::

                n1       : the number of vertices in the first class.
                           We must have n1=1..30.
                n2       : the number of vertices in the second class.
                           We must have n2=0..64 and n1+n2=1..64.
            mine:maxe    : <int>:<int> a range for the number of edges
                            <int>:0 means '<int> or more' except in the case 0:0
              res/mod    : only generate subset res out of subsets 0..mod-1
                -c       : only write connected graphs
                -z       : all the vertices in the second class must have
                           different neighbourhoods
                -F       : the vertices in the second class must have at least
                           two neighbours of degree at least 2
                -L       : there is no vertex in the first class whose removal
                           leaves the vertices in the second class unreachable
                           from each other
                -Y<int>  : two vertices in the second class must have at least
                           <int> common neighbours
                -Z<int>  : two vertices in the second class must have at most
                           <int> common neighbours
                -A       : no vertex in the second class has a neighbourhood
                           which is a subset of another vertex's neighbourhood
                           in the second class
                -D<int>  : specify an upper bound for the maximum degree
                           Example: -D6. You can also give separate maxima for
                           the two parts, for example: -D5:6
                -d<int>  : specify a lower bound for the minimum degree
                           Again, you can specify it separately for the two parts,
                           for example -d1:2
                -v       : display counts by number of edges to stderr
                -l       : canonically label output graphs

        Options which cause ``genbgL`` to use an output format different than
        the ``graph6`` format are not listed above (``-s``, ``-a``) as they will
        confuse the creation of a Sage graph. Option ``-q`` which suppress
        auxiliary output (except from ``-v``) should never be used as we are
        unable to recover the partition of the vertices of the bipartite graph
        without the auxiliary output. Hence the partition of the vertices of
        returned bipartite graphs might not respect the requirement.

        The res/mod option can be useful when using the output in a routine run
        several times in parallel.

        OUTPUT:

        A generator which will produce the graphs as
        :class:`~sage/graphs.bipartite_graph.BipartiteGraph`. These will be
        simple bipartite graphs: no loops, no multiple edges, no directed edges.

        EXAMPLES:

        The generator can be used to construct biparrtite graphs for testing,
        one at a time (usually inside a loop).  Or it can be used to
        create an entire list all at once if there is sufficient memory
        to contain it::

            sage: gen = graphs.nauty_genbg("1 1")
            sage: next(gen)
            Bipartite graph on 2 vertices
            sage: next(gen)
            Bipartite graph on 2 vertices
            sage: next(gen)
            Traceback (most recent call last):
            ...
            StopIteration

        Connected bipartite graphs of order 6 with different number of vertices
        in each side::

            sage: gen = graphs.nauty_genbg("1 5 -c")
            sage: len(list(gen))
            1
            sage: gen = graphs.nauty_genbg("2 4 -c")
            sage: len(list(gen))
            6
            sage: gen = graphs.nauty_genbg("3 3 -c")
            sage: len(list(gen))
            13

        Use :meth:`nauty_geng` instead if you want the list of all bipartite
        graphs of order `n`. For instance, the list of all connected bipartite
        graphs of order 6, which agrees with :oeis:`A005142`::

            sage: gen = graphs.nauty_geng("-b -c 6")
            sage: len(list(gen))
            17

        The ``debug`` switch can be used to examine ``genbgL``'s reaction to the
        input in the ``options`` string. A message starting with ">A" indicates
        success and a message starting with ">E" indicates a failure::

            sage: gen = graphs.nauty_genbg("2 3", debug=True)
            sage: print(next(gen))
            >A ...genbg... n=2+3 e=0:6 d=0:0 D=3:2
            sage: gen = graphs.nauty_genbg("-c2 3", debug=True)
            sage: next(gen)
            '>E Usage: ...genbg... [-c -ugs -vq -lzF] [-Z#] [-D#] [-A] [-d#|-d#:#] [-D#|-D#:#] n1 n2...

        Check that the partition of the bipartite graph is consistent::

            sage: gen = graphs.nauty_genbg("3 3")
            sage: left = set(range(3))
            sage: for g in gen:
            ....:     if g.left != left:
            ....:         raise ValueError('wrong partition')

        TESTS:

        Wrong input::

            sage: list(graphs.nauty_genbg("-c1 2", debug=False))
            Traceback (most recent call last):
            ...
            ValueError: wrong format of parameter options
            sage: list(graphs.nauty_genbg("-c1 2", debug=True))
            ['>E Usage: ...genbg... [-c -ugs -vq -lzF] [-Z#] [-D#] [-A] [-d#|-d#:#] [-D#|-D#:#] n1 n2...
            sage: list(graphs.nauty_genbg("-c 1 2", debug=True))
            ['>A ...genbg... n=1+2 e=2:2 d=1:1 D=2:1 c...\n', Bipartite graph on 3 vertices]

        We must have n1=1..30, n2=0..64 and n1+n2=1..64 (:issue:`34179`,
        :issue:`38618`)::

            sage: next(graphs.nauty_genbg("31 1", debug=False))
            Traceback (most recent call last):
            ...
            ValueError: wrong format of parameter options
            sage: next(graphs.nauty_genbg("31 1", debug=True))
            '>E ...genbg...: must have n1=1..30, n1+n2=1..64...
            sage: next(graphs.nauty_genbg("30 40", debug=True))
            '>E ...genbg...: must have n1=1..30, n1+n2=1..64...
            sage: next(graphs.nauty_genbg("1 63", debug=False))
            Bipartite graph on 64 vertices
            sage: next(graphs.nauty_genbg("1 64", debug=True))
            '>E ...genbg...: must have n1=1..30, n1+n2=1..64...
            sage: next(graphs.nauty_genbg("0 2", debug=True))
            '>E ...genbg...: must have n1=1..30, n1+n2=1..64...
            sage: next(graphs.nauty_genbg("2 0", debug=False))
            Bipartite graph on 2 vertices
            sage: next(graphs.nauty_genbg("2 -1", debug=True))
            '>E Usage: ...genbg... [-c -ugs -vq -lzF] [-Z#] [-D#] [-A] [-d#|-d#:#] [-D#|-D#:#] n1 n2...
        """
        import shlex
        from sage.features.nauty import NautyExecutable
        genbg_path = NautyExecutable("genbgL").absolute_filename()
        sp = subprocess.Popen(shlex.quote(genbg_path) + " {0}".format(options), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True,
                              encoding='latin-1')
        msg = sp.stderr.readline()
        if debug:
            yield msg
        elif msg.startswith('>E'):
            raise ValueError('wrong format of parameter options')

        if msg.startswith('>A'):
            # We extract the partition of the vertices from the msg string
            for s in msg.split(' '):
                if s.startswith('n='):
                    from sage.rings.integer import Integer
                    n1, n2 = (Integer(t) for t in s[2:].split('+') if t.isdigit())
                    partition = [set(range(n1)), set(range(n1, n1 + n2))]
                    break
            else:
                # should never happen
                raise ValueError('unable to recover the partition')
        else:
            # Either msg starts with >E or option -q has been given
            partition = None

        gen = sp.stdout
        from sage.graphs.bipartite_graph import BipartiteGraph
        while True:
            try:
                s = next(gen)
            except StopIteration:
                # Exhausted list of bipartite graphs from nauty genbgL
                return
            G = BipartiteGraph(s[:-1], format='graph6', partition=partition)
            yield G

    def nauty_genktreeg(self, options='', debug=False):
        r"""
        Return a generator which creates all `k`-trees using nauty..

        A `k`-tree is an undirected graph formed by starting with a complete
        graph on `k + 1` vertices and then repeatedly add vertices in such a
        way that each added vertex `v` has exactly `k` neighbors `U` such that,
        together, the `k + 1` vertices formed by `v` and `U` form a clique.
        See the :wikipedia:`K-tree` for more details.

        INPUT:

        - ``options`` -- string (default: ``""``); a string passed to
          ``genktreeg`` as if it was run at a system command line. At a minimum,
          you *must* pass the number of vertices you desire. Sage expects the
          graphs to be in nauty's "graph6" format, do not set an option to
          change this default or results will be unpredictable.

        - ``debug`` -- boolean (default: ``False``); if ``True`` the first line
          of ``genktreeg``'s output to standard error is captured and the first
          call to the generator's ``next()`` function will return this line as a
          string. A line leading with ">A" indicates a successful initiation of
          the program with some information on the arguments, while a line
          beginning with ">E" indicates an error with the input.

        The possible options, obtained as output of ``genktreeg --help``::

                 n       : the number of vertices
                -k<int>  : the value of `k`(default: 2)
              res/mod    : only generate subset res out of subsets 0..mod-1
                -l       : canonically label output graphs

        Options which cause ``genktreeg`` to use an output format different than
        the graph6 format are not listed above (-u, -s, -h) as they will confuse
        the creation of a Sage graph. The res/mod option can be useful when
        using the output in a routine run several times in parallel.

        OUTPUT:

        A generator which will produce the graphs as Sage graphs.
        These will be simple graphs: no loops, no multiple edges, no
        directed edges.

        EXAMPLES:

        A `k`-tree is a maximal graph with treewidth `k`::

            sage: # needs nauty
            sage: gen = graphs.nauty_genktreeg("10 -k4")
            sage: G = next(gen); G
            Graph on 10 vertices
            sage: G.treewidth()
            4

        A list of all 2-trees with 6, 7 and 8 vertices. This agrees with
        :oeis:`A054581`::

            sage: # needs nauty
            sage: gen = graphs.nauty_genktreeg("6")
            sage: len(list(gen))
            5
            sage: gen = graphs.nauty_genktreeg("7")
            sage: len(list(gen))
            12
            sage: gen = graphs.nauty_genktreeg("8")
            sage: len(list(gen))
            39

        The ``debug`` switch can be used to examine ``geng``'s reaction to the
        input in the ``options`` string.  We illustrate success.  (A failure
        will be a string beginning with ">E".)  Passing the "-q" switch to
        ``geng`` will suppress the indicator of a successful initiation, and so
        the first returned value might be an empty string if ``debug`` is
        ``True``::

            sage: gen = graphs.nauty_genktreeg("7", debug=True)                         # needs nauty
            sage: print(next(gen))                                                      # needs nauty
            >A ...genktreeg k=2 n=7

        TESTS:

        Wrong input::

            sage: # needs nauty
            sage: list(graphs.nauty_genktreeg("4 -k5", debug=True))
            ['>E genktreeg: n cannot be less than k\n']
            sage: list(graphs.nauty_genktreeg("10 -k 4", debug=True))
            ['>E genktreeg -k: missing argument value\n']
            sage: list(graphs.nauty_genktreeg("-c3", debug=False))
            Traceback (most recent call last):
            ...
            ValueError: wrong format of parameter option
        """
        import shlex
        from sage.features.nauty import NautyExecutable
        geng_path = NautyExecutable("genktreeg").absolute_filename()
        sp = subprocess.Popen(shlex.quote(geng_path) + " {0}".format(options), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True,
                              encoding='latin-1')
        msg = sp.stderr.readline()
        if debug:
            yield msg
        elif msg.startswith('>E'):
            raise ValueError('wrong format of parameter option')
        gen = sp.stdout
        while True:
            try:
                s = next(gen)
            except StopIteration:
                # Exhausted list of graphs from nauty geng
                return
            G = graph.Graph(s[:-1], format='graph6')
            yield G

    def cospectral_graphs(self, vertices, matrix_function=None, graphs=None):
        r"""
        Find all sets of graphs on ``vertices`` vertices (with
        possible restrictions) which are cospectral with respect to a
        constructed matrix.

        INPUT:

        - ``vertices`` -- the number of vertices in the graphs to be tested

        - ``matrix_function`` -- a function taking a graph and giving back
          a matrix.  This defaults to the adjacency matrix.  The spectra
          examined are the spectra of these matrices.

        - ``graphs`` -- one of three things:

           - ``None`` -- default; test all graphs having ``vertices``
             vertices

           - a function taking a graph and returning ``True`` or ``False``
             - test only the graphs on ``vertices`` vertices for which
             the function returns ``True``

           - a list of graphs (or other iterable object) -- these graphs
             are tested for cospectral sets.  In this case,
             ``vertices`` is ignored.

        OUTPUT:

           A list of lists of graphs.  Each sublist will be a list of
           cospectral graphs (lists of cardinality 1 being omitted).

        .. SEEALSO::

            :meth:`Graph.is_strongly_regular` -- tests whether a graph is
            strongly regular and/or returns its parameters.

        EXAMPLES::

            sage: g = graphs.cospectral_graphs(5)                                       # needs sage.modules
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)       # needs sage.modules
            [['Dr?', 'Ds_']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()                      # needs sage.modules
            True

        There are two sets of cospectral graphs on six vertices with no isolated vertices::

            sage: # needs sage.modules
            sage: g = graphs.cospectral_graphs(6, graphs=lambda x: min(x.degree())>0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)
            [['Ep__', 'Er?G'], ['ExGg', 'ExoG']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()
            True
            sage: g[1][1].am().charpoly()==g[1][1].am().charpoly()
            True

        There is one pair of cospectral trees on eight vertices::

            sage: g = graphs.cospectral_graphs(6, graphs=graphs.trees(8))               # needs sage.modules
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)       # needs sage.modules
            [['GiPC?C', 'GiQCC?']]
            sage: g[0][1].am().charpoly()==g[0][1].am().charpoly()                      # needs sage.modules
            True

        There are two sets of cospectral graphs (with respect to the
        Laplacian matrix) on six vertices::

            sage: # needs sage.modules
            sage: g = graphs.cospectral_graphs(6, matrix_function=lambda g: g.laplacian_matrix())
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
            ....:   A = g.adjacency_matrix().change_ring(QQ)
            ....:   for i in range(g.order()):
            ....:       A.rescale_row(i, 1 / len(A.nonzero_positions_in_row(i)))
            ....:   return A
            sage: g = graphs.cospectral_graphs(5, matrix_function=DinverseA,            # needs sage.libs.pari sage.modules
            ....:                              graphs=lambda g: min(g.degree()) > 0)
            sage: sorted(sorted(g.graph6_string() for g in glist) for glist in g)       # needs sage.modules
            [['Dlg', 'Ds_']]
            sage: (g[0][1].laplacian_matrix(normalized=True).charpoly()                 # needs sage.modules sage.symbolic
            ....:   == g[0][1].laplacian_matrix(normalized=True).charpoly())
            True
        """
        if matrix_function is None:
            matrix_function = lambda g: g.adjacency_matrix()

        from sage.graphs.graph_generators import graphs as graph_gen
        if graphs is None:
            graph_list = graph_gen(vertices, property=lambda _: True)
        elif callable(graphs):
            graph_list = iter(g for g in graph_gen(vertices, property=lambda _: True) if graphs(g))
        else:
            graph_list = iter(graphs)

        from collections import defaultdict
        charpolys = defaultdict(list)
        for g in graph_list:
            cp = matrix_function(g).charpoly()
            charpolys[cp].append(g)

        cospectral_graphs = []
        for cp, g_list in charpolys.items():
            if len(g_list) > 1:
                cospectral_graphs.append(g_list)

        return cospectral_graphs

    def _read_planar_code(self, code_input):
        r"""
        Return a generator for the plane graphs in planar code format in
        the file code_input (see [BM2016]_).

        A file with planar code starts with a header ``>>planar_code<<``.
        After the header each graph is stored in the following way :

        The first character is the number of vertices, followed by
        n11,...,n1k,null character,n21,...,n2k',null character, ...

        where the n1* are all neighbors of n1 and all n2* are the
        neighbors of n2, ...
        Besides, these neighbors are enumerated in clockwise order.

        INPUT:

        - ``code_input`` -- a file containing valid planar code data

        OUTPUT:

        A generator which will produce the plane graphs as Sage graphs
        with an embedding set. These will be simple graphs: no loops, no
        multiple edges, no directed edges (unless plantri is asked to give
        the dual graphs instead).

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        The following example creates a small planar code file in memory and
        reads it using the ``_read_planar_code`` method::

            sage: from io import StringIO
            sage: code_input = StringIO('>>planar_code<<')
            sage: _ = code_input.write('>>planar_code<<')
            sage: for c in [4,2,3,4,0,1,4,3,0,1,2,4,0,1,3,2,0]:
            ....:     _ = code_input.write('{:c}'.format(c))
            sage: _ = code_input.seek(0)
            sage: gen = graphs._read_planar_code(code_input)
            sage: l = list(gen)
            sage: l
            [Graph on 4 vertices]
            sage: l[0].is_isomorphic(graphs.CompleteGraph(4))
            True
            sage: l[0].get_embedding()
            {1: [2, 3, 4],
             2: [1, 4, 3],
             3: [1, 2, 4],
             4: [1, 3, 2]}
        """
        # start of code to read planar code
        header = code_input.read(15)
        assert header == '>>planar_code<<', 'Not a valid planar code header'

        # read graph per graph
        while True:
            c = code_input.read(1)
            if not c:
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
                c = code_input.read(1)
                if ord(c) == 0:
                    zeroCount += 1
                else:
                    g[zeroCount].append(ord(c))

            # construct graph based on g

            # first taking care that every edge is given twice
            edges_g = {i + 1: [j for j in di if j < i + 1]
                       for i, di in enumerate(g)}

            # then adding half of the loops (if any)
            has_loops = False
            for i, di in enumerate(g):
                Ni = di.count(i + 1)
                if Ni > 1:
                    edges_g[i + 1] += [i + 1] * (Ni // 2)
                    has_loops = True
            G = graph.Graph(edges_g, loops=has_loops)

            if not (G.has_multiple_edges() or has_loops):
                embed_g = {i + 1: di for i, di in enumerate(g)}
                G.set_embedding(embed_g)
            yield G

    def fullerenes(self, order, ipr=False):
        r"""
        Return a generator which creates fullerene graphs using
        the buckygen generator (see [BGM2012]_).

        INPUT:

        - ``order`` -- a positive even integer smaller than or equal to 254
          This specifies the number of vertices in the generated fullerenes

        - ``ipr`` -- (default: ``False``) if ``True`` only fullerenes that
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
        on 60 vertices::

            sage: gen = graphs.fullerenes(60)  # optional - buckygen
            sage: len(list(gen))               # optional - buckygen
            1812

        However, there is only one IPR fullerene graph on 60 vertices: the famous
        Buckminster Fullerene::

            sage: gen = graphs.fullerenes(60, ipr=True)  # optional - buckygen
            sage: next(gen)                              # optional - buckygen
            Graph on 60 vertices
            sage: next(gen)                              # optional - buckygen
            Traceback (most recent call last):
            ...
            StopIteration

        The unique fullerene graph on 20 vertices is isomorphic to the dodecahedron
        graph. ::

            sage: # optional - buckygen
            sage: gen = graphs.fullerenes(20)
            sage: g = next(gen)
            sage: g.is_isomorphic(graphs.DodecahedralGraph())
            True
            sage: g.get_embedding()
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
            sage: g.plot3d(layout='spring')
            Graphics3d Object
        """
        # number of vertices should be positive
        if order < 0:
            raise ValueError("number of vertices should be nonnegative")

        # buckygen can only output fullerenes on up to 254 vertices
        if order > 254:
            raise ValueError("number of vertices should be at most 254")

        # fullerenes only exist for an even number of vertices, larger than 20
        # and different from 22
        if order % 2 == 1 or order < 20 or order == 22:
            return

        from sage.features.graph_generators import Buckygen
        Buckygen().require()

        import shlex
        command = shlex.quote(Buckygen().absolute_filename())
        command += ' -' + ('I' if ipr else '') + 'd {0}d'.format(order)

        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True,
                              encoding='latin-1')

        sp.stdout.reconfigure(newline='')

        yield from graphs._read_planar_code(sp.stdout)

    def fusenes(self, hexagon_count, benzenoids=False):
        r"""
        Return a generator which creates fusenes and benzenoids using
        the benzene generator (see [BCH2002]_). Fusenes are planar
        polycyclic hydrocarbons with all bounded faces hexagons. Benzenoids
        are fusenes that are subgraphs of the hexagonal lattice.

        INPUT:

        - ``hexagon_count`` -- positive integer smaller than or equal to 30;
          this specifies the number of hexagons in the generated benzenoids

        - ``benzenoids`` -- (default: ``False``) if ``True`` only benzenoids are
          generated

        OUTPUT:

        A generator which will produce the fusenes as Sage graphs
        with an embedding set. These will be simple graphs: no loops, no
        multiple edges, no directed edges.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        There is a unique fusene with 2 hexagons::

            sage: gen = graphs.fusenes(2)  # optional - benzene
            sage: len(list(gen))           # optional - benzene
            1

        This fusene is naphthalene (`\textrm{C}_{10}\textrm{H}_{8}`).
        In the fusene graph the H-atoms are not stored, so this is
        a graph on just 10 vertices::

            sage: gen = graphs.fusenes(2)  # optional - benzene
            sage: next(gen)                # optional - benzene
            Graph on 10 vertices
            sage: next(gen)                # optional - benzene
            Traceback (most recent call last):
            ...
            StopIteration

        There are 6505 benzenoids with 9 hexagons::

            sage: gen = graphs.fusenes(9, benzenoids=True)  # optional - benzene
            sage: len(list(gen))                            # optional - benzene
            6505
        """
        if hexagon_count < 0:
            raise ValueError("number of hexagons should be nonnegative")

        # benzene is only built for fusenes with up to 30 hexagons
        if hexagon_count > 30:
            raise ValueError("number of hexagons should be at most 30")

        # there are no fusenes with 0 hexagons
        if hexagon_count == 0:
            return

        # there is only one unique fusene with 1 hexagon (and benzene doesn't generate it)
        if hexagon_count == 1:
            g = {1: [6, 2], 2: [1, 3], 3: [2, 4], 4: [3, 5], 5: [4, 6], 6: [5, 1]}
            G = graph.Graph(g)
            G.set_embedding(g)
            yield G
            return

        from sage.features.graph_generators import Benzene
        Benzene().require()

        import shlex
        command = shlex.quote(Benzene().absolute_filename())
        command += (' b' if benzenoids else '') + ' {0} p'.format(hexagon_count)

        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True,
                              encoding='latin-1')

        sp.stdout.reconfigure(newline='')

        for G in graphs._read_planar_code(sp.stdout):
            yield G

    def plantri_gen(self, options=""):
        r"""
        Iterator over planar graphs created using the ``plantri`` generator.

        ``plantri`` is a (optional) program that generates certain types of
        graphs that are embedded on the sphere. It outputs exactly one member of
        each isomorphism class, using an amount of memory almost independent of
        the number of graphs produced. Isomorphisms are defined with respect to
        the embeddings, so in some cases outputs may be isomorphic as abstract
        graphs.

        This method allows for passing command directly to ``plantry``,
        similarly to method :meth:`nauty_geng`, provide that the output format
        is not changed.

        INPUT:

        - ``options`` -- string (default: ``""``); a string passed to
          ``plantri`` as if it was run at a system command line. At a minimum,
          you *must* pass the number of vertices you desire. Sage expects the
          output of plantri to be in "planar code" format, so do not set an
          option to change this default or results will be unpredictable.

        The possible options are::

            n       : the number of vertices (the only compulsory parameter).
                      This number must be in range `3\cdots 64`.
                      It can also be given as "nd", where the suffix "d" means
                      "dual", in which case it is converted by adding 4 then
                      dividing by 2, i.e., `(28+4)/2 = 16`. In the case of
                      triangulations, this calculation yields the number of
                      faces, which is the number of vertices in the dual cubic
                      graph.

            -d      : output the dual instead of the original graph.
                      Note that it is applied only at the output stage. All
                      other switches refer to the original graph before the dual
                      is taken.

            -o      : Normally, one member of each isomorphism class is written.
                      If this switch is given, one member of each O-P
                      isomorphism class is written.

            -V      : output only graphs with non-trivial group. If -o is
                      given the O-P group is used, the full group otherwise.

            -m<int> : lower bound on the minimum degree. The default is -m3.
                      In the dual graph, this means a lower bound on the minimum
                      face size.

            -c<int> : lower bound on the connectivity. The default is -c3.

            -x      : when used in combination with -cN, the connectivity must
                      be exactly N rather than at least N.

            -e      : used to specify bounds on the number of edges.
                      There are four possible forms:
                          -e<int>        exactly <int> edges
                          -e:<int>       at most <int> edges
                          -e<int>:       at least <int> edges
                          -e<int>:<int>  between <int> and <int> edges

            -f<int> : upper bound on the size of a face, and so on the maximum
                      degree of the dual.

            -b but not -p : select eulerian triangulations, where "eulerian"
                            means that every vertex has even degree.
                            This parameter can be used in combination with
                            parameters -c and -x.

            -p but not -b : select general planar simple graphs.
                            This parameter can be used in combination with
                            parameters -m, -c, -x, -e and -f.

            -bp or -pb    : select general planar simple bipartite graphs.
                            This parameter can be used in combination with
                            parameters -m, -c, -x, -e and -f, except -c4, -m4,
                            -m5 and -f3.

            -P<int> : select triangulations of a disk. These are embedded simple
                      graphs with a distinguished "outer" face. The outer face
                      can be of any size (here called the disk size) but the
                      other faces must be triangles.  The argument <int> to -P
                      is the disk size. If no argument (or 0) is given, all disk
                      sizes are permitted.
                      This parameter can be used in combination with
                      parameters -m, -c, and -x.

            -q      : select simple quadrangulations. These are planar simple
                      graphs for which every face has length 4.
                      This parameter can be used in combination with parameters
                      -c and -m.

            -A      : select Appolonian networks. These are simple planar
                      triangulations that can be formed starting with `K_4` then
                      repeatedly dividing a face into three by addition of a new
                      vertex. They all have minimum degree and connectivity
                      equal to 3.

            res/mod : only generate subset res out of subsets 0..mod-1.
                      The set of objects is divided into mod disjoint classes
                      and only the res-th class is generated.

        If -b, -q, -p, -P and -A are absent, the graphs found are triangulations
        only restricted by connectivity and minimum degree. In this case,
        there is the possibility of connectivity lower than 3.

        Other options listed in the ``plantri`` guide might cause unpredictable
        behavior, in particular those changing the output format of ``plantri``
        as they will confuse the creation of a Sage graph.

        OUTPUT:

        An iterator which yields the graphs generated by ``plantri`` as Sage
        :class:`~sage.graphs.graph.Graph`.

        .. SEEALSO::

            - :meth:`planar_graphs` -- iterator over connected planar graphs
              using the ``plantri`` generator
            - :meth:`triangulations` -- iterator over connected planar
              triangulations using the ``plantri`` generator
            - :meth:`quadrangulations` -- iterator over connected planar
              quadrangulations using the ``plantri`` generator

        EXAMPLES:

        The generator can be used to construct graphs for testing, one at a time
        (usually inside a loop). Or it can be used to create an entire list all
        at once if there is sufficient memory to contain it::

            sage: # optional - plantri
            sage: gen = graphs.plantri_gen("6")
            sage: next(gen)
            Graph on 6 vertices
            sage: next(gen)
            Graph on 6 vertices
            sage: next(gen)
            Traceback (most recent call last):
            ...
            StopIteration

        An overview of the number of quadrangulations on up to 12 vertices. This
        agrees with :oeis:`A113201`::

            sage: for i in range(4, 13):                        # optional - plantri
            ....:     cmd = '-qm2c2 {}'.format(i)
            ....:     L = len(list(graphs.plantri_gen(cmd)))
            ....:     print("{:2d}   {:3d}".format(i, L))
             4     1
             5     1
             6     2
             7     3
             8     9
             9    18
            10    62
            11   198
            12   803

        TESTS:

        Wrong input, ``"-c=3"`` instead of ``"-c3"``::

            sage: list(graphs.plantri_gen("6 -c3"))  # optional - plantri
            [Graph on 6 vertices, Graph on 6 vertices]
            sage: list(graphs.plantri_gen("6 -c=3"))  # optional - plantri
            Traceback (most recent call last):
            ...
            AttributeError: invalid options '6 -c=3'
        """
        from sage.features.graph_generators import Plantri
        Plantri().require()

        import shlex
        command = '{} {}'.format(shlex.quote(Plantri().absolute_filename()),
                                 options)
        sp = subprocess.Popen(command, shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True,
                              encoding='latin-1')

        sp.stdout.reconfigure(newline='')

        try:
            yield from graphs._read_planar_code(sp.stdout)
        except AssertionError:
            raise AttributeError("invalid options '{}'".format(options))

    def planar_graphs(self, order, minimum_degree=None,
                      minimum_connectivity=None,
                      exact_connectivity=False,
                      minimum_edges=None,
                      maximum_edges=None,
                      maximum_face_size=None,
                      only_bipartite=False,
                      dual=False):
        r"""
        An iterator over connected planar graphs using the plantri generator.

        This uses the plantri generator (see [BM2007]_) which is available
        through the optional package plantri.

        .. NOTE::

            The non-3-connected graphs will be returned several times, with all
            its possible embeddings.

        INPUT:

        - ``order`` -- positive integer smaller than or equal to 64;
          this specifies the number of vertices in the generated graphs

        - ``minimum_degree`` -- (default: ``None``) a value `\geq 1` and `\leq
          5`, or ``None``. This specifies the minimum degree of the generated
          graphs. If this is ``None`` and the order is 1, then this is set to
          0. If this is ``None`` and the minimum connectivity is specified, then
          this is set to the same value as the minimum connectivity.  If the
          minimum connectivity is also equal to ``None``, then this is set to 1.

        - ``minimum_connectivity`` -- (default: ``None``) a value `\geq 1`
          and `\leq 3`, or ``None``. This specifies the minimum connectivity of the
          generated graphs. If this is ``None`` and the minimum degree is
          specified, then this is set to the minimum of the minimum degree
          and 3. If the minimum degree is also equal to ``None``, then this
          is set to 1.

        - ``exact_connectivity`` -- (default: ``False``) if ``True`` only
          graphs with exactly the specified connectivity will be generated.
          This option cannot be used with ``minimum_connectivity=3``, or if
          the minimum connectivity is not explicitly set.

        - ``minimum_edges`` -- integer (default: ``None``); lower bound on the
          number of edges

        - ``maximum_edges`` -- integer (default: ``None``); upper bound on the
          number of edges

        - ``maximum_face_size`` -- integer (default: ``None``); upper bound on
          the size of a face and so on the maximum degree of the dual graph

        - ``only_bipartite`` -- (default: ``False``) if ``True`` only bipartite
          graphs will be generated. This option cannot be used for graphs with
          a minimum degree larger than 3.

        - ``dual`` -- (default: ``False``) if ``True`` return instead the
          planar duals of the generated graphs

        OUTPUT:

        An iterator which will produce all planar graphs with the given
        number of vertices as Sage graphs with an embedding set. These will be
        simple graphs (no loops, no multiple edges, no directed edges)
        unless the option ``dual=True`` is used.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        There are 6 planar graphs on 4 vertices::

            sage: gen = graphs.planar_graphs(4)  # optional - plantri
            sage: len(list(gen))                 # optional - plantri
            6

        Three of these planar graphs are bipartite::

            sage: gen = graphs.planar_graphs(4, only_bipartite=True)  # optional - plantri
            sage: len(list(gen))                                      # optional - plantri
            3

        Setting ``dual=True`` gives the planar dual graphs::

            sage: gen = graphs.planar_graphs(4, dual=True)  # optional - plantri
            sage: [u for u in list(gen)]                    # optional - plantri
            [Graph on 4 vertices,
            Multi-graph on 3 vertices,
            Multi-graph on 2 vertices,
            Looped multi-graph on 2 vertices,
            Looped multi-graph on 1 vertex,
            Looped multi-graph on 1 vertex]

        The cycle of length 4 is the only 2-connected bipartite planar graph
        on 4 vertices::

            sage: l = list(graphs.planar_graphs(4, minimum_connectivity=2, only_bipartite=True))  # optional - plantri
            sage: l[0].get_embedding()                                                            # optional - plantri
            {1: [2, 3],
             2: [1, 4],
             3: [1, 4],
             4: [2, 3]}

        There is one planar graph with one vertex. This graph obviously has
        minimum degree equal to 0::

            sage: list(graphs.planar_graphs(1))                    # optional - plantri
            [Graph on 1 vertex]
            sage: list(graphs.planar_graphs(1, minimum_degree=1))  # optional - plantri
            []

        Specifying lower and upper bounds on the number of edges::

            sage: # optional - plantri
            sage: len(list(graphs.planar_graphs(4)))
            6
            sage: len(list(graphs.planar_graphs(4, minimum_edges=4)))
            4
            sage: len(list(graphs.planar_graphs(4, maximum_edges=4)))
            4
            sage: len(list(graphs.planar_graphs(4, minimum_edges=4, maximum_edges=4)))
            2

        Specifying the maximum size of a face::

            sage: len(list(graphs.planar_graphs(4, maximum_face_size=3)))  # optional - plantri
            1
            sage: len(list(graphs.planar_graphs(4, maximum_face_size=4)))  # optional - plantri
            3

        TESTS:

        The number of edges in a planar graph is equal to the number of edges in
        its dual::

            sage: # optional - plantri
            sage: planar      = list(graphs.planar_graphs(5,dual=True))
            sage: dual_planar = list(graphs.planar_graphs(5,dual=False))
            sage: planar_sizes      = [g.size() for g in planar]
            sage: dual_planar_sizes = [g.size() for g in dual_planar]
            sage: planar_sizes == dual_planar_sizes
            True
        """
        if order < 0:
            raise ValueError("number of vertices should be nonnegative")

        # plantri can only output general planar graphs on up to 64 vertices
        if order > 64:
            raise ValueError("number of vertices should be at most 64")

        if exact_connectivity and minimum_connectivity is None:
            raise ValueError("Minimum connectivity must be specified to use the exact_connectivity option.")

        if minimum_connectivity is not None and not (1 <= minimum_connectivity <= 3):
            raise ValueError("Minimum connectivity should be a number between 1 and 3.")

        # minimum degree should be None or a number between 1 and 5
        if minimum_degree == 0:
            if order != 1:
                raise ValueError("Minimum degree equal to 0 is only possible if the graphs have 1 vertex.")
        elif minimum_degree is not None and not (1 <= minimum_degree <= 5):
            raise ValueError("Minimum degree should be a number between 1 and 5 if the order is greater than 1.")
        elif minimum_degree is None and order == 1:
            minimum_degree = 0

        # check combination of values of minimum degree and minimum connectivity
        if minimum_connectivity is None:
            if minimum_degree is not None:
                minimum_connectivity = min(3, minimum_degree)
            elif minimum_degree is None:
                minimum_degree, minimum_connectivity = 1, 1
        else:
            if minimum_degree is None:
                minimum_degree = minimum_connectivity
            elif (minimum_degree < minimum_connectivity and
                  minimum_degree > 0):
                raise ValueError("Minimum connectivity can be at most the minimum degree.")

        # exact connectivity is not implemented for minimum connectivity 3
        if exact_connectivity and minimum_connectivity == 3:
            raise NotImplementedError("Generation of planar graphs with connectivity exactly 3 is not implemented.")

        if only_bipartite and minimum_degree > 3:
            raise NotImplementedError("Generation of bipartite planar graphs with minimum degree 4 or 5 is not implemented.")

        edges = ''
        if minimum_edges is None:
            if maximum_edges is not None:
                if maximum_edges < order - 1:
                    raise ValueError("the number of edges cannot be less than order - 1")
                edges = '-e:{}'.format(maximum_edges)
        else:
            if minimum_edges > 3*order - 6:
                raise ValueError("the number of edges cannot be more than 3*order - 6")
            if maximum_edges is None:
                edges = '-e{}:'.format(minimum_edges)
            elif minimum_edges > maximum_edges:
                raise ValueError("the maximum number of edges must be larger "
                                 "or equal to the minimum number of edges")
            elif minimum_edges == maximum_edges:
                edges = '-e{}'.format(minimum_edges)
            else:
                edges = '-e{}:{}'.format(minimum_edges, maximum_edges)

        faces = ''
        if maximum_face_size is not None:
            if maximum_face_size < 3:
                raise ValueError("the upper bound on the size of a face must be at least 3")
            faces = '-f{}'.format(maximum_face_size)

        if order == 0:
            return

        minimum_order = {0: 1, 1: 2, 2: 3, 3: 4, 4: 6, 5: 12}[minimum_degree]

        if order < minimum_order:
            return

        if order == 1:
            if minimum_degree == 0:
                G = graph.Graph(1)
                G.set_embedding({0: []})
                yield G
            return

        cmd = '-p{}m{}c{}{}{} {} {} {}'
        command = cmd.format('b' if only_bipartite else '',
                             minimum_degree,
                             minimum_connectivity,
                             'x' if exact_connectivity else '',
                             'd' if dual else '',
                             edges, faces,
                             order)

        yield from graphs.plantri_gen(command)

    def triangulations(self, order, minimum_degree=None, minimum_connectivity=None,
                       exact_connectivity=False, only_eulerian=False, dual=False):
        r"""
        An iterator over connected planar triangulations using the plantri generator.

        This uses the plantri generator (see [BM2007]_) which is available
        through the optional package plantri.

        INPUT:

        - ``order`` -- positive integer smaller than or equal to 64;
          this specifies the number of vertices in the generated triangulations

        - ``minimum_degree`` -- (default: ``None``) a value `\geq 3` and `\leq 5`,
          or ``None``. This specifies the minimum degree of the generated
          triangulations. If this is ``None`` and the minimum connectivity
          is specified, then this is set to the same value as the minimum
          connectivity. If the minimum connectivity is also equal to ``None``,
          then this is set to 3.

        - ``minimum_connectivity`` -- (default: ``None``) a value `\geq 3` and
          `\leq 5`, or ``None``. This specifies the minimum connectivity of the
          generated triangulations. If this is ``None`` and the minimum degree
          is specified, then this is set to the minimum of the minimum degree
          and 3. If the minimum degree is also equal to ``None``, then this is
          set to 3.

        - ``exact_connectivity`` -- (default: ``False``) if ``True`` only
          triangulations with exactly the specified connectivity will be generated.
          This option cannot be used with ``minimum_connectivity=3``, or if
          the minimum connectivity is not explicitly set.

        - ``only_eulerian`` -- (default: ``False``) if ``True`` only Eulerian
          triangulations will be generated. This option cannot be used if the
          minimum degree is explicitly set to anything else than 4.

        - ``dual`` -- (default: ``False``) if ``True`` return instead the
          planar duals of the generated graphs

        OUTPUT:

        An iterator which will produce all planar triangulations with the given
        number of vertices as Sage graphs with an embedding set. These will be
        simple graphs (no loops, no multiple edges, no directed edges).

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

            - :meth:`~sage.graphs.graph_generators.GraphGenerators.RandomTriangulation`
              -- build a random triangulation.

        EXAMPLES:

        The unique planar embedding of the `K_4` is the only planar triangulations
        on 4 vertices::

            sage: gen = graphs.triangulations(4)    # optional - plantri
            sage: [g.get_embedding() for g in gen]  # optional - plantri
            [{1: [2, 3, 4], 2: [1, 4, 3], 3: [1, 2, 4], 4: [1, 3, 2]}]

        but, of course, this graph is not Eulerian::

            sage: gen = graphs.triangulations(4, only_eulerian=True)  # optional - plantri
            sage: len(list(gen))                                      # optional - plantri
            0

        The unique Eulerian triangulation on 6 vertices is isomorphic to the octahedral
        graph. ::

            sage: gen = graphs.triangulations(6, only_eulerian=True)  # optional - plantri
            sage: g = next(gen)                                       # optional - plantri
            sage: g.is_isomorphic(graphs.OctahedralGraph())           # optional - plantri
            True

        The minimum degree of a triangulation is 3, so the method can not output
        a triangle::

            sage: list(graphs.triangulations(3))                      # optional - plantri
            []

        An overview of the number of 5-connected triangulations on up to 22 vertices. This
        agrees with :oeis:`A081621`::

            sage: for i in range(12, 23):                                             # optional - plantri
            ....:     L = len(list(graphs.triangulations(i, minimum_connectivity=5)))
            ....:     print("{}   {:3d}".format(i,L))
            12     1
            13     0
            14     1
            15     1
            16     3
            17     4
            18    12
            19    23
            20    71
            21   187
            22   627

        The minimum connectivity can be at most the minimum degree::

            sage: gen = next(graphs.triangulations(10, minimum_degree=3,     # optional - plantri
            ....:                                  minimum_connectivity=5))
            Traceback (most recent call last):
            ...
            ValueError: Minimum connectivity can be at most the minimum degree.

        There are 5 triangulations with 9 vertices and minimum degree equal to 4
        that are 3-connected, but only one of them is not 4-connected::

            sage: len([g for g in graphs.triangulations(9, minimum_degree=4,        # optional - plantri
            ....:                                       minimum_connectivity=3)])
            5
            sage: len([g for g in graphs.triangulations(9, minimum_degree=4,        # optional - plantri
            ....:                                       minimum_connectivity=3,
            ....:                                       exact_connectivity=True)])
            1

        Setting ``dual=True`` gives the planar dual graphs::

            sage: [len(g) for g in graphs.triangulations(9, minimum_degree=4,       # optional plantri
            ....:                                        minimum_connectivity=3, dual=True)]
            [14, 14, 14, 14, 14]

        TESTS::

            sage: [g.size() for g in graphs.triangulations(6, minimum_connectivity=3)]  # optional - plantri
            [12, 12]
        """
        if order < 0:
            raise ValueError("number of vertices should be nonnegative")

        # plantri can only output planar triangulations on up to 64 vertices
        if order > 64:
            raise ValueError("number of vertices should be at most 64")

        if exact_connectivity and minimum_connectivity is None:
            raise ValueError("Minimum connectivity must be specified to use the exact_connectivity option.")

        if minimum_connectivity is not None and not (3 <= minimum_connectivity <= 5):
            raise ValueError("Minimum connectivity should be None or a number between 3 and 5.")

        if minimum_degree is not None and not (3 <= minimum_degree <= 5):
            raise ValueError("Minimum degree should be None or a number between 3 and 5.")

        # for Eulerian triangulations the minimum degree is set to 4 (unless it was already specifically set)
        if only_eulerian and minimum_degree is None:
            minimum_degree = 4

        # check combination of values of minimum degree and minimum connectivity
        if minimum_connectivity is None:
            if minimum_degree is not None:
                minimum_connectivity = min(3, minimum_degree)
            else:
                minimum_degree, minimum_connectivity = 3, 3
        else:
            if minimum_degree is None:
                minimum_degree = minimum_connectivity
            elif minimum_degree < minimum_connectivity:
                raise ValueError("Minimum connectivity can be at most the minimum degree.")

        # exact connectivity is not implemented for minimum connectivity equal
        # to minimum degree
        if exact_connectivity and minimum_connectivity == minimum_degree:
            raise NotImplementedError("Generation of triangulations with minimum connectivity equal to minimum degree is not implemented.")

        minimum_order = {3: 4, 4: 6, 5: 12}[minimum_degree]

        if order < minimum_order:
            return

        if only_eulerian and order < 6:
            return

        cmd = '-{}m{}c{}{}{} {}'
        command = cmd.format('b' if only_eulerian else '',
                             minimum_degree,
                             minimum_connectivity,
                             'x' if exact_connectivity else '',
                             'd' if dual else '',
                             order)

        yield from graphs.plantri_gen(command)

    def quadrangulations(self, order, minimum_degree=None, minimum_connectivity=None,
                         no_nonfacial_quadrangles=False, dual=False):
        r"""
        An iterator over planar quadrangulations using the plantri generator.

        This uses the plantri generator (see [BM2007]_) which is available
        through the optional package plantri.

        INPUT:

        - ``order`` -- positive integer smaller than or equal to 64;
          this specifies the number of vertices in the generated quadrangulations

        - ``minimum_degree`` -- (default: ``None``) a value `\geq 2` and `\leq
          3`, or ``None``. This specifies the minimum degree of the generated
          quadrangulations. If this is ``None`` and the minimum connectivity is
          specified, then this is set to the same value as the minimum
          connectivity. If the minimum connectivity is also equal to ``None``,
          then this is set to 2.

        - ``minimum_connectivity`` -- (default: ``None``) a value `\geq 2` and
          `\leq 3`, or ``None``. This specifies the minimum connectivity of the
          generated quadrangulations. If this is ``None`` and the option
          ``no_nonfacial_quadrangles`` is set to ``True``, then this is set to
          3. Otherwise if this is ``None`` and the minimum degree is specified,
          then this is set to the minimum degree. If the minimum degree is also
          equal to ``None``, then this is set to 3.

        - ``no_nonfacial_quadrangles`` -- (default: ``False``) if ``True`` only
          quadrangulations with no non-facial quadrangles are generated. This
          option cannot be used if ``minimum_connectivity`` is set to 2.

        - ``dual`` -- (default: ``False``) if ``True`` return instead the
          planar duals of the generated graphs

        OUTPUT:

        An iterator which will produce all planar quadrangulations with the given
        number of vertices as Sage graphs with an embedding set. These will be
        simple graphs (no loops, no multiple edges, no directed edges).

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.set_embedding`,
              :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding` --
              get/set methods for embeddings.

        EXAMPLES:

        The cube is the only 3-connected planar quadrangulation on 8 vertices::

            sage: # optional - plantri
            sage: gen = graphs.quadrangulations(8, minimum_connectivity=3)
            sage: g = next(gen)
            sage: g.is_isomorphic(graphs.CubeGraph(3))
            True
            sage: next(gen)
            Traceback (most recent call last):
            ...
            StopIteration

        An overview of the number of quadrangulations on up to 12 vertices. This
        agrees with :oeis:`A113201`::

            sage: for i in range(4,13):                          # optional - plantri
            ....:     L =  len(list(graphs.quadrangulations(i)))
            ....:     print("{:2d}   {:3d}".format(i,L))
             4     1
             5     1
             6     2
             7     3
             8     9
             9    18
            10    62
            11   198
            12   803

        There are 2 planar quadrangulation on 12 vertices that do not have a
        non-facial quadrangle::

            sage: len([g for g in graphs.quadrangulations(12, no_nonfacial_quadrangles=True)])  # optional - plantri
            2

        Setting ``dual=True`` gives the planar dual graphs::

            sage: [len(g) for g in graphs.quadrangulations(12, no_nonfacial_quadrangles=True, dual=True)]  # optional - plantri
            [10, 10]
        """
        if order < 0:
            raise ValueError("number of vertices should be nonnegative")

        # plantri can only output planar quadrangulations on up to 64 vertices
        if order > 64:
            raise ValueError("number of vertices should be at most 64")

        if minimum_connectivity not in {None, 2, 3}:
            raise ValueError("Minimum connectivity should be None, 2 or 3.")

        if minimum_degree not in {None, 2, 3}:
            raise ValueError("Minimum degree should be None, 2 or 3.")

        if (no_nonfacial_quadrangles and
                minimum_connectivity == 2):
            raise NotImplementedError("Generation of no non-facial quadrangles "
                                      "and minimum connectivity 2 is not implemented")

        # check combination of values of minimum degree and minimum connectivity
        if minimum_connectivity is None:
            if minimum_degree is not None:
                minimum_connectivity = min(2, minimum_degree)
            else:
                minimum_degree, minimum_connectivity = 2, 2
        else:
            if minimum_degree is None:
                minimum_degree = minimum_connectivity
            elif minimum_degree < minimum_connectivity:
                raise ValueError("Minimum connectivity can be at most the minimum degree.")

        minimum_order = {2: 4, 3: 8}[minimum_degree]

        if order < minimum_order:
            return

        if no_nonfacial_quadrangles:
            # for plantri -q the option -c4 means 3-connected with no non-facial quadrangles
            minimum_connectivity = 4

        cmd = '-qm{}c{}{} {}'
        command = cmd.format(minimum_degree,
                             minimum_connectivity,
                             'd' if dual else '',
                             order)

        yield from graphs.plantri_gen(command)

###########################################################################
# Basic Graphs
###########################################################################
    from .generators import basic
    BullGraph = staticmethod(basic.BullGraph)
    ButterflyGraph = staticmethod(basic.ButterflyGraph)
    CircularLadderGraph = staticmethod(basic.CircularLadderGraph)
    ClawGraph = staticmethod(basic.ClawGraph)
    CycleGraph = staticmethod(basic.CycleGraph)
    CompleteGraph = staticmethod(basic.CompleteGraph)
    CompleteBipartiteGraph = staticmethod(basic.CompleteBipartiteGraph)
    CompleteMultipartiteGraph = staticmethod(basic.CompleteMultipartiteGraph)
    CorrelationGraph = staticmethod(basic.CorrelationGraph)
    DiamondGraph = staticmethod(basic.DiamondGraph)
    GemGraph = staticmethod(basic.GemGraph)
    DartGraph = staticmethod(basic.DartGraph)
    ForkGraph = staticmethod(basic.ForkGraph)
    EmptyGraph = staticmethod(basic.EmptyGraph)
    Grid2dGraph = staticmethod(basic.Grid2dGraph)
    GridGraph = staticmethod(basic.GridGraph)
    HouseGraph = staticmethod(basic.HouseGraph)
    HouseXGraph = staticmethod(basic.HouseXGraph)
    LadderGraph = staticmethod(basic.LadderGraph)
    MoebiusLadderGraph = staticmethod(basic.MoebiusLadderGraph)
    PathGraph = staticmethod(basic.PathGraph)
    StarGraph = staticmethod(basic.StarGraph)
    Toroidal6RegularGrid2dGraph = staticmethod(basic.Toroidal6RegularGrid2dGraph)
    ToroidalGrid2dGraph = staticmethod(basic.ToroidalGrid2dGraph)

###########################################################################
# Small Graphs
###########################################################################
    from .generators import smallgraphs, distance_regular
    Balaban10Cage = staticmethod(smallgraphs.Balaban10Cage)
    Balaban11Cage = staticmethod(smallgraphs.Balaban11Cage)
    BidiakisCube = staticmethod(smallgraphs.BidiakisCube)
    BiggsSmithGraph = staticmethod(smallgraphs.BiggsSmithGraph)
    BlanusaFirstSnarkGraph = staticmethod(smallgraphs.BlanusaFirstSnarkGraph)
    BlanusaSecondSnarkGraph = staticmethod(smallgraphs.BlanusaSecondSnarkGraph)
    BrinkmannGraph = staticmethod(smallgraphs.BrinkmannGraph)
    BrouwerHaemersGraph = staticmethod(smallgraphs.BrouwerHaemersGraph)
    BuckyBall = staticmethod(smallgraphs.BuckyBall)
    CameronGraph = staticmethod(smallgraphs.CameronGraph)
    Cell600 = staticmethod(smallgraphs.Cell600)
    Cell120 = staticmethod(smallgraphs.Cell120)
    ChvatalGraph = staticmethod(smallgraphs.ChvatalGraph)
    ClebschGraph = staticmethod(smallgraphs.ClebschGraph)
    cocliques_HoffmannSingleton = staticmethod(distance_regular.cocliques_HoffmannSingleton)
    ConwaySmith_for_3S7 = staticmethod(distance_regular.ConwaySmith_for_3S7)
    CoxeterGraph = staticmethod(smallgraphs.CoxeterGraph)
    CubeplexGraph = staticmethod(smallgraphs.CubeplexGraph)
    DejterGraph = staticmethod(smallgraphs.DejterGraph)
    DesarguesGraph = staticmethod(smallgraphs.DesarguesGraph)
    distance_3_doubly_truncated_Golay_code_graph = staticmethod(distance_regular.distance_3_doubly_truncated_Golay_code_graph)
    DoubleStarSnark = staticmethod(smallgraphs.DoubleStarSnark)
    DoublyTruncatedWittGraph = staticmethod(distance_regular.DoublyTruncatedWittGraph)
    DurerGraph = staticmethod(smallgraphs.DurerGraph)
    DyckGraph = staticmethod(smallgraphs.DyckGraph)
    EllinghamHorton54Graph = staticmethod(smallgraphs.EllinghamHorton54Graph)
    EllinghamHorton78Graph = staticmethod(smallgraphs.EllinghamHorton78Graph)
    ErreraGraph = staticmethod(smallgraphs.ErreraGraph)
    F26AGraph = staticmethod(smallgraphs.F26AGraph)
    FlowerSnark = staticmethod(smallgraphs.FlowerSnark)
    FolkmanGraph = staticmethod(smallgraphs.FolkmanGraph)
    FosterGraph = staticmethod(smallgraphs.FosterGraph)
    FosterGraph3S6 = staticmethod(distance_regular.FosterGraph3S6)
    FranklinGraph = staticmethod(smallgraphs.FranklinGraph)
    FruchtGraph = staticmethod(smallgraphs.FruchtGraph)
    GoldnerHararyGraph = staticmethod(smallgraphs.GoldnerHararyGraph)
    GolombGraph = staticmethod(smallgraphs.GolombGraph)
    GossetGraph = staticmethod(smallgraphs.GossetGraph)
    graph_3O73 = staticmethod(distance_regular.graph_3O73)
    GrayGraph = staticmethod(smallgraphs.GrayGraph)
    GritsenkoGraph = staticmethod(smallgraphs.GritsenkoGraph)
    GrotzschGraph = staticmethod(smallgraphs.GrotzschGraph)
    HallJankoGraph = staticmethod(smallgraphs.HallJankoGraph)
    WellsGraph = staticmethod(smallgraphs.WellsGraph)
    HarborthGraph = staticmethod(smallgraphs.HarborthGraph)
    HarriesGraph = staticmethod(smallgraphs.HarriesGraph)
    HarriesWongGraph = staticmethod(smallgraphs.HarriesWongGraph)
    HeawoodGraph = staticmethod(smallgraphs.HeawoodGraph)
    HerschelGraph = staticmethod(smallgraphs.HerschelGraph)
    HigmanSimsGraph = staticmethod(smallgraphs.HigmanSimsGraph)
    HoffmanGraph = staticmethod(smallgraphs.HoffmanGraph)
    HoffmanSingletonGraph = staticmethod(smallgraphs.HoffmanSingletonGraph)
    HoltGraph = staticmethod(smallgraphs.HoltGraph)
    HortonGraph = staticmethod(smallgraphs.HortonGraph)
    IoninKharaghani765Graph = staticmethod(smallgraphs.IoninKharaghani765Graph)
    IvanovIvanovFaradjevGraph = staticmethod(distance_regular.IvanovIvanovFaradjevGraph)
    J2Graph = staticmethod(distance_regular.J2Graph)
    JankoKharaghaniGraph = staticmethod(smallgraphs.JankoKharaghaniGraph)
    JankoKharaghaniTonchevGraph = staticmethod(smallgraphs.JankoKharaghaniTonchevGraph)
    KittellGraph = staticmethod(smallgraphs.KittellGraph)
    KrackhardtKiteGraph = staticmethod(smallgraphs.KrackhardtKiteGraph)
    Klein3RegularGraph = staticmethod(smallgraphs.Klein3RegularGraph)
    Klein7RegularGraph = staticmethod(smallgraphs.Klein7RegularGraph)
    LargeWittGraph = staticmethod(distance_regular.LargeWittGraph)
    LeonardGraph = staticmethod(distance_regular.LeonardGraph)
    LjubljanaGraph = staticmethod(smallgraphs.LjubljanaGraph)
    vanLintSchrijverGraph = staticmethod(distance_regular.vanLintSchrijverGraph)
    LivingstoneGraph = staticmethod(smallgraphs.LivingstoneGraph)
    locally_GQ42_distance_transitive_graph = staticmethod(distance_regular.locally_GQ42_distance_transitive_graph)
    LocalMcLaughlinGraph = staticmethod(smallgraphs.LocalMcLaughlinGraph)
    M22Graph = staticmethod(smallgraphs.M22Graph)
    MarkstroemGraph = staticmethod(smallgraphs.MarkstroemGraph)
    MathonStronglyRegularGraph = staticmethod(smallgraphs.MathonStronglyRegularGraph)
    McGeeGraph = staticmethod(smallgraphs.McGeeGraph)
    McLaughlinGraph = staticmethod(smallgraphs.McLaughlinGraph)
    MeredithGraph = staticmethod(smallgraphs.MeredithGraph)
    MoebiusKantorGraph = staticmethod(smallgraphs.MoebiusKantorGraph)
    MoserSpindle = staticmethod(smallgraphs.MoserSpindle)
    MurtyGraph = staticmethod(smallgraphs.MurtyGraph)
    NauruGraph = staticmethod(smallgraphs.NauruGraph)
    PappusGraph = staticmethod(smallgraphs.PappusGraph)
    PoussinGraph = staticmethod(smallgraphs.PoussinGraph)
    PerkelGraph = staticmethod(smallgraphs.PerkelGraph)
    PetersenGraph = staticmethod(smallgraphs.PetersenGraph)
    RobertsonGraph = staticmethod(smallgraphs.RobertsonGraph)
    SchlaefliGraph = staticmethod(smallgraphs.SchlaefliGraph)
    shortened_00_11_binary_Golay_code_graph = staticmethod(distance_regular.shortened_00_11_binary_Golay_code_graph)
    shortened_000_111_extended_binary_Golay_code_graph = staticmethod(distance_regular.shortened_000_111_extended_binary_Golay_code_graph)
    ShrikhandeGraph = staticmethod(smallgraphs.ShrikhandeGraph)
    SimsGewirtzGraph = staticmethod(smallgraphs.SimsGewirtzGraph)
    SousselierGraph = staticmethod(smallgraphs.SousselierGraph)
    SylvesterGraph = staticmethod(smallgraphs.SylvesterGraph)
    SzekeresSnarkGraph = staticmethod(smallgraphs.SzekeresSnarkGraph)
    ThomsenGraph = staticmethod(smallgraphs.ThomsenGraph)
    TietzeGraph = staticmethod(smallgraphs.TietzeGraph)
    TricornGraph = staticmethod(smallgraphs.TricornGraph)
    Tutte12Cage = staticmethod(smallgraphs.Tutte12Cage)
    TruncatedIcosidodecahedralGraph = staticmethod(smallgraphs.TruncatedIcosidodecahedralGraph)
    TruncatedTetrahedralGraph = staticmethod(smallgraphs.TruncatedTetrahedralGraph)
    TruncatedWittGraph = staticmethod(distance_regular.TruncatedWittGraph)
    TutteCoxeterGraph = staticmethod(smallgraphs.TutteCoxeterGraph)
    TutteGraph = staticmethod(smallgraphs.TutteGraph)
    TwinplexGraph = staticmethod(smallgraphs.TwinplexGraph)
    U42Graph216 = staticmethod(smallgraphs.U42Graph216)
    U42Graph540 = staticmethod(smallgraphs.U42Graph540)
    WagnerGraph = staticmethod(smallgraphs.WagnerGraph)
    WatkinsSnarkGraph = staticmethod(smallgraphs.WatkinsSnarkGraph)
    WienerArayaGraph = staticmethod(smallgraphs.WienerArayaGraph)
    SuzukiGraph = staticmethod(smallgraphs.SuzukiGraph)

###########################################################################
# Platonic Solids
###########################################################################
    from .generators import platonic_solids
    DodecahedralGraph = staticmethod(platonic_solids.DodecahedralGraph)
    HexahedralGraph = staticmethod(platonic_solids.HexahedralGraph)
    IcosahedralGraph = staticmethod(platonic_solids.IcosahedralGraph)
    OctahedralGraph = staticmethod(platonic_solids.OctahedralGraph)
    TetrahedralGraph = staticmethod(platonic_solids.TetrahedralGraph)

###########################################################################
# Families
###########################################################################
    from . import cographs as cographs_module
    from .generators import families
    from . import strongly_regular_db
    AlternatingFormsGraph = staticmethod(distance_regular.AlternatingFormsGraph)
    AztecDiamondGraph = staticmethod(families.AztecDiamondGraph)
    BalancedTree = staticmethod(families.BalancedTree)
    BarbellGraph = staticmethod(families.BarbellGraph)
    BilinearFormsGraph = staticmethod(distance_regular.BilinearFormsGraph)
    BiwheelGraph = staticmethod(families.BiwheelGraph)
    BubbleSortGraph = staticmethod(families.BubbleSortGraph)
    CaiFurerImmermanGraph = staticmethod(families.CaiFurerImmermanGraph)
    chang_graphs = staticmethod(families.chang_graphs)
    CirculantGraph = staticmethod(families.CirculantGraph)
    cographs = staticmethod(cographs_module.cographs)
    CubeGraph = staticmethod(families.CubeGraph)
    CubeConnectedCycle = staticmethod(families.CubeConnectedCycle)
    DipoleGraph = staticmethod(families.DipoleGraph)
    distance_regular_graph = staticmethod(distance_regular.distance_regular_graph)
    DorogovtsevGoltsevMendesGraph = staticmethod(families.DorogovtsevGoltsevMendesGraph)
    DoubleGeneralizedPetersenGraph = staticmethod(families.DoubleGeneralizedPetersenGraph)
    DoubleGrassmannGraph = staticmethod(distance_regular.DoubleGrassmannGraph)
    DoubleOddGraph = staticmethod(distance_regular.DoubleOddGraph)
    EgawaGraph = staticmethod(families.EgawaGraph)
    FibonacciTree = staticmethod(families.FibonacciTree)
    FoldedCubeGraph = staticmethod(families.FoldedCubeGraph)
    FriendshipGraph = staticmethod(families.FriendshipGraph)
    FurerGadget = staticmethod(families.FurerGadget)
    FuzzyBallGraph = staticmethod(families.FuzzyBallGraph)
    GeneralisedDodecagonGraph = staticmethod(distance_regular.GeneralisedDodecagonGraph)
    GeneralisedHexagonGraph = staticmethod(distance_regular.GeneralisedHexagonGraph)
    GeneralisedOctagonGraph = staticmethod(distance_regular.GeneralisedOctagonGraph)
    GeneralizedPetersenGraph = staticmethod(families.GeneralizedPetersenGraph)
    GeneralizedSierpinskiGraph = staticmethod(families.GeneralizedSierpinskiGraph)
    GoethalsSeidelGraph = staticmethod(families.GoethalsSeidelGraph)
    GrassmannGraph = staticmethod(distance_regular.GrassmannGraph)
    HalfCube = staticmethod(distance_regular.HalfCube)
    HammingGraph = staticmethod(families.HammingGraph)
    HanoiTowerGraph = staticmethod(families.HanoiTowerGraph)
    HararyGraph = staticmethod(families.HararyGraph)
    HermitianFormsGraph = staticmethod(distance_regular.HermitianFormsGraph)
    HyperStarGraph = staticmethod(families.HyperStarGraph)
    IGraph = staticmethod(families.IGraph)
    JohnsonGraph = staticmethod(families.JohnsonGraph)
    KneserGraph = staticmethod(families.KneserGraph)
    LCFGraph = staticmethod(families.LCFGraph)
    line_graph_forbidden_subgraphs = staticmethod(families.line_graph_forbidden_subgraphs)
    LollipopGraph = staticmethod(families.LollipopGraph)
    MathonPseudocyclicMergingGraph = staticmethod(families.MathonPseudocyclicMergingGraph)
    MathonPseudocyclicStronglyRegularGraph = staticmethod(families.MathonPseudocyclicStronglyRegularGraph)
    MuzychukS6Graph = staticmethod(families.MuzychukS6Graph)
    MycielskiGraph = staticmethod(families.MycielskiGraph)
    MycielskiStep = staticmethod(families.MycielskiStep)
    NKStarGraph = staticmethod(families.NKStarGraph)
    NStarGraph = staticmethod(families.NStarGraph)
    OddGraph = staticmethod(families.OddGraph)
    PaleyGraph = staticmethod(families.PaleyGraph)
    PasechnikGraph = staticmethod(families.PasechnikGraph)
    petersen_family = staticmethod(families.petersen_family)
    RingedTree = staticmethod(families.RingedTree)
    RoseWindowGraph = staticmethod(families.RoseWindowGraph)
    SierpinskiGasketGraph = staticmethod(families.SierpinskiGasketGraph)
    SquaredSkewHadamardMatrixGraph = staticmethod(families.SquaredSkewHadamardMatrixGraph)
    SwitchedSquaredSkewHadamardMatrixGraph = staticmethod(families.SwitchedSquaredSkewHadamardMatrixGraph)
    StaircaseGraph = staticmethod(families.StaircaseGraph)
    strongly_regular_graph = staticmethod(strongly_regular_db.strongly_regular_graph)
    TabacjnGraph = staticmethod(families.TabacjnGraph)
    TadpoleGraph = staticmethod(families.TadpoleGraph)
    trees = staticmethod(families.trees)
    TruncatedBiwheelGraph = staticmethod(families.TruncatedBiwheelGraph)
    nauty_gentreeg = staticmethod(families.nauty_gentreeg)
    TuranGraph = staticmethod(families.TuranGraph)
    UstimenkoGraph = staticmethod(distance_regular.UstimenkoGraph)
    WheelGraph = staticmethod(families.WheelGraph)
    WindmillGraph = staticmethod(families.WindmillGraph)

###########################################################################
# Graphs from classical geometries over `F_q`
###########################################################################
    from .generators import classical_geometries
    AffineOrthogonalPolarGraph = staticmethod(classical_geometries.AffineOrthogonalPolarGraph)
    AhrensSzekeresGeneralizedQuadrangleGraph = staticmethod(classical_geometries.AhrensSzekeresGeneralizedQuadrangleGraph)
    NonisotropicOrthogonalPolarGraph = staticmethod(classical_geometries.NonisotropicOrthogonalPolarGraph)
    NonisotropicUnitaryPolarGraph = staticmethod(classical_geometries.NonisotropicUnitaryPolarGraph)
    OrthogonalDualPolarGraph = staticmethod(classical_geometries.OrthogonalDualPolarGraph)
    OrthogonalPolarGraph = staticmethod(classical_geometries.OrthogonalPolarGraph)
    SymplecticDualPolarGraph = staticmethod(classical_geometries.SymplecticDualPolarGraph)
    SymplecticPolarGraph = staticmethod(classical_geometries.SymplecticPolarGraph)
    TaylorTwographDescendantSRG = staticmethod(classical_geometries.TaylorTwographDescendantSRG)
    TaylorTwographSRG = staticmethod(classical_geometries.TaylorTwographSRG)
    T2starGeneralizedQuadrangleGraph = staticmethod(classical_geometries.T2starGeneralizedQuadrangleGraph)
    Nowhere0WordsTwoWeightCodeGraph = staticmethod(classical_geometries.Nowhere0WordsTwoWeightCodeGraph)
    HaemersGraph = staticmethod(classical_geometries.HaemersGraph)
    CossidentePenttilaGraph = staticmethod(classical_geometries.CossidentePenttilaGraph)
    UnitaryDualPolarGraph = staticmethod(classical_geometries.UnitaryDualPolarGraph)
    UnitaryPolarGraph = staticmethod(classical_geometries.UnitaryPolarGraph)

###########################################################################
# Chessboard Graphs
###########################################################################
    from .generators import chessboard
    ChessboardGraphGenerator = staticmethod(chessboard.ChessboardGraphGenerator)
    BishopGraph = staticmethod(chessboard.BishopGraph)
    KingGraph = staticmethod(chessboard.KingGraph)
    KnightGraph = staticmethod(chessboard.KnightGraph)
    QueenGraph = staticmethod(chessboard.QueenGraph)
    RookGraph = staticmethod(chessboard.RookGraph)

###########################################################################
# Intersection graphs
###########################################################################
    from .generators import intersection
    IntervalGraph = staticmethod(intersection.IntervalGraph)
    IntersectionGraph = staticmethod(intersection.IntersectionGraph)
    PermutationGraph = staticmethod(intersection.PermutationGraph)
    OrthogonalArrayBlockGraph = staticmethod(intersection.OrthogonalArrayBlockGraph)
    ToleranceGraph = staticmethod(intersection.ToleranceGraph)

###########################################################################
# Random Graphs
###########################################################################
    from .generators import random
    RandomBarabasiAlbert = staticmethod(random.RandomBarabasiAlbert)
    RandomBipartite = staticmethod(random.RandomBipartite)
    RandomRegularBipartite = staticmethod(random.RandomRegularBipartite)
    RandomBicubicPlanar = staticmethod(random.RandomBicubicPlanar)
    RandomBlockGraph = staticmethod(random.RandomBlockGraph)
    RandomBoundedToleranceGraph = staticmethod(random.RandomBoundedToleranceGraph)
    RandomChordalGraph = staticmethod(random.RandomChordalGraph)
    RandomGNM = staticmethod(random.RandomGNM)
    RandomGNP = staticmethod(random.RandomGNP)
    RandomHolmeKim = staticmethod(random.RandomHolmeKim)
    RandomIntervalGraph = staticmethod(random.RandomIntervalGraph)
    RandomLobster = staticmethod(random.RandomLobster)
    RandomNewmanWattsStrogatz = staticmethod(random.RandomNewmanWattsStrogatz)
    RandomProperIntervalGraph = staticmethod(random.RandomProperIntervalGraph)
    RandomRegular = staticmethod(random.RandomRegular)
    RandomShell = staticmethod(random.RandomShell)
    RandomKTree = staticmethod(random.RandomKTree)
    RandomPartialKTree = staticmethod(random.RandomPartialKTree)
    RandomToleranceGraph = staticmethod(random.RandomToleranceGraph)
    RandomTreePowerlaw = staticmethod(random.RandomTreePowerlaw)
    RandomTree = staticmethod(random.RandomTree)
    RandomTriangulation = staticmethod(random.RandomTriangulation)
    RandomUnitDiskGraph = staticmethod(random.RandomUnitDiskGraph)

###########################################################################
# Maps
###########################################################################
    from .generators import world_map
    WorldMap = staticmethod(world_map.WorldMap)
    EuropeMap = staticmethod(world_map.EuropeMap)
    AfricaMap = staticmethod(world_map.AfricaMap)
    USAMap = staticmethod(world_map.USAMap)

###########################################################################
# Degree Sequence
###########################################################################
    from .generators import degree_sequence
    DegreeSequence = staticmethod(degree_sequence.DegreeSequence)
    DegreeSequenceBipartite = staticmethod(degree_sequence.DegreeSequenceBipartite)
    DegreeSequenceConfigurationModel = staticmethod(degree_sequence.DegreeSequenceConfigurationModel)
    DegreeSequenceTree = staticmethod(degree_sequence.DegreeSequenceTree)
    DegreeSequenceExpected = staticmethod(degree_sequence.DegreeSequenceExpected)


def canaug_traverse_vert(g, aut_gens, max_verts, property, dig=False, loops=False, sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free (di)graphs satisfying a
    given property.

    INPUT:

    - ``g`` -- current position on the tree

    - ``aut_gens`` -- list of generators of Aut(g), in list notation

    - ``max_verts`` -- when to retreat

    - ``property`` -- check before traversing below g

    - ``degree_sequence`` -- specify a degree sequence to try to obtain

    EXAMPLES::

        sage: from sage.graphs.graph_generators import canaug_traverse_vert
        sage: list(canaug_traverse_vert(Graph(), [], 3, lambda x: True))
        [Graph on 0 vertices, ... Graph on 3 vertices]

    The best way to access this function is through the graphs()
    iterator:

    Print graphs on 3 or less vertices.

    ::

        sage: for G in graphs(3, augment='vertices'):
        ....:    print(G)
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
        ....:     print(D)
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
            for i in range(len(children)):
                k = 0
                for j in range(possibilities):
                    if (1 << j) & i:
                        if dig and j >= n:
                            k += (1 << (gen[j - n] + n))
                        else:
                            k += (1 << gen[j])
                while children[k] != -1:
                    k = children[k]
                while children[i] != -1:
                    i = children[i]
                if i != k:
                    # union i & k
                    smaller, larger = sorted([i, k])
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
            z = g.copy(sparse=sparse)
            z.add_vertex(n)
            edges = []
            if dig:
                index = 0
                while 2 * index < possibilities:
                    if (1 << index) & i:
                        edges.append((index, n))
                    index += 1
                while index < possibilities:
                    if (1 << index) & i:
                        edges.append((n, index - n))
                    index += 1
            else:
                index = 0
                while (1 << index) <= i:
                    if (1 << index) & i:
                        edges.append((index, n))
                    index += 1
            z.add_edges(edges)
            z_s = []
            if property(z):
                z_s.append(z)
            if loops:
                z = z.copy(sparse=sparse)
                z.add_edge((n, n))
                if property(z):
                    z_s.append(z)
            for z in z_s:
                z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices(sort=True)], certificate=True, dig=(dig or loops))
                cut_vert = 0
                while canonical_relabeling[cut_vert] != n:
                    cut_vert += 1
                sub_verts = [v for v in z if v != cut_vert]
                m_z = z.subgraph(sub_verts)

                if m_z == g:
                    for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, sparse=sparse):
                        yield a
                else:
                    for possibility in check_aut(z_aut_gens, cut_vert, n):
                        if m_z.relabel(dict(enumerate(possibility)), check_input=False, inplace=False) == g:
                            for a in canaug_traverse_vert(z, z_aut_gens, max_verts, property, dig=dig, loops=loops, sparse=sparse):
                                yield a
                            break


def check_aut(aut_gens, cut_vert, n):
    """
    Helper function for exhaustive generation.

    At the start, check_aut is given a set of generators for the
    automorphism group, aut_gens. We already know we are looking for
    an element of the auto- morphism group that sends cut_vert to n,
    and check_aut generates these for the canaug_traverse function.

    EXAMPLES:

    Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = list(range(n + 1))
    seen_perms = [perm]
    unchecked_perms = [perm]
    while unchecked_perms:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for i in range(len(perm)):
                new_perm[i] = gen[perm[i]]
            if new_perm not in seen_perms:
                seen_perms.append(new_perm)
                unchecked_perms.append(new_perm)
                if new_perm[cut_vert] == n:
                    yield new_perm


def canaug_traverse_edge(g, aut_gens, property, dig=False, loops=False, sparse=True):
    """
    Main function for exhaustive generation. Recursive traversal of a
    canonically generated tree of isomorph free graphs satisfying a
    given property.

    INPUT:

    - ``g`` -- current position on the tree

    - ``aut_gens`` -- list of generators of Aut(g), in list notation

    - ``property`` -- check before traversing below g

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
        ....:     print(G)
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices
        Graph on 3 vertices

    Print digraphs on 3 or less vertices.

    ::

        sage: for G in digraphs(3):
        ....:     print(G)
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
        max_size = n * (n - 1)
    else:
        max_size = (n * (n - 1)) >> 1  # >> 1 is just / 2 (this is n choose 2)
    if loops:
        max_size += n
    if g.size() < max_size:
        # build a list representing C(g) - the edge to be added
        # is one of max_size choices
        if dig:
            children = [[(j, i) for i in range(n)] for j in range(n)]
        else:
            children = [[(j, i) for i in range(j)] for j in range(n)]
        # union-find C(g) under Aut(g)
        orbits = list(range(n))
        for gen in aut_gens:
            for iii in range(n):
                if orbits[gen[iii]] != orbits[iii]:
                    temp = orbits[gen[iii]]
                    for jjj in range(n):
                        if orbits[jjj] == temp:
                            orbits[jjj] = orbits[iii]
                if dig:
                    jjj_range = list(range(iii)) + list(range(iii + 1, n))
                else:
                    jjj_range = list(range(iii))  # iii > jjj
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
                j_range = list(range(i)) + list(range(i + 1, n))
            else:
                j_range = list(range(i))
            for j in j_range:
                if children[i][j] == (i, j):
                    roots.append((i, j))
        if loops:
            seen = []
            for i in range(n):
                if orbits[i] not in seen:
                    roots.append((i, i))
                    seen.append(orbits[i])
        for i, j in roots:
            if g.has_edge(i, j):
                continue
            # construct a z for each edge in roots...
            z = g.copy(sparse=sparse)
            z.add_edge(i, j)
            if not property(z):
                continue
            z_aut_gens, _, canonical_relabeling = search_tree(z, [z.vertices(sort=True)], certificate=True, dig=(dig or loops))
            relabel_inverse = [0]*n
            for ii in range(n):
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
                for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, sparse=sparse):
                    yield a
            else:
                for possibility in check_aut_edge(z_aut_gens, cut_edge, i, j, n, dig=dig):
                    if m_z.relabel(possibility, inplace=False) == g:
                        for a in canaug_traverse_edge(z, z_aut_gens, property, dig=dig, loops=loops, sparse=sparse):
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

    EXAMPLES:

    Note that the last two entries indicate that none of the
    automorphism group has yet been searched - we are starting at the
    identity [0, 1, 2, 3] and so far that is all we have seen. We
    return automorphisms mapping 2 to 3::

        sage: from sage.graphs.graph_generators import check_aut
        sage: list( check_aut( [ [0, 3, 2, 1], [1, 0, 3, 2], [2, 1, 0, 3] ], 2, 3))
        [[1, 0, 3, 2], [1, 2, 3, 0]]
    """
    from copy import copy
    perm = list(range(n))
    seen_perms = [perm]
    unchecked_perms = [perm]
    while unchecked_perms:
        perm = unchecked_perms.pop(0)
        for gen in aut_gens:
            new_perm = copy(perm)
            for ii in range(n):
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
