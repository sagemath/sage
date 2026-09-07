import pytest
import sys

from sage.graphs.generators.distance_regular import graph_from_GQ_spread
from sage.graphs.generators.smallgraphs import _EllipticLinesProjectivePlaneScheme as ES
from sage.graphs.graph_generators import graphs
from sage.misc.randstate import current_randstate
from sage.misc.randstate import set_random_seed

# This is certainly not the right method to get the seed
seed = int(current_randstate().long_seed() % sys.maxsize)
print(f"random seed used for these tests: {seed}")


@pytest.mark.longlong
def test_shortened_000_111_extended_binary_Golay_code_graph():
    r"""
    Test that Sage produces a graph equal to the one that we get
    from this construction.

    The construction itself takes a long time.
    """
    from sage.coding import codes_catalog
    from sage.coding.linear_code import LinearCode
    from sage.graphs.generators.distance_regular import (
        shortened_000_111_extended_binary_Golay_code_graph
    )
    from sage.matrix.constructor import matrix
    from sage.rings.finite_rings.finite_field_constructor import FiniteField

    code = codes_catalog.GolayCode(FiniteField(2))
    C_basis = code.basis()

    # now special shortening
    v = C_basis[0] + C_basis[1] + C_basis[2]  # v has 111 at the start
    C_basis = C_basis[3:]
    C_basis.append(v)
    C_basis = [x[3:] for x in C_basis]

    code = LinearCode(matrix(FiniteField(2), C_basis))
    G = code.cosetGraph()
    G.name("Shortened 000 111 extended binary Golay code")
    assert G.is_distance_regular()

    H = shortened_000_111_extended_binary_Golay_code_graph()
    assert G == H


# ****************************************************************************
# Methods for testing the behavior of parameter immutable
# ****************************************************************************

graph_constructors = [
    # basic.py
    (graphs.BullGraph, (), {}),
    (graphs.ButterflyGraph, (), {}),
    (graphs.CircularLadderGraph, (3,), {}),
    (graphs.ClawGraph, (), {}),
    (graphs.CompleteBipartiteGraph, (2, 3), {'set_position': False}),
    (graphs.CompleteGraph, (4,), {}),
    (graphs.CompleteMultipartiteGraph, ([2, 2, 3],), {}),
    (graphs.CorrelationGraph, ([[1,2,3], [4,5,6], [7,8,9999]], 0.9, False), {}),
    (graphs.CycleGraph, (4,), {}),
    (graphs.DartGraph, (), {}),
    (graphs.DiamondGraph, (), {}),
    (graphs.EmptyGraph, (), {}),
    (graphs.ForkGraph, (), {}),
    (graphs.GemGraph, (), {}),
    (graphs.Grid2dGraph, (2, 3), {'set_positions': False}),
    (graphs.GridGraph, ([2, 2, 3],), {}),
    (graphs.HouseGraph, (), {}),
    (graphs.HouseXGraph, (), {}),
    (graphs.LadderGraph, (3,), {}),
    (graphs.MoebiusLadderGraph, (3,), {}),
    (graphs.PathGraph, (3,), {}),
    (graphs.StarGraph, (3,), {}),
    (graphs.Toroidal6RegularGrid2dGraph, (4, 4), {}),
    (graphs.ToroidalGrid2dGraph, (3, 3), {}),
    # chessboard.py
    (graphs.BishopGraph, ([2, 2],), {}),
    (graphs.KingGraph, ([2, 2],), {}),
    (graphs.KnightGraph, ([2, 2],), {}),
    (graphs.QueenGraph, ([2, 2],), {'relabel': True}),
    (graphs.RookGraph, ([2, 2],), {}),
    # degree_sequence.py
    (graphs.DegreeSequence, ([3,3,3,3],), {}),
    (graphs.DegreeSequenceBipartite, ([1, 1, 2, 3, 3], [2, 4, 4]), {}),
    (graphs.DegreeSequenceConfigurationModel, ([2, 2, 2],), {}),
    (graphs.DegreeSequenceExpected, ([1,2,3,2,3],), {}),
    (graphs.DegreeSequenceTree, ([3,1,3,3,1,1,1,2,1],), {}),
    # distance_regular.py
    (graph_from_GQ_spread, (2, 4), {}),
    (graphs.AlternatingFormsGraph, (2, 2), {}),
    (graphs.BilinearFormsGraph, (1, 1, 2), {}),
    (graphs.ConwaySmith_for_3S7, (), {}),
    (graphs.DoubleGrassmannGraph, (2, 0), {}),
    (graphs.DoubleOddGraph, (1,), {}),
    (graphs.DoublyTruncatedWittGraph, (), {}),
    (graphs.FosterGraph3S6, (), {}),
    (graphs.GeneralisedHexagonGraph, (1, 2), {}),
    (graphs.GeneralisedOctagonGraph, (1, 1), {}),
    (graphs.GrassmannGraph, (2, 3, 1), {}),
    (graphs.HalfCube, (2,), {}),
    (graphs.HermitianFormsGraph, (1, 2), {}),
    (graphs.LargeWittGraph, (), {}),
    (graphs.LeonardGraph, (), {}),
    (graphs.TruncatedWittGraph, (), {}),
    (graphs.UstimenkoGraph, (2, 2), {}),
    (graphs.cocliques_HoffmannSingleton, (), {}),
    (graphs.distance_3_doubly_truncated_Golay_code_graph, (), {}),
    (graphs.shortened_000_111_extended_binary_Golay_code_graph, (), {}),
    (graphs.shortened_00_11_binary_Golay_code_graph, (), {}),
    (graphs.vanLintSchrijverGraph, (), {}),
    # (graphs.GeneralisedDodecagonGraph, (1, 2), {}),           # optional - internet gap_package_atlasrep
    # (graphs.graph_3O73, (), {}),                              # optional - internet gap_package_atlasrep
    # (graphs.IvanovIvanovFaradjevGraph, (), {}),               # optional - internet gap_package_atlasrep
    # (graphs.J2Graph, (), {}),                                 # optional - internet gap_package_atlasrep
    # (graphs.locally_GQ42_distance_transitive_graph, (), {}),  # optional - internet gap_package_atlasrep
    # families.py
    (graphs.AztecDiamondGraph, (2,), {}),
    (graphs.BarbellGraph, (3, 2), {}),
    (graphs.BiwheelGraph, (4,), {}),
    (graphs.BubbleSortGraph, (2,), {}),
    (graphs.CirculantGraph, (5, 2), {}),
    (graphs.CubeConnectedCycle, (3,), {}),
    (graphs.CubeGraph, (2,), {'embedding': 0}),
    (graphs.DipoleGraph, (1,), {}),
    (graphs.DorogovtsevGoltsevMendesGraph, (0,), {}),
    (graphs.DoubleGeneralizedPetersenGraph, (5, 2), {}),
    (graphs.EgawaGraph, (1, 2), {}),
    (graphs.FoldedCubeGraph, (2,), {}),
    (graphs.FriendshipGraph, (2,), {}),
    (graphs.FuzzyBallGraph, ([3, 1], 2), {}),
    (graphs.GeneralizedPetersenGraph, (5, 2), {}),
    (graphs.GeneralizedSierpinskiGraph, (graphs.CycleGraph(3), 2), {'stretch': 1}),
    (graphs.GoethalsSeidelGraph, (2, 3), {}),
    (graphs.HammingGraph, (1, 2), {}),
    (graphs.HanoiTowerGraph, (2, 3), {'labels': False, 'positions': False}),
    (graphs.HararyGraph, (3, 5), {}),
    (graphs.HyperStarGraph, (3, 2), {}),
    (graphs.IGraph, (5, 1, 2), {}),
    (graphs.JohnsonGraph, (5, 2), {}),
    (graphs.KneserGraph, (5, 2), {}),
    (graphs.LCFGraph, (4, [2,-2], 2), {}),
    (graphs.LollipopGraph, (3, 3), {}),
    (graphs.MathonPseudocyclicMergingGraph, (ES(3), 0), {}),
    (graphs.MathonPseudocyclicStronglyRegularGraph, (1,), {}),
    (graphs.MuzychukS6Graph, (4, 2), {}),
    (graphs.MycielskiGraph, (), {'k': 3, 'relabel': False}),
    (graphs.MycielskiStep, (graphs.CycleGraph(3),), {}),
    (graphs.NKStarGraph, (2, 1), {}),
    (graphs.NStarGraph, (2,), {}),
    (graphs.OddGraph, (2,), {}),
    (graphs.PaleyGraph, (5,), {}),
    (graphs.PasechnikGraph, (2,), {}),
    (graphs.RingedTree, (2,), {'vertex_labels': False}),
    (graphs.RingedTree, (2,), {'vertex_labels': True}),
    (graphs.RoseWindowGraph, (3, 2, 1), {}),
    (graphs.SierpinskiGasketGraph, (1,), {}),
    (graphs.SquaredSkewHadamardMatrixGraph, (1,), {}),
    (graphs.StaircaseGraph, (3,), {}),
    (graphs.SwitchedSquaredSkewHadamardMatrixGraph, (1,), {}),
    (graphs.TabacjnGraph, (3, 1, 2, 1), {}),
    (graphs.TadpoleGraph, (3, 1), {}),
    (graphs.TruncatedBiwheelGraph, (3,), {}),
    (graphs.TuranGraph, (2, 1), {}),
    (graphs.WheelGraph, (4,), {}),
    (graphs.WindmillGraph, (3, 3), {}),
    # intersection.py
    (graphs.IntersectionGraph, ([(1, 2), (2, 3)],), {}),
    (graphs.IntervalGraph, ([(1, 2), (1, 2), (2, 3)],), {}),
    (graphs.OrthogonalArrayBlockGraph, (2, 2), {}),
    (graphs.PermutationGraph, ([3, 4, 5, 1, 2],), {}),
    (graphs.ToleranceGraph, ([(1, 4, 3), (1, 2, 1)], ), {}),
    # platonic_solids.py
    (graphs.DodecahedralGraph, (), {}),
    (graphs.HexahedralGraph, (), {}),
    (graphs.IcosahedralGraph, (), {}),
    (graphs.OctahedralGraph, (), {}),
    (graphs.TetrahedralGraph, (), {}),
    # random.py
    (graphs.RandomBarabasiAlbert, (6, 2), {}),
    (graphs.RandomBicubicPlanar, (6,), {}),
    (graphs.RandomBipartite, (5, 2, .5), {}),
    (graphs.RandomBlockGraph, (6, 4), {}),
    (graphs.RandomBoundedToleranceGraph, (6,), {}),
    (graphs.RandomChordalGraph, (5,), {}),
    (graphs.RandomGNM, (3, 2), {'dense': False}),
    (graphs.RandomGNM, (3, 2), {'dense': True}),
    (graphs.RandomGNP, (6, .4), {'algorithm': 'Sage'}),
    (graphs.RandomGNP, (6, .4), {'algorithm': 'networkx'}),
    (graphs.RandomHolmeKim, (6, 2, .3), {}),
    (graphs.RandomIntervalGraph, (5,), {}),
    (graphs.RandomKTree, (8, 3), {}),
    (graphs.RandomNewmanWattsStrogatz, (7, 2, .2), {}),
    (graphs.RandomPartialKTree, (5, 2, 1), {}),
    (graphs.RandomProperIntervalGraph, (5,), {}),
    (graphs.RandomRegular, (3, 8), {}),
    (graphs.RandomRegularBipartite, (1, 2, 2), {}),
    (graphs.RandomShell, ([(10, 20, 0.8), (20, 40, 0.8)], ), {}),
    (graphs.RandomToleranceGraph, (8,), {}),
    (graphs.RandomTriangulation, (6,), {}),
    (graphs.RandomUnitDiskGraph, (6,), {}),
    # smallgraphs.py
    (graphs.Balaban10Cage, (), {'embedding': 2}),
    (graphs.Balaban11Cage, (), {'embedding': 3}),
    (graphs.BidiakisCube, (), {}),
    (graphs.BiggsSmithGraph, (), {'embedding': 2}),
    (graphs.BlanusaFirstSnarkGraph, (), {}),
    (graphs.BlanusaSecondSnarkGraph, (), {}),
    (graphs.BrinkmannGraph, (), {}),
    (graphs.BrouwerHaemersGraph, (), {}),
    (graphs.BuckyBall, (), {}),
    (graphs.CameronGraph, (), {}),
    (graphs.Cell120, (), {}),
    (graphs.Cell600, (), {'embedding': 2}),
    (graphs.ChvatalGraph, (), {}),
    (graphs.ClebschGraph, (), {}),
    (graphs.CoxeterGraph, (), {}),
    (graphs.CubeplexGraph, (), {'embedding': 'NT'}),
    (graphs.DejterGraph, (), {}),
    (graphs.DesarguesGraph, (), {}),
    (graphs.DoubleStarSnark, (), {}),
    (graphs.DurerGraph, (), {}),
    (graphs.DyckGraph, (), {}),
    (graphs.EllinghamHorton54Graph, (), {}),
    (graphs.EllinghamHorton78Graph, (), {}),
    (graphs.ErreraGraph, (), {}),
    (graphs.F26AGraph, (), {}),
    (graphs.FlowerSnark, (), {}),
    (graphs.FolkmanGraph, (), {}),
    (graphs.FosterGraph, (), {}),
    (graphs.FranklinGraph, (), {}),
    (graphs.FruchtGraph, (), {}),
    (graphs.GoldnerHararyGraph, (), {}),
    (graphs.GolombGraph, (), {}),
    (graphs.GossetGraph, (), {}),
    (graphs.GrayGraph, (), {'embedding': 2}),
    (graphs.GritsenkoGraph, (), {}),
    (graphs.GrotzschGraph, (), {}),
    (graphs.HallJankoGraph, (), {'from_string': True}),
    (graphs.HarborthGraph, (), {}),
    (graphs.HarriesGraph, (), {'embedding': 2}),
    (graphs.HarriesWongGraph, (), {'embedding': 2}),
    (graphs.HeawoodGraph, (), {}),
    (graphs.HerschelGraph, (), {}),
    (graphs.HigmanSimsGraph, (), {'relabel': False}),
    (graphs.HoffmanGraph, (), {}),
    (graphs.HoffmanSingletonGraph, (), {}),
    (graphs.HoltGraph, (), {}),
    (graphs.HortonGraph, (), {}),
    (graphs.IoninKharaghani765Graph, (), {}),
    (graphs.JankoKharaghaniGraph, (936,), {}),  # long time
    (graphs.JankoKharaghaniTonchevGraph, (), {}),
    (graphs.KittellGraph, (), {}),
    (graphs.Klein3RegularGraph, (), {}),
    (graphs.Klein7RegularGraph, (), {}),
    (graphs.KrackhardtKiteGraph, (), {}),
    (graphs.LjubljanaGraph, (), {'embedding': 2}),
    (graphs.M22Graph, (), {}),
    (graphs.MarkstroemGraph, (), {}),
    (graphs.MathonStronglyRegularGraph, (0,), {}),
    (graphs.McGeeGraph, (), {'embedding': 1}),
    (graphs.MeredithGraph, (), {}),
    (graphs.MoebiusKantorGraph, (), {}),
    (graphs.MoserSpindle, (), {}),
    (graphs.MurtyGraph, (), {}),
    (graphs.NauruGraph, (), {'embedding': 2}),
    (graphs.PappusGraph, (), {}),
    (graphs.PerkelGraph, (), {}),
    (graphs.PetersenGraph, (), {}),
    (graphs.PoussinGraph, (), {}),
    (graphs.RobertsonGraph, (), {}),
    (graphs.SchlaefliGraph, (), {}),
    (graphs.ShrikhandeGraph, (), {}),
    (graphs.SimsGewirtzGraph, (), {}),
    (graphs.SousselierGraph, (), {}),
    (graphs.SylvesterGraph, (), {}),
    (graphs.SzekeresSnarkGraph, (), {}),
    (graphs.ThomsenGraph, (), {}),
    (graphs.TietzeGraph, (), {}),
    (graphs.TricornGraph, (), {}),
    (graphs.TruncatedTetrahedralGraph, (), {}),
    (graphs.Tutte12Cage, (), {}),
    (graphs.TutteCoxeterGraph, (), {'embedding': 2}),
    (graphs.TutteGraph, (), {}),
    (graphs.TwinplexGraph, (), {'embedding': 'NT'}),
    (graphs.WagnerGraph, (), {}),
    (graphs.WatkinsSnarkGraph, (), {}),
    (graphs.WellsGraph, (), {}),
    (graphs.WienerArayaGraph, (), {}),
    # (graphs.LivingstoneGraph, (), {}),      # optional - internet # not tested
    # (graphs.LocalMcLaughlinGraph, (), {}),  # optional - gap_package_design
    # (graphs.McLaughlinGraph, (), {}),       # optional - gap_package_design
    # (graphs.SuzukiGraph, (), {}),           # optional internet # not tested
    # (graphs.U42Graph216, (), {}),           # optional - gap_package
    # (graphs.U42Graph540, (), {}),           # optional - gap_package
    # Can not be constructed currently, due to numerical issues
    # (graphs.TruncatedIcosidodecahedralGraph, (), {}),
    # trees.pyx
    (graphs.BalancedTree, (2, 1), {}),
    (graphs.Caterpillar, ([0, 1, 1, 0], ), {}),
    (graphs.FibonacciTree, (3,), {}),
    (graphs.RandomLobster, (9, .6, .3), {}),
    (graphs.RandomTree, (5,), {}),
    (graphs.RandomTreePowerlaw, (10,), {}),
    # world_maps.py
    (graphs.AfricaMap, (), {'continental': False}),
    (graphs.AfricaMap, (), {'continental': True}),
    (graphs.EuropeMap, (), {'continental': False}),
    (graphs.EuropeMap, (), {'continental': True}),
    (graphs.USAMap, (), {'continental': False}),
    (graphs.USAMap, (), {'continental': True}),
    (graphs.WorldMap, (), {}),
    ]


def _compare_graphs(mu, im):
    r"""
    Helper method to test the behavior of parameter ``immutable``.

    The tests are robust to vertex labels of different types.

    INPUT:

    - ``mu`` -- a mutable graph

    - ``im`` -- an immutable graph

    EXAMPLES::

        sage: from sage.graphs.generators.generators_test import _compare_graphs
        random seed used for these tests: ...
        sage: mu = graphs.CycleGraph(3, immutable=False)
        sage: im = graphs.CycleGraph(3, immutable=True)
        sage: _compare_graphs(mu, im)
        sage: _compare_graphs(mu, mu)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: _compare_graphs(im, im)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: im = graphs.CycleGraph(4, immutable=True)
        sage: _compare_graphs(mu, im)
        Traceback (most recent call last):
        ...
        AssertionError
    """
    assert mu.is_immutable() is False
    assert im.is_immutable() is True
    assert mu.order() == im.order()
    assert mu.size() == im.size()
    assert set(mu) == set(im)
    edges_mu = {(frozenset((u, v)), label) for u, v, label in mu.edges()}
    edges_im = {(frozenset((u, v)), label) for u, v, label in im.edges()}
    assert edges_mu == edges_im


@pytest.mark.parametrize("graph_constructor, args, kwds", graph_constructors)
def test_parameter_immutable(graph_constructor, args, kwds):
    r"""
    Check the behavior of parameter ``immutable`` on input graph generator.

    Check that the graphs returned by the generator when parameter ``immutable``
    is set to ``True`` or ``False`` have the same sets of vertices and edges.

    This method is applied on all entries defined in ``graph_constructors``.

    INPUT:

    - ``generator`` -- a graph generator

    - ``args`` -- ordered list of the arguments of the generator without default
      values

    - ``kwds`` -- dictionary mapping parameters with default values to the
      desired value
    """
    set_random_seed(seed)
    Gmu = graph_constructor(*args, **kwds, immutable=False)
    set_random_seed(seed)
    Gim = graph_constructor(*args, **kwds, immutable=True)
    _compare_graphs(Gmu, Gim)


def test_parameter_immutable_other():
    r"""
    Check the behavior of parameter ``immutable`` on other graph generators.

    This method is devoted to graph generators returning a pair of items, a list
    of graphs or an iterator.
    """
    def comp1(graph_constructor, args):
        set_random_seed(seed)
        Gmu = graph_constructor(*args, immutable=False)[0]
        set_random_seed(seed)
        Gim = graph_constructor(*args, immutable=True)[0]
        _compare_graphs(Gmu, Gim)

    comp1(graphs.CaiFurerImmermanGraph, (graphs.CycleGraph(3),))
    comp1(graphs.FurerGadget, (1,))

    # Generators returning lists of graphs
    gens = (graphs.chang_graphs, graphs.line_graph_forbidden_subgraphs,
            graphs.p2_forbidden_minors, graphs.petersen_family)
    for gen in gens:
        for Gmu, Gim in zip(gen(immutable=False), gen(immutable=True)):
            _compare_graphs(Gmu, Gim)

    # Check that iterators return trees in the same order
    gen_mu = graphs.trees(4, immutable=False)
    gen_im = graphs.trees(4, immutable=True)
    for Gmu, Gim in zip(gen_mu, gen_im):
        _compare_graphs(Gmu, Gim)

    gen_mu = graphs.nauty_gentreeg("4", immutable=False)
    gen_im = graphs.nauty_gentreeg("4", immutable=True)
    for Gmu, Gim in zip(gen_mu, gen_im):
        _compare_graphs(Gmu, Gim)
