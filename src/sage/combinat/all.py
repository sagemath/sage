r"""
Combinatorics

Introduction
------------

- :ref:`sage.combinat.quickref`
- :ref:`sage.combinat.tutorial`

Topics
------

- :ref:`sage.combinat.counting`
- :ref:`sage.combinat.enumerated_sets`
- :ref:`sage.combinat.catalog_partitions`
- :ref:`sage.combinat.finite_state_machine`
- :ref:`sage.combinat.posets.all`
- :ref:`sage.combinat.species.all`
- :ref:`sage.combinat.designs.all`
- :ref:`sage.combinat.words.all`
- :ref:`sage.combinat.bijectionist`
- :ref:`sage.combinat.chas.all`
- :ref:`sage.combinat.cluster_algebra_quiver.all`
- :ref:`sage.combinat.crystals.all`
- :ref:`sage.combinat.root_system.all`
- :ref:`sage.combinat.sf.all`
- :ref:`sage.combinat.fully_commutative_elements`

Utilities
---------

- :ref:`sage.combinat.output`
- :ref:`sage.combinat.ranker`
- :ref:`sage.combinat.combinatorial_map`
- :ref:`sage.combinat.misc`
"""
from sage.misc.namespace_package import install_doc, install_dict
# install the docstring of this module to the containing package
install_doc(__package__, __doc__)

# install modules quickref and tutorial to the containing package
from sage.combinat import quickref, tutorial
install_dict(__package__, {'quickref': quickref, 'tutorial': tutorial})
del quickref, tutorial

from sage.misc.lazy_import import lazy_import

from sage.combinat.combinat import (CombinatorialObject,
                       bell_number, bell_polynomial, bernoulli_polynomial,
                       catalan_number, euler_number,
                       fibonacci, fibonacci_sequence, fibonacci_xrange,
                       lucas_number1, lucas_number2,
                       number_of_tuples, number_of_unordered_tuples,
                       polygonal_number, stirling_number1, stirling_number2,
                       tuples, unordered_tuples)

from sage.combinat.expnums import expnums

from sage.combinat.chas.all import *
from sage.combinat.crystals.all import *
from sage.combinat.rigged_configurations.all import *

from sage.combinat.dlx import DLXMatrix, AllExactCovers, OneExactCover

# block designs, etc.
from sage.combinat.designs.all import *

# Free modules and friends
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.debruijn_sequence import DeBruijnSequences

from sage.combinat.schubert_polynomial import SchubertPolynomialRing
lazy_import('sage.combinat.key_polynomial', 'KeyPolynomialBasis', as_='KeyPolynomials')
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra, HeckeAlgebraSymmetricGroupT
from sage.combinat.symmetric_group_representations import SymmetricGroupRepresentation, SymmetricGroupRepresentations
from sage.combinat.yang_baxter_graph import YangBaxterGraph

# Permutations
from sage.combinat.permutation import Permutation, Permutations, Arrangements, CyclicPermutations, CyclicPermutationsOfPartition
from sage.combinat.affine_permutation import AffinePermutationGroup
lazy_import('sage.combinat.colored_permutations', ['ColoredPermutations',
                                                   'SignedPermutation',
                                                   'SignedPermutations'])
from sage.combinat.derangements import Derangements
lazy_import('sage.combinat.baxter_permutations', ['BaxterPermutations'])

# RSK
from sage.combinat.rsk import RSK, RSK_inverse, robinson_schensted_knuth, robinson_schensted_knuth_inverse, InsertionRules

# HillmanGrassl
lazy_import("sage.combinat.hillman_grassl", ["WeakReversePlanePartition", "WeakReversePlanePartitions"])

# PerfectMatchings
from sage.combinat.perfect_matching import PerfectMatching, PerfectMatchings

# Integer lists
from sage.combinat.integer_lists import IntegerListsLex

# Compositions
from sage.combinat.composition import Composition, Compositions
from sage.combinat.composition_signed import SignedCompositions

# Partitions
from sage.combinat.partition import (Partition, Partitions, PartitionsInBox,
                        OrderedPartitions, PartitionsGreatestLE,
                        PartitionsGreatestEQ, number_of_partitions)

lazy_import('sage.combinat.partition_tuple', ['PartitionTuple', 'PartitionTuples'])
lazy_import('sage.combinat.partition_kleshchev', ['KleshchevPartitions'])
lazy_import('sage.combinat.skew_partition', ['SkewPartition', 'SkewPartitions'])

# Partition algebra
from sage.combinat.partition_algebra import SetPartitionsAk, SetPartitionsPk, SetPartitionsTk, SetPartitionsIk, SetPartitionsBk, SetPartitionsSk, SetPartitionsRk, SetPartitionsPRk

# Raising operators
lazy_import('sage.combinat.partition_shifting_algebras', 'ShiftingOperatorAlgebra')

# Diagram algebra
from sage.combinat.diagram_algebras import PartitionAlgebra, BrauerAlgebra, TemperleyLiebAlgebra, PlanarAlgebra, PropagatingIdeal

# Descent algebra
lazy_import('sage.combinat.descent_algebra', 'DescentAlgebra')

# Vector Partitions
lazy_import('sage.combinat.vector_partition',
            ['VectorPartition', 'VectorPartitions'])

# Similarity class types
from sage.combinat.similarity_class_type import PrimarySimilarityClassType, PrimarySimilarityClassTypes, SimilarityClassType, SimilarityClassTypes

# Cores
from sage.combinat.core import Core, Cores

# Tableaux
lazy_import('sage.combinat.tableau',
            ["Tableau", "SemistandardTableau", "StandardTableau", "RowStandardTableau", "IncreasingTableau",
             "Tableaux", "SemistandardTableaux", "StandardTableaux", "RowStandardTableaux", "IncreasingTableaux"])
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux, StandardSkewTableaux, SemistandardSkewTableaux
from sage.combinat.ribbon_shaped_tableau import RibbonShapedTableau, RibbonShapedTableaux, StandardRibbonShapedTableaux
from sage.combinat.ribbon_tableau import RibbonTableaux, RibbonTableau, MultiSkewTableaux, MultiSkewTableau, SemistandardMultiSkewTableaux
from sage.combinat.composition_tableau import CompositionTableau, CompositionTableaux

lazy_import('sage.combinat.tableau_tuple',
            ['TableauTuple', 'StandardTableauTuple', 'RowStandardTableauTuple',
             'TableauTuples', 'StandardTableauTuples', 'RowStandardTableauTuples'])
from sage.combinat.k_tableau import WeakTableau, WeakTableaux, StrongTableau, StrongTableaux
lazy_import('sage.combinat.lr_tableau', ['LittlewoodRichardsonTableau',
                                         'LittlewoodRichardsonTableaux'])
lazy_import('sage.combinat.shifted_primed_tableau', ['ShiftedPrimedTableaux',
                                                     'ShiftedPrimedTableau'])

# SuperTableaux
lazy_import('sage.combinat.super_tableau',
            ["StandardSuperTableau", "SemistandardSuperTableau", "StandardSuperTableaux", "SemistandardSuperTableaux"])

# Words
from sage.combinat.words.all import *

lazy_import('sage.combinat.subword', 'Subwords')

from sage.combinat.graph_path import GraphPaths

# Tuples
from sage.combinat.tuple import Tuples, UnorderedTuples

# Alternating sign matrices
lazy_import('sage.combinat.alternating_sign_matrix', ('AlternatingSignMatrix',
                                                      'AlternatingSignMatrices',
                                                      'MonotoneTriangles',
                                                      'ContreTableaux',
                                                      'TruncatedStaircases'))

# Decorated Permutations
lazy_import('sage.combinat.decorated_permutation', ('DecoratedPermutation',
                                                    'DecoratedPermutations'))

# Plane Partitions
lazy_import('sage.combinat.plane_partition', ('PlanePartition',
                                              'PlanePartitions'))

# Parking Functions
lazy_import('sage.combinat.non_decreasing_parking_function',
            ['NonDecreasingParkingFunctions', 'NonDecreasingParkingFunction'])
lazy_import('sage.combinat.parking_functions',
            ['ParkingFunctions', 'ParkingFunction'])

# Trees and Tamari interval posets
from sage.combinat.ordered_tree import (OrderedTree, OrderedTrees,
                          LabelledOrderedTree, LabelledOrderedTrees)
from sage.combinat.binary_tree import (BinaryTree, BinaryTrees,
                         LabelledBinaryTree, LabelledBinaryTrees)
lazy_import('sage.combinat.interval_posets', ['TamariIntervalPoset', 'TamariIntervalPosets'])
lazy_import('sage.combinat.rooted_tree', ('RootedTree', 'RootedTrees',
                         'LabelledRootedTree', 'LabelledRootedTrees'))

from sage.combinat.combination import Combinations

from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.combinat.set_partition_ordered import OrderedSetPartition, OrderedSetPartitions
lazy_import('sage.combinat.multiset_partition_into_sets_ordered',
            ['OrderedMultisetPartitionIntoSets',
             'OrderedMultisetPartitionsIntoSets'])
from sage.combinat.subset import Subsets, subsets, powerset, uniq
from sage.combinat.necklace import Necklaces
lazy_import('sage.combinat.dyck_word', ('DyckWords', 'DyckWord'))
lazy_import('sage.combinat.nu_dyck_word', ('NuDyckWords', 'NuDyckWord'))
from sage.combinat.sloane_functions import sloane
lazy_import('sage.combinat.superpartition', ('SuperPartition',
                                             'SuperPartitions'))

lazy_import('sage.combinat.parallelogram_polyomino',
            ['ParallelogramPolyomino', 'ParallelogramPolyominoes'])

from sage.combinat.root_system.all import *
from sage.combinat.sf.all import *
from sage.combinat.ncsf_qsym.all import *
from sage.combinat.ncsym.all import *
lazy_import('sage.combinat.fqsym', 'FreeQuasisymmetricFunctions')
from sage.combinat.matrices.all import *
# Posets
from sage.combinat.posets.all import *

# Cluster Algebras and Quivers
from sage.combinat.cluster_algebra_quiver.all import *

from sage.combinat import ranker

from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.combinat.integer_vectors_mod_permgroup import IntegerVectorsModPermutationGroup

lazy_import('sage.combinat.q_analogues', ['gaussian_binomial', 'q_binomial', 'number_of_irreducible_polynomials'])

from sage.combinat.species.all import *

lazy_import('sage.combinat.kazhdan_lusztig', 'KazhdanLusztigPolynomial')

lazy_import('sage.combinat.degree_sequences', 'DegreeSequences')

lazy_import('sage.combinat.cyclic_sieving_phenomenon',
            ['CyclicSievingPolynomial', 'CyclicSievingCheck'])

lazy_import('sage.combinat.sidon_sets', 'sidon_sets')

# Puzzles
lazy_import('sage.combinat.knutson_tao_puzzles', 'KnutsonTaoPuzzleSolver')

# Gelfand-Tsetlin patterns
lazy_import('sage.combinat.gelfand_tsetlin_patterns',
            ['GelfandTsetlinPattern', 'GelfandTsetlinPatterns'])

# Finite State Machines (Automaton, Transducer)
lazy_import('sage.combinat.finite_state_machine',
            ['Automaton', 'Transducer', 'FiniteStateMachine'])
lazy_import('sage.combinat.finite_state_machine_generators',
            ['automata', 'transducers'])

# Sequences
lazy_import('sage.combinat.binary_recurrence_sequences',
            'BinaryRecurrenceSequence')
lazy_import('sage.combinat.recognizable_series', 'RecognizableSeriesSpace')
lazy_import('sage.combinat.regular_sequence', 'RegularSequenceRing')

# Six Vertex Model
lazy_import('sage.combinat.six_vertex_model', 'SixVertexModel')

# sine-Gordon Y-systems
lazy_import('sage.combinat.sine_gordon', 'SineGordonYsystem')

# Fully Packed Loop
lazy_import('sage.combinat.fully_packed_loop', ['FullyPackedLoop', 'FullyPackedLoops'])

# Subword complex and cluster complex
lazy_import('sage.combinat.subword_complex', 'SubwordComplex')
lazy_import("sage.combinat.cluster_complex", "ClusterComplex")

# Constellations
lazy_import('sage.combinat.constellation', ['Constellation', 'Constellations'])

# Growth diagrams
lazy_import('sage.combinat.growth', 'GrowthDiagram')

# Path Tableaux
lazy_import('sage.combinat.path_tableaux', 'catalog', as_='path_tableaux')

# Bijectionist
lazy_import('sage.combinat.bijectionist', 'Bijectionist')
