r"""
Combinatorics

Introductory material
---------------------

- :ref:`sage.combinat.quickref`
- :ref:`sage.combinat.tutorial`

Thematic indexes
----------------

- :ref:`sage.combinat.algebraic_combinatorics`

  - :ref:`sage.combinat.chas.all`
  - :ref:`sage.combinat.cluster_algebra_quiver.all`
  - :ref:`sage.combinat.crystals.all`
  - :ref:`sage.combinat.root_system.all`
  - :ref:`sage.combinat.sf.all`
  - :class:`~sage.combinat.fully_commutative_elements.FullyCommutativeElements`

- :ref:`sage.combinat.counting`
- :ref:`sage.combinat.enumerated_sets`
- :ref:`sage.combinat.catalog_partitions`
- :ref:`sage.combinat.finite_state_machine`
- :ref:`sage.combinat.species.all`
- :ref:`sage.combinat.designs.all`
- :ref:`sage.combinat.posets.all`
- :ref:`sage.combinat.words`

Utilities
---------

- :ref:`sage.combinat.output`
- :ref:`sage.combinat.ranker`
- :func:`Combinatorial maps <sage.combinat.combinatorial_map.combinatorial_map>`
- :ref:`sage.combinat.misc`

Related topics
--------------

- :ref:`sage.coding`
- :ref:`sage.dynamics`
- :ref:`sage.graphs`

"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from . import quickref, tutorial

from sage.misc.lazy_import import lazy_import

from .combinat import bell_number, catalan_number, euler_number, fibonacci, \
        lucas_number1, lucas_number2, stirling_number1, stirling_number2, \
        polygonal_number, CombinatorialObject, CombinatorialClass, \
        MapCombinatorialClass, \
        tuples, number_of_tuples, unordered_tuples, number_of_unordered_tuples, \
        bell_polynomial, fibonacci_sequence, fibonacci_xrange, bernoulli_polynomial

lazy_import('sage.combinat.combinat',
            ('InfiniteAbstractCombinatorialClass', 'UnionCombinatorialClass',
             'FilteredCombinatorialClass'),
            deprecation=(31545, 'this class is deprecated, do not use'))


from .expnums import expnums

from sage.combinat.chas.all import *
from sage.combinat.crystals.all import *
from .rigged_configurations.all import *

from sage.combinat.dlx import DLXMatrix, AllExactCovers, OneExactCover

# block designs, etc
from sage.combinat.designs.all import *

# Free modules and friends
from .free_module import CombinatorialFreeModule
from .debruijn_sequence import DeBruijnSequences

from .schubert_polynomial import SchubertPolynomialRing
from .symmetric_group_algebra import SymmetricGroupAlgebra, HeckeAlgebraSymmetricGroupT
from .symmetric_group_representations import SymmetricGroupRepresentation, SymmetricGroupRepresentations
from .yang_baxter_graph import YangBaxterGraph
#from hall_littlewood import HallLittlewood_qp, HallLittlewood_q, HallLittlewood_p

#Permutations
from .permutation import Permutation, Permutations, Arrangements, CyclicPermutations, CyclicPermutationsOfPartition
from .affine_permutation import AffinePermutationGroup
lazy_import('sage.combinat.colored_permutations', ['ColoredPermutations',
                                                   'SignedPermutations'])
from .derangements import Derangements
lazy_import('sage.combinat.baxter_permutations', ['BaxterPermutations'])

#RSK
from .rsk import RSK, RSK_inverse, robinson_schensted_knuth, robinson_schensted_knuth_inverse, InsertionRules

#HillmanGrassl
lazy_import("sage.combinat.hillman_grassl", ["WeakReversePlanePartition", "WeakReversePlanePartitions"])

#PerfectMatchings
from .perfect_matching import PerfectMatching, PerfectMatchings

# Integer lists
from .integer_lists import IntegerListsLex

#Compositions
from .composition import Composition, Compositions
from .composition_signed import SignedCompositions

#Partitions
from .partition import Partition, Partitions, PartitionsInBox,\
     OrderedPartitions, PartitionsGreatestLE, PartitionsGreatestEQ,\
     number_of_partitions

lazy_import('sage.combinat.partition_tuple', ['PartitionTuple', 'PartitionTuples'])
lazy_import('sage.combinat.partition_kleshchev', ['KleshchevPartitions'])
lazy_import('sage.combinat.skew_partition', ['SkewPartition', 'SkewPartitions'])

#Partition algebra
from .partition_algebra import SetPartitionsAk, SetPartitionsPk, SetPartitionsTk, SetPartitionsIk, SetPartitionsBk, SetPartitionsSk, SetPartitionsRk, SetPartitionsPRk

#Raising operators
lazy_import('sage.combinat.partition_shifting_algebras', 'ShiftingOperatorAlgebra')

#Diagram algebra
from .diagram_algebras import PartitionAlgebra, BrauerAlgebra, TemperleyLiebAlgebra, PlanarAlgebra, PropagatingIdeal

#Descent algebra
lazy_import('sage.combinat.descent_algebra', 'DescentAlgebra')

#Vector Partitions
lazy_import('sage.combinat.vector_partition',
            ['VectorPartition', 'VectorPartitions'])

#Similarity class types
from .similarity_class_type import PrimarySimilarityClassType, PrimarySimilarityClassTypes, SimilarityClassType, SimilarityClassTypes

#Cores
from .core import Core, Cores

#Tableaux
lazy_import('sage.combinat.tableau',["Tableau", "SemistandardTableau", "StandardTableau", "RowStandardTableau", "IncreasingTableau",
                                     "Tableaux","SemistandardTableaux","StandardTableaux","RowStandardTableaux", "IncreasingTableaux"])
from .skew_tableau import SkewTableau, SkewTableaux, StandardSkewTableaux, SemistandardSkewTableaux
from .ribbon_shaped_tableau import RibbonShapedTableau, RibbonShapedTableaux, StandardRibbonShapedTableaux
from .ribbon_tableau import RibbonTableaux, RibbonTableau, MultiSkewTableaux, MultiSkewTableau, SemistandardMultiSkewTableaux
from .composition_tableau import CompositionTableau, CompositionTableaux

lazy_import('sage.combinat.tableau_tuple',['TableauTuple', 'StandardTableauTuple', 'RowStandardTableauTuple',
                                           'TableauTuples', 'StandardTableauTuples', 'RowStandardTableauTuples'])
from .k_tableau import WeakTableau, WeakTableaux, StrongTableau, StrongTableaux
lazy_import('sage.combinat.lr_tableau', ['LittlewoodRichardsonTableau',
                                         'LittlewoodRichardsonTableaux'])
lazy_import('sage.combinat.shifted_primed_tableau', ['ShiftedPrimedTableaux',
                                                     'ShiftedPrimedTableau'])

#SuperTableaux
lazy_import('sage.combinat.super_tableau',["StandardSuperTableau", "SemistandardSuperTableau", "StandardSuperTableaux", "SemistandardSuperTableaux"])

#Words
from .words.all import *

from .subword import Subwords

from .graph_path import GraphPaths

#Tuples
from .tuple import Tuples, UnorderedTuples

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
from .ordered_tree import (OrderedTree, OrderedTrees,
                          LabelledOrderedTree, LabelledOrderedTrees)
from .binary_tree import (BinaryTree, BinaryTrees,
                         LabelledBinaryTree, LabelledBinaryTrees)
lazy_import('sage.combinat.interval_posets', ['TamariIntervalPoset', 'TamariIntervalPosets'])
lazy_import('sage.combinat.rooted_tree', ('RootedTree', 'RootedTrees',
                         'LabelledRootedTree', 'LabelledRootedTrees'))

from .combination import Combinations

from .set_partition import SetPartition, SetPartitions
from .set_partition_ordered import OrderedSetPartition, OrderedSetPartitions
lazy_import('sage.combinat.multiset_partition_into_sets_ordered', ['OrderedMultisetPartitionIntoSets',
                                                         'OrderedMultisetPartitionsIntoSets'])
from .subset import Subsets
#from subsets_pairwise import PairwiseCompatibleSubsets
from .necklace import Necklaces
lazy_import('sage.combinat.dyck_word', ('DyckWords', 'DyckWord'))
from .sloane_functions import sloane
lazy_import('sage.combinat.superpartition', ('SuperPartition',
                                             'SuperPartitions'))

lazy_import('sage.combinat.parallelogram_polyomino',
            ['ParallelogramPolyomino', 'ParallelogramPolyominoes'])

from .root_system.all import *
from .sf.all import *
from .ncsf_qsym.all import *
from .ncsym.all import *
lazy_import('sage.combinat.fqsym', 'FreeQuasisymmetricFunctions')
from .matrices.all import *
# Posets
from .posets.all import *

# Cluster Algebras and Quivers
from .cluster_algebra_quiver.all import *

from . import ranker

from .integer_vector import IntegerVectors
from .integer_vector_weighted import WeightedIntegerVectors
from .integer_vectors_mod_permgroup import IntegerVectorsModPermutationGroup

lazy_import('sage.combinat.q_analogues', ['gaussian_binomial', 'q_binomial'])

from .species.all import *

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
lazy_import('sage.combinat.k_regular_sequence', 'kRegularSequenceSpace')

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
