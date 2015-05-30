"""
Combinatorics features that are imported by default in the interpreter namespace
"""
from combinat import bell_number, catalan_number, euler_number, fibonacci, \
        lucas_number1, lucas_number2, stirling_number1, stirling_number2, \
        CombinatorialObject, CombinatorialClass, FilteredCombinatorialClass, \
        UnionCombinatorialClass, MapCombinatorialClass, \
        InfiniteAbstractCombinatorialClass, \
        tuples, number_of_tuples, \
        unordered_tuples, number_of_unordered_tuples, \
        permutations, cyclic_permutations, \
        cyclic_permutations_iterator, bell_polynomial, fibonacci_sequence, \
        fibonacci_xrange, bernoulli_polynomial

from expnums import expnums

from sage.combinat.crystals.all import *
from rigged_configurations.all import *

from sage.combinat.dlx import DLXMatrix, AllExactCovers, OneExactCover

# block designs, etc
from sage.combinat.designs.all import *

# Free modules and friends
from free_module import CombinatorialFreeModule
from combinatorial_algebra import CombinatorialAlgebra
from debruijn_sequence import DeBruijnSequences

from schubert_polynomial import SchubertPolynomialRing
from symmetric_group_algebra import SymmetricGroupAlgebra, HeckeAlgebraSymmetricGroupT
from symmetric_group_representations import SymmetricGroupRepresentation, SymmetricGroupRepresentations
from yang_baxter_graph import YangBaxterGraph
#from hall_littlewood import HallLittlewood_qp, HallLittlewood_q, HallLittlewood_p

#Permutations
from permutation import Permutation, Permutations, Arrangements, PermutationOptions, CyclicPermutations, CyclicPermutationsOfPartition
from affine_permutation import AffinePermutationGroup
from derangements import Derangements
lazy_import('sage.combinat.baxter_permutations', ['BaxterPermutations'])

#RSK
from rsk import RSK, RSK_inverse, robinson_schensted_knuth, robinson_schensted_knuth_inverse,\
                RobinsonSchenstedKnuth, RobinsonSchenstedKnuth_inverse

#PerfectMatchings
from perfect_matching import PerfectMatching, PerfectMatchings

# Integer lists lex

from integer_list import IntegerListsLex as IntegerListsLex

#Compositions
from composition import Composition, Compositions
from composition_signed import SignedCompositions

#Partitions
from partition import Partition, Partitions, PartitionsInBox,\
     OrderedPartitions, PartitionsGreatestLE, PartitionsGreatestEQ,\
     PartitionsGreatestLE, PartitionsGreatestEQ, number_of_partitions
#Functions being deprecated from partition
from partition import RestrictedPartitions

from sage.combinat.partition_tuple import PartitionTuple, PartitionTuples
from skew_partition import SkewPartition, SkewPartitions

#Partition algebra
from partition_algebra import SetPartitionsAk, SetPartitionsPk, SetPartitionsTk, SetPartitionsIk, SetPartitionsBk, SetPartitionsSk, SetPartitionsRk, SetPartitionsRk, SetPartitionsPRk

#Diagram algebra
from diagram_algebras import PartitionAlgebra, BrauerAlgebra, TemperleyLiebAlgebra, PlanarAlgebra, PropagatingIdeal

#Descent algebra
from descent_algebra import DescentAlgebra

#Vector Partitions
from vector_partition import VectorPartition, VectorPartitions

#Similarity class types
from similarity_class_type import PrimarySimilarityClassType, PrimarySimilarityClassTypes, SimilarityClassType, SimilarityClassTypes

#Cores
from core import Core, Cores

#Tableaux
from tableau import Tableau, SemistandardTableau, StandardTableau, \
        Tableaux, StandardTableaux, SemistandardTableaux
from skew_tableau import SkewTableau, SkewTableaux, StandardSkewTableaux, SemistandardSkewTableaux
from ribbon_shaped_tableau import RibbonShapedTableau, RibbonShapedTableaux, StandardRibbonShapedTableaux
from ribbon_tableau import RibbonTableaux, RibbonTableau, MultiSkewTableaux, MultiSkewTableau, SemistandardMultiSkewTableaux
from composition_tableau import CompositionTableau, CompositionTableaux

from sage.combinat.tableau_tuple import TableauTuple, StandardTableauTuple, TableauTuples, StandardTableauTuples
from k_tableau import WeakTableau, WeakTableaux, StrongTableau, StrongTableaux

#Words
from words.all import *

from subword import Subwords

from graph_path import GraphPaths

#Tuples
from tuple import Tuples, UnorderedTuples

#Alternating sign matrices
from alternating_sign_matrix import AlternatingSignMatrix, AlternatingSignMatrices, MonotoneTriangles, ContreTableaux, TruncatedStaircases

# Parking Functions
from non_decreasing_parking_function import NonDecreasingParkingFunctions, NonDecreasingParkingFunction
from parking_functions import ParkingFunctions, ParkingFunction

# Trees and Tamari interval posets
from sage.misc.lazy_import import lazy_import
from ordered_tree import (OrderedTree, OrderedTrees,
                          LabelledOrderedTree, LabelledOrderedTrees)
from binary_tree import (BinaryTree, BinaryTrees,
                         LabelledBinaryTree, LabelledBinaryTrees)

lazy_import('sage.combinat.interval_posets', ['TamariIntervalPoset', 'TamariIntervalPosets'])
from rooted_tree import (RootedTree, RootedTrees,
                         LabelledRootedTree, LabelledRootedTrees)

from combination import Combinations
from cartesian_product import CartesianProduct

from set_partition import SetPartition, SetPartitions
from set_partition_ordered import OrderedSetPartition, OrderedSetPartitions
from subset import Subsets
#from subsets_pairwise import PairwiseCompatibleSubsets
from necklace import Necklaces
from lyndon_word import LyndonWord, LyndonWords, StandardBracketedLyndonWords
from dyck_word import DyckWords, DyckWord
from sloane_functions import sloane

from root_system.all import *
from sf.all import *
from ncsf_qsym.all import *
from ncsym.all import *
from matrices.all import *
# Posets
from posets.all import *

from sage.misc.superseded import deprecated_callable_import
deprecated_callable_import(6637,
                           'sage.combinat.backtrack',
                           globals(),
                           locals(),
                           ["SearchForest", 
                            "TransitiveIdeal", 
                            "TransitiveIdealGraded"],
                           ("This class soon will not be available in that "
                            "way anymore. Use RecursivelyEnumeratedSet "
                            "instead."))

# Cluster Algebras and Quivers
from cluster_algebra_quiver.all import *

#import lrcalc

import ranker

from integer_vector import IntegerVectors
from integer_vector_weighted import WeightedIntegerVectors
from integer_vectors_mod_permgroup import IntegerVectorsModPermutationGroup

from finite_class import FiniteCombinatorialClass

from q_analogues import gaussian_binomial, q_binomial

from species.all import *

from multichoose_nk import MultichooseNK

from kazhdan_lusztig import KazhdanLusztigPolynomial

from degree_sequences import DegreeSequences

from cyclic_sieving_phenomenon import CyclicSievingPolynomial, CyclicSievingCheck

from sidon_sets import sidon_sets

# Puzzles
from knutson_tao_puzzles import KnutsonTaoPuzzleSolver

# Gelfand-Tsetlin patterns
from gelfand_tsetlin_patterns import GelfandTsetlinPattern, GelfandTsetlinPatterns

# Finite State Machines (Automaton, Transducer)
lazy_import('sage.combinat.finite_state_machine',
            ['Automaton', 'Transducer', 'FiniteStateMachine'])
lazy_import('sage.combinat.finite_state_machine_generators',
            ['transducers'])
# Binary Recurrence Sequences
from binary_recurrence_sequences import BinaryRecurrenceSequence

# Six Vertex Model
lazy_import('sage.combinat.six_vertex_model', 'SixVertexModel')

