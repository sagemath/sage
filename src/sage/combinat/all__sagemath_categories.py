# Integer lists
from .matrices.all__sagemath_categories import *

from .integer_lists import IntegerListsLex
from .integer_vector import IntegerVectors

from .combinat import (CombinatorialClass, CombinatorialObject,
                       MapCombinatorialClass,
                       bell_number, bell_polynomial, bernoulli_polynomial,
                       catalan_number, euler_number,
                       fibonacci, fibonacci_sequence, fibonacci_xrange,
                       lucas_number1, lucas_number2,
                       number_of_tuples, number_of_unordered_tuples,
                       polygonal_number, stirling_number1, stirling_number2,
                       tuples, unordered_tuples)

lazy_import('sage.combinat.combinat',
            ('InfiniteAbstractCombinatorialClass', 'UnionCombinatorialClass',
             'FilteredCombinatorialClass'),
            deprecation=(31545, 'this class is deprecated, do not use'))

from .combination import Combinations
from .composition import Composition, Compositions
from .permutation import Permutation, Permutations, Arrangements, CyclicPermutations, CyclicPermutationsOfPartition
from .subset import Subsets, subsets, powerset, uniq
from .tuple import Tuples, UnorderedTuples


from .dlx import DLXMatrix, AllExactCovers, OneExactCover
