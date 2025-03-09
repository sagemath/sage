"""
Examples of Combinatorial Species
"""
# ****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .set_species import SetSpecies
from .partition_species import PartitionSpecies
from .subset_species import SubsetSpecies
from .recursive_species import CombinatorialSpecies
from .characteristic_species import CharacteristicSpecies, SingletonSpecies, EmptySetSpecies
from .cycle_species import CycleSpecies
from .linear_order_species import LinearOrderSpecies
from .permutation_species import PermutationSpecies
from .empty_species import EmptySpecies
from .sum_species import SumSpecies
from .product_species import ProductSpecies
from .composition_species import CompositionSpecies
from .functorial_composition_species import FunctorialCompositionSpecies

from sage.misc.cachefunc import cached_function


@cached_function
def SimpleGraphSpecies():
    """
    Return the species of simple graphs.

    EXAMPLES::

        sage: S = species.SimpleGraphSpecies()
        sage: S.generating_series().counts(10)
        [1, 1, 2, 8, 64, 1024, 32768, 2097152, 268435456, 68719476736]
        sage: S.cycle_index_series()[:5]                                                # needs sage.modules
        [p[],
         p[1],
         p[1, 1] + p[2],
         4/3*p[1, 1, 1] + 2*p[2, 1] + 2/3*p[3],
         8/3*p[1, 1, 1, 1] + 4*p[2, 1, 1] + 2*p[2, 2] + 4/3*p[3, 1] + p[4]]
        sage: S.isotype_generating_series()[:6]                                         # needs sage.modules
        [1, 1, 2, 4, 11, 34]

    TESTS::

        sage: seq = S.isotype_generating_series().counts(6)[1:]                         # needs sage.modules
        sage: oeis(seq)[0]                              # optional - internet           # needs sage.modules
        A000088: ...

    ::

        sage: seq = S.generating_series().counts(10)[1:]                                # needs sage.modules
        sage: oeis(seq)[0]                              # optional - internet           # needs sage.modules
        A006125: a(n) = 2^(n*(n-1)/2).
    """
    E = SetSpecies()
    E2 = SetSpecies(size=2)
    WP = SubsetSpecies()
    P2 = E2 * E
    return WP.functorial_composition(P2)


@cached_function
def BinaryTreeSpecies():
    r"""
    Return the species of binary trees on `n` leaves.

    The species of binary trees `B` is defined by `B = X + B \cdot B`,
    where `X` is the singleton species.

    EXAMPLES::

        sage: B = species.BinaryTreeSpecies()
        sage: B.generating_series().counts(10)
        [0, 1, 2, 12, 120, 1680, 30240, 665280, 17297280, 518918400]
        sage: B.isotype_generating_series().counts(10)
        [0, 1, 1, 2, 5, 14, 42, 132, 429, 1430]
        sage: B._check()
        True

    ::

        sage: B = species.BinaryTreeSpecies()
        sage: a = B.structures([1,2,3,4,5])[187]; a
        2*((5*3)*(4*1))
        sage: a.automorphism_group()                                                    # needs sage.groups
        Permutation Group with generators [()]

    TESTS::

        sage: seq = B.isotype_generating_series().counts(10)[1:]
        sage: oeis(seq)[0]                              # optional -- internet
        A000108: Catalan numbers: ...
    """
    B = CombinatorialSpecies(min=1)
    X = SingletonSpecies()
    B.define(X + B * B)
    return B


@cached_function
def BinaryForestSpecies():
    """
    Return the species of binary forests.

    Binary forests are defined as sets of binary trees.

    EXAMPLES::

        sage: F = species.BinaryForestSpecies()
        sage: F.generating_series().counts(10)
        [1, 1, 3, 19, 193, 2721, 49171, 1084483, 28245729, 848456353]
        sage: F.isotype_generating_series().counts(10)                                  # needs sage.modules
        [1, 1, 2, 4, 10, 26, 77, 235, 758, 2504]
        sage: F.cycle_index_series()[:7]                                                # needs sage.modules
        [p[],
         p[1],
         3/2*p[1, 1] + 1/2*p[2],
         19/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3],
         193/24*p[1, 1, 1, 1] + 3/4*p[2, 1, 1] + 5/8*p[2, 2] + 1/3*p[3, 1] + 1/4*p[4],
         907/40*p[1, 1, 1, 1, 1] + 19/12*p[2, 1, 1, 1] + 5/8*p[2, 2, 1] + 1/2*p[3, 1, 1] + 1/6*p[3, 2] + 1/4*p[4, 1] + 1/5*p[5],
         49171/720*p[1, 1, 1, 1, 1, 1] + 193/48*p[2, 1, 1, 1, 1] + 15/16*p[2, 2, 1, 1] + 61/48*p[2, 2, 2] + 19/18*p[3, 1, 1, 1] + 1/6*p[3, 2, 1] + 7/18*p[3, 3] + 3/8*p[4, 1, 1] + 1/8*p[4, 2] + 1/5*p[5, 1] + 1/6*p[6]]

    TESTS::

        sage: seq = F.isotype_generating_series().counts(10)[1:]                        # needs sage.modules
        sage: oeis(seq)[0]                              # optional - internet           # needs sage.modules
        A052854: Number of forests of ordered trees on n total nodes.
    """
    B = BinaryTreeSpecies()
    S = SetSpecies()
    F = S(B)
    return F


del cached_function  # so it doesn't get picked up by tab completion
