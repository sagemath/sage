# -*- coding: utf-8 -*-
r"""
Permutations

The Permutations module. Use ``Permutation?`` to get information about
the Permutation class, and ``Permutations?`` to get information about
the combinatorial class of permutations.

.. WARNING::

   This file defined :class:`Permutation` which depends upon
   :class:`CombinatorialElement` despite it being deprecated (see
   :trac:`13742`). This is dangerous. In particular, the
   :meth:`Permutation._left_to_right_multiply_on_right` method (which can
   be called through multiplication) disables the input checks (see
   :meth:`Permutation`). This should not happen. Do not trust the results.

What does this file define ?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The main part of this file consists in the definition of permutation objects,
i.e. the :meth:`Permutation` method and the
:class:`~sage.combinat.permutation.Permutation` class. Global options for
elements of the permutation class can be set through the
:meth:`Permutations.options` object.

Below are listed all methods and classes defined in this file.

**Methods of Permutations objects**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.combinat.permutation.Permutation.left_action_product` | Returns the product of ``self`` with another permutation, in which the other permutation is applied first.
    :meth:`~sage.combinat.permutation.Permutation.right_action_product` | Returns the product of ``self`` with another permutation, in which ``self`` is applied first.
    :meth:`~sage.combinat.permutation.Permutation.size` | Returns the size of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.cycle_string` | Returns the disjoint-cycles representation of ``self`` as string.
    :meth:`~sage.combinat.permutation.Permutation.next` | Returns the permutation that follows ``self`` in lexicographic order (in the same symmetric group as ``self``).
    :meth:`~sage.combinat.permutation.Permutation.prev` | Returns the permutation that comes directly before ``self`` in lexicographic order (in the same symmetric group as ``self``).
    :meth:`~sage.combinat.permutation.Permutation.to_tableau_by_shape` | Returns a tableau of shape ``shape`` with the entries in ``self``.
    :meth:`~sage.combinat.permutation.Permutation.to_cycles` | Returns the permutation ``self`` as a list of disjoint cycles.
    :meth:`~sage.combinat.permutation.Permutation.forget_cycles` | Return ``self`` under the forget cycle map.
    :meth:`~sage.combinat.permutation.Permutation.to_permutation_group_element` | Returns a ``PermutationGroupElement`` equal to ``self``.
    :meth:`~sage.combinat.permutation.Permutation.signature` | Returns the signature of the permutation ``sef``.
    :meth:`~sage.combinat.permutation.Permutation.is_even` | Returns ``True`` if the permutation ``self`` is even, and ``False`` otherwise.
    :meth:`~sage.combinat.permutation.Permutation.to_matrix` | Returns a matrix representing the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.rank` | Returns the rank of ``self`` in lexicographic ordering (on the symmetric group containing ``self``).
    :meth:`~sage.combinat.permutation.Permutation.to_inversion_vector` | Returns the inversion vector of a permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.inversions` | Returns a list of the inversions of permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.stack_sort` | Returns the permutation obtained by sorting ``self`` through one stack.
    :meth:`~sage.combinat.permutation.Permutation.to_digraph` | Return a digraph representation of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.show` | Displays the permutation as a drawing.
    :meth:`~sage.combinat.permutation.Permutation.number_of_inversions` | Returns the number of inversions in the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.noninversions` | Returns the ``k``-noninversions in the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.number_of_noninversions` | Returns the number of ``k``-noninversions in the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.length` | Returns the Coxeter length of a permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.inverse` | Returns the inverse of a permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.ishift` | Returns the ``i``-shift of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.iswitch` | Returns the ``i``-switch of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.runs` | Returns a list of the runs in the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.longest_increasing_subsequence_length` | Returns the length of the longest increasing subsequences of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.longest_increasing_subsequences` | Returns the list of the longest increasing subsequences of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.cycle_type` | Returns the cycle type of ``self`` as a partition of ``len(self)``.
    :meth:`~sage.combinat.permutation.Permutation.foata_bijection` | Returns the image of the permutation ``self`` under the Foata bijection `\phi`.
    :meth:`~sage.combinat.permutation.Permutation.foata_bijection_inverse` | Returns the image of the permutation ``self`` under the inverse of the Foata bijection `\phi`.
    :meth:`~sage.combinat.permutation.Permutation.fundamental_transformation` | Returns the image of the permutation ``self`` under the Renyi-Foata-Schuetzenberger fundamental transformation.
    :meth:`~sage.combinat.permutation.Permutation.fundamental_transformation_inverse` | Returns the image of the permutation ``self`` under the inverse of the Renyi-Foata-Schuetzenberger fundamental transformation.
    :meth:`~sage.combinat.permutation.Permutation.destandardize` |  Return destandardization of ``self`` with respect to ``weight`` and ``ordered_alphabet``.
    :meth:`~sage.combinat.permutation.Permutation.to_lehmer_code` | Returns the Lehmer code of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.to_lehmer_cocode` | Returns the Lehmer cocode of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.reduced_word` | Returns the reduced word of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.reduced_words` | Returns a list of the reduced words of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.reduced_words_iterator` | An iterator for the reduced words of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.reduced_word_lexmin` | Returns a lexicographically minimal reduced word of a permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.fixed_points` | Returns a list of the fixed points of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.number_of_fixed_points` | Returns the number of fixed points of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.recoils` | Returns the list of the positions of the recoils of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.number_of_recoils` | Returns the number of recoils of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.recoils_composition` | Returns the composition corresponding to the recoils of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.descents` | Returns the list of the descents of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.idescents` | Returns a list of the idescents of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.idescents_signature` | Returns the list obtained by mapping each position in ``self`` to `-1` if it is an idescent and `1` if it is not an idescent.
    :meth:`~sage.combinat.permutation.Permutation.number_of_descents` | Returns the number of descents of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.number_of_idescents` | Returns the number of idescents of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.descents_composition` | Returns the composition corresponding to the descents of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.descent_polynomial` | Returns the descent polynomial of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.major_index` | Returns the major index of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.imajor_index` | Returns the inverse major index of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.to_major_code` | Returns the major code of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.peaks` | Returns a list of the peaks of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.number_of_peaks` | Returns the number of peaks of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.saliances` | Returns a list of the saliances of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.number_of_saliances` | Returns the number of saliances of the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_lequal` | Returns ``True`` if self is less or equal to ``p2`` in the Bruhat order.
    :meth:`~sage.combinat.permutation.Permutation.weak_excedences` | Returns all the numbers ``self[i]`` such that ``self[i] >= i+1``.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_inversions` | Returns the list of inversions of ``self`` such that the application of this inversion to ``self`` decrements its number of inversions.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_inversions_iterator` | Returns an iterator over Bruhat inversions of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_succ` | Returns a list of the permutations covering ``self`` in the Bruhat order.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_succ_iterator` | An iterator for the permutations covering ``self`` in the Bruhat order.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_pred` | Returns a list of the permutations covered by ``self`` in the Bruhat order.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_pred_iterator` | An iterator for the permutations covered by ``self`` in the Bruhat order.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_smaller` | Returns the combinatorial class of permutations smaller than or equal to ``self`` in the Bruhat order.
    :meth:`~sage.combinat.permutation.Permutation.bruhat_greater` | Returns the combinatorial class of permutations greater than or equal to ``self`` in the Bruhat order.
    :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal` | Returns ``True`` if ``self`` is less or equal to ``p2`` in the permutohedron order.
    :meth:`~sage.combinat.permutation.Permutation.permutohedron_succ` | Returns a list of the permutations covering ``self`` in the permutohedron order.
    :meth:`~sage.combinat.permutation.Permutation.permutohedron_pred` | Returns a list of the permutations covered by ``self`` in the permutohedron order.
    :meth:`~sage.combinat.permutation.Permutation.permutohedron_smaller` | Returns a list of permutations smaller than or equal to ``self`` in the permutohedron order.
    :meth:`~sage.combinat.permutation.Permutation.permutohedron_greater` | Returns a list of permutations greater than or equal to ``self`` in the permutohedron order.
    :meth:`~sage.combinat.permutation.Permutation.right_permutohedron_interval_iterator` | Returns an iterator over permutations in an interval of the permutohedron order.
    :meth:`~sage.combinat.permutation.Permutation.right_permutohedron_interval` | Returns a list of permutations in an interval of the permutohedron order.
    :meth:`~sage.combinat.permutation.Permutation.has_pattern` | Tests whether the permutation ``self`` matches the pattern.
    :meth:`~sage.combinat.permutation.Permutation.avoids` | Tests whether the permutation ``self`` avoids the pattern.
    :meth:`~sage.combinat.permutation.Permutation.pattern_positions` | Returns the list of positions where the pattern ``patt`` appears in ``self``.
    :meth:`~sage.combinat.permutation.Permutation.reverse` | Returns the permutation obtained by reversing the 1-line notation of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.complement` | Returns the complement of the permutation which is obtained by replacing each value `x` in the 1-line notation of ``self`` with `n - x + 1`.
    :meth:`~sage.combinat.permutation.Permutation.permutation_poset` | Returns the permutation poset of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.dict` | Returns a dictionary corresponding to the permutation ``self``.
    :meth:`~sage.combinat.permutation.Permutation.action` | Returns the action of the permutation ``self`` on a list.
    :meth:`~sage.combinat.permutation.Permutation.robinson_schensted` | Returns the pair of standard tableaux obtained by running the Robinson-Schensted Algorithm on ``self``.
    :meth:`~sage.combinat.permutation.Permutation.left_tableau` | Returns the left standard tableau after performing the RSK algorithm.
    :meth:`~sage.combinat.permutation.Permutation.right_tableau` | Returns the right standard tableau after performing the RSK algorithm.
    :meth:`~sage.combinat.permutation.Permutation.increasing_tree` | Returns the increasing tree of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.increasing_tree_shape` | Returns the shape of the increasing tree of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.binary_search_tree` | Returns the binary search tree of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.sylvester_class` | Iterates over the equivalence class of ``self`` under sylvester congruence
    :meth:`~sage.combinat.permutation.Permutation.RS_partition` | Returns the shape of the tableaux obtained by the RSK algorithm.
    :meth:`~sage.combinat.permutation.Permutation.remove_extra_fixed_points` | Returns the permutation obtained by removing any fixed points at the end of ``self``.
    :meth:`~sage.combinat.permutation.Permutation.retract_plain` | Returns the plain retract of ``self`` to a smaller symmetric group `S_m`.
    :meth:`~sage.combinat.permutation.Permutation.retract_direct_product` | Returns the direct-product retract of ``self`` to a smaller symmetric group `S_m`.
    :meth:`~sage.combinat.permutation.Permutation.retract_okounkov_vershik` | Returns the Okounkov-Vershik retract of ``self`` to a smaller symmetric group `S_m`.
    :meth:`~sage.combinat.permutation.Permutation.hyperoctahedral_double_coset_type` | Returns the coset-type of ``self`` as a partition.
    :meth:`~sage.combinat.permutation.Permutation.binary_search_tree_shape` | Returns the shape of the binary search tree of ``self`` (a non labelled binary tree).
    :meth:`~sage.combinat.permutation.Permutation.shifted_concatenation` | Returns the right (or left) shifted concatenation of ``self`` with a permutation ``other``.
    :meth:`~sage.combinat.permutation.Permutation.shifted_shuffle` | Returns the shifted shuffle of ``self`` with a permutation ``other``.

**Other classes defined in this file**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :class:`Permutations` |
    :class:`Permutations_nk` |
    :class:`Permutations_mset` |
    :class:`Permutations_set` |
    :class:`Permutations_msetk` |
    :class:`Permutations_setk` |
    :class:`Arrangements` |
    :class:`Arrangements_msetk` |
    :class:`Arrangements_setk` |
    :class:`StandardPermutations_all` |
    :class:`StandardPermutations_n_abstract` |
    :class:`StandardPermutations_n` |
    :class:`StandardPermutations_descents` |
    :class:`StandardPermutations_recoilsfiner` |
    :class:`StandardPermutations_recoilsfatter` |
    :class:`StandardPermutations_recoils` |
    :class:`StandardPermutations_bruhat_smaller` |
    :class:`StandardPermutations_bruhat_greater` |
    :class:`CyclicPermutations` |
    :class:`CyclicPermutationsOfPartition` |
    :class:`StandardPermutations_avoiding_12` |
    :class:`StandardPermutations_avoiding_21` |
    :class:`StandardPermutations_avoiding_132` |
    :class:`StandardPermutations_avoiding_123` |
    :class:`StandardPermutations_avoiding_321` |
    :class:`StandardPermutations_avoiding_231` |
    :class:`StandardPermutations_avoiding_312` |
    :class:`StandardPermutations_avoiding_213` |
    :class:`StandardPermutations_avoiding_generic` |
    :class:`PatternAvoider` |

**Functions defined in this file**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`from_major_code` | Returns the permutation corresponding to major code ``mc``.
    :meth:`from_permutation_group_element` | Returns a Permutation give a ``PermutationGroupElement`` ``pge``.
    :meth:`from_rank` | Returns the permutation with the specified lexicographic rank.
    :meth:`from_inversion_vector` | Returns the permutation corresponding to inversion vector ``iv``.
    :meth:`from_cycles` | Returns the permutation with given disjoint-cycle representation ``cycles``.
    :meth:`from_lehmer_code` | Returns the permutation with Lehmer code ``lehmer``.
    :meth:`from_reduced_word` | Returns the permutation corresponding to the reduced word ``rw``.
    :meth:`bistochastic_as_sum_of_permutations` | Returns a given bistochastic matrix as a nonnegative linear combination of permutations.
    :meth:`bounded_affine_permutation` | Returns a partial permutation representing the bounded affine permutation of a matrix.
    :meth:`descents_composition_list` | Returns a list of all the permutations in a given descent class (i. e., having a given descents composition).
    :meth:`descents_composition_first` | Returns the smallest element of a descent class.
    :meth:`descents_composition_last` | Returns the largest element of a descent class.
    :meth:`bruhat_lequal` | Returns ``True`` if ``p1`` is less or equal to ``p2`` in the Bruhat order.
    :meth:`permutohedron_lequal` | Returns ``True`` if ``p1`` is less or equal to ``p2`` in the permutohedron order.
    :meth:`to_standard` | Returns a standard permutation corresponding to the permutation ``self``.

AUTHORS:

- Mike Hansen

- Dan Drake (2008-04-07): allow Permutation() to take lists of tuples

- Sébastien Labbé (2009-03-17): added robinson_schensted_inverse

- Travis Scrimshaw:

  * (2012-08-16): ``to_standard()`` no longer modifies input
  * (2013-01-19): Removed RSK implementation and moved to
    :mod:`~sage.combinat.rsk`.
  * (2013-07-13): Removed ``CombinatorialClass`` and moved permutations to the
    category framework.

- Darij Grinberg (2013-09-07): added methods; ameliorated :trac:`14885` by
  exposing and documenting methods for global-independent
  multiplication.

- Travis Scrimshaw (2014-02-05): Made :class:`StandardPermutations_n` a
  finite Weyl group to make it more uniform with :class:`SymmetricGroup`.
  Added ability to compute the conjugacy classes.

Classes and methods
===================
"""

# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from typing import Iterator

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.finite_weyl_groups import FiniteWeylGroups
from sage.categories.finite_permutation_groups import FinitePermutationGroups
from sage.structure.list_clone import ClonableArray
from sage.structure.global_options import GlobalOptions
from sage.libs.gap.libgap import libgap
from sage.rings.all import ZZ, Integer, PolynomialRing
from sage.arith.all import factorial, multinomial
from sage.matrix.matrix_space import MatrixSpace
from sage.combinat.tools import transitive_ideal
from sage.combinat.composition import Composition
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.misc.prandom import sample
from sage.graphs.digraph import DiGraph
import itertools
from .combinat import CombinatorialElement, catalan_number
from sage.misc.cachefunc import cached_method
from .backtrack import GenericBacktracker
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.rsk import RSK, RSK_inverse
from sage.combinat.permutation_cython import (left_action_product,
             right_action_product, left_action_same_n, right_action_same_n,
             map_to_list, next_perm)


class Permutation(CombinatorialElement):
    r"""
    A permutation.

    Converts ``l`` to a permutation on `\{1, 2, \ldots, n\}`.

    INPUT:

    - ``l`` -- Can be any one of the following:

      - an instance of :class:`Permutation`,

      - list of integers, viewed as one-line permutation notation. The
        construction checks that you give an acceptable entry. To avoid
        the check, use the ``check_input`` option.

      - string, expressing the permutation in cycle notation.

      - list of tuples of integers, expressing the permutation in cycle
        notation.

      - a :class:`PermutationGroupElement`

      - a pair of two standard tableaux of the same shape. This yields
        the permutation obtained from the pair using the inverse of the
        Robinson-Schensted algorithm.

    - ``check_input`` (boolean) -- whether to check that input is correct. Slows
       the function down, but ensures that nothing bad happens. This is set to
       ``True`` by default.

    .. WARNING::

        Since :trac:`13742` the input is checked for correctness : it is not
        accepted unless it actually is a permutation on `\{1, \ldots, n\}`. It
        means that some :meth:`Permutation` objects cannot be created anymore
        without setting ``check_input = False``, as there is no certainty that
        its functions can handle them, and this should be fixed in a much
        better way ASAP (the functions should be rewritten to handle those
        cases, and new tests be added).

    .. WARNING::

        There are two possible conventions for multiplying permutations, and
        the one currently enabled in Sage by default is the one which has
        `(pq)(i) = q(p(i))` for any permutations `p \in S_n` and `q \in S_n`
        and any `1 \leq i \leq n`. (This equation looks less strange when
        the action of permutations on numbers is written from the right:
        then it takes the form `i^{pq} = (i^p)^q`, which is an associativity
        law). There is an alternative convention, which has
        `(pq)(i) = p(q(i))` instead. The conventions can be switched at
        runtime using
        :meth:`sage.combinat.permutation.Permutations.options`.
        It is best for code not to rely on this setting being set to a
        particular standard, but rather use the methods
        :meth:`left_action_product` and :meth:`right_action_product` for
        multiplying permutations (these methods don't depend on the setting).
        See :trac:`14885` for more details.

    .. NOTE::

        The ``bruhat*`` methods refer to the *strong* Bruhat order. To use
        the *weak* Bruhat order, look under ``permutohedron*``.

    EXAMPLES::

        sage: Permutation([2,1])
        [2, 1]
        sage: Permutation([2, 1, 4, 5, 3])
        [2, 1, 4, 5, 3]
        sage: Permutation('(1,2)')
        [2, 1]
        sage: Permutation('(1,2)(3,4,5)')
        [2, 1, 4, 5, 3]
        sage: Permutation( ((1,2),(3,4,5)) )
        [2, 1, 4, 5, 3]
        sage: Permutation( [(1,2),(3,4,5)] )
        [2, 1, 4, 5, 3]
        sage: Permutation( ((1,2)) )
        [2, 1]
        sage: Permutation( (1,2) )
        [2, 1]
        sage: Permutation( ((1,2),) )
        [2, 1]
        sage: Permutation( ((1,),) )
        [1]
        sage: Permutation( (1,) )
        [1]
        sage: Permutation( () )
        []
        sage: Permutation( ((),) )
        []
        sage: p = Permutation((1, 2, 5)); p
        [2, 5, 3, 4, 1]
        sage: type(p)
        <class 'sage.combinat.permutation.StandardPermutations_n_with_category.element_class'>

    Construction from a string in cycle notation::

        sage: p = Permutation( '(4,5)' ); p
        [1, 2, 3, 5, 4]

    The size of the permutation is the maximum integer appearing; add
    a 1-cycle to increase this::

        sage: p2 = Permutation( '(4,5)(10)' ); p2
        [1, 2, 3, 5, 4, 6, 7, 8, 9, 10]
        sage: len(p); len(p2)
        5
        10

    We construct a :class:`Permutation` from a
    :class:`PermutationGroupElement`::

        sage: g = PermutationGroupElement([2,1,3])
        sage: Permutation(g)
        [2, 1, 3]

    From a pair of tableaux of the same shape. This uses the inverse
    of the Robinson-Schensted algorithm::

        sage: p = [[1, 4, 7], [2, 5], [3], [6]]
        sage: q = [[1, 2, 5], [3, 6], [4], [7]]
        sage: P = Tableau(p)
        sage: Q = Tableau(q)
        sage: Permutation( (p, q) )
        [3, 6, 5, 2, 7, 4, 1]
        sage: Permutation( [p, q] )
        [3, 6, 5, 2, 7, 4, 1]
        sage: Permutation( (P, Q) )
        [3, 6, 5, 2, 7, 4, 1]
        sage: Permutation( [P, Q] )
        [3, 6, 5, 2, 7, 4, 1]

    TESTS::

        sage: Permutation([()])
        []
        sage: Permutation('()')
        []
        sage: Permutation(())
        []
        sage: Permutation( [1] )
        [1]

    From a pair of empty tableaux ::

        sage: Permutation( ([], []) )
        []
        sage: Permutation( [[], []] )
        []
    """
    @staticmethod
    def __classcall_private__(cls, l, check_input = True):
        """
        Return a permutation in the general permutations parent.

        EXAMPLES::

            sage: P = Permutation([2,1]); P
            [2, 1]
            sage: P.parent()
            Standard permutations
        """
        import sage.combinat.tableau as tableau
        if isinstance(l, Permutation):
            return l
        elif isinstance(l, PermutationGroupElement):
            l = l.domain()
        #if l is a string, then assume it is in cycle notation
        elif isinstance(l, str):
            if l == "()" or l == "":
                return from_cycles(0, [])
            cycles = l.split(")(")
            cycles[0] = cycles[0][1:]
            cycles[-1] = cycles[-1][:-1]
            cycle_list = []
            for c in cycles:
                cycle_list.append([int(_) for _ in c.split(",")])

            return from_cycles(max(max(c) for c in cycle_list), cycle_list)

        #if l is a pair of standard tableaux or a pair of lists
        elif isinstance(l, (tuple, list)) and len(l) == 2 and \
            all(isinstance(x, tableau.Tableau) for x in l):
            return RSK_inverse(*l, output='permutation')
        elif isinstance(l, (tuple, list)) and len(l) == 2 and \
            all(isinstance(x, list) for x in l):
            P,Q = [tableau.Tableau(_) for _ in l]
            return RSK_inverse(P, Q, 'permutation')
        # if it's a tuple or nonempty list of tuples, also assume cycle
        # notation
        elif isinstance(l, tuple) or \
             (isinstance(l, list) and l and
             all(isinstance(x, tuple) for x in l)):
            if l and (isinstance(l[0], (int,Integer)) or len(l[0]) > 0):
                if isinstance(l[0], tuple):
                    n = max(max(x) for x in l)
                    return from_cycles(n, [list(x) for x in l])
                else:
                    n = max(l)
                    return from_cycles(n, [list(l)])
            elif len(l) <= 1:
                return Permutations()([])
            else:
                raise ValueError("cannot convert l (= %s) to a Permutation"%l)

        # otherwise, it gets processed by CombinatorialElement's __init__.
        return Permutations()(l, check_input=check_input)

    def __init__(self, parent, l, check_input=True):
        """
        Constructor. Checks that INPUT is not a mess, and calls
        :class:`CombinatorialElement`. It should not, because
        :class:`CombinatorialElement` is deprecated.

        INPUT:

        - ``l`` -- a list of ``int`` variables

        - ``check_input`` (boolean) -- whether to check that input is
          correct. Slows the function down, but ensures that nothing bad
          happens.

          This is set to ``True`` by default.

        TESTS::

            sage: Permutation([1,2,3])
            [1, 2, 3]
            sage: Permutation([1,2,2,4])
            Traceback (most recent call last):
            ...
            ValueError: An element appears twice in the input. It should not.
            sage: Permutation([1,2,4,-1])
            Traceback (most recent call last):
            ...
            ValueError: the elements must be strictly positive integers
            sage: Permutation([1,2,4,5])
            Traceback (most recent call last):
            ...
            ValueError: The permutation has length 4 but its maximal element is
            5. Some element may be repeated, or an element is missing, but there
            is something wrong with its length.
        """
        l = list(l)

        if check_input and len(l) > 0:
            # Make a copy to sort later
            lst = list(l)

            # Is input a list of positive integers ?
            for i in lst:
                try:
                    i = int(i)
                except TypeError:
                    raise ValueError("the elements must be integer variables")
                if i < 1:
                    raise ValueError("the elements must be strictly positive integers")

            lst.sort()

            # Is the maximum element of the permutation the length of input,
            # or is some integer missing ?
            if int(lst[-1]) != len(lst):
                raise ValueError("The permutation has length "+str(len(lst))+
                                 " but its maximal element is "+
                                 str(int(lst[-1]))+". Some element "+
                                 "may be repeated, or an element is missing"+
                                 ", but there is something wrong with its length.")

            # Do the elements appear only once ?
            previous = lst[0]-1

            for i in lst:
                if i == previous:
                    raise ValueError("An element appears twice in the input. It should not.")
                previous = i

        CombinatorialElement.__init__(self, parent, l)

    def __setstate__(self, state):
        r"""
        In order to maintain backwards compatibility and be able to unpickle a
        old pickle from ``Permutation_class`` we have to override the default
        ``__setstate__``.

        EXAMPLES::

            sage: loads(b'x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+H-\xca--I,\xc9\xcc\xcf\xe3\n@\xb0\xe3\x93s\x12\x8b\x8b\xb9\n\x195\x1b'
            ....:       b'\x0b\x99j\x0b\x995BY\xe33\x12\x8b3\nY\xfc\x80\xac\x9c\xcc\xe2\x92B\xd6\xd8B6\r\x88iE\x99y\xe9\xc5z\x99y%\xa9\xe9'
            ....:       b'\xa9E\\\xb9\x89\xd9\xa9\xf10N!{(\xa3qkP!G\x06\x90a\x04dp\x82\x18\x86@\x06Wji\x92\x1e\x00i\x8d0q')
            [3, 2, 1]
            sage: loads(dumps( Permutation([3,2,1]) )) # indirect doctest
            [3, 2, 1]
        """
        if isinstance(state, dict):   # for old pickles from Permutation_class
            self._set_parent(Permutations())
            self.__dict__ = state
        else:
            self._set_parent(state[0])
            self.__dict__ = state[1]

    @cached_method
    def __hash__(self):
        """
        TESTS::

            sage: d = {}
            sage: p = Permutation([1,2,3])
            sage: d[p] = 1
            sage: d[p]
            1
        """
        try:
            return hash(tuple(self._list))
        except Exception:
            return hash(str(self._list))

    def __str__(self) -> str:
        """
        TESTS::

            sage: Permutations.options.display='list'
            sage: p = Permutation([2,1,3])
            sage: str(p)
            '[2, 1, 3]'
            sage: Permutations.options.display='cycle'
            sage: str(p)
            '(1,2)'
            sage: Permutations.options.display='singleton'
            sage: str(p)
            '(1,2)(3)'
            sage: Permutations.options.display='list'
        """
        return repr(self)

    def _repr_(self) -> str:
        """
        TESTS::

            sage: p = Permutation([2,1,3])
            sage: p
            [2, 1, 3]
            sage: Permutations.options.display='cycle'
            sage: p
            (1,2)
            sage: Permutations.options.display='singleton'
            sage: p
            (1,2)(3)
            sage: Permutations.options.display='reduced_word'
            sage: p
            [1]
            sage: Permutations.options._reset()
        """
        display = self.parent().options.display
        if display == 'list':
            return repr(self._list)
        elif display == 'cycle':
            return self.cycle_string()
        elif display == 'singleton':
            return self.cycle_string(singletons=True)
        elif display == 'reduced_word':
            return repr(self.reduced_word())
        raise ValueError('unknown display option')

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: p = Permutation([2,1,3])
            sage: latex(p)
            [2, 1, 3]
            sage: Permutations.options.latex='cycle'
            sage: latex(p)
            (1 \; 2)
            sage: Permutations.options.latex='singleton'
            sage: latex(p)
            (1 \; 2)(3)
            sage: Permutations.options.latex='reduced_word'
            sage: latex(p)
            s_{1}
            sage: latex(Permutation([1,2,3]))
            1
            sage: Permutations.options.latex_empty_str="e"
            sage: latex(Permutation([1,2,3]))
            e
            sage: Permutations.options.latex='twoline'
            sage: latex(p)
            \begin{pmatrix} 1 & 2 & 3 \\ 2 & 1 & 3 \end{pmatrix}
            sage: Permutations.options._reset()
        """
        display = self.parent().options.latex
        if display == "reduced_word":
            let = self.parent().options.generator_name
            redword = self.reduced_word()
            if not redword:
                return self.parent().options.latex_empty_str
            return " ".join("%s_{%s}" % (let, i) for i in redword)
        if display == "twoline":
            return "\\begin{pmatrix} %s \\\\ %s \\end{pmatrix}" % (
                    " & ".join("%s" % i for i in range(1, len(self._list)+1)),
                    " & ".join("%s" % i for i in self._list))
        if display == "list":
            return repr(self._list)
        if display == "cycle":
            ret = self.cycle_string()
        else: # Must be cycles with singletons
            ret = self.cycle_string(singletons=True)
        return ret.replace(",", " \\; ")

    def _gap_(self, gap):
        """
        Return a GAP version of this permutation.

        EXAMPLES::

            sage: gap(Permutation([1,2,3]))
            ()
            sage: gap(Permutation((1,2,3)))
            (1,2,3)
            sage: type(_)
            <class 'sage.interfaces.gap.GapElement'>
        """
        return self.to_permutation_group_element()._gap_(gap)

    def size(self) -> Integer:
        """
        Return the size of ``self``.

        EXAMPLES::

            sage: Permutation([3,4,1,2,5]).size()
            5
        """
        return len(self)

    grade = size  # for the category SetsWithGrading()

    def cycle_string(self, singletons=False) -> str:
        """
        Return a string of the permutation in cycle notation.

        If ``singletons=True``, it includes 1-cycles in the string.

        EXAMPLES::

            sage: Permutation([1,2,3]).cycle_string()
            '()'
            sage: Permutation([2,1,3]).cycle_string()
            '(1,2)'
            sage: Permutation([2,3,1]).cycle_string()
            '(1,2,3)'
            sage: Permutation([2,1,3]).cycle_string(singletons=True)
            '(1,2)(3)'
        """
        cycles = self.to_cycles(singletons=singletons)
        if not cycles:
            return "()"
        else:
            return "".join("("+",".join(str(l) for l in x)+")" for x in cycles)

    def __next__(self):
        r"""
        Return the permutation that follows ``self`` in lexicographic order on
        the symmetric group containing ``self``. If ``self`` is the last
        permutation, then ``next`` returns ``False``.

        EXAMPLES::

            sage: p = Permutation([1, 3, 2])
            sage: next(p)
            [2, 1, 3]
            sage: p = Permutation([4,3,2,1])
            sage: next(p)
            False

        TESTS::

            sage: p = Permutation([])
            sage: next(p)
            False
        """
        p = self[:]
        n = len(self)
        first = -1

        #Starting from the end, find the first o such that
        #p[o] < p[o+1]
        for i in reversed(range(n - 1)):
            if p[i] < p[i + 1]:
                first = i
                break

        #If first is still -1, then we are already at the last permutation
        if first == -1:
            return False

        #Starting from the end, find the first j such that p[j] > p[first]
        j = n - 1
        while p[j] < p[first]:
            j -= 1

        #Swap positions first and j
        (p[j], p[first]) = (p[first], p[j])

        #Reverse the list between first and the end
        first_half = p[:first+1]
        last_half  = p[first+1:]
        last_half.reverse()
        p = first_half + last_half

        return Permutations()(p)

    next = __next__

    def prev(self):
        r"""
        Return the permutation that comes directly before ``self`` in
        lexicographic order on the symmetric group containing ``self``.
        If ``self`` is the first permutation, then it returns ``False``.

        EXAMPLES::

            sage: p = Permutation([1,2,3])
            sage: p.prev()
            False
            sage: p = Permutation([1,3,2])
            sage: p.prev()
            [1, 2, 3]

        TESTS::

            sage: p = Permutation([])
            sage: p.prev()
            False

        Check that :trac:`16913` is fixed::

            sage: Permutation([1,4,3,2]).prev()
            [1, 4, 2, 3]
        """

        p = self[:]
        n = len(self)
        first = -1

        #Starting from the end, find the first o such that
        #p[o] > p[o+1]
        for i in reversed(range(n - 1)):
            if p[i] > p[i + 1]:
                first = i
                break

        #If first is still -1, that is we didn't find any descents,
        #then we are already at the last permutation
        if first == -1:
            return False

        #Starting from the end, find the first j such that p[j] < p[first]
        j = n - 1
        while p[j] > p[first]:
            j -= 1

        #Swap positions first and j
        (p[j], p[first]) = (p[first], p[j])

        #Reverse the list between first+1 and end
        first_half = p[:first+1]
        last_half  = p[first+1:]
        last_half.reverse()
        p = first_half + last_half

        return Permutations()(p)

    def to_tableau_by_shape(self, shape):
        """
        Return a tableau of shape ``shape`` with the entries
        in ``self``. The tableau is such that the reading word (i. e.,
        the word obtained by reading the tableau row by row, starting
        from the top row in English notation, with each row being
        read from left to right) is ``self``.

        EXAMPLES::

            sage: Permutation([3,4,1,2,5]).to_tableau_by_shape([3,2])
            [[1, 2, 5], [3, 4]]
            sage: Permutation([3,4,1,2,5]).to_tableau_by_shape([3,2]).reading_word_permutation()
            [3, 4, 1, 2, 5]
        """
        import sage.combinat.tableau as tableau
        if sum(shape) != len(self):
            raise ValueError("the size of the partition must be the size of self")

        t = []
        w = list(self)
        for i in reversed(shape):
            t = [ w[:i] ] + t
            w = w[i:]
        return tableau.Tableau(t)

    def to_cycles(self, singletons=True, use_min=True):
        """
        Return the permutation ``self`` as a list of disjoint cycles.

        The cycles are returned in the order of increasing smallest
        elements, and each cycle is returned as a tuple which starts
        with its smallest element.

        If ``singletons=False`` is given, the list does not contain the
        singleton cycles.

        If ``use_min=False`` is given, the cycles are returned in the
        order of increasing *largest* (not smallest) elements, and
        each cycle starts with its largest element.

        EXAMPLES::

            sage: Permutation([2,1,3,4]).to_cycles()
            [(1, 2), (3,), (4,)]
            sage: Permutation([2,1,3,4]).to_cycles(singletons=False)
            [(1, 2)]
            sage: Permutation([2,1,3,4]).to_cycles(use_min=True)
            [(1, 2), (3,), (4,)]
            sage: Permutation([2,1,3,4]).to_cycles(use_min=False)
            [(4,), (3,), (2, 1)]
            sage: Permutation([2,1,3,4]).to_cycles(singletons=False, use_min=False)
            [(2, 1)]

            sage: Permutation([4,1,5,2,6,3]).to_cycles()
            [(1, 4, 2), (3, 5, 6)]
            sage: Permutation([4,1,5,2,6,3]).to_cycles(use_min=False)
            [(6, 3, 5), (4, 2, 1)]

            sage: Permutation([6, 4, 5, 2, 3, 1]).to_cycles()
            [(1, 6), (2, 4), (3, 5)]
            sage: Permutation([6, 4, 5, 2, 3, 1]).to_cycles(use_min=False)
            [(6, 1), (5, 3), (4, 2)]

        The algorithm is of complexity `O(n)` where `n` is the size of the
        given permutation.

        TESTS::

            sage: from sage.combinat.permutation import from_cycles
            sage: for n in range(1,6):
            ....:    for p in Permutations(n):
            ....:       if from_cycles(n, p.to_cycles()) != p:
            ....:          print("There is a problem with {}".format(p))
            ....:          break
            sage: size = 10000
            sage: sample = (Permutations(size).random_element() for i in range(5))
            sage: all(from_cycles(size, p.to_cycles()) == p for p in sample)
            True

        Note: there is an alternative implementation called ``_to_cycle_set``
        which could be slightly (10%) faster for some input (typically for
        permutations of size in the range [100, 10000]). You can run the
        following benchmarks. For small permutations::

            sage: for size in range(9): # not tested
            ....:  print(size)
            ....:  lp = Permutations(size).list()
            ....:  timeit('[p.to_cycles(False) for p in lp]')
            ....:  timeit('[p._to_cycles_set(False) for p in lp]')
            ....:  timeit('[p._to_cycles_list(False) for p in lp]')
            ....:  timeit('[p._to_cycles_orig(False) for p in lp]')

        and larger ones::

            sage: for size in [10, 20, 50, 75, 100, 200, 500, 1000, # not tested
            ....:       2000, 5000, 10000, 15000, 20000, 30000,
            ....:       50000, 80000, 100000]:
            ....:    print(size)
            ....:    lp = [Permutations(size).random_element() for i in range(20)]
            ....:    timeit("[p.to_cycles() for p in lp]")
            ....:    timeit("[p._to_cycles_set() for p in lp]")
            ....:    timeit("[p._to_cycles_list() for p in lp]")
        """
        cycles = []

        l = self[:]

        if use_min:
            groundset = range(len(l))
        else:
            groundset = reversed(range(len(l)))

        # Go through until we've considered every number between 1 and len(l)
        for i in groundset:
            if not l[i]:
                continue
            cycleFirst = i + 1
            cycle = [cycleFirst]
            l[i], next = False, l[i]
            while next != cycleFirst:
                cycle.append( next )
                l[next - 1], next  = False, l[next - 1]
            # Add the cycle to the list of cycles
            if singletons or len(cycle) > 1:
                cycles.append(tuple(cycle))
        return cycles

    cycle_tuples = to_cycles

    def _to_cycles_orig(self, singletons=True):
        r"""
        Return the permutation ``self`` as a list of disjoint cycles.

        EXAMPLES::

            sage: Permutation([2,1,3,4])._to_cycles_orig()
            [(1, 2), (3,), (4,)]
            sage: Permutation([2,1,3,4])._to_cycles_orig(singletons=False)
            [(1, 2)]
        """
        p = self[:]
        cycles = []
        toConsider = -1

        #Create the list [1,2,...,len(p)]
        l = [ i+1 for i in range(len(p))]
        cycle = []

        #Go through until we've considered every number between
        #1 and len(p)
        while l:
            #If we are at the end of a cycle
            #then we want to add it to the cycles list
            if toConsider == -1:
                #Add the cycle to the list of cycles
                if singletons:
                    if cycle:
                        cycles.append(tuple(cycle))
                else:
                    if len(cycle) > 1:
                        cycles.append(tuple(cycle))
                #Start with the first element in the list
                toConsider = l[0]
                l.remove(toConsider)
                cycle = [ toConsider ]
                cycleFirst = toConsider

            #Figure out where the element under consideration
            #gets mapped to.
            next = p[toConsider - 1]

            #If the next element is the first one in the list
            #then we've reached the end of the cycle
            if next == cycleFirst:
                toConsider = -1
            else:
                cycle.append( next )
                l.remove( next )
                toConsider = next

        # When we're finished, add the last cycle
        if singletons:
            if cycle:
                cycles.append(tuple(cycle))
        else:
            if len(cycle) > 1:
                cycles.append(tuple(cycle))
        return cycles

    def _to_cycles_set(self, singletons=True):
        r"""
        Return the permutation ``self`` as a list of disjoint cycles.

        EXAMPLES::

            sage: Permutation([2,1,3,4])._to_cycles_set()
            [(1, 2), (3,), (4,)]
            sage: Permutation([2,1,3,4])._to_cycles_set(singletons=False)
            [(1, 2)]

        TESTS::

            sage: all((p._to_cycles_set(False) == p._to_cycles_orig(False)
            ....:          for i in range(7) for p in Permutations(i)))
            True
        """
        p = self[:]
        cycles = []

        if not singletons:
            #remove the fixed points
            L = set( i+1 for i,pi in enumerate(p) if pi != i+1 )
        else:
            L = set(range(1,len(p)+1))

        #Go through until we've considered every remaining number
        while L:
            # take the first remaining element
            cycleFirst = L.pop()
            next = p[cycleFirst-1]
            cycle = [cycleFirst]
            while next != cycleFirst:
                cycle.append(next)
                L.remove(next)
                next = p[next-1]
            # add the cycle
            cycles.append(tuple(cycle))

        return cycles

    def _to_cycles_list(self, singletons=True):
        r"""
        Return the permutation ``self`` as a list of disjoint cycles.

        EXAMPLES::

            sage: Permutation([2,1,3,4])._to_cycles_list()
            [(1, 2), (3,), (4,)]
            sage: Permutation([2,1,3,4])._to_cycles_list(singletons=False)
            [(1, 2)]

        TESTS::

            sage: all((p._to_cycles_list(False) == p._to_cycles_orig(False)
            ....:          for i in range(7) for p in Permutations(i)))
            True
        """
        p = self[:]
        cycles = []

        if not singletons:
            #remove the fixed points
            L = [i+1 for i,pi in enumerate(p) if pi != i+1]
        else:
            L = list(range(1, len(p) + 1))

        from bisect import bisect_left

        #Go through until we've considered every remaining number
        while L:
            # take the first remaining element
            cycleFirst = L.pop(0)
            next = p[cycleFirst-1]
            cycle = [cycleFirst]
            while next != cycleFirst:
                cycle.append(next)
                # remove next from L
                # we use a binary search to find it
                L.pop(bisect_left(L,next))
                next = p[next-1]
            # add the cycle
            cycles.append(tuple(cycle))

        return cycles


    def to_permutation_group_element(self):
        """
        Return a PermutationGroupElement equal to ``self``.

        EXAMPLES::

            sage: Permutation([2,1,4,3]).to_permutation_group_element()
            (1,2)(3,4)
            sage: Permutation([1,2,3]).to_permutation_group_element()
            ()
        """
        grp = SymmetricGroup(len(self))
        return grp.element_class(self.to_cycles(singletons=False), grp, check=False)

    def signature(self) -> Integer:
        r"""
        Return the signature of the permutation ``self``. This is
        `(-1)^l`, where `l` is the number of inversions of ``self``.

        .. NOTE::

            :meth:`sign` can be used as an alias for :meth:`signature`.

        EXAMPLES::

            sage: Permutation([4, 2, 3, 1, 5]).signature()
            -1
            sage: Permutation([1,3,2,5,4]).sign()
            1
            sage: Permutation([]).sign()
            1
        """
        return (-1)**(len(self)-len(self.to_cycles()))

    #one can also use sign as an alias for signature
    sign = signature

    def is_even(self) -> bool:
        r"""
        Return ``True`` if the permutation ``self`` is even and
        ``False`` otherwise.

        EXAMPLES::

            sage: Permutation([1,2,3]).is_even()
            True
            sage: Permutation([2,1,3]).is_even()
            False
        """
        return self.signature() == 1

    def to_matrix(self):
        r"""
        Return a matrix representing the permutation.

        EXAMPLES::

            sage: Permutation([1,2,3]).to_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        Alternatively::

            sage: matrix(Permutation([1,3,2]))
            [1 0 0]
            [0 0 1]
            [0 1 0]

        Notice that matrix multiplication corresponds to permutation
        multiplication only when the permutation option mult='r2l'

        ::

            sage: Permutations.options.mult='r2l'
            sage: p = Permutation([2,1,3])
            sage: q = Permutation([3,1,2])
            sage: (p*q).to_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: p.to_matrix()*q.to_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: Permutations.options.mult='l2r'
            sage: (p*q).to_matrix()
            [1 0 0]
            [0 0 1]
            [0 1 0]
        """
        # build the dictionary of entries since the matrix is
        # extremely sparse
        entries = { (v-1, i): 1 for i, v in enumerate(self) }
        M = MatrixSpace(ZZ, len(self), sparse=True)
        return M(entries)

    _matrix_ = to_matrix

    @combinatorial_map(name='to alternating sign matrix')
    def to_alternating_sign_matrix(self):
        r"""
        Return a matrix representing the permutation in the
        :class:`AlternatingSignMatrix` class.

        EXAMPLES::

            sage: m = Permutation([1,2,3]).to_alternating_sign_matrix(); m
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: parent(m)
            Alternating sign matrices of size 3
        """
        from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix
        return AlternatingSignMatrix(self.to_matrix().rows())

    def __mul__(self, rp) -> Permutation:
        """
        TESTS::

            sage: p213 = Permutation([2,1,3])
            sage: p312 = Permutation([3,1,2])
            sage: Permutations.options.mult='l2r'
            sage: p213*p312
            [1, 3, 2]
            sage: Permutations.options.mult='r2l'
            sage: p213*p312
            [3, 2, 1]
            sage: Permutations.options.mult='l2r'
        """
        if self.parent().options.mult == 'l2r':
            return self._left_to_right_multiply_on_right(rp)
        else:
            return self._left_to_right_multiply_on_left(rp)

    _mul_ = __mul__ # For ``prod()``

    def __rmul__(self, lp) -> Permutation:
        """
        TESTS::

            sage: p213 = Permutation([2,1,3])
            sage: p312 = Permutation([3,1,2])
            sage: Permutations.options.mult='l2r'
            sage: p213*p312
            [1, 3, 2]
            sage: Permutations.options.mult='r2l'
            sage: p213*p312
            [3, 2, 1]
            sage: Permutations.options.mult='l2r'
        """
        if self.parent().options.mult == 'l2r':
            return self._left_to_right_multiply_on_left(lp)
        else:
            return self._left_to_right_multiply_on_right(lp)

    def left_action_product(self, lp):
        r"""
        Return the permutation obtained by composing ``self`` with
        ``lp`` in such an order that ``lp`` is applied first and
        ``self`` is applied afterwards.

        This is usually denoted by either ``self * lp`` or ``lp * self``
        depending on the conventions used by the author. If the value
        of a permutation `p \in S_n` on an integer
        `i \in \{ 1, 2, \cdots, n \}` is denoted by `p(i)`, then this
        should be denoted by ``self * lp`` in order to have
        associativity (i.e., in order to have
        `(p \cdot q)(i) = p(q(i))` for all `p`, `q` and `i`). If, on
        the other hand, the value of a permutation `p \in S_n` on an
        integer `i \in \{ 1, 2, \cdots, n \}` is denoted by `i^p`, then
        this should be denoted by ``lp * self`` in order to have
        associativity (i.e., in order to have
        `i^{p \cdot q} = (i^p)^q` for all `p`, `q` and `i`).

        EXAMPLES::

            sage: p = Permutation([2,1,3])
            sage: q = Permutation([3,1,2])
            sage: p.left_action_product(q)
            [3, 2, 1]
            sage: q.left_action_product(p)
            [1, 3, 2]
        """
        return Permutations()(left_action_product(self._list, lp[:]))

    _left_to_right_multiply_on_left = left_action_product

    def right_action_product(self, rp):
        r"""
        Return the permutation obtained by composing ``self`` with
        ``rp`` in such an order that ``self`` is applied first and
        ``rp`` is applied afterwards.

        This is usually denoted by either ``self * rp`` or ``rp * self``
        depending on the conventions used by the author. If the value
        of a permutation `p \in S_n` on an integer
        `i \in \{ 1, 2, \cdots, n \}` is denoted by `p(i)`, then this
        should be denoted by ``rp * self`` in order to have
        associativity (i.e., in order to have
        `(p \cdot q)(i) = p(q(i))` for all `p`, `q` and `i`). If, on
        the other hand, the value of a permutation `p \in S_n` on an
        integer `i \in \{ 1, 2, \cdots, n \}` is denoted by `i^p`, then
        this should be denoted by ``self * rp`` in order to have
        associativity (i.e., in order to have
        `i^{p \cdot q} = (i^p)^q` for all `p`, `q` and `i`).

        EXAMPLES::

            sage: p = Permutation([2,1,3])
            sage: q = Permutation([3,1,2])
            sage: p.right_action_product(q)
            [1, 3, 2]
            sage: q.right_action_product(p)
            [3, 2, 1]
        """
        return Permutations()(right_action_product(self._list, rp[:]))

    _left_to_right_multiply_on_right = right_action_product

    def __call__(self, i):
        r"""
        Return the image of the integer `i` under ``self``.

        EXAMPLES::

            sage: p = Permutation([2, 1, 4, 5, 3])
            sage: p(1)
            2
            sage: p = Permutation(((1,2),(4,3,5)))
            sage: p(4)
            3
            sage: p(2)
            1
            sage: p = Permutation([5,2,1,6,3,7,4])
            sage: list(map(p, range(1,8)))
            [5, 2, 1, 6, 3, 7, 4]

        TESTS::

            sage: p = Permutation([5,2,1,6,3,7,4])
            sage: p(-1)
            Traceback (most recent call last):
            ...
            TypeError: i (= -1) must be between 1 and 7
            sage: p(10)
            Traceback (most recent call last):
            ...
            TypeError: i (= 10) must be between 1 and 7
        """
        if isinstance(i,(int,Integer)) and 1 <= i <= len(self):
            return self[i-1]
        else:
            raise TypeError("i (= %s) must be between %s and %s" % (i,1,len(self)))

    ########
    # Rank #
    ########

    def rank(self) -> Integer:
        r"""
        Return the rank of ``self`` in the lexicographic ordering on the
        symmetric group to which ``self`` belongs.

        EXAMPLES::

            sage: Permutation([1,2,3]).rank()
            0
            sage: Permutation([1, 2, 4, 6, 3, 5]).rank()
            10
            sage: perms = Permutations(6).list()
            sage: [p.rank() for p in perms] == list(range(factorial(6)))
            True
        """
        n = len(self)
        factoradic = self.to_lehmer_code()
        return sum(factoradic[n - 1 - i] * factorial(i)
                   for i in reversed(range(n)))

    ##############
    # Inversions #
    ##############

    def to_inversion_vector(self):
        r"""
        Return the inversion vector of ``self``.

        The inversion vector of a permutation `p \in S_n` is defined as
        the vector `(v_1, v_2, \ldots, v_n)`, where `v_i` is the
        number of elements larger than `i` that appear to the left
        of `i` in the permutation `p`.

        The algorithm is of complexity `O(n\log(n))` where `n` is the size of
        the given permutation.

        EXAMPLES::

            sage: Permutation([5,9,1,8,2,6,4,7,3]).to_inversion_vector()
            [2, 3, 6, 4, 0, 2, 2, 1, 0]
            sage: Permutation([8,7,2,1,9,4,6,5,10,3]).to_inversion_vector()
            [3, 2, 7, 3, 4, 3, 1, 0, 0, 0]
            sage: Permutation([3,2,4,1,5]).to_inversion_vector()
            [3, 1, 0, 0, 0]

        TESTS::

            sage: from sage.combinat.permutation import from_inversion_vector
            sage: all(from_inversion_vector(p.to_inversion_vector()) == p
            ....:   for n in range(6) for p in Permutations(n))
            True

            sage: P = Permutations(1000)
            sage: sample = (P.random_element() for i in range(5))
            sage: all(from_inversion_vector(p.to_inversion_vector()) == p
            ....:   for p in sample)
            True
        """
        p = self._list
        l = len(p)
        # lightning fast if the length is less than 3
        # (is it really useful?)
        if l<4:
            if l==0:
                return []
            if l==1:
                return [0]
            if l==2:
                return [p[0]-1,0]
            if l==3:
                if p[0]==1:
                    return [0,p[1]-2,0]
                if p[0]==2:
                    if p[1]==1:
                        return [1,0,0]
                    return [2,0,0]
                return [p[1],1,0]
        # choose the best one
        if l<411:
            return self._to_inversion_vector_small()
        else:
            return self._to_inversion_vector_divide_and_conquer()

    def _to_inversion_vector_orig(self):
        r"""
        Return the inversion vector of ``self``.

        The inversion vector of a permutation `p \in S_n` is defined as
        the vector `(v_1 , v_2 , \ldots , v_n)`, where `v_i` is the
        number of elements larger than `i` that appear to the left
        of `i` in the permutation `p`.

        (This implementation is probably not the most efficient one.)

        EXAMPLES::

            sage: p = Permutation([5,9,1,8,2,6,4,7,3])
            sage: p._to_inversion_vector_orig()
            [2, 3, 6, 4, 0, 2, 2, 1, 0]
        """
        p = self._list
        iv = [0]*len(p)
        for i in range(len(p)):
            for pj in p:
                if pj>i+1:
                    iv[i]+=1
                elif pj == i+1:
                    break
        return iv

    def _to_inversion_vector_small(self):
        r"""
        Return the inversion vector of ``self``.

        The inversion vector of a permutation `p \in S_n` is defined as
        the vector `(v_1, v_2, \ldots, v_n)`, where `v_i` is the
        number of elements larger than `i` that appear to the left
        of `i` in the permutation `p`.

        (This implementation is the best choice for ``5 < size < 420``
        approximately.)

        EXAMPLES::

            sage: p = Permutation([5,9,1,8,2,6,4,7,3])
            sage: p._to_inversion_vector_small()
            [2, 3, 6, 4, 0, 2, 2, 1, 0]
        """
        p = self._list
        l = len(p)+1
        iv = [0]*l
        checked = [1]*l
        for pi in reversed(p):
            checked[pi] = 0
            iv[pi] = sum(checked[pi:])
        return iv[1:]

    def _to_inversion_vector_divide_and_conquer(self):
        r"""
        Return the inversion vector of a permutation ``self``.

        The inversion vector of a permutation `p \in S_n` is defined as
        the vector `(v_1, v_2, \ldots, v_n)`, where `v_i` is the
        number of elements larger than `i` that appear to the left
        of `i` in the permutation `p`.

        (This implementation is the best choice for ``size > 410``
        approximately.)

        EXAMPLES::

            sage: p = Permutation([5,9,1,8,2,6,4,7,3])
            sage: p._to_inversion_vector_divide_and_conquer()
            [2, 3, 6, 4, 0, 2, 2, 1, 0]
        """
        # for big permutations,
        # we use a divide-and-conquer strategy
        # it's a merge sort, plus counting inversions
        def merge_and_countv(ivA_A, ivB_B):
            # iv* is the inversion vector of *
            (ivA, A) = ivA_A
            (ivB, B) = ivB_B
            C = []
            i, j = 0, 0
            ivC = []
            lA, lB = len(A), len(B)
            while(i < lA and j < lB):
                if B[j] < A[i]:
                    C.append(B[j])
                    ivC.append(ivB[j] + lA - i)
                    j += 1
                else:
                    C.append(A[i])
                    ivC.append(ivA[i])
                    i += 1
            if i < lA:
                C.extend(A[i:])
                ivC.extend(ivA[i:])
            else:
                C.extend(B[j:])
                ivC.extend(ivB[j:])
            return ivC,C

        def base_case(L):
            s = sorted(L)
            d = dict((j,i) for (i,j) in enumerate(s))
            iv = [0]*len(L)
            checked = [1]*len(L)
            for pi in reversed(L):
                dpi = d[pi]
                checked[dpi] = 0
                iv[dpi] = sum(checked[dpi:])
            return iv,s

        def sort_and_countv(L):
            if len(L)<250:
                return base_case(L)
            l = len(L)//2
            return merge_and_countv( sort_and_countv(L[:l]),
                                     sort_and_countv(L[l:]) )

        return sort_and_countv(self._list)[0]

    def inversions(self) -> list:
        r"""
        Return a list of the inversions of ``self``.

        An inversion of a permutation `p` is a pair `(i, j)` such that
        `i < j` and `p(i) > p(j)`.

        EXAMPLES::

            sage: Permutation([3,2,4,1,5]).inversions()
            [(1, 2), (1, 4), (2, 4), (3, 4)]
        """
        p = self[:]
        n = len(p)
        return [tuple([i+1,j+1]) for i in range(n-1) for j in range(i+1,n)
                if p[i]>p[j]]

    def stack_sort(self) -> Permutation:
        """
        Return the stack sort of a permutation.

        This is another permutation obtained through the
        process of sorting using one stack. If the result is the identity
        permutation, the original permutation is *stack-sortable*.

        See :wikipedia:`Stack-sortable_permutation`

        EXAMPLES::

            sage: p = Permutation([2,1,5,3,4,9,7,8,6])
            sage: p.stack_sort()
            [1, 2, 3, 4, 5, 7, 6, 8, 9]

            sage: S5 = Permutations(5)
            sage: len([1 for s in S5 if s.stack_sort() == S5.one()])
            42

        TESTS::

            sage: p = Permutation([])
            sage: p.stack_sort()
            []
            sage: p = Permutation([1])
            sage: p.stack_sort()
            [1]
        """
        stack: list[int] = []
        sorted_p: list[int] = []
        for j in self:
            if stack:
                for i in reversed(stack):
                    if i < j:
                        sorted_p.append(i)
                        stack.pop()
                    else:
                        break
            stack.append(j)
        sorted_p.extend(reversed(stack))
        return Permutation(sorted_p)  # type: ignore

    def to_digraph(self) -> DiGraph:
        r"""
        Return a digraph representation of ``self``.

        EXAMPLES::

            sage: d = Permutation([3, 1, 2]).to_digraph()
            sage: d.edges(labels=False)
            [(1, 3), (2, 1), (3, 2)]
            sage: P = Permutations(range(1, 10))
            sage: d = Permutation(P.random_element()).to_digraph()
            sage: all(c.is_cycle() for c in d.strongly_connected_components_subgraphs())
            True

        TESTS::

            sage: d = Permutation([1]).to_digraph()
            sage: d.edges(labels=False)
            [(1, 1)]
        """
        return DiGraph([self, enumerate(self, start=1)],
                       format='vertices_and_edges', loops=True)

    def show(self, representation="cycles", orientation="landscape", **args):
        r"""
        Display the permutation as a drawing.

        INPUT:

        - ``representation`` -- different kinds of drawings are available

          - ``"cycles"`` (default) -- the permutation is displayed as a
            collection of directed cycles

          - ``"braid"`` -- the permutation is displayed as segments linking
            each element `1, ..., n` to its image on a parallel line.

            When using this drawing, it is also possible to display the
            permutation horizontally (``orientation = "landscape"``, default
            option) or vertically (``orientation = "portrait"``).

          - ``"chord-diagram"`` -- the permutation is displayed as a directed
            graph, all of its vertices being located on a circle.

        All additional arguments are forwarded to the ``show`` subcalls.

        EXAMPLES::

            sage: Permutations(20).random_element().show(representation = "cycles")
            sage: Permutations(20).random_element().show(representation = "chord-diagram")
            sage: Permutations(20).random_element().show(representation = "braid")
            sage: Permutations(20).random_element().show(representation = "braid", orientation='portrait')

        TESTS::

            sage: Permutations(20).random_element().show(representation = "modern_art")
            Traceback (most recent call last):
            ...
            ValueError: The value of 'representation' must be equal to 'cycles', 'chord-diagram' or 'braid'
        """
        if representation == "cycles" or representation == "chord-diagram":
            d = self.to_digraph()

            if representation == "cycles":
                d.show(**args)
            else:
                d.show(layout = "circular", **args)

        elif representation == "braid":
            from sage.plot.line import line
            from sage.plot.text import text

            if orientation == "landscape":
                r = lambda x,y : (x,y)
            elif orientation == "portrait":
                r = lambda x,y : (-y,x)
            else:
                raise ValueError("The value of 'orientation' must be either "+
                                 "'landscape' or 'portrait'.")

            p = self[:]

            L = line([r(1,1)])
            for i in range(len(p)):
                L += line([r(i,1.0), r(p[i]-1,0)])
                L += text(str(i), r(i,1.05)) + text(str(i), r(p[i]-1,-.05))

            return L.show(axes = False, **args)

        else:
            raise ValueError("The value of 'representation' must be equal to "+
                             "'cycles', 'chord-diagram' or 'braid'")


    def number_of_inversions(self) -> Integer:
        r"""
        Return the number of inversions in ``self``.

        An inversion of a permutation is a pair of elements `(i, j)`
        with `i < j` and `p(i) > p(j)`.

        REFERENCES:

        - http://mathworld.wolfram.com/PermutationInversion.html

        EXAMPLES::

            sage: Permutation([3, 2, 4, 1, 5]).number_of_inversions()
            4
            sage: Permutation([1, 2, 6, 4, 7, 3, 5]).number_of_inversions()
            6
        """
        return sum(self.to_inversion_vector())

    def noninversions(self, k) -> list:
        r"""
        Return the list of all ``k``-noninversions in ``self``.

        If `k` is an integer and `p \in S_n` is a permutation, then
        a `k`-noninversion in `p` is defined as a strictly increasing
        sequence `(i_1, i_2, \ldots, i_k)` of elements of
        `\{ 1, 2, \ldots, n \}` satisfying
        `p(i_1) < p(i_2) < \cdots < p(i_k)`. (In other words, a
        `k`-noninversion in `p` can be regarded as a `k`-element
        subset of `\{ 1, 2, \ldots, n \}` on which `p` restricts
        to an increasing map.)

        EXAMPLES::

            sage: p = Permutation([3, 2, 4, 1, 5])
            sage: p.noninversions(1)
            [[3], [2], [4], [1], [5]]
            sage: p.noninversions(2)
            [[3, 4], [3, 5], [2, 4], [2, 5], [4, 5], [1, 5]]
            sage: p.noninversions(3)
            [[3, 4, 5], [2, 4, 5]]
            sage: p.noninversions(4)
            []
            sage: p.noninversions(5)
            []

        TESTS::

            sage: q = Permutation([])
            sage: q.noninversions(1)
            []
        """
        if k > len(self):
            return []
        return [list(pos) for pos in itertools.combinations(self, k) if all( pos[i] < pos[i+1] for i in range(k-1) )]

    def number_of_noninversions(self, k) -> Integer:
        r"""
        Return the number of ``k``-noninversions in ``self``.

        If `k` is an integer and `p \in S_n` is a permutation, then
        a `k`-noninversion in `p` is defined as a strictly increasing
        sequence `(i_1, i_2, \ldots, i_k)` of elements of
        `\{ 1, 2, \ldots, n \}` satisfying
        `p(i_1) < p(i_2) < \cdots < p(i_k)`. (In other words, a
        `k`-noninversion in `p` can be regarded as a `k`-element
        subset of `\{ 1, 2, \ldots, n \}` on which `p` restricts
        to an increasing map.)

        The number of `k`-noninversions in `p` has been denoted by
        `\mathrm{noninv}_k(p)` in [RSW2011]_, where conjectures
        and results regarding this number have been stated.

        EXAMPLES::

            sage: p = Permutation([3, 2, 4, 1, 5])
            sage: p.number_of_noninversions(1)
            5
            sage: p.number_of_noninversions(2)
            6
            sage: p.number_of_noninversions(3)
            2
            sage: p.number_of_noninversions(4)
            0
            sage: p.number_of_noninversions(5)
            0

        The number of `2`-noninversions of a permutation `p \in S_n`
        is `\binom{n}{2}` minus its number of inversions::

            sage: b = binomial(5, 2)
            sage: all( x.number_of_noninversions(2) == b - x.number_of_inversions()
            ....:      for x in Permutations(5) )
            True

        We also check some corner cases::

            sage: all( x.number_of_noninversions(1) == 5 for x in Permutations(5) )
            True
            sage: all( x.number_of_noninversions(0) == 1 for x in Permutations(5) )
            True
            sage: Permutation([]).number_of_noninversions(1)
            0
            sage: Permutation([]).number_of_noninversions(0)
            1
            sage: Permutation([2, 1]).number_of_noninversions(3)
            0
        """
        if k > len(self):
            return 0
        return sum(1 for pos in itertools.combinations(self, k)
                   if all(pos[i] < pos[i + 1] for i in range(k - 1)))

    def length(self) -> Integer:
        r"""
        Return the Coxeter length of ``self``.

        The length of a permutation `p` is given by the number of inversions
        of `p`.

        EXAMPLES::

            sage: Permutation([5, 1, 3, 4, 2]).length()
            6
        """
        return self.number_of_inversions()

    def absolute_length(self) -> Integer:
        """
        Return the absolute length of ``self``

        The absolute length is the length of the shortest expression
        of the element as a product of reflections.

        For permutations in the symmetric groups, the absolute
        length is the size minus the number of its disjoint
        cycles.

        EXAMPLES::

            sage: Permutation([4,2,3,1]).absolute_length()
            1
        """
        return self.size() - len(self.cycle_type())

    @combinatorial_map(order=2,name='inverse')
    def inverse(self) -> Permutation:
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: Permutation([3,8,5,10,9,4,6,1,7,2]).inverse()
            [8, 10, 1, 6, 3, 7, 9, 2, 5, 4]
            sage: Permutation([2, 4, 1, 5, 3]).inverse()
            [3, 1, 5, 2, 4]
            sage: ~Permutation([2, 4, 1, 5, 3])
            [3, 1, 5, 2, 4]
        """
        w = list(range(len(self)))
        for i,j in enumerate(self):
            w[j-1] = i+1
        return Permutations()(w)

    __invert__ = inverse

    def _icondition(self, i):
        """
        Return a string which shows the relative positions of `i-1,i,i+1` in
        ``self``, along with the actual positions of these three letters in
        ``self``. The string represents the letters `i-1,i,i+1` by `1,2,3`,
        respectively.

        .. NOTE::

            An imove (that is, an iswitch or an ishift) can only be applied
            when the relative positions of `i-1,i,i+1` are one of '213',
            '132', '231', or '312'. ``None`` is returned in the other cases
            to signal that an imove cannot be applied.

        EXAMPLES::

            sage: Permutation([2,1,3])._icondition(2)
            ('213', 1, 0, 2)
            sage: Permutation([1,3,2])._icondition(2)
            ('132', 0, 2, 1)
            sage: Permutation([2,3,1])._icondition(2)
            ('231', 2, 0, 1)
            sage: Permutation([3,1,2])._icondition(2)
            ('312', 1, 2, 0)
            sage: Permutation([1,2,3])._icondition(2)
            (None, 0, 1, 2)
            sage: Permutation([1,3,2,4])._icondition(3)
            ('213', 2, 1, 3)
            sage: Permutation([2,1,3])._icondition(3)
            Traceback (most recent call last):
            ...
            ValueError: i (= 3) must be between 2 and n-1

        .. SEEALSO::

            :meth:`ishift`, :meth:`iswitch`
        """
        if i not in range(2, len(self)):
            raise ValueError("i (= %s) must be between 2 and n-1" % i)
        pos_i   = self.index(i)
        pos_ip1 = self.index(i+1)
        pos_im1 = self.index(i-1)

        if pos_i < pos_im1 and pos_im1 < pos_ip1:
            state = '213'
        elif pos_im1 < pos_ip1 and pos_ip1 < pos_i:
            state = '132'
        elif pos_i < pos_ip1 and pos_ip1 < pos_im1:
            state = '231'
        elif pos_ip1 < pos_im1 and pos_im1 < pos_i:
            state = '312'
        else:
            state = None

        return (state, pos_im1, pos_i, pos_ip1)

    def ishift(self, i):
        """
        Return the ``i``-shift of ``self``. If an ``i``-shift of ``self``
        can't be performed, then ``self`` is returned.

        An `i`-shift can be applied when `i` is not inbetween `i-1` and
        `i+1`. The `i`-shift moves `i` to the other side, and leaves the
        relative positions of `i-1` and `i+1` in place. All other entries
        of the permutations are also left in place.

        EXAMPLES:

        Here, `2` is to the left of both `1` and `3`. A `2`-shift
        can be applied which moves the `2` to the right and leaves `1` and
        `3` in their same relative order::

            sage: Permutation([2,1,3]).ishift(2)
            [1, 3, 2]

        All entries other than `i`, `i-1` and `i+1` are unchanged::

            sage: Permutation([2,4,1,3]).ishift(2)
            [1, 4, 3, 2]

        Since `2` is between `1` and `3` in ``[1,2,3]``, a `2`-shift cannot
        be applied to ``[1,2,3]`` ::

            sage: Permutation([1,2,3]).ishift(2)
            [1, 2, 3]
        """
        state = self._icondition(i)
        if state[0] is None:
            return self

        state, pos_im1, pos_i, pos_ip1 = state
        l = list(self)

        if state == '213':   #goes to 132
            l[pos_i]   = i-1
            l[pos_im1] = i+1
            l[pos_ip1] = i
        elif state == '132': #goes to 213
            l[pos_im1] = i
            l[pos_ip1] = i-1
            l[pos_i]   = i+1
        elif state == '231': #goes to 312
            l[pos_i]   = i+1
            l[pos_ip1] = i-1
            l[pos_im1] = i
        elif state == '312': #goes to 231
            l[pos_ip1] = i
            l[pos_im1] = i+1
            l[pos_i]   = i-1
        else:
            # This branch should never occur, no matter what the user does.
            raise ValueError("invalid state")

        return Permutations()(l)

    def iswitch(self, i):
        """
        Return the ``i``-switch of ``self``. If an ``i``-switch of ``self``
        can't be performed, then ``self`` is returned.

        An `i`-switch can be applied when the subsequence of ``self`` formed
        by the entries `i-1`, `i` and `i+1` is neither increasing nor
        decreasing. In this case, this subsequence is reversed (i. e., its
        leftmost element and its rightmost element switch places), while all
        other letters of ``self`` are kept in place.

        EXAMPLES:

        Here, `2` is to the left of both `1` and `3`. A `2`-switch can be
        applied which moves the `2` to the right and switches the relative
        order between `1` and `3`::

            sage: Permutation([2,1,3]).iswitch(2)
            [3, 1, 2]

        All entries other than `i-1`, `i` and `i+1` are unchanged::

            sage: Permutation([2,4,1,3]).iswitch(2)
            [3, 4, 1, 2]

        Since `2` is between `1` and `3` in ``[1,2,3]``, a `2`-switch
        cannot be applied to ``[1,2,3]`` ::

            sage: Permutation([1,2,3]).iswitch(2)
            [1, 2, 3]
        """
        if i not in range(2, len(self)):
            raise ValueError("i (= %s) must between 2 and n-1"%i)

        state = self._icondition(i)
        if state[0] is None:
            return self

        state, pos_im1, pos_i, pos_ip1 = state
        l = list(self)

        if state == '213':    #goes to 312
            l[pos_i]   = i+1
            l[pos_ip1] = i
        elif state == '132':  #goes to 231
            l[pos_im1] = i
            l[pos_i]   = i-1
        elif state == '231':  #goes to 132
            l[pos_i]   = i-1
            l[pos_im1] = i
        elif state == '312':  #goes to 213
            l[pos_ip1] = i
            l[pos_i]   = i+1
        else:
            # This branch should never occur, no matter what the user does.
            raise ValueError("invalid state")

        return Permutations()(l)

    def runs(self, as_tuple=False):
        r"""
        Return a list of the runs in the nonempty permutation
        ``self``.

        A run in a permutation is defined to be a maximal (with
        respect to inclusion) nonempty increasing substring (i. e.,
        contiguous subsequence). For instance, the runs in the
        permutation ``[6,1,7,3,4,5,2]`` are ``[6]``, ``[1,7]``,
        ``[3,4,5]`` and ``[2]``.

        Runs in an empty permutation are not defined.

        INPUT:

        - ``as_tuple`` -- boolean (default: ``False``) choice of
          output format

        OUTPUT:

        a list of lists or a tuple of tuples

        REFERENCES:

        - http://mathworld.wolfram.com/PermutationRun.html

        EXAMPLES::

            sage: Permutation([1,2,3,4]).runs()
            [[1, 2, 3, 4]]
            sage: Permutation([4,3,2,1]).runs()
            [[4], [3], [2], [1]]
            sage: Permutation([2,4,1,3]).runs()
            [[2, 4], [1, 3]]
            sage: Permutation([1]).runs()
            [[1]]

        The example from above::

            sage: Permutation([6,1,7,3,4,5,2]).runs()
            [[6], [1, 7], [3, 4, 5], [2]]
            sage: Permutation([6,1,7,3,4,5,2]).runs(as_tuple=True)
            ((6,), (1, 7), (3, 4, 5), (2,))

        The number of runs in a nonempty permutation equals its
        number of descents plus 1::

            sage: all( len(p.runs()) == p.number_of_descents() + 1
            ....:      for p in Permutations(6) )
            True
        """
        p = self[:]
        runs = []
        current_value = p[0]
        current_run = [p[0]]
        for i in p[1:]:
            if i < current_value:
                runs.append(current_run)
                current_run = [i]
            else:
                current_run.append(i)

            current_value = i
        runs.append(current_run)

        if as_tuple:
            return tuple(tuple(r) for r in runs)

        return runs

    def decreasing_runs(self, as_tuple=False):
        """
        Decreasing runs of the permutation.

        INPUT:

        - ``as_tuple`` -- boolean (default: ``False``) choice of output
          format

        OUTPUT:

        a list of lists or a tuple of tuples

        .. SEEALSO::

            :meth:`runs`

        EXAMPLES::

            sage: s = Permutation([2,8,3,9,6,4,5,1,7])
            sage: s.decreasing_runs()
            [[2], [8, 3], [9, 6, 4], [5, 1], [7]]
            sage: s.decreasing_runs(as_tuple=True)
            ((2,), (8, 3), (9, 6, 4), (5, 1), (7,))
        """
        n = len(self)
        s_bar = self.complement()
        if as_tuple:
            return tuple(tuple(n + 1 - i for i in r) for r in s_bar.runs())
        return [[n + 1 - i for i in r] for r in s_bar.runs()]

    def longest_increasing_subsequence_length(self) -> Integer:
        r"""
        Return the length of the longest increasing subsequences of ``self``.

        EXAMPLES::

            sage: Permutation([2,3,1,4]).longest_increasing_subsequence_length()
            3
            sage: all(i.longest_increasing_subsequence_length() == len(RSK(i)[0][0]) for i in Permutations(5))
            True
            sage: Permutation([]).longest_increasing_subsequence_length()
            0
        """
        r: list[int] = []
        for x in self:
            if max(r+[0]) > x:
                y = min(z for z in r if z > x)
                r[r.index(y)] = x
            else:
                r.append(x)
        return len(r)

    def longest_increasing_subsequences(self):
        r"""
        Return the list of the longest increasing subsequences of ``self``.

        .. NOTE::

            The algorithm is not optimal.

        EXAMPLES::

            sage: Permutation([2,3,4,1]).longest_increasing_subsequences()
            [[2, 3, 4]]
            sage: Permutation([5, 7, 1, 2, 6, 4, 3]).longest_increasing_subsequences()
            [[1, 2, 6], [1, 2, 4], [1, 2, 3]]
        """
        patt = list(range(1,self.longest_increasing_subsequence_length()+1))
        return [[self[i] for i in m] for m in self.pattern_positions(patt)]

    def cycle_type(self):
        r"""
        Return a partition of ``len(self)`` corresponding to the cycle
        type of ``self``.

        This is a non-increasing sequence of the cycle lengths of ``self``.

        EXAMPLES::

            sage: Permutation([3,1,2,4]).cycle_type()
            [3, 1]
        """
        cycle_type = [len(c) for c in self.to_cycles()]
        cycle_type.sort(reverse=True)
        from sage.combinat.partition import Partition
        return Partition(cycle_type)

    @combinatorial_map(name='forget cycles')
    def forget_cycles(self):
        r"""
        Return the image of ``self`` under the map which forgets cycles.

        Consider a permutation `\sigma` written in standard cyclic form:

        .. MATH::

            \sigma
            = (a_{1,1}, \ldots, a_{1,k_1})
              (a_{2,1}, \ldots, a_{2,k_2})
              \cdots
              (a_{m,1}, \ldots, a_{m,k_m}),

        where `a_{1,1} < a_{2,1} < \cdots < a_{m,1}` and `a_{j,1} < a_{j,i}`
        for all `1 \leq j \leq m` and `2 \leq i \leq k_j` where we include
        cycles of length 1 as well. The image of the forget cycle map `\phi`
        is given by

        .. MATH::

            \phi(\sigma) = [a_{1,1}, \ldots, a_{1,k_1}, a_{2,1}, \ldots,
            a_{2,k_2}, \ldots, a_{m,1}, \ldots, a_{m,k_m}],

        considered as a permutation in 1-line notation.

        .. SEEALSO::

            :meth:`fundamental_transformation`, which is a similar map that
            is defined by instead taking `a_{j,1} > a_{j,i}` and is
            bijective.

        EXAMPLES::

            sage: P = Permutations(5)
            sage: x = P([1, 5, 3, 4, 2])
            sage: x.forget_cycles()
            [1, 2, 5, 3, 4]

        We select all permutations with a cycle composition of `[2, 3, 1]`
        in `S_6`::

            sage: P = Permutations(6)
            sage: l = [p for p in P if [len(t) for t in p.to_cycles()] == [1,3,2]]

        Next we apply `\phi` and then take the inverse, and then view the
        results as a poset under the Bruhat order::

            sage: l = [p.forget_cycles().inverse() for p in l]
            sage: B = Poset([l, lambda x,y: x.bruhat_lequal(y)])
            sage: R.<q> = QQ[]
            sage: sum(q^B.rank_function()(x) for x in B)
            q^5 + 2*q^4 + 3*q^3 + 3*q^2 + 2*q + 1

        We check the statement in [CC2013]_ that the posets
        `C_{[1,3,1,1]}` and `C_{[1,3,2]}` are isomorphic::

            sage: l2 = [p for p in P if [len(t) for t in p.to_cycles()] == [1,3,1,1]]
            sage: l2 = [p.forget_cycles().inverse() for p in l2]
            sage: B2 = Poset([l2, lambda x,y: x.bruhat_lequal(y)])
            sage: B.is_isomorphic(B2)
            True

        .. SEEALSO::

            :meth:`fundamental_transformation`.
        """
        ret = []
        for t in self.to_cycles():
            ret += list(t)
        return Permutations()(ret)

    @combinatorial_map(name='foata_bijection')
    def foata_bijection(self) -> Permutation:
        r"""
        Return the image of the permutation ``self`` under the Foata
        bijection `\phi`.

        The bijection shows that `\mathrm{maj}` (the major index)
        and `\mathrm{inv}` (the number of inversions) are
        equidistributed: if `\phi(P) = Q`, then `\mathrm{maj}(P) =
        \mathrm{inv}(Q)`.

        The Foata bijection `\phi` is a bijection on the set of words with
        no two equal letters. It can be defined by induction on the size
        of the word: Given a word `w_1 w_2 \cdots w_n`, start with
        `\phi(w_1) = w_1`. At the `i`-th step, if
        `\phi(w_1 w_2 \cdots w_i) = v_1 v_2 \cdots v_i`, we define
        `\phi(w_1 w_2 \cdots w_i w_{i+1})` by placing `w_{i+1}` on the end of
        the word `v_1 v_2 \cdots v_i` and breaking the word up into blocks
        as follows. If `w_{i+1} > v_i`, place a vertical line to the right
        of each `v_k` for which `w_{i+1} > v_k`. Otherwise, if
        `w_{i+1} < v_i`, place a vertical line to the right of each `v_k`
        for which `w_{i+1} < v_k`. In either case, place a vertical line at
        the start of the word as well. Now, within each block between
        vertical lines, cyclically shift the entries one place to the
        right.

        For instance, to compute `\phi([1,4,2,5,3])`, the sequence of
        words is

        * `1`,
        * `|1|4 \to 14`,
        * `|14|2 \to 412`,
        * `|4|1|2|5 \to 4125`,
        * `|4|125|3 \to 45123`.

        So `\phi([1,4,2,5,3]) = [4,5,1,2,3]`.

        See section 2 of [FS1978]_, and the proof of Proposition 1.4.6
        in [EnumComb1]_.

        .. SEEALSO::

            :meth:`foata_bijection_inverse` for the inverse map.

        EXAMPLES::

            sage: Permutation([1,2,4,3]).foata_bijection()
            [4, 1, 2, 3]
            sage: Permutation([2,5,1,3,4]).foata_bijection()
            [2, 1, 3, 5, 4]

            sage: P = Permutation([2,5,1,3,4])
            sage: P.major_index() == P.foata_bijection().number_of_inversions()
            True

            sage: all( P.major_index() == P.foata_bijection().number_of_inversions()
            ....:      for P in Permutations(4) )
            True

        The example from [FS1978]_::

            sage: Permutation([7,4,9,2,6,1,5,8,3]).foata_bijection()
            [4, 7, 2, 6, 1, 9, 5, 8, 3]

        Border cases::

            sage: Permutation([]).foata_bijection()
            []
            sage: Permutation([1]).foata_bijection()
            [1]
        """
        M: list[int] = []
        for e in self:
            k = len(M)
            if k <= 1:
                M.append(e)
                continue

            a = M[-1]
            M_prime = [0]*(k + 1)
            # Locate the positions of the vertical lines.
            if a > e:
                index_list = [-1] + [i for i, val in enumerate(M) if val > e]
            else:
                index_list = [-1] + [i for i, val in enumerate(M) if val < e]

            for j in range(1, len(index_list)):
                start = index_list[j-1] + 1
                end = index_list[j]
                M_prime[start] = M[end]
                for x in range(start + 1, end + 1):
                    M_prime[x] = M[x-1]
            M_prime[k] = e
            M = M_prime
        return Permutations()(M)

    @combinatorial_map(name='foata_bijection_inverse')
    def foata_bijection_inverse(self) -> Permutation:
        r"""
        Return the image of the permutation ``self`` under the inverse
        of the Foata bijection `\phi`.

        See :meth:`foata_bijection` for the definition of the Foata
        bijection.

        EXAMPLES::

            sage: Permutation([4, 1, 2, 3]).foata_bijection()
            [1, 2, 4, 3]

        TESTS::

            sage: all( P.foata_bijection().foata_bijection_inverse() == P
            ....:      for P in Permutations(5) )
            True

        Border cases::

            sage: Permutation([]).foata_bijection_inverse()
            []
            sage: Permutation([1]).foata_bijection_inverse()
            [1]
        """
        L = list(self)
        Mrev = [] # The resulting permutation, in reverse.
        while L:
            e = L.pop()
            Mrev.append(e)
            k = len(L)
            if k <= 1:
                continue
            L_prime = [0]*(k)
            a = L[0]
            # Locate the positions of the vertical lines.
            if a > e:
                index_list = [i for i, val in enumerate(L) if val > e]
            else:
                index_list = [i for i, val in enumerate(L) if val < e]
            index_list.append(k)

            for j in range(1, len(index_list)):
                start = index_list[j-1]
                end = index_list[j] - 1
                L_prime[end] = L[start]
                for x in range(start, end):
                    L_prime[x] = L[x+1]
            L = L_prime
        return Permutations()(reversed(Mrev))

    @combinatorial_map(name='fundamental_transformation')
    def fundamental_transformation(self) -> Permutation:
        r"""
        Return the image of the permutation ``self`` under the
        Renyi-Foata-Schuetzenberger fundamental transformation.

        The fundamental transformation is a bijection from the
        set of all permutations of `\{1, 2, \ldots, n\}` to
        itself, which transforms any such permutation `w`
        as follows:
        Write `w` in cycle form, with each cycle starting with
        its highest element, and the cycles being sorted in
        increasing order of their highest elements.
        Drop the parentheses in the resulting expression, thus
        reading it as a one-line notation of a new permutation
        `u`.
        Then, `u` is the image of `w` under the fundamental
        transformation.

        See [EnumComb1]_, Proposition 1.3.1.

        .. SEEALSO::

            :meth:`fundamental_transformation_inverse`
            for the inverse map.

            :meth:`forget_cycles` for a similar (but non-bijective)
            map where each cycle is starting from its lowest element.

        EXAMPLES::

            sage: Permutation([5, 1, 3, 4, 2]).fundamental_transformation()
            [3, 4, 5, 2, 1]
            sage: Permutations(5)([1, 5, 3, 4, 2]).fundamental_transformation()
            [1, 3, 4, 5, 2]
            sage: Permutation([8, 4, 7, 2, 9, 6, 5, 1, 3]).fundamental_transformation()
            [4, 2, 6, 8, 1, 9, 3, 7, 5]

        Comparison with :meth:`forget_cycles`::

            sage: P = Permutation([(1,3,4),(2,5)])
            sage: P
            [3, 5, 4, 1, 2]
            sage: P.forget_cycles()
            [1, 3, 4, 2, 5]
            sage: P.fundamental_transformation()
            [4, 1, 3, 5, 2]

        TESTS:

        Border cases::

            sage: Permutation([]).fundamental_transformation()
            []
            sage: Permutation([1]).fundamental_transformation()
            [1]
        """
        cycles = self.to_cycles(use_min=False)
        return Permutations()([a for c in reversed(cycles) for a in c])

    @combinatorial_map(name='fundamental_transformation_inverse')
    def fundamental_transformation_inverse(self):
        r"""
        Return the image of the permutation ``self`` under the
        inverse of the Renyi-Foata-Schuetzenberger fundamental
        transformation.

        The inverse of the fundamental transformation is a
        bijection from the set of all permutations of
        `\{1, 2, \ldots, n\}` to itself, which transforms any
        such permutation `w` as follows:
        Let `I = \{ i_1 < i_2 < \cdots < i_k \}` be the set of
        all left-to-right maxima of `w` (that is, of all indices
        `j` such that `w(j)` is bigger than each of
        `w(1), w(2), \ldots, w(j-1)`).
        The image of `w` under the inverse of the fundamental
        transformation is the permutation `u` that sends
        `w(i-1)` to `w(i)` for all `i \notin I` (notice that
        this makes sense, since `1 \in I` whenever `n > 0`),
        while sending each `w(i_p - 1)` (with `p \geq 2`)
        to `w(i_{p-1})`. Here, we set `i_{k+1} = n+1`.

        See [EnumComb1]_, Proposition 1.3.1.

        .. SEEALSO::

            :meth:`fundamental_transformation`
            for the inverse map.

        EXAMPLES::

            sage: Permutation([3, 4, 5, 2, 1]).fundamental_transformation_inverse()
            [5, 1, 3, 4, 2]
            sage: Permutation([4, 2, 6, 8, 1, 9, 3, 7, 5]).fundamental_transformation_inverse()
            [8, 4, 7, 2, 9, 6, 5, 1, 3]

        TESTS::

            sage: all(P.fundamental_transformation_inverse().
            ....:     fundamental_transformation() == P
            ....:     for P in Permutations(4))
            True

            sage: all(P.fundamental_transformation().
            ....:     fundamental_transformation_inverse() == P
            ....:     for P in Permutations(3))
            True

        Border cases::

            sage: Permutation([]).fundamental_transformation_inverse()
            []
            sage: Permutation([1]).fundamental_transformation_inverse()
            [1]
        """
        n = len(self)
        # We shall iterate through ``self`` from left to right.
        # At each step,
        # ``entry`` will be the entry we are currently seeing;
        # ``previous_entry`` will be the previous entry;
        # ``record_value`` will be the largest entry
        #                  encountered so far;
        # ``previous_record`` will be the largest entry
        #                     encountered before ``record_value``.
        record_value = 0
        previous_record = None
        previous_entry = None
        res = [0] * (n+1) # We'll use res[1], res[2], ..., res[n] only.
        for entry in self:
            if entry > record_value:
                record_value = entry
                if previous_record is not None:
                    res[previous_entry] = previous_record
                previous_record = entry
            else:
                res[previous_entry] = entry
            previous_entry = entry
        if n > 0:
            res[previous_entry] = previous_record
        return Permutations()(res[1:])

    def destandardize(self, weight, ordered_alphabet=None):
        r"""
        Return destandardization of ``self`` with respect to ``weight`` and ``ordered_alphabet``.

        INPUT:

        - ``weight`` -- list or tuple of nonnegative integers that sum to `n` if ``self``
          is a permutation in `S_n`.

        - ``ordered_alphabet`` -- (default: None) a list or tuple specifying the ordered alphabet the
          destandardized word is over

        OUTPUT: word over the ``ordered_alphabet`` which standardizes to ``self``

        Let `weight = (w_1,w_2,\ldots,w_\ell)`. Then this methods looks for an increasing
        sequence of `1,2,\ldots, w_1` and labels all letters in it by 1, then an increasing
        sequence of `w_1+1,w_1+2,\ldots,w_1+w_2` and labels all these letters by 2, etc..
        If an increasing sequence for the specified ``weight`` does not exist, an error is
        returned. The output is a word ``w`` over the specified ordered alphabet with
        evaluation ``weight`` such that ``w.standard_permutation()`` is ``self``.

        EXAMPLES::

            sage: p = Permutation([1,2,5,3,6,4])
            sage: p.destandardize([3,1,2])
            word: 113132
            sage: p = Permutation([2,1,3])
            sage: p.destandardize([2,1])
            Traceback (most recent call last):
            ...
            ValueError: Standardization with weight [2, 1] is not possible!

        TESTS::

            sage: p = Permutation([4,1,2,3,5,6])
            sage: p.destandardize([2,1,3], ordered_alphabet = [1,'a',3])
            word: 311a33
            sage: p.destandardize([2,1,3], ordered_alphabet = [1,'a'])
            Traceback (most recent call last):
            ...
            ValueError: Not enough letters in the alphabet are specified compared to the weight
        """
        ides = self.idescents()
        partial = [0]
        for a in weight:
            partial.append(partial[-1]+a)
        if not set(ides).issubset(set(partial)):
            raise ValueError("Standardization with weight {} is not possible!".format(weight))
        if ordered_alphabet is None:
            ordered_alphabet = list(range(1,len(weight)+1))
        else:
            if len(weight) > len(ordered_alphabet):
                raise ValueError("Not enough letters in the alphabet are specified compared to the weight")
        q = self.inverse()
        s = [0]*len(self)
        for i in range(len(partial)-1):
            for j in range(partial[i],partial[i+1]):
                s[q[j]-1] = ordered_alphabet[i]
        from sage.combinat.words.word import Word
        return Word(s)

    def to_lehmer_code(self) -> list:
        r"""
        Return the Lehmer code of the permutation ``self``.

        The Lehmer code of a permutation `p` is defined as the
        list `[c[1],c[2],...,c[n]]`, where `c[i]` is the number of
        `j>i` such that `p(j)<p(i)`.

        EXAMPLES::

            sage: p = Permutation([2,1,3])
            sage: p.to_lehmer_code()
            [1, 0, 0]
            sage: q = Permutation([3,1,2])
            sage: q.to_lehmer_code()
            [2, 0, 0]

            sage: Permutation([1]).to_lehmer_code()
            [0]
            sage: Permutation([]).to_lehmer_code()
            []

        TESTS::

            sage: from sage.combinat.permutation import from_lehmer_code
            sage: all(from_lehmer_code(p.to_lehmer_code()) == p
            ....:   for n in range(6) for p in Permutations(n))
            True

            sage: P = Permutations(1000)
            sage: sample = (P.random_element() for i in range(5))
            sage: all(from_lehmer_code(p.to_lehmer_code()) == p
            ....:   for p in sample)
            True

        """
        l = len(self._list)
        # choose the best implementations
        if l<577:
            return self._to_lehmer_code_small()
        else:
            return self.inverse().to_inversion_vector()

    def _to_lehmer_code_small(self) -> list:
        r"""
        Return the Lehmer code of the permutation ``self``.

        The Lehmer code of a permutation `p` is defined as the
        list `(c_1, c_2, \ldots, c_n)`, where `c_i` is the number
        of `j > i` such that `p(j) < p(i)`.

        (best choice for `size<577` approximately)

        EXAMPLES::

            sage: p = Permutation([7, 6, 10, 2, 3, 4, 8, 1, 9, 5])
            sage: p._to_lehmer_code_small()
            [6, 5, 7, 1, 1, 1, 2, 0, 1, 0]
        """
        p = self._list
        l = len(p)
        lehmer = []
        checked = [1]*l
        for pi in p:
            checked[pi-1] = 0
            lehmer.append(sum(checked[:pi]))
        return lehmer

    def to_lehmer_cocode(self) -> list:
        r"""
        Return the Lehmer cocode of the permutation ``self``.

        The Lehmer cocode of a permutation `p` is defined as the
        list `(c_1, c_2, \ldots, c_n)`, where `c_i` is the number
        of `j < i` such that `p(j) > p(i)`.

        EXAMPLES::

            sage: p = Permutation([2,1,3])
            sage: p.to_lehmer_cocode()
            [0, 1, 0]
            sage: q = Permutation([3,1,2])
            sage: q.to_lehmer_cocode()
            [0, 1, 1]
        """
        p = self[:]
        n = len(p)
        cocode = [0] * n
        for i in range(1, n):
            for j in range(i):
                if p[j] > p[i]:
                    cocode[i] += 1
        return cocode

    #################
    # Reduced Words #
    #################

    def reduced_word(self):
        r"""
        Return a reduced word of the permutation ``self``.

        See :meth:`reduced_words` for the definition of reduced words and
        a way to compute them all.

        .. WARNING::

            This does not respect the multiplication convention.

        EXAMPLES::

            sage: Permutation([3,5,4,6,2,1]).reduced_word()
            [2, 1, 4, 3, 2, 4, 3, 5, 4, 5]

            Permutation([1]).reduced_word_lexmin()
            []
            Permutation([]).reduced_word_lexmin()
            []
        """
        code = self.to_lehmer_code()
        return [i + ci - j for i, ci in enumerate(code) for j in range(ci)]

    def reduced_words_iterator(self):
        r"""
        Return an iterator for the reduced words of ``self``.

        EXAMPLES::

            sage: next(Permutation([5,2,3,4,1]).reduced_words_iterator())
            [1, 2, 3, 4, 3, 2, 1]

        """
        def aux(p):
            is_identity = True
            for d in range(len(p)-1):
                e = d+1
                if p[d] > p[e]:
                    is_identity = False
                    p[d], p[e] = p[e], p[d]

                    for x in aux(p):
                        x.append(e)
                        yield x

                    p[d], p[e] = p[e], p[d]
            if is_identity:
                yield []

        return aux(self[:])

    def reduced_words(self):
        r"""
        Return a list of the reduced words of ``self``.

        The notion of a reduced word is based on the well-known fact
        that every permutation can be written as a product of adjacent
        transpositions. In more detail: If `n` is a nonnegative integer,
        we can define the transpositions `s_i = (i, i+1) \in S_n`
        for all `i \in \{ 1, 2, \ldots, n-1 \}`, and every `p \in S_n`
        can then be written as a product `s_{i_1} s_{i_2} \cdots s_{i_k}`
        for some sequence `(i_1, i_2, \ldots, i_k)` of elements of
        `\{ 1, 2, \ldots, n-1 \}` (here `\{ 1, 2, \ldots, n-1 \}` denotes
        the empty set when `n \leq 1`). Fixing a `p`, the sequences
        `(i_1, i_2, \ldots, i_k)` of smallest length satisfying
        `p = s_{i_1} s_{i_2} \cdots s_{i_k}` are called the reduced words
        of `p`. (Their length is the Coxeter length of `p`, and can be
        computed using :meth:`length`.)

        Note that the product of permutations is defined here in such
        a way that `(pq)(i) = p(q(i))` for all permutations `p` and `q`
        and each `i \in \{ 1, 2, \ldots, n \}` (this is the same
        convention as in :meth:`left_action_product`, but not the
        default semantics of the `*` operator on permutations in Sage).
        Thus, for instance, `s_2 s_1` is the permutation obtained by
        first transposing `1` with `2` and then transposing `2` with `3`.

        .. SEEALSO::

            :meth:`reduced_word`, :meth:`reduced_word_lexmin`

        EXAMPLES::

            sage: Permutation([2,1,3]).reduced_words()
            [[1]]
            sage: Permutation([3,1,2]).reduced_words()
            [[2, 1]]
            sage: Permutation([3,2,1]).reduced_words()
            [[1, 2, 1], [2, 1, 2]]
            sage: Permutation([3,2,4,1]).reduced_words()
            [[1, 2, 3, 1], [1, 2, 1, 3], [2, 1, 2, 3]]

            Permutation([1]).reduced_words()
            [[]]
            Permutation([]).reduced_words()
            [[]]
        """
        return list(self.reduced_words_iterator())

    def reduced_word_lexmin(self):
        r"""
        Return a lexicographically minimal reduced word of the permutation
        ``self``.

        See :meth:`reduced_words` for the definition of reduced words and
        a way to compute them all.

        EXAMPLES::

            sage: Permutation([3,4,2,1]).reduced_word_lexmin()
            [1, 2, 1, 3, 2]

            Permutation([1]).reduced_word_lexmin()
            []
            Permutation([]).reduced_word_lexmin()
            []
        """
        cocode = self.inverse().to_lehmer_cocode()

        rw = []
        for i in range(len(cocode)):
            piece = [j+1 for j in range(i-cocode[i],i)]
            piece.reverse()
            rw += piece

        return rw


    ################
    # Fixed Points #
    ################

    def fixed_points(self) -> list:
        r"""
        Return a list of the fixed points of ``self``.

        EXAMPLES::

            sage: Permutation([1,3,2,4]).fixed_points()
            [1, 4]
            sage: Permutation([1,2,3,4]).fixed_points()
            [1, 2, 3, 4]
        """
        return [i for i, v in enumerate(self, start=1) if i == v]

    def number_of_fixed_points(self) -> Integer:
        r"""
        Return the number of fixed points of ``self``.

        EXAMPLES::

            sage: Permutation([1,3,2,4]).number_of_fixed_points()
            2
            sage: Permutation([1,2,3,4]).number_of_fixed_points()
            4
        """

        return len(self.fixed_points())

    def is_derangement(self) -> bool:
        r"""
        Return if ``self`` is a derangement.

        A permutation `\sigma` is a derangement if `\sigma` has no
        fixed points.

        EXAMPLES::

            sage: P = Permutation([1,4,2,3])
            sage: P.is_derangement()
            False
            sage: P = Permutation([2,3,1])
            sage: P.is_derangement()
            True
        """
        return not self.fixed_points()

    ############
    # Recoils  #
    ############

    def recoils(self) -> list:
        r"""
        Return the list of the positions of the recoils of ``self``.

        A recoil of a permutation `p` is an integer `i` such that `i+1`
        appears to the left of `i` in `p`.
        Here, the positions are being counted starting at `0`.
        (Note that it is the positions, not the recoils themselves, which
        are being listed.)

        EXAMPLES::

            sage: Permutation([1,4,3,2]).recoils()
            [2, 3]
            sage: Permutation([]).recoils()
            []
        """
        p = self
        recoils  = []
        for i in range(len(p)):
            if p[i] != len(self) and self.index(p[i]+1) < i:
                recoils.append(i)

        return recoils

    def number_of_recoils(self) -> Integer:
        r"""
        Return the number of recoils of the permutation ``self``.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).number_of_recoils()
            2
        """
        return len(self.recoils())

    def recoils_composition(self) -> Composition:
        r"""
        Return the recoils composition of ``self``.

        The recoils composition of a permutation `p \in S_n` is the
        composition of `n` whose descent set is the set of the recoils
        of `p` (not their positions). In other words, this is the
        descents composition of `p^{-1}`.

        EXAMPLES::

            sage: Permutation([1,3,2,4]).recoils_composition()
            [2, 2]
            sage: Permutation([]).recoils_composition()
            []
        """
        return self.inverse().descents_composition()

    ############
    # Descents #
    ############

    def descents(self, final_descent=False, side='right', positive=False,
                 from_zero=False, index_set=None):
        r"""
        Return the list of the descents of ``self``.

        A descent of a permutation `p` is an integer `i` such that
        `p(i) > p(i+1)`.

        .. WARNING::

            By default, the descents are returned as elements in the
            index set, i.e., starting at `1`. If you want them to
            start at `0`, set the keyword ``from_zero`` to ``True``.

        INPUT:

        - ``final_descent`` -- boolean (default ``False``);
          if ``True``, the last position of a non-empty
          permutation is also considered as a descent

        - ``side`` -- ``'right'`` (default) or ``'left'``;
          if ``'left'``, return the descents of the inverse permutation

        - ``positive`` -- boolean (default ``False``);
          if ``True``, return the positions that are not descents

        - ``from_zero`` -- boolean (default ``False``);
          if ``True``, return the positions starting from `0`

        - ``index_set`` -- list (default: ``[1, ..., n-1]`` where ``self``
          is a permutation of ``n``); the index set to check for descents

        EXAMPLES::

            sage: Permutation([3,1,2]).descents()
            [1]
            sage: Permutation([1,4,3,2]).descents()
            [2, 3]
            sage: Permutation([1,4,3,2]).descents(final_descent=True)
            [2, 3, 4]
            sage: Permutation([1,4,3,2]).descents(index_set=[1,2])
            [2]
            sage: Permutation([1,4,3,2]).descents(from_zero=True)
            [1, 2]

        TESTS:

        Check that the original error of :trac:`23891` is fixed::

            sage: Permutations(4)([1,4,3,2]).weak_covers()
            [[1, 3, 4, 2], [1, 4, 2, 3]]
        """
        if index_set is None:
            index_set = range(1, len(self))

        if side == 'right':
            p = self
        else:
            p = self.inverse()
        descents = []
        for i in index_set:
            if p[i-1] > p[i]:
                if not positive:
                    descents.append(i)
            else:
                if positive:
                    descents.append(i)

        if final_descent:
            descents.append(len(p))

        if from_zero:
            return [i - 1 for i in descents]

        return descents

    def idescents(self, final_descent=False, from_zero=False):
        """
        Return a list of the idescents of ``self``, that is the list of
        the descents of ``self``'s inverse.

        A descent of a permutation ``p`` is an integer ``i`` such that
        ``p(i) > p(i+1)``.

        .. WARNING::

            By default, the descents are returned as elements in the
            index set, i.e., starting at `1`. If you want them to
            start at `0`, set the keyword ``from_zero`` to ``True``.

        INPUT:

        - ``final_descent`` -- boolean (default ``False``);
          if ``True``, the last position of a non-empty
          permutation is also considered as a descent

        - ``from_zero`` -- optional boolean (default ``False``);
          if ``False``, return the positions starting from `1`

        EXAMPLES::

            sage: Permutation([2,3,1]).idescents()
            [1]
            sage: Permutation([1,4,3,2]).idescents()
            [2, 3]
            sage: Permutation([1,4,3,2]).idescents(final_descent=True)
            [2, 3, 4]
            sage: Permutation([1,4,3,2]).idescents(from_zero=True)
            [1, 2]
        """
        return self.inverse().descents(final_descent=final_descent,
                                       from_zero=from_zero)

    def idescents_signature(self, final_descent=False):
        """
        Return the list obtained as follows: Each position in ``self``
        is mapped to `-1` if it is an idescent and `1` if it is not an
        idescent.

        See :meth:`idescents` for a definition of idescents.

        With the ``final_descent`` option, the last position of a
        non-empty permutation is also considered as a descent.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).idescents()
            [2, 3]
            sage: Permutation([1,4,3,2]).idescents_signature()
            [1, -1, -1, 1]
        """
        idescents = self.idescents(final_descent=final_descent)
        d = {True: -1, False: 1}
        return [d[(i + 1) in idescents] for i in range(len(self))]

    def number_of_descents(self, final_descent=False) -> Integer:
        r"""
        Return the number of descents of ``self``.

        With the ``final_descent`` option, the last position of a
        non-empty permutation is also considered as a descent.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).number_of_descents()
            2
            sage: Permutation([1,4,3,2]).number_of_descents(final_descent=True)
            3
        """
        return len(self.descents(final_descent))

    def number_of_idescents(self, final_descent=False) -> Integer:
        r"""
        Return the number of idescents of ``self``.

        See :meth:`idescents` for a definition of idescents.

        With the ``final_descent`` option, the last position of a
        non-empty permutation is also considered as a descent.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).number_of_idescents()
            2
            sage: Permutation([1,4,3,2]).number_of_idescents(final_descent=True)
            3
        """
        return len(self.idescents(final_descent))

    @combinatorial_map(name='descent composition')
    def descents_composition(self) -> Composition:
        r"""
        Return the descent composition of ``self``.

        The descent composition of a permutation `p \in S_n` is defined
        as the composition of `n` whose descent set equals the descent
        set of `p`. Here, the descent set of `p` is defined as the set
        of all `i \in \{ 1, 2, \ldots, n-1 \}` satisfying
        `p(i) > p(i+1)`. The descent set
        of a composition `c = (i_1, i_2, \ldots, i_k)` is defined as
        the set `\{ i_1, i_1 + i_2, i_1 + i_2 + i_3, \ldots,
        i_1 + i_2 + \cdots + i_{k-1} \}`.

        EXAMPLES::

            sage: Permutation([1,3,2,4]).descents_composition()
            [2, 2]
            sage: Permutation([4,1,6,7,2,3,8,5]).descents_composition()
            [1, 3, 3, 1]
            sage: Permutation([]).descents_composition()
            []
        """
        if len(self) == 0:
            return Composition([])
        d = [0] + self.descents(final_descent=True)
        return Composition([d[i + 1] - d[i] for i in range(len(d) - 1)])

    def descent_polynomial(self):
        r"""
        Return the descent polynomial of the permutation ``self``.

        The descent polynomial of a permutation `p` is the product of
        all the ``z[p(i)]`` where ``i`` ranges over the descents of
        ``p``.

        A descent of a permutation ``p`` is an integer ``i`` such that
        ``p(i) > p(i+1)``.

        REFERENCES:

        - [GS1984]_

        EXAMPLES::

            sage: Permutation([2,1,3]).descent_polynomial()
            z1
            sage: Permutation([4,3,2,1]).descent_polynomial()
            z1*z2^2*z3^3

        .. TODO::

            This docstring needs to be fixed. First, the definition
            does not match the implementation (or the examples).
            Second, this doesn't seem to be defined in [GS1984]_
            (the descent monomial in their (7.23) is different).
        """
        p = self
        z = []
        P = PolynomialRing(ZZ, len(p), 'z')
        z = P.gens()
        result = 1
        pol = 1
        for i in range(len(p)-1):
            pol *= z[p[i]-1]
            if p[i] > p[i+1]:
                result *= pol

        return result

    ##############
    # Major Code #
    ##############

    def major_index(self, final_descent=False) -> Integer:
        r"""
        Return the major index of ``self``.

        The major index of a permutation `p` is the sum of the descents of `p`.
        Since our permutation indices are 0-based, we need to add the
        number of descents.

        With the ``final_descent`` option, the last position of a
        non-empty permutation is also considered as a descent.

        EXAMPLES::

            sage: Permutation([2,1,3]).major_index()
            1
            sage: Permutation([3,4,1,2]).major_index()
            2
            sage: Permutation([4,3,2,1]).major_index()
            6
        """
        descents = self.descents(final_descent)
        return sum(descents)

    def imajor_index(self, final_descent=False) -> Integer:
        """
        Return the inverse major index of the permutation ``self``, which is
        the major index of the inverse of ``self``.

        The major index of a permutation `p` is the sum of the descents of `p`.
        Since our permutation indices are 0-based, we need to add the
        number of descents.

        With the ``final_descent`` option, the last position of a
        non-empty permutation is also considered as a descent.

        EXAMPLES::

            sage: Permutation([2,1,3]).imajor_index()
            1
            sage: Permutation([3,4,1,2]).imajor_index()
            2
            sage: Permutation([4,3,2,1]).imajor_index()
            6
        """
        return sum(self.idescents(final_descent))

    def to_major_code(self, final_descent=False):
        r"""
        Return the major code of the permutation ``self``.

        The major code of a permutation `p` is defined as the sequence
        `(m_1-m_2, m_2-m_3, \ldots, m_n)`, where `m_i` is the major
        index of the permutation obtained by erasing all letters smaller than
        `i` from `p`.

        With the ``final_descent`` option, the last position of a
        non-empty permutation is also considered as a descent.
        This has an effect on the computation of major indices.

        REFERENCES:

        - Carlitz, L. *q-Bernoulli and Eulerian Numbers*.
          Trans. Amer. Math. Soc. 76 (1954) 332-350.
          http://www.ams.org/journals/tran/1954-076-02/S0002-9947-1954-0060538-2/

        - Skandera, M. *An Eulerian Partner for Inversions*.
          Sém. Lothar. Combin. 46 (2001) B46d.
          http://www.lehigh.edu/~mas906/papers/partner.ps

        EXAMPLES::

            sage: Permutation([9,3,5,7,2,1,4,6,8]).to_major_code()
            [5, 0, 1, 0, 1, 2, 0, 1, 0]
            sage: Permutation([2,8,4,3,6,7,9,5,1]).to_major_code()
            [8, 3, 3, 1, 4, 0, 1, 0, 0]
        """
        p = self
        n = len(p)
        major_indices = [0]*(n+1)
        smaller = p[:]
        P = Permutations()
        for i in range(n):
            major_indices[i] = P(smaller).major_index(final_descent)
            #Create the permutation that "erases" all the numbers
            #smaller than i+1
            smaller.remove(1)
            smaller = [i-1 for i in smaller]

        major_code = [ major_indices[i] - major_indices[i+1] for i in range(n) ]
        return major_code

    #########
    # Peaks #
    #########

    def peaks(self) -> list:
        r"""
        Return a list of the peaks of the permutation ``self``.

        A peak of a permutation `p` is an integer `i` such that
        `p(i-1) < p(i)` and `p(i) > p(i+1)`.

        EXAMPLES::

            sage: Permutation([1,3,2,4,5]).peaks()
            [1]
            sage: Permutation([4,1,3,2,6,5]).peaks()
            [2, 4]
            sage: Permutation([]).peaks()
            []
        """
        p = self
        peaks = []
        for i in range(1,len(p)-1):
            if p[i-1] <= p[i] and p[i] > p[i+1]:
                peaks.append(i)

        return peaks

    def number_of_peaks(self) -> Integer:
        r"""
        Return the number of peaks of the permutation ``self``.

        A peak of a permutation `p` is an integer `i` such that
        `p(i-1) < p(i)` and `p(i) > p(i+1)`.

        EXAMPLES::

            sage: Permutation([1,3,2,4,5]).number_of_peaks()
            1
            sage: Permutation([4,1,3,2,6,5]).number_of_peaks()
            2
        """
        return len(self.peaks())

    #############
    # Saliances #
    #############

    def saliances(self) -> list:
        r"""
        Return a list of the saliances of the permutation ``self``.

        A saliance of a permutation `p` is an integer `i` such that
        `p(i) > p(j)` for all `j > i`.

        EXAMPLES::

            sage: Permutation([2,3,1,5,4]).saliances()
            [3, 4]
            sage: Permutation([5,4,3,2,1]).saliances()
            [0, 1, 2, 3, 4]
        """
        p = self
        saliances = []
        for i in range(len(p)):
            is_saliance = True
            for j in range(i+1, len(p)):
                if p[i] <= p[j]:
                    is_saliance = False
            if is_saliance:
                saliances.append(i)

        return saliances

    def number_of_saliances(self) -> Integer:
        r"""
        Return the number of saliances of ``self``.

        A saliance of a permutation `p` is an integer `i` such that
        `p(i) > p(j)` for all `j > i`.

        EXAMPLES::

            sage: Permutation([2,3,1,5,4]).number_of_saliances()
            2
            sage: Permutation([5,4,3,2,1]).number_of_saliances()
            5
        """
        return len(self.saliances())

    ################
    # Bruhat Order #
    ################

    def bruhat_lequal(self, p2) -> bool:
        r"""
        Return ``True`` if ``self`` is less or equal to ``p2`` in
        the Bruhat order.

        The Bruhat order (also called strong Bruhat order or Chevalley
        order) on the symmetric group `S_n` is the partial order on `S_n`
        determined by the following condition: If `p` is a permutation,
        and `i` and `j` are two indices satisfying `p(i) > p(j)` and
        `i < j` (that is, `(i, j)` is an inversion of `p` with `i < j`),
        then `p \circ (i, j)` (the permutation obtained by first
        switching `i` with `j` and then applying `p`) is smaller than `p`
        in the Bruhat order.

        One can show that a permutation `p \in S_n` is less or equal to
        a permutation `q \in S_n` in the Bruhat order if and only if
        for every `i \in \{ 0, 1, \cdots , n \}` and
        `j \in \{ 1, 2, \cdots , n \}`, the number of the elements among
        `p(1), p(2), \cdots, p(j)` that are greater than `i` is `\leq`
        to the number of the elements among `q(1), q(2), \cdots, q(j)`
        that are greater than `i`.

        This method assumes that ``self`` and ``p2`` are permutations
        of the same integer `n`.

        EXAMPLES::

            sage: Permutation([2,4,3,1]).bruhat_lequal(Permutation([3,4,2,1]))
            True

            sage: Permutation([2,1,3]).bruhat_lequal(Permutation([2,3,1]))
            True
            sage: Permutation([2,1,3]).bruhat_lequal(Permutation([3,1,2]))
            True
            sage: Permutation([2,1,3]).bruhat_lequal(Permutation([1,2,3]))
            False
            sage: Permutation([1,3,2]).bruhat_lequal(Permutation([2,1,3]))
            False
            sage: Permutation([1,3,2]).bruhat_lequal(Permutation([2,3,1]))
            True
            sage: Permutation([2,3,1]).bruhat_lequal(Permutation([1,3,2]))
            False
            sage: sorted( [len([b for b in Permutations(3) if a.bruhat_lequal(b)])
            ....:          for a in Permutations(3)] )
            [1, 2, 2, 4, 4, 6]

            sage: Permutation([]).bruhat_lequal(Permutation([]))
            True
        """
        p1 = self
        n1 = len(p1)

        if n1 == 0:
            return True

        if p1[0] > p2[0] or p1[n1-1] < p2[n1-1]:
            return False

        for i in range(1, n1):
            c = 0
            for j in range(n1 - 2):
                # We should really check this for all j in range(n1), but for
                # j == n1 - 1 it is tautological and for j == n1 - 2 the check
                # is contained in the check p1[n1-1] < p2[n1-1] already made.
                if p2[j] > i:
                    c += 1
                if p1[j] > i:
                    c -= 1
                if c < 0:
                    return False

        return True

    def weak_excedences(self) -> list:
        """
        Return all the numbers ``self[i]`` such that ``self[i] >= i+1``.

        EXAMPLES::

            sage: Permutation([1,4,3,2,5]).weak_excedences()
            [1, 4, 3, 5]
        """
        res = []
        for i in range(len(self)):
            if self[i] >= i + 1:
                res.append(self[i])
        return res

    def bruhat_inversions(self) -> list:
        r"""
        Return the list of inversions of ``self`` such that the application of
        this inversion to ``self`` decreases its number of inversions by
        exactly 1.

        Equivalently, it returns the list of pairs `(i,j)` such that `i < j`,
        such that `p(i) > p(j)` and such that there exists no `k` (strictly)
        between `i` and `j` satisfying `p(i) > p(k) > p(j)`.

        EXAMPLES::

            sage: Permutation([5,2,3,4,1]).bruhat_inversions()
            [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
            sage: Permutation([6,1,4,5,2,3]).bruhat_inversions()
            [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [3, 4], [3, 5]]
        """
        return list(self.bruhat_inversions_iterator())

    def bruhat_inversions_iterator(self):
        """
        Return the iterator for the inversions of ``self`` such that the
        application of this inversion to ``self`` decreases its number of
        inversions by exactly 1.

        EXAMPLES::

            sage: list(Permutation([5,2,3,4,1]).bruhat_inversions_iterator())
            [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
            sage: list(Permutation([6,1,4,5,2,3]).bruhat_inversions_iterator())
            [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [3, 4], [3, 5]]
        """
        p = self
        n = len(p)

        for i in range(n-1):
            for j in range(i+1,n):
                if p[i] > p[j]:
                    ok = True
                    for k in range(i+1, j):
                        if p[i] > p[k] and p[k] > p[j]:
                            ok = False
                            break
                    if ok:
                        yield [i,j]

    def bruhat_succ(self) -> list:
        r"""
        Return a list of the permutations strictly greater than ``self`` in
        the Bruhat order (on the symmetric group containing ``self``) such
        that there is no permutation between one of those and ``self``.

        See :meth:`bruhat_lequal` for the definition of the Bruhat order.

        EXAMPLES::

            sage: Permutation([6,1,4,5,2,3]).bruhat_succ()
            [[6, 4, 1, 5, 2, 3],
             [6, 2, 4, 5, 1, 3],
             [6, 1, 5, 4, 2, 3],
             [6, 1, 4, 5, 3, 2]]
        """
        return list(self.bruhat_succ_iterator())

    def bruhat_succ_iterator(self):
        """
        An iterator for the permutations that are strictly greater than
        ``self`` in the Bruhat order (on the symmetric group containing
        ``self``) such that there is no permutation between one
        of those and ``self``.

        See :meth:`bruhat_lequal` for the definition of the Bruhat order.

        EXAMPLES::

            sage: [x for x in Permutation([6,1,4,5,2,3]).bruhat_succ_iterator()]
            [[6, 4, 1, 5, 2, 3],
             [6, 2, 4, 5, 1, 3],
             [6, 1, 5, 4, 2, 3],
             [6, 1, 4, 5, 3, 2]]
        """
        p = self
        n = len(p)
        P = Permutations()

        for z in P([n+1-x for x in p]).bruhat_inversions_iterator():
            pp = p[:]
            pp[z[0]] = p[z[1]]
            pp[z[1]] = p[z[0]]
            yield P(pp)

    def bruhat_pred(self) -> list:
        r"""
        Return a list of the permutations strictly smaller than ``self``
        in the Bruhat order (on the symmetric group containing ``self``) such
        that there is no permutation between one of those and ``self``.

        See :meth:`bruhat_lequal` for the definition of the Bruhat order.

        EXAMPLES::

            sage: Permutation([6,1,4,5,2,3]).bruhat_pred()
            [[1, 6, 4, 5, 2, 3],
             [4, 1, 6, 5, 2, 3],
             [5, 1, 4, 6, 2, 3],
             [6, 1, 2, 5, 4, 3],
             [6, 1, 3, 5, 2, 4],
             [6, 1, 4, 2, 5, 3],
             [6, 1, 4, 3, 2, 5]]
        """
        return list(self.bruhat_pred_iterator())

    def bruhat_pred_iterator(self):
        """
        An iterator for the permutations strictly smaller than ``self`` in
        the Bruhat order (on the symmetric group containing ``self``) such
        that there is no permutation between one of those and ``self``.

        See :meth:`bruhat_lequal` for the definition of the Bruhat order.

        EXAMPLES::

            sage: [x for x in Permutation([6,1,4,5,2,3]).bruhat_pred_iterator()]
            [[1, 6, 4, 5, 2, 3],
             [4, 1, 6, 5, 2, 3],
             [5, 1, 4, 6, 2, 3],
             [6, 1, 2, 5, 4, 3],
             [6, 1, 3, 5, 2, 4],
             [6, 1, 4, 2, 5, 3],
             [6, 1, 4, 3, 2, 5]]
        """
        p = self
        P = Permutations()
        for z in p.bruhat_inversions_iterator():
            pp = p[:]
            pp[z[0]] = p[z[1]]
            pp[z[1]] = p[z[0]]
            yield P(pp)

    def bruhat_smaller(self):
        r"""
        Return the combinatorial class of permutations smaller than or
        equal to ``self`` in the Bruhat order (on the symmetric group
        containing ``self``).

        See :meth:`bruhat_lequal` for the definition of the Bruhat order.

        EXAMPLES::

            sage: Permutation([4,1,2,3]).bruhat_smaller().list()
            [[1, 2, 3, 4],
             [1, 2, 4, 3],
             [1, 3, 2, 4],
             [1, 4, 2, 3],
             [2, 1, 3, 4],
             [2, 1, 4, 3],
             [3, 1, 2, 4],
             [4, 1, 2, 3]]
        """
        return StandardPermutations_bruhat_smaller(self)

    def bruhat_greater(self):
        r"""
        Return the combinatorial class of permutations greater than or
        equal to ``self`` in the Bruhat order (on the symmetric group
        containing ``self``).

        See :meth:`bruhat_lequal` for the definition of the Bruhat order.

        EXAMPLES::

            sage: Permutation([4,1,2,3]).bruhat_greater().list()
            [[4, 1, 2, 3],
             [4, 1, 3, 2],
             [4, 2, 1, 3],
             [4, 2, 3, 1],
             [4, 3, 1, 2],
             [4, 3, 2, 1]]
        """
        return StandardPermutations_bruhat_greater(self)

    ########################
    # Permutohedron  Order #
    ########################

    def permutohedron_lequal(self, p2, side="right") -> bool:
        r"""
        Return ``True`` if ``self`` is less or equal to ``p2`` in the
        permutohedron order.

        By default, the computations are done in the right permutohedron.
        If you pass the option ``side='left'``, then they will be done in
        the left permutohedron.

        For every nonnegative integer `n`, the right (resp. left)
        permutohedron order (also called the right (resp. left) weak
        order, or the right (resp. left) weak Bruhat order) is a partial
        order on the symmetric group `S_n`. It can be defined in various
        ways, including the following ones:

        - Two permutations `u` and `v` in `S_n` satisfy `u \leq v` in
          the right (resp. left) permutohedron order if and only if
          the (Coxeter) length of the permutation `v^{-1} \circ u`
          (resp. of the permutation `u \circ v^{-1}`) equals the
          length of `v` minus the length of `u`. Here, `p \circ q` means
          the permutation obtained by applying `q` first and then `p`.
          (Recall that the Coxeter length of a permutation is its number
          of inversions.)

        - Two permutations `u` and `v` in `S_n` satisfy `u \leq v` in
          the right (resp. left) permutohedron order if and only if
          every pair `(i, j)` of elements of `\{ 1, 2, \cdots, n \}`
          such that `i < j` and `u^{-1}(i) > u^{-1}(j)` (resp.
          `u(i) > u(j)`) also satisfies `v^{-1}(i) > v^{-1}(j)`
          (resp. `v(i) > v(j)`).

        - A permutation `v \in S_n` covers a permutation `u \in S_n` in
          the right (resp. left) permutohedron order if and only if we
          have `v = u \circ (i, i + 1)` (resp. `v = (i, i + 1) \circ u`)
          for some `i \in \{ 1, 2, \cdots, n - 1 \}` satisfying
          `u(i) < u(i + 1)` (resp. `u^{-1}(i) < u^{-1}(i + 1)`). Here,
          again, `p \circ q` means the permutation obtained by applying
          `q` first and then `p`.

        The right and the left permutohedron order are mutually
        isomorphic, with the isomorphism being the map sending every
        permutation to its inverse. Each of these orders endows the
        symmetric group `S_n` with the structure of a graded poset
        (the rank function being the Coxeter length).

        .. WARNING::

            The permutohedron order is not to be mistaken for the
            strong Bruhat order (:meth:`bruhat_lequal`), despite both
            orders being occasionally referred to as the Bruhat order.

        EXAMPLES::

            sage: p = Permutation([3,2,1,4])
            sage: p.permutohedron_lequal(Permutation([4,2,1,3]))
            False
            sage: p.permutohedron_lequal(Permutation([4,2,1,3]), side='left')
            True
            sage: p.permutohedron_lequal(p)
            True

            sage: Permutation([2,1,3]).permutohedron_lequal(Permutation([2,3,1]))
            True
            sage: Permutation([2,1,3]).permutohedron_lequal(Permutation([3,1,2]))
            False
            sage: Permutation([2,1,3]).permutohedron_lequal(Permutation([1,2,3]))
            False
            sage: Permutation([1,3,2]).permutohedron_lequal(Permutation([2,1,3]))
            False
            sage: Permutation([1,3,2]).permutohedron_lequal(Permutation([2,3,1]))
            False
            sage: Permutation([2,3,1]).permutohedron_lequal(Permutation([1,3,2]))
            False
            sage: Permutation([2,1,3]).permutohedron_lequal(Permutation([2,3,1]), side='left')
            False
            sage: sorted( [len([b for b in Permutations(3) if a.permutohedron_lequal(b)])
            ....:          for a in Permutations(3)] )
            [1, 2, 2, 3, 3, 6]
            sage: sorted( [len([b for b in Permutations(3) if a.permutohedron_lequal(b, side="left")])
            ....:          for a in Permutations(3)] )
            [1, 2, 2, 3, 3, 6]

            sage: Permutation([]).permutohedron_lequal(Permutation([]))
            True
        """
        p1 = self
        l1 = p1.number_of_inversions()
        l2 = p2.number_of_inversions()

        if l1 > l2:
            return False

        if side == "right":
            prod = p1._left_to_right_multiply_on_right(p2.inverse())
        else:
            prod = p1._left_to_right_multiply_on_left(p2.inverse())

        return prod.number_of_inversions() == l2 - l1

    def permutohedron_succ(self, side="right"):
        r"""
        Return a list of the permutations strictly greater than ``self``
        in the permutohedron order such that there is no permutation
        between any of those and ``self``.

        By default, the computations are done in the right permutohedron.
        If you pass the option ``side='left'``, then they will be done in
        the left permutohedron.

        See :meth:`permutohedron_lequal` for the definition of the
        permutohedron orders.

        EXAMPLES::

            sage: p = Permutation([4,2,1,3])
            sage: p.permutohedron_succ()
            [[4, 2, 3, 1]]
            sage: p.permutohedron_succ(side='left')
            [[4, 3, 1, 2]]
        """
        p = self
        n = len(p)
        P = Permutations()
        succ = []
        if side == "right":
            rise = lambda perm: [i for i in range(n - 1) if perm[i] < perm[i+1]]
            for i in rise(p):
                pp = p[:]
                pp[i] = p[i+1]
                pp[i+1] = p[i]
                succ.append(P(pp))
        else:
            advance = lambda perm: [i for i in range(1,n) if perm.index(i) < perm.index(i+1)]
            for i in advance(p):
                pp = p[:]
                pp[p.index(i)] = i+1
                pp[p.index(i+1)] = i
                succ.append(P(pp))
        return succ

    def permutohedron_pred(self, side="right") -> list:
        r"""
        Return a list of the permutations strictly smaller than ``self``
        in the permutohedron order such that there is no permutation
        between any of those and ``self``.

        By default, the computations are done in the right permutohedron.
        If you pass the option ``side='left'``, then they will be done in
        the left permutohedron.

        See :meth:`permutohedron_lequal` for the definition of the
        permutohedron orders.

        EXAMPLES::

            sage: p = Permutation([4,2,1,3])
            sage: p.permutohedron_pred()
            [[2, 4, 1, 3], [4, 1, 2, 3]]
            sage: p.permutohedron_pred(side='left')
            [[4, 1, 2, 3], [3, 2, 1, 4]]
        """
        p = self
        n = len(p)
        P = Permutations()
        pred = []
        if side == "right":
            for d in p.descents():
                pp = p[:]
                pp[d - 1] = p[d]
                pp[d] = p[d - 1]
                pred.append(P(pp))
        else:
            recoil = lambda perm: [i for i in range(1,n) if perm.index(i) > perm.index(i+1)]
            for i in recoil(p):
                pp = p[:]
                pp[p.index(i)] = i+1
                pp[p.index(i+1)] = i
                pred.append(P(pp))
        return pred

    def permutohedron_smaller(self, side="right") -> list:
        r"""
        Return a list of permutations smaller than or equal to ``self``
        in the permutohedron order.

        By default, the computations are done in the right permutohedron.
        If you pass the option ``side='left'``, then they will be done in
        the left permutohedron.

        See :meth:`permutohedron_lequal` for the definition of the
        permutohedron orders.

        EXAMPLES::

            sage: Permutation([4,2,1,3]).permutohedron_smaller()
            [[1, 2, 3, 4],
             [1, 2, 4, 3],
             [1, 4, 2, 3],
             [2, 1, 3, 4],
             [2, 1, 4, 3],
             [2, 4, 1, 3],
             [4, 1, 2, 3],
             [4, 2, 1, 3]]

        ::

            sage: Permutation([4,2,1,3]).permutohedron_smaller(side='left')
            [[1, 2, 3, 4],
             [1, 3, 2, 4],
             [2, 1, 3, 4],
             [2, 3, 1, 4],
             [3, 1, 2, 4],
             [3, 2, 1, 4],
             [4, 1, 2, 3],
             [4, 2, 1, 3]]
        """
        return transitive_ideal(lambda x: x.permutohedron_pred(side), self)

    def permutohedron_greater(self, side="right") -> list:
        r"""
        Return a list of permutations greater than or equal to ``self``
        in the permutohedron order.

        By default, the computations are done in the right permutohedron.
        If you pass the option ``side='left'``, then they will be done in
        the left permutohedron.

        See :meth:`permutohedron_lequal` for the definition of the
        permutohedron orders.

        EXAMPLES::

            sage: Permutation([4,2,1,3]).permutohedron_greater()
            [[4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 2, 1]]
            sage: Permutation([4,2,1,3]).permutohedron_greater(side='left')
            [[4, 2, 1, 3], [4, 3, 1, 2], [4, 3, 2, 1]]
        """
        return transitive_ideal(lambda x: x.permutohedron_succ(side), self)

    def right_permutohedron_interval_iterator(self, other):
        r"""
        Return an iterator on the permutations (represented as integer
        lists) belonging to the right permutohedron interval where
        ``self`` is the minimal element and ``other`` the maximal element.

        See :meth:`permutohedron_lequal` for the definition of the
        permutohedron orders.

        EXAMPLES::

            sage: Permutation([2, 1, 4, 5, 3]).right_permutohedron_interval(Permutation([2, 5, 4, 1, 3])) # indirect doctest
            [[2, 4, 5, 1, 3], [2, 4, 1, 5, 3], [2, 1, 4, 5, 3], [2, 1, 5, 4, 3], [2, 5, 1, 4, 3], [2, 5, 4, 1, 3]]
        """
        if len(self) != len(other):
            raise ValueError("len({}) and len({}) must be equal".format(self, other))
        if not self.permutohedron_lequal(other):
            raise ValueError("{} must be lower or equal than {} for the right permutohedron order".format(self, other))
        d = DiGraph()
        d.add_vertices(range(1, len(self) + 1))
        d.add_edges([(j, i) for (i, j) in self.inverse().inversions()])
        d.add_edges([(other[i], other[j]) for i in range(len(other) - 1)
                     for j in range(i, len(other)) if other[i] < other[j]])
        return d.topological_sort_generator()

    def right_permutohedron_interval(self, other):
        r"""
        Return the list of the permutations belonging to the right
        permutohedron interval where ``self`` is the minimal element and
        ``other`` the maximal element.

        See :meth:`permutohedron_lequal` for the definition of the
        permutohedron orders.

        EXAMPLES::

            sage: Permutation([2, 1, 4, 5, 3]).right_permutohedron_interval(Permutation([2, 5, 4, 1, 3]))
            [[2, 4, 5, 1, 3], [2, 4, 1, 5, 3], [2, 1, 4, 5, 3], [2, 1, 5, 4, 3], [2, 5, 1, 4, 3], [2, 5, 4, 1, 3]]

        TESTS::

            sage: Permutation([]).right_permutohedron_interval(Permutation([]))
            [[]]
            sage: Permutation([3, 1, 2]).right_permutohedron_interval(Permutation([3, 1, 2]))
            [[3, 1, 2]]
            sage: Permutation([1, 3, 2, 4]).right_permutohedron_interval(Permutation([3, 4, 2, 1]))
            [[3, 1, 4, 2], [3, 4, 1, 2], [3, 4, 2, 1], [1, 3, 4, 2], [1, 3, 2, 4], [3, 2, 4, 1], [3, 2, 1, 4], [3, 1, 2, 4]]
            sage: Permutation([2, 1, 4, 5, 3]).right_permutohedron_interval(Permutation([2, 5, 4, 1, 3]))
            [[2, 4, 5, 1, 3], [2, 4, 1, 5, 3], [2, 1, 4, 5, 3], [2, 1, 5, 4, 3], [2, 5, 1, 4, 3], [2, 5, 4, 1, 3]]
            sage: Permutation([2, 5, 4, 1, 3]).right_permutohedron_interval(Permutation([2, 1, 4, 5, 3]))
            Traceback (most recent call last):
            ...
            ValueError: [2, 5, 4, 1, 3] must be lower or equal than [2, 1, 4, 5, 3] for the right permutohedron order
            sage: Permutation([2, 4, 1, 3]).right_permutohedron_interval(Permutation([2, 1, 4, 5, 3]))
            Traceback (most recent call last):
            ...
            ValueError: len([2, 4, 1, 3]) and len([2, 1, 4, 5, 3]) must be equal
        """
        P = Permutations()
        return [P(p) for p in self.right_permutohedron_interval_iterator(other)]

    def permutohedron_join(self, other, side="right") -> Permutation:
        r"""
        Return the join of the permutations ``self`` and ``other``
        in the right permutohedron order (or, if ``side`` is set to
        ``'left'``, in the left permutohedron order).

        The permutohedron orders (see :meth:`permutohedron_lequal`)
        are lattices; the join operation refers to this lattice
        structure. In more elementary terms, the join of two
        permutations `\pi` and `\psi` in the symmetric group `S_n`
        is the permutation in `S_n` whose set of inversion is the
        transitive closure of the union of the set of inversions of
        `\pi` with the set of inversions of `\psi`.

        .. SEEALSO::

            :meth:`permutohedron_lequal`, :meth:`permutohedron_meet`.

        ALGORITHM:

        It is enough to construct the join of any two permutations
        `\pi` and `\psi` in `S_n` with respect to the right weak
        order. (The join of `\pi` and `\psi` with respect to the
        left weak order is the inverse of the join of `\pi^{-1}`
        and `\psi^{-1}` with respect to the right weak order.)
        Start with an empty list `l` (denoted ``xs`` in the actual
        code). For `i = 1, 2, \ldots, n` (in this order), we insert
        `i` into this list in the rightmost possible position such
        that any letter in `\{ 1, 2, ..., i-1 \}` which appears
        further right than `i` in either `\pi` or `\psi` (or both)
        must appear further right than `i` in the resulting list.
        After all numbers are inserted, we are left with a list
        which is precisely the join of `\pi` and `\psi` (in
        one-line notation). This algorithm is due to Markowsky,
        [Mar1994]_ (Theorem 1 (a)).

        AUTHORS:

        Viviane Pons and Darij Grinberg, 18 June 2014.

        EXAMPLES::

            sage: p = Permutation([3,1,2])
            sage: q = Permutation([1,3,2])
            sage: p.permutohedron_join(q)
            [3, 1, 2]
            sage: r = Permutation([2,1,3])
            sage: r.permutohedron_join(p)
            [3, 2, 1]

        ::

            sage: p = Permutation([3,2,4,1])
            sage: q = Permutation([4,2,1,3])
            sage: p.permutohedron_join(q)
            [4, 3, 2, 1]
            sage: r = Permutation([3,1,2,4])
            sage: p.permutohedron_join(r)
            [3, 2, 4, 1]
            sage: q.permutohedron_join(r)
            [4, 3, 2, 1]
            sage: s = Permutation([1,4,2,3])
            sage: s.permutohedron_join(r)
            [4, 3, 1, 2]

        The universal property of the join operation is
        satisfied::

            sage: def test_uni_join(p, q):
            ....:     j = p.permutohedron_join(q)
            ....:     if not p.permutohedron_lequal(j):
            ....:         return False
            ....:     if not q.permutohedron_lequal(j):
            ....:         return False
            ....:     for r in p.permutohedron_greater():
            ....:         if q.permutohedron_lequal(r) and not j.permutohedron_lequal(r):
            ....:             return False
            ....:     return True
            sage: all( test_uni_join(p, q) for p in Permutations(3) for q in Permutations(3) )
            True
            sage: test_uni_join(Permutation([6, 4, 7, 3, 2, 5, 8, 1]), Permutation([7, 3, 1, 2, 5, 4, 6, 8]))
            True

        Border cases::

            sage: p = Permutation([])
            sage: p.permutohedron_join(p)
            []
            sage: p = Permutation([1])
            sage: p.permutohedron_join(p)
            [1]

        The left permutohedron:

            sage: p = Permutation([3,1,2])
            sage: q = Permutation([1,3,2])
            sage: p.permutohedron_join(q, side="left")
            [3, 2, 1]
            sage: r = Permutation([2,1,3])
            sage: r.permutohedron_join(p, side="left")
            [3, 1, 2]
        """
        if side == "left":
            return self.inverse().permutohedron_join(other.inverse()).inverse()
        n = self.size()
        xs = []
        for i in range(1, n + 1):
            u = self.index(i)
            must_be_right = [f for f in self[u + 1:] if f < i]
            v = other.index(i)
            must_be_right += [f for f in other[v + 1:] if f < i]
            must_be_right = sorted(set(must_be_right))
            for j, q in enumerate(xs):
                if q in must_be_right:
                    xs = xs[:j] + [i] + xs[j:]
                    break
            else:
                xs.append(i)
        return Permutations(n)(xs)

    def permutohedron_meet(self, other, side="right") -> Permutation:
        r"""
        Return the meet of the permutations ``self`` and ``other``
        in the right permutohedron order (or, if ``side`` is set to
        ``'left'``, in the left permutohedron order).

        The permutohedron orders (see :meth:`permutohedron_lequal`)
        are lattices; the meet operation refers to this lattice
        structure. It is connected to the join operation by the
        following simple symmetry property: If `\pi` and `\psi`
        are two permutations `\pi` and `\psi` in the symmetric group
        `S_n`, and if `w_0` denotes the permutation
        `(n, n-1, \ldots, 1) \in S_n`, then

        .. MATH::

            \pi \wedge \psi = w_0 \circ ((w_0 \circ \pi) \vee (w_0 \circ \psi))
            = ((\pi \circ w_0) \vee (\psi \circ w_0)) \circ w_0

        and

        .. MATH::

            \pi \vee \psi = w_0 \circ ((w_0 \circ \pi) \wedge (w_0 \circ \psi))
            = ((\pi \circ w_0) \wedge (\psi \circ w_0)) \circ w_0,

        where `\wedge` means meet and `\vee` means join.

        .. SEEALSO::

            :meth:`permutohedron_lequal`, :meth:`permutohedron_join`.

        AUTHORS:

        Viviane Pons and Darij Grinberg, 18 June 2014.

        EXAMPLES::

            sage: p = Permutation([3,1,2])
            sage: q = Permutation([1,3,2])
            sage: p.permutohedron_meet(q)
            [1, 3, 2]
            sage: r = Permutation([2,1,3])
            sage: r.permutohedron_meet(p)
            [1, 2, 3]

        ::

            sage: p = Permutation([3,2,4,1])
            sage: q = Permutation([4,2,1,3])
            sage: p.permutohedron_meet(q)
            [2, 1, 3, 4]
            sage: r = Permutation([3,1,2,4])
            sage: p.permutohedron_meet(r)
            [3, 1, 2, 4]
            sage: q.permutohedron_meet(r)
            [1, 2, 3, 4]
            sage: s = Permutation([1,4,2,3])
            sage: s.permutohedron_meet(r)
            [1, 2, 3, 4]

        The universal property of the meet operation is
        satisfied::

            sage: def test_uni_meet(p, q):
            ....:     m = p.permutohedron_meet(q)
            ....:     if not m.permutohedron_lequal(p):
            ....:         return False
            ....:     if not m.permutohedron_lequal(q):
            ....:         return False
            ....:     for r in p.permutohedron_smaller():
            ....:         if r.permutohedron_lequal(q) and not r.permutohedron_lequal(m):
            ....:             return False
            ....:     return True
            sage: all( test_uni_meet(p, q) for p in Permutations(3) for q in Permutations(3) )
            True
            sage: test_uni_meet(Permutation([6, 4, 7, 3, 2, 5, 8, 1]), Permutation([7, 3, 1, 2, 5, 4, 6, 8]))
            True

        Border cases::

            sage: p = Permutation([])
            sage: p.permutohedron_meet(p)
            []
            sage: p = Permutation([1])
            sage: p.permutohedron_meet(p)
            [1]

        The left permutohedron:

            sage: p = Permutation([3,1,2])
            sage: q = Permutation([1,3,2])
            sage: p.permutohedron_meet(q, side="left")
            [1, 2, 3]
            sage: r = Permutation([2,1,3])
            sage: r.permutohedron_meet(p, side="left")
            [2, 1, 3]
        """
        return self.reverse().permutohedron_join(other.reverse(), side=side).reverse()

    ############
    # Patterns #
    ############

    def has_pattern(self, patt) -> bool:
        r"""
        Test whether the permutation ``self`` contains the pattern
        ``patt``.

        EXAMPLES::

            sage: Permutation([3,5,1,4,6,2]).has_pattern([1,3,2])
            True
        """
        p = self
        n = len(p)
        l = len(patt)
        if l > n:
            return False
        for pos in itertools.combinations(range(n), l):
            if to_standard([p[z] for z in pos]) == patt:
                return True
        return False

    def avoids(self, patt) -> bool:
        """
        Test whether the permutation ``self`` avoids the pattern
        ``patt``.

        EXAMPLES::

            sage: Permutation([6,2,5,4,3,1]).avoids([4,2,3,1])
            False
            sage: Permutation([6,1,2,5,4,3]).avoids([4,2,3,1])
            True
            sage: Permutation([6,1,2,5,4,3]).avoids([3,4,1,2])
            True
        """
        return not self.has_pattern(patt)

    def pattern_positions(self, patt) -> list:
        r"""
        Return the list of positions where the pattern ``patt`` appears
        in the permutation ``self``.

        EXAMPLES::

            sage: Permutation([3,5,1,4,6,2]).pattern_positions([1,3,2])
            [[0, 1, 3], [2, 3, 5], [2, 4, 5]]
        """
        p = self

        return [list(pos) for pos in itertools.combinations(range(len(p)), len(patt))
                if to_standard([p[z] for z in pos]) == patt]

    @combinatorial_map(name='Simion-Schmidt map')
    def simion_schmidt(self, avoid=[1,2,3]):
        r"""
        Implements the Simion-Schmidt map which sends an arbitrary permutation
        to a pattern avoiding permutation, where the permutation pattern is one
        of four length-three patterns.  This method also implements the bijection
        between (for example) ``[1,2,3]``- and ``[1,3,2]``-avoiding permutations.

        INPUT:

        - ``avoid`` -- one of the patterns ``[1,2,3]``, ``[1,3,2]``, ``[3,1,2]``, ``[3,2,1]``.

        EXAMPLES::

            sage: P = Permutations(6)
            sage: p = P([4,5,1,6,3,2])
            sage: pl = [ [1,2,3], [1,3,2], [3,1,2], [3,2,1] ]
            sage: for q in pl:
            ....:     s = p.simion_schmidt(q)
            ....:     print("{} {}".format(s, s.has_pattern(q)))
            [4, 6, 1, 5, 3, 2] False
            [4, 2, 1, 3, 5, 6] False
            [4, 5, 3, 6, 2, 1] False
            [4, 5, 1, 6, 2, 3] False
        """
        if len(list(self)) <= 2:
            return self
        targetPermutation = [self[0]]
        extreme = self[0]
        nonMinima = []
        if avoid == [1,2,3] or avoid == [1,3,2]:
            for i in range(1, len(list(self))):
                if self[i] < extreme:
                    targetPermutation.append(self[i])
                    extreme = self[i]
                else:
                    targetPermutation.append(None)
                    nonMinima.append(self[i])
            nonMinima.sort()
            if avoid == [1,3,2]:
                nonMinima.reverse()
        if avoid == [3,2,1] or avoid == [3,1,2]:
            for i in range(1, len(list(self))):
                if self[i] > extreme:
                    targetPermutation.append(self[i])
                    extreme = self[i]
                else:
                    targetPermutation.append(None)
                    nonMinima.append(self[i])
            nonMinima.sort()
            if avoid == [3,2,1]:
                nonMinima.reverse()

        for i in range(1, len(list(self))):
            if targetPermutation[i] is None:
                targetPermutation[i] = nonMinima.pop()
        return Permutations()(targetPermutation)

    @combinatorial_map(order=2,name='reverse')
    def reverse(self):
        """
        Return the permutation obtained by reversing the list.

        EXAMPLES::

            sage: Permutation([3,4,1,2]).reverse()
            [2, 1, 4, 3]
            sage: Permutation([1,2,3,4,5]).reverse()
            [5, 4, 3, 2, 1]
        """
        return self.__class__(self.parent(), [i for i in reversed(self)] )

    @combinatorial_map(order=2,name='complement')
    def complement(self):
        r"""
        Return the complement of the permutation ``self``.

        The complement of a permutation `w \in S_n` is defined as the
        permutation in `S_n` sending each `i` to `n + 1 - w(i)`.

        EXAMPLES::

            sage: Permutation([1,2,3]).complement()
            [3, 2, 1]
            sage: Permutation([1, 3, 2]).complement()
            [3, 1, 2]
        """
        n = len(self)
        return self.__class__(self.parent(), [n - x + 1 for x in self] )

    @combinatorial_map(name='permutation poset')
    def permutation_poset(self):
        r"""
        Return the permutation poset of ``self``.

        The permutation poset of a permutation `p` is the poset with
        vertices `(i, p(i))` for `i = 1, 2, \ldots, n` (where `n` is the
        size of `p`) and order inherited from `\ZZ \times \ZZ`.

        EXAMPLES::

            sage: Permutation([3,1,5,4,2]).permutation_poset().cover_relations()
            [[(2, 1), (5, 2)],
             [(2, 1), (3, 5)],
             [(2, 1), (4, 4)],
             [(1, 3), (3, 5)],
             [(1, 3), (4, 4)]]
            sage: Permutation([]).permutation_poset().cover_relations()
            []
            sage: Permutation([1,3,2]).permutation_poset().cover_relations()
            [[(1, 1), (2, 3)], [(1, 1), (3, 2)]]
            sage: Permutation([1,2]).permutation_poset().cover_relations()
            [[(1, 1), (2, 2)]]
            sage: P = Permutation([1,5,2,4,3])
            sage: P.permutation_poset().greene_shape() == P.RS_partition()   # This should hold for any P.
            True
        """
        from sage.combinat.posets.posets import Poset
        n = len(self)
        posetdict = {}
        for i in range(n):
            u = self[i]
            posetdict[(i + 1, u)] = [(j + 1, self[j]) for j in range(i + 1, n) if u < self[j]]
        return Poset(posetdict)

    def dict(self):
        """
        Return a dictionary corresponding to the permutation.

        EXAMPLES::

            sage: p = Permutation([2,1,3])
            sage: d = p.dict()
            sage: d[1]
            2
            sage: d[2]
            1
            sage: d[3]
            3
        """
        return dict(enumerate(self, start=1))

    def action(self, a):
        r"""
        Return the action of the permutation ``self`` on a list ``a``.

        The action of a permutation `p \in S_n` on an `n`-element list
        `(a_1, a_2, \ldots, a_n)` is defined to be
        `(a_{p(1)}, a_{p(2)}, \ldots, a_{p(n)})`.

        EXAMPLES::

            sage: p = Permutation([2,1,3])
            sage: a = list(range(3))
            sage: p.action(a)
            [1, 0, 2]
            sage: b = [1,2,3,4]
            sage: p.action(b)
            Traceback (most recent call last):
            ...
            ValueError: len(a) must equal len(self)

            sage: q = Permutation([2,3,1])
            sage: a = list(range(3))
            sage: q.action(a)
            [1, 2, 0]
        """
        if len(a) != len(self):
            raise ValueError("len(a) must equal len(self)")
        return [a[i-1] for i in self]

    ######################
    # Robinson-Schensted #
    ######################

    def robinson_schensted(self):
        """
        Return the pair of standard tableaux obtained by running the
        Robinson-Schensted algorithm on ``self``.

        This can also be done by running
        :func:`~sage.combinat.rsk.RSK` on ``self`` (with the optional argument
        ``check_standard=True`` to return standard Young tableaux).

        EXAMPLES::

            sage: Permutation([6,2,3,1,7,5,4]).robinson_schensted()
            [[[1, 3, 4], [2, 5], [6, 7]], [[1, 3, 5], [2, 6], [4, 7]]]
        """
        return RSK(self, check_standard=True)

    def _rsk_iter(self):
        r"""
        An iterator for RSK.

        Yields pairs ``[i, p(i)]`` for a permutation ``p``.

        EXAMPLES::

            sage: for x in Permutation([6,2,3,1,7,5,4])._rsk_iter(): x
            (1, 6)
            (2, 2)
            (3, 3)
            (4, 1)
            (5, 7)
            (6, 5)
            (7, 4)
        """
        return zip(range(1, len(self) + 1), self)

    @combinatorial_map(name='Robinson-Schensted insertion tableau')
    def left_tableau(self):
        """
        Return the left standard tableau after performing the RSK
        algorithm on ``self``.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).left_tableau()
            [[1, 2], [3], [4]]
        """
        return RSK(self, check_standard=True)[0]

    @combinatorial_map(name='Robinson-Schensted recording tableau')
    def right_tableau(self):
        """
        Return the right standard tableau after performing the RSK
        algorithm on ``self``.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).right_tableau()
            [[1, 2], [3], [4]]
        """
        return RSK(self, check_standard=True)[1]

    def increasing_tree(self, compare=min):
        """
        Return the increasing tree associated to ``self``.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).increasing_tree()
            1[., 2[3[4[., .], .], .]]
            sage: Permutation([4,1,3,2]).increasing_tree()
            1[4[., .], 2[3[., .], .]]

        By passing the option ``compare=max`` one can have the decreasing
        tree instead::

            sage: Permutation([2,3,4,1]).increasing_tree(max)
            4[3[2[., .], .], 1[., .]]
            sage: Permutation([2,3,1,4]).increasing_tree(max)
            4[3[2[., .], 1[., .]], .]
        """
        from sage.combinat.binary_tree import LabelledBinaryTree as LBT

        def rec(perm):
            if len(perm) == 0:
                return LBT(None)
            mn = compare(perm)
            k = perm.index(mn)
            return LBT([rec(perm[:k]), rec(perm[k + 1:])], label = mn)
        return rec(self)

    @combinatorial_map(name="Increasing tree")
    def increasing_tree_shape(self, compare=min):
        r"""
        Return the shape of the increasing tree associated with the
        permutation.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).increasing_tree_shape()
            [., [[[., .], .], .]]
            sage: Permutation([4,1,3,2]).increasing_tree_shape()
            [[., .], [[., .], .]]

        By passing the option ``compare=max`` one can have the decreasing
        tree instead::

            sage: Permutation([2,3,4,1]).increasing_tree_shape(max)
            [[[., .], .], [., .]]
            sage: Permutation([2,3,1,4]).increasing_tree_shape(max)
            [[[., .], [., .]], .]
        """
        return self.increasing_tree(compare).shape()

    def binary_search_tree(self, left_to_right=True):
        """
        Return the binary search tree associated to ``self``.

        If `w` is a word, then the binary search tree associated to `w`
        is defined as the result of starting with an empty binary tree,
        and then inserting the letters of `w` one by one into this tree.
        Here, the insertion is being done according to the method
        :meth:`~sage.combinat.binary_tree.LabelledBinaryTree.binary_search_insert`,
        and the word `w` is being traversed from left to right.

        A permutation is regarded as a word (using one-line notation),
        and thus a binary search tree associated to a permutation is
        defined.

        If the optional keyword variable ``left_to_right`` is set to
        ``False``, the word `w` is being traversed from right to left
        instead.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).binary_search_tree()
            1[., 4[3[2[., .], .], .]]
            sage: Permutation([4,1,3,2]).binary_search_tree()
            4[1[., 3[2[., .], .]], .]

        By passing the option ``left_to_right=False`` one can have
        the insertion going from right to left::

            sage: Permutation([1,4,3,2]).binary_search_tree(False)
            2[1[., .], 3[., 4[., .]]]
            sage: Permutation([4,1,3,2]).binary_search_tree(False)
            2[1[., .], 3[., 4[., .]]]

        TESTS::

            sage: Permutation([]).binary_search_tree()
            .
        """
        from sage.combinat.binary_tree import LabelledBinaryTree as LBT
        res = LBT(None)
        if left_to_right:
            gen = self
        else:
            gen = self[::-1]
        for i in gen:
            res = res.binary_search_insert(i)
        return res

    @combinatorial_map(name = "Binary search tree (left to right)")
    def binary_search_tree_shape(self, left_to_right=True):
        r"""
        Return the shape of the binary search tree of the permutation
        (a non labelled binary tree).

        EXAMPLES::

            sage: Permutation([1,4,3,2]).binary_search_tree_shape()
            [., [[[., .], .], .]]
            sage: Permutation([4,1,3,2]).binary_search_tree_shape()
            [[., [[., .], .]], .]

        By passing the option ``left_to_right=False`` one can have
        the insertion going from right to left::

            sage: Permutation([1,4,3,2]).binary_search_tree_shape(False)
            [[., .], [., [., .]]]
            sage: Permutation([4,1,3,2]).binary_search_tree_shape(False)
            [[., .], [., [., .]]]
        """
        from sage.combinat.binary_tree import binary_search_tree_shape
        return binary_search_tree_shape(list(self), left_to_right)

    def sylvester_class(self, left_to_right=False):
        r"""
        Iterate over the equivalence class of the permutation ``self``
        under sylvester congruence.

        Sylvester congruence is an equivalence relation on the set `S_n`
        of all permutations of `n`. It is defined as the smallest
        equivalence relation such that every permutation of the form
        `uacvbw` with `u`, `v` and `w` being words and `a`, `b` and `c`
        being letters satisfying `a \leq b < c` is equivalent to the
        permutation `ucavbw`. (Here, permutations are regarded as words
        by way of one-line notation.) This definition comes from [HNT2005]_,
        Definition 8, where it is more generally applied to arbitrary
        words.

        The equivalence class of a permutation `p \in S_n` under sylvester
        congruence is called the *sylvester class* of `p`. It is an
        interval in the right permutohedron order (see
        :meth:`permutohedron_lequal`) on `S_n`.

        This is related to the
        :meth:`~sage.combinat.binary_tree.LabelledBinaryTree.sylvester_class`
        method in that the equivalence class of a permutation `\pi` under
        sylvester congruence is the sylvester class of the right-to-left
        binary search tree of `\pi`. However, the present method
        yields permutations, while the method on labelled binary trees
        yields plain lists.

        If the variable ``left_to_right`` is set to ``True``, the method
        instead iterates over the equivalence class of ``self`` with
        respect to the *left* sylvester congruence. The left sylvester
        congruence is easiest to define by saying that two permutations
        are equivalent under it if and only if their reverses
        (:meth:`reverse`) are equivalent under (standard) sylvester
        congruence.

        EXAMPLES:

        The sylvester class of a permutation in `S_5`::

            sage: p = Permutation([3, 5, 1, 2, 4])
            sage: sorted(p.sylvester_class())
            [[1, 3, 2, 5, 4],
             [1, 3, 5, 2, 4],
             [1, 5, 3, 2, 4],
             [3, 1, 2, 5, 4],
             [3, 1, 5, 2, 4],
             [3, 5, 1, 2, 4],
             [5, 1, 3, 2, 4],
             [5, 3, 1, 2, 4]]

        The sylvester class of a permutation `p` contains `p`::

            sage: all( p in p.sylvester_class() for p in Permutations(4) )
            True

        Small cases::

            sage: list(Permutation([]).sylvester_class())
            [[]]

            sage: list(Permutation([1]).sylvester_class())
            [[1]]

        The sylvester classes in `S_3`::

            sage: [sorted(p.sylvester_class()) for p in Permutations(3)]
            [[[1, 2, 3]],
             [[1, 3, 2], [3, 1, 2]],
             [[2, 1, 3]],
             [[2, 3, 1]],
             [[1, 3, 2], [3, 1, 2]],
             [[3, 2, 1]]]

        The left sylvester classes in `S_3`::

            sage: [sorted(p.sylvester_class(left_to_right=True)) for p in Permutations(3)]
            [[[1, 2, 3]],
             [[1, 3, 2]],
             [[2, 1, 3], [2, 3, 1]],
             [[2, 1, 3], [2, 3, 1]],
             [[3, 1, 2]],
             [[3, 2, 1]]]

        A left sylvester class in `S_5`::

            sage: p = Permutation([4, 2, 1, 5, 3])
            sage: sorted(p.sylvester_class(left_to_right=True))
            [[4, 2, 1, 3, 5],
             [4, 2, 1, 5, 3],
             [4, 2, 3, 1, 5],
             [4, 2, 3, 5, 1],
             [4, 2, 5, 1, 3],
             [4, 2, 5, 3, 1],
             [4, 5, 2, 1, 3],
             [4, 5, 2, 3, 1]]
        """
        parself = self.parent()
        t = self.binary_search_tree(left_to_right=left_to_right)
        for u in t.sylvester_class(left_to_right=left_to_right):
            yield parself(u)

    @combinatorial_map(name='Robinson-Schensted tableau shape')
    def RS_partition(self):
        """
        Return the shape of the tableaux obtained by applying the RSK
        algorithm to ``self``.

        EXAMPLES::

            sage: Permutation([1,4,3,2]).RS_partition()
            [2, 1, 1]
        """
        return RSK(self)[1].shape()

    def remove_extra_fixed_points(self):
        """
        Return the permutation obtained by removing any fixed points at
        the end of ``self``.

        EXAMPLES::

            sage: Permutation([2,1,3]).remove_extra_fixed_points()
            [2, 1]
            sage: Permutation([1,2,3,4]).remove_extra_fixed_points()
            [1]

        .. SEEALSO::

            :meth:`retract_plain`
        """
        #Strip off all extra fixed points at the end of
        #the permutation.
        i = len(self)-1
        while i >= 1:
            if i != self[i] - 1:
                break
            i -= 1
        return Permutations()(self[:i+1])

    def retract_plain(self, m):
        r"""
        Return the plain retract of the permutation ``self`` in `S_n`
        to `S_m`, where `m \leq n`. If this retract is undefined, then
        ``None`` is returned.

        If `p \in S_n` is a permutation, and `m` is a nonnegative integer
        less or equal to `n`, then the plain retract of `p` to `S_m` is
        defined only if every `i > m` satisfies `p(i) = i`. In this case,
        it is defined as the permutation written
        `(p(1), p(2), \ldots, p(m))` in one-line notation.

        EXAMPLES::

            sage: Permutation([4,1,2,3,5]).retract_plain(4)
            [4, 1, 2, 3]
            sage: Permutation([4,1,2,3,5]).retract_plain(3)

            sage: Permutation([1,3,2,4,5,6]).retract_plain(3)
            [1, 3, 2]
            sage: Permutation([1,3,2,4,5,6]).retract_plain(2)

            sage: Permutation([1,2,3,4,5]).retract_plain(1)
            [1]
            sage: Permutation([1,2,3,4,5]).retract_plain(0)
            []

            sage: all( p.retract_plain(3) == p for p in Permutations(3) )
            True

        .. SEEALSO::

            :meth:`retract_direct_product`, :meth:`retract_okounkov_vershik`,
            :meth:`remove_extra_fixed_points`
        """
        n = len(self)
        p = list(self)
        for i in range(m, n):
            if p[i] != i + 1:
                return None
        return Permutations(m)(p[:m])

    def retract_direct_product(self, m):
        r"""
        Return the direct-product retract of the permutation
        ``self`` `\in S_n` to `S_m`, where `m \leq n`. If this retract
        is undefined, then ``None`` is returned.

        If `p \in S_n` is a permutation, and `m` is a nonnegative integer
        less or equal to `n`, then the direct-product retract of `p` to
        `S_m` is defined only if `p([m]) = [m]`, where `[m]` denotes the
        interval `\{1, 2, \ldots, m\}`. In this case, it is defined as the
        permutation written `(p(1), p(2), \ldots, p(m))` in one-line
        notation.

        EXAMPLES::

            sage: Permutation([4,1,2,3,5]).retract_direct_product(4)
            [4, 1, 2, 3]
            sage: Permutation([4,1,2,3,5]).retract_direct_product(3)

            sage: Permutation([1,4,2,3,6,5]).retract_direct_product(5)
            sage: Permutation([1,4,2,3,6,5]).retract_direct_product(4)
            [1, 4, 2, 3]
            sage: Permutation([1,4,2,3,6,5]).retract_direct_product(3)
            sage: Permutation([1,4,2,3,6,5]).retract_direct_product(2)
            sage: Permutation([1,4,2,3,6,5]).retract_direct_product(1)
            [1]
            sage: Permutation([1,4,2,3,6,5]).retract_direct_product(0)
            []

            sage: all( p.retract_direct_product(3) == p for p in Permutations(3) )
            True

        .. SEEALSO::

            :meth:`retract_plain`, :meth:`retract_okounkov_vershik`
        """
        n = len(self)
        p = list(self)
        for i in range(m, n):
            if p[i] <= m:
                return None
        return Permutations(m)(p[:m])

    def retract_okounkov_vershik(self, m):
        r"""
        Return the Okounkov-Vershik retract of the permutation
        ``self`` `\in S_n` to `S_m`, where `m \leq n`.

        If `p \in S_n` is a permutation, and `m` is a nonnegative integer
        less or equal to `n`, then the Okounkov-Vershik retract of `p` to
        `S_m` is defined as the permutation in `S_m` which sends every
        `i \in \{1, 2, \ldots, m\}` to `p^{k_i}(i)`, where `k_i` is the
        smallest positive integer `k` satisfying `p^k(i) \leq m`.

        In other words, the Okounkov-Vershik retract of `p` is the
        permutation whose disjoint cycle decomposition is obtained by
        removing all letters strictly greater than `m` from the
        decomposition of `p` into disjoint cycles (and removing all
        cycles which are emptied in the process).

        When `m = n-1`, the Okounkov-Vershik retract (as a map
        `S_n \to S_{n-1}`) is the map `\widetilde{p}_n` introduced in
        Section 7 of [VO2005]_, and appears as (3.20) in
        [CST2010]_. In the general case, the Okounkov-Vershik retract
        of a permutation in `S_n` to `S_m` can be obtained by first
        taking its Okounkov-Vershik retract to `S_{n-1}`, then that
        of the resulting permutation to `S_{n-2}`, etc. until arriving
        in `S_m`.

        EXAMPLES::

            sage: Permutation([4,1,2,3,5]).retract_okounkov_vershik(4)
            [4, 1, 2, 3]
            sage: Permutation([4,1,2,3,5]).retract_okounkov_vershik(3)
            [3, 1, 2]
            sage: Permutation([4,1,2,3,5]).retract_okounkov_vershik(2)
            [2, 1]
            sage: Permutation([4,1,2,3,5]).retract_okounkov_vershik(1)
            [1]
            sage: Permutation([4,1,2,3,5]).retract_okounkov_vershik(0)
            []

            sage: Permutation([1,4,2,3,6,5]).retract_okounkov_vershik(5)
            [1, 4, 2, 3, 5]
            sage: Permutation([1,4,2,3,6,5]).retract_okounkov_vershik(4)
            [1, 4, 2, 3]
            sage: Permutation([1,4,2,3,6,5]).retract_okounkov_vershik(3)
            [1, 3, 2]
            sage: Permutation([1,4,2,3,6,5]).retract_okounkov_vershik(2)
            [1, 2]
            sage: Permutation([1,4,2,3,6,5]).retract_okounkov_vershik(1)
            [1]
            sage: Permutation([1,4,2,3,6,5]).retract_okounkov_vershik(0)
            []

            sage: Permutation([6,5,4,3,2,1]).retract_okounkov_vershik(5)
            [1, 5, 4, 3, 2]
            sage: Permutation([6,5,4,3,2,1]).retract_okounkov_vershik(4)
            [1, 2, 4, 3]

            sage: Permutation([1,5,2,6,3,7,4,8]).retract_okounkov_vershik(4)
            [1, 3, 2, 4]

            sage: all( p.retract_direct_product(3) == p for p in Permutations(3) )
            True

        .. SEEALSO::

            :meth:`retract_plain`, :meth:`retract_direct_product`
        """
        res = []
        for i in range(1, m + 1):
            j = self(i)
            while j > m:
                j = self(j)
            res.append(j)
        return Permutations(m)(res)

    def hyperoctahedral_double_coset_type(self):
        r"""
        Return the coset-type of ``self`` as a partition.

        ``self`` must be a permutation of even size `2n`.  The coset-type
        determines the double class of the permutation, that is its image in
        `H_n \backslash S_{2n} / H_n`, where `H_n` is the `n`-th
        hyperoctahedral group.

        The coset-type is determined as follows. Consider the perfect matching
        `\{\{1,2\},\{3,4\},\dots,\{2n-1,2n\}\}` and its image by ``self``, and
        draw them simultaneously as edges of a graph whose vertices are labeled
        by `1,2,\dots,2n`. The coset-type is the ordered sequence of the
        semi-lengths of the cycles of this graph (see Chapter VII of [Mac1995]_ for
        more details, particularly Section VII.2).

        EXAMPLES::

            sage: Permutation([3, 4, 6, 1, 5, 7, 2, 8]).hyperoctahedral_double_coset_type()
            [3, 1]
            sage: all(p.hyperoctahedral_double_coset_type() ==
            ....:     p.inverse().hyperoctahedral_double_coset_type()
            ....:     for p in Permutations(4))
            True
            sage: Permutation([]).hyperoctahedral_double_coset_type()
            []
            sage: Permutation([3,1,2]).hyperoctahedral_double_coset_type()
            Traceback (most recent call last):
            ...
            ValueError: [3, 1, 2] is a permutation of odd size and has no coset-type
        """
        from sage.combinat.perfect_matching import PerfectMatchings
        n = len(self)
        if n % 2 == 1:
            raise ValueError("%s is a permutation of odd size and has no coset-type"%self)
        S = PerfectMatchings(n)([(2*i+1,2*i+2) for i in range(n//2)])
        return S.loop_type(S.apply_permutation(self))

    #####################
    # Binary operations #
    #####################

    def shifted_concatenation(self, other, side = "right"):
        r"""
        Return the right (or left) shifted concatenation of ``self``
        with a permutation ``other``. These operations are also known
        as the Loday-Ronco over and under operations.

        INPUT:

        - ``other`` -- a permutation, a list, a tuple, or any iterable
          representing a permutation.

        - ``side`` -- (default: ``"right"``) the string "left" or "right".

        OUTPUT:

        If ``side`` is ``"right"``, the method returns the permutation
        obtained by concatenating ``self`` with the letters of ``other``
        incremented by the size of ``self``. This is what is called
        ``side / other`` in [LR0102066]_, and denoted as the "over"
        operation.
        Otherwise, i. e., when ``side`` is ``"left"``, the method
        returns the permutation obtained by concatenating the letters
        of ``other`` incremented by the size of ``self`` with ``self``.
        This is what is called ``side \ other`` in [LR0102066]_
        (which seems to use the `(\sigma \pi)(i) = \pi(\sigma(i))`
        convention for the product of permutations).

        EXAMPLES::

            sage: Permutation([]).shifted_concatenation(Permutation([]), "right")
            []
            sage: Permutation([]).shifted_concatenation(Permutation([]), "left")
            []
            sage: Permutation([2, 4, 1, 3]).shifted_concatenation(Permutation([3, 1, 2]), "right")
            [2, 4, 1, 3, 7, 5, 6]
            sage: Permutation([2, 4, 1, 3]).shifted_concatenation(Permutation([3, 1, 2]), "left")
            [7, 5, 6, 2, 4, 1, 3]
        """
        if side == "right":
            return Permutations()(list(self) + [a + len(self) for a in other])
        elif side == "left":
            return Permutations()([a + len(self) for a in other] + list(self))
        else:
            raise ValueError("%s must be \"left\" or \"right\"" % side)

    def shifted_shuffle(self, other):
        r"""
        Return the shifted shuffle of two permutations ``self`` and ``other``.

        INPUT:

        - ``other`` -- a permutation, a list, a tuple, or any iterable
          representing a permutation.

        OUTPUT:

        The list of the permutations appearing in the shifted
        shuffle of the permutations ``self`` and ``other``.

        EXAMPLES::

            sage: Permutation([]).shifted_shuffle(Permutation([]))
            [[]]
            sage: Permutation([1, 2, 3]).shifted_shuffle(Permutation([1]))
            [[4, 1, 2, 3], [1, 2, 3, 4], [1, 2, 4, 3], [1, 4, 2, 3]]
            sage: Permutation([1, 2]).shifted_shuffle(Permutation([2, 1]))
            [[4, 1, 3, 2], [4, 3, 1, 2], [1, 4, 3, 2], [1, 4, 2, 3], [1, 2, 4, 3], [4, 1, 2, 3]]
            sage: Permutation([1]).shifted_shuffle([1])
            [[2, 1], [1, 2]]
            sage: len(Permutation([3, 1, 5, 4, 2]).shifted_shuffle(Permutation([2, 1, 4, 3])))
            126

        The shifted shuffle product is associative. We can test this on an
        admittedly toy example::

            sage: all( all( all( sorted(flatten([abs.shifted_shuffle(c)
            ....:                                for abs in a.shifted_shuffle(b)]))
            ....:                == sorted(flatten([a.shifted_shuffle(bcs)
            ....:                                   for bcs in b.shifted_shuffle(c)]))
            ....:                for c in Permutations(2) )
            ....:           for b in Permutations(2) )
            ....:      for a in Permutations(2) )
            True

        The ``shifted_shuffle`` method on permutations gives the same
        permutations as the ``shifted_shuffle`` method on words (but is
        faster)::

            sage: all( all( sorted(p1.shifted_shuffle(p2))
            ....:           == sorted([Permutation(p) for p in
            ....:                      Word(p1).shifted_shuffle(Word(p2))])
            ....:           for p2 in Permutations(3) )
            ....:      for p1 in Permutations(2) )
            True
        """
        return self.shifted_concatenation(other, "right").\
        right_permutohedron_interval(self.shifted_concatenation(other, "left"))

################################################################
# Parent classes
################################################################

# Base class for permutations
class Permutations(UniqueRepresentation, Parent):
    r"""
    Permutations.

    ``Permutations(n)`` returns the class of permutations of ``n``, if ``n``
    is an integer, list, set, or string.

    ``Permutations(n, k)`` returns the class of length-``k`` partial
    permutations of ``n`` (where ``n`` is any of the above things); ``k``
    must be a nonnegative integer. A length-`k` partial permutation of `n`
    is defined as a `k`-tuple of pairwise distinct elements of
    `\{ 1, 2, \ldots, n \}`.

    Valid keyword arguments are: 'descents', 'bruhat_smaller',
    'bruhat_greater', 'recoils_finer', 'recoils_fatter', 'recoils',
    and 'avoiding'. With the exception of 'avoiding', you cannot
    specify ``n`` or ``k`` along with a keyword.

    ``Permutations(descents=(list,n))`` returns the class of permutations of
    `n` with descents in the positions specified by ``list``. This uses the
    slightly nonstandard convention that the images of `1,2,...,n` under the
    permutation are regarded as positions `0,1,...,n-1`, so for example the
    presence of `1` in ``list`` signifies that the permutations `\pi` should
    satisfy `\pi(2) > \pi(3)`.
    Note that ``list`` is supposed to be a list of positions of the descents,
    not the descents composition.  It does *not* return the class of
    permutations with descents composition ``list``.

    ``Permutations(bruhat_smaller=p)`` and ``Permutations(bruhat_greater=p)``
    return the class of permutations smaller-or-equal or greater-or-equal,
    respectively, than the given permutation ``p`` in the Bruhat order.
    (The Bruhat order is defined in
    :meth:`~sage.combinat.permutation.Permutation.bruhat_lequal`.
    It is also referred to as the *strong* Bruhat order.)

    ``Permutations(recoils=p)`` returns the class of permutations whose
    recoils composition is ``p``. Unlike the ``descents=(list, n)`` syntax,
    this actually takes a *composition* as input.

    ``Permutations(recoils_fatter=p)`` and ``Permutations(recoils_finer=p)``
    return the class of permutations whose recoils composition is fatter or
    finer, respectively, than the given composition ``p``.

    ``Permutations(n, avoiding=P)`` returns the class of permutations of ``n``
    avoiding ``P``. Here ``P`` may be a single permutation or a list of
    permutations; the returned class will avoid all patterns in ``P``.

    EXAMPLES::

        sage: p = Permutations(3); p
        Standard permutations of 3
        sage: p.list()
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

    ::

        sage: p = Permutations(3, 2); p
        Permutations of {1,...,3} of length 2
        sage: p.list()
        [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]

    ::

        sage: p = Permutations(['c', 'a', 't']); p
        Permutations of the set ['c', 'a', 't']
        sage: p.list()
        [['c', 'a', 't'],
         ['c', 't', 'a'],
         ['a', 'c', 't'],
         ['a', 't', 'c'],
         ['t', 'c', 'a'],
         ['t', 'a', 'c']]

    ::

        sage: p = Permutations(['c', 'a', 't'], 2); p
        Permutations of the set ['c', 'a', 't'] of length 2
        sage: p.list()
        [['c', 'a'], ['c', 't'], ['a', 'c'], ['a', 't'], ['t', 'c'], ['t', 'a']]

    ::

        sage: p = Permutations([1,1,2]); p
        Permutations of the multi-set [1, 1, 2]
        sage: p.list()
        [[1, 1, 2], [1, 2, 1], [2, 1, 1]]

    ::

        sage: p = Permutations([1,1,2], 2); p
        Permutations of the multi-set [1, 1, 2] of length 2
        sage: p.list()
        [[1, 1], [1, 2], [2, 1]]

    ::

        sage: p = Permutations(descents=([1], 4)); p
        Standard permutations of 4 with descents [1]
        sage: p.list()
        [[2, 4, 1, 3], [3, 4, 1, 2], [1, 4, 2, 3], [1, 3, 2, 4], [2, 3, 1, 4]]

    ::

        sage: p = Permutations(bruhat_smaller=[1,3,2,4]); p
        Standard permutations that are less than or equal to [1, 3, 2, 4] in the Bruhat order
        sage: p.list()
        [[1, 2, 3, 4], [1, 3, 2, 4]]

    ::

        sage: p = Permutations(bruhat_greater=[4,2,3,1]); p
        Standard permutations that are greater than or equal to [4, 2, 3, 1] in the Bruhat order
        sage: p.list()
        [[4, 2, 3, 1], [4, 3, 2, 1]]

    ::

        sage: p = Permutations(recoils_finer=[2,1]); p
        Standard permutations whose recoils composition is finer than [2, 1]
        sage: p.list()
        [[3, 1, 2], [1, 2, 3], [1, 3, 2]]

    ::

        sage: p = Permutations(recoils_fatter=[2,1]); p
        Standard permutations whose recoils composition is fatter than [2, 1]
        sage: p.list()
        [[3, 1, 2], [3, 2, 1], [1, 3, 2]]

    ::

        sage: p = Permutations(recoils=[2,1]); p
        Standard permutations whose recoils composition is [2, 1]
        sage: p.list()
        [[3, 1, 2], [1, 3, 2]]

    ::

        sage: p = Permutations(4, avoiding=[1,3,2]); p
        Standard permutations of 4 avoiding [[1, 3, 2]]
        sage: p.list()
        [[4, 1, 2, 3],
         [4, 2, 1, 3],
         [4, 2, 3, 1],
         [4, 3, 1, 2],
         [4, 3, 2, 1],
         [3, 4, 1, 2],
         [3, 4, 2, 1],
         [2, 3, 4, 1],
         [3, 2, 4, 1],
         [1, 2, 3, 4],
         [2, 1, 3, 4],
         [2, 3, 1, 4],
         [3, 1, 2, 4],
         [3, 2, 1, 4]]

    ::

        sage: p = Permutations(5, avoiding=[[3,4,1,2], [4,2,3,1]]); p
        Standard permutations of 5 avoiding [[3, 4, 1, 2], [4, 2, 3, 1]]
        sage: p.cardinality()
        88
        sage: p.random_element().parent() is p
        True
    """
    @staticmethod
    def __classcall_private__(cls, n=None, k=None, **kwargs):
        """
        Return the correct parent based upon input.

        EXAMPLES::

            sage: Permutations()
            Standard permutations
            sage: Permutations(5, 3)
            Permutations of {1,...,5} of length 3
            sage: Permutations([1,2,3,4,6])
            Permutations of the set [1, 2, 3, 4, 6]
            sage: Permutations([1,2,3,4,5])
            Standard permutations of 5
        """
        valid_args = ['descents', 'bruhat_smaller', 'bruhat_greater',
                      'recoils_finer', 'recoils_fatter', 'recoils', 'avoiding']

        number_of_arguments = 0
        if n is not None:
            number_of_arguments += 1
        elif k is not None:
            number_of_arguments += 1

        #Make sure that exactly one keyword was passed
        for key in kwargs:
            if key not in valid_args:
                raise ValueError("unknown keyword argument: %s"%key)
            if key not in [ 'avoiding' ]:
                number_of_arguments += 1

        if number_of_arguments == 0:
            if 'avoiding' in kwargs:
                a = kwargs['avoiding']
                if len(a) == 0:
                    return StandardPermutations_all()
                if a in StandardPermutations_all():
                    a = (a,)
                return StandardPermutations_all_avoiding(a)
            return StandardPermutations_all()

        if number_of_arguments != 1:
            raise ValueError("you must specify exactly one argument")

        if n is not None:
            if isinstance(n, (int, Integer)):
                if k is None:
                    if 'avoiding' in kwargs:
                        a = kwargs['avoiding']
                        if not a:
                            return StandardPermutations_n(n)
                        if len(a) == 1 and a[0] != 1:
                            a = a[0]
                        if a in StandardPermutations_all():
                            if a == [1,2]:
                                return StandardPermutations_avoiding_12(n)
                            elif a == [2,1]:
                                return StandardPermutations_avoiding_21(n)
                            elif a == [1,2,3]:
                                return StandardPermutations_avoiding_123(n)
                            elif a == [1,3,2]:
                                return StandardPermutations_avoiding_132(n)
                            elif a == [2,1,3]:
                                return StandardPermutations_avoiding_213(n)
                            elif a == [2,3,1]:
                                return StandardPermutations_avoiding_231(n)
                            elif a == [3,1,2]:
                                return StandardPermutations_avoiding_312(n)
                            elif a == [3,2,1]:
                                return StandardPermutations_avoiding_321(n)
                            else:
                                return StandardPermutations_avoiding_generic(n, (a,))
                        elif isinstance(a, (list, tuple)):
                            a = tuple(map(Permutation, a))
                            return StandardPermutations_avoiding_generic(n, a)
                        else:
                            raise ValueError("do not know how to avoid %s"%a)
                    else:
                        return StandardPermutations_n(n)
                else:
                    return Permutations_nk(n,k)
            else:
                # In this case, we have that n is a list
                # Because of UniqueRepresentation, we require the elements
                #   to be hashable
                if len(set(n)) == len(n):
                    if list(n) == list(range(1, len(n)+1)):
                        if k is None:
                            return StandardPermutations_n(len(n))
                        else:
                            return Permutations_nk(len(n), k)
                    else:
                        if k is None:
                            return Permutations_set(n)
                        else:
                            return Permutations_setk(n,k)
                else:
                    if k is None:
                        return Permutations_mset(n)
                    else:
                        return Permutations_msetk(n, k)
        elif 'descents' in kwargs:
            #Descent positions specified
            if isinstance(kwargs['descents'], tuple):
                #Descent positions and size specified
                args = kwargs['descents']
                return StandardPermutations_descents(tuple(args[0]), args[1])
            else:
                #Size not specified
                return StandardPermutations_descents(kwargs['descents'])
        elif 'bruhat_smaller' in kwargs:
            return StandardPermutations_bruhat_smaller(Permutation(kwargs['bruhat_smaller']))
        elif 'bruhat_greater' in kwargs:
            return StandardPermutations_bruhat_greater(Permutation(kwargs['bruhat_greater']))
        elif 'recoils_finer' in kwargs:
            return StandardPermutations_recoilsfiner(Composition(kwargs['recoils_finer']))
        elif 'recoils_fatter' in kwargs:
            return StandardPermutations_recoilsfatter(Composition(kwargs['recoils_fatter']))
        elif 'recoils' in kwargs:
            return StandardPermutations_recoils(Composition(kwargs['recoils']))

    Element = Permutation

    # add options to class
    class options(GlobalOptions):
        r"""
        Set the global options for elements of the permutation class. The
        defaults are for permutations to be displayed in list notation and
        the multiplication done from left to right (like in GAP) -- that
        is, `(\pi \psi)(i) = \psi(\pi(i))` for all `i`.

        .. NOTE::

            These options have no effect on permutation group elements.

        @OPTIONS@

        EXAMPLES::

            sage: p213 = Permutation([2,1,3])
            sage: p312 = Permutation([3,1,2])
            sage: Permutations.options(mult='l2r', display='list')
            sage: Permutations.options.display
            list
            sage: p213
            [2, 1, 3]
            sage: Permutations.options.display='cycle'
            sage: p213
            (1,2)
            sage: Permutations.options.display='singleton'
            sage: p213
            (1,2)(3)
            sage: Permutations.options.display='list'

        ::

            sage: Permutations.options.mult
            l2r
            sage: p213*p312
            [1, 3, 2]
            sage: Permutations.options.mult='r2l'
            sage: p213*p312
            [3, 2, 1]
            sage: Permutations.options._reset()
        """
        NAME = 'Permutations'
        module = 'sage.combinat.permutation'
        display = dict(default="list",
                     description="Specifies how the permutations should be printed",
                     values=dict(list="the permutations are displayed in list notation"
                                      " (aka 1-line notation)",
                                 cycle="the permutations are displayed in cycle notation"
                                      " (i. e., as products of disjoint cycles)",
                                 singleton="the permutations are displayed in cycle notation"
                                           " with singleton cycles shown as well",
                                 reduced_word="the permutations are displayed as reduced words"),
                     alias=dict(word="reduced_word", reduced_expression="reduced_word"),
                     case_sensitive=False)
        latex = dict(default="list",
                   description="Specifies how the permutations should be latexed",
                   values=dict(list="latex as a list in one-line notation",
                               twoline="latex in two-line notation",
                               cycle="latex in cycle notation",
                               singleton="latex in cycle notation with singleton cycles shown as well",
                               reduced_word="latex as reduced words"),
                   alias=dict(word="reduced_word", reduced_expression="reduced_word", oneline="list"),
                   case_sensitive=False)
        latex_empty_str = dict(default="1",
                             description='The LaTeX representation of a reduced word when said word is empty',
                             checker=lambda char: isinstance(char,str))
        generator_name = dict(default="s",
                            description="the letter used in latexing the reduced word",
                            checker=lambda char: isinstance(char,str))
        mult = dict(default="l2r",
                  description="The multiplication of permutations",
                  values=dict(l2r=r"left to right: `(p_1 \cdot p_2)(x) = p_2(p_1(x))`",
                              r2l=r"right to left: `(p_1 \cdot p_2)(x) = p_1(p_2(x))`"),
                  case_sensitive=False)


class Permutations_nk(Permutations):
    r"""
    Length-`k` partial permutations of `\{1, 2, \ldots, n\}`.
    """
    def __init__(self, n, k):
        """
        TESTS::

            sage: P = Permutations(3,2)
            sage: TestSuite(P).run()
        """
        self.n = ZZ(n)
        self._k = ZZ(k)
        Permutations.__init__(self, category=FiniteEnumeratedSets())

    class Element(ClonableArray):
        """
        A length-`k` partial permutation of `[n]`.
        """
        def check(self):
            """
            Verify that ``self`` is a valid length-`k` partial
            permutation of `[n]`.

            EXAMPLES::

                sage: S = Permutations(4, 2)
                sage: elt = S([3, 1])
                sage: elt.check()
            """
            if self not in self.parent():
                raise ValueError("Invalid permutation")

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: [1,2] in Permutations(3,2)
            True
            sage: [1,1] in Permutations(3,2)
            False
            sage: [3,2,1] in Permutations(3,2)
            False
            sage: [3,1] in Permutations(3,2)
            True
        """
        if len(x) != self._k:
            return False

        r = list(range(1, self.n + 1))
        for i in x:
            if i in r:
                r.remove(i)
            else:
                return False

        return True

    def _repr_(self) -> str:
        """
        TESTS::

            sage: Permutations(3,2)
            Permutations of {1,...,3} of length 2
        """
        return "Permutations of {1,...,%s} of length %s" % (self.n, self._k)

    def __iter__(self) -> Iterator[Permutation]:
        """
        EXAMPLES::

            sage: [p for p in Permutations(3,2)]
            [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]
            sage: [p for p in Permutations(3,0)]
            [[]]
            sage: [p for p in Permutations(3,4)]
            []
        """
        for x in itertools.permutations(range(1,self.n+1), int(self._k)):
            yield self.element_class(self, x)

    def cardinality(self) -> Integer:
        """
        EXAMPLES::

            sage: Permutations(3,0).cardinality()
            1
            sage: Permutations(3,1).cardinality()
            3
            sage: Permutations(3,2).cardinality()
            6
            sage: Permutations(3,3).cardinality()
            6
            sage: Permutations(3,4).cardinality()
            0
        """
        if 0 <= self._k <= self.n:
            return factorial(self.n) // factorial(self.n - self._k)
        return ZZ.zero()

    def random_element(self):
        """
        EXAMPLES::

            sage: s = Permutations(3,2).random_element()
            sage: s in Permutations(3,2)
            True
        """
        return sample(range(1, self.n+1), self._k)


class Permutations_mset(Permutations):
    r"""
    Permutations of a multiset `M`.

    A permutation of a multiset `M` is represented by a list that
    contains exactly the same elements as `M` (with the same
    multiplicities), but possibly in different order. If `M` is
    a proper set there are `|M| !` such permutations.
    Otherwise, if the first element appears `k_1` times, the
    second element appears `k_2` times and so on, the number
    of permutations is `|M|! / (k_1! k_2! \ldots)`, which
    is sometimes called a multinomial coefficient.

    EXAMPLES::

        sage: mset = [1,1,2,2,2]
        sage: from sage.combinat.permutation import Permutations_mset
        sage: P = Permutations_mset(mset); P
        Permutations of the multi-set [1, 1, 2, 2, 2]
        sage: sorted(P)
        [[1, 1, 2, 2, 2],
         [1, 2, 1, 2, 2],
         [1, 2, 2, 1, 2],
         [1, 2, 2, 2, 1],
         [2, 1, 1, 2, 2],
         [2, 1, 2, 1, 2],
         [2, 1, 2, 2, 1],
         [2, 2, 1, 1, 2],
         [2, 2, 1, 2, 1],
         [2, 2, 2, 1, 1]]
        sage: MS = MatrixSpace(GF(2),2,2)
        sage: A = MS([1,0,1,1])
        sage: rows = A.rows()
        sage: rows[0].set_immutable()
        sage: rows[1].set_immutable()
        sage: P = Permutations_mset(rows); P
        Permutations of the multi-set [(1, 0), (1, 1)]
        sage: sorted(P)
        [[(1, 0), (1, 1)], [(1, 1), (1, 0)]]
    """
    @staticmethod
    def __classcall_private__(cls, mset):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(['c','a','c'])
            sage: S2 = Permutations(('c','a','c'))
            sage: S1 is S2
            True
        """
        return super(Permutations_mset, cls).__classcall__(cls, tuple(mset))

    def __init__(self, mset):
        """
        TESTS::

            sage: S = Permutations(['c','a','c'])
            sage: TestSuite(S).run()
        """
        self.mset = mset
        Permutations.__init__(self, category=FiniteEnumeratedSets())

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: p = Permutations([1,2,2])
            sage: [1,2,2] in p
            True
            sage: [] in p
            False
            sage: [2,2] in p
            False
            sage: [1,1] in p
            False
            sage: [2,1] in p
            False
            sage: [2,1,2] in p
            True
        """
        s = list(self.mset)
        if len(x) != len(s):
            return False
        for i in x:
            if i in s:
                s.remove(i)
            else:
                return False
        return True

    class Element(ClonableArray):
        """
        A permutation of an arbitrary multiset.
        """
        def check(self):
            """
            Verify that ``self`` is a valid permutation of the underlying
            multiset.

            EXAMPLES::

                sage: S = Permutations(['c','a','c'])
                sage: elt = S(['c','c','a'])
                sage: elt.check()
            """
            if self not in self.parent():
                raise ValueError("Invalid permutation")

    def _repr_(self) -> str:
        """
        TESTS::

            sage: Permutations(['c','a','c'])
            Permutations of the multi-set ['c', 'a', 'c']
        """
        return "Permutations of the multi-set %s"%list(self.mset)

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: [ p for p in Permutations(['c','t','t'])] # indirect doctest
            [['c', 't', 't'], ['t', 'c', 't'], ['t', 't', 'c']]

        TESTS:

        The empty multiset::

            sage: list(sage.combinat.permutation.Permutations_mset([]))
            [[]]
        """
        mset = self.mset
        n = len(mset)
        from array import array
        mset_list = array('I', sorted(mset.index(x) for x in mset))

        yield self.element_class(self, map_to_list(mset_list, mset, n), check=False)

        if n <= 1:
            return

        while next_perm(mset_list):
            #Yield the permutation
            yield self.element_class(self, map_to_list(mset_list, mset, n), check=False)

    def cardinality(self):
        """
        Return the cardinality of the set.

        EXAMPLES::

            sage: Permutations([1,2,2]).cardinality()
            3
            sage: Permutations([1,1,2,2,2]).cardinality()
            10
        """
        lmset = list(self.mset)
        mset_list = [lmset.index(x) for x in lmset]
        d = {}
        for i in mset_list:
            d[i] = d.get(i, 0) + 1

        return ZZ(multinomial(d.values()))

    def rank(self, p):
        r"""
        Return the rank of ``p`` in lexicographic order.

        INPUT:

        - ``p`` -- a permutation of `M`

        ALGORITHM:

        The algorithm uses the recurrence from the solution to exercise 4 in
        [Knu2011]_, Section 7.2.1.2:

        .. MATH::

            \mathrm{rank}(p_1 \ldots p_n) =
            \mathrm{rank}(p_2 \ldots p_n)
            + \frac{1}{n} \genfrac{(}{)}{0pt}{0}{n}{n_1, \ldots, n_t}
                \sum_{j=1}^t n_j \left[ x_j < p_1 \right],

        where `x_j, n_j` are the distinct elements of `p` with their
        multiplicities, `n` is the sum of `n_1, \ldots, n_t`,
        `\genfrac{(}{)}{0pt}{1}{n}{n_1, \ldots, n_t}` is the multinomial
        coefficient `\frac{n!}{n_1! \ldots n_t!}`, and
        `\sum_{j=1}^t n_j \left[ x_j < p_1 \right]` means "the number of
        elements to the right of the first element that are less than the first
        element".

        EXAMPLES::

            sage: mset = [1, 1, 2, 3, 4, 5, 5, 6, 9]
            sage: p = Permutations(mset)
            sage: p.rank(list(sorted(mset)))
            0
            sage: p.rank(list(reversed(sorted(mset)))) == p.cardinality() - 1
            True
            sage: p.rank([3, 1, 4, 1, 5, 9, 2, 6, 5])
            30991

        TESTS::

            sage: from sage.combinat.permutation import Permutations_mset
            sage: p = Permutations_mset([])
            sage: p.rank([])
            0
            sage: p = Permutations_mset([1, 1, 2, 3, 4, 5, 5, 6, 9])
            sage: p.rank([1, 2, 3, 4, 5, 6, 9])
            Traceback (most recent call last):
            ...
            ValueError: Invalid permutation

        Try with greater-than-unity multiplicity in the least and greatest
        elements::

            sage: mset = list(range(10)) * 3
            sage: p = Permutations_mset(mset)
            sage: p.rank(list(sorted(mset)))
            0
            sage: p.rank(list(reversed(sorted(mset))))
            4386797336285844479999999
            sage: p.cardinality()
            4386797336285844480000000

        Should match ``StandardPermutations_n`` when `M` is the set
        `\{1, 2, \ldots, n\}`::

            sage: ps = Permutations(4)
            sage: pm = Permutations_mset(list(range(1, 5)))
            sage: ps.rank([2, 3, 1, 4]) == pm.rank([2, 3, 1, 4])
            True
        """
        self(p).check()
        m = {}
        r = 0
        for n in range(1, len(p)+1):
            # ``p1`` is the first element of ``p[-n:]`` (i.e., the last ``n``
            # elements of ``p``). ``m`` represents the multiset of ``p[-n:]`` in
            # the form element→count.
            p1 = p[-n]
            m[p1] = m.get(p1, 0) + 1
            r += multinomial(m.values()) * sum(nj for xj, nj in m.items() if xj < p1) // n
        return r

    def unrank(self, r):
        r"""
        Return the permutation of `M` having lexicographic rank ``r``.

        INPUT:

        - ``r`` -- an integer between ``0`` and ``self.cardinality()-1``
          inclusive

        ALGORITHM:

        The algorithm is adapted from the solution to exercise 4 in [Knu2011]_,
        Section 7.2.1.2.

        EXAMPLES::

            sage: mset = [1, 1, 2, 3, 4, 5, 5, 6, 9]
            sage: p = Permutations(mset)
            sage: p.unrank(30991)
            [3, 1, 4, 1, 5, 9, 2, 6, 5]
            sage: p.rank(p.unrank(10))
            10
            sage: p.unrank(0) == list(sorted(mset))
            True
            sage: p.unrank(p.cardinality()-1) == list(reversed(sorted(mset)))
            True

        TESTS::

            sage: from sage.combinat.permutation import Permutations_mset
            sage: p = Permutations_mset([])
            sage: p.unrank(0)
            []
            sage: p.unrank(1)
            Traceback (most recent call last):
            ...
            ValueError: r must be between 0 and 0 inclusive
            sage: p = Permutations_mset([1, 1, 2, 3, 4, 5, 5, 6, 9])
            sage: p.unrank(-1)
            Traceback (most recent call last):
            ...
            ValueError: r must be between 0 and 90719 inclusive
            sage: p.unrank(p.cardinality())
            Traceback (most recent call last):
            ...
            ValueError: r must be between 0 and 90719 inclusive

        Try with a cardinality that exceeds the precise integer range of a
        float::

            sage: mset = list(range(10)) * 3
            sage: p = Permutations_mset(mset)
            sage: p.unrank(p.rank(mset)) == mset
            True
            sage: p.unrank(p.cardinality()-1) == list(reversed(sorted(mset)))
            True

        Exhaustive check of roundtrip and lexicographic order for a single
        multiset::

            sage: p = Permutations_mset([2, 2, 3, 3, 3, 5, 5, 5, 5, 5])
            sage: prev = None
            sage: for rank, perm in enumerate(p):
            ....:     perm = tuple(perm)
            ....:     assert p.rank(perm) == rank, (rank, perm, p.rank(perm))
            ....:     assert tuple(p.unrank(rank)) == perm, (rank, perm, p.unrank(rank))
            ....:     assert prev is None or prev < perm
            ....:     prev = perm

        Should match ``StandardPermutations_n`` when `M` is the set
        `\{1, 2, \ldots, n\}`::

            sage: ps = Permutations(4)
            sage: pm = Permutations_mset(list(range(1, 5)))
            sage: ps.unrank(5) == pm.unrank(5)
            True
        """
        range_error = ValueError("r must be between %d and %d inclusive" % (0, self.cardinality()-1))
        if r < 0:
            raise range_error

        # ``mset`` is a working copy of the remaining multiset, initially set to
        # ``self.mset`` and always kept sorted. ``m`` represents ``mset`` in
        # element→count form. One element is removed from both representations
        # of the multiset in each loop iteration.
        mset = list(sorted(self.mset))
        m = {}
        for x in mset:
            m[x] = m.get(x, 0) + 1

        from bisect import bisect_left

        # We will build up ``p`` from left to right. For each element ``pi``, we
        # use the residual rank to compute how many elements to the right of
        # position ``i`` have a value less than ``pi``. We use that count of
        # elements to index what remains of ``mset`` to get the value of ``pi``.
        p = []
        while mset:
            n = len(mset)
            mult = multinomial(m.values())
            # Now ``r = ci * mult / n + rank(p[i+1:])``, where ``ci`` is the
            # number of elements remaining in ``mset`` that are less than
            # ``pi``. We know ``r`` and ``mult / n``, but not ``rank(p[i+1:])``
            # nor ``ci``, except that we know both are integers.
            ti = r * n // mult
            if ti >= n:
                raise range_error
            # ``ti`` is an upper bound for ``ci``. The true value is
            # ``ci = (r - rank(p[i+1:])) * n / mult``, but we do not know
            # ``rank(p[i+1:])`` yet. When ``pi`` has multiplicity 1, ``ti = ci``
            # exactly, because in that case ``rank(p[i+1:]) < mult / n``. When
            # ``pi`` has multiplicity greater than 1, ``ti`` is equal to the
            # number of elements that are less than ``pi``, *plus* some fraction
            # of those that are equal to ``pi``. We need to search left from
            # ``ti`` to find the lowest-indexed element equal to ``mset[ti]``,
            # to find ``ci`` which counts only elements that are strictly less
            # than ``pi``. We do a binary search in ``mset`` to the left of
            # position ``ti`` to find the leftmost index whose value is equal to
            # ``mset[ti]``. We do not have to look more places to the left than
            # the multiplicity of ``pi`` (``m[mset[ti]]``).
            ci = bisect_left(mset, mset[ti], max(0, ti + 1 - m[mset[ti]]), ti + 1)
            # Now that we know ``ci``, we can set ``r = rank(p[i+1:])`` for the
            # next iteration.
            r -= ci * mult // n
            # Because ``mset`` is sorted, the element that is greater than
            # ``ci`` other elements is at index ``ci``. Add that element to the
            # permutation and remove it from the multiset for the next
            # iteration.
            pi = mset.pop(ci)
            m[pi] -= 1
            if m[pi] == 0:
                del m[pi]
            p.append(pi)

        if r > 0:
            raise range_error

        return p


class Permutations_set(Permutations):
    """
    Permutations of an arbitrary given finite set.

    Here, a "permutation of a finite set `S`" means a list of the
    elements of `S` in which every element of `S` occurs exactly
    once. This is not to be confused with bijections from `S` to
    `S`, which are also often called permutations in literature.
    """
    @staticmethod
    def __classcall_private__(cls, s):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(['c','a','t'])
            sage: S2 = Permutations(('c','a','t'))
            sage: S1 is S2
            True
        """
        return super(Permutations_set, cls).__classcall__(cls, tuple(s))

    def __init__(self, s):
        """
        TESTS::

            sage: S = Permutations(['c','a','t'])
            sage: TestSuite(S).run()
        """
        Permutations.__init__(self, category=FiniteEnumeratedSets())
        self._set = s

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: p = Permutations([4,-1,'cthulhu'])
            sage: [4,-1,'cthulhu'] in p
            True
            sage: [] in p
            False
            sage: [4,'cthulhu','fhtagn'] in p
            False
            sage: [4,'cthulhu',4,-1] in p
            False
            sage: [-1,'cthulhu',4] in p
            True
        """
        s = list(self._set)
        if len(x) != len(s):
            return False
        for i in x:
            if i in s:
                s.remove(i)
            else:
                return False
        return True

    def _repr_(self) -> str:
        """
        TESTS::

            sage: Permutations(['c','a','t'])
            Permutations of the set ['c', 'a', 't']
        """
        return "Permutations of the set %s"%list(self._set)

    class Element(ClonableArray):
        """
        A permutation of an arbitrary set.
        """
        def check(self):
            """
            Verify that ``self`` is a valid permutation of the underlying
            set.

            EXAMPLES::

                sage: S = Permutations(['c','a','t'])
                sage: elt = S(['t','c','a'])
                sage: elt.check()
            """
            if self not in self.parent():
                raise ValueError("Invalid permutation")

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: S = Permutations(['c','a','t'])
            sage: S.list()
            [['c', 'a', 't'],
             ['c', 't', 'a'],
             ['a', 'c', 't'],
             ['a', 't', 'c'],
             ['t', 'c', 'a'],
             ['t', 'a', 'c']]
        """
        for p in itertools.permutations(self._set, len(self._set)):
            yield self.element_class(self, p)

    def cardinality(self) -> Integer:
        """
        Return the cardinality of the set.

        EXAMPLES::

            sage: Permutations([1,2,3]).cardinality()
            6
        """
        return factorial(len(self._set))

    def random_element(self):
        """
        EXAMPLES::

            sage: s = Permutations([1,2,3]).random_element()
            sage: s.parent() is Permutations([1,2,3])
            True
        """
        return sample(self._set, len(self._set))


class Permutations_msetk(Permutations_mset):
    """
    Length-`k` partial permutations of a multiset.

    A length-`k` partial permutation of a multiset `M` is
    represented by a list of length `k` whose entries are
    elements of `M`, appearing in the list with a multiplicity
    not higher than their respective multiplicity in `M`.
    """
    @staticmethod
    def __classcall__(cls, mset, k):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(['c','a','c'], 2)
            sage: S2 = Permutations(('c','a','c'), 2)
            sage: S1 is S2
            True
        """
        return super(Permutations_msetk, cls).__classcall__(cls, tuple(mset), k)

    def __init__(self, mset, k):
        """
        TESTS::

            sage: P = Permutations([1,2,2],2)
            sage: TestSuite(P).run()
        """
        Permutations_mset.__init__(self, mset)
        self._k = k

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: p = Permutations([1,2,2],2)
            sage: [1,2,2] in p
            False
            sage: [2,2] in p
            True
            sage: [1,1] in p
            False
            sage: [2,1] in p
            True
        """
        if len(x) != self._k:
            return False
        s = list(self.mset)
        for i in x:
            if i in s:
                s.remove(i)
            else:
                return False
        return True

    def _repr_(self):
        """
        TESTS::

            sage: Permutations([1,2,2],2)
            Permutations of the multi-set [1, 2, 2] of length 2
        """
        return "Permutations of the multi-set %s of length %s" % (list(self.mset), self._k)

    def cardinality(self):
        """
        Return the cardinality of the set.

        EXAMPLES::

            sage: Permutations([1,2,2],2).cardinality()
            3
        """
        return ZZ.sum(1 for z in self)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations([1,2,2],2).list()
            [[1, 2], [2, 1], [2, 2]]
        """
        mset = self.mset
        lmset = list(mset)
        mset_list = [lmset.index(x) for x in lmset]
        indices = libgap.Arrangements(mset_list, self._k).sage()
        for ktuple in indices:
            yield self.element_class(self, [lmset[x] for x in ktuple])


class Permutations_setk(Permutations_set):
    """
    Length-`k` partial permutations of an arbitrary given finite set.

    Here, a "length-`k` partial permutation of a finite set `S`" means
    a list of length `k` whose entries are pairwise distinct and all
    belong to `S`.
    """
    @staticmethod
    def __classcall_private__(cls, s, k):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(['c','a','t'], 2)
            sage: S2 = Permutations(('c','a','t'), 2)
            sage: S1 is S2
            True
        """
        return super(Permutations_setk, cls).__classcall__(cls, tuple(s), k)

    def __init__(self, s, k):
        """
        TESTS::

            sage: P = Permutations([1,2,4],2)
            sage: TestSuite(P).run()
        """
        Permutations_set.__init__(self, s)
        self._k = k

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: p = Permutations([1,2,4],2)
            sage: [1,2,4] in p
            False
            sage: [2,2] in p
            False
            sage: [1,4] in p
            True
            sage: [2,1] in p
            True
        """
        if len(x) != self._k:
            return False
        s = list(self._set)
        return all(i in s for i in x) and len(set(x)) == len(x)

    def _repr_(self):
        """
        TESTS::

            sage: repr(Permutations([1,2,4],2))
            'Permutations of the set [1, 2, 4] of length 2'
        """
        return "Permutations of the set %s of length %s"%(list(self._set), self._k)

    def __iter__(self):
        """
        EXAMPLES::

            sage: [i for i in Permutations([1,2,4],2)]
            [[1, 2], [1, 4], [2, 1], [2, 4], [4, 1], [4, 2]]
        """
        for perm in itertools.permutations(self._set, int(self._k)):
            yield self.element_class(self, perm)

    def random_element(self):
        """
        EXAMPLES::

            sage: s = Permutations([1,2,4], 2).random_element()
            sage: s in Permutations([1,2,4], 2)
            True
        """
        return sample(self._set, self._k)

##################################
# Arrangements

class Arrangements(Permutations):
    r"""
    An arrangement of a multiset ``mset`` is an ordered selection
    without repetitions. It is represented by a list that contains
    only elements from ``mset``, but maybe in a different order.

    ``Arrangements`` returns the combinatorial class of
    arrangements of the multiset ``mset`` that contain ``k`` elements.

    EXAMPLES::

        sage: mset = [1,1,2,3,4,4,5]
        sage: Arrangements(mset,2).list()
        [[1, 1],
         [1, 2],
         [1, 3],
         [1, 4],
         [1, 5],
         [2, 1],
         [2, 3],
         [2, 4],
         [2, 5],
         [3, 1],
         [3, 2],
         [3, 4],
         [3, 5],
         [4, 1],
         [4, 2],
         [4, 3],
         [4, 4],
         [4, 5],
         [5, 1],
         [5, 2],
         [5, 3],
         [5, 4]]
         sage: Arrangements(mset,2).cardinality()
         22
         sage: Arrangements( ["c","a","t"], 2 ).list()
         [['c', 'a'], ['c', 't'], ['a', 'c'], ['a', 't'], ['t', 'c'], ['t', 'a']]
         sage: Arrangements( ["c","a","t"], 3 ).list()
         [['c', 'a', 't'],
          ['c', 't', 'a'],
          ['a', 'c', 't'],
          ['a', 't', 'c'],
          ['t', 'c', 'a'],
          ['t', 'a', 'c']]
    """
    @staticmethod
    def __classcall_private__(cls, mset, k):
        """
        Return the correct parent.

        EXAMPLES::

            sage: A1 = Arrangements( ["c","a","t"], 2)
            sage: A2 = Arrangements( ("c","a","t"), 2)
            sage: A1 is A2
            True
        """
        mset = tuple(mset)
        if [mset.index(_) for _ in mset] == list(range(len(mset))):
            return Arrangements_setk(mset, k)
        return Arrangements_msetk(mset, k)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: A = Arrangements([1,1,2,3,4,4,5], 2)
            sage: A.cardinality()
            22
        """
        one = ZZ.one()
        return sum(one for p in self)


class Arrangements_msetk(Arrangements, Permutations_msetk):
    r"""
    Arrangements of length `k` of a multiset `M`.
    """
    def _repr_(self):
        """
        TESTS::

            sage: Arrangements([1,2,2],2)
            Arrangements of the multi-set [1, 2, 2] of length 2
        """
        return "Arrangements of the multi-set %s of length %s"%(list(self.mset),self._k)


class Arrangements_setk(Arrangements, Permutations_setk):
    r"""
    Arrangements of length `k` of a set `S`.
    """
    def _repr_(self):
        """
        TESTS::

            sage: Arrangements([1,2,3],2)
            Arrangements of the set [1, 2, 3] of length 2
        """
        return "Arrangements of the set %s of length %s"%(list(self._set),self._k)

###############################################################
# Standard permutations

class StandardPermutations_all(Permutations):
    """
    All standard permutations.
    """
    def __init__(self):
        """
        TESTS::

            sage: SP = Permutations()
            sage: TestSuite(SP).run()
        """
        cat = InfiniteEnumeratedSets() & SetsWithGrading()
        Permutations.__init__(self, category=cat)

    def _repr_(self):
        """
        TESTS::

            sage: Permutations()
            Standard permutations
        """
        return "Standard permutations"

    def __contains__(self, x):
        """
        TESTS::

            sage: [] in Permutations()
            True
            sage: [1] in Permutations()
            True
            sage: [2] in Permutations()
            False
            sage: [1,2] in Permutations()
            True
            sage: [2,1] in Permutations()
            True
            sage: [1,2,2] in Permutations()
            False
            sage: [3,1,5,2] in Permutations()
            False
            sage: [3,4,1,5,2] in Permutations()
            True
        """
        if isinstance(x, Permutation):
            return True
        elif isinstance(x, list):
            s = sorted(x[:])
            if s != list(range(1, len(x)+1)):
                return False
            return True
        else:
            return False

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: it = iter(Permutations())
            sage: [next(it) for i in range(10)]
            [[], [1], [1, 2], [2, 1], [1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        n = 0
        while True:
            for p in itertools.permutations(range(1, n + 1)):
                yield self.element_class(self, p)
            n += 1

    def graded_component(self, n):
        """
        Return the graded component.

        EXAMPLES::

            sage: P = Permutations()
            sage: P.graded_component(4) == Permutations(4)
            True
        """
        return StandardPermutations_n(n)


class StandardPermutations_n_abstract(Permutations):
    r"""
    Abstract base class for subsets of permutations of the
    set `\{1, 2, \ldots, n\}`.

    .. WARNING::

        Anything inheriting from this class should override the
        ``__contains__`` method.
    """
    def __init__(self, n, category=None):
        """
        TESTS:

        We skip the descent and reduced word methods because they do
        not respect the ordering for multiplication::

            sage: SP = Permutations(3)
            sage: TestSuite(SP).run(skip=['_test_reduced_word', '_test_has_descent'])

            sage: SP.options.mult='r2l'
            sage: TestSuite(SP).run()
            sage: SP.options._reset()
        """
        self.n = ZZ(n)
        if category is None:
            category = FiniteEnumeratedSets()
        Permutations.__init__(self, category=category)

    def _element_constructor_(self, x, check_input=True):
        """
        Construct an element of ``self`` from ``x``.

        TESTS::

            sage: P = Permutations(5)
            sage: P([2,3,1])
            [2, 3, 1, 4, 5]

            sage: G = SymmetricGroup(4)
            sage: P = Permutations(4)
            sage: x = G([4,3,1,2]); x
            (1,4,2,3)
            sage: P(x)
            [4, 3, 1, 2]
            sage: G(P(x))
            (1,4,2,3)

            sage: P = PermutationGroup([[(1,3,5),(2,4)],[(1,3)]])
            sage: x = P([(3,5),(2,4)]); x
            (2,4)(3,5)
            sage: Permutations(6)(SymmetricGroup(6)(x))
            [1, 4, 5, 2, 3, 6]
            sage: Permutations(6)(x)  # known bug
            [1, 4, 5, 2, 3, 6]
        """
        if len(x) < self.n:
            x = list(x) + list(range(len(x) + 1, self.n + 1))
        return self.element_class(self, x, check_input=check_input)

    def __contains__(self, x):
        """
        TESTS::

            sage: [] in Permutations(0)
            True
            sage: [1,2] in Permutations(2)
            True
            sage: [1,2] in Permutations(3)
            False
            sage: [3,2,1] in Permutations(3)
            True
        """
        return Permutations.__contains__(self, x) and len(x) == self.n

class StandardPermutations_n(StandardPermutations_n_abstract):
    r"""
    Permutations of the set `\{1, 2, \ldots, n\}`.

    These are also called permutations of size `n`, or the elements
    of the `n`-th symmetric group.

    .. TODO::

        Have a :meth:`reduced_word` which works in both multiplication
        conventions.
    """
    def __init__(self, n):
        """
        Initialize ``self``.

        TESTS::

            sage: P = Permutations(5)
            sage: P.options.mult = 'r2l'
            sage: TestSuite(P).run(skip='_test_descents')
            sage: P.options._reset()
        """
        cat = FiniteWeylGroups().Irreducible() & FinitePermutationGroups()
        StandardPermutations_n_abstract.__init__(self, n, category=cat)

    def _repr_(self):
        """
        TESTS::

            sage: Permutations(3)
            Standard permutations of 3
        """
        return "Standard permutations of %s"%self.n

    def __iter__(self):
        """
        EXAMPLES::

            sage: [p for p in Permutations(0)]
            [[]]
            sage: [p for p in Permutations(3)]
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        for p in itertools.permutations(range(1, self.n + 1)):
            yield self.element_class(self, p)

    def _coerce_map_from_(self, G):
        """
        Return a coerce map or ``True`` if there exists a coerce map
        from ``G``.

        .. WARNING::

            The coerce maps between ``Permutations(n)`` and
            ``SymmetricGroup(n)`` exist, but do not respect the
            multiplication when the global variable
            ``Permutations.options.mult`` (see
            :meth:`sage.combinat.permutation.Permutations.options` )
            is set to ``'r2l'``. (Indeed, the multiplication
            in ``Permutations(n)`` depends on this global
            variable, while the multiplication in
            ``SymmetricGroup(n)`` does not.)

        EXAMPLES::

            sage: P = Permutations(6)
            sage: P.has_coerce_map_from(SymmetricGroup(6))
            True
            sage: P.has_coerce_map_from(SymmetricGroup(5))
            True
            sage: P.has_coerce_map_from(SymmetricGroup(7))
            False
            sage: P.has_coerce_map_from(Permutations(5))
            True
            sage: P.has_coerce_map_from(Permutations(7))
            False
        """
        if isinstance(G, SymmetricGroup):
            D = G.domain()
            if len(D) > self.n or list(D) != list(range(1, len(D)+1)):
                return False
            return self._from_permutation_group_element
        if isinstance(G, StandardPermutations_n) and G.n <= self.n:
            return True
        return super(StandardPermutations_n, self)._coerce_map_from_(G)

    def _from_permutation_group_element(self, x):
        """
        Return an element of ``self`` from a permutation group element.

        TESTS::

            sage: P = Permutations(4)
            sage: G = SymmetricGroup(4)
            sage: x = G([4,3,1,2])
            sage: P._from_permutation_group_element(x)
            [4, 3, 1, 2]
        """
        return self(x.domain())

    def as_permutation_group(self):
        """
        Return ``self`` as a permutation group.

        EXAMPLES::

            sage: P = Permutations(4)
            sage: PG = P.as_permutation_group()
            sage: PG
            Symmetric group of order 4! as a permutation group

            sage: G = SymmetricGroup(4)
            sage: PG is G
            True
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        return SymmetricGroup(self.n)

    def identity(self):
        r"""
        Return the identity permutation of size `n`.

        EXAMPLES::

            sage: Permutations(4).identity()
            [1, 2, 3, 4]
            sage: Permutations(0).identity()
            []
        """
        return self.element_class(self, range(1,self.n+1))

    one = identity

    def unrank(self, r):
        """
        EXAMPLES::

            sage: SP3 = Permutations(3)
            sage: l = list(map(SP3.unrank, range(6)))
            sage: l == SP3.list()
            True
            sage: SP0 = Permutations(0)
            sage: l = list(map(SP0.unrank, range(1)))
            sage: l == SP0.list()
            True
        """
        if r >= factorial(self.n) or r < 0:
            raise ValueError
        else:
            return from_rank(self.n, r)

    def rank(self, p=None):
        """
        Return the rank of ``self`` or ``p`` depending on input.

        If a permutation ``p`` is given, return the rank of ``p``
        in ``self``. Otherwise return the dimension of the
        underlying vector space spanned by the (simple) roots.

        EXAMPLES::

            sage: P = Permutations(5)
            sage: P.rank()
            4

            sage: SP3 = Permutations(3)
            sage: list(map(SP3.rank, SP3))
            [0, 1, 2, 3, 4, 5]
            sage: SP0 = Permutations(0)
            sage: list(map(SP0.rank, SP0))
            [0]
        """
        if p is None:
            return self.n - 1
        if p in self:
            return Permutation(p).rank()
        raise ValueError("p not in self")

    def random_element(self):
        """
        EXAMPLES::

            sage: s = Permutations(4).random_element(); s  # random
            [1, 2, 4, 3]
            sage: s in Permutations(4)
            True
        """
        return self.element_class(self, sample(range(1,self.n+1), self.n))

    def cardinality(self):
        """
        Return the number of permutations of size `n`, which is `n!`.

        EXAMPLES::

            sage: Permutations(0).cardinality()
            1
            sage: Permutations(3).cardinality()
            6
            sage: Permutations(4).cardinality()
            24
        """
        return factorial(self.n)

    def degree(self):
        """
        Return the degree of ``self``.

        This is the cardinality `n` of the set ``self`` acts on.

        EXAMPLES::

            sage: Permutations(0).degree()
            0
            sage: Permutations(1).degree()
            1
            sage: Permutations(5).degree()
            5
        """
        return self.n

    def degrees(self):
        """
        Return the degrees of ``self``.

        These are the degrees of the fundamental invariants of the
        ring of polynomial invariants.

        EXAMPLES::

            sage: Permutations(3).degrees()
            (2, 3)
            sage: Permutations(7).degrees()
            (2, 3, 4, 5, 6, 7)
        """
        return tuple(Integer(i) for i in range(2, self.n+1))

    def codegrees(self):
        """
        Return the codegrees of ``self``.

        EXAMPLES::

            sage: Permutations(3).codegrees()
            (0, 1)
            sage: Permutations(7).codegrees()
            (0, 1, 2, 3, 4, 5)
        """
        return tuple(Integer(i) for i in range(self.n-1))

    def element_in_conjugacy_classes(self, nu):
        r"""
        Return a permutation with cycle type ``nu``.

        If the size of ``nu`` is smaller than the size of permutations in
        ``self``, then some fixed points are added.

        EXAMPLES::

            sage: PP = Permutations(5)
            sage: PP.element_in_conjugacy_classes([2,2])
            [2, 1, 4, 3, 5]
        """
        from sage.combinat.partition import Partition
        nu = Partition(nu)
        if nu.size() > self.n:
            raise ValueError("The size of the partition (=%s) should be lower or equal"
                             " to the size of the permutations (=%s)"%(nu.size,self.n))
        l = []
        i = 0
        for nui in nu:
            for j in range(nui-1):
                l.append(i+j+2)
            l.append(i+1)
            i += nui
        for i in range(nu.size(), self.n):
            l.append(i+1)
        return self.element_class(self, l)

    def conjugacy_classes_representatives(self):
        r"""
        Return a complete list of representatives of conjugacy classes
        in ``self``.

        Let `S_n` be the symmetric group on `n` letters. The conjugacy
        classes are indexed by partitions `\lambda` of `n`. The ordering
        of the conjugacy classes is reverse lexicographic order of
        the partitions.

        EXAMPLES::

            sage: G = Permutations(5)
            sage: G.conjugacy_classes_representatives()
            [[1, 2, 3, 4, 5],
             [2, 1, 3, 4, 5],
             [2, 1, 4, 3, 5],
             [2, 3, 1, 4, 5],
             [2, 3, 1, 5, 4],
             [2, 3, 4, 1, 5],
             [2, 3, 4, 5, 1]]

        TESTS:

        Check some border cases::

            sage: S = Permutations(0)
            sage: S.conjugacy_classes_representatives()
            [[]]
            sage: S = Permutations(1)
            sage: S.conjugacy_classes_representatives()
            [[1]]
        """
        from sage.combinat.partition import Partitions_n
        return [ self.element_in_conjugacy_classes(la)
                 for la in reversed(Partitions_n(self.n)) ]

    def conjugacy_classes_iterator(self):
        """
        Iterate over the conjugacy classes of ``self``.

        EXAMPLES::

            sage: G = Permutations(4)
            sage: list(G.conjugacy_classes_iterator()) == G.conjugacy_classes()
            True
        """
        from sage.combinat.partition import Partitions_n
        from sage.groups.perm_gps.symgp_conjugacy_class import PermutationsConjugacyClass
        for la in reversed(Partitions_n(self.n)):
            yield PermutationsConjugacyClass(self, la)

    def conjugacy_classes(self):
        """
        Return a list of the conjugacy classes of ``self``.

        EXAMPLES::

            sage: G = Permutations(4)
            sage: G.conjugacy_classes()
            [Conjugacy class of cycle type [1, 1, 1, 1] in Standard permutations of 4,
             Conjugacy class of cycle type [2, 1, 1] in Standard permutations of 4,
             Conjugacy class of cycle type [2, 2] in Standard permutations of 4,
             Conjugacy class of cycle type [3, 1] in Standard permutations of 4,
             Conjugacy class of cycle type [4] in Standard permutations of 4]
        """
        return list(self.conjugacy_classes_iterator())

    def conjugacy_class(self, g):
        r"""
        Return the conjugacy class of ``g`` in ``self``.

        INPUT:

        - ``g`` -- a partition or an element of ``self``

        EXAMPLES::

            sage: G = Permutations(5)
            sage: g = G([2,3,4,1,5])
            sage: G.conjugacy_class(g)
            Conjugacy class of cycle type [4, 1] in Standard permutations of 5
            sage: G.conjugacy_class(Partition([2, 1, 1, 1]))
            Conjugacy class of cycle type [2, 1, 1, 1] in Standard permutations of 5
        """
        from sage.groups.perm_gps.symgp_conjugacy_class import PermutationsConjugacyClass
        return PermutationsConjugacyClass(self, g)

    def algebra(self, base_ring, category=None):
        """
        Return the symmetric group algebra associated to ``self``.

        INPUT:

        - ``base_ring`` -- a ring
        - ``category`` -- a category (default: the category of ``self``)

        EXAMPLES::

            sage: P = Permutations(4)
            sage: A = P.algebra(QQ); A
            Symmetric group algebra of order 4 over Rational Field

            sage: A.category()
            Join of Category of coxeter group algebras over Rational Field
                and Category of finite group algebras over Rational Field
                and Category of finite dimensional cellular algebras with basis over Rational Field
            sage: A = P.algebra(QQ, category=Monoids())
            sage: A.category()
            Category of finite dimensional cellular monoid algebras over Rational Field
        """
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        return SymmetricGroupAlgebra(base_ring, self, category=category)

    @cached_method
    def index_set(self):
        r"""
        Return the index set for the descents of the symmetric group ``self``.

        This is `\{ 1, 2, \ldots, n-1 \}`, where ``self`` is `S_n`.

        EXAMPLES::

            sage: P = Permutations(8)
            sage: P.index_set()
            (1, 2, 3, 4, 5, 6, 7)
        """
        return tuple(range(1, self.n))

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        The symmetric group `S_n` is a Coxeter group of type `A_{n-1}`.

        EXAMPLES::

            sage: A = SymmetricGroup([2,3,7]); A.cartan_type()
            ['A', 2]
            sage: A = SymmetricGroup([]); A.cartan_type()
            ['A', 0]
        """
        from sage.combinat.root_system.cartan_type import CartanType
        return CartanType(['A', max(self.n - 1,0)])

    def simple_reflection(self, i):
        r"""
        For `i` in the index set of ``self`` (that is, for `i` in
        `\{ 1, 2, \ldots, n-1 \}`, where ``self`` is `S_n`), this
        returns the elementary transposition `s_i = (i,i+1)`.

        EXAMPLES::

            sage: P = Permutations(4)
            sage: P.simple_reflection(2)
            [1, 3, 2, 4]
            sage: P.simple_reflections()
            Finite family {1: [2, 1, 3, 4], 2: [1, 3, 2, 4], 3: [1, 2, 4, 3]}
        """
        g = list(range(1, self.n+1))
        g[i-1] = i+1
        g[i] = i
        return self.element_class(self, g)

    class Element(Permutation):
        def has_left_descent(self, i, mult=None):
            r"""
            Check if ``i`` is a left descent of ``self``.

            A *left descent* of a permutation `\pi \in S_n` means an index
            `i \in \{ 1, 2, \ldots, n-1 \}` such that
            `s_i \circ \pi` has smaller length than `\pi`. Thus, a left
            descent of `\pi` is an index `i \in \{ 1, 2, \ldots, n-1 \}`
            satisfying `\pi^{-1}(i) > \pi^{-1}(i+1)`.

            .. WARNING::

                The methods :meth:`descents` and :meth:`idescents` behave
                differently than their Weyl group counterparts. In
                particular, the indexing is 0-based. This could lead to
                errors. Instead, construct the descent set as in the example.

            .. WARNING::

                This ignores the multiplication convention in order
                to be consistent with other Coxeter operations in
                permutations (e.g., computing :meth:`reduced_word`).

            EXAMPLES::

                sage: P = Permutations(4)
                sage: x = P([3, 2, 4, 1])
                sage: (~x).descents()
                [1, 2]
                sage: [i for i in P.index_set() if x.has_left_descent(i)]
                [1, 2]

            TESTS::

                sage: P = Permutations(4)
                sage: x = P([3, 2, 4, 1])
                sage: x.has_left_descent(2, mult='l2r')
                doctest:warning
                ...
                DeprecationWarning: The mult option is deprecated and ignored.
                See https://trac.sagemath.org/27467 for details.
                True
                sage: x.has_left_descent(2, mult='r2l')
                True
            """
            if mult is not None:
                from sage.misc.superseded import deprecation
                deprecation(27467, "The mult option is deprecated and ignored.")
            for val in self._list:
                if val == i:
                    return False
                if val == i + 1:
                    return True

        def has_right_descent(self, i, mult=None):
            r"""
            Check if ``i`` is a right descent of ``self``.

            A *right descent* of a permutation `\pi \in S_n` means an index
            `i \in \{ 1, 2, \ldots, n-1 \}` such that
            `\pi \circ s_i` has smaller length than `\pi`. Thus, a right
            descent of `\pi` is an index `i \in \{ 1, 2, \ldots, n-1 \}`
            satisfying `\pi(i) > \pi(i+1)`.

            .. WARNING::

                The methods :meth:`descents` and :meth:`idescents` behave
                differently than their Weyl group counterparts. In
                particular, the indexing is 0-based. This could lead to
                errors. Instead, construct the descent set as in the example.

            .. WARNING::

                This ignores the multiplication convention in order
                to be consistent with other Coxeter operations in
                permutations (e.g., computing :meth:`reduced_word`).

            EXAMPLES::

                sage: P = Permutations(4)
                sage: x = P([3, 2, 4, 1])
                sage: x.descents()
                [1, 3]
                sage: [i for i in P.index_set() if x.has_right_descent(i)]
                [1, 3]

            TESTS::

                sage: P = Permutations(4)
                sage: x = P([3, 2, 4, 1])
                sage: x.has_right_descent(3, mult='l2r')
                doctest:warning
                ...
                DeprecationWarning: The mult option is deprecated and ignored.
                See https://trac.sagemath.org/27467 for details.
                True
                sage: x.has_right_descent(3, mult='r2l')
                True
            """
            if mult is not None:
                from sage.misc.superseded import deprecation
                deprecation(27467, "The mult option is deprecated and ignored.")
            return self[i-1] > self[i]

        def __mul__(self, other):
            r"""
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: P = Permutations(4)
                sage: P.simple_reflection(1) * Permutation([6,5,4,3,2,1])
                [5, 6, 4, 3, 2, 1]
                sage: Permutations.options.mult='r2l'
                sage: P.simple_reflection(1) * Permutation([6,5,4,3,2,1])
                [6, 5, 4, 3, 1, 2]
                sage: Permutations.options.mult='l2r'
            """
            if not isinstance(other, StandardPermutations_n.Element):
                return Permutation.__mul__(self, other)
            if other.parent() is not self.parent():
                # They have different parents (but both are (like) Permutations of n)
                mul_order = self.parent().options.mult
                if mul_order == 'l2r':
                    p = right_action_product(self._list, other._list)
                elif mul_order == 'r2l':
                    p = left_action_product(self._list, other._list)
                return Permutations(len(p))(p)
            # They have the same parent
            return self._mul_(other)

        def _mul_(self, other):
            r"""
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: P = Permutations(4)
                sage: P.prod(P.group_generators()).parent() is P
                True
            """
            mul_order = self.parent().options.mult
            if mul_order == 'l2r':
                p = right_action_same_n(self._list, other._list)
            elif mul_order == 'r2l':
                p = left_action_same_n(self._list, other._list)
            return self.__class__(self.parent(), p)

        @combinatorial_map(order=2, name='inverse')
        def inverse(self):
            r"""
            Return the inverse of ``self``.

            EXAMPLES::

                sage: P = Permutations(4)
                sage: w0 = P([4,3,2,1])
                sage: w0.inverse() == w0
                True
                sage: w0.inverse().parent() is P
                True
                sage: P([3,2,4,1]).inverse()
                [4, 2, 1, 3]
            """
            w = list(range(len(self)))
            for i,j in enumerate(self):
                w[j-1] = i+1
            return self.__class__(self.parent(), w)

        __invert__ = inverse

        def apply_simple_reflection_left(self, i):
            r"""
            Return ``self`` multiplied by the simple reflection ``s[i]``
            on the left.

            This acts by switching the entries in positions `i` and `i+1`.

            .. WARNING::

                This ignores the multiplication convention in order
                to be consistent with other Coxeter operations in
                permutations (e.g., computing :meth:`reduced_word`).

            EXAMPLES::

                sage: W = Permutations(3)
                sage: w = W([2,3,1])
                sage: w.apply_simple_reflection_left(1)
                [1, 3, 2]
                sage: w.apply_simple_reflection_left(2)
                [3, 2, 1]
            """
            s = self.parent().simple_reflections()[i]
            p = right_action_same_n(self._list, s._list)
            return self.__class__(self.parent(), p)

        def apply_simple_reflection_right(self, i):
            r"""
            Return ``self`` multiplied by the simple reflection ``s[i]``
            on the right.

            This acts by switching the entries `i` and `i+1`.

            .. WARNING::

                This ignores the multiplication convention in order
                to be consistent with other Coxeter operations in
                permutations (e.g., computing :meth:`reduced_word`).

            EXAMPLES::

                sage: W = Permutations(3)
                sage: w = W([2,3,1])
                sage: w.apply_simple_reflection_right(1)
                [3, 2, 1]
                sage: w.apply_simple_reflection_right(2)
                [2, 1, 3]
            """
            s = self.parent().simple_reflections()[i]
            p = left_action_same_n(self._list, s._list)
            return self.__class__(self.parent(), p)

#############################
# Constructing Permutations #
#############################

# TODO: Make this a coercion
def from_permutation_group_element(pge, parent=None):
    """
    Return a :class:`Permutation` given a :class:`PermutationGroupElement`
    ``pge``.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: pge = PermutationGroupElement([(1,2),(3,4)])
        sage: permutation.from_permutation_group_element(pge)
        [2, 1, 4, 3]
    """
    if not isinstance(pge, PermutationGroupElement):
        raise TypeError("pge (= %s) must be a PermutationGroupElement"%pge)

    if parent is None:
        parent = Permutations( len(pge.domain()) )

    return parent(pge.domain())

def from_rank(n, rank):
    r"""
    Return the permutation of the set `\{1,...,n\}` with lexicographic
    rank ``rank``. This is the permutation whose Lehmer code is the
    factoradic representation of ``rank``. In particular, the
    permutation with rank `0` is the identity permutation.

    The permutation is computed without iterating through all of the
    permutations with lower rank. This makes it efficient for large
    permutations.

    .. NOTE::

        The variable ``rank`` is not checked for being in the interval
        from `0` to `n! - 1`. When outside this interval, it acts as
        its residue modulo `n!`.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: Permutation([3, 6, 5, 4, 2, 1]).rank()
        359
        sage: [permutation.from_rank(3, i) for i in range(6)]
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        sage: Permutations(6)[10]
        [1, 2, 4, 6, 3, 5]
        sage: permutation.from_rank(6,10)
        [1, 2, 4, 6, 3, 5]
    """
    #Find the factoradic of rank
    factoradic = [None] * n
    for j in range(1,n+1):
        factoradic[n-j] = Integer(rank % j)
        rank = int(rank) // j

    return from_lehmer_code(factoradic, Permutations(n))

def from_inversion_vector(iv, parent=None):
    r"""
    Return the permutation corresponding to inversion vector ``iv``.

    See `~sage.combinat.permutation.Permutation.to_inversion_vector`
    for a definition of the inversion vector of a permutation.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_inversion_vector([3,1,0,0,0])
        [3, 2, 4, 1, 5]
        sage: permutation.from_inversion_vector([2,3,6,4,0,2,2,1,0])
        [5, 9, 1, 8, 2, 6, 4, 7, 3]
        sage: permutation.from_inversion_vector([0])
        [1]
        sage: permutation.from_inversion_vector([])
        []
    """
    p = iv[:]
    open_spots = list(range(len(iv)))
    for i,ivi in enumerate(iv):
        p[open_spots.pop(ivi)] = i+1

    if parent is None:
        parent = Permutations()
    return parent(p)

def from_cycles(n, cycles, parent=None):
    r"""
    Return the permutation in the `n`-th symmetric group whose decomposition
    into disjoint cycles is ``cycles``.

    This function checks that its input is correct (i.e. that the cycles are
    disjoint and their elements integers among `1...n`). It raises an exception
    otherwise.

    .. WARNING::

        It assumes that the elements are of ``int`` type.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_cycles(4, [[1,2]])
        [2, 1, 3, 4]
        sage: permutation.from_cycles(4, [[1,2,4]])
        [2, 4, 3, 1]
        sage: permutation.from_cycles(10, [[3,1],[4,5],[6,8,9]])
        [3, 2, 1, 5, 4, 8, 7, 9, 6, 10]
        sage: permutation.from_cycles(10, ((2, 5), (6, 1, 3)))
        [3, 5, 6, 4, 2, 1, 7, 8, 9, 10]
        sage: permutation.from_cycles(4, [])
        [1, 2, 3, 4]
        sage: permutation.from_cycles(4, [[]])
        [1, 2, 3, 4]
        sage: permutation.from_cycles(0, [])
        []

    Bad input (see :trac:`13742`)::

        sage: Permutation("(-12,2)(3,4)")
        Traceback (most recent call last):
        ...
        ValueError: All elements should be strictly positive integers, and I just found a non-positive one.
        sage: Permutation("(1,2)(2,4)")
        Traceback (most recent call last):
        ...
        ValueError: An element appears twice. It should not.
        sage: permutation.from_cycles(4, [[1,18]])
        Traceback (most recent call last):
        ...
        ValueError: You claimed that this was a permutation on 1...4 but it contains 18
    """
    if parent is None:
        parent = Permutations(n)

    p = list(range(1, n+1))

    # Is it really a permutation on 1...n ?
    flattened_and_sorted = []
    for c in cycles:
        flattened_and_sorted.extend(c)
    flattened_and_sorted.sort()

    # Empty input
    if not flattened_and_sorted:
        return parent(p, check_input=False)

    # Only positive elements
    if int(flattened_and_sorted[0]) < 1:
        raise ValueError("All elements should be strictly positive "
                         "integers, and I just found a non-positive one.")

    # Really smaller or equal to n ?
    if flattened_and_sorted[-1] > n:
        raise ValueError("You claimed that this was a permutation on 1..."+
                         str(n)+" but it contains "+str(flattened_and_sorted[-1]))

    # Disjoint cycles ?
    previous = flattened_and_sorted[0] - 1
    for i in flattened_and_sorted:
        if i == previous:
            raise ValueError("An element appears twice. It should not.")
        else:
            previous = i

    for cycle in cycles:
        if not cycle:
            continue
        first = cycle[0]
        for i in range(len(cycle)-1):
            p[cycle[i]-1] = cycle[i+1]
        p[cycle[-1]-1] = first

    return parent(p, check_input=False)

def from_lehmer_code(lehmer, parent=None):
    r"""
    Return the permutation with Lehmer code ``lehmer``.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: lc = Permutation([2,1,5,4,3]).to_lehmer_code(); lc
        [1, 0, 2, 1, 0]
        sage: permutation.from_lehmer_code(lc)
        [2, 1, 5, 4, 3]
    """
    p = []
    open_spots = list(range(1,len(lehmer)+1))
    for ivi in lehmer:
        p.append(open_spots.pop(ivi))

    if parent is None:
        parent = Permutations()
    return parent(p)

def from_reduced_word(rw, parent=None):
    r"""
    Return the permutation corresponding to the reduced word ``rw``.

    See
    :meth:`~sage.combinat.permutation.Permutation.reduced_words` for
    a definition of reduced words and the convention on the order of
    multiplication used.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_reduced_word([3,2,3,1,2,3,1])
        [3, 4, 2, 1]
        sage: permutation.from_reduced_word([])
        []
    """
    if parent is None:
        parent = Permutations()

    if not rw:
        return parent([])

    p = [i+1 for i in range(max(rw)+1)]

    for i in rw:
        (p[i-1], p[i]) = (p[i], p[i-1])

    return parent(p)

def bistochastic_as_sum_of_permutations(M, check = True):
    r"""
    Return the positive sum of permutations corresponding to
    the bistochastic matrix ``M``.

    A stochastic matrix is a matrix with nonnegative real entries such that the
    sum of the elements of any row is equal to `1`. A bistochastic matrix is a
    stochastic matrix whose transpose matrix is also stochastic ( there are
    conditions both on the rows and on the columns ).

    According to the Birkhoff-von Neumann Theorem, any bistochastic matrix
    can be written as a convex combination of permutation matrices, which also
    means that the polytope of bistochastic matrices is integer.

    As a non-bistochastic matrix can obviously not be written as a convex
    combination of permutations, this theorem is an equivalence.

    This function, given a bistochastic matrix, returns the corresponding
    decomposition.

    INPUT:

    - ``M`` -- A bistochastic matrix

    - ``check`` (boolean) -- set to ``True`` (default) to check
      that the matrix is indeed bistochastic

    OUTPUT:

    - An element of ``CombinatorialFreeModule``, which is a free `F`-module
      ( where `F` is the ground ring of the given matrix ) whose basis is
      indexed by the permutations.

    .. NOTE::

        - In this function, we just assume 1 to be any constant : for us a
          matrix `M` is bistochastic if there exists `c>0` such that `M/c`
          is bistochastic.

        - You can obtain a sequence of pairs ``(permutation,coeff)``, where
          ``permutation`` is a Sage ``Permutation`` instance, and ``coeff``
          its corresponding coefficient from the result of this function
          by applying the ``list`` function.

        - If you are interested in the matrix corresponding to a ``Permutation``
          you will be glad to learn about the ``Permutation.to_matrix()`` method.

        - The base ring of the matrix can be anything that can be coerced to ``RR``.

    .. SEEALSO::

        - :meth:`~sage.matrix.matrix2.as_sum_of_permutations`
          to use this method through the ``Matrix`` class.

    EXAMPLES:

    We create a bistochastic matrix from a convex sum of permutations, then
    try to deduce the decomposition from the matrix::

        sage: from sage.combinat.permutation import bistochastic_as_sum_of_permutations
        sage: L = []
        sage: L.append((9,Permutation([4, 1, 3, 5, 2])))
        sage: L.append((6,Permutation([5, 3, 4, 1, 2])))
        sage: L.append((3,Permutation([3, 1, 4, 2, 5])))
        sage: L.append((2,Permutation([1, 4, 2, 3, 5])))
        sage: M = sum([c * p.to_matrix() for (c,p) in L])
        sage: decomp = bistochastic_as_sum_of_permutations(M)
        sage: print(decomp)
        2*B[[1, 4, 2, 3, 5]] + 3*B[[3, 1, 4, 2, 5]] + 9*B[[4, 1, 3, 5, 2]] + 6*B[[5, 3, 4, 1, 2]]

    An exception is raised when the matrix is not positive and bistochastic::

        sage: M = Matrix([[2,3],[2,2]])
        sage: decomp = bistochastic_as_sum_of_permutations(M)
        Traceback (most recent call last):
        ...
        ValueError: The matrix is not bistochastic

        sage: bistochastic_as_sum_of_permutations(Matrix(GF(7), 2, [2,1,1,2]))
        Traceback (most recent call last):
        ...
        ValueError: The base ring of the matrix must have a coercion map to RR

        sage: bistochastic_as_sum_of_permutations(Matrix(ZZ, 2, [2,-1,-1,2]))
        Traceback (most recent call last):
        ...
        ValueError: The matrix should have nonnegative entries
    """
    from sage.graphs.bipartite_graph import BipartiteGraph
    from sage.combinat.free_module import CombinatorialFreeModule
    from sage.rings.all import RR

    n = M.nrows()

    if n != M.ncols():
        raise ValueError("The matrix is expected to be square")

    if not all(x >= 0 for x in M.list()):
        raise ValueError("The matrix should have nonnegative entries")

    if check and not M.is_bistochastic(normalized = False):
        raise ValueError("The matrix is not bistochastic")

    if not RR.has_coerce_map_from(M.base_ring()):
        raise ValueError("The base ring of the matrix must have a coercion map to RR")

    CFM = CombinatorialFreeModule(M.base_ring(), Permutations(n))
    value = 0

    G = BipartiteGraph(M, weighted=True)
    P = Permutations()

    while G.size() > 0:
        matching = G.matching(use_edge_labels=True)
        matching = [(min(u,v), max(u, v), w) for u,v,w in matching]

        # This minimum is strictly larger than 0
        minimum = min([x[2] for x in matching])

        for (u,v,l) in matching:
            if minimum == l:
                G.delete_edge((u,v,l))
            else:
                G.set_edge_label(u,v,l-minimum)

        matching.sort(key=lambda x: x[0])
        value += minimum * CFM(P([x[1]-n+1 for x in matching]))

    return value


def bounded_affine_permutation(A):
    r"""
    Return the bounded affine permutation of a matrix.

    The *bounded affine permutation* of a matrix `A` with entries in `R`
    is a partial permutation of length `n`, where `n` is the number of
    columns of `A`. The entry in position `i` is the smallest value `j`
    such that column `i` is in the span of columns `i+1, \ldots, j`,
    over `R`, where column indices are taken modulo `n`.
    If column `i` is the zero vector, then the permutation has a
    fixed point at `i`.

    INPUT:

    - ``A`` -- matrix with entries in a ring `R`

    EXAMPLES::

        sage: from sage.combinat.permutation import bounded_affine_permutation
        sage: A = Matrix(ZZ, [[1,0,0,0], [0,1,0,0]])
        sage: bounded_affine_permutation(A)
        [5, 6, 3, 4]

        sage: A = Matrix(ZZ, [[0,1,0,1,0], [0,0,1,1,0]])
        sage: bounded_affine_permutation(A)
        [1, 4, 7, 8, 5]

    REFERENCES:

    - [KLS2013]_
    """
    n = A.ncols()
    R = A.base_ring()
    from sage.modules.free_module import FreeModule
    from sage.modules.free_module import span
    z = FreeModule(R, A.nrows()).zero()
    v = A.columns()
    perm = []
    for j in range(n):
        if not v[j]:
            perm.append(j + 1)
            continue
        V = span([z], R)
        for i in range(j + 1, j + n + 1):
            index = i % n
            V = V + span([v[index]], R)
            if not V.dimension():
                continue
            if v[j] in V:
                perm.append(i + 1)
                break
    S = Permutations(2 * n, n)
    return S(perm)


class StandardPermutations_descents(StandardPermutations_n_abstract):
    r"""
    Permutations of `\{1, \ldots, n\}` with a fixed set of descents.
    """
    @staticmethod
    def __classcall_private__(cls, d, n):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = Permutations(descents=([1,0,4,8],12))
            sage: P2 = Permutations(descents=((0,1,4,8),12))
            sage: P1 is P2
            True
        """
        return super(StandardPermutations_descents, cls).__classcall__(cls, tuple(sorted(d)), n)

    def __init__(self, d, n):
        r"""
        The class of all permutations of `\{1, 2, ..., n\}`
        with set of descent positions `d` (where the descent positions
        are being counted from `0`, so that `i` lies in this set if
        and only if the permutation takes a larger value at `i + 1` than
        at `i + 2`).

        TESTS::

            sage: P = Permutations(descents=([1,0,2], 5))
            sage: TestSuite(P).run()
        """
        StandardPermutations_n_abstract.__init__(self, n)
        self._d = d

    def _repr_(self):
        """
        TESTS::

            sage: Permutations(descents=([1,0,4,8],12))
            Standard permutations of 12 with descents [0, 1, 4, 8]
        """
        return "Standard permutations of %s with descents %s"%(self.n, list(self._d))

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        ALGORITHM:

        The algorithm described in [Vie1979]_ is implemented naively.

        EXAMPLES::

            sage: P = Permutations(descents=([1,0,2], 5))
            sage: P.cardinality()
            4

        TESTS::

            sage: Permutations(descents=([], 1)).cardinality()
            1

            sage: Permutations(descents=([1,4], 6)).cardinality()
            40

            sage: P = lambda D, n: Permutations(descents=(D, n+1))
            sage: all(P(D, n).cardinality() == len(P(D, n).list())
            ....:     for n in range(5) for D in subsets(range(n)))
            True

            sage: n = 20;
            sage: D = [6, 8, 10, 11, 12, 13, 14, 15, 17, 19]
            sage: P(D, n).cardinality()
            125291047596

        """
        def m(l):
            s = 0
            partial_sums = [0]
            for i in l:
                s += i
                partial_sums.append(s)
            return partial_sums

        def d(l):
            return m(reversed(l))[::-1]

        one = ZZ.one()
        if not self._d:
            return one

        l_ops = [1]*(self.n-1)
        for i in self._d:
            l_ops[i] = 0
        l = [one]
        for op in reversed(l_ops):
            if op:
                l = m(l)
            else:
                l = d(l)

        return sum(l)

    def first(self):
        """
        Return the first permutation with descents `d`.

        EXAMPLES::

            sage: Permutations(descents=([1,0,4,8],12)).first()
            [3, 2, 1, 4, 6, 5, 7, 8, 10, 9, 11, 12]
        """
        return descents_composition_first(Composition(descents=(self._d, self.n)))

    def last(self):
        """
        Return the last permutation with descents `d`.

        EXAMPLES::

            sage: Permutations(descents=([1,0,4,8],12)).last()
            [12, 11, 8, 9, 10, 4, 5, 6, 7, 1, 2, 3]
        """
        return descents_composition_last(Composition(descents=(self._d, self.n)))

    def __iter__(self):
        """
        Iterate over all the permutations that have the descents `d`.

        EXAMPLES::

            sage: Permutations(descents=([2,0],5)).list()
            [[5, 2, 4, 1, 3],
             [5, 3, 4, 1, 2],
             [4, 3, 5, 1, 2],
             [4, 2, 5, 1, 3],
             [3, 2, 5, 1, 4],
             [2, 1, 5, 3, 4],
             [3, 1, 5, 2, 4],
             [4, 1, 5, 2, 3],
             [5, 1, 4, 2, 3],
             [5, 1, 3, 2, 4],
             [4, 1, 3, 2, 5],
             [3, 1, 4, 2, 5],
             [2, 1, 4, 3, 5],
             [3, 2, 4, 1, 5],
             [4, 2, 3, 1, 5],
             [5, 2, 3, 1, 4]]
        """
        return iter(descents_composition_list(Composition(descents=(self._d, self.n))))

def descents_composition_list(dc):
    """
    Return a list of all the permutations that have the descent
    composition ``dc``.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.descents_composition_list([1,2,2])
        [[5, 2, 4, 1, 3],
         [5, 3, 4, 1, 2],
         [4, 3, 5, 1, 2],
         [4, 2, 5, 1, 3],
         [3, 2, 5, 1, 4],
         [2, 1, 5, 3, 4],
         [3, 1, 5, 2, 4],
         [4, 1, 5, 2, 3],
         [5, 1, 4, 2, 3],
         [5, 1, 3, 2, 4],
         [4, 1, 3, 2, 5],
         [3, 1, 4, 2, 5],
         [2, 1, 4, 3, 5],
         [3, 2, 4, 1, 5],
         [4, 2, 3, 1, 5],
         [5, 2, 3, 1, 4]]
    """
    return [p.inverse() for p in StandardPermutations_recoils(dc)]

def descents_composition_first(dc):
    r"""
    Compute the smallest element of a descent class having a descent
    composition ``dc``.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.descents_composition_first([1,1,3,4,3])
        [3, 2, 1, 4, 6, 5, 7, 8, 10, 9, 11, 12]
    """
    if not isinstance(dc, Composition):
        try:
            dc = Composition(dc)
        except TypeError:
            raise TypeError("The argument must be of type Composition")

    cpl = [x for x in reversed(dc.conjugate())]
    res = []
    s = 0
    for i in range(len(cpl)):
        res += [s + cpl[i]-j for j in range(cpl[i])]
        s   += cpl[i]

    return Permutations()(res)

def descents_composition_last(dc):
    r"""
    Return the largest element of a descent class having a descent
    composition ``dc``.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.descents_composition_last([1,1,3,4,3])
        [12, 11, 8, 9, 10, 4, 5, 6, 7, 1, 2, 3]
    """
    if not isinstance(dc, Composition):
        try:
            dc = Composition(dc)
        except TypeError:
            raise TypeError("The argument must be of type Composition")
    s = 0
    res = []
    for i in reversed(range(len(dc))):
        res = [j for j in range(s+1,s+dc[i]+1)] + res
        s += dc[i]

    return Permutations()(res)

class StandardPermutations_recoilsfiner(Permutations):
    @staticmethod
    def __classcall_private__(cls, recoils):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(recoils_finer=[2,2])
            sage: S2 = Permutations(recoils_finer=(2,2))
            sage: S1 is S2
            True
        """
        return super(StandardPermutations_recoilsfiner, cls).__classcall__(cls, Composition(recoils))

    def __init__(self, recoils):
        """
        TESTS::

            sage: P = Permutations(recoils_finer=[2,2])
            sage: TestSuite(P).run()
        """
        Permutations.__init__(self, category=FiniteEnumeratedSets())
        self.recoils = recoils

    def _repr_(self):
        """
        TESTS::

            sage: Permutations(recoils_finer=[2,2])
            Standard permutations whose recoils composition is finer than [2, 2]
        """
        return "Standard permutations whose recoils composition is finer than %s"%self.recoils

    def __iter__(self):
        """
        Iterate over of all of the permutations whose recoils composition
        is finer than ``self.recoils``.

        EXAMPLES::

            sage: Permutations(recoils_finer=[2,2]).list()
            [[3, 1, 4, 2],
             [3, 4, 1, 2],
             [1, 3, 4, 2],
             [1, 3, 2, 4],
             [1, 2, 3, 4],
             [3, 1, 2, 4]]
        """
        recoils = self.recoils
        dag = DiGraph()

        #Add the nodes
        for i in range(1, sum(recoils)+1):
            dag.add_vertex(i)

        #Add the edges to guarantee a finer recoil composition
        pos = 1
        for part in recoils:
            for i in range(part-1):
                dag.add_edge(pos, pos+1)
                pos += 1
            pos += 1

        for le in dag.topological_sort_generator():
            yield self.element_class(self, le)

class StandardPermutations_recoilsfatter(Permutations):
    @staticmethod
    def __classcall_private__(cls, recoils):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(recoils_fatter=[2,2])
            sage: S2 = Permutations(recoils_fatter=(2,2))
            sage: S1 is S2
            True
        """
        return super(StandardPermutations_recoilsfatter, cls).__classcall__(cls, Composition(recoils))

    def __init__(self, recoils):
        """
        TESTS::

            sage: P = Permutations(recoils_fatter=[2,2])
            sage: TestSuite(P).run()
        """
        Permutations.__init__(self, category=FiniteEnumeratedSets())
        self.recoils = recoils

    def _repr_(self):
        """
        TESTS::

            sage: Permutations(recoils_fatter=[2,2])
            Standard permutations whose recoils composition is fatter than [2, 2]
        """
        return "Standard permutations whose recoils composition is fatter than %s"%self.recoils

    def __iter__(self):
        """
        Iterate over of all of the permutations whose recoils composition
        is fatter than ``self.recoils``.

        EXAMPLES::

            sage: Permutations(recoils_fatter=[2,2]).list()
            [[4, 3, 2, 1],
             [3, 2, 1, 4],
             [3, 2, 4, 1],
             [3, 4, 2, 1],
             [3, 4, 1, 2],
             [3, 1, 4, 2],
             [1, 3, 4, 2],
             [1, 3, 2, 4],
             [3, 1, 2, 4],
             [1, 4, 3, 2],
             [4, 1, 3, 2],
             [4, 3, 1, 2]]
        """
        recoils = self.recoils
        dag = DiGraph()

        #Add the nodes
        for i in range(1, sum(recoils)+1):
            dag.add_vertex(i)

        #Add the edges to guarantee a fatter recoil composition
        pos = 0
        for i in range(len(recoils)-1):
            pos += recoils[i]
            dag.add_edge(pos+1, pos)

        for le in dag.topological_sort_generator():
            yield self.element_class(self, le)

class StandardPermutations_recoils(Permutations):
    r"""
    Permutations of `\{1, \ldots, n\}` with a fixed recoils composition.
    """
    @staticmethod
    def __classcall_private__(cls, recoils):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(recoils=[2,2])
            sage: S2 = Permutations(recoils=(2,2))
            sage: S1 is S2
            True
        """
        return super(StandardPermutations_recoils, cls).__classcall__(cls, Composition(recoils))

    def __init__(self, recoils):
        """
        TESTS::

            sage: P = Permutations(recoils=[2,2])
            sage: TestSuite(P).run()
        """
        Permutations.__init__(self, category=FiniteEnumeratedSets())
        self.recoils = recoils

    def _repr_(self):
        """
        TESTS::

            sage: Permutations(recoils=[2,2])
            Standard permutations whose recoils composition is [2, 2]
        """
        return "Standard permutations whose recoils composition is %s"%self.recoils

    def __iter__(self):
        """
        Iterate over of all of the permutations whose recoils composition
        is equal to ``self.recoils``.

        EXAMPLES::

            sage: Permutations(recoils=[2,2]).list()
            [[3, 1, 4, 2], [3, 4, 1, 2], [1, 3, 4, 2], [1, 3, 2, 4], [3, 1, 2, 4]]
        """
        recoils = self.recoils
        dag = DiGraph()

        #Add all the nodes
        for i in range(1, sum(recoils)+1):
            dag.add_vertex(i)

        #Add the edges which guarantee a finer recoil comp.
        pos = 1
        for part in recoils:
            for i in range(part-1):
                dag.add_edge(pos, pos+1)
                pos += 1
            pos += 1

        #Add the edges which guarantee a fatter recoil comp.
        pos = 0
        for i in range(len(recoils)-1):
            pos += recoils[i]
            dag.add_edge(pos+1, pos)

        for le in dag.topological_sort_generator():
            yield self.element_class(self, le)

def from_major_code(mc, final_descent=False):
    r"""
    Return the permutation with major code ``mc``.

    The major code of a permutation is defined in
    :meth:`~sage.combinat.permutation.Permutation.to_major_code`.

    .. WARNING::

        This function creates illegal permutations (i.e. ``Permutation([9])``,
        and this is dangerous as the :meth:`Permutation` class is only designed
        to handle permutations on `1...n`. This will have to be changed when Sage
        permutations will be able to handle anything, but right now this should
        be fixed. Be careful with the results.

    .. WARNING::

        If ``mc`` is not a major index of a permutation, then the return
        value of this method can be anything. Garbage in, garbage out!

    REFERENCES:

    - Skandera, M. *An Eulerian Partner for Inversions*. Sem.
      Lothar. Combin. 46 (2001) B46d.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_major_code([5, 0, 1, 0, 1, 2, 0, 1, 0])
        [9, 3, 5, 7, 2, 1, 4, 6, 8]
        sage: permutation.from_major_code([8, 3, 3, 1, 4, 0, 1, 0, 0])
        [2, 8, 4, 3, 6, 7, 9, 5, 1]
        sage: Permutation([2,1,6,4,7,3,5]).to_major_code()
        [3, 2, 0, 2, 2, 0, 0]
        sage: permutation.from_major_code([3, 2, 0, 2, 2, 0, 0])
        [2, 1, 6, 4, 7, 3, 5]

    TESTS::

        sage: permutation.from_major_code([])
        []

        sage: all( permutation.from_major_code(p.to_major_code()) == p
        ....:      for p in Permutations(5) )
        True
    """
    if not mc:
        w = []
    else:
        #define w^(n) to be the one-letter word n
        w = [len(mc)]

    #for i=n-1,..,1 let w^i be the unique word obtained by inserting
    #the letter i into the word w^(i+1) in such a way that
    #maj(w^i)-maj(w^(i+1)) = mc[i]

    for i in reversed(range(1,len(mc))):
        #Lemma 2.2 in Skandera

        #Get the descents of w and place them in reverse order
        d = Permutation(w, check_input=False).descents(final_descent=final_descent)
        d.reverse()

        #a is the list of all positions which are not descents
        a = [x for x in range(1, len(w) + 1) if x not in d]

        #d_k = -1    -- 0 in the lemma, but -1 due to 0-based indexing
        d.append(0)
        l = mc[i-1]
        indices = d + a
        w.insert(indices[l], i)

    return Permutation(w, check_input = False)

################
# Bruhat Order #
################


class StandardPermutations_bruhat_smaller(Permutations):
    r"""
    Permutations of `\{1, \ldots, n\}` that are less than or equal to a
    permutation `p` in the Bruhat order.
    """
    @staticmethod
    def __classcall_private__(cls, p):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(bruhat_smaller=[2,3,1])
            sage: S2 = Permutations(bruhat_smaller=(1,2,3))
            sage: S1 is S2
            True
        """
        return super(StandardPermutations_bruhat_smaller, cls).__classcall__(cls, Permutation(p))

    def __init__(self, p):
        """
        TESTS::

            sage: P = Permutations(bruhat_smaller=[3,2,1])
            sage: TestSuite(P).run()
        """
        Permutations.__init__(self, category=FiniteEnumeratedSets())
        self.p = p

    def _repr_(self):
        """
        TESTS::

            sage: Permutations(bruhat_smaller=[3,2,1])
            Standard permutations that are less than or equal to [3, 2, 1] in the Bruhat order
        """
        return "Standard permutations that are less than or equal to %s in the Bruhat order"%self.p

    def __iter__(self):
        r"""
        Iterate through a list of permutations smaller than or equal to ``p``
        in the Bruhat order.

        EXAMPLES::

            sage: Permutations(bruhat_smaller=[4,1,2,3]).list()
            [[1, 2, 3, 4],
             [1, 2, 4, 3],
             [1, 3, 2, 4],
             [1, 4, 2, 3],
             [2, 1, 3, 4],
             [2, 1, 4, 3],
             [3, 1, 2, 4],
             [4, 1, 2, 3]]
        """
        return iter(transitive_ideal(lambda x: x.bruhat_pred(), self.p))


class StandardPermutations_bruhat_greater(Permutations):
    r"""
    Permutations of `\{1, \ldots, n\}` that are greater than or equal to a
    permutation `p` in the Bruhat order.
    """
    @staticmethod
    def __classcall_private__(cls, p):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: S1 = Permutations(bruhat_greater=[2,3,1])
            sage: S2 = Permutations(bruhat_greater=(1,2,3))
            sage: S1 is S2
            True
        """
        return super(StandardPermutations_bruhat_greater, cls).__classcall__(cls, Permutation(p))

    def __init__(self, p):
        """
        TESTS::

            sage: P = Permutations(bruhat_greater=[3,2,1])
            sage: TestSuite(P).run()
        """
        Permutations.__init__(self, category=FiniteEnumeratedSets())
        self.p = p

    def _repr_(self):
        """
        TESTS::

            sage: Permutations(bruhat_greater=[3,2,1])
            Standard permutations that are greater than or equal to [3, 2, 1] in the Bruhat order
        """
        return "Standard permutations that are greater than or equal to %s in the Bruhat order"%self.p

    def __iter__(self):
        r"""
        Iterate through a list of permutations greater than or equal to ``p``
        in the Bruhat order.

        EXAMPLES::

            sage: Permutations(bruhat_greater=[4,1,2,3]).list()
            [[4, 1, 2, 3],
             [4, 1, 3, 2],
             [4, 2, 1, 3],
             [4, 2, 3, 1],
             [4, 3, 1, 2],
             [4, 3, 2, 1]]
        """
        return iter(transitive_ideal(lambda x: x.bruhat_succ(), self.p))

def bruhat_lequal(p1, p2):
    r"""
    Return ``True`` if ``p1`` is less than ``p2`` in the Bruhat order.

    Algorithm from mupad-combinat.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.bruhat_lequal([2,4,3,1],[3,4,2,1])
        True
    """
    n1 = len(p1)

    if n1 == 0:
        return True

    if p1[0] > p2[0] or p1[n1-1] < p2[n1-1]:
        return False

    for i in range(n1):
        c = 0
        for j in range(n1):
            if p2[j] > i+1:
                c += 1
            if p1[j] > i+1:
                c -= 1
            if c < 0:
                return False

    return True


#################
# Permutohedron #
#################

def permutohedron_lequal(p1, p2, side="right"):
    r"""
    Return ``True`` if ``p1`` is less than or equal to ``p2`` in the
    permutohedron order.

    By default, the computations are done in the right permutohedron.
    If you pass the option ``side='left'``, then they will be done in the
    left permutohedron.

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.permutohedron_lequal(Permutation([3,2,1,4]),Permutation([4,2,1,3]))
        False
        sage: permutation.permutohedron_lequal(Permutation([3,2,1,4]),Permutation([4,2,1,3]), side='left')
        True
    """
    l1 = p1.number_of_inversions()
    l2 = p2.number_of_inversions()

    if l1 > l2:
        return False

    if side == "right":
        prod = p1._left_to_right_multiply_on_right(p2.inverse())
    else:
        prod = p1._left_to_right_multiply_on_left(p2.inverse())

    return prod.number_of_inversions() == l2 - l1


############
# Patterns #
############
from sage.combinat.words.finite_word import evaluation_dict

def to_standard(p, key=None):
    r"""
    Return a standard permutation corresponding to the iterable ``p``.

    INPUT:

    - ``p`` -- an iterable
    - ``key`` -- (optional) a comparison key for the element
      ``x`` of ``p``

    EXAMPLES::

        sage: import sage.combinat.permutation as permutation
        sage: permutation.to_standard([4,2,7])
        [2, 1, 3]
        sage: permutation.to_standard([1,2,3])
        [1, 2, 3]
        sage: permutation.to_standard([])
        []
        sage: permutation.to_standard([1,2,3], key=lambda x: -x)
        [3, 2, 1]
        sage: permutation.to_standard([5,8,2,5], key=lambda x: -x)
        [2, 1, 4, 3]

    TESTS:

    Does not mutate the list::

        sage: a = [1,2,4]
        sage: permutation.to_standard(a)
        [1, 2, 3]
        sage: a
        [1, 2, 4]

    We check against the naive method::

        sage: def std(p):
        ....:     s = [0]*len(p)
        ....:     c = p[:]
        ....:     biggest = max(p) + 1
        ....:     i = 1
        ....:     for _ in range(len(c)):
        ....:         smallest = min(c)
        ....:         smallest_index = c.index(smallest)
        ....:         s[smallest_index] = i
        ....:         i += 1
        ....:         c[smallest_index] = biggest
        ....:     return Permutations()(s)
        sage: p = list(Words(100, 1000).random_element())
        sage: std(p) == permutation.to_standard(p)
        True
    """
    ev_dict = evaluation_dict(p)
    ordered_alphabet = sorted(ev_dict, key=key)
    offset = 0
    for k in ordered_alphabet:
        temp = ev_dict[k]
        ev_dict[k] = offset
        offset += temp
    result = []
    for l in p:
        ev_dict[l] += 1
        result.append(ev_dict[l])
    return Permutations(len(result))(result)


##########################################################

class CyclicPermutations(Permutations_mset):
    """
    Return the class of all cyclic permutations of ``mset`` in cycle notation.
    These are the same as necklaces.

    INPUT:

    - ``mset`` -- A multiset

    EXAMPLES::

        sage: CyclicPermutations(range(4)).list()
        [[0, 1, 2, 3],
         [0, 1, 3, 2],
         [0, 2, 1, 3],
         [0, 2, 3, 1],
         [0, 3, 1, 2],
         [0, 3, 2, 1]]
        sage: CyclicPermutations([1,1,1]).list()
        [[1, 1, 1]]
    """
    @staticmethod
    def __classcall_private__(cls, mset):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: CP1 = CyclicPermutations([1,1,1])
            sage: CP2 = CyclicPermutations((1,1,1))
            sage: CP1 is CP2
            True

            sage: CP = CyclicPermutations([1,2,3,3])
            sage: CP
            Cyclic permutations of [1, 2, 3, 3]
            sage: TestSuite(CP).run() # not tested -- broken
        """
        return super(CyclicPermutations, cls).__classcall__(cls, tuple(mset))

    def _repr_(self):
        """
        TESTS::

            sage: CyclicPermutations(range(4))
            Cyclic permutations of [0, 1, 2, 3]
        """
        return "Cyclic permutations of %s"%list(self.mset)

    def __iter__(self, distinct=False):
        """
        EXAMPLES::

            sage: CyclicPermutations(range(4)).list() # indirect doctest
            [[0, 1, 2, 3],
             [0, 1, 3, 2],
             [0, 2, 1, 3],
             [0, 2, 3, 1],
             [0, 3, 1, 2],
             [0, 3, 2, 1]]
             sage: CyclicPermutations([1,1,1]).list()
             [[1, 1, 1]]
             sage: CyclicPermutations([1,1,1]).list(distinct=True)
             [[1, 1, 1], [1, 1, 1]]
        """
        if distinct:
            content = [1] * len(self.mset)
        else:
            content = [0] * len(self.mset)
            index_list = map(self.mset.index, self.mset)
            for i in index_list:
                content[i] += 1

        from .necklace import Necklaces
        for necklace in Necklaces(content):
            yield [self.mset[x-1] for x in necklace]

    iterator = __iter__

    def list(self, distinct=False):
        """
        EXAMPLES::

            sage: CyclicPermutations(range(4)).list()
            [[0, 1, 2, 3],
             [0, 1, 3, 2],
             [0, 2, 1, 3],
             [0, 2, 3, 1],
             [0, 3, 1, 2],
             [0, 3, 2, 1]]
        """
        return list(self.__iter__(distinct=distinct))

#################################################

class CyclicPermutationsOfPartition(Permutations):
    """
    Combinations of cyclic permutations of each cell of a given partition.

    This is the same as a Cartesian product of necklaces.

    EXAMPLES::

        sage: CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]]).list()
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]

    ::

        sage: CyclicPermutationsOfPartition([[1,2,3,4],[4,4,4]]).list()
        [[[1, 2, 3, 4], [4, 4, 4]],
         [[1, 2, 4, 3], [4, 4, 4]],
         [[1, 3, 2, 4], [4, 4, 4]],
         [[1, 3, 4, 2], [4, 4, 4]],
         [[1, 4, 2, 3], [4, 4, 4]],
         [[1, 4, 3, 2], [4, 4, 4]]]

    ::

        sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list()
        [[[1, 2, 3], [4, 4, 4]], [[1, 3, 2], [4, 4, 4]]]

    ::

        sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list(distinct=True)
        [[[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]],
         [[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]]]
    """
    @staticmethod
    def __classcall_private__(cls, partition):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: CP1 = CyclicPermutationsOfPartition([[1,2,3],[4,4,4]])
            sage: CP2 = CyclicPermutationsOfPartition([[1,2,3],[4,4,4]])
            sage: CP1 is CP2
            True
        """
        partition = tuple(map(tuple, partition))
        return super(CyclicPermutationsOfPartition, cls).__classcall__(cls, partition)

    def __init__(self, partition):
        """
        TESTS::

            sage: CP = CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]])
            sage: CP
            Cyclic permutations of partition [[1, 2, 3, 4], [5, 6, 7]]
            sage: TestSuite(CP).run()
        """
        self.partition = partition
        Permutations.__init__(self, category=FiniteEnumeratedSets())

    class Element(ClonableArray):
        """
        A cyclic permutation of a partition.
        """
        def check(self):
            """
            Check that ``self`` is a valid element.

            EXAMPLES::

                sage: CP = CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]])
                sage: elt = CP[0]
                sage: elt.check()
            """
            if [sorted(_) for _ in self] != [sorted(_) for _ in self.parent().partition]:
                raise ValueError("Invalid cyclic permutation of the partition"%self.parent().partition)

    def _repr_(self):
        """
        TESTS::

            sage: CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]])
            Cyclic permutations of partition [[1, 2, 3, 4], [5, 6, 7]]
        """
        return "Cyclic permutations of partition {}".format(
            [list(_) for _ in self.partition])

    def __iter__(self, distinct=False):
        """
        AUTHORS:

        - Robert Miller

        EXAMPLES::

            sage: CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]]).list() # indirect doctest
            [[[1, 2, 3, 4], [5, 6, 7]],
             [[1, 2, 4, 3], [5, 6, 7]],
             [[1, 3, 2, 4], [5, 6, 7]],
             [[1, 3, 4, 2], [5, 6, 7]],
             [[1, 4, 2, 3], [5, 6, 7]],
             [[1, 4, 3, 2], [5, 6, 7]],
             [[1, 2, 3, 4], [5, 7, 6]],
             [[1, 2, 4, 3], [5, 7, 6]],
             [[1, 3, 2, 4], [5, 7, 6]],
             [[1, 3, 4, 2], [5, 7, 6]],
             [[1, 4, 2, 3], [5, 7, 6]],
             [[1, 4, 3, 2], [5, 7, 6]]]

        ::

            sage: CyclicPermutationsOfPartition([[1,2,3,4],[4,4,4]]).list()
            [[[1, 2, 3, 4], [4, 4, 4]],
             [[1, 2, 4, 3], [4, 4, 4]],
             [[1, 3, 2, 4], [4, 4, 4]],
             [[1, 3, 4, 2], [4, 4, 4]],
             [[1, 4, 2, 3], [4, 4, 4]],
             [[1, 4, 3, 2], [4, 4, 4]]]

        ::

            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list()
            [[[1, 2, 3], [4, 4, 4]], [[1, 3, 2], [4, 4, 4]]]

        ::

            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list(distinct=True)
            [[[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]],
             [[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]]]
        """
        if len(self.partition) == 1:
            for i in CyclicPermutations(self.partition[0]).iterator(distinct=distinct):
                yield self.element_class(self, [i])
        else:
            for right in CyclicPermutationsOfPartition(self.partition[1:]).iterator(distinct=distinct):
                for perm in CyclicPermutations(self.partition[0]).iterator(distinct=distinct):
                    yield self.element_class(self, [perm] + list(right))

    iterator = __iter__

    def list(self, distinct=False):
        """
        EXAMPLES::

            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list()
            [[[1, 2, 3], [4, 4, 4]], [[1, 3, 2], [4, 4, 4]]]
            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list(distinct=True)
            [[[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]],
             [[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]]]
        """
        return list(self.iterator(distinct=distinct))


###############################################
#Avoiding

class StandardPermutations_all_avoiding(StandardPermutations_all):
    """
    All standard permutations avoiding a set of patterns.
    """
    @staticmethod
    def __classcall_private__(cls, a):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: P1 = Permutations(avoiding=([2,1,3],[1,2,3]))
            sage: P2 = Permutations(avoiding=[[2,1,3],[1,2,3]])
            sage: P1 is P2
            True
        """
        a = tuple(map(Permutation, a))
        return super(StandardPermutations_all_avoiding, cls).__classcall__(cls, a)

    def __init__(self, a):
        """
        TESTS::

            sage: P = Permutations(avoiding=[[2,1,3],[1,2,3]])
            sage: TestSuite(P).run(max_runs=25)
        """
        Permutations.__init__(self, category=InfiniteEnumeratedSets())
        self._a = a

    def patterns(self):
        """
        Return the patterns avoided by this class of permutations.

        EXAMPLES::

            sage: P = Permutations(avoiding=[[2,1,3],[1,2,3]])
            sage: P.patterns()
            ([2, 1, 3], [1, 2, 3])
        """
        return self._a

    def _repr_(self):
        """
        EXAMPLES::

            sage: Permutations(avoiding=[[2,1,3],[1,2,3]])
            Standard permutations avoiding [[2, 1, 3], [1, 2, 3]]
        """
        return "Standard permutations avoiding %s"%(list(self._a))

    def __contains__(self, x):
        """
        TESTS::

            sage: [1,3,2] in Permutations(avoiding=[1,3,2])
            False
            sage: [1,3,2] in Permutations(avoiding=[[1,3,2]])
            False
            sage: [2,1,3] in Permutations(avoiding=[[1,3,2],[1,2,3]])
            True
            sage: [2,1,3] in Permutations(avoiding=[])
            True
        """
        if not super(StandardPermutations_all_avoiding, self).__contains__(x):
            return False
        x = Permutations()(x)
        return all(x.avoids(p) for p in self._a)

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: it = iter(Permutations(avoiding=[[2,1,3],[1,2,3]]))
            sage: [next(it) for i in range(10)]
            [[],
             [1],
             [1, 2],
             [2, 1],
             [1, 3, 2],
             [2, 3, 1],
             [3, 1, 2],
             [3, 2, 1],
             [1, 4, 3, 2],
             [2, 4, 3, 1]]
        """
        n = 0
        while True:
            for x in itertools.permutations(range(1, n + 1)):
                x = self.element_class(self, x)
                if all(x.avoids(p) for p in self._a):
                    yield x
            n += 1

class StandardPermutations_avoiding_generic(StandardPermutations_n_abstract):
    """
    Generic class for subset of permutations avoiding a set of patterns.
    """
    @staticmethod
    def __classcall_private__(cls, n, a):
        """
        Normalize arguments to ensure a unique representation.

        TESTS::

            sage: P1 = Permutations(3, avoiding=([2,1,3],[1,2,3]))
            sage: P2 = Permutations(3, avoiding=[[2,1,3],[1,2,3]])
            sage: P1 is P2
            True
        """
        a = tuple(map(Permutation, a))
        return super(StandardPermutations_avoiding_generic, cls).__classcall__(cls, n, a)

    def __init__(self, n, a):
        """
        EXAMPLES::

            sage: P = Permutations(3, avoiding=[[2,1,3],[1,2,3]])
            sage: TestSuite(P).run()
            sage: type(P)
            <class 'sage.combinat.permutation.StandardPermutations_avoiding_generic_with_category'>
        """
        StandardPermutations_n_abstract.__init__(self, n)
        self._a = a

    @property
    def a(self):
        r"""
        ``self.a`` is deprecated; use :meth:`patterns` instead.

        TESTS::

            sage: P = Permutations(3, avoiding=[[2,1,3],[1,2,3]])
            sage: P.a
            doctest:...: DeprecationWarning: The attribute a for the list of patterns to avoid is deprecated, use the method patterns instead.
            See https://trac.sagemath.org/26810 for details.
            ([2, 1, 3], [1, 2, 3])
        """
        from sage.misc.superseded import deprecation
        deprecation(26810, "The attribute a for the list of patterns to avoid is "
                           "deprecated, use the method patterns instead.")
        return self.patterns()

    def patterns(self):
        """
        Return the patterns avoided by this class of permutations.

        EXAMPLES::

            sage: P = Permutations(3, avoiding=[[2,1,3],[1,2,3]])
            sage: P.patterns()
            ([2, 1, 3], [1, 2, 3])
        """
        return self._a

    def __contains__(self, x):
        """
        TESTS::

            sage: [1,3,2] in Permutations(3, avoiding=[1,3,2])
            False
            sage: [1,3,2] in Permutations(3, avoiding=[[1,3,2]])
            False
            sage: [2,1,3] in Permutations(3, avoiding=[[1,3,2],[1,2,3]])
            True
            sage: [2,1,3] in Permutations(3, avoiding=[])
            True
        """
        if not super(StandardPermutations_avoiding_generic, self).__contains__(x):
            return False
        x = Permutations()(x)
        return all(x.avoids(p) for p in self._a)

    def _repr_(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[[2, 1, 3],[1,2,3]])
            Standard permutations of 3 avoiding [[2, 1, 3], [1, 2, 3]]
        """
        return "Standard permutations of %s avoiding %s"%(self.n, list(self._a))

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[[2, 1, 3],[1,2,3]]).list()
            [[1, 3, 2], [3, 1, 2], [2, 3, 1], [3, 2, 1]]
            sage: Permutations(0, avoiding=[[2, 1, 3],[1,2,3]]).list()
            [[]]
        """
        if self.n > 0:
            return iter(PatternAvoider(self, self._a))
        return iter([self.element_class(self, [])])

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: P = Permutations(3, avoiding=[[2, 1, 3],[1,2,3]])
            sage: P.cardinality()
            4
        """
        one = ZZ.one()
        return sum(one for p in self)

class StandardPermutations_avoiding_12(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[1, 2])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_12, self).__init__(n, (Permutations()([1, 2]),))

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[1,2]).list()
            [[3, 2, 1]]
        """
        yield self.element_class(self, range(self.n, 0, -1))

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: P = Permutations(3, avoiding=[1, 2])
            sage: P.cardinality()
            1
        """
        return ZZ.one()

class StandardPermutations_avoiding_21(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[2, 1])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_21, self).__init__(n, (Permutations()([2, 1]),))

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[2,1]).list()
            [[1, 2, 3]]
        """
        yield self.element_class(self, range(1, self.n+1))

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: P = Permutations(3, avoiding=[2, 1])
            sage: P.cardinality()
            1
        """
        return ZZ.one()

class StandardPermutations_avoiding_132(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[1, 3, 2])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_132, self).__init__(n, (Permutations()([1, 3, 2]),))

    def cardinality(self):
        """
        EXAMPLES::

            sage: Permutations(5, avoiding=[1, 3, 2]).cardinality()
            42
            sage: len( Permutations(5, avoiding=[1, 3, 2]).list() )
            42
        """
        return catalan_number(self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[1,3,2]).list() # indirect doctest
            [[1, 2, 3], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
            sage: Permutations(4, avoiding=[1,3,2]).list()
            [[4, 1, 2, 3],
             [4, 2, 1, 3],
             [4, 2, 3, 1],
             [4, 3, 1, 2],
             [4, 3, 2, 1],
             [3, 4, 1, 2],
             [3, 4, 2, 1],
             [2, 3, 4, 1],
             [3, 2, 4, 1],
             [1, 2, 3, 4],
             [2, 1, 3, 4],
             [2, 3, 1, 4],
             [3, 1, 2, 4],
             [3, 2, 1, 4]]
        """
        if self.n == 0:
            return

        elif self.n < 3:
            for p in itertools.permutations(range(1, self.n + 1)):
                yield self.element_class(self, p)
            return

        elif self.n == 3:
            for p in itertools.permutations(range(1, self.n + 1)):
                if p != (1, 3, 2):
                    yield self.element_class(self, p)
            return

        #Yield all the 132 avoiding permutations to the right.
        for right in StandardPermutations_avoiding_132(self.n - 1):
            yield self.element_class(self, [self.n] + list(right))

        #yi
        for i in range(1, self.n-1):
            for left in StandardPermutations_avoiding_132(i):
                for right in StandardPermutations_avoiding_132(self.n-i-1):
                    yield self.element_class(self, [x+(self.n-i-1) for x in left] + [self.n] + list(right) )


        #Yield all the 132 avoiding permutations to the left
        for left in StandardPermutations_avoiding_132(self.n - 1):
            yield self.element_class(self, list(left) + [self.n])

class StandardPermutations_avoiding_123(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[2, 1, 3])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_123, self).__init__(n, (Permutations()([1, 2, 3]),))

    def cardinality(self) -> Integer:
        """
        EXAMPLES::

            sage: Permutations(5, avoiding=[1, 2, 3]).cardinality()
            42
            sage: len( Permutations(5, avoiding=[1, 2, 3]).list() )
            42
        """
        return catalan_number(self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[1, 2, 3]).list() # indirect doctest
             [[1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
            sage: Permutations(2, avoiding=[1, 2, 3]).list()
            [[1, 2], [2, 1]]
            sage: Permutations(3, avoiding=[1, 2, 3]).list()
            [[1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        if self.n == 0:
            return

        elif self.n < 3:
            for p in itertools.permutations(range(1, self.n + 1)):
                yield self.element_class(self, p)
            return

        elif self.n == 3:
            for p in itertools.permutations(range(1, self.n + 1)):
                if p != (1, 2, 3):
                    yield self.element_class(self, p)
            return

        for p in StandardPermutations_avoiding_132(self.n):
            #Convert p to a 123 avoiding permutation by
            m = self.n+1
            minima_pos = []
            minima = []
            for i in range(self.n):
                if p[i] < m:
                    minima_pos.append(i)
                    minima.append(p[i])
                    m = p[i]

            new_p = []
            non_minima = [x for x in range(self.n, 0, -1) if x not in minima]
            a = 0
            b = 0
            for i in range(self.n):
                if i in minima_pos:
                    new_p.append( minima[a] )
                    a += 1
                else:
                    new_p.append( non_minima[b] )
                    b += 1

            yield self.element_class(self, new_p)

class StandardPermutations_avoiding_321(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[3, 2, 1])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_321, self).__init__(n, (Permutations()([3, 2, 1]),))

    def cardinality(self):
        """
        EXAMPLES::

            sage: Permutations(5, avoiding=[3, 2, 1]).cardinality()
            42
            sage: len( Permutations(5, avoiding=[3, 2, 1]).list() )
            42
        """
        return catalan_number(self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[3, 2, 1]).list() #indirect doctest
            [[2, 3, 1], [3, 1, 2], [1, 3, 2], [2, 1, 3], [1, 2, 3]]
        """
        for p in StandardPermutations_avoiding_123(self.n):
            yield self.element_class(self, p.reverse())

class StandardPermutations_avoiding_231(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[2, 3, 1])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_231, self).__init__(n, (Permutations()([2, 3, 1]),))

    def cardinality(self):
        """
        EXAMPLES::

            sage: Permutations(5, avoiding=[2, 3, 1]).cardinality()
            42
            sage: len( Permutations(5, avoiding=[2, 3, 1]).list() )
            42
        """
        return catalan_number(self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[2, 3, 1]).list()
            [[3, 2, 1], [3, 1, 2], [1, 3, 2], [2, 1, 3], [1, 2, 3]]
        """
        for p in StandardPermutations_avoiding_132(self.n):
            yield self.element_class(self, p.reverse())


class StandardPermutations_avoiding_312(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[3, 1, 2])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_312, self).__init__(n, (Permutations()([3, 1, 2]),))

    def cardinality(self):
        """
        EXAMPLES::

            sage: Permutations(5, avoiding=[3, 1, 2]).cardinality()
            42
            sage: len( Permutations(5, avoiding=[3, 1, 2]).list() )
            42
        """
        return catalan_number(self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[3, 1, 2]).list()
            [[3, 2, 1], [2, 3, 1], [2, 1, 3], [1, 3, 2], [1, 2, 3]]
        """
        for p in StandardPermutations_avoiding_132(self.n):
            yield self.element_class(self, p.complement())


class StandardPermutations_avoiding_213(StandardPermutations_avoiding_generic):
    def __init__(self, n):
        """
        TESTS::

            sage: P = Permutations(3, avoiding=[2, 1, 3])
            sage: TestSuite(P).run()
        """
        super(StandardPermutations_avoiding_213, self).__init__(n, (Permutations()([2, 1, 3]),))

    def cardinality(self):
        """
        EXAMPLES::

            sage: Permutations(5, avoiding=[2, 1, 3]).cardinality()
            42
            sage: len( Permutations(5, avoiding=[2, 1, 3]).list() )
            42
        """
        return catalan_number(self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Permutations(3, avoiding=[2, 1, 3]).list()
            [[1, 2, 3], [1, 3, 2], [3, 1, 2], [2, 3, 1], [3, 2, 1]]
        """
        for p in StandardPermutations_avoiding_132(self.n):
            yield p.complement().reverse()


class PatternAvoider(GenericBacktracker):
    def __init__(self, parent, patterns):
        """
        EXAMPLES::

            sage: from sage.combinat.permutation import PatternAvoider
            sage: P = Permutations(4)
            sage: p = PatternAvoider(P, [[1,2,3]])
            sage: loads(dumps(p))
            <sage.combinat.permutation.PatternAvoider object at 0x...>
        """
        GenericBacktracker.__init__(self, [], 1)
        self._patterns = patterns
        self._parent = parent

    def _rec(self, obj, state):
        """
        EXAMPLES::

            sage: from sage.combinat.permutation import PatternAvoider
            sage: P = Permutations(4)
            sage: p = PatternAvoider(P, [[1,2]])
            sage: list(p._rec([1], 2))
            [([2, 1], 3, False)]
        """
        i = state

        if state != self._parent.n:
            new_state = state + 1
            yld = False
        else:
            new_state = None
            yld = True

        for pos in reversed(range(len(obj)+1)):
            new_obj = self._parent.element_class(self._parent, obj[:pos] + [i] + obj[pos:])
            if all( not new_obj.has_pattern(p) for p in self._patterns):
                yield new_obj, new_state, yld


class PermutationsNK(Permutations_setk):
    """
    This exists solely for unpickling ``PermutationsNK`` objects created
    with Sage <= 6.3.
    """
    def __setstate__(self, state):
        r"""
        For unpickling old ``PermutationsNK`` objects.

        EXAMPLES::

            sage: loads(b"x\x9cM\x90\xcdN\xc30\x10\x84\xd5B\x0bM\x81\xf2\xd3\x1ex"
            ....:   b"\x03\xb8\xe4\x80x\x8bJ\x16B\xf2y\xb5qV\xa9\x95\xd8\xce"
            ....:   b"\xda[$\x0eHp\xe0\xc0[\xe3\xb4j\xe1bi\xfd\xcd\x8cg\xfd96"
            ....:   b"\t\x1b*Mp\x95\xf5(eO\xd1m\x05\xc5\x06\x0f\xbe-^\xfe\xc6"
            ....:   b"\xa4\xd6\x05\x8f\x1e\xbfx\xfc\xc1'\x0f\xba\x00r\x15\xd5"
            ....:   b"\xb5\xf5\r\x9f*\xbd\x04\x13\xfc\x1bE\x01G\xb2\t5xt\xc4"
            ....:   b"\x13\xa5\xa7`j\x14\xe4\xa9\xd230(\xd4\x84\xf8\xceg\x03"
            ....:   b"\x18$\x89\xcf\x95\x1e\x83\xe7\xd9\xbeH\xccy\xa9\xb4>\xeb"
            ....:   b"(\x16\x0e[\x82\xc3\xc0\x85\x1e=\x7f\xbf\xf2\\\xcf\xa1!O"
            ....:   b"\x11%\xc4\xc4\x17\x83\xbf\xe5\xcbM\xc6O\x19_\xe9\tT\x98"
            ....:   b"\x88\x17J/\xa0\xb7\xa6\xed\x08r\xb3\x94w\xe0\xeb\xf5(W"
            ....:   b"\xa5\x8e\x1cy\x19*'\x89[\x93s\xf8F\xe9U~\xca\x8a\xc5\xee"
            ....:   b"\xb8Kg\x93\xf0\xad\xd2\xf7G\xcb\xa0\x80\x1eS\xcaG\xcc\x17"
            ....:   b"|\xf7\x93\x03\x0f>4\xbb\x8f\xdb\xd9\x96\xea\x1f0\x81\xa2"
            ....:   b"\xa1=X\xa9mU\xfe\x02=\xaa\x87\x14")
            Permutations of the set [0, 1, 2, 3] of length 2
        """
        self.__class__ = Permutations_setk
        self.__init__(tuple(range(state['_n'])), state['_k'])

from sage.misc.persist import register_unpickle_override
register_unpickle_override("sage.combinat.permutation", "Permutation_class", Permutation)
register_unpickle_override("sage.combinat.permutation", "CyclicPermutationsOfPartition_partition", CyclicPermutationsOfPartition)
register_unpickle_override("sage.combinat.permutation", "CyclicPermutations_mset", CyclicPermutations)
register_unpickle_override('sage.combinat.permutation_nk', 'PermutationsNK', PermutationsNK)
