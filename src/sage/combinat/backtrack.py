r"""
Backtracking

This library contains a generic tool for constructing large sets whose
elements can be enumerated by exploring a search space with a (lazy)
tree or graph structure.

- :class:`GenericBacktracker`: Depth first search through a tree
  described by a ``children`` function, with branch pruning, etc.

This module has mostly been superseded by ``RecursivelyEnumeratedSet``.
"""
# ****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                     2010 Nicolas Borie <nicolas.borie at math.u-psud.fr>
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
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.monoids import Monoids
from sage.categories.commutative_additive_semigroups import (
    CommutativeAdditiveSemigroups)
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest


class GenericBacktracker:
    r"""
    A generic backtrack tool for exploring a search space organized as a tree,
    with branch pruning, etc.

    See also :class:`RecursivelyEnumeratedSet_forest` for
    handling simple special cases.
    """

    def __init__(self, initial_data, initial_state):
        r"""
        EXAMPLES::

            sage: from sage.combinat.backtrack import GenericBacktracker
            sage: p = GenericBacktracker([], 1)
            sage: loads(dumps(p))
            <sage.combinat.backtrack.GenericBacktracker object at 0x...>
        """
        self._initial_data = initial_data
        self._initial_state = initial_state

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.permutation import PatternAvoider
            sage: p = PatternAvoider(Permutations(4), [[1,3,2]])
            sage: len(list(p))                                                          # needs sage.combinat
            14
        """
        # Initialize the stack of generators with the initial data.
        # The generator in stack[i] is a generator for the i^th level
        # of the search tree.
        stack = []
        stack.append(self._rec(self._initial_data, self._initial_state))

        done = False
        while not done:
            # Try to get the next object in this level
            try:
                obj, state, yld = next(stack[-1])
            except StopIteration:
                # If there are no more, go back up the tree
                # We also need to check if we've exhausted all
                # possibilities
                stack.pop()
                done = len(stack) == 0
                continue

            # If the return state is None, then obj is a leaf
            # of the search tree.  If yld is True, then obj
            # should be yielded.
            if yld is True:
                yield obj
            if state is not None:
                stack.append(self._rec(obj, state))


class PositiveIntegerSemigroup(UniqueRepresentation, RecursivelyEnumeratedSet_forest):
    r"""
    The commutative additive semigroup of positive integers.

    This class provides an example of algebraic structure which
    inherits from :class:`RecursivelyEnumeratedSet_forest`. It builds the positive
    integers a la Peano, and endows it with its natural commutative
    additive semigroup structure.

    EXAMPLES::

        sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
        sage: PP = PositiveIntegerSemigroup()
        sage: PP.category()
        Join of Category of monoids and Category of commutative additive semigroups and Category of infinite enumerated sets and Category of facade sets
        sage: PP.cardinality()
        +Infinity
        sage: PP.one()
        1
        sage: PP.an_element()
        1
        sage: some_elements = list(PP.some_elements()); some_elements
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]

    TESTS::

        sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
        sage: PP = PositiveIntegerSemigroup()

    We factor out the long test from the ``TestSuite``::

        sage: TestSuite(PP).run(skip='_test_enumerated_set_contains')
        sage: PP._test_enumerated_set_contains()  # long time
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
        """
        RecursivelyEnumeratedSet_forest.__init__(self, facade=ZZ, category=(InfiniteEnumeratedSets(), CommutativeAdditiveSemigroups(), Monoids()))

    def roots(self):
        r"""
        Return the single root of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: list(PP.roots())
            [1]
        """
        return [ZZ(1)]

    def children(self, x):
        r"""
        Return the single child ``x+1`` of the integer ``x``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: list(PP.children(1))
            [2]
            sage: list(PP.children(42))
            [43]
        """
        return [ZZ(x + 1)]

    def one(self):
        r"""
        Return the unit of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: PP.one()
            1
        """
        return self.first()
