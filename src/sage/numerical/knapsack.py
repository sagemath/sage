r"""
Knapsack Problems

This module implements a number of solutions to various knapsack problems,
otherwise known as linear integer programming problems. Solutions to the
following knapsack problems are implemented:

- Solving the subset sum problem for super-increasing sequences.
- General case using Linear Programming

AUTHORS:

- Minh Van Nguyen (2009-04): initial version
- Nathann Cohen (2009-08): Linear Programming version

Definition of Knapsack problems
-------------------------------

You have already had a knapsack problem, so you should know, but in case you do
not, a knapsack problem is what happens when you have hundred of items to put
into a bag which is too small, and you want to pack the most useful of them.

When you formally write it, here is your problem:

* Your bag can contain a weight of at most `W`.
* Each item `i` has a weight `w_i`.
* Each item `i` has a usefulness `u_i`.

You then want to maximize the total usefulness of the items you will store into
your bag, while keeping sure the weight of the bag will not go over `W`.

As a linear program, this problem can be represented this way (if you define
`b_i` as the binary variable indicating whether the item `i` is to be included
in your bag):

.. MATH::

    \mbox{Maximize: }\sum_i b_i u_i \\
    \mbox{Such that: }
    \sum_i b_i w_i \leq W \\
    \forall i, b_i \mbox{ binary variable} \\

(For more information, see the :wikipedia:`Knapsack_problem`)

Examples
--------

If your knapsack problem is composed of three items (weight, value)
defined by (1,2), (1.5,1), (0.5,3), and a bag of maximum weight 2,
you can easily solve it this way::

    sage: from sage.numerical.knapsack import knapsack
    sage: knapsack( [(1,2), (1.5,1), (0.5,3)], max=2)
    [5.0, [(1, 2), (0.500000000000000, 3)]]

Super-increasing sequences
--------------------------

We can test for whether or not a sequence is super-increasing::

    sage: from sage.numerical.knapsack import Superincreasing
    sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
    sage: seq = Superincreasing(L)
    sage: seq
    Super-increasing sequence of length 8
    sage: seq.is_superincreasing()
    True
    sage: Superincreasing().is_superincreasing([1,3,5,7])
    False

Solving the subset sum problem for a super-increasing sequence
and target sum::

    sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
    sage: Superincreasing(L).subset_sum(98)
    [69, 21, 5, 2, 1]
"""

#*****************************************************************************
# Copyright (C) 2009 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
#
# http://www.gnu.org/licenses/
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#*****************************************************************************

from sage.misc.latex import latex
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject


class Superincreasing(SageObject):
    r"""
    A class for super-increasing sequences.

    Let `L = (a_1, a_2, a_3, \dots, a_n)` be a non-empty sequence of
    nonnegative integers. Then `L` is said to be super-increasing if
    each `a_i` is strictly greater than the sum of all previous values.
    That is, for each `a_i \in L` the sequence `L` must satisfy the property

    .. MATH::

        a_i > \sum_{k=1}^{i-1} a_k

    in order to be called a super-increasing sequence, where `|L| \geq 2`.
    If `L` has only one element, it is also defined to be a
    super-increasing sequence.

    If ``seq`` is ``None``, then construct an empty sequence. By
    definition, this empty sequence is not super-increasing.

    INPUT:

    - ``seq`` -- (default: ``None``) a non-empty sequence

    EXAMPLES::

        sage: from sage.numerical.knapsack import Superincreasing
        sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
        sage: Superincreasing(L).is_superincreasing()
        True
        sage: Superincreasing().is_superincreasing([1,3,5,7])
        False
        sage: seq = Superincreasing(); seq
        An empty sequence.
        sage: seq = Superincreasing([1, 3, 6]); seq
        Super-increasing sequence of length 3
        sage: seq = Superincreasing([1, 2, 5, 21, 69, 189, 376, 919]); seq
        Super-increasing sequence of length 8
    """

    def __init__(self, seq=None):
        r"""
        Constructing a super-increasing sequence object from ``seq``.

        If ``seq`` is ``None``, then construct an empty sequence. By
        definition, this empty sequence is not super-increasing.


        INPUT:

        - ``seq`` -- (default: ``None``) a non-empty sequence

        EXAMPLES::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
            sage: Superincreasing(L).is_superincreasing()
            True
            sage: Superincreasing().is_superincreasing([1,3,5,7])
            False
        """
        # argument seq is None, so construct an empty sequence
        if seq is None:
            self._seq = None
        # now seq is known to be not None, so try to construct a
        # super-increasing sequence from seq
        else:
            if self.is_superincreasing(seq):
                self._seq = seq
            else:
                raise ValueError("seq must be a super-increasing sequence")

    def __eq__(self, other):
        r"""
        Comparing ``self`` to ``other``.

        TESTS::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
            sage: seq = Superincreasing(L)
            sage: seq == loads(dumps(seq))
            True
        """
        return self._seq == other._seq

    def __ne__(self, other):
        r"""
        Comparing ``self`` to ``other``.

        TESTS::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
            sage: seq = Superincreasing(L)
            sage: seq != seq
            False
        """
        return not self.__eq__(other)

    def __repr__(self):
        r"""
        Return a string representation of this super-increasing
        sequence object.

        EXAMPLES::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: seq = Superincreasing(); seq
            An empty sequence.
            sage: seq = Superincreasing([1, 3, 6]); seq
            Super-increasing sequence of length 3
            sage: seq = Superincreasing([1, 2, 5, 21, 69, 189, 376, 919]); seq
            Super-increasing sequence of length 8
        """
        if self._seq is None:
            return "An empty sequence."
        else:
            return "Super-increasing sequence of length %s" % len(self._seq)

    def largest_less_than(self, N):
        r"""
        Return the largest integer in the sequence ``self`` that is less than
        or equal to ``N``.

        This function narrows down the candidate solution using a binary trim,
        similar to the way binary search halves the sequence at each iteration.


        INPUT:

        - ``N`` -- integer; the target value to search for

        OUTPUT:

        The largest integer in ``self`` that is less than or equal to ``N``. If
        no solution exists, then return ``None``.

        EXAMPLES:

        When a solution is found, return it::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [2, 3, 7, 25, 67, 179, 356, 819]
            sage: Superincreasing(L).largest_less_than(207)
            179
            sage: L = (2, 3, 7, 25, 67, 179, 356, 819)
            sage: Superincreasing(L).largest_less_than(2)
            2

        But if no solution exists, return ``None``::

            sage: L = [2, 3, 7, 25, 67, 179, 356, 819]
            sage: Superincreasing(L).largest_less_than(-1) is None
            True

        TESTS:

        The target ``N`` must be an integer::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [2, 3, 7, 25, 67, 179, 356, 819]
            sage: Superincreasing(L).largest_less_than(2.30)
            Traceback (most recent call last):
            ...
            TypeError: N (= 2.30000000000000) must be an integer.

        The sequence that ``self`` represents must also be non-empty::

            sage: Superincreasing([]).largest_less_than(2)
            Traceback (most recent call last):
            ...
            ValueError: seq must be a super-increasing sequence
            sage: Superincreasing(list()).largest_less_than(2)
            Traceback (most recent call last):
            ...
            ValueError: seq must be a super-increasing sequence
        """
        from sage.functions.other import Function_floor
        floor = Function_floor()
        # input error handling
        if len(self._seq) == 0:
            raise TypeError("self must be a non-empty list of integers.")
        if (not isinstance(N, Integer)) and (not isinstance(N, int)):
            raise TypeError("N (= %s) must be an integer." % N)

        # halving the list at each iteration, just like binary search
        # TODO: some error handling to ensure that self only contains integers?
        low = 0
        high = len(self._seq) - 1
        while low <= high:
            mid = floor((low + high) / 2)
            if N == self._seq[mid]:
                return self._seq[mid]
            if N < self._seq[mid]:
                high = mid - 1
            else:
                low = mid + 1
        if N >= self._seq[high]:
            return self._seq[high]
        else:
            return None

    def _latex_(self):
        r"""
        Return LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: latex(Superincreasing())
            \left[\right]
            sage: seq = Superincreasing([1, 2, 5, 21, 69, 189, 376, 919])
            sage: latex(seq)
            <BLANKLINE>
            \left[1,
            2,
            5,
            21,
            69,
            189,
            376,
            919\right]
        """
        if self._seq is None:
            return latex([])
        else:
            return latex(self._seq)

    def is_superincreasing(self, seq=None):
        r"""
        Determine whether or not ``seq`` is super-increasing.

        If ``seq=None`` then determine whether or not ``self`` is
        super-increasing.

        Let `L = (a_1, a_2, a_3, \dots, a_n)` be a non-empty sequence of
        nonnegative integers. Then `L` is said to be super-increasing if
        each `a_i` is strictly greater than the sum of all previous values.
        That is, for each `a_i \in L` the sequence `L` must satisfy the
        property

        .. MATH::

                  a_i > \sum_{k=1}^{i-1} a_k

        in order to be called a super-increasing sequence, where `|L| \geq 2`.
        If `L` has exactly one element, then it is also defined to be a
        super-increasing sequence.

        INPUT:

        - ``seq`` -- (default: ``None``) a sequence to test

        OUTPUT:

        - If ``seq`` is ``None``, then test ``self`` to determine whether or
          not it is super-increasing. In that case, return ``True`` if
          ``self`` is super-increasing; ``False`` otherwise.

        - If ``seq`` is not ``None``, then test ``seq`` to determine whether
          or not it is super-increasing. Return ``True`` if ``seq`` is
          super-increasing; ``False`` otherwise.

        EXAMPLES:

        By definition, an empty sequence is not super-increasing::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: Superincreasing().is_superincreasing([])
            False
            sage: Superincreasing().is_superincreasing()
            False
            sage: Superincreasing().is_superincreasing(tuple())
            False
            sage: Superincreasing().is_superincreasing(())
            False

        But here is an example of a super-increasing sequence::

            sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
            sage: Superincreasing(L).is_superincreasing()
            True
            sage: L = (1, 2, 5, 21, 69, 189, 376, 919)
            sage: Superincreasing(L).is_superincreasing()
            True

        A super-increasing sequence can have zero as one of its elements::

            sage: L = [0, 1, 2, 4]
            sage: Superincreasing(L).is_superincreasing()
            True

        A super-increasing sequence can be of length 1::

            sage: Superincreasing([randint(0, 100)]).is_superincreasing()
            True


        TESTS:

        The sequence must contain only integers::

            sage: # needs sage.symbolic
            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [1.0, 2.1, pi, 21, 69, 189, 376, 919]
            sage: Superincreasing(L).is_superincreasing()
            Traceback (most recent call last):
            ...
            TypeError: Element e (= 1.00000000000000) of seq must be a nonnegative integer.
            sage: L = [1, 2.1, pi, 21, 69, 189, 376, 919]
            sage: Superincreasing(L).is_superincreasing()
            Traceback (most recent call last):
            ...
            TypeError: Element e (= 2.10000000000000) of seq must be a nonnegative integer.
        """
        # argument seq is None, so test self for super-increasing
        if seq is None:
            # self must be a non-empty sequence
            if (self._seq is None) or len(self._seq) == 0:
                return False
            # so now self is known to represent a non-empty sequence
            if (not isinstance(self._seq[0], Integer)) and (not isinstance(self._seq[0], int)):
                raise TypeError("Element e (= %s) of self must be a non-negative integer." % self._seq[0])
            if self._seq[0] < 0:
                raise TypeError("Element e (= %s) of self must be a non-negative integer." % self._seq[0])
            cumSum = self._seq[0]  # the cumulative sum of the sequence represented by self
            for e in self._seq[1:]:
                if (not isinstance(e, Integer)) and (not isinstance(e, int)):
                    raise TypeError("Element e (= %s) of self must be a non-negative integer." % e)
                if e < 0:
                    raise TypeError("Element e (= %s) of self must be a non-negative integer." % e)
                if e <= cumSum:
                    return False
                cumSum += e
            return True
        # now we know that seq is not None, so test seq for super-increasing
        else:
            # seq must be a non-empty sequence
            if len(seq) == 0:
                return False
            # so now seq is known to represent a non-empty sequence
            if (not isinstance(seq[0], Integer)) and (not isinstance(seq[0], int)):
                raise TypeError("Element e (= %s) of seq must be a nonnegative integer." % seq[0])
            if seq[0] < 0:
                raise TypeError("Element e (= %s) of seq must be a nonnegative integer." % seq[0])
            cumSum = seq[0]  # the cumulative sum of the sequence seq
            for e in seq[1:]:
                if (not isinstance(e, Integer)) and (not isinstance(e, int)):
                    raise TypeError("Element e (= %s) of seq must be a nonnegative integer." % e)
                if e < 0:
                    raise TypeError("Element e (= %s) of seq must be a nonnegative integer." % e)
                if e <= cumSum:
                    return False
                cumSum += e
            return True

    def subset_sum(self, N):
        r"""
        Solving the subset sum problem for a super-increasing sequence.

        Let `S = (s_1, s_2, s_3, \dots, s_n)` be a non-empty sequence of
        nonnegative integers, and let `N \in \ZZ` be nonnegative. The
        subset sum problem asks for a  subset `A \subseteq S` all of whose
        elements sum to `N`. This method specializes the subset sum problem
        to the case of super-increasing sequences. If a solution exists, then
        it is also a super-increasing sequence.

        .. NOTE::

            This method only solves the subset sum problem for
            super-increasing sequences. In general, solving the subset sum
            problem for an arbitrary sequence is known to be computationally
            hard.

        INPUT:

        - ``N`` -- nonnegative integer

        OUTPUT:

        - A non-empty subset of ``self`` whose elements sum to ``N``. This
          subset is also a super-increasing sequence. If no such subset
          exists, then return the empty list.

        ALGORITHM:

        The algorithm used is adapted from page 355 of [HPS2008]_.

        EXAMPLES:

        Solving the subset sum problem for a super-increasing sequence
        and target sum::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [1, 2, 5, 21, 69, 189, 376, 919]
            sage: Superincreasing(L).subset_sum(98)
            [69, 21, 5, 2, 1]


        TESTS:

        The target ``N`` must be a nonnegative integer::

            sage: from sage.numerical.knapsack import Superincreasing
            sage: L = [0, 1, 2, 4]
            sage: Superincreasing(L).subset_sum(-6)
            Traceback (most recent call last):
            ...
            TypeError: N (= -6) must be a nonnegative integer.
            sage: Superincreasing(L).subset_sum(-6.2)
            Traceback (most recent call last):
            ...
            TypeError: N (= -6.20000000000000) must be a nonnegative integer.

        The sequence that ``self`` represents must only contain nonnegative
        integers::

            sage: L = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1]
            sage: Superincreasing(L).subset_sum(1)
            Traceback (most recent call last):
            ...
            TypeError: Element e (= -10) of seq must be a nonnegative integer.
        """
        # input error handling
        if not self.is_superincreasing():
            raise TypeError("self is not super-increasing. Only super-increasing sequences are currently supported.")
        if (not isinstance(N, Integer)) and (not isinstance(N, int)):
            raise TypeError("N (= %s) must be a nonnegative integer." % N)
        if N < 0:
            raise TypeError("N (= %s) must be a nonnegative integer." % N)

        # solve subset sum problem for super-increasing sequence
        candidates = []
        a = self.largest_less_than(N)
        while a is not None:
            candidates.append(a)
            a = self.largest_less_than(N - sum(candidates))

        lst = list(set(candidates))  # removing any duplicate elements
        if len(lst) != len(candidates):
            return []
        if sum(candidates) == N:
            return candidates
        else:
            return []


def knapsack(seq, binary=True, max=1, value_only=False, solver=None, verbose=0,
             *, integrality_tolerance=1e-3):
    r"""
    Solve the knapsack problem.

    For more information on the knapsack problem, see the documentation of the
    :mod:`knapsack module <sage.numerical.knapsack>` or the
    :wikipedia:`Knapsack_problem`.

    INPUT:

    - ``seq`` -- two different possible types:

      - A sequence of tuples ``(weight, value, something1, something2,
        ...)``. Note that only the first two coordinates (``weight`` and
        ``values``) will be taken into account. The rest (if any) will be
        ignored. This can be useful if you need to attach some information to
        the items.

      - A sequence of reals (a value of 1 is assumed).

    - ``binary`` -- when set to ``True``, an item can be taken 0 or 1 time
      When set to ``False``, an item can be taken any amount of times (while
      staying integer and positive).

    - ``max`` -- maximum admissible weight

    - ``value_only`` -- when set to ``True``, only the maximum useful value is
      returned. When set to ``False``, both the maximum useful value and an
      assignment are returned.

    - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear Programming
      (MILP) solver to be used. If set to ``None``, the default one is used. For
      more information on MILP solvers and which default solver is used, see
      the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
      of the class
      :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- parameter for use with MILP solvers over an
      inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`

    OUTPUT:

    If ``value_only`` is set to ``True``, only the maximum useful value is
    returned. Else (the default), the function returns a pair ``[value,list]``,
    where ``list`` can be of two types according to the type of ``seq``:

    - The list of tuples `(w_i, u_i, ...)` occurring in the solution.

    - A list of reals where each real is repeated the number of times it is
      taken into the solution.

    EXAMPLES:

    If your knapsack problem is composed of three items ``(weight, value)``
    defined by ``(1,2), (1.5,1), (0.5,3)``, and a bag of maximum weight `2`, you
    can easily solve it this way::

        sage: from sage.numerical.knapsack import knapsack
        sage: knapsack( [(1,2), (1.5,1), (0.5,3)], max=2)
        [5.0, [(1, 2), (0.500000000000000, 3)]]

        sage: knapsack( [(1,2), (1.5,1), (0.5,3)], max=2, value_only=True)
        5.0

    Besides weight and value, you may attach any data to the items::

        sage: from sage.numerical.knapsack import knapsack
        sage: knapsack( [(1, 2, 'spam'), (0.5, 3, 'a', 'lot')])
        [3.0, [(0.500000000000000, 3, 'a', 'lot')]]

    In the case where all the values (usefulness) of the items are equal to one,
    you do not need embarrass yourself with the second values, and you can just
    type for items `(1,1), (1.5,1), (0.5,1)` the command::

        sage: from sage.numerical.knapsack import knapsack
        sage: knapsack([1,1.5,0.5], max=2, value_only=True)
        2.0
    """
    reals = not isinstance(seq[0], tuple)
    if reals:
        seq = [(x,1) for x in seq]

    from sage.numerical.mip import MixedIntegerLinearProgram
    from sage.rings.integer_ring import ZZ

    p = MixedIntegerLinearProgram(solver=solver, maximization=True)

    if binary:
        present = p.new_variable(binary=True)
    else:
        present = p.new_variable(integer=True)

    p.set_objective(p.sum([present[i] * seq[i][1] for i in range(len(seq))]))
    p.add_constraint(p.sum([present[i] * seq[i][0] for i in range(len(seq))]), max=max)

    if value_only:
        return p.solve(objective_only=True, log=verbose)

    else:
        objective = p.solve(log=verbose)
        present = p.get_values(present, convert=ZZ, tolerance=integrality_tolerance)

        val = []

        if reals:
            [val.extend([seq[i][0]] * present[i]) for i in range(len(seq))]
        else:
            [val.extend([seq[i]] * present[i]) for i in range(len(seq))]

        return [objective,val]
