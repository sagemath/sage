r"""
The Steinhaus-Johnson-Trotter algorithm generates all permutations of a list in
an order such that each permutation is obtained by transposing two adjacent
elements from the previous permutation.

Each element of the list has a direction (initialized at -1) that changes at
each permutation and that is used to determine which elements to transpose. Thus
in addition to the permutation itself, the direction of each element is also
stored.

Note that the permutations are not generated in lexicographic order.

AUTHORS:

- Martin Grenouilloux (2024-05-22): initial version
"""

# ****************************************************************************
#       Copyright (C) 2024 Martin Grenouilloux <martin.grenouilloux@uib.no>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.combinat.combinat import CombinatorialElement

class SJT(CombinatorialElement):
    r"""
    A representation of a list permuted using the Steinhaus-Johnson-Trotter
    algorithm.

    Each element of the list has a direction (initialized at -1) that changes at
    each permutation and that is used to determine which elements to transpose.
    The directions have three possible values:

    - ``-1``: element tranposes to the left

    - ``1``: element transposes to the right

    - ``0``: element does not move

    Thus in addition to the permutation itself, the direction of each element is
    also stored.

    Note that the permutations are not generated in lexicographic order.

    .. WARNING::

        An ``SJT`` object should always be created with identity permutation for
        the algorithm to behave properly. If the identity permutation is not
        provided, it expects a coherent list of directions according to the
        provided input. This list is not checked.

    .. TODO::

        Implement the previous permutation for the Steinhaus-Johnson-Trotter
        algorithm.

    EXAMPLES::

        sage: from sage.combinat.SJT import SJT
        sage: s = SJT([1, 2, 3, 4]); s
        [1, 2, 3, 4]
        sage: s = s.next(); s
        [1, 2, 4, 3]
        sage: p = Permutation(s._list, algorithm='sjt', sjt=s)
        sage: p
        [1, 2, 4, 3]
        sage: p.next()
        [1, 4, 2, 3]

    TESTS::

        sage: from sage.combinat.SJT import SJT
        sage: s = SJT([1, 2, 3, 4]); s
        [1, 2, 3, 4]
        sage: s = SJT([1]); s
        [1]
        sage: s = s.next(); s
        False
        sage: s = SJT([]); s
        []
        sage: s = s.next(); s
        False
    """
    def __init__(self, l, directions=None) -> None:
        r"""
        Transpose two elements at positions ``a`` and ``b`` in ``perm`` and
        their corresponding directions as well following the
        Steinhaus-Johnson-Trotter algorithm.

        Each permutation is obtained by transposing two adjacent elements from
        the previous permutation.

        INPUT:

        - ``l`` -- list; a list of ordered ``int``.

        - ``directions`` -- list (default: ``None``); a list of directions for
          each element in the permuted list. Used when constructing permutations
          from a pre-defined internal state.

        EXAMPLES::

            sage: from sage.combinat.SJT import SJT
            sage: s = SJT([1, 2, 3, 4]); s
            [1, 2, 3, 4]
            sage: s = s.next(); s
            [1, 2, 4, 3]
            sage: p = Permutation(s._list, algorithm='sjt', sjt=s)
            sage: p
            [1, 2, 4, 3]
            sage: p.next()
            [1, 4, 2, 3]

        TESTS::

            sage: from sage.combinat.SJT import SJT
            sage: s = SJT([1, 3, 2, 4])
            Traceback (most recent call last):
            ...
            ValueError: no internal state directions were given for non-identity
            starting permutation for Steinhaus-Johnson-Trotter algorithm
            sage: s = SJT([]); s
            []
            sage: s = s.next(); s
            False
        """
        # The permuted list.
        self._list = l

        # The length of the permuted list. Return early on empty list.
        self._n = len(l)
        if self._n == 0:
            return

        if directions is None:
            if not all(l[i] <= l[i+1] for i in range(self._n - 1)):
                raise ValueError("no internal state directions were given for "
                "non-identity starting permutation for "
                "Steinhaus-Johnson-Trotter algorithm")
            self._directions = [-1] * self._n

            # The first element has null direction.
            self._directions[0] = 0
        else:
            self._directions = directions

    def __idx_largest_element_non_zero_direction(self, perm, directions):
        r"""
        Find the largest element in ``perm`` with a non null direction.
        """
        largest = 0
        index = None
        for i in range(self._n):
            if directions[i] != 0:
                e = perm[i]
                if e > largest:
                    index = i
                    largest = e

        return index

    def next(self):
        r"""
        Produce the next permutation of ``self`` following the
        Steinhaus-Johnson-Trotter algorithm.

        OUTPUT: the list of the next permutation

        EXAMPLES::

            sage: from sage.combinat.SJT import SJT
            sage: s = SJT([1, 2, 3, 4])
            sage: s = s.next(); s
            [1, 2, 4, 3]
            sage: s = s.next(); s
            [1, 4, 2, 3]

        TESTS::

            sage: from sage.combinat.SJT import SJT
            sage: s = SJT([1, 2, 3])
            sage: s.next()
            [1, 3, 2]

            sage: s = SJT([1])
            sage: s.next()
            False
        """
        # Return on empty list.
        if self._n == 0:
            return False

        # Copying lists of permutation and directions to avoid changing internal
        # state of the algorithm if ``next()`` is called without reassigning.
        perm = self._list[:]
        directions = self._directions[:]

        # Assume that the element to move is n (which will be in most cases).
        selected_elt = self._n
        xi = perm.index(selected_elt)
        direction = directions[xi]

        # If this element has null direction, find the largest whose is
        # non-null.
        if direction == 0:
            xi = self.__idx_largest_element_non_zero_direction(perm, directions)
            if xi is None:
                # We have created every permutation. Detected when all elements
                # have null direction.
                return False
            direction = directions[xi]
            selected_elt = perm[xi]

        new_pos = xi + direction

        # Proceed to transpose elements and corresponding directions.
        perm[xi], perm[new_pos] = perm[new_pos], perm[xi]
        directions[xi], directions[new_pos] = \
            directions[new_pos], directions[xi]

        # If the transposition results in the largest element being on one edge
        # or if the following element in its direction is greater than it, then
        # then set its direction to 0
        if new_pos == 0 or new_pos == self._n - 1 or \
            perm[new_pos + direction] > selected_elt:
            directions[new_pos] = 0

        # After each permutation, update each element's direction. If one
        # element is greater than selected element, change its direction towards
        # the selected element. This loops has no reason to be if selected
        # element is n and this will be the case most of the time.
        if selected_elt != self._n:
            for i in range(self._n):
                if perm[i] > selected_elt:
                    if i < new_pos:
                        directions[i] = 1
                    if i > new_pos:
                        directions[i] = -1

        return SJT(perm, directions)

    __next__ = next
