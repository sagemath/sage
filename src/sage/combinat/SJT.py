r"""
Steinhaus-Johnson-Trotter algorithm.

The Steinhaus-Johnson-Trotter algorithm generates permutations in a specific
order by transposing only two elements from the list at each operation. The
algorithm stores an internal state for every element in the permutated list
which corresponds to their direction. To know which elements to move, the
internal state table is accessed and the transposition of two elements is then
done.

It is important to notice that the permutations are generated in a different
order than the default lexicographic algorithm.

The class defined here is meant to be used in the ``Permutation`` class when
called with the parameter ``algorithm='sjt'``.

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
class SJT:
    def __init__(self, l, directions=None) -> None:
        r"""
        Transpose two elements at positions ``a`` and ``b`` in ``perm`` and
        their corresponding directions as well following the
        Steinhaus-Johnson-Trotter algorithm.

        INPUT:

        - ``l`` -- list: a list of ordered ``int``.

        - ``directions`` -- list (default: ``None`` ): a list of directions for
          each element in the permuted list. Used when constructing permutations
          from a pre-defined internal state. There are three possible values:

          - ``-1`` -> element tranposes to the left

          - ``1``  -> element transposes to the right

          - ``0``  -> element does not move

        """
        # The permuted list.
        self.__perm = l

        # The length of the permuted list.
        self.__n = len(l)

        if directions is None:
            if not all(l[i] <= l[i+1] for i in range(len(l) - 1)):
                raise ValueError("No internal state directions were given for "
                "non-identity starting permutation for "
                "Steinhaus-Johnson-Trotter algorithm. Expected identity "
                "permutation.")
            self.__directions = [-1] * self.__n

            # The first element has null direction.
            self.__directions[0] = 0
        else:
            self.__directions = directions

    def __idx_largest_element_non_zero_direction(self, perm, directions):
        r"""
        Find the largest element in ``perm`` with a non null direction.
        """
        largest = 0
        index = None
        for i in range(self.__n):
            if directions[i] != 0:
                e = perm[i]
                if e > largest:
                    index = i
                    largest = e

        return index

    def next(self):
        r"""
        Produce the next permutation following the Steinhaus-Johnson-Trotter
        algorithm.

        OUTPUT: the list of the next permutation.
        """
        # Copying lists of permutation and directions to avoid changing internal
        # state of the algorithm if ``next()`` is called without reassigning.
        perm = self.__perm[:]
        directions = self.__directions[:]

        # Assume that the element to move is n (which will be in most cases).
        selected_elt = self.__n
        xi = perm.index(selected_elt)
        direction = directions[xi]

        # If this element has null direction, find the largest whose is
        # non-null.
        if direction == 0:
            xi = self.__idx_largest_element_non_zero_direction(perm, directions)
            if xi is None:
                # We have created every permutation. Detected when all elements
                # have null direction.
                return False, None
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
        if new_pos == 0 or new_pos == self.__n - 1 or \
            perm[new_pos + direction] > selected_elt:
            directions[new_pos] = 0

        # After each permutation, update each element's direction. If one
        # element is greater than selected element, change its direction towards
        # the selected element. This loops has no reason to be if selected
        # element is n and this will be the case most of the time.
        if selected_elt != self.__n:
            for i in range(self.__n):
                if perm[i] > selected_elt:
                    if i < new_pos:
                        directions[i] = 1
                    if i > new_pos:
                        directions[i] = -1

        return perm, directions
