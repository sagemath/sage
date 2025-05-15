"""
Defines the internal CustomSwap class.

The CustomSwap class allows fast creation of permutations from a
list of transpositions.
"""

from typing import Any
from sage.graphs.maps.map_permutation import MapPermutation


class CustomSwap (MapPermutation):
    """An internal class which allows to initialize permutations from a list of transpositions."""

    def __init__(self, lst: list[tuple[int, int]]):
        """
        Initialize the CustomSwap structure.

        INPUT:
        - ``lst`` -- list[tuple[int,int]]: a list of transpositions of the form [(a,b)], where a,b are integers >=1

        EXAMPLES::

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: CustomSwap([(12,8)])
            [1, 2, 3, 4, 5, 6, 7, 12, 9, 10, 11, 8]

        NOTE:
            O(1); note that as of now, only arguments of the form [(i,j)] (of length 1) are supported.
            This class is not intended to be used by the user.
        """
        try:
            assert len(lst) == 1
            assert len(lst[0]) == 2

            self.a = min(lst[0][0], lst[0][1])
            self.b = max(lst[0][0], lst[0][1])

        except BaseException:
            raise

    def __eq__(self, other: Any) -> bool:
        """
        Return whether self is structurally equal to the given Permutation.
        
        INPUT:
        - ``other`` -- Any

        EXAMPLES:: 

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: c = CustomSwap([(12,8)])
            sage: c == c
            True

        NOTE:
            O(n), where n is the size of the permutation.

        """
        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False

    def size(self) -> int:
        """
        Return the size of self.

        EXAMPLES::

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: CustomSwap([(12,8)]).size()
            12

        NOTE:
            O(1)
        """
        return self.b

    def apply(self, i: int) -> int:
        """
        Return self(i).

        INPUT:
        - ``i`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: CustomSwap([(12,8)]).apply(8)
            12

        NOTE:
            O(1)
        """
        if i != self.a and i != self.b:
            return i
        return self.a + self.b - i

    def inverseApply(self, i: int) -> int:
        """
        Return self^-1(i) (ie. j such that self(j) = i).

        INPUT:
        - ``i`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: CustomSwap([(12,8)]).inverseApply(8)
            12

        NOTE:
            O(1)
        """
        return self(i)

    def number_of_fixed_points(self) -> int:
        """
        Return the number of fixed point (only considering those lower than self.size()).

        EXAMPLES::

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: CustomSwap([(12,8)]).number_of_fixed_points()
            10

        NOTE:
            O(1)
        """
        return self.b - (self.a != self.b)*2

    def to_cycles(self):
        """
        A list of tuple representing the cycle decomposition of self.

        EXAMPLES::

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: CustomSwap([(12,8)]).to_cycles()
            [(1,), (2,), (3,), (4,), (5,), (6,), (7,), (8, 12), (8,), (9,), (10,), (11,)]

        NOTE:
            O(n) where n is the size of the permutation
        """
        if self.a == self.b:
            return [(i,) for i in range(1, self.b + 1)]

        return [(i,) for i in range(1, self.a)] + [(self.a, self.b)] + [(j,)
                                                                        for j in range(self.a, self.b)]

    def inverse(self):
        """ 
        OUTPUT: 
            The inverse of self

        EXAMPLES::

            sage: from sage.graphs.maps.custom_swap import CustomSwap
            sage: CustomSwap([(12,8)]).inverse()
            [1, 2, 3, 4, 5, 6, 7, 12, 9, 10, 11, 8]

        NOTE:
            O(1)
        """
        return self
