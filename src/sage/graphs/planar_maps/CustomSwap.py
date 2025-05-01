from sage.graphs.planar_maps.MapPermutation import *


class CustomSwap (MapPermutation):
    """
    A custom class that permit to Initializes faster simple transposition permutation
    """

    def __init__(self, lst) -> None:
        try:
            assert len(lst) == 1
            assert len(lst[0]) == 2

            self.a = min(lst[0][0], lst[0][1])
            self.b = max(lst[0][0], lst[0][1])

        except BaseException:
            raise

    def __eq__(self, other):
        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False

    def size(self):
        """
        Returns:
            The size of self
        ----
        O(1)
        """
        return self.b

    def apply(self, i):
        """
        Args:
            i
        ---
        Returns:
            self(i)
        ---
        O(1)
        """
        if i != self.a and i != self.b:
            return i
        return self.a + self.b - i

    def inverseApply(self, i):
        """
        Args:
            i an index
        ----
        Returns:
            j such that self(j) = i
        ----
        O(1)
        """
        return self(i)

    def number_of_fixed_points(self):
        """
        Returns: the number of fixed point ( we only consider i such that i<=self.size())
        ----
        O(1)
        """

        return self.b - (1 + self.a != self.b)

    def to_cycles(self):
        """
        Returns:
            A list of tuple representing the cycle decomposition of self
        ---
        O(n)
        """
        if self.a == self.b:
            return [(i,) for i in range(1, self.b + 1)]

        return [(i,) for i in range(1, self.a)] + [(self.a, self.b)] + [(j,)
                                                                        for j in range(self.a, self.b)]

    def inverse(self):
        """ 
        Returns: 
            The inverse of self
        ---
        O(1)
        """
        return self
