from sage.graphs.maps.map_permutation import *


class CustomSwap (MapPermutation):
    """
    A custom class that permit to Initializes faster simple transposition permutation
    """

    def __init__(self, lst) -> None:
        """

        Initialize the CustomSwap 

        INPUT:
            lst a list of the form [(a,b)] where a,b are integer a,b>=1

        EXAMPLES::

            sage: CustomSwap([(12,8)])
            [1, 2, 3, 4, 5, 6, 7, 12, 9, 10, 11, 8]

        .. NOTE::
            O(1)
        """
        try:
            assert len(lst) == 1
            assert len(lst[0]) == 2

            self.a = min(lst[0][0], lst[0][1])
            self.b = max(lst[0][0], lst[0][1])

        except BaseException:
            raise

    def __eq__(self, other):
        """
        INPUT:
            other,another object

        OUTPUT:
            A boolean indicating if self equal other structurally

        EXAMPLES:: 

            sage: c = CustomSwap([(12,8)])
            sage: c == c
            True

        .. NOTE::

            O(n) where n is the size of the permutation.

        """
        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False

    def size(self):
        """
        OUTPUT:
            The size of self

        EXAMPLES::

            sage: CustomSwap([(12,8)]).size()
            12

        .. NOTE::
            O(1)
        """
        return self.b

    def apply(self, i):
        """
        INPUT:
            i an index

        OUTPUT:
            self(i)

        EXAMPLES::

            sage: CustomSwap([(12,8)]).apply(8)
            12

        .. NOTE::
            O(1)
        """
        if i != self.a and i != self.b:
            return i
        return self.a + self.b - i

    def inverseApply(self, i):
        """
        INPUT:
            i an index

        OUTPUT:
            j such that self(j) = i

        EXAMPLES::

            sage: CustomSwap([(12,8)]).inverseApply(8)
            12 

        .. NOTE::
            O(1)
        """
        return self(i)

    def number_of_fixed_points(self):
        """
        OUTPUT: 
            the number of fixed point ( we only consider i such that i<=self.size())


        EXAMPLES::

            sage: CustomSwap([(12,8)]).number_of_fixed_points()
            10 

        .. NOTE::
            O(1)
        """

        return self.b - (self.a != self.b)*2

    def to_cycles(self):
        """
        OUTPUT:
            A list of tuple representing the cycle decomposition of self

        EXAMPLES::

            sage: CustomSwap([(12,8)]).to_cycles()
            [(1,), (2,), (3,), (4,), (5,), (6,), (7,), (8, 12), (8,), (9,), (10,), (11,)]

        .. NOTE::
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

            sage: CustomSwap([(12,8)]).inverse()
            [1, 2, 3, 4, 5, 6, 7, 12, 9, 10, 11, 8]

        .. NOTE::
            O(1)
        """
        return self
