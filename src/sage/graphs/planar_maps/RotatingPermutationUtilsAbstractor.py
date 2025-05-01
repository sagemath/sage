
from sage.graphs.planar_maps.PermutationUtilsAbstractor import *


class RotatingPermutationUtilsAbstractor(PermutationUtilsAbstractor):
    def __init__(self, rpermutation) -> None:
        self.rpermutation = rpermutation

    def numberInCycle(self, index):
        """
        Returns:
            The size of the cycle containing index
        ---
        O(log(m)) where m is the size of the permutation 
        """
        return self.rpermutation.numberInCycle(index)

    def sameCycle(self, i, j):
        """
        Returns:
            A boolean indicating if i and j are on the same cycle
        ---
        O(log(m)) where m is the size of the permutation
        """
        return self.rpermutation.sameCycle(i, j)

    def numberOfCycles(self):
        """
        Returns:
            The number of cycles of the permutation
        ---
        O(1)
        """

        return self.rpermutation.numberOfCycles()

    def numberOfFixedPoint(self):
        """
        Returns:
            The number of fixed point of the permutation
        ---
        O(1)
        """

        return self.rpermutation.number_of_fixed_points()

    def checkTwoInTheSameCycle(self, listIndexes):
        """
        Args:
            listIndexes a list of indexes
        Returns:
            A boolean indicating if there are two indices in listIndexes on the sameCycle
        ---
        O(len(listIndexes)log(m)) where m is the permutation size
        """

        return self.rpermutation.checkTwoInTheSameCycle(listIndexes)
