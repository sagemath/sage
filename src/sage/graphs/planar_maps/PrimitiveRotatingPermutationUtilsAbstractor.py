from sage.graphs.planar_maps.PermutationUtilsAbstractor import *
from sage.graphs.planar_maps.MapError import NotImplementedError


class PrimitiveRotatingPermutationUtilsAbstractor(PermutationUtilsAbstractor):
    def __init__(self, rpermutation) -> None:
        self.rpermutation = rpermutation

    def numberInCycle(self, index):
        """
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor
        """
        raise NotImplementedError(self)

    def sameCycle(self, i, j):
        """
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor
        """
        raise NotImplementedError(self)

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
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor
        """
        raise NotImplementedError(self)
