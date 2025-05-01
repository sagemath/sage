import numpy as np


class PermutationUtilsAbstractor:
    """
    This class abstract some utils use in LabelledMap and MutableLabelledMap so that they can use the same
    apis but with different underlying implementation hence the version for MutableLabelledMap inherit from
    this class and is call RotatingUtilsAbstractor
    """

    def __init__(self, permutation) -> None:
        cycles = permutation.to_cycles()

        self._numberOfCycles = len(cycles)

        self._numberOfFixedPoint = sum(len(c) == 1 for c in cycles)

        self._cyclesLength = [len(c) for c in cycles]

        self._cycleIndexes = np.zeros(permutation.size() + 1, dtype=int)

        for j, c in enumerate(cycles):
            for i in c:
                self._cycleIndexes[i] = j

    def numberInCycle(self, index):
        """
        Returns:
            The size of the cycle containing index
        ---
        O(1)
        """
        return self._cyclesLength[self._cycleIndexes[index]]

    def sameCycle(self, i, j):
        """
        Returns:
            A boolean indicating if i and j are on the same cycle
        ---
        O(1)
        """
        return self._cycleIndexes[i] == self._cycleIndexes[j]

    def numberOfCycles(self):
        """
        Returns:
            The number of cycles of the permutation
        ---
        O(1)
        """
        return self._numberOfCycles

    def numberOfFixedPoint(self):
        """
        Returns:
            The number of fixed point of the permutation
        ---
        O(1)
        """
        self._numberOfFixedPoint

    def checkTwoInTheSameCycle(self, listIndexes):
        """
        Args:
            listIndexes a list of indexes
        Returns:
            A boolean indicating if there are two indices in listIndexes on the sameCycle
        ---
        O(len(listIndexes))
        """
        checkSet = set()
        for i in listIndexes:
            if self._cycleIndexes[i] in checkSet:
                return True
            checkSet.add(self._cycleIndexes[i])
        return False
