"""Define the internal class RotatingPermutationUtilsAbstractor."""

from sage.graphs.maps.rotating_permutation import RotatingPermutation
from sage.graphs.maps.permutation_utils_abstractor import PermutationUtilsAbstractor


class RotatingPermutationUtilsAbstractor(PermutationUtilsAbstractor):
    """Internal class."""

    def __init__(self, rpermutation: RotatingPermutation) -> None:
        """
        Init the RotatingPermutationUtilsAbstractor

        INPUT:

        - ``permutation`` -- RotatingPermutation


        EXAMPLES::

            sage: from sage.graphs.maps.rotating_permutation import RotatingPermutation
            sage: from sage.graphs.maps.rotating_permutation_utils_abstractor import RotatingPermutationUtilsAbstractor
            sage: rperm = RotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rAbstractor = RotatingPermutationUtilsAbstractor(rperm)

        NOTE:

        O(1)
        """
        self.rpermutation = rpermutation

    def numberInCycle(self, index: int) -> int:
        """
        INPUT:

        - ``index`` -- int

        OUTPUT:

        The size of the cycle containing index

        EXAMPLES::

            sage: from sage.graphs.maps.rotating_permutation import RotatingPermutation
            sage: from sage.graphs.maps.rotating_permutation_utils_abstractor import RotatingPermutationUtilsAbstractor
            sage: rperm = RotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rAbstractor = RotatingPermutationUtilsAbstractor(rperm)
            sage: rAbstractor.numberInCycle(5)
            4

        NOTE:

        O(log(m)) where m is the size of the permutation
        """
        return self.rpermutation.numberInCycle(index)

    def sameCycle(self, i: int, j: int) -> bool:
        """
        INPUT:

        - ``i`` -- int
        - ``j`` -- int

        OUTPUT:

        A boolean indicating if i and j are on the same cycle

        NOTE:

        O(log(m)) where m is the size of the permutation
        """
        return self.rpermutation.sameCycle(i, j)

    def numberOfCycles(self) -> int:
        """
        OUTPUT:

        The number of cycles of the permutation

        EXAMPLES::

            sage: from sage.graphs.maps.rotating_permutation import RotatingPermutation
            sage: from sage.graphs.maps.rotating_permutation_utils_abstractor import RotatingPermutationUtilsAbstractor
            sage: rperm = RotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rAbstractor = RotatingPermutationUtilsAbstractor(rperm)
            sage: rAbstractor.numberOfCycles()
            3

        NOTE:

        O(1)
        """

        return self.rpermutation.numberOfCycles()

    def numberOfFixedPoint(self) -> int:
        """
        OUTPUT:

        The number of fixed point of the permutation

        EXAMPLES::

            sage: from sage.graphs.maps.rotating_permutation import RotatingPermutation
            sage: from sage.graphs.maps.rotating_permutation_utils_abstractor import RotatingPermutationUtilsAbstractor
            sage: rperm = RotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rAbstractor = RotatingPermutationUtilsAbstractor(rperm)
            sage: rAbstractor.numberOfFixedPoint()
            1

        NOTE:

        O(1)
        """

        return self.rpermutation.number_of_fixed_points()

    def checkTwoInTheSameCycle(self, listIndexes: list[int]) -> bool:
        """
        INPUT:

        - ``listIndexes`` -- List[int]

        OUTPUT:

        A boolean indicating if there are two indices in listIndexes on the sameCycle

        EXAMPLES::

            sage: from sage.graphs.maps.rotating_permutation import RotatingPermutation
            sage: from sage.graphs.maps.rotating_permutation_utils_abstractor import RotatingPermutationUtilsAbstractor
            sage: rperm = RotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rAbstractor = RotatingPermutationUtilsAbstractor(rperm)
            sage: rAbstractor.checkTwoInTheSameCycle([1,6,3])
            True
            sage: rAbstractor.checkTwoInTheSameCycle([1,6,7])
            False

        NOTE:

        O(len(listIndexes)log(m)) where m is the permutation size
        """

        return self.rpermutation.checkTwoInTheSameCycle(listIndexes)
