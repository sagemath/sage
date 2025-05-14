from sage.graphs.maps.permutation_utils_abstractor import *
from sage.graphs.maps.map_error import NotImplementedError


class PrimitiveRotatingPermutationUtilsAbstractor(PermutationUtilsAbstractor):
    def __init__(self, rpermutation) -> None:
        """
        Init the PrimitiveRotatingPermutationUtilsAbstractor

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: from sage.graphs.maps.primitive_rotating_permutation_utils_abstractor import PrimitiveRotatingPermutationUtilsAbstractor
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: tAbstr = PrimitiveRotatingPermutationUtilsAbstractor(rperm)

        NOTE:
            O(1)

        """
        self.rpermutation = rpermutation

    def numberInCycle(self, index):
        """
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor


        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: from sage.graphs.maps.primitive_rotating_permutation_utils_abstractor import PrimitiveRotatingPermutationUtilsAbstractor
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: tAbstr = PrimitiveRotatingPermutationUtilsAbstractor(rperm)
            sage: try:
            ....:     tAbstr.numberInCycle(1)
            ....: except:
            ....:     print("OK")
            ....:
            OK

        """
        raise NotImplementedError(self)

    def sameCycle(self, i, j):
        """
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: from sage.graphs.maps.primitive_rotating_permutation_utils_abstractor import PrimitiveRotatingPermutationUtilsAbstractor
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: tAbstr = PrimitiveRotatingPermutationUtilsAbstractor(rperm)
            sage: try:
            ....:     tAbstr.sameCycle(1,2)
            ....: except:
            ....:     print("OK")
            ....:
            OK

        """
        raise NotImplementedError(self)

    def numberOfCycles(self):
        """
        OUTPUT:

            The number of cycles of the permutation


        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: from sage.graphs.maps.primitive_rotating_permutation_utils_abstractor import PrimitiveRotatingPermutationUtilsAbstractor
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: tAbstr = PrimitiveRotatingPermutationUtilsAbstractor(rperm)
            sage: tAbstr.numberOfCycles()
            3

        NOTE:
            O(1)
        """
        return self.rpermutation.numberOfCycles()

    def numberOfFixedPoint(self):
        """
        OUTPUT:
            The number of fixed point of the permutation

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: from sage.graphs.maps.primitive_rotating_permutation_utils_abstractor import PrimitiveRotatingPermutationUtilsAbstractor
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: tAbstr = PrimitiveRotatingPermutationUtilsAbstractor(rperm)
            sage: tAbstr.numberOfFixedPoint()
            1

        NOTE:
            O(1)
        """
        return self.rpermutation.number_of_fixed_points()

    def checkTwoInTheSameCycle(self, listIndexes):
        """ 
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: from sage.graphs.maps.primitive_rotating_permutation_utils_abstractor import PrimitiveRotatingPermutationUtilsAbstractor
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: tAbstr = PrimitiveRotatingPermutationUtilsAbstractor(rperm)
            sage: try:
            ....:     tAbstr.checkTwoInTheSameCycle(1)
            ....: except:
            ....:     print("OK")
            ....:
            OK

        """
        raise NotImplementedError(self)
