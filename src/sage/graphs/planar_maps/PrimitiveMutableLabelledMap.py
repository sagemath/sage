
from sage.graphs.planar_maps.LabelledMap import *
from sage.all import Permutation  # Import sage library
from sage.graphs.planar_maps.PrimitiveRotatingPermutation import PrimitiveRotatingPermutation
from sage.graphs.planar_maps.MapError import NotImplementedError
from sage.graphs.planar_maps.PrimitiveRotatingPermutationUtilsAbstractor import PrimitiveRotatingPermutationUtilsAbstractor
from sage.graphs.planar_maps.PrimitiveMutableTopologicalDemiEdge import *


class PrimitiveMutableLabelledMap(LabelledMap):
    """
    This class represent a more primitive version of MutableLabelledMap.
    It implements only the following operations for now ( addAfter, addBefore,deleteEdge,contractEdge)
    and place more responsability on the user and remove some of the other methods( such as checking if two demi edge are on the same node or face), 
    the only advantage is the fact that the operations listed are in O(1) instead of O(log(m)).
    """

    def __init__(
        self,
        sigma: Permutation = None,
        alpha: Permutation = None,
        adj=None,
        lmap=None
    ):

        if isinstance(lmap, LabelledMap):
            self.__init__(lmap.sigma, lmap.alpha)
            return

        super().__init__(sigma, alpha, adj)
        self.sigma = PrimitiveRotatingPermutation(self.sigma)
        self.alpha = PrimitiveRotatingPermutation(self.alpha)
        self.phi = PrimitiveRotatingPermutation(self.phi)
        self.phiUtilsAbstractor = PrimitiveRotatingPermutationUtilsAbstractor(
            self.phi)
        self.sigmaUtilsAbstractor = PrimitiveRotatingPermutationUtilsAbstractor(
            self.sigma)

        for e in range(1, self.q + 1):
            self.topologicalMap[e] = PrimitiveMutableTopologicalDemiEdge(
                self, e)

    def addEdge(self, startDemiEdge, endDemiEdge):
        """
        This will add an edge between the node of startDemiEdge to endDemiEdge(note that they
        need to be on the same node otherwise this will raise an error), the edge will be added as follow ,
        let denote (A,B) the demi edges composing this new edge A will be on the same node as startDemiEdge but before it
        and B on the same node as endDemiEdge but before it.
        It will return two MutableLabelledMap topoDemiEdgeA,topoDemiEdgeB corresponding to the new demi edge A and B

        WARNING: startDemiEdge and endDemiEdge must be on the same face otherwise there is no guarantee on the fact that self will stay stable
        ----
        Args:
            startDemiEdge: A demi edge of self
            endDemiEdge: A demi edge of self
            startDemiEdge,endDemiEdge must be on the same face or it will raise an error
        Returns:
            topoDemiEdgeA,topoDemiEdgeB as described above
        ----
        O(1)
        """

        # This will go before endDemiEdge
        newIndexEnd = self.q + 1

        # This will go before startDemiEdge
        newIndexStart = self.q + 2

        # We add the new edge to alpha
        self._addEdgeToAlpha()

        # We add the new edge to sigma
        # The order must be coherent with the choice made
        # newIndexStart,newIndexEnd
        self.sigma.addBefore(endDemiEdge)
        self.sigma.addBefore(startDemiEdge)

        # We add  the new demi edge to phi
        self.phi.cutAdd(startDemiEdge, endDemiEdge, newIndexStart, newIndexEnd)

        self.topologicalMap[newIndexStart] = PrimitiveMutableTopologicalDemiEdge(
            self, newIndexStart)
        self.topologicalMap[newIndexEnd] = PrimitiveMutableTopologicalDemiEdge(
            self, newIndexEnd)

        return self.X(newIndexStart), self.X(newIndexEnd)

    def _swapTopologicalDemiEdgeValue(self, demiEdge, otherDemiEdge):
        """
        This will change the index of the topologica demi edge associate to the demi edge of label demi edge to otherDemiEdge
        and same thing but swapping the role or demiEdge and otherDemiEdge
        ----
        O(1)
        """
        topoDemiEdge = self.X(demiEdge)
        otherTopoDemiEdge = self.X(otherDemiEdge)
        topoDemiEdge._swapIndex(otherTopoDemiEdge)
        self.topologicalMap[demiEdge] = otherTopoDemiEdge
        self.topologicalMap[otherDemiEdge] = topoDemiEdge

    def simpleSwap(self, demiEdge, otherDemiEdge):
        """
        This function relabel demiEdge into otherDemiEdge and otherDemiEdge into demiEdge
        -----
        Args:
            demiEdge,otherDemiEdge : two valid index
        -----
        O(1)
        """
        self._swapTopologicalDemiEdgeValue(demiEdge, otherDemiEdge)

        self.phi.swapIndex(demiEdge, otherDemiEdge)
        self.alpha.swapIndex(demiEdge, otherDemiEdge)
        self.sigma.swapIndex(demiEdge, otherDemiEdge)

    def labelToTheEnd(self, listIndexes):
        """
        This function relabel the indexes in listIndexes to take
        if there is say k indexes in listIndexes (and they are valid) it
        will relabel them into value in m,...,m-k+1
        -----
        listIndexes:
            A list of valid demi edge in self( otherwise it will raise an error)
        -----
        O(k) where k = len(listIndexes)
        """
        for i in listIndexes:
            if i != int(i) or i > self.q or i <= 0:
                raise ValueError("The list of indexes isn't correct")

        corres = self.phi.labelToTheEnd(listIndexes)
        self.alpha.labelToTheEnd(listIndexes)
        self.sigma.labelToTheEnd(listIndexes)
        for index in listIndexes:
            if index in corres:
                self._swapTopologicalDemiEdgeValue(index, corres[index])

    def _BruteDeleteEdge(self, demiEdge, sameFace):
        """
        This is an helper method it delete the demiEdge but don't check for
        connectivity before , if use in the wrong way it can break the connectivity invariant
        and thus many other method
        ----
        Args:
            demiEdge index corresponding to the edge to delete
            sameFace: a boolean indicating if demiEdge and alpha(demiEdge) are on the same face or not
        ----
        O(1)
        """
        otherDemiEdge = self.alpha(demiEdge)

        if demiEdge != self.q:
            self.labelToTheEnd([demiEdge, otherDemiEdge])
            self._BruteDeleteEdge(self.q, sameFace=sameFace)
            return

        # TopologicalDemiEdge thing
        topoDemiEdge = self.X(demiEdge)
        otherTopoDemiEdge = topoDemiEdge.c
        self._removeTopologicalDemiEdge(topoDemiEdge)
        self._removeTopologicalDemiEdge(otherTopoDemiEdge)

        # Actual removing
        self.sigma.deleteLastKIndex(2)
        self.alpha.deleteLastKIndex(2)
        if not sameFace:
            self.phi.mergeDelete(demiEdge, otherDemiEdge)
        else:
            # Because demiEdge doesn't break the connectivity
            # demiEdge and otherDemiEdge are on the same face
            # And the face will be cut
            # Look Strange  cause it mean that deleting a edge will increase the
            # Number of face but not so much because it won't increase the number of
            # On the current surface just the map obteined has more face on
            # The new lower genus surface
            # Note that this case only happen in genus > 0
            self.phi.cutDelete(demiEdge, otherDemiEdge)

    def deleteEdge(self, demiEdge, sameFace):
        """
        Args:
            demiEdge: a demi edge corresponding to the edge to delete
            sameFace: a boolean indicating if demiEdge and alpha(demiEdge) are on the same face or not
        ----
        O(1) 
        """
        self._BruteDeleteEdge(demiEdge, sameFace=sameFace)

    def _addEdgeToAlpha(self):
        """
        Add one edge compose of demi edges self.size() and self.size()+1 to alpha toward the end
        ----
        O(1)
        """
        self.alpha.stretch(2)
        self.alpha.addAfterGeneral(self.alpha.size(), self.alpha.size() - 1)

    def _addTopologicalDemiEdge(self, demiEdge):
        """
        Add a new MutableTopologicalDemiEdge to self with index associated to demiEdge
        ----
        O(1)
        """
        self.topologicalMap[demiEdge] = PrimitiveMutableTopologicalDemiEdge(
            self, demiEdge)

    def _removeTopologicalDemiEdge(self, topoDemiEdge):
        """
        Remove and make invalid topoDemiEdge
        ----
        O(1)
        """
        self.topologicalMap.pop(topoDemiEdge.raw)
        topoDemiEdge._invalidate()

    def addEdgeAfter(self, demiEdge):
        """
        This function will create an new edge attached to the same node as demi edge such that it is after
        demiEdge in the trigonometric order.
        ----
        Args:
            demiEdge
        ----
        O(1)
        """
        newDemiEdge = self.q + 1
        otherNewDemiEdge = self.q + 2

        # TopologicalDemiEdge things
        self._addTopologicalDemiEdge(newDemiEdge)
        self._addTopologicalDemiEdge(otherNewDemiEdge)

        # Updating sigma
        self.sigma.addAfter(demiEdge)
        self.sigma.stretch(1)

        # Updating alpha
        self._addEdgeToAlpha()

        # Updating phi
        self.phi.addAfter(self.alpha(demiEdge))
        self.phi.addAfter(newDemiEdge)

        return self.X(otherNewDemiEdge)

    def addEdgeBefore(self, demiEdge):
        """
        This function will create an new edge attached to the same node as demi edge such that it is before
        demiEdge in the trigonometric order.
        ----
        Args:
            demiEdge
        ----
        O(1)
        """

        return self.addEdgeAfter(self.sigma.inverseApply(demiEdge))

    def copy(self):
        """
        Returns: A copy of self
        -----
        O(m)
        """
        return PrimitiveMutableLabelledMap(sigma=self.sigma, alpha=self.alpha)

    def contractEdge(self, demiEdge):
        """
        Contract in self the edge corresponding to demiEdge.
        WARNING: For this version you demiEdge cannot be on a loop anyway because 
        contracting an loop is defined(in this library) as deleting an edge just
        call deleteEdge.
        -----
        Args:
            demiEdge the demiEdge to contract cannot be on a loop
        -----
        O(1)
        """
        if self.m == 1:
            raise ValueError(
                "Cannot contract an edge  in a map with only one edge")

        if demiEdge != self.q:
            self.labelToTheEnd([demiEdge, self.alpha(demiEdge)])
            self.contractEdge(self.q)
            return

        otherDemiEdge = self.alpha(demiEdge)

        self.sigma.mergeDelete(demiEdge, otherDemiEdge)
        self.phi.deleteLastKIndex(2)
        self.alpha.deleteLastKIndex(2)

    def areOnTheSameNode(self, demiEdgeA, demiEdgeB):
        """
        Not implemented for PrimitiveMutableLabelledMap
        """
        raise NotImplementedError(self)

    def areOnTheSameFace(self, demiEdgeA, demiEdgeB):
        """
        Not implemented for PrimitiveMutableLabelledMap
        """
        raise NotImplementedError(self)

    def numberInTheSameFace(self, demiEdge):
        """
        Not implemented for PrimitiveMutableLabelledMap
        """
        raise NotImplementedError(self)

    def numberInTheSameNode(self, demiEdge):
        """
        Not implemented for PrimitiveMutableLabelledMap
        """
        raise NotImplementedError(self)

    def checkTwoInTheSameFace(self, listDemiEdges):
        """
        Not implemented for PrimitiveMutableLabelledMap
        """
        raise NotImplementedError(self)

    def checkTwoInTheSameNode(self, listDemiEdges):
        """
        Not implemented for PrimitiveMutableLabelledMap
        """
        raise NotImplementedError(self)
