from numpy import delete
from sage.graphs.planar_maps.LabelledMap import *
from sage.all import Permutation  # Import sage library
from sage.graphs.planar_maps.MutableTopologicalDemiEdge import *
from sage.graphs.planar_maps.RotatingPermutationUtilsAbstractor import RotatingPermutationUtilsAbstractor
from sage.graphs.planar_maps.RotatingPermutation import RotatingPermutation


class MutableLabelledMap(LabelledMap):
    """
    This class represents a MutableLabelledMap , note that the indexes aren't fixed when changing the map so when want to keep
    an reference to a demi edge , get his MutableTopologicalDemiEdge with self.X(index), note that it  is guaranteed that when getting a MutableTopologicalDemiEdge from
    self that when calling the method public method of self  , the topological demi edge will always point to the correct label of the demi edge until the demi edge is deleted.
    .All the method returning map will return LabelledMap(not MutableLabelledMap) you will need to make them mutable
    by calling the constructor (e.g if your map is myMap call myMapMutable = MutableLabelledMap(lmap = myMap))
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
        self.sigma = RotatingPermutation(self.sigma)
        self.alpha = RotatingPermutation(self.alpha)
        self.phi = RotatingPermutation(self.phi)
        self.phiUtilsAbstractor = RotatingPermutationUtilsAbstractor(self.phi)
        self.sigmaUtilsAbstractor = RotatingPermutationUtilsAbstractor(
            self.sigma)

        for e in range(1, self.q + 1):
            self.topologicalMap[e] = MutableTopologicalDemiEdge(self, e)

    def willStillBeConnectedAfterEdgeDeletion(self, demiEdge):
        """
        This method return a boolean indicating if self will still be connected after deleting the edge corresponding to demiEdge
        If self of is a planar map it is efficient O(log(m)) otherwise it is O(m)
        -----
        Args:
            demiEdge: An index representing the demi edge corresponding to the edge we want to delete
        Returns:
            A boolean indicating whether or not self will still be connected if the edge is deleted
        -------
        O(log(m)) if self is a Planar map otherwise O(m)
        """
        # In case of higher genus than 0
        if self.genus() > 0:
            return self._willStillBeConnectedAfterEdgeDeletionHighGenus(
                demiEdge)

        if self.m == 1:
            return False
        otherHalf = self.alpha(demiEdge)

        # If they are not on the same face
        # Or one of them is attached to a node with degree 1
        # It means that it will still be connected after

        return (not self.areOnTheSameFace(demiEdge, otherHalf)) or self.sigma(
            otherHalf) == otherHalf or self.sigma(demiEdge) == demiEdge

    def willStillBeConnectedAfterNodeDeletion(self, demiEdge):
        """
        This method return a boolean indicating if self will still be connected after deleting the node corresponding to demiEdge
        If self of is a planar map it is efficient O(deg(node)*log(m)) otherwise it is O(m)
        ------
        Args:
            demiEdge: An index representing the demi edge corresponding to the node we want to delete
        Returns:
            A boolean indicating whether or not self will still be connected if the node is deleted
        -------
        O(deg(node)log(m)) where node is the node attached to demiEdge if self is planar map otherwise O(m+deg(node)log(m))
        """
        # In case of higher genus than 0
        if self.genus() > 0:
            return self._willStillBeConnectedAfterNodeDeletionHighGenus(
                demiEdge)

        curDemiEdge = demiEdge
        N = self.sigma.numberInCycle(demiEdge)

        # If the node has only one Edge
        if N == 1:
            return True

        lst = []

        # We're looking at each demiEdge on the node
        for _ in range(N):
            lst.append(curDemiEdge)

            # Here we see that  our node is link to a node with only one edge
            # We thus returned false
            if self.phi(curDemiEdge) == self.alpha(curDemiEdge):
                return False
            curDemiEdge = self.sigma(curDemiEdge)

        # WE MUST NOT have two demiEdges on the node on the same face
        # for our map to stay connected
        return not self.checkTwoInTheSameFace(lst)

    def _willStillBeConnectedAfterEdgeDeletionHighGenus(self, demiEdge):
        """
        A helper function use in case of high genus (>0) to check if after deleting the edge associated to demiEdge self
        will still be connected.
        -----
        O(m)
        """
        try:
            myCopy = self.copy()
            myCopy._BruteDeleteEdge(demiEdge)
            return transitiveCouplePermutation(myCopy.sigma, myCopy.alpha)
        except BaseException:
            return False

    def _willStillBeConnectedAfterNodeDeletionHighGenus(self, demiEdge):
        """
        A helper function use in case of high genus (>0) to check if after deleting the node associated to demiEdge self
        will still be connected.
        -----
        O(m+deg(node)log(m))
        """
        try:
            myCopy = self.copy()
            myCopy._BruteDeleteNode(demiEdge)
            return transitiveCouplePermutation(myCopy.sigma, myCopy.alpha)
        except BaseException:
            return False

    def addEdge(self, startDemiEdge, endDemiEdge):
        """
        This will add an edge between the node of startDemiEdge to endDemiEdge(note that they
        need to be on the same node otherwise this will raise an error), the edge will be added as follow ,
        let denote (A,B) the demi edges composing this new edge A will be on the same node as startDemiEdge but before it
        and B on the same node as endDemiEdge but before it.
        It will return two MutableLabelledMap topoDemiEdgeA,topoDemiEdgeB corresponding to the new demi edge A and B
        ----
        Args:
            startDemiEdge: A demi edge of self
            endDemiEdge: A demi edge of self
            startDemiEdge,endDemiEdge must be on the same face or it will raise an error
        Returns:
            topoDemiEdgeA,topoDemiEdgeB as described above
        """
        if not self.areOnTheSameFace(startDemiEdge, endDemiEdge):
            raise ValueError(
                f"{startDemiEdge} and {endDemiEdge} aren't on the same face")

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

        self.topologicalMap[newIndexStart] = MutableTopologicalDemiEdge(
            self, newIndexStart)
        self.topologicalMap[newIndexEnd] = MutableTopologicalDemiEdge(
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
        O(log(m))
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
        O(klog(m)) where k = len(listIndexes)
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

    def _BruteDeleteEdge(self, demiEdge):
        """
        This is an helper method it delete the demiEdge but don't check for
        connectivity before , if use in the wrong way it can break the connectivity invariant
        and thus many other method
        ----
        Args:
            demiEdge index corresponding to the edge to delete
        ----
        O(log(m))
        """
        otherDemiEdge = self.alpha(demiEdge)

        if demiEdge != self.q:
            self.labelToTheEnd([demiEdge, otherDemiEdge])
            self._BruteDeleteEdge(self.q)
            return

        # TopologicalDemiEdge thing
        topoDemiEdge = self.X(demiEdge)
        otherTopoDemiEdge = topoDemiEdge.c
        self._removeTopologicalDemiEdge(topoDemiEdge)
        self._removeTopologicalDemiEdge(otherTopoDemiEdge)

        # Actual removing
        self.sigma.deleteLastKIndex(2)
        self.alpha.deleteLastKIndex(2)
        if not self.areOnTheSameFace(demiEdge, otherDemiEdge):
            # We are merging the face
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

    def deleteEdge(self, demiEdge, trust=False):
        """
        This function delete demiEdge from self in case of planar map it is efficient O(log(m)) otherwise O(m) for higher genus
        map if trust = False. If trust = True as a parameter it is always of complexity O(log(m)), but it can break the invariant
        saying that the map is connected thus make all the other method unsafe by default trust = False.
        ----
        Args:
            demiEdge: a demi edge corresponding to the edge to delete
            trust: a parameter telling the function if it should trust the fact that the map will stay connected after
                deleting demiEdge default is False
        ----
        O(log(m)) if self is planar or trust = True otherwise O(m)
        """
        if not trust and not self.willStillBeConnectedAfterEdgeDeletion(
                demiEdge):
            raise ValueError(
                f"After deleting the edge associated to {demiEdge} the graph won't be connected anymore")

        self._BruteDeleteEdge(demiEdge)

    def _addEdgeToAlpha(self):
        """
        Add one edge compose of demi edges self.size() and self.size()+1 to alpha toward the end
        ----
        O(log(m))
        """
        self.alpha.stretch(2)
        self.alpha.addAfterGeneral(self.alpha.size(), self.alpha.size() - 1)

    def _addTopologicalDemiEdge(self, demiEdge):
        """
        Add a new MutableTopologicalDemiEdge to self with index associated to demiEdge
        ----
        O(1)
        """
        self.topologicalMap[demiEdge] = MutableTopologicalDemiEdge(
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
        O(log(m))
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
        O(log(m))
        """

        return self.addEdgeAfter(self.sigma.inverseApply(demiEdge))

    def deleteNode(self, demiEdge, trust=False):
        """
        This method will delete the node attached to demiEdge, when trust = False (default) it is efficient in case of planar map O(deg(node)*log(m))
        but for higher genus map it is of complexity O(m+deg(node)*log(m)). If trust = True as a parameter it is always of complexity O(log(m)), but it can break the invariant
        saying that the map is connected thus make all the other method unsafe by default trust = False.

        It will raise an error if the graph isn't connected after the operation and trust is set to False.
        ----
        Args:
            demiEdge : The demi edge on the node to delete
            trust: a parameter telling the method if it should trust the fact that the map will stay connected after
                deleting demiEdge default is False
        ----
        O(deg(node)*log(m)) if self is planar and O(m+deg(node)*log(m)) otherwise
        """
        if not trust and not self.willStillBeConnectedAfterNodeDeletion(
                demiEdge):
            raise ValueError(
                f"After deleting the node associated to {demiEdge} the graph won't be connected anymore")

        self._BruteDeleteNode(demiEdge)

    def _BruteDeleteNode(self, demiEdge):
        """
        This function will delete the node attached to self without checking if it will break
        the connectivity invariant thus it i dangerous if not used correctly
        ----
        Args:
            demiEdge : The demi edge on the node to delete
        ----
        O(deg(node)*log(m))
        """

        curDemiEdge = demiEdge
        N = self.numberInTheSameNode(curDemiEdge)
        while N > 0:

            self.simpleSwap(curDemiEdge, self.q)
            curDemiEdge = self.q
            self.simpleSwap(self.alpha(curDemiEdge), self.q - 1)
            otherDemiEdge = self.q - 1

            nxtDemiEdge = self.sigma(curDemiEdge)

            if self.areOnTheSameNode(otherDemiEdge, curDemiEdge):
                if nxtDemiEdge == otherDemiEdge:
                    nxtDemiEdge = self.sigma(otherDemiEdge)
                N -= 2
            else:
                N -= 1

            self._BruteDeleteEdge(curDemiEdge)
            curDemiEdge = nxtDemiEdge

    def contractFace(self, demiEdge):
        """
        Contract the face on which demiEdge is in
        ----
        Args:
            demiEdge
        ----
        O(tlog(m)) where t is the number of edge on the face containing demiEdge
        """
        if self.numberOfFaces() <= 2:
            raise ValueError("Cannot contract a face when there is <= 2 faces")

        curDemiEdge = demiEdge
        N = self.numberInTheSameFace(curDemiEdge)

        while N > 0:
            self.simpleSwap(curDemiEdge, self.q)
            curDemiEdge = self.q

            self.simpleSwap(self.alpha(curDemiEdge), self.q - 1)
            otherDemiEdge = self.q - 1

            nxtDemiEdge = self.phi(curDemiEdge)

            if self.areOnTheSameFace(otherDemiEdge, curDemiEdge):
                if nxtDemiEdge == otherDemiEdge:
                    nxtDemiEdge = self.phi(otherDemiEdge)
                N -= 2
            else:
                N -= 1

            self.contractEdge(curDemiEdge)

            curDemiEdge = nxtDemiEdge

    def copy(self):
        """
        Returns: A copy of self
        -----
        O(m)
        """
        return MutableLabelledMap(sigma=self.sigma, alpha=self.alpha)

    def contractEdge(self, demiEdge):
        """
        Contract in self the edge corresponding to demiEdge,demiEdge is on a loop edge it will just delete the edge
        -----
        O(log(m))
        """
        if self.m == 1:
            raise ValueError(
                "Cannot contract an edge  in a map with only one edge")

        # If this is a loop
        if self.areOnTheSameNode(demiEdge, self.alpha(demiEdge)):
            # We delete the edge
            # We can trust here cause a loop can't make the graph not connected
            self.deleteEdge(demiEdge, trust=True)
            return

        if demiEdge != self.q:
            self.labelToTheEnd([demiEdge, self.alpha(demiEdge)])
            self.contractEdge(self.q)
            return

        otherDemiEdge = self.alpha(demiEdge)

        # Merging the nodes here
        self.sigma.mergeDelete(demiEdge, otherDemiEdge)

        self.phi.deleteLastKIndex(2)
        self.alpha.deleteLastKIndex(2)

    def copyOnDemiEdge(self, demiEdge, otherMap, otherDemiEdge):
        """
        Given that demiEdge is such that it is attached to a node of degree one (otherwise this function will raise
        an error),this function will attach demiEdge before otherDemiEdge and then copy otherMap from this point on,
        I will return a MutableTopologicalDemiEdge corresponding to demiEdge in self but also a list of MutableTopologicalDemiEdge call it topoDemiEdgeList
        corresponding to the new MutableTopologicalDemiEdge of all the demi edge added from otherMap, such if a demi edge is of index i
        in otherMap than topoDemiEdgeList[i] is the MutableTopologicalDemiEdge corresponding to the copy of the demi edge inself.Note that this function
        won't modify otherMap(if otherMap isn't strictly the same object in memory as self).
        ----
        Args:
            demiEdge: A demi edge on a node of degree one in self
            otherMap: Another map it can be LabelledMap,RootedMap or MutableLabelledMap
            otherDemiEdge: A demi edge of otherMap
        ----
        O(p(log(m)+log(p))) where p = otherMap.m and m is the number of edge of self,
        note that it is much more efficient than O(p+m) mainly when m>>p
        """
        if not self.sigma(demiEdge) == demiEdge:
            raise ValueError(
                f"{demiEdge} must be on a node with degree one to call copyOnDemiEdge on it")

        C = self.q

        otherMap = otherMap.copy()

        # Alpha modification
        alphaCyclesToAdd = otherMap.alpha.to_cycles()

        for i in range(len(alphaCyclesToAdd)):
            a, b = alphaCyclesToAdd[i]
            alphaCyclesToAdd[i] = [a + C, b + C]
        self.alpha.addCycles(alphaCyclesToAdd)

        # Sigma modification
        sigmaCyclesToAdd = otherMap.sigma.to_cycles()

        for i in range(len(sigmaCyclesToAdd)):
            sigmaCyclesToAdd[i] = list(sigmaCyclesToAdd[i])
            for j in range(len(sigmaCyclesToAdd[i])):
                sigmaCyclesToAdd[i][j] += C

        self.sigma.addCycles(sigmaCyclesToAdd)

        self.sigma.addBeforeGeneral(otherDemiEdge + C, demiEdge)

        # Phi modification
        self.phi.stretch(otherMap.q)

        startDemiEdge = otherDemiEdge

        curDemiEdge = startDemiEdge

        self.phi.addAfterGeneral(self.alpha(demiEdge), startDemiEdge + C)

        while otherMap.phi(curDemiEdge) != startDemiEdge:
            self.phi.addAfterGeneral(
                curDemiEdge + C, otherMap.phi(curDemiEdge) + C)
            curDemiEdge = otherMap.phi(curDemiEdge)

        otherMapPhiCycles = otherMap.phi.to_cycles()
        phiCyclesToAdd = []

        for c in otherMapPhiCycles:
            if otherMap.areOnTheSameFace(otherDemiEdge, c[0]):
                continue
            phiCyclesToAdd.append(list(c))

        for i in range(len(phiCyclesToAdd)):
            for j in range(len(phiCyclesToAdd[i])):
                phiCyclesToAdd[i][j] += C
        self.phi.bruteAddCycles(phiCyclesToAdd)

        topoDemiEdgeList = []
        for i in range(1, otherMap.q + 1):
            self.topologicalMap[i +
                                C] = MutableTopologicalDemiEdge(self, i + C)
            topoDemiEdgeList.append(self.X(i + C))

        return self.X(demiEdge), topoDemiEdgeList

    def merge(self, demiEdge, otherMap, otherDemiEdge):
        """
        This will merge in self without modifying otherMap(if otherMap isn't the same object as self),
        what we mean by merging is the following draw an edge between the two map.
        It will return a couple (topoDemiEdge,topoDemiEdgeList) where topoDemiEdge is a MutableTopologicalDemiEdge associated
        to the edge that was added on the side of self, and topoDemiEdgeList which will contain a list of MutableTopologicalDemiEdge
        such that if i was a demi edge in  otherMap topoDemiEdgeList[i] is the MutableTopologicalDemiEdge associated  in self
        to the copy of i.
        ----
        Args:
            demiEdge: The demiEdge attached to the node on self on which to draw the edge
            otherMap: The other map (it can be a LabelledMap or RootedMap or MutableLabelledMap)
            otherDemiEdge: The otherDemiEdge attached to the node on otherMap on which to draw the new edge
        Returns:

            newTopoDemiEdge:
                A MutableTopologicalDemiEdge such that it corresponds to a new demi edge attached
                after demiEdge.

            topoDemiEdgeList:
                A list of MutableTopologicalDemiEdge corresponding to the description above
        ----
        O(p(log(m)+log(p))) where p = otherMap.m and m is the number of edge of self,
        note that it is much more efficient than O(p+m) mainly when m>>p
        """
        copyOtherMap = otherMap.copy()
        topoNewEdge = self.addEdgeBefore(demiEdge)
        return self.copyOnDemiEdge(
            topoNewEdge.raw, copyOtherMap, otherDemiEdge)

    def mergeNode(self, demiEdge, otherDemiEdge):
        """
        Merge the node attached to demiEdge and otherDemiEdge, they need to be on the same face
        ----
        Args:
            demiEdge,otherDemiEdge two demiEdge
        ----
        O(log(m)) where m is the number of edge of self
        """
        topoDemiEdge, _ = self.addEdge(demiEdge, otherDemiEdge)
        topoDemiEdge.contract()

    def areOnTheSameNode(self, demiEdgeA, demiEdgeB):
        """
        Args:
            demiEdgeA,demiEdgeB: Index representing demi edges of self
        Returns:
            A boolean indicating whether or note demiEdgeA and demiEdgeB are on the node
        -------
        O(log(m))
        """
        return self.sigmaUtilsAbstractor.sameCycle(demiEdgeA, demiEdgeB)

    def areOnTheSameFace(self, demiEdgeA, demiEdgeB):
        """
        Args:
            demiEdgeA,demiEdgeB: Index representing demi edges of self
        Returns:
            A boolean indicating whether or note demiEdgeA and demiEdgeB are on the same face
        -------
        O(log(m))
        """
        return self.phiUtilsAbstractor.sameCycle(demiEdgeA, demiEdgeB)

    def numberInTheSameFace(self, demiEdge):
        """
        Args:
            demiEdge
        Returns:
            The number of  demi edge on the same face as demi edge
        ----
        O(log(m))
        """
        return self.phiUtilsAbstractor.numberInCycle(demiEdge)

    def numberInTheSameNode(self, demiEdge):
        """
        Args:
            demiEdge
        Returns:
            The number of  demi edge on the same node as demi edge
        ----
        O(log(m))
        """

        return self.sigmaUtilsAbstractor.numberInCycle(demiEdge)

    def checkTwoInTheSameFace(self, listDemiEdges):
        """
        A method that will return a boolean indicating whether or not
        two demiEdge are on the same face
        -----
        Args:
            listDemiEdges: A list of demi edges
        OUTPUT:
            a boolean indicating whether or not there is two demi edge on the
            same face
        -----
        O(len(listDemiEdges)*log(m))
        """
        return self.phiUtilsAbstractor.checkTwoInTheSameCycle(listDemiEdges)

    def checkTwoInTheSameNode(self, listDemiEdges):
        """
        A method that will return a boolean indicating whether or not
        two demiEdge are on the same node
        -----
        Args:
            listDemiEdges: A list of demi edges
        OUTPUT:
            a boolean indicating whether or not there is two demi edge on the
            same node
        -----
        O(len(listDemiEdges)*log(m))
        """

        return self.sigmaUtilsAbstractor.checkTwoInTheSameCycle(listDemiEdges)
