from sage.all import Permutation  # Import sage library
from sage.graphs.maps.labelled_map import *
from sage.graphs.maps.mutable_topological_demi_edge import *
from sage.graphs.maps.rotating_permutation_utils_abstractor import RotatingPermutationUtilsAbstractor
from sage.graphs.maps.rotating_permutation import RotatingPermutation


class MutableLabelledMap(LabelledMap):
    """
    This class represents a MutableLabelledMap , note that the indexes aren't fixed when changing the map so when want to keep
    an reference to a demi edge , get his MutableTopologicalDemiEdge with self.X(index), note that it  is guaranteed that when getting a MutableTopologicalDemiEdge from
    self that when calling the method public method of self  , the topological demi edge will always point to the correct label of the demi edge until the demi edge is deleted.
    .All the method returning map will return LabelledMap(not MutableLabelledMap) you will need to make them mutable
    by calling the constructor (e.g if your map is myMap call myMapMutable = MutableLabelledMap(lmap = myMap))

    Attributes
    ----------
    sigma :  Permutation or MapPermutation
        Permutation that maps a half-edge to the half-edge incident to
        it in a anti-clockwise direction around the vertex it belongs to.

    alpha : Permutation or MapPermutation
        Fixed-point free involution whose cycles are given by the edges.

    phi: Permutation or MapPermutation
        Permutation that maps a half-edges to the half-edge next to it
        in his face in the clockwise orde.

    m: The number of edge of the map

    g: The genus of the map

    f: The number of face of the map

    q: The number of demi edges of the map

    """

    def __init__(
        self,
        sigma: Permutation = None,
        alpha: Permutation = None,
        adj=None,
        lmap=None
    ):
        r"""
        Init the MutableLabelledMap

        INPUT:

        - ``sigma`` -- Permutation | MapPermutation | None ; Permutation that maps a half-edge
          to the half-edge incident to it in anti-clockwise direction around
          the vertex it belongs to.
        - ``alpha`` -- Permutation | MapPermutation | None ;Fixed-point free involution whose cycles are given by the edges.
        - ``ajd``-- List[Tuples] | None ; and adjacency list be careful the order of the
            node in your adjaceny will be used to choose the embedding
        - ``lmap`` -- LabelledMap | None

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)

        NOTE:
            O(mlog(m)) where m is the size of the map
        """
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

        INPUT:
        - ``demiEdge`` -- int; An index representing the demi edge corresponding to the edge we want to delete

        OUTPUT:
            A boolean indicating whether or not self will still be connected if the edge is deleted

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.willStillBeConnectedAfterEdgeDeletion(2)
            False

        NOTE:
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
        # It means that it will still be connected after

        return (not self.areOnTheSameFace(demiEdge, otherHalf))

    def willStillBeConnectedAfterNodeDeletion(self, demiEdge):
        """
        This method return a boolean indicating if self will still be connected after deleting the node corresponding to demiEdge
        If self of is a planar map it is efficient O(deg(node)*log(m)) otherwise it is O(m)

        INPUT:
        - ``demiEdge`` -- int ;  An index representing the demi edge corresponding to the node we want to delete

        OUTPUT:
            A boolean indicating whether or not self will still be connected if the node is deleted

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.willStillBeConnectedAfterNodeDeletion(2)
            False

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int

        OUTPUT:
            A boolean indicating if self will still be connected after demiEdge deletion

        EXAMPLES::

            sage: alphaH = Permutation([(1,2),(3,4),(5,6)])
            sage: sigmaH = Permutation([(1,4,5,3),(6,2)])
            sage: mmH = MutableLabelledMap(sigma=sigmaH,alpha=alphaH)
            sage: mmH.g
            1
            sage: mmH._willStillBeConnectedAfterEdgeDeletionHighGenus(1)
            True

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int

        OUTPUT:
            A boolean indicating if self will still be connected after the deletion of the node on which demiEdge
            is attached

        EXAMPLES::

            sage: alphaH = Permutation([(1,2),(3,4),(5,6)])
            sage: sigmaH = Permutation([(1,4,5,3),(6,2)])
            sage: mmH = MutableLabelledMap(sigma=sigmaH,alpha=alphaH)
            sage: mmH.g
            1
            sage: mmH._willStillBeConnectedAfterNodeDeletionHighGenus(1)
            False

        NOTE:
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
        and B on the same node as endDemiEdge but before it.It will return two MutableLabelledMap topoDemiEdgeA,topoDemiEdgeB
        corresponding to the new demi edge A and B

        INPUT:
        - ``startDemiEdge`` -- int ;  A demi edge of self
        - ``endDemiEdge`` -- int ; A demi edge of self
        - ``startDemiEdge`` -- int
        - ``endDemiEdge`` -- int ;must be on the same face as startDemiEdge or it will raise an error

        OUTPUT:
            topoDemiEdgeA,topoDemiEdgeB as described above

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.addEdge(3,7)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 22, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 4, 21, 3)]

        NOTE:
            O(log(m))
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

        INPUT:
        - ``demiEdge`` -- int
        - ``,otherDemiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(1)
            sage: B = mm.X(2)
            sage: A,B
            (X(1), X(2))
            sage: mm._swapTopologicalDemiEdgeValue(1,2)
            sage: A,B
            (X(2), X(1))

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int
        - ``otherDemiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(1)
            sage: B = mm.X(2)
            sage: A,B
            (X(1), X(2))
            sage: mm.simpleSwap(2,1)
            sage: A,B
            (X(2), X(1))

        NOTE:
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

        INPUT:
        - ``listIndexes`` -- List[int]; A list of valid demi edge in self( otherwise it will raise an error)

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.labelToTheEnd([7,5,3])
            sage: mm.faces()
            [(1, 18, 2, 19, 4, 20, 11, 16, 3, 13, 12, 15, 14, 5, 7, 17, 9, 8, 10, 6)]

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int ; demi-edge on the edge to delete


        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.addEdge(3,11)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 22, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 4, 7, 21, 3)]
            sage: mm.deleteEdge(22)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]

        NOTE:
            O(log(m))
        """
        otherDemiEdge = self.alpha(demiEdge)

        if demiEdge != self.q or otherDemiEdge != self.q-1:
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

        INPUT:
        - ``demiEdge`` -- int ; a demi edge corresponding on the edge to delete
        - ``trust`` -- boolean ; a parameter telling the function if it should trust the fact that the map will stay connected after
        deleting demiEdge default is False

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.addEdge(3,11)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 22, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 4, 7, 21, 3)]
            sage: mm.deleteEdge(22)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]

        NOTE:
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

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.alpha.pretty_print()
            Rotating permutation: [(1, 3), (2, 5), (4, 6), (7, 9), (8, 10), (11, 13), (12, 15), (14, 17), (16, 18), (19, 20)]
            sage: mm._addEdgeToAlpha()
            sage: mm.alpha.pretty_print()
            Rotating permutation: [(1, 3), (2, 5), (4, 6), (7, 9), (8, 10), (11, 13), (12, 15), (14, 17), (16, 18), (19, 20), (21, 22)]

        NOTE:
            O(log(m))
        """
        self.alpha.stretch(2)
        self.alpha.addAfterGeneral(self.alpha.size(), self.alpha.size() - 1)

    def _addTopologicalDemiEdge(self, demiEdge):
        """
        Add a new MutableTopologicalDemiEdge to self with index associated to demiEdge

        INPUT:
        - ``demiEdge`` -- int


        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.X(42)
            ....: except:
            ....:    print("NOT TOPO")
            ....:
            NOT TOPO
            sage: mm._addTopologicalDemiEdge(42)
            sage: mm.X(42)
            X(42)

        ..NOTE:
            O(1)
        """
        self.topologicalMap[demiEdge] = MutableTopologicalDemiEdge(
            self, demiEdge)

    def _removeTopologicalDemiEdge(self, topoDemiEdge):
        """
        Remove and make invalid topoDemiEdge

        INPUT:
        - ``topoDemiEdge`` -- TopologicalDemiEdge


        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(10)
            sage: mm._removeTopologicalDemiEdge(A)
            sage: try:
            ....:    mm.X(10)
            ....: except:
            ....:    print("NOT TOPO")
            ....:
            NOT TOPO

        NOTE:
            O(1)
        """
        self.topologicalMap.pop(topoDemiEdge.raw)
        topoDemiEdge._invalidate()

    def addEdgeAfter(self, demiEdge):
        """
        This function will create an new edge attached to the same node as demi edge such that it is after
        demiEdge in the trigonometric order.

        INPUT:
        - ``demiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(10)
            sage: A.n,A.pn
            (X(10), X(10))
            sage: mm.addEdgeAfter(10)
            X(22)
            sage: mm.addEdgeBefore(10)
            X(24)
            sage: A.n,A.pn
            (X(21), X(23))

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(10)
            sage: A.n,A.pn
            (X(10), X(10))
            sage: mm.addEdgeAfter(10)
            X(22)
            sage: mm.addEdgeBefore(10)
            X(24)
            sage: A.n,A.pn
            (X(21), X(23))

        NOTE:
            O(log(m))
        """

        return self.addEdgeAfter(self.sigma.inverseApply(demiEdge))

    def deleteNode(self, demiEdge, trust=False):
        """
        This method will delete the node attached to demiEdge, when trust = False (default) it is efficient in case of planar map O(deg(node)*log(m))
        but for higher genus map it is of complexity O(m+deg(node)*log(m)). If trust = True as a parameter it is always of complexity O(log(m)), but it can break the invariant
        saying that the map is connected thus make all the other method unsafe by default trust = False.

        It will raise an error if the graph isn't connected after the operation and trust is set to False.

        INPUT:
        - ``demiEdge`` -- int ;  The demi edge on the node to delete
        - ``trust`` -- bool ; a parameter telling the method if it should trust the fact that the map will stay connected after
        deleting demiEdge default is False

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.nodes()
            [(1, 2, 4),
             (3,),
             (5,),
             (6, 7, 8),
             (9, 11, 12, 14),
             (10,),
             (13, 16),
             (15,),
             (17, 19),
             (18,),
             (20,)]
            sage: mm.deleteNode(3)
            sage: mm.nodes()
            [(1, 17),
             (2, 4),
             (3,),
             (5,),
             (6, 7, 8),
             (9, 11, 12, 14),
             (10,),
             (13, 16),
             (15,),
             (18,)]

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int ; The demi edge on the node to delete

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.nodes()
            [(1, 2, 4),
             (3,),
             (5,),
             (6, 7, 8),
             (9, 11, 12, 14),
             (10,),
             (13, 16),
             (15,),
             (17, 19),
             (18,),
             (20,)]
            sage: mm._BruteDeleteNode(3)
            sage: mm.nodes()
            [(1, 17),
             (2, 4),
             (3,),
             (5,),
             (6, 7, 8),
             (9, 11, 12, 14),
             (10,),
             (13, 16),
             (15,),
             (18,)]


        NOTE:
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
        INPUT:
        - ``demiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.addEdge(3,4)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 22, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 21, 3)]
            sage: mm.addEdge(4,7)
            (X(24), X(23))
            sage: mm.faces()
            [(1, 22, 24, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 21, 3),
             (4, 23)]
            sage: mm.contractFace(5)
            sage: mm.faces()
            [(1, 3, 17, 9, 8, 10, 6, 5, 7, 11, 16, 18, 13, 12, 15, 14), (2, 4)]

        NOTE:
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
        OUTPUT:
        Returns A copy of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.copy() == mm
            True

        NOTE:
            O(m)
        """
        return MutableLabelledMap(sigma=self.sigma, alpha=self.alpha)

    def contractEdge(self, demiEdge):
        """
        Contract in self the edge corresponding to demiEdge,demiEdge is on a loop edge it will just delete the edge

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.contractEdge(3)
            sage: mm.faces()
            [(1, 3, 17, 9, 8, 10, 6, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14)]

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int ; A demi edge on a node of degree one in self
        - ``otherMap`` -- LabelledMap; Another map it can be LabelledMap,RootedMap or MutableLabelledMap
        otherDemiEdge: A demi edge of otherMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.n
            11
            sage: lst = mm.copyOnDemiEdge(3,mm,4)
            sage: mm.n
            21

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int ; The demiEdge attached to the node on self on which to draw the edge
        - ``otherMap`` -- LabelledMap: The other map (it can be a LabelledMap or RootedMap or MutableLabelledMap)
        - ``otherDemiEdge`` -- int ; The otherDemiEdge attached to the node on otherMap on which to draw the new edge

        OUTPUT:

            newTopoDemiEdge:A MutableTopologicalDemiEdge such that it corresponds to a new demi edge attached
            after demiEdge.

            topoDemiEdgeList: A list of MutableTopologicalDemiEdge corresponding to the description above

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.n
            11
            sage: lst = mm.copyOnDemiEdge(3,mm,4)
            sage: mm.n
            21

        NOTE:
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

        INPUT:
        - ``demiEdge`` -- int
        - ``otherDemiEdge`` -- int


        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.mergeNode(1,11)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7), (6, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10)]

        NOTE:
            O(log(m)) where m is the number of edge of self
        """
        topoDemiEdge, _ = self.addEdge(demiEdge, otherDemiEdge)
        topoDemiEdge.contract()

    def areOnTheSameNode(self, demiEdgeA, demiEdgeB):
        """
        INPUT:
        - ``demiEdgeA`` -- int
        - ``demiEdgeB``-- int

        OUTPUT:
            A boolean indicating whether or note demiEdgeA and demiEdgeB are on the node

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.numberInTheSameNode(1)
            3
            sage: mm.numberInTheSameFace(1)
            20
            sage: mm.areOnTheSameNode(1,7)
            False
            sage: mm.areOnTheSameFace(1,7)
            True
            sage: mm.checkTwoInTheSameNode([1,7,8,9])
            True
            sage: mm.checkTwoInTheSameFace([1,7,8,9])
            True

        NOTE:
            O(log(m))
        """
        return self.sigmaUtilsAbstractor.sameCycle(demiEdgeA, demiEdgeB)

    def areOnTheSameFace(self, demiEdgeA, demiEdgeB):
        """
        INPUT:
        - ``demiEdgeA`` -- int
        - ``demiEdgeB``-- int

        OUTPUT:
            A boolean indicating whether or note demiEdgeA and demiEdgeB are on the same face


        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.numberInTheSameNode(1)
            3
            sage: mm.numberInTheSameFace(1)
            20
            sage: mm.areOnTheSameNode(1,7)
            False
            sage: mm.areOnTheSameFace(1,7)
            True
            sage: mm.checkTwoInTheSameNode([1,7,8,9])
            True
            sage: mm.checkTwoInTheSameFace([1,7,8,9])
            True

        NOTE:
            O(log(m))
        """
        return self.phiUtilsAbstractor.sameCycle(demiEdgeA, demiEdgeB)

    def numberInTheSameFace(self, demiEdge):
        """
        INPUT:
        - ``demiEdge`` -- int

        OUTPUT:
            The number of  demi edge on the same face as demi edge

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.numberInTheSameNode(1)
            3
            sage: mm.numberInTheSameFace(1)
            20
            sage: mm.areOnTheSameNode(1,7)
            False
            sage: mm.areOnTheSameFace(1,7)
            True
            sage: mm.checkTwoInTheSameNode([1,7,8,9])
            True
            sage: mm.checkTwoInTheSameFace([1,7,8,9])
            True

        NOTE:
            O(log(m))
        """
        return self.phiUtilsAbstractor.numberInCycle(demiEdge)

    def numberInTheSameNode(self, demiEdge):
        """
        INPUT:
        - ``demiEdge`` -- int

        OUTPUT:
            The number of  demi edge on the same node as demi edge

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.numberInTheSameNode(1)
            3
            sage: mm.numberInTheSameFace(1)
            20
            sage: mm.areOnTheSameNode(1,7)
            False
            sage: mm.areOnTheSameFace(1,7)
            True
            sage: mm.checkTwoInTheSameNode([1,7,8,9])
            True
            sage: mm.checkTwoInTheSameFace([1,7,8,9])
            True

        NOTE:
            O(log(m))
        """

        return self.sigmaUtilsAbstractor.numberInCycle(demiEdge)

    def checkTwoInTheSameFace(self, listDemiEdges):
        """
        A method that will return a boolean indicating whether or not
        two demiEdge are on the same face
        INPUT:
        - ``listDemiEdges`` -- List[int]; A list of demi edges

        OUTPUT:
            a boolean indicating whether or not there is two demi edge on the
            same face

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.numberInTheSameNode(1)
            3
            sage: mm.numberInTheSameFace(1)
            20
            sage: mm.areOnTheSameNode(1,7)
            False
            sage: mm.areOnTheSameFace(1,7)
            True
            sage: mm.checkTwoInTheSameNode([1,7,8,9])
            True
            sage: mm.checkTwoInTheSameFace([1,7,8,9])
            True

        NOTE:
            O(len(listDemiEdges)*log(m))
        """
        return self.phiUtilsAbstractor.checkTwoInTheSameCycle(listDemiEdges)

    def checkTwoInTheSameNode(self, listDemiEdges):
        """
        A method that will return a boolean indicating whether or not
        two demiEdge are on the same node

        INPUT:
        - ``listDemiEdges`` -- List[int]; A list of demi edges

        OUTPUT:
            a boolean indicating whether or not there is two demi edge on the
            same node

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.numberInTheSameNode(1)
            3
            sage: mm.numberInTheSameFace(1)
            20
            sage: mm.areOnTheSameNode(1,7)
            False
            sage: mm.areOnTheSameFace(1,7)
            True
            sage: mm.checkTwoInTheSameNode([1,7,8,9])
            True
            sage: mm.checkTwoInTheSameFace([1,7,8,9])
            True

        NOTE:
            O(len(listDemiEdges)*log(m))
        """

        return self.sigmaUtilsAbstractor.checkTwoInTheSameCycle(listDemiEdges)
