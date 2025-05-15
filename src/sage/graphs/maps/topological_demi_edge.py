"""Define the TopologicalDemiEdge class, an abstraction meant to represent a demi-edge of a map in a more user-friendly way than raw indices."""

# from sage.graphs.maps.labelled_map import LabelledMap


class TopologicalDemiEdge():
    """This class is a an abstraction meant to represent the demi edge of a map in a more user-friendly way than simple indexes. It is more related to the "topological structure" of the map than the raw index."""

    _lmap: "LabelledMap"

    def __init__(self, lmap: "LabelledMap", index: int):
        r"""
        This class is a an abstraction meant to represent the demi edge of a map in a more user-friendly way
        than simple indexes. It is more related to the "topological structure" of the map than the raw index.

        INPUT:
        - ``lmap`` -- LabelledMap: The map to which the demi-edge is bind
        - ``index`` -- int: The index associated to the demi-edge

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3)
            X(3)

        NOTE:
            O(1)
        """
        self._lmap = lmap
        self._index = index
        self._isValid = True

    def _checkValid(self) -> None:
        """Raise an error is self isn't valid anymore."""
        if not self._isValid:
            raise ValueError("This TopologicalDemiEdge isn't valid anymore.")

    @property
    def raw(self) -> int:
        """
        The raw index associated to the demi-edge to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).raw
            3

        NOTE:
            O(1)
        """
        return self.getIndex()

    @property
    def map(self) -> "LabelledMap":
        """
        The map to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).map
            Labelled map | Sigma : [3, 2, 1, 4], Alpha : [2, 1, 4, 3]

        NOTE:
            O(1)
        """

        return self.getMap()

    @property
    def c(self) -> "TopologicalDemiEdge":
        """
        The other TopologicalDemiEdge on the edge on which self is on (equivalent to alpha(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).c
            X(4)

        NOTE:
            O(1)
        """
        return self.nextOnEdge()

    @property
    def f(self) -> "TopologicalDemiEdge":
        """
        The next TopologicalDemiEdge on the face on which self is on (equivalent to phi(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 2).f
            X(3)

        NOTE:
            O(1)
        """
        return self.nextOnFace()

    @property
    def n(self) -> "TopologicalDemiEdge":
        """
        The next TopologicalDemiEdge on the node on which self is on (equivalent to sigma(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).n
            X(1)

        NOTE:
            O(1)
        """
        return self.nextOnNode()

    @property
    def pf(self) -> "TopologicalDemiEdge":
        """
        The previous TopologicalDemiEdge on the face on which self is on (equivalent to phi^-1(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).pf
            X(2)

        NOTE:
            O(1)
        """
        return self.prevOnFace()

    @property
    def pn(self) -> "TopologicalDemiEdge":
        """
        The previous TopologicalDemiEdge on the node on which self is on (equivalent to sigma(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 1).n
            X(3)

        NOTE:
            O(1)
        """
        return self.prevOnNode()

    def getIndex(self) -> int:
        """
        Returns the raw index associated to the demi-edge to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).getIndex()
            3

        NOTE:
            O(1)
        """
        self._checkValid()
        return self._index

    def getMap(self) -> "LabelledMap":
        """
        Returns the map to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).getMap()
            Labelled map | Sigma : [3, 2, 1, 4], Alpha : [2, 1, 4, 3]

        NOTE:
            O(1)
        """
        self._checkValid()
        return self._lmap

    def nextOnEdge(self) -> "TopologicalDemiEdge":
        """
        Returns the other TopologicalDemiEdge on the edge on which self is on (equivalent to alpha(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).nextOnEdge()
            X(4)

        NOTE:
            O(1)
        """
        self._checkValid()

        return self.map.getTopologicalDemiEdge(self.map.alpha(self.raw))

    def nextOnFace(self) -> "TopologicalDemiEdge":
        """
        Return the next TopologicalDemiEdge on the face on which self is on (equivalent to phi(self)).

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 2).nextOnFace()
            X(3)

        NOTE:
            O(1)
        """
        self._checkValid()
        return self.map.getTopologicalDemiEdge(self.map.phi(self.raw))

    def nextOnNode(self) -> "TopologicalDemiEdge":
        """
        Return the next TopologicalDemiEdge on the node on which self is on (equivalent to sigma(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).nextOnNode()
            X(1)

        NOTE:
            O(1)
        """
        self._checkValid()
        return self.map.getTopologicalDemiEdge(self.map.sigma(self.raw))

    def prevOnFace(self) -> "TopologicalDemiEdge":
        """
        Return the previous TopologicalDemiEdge on the face on which self is on (equivalent to phi^-1(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 3).prevOnFace()
            X(2)

        NOTE:
            O(1)
        """
        self._checkValid()
        return self.map.getTopologicalDemiEdge(
            self.map.phi.inverseApply(self.raw))

    def prevOnNode(self) -> "TopologicalDemiEdge":
        """
        The previous TopologicalDemiEdge on the node on which self is on (equivalent to sigma(self))

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: TopologicalDemiEdge(lm, 1).prevOnNode()
            X(3)

        NOTE:
            O(1)
        """
        self._checkValid()
        return self.map.getTopologicalDemiEdge(
            self.map.sigma.inverseApply(self.raw))

    def face(self) -> list["TopologicalDemiEdge"]:
        """
        Return the list containing all the demi-edges on the same face as self.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: print (A.face())
            [X(1), X(3), X(2), X(5), X(4), X(7), X(11), X(16), X(18), X(13), X(12), X(15), X(14), X(19), X(20), X(17), X(9), X(8), X(10), X(6)]

        NOTE:
            O(k), where k is the number of demi-edges on the face
        """
        self._checkValid()
        lst = []
        for e in self.map.demiEdgesOnTheSameFace(self.raw):
            lst.append(self.map.getTopologicalDemiEdge(e))
        return lst

    def node(self) -> list["TopologicalDemiEdge"]:
        """
        Return the list containing all the demi-edges on the same node as self.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: print (A.node())
            [X(1), X(2), X(4)]

        NOTE:
            O(k), where k is the number of demi-edges on the node
        """
        self._checkValid()

        lst = []
        for e in self.map.demiEdgesOnTheSameNode(self.raw):
            lst.append(self.map.getTopologicalDemiEdge(e))

        return lst

    def isOnSameFace(self, otherTopologicalDemiEdge: "TopologicalDemiEdge") -> bool:
        """
        Return whether self and otherTopologicalDemiEdge are on the same face.

        INPUT:
        - ``otherTopologicalDemiEdge`` -- TopologicalDemiEdge: a demi-edge on the same map as self

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: lm = LabelledMap(sigma=sigma, alpha=alpha)
            sage: A = lm.X(1)
            sage: B = lm.X(9)
            sage: A.isOnSameFace(B)
            True

        NOTE:
            O(1)
        """
        self._checkValid()
        return bool(self.map.areOnTheSameFace(
            self.raw, otherTopologicalDemiEdge.raw))

    def isOnSameNode(self, otherTopologicalDemiEdge: "TopologicalDemiEdge") -> bool:
        """
        Return whether self and otherTopologicalDemiEdge are on the same node.

        INPUT:
        - ``otherTopologicalDemiEdge`` -- TopologicalDemiEdge: a demi-edge on the same map as self

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: lm = LabelledMap(sigma=sigma, alpha=alpha)
            sage: A = lm.X(1)
            sage: B = lm.X(9)
            sage: B.isOnSameNode(A)
            False

    .. NOTE:
            O(1)
        """
        self._checkValid()

        return bool(self.map.areOnTheSameNode(
            self.raw, otherTopologicalDemiEdge.raw))

    def _invalidate(self) -> None:
        """
        Make self invalid.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: TopologicalDemiEdge(None,5)._invalidate()

    .. NOTE:
            O(1)
        """
        self._isValid = False

    def __repr__(self) -> str:
        """
        Return a string representing self.

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: TopologicalDemiEdge(None,5)
            X(5)
        """
        return f"X({self.raw})"

    def _setIndex(self, newIndex: int) -> None:
        """
        Change the index of self

        INPUT:
        - ``newIndex`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: U  = TopologicalDemiEdge(None,5)
            sage: U
            X(5)
            sage: U._setIndex(22)
            sage: U
            X(22)

        NOTE:
            O(1)
        """
        self._index = newIndex

    def _swapIndex(self, otherTopologicalDemiEdge: "TopologicalDemiEdge") -> None:
        """
        Swap indexes with otherTopologicalDemiEdge.

        INPUT:
        - ``otherTopologicalDemiEdge`` -- TopologicalDemiEdge

        EXAMPLES::

            sage: from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
            sage: U  = TopologicalDemiEdge(None,5)
            sage: V  = TopologicalDemiEdge(None,10)
            sage: U,V
            (X(5), X(10))
            sage: U._swapIndex(V)
            sage: U,V
            (X(10), X(5))

        NOTE:
            O(1)
        """
        tmp = self.raw
        self._setIndex(otherTopologicalDemiEdge.raw)
        otherTopologicalDemiEdge._setIndex(tmp)
