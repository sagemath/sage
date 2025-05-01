from sage.graphs.planar_maps.MapDecorator import CheckValid


@CheckValid
class TopologicalDemiEdge():
    def __init__(self, lmap, index):
        """
        This class is a an abstraction mean to represent the demi edge of a map in a more user
        friendly way than simple index, it is more link to the "topology" of the map than the raw index
        -----
        lmap: The map on which the self is bind
        index: The index associated to the demi edge it is bind
        -----
        O(1)
        """
        self._lmap = lmap
        self._index = index
        self._isValid = True

    @property
    def raw(self):
        """
        The raw index associated to demi edge on which self is bind to
        ----
        O(1)
        """
        return self.getIndex()

    @property
    def map(self):
        """
        The map on which self is bind to
        ----
        O(1)
        """
        return self.getMap()

    @property
    def c(self):
        """
        The next TopologicalDemiEdge  on the edge on which self is on
        ----
        O(1)
        """
        return self.nextOnEdge()

    @property
    def f(self):
        """
        The next TopologicalDemiEdge on the face on which self is on
        ----
        O(1)
        """
        return self.nextOnFace()

    @property
    def n(self):
        """
        The next TopologicalDemiEdge on the node on which self is on
        ----
        O(1)
        """
        return self.nextOnNode()

    @property
    def pf(self):
        """
        The previous TopologicalDemiEdge on the face on which self is on
        ----
        O(1)
        """
        return self.prevOnFace()

    @property
    def pn(self):
        """
        The previous TopologicalDemiEdge on the node on which self is on
        ----
        O(1)
        """
        return self.prevOnNode()

    @CheckValid
    def getIndex(self):
        """
        The raw index associated to demi edge on which self is bind to
        ----
        O(1)
        """
        return self._index

    @CheckValid
    def getMap(self):
        """
        Returns: The map on which self is bind
        ----
        O(1)
        """
        return self._lmap

    @CheckValid
    def nextOnEdge(self):
        """
        Returns: The next TopologicalDemiEdge  on the edge on which self is on
        ----
        O(1)
        """

        return self.map.getTopologicalDemiEdge(self.map.alpha(self.raw))

    @CheckValid
    def nextOnFace(self):
        """
        Returns: The next TopologicalDemiEdge on the face on which self is on
        ----
        O(1)
        """
        return self.map.getTopologicalDemiEdge(self.map.phi(self.raw))

    @CheckValid
    def nextOnNode(self):
        """
        Returns: The next TopologicalDemiEdge on the node on which self is on
        ----
        O(1)
        """
        return self.map.getTopologicalDemiEdge(self.map.sigma(self.raw))

    @CheckValid
    def prevOnFace(self):
        """
        Returns: The previous TopologicalDemiEdge on the face on which self is on
        ----
        O(1)
        """
        return self.map.getTopologicalDemiEdge(
            self.map.phi.inverseApply(self.raw))

    @CheckValid
    def prevOnNode(self):
        """
        Returns: The previous TopologicalDemiEdge on the node on which self is on
        ----
        O(1)
        """
        return self.map.getTopologicalDemiEdge(
            self.map.sigma.inverseApply(self.raw))

    @CheckValid
    def face(self):
        """
        Returns: A list containing all the TopologicalDemiEdge of demi edge on the
        same face as self
        ----
        O(f) where f is the number of demi edge on the face
        """
        lst = []
        for e in self.map.demiEdgesOnTheSameFace(self.raw):
            lst.append(self.map.getTopologicalDemiEdge(e))
        return lst

    @CheckValid
    def node(self):
        """
        Returns: A list containing all the TopologicalDemiEdge of demi edge on the
        same node as self
        ----
        O(n) where n is the number of demi edge on the node
        """

        lst = []
        for e in self.map.demiEdgesOnTheSameNode(self.raw):
            lst.append(self.map.getTopologicalDemiEdge(e))

        return lst

    @CheckValid
    def isOnSameFace(self, otherTopologicalDemiEdge):
        """
        Args:
            otherTopologicalDemiEdge a TopologicalDemiEdge on the same map as self
        Returns:
            A boolean indicating if they are on the same face
        ----
        O(1)
        """
        return self.map.areOnTheSameFace(
            self.raw, otherTopologicalDemiEdge.raw)

    @CheckValid
    def isOnSameNode(self, otherTopologicalDemiEdge):
        """
        Args:
            otherTopologicalDemiEdge a TopologicalDemiEdge on the same map as self
        Returns:
            A boolean indicating if they are on the same node
        ----
        O(1)
        """

        return self.map.areOnTheSameNode(
            self.raw, otherTopologicalDemiEdge.raw)

    def _invalidate(self):
        """
        Make self invalid
        -----
        O(1)
        """
        self._isValid = False

    def __repr__(self):
        return f"X({self.raw})"

    def _setIndex(self, newIndex):
        """
        Change the index of self
        ----
        Args: newIndex
        ----
        O(1)
        """
        self._index = newIndex

    def _swapIndex(self, otherTopoDemiEdge):
        """
        Swap indexes with otherTopoDemiEdge
        ----
        Args: otherTopoDemiEdge
        ----
        O(1)
        """
        tmp = self.raw
        self._setIndex(otherTopoDemiEdge.raw)
        otherTopoDemiEdge._setIndex(tmp)
