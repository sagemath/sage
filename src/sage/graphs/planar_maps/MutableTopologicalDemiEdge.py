from sage.graphs.planar_maps.TopologicalDemiEdge import TopologicalDemiEdge
from sage.graphs.planar_maps.MapDecorator import CheckValid


class MutableTopologicalDemiEdge(TopologicalDemiEdge):
    """
    This class implement a splecial version of TopologicalDemiEdge used in MutableLabelledMap
    their specificity is that they have method to modify the map for instance you can delete them
    ,add between them etc.
    """

    @CheckValid
    def delete(self, trust=False):
        """
        This function delete self from self.map in case of planar map it is efficient O(log(m)) otherwise O(m) for higher genus
        map if trust = False. If trust = True as a parameter it is always of complexity O(log(m)), but it can break the invariant
        saying that the self.map is connected thus make all the other method unsafe by default trust = False.
        ----
        Args:
            trust: a parameter telling the function if it should trust the fact that the map will stay connected after
                deleting demiEdge default is False
        ----
        O(log(m)) if self is planar or trust = True otherwise O(m)
        """

        self.map.deleteEdge(self.raw, trust=trust)

    @CheckValid
    def link(self, otherTopoDemiEdge):
        """
        This will add an edge between the node of self to otherTopoDemiEdge(note that they
        need to be on the same node otherwise this will raise an error), the edge will be added as follow ,
        let denote (A,B) the demi edges composing this new edge A will be on the same node as self but before it
        and B on the same node as otherTopoDemiEdge but before it.
        It will return two MutableLabelledMap topoDemiEdgeA,topoDemiEdgeB corresponding to the new demi edge A and B
        ----
        Args:
            otherTopoDemiEdge: Another TopologicalDemiEdge on the same facee as self
        Returns:
            topoDemiEdgeA,topoDemiEdgeB as described above
        """

        return self.map.addEdge(self.raw, otherTopoDemiEdge.raw)

    @CheckValid
    def addEdgeAfter(self):
        """
        This method will create a new edge, such that it is added on the same node as self but after it in the
        trigonometric order
        -----
        Returns: The TopologicalDemiEdge associate to demi edge of the new edge which is on the newly added node
        -----
        O(log(m))
        """
        return self.map.addEdgeAfter(self.raw)

    @CheckValid
    def addEdgeBefore(self):
        """
        This method will create a new edge, such that it is added on the same node as self but before it
        in the trigonometric order
        -----
        Returns: The TopologicalDemiEdge associate to demi edge of the new edge which is on the newly added node
        -----
        O(log(m))
        """

        return self.map.addEdgeBefore(self.raw)

    @CheckValid
    def deleteNode(self, trust=False):
        """
        This method will delete the node attached to self, when trust = False (default) it is efficient in case of planar map O(deg(node)*log(m))
        but for higher genus map it is of complexity O(m+deg(node)*log(m)) where m is the number of edge in self.map.
        If trust = True as a parameter it is always of complexity O(log(m)), but it can break the invariant saying that
        the self.map is connected thus make all the other method unsafe by default trust = False.
        It will raise an error if the self.map isn't connected after the operation and trust is set to False.
        ----
        Args:
            trust: a parameter telling the method if it should trust the fact that the self.map will stay connected after
                deleting self default is False
        ----
        O(deg(node)*log(m)) if self.map is planar and O(m+deg(node)*log(m)) otherwise
        """

        self.map.deleteNode(self.raw, trust=trust)

    @CheckValid
    def contract(self):
        """
        Contract the edge bind to self
        ----
        O(log(m))
        """
        self.map.contractEdge(self.raw)

    @CheckValid
    def contractFace(self):
        """
        Contract the face on which self is in
        ----
        O(tlog(m)) where t is the number of edge on the face containing self
        """
        self.map.contractFace(self.raw)

    @CheckValid
    def mergeMap(self, otherTopoDemiEdge):
        """
        This will merge in self.map without modifying otherTopoDemiEdge.map(if otherTopoDemiEdge.map isn't the same object as self.map),
        what we mean by merging is the following draw an edge between the two map.
        It will return a couple (topoDemiEdge,topoDemiEdgeList) where topoDemiEdge is a MutableTopologicalDemiEdge associated
        to the edge that was added on the side of self.map, and topoDemiEdgeList which will contain a list of MutableTopologicalDemiEdge
        such that if i was a demi edge in otherDemiEdge.map , topoDemiEdgeList[i] is the MutableTopologicalDemiEdge associated  in self.map
        to the copy of i.
        ----
        Args:
            otherTopoDemiEdge: The otherTopoDemiEdge attached to the node on otherTopoDemiEdge.map which to draw the new edge
        Returns:

            newTopoDemiEdge:
                A MutableTopologicalDemiEdge such that it corresponds to a new demi edge attached
                after demiEdge.

            topoDemiEdgeList:
                A list of MutableTopologicalDemiEdge corresponding to the description above
        ----
        O(p(log(m)+log(p))) where p = otherTopoDemiEdge.map.m and m is the number of edge of self.map,
        note that it is much more efficient than O(p+m) mainly when m>>p
        """

        return self.map.merge(
            self.raw, otherTopoDemiEdge.map, otherTopoDemiEdge.raw)

    @CheckValid
    def copyOn(self, otherTopoDemiEdge):
        """
        Given that self is such that it is attached to a node of degree one (otherwise this function will raise
        an error),this function will attach self before otherTopoDemiEdge and then copy otherTopoDemiEdge.map from this point on,
        I will return a MutableTopologicalDemiEdge corresponding to demiEdge in self.map but also a list of MutableTopologicalDemiEdge call it topoDemiEdgeList
        corresponding to the new MutableTopologicalDemiEdge of all the demi edge added from otherTopoDemiEdge.map, such if a demi edge is of index i
        in otherTopoDemiEdge.map than topoDemiEdgeList[i] is the MutableTopologicalDemiEdge corresponding to the copy of the demi edge in self.map.Note that this function
        won't modify otherTopoDemiEdge.map(if otherTopoDemiEdge.map isn't strictly the same object in memory as self.map).
        ----
        Args:
            otherTopoDemiEdge: Another TopologicalDemiEdge
        ----
        O(p(log(m)+log(p))) where p = otherTopoDemiEdge.map.m and m is the number of edge of self.map,
        note that it is much more efficient than O(p+m) mainly when m>>p
        """
        if self.map.numberInTheSameNode(self.raw) > 1:
            raise ValueError(
                "Self isn't attached to a node of degree one cannot copyOn")
        return self.map.copyOnDemiEdge(
            self.raw, otherTopoDemiEdge.map, otherTopoDemiEdge.raw)

    @CheckValid
    def mergeNode(self, otherTopoDemiEdge):
        """
        Merge the node attached to self and otherTopoDemiEdge, they need to be on the same face
        ----
        Args:
            otherTopoDemiEdge: the other TopologicalDemiEdge on the same face as sel
        ----
        O(log(m)) where m is the number of edge of self.map
        """

        if self.map != otherTopoDemiEdge.map:
            raise ValueError(
                "Cannot mergeNode between two demi edge on different map")
        self.map.mergeNode(self.raw, otherTopoDemiEdge.raw)

    @CheckValid
    def isOnSameFace(self, otherTopologicalDemiEdge):
        """
        Args:
            otherTopologicalDemiEdge a TopologicalDemiEdge on the same map as self
        Returns:
            A boolean indicating if they are on the same face
        ----
        O(log(m)) where m is the number of edge of self.map
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
        O(log(m)) where m is the number of edge of self.map
        """

        return self.map.areOnTheSameNode(
            self.raw, otherTopologicalDemiEdge.raw)
