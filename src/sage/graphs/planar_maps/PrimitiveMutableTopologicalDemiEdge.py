
from sage.graphs.planar_maps.TopologicalDemiEdge import TopologicalDemiEdge
from sage.graphs.planar_maps.MapDecorator import CheckValid
from sage.graphs.planar_maps.MapError import NotImplementedError


class PrimitiveMutableTopologicalDemiEdge(TopologicalDemiEdge):
    """
    This class implement a splecial version of TopologicalDemiEdge used in MutableLabelledMap
    their specificity is that they have method to modify the map for instance you can delete them
    ,add between them etc.
    """

    @CheckValid
    def delete(self, sameFace):
        """
        This function delete self from self.map         
        ----
        Args:
            sameFace: a indicating if self and self.c are on the same face or not
        WARNING: If this break the connectivity, no error will be raised and the map won't be stable anymore
        ----
        O(1)
        """

        self.map.deleteEdge(self.raw, sameFace)

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
        ----
        O(1)
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
        O(1)
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
        O(1)
        """

        return self.map.addEdgeBefore(self.raw)

    @CheckValid
    def contract(self):
        """
        Contract the edge bind to self , note that self cannot be on a loop 
        i.e self and self.c cannot be on the same node otherwise nothing is guaranted
        ----
        O(1)
        """
        self.map.contractEdge(self.raw)

    @CheckValid
    def isOnSameFace(self, otherTopologicalDemiEdge):
        """
        Not implemented for PrimitiveMutableTopologicalDemiEdge
        """
        raise NotImplementedError(self)

    @CheckValid
    def isOnSameNode(self, otherTopologicalDemiEdge):
        """
        Not implemented for PrimitiveMutableTopologicalDemiEdge
        """
        raise NotImplementedError(self)
