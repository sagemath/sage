
from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
from sage.graphs.maps.map_decorator import CheckValid
from sage.graphs.maps.map_error import NotImplemented


class PrimitiveMutableTopologicalDemiEdge(TopologicalDemiEdge):
    """
    This class implement a special version of TopologicalDemiEdge used in PrimitiveMutableLabelledMap
    their specificity is that they have method to modify the map for instance you can delete them
    ,add between them etc.
    """

    @CheckValid
    def delete(self, sameFace):
        """
        This function delete self from self.map         

        INPUT:
            sameFace: a indicating if self and self.c are on the same face or not


        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(sigma=sigma,alpha=alpha)
            sage: cp = mm.copy()
            sage: A = mm.X(5)
            sage: B = mm.X(13)
            sage: U,V = A.link(B)
            sage: cp == mm
            False
            sage: U.delete(sameFace=False)
            sage: cp == mm
            True

        .. NOTE::
            O(1),If this break the connectivity, no error will be raised and the map won't be stable anymore
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

        INPUT:
            otherTopoDemiEdge: Another TopologicalDemiEdge on the same facee as self

        OUTPUT:
            topoDemiEdgeA,topoDemiEdgeB as described above

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(sigma=sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: B = mm.X(13)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: A.link(B)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 3, 2, 22, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (4, 7, 11, 16, 18, 21, 5)]

        .. NOTE::
            O(1)
        """

        return self.map.addEdge(self.raw, otherTopoDemiEdge.raw)

    @CheckValid
    def addEdgeAfter(self):
        """
        This method will create a new edge, such that it is added on the same node as self but after it in the
        trigonometric order

        OUTPUT: 
            The TopologicalDemiEdge associate to demi edge of the new edge which is on the newly added node

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(sigma=sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: A.addEdgeAfter()
            X(22)
            sage: A.addEdgeBefore()
            X(24)

        .. NOTE::
            O(1)
        """
        return self.map.addEdgeAfter(self.raw)

    @CheckValid
    def addEdgeBefore(self):
        """
        This method will create a new edge, such that it is added on the same node as self but before it
        in the trigonometric order

        OUTPUT: 
            The TopologicalDemiEdge associate to demi edge of the new edge which is on the newly added node

        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(sigma=sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: A.addEdgeAfter()
            X(22)
            sage: A.addEdgeBefore()
            X(24)

        .. NOTE::
            O(1)
        """

        return self.map.addEdgeBefore(self.raw)

    @CheckValid
    def contract(self):
        """
        Contract the edge bind to self , note that self cannot be on a loop 
        i.e self and self.c cannot be on the same node otherwise nothing is guaranted

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(sigma=sigma,alpha=alpha)
            sage: A = mm.X(5)
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
            sage: A.contract()
            sage: mm.nodes()
            [(1, 4),
             (2, 17),
             (3,),
             (5,),
             (6, 7, 8),
             (9, 11, 12, 14),
             (10,),
             (13, 16),
             (15,),
             (18,)]

        .. NOTE::
            O(1)
        """
        self.map.contractEdge(self.raw)

    @CheckValid
    def isOnSameFace(self, otherTopologicalDemiEdge):
        """
        Not implemented for PrimitiveMutableTopologicalDemiEdge

        EXAMPLES::

            sage: try:
            ....:     PrimitiveMutableTopologicalDemiEdge(None,2).isOnSameFace(1)
            ....: except:
            ....:     print("OK")
            ....:
            OK

        """
        raise NotImplemented(self)

    @CheckValid
    def isOnSameNode(self, otherTopologicalDemiEdge):
        """
        Not implemented for PrimitiveMutableTopologicalDemiEdge


        EXAMPLES::
            sage: try:
            ....:     PrimitiveMutableTopologicalDemiEdge(None,2).isOnSameNode(1)
            ....: except:
            ....:     print("OK")
            ....:
            OK

        """
        raise NotImplemented(self)
