from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge
from sage.graphs.maps.map_decorator import CheckValid


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

        INPUT:
            trust: a parameter telling the function if it should trust the fact that the map will stay connected after
                deleting demiEdge default is False

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = MutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: cp = mm.copy()
            sage: A = mm.X(5)
            sage: B = mm.X(13)
            sage: U,V = A.link(B)
            sage: mm == cp
            False
            sage: U.delete()
            sage: mm == cp
            True

        .. NOTE::
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

        INPUT:
            otherTopoDemiEdge: Another TopologicalDemiEdge on the same facee as self


        OUTPUT:
            topoDemiEdgeA,topoDemiEdgeB as described above


        EXAMPLES::

            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(1)
            sage: B = mm.X(17)
            sage: A.isOnSameFace(B)
            True
            sage: A.link(B)
            (X(22), X(21))
            sage: A.isOnSameFace(B)
            False

        .. NOTE:
        O(log(m))
        """

        return self.map.addEdge(self.raw, otherTopoDemiEdge.raw)

    @CheckValid
    def addEdgeAfter(self):
        """
        This method will create a new edge, such that it is added on the same node as self but after it in the
        trigonometric order

        INPUT: The TopologicalDemiEdge associate to demi edge of the new edge which is on the newly added node


        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: A.addEdgeBefore()
            X(22)
            sage: A.addEdgeAfter()
            X(24)

        .. NOTE::
        O(log(m))
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
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: A.addEdgeBefore()
            X(22)
            sage: A.addEdgeAfter()
            X(24)   

        .. NOTE:
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

        INPUT:
            trust: a parameter telling the method if it should trust the fact that the self.map will stay connected after
                deleting self default is False

        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
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
            sage: A = mm.X(5)
            sage: A.deleteNode()
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
        O(deg(node)*log(m)) if self.map is planar and O(m+deg(node)*log(m)) otherwise
        """

        self.map.deleteNode(self.raw, trust=trust)

    @CheckValid
    def contract(self):
        """
        Contract the edge bind to self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: A = mm.X(5)
            sage: mm.faces()
            [(1, 3, 4, 7, 11, 16, 18, 13, 12, 15, 14, 2, 5, 17, 9, 8, 10, 6)]
            sage: A.contract()
            sage: mm.faces()
            [(1, 3, 4, 7, 11, 16, 5, 13, 12, 15, 14, 2, 9, 8, 10, 6)]

        .. NOTE::
        O(log(m))
        """
        self.map.contractEdge(self.raw)

    @CheckValid
    def contractFace(self):
        """
        Contract the face on which self is in

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: B = mm.X(13)
            sage: A.link(B)
            (X(22), X(21))
            sage: mm
            Labelled map | Sigma : [2, 4, 3, 1, 22, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 21, 19, 18, 17, 20, 13, 5], Alpha : [3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19, 22, 21]
            sage: mm.faces()
            [(1, 3, 2, 22, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (4, 7, 11, 16, 18, 21, 5)]
            sage: C = mm.X(4)
            sage: A.link(C)
            (X(24), X(23))
            sage: mm.faces()
            [(1, 3, 2, 22, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (4, 7, 11, 16, 18, 21, 24),
             (5, 23)]
            sage: A.contractFace()
            sage: mm.faces()
            [(1, 3, 2, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6), (4, 7, 11, 16, 18, 5)]

        .. NOTE::
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

        INPUT:
            otherTopoDemiEdge: The otherTopoDemiEdge attached to the node on otherTopoDemiEdge.map which to draw the new edge

        OUTPUT:

            newTopoDemiEdge:
                A MutableTopologicalDemiEdge such that it corresponds to a new demi edge attached
                after demiEdge.

            topoDemiEdgeList:
                A list of MutableTopologicalDemiEdge corresponding to the description above

        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: mm.m
            10
            sage: A.mergeMap(A)
            (X(22),
             [X(23),
              X(24),
              X(25),
              X(26),
              X(27),
              X(28),
              X(29),
              X(30),
              X(31),
              X(32),
              X(33),
              X(34),
              X(35),
              X(36),
              X(37),
              X(38),
              X(39),
              X(40),
              X(41),
              X(42)])
            sage: mm.m
            21
        .. NOTE::
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

        INPUT:
            otherTopoDemiEdge: Another TopologicalDemiEdge

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: cp = mm.copy()
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
            sage: A = mm.X(5)
            sage: B = cp.X(1)
            sage: lst = A.copyOn(B)
            sage: mm.m,cp.m
            (20, 10)

        .. NOTE::
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
        INPUT:
            otherTopoDemiEdge: the other TopologicalDemiEdge on the same face as sel

        EXAMPLES::
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(1)
            sage: B = mm.X(17)
            sage: A.isOnSameFace(B)
            True
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
            sage: A.mergeNode(B)
            sage: mm.nodes()
            [(1, 2, 4, 17, 19),
             (3,),
             (5,),
             (6, 7, 8),
             (9, 11, 12, 14),
             (10,),
             (13, 16),
             (15,),
             (18,),
             (20,)]

        .. NOTE::
        O(log(m)) where m is the number of edge of self.map
        """

        if self.map != otherTopoDemiEdge.map:
            raise ValueError(
                "Cannot mergeNode between two demi edge on different map")
        self.map.mergeNode(self.raw, otherTopoDemiEdge.raw)

    @CheckValid
    def isOnSameFace(self, otherTopologicalDemiEdge):
        """
        INPUT:
            otherTopologicalDemiEdge a TopologicalDemiEdge on the same map as self

        OUTPUT:
            A boolean indicating if they are on the same face

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(1)
            sage: B = mm.X(17)
            sage: A.isOnSameFace(B)
            True
            sage: A.isOnSameNode(B)
            False

        .. NOTE::
        O(log(m)) where m is the number of edge of self.map
        """
        return self.map.areOnTheSameFace(
            self.raw, otherTopologicalDemiEdge.raw)

    @CheckValid
    def isOnSameNode(self, otherTopologicalDemiEdge):
        """
        INPUT:
            otherTopologicalDemiEdge a TopologicalDemiEdge on the same map as self

        OUTPUT:
            A boolean indicating if they are on the same node

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(1)
            sage: B = mm.X(17)
            sage: A.isOnSameFace(B)
            True
            sage: A.isOnSameNode(B)
            False

        .. NOTE::
        O(log(m)) where m is the number of edge of self.map
        """

        return self.map.areOnTheSameNode(
            self.raw, otherTopologicalDemiEdge.raw)
