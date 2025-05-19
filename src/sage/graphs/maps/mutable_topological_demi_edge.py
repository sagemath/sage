"""Define the MutableTopologicalDemiEdge class, an abstraction an abstraction meant to represent a demi-edge of a map in a more user-friendly way than raw indices used in MutableLabelledMap."""

from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sage.graphs.maps.mutable_labelled_map import MutableLabelledMap        # see topological_demi_edge.py for explanations


class MutableTopologicalDemiEdge(TopologicalDemiEdge):
    """
    This class implement a splecial version of TopologicalDemiEdge used in MutableLabelledMap. Their specificity is that they have methods to modify the map: for instance, you can delete them, add edges, etc.
    """

    def __init__(self, lmap: "MutableLabelledMap", index: int):
        """This class is an abstraction meant to represent the demi edge of a mutable map in a more user-friendly way
        than simple indexes. It is more related to the "topological structure" of the map than the raw index.

        INPUT:

        - ``lmap`` -- MutableLabelledMap: The map to which the demi-edge is bind
        - ``index`` -- int: The index associated to the demi-edge

        EXAMPLES::

            sage: from sage.graphs.maps.mutable_topological_demi_edge import MutableTopologicalDemiEdge
            sage: lm = MutableLabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: MutableTopologicalDemiEdge(lm, 3)
            X(3)

        NOTE:

            O(1)
        """

        super().__init__(lmap, index)

    @property
    def map(self) -> "MutableLabelledMap":
        """
        The map to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.mutable_topological_demi_edge import MutableTopologicalDemiEdge

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: MutableTopologicalDemiEdge(lm, 3).map
            Labelled map | Sigma : [3, 2, 1, 4], Alpha : [2, 1, 4, 3]

        NOTE:

            O(1)
        """

        return self.getMap()

    def getMap(self) -> "MutableLabelledMap":
        """
        Returns the map to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.mutable_topological_demi_edge import MutableTopologicalDemiEdge

            sage: lm = LabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: MutableTopologicalDemiEdge(lm, 3).getMap()
            Labelled map | Sigma : [3, 2, 1, 4], Alpha : [2, 1, 4, 3]

        NOTE:

            O(1)
        """
        self._checkValid()
        return self._lmap

    def delete(self, trust=False) -> None:
        """
        This function delete self from self.map in case of planar map it is efficient O(log(m)) otherwise O(m) for higher genus
        map if trust = False. If trust = True as a parameter it is always of complexity O(log(m)), but it can break the invariant
        saying that the self.map is connected thus make all the other method unsafe by default trust = False.

        INPUT:

        -``trust`` -- bool ; a parameter telling the function if it should trust the fact that the map will stay connected after
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

        NOTE:

            O(log(m)) if self is planar or trust = True otherwise O(m)
        """
        self._checkValid()

        self.map.deleteEdge(self.raw, trust=trust)

    def link(self, otherTopoDemiEdge: "MutableTopologicalDemiEdge") -> "tuple[MutableTopologicalDemiEdge, MutableTopologicalDemiEdge]":
        """
        This will add an edge between the node of self to otherTopoDemiEdge(note that they
        need to be on the same node otherwise this will raise an error), the edge will be added as follow ,
        let denote (A,B) the demi edges composing this new edge A will be on the same node as self but before it
        and B on the same node as otherTopoDemiEdge but before it.
        It will return two MutableLabelledMap topoDemiEdgeA,topoDemiEdgeB corresponding to the new demi edge A and B

        INPUT:

        - ``otherTopoDemiEdge`` -- MutableTopologicalDemiEdge ; Another MutableTopologicalDemiEdge on the same facee as self

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

        NOTE:

            O(log(m))
        """
        self._checkValid()

        return self.map.addEdge(self.raw, otherTopoDemiEdge.raw)

    def addEdgeAfter(self) -> "MutableTopologicalDemiEdge":
        """
        This method will create a new edge, such that it is added on the same node as self but after it in the
        trigonometric order

        OUTPUT:

            The MutableTopologicalDemiEdge associate to demi edge of the new edge which is on the newly added node

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: A.addEdgeBefore()
            X(22)
            sage: A.addEdgeAfter()
            X(24)

        NOTE:

            O(log(m))
        """
        self._checkValid()
        return self.map.addEdgeAfter(self.raw)

    def addEdgeBefore(self) -> "MutableTopologicalDemiEdge":
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

        NOTE:

            O(log(m))
        """
        self._checkValid()

        return self.map.addEdgeBefore(self.raw)

    def deleteNode(self, trust=False) -> None:
        """
        This method will delete the node attached to self, when trust = False (default) it is efficient in case of planar map O(deg(node)*log(m))
        but for higher genus map it is of complexity O(m+deg(node)*log(m)) where m is the number of edge in self.map.
        If trust = True as a parameter it is always of complexity O(log(m)), but it can break the invariant saying that
        the self.map is connected thus make all the other method unsafe by default trust = False.
        It will raise an error if the self.map isn't connected after the operation and trust is set to False.

        INPUT:

        - ``trust`` -- bool ; a parameter telling the method if it should trust the fact that the self.map will stay connected after deleting self default is False

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

        NOTE:

            O(deg(node)*log(m)) if self.map is planar and O(m+deg(node)*log(m)) otherwise
        """
        self._checkValid()

        self.map.deleteNode(self.raw, trust=trust)

    def contract(self) -> None:
        """
        Contract the edge bind to self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm  = MutableLabelledMap(sigma = sigma,alpha=alpha)
            sage: A = mm.X(5)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: A.contract()
            sage: mm.faces()
            [(1, 3, 4, 7, 11, 16, 18, 13, 12, 15, 14, 2, 5, 17, 9, 8, 10, 6)]

        NOTE:

            O(log(m))
        """
        self._checkValid()
        self.map.contractEdge(self.raw)

    def contractFace(self) -> None:
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

        NOTE:

            O(tlog(m)) where t is the number of edge on the face containing self
        """
        self._checkValid()
        self.map.contractFace(self.raw)

    def mergeMap(self, otherTopoDemiEdge: "MutableTopologicalDemiEdge") -> "tuple[MutableTopologicalDemiEdge, list[MutableTopologicalDemiEdge]]":
        """
        This will merge in self.map without modifying otherTopoDemiEdge.map(if otherTopoDemiEdge.map isn't the same object as self.map),
        what we mean by merging is the following draw an edge between the two map.
        It will return a couple (topoDemiEdge,topoDemiEdgeList) where topoDemiEdge is a MutableTopologicalDemiEdge associated
        to the edge that was added on the side of self.map, and topoDemiEdgeList which will contain a list of MutableTopologicalDemiEdge
        such that if i was a demi edge in otherDemiEdge.map , topoDemiEdgeList[i] is the MutableTopologicalDemiEdge associated  in self.map
        to the copy of i.

        INPUT:

        -``otherTopoDemiEdge`` -- MutableTopologicalDemiEdge: The otherTopoDemiEdge attached to the node on otherTopoDemiEdge.map which to draw the new edge

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

        NOTE:

            O(p(log(m)+log(p))) where p = otherTopoDemiEdge.map.m and m is the number of edge of self.map,
            note that it is much more efficient than O(p+m) mainly when m>>p
        """
        self._checkValid()

        return self.map.merge(
            self.raw, otherTopoDemiEdge.map, otherTopoDemiEdge.raw)

    def copyOn(self, otherTopoDemiEdge: "MutableTopologicalDemiEdge") -> "tuple[MutableTopologicalDemiEdge, list[MutableTopologicalDemiEdge]]":
        """
        Given that self is such that it is attached to a node of degree one (otherwise this function will raise
        an error),this function will attach self before otherTopoDemiEdge and then copy otherTopoDemiEdge.map from this point on,
        I will return a MutableTopologicalDemiEdge corresponding to demiEdge in self.map but also a list of MutableTopologicalDemiEdge call it topoDemiEdgeList
        corresponding to the new MutableTopologicalDemiEdge of all the demi edge added from otherTopoDemiEdge.map, such if a demi edge is of index i
        in otherTopoDemiEdge.map than topoDemiEdgeList[i] is the MutableTopologicalDemiEdge corresponding to the copy of the demi edge in self.map.Note that this function
        won't modify otherTopoDemiEdge.map(if otherTopoDemiEdge.map isn't strictly the same object in memory as self.map).

        INPUT:

        - ``otherTopoDemiEdge`` -- MutableTopologicalDemiEdge: Another TopologicalDemiEdge

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

        NOTE:

            O(p(log(m)+log(p))) where p = otherTopoDemiEdge.map.m and m is the number of edge of self.map,
            note that it is much more efficient than O(p+m) mainly when m>>p
        """
        self._checkValid()
        if self.map.numberInTheSameNode(self.raw) > 1:
            raise ValueError(
                "Self isn't attached to a node of degree one cannot copyOn")
        return self.map.copyOnDemiEdge(
            self.raw, otherTopoDemiEdge.map, otherTopoDemiEdge.raw)

    def mergeNode(self, otherTopoDemiEdge: "MutableTopologicalDemiEdge") -> None:
        """
        Merge the node attached to self and otherTopoDemiEdge, they need to be on the same face

        INPUT:

        - ``otherTopoDemiEdge`` -- MutableTopologicalDemiEdge;the other TopologicalDemiEdge on the same face as self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
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

        NOTE:

            O(log(m)) where m is the number of edge of self.map
        """
        self._checkValid()

        if self.map is not otherTopoDemiEdge.map:
            raise ValueError(
                "Cannot mergeNode between two demi edge on different map")

        self.map.mergeNode(self.raw, otherTopoDemiEdge.raw)
        return
