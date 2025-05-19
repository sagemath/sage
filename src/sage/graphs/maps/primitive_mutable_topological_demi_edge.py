"""Define the PrimitiveMutableTopologicalDemiEdge class, an abstraction meant to represent a demi-edge of a map in a more user-friendly way than raw indices used in MutableLabelledMap."""

from sage.graphs.maps.map_error import NotImplementedError
from sage.graphs.maps.topological_demi_edge import TopologicalDemiEdge

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sage.graphs.maps.primitive_mutable_labelled_map import PrimitiveMutableLabelledMap        # see topological_demi_edge.py for explanations


class PrimitiveMutableTopologicalDemiEdge(TopologicalDemiEdge):
    """
    This class implement a special version of TopologicalDemiEdge used in PrimitiveMutableLabelledMap
    their specificity is that they have method to modify the map for instance you can delete them
    ,add between them etc.
    """

    def __init__(self, lmap: "PrimitiveMutableLabelledMap", index: int):
        """This class is an abstraction meant to represent the demi edge of a mutable map in a more user-friendly way
        than simple indexes. It is more related to the "topological structure" of the map than the raw index.

        INPUT:

        - ``lmap`` -- PrimitiveMutableLabelledMap: The map to which the demi-edge is bind
        - ``index`` -- int: The index associated to the demi-edge

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_mutable_topological_demi_edge import PrimitiveMutableTopologicalDemiEdge
            sage: lm = PrimitiveMutableLabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: PrimitiveMutableTopologicalDemiEdge(lm, 3)
            X(3)

        NOTE:

            O(1)"""

        super().__init__(lmap, index)

    @property
    def map(self) -> "PrimitiveMutableLabelledMap":
        """
        The map to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_mutable_topological_demi_edge import PrimitiveMutableTopologicalDemiEdge

            sage: lm = PrimitiveMutableLabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: PrimitiveMutableTopologicalDemiEdge(lm, 3).map
            Labelled map | Sigma : [3, 2, 1, 4], Alpha : [2, 1, 4, 3]

        NOTE:

            O(1)
        """

        return self.getMap()

    def getMap(self) -> "PrimitiveMutableLabelledMap":
        """
        Returns the map to which self is bind.

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_mutable_topological_demi_edge import PrimitiveMutableTopologicalDemiEdge

            sage: lm = PrimitiveMutableLabelledMap(Permutation([(1,3), (2,), (4,)]), Permutation([(1,2), (3,4)]))
            sage: PrimitiveMutableTopologicalDemiEdge(lm, 3).getMap()
            Labelled map | Sigma : [3, 2, 1, 4], Alpha : [2, 1, 4, 3]

        NOTE:

            O(1)
        """
        self._checkValid()
        return self._lmap

    def delete(self, sameFace: bool) -> None:
        """
        This function delete self from self.map

        INPUT:

        - ``sameFace`` -- bool ; a boolean indicating if self and self.c are on the same face or not

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

        NOTE:

            O(1),If this break the connectivity, no error will be raised and the map won't be stable anymore
        """
        self._checkValid()

        self.map.deleteEdge(self.raw, sameFace)

    def link(self, otherTopoDemiEdge: "PrimitiveMutableTopologicalDemiEdge") -> tuple["PrimitiveMutableTopologicalDemiEdge", "PrimitiveMutableTopologicalDemiEdge"]:
        """
        This will add an edge between the node of self to otherTopoDemiEdge(note that they
        need to be on the same node otherwise nothing is guaranteed), the edge will be added as follow ,
        let denote (A,B) the demi edges composing this new edge A will be on the same node as self but before it
        and B on the same node as otherTopoDemiEdge but before it.
        It will return two MutableLabelledMap topoDemiEdgeA,topoDemiEdgeB corresponding to the new demi edge A and B

        INPUT:

        - ``otherTopoDemiEdge`` -- PrimitiveMutableTopologicalDemiEdge ; Another PrimitiveMutableTopologicalDemiEdge on the same face as self

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

        NOTE:

            O(1)
        """
        self._checkValid()

        return self.map.addEdge(self.raw, otherTopoDemiEdge.raw)

    def addEdgeAfter(self) -> "PrimitiveMutableTopologicalDemiEdge":
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

        NOTE:

            O(1)
        """
        self._checkValid()
        return self.map.addEdgeAfter(self.raw)

    def addEdgeBefore(self) -> "PrimitiveMutableTopologicalDemiEdge":
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

        NOTE:

            O(1)
        """
        self._checkValid()

        return self.map.addEdgeBefore(self.raw)

    def contract(self) -> None:
        """
        Contract the edge bind to self , note that self cannot be on a loop
        i.e self and self.c cannot be on the same node otherwise nothing is guaranteed

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

        NOTE:

            O(1)
        """
        self._checkValid()
        self.map.contractEdge(self.raw)

    def isOnSameFace(self, otherTopologicalDemiEdge: TopologicalDemiEdge) -> bool:
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
        self._checkValid()
        raise NotImplementedError(self)

    def isOnSameNode(self, otherTopologicalDemiEdge: TopologicalDemiEdge) -> bool:
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
        self._checkValid()
        raise NotImplementedError(self)
