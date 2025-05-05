from sage.graphs.maps.map_decorator import CheckValid


@CheckValid
class TopologicalDemiEdge():
    def __init__(self, lmap, index):
        """
        This class is a an abstraction mean to represent the demi edge of a map in a more user 
        friendly way than simple index, it is more link to the "topology" of the map than the raw index

        INPUT:
            lmap: The map on which the self is bind
            index: The index associated to the demi edge it is bind

        EXAMPLES::
            sage: TopologicalDemiEdge(None,None)
            X(None)

        .. NOTE::
        O(1)
        """
        self._lmap = lmap
        self._index = index
        self._isValid = True

    @property
    def raw(self):
        """
        OUTPUT:

        The raw index associated to demi edge on which self is bind to

        EXAMPLES::

            sage: TopologicalDemiEdge(None,5).raw
            5

        .. NOTE::
        O(1)
        """
        return self.getIndex()

    @property
    def map(self):
        """
        The map on which self is bind to

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A.map == lm
            True

        .. NOTE::
        O(1)
        """
        return self.getMap()

    @property
    def c(self):
        """
        OUTPUT:
        The next TopologicalDemiEdge  on the edge on which self is on

        EXAMPLES:
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: A = lm.X(1)
            sage: A.c
            X(3)

        .. NOTE::
        O(1)
        """
        return self.nextOnEdge()

    @property
    def f(self):
        """
        OUTPUT:

        The next TopologicalDemiEdge on the face on which self is on


        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: A = lm.X(1)
            sage: A.c
            X(3)
            sage: A.f
            X(3)
            sage: A.pf
            X(6)
            sage: A.n
            X(2)
            sage: A.pn
            X(4)

        .. NOTE::
        O(1)
        """
        return self.nextOnFace()

    @property
    def n(self):
        """
        OUTPUT:

        The next TopologicalDemiEdge on the node on which self is on


        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: A = lm.X(1)
            sage: A.c
            X(3)
            sage: A.f
            X(3)
            sage: A.pf
            X(6)
            sage: A.n
            X(2)
            sage: A.pn
            X(4)

        .. NOTE::
        O(1)
        """
        return self.nextOnNode()

    @property
    def pf(self):
        """
        OUTPUT:
        The previous TopologicalDemiEdge on the face on which self is on

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: A = lm.X(1)
            sage: A.c
            X(3)
            sage: A.f
            X(3)
            sage: A.pf
            X(6)
            sage: A.n
            X(2)
            sage: A.pn
            X(4)

        .. NOTE::
        O(1)
        """
        return self.prevOnFace()

    @property
    def pn(self):
        """
        OUTPUT:
        The previous TopologicalDemiEdge on the node on which self is on

        EXAMPLES:

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: A = lm.X(1)
            sage: A.c
            X(3)
            sage: A.f
            X(3)
            sage: A.pf
            X(6)
            sage: A.n
            X(2)
            sage: A.pn
            X(4)

        .. NOTE::
        O(1)
        """
        return self.prevOnNode()

    @CheckValid
    def getIndex(self):
        """
        OUTPUT:
        The raw index associated to demi edge on which self is bind to

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.getIndex()
            1 

        .. NOTE::
        O(1)
        """
        return self._index

    @CheckValid
    def getMap(self):
        """
        OUTPUT: 
            The map on which self is bind

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.getIndex()
            1
            sage: A.getMap() == lm
            True

        .. NOTE::
            O(1)
        """
        return self._lmap

    @CheckValid
    def nextOnEdge(self):
        """
        OUTPUT: 
        The next TopologicalDemiEdge  on the edge on which self is on

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.getIndex()
            1
            sage: A.getMap() == lm
            True
            sage: A.nextOnEdge()
            X(3)

        .. NOTE::
        O(1)
        """

        return self.map.getTopologicalDemiEdge(self.map.alpha(self.raw))

    @CheckValid
    def nextOnFace(self):
        """
        OUTPUT: 
            The next TopologicalDemiEdge on the face on which self is on

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.getIndex()
            1
            sage: A.getMap() == lm
            True
            sage: A.nextOnEdge()
            X(3)
            sage: A.nextOnFace()
            X(3) 

        .. NOTE::
        O(1)
        """
        return self.map.getTopologicalDemiEdge(self.map.phi(self.raw))

    @CheckValid
    def nextOnNode(self):
        """
        OUTPUT: 

            The next TopologicalDemiEdge on the node on which self is on

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.getIndex()
            1
            sage: A.getMap() == lm
            True
            sage: A.nextOnEdge()
            X(3)
            sage: A.nextOnFace()
            X(3)
            sage: A.nextOnNode()
            X(2)


        .. NOTE::
        O(1)
        """
        return self.map.getTopologicalDemiEdge(self.map.sigma(self.raw))

    @CheckValid
    def prevOnFace(self):
        """
        OUTPUT: 
            The previous TopologicalDemiEdge on the face on which self is on

        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.getIndex()
            1
            sage: A.getMap() == lm
            True
            sage: A.nextOnEdge()
            X(3)
            sage: A.nextOnFace()
            X(3)
            sage: A.nextOnNode()
            X(2)
            sage: A.prevOnFace()
            X(6)

        .. NOTE::
        O(1)
        """
        return self.map.getTopologicalDemiEdge(
            self.map.phi.inverseApply(self.raw))

    @CheckValid
    def prevOnNode(self):
        """
        OUTPUT: 
            The previous TopologicalDemiEdge on the node on which self is on

        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.getIndex()
            1
            sage: A.getMap() == lm
            True
            sage: A.nextOnEdge()
            X(3)
            sage: A.nextOnFace()
            X(3)
            sage: A.nextOnNode()
            X(2)
            sage: A.prevOnFace()
            X(6)
            sage: A.prevOnNode()
            X(4)

        .. NOTE::
        O(1)
        """
        return self.map.getTopologicalDemiEdge(
            self.map.sigma.inverseApply(self.raw))

    @CheckValid
    def face(self):
        """
        OUTPUT: 
        A list containing all the TopologicalDemiEdge of demi edge on the
        same face as self

        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.face()
            [X(1),
             X(3),
             X(2),
             X(5),
             X(4),
             X(7),
             X(11),
             X(16),
             X(18),
             X(13),
             X(12),
             X(15),
             X(14),
             X(19),
             X(20),
             X(17),
             X(9),
             X(8),
             X(10),
             X(6)]

        .. NOTE::
        O(f) where f is the number of demi edge on the face
        """
        lst = []
        for e in self.map.demiEdgesOnTheSameFace(self.raw):
            lst.append(self.map.getTopologicalDemiEdge(e))
        return lst

    @CheckValid
    def node(self):
        """
        OUTPUT: A list containing all the TopologicalDemiEdge of demi edge on the
        same node as self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: lm = LabelledMap(sigma=sigma,alpha=alpha)
            sage: A = lm.X(1)
            sage: A.node()
            [X(1), X(2), X(4)]

        .. NOTE::
        O(n) where n is the number of demi edge on the node
        """

        lst = []
        for e in self.map.demiEdgesOnTheSameNode(self.raw):
            lst.append(self.map.getTopologicalDemiEdge(e))

        return lst

    @CheckValid
    def isOnSameFace(self, otherTopologicalDemiEdge):
        """
        INPUT:
            otherTopologicalDemiEdge a TopologicalDemiEdge on the same map as self

        OUTPUT:
            A boolean indicating if they are on the same face

        EXAMPLES::

            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: A = lm.X(1)
            sage: B = lm.X(9)
            sage: A.isOnSameFace(B)
            True
            sage: B.isOnSameNode(A)
            False

        .. NOTE::
        O(1)
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

            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: A = lm.X(1)
            sage: B = lm.X(9)
            sage: A.isOnSameFace(B)
            True
            sage: B.isOnSameNode(A)
            False

        .. NOTE::
        O(1)
        """

        return self.map.areOnTheSameNode(
            self.raw, otherTopologicalDemiEdge.raw)

    def _invalidate(self):
        """
        Make self invalid

        EXAMPLES::
            sage: TopologicalDemiEdge(None,5)._invalidate()

        .. NOTE::
        O(1)
        """
        self._isValid = False

    def __repr__(self):
        """
        OUTPUT:
            A string representation of self

        EXAMPLES::

            sage: TopologicalDemiEdge(None,5)
            X(5) 

        """
        return f"X({self.raw})"

    def _setIndex(self, newIndex):
        """
        Change the index of self

        INPUT: 
            newIndex

        EXAMPLES::

            sage: U  = TopologicalDemiEdge(None,5)
            sage: U
            X(5)
            sage: U._setIndex(22)
            sage: U
            X(22) 

        .. NOTE::
        O(1)
        """
        self._index = newIndex

    def _swapIndex(self, otherTopoDemiEdge):
        """
        Swap indexes with otherTopoDemiEdge

        INPUT: 
            otherTopoDemiEdge

        EXAMPLES::

            sage: U  = TopologicalDemiEdge(None,5)
            sage: V  = TopologicalDemiEdge(None,10)
            sage: U,V
            (X(5), X(10))
            sage: U._swapIndex(V)
            sage: U,V
            (X(10), X(5))

        .. NOTE::
        O(1)
        """
        tmp = self.raw
        self._setIndex(otherTopoDemiEdge.raw)
        otherTopoDemiEdge._setIndex(tmp)
