from sage.graphs.maps.labelled_map import LabelledMap


class RootedMap(LabelledMap):
    """
    This class represents a rooted map,it inherits from LabelledMap.
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
        labelledMap=None,
        sigma=None,
        alpha=None,
        adj=None,
        isAlreadyCanonical=False,
        trust=False,
    ):
        r"""
        Initializes the rooted map from either the permutations alpha
        and sigma, or an adjacency list or a Label

        INPUT:
        - ``sigma`` -- Permutation | MapPermutation | None; Permutation ; Permutation that maps a half-edge
          to the half-edge incident to it in anti-clockwise direction around
          the vertex it belongs to.
        - ``alpha`` -- Permutation | MapPermutation | None; Permutation that maps a half-edge
            Fixed-point free involution whose cycles are given by the edges.
        - ``ajd``-- List[Tuple] | None ; and adjacency list be careful the order of the
            node in your adjaceny will be used to choose the embedding

        - ``isAlreadyCanonical`` -- bool ; A parameter that indicates sigma,alpha given are already in canonical form
          i.e represent a the canonical representant of the rooted map

        - ``trust`` -- bool ; A parameter that indicates whether the validity check (i.e., whether the map is connex, etc.)
          should be skipped when initializing the map. It makes initialization faster but can be dangerous because
          if the map isn't well-formed, all the other methods become unsafe. You should be absolutely sure of your
          map's validity if you set this to true.

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)


        """

        self._production = True
        if labelledMap is None:
            labelledMap = LabelledMap(
                sigma=sigma, alpha=alpha, adj=adj, trust=trust)
        canonicalRepresentant = labelledMap
        if not isAlreadyCanonical:
            canonicalRepresentant = labelledMap.canonicalRepresentant()
        super().__init__(
            canonicalRepresentant.sigma, canonicalRepresentant.alpha, trust=self._production
        )

    def tetravalance(self):
        """
        This method provides a bijection between rooted maps
        with m edges of genus g and face-bicolorable tetravalent
        rooted maps of genus g with m vertices. This function returns
        the rooted face-bicolorable tetravalance associated with self.

        OUTPUT:

        A tetravalent face-bicolorable rooted map
        associated with self.

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19]) 
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])  
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.tetravalance().nodes()
            [(1, 2, 4, 7),
             (3, 6, 9, 11),
             (5, 8, 10, 12),
             (13, 15, 17, 20),
             (14, 16, 19, 23),
             (18, 22, 25, 29),
             (21, 24, 27, 32),
             (26, 31, 35, 30),
             (28, 34, 38, 40),
             (33, 37, 39, 36)]

        NOTE:

            Complexity is O(m), where m is the number of edges.
        """

        return RootedMap(
            labelledMap=super().tetravalance(), isAlreadyCanonical=True, trust=self._production
        )

    def edgeMap(self):
        """
        A method that return the edge map of this map

        OUTPUT:

        The edge map of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.edgeMap().nodes()
            [(1, 2, 4, 7),
             (3, 6, 9, 11),
             (5, 8, 10, 12),
             (13, 15, 17, 20),
             (14, 16, 19, 23),
             (18, 22, 25, 29),
             (21, 24, 27, 32),
             (26, 31, 35, 30),
             (28, 34, 38, 40),
             (33, 37, 39, 36)]

        NOTE:

            Complexity is O(m), where m is the number of edges.


        """

        return RootedMap(
            labelledMap=super().edgeMap(), isAlreadyCanonical=True, trust=self._production,
        )

    def incidenceMap(self):
        """
        A method that return the incidence map of this map

        OUTPUT:

        The incidence map of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.incidenceMap().faces()
            [(1, 6, 16, 20),
             (2, 9, 15, 3),
             (4, 8, 7, 5),
             (10, 19, 24, 11),
             (12, 23, 14, 13),
             (17, 27, 30, 18),
             (21, 26, 25, 22),
             (28, 34, 40, 29),
             (31, 33, 35, 32),
             (36, 39, 38, 37)]


        NOTE:

            Complexity is O(m), where m is the number of edges.

        """

        return RootedMap(
            labelledMap=super().incidenceMap(), isAlreadyCanonical=True, trust=self._production
        )

    def quadrangulation(self):
        """
        This function provides a bijection between rooted maps
        of genus g with m edges and bipartite rooted quadrangulations
        of genus g with m faces.


        OUTPUT:

        The bipartite quadrangulation associated to self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.quadrangulation().faces()
            [(1, 6, 16, 20),
             (2, 9, 15, 3),
             (4, 8, 7, 5),
             (10, 19, 24, 11),
             (12, 23, 14, 13),
             (17, 27, 30, 18),
             (21, 26, 25, 22),
             (28, 34, 40, 29),
             (31, 33, 35, 32),
             (36, 39, 38, 37)] 

        NOTE:

            Complexity is O(m), where m is the number of edges.

        """
        return RootedMap(
            labelledMap=super().quadrangulation(), isAlreadyCanonical=True, trust=self._production
        )

    def derivedMap(self):
        """
        OUTPUT:

        The derived map of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: dm = m.derivedMap()
            sage: dm.m,dm.g,dm.f
            (40, 0, 20)

        NOTE:

            Complexity is O(m), where m is the number of edges.

        """
        return RootedMap(
            labelledMap=super().derivedMap(), isAlreadyCanonical=True, trust=self._production
        )

    def dual(self):
        """
        OUTPUT:

        The dual  of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.dual().pretty_print()
                        Alpha: [(1, 3), (2, 5), (4, 6), (7, 9), (8, 11), (10, 13), (12, 14), (15, 16), (17, 19), (18, 20)]
                        Sigma (Node): [(1, 2, 4, 6, 8, 10, 12, 14, 13, 15, 16, 17, 18, 20, 19, 11, 5, 7, 9, 3)]
                        Phi (Face): [(1,), (2, 7, 3), (4, 8, 5), (6,), (9,), (10, 15, 17, 11), (12, 13), (14,), (16,), (18, 19), (20,)]

        NOTE:

            Complexity is O(m), where m is the number of edges.

        """

        return RootedMap(labelledMap=super().dual(), trust=self._production)

    def __repr__(self):
        r"""
        Return string representation of self

        OUTPUT: The string representation of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m
            Rooted map | Sigma : [2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20] Alpha : [3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19]

        """

        return (
            "Rooted map | Sigma : "
            + str(self.sigma)
            + " Alpha : "
            + str(self.alpha)
        )

    def inverseQuadrangulation(self):
        """
        This function is the inverse of quadrangulation given that self
        is a bipartite rooted quadrangulation. It returns the only
        rooted map M such that M.quadrangulation() = self. If self
        isn't a bipartite quadrangulation, it will raise an error.

        OUTPUT:
            The inverse of self from quadrangulation if self is a rooted
            bipartite quadrangulation; otherwise, raises an error.

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.quadrangulation().inverseQuadrangulation() == m
            True

        NOTE:

            Complexity is O(m), where m is the number of edges.

        """
        return RootedMap(
            labelledMap=super().inverseQuadrangulation(),
            isAlreadyCanonical=True,
            trust=self._production,
        )

    def relabel(self, tau):
        """
        This method, inherited from LabelledMap, is not applicable to
        RootedMap. It will simply return a copy of self.

        OUTPUT:
        A copy of self

        EXAMPLES::
            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.relabel(tau = Permutation([(1,2)])) == m
            True

        NOTE:

            Complexity is O(m), where m is the number of edges.

        """
        return RootedMap(labelledMap=self, isAlreadyCanonical=True, trust=self._production)

    def schaefferTree(self, markedDemiEdge):
        """
        The Schaeffer surjection from rooted bipartite quadrangulation of
        genus g with k faces and a marked node to a rooted one-face map
        (tree if g=0) of genus g with k edges and a labelling of its nodes.

        The function returns a rooted one-face map associated with self
        and a labelling on its demi-edges such that f(node) is the common
        value of all its demi-edges. If self isn't a bipartite
        quadrangulation, this function will raise an error.

        INPUT:
        - ``markedDemiEdge`` -- int ; A demi-edge on the node which is marked.

        OUTPUT:
            A couple (tree,labelling),
            tree: A rooted tree corresponding to the above description.
            labelling:A list of labellings on the demi-edges of tree.

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: quad = m.quadrangulation()
            sage: quad.schaefferTree(1)
            (Rooted map | Sigma : [1, 3, 4, 2, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20] Alpha : [2, 1, 5, 6, 3, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19],
             [-1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

        NOTE:

            Complexity is O(m), where m is the number of edges.

        """
        tree, labelled = super().schaefferTree(markedDemiEdge=markedDemiEdge)
        return RootedMap(labelledMap=tree, isAlreadyCanonical=True, trust=self._production), labelled

    def inverseShaefferTree(self, labelled, returnMarkedDemiEdge=True):
        """
        This method is the inverse of schaefferTree. Given that self is a
        one-face map, it returns a quadruple (quadA, quadB, markedDemiEdgeA,
        markedDemiEdgeB) where quadA and quadB are rooted quadrangulations.

        If returnMarkedDemiEdge=False, it will return only (quadA, quadB).
        If self isn't a one-face map, it will raise an error.

        INPUT:
        - ``labelled`` -- List[int] ; A list of size 2*m+1 such that for the demi-edge i,
        labelled[i] is the label of its attached node.
        - ``returnMarkedDemiEdge`` -- bool ;  Whether or not to return marked demi-edges
        (default is True).

        OUTPUT:
            (quadA, quadB, markedDemiEdgeA, markedDemiEdgeB) if
            returnMarkedDemiEdge=True; otherwise, (quadA, quadB).

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: quad = m.quadrangulation()
            sage: tree,labelled = quad.schaefferTree(1)
            sage: quadA,quadB,markedDemiEdgeA,markedDemiEdgeB = tree.inverseShaefferTree(labelled)
            sage: quadA == quad
            False
            sage: quadB == quad
            True

        NOTE:

            Complexity is O(m), where m is the number of edges.


        """
        if returnMarkedDemiEdge:
            (
                quadA,
                quadB,
                markedDemiEdgeA,
                markedDemiEdgeB,
            ) = super().inverseShaefferTree(
                labelled, returnMarkedDemiEdge=returnMarkedDemiEdge
            )

            return (
                RootedMap(labelledMap=quadA,
                          isAlreadyCanonical=True, trust=self._production),
                RootedMap(labelledMap=quadB,
                          isAlreadyCanonical=True, trust=self._production),
                markedDemiEdgeA,
                markedDemiEdgeB,
            )
        quadA, quadB = super().inverseShaefferTree(
            labelled, returnMarkedDemiEdge=returnMarkedDemiEdge
        )
        return RootedMap(
            labelledMap=quadA, isAlreadyCanonical=True, trust=self._production
        ), RootedMap(labelledMap=quadB, isAlreadyCanonical=True, trust=self._production)

    def copy(self):
        """
        OUTPUT:

            A copy of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.copy() == m
            True

        NOTE:

            Complexity is O(m), where m is the number of edges.

        """
        return RootedMap(sigma=self.sigma, alpha=self.alpha, trust=self._production, isAlreadyCanonical=True)

    @property
    def root(self):
        """
        OUTPUT:
            The TopologicalDemiEdge associated to the root

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = RootedMap(alpha = alpha,sigma=sigma)
            sage: m.root
            X(1) 

        NOTE:

            Complexity is O(1)

        """
        return self.X(1)
