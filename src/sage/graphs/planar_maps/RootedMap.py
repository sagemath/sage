from sage.graphs.planar_maps.LabelledMap import LabelledMap


class RootedMap(LabelledMap):
    """
    This class represents a rooted map.
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
        -------
        Returns:
            A tetravalent face-bicolorable rooted map
            associated with self.
        -------
        O(m), where m is the number of edges.
        """
        return RootedMap(
            labelledMap=super().tetravalance(), isAlreadyCanonical=True, trust=self._production
        )

    def edgeMap(self):
        """
        Returns the edge map of the rooted map.
        -------
        Returns:
            The edge map of self.
        -------
        O(m), where m is the number of edges.
        """
        return RootedMap(
            labelledMap=super().edgeMap(), isAlreadyCanonical=True, trust=self._production,
        )

    def incidenceMap(self):
        """
        Returns the incidence map of the rooted map.
        -------
        Returns:
            The incidence map of self.
        -------
        O(m), where m is the number of edges.
        """
        return RootedMap(
            labelledMap=super().incidenceMap(), isAlreadyCanonical=True, trust=self._production
        )

    def quadrangulation(self):
        """
        This function provides a bijection between rooted maps
        of genus g with m edges and bipartite rooted quadrangulations
        of genus g with m faces.
        -------
        Returns:
            A tetravalent bi-colorable rooted map associated with self.
        -------
        O(m), where m is the number of edges.
        """
        return RootedMap(
            labelledMap=super().quadrangulation(), isAlreadyCanonical=True, trust=self._production
        )

    def derivedMap(self):
        """
        Returns the derived map of self.
        -------
        Returns:
            The derived map of self.
        -------
        O(m)
        """
        return RootedMap(
            labelledMap=super().derivedMap(), isAlreadyCanonical=True, trust=self._production
        )

    def dual(self):
        """
        Returns the dual of the rooted map.
        -------
        Returns:
            The dual of self.
        -------
        O(m), where m is the number of edges.
        """
        return RootedMap(labelledMap=super().dual(), trust=self._production)

    def __repr__(self):
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
        -------
        Returns:
            The inverse of self from quadrangulation if self is a rooted
            bipartite quadrangulation; otherwise, raises an error.
        -------
        O(m), where m is the number of edges.
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
        -------
        Args:
            markedDemiEdge : A demi-edge on the node which is marked.
        Returns:
            - tree: A rooted tree corresponding to the above description.
            - labelling: A list of labellings on the demi-edges of tree.
        -------
        O(m), where m is the number of edges.
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
        -------
        Args:
            labelled : A list of size 2*m+1 such that for the demi-edge i,
                       labelled[i] is the label of its attached node.
            returnMarkedDemiEdge : Whether or not to return marked demi-edges
                                   (default is True).
        Returns:
            (quadA, quadB, markedDemiEdgeA, markedDemiEdgeB) if
            returnMarkedDemiEdge=True; otherwise, (quadA, quadB).
        -------
        O(m), where m is the number of edges.
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
        Returns:
            A copy of self
        ----
        O(m)
        """
        return RootedMap(sigma=self.sigma, alpha=self.alpha, trust=self._production)

    @property
    def root(self):
        """
        Returns:
            The TopologicalDemiEdge associated to the root
        ----
        O(1)
        """
        return self.X(1)
