r"""
Cluster complex (or generalized dual associahedron)

EXAMPLES:

A first example of a cluster complex::

    sage: C = ClusterComplex(['A', 2]); C
    Cluster complex of type ['A', 2] with 5 vertices and 5 facets

Its vertices, facets, and minimal non-faces::

    sage: C.vertices()
    (0, 1, 2, 3, 4)

    sage: C.facets()
    [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

    sage: for F in C.facets(): F.cluster()
    [(-1, 0), (0, -1)]
    [(-1, 0), (0, 1)]
    [(0, -1), (1, 0)]
    [(1, 0), (1, 1)]
    [(1, 1), (0, 1)]

    sage: C.minimal_nonfaces()
    [[0, 2], [0, 3], [1, 3], [1, 4], [2, 4]]

We can do everything we can do on simplicial complexes,
e.g. computing its homology::

    sage: C.homology()
    {0: 0, 1: Z}

AUTHORS:

- Christian Stump (2011) Initial version
"""

# ****************************************************************************
#       Copyright (C) 2011      Christian Stump <christian.stump@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.coxeter_groups import CoxeterGroups
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.combinat.subword_complex import SubwordComplex, SubwordComplexFacet
from sage.rings.semirings.non_negative_integer_semiring import NN


class ClusterComplexFacet(SubwordComplexFacet):
    r"""
    A cluster (i.e., a facet) of a cluster complex.
    """
    def cluster(self):
        """
        Return this cluster as a set of almost positive roots.

        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1))
            sage: F.cluster()
            [(-1, 0), (0, -1)]
        """
        if self.parent().k() != 1:
            raise NotImplementedError("not working for multi-cluster complexes")
        F = self.parent().greedy_facet(side="positive")
        R = F.extended_root_configuration()
        N = len(list(F))
        return [-R[i] if i < N else R[i] for i in self]

    def upper_cluster(self):
        """
        Return the part of the cluster that contains positive roots

        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1))
            sage: F.upper_cluster()
            []
        """
        return [beta for beta in self.cluster() if sum(beta) > 0]

    def product_of_upper_cluster(self):
        """
        Return the product of the upper cluster in reversed order.

        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: for F in C: F.product_of_upper_cluster().reduced_word()
            []
            [2]
            [1]
            [1, 2]
            [1, 2]
        """
        W = self.parent().group()
        return W.prod(W.reflections()[beta]
                      for beta in reversed(self.upper_cluster()))


class ClusterComplex(SubwordComplex):
    r"""
    A cluster complex (or generalized dual associahedron).

    The cluster complex (or generalized dual associahedron) is a
    simplicial complex constructed from a cluster algebra.  Its
    vertices are the cluster variables and its facets are the
    clusters, i.e., maximal subsets of compatible cluster variables.

    The cluster complex of type `A_n` is the simplicial complex with
    vertices being (proper) diagonals in a convex `(n+3)`-gon and with
    facets being triangulations.

    The implementation of the cluster complex depends on its
    connection to subword complexes, see [CLS2014]_. Let `c` be a Coxeter
    element with reduced word `{\bf c}` in a finite Coxeter group `W`,
    and let `{\bf w}_\circ` be the `c`-sorting word for the longest
    element `w_\circ \in W`.

    The ``multi-cluster complex`` `\Delta(W,k)` has vertices in
    one-to-one correspondence with letters in the word
    `Q = {\bf c^k w}_\circ` and with facets being complements
    in `Q` of reduced expressions for `w_\circ`.

    For `k = 1`, the multi-cluster complex is isomorphic to the
    cluster complex as defined above.

    EXAMPLES:

    A first example of a cluster complex::

        sage: C = ClusterComplex(['A', 2]); C
        Cluster complex of type ['A', 2] with 5 vertices and 5 facets

    Its vertices, facets, and minimal non-faces::

        sage: C.vertices()
        (0, 1, 2, 3, 4)

        sage: C.facets()
        [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

        sage: C.minimal_nonfaces()
        [[0, 2], [0, 3], [1, 3], [1, 4], [2, 4]]

    We can do everything we can do on simplicial complexes,
    e.g. computing its homology::

        sage: C.homology()
        {0: 0, 1: Z}

    We can also create a multi-cluster complex::

        sage: ClusterComplex(['A', 2], k=2)
        Multi-cluster complex of type ['A', 2] with 7 vertices and 14 facets

    REFERENCES:

    - [CLS2014]_
    """

    @staticmethod
    def __classcall__(cls, W, k=1, coxeter_element=None, algorithm="inductive"):
        r"""
        Standardize input to ensure a unique representation.

        TESTS::

            sage: S1 = ClusterComplex(['B',2])
            sage: W = CoxeterGroup(['B',2])
            sage: S2 = ClusterComplex(W)
            sage: S3 = ClusterComplex(CoxeterMatrix('B2'), coxeter_element=(1,2))
            sage: w = W.from_reduced_word([1,2])
            sage: S4 = ClusterComplex('B2', coxeter_element=w, algorithm="inductive")
            sage: S1 is S2 and S2 is S3 and S3 is S4
            True
        """
        if k not in NN:
            raise ValueError("the additional parameter must be a "
                             "nonnegative integer")

        if W not in CoxeterGroups:
            W = CoxeterGroup(W)

        if not W.is_finite():
            raise ValueError("the Coxeter group must be finite")

        if coxeter_element is None:
            coxeter_element = W.index_set()
        elif hasattr(coxeter_element, "reduced_word"):
            coxeter_element = coxeter_element.reduced_word()
        coxeter_element = tuple(coxeter_element)

        return super(SubwordComplex, cls).__classcall__(cls, W=W, k=k,
                                                        coxeter_element=coxeter_element,
                                                        algorithm=algorithm)

    def __init__(self, W, k, coxeter_element, algorithm):
        """
        Initialize ``self``.

        TESTS::

            sage: S = ClusterComplex(['A', 2])
            sage: TestSuite(S).run()
            sage: S = ClusterComplex(['A', 2], k=2)
            sage: TestSuite(S).run()
        """
        w = W.w0
        Q = coxeter_element * k + tuple(w.coxeter_sorting_word(coxeter_element))
        SubwordComplex.__init__(self, Q, w, algorithm=algorithm)
        self._W = W
        self._w0 = w
        self._k = k
        if k == 1:
            self.__custom_name = 'Cluster complex'
        else:
            self.__custom_name = 'Multi-cluster complex'

        self.set_immutable()

    def __call__(self, F, facet_test=True):
        r"""
        Create a facet of ``self``.

        INPUT:

        - ``F`` -- an iterable of positions
        - ``facet_test`` -- boolean (default: ``True``); tells whether
          or not the facet ``F`` should be tested before creation

        OUTPUT:

        The facet of ``self`` at positions given by ``F``.

        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1)); F
            (0, 1)

        TESTS::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1))
            sage: C(F) is F
            True
        """
        if isinstance(F, ClusterComplexFacet) and F.parent() is self:
            return F
        return self.element_class(self, F, facet_test=facet_test)

    Element = ClusterComplexFacet

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ClusterComplex(['A', 2])._repr_()
            "Cluster complex of type ['A', 2] with 5 vertices and 5 facets"
        """
        name = self.__custom_name
        name += (' of type %s with %s vertices and %s facets'
                 % (self.cartan_type(), len(self.vertices()),
                    len(self._facets)))
        return name

    def k(self):
        r"""
        Return the index `k` of ``self``.

        EXAMPLES::

            sage: ClusterComplex(['A', 2]).k()
            1
        """
        return self._k

    def minimal_nonfaces(self):
        """
        Return the minimal non-faces of ``self``.

        EXAMPLES::

            sage: ClusterComplex(['A', 2]).minimal_nonfaces()
            [[0, 2], [0, 3], [1, 3], [1, 4], [2, 4]]
        """
        from sage.combinat.combination import Combinations
        return [X for X in Combinations(self.vertices(), self.k() + 1)
                if not any(set(X).issubset(F) for F in self.facets())]

    def cyclic_rotation(self):
        """
        Return the operation on the facets of ``self`` obtained by the
        cyclic rotation as defined in [CLS2014]_.

        EXAMPLES::

            sage: ClusterComplex(['A', 2]).cyclic_rotation()
            <function ...act at ...>
        """
        W = self._W
        w = self._w0
        Q = self._Q
        l = len(Q)
        S = W.simple_reflections()
        S_inv = {S[j]: j for j in W.index_set()}
        Q = Q + tuple(S_inv[w * S[k] * w] for k in Q)
        D = {i: (Q[i + 1:].index(Q[i]) + i + 1) % l for i in range(l)}

        def act(F):
            return self.parent().element_class(sorted([D[i] for i in F]))
        return act
