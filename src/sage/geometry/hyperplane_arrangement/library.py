r"""
Library of Hyperplane Arrangements

A collection of useful or interesting hyperplane arrangements. See
:mod:`sage.geometry.hyperplane_arrangement.arrangement` for details
about how to construct your own hyperplane arrangements.
"""
# ****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import binomial
from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
from sage.matrix.constructor import matrix, random_matrix
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.rational_field import QQ
from sage.rings.semirings.non_negative_integer_semiring import NN


def make_parent(base_ring, dimension, names=None):
    """
    Construct the parent for the hyperplane arrangements.

    For internal use only.

    INPUT:

    - ``base_ring`` -- a ring

    - ``dimension`` -- integer

    - ``names`` -- ``None`` (default) or a list/tuple/iterable of
      strings

    OUTPUT:

    A new
    :class:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangements`
    instance.

    EXAMPLES::

        sage: from sage.geometry.hyperplane_arrangement.library import make_parent
        sage: make_parent(QQ, 3)
        Hyperplane arrangements in 3-dimensional linear space over
        Rational Field with coordinates t0, t1, t2
    """
    if names is None:
        names = tuple('t'+str(i) for i in range(dimension))
    else:
        names = tuple(map(str, names))
        if len(names) != dimension:
            raise ValueError('number of variable names does not match dimension')
    return HyperplaneArrangements(base_ring, names=names)


class HyperplaneArrangementLibrary:
    """
    The library of hyperplane arrangements.
    """

    def braid(self, n, K=QQ, names=None):
        r"""
        The braid arrangement.

        INPUT:

        - ``n`` -- integer

        - ``K`` -- field (default: ``QQ``)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The hyperplane arrangement consisting of the `n(n-1)/2`
        hyperplanes `\{ x_i - x_j = 0 : 1 \leq i \leq j \leq n \}`.

        EXAMPLES::

            sage: hyperplane_arrangements.braid(4)                                      # needs sage.graphs
            Arrangement of 6 hyperplanes of dimension 4 and rank 3
        """
        from sage.graphs.graph_generators import graphs

        x = polygen(QQ, 'x')
        A = self.graphical(graphs.CompleteGraph(n), K, names=names)
        charpoly = prod(x - i for i in range(n))
        A.characteristic_polynomial.set_cache(charpoly)
        return A

    def bigraphical(self, G, A=None, K=QQ, names=None):
        r"""
        Return a bigraphical hyperplane arrangement.

        INPUT:

        - ``G`` -- graph

        - ``A`` -- list, matrix, dictionary (default: ``None``
          gives semiorder), or the string 'generic'

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The hyperplane arrangement with hyperplanes `x_i - x_j =
        A[i,j]` and `x_j - x_i = A[j,i]` for each edge `v_i, v_j` of
        ``G``.  The indices `i,j` are the indices of elements of
        ``G.vertices()``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: G = graphs.CycleGraph(4)
            sage: G.edges(sort=True)
            [(0, 1, None), (0, 3, None), (1, 2, None), (2, 3, None)]
            sage: G.edges(sort=True, labels=False)
            [(0, 1), (0, 3), (1, 2), (2, 3)]
            sage: A = {0:{1:1, 3:2}, 1:{0:3, 2:0}, 2:{1:2, 3:1}, 3:{2:0, 0:2}}
            sage: HA = hyperplane_arrangements.bigraphical(G, A)
            sage: HA.n_regions()
            63
            sage: hyperplane_arrangements.bigraphical(G, 'generic').n_regions()
            65
            sage: hyperplane_arrangements.bigraphical(G).n_regions()
            59

        REFERENCES:

        - [HP2016]_
        """
        n = G.num_verts()
        if A is None:  # default to G-semiorder arrangement
            A = matrix(K, n, lambda i, j: 1)
        elif A == 'generic':
            A = random_matrix(ZZ, n, x=10000)
            A = matrix(K, A)
        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        vertex_to_int = {u: i for i, u in enumerate(G)}
        for u, v in G.edge_iterator(labels=False, sort_vertices=False):
            i = vertex_to_int[u]
            j = vertex_to_int[v]
            hyperplanes.append(x[i] - x[j] - A[i][j])
            hyperplanes.append(-x[i] + x[j] - A[j][i])
        return H(*hyperplanes)

    def Catalan(self, n, K=QQ, names=None):
        r"""
        Return the Catalan arrangement.

        INPUT:

        - ``n`` -- integer

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The arrangement of `3n(n-1)/2` hyperplanes `\{ x_i - x_j =
        -1,0,1 : 1 \leq i \leq j \leq n \}`.

        EXAMPLES::

            sage: hyperplane_arrangements.Catalan(5)
            Arrangement of 30 hyperplanes of dimension 5 and rank 4

        TESTS::

            sage: h = hyperplane_arrangements.Catalan(5)
            sage: h.characteristic_polynomial()
            x^5 - 30*x^4 + 335*x^3 - 1650*x^2 + 3024*x
            sage: h.characteristic_polynomial.clear_cache()  # long time
            sage: h.characteristic_polynomial()              # long time
            x^5 - 30*x^4 + 335*x^3 - 1650*x^2 + 3024*x
        """
        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                for k in [-1, 0, 1]:
                    hyperplanes.append(x[i] - x[j] - k)
        Cn = H(*hyperplanes)
        x = polygen(QQ, 'x')
        charpoly = x*prod([x-n-i for i in range(1, n)])
        Cn.characteristic_polynomial.set_cache(charpoly)
        return Cn

    def coordinate(self, n, K=QQ, names=None):
        r"""
        Return the coordinate hyperplane arrangement.

        INPUT:

        - ``n`` -- integer

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The coordinate hyperplane arrangement, which is the central
        hyperplane arrangement consisting of the coordinate
        hyperplanes `x_i = 0`.

        EXAMPLES::

            sage: hyperplane_arrangements.coordinate(5)
            Arrangement of 5 hyperplanes of dimension 5 and rank 5
        """
        H = make_parent(K, n, names)
        x = H.gens()
        return H(x)

    def Coxeter(self, data, K=QQ, names=None):
        r"""
        Return the Coxeter arrangement.

        This generalizes the braid arrangements to crystallographic
        root systems.

        INPUT:

        - ``data`` -- either an integer or a Cartan type (or coercible
          into; see "CartanType")

        - ``K`` -- field (default: ``QQ``)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        - If ``data`` is an integer `n`, return the braid arrangement in
          dimension `n`, i.e. the set of `n(n-1)` hyperplanes:
          `\{ x_i - x_j = 0,1 : 1 \leq i \leq j \leq n \}`. This corresponds
          to the Coxeter arrangement of Cartan type `A_{n-1}`.

        - If ``data`` is a Cartan type, return the Coxeter arrangement of given
          type.

        The Coxeter arrangement of a given crystallographic
        Cartan type is defined by the inner products
        `\langle a,x \rangle = 0` where
        `a \in \Phi^+` runs over positive roots of the root system `\Phi`.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: hyperplane_arrangements.Coxeter(4)
            Arrangement of 6 hyperplanes of dimension 4 and rank 3
            sage: hyperplane_arrangements.Coxeter("B4")
            Arrangement of 16 hyperplanes of dimension 4 and rank 4
            sage: hyperplane_arrangements.Coxeter("A3")
            Arrangement of 6 hyperplanes of dimension 4 and rank 3

        If the Cartan type is not crystallographic, the Coxeter arrangement
        is not implemented yet::

            sage: hyperplane_arrangements.Coxeter("H3")                                 # needs sage.libs.gap
            Traceback (most recent call last):
            ...
            NotImplementedError: Coxeter arrangements are not implemented
            for non crystallographic Cartan types

        The characteristic polynomial is pre-computed using the results
        of Terao, see [Ath2000]_::

            sage: # needs sage.combinat
            sage: hyperplane_arrangements.Coxeter("A3").characteristic_polynomial()
            x^3 - 6*x^2 + 11*x - 6
        """
        from sage.combinat.root_system.cartan_type import CartanType
        from sage.combinat.root_system.root_system import RootSystem
        from sage.combinat.root_system.weyl_group import WeylGroup

        if data in NN:
            cartan_type = CartanType(["A", data - 1])
        else:
            cartan_type = CartanType(data)
        if not cartan_type.is_crystallographic():
            raise NotImplementedError("Coxeter arrangements are not implemented for non crystallographic Cartan types")
        W = WeylGroup(cartan_type)
        Ra = RootSystem(cartan_type).ambient_space()
        PR = Ra.positive_roots()
        d = Ra.dimension()
        H = make_parent(K, d, names)
        x = H.gens()
        hyperplanes = []

        for a in PR:
            hyperplanes.append(sum(a[j] * x[j] for j in range(d)))
        A = H(*hyperplanes)
        x = polygen(QQ, 'x')
        charpoly = prod(x - d + 1 for d in W.degrees())
        A.characteristic_polynomial.set_cache(charpoly)
        return A

    def G_semiorder(self, G, K=QQ, names=None):
        r"""
        Return the semiorder hyperplane arrangement of a graph.

        INPUT:

        - ``G`` -- graph

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The semiorder hyperplane arrangement of a graph G is the
        arrangement `\{ x_i - x_j = -1,1 \}` where `ij` is an edge of
        ``G``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.G_semiorder(G)
            Arrangement of 20 hyperplanes of dimension 5 and rank 4
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.G_semiorder(g)
            Arrangement of 12 hyperplanes of dimension 5 and rank 4
        """
        n = G.num_verts()
        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        vertex_to_int = {u: i for i, u in enumerate(G.vertices(sort=True))}
        for u, v in G.edge_iterator(labels=False):
            i = vertex_to_int[u]
            j = vertex_to_int[v]
            hyperplanes.append(x[i] - x[j] - 1)
            hyperplanes.append(x[i] - x[j] + 1)
        return H(*hyperplanes)

    def G_Shi(self, G, K=QQ, names=None):
        r"""
        Return the Shi hyperplane arrangement of a graph `G`.

        INPUT:

        - ``G`` -- graph

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT: the Shi hyperplane arrangement of the given graph ``G``

        EXAMPLES::

            sage: # needs sage.graphs
            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.G_Shi(G)
            Arrangement of 20 hyperplanes of dimension 5 and rank 4
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.G_Shi(g)
            Arrangement of 12 hyperplanes of dimension 5 and rank 4
            sage: a = hyperplane_arrangements.G_Shi(graphs.WheelGraph(4)); a
            Arrangement of 12 hyperplanes of dimension 4 and rank 3
        """
        n = G.num_verts()
        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        vertex_to_int = {u: i for i, u in enumerate(G.vertices(sort=True))}
        for u, v in G.edge_iterator(labels=False):
            i = vertex_to_int[u]
            j = vertex_to_int[v]
            hyperplanes.append(x[i] - x[j])
            hyperplanes.append(x[i] - x[j] - 1)
        return H(*hyperplanes)

    def graphical(self, G, K=QQ, names=None):
        r"""
        Return the graphical hyperplane arrangement of a graph ``G``.

        INPUT:

        - ``G`` -- graph

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The graphical hyperplane arrangement of a graph G, which is
        the arrangement `\{ x_i - x_j = 0 \}` for all edges `ij` of the
        graph ``G``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.graphical(G)
            Arrangement of 10 hyperplanes of dimension 5 and rank 4
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.graphical(g)
            Arrangement of 6 hyperplanes of dimension 5 and rank 4

        TESTS::

            sage: # needs sage.graphs
            sage: h = hyperplane_arrangements.graphical(g)
            sage: h.characteristic_polynomial()
            x^5 - 6*x^4 + 14*x^3 - 15*x^2 + 6*x
            sage: h.characteristic_polynomial.clear_cache()     # long time
            sage: h.characteristic_polynomial()         # long time
            x^5 - 6*x^4 + 14*x^3 - 15*x^2 + 6*x
        """
        n = G.num_verts()
        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        vertex_to_int = {u: i for i, u in enumerate(G.vertices(sort=True))}
        for u, v in G.edge_iterator(labels=False):
            i = vertex_to_int[u]
            j = vertex_to_int[v]
            hyperplanes.append(x[i] - x[j])
        A = H(*hyperplanes)
        charpoly = G.chromatic_polynomial()
        A.characteristic_polynomial.set_cache(charpoly)
        return A

    def Ish(self, n, K=QQ, names=None):
        r"""
        Return the Ish arrangement.

        INPUT:

        - ``n`` -- integer

        - ``K`` -- field (default: ``QQ``)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The Ish arrangement, which is the set of `n(n-1)` hyperplanes.

        .. MATH::

            \{ x_i - x_j = 0 : 1 \leq i \leq j \leq n \}
            \cup
            \{ x_1 - x_j = i : 1 \leq i \leq j \leq n \}.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: a = hyperplane_arrangements.Ish(3); a
            Arrangement of 6 hyperplanes of dimension 3 and rank 2
            sage: a.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x
            sage: b = hyperplane_arrangements.Shi(3)
            sage: b.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x

        TESTS::

            sage: a.characteristic_polynomial.clear_cache()     # long time             # needs sage.combinat
            sage: a.characteristic_polynomial()         # long time                     # needs sage.combinat
            x^3 - 6*x^2 + 9*x

        REFERENCES:

        - [AR2012]_
        """
        from sage.combinat.combinat import stirling_number2

        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                hyperplanes.append(x[i] - x[j])
                hyperplanes.append(x[0] - x[j] - (i+1))
        A = H(*hyperplanes)
        x = polygen(QQ, 'x')
        charpoly = x * sum([(-1)**k * stirling_number2(n, n-k) *
                            prod([(x - 1 - j) for j in range(k, n-1)])
                            for k in range(n)])
        A.characteristic_polynomial.set_cache(charpoly)
        return A

    def IshB(self, n, K=QQ, names=None):
        r"""
        Return the type B Ish arrangement.

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default: ``QQ``)
        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT: the type `B` Ish arrangement, which is the set of `2n^2` hyperplanes

        .. MATH::

            \{ x_i \pm x_j = 0 : 1 \leq i < j \leq n \}
            \cup
            \{ x_i = a : 1 \leq i\leq n, \quad i - n \leq a \leq n - i + 1 \}.

        EXAMPLES::

            sage: a = hyperplane_arrangements.IshB(2)
            sage: a
            Arrangement of 8 hyperplanes of dimension 2 and rank 2
            sage: a.hyperplanes()
            (Hyperplane 0*t0 + t1 - 1,
             Hyperplane 0*t0 + t1 + 0,
             Hyperplane t0 - t1 + 0,
             Hyperplane t0 + 0*t1 - 2,
             Hyperplane t0 + 0*t1 - 1,
             Hyperplane t0 + 0*t1 + 0,
             Hyperplane t0 + 0*t1 + 1,
             Hyperplane t0 + t1 + 0)
            sage: a.cone().is_free()                                                    # needs sage.libs.singular
            True

        .. PLOT::
            :width: 500 px

            a = hyperplane_arrangements.IshB(2)
            sphinx_plot(a.plot())

        ::

            sage: a = hyperplane_arrangements.IshB(3); a
            Arrangement of 18 hyperplanes of dimension 3 and rank 3
            sage: a.characteristic_polynomial()
            x^3 - 18*x^2 + 108*x - 216
            sage: b = hyperplane_arrangements.Shi(['B', 3])
            sage: b.characteristic_polynomial()
            x^3 - 18*x^2 + 108*x - 216

        TESTS::

            sage: a.characteristic_polynomial.clear_cache()
            sage: a.characteristic_polynomial()
            x^3 - 18*x^2 + 108*x - 216

        REFERENCES:

        - [TT2023]_
        """
        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                hyperplanes.append(x[i] - x[j])
                hyperplanes.append(x[i] + x[j])
            for a in range(i+1-n, n-i+1):
                hyperplanes.append(x[i] - a)
        A = H(*hyperplanes)
        x = polygen(QQ, 'x')
        charpoly = (x - 2*n) ** n
        A.characteristic_polynomial.set_cache(charpoly)
        return A

    def linial(self, n, K=QQ, names=None):
        r"""
        Return the linial hyperplane arrangement.

        INPUT:

        - ``n`` -- integer

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The linial hyperplane arrangement is the set of hyperplanes
        `\{x_i - x_j = 1 : 1\leq i < j \leq n\}`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.linial(4);  a
            Arrangement of 6 hyperplanes of dimension 4 and rank 3
            sage: a.characteristic_polynomial()
            x^4 - 6*x^3 + 15*x^2 - 14*x

        TESTS::

            sage: h = hyperplane_arrangements.linial(5)
            sage: h.characteristic_polynomial()
            x^5 - 10*x^4 + 45*x^3 - 100*x^2 + 90*x
            sage: h.characteristic_polynomial.clear_cache()  # long time
            sage: h.characteristic_polynomial()              # long time
            x^5 - 10*x^4 + 45*x^3 - 100*x^2 + 90*x
        """
        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                hyperplanes.append(x[i] - x[j] - 1)
        A = H(*hyperplanes)
        x = polygen(QQ, 'x')
        charpoly = x * sum(binomial(n, k)*(x - k)**(n - 1) for k in range(n + 1)) / 2**n
        A.characteristic_polynomial.set_cache(charpoly)
        return A

    def semiorder(self, n, K=QQ, names=None):
        r"""
        Return the semiorder arrangement.

        INPUT:

        - ``n`` -- integer

        - ``K`` -- field (default: `\QQ`)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        OUTPUT:

        The semiorder arrangement, which is the set of `n(n-1)`
        hyperplanes `\{ x_i - x_j = -1,1 : 1 \leq i \leq j \leq n\}`.

        EXAMPLES::

            sage: hyperplane_arrangements.semiorder(4)
            Arrangement of 12 hyperplanes of dimension 4 and rank 3

        TESTS::

            sage: # needs sage.combinat
            sage: h = hyperplane_arrangements.semiorder(5)
            sage: h.characteristic_polynomial()
            x^5 - 20*x^4 + 180*x^3 - 790*x^2 + 1380*x
            sage: h.characteristic_polynomial.clear_cache()     # long time
            sage: h.characteristic_polynomial()         # long time
            x^5 - 20*x^4 + 180*x^3 - 790*x^2 + 1380*x
        """
        from sage.combinat.combinat import stirling_number2

        H = make_parent(K, n, names)
        x = H.gens()
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                for k in [-1, 1]:
                    hyperplanes.append(x[i] - x[j] - k)
        A = H(*hyperplanes)
        x = polygen(QQ, 'x')
        charpoly = x * sum([stirling_number2(n, k) * prod([x - k - i for i in range(1, k)])
                            for k in range(1, n+1)])
        A.characteristic_polynomial.set_cache(charpoly)
        return A

    def Shi(self, data, K=QQ, names=None, m=1):
        r"""
        Return the Shi arrangement.

        INPUT:

        - ``data`` -- either an integer or a Cartan type (or coercible
          into; see "CartanType")

        - ``K`` -- field (default: ``QQ``)

        - ``names`` -- tuple of strings or ``None`` (default); the
          variable names for the ambient space

        - ``m`` -- integer (default: 1)

        OUTPUT:

        - If ``data`` is an integer `n`, return the Shi arrangement in
          dimension `n`, i.e. the set of `n(n-1)` hyperplanes:
          `\{ x_i - x_j = 0,1 : 1 \leq i \leq j \leq n \}`. This corresponds
          to the Shi arrangement of Cartan type `A_{n-1}`.

        - If ``data`` is a Cartan type, return the Shi arrangement of given
          type.

        - If `m > 1`, return the `m`-extended Shi arrangement of given type.

        The `m`-extended Shi arrangement of a given crystallographic
        Cartan type is defined by the inner product
        `\langle a,x \rangle = k` for `-m < k \leq m` and
        `a \in \Phi^+` is a positive root of the root system `\Phi`.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: hyperplane_arrangements.Shi(4)
            Arrangement of 12 hyperplanes of dimension 4 and rank 3
            sage: hyperplane_arrangements.Shi("A3")
            Arrangement of 12 hyperplanes of dimension 4 and rank 3
            sage: hyperplane_arrangements.Shi("A3", m=2)
            Arrangement of 24 hyperplanes of dimension 4 and rank 3
            sage: hyperplane_arrangements.Shi("B4")
            Arrangement of 32 hyperplanes of dimension 4 and rank 4
            sage: hyperplane_arrangements.Shi("B4", m=3)
            Arrangement of 96 hyperplanes of dimension 4 and rank 4
            sage: hyperplane_arrangements.Shi("C3")
            Arrangement of 18 hyperplanes of dimension 3 and rank 3
            sage: hyperplane_arrangements.Shi("D4", m=3)
            Arrangement of 72 hyperplanes of dimension 4 and rank 4
            sage: hyperplane_arrangements.Shi("E6")
            Arrangement of 72 hyperplanes of dimension 8 and rank 6
            sage: hyperplane_arrangements.Shi("E6", m=2)
            Arrangement of 144 hyperplanes of dimension 8 and rank 6

        If the Cartan type is not crystallographic, the Shi arrangement
        is not defined::

            sage: hyperplane_arrangements.Shi("H4")
            Traceback (most recent call last):
            ...
            NotImplementedError: Shi arrangements are not defined for non crystallographic Cartan types

        The characteristic polynomial is pre-computed using the results
        of [Ath1996]_::

            sage: # needs sage.combinat
            sage: hyperplane_arrangements.Shi("A3").characteristic_polynomial()
            x^4 - 12*x^3 + 48*x^2 - 64*x
            sage: hyperplane_arrangements.Shi("A3", m=2).characteristic_polynomial()
            x^4 - 24*x^3 + 192*x^2 - 512*x
            sage: hyperplane_arrangements.Shi("C3").characteristic_polynomial()
            x^3 - 18*x^2 + 108*x - 216
            sage: hyperplane_arrangements.Shi("E6").characteristic_polynomial()
            x^8 - 72*x^7 + 2160*x^6 - 34560*x^5 + 311040*x^4 - 1492992*x^3 + 2985984*x^2
            sage: hyperplane_arrangements.Shi("B4", m=3).characteristic_polynomial()
            x^4 - 96*x^3 + 3456*x^2 - 55296*x + 331776

        TESTS::

            sage: # needs sage.combinat
            sage: h = hyperplane_arrangements.Shi(4)
            sage: h.characteristic_polynomial()
            x^4 - 12*x^3 + 48*x^2 - 64*x
            sage: h.characteristic_polynomial.clear_cache()     # long time
            sage: h.characteristic_polynomial()         # long time
            x^4 - 12*x^3 + 48*x^2 - 64*x
            sage: h = hyperplane_arrangements.Shi("A3", m=2)
            sage: h.characteristic_polynomial()
            x^4 - 24*x^3 + 192*x^2 - 512*x
            sage: h.characteristic_polynomial.clear_cache()
            sage: h.characteristic_polynomial()
            x^4 - 24*x^3 + 192*x^2 - 512*x
            sage: h = hyperplane_arrangements.Shi("B3", m=3)
            sage: h.characteristic_polynomial()
            x^3 - 54*x^2 + 972*x - 5832
            sage: h.characteristic_polynomial.clear_cache()
            sage: h.characteristic_polynomial()
            x^3 - 54*x^2 + 972*x - 5832
        """
        from sage.combinat.root_system.cartan_type import CartanType
        from sage.combinat.root_system.root_system import RootSystem

        if data in NN:
            cartan_type = CartanType(["A", data - 1])
        else:
            cartan_type = CartanType(data)
        if not cartan_type.is_crystallographic():
            raise NotImplementedError("Shi arrangements are not defined for non crystallographic Cartan types")
        n = cartan_type.rank()
        h = cartan_type.coxeter_number()
        Ra = RootSystem(cartan_type).ambient_space()
        PR = Ra.positive_roots()
        d = Ra.dimension()
        H = make_parent(K, d, names)
        x = H.gens()
        hyperplanes = []

        for a in PR:
            for const in range(-m + 1, m + 1):
                hyperplanes.append(sum(a[j]*x[j] for j in range(d))-const)
        A = H(*hyperplanes)
        x = polygen(QQ, 'x')
        charpoly = x**(d-n) * (x-m*h)**n
        A.characteristic_polynomial.set_cache(charpoly)
        return A


hyperplane_arrangements = HyperplaneArrangementLibrary()
