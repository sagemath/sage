# sage_setup: distribution = sagemath-modules
"""
Root system data for super type A
"""
# ****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from . import ambient_space
from .cartan_type import SuperCartanType_standard


class AmbientSpace(ambient_space.AmbientSpace):
    r"""
    The ambient space for (super) type `A(m|n)`.

    EXAMPLES::

        sage: R = RootSystem(['A', [2,1]])
        sage: AL = R.ambient_space(); AL
        Ambient space of the Root system of type ['A', [2, 1]]
        sage: AL.basis()
        Finite family {-3: (1, 0, 0, 0, 0),
         -2: (0, 1, 0, 0, 0),
         -1: (0, 0, 1, 0, 0),
         1: (0, 0, 0, 1, 0),
         2: (0, 0, 0, 0, 1)}
    """

    def __init__(self, root_system, base_ring, index_set=None):
        """
        Initialize ``self``.

        TESTS::

            sage: R = RootSystem(['A', [4,2]])
            sage: AL = R.ambient_space(); AL
            Ambient space of the Root system of type ['A', [4, 2]]
            sage: TestSuite(AL).run(skip='_test_norm_of_simple_roots')
        """
        ct = root_system.cartan_type()
        if index_set is None:
            index_set = tuple(list(range(-ct.m - 1, 0)) +
                              list(range(1, ct.n + 2)))
        ambient_space.AmbientSpace.__init__(self, root_system, base_ring,
                                            index_set=index_set)

    @classmethod
    def smallest_base_ring(cls, cartan_type=None):
        """
        Return the smallest base ring the ambient space can be defined upon.

        .. SEEALSO::

            :meth:`~sage.combinat.root_system.ambient_space.AmbientSpace.smallest_base_ring`

        EXAMPLES::

            sage: e = RootSystem(['A', [3,1]]).ambient_space()
            sage: e.smallest_base_ring()
            Integer Ring
        """
        return ZZ

    def dimension(self):
        """
        Return the dimension of this ambient space.

        EXAMPLES::

            sage: e = RootSystem(['A', [4,2]]).ambient_space()
            sage: e.dimension()
            8
        """
        ct = self.root_system.cartan_type()
        return ct.m + ct.n + 2

    def simple_root(self, i):
        """
        Return the `i`-th simple root of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: list(e.simple_roots())
            [(1, -1, 0, 0, 0), (0, 1, -1, 0, 0),
             (0, 0, 1, -1, 0), (0, 0, 0, 1, -1)]
        """
        if i < 0:
            return self.monomial(i-1) - self.monomial(i)
        if i == 0:
            return self.monomial(-1) - self.monomial(1)
        return self.monomial(i) - self.monomial(i+1)

    def positive_roots(self):
        """
        Return the positive roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.positive_roots()
            [(0, 1, -1, 0, 0),
             (1, 0, -1, 0, 0),
             (1, -1, 0, 0, 0),
             (0, 0, 0, 1, -1),
             (0, 0, 1, -1, 0),
             (0, 0, 1, 0, -1),
             (0, 1, 0, -1, 0),
             (0, 1, 0, 0, -1),
             (1, 0, 0, -1, 0),
             (1, 0, 0, 0, -1)]
        """
        return self.positive_even_roots() + self.positive_odd_roots()

    def positive_even_roots(self):
        """
        Return the positive even roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.positive_even_roots()
            [(0, 1, -1, 0, 0), (1, 0, -1, 0, 0),
             (1, -1, 0, 0, 0), (0, 0, 0, 1, -1)]
        """
        ct = self.root_system.cartan_type()
        ret = []
        ret += [self.monomial(-j) - self.monomial(-i)
                for i in range(1, ct.m + 2)
                for j in range(i + 1, ct.m + 2)]
        ret += [self.monomial(i) - self.monomial(j)
                for i in range(1, ct.n + 2)
                for j in range(i + 1, ct.n + 2)]
        return ret

    def positive_odd_roots(self):
        """
        Return the positive odd roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.positive_odd_roots()
            [(0, 0, 1, -1, 0),
             (0, 0, 1, 0, -1),
             (0, 1, 0, -1, 0),
             (0, 1, 0, 0, -1),
             (1, 0, 0, -1, 0),
             (1, 0, 0, 0, -1)]
        """
        ct = self.root_system.cartan_type()
        return [self.monomial(-i) - self.monomial(j)
                for i in range(1, ct.m + 2)
                for j in range(1, ct.n + 2)]

    def highest_root(self):
        """
        Return the highest root of ``self``.

        EXAMPLES::

           sage: e = RootSystem(['A', [4,2]]).ambient_lattice()
           sage: e.highest_root()
           (1, 0, 0, 0, 0, 0, 0, -1)
        """
        ct = self.root_system.cartan_type()
        return self.monomial(-ct.m-1) - self.monomial(ct.n+1)

    def negative_roots(self):
        """
        Return the negative roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.negative_roots()
            [(0, -1, 1, 0, 0),
             (-1, 0, 1, 0, 0),
             (-1, 1, 0, 0, 0),
             (0, 0, 0, -1, 1),
             (0, 0, -1, 1, 0),
             (0, 0, -1, 0, 1),
             (0, -1, 0, 1, 0),
             (0, -1, 0, 0, 1),
             (-1, 0, 0, 1, 0),
             (-1, 0, 0, 0, 1)]
        """
        return self.negative_even_roots() + self.negative_odd_roots()

    def negative_even_roots(self):
        """
        Return the negative even roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.negative_even_roots()
            [(0, -1, 1, 0, 0), (-1, 0, 1, 0, 0),
             (-1, 1, 0, 0, 0), (0, 0, 0, -1, 1)]
        """
        ct = self.root_system.cartan_type()
        ret = []
        ret += [self.monomial(-i) - self.monomial(-j)
                for i in range(1, ct.m + 2)
                for j in range(i + 1, ct.m + 2)]
        ret += [self.monomial(j) - self.monomial(i)
                for i in range(1, ct.n + 2)
                for j in range(i + 1, ct.n + 2)]
        return ret

    def negative_odd_roots(self):
        """
        Return the negative odd roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.negative_odd_roots()
            [(0, 0, -1, 1, 0),
             (0, 0, -1, 0, 1),
             (0, -1, 0, 1, 0),
             (0, -1, 0, 0, 1),
             (-1, 0, 0, 1, 0),
             (-1, 0, 0, 0, 1)]
        """
        ct = self.root_system.cartan_type()
        return [self.monomial(j) - self.monomial(-i)
                for i in range(1, ct.m + 2)
                for j in range(1, ct.n + 2)]

    def fundamental_weight(self, i):
        r"""
        Return the fundamental weight `\Lambda_i` of ``self``.

        EXAMPLES::

            sage: L = RootSystem(['A', [3,2]]).ambient_space()
            sage: L.fundamental_weight(-1)
            (1, 1, 1, 0, 0, 0, 0)
            sage: L.fundamental_weight(0)
            (1, 1, 1, 1, 0, 0, 0)
            sage: L.fundamental_weight(2)
            (1, 1, 1, 1, -1, -1, -2)
            sage: list(L.fundamental_weights())
            [(1, 0, 0, 0, 0, 0, 0),
             (1, 1, 0, 0, 0, 0, 0),
             (1, 1, 1, 0, 0, 0, 0),
             (1, 1, 1, 1, 0, 0, 0),
             (1, 1, 1, 1, -1, -2, -2),
             (1, 1, 1, 1, -1, -1, -2)]

        ::

            sage: L = RootSystem(['A', [2,3]]).ambient_space()
            sage: La = L.fundamental_weights()
            sage: al = L.simple_roots()
            sage: I = L.index_set()
            sage: matrix([[al[i].scalar(La[j]) for i in I] for j in I])
            [ 1  0  0  0  0  0]
            [ 0  1  0  0  0  0]
            [ 0  0  1  0  0  0]
            [ 0  0  0 -1  0  0]
            [ 0  0  0  0 -1  0]
            [ 0  0  0  0  0 -1]
        """
        m = self.root_system.cartan_type().m
        n = self.root_system.cartan_type().n
        if i <= 0:
            return self.sum(self.monomial(j) for j in range(-m-1,i))
        return (self.sum(self.monomial(j) for j in range(-m-1,1))
                - self.sum(self.monomial(j) for j in range(i+1))
                - 2*self.sum(self.monomial(j) for j in range(i+1,n+2)))

    def simple_coroot(self, i):
        """
        Return the simple coroot `h_i` of ``self``.

        EXAMPLES::

            sage: L = RootSystem(['A', [3,2]]).ambient_space()
            sage: L.simple_coroot(-2)
            (0, 1, -1, 0, 0, 0, 0)
            sage: L.simple_coroot(0)
            (0, 0, 0, 1, -1, 0, 0)
            sage: L.simple_coroot(2)
            (0, 0, 0, 0, 0, -1, 1)
            sage: list(L.simple_coroots())
            [(1, -1, 0, 0, 0, 0, 0),
             (0, 1, -1, 0, 0, 0, 0),
             (0, 0, 1, -1, 0, 0, 0),
             (0, 0, 0, 1, -1, 0, 0),
             (0, 0, 0, 0, -1, 1, 0),
             (0, 0, 0, 0, 0, -1, 1)]
        """
        if i <= 0:
            return self.simple_root(i)
        return -self.simple_root(i)

    class Element(ambient_space.AmbientSpaceElement):
        def inner_product(self, lambdacheck):
            """
            The scalar product with elements of the coroot lattice
            embedded in the ambient space.

            EXAMPLES::

                sage: L = RootSystem(['A', [2,1]]).ambient_space()
                sage: a = L.simple_roots()
                sage: matrix([[a[i].inner_product(a[j]) for j in L.index_set()] for i in L.index_set()])
                [ 2 -1  0  0]
                [-1  2 -1  0]
                [ 0 -1  0  1]
                [ 0  0  1 -2]
            """
            self_mc = self._monomial_coefficients
            lambdacheck_mc = lambdacheck._monomial_coefficients

            result = self.parent().base_ring().zero()
            for t,c in lambdacheck_mc.items():
                if t not in self_mc:
                    continue
                if t > 0:
                    result -= c*self_mc[t]
                else:
                    result += c*self_mc[t]
            return result

        scalar = inner_product
        dot_product = inner_product

        def associated_coroot(self):
            """
            Return the coroot associated to ``self``.

            EXAMPLES::

                sage: L = RootSystem(['A', [3,2]]).ambient_space()
                sage: al = L.simple_roots()
                sage: al[-1].associated_coroot()
                (0, 0, 1, -1, 0, 0, 0)
                sage: al[0].associated_coroot()
                (0, 0, 0, 1, -1, 0, 0)
                sage: al[1].associated_coroot()
                (0, 0, 0, 0, -1, 1, 0)

                sage: a = al[-1] + al[0] + al[1]; a
                (0, 0, 1, 0, 0, -1, 0)
                sage: a.associated_coroot()
                (0, 0, 1, 0, -2, 1, 0)
                sage: h = L.simple_coroots()
                sage: h[-1] + h[0] + h[1]
                (0, 0, 1, 0, -2, 1, 0)

                sage: (al[-1] + al[0] + al[2]).associated_coroot()
                (0, 0, 1, 0, -1, -1, 1)
            """
            P = self.parent()
            al = P.simple_roots()
            h = P.simple_coroots()
            try:
                return h[al.inverse_family()[self]]
            except KeyError:
                pass
            V = P._dense_free_module()
            dep = V.linear_dependence([self._vector_()] +
                                      [al[i]._vector_() for i in P.index_set()])[0]
            I = P.index_set()
            return P.sum((-c/dep[0]) * h[I[i]] for i,c in dep[1:].items())

        def has_descent(self, i, positive=False):
            """
            Test if ``self`` has a descent at position `i`, that is
            if ``self`` is on the strict negative side of the `i`-th
            simple reflection hyperplane.

            If ``positive`` is ``True``, tests if it is on the strict
            positive side instead.

            EXAMPLES::

                sage: L = RootSystem(['A', [2,1]]).ambient_space()
                sage: al = L.simple_roots()
                sage: [al[i].has_descent(1) for i in L.index_set()]
                [False, False, True, False]
                sage: [(-al[i]).has_descent(1) for i in L.index_set()]
                [False, False, False, True]
                sage: [al[i].has_descent(1, True) for i in L.index_set()]
                [False, False, False, True]
                sage: [(-al[i]).has_descent(1, True) for i in L.index_set()]
                [False, False, True, False]
                sage: (al[-2] + al[0] + al[1]).has_descent(-1)
                True
                sage: (al[-2] + al[0] + al[1]).has_descent(1)
                False
                sage: (al[-2] + al[0] + al[1]).has_descent(1, positive=True)
                True
                sage: all(all(not la.has_descent(i) for i in L.index_set())
                ....:     for la in L.fundamental_weights())
                True
            """
            s = self.scalar(self.parent().simple_roots()[i])
            if i > 0:
                s = -s
            if positive:
                return s > 0
            else:
                return s < 0

        def is_dominant_weight(self):
            """
            Test whether ``self`` is a dominant element of the weight lattice.

            EXAMPLES::

                sage: L = RootSystem(['A',2]).ambient_lattice()
                sage: Lambda = L.fundamental_weights()
                sage: [x.is_dominant() for x in Lambda]
                [True, True]
                sage: (3*Lambda[1]+Lambda[2]).is_dominant()
                True
                sage: (Lambda[1]-Lambda[2]).is_dominant()
                False
                sage: (-Lambda[1]+Lambda[2]).is_dominant()
                False

            Tests that the scalar products with the coroots are all
            nonnegative integers. For example, if `x` is the sum of a
            dominant element of the weight lattice plus some other element
            orthogonal to all coroots, then the implementation correctly
            reports `x` to be a dominant weight::

               sage: x = Lambda[1] + L([-1,-1,-1])
               sage: x.is_dominant_weight()
               True
            """
            alpha = self.parent().simple_roots()
            l = self.parent().cartan_type().symmetrizer()
            from sage.rings.semirings.non_negative_integer_semiring import NN
            return all(l[i] * self.inner_product(alpha[i]) in NN
                       for i in self.parent().index_set())


class CartanType(SuperCartanType_standard):
    """
    Cartan Type `A(m|n)`.

    .. SEEALSO:: :func:`~sage.combinat.root_systems.cartan_type.CartanType`
    """

    def __init__(self, m, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A', [4,2]])
            sage: ct
            ['A', [4, 2]]
            sage: ct._repr_(compact=True)
            'A4|2'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.affine() # Not tested -- to be implemented
            ['A', [4, 2], 1]
            sage: ct.dual()
            ['A', [4, 2]]

        TESTS::

            sage: TestSuite(ct).run()
        """
        self.m = m
        self.n = n
        self.letter = "A"

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['A',[4,3]]))
            A(4|3)
        """
        return "A(%s|%s)" % (self.m, self.n)

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: CartanType(['A', [2,3]]).index_set()
            (-2, -1, 0, 1, 2, 3)
        """
        return tuple(range(-self.m, self.n + 1))

    AmbientSpace = AmbientSpace

    def is_irreducible(self):
        """
        Return whether ``self`` is irreducible, which is ``True``.

        EXAMPLES::

            sage: CartanType(['A', [3,4]]).is_irreducible()
            True
        """
        return True

    # A lot of these methods should be implemented by the ABCs of CartanType

    def is_affine(self):
        """
        Return whether ``self`` is affine or not.

        EXAMPLES::

            sage: CartanType(['A', [2,3]]).is_affine()
            False
        """
        return False

    def is_finite(self):
        """
        Return whether ``self`` is finite or not.

        EXAMPLES::

            sage: CartanType(['A', [2,3]]).is_finite()
            True
        """
        return True

    def dual(self):
        """
        Return dual of ``self``.

        EXAMPLES::

            sage: CartanType(['A', [2,3]]).dual()
            ['A', [2, 3]]
        """
        return self

    def type(self):
        """
        Return type of ``self``.

        EXAMPLES::

            sage: CartanType(['A', [2,3]]).type()
            'A'
        """
        return 'A'

    def root_system(self):
        """
        Return root system of ``self``.

        EXAMPLES::

            sage: CartanType(['A', [2,3]]).root_system()
            Root system of type ['A', [2, 3]]
        """
        from sage.combinat.root_system.root_system import RootSystem
        return RootSystem(self)

    @cached_method
    def symmetrizer(self):
        """
        Return symmetrizing matrix for ``self``.

        EXAMPLES::

            sage: CartanType(['A', [2,3]]).symmetrizer()
            Finite family {-2: 1, -1: 1, 0: 1, 1: -1, 2: -1, 3: -1}
        """
        from sage.sets.family import Family

        def ell(i):
            return ZZ.one() if i <= 0 else -ZZ.one()
        return Family(self.index_set(), ell)

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram of super type A.

        EXAMPLES::

            sage: a = CartanType(['A', [4,2]]).dynkin_diagram(); a                      # needs sage.graphs
            O---O---O---O---X---O---O
            -4  -3  -2  -1  0   1   2
            A4|2
            sage: a.edges(sort=True)                                                    # needs sage.graphs
            [(-4, -3, 1), (-3, -4, 1), (-3, -2, 1), (-2, -3, 1),
             (-2, -1, 1), (-1, -2, 1), (-1, 0, 1), (0, -1, 1),
             (0, 1, 1), (1, 0, -1), (1, 2, 1), (2, 1, 1)]

        TESTS::

            sage: a = DynkinDiagram(['A', [0,0]]); a                                    # needs sage.graphs
            X
            0
            A0|0
            sage: a.vertices(sort=False), a.edges(sort=False)                           # needs sage.graphs
            ([0], [])

            sage: a = DynkinDiagram(['A', [1,0]]); a                                    # needs sage.graphs
            O---X
            -1  0
            A1|0
            sage: a.vertices(sort=True), a.edges(sort=True)                             # needs sage.graphs
            ([-1, 0], [(-1, 0, 1), (0, -1, 1)])

            sage: a = DynkinDiagram(['A', [0,1]]); a                                    # needs sage.graphs
            X---O
            0   1
            A0|1
            sage: a.vertices(sort=True), a.edges(sort=True)                             # needs sage.graphs
            ([0, 1], [(0, 1, 1), (1, 0, -1)])
        """
        from .dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self, odd_isotropic_roots=[0])
        for i in range(self.m):
            g.add_edge(-i-1, -i)
        for i in range(1, self.n):
            g.add_edge(i, i+1)
        g.add_vertex(0)  # Usually there, but not when m == n == 0
        if self.m > 0:
            g.add_edge(-1, 0)
        if self.n > 0:
            g.add_edge(1, 0, -1)
        return g

    def cartan_matrix(self):
        """
        Return the Cartan matrix associated to ``self``.

        EXAMPLES::

            sage: ct = CartanType(['A', [2,3]])
            sage: ct.cartan_matrix()                                                    # needs sage.graphs
            [ 2 -1  0  0  0  0]
            [-1  2 -1  0  0  0]
            [ 0 -1  0  1  0  0]
            [ 0  0 -1  2 -1  0]
            [ 0  0  0 -1  2 -1]
            [ 0  0  0  0 -1  2]

        TESTS::

            sage: ct = CartanType(['A', [0,0]])
            sage: ct.cartan_matrix()                                                    # needs sage.graphs
            [0]

            sage: ct = CartanType(['A', [1,0]])
            sage: ct.cartan_matrix()                                                    # needs sage.graphs
            [ 2 -1]
            [-1  0]

            sage: ct = CartanType(['A', [0,1]])
            sage: ct.cartan_matrix()                                                    # needs sage.graphs
            [ 0  1]
            [-1  2]
        """
        return self.dynkin_diagram().cartan_matrix()

    def relabel(self, relabelling):
        """
        Return a relabelled copy of this Cartan type.

        INPUT:

        - ``relabelling`` -- a function (or a list or dictionary)

        OUTPUT:

        an isomorphic Cartan type obtained by relabelling the nodes of
        the Dynkin diagram. Namely, the node with label ``i`` is
        relabelled ``f(i)`` (or, by ``f[i]`` if ``f`` is a list or
        dictionary).

        EXAMPLES::

            sage: ct = CartanType(['A', [1,2]])
            sage: ct.dynkin_diagram()                                                   # needs sage.graphs
            O---X---O---O
            -1  0   1   2
            A1|2
            sage: f = {1:2, 2:1, 0:0, -1:-1}
            sage: ct.relabel(f)
            ['A', [1, 2]] relabelled by {-1: -1, 0: 0, 1: 2, 2: 1}
            sage: ct.relabel(f).dynkin_diagram()                                        # needs sage.graphs
            O---X---O---O
            -1  0   2   1
            A1|2 relabelled by {-1: -1, 0: 0, 1: 2, 2: 1}
        """
        from . import type_relabel
        return type_relabel.CartanType(self, relabelling)

    def _latex_draw_node(self, x, y, label, position='below=4pt'):
        r"""
        Draw (possibly marked [crossed out]) circular node ``i`` at the
        position ``(x,y)`` with node label ``label`` .

        - ``position`` -- position of the label relative to the node
        - ``anchor`` -- (optional) the anchor point for the label

        EXAMPLES::

            sage: t = CartanType(['A', [3,2]])
            sage: print(t._latex_draw_node(0, 0, 0))
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            \draw[-,thick] (0.17 cm, 0.17 cm) -- (-0.17 cm, -0.17 cm);
            \draw[-,thick] (0.17 cm, -0.17 cm) -- (-0.17 cm, 0.17 cm);
            sage: print(t._latex_draw_node(0, 0, 1))
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
        """
        ret = "\\draw[fill={}] ({} cm, {} cm) circle (.25cm) node[{}]{{${}$}};\n".format(
              'white', x, y, position, label)
        if label == 0:
            ret += "\\draw[-,thick] ({} cm, {} cm) -- ({} cm, {} cm);\n".format(
                                    x+.17, y+.17, x-.17, y-.17)
            ret += "\\draw[-,thick] ({} cm, {} cm) -- ({} cm, {} cm);\n".format(
                                    x+.17, y-.17, x-.17, y+.17)
        return ret

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A', [3,2]])._latex_dynkin_diagram())
            \draw (0 cm, 0 cm) -- (10 cm, 0 cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$-3$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$-2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$-1$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            \draw[-,thick] (6.17 cm, 0.17 cm) -- (5.83 cm, -0.17 cm);
            \draw[-,thick] (6.17 cm, -0.17 cm) -- (5.83 cm, 0.17 cm);
            \draw[fill=white] (8 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (10 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};

            sage: print(CartanType(['A', [0,2]])._latex_dynkin_diagram())
            \draw (0 cm, 0 cm) -- (4 cm, 0 cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            \draw[-,thick] (0.17 cm, 0.17 cm) -- (-0.17 cm, -0.17 cm);
            \draw[-,thick] (0.17 cm, -0.17 cm) -- (-0.17 cm, 0.17 cm);
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};

            sage: print(CartanType(['A', [2,0]])._latex_dynkin_diagram())
            \draw (0 cm, 0 cm) -- (4 cm, 0 cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$-2$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$-1$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            \draw[-,thick] (4.17 cm, 0.17 cm) -- (3.83 cm, -0.17 cm);
            \draw[-,thick] (4.17 cm, -0.17 cm) -- (3.83 cm, 0.17 cm);

            sage: print(CartanType(['A', [0,0]])._latex_dynkin_diagram())
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            \draw[-,thick] (0.17 cm, 0.17 cm) -- (-0.17 cm, -0.17 cm);
            \draw[-,thick] (0.17 cm, -0.17 cm) -- (-0.17 cm, 0.17 cm);
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        if self.n + self.m > 1:
            ret = "\\draw (0 cm, 0 cm) -- ({} cm, 0 cm);\n".format((self.n+self.m)*node_dist)
        else:
            ret = ""
        return ret + "".join(node((self.m+i)*node_dist, 0, label(i))
                             for i in self.index_set())

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of the Dynkin diagram.

        EXAMPLES::

            sage: t = CartanType(['A', [3,2]])
            sage: print(t.ascii_art())
            O---O---O---X---O---O
            -3  -2  -1  0   1   2
            sage: t = CartanType(['A', [3,7]])
            sage: print(t.ascii_art())
            O---O---O---X---O---O---O---O---O---O---O
            -3  -2  -1  0   1   2   3   4   5   6   7

            sage: t = CartanType(['A', [0,7]])
            sage: print(t.ascii_art())
            X---O---O---O---O---O---O---O
            0   1   2   3   4   5   6   7
            sage: t = CartanType(['A', [0,0]])
            sage: print(t.ascii_art())
            X
            0
            sage: t = CartanType(['A', [5,0]])
            sage: print(t.ascii_art())
            O---O---O---O---O---X
            -5  -4  -3  -2  -1  0
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = lambda i: 'O'
        ret = "---".join(node(label(i)) for i in range(1,self.m+1))
        if self.m == 0:
            if self.n == 0:
                ret = "X"
            else:
                ret += "X---"
        else:
            if self.n == 0:
                ret += "---X"
            else:
                ret += "---X---"
        ret += "---".join(node(label(i)) for i in range(1,self.n+1)) + "\n"
        ret += "".join("{!s:4}".format(label(-i)) for i in reversed(range(1,self.m+1)))
        ret += "{!s:4}".format(label(0))
        ret += "".join("{!s:4}".format(label(i)) for i in range(1,self.n+1))
        return ret
