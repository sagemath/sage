"""
Root system data for type D
"""
# ****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from . import ambient_space


class AmbientSpace(ambient_space.AmbientSpace):
    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['D',3]).ambient_space()
            sage: e.dimension()
            3
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j, p1, p2):
        """
        Note that indexing starts at 0.

        EXAMPLES::

            sage: e = RootSystem(['D',3]).ambient_space()
            sage: e.root(0, 1, 1, 1)
            (-1, -1, 0)
            sage: e.root(0, 0, 1, 1)
            (-1, 0, 0)
        """
        if i != j:
            return (-1)**p1 * self.monomial(i) + (-1)**p2 * self.monomial(j)
        return (-1)**p1 * self.monomial(i)

    def simple_root(self, i):
        """
        EXAMPLES::

            sage: RootSystem(['D',4]).ambient_space().simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1), 4: (0, 0, 1, 1)}
        """
        if i not in self.index_set():
            raise ValueError("{} is not in the index set".format(i))
        return self.root(i-1, i, 0, 1) if i < self.n else self.root(self.n-2, self.n-1, 0, 0)

    def positive_roots(self):
        """
        EXAMPLES::

            sage: RootSystem(['D',4]).ambient_space().positive_roots()
            [(1, 1, 0, 0),
             (1, 0, 1, 0),
             (0, 1, 1, 0),
             (1, 0, 0, 1),
             (0, 1, 0, 1),
             (0, 0, 1, 1),
             (1, -1, 0, 0),
             (1, 0, -1, 0),
             (0, 1, -1, 0),
             (1, 0, 0, -1),
             (0, 1, 0, -1),
             (0, 0, 1, -1)]
        """
        res = []
        for p in [0, 1]:
            for j in range(self.n):
                res.extend([self.root(i, j, 0, p) for i in range(j)])
        return res

    def negative_roots(self):
        """
        EXAMPLES::

            sage: RootSystem(['D',4]).ambient_space().negative_roots()
            [(-1, 1, 0, 0),
             (-1, 0, 1, 0),
             (0, -1, 1, 0),
             (-1, 0, 0, 1),
             (0, -1, 0, 1),
             (0, 0, -1, 1),
             (-1, -1, 0, 0),
             (-1, 0, -1, 0),
             (0, -1, -1, 0),
             (-1, 0, 0, -1),
             (0, -1, 0, -1),
             (0, 0, -1, -1)]
        """
        res = []
        for p in [0, 1]:
            for j in range(self.n):
                res.extend([self.root(i, j, 1, p) for i in range(j)])
        return res

    def fundamental_weight(self, i):
        """
        EXAMPLES::

            sage: RootSystem(['D',4]).ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1/2, 1/2, 1/2, -1/2), 4: (1/2, 1/2, 1/2, 1/2)}
        """
        if i not in self.index_set():
            raise ValueError("{} is not in the index set".format(i))
        n = self.dimension()
        if i == n:
            return self.sum(self.monomial(j) for j in range(n)) / 2
        elif i == n - 1:
            return (self.sum(self.monomial(j) for j in range(n-1)) - self.monomial(n-1)) / 2
        else:
            return self.sum(self.monomial(j) for j in range(i))


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_A', 'ambient_space',  AmbientSpace)

from sage.misc.cachefunc import cached_method
from .cartan_type import CartanType_standard_finite, CartanType_simply_laced, CartanType_simple


class CartanType(CartanType_standard_finite, CartanType_simply_laced):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['D',4])
            sage: ct
            ['D', 4]
            sage: ct._repr_(compact = True)
            'D4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.dual()
            ['D', 4]
            sage: ct.affine()
            ['D', 4, 1]

            sage: ct = CartanType(['D',2])
            sage: ct.is_irreducible()
            False
            sage: ct.dual()
            ['D', 2]
            sage: ct.affine()
            Traceback (most recent call last):
            ...
            ValueError: ['D', 2, 1] is not a valid Cartan type


        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 2
        CartanType_standard_finite.__init__(self, "D", n)
        if n >= 3:
            self._add_abstract_superclass(CartanType_simple)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['D',4]))
            D_{4}
        """
        return "D_{%s}" % self.n

    AmbientSpace = AmbientSpace

    def is_atomic(self):
        """
        Implement :meth:`CartanType_abstract.is_atomic`.

        `D_2` is atomic, like all `D_n`, despite being non irreducible.

        EXAMPLES::

            sage: CartanType(["D",2]).is_atomic()
            True
            sage: CartanType(["D",2]).is_irreducible()
            False
        """
        return True

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['D',4]).coxeter_number()
            6
        """
        return 2*self.n - 2

    def dual_coxeter_number(self):
        """
        Return the dual Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['D',4]).dual_coxeter_number()
            6
        """
        return 2*self.n - 2

    @cached_method
    def dynkin_diagram(self):
        """
        Return a Dynkin diagram for type D.

        EXAMPLES::

            sage: d = CartanType(['D',5]).dynkin_diagram(); d                           # needs sage.graphs
                    O 5
                    |
                    |
            O---O---O---O
            1   2   3   4
            D5
            sage: d.edges(sort=True)                                                    # needs sage.graphs
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 1),
             (3, 4, 1), (3, 5, 1), (4, 3, 1), (5, 3, 1)]

            sage: d = CartanType(['D',4]).dynkin_diagram(); d                           # needs sage.graphs
                O 4
                |
                |
            O---O---O
            1   2   3
            D4
            sage: d.edges(sort=True)                                                    # needs sage.graphs
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (2, 4, 1), (3, 2, 1), (4, 2, 1)]

            sage: d = CartanType(['D',3]).dynkin_diagram(); d                           # needs sage.graphs
            O 3
            |
            |
            O---O
            1   2
            D3
            sage: d.edges(sort=True)                                                    # needs sage.graphs
            [(1, 2, 1), (1, 3, 1), (2, 1, 1), (3, 1, 1)]


            sage: d = CartanType(['D',2]).dynkin_diagram(); d                           # needs sage.graphs
            O   O
            1   2
            D2
            sage: d.edges(sort=True)                                                    # needs sage.graphs
            []
        """
        from .dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        n = self.n
        if n >= 3:
            for i in range(1, n-1):
                g.add_edge(i, i+1)
            g.add_edge(n-2, n)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['D',4])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (2 cm,0);
            \draw (2 cm,0) -- (4 cm,0.7 cm);
            \draw (2 cm,0) -- (4 cm,-0.7 cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0.7 cm) circle (.25cm) node[right=3pt]{$4$};
            \draw[fill=white] (4 cm, -0.7 cm) circle (.25cm) node[right=3pt]{$3$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        if self.n == 2:
            ret = node(0, 0, label(1))
            ret += node(node_dist, 0, label(2))
            return ret
        rt_most = (self.n-2) * node_dist
        center_point = rt_most - node_dist
        ret = "\\draw (0 cm,0) -- (%s cm,0);\n" % center_point
        ret += "\\draw (%s cm,0) -- (%s cm,0.7 cm);\n" % (center_point, rt_most)
        ret += "\\draw (%s cm,0) -- (%s cm,-0.7 cm);\n" % (center_point, rt_most)
        for i in range(self.n-2):
            ret += node(i*node_dist, 0, label(i+1))
        ret += node(rt_most, 0.7, label(self.n), 'right=3pt')
        ret += node(rt_most, -0.7, label(self.n-1), 'right=3pt')
        return ret

    def ascii_art(self, label=None, node=None):
        """
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['D',3]).ascii_art())
            O 3
            |
            |
            O---O
            1   2
            sage: print(CartanType(['D',4]).ascii_art())
                O 4
                |
                |
            O---O---O
            1   2   3
            sage: print(CartanType(['D',4]).ascii_art(label = lambda x: x+2))
                O 6
                |
                |
            O---O---O
            3   4   5
            sage: print(CartanType(['D',6]).ascii_art(label = lambda x: x+2))
                        O 8
                        |
                        |
            O---O---O---O---O
            3   4   5   6   7
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        n = self.n
        if n == 2:
            ret = "{}   {}\n".format(node(label(1)), node(label(2)))
            return ret + "{!s:4}{!s:4}".format(label(1), label(2))
        ret = (4*(n-3))*" "+"{} {}\n".format(node(label(n)), label(n))
        ret += ((4*(n-3))*" " + "|\n")*2
        ret += "---".join(node(label(i)) for i in range(1, n)) + "\n"
        ret += "".join("{!s:4}".format(label(i)) for i in range(1,n))
        return ret


# For unpickling backward compatibility (Sage <= 4.1)
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_D', 'ambient_space',  AmbientSpace)
