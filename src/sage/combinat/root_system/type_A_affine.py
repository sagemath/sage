"""
Root system data for (untwisted) type A affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .cartan_type import CartanType_standard_untwisted_affine, CartanType_simply_laced


class CartanType(CartanType_standard_untwisted_affine):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4,1])
            sage: ct
            ['A', 4, 1]
            sage: ct._repr_(compact = True)
            'A4~'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            False
            sage: ct.is_affine()
            True
            sage: ct.is_untwisted_affine()
            True
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.classical()
            ['A', 4]
            sage: ct.dual()
            ['A', 4, 1]

            sage: ct = CartanType(['A', 1, 1])
            sage: ct.is_simply_laced()
            False
            sage: ct.dual()
            ['A', 1, 1]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 1
        CartanType_standard_untwisted_affine.__init__(self, "A", n)
        if n >= 2:
            self._add_abstract_superclass(CartanType_simply_laced)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: ct = CartanType(['A',4,1])
            sage: latex(ct)
            A_{4}^{(1)}
        """
        return "A_{%s}^{(1)}" % self.n

    def dynkin_diagram(self):
        """
        Return the extended Dynkin diagram for affine type A.

        EXAMPLES::

            sage: a = CartanType(['A',3,1]).dynkin_diagram(); a                         # needs sage.graphs
             0
             O-------+
             |       |
             |       |
             O---O---O
             1   2   3
             A3~
            sage: a.edges(sort=True)                                                    # needs sage.graphs
            [(0, 1, 1),
             (0, 3, 1),
             (1, 0, 1),
             (1, 2, 1),
             (2, 1, 1),
             (2, 3, 1),
             (3, 0, 1),
             (3, 2, 1)]

            sage: a = DynkinDiagram(['A',1,1]); a                                       # needs sage.graphs
            O<=>O
            0   1
            A1~
            sage: a.edges(sort=True)                                                    # needs sage.graphs
            [(0, 1, 2), (1, 0, 2)]
        """
        from .dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)

        if n == 1:
            g.add_edge(0, 1, 2)
            g.add_edge(1, 0, 2)
        else:
            for i in range(1, n):
                g.add_edge(i, i+1)
            g.add_edge(0, 1)
            g.add_edge(0, n)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A',4,1])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (6 cm,0);
            \draw (0 cm,0) -- (3.0 cm, 1.2 cm);
            \draw (3.0 cm, 1.2 cm) -- (6 cm, 0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            \draw[fill=white] (3.0 cm, 1.2 cm) circle (.25cm) node[anchor=south east]{$0$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        if self.n == 1:
            ret = "\\draw (0, 0.1 cm) -- +(%s cm,0);\n" % node_dist
            ret += "\\draw (0, -0.1 cm) -- +(%s cm,0);\n" % node_dist
            ret += self._latex_draw_arrow_tip(0.33*node_dist-0.2, 0, 180)
            ret += self._latex_draw_arrow_tip(0.66*node_dist+0.2, 0, 0)
            ret += node(0, 0, label(0))
            ret += node(node_dist, 0, label(1))
            return ret
        rt_most = (self.n-1)*node_dist
        mid = 0.5 * rt_most
        ret = "\\draw (0 cm,0) -- (%s cm,0);\n" % rt_most
        ret += "\\draw (0 cm,0) -- (%s cm, 1.2 cm);\n" % mid
        ret += "\\draw (%s cm, 1.2 cm) -- (%s cm, 0);\n" % (mid, rt_most)
        for i in range(self.n):
            ret += node(i*node_dist, 0, label(i+1))
        ret += node(mid, 1.2, label(0), 'anchor=south east')
        return ret

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A',3,1]).ascii_art())
            0
            O-------+
            |       |
            |       |
            O---O---O
            1   2   3

            sage: print(CartanType(['A',5,1]).ascii_art(label = lambda x: x+2))
            2
            O---------------+
            |               |
            |               |
            O---O---O---O---O
            3   4   5   6   7

            sage: print(CartanType(['A',1,1]).ascii_art())
            O<=>O
            0   1

            sage: print(CartanType(['A',1,1]).ascii_art(label = lambda x: x+2))
            O<=>O
            2   3
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        n = self.n
        if n == 1:
            l0 = label(0)
            l1 = label(1)
            return "{}<=>{}\n{!s:4}{}".format(node(l0), node(l1), l0, l1)
        ret = "{}\n{}".format(label(0), node(label(0)))
        ret += "----"*(n-2) + "---+\n|" + "    "*(n-2) + "   |\n|" + "    "*(n-2) + "   |\n"
        ret += "---".join(node(label(i)) for i in range(1,n+1)) + "\n"
        ret += "".join("{!s:4}".format(label(i)) for i in range(1,n+1))
        return ret

    def dual(self):
        """
        Type `A_1^1` is self dual despite not being simply laced.

        EXAMPLES::

            sage: CartanType(['A',1,1]).dual()
            ['A', 1, 1]
        """
        return self

    def _default_folded_cartan_type(self):
        r"""
        Return the default folded Cartan type.

        In general, this just returns ``self`` in ``self`` with `\sigma` as
        the identity map.

        EXAMPLES::

            sage: CartanType(['A',1,1])._default_folded_cartan_type()
            ['A', 1, 1] as a folding of ['A', 3, 1]
            sage: CartanType(['A',3,1])._default_folded_cartan_type()
            ['A', 3, 1] as a folding of ['A', 3, 1]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        if self.n == 1:
            return CartanTypeFolded(self, ['A', 3, 1], [[0,2], [1,3]])
        return CartanTypeFolded(self, self, [[i] for i in self.index_set()])
