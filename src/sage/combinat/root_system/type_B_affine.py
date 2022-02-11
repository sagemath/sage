"""
Root system data for (untwisted) type B affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .cartan_type import CartanType_standard_untwisted_affine
class CartanType(CartanType_standard_untwisted_affine):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['B',4,1])
            sage: ct
            ['B', 4, 1]
            sage: ct._repr_(compact = True)
            'B4~'

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
            False
            sage: ct.classical()
            ['B', 4]
            sage: ct.dual()
            ['B', 4, 1]^*
            sage: ct.dual().is_untwisted_affine()
            False

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 1
        CartanType_standard_untwisted_affine.__init__(self, "B", n)

    def dynkin_diagram(self):
        """
        Return the extended Dynkin diagram for affine type `B`.

        EXAMPLES::

            sage: b = CartanType(['B',3,1]).dynkin_diagram()
            sage: b
                O 0
                |
                |
            O---O=>=O
            1   2   3
            B3~
            sage: sorted(b.edges())
            [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]

            sage: b = CartanType(['B',2,1]).dynkin_diagram(); b
            O=>=O=<=O
            0   2   1
            B2~
            sage: sorted(b.edges())
            [(0, 2, 2), (1, 2, 2), (2, 0, 1), (2, 1, 1)]

            sage: b = CartanType(['B',1,1]).dynkin_diagram(); b
            O<=>O
            0   1
            B1~
            sage: sorted(b.edges())
            [(0, 1, 2), (1, 0, 2)]

        """
        from . import cartan_type
        n = self.n
        if n == 1:
            res = cartan_type.CartanType(["A",1,1]).dynkin_diagram()
            res._cartan_type = self
            return res
        if n == 2:
            res = cartan_type.CartanType(["C",2,1]).relabel({0:0, 1:2, 2:1}).dynkin_diagram()
            res._cartan_type = self
            return res
        from .dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        for i in range(1, n):
            g.add_edge(i, i+1)
        g.set_edge_label(n-1, n, 2)
        g.add_edge(0,2)
        return g

    def _latex_dynkin_diagram(self, label=lambda i: i, node=None, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['B',4,1])._latex_dynkin_diagram())
            \draw (0,0.7 cm) -- (2 cm,0);
            \draw (0,-0.7 cm) -- (2 cm,0);
            \draw (2 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(5.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0.7 cm) circle (.25cm) node[left=3pt]{$0$};
            \draw[fill=white] (0 cm, -0.7 cm) circle (.25cm) node[left=3pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>

            sage: print(CartanType(['B',4,1]).dual()._latex_dynkin_diagram())
            \draw (0,0.7 cm) -- (2 cm,0);
            \draw (0,-0.7 cm) -- (2 cm,0);
            \draw (2 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(4.8, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0.7 cm) circle (.25cm) node[left=3pt]{$0$};
            \draw[fill=white] (0 cm, -0.7 cm) circle (.25cm) node[left=3pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>
        """
        if node is None:
            node = self._latex_draw_node
        if self.n == 1:
            from . import cartan_type
            return cartan_type.CartanType(["A",1,1])._latex_dynkin_diagram(label, node, node_dist)
        elif self.n == 2:
            from . import cartan_type
            return cartan_type.CartanType(["C",2,1])._latex_dynkin_diagram(label, node, node_dist, dual)
        n = self.n
        single_end = (n-2)*node_dist # Where the single line ends
        ret = "\\draw (0,0.7 cm) -- (%s cm,0);\n"%node_dist
        ret += "\\draw (0,-0.7 cm) -- (%s cm,0);\n"%node_dist
        ret += "\\draw (%s cm,0) -- (%s cm,0);\n"%(node_dist, single_end)
        ret += "\\draw (%s cm, 0.1 cm) -- +(%s cm,0);\n"%(single_end, node_dist)
        ret += "\\draw (%s cm, -0.1 cm) -- +(%s cm,0);\n"%(single_end, node_dist)
        if dual:
            ret += self._latex_draw_arrow_tip(single_end+0.5*node_dist-0.2, 0, 180)
        else:
            ret += self._latex_draw_arrow_tip(single_end+0.5*node_dist+0.2, 0, 0)
        ret += node(0, 0.7, label(0), 'left=3pt')
        ret +=node(0, -0.7, label(1), 'left=3pt')
        for i in range(1, n):
            ret += node(i*node_dist, 0, label(i+1))
        return ret

    def ascii_art(self, label=lambda i: i, node=None):
        """
        Return an ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['B',3,1]).ascii_art())
                O 0
                |
                |
            O---O=>=O
            1   2   3

            sage: print(CartanType(['B',5,1]).ascii_art(label = lambda x: x+2))
                O 2
                |
                |
            O---O---O---O=>=O
            3   4   5   6   7

            sage: print(CartanType(['B',2,1]).ascii_art(label = lambda x: x+2))
            O=>=O=<=O
            2   4   3
            sage: print(CartanType(['B',1,1]).ascii_art(label = lambda x: x+2))
            O<=>O
            2   3
        """
        n = self.n
        from .cartan_type import CartanType
        if node is None:
            node = self._ascii_art_node
        if n == 1:
            return CartanType(["A",1,1]).ascii_art(label, node)
        if n == 2:
            return CartanType(["C",2,1]).relabel({0:0, 1:2, 2:1}).ascii_art(label, node)
        ret  = "    {} {}\n    |\n    |\n".format(node(label(0)), label(0))
        ret += "---".join(node(label(i)) for i in range(1,n)) + "=>={}\n".format(node(label(n)))
        ret += "".join("{!s:4}".format(label(i)) for i in range(1,n+1))
        return ret

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['B', 4, 1])._default_folded_cartan_type()
            ['B', 4, 1] as a folding of ['D', 5, 1]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        n = self.n
        if n == 1:
            return CartanTypeFolded(self, ['A', 1, 1], [[0], [1]])
        return CartanTypeFolded(self, ['D', n + 1, 1],
            [[i] for i in range(n)] + [[n, n+1]])

