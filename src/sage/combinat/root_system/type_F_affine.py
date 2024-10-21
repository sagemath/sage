# sage_setup: distribution = sagemath-modules
"""
Root system data for (untwisted) type F affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .cartan_type import CartanType_standard_untwisted_affine


class CartanType(CartanType_standard_untwisted_affine):
    def __init__(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['F',4,1])
            sage: ct
            ['F', 4, 1]
            sage: ct._repr_(compact = True)
            'F4~'

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
            ['F', 4]
            sage: ct.dual()
            ['F', 4, 1]^*
            sage: ct.dual().is_untwisted_affine()
            False

        TESTS::

            sage: TestSuite(ct).run()
        """
        CartanType_standard_untwisted_affine.__init__(self, "F", 4)

    def dynkin_diagram(self):
        """
        Return the extended Dynkin diagram for affine type F.

        EXAMPLES::

            sage: f = CartanType(['F', 4, 1]).dynkin_diagram(); f                       # needs sage.graphs
            O---O---O=>=O---O
            0   1   2   3   4
            F4~
            sage: f.edges(sort=True)                                                    # needs sage.graphs
            [(0, 1, 1), (1, 0, 1), (1, 2, 1), (2, 1, 1),
             (2, 3, 2), (3, 2, 1), (3, 4, 1), (4, 3, 1)]
        """
        from .dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        for i in range(1, 4):
            g.add_edge(i, i+1)
        g.set_edge_label(2,3,2)
        g.add_edge(0, 1)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['F',4,1])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (2 cm,0);
            {
            \pgftransformxshift{2 cm}
            \draw (0 cm,0) -- (2 cm,0);
            \draw (2 cm, 0.1 cm) -- +(2 cm,0);
            \draw (2 cm, -0.1 cm) -- +(2 cm,0);
            \draw (4.0 cm,0) -- +(2 cm,0);
            \draw[shift={(3.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            }
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        ret = "\\draw (0 cm,0) -- (%s cm,0);\n" % node_dist
        ret += "{\n\\pgftransformxshift{%s cm}\n" % node_dist
        ret += self.classical()._latex_dynkin_diagram(label, node, node_dist, dual)
        ret += "}\n" + node(0, 0, label(0))
        return ret

    def ascii_art(self, label=None, node=None):
        """
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['F',4,1]).ascii_art(label = lambda x: x+2))
            O---O---O=>=O---O
            2   3   4   5   6
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        ret = "{}---{}---{}=>={}---{}\n".format(node(label(0)), node(label(1)),
                             node(label(2)), node(label(3)), node(label(4)))
        ret += ("{!s:4}"*5 + "\n").format(label(0), label(1), label(2), label(3), label(4))
        return ret

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['F', 4, 1])._default_folded_cartan_type()
            ['F', 4, 1] as a folding of ['E', 6, 1]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        return CartanTypeFolded(self, ['E', 6, 1], [[0], [2], [4], [3, 5], [1, 6]])
