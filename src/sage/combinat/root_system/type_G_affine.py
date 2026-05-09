"""
Root system data for (untwisted) type G affine
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

            sage: ct = CartanType(['G',2,1])
            sage: ct
            ['G', 2, 1]
            sage: ct._repr_(compact = True)
            'G2~'

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
            ['G', 2]
            sage: ct.dual()
            ['G', 2, 1]^*
            sage: ct.dual().is_untwisted_affine()
            False

        TESTS::

            sage: TestSuite(ct).run()
        """
        CartanType_standard_untwisted_affine.__init__(self, "G",2)

    def dynkin_diagram(self):
        """
        Return the extended Dynkin diagram for type G.

        EXAMPLES::

            sage: g = CartanType(['G',2,1]).dynkin_diagram(); g                         # needs sage.graphs
              3
            O=<=O---O
            1   2   0
            G2~
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 3)]
        """
        from .dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        g.add_edge(1, 2)
        g.set_edge_label(2,1,3)
        g.add_edge(0, 2)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['G',2,1])._latex_dynkin_diagram())
            \draw (2 cm,0) -- (4.0 cm,0);
            \draw (0, 0.15 cm) -- +(2 cm,0);
            \draw (0, -0.15 cm) -- +(2 cm,0);
            \draw (0,0) -- (2 cm,0);
            \draw (0, 0.15 cm) -- +(2 cm,0);
            \draw (0, -0.15 cm) -- +(2 cm,0);
            \draw[shift={(0.8, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda x: x
        if node is None:
            node = self._latex_draw_node
        ret = "\\draw (%s cm,0) -- (%s cm,0);\n" % (node_dist, node_dist*2.0)
        ret += "\\draw (0, 0.15 cm) -- +(%s cm,0);\n" % node_dist
        ret += "\\draw (0, -0.15 cm) -- +(%s cm,0);\n" % node_dist
        ret += self.classical()._latex_dynkin_diagram(label, node, node_dist, dual)
        ret += node(2*node_dist, 0, label(0))
        return ret

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['G',2,1]).ascii_art(label = lambda x: x+2))
              3
            O=<=O---O
            3   4   2
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        ret = "  3\n{}=<={}---{}".format(node(label(1)), node(label(2)), node(label(0)))
        return ret + "\n{!s:4}{!s:4}{!s:4}".format(label(1), label(2), label(0))

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['G', 2, 1])._default_folded_cartan_type()
            ['G', 2, 1] as a folding of ['D', 4, 1]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        return CartanTypeFolded(self, ['D', 4, 1], [[0], [1, 3, 4], [2]])
