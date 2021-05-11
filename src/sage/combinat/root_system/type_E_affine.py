"""
Root system data for (untwisted) type E affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .cartan_type import CartanType_standard_untwisted_affine, CartanType_simply_laced
class CartanType(CartanType_standard_untwisted_affine, CartanType_simply_laced):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['E',6,1])
            sage: ct
            ['E', 6, 1]
            sage: ct._repr_(compact = True)
            'E6~'

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
            ['E', 6]
            sage: ct.dual()
            ['E', 6, 1]

        TESTS::

            sage: TestSuite(ct).run()
        """
        if n < 6 or n > 8:
            raise ValueError("Invalid Cartan Type for Type E")
        CartanType_standard_untwisted_affine.__init__(self, "E", n)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['E',7,1]))
            E_7^{(1)}
        """
        return "E_%s^{(1)}"%self.n

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type E.

        EXAMPLES::

            sage: e = CartanType(['E', 6, 1]).dynkin_diagram()
            sage: e
                    O 0
                    |
                    |
                    O 2
                    |
                    |
            O---O---O---O---O
            1   3   4   5   6
            E6~
            sage: sorted(e.edges())
            [(0, 2, 1),
             (1, 3, 1),
             (2, 0, 1),
             (2, 4, 1),
             (3, 1, 1),
             (3, 4, 1),
             (4, 2, 1),
             (4, 3, 1),
             (4, 5, 1),
             (5, 4, 1),
             (5, 6, 1),
             (6, 5, 1)]

            sage: e = CartanType(['E', 7, 1]).dynkin_diagram()
            sage: e
                        O 2
                        |
                        |
            O---O---O---O---O---O---O
            0   1   3   4   5   6   7
            E7~
            sage: sorted(e.edges())
            [(0, 1, 1), (1, 0, 1), (1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1),
             (4, 2, 1), (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1),
             (6, 5, 1), (6, 7, 1), (7, 6, 1)]
            sage: e = CartanType(['E', 8, 1]).dynkin_diagram()
            sage: e
                    O 2
                    |
                    |
            O---O---O---O---O---O---O---O
            1   3   4   5   6   7   8   0
            E8~
            sage: sorted(e.edges())
            [(0, 8, 1), (1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1),
             (4, 2, 1), (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1),
             (6, 5, 1), (6, 7, 1), (7, 6, 1), (7, 8, 1), (8, 0, 1), (8, 7, 1)]

        """
        from .dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        g.add_edge(1,3)
        g.add_edge(2,4)
        for i in range(3,n):
            g.add_edge(i, i+1)
        if n == 6:
            g.add_edge(0, 2)
        elif n == 7:
            g.add_edge(0, 1)
        elif n == 8:
            g.add_edge(0, 8)
        else:
            raise ValueError("Invalid Cartan Type for Type E affine")
        return g

    def _latex_dynkin_diagram(self, label=lambda i: i, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['E',7,1])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (12 cm,0);
            \draw (6 cm, 0 cm) -- +(0,2 cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            \draw[fill=white] (8 cm, 0 cm) circle (.25cm) node[below=4pt]{$5$};
            \draw[fill=white] (10 cm, 0 cm) circle (.25cm) node[below=4pt]{$6$};
            \draw[fill=white] (12 cm, 0 cm) circle (.25cm) node[below=4pt]{$7$};
            \draw[fill=white] (6 cm, 2 cm) circle (.25cm) node[right=3pt]{$2$};
            <BLANKLINE>
        """
        n = self.n
        if node is None:
            node = self._latex_draw_node

        if n == 7:
            ret = "\\draw (0 cm,0) -- (%s cm,0);\n"%((n-1)*node_dist)
            ret += "\\draw (%s cm, 0 cm) -- +(0,%s cm);\n"%(3*node_dist, node_dist)
            ret += node(0, 0, label(0))
            ret += node(node_dist, 0, label(1))
            for i in range(2, n):
                ret += node(i*node_dist, 0, label(i+1))
            ret += node(3*node_dist, node_dist, label(2), "right=3pt")
            return ret

        ret = "\\draw (0 cm,0) -- (%s cm,0);\n"%((n-2)*node_dist)
        ret += "\\draw (%s cm, 0 cm) -- +(0,%s cm);\n"%(2*node_dist, node_dist)

        if n == 6:
            ret += "\\draw (%s cm, %s cm) -- +(0,%s cm);\n"%(2*node_dist, node_dist, node_dist)
            ret += node(2*node_dist, 2*node_dist, label(0), "right=3pt")
        else: # n == 8
            ret += "\\draw (%s cm,0) -- +(%s cm,0);\n"%((n-2)*node_dist, node_dist)
            ret += node((n-1)*node_dist, 0, label(0))

        ret += node(0, 0, label(1))
        for i in range(1, n-1):
            ret += node(i*node_dist, 0, label(i+2))
        ret += node(2*node_dist, node_dist, label(2), "right=3pt")
        return ret

    def ascii_art(self, label=lambda x: x, node=None):
        """
        Return an ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['E',6,1]).ascii_art(label = lambda x: x+2))
                    O 2
                    |
                    |
                    O 4
                    |
                    |
            O---O---O---O---O
            3   5   6   7   8
            sage: print(CartanType(['E',7,1]).ascii_art(label = lambda x: x+2))
                        O 4
                        |
                        |
            O---O---O---O---O---O---O
            2   3   5   6   7   8   9
            sage: print(CartanType(['E',8,1]).ascii_art(label = lambda x: x-3))
                    O -1
                    |
                    |
            O---O---O---O---O---O---O---O
            -2  0   1   2   3   4   5   -3
        """
        n = self.n
        if node is None:
            node = self._ascii_art_node
        if n == 6:
            ret = "        {} {}\n        |\n        |\n".format(node(label(0)), label(0))
            return ret + self.classical().ascii_art(label, node)
        elif n == 7:
            ret = "            {} {}\n            |\n            |\n".format(node(label(2)), label(2))
            labels = [label(_) for _ in [0,1,3,4,5,6,7]]
            nodes = [node(_) for _ in labels]
            return ret + '---'.join(n for n in nodes) + '\n' + "".join("{!s:4}".format(i) for i in labels)
        elif n == 8:
            ret = "        {} {}\n        |\n        |\n".format(node(label(2)), label(2))
            labels = [label(_) for _ in [1,3,4,5,6,7,8,0]]
            nodes = [node(_) for _ in labels]
            return ret + '---'.join(n for n in nodes) + '\n' + "".join("{!s:4}".format(i) for i in labels)

