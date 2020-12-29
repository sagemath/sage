"""
Root system data for (untwisted) type D affine
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

            sage: ct = CartanType(['D',4,1])
            sage: ct
            ['D', 4, 1]
            sage: ct._repr_(compact = True)
            'D4~'

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
            ['D', 4]
            sage: ct.dual()
            ['D', 4, 1]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 3
        CartanType_standard_untwisted_affine.__init__(self, "D", n)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type D.

        EXAMPLES::

           sage: d = CartanType(['D', 6, 1]).dynkin_diagram()
           sage: d
              0 O       O 6
                |       |
                |       |
            O---O---O---O---O
            1   2   3   4   5
            D6~
           sage: sorted(d.edges())
           [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 1), (2, 3, 1),
            (3, 2, 1), (3, 4, 1), (4, 3, 1), (4, 5, 1), (4, 6, 1), (5, 4, 1), (6, 4, 1)]

           sage: d = CartanType(['D', 4, 1]).dynkin_diagram()
           sage: d
               O 4
               |
               |
           O---O---O
           1   |2  3
               |
               O 0
           D4~
           sage: sorted(d.edges())
           [(0, 2, 1),
            (1, 2, 1),
            (2, 0, 1),
            (2, 1, 1),
            (2, 3, 1),
            (2, 4, 1),
            (3, 2, 1),
            (4, 2, 1)]

           sage: d = CartanType(['D', 3, 1]).dynkin_diagram()
           sage: d
           0
           O-------+
           |       |
           |       |
           O---O---O
           3   1   2
           D3~
           sage: sorted(d.edges())
           [(0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1), (2, 0, 1), (2, 1, 1), (3, 0, 1), (3, 1, 1)]

        """
        from .dynkin_diagram import DynkinDiagram_class
        n = self.n
        if n == 3:
            from . import cartan_type
            res = cartan_type.CartanType(["A",3,1]).relabel({0:0, 1:3, 2:1, 3: 2}).dynkin_diagram()
            res._cartan_type = self
            return res
        g = DynkinDiagram_class(self)
        for i in range(1, n-1):
            g.add_edge(i, i+1)
        g.add_edge(n-2,n)
        g.add_edge(0,2)
        return g

    def _latex_dynkin_diagram(self, label=lambda i: i, node=None, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['D',4,1])._latex_dynkin_diagram())
            \draw (0,0.7 cm) -- (2 cm,0);
            \draw (0,-0.7 cm) -- (2 cm,0);
            \draw (2 cm,0) -- (2 cm,0);
            \draw (2 cm,0) -- (4 cm,0.7 cm);
            \draw (2 cm,0) -- (4 cm,-0.7 cm);
            \draw[fill=white] (0 cm, 0.7 cm) circle (.25cm) node[left=3pt]{$0$};
            \draw[fill=white] (0 cm, -0.7 cm) circle (.25cm) node[left=3pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0.7 cm) circle (.25cm) node[right=3pt]{$4$};
            \draw[fill=white] (4 cm, -0.7 cm) circle (.25cm) node[right=3pt]{$3$};
            <BLANKLINE>
        """
        if node is None:
            node = self._latex_draw_node
        n = self.n
        if n == 3:
            from . import cartan_type
            relabel = {0:label(0), 1:label(3), 2:label(1), 3:label(2)}
            return cartan_type.CartanType(["A",3,1]).relabel(relabel)._latex_dynkin_diagram(node_dist=node_dist)
        rt_most = (n - 2) * node_dist
        center_point = rt_most - node_dist
        ret = "\\draw (0,0.7 cm) -- (%s cm,0);\n"%node_dist
        ret += "\\draw (0,-0.7 cm) -- (%s cm,0);\n"%node_dist
        ret += "\\draw (%s cm,0) -- (%s cm,0);\n"%(node_dist, center_point)
        ret += "\\draw (%s cm,0) -- (%s cm,0.7 cm);\n"%(center_point, rt_most)
        ret += "\\draw (%s cm,0) -- (%s cm,-0.7 cm);\n"%(center_point, rt_most)
        ret += node(0, 0.7, label(0), "left=3pt")
        ret += node(0, -0.7, label(1), "left=3pt")
        for i in range(1, self.n-2):
            ret += node(i*node_dist, 0, label(i+1))
        ret += node(rt_most, 0.7, label(n), "right=3pt")
        ret += node(rt_most, -0.7, label(n-1), "right=3pt")
        return ret

    def ascii_art(self, label=lambda i: i, node=None):
        """
        Return an ascii art representation of the extended Dynkin diagram.

        TESTS::

            sage: print(CartanType(['D',6,1]).ascii_art(label = lambda x: x+2))
              2 O       O 8
                |       |
                |       |
            O---O---O---O---O
            3   4   5   6   7

            sage: print(CartanType(['D',4,1]).ascii_art(label = lambda x: x+2))
                O 6
                |
                |
            O---O---O
            3   |4  5
                |
                O 2

            sage: print(CartanType(['D',3,1]).ascii_art(label = lambda x: x+2))
            2
            O-------+
            |       |
            |       |
            O---O---O
            5   3   4
        """
        if node is None:
            node = self._ascii_art_node
        n = self.n
        if n == 3:
            from . import cartan_type
            return cartan_type.CartanType(["A",3,1]).relabel({0:0, 1:3, 2:1, 3: 2}).ascii_art(label, node)
        if n == 4:
            ret = "    {} {}\n".format(node(label(4)), label(4)) + "    |\n    |\n"
            ret += "{}---{}---{}\n".format(node(label(1)), node(label(2)), node(label(3)))
            ret += "{!s:4}|{!s:3}{!s:4}\n".format(label(1), label(2), label(3))
            ret += "    |\n    {} {}".format(node(label(0)), label(0))
            return ret

        ret = "{!s:>3} {}".format(label(0), node(label(0)))
        ret += (4*(n-4)-1)*" "+"{} {}\n".format(node(label(n)), label(n))
        ret += "    |" + (4*(n-4)-1)*" " + "|\n"
        ret += "    |" + (4*(n-4)-1)*" " + "|\n"
        ret += "---".join(node(label(i)) for i in range(1, n))
        ret += '\n' + "".join("{!s:4}".format(label(i)) for i in range(1,n))
        return ret

