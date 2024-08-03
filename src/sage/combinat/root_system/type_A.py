"""
Root system data for type A
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations
from . import ambient_space


class AmbientSpace(ambient_space.AmbientSpace):
    r"""
    EXAMPLES::

        sage: R = RootSystem(["A",3])
        sage: e = R.ambient_space(); e
        Ambient space of the Root system of type ['A', 3]
        sage: TestSuite(e).run()                                                        # needs sage.graphs

    By default, this ambient space uses the barycentric projection for plotting::

        sage: # needs sage.symbolic
        sage: L = RootSystem(["A",2]).ambient_space()
        sage: e = L.basis()
        sage: L._plot_projection(e[0])
        (1/2, 989/1142)
        sage: L._plot_projection(e[1])
        (-1, 0)
        sage: L._plot_projection(e[2])
        (1/2, -989/1142)
        sage: L = RootSystem(["A",3]).ambient_space()
        sage: l = L.an_element(); l
        (2, 2, 3, 0)
        sage: L._plot_projection(l)
        (0, -1121/1189, 7/3)

    .. SEEALSO::

        - :meth:`sage.combinat.root_system.root_lattice_realizations.RootLatticeRealizations.ParentMethods._plot_projection`
    """
    @classmethod
    def smallest_base_ring(cls, cartan_type=None):
        """
        Return the smallest base ring the ambient space can be defined upon.

        .. SEEALSO:: :meth:`~sage.combinat.root_system.ambient_space.AmbientSpace.smallest_base_ring`

        EXAMPLES::

            sage: e = RootSystem(["A",3]).ambient_space()
            sage: e.smallest_base_ring()
            Integer Ring
        """
        return ZZ

    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(["A",3]).ambient_space()
            sage: e.dimension()
            4
        """
        return self.root_system.cartan_type().rank()+1

    def root(self, i, j):
        """
        Note that indexing starts at 0.

        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.root(0,1)
            (1, -1, 0, 0)
        """
        return self.monomial(i) - self.monomial(j)

    def simple_root(self, i):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1)}
        """
        return self.root(i-1, i)

    def negative_roots(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.negative_roots()
            [(-1, 1, 0, 0),
             (-1, 0, 1, 0),
             (-1, 0, 0, 1),
             (0, -1, 1, 0),
             (0, -1, 0, 1),
             (0, 0, -1, 1)]
        """
        res = []
        for j in range(self.n-1):
            for i in range(j+1,self.n):
                res.append(  self.root(i,j) )
        return res

    def positive_roots(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.positive_roots()
            [(1, -1, 0, 0),
             (1, 0, -1, 0),
             (0, 1, -1, 0),
             (1, 0, 0, -1),
             (0, 1, 0, -1),
             (0, 0, 1, -1)]
        """
        res = []
        for j in range(self.n):
            for i in range(j):
                res.append(  self.root(i,j) )
        return res

    def highest_root(self):
        """
        EXAMPLES::

           sage: e = RootSystem(['A',3]).ambient_lattice()
           sage: e.highest_root()
           (1, 0, 0, -1)
        """
        return self.root(0,self.n-1)

    def fundamental_weight(self, i):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1, 1, 1, 0)}
        """
        return self.sum(self.monomial(j) for j in range(i))

    def det(self, k=1):
        """
        Return the vector (1, ... ,1) which in the ['A',r]
        weight lattice, interpreted as a weight of GL(r+1,CC)
        is the determinant. If the optional parameter k is
        given, returns (k, ... ,k), the `k`-th power of the
        determinant.

        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_space()
            sage: e.det(1/2)
            (1/2, 1/2, 1/2, 1/2)
        """
        return self.sum(self.monomial(j)*k for j in range(self.n))

    _plot_projection = RootLatticeRealizations.ParentMethods.__dict__['_plot_projection_barycentric']


from .cartan_type import CartanType_standard_finite, CartanType_simply_laced, CartanType_simple


class CartanType(CartanType_standard_finite, CartanType_simply_laced, CartanType_simple):
    """
    Cartan Type `A_n`.

    .. SEEALSO:: :func:`~sage.combinat.root_systems.cartan_type.CartanType`
    """

    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4])
            sage: ct
            ['A', 4]
            sage: ct._repr_(compact = True)
            'A4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['A', 4, 1]
            sage: ct.dual()
            ['A', 4]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 0
        CartanType_standard_finite.__init__(self, "A", n)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['A',4]))
            A_{4}
        """
        return "A_{%s}" % self.n

    AmbientSpace = AmbientSpace

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['A',4]).coxeter_number()
            5
        """
        return self.n + 1

    def dual_coxeter_number(self):
        """
        Return the dual Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['A',4]).dual_coxeter_number()
            5
        """
        return self.n + 1

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram of type A.

        EXAMPLES::

            sage: a = CartanType(['A',3]).dynkin_diagram(); a                           # needs sage.graphs
            O---O---O
            1   2   3
            A3
            sage: a.edges(sort=True)                                                    # needs sage.graphs
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 1)]

        TESTS::

            sage: a = DynkinDiagram(['A',1]); a                                         # needs sage.graphs
            O
            1
            A1
            sage: a.vertices(sort=False), a.edges(sort=False)                           # needs sage.graphs
            ([1], [])
        """
        from .dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        for i in range(1, n):
            g.add_edge(i, i+1)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A',4])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (6 cm,0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>

            sage: print(CartanType(['A',0])._latex_dynkin_diagram())
            <BLANKLINE>
            sage: print(CartanType(['A',1])._latex_dynkin_diagram())
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        if self.n > 1:
            ret = "\\draw (0 cm,0) -- ({} cm,0);\n".format((self.n-1)*node_dist)
        else:
            ret = ""
        return ret + "".join(node((i-1)*node_dist, 0, label(i))
                             for i in self.index_set())

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A',0]).ascii_art())
            sage: print(CartanType(['A',1]).ascii_art())
            O
            1
            sage: print(CartanType(['A',3]).ascii_art())
            O---O---O
            1   2   3
            sage: print(CartanType(['A',12]).ascii_art())
            O---O---O---O---O---O---O---O---O---O---O---O
            1   2   3   4   5   6   7   8   9   10  11  12
            sage: print(CartanType(['A',5]).ascii_art(label = lambda x: x+2))
            O---O---O---O---O
            3   4   5   6   7
            sage: print(CartanType(['A',5]).ascii_art(label = lambda x: x-2))
            O---O---O---O---O
            -1  0   1   2   3
        """
        n = self.n
        if n == 0:
            return ""
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        ret = "---".join(node(label(i)) for i in range(1,n+1)) + "\n"
        ret += "".join("{!s:4}".format(label(i)) for i in range(1,n+1))
        return ret


# For unpickling backward compatibility (Sage <= 4.1)
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_A', 'ambient_space',  AmbientSpace)
