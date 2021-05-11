"""
Root system data for type C
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import ambient_space

class AmbientSpace(ambient_space.AmbientSpace):
    """
    EXAMPLES::

        sage: e = RootSystem(['C',2]).ambient_space(); e
        Ambient space of the Root system of type ['C', 2]

    One cannot construct the ambient lattice because the fundamental
    coweights have rational coefficients::

        sage: e.smallest_base_ring()
        Rational Field

        sage: RootSystem(['B',2]).ambient_space().fundamental_weights()
        Finite family {1: (1, 0), 2: (1/2, 1/2)}

    TESTS::

        sage: TestSuite(e).run()
    """


    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['C',3]).ambient_space()
            sage: e.dimension()
            3
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j, p1, p2):
        """
        Note that indexing starts at 0.

        EXAMPLES::

            sage: e = RootSystem(['C',3]).ambient_space()
            sage: e.root(0, 1, 1, 1)
            (-1, -1, 0)
        """
        return (-1)**p1 * self.monomial(i) + (-1)**p2 * self.monomial(j)

    def simple_root(self, i):
        """
        EXAMPLES::

            sage: RootSystem(['C',3]).ambient_space().simple_roots()
            Finite family {1: (1, -1, 0), 2: (0, 1, -1), 3: (0, 0, 2)}
        """
        if i not in self.index_set():
            raise ValueError("{} is not in the index set".format(i))
        return self.root(i-1, i,0,1) if i < self.n else self.root(self.n-1, self.n-1, 0, 0)

    def positive_roots(self):
        """
        EXAMPLES::

            sage: RootSystem(['C',3]).ambient_space().positive_roots()
            [(1, 1, 0),
             (1, 0, 1),
             (0, 1, 1),
             (1, -1, 0),
             (1, 0, -1),
             (0, 1, -1),
             (2, 0, 0),
             (0, 2, 0),
             (0, 0, 2)]
        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend([self.root(i,j,0,p) for i in range(j)])
        res.extend([self.root(i,i,0,0) for i in range(self.n)])
        return res

    def negative_roots(self):
        """
        EXAMPLES::

            sage: RootSystem(['C',3]).ambient_space().negative_roots()
            [(-1, 1, 0),
             (-1, 0, 1),
             (0, -1, 1),
             (-1, -1, 0),
             (-1, 0, -1),
             (0, -1, -1),
             (-2, 0, 0),
             (0, -2, 0),
             (0, 0, -2)]
        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend( [self.root(i,j,1,p) for i in range(j) ] )
        res.extend( [ self.root(i,i,1,1) for i in range(self.n) ] )
        return res


    def fundamental_weight(self, i):
        """
        EXAMPLES::

            sage: RootSystem(['C',3]).ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1)}
        """
        return self.sum(self.monomial(j) for j in range(i))

from .cartan_type import CartanType_standard_finite, CartanType_simple, CartanType_crystallographic, CartanType_simply_laced
class CartanType(CartanType_standard_finite, CartanType_simple, CartanType_crystallographic):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['C',4])
            sage: ct
            ['C', 4]
            sage: ct._repr_(compact = True)
            'C4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            False
            sage: ct.affine()
            ['C', 4, 1]
            sage: ct.dual()
            ['B', 4]

            sage: ct = CartanType(['C',1])
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['C', 1, 1]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 1
        CartanType_standard_finite.__init__(self, "C", n)
        if n == 1:
            self._add_abstract_superclass(CartanType_simply_laced)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['C',4]))
            C_{4}
        """
        return "C_{%s}"%self.n

    AmbientSpace = AmbientSpace

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['C',4]).coxeter_number()
            8
        """
        return 2*self.n

    def dual_coxeter_number(self):
        """
        Return the dual Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['C',4]).dual_coxeter_number()
            5
        """
        return self.n + 1

    def dual(self):
        """
        Types B and C are in duality:

        EXAMPLES::

            sage: CartanType(["C", 3]).dual()
            ['B', 3]
        """
        from . import cartan_type
        return cartan_type.CartanType(["B", self.n])

    def dynkin_diagram(self):
        """
        Returns a Dynkin diagram for type C.

        EXAMPLES::

            sage: c = CartanType(['C',3]).dynkin_diagram()
            sage: c
            O---O=<=O
            1   2   3
            C3
            sage: c.edges(sort=True)
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

             sage: b = CartanType(['C',1]).dynkin_diagram()
             sage: b
             O
             1
             C1
             sage: sorted(b.edges())
             []
        """
        return self.dual().dynkin_diagram().dual()

    def _latex_dynkin_diagram(self, label=lambda x: x, node=None, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['C',4])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(4.8, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>

        When ``dual=True``, the Dynkin diagram for the dual Cartan
        type `B_n` is returned::

            sage: print(CartanType(['C',4])._latex_dynkin_diagram(dual=True))
            \draw (0 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(5.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>

        .. SEEALSO::

            - :meth:`sage.combinat.root_system.type_C.CartanType._latex_dynkin_diagram`
            - :meth:`sage.combinat.root_system.type_BC_affine.CartanType._latex_dynkin_diagram`
        """
        return self.dual()._latex_dynkin_diagram(label=label, node=node, node_dist=node_dist, dual=not dual)

    def ascii_art(self, label=lambda i: i, node=None):
        """
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['C',1]).ascii_art())
            O
            1
            sage: print(CartanType(['C',2]).ascii_art())
            O=<=O
            1   2
            sage: print(CartanType(['C',3]).ascii_art())
            O---O=<=O
            1   2   3
            sage: print(CartanType(['C',5]).ascii_art(label = lambda x: x+2))
            O---O---O---O=<=O
            3   4   5   6   7
        """
        return self.dual().ascii_art(label=label, node=node).replace("=>=", "=<=")

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['C', 3])._default_folded_cartan_type()
            ['C', 3] as a folding of ['A', 5]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        n = self.n
        return CartanTypeFolded(self, ['A', 2*n-1],
            [[i, 2*n-i] for i in range(1, n)] + [[n]])

# For unpickling backward compatibility (Sage <= 4.1)
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_C', 'ambient_space',  AmbientSpace)
