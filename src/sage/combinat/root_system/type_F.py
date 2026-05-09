"""
Root system data for type F
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
from sage.rings.integer_ring import ZZ
from sage.sets.family import Family

# TODO: double check that this can't be defined over ZZ


class AmbientSpace(ambient_space.AmbientSpace):
    """
    The lattice behind `F_4`.  The computations are based on Bourbaki,
    Groupes et Algèbres de Lie, Ch. 4,5,6 (planche VIII).
    """

    def __init__(self, root_system, base_ring):
        r"""
        Initialize the ambient lattice for the root system of type `F_4`.

        This essentially initializes ``Base`` with the coordinates of
        the simple roots in the canonical basis for `\RR^4`.

        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()

        TESTS::

            sage: TestSuite(e).run()                                                    # needs sage.graphs
        """
        ambient_space.AmbientSpace.__init__(self, root_system, base_ring)
        v = ZZ(1)/ZZ(2)
        self.Base = [self.root(1,2,p2=1),
                     self.root(2,3,p2=1),
                     self.root(3),
                     v*(self.root(0)-self.root(1)-self.root(2)-self.root(3))]

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()
            sage: e.dimension()
            4
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j=None, k=None, l=None, p1=0, p2=0, p3=0, p4=0):
        """
        Compute a root from base elements of the underlying lattice.
        The arguments specify the basis elements and the signs.
        Sadly, the base elements are indexed zero-based.
        We assume that if one of the indices is not given, the rest are not as well.

        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()
            sage: [ e.root(i,j,p2=1) for i in range(e.n) for j in range(i+1,e.n) ]
            [(1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1), (0, 1, -1, 0), (0, 1, 0, -1), (0, 0, 1, -1)]
        """
        if i == j or j is None:
            return (-1)**p1*self.monomial(i)
        if k is None:
            return (-1)**p1*self.monomial(i) + (-1)**p2*self.monomial(j)
        if l is None:
            return (-1)**p1*self.monomial(i) + (-1)**p2*self.monomial(j)+(-1)**p3*self.monomial(k)
        return (-1)**p1*self.monomial(i) + (-1)**p2*self.monomial(j)+(-1)**p3*self.monomial(k)+(-1)**p4*self.monomial(l)

    def simple_root(self, i):
        r"""
        Return the `i`-th simple root.

        It is computed according to what Bourbaki calls the Base:

        .. MATH::

            \alpha_1 = \epsilon_2-\epsilon_3,
            \alpha_2 = \epsilon_3-\epsilon_4,
            \alpha_3 = \epsilon_4,
            \alpha_4 = \frac{1}{2} \left( \epsilon_1-\epsilon_2-\epsilon_3-\epsilon_4 \right).

        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()
            sage: e.simple_roots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 1), 4: (1/2, -1/2, -1/2, -1/2)}
        """
        return self.Base[i-1]

    def negative_roots(self):
        """
        Return the negative roots.

        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()
            sage: e.negative_roots()
            [(-1, 0, 0, 0),
            (0, -1, 0, 0),
            (0, 0, -1, 0),
            (0, 0, 0, -1),
            (-1, -1, 0, 0),
            (-1, 0, -1, 0),
            (-1, 0, 0, -1),
            (0, -1, -1, 0),
            (0, -1, 0, -1),
            (0, 0, -1, -1),
            (-1, 1, 0, 0),
            (-1, 0, 1, 0),
            (-1, 0, 0, 1),
            (0, -1, 1, 0),
            (0, -1, 0, 1),
            (0, 0, -1, 1),
            (-1/2, -1/2, -1/2, -1/2),
            (-1/2, -1/2, -1/2, 1/2),
            (-1/2, -1/2, 1/2, -1/2),
            (-1/2, -1/2, 1/2, 1/2),
            (-1/2, 1/2, -1/2, -1/2),
            (-1/2, 1/2, -1/2, 1/2),
            (-1/2, 1/2, 1/2, -1/2),
            (-1/2, 1/2, 1/2, 1/2)]
        """
        return [ -a for a in self.positive_roots()]

    def positive_roots(self):
        r"""
        Return the positive roots.

        These are the roots which are positive with respect to the
        lexicographic ordering of the basis elements
        (`\epsilon_1<\epsilon_2<\epsilon_3<\epsilon_4`).

        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()
            sage: e.positive_roots()
            [(1, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, 1, 0),
            (0, 0, 0, 1),
            (1, 1, 0, 0),
            (1, 0, 1, 0),
            (1, 0, 0, 1),
            (0, 1, 1, 0),
            (0, 1, 0, 1),
            (0, 0, 1, 1),
            (1, -1, 0, 0),
            (1, 0, -1, 0),
            (1, 0, 0, -1),
            (0, 1, -1, 0),
            (0, 1, 0, -1),
            (0, 0, 1, -1),
            (1/2, 1/2, 1/2, 1/2),
            (1/2, 1/2, 1/2, -1/2),
            (1/2, 1/2, -1/2, 1/2),
            (1/2, 1/2, -1/2, -1/2),
            (1/2, -1/2, 1/2, 1/2),
            (1/2, -1/2, 1/2, -1/2),
            (1/2, -1/2, -1/2, 1/2),
            (1/2, -1/2, -1/2, -1/2)]
            sage: e.rho()
            (11/2, 5/2, 3/2, 1/2)
        """
        v = ZZ(1)/ZZ(2)
        if not hasattr(self, 'PosRoots'):
            self.PosRoots = ([ self.monomial(i) for i in range(self.n) ] +
                            [ self.root(i,j,p2=0) for i in range(self.n) for j in range(i+1,self.n) ] +
                            [ self.root(i,j,p2=1) for i in range(self.n) for j in range(i+1,self.n) ] +
                            [ v*self.root(0,1,2,3,0,p2,p3,p4) for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] ])
        return self.PosRoots

    def fundamental_weights(self):
        """
        Return the fundamental weights of ``self``.

        EXAMPLES::

            sage: e =  RootSystem(['F',4]).ambient_space()
            sage: e.fundamental_weights()
            Finite family {1: (1, 1, 0, 0), 2: (2, 1, 1, 0), 3: (3/2, 1/2, 1/2, 1/2), 4: (1, 0, 0, 0)}
        """
        v = ZZ(1)/ZZ(2)
        return Family({ 1: self.monomial(0)+self.monomial(1),
                        2: 2*self.monomial(0)+self.monomial(1)+self.monomial(2),
                        3: v*(3*self.monomial(0)+self.monomial(1)+self.monomial(2)+self.monomial(3)),
                        4: self.monomial(0)})


from .cartan_type import CartanType_standard_finite, CartanType_simple, CartanType_crystallographic


class CartanType(CartanType_standard_finite, CartanType_simple, CartanType_crystallographic):
    def __init__(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['F',4])
            sage: ct
            ['F', 4]
            sage: ct._repr_(compact = True)
            'F4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            False
            sage: ct.dual()
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
            sage: ct.affine()
            ['F', 4, 1]

        TESTS::

            sage: TestSuite(ct).run()
        """
        CartanType_standard_finite.__init__(self, "F", 4)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['F',4]))
            F_4
            sage: latex(CartanType(['F',4]).dual())
            F_4 \text{ relabelled by } \left\{1 : 4, 2 : 3, 3 : 2, 4 : 1\right\}
        """
        return "F_4"

    AmbientSpace = AmbientSpace

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['F',4]).coxeter_number()
            12
        """
        return 12

    def dual_coxeter_number(self):
        """
        Return the dual Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['F',4]).dual_coxeter_number()
            9
        """
        return 9

    def dynkin_diagram(self):
        """
        Return a Dynkin diagram for type F.

        EXAMPLES::

            sage: f = CartanType(['F',4]).dynkin_diagram(); f                           # needs sage.graphs
            O---O=>=O---O
            1   2   3   4
            F4
            sage: f.edges(sort=True)                                                    # needs sage.graphs
            [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1), (3, 4, 1), (4, 3, 1)]
        """
        from .dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        for i in range(1, 4):
            g.add_edge(i, i+1)
        g.set_edge_label(2,3,2)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['F',4])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (2 cm,0);
            \draw (2 cm, 0.1 cm) -- +(2 cm,0);
            \draw (2 cm, -0.1 cm) -- +(2 cm,0);
            \draw (4.0 cm,0) -- +(2 cm,0);
            \draw[shift={(3.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        ret = "\\draw (0 cm,0) -- (%s cm,0);\n" % node_dist
        ret += "\\draw (%s cm, 0.1 cm) -- +(%s cm,0);\n" % (node_dist, node_dist)
        ret += "\\draw (%s cm, -0.1 cm) -- +(%s cm,0);\n" % (node_dist, node_dist)
        ret += "\\draw (%s cm,0) -- +(%s cm,0);\n" % (node_dist*2.0, node_dist)
        if dual:
            ret += self._latex_draw_arrow_tip(1.5*node_dist-0.2, 0, 180)
        else:
            ret += self._latex_draw_arrow_tip(1.5*node_dist+0.2, 0, 0)
        for i in range(4):
            ret += node(i*node_dist, 0, label(i+1))
        return ret

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['F',4]).ascii_art(label = lambda x: x+2))
            O---O=>=O---O
            3   4   5   6
            sage: print(CartanType(['F',4]).ascii_art(label = lambda x: x-2))
            O---O=>=O---O
            -1  0   1   2
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        ret = "{}---{}=>={}---{}\n".format(node(label(1)), node(label(2)),
                                           node(label(3)), node(label(4)))
        ret += ("{!s:4}"*4).format(label(1), label(2), label(3), label(4))
        return ret

    def dual(self):
        r"""
        Return the dual Cartan type.

        This uses that `F_4` is self-dual up to relabelling.

        EXAMPLES::

            sage: F4 = CartanType(['F',4])
            sage: F4.dual()
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}

            sage: F4.dynkin_diagram()                                                   # needs sage.graphs
            O---O=>=O---O
            1   2   3   4
            F4
            sage: F4.dual().dynkin_diagram()                                            # needs sage.graphs
            O---O=>=O---O
            4   3   2   1
            F4 relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
        """
        return self.relabel({1: 4, 2: 3, 3: 2, 4: 1})

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['F', 4])._default_folded_cartan_type()
            ['F', 4] as a folding of ['E', 6]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        return CartanTypeFolded(self, ['E', 6], [[2], [4], [3, 5], [1, 6]])


# For unpickling backward compatibility (Sage <= 4.1)
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_F',
                           'ambient_space', AmbientSpace)
