"""
Root system data for Cartan types with marked nodes
"""
#*****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.root_system import cartan_type
from sage.combinat.root_system import ambient_space
from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

class CartanType(cartan_type.CartanType_decorator):
    r"""
    A class for Cartan types with marked nodes.

    INPUT:

    - ``ct`` -- a Cartan type

    - ``marked_nodes`` -- a list of marked nodes

    EXAMPLES:

    We take the Cartan type `B_4`::

        sage: T = CartanType(['B',4])
        sage: T.dynkin_diagram()
        O---O---O=>=O
        1   2   3   4
        B4

    And mark some of its nodes::

        sage: T = T.marked_nodes([2,3])
        sage: T.dynkin_diagram()
        O---X---X=>=O
        1   2   3   4
        B4 with nodes (2, 3) marked

    Markings are not additive::

        sage: T.marked_nodes([1,4]).dynkin_diagram()
        X---O---O=>=X
        1   2   3   4
        B4 with nodes (1, 4) marked

    And trivial relabelling are honoured nicely::

        sage: T = T.marked_nodes([])
        sage: T.dynkin_diagram()
        O---O---O=>=O
        1   2   3   4
        B4
    """
    @staticmethod
    def __classcall__(cls, ct, marked_nodes):
        """
        This standardizes the input of the constructor to ensure
        unique representation.

        EXAMPLES::

            sage: ct1 = CartanType(['B',2]).marked_nodes([1,2])
            sage: ct2 = CartanType(['B',2]).marked_nodes([2,1])
            sage: ct3 = CartanType(['B',2]).marked_nodes((1,2))
            sage: ct4 = CartanType(['D',4]).marked_nodes([1,2])
            sage: ct1 is ct2
            True
            sage: ct1 is ct3
            True
            sage: ct1 == ct4
            False
        """
        ct = cartan_type.CartanType(ct)
        if not marked_nodes:
            return ct
        if any(node not in ct.index_set() for node in marked_nodes):
            raise ValueError("invalid marked node")
        marked_nodes = tuple(sorted(marked_nodes))
        return super(CartanType, cls).__classcall__(cls, ct, marked_nodes)

    def __init__(self, ct, marked_nodes):
        """
        Return an isomorphic Cartan type obtained by marking the
        nodes of the Dynkin diagram.

        TESTS:

        Test that the produced Cartan type is in the appropriate
        abstract classes::

            sage: ct = CartanType(['B',4]).marked_nodes([1,2])
            sage: TestSuite(ct).run()
            sage: from sage.combinat.root_system import cartan_type
            sage: isinstance(ct, cartan_type.CartanType_finite)
            True
            sage: isinstance(ct, cartan_type.CartanType_simple)
            True
            sage: isinstance(ct, cartan_type.CartanType_affine)
            False
            sage: isinstance(ct, cartan_type.CartanType_crystallographic)
            True
            sage: isinstance(ct, cartan_type.CartanType_simply_laced)
            False

            sage: ct = CartanType(['A',3,1]).marked_nodes([1,2])
            sage: TestSuite(ct).run()
            sage: isinstance(ct, cartan_type.CartanType_simple)
            True
            sage: isinstance(ct, cartan_type.CartanType_finite)
            False
            sage: isinstance(ct, cartan_type.CartanType_affine)
            True
            sage: isinstance(ct, cartan_type.CartanType_crystallographic)
            True
            sage: isinstance(ct, cartan_type.CartanType_simply_laced)
            True
        """
        cartan_type.CartanType_decorator.__init__(self, ct)
        self._marked_nodes = marked_nodes
        # TODO: design an appropriate infrastructure to handle this
        # automatically? Maybe using categories and axioms?
        # See also type_dual.CartanType.__init__
        if ct.is_finite():
            self.__class__ = CartanType_finite
        elif ct.is_affine():
            self.__class__ = CartanType_affine
        abstract_classes = tuple(cls
                                 for cls in self._stable_abstract_classes
                                 if isinstance(ct, cls))
        if abstract_classes:
            self._add_abstract_superclass(abstract_classes)

    # For each class cls in _stable_abstract_classes, if ct is an
    # instance of A then ct.relabel(...) is put in this class as well.
    # The order is relevant to avoid MRO issues!
    _stable_abstract_classes = [
        cartan_type.CartanType_finite,
        cartan_type.CartanType_affine,
        cartan_type.CartanType_simple,
        cartan_type.CartanType_simply_laced,
        cartan_type.CartanType_crystallographic]

    def _repr_(self, compact=False):
        """
        EXAMPLES::

           sage: CartanType(['F', 4]).marked_nodes([2])
           ['F', 4] with node 2 marked

           sage: CartanType(['F', 4, 1]).dual().marked_nodes([0, 2])
           ['F', 4, 1]^* with nodes (0, 2) marked

           sage: CartanType(['F', 4, 1]).marked_nodes([0, 2])._repr_(compact = True)
           'F4~ with nodes (0, 2) marked'

            sage: D = DynkinDiagram("A2")
            sage: D.marked_nodes([1])
            O---O
            1   2
            A2 with node 1 marked

            sage: CM = CartanMatrix([[2,-4],[-5,2]])
            sage: CM.marked_nodes([1])
            [ 2 -4]
            [-5  2] with node 1 marked
        """
        if not compact:
            base = repr(self._type)
        else:
            try:
                base = self._type._repr_(compact=True)
            except TypeError:
                base = repr(self._type)
        if len(self._marked_nodes) == 1:
            return base + " with node {} marked".format(self._marked_nodes[0])
        return base + " with nodes {} marked".format(self._marked_nodes)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: ct = CartanType(['A',4]).marked_nodes([1, 3])
            sage: latex(ct)
            A_{4} \text{ with nodes $\left(1, 3\right)$ marked}

        A more compact, but potentially confusing, representation can
        be obtained using the ``latex_marked`` global option::

            sage: CartanType.options['latex_marked'] = False
            sage: latex(ct)
            A_{4}
            sage: CartanType.options['latex_marked'] = True

        Kac's notations are implemented::

            sage: CartanType.options['notation'] = 'Kac'
            sage: latex(CartanType(['D',4,3]).marked_nodes([0]))
            D_4^{(3)} \text{ with node $0$ marked}
            sage: CartanType.options._reset()
        """
        from sage.misc.latex import latex
        ret = self._type._latex_()
        if self.options('latex_marked'):
            if len(self._marked_nodes) == 1:
                ret += " \\text{{ with node ${}$ marked}} ".format(latex(self._marked_nodes[0]))
            else:
                ret += " \\text{{ with nodes ${}$ marked}} ".format(latex(self._marked_nodes))
        return ret

    def _ascii_art_node(self, label):
        """
        Return the ascii art for the node labeled by ``label``.

        EXAMPLES::

            sage: ct = CartanType(['A',4]).marked_nodes([1, 2])
            sage: ct._ascii_art_node(1)
            'X'
            sage: ct._ascii_art_node(3)
            'O'
            sage: CartanType.options._reset()
        """
        if label in self._marked_nodes:
            return self.options('marked_node_str')
        return 'O'

    def _latex_draw_node(self, x, y, label, position="below=4pt", fill='white'):
        r"""
        Draw (possibly marked [crossed out]) circular node ``i`` at the
        position ``(x,y)`` with node label ``label`` .

        - ``position`` -- position of the label relative to the node
        - ``anchor`` -- (optional) the anchor point for the label

        EXAMPLES::

            sage: CartanType.options(mark_special_node='both')
            sage: CartanType(['A',3,1]).marked_nodes([1,3])._latex_draw_node(0, 0, 0)
            '\\draw[fill=black] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};\n'
            sage: CartanType.options._reset()
        """
        ret = cartan_type.CartanType_abstract._latex_draw_node(self, x, y, label, position, fill)
        if label in self._marked_nodes:
            ret += self._latex_draw_mark(x, y)
        return ret

    def _latex_draw_mark(self, x, y, color='black', thickness='thin'):
        r"""
        Draw a mark as a cross `\times` at the point ``(x, y)``.

        INPUT:

        - ``(x, y)`` -- the coordinates of a point, in cm

        - ``color`` -- the color of the mark

        This is an internal function used to assist drawing marked points of
        the Dynkin diagrams. See e.g.
        :meth:`~sage.combinat.root_system.type_marked.CartanType._latex_dynkin_diagram`.

        EXAMPLES::

            sage: print(CartanType(['B',2]).marked_nodes([1,2])._latex_draw_mark(1, 0))
            \draw[shift={(1, 0)}, black, thin] (0.25cm, 0.25cm) -- (-0.25cm, -0.25cm);
            \draw[shift={(1, 0)}, black, thin] (0.25cm, -0.25cm) -- (-0.25cm, 0.25cm);
            <BLANKLINE>
        """
        ret = "\\draw[shift={{({}, {})}}, {}, {}] (0.25cm, 0.25cm) -- (-0.25cm, -0.25cm);\n".format(x, y, color, thickness)
        ret += "\\draw[shift={{({}, {})}}, {}, {}] (0.25cm, -0.25cm) -- (-0.25cm, 0.25cm);\n".format(x, y, color, thickness)
        return ret

    def _latex_dynkin_diagram(self, label=lambda i: i, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A',4]).marked_nodes([1,3])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (6 cm,0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[shift={(0, 0)}, black, thin] (0.25cm, 0.25cm) -- (-0.25cm, -0.25cm);
            \draw[shift={(0, 0)}, black, thin] (0.25cm, -0.25cm) -- (-0.25cm, 0.25cm);
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[shift={(4, 0)}, black, thin] (0.25cm, 0.25cm) -- (-0.25cm, -0.25cm);
            \draw[shift={(4, 0)}, black, thin] (0.25cm, -0.25cm) -- (-0.25cm, 0.25cm);
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>
        """
        if node is None:
            node = self._latex_draw_node
        return self._type._latex_dynkin_diagram(label, node, node_dist)

    def ascii_art(self, label=lambda i: i, node=None):
        """
        Return an ascii art representation of this Cartan type.

        EXAMPLES::

            sage: print(CartanType(["G", 2]).marked_nodes([2]).ascii_art())
              3
            O=<=X
            1   2
            sage: print(CartanType(["B", 3, 1]).marked_nodes([0, 3]).ascii_art())
                X 0
                |
                |
            O---O=>=X
            1   2   3
            sage: print(CartanType(["F", 4, 1]).marked_nodes([0, 2]).ascii_art())
            X---O---X=>=O---O
            0   1   2   3   4
        """
        if node is None:
            node = self._ascii_art_node
        return self._type.ascii_art(label, node)

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram for this Cartan type.

        EXAMPLES::

            sage: CartanType(["G", 2]).marked_nodes([2]).dynkin_diagram()
              3
            O=<=X
            1   2
            G2 with node 2 marked

        TESTS:

        To be compared with the examples in :meth:`ascii_art`::

            sage: sorted(CartanType(["G", 2]).relabel({1:2,2:1}).dynkin_diagram().edges())
            [(1, 2, 3), (2, 1, 1)]
            sage: sorted(CartanType(["B", 3, 1]).relabel([1,3,2,0]).dynkin_diagram().edges())
            [(0, 2, 1), (1, 2, 1), (2, 0, 2), (2, 1, 1), (2, 3, 1), (3, 2, 1)]
            sage: sorted(CartanType(["F", 4, 1]).relabel(lambda n: 4-n).dynkin_diagram().edges())
            [(0, 1, 1), (1, 0, 1), (1, 2, 1), (2, 1, 2), (2, 3, 1), (3, 2, 1), (3, 4, 1), (4, 3, 1)]
        """
        result = self._type.dynkin_diagram().copy()
        result._cartan_type = self
        return result

    def dual(self):
        """
        Implements
        :meth:`sage.combinat.root_system.cartan_type.CartanType_abstract.dual`,
        using that taking the dual and marking nodes are commuting operations.

        EXAMPLES::

            sage: T = CartanType(["BC",3, 2])
            sage: T.marked_nodes([1,3]).dual().dynkin_diagram()
            O=>=X---O=>=X
            0   1   2   3
            BC3~* with nodes (1, 3) marked
            sage: T.dual().marked_nodes([1,3]).dynkin_diagram()
            O=>=X---O=>=X
            0   1   2   3
            BC3~* with nodes (1, 3) marked
        """
        return self._type.dual().marked_nodes(self._marked_nodes)

    def relabel(self, relabelling):
        """
        Return the relabelling of ``self``.

        EXAMPLES::

            sage: T = CartanType(["BC",3, 2])
            sage: T.marked_nodes([1,3]).relabel(lambda x: x+2).dynkin_diagram()
            O=<=X---O=<=X
            2   3   4   5
            BC3~ relabelled by {0: 2, 1: 3, 2: 4, 3: 5} with nodes (3, 5) marked
            sage: T.relabel(lambda x: x+2).marked_nodes([3,5]).dynkin_diagram()
            O=<=X---O=<=X
            2   3   4   5
            BC3~ relabelled by {0: 2, 1: 3, 2: 4, 3: 5} with nodes (3, 5) marked
        """
        rct = self._type.relabel(relabelling)
        rd = rct._relabelling
        marked_nodes = [rd[node] for node in self._marked_nodes]
        return rct.marked_nodes(marked_nodes)

    def marked_nodes(self, marked_nodes):
        """
        Return ``self`` with nodes ``marked_nodes`` marked.

        EXAMPLES::

            sage: ct = CartanType(['A',12])
            sage: m = ct.marked_nodes([1,4,6,7,8,12]); m
            ['A', 12] with nodes (1, 4, 6, 7, 8, 12) marked
            sage: m.marked_nodes([2])
            ['A', 12] with node 2 marked
            sage: m.marked_nodes([]) is ct
            True
        """
        if not marked_nodes:
            return self._type
        return CartanType(self._type, marked_nodes)

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        If a node `a` is marked, then all nodes in the orbit of `a` are marked
        in the ambient type.

        EXAMPLES::

            sage: fct = CartanType(['D', 4, 3])._default_folded_cartan_type(); fct
            ['G', 2, 1]^* relabelled by {0: 0, 1: 2, 2: 1} as a folding of ['D', 4, 1]
            sage: fct.folding_orbit()
            Finite family {0: (0,), 1: (2,), 2: (1, 3, 4)}
            sage: CartanType(['G',2,1]).dual()._default_folded_cartan_type().folding_orbit()
            Finite family {0: (0,), 1: (1, 3, 4), 2: (2,)}
            sage: CartanType(['C',3,1]).relabel({0:1, 1:0, 2:3, 3:2}).as_folding().scaling_factors()
            Finite family {0: 1, 1: 2, 2: 2, 3: 1}
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        vct = self._type._default_folded_cartan_type()
        sigma = vct.folding_orbit()
        marked_nodes = sum([sigma[i] for i in self._marked_nodes], ())
        folding = vct._folding.marked_nodes(marked_nodes)
        return CartanTypeFolded(self, folding, sigma)

    def type(self):
        """
        Return the type of ``self`` or ``None`` if unknown.

        EXAMPLES::

            sage: ct = CartanType(['F', 4]).marked_nodes([1,3])
            sage: ct.type()
            'F'
        """
        return self._type.type()

###########################################################################

class AmbientSpace(ambient_space.AmbientSpace):
    """
    Ambient space for a marked finite Cartan type.

    It is constructed in the canonical way from the ambient space of
    the original Cartan type.

    EXAMPLES::

        sage: L = CartanType(["F",4]).marked_nodes([1,3]).root_system().ambient_space(); L
        Ambient space of the Root system of type ['F', 4] with nodes (1, 3) marked
        sage: TestSuite(L).run()
    """
    @lazy_attribute
    def _space(self):
        """
        The ambient space this is a marking of.

        EXAMPLES::

            sage: L = CartanType(["F",4]).marked_nodes([1,3]).root_system().ambient_space()
            sage: L._space
            Ambient space of the Root system of type ['F', 4]
        """
        K = self.base_ring()
        return self.cartan_type()._type.root_system().ambient_space(K)

    def dimension(self):
        """
        Return the dimension of this ambient space.

        .. SEEALSO:: :meth:`sage.combinat.root_system.ambient_space.AmbientSpace.dimension`

        EXAMPLES::

            sage: L = CartanType(["F",4]).marked_nodes([1,3]).root_system().ambient_space()
            sage: L.dimension()
            4
        """
        # Can't yet use _dual_space for the base ring (and cartan_type?) is not yet initialized
        return self.root_system.cartan_type()._type.root_system().ambient_space().dimension()

    @cached_method
    def simple_root(self, i):
        """
        Return the ``i``-th simple root.

        It is constructed by looking up the corresponding simple
        coroot in the ambient space for the original Cartan type.

        EXAMPLES::

            sage: L = CartanType(["F",4]).marked_nodes([1,3]).root_system().ambient_space()
            sage: L.simple_root(1)
            (0, 1, -1, 0)
            sage: L.simple_roots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1),
                           3: (0, 0, 0, 1), 4: (1/2, -1/2, -1/2, -1/2)}
            sage: L.simple_coroots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1),
                           3: (0, 0, 0, 2), 4: (1, -1, -1, -1)}
        """
        return self.sum_of_terms(self._space.simple_root(i))

    @cached_method
    def fundamental_weight(self, i):
        """
        Return the ``i``-th fundamental weight.

        It is constructed by looking up the corresponding simple
        coroot in the ambient space for the original Cartan type.

        EXAMPLES::

            sage: L = CartanType(["F",4]).marked_nodes([1,3]).root_system().ambient_space()
            sage: L.fundamental_weight(1)
            (1, 1, 0, 0)
            sage: L.fundamental_weights()
            Finite family {1: (1, 1, 0, 0), 2: (2, 1, 1, 0),
                           3: (3/2, 1/2, 1/2, 1/2), 4: (1, 0, 0, 0)}
        """
        return self.sum_of_terms(self._space.fundamental_weight(i))

    @lazy_attribute
    def _plot_projection(self):
        """
        A hack so that if an ambient space uses barycentric projection,
        then so does its dual.

        EXAMPLES::

            sage: L = CartanType(["G",2]).marked_nodes([1]).root_system().ambient_space()
            sage: L._plot_projection == L._plot_projection_barycentric
            True

            sage: L = CartanType(["F",4]).marked_nodes([1,3]).root_system().ambient_space()
            sage: L._plot_projection == L._plot_projection_barycentric
            False
        """
        if self._space._plot_projection == self._space._plot_projection_barycentric:
            return self._plot_projection_barycentric
        else:
            RootLatticeRealizations.ParentMethods.__dict__["_plot_projection"]

###########################################################################

class CartanType_finite(CartanType, cartan_type.CartanType_finite):
    AmbientSpace = AmbientSpace

    def affine(self):
        """
        Return the affine Cartan type associated with ``self``.

        EXAMPLES::

            sage: B4 = CartanType(['B',4]).marked_nodes([1,3])
            sage: B4.dynkin_diagram()
            X---O---X=>=O
            1   2   3   4
            B4 with nodes (1, 3) marked
            sage: B4.affine().dynkin_diagram()
                O 0
                |
                |
            X---O---X=>=O
            1   2   3   4
            B4~ with nodes (1, 3) marked

        TESTS:

        Check that we don't inadvertently change the internal
        marking of ``ct``::

            sage: ct = CartanType(['F', 4]).marked_nodes([1,3])
            sage: ct.affine()
            ['F', 4, 1] with nodes (1, 3) marked
            sage: ct
            ['F', 4] with nodes (1, 3) marked
        """
        return self._type.affine().marked_nodes(self._marked_nodes)

###########################################################################
class CartanType_affine(CartanType, cartan_type.CartanType_affine):
    """
    TESTS::

        sage: ct = CartanType(['B',3,1]).marked_nodes([1,3])
        sage: ct
        ['B', 3, 1] with nodes (1, 3) marked

        sage: L = ct.root_system().ambient_space(); L
        Ambient space of the Root system of type ['B', 3, 1] with nodes (1, 3) marked
        sage: L.classical()
        Ambient space of the Root system of type ['B', 3] with nodes (1, 3) marked
        sage: TestSuite(L).run()
    """
    def _latex_draw_node(self, x, y, label, position="below=4pt"):
        r"""
        Draw the possibly marked (crossed out) circular node ``i`` at the
        position ``(x,y)`` with node label ``label`` .

        - ``position`` -- position of the label relative to the node
        - ``anchor`` -- (optional) the anchor point for the label

        EXAMPLES::

            sage: CartanType.options(mark_special_node='both')
            sage: print(CartanType(['A',3,1]).marked_nodes([0,1,3])._latex_draw_node(0, 0, 0))
            \draw[fill=black] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
            \draw[shift={(0, 0)}, lightgray, very thick] (0.25cm, 0.25cm) -- (-0.25cm, -0.25cm);
            \draw[shift={(0, 0)}, lightgray, very thick] (0.25cm, -0.25cm) -- (-0.25cm, 0.25cm);
            <BLANKLINE>
            sage: CartanType.options._reset()
        """
        mark_special = (label == self.special_node()
                        and self.options('mark_special_node') in ['latex', 'both'])
        if mark_special:
            fill = 'black'
        else:
            fill = 'white'

        ret = cartan_type.CartanType_abstract._latex_draw_node(self, x, y, label, position, fill)

        if label in self._marked_nodes:
            if mark_special:
                ret += self._latex_draw_mark(x, y, 'lightgray', 'very thick')
            else:
                ret += self._latex_draw_mark(x, y)
        return ret

    def _ascii_art_node(self, label):
        """
        Return the ascii art for the node labeled by ``label``.

        EXAMPLES::

            sage: ct = CartanType(['A',4, 1]).marked_nodes([0, 2])
            sage: CartanType.options(mark_special_node='both')
            sage: ct._ascii_art_node(0)
            '#'
            sage: CartanType.options._reset()
        """
        if label in self._marked_nodes:
            if (label == self.special_node()
                    and self.options('mark_special_node') in ['printing', 'both']):
                return '#'
            return self.options('marked_node_str')
        return 'O'

    def classical(self):
        """
        Return the classical Cartan type associated with ``self``.

        EXAMPLES::

            sage: T = CartanType(['A',4,1]).marked_nodes([0,2,4])
            sage: T.dynkin_diagram()
            0
            X-----------+
            |           |
            |           |
            O---X---O---X
            1   2   3   4
            A4~ with nodes (0, 2, 4) marked

            sage: T0 = T.classical()
            sage: T0
            ['A', 4] with nodes (2, 4) marked
            sage: T0.dynkin_diagram()
            O---X---O---X
            1   2   3   4
            A4 with nodes (2, 4) marked
        """
        if self._type.special_node() in self._marked_nodes:
            marked_nodes = list(self._marked_nodes)
            marked_nodes.remove(self._type.special_node())
            return self._type.classical().marked_nodes(marked_nodes)
        return self._type.classical().marked_nodes(self._marked_nodes)

    def basic_untwisted(self):
        r"""
        Return the basic untwisted Cartan type associated with this affine
        Cartan type.

        Given an affine type `X_n^{(r)}`, the basic untwisted type is `X_n`.
        In other words, it is the classical Cartan type that is twisted to
        obtain ``self``.

        EXAMPLES::

            sage: CartanType(['A', 7, 2]).marked_nodes([1,3]).basic_untwisted()
            ['A', 7] with nodes (1, 3) marked
            sage: CartanType(['D', 4, 3]).marked_nodes([0,2]).basic_untwisted()
            ['D', 4] with node 2 marked
        """
        if self._type.special_node() in self._marked_nodes:
            marked_nodes = list(self._marked_nodes)
            marked_nodes.remove(self._type.special_node())
            return self._type.basic_untwisted().marked_nodes(marked_nodes)
        return self._type.basic_untwisted().marked_nodes(self._marked_nodes)

    def special_node(self):
        r"""
        Return the special node of the Cartan type.

        .. SEEALSO:: :meth:`~sage.combinat.root_system.CartanType_affine.special_node`

        It is the special node of the non-marked Cartan type..

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).marked_nodes([1,3]).special_node()
            0
        """
        return self._type.special_node()

    def is_untwisted_affine(self):
        """
        Implement :meth:`CartanType_affine.is_untwisted_affine`.

        A marked Cartan type is untwisted affine if the original is.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).marked_nodes([1,3]).is_untwisted_affine()
            True
        """
        return self._type.is_untwisted_affine()
