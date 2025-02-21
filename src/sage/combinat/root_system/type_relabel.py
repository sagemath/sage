"""
Root system data for relabelled Cartan types
"""
#*****************************************************************************
#       Copyright (C) 2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.family import Family, FiniteFamily
from sage.combinat.root_system import cartan_type
from sage.combinat.root_system import ambient_space
from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations


class CartanType(cartan_type.CartanType_decorator):
    r"""
    A class for relabelled Cartan types.
    """
    @staticmethod
    def __classcall__(cls, type, relabelling):
        """
        This standardizes the input of the constructor to ensure
        unique representation.

        EXAMPLES::

            sage: ct1 = CartanType(['B',2]).relabel({1:2, 2:1})    # indirect doctest
            sage: ct2 = CartanType(['B',2]).relabel(lambda x: 3-x)
            sage: ct3 = CartanType(['B',2]).relabel({1:3, 2: 4})
            sage: ct4 = CartanType(['D',4]).relabel(lambda x: 3-x)
            sage: ct1 == ct2
            True
            sage: ct1 == ct3
            False
            sage: ct1 == ct4
            False
        """
        if isinstance(relabelling, (list, tuple, dict, FiniteFamily)):
            # allows for using relabellings with more entries than in the index_set
            # and by the way makes a copy of relabelling
            relabelling = {i: relabelling[i] for i in type.index_set()}
        else:
            relabelling = {i: relabelling(i) for i in type.index_set()}

        if isinstance(type, CartanType): # type is already a relabelled type
            relabelling = {i: relabelling[type._relabelling[i]]
                           for i in type._type.index_set()}
            type = type._type

        if all( relabelling[i] == i for i in type.index_set() ):
            return type

        relabelling = FiniteFamily(relabelling) # Hack to emulate a frozendict which would be hashable!!!!
        return super().__classcall__(cls, type, relabelling)

    def __init__(self, type, relabelling):
        """
        INPUT:

        - ``type`` -- a Cartan type

        - ``relabelling`` -- a function (or a list, or a dictionary)

        Returns an isomorphic Cartan type obtained by relabelling the
        nodes of the Dynkin diagram. Namely the node with label ``i``
        is relabelled ``f(i)`` (or, by ``f[i]`` if ``f`` is a list or
        dictionary).

        EXAMPLES:

        We take the Cartan type `B_4`::

            sage: T = CartanType(['B',4])
            sage: T.dynkin_diagram()                                                    # needs sage.graphs
            O---O---O=>=O
            1   2   3   4
            B4

        And relabel its nodes::

            sage: cycle = {1:2, 2:3, 3:4, 4:1}

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()                                                    # needs sage.graphs
            O---O---O=>=O
            2   3   4   1
            B4 relabelled by {1: 2, 2: 3, 3: 4, 4: 1}
            sage: T.dynkin_diagram().edges(sort=True)                                   # needs sage.graphs
            [(1, 4, 1), (2, 3, 1), (3, 2, 1), (3, 4, 1), (4, 1, 2), (4, 3, 1)]

        Multiple relabelling are recomposed into a single one::

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()                                                    # needs sage.graphs
            O---O---O=>=O
            3   4   1   2
            B4 relabelled by {1: 3, 2: 4, 3: 1, 4: 2}

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()                                                    # needs sage.graphs
            O---O---O=>=O
            4   1   2   3
            B4 relabelled by {1: 4, 2: 1, 3: 2, 4: 3}

        And trivial relabelling are honoured nicely::

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()                                                    # needs sage.graphs
            O---O---O=>=O
            1   2   3   4
            B4

        TESTS:

        Test that the produced Cartan type is in the appropriate
        abstract classes (see :issue:`13724`)::

            sage: ct = CartanType(['B',4]).relabel(cycle)
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

            sage: ct = CartanType(['A',3,1]).relabel({0:3,1:2, 2:1,3:0})
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

        Check for the original issues of :issue:`13724`::

            sage: A3 = CartanType("A3")
            sage: A3.cartan_matrix()                                                    # needs sage.graphs
            [ 2 -1  0]
            [-1  2 -1]
            [ 0 -1  2]
            sage: A3r = A3.relabel({1:2,2:3,3:1})
            sage: A3r.cartan_matrix()                                                   # needs sage.graphs
            [ 2  0 -1]
            [ 0  2 -1]
            [-1 -1  2]

            sage: ct = CartanType(["D",4,3]).classical(); ct
            ['G', 2]
            sage: ct.symmetrizer()                                                      # needs sage.graphs
            Finite family {1: 1, 2: 3}

        Check the underlying issue of :issue:`24892`, that the root system
        of a relabelled non-crystallographic Cartan type has an
        ``ambient_space()`` that does not result in an error (note that
        this should actually return a valid ambient space, which requires
        the non-crystallographic finite types to have them implemented)::

            sage: rI5 = CartanType(['I',5]).relabel({1:0,2:1})
            sage: rI5.root_system().ambient_space()
        """
        cartan_type.CartanType_decorator.__init__(self, type)
        relabelling = Family(relabelling)
        self._relabelling = dict(relabelling.items())
        self._relabelling_inverse = dict(relabelling.inverse_family().items())
        self._index_set = tuple(sorted(relabelling[i] for i in type.index_set()))
        # TODO: design an appropriate infrastructure to handle this
        # automatically? Maybe using categories and axioms?
        # See also type_dual.CartanType.__init__
        if type.is_finite() and (isinstance(type, cartan_type.SuperCartanType_standard)
                                 or type.is_crystallographic()):
            # FIXME: Remove the is_crystallographic (and the short-circuiting
            #   super) check once the non-crystallographic finite types
            #   (i.e., H_3, H_4, I_2(p)) have an implementation of an
            #   ambient space. See issue #24892.
            self.__class__ = CartanType_finite
        elif type.is_affine():
            self.__class__ = CartanType_affine
        abstract_classes = tuple(cls
                                 for cls in self._stable_abstract_classes
                                 if isinstance(type, cls))
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

           sage: CartanType(['F', 4]).relabel(lambda x: 5-x)
           ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}

           sage: CartanType(['F', 4]).relabel(lambda x: 5-x)._repr_(compact = True)
           'F4 relabelled by {1: 4, 2: 3, 3: 2, 4: 1}'

        TESTS::

            sage: CoxeterType(['I',5]).relabel({1:0,2:1})
            Coxeter type of ['I', 5] relabelled by {1: 0, 2: 1}
        """
        from pprint import pformat
        # Special case for type D_4^3
        if (self._type.is_affine() and self._type.dual().type() == 'G'
                and self.options("notation") == "Kac"):
            if compact:
                return 'D4^3'
            return "['D', 4, 3]"
        relab = pformat(self._relabelling)
        return self._type._repr_(compact=compact) + " relabelled by {}".format(relab)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: ct = CartanType(['A',4]).relabel(lambda x: (x+1)%4+1)
            sage: latex(ct)
            A_{4} \text{ relabelled by } \left\{1 : 3, 2 : 4, 3 : 1, 4 : 2\right\}

        A more compact, but potentially confusing, representation can
        be obtained using the ``latex_relabel`` global option::

            sage: CartanType.options['latex_relabel'] = False
            sage: latex(ct)
            A_{4}
            sage: CartanType.options['latex_relabel'] = True

        Kac's notations are implemented::

            sage: CartanType.options['notation'] = 'Kac'
            sage: latex(CartanType(['D',4,3]))
            D_4^{(3)}
            sage: CartanType.options._reset()

        TESTS::

            sage: latex(CoxeterType(['I',5]).relabel({1:0,2:1}))
            I_2(5) \text{ relabelled by } \left\{1 : 0, 2 : 1\right\}
        """
        from sage.misc.latex import latex
        # Special case for type D_4^{(3)}
        if (self._type.is_affine() and self._type.dual().type() == 'G'
                and self.options("notation") == "Kac"):
            return 'D_4^{(3)}'
        ret = self._type._latex_()
        if self.options('latex_relabel'):
            ret += " \\text{ relabelled by } " + latex(self._relabelling)
        return ret

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A',4]).relabel(lambda x: (x+1)%4+1)._latex_dynkin_diagram())
            \draw (0 cm,0) -- (6 cm,0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda i: i
        return self._type._latex_dynkin_diagram(lambda i: label(self._relabelling[i]), node, node_dist)

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of this Cartan type.

        EXAMPLES::

            sage: print(CartanType(["G", 2]).relabel({1:2,2:1}).ascii_art())
              3
            O=<=O
            2   1
            sage: print(CartanType(["B", 3, 1]).relabel([1,3,2,0]).ascii_art())
                O 1
                |
                |
            O---O=>=O
            3   2   0
            sage: print(CartanType(["F", 4, 1]).relabel(lambda n: 4-n).ascii_art())
            O---O---O=>=O---O
            4   3   2   1   0
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        return self._type.ascii_art(lambda i: label(self._relabelling[i]), node)

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram for this Cartan type.

        EXAMPLES::

            sage: CartanType(["G", 2]).relabel({1:2,2:1}).dynkin_diagram()              # needs sage.graphs
              3
            O=<=O
            2   1
            G2 relabelled by {1: 2, 2: 1}

        TESTS:

        To be compared with the examples in :meth:`ascii_art`::

            sage: CartanType(["G", 2]).relabel({1:2,2:1}).dynkin_diagram().edges(sort=True)         # needs sage.graphs
            [(1, 2, 3), (2, 1, 1)]
            sage: CartanType(["B", 3, 1]).relabel([1,3,2,0]).dynkin_diagram().edges(sort=True)      # needs sage.graphs
            [(0, 2, 1), (1, 2, 1), (2, 0, 2), (2, 1, 1), (2, 3, 1), (3, 2, 1)]
            sage: CartanType(["F", 4, 1]).relabel(lambda n: 4-n).dynkin_diagram().edges(sort=True)  # needs sage.graphs
            [(0, 1, 1), (1, 0, 1), (1, 2, 1), (2, 1, 2), (2, 3, 1), (3, 2, 1), (3, 4, 1), (4, 3, 1)]
        """
        # Maybe we want to move this up as a relabel method for Dynkin diagram
        # We will have to be careful setting the Cartan type of the result though
        result = self._type.dynkin_diagram().copy()
        # relabelling in place allows to keep the extra Dynkin diagram structure
        super(result.__class__, result).relabel(self._relabelling, inplace=True)
        result._cartan_type = self
        return result

    def index_set(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: ct.index_set()
            (1, 2)
        """
        return self._index_set

    def dual(self):
        """
        Implement :meth:`sage.combinat.root_system.cartan_type.CartanType_abstract.dual`,
        using that taking the dual and relabelling are commuting operations.

        EXAMPLES::

            sage: T = CartanType(["BC",3, 2])
            sage: cycle = {1:2, 2:3, 3:0, 0:1}
            sage: T.relabel(cycle).dual().dynkin_diagram()                              # needs sage.graphs
            O=>=O---O=>=O
            1   2   3   0
            BC3~* relabelled by {0: 1, 1: 2, 2: 3, 3: 0}
            sage: T.dual().relabel(cycle).dynkin_diagram()                              # needs sage.graphs
            O=>=O---O=>=O
            1   2   3   0
            BC3~* relabelled by {0: 1, 1: 2, 2: 3, 3: 0}
        """
        return self._type.dual().relabel(self._relabelling)

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: fct = CartanType(['D', 4, 3])._default_folded_cartan_type(); fct
            ['G', 2, 1]^* relabelled by {0: 0, 1: 2, 2: 1} as a folding of ['D', 4, 1]
            sage: fct.folding_orbit()
            Finite family {0: (0,), 1: (2,), 2: (1, 3, 4)}
            sage: CartanType(['G',2,1]).dual()._default_folded_cartan_type().folding_orbit()
            Finite family {0: (0,), 1: (1, 3, 4), 2: (2,)}
            sage: CartanType(['C',3,1]).relabel({0:1, 1:0, 2:3, 3:2}).as_folding().scaling_factors()                    # needs sage.graphs
            Finite family {0: 1, 1: 2, 2: 2, 3: 1}
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        vct = self._type._default_folded_cartan_type()
        sigma = vct.folding_orbit()
        return CartanTypeFolded(self, vct._folding,
            {self._relabelling[i]: sigma[i] for i in self._type.index_set()})

    def type(self):
        """
        Return the type of ``self`` or ``None`` if unknown.

        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: ct.type()
            'G'
        """
        return self._type.type()

    @cached_method
    def coxeter_diagram(self):
        """
        Return the Coxeter diagram for ``self``.

        EXAMPLES::

            sage: ct = CartanType(['H', 3]).relabel({1:3,2:2,3:1})
            sage: G = ct.coxeter_diagram(); G                                           # needs sage.graphs
            Graph on 3 vertices
            sage: G.edges(sort=True)                                                    # needs sage.graphs
            [(1, 2, 5), (2, 3, 3)]
        """
        return self._type.coxeter_diagram().relabel(self._relabelling, inplace=False, immutable=True)

###########################################################################


class AmbientSpace(ambient_space.AmbientSpace):
    """
    Ambient space for a relabelled finite Cartan type.

    It is constructed in the canonical way from the ambient space of
    the original Cartan type, by relabelling the simple roots,
    fundamental weights, etc.

    EXAMPLES::

        sage: cycle = {1:2, 2:3, 3:4, 4:1}
        sage: L = CartanType(["F",4]).relabel(cycle).root_system().ambient_space(); L
        Ambient space of the Root system of type ['F', 4] relabelled by {1: 2, 2: 3, 3: 4, 4: 1}
        sage: TestSuite(L).run()                                                        # needs sage.graphs
    """

    @lazy_attribute
    def _space(self):
        """
        The ambient space this is a relabelling of.

        EXAMPLES::

            sage: cycle = {1:2, 2:3, 3:4, 4:1}
            sage: L = CartanType(["F",4]).relabel(cycle).root_system().ambient_space()
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

            sage: cycle = {1:2, 2:3, 3:4, 4:1}
            sage: L = CartanType(["F",4]).relabel(cycle).root_system().ambient_space()
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

            sage: cycle = {1:2, 2:3, 3:4, 4:1}
            sage: L = CartanType(["F",4]).relabel(cycle).root_system().ambient_space()
            sage: K = CartanType(["F",4]).root_system().ambient_space()
            sage: K.simple_roots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 1), 4: (1/2, -1/2, -1/2, -1/2)}
            sage: K.simple_coroots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 2), 4: (1, -1, -1, -1)}
            sage: L.simple_root(1)
            (1/2, -1/2, -1/2, -1/2)

            sage: L.simple_roots()
            Finite family {1: (1/2, -1/2, -1/2, -1/2), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1), 4: (0, 0, 0, 1)}

            sage: L.simple_coroots()
            Finite family {1: (1, -1, -1, -1), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1), 4: (0, 0, 0, 2)}
        """
        i = self.cartan_type()._relabelling_inverse[i]
        return self.sum_of_terms(self._space.simple_root(i))

    @cached_method
    def fundamental_weight(self, i):
        """
        Return the ``i``-th fundamental weight.

        It is constructed by looking up the corresponding simple
        coroot in the ambient space for the original Cartan type.

        EXAMPLES::

            sage: cycle = {1:2, 2:3, 3:4, 4:1}
            sage: L = CartanType(["F",4]).relabel(cycle).root_system().ambient_space()
            sage: K = CartanType(["F",4]).root_system().ambient_space()
            sage: K.fundamental_weights()
            Finite family {1: (1, 1, 0, 0), 2: (2, 1, 1, 0), 3: (3/2, 1/2, 1/2, 1/2), 4: (1, 0, 0, 0)}
            sage: L.fundamental_weight(1)
            (1, 0, 0, 0)
            sage: L.fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (2, 1, 1, 0), 4: (3/2, 1/2, 1/2, 1/2)}
        """
        i = self.cartan_type()._relabelling_inverse[i]
        return self.sum_of_terms(self._space.fundamental_weight(i))

    @lazy_attribute
    def _plot_projection(self):
        """
        A hack so that if an ambient space uses barycentric projection, then so does its dual.

        EXAMPLES::

            sage: cycle = {1:2, 2:1}
            sage: L = CartanType(["G",2]).relabel(cycle).root_system().ambient_space()
            sage: L._plot_projection == L._plot_projection_barycentric
            True

            sage: cycle = {1:2, 2:3, 3:4, 4:1}
            sage: L = CartanType(["F",4]).relabel(cycle).root_system().ambient_space()
            sage: L._plot_projection == L._plot_projection_barycentric
            False
        """
        if self._space._plot_projection == self._space._plot_projection_barycentric:
            return self._plot_projection_barycentric
        else:
            RootLatticeRealizations.ParentMethods.__dict__["_plot_projection"]


class CartanType_finite(CartanType, cartan_type.CartanType_finite):
    AmbientSpace = AmbientSpace

    def affine(self):
        """
        Return the affine Cartan type associated with ``self``.

        EXAMPLES::

            sage: B4 = CartanType(['B',4])
            sage: B4.dynkin_diagram()                                                   # needs sage.graphs
            O---O---O=>=O
            1   2   3   4
            B4
            sage: B4.affine().dynkin_diagram()                                          # needs sage.graphs
                O 0
                |
                |
            O---O---O=>=O
            1   2   3   4
            B4~

        If possible, this reuses the original label for the special node::

            sage: T = B4.relabel({1:2, 2:3, 3:4, 4:1}); T.dynkin_diagram()              # needs sage.graphs
            O---O---O=>=O
            2   3   4   1
            B4 relabelled by {1: 2, 2: 3, 3: 4, 4: 1}
            sage: T.affine().dynkin_diagram()                                           # needs sage.graphs
                O 0
                |
                |
            O---O---O=>=O
            2   3   4   1
            B4~ relabelled by {0: 0, 1: 2, 2: 3, 3: 4, 4: 1}

        Otherwise, it chooses a label for the special_node in `0,1,...`::

            sage: T = B4.relabel({1:0, 2:1, 3:2, 4:3}); T.dynkin_diagram()              # needs sage.graphs
            O---O---O=>=O
            0   1   2   3
            B4 relabelled by {1: 0, 2: 1, 3: 2, 4: 3}
            sage: T.affine().dynkin_diagram()                                           # needs sage.graphs
                O 4
                |
                |
            O---O---O=>=O
            0   1   2   3
            B4~ relabelled by {0: 4, 1: 0, 2: 1, 3: 2, 4: 3}

        This failed before :issue:`13724`::

            sage: ct = CartanType(["G",2]).dual(); ct
            ['G', 2] relabelled by {1: 2, 2: 1}
            sage: ct.affine()
            ['G', 2, 1] relabelled by {0: 0, 1: 2, 2: 1}

            sage: ct = CartanType(["F",4]).dual(); ct
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
            sage: ct.affine()
            ['F', 4, 1] relabelled by {0: 0, 1: 4, 2: 3, 3: 2, 4: 1}

        Check that we don't inadvertently change the internal
        relabelling of ``ct``::

            sage: ct
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
        """
        affine = self._type.affine()
        relabelling = self._relabelling.copy()
        for special_node in [affine.special_node()] + list(range(affine.rank())):
            if special_node not in self._relabelling_inverse:
                relabelling[affine.special_node()] = special_node
                break
        return self._type.affine().relabel(relabelling)

###########################################################################


class CartanType_affine(CartanType, cartan_type.CartanType_affine):
    """
    TESTS::

        sage: ct = CartanType(['D',4,3]); ct
        ['G', 2, 1]^* relabelled by {0: 0, 1: 2, 2: 1}

        sage: L = ct.root_system().ambient_space(); L
        Ambient space of the Root system of type ['G', 2, 1]^* relabelled by {0: 0, 1: 2, 2: 1}
        sage: L.classical()
        Ambient space of the Root system of type ['G', 2]
        sage: TestSuite(L).run()
    """

    def classical(self):
        """
        Return the classical Cartan type associated with ``self``.

        EXAMPLES::

            sage: A41 = CartanType(['A',4,1])
            sage: A41.dynkin_diagram()                                                  # needs sage.graphs
            0
            O-----------+
            |           |
            |           |
            O---O---O---O
            1   2   3   4
            A4~

            sage: T = A41.relabel({0:1, 1:2, 2:3, 3:4, 4:0})
            sage: T
            ['A', 4, 1] relabelled by {0: 1, 1: 2, 2: 3, 3: 4, 4: 0}
            sage: T.dynkin_diagram()                                                    # needs sage.graphs
            1
            O-----------+
            |           |
            |           |
            O---O---O---O
            2   3   4   0
            A4~ relabelled by {0: 1, 1: 2, 2: 3, 3: 4, 4: 0}

            sage: T0 = T.classical()
            sage: T0
            ['A', 4] relabelled by {1: 2, 2: 3, 3: 4, 4: 0}
            sage: T0.dynkin_diagram()                                                   # needs sage.graphs
            O---O---O---O
            2   3   4   0
            A4 relabelled by {1: 2, 2: 3, 3: 4, 4: 0}
        """
        return self._type.classical().relabel(self._relabelling)

    def basic_untwisted(self):
        r"""
        Return the basic untwisted Cartan type associated with this affine
        Cartan type.

        Given an affine type `X_n^{(r)}`, the basic untwisted type is `X_n`.
        In other words, it is the classical Cartan type that is twisted to
        obtain ``self``.

        EXAMPLES::

            sage: ct = CartanType(['A', 5, 2]).relabel({0:1, 1:0, 2:2, 3:3})
            sage: ct.basic_untwisted()
            ['A', 5]
        """
        return self._type.basic_untwisted()

    def special_node(self):
        r"""
        Return a special node of the Dynkin diagram.

        .. SEEALSO:: :meth:`~sage.combinat.root_system.CartanType_affine.special_node`

        It is obtained by relabelling of the special node of the non
        relabelled Dynkin diagram.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).special_node()
            0
            sage: CartanType(['B', 3, 1]).relabel({1:2, 2:3, 3:0, 0:1}).special_node()
            1
        """
        return self._relabelling[self._type.special_node()]

    def is_untwisted_affine(self):
        """
        Implement :meth:`CartanType_affine.is_untwisted_affine`.

        A relabelled Cartan type is untwisted affine if the original is.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).relabel({1:2, 2:3, 3:0, 0:1}).is_untwisted_affine()
            True
        """
        return self._type.is_untwisted_affine()
