"""
Root system data for dual Cartan types
"""
# ****************************************************************************
#       Copyright (C) 2008-2009 Anne Schilling <anne at math.ucdavis.edu>
#       Copyright (C) 2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.call import attrcall
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.root_system import cartan_type
from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations
from sage.combinat.root_system import ambient_space


class CartanType(cartan_type.CartanType_decorator, cartan_type.CartanType_crystallographic):
    r"""
    A class for dual Cartan types.

    The dual of a (crystallographic) Cartan type is a Cartan type with
    the same index set, but all arrows reversed in the Dynkin diagram
    (otherwise said, the Cartan matrix is transposed). It shares a lot
    of properties in common with its dual. In particular, the Weyl
    group is isomorphic to that of the dual as a Coxeter group.

    EXAMPLES:

    For all finite Cartan types, and in particular the simply laced
    ones, the dual Cartan type is given by another preexisting Cartan
    type::

        sage: CartanType(['A',4]).dual()
        ['A', 4]
        sage: CartanType(['B',4]).dual()
        ['C', 4]
        sage: CartanType(['C',4]).dual()
        ['B', 4]
        sage: CartanType(['F',4]).dual()
        ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}

    So to exercise this class we consider some non simply laced affine
    Cartan types and also create explicitly `F_4^*` as a dual cartan
    type::

        sage: from sage.combinat.root_system.type_dual import CartanType as CartanTypeDual
        sage: F4d = CartanTypeDual(CartanType(['F',4])); F4d
        ['F', 4]^*
        sage: G21d = CartanType(['G',2,1]).dual(); G21d
        ['G', 2, 1]^*

    They share many properties with their original Cartan types::

        sage: F4d.is_irreducible()
        True
        sage: F4d.is_crystallographic()
        True
        sage: F4d.is_simply_laced()
        False
        sage: F4d.is_finite()
        True
        sage: G21d.is_finite()
        False
        sage: F4d.is_affine()
        False
        sage: G21d.is_affine()
        True

    TESTS::

        sage: TestSuite(F4d).run(skip=["_test_pickling"])
        sage: TestSuite(G21d).run()

    .. NOTE:: F4d is pickled by construction as F4.dual() hence the above failure.
    """

    def __init__(self, type):
        """
        INPUT:

        - ``type`` -- a Cartan type

        EXAMPLES::

           sage: ct = CartanType(['F',4,1]).dual()
           sage: TestSuite(ct).run()

        TESTS::

            sage: ct1 = CartanType(['B',3,1]).dual()
            sage: ct2 = CartanType(['B',3,1]).dual()
            sage: ct3 = CartanType(['D',4,1]).dual()
            sage: ct1 == ct2
            True
            sage: ct1 == ct3
            False

        Test that the produced Cartan type is in the appropriate
        abstract classes (see :issue:`13724`)::

            sage: from sage.combinat.root_system import cartan_type
            sage: ct = CartanType(['B',3,1]).dual()
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
            False

        By default, the dual of a reducible and finite type is not
        constructed as such::

            sage: ct = CartanType([['B',4],['A',2]]).dual(); ct
            C4xA2

        In order to exercise the dual infrastructure we force the
        construction as a dual::

            sage: from sage.combinat.root_system import type_dual
            sage: ct = type_dual.CartanType(CartanType([['B',4],['A',2]])); ct
            B4xA2^*
            sage: isinstance(ct, type_dual.CartanType)
            True
            sage: TestSuite(ct).run(skip=["_test_pickling"])
            sage: isinstance(ct, cartan_type.CartanType_finite)
            True
            sage: isinstance(ct, cartan_type.CartanType_simple)
            False
            sage: isinstance(ct, cartan_type.CartanType_affine)
            False
            sage: isinstance(ct, cartan_type.CartanType_crystallographic)
            True
            sage: isinstance(ct, cartan_type.CartanType_simply_laced)
            False
        """
        if not type.is_crystallographic():
            raise NotImplementedError("only implemented for crystallographic Cartan types")
        cartan_type.CartanType_decorator.__init__(self, type)
        # TODO: design an appropriate infrastructure to handle this
        # automatically? Maybe using categories and axioms?
        # See also type_relabel.CartanType.__init__
        if type.is_finite():
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
        cartan_type.CartanType_simple]

    def _repr_(self, compact=False):
        """
        EXAMPLES::

           sage: CartanType(['F', 4, 1]).dual()
           ['F', 4, 1]^*

           sage: CartanType(['F', 4, 1]).dual()._repr_(compact = True)
           'F4~*'
        """
        dual_str = self.options.dual_str
        if self.is_affine() and self.options.notation == "Kac":
            if self._type.type() == 'B':
                if compact:
                    return 'A%s^2' % (self.classical().rank()*2-1)
                return "['A', %s, 2]" % (self.classical().rank()*2-1)
            elif self._type.type() == 'BC':
                dual_str = '+'
            elif self._type.type() == 'C':
                if compact:
                    return 'D%s^2' % (self.rank())
                return "['D', %s, 2]" % (self.rank())
            elif self._type.type() == 'F':
                if compact:
                    return 'E6^2'
                return "['E', 6, 2]"
        return self.dual()._repr_(compact)+(dual_str if compact else "^"+dual_str)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(CartanType(['F', 4, 1]).dual())
            F_4^{(1)\vee}
        """
        return self._type._latex_()+"^"+self.options.dual_latex

    def __reduce__(self):
        """
        TESTS::

            sage: CartanType(['F', 4, 1]).dual().__reduce__()
            (*.dual(), (['F', 4, 1],))
        """
        return (attrcall("dual"), (self._type,))

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        EXAMPLES::

            sage: print(CartanType(['F',4,1]).dual()._latex_dynkin_diagram())
            \draw (0 cm,0) -- (2 cm,0);
            {
            \pgftransformxshift{2 cm}
            \draw (0 cm,0) -- (2 cm,0);
            \draw (2 cm, 0.1 cm) -- +(2 cm,0);
            \draw (2 cm, -0.1 cm) -- +(2 cm,0);
            \draw (4.0 cm,0) -- +(2 cm,0);
            \draw[shift={(2.8, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            }
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        return self._type._latex_dynkin_diagram(label, node, node_dist, dual=True)

    def ascii_art(self, label=None, node=None):
        """
        Return an ascii art representation of this Cartan type.

        (by hacking the ascii art representation of the dual Cartan type)

        EXAMPLES::

            sage: print(CartanType(["B", 3, 1]).dual().ascii_art())
                O 0
                |
                |
            O---O=<=O
            1   2   3
            sage: print(CartanType(["C", 4, 1]).dual().ascii_art())
            O=<=O---O---O=>=O
            0   1   2   3   4
            sage: print(CartanType(["G", 2, 1]).dual().ascii_art())
              3
            O=>=O---O
            1   2   0
            sage: print(CartanType(["F", 4, 1]).dual().ascii_art())
            O---O---O=<=O---O
            0   1   2   3   4
            sage: print(CartanType(["BC", 4, 2]).dual().ascii_art())
            O=>=O---O---O=>=O
            0   1   2   3   4
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        res = self._type.ascii_art(label, node)
        # swap, like a computer science freshman!
        # This assumes that the oriented multiple arrows are always ascii arted as =<= or =>=
        res = res.replace("=<=", "=?=")
        res = res.replace("=>=", "=<=")
        res = res.replace("=?=", "=>=")
        return res

    def __eq__(self, other):
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: B41 = CartanType(['B', 4, 1])
            sage: B41dual = CartanType(['B', 4, 1]).dual()
            sage: F41dual = CartanType(['F', 4, 1]).dual()

            sage: F41dual == F41dual
            True
            sage: F41dual == B41dual
            False
            sage: B41dual == B41
            False
        """
        if not isinstance(other, CartanType):
            return False
        return self._type == other._type

    def __ne__(self, other):
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: B41 = CartanType(['B', 4, 1])
            sage: B41dual = CartanType(['B', 4, 1]).dual()
            sage: F41dual = CartanType(['F', 4, 1]).dual()

            sage: F41dual != F41dual
            False
            sage: F41dual != B41dual
            True
            sage: B41dual != B41
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Compute the hash of ``self``.

        EXAMPLES::

            sage: B41 = CartanType(['B', 4, 1])
            sage: B41dual = CartanType(['B', 4, 1]).dual()
            sage: h = hash(B41dual)
        """
        return hash(self._type)

    def dual(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4, 1]).dual()
           sage: ct.dual()
           ['F', 4, 1]
        """
        return self._type

    def dynkin_diagram(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['F', 4, 1]).dual()
            sage: ct.dynkin_diagram()                                                   # needs sage.graphs
            O---O---O=<=O---O
            0   1   2   3   4
            F4~*
        """
        return self._type.dynkin_diagram().dual()

###########################################################################


class AmbientSpace(ambient_space.AmbientSpace):
    """
    Ambient space for a dual finite Cartan type.

    It is constructed in the canonical way from the ambient space of
    the original Cartan type by switching the roles of simple roots,
    fundamental weights, etc.

    .. NOTE::

        Recall that, for any finite Cartan type, and in particular the
        a simply laced one, the dual Cartan type is constructed as
        another preexisting Cartan type. Furthermore the ambient space
        for an affine type is constructed from the ambient space for
        its classical type. Thus this code is not actually currently
        used.

        It is kept for cross-checking and for reference in case it
        could become useful, e.g., for dual of general Kac-Moody
        types.

        For the doctests, we need to explicitly create a dual type.
        Subsequently, since reconstruction of the dual of type `F_4`
        is the relabelled Cartan type, pickling fails on the
        ``TestSuite`` run.

    EXAMPLES::

        sage: ct = sage.combinat.root_system.type_dual.CartanType(CartanType(['F',4]))
        sage: L = ct.root_system().ambient_space(); L
        Ambient space of the Root system of type ['F', 4]^*
        sage: TestSuite(L).run(skip=["_test_elements","_test_pickling"])                # needs sage.graphs
    """

    @lazy_attribute
    def _dual_space(self):
        """
        The dual of this ambient space.

        EXAMPLES::

            sage: ct = sage.combinat.root_system.type_dual.CartanType(CartanType(['F',4]))
            sage: L = ct.root_system().ambient_space(); L
            Ambient space of the Root system of type ['F', 4]^*
            sage: L._dual_space
            Ambient space of the Root system of type ['F', 4]

        The basic data for this space is fetched from the dual space::

            sage: L._dual_space.simple_root(1)
            (0, 1, -1, 0)
            sage: L.simple_root(1)
            (0, 1, -1, 0)
        """
        K = self.base_ring()
        return self.cartan_type().dual().root_system().ambient_space(K)
        #return self.root_system.dual.ambient_space()

    def dimension(self):
        """
        Return the dimension of this ambient space.

        .. SEEALSO:: :meth:`sage.combinat.root_system.ambient_space.AmbientSpace.dimension`

        EXAMPLES::

            sage: ct = sage.combinat.root_system.type_dual.CartanType(CartanType(['F',4]))
            sage: L = ct.root_system().ambient_space()
            sage: L.dimension()
            4
        """
        # Can't yet use _dual_space for the base ring (and the Cartan type?) is not yet initialized
        return self.root_system.dual.ambient_space().dimension()

    @cached_method
    def simple_root(self, i):
        """
        Return the ``i``-th simple root.

        It is constructed by looking up the corresponding simple
        coroot in the ambient space for the dual Cartan type.

        EXAMPLES::

            sage: ct = sage.combinat.root_system.type_dual.CartanType(CartanType(['F',4]))
            sage: ct.root_system().ambient_space().simple_root(1)
            (0, 1, -1, 0)

            sage: ct.root_system().ambient_space().simple_roots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 2), 4: (1, -1, -1, -1)}

            sage: ct.dual().root_system().ambient_space().simple_coroots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 2), 4: (1, -1, -1, -1)}

        Note that this ambient space is isomorphic, but not equal, to
        that obtained by constructing `F_4` dual by relabelling::

            sage: ct = CartanType(['F',4]).dual(); ct
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
            sage: ct.root_system().ambient_space().simple_roots()
            Finite family {1: (1/2, -1/2, -1/2, -1/2), 2: (0, 0, 0, 1), 3: (0, 0, 1, -1), 4: (0, 1, -1, 0)}
        """
        dual_coroot = self._dual_space.simple_coroot(i)
        return self.sum_of_terms(dual_coroot)

    @cached_method
    def fundamental_weights(self):
        """
        Return the fundamental weights.

        They are computed from the simple roots by inverting the
        Cartan matrix. This is acceptable since this is only about
        ambient spaces for finite Cartan types. Also, we do not have
        to worry about the usual `GL_n` vs `SL_n` catch because type
        `A` is self dual.

        An alternative would have been to start from the fundamental
        coweights in the dual ambient space, but those are not yet
        implemented.

        EXAMPLES::

            sage: ct = sage.combinat.root_system.type_dual.CartanType(CartanType(['F',4]))
            sage: L = ct.root_system().ambient_space()
            sage: L.fundamental_weights()                                               # needs sage.graphs
            Finite family {1: (1, 1, 0, 0), 2: (2, 1, 1, 0), 3: (3, 1, 1, 1), 4: (2, 0, 0, 0)}

        Note that this ambient space is isomorphic, but not equal, to
        that obtained by constructing `F_4` dual by relabelling::

            sage: ct = CartanType(['F',4]).dual(); ct
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
            sage: ct.root_system().ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (3/2, 1/2, 1/2, 1/2), 3: (2, 1, 1, 0), 4: (1, 1, 0, 0)}
        """
        return self.fundamental_weights_from_simple_roots()

    @lazy_attribute
    def _plot_projection(self):
        """
        Return the default plot projection for ``self``.

        If an ambient space uses barycentric projection, then so does
        its dual.

        .. SEEALSO::

            - :meth:`sage.combinat.root_system.root_lattice_realizations.RootLatticeRealizations.ParentMethods._plot_projection`

        EXAMPLES::

            sage: ct = sage.combinat.root_system.type_dual.CartanType(CartanType(['G',2]))
            sage: L = ct.root_system().ambient_space()
            sage: L._plot_projection == L._plot_projection_barycentric
            True

            sage: L = RootSystem(['G',2]).coambient_space()
            sage: L._plot_projection == L._plot_projection_barycentric
            True
        """
        dual_space = self.cartan_type().dual().root_system().ambient_space(self.base_ring())
        if dual_space._plot_projection == dual_space._plot_projection_barycentric:
            return self._plot_projection_barycentric
        else:
            RootLatticeRealizations.ParentMethods.__dict__["_plot_projection"]


class CartanType_finite(CartanType, cartan_type.CartanType_finite):
    AmbientSpace = AmbientSpace

###########################################################################


class CartanType_affine(CartanType, cartan_type.CartanType_affine):
    def classical(self):
        """
        Return the classical Cartan type associated with ``self`` (which should
        be affine).

        EXAMPLES::

            sage: CartanType(['A',3,1]).dual().classical()
            ['A', 3]
            sage: CartanType(['B',3,1]).dual().classical()
            ['C', 3]
            sage: CartanType(['F',4,1]).dual().classical()
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
            sage: CartanType(['BC',4,2]).dual().classical()
            ['B', 4]
        """
        return self.dual().classical().dual()

    def basic_untwisted(self):
        r"""
        Return the basic untwisted Cartan type associated with this affine
        Cartan type.

        Given an affine type `X_n^{(r)}`, the basic untwisted type is `X_n`.
        In other words, it is the classical Cartan type that is twisted to
        obtain ``self``.

        EXAMPLES::

            sage: CartanType(['A', 7, 2]).basic_untwisted()
            ['A', 7]
            sage: CartanType(['E', 6, 2]).basic_untwisted()
            ['E', 6]
            sage: CartanType(['D', 4, 3]).basic_untwisted()
            ['D', 4]
        """
        from . import cartan_type
        if self.dual().type() == 'B':
            return cartan_type.CartanType(['A', self.classical().rank()*2-1])
        elif self.dual().type() == 'BC':
            return cartan_type.CartanType(['A', self.classical().rank()*2])
        elif self.dual().type() == 'C':
            return cartan_type.CartanType(['D', self.classical().rank()+1])
        elif self.dual().type() == 'F':
            return cartan_type.CartanType(['E', 6])
        elif self.dual().type() == 'G':
            return cartan_type.CartanType(['D', 4])

    def special_node(self):
        """
        Implement :meth:`CartanType_affine.special_node`.

        The special node of the dual of an affine type `T` is the
        special node of `T`.

        EXAMPLES::

            sage: CartanType(['A',3,1]).dual().special_node()
            0
            sage: CartanType(['B',3,1]).dual().special_node()
            0
            sage: CartanType(['F',4,1]).dual().special_node()
            0
            sage: CartanType(['BC',4,2]).dual().special_node()
            0
        """
        return self.dual().special_node()

    def _repr_(self, compact=False):
        """
        EXAMPLES::

           sage: CartanType(['F', 4, 1]).dual()
           ['F', 4, 1]^*

           sage: CartanType(['F', 4, 1]).dual()._repr_(compact = True)
           'F4~*'
        """
        if self.options.notation == "Kac":
            if self._type.type() == 'B':
                if compact:
                    return 'A%s^2' % (self.classical().rank()*2-1)
                return "['A', %s, 2]" % (self.classical().rank()*2-1)
            elif self._type.type() == 'BC':
                pass
            elif self._type.type() == 'C':
                if compact:
                    return 'D%s^2' % (self.rank())
                return "['D', %s, 2]" % (self.rank())
            elif self._type.type() == 'F':
                if compact:
                    return 'E6^2'
                return "['E', 6, 2]"
        return CartanType._repr_(self, compact)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['B',4,1]).dual())
            B_{4}^{(1)\vee}
            sage: latex(CartanType(['BC',4,2]).dual())
            BC_{4}^{(2)\vee}
            sage: latex(CartanType(['G',2,1]).dual())
            G_2^{(1)\vee}

            sage: CartanType.options['notation'] = 'Kac'
            sage: latex(CartanType(['A',7,2]))
            A_{7}^{(2)}
            sage: latex(CartanType(['B',4,1]).dual())
            A_{7}^{(2)}
            sage: latex(CartanType(['A',8,2]))
            A_{8}^{(2)}
            sage: latex(CartanType(['A',8,2]).dual())
            A_{8}^{(2)\dagger}
            sage: latex(CartanType(['E',6,2]))
            E_6^{(2)}
            sage: latex(CartanType(['D',5,2]))
            D_{5}^{(2)}
            sage: CartanType.options._reset()
        """
        if self.options('notation') == "Kac":
            if self._type.type() == 'B':
                return "A_{%s}^{(2)}" % (self.classical().rank()*2-1)
            elif self._type.type() == 'BC':
                return "A_{%s}^{(2)\\dagger}" % (2*self.classical().rank())
            elif self._type.type() == 'C':
                return "D_{%s}^{(2)}" % (self.rank)()
            elif self._type.type() == 'F':
                return "E_6^{(2)}"
        result = self._type._latex_()
        import re
        if re.match(r".*\^{\(\d\)}$", result):
            return "%s%s}" % (result[:-1], self.options('dual_latex'))
        else:
            return "{%s}^%s" % (result, self.options('dual_latex'))

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['A', 6, 2]).dual()._default_folded_cartan_type()
            ['BC', 3, 2]^* as a folding of ['A', 5, 1]
            sage: CartanType(['A', 5, 2])._default_folded_cartan_type()
            ['B', 3, 1]^* as a folding of ['D', 4, 1]
            sage: CartanType(['D', 4, 2])._default_folded_cartan_type()
            ['C', 3, 1]^* as a folding of ['A', 5, 1]
            sage: CartanType(['E', 6, 2])._default_folded_cartan_type()
            ['F', 4, 1]^* as a folding of ['E', 6, 1]
            sage: CartanType(['G', 2, 1]).dual()._default_folded_cartan_type()
            ['G', 2, 1]^* as a folding of ['D', 4, 1]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        letter = self._type.type()
        if letter == 'BC':  # A_{2n}^{(2)\dagger}
            n = self._type.classical().rank()
            return CartanTypeFolded(self, ['A', 2*n - 1, 1],
                [[0]] + [[i, 2*n-i] for i in range(1, n)] + [[n]])
        if letter == 'B':  # A_{2n-1}^{(2)}
            n = self._type.classical().rank()
            return CartanTypeFolded(self, ['D', n + 1, 1],
                [[i] for i in range(n)] + [[n, n+1]])
        if letter == 'C':  # D_{n+1}^{(2)}
            n = self._type.classical().rank()
            return CartanTypeFolded(self, ['A', 2*n-1, 1],
                [[0]] + [[i, 2*n-i] for i in range(1, n)] + [[n]])
        if letter == 'F':  # E_6^{(2)}
            return CartanTypeFolded(self, ['E', 6, 1], [[0], [2], [4], [3, 5], [1, 6]])
        if letter == 'G':  # D_4^{(3)}
            return CartanTypeFolded(self, ['D', 4, 1], [[0], [1, 3, 4], [2]])
        return super()._default_folded_cartan_type()
