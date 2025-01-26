# sage.doctest: needs sage.graphs sage.groups
r"""
Fundamental group of an extended affine Weyl group

AUTHORS:

- Mark Shimozono (2013) initial version
"""
# ****************************************************************************
#       Copyright (C) 2013 Mark Shimozono <mshimo at math.vt.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.combinat.root_system.cartan_type import CartanType
from sage.categories.groups import Groups
from sage.misc.cachefunc import cached_method
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.family import Family
from sage.combinat.root_system.root_system import RootSystem
from sage.categories.category import Category
from sage.categories.enumerated_sets import EnumeratedSets
from sage.rings.integer_ring import ZZ
from sage.sets.family import LazyFamily


def FundamentalGroupOfExtendedAffineWeylGroup(cartan_type, prefix='pi',
                                              general_linear=None):
    r"""
    Factory for the fundamental group of an extended affine Weyl group.

    INPUT:

    - ``cartan_type`` -- a Cartan type that is either affine or finite, with the latter being a
      shorthand for the untwisted affinization
    - ``prefix`` -- (default: ``'pi'``) string that labels the elements of the group
    - ``general_linear`` -- (default: ``None``, meaning ``False``) In untwisted
      type `A`, if ``True``, use the universal central extension

    .. RUBRIC:: Fundamental group

    Associated to each affine Cartan type `\tilde{X}` is an extended affine Weyl group `E`.
    Its subgroup of length-zero elements is called the fundamental group `F`.
    The group `F` can be identified with a subgroup of the group of automorphisms of the
    affine Dynkin diagram. As such, every element of `F` can be viewed as a permutation of the
    set `I` of affine Dynkin nodes.

    Let `0 \in I` be the distinguished affine node; it is the one whose removal produces the
    associated finite Cartan type (call it `X`). A node `i \in I` is called *special*
    if some automorphism of the affine Dynkin diagram, sends `0` to `i`.
    The node `0` is always special due to the identity automorphism.
    There is a bijection of the set of special nodes with the fundamental group. We denote the
    image of `i` by `\pi_i`. The structure of `F` is determined as follows.

    - `\tilde{X}` is untwisted -- `F` is isomorphic to `P^\vee/Q^\vee` where `P^\vee` and `Q^\vee` are the
      coweight and coroot lattices of type `X`. The group `P^\vee/Q^\vee` consists of the cosets `\omega_i^\vee + Q^\vee`
      for special nodes `i`, where `\omega_0^\vee = 0` by convention. In this case the special nodes `i`
      are the *cominuscule* nodes, the ones such that `\omega_i^\vee(\alpha_j)` is `0` or `1` for all `j\in I_0 = I \setminus \{0\}`.
      For `i` special, addition by `\omega_i^\vee+Q^\vee` permutes `P^\vee/Q^\vee` and therefore permutes the set of special nodes.
      This permutation extends uniquely to an automorphism of the affine Dynkin diagram.
    - `\tilde{X}` is dual untwisted -- (that is, the dual of `\tilde{X}` is untwisted) `F` is isomorphic to `P/Q`
      where `P` and `Q` are the weight and root lattices of type `X`. The group `P/Q` consists of the cosets
      `\omega_i + Q` for special nodes `i`, where `\omega_0 = 0` by convention. In this case the special nodes `i`
      are the *minuscule* nodes, the ones such that `\alpha_j^\vee(\omega_i)` is `0` or `1` for all `j \in I_0`.
      For `i` special, addition by `\omega_i+Q` permutes `P/Q` and therefore permutes the set of special nodes.
      This permutation extends uniquely to an automorphism of the affine Dynkin diagram.
    - `\tilde{X}` is mixed -- (that is, not of the above two types) `F` is the trivial group.

    EXAMPLES::

        sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
        sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1]); F
        Fundamental group of type ['A', 3, 1]
        sage: F.cartan_type().dynkin_diagram()
        0
        O-------+
        |       |
        |       |
        O---O---O
        1   2   3
        A3~
        sage: F.special_nodes()
        (0, 1, 2, 3)
        sage: F(1)^2
        pi[2]
        sage: F(1)*F(2)
        pi[3]
        sage: F(3)^(-1)
        pi[1]

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("B3"); F
        Fundamental group of type ['B', 3, 1]
        sage: F.cartan_type().dynkin_diagram()
            O 0
            |
            |
        O---O=>=O
        1   2   3
        B3~
        sage: F.special_nodes()
        (0, 1)

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("C2"); F
        Fundamental group of type ['C', 2, 1]
        sage: F.cartan_type().dynkin_diagram()
        O=>=O=<=O
        0   1   2
        C2~
        sage: F.special_nodes()
        (0, 2)

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("D4"); F
        Fundamental group of type ['D', 4, 1]
        sage: F.cartan_type().dynkin_diagram()
            O 4
            |
            |
        O---O---O
        1   |2  3
            |
            O 0
        D4~
        sage: F.special_nodes()
        (0, 1, 3, 4)
        sage: (F(4), F(4)^2)
        (pi[4], pi[0])

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("D5"); F
        Fundamental group of type ['D', 5, 1]
        sage: F.cartan_type().dynkin_diagram()
          0 O   O 5
            |   |
            |   |
        O---O---O---O
        1   2   3   4
        D5~
        sage: F.special_nodes()
        (0, 1, 4, 5)
        sage: (F(5), F(5)^2, F(5)^3, F(5)^4)
        (pi[5], pi[1], pi[4], pi[0])
        sage: F = FundamentalGroupOfExtendedAffineWeylGroup("E6"); F
        Fundamental group of type ['E', 6, 1]
        sage: F.cartan_type().dynkin_diagram()
                O 0
                |
                |
                O 2
                |
                |
        O---O---O---O---O
        1   3   4   5   6
        E6~
        sage: F.special_nodes()
        (0, 1, 6)
        sage: F(1)^2
        pi[6]

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['D',4,2]); F
        Fundamental group of type ['C', 3, 1]^*
        sage: F.cartan_type().dynkin_diagram()
        O=<=O---O=>=O
        0   1   2   3
        C3~*
        sage: F.special_nodes()
        (0, 3)

    We also implement a fundamental group for `GL_n`. It is defined to be the group of integers, which is the
    covering group of the fundamental group Z/nZ for affine `SL_n`::

        sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True); F
        Fundamental group of GL(3)
        sage: x = F.an_element(); x
        pi[5]
        sage: x*x
        pi[10]
        sage: x.inverse()
        pi[-5]
        sage: wt = F.cartan_type().classical().root_system().ambient_space().an_element(); wt
        (2, 2, 3)
        sage: x.act_on_classical_ambient(wt)
        (2, 3, 2)
        sage: w = WeylGroup(F.cartan_type(),prefix='s').an_element(); w
        s0*s1*s2
        sage: x.act_on_affine_weyl(w)
        s2*s0*s1
    """
    cartan_type = CartanType(cartan_type)
    if cartan_type.is_finite():
        cartan_type = cartan_type.affine()
    if not cartan_type.is_affine():
        raise NotImplementedError("Cartan type is not affine")
    if general_linear is True:
        if cartan_type.is_untwisted_affine() and cartan_type.type() == "A":
            return FundamentalGroupGL(cartan_type, prefix)
        else:
            raise ValueError("General Linear Fundamental group is untwisted type A")
    return FundamentalGroupOfExtendedAffineWeylGroup_Class(cartan_type, prefix,
                                                           finite=True)


class FundamentalGroupElement(MultiplicativeGroupElement):
    def __init__(self, parent, x):
        r"""
        This should not be called directly.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: x = FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1], prefix='f').an_element()
            sage: TestSuite(x).run()
        """
        if x not in parent.special_nodes():
            raise ValueError("%s is not a special node" % x)
        self._value = x
        MultiplicativeGroupElement.__init__(self, parent)

    def value(self):
        r"""
        Return the special node which indexes the special automorphism ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1], prefix='f')
            sage: F.special_nodes()
            (0, 1, 2, 3, 4)
            sage: x = F(4); x
            f[4]
            sage: x.value()
            4
        """
        return self._value

    def _repr_(self):
        r"""
        Return a string representing ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1], prefix='f')
            sage: F(2)^3 # indirect doctest
            f[1]
        """
        return self.parent()._prefix + "[" + repr(self.value()) + "]"

    def __invert__(self):
        r"""
        Return the inverse element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: F(1).inverse()   # indirect doctest
            pi[3]
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['E',6,1], prefix='f')
            sage: F(1).inverse()
            f[6]
        """
        par = self.parent()
        return self.__class__(par, par.dual_node(self.value()))

    def _richcmp_(self, x, op):
        r"""
        Compare ``self`` with `x`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: x = F(0); y = F(2)
            sage: y > x
            True
            sage: y == y
            True
            sage: y != y
            False
            sage: x <= y
            True
        """
        return richcmp(self.value(), x.value(), op)

    def act_on_affine_weyl(self, w):
        r"""
        Act by ``self`` on the element `w` of the affine Weyl group.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: W = WeylGroup(F.cartan_type(),prefix='s')
            sage: w = W.from_reduced_word([2,3,0])
            sage: F(1).act_on_affine_weyl(w).reduced_word()
            [3, 0, 1]
        """
        par = self.parent()
        if self == par.one():
            return w
        action = par.action(self.value())
        return w.parent().from_reduced_word([action(j) for j in w.reduced_word()])

    def act_on_affine_lattice(self, wt):
        r"""
        Act by ``self`` on the element ``wt`` of an affine root/weight lattice realization.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: wt = RootSystem(F.cartan_type()).weight_lattice().an_element(); wt
            2*Lambda[0] + 2*Lambda[1] + 3*Lambda[2]
            sage: F(3).act_on_affine_lattice(wt)
            2*Lambda[0] + 3*Lambda[1] + 2*Lambda[3]

        .. WARNING::

            Doesn't work on ambient spaces.
        """
        return wt.map_support(self.parent().action(self.value()))

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        TESTS:

        Check that :issue:`36121` is fixed::

            sage: W = ExtendedAffineWeylGroup('A4')
            sage: fun = W.fundamental_group().an_element()
            sage: hash(fun) == hash(fun.value())
            True
        """
        return hash(self.value())


class FundamentalGroupOfExtendedAffineWeylGroup_Class(UniqueRepresentation,
                                                      Parent):
    r"""
    The group of length zero elements in the extended affine Weyl group.
    """
    Element = FundamentalGroupElement

    def __init__(self, cartan_type, prefix, finite=True):
        r"""

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: F in Groups().Commutative().Finite()
            True
            sage: TestSuite(F).run()
        """
        def leading_support(beta):
            r"""
            Given a dictionary with one key, return this key
            """
            supp = beta.support()
            assert len(supp) == 1
            return supp[0]

        self._cartan_type = cartan_type
        self._prefix = prefix
        special_node = cartan_type.special_node()
        self._special_nodes = cartan_type.special_nodes()

        # initialize dictionaries with the entries for the
        # distinguished special node

        # dictionary of inverse elements
        inverse_dict = {}
        inverse_dict[special_node] = special_node
        # dictionary for the action of special automorphisms by
        # permutations of the affine Dynkin nodes
        auto_dict = {}
        for i in cartan_type.index_set():
            auto_dict[special_node,i] = i
        # dictionary for the finite Weyl component of the special automorphisms
        reduced_words_dict = {}
        reduced_words_dict[0] = tuple()

        if cartan_type.dual().is_untwisted_affine():
            # this combines the computations for an untwisted affine
            # type and its affine dual
            cartan_type = cartan_type.dual()
        if cartan_type.is_untwisted_affine():
            cartan_type_classical = cartan_type.classical()
            I = list(cartan_type_classical.index_set())
            Q = RootSystem(cartan_type_classical).root_lattice()
            alpha = Q.simple_roots()
            omega = RootSystem(cartan_type_classical).weight_lattice().fundamental_weights()
            W = Q.weyl_group(prefix='s')
            for i in self._special_nodes:
                if i == special_node:
                    continue
                antidominant_weight, reduced_word = omega[i].to_dominant_chamber(reduced_word=True, positive=False)
                reduced_words_dict[i] = tuple(reduced_word)
                w0i = W.from_reduced_word(reduced_word)
                idual = leading_support(-antidominant_weight)
                inverse_dict[i] = idual
                auto_dict[i,special_node] = i
                for j in I:
                    if j == idual:
                        auto_dict[i, j] = special_node
                    else:
                        auto_dict[i, j] = leading_support(w0i.action(alpha[j]))

        self._action = Family(self._special_nodes, lambda i: Family(cartan_type.index_set(), lambda j: auto_dict[i, j]))
        self._dual_node = Family(self._special_nodes, inverse_dict.__getitem__)
        self._reduced_words = Family(self._special_nodes, reduced_words_dict.__getitem__)

        if finite:
            cat = Category.join((Groups().Commutative().Finite(),
                                 EnumeratedSets()))
        else:
            cat = Groups().Commutative().Infinite()
        Parent.__init__(self, category=cat)

    @cached_method
    def one(self):
        r"""
        Return the identity element of the fundamental group.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: F.one()
            pi[0]
        """
        return self(self.cartan_type().special_node())

    def product(self, x, y):
        r"""
        Return the product of `x` and `y`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: F.special_nodes()
            (0, 1, 2, 3)
            sage: F(2)*F(3)
            pi[1]
            sage: F(1)*F(3)^(-1)
            pi[2]
        """
        return self(self.action(x.value())(y.value()))

    def cartan_type(self):
        r"""
        The Cartan type of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1]).cartan_type()
            ['A', 3, 1]
        """
        return self._cartan_type

    def _repr_(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1]) # indirect doctest
            Fundamental group of type ['A', 3, 1]
        """
        return "Fundamental group of type %s" % self.cartan_type()

    def special_nodes(self):
        r"""
        Return the special nodes of ``self``.

        See :meth:`sage.combinat.root_system.cartan_type.special_nodes()`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['D',4,1]).special_nodes()
            (0, 1, 3, 4)
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1]).special_nodes()
            (0, 1, 2)
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['C',3,1]).special_nodes()
            (0, 3)
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['D',4,2]).special_nodes()
            (0, 3)
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).special_nodes()
            Integer Ring
        """
        return self._special_nodes

    def group_generators(self):
        r"""
        Return a tuple of generators of the fundamental group.

        .. WARNING::

            This returns the entire group, a necessary behavior because it
            is used in :meth:`__iter__`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['E',6,1],prefix='f').group_generators()
            Finite family {0: f[0], 1: f[1], 6: f[6]}
        """
        return Family(self.special_nodes(), self)

    def __iter__(self):
        r"""
        Return the iterator for ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['E',6,1],prefix='f')
            sage: [x for x in F] # indirect doctest
            [f[0], f[1], f[6]]
        """
        return iter(self.group_generators())

    @cached_method
    def an_element(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1],prefix='f').an_element()
            f[4]
        """
        return self.last()

    @cached_method
    def index_set(self):
        r"""
        The node set of the affine Cartan type of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1]).index_set()
            (0, 1, 2)
        """
        return self.cartan_type().index_set()

    def action(self, i):
        r"""
        Return a function which permutes the affine Dynkin node set by the `i`-th special automorphism.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1])
            sage: [[(i, j, F.action(i)(j)) for j in F.index_set()] for i in F.special_nodes()]
            [[(0, 0, 0), (0, 1, 1), (0, 2, 2)], [(1, 0, 1), (1, 1, 2), (1, 2, 0)], [(2, 0, 2), (2, 1, 0), (2, 2, 1)]]
            sage: G = FundamentalGroupOfExtendedAffineWeylGroup(['D',4,1])
            sage: [[(i, j, G.action(i)(j)) for j in G.index_set()] for i in G.special_nodes()]
            [[(0, 0, 0), (0, 1, 1), (0, 2, 2), (0, 3, 3), (0, 4, 4)], [(1, 0, 1), (1, 1, 0), (1, 2, 2), (1, 3, 4), (1, 4, 3)], [(3, 0, 3), (3, 1, 4), (3, 2, 2), (3, 3, 0), (3, 4, 1)], [(4, 0, 4), (4, 1, 3), (4, 2, 2), (4, 3, 1), (4, 4, 0)]]
        """
        return lambda j: self._action[i][j]

    def dual_node(self, i):
        r"""
        Return the node that indexes the inverse of the `i`-th element.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',4,1])
            sage: [(i, F.dual_node(i)) for i in F.special_nodes()]
            [(0, 0), (1, 4), (2, 3), (3, 2), (4, 1)]
            sage: G = FundamentalGroupOfExtendedAffineWeylGroup(['E',6,1])
            sage: [(i, G.dual_node(i)) for i in G.special_nodes()]
            [(0, 0), (1, 6), (6, 1)]
            sage: H = FundamentalGroupOfExtendedAffineWeylGroup(['D',5,1])
            sage: [(i, H.dual_node(i)) for i in H.special_nodes()]
            [(0, 0), (1, 1), (4, 5), (5, 4)]
        """
        return self._dual_node[i]

    def reduced_word(self, i):
        r"""
        Return a reduced word for the finite Weyl group element associated with the `i`-th special automorphism.

        More precisely, for each special node `i`, ``self.reduced_word(i)`` is a reduced word for
        the element `v` in the finite Weyl group such that in the extended affine Weyl group,
        the `i`-th special automorphism is equal to `t v` where `t` is a translation element.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',3,1])
            sage: [(i, F.reduced_word(i)) for i in F.special_nodes()]
            [(0, ()), (1, (1, 2, 3)), (2, (2, 1, 3, 2)), (3, (3, 2, 1))]
        """
        return self._reduced_words[i]


class FundamentalGroupGLElement(FundamentalGroupElement):
    def act_on_classical_ambient(self, wt):
        r"""
        Act by ``self`` on the classical ambient weight vector ``wt``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: f = F.an_element(); f
            pi[5]
            sage: la = F.cartan_type().classical().root_system().ambient_space().an_element(); la
            (2, 2, 3)
            sage: f.act_on_classical_ambient(la)
            (2, 3, 2)
        """
        return wt.map_support(self.parent().action(self.value()))


class FundamentalGroupGL(FundamentalGroupOfExtendedAffineWeylGroup_Class):
    r"""
    Fundamental group of `GL_n`. It is just the integers with extra privileges.
    """
    Element = FundamentalGroupGLElement

    def __init__(self, cartan_type, prefix='pi'):
        r"""

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F in Groups().Commutative().Infinite()
            True
            sage: TestSuite(F).run()
        """
        FundamentalGroupOfExtendedAffineWeylGroup_Class.__init__(self, cartan_type, prefix, finite=False)
        self._special_nodes = ZZ
        self._n = cartan_type.n + 1

    @cached_method
    def one(self):
        r"""
        Return the identity element of the fundamental group.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).one()
            pi[0]
        """
        return self(ZZ.zero())

    def product(self, x, y):
        r"""
        Return the product of `x` and `y`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.special_nodes()
            Integer Ring
            sage: F(2)*F(3)
            pi[5]
            sage: F(1)*F(3)^(-1)
            pi[-2]
        """
        return self(x.value() + y.value())

    def _repr_(self):
        r"""
        Return a string representing the fundamental group.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True) # indirect doctest
            Fundamental group of GL(3)
        """
        return "Fundamental group of GL(%s)" % self._n

    def family(self):
        r"""
        The family associated with the set of special nodes.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: fam = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).family() # indirect doctest
            sage: fam
            Lazy family (<lambda>(i))_{i in Integer Ring}
            sage: fam[-3]
            -3
        """
        return LazyFamily(ZZ, lambda i: i)

    @cached_method
    def an_element(self):
        r"""
        An element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).an_element()
            pi[5]
        """
        return self(ZZ(5))

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).some_elements()
            [pi[-2], pi[2], pi[5]]
        """
        return [self(ZZ(i)) for i in [-2, 2, 5]]

    def group_generators(self):
        r"""
        Return group generators for ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True).group_generators()
            (pi[1],)
        """
        return (self(ZZ.one()),)

    def action(self, i):
        r"""
        The action of the `i`-th automorphism on the affine Dynkin node set.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.action(4)(2)
            0
            sage: F.action(-4)(2)
            1
        """
        return lambda j: (i + j) % self._n

    def dual_node(self, i):
        r"""
        The node whose special automorphism is inverse to that of `i`.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.dual_node(2)
            -2
        """
        return -i

    @cached_method
    def reduced_word(self, i):
        r"""
        A reduced word for the finite permutation part of the
        special automorphism indexed by `i`.

        More precisely, return a reduced word for the finite Weyl group element `u`
        where `i`-th automorphism (expressed in the extended affine Weyl group)
        has the form `t u` where `t` is a translation element.

        EXAMPLES::

            sage: from sage.combinat.root_system.fundamental_group import FundamentalGroupOfExtendedAffineWeylGroup
            sage: F = FundamentalGroupOfExtendedAffineWeylGroup(['A',2,1], general_linear=True)
            sage: F.reduced_word(10)
            (1, 2)
        """
        i = i % self._n
        if i == 0:
            return tuple()
        om = self.cartan_type().classical().root_system().weight_lattice().fundamental_weight(i)
        return tuple((-om).reduced_word())
