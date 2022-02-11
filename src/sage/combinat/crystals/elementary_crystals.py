r"""
Elementary Crystals

Let `\lambda` be a weight. The crystals `T_{\lambda}`, `R_{\lambda}`, `B_i`,
and `C` are important objects in the tensor category of crystals.
For example, the crystal `T_0` is the neutral object in this category; i.e.,
`T_0 \otimes B \cong B \otimes T_0 \cong B` for any crystal `B`.  We list
some other properties of these crystals:

- The crystal `T_{\lambda} \otimes B(\infty)` is the crystal of the Verma
  module with highest weight `\lambda`, where `\lambda` is a dominant integral
  weight.

- Let `u_{\infty}` be the highest weight vector of `B(\infty)` and `\lambda`
  be a dominant integral weight. There is an embedding of crystals `B(\lambda)
  \longrightarrow T_{\lambda} \otimes B(\infty)` sending `u_{\lambda} \mapsto
  t_{\lambda} \otimes u_{\infty}` which is not strict, but the embedding
  `B(\lambda) \longrightarrow C \otimes T_{\lambda} \otimes B(\infty)` by
  `u_{\lambda} \mapsto c \otimes t_{\lambda} \otimes u_{\infty}` is a strict
  embedding.

- For any dominant integral weight `\lambda`, there is a surjective crystal
  morphism `\Psi_{\lambda} \colon R_{\lambda} \otimes B(\infty) \longrightarrow
  B(\lambda)`.  More precisely, if `B = \{r_{\lambda} \otimes b \in R_{\lambda}
  \otimes B(\infty) : \Psi_{\lambda}(r_{\lambda} \otimes b) \neq 0 \}`, then
  `B \cong B(\lambda)` as crystals.

- For all Cartan types and all weights `\lambda`, we have `R_{\lambda} \cong C
  \otimes T_{\lambda}` as crystals.

- For each `i`, there is a strict crystal morphism `\Psi_i \colon B(\infty)
  \longrightarrow B_i \otimes B(\infty)` defined by `u_{\infty} \mapsto
  b_i(0) \otimes u_{\infty}`, where `u_\infty` is the highest weight vector
  of `B(\infty)`.

For more information on `B(\infty)`, see
:class:`~sage.combinat.crystals.infinity_crystals.InfinityCrystalOfTableaux`.

.. NOTE::

    As with
    :class:`~sage.combinat.crystals.tensor_product.TensorProductOfCrystals`,
    we are using the opposite of Kashiwara's convention.

AUTHORS:

- Ben Salisbury: Initial version

REFERENCES:

- [Ka1993]_

- [NZ1997]_
"""

# ****************************************************************************
#  Copyright (C) 2013 Ben Salisbury <benjamin_salisbury at brown dot edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.categories.crystals import Crystals
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.supercrystals import SuperCrystals
from sage.categories.regular_supercrystals import RegularSuperCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


class AbstractSingleCrystalElement(Element):
    r"""
    Abstract base class for elements in crystals with a single element.
    """
    def __lt__(self, other):
        r"""
        EXAMPLES::

            sage: La = RootSystem("D4").ambient_space().fundamental_weights()
            sage: T = crystals.elementary.T("D4",La[3]+La[4])
            sage: t = T.highest_weight_vector()
            sage: t < t.e(1)
            False
            sage: t < t
            False
        """
        return False

    def __hash__(self):
        r"""
        TESTS::

            sage: C = crystals.elementary.Component("D7")
            sage: c = C.highest_weight_vector()
            sage: hash(c) # random
            879
        """
        return hash(self.parent())

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: La = RootSystem("A2").weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T("A2",La[1])
            sage: U = crystals.elementary.T("A2",La[2])
            sage: la = RootSystem("B2").weight_lattice().fundamental_weights()
            sage: V = crystals.elementary.T("B2",la[1])
            sage: t = T.highest_weight_vector()
            sage: u = U.highest_weight_vector()
            sage: v = V.highest_weight_vector()
            sage: [t == t, u == u, v == v]
            [True, True, True]
            sage: [t == u, u == v, t == v]
            [False, False, False]

            sage: C = crystals.elementary.Component("D7")
            sage: c = C.highest_weight_vector()
            sage: c == c
            True
            sage: c == c.f(7)
            False
        """
        if isinstance(other, AbstractSingleCrystalElement):
            return self.parent() is other.parent()
        return False

    def __ne__(self, other):
        r"""
        EXAMPLES::

            sage: La = RootSystem("A2").weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T("A2",La[1])
            sage: T.highest_weight_vector() != T.highest_weight_vector()
            False
            sage: T.highest_weight_vector() != T.highest_weight_vector().e(1)
            True
        """
        return not self == other

    def e(self, i):
        r"""
        Return `e_i` of ``self``, which is ``None`` for all `i`.

        INPUT:

        - ``i`` -- An element of the index set

        EXAMPLES::

            sage: ct = CartanType(['A',2])
            sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T(ct,la[1])
            sage: t = T.highest_weight_vector()
            sage: t.e(1)
            sage: t.e(2)
        """
        return None

    def f(self, i):
        r"""
        Return `f_i` of ``self``, which is ``None`` for all `i`.

        INPUT:

        - ``i`` -- An element of the index set

        EXAMPLES::

            sage: ct = CartanType(['A',2])
            sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T(ct,la[1])
            sage: t = T.highest_weight_vector()
            sage: t.f(1)
            sage: t.f(2)
        """
        return None


class TCrystal(UniqueRepresentation, Parent):
    r"""
    The crystal `T_{\lambda}`.

    Let `\lambda` be a weight. As defined in [Ka1993]_ the crystal
    `T_{\lambda} = \{ t_{\lambda} \}` is a single element crystal with the
    crystal structure defined by

    .. MATH::

        \mathrm{wt}(t_\lambda) = \lambda, \quad
        e_i t_{\lambda} = f_i t_{\lambda} = 0, \quad
        \varepsilon_i(t_{\lambda}) = \varphi_i(t_{\lambda}) = -\infty.

    The crystal `T_{\lambda}` shifts the weights of the vertices in a crystal
    `B` by `\lambda` when tensored with `B`, but leaves the graph structure of
    `B` unchanged. That is to say, for all `b \in B`, we have `\mathrm{wt}(b
    \otimes t_{\lambda}) = \mathrm{wt}(b) + \lambda`.

    INPUT:

    - ``cartan_type`` -- A Cartan type

    - ``weight`` -- An element of the weight lattice of type ``cartan_type``

    EXAMPLES::

        sage: ct = CartanType(['A',2])
        sage: C = crystals.Tableaux(ct, shape=[1])
        sage: for x in C: x.weight()
        (1, 0, 0)
        (0, 1, 0)
        (0, 0, 1)
        sage: La = RootSystem(ct).ambient_space().fundamental_weights()
        sage: TLa = crystals.elementary.T(ct, 3*(La[1] + La[2]))
        sage: TP = crystals.TensorProduct(TLa, C)
        sage: for x in TP: x.weight()
        (7, 3, 0)
        (6, 4, 0)
        (6, 3, 1)
        sage: G = C.digraph()
        sage: H = TP.digraph()
        sage: G.is_isomorphic(H,edge_labels=True)
        True
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, weight=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: ct = CartanType(['A',3])
            sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
            sage: wts = RootSystem(ct).ambient_space().fundamental_weights()
            sage: X = crystals.elementary.T(['A',3], la[1])
            sage: Y = crystals.elementary.T(la[1])
            sage: X is Y
            True
        """
        if weight is None:
            weight = cartan_type
            cartan_type = weight.parent().cartan_type()
        cartan_type = CartanType(cartan_type)
        return super(TCrystal, cls).__classcall__(cls, cartan_type, weight)

    def __init__(self, cartan_type, weight):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: la = RootSystem("A2").weight_lattice().fundamental_weights()
            sage: B = crystals.elementary.T("A2", 5*la[2])
            sage: TestSuite(B).run()
        """
        self._weight = weight
        self._cartan_type = cartan_type
        if self._cartan_type.type() == 'Q':
            category = SuperCrystals().Finite()
        else:
            category = (FiniteCrystals(), HighestWeightCrystals())
        Parent.__init__(self, category=category)
        self.module_generators = (self.element_class(self),)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: la = RootSystem(['E',6]).weight_lattice().fundamental_weights()
            sage: B = crystals.elementary.T(['E',6], la[6])
            sage: B
            The T crystal of type ['E', 6] and weight Lambda[6]
        """
        return "The T crystal of type {1!s} and weight {0!s}".format(self._weight, self._cartan_type)

    def _element_constructor_(self, weight):
        r"""
        Construct an element of ``self`` from ``weight``.

        INPUT:

        - ``weight`` -- An element of the weight lattice

        EXAMPLES::

            sage: la = RootSystem("E8").weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T("E8",la[7]+la[8])
            sage: T(la[7]+la[8])
            Lambda[7] + Lambda[8]
        """
        if weight != self._weight:
            raise ValueError("Only element is t(%s)" % self._weight)
        return self.element_class(self)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``, which is always `1`.

        EXAMPLES::

            sage: La = RootSystem(['C',12]).weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T(['C',12], La[9])
            sage: T.cardinality()
            1
        """
        return ZZ.one()

    def weight_lattice_realization(self):
        """
        Return a realization of the lattice containing the weights
        of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['C',12]).weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T(['C',12], La[9])
            sage: T.weight_lattice_realization()
            Weight lattice of the Root system of type ['C', 12]

            sage: ct = CartanMatrix([[2, -4], [-5, 2]])
            sage: La = RootSystem(ct).weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T(ct, La[1])
            sage: T.weight_lattice_realization()
            Weight lattice of the Root system of type
            [ 2 -4]
            [-5  2]
        """
        return self._weight.parent()

    class Element(AbstractSingleCrystalElement):
        r"""
        Element of a `T_{\lambda}` crystal.
        """
        def _repr_(self):
            r"""
            EXAMPLES::

                sage: ct = CartanType(['F',4])
                sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
                sage: T = crystals.elementary.T(ct,2*la[1]-3*la[3])
                sage: t = T.highest_weight_vector()
                sage: t
                2*Lambda[1] - 3*Lambda[3]
            """
            return repr(self.parent()._weight)

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            EXAMPLES::

                sage: ct = CartanType(['B',5,1])
                sage: la = RootSystem(ct).ambient_space().fundamental_weights()
                sage: T = crystals.elementary.T(ct, 2*la[1]-3*la[3]+la[0])
                sage: t = T.highest_weight_vector()
                sage: latex(t)
                {t_{-e_{0} - 3 e_{1} - 3 e_{2} - 3 e_{deltacheck}}}
            """
            return "{t_{" + self.parent()._weight._latex_() + "}}"

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``, which is `-\infty` for all `i`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: ct = CartanType(['C',5])
                sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
                sage: T = crystals.elementary.T(ct,la[4]+la[5]-la[1]-la[2])
                sage: t = T.highest_weight_vector()
                sage: [t.epsilon(i) for i in T.index_set()]
                [-inf, -inf, -inf, -inf, -inf]
            """
            return float("-inf")

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``, which is `-\infty` for all `i`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: ct = CartanType(['C',5])
                sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
                sage: T = crystals.elementary.T(ct,la[4]+la[5]-la[1]-la[2])
                sage: t = T.highest_weight_vector()
                sage: [t.phi(i) for i in T.index_set()]
                [-inf, -inf, -inf, -inf, -inf]
            """
            return float("-inf")

        def weight(self):
            r"""
            Return the weight of ``self``, which is always `\lambda`.

            EXAMPLES::

                sage: ct = CartanType(['C',5])
                sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
                sage: T = crystals.elementary.T(ct,la[4]+la[5]-la[1]-la[2])
                sage: t = T.highest_weight_vector()
                sage: t.weight()
                -Lambda[1] - Lambda[2] + Lambda[4] + Lambda[5]
            """
            return self.parent()._weight


class RCrystal(UniqueRepresentation, Parent):
    r"""
    The crystal `R_{\lambda}`.

    For a fixed weight `\lambda`, the crystal `R_{\lambda} = \{ r_{\lambda} \}`
    is a single element crystal with the crystal structure defined by

    .. MATH::

        \mathrm{wt}(r_{\lambda}) = \lambda, \quad
        e_i r_{\lambda} = f_i r_{\lambda} = 0, \quad
        \varepsilon_i(r_{\lambda}) = -\langle h_i, \lambda\rangle, \quad
        \varphi_i(r_{\lambda}) = 0,

    where `\{h_i\}` are the simple coroots.

    Tensoring `R_{\lambda}` with a crystal `B` results in shifting the weights
    of the vertices in `B` by `\lambda` and may also cut a subset out of the
    original graph of `B`.  That is, `\mathrm{wt}(r_{\lambda} \otimes b) =
    \mathrm{wt}(b) + \lambda`, where `b \in B`, provided `r_{\lambda} \otimes
    b \neq 0`. For example, the crystal graph of `B(\lambda)` is the same as
    the crystal graph of `R_{\lambda} \otimes B(\infty)` generated from the
    component `r_{\lambda} \otimes u_{\infty}`.

    There is also a dual version of this crystal given by
    `R^{\vee}_{\lambda} = \{ r^{\vee}_{\lambda} \}` with the crystal
    structure defined by

    .. MATH::

        \mathrm{wt}(r^{\vee}_{\lambda}) = \lambda, \quad
        e_i r^{\vee}_{\lambda} = f_i r^{\vee}_{\lambda} = 0, \quad
        \varepsilon_i(r^{\vee}_{\lambda}) = 0, \quad
        \varphi_i(r^{\vee}_{\lambda}) = \langle h_i, \lambda\rangle.

    INPUT:

    - ``cartan_type`` -- a Cartan type
    - ``weight`` -- an element of the weight lattice of type ``cartan_type``
    - ``dual`` -- (default: ``False``) boolean

    EXAMPLES:

    We check by tensoring `R_{\lambda}` with `B(\infty)` results in a
    component of `B(\lambda)`::

        sage: B = crystals.infinity.Tableaux("A2")
        sage: R = crystals.elementary.R("A2", B.Lambda()[1]+B.Lambda()[2])
        sage: T = crystals.TensorProduct(R, B)
        sage: mg = T(R.highest_weight_vector(), B.highest_weight_vector())
        sage: S = T.subcrystal(generators=[mg])
        sage: sorted([x.weight() for x in S], key=str)
        [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 1, 1),
         (1, 1, 1), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        sage: C = crystals.Tableaux("A2", shape=[2,1])
        sage: sorted([x.weight() for x in C], key=str)
        [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 1, 1),
         (1, 1, 1), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        sage: GT = T.digraph(subset=S)
        sage: GC = C.digraph()
        sage: GT.is_isomorphic(GC, edge_labels=True)
        True
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, weight=None, dual=False):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: ct = CartanType(['A',3])
            sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
            sage: X = crystals.elementary.R(['A',3], la[1])
            sage: Y = crystals.elementary.R(la[1])
            sage: X is Y
            True
        """
        if weight is None:
            weight = cartan_type
            cartan_type = weight.parent().cartan_type()
        cartan_type = CartanType(cartan_type)
        return super(RCrystal, cls).__classcall__(cls, cartan_type, weight, dual)

    def __init__(self, cartan_type, weight, dual):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: la = RootSystem("A2").weight_lattice().fundamental_weights()
            sage: B = crystals.elementary.R("A2",5*la[2])
            sage: TestSuite(B).run()
        """
        self._weight = weight
        self._cartan_type = cartan_type
        self._dual = dual
        if self._cartan_type.type() == 'Q':
            category = SuperCrystals().Finite()
        else:
            category = (FiniteCrystals(), HighestWeightCrystals())
        Parent.__init__(self, category=category)
        self.module_generators = (self.element_class(self),)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: la = RootSystem(['E',6]).weight_lattice().fundamental_weights()
            sage: B = crystals.elementary.R(['E',6], la[6])
            sage: B
            The R crystal of weight Lambda[6] and type ['E', 6]

            sage: crystals.elementary.R(['E',6], la[1], dual=True)
            The dual R crystal of weight Lambda[1] and type ['E', 6]
        """
        dual_str = " dual" if self._dual else ""
        return "The{} R crystal of weight {} and type {}".format(
                    dual_str, self._weight, self._cartan_type)

    def _element_constructor_(self, weight):
        r"""
        Construct an element of ``self`` from ``weight``.

        INPUT:

        - ``weight`` -- An element of the weight lattice

        EXAMPLES::

            sage: la = RootSystem("E8").weight_lattice().fundamental_weights()
            sage: R = crystals.elementary.R("E8",la[7]+la[8])
            sage: R(la[7]+la[8])
            Lambda[7] + Lambda[8]
        """
        if weight != self._weight:
            raise ValueError("Only element is r(%s)" % self._weight)
        return self.element_class(self)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``, which is always `1`.

        EXAMPLES::

            sage: La = RootSystem(['C',12]).weight_lattice().fundamental_weights()
            sage: R = crystals.elementary.R(['C',12],La[9])
            sage: R.cardinality()
            1
        """
        return ZZ.one()

    def weight_lattice_realization(self):
        """
        Return a realization of the lattice containing the weights
        of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['C',12]).weight_lattice().fundamental_weights()
            sage: R = crystals.elementary.R(['C',12], La[9])
            sage: R.weight_lattice_realization()
            Weight lattice of the Root system of type ['C', 12]

            sage: ct = CartanMatrix([[2, -4], [-5, 2]])
            sage: La = RootSystem(ct).weight_lattice().fundamental_weights()
            sage: R = crystals.elementary.R(ct, La[1])
            sage: R.weight_lattice_realization()
            Weight lattice of the Root system of type
            [ 2 -4]
            [-5  2]
        """
        return self._weight.parent()

    class Element(AbstractSingleCrystalElement):
        r"""
        Element of a `R_{\lambda}` crystal.
        """
        def _repr_(self):
            r"""
            EXAMPLES::

                sage: ct = CartanType(['F',4])
                sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
                sage: T = crystals.elementary.T(ct,2*la[1]-3*la[3])
                sage: t = T.highest_weight_vector()
                sage: t
                2*Lambda[1] - 3*Lambda[3]
            """
            return repr(self.parent()._weight)

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            EXAMPLES::

                sage: la = RootSystem("G2").weight_lattice().fundamental_weights()
                sage: R = crystals.elementary.R("G2",la[1])
                sage: r = R.highest_weight_vector()
                sage: latex(r)
                {r_{\Lambda_{1}}}

                sage: R = crystals.elementary.R("G2", la[1], dual=True)
                sage: latex(R.highest_weight_vector())
                {r^{\vee}_{\Lambda_{1}}}
            """
            if self.parent()._dual:
                return r"{r^{\vee}_{" + self.parent()._weight._latex_() + "}}"
            return "{r_{" + self.parent()._weight._latex_() + "}}"

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            We have `\varepsilon_i(r_{\lambda}) = -\langle h_i, \lambda
            \rangle` for all `i`, where `h_i` is a simple coroot.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: la = RootSystem(['A',2]).weight_lattice().fundamental_weights()
                sage: R = crystals.elementary.R("A2", la[1])
                sage: r = R.highest_weight_vector()
                sage: [r.epsilon(i) for i in R.index_set()]
                [-1, 0]

                sage: R = crystals.elementary.R("A2", la[1], dual=True)
                sage: r = R.highest_weight_vector()
                sage: [r.epsilon(i) for i in R.index_set()]
                [0, 0]
            """
            if self.parent()._dual:
                return ZZ.zero()
            else:
                P = self.parent().weight_lattice_realization()
                h = P.simple_coroots()
                return -P(self.weight()).scalar(h[i])

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``, which is `0` for all `i`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: la = RootSystem("C5").weight_lattice().fundamental_weights()
                sage: R = crystals.elementary.R("C5", la[4]+la[5])
                sage: r = R.highest_weight_vector()
                sage: [r.phi(i) for i in R.index_set()]
                [0, 0, 0, 0, 0]

                sage: R = crystals.elementary.R("C5", la[4]+la[5], dual=True)
                sage: r = R.highest_weight_vector()
                sage: [r.phi(i) for i in R.index_set()]
                [0, 0, 0, 1, 1]
            """
            if self.parent()._dual:
                P = self.parent().weight_lattice_realization()
                h = P.simple_coroots()
                return P(self.weight()).scalar(h[i])
            else:
                return ZZ.zero()

        def weight(self):
            r"""
            Return the weight of ``self``, which is always `\lambda`.

            EXAMPLES::

                sage: ct = CartanType(['C',5])
                sage: la = RootSystem(ct).weight_lattice().fundamental_weights()
                sage: T = crystals.elementary.T(ct,la[4]+la[5]-la[1]-la[2])
                sage: t = T.highest_weight_vector()
                sage: t.weight()
                -Lambda[1] - Lambda[2] + Lambda[4] + Lambda[5]
            """
            return self.parent()._weight


class ElementaryCrystal(UniqueRepresentation, Parent):
    r"""
    The elementary crystal `B_i`.

    For `i` an element of the index set of type `X`, the crystal `B_i` of type
    `X` is the set

    .. MATH::

        B_i = \{ b_i(m) : m \in \ZZ \},

    where the crystal structure is given by

    .. MATH::

        \begin{aligned}
        \mathrm{wt}\bigl(b_i(m)\bigr) &= m\alpha_i \\
        \varphi_j\bigl(b_i(m)\bigr) &= \begin{cases}
            m & \text{ if } j=i, \\
            -\infty & \text{ if } j\neq i,
        \end{cases} \\
        \varepsilon_j\bigl(b_i(m)\bigr) &= \begin{cases}
            -m & \text{ if } j=i, \\
            -\infty & \text{ if } j\neq i,
        \end{cases} \\
        e_j b_i(m) &= \begin{cases}
            b_i(m+1) & \text{ if } j=i, \\
            0 & \text{ if } j\neq i,
        \end{cases} \\
        f_j b_i(m) &= \begin{cases}
            b_i(m-1) & \text{ if } j=i, \\
            0 & \text{ if } j\neq i.
        \end{cases}
        \end{aligned}

    The *Kashiwara embedding theorem* asserts there is a unique strict crystal
    embedding of crystals

    .. MATH::

        B(\infty) \hookrightarrow B_i \otimes B(\infty),

    satisfying certain properties (see [Ka1993]_).  The above embedding
    may be iterated to obtain a new embedding

    .. MATH::

        B(\infty) \hookrightarrow B_{i_N} \otimes B_{i_{N-1}}
        \otimes \cdots \otimes B_{i_2} \otimes B_{i_1} \otimes B(\infty),

    which is a foundational object in the study of *polyhedral realizations of
    crystals* (see, for example, [NZ1997]_).
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, i):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: B = crystals.elementary.Elementary(['A',4], 3)
            sage: C = crystals.elementary.Elementary(CartanType("A4"), int(3))
            sage: B is C
            True
        """
        cartan_type = CartanType(cartan_type)
        if i not in cartan_type.index_set():
            raise ValueError('i must an element of the index set.')
        return super(ElementaryCrystal, cls).__classcall__(cls, cartan_type, i)

    def __init__(self, cartan_type, i):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.elementary.Elementary("D4",3)
            sage: TestSuite(B).run()
        """
        Parent.__init__(self, category=(Crystals(), InfiniteEnumeratedSets()))
        self._i = i
        self._cartan_type = cartan_type
        self.module_generators = (self.element_class(self, 0),)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.elementary.Elementary(['B',5,1], 4)
            sage: B
            The 4-elementary crystal of type ['B', 5, 1]
        """
        return "The {0!s}-elementary crystal of type {1!s}".format(self._i, self._cartan_type)

    def _element_constructor_(self, m):
        r"""
        Construct an element of ``self`` from ``weight``.

        INPUT:

        - ``m`` -- An integer

        EXAMPLES::

            sage: B = crystals.elementary.Elementary(['F',4], 2)
            sage: B(0)
            0
            sage: B(-15)
            -15
            sage: B(721)
            721
        """
        return self.element_class(self, ZZ(m))

    def weight_lattice_realization(self):
        """
        Return a realization of the lattice containing the weights
        of ``self``.

        EXAMPLES::

            sage: B = crystals.elementary.Elementary(['A',4, 1], 2)
            sage: B.weight_lattice_realization()
            Root lattice of the Root system of type ['A', 4, 1]
        """
        return self.cartan_type().root_system().root_lattice()

    class Element(Element):
        r"""
        Element of a `B_i` crystal.
        """
        def __init__(self, parent, m):
            r"""
            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['B',7],7)
                sage: elt = B(17); elt
                17
            """
            self._m = m
            Element.__init__(self, parent)

        def __hash__(self):
            r"""
            TESTS::

                sage: B = crystals.elementary.Elementary(['B',7],7)
                sage: hash(B(17))
                17
            """
            return hash(self._m)

        def _repr_(self):
            r"""
            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['A',4],3)
                sage: B(-47)
                -47
            """
            return repr(self._m)

        def __lt__(self, other):
            r"""
            EXAMPLES::

                sage: B = crystals.elementary.Elementary("D4",3)
                sage: b = B(1)
                sage: c = B(-1)
                sage: b.__lt__(c)
                False
                sage: c.__lt__(b)
                True
            """
            if self.parent() is not other.parent():
                return False
            return Integer(self._m) < Integer(other._m)

        def __eq__(self, other):
            r"""
            EXAMPLES::

                sage: B = crystals.elementary.Elementary("A2",1)
                sage: C = crystals.elementary.Elementary("A2",2)
                sage: D = crystals.elementary.Elementary("B2",1)
                sage: [B(0) == B(1), B(0) == C(0), B(0) == D(0), C(0) == D(0)]
                [False, False, False, False]
                sage: [B(1) == B(1), C(12) == C(12), D(-1) == D(-1)]
                [True, True, True]
            """
            if isinstance(other, ElementaryCrystal.Element):
                return self.parent() is other.parent() and self._m == other._m
            return False

        def __ne__(self, other):
            r"""
            EXAMPLES::

                sage: B = crystals.elementary.Elementary("A2",1)
                sage: B(0) != B(2)
                True
                sage: B(0) != B(0)
                False
            """
            return not self == other

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['B',11,1],6)
                sage: latex(B(26))
                {b_{6}(26)}
            """
            return "{b_{%s}(%s)}" % (self.parent()._i, self._m)

        def e(self, i):
            r"""
            Return the action of `e_i` on ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['E',7],1)
                sage: B(3).e(1)
                4
                sage: B(172).e_string([1]*171)
                343
                sage: B(0).e(2)
            """
            if i == self.parent()._i:
                return self.__class__(self.parent(), self._m + 1)
            else:
                return None

        def f(self, i):
            r"""
            Return the action of `f_i` on ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['E',7],1)
                sage: B(3).f(1)
                2
                sage: B(172).f_string([1]*171)
                1
                sage: B(0).e(2)
            """
            if i == self.parent()._i:
                return self.__class__(self.parent(), self._m - 1)
            else:
                return None

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['F',4],3)
                sage: [[B(j).epsilon(i) for i in B.index_set()] for j in range(5)]
                [[-inf, -inf, 0, -inf],
                 [-inf, -inf, -1, -inf],
                 [-inf, -inf, -2, -inf],
                 [-inf, -inf, -3, -inf],
                 [-inf, -inf, -4, -inf]]
            """
            if i == self.parent()._i:
                return -self._m
            else:
                return float("-inf")

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['E',8,1],4)
                sage: [[B(m).phi(j) for j in B.index_set()] for m in range(44,49)]
                [[-inf, -inf, -inf, -inf, 44, -inf, -inf, -inf, -inf],
                 [-inf, -inf, -inf, -inf, 45, -inf, -inf, -inf, -inf],
                 [-inf, -inf, -inf, -inf, 46, -inf, -inf, -inf, -inf],
                 [-inf, -inf, -inf, -inf, 47, -inf, -inf, -inf, -inf],
                 [-inf, -inf, -inf, -inf, 48, -inf, -inf, -inf, -inf]]
            """
            if i == self.parent()._i:
                return self._m
            else:
                return float("-inf")

        def weight(self):
            r"""
            Return the weight of ``self``.

            EXAMPLES::

                sage: B = crystals.elementary.Elementary(['C',14],12)
                sage: B(-385).weight()
                -385*alpha[12]
            """
            Q = self.parent().weight_lattice_realization()
            return self._m * Q.simple_root(self.parent()._i)


class ComponentCrystal(UniqueRepresentation, Parent):
    r"""
    The component crystal.

    Defined in [Ka1993]_, the component crystal `C = \{c\}` is the single
    element crystal whose crystal structure is defined by

    .. MATH::

        \mathrm{wt}(c) = 0, \quad
        e_i c = f_i c = 0, \quad
        \varepsilon_i(c) = \varphi_i(c) = 0.

    Note `C \cong B(0)`, where `B(0)` is the highest weight crystal of highest
    weight `0`.

    INPUT:

    - ``cartan_type`` -- a Cartan type
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, P=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: C = crystals.elementary.Component("A2")
            sage: D = crystals.elementary.Component(CartanType(['A',2]))
            sage: C is D
            True
            sage: AS = RootSystem(['A',2]).ambient_space()
            sage: E = crystals.elementary.Component(AS)
            sage: F = crystals.elementary.Component(CartanType(['A',2]), AS)
            sage: C is E and C is F
            True
        """
        if cartan_type in RootLatticeRealizations:
            P = cartan_type
            cartan_type = P.cartan_type()
        elif P is None:
            cartan_type = CartanType(cartan_type)
            P = cartan_type.root_system().ambient_space()
            if P is None:
                P = cartan_type.root_system().weight_lattice()
        return super(ComponentCrystal, cls).__classcall__(cls, cartan_type, P)

    def __init__(self, cartan_type, P):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.elementary.Component("D4")
            sage: TestSuite(B).run()
        """
        self._weight_lattice_realization = P
        self._cartan_type = cartan_type
        if self._cartan_type.type() == 'Q':
            category = RegularSuperCrystals()
        else:
            category = (RegularCrystals().Finite(), HighestWeightCrystals())
        Parent.__init__(self, category=category)
        self.module_generators = (self.element_class(self),)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.elementary.Component("D4")
            sage: C
            The component crystal of type ['D', 4]
        """
        return "The component crystal of type {0!s}".format(self._cartan_type)

    def _element_constructor_(self, weight):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: C = crystals.elementary.Component("E6")
            sage: c = C.highest_weight_vector()
            sage: c
            c
        """
        if weight != self.weight_lattice_realization().zero():
            raise ValueError("only element is c")
        return self.element_class(self)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``, which is always `1`.

        EXAMPLES::

            sage: C = crystals.elementary.Component("E6")
            sage: c = C.highest_weight_vector()
            sage: C.cardinality()
            1
        """
        return ZZ.one()

    def weight_lattice_realization(self):
        """
        Return the weight lattice realization of ``self``.

        EXAMPLES::

            sage: C = crystals.elementary.Component("A2")
            sage: C.weight_lattice_realization()
            Ambient space of the Root system of type ['A', 2]

            sage: P = RootSystem(['A',2]).weight_lattice()
            sage: C = crystals.elementary.Component(P)
            sage: C.weight_lattice_realization() is P
            True
        """
        return self._weight_lattice_realization

    class Element(AbstractSingleCrystalElement):
        r"""
        Element of a component crystal.
        """
        def _repr_(self):
            r"""
            EXAMPLES::

                sage: C = crystals.elementary.Component("F4")
                sage: c = C.highest_weight_vector()
                sage: c
                c
            """
            return 'c'

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            EXAMPLES::

                sage: C = crystals.elementary.Component("E7")
                sage: c = C.highest_weight_vector()
                sage: latex(c)
                {c}
            """
            return "{c}"

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``, which is `0` for all `i`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: C = crystals.elementary.Component("C5")
                sage: c = C.highest_weight_vector()
                sage: [c.epsilon(i) for i in C.index_set()]
                [0, 0, 0, 0, 0]
            """
            return 0

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``, which is `0` for all `i`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: C = crystals.elementary.Component("C5")
                sage: c = C.highest_weight_vector()
                sage: [c.phi(i) for i in C.index_set()]
                [0, 0, 0, 0, 0]
            """
            return 0

        def weight(self):
            r"""
            Return the weight of ``self``, which is always `0`.

            EXAMPLES::

                sage: C = crystals.elementary.Component("F4")
                sage: c = C.highest_weight_vector()
                sage: c.weight()
                (0, 0, 0, 0)
            """
            return self.parent().weight_lattice_realization().zero()
