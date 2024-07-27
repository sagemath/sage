# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Kyoto Path Model for Affine Highest Weight Crystals
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
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
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.parent import Parent
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals, \
    TensorProductOfRegularCrystalsElement


class KyotoPathModel(TensorProductOfCrystals):
    r"""
    The Kyoto path model for an affine highest weight crystal.

    .. NOTE::

        Here we are using anti-Kashiwara notation and might differ from
        some of the literature.

    Consider a Kac--Moody algebra `\mathfrak{g}` of affine Cartan type `X`,
    and we want to model the `U_q'(\mathfrak{g})`-crystal `B(\lambda)`.
    First we consider the set of fundamental weights `\{\Lambda_i\}_{i \in I}`
    of `\mathfrak{g}` and let `\{\overline{\Lambda}_i\}_{i \in I_0}` be the
    corresponding fundamental weights of the corresponding classical Lie
    algebra `\mathfrak{g}_0`. To model `B(\lambda)`, we start with a sequence
    of perfect `U_q'(\mathfrak{g})`-crystals `(B^{(i)})_i` of level
    `l` such that

    .. MATH::

        \lambda \in \overline{P}_l^+ = \left\{ \mu \in \overline{P}^+ \mid
        \langle c, \mu \rangle = l \right\}

    where `c` is the canonical central element of `U_q'(\mathfrak{g})`
    and `\overline{P}^+` is the nonnegative weight lattice spanned by
    `\{ \overline{\Lambda}_i \mid i \in I \}`.

    Next we consider the crystal isomorphism `\Phi_0 : B(\lambda_0) \to B^{(0)}
    \otimes B(\lambda_1)` defined by `u_{\lambda_0} \mapsto b^{(0)}_{\lambda_0}
    \otimes u_{\lambda_1}` where `b^{(0)}_{\lambda_0}` is the unique element in
    `B^{(0)}` such that `\varphi\left( b^{(0)}_{\lambda_0} \right) = \lambda_0`
    and `\lambda_1 = \varepsilon\left( b^{(0)}_{\lambda_0} \right)` and
    `u_{\mu}` is the highest weight element in `B(\mu)`. Iterating this, we
    obtain the following isomorphism:

    .. MATH::

        \Phi_n : B(\lambda) \to B^{(0)} \otimes B^{(1)} \otimes \cdots
        \otimes B^{(N)} \otimes B(\lambda_{N+1}).

    We note by Lemma 10.6.2 in [HK2002]_ that for any `b \in B(\lambda)` there
    exists a finite `N` such that

    .. MATH::

        \Phi_N(b) = \left( \bigotimes_{k=0}^{N-1} b^{(k)} \right)
        \otimes u_{\lambda_N}.

    Therefore we can model elements `b \in B(\lambda)` as a
    `U_q'(\mathfrak{g})`-crystal by considering an infinite list of
    elements `b^{(k)} \in B^{(k)}` and defining the crystal structure by:

    .. MATH::

        \begin{aligned}
        \overline{\mathrm{wt}}(b) & = \lambda_N + \sum_{k=0}^{N-1}
        \overline{\mathrm{wt}}\left( b^{(k)} \right)
        \\ e_i(b) & = e_i\left( b^{\prime} \otimes b^{(N)} \right) \otimes
        u_{\lambda_N},
        \\ f_i(b) & = f_i\left( b^{\prime} \otimes b^{(N)} \right) \otimes
        u_{\lambda_N},
        \\ \varepsilon_i(b) & = \max\left( \varepsilon_i(b^{\prime}) -
        \varphi_i\left( b^{(N)} \right), 0 \right),
        \\ \varphi_i(b) & = \varphi_i(b^{\prime}) + \max\left(
        \varphi_i\left( b^{(N)} \right) - \varepsilon_i(b^{\prime}), 0 \right),
        \end{aligned}

    where `b^{\prime} = b^{(0)} \otimes \cdots \otimes b^{(N-1)}`. To
    translate this into a finite list, we consider a finite sequence
    `b^{(0)} \otimes \cdots \otimes b^{(N-1)} \otimes b^{(N)}_{\lambda_N}`
    and if

    .. MATH::

        f_i\left( b^{(0)} \otimes \cdots b^{(N-1)} \otimes
        b^{(N)}_{\lambda_N} \right) = b_0 \otimes \cdots \otimes b^{(N-1)}
        \otimes f_i\left( b^{(N)}_{\lambda_N} \right),

    then we take the image as `b^{(0)} \otimes \cdots \otimes f_i\left(
    b^{(N)}_{\lambda_N}\right) \otimes b^{(N+1)}_{\lambda_{N+1}}`. Similarly
    we remove `b^{(N)}_{\lambda_{N}}` if we have `b_0 \otimes \cdots
    \otimes b^{(N-1)} \otimes b^{(N-1)}_{\lambda_{N-1}} \otimes
    b^{(N)}_{\lambda_N}`. Additionally if

    .. MATH::

        e_i\left( b^{(0)} \otimes \cdots \otimes b^{(N-1)} \otimes
        b^{(N)}_{\lambda_N} \right) = b^{(0)} \otimes \cdots \otimes
        b^{(N-1)} \otimes e_i\left( b^{(N)}_{\lambda_N} \right),

    then we consider this to be `0`.

    We can then lift the `U_q'(\mathfrak{g})`-crystal structure to a
    `U_q(\mathfrak{g})`-crystal structure by using a tensor product of
    the :class:`affinization
    <sage.combinat.crystals.affinization.AffinizationOfCrystal>` of the
    of crystals `B^{(i)}` for all `i`.

    INPUT:

    - ``B`` -- a single or list of `U_q^{\prime}` perfect crystal(s) of
      level `l`
    - ``weight`` -- a weight in `\overline{P}_l^+`

    EXAMPLES::

        sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
        sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
        sage: C = crystals.KyotoPathModel(B, La[0])
        sage: mg = C.module_generators[0]; mg
        [[[3]]]
        sage: mg.f_string([0,1,2,2])
        [[[3]], [[3]], [[1]]]
        sage: x = mg.f_string([0,1,2]); x
        [[[2]], [[3]], [[1]]]
        sage: x.weight()
        Lambda[0]

    An example of type `A_5^{(2)}`::

        sage: B = crystals.KirillovReshetikhin(['A',5,2], 1,1)
        sage: La = RootSystem(['A',5,2]).weight_lattice().fundamental_weights()
        sage: C = crystals.KyotoPathModel(B, La[0])
        sage: mg = C.module_generators[0]; mg
        [[[-1]]]
        sage: mg.f_string([0,2,1,3])
        [[[-3]], [[2]], [[-1]]]
        sage: mg.f_string([0,2,3,1])
        [[[-3]], [[2]], [[-1]]]

    An example of type `D_3^{(2)}`::

        sage: B = crystals.KirillovReshetikhin(['D',3,2], 1,1)
        sage: La = RootSystem(['D',3,2]).weight_lattice().fundamental_weights()
        sage: C = crystals.KyotoPathModel(B, La[0])
        sage: mg = C.module_generators[0]; mg
        [[]]
        sage: mg.f_string([0,1,2,0])
        [[[0]], [[1]], []]

    An example using multiple crystals of the same level::

        sage: B1 = crystals.KirillovReshetikhin(['A',2,1], 1,1)
        sage: B2 = crystals.KirillovReshetikhin(['A',2,1], 2,1)
        sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
        sage: C = crystals.KyotoPathModel([B1, B2, B1], La[0])
        sage: mg = C.module_generators[0]; mg
        [[[3]]]
        sage: mg.f_string([0,1,2,2])
        [[[3]], [[1], [3]], [[3]]]
        sage: mg.f_string([0,1,2,2,2])
        sage: mg.f_string([0,1,2,2,1,0])
        [[[3]], [[2], [3]], [[1]], [[2]]]
        sage: mg.f_string([0,1,2,2,1,0,0,2])
        [[[3]], [[1], [2]], [[1]], [[3]], [[1], [3]]]

    By using the extended weight lattice, the Kyoto path model lifts
    the perfect crystals to their affinizations::

        sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
        sage: P = RootSystem(['A',2,1]).weight_lattice(extended=True)
        sage: La = P.fundamental_weights()
        sage: C = crystals.KyotoPathModel(B, La[0])
        sage: mg = C.module_generators[0]; mg
        [[[3]](0)]
        sage: x = mg.f_string([0,1,2]); x
        [[[2]](-1), [[3]](0), [[1]](0)]
        sage: x.weight()
        Lambda[0] - delta
    """
    @staticmethod
    def __classcall_private__(cls, crystals, weight, P=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
            sage: P = RootSystem(['A',2,1]).weight_lattice()
            sage: La = P.fundamental_weights()
            sage: C = crystals.KyotoPathModel(B, La[0])
            sage: C2 = crystals.KyotoPathModel((B,), La[0])
            sage: C3 = crystals.KyotoPathModel([B], La[0], P)
            sage: C is C2 and C2 is C3
            True
        """
        if isinstance(crystals, list):
            crystals = tuple(crystals)
        elif not isinstance(crystals, tuple):
            crystals = (crystals,)

        if any(not B.is_perfect() for B in crystals):
            raise ValueError("all crystals must be perfect")
        level = crystals[0].level()
        if any(B.level() != level for B in crystals[1:]):
            raise ValueError("all crystals must have the same level")
        ct = crystals[0].cartan_type()
        if P is None:
            P = weight.parent()
        if sum( ct.dual().c()[i] * weight.scalar(h) for i,h in
                enumerate(P.simple_coroots()) ) != level:
            raise ValueError( "{} is not a level {} weight".format(weight, level) )

        return super().__classcall__(cls, crystals, weight, P)

    def __init__(self, crystals, weight, P):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
            sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
            sage: C = crystals.KyotoPathModel(B, La[0])
            sage: TestSuite(C).run() # long time
        """
        Parent.__init__(self, category=(HighestWeightCrystals(), InfiniteEnumeratedSets()))

        self._cartan_type = crystals[0].cartan_type()
        self._weight = weight
        if weight.parent().is_extended():
            # public for TensorProductOfCrystals
            self.crystals = tuple([C.affinization() for C in crystals])
            self._epsilon_dicts = [{b.Epsilon(): self.crystals[i](b, 0) for b in B}
                                   for i,B in enumerate(crystals)]
            self._phi_dicts = [{b.Phi(): self.crystals[i](b, 0) for b in B}
                               for i,B in enumerate(crystals)]
        else:
            # public for TensorProductOfCrystals
            self.crystals = tuple(crystals)
            self._epsilon_dicts = [{b.Epsilon(): b for b in B}
                                   for B in crystals]
            self._phi_dicts = [{b.Phi(): b for b in B}
                               for B in crystals]
        self.module_generators = (self.element_class(self, [self._phi_dicts[0][weight]]),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
            sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
            sage: crystals.KyotoPathModel(B, La[0])
            Kyoto path realization of B(Lambda[0]) using
             [Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1)]
        """
        return "Kyoto path realization of B({}) using {}".format(self._weight, list(self.crystals))

    def finite_tensor_product(self, k):
        """
        Return the finite tensor product of crystals of length ``k``
        from truncating ``self``.

        EXAMPLES::

            sage: B1 = crystals.KirillovReshetikhin(['A',2,1], 1,1)
            sage: B2 = crystals.KirillovReshetikhin(['A',2,1], 2,1)
            sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
            sage: C = crystals.KyotoPathModel([B1,B2,B1], La[0])
            sage: C.finite_tensor_product(5)
            Full tensor product of the crystals
             [Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1),
              Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(2,1),
              Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1),
              Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1),
              Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(2,1)]
        """
        N = len(self.crystals)
        crystals = [self.crystals[i % N] for i in range(k)]
        return TensorProductOfCrystals(*crystals)

    def weight_lattice_realization(self):
        """
        Return the weight lattice realization used to express weights.

        EXAMPLES::

            sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
            sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
            sage: C = crystals.KyotoPathModel(B, La[0])
            sage: C.weight_lattice_realization()
            Weight lattice of the Root system of type ['A', 2, 1]

            sage: P = RootSystem(['A',2,1]).weight_lattice(extended=True)
            sage: C = crystals.KyotoPathModel(B, P.fundamental_weight(0))
            sage: C.weight_lattice_realization()
            Extended weight lattice of the Root system of type ['A', 2, 1]
        """
        return self._weight.parent()

    class Element(TensorProductOfRegularCrystalsElement):
        """
        An element in the Kyoto path model.
        """
        # For simplicity (and safety), we use the regular crystals implementation

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
                sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
                sage: C = crystals.KyotoPathModel(B, La[0])
                sage: mg = C.module_generators[0]
                sage: [mg.epsilon(i) for i in C.index_set()]
                [0, 0, 0]
                sage: elt = mg.f(0)
                sage: [elt.epsilon(i) for i in C.index_set()]
                [1, 0, 0]
                sage: elt = mg.f_string([0,1,2])
                sage: [elt.epsilon(i) for i in C.index_set()]
                [0, 0, 1]
                sage: elt = mg.f_string([0,1,2,2])
                sage: [elt.epsilon(i) for i in C.index_set()]
                [0, 0, 2]
            """
            x = self.e(i)
            eps = 0
            while x is not None:
                x = x.e(i)
                eps = eps + 1
            return eps

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
                sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
                sage: C = crystals.KyotoPathModel(B, La[0])
                sage: mg = C.module_generators[0]
                sage: [mg.phi(i) for i in C.index_set()]
                [1, 0, 0]
                sage: elt = mg.f(0)
                sage: [elt.phi(i) for i in C.index_set()]
                [0, 1, 1]
                sage: elt = mg.f_string([0,1])
                sage: [elt.phi(i) for i in C.index_set()]
                [0, 0, 2]
            """
            x = self.f(i)
            phi = 0
            while x is not None:
                x = x.f(i)
                phi = phi + 1
            return phi

        def e(self, i):
            """
            Return the action of `e_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
                sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
                sage: C = crystals.KyotoPathModel(B, La[0])
                sage: mg = C.module_generators[0]
                sage: all(mg.e(i) is None for i in C.index_set())
                True
                sage: mg.f(0).e(0) == mg
                True
            """
            k = self.position_of_first_unmatched_plus(i)
            if k is None:
                return None
            if k == len(self)-1:
                return None
            crystal = self[k].e(i)
            if k == len(self)-2 and crystal.Epsilon() == self[-1].Phi():
                l = self[:-1]
                l[-1] = crystal
                return self.__class__(self.parent(), l)
            return self._set_index(k, crystal)

        def f(self, i):
            """
            Return the action of `f_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
                sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
                sage: C = crystals.KyotoPathModel(B, La[0])
                sage: mg = C.module_generators[0]
                sage: mg.f(2)
                sage: mg.f(0)
                [[[1]], [[2]]]
                sage: mg.f_string([0,1,2])
                [[[2]], [[3]], [[1]]]
            """
            k = self.position_of_last_unmatched_minus(i)
            if k is None:
                return None
            if k == len(self)-1:
                l = list(self)
                k = len(l) % len(self.parent().crystals)
                l.append(self.parent()._phi_dicts[k][ l[-1].Epsilon() ])
                l[-2] = l[-2].f(i)
                return self.__class__(self.parent(), l)
            return self._set_index(k, self[k].f(i))

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: B = crystals.KirillovReshetikhin(['A',2,1], 1,1)
                sage: P = RootSystem(['A',2,1]).weight_lattice(extended=True)
                sage: La = P.fundamental_weights()
                sage: C = crystals.KyotoPathModel(B, La[0])
                sage: mg = C.module_generators[0]
                sage: mg.weight()
                Lambda[0]
                sage: mg.f_string([0,1,2]).weight()
                Lambda[0] - delta
            """
            wt = TensorProductOfRegularCrystalsElement.weight(self)
            return wt + self[-1].Epsilon()

        def truncate(self, k=None):
            r"""
            Truncate ``self`` to have length ``k`` and return as an element
            in a (finite) tensor product of crystals.

            INPUT:

            - ``k`` -- (optional) the length to truncate to; if not specified,
              then returns one more than the current non-ground-state elements
              (i.e. the current list in ``self``)

            EXAMPLES::

                sage: B1 = crystals.KirillovReshetikhin(['A',2,1], 1,1)
                sage: B2 = crystals.KirillovReshetikhin(['A',2,1], 2,1)
                sage: La = RootSystem(['A',2,1]).weight_lattice().fundamental_weights()
                sage: C = crystals.KyotoPathModel([B1,B2,B1], La[0])
                sage: mg = C.highest_weight_vector()
                sage: elt = mg.f_string([0,1,2,2,1,0]); elt
                [[[3]], [[2], [3]], [[1]], [[2]]]
                sage: t = elt.truncate(); t
                [[[3]], [[2], [3]], [[1]], [[2]]]
                sage: t.parent() is C.finite_tensor_product(4)
                True
                sage: elt.truncate(2)
                [[[3]], [[2], [3]]]
                sage: elt.truncate(10)
                [[[3]], [[2], [3]], [[1]], [[2]], [[1], [3]],
                 [[2]], [[1]], [[2], [3]], [[1]], [[3]]]
            """
            if k is None:
                k = len(self)

            P = self.parent().finite_tensor_product(k)
            if k <= len(self):
                l = self[:k]
            else:
                l = list(self)
                N = len(self.parent().crystals)
                while len(l) < k:
                    i = len(l) % N
                    l.append(self.parent()._phi_dicts[i][ l[-1].Epsilon() ])
            return P(*l)
