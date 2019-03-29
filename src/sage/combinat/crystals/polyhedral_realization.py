r"""
Polyhedral Realization of `B(\infty)`
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
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
    TensorProductOfCrystalsElement
from sage.combinat.crystals.elementary_crystals import ElementaryCrystal
from sage.combinat.root_system.cartan_type import CartanType

class InfinityCrystalAsPolyhedralRealization(TensorProductOfCrystals):
    r"""
    The polyhedral realization of `B(\infty)`.

    .. NOTE::

        Here we are using anti-Kashiwara notation and might differ from
        some of the literature.

    Consider a Kac-Moody algebra `\mathfrak{g}` of Cartan type `X` with
    index set `I`,  and consider a finite sequence `J = (j_1, j_2, \ldots, j_m)`
    whose support equals `I`.  We extend this to an infinite sequence
    by taking `\bar{J} = J \cdot J \cdot J \cdots`, where `\cdot` denotes
    concatenation of sequences. Let

    .. MATH::

        B_J = B_{j_m} \otimes \cdots \otimes B_{j_2} \otimes B_{j_1},

    where `B_i` is an
    :class:`~sage.combinat.crystals.elementary_crystals.ElementaryCrystal`.

    As given in Theorem 2.1.1 of [Ka1993]_, there exists a strict crystal embedding
    `\Psi_i \colon B(\infty) \to B_i \otimes B(\infty)` defined by `u_{\infty}
    \mapsto b_i(0) \otimes u_{\infty}`, where `b_i(0) \in B_i` and `u_{\infty}`
    is the (unique) highest weight element in `B(\infty)`. This is sometimes
    known as the *Kashiwara embedding* [NZ1997]_ (though, in [NZ1997]_, the target
    of this map is denoted by `\ZZ_J^\infty`).  By iterating this embedding by
    taking `\Psi_J = \Psi_{j_n} \circ \Psi_{j_{n-1}} \circ \cdots \circ
    \Psi_{j_1}`, we obtain the following strict crystal embedding:

    .. MATH::

        \Psi_J^n \colon B(\infty) \to B_J^{\otimes n} \otimes B(\infty).

    We note there is a natural analog of Lemma 10.6.2 in [HK2002]_ that
    for any `b \in B(\infty)`, there exists a positive integer `N` such that

    .. MATH::

        \Psi^N_J(b) = \left( \bigotimes_{k=1}^N b^{(k)} \right)
        \otimes u_{\infty}.

    Therefore we can model elements `b \in B(\infty)` by considering
    an infinite list of elements `b^{(k)} \in B_J` and defining the crystal
    structure by:

    .. MATH::

        \begin{aligned}
        \mathrm{wt}(b) & = \sum_{k=1}^N \mathrm{wt}(b^{(k)})
        \\ e_i(b) & = e_i\left( \left( \bigotimes_{k=1}^N b^{(k)} \right)
        \right) \otimes u_{\infty},
        \\ f_i(b) & = f_i\left( \left( \bigotimes_{k=1}^N b^{(k)} \right)
        \right) \otimes u_{\infty},
        \\ \varepsilon_i(b) & = \max_{ e_i^k(b) \neq 0 } k,
        \\ \varphi_i(b) & = \varepsilon_i(b) - \langle \mathrm{wt}(b),
        h_i^{\vee} \rangle.
        \end{aligned}

    To translate this into a finite list, we consider a finite sequence
    `b_1 \otimes \cdots \otimes b_N` and if

    .. MATH::

        f_i\left( b^{(1)} \otimes \cdots b^{(N-1)} \otimes b^{(N)} \right)
        = b^{(1)} \otimes \cdots \otimes b^{(N-1)} \otimes
        f_i\left( b^{(N)} \right),

    then we take the image as `b^{(1)} \otimes \cdots \otimes f_i\left(
    b^{(N)} \right) \otimes b^{(N+1)}`. Similarly we remove `b^{(N)}` if
    we have `b^{(N)} = \bigotimes_{k=1}^m b_{j_k}(0)`. Additionally if

    .. MATH::

        e_i\left( b^{(1)} \otimes \cdots \otimes b^{(N-1)} \otimes
        b^{(N)} \right) = b^{(1)} \otimes \cdots \otimes b^{(N-1)}
        \otimes e_i\left( b^{(N)} \right),

    then we consider this to be `0`.

    INPUT:

    - ``cartan_type`` -- a Cartan type
    - ``seq`` -- (default: ``None``) a finite sequence whose support
      equals the index set of the Cartan type; if ``None``, then this
      is the index set

    EXAMPLES::

        sage: B = crystals.infinity.PolyhedralRealization(['A',2])
        sage: mg = B.module_generators[0]; mg
        [0, 0]
        sage: mg.f_string([2,1,2,2])
        [0, -3, -1, 0, 0, 0]

    An example of type `B_2`::

        sage: B = crystals.infinity.PolyhedralRealization(['B',2])
        sage: mg = B.module_generators[0]; mg
        [0, 0]
        sage: mg.f_string([2,1,2,2])
        [0, -2, -1, -1, 0, 0]

    An example of type `G_2`::

        sage: B = crystals.infinity.PolyhedralRealization(['G',2])
        sage: mg = B.module_generators[0]; mg
        [0, 0]
        sage: mg.f_string([2,1,2,2])
        [0, -3, -1, 0, 0, 0]
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, seq=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: B1 = crystals.infinity.PolyhedralRealization(['A',2])
            sage: B2 = crystals.infinity.PolyhedralRealization(['A',2], [1,2])
            sage: B1 is B2
            True
        """
        cartan_type = CartanType(cartan_type)
        if seq is None:
            seq = cartan_type.index_set()
        else:
            seq = tuple(seq)
        if set(seq) != set(cartan_type.index_set()):
            raise ValueError("the support of seq is not the index set")
        return super(InfinityCrystalAsPolyhedralRealization, cls).__classcall__(cls, cartan_type, seq)

    def __init__(self, cartan_type, seq):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PolyhedralRealization(['A',2])
            sage: TestSuite(B).run() # long time
        """
        cat = (HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=cat)
        self._cartan_type = cartan_type
        self._seq = seq
        # These are the additional factors we add as necessary
        self._factors = tuple([ElementaryCrystal(cartan_type, i) for i in seq])
        # public for TensorProductOfCrystals
        self.crystals = self._factors
        self._tp = [C.module_generators[0] for C in self.crystals]
        self.module_generators = (self.element_class(self, self._tp),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.infinity.PolyhedralRealization(['A',2])
            Polyhedral realization of B(oo) of type ['A', 2] using (1, 2)
        """
        return "Polyhedral realization of B(oo) of type {} using {}".format(self._cartan_type, self._seq)

    def finite_tensor_product(self, k):
        """
        Return the finite tensor product of crystals of length ``k``
        by truncating ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PolyhedralRealization(['A',2])
            sage: B.finite_tensor_product(5)
            Full tensor product of the crystals
             [The 1-elementary crystal of type ['A', 2],
              The 2-elementary crystal of type ['A', 2],
              The 1-elementary crystal of type ['A', 2],
              The 2-elementary crystal of type ['A', 2],
              The 1-elementary crystal of type ['A', 2]]
        """
        N = len(self._factors)
        crystals = [self._factors[i % N] for i in range(k)]
        return TensorProductOfCrystals(*crystals)

    class Element(TensorProductOfCrystalsElement):
        r"""
        An element in the polyhedral realization of `B(\infty)`.
        """
        # For simplicity (and safety), we use the regular crystals implementation
        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: B = crystals.infinity.PolyhedralRealization(['A',2,1])
                sage: mg = B.module_generators[0]
                sage: [mg.epsilon(i) for i in B.index_set()]
                [0, 0, 0]
                sage: elt = mg.f(0)
                sage: [elt.epsilon(i) for i in B.index_set()]
                [1, 0, 0]
                sage: elt = mg.f_string([0,1,2])
                sage: [elt.epsilon(i) for i in B.index_set()]
                [0, 0, 1]
                sage: elt = mg.f_string([0,1,2,2])
                sage: [elt.epsilon(i) for i in B.index_set()]
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

                sage: B = crystals.infinity.PolyhedralRealization(['A',2,1])
                sage: mg = B.module_generators[0]
                sage: [mg.phi(i) for i in B.index_set()]
                [0, 0, 0]
                sage: elt = mg.f(0)
                sage: [elt.phi(i) for i in B.index_set()]
                [-1, 1, 1]
                sage: elt = mg.f_string([0,1])
                sage: [elt.phi(i) for i in B.index_set()]
                [-1, 0, 2]
                sage: elt = mg.f_string([0,1,2,2])
                sage: [elt.phi(i) for i in B.index_set()]
                [1, 1, 0]
            """
            P = self.parent().weight_lattice_realization()
            h = P.simple_coroots()
            omega = P(self.weight()).scalar(h[i])
            return self.epsilon(i) + omega

        def e(self, i):
            """
            Return the action of `e_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.infinity.PolyhedralRealization(['A',2])
                sage: mg = B.module_generators[0]
                sage: all(mg.e(i) is None for i in B.index_set())
                True
                sage: mg.f(1).e(1) == mg
                True
            """
            N = len(self) + 1
            pos = None
            for k in range(1, N):
                if all(self._sig(i,k) > self._sig(i,j) for j in range(1, k)) and \
                        all(self._sig(i,k) >= self._sig(i,j) for j in range(k+1, N)):
                    crystal = self[-k].e(i)
                    pos = k
                    break

            nf = len(self.parent()._factors)
            if pos is None or pos <= nf:
                return None

            l = list(self)
            l[-pos] = crystal
            if pos <= 2*nf and all(b._m == 0 for b in l[-2*nf:-nf]):
                return self.__class__(self.parent(), l[:-nf])
            return self.__class__(self.parent(), l)

        def f(self, i):
            """
            Return the action of `f_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.infinity.PolyhedralRealization(['A',2])
                sage: mg = B.module_generators[0]
                sage: mg.f(1)
                [-1, 0, 0, 0]
                sage: mg.f_string([1,2,2,1])
                [-1, -2, -1, 0, 0, 0]
            """
            N = len(self) + 1
            pos = None
            for k in range(1, N):
                if all(self._sig(i,k) >= self._sig(i,j) for j in range(1, k)) and \
                        all(self._sig(i,k) > self._sig(i,j) for j in range(k+1, N)):
                    crystal = self[-k].f(i)
                    pos = k
                    break

            nf = len(self.parent()._factors)
            if pos <= nf:
                l = list(self)
                l[-pos] = l[-pos].f(i)
                return self.__class__(self.parent(), l + self.parent()._tp)
            return self._set_index(-pos, crystal)

        def truncate(self, k=None):
            r"""
            Truncate ``self`` to have length ``k`` and return as an element
            in a (finite) tensor product of crystals.

            INPUT:

            - ``k`` -- (optional) the length of the truncation; if not
              specified, then returns one more than the current non-ground-state
              elements (i.e. the current list in ``self``)

            EXAMPLES::

                sage: B = crystals.infinity.PolyhedralRealization(['A',2])
                sage: mg = B.module_generators[0]
                sage: elt = mg.f_string([1,2,2,1]); elt
                [-1, -2, -1, 0, 0, 0]
                sage: t = elt.truncate(); t
                [-1, -2, -1, 0, 0, 0]
                sage: t.parent() is B.finite_tensor_product(6)
                True
                sage: elt.truncate(2)
                [-1, -2]
                sage: elt.truncate(10)
                [-1, -2, -1, 0, 0, 0, 0, 0, 0, 0]
            """
            if k is None:
                k = len(self)

            P = self.parent().finite_tensor_product(k)
            if k <= len(self):
                l = self[:k]
            else:
                l = list(self)
                N = len(self.parent()._tp)
                while len(l) < k:
                    i = len(l) % N
                    l.append(self.parent()._tp[i])
            return P(*l)

