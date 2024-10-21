# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Crystals of Modified Nakajima Monomials

AUTHORS:

- Arthur Lubovsky: Initial version

- Ben Salisbury: Initial version

Let `Y_{i,k}`, for `i \in I` and `k \in \ZZ`, be a commuting set of
variables, and let `\boldsymbol{1}` be a new variable which commutes with
each `Y_{i,k}`.  (Here, `I` represents the index set of a Cartan datum.)  One
may endow the structure of a crystal on the set `\widehat{\mathcal{M}}` of
monomials of the form

.. MATH::

    M = \prod_{(i,k) \in I\times \ZZ_{\ge0}} Y_{i,k}^{y_i(k)}\boldsymbol{1}.

Elements of `\widehat{\mathcal{M}}` are called  *modified Nakajima monomials*.
We will omit the `\boldsymbol{1}` from the end of a monomial if there exists
at least one `y_i(k) \neq 0`.  The crystal structure on this set is defined by

.. MATH::

    \begin{aligned}
    \mathrm{wt}(M) &= \sum_{i\in I} \Bigl( \sum_{k\ge 0} y_i(k) \Bigr) \Lambda_i, \\
    \varphi_i(M) &= \max\Bigl\{ \sum_{0\le j \le k} y_i(j) : k\ge 0 \Bigr\}, \\
    \varepsilon_i(M) &= \varphi_i(M) - \langle h_i, \mathrm{wt}(M) \rangle, \\
    k_f = k_f(M) &= \min\Bigl\{ k\ge 0 : \varphi_i(M) = \sum_{0\le j\le k} y_i(j) \Bigr\}, \\
    k_e = k_e(M) &= \max\Bigl\{ k\ge 0 : \varphi_i(M) = \sum_{0\le j\le k} y_i(j) \Bigr\},
    \end{aligned}

where `\{h_i : i \in I\}` and `\{\Lambda_i : i \in I \}` are the simple
coroots and fundamental weights, respectively.  With a chosen set of integers
`C = (c_{ij})_{i\neq j}` such that `c_{ij}+c_{ji} =1`, one defines

.. MATH::

    A_{i,k} = Y_{i,k} Y_{i,k+1} \prod_{j\neq i} Y_{j,k+c_{ji}}^{a_{ji}},

where `(a_{ij})` is a Cartan matrix.  Then

.. MATH::

    \begin{aligned}
    e_iM &= \begin{cases} 0 & \text{if } \varepsilon_i(M) = 0, \\
    A_{i,k_e}M & \text{if } \varepsilon_i(M) > 0, \end{cases} \\
    f_iM &= A_{i,k_f}^{-1} M.
    \end{aligned}

It is shown in [KKS2007]_ that the connected component of `\widehat{\mathcal{M}}`
containing the element `\boldsymbol{1}`, which we denote by
`\mathcal{M}(\infty)`, is crystal isomorphic to the crystal `B(\infty)`.

Let `\widetilde{\mathcal{M}}` be `\widehat{\mathcal{M}}` as a set, and with
crystal structure defined as on `\widehat{\mathcal{M}}` with the exception
that

.. MATH::

    f_iM = \begin{cases} 0 & \text{if } \varphi_i(M) = 0, \\
    A_{i,k_f}^{-1}M & \text{if } \varphi_i(M) > 0. \end{cases}

Then Kashiwara [Ka2003]_ showed that the connected component in
`\widetilde{\mathcal{M}}` containing a monomial `M` such that `e_iM = 0`, for
all `i \in I`, is crystal isomorphic to the irreducible highest weight
crystal `B(\mathrm{wt}(M))`.

WARNING:

    Monomial crystals depend on the choice of positive integers
    `C = (c_{ij})_{i\neq j}` satisfying the condition `c_{ij}+c_{ji}=1`.
    We have chosen such integers uniformly such that `c_{ij} = 1` if
    `i < j` and `c_{ij} = 0` if `i>j`.
"""

# *****************************************************************************
#  Copyright (C) 2013
#
#  Arthur Lubovsky (alubovsky at albany dot edu)
#  Ben Salisbury (ben dot salisbury at cmich dot edu)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from copy import copy
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_space import MatrixSpace


class NakajimaMonomial(Element):
    r"""
    An element of the monomial crystal.

    Monomials of the form `Y_{i_1,k_1}^{y_1} \cdots Y_{i_t,k_t}^{y_t}`,
    where `i_1, \dots, i_t` are elements of the index set, `k_1, \dots, k_t`
    are nonnegative integers, and `y_1, \dots, y_t` are integers.

    EXAMPLES::

        sage: M = crystals.infinity.NakajimaMonomials(['B',3,1])
        sage: mg = M.module_generators[0]
        sage: mg
        1
        sage: mg.f_string([1,3,2,0,1,2,3,0,0,1])
        Y(0,0)^-1 Y(0,1)^-1 Y(0,2)^-1 Y(0,3)^-1 Y(1,0)^-3
         Y(1,1)^-2 Y(1,2) Y(2,0)^3 Y(2,2) Y(3,0) Y(3,2)^-1

    An example using the `A` variables::

        sage: M = crystals.infinity.NakajimaMonomials("A3")
        sage: M.set_variables('A')
        sage: mg = M.module_generators[0]
        sage: mg.f_string([1,2,3,2,1])
        A(1,0)^-1 A(1,1)^-1 A(2,0)^-2 A(3,0)^-1
        sage: mg.f_string([3,2,1])
        A(1,2)^-1 A(2,1)^-1 A(3,0)^-1
        sage: M.set_variables('Y')
    """

    def __init__(self, parent, Y, A):
        r"""
        INPUT:

        - ``d`` -- dictionary of with pairs of the form ``{(i,k): y}``

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials("C5")
            sage: mg = M.module_generators[0]
            sage: TestSuite(mg).run()
        """
        self._Y = Y
        self._A = A
        Element.__init__(self, parent)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A',5,2])
            sage: x = M({(1,0):1, (2,2):-2, (0,5):10}); x
            Y(0,5)^10 Y(1,0) Y(2,2)^-2
            sage: M.set_variables('A')
            sage: x
            A(1,0)^-2 A(1,1)^-2 A(2,0)^-4 A(2,1)^-2 A(3,0)^-2
            sage: M.set_variables('Y')
        """
        return getattr(self, '_repr_' + self.parent()._variable)()

    def _repr_Y(self):
        r"""
        Return a string representation of ``self`` in the `Y` variables.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A',5,2])
            sage: M({(1,0):1,(2,2):-2,(0,5):10})
            Y(0,5)^10 Y(1,0) Y(2,2)^-2
        """
        if not self._Y:
            return "1"

        L = sorted(self._Y.items(), key=lambda x: (x[0][0], x[0][1]))
        exp = lambda e: "^{}".format(e) if e != 1 else ""
        return ' '.join("Y({},{})".format(mon[0][0], mon[0][1]) + exp(mon[1])
                        for mon in L)

    def _repr_A(self):
        r"""
        Return a string representation of ``self`` in the `A` variables.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['B',4,1])
            sage: m = M.module_generators[0].f_string([4,2,1])
            sage: m._repr_A()
            'A(1,1)^-1 A(2,0)^-1 A(4,0)^-1'
        """
        try:
            Y = {(i,0): c for i,c in self.parent().hw}
        except Exception:
            Y = {}

        if not Y and not self._A:
            return "1"

        L = sorted(Y.items(), key=lambda x: (x[0][0], x[0][1]))
        exp = lambda e: "^{}".format(e) if e != 1 else ""
        ret = ' '.join("Y({},{})".format(mon[0][0], mon[0][1]) + exp(mon[1])
                        for mon in L)
        if not self._A:
            return ret
        if Y:
            ret += ' '
        L = sorted(self._A.items(), key=lambda x: (x[0][0], x[0][1]))
        return ret + ' '.join("A({},{})".format(mon[0][0], mon[0][1]) + exp(mon[1])
                              for mon in L)

    def __hash__(self):
        r"""
        TESTS::

            sage: M = crystals.infinity.NakajimaMonomials(['C',5])
            sage: m1 = M.module_generators[0].f(1)
            sage: m2 = M.module_generators[0].f(2)
            sage: hash(m1) != hash(m2)
            True
        """
        return hash(frozenset(tuple(self._Y.items())))

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['C',5])
            sage: m1 = M.module_generators[0].f(1)
            sage: m2 = M.module_generators[0].f(2)
            sage: m1.__eq__(m2)
            False
            sage: m1.__eq__(m1)
            True
        """
        if isinstance(other, NakajimaMonomial):
            return self._Y == other._Y
        return self._Y == other

    def __ne__(self, other):
        r"""
        EXAMPLES::

            sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['A',2],La[1]+La[2])
            sage: m0 = M.module_generators[0]
            sage: m = M.module_generators[0].f(1).f(2).f(2).f(1)
            sage: m.__ne__(m0)
            True
            sage: m.__ne__(m)
            False
        """
        return not self == other

    def __lt__(self, other):
        r"""
        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['F',4])
            sage: mg = M.module_generators[0]
            sage: m = mg.f(4)
            sage: m.__lt__(mg)
            False
            sage: mg.__lt__(m)
            False
        """
        return False

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['G',2,1])
            sage: x = M.module_generators[0].f_string([1,0,2])
            sage: latex(x)
            Y_{0,0}^{-1} Y_{1,0}^{-1} Y_{1,1}^{2} Y_{2,0} Y_{2,1}^{-1}
            sage: M.set_variables('A')
            sage: latex(x)
            A_{0,0}^{-1} A_{1,0}^{-1} A_{2,0}^{-1}
            sage: M.set_variables('Y')
        """
        return getattr(self, '_latex_' + self.parent()._variable)()

    def _latex_Y(self):
        r"""
        Return a `\LaTeX` representation of ``self`` in the `Y` variables.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['G',2,1])
            sage: M.module_generators[0].f_string([1,0,2])._latex_Y()
            'Y_{0,0}^{-1} Y_{1,0}^{-1} Y_{1,1}^{2} Y_{2,0} Y_{2,1}^{-1} '
        """
        if not self._Y:
            return "\\boldsymbol{1}"

        L = sorted(self._Y.items(), key=lambda x:(x[0][0],x[0][1]))
        return_str = ''
        for x in L:
            if x[1] != 1:
                return_str += "Y_{%s,%s}" % (x[0][0],x[0][1]) + "^{%s} " % x[1]
            else:
                return_str += "Y_{%s,%s} " % (x[0][0],x[0][1])
        return return_str

    def _latex_A(self):
        r"""
        Return a `\LaTeX` representation of ``self`` in the `A` variables.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['C',4,1])
            sage: m = M.module_generators[0].f_string([4,2,3])
            sage: m._latex_A()
            'A_{2,0}^{-1} A_{3,1}^{-1} A_{4,0}^{-1} '
        """
        try:
            Y = {(i,0): c for i,c in self.parent().hw}
        except Exception:
            Y = {}

        if not Y and not self._A:
            return "\\boldsymbol{1}"

        L = sorted(Y.items(), key=lambda x:(x[0][0],x[0][1]))
        return_str = ''
        for x in L:
            if x[1] != 1:
                return_str += "Y_{%s,%s}" % (x[0][0],x[0][1]) + "^{%s} " % x[1]
            else:
                return_str += "Y_{%s,%s} " % (x[0][0],x[0][1])
        L = sorted(self._A.items(), key=lambda x:(x[0][0],x[0][1]))
        for x in L:
            if x[1] != 1:
                return_str += "A_{%s,%s}" % (x[0][0],x[0][1]) + "^{%s} " % x[1]
            else:
                return_str += "A_{%s,%s} " % (x[0][0],x[0][1])
        return return_str

    def _classical_weight(self):
        r"""
        Return the weight of ``self`` as an element of the classical version of
        ``self.parent().weight_lattice_realization``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,2])
            sage: m = M.module_generators[0].f_string([0,3,2,0,1])
            sage: m._classical_weight()
            -2*Lambda[0] + Lambda[1]

            sage: M = crystals.infinity.NakajimaMonomials(['E',6])
            sage: m = M.module_generators[0].f_string([1,5,2,6,3])
            sage: m._classical_weight()
            (-1/2, -3/2, 3/2, 1/2, -1/2, 1/2, 1/2, -1/2)
        """
        P = self.parent().weight_lattice_realization()
        La = P.fundamental_weights()
        return P(sum(v*La[k[0]] for k,v in self._Y.items()))

    def weight_in_root_lattice(self):
        r"""
        Return the weight of ``self`` as an element of the root lattice.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['F',4])
            sage: m = M.module_generators[0].f_string([3,3,1,2,4])
            sage: m.weight_in_root_lattice()
            -alpha[1] - alpha[2] - 2*alpha[3] - alpha[4]

            sage: M = crystals.infinity.NakajimaMonomials(['B',3,1])
            sage: mg = M.module_generators[0]
            sage: m = mg.f_string([1,3,2,0,1,2,3,0,0,1])
            sage: m.weight_in_root_lattice()
            -3*alpha[0] - 3*alpha[1] - 2*alpha[2] - 2*alpha[3]

            sage: M = crystals.infinity.NakajimaMonomials(['C',3,1])
            sage: m = M.module_generators[0].f_string([3,0,1,2,0])
            sage: m.weight_in_root_lattice()
            -2*alpha[0] - alpha[1] - alpha[2] - alpha[3]
        """
        Q = RootSystem(self.parent().cartan_type()).root_lattice()
        al = Q.simple_roots()
        return Q.sum(e*al[k[0]] for k,e in self._A.items())

    def weight(self):
        r"""
        Return the weight of ``self`` as an element of the weight lattice.

        EXAMPLES::

            sage: C = crystals.infinity.NakajimaMonomials(['A',1,1])
            sage: v = C.highest_weight_vector()
            sage: v.f(1).weight() + v.f(0).weight()
            -delta

            sage: M = crystals.infinity.NakajimaMonomials(['A',4,2])
            sage: m = M.highest_weight_vector().f_string([1,2,0,1])
            sage: m.weight()
            2*Lambda[0] - Lambda[1] - delta
        """
        P = self.parent().weight_lattice_realization()
        return P(self.weight_in_root_lattice())

    def epsilon(self, i):
        r"""
        Return the value of `\varepsilon_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['G',2,1])
            sage: m = M.module_generators[0].f(2)
            sage: [m.epsilon(i) for i in M.index_set()]
            [0, 0, 1]

            sage: M = crystals.infinity.NakajimaMonomials(['C',4,1])
            sage: m = M.module_generators[0].f_string([4,2,3])
            sage: [m.epsilon(i) for i in M.index_set()]
            [0, 0, 0, 1, 0]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        h = self.parent().weight_lattice_realization().simple_coroots()
        return self.phi(i) - self._classical_weight().scalar(h[i])

    def phi(self, i):
        r"""
        Return the value of `\varphi_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,3])
            sage: m = M.module_generators[0].f(1)
            sage: [m.phi(i) for i in M.index_set()]
            [1, -1, 1]

            sage: M = crystals.infinity.NakajimaMonomials(['C',4,1])
            sage: m = M.module_generators[0].f_string([4,2,3])
            sage: [m.phi(i) for i in M.index_set()]
            [0, 1, -1, 2, -1]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        if not self._Y or all(x[0] != i for x in self._Y):
            return ZZ.zero()

        d = copy(self._Y)
        K = max(x[1] for x in d if x[0] == i)
        for a in range(K):
            if (i,a) in d:
                continue
            else:
                d[(i,a)] = 0
        S = sorted((x for x in d.items() if x[0][0] == i), key=lambda x: x[0][1])
        return max(sum(S[k][1] for k in range(s)) for s in range(1,len(S)+1))

    def _ke(self, i):
        r"""
        Return the value `k_e` with respect to ``i`` and ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,3])
            sage: m = M.module_generators[0].f(1)
            sage: [m._ke(i) for i in M.index_set()]
            [+Infinity, 0, +Infinity]
        """
        h = self.parent().weight_lattice_realization().simple_coroots()
        phi = self.phi(i)
        if phi == self._classical_weight().scalar(h[i]): # self.epsilon(i) == 0
            return Infinity

        d = copy(self._Y)
        K = max(x[1] for x in d if x[0] == i)
        for a in range(K):
            if (i,a) in d:
                continue
            else:
                d[(i,a)] = 0
        total = ZZ.zero()
        L = []
        S = sorted((x for x in d.items() if x[0][0] == i), key=lambda x: x[0][1])
        for var,exp in S:
            total += exp
            if total == phi:
                L.append(var[1])

        return max(L) if L else ZZ.zero()

    def _kf(self, i):
        r"""
        Return the value `k_f` with respect to ``i`` and ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['F',4,1])
            sage: m = M.module_generators[0].f_string([0,1,4,3])
            sage: [m._kf(i) for i in M.index_set()]
            [0, 0, 2, 0, 0]
        """
        if all(i != x[0] for x in self._Y):
            return ZZ.zero()

        d = copy(self._Y)
        K = max(key[1] for key in d if key[0] == i)
        for a in range(K):
            if (i,a) in d:
                continue
            else:
                d[(i,a)] = 0
        S = sorted((x for x in d.items() if x[0][0] == i), key=lambda x: x[0][1])
        sum = 0
        phi = self.phi(i)
        for var,exp in S:
            sum += exp
            if sum == phi:
                return var[1]

    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['E',7,1])
            sage: m = M.module_generators[0].f_string([0,1,4,3])
            sage: [m.e(i) for i in M.index_set()]
            [None,
             None,
             None,
             Y(0,0)^-1 Y(1,1)^-1 Y(2,1) Y(3,0) Y(3,1) Y(4,0)^-1 Y(4,1)^-1 Y(5,0),
             None,
             None,
             None,
             None]

            sage: M = crystals.infinity.NakajimaMonomials("C5")
            sage: m = M.module_generators[0].f_string([1,3])
            sage: [m.e(i) for i in M.index_set()]
            [Y(2,1) Y(3,0)^-1 Y(3,1)^-1 Y(4,0),
             None,
             Y(1,0)^-1 Y(1,1)^-1 Y(2,0),
             None,
             None]

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,1])
            sage: M.set_variables('A')
            sage: m = M.module_generators[0].f_string([4,2,3,0])
            sage: [m.e(i) for i in M.index_set()]
            [A(2,1)^-1 A(3,1)^-1 A(4,0)^-1,
             None,
             None,
             A(0,2)^-1 A(2,1)^-1 A(4,0)^-1,
             None]
            sage: M.set_variables('Y')
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        if self.epsilon(i) == 0:
            return None

        newdict = copy(self._Y)
        ke = self._ke(i)
        Aik = {(i, ke): 1, (i, ke+1): 1}
        ct = self.parent().cartan_type()
        cm = ct.cartan_matrix()
        shift = 0
        if self.parent().cartan_type().is_finite():
            shift = 1
        for j_index,j in enumerate(self.parent().index_set()):
            if i == j:
                continue
            c = self.parent()._c[j_index,i-shift]
            if cm[j_index,i-shift] != 0:
                Aik[(j, ke+c)] = cm[j_index,i-shift]
        # Multiply by Aik
        for key,value in Aik.items():
            if key in newdict:
                if newdict[key] == -value: # The result would be a 0 exponent
                    del newdict[key]
                else:
                    newdict[key] += value
            else:
                newdict[key] = value
        A = copy(self._A)
        A[(i,ke)] = A.get((i,ke),0) + 1
        if not A[(i,ke)]:
            del A[(i,ke)]
        return self.__class__(self.parent(), newdict, A)

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials("B4")
            sage: m = M.module_generators[0].f_string([1,3,4])
            sage: [m.f(i) for i in M.index_set()]
            [Y(1,0)^-2 Y(1,1)^-2 Y(2,0)^2 Y(2,1) Y(3,0)^-1 Y(4,0) Y(4,1)^-1,
             Y(1,0)^-1 Y(1,1)^-1 Y(1,2) Y(2,0) Y(2,2)^-1 Y(3,0)^-1 Y(3,1) Y(4,0) Y(4,1)^-1,
             Y(1,0)^-1 Y(1,1)^-1 Y(2,0) Y(2,1)^2 Y(3,0)^-2 Y(3,1)^-1 Y(4,0)^3 Y(4,1)^-1,
             Y(1,0)^-1 Y(1,1)^-1 Y(2,0) Y(2,1) Y(3,0)^-1 Y(3,1) Y(4,1)^-2]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        newdict = copy(self._Y)
        kf = self._kf(i)
        Aik = {(i, kf): -1, (i, kf+1): -1}
        ct = self.parent().cartan_type()
        cm = ct.cartan_matrix()
        shift = 0
        if ct.is_finite():
            shift = 1
        for j_index,j in enumerate(self.parent().index_set()):
            if i == j:
                continue
            c = self.parent()._c[j_index,i-shift]
            if cm[j_index,i-shift] != 0:
                Aik[(j, kf+c)] = -cm[j_index,i-shift]
        # Multiply by Aik
        for key,value in Aik.items():
            if key in newdict:
                if newdict[key] == -value: # The result would be a 0 exponent
                    del newdict[key]
                else:
                    newdict[key] += value
            else:
                newdict[key] = value
        A = copy(self._A)
        A[(i,kf)] = A.get((i,kf),0) - 1
        if not A[(i,kf)]:
            del A[(i,kf)]
        return self.__class__(self.parent(), newdict, A)


class InfinityCrystalOfNakajimaMonomials(UniqueRepresentation, Parent):
    r"""
    Crystal `B(\infty)` in terms of (modified) Nakajima monomials.

    Let `Y_{i,k}`, for `i \in I` and `k \in \ZZ`, be a commuting set of
    variables, and let `\boldsymbol{1}` be a new variable which commutes
    with each `Y_{i,k}`.  (Here, `I` represents the index set of a Cartan
    datum.)  One may endow the structure of a crystal on the
    set `\widehat{\mathcal{M}}` of monomials of the form

    .. MATH::

        M = \prod_{(i,k) \in I\times \ZZ_{\ge0}} Y_{i,k}^{y_i(k)}\boldsymbol{1}.

    Elements of `\widehat{\mathcal{M}}` are called
    *modified Nakajima monomials*. We will omit the `\boldsymbol{1}`
    from the end of a monomial if there exists at least one `y_i(k) \neq 0`.
    The crystal structure on this set is defined by

    .. MATH::

        \begin{aligned}
        \mathrm{wt}(M) & = \sum_{i\in I} \Bigl( \sum_{k \ge 0}
        y_i(k) \Bigr) \Lambda_i, \\
        \varphi_i(M) & = \max\Bigl\{ \sum_{0 \le j \le k} y_i(j) :
        k \ge 0 \Bigr\}, \\
        \varepsilon_i(M) & = \varphi_i(M) -
        \langle h_i, \mathrm{wt}(M) \rangle, \\
        k_f = k_f(M) & = \min\Bigl\{ k \ge 0 :
        \varphi_i(M) = \sum_{0 \le j \le k} y_i(j) \Bigr\}, \\
        k_e = k_e(M) & = \max\Bigl\{ k \ge 0 :
        \varphi_i(M) = \sum_{0 \le j \le k} y_i(j) \Bigr\},
        \end{aligned}

    where `\{h_i : i \in I\}` and `\{\Lambda_i : i \in I \}` are the simple
    coroots and fundamental weights, respectively.  With a chosen set of
    nonnegative integers `C = (c_{ij})_{i\neq j}` such that
    `c_{ij} + c_{ji} = 1`, one defines

    .. MATH::

        A_{i,k} = Y_{i,k} Y_{i,k+1} \prod_{j\neq i} Y_{j,k+c_{ji}}^{a_{ji}},

    where `(a_{ij})_{i,j \in I}` is a Cartan matrix.  Then

    .. MATH::

        \begin{aligned}
        e_iM &= \begin{cases} 0 & \text{if } \varepsilon_i(M) = 0, \\
        A_{i,k_e}M & \text{if } \varepsilon_i(M) > 0, \end{cases} \\
        f_iM &= A_{i,k_f}^{-1} M.
        \end{aligned}

    It is shown in [KKS2007]_ that the connected component of
    `\widehat{\mathcal{M}}` containing the element `\boldsymbol{1}`,
    which we denote by `\mathcal{M}(\infty)`, is crystal isomorphic
    to the crystal `B(\infty)`.

    INPUT:

    - ``cartan_type`` -- a Cartan type

    - ``c`` -- (optional) the matrix `(c_{ij})_{i,j \in I}` such that
      `c_{ii} = 0` for all `i \in I`, `c_{ij} \in \ZZ_{>0}` for all
      `i,j \in I`, and `c_{ij} + c_{ji} = 1` for all `i \neq j`; the
      default is `c_{ij} = 0` if `i < j` and `0` otherwise

    EXAMPLES::

        sage: B = crystals.infinity.Tableaux("C3")
        sage: S = B.subcrystal(max_depth=4)
        sage: G = B.digraph(subset=S) # long time
        sage: M = crystals.infinity.NakajimaMonomials("C3") # long time
        sage: T = M.subcrystal(max_depth=4) # long time
        sage: H = M.digraph(subset=T) # long time
        sage: G.is_isomorphic(H,edge_labels=True) # long time
        True

        sage: M = crystals.infinity.NakajimaMonomials(['A',2,1])
        sage: T = M.subcrystal(max_depth=3)
        sage: H = M.digraph(subset=T) # long time
        sage: Y = crystals.infinity.GeneralizedYoungWalls(2)
        sage: YS = Y.subcrystal(max_depth=3)
        sage: YG = Y.digraph(subset=YS) # long time
        sage: YG.is_isomorphic(H,edge_labels=True) # long time
        True

        sage: M = crystals.infinity.NakajimaMonomials("D4")
        sage: B = crystals.infinity.Tableaux("D4")
        sage: MS = M.subcrystal(max_depth=3)
        sage: BS = B.subcrystal(max_depth=3)
        sage: MG = M.digraph(subset=MS) # long time
        sage: BG = B.digraph(subset=BS) # long time
        sage: BG.is_isomorphic(MG,edge_labels=True) # long time
        True
    """
    @staticmethod
    def _normalize_c(c, n):
        """
        Normalize the input ``c``.

        EXAMPLES::

            sage: from sage.combinat.crystals.monomial_crystals import InfinityCrystalOfNakajimaMonomials
            sage: InfinityCrystalOfNakajimaMonomials._normalize_c(None, 4)
            [0 1 1 1]
            [0 0 1 1]
            [0 0 0 1]
            [0 0 0 0]
            sage: c = matrix([[0,1,1],[0,0,0],[0,1,0]]); c
            [0 1 1]
            [0 0 0]
            [0 1 0]
            sage: c.is_mutable()
            True
            sage: C = InfinityCrystalOfNakajimaMonomials._normalize_c(c, 3); C
            [0 1 1]
            [0 0 0]
            [0 1 0]
            sage: C.is_mutable()
            False

        TESTS::

            sage: c = matrix([[0,1],[0,1]])
            sage: C = InfinityCrystalOfNakajimaMonomials._normalize_c(c, 2)
            Traceback (most recent call last):
            ...
            ValueError: the c matrix must have 0s on the diagonal
            sage: c = matrix([[0,2],[-1,0]])
            sage: C = InfinityCrystalOfNakajimaMonomials._normalize_c(c, 2)
            Traceback (most recent call last):
            ...
            ValueError: the c matrix must have nonnegative entries
            sage: c = matrix([[0,1],[1,0]])
            sage: C = InfinityCrystalOfNakajimaMonomials._normalize_c(c, 2)
            Traceback (most recent call last):
            ...
            ValueError: transpose entries do not sum to 1
        """
        if c is None:
            # Default is i < j <=> c_{ij} = 1 (0 otherwise)
            c = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
        MS = MatrixSpace(ZZ, n, n)
        c = MS(c)
        c.set_immutable()
        if any(c[i,i] != 0 for i in range(n)):
            raise ValueError("the c matrix must have 0s on the diagonal")
        if any(c[i,j] + c[j,i] != 1 for i in range(n) for j in range(i)):
            raise ValueError("transpose entries do not sum to 1")
        if any(c[i,j] < 0 or c[j,i] < 0 for i in range(n) for j in range(i)):
            raise ValueError("the c matrix must have nonnegative entries")
        return c

    @staticmethod
    def __classcall_private__(cls, ct, c=None):
        r"""
        Normalize input to ensure a unique representation.

        INPUT:

        - ``ct`` -- a Cartan type

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials("E8")
            sage: M1 = crystals.infinity.NakajimaMonomials(['E',8])
            sage: M2 = crystals.infinity.NakajimaMonomials(CartanType(['E',8]))
            sage: M is M1 is M2
            True
        """
        cartan_type = CartanType(ct)
        n = len(cartan_type.index_set())
        c = InfinityCrystalOfNakajimaMonomials._normalize_c(c, n)
        M = super().__classcall__(cls, cartan_type, c)
        M.set_variables('Y')
        return M

    def __init__(self, ct, c, category=None):
        r"""
        EXAMPLES::

            sage: Minf = crystals.infinity.NakajimaMonomials(['A',3])
            sage: TestSuite(Minf).run() # long time
        """
        self._cartan_type = ct
        self._c = c
        self._variable = 'Y'

        if category is None:
            category = (HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=category)
        self.module_generators = (self.element_class(self, {}, {}),)

    def _element_constructor_(self, Y=None, A=None):
        r"""
        Construct an element of ``self`` from ``Y``.

        INPUT:

        - ``Y`` -- dictionary whose key is a pair and whose value
          is an integer
        - ``A`` -- dictionary whose key is a pair and whose value
          is an integer

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,1])
            sage: m = M({(1,0):-1,(1,1):-1,(2,0):1})
            sage: m
            Y(1,0)^-1 Y(1,1)^-1 Y(2,0)

            sage: M = crystals.infinity.NakajimaMonomials(['A',2,1])
            sage: m = M(A={(0,1): -1, (1,1): -2, (2,0): -1, (2,1): -1})
            sage: m._repr_A()
            'A(0,1)^-1 A(1,1)^-2 A(2,0)^-1 A(2,1)^-1'
            sage: m
            Y(0,2)^2 Y(1,2)^-1 Y(2,0)^-1 Y(2,1) Y(2,2)^-1
            sage: m == M.highest_weight_vector().f_string([2,0,1,2,1])
            True
        """
        if A is None:
            if Y is None:
                return self.module_generators[0]
            # This is a crude way to determine the A, but it works
            hw,path = self.element_class(self, Y, {}).to_highest_weight()
            hw._A = {}
            return hw.f_string(reversed(path))
        elif Y is None or Y == 0:
            # The Y == 0 check is because the parent's __call__ has that
            #   as the first default value
            ct = self.cartan_type()
            cm = ct.cartan_matrix()
            I = self.index_set()
            shift = 0
            if ct.is_finite():
                shift = 1
            Y = {}
            for k,v in A.items():
                Y[k] = Y.get(k, 0) + v
                Y[(k[0],k[1]+1)] = Y.get((k[0],k[1]+1), 0) + v
                for j_index,j in enumerate(I):
                    if k[0] == j:
                        continue
                    c = self._c[j_index,k[0]-shift]
                    if cm[j_index,k[0]-shift] != 0:
                        Y[(j,k[1]+c)] = Y.get((j,k[1]+c), 0) + v*cm[j_index,k[0]-shift]
            for k in list(Y):
                if Y[k] == 0:
                    del Y[k]
        return self.element_class(self, Y, A)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,1])
            sage: m = M({(1,0):-1,(1,1):-1,(2,0):1})
            sage: m
            Y(1,0)^-1 Y(1,1)^-1 Y(2,0)
        """
        return "Infinity Crystal of modified Nakajima monomials of type {}".format(self._cartan_type)

    def c(self):
        """
        Return the matrix `c_{ij}` of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['B',3]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(La[1]+La[2])
            sage: M.c()
            [0 1 1]
            [0 0 1]
            [0 0 0]

            sage: c = Matrix([[0,0,1],[1,0,0],[0,1,0]])
            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: M = crystals.NakajimaMonomials(2*La[1], c=c)
            sage: M.c() == c
            True
        """
        return self._c

    def cardinality(self):
        r"""
        Return the cardinality of ``self``, which is always `\infty`.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A',5,2])
            sage: M.cardinality()
            +Infinity
        """
        return Infinity

    def set_variables(self, letter):
        r"""
        Set the type of monomials to use for the element output.

        If the `A` variables are used, the output is written as
        `\prod_{i\in I} Y_{i,0}^{\lambda_i} \prod_{i,k} A_{i,k}^{c_{i,k}}`, where
        `\sum_{i \in I} \lambda_i \Lambda_i` is the corresponding
        dominant weight.

        INPUT:

        - ``letter`` -- can be one of the following:

          * ``'Y'`` -- use `Y_{i,k}`, corresponds to fundamental weights
          * ``'A'`` -- use `A_{i,k}`, corresponds to simple roots

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A', 4])
            sage: elt = M.highest_weight_vector().f_string([2,1,3,2,3,2,4,3])
            sage: elt
            Y(1,2) Y(2,0)^-1 Y(2,2)^-1 Y(3,0)^-1 Y(3,2)^-1 Y(4,0)
            sage: M.set_variables('A')
            sage: elt
            A(1,1)^-1 A(2,0)^-1 A(2,1)^-2 A(3,0)^-2 A(3,1)^-1 A(4,0)^-1
            sage: M.set_variables('Y')

        ::

            sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(La[1]+La[2])
            sage: lw = M.lowest_weight_vectors()[0]
            sage: lw
            Y(1,2)^-1 Y(2,1)^-1
            sage: M.set_variables('A')
            sage: lw
            Y(1,0) Y(2,0) A(1,0)^-1 A(1,1)^-1 A(2,0)^-2
            sage: M.set_variables('Y')
        """
        if letter not in ['Y', 'A']:
            raise ValueError("invalid monomial type")
        self._variable = letter

    def get_variables(self):
        """
        Return the type of monomials to use for the element output.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A', 4])
            sage: M.get_variables()
            'Y'
        """
        return self._variable

    Element = NakajimaMonomial


class CrystalOfNakajimaMonomialsElement(NakajimaMonomial):
    r"""
    Element class for
    :class:`~sage.combinat.crystals.monomial_crystals.CrystalOfNakajimaMonomials`.

    The `f_i` operators need to be modified from the version in
    :class:`~sage.combinat.crystals.monomial_crystalsNakajimaMonomial`
    in order to create irreducible highest weight realizations.
    This modified `f_i` is defined as

    .. MATH::

        f_iM = \begin{cases} 0 & \text{if } \varphi_i(M) = 0, \\
        A_{i,k_f}^{-1}M & \text{if } \varphi_i(M) > 0. \end{cases}

    EXAMPLES::

        sage: La = RootSystem(['A',5,2]).weight_lattice(extended=True).fundamental_weights()
        sage: M = crystals.NakajimaMonomials(['A',5,2],3*La[0])
        sage: m = M.module_generators[0].f(0); m
        Y(0,0)^2 Y(0,1)^-1 Y(2,0)
        sage: TestSuite(m).run()
    """

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: La = RootSystem(['A',5,2]).weight_lattice(extended=True).fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['A',5,2],3*La[0])
            sage: m = M.module_generators[0]
            sage: [m.f(i) for i in M.index_set()]
            [Y(0,0)^2 Y(0,1)^-1 Y(2,0), None, None, None]

        ::

            sage: M = crystals.infinity.NakajimaMonomials("E8")
            sage: M.set_variables('A')
            sage: m = M.module_generators[0].f_string([4,2,3,8])
            sage: m
            A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1
            sage: [m.f(i) for i in M.index_set()]
            [A(1,2)^-1 A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1,
             A(2,0)^-1 A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1,
             A(2,1)^-1 A(3,0)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(4,1)^-1 A(8,0)^-1,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(5,0)^-1 A(8,0)^-1,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(6,0)^-1 A(8,0)^-1,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(7,1)^-1 A(8,0)^-1,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-2]
            sage: M.set_variables('Y')
        """
        if self.phi(i) == 0:
            return None
        return super().f(i)

    def weight(self):
        r"""
        Return the weight of ``self`` as an element of the weight lattice.

        EXAMPLES::

            sage: La = RootSystem("A2").weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials("A2",La[1]+La[2])
            sage: M.module_generators[0].weight()
            (2, 1, 0)
        """
        P = self.parent().weight_lattice_realization()
        return P(self.weight_in_root_lattice()) + P(self.parent().hw)


class CrystalOfNakajimaMonomials(InfinityCrystalOfNakajimaMonomials):
    r"""
    Let `\widetilde{\mathcal{M}}` be `\widehat{\mathcal{M}}` as a set, and with
    crystal structure defined as on `\widehat{\mathcal{M}}` with the exception
    that

    .. MATH::

        f_iM = \begin{cases} 0 & \text{if } \varphi_i(M) = 0, \\
        A_{i,k_f}^{-1}M & \text{if } \varphi_i(M) > 0. \end{cases}

    Then Kashiwara [Ka2003]_ showed that the connected component in
    `\widetilde{\mathcal{M}}` containing a monomial `M` such that `e_iM = 0`,
    for all `i \in I`, is crystal isomorphic to the irreducible highest weight
    crystal `B(\mathrm{wt}(M))`.

    INPUT:

    - ``ct`` -- a Cartan type

    - ``La`` -- an element of the weight lattice

    EXAMPLES::

        sage: La = RootSystem("A2").weight_lattice().fundamental_weights()
        sage: M = crystals.NakajimaMonomials("A2",La[1]+La[2])
        sage: B = crystals.Tableaux("A2",shape=[2,1])
        sage: GM = M.digraph()
        sage: GB = B.digraph()
        sage: GM.is_isomorphic(GB,edge_labels=True)
        True

        sage: La = RootSystem("G2").weight_lattice().fundamental_weights()
        sage: M = crystals.NakajimaMonomials("G2",La[1]+La[2])
        sage: B = crystals.Tableaux("G2",shape=[2,1])
        sage: GM = M.digraph()
        sage: GB = B.digraph()
        sage: GM.is_isomorphic(GB,edge_labels=True)
        True

        sage: La = RootSystem("B2").weight_lattice().fundamental_weights()
        sage: M = crystals.NakajimaMonomials(['B',2],La[1]+La[2])
        sage: B = crystals.Tableaux("B2",shape=[3/2,1/2])
        sage: GM = M.digraph()
        sage: GB = B.digraph()
        sage: GM.is_isomorphic(GB,edge_labels=True)
        True

        sage: La = RootSystem(['A',3,1]).weight_lattice(extended=True).fundamental_weights()
        sage: M = crystals.NakajimaMonomials(['A',3,1],La[0]+La[2])
        sage: B = crystals.GeneralizedYoungWalls(3,La[0]+La[2])
        sage: SM = M.subcrystal(max_depth=4)
        sage: SB = B.subcrystal(max_depth=4)
        sage: GM = M.digraph(subset=SM) # long time
        sage: GB = B.digraph(subset=SB) # long time
        sage: GM.is_isomorphic(GB,edge_labels=True) # long time
        True

        sage: La = RootSystem(['A',5,2]).weight_lattice(extended=True).fundamental_weights()
        sage: LA = RootSystem(['A',5,2]).weight_space().fundamental_weights()
        sage: M = crystals.NakajimaMonomials(['A',5,2],3*La[0])
        sage: B = crystals.LSPaths(3*LA[0])
        sage: SM = M.subcrystal(max_depth=4)
        sage: SB = B.subcrystal(max_depth=4)
        sage: GM = M.digraph(subset=SM)
        sage: GB = B.digraph(subset=SB)
        sage: GM.is_isomorphic(GB,edge_labels=True)
        True

        sage: c = matrix([[0,1,0],[0,0,1],[1,0,0]])
        sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
        sage: M = crystals.NakajimaMonomials(2*La[1], c=c)
        sage: sorted(M.subcrystal(max_depth=3), key=str)
        [Y(0,0) Y(0,1) Y(1,0) Y(2,1)^-1,
         Y(0,0) Y(0,1)^2 Y(1,1)^-1 Y(2,0) Y(2,1)^-1,
         Y(0,0) Y(0,2)^-1 Y(1,0) Y(1,1) Y(2,1)^-1 Y(2,2),
         Y(0,1) Y(0,2)^-1 Y(1,1)^-1 Y(2,0)^2 Y(2,2),
         Y(0,1) Y(1,0) Y(1,1)^-1 Y(2,0),
         Y(0,1)^2 Y(1,1)^-2 Y(2,0)^2,
         Y(0,2)^-1 Y(1,0) Y(2,0) Y(2,2),
         Y(1,0) Y(1,3) Y(2,0) Y(2,3)^-1,
         Y(1,0)^2]
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, La=None, c=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: La = RootSystem(['E',8,1]).weight_lattice(extended=True).fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['E',8,1],La[0]+La[8])
            sage: M1 = crystals.NakajimaMonomials(CartanType(['E',8,1]),La[0]+La[8])
            sage: M2 = crystals.NakajimaMonomials(['E',8,1],M.Lambda()[0] + M.Lambda()[8])
            sage: M is M1 is M2
            True
        """
        if La is None:
            La = cartan_type
            cartan_type = La.parent().cartan_type()
        cartan_type = CartanType(cartan_type)
        if cartan_type.is_affine():
            La = RootSystem(cartan_type).weight_lattice(extended=True)(La)
        else:
            La = RootSystem(cartan_type).weight_lattice()(La)
        n = len(cartan_type.index_set())
        c = InfinityCrystalOfNakajimaMonomials._normalize_c(c, n)
        return super().__classcall__(cls, cartan_type, La, c)

    def __init__(self, ct, La, c):
        r"""
        EXAMPLES::

            sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['A',2], La[1]+La[2])
            sage: TestSuite(M).run()

            sage: La = RootSystem(['C',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['C',2,1], La[0])
            sage: TestSuite(M).run(max_runs=100)
        """
        if ct.is_finite():
            cat = ClassicalCrystals()
        else:
            cat = (RegularCrystals(), HighestWeightCrystals(), InfiniteEnumeratedSets())
        InfinityCrystalOfNakajimaMonomials.__init__(self, ct, c, cat)
        self._cartan_type = ct
        self.hw = La
        gen = {(i,0): c for i,c in La}
        self.module_generators = (self.element_class(self, gen, {}),)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['C',3,1]).weight_lattice(extended=True).fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['C',3,1],La[0]+5*La[3])
            sage: M
            Highest weight crystal of modified Nakajima monomials of Cartan type ['C', 3, 1] and highest weight Lambda[0] + 5*Lambda[3]
        """
        return "Highest weight crystal of modified Nakajima monomials of Cartan type {1!s} and highest weight {0!s}".format(self.hw, self._cartan_type)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['A',2], La[1])
            sage: M.cardinality()
            3

            sage: La = RootSystem(['D',4,2]).weight_lattice(extended=True).fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['D',4,2], La[1])
            sage: M.cardinality()
            +Infinity
        """
        if not self.cartan_type().is_finite():
            return Infinity
        return super(InfinityCrystalOfNakajimaMonomials, self).cardinality()

    Element = CrystalOfNakajimaMonomialsElement
