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
`C = (c_{ij})_{i\neq j}` such that `c_{ij}+c{ji} =1`, one defines

.. MATH::

    A_{i,k} = Y_{i,k} Y_{i,k+1} \prod_{j\neq i} Y_{j,k+c_{ji}}^{a_{ji}},

where `(a_{ij})` is a Cartan matrix.  Then

.. MATH::

    \begin{aligned}
    e_iM &= \begin{cases} 0 & \text{if } \varepsilon_i(M) = 0, \\
    A_{i,k_e}M & \text{if } \varepsilon_i(M) > 0, \end{cases} \\
    f_iM &= A_{i,k_f}^{-1} M.
    \end{aligned}

It is shown in [KKS07]_ that the connected component of `\widehat{\mathcal{M}}`
containing the element `\boldsymbol{1}`, which we denote by
`\mathcal{M}(\infty)`, is crystal isomorphic to the crystal `B(\infty)`.

Let `\widetilde{\mathcal{M}}` be `\widehat{\mathcal{M}}` as a set, and with
crystal structure defined as on `\widehat{\mathcal{M}}` with the exception
that

.. MATH::

    f_iM = \begin{cases} 0 & \text{if } \varphi_i(M) = 0, \\
    A_{i,k_f}^{-1}M & \text{if } \varphi_i(M) > 0. \end{cases}

Then Kashiwara [Kash03]_ showed that the connected component in
`\widetilde{\mathcal{M}}` containing a monomial `M` such that `e_iM = 0`, for
all `i \in I`, is crystal isomorphic to the irreducible highest weight
crystal `B(\mathrm{wt}(M))`.

WARNING:

    Monomial crystals depend on the choice of positive integers
    `C = (c_{ij})_{i\neq j}` satisfying the condition `c_{ij}+c_{ji}=1`.
    We have chosen such integers uniformly such that `c_{ij} = 1` if
    `i < j` and `c_{ij} = 0` if `i>j`.

REFERENCES:

.. [KKS07] S.-J. Kang, J.-A. Kim, and D.-U. Shin.
   Modified Nakajima Monomials and the Crystal `B(\infty)`.
   J. Algebra **308**, pp. 524--535, 2007.

.. [Kash03] M. Kashiwara.
   Realizations of Crystals.
   Combinatorial and geometric representation theory (Seoul, 2001),
   Contemp. Math. **325**, Amer. Math. Soc., pp. 133--139, 2003.
"""

#******************************************************************************
#  Copyright (C) 2013
#
#  Arthur Lubovsky (alubovsky at albany dot edu)
#  Ben Salisbury (ben dot salisbury at cmich dot edu)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from copy import copy
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.rings.infinity import Infinity

class NakajimaYMonomial(Element):
    r"""
    Monomials of the form `Y_{i_1,k_1}^{a_1}\cdots Y_{i_t,k_t}^{y_t}`, where
    `i_1,\dots,i_t` are elements of the index set, `k_1,\dots,k_t` are
    nonnegative integers, and `y_1,\dots,y_t` are integers.

    EXAMPLES::

        sage: M = crystals.infinity.NakajimaMonomials(['B',3,1])
        sage: mg = M.module_generators[0]
        sage: mg
        1
        sage: mg.f_string([1,3,2,0,1,2,3,0,0,1])
        Y(0,0)^-1 Y(0,1)^-1 Y(0,2)^-1 Y(0,3)^-1 Y(1,0)^-3 Y(1,1)^-2 Y(1,2) Y(2,0)^3 Y(2,2) Y(3,0) Y(3,2)^-1
    """

    def __init__(self,parent,dict):
        r"""
        INPUT:

        - ``dict`` -- a dictionary of with pairs of the form ``{(i,k):y}``

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials("C5")
            sage: mg = M.module_generators[0]
            sage: TestSuite(mg).run()
        """
        self._dict = dict
        Element.__init__(self, parent)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A',5,2])
            sage: M({(1,0):1,(2,2):-2,(0,5):10})
            Y(0,5)^10 Y(1,0) Y(2,2)^-2
        """
        if self._dict == {}:
            return "1"
        else:
            L = sorted(self._dict.iteritems(), key=lambda x:(x[0][0],x[0][1]))
            return_str = ''
            for x in range(len(L)):
                if L[x][1] != 1:
                    return_str += "Y(%s,%s)"%(L[x][0][0],L[x][0][1]) + "^%s "%L[x][1]
                else:
                    return_str += "Y(%s,%s) "%(L[x][0][0],L[x][0][1])
            return return_str

    def __eq__(self,other):
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
        if isinstance(other, NakajimaYMonomial):
            return self._dict == other._dict
        return self._dict == other

    def __ne__(self,other):
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
        return not self.__eq__(other)

    def __lt__(self,other):
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
            sage: M.module_generators[0].f_string([1,0,2])._latex_()
            'Y_{0,0}^{-1} Y_{1,0}^{-1} Y_{1,1}^{2} Y_{2,0} Y_{2,1}^{-1} '
        """
        if self._dict == {}:
            return "\\boldsymbol{1}"
        else:
            L = sorted(self._dict.iteritems(), key=lambda x:(x[0][0],x[0][1]))
            return_str = ''
            for x in range(len(L)):
                if L[x][1] != 1:
                    return_str += "Y_{%s,%s}"%(L[x][0][0],L[x][0][1]) + "^{%s} "%L[x][1]
                else:
                    return_str += "Y_{%s,%s} "%(L[x][0][0],L[x][0][1])
            return return_str

    def weight(self):
        r"""
        Return the weight of ``self`` as an element of
        ``self.parent().weight_lattice_realization``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,2])
            sage: m = M.module_generators[0].f_string([0,3,2,0,1])
            sage: m.weight()
            -2*Lambda[0] + Lambda[1]

            sage: M = crystals.infinity.NakajimaMonomials(['E',6])
            sage: m = M.module_generators[0].f_string([1,5,2,6,3])
            sage: m.weight()
            (-1/2, -3/2, 3/2, 1/2, -1/2, 1/2, 1/2, -1/2)
        """
        P = self.parent().weight_lattice_realization()
        La = P.fundamental_weights()
        return P(sum(v*La[k[0]] for k,v in self._dict.iteritems()))

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
        """
        Q = RootSystem(self.parent().cartan_type()).root_lattice()
        alpha = Q.simple_roots()
        path = self.to_highest_weight()
        return Q(sum(-alpha[j] for j in path[1]))

    def epsilon(self,i):
        r"""
        Return the value of `\varepsilon_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['G',2,1])
            sage: m = M.module_generators[0].f(2)
            sage: [m.epsilon(i) for i in M.index_set()]
            [0, 0, 1]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        h = self.parent().weight_lattice_realization().simple_coroots()
        return self.phi(i) - self.weight().scalar(h[i])

    def phi(self,i):
        r"""
        Return the value of `\varphi_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,3])
            sage: m = M.module_generators[0].f(1)
            sage: [m.phi(i) for i in M.index_set()]
            [1, -1, 1]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        dict = self._dict
        if dict == {}:
            return 0
        else:
            L = [x[0] for x in dict.keys()]
            if i not in L:
                return 0
            else:
                d = copy(dict)
                K = max(x[1] for x in list(d) if x[0] ==i)
                for a in range(K):
                    if (i,a) in d:
                        continue
                    else:
                        d[(i,a)] = 0
                S = sorted((x for x in d.iteritems() if x[0][0]==i), key=lambda x: x[0][1])
                return max(sum(S[k][1] for k in range(s)) for s in range(1,len(S)+1))

    def _ke(self,i):
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
        dict = self._dict
        sum = 0
        L = []
        phi = self.phi(i)
        if self.epsilon(i) == 0:
            return Infinity
        else:
            d = copy(dict)
            K = max(x[1] for x in list(d) if x[0] ==i)
            for a in range(K):
                if (i,a) in d:
                    continue
                else:
                    d[(i,a)] = 0
            S = sorted((x for x in d.iteritems() if x[0][0]==i), key=lambda x: x[0][1])
            for var,exp in S:
                sum += exp
                if sum == phi:
                    L.append(var[1])
            if L == []:
                return 0
            else:
                return max(L)

    def _kf(self,i):
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
        d = copy(self._dict)
        I = [x[0] for x in d]
        if i not in I:
            return 0
        else:
            K = max(x[1] for x in list(d) if x[0] ==i)
            for a in range(K):
                if (i,a) in d:
                    continue
                else:
                    d[(i,a)] = 0
            S = sorted((x for x in d.iteritems() if x[0][0]==i), key=lambda x: x[0][1])
            sum = 0
            phi = self.phi(i)
            for var,exp in S:
                sum += exp
                if sum == phi:
                    return var[1]

    def e(self,i):
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
             Y(0,0)^-1 Y(1,1)^-1 Y(2,1) Y(3,0) Y(3,1) Y(4,0)^-1 Y(4,1)^-1 Y(5,0) ,
             None,
             None,
             None,
             None]

            sage: M = crystals.infinity.NakajimaMonomials("C5")
            sage: m = M.module_generators[0].f_string([1,3])
            sage: [m.e(i) for i in M.index_set()]
            [Y(2,1) Y(3,0)^-1 Y(3,1)^-1 Y(4,0) ,
             None,
             Y(1,0)^-1 Y(1,1)^-1 Y(2,0) ,
             None,
             None]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        if self.epsilon(i) == 0:
            return None

        newdict = copy(self._dict)
        ke = self._ke(i)
        Aik = {(i, ke):1, (i, ke+1):1}
        ct = self.parent().cartan_type()
        cm = ct.cartan_matrix()
        shift = 0
        if self.parent().cartan_type().is_finite():
            shift = 1
        for j in self.parent().index_set():
            if i == j:
                continue
            c = 0
            if i > j:
                c = 1
            if ct.is_affine() and ct.type() == 'A' and abs(i-j) == ct.rank() - 1:
                c = 1 - c
            if cm[j-shift][i-shift] != 0:
                Aik[(j, ke+c)] = cm[j-shift][i-shift]
        for key,value in Aik.iteritems():
            if key in newdict:
                newdict[key] +=value
            else:
                newdict[key] = value
        for k in list(newdict):
            if newdict[k] == 0:
                newdict.pop(k)
        return self.__class__(self.parent(),newdict)

    def f(self,i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials("B4")
            sage: m = M.module_generators[0].f_string([1,3,4])
            sage: [m.f(i) for i in M.index_set()]
            [Y(1,0)^-2 Y(1,1)^-2 Y(2,0)^2 Y(2,1) Y(3,0)^-1 Y(4,0) Y(4,1)^-1 ,
             Y(1,0)^-1 Y(1,1)^-1 Y(1,2) Y(2,0) Y(2,2)^-1 Y(3,0)^-1 Y(3,1) Y(4,0) Y(4,1)^-1 ,
             Y(1,0)^-1 Y(1,1)^-1 Y(2,0) Y(2,1)^2 Y(3,0)^-2 Y(3,1)^-1 Y(4,0)^3 Y(4,1)^-1 ,
             Y(1,0)^-1 Y(1,1)^-1 Y(2,0) Y(2,1) Y(3,0)^-1 Y(3,1) Y(4,1)^-2 ]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        newdict = copy(self._dict)
        kf = self._kf(i)
        Aik = {(i, kf):-1, (i, kf+1):-1}
        ct = self.parent().cartan_type()
        cm = ct.cartan_matrix()
        shift = 0
        if ct.is_finite():
            shift = 1
        for j in self.parent().index_set():
            if i == j:
                continue
            c = 0
            if i > j:
                c = 1
            if ct.is_affine() and ct.type() == 'A' and abs(i-j) == ct.rank() - 1:
                c = 1 - c
            if cm[j-shift][i-shift] != 0:
                Aik[(j, kf+c)] = -cm[j-shift][i-shift]
        for key,value in Aik.iteritems():
            if key in newdict:
                newdict[key] +=value
            else:
                newdict[key] = value
        for k in list(newdict):
            if newdict[k] == 0:
                newdict.pop(k)
        return self.__class__(self.parent(),newdict)

class NakajimaAMonomial(NakajimaYMonomial):
    r"""
    Monomials of the form `A_{i_1,k_1}^{a_1}\cdots A_{i_t,k_t}^{a_t}`, where
    `i_1,\dots,i_t` are elements of the index set, `k_1,\dots,k_t` are
    nonnegative integers, and `a_1,\dots,a_t` are integers.

    EXAMPLES::

        sage: M = crystals.infinity.NakajimaMonomials("A3",use_Y=False)
        sage: mg = M.module_generators[0]
        sage: mg.f_string([1,2,3,2,1])
        A(1,0)^-1 A(1,1)^-1 A(2,0)^-2 A(3,0)^-1
        sage: mg.f_string([3,2,1])
        A(1,2)^-1 A(2,1)^-1 A(3,0)^-1
    """

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['B',4,1],use_Y=False)
            sage: m = M.module_generators[0].f_string([4,2,1])
            sage: m
            A(1,1)^-1 A(2,0)^-1 A(4,0)^-1
        """
        if self._dict == {}:
            return "1"
        else:
            L = sorted(self._dict.iteritems(), key=lambda x:(x[0][0],x[0][1]))
            return_str = ''
            for x in range(len(L)):
                if L[x][1] != 1:
                    return_str += "A(%s,%s)"%(L[x][0][0],L[x][0][1]) + "^%s "%L[x][1]
                else:
                    return_str += "A(%s,%s) "%(L[x][0][0],L[x][0][1])
            return return_str

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['C',4,1],use_Y=False)
            sage: m = M.module_generators[0].f_string([4,2,3])
            sage: m._latex_()
            'A_{2,0}^{-1} A_{3,1}^{-1} A_{4,0}^{-1} '
        """
        if self._dict == {}:
            return "\\boldsymbol{1}"
        else:
            L = sorted(self._dict.iteritems(), key=lambda x:(x[0][0],x[0][1]))
            return_str = ''
            for x in range(len(L)):
                if L[x][1] !=1:
                    return_str += "A_{%s,%s}"%(L[x][0][0],L[x][0][1]) + "^{%s} "%L[x][1]
                else:
                    return_str += "A_{%s,%s} "%(L[x][0][0],L[x][0][1])
            return return_str

    def to_Y_monomial(self):
        r"""
        Represent `\prod_{(i,k)} A_{i,k}^{a_{i}(k)}` in the form
        `\prod_{(i,k)} Y_{i,k}^{y_i(k)}` using the formula

        .. MATH::

            A_{i,k} = Y_{i,k} Y_{i,k+1} \prod_{\substack{j \in I \\ j\neq i}}
            Y_{i,k+c_{ji}}^{a_{ji}}.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A',2,1],use_Y=False)
            sage: m = M.module_generators[0].f_string([2,0,1,2,1])
            sage: m
            A(0,0)^-1 A(1,0)^-1 A(1,1)^-1 A(2,0)^-1 A(2,1)^-1
            sage: m.to_Y_monomial()
            Y(0,1) Y(0,2) Y(1,1)^-1 Y(2,2)^-1
        """
        Y = {}
        d = self._dict
        ct = self.parent().cartan_type()
        cm = ct.cartan_matrix()
        for k,v in d.iteritems():
            Y[k] = Y.get(k,0) + v
            Y[(k[0],k[1]+1)] = Y.get((k[0],k[1]+1), 0) + v
            shift = 0
            if ct.is_finite():
                shift = 1
            for j in self.parent().index_set():
                if k[0] == j:
                    continue
                c = 0
                if k[0] > j:
                    c = 1
                if ct.is_affine() and ct.type() == 'A' and abs(k[0]-j) == ct.rank() - 1:
                    c = 1 - c
                if cm[j-shift][k[0]-shift] != 0:
                    Y[(j, k[1]+c)] = Y.get((j,k[1]+c),0) + v*cm[j-shift][k[0]-shift]
        for k in Y.keys():
            if Y[k] == 0:
                Y.pop(k)
        return NakajimaYMonomial(self.parent(), Y)

    def weight(self):
        r"""
        Return the weight of ``self`` as an element of
        ``self.parent().weight_lattice_realization()``.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A',4,2],use_Y=False)
            sage: m = M.module_generators[0].f_string([1,2,0,1])
            sage: m.weight()
            2*Lambda[0] - Lambda[1]
        """
        return self.to_Y_monomial().weight()

    def weight_in_root_lattice(self):
        r"""
        Return the weight of ``self`` as an element of the root lattice.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['C',3,1],use_Y=False)
            sage: m = M.module_generators[0].f_string([3,0,1,2,0])
            sage: m.weight_in_root_lattice()
            -2*alpha[0] - alpha[1] - alpha[2] - alpha[3]
        """
        return self.to_Y_monomial().weight_in_root_lattice()

    def epsilon(self,i):
        r"""
        Return the action of `\varepsilon_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['C',4,1],use_Y=False)
            sage: m = M.module_generators[0].f_string([4,2,3])
            sage: [m.epsilon(i) for i in M.index_set()]
            [0, 0, 0, 1, 0]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        return self.to_Y_monomial().epsilon(i)

    def phi(self,i):
        r"""
        Return the action of `\varphi_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['C',4,1],use_Y=False)
            sage: m = M.module_generators[0].f_string([4,2,3])
            sage: [m.phi(i) for i in M.index_set()]
            [0, 1, -1, 2, -1]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        return self.to_Y_monomial().phi(i)

    def e(self,i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,1],use_Y=False)
            sage: m = M.module_generators[0].f_string([4,2,3,0])
            sage: [m.e(i) for i in M.index_set()]
            [A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 ,
             None,
             None,
             A(0,2)^-1 A(2,1)^-1 A(4,0)^-1 ,
             None]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        if self.epsilon(i) == 0:
            return None
        ke = self.to_Y_monomial()._ke(i)
        d = copy(self._dict)
        d[(i,ke)] = d.get((i,ke),0)+1
        for k in list(d):
            if d[k] == 0:
                d.pop(k)
        return self.__class__(self.parent(), d)

    def f(self,i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials("E8",use_Y=False)
            sage: m = M.module_generators[0].f_string([4,2,3,8])
            sage: m
            A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1
            sage: [m.f(i) for i in M.index_set()]
            [A(1,2)^-1 A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1 ,
             A(2,0)^-1 A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1 ,
             A(2,1)^-1 A(3,0)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-1 ,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(4,1)^-1 A(8,0)^-1 ,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(5,0)^-1 A(8,0)^-1 ,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(6,0)^-1 A(8,0)^-1 ,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(7,1)^-1 A(8,0)^-1 ,
             A(2,1)^-1 A(3,1)^-1 A(4,0)^-1 A(8,0)^-2 ]
        """
        if i not in self.parent().index_set():
            raise ValueError("i must be an element of the index set")
        kf = self.to_Y_monomial()._kf(i)
        d = copy(self._dict)
        d[(i,kf)] = d.get((i,kf),0) - 1
        return self.__class__(self.parent(), d)

class InfinityCrystalOfNakajimaMonomials(Parent,UniqueRepresentation):
    r"""
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
    `C = (c_{ij})_{i\neq j}` such that `c_{ij}+c{ji} =1`, one defines

    .. MATH::

        A_{i,k} = Y_{i,k} Y_{i,k+1} \prod_{j\neq i} Y_{j,k+c_{ji}}^{a_{ji}},

    where `(a_{ij})` is a Cartan matrix.  Then

    .. MATH::

        \begin{aligned}
        e_iM &= \begin{cases} 0 & \text{if } \varepsilon_i(M) = 0, \\
        A_{i,k_e}M & \text{if } \varepsilon_i(M) > 0, \end{cases} \\
        f_iM &= A_{i,k_f}^{-1} M.
        \end{aligned}

    It is shown in [KKS07]_ that the connected component of
    `\widehat{\mathcal{M}}` containing the element `\boldsymbol{1}`, which we
    denote by `\mathcal{M}(\infty)`, is crystal isomorphic to the crystal
    `B(\infty)`.

    INPUT:

    - ``cartan_type`` -- a Cartan type

    - ``use_Y`` -- choice of monomials in terms of `A` or `Y`

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
    def __classcall_private__(cls, ct, category=None, use_Y=True):
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
        if isinstance(use_Y, bool):
            if use_Y:
                elt_class = NakajimaYMonomial
            else:
                elt_class = NakajimaAMonomial
        else:
            elt_class = use_Y
        cartan_type = CartanType(ct)
        return super(InfinityCrystalOfNakajimaMonomials,cls).__classcall__(cls,cartan_type,category,elt_class)

    def __init__(self, ct, category, elt_class):
        r"""
        EXAMPLES::

            sage: Minf = crystals.infinity.NakajimaMonomials(['A',3])
            sage: TestSuite(Minf).run() # long time
        """
        self._cartan_type = ct

        self.Element = elt_class
        if category is None:
            category = (HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=category)
        self.module_generators = (self.element_class(self,{}),)

    def _element_constructor_(self,dict):
        r"""
        Construct an element of ``self`` from ``dict``.

        INPUT:

        - ``dict`` -- a dictionary whose key is a pair and whose value is an integer

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,1])
            sage: m = M({(1,0):-1,(1,1):-1,(2,0):1})
            sage: m
            Y(1,0)^-1 Y(1,1)^-1 Y(2,0)
        """
        return self.element_class(self,dict)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['D',4,1])
            sage: m = M({(1,0):-1,(1,1):-1,(2,0):1})
            sage: m
            Y(1,0)^-1 Y(1,1)^-1 Y(2,0)
        """
        return "Infinity Crystal of modified Nakajima monomials of type {}".format(self._cartan_type)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``, which is always `\infty`.

        EXAMPLES::

            sage: M = crystals.infinity.NakajimaMonomials(['A',5,2])
            sage: M.cardinality()
            +Infinity
        """
        return Infinity

class CrystalOfNakajimaMonomialsElement(NakajimaYMonomial):
    r"""
    Element class for
    :class:`~sage.combinat.crystals.monomial_crystals.CrystalOfNakajimaMonomials`.

    The `f_i` operators need to be modified from the version in
    :class:`~sage.combinat.crystals.monomial_crystalsNakajimaYMonomial`
    in order to create irreducible highest weight realizations.
    This modified `f_i` is defined as

    .. MATH::

        f_iM = \begin{cases} 0 & \text{if } \varphi_i(M) = 0, \\
        A_{i,k_f}^{-1}M & \text{if } \varphi_i(M) > 0. \end{cases}

    EXAMPLES::

        sage: La = RootSystem(['A',5,2]).weight_lattice().fundamental_weights()
        sage: M = crystals.NakajimaMonomials(['A',5,2],3*La[0])
        sage: m = M.module_generators[0].f(0); m
        Y(0,0)^2 Y(0,1)^-1 Y(2,0)
        sage: TestSuite(m).run()
    """
    def f(self,i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: La = RootSystem(['A',5,2]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['A',5,2],3*La[0])
            sage: m = M.module_generators[0]
            sage: [m.f(i) for i in M.index_set()]
            [Y(0,0)^2 Y(0,1)^-1 Y(2,0) , None, None, None]
        """
        if self.phi(i) == 0:
            return None
        return super(CrystalOfNakajimaMonomialsElement, self).f(i)

class CrystalOfNakajimaMonomials(InfinityCrystalOfNakajimaMonomials):
    r"""
    Let `\widetilde{\mathcal{M}}` be `\widehat{\mathcal{M}}` as a set, and with
    crystal structure defined as on `\widehat{\mathcal{M}}` with the exception
    that

    .. MATH::

        f_iM = \begin{cases} 0 & \text{if } \varphi_i(M) = 0, \\
        A_{i,k_f}^{-1}M & \text{if } \varphi_i(M) > 0. \end{cases}

    Then Kashiwara [Kash03]_ showed that the connected component in
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

        sage: La = RootSystem(['A',3,1]).weight_lattice().fundamental_weights()
        sage: M = crystals.NakajimaMonomials(['A',3,1],La[0]+La[2])
        sage: B = crystals.GeneralizedYoungWalls(3,La[0]+La[2])
        sage: SM = M.subcrystal(max_depth=4)
        sage: SB = B.subcrystal(max_depth=4)
        sage: GM = M.digraph(subset=SM) # long time
        sage: GB = B.digraph(subset=SB) # long time
        sage: GM.is_isomorphic(GB,edge_labels=True) # long time
        True

        sage: La = RootSystem(['A',5,2]).weight_lattice().fundamental_weights()
        sage: LA = RootSystem(['A',5,2]).weight_space().fundamental_weights()
        sage: M = crystals.NakajimaMonomials(['A',5,2],3*La[0])
        sage: B = crystals.LSPaths(3*LA[0])
        sage: SM = M.subcrystal(max_depth=4)
        sage: SB = B.subcrystal(max_depth=4)
        sage: GM = M.digraph(subset=SM)
        sage: GB = B.digraph(subset=SB)
        sage: GM.is_isomorphic(GB,edge_labels=True)
        True
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, La):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: La = RootSystem(['E',8,1]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['E',8,1],La[0]+La[8])
            sage: M1 = crystals.NakajimaMonomials(CartanType(['E',8,1]),La[0]+La[8])
            sage: M2 = crystals.NakajimaMonomials(['E',8,1],M.Lambda()[0] + M.Lambda()[8])
            sage: M is M1 is M2
            True
        """
        cartan_type = CartanType(cartan_type)
        La = RootSystem(cartan_type).weight_lattice()(La)
        return super(CrystalOfNakajimaMonomials, cls).__classcall__(cls, cartan_type, La)

    def __init__(self, ct, La):
        r"""
        EXAMPLES::

            sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['A',2],La[1]+La[2])
            sage: TestSuite(M).run()
        """
        InfinityCrystalOfNakajimaMonomials.__init__( self, ct,
                (RegularCrystals(), HighestWeightCrystals()), CrystalOfNakajimaMonomialsElement )
        self._cartan_type = ct
        self.hw = La
        gen = {}
        for j in range(len(La.support())):
            gen[(La.support()[j],0)] = La.coefficients()[j]
        self.module_generators = (self.element_class(self,gen),)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: La = RootSystem(['C',3,1]).weight_lattice().fundamental_weights()
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
            sage: M = crystals.NakajimaMonomials(['A',2],La[1])
            sage: M.cardinality()
            3

            sage: La = RootSystem(['D',4,2]).weight_lattice().fundamental_weights()
            sage: M = crystals.NakajimaMonomials(['D',4,2],La[1])
            sage: M.cardinality()
            +Infinity
        """
        if not self.cartan_type().is_finite():
            return Infinity
        return len(list(self.subcrystal(generators=[self.module_generators[0]])))

