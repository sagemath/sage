# sage.doctest: needs sage.combinat sage.graphs sage.modules
r"""
Quantum Group Representations

AUTHORS:

- Travis Scrimshaw (2018): initial version
"""

#*****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.algebras.quantum_groups.q_numbers import q_int
from sage.categories.crystals import Crystals
from sage.categories.quantum_group_representations import QuantumGroupRepresentations


class QuantumGroupRepresentation(CombinatorialFreeModule):
    """
    A representation of a quantum group whose basis is indexed
    by the corresponding (combinatorial) crystal.

    INPUT:

    - ``C`` -- the crystal corresponding to the representation
    - ``R`` -- the base ring
    - ``q`` -- (default: the generator of ``R``) the parameter `q`
      of the quantum group
    """
    @staticmethod
    def __classcall__(cls, R, C, q=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
            sage: C = crystals.Tableaux(['A',3], shape=[1,1])
            sage: R = ZZ['q'].fraction_field()
            sage: V1 = MinusculeRepresentation(R, C)
            sage: V2 = MinusculeRepresentation(R, C, R.gen())
            sage: V1 is V2
            True
        """
        if q is None:
            q = R.gen()
        return super().__classcall__(cls, R, C, q)

    def __init__(self, R, C, q):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
            sage: C = crystals.Tableaux(['A',3], shape=[1,1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = MinusculeRepresentation(R, C)
            sage: TestSuite(V).run()
        """
        self._q = q
        self._d = C.cartan_type().symmetrizer()
        cat = QuantumGroupRepresentations(R).WithBasis()
        if C in Crystals().Finite():
            cat = cat.FiniteDimensional()
        CombinatorialFreeModule.__init__(self, R, C, category=cat)

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: C = crystals.Tableaux(['C',3], shape=[1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = AdjointRepresentation(R, C)
            sage: V.cartan_type()
            ['C', 3]
        """
        return self.basis().keys().cartan_type()

    def K_on_basis(self, i, b, power=1):
        r"""
        Return the action of `K_i` on the basis element indexed by ``b``
        to the power ``power``.

        INPUT:

        - ``i`` -- an element of the index set
        - ``b`` -- an element of basis keys
        - ``power`` -- (default: 1) the power of `K_i`

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
            sage: C = crystals.Tableaux(['A',3], shape=[1,1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = MinusculeRepresentation(R, C)
            sage: [[V.K_on_basis(i, b) for i in V.index_set()] for b in C]
            [[B[[[1], [2]]], q*B[[[1], [2]]], B[[[1], [2]]]],
             [q*B[[[1], [3]]], 1/q*B[[[1], [3]]], q*B[[[1], [3]]]],
             [1/q*B[[[2], [3]]], B[[[2], [3]]], q*B[[[2], [3]]]],
             [q*B[[[1], [4]]], B[[[1], [4]]], 1/q*B[[[1], [4]]]],
             [1/q*B[[[2], [4]]], q*B[[[2], [4]]], 1/q*B[[[2], [4]]]],
             [B[[[3], [4]]], 1/q*B[[[3], [4]]], B[[[3], [4]]]]]
            sage: [[V.K_on_basis(i, b, -1) for i in V.index_set()] for b in C]
            [[B[[[1], [2]]], 1/q*B[[[1], [2]]], B[[[1], [2]]]],
             [1/q*B[[[1], [3]]], q*B[[[1], [3]]], 1/q*B[[[1], [3]]]],
             [q*B[[[2], [3]]], B[[[2], [3]]], 1/q*B[[[2], [3]]]],
             [1/q*B[[[1], [4]]], B[[[1], [4]]], q*B[[[1], [4]]]],
             [q*B[[[2], [4]]], 1/q*B[[[2], [4]]], q*B[[[2], [4]]]],
             [B[[[3], [4]]], q*B[[[3], [4]]], B[[[3], [4]]]]]
        """
        WLR = self.basis().keys().weight_lattice_realization()
        alc = WLR.simple_coroots()
        return self.term( b, self._q**(b.weight().scalar(alc[i]) * self._d[i] * power) )


class CyclicRepresentation(QuantumGroupRepresentation):
    """
    A cyclic quantum group representation that is indexed by either a
    highest weight crystal or Kirillov-Reshetikhin crystal.

    The crystal ``C`` must either allow ``C.module_generator()``,
    otherwise it is assumed to be generated by ``C.module_generators[0]``.

    This is meant as an abstract base class for
    :class:`~sage.algebras.quantum_groups.representation.AdjointRepresentation`
    and
    :class:`~sage.algebras.quantum_groups.representation.MinusculeRepresentation`.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: C = crystals.Tableaux(['C',3], shape=[1])
            sage: R = ZZ['q'].fraction_field()
            sage: AdjointRepresentation(R, C)
            V((1, 0, 0))
        """
        try:
            mg = self.basis().keys().module_generator()
        except (TypeError, AttributeError):
            mg = self.basis().keys().module_generators[0]
        return "V({})".format(mg.weight())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: C = crystals.Tableaux(['G',2], shape=[1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = AdjointRepresentation(R, C)
            sage: latex(V)
            V\left( e_{0} - e_{2} \right)

            sage: La = RootSystem(['E',7,1]).weight_space().fundamental_weights()
            sage: K = crystals.ProjectedLevelZeroLSPaths(La[1])
            sage: A = AdjointRepresentation(R, K)
            sage: latex(A)
            V\left( -2 \Lambda_{0} + \Lambda_{1} \right)
        """
        try:
            mg = self.basis().keys().module_generator()
        except (TypeError, AttributeError):
            mg = self.basis().keys().module_generators[0]
        from sage.misc.latex import latex
        return r"V\left( {} \right)".format(latex(mg.weight()))

    @cached_method
    def module_generator(self):
        """
        Return the module generator of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: C = crystals.Tableaux(['G',2], shape=[1,1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = AdjointRepresentation(R, C)
            sage: V.module_generator()
            B[[[1], [2]]]

            sage: K = crystals.KirillovReshetikhin(['D',4,2], 1,1)
            sage: A = AdjointRepresentation(R, K)
            sage: A.module_generator()
            B[[[1]]]
        """
        try:
            mg = self.basis().keys().module_generator()
        except (TypeError, AttributeError):
            mg = self.basis().keys().module_generators[0]
        return self.monomial(mg)


class AdjointRepresentation(CyclicRepresentation):
    r"""
    An (generalized) adjoint representation of a quantum group.

    We define an *(generalized) adjoint representation* `V` of a
    quantum group `U_q` to be a cyclic `U_q`-module with a weight
    space decomposition `V = \bigoplus_{\mu} V_{\mu}` such that
    `\dim V_{\mu} \leq 1` unless `\mu = 0`. Moreover, we require
    that there exists a basis `\{y_j | j \in J\}` for `V_0` such
    that `e_i y_j = 0` for all `j \neq i \in I`.

    For a base ring `R`, we construct an adjoint representation from
    its (combinatorial) crystal `B` by `V = R \{v_b | b \in B\}` with

    .. MATH::

        \begin{aligned}
        e_i v_b & = \begin{cases}
            v_{e_i b} / [\varphi_i(e_i b)]_{q_i},
                & \text{if } \operatorname{wt}(b) \neq 0, \\
            v_{e_i b} + \sum_{j \neq i} [-A_{ij}]_{q_i} / [2]_{q_i} v_{y_j}
                & \text{otherwise}
            \end{cases} \\
        f_i v_b & = \begin{cases}
            v_{f_i b} / [\varepsilon_i(f_i b)]_{q_i},
                & \text{if } \operatorname{wt}(b) \neq 0, \\
            v_{f_i b} + \sum_{j \neq i} [-A_{ij}]_{q_i} / [2]_{q_i} v_{y_j}
                & \text{otherwise}
            \end{cases} \\
        K_i v_b & = q^{\langle h_i, \operatorname{wt}(b) \rangle} v_b,
        \end{aligned}

    where `(A_{ij})_{i,j \in I}` is the Cartan matrix, and we
    consider `v_0 := 0`.

    INPUT:

    - ``C`` -- the crystal corresponding to the representation
    - ``R`` -- the base ring
    - ``q`` -- (default: the generator of ``R``) the parameter `q`
      of the quantum group

    .. WARNING::

        This assumes that `q` is generic.

    EXAMPLES::

        sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
        sage: R = ZZ['q'].fraction_field()
        sage: C = crystals.Tableaux(['D',4], shape=[1,1])
        sage: V = AdjointRepresentation(R, C)
        sage: V
        V((1, 1, 0, 0))
        sage: v = V.an_element(); v
        2*B[[[1], [2]]] + 2*B[[[1], [3]]] + 3*B[[[2], [3]]]
        sage: v.e(2)
        2*B[[[1], [2]]]
        sage: v.f(2)
        2*B[[[1], [3]]]
        sage: v.f(4)
        2*B[[[1], [-4]]] + 3*B[[[2], [-4]]]
        sage: v.K(3)
        2*B[[[1], [2]]] + 2*q*B[[[1], [3]]] + 3*q*B[[[2], [3]]]
        sage: v.K(2,-2)
        2/q^2*B[[[1], [2]]] + 2*q^2*B[[[1], [3]]] + 3*B[[[2], [3]]]

        sage: La = RootSystem(['F',4,1]).weight_space().fundamental_weights()
        sage: K = crystals.ProjectedLevelZeroLSPaths(La[4])
        sage: A = AdjointRepresentation(R, K)
        sage: A
        V(-Lambda[0] + Lambda[4])

    Sort the summands uniformly in Python 2 and Python 3::

        sage: A.print_options(sorting_key=lambda x: str(x))
        sage: v = A.an_element(); v
        2*B[(-Lambda[0] + Lambda[3] - Lambda[4],)]
         + 2*B[(-Lambda[0] + Lambda[4],)]
         + 3*B[(Lambda[0] - Lambda[1] + Lambda[4],)]
        sage: v.e(0)
        2*B[(Lambda[0] - Lambda[1] + Lambda[3] - Lambda[4],)]
         + 2*B[(Lambda[0] - Lambda[1] + Lambda[4],)]
        sage: v.f(0)
        3*B[(-Lambda[0] + Lambda[4],)]

    REFERENCES:

    - [OS2018]_
    """
    def __init__(self, R, C, q):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: C = crystals.Tableaux(['B',3], shape=[1,1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = AdjointRepresentation(R, C)
            sage: TestSuite(V).run()

            sage: A = crystals.Tableaux(['A',2], shape=[2,1])
            sage: VA = AdjointRepresentation(R, A)
            sage: TestSuite(VA).run()

            sage: K1 = crystals.KirillovReshetikhin(['C',3,1], 1,1)
            sage: A1 = AdjointRepresentation(R, K1)
            sage: TestSuite(A1).run()
            sage: K2 = crystals.KirillovReshetikhin(['C',2,1], 1,2)
            sage: A2 = AdjointRepresentation(R, K2)
            sage: TestSuite(A2).run()
        """
        self._WLR_zero = C.weight_lattice_realization().zero()
        CyclicRepresentation.__init__(self, R, C, q)
        ct = C.cartan_type()
        if ct.is_finite() and ct.type() == 'A':
            def test_zero(x):
                wt = x.weight()
                return all(wt.scalar(ac) == 0
                           for ac in self._WLR_zero.parent().simple_coroots())
            self._check_zero_wt = test_zero
        else:
            self._check_zero_wt = lambda x: x.weight() == self._WLR_zero

    @lazy_attribute
    def _zero_elts(self):
        r"""
        Find all of the elements of weight `0` in the basis keys.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,1)
            sage: R = ZZ['q'].fraction_field()
            sage: V = AdjointRepresentation(R, K)
            sage: V._zero_elts
            {0: [], 1: [[2], [-2]], 2: [[3], [-3]],
             3: [[4], [-4]], 4: [[-4], [4]]}
        """
        C = self.basis().keys()
        ret = {}
        for x in C:
            if self._check_zero_wt(x):
                for i in C.index_set():
                    if x.epsilon(i) > 0:
                        ret[i] = x
                        break
        return ret

    def e_on_basis(self, i, b):
        r"""
        Return the action of `e_i` on the basis element indexed by ``b``.

        INPUT:

        - ``i`` -- an element of the index set
        - ``b`` -- an element of basis keys

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: K = crystals.KirillovReshetikhin(['D',3,2], 1,1)
            sage: R = ZZ['q'].fraction_field()
            sage: V = AdjointRepresentation(R, K)
            sage: mg0 = K.module_generators[0]; mg0
            []
            sage: mg1 = K.module_generators[1]; mg1
            [[1]]
            sage: V.e_on_basis(0, mg0)
            ((q^2+1)/q)*B[[[-1]]]
            sage: V.e_on_basis(0, mg1)
            B[[]]
            sage: V.e_on_basis(1, mg0)
            0
            sage: V.e_on_basis(1, mg1)
            0
            sage: V.e_on_basis(2, mg0)
            0
            sage: V.e_on_basis(2, mg1)
            0

            sage: K = crystals.KirillovReshetikhin(['D',4,3], 1,1)
            sage: V = AdjointRepresentation(R, K)
            sage: V.e_on_basis(0, K.module_generator())
            B[[]] + (q/(q^2+1))*B[[[0]]]
        """
        C = self.basis().keys()
        x = b.e(i)
        if x is None:
            return self.zero()
        I = {j: pos for pos,j in enumerate(C.index_set())}
        if self._check_zero_wt(x):
            A = C.cartan_type().cartan_matrix()
            return self.monomial(x) + sum(self.term(self._zero_elts[j],
                                                    q_int(-A[I[i],I[j]], self._q**self._d[i])
                                                    / q_int(2, self._q**self._d[j]))
                                          for j in C.index_set()
                                          if A[I[i],I[j]] < 0 and j in self._zero_elts)
        return self.term(x, q_int(x.phi(i), self._q**self._d[i]))

    def f_on_basis(self, i, b):
        r"""
        Return the action of `f_i` on the basis element indexed by ``b``.

        INPUT:

        - ``i`` -- an element of the index set
        - ``b`` -- an element of basis keys

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
            sage: K = crystals.KirillovReshetikhin(['D',3,2], 1,1)
            sage: R = ZZ['q'].fraction_field()
            sage: V = AdjointRepresentation(R, K)
            sage: mg0 = K.module_generators[0]; mg0
            []
            sage: mg1 = K.module_generators[1]; mg1
            [[1]]
            sage: V.f_on_basis(0, mg0)
            ((q^2+1)/q)*B[[[1]]]
            sage: V.f_on_basis(0, mg1)
            0
            sage: V.f_on_basis(1, mg0)
            0
            sage: V.f_on_basis(1, mg1)
            B[[[2]]]
            sage: V.f_on_basis(2, mg0)
            0
            sage: V.f_on_basis(2, mg1)
            0

            sage: K = crystals.KirillovReshetikhin(['D',4,3], 1,1)
            sage: V = AdjointRepresentation(R, K)
            sage: lw = K.module_generator().to_lowest_weight([1,2])[0]
            sage: V.f_on_basis(0, lw)
            B[[]] + (q/(q^2+1))*B[[[0]]]
        """
        C = self.basis().keys()
        x = b.f(i)
        if x is None:
            return self.zero()
        I = {j: pos for pos,j in enumerate(C.index_set())}
        if self._check_zero_wt(x):
            A = C.cartan_type().cartan_matrix()
            return self.monomial(x) + sum(self.term(self._zero_elts[j],
                                                    q_int(-A[I[i],I[j]],
                                                    self._q**self._d[i])
                                                    / q_int(2, self._q**self._d[j]))
                                          for j in C.index_set()
                                          if A[I[i],I[j]] < 0 and j in self._zero_elts)
        return self.term(x, q_int(x.epsilon(i), self._q**self._d[i]))


class MinusculeRepresentation(CyclicRepresentation):
    r"""
    A minuscule representation of a quantum group.

    A quantum group representation `V` is *minuscule* if it is
    cyclic, there is a weight space decomposition
    `V = \bigoplus_{\mu} V_{\mu}` with `\dim V_{\mu} \leq 1`,
    and `e_i^2 V = 0` and `f_i^2 V = 0`.

    For a base ring `R`, we construct a minuscule representation from
    its (combinatorial) crystal `B` by `V = R \{v_b | b \in B\}` with
    `e_i v_b = v_{e_i b}`, `f_i v_b = v_{f_i b}`, and
    `K_i v_b = q^{\langle h_i, \operatorname{wt}(b) \rangle} v_b`,
    where we consider `v_0 := 0`.

    INPUT:

    - ``C`` -- the crystal corresponding to the representation
    - ``R`` -- the base ring
    - ``q`` -- (default: the generator of ``R``) the parameter `q`
      of the quantum group

    .. WARNING::

        This assumes that `q` is generic.

    EXAMPLES::

        sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
        sage: R = ZZ['q'].fraction_field()
        sage: C = crystals.Tableaux(['B',3], shape=[1/2,1/2,1/2])
        sage: V = MinusculeRepresentation(R, C)
        sage: V
        V((1/2, 1/2, 1/2))
        sage: v = V.an_element(); v
        2*B[[+++, []]] + 2*B[[++-, []]] + 3*B[[+-+, []]]
        sage: v.e(3)
        2*B[[+++, []]]
        sage: v.f(1)
        3*B[[-++, []]]
        sage: v.f(3)
        2*B[[++-, []]] + 3*B[[+--, []]]
        sage: v.K(2)
        2*B[[+++, []]] + 2*q^2*B[[++-, []]] + 3/q^2*B[[+-+, []]]
        sage: v.K(3, -2)
        2/q^2*B[[+++, []]] + 2*q^2*B[[++-, []]] + 3/q^2*B[[+-+, []]]

        sage: K = crystals.KirillovReshetikhin(['D',4,2], 3,1)
        sage: A = MinusculeRepresentation(R, K)
        sage: A
        V(-Lambda[0] + Lambda[3])
        sage: v = A.an_element(); v
        2*B[[+++, []]] + 2*B[[++-, []]] + 3*B[[+-+, []]]
        sage: v.f(0)
        0
        sage: v.e(0)
        2*B[[-++, []]] + 2*B[[-+-, []]] + 3*B[[--+, []]]

    REFERENCES:

    - [OS2018]_
    """
    def e_on_basis(self, i, b):
        r"""
        Return the action of `e_i` on the basis element indexed by ``b``.

        INPUT:

        - ``i`` -- an element of the index set
        - ``b`` -- an element of basis keys

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
            sage: C = crystals.Tableaux(['A',3], shape=[1,1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = MinusculeRepresentation(R, C)
            sage: lw = C.lowest_weight_vectors()[0]
            sage: V.e_on_basis(1, lw)
            0
            sage: V.e_on_basis(2, lw)
            B[[[2], [4]]]
            sage: V.e_on_basis(3, lw)
            0
            sage: hw = C.highest_weight_vectors()[0]
            sage: all(V.e_on_basis(i, hw) == V.zero() for i in V.index_set())
            True
        """
        x = b.e(i)
        if x is None:
            return self.zero()
        return self.monomial(x)

    def f_on_basis(self, i, b):
        r"""
        Return the action of `f_i` on the basis element indexed by ``b``.

        INPUT:

        - ``i`` -- an element of the index set
        - ``b`` -- an element of basis keys

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
            sage: C = crystals.Tableaux(['A',3], shape=[1,1])
            sage: R = ZZ['q'].fraction_field()
            sage: V = MinusculeRepresentation(R, C)
            sage: hw = C.highest_weight_vectors()[0]
            sage: V.f_on_basis(1, hw)
            0
            sage: V.f_on_basis(2, hw)
            B[[[1], [3]]]
            sage: V.f_on_basis(3, hw)
            0
            sage: lw = C.lowest_weight_vectors()[0]
            sage: all(V.f_on_basis(i, lw) == V.zero() for i in V.index_set())
            True
        """
        x = b.f(i)
        if x is None:
            return self.zero()
        return self.monomial(x)
