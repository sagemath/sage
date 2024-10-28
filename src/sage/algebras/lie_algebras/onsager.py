"""
Onsager Algebra

AUTHORS:

- Travis Scrimshaw (2017-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.structure.indexed_generators import IndexedGenerators
from sage.sets.family import Family
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import LieAlgebraWithGenerators, InfinitelyGeneratedLieAlgebra


class OnsagerAlgebra(LieAlgebraWithGenerators, IndexedGenerators):
    r"""
    The Onsager (Lie) algebra.

    The Onsager (Lie) algebra `\mathcal{O}` is a Lie algebra with
    generators `A_0, A_1` that satisfy

    .. MATH::

        [A_0, [A_0, [A_0, A_1]]] = -4 [A_0, A_1],
        \qquad
        [A_1, [A_1, [A_1, A_0]]] = -4 [A_1, A_0].

    .. NOTE::

        We are using a rescaled version of the usual defining generators.

    There exist a basis `\{A_m, G_n \mid m \in \ZZ, n \in \ZZ_{>0}\}`
    for `\mathcal{O}` with structure coefficients

    .. MATH::

        [A_m, A_{m'}] = G_{m-m'},
        \qquad
        [G_n, G_{n'}] = 0,
        \qquad
        [G_n, A_m] = 2A_{m-n} - 2A_{m+n},

    where `m > m'`.

    The Onsager algebra is isomorphic to the subalgebra of the affine
    Lie algebra `\widehat{\mathfrak{sl}}_2 = \mathfrak{sl}_2 \otimes
    \CC[t,t^{-1}] \oplus \CC K \oplus \CC d` that is invariant under
    the Chevalley involution. In particular, we have

    .. MATH::

        A_i \mapsto f \otimes t^i - e \otimes t^{-i},
        \qquad
        G_i \mapsto h \otimes t^{-i} - h \otimes t^i.

    where `e,f,h` are the Chevalley generators of `\mathfrak{sl}_2`.

    EXAMPLES:

    We construct the Onsager algebra and do some basic computations::

        sage: O = lie_algebras.OnsagerAlgebra(QQ)
        sage: O.inject_variables()
        Defining A0, A1

    We verify the defining relations::

        sage: O([A0, [A0, [A0, A1]]]) == -4 * O([A0, A1])
        True
        sage: O([A1, [A1, [A1, A0]]]) == -4 * O([A1, A0])
        True

    We check the embedding into `\widehat{\mathfrak{sl}}_2`::

        sage: L = LieAlgebra(QQ, cartan_type=['A',1,1])
        sage: B = L.basis()
        sage: al = RootSystem(['A',1]).root_lattice().simple_root(1)
        sage: ac = al.associated_coroot()
        sage: def emb_A(i): return B[-al,i] - B[al,-i]
        sage: def emb_G(i): return B[ac,i] - B[ac,-i]
        sage: a0 = emb_A(0)
        sage: a1 = emb_A(1)
        sage: L([a0, [a0, [a0, a1]]]) == -4 * L([a0, a1])
        True
        sage: L([a1, [a1, [a1, a0]]]) == -4 * L([a1, a0])
        True

        sage: all(emb_G(n).bracket(emb_A(m)) == 2*emb_A(m-n) - 2*emb_A(m+n)
        ....:     for m in range(-10, 10) for n in range(1,10))
        True
        sage: all(emb_A(m).bracket(emb_A(mp)) == emb_G(m-mp)
        ....:     for m in range(-10,10) for mp in range(m-10, m))
        True

    REFERENCES:

    - [Onsager1944]_
    - [DG1982]_
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: TestSuite(O).run()
        """
        cat = LieAlgebras(R).WithBasis()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        IndexedGenerators.__init__(self, FiniteEnumeratedSet([0,1]))
        LieAlgebraWithGenerators.__init__(self, R, index_set=self._indices,
                                          names=('A0', 'A1'), category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.OnsagerAlgebra(QQ)
            Onsager algebra over Rational Field
        """
        return "Onsager algebra over {}".format(self.base_ring())

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: latex(O)
            \mathcal{O}_{\Bold{Q}}
        """
        from sage.misc.latex import latex
        return "\\mathcal{{O}}_{{{}}}".format(latex(self.base_ring()))

    def _repr_generator(self, m):
        """
        Return a string representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O._repr_generator((0,-2))
            'A[-2]'
            sage: O._repr_generator((1,4))
            'G[4]'
        """
        if m[0] == 0:
            return 'A[{}]'.format(m[1])
        return 'G[{}]'.format(m[1])

    def _latex_generator(self, m):
        r"""
        Return a LaTeX representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O._latex_generator((0,-2))
            'A_{-2}'
            sage: O._latex_generator((1,4))
            'G_{4}'
        """
        if m[0] == 0:
            return 'A_{{{}}}'.format(m[1])
        return 'G_{{{}}}'.format(m[1])

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = _repr_generator
    _latex_term = _latex_generator

    @cached_method
    def basis(self):
        r"""
        Return the basis of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.basis()
            Lazy family (Onsager monomial(i))_{i in
             Disjoint union of Family (Integer Ring, Positive integers)}
        """
        from sage.rings.integer_ring import ZZ
        from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
        from sage.sets.positive_integers import PositiveIntegers
        I = DisjointUnionEnumeratedSets([ZZ, PositiveIntegers()],
                                        keepkey=True, facade=True)
        return Family(I, self.monomial, name='Onsager monomial')

    @cached_method
    def lie_algebra_generators(self):
        r"""
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.lie_algebra_generators()
            Finite family {'A0': A[0], 'A1': A[1]}
        """
        d = {"A0": self.basis()[0,0], "A1": self.basis()[0,1]}
        return Family(self._names, d.__getitem__)

    def bracket_on_basis(self, x, y):
        r"""
        Return the bracket of basis elements indexed by ``x`` and ``y``
        where ``x < y``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.bracket_on_basis((1,3), (1,9))  # [G, G]
            0
            sage: O.bracket_on_basis((0,8), (1,13))  # [A, G]
            -2*A[-5] + 2*A[21]
            sage: O.bracket_on_basis((0,-9), (0, 7))  # [A, A]
            -G[16]
        """
        if x[0] == 1:
            # From < property, we have y[0] == 1
            # Therefore, we have [G_n, G_{n'}] = 0
            return self.zero()
        R = self.base_ring()
        if y[0] == 1: # [A_m, G_n] = -(2A_{m-n} - 2A_{m+n})
            d = {(0, x[1]-y[1]): R(-2), (0, x[1]+y[1]): R(2)}
            return self.element_class(self, d)
        # [A_m, A_{m'}] = -G_{m' - m}, where m < m'
        return self.element_class(self, {(1, y[1]-x[1]): -R.one()})

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.an_element()
            -2*A[-3] + A[2] + 3*G[2]
        """
        B = self.basis()
        return B[0,2] - 2*B[0,-3] + 3*B[1,2]

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.some_elements()
            [A[0], A[2], A[-1], G[4], -2*A[-3] + A[2] + 3*G[2]]
        """
        B = self.basis()
        return [B[0,0], B[0,2], B[0,-1], B[1,4], self.an_element()]

    def quantum_group(self, q=None, c=None):
        r"""
        Return the quantum group of ``self``.

        The corresponding quantum group is the
        :class:`~sage.algebras.lie_algebras.onsager.QuantumOnsagerAlgebra`.
        The parameter `c` must be such that `c(1) = 1`

        INPUT:

        - ``q`` -- (optional) the quantum parameter; the default
          is `q \in R(q)`, where `R` is the base ring of ``self``
        - ``c`` -- (optional) the parameter `c`; the default is ``q``

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: Q
            q-Onsager algebra with c=q over Fraction Field of
             Univariate Polynomial Ring in q over Rational Field
        """
        if q is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            q = PolynomialRing(self.base_ring(), 'q').fraction_field().gen()
        if c is None:
            c = q
        else:
            c = q.parent()(c)
        return QuantumOnsagerAlgebra(self, q, c)

    def alternating_central_extension(self):
        r"""
        Return the alternating central extension of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: ACE = O.alternating_central_extension()
            sage: ACE
            Alternating central extension of the Onsager algebra over Rational Field
        """
        return OnsagerAlgebraACE(self.base_ring())

    Element = LieAlgebraElement

#####################################################################
# q-Onsager algebra (the quantum group)


class QuantumOnsagerAlgebra(CombinatorialFreeModule):
    r"""
    The quantum Onsager algebra.

    The *quantum Onsager algebra*, or `q`-Onsager algebra, is a
    quantum group analog of the Onsager algebra. It is the left
    (or right) coideal subalgebra of the quantum group
    `U_q(\widehat{\mathfrak{sl}}_2)` and is the simplest example
    of a quantum symmetric pair coideal subalgebra of affine type.

    The `q`-Onsager algebra depends on a parameter `c` such that
    `c(1) = 1`. The `q`-Onsager algebra with parameter `c` is denoted
    `U_q(\mathcal{O}_R)_c`, where `R` is the base ring of the
    defining Onsager algebra.

    EXAMPLES:

    We create the `q`-Onsager algebra and its generators::

        sage: O = lie_algebras.OnsagerAlgebra(QQ)
        sage: Q = O.quantum_group()
        sage: G = Q.algebra_generators()

    The generators are given as pairs, where `G[0,n]` is the generator
    `B_{n\delta+\alpha_1}` and `G[1,n]` is the generator `B_{n\delta}`.
    We use the convention that
    `n\delta + \alpha_1 \equiv (-n-1)\delta + \alpha_0`. ::

        sage: G[0,5]
        B[5d+a1]
        sage: G[0,-5]
        B[4d+a0]
        sage: G[1,5]
        B[5d]
        sage: (G[0,5] + G[0,-3]) * (G[1,2] - G[0,3])
        B[2d+a0]*B[2d] - B[2d+a0]*B[3d+a1]
         + ((-q^4+1)/q^2)*B[1d]*B[6d+a1]
         + ((q^4-1)/q^2)*B[1d]*B[4d+a1] + B[2d]*B[5d+a1]
         - B[5d+a1]*B[3d+a1] + ((q^2+1)/q^2)*B[7d+a1]
         + ((q^6+q^4-q^2-1)/q^2)*B[5d+a1] + (-q^4-q^2)*B[3d+a1]
        sage: (G[0,5] + G[0,-3] + G[1,4]) * (G[0,2] - G[1,3])
        -B[2d+a0]*B[3d] + B[2d+a0]*B[2d+a1]
         + ((q^4-1)/q^4)*B[1d]*B[7d+a1]
         + ((q^8-2*q^4+1)/q^4)*B[1d]*B[5d+a1]
         + (-q^4+1)*B[1d]*B[3d+a1] + ((q^4-1)/q^2)*B[2d]*B[6d+a1]
         + ((-q^4+1)/q^2)*B[2d]*B[4d+a1] - B[3d]*B[4d]
         - B[3d]*B[5d+a1] + B[4d]*B[2d+a1] + B[5d+a1]*B[2d+a1]
         + ((-q^2-1)/q^4)*B[8d+a1] + ((-q^6-q^4+q^2+1)/q^4)*B[6d+a1]
         + (-q^6-q^4+q^2+1)*B[4d+a1] + (q^6+q^4)*B[2d+a1]

    We check the `q`-Dolan-Grady relations::

        sage: def q_dolan_grady(a, b, q):
        ....:     x = q*a*b - ~q*b*a
        ....:     y = ~q*a*x - q*x*a
        ....:     return a*y - y*a
        sage: A0, A1 = G[0,-1], G[0,0]
        sage: q = Q.q()
        sage: q_dolan_grady(A1, A0, q) == (q^4 + 2*q^2 + 1) * (A0*A1 - A1*A0)
        True
        sage: q_dolan_grady(A0, A1, q) == (q^4 + 2*q^2 + 1) * (A1*A0 - A0*A1)
        True

    REFERENCES:

    - [BK2017]_
    """
    def __init__(self, g, q, c):
        """
        Initialize ``self``.

        TESTS::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: TestSuite(Q).run()  # long time
        """
        self._g = g
        self._q = q
        self._c = c
        self._q_two = q + ~q
        R = self._q_two.parent()
        from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
        monomials = IndexedFreeAbelianMonoid(g.basis().keys(),
                                             prefix='B', bracket=False,
                                             sorting_key=self._monoid_key)
        CombinatorialFreeModule.__init__(self, R, monomials,
                                         prefix='', bracket=False, latex_bracket=False,
                                         sorting_key=self._monomial_key,
                                         category=Algebras(R).WithBasis().Filtered())

    def _basis_key(self, k):
        r"""
        Key for ordering the basis elements of ``self._g``.

        We choose a key in order to obtain the ordering from [BK2017]_
        in the quantum group.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: Q._basis_key((0,2))
            (1, -2)
            sage: Q._basis_key((0,-2))
            (-1, 2)
            sage: Q._basis_key((1,2))
            (0, 2)
        """
        if k[0] == 0: # B_{m\delta + \alpha_1}
            if k[1] < 0:
                return (-1, -k[1])
            else:
                return (1, -k[1])
        # B_{n\delta}
        return (0, k[1])

    def _monoid_key(self, x):
        r"""
        Key function for the underlying monoid of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: G = Q.algebra_generators()
            sage: I = Q._indices.gens()
            sage: I[0,1] * I[1,3] * I[1,2] * I[0,-4]^3 # indirect doctest
            B(0, -4)^3*B(1, 2)*B(1, 3)*B(0, 1)
        """
        return self._basis_key(x[0])

    def _monomial_key(self, x):
        r"""
        Compute the key for ``x`` so that the comparison is done by
        reverse degree lexicographic order.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: G = Q.algebra_generators()
            sage: G[0,0] * G[1,1] * G[0,-2]  # indirect doctest
            (q^2-1)*B[a0]^2*B[1d] + q^2*B[1d+a0]*B[1d]*B[a1]
             + ((q^6-2*q^2-1)/q^2)*B[a0]*B[1d+a0] + (-q^4-q^2)*B[a0]*B[a1]
             + (q^4+q^2)*B[1d+a0]*B[1d+a1] + (q^4+q^2)*B[2d+a0]*B[a1]
             + q^2*B[1d]*B[2d] + (-q^4+1)*B[1d] + (q^4+q^2)*B[3d]
        """
        return (-len(x), [self._basis_key(l) for l in x.to_word_list()])

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.quantum_group()
            q-Onsager algebra with c=q over Fraction Field of
             Univariate Polynomial Ring in q over Rational Field
        """
        return "{}-Onsager algebra with c={} over {}".format(self._q, self._c,
                                                             self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group(q=-1)
            sage: latex(Q)
            U_{-1}(\mathcal{O}_{\Bold{Q}})_{-1}
        """
        from sage.misc.latex import latex
        return "U_{{{}}}(\\mathcal{{O}}_{{{}}})_{{{}}}".format(latex(self._q),
                            latex(self._g.base_ring()), latex(self._c))

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: I = Q._indices.gens()
            sage: Q._repr_term(I[0,3])
            'B[3d+a1]'
            sage: Q._repr_term(I[0,-3])
            'B[2d+a0]'
            sage: Q._repr_term(I[1,3])
            'B[3d]'
            sage: Q._repr_term(I[0,-1]^2 * I[1,3]^13 * I[0,3])
            'B[a0]^2*B[3d]^13*B[3d+a1]'
        """
        def to_str(x):
            k,e = x
            if k[0] == 0:
                if k[1] == -1:
                    ret = 'B[a0]'
                elif k[1] == 0:
                    ret = 'B[a1]'
                elif k[1] < -1:
                    ret = 'B[{}d+a0]'.format(-k[1]-1)
                elif k[1] > 0:
                    ret = 'B[{}d+a1]'.format(k[1])
            else:
                ret = 'B[{}d]'.format(k[1])
            if e > 1:
                ret = ret + '^{}'.format(e)
            return ret
        return '*'.join(to_str(x) for x in m._sorted_items())

    def _latex_term(self, m):
        r"""
        Return a latex representation of the term indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: I = Q._indices.gens()
            sage: Q._latex_term(I[0,3])
            'B_{3\\delta+\\alpha_1}'
            sage: Q._latex_term(I[0,-3])
            'B_{2\\delta+\\alpha_0}'
            sage: Q._latex_term(I[1,3])
            'B_{3\\delta}'
            sage: Q._latex_term(I[0,-1]^2 * I[1,3]^13 * I[0,3])
            'B_{\\alpha_0}^{2} B_{3\\delta}^{13} B_{3\\delta+\\alpha_1}'
        """
        def to_str(x):
            k,e = x
            if k[0] == 0:
                if k[1] == -1:
                    ret = 'B_{\\alpha_0}'
                elif k[1] == 0:
                    ret = 'B_{\\alpha_1}'
                elif k[1] < -1:
                    ret = 'B_{{{}\\delta+\\alpha_0}}'.format(-k[1]-1)
                elif k[1] > 0:
                    ret = 'B_{{{}\\delta+\\alpha_1}}'.format(k[1])
            else:
                ret = 'B_{{{}\\delta}}'.format(k[1])
            if e > 1:
                ret = ret + '^{{{}}}'.format(e)
            return ret
        return ' '.join(to_str(x) for x in m._sorted_items())

    def lie_algebra(self):
        r"""
        Return the underlying Lie algebra of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: Q.lie_algebra()
            Onsager algebra over Rational Field
            sage: Q.lie_algebra() is O
            True
        """
        return self._g

    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: Q.algebra_generators()
            Lazy family (generator map(i))_{i in Disjoint union of
             Family (Integer Ring, Positive integers)}
        """
        G = self._indices.gens()
        return Family(self._indices._indices, lambda x: self.monomial(G[x]),
                      name="generator map")

    gens = algebra_generators

    def q(self):
        """
        Return the parameter `q` of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: Q.q()
            q
        """
        return self._q

    def c(self):
        """
        Return the parameter `c` of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group(c=-3)
            sage: Q.c()
            -3
        """
        return self._c

    @cached_method
    def one_basis(self):
        """
        Return the basis element indexing `1`.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: ob = Q.one_basis(); ob
            1
            sage: ob.parent()
            Free abelian monoid indexed by
             Disjoint union of Family (Integer Ring, Positive integers)
        """
        return self._indices.one()

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: Q.an_element()
            -2*B[2d+a0] + q*B[2d] + B[2d+a1]
        """
        G = self.algebra_generators()
        return G[0,2] - 2*G[0,-3] + self.base_ring().an_element()*G[1,2]

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: Q.some_elements()
            [B[a1], B[3d+a1], B[a0], B[1d], B[4d]]
        """
        G = self.algebra_generators()
        return [G[0,0], G[0,3], G[0,-1], G[1,1], G[1,4]]

    def degree_on_basis(self, m):
        r"""
        Return the degree of the basis element indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: G = Q.algebra_generators()
            sage: B0 = G[0,0]
            sage: B1 = G[0,-1]
            sage: Q.degree_on_basis(B0.leading_support())
            1
            sage: Q.degree_on_basis((B1^10 * B0^10).leading_support())
            20
            sage: ((B0 * B1)^3).maximal_degree()
            6
        """
        return m.length()

    @cached_method
    def product_on_basis(self, lhs, rhs):
        r"""
        Return the product of the two basis elements ``lhs`` and ``rhs``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: I = Q._indices.gens()
            sage: Q.product_on_basis(I[1,21]^2, I[1,31]^3)
            B[21d]^2*B[31d]^3
            sage: Q.product_on_basis(I[1,31]^3, I[1,21]^2)
            B[21d]^2*B[31d]^3
            sage: Q.product_on_basis(I[0,8], I[0,6])
            B[8d+a1]*B[6d+a1]
            sage: Q.product_on_basis(I[0,-8], I[0,6])
            B[7d+a0]*B[6d+a1]
            sage: Q.product_on_basis(I[0,-6], I[0,-8])
            B[5d+a0]*B[7d+a0]
            sage: Q.product_on_basis(I[0,-6], I[1,2])
            B[5d+a0]*B[2d]
            sage: Q.product_on_basis(I[1,6], I[0,2])
            B[6d]*B[2d+a1]

            sage: Q.product_on_basis(I[0,1], I[0,2])
            1/q^2*B[2d+a1]*B[1d+a1] - B[1d]
            sage: Q.product_on_basis(I[0,-3], I[0,-1])
            1/q^2*B[a0]*B[2d+a0] + ((-q^2+1)/q^2)*B[1d+a0]^2 - B[2d]
            sage: Q.product_on_basis(I[0,2], I[0,-1])
            q^2*B[a0]*B[2d+a1] + ((q^4-1)/q^2)*B[1d+a1]*B[a1]
             + (-q^2+1)*B[1d] + q^2*B[3d]
            sage: Q.product_on_basis(I[0,2], I[1,1])
            B[1d]*B[2d+a1] + (q^2+1)*B[3d+a1] + (-q^2-1)*B[1d+a1]
            sage: Q.product_on_basis(I[0,1], I[1,2])
            ((-q^4+1)/q^2)*B[1d]*B[2d+a1] + ((q^4-1)/q^2)*B[1d]*B[a1]
             + B[2d]*B[1d+a1] + (-q^4-q^2)*B[a0]
             + ((q^2+1)/q^2)*B[3d+a1] + ((q^6+q^4-q^2-1)/q^2)*B[1d+a1]
            sage: Q.product_on_basis(I[1,2], I[0,-1])
            B[a0]*B[2d] + ((-q^4+1)/q^2)*B[1d+a0]*B[1d]
             + ((q^4-1)/q^2)*B[1d]*B[a1] + ((q^2+1)/q^2)*B[2d+a0]
             + ((-q^2-1)/q^2)*B[1d+a1]
            sage: Q.product_on_basis(I[1,2], I[0,-4])
            ((q^4-1)/q^2)*B[2d+a0]*B[1d] + B[3d+a0]*B[2d]
             + ((-q^4+1)/q^2)*B[4d+a0]*B[1d] + (-q^4-q^2)*B[1d+a0]
             + ((q^6+q^4-q^2-1)/q^2)*B[3d+a0] + ((q^2+1)/q^2)*B[5d+a0]

        TESTS::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: Q = O.quantum_group()
            sage: G = Q.gens()
            sage: G[0,2]*(G[0,1]*G[0,3]) - (G[0,2]*G[0,1])*G[0,3]
            0
            sage: G[0,-2]*(G[0,-1]*G[0,-3]) - (G[0,-2]*G[0,-1])*G[0,-3]
            0
            sage: G[0,1]*(G[0,3]*G[0,-2]) - (G[0,1]*G[0,3])*G[0,-2]
            0
            sage: G[0,2]*(G[0,1]*G[1,3]) - (G[0,2]*G[0,1])*G[1,3]
            0
            sage: G[0,-2]*(G[0,1]*G[1,3]) - (G[0,-2]*G[0,1])*G[1,3]
            0
            sage: G[0,-2]*(G[1,1]*G[1,3]) - (G[0,-2]*G[1,1])*G[1,3]
            0
        """
        # Some trivial base cases
        if lhs == self.one_basis():
            return self.monomial(rhs)
        if rhs == self.one_basis():
            return self.monomial(lhs)

        I = self._indices
        B = I.gens()
        q = self._q
        kl = lhs.trailing_support()
        kr = rhs.leading_support()
        if self._basis_key(kl) <= self._basis_key(kr):
            return self.monomial(lhs * rhs)

        # Create the commutator
        # We have xy - yx = [x, y] -> xy = yx + LOT for x > y
        if kl[0] == 1 and kr[0] == 1:
            # [B[rd], B[md]] == 0
            return self.monomial(lhs * B[kr]) * self.monomial(rhs // B[kr])

        if kl[0] == 0 and kr[0] == 0:
            def a(m, p):
                if p <= (m - 1) // 2:
                    return q**(-2*(p-1)) * (1 + q**-2)
                # Assume m is even and p == m/2
                assert p == m // 2 and m % 2 == 0
                return q**(-m+2)
            if kl[1] * kr[1] > 0 or (kl[1] == 0 and kr[1] > 0):
                # Same sign
                # [B[rd+a1], B[(r+m)d+a1]]
                m = kr[1] - kl[1]
                assert m > 0
                terms = q**-2 * self.monomial(B[kr] * B[kl])
                terms -= self.monomial(B[1,m])
                temp = (-sum(q**(-2*(p-1)) * self.monomial(B[1,m-2*p])
                             for p in range(1, (m - 1) // 2 + 1))
                         + sum(a(m,p) * self.monomial(B[0,kr[1]-p]) * self.monomial(B[0,p+kl[1]])
                               for p in range(1, m // 2 + 1)))
                terms += (q**-2 - 1) * temp
            else:
                r = -kr[1] - 1
                # s = kl[1]
                if r <= kl[1]:
                    # [B[rd+a0], B[sd+a1]] r <= s
                    terms = -self.monomial(B[kr] * B[kl])
                    terms -= self.monomial(B[1,r+kl[1]+1])
                    terms -= (q**2-1) * sum(q**(2*k) * self.monomial(B[1,r+kl[1]-1-2*k])
                                            for k in range(r))
                    terms -= (q**2-q**-2) * sum(q**(2*(r-1-k)) * self.monomial(B[0,-(k+1)]) * self.monomial(B[0,-r+kl[1]+k])
                                                for k in range(r))
                    m = -r + kl[1] + 1
                    temp = (-sum(q**(-2*(p-1)) * self.monomial(B[1,m-2*p])
                                 for p in range(1, (m - 1) // 2 + 1))
                             + sum(a(m,p) * self.monomial(B[0,m-p-1]) * self.monomial(B[0,p-1])
                                   for p in range(1, m // 2 + 1)))
                    terms += (q**-2 - 1) * q**(2*r) * temp
                else:
                    # [B[rd+a0], B[sd+a1]] r > s
                    terms = -self.monomial(B[kr] * B[kl])
                    terms -= self.monomial(B[1,r+kl[1]+1])
                    terms -= (q**2-1) * sum(q**(2*k) * self.monomial(B[1,r+kl[1]-1-2*k])
                                            for k in range(kl[1]))
                    terms -= (q**2-q**-2) * sum(q**(2*(kl[1]-1-k)) * self.monomial(B[0,-(r-kl[1]+k+1)]) * self.monomial(B[0,k])
                                                for k in range(kl[1]))
                    m = r - kl[1] + 1
                    temp = (-sum(q**(-2*(p-1)) * self.monomial(B[1,m-2*p])
                                 for p in range(1, (m - 1) // 2 + 1))
                             + sum(a(m,p) * self.monomial(B[0,-p]) * self.monomial(B[0,p-m])
                                   for p in range(1, m // 2 + 1)))
                    terms += (q**-2 - 1) * q**(2*kl[1]) * temp
                terms = -q**2 * terms
        elif kl[0] == 1 and kr[0] == 0:
            terms = self.monomial(B[kr] * B[kl])
            # We have kr[1] < 0
            assert kr[1] < 0
            p = -kr[1] - 1
            if p < kl[1]:
                # [B[md], B[pd+a0]] with p < m
                # m = kl[1]
                terms += self._c * self._q_two * (
                           q**(-2*(kl[1]-1)) * self.monomial(B[0,-(kl[1]+p+1)])
                         + (q**2 - q**-2) * sum(q**(-2*(kl[1]-2*p+2*h))
                                                * self.monomial(B[0,-(kl[1]-p+2*h+1)])
                                                for h in range(p))
                         - q**(-2*(kl[1]-2*p-1)) * self.monomial(B[0,kl[1]-p-1])
                         )
                terms -= (q**2 - q**-2) * sum(
                           q**(-2*(ell-1)) * self.monomial(B[0,-(ell+p+1)] * B[1,kl[1]-ell])
                         + (q**2 - q**-2) * sum(q**(-2*(ell-2*h)) * self.monomial(B[0,-(ell+p-2*h+1)] * B[1,kl[1]-ell])
                                                for h in range(1, ell))
                         - q**(2*(ell-1)) * self.monomial(B[0,-(p-ell+1)] * B[1,kl[1]-ell])
                         for ell in range(1, p+1))
                terms -= (q**2 - q**-2) * sum(
                           q**(-2*(ell-1)) * self.monomial(B[0,-(ell+p+1)] * B[1,kl[1]-ell])
                         + (q**2 - q**-2) * sum(q**(-2*(ell-2*h)) * self.monomial(B[0,-(ell+p-2*h+1)] * B[1,kl[1]-ell])
                                                for h in range(1, p+1))
                         for ell in range(p+1, kl[1]))
                terms += (q**2 - q**-2) * sum(
                           q**(-2*(ell-2*p-1)) * self.monomial(B[1,kl[1]-ell] * B[0,ell-p-1])
                         for ell in range(p+1, kl[1]))
            else:
                # [B[md], B[pd+a0]] with p >= m
                # m = kl[1]
                terms += self._c * self._q_two * (
                           q**(-2*(kl[1]-1)) * self.monomial(B[0,-(p+kl[1]+1)])
                         + (q**2 - q**-2) * sum(q**(2*(kl[1]-2-2*h))
                                                * self.monomial(B[0,-(p-kl[1]+2+2*h+1)])
                                                for h in range(kl[1]-1))
                         - q**(2*(kl[1]-1)) * self.monomial(B[0,-(p-kl[1]+1)])
                         )
                terms -= (q**2 - q**-2) * sum(
                           q**(-2*(ell-1)) * self.monomial(B[0,-(p+ell+1)] * B[1,kl[1]-ell])
                         + (q**2 - q**-2) * sum(q**(-2*(ell-2*h)) * self.monomial(B[0,-(p+ell-2*h+1)] * B[1,kl[1]-ell])
                                                for h in range(1, ell))
                         - q**(2*(ell-1)) * self.monomial(B[0,-(p-ell+1)] * B[1,kl[1]-ell])
                         for ell in range(1, kl[1]))
        else:  # kl[0] == 0 and kr[0] == 1:
            terms = self.monomial(B[kr] * B[kl])
            if kl[1] < kr[1]:
                # [B[pd+a1], B[md]] with p < m
                # p = kl[1], m = kr[1]
                terms += self._c * self._q_two * (
                           q**(-2*(kr[1]-1)) * self.monomial(B[0,kr[1]+kl[1]])
                         + (q**2 - q**-2) * sum(q**(-2*(kr[1]-2*kl[1]+2*h))
                                                * self.monomial(B[0,kr[1]-kl[1]+2*h])
                                                for h in range(kl[1]))
                         - q**(-2*(kr[1]-2*kl[1]-1)) * self.monomial(B[0,kl[1]-kr[1]])
                         )
                terms -= (q**2 - q**-2) * sum(
                           q**(-2*(ell-1)) * self.monomial(B[1,kr[1]-ell] * B[0,ell+kl[1]])
                         + (q**2 - q**-2) * sum(q**(-2*(ell-2*h)) * self.monomial(B[1,kr[1]-ell] * B[0,ell+kl[1]-2*h])
                                                for h in range(1, ell))
                         - q**(2*(ell-1)) * self.monomial(B[1,kr[1]-ell] * B[0,kl[1]-ell])
                         for ell in range(1, kl[1]+1))
                terms -= (q**2 - q**-2) * sum(
                           q**(-2*(ell-1)) * self.monomial(B[1,kr[1]-ell] * B[0,ell+kl[1]])
                         + (q**2 - q**-2) * sum(q**(-2*(ell-2*h)) * self.monomial(B[1,kr[1]-ell] * B[0,ell+kl[1]-2*h])
                                                for h in range(1, kl[1]+1))
                         for ell in range(kl[1]+1, kr[1]))
                terms += (q**2 - q**-2) * sum(
                           q**(-2*(ell-2*kl[1]-1)) * self.monomial(B[0,kl[1]-ell] * B[1,kr[1]-ell])
                         for ell in range(kl[1]+1, kr[1]))
            else:
                # [B[pd+a1], B[md]] with p >= m
                # p = kl[1], m = kr[1]
                terms += self._c * self._q_two * (
                           q**(-2*(kr[1]-1)) * self.monomial(B[0,kl[1]+kr[1]])
                         + (q**2 - q**-2) * sum(q**(2*(kr[1]-2-2*h))
                                                * self.monomial(B[0,kl[1]-kr[1]+2+2*h])
                                                for h in range(kr[1]-1))
                         - q**(2*(kr[1]-1)) * self.monomial(B[0,kl[1]-kr[1]])
                         )
                terms -= (q**2 - q**-2) * sum(
                           q**(-2*(ell-1)) * self.monomial(B[1,kr[1]-ell] * B[0,kl[1]+ell])
                         + (q**2 - q**-2) * sum(q**(-2*(ell-2*h)) * self.monomial(B[1,kr[1]-ell] * B[0,kl[1]+ell-2*h])
                                                for h in range(1, ell))
                         - q**(2*(ell-1)) * self.monomial(B[1,kr[1]-ell] * B[0,kl[1]-ell])
                         for ell in range(1, kr[1]))

        return self.monomial(lhs // B[kl]) * terms * self.monomial(rhs // B[kr])

#####################################################################
# ACE of the Onsager algebra


class OnsagerAlgebraACE(InfinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    The alternating central extension of the Onsager algebra.

    The *alternating central extension* of the :class:`Onsager algebra
    <sage.algebras.lie_algebras.onsager.OnsagerAlgebra>` is the Lie algebra
    with basis elements `\{\mathcal{A}_k, \mathcal{B}_k\}_{k \in \ZZ}`
    that satisfy the relations

    .. MATH::

        \begin{aligned}
        [\mathcal{A}_k, \mathcal{A}_m] & = \mathcal{B}_{k-m} - \mathcal{B}_{m-k},
        \\ [\mathcal{A}_k, \mathcal{B}_m] & = \mathcal{A}_{k+m} - \mathcal{A}_{k-m},
        \\ [\mathcal{B}_k, \mathcal{B}_m] & = 0.
        \end{aligned}

    This has a natural injection from the Onsager algebra by the map `\iota`
    defined by

    .. MATH::

        \iota(A_k) = \mathcal{A}_k,
        \qquad\qquad
        \iota(B_k) = \mathcal{B}_k - \mathcal{B}_{-k}.

    Note that the map `\iota` differs slightly from Lemma 4.18 in [Ter2021b]_
    due to our choice of basis of the Onsager algebra.

    .. WARNING::

        We have added an extra basis vector `\mathcal{B}_0`, which would
        be `0` in the definition given in [Ter2021b]_.

    EXAMPLES:

    We begin by constructing the ACE and doing some sample computations::

        sage: O = lie_algebras.OnsagerAlgebra(QQ)
        sage: ACE = O.alternating_central_extension()
        sage: ACE
        Alternating central extension of the Onsager algebra over Rational Field

        sage: B = ACE.basis()
        sage: A1, A2, Am2 = B[0,1], B[0,2], B[0,-2]
        sage: B1, B2, Bm2 = B[1,1], B[1,2], B[1,-2]
        sage: A1.bracket(Am2)
        -B[-3] + B[3]
        sage: A1.bracket(A2)
        B[-1] - B[1]
        sage: A1.bracket(B2)
        -A[-1] + A[3]
        sage: A1.bracket(Bm2)
        A[-1] - A[3]
        sage: B2.bracket(B1)
        0
        sage: Bm2.bracket(B2)
        0
        sage: (A2 + Am2).bracket(B1 + A2 + B2 + Bm2)
        -A[-3] + A[-1] - A[1] + A[3] + B[-4] - B[4]

    The natural inclusion map `\iota` is implemented as a coercion map::

        sage: iota = ACE.coerce_map_from(O)
        sage: b = O.basis()
        sage: am1, a2, b4 = b[0,-1], b[0,2], b[1,4]
        sage: iota(am1.bracket(a2)) == iota(am1).bracket(iota(a2))
        True
        sage: iota(am1.bracket(b4)) == iota(am1).bracket(iota(b4))
        True
        sage: iota(b4.bracket(a2)) == iota(b4).bracket(iota(a2))
        True

        sage: am1 + B2
        A[-1] + B[2]
        sage: am1.bracket(B2)
        -A[-3] + A[1]
        sage: Bm2.bracket(a2)
        -A[0] + A[4]

    We have the projection map `\rho` from Lemma 4.19 in [Ter2021b]_:

    .. MATH::

        \rho(\mathcal{A}_k) = A_k,
        \qquad\qquad
        \rho(\mathcal{B}_k) = \mathrm{sgn}(k) B_{|k|}.

    The kernel of `\rho` is the center `\mathcal{Z}`, which has a basis
    `\{B_k + B_{-k}\}_{k \in \ZZ}`::

        sage: rho = ACE.projection()
        sage: rho(A1)
        A[1]
        sage: rho(Am2)
        A[-2]
        sage: rho(B1)
        1/2*G[1]
        sage: rho(Bm2)
        -1/2*G[2]
        sage: all(rho(B[1,k] + B[1,-k]) == 0 for k in range(-6,6))
        True
        sage: all(B[0,m].bracket(B[1,k] + B[1,-k]) == 0
        ....:     for k in range(-4,4) for m in range(-4,4))
        True
    """
    def __init__(self, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: ACE = lie_algebras.AlternatingCentralExtensionOnsagerAlgebra(QQ)
            sage: TestSuite(ACE).run()

            sage: B = ACE.basis()
            sage: A1, A2, Am2 = B[0,1], B[0,2], B[0,-2]
            sage: B1, B2, Bm2 = B[1,1], B[1,2], B[1,-2]
            sage: TestSuite(ACE).run(elements=[A1,A2,Am2,B1,B2,Bm2,ACE.an_element()])
        """
        cat = LieAlgebras(R).WithBasis()
        from sage.rings.integer_ring import ZZ
        from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
        I = DisjointUnionEnumeratedSets([ZZ, ZZ], keepkey=True, facade=True)
        IndexedGenerators.__init__(self, I)
        InfinitelyGeneratedLieAlgebra.__init__(self, R, index_set=I, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            Alternating central extension of the Onsager algebra over Rational Field
        """
        return "Alternating central extension of the Onsager algebra over {}".format(self.base_ring())

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: latex(O)
            \mathcal{O}_{\Bold{Q}}
        """
        from sage.misc.latex import latex
        return "\\mathcal{{O}}_{{{}}}".format(latex(self.base_ring()))

    def _repr_generator(self, m):
        """
        Return a string representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: O._repr_generator((0,-2))
            'A[-2]'
            sage: O._repr_generator((1,4))
            'B[4]'
        """
        if m[0] == 0:
            return 'A[{}]'.format(m[1])
        return 'B[{}]'.format(m[1])

    def _latex_generator(self, m):
        r"""
        Return a LaTeX representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: O._latex_generator((0,-2))
            '\\mathcal{A}_{-2}'
            sage: O._latex_generator((1,4))
            '\\mathcal{B}_{4}'
        """
        if m[0] == 0:
            return '\\mathcal{{A}}_{{{}}}'.format(m[1])
        return '\\mathcal{{B}}_{{{}}}'.format(m[1])

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = _repr_generator
    _latex_term = _latex_generator

    @cached_method
    def basis(self):
        r"""
        Return the basis of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: O.basis()
            Lazy family (Onsager ACE monomial(i))_{i in
             Disjoint union of Family (Integer Ring, Integer Ring)}
        """
        return Family(self._indices, self.monomial, name='Onsager ACE monomial')

    @cached_method
    def lie_algebra_generators(self):
        r"""
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: O.lie_algebra_generators()
            Lazy family (Onsager ACE monomial(i))_{i in
             Disjoint union of Family (Integer Ring, Integer Ring)}
        """
        return self.basis()

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: O.an_element()
            -2*A[-3] + A[2] + B[-1] + 3*B[2]
        """
        B = self.basis()
        return B[0,2] - 2*B[0,-3] + 3*B[1,2] + B[1,-1]

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: O.some_elements()
            [A[0], A[2], A[-1], B[4], B[-3], -2*A[-3] + A[2] + B[-1] + 3*B[2]]
        """
        B = self.basis()
        return [B[0,0], B[0,2], B[0,-1], B[1,4], B[1,-3], self.an_element()]

    def bracket_on_basis(self, x, y):
        r"""
        Return the bracket of basis elements indexed by ``x`` and ``y``
        where ``x < y``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ).alternating_central_extension()
            sage: O.bracket_on_basis((1,3), (1,9))  # [B, B]
            0
            sage: O.bracket_on_basis((0,8), (1,13))  # [A, B]
            -A[-5] + A[21]
            sage: O.bracket_on_basis((0,-9), (0, 7))  # [A, A]
            B[-16] - B[16]
        """
        if x[0] == 1:
            # From < property, we have y[0] == 1
            # Therefore, we have [B_k, B_m] = 0
            return self.zero()
        R = self.base_ring()
        one = R.one()
        if y[0] == 1: # [A_k, B_m] = A_{k+m} - A_{k-m}
            if y[1] == 0:  # special case for m = 0, as A_k - A_k = 0
                return self.zero()
            d = {(0, x[1]-y[1]): -one, (0, y[1]+x[1]): one}
        else:
            # [A_k, A_m] = B_{k-m} - B_{m-k}
            d = {(1, x[1]-y[1]): one, (1, y[1]-x[1]): -one}
        return self.element_class(self, d)

    def _coerce_map_from_(self, R):
        r"""
        Return if there is a coercion map from ``R``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: ACE = O.alternating_central_extension()
            sage: ACE.has_coerce_map_from(O)  # indirect doctest
            True
        """
        if isinstance(R, OnsagerAlgebra):
            if R.base_ring().has_coerce_map_from(self.base_ring()):
                return R.module_morphism(self._from_onsager_on_basis, codomain=self)
        return super()._coerce_map_from_(R)

    def _from_onsager_on_basis(self, x):
        r"""
        Map the basis element indexed by ``x`` from the corresponding
        Onsager algebra to ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: ACE = O.alternating_central_extension()
            sage: ACE._from_onsager_on_basis((0, 2))
            A[2]
            sage: ACE._from_onsager_on_basis((1, 4))
            -B[-4] + B[4]

            sage: phi = ACE.coerce_map_from(O)
            sage: a1 = O.basis()[0,1]
            sage: a3 = O.basis()[0,3]
            sage: b2 = O.basis()[1,2]
            sage: phi(a3)
            A[3]
            sage: phi(b2)
            -B[-2] + B[2]
            sage: b2.bracket(a3)
            2*A[1] - 2*A[5]
            sage: phi(b2).bracket(phi(a3))
            2*A[1] - 2*A[5]
            sage: phi(b2.bracket(a3))
            2*A[1] - 2*A[5]

            sage: a1.bracket(a3)
            -G[2]
            sage: phi(a1).bracket(phi(a3))
            B[-2] - B[2]
            sage: phi(a1.bracket(a3))
            B[-2] - B[2]
        """
        one = self.base_ring().one()
        if x[0] == 0:
            return self._from_dict({x: one}, remove_zeros=False)
        return self._from_dict({(1, x[1]): one, (1, -x[1]): -one}, remove_zeros=False)

    def projection(self):
        r"""
        Return the projection map `\rho` from Lemma 4.19 in [Ter2021b]_
        to the Onsager algebra.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: ACE = O.alternating_central_extension()
            sage: rho = ACE.projection()
            sage: B = ACE.basis()
            sage: A1, A2, Am2 = B[0,1], B[0,2], B[0,-2]
            sage: B1, B2, Bm2 = B[1,1], B[1,2], B[1,-2]

            sage: rho(A1)
            A[1]
            sage: rho(Am2)
            A[-2]
            sage: rho(B1)
            1/2*G[1]
            sage: rho(B2)
            1/2*G[2]
            sage: rho(Bm2)
            -1/2*G[2]

            sage: rho(A1.bracket(A2))
            -G[1]
            sage: rho(A1).bracket(rho(A2))
            -G[1]
            sage: rho(B1.bracket(Am2))
            A[-3] - A[-1]
            sage: rho(B1).bracket(rho(Am2))
            A[-3] - A[-1]
        """
        O = OnsagerAlgebra(self.base_ring())
        return self.module_morphism(self._projection_on_basis, codomain=O)

    def _projection_on_basis(self, x):
        r"""
        Compute the projection map `\rho` on the basis element ``x``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: ACE = O.alternating_central_extension()
            sage: ACE._projection_on_basis((0,2))
            A[2]
            sage: ACE._projection_on_basis((1,4))
            1/2*G[4]
            sage: ACE._projection_on_basis((1,-4))
            -1/2*G[4]
        """
        R = self.base_ring()
        O = OnsagerAlgebra(R)
        if x[0] == 0: # A_k
            return O._from_dict({x: R.one()}, remove_zeros=False)
        # Otherwise B_k
        c = R.one() / 2
        if x[1] < 0:
            return O._from_dict({(1, -x[1]): -c}, remove_zeros=False)
        elif x[1] == 0:
            return O.zero()
        else:
            return O._from_dict({x: c}, remove_zeros=False)

    Element = LieAlgebraElement
