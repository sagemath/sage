"""
Alternating Central Extension Quantum Onsager Algebra

AUTHORS:

- Travis Scrimshaw (2021-03): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.sets.positive_integers import PositiveIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.rings.integer_ring import ZZ


class ACEQuantumOnsagerAlgebra(CombinatorialFreeModule):
    r"""
    The alternating central extension of the `q`-Onsager algebra.

    The *alternating central extension* `\mathcal{A}_q` of the `q`-Onsager
    algebra `O_q` is a current algebra of `O_q` introduced by Baseilhac
    and Koizumi [BK2005]_. A presentation was given by Baseilhac
    and Shigechi [BS2010]_, which was then reformulated in terms of currents
    in [Ter2021]_ and then used to prove that the generators form a PBW basis.

    .. NOTE::

        This is only for the `q`-Onsager algebra with parameter
        `c = q^{-1} (q - q^{-1})^2`.

    EXAMPLES::

        sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
        sage: AG = A.algebra_generators()

    We construct the generators `\mathcal{G}_3`, `\mathcal{W}_{-5}`,
    `\mathcal{W}_2`, and `\widetilde{\mathcal{G}}_{4}` and perform
    some computations::

        sage: G3 = AG[0,3]
        sage: Wm5 = AG[1,-5]
        sage: W2 = AG[1,2]
        sage: Gt4 = AG[2,4]
        sage: [G3, Wm5, W2, Gt4]
        [G[3], W[-5], W[2], Gt[4]]
        sage: Gt4 * G3
        G[3]*Gt[4] + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[-6]*W[1]
         + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[-5]*W[2]
         + ((q^12-3*q^8+3*q^4-1)/q^6)*W[-4]*W[1]
         + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[-4]*W[3]
         + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[-3]*W[-2]
         + ((q^12-3*q^8+3*q^4-1)/q^6)*W[-3]*W[2]
         + ((q^12-3*q^8+3*q^4-1)/q^6)*W[-2]*W[5]
         + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[-1]*W[4]
         + ((q^12-3*q^8+3*q^4-1)/q^6)*W[-1]*W[6]
         + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[0]*W[5]
         + ((q^12-3*q^8+3*q^4-1)/q^6)*W[0]*W[7]
         + ((q^12-3*q^8+3*q^4-1)/q^6)*W[3]*W[4]
        sage: Wm5 * G3
        ((q^2-1)/q^2)*G[1]*W[-7] + ((-q^2+1)/q^2)*G[1]*W[7]
         + ((q^2-1)/q^2)*G[2]*W[-6] + ((-q^2+1)/q^2)*G[2]*W[6] + G[3]*W[-5]
         + ((-q^2+1)/q^2)*G[6]*W[-2] + ((q^2-1)/q^2)*G[6]*W[2]
         + ((-q^2+1)/q^2)*G[7]*W[-1] + ((q^2-1)/q^2)*G[7]*W[1]
         + ((-q^2+1)/q^2)*G[8]*W[0] + ((-q^8+2*q^4-1)/q^5)*W[-8]
         + ((q^8-2*q^4+1)/q^5)*W[8]
        sage: W2 * G3
        (q^2-1)*G[1]*W[-2] + (-q^2+1)*G[1]*W[4] + (-q^2+1)*G[3]*W[0]
         + q^2*G[3]*W[2] + (q^2-1)*G[4]*W[1] + ((-q^8+2*q^4-1)/q^3)*W[-3]
         + ((q^8-2*q^4+1)/q^3)*W[5]
        sage: W2 * Wm5
        (q^4/(q^8+2*q^6-2*q^2-1))*G[1]*Gt[6] + (-q^4/(q^8+2*q^6-2*q^2-1))*G[6]*Gt[1]
         + W[-5]*W[2] + (q/(q^2+1))*G[7] + (-q/(q^2+1))*Gt[7]
        sage: Gt4 * Wm5
        ((q^2-1)/q^2)*W[-8]*Gt[1] + ((q^2-1)/q^2)*W[-7]*Gt[2]
         + ((q^2-1)/q^2)*W[-6]*Gt[3] + W[-5]*Gt[4] + ((-q^2+1)/q^2)*W[-3]*Gt[6]
         + ((-q^2+1)/q^2)*W[-2]*Gt[7] + ((-q^2+1)/q^2)*W[-1]*Gt[8]
         + ((-q^2+1)/q^2)*W[0]*Gt[9] + ((q^2-1)/q^2)*W[1]*Gt[8]
         + ((q^2-1)/q^2)*W[2]*Gt[7] + ((q^2-1)/q^2)*W[3]*Gt[6]
         + ((-q^2+1)/q^2)*W[6]*Gt[3] + ((-q^2+1)/q^2)*W[7]*Gt[2]
         + ((-q^2+1)/q^2)*W[8]*Gt[1] + ((-q^8+2*q^4-1)/q^5)*W[-9]
         + ((q^8-2*q^4+1)/q^5)*W[9]
        sage: Gt4 * W2
        (q^2-1)*W[-3]*Gt[1] + (-q^2+1)*W[0]*Gt[4] + (q^2-1)*W[1]*Gt[5]
         + q^2*W[2]*Gt[4] + (-q^2+1)*W[5]*Gt[1] + ((-q^8+2*q^4-1)/q^3)*W[-4]
         + ((q^8-2*q^4+1)/q^3)*W[6]

    REFERENCES:

    - [BK2005]_
    - [BS2010]_
    - [Ter2021]_
    """
    @staticmethod
    def __classcall_private__(cls, R=None, q=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: A1 = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: R.<q> = QQ[]
            sage: q = R.gen()
            sage: A2 = algebras.AlternatingCentralExtensionQuantumOnsager(R.fraction_field(), q)
            sage: A1 is A2
            True
            sage: A2.q().parent() is R.fraction_field()
            True
            sage: q = R.fraction_field().gen()
            sage: A3 = algebras.AlternatingCentralExtensionQuantumOnsager(q=q)
            sage: A1 is A3
            True
        """
        if q is None:
            if R is None:
                raise ValueError("either base ring or q must be specified")
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            q = PolynomialRing(R, 'q').fraction_field().gen()
            R = q.parent()
        else:
            if R is None:
                R = q.parent()
            else:
                q = R(q)
        return super(ACEQuantumOnsagerAlgebra, cls).__classcall__(cls, R, q)

    def __init__(self, R, q):
        r"""
        Initialize ``self``.

        TESTS::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: TestSuite(A).run()  # long time
        """
        I = DisjointUnionEnumeratedSets([PositiveIntegers(), ZZ, PositiveIntegers()],
                                        keepkey=True, facade=True)
        monomials = IndexedFreeAbelianMonoid(I, prefix='A', bracket=False)
        self._q = q
        CombinatorialFreeModule.__init__(self, R, monomials,
                                         prefix='', bracket=False, latex_bracket=False,
                                         sorting_key=self._monomial_key,
                                         category=Algebras(R).WithBasis().Filtered())

    def _monomial_key(self, x):
        r"""
        Compute the key for ``x`` so that the comparison is done by
        reverse degree lexicographic order.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: AG = A.algebra_generators()
            sage: AG[1,1] * AG[1,0] * AG[0,1]               # indirect doctest
            G[1]*W[0]*W[1] + (q/(q^2+1))*G[1]^2 + (-q/(q^2+1))*G[1]*Gt[1]
             + ((-q^8+2*q^4-1)/q^5)*W[-1]*W[1] + ((-q^8+2*q^4-1)/q^5)*W[0]^2
             + ((q^8-2*q^4+1)/q^5)*W[0]*W[2] + ((q^8-2*q^4+1)/q^5)*W[1]^2
        """
        return (-len(x), x.to_word_list())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: A
            Alternating Central Extension of q-Onsager algebra over Fraction
             Field of Univariate Polynomial Ring in q over Rational Field
        """
        return "Alternating Central Extension of {}-Onsager algebra over {}".format(
                                                     self._q, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: latex(A)
            \mathcal{A}_{q,\mathrm{Frac}(\Bold{Q}[q])}
        """
        from sage.misc.latex import latex
        return "\\mathcal{{A}}_{{{},{}}}".format(latex(self._q),
                                                   latex(self.base_ring()))

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: I = A._indices.gens()
            sage: A._repr_term(I[0,3])
            'G[3]'
            sage: A._repr_term(I[1,-2])
            'W[-2]'
            sage: A._repr_term(I[1,3])
            'W[3]'
            sage: A._repr_term(I[2,5])
            'Gt[5]'
            sage: A._repr_term(I[0,1]^2 * I[1,0] * I[1,3]^13 * I[2,3])
            'G[1]^2*W[0]*W[3]^13*Gt[3]'
        """
        def to_str(x):
            k, e = x
            if k[0] == 0:
                ret = "G[{}]".format(k[1])
            elif k[0] == 1:
                ret = "W[{}]".format(k[1])
            elif k[0] == 2:
                ret = "Gt[{}]".format(k[1])
            if e > 1:
                ret = ret + "^{}".format(e)
            return ret
        return '*'.join(to_str(x) for x in m._sorted_items())

    def _latex_term(self, m):
        r"""
        Return a latex representation of the term indexed by ``m``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: I = A._indices.gens()
            sage: A._latex_term(I[0,3])
            '\\mathcal{G}_{3}'
            sage: A._latex_term(I[1,-2])
            '\\mathcal{W}_{-2}'
            sage: A._latex_term(I[1,3])
            '\\mathcal{W}_{3}'
            sage: A._latex_term(I[2,5])
            '\\widetilde{\\mathcal{G}}_{5}'
            sage: A._latex_term(I[0,1]^2 * I[1,0] * I[1,3]^13 * I[2,3])
            '\\mathcal{G}_{1}^{2} \\mathcal{W}_{0} \\mathcal{W}_{3}^{13} \\widetilde{\\mathcal{G}}_{3}'
        """
        def to_str(x):
            k, e = x
            if k[0] == 0:
                ret = "\\mathcal{{G}}_{{{}}}".format(k[1])
            elif k[0] == 1:
                ret = "\\mathcal{{W}}_{{{}}}".format(k[1])
            elif k[0] == 2:
                ret = "\\widetilde{{\\mathcal{{G}}}}_{{{}}}".format(k[1])
            if e > 1:
                ret = ret + '^{{{}}}'.format(e)
            return ret
        return ' '.join(to_str(x) for x in m._sorted_items())

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: A.algebra_generators()
            Lazy family (generator map(i))_{i in Disjoint union of
             Family (Positive integers, Integer Ring, Positive integers)}
        """
        G = self._indices.gens()
        q = self._q

        def monomial_map(x):
            if x[0] != 1 and x[1] == 0:
                return self.term(self.one_basis(), -(q-~q)*(q+~q)**2)
            return self.monomial(G[x])
        return Family(self._indices._indices, monomial_map,
                      name="generator map")

    gens = algebra_generators

    def q(self):
        r"""
        Return the parameter `q` of ``self``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: A.q()
            q
        """
        return self._q

    @cached_method
    def one_basis(self):
        r"""
        Return the basis element indexing `1`.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: ob = A.one_basis(); ob
            1
            sage: ob.parent()
            Free abelian monoid indexed by Disjoint union of
             Family (Positive integers, Integer Ring, Positive integers)
        """
        return self._indices.one()

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: A.an_element()
            q*G[2] - 2*W[-3] + W[2] - q*Gt[1]
        """
        G = self.algebra_generators()
        return G[1,2] - 2*G[1,-3] + self.base_ring().an_element()*(G[0,2] - G[2,1])

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: A.some_elements()
            [W[0], W[3], W[-1], W[1], W[-2], G[1], G[2], Gt[1], Gt[2]]
        """
        G = self.algebra_generators()
        return [G[1,0], G[1,3], G[1,-1], G[1,1], G[1,-2], G[0,1], G[0,2], G[2,1], G[2,2]]

    def degree_on_basis(self, m):
        r"""
        Return the degree of the basis element indexed by ``m``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: G = A.algebra_generators()
            sage: A.degree_on_basis(G[0,1].leading_support())
            2
            sage: A.degree_on_basis(G[0,2].leading_support())
            4
            sage: A.degree_on_basis(G[1,-1].leading_support())
            3
            sage: A.degree_on_basis(G[1,0].leading_support())
            1
            sage: A.degree_on_basis(G[1,1].leading_support())
            1
            sage: A.degree_on_basis(G[2,1].leading_support())
            2
            sage: A.degree_on_basis(G[2,2].leading_support())
            4
            sage: [x.degree() for x in A.some_elements()]
            [1, 5, 3, 1, 5, 2, 4, 2, 4]
        """
        def deg(k):
            if k[0] != 1:
                return 2*k[1]
            return -2*k[1]+1 if k[1] <= 0 else 2*k[1] - 1
        return ZZ.sum(deg(k) * c for k, c in m._monomial.items())

    @cached_method
    def quantum_onsager_pbw_generator(self, i):
        r"""
        Return the image of the PBW generator of the `q`-Onsager algebra
        in ``self``.

        INPUT:

        - ``i`` -- a pair ``(k, m)`` such that

          * ``k=0`` and ``m`` is an integer
          * ``k=1`` and ``m`` is a positive integer

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: A.quantum_onsager_pbw_generator((0,0))
            W[1]
            sage: A.quantum_onsager_pbw_generator((0,1))
            (q^3/(q^4-1))*W[1]*Gt[1] - q^2*W[0] + (q^2+1)*W[2]
            sage: A.quantum_onsager_pbw_generator((0,2))
            (q^6/(q^8-2*q^4+1))*W[1]*Gt[1]^2 + (-q^5/(q^4-1))*W[0]*Gt[1]
             + (q^3/(q^2-1))*W[1]*Gt[2] + (q^3/(q^2-1))*W[2]*Gt[1]
             + (-q^4-q^2)*W[-1] - q^2*W[1] + (q^4+2*q^2+1)*W[3]
            sage: A.quantum_onsager_pbw_generator((0,-1))
            W[0]
            sage: A.quantum_onsager_pbw_generator((0,-2))
            (q/(q^4-1))*W[0]*Gt[1] + ((q^2+1)/q^2)*W[-1] - 1/q^2*W[1]
            sage: A.quantum_onsager_pbw_generator((0,-3))
            (q^2/(q^8-2*q^4+1))*W[0]*Gt[1]^2 + (1/(q^3-q))*W[-1]*Gt[1]
             + (1/(q^3-q))*W[0]*Gt[2] - (1/(q^5-q))*W[1]*Gt[1]
             + ((q^4+2*q^2+1)/q^4)*W[-2] - 1/q^2*W[0] + ((-q^2-1)/q^4)*W[2]
            sage: A.quantum_onsager_pbw_generator((1,1))
            ((-q^2+1)/q^2)*W[0]*W[1] + (1/(q^3+q))*G[1] - (1/(q^3+q))*Gt[1]
            sage: A.quantum_onsager_pbw_generator((1,2))
            -1/q*W[0]*W[1]*Gt[1] + (1/(q^6+q^4-q^2-1))*G[1]*Gt[1]
             + ((-q^4+1)/q^4)*W[-1]*W[1] + (q^2-1)*W[0]^2
             + ((-q^4+1)/q^2)*W[0]*W[2] + ((q^2-1)/q^4)*W[1]^2
             - (1/(q^6+q^4-q^2-1))*Gt[1]^2 + 1/q^3*G[2] - 1/q^3*Gt[2]
        """
        W0 = self.algebra_generators()[1,0]
        W1 = self.algebra_generators()[1,1]
        q = self._q
        if i[0] == 0:
            if i[1] < 0:
                if i[1] == -1:
                    return W0
                Bd = self.quantum_onsager_pbw_generator((1, 1))
                Bm1 = self.quantum_onsager_pbw_generator((0, i[1]+1))
                Bm2 = self.quantum_onsager_pbw_generator((0, i[1]+2))
                return Bm2 + q/(q**-3-~q-q+q**3) * (Bd * Bm1 - Bm1 * Bd)
            else:
                if i[1] == 0:
                    return W1
                Bd = self.quantum_onsager_pbw_generator((1, 1))
                Bm1 = self.quantum_onsager_pbw_generator((0, i[1]-1))
                Bm2 = self.quantum_onsager_pbw_generator((0, i[1]-2))
                return Bm2 - q/(q**-3-~q-q+q**3) * (Bd * Bm1 - Bm1 * Bd)
        elif i[0] == 1:
            if i[1] == 1:
                return q**-2 * W1 * W0 - W0 * W1
            if i[1] <= 0:
                raise ValueError("not an index of a PBW basis element")
            B = self.quantum_onsager_pbw_generator
            n = i[1]
            Bm1 = self.quantum_onsager_pbw_generator((0, n-1))
            return (q**-2 * Bm1 * W0 - W0 * Bm1
                    + (q**-2 - 1) * sum(B((0,ell)) * B((0,n-ell-2))
                                        for ell in range(n-1)))
        raise ValueError("not an index of a PBW basis element")

    @cached_method
    def product_on_basis(self, lhs, rhs):
        r"""
        Return the product of the two basis elements ``lhs`` and ``rhs``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: G = A.algebra_generators()
            sage: q = A.q()
            sage: rho = -(q^2 - q^-2)^2

        We verify the PBW ordering::

            sage: G[0,1] * G[1,1]  # indirect doctest
            G[1]*W[1]
            sage: G[1,1] * G[0,1]
            q^2*G[1]*W[1] + ((-q^8+2*q^4-1)/q^3)*W[0] + ((q^8-2*q^4+1)/q^3)*W[2]
            sage: G[1,-1] * G[1,1]
            W[-1]*W[1]
            sage: G[1,1] * G[1,-1]
            W[-1]*W[1] + (q/(q^2+1))*G[2] + (-q/(q^2+1))*Gt[2]
            sage: G[1,1] * G[2,1]
            W[1]*Gt[1]
            sage: G[2,1] * G[1,1]
            q^2*W[1]*Gt[1] + ((-q^8+2*q^4-1)/q^3)*W[0] + ((q^8-2*q^4+1)/q^3)*W[2]
            sage: G[0,1] * G[2,1]
            G[1]*Gt[1]
            sage: G[2,1] * G[0,1]
            G[1]*Gt[1] + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[-1]*W[1]
             + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[0]^2
             + ((q^12-3*q^8+3*q^4-1)/q^6)*W[0]*W[2]
             + ((q^12-3*q^8+3*q^4-1)/q^6)*W[1]^2

        We verify some of the defining relations (see Equations (3-14)
        in [Ter2021]_), which are used to construct the PBW basis::

            sage: G[0,1] * G[0,2] == G[0,2] * G[0,1]
            True
            sage: G[1,-1] * G[1,-2] == G[1,-2] * G[1,-1]
            True
            sage: G[1,1] * G[1,2] == G[1,2] * G[1,1]
            True
            sage: G[2,1] * G[2,2] == G[2,2] * G[2,1]
            True
            sage: G[1,0] * G[1,2] - G[1,2] * G[1,0] == G[1,-1] * G[1,1] - G[1,1] * G[1,-1]
            True
            sage: G[1,0] * G[1,2] - G[1,2] * G[1,0] == (G[2,2] - G[0,2]) / (q + ~q)
            True
            sage: q * G[1,0] * G[0,2] - ~q * G[0,2] * G[1,0] == q * G[2,2] * G[1,0] - ~q * G[1,0] * G[2,2]
            True
            sage: q * G[1,0] * G[0,2] - ~q * G[0,2] * G[1,0] == rho * G[1,-2] - rho * G[1,2]
            True
            sage: q * G[0,2] * G[1,1] - ~q * G[1,1] * G[0,2] == q * G[1,1] * G[2,2] - ~q * G[2,2] * G[1,1]
            True
            sage: q * G[0,2] * G[1,1] - ~q * G[1,1] * G[0,2] == rho * G[1,3] - rho * G[1,-1]
            True
            sage: G[1,-2] * G[1,2] - G[1,2] * G[1,-2] == G[1,-1] * G[1,3] - G[1,3] * G[1,-1]
            True
            sage: G[1,-2] * G[0,2] - G[0,2] * G[1,-2] == G[1,-1] * G[0,3] - G[0,3] * G[1,-1]
            True
            sage: G[1,1] * G[0,2] - G[0,2] * G[1,1] == G[1,2] * G[0,1] - G[0,1] * G[1,2]
            True
            sage: G[1,-2] * G[2,2] - G[2,2] * G[1,-2] == G[1,-1] * G[2,3] - G[2,3] * G[1,-1]
            True
            sage: G[1,1] * G[2,2] - G[2,2] * G[1,1] == G[1,2] * G[2,1] - G[2,1] * G[1,2]
            True
            sage: G[0,1] * G[2,2] - G[2,2] * G[0,1] == G[0,2] * G[2,1] - G[2,1] * G[0,2]
            True
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
        if kl <= kr:
            return self.monomial(lhs * rhs)

        A = self.algebra_generators()

        # Create the commutator
        # We have xy - yx = [x, y] -> xy = yx + LOT for x > y
        if kl[0] == kr[0]:
            # Commuting elements
            if kl[0] != 1 or (kl[1] > 0 and kr[1] > 0) or (kl[1] <= 0 and kr[1] <= 0):
                return self.monomial(lhs * B[kr]) * self.monomial(rhs // B[kr])

            # relation (ii)
            i = kl[1] - 1
            j = -kr[1]
            denom = (q**2 - q**-2) * (q + q**-1)**2
            terms = A[1,-j]*A[1,i+1] + self.sum(A[0,ell]*A[2,i+j+1-ell] - A[0,i+j+1-ell]*A[2,ell] for ell in range(min(i,j)+1)) / denom
        elif kl[0] == 2 and kr[0] == 0:
            # relation (iii)
            i = kl[1] - 1
            j = kr[1] - 1
            coeff = (q**2 - q**-2)**3
            terms = (A[0,j+1]*A[2,i+1] - coeff * A[1,-i]*A[1,-j] + coeff * A[1,i+1]*A[1,j+1]
                                       + coeff * sum(A[1,-ell]*A[1,i+j+2-ell] - A[1,ell-1-i-j]*A[1,ell+1] for ell in range(min(i,j)+1))
                                       - coeff * sum(A[1,1-ell]*A[1,i+j+1-ell] - A[1,ell-i-j]*A[1,ell] for ell in range(1, min(i,j)+1)))
        elif kl[0] == 1 and kr[0] == 0:
            if kl[1] > 0:
                # relation (vi)
                i = kl[1] - 1
                j = kr[1] - 1
                coeff = q * (q - ~q)
                terms = (A[0,j+1]*A[1,i+1] + coeff * sum(A[0,ell]*A[1,ell-i-j] for ell in range(min(i,j)+1))
                                           + coeff * sum(A[0,i+j+1-ell]*A[1,ell+1] - A[0,ell]*A[1,i+j+2-ell] for ell in range(min(i,j)+1))
                                           - coeff * sum(A[0,i+j+1-ell]*A[1,1-ell] for ell in range(1, min(i,j)+1)))
            else:
                # relation (v)
                i = -kl[1]
                j = kr[1] - 1
                coeff = ~q * (q - ~q)
                terms = (A[0,j+1]*A[1,-i] - coeff * sum(A[0,ell]*A[1,i+j+1-ell] for ell in range(min(i,j)+1))
                                          + coeff * sum(A[0,ell]*A[1,ell-1-i-j] - A[0,i+j+1-ell]*A[1,-ell] for ell in range(min(i,j)+1))
                                          + coeff * sum(A[0,i+j+1-ell]*A[1,ell] for ell in range(1, min(i,j)+1)))
        elif kl[0] == 2 and kr[0] == 1:
            if kr[1] > 0:
                # relation (vi)
                i = kl[1] - 1
                j = kr[1] - 1
                coeff = q * (q - ~q)
                terms = (A[1,j+1]*A[2,i+1] + coeff * sum(A[1,ell-i-j]*A[2,ell] for ell in range(min(i,j)+1))
                                           + coeff * sum(A[1,ell+1]*A[2,i+j+1-ell] - A[1,i+j+2-ell]*A[2,ell] for ell in range(min(i,j)+1))
                                           - coeff * sum(A[1,1-ell]*A[2,i+j+1-ell] for ell in range(1, min(i,j)+1)))
            else:
                # relation (vii)
                i = kl[1] - 1
                j = -kr[1]
                coeff = ~q * (q - ~q)
                terms = (A[1,-j]*A[2,i+1] - coeff * sum(A[1,i+j+1-ell]*A[2,ell] for ell in range(min(i,j)+1))
                                          + coeff * sum(A[1,ell-1-i-j]*A[2,ell] - A[1,-ell]*A[2,i+j+1-ell] for ell in range(min(i,j)+1))
                                          + coeff * sum(A[1,ell]*A[2,i+j+1-ell] for ell in range(1, min(i,j)+1)))

        return self.monomial(lhs // B[kl]) * terms * self.monomial(rhs // B[kr])

    def _sigma_on_basis(self, x):
        r"""
        Return the action of the `\sigma` automorphism on the basis element
        indexed by ``x``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: I = A._indices.monoid_generators()
            sage: A._sigma_on_basis(I[1,-1] * I[1,1])
            W[0]*W[2] + (q/(q^2+1))*G[2] + (-q/(q^2+1))*Gt[2]
            sage: A._sigma_on_basis(I[0,1] * I[2,1])
            G[1]*Gt[1] + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[-1]*W[1]
             + ((-q^12+3*q^8-3*q^4+1)/q^6)*W[0]^2
             + ((q^12-3*q^8+3*q^4-1)/q^6)*W[0]*W[2]
             + ((q^12-3*q^8+3*q^4-1)/q^6)*W[1]^2
            sage: [(x, A.sigma(x)) for x in A.some_elements()]
            [(W[0], W[1]), (W[3], W[-2]), (W[-1], W[2]), (W[1], W[0]),
             (W[-2], W[3]), (G[1], Gt[1]), (G[2], Gt[2]), (Gt[1], G[1]),
             (Gt[2], G[2])]
        """
        def tw(m):
            if m[0] == 0:
                return (2, m[1])
            if m[0] == 1:
                return (1, -m[1]+1)
            if m[0] == 2:
                return (0, m[1])
        A = self.algebra_generators()
        return self.prod(A[tw(m)]**e for m,e in x._sorted_items())

    def _dagger_on_basis(self, x):
        r"""
        Return the action of the `\dagger` antiautomorphism on the basis element
        indexed by ``x``.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: I = A._indices.monoid_generators()
            sage: A._dagger_on_basis(I[0,1] * I[1,-1] * I[2,1])
            G[1]*W[-1]*Gt[1]
            sage: A._dagger_on_basis(I[1,-1] * I[1,1])
            W[-1]*W[1] + (q/(q^2+1))*G[2] + (-q/(q^2+1))*Gt[2]
            sage: A._dagger_on_basis(I[0,1] * I[1,-1] * I[1,2] * I[2,1])
            (q^4/(q^8+2*q^6-2*q^2-1))*G[1]^2*Gt[1]*Gt[2]
             + (-q^4/(q^8+2*q^6-2*q^2-1))*G[1]*G[2]*Gt[1]^2
             + G[1]*W[-1]*W[2]*Gt[1] + (q/(q^2+1))*G[1]*G[3]*Gt[1]
             + (-q/(q^2+1))*G[1]*Gt[1]*Gt[3]
            sage: [(x, A.dagger(x)) for x in A.some_elements()]
            [(W[0], W[0]), (W[3], W[3]), (W[-1], W[-1]), (W[1], W[1]),
             (W[-2], W[-2]), (G[1], Gt[1]), (G[2], Gt[2]), (Gt[1], G[1]),
             (Gt[2], G[2])]
        """
        def tw(m):
            if m[0] == 0:
                return (2, m[1])
            if m[0] == 1:
                return (1, m[1])
            if m[0] == 2:
                return (0, m[1])
        A = self.algebra_generators()
        return self.prod(A[tw(m)]**e for m,e in reversed(x._sorted_items()))

    @lazy_attribute
    def sigma(self):
        r"""
        The automorphism `\sigma`.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: G = A.algebra_generators()
            sage: x = A.an_element()^2
            sage: A.sigma(A.sigma(x)) == x
            True
            sage: A.sigma(G[1,-1] * G[1,1]) == A.sigma(G[1,-1]) * A.sigma(G[1,1])
            True
            sage: A.sigma(G[0,2] * G[1,3]) == A.sigma(G[0,2]) * A.sigma(G[1,3])
            True
        """
        return self.module_morphism(self._sigma_on_basis, codomain=self, category=self.category())

    @lazy_attribute
    def dagger(self):
        r"""
        The antiautomorphism `\dagger`.

        EXAMPLES::

            sage: A = algebras.AlternatingCentralExtensionQuantumOnsager(QQ)
            sage: G = A.algebra_generators()
            sage: x = A.an_element()^2
            sage: A.dagger(A.dagger(x)) == x
            True
            sage: A.dagger(G[1,-1] * G[1,1]) == A.dagger(G[1,1]) * A.dagger(G[1,-1])
            True
            sage: A.dagger(G[0,2] * G[1,3]) == A.dagger(G[1,3]) * A.dagger(G[0,2])
            True
            sage: A.dagger(G[2,2] * G[1,3]) == A.dagger(G[1,3]) * A.dagger(G[2,2])
            True
        """
        return self.module_morphism(self._dagger_on_basis, codomain=self)
