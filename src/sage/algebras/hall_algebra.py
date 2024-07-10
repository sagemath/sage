# sage.doctest: needs sage.combinat sage.modules
r"""
Hall Algebras

AUTHORS:

- Travis Scrimshaw (2013-10-17): Initial version
"""
# ****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_method
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.combinat.partition import Partition, Partitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.hall_polynomial import hall_polynomial
from sage.combinat.sf.sf import SymmetricFunctions
from sage.rings.integer_ring import ZZ
from functools import cmp_to_key, reduce


def transpose_cmp(x, y):
    r"""
    Compare partitions ``x`` and ``y`` in transpose dominance order.

    We say partitions `\mu` and `\lambda` satisfy `\mu \prec \lambda`
    in transpose dominance order if for all `i \geq 1` we have:

    .. MATH::

        l_1 + 2 l_2 + \cdots + (i-1) l_{i-1} + i(l_i + l_{i+1} + \cdots) \leq
        m_1 + 2 m_2 + \cdots + (i-1) m_{i-1} + i(m_i + m_{i+1} + \cdots),

    where `l_k` denotes the number of appearances of `k` in
    `\lambda`, and `m_k` denotes the number of appearances of `k`
    in `\mu`.

    Equivalently, `\mu \prec \lambda` if the conjugate of the
    partition `\mu` dominates the conjugate of the partition
    `\lambda`.

    Since this is a partial ordering, we fallback to lex ordering
    `\mu <_L \lambda` if we cannot compare in the transpose order.

    EXAMPLES::

        sage: from sage.algebras.hall_algebra import transpose_cmp
        sage: transpose_cmp(Partition([4,3,1]), Partition([3,2,2,1]))
        -1
        sage: transpose_cmp(Partition([2,2,1]), Partition([3,2]))
        1
        sage: transpose_cmp(Partition([4,1,1]), Partition([4,1,1]))
        0
    """
    if x == y:
        return 0
    xexp = x.to_exp()
    yexp = y.to_exp()
    n = min(len(xexp), len(yexp))

    def check(m, l):
        s1 = 0
        s2 = 0
        for i in range(n):
            s1 += sum(l[i:])
            s2 += sum(m[i:])
            if s1 > s2:
                return False
        return sum(l) <= sum(m)
    if check(xexp, yexp):
        return 1
    if check(yexp, xexp):
        return -1
    if x < y:
        return -1
    return 1


class HallAlgebra(CombinatorialFreeModule):
    r"""
    The (classical) Hall algebra.

    The *(classical) Hall algebra* over a commutative ring `R` with a
    parameter `q \in R` is defined to be the free `R`-module with
    basis `(I_\lambda)`, where `\lambda` runs over all integer
    partitions. The algebra structure is given by a product defined by

    .. MATH::

        I_\mu \cdot I_\lambda = \sum_\nu P^{\nu}_{\mu, \lambda}(q) I_\nu,

    where `P^{\nu}_{\mu, \lambda}` is a Hall polynomial (see
    :func:`~sage.combinat.hall_polynomial.hall_polynomial`). The
    unity of this algebra is `I_{\emptyset}`.

    The (classical) Hall algebra is also known as the Hall-Steinitz
    algebra.

    We can define an `R`-algebra isomorphism `\Phi` from the
    `R`-algebra of symmetric functions (see
    :class:`~sage.combinat.sf.sf.SymmetricFunctions`) to the
    (classical) Hall algebra by sending the `r`-th elementary
    symmetric function `e_r` to `q^{r(r-1)/2} I_{(1^r)}` for every
    positive integer `r`. This isomorphism used to transport the
    Hopf algebra structure from the `R`-algebra of symmetric functions
    to the Hall algebra, thus making the latter a connected graded
    Hopf algebra. If `\lambda` is a partition, then the preimage
    of the basis element `I_{\lambda}` under this isomorphism is
    `q^{n(\lambda)} P_{\lambda}(x; q^{-1})`, where `P_{\lambda}` denotes
    the `\lambda`-th Hall-Littlewood `P`-function, and where
    `n(\lambda) = \sum_i (i - 1) \lambda_i`.

    See section 2.3 in [Sch2006]_, and sections II.2 and III.3
    in [Mac1995]_ (where our `I_{\lambda}` is called `u_{\lambda}`).

    EXAMPLES::

        sage: R.<q> = ZZ[]
        sage: H = HallAlgebra(R, q)
        sage: H[2,1]*H[1,1]
        H[3, 2] + (q+1)*H[3, 1, 1] + (q^2+q)*H[2, 2, 1] + (q^4+q^3+q^2)*H[2, 1, 1, 1]
        sage: H[2]*H[2,1]
        H[4, 1] + q*H[3, 2] + (q^2-1)*H[3, 1, 1] + (q^3+q^2)*H[2, 2, 1]
        sage: H[3]*H[1,1]
        H[4, 1] + q^2*H[3, 1, 1]
        sage: H[3]*H[2,1]
        H[5, 1] + q*H[4, 2] + (q^2-1)*H[4, 1, 1] + q^3*H[3, 2, 1]

    We can rewrite the Hall algebra in terms of monomials of
    the elements `I_{(1^r)}`::

        sage: I = H.monomial_basis()
        sage: H(I[2,1,1])
        H[3, 1] + (q+1)*H[2, 2] + (2*q^2+2*q+1)*H[2, 1, 1]
         + (q^5+2*q^4+3*q^3+3*q^2+2*q+1)*H[1, 1, 1, 1]
        sage: I(H[2,1,1])
        I[3, 1] + (-q^3-q^2-q-1)*I[4]

    The isomorphism between the Hall algebra and the symmetric
    functions described above is implemented as a coercion::

        sage: R = PolynomialRing(ZZ, 'q').fraction_field()
        sage: q = R.gen()
        sage: H = HallAlgebra(R, q)
        sage: e = SymmetricFunctions(R).e()
        sage: e(H[1,1,1])
        1/q^3*e[3]

    We can also do computations with any special value of ``q``,
    such as `0` or `1` or (most commonly) a prime power. Here
    is an example using a prime::

        sage: H = HallAlgebra(ZZ, 2)
        sage: H[2,1]*H[1,1]
        H[3, 2] + 3*H[3, 1, 1] + 6*H[2, 2, 1] + 28*H[2, 1, 1, 1]
        sage: H[3,1]*H[2]
        H[5, 1] + H[4, 2] + 6*H[3, 3] + 3*H[4, 1, 1] + 8*H[3, 2, 1]
        sage: H[2,1,1]*H[3,1]
        H[5, 2, 1] + 2*H[4, 3, 1] + 6*H[4, 2, 2] + 7*H[5, 1, 1, 1]
         + 19*H[4, 2, 1, 1] + 24*H[3, 3, 1, 1] + 48*H[3, 2, 2, 1]
         + 105*H[4, 1, 1, 1, 1] + 224*H[3, 2, 1, 1, 1]
        sage: I = H.monomial_basis()
        sage: H(I[2,1,1])
        H[3, 1] + 3*H[2, 2] + 13*H[2, 1, 1] + 105*H[1, 1, 1, 1]
        sage: I(H[2,1,1])
        I[3, 1] - 15*I[4]

    If `q` is set to `1`, the coercion to the symmetric functions
    sends `I_{\lambda}` to `m_{\lambda}`::

        sage: H = HallAlgebra(QQ, 1)
        sage: H[2,1] * H[2,1]
        H[4, 2] + 2*H[3, 3] + 2*H[4, 1, 1] + 2*H[3, 2, 1] + 6*H[2, 2, 2] + 4*H[2, 2, 1, 1]
        sage: m = SymmetricFunctions(QQ).m()
        sage: m[2,1] * m[2,1]
        4*m[2, 2, 1, 1] + 6*m[2, 2, 2] + 2*m[3, 2, 1] + 2*m[3, 3] + 2*m[4, 1, 1] + m[4, 2]
        sage: m(H[3,1])
        m[3, 1]

    We can set `q` to `0` (but should keep in mind that we don't get
    the Schur functions this way)::

        sage: H = HallAlgebra(QQ, 0)
        sage: H[2,1] * H[2,1]
        H[4, 2] + H[3, 3] + H[4, 1, 1] - H[3, 2, 1] - H[3, 1, 1, 1]

    TESTS:

    The coefficients are actually Laurent polynomials in general, so we don't
    have to work over the fraction field of `\ZZ[q]`. This didn't work before
    :issue:`15345`::

        sage: R.<q> = LaurentPolynomialRing(ZZ)
        sage: H = HallAlgebra(R, q)
        sage: I = H.monomial_basis()
        sage: hi = H(I[2,1]); hi
        H[2, 1] + (1+q+q^2)*H[1, 1, 1]
        sage: hi.parent() is H
        True
        sage: h22 = H[2]*H[2]; h22
        H[4] - (1-q)*H[3, 1] + (q+q^2)*H[2, 2]
        sage: h22.parent() is H
        True
        sage: e = SymmetricFunctions(R).e()
        sage: e(H[1,1,1])
        (q^-3)*e[3]
    """
    def __init__(self, base_ring, q, prefix='H'):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: H = HallAlgebra(R, q)
            sage: TestSuite(H).run()
            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: H = HallAlgebra(R, q)
            sage: TestSuite(H).run() # long time
            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: H = HallAlgebra(R, q)
            sage: TestSuite(H).run() # long time
        """
        self._q = q
        try:
            q_inverse = q**-1
            if q_inverse not in base_ring:
                hopf_structure = False
            else:
                hopf_structure = True
        except Exception:
            hopf_structure = False
        if hopf_structure:
            category = HopfAlgebrasWithBasis(base_ring)
        else:
            category = AlgebrasWithBasis(base_ring)
        CombinatorialFreeModule.__init__(self, base_ring, Partitions(),
                                         prefix=prefix, bracket=False,
                                         sorting_key=cmp_to_key(transpose_cmp),
                                         category=category)

        # Coercions
        I = self.monomial_basis()
        M = I.module_morphism(I._to_natural_on_basis, codomain=self,
                              triangular='upper', unitriangular=True,
                              inverse_on_support=lambda x: x.conjugate(),
                              invertible=True)
        M.register_as_coercion()
        (~M).register_as_coercion()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: HallAlgebra(R, q)
            Hall algebra with q=q over Univariate Polynomial Ring in q over Integer Ring
        """
        return "Hall algebra with q={} over {}".format(self._q, self.base_ring())

    def one_basis(self):
        """
        Return the index of the basis element `1`.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: H = HallAlgebra(R, q)
            sage: H.one_basis()
            []
        """
        return Partition([])

    def product_on_basis(self, mu, la):
        """
        Return the product of the two basis elements indexed by ``mu``
        and ``la``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: H = HallAlgebra(R, q)
            sage: H.product_on_basis(Partition([1,1]), Partition([1]))
            H[2, 1] + (q^2+q+1)*H[1, 1, 1]
            sage: H.product_on_basis(Partition([2,1]), Partition([1,1]))
            H[3, 2] + (q+1)*H[3, 1, 1] + (q^2+q)*H[2, 2, 1] + (q^4+q^3+q^2)*H[2, 1, 1, 1]
            sage: H.product_on_basis(Partition([3,2]), Partition([2,1]))
            H[5, 3] + (q+1)*H[4, 4] + q*H[5, 2, 1] + (2*q^2-1)*H[4, 3, 1]
             + (q^3+q^2)*H[4, 2, 2] + (q^4+q^3)*H[3, 3, 2]
             + (q^4-q^2)*H[4, 2, 1, 1] + (q^5+q^4-q^3-q^2)*H[3, 3, 1, 1]
             + (q^6+q^5)*H[3, 2, 2, 1]
            sage: H.product_on_basis(Partition([3,1,1]), Partition([2,1]))
            H[5, 2, 1] + q*H[4, 3, 1] + (q^2-1)*H[4, 2, 2]
             + (q^3+q^2)*H[3, 3, 2] + (q^2+q+1)*H[5, 1, 1, 1]
             + (2*q^3+q^2-q-1)*H[4, 2, 1, 1] + (q^4+2*q^3+q^2)*H[3, 3, 1, 1]
             + (q^5+q^4)*H[3, 2, 2, 1] + (q^6+q^5+q^4-q^2-q-1)*H[4, 1, 1, 1, 1]
             + (q^7+q^6+q^5)*H[3, 2, 1, 1, 1]
        """
        # Check conditions for multiplying by 1
        if len(mu) == 0:
            return self.monomial(la)
        if len(la) == 0:
            return self.monomial(mu)

        if all(x == 1 for x in la):
            return self.sum_of_terms([(p, hall_polynomial(p, mu, la, self._q))
                                      for p in Partitions(sum(mu) + len(la))],
                                     distinct=True)

        I = HallAlgebraMonomials(self.base_ring(), self._q)
        mu = self.monomial(mu)
        la = self.monomial(la)
        return self(I(mu) * I(la))

    def coproduct_on_basis(self, la):
        """
        Return the coproduct of the basis element indexed by ``la``.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: H = HallAlgebra(R, q)
            sage: H.coproduct_on_basis(Partition([1,1]))
            H[] # H[1, 1] + 1/q*H[1] # H[1] + H[1, 1] # H[]
            sage: H.coproduct_on_basis(Partition([2]))
            H[] # H[2] + ((q-1)/q)*H[1] # H[1] + H[2] # H[]
            sage: H.coproduct_on_basis(Partition([2,1]))
            H[] # H[2, 1] + ((q^2-1)/q^2)*H[1] # H[1, 1] + 1/q*H[1] # H[2]
             + ((q^2-1)/q^2)*H[1, 1] # H[1] + 1/q*H[2] # H[1] + H[2, 1] # H[]

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: H = HallAlgebra(R, q)
            sage: H.coproduct_on_basis(Partition([2]))
            H[] # H[2] - (q^-1-1)*H[1] # H[1] + H[2] # H[]
            sage: H.coproduct_on_basis(Partition([2,1]))
            H[] # H[2, 1] - (q^-2-1)*H[1] # H[1, 1] + (q^-1)*H[1] # H[2]
             - (q^-2-1)*H[1, 1] # H[1] + (q^-1)*H[2] # H[1] + H[2, 1] # H[]
        """
        S = self.tensor_square()
        if all(x == 1 for x in la):
            n = len(la)
            return S.sum_of_terms([( (Partition([1]*r), Partition([1]*(n-r))), self._q**(-r*(n-r)) )
                                   for r in range(n+1)], distinct=True)

        I = HallAlgebraMonomials(self.base_ring(), self._q)
        la = self.monomial(la)
        return S(I(la).coproduct())

    def antipode_on_basis(self, la):
        """
        Return the antipode of the basis element indexed by ``la``.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: H = HallAlgebra(R, q)
            sage: H.antipode_on_basis(Partition([1,1]))
            1/q*H[2] + 1/q*H[1, 1]
            sage: H.antipode_on_basis(Partition([2]))
            -1/q*H[2] + ((q^2-1)/q)*H[1, 1]

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: H = HallAlgebra(R, q)
            sage: H.antipode_on_basis(Partition([1,1]))
            (q^-1)*H[2] + (q^-1)*H[1, 1]
            sage: H.antipode_on_basis(Partition([2]))
            -(q^-1)*H[2] - (q^-1-q)*H[1, 1]
        """
        if all(x == 1 for x in la):
            r = len(la)
            q = (-1) ** r * self._q ** (-(r * (r - 1)) // 2)
            return self._from_dict({p: q for p in Partitions(r)})

        I = HallAlgebraMonomials(self.base_ring(), self._q)
        return self(I(self.monomial(la)).antipode())

    def counit(self, x):
        """
        Return the counit of the element ``x``.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: H = HallAlgebra(R, q)
            sage: H.counit(H.an_element())
            2
        """
        return x.coefficient(self.one_basis())

    def monomial_basis(self):
        """
        Return the basis of the Hall algebra given by monomials in the
        `I_{(1^r)}`.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: H = HallAlgebra(R, q)
            sage: H.monomial_basis()
            Hall algebra with q=q over Univariate Polynomial Ring in q over
             Integer Ring in the monomial basis
        """
        return HallAlgebraMonomials(self.base_ring(), self._q)

    def __getitem__(self, la):
        """
        Return the basis element indexed by ``la``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: H = HallAlgebra(R, q)
            sage: H[[]]
            H[]
            sage: H[2]
            H[2]
            sage: H[[2]]
            H[2]
            sage: H[2,1]
            H[2, 1]
            sage: H[Partition([2,1])]
            H[2, 1]
            sage: H[(2,1)]
            H[2, 1]
        """
        if la in ZZ:
            return self.monomial(Partition([la]))
        return self.monomial(Partition(la))

    class Element(CombinatorialFreeModule.Element):
        def scalar(self, y):
            r"""
            Return the scalar product of ``self`` and ``y``.

            The scalar product is given by

            .. MATH::

                (I_{\lambda}, I_{\mu}) = \delta_{\lambda,\mu}
                \frac{1}{a_{\lambda}},

            where `a_{\lambda}` is given by

            .. MATH::

                a_{\lambda} = q^{|\lambda| + 2 n(\lambda)} \prod_k
                \prod_{i=1}^{l_k} (1 - q^{-i})

            where `n(\lambda) = \sum_i (i - 1) \lambda_i` and
            `\lambda = (1^{l_1}, 2^{l_2}, \ldots, m^{l_m})`.

            Note that `a_{\lambda}` can be interpreted as the number
            of automorphisms of a certain object in a category
            corresponding to `\lambda`. See Lemma 2.8 in [Sch2006]_
            for details.

            EXAMPLES::

                sage: R.<q> = ZZ[]
                sage: H = HallAlgebra(R, q)
                sage: H[1].scalar(H[1])
                1/(q - 1)
                sage: H[2].scalar(H[2])
                1/(q^2 - q)
                sage: H[2,1].scalar(H[2,1])
                1/(q^5 - 2*q^4 + q^3)
                sage: H[1,1,1,1].scalar(H[1,1,1,1])
                1/(q^16 - q^15 - q^14 + 2*q^11 - q^8 - q^7 + q^6)
                sage: H.an_element().scalar(H.an_element())
                (4*q^2 + 9)/(q^2 - q)
            """
            q = self.parent()._q
            f = lambda la: ~( q**(sum(la) + 2*la.weighted_size())
                              * prod(prod((1 - q**-i) for i in range(1,k+1))
                                     for k in la.to_exp()) )
            y = self.parent()(y)
            ret = q.parent().zero()
            for mx, cx in self:
                cy = y.coefficient(mx)
                if cy != 0:
                    ret += cx * cy * f(mx)
            return ret


class HallAlgebraMonomials(CombinatorialFreeModule):
    r"""
    The classical Hall algebra given in terms of monomials in the
    `I_{(1^r)}`.

    We first associate a monomial `I_{(1^{r_1})} I_{(1^{r_2})} \cdots
    I_{(1^{r_k})}` with the composition `(r_1, r_2, \ldots, r_k)`. However
    since `I_{(1^r)}` commutes with `I_{(1^s)}`, the basis is indexed
    by partitions.

    EXAMPLES:

    We use the fraction field of `\ZZ[q]` for our initial example::

        sage: R = PolynomialRing(ZZ, 'q').fraction_field()
        sage: q = R.gen()
        sage: H = HallAlgebra(R, q)
        sage: I = H.monomial_basis()

    We check that the basis conversions are mutually inverse::

        sage: all(H(I(H[p])) == H[p] for i in range(7) for p in Partitions(i))
        True
        sage: all(I(H(I[p])) == I[p] for i in range(7) for p in Partitions(i))
        True

    Since Laurent polynomials are sufficient, we run the same check with
    the Laurent polynomial ring `\ZZ[q, q^{-1}]`::

        sage: R.<q> = LaurentPolynomialRing(ZZ)
        sage: H = HallAlgebra(R, q)
        sage: I = H.monomial_basis()
        sage: all(H(I(H[p])) == H[p] for i in range(6) for p in Partitions(i)) # long time
        True
        sage: all(I(H(I[p])) == I[p] for i in range(6) for p in Partitions(i)) # long time
        True

    We can also convert to the symmetric functions. The natural basis
    corresponds to the Hall-Littlewood basis (up to a renormalization and
    an inversion of the `q` parameter), and this basis corresponds
    to the elementary basis (up to a renormalization)::

        sage: Sym = SymmetricFunctions(R)
        sage: e = Sym.e()
        sage: e(I[2,1])
        (q^-1)*e[2, 1]
        sage: e(I[4,2,2,1])
        (q^-8)*e[4, 2, 2, 1]
        sage: HLP = Sym.hall_littlewood(q).P()
        sage: H(I[2,1])
        H[2, 1] + (1+q+q^2)*H[1, 1, 1]
        sage: HLP(e[2,1])
        (1+q+q^2)*HLP[1, 1, 1] + HLP[2, 1]
        sage: all( e(H[lam]) == q**-sum([i * x for i, x in enumerate(lam)])
        ....:          * e(HLP[lam]).map_coefficients(lambda p: p(q**(-1)))
        ....:      for lam in Partitions(4) )
        True

    We can also do computations using a prime power::

        sage: H = HallAlgebra(ZZ, 3)
        sage: I = H.monomial_basis()
        sage: i_elt = I[2,1]*I[1,1]; i_elt
        I[2, 1, 1, 1]
        sage: H(i_elt)
        H[4, 1] + 7*H[3, 2] + 37*H[3, 1, 1] + 136*H[2, 2, 1]
         + 1495*H[2, 1, 1, 1] + 62920*H[1, 1, 1, 1, 1]
    """
    def __init__(self, base_ring, q, prefix='I'):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: TestSuite(I).run()
            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: TestSuite(I).run()
            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: TestSuite(I).run()
        """
        self._q = q
        try:
            q_inverse = q**-1
            if q_inverse not in base_ring:
                hopf_structure = False
            else:
                hopf_structure = True
        except Exception:
            hopf_structure = False
        if hopf_structure:
            category = HopfAlgebrasWithBasis(base_ring)
        else:
            category = AlgebrasWithBasis(base_ring)
        CombinatorialFreeModule.__init__(self, base_ring, Partitions(),
                                         prefix=prefix, bracket=False,
                                         category=category)

        # Coercions
        if hopf_structure:
            e = SymmetricFunctions(base_ring).e()
            f = lambda la: q ** sum(-((r * (r - 1)) // 2) for r in la)
            M = self.module_morphism(diagonal=f, codomain=e)
            M.register_as_coercion()
            (~M).register_as_coercion()

    @cached_method
    def _to_natural_on_basis(self, a):
        """
        Return the basis element indexed by ``a`` converted into
        the partition basis.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I._to_natural_on_basis(Partition([3]))
            H[1, 1, 1]
            sage: I._to_natural_on_basis(Partition([2,1,1]))
            H[3, 1] + (q+1)*H[2, 2] + (2*q^2+2*q+1)*H[2, 1, 1]
             + (q^5+2*q^4+3*q^3+3*q^2+2*q+1)*H[1, 1, 1, 1]
        """
        H = HallAlgebra(self.base_ring(), self._q)
        return reduce(lambda cur,r: cur * H.monomial(Partition([1]*r)), a, H.one())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: HallAlgebra(R, q).monomial_basis()
            Hall algebra with q=q over Univariate Polynomial Ring in q over
             Integer Ring in the monomial basis
        """
        return "Hall algebra with q={} over {} in the monomial basis".format(self._q, self.base_ring())

    def one_basis(self):
        """
        Return the index of the basis element `1`.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I.one_basis()
            []
        """
        return Partition([])

    def product_on_basis(self, a, b):
        """
        Return the product of the two basis elements indexed by ``a``
        and ``b``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I.product_on_basis(Partition([4,2,1]), Partition([3,2,1]))
            I[4, 3, 2, 2, 1, 1]
        """
        return self.monomial(Partition(sorted(list(a) + list(b), reverse=True)))

    def coproduct_on_basis(self, a):
        """
        Return the coproduct of the basis element indexed by ``a``.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I.coproduct_on_basis(Partition([1]))
            I[] # I[1] + I[1] # I[]
            sage: I.coproduct_on_basis(Partition([2]))
            I[] # I[2] + 1/q*I[1] # I[1] + I[2] # I[]
            sage: I.coproduct_on_basis(Partition([2,1]))
            I[] # I[2, 1] + 1/q*I[1] # I[1, 1] + I[1] # I[2]
             + 1/q*I[1, 1] # I[1] + I[2] # I[1] + I[2, 1] # I[]

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I.coproduct_on_basis(Partition([2,1]))
            I[] # I[2, 1] + (q^-1)*I[1] # I[1, 1] + I[1] # I[2]
             + (q^-1)*I[1, 1] # I[1] + I[2] # I[1] + I[2, 1] # I[]
        """
        S = self.tensor_square()
        return S.prod(S.sum_of_terms([( (Partition([r]), Partition([n-r]) ), self._q**(-r*(n-r)) )
                                      for r in range(n+1)], distinct=True) for n in a)

    def antipode_on_basis(self, a):
        """
        Return the antipode of the basis element indexed by ``a``.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I.antipode_on_basis(Partition([1]))
            -I[1]
            sage: I.antipode_on_basis(Partition([2]))
            1/q*I[1, 1] - I[2]
            sage: I.antipode_on_basis(Partition([2,1]))
            -1/q*I[1, 1, 1] + I[2, 1]

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I.antipode_on_basis(Partition([2,1]))
            -(q^-1)*I[1, 1, 1] + I[2, 1]
        """
        H = HallAlgebra(self.base_ring(), self._q)
        cur = self.one()
        for r in a:
            q = (-1) ** r * self._q ** (-(r * (r - 1)) // 2)
            cur *= self(H._from_dict({p: q for p in Partitions(r)}))
        return cur

    def counit(self, x):
        """
        Return the counit of the element ``x``.

        EXAMPLES::

            sage: R = PolynomialRing(ZZ, 'q').fraction_field()
            sage: q = R.gen()
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I.counit(I.an_element())
            2
        """
        return x.coefficient(self.one_basis())

    def __getitem__(self, a):
        """
        Return the basis element indexed by ``a``.

        EXAMPLES::

            sage: R.<q> = ZZ[]
            sage: I = HallAlgebra(R, q).monomial_basis()
            sage: I[3,1,1] + 3*I[1,1]
            3*I[1, 1] + I[3, 1, 1]
            sage: I[Partition([3,2,2])]
            I[3, 2, 2]
            sage: I[2]
            I[2]
            sage: I[[2]]
            I[2]
            sage: I[[]]
            I[]
        """
        if a in ZZ:
            return self.monomial(Partition([a]))
        return self.monomial(Partition(a))

    class Element(CombinatorialFreeModule.Element):
        def scalar(self, y):
            r"""
            Return the scalar product of ``self`` and ``y``.

            The scalar product is computed by converting into the
            natural basis.

            EXAMPLES::

                sage: R.<q> = ZZ[]
                sage: I = HallAlgebra(R, q).monomial_basis()
                sage: I[1].scalar(I[1])
                1/(q - 1)
                sage: I[2].scalar(I[2])
                1/(q^4 - q^3 - q^2 + q)
                sage: I[2,1].scalar(I[2,1])
                (2*q + 1)/(q^6 - 2*q^5 + 2*q^3 - q^2)
                sage: I[1,1,1,1].scalar(I[1,1,1,1])
                24/(q^4 - 4*q^3 + 6*q^2 - 4*q + 1)
                sage: I.an_element().scalar(I.an_element())
                (4*q^4 - 4*q^2 + 9)/(q^4 - q^3 - q^2 + q)
            """
            H = HallAlgebra(self.parent().base_ring(), self.parent()._q)
            return H(self).scalar(H(y))
