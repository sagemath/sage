# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Hecke Character Basis

The basis of symmetric functions given by characters of the
Hecke algebra (of type `A`).

AUTHORS:

- Travis Scrimshaw (2017-08): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
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
# ****************************************************************************

from sage.combinat.partition import _Partitions, Partitions
from sage.combinat.sf.multiplicative import SymmetricFunctionAlgebra_multiplicative


class HeckeCharacter(SymmetricFunctionAlgebra_multiplicative):
    r"""
    Basis of the symmetric functions that gives the characters of the
    Hecke algebra in analogy to the Frobenius formula for the
    symmetric group.

    Consider the Hecke algebra `H_n(q)` with quadratic relations

    .. MATH::

        T_i^2 = (q - 1) T_i + q.

    Let `\mu` be a partition of `n` with length `\ell`. The character
    `\chi` of a `H_n(q)`-representation is completely determined by
    the elements `T_{\gamma_{\mu}}`, where

    .. MATH::

        \gamma_{\mu} = (\mu_1, \ldots, 1) (\mu_2 + \mu_1, \ldots, 1 + \mu_1)
        \cdots (n, \ldots, 1 + \sum_{i < \ell} \mu_i),

    (written in cycle notation). We define a basis of the symmetric
    functions by

    .. MATH::

        \bar{q}_{\mu} = \sum_{\lambda \vdash n}
        \chi^{\lambda}(T_{\gamma_{\mu}}) s_{\lambda}.

    INPUT:

    - ``sym`` -- the ring of symmetric functions
    - ``q`` -- (default: ``'q'``) the parameter `q`

    EXAMPLES::

        sage: q = ZZ['q'].fraction_field().gen()
        sage: Sym = SymmetricFunctions(q.parent())
        sage: qbar = Sym.hecke_character(q)
        sage: qbar[2] * qbar[3] * qbar[3,1]
        qbar[3, 3, 2, 1]

        sage: s = Sym.s()
        sage: s(qbar([2]))
        -s[1, 1] + q*s[2]
        sage: s(qbar([4]))
        -s[1, 1, 1, 1] + q*s[2, 1, 1] - q^2*s[3, 1] + q^3*s[4]
        sage: qbar(s[2])
        (1/(q+1))*qbar[1, 1] + (1/(q+1))*qbar[2]
        sage: qbar(s[1,1])
        (q/(q+1))*qbar[1, 1] - (1/(q+1))*qbar[2]

        sage: s(qbar[2,1])
        -s[1, 1, 1] + (q-1)*s[2, 1] + q*s[3]
        sage: qbar(s[2,1])
        (q/(q^2+q+1))*qbar[1, 1, 1] + ((q-1)/(q^2+q+1))*qbar[2, 1]
         - (1/(q^2+q+1))*qbar[3]

    We compute character tables for Hecke algebras, which correspond
    to the transition matrix from the `\bar{q}` basis to the Schur
    basis::

        sage: qbar.transition_matrix(s, 1)
        [1]
        sage: qbar.transition_matrix(s, 2)
        [ q -1]
        [ 1  1]
        sage: qbar.transition_matrix(s, 3)
        [  q^2    -q     1]
        [    q q - 1    -1]
        [    1     2     1]
        sage: qbar.transition_matrix(s, 4)
        [      q^3      -q^2         0         q        -1]
        [      q^2   q^2 - q        -q    -q + 1         1]
        [      q^2 q^2 - 2*q   q^2 + 1  -2*q + 1         1]
        [        q   2*q - 1     q - 1     q - 2        -1]
        [        1         3         2         3         1]

    We can do computations with a specialized `q` to a generic element
    of the base ring. We compute some examples with `q = 2`::

        sage: qbar = Sym.qbar(q=2)
        sage: s = Sym.schur()
        sage: qbar(s[2,1])
        2/7*qbar[1, 1, 1] + 1/7*qbar[2, 1] - 1/7*qbar[3]
        sage: s(qbar[2,1])
        -s[1, 1, 1] + s[2, 1] + 2*s[3]

    REFERENCES:

    - [Ram1991]_
    - [RR1997]_
    """
    @staticmethod
    def __classcall__(cls, Sym, q='q'):
        """
        Normalize the arguments.

        TESTS::

            sage: R.<q, t> = QQ[]
            sage: B1 = SymmetricFunctions(R).qbar()
            sage: B2 = SymmetricFunctions(R).qbar(q)
            sage: B3 = SymmetricFunctions(R).qbar(t)
            sage: B1 is B2
            True
            sage: B1 == B3
            False
        """
        return super().__classcall__(cls, Sym, Sym.base_ring()(q))

    def __init__(self, sym, q):
        r"""
        Initialize ``self``.

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(ZZ['q']))
            sage: qbar = Sym.qbar()
            sage: TestSuite(qbar).run()

            sage: Sym = SymmetricFunctions(QQ)
            sage: qbar = Sym.qbar(q=2)
            sage: TestSuite(qbar).run()

        Check that the conversion `q \to p \to s` agrees with
        the definition of `q \to s` from [Ram1991]_::

            sage: Sym = SymmetricFunctions(FractionField(ZZ['q']))
            sage: qbar = Sym.qbar()
            sage: s = Sym.s()
            sage: q = qbar.q()
            sage: def to_schur(mu):
            ....:     if not mu:
            ....:        return s.one()
            ....:     mone = -qbar.base_ring().one()
            ....:     return s.prod(sum(mone**(r-m) * q**(m-1)
            ....:                       * s[Partition([m] + [1]*(r-m))]
            ....:                       for m in range(1, r+1))
            ....:                   for r in mu)
            sage: phi = qbar.module_morphism(to_schur, codomain=s)
            sage: all(phi(qbar[mu]) == s(qbar[mu]) for n in range(6)
            ....:     for mu in Partitions(n))
            True
        """
        self.q = q
        SymmetricFunctionAlgebra_multiplicative.__init__(self, sym,
            basis_name="Hecke character with q={}".format(self.q),
            prefix='qbar')
        self._p = sym.power()

        # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
        # category = ModulesWithBasis(self._sym.base_ring())
        self.register_coercion(self._p._module_morphism(self._p_to_qbar_on_basis,
                                                           codomain=self))
        self._p.register_coercion(self._module_morphism(self._qbar_to_p_on_basis,
                                                        codomain=self._p))

    def construction(self):
        """
        Return a pair ``(F, R)``, where ``F`` is a
        :class:`SymmetricFunctionsFunctor` and `R` is a ring, such
        that ``F(R)`` returns ``self``.

        EXAMPLES::

            sage: qbar = SymmetricFunctions(QQ['q']).qbar('q')
            sage: qbar.construction()
            (SymmetricFunctionsFunctor[Hecke character with q=q],
             Univariate Polynomial Ring in q over Rational Field)
        """

        from sage.combinat.sf.sfa import SymmetricFunctionsFunctor
        return (SymmetricFunctionsFunctor(self, self.basis_name(), self.q),
                self.base_ring())

    def _p_to_qbar_on_generator(self, n):
        r"""
        Convert `p_n` to ``self``.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: qbar = SymmetricFunctions(QQ['q'].fraction_field()).qbar('q')
            sage: qbar._p_to_qbar_on_generator(3)
            ((q^2-2*q+1)/(q^2+q+1))*qbar[1, 1, 1]
             + ((-3*q+3)/(q^2+q+1))*qbar[2, 1]
             + (3/(q^2+q+1))*qbar[3]
            sage: qbar = SymmetricFunctions(QQ).qbar(-1)
            sage: qbar._p_to_qbar_on_generator(3)
            2*qbar[2, 1] + 3*qbar[3]
        """
        if n == 1:
            return self([1])
        q = self.q
        if q**n == self.base_ring().one():
            raise ValueError("the parameter q=%s must not be a %s root of unity" % (q,n))
        out = n * self([n]) - sum((q**i-1) * self._p_to_qbar_on_generator(i)
                                  * self([n-i]) for i in range(1,n) if q**i != 1)
        return out*(q-1) / (q**n-1)

    def _p_to_qbar_on_basis(self, mu):
        r"""
        Convert the power sum basis element indexed by ``mu`` to ``self``.

        INPUT:

        - ``mu`` -- a partition or a list of nonnegative integers

        EXAMPLES::

            sage: qbar = SymmetricFunctions(QQ['q'].fraction_field()).qbar('q')
            sage: qbar._p_to_qbar_on_basis([3,1])
            ((q^2-2*q+1)/(q^2+q+1))*qbar[1, 1, 1, 1]
             + ((-3*q+3)/(q^2+q+1))*qbar[2, 1, 1]
             + (3/(q^2+q+1))*qbar[3, 1]
            sage: qbar = SymmetricFunctions(QQ).qbar(2)
            sage: qbar._p_to_qbar_on_basis([3,1])
            1/7*qbar[1, 1, 1, 1] - 3/7*qbar[2, 1, 1] + 3/7*qbar[3, 1]
        """
        return self.prod(self._p_to_qbar_on_generator(p) for p in mu)

    def _qbar_to_p_on_generator(self, n):
        r"""
        Convert a generator of the basis indexed by ``n`` to the
        power sum basis.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: qbar = SymmetricFunctions(QQ['q'].fraction_field()).qbar('q')
            sage: qbar._qbar_to_p_on_generator(3)
            (1/6*q^2-1/3*q+1/6)*p[1, 1, 1]
             + (1/2*q^2-1/2)*p[2, 1]
             + (1/3*q^2+1/3*q+1/3)*p[3]
            sage: qbar = SymmetricFunctions(QQ).qbar(-1)
            sage: qbar._qbar_to_p_on_generator(3)
            2/3*p[1, 1, 1] + 1/3*p[3]
        """
        if n == 1:
            return self._p([1])
        q = self.q
        BR = self.base_ring()
        return q**(n-1) * self._p.sum(sum(q**(-i) for i in range(mu[0]))
                                      * BR.prod(1 - q**(-p) for p in mu[1:])
                                      * self._p(mu) / mu.centralizer_size()
                                      for mu in Partitions(n)
                                      if not any(q**p == 1 for p in mu[1:]))

    def _qbar_to_p_on_basis(self, mu):
        r"""
        Convert a basis element indexed by the partition ``mu``
        to the power basis.

        INPUT:

        - ``mu`` -- a partition or a list of nonnegative integers

        EXAMPLES::

            sage: qbar = SymmetricFunctions(QQ['q'].fraction_field()).qbar('q')
            sage: qbar._qbar_to_p_on_basis([3,1])
            (1/6*q^2-1/3*q+1/6)*p[1, 1, 1, 1]
             + (1/2*q^2-1/2)*p[2, 1, 1]
             + (1/3*q^2+1/3*q+1/3)*p[3, 1]
            sage: qbar = SymmetricFunctions(QQ).qbar(-1)
            sage: qbar._qbar_to_p_on_basis([3,1])
            2/3*p[1, 1, 1, 1] + 1/3*p[3, 1]
        """
        return self._p.prod(self._qbar_to_p_on_generator(p) for p in mu)

    def coproduct_on_generators(self, r):
        r"""
        Return the coproduct on the generator `\bar{q}_r` of ``self``.

        Define the coproduct on `\bar{q}_r` by

        .. MATH::

            \Delta(\bar{q}_r) = \bar{q}_0 \otimes \bar{q}_r
            + (q - 1) \sum_{j=1}^{r-1} \bar{q}_j \otimes \bar{q}_{r-j}
            + \bar{q}_r \otimes \bar{q}_0.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: Sym = SymmetricFunctions(q.parent())
            sage: qbar = Sym.hecke_character()
            sage: s = Sym.s()
            sage: qbar[2].coproduct()
            qbar[] # qbar[2] + (q-1)*qbar[1] # qbar[1] + qbar[2] # qbar[]
        """
        def P(i):
            return _Partitions([i]) if i else _Partitions([])
        T = self.tensor_square()
        one = self.base_ring().one()
        q = self.q
        return T.sum_of_terms(((P(j), P(r-j)), one if j in [0,r] else q-one)
                              for j in range(r+1))
