# sage.doctest: needs sage.combinat sage.modules
"""
Elementary symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#                     2012 Anne Schilling <anne@math.ucdavis.edu>
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
#*****************************************************************************
from sage.arith.misc import binomial, factorial
from sage.combinat.partition import Partition
from sage.combinat.sf import multiplicative, classical
from sage.misc.misc_c import prod
from sage.rings.infinity import infinity


###################################
#                                 #
# Elementary Symmetric Functions  #
#                                 #
###################################


class SymmetricFunctionAlgebra_elementary(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, Sym):
        """
        A class for methods for the elementary basis of the symmetric functions.

        INPUT:

        - ``self`` -- an elementary basis of the symmetric functions
        - ``Sym`` -- an instance of the ring of symmetric functions

        TESTS::

            sage: e = SymmetricFunctions(QQ).e()
            sage: e == loads(dumps(e))
            True
            sage: TestSuite(e).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(e).run(elements = [e[1,1]+e[2], e[1]+2*e[1,1]])
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, Sym, "elementary", 'e')

    def _dual_basis_default(self):
        """
        Return the default value for ``self.dual_basis()``.

        This method returns the dual basis to the elementary basis
        with respect to the standard scalar product, that is the
        forgotten basis.

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: e.dual_basis()
            Symmetric Functions over Rational Field in the forgotten basis

        TESTS::

            sage: e._dual_basis_default() is e.dual_basis()
            True
        """
        return self.dual_basis(scalar=None, prefix='f', basis_name='forgotten')

    def coproduct_on_generators(self, i):
        r"""
        Return the coproduct on ``self[i]``.

        INPUT:

        - ``self`` -- an elementary basis of the symmetric functions
        - ``i`` -- nonnegative integer

        OUTPUT: the coproduct on the elementary generator `e(i)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: e = Sym.elementary()
            sage: e.coproduct_on_generators(2)
            e[] # e[2] + e[1] # e[1] + e[2] # e[]
            sage: e.coproduct_on_generators(0)
            e[] # e[]
        """
        def P(i):
            return Partition([i]) if i else Partition([])
        T = self.tensor_square()
        return T.sum_of_monomials( (P(j), P(i-j)) for j in range(i+1) )

    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        def omega(self):
            r"""
            Return the image of ``self`` under the omega automorphism.

            The *omega automorphism* is defined to be the unique algebra
            endomorphism `\omega` of the ring of symmetric functions that
            satisfies `\omega(e_k) = h_k` for all positive integers `k`
            (where `e_k` stands for the `k`-th elementary symmetric
            function, and `h_k` stands for the `k`-th complete homogeneous
            symmetric function). It furthermore is a Hopf algebra
            endomorphism and an involution, and it is also known as the
            *omega involution*. It sends the power-sum symmetric function
            `p_k` to `(-1)^{k-1} p_k` for every positive integer `k`.

            The images of some bases under the omega automorphism are given by

            .. MATH::

                \omega(e_{\lambda}) = h_{\lambda}, \qquad
                \omega(h_{\lambda}) = e_{\lambda}, \qquad
                \omega(p_{\lambda}) = (-1)^{|\lambda| - \ell(\lambda)}
                p_{\lambda}, \qquad
                \omega(s_{\lambda}) = s_{\lambda^{\prime}},

            where `\lambda` is any partition, where `\ell(\lambda)` denotes
            the length (:meth:`~sage.combinat.partition.Partition.length`)
            of the partition `\lambda`, where `\lambda^{\prime}` denotes the
            conjugate partition
            (:meth:`~sage.combinat.partition.Partition.conjugate`) of
            `\lambda`, and where the usual notations for bases are used
            (`e` = elementary, `h` = complete homogeneous, `p` = powersum,
            `s` = Schur).

            :meth:`omega_involution` is a synonym for the :meth:`omega`
            method.

            EXAMPLES::

                sage: e = SymmetricFunctions(QQ).e()
                sage: a = e([2,1]); a
                e[2, 1]
                sage: a.omega()
                e[1, 1, 1] - e[2, 1]

            ::

                sage: h = SymmetricFunctions(QQ).h()
                sage: h(e([2,1]).omega())
                h[2, 1]
            """
            e = self.parent()
            h = e.realization_of().h()
            return e( h._from_element(self) )

        omega_involution = omega

        def verschiebung(self, n):
            r"""
            Return the image of the symmetric function ``self`` under the
            `n`-th Verschiebung operator.

            The `n`-th Verschiebung operator `\mathbf{V}_n` is defined to be
            the unique algebra endomorphism `V` of the ring of symmetric
            functions that satisfies `V(h_r) = h_{r/n}` for every positive
            integer `r` divisible by `n`, and satisfies `V(h_r) = 0` for
            every positive integer `r` not divisible by `n`. This operator
            `\mathbf{V}_n` is a Hopf algebra endomorphism. For every
            nonnegative integer `r` with `n \mid r`, it satisfies

            .. MATH::

                \mathbf{V}_n(h_r) = h_{r/n},
                \quad \mathbf{V}_n(p_r) = n p_{r/n},
                \quad \mathbf{V}_n(e_r) = (-1)^{r - r/n} e_{r/n}

            (where `h` is the complete homogeneous basis, `p` is the
            powersum basis, and `e` is the elementary basis). For every
            nonnegative integer `r` with `n \nmid r`, it satisfes

            .. MATH::

                \mathbf{V}_n(h_r) = \mathbf{V}_n(p_r) = \mathbf{V}_n(e_r) = 0.

            The `n`-th Verschiebung operator is also called the `n`-th
            Verschiebung endomorphism. Its name derives from the Verschiebung
            (German for "shift") endomorphism of the Witt vectors.

            The `n`-th Verschiebung operator is adjoint to the `n`-th
            Frobenius operator (see :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.frobenius`
            for its definition) with respect to the Hall scalar product
            (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.scalar`).

            The action of the `n`-th Verschiebung operator on the Schur basis
            can also be computed explicitly. The following (probably clumsier
            than necessary) description can be obtained by solving exercise
            7.61 in Stanley [STA]_.

            Let `\lambda` be a partition. Let `n` be a positive integer. If
            the `n`-core of `\lambda` is nonempty, then
            `\mathbf{V}_n(s_\lambda) = 0`. Otherwise, the following method
            computes `\mathbf{V}_n(s_\lambda)`: Write the partition `\lambda`
            in the form `(\lambda_1, \lambda_2, ..., \lambda_{ns})` for some
            nonnegative integer `s`. (If `n` does not divide the length of
            `\lambda`, then this is achieved by adding trailing zeroes to
            `\lambda`.) Set `\beta_i = \lambda_i + ns - i` for every
            `s \in \{ 1, 2, \ldots, ns \}`. Then,
            `(\beta_1, \beta_2, ..., \beta_{ns})` is a strictly decreasing
            sequence of nonnegative integers. Stably sort the list
            `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
            `-1 - \beta_i` modulo `n`. Let `\xi` be the sign of the
            permutation that is used for this sorting. Let `\psi` be the sign
            of the permutation that is used to stably sort the list
            `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
            `i - 1` modulo `n`. (Notice that `\psi = (-1)^{n(n-1)s(s-1)/4}`.)
            Then, `\mathbf{V}_n(s_\lambda) = \xi \psi \prod_{i=0}^{n-1}
            s_{\lambda^{(i)}}`, where
            `(\lambda^{(0)}, \lambda^{(1)}, \ldots, \lambda^{(n - 1)})`
            is the `n`-quotient of `\lambda`.

            INPUT:

            - ``n`` -- positive integer

            OUTPUT:

            The result of applying the `n`-th Verschiebung operator (on the
            ring of symmetric functions) to ``self``.

            EXAMPLES::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: e = Sym.e()
                sage: e[3].verschiebung(2)
                0
                sage: e[4].verschiebung(4)
                -e[1]

            The Verschiebung endomorphisms are multiplicative::

                sage: all( all( e(lam).verschiebung(2) * e(mu).verschiebung(2)
                ....:           == (e(lam) * e(mu)).verschiebung(2)
                ....:           for mu in Partitions(4) )
                ....:      for lam in Partitions(4) )
                True

            TESTS:

            Let us check that this method on the elementary basis gives the
            same result as the implementation in :mod:`sage.combinat.sf.sfa`
            on the complete homogeneous basis::

                sage: Sym = SymmetricFunctions(QQ)
                sage: e = Sym.e(); h = Sym.h()
                sage: all( h(e(lam)).verschiebung(3) == h(e(lam).verschiebung(3))
                ....:      for lam in Partitions(6) )
                True
                sage: all( e(h(lam)).verschiebung(2) == e(h(lam).verschiebung(2))
                ....:      for lam in Partitions(4) )
                True
            """
            parent = self.parent()
            e_coords_of_self = self.monomial_coefficients().items()
            dct = {Partition([i // n for i in lam]):
                   (-1) ** (sum(lam) - (sum(lam) // n)) * coeff
                   for (lam, coeff) in e_coords_of_self
                   if all( i % n == 0 for i in lam )}
            result_in_e_basis = parent._from_dict(dct)
            return parent(result_in_e_basis)

        def expand(self, n, alphabet='x'):
            """
            Expand the symmetric function ``self`` as a symmetric polynomial
            in ``n`` variables.

            INPUT:

            - ``n`` -- nonnegative integer

            - ``alphabet`` -- (default: ``'x'``) a variable for the expansion

            OUTPUT:

            A monomial expansion of ``self`` in the `n` variables
            labelled by ``alphabet``.

            EXAMPLES::

                sage: e = SymmetricFunctions(QQ).e()
                sage: e([2,1]).expand(3)
                x0^2*x1 + x0*x1^2 + x0^2*x2 + 3*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
                sage: e([1,1,1]).expand(2)
                x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3
                sage: e([3]).expand(2)
                0
                sage: e([2]).expand(3)
                x0*x1 + x0*x2 + x1*x2
                sage: e([3]).expand(4,alphabet='x,y,z,t')
                x*y*z + x*y*t + x*z*t + y*z*t
                sage: e([3]).expand(4,alphabet='y')
                y0*y1*y2 + y0*y1*y3 + y0*y2*y3 + y1*y2*y3
                sage: e([]).expand(2)
                1
                sage: e([]).expand(0)
                1
                sage: (3*e([])).expand(0)
                3
            """
            condition = lambda part: max(part) > n
            return self._expand(condition, n, alphabet)

        def principal_specialization(self, n=infinity, q=None):
            r"""
            Return the principal specialization of a symmetric function.

            The *principal specialization* of order `n` at `q`
            is the ring homomorphism `ps_{n,q}` from the ring of
            symmetric functions to another commutative ring `R`
            given by `x_i \mapsto q^{i-1}` for `i \in \{1,\dots,n\}`
            and `x_i \mapsto 0` for `i > n`.
            Here, `q` is a given element of `R`, and we assume that
            the variables of our symmetric functions are
            `x_1, x_2, x_3, \ldots`.
            (To be more precise, `ps_{n,q}` is a `K`-algebra
            homomorphism, where `K` is the base ring.)
            See Section 7.8 of [EnumComb2]_.

            The *stable principal specialization* at `q` is the ring
            homomorphism `ps_q` from the ring of symmetric functions
            to another commutative ring `R` given by
            `x_i \mapsto q^{i-1}` for all `i`.
            This is well-defined only if the resulting infinite sums
            converge; thus, in particular, setting `q = 1` in the
            stable principal specialization is an invalid operation.

            INPUT:

            - ``n`` -- (default: ``infinity``) a nonnegative integer or
              ``infinity``, specifying whether to compute the principal
              specialization of order ``n`` or the stable principal
              specialization.

            - ``q`` -- (default: ``None``) the value to use for `q`; the
              default is to create a ring of polynomials in ``q``
              (or a field of rational functions in ``q``) over the
              given coefficient ring.

            We use the formulas from Proposition 7.8.3 of [EnumComb2]_
            (using Gaussian binomial coefficients `\binom{u}{v}_q`):

            .. MATH::

                ps_{n,q}(e_\lambda) = \prod_i q^{\binom{\lambda_i}{2}} \binom{n}{\lambda_i}_q,

                ps_{n,1}(e_\lambda) = \prod_i \binom{n}{\lambda_i},

                ps_q(e_\lambda) = \prod_i q^{\binom{\lambda_i}{2}} / \prod_{j=1}^{\lambda_i} (1-q^j).

            EXAMPLES::

                sage: e = SymmetricFunctions(QQ).e()
                sage: x = e[3,1]
                sage: x.principal_specialization(3)
                q^5 + q^4 + q^3
                sage: x = 5*e[1,1,1] + 3*e[2,1] + 1
                sage: x.principal_specialization(3)
                5*q^6 + 18*q^5 + 36*q^4 + 44*q^3 + 36*q^2 + 18*q + 6

            By default, we return a rational functions in `q`.  Sometimes
            it is better to obtain an element of the symbolic ring::

                sage: x.principal_specialization(q=var("q"))                            # needs sage.symbolic
                -3*q/((q^2 - 1)*(q - 1)^2) - 5/(q - 1)^3 + 1

            TESTS::

                sage: e.zero().principal_specialization(3)
                0
            """
            from sage.combinat.q_analogues import q_binomial

            def get_variable(ring, name):
                try:
                    ring(name)
                except TypeError:
                    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                    return PolynomialRing(ring, name).gen()
                else:
                    raise ValueError("the variable %s is in the base ring, pass it explicitly" % name)

            if q is None:
                q = get_variable(self.base_ring(), "q")

            if q == 1:
                if n == infinity:
                    raise ValueError("the stable principal specialization at q=1 is not defined")
                f = lambda partition: prod(binomial(n, part) for part in partition)
            elif n == infinity:
                f = lambda partition: prod(q**binomial(part, 2)/prod((1-q**i)
                                                                     for i in range(1,part+1))
                                           for part in partition)
            else:
                f = lambda partition: prod(q**binomial(part, 2)*q_binomial(n, part, q=q)
                                           for part in partition)

            return self.parent()._apply_module_morphism(self, f, q.parent())

        def exponential_specialization(self, t=None, q=1):
            r"""
            Return the exponential specialization of a
            symmetric function (when `q = 1`), or the
            `q`-exponential specialization (when `q \neq 1`).

            The *exponential specialization* `ex` at `t` is a
            `K`-algebra homomorphism from the `K`-algebra of
            symmetric functions to another `K`-algebra `R`.
            It is defined whenever the base ring `K` is a
            `\QQ`-algebra and `t` is an element of `R`.
            The easiest way to define it is by specifying its
            values on the powersum symmetric functions to be
            `p_1 = t` and `p_n = 0` for `n > 1`.
            Equivalently, on the homogeneous functions it is
            given by `ex(h_n) = t^n / n!`; see Proposition 7.8.4 of
            [EnumComb2]_.

            By analogy, the `q`-exponential specialization is a
            `K`-algebra homomorphism from the `K`-algebra of
            symmetric functions to another `K`-algebra `R` that
            depends on two elements `t` and `q` of `R` for which
            the elements `1 - q^i` for all positive integers `i`
            are invertible.
            It can be defined by specifying its values on the
            complete homogeneous symmetric functions to be

            .. MATH::

                ex_q(h_n) = t^n / [n]_q!,

            where `[n]_q!` is the `q`-factorial.  Equivalently, for
            `q \neq 1` and a homogeneous symmetric function `f` of
            degree `n`, we have

            .. MATH::

                ex_q(f) = (1-q)^n t^n ps_q(f),

            where `ps_q(f)` is the stable principal specialization of `f`
            (see :meth:`principal_specialization`).
            (See (7.29) in [EnumComb2]_.)

            The limit of `ex_q` as `q \to 1` is `ex`.

            INPUT:

            - ``t`` -- (default: ``None``) the value to use for `t`.
              The default is to create a ring of polynomials in `t`.

            - ``q`` -- (default: `1`) the value to use for `q`.  If
              ``q`` is ``None``, then a ring (or fraction field) of
              polynomials in ``q`` is created.

            EXAMPLES::

                sage: e = SymmetricFunctions(QQ).e()
                sage: x = e[3,2]
                sage: x.exponential_specialization()
                1/12*t^5
                sage: x = 5*e[2] + 3*e[1] + 1
                sage: x.exponential_specialization(t=var("t"), q=var("q"))              # needs sage.symbolic
                5*q*t^2/(q + 1) + 3*t + 1

            TESTS::

                sage: e.zero().exponential_specialization()
                0
            """
            from sage.combinat.q_analogues import q_factorial

            def get_variable(ring, name):
                try:
                    ring(name)
                except TypeError:
                    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                    return PolynomialRing(ring, name).gen()
                else:
                    raise ValueError("the variable %s is in the base ring, pass it explicitly" % name)

            if q == 1:
                if t is None:
                    t = get_variable(self.base_ring(), 't')

                def f(partition):
                    n = 0
                    m = 1
                    for part in partition:
                        n += part
                        m *= factorial(part)
                    return t**n/m

                return self.parent()._apply_module_morphism(self, f, t.parent())

            if q is None and t is None:
                q = get_variable(self.base_ring(), 'q')
                t = get_variable(q.parent(), 't')
            elif q is None:
                q = get_variable(t.parent(), 'q')
            elif t is None:
                t = get_variable(q.parent(), 't')

            def f(partition):
                n = 0
                m = 1
                for part in partition:
                    n += part
                    m *= q**binomial(part, 2)/q_factorial(part, q=q)

                return t**n * m

            return self.parent()._apply_module_morphism(self, f, t.parent())


# Backward compatibility for unpickling
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.sf.elementary',
                           'SymmetricFunctionAlgebraElement_elementary',
                           SymmetricFunctionAlgebra_elementary.Element)
