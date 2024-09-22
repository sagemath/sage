# sage.doctest: needs sage.combinat sage.modules
r"""
Homogeneous symmetric functions

By this we mean the basis formed of the complete homogeneous
symmetric functions `h_\lambda`, not an arbitrary graded basis.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
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

####################################
#                                  #
# Homogeneous Symmetric Functions  #
#                                  #
####################################
from sage.arith.misc import binomial, factorial
from sage.combinat.partition import Partition
from sage.combinat.sf import multiplicative, classical
from sage.misc.misc_c import prod
from sage.rings.infinity import infinity


class SymmetricFunctionAlgebra_homogeneous(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, Sym):
        """
        A class of methods specific to the homogeneous basis of
        symmetric functions.

        INPUT:

        - ``self`` -- a homogeneous basis of symmetric functions
        - ``Sym`` -- an instance of the ring of symmetric functions

        TESTS::

            sage: h = SymmetricFunctions(QQ).e()
            sage: h == loads(dumps(h))
            True
            sage: TestSuite(h).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(h).run(elements = [h[1,1]+h[2], h[1]+2*h[1,1]])
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, Sym, "homogeneous", 'h')

    def _dual_basis_default(self):
        r"""
        Return the dual basis to ``self``.

        INPUT:

        - ``self`` -- a homogeneous basis of symmetric functions
        - ``scalar`` -- (optional) input which specifies a function ``zee``
          on partitions. The function ``zee`` determines the scalar
          product on the power sum basis with normalization
          `\langle p_\mu, p_\mu \rangle = \mathrm{zee}(mu)`.
          (default: uses standard ``zee`` function)
        - ``scalar_name`` -- specifies the name of the scalar function
          (optional)
        - ``prefix`` -- (optional) input which specifies the prefix to be
          used to display the basis

        OUTPUT:

        The dual basis of the homogeneous basis with respect to the
        standard scalar product (the monomial basis).  If a function
        ``zee`` is specified, the dual basis is with respect to the
        modified scalar product.

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: h = SymmetricFunctions(QQ).h()
            sage: h.dual_basis() == m
            True

            sage: zee = lambda x : 2
            sage: hh = h.dual_basis(zee); hh
            Dual basis to Symmetric Functions over Rational Field in the homogeneous basis
            sage: hh[2,1].scalar(h[2,1])
            1
            sage: hh[2,2].scalar(h[2,2])
            4

        TESTS::

            sage: h._dual_basis_default() is h.dual_basis()
            True
        """
        return self.realization_of().m()

    def coproduct_on_generators(self, i):
        r"""
        Return the coproduct on `h_i`.

        INPUT:

        - ``self`` -- a homogeneous basis of symmetric functions
        - ``i`` -- nonnegative integer

        OUTPUT: the sum `\sum_{r=0}^i h_r \otimes h_{i-r}`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: h = Sym.homogeneous()
            sage: h.coproduct_on_generators(2)
            h[] # h[2] + h[1] # h[1] + h[2] # h[]
            sage: h.coproduct_on_generators(0)
            h[] # h[]
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

            :meth:`omega_involution()` is a synonym for the :meth:`omega()`
            method.

            OUTPUT: the image of ``self`` under the omega automorphism

            EXAMPLES::

                sage: h = SymmetricFunctions(QQ).h()
                sage: a = h([2,1]); a
                h[2, 1]
                sage: a.omega()
                h[1, 1, 1] - h[2, 1]
                sage: e = SymmetricFunctions(QQ).e()
                sage: e(h([2,1]).omega())
                e[2, 1]
            """
            e = self.parent().realization_of().e()
            return self.parent()(e._from_element(self))

        omega_involution = omega

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

                sage: h = SymmetricFunctions(QQ).h()
                sage: h([3]).expand(2)
                x0^3 + x0^2*x1 + x0*x1^2 + x1^3
                sage: h([1,1,1]).expand(2)
                x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3
                sage: h([2,1]).expand(3)
                x0^3 + 2*x0^2*x1 + 2*x0*x1^2 + x1^3 + 2*x0^2*x2 + 3*x0*x1*x2 + 2*x1^2*x2 + 2*x0*x2^2 + 2*x1*x2^2 + x2^3
                sage: h([3]).expand(2,alphabet='y')
                y0^3 + y0^2*y1 + y0*y1^2 + y1^3
                sage: h([3]).expand(2,alphabet='x,y')
                x^3 + x^2*y + x*y^2 + y^3
                sage: h([3]).expand(3,alphabet='x,y,z')
                x^3 + x^2*y + x*y^2 + y^3 + x^2*z + x*y*z + y^2*z + x*z^2 + y*z^2 + z^3
                sage: (h([]) + 2*h([1])).expand(3)
                2*x0 + 2*x1 + 2*x2 + 1
                sage: h([1]).expand(0)
                0
                sage: (3*h([])).expand(0)
                3
            """
            if n == 0:   # Symmetrica crashes otherwise...
                return self.counit()
            condition = lambda part: False
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

                ps_{n,q}(h_\lambda) = \prod_i \binom{n+\lambda_i-1}{\lambda_i}_q,

                ps_{n,1}(h_\lambda) = \prod_i \binom{n+\lambda_i-1}{\lambda_i},

                ps_q(h_\lambda) = 1 / \prod_i \prod_{j=1}^{\lambda_i} (1-q^j).

            EXAMPLES::

                sage: h = SymmetricFunctions(QQ).h()
                sage: x = h[2,1]
                sage: x.principal_specialization(3)
                q^6 + 2*q^5 + 4*q^4 + 4*q^3 + 4*q^2 + 2*q + 1
                sage: x = 3*h[2] + 2*h[1] + 1
                sage: x.principal_specialization(3, q=var("q"))                         # needs sage.symbolic
                2*(q^3 - 1)/(q - 1) + 3*(q^4 - 1)*(q^3 - 1)/((q^2 - 1)*(q - 1)) + 1

            TESTS::

                sage: x = h.zero()
                sage: s = x.principal_specialization(3); s
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
                q = get_variable(self.base_ring(), 'q')

            if q == 1:
                if n == infinity:
                    raise ValueError("the stable principal specialization at q=1 is not defined")
                f = lambda partition: prod(binomial(n+part-1, part) for part in partition)
            elif n == infinity:
                f = lambda partition: prod(1/prod((1-q**i) for i in range(1, part+1)) for part in partition)
            else:
                f = lambda partition: prod(q_binomial(n+part-1, part, q=q) for part in partition)

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

            - ``t`` -- (default: ``None``) the value to use for `t`;
              the default is to create a ring of polynomials in ``t``

            - ``q`` -- (default: `1`) the value to use for `q`;  if
              ``q`` is ``None``, then a ring (or fraction field) of
              polynomials in ``q`` is created

            EXAMPLES::

                sage: h = SymmetricFunctions(QQ).h()
                sage: x = h[5,3]
                sage: x.exponential_specialization()
                1/720*t^8
                sage: factorial(5)*factorial(3)
                720

                sage: x = 5*h[1,1,1] + 3*h[2,1] + 1
                sage: x.exponential_specialization()
                13/2*t^3 + 1

            We also support the `q`-exponential_specialization::

                sage: factor(h[3].exponential_specialization(q=var("q"), t=var("t")))   # needs sage.symbolic
                t^3/((q^2 + q + 1)*(q + 1))

            TESTS::

                sage: x = h.zero()
                sage: s = x.exponential_specialization(); s
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
                    m *= q_factorial(part, q=q)
                return t**n/m

            return self.parent()._apply_module_morphism(self, f, t.parent())


# Backward compatibility for unpickling
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.sf.homogeneous',
                           'SymmetricFunctionAlgebraElement_homogeneous',
                           SymmetricFunctionAlgebra_homogeneous.Element)
