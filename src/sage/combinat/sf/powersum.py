# sage.doctest: needs sage.combinat sage.modules
"""
Power sum symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#                     2012 Anne Schilling <anne at math.ucdavis.edu>
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
from . import sfa, multiplicative, classical
from sage.combinat.partition import Partition
from sage.arith.misc import divisors
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod
from sage.misc.superseded import deprecated_function_alias


class SymmetricFunctionAlgebra_power(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, Sym):
        """
        A class for methods associated to the power sum basis of the symmetric functions.

        INPUT:

        - ``self`` -- the power sum basis of the symmetric functions
        - ``Sym`` -- an instance of the ring of symmetric functions

        TESTS::

            sage: p = SymmetricFunctions(QQ).p()
            sage: p == loads(dumps(p))
            True
            sage: TestSuite(p).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(p).run(elements = [p[1,1]+p[2], p[1]+2*p[1,1]])
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, Sym, "powersum", 'p')

    def coproduct_on_generators(self, i):
        r"""
        Return coproduct on generators for power sums `p_i`
        (for `i > 0`).

        The elements `p_i` are primitive elements.

        INPUT:

        - ``self`` -- the power sum basis of the symmetric functions
        - ``i`` -- positive integer

        OUTPUT: the result of the coproduct on the generator `p(i)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.powersum()
            sage: p.coproduct_on_generators(2)
            p[] # p[2] + p[2] # p[]
        """
        Pi = Partition([i])
        P0 = Partition([])
        T = self.tensor_square()
        return T.sum_of_monomials( [(Pi, P0), (P0, Pi)] )

    def antipode_on_basis(self, partition):
        r"""
        Return the antipode of ``self[partition]``.

        The antipode on the generator `p_i` (for `i > 0`) is `-p_i`,
        and the antipode on `p_\mu` is `(-1)^{length(\mu)} p_\mu`.

        INPUT:

        - ``self`` -- the power sum basis of the symmetric functions
        - ``partition`` -- a partition

        OUTPUT:

        - the result of the antipode on ``self(partition)``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.p()
            sage: p.antipode_on_basis([2])
            -p[2]
            sage: p.antipode_on_basis([3])
            -p[3]
            sage: p.antipode_on_basis([2,2])
            p[2, 2]
            sage: p.antipode_on_basis([])
            p[]
        """
        if len(partition) % 2 == 0:
            return self[partition]
        return -self[partition]
        #This is slightly faster than: return (-1)**len(partition) * self[partition]

    def bottom_schur_function(self, partition, degree=None):
        r"""
        Return the least-degree component of ``s[partition]``,
        where ``s`` denotes the Schur basis of the symmetric
        functions, and the grading is not the usual grading on the
        symmetric functions but rather the grading which gives
        every `p_i` degree `1`.

        This least-degree component has its degree equal to the
        Frobenius rank of ``partition``, while the degree with respect
        to the usual grading is still the size of ``partition``.

        This method requires the base ring to be a (commutative)
        `\QQ`-algebra. This restriction is unavoidable, since
        the least-degree component (in general) has noninteger
        coefficients in all classical bases of the symmetric
        functions.

        The optional keyword ``degree`` allows taking any
        homogeneous component rather than merely the least-degree
        one. Specifically, if ``degree`` is set, then the
        ``degree``-th component will be returned.

        REFERENCES:

        .. [ClSt03] Peter Clifford, Richard P. Stanley,
           *Bottom Schur functions*.
           :arxiv:`math/0311382v2`.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.p()
            sage: p.bottom_schur_function([2,2,1])
            -1/6*p[3, 2] + 1/4*p[4, 1]
            sage: p.bottom_schur_function([2,1])
            -1/3*p[3]
            sage: p.bottom_schur_function([3])
            1/3*p[3]
            sage: p.bottom_schur_function([1,1,1])
            1/3*p[3]
            sage: p.bottom_schur_function(Partition([1,1,1]))
            1/3*p[3]
            sage: p.bottom_schur_function([2,1], degree=1)
            -1/3*p[3]
            sage: p.bottom_schur_function([2,1], degree=2)
            0
            sage: p.bottom_schur_function([2,1], degree=3)
            1/3*p[1, 1, 1]
            sage: p.bottom_schur_function([2,2,1], degree=3)
            1/8*p[2, 2, 1] - 1/6*p[3, 1, 1]
        """
        from sage.combinat.partition import _Partitions
        s = self.realization_of().schur()
        partition = _Partitions(partition)
        if degree is None:
            degree = partition.frobenius_rank()
        s_partition = self(s[partition])
        return self.sum_of_terms([(p, coeff) for p, coeff
                                  in s_partition if len(p) == degree],
                                 distinct=True)

    def eval_at_permutation_roots_on_generators(self, k, rho):
        r"""
        Evaluate `p_k` at eigenvalues of permutation matrix.

        This function evaluates a symmetric function ``p([k])``
        at the eigenvalues of a permutation matrix with cycle
        structure ``\rho``.

        This function evaluates a `p_k` at the roots of unity

        .. MATH::

            \Xi_{\rho_1},\Xi_{\rho_2},\ldots,\Xi_{\rho_\ell}

        where

        .. MATH::

            \Xi_{m} = 1,\zeta_m,\zeta_m^2,\ldots,\zeta_m^{m-1}

        and `\zeta_m` is an `m` root of unity.
        This is characterized by `p_k[ A , B ] = p_k[A] + p_k[B]` and
        `p_k[ \Xi_m ] = 0` unless `m` divides `k` and `p_{rm}[\Xi_m]=m`.

        INPUT:

        - ``k`` -- nonnegative integer
        - ``rho`` -- a partition or a list of nonnegative integers

        OUTPUT: an element of the base ring

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: p.eval_at_permutation_roots_on_generators(3, [6])
            0
            sage: p.eval_at_permutation_roots_on_generators(3, [3])
            3
            sage: p.eval_at_permutation_roots_on_generators(3, [1])
            1
            sage: p.eval_at_permutation_roots_on_generators(3, [3,3])
            6
            sage: p.eval_at_permutation_roots_on_generators(3, [1,1,1,1,1])
            5
        """
        return self.base_ring().sum(d*list(rho).count(d) for d in divisors(k))

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

                sage: p = SymmetricFunctions(QQ).p()
                sage: a = p([2,1]); a
                p[2, 1]
                sage: a.omega()
                -p[2, 1]
                sage: p([]).omega()
                p[]
                sage: p(0).omega()
                0
                sage: p = SymmetricFunctions(ZZ).p()
                sage: (p([3,1,1]) - 2 * p([2,1])).omega()
                2*p[2, 1] + p[3, 1, 1]
            """
            f = lambda part, coeff: (part, (-1)**(sum(part)-len(part)) * coeff)
            return self.map_item(f)

        omega_involution = omega

        def scalar(self, x, zee=None):
            r"""
            Return the standard scalar product of ``self`` and ``x``.

            INPUT:

            - ``x`` -- a power sum symmetric function
            - ``zee`` -- (default: uses standard ``zee`` function) optional
              input specifying the scalar product on the power sum basis with
              normalization `\langle p_{\mu}, p_{\mu} \rangle =
              \mathrm{zee}(\mu)`. ``zee`` should be a function on partitions.

            Note that the power-sum symmetric functions are orthogonal under
            this scalar product. With the default value of ``zee``, the value
            of `\langle p_{\lambda}, p_{\lambda} \rangle` is given by the
            size of the centralizer in `S_n` of a permutation of cycle
            type `\lambda`.

            OUTPUT:

            - the standard scalar product between ``self`` and ``x``, or, if
              the optional parameter ``zee`` is specified, then the scalar
              product with respect to the normalization `\langle p_{\mu},
              p_{\mu} \rangle = \mathrm{zee}(\mu)` with the power sum basis
              elements being orthogonal

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: p4 = Partitions(4)
                sage: matrix([ [p(a).scalar(p(b)) for a in p4] for b in p4])
                [ 4  0  0  0  0]
                [ 0  3  0  0  0]
                [ 0  0  8  0  0]
                [ 0  0  0  4  0]
                [ 0  0  0  0 24]
                sage: p(0).scalar(p(1))
                0
                sage: p(1).scalar(p(2))
                2

                sage: zee = lambda x : 1
                sage: matrix( [[p[la].scalar(p[mu], zee) for la in Partitions(3)] for mu in Partitions(3)])
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            parent = self.parent()
            x = parent(x)
            if zee is None:
                f = lambda part1, part2: sfa.zee(part1)
            else:
                f = lambda part1, part2: zee(part1)
            return parent._apply_multi_module_morphism(self, x, f, orthogonal=True)

        def _derivative(self, part):
            """
            Return the 'derivative' of ``p([part])`` with respect to ``p([1])``
            (where ``p([part])`` is regarded as a polynomial in the
            indeterminates ``p([i])``).

            INPUT:

            - ``part`` -- a partition

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: a = p([2,1])
                sage: a._derivative(Partition([2,1]))
                p[2]
                sage: a._derivative(Partition([1,1,1]))
                3*p[1, 1]
            """
            p = self.parent()
            if 1 not in part:
                return p.zero()
            else:
                return len([i for i in part if i == 1]) * p(part[:-1])

        def _derivative_with_respect_to_p1(self):
            """
            Return the 'derivative' of a symmetric function in the power sum
            basis with respect to ``p([1])`` (where ``p([part])`` is regarded
            as a polynomial in the indeterminates ``p([i])``).

            On the Frobenius image of an `S_n`-module, the resulting character
            is the Frobenius image of the restriction of this module
            to `S_{n-1}`.

            OUTPUT:

            - a symmetric function (in the power sum basis) of degree one
              smaller than ``self``, obtained by differentiating ``self``
              by `p_1`.

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: a = p([2,1,1,1])
                sage: a._derivative_with_respect_to_p1()
                3*p[2, 1, 1]
                sage: a = p([3,2])
                sage: a._derivative_with_respect_to_p1()
                0
                sage: p(0)._derivative_with_respect_to_p1()
                0
                sage: p(1)._derivative_with_respect_to_p1()
                0
                sage: p([1])._derivative_with_respect_to_p1()
                p[]
                sage: f = p[1] + p[2,1]
                sage: f._derivative_with_respect_to_p1()
                p[] + p[2]
            """
            p = self.parent()
            if self == p.zero():
                return self
            return p._apply_module_morphism(self, self._derivative)

        def adams_operator(self, n):
            r"""
            Return the image of the symmetric function ``self`` under the
            `n`-th Adams operator.

            The `n`-th Adams operator `\mathbf{f}_n` is defined to be the
            map from the ring of symmetric functions to itself that sends
            every symmetric function `P(x_1, x_2, x_3, \ldots)` to
            `P(x_1^n, x_2^n, x_3^n, \ldots)`. This operator `\mathbf{f}_n`
            is a Hopf algebra endomorphism, and satisfies

            .. MATH::

                \mathbf{f}_n m_{(\lambda_1, \lambda_2, \lambda_3, \ldots)} =
                m_{(n\lambda_1, n\lambda_2, n\lambda_3, \ldots)}

            for every partition `(\lambda_1, \lambda_2, \lambda_3, \ldots)`
            (where `m` means the monomial basis). Moreover,
            `\mathbf{f}_n (p_r) = p_{nr}` for every positive integer `r` (where
            `p_k` denotes the `k`-th powersum symmetric function).

            The `n`-th Adams operator is also called the `n`-th
            Frobenius endomorphism. It is not related to the Frobenius map
            which connects the ring of symmetric functions with the
            representation theory of the symmetric group.

            The `n`-th Adams operator is the `n`-th Adams operator
            of the `\Lambda`-ring of symmetric functions over the integers.

            The `n`-th Adams operator can also be described via plethysm:
            Every symmetric function `P` satisfies
            `\mathbf{f}_n(P) = p_n \circ P = P \circ p_n`,
            where `p_n` is the `n`-th powersum symmetric function, and `\circ`
            denotes (outer) plethysm.

            INPUT:

            - ``n`` -- positive integer

            OUTPUT:

            The result of applying the `n`-th Adams operator (on the ring
            of symmetric functions) to ``self``.

            EXAMPLES::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: p = Sym.p()
                sage: p[3].adams_operator(2)
                p[6]
                sage: p[4,2,1].adams_operator(3)
                p[12, 6, 3]
                sage: p([]).adams_operator(4)
                p[]
                sage: p[3].adams_operator(1)
                p[3]
                sage: (p([3]) - p([2]) + p([])).adams_operator(3)
                p[] - p[6] + p[9]

            TESTS:

            Let us check that this method on the powersum basis gives the
            same result as the implementation in :mod:`sage.combinat.sf.sfa`
            on the complete homogeneous basis::

                sage: Sym = SymmetricFunctions(QQ)
                sage: p = Sym.p(); h = Sym.h()
                sage: all( h(p(lam)).adams_operator(3) == h(p(lam).adams_operator(3))
                ....:      for lam in Partitions(3) )
                True
                sage: all( p(h(lam)).adams_operator(2) == p(h(lam).adams_operator(2))
                ....:      for lam in Partitions(4) )
                True

            .. SEEALSO::

                :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.plethysm`
            """
            dct = {lam.stretch(n): coeff
                   for lam, coeff in self.monomial_coefficients().items()}
            return self.parent()._from_dict(dct)

        frobenius = deprecated_function_alias(36396, adams_operator)

        adams_operation = deprecated_function_alias(36396, adams_operator)

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
            Adams operator (see :meth:`adams_operator` for its definition)
            with respect to the Hall scalar product (:meth:`scalar`).

            The action of the `n`-th Verschiebung operator on the Schur basis
            can also be computed explicitly. The following (probably clumsier
            than necessary) description can be obtained by solving exercise
            7.61 in Stanley's [STA]_.

            Let `\lambda` be a partition. Let `n` be a positive integer. If
            the `n`-core of `\lambda` is nonempty, then
            `\mathbf{V}_n(s_\lambda) = 0`. Otherwise, the following method
            computes `\mathbf{V}_n(s_\lambda)`: Write the partition `\lambda`
            in the form `(\lambda_1, \lambda_2, \ldots, \lambda_{ns})` for some
            nonnegative integer `s`. (If `n` does not divide the length of
            `\lambda`, then this is achieved by adding trailing zeroes to
            `\lambda`.) Set `\beta_i = \lambda_i + ns - i` for every
            `s \in \{ 1, 2, \ldots, ns \}`. Then,
            `(\beta_1, \beta_2, \ldots, \beta_{ns})` is a strictly decreasing
            sequence of nonnegative integers. Stably sort the list
            `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
            `-1 - \beta_i` modulo `n`. Let `\xi` be the sign of the
            permutation that is used for this sorting. Let `\psi` be the sign
            of the permutation that is used to stably sort the list
            `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
            `i - 1` modulo `n`. (Notice that `\psi = (-1)^{n(n-1)s(s-1)/4}`.)
            Then, `\mathbf{V}_n(s_\lambda) = \xi \psi \prod_{i = 0}^{n - 1}
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
                sage: p = Sym.p()
                sage: p[3].verschiebung(2)
                0
                sage: p[4].verschiebung(4)
                4*p[1]

            The Verschiebung endomorphisms are multiplicative::

                sage: all( all( p(lam).verschiebung(2) * p(mu).verschiebung(2)
                ....:           == (p(lam) * p(mu)).verschiebung(2)
                ....:           for mu in Partitions(4) )
                ....:      for lam in Partitions(4) )
                True

            Testing the adjointness between the Adams operators
            `\mathbf{f}_n` and the Verschiebung operators
            `\mathbf{V}_n`::

                sage: Sym = SymmetricFunctions(QQ)
                sage: p = Sym.p()
                sage: all( all( p(lam).verschiebung(2).scalar(p(mu))
                ....:           == p(lam).scalar(p(mu).adams_operator(2))
                ....:           for mu in Partitions(2) )
                ....:      for lam in Partitions(4) )
                True

            TESTS:

            Let us check that this method on the powersum basis gives the
            same result as the implementation in :mod:`sage.combinat.sf.sfa`
            on the monomial basis::

                sage: Sym = SymmetricFunctions(QQ)
                sage: p = Sym.p(); m = Sym.m()
                sage: all( m(p(lam)).verschiebung(3) == m(p(lam).verschiebung(3))
                ....:      for lam in Partitions(6) )
                True
                sage: all( p(m(lam)).verschiebung(2) == p(m(lam).verschiebung(2))
                ....:      for lam in Partitions(4) )
                True
            """
            parent = self.parent()
            p_coords_of_self = self.monomial_coefficients().items()
            dct = {Partition([i // n for i in lam]): coeff * (n ** len(lam))
                   for (lam, coeff) in p_coords_of_self
                   if all( i % n == 0 for i in lam )}
            result_in_p_basis = parent._from_dict(dct)
            return parent(result_in_p_basis)

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

                sage: p = SymmetricFunctions(QQ).p()
                sage: a = p([2])
                sage: a.expand(2)
                x0^2 + x1^2
                sage: a.expand(3, alphabet=['a','b','c'])
                a^2 + b^2 + c^2
                sage: p([2,1,1]).expand(2)
                x0^4 + 2*x0^3*x1 + 2*x0^2*x1^2 + 2*x0*x1^3 + x1^4
                sage: p([7]).expand(4)
                x0^7 + x1^7 + x2^7 + x3^7
                sage: p([7]).expand(4,alphabet='t')
                t0^7 + t1^7 + t2^7 + t3^7
                sage: p([7]).expand(4,alphabet='x,y,z,t')
                x^7 + y^7 + z^7 + t^7
                sage: p(1).expand(4)
                1
                sage: p(0).expand(4)
                0
                sage: (p([]) + 2*p([1])).expand(3)
                2*x0 + 2*x1 + 2*x2 + 1
                sage: p([1]).expand(0)
                0
                sage: (3*p([])).expand(0)
                3
            """
            if n == 0:   # Symmetrica crashes otherwise...
                return self.counit()
            condition = lambda part: False
            return self._expand(condition, n, alphabet)

        def eval_at_permutation_roots(self, rho):
            r"""
            Evaluate at eigenvalues of a permutation matrix.

            Evaluate an element of the power sum basis at the eigenvalues
            of a permutation matrix with cycle structure `\rho`.

            This function evaluates an element at the roots of unity

            .. MATH::

                \Xi_{\rho_1},\Xi_{\rho_2},\ldots,\Xi_{\rho_\ell}

            where

            .. MATH::

                \Xi_{m} = 1,\zeta_m,\zeta_m^2,\ldots,\zeta_m^{m-1}

            and `\zeta_m` is an `m` root of unity.
            These roots of unity represent the eigenvalues of permutation
            matrix with cycle structure `\rho`.

            INPUT:

            - ``rho`` -- a partition or a list of nonnegative integers

            OUTPUT: an element of the base ring

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: p([3,3]).eval_at_permutation_roots([6])
                0
                sage: p([3,3]).eval_at_permutation_roots([3])
                9
                sage: p([3,3]).eval_at_permutation_roots([1])
                1
                sage: p([3,3]).eval_at_permutation_roots([3,3])
                36
                sage: p([3,3]).eval_at_permutation_roots([1,1,1,1,1])
                25
                sage: (p[1]+p[2]+p[3]).eval_at_permutation_roots([3,2])
                5
            """
            p = self.parent()
            R = self.base_ring()
            on_basis = lambda lam: R.prod(
                p.eval_at_permutation_roots_on_generators(k, rho) for k in lam)
            return p._apply_module_morphism(self, on_basis, R)

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

            We use the formulas from Proposition 7.8.3 of [EnumComb2]_:

            .. MATH::

                ps_{n,q}(p_\lambda) = \prod_i (1-q^{n\lambda_i}) / (1-q^{\lambda_i}),

                ps_{n,1}(p_\lambda) = n^{\ell(\lambda)},

                ps_q(p_\lambda) = 1 / \prod_i (1-q^{\lambda_i}),

            where `\ell(\lambda)` denotes the length of `\lambda`,
            and where the products range from `i=1` to `i=\ell(\lambda)`.

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: x = p[8,7,3,1]
                sage: x.principal_specialization(3, q=var("q"))                         # needs sage.symbolic
                (q^24 - 1)*(q^21 - 1)*(q^9 - 1)/((q^8 - 1)*(q^7 - 1)*(q - 1))

                sage: x = 5*p[1,1,1] + 3*p[2,1] + 1
                sage: x.principal_specialization(3, q=var("q"))                         # needs sage.symbolic
                5*(q^3 - 1)^3/(q - 1)^3 + 3*(q^6 - 1)*(q^3 - 1)/((q^2 - 1)*(q - 1)) + 1

            By default, we return a rational function in `q`::

                sage: x.principal_specialization(3)
                8*q^6 + 18*q^5 + 36*q^4 + 38*q^3 + 36*q^2 + 18*q + 9

            If ``n`` is not given we return the stable principal specialization::

                sage: x.principal_specialization(q=var("q"))                            # needs sage.symbolic
                3/((q^2 - 1)*(q - 1)) - 5/(q - 1)^3 + 1

            TESTS::

                sage: p.zero().principal_specialization(3)
                0
            """
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
                f = lambda partition: n**len(partition)
            elif n == infinity:
                f = lambda partition: prod(1/(1-q**part) for part in partition)
            else:
                from sage.rings.integer_ring import ZZ
                ZZq = PolynomialRing(ZZ, "q")
                q_lim = ZZq.gen()

                def f(partition):
                    denom = prod((1 - q**part) for part in partition)
                    try:
                        ~denom
                        rational = prod((1 - q**(n*part)) for part in partition) / denom
                        return q.parent()(rational)
                    except (ZeroDivisionError, NotImplementedError, TypeError):
                        # If denom is not invertible, we need to do the
                        # computation with universal coefficients instead:
                        quotient = ZZq(prod((1-q_lim**(n*part))/(1-q_lim**part) for part in partition))
                        return quotient.subs({q_lim: q})

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

            - ``q`` -- (default: `1`) the value to use for `q`.  If
              ``q`` is ``None``, then a ring (or fraction field) of
              polynomials in ``q`` is created.

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: x = p[8,7,3,1]
                sage: x.exponential_specialization()
                0
                sage: x = p[3] + 5*p[1,1] + 2*p[1] + 1
                sage: x.exponential_specialization(t=var("t"))                          # needs sage.symbolic
                5*t^2 + 2*t + 1

            We also support the `q`-exponential_specialization::

                sage: factor(p[3].exponential_specialization(q=var("q"), t=var("t")))   # needs sage.symbolic
                (q - 1)^2*t^3/(q^2 + q + 1)

            TESTS::

                sage: p.zero().exponential_specialization()
                0
            """
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
                    for part in partition:
                        if part != 1:
                            return 0
                        n += 1
                    return t**n

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
                    m *= 1-q**part
                return (1-q)**n * t**n / m

            return self.parent()._apply_module_morphism(self, f, t.parent())


# Backward compatibility for unpickling
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.sf.powersum',
                           'SymmetricFunctionAlgebraElement_power',
                           SymmetricFunctionAlgebra_power.Element)
