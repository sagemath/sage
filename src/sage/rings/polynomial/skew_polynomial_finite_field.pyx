# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.combinat sage.rings.finite_rings
r"""
Univariate dense skew polynomials over finite fields

This module provides the
class:`~sage.rings.polynomial.skew_polynomial_finite_field.SkewPolynomial_finite_field_dense`,
which constructs a single univariate skew polynomial over a finite field
equipped with the Frobenius endomorphism. Among other things, it implements
the fast factorization algorithm designed in [CL2017]_.

AUTHOR::

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests and
  refactored classes and methods
"""

# ***************************************************************************
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.element cimport parent
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix2 import NotFullRankError

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.skew_polynomial_finite_order cimport SkewPolynomial_finite_order_dense

from sage.structure.factorization import Factorization
from sage.misc.mrange import xmrange_iter
from sage.misc.prandom import sample
from sage.arith.misc import factorial
from sage.combinat.q_analogues import q_jordan


cdef class SkewPolynomial_finite_field_dense(SkewPolynomial_finite_order_dense):
    cdef inline _reduced_norm_factored(self):
        """
        Return the reduced norm of this polynomial factorized in the center.
        """
        if self._norm_factor is None:
            N = self._parent._working_center(self.reduced_norm(var=False))
            self._norm_factor = N.factor()
        return self._norm_factor

    def is_irreducible(self):
        r"""
        Return ``True`` if this skew polynomial is irreducible.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: a = x^2 + t*x + 1
            sage: a.is_irreducible()
            False
            sage: a.factor()
            (x + 4*t^2 + 4*t + 1) * (x + 3*t + 2)

            sage: a = x^2 + t*x + t + 1
            sage: a.is_irreducible()
            True
            sage: a.factor()
            x^2 + t*x + t + 1

        Skew polynomials of degree `1` are of course irreducible::

            sage: a = x + t
            sage: a.is_irreducible()
            True

        A random irreducible skew polynomial is irreducible::

            sage: a = S.random_irreducible(degree=4,monic=True); a   # random
            x^4 + (t + 1)*x^3 + (3*t^2 + 2*t + 3)*x^2 + 3*t*x + 3*t
            sage: a.is_irreducible()
            True

        By convention, constant skew polynomials are not irreducible::

            sage: S(1).is_irreducible()
            False
            sage: S(0).is_irreducible()
            False
        """
        if self._norm_factor is not None:
            return len(self._norm_factor) == 1 and self._norm_factor[0][1] == 1
        N = self._parent._working_center(self.reduced_norm(var=False))
        return N.is_irreducible()

    def type(self, N):
        r"""
        Return the `N`-type of this skew polynomial (see definition below).

        INPUT:

        - ``N`` -- an irreducible polynomial in the
          center of the underlying skew polynomial ring

        .. NOTE::

            The result is cached.

        DEFINITION:

        The `N`-type of a skew polynomial `a` is the Partition
        `(t_0, t_1, t_2, \ldots)` defined by

        .. MATH::

            t_0 + \cdots + t_i = \frac{\deg gcd(a,N^i)}{\deg N},

        where `\deg N` is the degree of `N` considered as an
        element in the center.

        This notion has an important mathematical interest because
        it corresponds to the Jordan type of the `N`-typical part
        of the associated Galois representation.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center(); x3 = Z.gen()

            sage: a = x^4 + x^3 + (4*t^2 + 4)*x^2 + (t^2 + 2)*x + 4*t^2
            sage: N = x3^2 + x3 + 1
            sage: a.type(N)
            [1]
            sage: N = x3 + 1
            sage: a.type(N)
            [2]

            sage: a = x^3 + (3*t^2 + 1)*x^2 + (3*t^2 + t + 1)*x + t + 1
            sage: N = x3 + 1
            sage: a.type(N)
            [2, 1]

        If `N` does not divide the reduced map of `a`, the type
        is empty::

            sage: N = x3 + 2
            sage: a.type(N)
            []

        If `a = N`, the type is just `[r]` where `r` is the order
        of the twisting morphism ``Frob``::

            sage: N = x3^2 + x3 + 1
            sage: S(N).type(N)
            [3]

        `N` must be irreducible::

            sage: N = (x3 + 1) * (x3 + 2)
            sage: a.type(N)
            Traceback (most recent call last):
            ...
            ValueError: N is not irreducible
        """
        cdef SkewPolynomial_finite_field_dense P, d, NS
        cdef Py_ssize_t deg, degN, m
        cdef list type
        skew_ring = self._parent
        if parent(N) is not skew_ring._working_center:
            N = skew_ring._working_center(N)
        if self._types is None:
            self._types = dict()
        elif N in self._types:
            return self._types[N]
        if not N.is_irreducible():
            raise ValueError("N is not irreducible")
        if self._norm_factor is None:
            m = -1
        else:
            for n, m in self._norm_factor:
                if n == N:
                    break
        NS = skew_ring(N)
        type = [ ]
        degN = N.degree()
        P = self
        d = P.right_gcd(NS)
        deg = d.degree() // degN
        while deg > 0:
            if m >= 0:
                if deg == 1:
                    type += m * [1]
                    break
                m -= deg
            P = P // d
            type.append(deg)
            d = P.right_gcd(NS)
            deg = d.degree() // degN
        self._types[N] = type
        return type

    # Finding divisors
    # ----------------

    cdef SkewPolynomial_finite_field_dense _rdivisor_c(self, N):
        r"""
        Return a right divisor of this skew polynomial whose
        reduced norm is `N`.

        .. WARNING::

            `N` needs to be an irreducible factor of the
            reduced norm. This function does not check
            this (and his behaviour is not defined if the
            require property doesn't hold).
        """
        cdef skew_ring = self._parent
        cdef SkewPolynomial_finite_field_dense NS = <SkewPolynomial_finite_field_dense>skew_ring(N)
        cdef SkewPolynomial_finite_field_dense P = self.right_gcd(NS)
        cdef Py_ssize_t d = N.degree()
        cdef Py_ssize_t e = P.degree() // d
        if e == 1:
            return P.right_monic()

        cdef SkewPolynomial_finite_field_dense D
        cdef SkewPolynomial_finite_field_dense Q = <SkewPolynomial_finite_field_dense>(NS // P)
        cdef SkewPolynomial_finite_field_dense R, X
        cdef Py_ssize_t i, j, t, r = skew_ring._order
        cdef Polynomial dd, xx, yy, zz
        cdef Integer exp
        cdef list lM, lV
        cdef bint char2

        center = parent(N)
        E = center.quo(N)
        PE = PolynomialRing(E, name='T')
        char2 = skew_ring.characteristic() != 2
        if not char2:
            exp = <Integer>((E.cardinality()-1) // 2)
        while True:
            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((e*r-1, e*r-1))
            R = Q*R
            X = <SkewPolynomial_finite_field_dense>Q._new_c(list(Q._coeffs), Q._parent)
            lM = [ None ] * (e**2)
            for j in range(e):
                for i in range(e):
                    coeffs = [skew_ring._retraction(X[t*r+i]) for t in range(d)]
                    value = E(coeffs)
                    lM[i*e+j] = value
                X = (R*X) % NS
            M = MatrixSpace(E,e,e)(lM)
            V = MatrixSpace(E,e,1)([ E([skew_ring._retraction(X[t*r+i]) for t in range(d)]) for i in range(e) ])
            try:
                W = M._solve_right_nonsingular_square(V)
            except NotFullRankError:
                skew_ring._new_retraction_map()
                continue
            xx = PE(<list>W.list() + [E(-1)])
            if char2:
                zz = yy = PE.gen()
                for i in range(1,d):
                    zz = (zz*zz) % xx
                    yy += zz
                dd = xx.gcd(yy)
                if dd.degree() != 1:
                    continue
            else:
                yy = PE.gen().__pow__(exp,xx) - 1
                dd = xx.gcd(yy)
                if dd.degree() != 1:
                    yy += 2
                    dd = xx.gcd(yy)
                    if dd.degree() != 1:
                        continue
            D = P.right_gcd(R + skew_ring(center((dd[0]/dd[1]).list())))
            if D.degree() == 0:
                continue
            return D

    def _reduced_norm_factor_uniform(self):
        r"""
        Return a factor of the reduced norm of this skew
        polynomial, the probability of a given factor to
        show up being proportional to the number of irreducible
        divisors of ``self`` having this norm.

        This method is an helper function; it is not supposed
        to be called directly.

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: Z.<z> = S.center()

        In the two examples below, the reduced norms of `P` and `Q`
        have a unique factor::

            sage: P = x + a
            sage: P.reduced_norm().factor()
            z + 2
            sage: P._reduced_norm_factor_uniform()
            z + 2

            sage: Q = x^3 + 1
            sage: Q.reduced_norm().factor()
            (z + 1)^3
            sage: Q._reduced_norm_factor_uniform()
            z + 1

        Now, we consider the product `R = P \cdot Q`; it admits `32`
        irreducible divisors but among them, only one has norm `z + 2`,
        the others having norm `z + 1`.
        Therefore this method outputs `z + 2` with probability `1 / 32`
        and `z + 1` with probability `31 / 32`::

            sage: R = P*Q
            sage: counts = { z+1: 0, z+2: 0 }
            sage: for _ in range(1000):
            ....:     N = R._reduced_norm_factor_uniform()
            ....:     counts[N] += 1
            sage: counts   # random
            {z + 1: 969, z + 2: 31}
        """
        skew_ring = self._parent
        F = self._reduced_norm_factored()
        center = F[0][0].parent()
        cardcenter = center.base_ring().cardinality()
        gencenter = center.gen()
        count = [ ]
        total = 0
        for n, _ in F:
            if n == gencenter:
                total += 1
            else:
                degn = n.degree()
                P = self.right_gcd(skew_ring(n))
                m = P.degree() // degn
                cardL = cardcenter**degn
                total += (cardL**m - 1) / (cardL - 1)
            count.append(total)
        random = ZZ.random_element(total)
        for i in range(len(F)):
            if random < count[i]:
                return F[i][0]

    def _irreducible_divisors(self, bint right):
        r"""
        Return an iterator over all irreducible monic
        divisors of this skew polynomial.

        Do not use this function. Use instead
        :meth:`right_irreducible_divisors` and
        :meth:`left_irreducible_divisors`.

        INPUT:

        - ``right`` -- boolean; if ``True``, return right divisors,
          otherwise, return left divisors

        TESTS::

            sage: k.<a> = GF(7^4)
            sage: Frob = k.frobenius_endomorphism(2)
            sage: S.<x> = k['x', Frob]

            sage: P = S.random_element(degree=10)
            sage: rightdiv = [ f for f in P.right_irreducible_divisors() ]   # indirect doctest
            sage: len(rightdiv) == P.count_irreducible_divisors()
            True
            sage: len(rightdiv) == Set(rightdiv).cardinality()  # check no duplicates
            True
            sage: for D in rightdiv:
            ....:     assert P.is_right_divisible_by(D), "not right divisible"
            ....:     assert D.is_irreducible(), "not irreducible"

            sage: P = S.random_element(degree=10)
            sage: leftdiv = [ f for f in P.left_irreducible_divisors() ]   # indirect doctest
            sage: len(leftdiv) == P.count_irreducible_divisors()
            True
            sage: len(leftdiv) == Set(leftdiv).cardinality()  # check no duplicates
            True
            sage: for D in leftdiv:
            ....:     assert P.is_left_divisible_by(D), "not left divisible"
            ....:     assert D.is_irreducible(), "not irreducible"
        """
        cdef SkewPolynomial_finite_field_dense NS, P, Q, R, P1, Q1, L, V, g, d
        cdef Py_ssize_t i, m, degrandom
        if not self:
            return
        if right:
            quo_rem = SkewPolynomial_finite_field_dense.right_quo_rem
            quo_rem2 = SkewPolynomial_finite_field_dense.left_quo_rem
            gcd = SkewPolynomial_finite_field_dense.right_gcd
            gcd2 = SkewPolynomial_finite_field_dense.left_gcd
            mul = SkewPolynomial_finite_field_dense.__mul__  # why _mul_ doesn't work?
        else:
            quo_rem = SkewPolynomial_finite_field_dense.left_quo_rem
            quo_rem2 = SkewPolynomial_finite_field_dense.right_quo_rem
            gcd = SkewPolynomial_finite_field_dense.left_gcd
            gcd2 = SkewPolynomial_finite_field_dense.right_gcd
            mul = mul_op
        skew_ring = self._parent
        center = skew_ring._working_center
        kfixed = center.base_ring()
        F = self._reduced_norm_factored()
        for N, _ in F:
            if N == center.gen():
                yield skew_ring.gen()
                continue
            degN = N.degree()
            NS = skew_ring(N)
            P = gcd(self, NS)
            m = P.degree()/degN
            if m == 1:
                yield P
                continue
            Q,_ = quo_rem(NS, P)
            if right:
                P1 = self._rdivisor_c(N)
            else:
                D = self._rdivisor_c(N)
                degrandom = NS.degree() - 1
                while True:
                    P1 = P // P.right_gcd(NS // D)
                    if P1.degree() == degN:
                        break
                    while True:
                        R = skew_ring.random_element((degrandom,degrandom))
                        if NS.right_gcd(R) == 1:
                            break
                    D = NS.right_gcd(D*R)
            Q1,_ = quo_rem(P, P1)
            degrandom = P.degree() - 1
            while True:
                R = skew_ring.random_element((degrandom, degrandom))
                _, g = quo_rem2(mul(R,Q), P)
                if gcd2(g,P) != 1:
                    continue
                L = Q1
                V = L
                for i in range(1,m):
                    L = gcd2(mul(g,L), P)
                    V = gcd2(V,L)
                if V == 1:
                    break
            rng = xmrange_iter([kfixed]*degN, center)
            for i in range(m):
                for pol in xmrange_iter([rng]*i):
                    f = skew_ring.one()
                    for j in range(i):
                        coeff = pol.pop()
                        _, f = quo_rem2(g*f + coeff, P)
                    d = gcd2(mul(f,Q1), P)
                    d, _ = quo_rem2(P, d)
                    yield d

    def right_irreducible_divisor(self, uniform=False):
        r"""
        Return a right irreducible divisor of this skew polynomial.

        INPUT:

        - ``uniform`` -- boolean (default: ``False``); whether the
          output irreducible divisor should be uniformly distributed
          among all possibilities

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^6 + 3*t*x^5 + (3*t + 1)*x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + 2)*x + 4*t^2 + 3*t + 3

            sage: dr = a.right_irreducible_divisor(); dr  # random
            x^3 + (2*t^2 + t + 4)*x^2 + (4*t + 1)*x + 4*t^2 + t + 1
            sage: a.is_right_divisible_by(dr)
            True

        Right divisors are cached. Hence, if we ask again for a
        right divisor, we will get the same answer::

            sage: a.right_irreducible_divisor()  # random
            x^3 + (2*t^2 + t + 4)*x^2 + (4*t + 1)*x + 4*t^2 + t + 1

        However the algorithm is probabilistic. Hence, if we first
        reinitialize `a`, we may get a different answer::

            sage: a = x^6 + 3*t*x^5 + (3*t + 1)*x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + 2)*x + 4*t^2 + 3*t + 3
            sage: a.right_irreducible_divisor()  # random
            x^3 + (t^2 + 3*t + 4)*x^2 + (t + 2)*x + 4*t^2 + t + 1

        We can also generate uniformly distributed irreducible monic
        divisors as follows::

            sage: a.right_irreducible_divisor(uniform=True)  # random
            x^3 + (4*t + 2)*x^2 + (2*t^2 + 2*t + 2)*x + 2*t^2 + 2
            sage: a.right_irreducible_divisor(uniform=True)  # random
            x^3 + (t^2 + 2)*x^2 + (3*t^2 + 1)*x + 4*t^2 + 2*t
            sage: a.right_irreducible_divisor(uniform=True)  # random
            x^3 + x^2 + (4*t^2 + 2*t + 4)*x + t^2 + 3

        By convention, the zero skew polynomial has no irreducible
        divisor::

            sage: S(0).right_irreducible_divisor()
            Traceback (most recent call last):
            ...
            ValueError: 0 has no irreducible divisor
        """
        if self.is_zero():
            raise ValueError("0 has no irreducible divisor")
        skew_ring = self._parent
        if uniform:
            N = self._reduced_norm_factor_uniform()
        else:
            N = self._reduced_norm_factored()[0][0]
        NS = skew_ring(N)
        degN = N.degree()
        D = self._rdivisor_c(N)
        if uniform:
            P1 = self.right_gcd(NS)
            if P1.degree() != degN:
                Q1 = NS // P1
                deg = P1.degree() - 1
                while True:
                    R = Q1 * skew_ring.random_element((deg,deg))
                    if P1.right_gcd(R) == 1:
                        break
                D = P1.right_gcd(D*R)
        return D

    def left_irreducible_divisor(self, uniform=False):
        r"""
        Return a left irreducible divisor of this skew polynomial.

        INPUT:

        - ``uniform`` -- boolean (default: ``False``); whether the
          output irreducible divisor should be uniformly distributed
          among all possibilities

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^6 + 3*t*x^5 + (3*t + 1)*x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + 2)*x + 4*t^2 + 3*t + 3
            sage: dl = a.left_irreducible_divisor(); dl  # random
            x^3 + (t^2 + t + 2)*x^2 + (t + 2)*x + 3*t^2 + t + 4
            sage: a.is_left_divisible_by(dl)
            True

        The algorithm is probabilistic. Hence, if we ask again for
        a left irreducible divisor of `a`, we may get a different
        answer::

            sage: a.left_irreducible_divisor()  # random
            x^3 + (4*t + 3)*x^2 + (2*t^2 + 3*t + 4)*x + 4*t^2 + 2*t

        We can also generate uniformly distributed irreducible monic
        divisors as follows::

            sage: a.left_irreducible_divisor(uniform=True)  # random
            x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + t + 3)*x + 2*t^2 + 3
            sage: a.left_irreducible_divisor(uniform=True)  # random
            x^3 + (2*t^2 + t + 4)*x^2 + (2*t^2 + 4*t + 4)*x + 2*t + 3
            sage: a.left_irreducible_divisor(uniform=True)  # random
            x^3 + (t^2 + t + 2)*x^2 + (3*t^2 + t)*x + 2*t + 1

        By convention, the zero skew polynomial has no irreducible
        divisor::

            sage: S(0).left_irreducible_divisor()
            Traceback (most recent call last):
            ...
            ValueError: 0 has no irreducible divisor
        """
        if self.is_zero():
            raise ValueError("0 has no irreducible divisor")
        if uniform:
            N = self._reduced_norm_factor_uniform()
        else:
            N = self._reduced_norm_factored()[0][0]
        skew_ring = self._parent
        NS = skew_ring(N)
        degN = N.degree()
        D = self._rdivisor_c(N)
        deg = NS.degree() - 1
        P1 = self.left_gcd(NS)
        if not uniform:
            LD = P1 // P1.right_gcd(NS // D)
            if LD.degree() == degN:
                return LD
        while True:
            while True:
                R = skew_ring.random_element((deg,deg))
                if NS.right_gcd(R) == 1:
                    break
            D = NS.right_gcd(D*R)
            LD = P1 // P1.right_gcd(NS // D)
            if LD.degree() == degN:
                return LD

    def right_irreducible_divisors(self):
        r"""
        Return an iterator over all irreducible monic right divisors
        of this skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + 2*t*x^3 + 3*t^2*x^2 + (t^2 + t + 1)*x + 4*t + 3
            sage: iter = a.right_irreducible_divisors(); iter
            <...generator object at 0x...>
            sage: next(iter)   # random
            x + 2*t^2 + 4*t + 4
            sage: next(iter)   # random
            x + 3*t^2 + 4*t + 1

        We can use this function to build the list of all monic
        irreducible divisors of `a`::

            sage: rightdiv = [ d for d in a.right_irreducible_divisors() ]

        Note that the algorithm is probabilistic. As a consequence, if we
        build again the list of right monic irreducible divisors of `a`, we
        may get a different ordering::

            sage: rightdiv2 = [ d for d in a.right_irreducible_divisors() ]
            sage: rightdiv == rightdiv2
            False
            sage: Set(rightdiv) == Set(rightdiv2)
            True
        """
        return self._irreducible_divisors(True)

    def left_irreducible_divisors(self):
        r"""
        Return an iterator over all irreducible monic left divisors
        of this skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + 2*t*x^3 + 3*t^2*x^2 + (t^2 + t + 1)*x + 4*t + 3
            sage: iter = a.left_irreducible_divisors(); iter
            <...generator object at 0x...>
            sage: next(iter)  # random
            x + 3*t + 3
            sage: next(iter)  # random
            x + 4*t + 2

        We can use this function to build the list of all monic
        irreducible divisors of `a`::

            sage: leftdiv = [ d for d in a.left_irreducible_divisors() ]

        Note that the algorithm is probabilistic. As a consequence, if we
        build again the list of left monic irreducible divisors of `a`, we
        may get a different ordering::

            sage: leftdiv2 = [ d for d in a.left_irreducible_divisors() ]
            sage: Set(leftdiv) == Set(leftdiv2)
            True
        """
        return self._irreducible_divisors(False)

    def count_irreducible_divisors(self):
        r"""
        Return the number of irreducible monic divisors of
        this skew polynomial.

        .. NOTE::

            One can prove that there are always as many left
            irreducible monic divisors as right irreducible
            monic divisors.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

        We illustrate that a skew polynomial may have a number of irreducible
        divisors greater than its degree::

            sage: a = x^4 + (4*t + 3)*x^3 + t^2*x^2 + (4*t^2 + 3*t)*x + 3*t
            sage: a.count_irreducible_divisors()
            12

        We illustrate that an irreducible polynomial in the center have
        in general a lot of irreducible divisors in the skew polynomial
        ring::

            sage: Z.<x3> = S.center()
            sage: N = x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3; N
            x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3
            sage: N.is_irreducible()
            True
            sage: S(N).count_irreducible_divisors()
            9768751
        """
        if self.is_zero():
            return 0
        skew_ring = self.parent()
        cardcenter = skew_ring._working_center.base_ring().cardinality()
        gencenter = skew_ring._working_center.gen()
        F = self._reduced_norm_factored()
        val = self.valuation()
        self >>= val
        count = 0
        if val > 0:
            count = 1
        for N,_ in F:
            if N == gencenter:
                continue
            degN = N.degree()
            P = self.right_gcd(skew_ring(N))
            m = P.degree() // degN
            cardL = cardcenter**degN
            count += (cardL**m - 1) // (cardL - 1)
        return count

    # Finding factorizations
    # ----------------------

    cdef _factor_c(self):
        r"""
        Compute a factorization of ``self``.

        This is the low level implementation of :meth:`factor`.
        """
        cdef skew_ring = self._parent
        cdef SkewPolynomial_finite_field_dense poly = self.right_monic()
        cdef list a = poly._coeffs
        cdef Py_ssize_t val = 0
        while not a[0]:
            del a[0]
            val += 1

        cdef Py_ssize_t degQ, degrandom, m, mP, i
        cdef N
        cdef list factors = [ (skew_ring.gen(), val) ]
        cdef SkewPolynomial_finite_field_dense P, Q, P1, NS, g, right, Pn
        cdef unit = self.leading_coefficient()
        cdef Polynomial gencenter = skew_ring._working_center.gen()
        cdef Py_ssize_t p = skew_ring.characteristic()
        cdef F = self._reduced_norm_factored()

        for N, m in F:
            if N == gencenter:
                continue
            degN = N.degree()
            if poly.degree() == degN:
                factors.append((poly, 1))
                break
            NS = <SkewPolynomial_finite_field_dense>skew_ring(N)
            P1 = None
            while True:
                P = <SkewPolynomial_finite_field_dense>poly.right_gcd(NS)
                mP = P.degree() / degN
                if mP == 0: break
                if mP == 1:
                    factors.append((P, 1))
                    poly = poly // P
                    for i in range(1, m):
                        if poly.degree() == degN:
                            factors.append((poly,1))
                            break
                        P = poly.right_gcd(NS)
                        factors.append((P, 1))
                        poly = poly // P
                    break
                if P1 is None:
                    P1 = P._rdivisor_c(N)
                Q = <SkewPolynomial_finite_field_dense>NS._new_c(NS._coeffs[:], NS._parent)
                Q = P1 * (Q // P)
                factors.append((P1, 1))
                right = <SkewPolynomial_finite_field_dense>P1._new_c(P1._coeffs[:], P1._parent)
                m -= (mP-1)
                degrandom = P.degree()
                while mP > 2:
                    while True:
                        g = <SkewPolynomial_finite_field_dense>skew_ring.random_element((degrandom, degrandom))
                        g = (Q*g).right_gcd(P)
                        Pn = right._left_lcm_cofactor(g)
                        if Pn.degree() == degN: break
                    Pn = Pn.right_monic()
                    factors.append((Pn, 1))
                    right = Pn * right
                    degrandom -= degN
                    mP -= 1
                poly = poly // right
                P1, _ = P.right_quo_rem(right)
        factors.reverse()
        return Factorization(factors, sort=False, unit=unit)

    cdef _factor_uniform_c(self):
        r"""
        Compute a uniformly distributed factorization of ``self``.

        This is the low level implementation of :meth:`factor`.
        """
        skew_ring = self._parent
        cdef Integer cardE, cardcenter = skew_ring._working_center.base_ring().cardinality()
        cdef gencenter = skew_ring._working_center.gen()
        cdef SkewPolynomial_finite_field_dense gen = <SkewPolynomial_finite_field_dense>skew_ring.gen()

        cdef list factorsN = [ ]
        cdef dict dict_divisor = { }
        cdef dict dict_type = { }
        cdef dict dict_right = { }
        cdef Py_ssize_t m
        cdef list type

        for N, m in self._reduced_norm_factored():
            factorsN += m * [N]
            if N == gencenter:
                continue
            type = list(self.type(N))
            dict_type[N] = type
            if <Py_ssize_t>(type[0]) > 1:
                dict_divisor[N] = self._rdivisor_c(N)
                dict_right[N] = skew_ring.one()
        m = len(factorsN)
        cdef list indices = <list>sample(range(1, m+1), m)

        cdef unit = self.leading_coefficient()
        cdef SkewPolynomial_finite_field_dense left = self._new_c(self._coeffs[:],skew_ring)
        left = left.right_monic()
        cdef SkewPolynomial_finite_field_dense right = <SkewPolynomial_finite_field_dense>skew_ring.one()
        cdef SkewPolynomial_finite_field_dense L, R
        cdef SkewPolynomial_finite_field_dense NS, P, Q, D, D1, D2, d
        cdef list factors = [ ]
        cdef list maxtype
        cdef Py_ssize_t i, j, degN, deg
        cdef count, maxcount

        for i in indices:
            N = factorsN[i-1]
            if N == gencenter:
                D1 = gen
            else:
                type = dict_type[N]
                NS = skew_ring(N)
                P = left.right_gcd(NS)
                if <Py_ssize_t>(type[0]) == 1:
                    D1 = P
                else:
                    R = right._new_c(right._coeffs[:],skew_ring)
                    R = R // dict_right[N]
                    D = R._left_lcm_cofactor(dict_divisor[N])
                    maxtype = list(type)
                    maxtype[-1] -= 1
                    degN = N.degree()
                    cardE = cardcenter ** degN
                    maxcount = q_jordan(maxtype, cardE)
                    Q = NS // P
                    deg = P.degree()-1
                    while True:
                        while True:
                            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((deg,deg))
                            R = Q * R
                            if P.right_gcd(R).degree() == 0:
                                break
                        D1 = P.right_gcd(D*R)

                        L = left._new_c(list(left._coeffs),skew_ring)
                        L = L // D1
                        degN = N.degree()
                        for j in range(len(type)):
                            if <Py_ssize_t>(type[j]) == 1:
                                newtype = type[:-1]
                                break
                            d = L.right_gcd(NS)
                            deg = d.degree() / degN
                            if deg < <Py_ssize_t>(type[j]):
                                newtype = list(type)
                                newtype[j] = deg
                                break
                            L = L // d
                        count = q_jordan(newtype, cardE)
                        if ZZ.random_element(maxcount) < count:
                            break
                    dict_type[N] = newtype

                    D2 = D._new_c(list(D._coeffs), skew_ring)
                    D2 = D2.right_monic()
                    while D2 == D1:
                        while True:
                            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((deg,deg))
                            R = Q * R
                            if P.right_gcd(R).degree() == 0:
                                break
                        D2 = P.right_gcd(D*R)
                    dict_divisor[N] = D1.left_lcm(D2)
            factors.append((D1,1))
            left = left // D1
            right = D1 * right
            dict_right[N] = right._new_c(list(right._coeffs), skew_ring)

        factors.reverse()
        return Factorization(factors, sort=False, unit=unit)

    def factor(self, uniform=False):
        r"""
        Return a factorization of this skew polynomial.

        INPUT:

        - ``uniform`` -- boolean (default: ``False``); whether the
          output irreducible divisor should be uniformly distributed
          among all possibilities

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^3 + (t^2 + 4*t + 2)*x^2 + (3*t + 3)*x + t^2 + 1
            sage: F = a.factor(); F   # random
            (x + t^2 + 4) * (x + t + 3) * (x + t)
            sage: F.value() == a
            True

        The result of the factorization is cached. Hence, if we try
        again to factor `a`, we will get the same answer::

            sage: a.factor()  # random
            (x + t^2 + 4) * (x + t + 3) * (x + t)

        However, the algorithm is probabilistic. Hence if we first
        reinitialiaze `a`, we may get a different answer::

            sage: a = x^3 + (t^2 + 4*t + 2)*x^2 + (3*t + 3)*x + t^2 + 1
            sage: F = a.factor(); F   # random
            (x + t^2 + t + 2) * (x + 2*t^2 + t + 4) * (x + t)
            sage: F.value() == a
            True

        There is a priori no guarantee on the distribution of the
        factorizations we get. Passing in the keyword ``uniform=True``
        ensures the output is uniformly distributed among all
        factorizations::

            sage: a.factor(uniform=True)   # random
            (x + t^2 + 4) * (x + t) * (x + t + 3)
            sage: a.factor(uniform=True)   # random
            (x + 2*t^2) * (x + t^2 + t + 1) * (x + t^2 + t + 2)
            sage: a.factor(uniform=True)   # random
            (x + 2*t^2 + 3*t) * (x + 4*t + 2) * (x + 2*t + 2)

        By convention, the zero skew polynomial has no factorization::

            sage: S(0).factor()
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined
        """
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        if uniform:
            F = self._factor_uniform_c()
            if self._factorization is None:
                self._factorization = F
        else:
            if self._factorization is None:
                self._factorization = self._factor_c()
            F = self._factorization
        return F

    def count_factorizations(self):
        r"""
        Return the number of factorizations (as a product of a
        unit and a product of irreducible monic factors) of this
        skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + (4*t + 3)*x^3 + t^2*x^2 + (4*t^2 + 3*t)*x + 3*t
            sage: a.count_factorizations()
            216

        We illustrate that an irreducible polynomial in the center have
        in general a lot of distinct factorizations in the skew polynomial
        ring::

            sage: Z.<x3> = S.center()
            sage: N = x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3
            sage: N.is_irreducible()
            True
            sage: S(N).count_factorizations()
            30537115626
        """
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        cardcenter = self._parent._working_center.base_ring().cardinality()
        gencenter = self._parent._working_center.gen()
        F = self._reduced_norm_factored()
        summ = 0
        count = 1
        for N, m in F:
            summ += m
            if m == 1:
                continue
            if N != gencenter:
                count *= q_jordan(self.type(N), cardcenter**N.degree())
            count /= factorial(m)
        return count * factorial(summ)

    # Not optimized:
    # many calls to reduced_norm, reduced_norm_factor, _rdivisor_c, which are slow

    def factorizations(self):
        r"""
        Return an iterator over all factorizations (as a product
        of a unit and a product of irreducible monic factors) of
        this skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^3 + (t^2 + 1)*x^2 + (2*t + 3)*x + t^2 + t + 2
            sage: iter = a.factorizations(); iter
            <...generator object at 0x...>
            sage: next(iter)   # random
            (x + 3*t^2 + 4*t) * (x + 2*t^2) * (x + 4*t^2 + 4*t + 2)
            sage: next(iter)   # random
            (x + 3*t^2 + 4*t) * (x + 3*t^2 + 2*t + 2) * (x + 4*t^2 + t + 2)

        We can use this function to build the list of factorizations
        of `a`::

            sage: factorizations = [ F for F in a.factorizations() ]

        We do some checks::

            sage: len(factorizations) == a.count_factorizations()
            True
            sage: len(factorizations) == Set(factorizations).cardinality()  # check no duplicates
            True
            sage: for F in factorizations:
            ....:     assert F.value() == a, "factorization has a different value"
            ....:     for d,_ in F:
            ....:         assert d.is_irreducible(), "a factor is not irreducible"

        Note that the algorithm used in this method is probabilistic.
        As a consequence, if we call it two times with the same input,
        we can get different orderings::

            sage: factorizations2 = [ F for F in a.factorizations() ]
            sage: factorizations == factorizations2  # random
            False
            sage: sorted(factorizations) == sorted(factorizations2)
            True
        """
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")

        def factorizations_rec(P):
            if P.is_irreducible():
                yield [(P, 1)]
            else:
                for div in P._irreducible_divisors(True):
                    Q = P // div
                    # Here, we should update Q._norm, Q._norm_factor
                    for factors in factorizations_rec(Q):
                        factors.append((div, 1))
                        yield factors
        unit = self.leading_coefficient()
        for factors in factorizations_rec(~unit * self):
            yield Factorization(factors, sort=False, unit=unit)
