# sage.doctest: optional - sage.rings.finite_rings
r"""
Drinfeld modules over rings of characteristic zero

This module provides the classes
:class:`sage.rings.function_fields.drinfeld_module.charzero_drinfeld_module.DrinfeldModule_charzero`,
:class:`sage.rings.function_fields.drinfeld_module.charzero_drinfeld_module.DrinfeldModule_rational`,
which inherits
:class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

AUTHORS:

- David Ayotte (2023-09)
- Xavier Caruso (2024-12) - computation of class polynomials and Taelman's units
"""

# *****************************************************************************
#        Copyright (C) 2022 David Ayotte <david.ayotte@mail.concordia.ca>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from .drinfeld_module import DrinfeldModule

from sage.rings.integer_ring import ZZ

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')


class DrinfeldModule_charzero(DrinfeldModule):
    r"""
    This class implements Drinfeld `\mathbb{F}_q[T]`-modules defined
    over fields of `\mathbb{F}_q[T]`-characteristic zero.

    Recall that the `\mathbb{F}_q[T]`-*characteristic* is defined as the
    kernel of the underlying structure morphism. For general definitions
    and help on Drinfeld modules, see class
    :class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

    .. RUBRIC:: Construction:

    The user does not ever need to directly call
    ``DrinfeldModule_charzero`` --- the metaclass ``DrinfeldModule`` is
    responsible for instantiating the right class depending on the
    input::

        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: phi = DrinfeldModule(A, [T, 1])
        sage: phi
        Drinfeld module defined by T |--> t + T

    ::

        sage: isinstance(phi, DrinfeldModule)
        True
        sage: from sage.rings.function_field.drinfeld_modules.charzero_drinfeld_module import DrinfeldModule_charzero
        sage: isinstance(phi, DrinfeldModule_charzero)
        True

    .. RUBRIC:: Logarithm and exponential

    It is possible to calculate the logarithm and the exponential of
    any Drinfeld modules of characteristic zero::

        sage: A = GF(2)['T']
        sage: K.<T> = Frac(A)
        sage: phi = DrinfeldModule(A, [T, 1])
        sage: phi.exponential()
        z + ((1/(T^2+T))*z^2) + ((1/(T^8+T^6+T^5+T^3))*z^4) + O(z^8)
        sage: phi.logarithm()
        z + ((1/(T^2+T))*z^2) + ((1/(T^6+T^5+T^3+T^2))*z^4) + O(z^8)

    .. RUBRIC:: Goss polynomials

    Goss polynomials are a sequence of polynomials related with the
    analytic theory of Drinfeld module. They provide a function field
    analogue of certain classical trigonometric functions::

        sage: A = GF(2)['T']
        sage: K.<T> = Frac(A)
        sage: phi = DrinfeldModule(A, [T, 1])
        sage: phi.goss_polynomial(1)
        X
        sage: phi.goss_polynomial(2)
        X^2
        sage: phi.goss_polynomial(3)
        X^3 + (1/(T^2 + T))*X^2

    .. RUBRIC:: Base fields of `\mathbb{F}_q[T]`-characteristic zero

    The base fields need not only be fraction fields of polynomials
    ring. In the following example, we construct a Drinfeld module over
    `\mathbb{F}_q((1/T))`, the completion of the rational function field
    at the place `1/T`::

        sage: A.<T> = GF(2)[]
        sage: L.<s> = LaurentSeriesRing(GF(2))  # s = 1/T
        sage: phi = DrinfeldModule(A, [1/s, s + s^2 + s^5 + O(s^6), 1+1/s])
        sage: phi(T)
        (s^-1 + 1)*t^2 + (s + s^2 + s^5 + O(s^6))*t + s^-1

    One can also construct Drinfeld modules over SageMath's global
    function fields::

        sage: A.<T> = GF(5)[]
        sage: K.<z> = FunctionField(GF(5))  # z = T
        sage: phi = DrinfeldModule(A, [z, 1, z^2])
        sage: phi(T)
        z^2*t^2 + t + z
    """
    @cached_method
    def _compute_coefficient_exp(self, k):
        r"""
        Return the `q^k`-th coefficient of the exponential of this
        Drinfeld module.

        INPUT:

        - ``k`` -- integer; the index of the coefficient

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: phi._compute_coefficient_exp(0)
            1
            sage: phi._compute_coefficient_exp(1)
            1/(T^2 + T)
            sage: phi._compute_coefficient_exp(2)
            1/(T^8 + T^6 + T^5 + T^3)
            sage: phi._compute_coefficient_exp(3)
            1/(T^24 + T^20 + T^18 + T^17 + T^14 + T^13 + T^11 + T^7)
        """
        k = ZZ(k)
        if k.is_zero():
            return self._base.one()
        q = self._Fq.cardinality()
        c = self._base.zero()
        for i in range(k):
            j = k - i
            c += self._compute_coefficient_exp(i)*self._compute_coefficient_log(j)**(q**i)
        return -c

    def exponential(self, name='z'):
        r"""
        Return the exponential of this Drinfeld module.

        Note that the exponential is only defined when the
        `\mathbb{F}_q[T]`-characteristic is zero.

        INPUT:

        - ``name`` -- string (default: ``'z'``); the name of the
          generator of the lazy power series ring

        OUTPUT: a lazy power series over the base field

        EXAMPLES::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: exp = phi.exponential(); exp
            z + ((1/(T^2+T))*z^2) + ((1/(T^8+T^6+T^5+T^3))*z^4) + O(z^8)

        The exponential is returned as a lazy power series, meaning that
        any of its coefficients can be computed on demands::

            sage: exp[2^4]
            1/(T^64 + T^56 + T^52 + ... + T^27 + T^23 + T^15)
            sage: exp[2^5]
            1/(T^160 + T^144 + T^136 + ... + T^55 + T^47 + T^31)

        Example in higher rank::

            sage: A = GF(5)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T^2, T + T^2 + T^4, 1])
            sage: exp = phi.exponential(); exp
            z + ((T/(T^4+4))*z^5) + O(z^8)

        The exponential is the compositional inverse of the logarithm
        (see :meth:`logarithm`)::

            sage: log = phi.logarithm(); log
            z + ((4*T/(T^4+4))*z^5) + O(z^8)
            sage: exp.compose(log)
            z + O(z^8)
            sage: log.compose(exp)
            z + O(z^8)

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: exp = phi.exponential()
            sage: exp[2] == 1/(T**q - T)  # expected value
            True
            sage: exp[2^2] == 1/((T**(q**2) - T)*(T**q - T)**q)  # expected value
            True
            sage: exp[2^3] == 1/((T**(q**3) - T)*(T**(q**2) - T)**q*(T**q - T)**(q**2))  # expected value
            True

        REFERENCE:

        See section 4.6 of [Gos1998]_ for the definition of the
        exponential.
        """
        L = LazyPowerSeriesRing(self._base, name)
        zero = self._base.zero()
        q = self._Fq.cardinality()

        def coeff_exp(k):
            # Return the k-th coefficient of the exponential.
            k = ZZ(k)
            v, u = k.val_unit(q)
            if u == 1:
                return self._compute_coefficient_exp(v)
            else:
                return zero
        return L(coeff_exp, valuation=1)

    @cached_method
    def _compute_coefficient_log(self, k):
        r"""
        Return the `q^k`-th coefficient of the logarithm of this
        Drinfeld module.

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: phi._compute_coefficient_log(0)
            1
            sage: phi._compute_coefficient_log(1)
            1/(T^2 + T)
            sage: phi._compute_coefficient_log(2)
            1/(T^6 + T^5 + T^3 + T^2)
            sage: phi._compute_coefficient_log(3)
            1/(T^14 + T^13 + T^11 + T^10 + T^7 + T^6 + T^4 + T^3)
        """
        k = ZZ(k)
        if k.is_zero():
            return self._base.one()
        r = self._gen.degree()
        T = self._gen[0]
        q = self._Fq.cardinality()
        c = self._base.zero()
        for i in range(k):
            j = k - i
            if j < r + 1:
                c += self._compute_coefficient_log(i)*self._gen[j]**(q**i)
        return c/(T - T**(q**k))

    def logarithm(self, name='z'):
        r"""
        Return the logarithm of the given Drinfeld module.

        By definition, the logarithm is the compositional inverse of the
        exponential (see :meth:`exponential`). Note that the logarithm
        is only defined when the `\mathbb{F}_q[T]`-characteristic is
        zero.

        INPUT:

        - ``name`` -- string (default: ``'z'``); the name of the
          generator of the lazy power series ring

        OUTPUT: a lazy power series over the base field

        EXAMPLES::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: log = phi.logarithm(); log
            z + ((1/(T^2+T))*z^2) + ((1/(T^6+T^5+T^3+T^2))*z^4) + O(z^8)

        The logarithm is returned as a lazy power series, meaning that
        any of its coefficients can be computed on demands::

            sage: log[2^4]
            1/(T^30 + T^29 + T^27 + ... + T^7 + T^5 + T^4)
            sage: log[2^5]
            1/(T^62 + T^61 + T^59 + ... + T^8 + T^6 + T^5)

        Example in higher rank::

            sage: A = GF(5)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T^2, T + T^2 + T^4, 1])
            sage: phi.logarithm()
            z + ((4*T/(T^4+4))*z^5) + O(z^8)

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = 2
            sage: log[2] == -1/((T**q - T))  # expected value
            True
            sage: log[2**2] == 1/((T**q - T)*(T**(q**2) - T))  # expected value
            True
            sage: log[2**3] == -1/((T**q - T)*(T**(q**2) - T)*(T**(q**3) - T))  # expected value
            True
        """
        L = LazyPowerSeriesRing(self._base, name)
        q = self._Fq.cardinality()

        def coeff_log(k):
            # Return the k-th coefficient of the logarithm
            k = ZZ(k)
            v, u = k.val_unit(q)
            if u == 1:
                return self._compute_coefficient_log(v)
            else:
                return self._base.zero()
        return L(coeff_log, valuation=1)

    @cached_method
    def _compute_goss_polynomial(self, n, q, poly_ring, X):
        r"""
        Utility function for computing the n-th Goss polynomial.

        The user should not call this method directly, but
        :meth:`goss_polynomial` instead.

        TESTS::

            sage: A = GF(2^2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T+1, T^2, 1])
            sage: poly_ring = phi.base()['X']
            sage: X = poly_ring.gen()
            sage: phi._compute_goss_polynomial(0, 2^2, poly_ring, X)
            0
            sage: phi._compute_goss_polynomial(3, 2^2, poly_ring, X)
            X^3
            sage: phi._compute_goss_polynomial(4*3, 2^2, poly_ring, X)
            X^12
            sage: phi._compute_goss_polynomial(9, 2^2, poly_ring, X)
            X^9 + (1/(T^3 + T^2 + T))*X^6 + (1/(T^6 + T^4 + T^2))*X^3
        """
        # Trivial cases
        if n.is_zero():
            return poly_ring.zero()
        if n <= q - 1:
            return X**n
        if n % q == 0:
            return self.goss_polynomial(n // q)**q
        # General case
        pol = poly_ring.zero()
        m = q
        i = 1
        while m < n:
            pol += self._compute_coefficient_exp(i) * self._compute_goss_polynomial(n - m, q, poly_ring, X)
            m *= q
            i += 1
        return X*(self._compute_goss_polynomial(n - 1, q, poly_ring, X) + pol)

    def goss_polynomial(self, n, var='X'):
        r"""
        Return the `n`-th Goss polynomial of the Drinfeld module.

        Note that Goss polynomials are only defined for Drinfeld modules
        of characteristic zero.

        INPUT:

        - ``n`` -- integer; the index of the Goss polynomial

        - ``var``-- string (default: ``'X'``); the name of polynomial
          variable

        OUTPUT: a univariate polynomial in ``var`` over the base `A`-field

        EXAMPLES::

            sage: A = GF(3)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])  # Carlitz module
            sage: phi.goss_polynomial(1)
            X
            sage: phi.goss_polynomial(2)
            X^2
            sage: phi.goss_polynomial(4)
            X^4 + (1/(T^3 + 2*T))*X^2
            sage: phi.goss_polynomial(5)
            X^5 + (2/(T^3 + 2*T))*X^3
            sage: phi.goss_polynomial(10)
            X^10 + (1/(T^3 + 2*T))*X^8 + (1/(T^6 + T^4 + T^2))*X^6 + (1/(T^9 + 2*T^3))*X^4 + (1/(T^18 + 2*T^12 + 2*T^10 + T^4))*X^2

        REFERENCE:

        Section 3 of [Gek1988]_ provides an exposition of Goss
        polynomials.
        """
        n = ZZ(n)
        K = self.base()
        poly_ring = K[var]
        X = poly_ring.gen()
        q = self._Fq.cardinality()
        return self._compute_goss_polynomial(n, q, poly_ring, X)


class DrinfeldModule_rational(DrinfeldModule_charzero):
    """
    A class for Drinfeld modules defined over the fraction
    field of the underlying function field
    """
    def _phiT_matrix(self, polynomial_part):
        r"""
        Return the matrix giving the action of `\phi_T` modulo `u^s`
        where `u = 1/T` is the uniformizer at infinity `s` is chosen
        such that `u^s` is in the domain of convergence of the logarithm.

        It is an helper function; do not call it directly.

        INPUT:

        - ``polynomial_part`` -- boolean; if ``False``, omit the
          part with negative powers of `u`; if ``True``, return this
          part as a polynomial vector in `T`

        TESTS::

            sage: q = 5
            sage: Fq = GF(q)
            sage: A = Fq['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T^20])
            sage: phi._phiT_matrix(False)
            [0 0 0 0 0]
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 1]
            sage: phi._phiT_matrix(True)
            (
            [0 0 0 0 0]
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 1], (T^15 + 1, T^10, T^5, 1, 0)
            )
        ::

            sage: psi = DrinfeldModule(A, [T, 1/T])
            sage: psi._phiT_matrix(False)
            Traceback (most recent call last):
            ...
            ValueError: the Drinfeld module must have polynomial coefficients
        """
        A = self.function_ring()
        Fq = A.base_ring()
        q = Fq.cardinality()
        r = self.rank()

        gs = []
        for g in self.coefficients(sparse=False):
            g = g.backend(force=True)
            if g.denominator().is_one():
                gs.append(A(g.numerator().list()))
            else:
                raise ValueError("the Drinfeld module must have polynomial coefficients")
        s = max(gs[i].degree() // (q**i - 1) for i in range(1, r+1))

        M = matrix(Fq, s)
        if polynomial_part:
            P = vector(A, s)
        qk = 1
        for k in range(r+1):
            for i in range(s):
                e = (i+1)*qk
                if polynomial_part:
                    P[i] += gs[k] >> e
                for j in range(s):
                    e -= 1
                    if e < 0:
                        break
                    M[i, j] += gs[k][e]
            qk *= q

        if polynomial_part:
            return M, P
        else:
            return M

    def class_polynomial(self):
        r"""
        Return the class polynomial, that is the Fitting ideal
        of the class module, of this Drinfeld module.

        EXAMPLES:

        We check that the class module of the Carlitz module
        is trivial::

            sage: q = 5
            sage: Fq = GF(q)
            sage: A = Fq['T']
            sage: K.<T> = Frac(A)
            sage: C = DrinfeldModule(A, [T, 1]); C
            Drinfeld module defined by T |--> t + T
            sage: C.class_polynomial()
            1

        When the coefficients of the Drinfeld module have small
        enough degrees, the class module is always trivial::

            sage: r = 4
            sage: phi = DrinfeldModule(A, [T] + [A.random_element(degree=q**i) for i in range(1, r+1)])
            sage: phi.class_polynomial()
            1

        Here is an example with a nontrivial class module::

            sage: phi = DrinfeldModule(A, [T, -T^(2*q-1) + 2*T^(q-1)])
            sage: phi.class_polynomial()
            T + 3

        TESTS:

        The Drinfeld module must have polynomial coefficients::

            sage: phi = DrinfeldModule(A, [T, 1/T])
            sage: phi.class_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: the Drinfeld module must have polynomial coefficients
        """
        A = self.function_ring()
        Fq = A.base_ring()
        M = self._phiT_matrix(False)
        s = M.nrows()
        if s == 0:
            # small case
            return A.one()

        v = vector(Fq, s)
        v[s-1] = 1
        vs = [v]
        for i in range(s-1):
            v = v*M
            vs.append(v)
        V = matrix(vs)
        V.echelonize()

        dim = V.rank()
        pivots = V.pivots()
        j = ip = 0
        for i in range(dim, s):
            while ip < dim and j == pivots[ip]:
                j += 1
                ip += 1
            V[i,j] = 1

        N = (V * M * ~V).submatrix(dim, dim)
        return A(N.charpoly())

    def taelman_exponential_unit(self):
        r"""
        Return the exponential of a fundamental Taelman's unit
        of this Drinfeld module.

        A Taelman's unit is by definition an element `x \in
        \FF_q((1/T))` whose exponential falls in `\FF_q[T]`.

        Taelman's units form a `\FF_q[T]`-line in `\FF_q((1/T))`;
        a fundamental unit is by definition a generator of this line.

        We refer to [Tae2012]_ for more details about this construction.

        EXAMPLES:

        The Taelman exponential unit of The Carlitz module is `1`::

            sage: q = 7
            sage: Fq = GF(q)
            sage: A = Fq['T']
            sage: K.<T> = Frac(A)
            sage: C = DrinfeldModule(A, [T, 1]); C
            Drinfeld module defined by T |--> t + T
            sage: C.taelman_exponential_unit()
            1

        The same occurs more generally when the coefficients of the
        Drinfeld module have small enough degrees::

            sage: r = 4
            sage: phi = DrinfeldModule(A, [T] + [A.random_element(degree=q**i) for i in range(1, r+1)])
            sage: phi.taelman_exponential_unit()
            1

        Usually, as soon as we leave the world of small Drinfeld modules,
        Taelman's exponential units are highly non trivial::

            sage: phi = DrinfeldModule(A, [T, T^(2*q+1), T^3])
            sage: phi.taelman_exponential_unit()
            T^52 + T^22 + T^8 + T^2 + 1
        """
        A = self.function_ring()
        Fq = A.base_ring()
        q = Fq.cardinality()
        M, P = self._phiT_matrix(True)
        s = M.nrows()
        if s == 0:
            # small case
            return self.base().one()

        gs = self.coefficients(sparse=False)
        v = vector(Fq, s)
        v[s-1] = 1
        p = A.zero()
        vs = [v]
        ps = [p]
        for i in range(s):
            pq = p
            p = v * P
            for j in range(len(gs) - 1):
                p += gs[j] * pq
                pq = pq ** q
            p += gs[-1] * pq
            v = v * M
            vs.append(v)
            ps.append(p)
        vs.reverse()
        ps.reverse()
        V = matrix(vs)

        unit = V.left_kernel().basis()[0]
        expunit = sum(unit[i]*ps[i] for i in range(s+1))
        expunit /= expunit.numerator().leading_coefficient()
        return expunit
