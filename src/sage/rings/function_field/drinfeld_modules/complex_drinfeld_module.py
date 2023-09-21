from .drinfeld_module import DrinfeldModule

from sage.rings.integer_ring import ZZ

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')

class DrinfeldModule_complex(DrinfeldModule):
    @cached_method
    def _compute_coefficient_exp(self, k):
        r"""
        Return the `q^k`-th coefficient of the exponential of this Drinfeld module.

        INPUT:

        - ``k`` (integer) -- the index of the coefficient

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

        - ``name`` (string, default: ``'z'``) -- the name of the
          generator of the lazy power series ring.

        OUTPUT:

        A lazy power series over the base field.

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
        if self.category()._characteristic:
            raise ValueError(f"characteristic must be zero (={self.characteristic()})")
        L = LazyPowerSeriesRing(self._base, name)
        zero = self._base.zero()
        q = self._Fq.cardinality()

        def coeff_exp(k):
            # Return the k-th coefficient of the exponential.
            k = ZZ(k)
            if k.is_power_of(q):
                return self._compute_coefficient_exp(k.log(q))
            else:
                return zero
        return L(coeff_exp, valuation=1)

    @cached_method
    def _compute_coefficient_log(self, k):
        r"""
        Return the `q^k`-th coefficient of the logarithm of this Drinfeld module.

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

        - ``name`` (string, default: ``'z'``) -- the name of the
          generator of the lazy power series ring.

        OUTPUT:

        A lazy power series over the base field.

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
        if self.category()._characteristic:
            raise ValueError(f"characteristic must be zero (={self.characteristic()})")
        L = LazyPowerSeriesRing(self._base, name)
        zero = self._base.zero()
        q = self._Fq.cardinality()

        def coeff_log(k):
            # Return the k-th coefficient of the logarithm
            k = ZZ(k)
            if k.is_power_of(q):
                return self._compute_coefficient_log(k.log(q))
            else:
                return self._base.zero()
        return L(coeff_log, valuation=1)
