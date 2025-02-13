r"""
Witt vector rings

This module provides the class :class:`WittVectorRing` of rings of truncated
Witt vectors.

AUTHORS:

- Jacob Dennerlein (2022-11-28): initial version
- Rubén Muñoz--Bertrand (2025-02-13): major refactoring and clean-up

"""

# ****************************************************************************
#       Copyright (C) 2025 Rubén Muñoz--Bertrand
#                          <ruben.munoz--bertrand@univ-fcomte.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from itertools import product

from sage.categories.fields import Fields
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Zp, Zq
from sage.rings.padics.witt_vector import WittVector
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ring import CommutativeRing
from sage.sets.primes import Primes
from sage.structure.unique_representation import UniqueRepresentation


def _fast_char_p_power(x, n, p=None):
    r"""
    Raise x^n power in characteristic p.

    If x is not an element of a ring of characteristic p, this throws an error.

    If x is an element of GF(p^k), this is already fast.
    However, is x is a polynomial, this seems to be slow?

    EXAMPLES::

        sage: t = GF(1913)(33)
        sage: sage.rings.padics.witt_vector_ring._fast_char_p_power(t, 77)
        1371
    """
    if n not in ZZ:
        raise ValueError(f'Exponent {n} is not an integer')
    if n == 0 or x == 1:
        return x.parent().one()
    if x.parent().characteristic() not in Primes():
        raise ValueError(f'{x} is not in a ring of prime characteristic')

    x_is_Polynomial = isinstance(x, Polynomial)
    x_is_MPolynomial = isinstance(x, MPolynomial)

    if not (x_is_Polynomial or x_is_MPolynomial):
        return x**n
    if x.is_gen():
        return x**n
    if n < 0:
        x = x**-1  # This may throw an error.
        n = -n

    P = x.parent()
    if p is None:
        p = P.characteristic()
    base_p_digits = ZZ(n).digits(base=p)

    x_to_the_n = 1

    for p_exp, digit in enumerate(base_p_digits):
        if digit == 0:
            continue
        inner_term = x**digit
        term_dict = {}
        for e_int_or_tuple, c in inner_term.dict().items():
            power = p**p_exp
            new_c = _fast_char_p_power(c, power)
            new_e_tuple = None
            if x_is_Polynomial:  # Then the dict keys are ints
                new_e_tuple = e_int_or_tuple * power
            elif x_is_MPolynomial:  # Then the dict keys are ETuples
                new_e_tuple = e_int_or_tuple.emul(power)
            term_dict[new_e_tuple] = new_c
        term = P(term_dict)
        x_to_the_n *= term

    return x_to_the_n


class WittVectorRing(CommutativeRing, UniqueRepresentation):
    r"""
    Return the appropriate p-typical truncated Witt vector ring, depending on
    the input.

    INPUT::

    -- ``base_ring`` -- commutative ring of coefficients

    -- ``prec`` -- integer (default: ``1``), length of the truncated Witt
       vectors in the ring

    -- ``p`` -- prime (default: ``None``) number; when it is not set, it
       defaults to the characteristic of ``base_ring`` when it is prime.

    -- ``algorithm`` -- the name of the algorithm to use for the ring laws
       (default: ``None``); must be either ``None``, ``standard``,
       ``p_invertible``, ``finotti`` or ``Zq_isomorphism``. When it is not
       set, the most adequate algorithm is chosen. The ``standard`` algorithm
       is the schoolbook algorithm, the ``p_invertible`` is the same one in
       the case where ``p`` is invertible in ``base_ring`` and uses some
       optimisations in that setting, ``finotti`` is Finotti's algorithm which
       can be used when ``base_ring`` has characteristic ``p``, and
       ``Zq_isomorphism`` computes the ring laws in `\mathbb Z_q` when
       ``base_ring`` is `\mathbb F_q` for `q` a power of ``p``.

    EXAMPLES::

        sage: WittVectorRing(QQ,p=5)
        Ring of truncated 5-typical Witt vectors of length 1 over Rational
        Field
        sage: WittVectorRing(GF(3))
        Ring of truncated 3-typical Witt vectors of length 1 over Finite Field
        of size 3
        sage: WittVectorRing(GF(3)['t'])
        Ring of truncated 3-typical Witt vectors of length 1 over Univariate
        Polynomial Ring in t over Finite Field of size 3
        sage: WittVectorRing(Qp(7), prec=30, p=5)
        Ring of truncated 5-typical Witt vectors of length 30 over 7-adic
        Field with capped relative precision 20

    TESTS::

        sage: A = SymmetricGroup(3).algebra(QQ)
        sage: WittVectorRing(A)
        Traceback (most recent call last):
        ...
        TypeError: Symmetric group algebra of order 3 over Rational Field
        is not a commutative ring
        sage: WittVectorRing(QQ)
        Traceback (most recent call last):
        ...
        ValueError: Rational Field has non-prime characteristic
        and no prime was supplied

        sage: WittVectorRing(QQ, p=5, algorithm='moon')
        Traceback (most recent call last):
        ...
        ValueError: algorithm must be one of None, 'standard',
        'p_invertible', 'finotti', 'Zq_isomorphism'
    """
    Element = WittVector

    def __classcall_private__(cls, base_ring, prec=1, p=None, algorithm=None):
        r"""
        Construct the ring of truncated Witt vectors from the parameters.
        """
        if not isinstance(base_ring, CommutativeRing):
            raise TypeError(f'{base_ring} is not a commutative ring')

        if not (isinstance(prec, int) or isinstance(prec, Integer)):
            raise TypeError(f'{prec} is not an integer')
        elif prec <= 0:
            raise ValueError(f'{prec} must be positive')

        char = base_ring.characteristic()
        if p is None:
            if char not in Primes():
                raise ValueError(f'{base_ring} has non-prime characteristic '
                                 'and no prime was supplied')
            else:
                prime = char
        else:
            prime = p

        if algorithm not in [None, 'standard', 'p_invertible', 'finotti',
                             'Zq_isomorphism']:
            raise ValueError("algorithm must be one of None, 'standard', "
                             "'p_invertible', 'finotti', 'Zq_isomorphism'")

        if prime == char:
            if algorithm == 'p_invertible':
                raise ValueError("The 'p_invertible' algorithm only works "
                                 "when p is a unit in the ring of "
                                 "coefficients.")
            elif base_ring in Fields().Finite():
                if algorithm is None:
                    algorithm = 'Zq_isomorphism'
            else:
                if algorithm is None:
                    algorithm = 'finotti'
                elif algorithm == 'Zq_isomorphism':
                    raise ValueError("The 'Zq_isomorphism' algorithm only "
                                     "works when the coefficient ring is a "
                                     "finite field of characteristic p.")
        else:
            if algorithm == 'finotti':
                raise ValueError("The 'finotti' algorithm only works for "
                                 "coefficients rings of characteristic p.")
            elif algorithm == 'Zq_isomorphism':
                raise ValueError("The 'Zq_isomorphism' algorithm only works "
                                 "when the coefficient ring is a finite "
                                 "field of characteristic p.")
            elif base_ring(prime).is_unit():
                if algorithm is None:
                    algorithm = 'p_invertible'
            else:
                if algorithm is None:
                    algorithm = 'standard'
                elif algorithm == 'p_invertible':
                    raise ValueError("The 'p_invertible' algorithm only "
                                     "works when p is a unit in the ring of "
                                     "coefficients.")

        return cls.__classcall__(cls, base_ring, prec, prime, algorithm)

    def __init__(self, base_ring, prec, prime, algorithm):
        r"""
        Initialises ``self``.

        EXAMPLES::

            sage: WittVectorRing(ZZ, p=5, prec=2)
            Ring of truncated 5-typical Witt vectors of length 2 over Integer
            Ring
        """
        self.prec = prec
        self.prime = prime

        self._algorithm = algorithm
        self.sum_polynomials = None
        self.prod_polynomials = None

        if algorithm == 'standard':
            self._generate_sum_and_product_polynomials(base_ring)

        elif algorithm == 'finotti':
            self._generate_binomial_table()

        CommutativeRing.__init__(self, base_ring)

    def __iter__(self):
        """
        Iterator for truncated Witt vector rings.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), p=3, prec=2)
            sage: [w for w in W]
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1),
            (2, 2)]
        """
        for t in product(self.base(), repeat=self.prec):
            yield self(t)

    def _repr_(self):
        """
        Return a string representation of the ring.

        EXAMPLES::

            sage: WittVectorRing(QQ, p=2, prec=5)
            Ring of truncated 2-typical Witt vectors of length 5 over Rational
            Field
        """
        return f"Ring of truncated {self.prime}-typical Witt vectors of "\
               f"length {self.prec} over {self.base()}"

    def _coerce_map_from_(self, S):
        """"
        Check whether there is a coerce map from ``S``.

        EXAMPLES::

            sage: W=WittVectorRing(GF(25), p=5,prec=2)
            sage: W.has_coerce_map_from(WittVectorRing(GF(5), p=5, prec=3))
            True
        """
        if isinstance(S, WittVectorRing):
            return (
                S.precision() >= self.precision()
                and self.base().has_coerce_map_from(S.base()))
        return S is ZZ

    def _generate_sum_and_product_polynomials(self, base):
        """
        Generates the sum and product polynomials defining the ring laws of
        truncated Witt vectors for the ``standard`` algorithm.

        EXAMPLES::

            sage: P.<X1,X2,Y1,Y2>=PolynomialRing(GF(3),'X1,X2,Y1,Y2')
            sage: W=WittVectorRing(P, p=3, prec=2)
            sage: W([X1,X2])+W([Y1,Y2])
            (X1 + Y1, -X1^2*Y1 - X1*Y1^2 + X2 + Y2)
            sage: W([X1,X2])*W([Y1,Y2])
            (X1*Y1, X2*Y1^3 + X1^3*Y2)
        """
        p = self.prime
        prec = self.prec
        x_var_names = [f'X{i}' for i in range(prec)]
        y_var_names = [f'Y{i}' for i in range(prec)]
        var_names = x_var_names + y_var_names

        # Okay, what's going on here? Sage, by default, relies on
        # Singular for Multivariate Polynomial Rings, but Singular uses
        # only SIXTEEN bits (unsigned) to store its exponents. So if we
        # want exponents larger than 2^16 - 1, we have to use the
        # generic implementation. However, after some experimentation,
        # it seems like the generic implementation is faster?
        #
        # After trying to compute S_4 for p=5, it looks like generic is
        # faster for  very small polys, and MUCH slower for large polys.
        # So we'll default to singular unless we can't use it.
        #
        # Remark: Since when is SIXTEEN bits sufficient for anyone???
        #
        if p**(prec-1) >= 2**16:
            implementation = 'generic'
        else:
            implementation = 'singular'

        # We first generate the "universal" polynomials and then project
        # to the base ring.
        R = PolynomialRing(ZZ, var_names, implementation=implementation)
        x_y_vars = R.gens()
        x_vars = x_y_vars[:prec]
        y_vars = x_y_vars[prec:]

        self.sum_polynomials = [0]*(self.prec)
        for n in range(prec):
            s_n = x_vars[n] + y_vars[n]
            for i in range(n):
                s_n += ((x_vars[i]**(p**(n-i)) + y_vars[i]**(p**(n-i))
                        - self.sum_polynomials[i]**(p**(n-i))) / p**(n-i))
            self.sum_polynomials[n] = R(s_n)

        self.prod_polynomials = [x_vars[0] * y_vars[0]] + [0]*(self.prec)
        for n in range(1, prec):
            x_poly = sum([p**i * x_vars[i]**(p**(n-i)) for i in range(n+1)])
            y_poly = sum([p**i * y_vars[i]**(p**(n-i)) for i in range(n+1)])
            p_poly = sum([p**i * self.prod_polynomials[i]**(p**(n-i))
                         for i in range(n)])
            p_n = (x_poly*y_poly - p_poly) // p**n
            self.prod_polynomials[n] = p_n

        # We have to use generic here, because Singular doesn't support
        # Polynomial Rings over Polynomial Rings. For example,
        # ``PolynomialRing(GF(5)['x'], ['X', 'Y'],
        #                  implementation='singular')``
        # will fail.
        S = PolynomialRing(base, x_y_vars, implementation='generic')
        for n in range(prec):
            self.sum_polynomials[n] = S(self.sum_polynomials[n])
            self.prod_polynomials[n] = S(self.prod_polynomials[n])

    def _generate_binomial_table(self):
        """
        Generates a table of binomial coefficients.

        TODO::

            Suprisingly, it seems that no function in SageMath does this;
            ``arith.misc.binomial`` and ``combinat.sloane_functions.A000984``
            only compute one value at a time, so it is not efficient. Ideally,
            this should be implemented elsewhere and in a low-level language.
        """
        import numpy as np
        p = self.prime
        R = Zp(p, prec=self.prec+1, type='fixed-mod')
        v_p = ZZ.valuation(p)
        table = [[0]]
        for k in range(1, self.prec+1):
            row = np.empty(p**k, dtype=int)
            row[0] = 0
            prev_bin = 1
            for i in range(1, p**k // 2 + 1):
                val = v_p(i)
                # Instead of calling binomial each time, we compute the
                # coefficients recursively. This is MUCH faster.
                next_bin = prev_bin * (p**k - (i-1)) // i
                prev_bin = next_bin
                series = R(-next_bin // p**(k-val))
                for _ in range(val):
                    temp = series % p
                    series = (series - R.teichmuller(temp)) // p
                row[i] = ZZ(series % p)
                row[p**k - i] = row[i]  # binomial coefficients are symmetric
            table.append(row)
        self.binomial_table = table

    def _eta_bar(self, vec, eta_index):
        r"""
        Generate the `\eta_i` for ``finotti``'s algorithm.
        """
        vec = tuple(x for x in vec if x != 0)  # strip zeroes

        # special cases
        if len(vec) <= 1:
            return 0
        if eta_index == 0:
            return sum(vec)

        # renaming to match notation in paper
        k = eta_index
        p = self.prime
        # if vec = (x,y), we know what to do: Theorem 8.6
        if len(vec) == 2:
            # Here we have to check if we've pre-computed already
            x, y = vec
            scriptN = [[None] for _ in range(k+1)]  # each list starts with
            # None, so that indexing matches paper

            # calculate first N_t scriptN's
            for t in range(1, k+1):
                for i in range(1, p**t):
                    scriptN[t].append(self.binomial_table[t][i]
                                      * _fast_char_p_power(x, i)
                                      * _fast_char_p_power(y, p**t - i))
            indexN = [p**i - 1 for i in range(k+1)]
            for t in range(2, k+1):
                for i in range(1, t):
                    # append scriptN_{t, N_t+l}
                    next_scriptN = self._eta_bar(
                        scriptN[t-i][1:indexN[t-i]+t-i], i
                    )
                    scriptN[t].append(next_scriptN)
            return sum(scriptN[k][1:])

        # if vec is longer, we split and recurse: Proposition 5.4
        # This is where we need to using multiprocessing.
        else:
            m = len(vec) // 2
            v_1 = vec[:m]
            v_2 = vec[m:]
            s_1 = sum(v_1)
            s_2 = sum(v_2)
            scriptM = [[] for _ in range(k+1)]
            for t in range(1, k+1):
                scriptM[t].append(self._eta_bar(v_1, t))
                scriptM[t].append(self._eta_bar(v_2, t))
                scriptM[t].append(self._eta_bar((s_1, s_2), t))
            for t in range(2, k+1):
                for s in range(1, t):
                    result = self._eta_bar(scriptM[t-s], s)
                    scriptM[t].append(result)
            return sum(scriptM[k])

    def _series_to_vector(self, series):
        r"""
        Computes the canonical bijection from `\mathbb Z_q` to
        `W(\mathbb F_q)`.
        """
        F = self.base()  # known to be finite
        R = Zq(F.cardinality(), prec=self.prec, type='fixed-mod',
               modulus=F.polynomial(), names=['z'])
        K = R.residue_field()
        p = self.prime

        series = R(series)
        witt_vector = []
        for i in range(self.prec):
            temp = K(series)
            elem = temp.polynomial()(F.gen())  # hack to convert to F
            # (K != F for some reason)
            witt_vector.append(elem**(p**i))
            series = (series - R.teichmuller(temp)) // p
        return witt_vector

    def _vector_to_series(self, vec):
        r"""
        Computes the canonical bijection from `W(\mathbb F_q)` to
        `\mathbb Z_q`.
        """
        F = self.base()
        R = Zq(F.cardinality(), prec=self.prec, type='fixed-mod',
               modulus=F.polynomial(), names=['z'])
        K = R.residue_field()
        p = self.prime

        series = R.zero()
        for i in range(self.prec):
            temp = vec[i].nth_root(p**i)
            elem = temp.polynomial()(K.gen())
            # hack to convert to K (F != K for some reason)

            series += p**i * R.teichmuller(elem)
        return series

    def characteristic(self):
        """
        Return the characteristic of ``self``.

        EXAMPLES::

            WittVectorRing(GF(25), p=5, prec=3).characteristic()
            125
            sage: WittVectorRing(ZZ, p=2, prec=4).characteristic()
            0
            sage: WittVectorRing(Integers(18), p=3, prec=3).characteristic()
            162
        """
        p = self.prime
        if self.base()(p).is_unit():
            # If p is invertible, W_n(R) is isomorphic to R^n.
            return self.base().characteristic()

        # This is Jacob Dennerlein's Corollary 3.3. in "Computational
        # Aspects of Mixed Characteristic Witt Vectors" (preprint)
        return p**(self.prec-1) * self.base().characteristic()

    def precision(self):
        """
        Return the length of the truncated Witt vectors in ``length``.

        EXAMPLES::
            sage: WittVectorRing(GF(9), p=3, prec=3).precision()
            3
        """
        return self.prec

    def random_element(self):
        """
        Return the length of the truncated Witt vectors in ``length``.

        EXAMPLES::
            sage: WittVectorRing(GF(27), prec=2).random_element()  # random
            (z3, 2*z3^2 + 1)
        """
        return self.element_class(self, tuple(self.base().random_element()
                                              for _ in range(self.prec)))

    def teichmuller_lift(self, x):
        """
        Return the Teichmüller lift of ``x`` in ``self``. This lift is
        sometimes known as the multiplicative lift of ``x``, in order to avoid
        refering to a nazi mathematician.

        EXAMPLES::
            sage: WittVectorRing(GF(125), prec=2).teichmuller_lift(3)
            (3, 0)
        """
        if x not in self.base():
            raise Exception(f'{x} not in {self.base()}')
        return self.element_class(self,
                                  (x,) + tuple(0 for _ in range(self.prec-1)))

    def is_finite(self):
        """
        Return whether ``self`` is a finite ring.

        EXAMPLES::
            sage: WittVectorRing(GF(23)).is_finite()
            True
            sage: WittVectorRing(ZZ, p=2).is_finite()
            False
        """
        return self.base().is_finite()

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::
            sage: WittVectorRing(GF(17), prec=2).cardinality()
            289
            sage: WittVectorRing(QQ, p=2).cardinality()
            +Infinity
        """
        return self.base().cardinality()**(self.prec)
