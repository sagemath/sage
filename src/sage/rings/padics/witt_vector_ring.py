r"""
Witt vector rings

This module provides the class :class:`WittVectorRing` of rings of truncated
Witt vectors.

AUTHORS:

- Jacob Dennerlein (2022-11-28): initial version
- Rubén Muñoz-\-Bertrand (2025-02-13, 2025-12-23): major refactoring and
  new features
"""

# ****************************************************************************
#       Copyright (C) 2025 Rubén Muñoz--Bertrand
#                          <ruben.munoz-bertrand@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from itertools import product
from math import prod

from collections.abc import Iterator

from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.categories.integral_domains import IntegralDomains
from sage.combinat.tuple import Tuples
from sage.misc.latex import latex
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.morphism import RingHomomorphism, RingMap
from sage.rings.padics.factory import Zp
from sage.rings.padics.witt_vector import (
    WittVector_finotti,
    WittVector_phantom,
    WittVector_pinvertible,
    WittVector_standard,
)
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.integer_range import IntegerRange
from sage.sets.primes import Primes
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


def fast_char_p_power(x, n, p=None):
    r"""
    Return `x^n` assuming that `x` lives in a ring of characteristic `p`.

    If `x` is not an element of a ring of characteristic `p`,
    this throws an error.

    EXAMPLES::

        sage: from sage.rings.padics.witt_vector_ring import fast_char_p_power
        sage: t = GF(1913)(33)
        sage: fast_char_p_power(t, 77)
        1371

    ::

        sage: K.<t> = GF(5^3)
        sage: fast_char_p_power(t, 385)
        4*t^2 + 1
        sage: t^385
        4*t^2 + 1

    ::

        sage: A.<x> = K[]
        sage: fast_char_p_power(x + 1, 10)
        x^10 + 2*x^5 + 1

    ::

        sage: B.<u,v> = K[]
        sage: fast_char_p_power(u + v, 1250)
        u^1250 + 2*u^625*v^625 + v^1250
    """
    x_is_Polynomial = isinstance(x, Polynomial)
    x_is_MPolynomial = isinstance(x, MPolynomial)

    if not (x_is_Polynomial or x_is_MPolynomial):
        return x**n
    if x.is_gen():
        return x**n
    if n < 0:
        x = ~x
        n = -n

    P = x.parent()
    if p is None:
        p = P.characteristic()
    base_p_digits = ZZ(n).digits(base=p)

    xn = 1

    for p_exp, digit in enumerate(base_p_digits):
        if digit == 0:
            continue
        inner_term = x**digit
        term_dict = {}
        for e_int_or_tuple, c in inner_term.dict().items():
            power = p**p_exp
            new_c = fast_char_p_power(c, power)
            new_e_tuple = None
            if x_is_Polynomial:  # Then the dict keys are ints
                new_e_tuple = e_int_or_tuple * power
            elif x_is_MPolynomial:  # Then the dict keys are ETuples
                new_e_tuple = e_int_or_tuple.emul(power)
            term_dict[new_e_tuple] = new_c
        term = P(term_dict)
        xn *= term

    return xn


class WittVectorRing(Parent, UniqueRepresentation):
    r"""
    Return the appropriate `p`-typical truncated Witt vector ring.

    INPUT:

    - ``coefficient_ring`` -- commutative ring of coefficients

    - ``prec`` -- integer (default: `1`), length of the truncated Witt
      vectors in the ring

    - ``p`` -- a prime number (default: ``None``); when it is not set, it
      defaults to the characteristic of ``coefficient_ring`` when it is prime.

    - ``algorithm`` -- the name of the algorithm to use for the ring laws
      (default: ``None``); when it is not set, the most adequate algorithm
      is chosen

    Available algorithms are:

    - ``standard`` -- the schoolbook algorithm;

    - ``finotti`` -- Finotti's algorithm; it can be used when the coefficient
      ring has characteristic `p`;

    - ``phantom`` -- computes the ring laws using the phantom components
      using a lift of ``coefficient_ring``, assuming that it is either
      `\mathbb F_q` for a power `q` of `p`, or a polynomial ring on that field;

    - ``p_invertible`` -- uses some optimisations when `p` is invertible
      in the coefficient ring.

    EXAMPLES::

        sage: WittVectorRing(QQ, p=5)
        Ring of truncated 5-typical Witt vectors of length 1 over
        Rational Field

    ::

        sage: WittVectorRing(GF(3))
        Ring of truncated 3-typical Witt vectors of length 1 over
        Finite Field of size 3

    ::

        sage: WittVectorRing(GF(3)['t'])
        Ring of truncated 3-typical Witt vectors of length 1 over
        Univariate Polynomial Ring in t over Finite Field of size 3

    ::

        sage: WittVectorRing(Qp(7), prec=30, p=5)
        Ring of truncated 5-typical Witt vectors of length 30 over
        7-adic Field with capped relative precision 20

    ::

        sage: p = 5
        sage: W = WittVectorRing(ZZ, prec=3, p=p)
        sage: W
        Ring of truncated 5-typical Witt vectors of length 3 over Integer Ring
        sage: V = W.verschiebung(extend=True)
        sage: WW = V.codomain()
        sage: WW
        Ring of truncated 5-typical Witt vectors of length 4 over Integer Ring
        sage: F = WW.frobenius_morphism()
        sage: w = W.random_element(x=-100, y=100)
        sage: w  # random
        (-13, 20, -33)
        sage: p*w == F(V(w))
        True

    ::

        sage: p = 7
        sage: W = WittVectorRing(GF(p)['x,y'], prec=5)
        sage: W
        Ring of truncated 7-typical Witt vectors of length 5 over Multivariate
        Polynomial Ring in x, y over Finite Field of size 7
        sage: V = W.verschiebung()
        sage: F = W.frobenius_morphism()
        sage: w = W.random_element()
        sage: w  # random
        (-3*x^2 + 2*y^2 - 3*y + 2, -2*x*y - y^2 - y - 2, -x^2 + y^2,
        x^2 + 2*y^2 + 3*x + y, 2*x^2 - 3*x*y + y^2 + y + 1)
        sage: p*w == F(V(w))
        True
        sage: p*w == V(F(w))
        True

    TESTS::

        sage: A = SymmetricGroup(3).algebra(QQ)
        sage: WittVectorRing(A)
        Traceback (most recent call last):
        ...
        TypeError: Symmetric group algebra of order 3 over Rational Field is not a commutative ring

        sage: WittVectorRing(QQ)
        Traceback (most recent call last):
        ...
        ValueError: Rational Field has non-prime characteristic and no prime was supplied

        sage: WittVectorRing(QQ, p=5, algorithm='moon')
        Traceback (most recent call last):
        ...
        ValueError: algorithm must be one of None, 'standard', 'p_invertible', 'finotti', 'phantom'

        sage: W = WittVectorRing(ZZ, p=53)
        sage: type(W)
        <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_standard_with_category'>

        sage: W = WittVectorRing(PolynomialRing(GF(13), 't'))
        sage: type(W)
        <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_phantom_with_category'>

        sage: W = WittVectorRing(PolynomialRing(GF(5), 't,u'))
        sage: type(W)
        <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_finotti_with_category'>
    """
    def __classcall_private__(cls, coefficient_ring, prec=1, p=None, algorithm=None):
        r"""
        Construct the ring of truncated Witt vectors from the parameters.

        TESTS::

            sage: W = WittVectorRing(QQ, p=5)
            sage: W
            Ring of truncated 5-typical Witt vectors of length 1 over Rational Field
        """
        if coefficient_ring not in CommutativeRings():
            raise TypeError(f"{coefficient_ring} is not a commutative ring")
        elif not isinstance(prec, (int, Integer)):
            raise TypeError(f"{prec} is not an integer")
        elif prec <= 0:
            raise ValueError(f"{prec} must be positive")

        prec = Integer(prec)
        char = coefficient_ring.characteristic()

        if p is None:
            if char not in Primes():
                raise ValueError(f"{coefficient_ring} has non-prime "
                                 "characteristic and no prime was supplied")
            p = char
        elif p not in Primes():
            raise ValueError(f"p must be a prime number, here {p} was given")

        match algorithm:
            case None:
                if p == char:
                    if (coefficient_ring in Fields().Finite()
                        or isinstance(coefficient_ring,
                                      PolynomialRing_generic)
                        and coefficient_ring.base()
                            in Fields().Finite()):
                        child = WittVectorRing_phantom
                    else:
                        child = WittVectorRing_finotti
                elif coefficient_ring(p).is_unit():
                    child = WittVectorRing_pinvertible
                else:
                    child = WittVectorRing_standard
            case 'finotti':
                child = WittVectorRing_finotti
            case 'phantom':
                child = WittVectorRing_phantom
            case 'p_invertible':
                child = WittVectorRing_pinvertible
            case 'standard':
                child = WittVectorRing_standard
            case _:
                raise ValueError("algorithm must be one of None, 'standard', "
                                 "'p_invertible', 'finotti', 'phantom'")

        return child.__classcall__(child, coefficient_ring, prec, p)

    def __init__(self, coefficient_ring, prec, prime, base) -> None:
        r"""
        Initialise ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(PolynomialRing(GF(5), 't'), prec=4); W
            Ring of truncated 5-typical Witt vectors of length 4 over Univariate Polynomial Ring in t over Finite Field of size 5
            sage: type(W)
            <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_phantom_with_category'>

            sage: TestSuite(W).run()
        """
        cring = coefficient_ring
        self._coefficient_ring = cring
        self._prec = prec
        self._prime = prime

        if prec.is_one() and cring in Fields():
            cat = Fields()
        elif prec.is_one() and cring in IntegralDomains():
            cat = IntegralDomains()
        else:
            cat = CommutativeRings()

        Parent.__init__(self, base=base, category=cat)

    def __iter__(self) -> Iterator:
        """
        Iterator for truncated Witt vector rings.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), p=3, prec=2)
            sage: [w for w in W]
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1),
            (2, 2)]
        """
        for t in product(self._coefficient_ring, repeat=self._prec):
            yield self(t)

    def _coerce_map_from_(self, S):
        """"
        Check whether there is a coerce map from ``S``.

        EXAMPLES::

            sage: K = GF(7)
            sage: WK = WittVectorRing(K, prec=3)
            sage: W = WittVectorRing(PolynomialRing(K, 't'), prec=3)
            sage: Wf = WittVectorRing(PolynomialRing(K, 't'), prec=3, algorithm='finotti')
            sage: WW = WittVectorRing(PolynomialRing(K, 't,u'), prec=3)
            sage: Wf.has_coerce_map_from(WK)  # indirect doctest
            True
            sage: Wf.has_coerce_map_from(W)  # indirect doctest
            False
            sage: WW.has_coerce_map_from(W)  # indirect doctest
            True
            sage: WW.has_coerce_map_from(Wf)  # indirect doctest
            True
            sage: WW.has_coerce_map_from(WK)  # indirect doctest
            True

            sage: W = WittVectorRing(GF(25), prec=2)
            sage: W.has_coerce_map_from(WittVectorRing(GF(5), prec=3))  # indirect doctest
            True
            sage: W.has_coerce_map_from(WittVectorRing(ZZ, p=5, prec=3))  # indirect doctest
            True
            sage: W.has_coerce_map_from(WittVectorRing(PolynomialRing(GF(5), 't,u'), prec=3))  # indirect doctest
            False
            sage: WW = WittVectorRing(PolynomialRing(GF(5), 't'), prec=3)
            sage: WW.has_coerce_map_from(W)  # indirect doctest
            False
            sage: W.has_coerce_map_from(WW)  # indirect doctest
            False

            sage: W = WittVectorRing(QQ, p=3, prec=3)
            sage: W.has_coerce_map_from(WittVectorRing(ZZ, p=3, prec=3))  # indirect doctest
            True
            sage: W.has_coerce_map_from(WittVectorRing(QQ, p=3, prec=2))  # indirect doctest
            False

            sage: W = WittVectorRing(PolynomialRing(ZZ, 'x'), p=5, prec=2)
            sage: W.has_coerce_map_from(WittVectorRing(ZZ, p=5, prec=3))  # indirect doctest
            True
            sage: W.has_coerce_map_from(WittVectorRing(ZZ, p=3, prec=3))  # indirect doctest
            False
        """
        if (isinstance(S, WittVectorRing)
            and S.precision() >= self._prec and S.prime() == self._prime
            and self._coefficient_ring.has_coerce_map_from(
                S.coefficient_ring())):
            return (any(isinstance(S, rng) for rng in self._always_coerce)
                    or (S.precision() != self._prec
                        or S.coefficient_ring() is not self._coefficient_ring)
                    and any(isinstance(S, rng)
                            for rng in self._coerce_when_different))
        if S is ZZ:
            return True

    def _generate_generators(self):
        """
        Generate a list of generators of ``self`` as a ring.

        EXAMPLES::

            sage: R.<x,y> = GF(2)[]
            sage: W = WittVectorRing(R, p=3, prec=3)
            sage: W.gens()  # indirect doctest
            ((x, 0, 0), (y, 0, 0), (0, 1, 0), (0, x, 0), (0, x^2, 0),
            (0, y, 0), (0, x*y, 0), (0, x^2*y, 0), (0, y^2, 0), (0, x*y^2, 0),
            (0, x^2*y^2, 0), (0, 0, 1), (0, 0, x), (0, 0, x^2), (0, 0, x^3),
            (0, 0, x^4), (0, 0, x^5), (0, 0, x^6), (0, 0, x^7), (0, 0, x^8),
            (0, 0, y), (0, 0, x*y), (0, 0, x^2*y), (0, 0, x^3*y),
            (0, 0, x^4*y), (0, 0, x^5*y), (0, 0, x^6*y), (0, 0, x^7*y),
            (0, 0, x^8*y), (0, 0, y^2), (0, 0, x*y^2), (0, 0, x^2*y^2),
            (0, 0, x^3*y^2), (0, 0, x^4*y^2), (0, 0, x^5*y^2), (0, 0, x^6*y^2),
            (0, 0, x^7*y^2), (0, 0, x^8*y^2), (0, 0, y^3), (0, 0, x*y^3),
            (0, 0, x^2*y^3), (0, 0, x^3*y^3), (0, 0, x^4*y^3), (0, 0, x^5*y^3),
            (0, 0, x^6*y^3), (0, 0, x^7*y^3), (0, 0, x^8*y^3), (0, 0, y^4),
            (0, 0, x*y^4), (0, 0, x^2*y^4), (0, 0, x^3*y^4), (0, 0, x^4*y^4),
            (0, 0, x^5*y^4), (0, 0, x^6*y^4), (0, 0, x^7*y^4), (0, 0, x^8*y^4),
            (0, 0, y^5), (0, 0, x*y^5), (0, 0, x^2*y^5), (0, 0, x^3*y^5),
            (0, 0, x^4*y^5), (0, 0, x^5*y^5), (0, 0, x^6*y^5), (0, 0, x^7*y^5),
            (0, 0, x^8*y^5), (0, 0, y^6), (0, 0, x*y^6), (0, 0, x^2*y^6),
            (0, 0, x^3*y^6), (0, 0, x^4*y^6), (0, 0, x^5*y^6), (0, 0, x^6*y^6),
            (0, 0, x^7*y^6), (0, 0, x^8*y^6), (0, 0, y^7), (0, 0, x*y^7),
            (0, 0, x^2*y^7), (0, 0, x^3*y^7), (0, 0, x^4*y^7), (0, 0, x^5*y^7),
            (0, 0, x^6*y^7), (0, 0, x^7*y^7), (0, 0, x^8*y^7), (0, 0, y^8),
            (0, 0, x*y^8), (0, 0, x^2*y^8), (0, 0, x^3*y^8), (0, 0, x^4*y^8),
            (0, 0, x^5*y^8), (0, 0, x^6*y^8), (0, 0, x^7*y^8), (0, 0, x^8*y^8))
        """
        coeff_ring_gens = self._coefficient_ring.gens()
        coeff_ring_names = self._coefficient_ring.variable_names()
        self._gens = tuple(self.teichmuller_lift(x) for x in coeff_ring_gens)
        names = tuple(f"T{x}" for x in coeff_ring_names)
        p = self._prime
        R = PolynomialRing(self._coefficient_ring, len(coeff_ring_gens), 'T')
        var = R.gens()
        vec = [self._coefficient_ring.zero()] * self._prec

        if (len(var) == 1 and
                coeff_ring_gens[0] == self._coefficient_ring.one() and
                (self._coefficient_ring.characteristic() == p
                 or self.base_ring().coefficient_ring()
                 is self._coefficient_ring)):
            self._assign_names(names)
            return

        for i in range(1, self._prec):
            p_i = p**i
            degrees = [p_i] * len(var)
            power = -1

            for powers in Tuples(IntegerRange(p_i), len(var)):
                power += 1
                if (self._coefficient_ring.characteristic() == p and
                        any((po % p).is_zero() and
                            not po.is_zero() for po in powers)):
                    continue

                vec[i] = prod(coeff_ring_gens[j]**powers[j] for j in range(len(var)))

                if (vec[i].is_one() and
                       self._coefficient_ring.characteristic() == p):
                    continue

                w = self(vec)

                if w not in self._gens:
                    self._gens += (w,)
                    name = f"V{i}T"
                    for j in range(len(var)):
                        if not powers[j].is_zero():
                            name += f"{coeff_ring_names[j]}"
                    name += f"_{power}"
                    names += (name,)

            vec[i] = self._coefficient_ring.zero()

        self._assign_names(names)

    def _generate_witt_polynomials(self, coefficient_ring, prec, p):
        """
        Generate the sum and product polynomials defining the ring laws of
        truncated Witt vectors for the ``standard`` algorithm.

        EXAMPLES::

            sage: P.<X1,X2,Y1,Y2> = PolynomialRing(GF(3),'X1,X2,Y1,Y2')
            sage: W = WittVectorRing(P, p=3, prec=2, algorithm='standard')
            sage: W([X1,X2]) + W([Y1,Y2])  # indirect doctest
            (X1 + Y1, -X1^2*Y1 - X1*Y1^2 + X2 + Y2)
            sage: W([X1,X2]) * W([Y1,Y2])  # indirect doctest
            (X1*Y1, X2*Y1^3 + X1^3*Y2)
            sage: W.frobenius_morphism()(W([X1,X2]))  # indirect doctest
            (X1^3, X2^3)
        """
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
        if p**(prec - 1) >= 2**16:
            implementation = 'generic'
        else:
            implementation = 'singular'

        # We first generate the "universal" polynomials and then project
        # to the coefficient ring.
        R = PolynomialRing(ZZ, var_names, implementation=implementation)
        x_y_vars = R.gens()
        x_vars = x_y_vars[:prec]
        y_vars = x_y_vars[prec:]

        self._sum_polynomials = [0]*(prec)
        for n in range(prec):
            s_n = x_vars[n] + y_vars[n]
            for i in range(n):
                s_n += ((x_vars[i]**(p**(n-i)) + y_vars[i]**(p**(n-i))
                        - self._sum_polynomials[i]**(p**(n-i))) / p**(n-i))
            self._sum_polynomials[n] = R(s_n)

        self._prod_polynomials = [x_vars[0] * y_vars[0]] + [0]*(prec-1)
        for n in range(1, prec):
            x_poly = sum([p**i * x_vars[i]**(p**(n-i)) for i in range(n+1)])
            y_poly = sum([p**i * y_vars[i]**(p**(n-i)) for i in range(n+1)])
            p_poly = sum([p**i * self._prod_polynomials[i]**(p**(n-i))
                         for i in range(n)])
            p_n = (x_poly*y_poly - p_poly) // p**n
            self._prod_polynomials[n] = p_n

        R = PolynomialRing(ZZ, x_var_names, implementation=implementation)
        x_vars = R.gens()

        self._frob_polynomials = []
        if not prec.is_one():
            self._frob_polynomials = [x_vars[0]**p + p*x_vars[1]]
            for n in range(2, prec):
                x_poly = sum([p**i * x_vars[i]**(p**(n-i)) for i in range(n+1)])
                p_poly = sum([p**i * self._frob_polynomials[i]**(p**(n-1-i))
                              for i in range(n-1)])
                self._frob_polynomials.append((x_poly - p_poly) / p**(n-1))

        S = PolynomialRing(coefficient_ring, x_y_vars)
        for n in range(prec):
            self._sum_polynomials[n] = S(self._sum_polynomials[n])
            self._prod_polynomials[n] = S(self._prod_polynomials[n])
        S = PolynomialRing(coefficient_ring, x_vars)
        for n in range(prec-1):
            self._frob_polynomials[n] = S(self._frob_polynomials[n])

    def _latex_(self) -> str:
        r"""
        Return a `\LaTeX` representation of ``self``.

        .. WARNING::

            This representation follows the standard representation in the
            literature which does not mention `p`.

        EXAMPLES::

            sage: W = WittVectorRing(PolynomialRing(GF(3),'t'))
            sage: latex(W)
            W_{1}\left(\Bold{F}_{3}[t]\right)
        """
        return "W_{%s}\\left(%s\\right)" % (latex(self._prec),
                                            latex(self._coefficient_ring))

    def _repr_(self) -> str:
        """
        Return a string representation of the ring.

        EXAMPLES::

            sage: WittVectorRing(QQ, p=2, prec=5)
            Ring of truncated 2-typical Witt vectors of length 5 over Rational Field
        """
        return f"Ring of truncated {self._prime}-typical Witt vectors of "\
               f"length {self._prec} over {self._coefficient_ring}"

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: WittVectorRing(GF(17), prec=2).cardinality()
            289
            sage: WittVectorRing(QQ, p=2).cardinality()
            +Infinity
        """
        return self._coefficient_ring.cardinality()**(self._prec)

    def characteristic(self):
        """
        Return the characteristic of ``self``.

        EXAMPLES::

            sage: WittVectorRing(GF(25), p=5, prec=3).characteristic()
            125
            sage: WittVectorRing(ZZ, p=2, prec=4).characteristic()
            0
            sage: WittVectorRing(Integers(18), p=3, prec=3).characteristic()
            162
        """
        p = self._prime
        if self._coefficient_ring(p).is_unit():
            # If p is invertible, W_n(R) is isomorphic to R^n.
            return self._coefficient_ring.characteristic()

        # This is Jacob Dennerlein's Corollary 3.3. in "Computational
        # Aspects of Mixed Characteristic Witt Vectors" (preprint)
        return p**(self._prec-1) * self._coefficient_ring.characteristic()

    def coefficient_ring(self):
        """
        Return the coefficient ring of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(Zp(5), p=5)
            sage: W.coefficient_ring()
            5-adic Ring with capped relative precision 20
        """
        return self._coefficient_ring

    def frobenius_morphism(self, truncate=False):
        """
        Return the Frobenius homomorphism of this ring of Witt vectors.

        INPUT:

        - ``truncate`` -- boolean (default: ``False``); when the coefficient
          ring has positive characteristic `p`, chose whether to lower the
          precision of the codomain by `1`. This is always the case when the
          characteristic of the coefficient ring is not `p`.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: W = WittVectorRing(R, p=5, prec=2)
            sage: w = W([x, y])
            sage: frob = W.frobenius_morphism()
            sage: frob(w)
            (x^5 + 5*y)

        TESTS::

            sage: TestSuite(frob).run(skip="_test_category")
        """
        return WittVectorFrobeniusMorphism(self, truncate=truncate)

    def frobenius_polynomials(self, variables=None):
        """
        Return the Witt Frobenius polynomials.

        INPUT:

        - ``variables`` -- names of the indeterminates (default: ``None``),
          given as a string, or as a list of strings, whose length must be the
          precision of the ring. When nothing is given, variables indexed by
          `X` are used.

        EXAMPLES::

            sage: W = WittVectorRing(GF(5), prec=4)
            sage: W.frobenius_polynomials()
            [X0^5, X1^5, X2^5]

            sage: W = WittVectorRing(ZZ, p=2, prec=3)
            sage: W.frobenius_polynomials('T0, T1, T2')
            [T0^2 + 2*T1, -2*T0^2*T1 - T1^2 + 2*T2]

            sage: W = WittVectorRing(ZZ, p=3, prec=3)
            sage: W.frobenius_polynomials('T0, T1, T2')
            [T0^3 + 3*T1, -3*T0^6*T1 - 9*T0^3*T1^2 - 8*T1^3 + 3*T2]
        """
        if not hasattr(self, '_frob_polynomials'):
            self._generate_witt_polynomials(self._coefficient_ring, self._prec,
                                            self._prime)
        if variables is None:
            return self._frob_polynomials.copy()
        R = PolynomialRing(self._coefficient_ring, variables)
        return [R(self._frob_polynomials[i]) for i in range(self._prec-1)]

    def gen(self, n=0):
        """
        Return the ``n``-th generator of ``self``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: W = WittVectorRing(R, p=3, prec=3)
            sage: W.gen()
            (x, 0, 0)
            sage: W.gen(n=12)
            (0, 0, x^8)
            sage: W.gen(n=13)
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.gens()[int(n)]

    def gens(self):
        """
        Return a tuple of generators of ``self`` as a ring.

        .. NOTE::

            This tuple is not necessarily minimal.

        EXAMPLES::

            sage: W = WittVectorRing(GF(29), prec=3)
            sage: W.gens()
            ((1, 0, 0),)

            sage: F.<a> = GF(27)
            sage: W = WittVectorRing(F, prec=3)
            sage: W.gens()
            ((a, 0, 0), (0, a, 0), (0, a^2, 0), (0, 0, a), (0, 0, a^2),
            (0, 0, a^2 + 2*a), (0, 0, 2*a^2 + a + 2), (0, 0, a^2 + 2*a + 2),
            (0, 0, 2*a^2 + 2))

            sage: R.<x> = ZZ[]
            sage: W = WittVectorRing(R, p=7, prec=2)
            sage: W.gens()
            ((x, 0), (0, 1), (0, x), (0, x^2), (0, x^3), (0, x^4), (0, x^5),
            (0, x^6))
        """
        if not hasattr(self, "_gens"):
            self._generate_generators()

        return self._gens

    def is_exact(self):
        """
        Return whether the elements of ``self`` are represented exactly.

        EXAMPLES::

            sage: WittVectorRing(ZZ, p=3, prec=3).is_exact()
            True
            sage: WittVectorRing(Qp(7), p=7, prec=2).is_exact()
            False
        """
        return self._coefficient_ring.is_exact()

    def is_field(self, proof=True):
        """
        Return ``True`` if this ring is a field.

        INPUT:

        - ``proof`` -- boolean (default: ``True``); if set to ``True``, the
          returned value is correct but the method might throw an exception.
          Otherwise, if it is set to ``False``, the method returns ``True`` if
          it can establish that ``self`` is a field, and ``False`` otherwise.

        EXAMPLES::

                sage: WittVectorRing(QQ, p=2).is_field()
                True
                sage: WittVectorRing(QQ, p=2, prec=2).is_field()
                False
                sage: WittVectorRing(GF(9, 'a')).is_field()
                True
                sage: WittVectorRing(GF(9, 'a'), prec=3).is_field()
                False
                sage: WittVectorRing(ZZ, p=5).is_field()
                False
                sage: L.<z> = LazyLaurentSeriesRing(ZZ)
                sage: W = WittVectorRing(L, p=7, algorithm="standard")
                sage: W.is_field(proof=True)
                Traceback (most recent call last):
                ...
                NotImplementedError: unable to determine whether or not Lazy
                Laurent Series Ring in z over Integer Ring is a field.
                sage: W.is_field(proof=False)
                False
        """
        return self._prec.is_one() and self._coefficient_ring.is_field(proof=proof)

    def is_finite(self) -> bool:
        """
        Return whether ``self`` is a finite ring.

        EXAMPLES::

            sage: WittVectorRing(GF(23)).is_finite()
            True
            sage: WittVectorRing(ZZ, p=2).is_finite()
            False
        """
        return self._coefficient_ring.is_finite()

    def is_integral_domain(self, proof=True):
        """
        Return ``True`` if this ring is an integral domain.

        INPUT:

        - ``proof`` -- boolean (default: ``True``); if set to ``True``, the
          returned value is correct but the method might throw an exception.
          Otherwise, if it is set to ``False``, the method returns ``True`` if
          it can establish that ``self`` is an integral domain, and ``False``
          otherwise.

        EXAMPLES::

                sage: WittVectorRing(ZZ, p=2).is_integral_domain()
                True
                sage: WittVectorRing(ZZ, p=2, prec=2).is_integral_domain()
                False
                sage: WittVectorRing(GF(9, 'a')).is_integral_domain()
                True
                sage: WittVectorRing(GF(9, 'a'), prec=3).is_integral_domain()
                False
                sage: R.<x> = ZZ[]
                sage: S.<y> = R.quo(x^2 + 7*x + 8)
                sage: T.<t> = S[]
                sage: Q.<z> = T.quo(t^2)
                sage: W = WittVectorRing(Q, p=7, algorithm="standard")
                sage: W.is_integral_domain(proof=True)
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot rewrite Univariate Quotient
                Polynomial Ring in y over Integer Ring with modulus
                x^2 + 7*x + 8 as an isomorphic ring
                sage: W.is_integral_domain(proof=False)
                False
        """
        return (self._prec.is_one() and
                self._coefficient_ring.is_integral_domain(proof=proof))

    def is_integrally_closed(self):
        """
        Return ``True`` if this ring is an integrally closed domain.

        EXAMPLES::

                sage: x = polygen(ZZ, 'x')
                sage: K.<a> = NumberField(x^2 + 189*x + 394)
                sage: R = K.order(2*a)
                sage: WittVectorRing(R, p=79).is_integrally_closed()
                False
                sage: WittVectorRing(ZZ, p=101).is_integrally_closed()
                True
                sage: WittVectorRing(ZZ, p=101, prec=2).is_integrally_closed()
                False
                sage: S.<t> = ZZ[]
                sage: WittVectorRing(S, p=101).is_integrally_closed()
                Traceback (most recent call last):
                ...
                NotImplementedError
        """
        return (self._prec.is_one() and
                self._coefficient_ring.is_integrally_closed())

    def is_prime_field(self):
        r"""
        Return ``True`` if ``self`` is isomorphic to one of the prime fields
        `\mathbb Q` or `\mathbb F_p`.

        EXAMPLES::

                sage: WittVectorRing(QQ, p=2).is_prime_field()
                True
                sage: WittVectorRing(QQ, p=2, prec=2).is_prime_field()
                False
                sage: WittVectorRing(GF(9, 'a')).is_prime_field()
                False
                sage: WittVectorRing(GF(7)).is_prime_field()
                True
                sage: WittVectorRing(GF(7), prec=3).is_prime_field()
                False
        """
        return self._prec.is_one() and self._coefficient_ring.is_prime_field()

    def ngens(self):
        """
        Return the number of generators of ``self``.

        .. NOTE::

            This number is not necessarily minimal.

        EXAMPLES::

            sage: W = WittVectorRing(GF(2), p=3, prec=3)
            sage: W.ngens()
            1
            sage: W = WittVectorRing(GF(2)['t'], p=3, prec=3)
            sage: W.ngens()
            13
            sage: W = WittVectorRing(GF(3), p=3, prec=3)
            sage: W.ngens()
            1
            sage: W = WittVectorRing(GF(3)['t'], p=3, prec=3)
            sage: W.ngens()
            9
            sage: W = WittVectorRing(GF(9,'a'), p=3, prec=3)
            sage: W.ngens()
            8
            sage: W = WittVectorRing(GF(9,'a')['t'], p=3, prec=3)
            sage: W.ngens()
            9
        """
        return len(self.gens())

    def precision(self):
        """
        Return the length of the truncated Witt vectors in ``self``.

        EXAMPLES::

            sage: WittVectorRing(GF(9), p=3, prec=3).precision()
            3
        """
        return self._prec

    def prime(self):
        """
        Return the prime from which the truncated Witt vector ring has been
        constructed.

        EXAMPLES::

            sage: W = WittVectorRing(GF(81), prec=3)
            sage: W.prime()
            3

            sage: W = WittVectorRing(ZZ, p=7, prec=2)
            sage: W.prime()
            7
        """
        return self._prime

    def prod_polynomials(self, variables=None):
        """
        Return the Witt product polynomials.

        INPUT:

        - ``variables`` -- names of the indeterminates (default: ``None``),
          given as a string, or as a list of strings, whose length must be the
          double of the precision of the ring. When nothing is given,
          variables indexed by `X` and `Y` are used.

        EXAMPLES::

            sage: W = WittVectorRing(GF(5), prec=3)
            sage: W.prod_polynomials()
            [X0*Y0,
            X1*Y0^5 + X0^5*Y1,
            -X0^5*X1^4*Y0^20*Y1 - 2*X0^10*X1^3*Y0^15*Y1^2 - 2*X0^15*X1^2*Y0^10*Y1^3 - X0^20*X1*Y0^5*Y1^4 + X2*Y0^25 + X0^25*Y2 + X1^5*Y1^5]

            sage: W = WittVectorRing(ZZ, p=2, prec=2)
            sage: W.prod_polynomials('T0, T1, U0, U1')
            [T0*U0, T1*U0^2 + T0^2*U1 + 2*T1*U1]
        """
        if not hasattr(self, '_prod_polynomials'):
            self._generate_witt_polynomials(self._coefficient_ring, self._prec,
                                            self._prime)
        if variables is None:
            return self._prod_polynomials.copy()
        R = PolynomialRing(self._coefficient_ring, variables)
        return [R(self._prod_polynomials[i]) for i in range(self._prec)]

    def random_element(self, *args, **kwds):
        """
        Return a random truncated Witt vector.

        Extra arguments are passed to
        the random generator of the coefficient ring.

        EXAMPLES::

            sage: WittVectorRing(GF(27,'t'), prec=2).random_element()  # random
            (z3, 2*z3^2 + 1)

            sage: W = WittVectorRing(PolynomialRing(ZZ,'x'), p=3, prec=3)
            sage: W.random_element(5)  # random
            (x^5 - 2*x^4 - 4*x^3 - 2*x^2 + 1, -x^5 + 2*x^4 - x - 1,
            -x^5 + 7*x^4 + 3*x^3 - 24*x^2 - 1)
        """
        return self(tuple(self._coefficient_ring.random_element(*args, **kwds)
                          for _ in range(self._prec)))

    def sum_polynomials(self, variables=None):
        """
        Return the Witt sum polynomials.

        INPUT:

        - ``variables`` -- names of the indeterminates (default: ``None``),
          given as a string, or as a list of strings, whose length must be the
          double of the precision of the ring. When nothing is given,
          variables indexed by `X` and `Y` are used.

        EXAMPLES::

            sage: W = WittVectorRing(GF(5), prec=2)
            sage: W.sum_polynomials(['T0', 'T1', 'U0', 'U1'])
            [T0 + U0, -T0^4*U0 - 2*T0^3*U0^2 - 2*T0^2*U0^3 - T0*U0^4 + T1 + U1]

            sage: W = WittVectorRing(ZZ, p=2, prec=3)
            sage: W.sum_polynomials()
            [X0 + Y0,
            -X0*Y0 + X1 + Y1,
            -X0^3*Y0 - 2*X0^2*Y0^2 - X0*Y0^3 + X0*X1*Y0 + X0*Y0*Y1 - X1*Y1 + X2 + Y2]
        """
        if not hasattr(self, '_sum_polynomials'):
            self._generate_witt_polynomials(self._coefficient_ring, self._prec,
                                            self._prime)
        if variables is None:
            return self._sum_polynomials.copy()
        R = PolynomialRing(self._coefficient_ring, variables)
        return [R(self._sum_polynomials[i]) for i in range(self._prec)]

    def teichmuller_lift(self, x):
        """
        Return the Teichmüller lift of ``x`` in ``self``.

        This lift is sometimes known as the multiplicative lift of ``x``.

        EXAMPLES::

            sage: WittVectorRing(GF(125,'t'), prec=2).teichmuller_lift(3)
            (3, 0)
        """
        if x not in self._coefficient_ring:
            raise TypeError(f"{x} not in {self._coefficient_ring}")
        return self((x,) + tuple(0 for _ in range(self._prec-1)))

    def verschiebung(self, extend=False):
        """
        Return the Verschiebung map of this ring of Witt vectors.

        - ``extend`` -- boolean (default: ``False``); whether the codomain of
          the map has precision `1` more than ``self``.

        INPUT:

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: W = WittVectorRing(R, p=5, prec=2)
            sage: w = W([x, y])
            sage: V = W.verschiebung()
            sage: V(w)
            (0, x)

        TESTS::

            sage: TestSuite(V).run(skip="_test_category")
        """
        return WittVectorVerschiebung(self, extend=extend)


class WittVectorRing_finotti(WittVectorRing):
    """
    Child class for truncated Witt vectors using Finotti's algorithm.

    .. WARNING::

        This class should never be called directly, use ``WittVectorRing``
        instead.

    EXAMPLES::

        sage: W = WittVectorRing(GF(49), prec=3, algorithm='finotti')
        sage: W
        Ring of truncated 7-typical Witt vectors of length 3 over Finite Field in z2 of size 7^2

        sage: W = WittVectorRing(ZZ, p=11, prec=3, algorithm='finotti')
        Traceback (most recent call last):
        ...
        ValueError: the 'finotti' algorithm only works for coefficients rings of characteristic p
    """
    Element = WittVector_finotti

    def __init__(self, coefficient_ring, prec, prime) -> None:
        r"""
        Initialise ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(PowerSeriesRing(GF(3), 't'), p=3, prec=3)
            sage: W
            Ring of truncated 3-typical Witt vectors of length 3 over Power Series Ring in t over Finite Field of size 3
            sage: type(W)
            <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_finotti_with_category'>

            sage: TestSuite(W).run()
        """
        if coefficient_ring.characteristic() != prime:
            raise ValueError("the 'finotti' algorithm only works for "
                             "coefficients rings of characteristic p")

        if isinstance(coefficient_ring, MPolynomialRing_base):
            self._always_coerce = [WittVectorRing_finotti,
                                   WittVectorRing_phantom,
                                   WittVectorRing_standard]
            self._coerce_when_different = []
        else:
            self._always_coerce = [WittVectorRing_finotti,
                                   WittVectorRing_standard]
            self._coerce_when_different = [WittVectorRing_phantom]

        import numpy as np
        R = Zp(prime, prec=prec+1, type='fixed-mod')
        v_p = ZZ.valuation(prime)
        table = [[0]]
        for k in range(1, prec+1):
            pk = prime**k
            row = np.empty(pk, dtype=int)
            row[0] = 0
            prev_bin = 1
            for i in range(1, pk // 2 + 1):
                val = v_p(i)
                # Instead of calling binomial each time, we compute the
                # coefficients recursively. This is MUCH faster.
                next_bin = prev_bin * (pk - (i-1)) // i
                prev_bin = next_bin
                series = R(-next_bin // prime**(k-val))
                for _ in range(val):
                    temp = series % prime
                    series = (series - R.teichmuller(temp)) // prime
                row[i] = ZZ(series % prime)
                row[pk - i] = row[i]  # binomial coefficients are symmetric
            table.append(row)
        self._binomial_table = table

        if coefficient_ring.base_ring() is coefficient_ring:
            base = self
        else:
            base = WittVectorRing(coefficient_ring.base_ring(), prec=prec,
                                  p=prime, algorithm="finotti")

        super().__init__(coefficient_ring, prec, prime, base)

    def _eta_bar(self, vec, eta_index):
        r"""
        Generate the `\eta_i` for ``finotti``'s algorithm.

        EXAMPLES::

            sage: R.<x,y,z,t> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=2, algorithm='finotti')
            sage: (W([x,y]) + W([z,t]))[1]  # indirect doctest
            -x^4*z - 2*x^3*z^2 - 2*x^2*z^3 - x*z^4 + y + t
        """
        vec = tuple(x for x in vec if x != 0)  # strip zeroes

        # special cases
        if len(vec) <= 1:
            return 0
        if eta_index == 0:
            return sum(vec)

        # renaming to match notation in Finotti's "Computations with Witt
        # vectors and the Greenberg transform", doi:10.1142/S1793042114500377
        k = eta_index
        p = self._prime
        # if vec = (x,y), we know what to do: Theorem 8.6
        if len(vec) == 2:
            # Here we have to check if we've pre-computed already
            x, y = vec
            scriptN = [[None] for _ in range(k+1)]  # each list starts with
            # None, so that indexing matches paper

            # calculate first N_t scriptN's
            for t in range(1, k+1):
                for i in range(1, p**t):
                    scriptN[t].append(self._binomial_table[t][i]
                                      * fast_char_p_power(x, i)
                                      * fast_char_p_power(y, p**t - i))
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


class WittVectorRing_phantom(WittVectorRing):
    """
    Child class for truncated Witt vectors using the ``phantom`` algorithm.

    .. WARNING::

        This class should never be called directly, use ``WittVectorRing``
        instead.

    EXAMPLES::

        sage: W = WittVectorRing(GF(19), prec=20)
        sage: W
        Ring of truncated 19-typical Witt vectors of length 20 over Finite Field of size 19

        sage: W = WittVectorRing(QQ, p=23, prec=3, algorithm='phantom')
        Traceback (most recent call last):
        ...
        ValueError: the 'phantom' algorithm only works when the coefficient ring is a finite field of char. p, or a polynomial ring on that field
    """
    Element = WittVector_phantom

    def __init__(self, coefficient_ring, prec, prime) -> None:
        r"""
        Initialise ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(23,'t'), p=23, prec=2)
            sage: W
            Ring of truncated 23-typical Witt vectors of length 2 over Finite Field of size 23
            sage: type(W)
            <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_phantom_with_category'>

            sage: TestSuite(W).run()
        """
        msg = "the 'phantom' algorithm only works when the coefficient ring is"\
            " a finite field of char. p, or a polynomial ring on that field"

        if coefficient_ring.characteristic() != prime:
            raise ValueError(msg)

        if not (coefficient_ring in Fields().Finite() or
                (isinstance(coefficient_ring, (PolynomialRing_generic,
                                               MPolynomialRing_base)) and
                 coefficient_ring.base() in Fields().Finite())):
            raise ValueError(msg)

        if (coefficient_ring in Fields().Finite()
            or isinstance(coefficient_ring,
                          PolynomialRing_generic)):
            self._always_coerce = [WittVectorRing_finotti,
                                   WittVectorRing_phantom,
                                   WittVectorRing_standard]
            self._coerce_when_different = []
        else:
            self._always_coerce = [WittVectorRing_phantom,
                                   WittVectorRing_standard]
            self._coerce_when_different = [WittVectorRing_finotti]

        if coefficient_ring.base_ring() is coefficient_ring:
            base = self
        else:
            base = WittVectorRing(coefficient_ring.base_ring(), prec=prec,
                                  p=prime, algorithm="phantom")

        super().__init__(coefficient_ring, prec, prime, base)


class WittVectorRing_pinvertible(WittVectorRing):
    """
    Child class for truncated Witt vectors using the ``p_invertible`` algorithm.

    .. WARNING::

        This class should never be called directly, use ``WittVectorRing``
        instead.

    EXAMPLES::

        sage: W = WittVectorRing(QQ, p=31, prec=20)
        sage: W
        Ring of truncated 31-typical Witt vectors of length 20 over Rational Field

        sage: W = WittVectorRing(GF(3), prec=3, algorithm='p_invertible')
        Traceback (most recent call last):
        ...
        ValueError: the 'p_invertible' algorithm only works when p is a unit in the ring of coefficients
    """
    Element = WittVector_pinvertible

    def __init__(self, coefficient_ring, prec, prime) -> None:
        r"""
        Initialise ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(QQ, p=11, prec=3)
            sage: W
            Ring of truncated 11-typical Witt vectors of length 3 over Rational Field
            sage: type(W)
            <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_pinvertible_with_category'>

            sage: TestSuite(W).run()
        """
        if not coefficient_ring(prime).is_unit():
            raise ValueError("the 'p_invertible' algorithm only works when p "
                             "is a unit in the ring of coefficients")

        self._always_coerce = [WittVectorRing_pinvertible,
                               WittVectorRing_standard]
        self._coerce_when_different = []

        if coefficient_ring.base_ring() is coefficient_ring:
            base = self
        else:
            base = WittVectorRing(coefficient_ring.base_ring(), prec=prec,
                                  p=prime, algorithm="p_invertible")

        super().__init__(coefficient_ring, prec, prime, base)


class WittVectorRing_standard(WittVectorRing):
    """
    Child class for truncated Witt vectors using the ``standard`` algorithm.

    .. WARNING::

        This class should never be called directly, use ``WittVectorRing``
        instead.

    EXAMPLES::

        sage: W = WittVectorRing(GF(3), prec=3, algorithm='standard')
        sage: W
        Ring of truncated 3-typical Witt vectors of length 3 over Finite Field of size 3
    """
    Element = WittVector_standard

    def __init__(self, coefficient_ring, prec, prime) -> None:
        r"""
        Initialise ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(ZZ, p=5, prec=2)
            sage: W
            Ring of truncated 5-typical Witt vectors of length 2 over Integer Ring
            sage: type(W)
            <class 'sage.rings.padics.witt_vector_ring.WittVectorRing_standard_with_category'>

            sage: TestSuite(W).run()
        """
        self._always_coerce = []
        self._coerce_when_different = [WittVectorRing]

        self._generate_witt_polynomials(coefficient_ring, prec, prime)

        if coefficient_ring.base_ring() is coefficient_ring:
            base = self
        else:
            base = WittVectorRing(coefficient_ring.base_ring(), prec=prec,
                                  p=prime, algorithm="standard")

        super().__init__(coefficient_ring, prec, prime, base)


class WittVectorFrobeniusMorphism(RingHomomorphism):
    """
    A class implementing Frobenius homomorphisms on Witt vector rings.

    .. WARNING::

        This class should not be called directly, use
        :meth:`WittVectorRing.frobenius_morphism` instead.

    EXAMPLES::

        sage: R.<x,y,z> = PolynomialRing(ZZ)
        sage: W = WittVectorRing(R, p=3, prec=3)
        sage: w = W([x, y, z])
        sage: frob = W.frobenius_morphism()
        sage: frob(w)
        (x^3 + 3*y, -3*x^6*y - 9*x^3*y^2 - 8*y^3 + 3*z)

        sage: W = WittVectorRing(R, p=3, prec=1)
        sage: W.frobenius_morphism()
        Traceback (most recent call last):
        ...
        ValueError: the ring of Witt vectors must have precision at least 2
        when truncate=True

    TESTS::

        sage: TestSuite(frob).run(skip="_test_category")
    """
    def __init__(self, domain, truncate=False):
        """
        INPUT:

        - ``domain`` -- the ring of Witt vectors whose Frobenius endomorphism
          we are computing.

        - ``truncate`` -- boolean (default: ``False``); when the coefficient
          ring of the domain has positive characteristic `p`, chose whether to
          lower the precision of the codomain by `1`. This is always the case
          when the characteristic of the coefficient ring is not `p`.

        EXAMPLES::

            sage: R.<x,y,z,t> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=4)
            sage: w = W([x, y, z, t])
            sage: frob = W.frobenius_morphism()
            sage: frob(w)
            (x^5, y^5, z^5, t^5)
            sage: frob = W.frobenius_morphism(truncate=True)
            sage: frob(w)
            (x^5, y^5, z^5)

        TESTS::

            sage: R.<x,y> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=2)
            sage: frob = W.frobenius_morphism()
            sage: TestSuite(frob).run(skip="_test_category")
            sage: frob = W.frobenius_morphism(truncate=True)
            sage: TestSuite(frob).run(skip="_test_category")
        """

        coeff_ring = domain.coefficient_ring()

        if not truncate:
            if (isinstance(domain, WittVectorRing_standard)
                    and coeff_ring.characteristic() != domain.prime()):
                truncate = True

            elif isinstance(domain, WittVectorRing_pinvertible):
                truncate = True

            else:
                codomain = domain
                self._call_ = self._call_char_p

        if truncate:
            prec = domain.precision()

            if prec.is_one():
                raise ValueError("the ring of Witt vectors must have "
                                 "precision at least 2 when truncate=True")

            if isinstance(domain, WittVectorRing_finotti):
                algorithm = "finotti"
                self._call_ = self._call_char_p_truncate
            elif isinstance(domain, WittVectorRing_phantom):
                algorithm = "phantom"
                self._call_ = self._call_phantom_truncate
            elif isinstance(domain, WittVectorRing_pinvertible):
                algorithm = "p_invertible"
                self._call_ = self._call_standard
            elif isinstance(domain, WittVectorRing_standard):
                algorithm = "standard"
                self._call_ = self._call_standard

            codomain = WittVectorRing(coeff_ring, prec=prec-1,
                                      p=domain.prime(), algorithm=algorithm)

        super().__init__(Hom(domain, codomain))

    def _call_standard(self, x):
        """
        Evaluate ``self`` at ``x`` using the general algorithm.

        EXAMPLES::

            sage: W = WittVectorRing(ZZ, p=11, prec=2)
            sage: w = W([1, 100])
            sage: frob = W.frobenius_morphism()
            sage: frob(w)
            (1101)

        TESTS::

            sage: TestSuite(frob).run(skip="_test_category")
        """
        W = self.domain()
        Wcod = self.codomain()
        return Wcod([W.frobenius_polynomials()[i](x.coordinates())
                     for i in range(Wcod.precision())])

    def _call_char_p_truncate(self, x):
        """
        Evaluate ``self`` at ``x`` when the coefficient ring has charaterisic `p`
        for the truncated Frobenius.

        EXAMPLES::

            sage: W = WittVectorRing(GF(7), prec=4, algorithm="finotti")
            sage: w = W([1, 2, 3, 4])
            sage: frob = W.frobenius_morphism(truncate=True)
            sage: frob(w)
            (1, 2, 3)

        TESTS::

            sage: TestSuite(frob).run(skip="_test_category")
        """
        Wcod = self.codomain()
        F = Wcod.coefficient_ring().frobenius_endomorphism()
        return Wcod([F(w) for w in x])

    def _call_char_p(self, x):
        """
        Evaluate ``self`` at ``x`` when the coefficient ring has charaterisic `p`
        for the non truncated Frobenius.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=3, algorithm="phantom")
            sage: w = W([x, x + y, x*y + y])
            sage: frob = W.frobenius_morphism()
            sage: frob(w)
            (x^5, x^5 + y^5, x^5*y^5 + y^5)

        TESTS::

            sage: TestSuite(frob).run(skip="_test_category")
        """
        W = self.domain()
        F = W.coefficient_ring().frobenius_endomorphism()
        return W([F(w) for w in x])

    def _call_phantom_truncate(self, x):
        """
        Evaluate ``self`` at ``x`` when the algorithm is ``phantom`` and the
        Frobenius is truncated.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=3, algorithm="phantom")
            sage: w = W([x, x + y, x*y + y])
            sage: frob = W.frobenius_morphism(truncate=True)
            sage: fw = frob(w)
            sage: fw.phantom(lift=True)
            (x^5 + 5*x + 5*y, x^25 + 5*x^5 + 5^2*x^4*y + 2*5^2*x^3*y^2 +
            2*5^2*x^2*y^3 + 5^2*x*y^4 + 5*y^5 + 5^2*x*y + 5^2*y)

        TESTS::

            sage: TestSuite(frob).run(skip="_test_category")
        """
        phantom = x.phantom(lift=True)
        W = self.domain()
        Wcod = self.codomain()
        return Wcod(phantom=[phantom[i] for i in range(1, W.precision())])

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(PolynomialRing(GF(3),'t'))
            sage: latex(W.frobenius_morphism())
            F
        """
        return "F"

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(11), prec=5)
            sage: W.frobenius_morphism()
            Frobenius endomorphism on the Ring of truncated 11-typical Witt
            vectors of length 5 over Finite Field of size 11
            sage: W.frobenius_morphism(truncate=True)
            Frobenius homomorphism from the Ring of truncated 11-typical Witt
            vectors of length 5 over Finite Field of size 11 to the Ring of
            truncated 11-typical Witt vectors of length 4 over Finite Field of
            size 11
        """
        if self.domain() is self.codomain():
            return f"Frobenius endomorphism on the {self.domain()}"
        return f"Frobenius homomorphism from the {self.domain()} to the "\
               f"{self.codomain()}"


class WittVectorVerschiebung(RingMap):
    """
    A class implementing Veschiebung maps on Witt vector rings.

    .. WARNING::

        This class should not be called directly, use
        :meth:`WittVectorRing.verschiebung` instead.

    EXAMPLES::

        sage: R.<x,y,z> = PolynomialRing(ZZ)
        sage: W = WittVectorRing(R, p=3, prec=3)
        sage: w = W([x, y, z])
        sage: V = W.verschiebung()
        sage: V(w)
        (0, x, y)

    TESTS::

        sage: TestSuite(V).run(skip="_test_category")
    """
    def __init__(self, domain, extend=False):
        """
        INPUT:

        - ``domain`` -- the ring of Witt vectors whose Verschiebung map we are
          computing.

        - ``extend`` -- boolean (default: ``False``); whether the codomain of
          the Verschiebung map is a ring of Witt vectors whose lengths are `1`
          more than the ones in the domain.

        EXAMPLES::

            sage: R.<x,y,z,t> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=4)
            sage: w = W([x, y, z, t])
            sage: V = W.verschiebung()
            sage: V(w)
            (0, x, y, z)
            sage: V = W.verschiebung(extend=True)
            sage: V(w)
            (0, x, y, z, t)

        TESTS::

            sage: R.<x,y> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=2)
            sage: V = W.verschiebung()
            sage: TestSuite(V).run(skip="_test_category")
            sage: V = W.verschiebung(extend=True)
            sage: TestSuite(V).run(skip="_test_category")
        """
        if extend:
            self._call_ = self._call_extend
            prec = domain.precision()

            if isinstance(domain, WittVectorRing_finotti):
                algorithm = "finotti"
            elif isinstance(domain, WittVectorRing_phantom):
                algorithm = "phantom"
            elif isinstance(domain, WittVectorRing_pinvertible):
                algorithm = "p_invertible"
            elif isinstance(domain, WittVectorRing_standard):
                algorithm = "standard"

            codomain = WittVectorRing(domain.coefficient_ring(), prec=prec+1,
                                      p=domain.prime(), algorithm=algorithm)

        else:
            self._call_ = self._call_no_extend
            codomain = domain

        super().__init__(domain.Hom(codomain, CommutativeAdditiveGroups()))

    def _call_extend(self, x):
        """
        Evaluate ``self`` at ``x`` for the extended Verschiebung.

        EXAMPLES::

            sage: W = WittVectorRing(ZZ, p=11, prec=2)
            sage: w = W([1, 100])
            sage: V_ZZ = W.verschiebung(extend=True)
            sage: V_ZZ(w)
            (0, 1, 100)

            sage: R.<x,y> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=3, algorithm="phantom")
            sage: w = W([x, x + y, x*y + y])
            sage: V_GF = W.verschiebung(extend=True)
            sage: V_GF(w)
            (0, x, x + y, x*y + y)

        TESTS::

            sage: TestSuite(V_ZZ).run(skip="_test_category")
            sage: TestSuite(V_GF).run(skip="_test_category")

        """
        Wcod = self.codomain()
        R = Wcod.coefficient_ring()
        return Wcod((R.zero(),) + x.coordinates())

    def _call_no_extend(self, x):
        """
        Evaluate ``self`` at ``x`` when the domain is the codomain.

        EXAMPLES::

            sage: W = WittVectorRing(GF(7), prec=4, algorithm="finotti")
            sage: w = W([1, 2, 3, 4])
            sage: V_GF = W.verschiebung()
            sage: V_GF(w)
            (0, 1, 2, 3)

            sage: R.<x,y> = PolynomialRing(GF(5))
            sage: W = WittVectorRing(R, prec=3, algorithm="phantom")
            sage: w = W([x, x + y, x*y + y])
            sage: V_poly = W.verschiebung()
            sage: V_poly(w)
            (0, x, x + y)

        TESTS::

            sage: TestSuite(V_GF).run(skip="_test_category")
            sage: TestSuite(V_poly).run(skip="_test_category")
        """
        W = self.domain()
        R = W.coefficient_ring()
        return W((R.zero(),) + x.coordinates())

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(PolynomialRing(GF(3),'t'))
            sage: latex(W.verschiebung())
            V
        """
        return "V"

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(11), prec=5)
            sage: W.verschiebung()
            Verschiebung map on the Ring of truncated 11-typical Witt vectors
            of length 5 over Finite Field of size 11
            sage: W.verschiebung(extend=True)
            Verschiebung map from the Ring of truncated 11-typical Witt vectors
            of length 5 over Finite Field of size 11 to the Ring of truncated
            11-typical Witt vectors of length 6 over Finite Field of size 11
        """
        if self.domain() is self.codomain():
            return f"Verschiebung map on the {self.domain()}"
        return f"Verschiebung map from the {self.domain()} to the "\
               f"{self.codomain()}"
