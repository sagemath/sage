r"""
Hyperelliptic curves (smooth model) over the rationals

AUTHORS:

- David Kohel (2006): initial version
- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                     2025 Sabrina Kunzweiler <sabrina.kunzweiler@math.u-bordeaux.fr>
#                     2025 Gareth Ma <grhkm21@gmail.com>
#                     2025 Giacomo Pope <giacomopope@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from multiprocessing import Pool

import sage.rings.abc
from sage.arith.misc import binomial, prime_range
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp as pAdicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.real_mpfr import RR
from sage.schemes.hyperelliptic_curves import hyperelliptic_generic


def _hypellfrob_euler_factor(args):
    p, N, coefficients, genus = args
    if N is None:
        from sage.rings.finite_rings.finite_field_constructor import GF
        from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve

        R = PolynomialRing(GF(p), "x")
        Q = R(coefficients)
        P = HyperellipticCurve(Q).frobenius_polynomial()
        return p, P.reverse()

    from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob

    R = PolynomialRing(ZZ, "x")
    Q = R(coefficients)
    M = hypellfrob(p, N, Q).change_ring(ZZ)
    coefficients = M.charpoly().list()[genus : 2 * genus + 1]
    modulus = p**N
    coefficients = [c % modulus for c in coefficients]
    coefficients = [c if 2 * c < modulus else c - modulus
                    for c in coefficients]
    for i in range(1, genus + 1):
        coefficient = coefficients[genus - i]
        bound = binomial(2 * genus, i)
        if coefficient**2 > bound**2 * p**i:
            raise ArithmeticError(
                "Frobenius coefficient is outside its Weil bound"
            )
    coefficients = [coefficients[genus - i] * p ** (genus - i)
                    for i in range(genus)] + coefficients
    P = PolynomialRing(ZZ, "T")(coefficients)
    return p, P.reverse()


def _hypellfrob_precision(p, genus):
    bound = 2 * binomial(2 * genus, genus) * RR(p).sqrt() ** genus
    N = ZZ(bound.ceil()).exact_log(p)
    if p**N <= bound:
        N += 1
    return N


class HyperellipticCurve_rational_field(
    hyperelliptic_generic.HyperellipticCurve_generic
):
    def __init__(
        self, projective_model, f, h, genus: Integer, names=["x", "y"]
    ) -> None:
        r"""
        Create a hyperelliptic curve over the rationals.

        TESTS::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(-x^2, x^3 + 1)
            sage: H
            Hyperelliptic Curve over Rational Field defined by y^2 + (x^3 + 1)*y = -x^2
        """
        super().__init__(projective_model, f, h, genus, names)

    def matrix_of_frobenius(self, p, prec=20):
        r"""
        Compute the matrix of Frobenius on Monsky-Washnitzer cohomology using
        the `p`-adic field with precision ``prec``.

        This function is essentially a wrapper function of
        :meth:`sage.schemes.hyperelliptic_curves.monsky_washnitzer.matrix_of_frobenius_hyperelliptic`.

        INPUT:

        - ``p`` (prime integer or pAdic ring / field ) -- if ``p`` is an integer,
          constructs a ``pAdicField`` with ``p`` to compute the matrix of
          Frobenius, otherwise uses the supplied pAdic ring or field.

        - ``prec`` (optional) -- if ``p`` is an prime integer, the `p`-adic
          precision of the coefficient ring constructed

        EXAMPLES::

            sage: K = pAdicField(5, prec=3)
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^5 - 2*x + 3)
            sage: H.matrix_of_frobenius(K)
            [            4*5 + O(5^3)       5 + 2*5^2 + O(5^3) 2 + 3*5 + 2*5^2 + O(5^3)     2 + 5 + 5^2 + O(5^3)]
            [      3*5 + 5^2 + O(5^3)             3*5 + O(5^3)             4*5 + O(5^3)         2 + 5^2 + O(5^3)]
            [    4*5 + 4*5^2 + O(5^3)     3*5 + 2*5^2 + O(5^3)       5 + 3*5^2 + O(5^3)     2*5 + 2*5^2 + O(5^3)]
            [            5^2 + O(5^3)       5 + 4*5^2 + O(5^3)     4*5 + 3*5^2 + O(5^3)             2*5 + O(5^3)]

        You can also pass directly a prime `p` with to construct a pAdic field with precision
        ``prec``::

            sage: H.matrix_of_frobenius(3, prec=2)
            [        O(3^2)     3 + O(3^2)         O(3^2)         O(3^2)]
            [    3 + O(3^2)         O(3^2)         O(3^2) 2 + 3 + O(3^2)]
            [  2*3 + O(3^2)         O(3^2)         O(3^2)    3^-1 + O(3)]
            [        O(3^2)         O(3^2)     3 + O(3^2)         O(3^2)]
        """
        from sage.schemes.hyperelliptic_curves import monsky_washnitzer

        if isinstance(p, (sage.rings.abc.pAdicField, sage.rings.abc.pAdicRing)):
            K = p
        else:
            K = pAdicField(p, prec)
        frob_p, _ = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(
            self.change_ring(K)
        )
        return frob_p

    def euler_factors(self, X, processes=None):
        r"""
        Return Euler factors at good primes up to ``X`` using ``hypellfrob``.

        The result is a list of pairs ``(p, L_p(T))`` in increasing order of
        ``p``, where ``L_p(T) = det(1 - Frob_p*T)``.  The minimal precision
        implied by the Weil bound is used separately for each prime.  The
        computations are distributed among ``processes`` independent
        processes.  If ``processes`` is ``None``, the multiprocessing default
        is used.

        This method currently requires an integral monic odd-degree model of
        the form ``y^2 = f(x)``.  Primes dividing the polynomial discriminant,
        as well as 2, are omitted.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: C = HyperellipticCurve(x^5 + 2*x^2 + x + 1)
            sage: [factor for p, factor in C.euler_factors(110, processes=2) if p == 37][0]
            1369*T^4 - 37*T^3 + 22*T^2 - T + 1

        The order of the result is independent of the number of processes.
        """
        if X < 2:
            return []

        f, h = self.hyperelliptic_polynomials()
        if h != 0:
            raise NotImplementedError("only implemented for y^2 = f(x)")
        if f.degree() < 3 or f.degree() % 2 == 0 or not f.is_monic():
            raise NotImplementedError(
                "requires a monic polynomial of odd degree at least 3"
            )
        if any(c not in ZZ for c in f.list()):
            raise ValueError("the defining polynomial must be integral")

        discriminant = ZZ(f.discriminant())
        coefficients = tuple(f.list())
        genus = self.genus()
        jobs = []
        for p in prime_range(3, X + 1):
            if discriminant % p:
                precision = _hypellfrob_precision(p, genus)
                if p <= f.degree() * (2 * precision - 1):
                    precision = None
                jobs.append((p, precision, coefficients, genus))
        with Pool(processes=processes) as pool:
            return pool.map(_hypellfrob_euler_factor, jobs)

    def lseries(self, prec=53):
        r"""
        Return the L-series of this hyperelliptic curve of genus 2.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: C = HyperellipticCurve(x^2+x, x^3+x^2+1)
            sage: C.lseries()
            PARI L-function associated to Hyperelliptic Curve
            over Rational Field defined by y^2 + (x^3 + x^2 + 1)*y = x^2 + x
        """
        from sage.lfunctions.pari import LFunction, lfun_genus2

        L = LFunction(lfun_genus2(self), prec=prec)
        L.rename("PARI L-function associated to %s" % self)
        return L
