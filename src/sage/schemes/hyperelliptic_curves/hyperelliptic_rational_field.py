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

import sage.rings.abc
from sage.rings.integer import Integer
from sage.rings.padics.factory import Qp as pAdicField
from sage.schemes.hyperelliptic_curves import hyperelliptic_generic


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
