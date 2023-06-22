r"""
Univariate differential polynomials.

AUTHORS:

- Xavier Caruso, Raphaël Pagès (2023-06): initial version

"""

# ***************************************************************************
#    Copyright (C) 2023 Xavier Caruso <xavier.caruso@normalesup.org>
#                  2023 Raphaël Pagès <raphael.pages@math.u-bordeaux.fr>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.rings.polynomial.ore_polynomial_element import OrePolynomial_generic_dense
from sage.matrix.constructor import matrix


class DifferentialPolynomial_generic_dense(OrePolynomial_generic_dense):
    r"""
    Generic implementation of dense differential polynomial supporting
    any valid base ring and derivation.

    TESTS::

    sage: R.<t> = ZZ[]
    sage: S.<D>=OrePolynomialRing(R, t*R.derivation())
    sage: type(D)
    <class 'sage.rings.polynomial.differential_polynomial_element.DifferentialPolynomial_generic_dense'>

    sage: D*t
    t*D + t

    sage: category(D)
    Category of elements of Ore Polynomial Ring in D over Univariate Polynomial Ring in t over Integer Ring twisted by t*d/dt
    """

    def p_curvature(self, algorithm=None):
        r"""
        Return the `p`-curvature of this differential polynomial.

        INPUT:

        - ``algorithm`` -- a string or ``None`` (default: ``None``);
          the algorithm to use to compute the p-curvature, currenly
          only Katz's algorithm is available

        EXAMPLES::

            sage: A.<t> = GF(5)[]
            sage: S.<d> = A['d', A.derivation()]
            sage: d.p_curvature()
            [0]
            sage: (d^2).p_curvature()
            [0 0]
            [0 0]

            sage: L = d^3 + t*d
            sage: M = L.p_curvature()
            sage: M
            [        0         0         0]
            [      t^2       4*t 4*t^3 + 4]
            [        3       t^2         t]

        We verify that the coefficients of characteristic polynomial of
        the `p`-curvature are polynomials in `t^p`::

            sage: M.charpoly()
            x^3 + t^5*x

        When the base ring has characteristic zero, the `p`-curvature is
        not defined and an error is raised::

            sage: A.<t> = QQ[]
            sage: S.<d> = A['d', A.derivation()]
            sage: d.p_curvature()
            Traceback (most recent call last):
            ...
            ValueError: p-curvature only makes sense in positive characteristic

        TESTS::

            sage: A.<t> = GF(5)[]
            sage: S.<d> = A['d', A.derivation()]
            sage: d.p_curvature(algorithm="fast")
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'fast' is not available

        """
        p = self.parent().characteristic()
        if p == 0:
            raise ValueError("p-curvature only makes sense in positive characteristic")
        if algorithm is None:
            algorithm = "katz"
        methodname = "_p_curvature_" + algorithm
        if hasattr(self, methodname):
            method = getattr(self, "_p_curvature_%s" % algorithm)
            return method()
        raise ValueError("algorithm '%s' is not available" % algorithm)

    def _p_curvature_katz(self):
        r"""
        Return the `p`-curvature of this differential polynomial
        computed using Katz' algorithm.

        TESTS::

            sage: A.<t> = GF(5)[]
            sage: K = A.fraction_field()
            sage: S.<d> = K['d', K.derivation()]
            sage: L = d^2 + t^2*d - 1/t
            sage: L.p_curvature(algorithm='katz')  # indirect doctest
            [                       (4*t^9 + 4*t^6 + 3*t^4 + 3*t^3 + t + 4)/t^4 (t^12 + 3*t^9 + 3*t^7 + 3*t^6 + 2*t^4 + 2*t^3 + t^2 + 3*t + 4)/t^5]
            [                        (t^11 + 3*t^8 + 3*t^6 + 2*t^3 + t + 1)/t^3                 (4*t^14 + t^9 + t^6 + 2*t^4 + 2*t^3 + 4*t + 1)/t^4]

        ::

            sage: S.<x> = K['x', t*K.derivation()]
            sage: x.p_curvature(algorithm='katz')  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of the p-curvature is only implemented when d^p = 0

        """
        KD = self.parent()
        d = KD.twisting_derivation()
        p = KD.characteristic()
        if d.pth_power() != 0:
            raise NotImplementedError("computation of the p-curvature is only implemented when d^p = 0")
        A = self.companion_matrix()
        B = A
        for _ in range(p-1):
            B = A*B + B.apply_morphism(d)
        return B

