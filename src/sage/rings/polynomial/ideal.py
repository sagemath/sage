"""
Ideals in Univariate Polynomial Rings

AUTHORS:

- David Roe (2009-12-14) -- initial version.
"""

#*****************************************************************************
#       Copyright (C) 2009 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ideal import Ideal_pid

class Ideal_1poly_field(Ideal_pid):
    """
    An ideal in a univariate polynomial ring over a field.
    """
    def residue_class_degree(self):
        """
        Return the degree of the generator of this ideal.

        This function is included for compatibility with ideals in rings of integers of number fields.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: P = R.ideal(t^4 + t + 1)
            sage: P.residue_class_degree()
            4
        """
        return self.gen().degree()

    def residue_field(self, names=None, check=True):
        r"""
        If this ideal is `P \subset F_p[t]`, return the quotient `F_p[t]/P`.

        EXAMPLES::

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + 2*t + 9)
            sage: k.<a> = P.residue_field(); k                                          # needs sage.rings.finite_rings
            Residue field in a of Principal ideal (t^3 + 2*t + 9) of
             Univariate Polynomial Ring in t over Finite Field of size 17
        """
        if check:
            if not self.ring().base_ring().is_finite():
                raise TypeError("residue fields only supported for polynomial rings over finite fields.")
            if not self.is_prime():
                raise ValueError("%s is not a prime ideal" % self)

        from sage.rings.finite_rings.residue_field import ResidueField
        return ResidueField(self, names, check=False)

    def groebner_basis(self, algorithm=None):
        """
        Return a Gröbner basis for this ideal.

        The Gröbner basis has 1 element, namely the generator of the
        ideal. This trivial method exists for compatibility with
        multi-variate polynomial rings.

        INPUT:

        - ``algorithm`` -- ignored

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: I = R.ideal([x^2 - 1, x^3 - 1])
            sage: G = I.groebner_basis(); G
            [x - 1]
            sage: type(G)
            <class 'sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic'>
            sage: list(G)
            [x - 1]
        """
        gb = self.gens_reduced()
        from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence_generic
        return PolynomialSequence_generic([gb], self.ring(), immutable=True)

    def change_ring(self, R):
        """
        Coerce an ideal into a new ring.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: I = R.ideal([q^2 + q - 1])
            sage: I.change_ring(RR['q'])                                                # needs sage.rings.real_mpfr
            Principal ideal (q^2 + q - 1.00000000000000) of
             Univariate Polynomial Ring in q over Real Field with 53 bits of precision
        """
        return R.ideal(self.gens())
