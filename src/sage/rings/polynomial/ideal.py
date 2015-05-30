"""
Ideals in Univariate Polynomial Rings.

AUTHORS:

- David Roe (2009-12-14) -- initial version.
"""

#*****************************************************************************
#
#   Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2009 DavidRoe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
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

from sage.rings.ideal import Ideal_pid

class Ideal_1poly_field(Ideal_pid):
    """
    An ideal in a univariate polynomial ring over a field.
    """
    def residue_class_degree(self):
        """
        Returns the degree of the generator of this ideal.

        This function is included for compatibility with ideals in rings of integers of number fields.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: P = R.ideal(t^4 + t + 1)
            sage: P.residue_class_degree()
            4
        """
        return self.gen().degree()

    def residue_field(self, names=None, check=True):
        """
        If this ideal is `P \subset F_p[t]`, returns the quotient `F_p[t]/P`.

        EXAMPLES::

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + 2*t + 9)
            sage: k.<a> = P.residue_field(); k
            Residue field in a of Principal ideal (t^3 + 2*t + 9) of Univariate Polynomial Ring in t over Finite Field of size 17
        """
        if check:
            if not self.ring().base_ring().is_finite():
                raise TypeError("residue fields only supported for polynomial rings over finite fields.")
            if not self.is_prime():
                raise ValueError("%s is not a prime ideal"%self)

        from sage.rings.finite_rings.residue_field import ResidueField
        return ResidueField(self, names, check=False)

