r"""
Projective plane conics over prime finite fields

AUTHORS:

- Marco Streng (2010-07-20)


"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2010 Marco Streng <marco.streng@gmail.com>
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

from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_prime_finite_field
from .con_finite_field import ProjectiveConic_finite_field

class ProjectiveConic_prime_finite_field(ProjectiveConic_finite_field, ProjectivePlaneCurve_prime_finite_field):
    r"""
    Create a projective plane conic curve over a prime finite field.
    See ``Conic`` for full documentation.

    EXAMPLES::

        sage: P.<X, Y, Z> = FiniteField(5)[]
        sage: Conic(X^2 + Y^2 - 2*Z^2)
        Projective Conic Curve over Finite Field of size 5 defined by X^2 + Y^2 - 2*Z^2

    TESTS::

        sage: Conic([FiniteField(7)(1), 1, -1])._test_pickling()
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        EXAMPLES ::

            sage: Conic([GF(3)(1), 1, 1])
            Projective Conic Curve over Finite Field of size 3 defined by x^2 + y^2 + z^2
        """
        ProjectiveConic_finite_field.__init__(self, A, f)
        ProjectivePlaneCurve_prime_finite_field.__init__(self, A, f)


