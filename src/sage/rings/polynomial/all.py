"""
Polynomials
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import

from sage.rings.polynomial.all__sagemath_categories import *

# Generic convolution
from sage.rings.polynomial.convolution import convolution

# Boolean Polynomial Rings
from sage.rings.polynomial.polynomial_ring_constructor import BooleanPolynomialRing_constructor as BooleanPolynomialRing

# Laurent Polynomial Rings
lazy_import('sage.rings.polynomial.omega', 'MacMahonOmega')

# Integer-valued Univariate Polynomial Ring
lazy_import('sage.rings.polynomial.integer_valued_polynomials',
            'IntegerValuedPolynomialRing')

# Ore Polynomial Rings
lazy_import('sage.rings.polynomial.ore_polynomial_ring', 'OrePolynomialRing')
SkewPolynomialRing = OrePolynomialRing

del lazy_import
