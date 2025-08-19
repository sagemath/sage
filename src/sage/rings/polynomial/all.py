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

# Quotient of polynomial ring
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing
from sage.rings.polynomial.polynomial_quotient_ring_element import PolynomialQuotientRingElement

# Univariate Polynomial Rings
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import polygen, polygens
from sage.rings.polynomial.polynomial_element import Polynomial

# Multivariate Polynomial Rings
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.multi_polynomial_element import degree_lowest_rational_function

# Generic convolution
from sage.rings.polynomial.convolution import convolution

# Boolean Polynomial Rings
from sage.rings.polynomial.polynomial_ring_constructor import BooleanPolynomialRing_constructor as BooleanPolynomialRing

# Laurent Polynomial Rings
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
lazy_import('sage.rings.polynomial.omega', 'MacMahonOmega')

# Infinite Polynomial Rings
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing

# Ore Polynomial Rings
lazy_import('sage.rings.polynomial.ore_polynomial_ring', 'OrePolynomialRing')
SkewPolynomialRing = OrePolynomialRing

# Evaluation of cyclotomic polynomials
from sage.rings.polynomial.cyclotomic import cyclotomic_value

# Integer-valued Univariate Polynomial Ring
lazy_import('sage.rings.polynomial.integer_valued_polynomials',
            'IntegerValuedPolynomialRing')
lazy_import('sage.rings.polynomial.q_integer_valued_polynomials',
            'QuantumValuedPolynomialRing')
