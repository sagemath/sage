r"""
Multivariate Tropical Polynomial Semirings

<Description>

AUTHORS:

- Verrel Rievaldo Wijaya

EXAMPLES:

Construct multivariate tropical polynomial semirings::

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R = PolynomialRing(T, 'a,b')
    sage: R
    Multivarite Tropical Polynomial Semiring in a, b over Rational Field
    
Create an element by inputting a dictionary::

    sage: dict1 = {(1,0):0, (0,1):-1, (1,1):3}
    sage: p1 = R(dict1); p1
    3*a*b + 0*a + (-1)*b

We can also create an element by converting from classical polynomial::

    sage: S.<a,b> = QQ[]
    sage: f = a + b + a*b
    sage: p2 = R(f); p2
    1*a*b + 1*a + 1*b

Some basic arithmetic operations::

    sage: p1 + p2
    3*a*b + 1*a + 1*b
    sage: p1 * p2
    4*a^2*b^2 + 4*a^2*b + 1*a^2 + 4*a*b^2 + 1*a*b + 0*b^2
    sage: T(2) * p1
    5*a*b + 2*a + 1*b

TESTS:

    sage: -p1
    Traceback (most recent call last):
    ...
    ArithmeticError: cannot negate any non-infinite element

REFERENCES:

    - [Fil2017]_

"""

# ****************************************************************************
#       Copyright (C) 2024 Verrel Rievaldo Wijaya <verrelrievaldo@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

class TropicalMPolynomial(MPolynomial_polydict):
    
    def roots(self):
        pass

class TropicalMPolynomialSemiring(UniqueRepresentation, Parent):
    def __init__(self, base_semiring, names):
        Parent.__init__(self, base=base_semiring, names=names)

    Element = TropicalMPolynomial

    def _element_constructor_(self, x):
        C = self.element_class
        new_dict = {}

        if isinstance(x, MPolynomial):
            x = x.dict()

        for key, value in x.items():
            new_dict[key] = self.base()(value)

        return C(self, new_dict)
    
    def _repr_(self):
        return (f"Multivarite Tropical Polynomial Semiring in {', '.join(self.variable_names())}"
            f" over {self.base_ring().base_ring()}")