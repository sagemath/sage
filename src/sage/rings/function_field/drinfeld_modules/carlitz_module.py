r"""
Carlitz module

AUTHORS:

- Xavier Caruso (2025-07): initial version
"""

# *****************************************************************************
#        Copyright (C) 2025 Xavier Caruso <xavier@caruso.ovh>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.finite_fields import FiniteFields

from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule


def CarlitzModule(A, base=None):
    r"""
    Return the Carlitz module over `A`.

    INPUT:

    - ``A`` -- a polynomial ring over a finite field

    - ``base`` -- a field or an element in a field
      (default: the fraction field of ``A``)

    EXAMPLES::

        sage: Fq = GF(7)
        sage: A.<T> = Fq[]
        sage: CarlitzModule(A)
        Drinfeld module defined by T |--> τ + T

    We can specify a different base.
    This is interesting for instance for having two different variable
    names::

        sage: R.<z> = Fq[]
        sage: K = Frac(R)
        sage: CarlitzModule(A, K)
        Drinfeld module defined by T |--> τ + z

    Using a similar syntax, we can construct the reduction over the
    Carlitz module modulo primes::

        sage: F.<a> = Fq.extension(z^2 + 1)
        sage: CarlitzModule(A, F)
        Drinfeld module defined by T |--> τ + a

    It is also possible to pass in any element in the base field
    (in this case, the result might not be strictly speaking the
    Carlitz module, but it is always a Drinfeld module of rank 1)::

        sage: CarlitzModule(A, z^2)
        Drinfeld module defined by T |--> τ + z^2

    TESTS::

        sage: CarlitzModule(Fq)
        Traceback (most recent call last):
        ...
        TypeError: the function ring must be defined over a finite field

    ::

        sage: S.<x,y> = QQ[]
        sage: CarlitzModule(A, S)
        Traceback (most recent call last):
        ...
        ValueError: function ring base must coerce into base field
    """
    if (not isinstance(A, PolynomialRing_generic)
     or A.base_ring() not in FiniteFields()):
        raise TypeError('the function ring must be defined over a finite field')
    if base is None:
        K = A.fraction_field()
        z = K.gen()
    elif isinstance(base, Parent):
        if base.has_coerce_map_from(A):
            z = base(A.gen())
        else:
            z = base.gen()
    elif isinstance(base, Element):
        z = base
    else:
        raise ValueError("cannot construct a Carlitz module from the given data")
    return DrinfeldModule(A, [z, 1])
