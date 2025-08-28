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

from sage.rings.infinity import Infinity

from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule


def CarlitzModule(A, base=None):
    r"""
    Return the Carlitz module over `A`.

    INPUT:

    - ``A`` -- a polynomial ring over a finite field

    - ``base`` -- a field, an element in a field or a
      string (default: the fraction field of ``A``)

    EXAMPLES::

        sage: Fq = GF(7)
        sage: A.<T> = Fq[]
        sage: CarlitzModule(A)
        Drinfeld module defined by T |--> τ + T

    We can specify a different base.
    This is interesting for instance for having two different variable
    names::

        sage: R.<z> = Fq[]
        sage: CarlitzModule(A, R)
        Drinfeld module defined by T |--> τ + z

    One can even use the following shortcut, which avoids the
    construction of `R`::

        sage: CarlitzModule(A, 'z')
        Drinfeld module defined by T |--> τ + z

    Using a similar syntax, we can construct the reduction of the
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
    elif isinstance(base, str):
        K = A.base_ring()[base]
        z = K.gen()
    else:
        raise ValueError("cannot construct a Carlitz module from the given data")
    return DrinfeldModule(A, [z, 1])


def carlitz_exponential(A, prec=+Infinity, name='z'):
    r"""
    Return the Carlitz exponential attached the ring `A`.

    INPUT:

    - ``prec`` -- an integer or ``Infinity`` (default: ``Infinity``);
      the precision at which the series is returned; if ``Infinity``,
      a lazy power series in returned, else, a classical power series
      is returned.

    - ``name`` -- string (default: ``'z'``); the name of the
      generator of the lazy power series ring

    EXAMPLES::

        sage: A.<T> = GF(2)[]

     When ``prec`` is ``Infinity`` (which is the default),
     the exponential is returned as a lazy power series, meaning
     that any of its coefficients can be computed on demands::

        sage: exp = carlitz_exponential(A)
        sage: exp
        z + ((1/(T^2+T))*z^2) + ((1/(T^8+T^6+T^5+T^3))*z^4) + O(z^8)
        sage: exp[2^4]
        1/(T^64 + T^56 + T^52 + ... + T^27 + T^23 + T^15)
        sage: exp[2^5]
        1/(T^160 + T^144 + T^136 + ... + T^55 + T^47 + T^31)

    On the contrary, when ``prec`` is a finite number, all the
    required coefficients are computed at once::

        sage: carlitz_exponential(A, prec=10)
        z + (1/(T^2 + T))*z^2 + (1/(T^8 + T^6 + T^5 + T^3))*z^4 + (1/(T^24 + T^20 + T^18 + T^17 + T^14 + T^13 + T^11 + T^7))*z^8 + O(z^10)

    We check that the Carlitz exponential is the compositional inverse
    of the Carlitz logarithm::

        sage: log = carlitz_logarithm(A)
        sage: exp(log)
        z + O(z^8)
        sage: log(exp)
        z + O(z^8)
    """
    C = CarlitzModule(A)
    return C.exponential(prec, name)


def carlitz_logarithm(A, prec=+Infinity, name='z'):
    r"""
    Return the Carlitz exponential attached the ring `A`.

    INPUT:

    - ``prec`` -- an integer or ``Infinity`` (default: ``Infinity``);
      the precision at which the series is returned; if ``Infinity``,
      a lazy power series in returned, else, a classical power series
      is returned.

    - ``name`` -- string (default: ``'z'``); the name of the
      generator of the lazy power series ring

    EXAMPLES::

        sage: A.<T> = GF(2)[]

     When ``prec`` is ``Infinity`` (which is the default),
     the exponential is returned as a lazy power series, meaning
     that any of its coefficients can be computed on demands::

        sage: log = carlitz_logarithm(A)
        sage: log
        z + ((1/(T^2+T))*z^2) + ((1/(T^6+T^5+T^3+T^2))*z^4) + O(z^8)
        sage: log[2^4]
        1/(T^30 + T^29 + T^27 + ... + T^7 + T^5 + T^4)
        sage: log[2^5]
        1/(T^62 + T^61 + T^59 + ... + T^8 + T^6 + T^5)

    On the contrary, when ``prec`` is a finite number, all the
    required coefficients are computed at once::

        sage: carlitz_logarithm(A, prec=10)
        z + (1/(T^2 + T))*z^2 + (1/(T^6 + T^5 + T^3 + T^2))*z^4 + (1/(T^14 + T^13 + T^11 + T^10 + T^7 + T^6 + T^4 + T^3))*z^8 + O(z^10)
    """
    C = CarlitzModule(A)
    return C.logarithm(prec, name)
