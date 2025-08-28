# sage.doctest: needs sage.libs.flint
r"""
Integer factorization using FLINT

AUTHORS:

- Michael Orlitzky (2023)
"""

#*****************************************************************************
#       Copyright (C) 2023 Michael Orlitzky
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_set_mpz
from sage.libs.flint.fmpz_factor cimport *
from sage.libs.flint.fmpz_factor_sage cimport *
from sage.rings.integer cimport Integer


def factor_using_flint(Integer n):
    r"""
    Factor the nonzero integer ``n`` using FLINT.

    This function returns a list of (factor, exponent) pairs. The
    factors will be of type ``Integer``, and the exponents will be of
    type ``int``.

    INPUT:

    - ``n`` -- a nonzero sage Integer; the number to factor

    OUTPUT:

    A list of ``(Integer, int)`` pairs representing the factors and
    their exponents.

    EXAMPLES::

        sage: from sage.rings.factorint_flint import factor_using_flint
        sage: n = ZZ(9962572652930382)
        sage: factors = factor_using_flint(n)
        sage: factors
        [(2, 1), (3, 1), (1660428775488397, 1)]
        sage: prod( f^e for (f,e) in factors ) == n
        True

     Negative numbers will have a leading factor of ``(-1)^1``::

        sage: n = ZZ(-1 * 2 * 3)
        sage: factor_using_flint(n)
        [(-1, 1), (2, 1), (3, 1)]

    The factorization of unity is empty::

        sage: factor_using_flint(ZZ.one())
        []

    While zero has a single factor, of... zero::

        sage: factor_using_flint(ZZ.zero())
        [(0, 1)]

    TESTS:

    Check that the integers [-10,000, 10,000] are factored correctly::

        sage: all(
        ....:   prod( f^e for (f,e) in factor_using_flint(ZZ(c*k)) ) == c*k
        ....:   for k in range(10000)
        ....:   for c in [-1, 1]
        ....: )
        True
    """
    if n.is_zero():
        return [(n, int(1))]

    cdef fmpz_t p
    fmpz_init(p)
    fmpz_set_mpz(p, (<Integer>n).value)

    cdef fmpz_factor_t factors
    fmpz_factor_init(factors)

    sig_on()
    fmpz_factor(factors, p)
    sig_off()

    pairs = fmpz_factor_to_pairlist(factors)

    fmpz_factor_clear(factors)
    return pairs
