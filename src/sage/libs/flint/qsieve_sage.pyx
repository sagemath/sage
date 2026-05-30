"""
Interface to FLINT's ``qsieve_factor()``. This used to interact
with an external "QuadraticSieve" program, but its functionality has
been absorbed into flint.
"""

from cysignals.signals cimport sig_on, sig_off
from .types cimport fmpz_t, fmpz_factor_t
from .fmpz cimport fmpz_init, fmpz_set_mpz
from .fmpz_factor cimport fmpz_factor_init, fmpz_factor_clear
from .fmpz_factor_sage cimport fmpz_factor_to_pairlist
from .qsieve cimport qsieve_factor
from sage.rings.integer cimport Integer


def qsieve(n):
    r"""
    Factor ``n`` using the quadratic sieve.

    INPUT:

    - ``n`` -- integer; neither prime nor a perfect power

    OUTPUT:

    A list of the factors of ``n``. There is no guarantee that the
    factors found will be prime, or distinct.

    EXAMPLES::

        sage: k = 19; n = next_prime(10^k)*next_prime(10^(k+1))
        sage: factor(n)  # (currently) uses PARI
        10000000000000000051 * 100000000000000000039
        sage: qsieve(n)
        [(10000000000000000051, 1), (100000000000000000039, 1)]

    TESTS:

    The factorization of zero is undefined, to match the behavior of
    ``ZZ.zero().factor()``::

        sage: qsieve(ZZ.zero())
        Traceback (most recent call last):
        ...
        ArithmeticError: factorization of 0 is not defined
    """
    n = Integer(n)

    if n.is_zero():
        raise ArithmeticError("factorization of 0 is not defined")

    cdef fmpz_t p
    fmpz_init(p)
    fmpz_set_mpz(p, (<Integer>n).value)

    cdef fmpz_factor_t factors
    fmpz_factor_init(factors)
    sig_on()
    qsieve_factor(factors, p)
    sig_off()

    pairs = fmpz_factor_to_pairlist(factors)

    fmpz_factor_clear(factors)
    return pairs
