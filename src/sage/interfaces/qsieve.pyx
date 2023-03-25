"""
Interface to FLINT's ``qsieve_factor()``. This used to interact
with an external "QuadraticSieve" program, but its functionality has
been absorbed into flint.
"""

from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_set_mpz, fmpz_get_mpz
from sage.libs.flint.fmpz_factor cimport *
from sage.libs.flint.qsieve cimport *
from sage.libs.gmp.mpz cimport mpz_t, mpz_init, mpz_set, mpz_clear
from sage.rings.integer cimport Integer

def qsieve(n):
    r"""
    Factor ``n`` using the quadratic sieve.

    INPUT:

    - ``n`` -- an integer; neither prime nor a perfect power.

    OUTPUT:

    A list of the factors of ``n``. There is no guarantee that the
    factors found will be prime, or distinct.

    EXAMPLES::

        sage: k = 19; n = next_prime(10^k)*next_prime(10^(k+1))
        sage: factor(n)  # (currently) uses PARI
        10000000000000000051 * 100000000000000000039
        sage: qsieve(n)
        [(10000000000000000051, 1), (100000000000000000039, 1)]

    """
    if n.is_zero():
        return [n]

    n = Integer(n)

    sig_on()
    cdef fmpz_t p
    fmpz_init(p)
    fmpz_set_mpz(p, (<Integer>n).value)

    cdef fmpz_factor_t factors
    fmpz_factor_init(factors)
    qsieve_factor(factors,p)
    sig_off()

    pairs = []
    if factors.sign < 0:
        # FLINT doesn't return the plus/minus one factor.
        pairs.append( (Integer(-1), int(1)) )

    cdef mpz_t mpz_factor;
    for i in range(factors.num):
        mpz_init(mpz_factor)
        fmpz_get_mpz(mpz_factor, &factors.p[i])
        f = Integer()
        mpz_set(f.value, mpz_factor)
        mpz_clear(mpz_factor)
        e = int(factors.exp[i])
        pairs.append( (f,e) )
        sig_check()

    fmpz_factor_clear(factors)
    return pairs
