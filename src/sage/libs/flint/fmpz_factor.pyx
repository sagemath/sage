from cysignals.signals cimport sig_check
from sage.libs.flint.fmpz cimport fmpz_get_mpz
from sage.libs.gmp.mpz cimport mpz_t, mpz_init, mpz_set, mpz_clear
from sage.rings.integer cimport Integer

cdef fmpz_factor_to_pairlist(const fmpz_factor_t factors):
    r"""
    Helper function that converts a fmpz_factor_t into a list of
    (factor, exponent) pairs. The factors are Integers, and the
    exponents are Python ints. This is used and indirectly tested by
    both :func:`qsieve` and
    :func:`sage.rings.factorint.factor_using_flint`.  The output
    format was ultimately based on that of
    :func:`sage.rings.factorint.factor_using_pari`.
    """
    cdef list pairs = []

    if factors.sign < 0:
        # FLINT doesn't return the plus/minus one factor.
        pairs.append( (Integer(-1), int(1)) )

    cdef mpz_t mpz_factor
    for i in range(factors.num):
        mpz_init(mpz_factor)
        fmpz_get_mpz(mpz_factor, &factors.p[i])
        f = Integer()
        mpz_set(f.value, mpz_factor)
        mpz_clear(mpz_factor)
        e = int(factors.exp[i])
        pairs.append( (f,e) )
        sig_check()

    return pairs
