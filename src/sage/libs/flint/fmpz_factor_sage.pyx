# sage_setup: distribution = sagemath-flint
from cysignals.signals cimport sig_check
from sage.libs.flint.fmpz cimport fmpz_get_mpz
from sage.rings.integer cimport Integer

cdef fmpz_factor_to_pairlist(const fmpz_factor_t factors):
    r"""
    Helper function that converts a fmpz_factor_t into a list of
    (factor, exponent) pairs. The factors are Integers, and the
    exponents are Python ints. This is used and indirectly tested by
    both :func:`qsieve` and
    :func:`sage.rings.factorint_flint.factor_using_flint`.  The output
    format was ultimately based on that of
    :func:`sage.rings.factorint_pari.factor_using_pari`.
    """
    cdef list pairs = []

    if factors.sign < 0:
        # FLINT doesn't return the plus/minus one factor.
        pairs.append((Integer(-1), int(1)))

    for i in range(factors.num):
        f = Integer()
        fmpz_get_mpz(f.value, &factors.p[i])
        e = int(factors.exp[i])
        pairs.append((f, e))
        sig_check()

    return pairs
