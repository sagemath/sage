# sage_setup: distribution = sagemath-mpmath
"""
Utility functions called by mpmath

Upstream mpmath uses these functions when it detects Sage.

Our vendored copy of mpmath, :mod:`sage.libs.mpmath._vendor.mpmath`,
uses these functions unconditionally,
see ``SAGE_ROOT/pkgs/sagemath-mpmath/pyproject.toml``.
"""

from sage.ext.stdsage cimport PY_NEW

from sage.rings.integer cimport Integer

from sage.libs.mpfr cimport *
from sage.libs.gmp.all cimport *

cpdef int bitcount(n) noexcept:
    """
    Bitcount of a Sage Integer or Python int/long.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import bitcount
        sage: bitcount(0)
        0
        sage: bitcount(1)
        1
        sage: bitcount(100)
        7
        sage: bitcount(-100)
        7
        sage: bitcount(2r)
        2
        sage: bitcount(2L)
        2
    """
    cdef Integer m
    if isinstance(n, Integer):
        m = <Integer>n
    else:
        m = Integer(n)
    if mpz_sgn(m.value) == 0:
        return 0
    return mpz_sizeinbase(m.value, 2)

cpdef isqrt(n):
    """
    Square root (rounded to floor) of a Sage Integer or Python int/long.
    The result is a Sage Integer.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import isqrt
        sage: isqrt(0)
        0
        sage: isqrt(100)
        10
        sage: isqrt(10)
        3
        sage: isqrt(10r)
        3
        sage: isqrt(10L)
        3
    """
    cdef Integer m, y
    if isinstance(n, Integer):
        m = <Integer>n
    else:
        m = Integer(n)
    if mpz_sgn(m.value) < 0:
        raise ValueError("square root of negative integer not defined.")
    y = PY_NEW(Integer)
    mpz_sqrt(y.value, m.value)
    return y

cpdef from_man_exp(man, exp, long prec = 0, str rnd = 'd'):
    """
    Create normalized mpf value tuple from mantissa and exponent.

    With prec > 0, rounds the result in the desired direction
    if necessary.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import from_man_exp
        sage: from_man_exp(-6, -1)
        (1, 3, 0, 2)
        sage: from_man_exp(-6, -1, 1, 'd')
        (1, 1, 1, 1)
        sage: from_man_exp(-6, -1, 1, 'u')
        (1, 1, 2, 1)
    """
    cdef Integer res
    cdef long bc
    res = Integer(man)
    bc = mpz_sizeinbase(res.value, 2)
    if not prec:
        prec = bc
    if mpz_sgn(res.value) < 0:
        mpz_neg(res.value, res.value)
        return normalize(1, res, exp, bc, prec, rnd)
    else:
        return normalize(0, res, exp, bc, prec, rnd)

cpdef normalize(long sign, Integer man, exp, long bc, long prec, str rnd):
    """
    Create normalized mpf value tuple from full list of components.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import normalize
        sage: normalize(0, 4, 5, 3, 53, 'n')
        (0, 1, 7, 1)
    """
    cdef long shift
    cdef Integer res
    cdef unsigned long trail
    if mpz_sgn(man.value) == 0:
        from sage.libs.mpmath._vendor.mpmath.libmp import fzero
        return fzero
    if bc <= prec and mpz_odd_p(man.value):
        return (sign, man, exp, bc)
    shift = bc - prec
    res = PY_NEW(Integer)
    if shift > 0:
        if rnd == 'n':
            if mpz_tstbit(man.value, shift-1) and (mpz_tstbit(man.value, shift)
                or (mpz_scan1(man.value, 0) < (shift-1))):
                mpz_cdiv_q_2exp(res.value, man.value, shift)
            else:
                mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'd':
            mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'f':
            if sign: mpz_cdiv_q_2exp(res.value, man.value, shift)
            else:    mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'c':
            if sign: mpz_fdiv_q_2exp(res.value, man.value, shift)
            else:    mpz_cdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'u':
            mpz_cdiv_q_2exp(res.value, man.value, shift)
        exp += shift
    else:
        mpz_set(res.value, man.value)
    # Strip trailing bits
    trail = mpz_scan1(res.value, 0)
    if 0 < trail < bc:
        mpz_tdiv_q_2exp(res.value, res.value, trail)
        exp += trail
    bc = mpz_sizeinbase(res.value, 2)
    return (sign, res, int(exp), bc)
