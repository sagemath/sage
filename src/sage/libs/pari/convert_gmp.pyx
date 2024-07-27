# sage_setup: distribution = sagemath-pari

"""
Convert PARI objects to/from GMP objects

Utility function to convert PARI ``GEN``s to/from the GMP types
``mpz_t`` and ``mpq_t``.

AUTHORS:

- Luca De Feo (2016-09-06): Separate Sage-specific components from
  generic C-interface in ``Pari`` (:issue:`20241`)
"""

#*****************************************************************************
#       Copyright (C) 2016 Luca De Feo <luca.defeo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on

from sage.libs.gmp.all cimport *

from cypari2.paridecl cimport *
from cypari2.stack cimport new_gen

cdef Gen new_gen_from_mpz_t(mpz_t value):
    """
    Create a new PARI Gen of type ``t_INT`` from a given
    GMP integer ``value``.

    EXAMPLES::

        sage: pari(42)       # indirect doctest
        42

    TESTS:

    Check that the hash of an integer does not depend on existing
    garbage on the stack (:issue:`11611`)::

        sage: foo = pari(2^(32*1024))  # Create large integer to put PARI stack in known state
        sage: a5 = pari(5)
        sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1))  # Dirty PARI stack
        sage: b5 = pari(5)
        sage: a5.__hash__() == b5.__hash__()
        True
    """
    sig_on()
    return new_gen(_new_GEN_from_mpz_t(value))


cdef inline GEN _new_GEN_from_mpz_t(mpz_t value) noexcept:
    r"""
    Create a new PARI ``t_INT`` from a ``mpz_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    if mpz_sgn(value) == 0:
        return gen_0

    cdef unsigned long limbs = mpz_size(value)

    cdef GEN z = cgeti(limbs + 2)
    # Set sign and "effective length"
    z[1] = evalsigne(mpz_sgn(value)) + evallgefint(limbs + 2)
    mpz_export(int_LSW(z), NULL, -1, sizeof(long), 0, 0, value)

    return z


cdef Gen new_gen_from_mpq_t(mpq_t value):
    """
    Create a new PARI Gen of type ``t_INT`` or ``t_FRAC`` from a given
    GMP rational ``value``.

    EXAMPLES::

        sage: pari(-2/3)
        -2/3
        sage: pari(QQ(42))
        42
        sage: pari(QQ(42)).type()
        't_INT'
        sage: pari(QQ(1/42)).type()
        't_FRAC'

    TESTS:

    Check that the hash of a rational does not depend on existing
    garbage on the stack (:issue:`11854`)::

        sage: foo = pari(2^(32*1024))  # Create large integer to put PARI stack in known state
        sage: a5 = pari(5/7)
        sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1))  # Dirty PARI stack
        sage: b5 = pari(5/7)
        sage: a5.__hash__() == b5.__hash__()
        True
    """
    sig_on()
    return new_gen(_new_GEN_from_mpq_t(value))


cdef inline GEN _new_GEN_from_mpq_t(mpq_t value) noexcept:
    r"""
    Create a new PARI ``t_INT`` or ``t_FRAC`` from a ``mpq_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN num = _new_GEN_from_mpz_t(mpq_numref(value))
    if mpz_cmpabs_ui(mpq_denref(value), 1) == 0:
        # Denominator is 1, return the numerator (an integer)
        return num
    cdef GEN denom = _new_GEN_from_mpz_t(mpq_denref(value))
    return mkfrac(num, denom)


cdef Gen new_gen_from_padic(long ordp, long relprec,
                            mpz_t prime, mpz_t p_pow, mpz_t unit):
    """
    Create a new PARI Gen of type ``t_PADIC`` from the given input data
    as GMP integers.
    """
    cdef GEN z
    sig_on()
    z = cgetg(5, t_PADIC)
    z[1] = evalprecp(relprec) + evalvalp(ordp)
    set_gel(z, 2, _new_GEN_from_mpz_t(prime))
    set_gel(z, 3, _new_GEN_from_mpz_t(p_pow))
    set_gel(z, 4, _new_GEN_from_mpz_t(unit))
    return new_gen(z)


cdef GEN _new_GEN_from_mpq_t_matrix(mpq_t** B, long nr, long nc) noexcept:
    """
    Create a new PARI ``t_MAT`` from a given
    2-dimensional array of GMP rationals ``mpq_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN x
    # Allocate zero matrix
    cdef GEN A = zeromatcopy(nr, nc)
    cdef long i, j
    for i in range(nr):
        for j in range(nc):
            x = _new_GEN_from_mpq_t(B[i][j])
            set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
    return A


cdef Gen rational_matrix(mpq_t** B, long nr, long nc):
    """
    Create a new PARI matrix of type ``t_MAT`` from a given
    array of GMP rationals ``mpq_t``.

    INPUT:

    - ``B`` -- a 2-dimensional array of ``mpq_t`` values. This array is
      accessed as ``B[i][j]``, where ``i`` is the row index and ``j``
      the column index.

    - ``nr`` -- number of rows of this matrix

    - ``nc`` -- number of columns of this matrix

    EXAMPLES::

        sage: matrix(QQ,2,[1..6]).__pari__()   # indirect doctest                       # needs sage.modules
        [1, 2, 3; 4, 5, 6]
    """
    sig_on()
    cdef GEN g = _new_GEN_from_mpq_t_matrix(B, nr, nc)
    return new_gen(g)


cdef inline void INT_to_mpz(mpz_ptr value, GEN g) noexcept:
    """
    Convert a PARI ``t_INT`` to a GMP integer, stored in ``value``.
    """
    if typ(g) != t_INT:
        pari_err(e_TYPE, <char*>"conversion to mpz", g)

    cdef long size = lgefint(g) - 2
    mpz_import(value, size, -1, sizeof(long), 0, 0, int_LSW(g))

    if signe(g) < 0:
        mpz_neg(value, value)


cdef void INTFRAC_to_mpq(mpq_ptr value, GEN g) noexcept:
    """
    Convert a PARI ``t_INT`` or ``t_FRAC`` to a GMP rational, stored in
    ``value``.
    """
    if typ(g) == t_FRAC:
        INT_to_mpz(mpq_numref(value), gel(g, 1))
        INT_to_mpz(mpq_denref(value), gel(g, 2))
    elif typ(g) == t_INT:
        INT_to_mpz(mpq_numref(value), g)
        mpz_set_ui(mpq_denref(value), 1)
    else:
        pari_err(e_TYPE, <char*>"conversion to mpq", g)
