"""
Convert PARI objects to/from Flint objects

Utility function to convert PARI ``GEN``s to/from flint types.

AUTHORS:

- Luca De Feo (2016-09-06): Separate Sage-specific components from
  generic C-interface in ``Pari`` (:issue:`20241`)
"""

# ****************************************************************************
#       Copyright (C) 2016 Luca De Feo <luca.defeo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on

from sage.libs.flint.fmpz cimport COEFF_IS_MPZ, COEFF_TO_PTR, fmpz_is_one
from sage.libs.flint.fmpq cimport fmpq_numref, fmpq_denref
from sage.libs.flint.fmpz_mat cimport fmpz_mat_nrows, fmpz_mat_ncols, fmpz_mat_entry
from sage.libs.flint.fmpq_mat cimport fmpq_mat_nrows, fmpq_mat_ncols, fmpq_mat_entry

from cypari2.paridecl cimport *
from cypari2.stack cimport new_gen
from sage.libs.pari.convert_gmp cimport _new_GEN_from_mpz_t


cdef inline GEN _new_GEN_from_fmpz_t(fmpz_t value) noexcept:
    r"""
    Create a new PARI ``t_INT`` from a ``fmpz_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    if COEFF_IS_MPZ(value[0]):
        return _new_GEN_from_mpz_t(COEFF_TO_PTR(value[0]))
    else:
        return stoi(value[0])


cdef inline GEN _new_GEN_from_fmpq_t(fmpq_t value) noexcept:
    r"""
    Create a new PARI ``t_RAT`` from a ``fmpq_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN num = _new_GEN_from_fmpz_t(fmpq_numref(value))
    if fmpz_is_one(fmpq_denref(value)):
        return num
    cdef GEN denom = _new_GEN_from_fmpz_t(fmpq_denref(value))
    return mkfrac(num, denom)


cdef GEN _new_GEN_from_fmpz_mat_t(fmpz_mat_t B) noexcept:
    r"""
    Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
    from a ``fmpz_mat_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN x
    cdef GEN A = zeromatcopy(fmpz_mat_nrows(B), fmpz_mat_ncols(B))
    cdef Py_ssize_t i, j
    for i in range(fmpz_mat_nrows(B)):
        for j in range(fmpz_mat_ncols(B)):
            x = _new_GEN_from_fmpz_t(fmpz_mat_entry(B, i, j))
            set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
    return A


cdef GEN _new_GEN_from_fmpq_mat_t(fmpq_mat_t B) noexcept:
    cdef GEN x
    cdef GEN A = zeromatcopy(fmpq_mat_nrows(B), fmpq_mat_ncols(B))
    cdef Py_ssize_t i, j
    for i in range(fmpq_mat_nrows(B)):
        for j in range(fmpq_mat_ncols(B)):
            x = _new_GEN_from_fmpq_t(fmpq_mat_entry(B, i, j))
            set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
    return A

cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(fmpz_mat_t B) noexcept:
    r"""
    Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
    from a ``fmpz_mat_t`` and rotate the matrix 90 degrees
    counterclockwise.  So the resulting matrix will have ``nc`` rows
    and ``nr`` columns.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN x
    cdef GEN A = zeromatcopy(fmpz_mat_ncols(B), fmpz_mat_nrows(B))
    cdef Py_ssize_t i, j
    for i in range(fmpz_mat_nrows(B)):
        for j in range(fmpz_mat_ncols(B)):
            x = _new_GEN_from_fmpz_t(fmpz_mat_entry(B, i, fmpz_mat_ncols(B) - j - 1))
            set_gcoeff(A, j+1, i+1, x)  # A[j+1, i+1] = x (using 1-based indexing)
    return A


cdef GEN _new_GEN_from_fmpq_mat_t_rotate90(fmpq_mat_t B) noexcept:
    r"""
    Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
    from a ``fmpq_mat_t`` and rotate the matrix 90 degrees
    counterclockwise.  So the resulting matrix will have ``nc`` rows
    and ``nr`` columns.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN x
    cdef GEN A = zeromatcopy(fmpq_mat_ncols(B), fmpq_mat_nrows(B))
    cdef Py_ssize_t i, j
    for i in range(fmpq_mat_nrows(B)):
        for j in range(fmpq_mat_ncols(B)):
            x = _new_GEN_from_fmpq_t(fmpq_mat_entry(B, i, fmpq_mat_ncols(B) - j - 1))
            set_gcoeff(A, j+1, i+1, x)  # A[j+1, i+1] = x (using 1-based indexing)
    return A


cdef Gen integer_matrix(fmpz_mat_t B, bint rotate):
    """
    EXAMPLES::

        sage: matrix(ZZ,2,[1..6]).__pari__()   # indirect doctest
        [1, 2, 3; 4, 5, 6]
    """
    cdef GEN g
    sig_on()
    if rotate:
        g = _new_GEN_from_fmpz_mat_t_rotate90(B)
    else:
        g = _new_GEN_from_fmpz_mat_t(B)
    return new_gen(g)


cdef Gen rational_matrix(fmpq_mat_t B, bint rotate):
    """
    EXAMPLES::

        sage: matrix(QQ, 2, [1/2, 2/3, 3/4, 4/5]).__pari__() # indirect doctest
        [1/2, 2/3; 3/4, 4/5]
    """
    cdef GEN g
    sig_on()
    if rotate:
        g = _new_GEN_from_fmpq_mat_t_rotate90(B)
    else:
        g = _new_GEN_from_fmpq_mat_t(B)
    return new_gen(g)
