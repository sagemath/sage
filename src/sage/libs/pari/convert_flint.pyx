"""
Convert PARI objects to/from Flint objects

Utility function to convert PARI ``GEN``s to/from flint types.

AUTHORS:

- Luca De Feo (2016-09-06): Separate Sage-specific components from
  generic C-interface in ``PariInstance`` (:trac:`20241`)
"""

#*****************************************************************************
#       Copyright (C) 2016 Luca De Feo <luca.defeo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, division, print_function

include "cysignals/signals.pxi"

from sage.libs.flint.fmpz cimport fmpz_get_mpz, COEFF_IS_MPZ, COEFF_TO_PTR
from sage.libs.flint.fmpz_mat cimport *

from sage.libs.cypari2.paridecl cimport *
from sage.libs.cypari2.stack cimport new_gen
from .convert_gmp cimport _new_GEN_from_mpz_t


cdef inline GEN _new_GEN_from_fmpz_t(fmpz_t value):
    r"""
    Create a new PARI ``t_INT`` from a ``fmpz_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    if COEFF_IS_MPZ(value[0]):
        return _new_GEN_from_mpz_t(COEFF_TO_PTR(value[0]))
    else:
        return stoi(value[0])


cdef GEN _new_GEN_from_fmpz_mat_t(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc):
    r"""
    Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
    from a ``fmpz_mat_t``.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN x
    cdef GEN A = zeromatcopy(nr, nc)
    cdef Py_ssize_t i, j
    for i in range(nr):
        for j in range(nc):
            x = _new_GEN_from_fmpz_t(fmpz_mat_entry(B,i,j))
            set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
    return A


cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc):
    r"""
    Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
    from a ``fmpz_mat_t`` and rotate the matrix 90 degrees
    counterclockwise.  So the resulting matrix will have ``nc`` rows
    and ``nr`` columns.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN x
    cdef GEN A = zeromatcopy(nc, nr)
    cdef Py_ssize_t i, j
    for i in range(nr):
        for j in range(nc):
            x = _new_GEN_from_fmpz_t(fmpz_mat_entry(B,i,nc-j-1))
            set_gcoeff(A, j+1, i+1, x)  # A[j+1, i+1] = x (using 1-based indexing)
    return A


cdef gen integer_matrix(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf):
    """
    EXAMPLES::

        sage: matrix(ZZ,2,[1..6])._pari_()   # indirect doctest
        [1, 2, 3; 4, 5, 6]
    """
    sig_on()
    cdef GEN g
    if permute_for_hnf:
        g = _new_GEN_from_fmpz_mat_t_rotate90(B, nr, nc)
    else:
        g = _new_GEN_from_fmpz_mat_t(B, nr, nc)
    return new_gen(g)
