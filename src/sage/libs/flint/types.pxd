"""
Declarations for FLINT types
"""

#*****************************************************************************
#       Copyright (C) 2014 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.types cimport *

cdef extern from "flint/flint.h":
    ctypedef mp_limb_t ulong
    ctypedef mp_limb_signed_t slong

cdef extern from "flint/fmpq.h":
    ctypedef struct fmpq:
        pass

    ctypedef fmpq fmpq_t[1]

cdef extern from "flint/fmpq_poly.h":
    ctypedef struct fmpq_poly_struct:
        pass

    ctypedef fmpq_poly_struct fmpq_poly_t[1]

cdef extern from "flint/fmpz.h":
    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]

    ctypedef struct fmpz_preinvn_struct:
        mp_ptr dinv
        long n
        mp_bitcnt_t norm

    ctypedef fmpz_preinvn_struct[1] fmpz_preinvn_t

cdef extern from "flint/fmpz_mat.h":
    ctypedef struct fmpz_mat_struct:
        pass

    ctypedef fmpz_mat_struct fmpz_mat_t[1]

cdef extern from "flint/fmpz_poly.h":
    ctypedef struct fmpz_poly_struct:
        pass

    ctypedef fmpz_poly_struct fmpz_poly_t[1]

cdef extern from "flint/nmod_poly.h":
    ctypedef struct nmod_t:
        mp_limb_t n
        mp_limb_t ninv
        mp_bitcnt_t norm

    ctypedef struct nmod_poly_struct:
        mp_limb_t *coeffs
        long alloc
        long length
        nmod_t mod

    ctypedef nmod_poly_struct nmod_poly_t[1]

    ctypedef struct nmod_poly_factor_struct:
        nmod_poly_t p
        long *exp
        long num
        long alloc

    ctypedef nmod_poly_factor_struct nmod_poly_factor_t[1]

cdef extern from "flint/ulong_extras.h":
    ctypedef struct n_factor_t:
        int num
        unsigned long exp[15]
        unsigned long p[15]
