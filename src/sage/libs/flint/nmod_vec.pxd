# distutils: libraries = flint
# distutils: depends = flint/nmod_vec.h

# This file was (manually) generated from FLINT's nmod_vec.h.
#*****************************************************************************
#       Copyright (C) 2010 William Hart
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.flint.types cimport *

# flint/nmod_vec.h
cdef extern from "flint_wrap.h":
    mp_limb_t _nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t _nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
    mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)
    void nmod_init(nmod_t * mod, mp_limb_t n)
    mp_ptr _nmod_vec_init(long len)
    void _nmod_vec_clear(mp_ptr vec)
    void _nmod_vec_randtest(mp_ptr vec, flint_rand_t state, long len, nmod_t mod)
    void _nmod_vec_zero(mp_ptr vec, long len)
    mp_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, long len)
    void _nmod_vec_set(mp_ptr res, mp_srcptr vec, long len)
    void _nmod_vec_swap(mp_ptr a, mp_ptr b, long length)
    int _nmod_vec_equal(mp_ptr vec, mp_srcptr vec2, long len)
    int _nmod_vec_is_zero(mp_srcptr vec, long len)
    void _nmod_vec_reduce(mp_ptr res, mp_srcptr vec, long len, nmod_t mod)
    void _nmod_vec_add(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, long len, nmod_t mod)
    void _nmod_vec_sub(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, long len, nmod_t mod)
    void _nmod_vec_neg(mp_ptr res, mp_srcptr vec, long len, nmod_t mod)
    void _nmod_vec_scalar_mul_nmod(mp_ptr res, mp_srcptr vec, long len, mp_limb_t c, nmod_t mod)
    void _nmod_vec_scalar_addmul_nmod(mp_ptr res, mp_srcptr vec, long len, mp_limb_t c, nmod_t mod)
    int _nmod_vec_dot_bound_limbs(long len, nmod_t mod)
    mp_limb_t _nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2, long len, nmod_t mod, int nlimbs)
    mp_limb_t _nmod_vec_dot_ptr(mp_srcptr vec1, mp_ptr * vec2, long offset, long len, nmod_t mod, int nlimbs)
