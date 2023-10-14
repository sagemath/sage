# distutils: depends = flint/fmpz_poly_mat.h flint/fmpq_mpoly.h flint/nmod_poly_mat.h flint/perm.h flint/nmod_vec.h flint/fmpzi.h flint/mpoly.h flint/mpfr_mat.h flint/nmod_poly.h flint/long_extras.h flint/partitions.h flint/fq_nmod_poly_factor.h flint/mpf_mat.h flint/fmpz_poly_q.h flint/ca_vec.h flint/double_extras.h flint/fq_nmod_mpoly_factor.h flint/fmpz_mpoly.h flint/fmpq_vec.h flint/fq_poly.h flint/fq_nmod_mat.h flint/calcium.h flint/fmpz_vec.h flint/fexpr_builtin.h flint/ca_poly.h flint/fq_embed.h flint/nmod_poly_factor.h flint/fmpz_extras.h flint/acb_mat.h flint/gr_generic.h flint/gr_poly.h flint/fq_zech_poly_factor.h flint/fmpz_mpoly_q.h flint/nmod.h flint/ca_mat.h flint/fq_zech_vec.h flint/arb.h flint/nf_elem.h flint/fmpz_poly.h flint/fq_nmod_embed.h flint/arb_poly.h flint/fmpq.h flint/fmpq_mpoly_factor.h flint/acb_calc.h flint/fq_vec.h flint/padic_mat.h flint/fmpz.h flint/fmpz_mat.h flint/fmpz_mod_mat.h flint/fq_zech.h flint/double_interval.h flint/fq_default_mat.h flint/fmpz_mod_mpoly.h flint/qsieve.h flint/qfb.h flint/thread_pool.h flint/fmpz_mod_mpoly_factor.h flint/acb_modular.h flint/fq_nmod_poly.h flint/profiler.h flint/acb_hypgeom.h flint/d_mat.h flint/fq_zech_poly.h flint/fmpz_mod.h flint/qqbar.h flint/hypgeom.h flint/acb_elliptic.h flint/acb_dft.h flint/d_vec.h flint/ulong_extras.h flint/dlog.h flint/bool_mat.h flint/fmpq_poly.h flint/fexpr.h flint/mag.h flint/fmpz_mpoly_factor.h flint/mpfr_vec.h flint/ca_ext.h flint/gr.h flint/acb_poly.h flint/fft.h flint/padic.h flint/dirichlet.h flint/fmpz_mod_poly_factor.h flint/fq_default.h flint/fq.h flint/gr_vec.h flint/fq_nmod.h flint/fq_zech_embed.h flint/nmod_mpoly_factor.h flint/flint.h flint/fmpz_factor.h flint/qadic.h flint/nf.h flint/gr_mpoly.h flint/arb_fpwrap.h flint/fq_default_poly.h flint/aprcl.h flint/acf.h flint/gr_special.h flint/arf.h flint/fq_zech_mat.h flint/gr_mat.h flint/arb_fmpz_poly.h flint/fq_nmod_vec.h flint/fmpz_mod_poly.h flint/bernoulli.h flint/fq_default_poly_factor.h flint/arb_mat.h flint/arb_hypgeom.h flint/fq_poly_factor.h flint/nmod_mpoly.h flint/acb.h flint/fmpz_lll.h flint/fmpz_poly_factor.h flint/ca_field.h flint/acb_dirichlet.h flint/arith.h flint/fmpz_mod_vec.h flint/fmpq_mat.h flint/fq_mat.h flint/mpn_extras.h flint/mpf_vec.h flint/ca.h flint/padic_poly.h flint/nmod_mat.h flint/fq_nmod_mpoly.h flint/arb_calc.h flint/nmod_types.h

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

# Use these typedefs in lieu of flint's ulong and slong macros
ctypedef mp_limb_t ulong
ctypedef mp_limb_signed_t slong
ctypedef mp_limb_t flint_bitcnt_t


cdef extern from "flint_wrap.h":
    # flint/fmpz.h
    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]

    bint COEFF_IS_MPZ(fmpz)
    mpz_ptr COEFF_TO_PTR(fmpz)

    ctypedef struct fmpz_preinvn_struct:
        mp_ptr dinv
        long n
        mp_bitcnt_t norm

    ctypedef fmpz_preinvn_struct[1] fmpz_preinvn_t

    ctypedef struct fmpz_comb_struct:
        pass

    ctypedef fmpz_comb_struct fmpz_comb_t[1]

    ctypedef struct fmpz_comb_temp_struct:
        slong Alen, Tlen
        fmpz * A, * T

    ctypedef fmpz_comb_temp_struct fmpz_comb_temp_t[1]

    ctypedef struct fmpz_multi_CRT_struct:
        pass

    ctypedef fmpz_multi_CRT_struct fmpz_multi_CRT_t[1]


    # flint/fmpq.h
    ctypedef struct fmpq:
        pass

    ctypedef fmpq fmpq_t[1]


    # flint/arith.h
    ctypedef struct trig_prod_struct:
        pass
    ctypedef trig_prod_struct trig_prod_t[1]


    # flint/ulong_extras.pxd
    ctypedef struct n_ecm_s:
        pass
    ctypedef n_ecm_s n_ecm_t[1];


    # flint/arf.h
    ctypedef struct arf_struct:
        pass
    ctypedef arf_struct arf_t[1]
    ctypedef arf_struct * arf_ptr
    ctypedef const arf_struct * arf_srcptr
    cdef enum arf_rnd_t:
        ARF_RND_DOWN
        ARF_RND_UP
        ARF_RND_FLOOR
        ARF_RND_CEIL
        ARF_RND_NEAR
    long ARF_PREC_EXACT


    # flint/arb_types.h
    ctypedef struct mag_struct:
        pass
    ctypedef mag_struct mag_t[1]
    ctypedef mag_struct * mag_ptr
    ctypedef const mag_struct * mag_srcptr

    ctypedef struct arb_struct:
        pass
    ctypedef arb_struct arb_t[1]
    ctypedef arb_struct * arb_ptr
    ctypedef const arb_struct * arb_srcptr

    ctypedef struct arb_mat_struct:
        arb_ptr entries
        slong r
        slong c
        arb_ptr * rows
    ctypedef arb_mat_struct arb_mat_t[1]

    ctypedef struct arb_poly_struct:
        pass
    ctypedef arb_poly_struct[1] arb_poly_t
    ctypedef arb_poly_struct * arb_poly_ptr
    ctypedef const arb_poly_struct * arb_poly_srcptr


    # flint/arb_calc.h
    ctypedef int (*arb_calc_func_t)(arb_ptr out, const arb_t inp,
                                    void * param, slong order, slong prec)


    # flint/acb.h
    ctypedef struct acb_struct:
        pass
    ctypedef acb_struct[1] acb_t
    ctypedef acb_struct * acb_ptr
    ctypedef const acb_struct * acb_srcptr



    # flint/acb_mat.h
    ctypedef struct acb_mat_struct:
        pass
    ctypedef acb_mat_struct[1] acb_mat_t


    # flint/acb_modular.h
    ctypedef struct psl2z_struct:
        fmpz a
        fmpz b
        fmpz c
        fmpz d
    ctypedef psl2z_struct psl2z_t[1];


    # flint/acb_poly.h
    ctypedef struct acb_poly_struct:
        pass
    ctypedef acb_poly_struct[1] acb_poly_t
    ctypedef acb_poly_struct * acb_poly_ptr
    ctypedef const acb_poly_struct * acb_poly_srcptr

    # flint/acb_calc.h
    ctypedef struct acb_calc_integrate_opt_struct:
        long deg_limit
        long eval_limit
        long depth_limit
        bint use_heap
        int verbose
    ctypedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1]
    ctypedef int (*acb_calc_func_t)(acb_ptr out,
            const acb_t inp, void * param, long order, long prec)

    # flint/d_mat.h
    ctypedef struct d_mat_struct:
        double * entries
        slong r
        slong c
        double ** rows

    ctypedef d_mat_struct d_mat_t[1]


    # flint/flint.h
    ctypedef struct flint_rand_s:
        pass
    ctypedef flint_rand_s flint_rand_t[1]

    ctypedef void* flint_rand_t
    cdef long FLINT_BITS
    cdef long FLINT_D_BITS


    # flint/limb_types.h
    ctypedef struct n_factor_t:
        int num
        unsigned long exp[15]
        unsigned long p[15]

    ctypedef struct n_primes_struct:
        pass

    ctypedef n_primes_struct n_primes_t[1]

    long FLINT_MAX_FACTORS_IN_LIMB


    # flint/fmpz_factor.h
    ctypedef struct ecm_s:
        pass

    ctypedef ecm_s ecm_t[1]


    # flint/fmpz_mod.h
    ctypedef struct fmpz_mod_discrete_log_pohlig_hellman_entry_struct:
        pass

    ctypedef struct fmpz_mod_discrete_log_pohlig_hellman_struct:
        pass

    ctypedef fmpz_mod_discrete_log_pohlig_hellman_struct fmpz_mod_discrete_log_pohlig_hellman_t[1]


    # flint/fmpz_types.h
    ctypedef struct fmpz_factor_struct:
        int sign
        fmpz* p
        ulong* exp
        slong alloc
        slong num

    ctypedef fmpz_factor_struct fmpz_factor_t[1]

    ctypedef struct fmpz_poly_struct:
        fmpz* coeffs
        long alloc
        long length

    ctypedef fmpz_poly_struct fmpz_poly_t[1]

    ctypedef struct fmpz_poly_factor_struct:
        pass

    ctypedef fmpz_poly_factor_struct fmpz_poly_factor_t[1]

    ctypedef struct fmpz_mat_struct:
        pass

    ctypedef fmpz_mat_struct fmpz_mat_t[1]

    ctypedef struct fmpz_poly_mat_struct:
        pass

    ctypedef fmpz_poly_mat_struct fmpz_poly_mat_t[1]

    ctypedef struct fmpz_mpoly_struct:
        pass

    ctypedef fmpz_mpoly_struct fmpz_mpoly_t[1]

    ctypedef struct fmpz_mpoly_factor_struct:
        pass

    ctypedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1];

    ctypedef struct fmpz_poly_q_struct:
        fmpz_poly_struct *num
        fmpz_poly_struct *den

    ctypedef fmpz_poly_q_struct fmpz_poly_q_t[1]

    ctypedef struct fmpz_mpoly_q_struct:
        pass

    ctypedef fmpz_mpoly_q_struct fmpz_mpoly_q_t[1]

    ctypedef struct fmpzi_struct:
        fmpz a
        fmpz b

    ctypedef fmpzi_struct fmpzi_t[1]


    # flint/fmpz_poly.h
    ctypedef struct fmpz_poly_powers_precomp_struct:
        fmpz ** powers
        slong len

    ctypedef fmpz_poly_powers_precomp_struct fmpz_poly_powers_precomp_t[1]

    ctypedef struct fmpz_poly_mul_precache_struct:
       mp_limb_t ** jj
       slong n
       slong len2
       slong loglen
       slong bits2
       slong limbs
       fmpz_poly_t poly2

    ctypedef fmpz_poly_mul_precache_struct fmpz_poly_mul_precache_t[1]


    # flint/fmpz_mod_types.h
    ctypedef struct fmpz_mod_ctx_struct:
        pass

    ctypedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1]

    ctypedef struct fmpz_mod_mat_struct:
        fmpz_mat_t mat
        fmpz_t mod

    ctypedef fmpz_mod_mat_struct fmpz_mod_mat_t[1]

    ctypedef struct fmpz_mod_poly_struct:
        fmpz * coeffs
        slong alloc
        slong length

    ctypedef fmpz_mod_poly_struct fmpz_mod_poly_t[1]

    ctypedef struct fmpz_mod_poly_factor_struct:
        fmpz_mod_poly_struct * poly
        slong *exp
        slong num
        slong alloc

    ctypedef fmpz_mod_poly_factor_struct fmpz_mod_poly_factor_t[1]

    ctypedef struct fmpz_mod_mpoly_struct:
        fmpz * coeffs
        ulong * exps
        slong length
        flint_bitcnt_t bits
        slong coeffs_alloc
        slong exps_alloc

    ctypedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1]

    ctypedef struct fmpz_mod_mpoly_factor_struct:
        fmpz_t constant
        fmpz_mod_mpoly_struct * poly
        fmpz * exp
        slong num
        slong alloc

    ctypedef fmpz_mod_mpoly_factor_struct fmpz_mod_mpoly_factor_t[1]


    # flint/fmpq_poly.h
    ctypedef struct fmpq_poly_struct:
        pass

    ctypedef fmpq_poly_struct fmpq_poly_t[1]

    ctypedef struct fmpq_poly_powers_precomp_struct:
        fmpq_poly_struct * powers
        slong len

    ctypedef fmpq_poly_powers_precomp_struct fmpq_poly_powers_precomp_t[1]


    # flint/fmpq_mat.h
    ctypedef struct fmpq_mat_struct:
        pass

    ctypedef fmpq_mat_struct fmpq_mat_t[1]


    # flint/fmpz_mod_poly.h:
    ctypedef struct fmpz_mod_poly_res_struct:
        pass

    ctypedef fmpz_mod_poly_res_struct fmpz_mod_poly_res_t[1]

    ctypedef struct fmpz_mod_poly_frobenius_powers_2exp_struct:
        pass

    ctypedef fmpz_mod_poly_frobenius_powers_2exp_struct fmpz_mod_poly_frobenius_powers_2exp_t[1]

    ctypedef struct fmpz_mod_poly_frobenius_powers_struct:
        pass

    ctypedef fmpz_mod_poly_frobenius_powers_struct fmpz_mod_poly_frobenius_powers_t[1]

    ctypedef struct fmpz_mod_poly_matrix_precompute_arg_t:
        pass

    ctypedef struct fmpz_mod_poly_compose_mod_precomp_preinv_arg_t:
        pass

    ctypedef struct fmpz_mod_poly_radix_struct:
        pass

    ctypedef fmpz_mod_poly_radix_struct fmpz_mod_poly_radix_t[1]

    ctypedef struct fmpz_mod_berlekamp_massey_struct:
        pass

    ctypedef fmpz_mod_berlekamp_massey_struct fmpz_mod_berlekamp_massey_t[1]


    # flint/nmod_poly.h
    ctypedef struct nmod_t:
        mp_limb_t n
        mp_limb_t ninv
        mp_bitcnt_t norm

    ctypedef struct nmod_poly_multi_crt_struct:
        pass

    ctypedef nmod_poly_multi_crt_struct nmod_poly_multi_crt_t[1]

    ctypedef struct nmod_berlekamp_massey_struct:
        pass

    ctypedef nmod_berlekamp_massey_struct nmod_berlekamp_massey_t[1]


    # flint/nmod_types.h
    ctypedef struct nmod_mat_struct:
        mp_limb_t * entries
        slong r
        slong c
        mp_limb_t ** rows
        nmod_t mod

    ctypedef nmod_mat_struct nmod_mat_t[1]

    ctypedef struct nmod_poly_struct:
        mp_ptr coeffs
        slong alloc
        slong length
        nmod_t mod

    ctypedef nmod_poly_struct nmod_poly_t[1]

    ctypedef struct nmod_poly_factor_struct:
        nmod_poly_struct * p
        slong *exp
        slong num
        slong alloc

    ctypedef nmod_poly_factor_struct nmod_poly_factor_t[1]

    ctypedef struct nmod_poly_mat_struct:
        nmod_poly_struct * entries
        slong r
        slong c
        nmod_poly_struct ** rows
        mp_limb_t modulus

    ctypedef nmod_poly_mat_struct nmod_poly_mat_t[1]

    ctypedef struct nmod_mpoly_struct:
        mp_limb_t * coeffs
        ulong * exps
        slong length
        flint_bitcnt_t bits
        slong coeffs_alloc
        slong exps_alloc

    ctypedef nmod_mpoly_struct nmod_mpoly_t[1]

    ctypedef struct nmod_mpoly_factor_struct:
        mp_limb_t constant
        nmod_mpoly_struct * poly
        fmpz * exp
        slong num
        slong alloc

    ctypedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1]


    # flint/fq.h
    ctypedef struct fq_ctx_struct:
        fmpz_mod_poly_t modulus

    ctypedef fq_ctx_struct fq_ctx_t[1]

    ctypedef fmpz_poly_struct fq_struct
    ctypedef fmpz_poly_t fq_t


    # flint/fq_nmod.h
    ctypedef struct fq_nmod_ctx_struct:
        nmod_poly_t modulus

    ctypedef fq_nmod_ctx_struct fq_nmod_ctx_t[1]

    ctypedef nmod_poly_struct fq_nmod_struct
    ctypedef nmod_poly_t fq_nmod_t


    # flint/padic.h
    ctypedef struct padic_struct:
        fmpz u
        long v

    ctypedef padic_struct padic_t[1]

    cdef enum padic_print_mode:
        PADIC_TERSE
        PADIC_SERIES
        PADIC_VAL_UNIT

    ctypedef struct padic_ctx_struct:
        fmpz_t p
        long N
        double pinv
        fmpz* pow
        long min
        long max

    ctypedef padic_ctx_struct padic_ctx_t[1]

    ctypedef struct padic_inv_struct:
        long n
        fmpz *pow
        fmpz *u

    ctypedef padic_inv_struct padic_inv_t[1]


    # flint/padic_poly.h
    ctypedef struct padic_poly_struct:
        fmpz *coeffs
        long alloc
        long length
        long val
        long N

    ctypedef padic_poly_struct padic_poly_t[1]


    # flint/qadic.h
    ctypedef struct qadic_ctx_struct:
        padic_ctx_struct pctx
        fmpz *a
        long *j
        long len
        char *var

    ctypedef qadic_ctx_struct qadic_ctx_t[1]

    ctypedef padic_poly_struct qadic_struct
    ctypedef padic_poly_t qadic_t


    # flint/qsieve.h
    ctypedef struct prime_t:
        pass

    ctypedef struct fac_t:
        pass

    ctypedef struct la_col_t:
        pass

    ctypedef struct hash_t:
        pass

    ctypedef struct relation_t:
        pass

    ctypedef struct qs_poly_s:
        pass
    ctypedef qs_poly_s qs_poly_t[1]

    ctypedef struct qs_s:
        pass
    ctypedef qs_s qs_t[1];


    # flint/thread_pool.h
    ctypedef struct thread_pool_entry_struct:
        pass

    ctypedef thread_pool_entry_struct thread_pool_entry_t[1]

    ctypedef struct thread_pool_struct:
        pass

    ctypedef thread_pool_struct thread_pool_t[1]
    ctypedef int thread_pool_handle


    # flint/bernoulli.h
    ctypedef struct bernoulli_rev_struct:
        pass
    ctypedef bernoulli_rev_struct bernoulli_rev_t[1]


    # flint/nf.h
    ctypedef struct nf_struct:
        pass
    ctypedef nf_struct nf_t[1]


    # flint/nf_elem.h
    ctypedef struct lnf_elem_struct:
        pass
    ctypedef lnf_elem_struct lnf_elem_t[1]

    ctypedef struct qnf_elem_struct:
        pass
    ctypedef qnf_elem_struct qnf_elem_t[1]

    ctypedef union nf_elem_struct:
        fmpq_poly_t elem
        lnf_elem_t lelem
        qnf_elem_t qelem
    ctypedef nf_elem_struct nf_elem_t[1]


    # flint/ca.h
    ctypedef union ca_elem_struct:
        fmpq q
        nf_elem_struct nf
        fmpz_mpoly_q_struct * mpoly_q

    ctypedef struct ca_struct:
        ulong field
        ca_elem_struct elem

    ctypedef ca_struct ca_t[1]
    ctypedef ca_struct * ca_ptr
    ctypedef const ca_struct * ca_srcptr

    ctypedef struct ca_ext_qqbar:
        pass

    ctypedef struct ca_ext_func_data:
        pass

    ctypedef struct ca_ext_struct:
        pass

    ctypedef ca_ext_struct ca_ext_t[1]
    ctypedef ca_ext_struct * ca_ext_ptr
    ctypedef const ca_ext_struct * ca_ext_srcptr

    ctypedef struct ca_ext_cache_struct:
        pass

    ctypedef ca_ext_cache_struct ca_ext_cache_t[1]

    ctypedef struct ca_field_struct:
        pass

    ctypedef ca_field_struct ca_field_t[1]
    ctypedef ca_field_struct * ca_field_ptr
    ctypedef const ca_field_struct * ca_field_srcptr

    ctypedef struct ca_field_cache_struct:
        pass

    ctypedef ca_field_cache_struct ca_field_cache_t[1];

    cdef enum:
        CA_OPT_VERBOSE
        CA_OPT_PRINT_FLAGS
        CA_OPT_MPOLY_ORD
        CA_OPT_PREC_LIMIT
        CA_OPT_QQBAR_DEG_LIMIT
        CA_OPT_LOW_PREC
        CA_OPT_SMOOTH_LIMIT
        CA_OPT_LLL_PREC
        CA_OPT_POW_LIMIT
        CA_OPT_USE_GROEBNER
        CA_OPT_GROEBNER_LENGTH_LIMIT
        CA_OPT_GROEBNER_POLY_LENGTH_LIMIT
        CA_OPT_GROEBNER_POLY_BITS_LIMIT
        CA_OPT_VIETA_LIMIT
        CA_OPT_TRIG_FORM
        CA_OPT_NUM_OPTIONS

    ctypedef struct ca_ctx_struct:
        pass

    ctypedef ca_ctx_struct ca_ctx_t[1]
