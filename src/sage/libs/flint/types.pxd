# distutils: depends = flint/acb.h flint/acb_calc.h flint/acb_dft.h flint/acb_dirichlet.h flint/acb_elliptic.h flint/acb_hypgeom.h flint/acb_mat.h flint/acb_modular.h flint/acb_poly.h flint/acf.h flint/aprcl.h flint/arb.h flint/arb_calc.h flint/arb_fmpz_poly.h flint/arb_fpwrap.h flint/arb_hypgeom.h flint/arb_mat.h flint/arb_poly.h flint/arf.h flint/arith.h flint/bernoulli.h flint/bool_mat.h flint/ca.h flint/ca_ext.h flint/ca_field.h flint/ca_mat.h flint/ca_poly.h flint/ca_vec.h flint/calcium.h flint/d_mat.h flint/d_vec.h flint/dirichlet.h flint/dlog.h flint/double_extras.h flint/double_interval.h flint/fexpr.h flint/fexpr_builtin.h flint/fft.h flint/flint.h flint/fmpq.h flint/fmpq_mat.h flint/fmpq_mpoly.h flint/fmpq_mpoly_factor.h flint/fmpq_poly.h flint/fmpq_vec.h flint/fmpz.h flint/fmpz_extras.h flint/fmpz_factor.h flint/fmpz_lll.h flint/fmpz_mat.h flint/fmpz_mod.h flint/fmpz_mod_mat.h flint/fmpz_mod_mpoly.h flint/fmpz_mod_mpoly_factor.h flint/fmpz_mod_poly.h flint/fmpz_mod_poly_factor.h flint/fmpz_mod_vec.h flint/fmpz_mpoly.h flint/fmpz_mpoly_factor.h flint/fmpz_mpoly_q.h flint/fmpz_poly.h flint/fmpz_poly_factor.h flint/fmpz_poly_mat.h flint/fmpz_poly_q.h flint/fmpz_vec.h flint/fmpzi.h flint/fq.h flint/fq_default.h flint/fq_default_mat.h flint/fq_default_poly.h flint/fq_default_poly_factor.h flint/fq_embed.h flint/fq_mat.h flint/fq_nmod.h flint/fq_nmod_embed.h flint/fq_nmod_mat.h flint/fq_nmod_mpoly.h flint/fq_nmod_mpoly_factor.h flint/fq_nmod_poly.h flint/fq_nmod_poly_factor.h flint/fq_nmod_vec.h flint/fq_poly.h flint/fq_poly_factor.h flint/fq_vec.h flint/fq_zech.h flint/fq_zech_embed.h flint/fq_zech_mat.h flint/fq_zech_poly.h flint/fq_zech_poly_factor.h flint/fq_zech_vec.h flint/gr.h flint/gr_generic.h flint/gr_mat.h flint/gr_mpoly.h flint/gr_poly.h flint/gr_special.h flint/gr_vec.h flint/hypgeom.h flint/long_extras.h flint/mag.h flint/mpf_mat.h flint/mpf_vec.h flint/mpfr_mat.h flint/mpfr_vec.h flint/mpn_extras.h flint/mpoly.h flint/nf.h flint/nf_elem.h flint/nmod.h flint/nmod_mat.h flint/nmod_mpoly.h flint/nmod_mpoly_factor.h flint/nmod_poly.h flint/nmod_poly_factor.h flint/nmod_poly_mat.h flint/nmod_types.h flint/nmod_vec.h flint/padic.h flint/padic_mat.h flint/padic_poly.h flint/partitions.h flint/perm.h flint/profiler.h flint/qadic.h flint/qfb.h flint/qqbar.h flint/qsieve.h flint/thread_pool.h flint/ulong_extras.h

# WARNING: src/sage/libs/flint/types.pxd is generated from
# src/sage_setup/autogen/flint/templates/types.pxd.template
# please make sure that you are modifying the correct file!
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
        fmpz * A
        fmpz * T

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
    ctypedef n_ecm_s n_ecm_t[1]


    # flint/arf.h
    ctypedef enum arf_rnd_t:
        ARF_RND_DOWN
        ARF_RND_UP
        ARF_RND_FLOOR
        ARF_RND_CEIL
        ARF_RND_NEAR
    long ARF_PREC_EXACT


    # flint/arf_types.h
    ctypedef struct mantissa_noptr_struct:
        pass

    ctypedef struct mantissa_ptr_struct:
        pass

    ctypedef union mantissa_struct:
        mantissa_noptr_struct noptr
        mantissa_ptr_struct ptr

    ctypedef struct arf_struct:
        pass
    ctypedef arf_struct arf_t[1]
    ctypedef arf_struct * arf_ptr
    ctypedef const arf_struct * arf_srcptr

    ctypedef struct arf_interval_struct:
        arf_struct a
        arf_struct b
    ctypedef arf_interval_struct arf_interval_t[1]
    ctypedef arf_interval_struct * arf_interval_ptr
    ctypedef const arf_interval_struct * arf_interval_srcptr


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
        arb_ptr coeffs
        long alloc
        long length
    ctypedef arb_poly_struct[1] arb_poly_t


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
    ctypedef psl2z_struct psl2z_t[1]


    # flint/acb_poly.h
    ctypedef struct acb_poly_struct:
        acb_ptr coeffs
        long alloc
        long length
    ctypedef acb_poly_struct[1] acb_poly_t

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


    # flint/acb_dft.h
    ctypedef struct crt_struct:
        pass
    ctypedef crt_struct crt_t[1]

    ctypedef struct acb_dft_step_struct:
        pass
    ctypedef acb_dft_step_struct * acb_dft_step_ptr

    ctypedef struct acb_dft_cyc_struct:
        pass
    ctypedef acb_dft_cyc_struct acb_dft_cyc_t[1]

    ctypedef struct acb_dft_rad2_struct:
        pass
    ctypedef acb_dft_rad2_struct acb_dft_rad2_t[1]

    ctypedef struct acb_dft_bluestein_struct:
        pass
    ctypedef acb_dft_bluestein_struct acb_dft_bluestein_t[1]

    ctypedef struct acb_dft_prod_struct:
        pass
    ctypedef acb_dft_prod_struct acb_dft_prod_t[1]

    ctypedef struct acb_dft_crt_struct:
        pass
    ctypedef acb_dft_crt_struct acb_dft_crt_t[1]

    ctypedef struct acb_dft_naive_struct:
        pass
    ctypedef acb_dft_naive_struct acb_dft_naive_t[1]

    ctypedef struct acb_dft_pre_struct:
        pass
    ctypedef acb_dft_pre_struct acb_dft_pre_t[1]


    # flint/acb_dirichlet.h
    ctypedef struct acb_dirichlet_hurwitz_precomp_struct:
        pass
    ctypedef acb_dirichlet_hurwitz_precomp_struct acb_dirichlet_hurwitz_precomp_t[1]

    ctypedef struct acb_dirichlet_roots_struct:
        pass
    ctypedef acb_dirichlet_roots_struct acb_dirichlet_roots_t[1]

    ctypedef struct acb_dirichlet_platt_c_precomp_struct:
        pass
    ctypedef acb_dirichlet_platt_c_precomp_struct acb_dirichlet_platt_c_precomp_t[1]

    ctypedef struct acb_dirichlet_platt_i_precomp_struct:
        pass
    ctypedef acb_dirichlet_platt_i_precomp_struct acb_dirichlet_platt_i_precomp_t[1]

    ctypedef struct acb_dirichlet_platt_ws_precomp_struct:
        pass
    ctypedef acb_dirichlet_platt_ws_precomp_struct acb_dirichlet_platt_ws_precomp_t[1]


    # flint/acb_theta.h
    cdef struct acb_theta_eld_struct:
        pass
    ctypedef acb_theta_eld_struct acb_theta_eld_t[1]

    ctypedef void (*acb_theta_naive_worker_t)(acb_ptr, acb_srcptr, acb_srcptr, const slong *,
        slong, const acb_t, const slong *, slong, slong, slong, slong)
    ctypedef int (*acb_theta_ql_worker_t)(acb_ptr, acb_srcptr, acb_srcptr,
        arb_srcptr, arb_srcptr, const acb_mat_t, slong, slong)

    # flint/d_mat.h
    ctypedef struct d_mat_struct:
        double * entries
        slong r
        slong c
        double ** rows

    ctypedef d_mat_struct d_mat_t[1]


    # flint/flint.h
    ctypedef struct flint_rand_struct:
        pass
    ctypedef flint_rand_struct flint_rand_t[1]

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

    ctypedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1]

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


    # flint/fmpq_types.h
    ctypedef struct fmpq_mat_struct:
        pass
    ctypedef fmpq_mat_struct fmpq_mat_t[1]

    ctypedef struct fmpq_poly_struct:
        pass
    ctypedef fmpq_poly_struct fmpq_poly_t[1]

    ctypedef struct fmpq_mpoly_struct:
        pass
    ctypedef fmpq_mpoly_struct fmpq_mpoly_t[1]

    ctypedef struct fmpq_mpoly_factor_struct:
        pass
    ctypedef fmpq_mpoly_factor_struct fmpq_mpoly_factor_t[1]


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


    # flint/nmod_mpoly.h
    ctypedef struct nmod_mpoly_univar_struct:
        pass
    ctypedef nmod_mpoly_univar_struct nmod_mpoly_univar_t[1]

    ctypedef struct nmod_mpolyu_struct:
        pass
    ctypedef nmod_mpolyu_struct nmod_mpolyu_t[1]

    ctypedef struct nmod_mpolyn_struct:
        pass
    ctypedef nmod_mpolyn_struct nmod_mpolyn_t[1]

    ctypedef struct nmod_mpolyun_struct:
        pass
    ctypedef nmod_mpolyun_struct nmod_mpolyun_t[1]

    ctypedef struct nmod_mpolyd_struct:
        pass
    ctypedef nmod_mpolyd_struct nmod_mpolyd_t[1]

    ctypedef struct nmod_poly_stack_struct:
        pass
    ctypedef nmod_poly_stack_struct nmod_poly_stack_t[1]


    # flint/fq_nmod_types.h
    ctypedef nmod_poly_t fq_nmod_t
    ctypedef nmod_poly_struct fq_nmod_struct

    ctypedef struct fq_nmod_ctx_struct:
        pass
    ctypedef fq_nmod_ctx_struct fq_nmod_ctx_t[1]

    ctypedef struct fq_nmod_mat_struct:
        pass
    ctypedef fq_nmod_mat_struct fq_nmod_mat_t[1]

    ctypedef struct fq_nmod_poly_struct:
        pass
    ctypedef fq_nmod_poly_struct fq_nmod_poly_t[1]

    ctypedef struct fq_nmod_poly_factor_struct:
        pass
    ctypedef fq_nmod_poly_factor_struct fq_nmod_poly_factor_t[1]


    # flint2/fq_nmod_mpoly.h
    ctypedef struct fq_nmod_mpoly_ctx_struct:
        pass
    ctypedef fq_nmod_mpoly_ctx_struct fq_nmod_mpoly_ctx_t[1]

    ctypedef struct fq_nmod_mpoly_struct:
        pass
    ctypedef fq_nmod_mpoly_struct fq_nmod_mpoly_t[1]

    ctypedef struct fq_nmod_mpoly_univar_struct:
        pass
    ctypedef fq_nmod_mpoly_univar_struct fq_nmod_mpoly_univar_t[1]

    ctypedef struct fq_nmod_mpolyu_struct:
        pass
    ctypedef fq_nmod_mpolyu_struct fq_nmod_mpolyu_t[1]

    ctypedef struct fq_nmod_mpolyn_struct:
        pass
    ctypedef fq_nmod_mpolyn_struct fq_nmod_mpolyn_t[1]

    ctypedef struct fq_nmod_mpolyun_struct:
        pass
    ctypedef fq_nmod_mpolyun_struct fq_nmod_mpolyun_t[1]

    ctypedef struct bad_fq_nmod_embed_struct:
        pass
    ctypedef bad_fq_nmod_embed_struct bad_fq_nmod_embed_t[1]


    # flint2/fq_nmod_mpoly_factor.h
    ctypedef struct fq_nmod_mpoly_factor_struct:
        pass
    ctypedef fq_nmod_mpoly_factor_struct fq_nmod_mpoly_factor_t[1]

    ctypedef struct fq_nmod_mpolyv_struct:
        pass
    ctypedef fq_nmod_mpolyv_struct fq_nmod_mpolyv_t[1]

    ctypedef struct fq_nmod_mpoly_pfrac_struct:
        pass
    ctypedef fq_nmod_mpoly_pfrac_struct fq_nmod_mpoly_pfrac_t[1]


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
    ctypedef struct fmpq_poly_powers_precomp_struct:
        fmpq_poly_struct * powers
        slong len

    ctypedef fmpq_poly_powers_precomp_struct fmpq_poly_powers_precomp_t[1]


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


    # flint/nmod.h
    ctypedef struct nmod_discrete_log_pohlig_hellman_table_entry_struct:
        pass

    ctypedef struct nmod_discrete_log_pohlig_hellman_entry_struct:
        pass

    ctypedef struct nmod_discrete_log_pohlig_hellman_struct:
        pass
    ctypedef nmod_discrete_log_pohlig_hellman_struct nmod_discrete_log_pohlig_hellman_t[1]


    # flint/fq_types.h
    ctypedef fmpz_poly_t fq_t
    ctypedef fmpz_poly_struct fq_struct

    ctypedef struct fq_ctx_struct:
        pass

    ctypedef fq_ctx_struct fq_ctx_t[1]

    ctypedef struct fq_mat_struct:
        pass
    ctypedef fq_mat_struct fq_mat_t[1]

    ctypedef struct fq_poly_struct:
        pass
    ctypedef fq_poly_struct fq_poly_t[1]

    ctypedef struct fq_poly_factor_struct:
        pass
    ctypedef fq_poly_factor_struct fq_poly_factor_t[1]


    # flint/fq_zech_types.h
    ctypedef struct fq_zech_struct:
        pass
    ctypedef fq_zech_struct fq_zech_t[1]

    ctypedef struct fq_zech_ctx_struct:
        pass
    ctypedef fq_zech_ctx_struct fq_zech_ctx_t[1]

    ctypedef struct fq_zech_mat_struct:
        pass
    ctypedef fq_zech_mat_struct fq_zech_mat_t[1]

    ctypedef struct fq_zech_poly_struct:
        pass
    ctypedef fq_zech_poly_struct fq_zech_poly_t[1]

    ctypedef struct fq_zech_poly_factor_struct:
        pass
    ctypedef fq_zech_poly_factor_struct fq_zech_poly_factor_t[1]


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


    # flint/fmpz_mod_mpoly.h
    ctypedef struct fmpz_mod_mpoly_univar_struct:
        pass
    ctypedef fmpz_mod_mpoly_univar_struct fmpz_mod_mpoly_univar_t[1]


    # flint/padic_poly.h
    ctypedef struct padic_poly_struct:
        fmpz *coeffs
        long alloc
        long length
        long val
        long N

    ctypedef padic_poly_struct padic_poly_t[1]


    # flint/padic_mat.h
    ctypedef struct padic_mat_struct:
        pass
    ctypedef padic_mat_struct padic_mat_t[1]


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
    ctypedef qs_s qs_t[1]


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


    # flint/mpoly_types.h
    ctypedef enum ordering_t:
        ORD_LEX
        ORD_DEGLEX
        ORD_DEGREVLEX

    ctypedef struct mpoly_ctx_struct:
        pass
    ctypedef mpoly_ctx_struct mpoly_ctx_t[1]

    ctypedef struct nmod_mpoly_ctx_struct:
        pass
    ctypedef nmod_mpoly_ctx_struct nmod_mpoly_ctx_t[1]

    ctypedef struct fmpz_mpoly_ctx_struct:
        pass
    ctypedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1]

    ctypedef struct fmpq_mpoly_ctx_struct:
        pass

    ctypedef fmpq_mpoly_ctx_struct fmpq_mpoly_ctx_t[1]

    ctypedef struct fmpz_mod_mpoly_ctx_struct:
        pass
    ctypedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1]


    # flint/mpoly.h
    ctypedef struct mpoly_heap_t:
        pass

    ctypedef struct mpoly_nheap_t:
        pass

    ctypedef struct mpoly_heap1_s:
        pass

    ctypedef struct mpoly_heap_s:
        pass

    ctypedef struct mpoly_rbnode_ui_struct:
        pass

    ctypedef struct mpoly_rbtree_ui_struct:
        pass

    ctypedef mpoly_rbtree_ui_struct mpoly_rbtree_ui_t[1]

    ctypedef struct mpoly_rbnode_fmpz_struct:
        pass

    ctypedef struct mpoly_rbtree_fmpz_struct:
        pass
    ctypedef mpoly_rbtree_fmpz_struct mpoly_rbtree_fmpz_t[1]

    ctypedef struct mpoly_gcd_info_struct:
        pass
    ctypedef mpoly_gcd_info_struct mpoly_gcd_info_t[1]

    ctypedef struct mpoly_compression_struct:
        pass
    ctypedef mpoly_compression_struct mpoly_compression_t[1]

    ctypedef struct mpoly_univar_struct:
        pass
    ctypedef mpoly_univar_struct mpoly_univar_t[1]

    ctypedef struct mpoly_void_ring_struct:
        pass
    ctypedef mpoly_void_ring_struct mpoly_void_ring_t[1]

    ctypedef struct string_with_length_struct:
        pass

    ctypedef struct mpoly_parse_struct:
        pass
    ctypedef mpoly_parse_struct mpoly_parse_t[1]


    # flint/fmpz_mpoly.h
    ctypedef struct fmpz_mpoly_univar_struct:
        pass
    ctypedef fmpz_mpoly_univar_struct fmpz_mpoly_univar_t[1]

    ctypedef struct fmpz_mpolyd_struct:
        pass
    ctypedef fmpz_mpolyd_struct fmpz_mpolyd_t[1]

    ctypedef struct fmpz_mpoly_vec_struct:
        pass
    ctypedef fmpz_mpoly_vec_struct fmpz_mpoly_vec_t[1]

    ctypedef struct fmpz_mpolyd_ctx_struct:
        pass
    ctypedef fmpz_mpolyd_ctx_struct fmpz_mpolyd_ctx_t[1]

    ctypedef struct fmpz_pow_cache_struct:
        pass
    ctypedef fmpz_pow_cache_struct fmpz_pow_cache_t[1]

    ctypedef struct fmpz_mpoly_geobucket_struct:
        pass
    ctypedef fmpz_mpoly_geobucket_struct fmpz_mpoly_geobucket_t[1]


    # flint/fmpq_mpoly.h
    ctypedef struct fmpq_mpoly_univar_struct:
        pass
    ctypedef fmpq_mpoly_univar_struct fmpq_mpoly_univar_t[1]


    # flint/fft_small.h
    ctypedef struct mpn_ctx_struct:
        pass
    ctypedef mpn_ctx_struct mpn_ctx_t[1]

    ctypedef struct mul_precomp_struct:
        pass

    ctypedef struct nmod_poly_divrem_precomp_struct:
        pass


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


    # flint/machine_vectors.h
    ctypedef struct vec1n:
        pass
    ctypedef struct vec2n:
        pass
    ctypedef struct vec4n:
        pass
    ctypedef struct vec8n:
        pass

    ctypedef struct vec1d:
        pass
    ctypedef struct vec2d:
        pass
    ctypedef struct vec4d:
        pass
    ctypedef struct vec8d:
        pass


    # flint/calcium.h
    ctypedef struct calcium_stream_struct:
        pass
    ctypedef calcium_stream_struct calcium_stream_t[1]


    ctypedef enum truth_t:
        T_TRUE
        T_FALSE
        T_UNKNOWN

    ctypedef enum calcium_func_code:
        CA_QQBar
        CA_Neg
        CA_Add
        CA_Sub
        CA_Mul
        CA_Div
        CA_Sqrt
        CA_Cbrt
        CA_Root
        CA_Floor
        CA_Ceil
        CA_Abs
        CA_Sign
        CA_Re
        CA_Im
        CA_Arg
        CA_Conjugate
        CA_Pi
        CA_Sin
        CA_Cos
        CA_Exp
        CA_Log
        CA_Pow
        CA_Tan
        CA_Cot
        CA_Cosh
        CA_Sinh
        CA_Tanh
        CA_Coth
        CA_Atan
        CA_Acos
        CA_Asin
        CA_Acot
        CA_Atanh
        CA_Acosh
        CA_Asinh
        CA_Acoth
        CA_Euler
        CA_Gamma
        CA_LogGamma
        CA_Psi
        CA_Erf
        CA_Erfc
        CA_Erfi
        CA_RiemannZeta
        CA_HurwitzZeta
        CA_FUNC_CODE_LENGTH


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

    ctypedef ca_field_cache_struct ca_field_cache_t[1]

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

    ctypedef struct ca_factor_struct:
        pass
    ctypedef ca_factor_struct ca_factor_t[1]


    # flint/ca_poly.h
    ctypedef struct ca_poly_struct:
        pass
    ctypedef ca_poly_struct ca_poly_t[1]

    ctypedef struct ca_poly_vec_struct:
        pass
    ctypedef ca_poly_vec_struct ca_poly_vec_t[1]


    # flint/ca_vec.h
    ctypedef struct ca_vec_struct:
        pass
    ctypedef ca_vec_struct ca_vec_t[1]


    # flint/ca_mat.h
    ctypedef struct ca_mat_struct:
        pass
    ctypedef ca_mat_struct ca_mat_t[1]


    # flint/mpf_vec.h
    ctypedef __mpf_struct mpf


    # flint/mpf_mat.h
    ctypedef struct mpf_mat_struct:
        pass
    ctypedef mpf_mat_struct mpf_mat_t[1]


    # flint/mpfr_mat.h
    ctypedef struct mpfr_mat_struct:
        pass
    ctypedef mpfr_mat_struct mpfr_mat_t[1]


    # flint/gr.h
    ctypedef int (*gr_funcptr)()

    ctypedef struct gr_stream_struct:
        pass
    ctypedef gr_stream_struct gr_stream_t[1]

    ctypedef void * gr_ptr
    ctypedef const void * gr_srcptr
    ctypedef void * gr_ctx_ptr

    ctypedef struct gr_vec_struct:
        pass
    ctypedef gr_vec_struct gr_vec_t[1]

    ctypedef enum gr_method:
        GR_METHOD_CTX_WRITE
        GR_METHOD_CTX_CLEAR
        GR_METHOD_CTX_IS_RING
        GR_METHOD_CTX_IS_COMMUTATIVE_RING
        GR_METHOD_CTX_IS_INTEGRAL_DOMAIN
        GR_METHOD_CTX_IS_FIELD
        GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN
        GR_METHOD_CTX_IS_FINITE
        GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC
        GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED
        GR_METHOD_CTX_IS_ORDERED_RING
        GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP
        GR_METHOD_CTX_IS_EXACT
        GR_METHOD_CTX_IS_CANONICAL
        GR_METHOD_CTX_IS_THREADSAFE
        GR_METHOD_CTX_HAS_REAL_PREC
        GR_METHOD_CTX_SET_REAL_PREC
        GR_METHOD_CTX_GET_REAL_PREC
        GR_METHOD_CTX_SET_GEN_NAME
        GR_METHOD_CTX_SET_GEN_NAMES
        GR_METHOD_INIT
        GR_METHOD_CLEAR
        GR_METHOD_SWAP
        GR_METHOD_SET_SHALLOW
        GR_METHOD_WRITE
        GR_METHOD_WRITE_N
        GR_METHOD_RANDTEST
        GR_METHOD_RANDTEST_NOT_ZERO
        GR_METHOD_RANDTEST_SMALL
        GR_METHOD_ZERO
        GR_METHOD_ONE
        GR_METHOD_NEG_ONE
        GR_METHOD_IS_ZERO
        GR_METHOD_IS_ONE
        GR_METHOD_IS_NEG_ONE
        GR_METHOD_EQUAL
        GR_METHOD_SET
        GR_METHOD_SET_UI
        GR_METHOD_SET_SI
        GR_METHOD_SET_FMPZ
        GR_METHOD_SET_FMPQ
        GR_METHOD_SET_D
        GR_METHOD_SET_OTHER
        GR_METHOD_SET_STR
        GR_METHOD_GET_SI
        GR_METHOD_GET_UI
        GR_METHOD_GET_FMPZ
        GR_METHOD_GET_FMPQ
        GR_METHOD_GET_D
        GR_METHOD_GET_FEXPR
        GR_METHOD_GET_FEXPR_SERIALIZE
        GR_METHOD_SET_FEXPR
        GR_METHOD_NEG
        GR_METHOD_ADD
        GR_METHOD_ADD_UI
        GR_METHOD_ADD_SI
        GR_METHOD_ADD_FMPZ
        GR_METHOD_ADD_FMPQ
        GR_METHOD_ADD_OTHER
        GR_METHOD_OTHER_ADD
        GR_METHOD_SUB
        GR_METHOD_SUB_UI
        GR_METHOD_SUB_SI
        GR_METHOD_SUB_FMPZ
        GR_METHOD_SUB_FMPQ
        GR_METHOD_SUB_OTHER
        GR_METHOD_OTHER_SUB
        GR_METHOD_MUL
        GR_METHOD_MUL_UI
        GR_METHOD_MUL_SI
        GR_METHOD_MUL_FMPZ
        GR_METHOD_MUL_FMPQ
        GR_METHOD_MUL_OTHER
        GR_METHOD_OTHER_MUL
        GR_METHOD_ADDMUL
        GR_METHOD_ADDMUL_UI
        GR_METHOD_ADDMUL_SI
        GR_METHOD_ADDMUL_FMPZ
        GR_METHOD_ADDMUL_FMPQ
        GR_METHOD_ADDMUL_OTHER
        GR_METHOD_SUBMUL
        GR_METHOD_SUBMUL_UI
        GR_METHOD_SUBMUL_SI
        GR_METHOD_SUBMUL_FMPZ
        GR_METHOD_SUBMUL_FMPQ
        GR_METHOD_SUBMUL_OTHER
        GR_METHOD_FMA
        GR_METHOD_FMS
        GR_METHOD_FMMA
        GR_METHOD_FMMS
        GR_METHOD_MUL_TWO
        GR_METHOD_SQR
        GR_METHOD_MUL_2EXP_SI
        GR_METHOD_MUL_2EXP_FMPZ
        GR_METHOD_SET_FMPZ_2EXP_FMPZ
        GR_METHOD_GET_FMPZ_2EXP_FMPZ
        GR_METHOD_IS_INVERTIBLE
        GR_METHOD_INV
        GR_METHOD_DIV
        GR_METHOD_DIV_UI
        GR_METHOD_DIV_SI
        GR_METHOD_DIV_FMPZ
        GR_METHOD_DIV_FMPQ
        GR_METHOD_DIV_OTHER
        GR_METHOD_OTHER_DIV
        GR_METHOD_DIVEXACT
        GR_METHOD_DIVEXACT_UI
        GR_METHOD_DIVEXACT_SI
        GR_METHOD_DIVEXACT_FMPZ
        GR_METHOD_DIVEXACT_FMPQ
        GR_METHOD_DIVEXACT_OTHER
        GR_METHOD_OTHER_DIVEXACT
        GR_METHOD_POW
        GR_METHOD_POW_UI
        GR_METHOD_POW_SI
        GR_METHOD_POW_FMPZ
        GR_METHOD_POW_FMPQ
        GR_METHOD_POW_OTHER
        GR_METHOD_OTHER_POW
        GR_METHOD_IS_SQUARE
        GR_METHOD_SQRT
        GR_METHOD_RSQRT
        GR_METHOD_HYPOT
        GR_METHOD_IS_PERFECT_POWER
        GR_METHOD_FACTOR_PERFECT_POWER
        GR_METHOD_ROOT_UI
        GR_METHOD_DIVIDES
        GR_METHOD_EUCLIDEAN_DIV
        GR_METHOD_EUCLIDEAN_REM
        GR_METHOD_EUCLIDEAN_DIVREM
        GR_METHOD_GCD
        GR_METHOD_LCM
        GR_METHOD_FACTOR
        GR_METHOD_NUMERATOR
        GR_METHOD_DENOMINATOR
        GR_METHOD_FLOOR
        GR_METHOD_CEIL
        GR_METHOD_TRUNC
        GR_METHOD_NINT
        GR_METHOD_CMP
        GR_METHOD_CMPABS
        GR_METHOD_CMP_OTHER
        GR_METHOD_CMPABS_OTHER
        GR_METHOD_MIN
        GR_METHOD_MAX
        GR_METHOD_ABS
        GR_METHOD_ABS2
        GR_METHOD_I
        GR_METHOD_CONJ
        GR_METHOD_RE
        GR_METHOD_IM
        GR_METHOD_SGN
        GR_METHOD_CSGN
        GR_METHOD_ARG
        GR_METHOD_POS_INF
        GR_METHOD_NEG_INF
        GR_METHOD_UINF
        GR_METHOD_UNDEFINED
        GR_METHOD_UNKNOWN
        GR_METHOD_IS_INTEGER
        GR_METHOD_IS_RATIONAL
        GR_METHOD_IS_REAL
        GR_METHOD_IS_POSITIVE_INTEGER
        GR_METHOD_IS_NONNEGATIVE_INTEGER
        GR_METHOD_IS_NEGATIVE_INTEGER
        GR_METHOD_IS_NONPOSITIVE_INTEGER
        GR_METHOD_IS_POSITIVE_REAL
        GR_METHOD_IS_NONNEGATIVE_REAL
        GR_METHOD_IS_NEGATIVE_REAL
        GR_METHOD_IS_NONPOSITIVE_REAL
        GR_METHOD_IS_ROOT_OF_UNITY
        GR_METHOD_ROOT_OF_UNITY_UI
        GR_METHOD_ROOT_OF_UNITY_UI_VEC
        GR_METHOD_PI
        GR_METHOD_EXP
        GR_METHOD_EXPM1
        GR_METHOD_EXP_PI_I
        GR_METHOD_EXP2
        GR_METHOD_EXP10
        GR_METHOD_LOG
        GR_METHOD_LOG1P
        GR_METHOD_LOG_PI_I
        GR_METHOD_LOG2
        GR_METHOD_LOG10
        GR_METHOD_SIN
        GR_METHOD_COS
        GR_METHOD_SIN_COS
        GR_METHOD_TAN
        GR_METHOD_COT
        GR_METHOD_SEC
        GR_METHOD_CSC
        GR_METHOD_SIN_PI
        GR_METHOD_COS_PI
        GR_METHOD_SIN_COS_PI
        GR_METHOD_TAN_PI
        GR_METHOD_COT_PI
        GR_METHOD_SEC_PI
        GR_METHOD_CSC_PI
        GR_METHOD_SINC
        GR_METHOD_SINC_PI
        GR_METHOD_SINH
        GR_METHOD_COSH
        GR_METHOD_SINH_COSH
        GR_METHOD_TANH
        GR_METHOD_COTH
        GR_METHOD_SECH
        GR_METHOD_CSCH
        GR_METHOD_ASIN
        GR_METHOD_ACOS
        GR_METHOD_ATAN
        GR_METHOD_ATAN2
        GR_METHOD_ACOT
        GR_METHOD_ASEC
        GR_METHOD_ACSC
        GR_METHOD_ASINH
        GR_METHOD_ACOSH
        GR_METHOD_ATANH
        GR_METHOD_ACOTH
        GR_METHOD_ASECH
        GR_METHOD_ACSCH
        GR_METHOD_ASIN_PI
        GR_METHOD_ACOS_PI
        GR_METHOD_ATAN_PI
        GR_METHOD_ACOT_PI
        GR_METHOD_ASEC_PI
        GR_METHOD_ACSC_PI
        GR_METHOD_FAC
        GR_METHOD_FAC_UI
        GR_METHOD_FAC_FMPZ
        GR_METHOD_FAC_VEC
        GR_METHOD_RFAC
        GR_METHOD_RFAC_UI
        GR_METHOD_RFAC_FMPZ
        GR_METHOD_RFAC_VEC
        GR_METHOD_BIN
        GR_METHOD_BIN_UI
        GR_METHOD_BIN_UIUI
        GR_METHOD_BIN_VEC
        GR_METHOD_BIN_UI_VEC
        GR_METHOD_RISING_UI
        GR_METHOD_RISING
        GR_METHOD_FALLING_UI
        GR_METHOD_FALLING
        GR_METHOD_GAMMA
        GR_METHOD_GAMMA_FMPZ
        GR_METHOD_GAMMA_FMPQ
        GR_METHOD_RGAMMA
        GR_METHOD_LGAMMA
        GR_METHOD_DIGAMMA
        GR_METHOD_BETA
        GR_METHOD_DOUBLEFAC
        GR_METHOD_DOUBLEFAC_UI
        GR_METHOD_BARNES_G
        GR_METHOD_LOG_BARNES_G
        GR_METHOD_HARMONIC
        GR_METHOD_HARMONIC_UI
        GR_METHOD_BERNOULLI_UI
        GR_METHOD_BERNOULLI_FMPZ
        GR_METHOD_BERNOULLI_VEC
        GR_METHOD_FIB_UI
        GR_METHOD_FIB_FMPZ
        GR_METHOD_FIB_VEC
        GR_METHOD_STIRLING_S1U_UIUI
        GR_METHOD_STIRLING_S1_UIUI
        GR_METHOD_STIRLING_S2_UIUI
        GR_METHOD_STIRLING_S1U_UI_VEC
        GR_METHOD_STIRLING_S1_UI_VEC
        GR_METHOD_STIRLING_S2_UI_VEC
        GR_METHOD_EULERNUM_UI
        GR_METHOD_EULERNUM_FMPZ
        GR_METHOD_EULERNUM_VEC
        GR_METHOD_BELLNUM_UI
        GR_METHOD_BELLNUM_FMPZ
        GR_METHOD_BELLNUM_VEC
        GR_METHOD_PARTITIONS_UI
        GR_METHOD_PARTITIONS_FMPZ
        GR_METHOD_PARTITIONS_VEC
        GR_METHOD_CHEBYSHEV_T_FMPZ
        GR_METHOD_CHEBYSHEV_U_FMPZ
        GR_METHOD_CHEBYSHEV_T
        GR_METHOD_CHEBYSHEV_U
        GR_METHOD_JACOBI_P
        GR_METHOD_GEGENBAUER_C
        GR_METHOD_LAGUERRE_L
        GR_METHOD_HERMITE_H
        GR_METHOD_LEGENDRE_P
        GR_METHOD_LEGENDRE_Q
        GR_METHOD_LEGENDRE_P_ROOT_UI
        GR_METHOD_SPHERICAL_Y_SI
        GR_METHOD_EULER
        GR_METHOD_CATALAN
        GR_METHOD_KHINCHIN
        GR_METHOD_GLAISHER
        GR_METHOD_LAMBERTW
        GR_METHOD_LAMBERTW_FMPZ
        GR_METHOD_ERF
        GR_METHOD_ERFC
        GR_METHOD_ERFCX
        GR_METHOD_ERFI
        GR_METHOD_ERFINV
        GR_METHOD_ERFCINV
        GR_METHOD_FRESNEL_S
        GR_METHOD_FRESNEL_C
        GR_METHOD_FRESNEL
        GR_METHOD_GAMMA_UPPER
        GR_METHOD_GAMMA_LOWER
        GR_METHOD_BETA_LOWER
        GR_METHOD_EXP_INTEGRAL
        GR_METHOD_EXP_INTEGRAL_EI
        GR_METHOD_SIN_INTEGRAL
        GR_METHOD_COS_INTEGRAL
        GR_METHOD_SINH_INTEGRAL
        GR_METHOD_COSH_INTEGRAL
        GR_METHOD_LOG_INTEGRAL
        GR_METHOD_DILOG
        GR_METHOD_BESSEL_J
        GR_METHOD_BESSEL_Y
        GR_METHOD_BESSEL_J_Y
        GR_METHOD_BESSEL_I
        GR_METHOD_BESSEL_I_SCALED
        GR_METHOD_BESSEL_K
        GR_METHOD_BESSEL_K_SCALED
        GR_METHOD_AIRY
        GR_METHOD_AIRY_AI
        GR_METHOD_AIRY_BI
        GR_METHOD_AIRY_AI_PRIME
        GR_METHOD_AIRY_BI_PRIME
        GR_METHOD_AIRY_AI_ZERO
        GR_METHOD_AIRY_BI_ZERO
        GR_METHOD_AIRY_AI_PRIME_ZERO
        GR_METHOD_AIRY_BI_PRIME_ZERO
        GR_METHOD_COULOMB
        GR_METHOD_COULOMB_F
        GR_METHOD_COULOMB_G
        GR_METHOD_COULOMB_HPOS
        GR_METHOD_COULOMB_HNEG
        GR_METHOD_ZETA
        GR_METHOD_ZETA_UI
        GR_METHOD_HURWITZ_ZETA
        GR_METHOD_LERCH_PHI
        GR_METHOD_STIELTJES
        GR_METHOD_DIRICHLET_ETA
        GR_METHOD_DIRICHLET_BETA
        GR_METHOD_RIEMANN_XI
        GR_METHOD_ZETA_ZERO
        GR_METHOD_ZETA_ZERO_VEC
        GR_METHOD_ZETA_NZEROS
        GR_METHOD_DIRICHLET_CHI_UI
        GR_METHOD_DIRICHLET_CHI_FMPZ
        GR_METHOD_DIRICHLET_L
        GR_METHOD_DIRICHLET_HARDY_THETA
        GR_METHOD_DIRICHLET_HARDY_Z
        GR_METHOD_BERNPOLY_UI
        GR_METHOD_EULERPOLY_UI
        GR_METHOD_POLYLOG
        GR_METHOD_POLYGAMMA
        GR_METHOD_HYPGEOM_0F1
        GR_METHOD_HYPGEOM_1F1
        GR_METHOD_HYPGEOM_2F0
        GR_METHOD_HYPGEOM_2F1
        GR_METHOD_HYPGEOM_U
        GR_METHOD_HYPGEOM_PFQ
        GR_METHOD_JACOBI_THETA
        GR_METHOD_JACOBI_THETA_1
        GR_METHOD_JACOBI_THETA_2
        GR_METHOD_JACOBI_THETA_3
        GR_METHOD_JACOBI_THETA_4
        GR_METHOD_JACOBI_THETA_Q
        GR_METHOD_JACOBI_THETA_Q_1
        GR_METHOD_JACOBI_THETA_Q_2
        GR_METHOD_JACOBI_THETA_Q_3
        GR_METHOD_JACOBI_THETA_Q_4
        GR_METHOD_MODULAR_J
        GR_METHOD_MODULAR_LAMBDA
        GR_METHOD_MODULAR_DELTA
        GR_METHOD_HILBERT_CLASS_POLY
        GR_METHOD_DEDEKIND_ETA
        GR_METHOD_DEDEKIND_ETA_Q
        GR_METHOD_EISENSTEIN_E
        GR_METHOD_EISENSTEIN_G
        GR_METHOD_EISENSTEIN_G_VEC
        GR_METHOD_AGM
        GR_METHOD_AGM1
        GR_METHOD_ELLIPTIC_K
        GR_METHOD_ELLIPTIC_E
        GR_METHOD_ELLIPTIC_PI
        GR_METHOD_ELLIPTIC_F
        GR_METHOD_ELLIPTIC_E_INC
        GR_METHOD_ELLIPTIC_PI_INC
        GR_METHOD_CARLSON_RF
        GR_METHOD_CARLSON_RC
        GR_METHOD_CARLSON_RJ
        GR_METHOD_CARLSON_RG
        GR_METHOD_CARLSON_RD
        GR_METHOD_ELLIPTIC_ROOTS
        GR_METHOD_ELLIPTIC_INVARIANTS
        GR_METHOD_WEIERSTRASS_P
        GR_METHOD_WEIERSTRASS_P_PRIME
        GR_METHOD_WEIERSTRASS_P_INV
        GR_METHOD_WEIERSTRASS_ZETA
        GR_METHOD_WEIERSTRASS_SIGMA
        GR_METHOD_GEN
        GR_METHOD_GENS
        GR_METHOD_CTX_FQ_PRIME
        GR_METHOD_CTX_FQ_DEGREE
        GR_METHOD_CTX_FQ_ORDER
        GR_METHOD_FQ_FROBENIUS
        GR_METHOD_FQ_MULTIPLICATIVE_ORDER
        GR_METHOD_FQ_NORM
        GR_METHOD_FQ_TRACE
        GR_METHOD_FQ_IS_PRIMITIVE
        GR_METHOD_FQ_PTH_ROOT
        GR_METHOD_VEC_INIT
        GR_METHOD_VEC_CLEAR
        GR_METHOD_VEC_SWAP
        GR_METHOD_VEC_SET
        GR_METHOD_VEC_ZERO
        GR_METHOD_VEC_EQUAL
        GR_METHOD_VEC_IS_ZERO
        GR_METHOD_VEC_NEG
        GR_METHOD_VEC_NORMALISE
        GR_METHOD_VEC_NORMALISE_WEAK
        GR_METHOD_VEC_ADD
        GR_METHOD_VEC_SUB
        GR_METHOD_VEC_MUL
        GR_METHOD_VEC_DIV
        GR_METHOD_VEC_DIVEXACT
        GR_METHOD_VEC_POW
        GR_METHOD_VEC_ADD_SCALAR
        GR_METHOD_VEC_SUB_SCALAR
        GR_METHOD_VEC_MUL_SCALAR
        GR_METHOD_VEC_DIV_SCALAR
        GR_METHOD_VEC_DIVEXACT_SCALAR
        GR_METHOD_VEC_POW_SCALAR
        GR_METHOD_SCALAR_ADD_VEC
        GR_METHOD_SCALAR_SUB_VEC
        GR_METHOD_SCALAR_MUL_VEC
        GR_METHOD_SCALAR_DIV_VEC
        GR_METHOD_SCALAR_DIVEXACT_VEC
        GR_METHOD_SCALAR_POW_VEC
        GR_METHOD_VEC_ADD_OTHER
        GR_METHOD_VEC_SUB_OTHER
        GR_METHOD_VEC_MUL_OTHER
        GR_METHOD_VEC_DIV_OTHER
        GR_METHOD_VEC_DIVEXACT_OTHER
        GR_METHOD_VEC_POW_OTHER
        GR_METHOD_OTHER_ADD_VEC
        GR_METHOD_OTHER_SUB_VEC
        GR_METHOD_OTHER_MUL_VEC
        GR_METHOD_OTHER_DIV_VEC
        GR_METHOD_OTHER_DIVEXACT_VEC
        GR_METHOD_OTHER_POW_VEC
        GR_METHOD_VEC_ADD_SCALAR_OTHER
        GR_METHOD_VEC_SUB_SCALAR_OTHER
        GR_METHOD_VEC_MUL_SCALAR_OTHER
        GR_METHOD_VEC_DIV_SCALAR_OTHER
        GR_METHOD_VEC_DIVEXACT_SCALAR_OTHER
        GR_METHOD_VEC_POW_SCALAR_OTHER
        GR_METHOD_SCALAR_OTHER_ADD_VEC
        GR_METHOD_SCALAR_OTHER_SUB_VEC
        GR_METHOD_SCALAR_OTHER_MUL_VEC
        GR_METHOD_SCALAR_OTHER_DIV_VEC
        GR_METHOD_SCALAR_OTHER_DIVEXACT_VEC
        GR_METHOD_SCALAR_OTHER_POW_VEC
        GR_METHOD_VEC_ADD_SCALAR_SI
        GR_METHOD_VEC_ADD_SCALAR_UI
        GR_METHOD_VEC_ADD_SCALAR_FMPZ
        GR_METHOD_VEC_ADD_SCALAR_FMPQ
        GR_METHOD_VEC_SUB_SCALAR_SI
        GR_METHOD_VEC_SUB_SCALAR_UI
        GR_METHOD_VEC_SUB_SCALAR_FMPZ
        GR_METHOD_VEC_SUB_SCALAR_FMPQ
        GR_METHOD_VEC_MUL_SCALAR_SI
        GR_METHOD_VEC_MUL_SCALAR_UI
        GR_METHOD_VEC_MUL_SCALAR_FMPZ
        GR_METHOD_VEC_MUL_SCALAR_FMPQ
        GR_METHOD_VEC_DIV_SCALAR_SI
        GR_METHOD_VEC_DIV_SCALAR_UI
        GR_METHOD_VEC_DIV_SCALAR_FMPZ
        GR_METHOD_VEC_DIV_SCALAR_FMPQ
        GR_METHOD_VEC_DIVEXACT_SCALAR_SI
        GR_METHOD_VEC_DIVEXACT_SCALAR_UI
        GR_METHOD_VEC_DIVEXACT_SCALAR_FMPZ
        GR_METHOD_VEC_DIVEXACT_SCALAR_FMPQ
        GR_METHOD_VEC_POW_SCALAR_SI
        GR_METHOD_VEC_POW_SCALAR_UI
        GR_METHOD_VEC_POW_SCALAR_FMPZ
        GR_METHOD_VEC_POW_SCALAR_FMPQ
        GR_METHOD_VEC_MUL_SCALAR_2EXP_SI
        GR_METHOD_VEC_MUL_SCALAR_2EXP_FMPZ
        GR_METHOD_VEC_ADDMUL_SCALAR
        GR_METHOD_VEC_SUBMUL_SCALAR
        GR_METHOD_VEC_ADDMUL_SCALAR_SI
        GR_METHOD_VEC_SUBMUL_SCALAR_SI
        GR_METHOD_VEC_SUM
        GR_METHOD_VEC_PRODUCT
        GR_METHOD_VEC_DOT
        GR_METHOD_VEC_DOT_REV
        GR_METHOD_VEC_DOT_UI
        GR_METHOD_VEC_DOT_SI
        GR_METHOD_VEC_DOT_FMPZ
        GR_METHOD_VEC_SET_POWERS
        GR_METHOD_VEC_RECIPROCALS
        GR_METHOD_POLY_MULLOW
        GR_METHOD_POLY_DIV
        GR_METHOD_POLY_DIVREM
        GR_METHOD_POLY_DIVEXACT
        GR_METHOD_POLY_TAYLOR_SHIFT
        GR_METHOD_POLY_INV_SERIES
        GR_METHOD_POLY_INV_SERIES_BASECASE
        GR_METHOD_POLY_DIV_SERIES
        GR_METHOD_POLY_DIV_SERIES_BASECASE
        GR_METHOD_POLY_RSQRT_SERIES
        GR_METHOD_POLY_SQRT_SERIES
        GR_METHOD_POLY_EXP_SERIES
        GR_METHOD_POLY_FACTOR
        GR_METHOD_POLY_ROOTS
        GR_METHOD_POLY_ROOTS_OTHER
        GR_METHOD_MAT_MUL
        GR_METHOD_MAT_DET
        GR_METHOD_MAT_EXP
        GR_METHOD_MAT_LOG
        GR_METHOD_MAT_FIND_NONZERO_PIVOT
        GR_METHOD_MAT_DIAGONALIZATION
        GR_METHOD_TAB_SIZE

    ctypedef struct gr_method_tab_input:
        pass

    ctypedef enum gr_which_structure:
        GR_CTX_FMPZ
        GR_CTX_FMPQ
        GR_CTX_FMPZI
        GR_CTX_FMPZ_MOD
        GR_CTX_NMOD
        GR_CTX_NMOD8
        GR_CTX_NMOD32
        GR_CTX_FQ
        GR_CTX_FQ_NMOD
        GR_CTX_FQ_ZECH
        GR_CTX_NF
        GR_CTX_REAL_ALGEBRAIC_QQBAR
        GR_CTX_COMPLEX_ALGEBRAIC_QQBAR
        GR_CTX_REAL_ALGEBRAIC_CA
        GR_CTX_COMPLEX_ALGEBRAIC_CA
        GR_CTX_RR_CA
        GR_CTX_CC_CA
        GR_CTX_COMPLEX_EXTENDED_CA
        GR_CTX_RR_ARB
        GR_CTX_CC_ACB
        GR_CTX_REAL_FLOAT_ARF
        GR_CTX_COMPLEX_FLOAT_ACF
        GR_CTX_FMPZ_POLY
        GR_CTX_FMPQ_POLY
        GR_CTX_GR_POLY
        GR_CTX_FMPZ_MPOLY
        GR_CTX_GR_MPOLY
        GR_CTX_FMPZ_MPOLY_Q
        GR_CTX_GR_SERIES
        GR_CTX_GR_SERIES_MOD
        GR_CTX_GR_MAT
        GR_CTX_GR_VEC
        GR_CTX_PSL2Z
        GR_CTX_DIRICHLET_GROUP
        GR_CTX_PERM
        GR_CTX_FEXPR
        GR_CTX_UNKNOWN_DOMAIN
        GR_CTX_WHICH_STRUCTURE_TAB_SIZE

    ctypedef struct gr_ctx_struct:
        pass
    ctypedef gr_ctx_struct gr_ctx_t[1]

    ctypedef void ((*gr_method_init_clear_op)(gr_ptr, gr_ctx_ptr))
    ctypedef void ((*gr_method_swap_op)(gr_ptr, gr_ptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_ctx)(gr_ctx_ptr))
    # NOTE: we removed an extra parenthesis so that Cython is less confused
    # see https://github.com/cython/cython/issues/5779
    ctypedef truth_t (*gr_method_ctx_predicate)(gr_ctx_ptr)
    ctypedef int ((*gr_method_ctx_set_si)(gr_ctx_ptr, slong))
    ctypedef int ((*gr_method_ctx_get_si)(slong *, gr_ctx_ptr))
    ctypedef int ((*gr_method_ctx_stream)(gr_stream_t, gr_ctx_ptr))
    ctypedef int ((*gr_method_ctx_set_str)(gr_ctx_ptr, const char *))
    ctypedef int ((*gr_method_ctx_set_strs)(gr_ctx_ptr, const char **))
    ctypedef int ((*gr_method_stream_in)(gr_stream_t, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_stream_in_si)(gr_stream_t, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_randtest)(gr_ptr, flint_rand_t state, gr_ctx_ptr))
    ctypedef int ((*gr_method_constant_op)(gr_ptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_constant_op_get_si)(slong *, gr_ctx_ptr))
    ctypedef int ((*gr_method_constant_op_get_fmpz)(fmpz_t, gr_ctx_ptr))
    ctypedef void ((*gr_method_void_unary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_si)(gr_ptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_ui)(gr_ptr, ulong, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_fmpz)(gr_ptr, const fmpz_t, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_fmpq)(gr_ptr, const fmpq_t, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_d)(gr_ptr, double, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_other)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_str)(gr_ptr, const char *, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_get_ui)(ulong *, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_get_si)(slong *, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_get_fmpz)(fmpz_t, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_get_fmpq)(fmpq_t, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_get_d)(double *, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_get_fmpz_fmpz)(fmpz_t, fmpz_t, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_unary_op_with_flag)(gr_ptr, gr_srcptr, int, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_unary_op)(gr_ptr, gr_ptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_unary_op_with_flag)(gr_ptr, gr_ptr, gr_srcptr, int, gr_ctx_ptr))
    ctypedef int ((*gr_method_quaternary_unary_op)(gr_ptr, gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_si)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_ui)(gr_ptr, gr_srcptr, ulong, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_fmpz)(gr_ptr, gr_srcptr, const fmpz_t, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_fmpz_fmpz)(gr_ptr, const fmpz_t, const fmpz_t, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_fmpz_si)(gr_ptr, const fmpz_t, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_fmpq)(gr_ptr, gr_srcptr, const fmpq_t, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_other)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_other_binary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_si_binary_op)(gr_ptr, slong, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_ui_binary_op)(gr_ptr, ulong, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_fmpz_binary_op)(gr_ptr, const fmpz_t, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_fmpq_binary_op)(gr_ptr, const fmpq_t, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_ui_ui)(gr_ptr, ulong, ulong, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_ui_si)(gr_ptr, ulong, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_get_int)(int *, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_other_get_int)(int *, gr_srcptr, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_binary_op)(gr_ptr, gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_op_with_flag)(gr_ptr, gr_srcptr, gr_srcptr, int, gr_ctx_ptr))
    ctypedef int ((*gr_method_binary_binary_op_ui_ui)(gr_ptr, gr_ptr, ulong, ulong, gr_ctx_ptr))
    ctypedef int ((*gr_method_ternary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_ternary_op_with_flag)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, int, gr_ctx_ptr))
    ctypedef int ((*gr_method_ternary_unary_op)(gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_quaternary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_quaternary_op_with_flag)(gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_srcptr, int, gr_ctx_ptr))
    ctypedef int ((*gr_method_quaternary_binary_op)(gr_ptr, gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_quaternary_ternary_op)(gr_ptr, gr_ptr, gr_ptr, gr_ptr, gr_srcptr, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_si_si_quaternary_op)(gr_ptr, slong, slong, gr_srcptr, gr_srcptr, gr_ctx_ptr))
    # NOTE: we removed an extra parenthesis so that Cython is less confused
    # see https://github.com/cython/cython/issues/5779
    ctypedef truth_t (*gr_method_unary_predicate)(gr_srcptr, gr_ctx_ptr)
    ctypedef truth_t (*gr_method_binary_predicate)(gr_srcptr, gr_srcptr, gr_ctx_ptr)
    ctypedef void ((*gr_method_vec_init_clear_op)(gr_ptr, slong, gr_ctx_ptr))
    ctypedef void ((*gr_method_vec_swap_op)(gr_ptr, gr_ptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_constant_op)(gr_ptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_op)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_vec_op)(gr_ptr, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_scalar_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_scalar_vec_op)(gr_ptr, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_op_other)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_other_op_vec)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_op_scalar_other)(gr_ptr, gr_srcptr, slong, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr))
    ctypedef int ((*gr_method_scalar_other_op_vec)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_scalar_op_si)(gr_ptr, gr_srcptr, slong, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_scalar_op_ui)(gr_ptr, gr_srcptr, slong, ulong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_scalar_op_fmpz)(gr_ptr, gr_srcptr, slong, const fmpz_t, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_scalar_op_fmpq)(gr_ptr, gr_srcptr, slong, const fmpq_t, gr_ctx_ptr))
    # NOTE: we removed an extra parenthesis so that Cython is less confused
    # see https://github.com/cython/cython/issues/5779
    ctypedef truth_t (*gr_method_vec_predicate)(gr_srcptr, slong, gr_ctx_ptr)
    ctypedef truth_t (*gr_method_vec_vec_predicate)(gr_srcptr, gr_srcptr, slong, gr_ctx_ptr)
    ctypedef int ((*gr_method_factor_op)(gr_ptr, gr_vec_t, gr_vec_t, gr_srcptr, int, gr_ctx_ptr))
    ctypedef int ((*gr_method_poly_unary_trunc_op)(gr_ptr, gr_srcptr, slong, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_poly_binary_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_poly_binary_binary_op)(gr_ptr, gr_ptr, gr_srcptr, slong, gr_srcptr, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_poly_binary_trunc_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, slong, slong, gr_ctx_ptr))
    ctypedef int ((*gr_method_vec_ctx_op)(gr_vec_t, gr_ctx_ptr))

    ctypedef struct polynomial_ctx_t:
        pass

    ctypedef struct vector_ctx_t:
        pass

    ctypedef struct matrix_ctx_t:
        pass


    # flint/gr_poly.h
    ctypedef struct gr_poly_struct:
        pass
    ctypedef gr_poly_struct gr_poly_t[1]


    # flint/gr_mat.h
    ctypedef struct gr_mat_struct:
        pass
    ctypedef gr_mat_struct gr_mat_t[1]


    # flint/gr_mpoly.h
    ctypedef struct gr_mpoly_struct:
        pass
    ctypedef gr_mpoly_struct gr_mpoly_t[1]


    # flint/double_interval.h
    ctypedef struct di_t:
        double a
        double b


    # flint/fq_default.h
    ctypedef union fq_default_struct:
        fq_t fq
        fq_nmod_t fq_nmod
        fq_zech_t fq_zech
        ulong nmod
        fmpz_t fmpz_mod
    ctypedef fq_default_struct fq_default_t[1]

    ctypedef struct fq_default_ctx_struct:
        pass
    ctypedef fq_default_ctx_struct fq_default_ctx_t[1]


    # flint/fq_default_poly.h
    ctypedef union fq_default_poly_struct:
        fq_poly_t fq
        fq_nmod_poly_t fq_nmod
        fq_zech_poly_t fq_zech
        nmod_poly_t nmod
        fmpz_mod_poly_t fmpz_mod
    ctypedef fq_default_poly_struct fq_default_poly_t[1]


    # flint/fq_poly_factor.h
    ctypedef union fq_default_poly_factor_struct:
        fq_poly_factor_t fq
        fq_nmod_poly_factor_t fq_nmod
        fq_zech_poly_factor_t fq_zech
        nmod_poly_factor_t nmod
        fmpz_mod_poly_factor_t fmpz_mod
    ctypedef fq_default_poly_factor_struct fq_default_poly_factor_t[1]


    # flint/fq_default_mat.h
    ctypedef union fq_default_mat_struct:
        pass
    ctypedef fq_default_mat_struct fq_default_mat_t[1]


    # flint/qfb.h
    ctypedef struct qfb:
        pass
    ctypedef qfb qfb_t[1]

    ctypedef struct qfb_hash_t:
        pass


    # flint/qqbar.h
    ctypedef struct qqbar_struct:
        pass
    ctypedef qqbar_struct qqbar_t[1]
    ctypedef qqbar_struct * qqbar_ptr
    ctypedef const qqbar_struct * qqbar_srcptr


    # flint/profiler.h
    ctypedef struct struct_meminfo:
        pass
    ctypedef struct_meminfo meminfo_t[1]

    ctypedef struct struct_timeit:
        pass
    ctypedef struct_timeit timeit_t[1]

    ctypedef void (*profile_target_t)(void* arg, ulong count)


    # flint/fexpr.h
    ctypedef struct fexpr_struct:
        pass
    ctypedef fexpr_struct fexpr_t[1]
    ctypedef fexpr_struct * fexpr_ptr
    ctypedef const fexpr_struct * fexpr_srcptr

    ctypedef struct fexpr_vec_struct:
        pass
    ctypedef fexpr_vec_struct fexpr_vec_t[1]


    # flint/hypgeom.h
    ctypedef struct hypgeom_struct:
        pass
    ctypedef hypgeom_struct hypgeom_t[1]


    # flint/dlog.h
    ctypedef struct dlog_precomp_struct:
        pass
    ctypedef dlog_precomp_struct * dlog_precomp_ptr
    ctypedef dlog_precomp_struct dlog_precomp_t[1]

    ctypedef struct dlog_1modpe_struct:
        pass
    ctypedef dlog_1modpe_struct dlog_1modpe_t[1]

    ctypedef struct dlog_modpe_struct:
        pass
    ctypedef dlog_modpe_struct dlog_modpe_t[1]

    ctypedef struct dlog_table_struct:
        pass
    ctypedef dlog_table_struct dlog_table_t[1]

    ctypedef struct apow_t:
        pass

    ctypedef struct dlog_bsgs_struct:
        pass
    ctypedef dlog_bsgs_struct dlog_bsgs_t[1]

    ctypedef struct dlog_rho_struct:
        pass
    ctypedef dlog_rho_struct dlog_rho_t[1]

    ctypedef struct dlog_crt_struct:
        pass
    ctypedef dlog_crt_struct dlog_crt_t[1]

    ctypedef struct dlog_power_struct:
        pass
    ctypedef dlog_power_struct dlog_power_t[1]

    ctypedef ulong dlog_order23_t[1]


    # flint/bool_mat.h
    ctypedef struct bool_mat_struct:
        pass
    ctypedef bool_mat_struct bool_mat_t[1]


    # flint/dirichlet.h
    ctypedef struct dirichlet_prime_group_struct:
        pass

    ctypedef struct dirichlet_group_struct:
        pass
    ctypedef dirichlet_group_struct dirichlet_group_t[1]

    ctypedef struct dirichlet_char_struct:
        pass
    ctypedef dirichlet_char_struct dirichlet_char_t[1]


    # flint/arb_fpwrap.h
    ctypedef struct complex_double:
        double real
        double imag


    # flint/aprcl.h
    ctypedef struct _aprcl_config:
        pass
    ctypedef _aprcl_config aprcl_config[1];

    ctypedef struct _unity_zpq:
        pass
    ctypedef _unity_zpq unity_zpq[1]

    ctypedef struct _unity_zp:
        pass
    ctypedef _unity_zp unity_zp[1]

    ctypedef enum primality_test_status:
        UNKNOWN
        PRIME
        COMPOSITE
        PROBABPRIME


    # flint/acf.h
    ctypedef struct acf_struct:
        pass

    ctypedef acf_struct acf_t[1]
    ctypedef acf_struct * acf_ptr
    ctypedef const acf_struct * acf_srcptr


    # flint/fmpz_lll.h
    ctypedef enum rep_type:
        GRAM
        Z_BASIS

    ctypedef enum gram_type:
        APPROX
        EXACT

    ctypedef struct fmpz_lll_struct:
        pass
    ctypedef fmpz_lll_struct fmpz_lll_t[1]

    ctypedef union fmpz_gram_union:
        d_mat_t appSP
        mpf_mat_t appSP2
        fmpz_mat_t exactSP
    ctypedef fmpz_gram_union fmpz_gram_t[1]
