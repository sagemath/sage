# distutils: libraries = flint
# distutils: depends = flint/dlog.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    unsigned long dlog_once(unsigned long b, unsigned long a, const nmod_t mod, unsigned long n)

    void dlog_precomp_n_init(dlog_precomp_t pre, unsigned long a, unsigned long mod, unsigned long n, unsigned long num)

    unsigned long dlog_precomp(const dlog_precomp_t pre, unsigned long b)

    void dlog_precomp_clear(dlog_precomp_t pre)

    void dlog_precomp_modpe_init(dlog_precomp_t pre, unsigned long a, unsigned long p, unsigned long e, unsigned long pe, unsigned long num)

    void dlog_precomp_p_init(dlog_precomp_t pre, unsigned long a, unsigned long mod, unsigned long p, unsigned long num)

    void dlog_precomp_pe_init(dlog_precomp_t pre, unsigned long a, unsigned long mod, unsigned long p, unsigned long e, unsigned long pe, unsigned long num)

    void dlog_precomp_small_init(dlog_precomp_t pre, unsigned long a, unsigned long mod, unsigned long n, unsigned long num)

    void dlog_vec_fill(unsigned long * v, unsigned long nv, unsigned long x)

    void dlog_vec_set_not_found(unsigned long * v, unsigned long nv, nmod_t mod)

    void dlog_vec(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    void dlog_vec_add(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    void dlog_vec_loop(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    void dlog_vec_loop_add(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    void dlog_vec_eratos(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    void dlog_vec_eratos_add(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    void dlog_vec_sieve_add(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    void dlog_vec_sieve(unsigned long * v, unsigned long nv, unsigned long a, unsigned long va, nmod_t mod, unsigned long na, nmod_t order)

    unsigned long dlog_table_init(dlog_table_t t, unsigned long a, unsigned long mod)

    void dlog_table_clear(dlog_table_t t)

    unsigned long dlog_table(const dlog_table_t t, unsigned long b)

    unsigned long dlog_bsgs_init(dlog_bsgs_t t, unsigned long a, unsigned long mod, unsigned long n, unsigned long m)

    void dlog_bsgs_clear(dlog_bsgs_t t)

    unsigned long dlog_bsgs(const dlog_bsgs_t t, unsigned long b)

    unsigned long dlog_modpe_init(dlog_modpe_t t, unsigned long a, unsigned long p, unsigned long e, unsigned long pe, unsigned long num)

    void dlog_modpe_clear(dlog_modpe_t t)

    unsigned long dlog_modpe(const dlog_modpe_t t, unsigned long b)

    unsigned long dlog_crt_init(dlog_crt_t t, unsigned long a, unsigned long mod, unsigned long n, unsigned long num)

    void dlog_crt_clear(dlog_crt_t t)

    unsigned long dlog_crt(const dlog_crt_t t, unsigned long b)

    unsigned long dlog_power_init(dlog_power_t t, unsigned long a, unsigned long mod, unsigned long p, unsigned long e, unsigned long num)

    void dlog_power_clear(dlog_power_t t)

    unsigned long dlog_power(const dlog_power_t t, unsigned long b)

    void dlog_rho_init(dlog_rho_t t, unsigned long a, unsigned long mod, unsigned long n)

    void dlog_rho_clear(dlog_rho_t t)

    unsigned long dlog_rho(const dlog_rho_t t, unsigned long b)
