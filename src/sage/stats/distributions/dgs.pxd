# sage_setup: distribution = sagemath-modules
"""
AUTHOR: Martin Albrecht <martinralbrecht@googlemail.com>
"""

from sage.libs.gmp.mpz cimport mpz_t
from sage.libs.gmp.random cimport gmp_randstate_t
from sage.libs.mpfr.types cimport mpfr_t
from libc.stdint cimport uint64_t

cdef extern from "dgs.h":
    int DGS_BERN_UNIFORM_DEFAULT_LENGTH
    int DGS_BERN_UNIFORM_MAX_LENGTH

    ctypedef struct dgs_bern_uniform_t:
        size_t length
        size_t count
        mpz_t tmp
        uint64_t pool

    dgs_bern_uniform_t *dgs_bern_uniform_init(int length)
    unsigned long dgs_bern_uniform_call(dgs_bern_uniform_t *self, gmp_randstate_t state)
    void dgs_bern_uniform_clear(dgs_bern_uniform_t *self)

    ctypedef struct dgs_bern_mp_t:
        mpfr_t c
        mpfr_t tmp

    dgs_bern_mp_t *dgs_bern_mp_init(mpfr_t c)
    unsigned long dgs_bern_mp_call(dgs_bern_mp_t *self, gmp_randstate_t state)
    void dgs_bern_mp_clear(dgs_bern_mp_t *self)

    ctypedef struct dgs_bern_exp_mp_t:
        size_t l
        mpfr_t *c
        dgs_bern_mp_t **B

    dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l)
    unsigned long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x, gmp_randstate_t state)
    void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self)

    ctypedef struct dgs_bern_dp_t:
        double p

    dgs_bern_dp_t *dgs_bern_dp_init(double p)
    long dgs_bern_dp_call(dgs_bern_dp_t *self)
    void dgs_bern_dp_clear(dgs_bern_dp_t *self)

    ctypedef struct dgs_bern_exp_dp_t:
        size_t l
        double *p
        dgs_bern_dp_t **B

    dgs_bern_exp_dp_t* dgs_bern_exp_dp_init(double f, size_t l)
    long dgs_bern_exp_dp_call(dgs_bern_exp_dp_t *self, long x)
    void dgs_bern_exp_dp_clear(dgs_bern_exp_dp_t *self)

    ctypedef enum dgs_disc_gauss_alg_t:
        DGS_DISC_GAUSS_UNIFORM_TABLE
        DGS_DISC_GAUSS_UNIFORM_ONLINE
        DGS_DISC_GAUSS_UNIFORM_LOGTABLE
        DGS_DISC_GAUSS_SIGMA2_LOGTABLE

    ctypedef struct dgs_disc_gauss_sigma2p_t:
        dgs_bern_uniform_t *B

    dgs_disc_gauss_sigma2p_t *dgs_disc_gauss_sigma2p_init()
    void dgs_disc_gauss_sigma2p_call(mpz_t rop, dgs_disc_gauss_sigma2p_t *self, gmp_randstate_t state)
    void dgs_disc_gauss_sigma2p_clear(dgs_disc_gauss_sigma2p_t *self)

    ctypedef struct dgs_disc_gauss_mp_t:
        mpfr_t sigma
        mpfr_t c
        size_t tailcut
        dgs_disc_gauss_alg_t algorithm
        dgs_bern_mp_t *B
        void call(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state)

        unsigned long upper_bound
        mpz_t x
        mpfr_t y
        mpfr_t *rho

    dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, mpfr_t c, size_t tailcut, dgs_disc_gauss_alg_t algorithm)
    void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state)
    void dgs_disc_gauss_mp_flush_cache(dgs_disc_gauss_mp_t *self)
    void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self)

    ctypedef struct dgs_disc_gauss_dp_t:
        double sigma
        double c
        size_t tailcut
        dgs_disc_gauss_alg_t algorithm
        dgs_bern_uniform_t *B
        dgs_bern_exp_dp_t *Bexp
        dgs_disc_gauss_sigma2p_t *D2
        long (*call)(dgs_disc_gauss_dp_t *self)

        double f
        long upper_bound
        long two_upper_bound_plus_one
        long k
        double *rho

    dgs_disc_gauss_dp_t *dgs_disc_gauss_dp_init(double sigma, double c, size_t tailcut, dgs_disc_gauss_alg_t algorithm)
    long dgs_disc_gauss_dp_call_uniform_table(dgs_disc_gauss_dp_t *self)
    long dgs_disc_gauss_dp_call_uniform_logtable(dgs_disc_gauss_dp_t *self)
    long dgs_disc_gauss_dp_call_uniform_online(dgs_disc_gauss_dp_t *self)
    long dgs_disc_gauss_dp_call_sigma2_logtable(dgs_disc_gauss_dp_t *self)
    void dgs_disc_gauss_dp_flush_cache(dgs_disc_gauss_dp_t *self)
    void dgs_disc_gauss_dp_clear(dgs_disc_gauss_dp_t *self)
