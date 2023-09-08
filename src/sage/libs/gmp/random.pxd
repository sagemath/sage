# sage_setup: distribution = sagemath-objects
# distutils: libraries = gmp

from sage.libs.gmp.types cimport *

cdef extern from "gmp.h":

    ### Random Number Functions ###

    # Random State Initialization
    void gmp_randinit_default (gmp_randstate_t state)
    int gmp_randinit_mt (gmp_randstate_t state)
    void gmp_randinit_lc_2exp (gmp_randstate_t state, mpz_t a, unsigned long c, unsigned long m2exp)
    int gmp_randinit_lc_2exp_size (gmp_randstate_t state, unsigned long size)
    int gmp_randinit_set (gmp_randstate_t rop, gmp_randstate_t op)
    # void gmp_randinit (gmp_randstate_t state, gmp_randalg_t alg, ...)
    void gmp_randclear (gmp_randstate_t state)

    # Random State Seeding
    void gmp_randseed (gmp_randstate_t state, mpz_t seed)
    void gmp_randseed_ui (gmp_randstate_t state, unsigned long int seed)

    # Random State Miscellaneous
    unsigned long gmp_urandomb_ui (gmp_randstate_t state, unsigned long n)
    unsigned long gmp_urandomm_ui (gmp_randstate_t state, unsigned long n)
