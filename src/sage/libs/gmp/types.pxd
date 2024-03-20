# sage_setup: distribution = sagemath-objects
from libc.stdio cimport FILE

cdef extern from "gmp.h":
    # GMP's configuration of how many bits are stuffed into a limb
    cdef unsigned int GMP_LIMB_BITS
    cdef int mp_bits_per_limb

    # Underlying typedefs
    ctypedef unsigned long mp_limb_t
    ctypedef long mp_limb_signed_t
    ctypedef unsigned long mp_bitcnt_t
    ctypedef long mp_size_t
    ctypedef long mp_exp_t

    ctypedef mp_limb_t* mp_ptr
    ctypedef mp_limb_t* mp_srcptr

    # This internal structure is not guaranteed to stay the same with
    # future releases of GMP or MPIR.
    ctypedef struct __mpz_struct:
        int _mp_alloc
        int _mp_size
        mp_ptr _mp_d

    ctypedef struct __mpq_struct:
        __mpz_struct _mp_num
        __mpz_struct _mp_den

    ctypedef struct __mpf_struct:
        int _mp_prec
        int _mp_size
        mp_exp_t _mp_exp
        mp_limb_t *_mp_d

    ctypedef struct __gmp_randstate_struct:
        pass

    # User facing types
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpq_struct mpq_t[1]
    ctypedef __mpf_struct mpf_t[1]
    ctypedef __gmp_randstate_struct gmp_randstate_t[1]

    ctypedef __mpz_struct *mpz_ptr
    ctypedef __mpq_struct *mpq_ptr
    ctypedef __mpf_struct *mpf_ptr
    ctypedef __mpz_struct *mpz_srcptr
    ctypedef __mpq_struct *mpq_srcptr
    ctypedef __mpf_struct *mpf_srcptr

    # Cython doesn't want to take the address of an mpz_t
    cdef mpz_t* address_of_mpz "&"(mpz_t x)
