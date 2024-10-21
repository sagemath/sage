# sage_setup: distribution = sagemath-ntl
# distutils: libraries = gmp mpfr mpfi

from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.mpfi.types cimport *

cdef extern from "mpfi.h":
    # Rounding
    int mpfi_round_prec(mpfi_ptr, mpfr_prec_t prec)

    # Initialization, destruction, and assignment
    # Initializations
    void mpfi_init(mpfi_ptr)
    void mpfi_init2(mpfi_ptr, mpfr_prec_t)
    void mpfi_clear(mpfi_ptr)

    # mpfi bounds have the same precision
    mpfr_prec_t mpfi_get_prec(mpfi_srcptr)
    void mpfi_set_prec(mpfi_ptr, mpfr_prec_t)

    # assignment functions
    int mpfi_set(mpfi_ptr, mpfi_srcptr)
    int mpfi_set_si(mpfi_ptr, long)
    int mpfi_set_ui(mpfi_ptr, unsigned long)
    int mpfi_set_d(mpfi_ptr, double)
    int mpfi_set_z(mpfi_ptr, mpz_t)
    int mpfi_set_q(mpfi_ptr, mpq_t)
    int mpfi_set_fr(mpfi_ptr, mpfr_srcptr)
    int mpfi_set_str(mpfi_ptr, char *, int)

    # combined initialization and assignment functions
    int mpfi_init_set(mpfi_ptr, mpfi_srcptr)
    int mpfi_init_set_si(mpfi_ptr, long)
    int mpfi_init_set_ui(mpfi_ptr, unsigned long)
    int mpfi_init_set_d(mpfi_ptr, double)
    int mpfi_init_set_z(mpfi_ptr, mpz_srcptr)
    int mpfi_init_set_q(mpfi_ptr, mpq_srcptr)
    int mpfi_init_set_fr(mpfi_ptr, mpfr_srcptr)
    int mpfi_init_set_str(mpfi_ptr, char *, int)

    # swapping two intervals
    void mpfi_swap(mpfi_ptr, mpfi_ptr)

    # Various useful interval functions
    # with scalar or interval results

    # absolute diameter
    int mpfi_diam_abs(mpfr_ptr, mpfi_srcptr)
    # relative diameter
    int mpfi_diam_rel(mpfr_ptr, mpfi_srcptr)
    # diameter: relative if the interval does not contain 0
    # absolute otherwise
    int mpfi_diam(mpfr_ptr, mpfi_srcptr)
    # magnitude: the largest absolute value of any element
    int mpfi_mag(mpfr_ptr, mpfi_srcptr)
    # magnitude: the smallest absolute value of any element
    int mpfi_mig(mpfr_ptr, mpfi_srcptr)
    # middle of y
    int mpfi_mid(mpfr_ptr, mpfi_srcptr)
    # picks randomly a point m in y
    void mpfi_alea(mpfr_ptr, mpfi_srcptr)

    # Conversions
    double mpfi_get_d(mpfi_srcptr)
    void mpfi_get_fr(mpfr_ptr, mpfi_srcptr)

    # Basic arithmetic operations

    # arithmetic operations between two interval operands
    int mpfi_add(mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int mpfi_sub(mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int mpfi_mul(mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int mpfi_div(mpfi_ptr, mpfi_srcptr, mpfi_srcptr)

    # arithmetic operations between an interval operand and a double prec. floating-point
    int mpfi_add_d(mpfi_ptr, mpfi_srcptr, double)
    int mpfi_sub_d(mpfi_ptr, mpfi_srcptr, double)
    int mpfi_d_sub(mpfi_ptr, double, mpfi_srcptr)
    int mpfi_mul_d(mpfi_ptr, mpfi_srcptr, double)
    int mpfi_div_d(mpfi_ptr, mpfi_srcptr, double)
    int mpfi_d_div(mpfi_ptr, double, mpfi_srcptr)

    # arithmetic operations between an interval operand and an unsigned long integer
    int mpfi_add_ui(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_sub_ui(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_ui_sub(mpfi_ptr, unsigned long, mpfi_srcptr)
    int mpfi_mul_ui(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_div_ui(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_ui_div(mpfi_ptr, unsigned long, mpfi_srcptr)

    # arithmetic operations between an interval operand and a long integer
    int mpfi_add_si(mpfi_ptr, mpfi_srcptr, long)
    int mpfi_sub_si(mpfi_ptr, mpfi_srcptr, long)
    int mpfi_si_sub(mpfi_ptr, long, mpfi_srcptr)
    int mpfi_mul_si(mpfi_ptr, mpfi_srcptr, long)
    int mpfi_div_si(mpfi_ptr, mpfi_srcptr, long)
    int mpfi_si_div(mpfi_ptr, long, mpfi_srcptr)

    # arithmetic operations between an interval operand and a multiple prec. integer
    int mpfi_add_z(mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int mpfi_sub_z(mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int mpfi_z_sub(mpfi_ptr, mpz_srcptr, mpfi_srcptr)
    int mpfi_mul_z(mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int mpfi_div_z(mpfi_ptr, mpfi_srcptr, mpz_srcptr)
    int mpfi_z_div(mpfi_ptr, mpz_srcptr, mpfi_srcptr)

    # arithmetic operations between an interval operand and a multiple prec. rational
    int mpfi_add_q(mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int mpfi_sub_q(mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int mpfi_q_sub(mpfi_ptr, mpq_srcptr, mpfi_srcptr)
    int mpfi_mul_q(mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int mpfi_div_q(mpfi_ptr, mpfi_srcptr, mpq_srcptr)
    int mpfi_q_div(mpfi_ptr, mpq_srcptr, mpfi_srcptr)

    # arithmetic operations between an interval operand and a mult. prec. floating-pt nb
    int mpfi_add_fr(mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int mpfi_sub_fr(mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int mpfi_fr_sub(mpfi_ptr, mpfr_srcptr, mpfi_srcptr)
    int mpfi_mul_fr(mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int mpfi_div_fr(mpfi_ptr, mpfi_srcptr, mpfr_srcptr)
    int mpfi_fr_div(mpfi_ptr, mpfr_srcptr, mpfi_srcptr)

    # arithmetic operations taking a single interval operand
    int mpfi_neg(mpfi_ptr, mpfi_srcptr)
    int mpfi_sqr(mpfi_ptr, mpfi_srcptr)
    # the inv function generates the whole real interval
    # if 0 is in the interval defining the divisor
    int mpfi_inv(mpfi_ptr, mpfi_srcptr)
    # the sqrt of a (partially) negative interval is a NaN
    int mpfi_sqrt(mpfi_ptr, mpfi_srcptr)
    # the first interval contains the absolute values of
    # every element of the second interval
    int mpfi_abs(mpfi_ptr, mpfi_srcptr)

    # various operations
    int mpfi_mul_2exp(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_mul_2ui(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_mul_2si(mpfi_ptr, mpfi_srcptr, long)
    int mpfi_div_2exp(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_div_2ui(mpfi_ptr, mpfi_srcptr, unsigned long)
    int mpfi_div_2si(mpfi_ptr, mpfi_srcptr, long)

    # Special functions
    int mpfi_log(mpfi_ptr, mpfi_srcptr)
    int mpfi_exp(mpfi_ptr, mpfi_srcptr)
    int mpfi_exp2(mpfi_ptr, mpfi_srcptr)

    int mpfi_cos(mpfi_ptr, mpfi_srcptr)
    int mpfi_sin(mpfi_ptr, mpfi_srcptr)
    int mpfi_tan(mpfi_ptr, mpfi_srcptr)
    int mpfi_acos(mpfi_ptr, mpfi_srcptr)
    int mpfi_asin(mpfi_ptr, mpfi_srcptr)
    int mpfi_atan(mpfi_ptr, mpfi_srcptr)

    int mpfi_cosh(mpfi_ptr, mpfi_srcptr)
    int mpfi_sinh(mpfi_ptr, mpfi_srcptr)
    int mpfi_tanh(mpfi_ptr, mpfi_srcptr)
    int mpfi_acosh(mpfi_ptr, mpfi_srcptr)
    int mpfi_asinh(mpfi_ptr, mpfi_srcptr)
    int mpfi_atanh(mpfi_ptr, mpfi_srcptr)

    int mpfi_log1p(mpfi_ptr, mpfi_srcptr)
    int mpfi_expm1(mpfi_ptr, mpfi_srcptr)

    int mpfi_log2(mpfi_ptr, mpfi_srcptr)
    int mpfi_log10(mpfi_ptr, mpfi_srcptr)

    int mpfi_const_log2(mpfi_ptr)
    int mpfi_const_pi(mpfi_ptr)
    int mpfi_const_euler(mpfi_ptr)

    # Comparison functions
    # Warning: the meaning of interval comparison is not clearly defined
    # customizable comparison functions

    int (*mpfi_cmp) (mpfi_srcptr, mpfi_srcptr)
    int (*mpfi_cmp_d) (mpfi_srcptr, double)
    int (*mpfi_cmp_ui) (mpfi_srcptr, unsigned long)
    int (*mpfi_cmp_si) (mpfi_srcptr, long)
    int (*mpfi_cmp_z) (mpfi_srcptr, mpz_srcptr)
    int (*mpfi_cmp_q) (mpfi_srcptr, mpq_srcptr)
    int (*mpfi_cmp_fr) (mpfi_srcptr, mpfr_srcptr)

    bint (*mpfi_is_pos) (mpfi_srcptr)
    bint (*mpfi_is_nonneg) (mpfi_srcptr)
    bint (*mpfi_is_neg) (mpfi_srcptr)
    bint (*mpfi_is_nonpos) (mpfi_srcptr)
    bint (*mpfi_is_zero) (mpfi_srcptr)
    bint (*mpfi_is_strictly_pos) (mpfi_srcptr)
    bint (*mpfi_is_strictly_neg) (mpfi_srcptr)

    bint mpfi_has_zero(mpfi_srcptr)

    bint mpfi_nan_p(mpfi_srcptr)
    bint mpfi_inf_p(mpfi_srcptr)
    bint mpfi_bounded_p(mpfi_srcptr)

    # Interval manipulation

    # operations related to the internal representation by endpoints

    # get left or right bound of the interval defined by the
    # second argument and put the result in the first one
    int mpfi_get_left(mpfr_ptr, mpfi_srcptr)
    int mpfi_get_right(mpfr_ptr, mpfi_srcptr)

    int mpfi_revert_if_needed(mpfi_ptr)

    # Set operations on intervals
    # "Convex hulls"
    # extends the interval defined by the first argument
    # so that it contains the second one

    int mpfi_put(mpfi_ptr, mpfi_srcptr)
    int mpfi_put_d(mpfi_ptr, double)
    int mpfi_put_si(mpfi_ptr, long)
    int mpfi_put_ui(mpfi_ptr, unsigned long)
    int mpfi_put_z(mpfi_ptr, mpz_srcptr)
    int mpfi_put_q(mpfi_ptr, mpq_srcptr)
    int mpfi_put_fr(mpfi_ptr, mpfr_srcptr)

    # builds an interval whose left bound is the lower (round -infty)
    # than the second argument and the right bound is greater
    # (round +infty) than the third one


    int mpfi_interv_d(mpfi_ptr, double,double)
    int mpfi_interv_si(mpfi_ptr, long,long)
    int mpfi_interv_ui(mpfi_ptr, unsigned long,unsigned long)
    int mpfi_interv_z(mpfi_ptr, mpz_srcptr,mpz_srcptr)
    int mpfi_interv_q(mpfi_ptr, mpq_srcptr,mpq_srcptr)
    int mpfi_interv_fr(mpfi_ptr, mpfr_srcptr,mpfr_srcptr)

    # Inclusion tests
    # tests if the first argument is inside the interval
    # defined by the second one
    bint mpfi_is_strictly_inside(mpfi_srcptr, mpfi_srcptr)
    bint mpfi_is_inside(mpfi_srcptr, mpfi_srcptr)
    bint mpfi_is_inside_d(double, mpfi_srcptr)
    bint mpfi_is_inside_ui(unsigned long, mpfi_srcptr)
    bint mpfi_is_inside_si(long, mpfi_srcptr)
    bint mpfi_is_inside_z(mpz_srcptr, mpfi_srcptr)
    bint mpfi_is_inside_q(mpq_srcptr, mpfi_srcptr)
    bint mpfi_is_inside_fr(mpfr_srcptr, mpfi_srcptr)

    # set operations
    bint mpfi_is_empty(mpfi_srcptr)
    int mpfi_intersect(mpfi_ptr, mpfi_srcptr, mpfi_srcptr)
    int mpfi_union(mpfi_ptr, mpfi_srcptr, mpfi_srcptr)

    # Miscellaneous

    # adds the second argument to the right bound of the first one
    # and subtracts the second argument to the left bound of
    # the first one
    int mpfi_increase(mpfi_ptr, mpfr_srcptr)
    # keeps the same center and multiply the radius by 2*(1+fact)
    int mpfi_blow(mpfi_ptr, mpfi_srcptr, double)
    # splits the interval into 2 halves
    int mpfi_bisect(mpfi_ptr, mpfi_ptr, mpfi_srcptr)

    char * mpfi_get_version()

    # Error handling
    void mpfi_reset_error()
    void mpfi_set_error(int)
    int mpfi_is_error()

    ctypedef enum mpfi_flags_exact:
        MPFI_FLAGS_BOTH_ENDPOINTS_EXACT
        MPFI_FLAGS_LEFT_ENDPOINT_INEXACT
        MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT
        MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT

    bint MPFI_BOTH_ARE_EXACT(mpfi_flags_exact)
    bint MPFI_LEFT_IS_INEXACT(mpfi_flags_exact)
    bint MPFI_RIGHT_IS_INEXACT(mpfi_flags_exact)
    bint MPFI_BOTH_ARE_INEXACT(mpfi_flags_exact)

    mpfi_flags_exact MPFI_REVERT_INEXACT_FLAGS(mpfi_flags_exact)

    bint MPFI_NAN_P(mpfi_srcptr)
    bint MPFI_INF_P(mpfi_srcptr)
    bint MPFI_IS_ZERO(mpfi_srcptr)

    void MPFI_CLEAR(mpfi_ptr)
