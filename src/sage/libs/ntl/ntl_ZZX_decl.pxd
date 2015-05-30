include "sage/ext/cdefs.pxi"

from sage.libs.ntl.ntl_ZZ_decl cimport vec_ZZ_c
from sage.libs.ntl.ntl_ZZ_decl cimport ZZ_c

cdef extern from "ntl_wrap.h":
    # really, this is from NTL/ZZX.h
    ctypedef struct ZZX_c "struct ZZX":
        vec_ZZ_c rep

    ctypedef struct pair_ZZX_long_c "pair_ZZX_long":
        ZZX_c a
        long b

    ctypedef struct vec_pair_ZZX_long_c "vec_pair_ZZX_long":
        pair_ZZX_long_c RawGet(long i)
        long length()

    # Some boiler-plate
    ZZX_c* ZZX_new "New<ZZX>"()
    ZZX_c* ZZX_construct "Construct<ZZX>"(void *mem)
    void ZZX_destruct "Destruct<ZZX>"(ZZX_c *mem)
    void ZZX_swap "swap"(ZZX_c x, ZZX_c y)
    void ZZX_delete "Delete<ZZX>"(ZZX_c *mem)
    void ZZX_from_str "_from_str<ZZX>"(ZZX_c* dest, char* s)
    object ZZX_to_PyString "_to_PyString<ZZX>"(ZZX_c *x)
    int ZZX_equal "_equal<ZZX>"(ZZX_c x, ZZX_c y)

    void ZZX_PseudoRem "PseudoRem"(ZZX_c x, ZZX_c a, ZZX_c b)
    ZZ_c ZZX_LeadCoeff "LeadCoeff"(ZZX_c x)
    ZZ_c ZZX_ConstTerm "ConstTerm"(ZZX_c x)
    ZZ_c ZZX_coeff "coeff"(ZZX_c a, long i)
    void ZZX_SetCoeff "SetCoeff"(ZZX_c x, long i, ZZ_c a)
    void ZZX_SetCoeff_long "SetCoeff"(ZZX_c x, long i, long a)
    long IsZero_ZZX "IsZero"(ZZX_c a)
    # f must be monic!
    void ZZX_MulMod "MulMod"(ZZX_c x, ZZX_c a, ZZX_c b, ZZX_c f)
    void ZZX_mul_long "mul"( ZZX_c x, ZZX_c a, long b)
    void ZZX_mul_ZZ "mul"( ZZX_c x, ZZX_c a, ZZ_c b)
    void ZZX_mul "mul"( ZZX_c x, ZZX_c a, ZZX_c b)
    void ZZX_add "add"( ZZX_c x, ZZX_c a, ZZX_c b)
    void ZZX_sub "sub"( ZZX_c x, ZZX_c a, ZZX_c b)
    void ZZX_negate "negate"( ZZX_c x, ZZX_c a)
    void ZZX_div_ZZ "div"( ZZX_c x, ZZX_c a, ZZ_c b)
    long ZZX_deg "deg"( ZZX_c x )
    void ZZX_rem "rem"(ZZX_c r, ZZX_c a, ZZX_c b)
    void ZZX_XGCD "XGCD"(ZZ_c r, ZZX_c s, ZZX_c t, ZZX_c a, ZZX_c b, long deterministic)
    void ZZX_content "content"(ZZ_c d, ZZX_c f)
    void ZZX_factor "factor"(ZZ_c c, vec_pair_ZZX_long_c factors, ZZX_c f, long verbose, long bnd)

    void ZZX_squarefree_decomposition(ZZX_c*** v, long** e, long* n, ZZX_c* x)

    char* ZZX_repr(ZZX_c* x)
    ## for cleaning up after ZZX_repr:
    void cpp_delete_array "delete []"(char *str)
    ZZX_c* ZZX_copy(ZZX_c* x)
    ZZX_c* ZZX_div(ZZX_c* x, ZZX_c* y, int* divisible)
    void ZZX_quo_rem(ZZX_c* x, ZZX_c* other, ZZX_c** r, ZZX_c** q)
    ZZX_c* ZZX_square(ZZX_c* x)
    int ZZX_IsZero "IsZero"(ZZX_c x)
    int ZZX_IsOne "IsOne"(ZZX_c x)
    int ZZX_is_monic(ZZX_c* x)
    ZZX_c* ZZX_neg(ZZX_c* x)
    ZZX_c* ZZX_left_shift(ZZX_c* x, long n)
    ZZX_c* ZZX_right_shift(ZZX_c* x, long n)
    ZZX_c* ZZX_primitive_part(ZZX_c* x)
    void ZZX_pseudo_quo_rem(ZZX_c* x, ZZX_c* y, ZZX_c** r, ZZX_c** q)
    ZZX_c* ZZX_gcd(ZZX_c* x, ZZX_c* y)
    ZZX_c* ZZX_xgcd(ZZX_c* x, ZZX_c* y, ZZ_c** r, ZZX_c** s, ZZX_c** t, int proof)
    void ZZX_set_x(ZZX_c* x)
    int ZZX_is_x(ZZX_c* x)
    ZZX_c* ZZX_derivative(ZZX_c* x)
    ZZX_c* ZZX_reverse(ZZX_c* x)
    ZZX_c* ZZX_reverse_hi(ZZX_c* x, long hi)
    ZZX_c* ZZX_truncate(ZZX_c* x, long m)
    ZZX_c* ZZX_multiply_and_truncate(ZZX_c* x, ZZX_c* y, long m)
    ZZX_c* ZZX_square_and_truncate(ZZX_c* x, long m)
    ZZX_c* ZZX_invert_and_truncate(ZZX_c* x, long m)
    ZZX_c* ZZX_multiply_mod(ZZX_c* x, ZZX_c* y,  ZZX_c* modulus)
    ZZ_c* ZZX_trace_mod(ZZX_c* x, ZZX_c* y)
    char* ZZX_trace_list(ZZX_c* x)
    ZZ_c* ZZX_resultant(ZZX_c* x, ZZX_c* y, int proof)
    ZZ_c* ZZX_norm_mod(ZZX_c* x, ZZX_c* y, int proof)
    ZZ_c* ZZX_discriminant(ZZX_c* x, int proof)
    ZZ_c* ZZX_polyeval(ZZX_c* x, ZZ_c* a)
    ZZX_c* ZZX_charpoly_mod(ZZX_c* x, ZZX_c* y, int proof)
    ZZX_c* ZZX_minpoly_mod(ZZX_c* x, ZZX_c* y)
    void ZZX_clear(ZZX_c* x)
    void ZZX_preallocate_space(ZZX_c* x, long n)
    void ZZX_getitem_as_mpz(mpz_t output, ZZX_c* x, long i)
    cdef void ZZX_setitem_from_int(ZZX_c* x, long i, int value)
    cdef int ZZX_getitem_as_int(ZZX_c* x, long i)
