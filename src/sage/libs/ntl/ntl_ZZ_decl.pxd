include "sage/ext/cdefs.pxi"

cdef extern from "ntl_wrap.h":
    #### ZZ_c
    ctypedef struct ZZ_c "struct ZZ":
        pass

    ctypedef struct vec_ZZ_c "vec_ZZ":
        ZZ_c RawGet(long i)
        ZZ_c *elts()
        long length()

    void del_charstar(char*)

    # Some boiler-plate
    ZZ_c* ZZ_new "New<ZZ>"()
    ZZ_c* ZZ_construct "Construct<ZZ>"(void *mem)
    void ZZ_destruct "Destruct<ZZ>"(ZZ_c *mem)
    void ZZ_delete "Delete<ZZ>"(ZZ_c *mem)
    void ZZ_from_str "_from_str<ZZ>"(ZZ_c* dest, char* s)
    object ZZ_to_PyString "_to_PyString<ZZ>"(ZZ_c *x)
    int ZZ_equal "_equal<ZZ>"(ZZ_c x, ZZ_c y)

    void ZZ_conv_from_int "conv"(ZZ_c x, int i)
    void ZZ_conv_to_int "conv"(int i, ZZ_c x)
    void ZZ_conv_from_long "conv"(ZZ_c x, long l)
    void ZZ_conv_to_long "conv"(long l, ZZ_c x)
    void ZZ_to_mpz(mpz_t output, ZZ_c* x)
    void mpz_to_ZZ(ZZ_c *output, mpz_t x)
    cdef int ZZ_to_int(ZZ_c* x)
    cdef ZZ_c* int_to_ZZ(int value)
    cdef void ZZ_set_from_int(ZZ_c* x, int value)

    long ZZ_sign "sign"(ZZ_c a)
    long ZZ_IsZero "IsZero"(ZZ_c a)
    long ZZ_IsOne "IsOne"(ZZ_c a)
    long ZZ_compare "compare"(ZZ_c a, ZZ_c b)
    void ZZ_add "add"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_add_long "add"(ZZ_c x, ZZ_c a, long b)
    void ZZ_sub "sub"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_sub_long "sub"(ZZ_c x, long a, ZZ_c b)
    void ZZ_negate "negate"(ZZ_c x, ZZ_c a)
    void ZZ_abs "abs"(ZZ_c x, ZZ_c a)

    void ZZ_mul "mul"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_mul_long "mul"( ZZ_c x, ZZ_c a, long b)
    void ZZ_sqr "sqr"(ZZ_c x, ZZ_c a)

    void ZZ_DivRem "DivRem"(ZZ_c q, ZZ_c r, ZZ_c a, ZZ_c b)
    void ZZ_div "div"( ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_rem "rem"( ZZ_c r, ZZ_c a, ZZ_c b)
    long ZZ_DivRem_long "DivRem"(ZZ_c q, ZZ_c a, long b)
    long ZZ_rem_long "rem"(ZZ_c a, long b)
    long ZZ_divide "divide"(ZZ_c q, ZZ_c a, ZZ_c b)
    long ZZ_divide_long "divide"(ZZ_c q, ZZ_c a, long b)
    long ZZ_divide_test "divide"(ZZ_c a, ZZ_c b)
    long ZZ_divide_test_long "divide"(ZZ_c a, long b)

    void ZZ_GCD "GCD"(ZZ_c d, ZZ_c a, ZZ_c b)
    void ZZ_XGCD "XGCD"(ZZ_c d, ZZ_c s, ZZ_c t, ZZ_c a, ZZ_c b)

    void ZZ_AddMod "AddMod"(ZZ_c x, ZZ_c a, ZZ_c b, ZZ_c n)
    void ZZ_SubMod "SubMod"(ZZ_c x, ZZ_c a, ZZ_c b, ZZ_c n)
    void ZZ_NegateMod "NegateMod"(ZZ_c x, ZZ_c a, ZZ_c n)
    void ZZ_MulMod "MulMod"(ZZ_c x, ZZ_c a, ZZ_c b, ZZ_c n)
    void ZZ_SqrMod "SqrMod"(ZZ_c x, ZZ_c a, ZZ_c n)
    void ZZ_InvMod "InvMod"(ZZ_c x, ZZ_c a, ZZ_c n)
    long ZZ_InvModStatus "InvModStatus"(ZZ_c x, ZZ_c a, ZZ_c n)
    void ZZ_PowerMod "PowerMod"(ZZ_c x, ZZ_c a, ZZ_c e, ZZ_c n)
    void ZZ_PowerMod_long "PowerMod"(ZZ_c x, ZZ_c a, long e, ZZ_c n)

    void ZZ_LeftShift "LeftShift"(ZZ_c x, ZZ_c a, long n)
    void ZZ_RightShift "RightShift"(ZZ_c x, ZZ_c a, long n)

    long ZZ_MakeOdd "MakeOdd"(ZZ_c x)
    long ZZ_NumTwos "NumTwos"(ZZ_c x)
    long ZZ_IsOdd "IsOdd"(ZZ_c a)
    long ZZ_NumBits "NumBits"(ZZ_c a)
    long ZZ_bit "bit"(ZZ_c a, long k)
    void ZZ_trunc "trunc"(ZZ_c x, ZZ_c a, long k)
    void ZZ_SetBit "SetBit"(ZZ_c x, long p)
    void ZZ_SwitchBit "SwitchBit"(ZZ_c x, long p)
    void ZZ_weight "weight"(ZZ_c a)
    void ZZ_bit_and "bit_and"(ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_bit_or "bit_or"(ZZ_c x, ZZ_c a, ZZ_c b)
    void ZZ_bit_xor "bit_xor"(ZZ_c x, ZZ_c a, ZZ_c b)

    void ZZ_SetSeed "SetSeed"(ZZ_c s)
    void ZZ_RandomBnd "RandomBnd"(ZZ_c x, ZZ_c n)
    void ZZ_RandomBits "RandomBits"(ZZ_c x, long l)

    long ZZ_CRT "CRT"(ZZ_c a, ZZ_c p, ZZ_c A, ZZ_c P)

    long ZZ_ReconstructRational "ReconstructRational"(ZZ_c a, ZZ_c b, ZZ_c x, ZZ_c m, ZZ_c a_bound, ZZ_c b_bound)

    void ZZ_GenPrime "GenPrime"(ZZ_c n, long l, long err)
    void ZZ_GenGermainPrime "GenGermainPrime"(ZZ_c n, long l, long err)
    long ZZ_ProbPrime "ProbPrime"(ZZ_c n, long NumTrials)
    void ZZ_RandomPrime "RandomPrime"(ZZ_c n, long l, long NumTrials)
    void ZZ_NextPrime "NextPrime"(ZZ_c n, ZZ_c m, long NumTrials)
    long ZZ_MillerWitness "MillerWitness"(ZZ_c n, ZZ_c w)

    void ZZ_power "power"(ZZ_c t, ZZ_c x, long e)
    void ZZ_power2 "power2"(ZZ_c x, long e)

    void ZZ_SqrRoot "SqrRoot"(ZZ_c x, ZZ_c a)

    long ZZ_Jacobi "Jacobi"(ZZ_c a, ZZ_c n)
    void ZZ_SqrRootMod "SqrRootMod"(ZZ_c x, ZZ_c a, ZZ_c n)

    long ZZ_remove(ZZ_c x, ZZ_c a, ZZ_c p) # a la mpz_remove.  Written in ntl_wrap.cpp.

    # Random-number generation
    #void setSeed(ZZ_c* x)
    #ZZ_c* ZZ_randomBnd(ZZ_c* x)
    #ZZ_c* ZZ_randomBits(long n)

