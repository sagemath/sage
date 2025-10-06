# distutils: depends = NTL/ZZ.h

from sage.libs.ntl.types cimport (ZZ_c, ZZX_c, ZZ_p_c, vec_ZZ_p_c, ZZ_pContext_c,
        ZZ_pX_c, ZZ_pX_Modulus_c, ZZ_pX_Multiplier_c)


cdef extern from "ntlwrap.h":
    long ZZ_pX_IsZero "IsZero"(ZZ_pX_c a)
    long ZZ_pX_IsOne "IsOne"(ZZ_pX_c a)

    void ZZ_pX_add "add"( ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_add_long "add"(ZZ_pX_c x, ZZ_pX_c a, long b)
    void ZZ_pX_sub "sub"( ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_sub_long "sub"(ZZ_pX_c x, long a, ZZ_pX_c b)
    void ZZ_pX_negate "NTL::negate"(ZZ_pX_c x, ZZ_pX_c a)

    void ZZ_pX_mul "mul"( ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_mul_long "mul"( ZZ_pX_c x, ZZ_pX_c a, long b)
    void ZZ_pX_mul_ZZ_p "mul"( ZZ_pX_c x, ZZ_pX_c a, ZZ_p_c b)
    void ZZ_pX_rmul "mul"( ZZ_pX_c x, ZZ_pX_c a, ZZ_p_c b)
    void ZZ_pX_sqr "sqr"( ZZ_pX_c x, ZZ_pX_c a)
    long ZZ_pX_power "NTL::power"( ZZ_pX_c x, ZZ_pX_c a, long e)

    void ZZ_pX_LeftShift "LeftShift"(ZZ_pX_c x, ZZ_pX_c a, long n)
    void ZZ_pX_RightShift "RightShift"(ZZ_pX_c x, ZZ_pX_c a, long n)

    void ZZ_pX_DivRem "DivRem"(ZZ_pX_c q, ZZ_pX_c r, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_div_ZZ "div"( ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_div "div"( ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_c b)
    long ZZ_pX_divide "divide"( ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_rem "rem"(ZZ_pX_c r, ZZ_pX_c a, ZZ_pX_c b)

    void ZZ_pX_GCD "GCD"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_PlainXGCD "PlainXGCD"(ZZ_pX_c d, ZZ_pX_c s, ZZ_pX_c t, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_XGCD "XGCD"(ZZ_pX_c d, ZZ_pX_c s, ZZ_pX_c t, ZZ_pX_c a, ZZ_pX_c b)

    long ZZ_pX_deg "deg"( ZZ_pX_c x )
    ZZ_p_c ZZ_pX_coeff "coeff"(ZZ_pX_c a, long i)
    ZZ_p_c ZZ_pX_LeadCoeff "LeadCoeff"(ZZ_pX_c x)
    ZZ_p_c ZZ_pX_ConstTerm "ConstTerm"(ZZ_pX_c x)
    void ZZ_pX_SetCoeff "SetCoeff"(ZZ_pX_c x, long i, ZZ_p_c a)
    void ZZ_pX_SetCoeff_long "SetCoeff"(ZZ_pX_c x, long i, long a)
    void ZZ_pX_SetX "SetX"(ZZ_pX_c x)
    long ZZ_pX_IsX "IsX"(ZZ_pX_c a)
    void ZZ_pX_diff "diff"(ZZ_pX_c x, ZZ_pX_c a)
    void ZZ_pX_MakeMonic "MakeMonic"(ZZ_pX_c x)
    void ZZ_pX_reverse_hi "reverse"(ZZ_pX_c x, ZZ_pX_c a, long hi)
    void ZZ_pX_reverse "reverse"(ZZ_pX_c x, ZZ_pX_c a)
    void ZZ_pX_VectorCopy "VectorCopy"(vec_ZZ_p_c x, ZZ_pX_c a, long n)

    void ZZ_pX_random "random"(ZZ_pX_c x, long n)

    void ZZ_pX_BuildFromRoots "BuildFromRoots"(ZZ_pX_c x, vec_ZZ_p_c a)
    void ZZ_pX_eval "eval"(ZZ_p_c fa, ZZ_pX_c f, ZZ_p_c a)
    void ZZ_pX_eval_vec "eval"(vec_ZZ_p_c fa, ZZ_pX_c f, vec_ZZ_p_c a)
    void ZZ_pX_interpolate "interpolate"(ZZ_pX_c f, vec_ZZ_p_c a, vec_ZZ_p_c b)

    void ZZ_pX_trunc "trunc"(ZZ_pX_c x, ZZ_pX_c a, long n)
    void ZZ_pX_MulTrunc "MulTrunc"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b, long n)
    void ZZ_pX_SqrTrunc "SqrTrunc"(ZZ_pX_c x, ZZ_pX_c a, long n)
    void ZZ_pX_InvTrunc "InvTrunc"(ZZ_pX_c x, ZZ_pX_c a, long n)

    void ZZ_pX_MulMod "MulMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b, ZZ_pX_c f)
    void ZZ_pX_SqrMod "SqrMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_MulByXMod "MulByXMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_InvMod "InvMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)
    long ZZ_pX_InvModStatus "InvModStatus"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c f)

    void ZZ_pX_Modulus_build "build"(ZZ_pX_Modulus_c F, ZZ_pX_c f) # MUST be called before using the modulus
    long ZZ_pX_Modulus_deg "deg"(ZZ_pX_Modulus_c F)

    void ZZ_pX_MulMod_pre "MulMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_c b, ZZ_pX_Modulus_c F)
    void ZZ_pX_SqrMod_pre "SqrMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerMod_pre "PowerMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_c e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerMod_long_pre "PowerMod"(ZZ_pX_c x, ZZ_pX_c a, long e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXMod_pre "PowerXMod"(ZZ_pX_c x, ZZ_c e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXMod_long_pre "PowerXMod"(ZZ_pX_c x, long e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXPlusAMod_pre "PowerXPlusAMod"(ZZ_pX_c x, ZZ_p_c a, ZZ_c e, ZZ_pX_Modulus_c F)
    void ZZ_pX_PowerXPlusAMod_long_pre "PowerXPlusAMod"(ZZ_pX_c x, ZZ_p_c a, long e, ZZ_pX_Modulus_c F)
    void ZZ_pX_rem_pre "rem"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_DivRem_pre "DivRem"(ZZ_pX_c q, ZZ_pX_c r, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_div_pre "div"(ZZ_pX_c q, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_InvMod_pre "InvMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)

    void ZZ_pX_Multiplier_build "build"(ZZ_pX_Multiplier_c F, ZZ_pX_c b, ZZ_pX_Modulus_c F) # MUST be called before using the multiplier
    void ZZ_pX_MulMod_premul "MulMod"(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Multiplier_c B, ZZ_pX_Modulus_c F)

    void ZZ_pX_CompMod "CompMod"(ZZ_pX_c x, ZZ_pX_c g, ZZ_pX_c h, ZZ_pX_Modulus_c F)
    # Add other composition functions here
    # Add power projection routines
    # Add minimum polynomials of recurrence sequences
    void ZZ_pX_MinPolyMod "MinPolyMod"(ZZ_pX_c h, ZZ_pX_c g, ZZ_pX_c f)
    void ZZ_pX_MinPolyMod_pre "MinPolyMod"(ZZ_pX_c h, ZZ_pX_c g, ZZ_pX_Modulus_c F)

    void ZZ_pX_TraceMod "TraceMod"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_TraceMod_pre "TraceMod"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F)
    void ZZ_pX_TraceVec "TraceVec"(vec_ZZ_p_c S, ZZ_pX_c f)
    void ZZ_pX_NormMod "NormMod"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_c f)
    void ZZ_pX_resultant "resultant"(ZZ_p_c x, ZZ_pX_c a, ZZ_pX_c b)
    void ZZ_pX_CharPolyMod "CharPolyMod"(ZZ_pX_c g, ZZ_pX_c a, ZZ_pX_c f)

    void ZZ_pX_clear "clear"(ZZ_pX_c x)
    void ZZ_pX_set "set"(ZZ_pX_c x)

    void ZZ_pX_to_ZZX "conv"(ZZX_c x, ZZ_pX_c a)
    void ZZX_to_ZZ_pX "conv"(ZZ_pX_c x, ZZX_c a)


cdef extern from "ntlwrap_impl.h":
    char* ZZ_pX_trace_list(ZZ_pX_c* x)
    void ZZ_pX_factor(ZZ_pX_c*** v, long** e, long* n, ZZ_pX_c* x, long verbose)
    void ZZ_pX_linear_roots(ZZ_p_c*** v, long* n, ZZ_pX_c* x)

    # The following are all used for padics.
    void ZZ_pX_conv_modulus(ZZ_pX_c fout, ZZ_pX_c fin, ZZ_pContext_c c)
    void ZZ_pX_min_val_coeff(long valuation, long index, ZZ_pX_c f, ZZ_c p)
    long ZZ_pX_get_val_coeff(ZZ_pX_c f, ZZ_c p, long i)
    void ZZ_pX_left_pshift(ZZ_pX_c x, ZZ_pX_c a, ZZ_c pn, ZZ_pContext_c c)
    void ZZ_pX_right_pshift(ZZ_pX_c x, ZZ_pX_c a, ZZ_c pn, ZZ_pContext_c c)
    void ZZ_pX_InvMod_newton_unram(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F, ZZ_pContext_c cpn, ZZ_pContext_c cp)
    void ZZ_pX_InvMod_newton_ram(ZZ_pX_c x, ZZ_pX_c a, ZZ_pX_Modulus_c F, ZZ_pContext_c cpn)
