# sage_setup: distribution = sagemath-ntl
# distutils: depends = NTL/ZZ.h

from sage.libs.ntl.types cimport (ZZ_c, ZZ_p_c, ZZ_pContext_c, ZZ_pE_c, vec_ZZ_p_c,
        vec_ZZ_pE_c, ZZ_pEX_c, ZZ_pEX_Modulus_c)


cdef extern from "ntlwrap.h":
    long ZZ_pEX_IsZero "IsZero"(ZZ_pEX_c a)
    long ZZ_pEX_IsOne "IsOne"(ZZ_pEX_c a)

    void ZZ_pEX_add "add"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_add_ZZ_p "add"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_add_ZZ_pE "add"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_sub "sub"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_sub_ZZ_p "sub"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_sub_ZZ_pE "sub"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_negate "NTL::negate"(ZZ_pEX_c x, ZZ_pEX_c a)

    void ZZ_pEX_mul "mul"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_mul_long "mul"( ZZ_pEX_c x, ZZ_pEX_c a, long b)
    void ZZ_pEX_mul_ZZ_p "mul"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_mul_ZZ_pE "mul"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_sqr "sqr"( ZZ_pEX_c x, ZZ_pEX_c a)
    void ZZ_pEX_power "NTL::power"( ZZ_pEX_c x, ZZ_pEX_c a, long e)

    void ZZ_pEX_LeftShift "LeftShift"(ZZ_pEX_c x, ZZ_pEX_c a, long n)
    void ZZ_pEX_RightShift "RightShift"(ZZ_pEX_c x, ZZ_pEX_c a, long n)

    void ZZ_pEX_DivRem "DivRem"(ZZ_pEX_c q, ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_div_ZZ_pEX "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_div_ZZ_pE "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_pE_c b)
    void ZZ_pEX_div_ZZ_p "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_p_c b)
    void ZZ_pEX_div_long "div"( ZZ_pEX_c x, ZZ_pEX_c a, long b)
    void ZZ_pEX_rem "rem"(ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_c b)
    long ZZ_pEX_divide "divide"( ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b)

    void ZZ_pEX_GCD "GCD"(ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_c b)
    void ZZ_pEX_XGCD "XGCD"(ZZ_pEX_c d, ZZ_pEX_c s, ZZ_pEX_c t, ZZ_pEX_c a, ZZ_pEX_c b)

    long ZZ_pEX_deg "deg"( ZZ_pEX_c x )
    ZZ_pE_c ZZ_pEX_coeff "coeff"(ZZ_pEX_c a, long i)
    ZZ_pE_c ZZ_pEX_LeadCoeff "LeadCoeff"(ZZ_pEX_c a)
    void ZZ_pEX_SetCoeff "SetCoeff"(ZZ_pEX_c x, long i, ZZ_pE_c a)
    void ZZ_pEX_SetCoeff_ZZ_p "SetCoeff"(ZZ_pEX_c x, long i, ZZ_p_c a)
    void ZZ_pEX_SetCoeff_long "SetCoeff"(ZZ_pEX_c x, long i, long a)
    void ZZ_pEX_SetCoeff_one "SetCoeff"(ZZ_pEX_c x, long i)
    void ZZ_pEX_SetX "SetX"(ZZ_pEX_c x)
    long ZZ_pEX_IsX "IsX"(ZZ_pEX_c x)
    void ZZ_pEX_diff "diff"(ZZ_pEX_c x, ZZ_pEX_c a)
    void ZZ_pEX_MakeMonic "MakeMonic"(ZZ_pEX_c x)
    void ZZ_pEX_reverse_hi "reverse"(ZZ_pEX_c x, ZZ_pEX_c a, long hi)
    void ZZ_pEX_reverse "reverse"(ZZ_pEX_c x, ZZ_pEX_c a)

    void ZZ_pEX_random "random"(ZZ_pEX_c x, long n)

    void ZZ_pEX_BuildFromRoots "BuildFromRoots"(ZZ_pEX_c x, vec_ZZ_pE_c a)
    void ZZ_pEX_eval "eval"(ZZ_pE_c b, ZZ_pEX_c f, ZZ_pE_c a)
    void ZZ_pEX_eval_vec "eval"(vec_ZZ_pE_c b, ZZ_pEX_c f, vec_ZZ_pE_c a)
    void ZZ_pEX_interpolate "interpolate"(ZZ_pEX_c f, vec_ZZ_pE_c a, vec_ZZ_pE_c b)

    void ZZ_pEX_trunc "trunc"(ZZ_pEX_c x, ZZ_pEX_c a, long n)
    void ZZ_pEX_MulTrunc "MulTrunc"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b, long n)
    void ZZ_pEX_SqrTrunc "SqrTrunc"(ZZ_pEX_c x, ZZ_pEX_c a, long n)
    void ZZ_pEX_InvTrunc "InvTrunc"(ZZ_pEX_c x, ZZ_pEX_c a, long n)

    void ZZ_pEX_MulMod "MulMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b, ZZ_pEX_c f)
    void ZZ_pEX_SqrMod "SqrMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_MulByXMod "MulByXMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_InvMod "InvMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    long ZZ_pEX_InvModStatus "InvModStatus"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c f)

    void ZZ_pEX_Modulus_build "build"(ZZ_pEX_Modulus_c F, ZZ_pEX_c f)
    long ZZ_pEX_Modulus_deg "deg"(ZZ_pEX_Modulus_c F)

    void ZZ_pEX_MulMod_pre "MulMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_c b, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_SqrMod_pre "SqrMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerMod_pre "PowerMod"(ZZ_pEX_c x, ZZ_pEX_c a, long e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerMod_ZZ_pre "PowerMod"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_c e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerXMod_pre "PowerXMod"(ZZ_pEX_c x, long e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_PowerXMod_ZZ_pre "PowerXMod"(ZZ_pEX_c x, ZZ_c e, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_rem_pre "rem"(ZZ_pEX_c x, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_DivRem_pre "DivRem"(ZZ_pEX_c q, ZZ_pEX_c r, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_div_pre "div"(ZZ_pEX_c q, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)

    void ZZ_pEX_CompMod "CompMod"(ZZ_pEX_c x, ZZ_pEX_c f, ZZ_pEX_c g, ZZ_pEX_Modulus_c F)

    void ZZ_pEX_MinPolyMod "MinPolyMod"(ZZ_pEX_c h, ZZ_pEX_c g, ZZ_pEX_c f)
    void ZZ_pEX_MinPolyMod_pre "MinPolyMod"(ZZ_pEX_c h, ZZ_pEX_c g, ZZ_pEX_Modulus_c F)

    void ZZ_pEX_TraceMod "TraceMod"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_TraceMod_pre "TraceMod"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_Modulus_c F)
    void ZZ_pEX_TraceVec "TraceVec"(vec_ZZ_pE_c x, ZZ_pEX_c f)
    void ZZ_pEX_NormMod "NormMod"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_c f)
    void ZZ_pEX_resultant "resultant"(ZZ_pE_c x, ZZ_pEX_c a, ZZ_pEX_c b)

    void ZZ_pEX_clear "clear"(ZZ_pEX_c x)
    void ZZ_pEX_set "set"(ZZ_pEX_c x)

    long ZZ_pEX_IterIrredTest "IterIrredTest"(ZZ_pEX_c x)
    long ZZ_pEX_DetIrredTest "DetIrredTest"(ZZ_pEX_c x)
    long ZZ_pEX_ProbIrredTest "ProbIrredTest"(ZZ_pEX_c x, long iter)


cdef extern from "ntlwrap_impl.h":
    void ZZ_pEX_conv_modulus(ZZ_pEX_c fout, ZZ_pEX_c fin, ZZ_pContext_c c)
