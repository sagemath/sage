# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport GF2X_c, GF2_c, GF2XModulus_c, vec_GF2_c, ZZ_c


cdef extern from "ntlwrap.h":
    long *GF2XHexOutput_c "(&GF2X::HexOutput)" # work-around for Cython bug

    int GF2X_IsOne "IsOne"(GF2X_c x)
    int GF2X_IsZero "IsZero"(GF2X_c x)
    int GF2X_IsX "IsX"(GF2X_c x)

    void GF2X_add "add"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_sub "sub"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_mul "mul"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_negate "NTL::negate"(GF2X_c x, GF2X_c a)
    void GF2X_power "NTL::power"(GF2X_c t, GF2X_c x, long e)
    long GF2X_deg "deg"(GF2X_c x)

    void GF2X_conv_long "conv" (GF2X_c x, long a)
    void GF2X_conv_GF2 "conv" (GF2X_c x, GF2_c a)

    void GF2X_LeftShift "LeftShift"( GF2X_c r, GF2X_c a, long offset)
    void GF2X_RightShift "RightShift"( GF2X_c r, GF2X_c a, long offset)

    void GF2X_trunc "trunc"(GF2X_c r, GF2X_c a, long n)
    void GF2X_DivRem "DivRem"(GF2X_c q, GF2X_c r, GF2X_c a, GF2X_c b)
    void GF2X_div "div" (GF2X_c q, GF2X_c a, GF2X_c b)
    void GF2X_rem "rem" (GF2X_c r, GF2X_c a, GF2X_c b)
    long GF2X_divide "divide"(GF2X_c q, GF2X_c a, GF2X_c b)

    void GF2X_GCD "GCD" (GF2X_c r, GF2X_c a, GF2X_c b)
    void GF2X_XGCD "XGCD" (GF2X_c r, GF2X_c s, GF2X_c t, GF2X_c a, GF2X_c b)

    void GF2XFromBytes(GF2X_c a, unsigned char *p, long n)
    void BytesFromGF2X "BytesFromGF2X" (unsigned char *p, GF2X_c a, long n)

    GF2_c GF2X_coeff "coeff"(GF2X_c a, long i)
    GF2_c GF2X_LeadCoeff "LeadCoeff"(GF2X_c a)
    GF2_c GF2X_ConstTerm "ConstTerm"(GF2X_c a)
    void GF2X_SetCoeff "SetCoeff"(GF2X_c x, long i, GF2_c a)
    void GF2X_SetCoeff_long "SetCoeff"(GF2X_c x, long i, long a)

    GF2X_c GF2X_diff "diff"(GF2X_c a)
    GF2X_c GF2X_reverse "reverse"(GF2X_c a, long hi)

    long GF2X_weight "weight"(GF2X_c a)
    long GF2X_NumBits "NumBits" (GF2X_c a)
    long GF2X_NumBytes "NumBytes"(GF2X_c a)

    void GF2X_MinPolySeq "MinPolySeq" (GF2X_c h, vec_GF2_c a, long m)

    #### GF2XFactoring
    long GF2X_IterIrredTest "IterIrredTest" (GF2X_c f)
    void GF2X_BuildSparseIrred "BuildSparseIrred" (GF2X_c f, long n)
    void GF2X_BuildRandomIrred "BuildRandomIrred" (GF2X_c f, GF2X_c g)
    void GF2X_BuildIrred "BuildIrred" (GF2X_c f, long n)

    #### GF2XModulus_c
    void GF2XModulus_build "build"(GF2XModulus_c F, GF2X_c f) # MUST be called before using the modulus
    long GF2XModulus_deg "deg"(GF2XModulus_c F)

    GF2X_c GF2XModulus_GF2X "GF2X" (GF2XModulus_c m)

    GF2X_c GF2X_IrredPolyMod "IrredPolyMod" (GF2X_c g, GF2XModulus_c F)

    void GF2X_MulMod_pre "MulMod"(GF2X_c x, GF2X_c a, GF2X_c b, GF2XModulus_c F)
    void GF2X_SqrMod_pre "SqrMod"(GF2X_c x, GF2X_c a, GF2XModulus_c F)
    void GF2X_PowerMod_pre "PowerMod"(GF2X_c x, GF2X_c a, ZZ_c e, GF2XModulus_c F)
    void GF2X_PowerMod_long_pre "PowerMod"(GF2X_c x, GF2X_c a, long e, GF2XModulus_c F)
    void GF2X_PowerXMod_pre "PowerXMod"(GF2X_c x, ZZ_c e, GF2XModulus_c F)
    void GF2X_PowerXMod_long_pre "PowerXMod"(GF2X_c x, long e, GF2XModulus_c F)
    void GF2X_PowerXPlusAMod_pre "PowerXPlusAMod"(GF2X_c x, GF2_c a, GF2_c e, GF2XModulus_c F)
    void GF2X_PowerXPlusAMod_long_pre "PowerXPlusAMod"(GF2X_c x, GF2_c a, long e, GF2XModulus_c F)

    # x = g(h) mod f; deg(h) < n
    void GF2X_CompMod "CompMod"(GF2X_c x, GF2X_c g, GF2X_c h, GF2XModulus_c F)
    # xi = gi(h) mod f (i=1,2), deg(h) < n.
    void GF2X_Comp2Mod "Comp2Mod"(GF2X_c x1, GF2X_c x2, GF2X_c g1, GF2X_c g2, GF2X_c h, GF2XModulus_c F)

    # xi = gi(h) mod f (i=1,2,3), deg(h) < n.
    void GF2X_CompMod3 "Comp2Mod"(GF2X_c x1, GF2X_c x2, GF2X_c x3, GF2X_c g1, GF2X_c g2, GF2X_c g3, GF2X_c h, GF2XModulus_c F)
