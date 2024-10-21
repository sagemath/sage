# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport GF2_c


cdef extern from "ntlwrap.h":
    int GF2_IsOne "IsOne"(GF2_c x)
    int GF2_IsZero "IsZero"(GF2_c x)

    void GF2_add "add"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_sub "sub"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_mul "mul"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_div "div"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_negate "NTL::negate"(GF2_c x, GF2_c a)
    void GF2_power "NTL::power"(GF2_c t, GF2_c x, long e)
    long GF2_deg "deg"(GF2_c x)

    void GF2_conv_long "conv" (GF2_c x, long i)
    long GF2_conv_to_long "rep" (GF2_c x)
