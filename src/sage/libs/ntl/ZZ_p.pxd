# distutils: depends = NTL/ZZ.h

from .types cimport ZZ_c, ZZ_p_c


cdef extern from "ntlwrap.cpp":
    char* ZZ_p_to_str(ZZ_p_c* x)
    void ZZ_p_add "add"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_sub "sub"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_mul "mul"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_mul_long "mul"( ZZ_p_c x, ZZ_p_c a, long b)
    void ZZ_p_div "div"( ZZ_p_c x, ZZ_p_c a, ZZ_p_c b)
    void ZZ_p_negate "NTL::negate"(ZZ_p_c x, ZZ_p_c a)
    void ZZ_p_power "NTL::power"(ZZ_p_c t, ZZ_p_c x, long e)
    int ZZ_p_IsOne "IsOne"(ZZ_p_c x)
    int ZZ_p_IsZero "IsZero"(ZZ_p_c x)
    ZZ_c ZZ_p_rep "rep"(ZZ_p_c z)
    ZZ_p_c* ZZ_p_neg(ZZ_p_c* x)
    void ntl_ZZ_set_modulus(ZZ_c* x)
    void ZZ_p_inv "inv"(ZZ_p_c r, ZZ_p_c x)
    void ZZ_p_random "random"(ZZ_p_c r)
    ZZ_p_c long_to_ZZ_p "to_ZZ_p"(long i)
    ZZ_p_c ZZ_to_ZZ_p "to_ZZ_p"(ZZ_c i)
    int ZZ_p_to_int(ZZ_p_c x)
    ZZ_p_c int_to_ZZ_p(int i)
    void ZZ_p_modulus(ZZ_c* mod, ZZ_p_c* x)
