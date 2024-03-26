# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport mat_GF2_c, vec_GF2_c, GF2_c


cdef extern from "ntlwrap.h":
    void mat_GF2_add "add"( mat_GF2_c x, mat_GF2_c a, mat_GF2_c b)
    void mat_GF2_sub "sub"( mat_GF2_c x, mat_GF2_c a, mat_GF2_c b)
    void mat_GF2_mul "mul"( mat_GF2_c x, mat_GF2_c a, mat_GF2_c b)
    void mat_GF2_negate "NTL::negate"(mat_GF2_c x, mat_GF2_c a)
    void mat_GF2_power "NTL::power"(mat_GF2_c t, mat_GF2_c x, long e)
    GF2_c mat_GF2_determinant "determinant"(mat_GF2_c m)
    void mat_GF2_transpose "transpose"(mat_GF2_c r, mat_GF2_c m)
    long mat_GF2_IsZero "IsZero"(mat_GF2_c x)

    long mat_GF2_gauss "gauss"(mat_GF2_c A, long w)
    void mat_GF2_solve "solve"(GF2_c d, vec_GF2_c X, mat_GF2_c A, vec_GF2_c b)
    void mat_GF2_inv "inv" (mat_GF2_c X, mat_GF2_c A)

    long mat_GF2_IsIdent "IsIdent"(mat_GF2_c A, long n)
    long mat_GF2_IsDiag "IsDiag"(mat_GF2_c A, long n, GF2_c d)

    void mat_GF2_image "image"(mat_GF2_c X, mat_GF2_c A)
    void mat_GF2_kernel "kernel" (mat_GF2_c X, mat_GF2_c A)

    void vec_GF2_conv_mat_GF2 "conv" (vec_GF2_c out, mat_GF2_c inp)
    void mat_GF2_conv_vec_GF2(mat_GF2_c out, vec_GF2_c inp)


cdef extern from "ntlwrap_impl.h":
    void mat_GF2_setitem(mat_GF2_c* x, int i, int j, GF2_c* z)
