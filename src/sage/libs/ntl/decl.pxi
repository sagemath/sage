include "sage/ext/cdefs.pxi"
include "sage/ext/python.pxi"

cdef extern from "ccobject.h":
    pass

cdef extern from "ntl_wrap.h":
    long NTL_OVFBND
    bint NTL_OVERFLOW(long, long, long)

from sage.libs.ntl.ntl_ZZ_decl cimport *
from sage.libs.ntl.ntl_lzz_pX_decl cimport *
from sage.libs.ntl.ntl_ZZ_pX_decl cimport *
from sage.libs.ntl.ntl_ZZ_pContext_decl cimport *
from sage.libs.ntl.ntl_ZZ_p_decl cimport *
from sage.libs.ntl.ntl_vec_ZZ_p_decl cimport *
from sage.libs.ntl.ntl_ZZX_decl cimport *
from sage.libs.ntl.ntl_lzz_pContext_decl cimport *
from sage.libs.ntl.ntl_ZZ_pEContext_decl cimport *
from sage.libs.ntl.ntl_ZZ_pE_decl cimport *
from sage.libs.ntl.ntl_vec_ZZ_pE_decl cimport *
from sage.libs.ntl.ntl_ZZ_pEX_decl cimport *

cdef extern from "ntl_wrap.h":
    #### mat_ZZ_c
    cdef cppclass mat_ZZ_c "mat_ZZ":
        pass

    # Some boiler-plate
    mat_ZZ_c* mat_ZZ_construct "Construct<mat_ZZ>"(void *mem)
    void mat_ZZ_destruct "Destruct<mat_ZZ>"(mat_ZZ_c *mem)
    void mat_ZZ_delete "Delete<mat_ZZ>"(mat_ZZ_c *mem)
    object mat_ZZ_to_PyString "_to_PyString<mat_ZZ>"(mat_ZZ_c *x)
    int mat_ZZ_equal "_equal<mat_ZZ>"(mat_ZZ_c x, mat_ZZ_c y)

    void mat_ZZ_mul "mul"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_add "add"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_sub "sub"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_power "power"( mat_ZZ_c x, mat_ZZ_c a, long e)
    void mat_ZZ_CharPoly "CharPoly"(ZZX_c r, mat_ZZ_c m)
    void mat_ZZ_SetDims(mat_ZZ_c* mZZ, long nrows, long ncols)

    mat_ZZ_c* mat_ZZ_pow(mat_ZZ_c* x, long e)
    long mat_ZZ_nrows(mat_ZZ_c* x)
    long mat_ZZ_ncols(mat_ZZ_c* x)
    void mat_ZZ_setitem(mat_ZZ_c* x, int i, int j, ZZ_c* z)
    ZZ_c* mat_ZZ_getitem(mat_ZZ_c* x, int i, int j)
    ZZ_c* mat_ZZ_determinant(mat_ZZ_c* x, long deterministic)
    mat_ZZ_c* mat_ZZ_HNF(mat_ZZ_c* A, ZZ_c* D)
    ZZX_c* mat_ZZ_charpoly(mat_ZZ_c* A)

    cdef long mat_ZZ_LLL(ZZ_c **det, mat_ZZ_c *x, long a, long b, long verbose)
    cdef long mat_ZZ_LLL_U(ZZ_c **det, mat_ZZ_c *x, mat_ZZ_c *U, long a, long b, long verbose)

    cdef long mat_ZZ_LLL_FP   "LLL_FP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_FP_U "LLL_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_QP   "LLL_QP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_QP_U "LLL_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_XD   "LLL_XD"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_XD_U "LLL_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_RR   "LLL_RR"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_RR_U "LLL_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)

    cdef long mat_ZZ_G_LLL_FP   "G_LLL_FP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_FP_U "G_LLL_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_QP   "G_LLL_QP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_QP_U "G_LLL_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_XD   "G_LLL_XD"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_XD_U "G_LLL_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_RR   "G_LLL_RR"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_RR_U "G_LLL_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)

    cdef long mat_ZZ_BKZ_FP     "BKZ_FP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_FP_U   "BKZ_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_XD     "BKZ_XD"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_XD_U   "BKZ_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP     "BKZ_QP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP_U   "BKZ_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP1    "BKZ_QP1"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP1_U  "BKZ_QP1"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_RR     "BKZ_RR"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_RR_U   "BKZ_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)

    cdef long mat_ZZ_G_BKZ_FP     "G_BKZ_FP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_FP_U   "G_BKZ_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_XD     "G_BKZ_XD"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_XD_U   "G_BKZ_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP     "G_BKZ_QP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP_U   "G_BKZ_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP1    "G_BKZ_QP1"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP1_U  "G_BKZ_QP1"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_RR     "G_BKZ_RR"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_RR_U   "G_BKZ_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)

    #### GF2_c
    ctypedef struct GF2_c "struct GF2":
        pass

    GF2_c* GF2_new "New<GF2>"()
    GF2_c* GF2_construct "Construct<GF2>"(void *mem)
    void GF2_destruct "Destruct<GF2>"(GF2_c *mem)
    void GF2_delete "Delete<GF2>"(GF2_c *mem)
    void GF2_from_str "_from_str<GF2>"(GF2_c* dest, char* s)
    object GF2_to_PyString "_to_PyString<GF2>"(GF2_c *x)
    int GF2_equal "_equal<GF2>"(GF2_c x, GF2_c y)
    int GF2_IsOne "IsOne"(GF2_c x)
    int GF2_IsZero "IsZero"(GF2_c x)

    void GF2_add "add"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_sub "sub"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_mul "mul"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_div "div"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_negate "negate"(GF2_c x, GF2_c a)
    void GF2_power "power"(GF2_c t, GF2_c x, long e)
    long GF2_deg "deg"(GF2_c x)

    void GF2_conv_long "conv" (GF2_c x, long i)
    long GF2_conv_to_long "rep" (GF2_c x)

    #### GF2X_c

    ctypedef struct GF2X_c "struct GF2X":
        pass

    long *GF2XHexOutput_c "(&GF2X::HexOutput)" # work-around for Cython bug

    GF2X_c* GF2X_new "New<GF2X>"()
    GF2X_c* GF2X_construct "Construct<GF2X>"(void *mem)
    void GF2X_destruct "Destruct<GF2X>"(GF2X_c *mem)
    void GF2X_delete "Delete<GF2X>"(GF2X_c *mem)
    void GF2X_from_str "_from_str<GF2X>"(GF2X_c* dest, char* s)
    object GF2X_to_PyString "_to_PyString<GF2X>"(GF2X_c *x)
    int GF2X_equal "_equal<GF2X>"(GF2X_c x, GF2X_c y)
    int GF2X_IsOne "IsOne"(GF2X_c x)
    int GF2X_IsZero "IsZero"(GF2X_c x)

    void GF2X_add "add"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_sub "sub"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_mul "mul"( GF2X_c x, GF2X_c a, GF2X_c b)
    void GF2X_negate "negate"(GF2X_c x, GF2X_c a)
    void GF2X_power "power"(GF2X_c t, GF2X_c x, long e)
    long GF2X_deg "deg"(GF2X_c x)

    void GF2X_conv_long "conv" (GF2X_c x, long a)
    void GF2X_conv_GF2 "conv" (GF2X_c x, GF2_c a)

    void GF2X_LeftShift "LeftShift"( GF2X_c r, GF2X_c a, long offset)
    void GF2X_RightShift "RightShift"( GF2X_c r, GF2X_c a, long offset)

    void GF2X_DivRem "DivRem"(GF2X_c q, GF2X_c r, GF2X_c a, GF2X_c b)
    void GF2X_div "div" (GF2X_c q, GF2X_c a, GF2X_c b)
    void GF2X_rem "rem" (GF2X_c r, GF2X_c a, GF2X_c b)
    long GF2X_divide "divide"(GF2X_c q, GF2X_c a, GF2X_c b)

    GF2X_c GF2X_GCD "GCD" (GF2X_c a, GF2X_c b)
    void GF2X_XGCD "XGCD" (GF2X_c r, GF2X_c s, GF2X_c t, GF2X_c a, GF2X_c b)

    void GF2XFromBytes(GF2X_c a, unsigned char *p, long n)
    void BytesFromGF2X "BytesFromGF2X" (unsigned char *p, GF2X_c a, long n)

    GF2_c GF2X_coeff "coeff"(GF2X_c a, long i)
    GF2_c GF2X_LeadCoeff "LeadCoeff"(GF2X_c a)
    GF2_c GF2X_ConstTerm "ConstTerm"(GF2X_c a)
    void GF2X_SetCoeff "SetCoeff"(GF2X_c x, long i, GF2_c a)

    GF2X_c GF2X_diff "diff"(GF2X_c a)
    GF2X_c GF2X_reverse "reverse"(GF2X_c a, long hi)

    long GF2X_weight "weight"(GF2X_c a)
    long GF2X_NumBits "NumBits" (GF2X_c a)
    long GF2X_NumBytes "NumBytes"(GF2X_c a)

    #### GF2XFactoring
    void GF2X_BuildSparseIrred "BuildSparseIrred" (GF2X_c f, long n)
    void GF2X_BuildRandomIrred "BuildRandomIrred" (GF2X_c f, GF2X_c g)
    void GF2X_BuildIrred "BuildIrred" (GF2X_c f, long n)

    #### GF2XModulus_c
    ctypedef struct GF2XModulus_c "struct GF2XModulus":
        pass

    GF2X_c GF2XModulus_GF2X "GF2X" (GF2XModulus_c m)

    GF2X_c GF2X_IrredPolyMod "IrredPolyMod" (GF2X_c g, GF2XModulus_c F)

    #### GF2EContext_c

    ctypedef struct GF2EContext_c "struct GF2EContext":
        void (*restore)()

    GF2EContext_c* GF2EContext_new "New<GF2EContext>"()
    GF2EContext_c* GF2EContext_construct "Construct<GF2EContext>"(void *mem)
    GF2EContext_c* GF2EContext_new_GF2X "GF2EContext_new"(GF2X_c* p)
    GF2EContext_c* GF2EContext_construct_GF2X "GF2EContext_construct"(void *mem, GF2X_c* p)
    void GF2EContext_destruct "Destruct<GF2EContext>"(GF2EContext_c *mem)
    void GF2EContext_delete "Delete<GF2EContext>"(GF2EContext_c *mem)

    #### GF2E_c
    ctypedef struct GF2E_c "struct GF2E":
        pass

    void GF2E_init "GF2E::init"(GF2X_c x)
    long GF2E_degree "GF2E::degree"()
    GF2XModulus_c GF2E_modulus "GF2E::modulus"()

    GF2E_c* GF2E_new "New<GF2E>"()
    GF2E_c* GF2E_construct "Construct<GF2E>"(void *mem)
    void GF2E_destruct "Destruct<GF2E>"(GF2E_c *mem)
    void GF2E_delete "Delete<GF2E>"(GF2E_c *mem)
    void GF2E_from_str "_from_str<GF2E>"(GF2E_c* dest, char* s)
    object GF2E_to_PyString "_to_PyString<GF2E>"(GF2E_c *x)
    int GF2E_equal "_equal<GF2E>"(GF2E_c x, GF2E_c y)
    int GF2E_IsOne "IsOne"(GF2E_c x)
    int GF2E_IsZero "IsZero"(GF2E_c x)

    void GF2E_add "add"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_sub "sub"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_mul "mul"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_div "div"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_power "power"(GF2E_c t, GF2E_c x, long e)
    long GF2E_deg "deg"(GF2E_c x)

    void GF2E_conv_GF2X "conv" (GF2E_c out, GF2X_c inp)
    void GF2E_conv_long "conv" (GF2E_c out, long inp)
    void GF2E_conv_ZZ "conv" (GF2E_c out, ZZ_c inp)
    void GF2E_conv_GF2 "conv" (GF2E_c out, GF2_c inp)
    GF2X_c GF2E_rep "rep"(GF2E_c x)

    GF2E_c GF2E_random "random_GF2E"()

    GF2_c GF2E_trace "trace"(GF2E_c x)

    #### GF2EX_c
    ctypedef struct GF2EX_c "struct GF2EX":
        pass

    GF2EX_c* GF2EX_new "New<GF2EX>"()
    GF2EX_c* GF2EX_construct "Construct<GF2EX>"(void *mem)
    void GF2EX_destruct "Destruct<GF2EX>"(GF2EX_c *mem)
    void GF2EX_delete "Delete<GF2EX>"(GF2EX_c *mem)
    void GF2EX_from_str "_from_str<GF2EX>"(GF2EX_c* dest, char* s)
    object GF2EX_to_PyString "_to_PyString<GF2EX>"(GF2EX_c *x)
    int GF2EX_equal "_equal<GF2EX>"(GF2EX_c x, GF2EX_c y)
    void GF2EX_add "add"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_sub "sub"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_mul "mul"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_negate "negate"(GF2EX_c x, GF2EX_c a)
    void GF2EX_power "power"(GF2EX_c t, GF2EX_c x, long e)
    int GF2EX_IsOne "IsOne"(GF2EX_c x)
    int GF2EX_IsZero "IsZero"(GF2EX_c x)

    #### vec_GF2E_c
    ctypedef struct vec_GF2E_c "struct vec_GF2E":
        pass

    vec_GF2E_c* vec_GF2E_new "New<vec_GF2E>"()
    vec_GF2E_c* vec_GF2E_construct "Construct<vec_GF2E>"(void *mem)
    void vec_GF2E_destruct "Destruct<vec_GF2E>"(vec_GF2E_c *mem)
    void vec_GF2E_delete "Delete<vec_GF2E>"(vec_GF2E_c *mem)
    void vec_GF2E_from_str "_from_str<vec_GF2E>"(vec_GF2E_c* dest, char* s)
    object vec_GF2E_to_PyString "_to_PyString<vec_GF2E>"(vec_GF2E_c *x)

    #### mat_GF2E_c
    cdef cppclass mat_GF2E_c "mat_GF2E":
        void (*SetDims)(long nrows, long ncols)
        long (*NumRows)()
        long (*NumCols)()
        GF2E_c (*get "operator()") (long i, long j)

    mat_GF2E_c* mat_GF2E_new "New<mat_GF2E>"()
    mat_GF2E_c* mat_GF2E_construct "Construct<mat_GF2E>"(void *mem)
    void mat_GF2E_destruct "Destruct<mat_GF2E>"(mat_GF2E_c *mem)
    void mat_GF2E_delete "Delete<mat_GF2E>"(mat_GF2E_c *mem)
    void mat_GF2E_from_str "_from_str<mat_GF2E>"(mat_GF2E_c* dest, char* s)
    object mat_GF2E_to_PyString "_to_PyString<mat_GF2E>"(mat_GF2E_c *x)
    int mat_GF2E_equal "_equal<mat_GF2E>"(mat_GF2E_c x, mat_GF2E_c y)
    void mat_GF2E_add "add"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_sub "sub"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_mul "mul"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_negate "negate"(mat_GF2E_c x, mat_GF2E_c a)
    void mat_GF2E_power "power"(mat_GF2E_c t, mat_GF2E_c x, long e)
    GF2E_c mat_GF2E_determinant "determinant"(mat_GF2E_c m)
    void mat_GF2E_transpose "transpose"(mat_GF2E_c r, mat_GF2E_c m)
    long mat_GF2E_IsZero "IsZero"(mat_GF2E_c x)
    void mat_GF2E_setitem(mat_GF2E_c* x, int i, int j, GF2E_c* z)

    long mat_GF2E_gauss "gauss"(mat_GF2E_c A, long w)
    void mat_GF2E_solve "solve"(GF2E_c d, vec_GF2E_c X, mat_GF2E_c A, vec_GF2E_c b)
    void mat_GF2E_inv "inv" (mat_GF2E_c X, mat_GF2E_c A)


    long mat_GF2E_IsIdent "IsIdent"(mat_GF2E_c A, long n)
    long mat_GF2E_IsDiag "IsDiag"(mat_GF2E_c A, long n, GF2E_c d)


    void mat_GF2E_image "image"(mat_GF2E_c X, mat_GF2E_c A)
    void mat_GF2E_kernel "kernel" (mat_GF2E_c X, mat_GF2E_c A)

    void vec_GF2E_conv_mat_GF2E "conv" (vec_GF2E_c out, mat_GF2E_c inp)
    void mat_GF2E_conv_vec_GF2E(mat_GF2E_c out, vec_GF2E_c inp)

    #### vec_GF2_c
    ctypedef struct vec_GF2_c "struct vec_GF2":
        pass

    vec_GF2_c* vec_GF2_new "New<vec_GF2>"()
    vec_GF2_c* vec_GF2_construct "Construct<vec_GF2>"(void *mem)
    void vec_GF2_destruct "Destruct<vec_GF2>"(vec_GF2_c *mem)
    void vec_GF2_delete "Delete<vec_GF2>"(vec_GF2_c *mem)
    void vec_GF2_from_str "_from_str<vec_GF2>"(vec_GF2_c* dest, char* s)
    object vec_GF2_to_PyString "_to_PyString<vec_GF2>"(vec_GF2_c *x)


    #### mat_GF2_c
    cdef cppclass mat_GF2_c "mat_GF2":
        void (*SetDims)(long nrows, long ncols)
        long (*NumRows)()
        long (*NumCols)()
        GF2_c (*get "operator()") (long i, long j)

    mat_GF2_c* mat_GF2_new "New<mat_GF2>"()
    mat_GF2_c* mat_GF2_construct "Construct<mat_GF2>"(void *mem)
    void mat_GF2_destruct "Destruct<mat_GF2>"(mat_GF2_c *mem)
    void mat_GF2_delete "Delete<mat_GF2>"(mat_GF2_c *mem)
    void mat_GF2_from_str "_from_str<mat_GF2>"(mat_GF2_c* dest, char* s)
    object mat_GF2_to_PyString "_to_PyString<mat_GF2>"(mat_GF2_c *x)
    int mat_GF2_equal "_equal<mat_GF2>"(mat_GF2_c x, mat_GF2_c y)
    void mat_GF2_add "add"( mat_GF2_c x, mat_GF2_c a, mat_GF2_c b)
    void mat_GF2_sub "sub"( mat_GF2_c x, mat_GF2_c a, mat_GF2_c b)
    void mat_GF2_mul "mul"( mat_GF2_c x, mat_GF2_c a, mat_GF2_c b)
    void mat_GF2_negate "negate"(mat_GF2_c x, mat_GF2_c a)
    void mat_GF2_power "power"(mat_GF2_c t, mat_GF2_c x, long e)
    GF2_c mat_GF2_determinant "determinant"(mat_GF2_c m)
    void mat_GF2_transpose "transpose"(mat_GF2_c r, mat_GF2_c m)
    long mat_GF2_IsZero "IsZero"(mat_GF2_c x)
    void mat_GF2_setitem(mat_GF2_c* x, int i, int j, GF2_c* z)

    long mat_GF2_gauss "gauss"(mat_GF2_c A, long w)
    void mat_GF2_solve "solve"(GF2_c d, vec_GF2_c X, mat_GF2_c A, vec_GF2_c b)
    void mat_GF2_inv "inv" (mat_GF2_c X, mat_GF2_c A)


    long mat_GF2_IsIdent "IsIdent"(mat_GF2_c A, long n)
    long mat_GF2_IsDiag "IsDiag"(mat_GF2_c A, long n, GF2_c d)


    void mat_GF2_image "image"(mat_GF2_c X, mat_GF2_c A)
    void mat_GF2_kernel "kernel" (mat_GF2_c X, mat_GF2_c A)

    void vec_GF2_conv_mat_GF2 "conv" (vec_GF2_c out, mat_GF2_c inp)
    void mat_GF2_conv_vec_GF2(mat_GF2_c out, vec_GF2_c inp)
