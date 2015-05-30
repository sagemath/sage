#ifdef __cplusplus
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/HNF.h>
#include <NTL/LLL.h>
#include <gmp.h>
using namespace NTL;
#endif

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

#include "Python.h"

EXTERN void del_charstar(char*);

EXTERN void setup_NTL_error_callback(void (*function)(const char*, void*), void* context);

////////  ZZ //////////

#ifndef __cplusplus
struct ZZ;
#endif

EXTERN int ZZ_to_int(const struct ZZ* x);
EXTERN struct ZZ* int_to_ZZ(int value);
EXTERN void ZZ_to_mpz(mpz_t output, const struct ZZ* x);
EXTERN void mpz_to_ZZ(struct ZZ *output, const mpz_t x);
EXTERN void ZZ_set_from_int(struct ZZ* x, int value);
/*Random-number generation */
//EXTERN void setSeed(const struct ZZ* n);
//EXTERN struct ZZ* ZZ_randomBnd(const struct ZZ* x);
//EXTERN struct ZZ* ZZ_randomBits(long n);

#ifdef __cplusplus
EXTERN long ZZ_remove(struct ZZ& x, const struct ZZ& a, const struct ZZ& p);
#endif

////////  ZZ_p //////////

#ifndef __cplusplus
struct ZZ_p;
#endif

#ifdef __cplusplus  // sorry, if you want a C version, feel free to add it
EXTERN int ZZ_p_to_int(const ZZ_p& x);
EXTERN ZZ_p int_to_ZZ_p(int value);
#endif
EXTERN void ZZ_p_set_from_int(struct ZZ_p* x, int value);
EXTERN struct ZZ_p* ZZ_p_pow(const struct ZZ_p* x, long e);
EXTERN void ntl_ZZ_set_modulus(struct ZZ* x);
EXTERN struct ZZ_p* ZZ_p_inv(struct ZZ_p* x);
EXTERN struct ZZ_p* ZZ_p_neg(struct ZZ_p* x);
EXTERN struct ZZ_p* ZZ_p_random(void);
EXTERN void ZZ_p_modulus(struct ZZ* mod, const struct ZZ_p* x);


EXTERN struct ZZ_pContext* ZZ_pContext_new(struct ZZ* p);
EXTERN struct ZZ_pContext* ZZ_pContext_construct(void* mem, struct ZZ* p);

//////// ZZX //////////
#ifndef __cplusplus
struct ZZX;
#endif

EXTERN char* ZZX_repr(struct ZZX* x);
EXTERN struct ZZX* ZZX_copy(struct ZZX* x);
EXTERN void ZZX_setitem_from_int(struct ZZX* x, long i, int value);
EXTERN int ZZX_getitem_as_int(struct ZZX* x, long i);
EXTERN void ZZX_getitem_as_mpz(mpz_t output, struct ZZX* x, long i);
EXTERN struct ZZX* ZZX_div(struct ZZX* x, struct ZZX* y, int* divisible);
EXTERN void ZZX_quo_rem(struct ZZX* x, struct ZZX* other, struct ZZX** r, struct ZZX** q);
EXTERN struct ZZX* ZZX_square(struct ZZX* x);
EXTERN int ZZX_equal(struct ZZX* x, struct ZZX* y);
EXTERN int ZZX_is_monic(struct ZZX* x);
EXTERN struct ZZX* ZZX_neg(struct ZZX* x);
EXTERN struct ZZX* ZZX_left_shift(struct ZZX* x, long n);
EXTERN struct ZZX* ZZX_right_shift(struct ZZX* x, long n);
EXTERN char* ZZX_content(struct ZZX* x);
EXTERN struct ZZX* ZZX_primitive_part(struct ZZX* x);
EXTERN void ZZX_pseudo_quo_rem(struct ZZX* x, struct ZZX* y, struct ZZX** r, struct ZZX** q);
EXTERN struct ZZX* ZZX_gcd(struct ZZX* x, struct ZZX* y);
EXTERN void ZZX_xgcd(struct ZZX* x, struct ZZX* y, struct ZZ** r, struct ZZX** s, struct ZZX** t, int proof);
EXTERN long ZZX_degree(struct ZZX* x);
EXTERN struct ZZ* ZZX_leading_coefficient(struct ZZX* x);
EXTERN char* ZZX_constant_term(struct ZZX* x);
EXTERN void ZZX_set_x(struct ZZX* x);
EXTERN int ZZX_is_x(struct ZZX* x);
EXTERN struct ZZX* ZZX_derivative(struct ZZX* x);
EXTERN struct ZZX* ZZX_reverse(struct ZZX* x);
EXTERN struct ZZX* ZZX_reverse_hi(struct ZZX* x, int hi);
EXTERN struct ZZX* ZZX_truncate(struct ZZX* x, long m);
EXTERN struct ZZX* ZZX_multiply_and_truncate(struct ZZX* x, struct ZZX* y, long m);
EXTERN struct ZZX* ZZX_square_and_truncate(struct ZZX* x, long m);
EXTERN struct ZZX* ZZX_invert_and_truncate(struct ZZX* x, long m);
EXTERN struct ZZX* ZZX_multiply_mod(struct ZZX* x, struct ZZX* y,  struct ZZX* modulus);
EXTERN struct ZZ* ZZX_trace_mod(struct ZZX* x, struct ZZX* y);
/* EXTERN struct ZZ* ZZX_polyeval(struct ZZX* f, struct ZZ* a); */
EXTERN char* ZZX_trace_list(struct ZZX* x);
EXTERN struct ZZ* ZZX_resultant(struct ZZX* x, struct ZZX* y, int proof);
EXTERN struct ZZ* ZZX_norm_mod(struct ZZX* x, struct ZZX* y, int proof);
EXTERN struct ZZ* ZZX_discriminant(struct ZZX* x, int proof);
EXTERN struct ZZX* ZZX_charpoly_mod(struct ZZX* x, struct ZZX* y, int proof);
EXTERN struct ZZX* ZZX_minpoly_mod(struct ZZX* x, struct ZZX* y);
EXTERN void ZZX_clear(struct ZZX* x);
EXTERN void ZZX_preallocate_space(struct ZZX* x, long n);

//////// ZZXFactoring //////////

// OUTPUT: v -- pointer to list of n ZZX elements (the squarefree factors)
//         e -- point to list of e longs (the exponents)
//         n -- length of above two lists
//  The lists v and e are mallocd, and must be freed by the calling code.
EXTERN void ZZX_squarefree_decomposition(struct ZZX*** v, long** e, long* n, struct ZZX* x);


//////// ZZ_pX //////////
#ifndef __cplusplus
struct ZZ_pX;
#endif

EXTERN struct ZZ_pX* ZZ_pX_init(void);
//EXTERN char* ZZ_pX_repr(struct ZZ_pX* x);
/* EXTERN struct ZZ_pX* ZZ_pX_copy(struct ZZ_pX* x); */
/* EXTERN void ZZ_pX_setitem_from_int(struct ZZ_pX* x, long i, int value); */
/* EXTERN int ZZ_pX_getitem_as_int(struct ZZ_pX* x, long i); */
/* EXTERN struct ZZ_pX* ZZ_pX_div(struct ZZ_pX* x, struct ZZ_pX* y, int* divisible); */
/* EXTERN struct ZZ_pX* ZZ_pX_mod(struct ZZ_pX* x, struct ZZ_pX* y); */
/* EXTERN void ZZ_pX_quo_rem(struct ZZ_pX* x, struct ZZ_pX* other, struct ZZ_pX** r, struct ZZ_pX** q); */
/* EXTERN struct ZZ_pX* ZZ_pX_square(struct ZZ_pX* x); */
/* EXTERN int ZZ_pX_equal(struct ZZ_pX* x, struct ZZ_pX* y); */
/* EXTERN int ZZ_pX_is_monic(struct ZZ_pX* x); */
/* EXTERN struct ZZ_pX* ZZ_pX_neg(struct ZZ_pX* x); */
/* EXTERN struct ZZ_pX* ZZ_pX_left_shift(struct ZZ_pX* x, long n); */
/* EXTERN struct ZZ_pX* ZZ_pX_right_shift(struct ZZ_pX* x, long n); */
/* EXTERN void ZZ_pX_quo_rem(struct ZZ_pX* x, struct ZZ_pX* y, struct ZZ_pX** r, struct ZZ_pX** q); */
/* EXTERN struct ZZ_pX* ZZ_pX_gcd(struct ZZ_pX* x, struct ZZ_pX* y); */
/* EXTERN void ZZ_pX_xgcd(struct ZZ_pX** d, struct ZZ_pX** s, struct ZZ_pX** t, struct ZZ_pX* a, struct ZZ_pX* b); */
/* EXTERN void ZZ_pX_plain_xgcd(struct ZZ_pX** d, struct ZZ_pX** s, struct ZZ_pX** t, struct ZZ_pX* a, struct ZZ_pX* b); */
/* EXTERN long ZZ_pX_degree(struct ZZ_pX* x); */
/* EXTERN void ZZ_pX_set_x(struct ZZ_pX* x); */
/* EXTERN int ZZ_pX_is_x(struct ZZ_pX* x); */
/* EXTERN struct ZZ_pX* ZZ_pX_derivative(struct ZZ_pX* x); */
/* EXTERN struct ZZ_pX* ZZ_pX_reverse(struct ZZ_pX* x); */
/* EXTERN struct ZZ_pX* ZZ_pX_reverse_hi(struct ZZ_pX* x, int hi); */
/* EXTERN struct ZZ_pX* ZZ_pX_truncate(struct ZZ_pX* x, long m); */
/* EXTERN struct ZZ_pX* ZZ_pX_multiply_and_truncate(struct ZZ_pX* x, struct ZZ_pX* y, long m); */
/* EXTERN struct ZZ_pX* ZZ_pX_square_and_truncate(struct ZZ_pX* x, long m); */
/* EXTERN struct ZZ_pX* ZZ_pX_invert_and_truncate(struct ZZ_pX* x, long m); */
/* EXTERN struct ZZ_pX* ZZ_pX_multiply_mod(struct ZZ_pX* x, struct ZZ_pX* y,  struct ZZ_pX* modulus); */
/* EXTERN struct ZZ_p* ZZ_pX_trace_mod(struct ZZ_pX* x, struct ZZ_pX* y); */
EXTERN char* ZZ_pX_trace_list(struct ZZ_pX* x);
/* EXTERN struct ZZ_p* ZZ_pX_resultant(struct ZZ_pX* x, struct ZZ_pX* y); */
/* EXTERN struct ZZ_p* ZZ_pX_norm_mod(struct ZZ_pX* x, struct ZZ_pX* y); */
/* EXTERN struct ZZ_pX* ZZ_pX_charpoly_mod(struct ZZ_pX* x, struct ZZ_pX* y); */
/* EXTERN struct ZZ_pX* ZZ_pX_minpoly_mod(struct ZZ_pX* x, struct ZZ_pX* y); */
/* EXTERN void ZZ_pX_clear(struct ZZ_pX* x); */
// EXTERN void ZZ_pX_preallocate_space(struct ZZ_pX* x, long n);

// Factoring elements of ZZ_pX:
// OUTPUT: v -- pointer to list of n ZZ_pX elements (the irred factors)
//         e -- point to list of e longs (the exponents)
//         n -- length of above two lists
//  The lists v and e are mallocd, and must be freed by the calling code.
EXTERN void ZZ_pX_factor(struct ZZ_pX*** v, long** e, long* n, struct ZZ_pX* x, long verbose);
EXTERN void ZZ_pX_linear_roots(struct ZZ_p*** v, long* n, struct ZZ_pX* f);

#ifdef __cplusplus
EXTERN void ZZ_pX_conv_modulus(struct ZZ_pX &fout, const struct ZZ_pX &fin, const struct ZZ_pContext &mod);
EXTERN void ZZ_pEX_conv_modulus(struct ZZ_pEX &fout, const struct ZZ_pEX &fin, const struct ZZ_pContext &mod);
EXTERN void ZZ_pX_min_val_coeff(long &valuation, long &index, const struct ZZ_pX &f, const struct ZZ &p);
EXTERN long ZZ_pX_get_val_coeff(const struct ZZ_pX &f, const struct ZZ &p, long i);
EXTERN void ZZ_pX_left_pshift(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ &pn, const struct ZZ_pContext &c);
EXTERN void ZZ_pX_right_pshift(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ &pn, const struct ZZ_pContext &c);
EXTERN void ZZ_pX_InvMod_newton_unram(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ_pXModulus &F, const struct ZZ_pContext &cpn, const struct ZZ_pContext &cp);
EXTERN void ZZ_pX_InvMod_newton_ram(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ_pXModulus &F, const struct ZZ_pContext &cpn);

#endif

//////// zz_p //////////

#ifndef __cplusplus
struct zz_p;
#endif

#define zz_p_set_from_long( obj1, obj2 )\
        (obj1) = (obj2)
#define NTL_zz_p_DOUBLE_EQUALS( obj1, obj2 )\
        ((obj1) == (obj2))

EXTERN struct zz_pContext* zz_pContext_new(long p);
EXTERN struct zz_pContext* zz_pContext_construct(void* mem, long p);
EXTERN void zz_pContext_restore(struct zz_pContext* ctx);

//////// zz_pX //////////

#ifndef __cplusplus
struct zz_pX;
#endif

#define NTL_zz_pX_DOUBLE_EQUALS( obj1, obj2 )\
        ((obj1) == (obj2))

//////// ZZ_pEContext ///////////////

#ifndef __cplusplus
struct ZZ_pEContext;
#endif

EXTERN struct ZZ_pEContext* ZZ_pEContext_new(struct ZZ_pX *f);
EXTERN struct ZZ_pEContext* ZZ_pEContext_construct(void* mem, struct ZZ_pX *f);
EXTERN void ZZ_pEContext_restore(struct ZZ_pEContext* ctx);

//////// ZZ_pE ////////////

#ifndef __cplusplus
struct ZZ_pE;
#endif

EXTERN struct ZZ_pX ZZ_pE_to_ZZ_pX(struct ZZ_pE x);

//////// ZZ_pEX /////////

#ifndef __cplusplus
struct ZZ_pEX;
#endif

//////// mat_ZZ //////////

#ifndef __cplusplus
typedef struct {} mat_ZZ;
#endif

EXTERN void mat_ZZ_SetDims(mat_ZZ* mZZ, long nrows, long ncols);
EXTERN mat_ZZ* mat_ZZ_pow(const mat_ZZ* x, long e);
EXTERN long mat_ZZ_nrows(const mat_ZZ* x);
EXTERN long mat_ZZ_ncols(const mat_ZZ* x);
EXTERN void mat_ZZ_setitem(mat_ZZ* x, int i, int j, const struct ZZ* z);
EXTERN struct ZZ* mat_ZZ_getitem(const mat_ZZ* x, int i, int j);
EXTERN struct ZZ* mat_ZZ_determinant(const mat_ZZ* x, long deterministic);
EXTERN mat_ZZ* mat_ZZ_HNF(const mat_ZZ* A, const struct ZZ* D);
EXTERN struct ZZX* mat_ZZ_charpoly(const mat_ZZ* A);
EXTERN long mat_ZZ_LLL(struct ZZ **det, mat_ZZ *x, long a, long b, long verbose);
EXTERN long mat_ZZ_LLL_U(struct ZZ **det, mat_ZZ *x, mat_ZZ *U, long a, long b, long verbose);

/* //////// ZZ_p ////////// */
/* #ifndef __cplusplus */
/* struct ZZ_p; */
/* #endif */

/* EXTERN void ZZ_p_set_modulus(const struct ZZ* p); */
/* EXTERN struct ZZ_p* new_ZZ_p(void); */
/* EXTERN void del_ZZ_p(struct ZZ_p* x); */
/* EXTERN struct ZZ_p* ZZ_p_add(const struct ZZ_p* x, const struct ZZ_p* y); */
/* EXTERN struct ZZ_p* ZZ_p_sub(const struct ZZ_p* x, const struct ZZ_p* y); */
/* EXTERN struct ZZ_p* ZZ_p_mul(const struct ZZ_p* x, const struct ZZ_p* y); */
/* EXTERN struct ZZ_p* ZZ_p_pow(const struct ZZ_p* x, long e); */
/* EXTERN int ZZ_p_is_zero(struct ZZ_p*x ); */
/* EXTERN int ZZ_p_is_one(struct ZZ_p*x ); */


//////// ZZ_pE //////////
#ifndef __cplusplus
struct ZZ_pE;
#endif

// EXTERN struct ZZ_pE* new_ZZ_pE



//////// ZZ_pEX //////////

//#ifndef __cplusplus
//struct ZZ_pEX;
//#endif

//EXTERN struct ZZ_pEX* new_ZZ_pEX

/////// GF2X ////////////////
#ifndef __cplusplus
struct GF2X;
#endif

/////// GF2EContext ////////////////

#ifndef __cplusplus
struct GF2EContext;
struct GF2X_c;
#endif

EXTERN struct GF2EContext* GF2EContext_new(struct GF2X_c* p);
EXTERN struct GF2EContext* GF2EContext_construct(void *mem, const struct GF2X *p);

//////// mat_GF2E //////////

#ifndef __cplusplus
typedef struct {} mat_GF2E;
struct GF2E;
#endif

EXTERN void mat_GF2E_setitem(mat_GF2E* x, int i, int j, const struct GF2E* z);

//////// mat_GF2 //////////

#ifndef __cplusplus
typedef struct {} mat_GF2;
struct GF2;
#endif

EXTERN void mat_GF2_setitem(mat_GF2* x, int i, int j, const struct GF2* z);
