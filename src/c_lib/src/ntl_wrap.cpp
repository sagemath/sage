#include <iostream>
#include <sstream>
using namespace std;

#include "ntl_wrap.h"
#include <NTL/mat_poly_ZZ.h>
#include <NTL/LLL.h>
#include <NTL/tools.h>

void del_charstar(char* a) {
    delete[] a;
}


void setup_NTL_error_callback(void (*function)(const char*, void*), void* context)
{
    NTL::SetErrorCallbackFunction(function, context);
}


//////// ZZ //////////

/* Return value is only valid if the result should fit into an int.
   AUTHOR: David Harvey (2008-06-08) */
int ZZ_to_int(const ZZ* x)
{
    return to_int(*x);
}

/* Returns a *new* ZZ object.
   AUTHOR: David Harvey (2008-06-08) */
struct ZZ* int_to_ZZ(int value)
{
    ZZ* output = new ZZ();
    conv(*output, value);
    return output;
}

/* Copies the ZZ into the mpz_t
   Assumes output has been mpz_init'd.
   AUTHOR: David Harvey
           Joel B. Mohler moved the ZZX_getitem_as_mpz code out to this function (2007-03-13) */
void ZZ_to_mpz(mpz_t output, const struct ZZ* x)
{
    unsigned char stack_bytes[4096];
    int use_heap;
    unsigned long size = NumBytes(*x);
    use_heap = (size > sizeof(stack_bytes));
    unsigned char* bytes = use_heap ? (unsigned char*) malloc(size) : stack_bytes;
    BytesFromZZ(bytes, *x, size);
    mpz_import(output, size, -1, 1, 0, 0, bytes);
    if (sign(*x) < 0)
        mpz_neg(output, output);
    if (use_heap)
        free(bytes);
}

// Ok, I know that this is obvious
// I just wanted to document the appearance of the magic number 8 in the code below
#define bits_in_byte        8

/* Copies the mpz_t into the ZZ
   AUTHOR: Joel B. Mohler (2007-03-15) */
void mpz_to_ZZ(struct ZZ* output, const mpz_t x)
{
    unsigned char stack_bytes[4096];
    int use_heap;
    size_t size = (mpz_sizeinbase(x, 2) + bits_in_byte-1) / bits_in_byte;
    use_heap = (size > sizeof(stack_bytes));
    void* bytes = use_heap ? malloc(size) : stack_bytes;
    size_t words_written;
    mpz_export(bytes, &words_written, -1, 1, 0, 0, x);
    clear(*output);
    ZZFromBytes(*output, (unsigned char *)bytes, words_written);
    if (mpz_sgn(x) < 0)
        NTL::negate(*output, *output);
    if (use_heap)
        free(bytes);
}

/* Sets given ZZ to value
   AUTHOR: David Harvey (2008-06-08) */
void ZZ_set_from_int(ZZ* x, int value)
{
    conv(*x, value);
}

long ZZ_remove(struct ZZ &dest, const struct ZZ &src, const struct ZZ &f)
{
    // Based on the code for mpz_remove
    ZZ fpow[40];            // inexaustible...until year 2020 or so
    ZZ x, rem;
    long pwr;
    int p;

    if (compare(f, 1) <= 0 && compare(f, -1) >= 0)
        Error("Division by zero");

    if (compare(src, 0) == 0)
    {
        if (src != dest)
           dest = src;
        return 0;
    }

    if (compare(f, 2) == 0)
    {
        dest = src;
        return MakeOdd(dest);
    }

    /* We could perhaps compute mpz_scan1(src,0)/mpz_scan1(f,0).  It is an
     upper bound of the result we're seeking.  We could also shift down the
     operands so that they become odd, to make intermediate values smaller.  */

    pwr = 0;
    fpow[0] = ZZ(f);
    dest = src;
    rem = ZZ();
    x = ZZ();

    /* Divide by f, f^2, ..., f^(2^k) until we get a remainder for f^(2^k).  */
    for (p = 0;;p++)
    {
        DivRem(x, rem, dest, fpow[p]);
        if (compare(rem, 0) != 0)
            break;
        fpow[p+1] = ZZ();
        mul(fpow[p+1], fpow[p], fpow[p]);
        dest = x;
    }

    pwr = (1 << p) - 1;

    /* Divide by f^(2^(k-1)), f^(2^(k-2)), ..., f for all divisors that give a
       zero remainder.  */
    while (--p >= 0)
    {
        DivRem(x, rem, dest, fpow[p]);
        if (compare(rem, 0) == 0)
        {
            pwr += 1 << p;
            dest = x;
        }
    }
    return pwr;
}

//////// ZZ_p //////////

/* Return value is only valid if the result should fit into an int.
   AUTHOR: David Harvey (2008-06-08) */
int ZZ_p_to_int(const ZZ_p& x )
{
    return ZZ_to_int(&rep(x));
}

/* Returns a *new* ZZ_p object.
   AUTHOR: David Harvey (2008-06-08) */
ZZ_p int_to_ZZ_p(int value)
{
    ZZ_p r;
    r = value;
    return r;
}

/* Sets given ZZ_p to value
   AUTHOR: David Harvey (2008-06-08) */
void ZZ_p_set_from_int(ZZ_p* x, int value)
{
    conv(*x, value);
}

void ZZ_p_modulus(struct ZZ* mod, const struct ZZ_p* x)
{
    (*mod) = x->modulus();
}

struct ZZ_p* ZZ_p_pow(const struct ZZ_p* x, long e)
{
    ZZ_p *z = new ZZ_p();
    power(*z, *x, e);
    return z;
}

void ntl_ZZ_set_modulus(ZZ* x)
{
    ZZ_p::init(*x);
}

ZZ_p* ZZ_p_inv(ZZ_p* x)
{
    ZZ_p *z = new ZZ_p();
    inv(*z, *x);
    return z;
}

ZZ_p* ZZ_p_random(void)
{
    ZZ_p *z = new ZZ_p();
    random(*z);
    return z;
}

struct ZZ_p* ZZ_p_neg(struct ZZ_p* x)
{
    return new ZZ_p(-(*x));
}



///////////////////////////////////////////////
//////// ZZX //////////
///////////////////////////////////////////////

char* ZZX_repr(struct ZZX* x)
{
    ostringstream instore;
    instore << (*x);
    int n = strlen(instore.str().data());
    char* buf = new char[n+1];
    strcpy(buf, instore.str().data());
    return buf;
}

struct ZZX* ZZX_copy(struct ZZX* x) {
    return new ZZX(*x);
}

/* Sets ith coefficient of x to value.
   AUTHOR: David Harvey (2006-06-08) */
void ZZX_setitem_from_int(struct ZZX* x, long i, int value)
{
    SetCoeff(*x, i, value);
}

/* Returns ith coefficient of x.
   Return value is only valid if the result should fit into an int.
   AUTHOR: David Harvey (2006-06-08) */
int ZZX_getitem_as_int(struct ZZX* x, long i)
{
    return ZZ_to_int(&coeff(*x, i));
}

/* Copies ith coefficient of x to output.
   Assumes output has been mpz_init'd.
   AUTHOR: David Harvey (2007-02) */
void ZZX_getitem_as_mpz(mpz_t output, struct ZZX* x, long i)
{
    const ZZ& z = coeff(*x, i);
    ZZ_to_mpz(output, &z);
}

struct ZZX* ZZX_div(struct ZZX* x, struct ZZX* y, int* divisible)
{
    struct ZZX* z = new ZZX();
    *divisible = divide(*z, *x, *y);
    return z;
}



void ZZX_quo_rem(struct ZZX* x, struct ZZX* other, struct ZZX** r, struct ZZX** q)
{
    struct ZZX *qq = new ZZX(), *rr = new ZZX();
    DivRem(*qq, *rr, *x, *other);
    *r = rr; *q = qq;
}


struct ZZX* ZZX_square(struct ZZX* x)
{
    struct ZZX* s = new ZZX();
    sqr(*s, *x);
    return s;
}


int ZZX_is_monic(struct ZZX* x)
{
    return IsOne(LeadCoeff(*x));
}


struct ZZX* ZZX_neg(struct ZZX* x)
{
    struct ZZX* y = new ZZX();
    *y = -*x;
    return y;
}


struct ZZX* ZZX_left_shift(struct ZZX* x, long n)
{
    struct ZZX* y = new ZZX();
    LeftShift(*y, *x, n);
    return y;
}


struct ZZX* ZZX_right_shift(struct ZZX* x, long n)
{
    struct ZZX* y = new ZZX();
    RightShift(*y, *x, n);
    return y;
}

struct ZZX* ZZX_primitive_part(struct ZZX* x)
{
    struct ZZX* p = new ZZX();
    PrimitivePart(*p, *x);
    return p;
}


void ZZX_pseudo_quo_rem(struct ZZX* x, struct ZZX* y, struct ZZX** r, struct ZZX** q)
{
    *r = new ZZX();
    *q = new ZZX();
    PseudoDivRem(**q, **r, *x, *y);
}


struct ZZX* ZZX_gcd(struct ZZX* x, struct ZZX* y)
{
    struct ZZX* g = new ZZX();
    GCD(*g, *x, *y);
    return g;
}


void ZZX_xgcd(struct ZZX* x, struct ZZX* y, struct ZZ** r, struct ZZX** s,
          struct ZZX** t, int proof)
{
    *r = new ZZ();
    *s = new ZZX();
    *t = new ZZX();
    XGCD(**r, **s, **t, *x, *y, proof);
}


long ZZX_degree(struct ZZX* x)
{
    return deg(*x);
}

void ZZX_set_x(struct ZZX* x)
{
    SetX(*x);
}


int ZZX_is_x(struct ZZX* x)
{
    return IsX(*x);
}


struct ZZX* ZZX_derivative(struct ZZX* x)
{
    ZZX* d = new ZZX();
    diff(*d, *x);
    return d;
}


struct ZZX* ZZX_reverse(struct ZZX* x)
{
    ZZX* r = new ZZX();
    reverse(*r, *x);
    return r;
}

struct ZZX* ZZX_reverse_hi(struct ZZX* x, int hi)
{
    ZZX* r = new ZZX();
    reverse(*r, *x, hi);
    return r;
}


struct ZZX* ZZX_truncate(struct ZZX* x, long m)
{
    ZZX* t = new ZZX();
    trunc(*t, *x, m);
    return t;
}


struct ZZX* ZZX_multiply_and_truncate(struct ZZX* x, struct ZZX* y, long m)
{
    ZZX* t = new ZZX();
    MulTrunc(*t, *x, *y, m);
    return t;
}


struct ZZX* ZZX_square_and_truncate(struct ZZX* x, long m)
{
    ZZX* t = new ZZX();
    SqrTrunc(*t, *x, m);
    return t;
}


struct ZZX* ZZX_invert_and_truncate(struct ZZX* x, long m)
{
    ZZX* t = new ZZX();
    InvTrunc(*t, *x, m);
    return t;
}


struct ZZX* ZZX_multiply_mod(struct ZZX* x, struct ZZX* y,  struct ZZX* modulus)
{
    ZZX* p = new ZZX();
    MulMod(*p, *x, *y, *modulus);
    return p;
}


struct ZZ* ZZX_trace_mod(struct ZZX* x, struct ZZX* y)
{
    ZZ* p = new ZZ();
    TraceMod(*p, *x, *y);
    return p;
}


char* ZZX_trace_list(struct ZZX* x)
{
    vec_ZZ v;
    TraceVec(v, *x);
    ostringstream instore;
    instore << v;
    int n = strlen(instore.str().data());
    char* buf = new char[n+1];
    strcpy(buf, instore.str().data());
    return buf;
}


struct ZZ* ZZX_resultant(struct ZZX* x, struct ZZX* y, int proof)
{
    ZZ* res = new ZZ();
    resultant(*res, *x, *y, proof);
    return res;
}


struct ZZ* ZZX_norm_mod(struct ZZX* x, struct ZZX* y, int proof)
{
    ZZ* res = new ZZ();
    NormMod(*res, *x, *y, proof);
    return res;
}


struct ZZ* ZZX_discriminant(struct ZZX* x, int proof)
{
    ZZ* d = new ZZ();
    discriminant(*d, *x, proof);
    return d;
}


struct ZZX* ZZX_charpoly_mod(struct ZZX* x, struct ZZX* y, int proof)
{
    ZZX* f = new ZZX();
    CharPolyMod(*f, *x, *y, proof);
    return f;
}


struct ZZX* ZZX_minpoly_mod(struct ZZX* x, struct ZZX* y)
{
    ZZX* f = new ZZX();
    MinPolyMod(*f, *x, *y);
    return f;
}


void ZZX_clear(struct ZZX* x)
{
    clear(*x);
}


void ZZX_preallocate_space(struct ZZX* x, long n)
{
    x->SetMaxLength(n);
}

/*
EXTERN struct ZZ* ZZX_polyeval(struct ZZX* f, struct ZZ* a)
{
    ZZ* b = new ZZ();
    *b = PolyEval(*f, *a);
    return b;
}
*/

void ZZX_squarefree_decomposition(struct ZZX*** v, long** e, long* n, struct ZZX* x)
{
    vec_pair_ZZX_long factors;
    SquareFreeDecomp(factors, *x);
    *n = factors.length();
    *v = (ZZX**) malloc(sizeof(ZZX*) * (*n));
    *e = (long*) malloc(sizeof(long) * (*n));
    for (long i = 0; i < (*n); i++) {
        (*v)[i] = new ZZX(factors[i].a);
        (*e)[i] = factors[i].b;
    }
}

///////////////////////////////////////////////
//////// ZZ_pX //////////
///////////////////////////////////////////////

// char* ZZ_pX_repr(struct ZZ_pX* x)
// {
//   ostringstream instore;
//   instore << (*x);
//   int n = strlen(instore.str().data());
//   char* buf = new char[n+1];
//   strcpy(buf, instore.str().data());
//   return buf;
// }

// void ZZ_pX_dealloc(struct ZZ_pX* x) {
//   delete x;
// }

// struct ZZ_pX* ZZ_pX_copy(struct ZZ_pX* x) {
//   return new ZZ_pX(*x);
// }

// /* Sets ith coefficient of x to value.
//    AUTHOR: David Harvey (2008-06-08) */
// void ZZ_pX_setitem_from_int(struct ZZ_pX* x, long i, int value)
// {
//   SetCoeff(*x, i, value);
// }

// /* Returns ith coefficient of x.
//    Return value is only valid if the result should fit into an int.
//    AUTHOR: David Harvey (2008-06-08) */
// int ZZ_pX_getitem_as_int(struct ZZ_pX* x, long i)
// {
//     return ZZ_to_int(&rep(coeff(*x, i)));
// }

// struct ZZ_pX* ZZ_pX_add(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_pX *z = new ZZ_pX();
//   add(*z, *x, *y);
//   return z;
// }

// struct ZZ_pX* ZZ_pX_sub(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_pX *z = new ZZ_pX();
//   sub(*z, *x, *y);
//   return z;
// }

// struct ZZ_pX* ZZ_pX_mul(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_pX *z = new ZZ_pX();
//   mul(*z, *x, *y);
//   return z;
// }


// struct ZZ_pX* ZZ_pX_div(struct ZZ_pX* x, struct ZZ_pX* y, int* divisible)
// {
//   struct ZZ_pX* z = new ZZ_pX();
//   *divisible = divide(*z, *x, *y);
//   return z;
// }


// struct ZZ_pX* ZZ_pX_mod(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   struct ZZ_pX* z = new ZZ_pX();
//   rem(*z, *x, *y);
//   return z;
// }



// void ZZ_pX_quo_rem(struct ZZ_pX* x, struct ZZ_pX* y, struct ZZ_pX** r, struct ZZ_pX** q)
// {
//   *r = new ZZ_pX();
//   *q = new ZZ_pX();
//   DivRem(**q, **r, *x, *y);
// }


// struct ZZ_pX* ZZ_pX_square(struct ZZ_pX* x)
// {
//   struct ZZ_pX* s = new ZZ_pX();
//   sqr(*s, *x);
//   return s;
// }



// int ZZ_pX_is_monic(struct ZZ_pX* x)
// {
//   IsOne(LeadCoeff(*x));
// }


// struct ZZ_pX* ZZ_pX_neg(struct ZZ_pX* x)
// {
//   struct ZZ_pX* y = new ZZ_pX();
//   *y = -*x;
//   return y;
// }


// struct ZZ_pX* ZZ_pX_left_shift(struct ZZ_pX* x, long n)
// {
//   struct ZZ_pX* y = new ZZ_pX();
//   LeftShift(*y, *x, n);
//   return y;
// }


// struct ZZ_pX* ZZ_pX_right_shift(struct ZZ_pX* x, long n)
// {
//   struct ZZ_pX* y = new ZZ_pX();
//   RightShift(*y, *x, n);
//   return y;
// }



// struct ZZ_pX* ZZ_pX_gcd(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   struct ZZ_pX* g = new ZZ_pX();
//   GCD(*g, *x, *y);
//   return g;
// }


// void ZZ_pX_xgcd(struct ZZ_pX** d, struct ZZ_pX** s, struct ZZ_pX** t, struct ZZ_pX* a, struct ZZ_pX* b)
// {
//   *d = new ZZ_pX();
//   *s = new ZZ_pX();
//   *t = new ZZ_pX();
//   XGCD(**d, **s, **t, *a, *b);
// }

// void ZZ_pX_plain_xgcd(struct ZZ_pX** d, struct ZZ_pX** s, struct ZZ_pX** t, struct ZZ_pX* a, struct ZZ_pX* b)
// {
//   *d = new ZZ_pX();
//   *s = new ZZ_pX();
//   *t = new ZZ_pX();
//   PlainXGCD(**d, **s, **t, *a, *b);
// }

// ZZ_p* ZZ_pX_leading_coefficient(struct ZZ_pX* x)
// {
//   return new ZZ_p(LeadCoeff(*x));
// }


// void ZZ_pX_set_x(struct ZZ_pX* x)
// {
//   SetX(*x);
// }


// int ZZ_pX_is_x(struct ZZ_pX* x)
// {
//   return IsX(*x);
// }


// struct ZZ_pX* ZZ_pX_derivative(struct ZZ_pX* x)
// {
//   ZZ_pX* d = new ZZ_pX();
//   diff(*d, *x);
//   return d;
// }


// struct ZZ_pX* ZZ_pX_reverse(struct ZZ_pX* x)
// {
//   ZZ_pX* r = new ZZ_pX();
//   reverse(*r, *x);
//   return r;
// }

// struct ZZ_pX* ZZ_pX_reverse_hi(struct ZZ_pX* x, int hi)
// {
//   ZZ_pX* r = new ZZ_pX();
//   reverse(*r, *x, hi);
//   return r;
// }


// struct ZZ_pX* ZZ_pX_truncate(struct ZZ_pX* x, long m)
// {
//   ZZ_pX* t = new ZZ_pX();
//   trunc(*t, *x, m);
//   return t;
// }


// struct ZZ_pX* ZZ_pX_multiply_and_truncate(struct ZZ_pX* x, struct ZZ_pX* y, long m)
// {
//   ZZ_pX* t = new ZZ_pX();
//   MulTrunc(*t, *x, *y, m);
//   return t;
// }


// struct ZZ_pX* ZZ_pX_square_and_truncate(struct ZZ_pX* x, long m)
// {
//   ZZ_pX* t = new ZZ_pX();
//   SqrTrunc(*t, *x, m);
//   return t;
// }


// struct ZZ_pX* ZZ_pX_invert_and_truncate(struct ZZ_pX* x, long m)
// {
//   ZZ_pX* t = new ZZ_pX();
//   InvTrunc(*t, *x, m);
//   return t;
// }


// struct ZZ_pX* ZZ_pX_multiply_mod(struct ZZ_pX* x, struct ZZ_pX* y,  struct ZZ_pX* modulus)
// {
//   ZZ_pX* p = new ZZ_pX();
//   MulMod(*p, *x, *y, *modulus);
//   return p;
// }


// struct ZZ_p* ZZ_pX_trace_mod(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_p* p = new ZZ_p();
//   TraceMod(*p, *x, *y);
//   return p;
// }


char* ZZ_pX_trace_list(struct ZZ_pX* x)
{
    vec_ZZ_p v;
    TraceVec(v, *x);
    ostringstream instore;
    instore << v;
    int n = strlen(instore.str().data());
    char* buf = new char[n+1];
    strcpy(buf, instore.str().data());
    return buf;
}


// struct ZZ_p* ZZ_pX_resultant(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_p* res = new ZZ_p();
//   resultant(*res, *x, *y);
//   return res;
// }


// struct ZZ_p* ZZ_pX_norm_mod(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_p* res = new ZZ_p();
//   NormMod(*res, *x, *y);
//   return res;
// }



// struct ZZ_pX* ZZ_pX_charpoly_mod(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_pX* f = new ZZ_pX();
//   CharPolyMod(*f, *x, *y);
//   return f;
// }


// struct ZZ_pX* ZZ_pX_minpoly_mod(struct ZZ_pX* x, struct ZZ_pX* y)
// {
//   ZZ_pX* f = new ZZ_pX();
//   MinPolyMod(*f, *x, *y);
//   return f;
// }


// void ZZ_pX_clear(struct ZZ_pX* x)
// {
//   clear(*x);
// }


// void ZZ_pX_preallocate_space(struct ZZ_pX* x, long n)
// {
//   x->SetMaxLength(n);
// }

void ZZ_pX_factor(struct ZZ_pX*** v, long** e, long* n, struct ZZ_pX* x, long verbose)
{
    long i;
    vec_pair_ZZ_pX_long factors;
    berlekamp(factors, *x, verbose);
    *n = factors.length();
    *v = (ZZ_pX**) malloc(sizeof(ZZ_pX*) * (*n));
    *e = (long*) malloc(sizeof(long)*(*n));
    for (i=0; i<(*n); i++) {
        (*v)[i] = new ZZ_pX(factors[i].a);
        (*e)[i] = factors[i].b;
    }
}

void ZZ_pX_linear_roots(struct ZZ_p*** v, long* n, struct ZZ_pX* f)
{
    long i;
    vec_ZZ_p w;
    FindRoots(w, *f);
    *n = w.length();
    (*v) = (ZZ_p**) malloc(sizeof(ZZ_p*)*(*n));
    for (i=0; i<(*n); i++) {
        (*v)[i] = new ZZ_p(w[i]);
    }
}

/////////// ZZ_pE //////////////

struct ZZ_pX ZZ_pE_to_ZZ_pX(struct ZZ_pE x)
{
    return ZZ_pX(rep(x));
}



//////// mat_ZZ //////////

void mat_ZZ_SetDims(mat_ZZ* mZZ, long nrows, long ncols){
    mZZ->SetDims(nrows, ncols);
}

mat_ZZ* mat_ZZ_pow(const mat_ZZ* x, long e) 
{
    mat_ZZ *z = new mat_ZZ();
    power(*z, *x, e);
    return z;
}

long mat_ZZ_nrows(const mat_ZZ* x)
{
    return x->NumRows();
}


long mat_ZZ_ncols(const mat_ZZ* x)
{
    return x->NumCols();
}

void mat_ZZ_setitem(mat_ZZ* x, int i, int j, const struct ZZ* z)
{
    (*x)[i][j] = *z;

}

struct ZZ* mat_ZZ_getitem(const mat_ZZ* x, int i, int j)
{
    return new ZZ((*x)(i,j));
}

struct ZZ* mat_ZZ_determinant(const mat_ZZ* x, long deterministic)
{
    ZZ* d = new ZZ();
    determinant(*d, *x, deterministic);
    return d;
}

mat_ZZ* mat_ZZ_HNF(const mat_ZZ* A, const struct ZZ* D)
{
    mat_ZZ* W = new mat_ZZ();
    HNF(*W, *A, *D);
    return W;
}

long mat_ZZ_LLL(struct ZZ **det, mat_ZZ *x, long a, long b, long verbose)
{
    *det = new ZZ();
    return LLL(**det,*x,a,b,verbose);
}

long mat_ZZ_LLL_U(struct ZZ **det, mat_ZZ *x, mat_ZZ *U, long a, long b, long verbose)
{
    *det = new ZZ();
    return LLL(**det,*x,*U,a,b,verbose);
}


struct ZZX* mat_ZZ_charpoly(const mat_ZZ* A)
{
    ZZX* f = new ZZX();
    CharPoly(*f, *A);
    return f;
}

/**
 * GF2EContext
 */

GF2EContext* GF2EContext_construct(void *mem, const GF2X *p)
{
    return new(mem) GF2EContext(*p);
}


GF2EContext* GF2EContext_new(const GF2X *p)
{
    return new GF2EContext(*p);
}


void mat_GF2E_setitem(mat_GF2E* x, int i, int j, const struct GF2E* z)
{
    (*x)[i][j] = *z;
}

void mat_GF2_setitem(mat_GF2* x, int i, int j, const struct GF2* z)
{
    (*x)[i][j] = *z;
}

/**
 * ZZ_pContext
 */

ZZ_pContext* ZZ_pContext_new(ZZ *p)
{
    return new ZZ_pContext(*p);
}

ZZ_pContext* ZZ_pContext_construct(void *mem, ZZ *p)
{
    return new(mem) ZZ_pContext(*p);
}

void ZZ_pContext_restore(ZZ_pContext *ctx)
{
    ctx->restore();
}

// Functions for using ZZ_pX's for p-adic extensions

void ZZ_pX_conv_modulus(ZZ_pX &fout, const ZZ_pX &fin, const ZZ_pContext &modout)
{
    // Changes the modulus of fin to modout, and puts the result in fout.
    long i, n;

    n = fin.rep.length();
    fout.rep.SetLength(n);

    ZZ_p* xp = fout.rep.elts();
    const ZZ_p* ap = fin.rep.elts();

    // I think it's enough to just restore modout once.
    // This should be true as long as the function rep taking a ZZ_p as an argument
    // and returning a ZZ works when the ZZ_p::modulus is incorrect.
    modout.restore();

    for (i = 0; i < n; i++)
    {
        conv(xp[i], rep(ap[i]));
    }

    // We may have set a leading coefficient to 0, so we have to normalize
    fout.normalize();
}

void ZZ_pEX_conv_modulus(ZZ_pEX &fout, const ZZ_pEX &fin, const ZZ_pContext &modout)
{
    // Changes the modulus of fin to modout, and puts the result in fout.
    long i, n, j, m;

    n = fin.rep.length();
    fout.rep.SetLength(n);

    ZZ_pE* outpe = fout.rep.elts();
    const ZZ_pE* inpe = fin.rep.elts();

    ZZ_p* xp;
    const ZZ_p* ap;
    // I think it's enough to just restore modout once
    // This should be true as long as Loophole() offers access to
    // the underlying ZZ_pX representations of ZZ_pEs,
    // and rep of a ZZ_p (giving a ZZ) works even if the ZZ_p::modulus is
    // incorrect
    modout.restore();

    for (i = 0; i < n; i++)
    {
        m = rep(inpe[i]).rep.length();
        outpe[i]._ZZ_pE__rep.rep.SetLength(m);

        xp = outpe[i]._ZZ_pE__rep.rep.elts();
        ap = rep(inpe[i]).rep.elts();

        for (j = 0; j < m; j++)
            conv(xp[j], rep(ap[j]));

        // We may have set a leading coefficient to 0, so we have to normalize
        outpe[i]._ZZ_pE__rep.normalize();
    }
    // We may have set a leading coefficient to 0, so we have to normalize
    fout.normalize();
}

void ZZ_pX_min_val_coeff(long & valuation, long &index, const struct ZZ_pX &f, const struct ZZ &p)
{
    // Sets index, where the indexth coefficient of f has the minimum p-adic valuation.
    // Sets valuation to be this valuation.
    // If there are ties, index will be the lowest of the tied indices
    // This only makes mathematical sense when p divides the modulus of f.
    long i, n, v;

    n = f.rep.length();
    if (n == 0)
    {
        index = -1;
        return;
    }

    const ZZ_p* fp = f.rep.elts();
    ZZ *u = new ZZ();

    valuation = -1;
    i = 0;

    while (valuation == -1)
    {
        if (rep(fp[i]) != 0)
        {
            index = i;
            valuation = ZZ_remove(*u, rep(fp[i]), p);
        }
        i++;
    }
    for (; i < n; i++)
    {
        if (rep(fp[i]) != 0)
        {
            v = ZZ_remove(*u, rep(fp[i]), p);
            if (v < valuation)
            {
                valuation = v;
                index = i;
            }
        }
    }
    delete u;
}

long ZZ_pX_get_val_coeff(const struct ZZ_pX &f, const struct ZZ &p, long i)
{
    // Gets the p-adic valuation of the ith coefficient of f.
    ZZ *u = new ZZ();
    long ans = ZZ_remove(*u, rep(coeff(f, i)), p);
    delete u;
    return ans;
}

void ZZ_pX_left_pshift(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ &pn, const struct ZZ_pContext &c)
{
    // Multiplies each coefficient by pn, and sets the context of the answer to c.

    long i, n;

    n = a.rep.length();
    x.rep.SetLength(n);

    ZZ_p* xp = x.rep.elts();
    const ZZ_p* ap = a.rep.elts();

    // I think it's enough to just restore modout once.
    // This should be true as long as the function rep taking a ZZ_p as an argument
    // and returning a ZZ works when the ZZ_p::modulus is incorrect.
    c.restore();

    for (i = 0; i < n; i++)
    {
        conv(xp[i], rep(ap[i]) * pn);
    }

    // We may have set a leading coefficient to 0, so we have to normalize
    x.normalize();
}

void ZZ_pX_right_pshift(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ &pn, const struct ZZ_pContext &c)
{
    // Divides each coefficient by pn, and sets the context of the answer to c.

    long i, n;

    n = a.rep.length();
    x.rep.SetLength(n);

    ZZ_p* xp = x.rep.elts();
    const ZZ_p* ap = a.rep.elts();

    // I think it's enough to just restore modout once.
    // This should be true as long as the function rep taking a ZZ_p as an argument
    // and returning a ZZ works when the ZZ_p::modulus is incorrect.
    c.restore();

    for (i = 0; i < n; i++)
    {
        conv(xp[i], rep(ap[i]) / pn);
    }

    // We may have set a leading coefficient to 0, so we have to normalize
    x.normalize();
}

void ZZ_pX_InvMod_newton_unram(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ_pXModulus &F, const struct ZZ_pContext &cpn, const struct ZZ_pContext &cp)
{
    //int j;
    cp.restore();
    ZZ_pX *amodp = new ZZ_pX();
    ZZ_pX *xmodp = new ZZ_pX();
    ZZ_pX *fmodp = new ZZ_pX();
    ZZ_pX_conv_modulus(*amodp, a, cp);
    ZZ_pX_conv_modulus(*fmodp, F.val(), cp);
    InvMod(*xmodp, *amodp, *fmodp);
    //cout << "xmodp: " << *xmodp << "\namodp: " << *amodp << "\nfmodp: " << *fmodp << "\n";
    cpn.restore();
    ZZ_pX *minusa = new ZZ_pX();
    ZZ_pX *xn = new ZZ_pX();
    ZZ_pX_conv_modulus(*xn, *xmodp, cpn);
    NTL::negate(*minusa, a);
    while (1 > 0)
    {
        // x_n = 2*x_{n-1} - a*x_{n-1}^2 = (2 - a*x_{n-1})*x_{n-1}
        MulMod(x, *minusa, *xn, F);
        SetCoeff(x, 0, ConstTerm(x) + 2);
        MulMod(x, x, *xn, F);
        if (x == *xn)
            break;
        *xn = x;
        //cout << "x: " << x << "\nxn: " << *xn << "\n";
        //cin >> j;
    }
    delete amodp;
    delete xmodp;
    delete fmodp;
    delete minusa;
    delete xn;
}

void ZZ_pX_InvMod_newton_ram(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ_pXModulus &F, const struct ZZ_pContext &cpn)
{
    //int j;
    cpn.restore();
    ZZ_pX *minusa = new ZZ_pX();
    ZZ_pX *xn = new ZZ_pX();
    SetCoeff(*xn, 0, inv(ConstTerm(a)));
    NTL::negate(*minusa, a);
    while (1 > 0)
    {
        // x_n = 2*x_{n-1} - a*x_{n-1}^2 = (2 - a*x_{n-1})*x_{n-1}
        MulMod(x, *minusa, *xn, F);
        SetCoeff(x, 0, ConstTerm(x) + 2);
        MulMod(x, x, *xn, F);
        //cout << "x: " << x << "\nxn: " << *xn << "\n";
        if (x == *xn)
            break;
        *xn = x;
        //cin >> j;
    }
    delete minusa;
    delete xn;
}

/**
 * ZZ_pEContext
 */

ZZ_pEContext* ZZ_pEContext_new(ZZ_pX *f)
{
    return new ZZ_pEContext(*f);
}

ZZ_pEContext* ZZ_pEContext_construct(void *mem, ZZ_pX *f)
{
    return new(mem) ZZ_pEContext(*f);
}

void ZZ_pEContext_restore(ZZ_pEContext *ctx)
{
    ctx->restore();
}

/**
 * zz_pContext
 */

zz_pContext* zz_pContext_new(long p)
{
    return new zz_pContext(p);
}

zz_pContext* zz_pContext_construct(void *mem, long p)
{
    return new(mem) zz_pContext(p);
}

void zz_pContext_restore(zz_pContext *ctx)
{
    ctx->restore();
}
