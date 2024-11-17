#ifndef _SAGE_NTLWRAP_CPP
#define _SAGE_NTLWRAP_CPP
#include "ntlwrap.h"

#ifdef __cplusplus

#include <iostream>
#include <sstream>
#include <gmp.h>

#include <NTL/ZZ_limbs.h>

//////// ZZ //////////

/* Return value is only valid if the result should fit into an int.
   AUTHOR: David Harvey (2008-06-08) */
static CYTHON_INLINE int ZZ_to_int(const ZZ* x)
{
    return to_int(*x);
}

/* Copies the ZZ into the mpz_t
   Assumes output has been mpz_init'd.
   AUTHOR: David Harvey
           Joel B. Mohler moved the ZZX_getitem_as_mpz code out to this function (2007-03-13) */
static void ZZ_to_mpz(mpz_t output, const struct ZZ* x)
{
    mpz_import(output, x->size(), -1, sizeof(mp_limb_t), 0, 0, ZZ_limbs_get(*x));
    if (sign(*x) < 0)
        mpz_neg(output, output);
}

/* Copies the mpz_t into the ZZ
   AUTHOR: Joel B. Mohler (2007-03-15) */
static void mpz_to_ZZ(struct ZZ* output, mpz_srcptr x)
{
    unsigned char stack_bytes[4096];
    size_t size = (mpz_sizeinbase(x, 2) + 7) / 8;
    int use_heap = (size > sizeof(stack_bytes));
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

static long ZZ_remove(struct ZZ &dest, const struct ZZ &src, const struct ZZ &f)
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
        NTL::mul(fpow[p+1], fpow[p], fpow[p]);
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
static CYTHON_INLINE int ZZ_p_to_int(const ZZ_p& x )
{
    return ZZ_to_int(&rep(x));
}


static CYTHON_INLINE void ZZ_p_modulus(struct ZZ* mod, const struct ZZ_p* x)
{
    (*mod) = x->modulus();
}


///////////////////////////////////////////////
//////// ZZX //////////
///////////////////////////////////////////////

static char* ZZX_repr(struct ZZX* x)
{
    std::ostringstream instore;
    instore << (*x);
    int n = strlen(instore.str().data());
    char* buf = new char[n+1];
    strcpy(buf, instore.str().data());
    return buf;
}

static CYTHON_INLINE struct ZZX* ZZX_copy(struct ZZX* x) {
    return new ZZX(*x);
}

/* Sets ith coefficient of x to value.
   AUTHOR: David Harvey (2006-06-08) */
static CYTHON_INLINE void ZZX_setitem_from_int(struct ZZX* x, long i, int value)
{
    SetCoeff(*x, i, value);
}

/* Returns ith coefficient of x.
   Return value is only valid if the result should fit into an int.
   AUTHOR: David Harvey (2006-06-08) */
static CYTHON_INLINE int ZZX_getitem_as_int(struct ZZX* x, long i)
{
    return ZZ_to_int(&coeff(*x, i));
}

/* Copies ith coefficient of x to output.
   Assumes output has been mpz_init'd.
   AUTHOR: David Harvey (2007-02) */
static CYTHON_INLINE void ZZX_getitem_as_mpz(mpz_t output, struct ZZX* x, long i)
{
    const ZZ& z = coeff(*x, i);
    ZZ_to_mpz(output, &z);
}

static CYTHON_INLINE struct ZZX* ZZX_div(struct ZZX* x, struct ZZX* y, int* divisible)
{
    struct ZZX* z = new ZZX();
    *divisible = divide(*z, *x, *y);
    return z;
}


static CYTHON_INLINE void ZZX_quo_rem(struct ZZX* x, struct ZZX* other, struct ZZX** r, struct ZZX** q)
{
    struct ZZX *qq = new ZZX(), *rr = new ZZX();
    DivRem(*qq, *rr, *x, *other);
    *r = rr; *q = qq;
}


static CYTHON_INLINE struct ZZX* ZZX_square(struct ZZX* x)
{
    struct ZZX* s = new ZZX();
    sqr(*s, *x);
    return s;
}


static CYTHON_INLINE struct ZZX* ZZX_neg(struct ZZX* x)
{
    struct ZZX* y = new ZZX();
    *y = -*x;
    return y;
}


static CYTHON_INLINE struct ZZX* ZZX_left_shift(struct ZZX* x, long n)
{
    struct ZZX* y = new ZZX();
    LeftShift(*y, *x, n);
    return y;
}


static CYTHON_INLINE struct ZZX* ZZX_right_shift(struct ZZX* x, long n)
{
    struct ZZX* y = new ZZX();
    RightShift(*y, *x, n);
    return y;
}

static CYTHON_INLINE struct ZZX* ZZX_primitive_part(struct ZZX* x)
{
    struct ZZX* p = new ZZX();
    PrimitivePart(*p, *x);
    return p;
}


static CYTHON_INLINE void ZZX_pseudo_quo_rem(struct ZZX* x, struct ZZX* y, struct ZZX** r, struct ZZX** q)
{
    *r = new ZZX();
    *q = new ZZX();
    PseudoDivRem(**q, **r, *x, *y);
}


static CYTHON_INLINE struct ZZX* ZZX_gcd(struct ZZX* x, struct ZZX* y)
{
    struct ZZX* g = new ZZX();
    GCD(*g, *x, *y);
    return g;
}


static CYTHON_INLINE void ZZX_xgcd(struct ZZX* x, struct ZZX* y, struct ZZ** r, struct ZZX** s,
          struct ZZX** t, int proof)
{
    *r = new ZZ();
    *s = new ZZX();
    *t = new ZZX();
    XGCD(**r, **s, **t, *x, *y, proof);
}


static CYTHON_INLINE void ZZX_set_x(struct ZZX* x)
{
    SetX(*x);
}


static CYTHON_INLINE int ZZX_is_x(struct ZZX* x)
{
    return IsX(*x);
}


static CYTHON_INLINE struct ZZX* ZZX_derivative(struct ZZX* x)
{
    ZZX* d = new ZZX();
    diff(*d, *x);
    return d;
}


static CYTHON_INLINE struct ZZX* ZZX_reverse(struct ZZX* x)
{
    ZZX* r = new ZZX();
    reverse(*r, *x);
    return r;
}

static CYTHON_INLINE struct ZZX* ZZX_reverse_hi(struct ZZX* x, int hi)
{
    ZZX* r = new ZZX();
    reverse(*r, *x, hi);
    return r;
}


static CYTHON_INLINE struct ZZX* ZZX_truncate(struct ZZX* x, long m)
{
    ZZX* t = new ZZX();
    trunc(*t, *x, m);
    return t;
}


static CYTHON_INLINE struct ZZX* ZZX_multiply_and_truncate(struct ZZX* x, struct ZZX* y, long m)
{
    ZZX* t = new ZZX();
    MulTrunc(*t, *x, *y, m);
    return t;
}


static CYTHON_INLINE struct ZZX* ZZX_square_and_truncate(struct ZZX* x, long m)
{
    ZZX* t = new ZZX();
    SqrTrunc(*t, *x, m);
    return t;
}


static CYTHON_INLINE struct ZZX* ZZX_invert_and_truncate(struct ZZX* x, long m)
{
    ZZX* t = new ZZX();
    InvTrunc(*t, *x, m);
    return t;
}


static CYTHON_INLINE struct ZZX* ZZX_multiply_mod(struct ZZX* x, struct ZZX* y,  struct ZZX* modulus)
{
    ZZX* p = new ZZX();
    MulMod(*p, *x, *y, *modulus);
    return p;
}


static CYTHON_INLINE struct ZZ* ZZX_trace_mod(struct ZZX* x, struct ZZX* y)
{
    ZZ* p = new ZZ();
    TraceMod(*p, *x, *y);
    return p;
}


static char* ZZX_trace_list(struct ZZX* x)
{
    vec_ZZ v;
    TraceVec(v, *x);
    std::ostringstream instore;
    instore << v;
    int n = strlen(instore.str().data());
    char* buf = new char[n+1];
    strcpy(buf, instore.str().data());
    return buf;
}


static CYTHON_INLINE struct ZZ* ZZX_resultant(struct ZZX* x, struct ZZX* y, int proof)
{
    ZZ* res = new ZZ();
    resultant(*res, *x, *y, proof);
    return res;
}


static CYTHON_INLINE struct ZZ* ZZX_norm_mod(struct ZZX* x, struct ZZX* y, int proof)
{
    ZZ* res = new ZZ();
    NormMod(*res, *x, *y, proof);
    return res;
}


static CYTHON_INLINE struct ZZ* ZZX_discriminant(struct ZZX* x, int proof)
{
    ZZ* d = new ZZ();
    discriminant(*d, *x, proof);
    return d;
}


static CYTHON_INLINE struct ZZX* ZZX_charpoly_mod(struct ZZX* x, struct ZZX* y, int proof)
{
    ZZX* f = new ZZX();
    CharPolyMod(*f, *x, *y, proof);
    return f;
}


static CYTHON_INLINE struct ZZX* ZZX_minpoly_mod(struct ZZX* x, struct ZZX* y)
{
    ZZX* f = new ZZX();
    MinPolyMod(*f, *x, *y);
    return f;
}


static CYTHON_INLINE void ZZX_clear(struct ZZX* x)
{
    clear(*x);
}


static CYTHON_INLINE void ZZX_preallocate_space(struct ZZX* x, long n)
{
    x->SetMaxLength(n);
}

static void ZZX_squarefree_decomposition(struct ZZX*** v, long** e, long* n, struct ZZX* x)
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

static char* ZZ_pX_trace_list(struct ZZ_pX* x)
{
    vec_ZZ_p v;
    TraceVec(v, *x);
    std::ostringstream instore;
    instore << v;
    int n = strlen(instore.str().data());
    char* buf = new char[n+1];
    strcpy(buf, instore.str().data());
    return buf;
}

static void ZZ_pX_factor(struct ZZ_pX*** v, long** e, long* n, struct ZZ_pX* x, long verbose)
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

static void ZZ_pX_linear_roots(struct ZZ_p*** v, long* n, struct ZZ_pX* f)
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

static CYTHON_INLINE struct ZZ_pX ZZ_pE_to_ZZ_pX(struct ZZ_pE x)
{
    return ZZ_pX(rep(x));
}

//////// mat_ZZ //////////

static CYTHON_INLINE void mat_ZZ_setitem(mat_ZZ* x, int i, int j, const struct ZZ* z)
{
    (*x)[i][j] = *z;

}

static CYTHON_INLINE struct ZZ* mat_ZZ_getitem(const mat_ZZ* x, int i, int j)
{
    return new ZZ((*x)(i,j));
}

static CYTHON_INLINE struct ZZ* mat_ZZ_determinant(const mat_ZZ* x, long deterministic)
{
    ZZ* d = new ZZ();
    determinant(*d, *x, deterministic);
    return d;
}

static CYTHON_INLINE mat_ZZ* mat_ZZ_HNF(const mat_ZZ* A, const struct ZZ* D)
{
    mat_ZZ* W = new mat_ZZ();
    HNF(*W, *A, *D);
    return W;
}

static CYTHON_INLINE long mat_ZZ_LLL(struct ZZ **det, mat_ZZ *x, long a, long b, long verbose)
{
    *det = new ZZ();
    return LLL(**det, *x, a, b, verbose);
}

static CYTHON_INLINE long mat_ZZ_LLL_U(struct ZZ **det, mat_ZZ *x, mat_ZZ *U, long a, long b, long verbose)
{
    *det = new ZZ();
    return LLL(**det, *x, *U, a, b, verbose);
}


static CYTHON_INLINE void mat_GF2E_setitem(mat_GF2E* x, int i, int j, const struct GF2E* z)
{
    (*x)[i][j] = *z;
}

static CYTHON_INLINE void mat_GF2_setitem(mat_GF2* x, int i, int j, const struct GF2* z)
{
    (*x)[i][j] = *z;
}

// Functions for using ZZ_pX's for p-adic extensions

static void ZZ_pX_conv_modulus(ZZ_pX &fout, const ZZ_pX &fin, const ZZ_pContext &modout)
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

static void ZZ_pEX_conv_modulus(ZZ_pEX &fout, const ZZ_pEX &fin, const ZZ_pContext &modout)
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

static void ZZ_pX_min_val_coeff(long & valuation, long &index, const struct ZZ_pX &f, const struct ZZ &p)
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

static CYTHON_INLINE long ZZ_pX_get_val_coeff(const struct ZZ_pX &f, const struct ZZ &p, long i)
{
    // Gets the p-adic valuation of the ith coefficient of f.
    ZZ *u = new ZZ();
    long ans = ZZ_remove(*u, rep(coeff(f, i)), p);
    delete u;
    return ans;
}

static void ZZ_pX_left_pshift(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ &pn, const struct ZZ_pContext &c)
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

static void ZZ_pX_right_pshift(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ &pn, const struct ZZ_pContext &c)
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

static void ZZ_pX_InvMod_newton_unram(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ_pXModulus &F, const struct ZZ_pContext &cpn, const struct ZZ_pContext &cp)
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

static void ZZ_pX_InvMod_newton_ram(struct ZZ_pX &x, const struct ZZ_pX &a, const struct ZZ_pXModulus &F, const struct ZZ_pContext &cpn)
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

#endif  /* #ifdef __cplusplus */
#endif  /* #ifndef  _SAGE_NTLWRAP_CPP */
