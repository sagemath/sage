/* sage_setup: distribution = sagemath-pari
 */
/*
 *  C helper functions for the computation
 *  of p-adic transcendantal functions
 *
 *********************************************/

#include <math.h>
#include <gmp.h>
#include <stdlib.h>
#include <macros.h>  /* cysignals library */

/* p-adic logarithm */
void padiclog(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo) {
    /*  Compute the p-adic logarithm of a,
        which is supposed to be congruent to 1 mod p

        Algorithm:
         1. we raise a at the power p^(v-1) (for a suitable v) in order
            to make it closer to 1
         2. we write the new a as a product
              1/a = (1 - a_0*p^v) (1 - a_1*p^(2*v) (1 - a_2*p^(4*v) ...
            with 0 <= a_i < p^(v*2^i).
         3. we compute each log(1 - a_i*p^(v*2^i)) using Taylor expansion
            and a binary splitting strategy.                                */

    unsigned long i, v, e, N, saveN, Np, tmp, trunc, step;
    double den = log(p);
    mpz_t f, arg, trunc_mod, h, hpow, mpz_tmp, mpz_tmp2, d, inv, mod2;
    mpz_t *num, *denom;

    mpz_init(mpz_tmp);
    mpz_init(mpz_tmp2);
    mpz_init(arg);
    mpz_set_ui(ans, 0);

    mpz_fdiv_r_ui(mpz_tmp, a, p);
    mpz_set(arg, a);

    /* First we make the argument closer to 1 by raising it to the p^(v-1) */
    if (prec < p) {
        v = 0; e = 1;
    } else {
        v = (unsigned long)(log(prec)/den);  // v here is v-1
        e = pow(p,v);
        mpz_mul_ui(mpz_tmp, modulo, e);
        mpz_powm_ui(arg, arg, e, mpz_tmp);
        prec += v;
    }

    /* Where do we need to truncate the Taylor expansion */
    N = prec+v; N /= ++v;                 // note the ++v
    Np = N;
    den *= v;
    while(1) {
        tmp = Np + (unsigned long)(log(N)/den);
        if (tmp == N) break;
        N = tmp;
    }

    /* We allocate memory and initialize variables */
    mpz_init(f); mpz_init(mod2);
    mpz_init(h); mpz_init(hpow);
    mpz_init(d); mpz_init(inv);
    sig_block();
    num = (mpz_t*)malloc(N*sizeof(mpz_t));
    denom = (mpz_t*)malloc(N*sizeof(mpz_t));
    sig_unblock();
    for (i = 0; i < N; i++) {
        mpz_init(num[i]);
        mpz_init(denom[i]);
    }

    trunc = v << 1;
    mpz_init(trunc_mod);
    mpz_ui_pow_ui(trunc_mod, p, trunc);
    while(1) {
        /* We compute f = 1 - a_i*p^((v+1)*2^i)
           trunc_mod is p^((v+1)*2^(i+1)) */
        mpz_fdiv_r(f, arg, trunc_mod);

        if (mpz_cmp_ui(f, 1) != 0) {

            mpz_ui_sub(f, 2, f);
            mpz_mul(arg, arg, f);

            /* We compute the Taylor expansion of log(f)
               For now, computations are carried out over the rationals */
            for (i = 0; i < N; i++) {
                mpz_set_ui(num[i], 1);
                mpz_set_ui(denom[i], i+1);
            }
            step = 1;
            mpz_ui_sub(h, 1, f);   // we write f = 1 - h, i.e. h = a_i*p^(2^i)
            mpz_set(hpow, h);
            while(step < N) {
                for (i = 0; i < N - step; i += step << 1) {
                    mpz_mul(mpz_tmp2, hpow, num[i+step]);
                    mpz_mul(mpz_tmp, mpz_tmp2, denom[i]);
                    mpz_mul(num[i], num[i], denom[i+step]);
                    mpz_add(num[i], num[i], mpz_tmp);
                    mpz_mul(denom[i], denom[i], denom[i+step]);
                }
                step <<= 1;
                mpz_mul(hpow, hpow, hpow);
            }

            /* We simplify the fraction */
            Np = N; tmp = 0;
            while(Np > 0) { Np /= p; tmp += Np; }
            mpz_ui_pow_ui(d, p, tmp);
            mpz_divexact(mpz_tmp, num[0], d);
            mpz_divexact(denom[0], denom[0], d);

            mpz_divexact_ui(h, h, e);
            mpz_mul(mpz_tmp, h, mpz_tmp);

            /* We coerce the result from Q to Zp */
            mpz_gcdext(d, inv, NULL, denom[0], modulo);
            mpz_mul(mpz_tmp, mpz_tmp, inv);

            /* We add this contribution to log(f) */
            mpz_add(ans, ans, mpz_tmp);

        }

        if (trunc > prec) break;

        /* We update the variables for the next step */
        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        trunc <<= 1;
        for (i = N >> 1; i < N; i++) {
            mpz_clear(num[i]);
            mpz_clear(denom[i]);
        }
        N >>= 1;
    }

    mpz_fdiv_r(ans, ans, modulo);

    /* We clear memory */
    mpz_clear(arg);
    mpz_clear(f);
    mpz_clear(trunc_mod);
    mpz_clear(h);
    mpz_clear(hpow);
    mpz_clear(mpz_tmp);
    mpz_clear(d);
    mpz_clear(inv);
    mpz_clear(mod2);
    for (i = 0; i < N; i++) {
        mpz_clear(num[i]);
        mpz_clear(denom[i]);
    }
    sig_block();
    free(num);
    free(denom);
    sig_unblock();
}


/* p-adic exponential */
void padicexp(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo) {
    /*  Compute the p-adic exponential of a, which is supposed
         - to be congruent to 0 mod p if p > 2
         - to be congruent to 0 mod 4 if p = 2

        Algorithm:
         1. we write a as a sum
              a = a_0*p + a_1*p^2 + a_2*p^4 + ...
            with 0 <= a_i < p^(2^i).
         2. we compute each exp(a_i*p^(2^i)) using Taylor expansion
            and a binary splitting strategy.                           */

    unsigned long i, N, saveN, Np, tmp, trunc, step;
    mpz_t f, arg, trunc_mod, h, hpow, mpz_tmp, d, inv;
    mpz_t denominator;
    mpz_t *num, *denom;

    mpz_init(mpz_tmp);
    mpz_init(arg);
    mpz_set_ui(ans, 1);
    mpz_init(denominator);
    mpz_set_ui(denominator, 1);
    mpz_set(arg,a);

    /* Where do we need to truncate the Taylor expansion */
    if (p == 2) {
        N = prec;
    } else {
        N = (prec*(p-1)) / (p-2);
    }
    saveN = N;

    /* We allocate memory and initialize variables */
    mpz_init(f);
    mpz_init(h); mpz_init(hpow);
    mpz_init(d); mpz_init(inv);
    sig_block();
    num = (mpz_t*)malloc((N+1)*sizeof(mpz_t));
    denom = (mpz_t*)malloc((N+1)*sizeof(mpz_t));
    sig_unblock();
    for (i = 0; i <= N; i++) {
        mpz_init(num[i]);
        mpz_init(denom[i]);
    }

    if (p == 2) {
        trunc = 4;
        mpz_init_set_ui(trunc_mod, p);
        mpz_mul_ui(trunc_mod, trunc_mod, p);
        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
    } else {
        trunc = 2;
        mpz_init_set_ui(trunc_mod, p);
        mpz_mul_ui(trunc_mod, trunc_mod, p);
    }
    while(1) {
        mpz_fdiv_r(f, arg, trunc_mod);
        mpz_sub(arg, arg, f);

        if (mpz_cmp_ui(f, 0) != 0) {

            /* We compute the Taylor expansion of exp(f)
               For now, computations are carried out over the rationals */
            mpz_set_ui(num[0], 1);
            mpz_set_ui(denom[0], 1);
            for (i = 1; i <= N; i++) {
                mpz_set_ui(num[i], 1);
                mpz_set_ui(denom[i], i);
            }
            step = 1;
            mpz_set(h, f);
            mpz_set(hpow, h);

            while(1) {
                for (i = 0; i <= N - step; i += step << 1) {
                    mpz_mul(mpz_tmp, hpow, num[i+step]);
                    mpz_mul(num[i], num[i], denom[i+step]);
                    mpz_add(num[i], num[i], mpz_tmp);
                    mpz_mul(denom[i], denom[i], denom[i+step]);
                }
                step <<= 1;
                if (step > N) break;
                mpz_mul(hpow, hpow, hpow);
            }

            /* We simplify the fraction */
            Np = N; tmp = 0;
            while(Np > 0) { Np /= p; tmp += Np; }
            mpz_ui_pow_ui(d, p, tmp);
            mpz_divexact(num[0], num[0], d);
            mpz_divexact(denom[0], denom[0], d);

            /* We add this contribution to exp(f) */
            mpz_mul(ans, ans, num[0]);
            mpz_fdiv_r(ans, ans, modulo);
            mpz_mul(denominator, denominator, denom[0]);
            mpz_fdiv_r(denominator, denominator, modulo);

        }

        if (trunc > prec) break;

        /* We update the variables for the next step */
        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        trunc <<= 1;
        N >>= 1;
    }

    /* We coerce the result from Q to Zp */
    mpz_gcdext(d, inv, NULL, denominator, modulo);
    mpz_mul(ans, ans, inv);
    mpz_fdiv_r(ans, ans, modulo);

    /* We clear memory */
    mpz_clear(arg);
    mpz_clear(denominator);
    mpz_clear(f);
    mpz_clear(trunc_mod);
    mpz_clear(h);
    mpz_clear(hpow);
    mpz_clear(mpz_tmp);
    mpz_clear(d);
    mpz_clear(inv);
    for (i = 0; i <= saveN; i++) {
        mpz_clear(num[i]);
        mpz_clear(denom[i]);
    }
    sig_block();
    free(num);
    free(denom);
    sig_unblock();
}


void padicexp_Newton(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, unsigned long precinit, const mpz_t modulo) {
    /*  Compute the p-adic exponential of a,
        which is supposed to be congruent to 0 mod p

        Algorithm:
        We solve the equation log(ans) = a using Newton iteration
        and apply a binary splitting strategy for computing the log   */

    unsigned long i, N, saveN, Np, tmp, trunc, step;
    mpz_t f, arg, li, bi, trunc_mod, h, hpow, mpz_tmp;
    mpz_t d, inv;
    mpz_t *num, *denom;

    N = 1 + prec;
    while(1) {
      tmp = 1 + prec + (unsigned long)(log(N)/log(p));
      if (tmp == N) break;
      N = tmp;
    }
    saveN = N;

    mpz_init(mpz_tmp);
    mpz_init(arg);
    mpz_set(arg, ans);
    mpz_set_ui(ans, 1);
    mpz_init(li);

    // We first compute l_1 = log(ans) but
    // stop the computation at some finite level
    // We also update ans so that ans = exp(l_1)

    trunc = 2;
    mpz_init_set_ui(trunc_mod, p);
    mpz_mul_ui(trunc_mod, trunc_mod, p);
    mpz_init(f);
    mpz_init(h); mpz_init(hpow);
    mpz_init(d); mpz_init(inv);
    sig_block();
    num = (mpz_t*)malloc(N*sizeof(mpz_t));
    denom = (mpz_t*)malloc(N*sizeof(mpz_t));
    sig_unblock();
    for (i = 0; i < N; i++) {
        mpz_init(num[i]);
        mpz_init(denom[i]);
    }

    while(1) {
        mpz_fdiv_r(f, arg, trunc_mod);

        if (mpz_cmp_ui(f, 1) != 0) {

            mpz_mul(ans, ans, f);
            mpz_fdiv_r(ans, ans, modulo);
            mpz_ui_sub(f, 2, f);
            mpz_mul(arg, arg, f);

            for (i = 0; i < N; i++) {
                mpz_set_ui(num[i], 1);
                mpz_set_ui(denom[i], i+1);
            }
            step = 1;
            mpz_ui_sub(h, 1, f);
            mpz_set(hpow, h);
            while(1) {
                for (i = 0; i < N - step; i += step << 1) {
                    mpz_mul(mpz_tmp, hpow, num[i+step]);
                    mpz_mul(mpz_tmp, mpz_tmp, denom[i]);
                    mpz_mul(num[i], num[i], denom[i+step]);
                    mpz_add(num[i], num[i], mpz_tmp);
                    mpz_mul(denom[i], denom[i], denom[i+step]);
                }
                step <<= 1;
                if (step >= N) break;
                mpz_mul(hpow, hpow, hpow);
            }

            Np = N; tmp = 0;
            while(Np > 0) { Np /= p; tmp += Np; }
            mpz_ui_pow_ui(d, p, tmp);
            mpz_divexact(mpz_tmp, num[0], d);
            mpz_mul(mpz_tmp, h, mpz_tmp);

            mpz_divexact(denom[0], denom[0], d);
            mpz_gcdext(d, inv, NULL, denom[0], modulo);
            mpz_mul(mpz_tmp, mpz_tmp, inv);

            mpz_add(li, li, mpz_tmp);

        }

        if (trunc > precinit) break;

        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        trunc <<= 1;
        N >>= 1;
    }

    mpz_gcdext(d, inv, NULL, ans, modulo);
    mpz_mul(ans, ans, inv);

    // Here comes the Newton iteration
    //   a_(i+1) = a_i * (1 + b_i)
    //   l_(i+1) = l_i + log(1 + b_i)
    //   b_(i+1) = (a - l_(i+1)) correctly truncated
    // NB: the value of a_i is stored in the variable ans

    // Initialization of b_1
    N = 1 + prec/precinit;
    while(1) {
        tmp = 1 + prec/precinit + (unsigned long)(log(N)/log(p));
        if (tmp == N) break;
        N = tmp;
    }
    if (p == 2) trunc = 2*precinit - 1;
    else trunc = precinit << 1;

    mpz_ui_pow_ui(trunc_mod, p, trunc);
    mpz_init(bi);
    mpz_sub(bi, a, li);
    mpz_fdiv_r(bi, bi, trunc_mod);

    while(1) {
        if (mpz_cmp_ui(bi, 0) != 0) {

            // We set a_(i+1) = a_i * (1 + b_i)
            mpz_add_ui(mpz_tmp, bi, 1);
            mpz_mul(ans, ans, mpz_tmp);
            mpz_fdiv_r(ans, ans, modulo);

            // We compute l_(i+1) = l_i + log(1 + b_i)
            for (i = 0; i < N; i++) {
                mpz_set_ui(num[i], 1);
                mpz_set_ui(denom[i], i+1);
            }
            step = 1;
            mpz_neg(h, bi);
            mpz_set(hpow, h);
            while(1) {
                for (i = 0; i < N - step; i += step << 1) {
                    mpz_mul(mpz_tmp, hpow, num[i+step]);
                    mpz_mul(mpz_tmp, mpz_tmp, denom[i]);
                    mpz_mul(num[i], num[i], denom[i+step]);
                    mpz_add(num[i], num[i], mpz_tmp);
                    mpz_mul(denom[i], denom[i], denom[i+step]);
                }
                step <<= 1;
                if (step >= N) break;
                mpz_mul(hpow, hpow, hpow);
            }

            Np = N; tmp = 0;
            while(Np > 0) { Np /= p; tmp += Np; }
            mpz_ui_pow_ui(d, p, tmp);
            mpz_divexact(mpz_tmp, num[0], d);
            mpz_mul(mpz_tmp, h, mpz_tmp);

            mpz_divexact(denom[0], denom[0], d);
            mpz_gcdext(d, inv, NULL, denom[0], modulo);
            mpz_mul(mpz_tmp, mpz_tmp, inv);

            mpz_sub(li, li, mpz_tmp);
        }

        if (trunc > prec) break;

        // We update N, trunc and trunc_mod
        if (p == 2) {
            N = 1 + prec/trunc;
            while(1) {
                tmp = 1 + prec/trunc + (unsigned long)(log(N)/log(p));
                if (tmp == N) break;
                N = tmp;
            }
            trunc = 2*trunc - 1;
            mpz_mul(trunc_mod, trunc_mod, trunc_mod);
            mpz_divexact_ui(trunc_mod, trunc_mod, p);
        } else {
            trunc <<= 1;
            N >>= 1;
            mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        }

        // We set b_(i+1) = a - l_(i+1)  mod  trunc_mod
        mpz_sub(bi, a, li);
        mpz_fdiv_r(bi, bi, trunc_mod);

    }

    mpz_fdiv_r(ans, ans, modulo);

    /* We clear memory */
    mpz_clear(arg);
    mpz_clear(f);
    mpz_clear(trunc_mod);
    mpz_clear(h);
    mpz_clear(hpow);
    mpz_clear(mpz_tmp);
    mpz_clear(d);
    mpz_clear(inv);
    mpz_clear(li);
    mpz_clear(bi);
    for (i = 0; i < saveN; i++) {
        mpz_clear(num[i]);
        mpz_clear(denom[i]);
    }
    sig_block();
    free(num);
    free(denom);
    sig_unblock();
}
