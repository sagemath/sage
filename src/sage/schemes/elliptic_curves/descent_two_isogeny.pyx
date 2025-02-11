r"""
Descent on elliptic curves over `\QQ` with a 2-isogeny
"""

# ***************************************************************************
#       Copyright (C) 2009 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
cdef object x_ZZ = polygen(ZZ)
from sage.rings.polynomial.real_roots import real_roots
from sage.arith.misc import prime_divisors
import sage.libs.ntl.all as ntl

from sage.rings.integer cimport Integer
from sage.libs.gmp.mpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpz_poly_sage cimport *
from sage.libs.flint.nmod_poly cimport *
from sage.libs.flint.nmod_poly_factor cimport *
from sage.libs.flint.ulong_extras cimport *

from cypari2.paridecl cimport (GEN, cgetg, t_POL, set_gel, gel, stoi, lg,
                               evalvarn, evalsigne, Z_issquare,
                               hyperellratpoints)
from cypari2.stack cimport clear_stack
from sage.libs.pari.convert_gmp cimport _new_GEN_from_mpz_t


DEF N_RES_CLASSES_BSD = 10


cdef unsigned long valuation(mpz_t a, mpz_t p) noexcept:
    """
    Return the number of times p divides a.
    """
    cdef mpz_t aa
    cdef unsigned long v
    mpz_init(aa)
    v = mpz_remove(aa, a, p)
    mpz_clear(aa)
    return v


def test_valuation(a, p):
    """
    Doctest function for cdef long valuation(mpz_t, mpz_t).

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.descent_two_isogeny import test_valuation as tv
        sage: for i in [1..20]:
        ....:     print('{:>10} {} {} {}'.format(str(factor(i)), tv(i,2), tv(i,3), tv(i,5)))
                 1 0 0 0
                 2 1 0 0
                 3 0 1 0
               2^2 2 0 0
                 5 0 0 1
             2 * 3 1 1 0
                 7 0 0 0
               2^3 3 0 0
               3^2 0 2 0
             2 * 5 1 0 1
                11 0 0 0
           2^2 * 3 2 1 0
                13 0 0 0
             2 * 7 1 0 0
             3 * 5 0 1 1
               2^4 4 0 0
                17 0 0 0
           2 * 3^2 1 2 0
                19 0 0 0
           2^2 * 5 2 0 1
    """
    cdef Integer A = Integer(a)
    cdef Integer P = Integer(p)
    return valuation(A.value, P.value)


cdef int padic_square(mpz_t a, mpz_t p) noexcept:
    """
    Test if a is a `p`-adic square.
    """
    cdef unsigned long v
    cdef mpz_t aa
    cdef int result

    if mpz_sgn(a) == 0:
        return 1

    v = valuation(a, p)
    if v & 1:
        return 0

    mpz_init_set(aa, a)
    while v:
        v -= 1
        mpz_divexact(aa, aa, p)
    if mpz_cmp_ui(p, 2) == 0:
        result = bool(mpz_fdiv_ui(aa, 8) == 1)
    else:
        result = bool(mpz_legendre(aa, p) == 1)
    mpz_clear(aa)
    return result


def test_padic_square(a, p):
    """
    Doctest function for cdef int padic_square(mpz_t, unsigned long).

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.descent_two_isogeny import test_padic_square as ps
        sage: for i in [1..300]:
        ....:     for p in prime_range(100):
        ....:          if Qp(p)(i).is_square() != bool(ps(i,p)):
        ....:              print(i, p)
    """
    cdef Integer A = Integer(a)
    cdef Integer P = Integer(p)
    return padic_square(A.value, P.value)


cdef int lemma6(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e,
                mpz_t x, mpz_t p, unsigned long nu) noexcept:
    """
    Implement Lemma 6 of BSD's "Notes on elliptic curves, I" for odd `p`.

    Returns -1 for insoluble, 0 for undecided, +1 for soluble.
    """
    cdef mpz_t g_of_x, g_prime_of_x
    cdef unsigned long lambd, mu
    cdef int result = -1

    mpz_init(g_of_x)
    mpz_mul(g_of_x, a, x)
    mpz_add(g_of_x, g_of_x, b)
    mpz_mul(g_of_x, g_of_x, x)
    mpz_add(g_of_x, g_of_x, c)
    mpz_mul(g_of_x, g_of_x, x)
    mpz_add(g_of_x, g_of_x, d)
    mpz_mul(g_of_x, g_of_x, x)
    mpz_add(g_of_x, g_of_x, e)

    if padic_square(g_of_x, p):
        mpz_clear(g_of_x)
        return +1  # soluble

    mpz_init_set(g_prime_of_x, x)
    mpz_mul(g_prime_of_x, a, x)
    mpz_mul_ui(g_prime_of_x, g_prime_of_x, 4)
    mpz_addmul_ui(g_prime_of_x, b, 3)
    mpz_mul(g_prime_of_x, g_prime_of_x, x)
    mpz_addmul_ui(g_prime_of_x, c, 2)
    mpz_mul(g_prime_of_x, g_prime_of_x, x)
    mpz_add(g_prime_of_x, g_prime_of_x, d)

    lambd = valuation(g_of_x, p)
    if mpz_sgn(g_prime_of_x) == 0:
        if lambd >= 2*nu:
            result = 0  # undecided
    else:
        mu = valuation(g_prime_of_x, p)
        if lambd > 2*mu:
            result = 1  # soluble
        elif lambd >= 2*nu and mu >= nu:
            result = 0  # undecided

    mpz_clear(g_prime_of_x)
    mpz_clear(g_of_x)
    return result


cdef int lemma7(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e,
                mpz_t x, mpz_t p, unsigned long nu) noexcept:
    """
    Implement Lemma 7 of BSD's "Notes on elliptic curves, I" for `p=2`.

    Returns -1 for insoluble, 0 for undecided, +1 for soluble.
    """
    cdef mpz_t g_of_x, g_prime_of_x, g_of_x_odd_part
    cdef unsigned long lambd, mu, g_of_x_odd_part_mod_4
    cdef int result = -1

    mpz_init(g_of_x)
    mpz_mul(g_of_x, a, x)
    mpz_add(g_of_x, g_of_x, b)
    mpz_mul(g_of_x, g_of_x, x)
    mpz_add(g_of_x, g_of_x, c)
    mpz_mul(g_of_x, g_of_x, x)
    mpz_add(g_of_x, g_of_x, d)
    mpz_mul(g_of_x, g_of_x, x)
    mpz_add(g_of_x, g_of_x, e)

    if padic_square(g_of_x, p):
        mpz_clear(g_of_x)
        return +1  # soluble

    mpz_init_set(g_prime_of_x, x)
    mpz_mul(g_prime_of_x, a, x)
    mpz_mul_ui(g_prime_of_x, g_prime_of_x, 4)
    mpz_addmul_ui(g_prime_of_x, b, 3)
    mpz_mul(g_prime_of_x, g_prime_of_x, x)
    mpz_addmul_ui(g_prime_of_x, c, 2)
    mpz_mul(g_prime_of_x, g_prime_of_x, x)
    mpz_add(g_prime_of_x, g_prime_of_x, d)

    lambd = valuation(g_of_x, p)
    mpz_init_set(g_of_x_odd_part, g_of_x)
    while mpz_even_p(g_of_x_odd_part):
        mpz_divexact_ui(g_of_x_odd_part, g_of_x_odd_part, 2)
    g_of_x_odd_part_mod_4 = mpz_fdiv_ui(g_of_x_odd_part, 4)
    if mpz_sgn(g_prime_of_x) == 0:
        if lambd >= 2*nu:
            result = 0  # undecided
        elif lambd == 2*nu-2 and g_of_x_odd_part_mod_4==1:
            result = 0  # undecided
    else:
        mu = valuation(g_prime_of_x, p)
        if lambd > 2*mu:
            result = 1  # soluble
        elif nu > mu:
            if lambd >= mu+nu:
                result = 1  # soluble
            elif lambd+1 == mu+nu and (lambd & 1) == 0:
                result = 1  # soluble
            elif lambd+2 == mu+nu and (lambd & 1) == 0 and g_of_x_odd_part_mod_4 == 1:
                result = 1  # soluble
        else:  # nu <= mu
            if lambd >= 2*nu:
                result = 0  # undecided
            elif lambd+2 == 2*nu and g_of_x_odd_part_mod_4==1:
                result = 0  # undecided

    mpz_clear(g_prime_of_x)
    mpz_clear(g_of_x)
    return result


cdef int Zp_soluble_BSD(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e,
                        mpz_t x_k, mpz_t p, unsigned long k) noexcept:
    """
    Uses the approach of BSD's "Notes on elliptic curves, I" to test for
    solubility of y^2 == ax^4 + bx^3 + cx^2 + dx + e over Zp, with
    x=x_k (mod p^k).
    """
    # returns solubility of y^2 = ax^4 + bx^3 + cx^2 + dx + e
    # in Zp with x=x_k (mod p^k)
    cdef int code
    cdef unsigned long t
    cdef mpz_t s

    if mpz_cmp_ui(p, 2) == 0:
        code = lemma7(a, b, c, d, e, x_k, p, k)
    else:
        code = lemma6(a, b, c, d, e, x_k, p, k)
    if code == 1:
        return 1
    if code == -1:
        return 0

    # now code == 0
    t = 0
    mpz_init(s)
    while code == 0 and mpz_cmp_ui(p, t) > 0 and t < N_RES_CLASSES_BSD:
        mpz_pow_ui(s, p, k)
        mpz_mul_ui(s, s, t)
        mpz_add(s, s, x_k)
        code = Zp_soluble_BSD(a, b, c, d, e, s, p, k+1)
        t += 1
    mpz_clear(s)
    return code


cdef bint Zp_soluble_siksek(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e,
                            mpz_t pp, unsigned long pp_ui,
                            nmod_poly_factor_t f_factzn, nmod_poly_t f,
                            fmpz_poly_t f1, fmpz_poly_t linear) noexcept:
    """
    Uses the approach of Algorithm 5.3.1 of Siksek's thesis to test for
    solubility of y^2 == ax^4 + bx^3 + cx^2 + dx + e over Zp.
    """
    cdef unsigned long v_min, v
    cdef unsigned long roots[4]
    cdef int i, j, has_roots, has_single_roots
    cdef bint result

    cdef mpz_t aa, bb, cc, dd, ee
    cdef mpz_t aaa, bbb, ccc, ddd, eee
    cdef unsigned long qq
    cdef unsigned long rr, ss
    cdef mpz_t tt

    # Step 0: divide out all common p from the quartic
    v_min = valuation(a, pp)
    if mpz_cmp_ui(b, 0) != 0:
        v = valuation(b, pp)
        if v < v_min:
            v_min = v
    if mpz_cmp_ui(c, 0) != 0:
        v = valuation(c, pp)
        if v < v_min:
            v_min = v
    if mpz_cmp_ui(d, 0) != 0:
        v = valuation(d, pp)
        if v < v_min:
            v_min = v
    if mpz_cmp_ui(e, 0) != 0:
        v = valuation(e, pp)
        if v < v_min:
            v_min = v
    for v in range(v_min):
        mpz_divexact(a, a, pp)
        mpz_divexact(b, b, pp)
        mpz_divexact(c, c, pp)
        mpz_divexact(d, d, pp)
        mpz_divexact(e, e, pp)

    if not v_min % 2:
        # Step I in Alg. 5.3.1 of Siksek's thesis
        nmod_poly_set_coeff_ui(f, 0, mpz_fdiv_ui(e, pp_ui))
        nmod_poly_set_coeff_ui(f, 1, mpz_fdiv_ui(d, pp_ui))
        nmod_poly_set_coeff_ui(f, 2, mpz_fdiv_ui(c, pp_ui))
        nmod_poly_set_coeff_ui(f, 3, mpz_fdiv_ui(b, pp_ui))
        nmod_poly_set_coeff_ui(f, 4, mpz_fdiv_ui(a, pp_ui))

        result = 0
        (<nmod_poly_factor_struct *>f_factzn)[0].num = 0  # reset data struct
        qq = nmod_poly_factor(f_factzn, f)
        for i in range(f_factzn.num):
            if f_factzn.exp[i]&1:
                result = 1
                break
        if result == 0 and n_jacobi(qq, pp_ui) == 1:
            result = 1
        if result:
            return 1

        nmod_poly_zero(f)
        nmod_poly_set_coeff_ui(f, 0, 1)
        for i in range(f_factzn.num):
            for j in range(f_factzn.exp[i]>>1):
                nmod_poly_mul(f, f, &f_factzn.p[i])

        (<nmod_poly_factor_struct *>f_factzn)[0].num = 0  # reset data struct
        nmod_poly_factor(f_factzn, f)
        has_roots = 0
        j = 0
        for i in range(f_factzn.num):
            if nmod_poly_degree(&f_factzn.p[i]) == 1 and 0 != nmod_poly_get_coeff_ui(&f_factzn.p[i], 1):
                has_roots = 1
                roots[j] = pp_ui - nmod_poly_get_coeff_ui(&f_factzn.p[i], 0)
                j += 1
        if not has_roots:
            return 0

        i = nmod_poly_degree(f)
        mpz_init(aaa)
        mpz_init(bbb)
        mpz_init(ccc)
        mpz_init(ddd)
        mpz_init(eee)

        if i == 0:  # g == 1
            mpz_set(aaa, a)
            mpz_set(bbb, b)
            mpz_set(ccc, c)
            mpz_set(ddd, d)
            mpz_sub_ui(eee, e, qq)
        elif i == 1:  # g == x + rr
            mpz_set(aaa, a)
            mpz_set(bbb, b)
            mpz_sub_ui(ccc, c, qq)
            rr = nmod_poly_get_coeff_ui(f, 0)
            ss = rr*qq
            mpz_set(ddd,d)
            mpz_sub_ui(ddd, ddd, ss*2)
            mpz_set(eee,e)
            mpz_sub_ui(eee, eee, ss*rr)
        elif i == 2:  # g == x^2 + rr*x + ss
            mpz_sub_ui(aaa, a, qq)
            rr = nmod_poly_get_coeff_ui(f, 1)
            mpz_init(tt)
            mpz_set_ui(tt, rr*qq)
            mpz_set(bbb,b)
            mpz_submul_ui(bbb, tt, 2)
            mpz_set(ccc,c)
            mpz_submul_ui(ccc, tt, rr)
            ss = nmod_poly_get_coeff_ui(f, 0)
            mpz_set_ui(tt, ss*qq)
            mpz_set(eee,e)
            mpz_submul_ui(eee, tt, ss)
            mpz_mul_ui(tt, tt, 2)
            mpz_sub(ccc, ccc, tt)
            mpz_set(ddd,d)
            mpz_submul_ui(ddd, tt, rr)
            mpz_clear(tt)
        mpz_divexact(aaa, aaa, pp)
        mpz_divexact(bbb, bbb, pp)
        mpz_divexact(ccc, ccc, pp)
        mpz_divexact(ddd, ddd, pp)
        mpz_divexact(eee, eee, pp)
        # now aaa, bbb, ccc, ddd, eee represents h(x)

        result = 0
        mpz_init(tt)
        for i in range(j):
            mpz_mul_ui(tt, aaa, roots[i])
            mpz_add(tt, tt, bbb)
            mpz_mul_ui(tt, tt, roots[i])
            mpz_add(tt, tt, ccc)
            mpz_mul_ui(tt, tt, roots[i])
            mpz_add(tt, tt, ddd)
            mpz_mul_ui(tt, tt, roots[i])
            mpz_add(tt, tt, eee)
            # tt == h(r) mod p
            mpz_mod(tt, tt, pp)
            if mpz_sgn(tt) == 0:
                fmpz_poly_zero(f1)
                fmpz_poly_zero(linear)
                fmpz_poly_set_coeff_mpz(f1, 0, e)
                fmpz_poly_set_coeff_mpz(f1, 1, d)
                fmpz_poly_set_coeff_mpz(f1, 2, c)
                fmpz_poly_set_coeff_mpz(f1, 3, b)
                fmpz_poly_set_coeff_mpz(f1, 4, a)
                fmpz_poly_set_coeff_ui(linear, 0, roots[i])
                fmpz_poly_set_coeff_mpz(linear, 1, pp)
                fmpz_poly_compose(f1, f1, linear)
                fmpz_poly_scalar_fdiv_ui(f1, f1, pp_ui)
                fmpz_poly_scalar_fdiv_ui(f1, f1, pp_ui)
                mpz_init(aa)
                mpz_init(bb)
                mpz_init(cc)
                mpz_init(dd)
                mpz_init(ee)
                fmpz_poly_get_coeff_mpz(aa, f1, 4)
                fmpz_poly_get_coeff_mpz(bb, f1, 3)
                fmpz_poly_get_coeff_mpz(cc, f1, 2)
                fmpz_poly_get_coeff_mpz(dd, f1, 1)
                fmpz_poly_get_coeff_mpz(ee, f1, 0)
                result = Zp_soluble_siksek(aa, bb, cc, dd, ee, pp, pp_ui, f_factzn, f, f1, linear)
                mpz_clear(aa)
                mpz_clear(bb)
                mpz_clear(cc)
                mpz_clear(dd)
                mpz_clear(ee)
                if result == 1:
                    break
        mpz_clear(aaa)
        mpz_clear(bbb)
        mpz_clear(ccc)
        mpz_clear(ddd)
        mpz_clear(eee)
        mpz_clear(tt)
        return result
    else:
        # Step II in Alg. 5.3.1 of Siksek's thesis
        nmod_poly_set_coeff_ui(f, 0, mpz_fdiv_ui(e, pp_ui))
        nmod_poly_set_coeff_ui(f, 1, mpz_fdiv_ui(d, pp_ui))
        nmod_poly_set_coeff_ui(f, 2, mpz_fdiv_ui(c, pp_ui))
        nmod_poly_set_coeff_ui(f, 3, mpz_fdiv_ui(b, pp_ui))
        nmod_poly_set_coeff_ui(f, 4, mpz_fdiv_ui(a, pp_ui))
        (<nmod_poly_factor_struct *>f_factzn)[0].num = 0  # reset data struct
        nmod_poly_factor(f_factzn, f)
        has_roots = 0
        has_single_roots = 0
        j = 0
        for i in range(f_factzn.num):
            if nmod_poly_degree(&f_factzn.p[i]) == 1 and 0 != nmod_poly_get_coeff_ui(&f_factzn.p[i], 1):
                has_roots = 1
                if f_factzn.exp[i] == 1:
                    has_single_roots = 1
                    break
                roots[j] = pp_ui - nmod_poly_get_coeff_ui(&f_factzn.p[i], 0)
                j += 1

        if not has_roots:
            return 0
        if has_single_roots:
            return 1

        result = 0
        if j > 0:
            mpz_init(aa)
            mpz_init(bb)
            mpz_init(cc)
            mpz_init(dd)
            mpz_init(ee)
        for i in range(j):
            fmpz_poly_zero(f1)
            fmpz_poly_zero(linear)
            fmpz_poly_set_coeff_mpz(f1, 0, e)
            fmpz_poly_set_coeff_mpz(f1, 1, d)
            fmpz_poly_set_coeff_mpz(f1, 2, c)
            fmpz_poly_set_coeff_mpz(f1, 3, b)
            fmpz_poly_set_coeff_mpz(f1, 4, a)
            fmpz_poly_set_coeff_ui(linear, 0, roots[i])
            fmpz_poly_set_coeff_mpz(linear, 1, pp)
            fmpz_poly_compose(f1, f1, linear)
            fmpz_poly_scalar_fdiv_ui(f1, f1, pp_ui)
            fmpz_poly_get_coeff_mpz(aa, f1, 4)
            fmpz_poly_get_coeff_mpz(bb, f1, 3)
            fmpz_poly_get_coeff_mpz(cc, f1, 2)
            fmpz_poly_get_coeff_mpz(dd, f1, 1)
            fmpz_poly_get_coeff_mpz(ee, f1, 0)
            result = Zp_soluble_siksek(aa, bb, cc, dd, ee, pp, pp_ui, f_factzn, f, f1, linear)
            if result == 1:
                break
        if j > 0:
            mpz_clear(aa)
            mpz_clear(bb)
            mpz_clear(cc)
            mpz_clear(dd)
            mpz_clear(ee)
        return result

cdef bint Zp_soluble_siksek_large_p(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e, mpz_t pp,
                                    fmpz_poly_t f1, fmpz_poly_t linear) noexcept:
    """
    Uses the approach of Algorithm 5.3.1 of Siksek's thesis to test for
    solubility of y^2 == ax^4 + bx^3 + cx^2 + dx + e over Zp.
    """
    cdef unsigned long v_min, v
    cdef mpz_t roots[4]
    cdef int i, j, has_roots, has_single_roots
    cdef bint result

    cdef mpz_t aa, bb, cc, dd, ee
    cdef mpz_t aaa, bbb, ccc, ddd, eee
    cdef mpz_t qq, rr, ss, tt
    cdef Integer A, B, C, D, E, P

    # Step 0: divide out all common p from the quartic
    v_min = valuation(a, pp)
    if mpz_cmp_ui(b, 0) != 0:
        v = valuation(b, pp)
        if v < v_min:
            v_min = v
    if mpz_cmp_ui(c, 0) != 0:
        v = valuation(c, pp)
        if v < v_min:
            v_min = v
    if mpz_cmp_ui(d, 0) != 0:
        v = valuation(d, pp)
        if v < v_min:
            v_min = v
    if mpz_cmp_ui(e, 0) != 0:
        v = valuation(e, pp)
        if v < v_min:
            v_min = v
    for v in range(v_min):
        mpz_divexact(a, a, pp)
        mpz_divexact(b, b, pp)
        mpz_divexact(c, c, pp)
        mpz_divexact(d, d, pp)
        mpz_divexact(e, e, pp)

    if not v_min % 2:
        # Step I in Alg. 5.3.1 of Siksek's thesis
        A = Integer(0)
        B = Integer(0)
        C = Integer(0)
        D = Integer(0)
        E = Integer(0)
        P = Integer(0)
        mpz_set(A.value, a)
        mpz_set(B.value, b)
        mpz_set(C.value, c)
        mpz_set(D.value, d)
        mpz_set(E.value, e)
        mpz_set(P.value, pp)
        f = ntl.ZZ_pX([E, D, C, B, A], P)
        f /= ntl.ZZ_pX([A], P)  # now f is monic, and we are done with A,B,C,D,E
        mpz_set(qq, A.value)  # qq is the leading coefficient of the polynomial
        f_factzn = f.factor()
        result = 0
        for factor, exponent in f_factzn:
            if exponent&1:
                result = 1
                break
        if result == 0 and mpz_legendre(qq, pp) == 1:
            result = 1
        if result:
            return 1

        f = ntl.ZZ_pX([1], P)
        for factor, exponent in f_factzn:
            for j in range(exponent // 2):
                f *= factor

        f /= f.leading_coefficient()
        f_factzn = f.factor()

        has_roots = 0
        j = 0
        for factor, exponent in f_factzn:
            if factor.degree() == 1:
                has_roots = 1
                A = P - Integer(factor[0])
                mpz_set(roots[j], A.value)
                j += 1
        if not has_roots:
            return 0

        i = f.degree()
        mpz_init(aaa)
        mpz_init(bbb)
        mpz_init(ccc)
        mpz_init(ddd)
        mpz_init(eee)

        if i == 0:  # g == 1
            mpz_set(aaa, a)
            mpz_set(bbb, b)
            mpz_set(ccc, c)
            mpz_set(ddd, d)
            mpz_sub(eee, e, qq)
        elif i == 1:  # g == x + rr
            mpz_set(aaa, a)
            mpz_set(bbb, b)
            mpz_sub(ccc, c, qq)
            A = Integer(f[0])
            mpz_set(rr, A.value)
            mpz_mul(ss, rr, qq)
            mpz_set(ddd,d)
            mpz_sub(ddd, ddd, ss)
            mpz_sub(ddd, ddd, ss)
            mpz_set(eee,e)
            mpz_mul(ss, ss, rr)
            mpz_sub(eee, eee, ss)
            mpz_divexact(ss, ss, rr)
        elif i == 2:  # g == x^2 + rr*x + ss
            mpz_sub(aaa, a, qq)
            A = Integer(f[1])
            mpz_set(rr, A.value)
            mpz_init(tt)
            mpz_mul(tt, rr, qq)
            mpz_set(bbb,b)
            mpz_submul_ui(bbb, tt, 2)
            mpz_set(ccc,c)
            mpz_submul(ccc, tt, rr)
            A = Integer(f[0])
            mpz_set(ss, A.value)
            mpz_mul(tt, ss, qq)
            mpz_set(eee,e)
            mpz_submul(eee, tt, ss)
            mpz_mul_ui(tt, tt, 2)
            mpz_sub(ccc, ccc, tt)
            mpz_set(ddd,d)
            mpz_submul(ddd, tt, rr)
            mpz_clear(tt)
        mpz_divexact(aaa, aaa, pp)
        mpz_divexact(bbb, bbb, pp)
        mpz_divexact(ccc, ccc, pp)
        mpz_divexact(ddd, ddd, pp)
        mpz_divexact(eee, eee, pp)
        # now aaa, bbb, ccc, ddd, eee represents h(x)

        result = 0
        mpz_init(tt)
        for i in range(j):
            mpz_mul(tt, aaa, roots[i])
            mpz_add(tt, tt, bbb)
            mpz_mul(tt, tt, roots[i])
            mpz_add(tt, tt, ccc)
            mpz_mul(tt, tt, roots[i])
            mpz_add(tt, tt, ddd)
            mpz_mul(tt, tt, roots[i])
            mpz_add(tt, tt, eee)
            # tt == h(r) mod p
            mpz_mod(tt, tt, pp)
            if mpz_sgn(tt) == 0:
                fmpz_poly_zero(f1)
                fmpz_poly_zero(linear)
                fmpz_poly_set_coeff_mpz(f1, 0, e)
                fmpz_poly_set_coeff_mpz(f1, 1, d)
                fmpz_poly_set_coeff_mpz(f1, 2, c)
                fmpz_poly_set_coeff_mpz(f1, 3, b)
                fmpz_poly_set_coeff_mpz(f1, 4, a)
                fmpz_poly_set_coeff_mpz(linear, 0, roots[i])
                fmpz_poly_set_coeff_mpz(linear, 1, pp)
                fmpz_poly_compose(f1, f1, linear)
                fmpz_poly_scalar_fdiv_mpz(f1, f1, pp)
                fmpz_poly_scalar_fdiv_mpz(f1, f1, pp)

                mpz_init(aa)
                mpz_init(bb)
                mpz_init(cc)
                mpz_init(dd)
                mpz_init(ee)
                fmpz_poly_get_coeff_mpz(aa, f1, 4)
                fmpz_poly_get_coeff_mpz(bb, f1, 3)
                fmpz_poly_get_coeff_mpz(cc, f1, 2)
                fmpz_poly_get_coeff_mpz(dd, f1, 1)
                fmpz_poly_get_coeff_mpz(ee, f1, 0)
                result = Zp_soluble_siksek_large_p(aa, bb, cc, dd, ee, pp, f1, linear)
                mpz_clear(aa)
                mpz_clear(bb)
                mpz_clear(cc)
                mpz_clear(dd)
                mpz_clear(ee)
                if result == 1:
                    break
        mpz_clear(aaa)
        mpz_clear(bbb)
        mpz_clear(ccc)
        mpz_clear(ddd)
        mpz_clear(eee)
        mpz_clear(tt)
        return result
    else:
        # Step II in Alg. 5.3.1 of Siksek's thesis
        A = Integer(0)
        B = Integer(0)
        C = Integer(0)
        D = Integer(0)
        E = Integer(0)
        P = Integer(0)
        mpz_set(A.value, a)
        mpz_set(B.value, b)
        mpz_set(C.value, c)
        mpz_set(D.value, d)
        mpz_set(E.value, e)
        mpz_set(P.value, pp)
        f = ntl.ZZ_pX([E, D, C, B, A], P)
        f /= ntl.ZZ_pX([A], P)  # now f is monic
        f_factzn = f.factor()

        has_roots = 0
        has_single_roots = 0
        j = 0
        for factor, exponent in f_factzn:
            if factor.degree() == 1:
                has_roots = 1
                if exponent == 1:
                    has_single_roots = 1
                    break
                A = P - Integer(factor[0])
                mpz_set(roots[j], A.value)
                j += 1

        if not has_roots:
            return 0
        if has_single_roots:
            return 1

        result = 0
        if j > 0:
            mpz_init(aa)
            mpz_init(bb)
            mpz_init(cc)
            mpz_init(dd)
            mpz_init(ee)
        for i in range(j):
            fmpz_poly_zero(f1)
            fmpz_poly_zero(linear)
            fmpz_poly_set_coeff_mpz(f1, 0, e)
            fmpz_poly_set_coeff_mpz(f1, 1, d)
            fmpz_poly_set_coeff_mpz(f1, 2, c)
            fmpz_poly_set_coeff_mpz(f1, 3, b)
            fmpz_poly_set_coeff_mpz(f1, 4, a)
            fmpz_poly_set_coeff_mpz(linear, 0, roots[i])
            fmpz_poly_set_coeff_mpz(linear, 1, pp)
            fmpz_poly_compose(f1, f1, linear)
            fmpz_poly_scalar_fdiv_mpz(f1, f1, pp)
            fmpz_poly_get_coeff_mpz(aa, f1, 4)
            fmpz_poly_get_coeff_mpz(bb, f1, 3)
            fmpz_poly_get_coeff_mpz(cc, f1, 2)
            fmpz_poly_get_coeff_mpz(dd, f1, 1)
            fmpz_poly_get_coeff_mpz(ee, f1, 0)
            result = Zp_soluble_siksek_large_p(aa, bb, cc, dd, ee, pp, f1, linear)
            if result == 1:
                break
        if j > 0:
            mpz_clear(aa)
            mpz_clear(bb)
            mpz_clear(cc)
            mpz_clear(dd)
            mpz_clear(ee)
        return result

cdef bint Qp_soluble_siksek(mpz_t A, mpz_t B, mpz_t C, mpz_t D, mpz_t E,
                            mpz_t p, unsigned long P,
                            nmod_poly_factor_t f_factzn, fmpz_poly_t f1,
                            fmpz_poly_t linear) noexcept:
    """
    Uses Samir Siksek's thesis results to determine whether the quartic is
    locally soluble at p.
    """
    cdef int result = 0
    cdef mpz_t a, b, c, d, e
    cdef nmod_poly_t f
    nmod_poly_init(f, P)

    mpz_init_set(a, A)
    mpz_init_set(b, B)
    mpz_init_set(c, C)
    mpz_init_set(d, D)
    mpz_init_set(e, E)

    if Zp_soluble_siksek(a, b, c, d, e,
                         p, P, f_factzn, f, f1, linear):
        result = 1
    else:
        mpz_set(a, A)
        mpz_set(b, B)
        mpz_set(c, C)
        mpz_set(d, D)
        mpz_set(e, E)
        if Zp_soluble_siksek(e, d, c, b, a,
                             p, P,f_factzn, f, f1, linear):
            result = 1

    mpz_clear(a)
    mpz_clear(b)
    mpz_clear(c)
    mpz_clear(d)
    mpz_clear(e)
    nmod_poly_clear(f)
    return result


cdef bint Qp_soluble_siksek_large_p(mpz_t A, mpz_t B, mpz_t C, mpz_t D, mpz_t E,
                                    mpz_t p, fmpz_poly_t f1, fmpz_poly_t linear) noexcept:
    """
    Uses Samir Siksek's thesis results to determine whether the quartic is
    locally soluble at p, when p is bigger than wordsize, and we can't use
    FLINT.
    """
    cdef int result = 0
    cdef mpz_t a, b, c, d, e

    mpz_init_set(a,A)
    mpz_init_set(b,B)
    mpz_init_set(c,C)
    mpz_init_set(d,D)
    mpz_init_set(e,E)

    if Zp_soluble_siksek_large_p(a, b, c, d, e,
                                 p, f1, linear):
        result = 1
    else:
        mpz_set(a,A)
        mpz_set(b,B)
        mpz_set(c,C)
        mpz_set(d,D)
        mpz_set(e,E)
        if Zp_soluble_siksek_large_p(e,d,c,b,a,p,f1,linear):
            result = 1

    mpz_clear(a)
    mpz_clear(b)
    mpz_clear(c)
    mpz_clear(d)
    mpz_clear(e)
    return result

cdef bint Qp_soluble_BSD(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e, mpz_t p) noexcept:
    """
    Uses the original test of Birch and Swinnerton-Dyer to test for local
    solubility of the quartic at p.
    """
    cdef mpz_t zero
    cdef int result = 0
    mpz_init_set_ui(zero, 0)
    if Zp_soluble_BSD(a, b, c, d, e,
                      zero, p, 0):
        result = 1
    elif Zp_soluble_BSD(e, d, c, b, a,
                        zero, p, 1):
        result = 1
    mpz_clear(zero)
    return result

cdef bint Qp_soluble(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e, mpz_t p) noexcept:
    """
    Try the BSD approach for a few residue classes and if no solution is found,
    switch to Siksek to try to prove insolubility.
    """
    cdef int bsd_sol, sik_sol
    cdef unsigned long pp
    cdef fmpz_poly_t f1, linear
    cdef nmod_poly_factor_t f_factzn
    bsd_sol = Qp_soluble_BSD(a, b, c, d, e, p)
    if mpz_cmp_ui(p,N_RES_CLASSES_BSD)>0 and not bsd_sol:
        fmpz_poly_init(f1)
        fmpz_poly_init(linear)
        if mpz_fits_ulong_p(p):
            nmod_poly_factor_init(f_factzn)
            pp = mpz_get_ui(p)
            sik_sol = Qp_soluble_siksek(a, b, c, d, e,
                                        p, pp, f_factzn, f1,linear)
            nmod_poly_factor_clear(f_factzn)
        else:
            sik_sol = Qp_soluble_siksek_large_p(a, b, c, d, e,
                                                p, f1, linear)
        fmpz_poly_clear(f1)
        fmpz_poly_clear(linear)
    else:
        sik_sol = bsd_sol
    return sik_sol


def test_qpls(a, b, c, d, e, p):
    """
    Testing function for Qp_soluble.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.descent_two_isogeny import test_qpls as tq
        sage: tq(1,2,3,4,5,7)
        1
    """
    cdef Integer A, B, C, D, E, P
    cdef int result
    cdef mpz_t aa, bb, cc, dd, ee, pp
    A = Integer(a)
    B = Integer(b)
    C = Integer(c)
    D = Integer(d)
    E = Integer(e)
    P = Integer(p)
    mpz_init_set(aa, A.value)
    mpz_init_set(bb, B.value)
    mpz_init_set(cc, C.value)
    mpz_init_set(dd, D.value)
    mpz_init_set(ee, E.value)
    mpz_init_set(pp, P.value)
    result = Qp_soluble(aa, bb, cc, dd, ee, pp)
    mpz_clear(aa)
    mpz_clear(bb)
    mpz_clear(cc)
    mpz_clear(dd)
    mpz_clear(ee)
    mpz_clear(pp)
    return result


cdef int everywhere_locally_soluble(mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t e) except -1:
    """
    Return whether the quartic has local solutions at all primes `p`.
    """
    cdef Integer A, B, C, D, E, Delta,p
    cdef mpz_t mpz_2
    A = Integer(0)
    B = Integer(0)
    C = Integer(0)
    D = Integer(0)
    E = Integer(0)
    mpz_set(A.value, a)
    mpz_set(B.value, b)
    mpz_set(C.value, c)
    mpz_set(D.value, d)
    mpz_set(E.value, e)
    f = (((A*x_ZZ + B)*x_ZZ + C)*x_ZZ + D)*x_ZZ + E

    # RR soluble:
    if mpz_sgn(a) != 1:
        if not real_roots(f):
            return 0

    # Q2 soluble:
    mpz_init_set_ui(mpz_2, 2)
    if not Qp_soluble(a, b, c, d, e, mpz_2):
        mpz_clear(mpz_2)
        return 0
    mpz_clear(mpz_2)

    # Odd finite primes
    Delta = f.discriminant()
    for p in prime_divisors(Delta):
        if p == 2:
            continue
        if not Qp_soluble(a, b, c, d, e, p.value):
            return 0

    return 1


def test_els(a, b, c, d, e):
    """
    Doctest function for cdef int everywhere_locally_soluble(mpz_t, mpz_t, mpz_t, mpz_t, mpz_t).

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.descent_two_isogeny import test_els
        sage: for _ in range(1000):
        ....:     a,b,c,d,e = randint(1,1000), randint(1,1000), randint(1,1000), randint(1,1000), randint(1,1000)
        ....:     if pari.Pol([a,b,c,d,e]).hyperellratpoints(1000, 1):
        ....:         try:
        ....:             if not test_els(a,b,c,d,e):
        ....:                 print("This never happened", a, b, c, d, e)
        ....:         except ValueError:
        ....:             continue
    """
    cdef Integer A, B, C, D, E
    A = Integer(a)
    B = Integer(b)
    C = Integer(c)
    D = Integer(d)
    E = Integer(e)
    return everywhere_locally_soluble(A.value, B.value, C.value, D.value, E.value)


cdef int count(mpz_t c_mpz, mpz_t d_mpz, mpz_t *p_list, unsigned long p_list_len,
               int global_limit_small, int global_limit_large,
               int verbosity, bint selmer_only, mpz_t n1, mpz_t n2) except -1:
    """
    Count the number of els/gls quartic 2-covers of E.
    """
    cdef unsigned long n_primes, i
    cdef bint found_global_points, els, check_negs, verbose = (verbosity > 4)
    cdef Integer a_Int, c_Int, e_Int
    cdef mpz_t c_sq_mpz, d_prime_mpz
    cdef mpz_t n_divisors, j

    mpz_init(c_sq_mpz)
    mpz_mul(c_sq_mpz, c_mpz, c_mpz)
    mpz_init_set(d_prime_mpz, c_sq_mpz)
    mpz_submul_ui(d_prime_mpz, d_mpz, 4)
    check_negs = 0
    if mpz_sgn(d_prime_mpz) > 0:
        if mpz_sgn(c_mpz) >= 0 or mpz_cmp(c_sq_mpz, d_prime_mpz) <= 0:
            check_negs = 1
    mpz_clear(c_sq_mpz)
    mpz_clear(d_prime_mpz)

    # Set up coefficient array, and static variables
    cdef mpz_t *coeffs = <mpz_t *> sig_malloc(5 * sizeof(mpz_t))
    for i in range(5):
        mpz_init(coeffs[i])
    mpz_set_ui(coeffs[1], 0)     #
    mpz_set(coeffs[2], c_mpz)    # These never change
    mpz_set_ui(coeffs[3], 0)     #

    # Get prime divisors, and put them in an mpz_t array
    # (this block, by setting check_negs, takes care of
    # local solubility over RR)
    cdef mpz_t *p_div_d_mpz = <mpz_t *> sig_malloc((p_list_len+1) * sizeof(mpz_t))
    n_primes = 0
    for i in range(p_list_len):
        if mpz_divisible_p(d_mpz, p_list[i]):
            mpz_init(p_div_d_mpz[n_primes])
            mpz_set(p_div_d_mpz[n_primes], p_list[i])
            n_primes += 1
    if check_negs:
        mpz_init(p_div_d_mpz[n_primes])
        mpz_set_si(p_div_d_mpz[n_primes], -1)
        n_primes += 1
    mpz_init_set_ui(n_divisors, 1)
    mpz_mul_2exp(n_divisors, n_divisors, n_primes)

    mpz_init_set_ui(j, 0)
    if not selmer_only:
        mpz_set_ui(n1, 0)
    mpz_set_ui(n2, 0)
    while mpz_cmp(j, n_divisors) < 0:
        mpz_set_ui(coeffs[4], 1)
        for i in range(n_primes):
            if mpz_tstbit(j, i):
                mpz_mul(coeffs[4], coeffs[4], p_div_d_mpz[i])
        if verbosity > 3:
            a_Int = Integer(0)
            mpz_set(a_Int.value, coeffs[4])
            print('\nSquarefree divisor:', a_Int)
        mpz_divexact(coeffs[0], d_mpz, coeffs[4])
        found_global_points = 0
        if not selmer_only:
            if verbose:
                print("\nCalling ratpoints for small point search")
            found_global_points = ratpoints_mpz_exists_only(coeffs, 4, global_limit_small)
            if found_global_points:
                if verbosity > 2:
                    a_Int = Integer(0)
                    mpz_set(a_Int.value, coeffs[4])
                    c_Int = Integer(0)
                    mpz_set(c_Int.value, coeffs[2])
                    e_Int = Integer(0)
                    mpz_set(e_Int.value, coeffs[0])
                    print('Found small global point, quartic (%d,%d,%d,%d,%d)' % (a_Int, 0, c_Int, 0, e_Int))
                mpz_add_ui(n1, n1, 1)
                mpz_add_ui(n2, n2, 1)
            if verbose:
                print("\nDone calling ratpoints for small point search")
        if not found_global_points:
            # Test whether the quartic is everywhere locally soluble:
            els = 1
            for i in range(p_list_len):
                if not Qp_soluble(coeffs[4], coeffs[3], coeffs[2], coeffs[1], coeffs[0], p_list[i]):
                    els = 0
                    break
            if els:
                if verbosity > 2:
                    a_Int = Integer(0)
                    mpz_set(a_Int.value, coeffs[4])
                    c_Int = Integer(0)
                    mpz_set(c_Int.value, coeffs[2])
                    e_Int = Integer(0)
                    mpz_set(e_Int.value, coeffs[0])
                    print('ELS without small global points, quartic (%d,%d,%d,%d,%d)' % (a_Int, 0, c_Int, 0, e_Int))
                mpz_add_ui(n2, n2, 1)
                if not selmer_only:
                    if verbose:
                        print("\nCalling ratpoints for large point search")
                    found_global_points = ratpoints_mpz_exists_only(coeffs, 4, global_limit_large)
                    if found_global_points:
                        if verbosity > 2:
                            print('  -- Found large global point.')
                        mpz_add_ui(n1, n1, 1)
                    if verbose:
                        print("\nDone calling ratpoints for large point search")
        mpz_add_ui(j, j, 1)
    mpz_clear(j)
    for i in range(n_primes):
        mpz_clear(p_div_d_mpz[i])
    sig_free(p_div_d_mpz)
    mpz_clear(n_divisors)
    for i in range(5):
        mpz_clear(coeffs[i])
    sig_free(coeffs)
    return 0


def two_descent_by_two_isogeny(E,
                               int global_limit_small=10,
                               int global_limit_large=10000,
                               int verbosity=0,
                               bint selmer_only=0, bint proof=1):
    """
    Given an elliptic curve E with a two-isogeny phi : E --> E' and dual isogeny
    phi', runs a two-isogeny descent on E, returning n1, n2, n1' and n2'. Here
    n1 is the number of quartic covers found with a rational point, and n2 is
    the number which are ELS.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
        sage: E = EllipticCurve('14a')
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny(E)
        sage: log(n1,2) + log(n1_prime,2) - 2 # the rank
        0
        sage: E = EllipticCurve('65a')
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny(E)
        sage: log(n1,2) + log(n1_prime,2) - 2 # the rank
        1

        sage: # needs sage.symbolic
        sage: x,y = var('x,y')
        sage: E = EllipticCurve(y^2 == x^3 + x^2 - 25*x + 39)
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny(E)
        sage: log(n1,2) + log(n1_prime,2) - 2  # the rank
        2
        sage: E = EllipticCurve(y^2 + x*y + y == x^3 - 131*x + 558)
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny(E)
        sage: log(n1,2) + log(n1_prime,2) - 2  # the rank
        3

    Using the verbosity option::

        sage: E = EllipticCurve('14a')
        sage: two_descent_by_two_isogeny(E, verbosity=1)
        2-isogeny
        Results:
        2 <= #E(Q)/phi'(E'(Q)) <= 2
        2 <= #E'(Q)/phi(E(Q)) <= 2
        #Sel^(phi')(E'/Q) = 2
        #Sel^(phi)(E/Q) = 2
        1 <= #Sha(E'/Q)[phi'] <= 1
        1 <= #Sha(E/Q)[phi] <= 1
        1 <= #Sha(E/Q)[2], #Sha(E'/Q)[2] <= 1
        0 <= rank of E(Q) = rank of E'(Q) <= 0
        (2, 2, 2, 2)

    Handling curves whose discriminants involve larger than wordsize primes::

        sage: E = EllipticCurve('14a')
        sage: E = E.quadratic_twist(next_prime(10^20))
        sage: E
        Elliptic Curve defined by y^2 = x^3 + x^2 + 716666666666666667225666666666666666775672*x - 391925925925925926384240370370370370549019837037037037060249356 over Rational Field
        sage: E.discriminant().factor()
        -1 * 2^18 * 7^3 * 100000000000000000039^6
        sage: log(100000000000000000039.0, 2.0)
        66.438...
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny(E)
        sage: log(n1,2) + log(n1_prime,2) - 2 # the rank
        0

    TESTS:

    Here we contrive an example to demonstrate that a keyboard interrupt
    is caught. Here we let `E` be the smallest optimal curve with two-torsion
    and nontrivial `Sha[2]`. This ensures that the two-descent will be looking
    for rational points which do not exist, and by setting global_limit_large
    to a very high bound, it will still be working when we simulate a ``CTRL-C``::

        sage: from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
        sage: E = EllipticCurve('960d'); E
        Elliptic Curve defined by y^2 = x^3 - x^2 - 900*x - 10098 over Rational Field
        sage: E.sha().an()
        4
        sage: from sage.doctest.util import ensure_interruptible_after
        sage: with ensure_interruptible_after(0.5): two_descent_by_two_isogeny(E, global_limit_large=10^8)
    """
    cdef Integer a1, a2, a3, a4, a6, s2, s4, s6
    cdef Integer c, d, x0
    cdef list x_list
    assert E.torsion_order() % 2 == 0, 'Need rational two-torsion for isogeny descent.'
    if verbosity > 0:
        print('\n2-isogeny')
        if verbosity > 1:
            print('\nchanging coordinates')
    a1 = Integer(E.a1())
    a2 = Integer(E.a2())
    a3 = Integer(E.a3())
    a4 = Integer(E.a4())
    a6 = Integer(E.a6())
    if a1 == 0 == a3:
        s2 = a2
        s4 = a4
        s6 = a6
    else:
        s2 = a1*a1+4*a2
        s4 = 8*(a1*a3+2*a4)
        s6 = 16*(a3*a3+4*a6)
    f = ((x_ZZ + s2)*x_ZZ + s4)*x_ZZ + s6
    x_list = f.roots()  # over ZZ -- use FLINT directly?
    x0 = x_list[0][0]
    c = 3*x0+s2
    d = (c+s2)*x0+s4
    return two_descent_by_two_isogeny_work(c, d,
                                           global_limit_small,
                                           global_limit_large, verbosity,
                                           selmer_only, proof)


def two_descent_by_two_isogeny_work(Integer c, Integer d,
                                    int global_limit_small=10,
                                    int global_limit_large=10000,
                                    int verbosity=0, bint selmer_only=0,
                                    bint proof=1):
    """
    Do all the work in doing a two-isogeny descent.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny_work
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny_work(13,128)
        sage: log(n1,2) + log(n1_prime,2) - 2 # the rank
        0
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny_work(1,-16)
        sage: log(n1,2) + log(n1_prime,2) - 2 # the rank
        1
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny_work(10,8)
        sage: log(n1,2) + log(n1_prime,2) - 2 # the rank
        2
        sage: n1, n2, n1_prime, n2_prime = two_descent_by_two_isogeny_work(85,320)
        sage: log(n1,2) + log(n1_prime,2) - 2 # the rank
        3
    """
    cdef mpz_t c_mpz, d_mpz, c_prime_mpz, d_prime_mpz
    cdef mpz_t *p_list_mpz
    cdef unsigned long i, j, p, p_list_len
    cdef Integer P, n1, n2, n1_prime, n2_prime, c_prime, d_prime
    cdef object PO
    cdef bint found, d_neg, d_prime_neg
    cdef n_factor_t fact
    cdef list primes
    mpz_init_set(c_mpz, c.value)
    mpz_init_set(d_mpz, d.value)
    mpz_init(c_prime_mpz)
    mpz_init(d_prime_mpz)
    mpz_mul_si(c_prime_mpz, c_mpz, -2)
    mpz_mul(d_prime_mpz, c_mpz, c_mpz)
    mpz_submul_ui(d_prime_mpz, d_mpz, 4)

    d_neg = 0
    d_prime_neg = 0
    if mpz_sgn(d_mpz) < 0:
        d_neg = 1
        mpz_neg(d_mpz, d_mpz)
    if mpz_sgn(d_prime_mpz) < 0:
        d_prime_neg = 1
        mpz_neg(d_prime_mpz, d_prime_mpz)
    if mpz_fits_ulong_p(d_mpz) and mpz_fits_ulong_p(d_prime_mpz):
        # Factor very quickly using FLINT.
        p_list_mpz = <mpz_t *> sig_malloc(20 * sizeof(mpz_t))
        mpz_init_set_ui(p_list_mpz[0], 2)
        p_list_len = 1
        n_factor_init(&fact)
        n_factor(&fact, mpz_get_ui(d_mpz), proof)
        for i in range(fact.num):
            p = fact.p[i]
            if p != 2:
                mpz_init_set_ui(p_list_mpz[p_list_len], p)
                p_list_len += 1
        n_factor(&fact, mpz_get_ui(d_prime_mpz), proof)
        for i in range(fact.num):
            p = fact.p[i]
            found = 0
            for j in range(p_list_len):
                if mpz_cmp_ui(p_list_mpz[j], p) == 0:
                    found = 1
                    break
            if not found:
                mpz_init_set_ui(p_list_mpz[p_list_len], p)
                p_list_len += 1
    else:
        # Factor more slowly using Pari via Python.
        from sage.libs.pari import pari
        d = Integer(0)
        mpz_set(d.value, d_mpz)
        primes = list(pari(d).factor()[0])
        d_prime = Integer(0)
        mpz_set(d_prime.value, d_prime_mpz)
        for PO in pari(d_prime).factor()[0]:
            if PO not in primes:
                primes.append(PO)
        P = Integer(2)
        if P not in primes:
            primes.append(P)
        p_list_len = len(primes)
        p_list_mpz = <mpz_t *> sig_malloc(p_list_len * sizeof(mpz_t))
        for i in range(p_list_len):
            P = Integer(primes[i])
            mpz_init_set(p_list_mpz[i], P.value)
    if d_neg:
        mpz_neg(d_mpz, d_mpz)
    if d_prime_neg:
        mpz_neg(d_prime_mpz, d_prime_mpz)

    if verbosity > 1:
        c_prime = -2*c
        d_prime = c*c-4*d
        print('\nnew curve is y^2 == x( x^2 + (%d)x + (%d) )' % (int(c), int(d)))
        print('new isogenous curve is' +
              ' y^2 == x( x^2 + (%d)x + (%d) )' % (int(c_prime), int(d_prime)))

    n1 = Integer(0)
    n2 = Integer(0)
    n1_prime = Integer(0)
    n2_prime = Integer(0)
    count(c.value, d.value, p_list_mpz, p_list_len,
          global_limit_small, global_limit_large, verbosity, selmer_only,
          n1.value, n2.value)
    count(c_prime_mpz, d_prime_mpz, p_list_mpz, p_list_len,
          global_limit_small, global_limit_large, verbosity, selmer_only,
          n1_prime.value, n2_prime.value)

    for i in range(p_list_len):
        mpz_clear(p_list_mpz[i])
    sig_free(p_list_mpz)

    if verbosity > 0:
        print("\nResults:")
        print(n1, "<= #E(Q)/phi'(E'(Q)) <=", n2)
        print(n1_prime, "<= #E'(Q)/phi(E(Q)) <=", n2_prime)
        print("#Sel^(phi')(E'/Q) =", n2)
        print("#Sel^(phi)(E/Q) =", n2_prime)
        print("1 <= #Sha(E'/Q)[phi'] <=", n2/n1)
        print("1 <= #Sha(E/Q)[phi] <=", n2_prime/n1_prime)
        print("1 <= #Sha(E/Q)[2], #Sha(E'/Q)[2] <=", (n2_prime/n1_prime)*(n2/n1))
        a = Integer(n1*n1_prime).log(Integer(2))
        e = Integer(n2*n2_prime).log(Integer(2))
        print(a - 2, "<= rank of E(Q) = rank of E'(Q) <=", e - 2)

    return n1, n2, n1_prime, n2_prime


cdef bint ratpoints_mpz_exists_only(mpz_t *coeffs, long degree, long H) except -1:
    """
    Search for projective points on the hyperelliptic curve
    ``y^2 = P(x)``.

    INPUT:

    - ``coeffs`` -- an array of length ``degree + 1`` giving the
      coefficients of ``P``, starting with the constant coefficient

    - ``degree`` -- degree of ``P``

    - ``H`` -- bound on the naive height for search

    OUTPUT: boolean, whether or not a projective point was found
    """
    sig_on()
    cdef GEN pol = cgetg(degree + 3, t_POL)
    pol[1] = evalvarn(0) + evalsigne(1)
    cdef long i
    for i in range(degree + 1):
        set_gel(pol, i+2, _new_GEN_from_mpz_t(coeffs[i]))

    # PARI checks only for affine points, so we manually check for
    # points at infinity (of the smooth model)
    cdef int r
    if degree % 2 or Z_issquare(gel(pol, degree+2)):
        r = 1
    else:
        R = hyperellratpoints(pol, stoi(H), 1)
        r = (lg(R) > 1)
    clear_stack()
    return r
