r"""
Fast calculation of cyclotomic polynomials

This module provides a function :func:`cyclotomic_coeffs`, which calculates the
coefficients of cyclotomic polynomials. This is not intended to be invoked
directly by the user, but it is called by the method
:meth:`~sage.rings.polynomial.polynomial_ring.PolynomialRing_generic.cyclotomic_polynomial`
method of univariate polynomial ring objects and the top-level
:func:`~sage.misc.functional.cyclotomic_polynomial` function.
"""

# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys

from libc.stdint cimport int64_t
from cysignals.memory cimport sig_malloc, check_calloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.arith.misc import factor, gcd
from sage.combinat.subset import subsets
from sage.misc.misc_c import prod
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.integer_ring import ZZ
from sage.structure.element cimport parent
from sage.structure.factorization import Factorization

cdef extern from *:
    """
    #include <stdint.h>
    #if defined(__SIZEOF_INT128__)
    #define DICKSON_WORD_BITS 62
    static inline int64_t mulmod_i64(int64_t a, int64_t b, int64_t m) {
        return (int64_t)(((unsigned __int128)a * (unsigned __int128)b) % (unsigned __int128)m);
    }
    #else
    #define DICKSON_WORD_BITS 31
    static inline int64_t mulmod_i64(int64_t a, int64_t b, int64_t m) {
        return (int64_t)(((uint64_t)a * (uint64_t)b) % (uint64_t)m);
    }
    #endif
    static inline int64_t invmod_i64(int64_t a, int64_t m) {
        int64_t t = 0, newt = 1;
        int64_t r = m, newr = a % m;
        if (newr < 0) newr += m;
        while (newr != 0) {
            int64_t q = r / newr;
            int64_t tmp = t - q * newt;
            t = newt; newt = tmp;
            tmp = r - q * newr;
            r = newr; newr = tmp;
        }
        if (r > 1) return 0;
        return t < 0 ? t + m : t;
    }
    static inline long gcd_long(long a, long b) {
        long tmp;
        while (b != 0) {
            tmp = b;
            b = a % b;
            a = tmp;
        }
        return a < 0 ? -a : a;
    }
    """
    # Largest bit size of moduli handled by the word-sized C code path
    int DICKSON_WORD_BITS
    # All arguments are expected to be nonnegative and reduced modulo ``m``
    int64_t mulmod_i64(int64_t a, int64_t b, int64_t m) noexcept nogil
    int64_t invmod_i64(int64_t a, int64_t m) noexcept nogil
    long gcd_long(long a, long b) noexcept nogil

try:
    from sage.libs.pari import pari
except ImportError:
    pass


def cyclotomic_coeffs(nn, sparse=None):
    """
    Return the coefficients of the `n`-th cyclotomic polynomial
    by using the formula

    .. MATH::

        \\Phi_n(x) = \\prod_{d|n} (1-x^{n/d})^{\\mu(d)}

    where `\\mu(d)` is the Möbius function that is 1 if `d` has an even
    number of distinct prime divisors, `-1` if it has an odd number of
    distinct prime divisors, and `0` if `d` is not squarefree.

    Multiplications and divisions by polynomials of the
    form `1-x^n` can be done very quickly in a single pass.

    If ``sparse`` is ``True``, the result is returned as a dictionary of
    the nonzero entries, otherwise the result is returned as a list
    of python ints.

    EXAMPLES::

        sage: from sage.rings.polynomial.cyclotomic import cyclotomic_coeffs
        sage: cyclotomic_coeffs(30)
        [1, 1, 0, -1, -1, -1, 0, 1, 1]
        sage: cyclotomic_coeffs(10^5)
        {0: 1, 10000: -1, 20000: 1, 30000: -1, 40000: 1}
        sage: R = QQ['x']
        sage: R(cyclotomic_coeffs(30))
        x^8 + x^7 - x^5 - x^4 - x^3 + x + 1

    Check that it has the right degree::

        sage: euler_phi(30)                                                             # needs sage.libs.pari
        8
        sage: R(cyclotomic_coeffs(14)).factor()                                         # needs sage.libs.pari
        x^6 - x^5 + x^4 - x^3 + x^2 - x + 1

    The coefficients are not always +/-1::

        sage: cyclotomic_coeffs(105)
        [1, 1, 1, 0, 0, -1, -1, -2, -1, -1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, -1,
         0, -1, 0, -1, 0, -1, 0, -1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -2,
         -1, -1, 0, 0, 1, 1, 1]

    In fact the height is not bounded by any polynomial in `n` (Erdős),
    although takes a while just to exceed linear::

        sage: v = cyclotomic_coeffs(1181895)
        sage: max(v)
        14102773

    The polynomial is a palindrome for any n::

        sage: n = ZZ.random_element(50000)
        sage: v = cyclotomic_coeffs(n, sparse=False)
        sage: v == list(reversed(v))
        True

    AUTHORS:

    - Robert Bradshaw (2007-10-27): initial version (inspired by work of Andrew
      Arnold and Michael Monagan)

    REFERENCE:

    - http://www.cecm.sfu.ca/~ada26/cyclotomic/
    """
    factors = factor(nn)
    if any(e != 1 for _, e in factors):
        # If there are primes that occur in the factorization with multiplicity
        # greater than one we use the fact that Phi_ar(x) = Phi_r(x^a) when all
        # primes dividing a divide r.
        rad = prod(p for p, _ in factors)
        rad_coeffs = cyclotomic_coeffs(rad, sparse=True)
        pow = int(nn // rad)
        if sparse is None or sparse:
            L = {}
        else:
            L = [0] * (1 + pow * prod(p - 1 for p, _ in factors))
        for mon, c in rad_coeffs.items():
            L[mon * pow] = c
        return L

    elif len(factors) == 1 and not sparse:
        # \Phi_p is easy to calculate for p prime.
        return [1] * factors[0][0]

    # The following bounds are from Michael Monagan:
    #    For all n < 169,828,113, the height of Phi_n(x) is less than 60 bits.
    #    At n = 169828113, we get a height of 31484567640915734951 which is 65 bits
    #    For n=10163195, the height of Phi_n(x) is 1376877780831,  40.32 bits.
    #    For n<10163195, the height of Phi_n(x) is <= 74989473, 26.16 bits.
    cdef long fits_long_limit = 169828113 if sizeof(long) >= 8 else 10163195
    if nn >= fits_long_limit and bateman_bound(nn) > sys.maxsize:
        # Do this to avoid overflow.
        print("Warning: using PARI (slow!)")
        return [int(a) for a in pari.polcyclo(nn).Vecrev()]

    cdef long d, max_deg = 0, n = nn
    primes = [int(p) for p, _ in factors]
    prime_subsets = list(subsets(primes))
    if n > 5000:
        prime_subsets.sort(key=lambda a: -prod(a))

    for s in prime_subsets:
        if len(s) % 2 == 0:
            d = prod(s)
            max_deg += n / d

    cdef long* coeffs = <long*>check_calloc(max_deg + 1, sizeof(long))
    coeffs[0] = 1

    cdef long k, dd, offset = 0, deg = 0
    for s in prime_subsets:
        if len(s) % 2 == 0:
            d = prod(s)
            dd = n / d
            # f *= (1-x^dd)
            sig_on()
            for k from deg + dd >= k >= dd:
                coeffs[k] -= coeffs[k - dd]
            deg += dd
            sig_off()

    prime_subsets.reverse()
    for s in prime_subsets:
        if len(s) % 2:
            d = prod(s)
            dd = n / d
            # f /= (1-x^dd)
            sig_on()
            for k from deg >= k > deg - dd:
                coeffs[k] = -coeffs[k]
            for k from deg - dd >= k >= offset:
                coeffs[k] = coeffs[k + dd] - coeffs[k]
            offset += dd
            sig_off()

    cdef long non_zero = 0
    if sparse is None:
        for k in range(offset, deg + 1):
            non_zero += coeffs[k] != 0
        sparse = (4 * non_zero) < (deg - offset)

    if sparse:
        L = {}
        for k in range(offset, deg + 1):
            if coeffs[k]:
                L[k - offset] = coeffs[k]
    else:
        L = [coeffs[k] for k in range(offset, deg + 1)]

    sig_free(coeffs)
    return L


def cyclotomic_value(n, x):
    r"""
    Return the value of the `n`-th cyclotomic polynomial evaluated at `x`.

    INPUT:

    - ``n`` -- an Integer, specifying which cyclotomic polynomial is to be
      evaluated

    - ``x`` -- an element of a ring

    OUTPUT: the value of the cyclotomic polynomial `\Phi_n` at `x`

    ALGORITHM:

    - Reduce to the case that `n` is squarefree: use the identity

    .. MATH::

        \Phi_n(x) = \Phi_q(x^{n/q})

    where `q` is the radical of `n`.

    - Use the identity

    .. MATH::

        \Phi_n(x) = \prod_{d | n} (x^d - 1)^{\mu(n / d)},

    where `\mu` is the Möbius function.

    - Handles the case that `x^d = 1` for some `d`, but not the case that
      `x^d - 1` is non-invertible: in this case polynomial evaluation is
      used instead.

    EXAMPLES::

        sage: cyclotomic_value(51, 3)
        1282860140677441
        sage: cyclotomic_polynomial(51)(3)
        1282860140677441

    It works for non-integral values as well::

        sage: cyclotomic_value(144, 4/3)
        79148745433504023621920372161/79766443076872509863361
        sage: cyclotomic_polynomial(144)(4/3)
        79148745433504023621920372161/79766443076872509863361

    TESTS::

        sage: elements = [-1, 0, 1, 2, 1/2, Mod(3, 8), Mod(3,11)]
        sage: R.<x> = QQ[]; elements += [x^2 + 2]
        sage: K.<i> = NumberField(x^2 + 1); elements += [i]                             # needs sage.rings.number_fields
        sage: elements += [GF(9,'a').gen()]                                             # needs sage.rings.finite_rings
        sage: elements += [Zp(3)(54)]                                                   # needs sage.rings.padics
        sage: for y in elements:
        ....:     for n in [1..60]:
        ....:         val1 = cyclotomic_value(n, y)
        ....:         val2 = cyclotomic_polynomial(n)(y)
        ....:         if val1 != val2:
        ....:             print("Wrong value for cyclotomic_value(%s, %s) in %s"%(n,y,parent(y)))
        ....:         if val1.parent() is not val2.parent():
        ....:             print("Wrong parent for cyclotomic_value(%s, %s) in %s"%(n,y,parent(y)))

        sage: cyclotomic_value(20, I)                                                   # needs sage.symbolic
        5
        sage: a = cyclotomic_value(10, mod(3, 11)); a
        6
        sage: a.parent()
        Ring of integers modulo 11
        sage: cyclotomic_value(30, -1.0)                                                # needs sage.rings.real_mpfr
        1.00000000000000

        sage: S.<t> = R.quotient(R.cyclotomic_polynomial(15))
        sage: cyclotomic_value(15, t)
        0
        sage: cyclotomic_value(30, t)
        2*t^7 - 2*t^5 - 2*t^3 + 2*t
        sage: S.<t> = R.quotient(x^10)
        sage: cyclotomic_value(2^128 - 1, t)
        -t^7 - t^6 - t^5 + t^2 + t + 1
        sage: cyclotomic_value(10, mod(3,4))
        1

    Check that the issue with symbolic element in :issue:`14982` is fixed::

        sage: a = cyclotomic_value(3, I)                                                # needs sage.rings.number_fields
        sage: parent(a)                                                                 # needs sage.rings.number_fields
        Number Field in I with defining polynomial x^2 + 1 with I = 1*I
    """
    n = ZZ(n)
    if n < 3:
        if n == 1:
            return x - ZZ.one()
        if n == 2:
            return x + ZZ.one()
        raise ValueError("n must be positive")

    P = parent(x)
    try:
        return P(pari.polcyclo(n, x).sage())
    except Exception:
        pass
    one = P(1)

    # The following is modeled on the implementation in PARI and is
    # used for cases for which PARI doesn't work. These are in
    # particular:
    # - n does not fit in a C long;
    # - x is some Sage type which cannot be converted to PARI;
    # - PARI's algorithm encounters a zero-divisor which is not zero.

    factors = n.factor()
    cdef Py_ssize_t i, j, ti, L, root_of_unity = -1
    primes = [p for p, _ in factors]
    L = len(primes)
    if any(e != 1 for _, e in factors):
        # If there are primes that occur in the factorization with multiplicity
        # greater than one we use the fact that Phi_ar(x) = Phi_r(x^a) when all
        # primes dividing a divide r.
        rad = prod(primes)
        pow = n // rad
        x = x ** pow
        n = rad
    if x == 1:
        # if n is prime, return n
        if L == 1:
            return n * x  # in case the parent of x has nonzero characteristic
        return x
    xd = [x]  # the x^d for d | n
    cdef char mu
    cdef char* md = <char*>sig_malloc(sizeof(char) * (1 << L))  # the mu(d) for d | n
    try:
        md[0] = 1
        if L & 1:
            mu = -1
            num = 1
            den = x - 1
        else:
            mu = 1
            num = x - 1
            den = 1
        for i in range(L):
            ti = 1 << i
            p = primes[i]
            for j in range(ti):
                xpow = xd[j]**p
                xd.append(xpow)
                md[ti + j] = -md[j]
                # if xpow = 1, we record such smallest index,
                # and deal with the corresponding factors at the end.
                if xpow == one:
                    if root_of_unity == -1:
                        root_of_unity = ti + j
                elif mu == md[ti + j]:
                    num *= xpow - one
                else:
                    den *= xpow - one
    finally:
        sig_free(md)
    try:
        ans = num / den
    except ZeroDivisionError:
        # We fall back on evaluation of the cyclotomic polynomial.
        # This case is triggered in cyclotomic_value(10, mod(3, 4)) for example.
        from sage.misc.functional import cyclotomic_polynomial
        return cyclotomic_polynomial(n)(x)
    if root_of_unity >= 0:
        # x is a root of unity.  If root_of_unity=2^L, x is a primitive
        # root of unity and the value is zero
        if root_of_unity == (1 << L) - 1:
            return x - x  # preserves the parent, as well as precision for p-adic x
        # x is a primitive d-th root of unity, where d|n and d<n.
        # If root_of_unity = (1<<L) - (1<<(i-1)) - 1 for some i < L,
        # then n/d == primes[i] and we need to multiply by primes[i],
        # otherwise n/d is composite and nothing more needs to be done.
        for i in range(L):
            if root_of_unity + (1 << i) + 1 == 1 << L:
                ans *= primes[i]
                break
    return x.parent()(ans)


def bateman_bound(nn):
    """
    Reference:

    Bateman, P. T.; Pomerance, C.; Vaughan, R. C.
    *On the size of the coefficients of the cyclotomic polynomial.*

    EXAMPLES::

        sage: from sage.rings.polynomial.cyclotomic import bateman_bound
        sage: bateman_bound(2**8 * 1234567893377)                                       # needs sage.libs.pari
        66944986927
    """
    _, n = nn.val_unit(2)
    primes = [p for p, _ in factor(n)]
    j = len(primes)
    return prod(primes[k] ** (2 ** (j - k - 2) - 1) for k in range(j - 2))


# ============================================================================
# Dickson Engine for factoring x^n - 1 and cyclotomic polynomials over Z/p^e Z
# ============================================================================

cdef list _poly_substitute_power(list poly, Py_ssize_t power):
    """
    Return the coefficient list of ``poly(x^power)``.
    """
    cdef Py_ssize_t deg = len(poly) - 1
    cdef Py_ssize_t i
    cdef list res = [0] * (deg * power + 1)
    for i in range(deg + 1):
        res[i * power] = poly[i]
    return res


cdef list _coeff_list(f):
    """
    Return the coefficients of the polynomial ``f`` as a list of Python ints.
    """
    cdef list L = [int(c) for c in f.list()]
    return L if L else [0]


cdef tuple _find_primitive_seed_fp(p_val, Py_ssize_t m_unram, Py_ssize_t m_order):
    r"""
    Return a pair ``(G, H_inv)`` of polynomials over `\GF{p}`.

    Here ``G`` is a monic irreducible factor of degree ``m_order`` of the
    cyclotomic polynomial `\Phi_m(x)`, where `m` is ``m_unram``, and ``H_inv``
    is `H(x)^{-1} \bmod G(x)` for the cofactor `H = (x^m - 1)/G`. The inverse
    is obtained without forming `H`, from the derivative identity
    `H(x)^{-1} \equiv m^{-1} x G'(x) \pmod{G(x)}`.
    """
    from sage.arith.misc import euler_phi, power_mod
    from sage.rings.finite_rings.finite_field_constructor import GF

    F_p = GF(p_val)
    R_fp = F_p['x']
    x = R_fp.gen()

    if m_order == euler_phi(m_unram):
        # Phi_m(x) is irreducible over GF(p)
        G = R_fp.cyclotomic_polynomial(m_unram)
    elif m_order == 1:
        G = x - F_p.zeta(m_unram)
    else:
        # Project a random element of GF(p^d) = GF(p)[x]/(P) onto the
        # subgroup of m-th roots of unity until a primitive one is found.
        p_int = int(p_val)
        P = R_fp.irreducible_element(m_order)
        exp = (p_int ** m_order - 1) // m_unram
        prime_divs = [int(q) for q, _ in factor(m_unram)]
        c = 0
        while True:
            if c < p_int and c < 32:
                alpha = x + F_p(c)
            else:
                alpha = R_fp.random_element(degree=(1, m_order - 1))
            c += 1
            zeta = power_mod(alpha, exp, P)
            if zeta and all(power_mod(zeta, m_unram // q, P) != 1
                            for q in prime_divs):
                break
        G = R_fp(pari(zeta).Mod(pari(P)).minpoly())

    H_inv = (~F_p(m_unram) * x * G.derivative()) % G
    return G, H_inv


cdef _hensel_lift_factor(G, H_inv, n, p_val, Py_ssize_t e, R_e):
    r"""
    Lift a monic factor ``G`` of `x^n - 1` over `\GF{p}`, with `p \nmid n`, to
    the unique monic factor of `x^n - 1` over `\ZZ / p^e \ZZ` reducing to ``G``.

    ``H_inv`` is the inverse modulo ``G`` over `\GF{p}` of the cofactor
    `(x^n - 1)/G`. It stays fixed during the lifting (single-seed,
    cofactor-free linear Hensel lifting), so each of the `e - 1` steps costs
    one modular exponentiation `x^n \bmod G` in ``R_e``, a polynomial ring over
    `\ZZ / p^e \ZZ` whose arithmetic is delegated to FLINT or NTL.
    """
    cdef list G_list = _coeff_list(G)
    if e <= 1:
        return R_e(G_list)
    cdef Py_ssize_t i
    cdef Py_ssize_t d = len(G_list) - 1
    R_fp = G.parent()
    X = R_e.gen()
    p_int = int(p_val)
    pk = p_int
    for _ in range(1, e):
        # x^n - 1 = 0 mod (p^h, G); extract the next p-adic digit of the error
        r_list = _coeff_list(pow(X, n, R_e(G_list)))
        r_list[0] -= 1
        E = R_fp([c // pk for c in r_list])
        D = (E * H_inv) % G
        for i in range(d):
            G_list[i] += pk * int(D[i])
        pk *= p_int
    return R_e(G_list)


cdef list _coset_factor_fallback(G, Py_ssize_t i, n, p_val, Py_ssize_t e, R_e):
    r"""
    Return the coefficient list of the monic irreducible factor of `x^n - 1`
    over `\ZZ / p^e \ZZ` whose roots are the `\zeta^{i p^j}`, where `\zeta` is
    a root of the seed ``G``.

    This is used when the Newton-Girard reconstruction would require dividing
    by a multiple of `p`: the factor is computed over `\GF{p}` as the minimal
    polynomial of `x^i \bmod G` and then Hensel lifted on its own.
    """
    R_fp = G.parent()
    x = R_fp.gen()
    n_i = n // int(gcd(i, n))
    beta = pow(x, i, G)
    if beta.degree() <= 0:
        g = x - beta[0]
    else:
        g = R_fp(pari(beta).Mod(pari(G)).minpoly())
    H_inv = (~R_fp.base_ring()(n_i) * x * g.derivative()) % g
    return _coeff_list(_hensel_lift_factor(g, H_inv, n_i, p_val, e, R_e))


cdef list _dickson_unramified_c(int64_t p_c, Py_ssize_t m_order,
                                list G_lifted, long m_unram, int64_t final_mod,
                                bint full, G_fp, Py_ssize_t e, R_e):
    """
    Fast C-level trace generation and MED coset factor reconstruction
    for word-sized prime-power moduli.
    """
    cdef long max_T = m_unram if m_unram > m_order else m_order
    cdef int64_t *A = NULL
    cdef int64_t *T = NULL
    cdef char *visited = NULL
    cdef int64_t *med_traces = NULL
    cdef int64_t *e_sym = NULL

    cdef long j, k, i, curr, m_coset, idx
    cdef int64_t coeff, param_sum, term, tail, next_T, scale, scale_inv, raw_trace, sum_val, k_inv
    cdef bint ok
    cdef list factors = []
    cdef list f_coeffs

    try:
        A = <int64_t*>check_calloc(m_order + 1, sizeof(int64_t))
        T = <int64_t*>check_calloc(max_T + 1, sizeof(int64_t))
        visited = <char*>check_calloc(m_unram, sizeof(char))
        med_traces = <int64_t*>check_calloc(m_order + 1, sizeof(int64_t))
        e_sym = <int64_t*>check_calloc(m_order + 1, sizeof(int64_t))

        for j in range(1, m_order + 1):
            coeff = <int64_t>G_lifted[m_order - j] if (m_order - j) < len(G_lifted) else 0
            coeff %= final_mod
            if (j & 1) and coeff:
                coeff = final_mod - coeff
            A[j] = coeff

        T[0] = m_order % final_mod
        if max_T >= 1:
            T[1] = A[1]

        sig_on()
        for k in range(2, m_order + 1):
            param_sum = 0
            for j in range(1, k):
                term = mulmod_i64(A[j], T[k - j], final_mod)
                if not (j & 1):
                    term = final_mod - term if term else 0
                param_sum += term
                if param_sum >= final_mod:
                    param_sum -= final_mod
            tail = mulmod_i64(k % final_mod, A[k], final_mod)
            if not (k & 1):
                tail = final_mod - tail if tail else 0
            T[k] = (param_sum + tail) % final_mod

        for k in range(m_order + 1, m_unram + 1):
            next_T = 0
            for j in range(1, m_order + 1):
                term = mulmod_i64(A[j], T[k - j], final_mod)
                if not (j & 1):
                    term = final_mod - term if term else 0
                next_T += term
                if next_T >= final_mod:
                    next_T -= final_mod
            T[k] = next_T
        sig_off()

        # Index 0 corresponds to (x - 1) = Phi_1(x)
        visited[0] = 1
        if full or m_unram == 1:
            factors.append([final_mod - 1, 1])

        for i in range(1, m_unram):
            if visited[i]:
                continue
            m_coset = 0
            curr = i
            while not visited[curr]:
                visited[curr] = 1
                curr = <long>mulmod_i64(curr, p_c % m_unram, m_unram)
                m_coset += 1

            if not full and gcd_long(i, m_unram) != 1:
                continue

            scale = m_order // m_coset
            scale_inv = invmod_i64(scale, final_mod) if p_c > m_coset else 0

            if p_c > m_coset and scale_inv != 0:
                for k in range(1, m_coset + 1):
                    idx = <long>mulmod_i64(i, k, m_unram)
                    raw_trace = (m_order % final_mod) if idx == 0 else T[idx]
                    med_traces[k] = mulmod_i64(raw_trace, scale_inv, final_mod)

                e_sym[0] = 1
                ok = True
                for k in range(1, m_coset + 1):
                    sum_val = 0
                    for j in range(1, k + 1):
                        term = mulmod_i64(e_sym[k - j], med_traces[j], final_mod)
                        if not (j & 1):
                            term = final_mod - term if term else 0
                        sum_val += term
                        if sum_val >= final_mod:
                            sum_val -= final_mod
                    k_inv = invmod_i64(k, final_mod)
                    if k_inv == 0:
                        ok = False
                        break
                    e_sym[k] = mulmod_i64(sum_val, k_inv, final_mod)

                if ok:
                    f_coeffs = [0] * (m_coset + 1)
                    f_coeffs[m_coset] = 1
                    for k in range(1, m_coset + 1):
                        if k & 1:
                            f_coeffs[m_coset - k] = (final_mod - e_sym[k]) if e_sym[k] else 0
                        else:
                            f_coeffs[m_coset - k] = e_sym[k]
                    factors.append(f_coeffs)
                    continue

            factors.append(_coset_factor_fallback(G_fp, i, m_unram, p_c, e, R_e))
    finally:
        sig_free(A)
        sig_free(T)
        sig_free(visited)
        sig_free(med_traces)
        sig_free(e_sym)

    return factors


cdef list _dickson_unramified_bignum(p_val, Py_ssize_t m_order,
                                     list G_lifted, m_unram, final_mod,
                                     bint full, G_fp, Py_ssize_t e, R_e):
    """
    Arbitrary-precision trace generation and MED coset factor reconstruction
    for prime-power moduli exceeding the word-sized C code path.
    """
    p_int = int(p_val)
    n_int = int(m_unram)
    mod_int = int(final_mod)

    cdef list A = [0] * (m_order + 1)
    cdef Py_ssize_t j, k, i, curr, m_coset, idx
    for j in range(1, m_order + 1):
        coeff = G_lifted[m_order - j] if (m_order - j) < len(G_lifted) else 0
        sign = -1 if (j & 1) else 1
        A[j] = (coeff * sign) % mod_int

    max_T = max(n_int, m_order)
    cdef list T = [0] * (max_T + 1)
    T[0] = m_order % mod_int
    if max_T >= 1:
        T[1] = A[1]

    for k in range(2, m_order + 1):
        param_sum = 0
        for j in range(1, k):
            term = (A[j] * T[k - j]) % mod_int
            if not (j & 1):
                term = -term
            param_sum = (param_sum + term) % mod_int
        tail = (k * A[k]) % mod_int
        if not (k & 1):
            tail = -tail
        T[k] = (param_sum + tail) % mod_int

    for k in range(m_order + 1, n_int + 1):
        next_T = 0
        for j in range(1, m_order + 1):
            term = (A[j] * T[k - j]) % mod_int
            if not (j & 1):
                term = -term
            next_T = (next_T + term) % mod_int
        T[k] = next_T

    cdef list factors = []
    cdef list visited = [False] * n_int
    visited[0] = True
    if full or n_int == 1:
        factors.append([(mod_int - 1) % mod_int, 1])

    for i in range(1, n_int):
        if visited[i]:
            continue
        m_coset = 0
        curr = i
        while not visited[curr]:
            visited[curr] = True
            curr = (curr * p_int) % n_int
            m_coset += 1

        if not full and int(gcd(i, n_int)) != 1:
            continue

        scale = m_order // m_coset
        if p_int > m_coset and scale % p_int != 0:
            scale_inv = pow(scale, -1, mod_int)
            med_traces = [0] * (m_coset + 1)
            for k in range(1, m_coset + 1):
                idx = (i * k) % n_int
                raw_trace = (m_order % mod_int) if idx == 0 else T[idx]
                med_traces[k] = (raw_trace * scale_inv) % mod_int

            e_sym = [0] * (m_coset + 1)
            e_sym[0] = 1
            for k in range(1, m_coset + 1):
                sum_val = 0
                for j in range(1, k + 1):
                    term = (e_sym[k - j] * med_traces[j]) % mod_int
                    if not (j & 1):
                        term = -term
                    sum_val = (sum_val + term) % mod_int
                k_inv = pow(k, -1, mod_int)
                e_sym[k] = (sum_val * k_inv) % mod_int

            f_coeffs = [0] * (m_coset + 1)
            f_coeffs[m_coset] = 1
            for k in range(1, m_coset + 1):
                sign = -1 if (k & 1) else 1
                f_coeffs[m_coset - k] = (sign * e_sym[k]) % mod_int
            factors.append(f_coeffs)
        else:
            factors.append(_coset_factor_fallback(G_fp, i, n_int, p_int, e, R_e))

    return factors


def dickson_cyclotomic_factors(n, p, e=1, bint full=True, ring=None, unit=None):
    r"""
    Return the canonical monic irreducible factorization of `x^n - 1` (if
    ``full=True``) or the `n`-th cyclotomic polynomial `\Phi_n(x)` (if
    ``full=False``) over `\ZZ / p^e \ZZ` using the Dickson Engine.

    Both the unramified regime (`\gcd(n, p) = 1`) and the ramified regime
    (`p \mid n`) are supported for arbitrary prime powers `p^e`.

    INPUT:

    - ``n`` -- positive integer
    - ``p`` -- prime number
    - ``e`` -- positive integer (default: ``1``)
    - ``full`` -- boolean (default: ``True``); if ``True``, factor `x^n - 1`;
      if ``False``, factor `\Phi_n(x)`
    - ``ring`` -- optional univariate polynomial ring over `\ZZ / p^e \ZZ`
      (or `\GF{p}`); if provided, a :class:`~sage.structure.factorization.Factorization`
      object in ``ring`` is returned, otherwise a list of pairs
      ``(coeffs, multiplicity)`` is returned
    - ``unit`` -- optional unit element in the base ring of ``ring`` for scaling
      the returned :class:`~sage.structure.factorization.Factorization`

    OUTPUT:

    A :class:`~sage.structure.factorization.Factorization` object in ``ring``
    if ``ring`` is given, or a list of ``(coeffs, multiplicity)`` tuples where
    ``coeffs`` is the coefficient list ``[c_0, c_1, ..., 1]`` of a monic factor.

    ALGORITHM:

    Write `n = p^k m` with `\gcd(m, p) = 1`:

    1. **Unramified core (`\gcd(m, p) = 1`):**
       Find a single primitive irreducible factor `G_1(x) \in \GF{p}[x]` of
       `\Phi_m(x)` of degree `d = \mathrm{ord}_m(p)` and lift it to
       `\ZZ / p^e \ZZ[x]` by cofactor-free single-seed linear Hensel lifting:
       the inverse `H_1(x)^{-1} \bmod G_1(x)` of the cofactor is computed once
       over `\GF{p}` as `m^{-1} x G_1'(x) \bmod G_1(x)`, and each of the
       `e - 1` lifting steps costs one modular exponentiation `x^m \bmod G_1`
       in degree `d` [WDYH2026a]_, [WDYH2026b]_. Next, generate the power-sum
       trace sequence `T_1, \dots, T_m \pmod{p^e}` via the order-`d` Dickson
       linear recurrence and reconstruct all irreducible factors across
       `p`-cyclotomic cosets using Multiple Equal-Difference (MED) coset
       partitioning [ZZLW2025]_ and Newton-Girard identities. When a
       Newton-Girard denominator is divisible by `p` (which can only happen
       if `p \le d`), the factor attached to the coset is instead computed
       over `\GF{p}` and Hensel lifted on its own.

    2. **Ramified extension (`k \ge 1`):**
       When `e \ge 2`, by the linear disjointness of unramified and totally
       ramified cyclotomic extensions of `\QQ_p`, each unramified factor
       `G_i(x)` of `x^m - 1` yields canonical monic irreducible factors
       `H_i^{(j)}(x) = G_i(x^{p^j}) / G_i(x^{p^{j-1}})` over `\ZZ / p^e \ZZ`
       for `1 \le j \le k` via exact polynomial division and substitution.
       When `e = 1` (over `\GF{p}`), `G_i(x)` appears with multiplicity `p^k`
       in `x^n - 1` and `\varphi(p^k)` in `\Phi_n(x)`.

    .. NOTE::

        Factorizations over `\ZZ / p^e \ZZ` with `e \ge 2` into irreducible
        polynomials are in general not unique. The factors returned here are
        the reductions modulo `p^e` of the monic irreducible factors of
        `x^n - 1` over the `p`-adic integers `\ZZ_p`.

    REFERENCES:

    - [ZZLW2025]_
    - [WDYH2026a]_
    - [WDYH2026b]_

    EXAMPLES:

    Unramified factorization of `x^n - 1` over `\ZZ / p^e \ZZ`::

        sage: from sage.rings.polynomial.cyclotomic import dickson_cyclotomic_factors
        sage: R.<x> = Zmod(5^3)[]
        sage: F = dickson_cyclotomic_factors(6, 5, 3, ring=R); F
        (x + 1) * (x + 124) * (x^2 + x + 1) * (x^2 + 124*x + 1)
        sage: F.expand() == x^6 - 1
        True

    Factoring a single cyclotomic polynomial `\Phi_n(x)` (``full=False``)::

        sage: R.<x> = Zmod(11^3)[]
        sage: F = dickson_cyclotomic_factors(12, 11, 3, full=False, ring=R); F
        (x^2 + 578*x + 1) * (x^2 + 753*x + 1)
        sage: F.expand() == R.cyclotomic_polynomial(12)
        True

    Small characteristic (`p \le \mathrm{ord}_n(p)`), where some factors are
    computed over `\GF{p}` and lifted separately::

        sage: R.<x> = Zmod(2^4)[]
        sage: F = dickson_cyclotomic_factors(15, 2, 4, ring=R)
        sage: len(F)
        5
        sage: F.expand() == x^15 - 1
        True

    Ramified regime (`p \mid n`) over `\ZZ / p^e \ZZ` via cyclotomic substitution::

        sage: R.<x> = Zmod(3^2)[]
        sage: F = dickson_cyclotomic_factors(12, 3, 2, ring=R); F
        (x + 1) * (x + 8) * (x^2 + 1) * (x^2 + x + 1) * (x^2 + 8*x + 1) * (x^4 + 8*x^2 + 1)
        sage: F.expand() == x^12 - 1
        True

    Large prime-power moduli exceeding 64-bit machine words::

        sage: R.<x> = Zmod(101^10)[]
        sage: F = dickson_cyclotomic_factors(102, 101, 10, ring=R)
        sage: len(F)
        52
        sage: F.expand() == x^102 - 1
        True

    Over `\ZZ / 8 \ZZ` we have `x^2 - 1 = (x - 1)(x + 1) = (x - 3)(x + 3)`;
    the factorization coming from `\ZZ_2` is returned::

        sage: R.<x> = Zmod(8)[]
        sage: dickson_cyclotomic_factors(2, 2, 3, ring=R)
        (x + 1) * (x + 7)
        sage: (x + 3) * (x + 5) == x^2 - 1
        True

    TESTS::

        sage: R.<x> = GF(3)[]
        sage: F = dickson_cyclotomic_factors(12, 3, 1, ring=R); F
        (x + 1)^3 * (x + 2)^3 * (x^2 + 1)^3
        sage: F.expand() == x^12 - 1
        True
        sage: F_cyclo = dickson_cyclotomic_factors(12, 3, 1, full=False, ring=R); F_cyclo
        (x^2 + 1)^2
        sage: F_cyclo.expand() == R.cyclotomic_polynomial(12)
        True

    Large multiplicative order `\mathrm{ord}_n(p)`::

        sage: R.<x> = Zmod(257^5)[]
        sage: F = dickson_cyclotomic_factors(509, 257, 5, ring=R)
        sage: [f.degree() for f, _ in F]
        [1, 508]
        sage: F.expand() == x^509 - 1
        True
        sage: R.<x> = Zmod(2^3)[]
        sage: F = dickson_cyclotomic_factors(4095, 2, 3, ring=R)
        sage: len(F), F.expand() == x^4095 - 1
        (351, True)

    Consistency checks::

        sage: all(dickson_cyclotomic_factors(n, p, e, ring=Zmod(p^e)['y']).expand()
        ....:     == Zmod(p^e)['y'].gen()^n - 1
        ....:     for p in primes(14) for e in [1, 2, 3] for n in [1..30])
        True
        sage: all(dickson_cyclotomic_factors(n, p, e, full=False, ring=Zmod(p^e)['y']).expand()
        ....:     == Zmod(p^e)['y'].cyclotomic_polynomial(n)
        ....:     for p in primes(14) for e in [1, 2, 3] for n in [1..30])
        True
        sage: dickson_cyclotomic_factors(6, 5, 2, ring=Zmod(5)['x'])
        Traceback (most recent call last):
        ...
        ValueError: ring must be a polynomial ring over Z/p^e Z
    """
    n_zz = ZZ(n)
    p_zz = ZZ(p)
    e_zz = ZZ(e)
    if n_zz < 1:
        raise ValueError("n must be a positive integer")
    if not p_zz.is_prime():
        raise ValueError("p must be a prime number")
    if e_zz < 1:
        raise ValueError("e must be a positive integer")

    cdef Py_ssize_t k_ram, j, e_int = <Py_ssize_t>e_zz
    k_ram, m_unram_zz = n_zz.val_unit(p_zz)
    final_mod_zz = p_zz ** e_zz
    final_mod_int = int(final_mod_zz)
    p_int = int(p_zz)
    m_unram_int = int(m_unram_zz)

    if ring is None:
        from sage.rings.finite_rings.integer_mod_ring import Zmod
        R_e = Zmod(final_mod_zz)['x']
    else:
        if ring.base_ring().characteristic() != final_mod_zz:
            raise ValueError("ring must be a polynomial ring over Z/p^e Z")
        R_e = ring

    cdef list unramified_factors
    cdef Py_ssize_t m_order
    if m_unram_int == 1:
        unramified_factors = [[final_mod_int - 1, 1]]
    else:
        m_order = <Py_ssize_t>int(Mod(p_zz, m_unram_zz).multiplicative_order())
        G_fp, H_inv = _find_primitive_seed_fp(p_zz, m_unram_int, m_order)
        G_lifted = _coeff_list(_hensel_lift_factor(G_fp, H_inv, m_unram_int,
                                                   p_zz, e_int, R_e))
        if final_mod_zz.nbits() <= DICKSON_WORD_BITS and m_unram_int < (1 << 30):
            unramified_factors = _dickson_unramified_c(
                <int64_t>p_int, m_order, G_lifted, <long>m_unram_int,
                <int64_t>final_mod_int, full, G_fp, e_int, R_e)
        else:
            unramified_factors = _dickson_unramified_bignum(
                p_int, m_order, G_lifted, m_unram_int, final_mod_int, full,
                G_fp, e_int, R_e)

    cdef list result_pairs = []
    cdef list G, H1
    if k_ram == 0:
        for G in unramified_factors:
            result_pairs.append((G, 1))
    elif e_int == 1:
        mult = p_int ** k_ram if full else (p_int ** (k_ram - 1)) * (p_int - 1)
        for G in unramified_factors:
            result_pairs.append((G, mult))
    else:
        for G in unramified_factors:
            # first ramified layer G(x^p) / G(x), by exact monic division
            H1 = _coeff_list(R_e(_poly_substitute_power(G, p_int)) // R_e(G))
            if full:
                result_pairs.append((G, 1))
                result_pairs.append((H1, 1))
                for j in range(2, k_ram + 1):
                    result_pairs.append((_poly_substitute_power(H1, p_int ** (j - 1)), 1))
            elif k_ram == 1:
                result_pairs.append((H1, 1))
            else:
                result_pairs.append((_poly_substitute_power(H1, p_int ** (k_ram - 1)), 1))

    if ring is None:
        return result_pairs

    if unit is None:
        unit = ring.base_ring().one()
    F = Factorization([(ring(coeffs), mult) for coeffs, mult in result_pairs], unit=unit)
    F.sort()
    return F
