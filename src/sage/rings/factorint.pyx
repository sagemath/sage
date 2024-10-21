# sage_setup: distribution = sagemath-categories
r"""
Integer factorization functions

AUTHORS:

- Andre Apitzsch (2011-01-13): initial version
"""

# ****************************************************************************
#       Copyright (C) 2010-2011 André Apitzsch <andre.apitzsch@st.ovgu.de>
#                     2012      Nils Bruin
#                     2014      David Roe
#                     2014      Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport *

from sage.rings.integer cimport Integer
from sage.structure.factorization_integer import IntegerFactorization
from sage.misc.misc_c import prod

cdef extern from "limits.h":
    long LONG_MAX

cpdef aurifeuillian(n, m, F=None, bint check=True):
    r"""
    Return the Aurifeuillian factors `F_n^\pm(m^2n)`.

    This is based off Theorem 3 of [Bre1993]_.

    INPUT:

    - ``n`` -- integer
    - ``m`` -- integer
    - ``F`` -- integer (default: ``None``)
    - ``check`` -- boolean (default: ``True``)

    OUTPUT: list of factors

    EXAMPLES::

        sage: from sage.rings.factorint import aurifeuillian

        sage: # needs sage.libs.pari sage.rings.real_interval_field
        sage: aurifeuillian(2, 2)
        [5, 13]
        sage: aurifeuillian(2, 2^5)
        [1985, 2113]
        sage: aurifeuillian(5, 3)
        [1471, 2851]
        sage: aurifeuillian(15, 1)
        [19231, 142111]

        sage: # needs sage.libs.pari
        sage: aurifeuillian(12, 3)
        Traceback (most recent call last):
        ...
        ValueError: n has to be square-free
        sage: aurifeuillian(1, 2)
        Traceback (most recent call last):
        ...
        ValueError: n has to be greater than 1
        sage: aurifeuillian(2, 0)
        Traceback (most recent call last):
        ...
        ValueError: m has to be positive

    .. NOTE::

        There is no need to set `F`. It's only for increasing speed
        of :meth:`factor_aurifeuillian()`.
    """
    if check:
        if not n.is_squarefree():
            raise ValueError("n has to be square-free")
        if n < 2:
            raise ValueError("n has to be greater than 1")
        if m < 1:
            raise ValueError("m has to be positive")

    from sage.arith.misc import euler_phi
    from sage.rings.real_mpfi import RealIntervalField

    x = m**2*n
    cdef Py_ssize_t y = euler_phi(2*n)//2
    if F is None:
        from sage.rings.polynomial.cyclotomic import cyclotomic_value
        if n%2:
            if n%4 == 3:
                s = -1
            else:
                s = 1
            F = cyclotomic_value(n, s*x)
        else:
            F = cyclotomic_value(n//2, -x**2)
            if n == 2:
                F = -F
    cdef Py_ssize_t j
    tmp = 0
    for j in range(y):
        tmp += n.kronecker(2*j + 1) / ((2*j + 1) * x**j)
    prec = Integer(150)
    R = RealIntervalField(prec)
    Fm = R(F).sqrt() * R(-1/m*tmp).exp()
    while Fm.upper().round() != Fm.lower().round():
        prec *= 2
        R = RealIntervalField(prec)
        Fm = R(F).sqrt() * R(-1/m*tmp).exp()
    Fm = Fm.upper().round()
    assert (not check or Fm.divides(F))
    return [Fm, F // Fm]

cpdef factor_aurifeuillian(n, check=True):
    r"""
    Return Aurifeuillian factors of `n` if `n = x^{(2k-1)x} \pm 1`
    (where the sign is '-' if x = 1 mod 4, and '+' otherwise) else `n`

    INPUT:

    - ``n`` -- integer

    OUTPUT: list of factors of `n` found by Aurifeuillian factorization

    EXAMPLES::

        sage: # needs sage.libs.pari sage.rings.real_interval_field
        sage: from sage.rings.factorint import factor_aurifeuillian as fa
        sage: fa(2^6 + 1)
        [5, 13]
        sage: fa(2^58 + 1)
        [536838145, 536903681]
        sage: fa(3^3 + 1)
        [4, 1, 7]
        sage: fa(5^5 - 1)
        [4, 11, 71]
        sage: prod(_) == 5^5 - 1
        True
        sage: fa(2^4 + 1)
        [17]
        sage: fa((6^2*3)^3 + 1)
        [109, 91, 127]

    TESTS::

        sage: # needs sage.libs.pari sage.rings.real_interval_field
        sage: for n in [2,3,5,6,30,31,33]:
        ....:     for m in [8,96,109201283]:
        ....:         s = -1 if n % 4 == 1 else 1
        ....:         y = (m^2*n)^n + s
        ....:         F = fa(y)
        ....:         assert(F and prod(F) == y)

    REFERENCES:

    - http://mathworld.wolfram.com/AurifeuilleanFactorization.html
    - [Bre1993]_ Theorem 3
    """
    if n in [-2, -1, 0, 1, 2]:
        return [n]
    cdef int exp = 1
    for s in [-1, 1]:
        x, exp = (n - s).perfect_power()
        if exp > 1:
            # factorization is possible if x = m^2 * n and exp = n, n squarefree
            # We have the freedom to replace exp by exp/a and x by x^a.  If a is even,
            # n = 1 and we've gotten nowhere.  So set n as the squarefree part of x,
            # and m^2 the remainder.  We can replace exp by exp/(2a+1) as long as we
            # replace (m^2 * n) by (m^2 * n)^(2a + 1) = (m^(2a + 1) * n^a)^2 * n.
            # In particular, n needs to be a divisor of both x and exp (as well as being squarefree).
            a_guess = x.gcd(exp)
            m = x // a_guess
            # n_guess is small, so we can factor it.
            a = 1
            for p, e in a_guess.factor():
                v = m.valuation(p)
                if v + e % 2:
                    # one factor of p remains in a.
                    a *= p
                    if e > 1:
                        m *= p**(e-1)
                else:
                    # all factors of p go to m
                    m *= p**e
            if a == 1 or (a % 4 == 1 and s == 1) or (a % 4 != 1 and s == -1): continue
            m, r = m.sqrtrem()
            if r: continue
            exp_adjust = exp // a
            if exp_adjust % 2 == 0: continue
            m = m**exp_adjust * a**(exp_adjust//2)
            F = aurifeuillian(a, m, check=False)
            rem = prod(F)
            if check and not rem.divides(n):
                raise RuntimeError("rem=%s, F=%s, n=%s, m=%s"%(rem, F, n, m))
            rem = n // rem
            if rem != 1:
                return [rem] + F
            return F
    return [n]


def factor_cunningham(m, proof=None):
    r"""
    Return factorization of ``self`` obtained using trial division
    for all primes in the so called Cunningham table. This is
    efficient if ``self`` has some factors of type `b^n+1` or `b^n-1`,
    with `b` in `\{2,3,5,6,7,10,11,12\}`.

    You need to install an optional package to use this method,
    this can be done with the following command line:
    ``sage -i cunningham_tables``.

    INPUT:

    - ``proof`` -- boolean (default: ``None``); whether or not to
      prove primality of each factor, this is only for factors
      not in the Cunningham table

    EXAMPLES::

        sage: from sage.rings.factorint import factor_cunningham
        sage: factor_cunningham(2^257-1) # optional - cunningham_tables
        535006138814359 * 1155685395246619182673033 * 374550598501810936581776630096313181393
        sage: factor_cunningham((3^101+1)*(2^60).next_prime(), proof=False) # optional - cunningham_tables
        2^2 * 379963 * 1152921504606847009 * 1017291527198723292208309354658785077827527
    """
    from sage.databases import cunningham_tables
    cunningham_prime_factors = cunningham_tables.cunningham_prime_factors()
    if m.nbits() < 100 or not cunningham_prime_factors:
        return m.factor(proof=proof)
    n = Integer(m)
    L = []
    for p in cunningham_prime_factors:
        if p > n:
            break
        if p.divides(n):
            v,n = n.val_unit(p)
            L.append( (p,v) )
    if n.is_one():
        return IntegerFactorization(L)
    else:
        return IntegerFactorization(L)*n.factor(proof=proof)


cpdef factor_trial_division(m, long limit=LONG_MAX):
    r"""
    Return partial factorization of ``self`` obtained using trial division
    for all primes up to ``limit``, where ``limit`` must fit in a C ``signed long``.

    INPUT:

    - ``limit`` -- integer (default: ``LONG_MAX``); that fits in a C ``signed long``

    EXAMPLES::

        sage: from sage.rings.factorint import factor_trial_division
        sage: n = 920384092842390423848290348203948092384082349082
        sage: factor_trial_division(n, 1000)
        2 * 11 * 41835640583745019265831379463815822381094652231
        sage: factor_trial_division(n, 2000)
        2 * 11 * 1531 * 27325696005058797691594630609938486205809701

    TESTS:

    Test that :issue:`13692` is solved::

        sage: from sage.rings.factorint import factor_trial_division
        sage: list(factor_trial_division(8))
        [(2, 3)]
    """
    cdef Integer n = PY_NEW(Integer), unit = PY_NEW(Integer), p = Integer(2)
    cdef long e

    n = Integer(m)
    if mpz_sgn(n.value) > 0:
        mpz_set_si(unit.value, 1)
    else:
        mpz_neg(n.value,n.value)
        mpz_set_si(unit.value, -1)

    F = []
    while mpz_cmpabs_ui(n.value, 1):
        p = n.trial_division(bound=limit,start=mpz_get_ui(p.value))
        e = mpz_remove(n.value, n.value, p.value)
        F.append((p,e))

    return IntegerFactorization(F, unit=unit, unsafe=True,
                                   sort=False, simplify=False)
