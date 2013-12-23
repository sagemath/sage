r"""
Integer factorization functions

AUTHORS:

- Andre Apitzsch (2011-01-13): initial version

"""

#*****************************************************************************
#       Copyright (C) 2010 André Apitzsch <andre.apitzsch@st.ovgu.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/libs/pari/decl.pxi"
include "sage/ext/gmp.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.integer cimport Integer
from sage.rings.fast_arith import prime_range
from sage.structure.factorization_integer import IntegerFactorization
from math import floor
from sage.misc.superseded import deprecated_function_alias

cdef extern from "limits.h":
    long LONG_MAX

cpdef aurifeuillian(n, m, F=None):
    r"""
    Return Aurifeuillian factors `F_n^\pm(m^2n)`

    INPUT:

    - ``n`` - integer
    - ``m`` - integer
    - ``F`` - integer (default: None)

    OUTPUT:

        List of factors

    EXAMPLES:

        sage: from sage.rings.factorint import aurifeuillian
        sage: aurifeuillian(2,2)
        [5, 13]
        sage: aurifeuillian(2,2^5)
        [1985, 2113]
        sage: aurifeuillian(5,3)
        [1471, 2851]
        sage: aurifeuillian(15,1)
        [19231, 142111]
        sage: aurifeuillian(12,3)
        Traceback (most recent call last):
        ...
        ValueError: n has to be square-free
        sage: aurifeuillian(1,2)
        Traceback (most recent call last):
        ...
        ValueError: n has to be greater than 1
        sage: aurifeuillian(2,0)
        Traceback (most recent call last):
        ...
        ValueError: m has to be positive

    .. NOTE::

    There is no need to set `F`. It's only for increasing speed
    of factor_aurifeuillian().

    REFERENCES:

    .. Brent, On computing factors of cyclotomic polynomials, Theorem 3
        arXiv:1004.5466v1 [math.NT]

    """
    from sage.rings.arith import euler_phi, kronecker_symbol
    from sage.functions.log import exp
    from sage.rings.real_mpfr import RealField
    if not n.is_squarefree():
        raise ValueError, "n has to be square-free"
    if n < 2:
        raise ValueError, "n has to be greater than 1"
    if m < 1:
        raise ValueError, "m has to be positive"
    x = m**2*n
    y = euler_phi(2*n)//2
    if F == None:
        from sage.misc.functional import cyclotomic_polynomial
        if n%2:
            if n%4 == 3:
                s = -1
            else:
                s = 1
            F = cyclotomic_polynomial(n)(s*x)
        else:
            F = (-1)**euler_phi(n//2)*cyclotomic_polynomial(n//2)(-x**2)
    tmp = sum([kronecker_symbol(n,2*j+1)/((2*j+1)*x**j) for j in range(y)])
    R = RealField(300)
    Fm = R(F.sqrt()*R(-1/m*tmp).exp()).round()
    return [Fm, Integer(round(F//Fm))]

base_exponent = deprecated_function_alias(12116, lambda n: n.perfect_power())

cpdef factor_aurifeuillian(n):
    r"""
    Return Aurifeuillian factors of `n` if `n = x^{(2k-1)x} \pw 1`
    (where the sign is '-' if x = 1 mod 4, and '+' otherwise) else `n`

    INPUT:

    - ``n`` - integer

    OUTPUT:

        List of factors of `n` found by Aurifeuillian factorization.

    EXAMPLES:

        sage: from sage.rings.factorint import factor_aurifeuillian as fa
        sage: fa(2^6+1)
        [5, 13]
        sage: fa(2^58+1)
        [536838145, 536903681]
        sage: fa(3^3+1)
        [4, 1, 7]
        sage: fa(5^5-1)
        [4, 11, 71]
        sage: prod(_) == 5^5-1
        True
        sage: fa(2^4+1)
        [17]

    REFERENCES:

    .. http://mathworld.wolfram.com/AurifeuilleanFactorization.html

    """
    if n in [-2, -1, 0, 1, 2]:
        return [n]
    cdef int exp = 1
    for x in [-1, 1]:
        b = n + x
        b, exp = b.perfect_power()
        if exp > 1:
            if not b.is_prime():
                continue
            m = b**((exp+b)/(2*b)-1)
            if m != floor(m):
                continue
            if b == 2:
                if x == -1:
                    return aurifeuillian(b, m, n)
                else:
                    return [2**(exp/2)-1, 2**(exp/2)+1]
            F = b**(exp/b)-x
            result = [F] + aurifeuillian(b, m, n/F)
            prod = 1
            for a in result:
                prod *= a
            if prod == n:
                return result
    return [n]

def factor_cunningham(m, proof=None):
    r"""
    Return factorization of self obtained using trial division
    for all primes in the so called Cunningham table. This is
    efficient if self has some factors of type $b^n+1$ or $b^n-1$,
    with $b$ in $\{2,3,5,6,7,10,11,12\}$.

    You need to install an optional package to use this method,
    this can be done with the following command line:
    ``sage -i cunningham_tables``

    INPUT:

     - ``proof`` - bool (default: None) whether or not to
       prove primality of each factor, this is only for factors
       not in the Cunningham table.

    EXAMPLES::

        sage: from sage.rings.factorint import factor_cunningham
        sage: factor_cunningham(2^257-1) # optional - cunningham
        535006138814359 * 1155685395246619182673033 * 374550598501810936581776630096313181393
        sage: factor_cunningham((3^101+1)*(2^60).next_prime(),proof=False) # optional - cunningham
        2^2 * 379963 * 1152921504606847009 * 1017291527198723292208309354658785077827527

    """
    from sage.databases import cunningham_tables
    cunningham_prime_factors = cunningham_tables.cunningham_prime_factors()
    if m.nbits() < 100 or len(cunningham_prime_factors) == 0:
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
    Return partial factorization of self obtained using trial division
    for all primes up to limit, where limit must fit in a C signed long.

    INPUT:

     - ``limit`` - integer (default: LONG_MAX) that fits in a C signed long

    EXAMPLES::

        sage: from sage.rings.factorint import factor_trial_division
        sage: n = 920384092842390423848290348203948092384082349082
        sage: factor_trial_division(n, 1000)
        2 * 11 * 41835640583745019265831379463815822381094652231
        sage: factor_trial_division(n, 2000)
        2 * 11 * 1531 * 27325696005058797691594630609938486205809701

    TESTS:

    Test that :trac:`13692` is solved::

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

cpdef factor_using_pari(n, int_=False, debug_level=0, proof=None):
    r"""
    Factors this (positive) integer using PARI.

    This method returns a list of pairs, not a ``Factorization``
    object. The first element of each pair is the factor, of type
    ``Integer`` if ``int_`` is ``False`` or ``int`` otherwise,
    the second element is the positive exponent, of type ``int``.

    INPUT:

     - ``int_`` -- (default: ``False``), whether the factors are
       of type ``int`` instead of ``Integer``

     - ``debug_level`` -- (default: 0), debug level of the call
       to PARI

     - ``proof`` -- (default: ``None``), whether the factors are
       required to be proven prime;  if ``None``, the global default
       is used

    OUTPUT:

        -  a list of pairs

    EXAMPLES::

        sage: factor(-2**72 + 3, algorithm='pari')  # indirect doctest
        -1 * 83 * 131 * 294971519 * 1472414939
    """
    from sage.libs.pari.all import pari

    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "arithmetic")

    prev = pari.get_debug_level()

    if prev != debug_level:
        pari.set_debug_level(debug_level)

    F = pari(n).factor(proof=proof)
    B, e = F
    if int_:
        v = [(int(B[i]), int(e[i])) for i in xrange(len(B))]
    else:
        v = [(Integer(B[i]), int(e[i])) for i in xrange(len(B))]

    if prev != debug_level:
        pari.set_debug_level(prev)

    return v
