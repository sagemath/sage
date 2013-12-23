"""
Miscellaneous arithmetic functions
"""
#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import math
import sys
import sage.misc.misc as misc
from sage.libs.pari.all import pari
import sage.libs.flint.arith as flint_arith

from sage.rings.rational_field import QQ
import sage.rings.rational
import sage.rings.complex_field
import sage.rings.complex_number
import sage.rings.real_mpfr
from sage.structure.element import parent
from sage.misc.misc import prod, union
from sage.rings.real_mpfi import RealIntervalField

import fast_arith

from integer_ring import ZZ
import integer

##################################################################
# Elementary Arithmetic
##################################################################

def algdep(z, degree, known_bits=None, use_bits=None, known_digits=None, use_digits=None, height_bound=None, proof=False):
    """
    Returns a polynomial of degree at most `degree` which is
    approximately satisfied by the number `z`. Note that the returned
    polynomial need not be irreducible, and indeed usually won't be if
    `z` is a good approximation to an algebraic number of degree less
    than `degree`.

    You can specify the number of known bits or digits of `z` with
    ``known_bits=k`` or ``known_digits=k``. PARI is then told to
    compute the result using `0.8k` of these bits/digits. Or, you can
    specify the precision to use directly with ``use_bits=k`` or
    ``use_digits=k``. If none of these are specified, then the precision
    is taken from the input value.

    A height bound may be specified to indicate the maximum coefficient
    size of the returned polynomial; if a sufficiently small polynomial
    is not found, then ``None`` will be returned. If ``proof=True`` then
    the result is returned only if it can be proved correct (i.e. the
    only possible minimal polynomial satisfying the height bound, or no
    such polynomial exists). Otherwise a ``ValueError`` is raised
    indicating that higher precision is required.

    ALGORITHM: Uses LLL for real/complex inputs, PARI C-library
    ``algdep`` command otherwise.

    Note that ``algebraic_dependency`` is a synonym for ``algdep``.

    INPUT:


    -  ``z`` - real, complex, or `p`-adic number

    -  ``degree`` - an integer

    -  ``height_bound`` - an integer (default: ``None``) specifying the maximum
                          coefficient size for the returned polynomial

    -  ``proof`` - a boolean (default: ``False``), requires height_bound to be set


    EXAMPLES::

        sage: algdep(1.888888888888888, 1)
        9*x - 17
        sage: algdep(0.12121212121212,1)
        33*x - 4
        sage: algdep(sqrt(2),2)
        x^2 - 2

    This example involves a complex number::

        sage: z = (1/2)*(1 + RDF(sqrt(3)) *CC.0); z
        0.500000000000000 + 0.866025403784439*I
        sage: p = algdep(z, 6); p
        x^3 + 1
        sage: p.factor()
        (x + 1) * (x^2 - x + 1)
        sage: z^2 - z + 1
        0.000000000000000

    This example involves a `p`-adic number::

        sage: K = Qp(3, print_mode = 'series')
        sage: a = K(7/19); a
        1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
        sage: algdep(a, 1)
        19*x - 7

    These examples show the importance of proper precision control. We
    compute a 200-bit approximation to `sqrt(2)` which is wrong in the
    33'rd bit::

        sage: z = sqrt(RealField(200)(2)) + (1/2)^33
        sage: p = algdep(z, 4); p
        227004321085*x^4 - 216947902586*x^3 - 99411220986*x^2 + 82234881648*x - 211871195088
        sage: factor(p)
        227004321085*x^4 - 216947902586*x^3 - 99411220986*x^2 + 82234881648*x - 211871195088
        sage: algdep(z, 4, known_bits=32)
        x^2 - 2
        sage: algdep(z, 4, known_digits=10)
        x^2 - 2
        sage: algdep(z, 4, use_bits=25)
        x^2 - 2
        sage: algdep(z, 4, use_digits=8)
        x^2 - 2

    Using the ``height_bound`` and ``proof`` parameters, we can see that
    `pi` is not the root of an integer polynomial of degree at most 5
    and coefficients bounded above by 10::

        sage: algdep(pi.n(), 5, height_bound=10, proof=True) is None
        True

    For stronger results, we need more precicion::

        sage: algdep(pi.n(), 5, height_bound=100, proof=True) is None
        Traceback (most recent call last):
        ...
        ValueError: insufficient precision for non-existence proof
        sage: algdep(pi.n(200), 5, height_bound=100, proof=True) is None
        True

        sage: algdep(pi.n(), 10, height_bound=10, proof=True) is None
        Traceback (most recent call last):
        ...
        ValueError: insufficient precision for non-existence proof
        sage: algdep(pi.n(200), 10, height_bound=10, proof=True) is None
        True

    We can also use ``proof=True`` to get positive results::

        sage: a = sqrt(2) + sqrt(3) + sqrt(5)
        sage: algdep(a.n(), 8, height_bound=1000, proof=True)
        Traceback (most recent call last):
        ...
        ValueError: insufficient precision for uniqueness proof
        sage: f = algdep(a.n(1000), 8, height_bound=1000, proof=True); f
        x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
        sage: f(a).expand()
        0

    TESTS::

        sage: algdep(complex("1+2j"), 4)
        x^2 - 2*x + 5
    """
    if proof and not height_bound:
        raise ValueError, "height_bound must be given for proof=True"

    x = ZZ['x'].gen()

    if isinstance(z, (int, long, integer.Integer)):
        if height_bound and abs(z) >= height_bound:
            return None
        return x - ZZ(z)

    degree = ZZ(degree)

    if isinstance(z, (sage.rings.rational.Rational)):
        if height_bound and max(abs(z.denominator()), abs(z.numerator())) >= height_bound:
            return None
        return z.denominator()*x - z.numerator()

    if isinstance(z, float):
        z = sage.rings.all.RR(z)
    elif isinstance(z, complex):
        z = sage.rings.all.CC(z)

    if isinstance(z, (sage.rings.real_mpfr.RealNumber,
                      sage.rings.complex_number.ComplexNumber)):

        log2_10 = 3.32192809488736

        prec = z.prec() - 6
        if known_digits is not None:
            known_bits = known_digits * log2_10
        if known_bits is not None:
            use_bits = known_bits * 0.8
        if use_digits is not None:
            use_bits = use_digits * log2_10
        if use_bits is not None:
            prec = int(use_bits)

        is_complex = isinstance(z, sage.rings.complex_number.ComplexNumber)
        n = degree+1
        from sage.matrix.all import matrix
        M = matrix(ZZ, n, n+1+int(is_complex))
        r = ZZ(1) << prec
        M[0, 0] = 1
        M[0,-1] = r
        for k in range(1, degree+1):
            M[k,k] = 1
            r *= z
            if is_complex:
                M[k, -1] = r.real().round()
                M[k, -2] = r.imag().round()
            else:
                M[k, -1] = r.round()
        LLL = M.LLL(delta=.75)
        coeffs = LLL[0][:n]
        if height_bound:
            def norm(v):
                # norm on an integer vector invokes Integer.sqrt() which tries to factor...
                from sage.rings.real_mpfi import RIF
                return v.change_ring(RIF).norm()
            if max(abs(a) for a in coeffs) > height_bound:
                if proof:
                    # Given an LLL reduced basis $b_1, ..., b_n$, we only
                    # know that $|b_1| <= 2^((n-1)/2) |x|$ for non-zero $x \in L$.
                    if norm(LLL[0]) <= 2**((n-1)/2) * n.sqrt() * height_bound:
                        raise ValueError, "insufficient precision for non-existence proof"
                return None
            elif proof and norm(LLL[1]) < 2**((n-1)/2) * max(norm(LLL[0]), n.sqrt()*height_bound):
                raise ValueError, "insufficient precision for uniqueness proof"
        if coeffs[degree] < 0:
            coeffs = -coeffs
        f = list(coeffs)

    elif proof or height_bound:
        raise NotImplementedError, "proof and height bound only implemented for real and complex numbers"

    else:
        y = pari(z)
        f = y.algdep(degree)

    return x.parent()(f)


algebraic_dependency = algdep

def bernoulli(n, algorithm='default', num_threads=1):
    r"""
    Return the n-th Bernoulli number, as a rational number.

    INPUT:

    -  ``n`` - an integer
    -  ``algorithm``:

       -  ``'default'`` - (default) use 'pari' for n <= 30000, and 'bernmm' for
          n > 30000 (this is just a heuristic, and not guaranteed to be optimal
          on all hardware)
       -  ``'pari'`` - use the PARI C library
       -  ``'gap'`` - use GAP
       -  ``'gp'`` - use PARI/GP interpreter
       -  ``'magma'`` - use MAGMA (optional)
       -  ``'bernmm'`` - use bernmm package (a multimodular algorithm)

    -  ``num_threads`` - positive integer, number of
       threads to use (only used for bernmm algorithm)

    EXAMPLES::

        sage: bernoulli(12)
        -691/2730
        sage: bernoulli(50)
        495057205241079648212477525/66

    We demonstrate each of the alternative algorithms::

        sage: bernoulli(12, algorithm='gap')
        -691/2730
        sage: bernoulli(12, algorithm='gp')
        -691/2730
        sage: bernoulli(12, algorithm='magma')           # optional - magma
        -691/2730
        sage: bernoulli(12, algorithm='pari')
        -691/2730
        sage: bernoulli(12, algorithm='bernmm')
        -691/2730
        sage: bernoulli(12, algorithm='bernmm', num_threads=4)
        -691/2730

    TESTS::

        sage: algs = ['gap','gp','pari','bernmm']
        sage: test_list = [ZZ.random_element(2, 2255) for _ in range(500)]
        sage: vals = [[bernoulli(i,algorithm = j) for j in algs] for i in test_list]  # long time (up to 21s on sage.math, 2011)
        sage: union([len(union(x))==1 for x in vals])  # long time (depends on previous line)
        [True]
        sage: algs = ['gp','pari','bernmm']
        sage: test_list = [ZZ.random_element(2256, 5000) for _ in range(500)]
        sage: vals = [[bernoulli(i,algorithm = j) for j in algs] for i in test_list]  # long time (up to 30s on sage.math, 2011)
        sage: union([len(union(x))==1 for x in vals])  # long time (depends on previous line)
        [True]

    AUTHOR:

    - David Joyner and William Stein
    """
    from sage.rings.all import Integer, Rational
    n = Integer(n)

    if algorithm == 'default':
        algorithm = 'pari' if n <= 30000 else 'bernmm'

    if algorithm == 'pari':
        x = pari(n).bernfrac()         # Use the PARI C library
        return Rational(x)
    elif algorithm == 'gap':
        import sage.interfaces.gap
        x = sage.interfaces.gap.gap('Bernoulli(%s)'%n)
        return Rational(x)
    elif algorithm == 'magma':
        import sage.interfaces.magma
        x = sage.interfaces.magma.magma('Bernoulli(%s)'%n)
        return Rational(x)
    elif algorithm == 'gp':
        import sage.interfaces.gp
        x = sage.interfaces.gp.gp('bernfrac(%s)'%n)
        return Rational(x)
    elif algorithm == 'bernmm':
        import sage.rings.bernmm
        return sage.rings.bernmm.bernmm_bern_rat(n, num_threads)
    else:
        raise ValueError, "invalid choice of algorithm"


def factorial(n, algorithm='gmp'):
    r"""
    Compute the factorial of `n`, which is the product
    `1\cdot 2\cdot 3 \cdots (n-1)\cdot n`.

    INPUT:

    -  ``n`` - an integer

    -  ``algorithm`` - string (default: 'gmp'):

       -  ``'gmp'`` - use the GMP C-library factorial function

       -  ``'pari'`` - use PARI's factorial function

    OUTPUT: an integer

    EXAMPLES::

        sage: from sage.rings.arith import factorial
        sage: factorial(0)
        1
        sage: factorial(4)
        24
        sage: factorial(10)
        3628800
        sage: factorial(1) == factorial(0)
        True
        sage: factorial(6) == 6*5*4*3*2
        True
        sage: factorial(1) == factorial(0)
        True
        sage: factorial(71) == 71* factorial(70)
        True
        sage: factorial(-32)
        Traceback (most recent call last):
        ...
        ValueError: factorial -- must be nonnegative

    PERFORMANCE: This discussion is valid as of April 2006. All timings
    below are on a Pentium Core Duo 2Ghz MacBook Pro running Linux with
    a 2.6.16.1 kernel.


    -  It takes less than a minute to compute the factorial of
       `10^7` using the GMP algorithm, and the factorial of
       `10^6` takes less than 4 seconds.

    -  The GMP algorithm is faster and more memory efficient than the
       PARI algorithm. E.g., PARI computes `10^7` factorial in 100
       seconds on the core duo 2Ghz.

    -  For comparison, computation in Magma `\leq` 2.12-10 of
       `n!` is best done using ``*[1..n]``. It takes
       113 seconds to compute the factorial of `10^7` and 6
       seconds to compute the factorial of `10^6`. Mathematica
       V5.2 compute the factorial of `10^7` in 136 seconds and the
       factorial of `10^6` in 7 seconds. (Mathematica is notably
       very efficient at memory usage when doing factorial
       calculations.)
    """
    if n < 0:
        raise ValueError, "factorial -- must be nonnegative"
    if algorithm == 'gmp':
        return ZZ(n).factorial()
    elif algorithm == 'pari':
        return pari.factorial(n)
    else:
        raise ValueError, 'unknown algorithm'

def is_prime(n):
    r"""
    Returns ``True`` if `n` is prime, and ``False`` otherwise.

    AUTHORS:

    - Kevin Stueve kstueve@uw.edu (2010-01-17):
      delegated calculation to ``n.is_prime()``

    INPUT:

    -  ``n`` - the object for which to determine primality

    OUTPUT:

    -  ``bool`` - ``True`` or ``False``

    EXAMPLES::

        sage: is_prime(389)
        True
        sage: is_prime(2000)
        False
        sage: is_prime(2)
        True
        sage: is_prime(-1)
        False
        sage: factor(-6)
        -1 * 2 * 3
        sage: is_prime(1)
        False
        sage: is_prime(-2)
        False

    ALGORITHM:

    Calculation is delegated to the ``n.is_prime()`` method, or in special
    cases (e.g., Python ``int``s) to ``Integer(n).is_prime()``.  If an
    ``n.is_prime()`` method is not available, it otherwise raises a
    ``TypeError``.
    """
    if isinstance(n, sage.symbolic.expression.Expression):
        try:
            n = n.pyobject()
        except TypeError:
            pass

    from sage.structure.proof.all import arithmetic
    proof = arithmetic()
    if isinstance(n, int) or isinstance(n, long):
        from sage.rings.integer import Integer
        if proof:
            return Integer(n).is_prime()
        else:
            return Integer(n).is_pseudoprime()
    try:
        if proof:
            return n.is_prime()
        else:
            return n.is_pseudoprime()
    except AttributeError:
        raise TypeError, "is_prime() is not written for this type"

def is_pseudoprime(n, flag=0):
    r"""
    Returns True if `x` is a pseudo-prime, and False otherwise.  The
    result is *NOT* proven correct - *this is a pseudo-primality
    test!*.

    INPUT:

        -  ``flag`` - int
        - ``0`` (default): checks whether x is a Baillie-Pomerance-
          Selfridge-Wagstaff pseudo prime (strong Rabin-Miller pseudo
          prime for base 2, followed by strong Lucas test for the
          sequence (P,-1), P smallest positive integer such that `P^2
          - 4` is not a square mod x).
        - ``>0``: checks whether x is a strong Miller-Rabin pseudo
          prime for flag randomly chosen bases (with end-matching to
          catch square roots of -1).

    OUTPUT:

        -  ``bool`` - True or False

    .. note::

       We do not consider negatives of prime numbers as prime.

    EXAMPLES::

        sage: is_pseudoprime(389)
        True
        sage: is_pseudoprime(2000)
        False
        sage: is_pseudoprime(2)
        True
        sage: is_pseudoprime(-1)
        False
        sage: factor(-6)
        -1 * 2 * 3
        sage: is_pseudoprime(1)
        False
        sage: is_pseudoprime(-2)
        False

    IMPLEMENTATION: Calls the PARI ispseudoprime function.
    """
    n = ZZ(n)
    return pari(n).ispseudoprime()

def is_prime_power(n, flag=0):
    r"""
    Returns True if `n` is a prime power, and False otherwise.  The
    result is proven correct - *this is NOT a pseudo-primality test!*.

    INPUT:

        -  ``n`` - an integer or rational number
        -  ``flag (for primality testing)`` - int
        - ``0`` (default): use a combination of algorithms.
        - ``1``: certify primality using the Pocklington-Lehmer Test.
        - ``2``: certify primality using the APRCL test.


    EXAMPLES::

        sage: is_prime_power(389)
        True
        sage: is_prime_power(2000)
        False
        sage: is_prime_power(2)
        True
        sage: is_prime_power(1024)
        True
        sage: is_prime_power(-1)
        False
        sage: is_prime_power(1)
        True
        sage: is_prime_power(997^100)
        True
        sage: is_prime_power(1/2197)
        True
        sage: is_prime_power(1/100)
        False
        sage: is_prime_power(2/5)
        False
    """
    try:
        n = ZZ(n)
    except TypeError:
        # n might be a nonintegral rational number, in which case it is a
        # prime power iff the integer 1/n is a prime power
        r = QQ(1/n)
        if not r.is_integral():
            return False
        n = ZZ(r)

    from sage.structure.proof.all import arithmetic
    proof = arithmetic()
    if proof:
        return n.is_prime_power(flag=flag)
    else:
        return is_pseudoprime_small_power(n)

def is_pseudoprime_small_power(n, bound=1024, get_data=False):
    r"""
    Return True if `n` is a small power of a pseudoprime, and False
    otherwise.  The result is *NOT* proven correct - *this IS a
    pseudo-primality test!*.

    If `get_data` is set to true and `n = p^d`, for a pseudoprime `p`
    and power `d`, return [(p, d)].


    INPUT:

        -  ``n`` - an integer
        -  ``bound (default: 1024)`` - int: highest power to test.
        -  ``get_data`` - boolean: return small pseudoprime and the power.

    EXAMPLES::

        sage: is_pseudoprime_small_power(389)
        True
        sage: is_pseudoprime_small_power(2000)
        False
        sage: is_pseudoprime_small_power(2)
        True
        sage: is_pseudoprime_small_power(1024)
        True
        sage: is_pseudoprime_small_power(-1)
        False
        sage: is_pseudoprime_small_power(1)
        True
        sage: is_pseudoprime_small_power(997^100)
        True

    The default bound is 1024::

        sage: is_pseudoprime_small_power(3^1024)
        True
        sage: is_pseudoprime_small_power(3^1025)
        False

    But it can be set higher or lower::

        sage: is_pseudoprime_small_power(3^1025, bound=2000)
        True
        sage: is_pseudoprime_small_power(3^100, bound=20)
        False

    Use of the get_data keyword::

        sage: is_pseudoprime_small_power(3^1024, get_data=True)
        [(3, 1024)]
        sage: is_pseudoprime_small_power(2^256, get_data=True)
        [(2, 256)]
        sage: is_pseudoprime_small_power(31, get_data=True)
        [(31, 1)]
        sage: is_pseudoprime_small_power(15, get_data=True)
        False
    """
    n = ZZ(n)
    if n == 1:
        # canonical way to write 1 as a prime power?
        return True
    if n <= 0:
        return False
    if n.is_pseudoprime():
        if get_data == True:
            return [(n, 1)]
        else:
            return True
    for i in xrange(2, bound + 1):
        p, boo = n.nth_root(i, truncate_mode=True)
        if boo:
            if p.is_pseudoprime():
                if get_data == True:
                    return [(p, i)]
                else:
                    return True
    return False


def valuation(m,*args1, **args2):
    """
    This actually just calls the m.valuation() method.
    See the documentation of m.valuation() for a more precise description.
    Use of this function by developers is discouraged. Use m.valuation() instead.

    .. NOTE::

        This is not always a valuation in the mathematical sense.
        For more information see:
        sage.rings.finite_rings.integer_mod.IntegerMod_int.valuation

    EXAMPLES::

        sage: valuation(512,2)
        9
        sage: valuation(1,2)
        0
        sage: valuation(5/9, 3)
        -2

    Valuation of 0 is defined, but valuation with respect to 0 is not::

        sage: valuation(0,7)
        +Infinity
        sage: valuation(3,0)
        Traceback (most recent call last):
        ...
        ValueError: You can only compute the valuation with respect to a integer larger than 1.

    Here are some other examples::

        sage: valuation(100,10)
        2
        sage: valuation(200,10)
        2
        sage: valuation(243,3)
        5
        sage: valuation(243*10007,3)
        5
        sage: valuation(243*10007,10007)
        1
        sage: y = QQ['y'].gen()
        sage: valuation(y^3, y)
        3
        sage: x = QQ[['x']].gen()
        sage: valuation((x^3-x^2)/(x-4))
        2
        sage: valuation(4r,2r)
        2
        sage: valuation(1r,1r)
        Traceback (most recent call last):
        ...
        ValueError: You can only compute the valuation with respect to a integer larger than 1.
    """
    if isinstance(m,(int,long)):
        m=ZZ(m)
    return m.valuation(*args1, **args2)

def prime_powers(start, stop=None):
    r"""
    List of all positive primes powers between ``start`` and
    ``stop``-1, inclusive. If the second argument is omitted, returns
    the prime powers up to the first argument.

    INPUT:

    - ``start`` - an integer. If two inputs are given, a lower bound
      for the returned set of prime powers. If this is the only input,
      then it is an upper bound.

    - ``stop`` - an integer (default: ``None``) An upper bound for the
      returned set of prime powers.

    OUTPUT:

    The set of all prime powers between ``start`` and ``stop`` or, if
    only one argument is passed, the set of all prime powers between 1
    and ``start``. Note that we will here say that the number `n` is a
    prime power if `n=p^k`, where `p` is a prime number and `k` is a
    nonnegative integer. Thus, `1` is a prime power, as `1 = 2^0`.

    EXAMPLES::

        sage: prime_powers(20)
        [1, 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19]
        sage: len(prime_powers(1000))
        194
        sage: len(prime_range(1000))
        168
        sage: a = [z for z in range(95,1234) if is_prime_power(z)]
        sage: b = prime_powers(95,1234)
        sage: len(b)
        194
        sage: len(a)
        194
        sage: a[:10]
        [97, 101, 103, 107, 109, 113, 121, 125, 127, 128]
        sage: b[:10]
        [97, 101, 103, 107, 109, 113, 121, 125, 127, 128]
        sage: a == b
        True
        sage: prime_powers(10,7)
        []
        sage: prime_powers(-5)
        []
        sage: prime_powers(-1,2)
        [1]

    TESTS::

        sage: v = prime_powers(10)
        sage: type(v[0])      # trac #922
        <type 'sage.rings.integer.Integer'>

        sage: prime_powers("foo")
        Traceback (most recent call last):
        ...
        TypeError: start must be an integer, foo is not an integer

        sage: prime_powers(6, "bar")
        Traceback (most recent call last):
        ...
        TypeError: stop must be an integer, bar is not an integer

    """
    from sage.rings.integer import Integer
    # check to ensure that both inputs are positive integers
    if not isinstance(start, (int, Integer)):
        raise TypeError("start must be an integer, {} is not an integer".format(start))
    if not (isinstance(stop, (int, Integer)) or stop is None):
        raise TypeError("stop must be an integer, {} is not an integer".format(stop))

    # coerce inputs that are ints into Integers for efficiency
    start = Integer(start)
    if(stop is not None):
        stop = Integer(stop)

    # deal with the case in which only one input is given
    if stop is None:
        start, stop = 1, Integer(start)

    # inserted to prevent an error from occurring
    if stop < 1:
        return [];

    # find all the primes in the given range
    from fast_arith import prime_range
    temp = prime_range(stop)
    output = [p for p in temp if p>=start]

    if start <= 1:
        output.append(Integer(1))

    s = stop.sqrt()
    for p in temp:
        # if p > the square root of stop, p^2 will be outside the given
        # range
        if p > s:
            break
        q = p*p
        # check if each power of p falls within the given range
        while q < stop:
            if start <= q:
                output.append(q)
            q *= p

    output.sort()
    return output

def primes_first_n(n, leave_pari=False):
    r"""
    Return the first `n` primes.

    INPUT:

    - `n` - a nonnegative integer

    OUTPUT:

    - a list of the first `n` prime numbers.

    EXAMPLES::

        sage: primes_first_n(10)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        sage: len(primes_first_n(1000))
        1000
        sage: primes_first_n(0)
        []
    """
    if n < 0:
        raise ValueError, "n must be nonnegative"
    if n < 1:
        return []
    return fast_arith.prime_range(pari.nth_prime(n) + 1)

#
# This is from
#    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/366178
# It's impressively fast given that it's in Pure Python.
#
def eratosthenes(n):
    r"""
    Return a list of the primes `\leq n`.

    This is extremely slow and is for educational purposes only.

    INPUT:

    -  ``n`` - a positive integer

    OUTPUT:

    - a list of primes less than or equal to n.


    EXAMPLES::

        sage: len(eratosthenes(100))
        25
        sage: eratosthenes(3)
        [2, 3]
    """
    n = int(n)
    if n == 2:
        return [2]
    elif n<2:
        return []
    s = range(3,n+3,2)
    mroot = n ** 0.5
    half = (n+1)/2
    i = 0
    m = 3
    while m <= mroot:
        if s[i]:
            j = (m*m-3)/2
            s[j] = 0
            while j < half:
                s[j] = 0
                j += m
        i = i+1
        m = 2*i+3
    return [ZZ(2)] + [ZZ(x) for x in s if x and x <= n]

# My old versions; not as fast as the above.
## def eratosthenes(n):
##     """
##     Returns a list of the primes up to n, computed
##     using the Sieve of Eratosthenes.
##     Input:
##         n -- a positive integer
##     Output:
##         list -- a list of the primes up to n
##     Examples:
##     sage: eratosthenes(7)
##     [2, 3, 5, 7]
##     sage: eratosthenes(45)
##     [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
##     """
##     if n <= 1: return []
##     X = [i for i in range(3,n+1) if i%2 != 0]
##     P = [2]
##     sqrt_n = sqrt(n)
##     while len(X) > 0 and X[0] <= sqrt_n:
##         p = X[0]
##         P.append(p)
##         X = [a for a in X if a%p != 0]
##     return P + X

def primes(start, stop=None, proof=None):
    r"""
    Returns an iterator over all primes between start and stop-1,
    inclusive. This is much slower than ``prime_range``, but
    potentially uses less memory.  As with :func:`next_prime`, the optional
    argument proof controls whether the numbers returned are
    guaranteed to be prime or not.

    This command is like the xrange command, except it only iterates
    over primes. In some cases it is better to use primes than
    ``prime_range``, because primes does not build a list of all primes in
    the range in memory all at once. However, it is potentially much
    slower since it simply calls the :func:`next_prime` function
    repeatedly, and :func:`next_prime` is slow.

    INPUT:

    - ``start`` - an integer - lower bound for the primes

    - ``stop`` - an integer (or infinity) optional argument -
      giving upper (open) bound for the primes

    - ``proof`` - bool or None (default: None)  If True, the function
      yields only proven primes.  If False, the function uses a
      pseudo-primality test, which is much faster for really big
      numbers but does not provide a proof of primality. If None,
      uses the global default (see :mod:`sage.structure.proof.proof`)

    OUTPUT:

    -  an iterator over primes from start to stop-1, inclusive


    EXAMPLES::

        sage: for p in primes(5,10):
        ...     print p
        ...
        5
        7
        sage: list(primes(13))
        [2, 3, 5, 7, 11]
        sage: list(primes(10000000000, 10000000100))
        [10000000019, 10000000033, 10000000061, 10000000069, 10000000097]
        sage: max(primes(10^100, 10^100+10^4, proof=False))
        10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000009631
        sage: next(p for p in primes(10^20, infinity) if is_prime(2*p+1))
        100000000000000001243


    TESTS::

        sage: for a in range(-10, 50):
        ...    for b in range(-10, 50):
        ...        assert list(primes(a,b)) == list(filter(is_prime, xrange(a,b)))
        ...
        sage: sum(primes(-10, 9973, proof=False)) == sum(filter(is_prime, range(-10, 9973)))
        True
        sage: for p in primes(10, infinity):
        ...    if p > 20: break
        ...    print p
        ...
        11
        13
        17
        19
        sage: next(p for p in primes(10,oo)) # checks alternate infinity notation
        11
    """
    from sage.rings.infinity import infinity

    start = ZZ(start)
    if stop == None:
        stop = start
        start = ZZ(2)
    elif stop != infinity:
        stop = ZZ(stop)
    n = start - 1
    while True:
        n = next_prime(n, proof)
        if n < stop:
            yield n
        else:
            return

def next_prime_power(n):
    """
    The next prime power greater than the integer n. If n is a prime
    power, then this function does not return n, but the next prime
    power after n.

    EXAMPLES::

        sage: next_prime_power(-10)
        1
        sage: is_prime_power(1)
        True
        sage: next_prime_power(0)
        1
        sage: next_prime_power(1)
        2
        sage: next_prime_power(2)
        3
        sage: next_prime_power(10)
        11
        sage: next_prime_power(7)
        8
        sage: next_prime_power(99)
        101
    """
    if n < 0:   # negatives are not prime.
        return ZZ(1)
    if n == 2:
        return ZZ(3)
    n = ZZ(n) + 1
    while not is_prime_power(n):  # pari isprime is provably correct
        n += 1
    return n

def next_probable_prime(n):
    """
    Returns the next probable prime after self, as determined by PARI.

    INPUT:


    -  ``n`` - an integer


    EXAMPLES::

        sage: next_probable_prime(-100)
        2
        sage: next_probable_prime(19)
        23
        sage: next_probable_prime(int(999999999))
        1000000007
        sage: next_probable_prime(2^768)
        1552518092300708935148979488462502555256886017116696611139052038026050952686376886330878408828646477950487730697131073206171580044114814391444287275041181139204454976020849905550265285631598444825262999193716468750892846853816058039
    """
    return ZZ(n).next_probable_prime()

def next_prime(n, proof=None):
    """
    The next prime greater than the integer n. If n is prime, then this
    function does not return n, but the next prime after n. If the
    optional argument proof is False, this function only returns a
    pseudo-prime, as defined by the PARI nextprime function. If it is
    None, uses the global default (see :mod:`sage.structure.proof.proof`)

    INPUT:


    -  ``n`` - integer

    -  ``proof`` - bool or None (default: None)


    EXAMPLES::

        sage: next_prime(-100)
        2
        sage: next_prime(1)
        2
        sage: next_prime(2)
        3
        sage: next_prime(3)
        5
        sage: next_prime(4)
        5

    Notice that the next_prime(5) is not 5 but 7.

    ::

        sage: next_prime(5)
        7
        sage: next_prime(2004)
        2011
    """
    try:
        return n.next_prime(proof)
    except AttributeError:
        return ZZ(n).next_prime(proof)

def previous_prime(n):
    """
    The largest prime < n. The result is provably correct. If n <= 1,
    this function raises a ValueError.

    EXAMPLES::

        sage: previous_prime(10)
        7
        sage: previous_prime(7)
        5
        sage: previous_prime(8)
        7
        sage: previous_prime(7)
        5
        sage: previous_prime(5)
        3
        sage: previous_prime(3)
        2
        sage: previous_prime(2)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime
        sage: previous_prime(1)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime
        sage: previous_prime(-20)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime
    """
    n = ZZ(n)-1
    if n <= 1:
        raise ValueError, "no previous prime"
    if n <= 3:
        return ZZ(n)
    if n%2 == 0:
        n -= 1
    while not is_prime(n):
        n -= 2
    return ZZ(n)

def previous_prime_power(n):
    r"""
    The largest prime power `< n`. The result is provably
    correct. If `n \leq 2`, this function returns `-x`,
    where `x` is prime power and `-x < n` and no larger
    negative of a prime power has this property.

    EXAMPLES::

        sage: previous_prime_power(2)
        1
        sage: previous_prime_power(10)
        9
        sage: previous_prime_power(7)
        5
        sage: previous_prime_power(127)
        125

    ::

        sage: previous_prime_power(0)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime power
        sage: previous_prime_power(1)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime power

    ::

        sage: n = previous_prime_power(2^16 - 1)
        sage: while is_prime(n):
        ...    n = previous_prime_power(n)
        sage: factor(n)
        251^2
    """
    n = ZZ(n)-1
    if n <= 0:
        raise ValueError, "no previous prime power"
    while not is_prime_power(n):
        n -= 1
    return n

def random_prime(n, proof=None, lbound=2):
    """
    Returns a random prime p between `lbound` and n (i.e. `lbound <= p <= n`).
    The returned prime is chosen uniformly at random from the set of prime
    numbers less than or equal to n.

    INPUT:

    -  ``n`` - an integer >= 2.

    -  ``proof`` - bool or None (default: None) If False, the function uses a
       pseudo-primality test, which is much faster for really big numbers but
       does not provide a proof of primality. If None, uses the global default
       (see :mod:`sage.structure.proof.proof`)

    - ``lbound`` - an integer >= 2
      lower bound for the chosen primes

    EXAMPLES::

        sage: random_prime(100000)
        88237
        sage: random_prime(2)
        2

    Here we generate a random prime between 100 and 200::

        sage: random_prime(200, lbound=100)
        149

    If all we care about is finding a pseudo prime, then we can pass
    in ``proof=False`` ::

        sage: random_prime(200, proof=False, lbound=100)
        149

    TESTS::

        sage: type(random_prime(2))
        <type 'sage.rings.integer.Integer'>
        sage: type(random_prime(100))
        <type 'sage.rings.integer.Integer'>
        sage: random_prime(1, lbound=-2)   #caused Sage hang #10112
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than or equal to 2
        sage: random_prime(126, lbound=114)
        Traceback (most recent call last):
        ...
        ValueError: There are no primes between 114 and 126 (inclusive)


    AUTHORS:

    - Jon Hanke (2006-08-08): with standard Stein cleanup

    - Jonathan Bober (2007-03-17)
    """
    # since we don't want current_randstate to get
    # pulled when you say "from sage.arith import *".
    from sage.misc.randstate import current_randstate
    from sage.structure.proof.proof import get_flag
    proof = get_flag(proof, "arithmetic")
    n = ZZ(n)
    if n < 2:
        raise ValueError, "n must be greater than or equal to 2"
    if n < lbound:
        raise ValueError, "n must be at least lbound: %s"%(lbound)
    elif n == 2:
        return n
    lbound = max(2, lbound)
    if lbound > 2:
        if lbound == 3 or n <= 2*lbound - 2:
        # check for Betrand's postulate (proved by Chebyshev)
            if lbound < 25 or n <= 6*lbound/5:
            # see J. Nagura, Proc. Japan Acad. 28, (1952). 177-181.
                if lbound < 2010760 or n <= 16598*lbound/16597:
                # see L. Schoenfeld, Math. Comp. 30 (1976), no. 134, 337-360.
                    if proof:
                        smallest_prime = ZZ(lbound-1).next_prime()
                    else:
                        smallest_prime = ZZ(lbound-1).next_probable_prime()
                    if smallest_prime > n:
                        raise ValueError, \
                              "There are no primes between %s and %s (inclusive)" % (lbound, n)

    if proof:
        prime_test = is_prime
    else:
        prime_test = is_pseudoprime
    randint = current_randstate().python_random().randint
    while True:
        # In order to ensure that the returned prime is chosen
        # uniformly from the set of primes it is necessary to
        # choose a random number and then test for primality.
        # The method of choosing a random number and then returning
        # the closest prime smaller than it would typically not,
        # for example, return the first of a pair of twin primes.
        p = randint(lbound, n)
        if prime_test(p):
            return ZZ(p)


def divisors(n):
    """
    Returns a list of all positive integer divisors of the nonzero
    integer n.

    INPUT:


    -  ``n`` - the element


    EXAMPLES::

        sage: divisors(-3)
        [1, 3]
        sage: divisors(6)
        [1, 2, 3, 6]
        sage: divisors(28)
        [1, 2, 4, 7, 14, 28]
        sage: divisors(2^5)
        [1, 2, 4, 8, 16, 32]
        sage: divisors(100)
        [1, 2, 4, 5, 10, 20, 25, 50, 100]
        sage: divisors(1)
        [1]
        sage: divisors(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be nonzero
        sage: divisors(2^3 * 3^2 * 17)
        [1, 2, 3, 4, 6, 8, 9, 12, 17, 18, 24, 34, 36, 51, 68, 72, 102, 136, 153, 204, 306, 408, 612, 1224]

    This function works whenever one has unique factorization::

        sage: K.<a> = QuadraticField(7)
        sage: divisors(K.ideal(7))
        [Fractional ideal (1), Fractional ideal (a), Fractional ideal (7)]
        sage: divisors(K.ideal(3))
        [Fractional ideal (1), Fractional ideal (3), Fractional ideal (-a + 2), Fractional ideal (-a - 2)]
        sage: divisors(K.ideal(35))
        [Fractional ideal (1), Fractional ideal (35), Fractional ideal (5*a), Fractional ideal (5), Fractional ideal (a), Fractional ideal (7)]

    TESTS::

        sage: divisors(int(300))
        [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 25, 30, 50, 60, 75, 100, 150, 300]
    """
    if not n:
        raise ValueError, "n must be nonzero"

    R = parent(n)
    if R in [int, long]:
        n = ZZ(n) # we have specialized code for this case, make sure it gets used
    try:
        return n.divisors()
    except AttributeError:
        pass
    f = factor(n)
    one = R(1)
    all = [one]
    for p, e in f:
        prev = all[:]
        pn = one
        for i in range(e):
            pn *= p
            all.extend([a*pn for a in prev])
    all.sort()
    return all

class Sigma:
    """
    Return the sum of the k-th powers of the divisors of n.

    INPUT:


    -  ``n`` - integer

    -  ``k`` - integer (default: 1)


    OUTPUT: integer

    EXAMPLES::

        sage: sigma(5)
        6
        sage: sigma(5,2)
        26

    The sigma function also has a special plotting method.

    ::

        sage: P = plot(sigma, 1, 100)

    This method also works with k-th powers.

    ::

        sage: P = plot(sigma, 1, 100, k=2)

    AUTHORS:

    - William Stein: original implementation

    - Craig Citro (2007-06-01): rewrote for huge speedup

    TESTS::

        sage: sigma(100,4)
        106811523
        sage: sigma(factorial(100),3).mod(144169)
        3672
        sage: sigma(factorial(150),12).mod(691)
        176
        sage: RR(sigma(factorial(133),20))
        2.80414775675747e4523
        sage: sigma(factorial(100),0)
        39001250856960000
        sage: sigma(factorial(41),1)
        229199532273029988767733858700732906511758707916800
    """
    def __repr__(self):
        """
        A description of this class, which computes the sum of the
        k-th powers of the divisors of n.

        EXAMPLES::

            sage: Sigma().__repr__()
            'Function that adds up (k-th powers of) the divisors of n'
        """
        return "Function that adds up (k-th powers of) the divisors of n"

    def __call__(self, n, k=1):
        """
        Computes the sum of (the k-th powers of) the divisors of n.

        EXAMPLES::

            sage: q = Sigma()
            sage: q(10)
            18
            sage: q(10,2)
            130
        """
        n = ZZ(n)
        k = ZZ(k)
        one = ZZ(1)

        if (k == ZZ(0)):
            return prod([ expt+one for p, expt in factor(n) ])
        elif (k == one):
            return prod([ (p**(expt+one) - one).divide_knowing_divisible_by(p - one)
                          for p, expt in factor(n) ])
        else:
            return prod([ (p**((expt+one)*k)-one).divide_knowing_divisible_by(p**k-one)
                          for p,expt in factor(n) ])

    def plot(self, xmin=1, xmax=50, k=1, pointsize=30, rgbcolor=(0,0,1), join=True,
             **kwds):
        """
        Plot the sigma (sum of k-th powers of divisors) function.

        INPUT:


        -  ``xmin`` - default: 1

        -  ``xmax`` - default: 50

        -  ``k`` - default: 1

        -  ``pointsize`` - default: 30

        -  ``rgbcolor`` - default: (0,0,1)

        -  ``join`` - default: True; whether to join the
           points.

        -  ``**kwds`` - passed on

        EXAMPLES::

            sage: p = Sigma().plot()
            sage: p.ymax()
            124.0
        """
        v = [(n,sigma(n,k)) for n in range(xmin,xmax + 1)]
        from sage.plot.all import list_plot
        P = list_plot(v, pointsize=pointsize, rgbcolor=rgbcolor, **kwds)
        if join:
            P += list_plot(v, plotjoined=True, rgbcolor=(0.7,0.7,0.7), **kwds)
        return P

sigma = Sigma()

def gcd(a, b=None, **kwargs):
    r"""
    The greatest common divisor of a and b, or if a is a list and b is
    omitted the greatest common divisor of all elements of a.

    INPUT:


    -  ``a,b`` - two elements of a ring with gcd or

    -  ``a`` - a list or tuple of elements of a ring with
       gcd


    Additional keyword arguments are passed to the respectively called
    methods.

    OUTPUT:

    The given elements are first coerced into a common parent. Then,
    their greatest common divisor *in that common parent* is returned.

    EXAMPLES::

        sage: GCD(97,100)
        1
        sage: GCD(97*10^15, 19^20*97^2)
        97
        sage: GCD(2/3, 4/5)
        2/15
        sage: GCD([2,4,6,8])
        2
        sage: GCD(srange(0,10000,10))  # fast  !!
        10

    Note that to take the gcd of `n` elements for `n \not= 2` you must
    put the elements into a list by enclosing them in ``[..]``.  Before
    #4988 the following wrongly returned 3 since the third parameter
    was just ignored::

        sage: gcd(3,6,2)
        Traceback (most recent call last):
        ...
        TypeError: gcd() takes at most 2 arguments (3 given)
        sage: gcd([3,6,2])
        1

    Similarly, giving just one element (which is not a list) gives an error::

        sage: gcd(3)
        Traceback (most recent call last):
        ...
        TypeError: 'sage.rings.integer.Integer' object is not iterable

    By convention, the gcd of the empty list is (the integer) 0::

        sage: gcd([])
        0
        sage: type(gcd([]))
        <type 'sage.rings.integer.Integer'>

    TESTS:

    The following shows that indeed coercion takes place before computing
    the gcd. This behaviour was introduced in trac ticket #10771::

        sage: R.<x>=QQ[]
        sage: S.<x>=ZZ[]
        sage: p = S.random_element()
        sage: q = R.random_element()
        sage: parent(gcd(1/p,q))
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: parent(gcd([1/p,q]))
        Fraction Field of Univariate Polynomial Ring in x over Rational Field

    Make sure we try QQ and not merely ZZ (:trac:`13014`)::

        sage: bool(gcd(2/5, 3/7) == gcd(SR(2/5), SR(3/7)))
        True

    Make sure that the gcd of Expressions stays symbolic::

        sage: parent(gcd(2, 4))
        Integer Ring
        sage: parent(gcd(SR(2), 4))
        Symbolic Ring
        sage: parent(gcd(2, SR(4)))
        Symbolic Ring
        sage: parent(gcd(SR(2), SR(4)))
        Symbolic Ring

    Verify that objects without gcd methods but which can't be
    coerced to ZZ or QQ raise an error::

        sage: F.<a,b> = FreeMonoid(2)
        sage: gcd(a,b)
        Traceback (most recent call last):
        ...
        TypeError: unable to find gcd

    """
    # Most common use case first:
    if b is not None:
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        a,b = cm.canonical_coercion(a,b)
        try:
            GCD = a.gcd
        except AttributeError:
            try:
                return ZZ(a).gcd(ZZ(b))
            except TypeError:
                raise TypeError("unable to find gcd")
        return GCD(b, **kwargs)

    from sage.structure.sequence import Sequence
    seq = Sequence(a)
    U = seq.universe()
    if U is ZZ or U is int or U is long:# ZZ.has_coerce_map_from(U):
        return sage.rings.integer.GCD_list(a)
    return __GCD_sequence(seq, **kwargs)

GCD = gcd

def __GCD_sequence(v, **kwargs):
    """
    Internal function returning the gcd of the elements of a sequence

    INPUT:


    -  ``v`` - A sequence (possibly empty)


    OUTPUT: The gcd of the elements of the sequence as an element of
    the sequence's universe, or the integer 0 if the sequence is
    empty.

    EXAMPLES::

        sage: from sage.rings.arith import __GCD_sequence
        sage: from sage.structure.sequence import Sequence
        sage: l = ()
        sage: __GCD_sequence(l)
        0
        sage: __GCD_sequence(Sequence(srange(10)))
        1
        sage: X=polygen(QQ)
        sage: __GCD_sequence(Sequence((2*X+4,2*X^2,2)))
        1
        sage: X=polygen(ZZ)
        sage: __GCD_sequence(Sequence((2*X+4,2*X^2,2)))
        2
    """
    if len(v) == 0:
        return ZZ(0)
    if hasattr(v,'universe'):
        g = v.universe()(0)
    else:
        g = ZZ(0)
    one = v.universe()(1)
    for vi in v:
        g = vi.gcd(g, **kwargs)
        if g == one:
            return g
    return g

def lcm(a, b=None):
    """
    The least common multiple of a and b, or if a is a list and b is
    omitted the least common multiple of all elements of a.

    Note that LCM is an alias for lcm.

    INPUT:


    -  ``a,b`` - two elements of a ring with lcm or

    -  ``a`` - a list or tuple of elements of a ring with
       lcm

    OUTPUT:

    First, the given elements are coerced into a common parent. Then,
    their least common multiple *in that parent* is returned.

    EXAMPLES::

        sage: lcm(97,100)
        9700
        sage: LCM(97,100)
        9700
        sage: LCM(0,2)
        0
        sage: LCM(-3,-5)
        15
        sage: LCM([1,2,3,4,5])
        60
        sage: v = LCM(range(1,10000))   # *very* fast!
        sage: len(str(v))
        4349


    TESTS:

    The following tests against a bug that was fixed in trac
    ticket #10771::

        sage: lcm(4/1,2)
        4

    The following shows that indeed coercion takes place before
    computing the least common multiple::

        sage: R.<x>=QQ[]
        sage: S.<x>=ZZ[]
        sage: p = S.random_element()
        sage: q = R.random_element()
        sage: parent(lcm([1/p,q]))
        Fraction Field of Univariate Polynomial Ring in x over Rational Field

    Make sure we try QQ and not merely ZZ (:trac:`13014`)::

        sage: bool(lcm(2/5, 3/7) == lcm(SR(2/5), SR(3/7)))
        True

    Make sure that the lcm of Expressions stays symbolic::

        sage: parent(lcm(2, 4))
        Integer Ring
        sage: parent(lcm(SR(2), 4))
        Symbolic Ring
        sage: parent(lcm(2, SR(4)))
        Symbolic Ring
        sage: parent(lcm(SR(2), SR(4)))
        Symbolic Ring

    Verify that objects without lcm methods but which can't be
    coerced to ZZ or QQ raise an error::

        sage: F.<a,b> = FreeMonoid(2)
        sage: lcm(a,b)
        Traceback (most recent call last):
        ...
        TypeError: unable to find lcm

    """
    # Most common use case first:
    if b is not None:
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        a,b = cm.canonical_coercion(a,b)
        try:
            LCM = a.lcm
        except AttributeError:
            try:
                return ZZ(a).lcm(ZZ(b))
            except TypeError:
                raise TypeError("unable to find lcm")
        return LCM(b)

    from sage.structure.sequence import Sequence
    seq = Sequence(a)
    U = seq.universe()
    if U is ZZ or U is int or U is long:
        return sage.rings.integer.LCM_list(a)
    return __LCM_sequence(seq)

LCM = lcm

def __LCM_sequence(v):
    """
    Internal function returning the lcm of the elements of a sequence

    INPUT:


    -  ``v`` - A sequence (possibly empty)


    OUTPUT: The lcm of the elements of the sequence as an element of
    the sequence's universe, or the integer 1 if the sequence is
    empty.

    EXAMPLES::

        sage: from sage.structure.sequence import Sequence
        sage: from sage.rings.arith import __LCM_sequence
        sage: l = Sequence(())
        sage: __LCM_sequence(l)
        1

    This is because lcm(0,x)=0 for all x (by convention)

    ::

        sage: __LCM_sequence(Sequence(srange(100)))
        0

    So for the lcm of all integers up to 10 you must do this::

        sage: __LCM_sequence(Sequence(srange(1,100)))
        69720375229712477164533808935312303556800

    Note that the following example did not work in QQ[] as of 2.11,
    but does in 3.1.4; the answer is different, though equivalent::

        sage: R.<X>=ZZ[]
        sage: __LCM_sequence(Sequence((2*X+4,2*X^2,2)))
        2*X^3 + 4*X^2
        sage: R.<X>=QQ[]
        sage: __LCM_sequence(Sequence((2*X+4,2*X^2,2)))
        X^3 + 2*X^2
    """
    if len(v) == 0:
        return ZZ(1)
    try:
        g = v.universe()(1)
    except AttributeError:
        g = ZZ(1)
    for vi in v:
        g = vi.lcm(g)
        if not g:
            return g
    return g

def xlcm(m,n):
    r"""
    Extended lcm function: given two positive integers `m,n`, returns
    a triple `(l,m_1,n_1)` such that `l=\mathop{\mathrm{lcm}}(m,n)=m_1
    \cdot n_1` where `m_1|m`, `n_1|n` and `\gcd(m_1,n_1)=1`, all with no
    factorization.

    Used to construct an element of order `l` from elements of orders `m,n`
    in any group: see sage/groups/generic.py for examples.

    EXAMPLES::

        sage: xlcm(120,36)
        (360, 40, 9)
    """
    g=gcd(m,n)
    l=m*n//g      # = lcm(m,n)
    g=gcd(m,n//g) # divisible by those primes which divide n to a
                  # higher power than m

    while not g==1:
        m//=g
        g=gcd(m,g)

    n=l//m;
    return (l,m,n)

def xgcd(a, b):
    r"""
    Return a triple ``(g,s,t)`` such that `g = s\cdot a+t\cdot b = \gcd(a,b)`.

    .. note::

       One exception is if `a` and `b` are not in a PID, e.g., they are
       both polynomials over the integers, then this function can't in
       general return ``(g,s,t)`` as above, since they need not exist.
       Instead, over the integers, we first multiply `g` by a divisor of
       the resultant of `a/g` and `b/g`, up to sign.

    INPUT:


    -  ``a, b`` - integers or univariate polynomials (or
       any type with an xgcd method).


    OUTPUT:

    -  ``g, s, t`` - such that `g = s\cdot a + t\cdot b`


    .. note::

       There is no guarantee that the returned cofactors (s and t) are
       minimal. In the integer case, see
       :meth:`sage.rings.integer.Integer._xgcd()` for minimal
       cofactors.

    EXAMPLES::

        sage: xgcd(56, 44)
        (4, 4, -5)
        sage: 4*56 + (-5)*44
        4
        sage: g, a, b = xgcd(5/1, 7/1); g, a, b
        (1, 3, -2)
        sage: a*(5/1) + b*(7/1) == g
        True
        sage: x = polygen(QQ)
        sage: xgcd(x^3 - 1, x^2 - 1)
        (x - 1, 1, -x)
        sage: K.<g> = NumberField(x^2-3)
        sage: R.<a,b> = K[]
        sage: S.<y> = R.fraction_field()[]
        sage: xgcd(y^2, a*y+b)
        (1, a^2/b^2, ((-a)/b^2)*y + 1/b)
        sage: xgcd((b+g)*y^2, (a+g)*y+b)
        (1, (a^2 + (2*g)*a + 3)/(b^3 + (g)*b^2), ((-a + (-g))/b^2)*y + 1/b)

    We compute an xgcd over the integers, where the linear combination
    is not the gcd but the resultant::

        sage: R.<x> = ZZ[]
        sage: gcd(2*x*(x-1), x^2)
        x
        sage: xgcd(2*x*(x-1), x^2)
        (2*x, -1, 2)
        sage: (2*(x-1)).resultant(x)
        2
    """
    try:
        return a.xgcd(b)
    except AttributeError:
        pass
    if not isinstance(a, sage.rings.integer.Integer):
        a = ZZ(a)
    return a.xgcd(ZZ(b))

XGCD = xgcd

## def XGCD_python(a, b):
##     """
##     Returns triple (g,p,q) such that g = p*a+b*q = GCD(a,b).
##     This function should behave exactly the same as XGCD,
##     but is implemented in pure python.
##     """
##     if a == 0 and b == 0:
##         return (0,0,1)
##     if a == 0:
##         return (abs(b), 0, b/abs(b))
##     if b == 0:
##         return (abs(a), a/abs(a), 0)
##     psign = 1
##     qsign = 1
##     if a < 0:
##         a = -a
##         psign = -1
##     if b < 0:
##         b = -b
##         qsign = -1
##     p = 1; q = 0; r = 0; s = 1
##     while b != 0:
##         c = a % b
##         quot = a/b
##         a = b; b = c
##         new_r = p - quot*r
##         new_s = q - quot*s
##         p = r; q = s
##         r = new_r; s = new_s
##     return (a, p*psign, q*qsign)

def inverse_mod(a, m):
    """
    The inverse of the ring element a modulo m.

    If no special inverse_mod is defined for the elements, it tries to
    coerce them into integers and perform the inversion there

    ::

        sage: inverse_mod(7,1)
        0
        sage: inverse_mod(5,14)
        3
        sage: inverse_mod(3,-5)
        2
    """
    try:
        return a.inverse_mod(m)
    except AttributeError:
        return integer.Integer(a).inverse_mod(m)

#######################################################
# Functions to find the fastest available commands
# for gcd and inverse_mod
#######################################################

def get_gcd(order):
    """
    Return the fastest gcd function for integers of size no larger than
    order.

    EXAMPLES::

        sage: sage.rings.arith.get_gcd(4000)
        <built-in method gcd_int of sage.rings.fast_arith.arith_int object at ...>
        sage: sage.rings.arith.get_gcd(400000)
        <built-in method gcd_longlong of sage.rings.fast_arith.arith_llong object at ...>
        sage: sage.rings.arith.get_gcd(4000000000)
        <function gcd at ...>
    """
    if order <= 46340:   # todo: don't hard code
        return fast_arith.arith_int().gcd_int
    elif order <= 2147483647:   # todo: don't hard code
        return fast_arith.arith_llong().gcd_longlong
    else:
        return gcd

def get_inverse_mod(order):
    """
    Return the fastest inverse_mod function for integers of size no
    larger than order.

    EXAMPLES::

        sage: sage.rings.arith.get_inverse_mod(6000)
        <built-in method inverse_mod_int of sage.rings.fast_arith.arith_int object at ...>
        sage: sage.rings.arith.get_inverse_mod(600000)
        <built-in method inverse_mod_longlong of sage.rings.fast_arith.arith_llong object at ...>
        sage: sage.rings.arith.get_inverse_mod(6000000000)
        <function inverse_mod at ...>
    """
    if order <= 46340:   # todo: don't hard code
        return fast_arith.arith_int().inverse_mod_int
    elif order <= 2147483647:   # todo: don't hard code
        return fast_arith.arith_llong().inverse_mod_longlong
    else:
        return inverse_mod

# def sqrt_mod(a, m):
#     """A square root of a modulo m."""

# def xxx_inverse_mod(a, m):
#     """The inverse of a modulo m."""
#     g,s,t = XGCD(a,m)
#     if g != 1:
#         raise "inverse_mod(a=%s,m=%s), error since GCD=%s"%(a,m,g)
#     return s

def power_mod(a,n,m):
    """
    The n-th power of a modulo the integer m.

    EXAMPLES::

        sage: power_mod(0,0,5)
        Traceback (most recent call last):
        ...
        ArithmeticError: 0^0 is undefined.
        sage: power_mod(2,390,391)
        285
        sage: power_mod(2,-1,7)
        4
        sage: power_mod(11,1,7)
        4
        sage: R.<x> = ZZ[]
        sage: power_mod(3*x, 10, 7)
        4*x^10

        sage: power_mod(11,1,0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: modulus must be nonzero.
    """
    if m==0:
        raise ZeroDivisionError, "modulus must be nonzero."
    if m==1:
        return 0
    if n < 0:
        ainv = inverse_mod(a,m)
        return power_mod(ainv, -n, m)
    if n==0:
        if a == 0:
            raise ArithmeticError, "0^0 is undefined."
        return 1

    apow = a % m
    while n&1 == 0:
        apow = (apow*apow) % m
        n = n >> 1
    power = apow
    n = n >> 1
    while n != 0:
        apow = (apow*apow) % m
        if n&1 != 0:
            power = (power*apow) % m
        n = n >> 1

    return power


def rational_reconstruction(a, m, algorithm='fast'):
    r"""
    This function tries to compute `x/y`, where `x/y` is a rational number in
    lowest terms such that the reduction of `x/y` modulo `m` is equal to `a` and
    the absolute values of `x` and `y` are both `\le \sqrt{m/2}`. If such `x/y`
    exists, that pair is unique and this function returns it. If no
    such pair exists, this function raises ZeroDivisionError.

    An efficient algorithm for computing rational reconstruction is
    very similar to the extended Euclidean algorithm. For more details,
    see Knuth, Vol 2, 3rd ed, pages 656-657.

    INPUT:

    - ``a`` - an integer

    - ``m`` - a modulus

    - ``algorithm`` - (default: 'fast')

      - ``'fast'`` - a fast compiled implementation

      - ``'python'`` - a slow pure python implementation

    OUTPUT:

    Numerator and denominator `n`, `d` of the unique rational number
    `r=n/d`, if it exists, with `n` and `|d| \le \sqrt{N/2}`. Return
    `(0,0)` if no such number exists.

    The algorithm for rational reconstruction is described (with a
    complete nontrivial proof) on pages 656-657 of Knuth, Vol 2, 3rd
    ed. as the solution to exercise 51 on page 379. See in particular
    the conclusion paragraph right in the middle of page 657, which
    describes the algorithm thus:

        This discussion proves that the problem can be solved
        efficiently by applying Algorithm 4.5.2X with `u=m` and `v=a`,
        but with the following replacement for step X2: If
        `v3 \le \sqrt{m/2}`, the algorithm terminates. The pair
        `(x,y)=(|v2|,v3*\mathrm{sign}(v2))` is then the unique
        solution, provided that `x` and `y` are coprime and
        `x \le \sqrt{m/2}`; otherwise there is no solution. (Alg 4.5.2X is
        the extended Euclidean algorithm.)

    Knuth remarks that this algorithm is due to Wang, Kornerup, and
    Gregory from around 1983.

    EXAMPLES::

        sage: m = 100000
        sage: (119*inverse_mod(53,m))%m
        11323
        sage: rational_reconstruction(11323,m)
        119/53

    ::

        sage: rational_reconstruction(400,1000)
        Traceback (most recent call last):
        ...
        ValueError: Rational reconstruction of 400 (mod 1000) does not exist.

    ::

        sage: rational_reconstruction(3,292393, algorithm='python')
        3
        sage: a = Integers(292393)(45/97); a
        204977
        sage: rational_reconstruction(a,292393, algorithm='python')
        45/97
        sage: a = Integers(292393)(45/97); a
        204977
        sage: rational_reconstruction(a,292393, algorithm='fast')
        45/97
        sage: rational_reconstruction(293048,292393, algorithm='fast')
        Traceback (most recent call last):
        ...
        ValueError: Rational reconstruction of 655 (mod 292393) does not exist.
        sage: rational_reconstruction(293048,292393, algorithm='python')
        Traceback (most recent call last):
        ...
        ValueError: Rational reconstruction of 655 (mod 292393) does not exist.

    TESTS:

    Check that ticket #9345 is fixed::

        sage: rational_reconstruction(1, 0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero
        sage: rational_reconstruction(randint(-10^6, 10^6), 0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero
    """
    if algorithm == 'fast':
        return ZZ(a).rational_reconstruction(m)
    elif algorithm == 'python':
        return _rational_reconstruction_python(a,m)
    else:
        raise ValueError("unknown algorithm")

def _rational_reconstruction_python(a,m):
    """
    Internal fallback function for rational_reconstruction; see
    documentation of that function for details.

    INPUT:

    - ``a`` - an integer

    - ``m`` - a modulus

    EXAMPLES::

        sage: from sage.rings.arith import _rational_reconstruction_python
        sage: _rational_reconstruction_python(20,31)
        -2/3
        sage: _rational_reconstruction_python(11323,100000)
        119/53
    """
    a = int(a); m = int(m)
    a %= m
    if a == 0 or m==0:
        return ZZ(0)/ZZ(1)
    if m < 0:
        m = -m
    if a < 0:
        a = m-a
    if a == 1:
        return ZZ(1)/ZZ(1)
    u = m
    v = a
    bnd = math.sqrt(m/2)
    U = (1,0,u)
    V = (0,1,v)
    while abs(V[2]) > bnd:
        q = U[2]/V[2]  # floor is implicit
        T = (U[0]-q*V[0], U[1]-q*V[1], U[2]-q*V[2])
        U = V
        V = T
    x = abs(V[1])
    y = V[2]
    if V[1] < 0:
        y *= -1
    if x <= bnd and GCD(x,y) == 1:
        return ZZ(y) / ZZ(x)
    raise ValueError, "Rational reconstruction of %s (mod %s) does not exist."%(a,m)

def mqrr_rational_reconstruction(u, m, T):
    r"""
    Maximal Quotient Rational Reconstruction.

    For research purposes only - this is pure Python, so slow.

    INPUT:

    - ``u, m, T`` -  integers such that `m > u \ge 0`, `T > 0`.

    OUTPUT:

    Either integers `n,d` such that `d>0`, `\mathop{\mathrm{gcd}}(n,d)=1`, `n/d=u \bmod m`, and
    `T \cdot d \cdot |n| < m`, or ``None``.

    Reference: Monagan, Maximal Quotient Rational Reconstruction: An
    Almost Optimal Algorithm for Rational Reconstruction (page 11)

    This algorithm is probabilistic.

    EXAMPLES::

        sage: mqrr_rational_reconstruction(21,3100,13)
        (21, 1)
    """
    if u == 0:
        if m > T:
            return (0,1)
        else:
            return None
    n, d = 0, 0
    t0, r0 = 0, m
    t1, r1 = 1, u
    while r1 != 0 and r0 > T:
        q = r0/r1   # C division implicit floor
        if q > T:
            n, d, T = r1, t1, q
        r0, r1 = r1, r0 - q*r1
        t0, t1 = t1, t0 - q*t1
    if d != 0 and GCD(n,d) == 1:
        return (n,d)
    return None


######################


def trial_division(n, bound=None):
    """
    Return the smallest prime divisor <= bound of the positive integer
    n, or n if there is no such prime. If the optional argument bound
    is omitted, then bound <= n.

    INPUT:

    -  ``n`` - a positive integer

    - ``bound`` - (optional) a positive integer

    OUTPUT:

    -  ``int`` - a prime p=bound that divides n, or n if
       there is no such prime.


    EXAMPLES::

        sage: trial_division(15)
        3
        sage: trial_division(91)
        7
        sage: trial_division(11)
        11
        sage: trial_division(387833, 300)
        387833
        sage: # 300 is not big enough to split off a
        sage: # factor, but 400 is.
        sage: trial_division(387833, 400)
        389
    """
    if bound is None:
        return ZZ(n).trial_division()
    else:
        return ZZ(n).trial_division(bound)

def factor(n, proof=None, int_=False, algorithm='pari', verbose=0, **kwds):
    """
    Returns the factorization of ``n``.  The result depends on the
    type of ``n``.

    If ``n`` is an integer, returns the factorization as an object
    of type ``Factorization``.

    If n is not an integer, ``n.factor(proof=proof, **kwds)`` gets called.
    See ``n.factor??`` for more documentation in this case.

    .. warning::

       This means that applying ``factor`` to an integer result of
       a symbolic computation will not factor the integer, because it is
       considered as an element of a larger symbolic ring.

       EXAMPLE::

           sage: f(n)=n^2
           sage: is_prime(f(3))
           False
           sage: factor(f(3))
           9

    INPUT:

    -  ``n`` - an nonzero integer

    -  ``proof`` - bool or None (default: None)

    -  ``int_`` - bool (default: False) whether to return
       answers as Python ints

    -  ``algorithm`` - string

       - ``'pari'`` - (default) use the PARI c library

       - ``'kash'`` - use KASH computer algebra system (requires the
         optional kash package be installed)

       - ``'magma'`` - use Magma (requires magma be installed)

    -  ``verbose`` - integer (default: 0); PARI's debug
       variable is set to this; e.g., set to 4 or 8 to see lots of output
       during factorization.

    OUTPUT:

    -  factorization of n

    The qsieve and ecm commands give access to highly optimized
    implementations of algorithms for doing certain integer
    factorization problems. These implementations are not used by the
    generic factor command, which currently just calls PARI (note that
    PARI also implements sieve and ecm algorithms, but they aren't as
    optimized). Thus you might consider using them instead for certain
    numbers.

    The factorization returned is an element of the class
    :class:`~sage.structure.factorization.Factorization`; see Factorization??
    for more details, and examples below for usage. A Factorization contains
    both the unit factor (+1 or -1) and a sorted list of (prime, exponent)
    pairs.

    The factorization displays in pretty-print format but it is easy to
    obtain access to the (prime,exponent) pairs and the unit, to
    recover the number from its factorization, and even to multiply two
    factorizations. See examples below.

    EXAMPLES::

        sage: factor(500)
        2^2 * 5^3
        sage: factor(-20)
        -1 * 2^2 * 5
        sage: f=factor(-20)
        sage: list(f)
        [(2, 2), (5, 1)]
        sage: f.unit()
        -1
        sage: f.value()
        -20
        sage: factor( -next_prime(10^2) * next_prime(10^7) )
        -1 * 101 * 10000019

    ::

        sage: factor(-500, algorithm='kash')      # optional - kash
        -1 * 2^2 * 5^3

    ::

        sage: factor(-500, algorithm='magma')     # optional - magma
        -1 * 2^2 * 5^3

    ::

        sage: factor(0)
        Traceback (most recent call last):
        ...
        ArithmeticError: Prime factorization of 0 not defined.
        sage: factor(1)
        1
        sage: factor(-1)
        -1
        sage: factor(2^(2^7)+1)
        59649589127497217 * 5704689200685129054721

    Sage calls PARI's factor, which has proof False by default.
    Sage has a global proof flag, set to True by default (see
    :mod:`sage.structure.proof.proof`, or proof.[tab]). To override
    the default, call this function with proof=False.

    ::

        sage: factor(3^89-1, proof=False)
        2 * 179 * 1611479891519807 * 5042939439565996049162197

    ::

        sage: factor(2^197 + 1)  # long time (2s)
        3 * 197002597249 * 1348959352853811313 * 251951573867253012259144010843

    Any object which has a factor method can be factored like this::

        sage: K.<i> = QuadraticField(-1)
        sage: factor(122 - 454*i)
        (-1) * (-i - 4) * (-3*i - 2) * (-i - 2)^3 * (i + 1)^3

    To access the data in a factorization::

        sage: f = factor(420); f
        2^2 * 3 * 5 * 7
        sage: [x for x in f]
        [(2, 2), (3, 1), (5, 1), (7, 1)]
        sage: [p for p,e in f]
        [2, 3, 5, 7]
        sage: [e for p,e in f]
        [2, 1, 1, 1]
        sage: [p^e for p,e in f]
        [4, 3, 5, 7]

    """
    if isinstance(n, (int, long)):
        n = ZZ(n)

    if isinstance(n, integer.Integer):
        return n.factor(proof=proof, algorithm=algorithm,
                        int_ = int_, verbose=verbose)
    else:
        # e.g. n = x**2 + y**2 + 2*x*y
        try:
            return n.factor(proof=proof, **kwds)
        except AttributeError:
            raise TypeError, "unable to factor n"
        except TypeError:
            # Just in case factor method doesn't have a proof option.
            try:
                return n.factor(**kwds)
            except AttributeError:
                raise TypeError, "unable to factor n"

def radical(n, *args, **kwds):
    """
    Return the product of the prime divisors of n.

    This calls ``n.radical(*args, **kwds)``.  If that doesn't work, it
    does ``n.factor(*args, **kwds)`` and returns the product of the prime
    factors in the resulting factorization.

    EXAMPLES::

        sage: radical(2 * 3^2 * 5^5)
        30
        sage: radical(0)
        Traceback (most recent call last):
        ...
        ArithmeticError: Radical of 0 not defined.
        sage: K.<i> = QuadraticField(-1)
        sage: radical(K(2))
        i + 1

    The next example shows how to compute the radical of a number,
    assuming no prime > 100000 has exponent > 1 in the factorization::

        sage: n = 2^1000-1; n / radical(n, limit=100000)
        125
    """
    try:
        return n.radical(*args, **kwds)
    except AttributeError:
        return n.factor(*args, **kwds).radical_value()

def prime_divisors(n):
    """
    The prime divisors of the integer n, sorted in increasing order. If
    n is negative, we do *not* include -1 among the prime divisors,
    since -1 is not a prime number.

    EXAMPLES::

        sage: prime_divisors(1)
        []
        sage: prime_divisors(100)
        [2, 5]
        sage: prime_divisors(-100)
        [2, 5]
        sage: prime_divisors(2004)
        [2, 3, 167]
    """
    return [p for p,_ in factor(n) if p != -1]

prime_factors = prime_divisors

def odd_part(n):
    r"""
    The odd part of the integer `n`. This is `n / 2^v`,
    where `v = \mathrm{valuation}(n,2)`.

    EXAMPLES::

        sage: odd_part(5)
        5
        sage: odd_part(4)
        1
        sage: odd_part(factorial(31))
        122529844256906551386796875
    """
    if not isinstance(n, integer.Integer):
        n = ZZ(n)
    return n.odd_part()

def prime_to_m_part(n,m):
    """
    Returns the prime-to-m part of n, i.e., the largest divisor of n
    that is coprime to m.

    INPUT:


    -  ``n`` - Integer (nonzero)

    -  ``m`` - Integer


    OUTPUT: Integer

    EXAMPLES::

        sage: z = 43434
        sage: z.prime_to_m_part(20)
        21717
    """
    if n == 0:
        raise ValueError, "n must be nonzero."
    if m == 0:
        return ZZ(1)
    n = ZZ(n)
    m = ZZ(m)
    while True:
        g = gcd(n,m)
        if g == 1:
            return n
        n = n // g


def is_square(n, root=False):
    """
    Returns whether or not n is square, and if n is a square also
    returns the square root. If n is not square, also returns None.

    INPUT:


    -  ``n`` - an integer

    -  ``root`` - whether or not to also return a square
       root (default: False)


    OUTPUT:


    -  ``bool`` - whether or not a square

    -  ``object`` - (optional) an actual square if found,
       and None otherwise.


    EXAMPLES::

        sage: is_square(2)
        False
        sage: is_square(4)
        True
        sage: is_square(2.2)
        True
        sage: is_square(-2.2)
        False
        sage: is_square(CDF(-2.2))
        True
        sage: is_square((x-1)^2)
        True

    ::

        sage: is_square(4, True)
        (True, 2)
    """
    try:
        if root:
            try:
                return n.is_square(root)
            except TypeError:
                if n.is_square():
                    return True, n.sqrt()
                else:
                    return False, None
        return n.is_square()
    except (AttributeError, NotImplementedError):
        pass
    t, x = pari(n).issquare(find_root=True)
    if root:
        if t:
            if hasattr(n, 'parent'):
                x = n.parent()(str(x))
            else:
                x = x.python()
        return t, x
    return t


def is_squarefree(n):
    """
    Returns True if and only if n is not divisible by the square of an
    integer > 1.

    EXAMPLES::

        sage: is_squarefree(100)
        False
        sage: is_squarefree(101)
        True
    """
    if n==0:
        return False
    for p, r in factor(n):
        if r>1:
            return False
    return True


#################################################################
# Euler phi function
#################################################################
class Euler_Phi:
    r"""
    Return the value of the Euler phi function on the integer n. We
    defined this to be the number of positive integers <= n that are
    relatively prime to n. Thus if n<=0 then
    ``euler_phi(n)`` is defined and equals 0.

    INPUT:


    -  ``n`` - an integer


    EXAMPLES::

        sage: euler_phi(1)
        1
        sage: euler_phi(2)
        1
        sage: euler_phi(3)
        2
        sage: euler_phi(12)
        4
        sage: euler_phi(37)
        36

    Notice that euler_phi is defined to be 0 on negative numbers and
    0.

    ::

        sage: euler_phi(-1)
        0
        sage: euler_phi(0)
        0
        sage: type(euler_phi(0))
        <type 'sage.rings.integer.Integer'>

    We verify directly that the phi function is correct for 21.

    ::

        sage: euler_phi(21)
        12
        sage: [i for i in range(21) if gcd(21,i) == 1]
        [1, 2, 4, 5, 8, 10, 11, 13, 16, 17, 19, 20]

    The length of the list of integers 'i' in range(n) such that the
    gcd(i,n) == 1 equals euler_phi(n).

    ::

        sage: len([i for i in range(21) if gcd(21,i) == 1]) == euler_phi(21)
        True

    The phi function also has a special plotting method.

    ::

        sage: P = plot(euler_phi, -3, 71)

    AUTHORS:

    - William Stein

    - Alex Clemesha (2006-01-10): some examples
    """
    def __repr__(self):
        """
        Returns a string describing this class.

        EXAMPLES::

            sage: Euler_Phi().__repr__()
            'Number of positive integers <=n but relatively prime to n'
        """
        return "Number of positive integers <=n but relatively prime to n"

    def __call__(self, n):
        """
        Calls the euler_phi function.

        EXAMPLES::

            sage: Euler_Phi()(10)
            4
            sage: Euler_Phi()(720)
            192
        """
        if n<=0:
            return ZZ(0)
        if n<=2:
            return ZZ(1)
        return ZZ(pari(n).phi())
        #return misc.mul([(p-1)*p**(r-1) for p, r in factor(n)])

    def plot(self, xmin=1, xmax=50, pointsize=30, rgbcolor=(0,0,1), join=True,
             **kwds):
        """
        Plot the Euler phi function.

        INPUT:


        -  ``xmin`` - default: 1

        -  ``xmax`` - default: 50

        -  ``pointsize`` - default: 30

        -  ``rgbcolor`` - default: (0,0,1)

        -  ``join`` - default: True; whether to join the
           points.

        -  ``**kwds`` - passed on

        EXAMPLES::

            sage: p = Euler_Phi().plot()
            sage: p.ymax()
            46.0
        """
        v = [(n,euler_phi(n)) for n in range(xmin,xmax + 1)]
        from sage.plot.all import list_plot
        P = list_plot(v, pointsize=pointsize, rgbcolor=rgbcolor, **kwds)
        if join:
            P += list_plot(v, plotjoined=True, rgbcolor=(0.7,0.7,0.7), **kwds)
        return P

euler_phi = Euler_Phi()

def crt(a,b,m=None,n=None):
    r"""
    Returns a solution to a Chinese Remainder Theorem problem.

    INPUT:

    - ``a``, ``b`` - two residues (elements of some ring for which
      extended gcd is available), or two lists, one of residues and
      one of moduli.

    - ``m``, ``n`` - (default: ``None``) two moduli, or ``None``.

    OUTPUT:

    If ``m``, ``n`` are not ``None``, returns a solution `x` to the
    simultaneous congruences `x\equiv a \bmod m` and `x\equiv b \bmod
    n`, if one exists. By the Chinese Remainder Theorem, a solution to the
    simultaneous congruences exists if and only if
    `a\equiv b\pmod{\gcd(m,n)}`. The solution `x` is only well-defined modulo
    `\text{lcm}(m,n)`.

    If ``a`` and ``b`` are lists, returns a simultaneous solution to
    the congruences `x\equiv a_i\pmod{b_i}`, if one exists.

    .. SEEALSO::

        - :func:`CRT_list`

    EXAMPLES:

    Using ``crt`` by giving it pairs of residues and moduli::

        sage: crt(2, 1, 3, 5)
        11
        sage: crt(13, 20, 100, 301)
        28013
        sage: crt([2, 1], [3, 5])
        11
        sage: crt([13, 20], [100, 301])
        28013

    You can also use upper case::

        sage: c = CRT(2,3, 3, 5); c
        8
        sage: c % 3 == 2
        True
        sage: c % 5 == 3
        True

    Note that this also works for polynomial rings::

        sage: K.<a> = NumberField(x^3 - 7)
        sage: R.<y> = K[]
        sage: f = y^2 + 3
        sage: g = y^3 - 5
        sage: CRT(1,3,f,g)
        -3/26*y^4 + 5/26*y^3 + 15/26*y + 53/26
        sage: CRT(1,a,f,g)
        (-3/52*a + 3/52)*y^4 + (5/52*a - 5/52)*y^3 + (15/52*a - 15/52)*y + 27/52*a + 25/52

    You can also do this for any number of moduli::

        sage: K.<a> = NumberField(x^3 - 7)
        sage: R.<x> = K[]
        sage: CRT([], [])
        0
        sage: CRT([a], [x])
        a
        sage: f = x^2 + 3
        sage: g = x^3 - 5
        sage: h = x^5 + x^2 - 9
        sage: k = CRT([1, a, 3], [f, g, h]); k
        (127/26988*a - 5807/386828)*x^9 + (45/8996*a - 33677/1160484)*x^8 + (2/173*a - 6/173)*x^7 + (133/6747*a - 5373/96707)*x^6 + (-6/2249*a + 18584/290121)*x^5 + (-277/8996*a + 38847/386828)*x^4 + (-135/4498*a + 42673/193414)*x^3 + (-1005/8996*a + 470245/1160484)*x^2 + (-1215/8996*a + 141165/386828)*x + 621/8996*a + 836445/386828
        sage: k.mod(f)
        1
        sage: k.mod(g)
        a
        sage: k.mod(h)
        3

    If the moduli are not coprime, a solution may not exist::

        sage: crt(4,8,8,12)
        20
        sage: crt(4,6,8,12)
        Traceback (most recent call last):
        ...
        ValueError: No solution to crt problem since gcd(8,12) does not divide 4-6

        sage: x = polygen(QQ)
        sage: crt(2,3,x-1,x+1)
        -1/2*x + 5/2
        sage: crt(2,x,x^2-1,x^2+1)
        -1/2*x^3 + x^2 + 1/2*x + 1
        sage: crt(2,x,x^2-1,x^3-1)
        Traceback (most recent call last):
        ...
        ValueError: No solution to crt problem since gcd(x^2 - 1,x^3 - 1) does not divide 2-x

        sage: crt(int(2), int(3), int(7), int(11))
        58
    """
    if isinstance(a, list):
        return CRT_list(a, b)
    if isinstance(a, (int, long)):
        a = integer.Integer(a) # otherwise we get an error at (b-a).quo_rem(g)
    g, alpha, beta = XGCD(m, n)
    q, r = (b - a).quo_rem(g)
    if r != 0:
        raise ValueError("No solution to crt problem since gcd(%s,%s) does not divide %s-%s" % (m, n, a, b))
    return (a + q*alpha*m) % lcm(m, n)

CRT = crt

def CRT_list(v, moduli):
    r""" Given a list ``v`` of elements and a list of corresponding
    ``moduli``, find a single element that reduces to each element of
    ``v`` modulo the corresponding moduli.

    .. SEEALSO::

        - :func:`crt`

    EXAMPLES::

        sage: CRT_list([2,3,2], [3,5,7])
        23
        sage: x = polygen(QQ)
        sage: c = CRT_list([3], [x]); c
        3
        sage: c.parent()
        Univariate Polynomial Ring in x over Rational Field

    It also works if the moduli are not coprime::

        sage: CRT_list([32,2,2],[60,90,150])
        452

    But with non coprime moduli there is not always a solution::

        sage: CRT_list([32,2,1],[60,90,150])
        Traceback (most recent call last):
        ...
        ValueError: No solution to crt problem since gcd(180,150) does not divide 92-1

    The arguments must be lists::

        sage: CRT_list([1,2,3],"not a list")
        Traceback (most recent call last):
        ...
        ValueError: Arguments to CRT_list should be lists
        sage: CRT_list("not a list",[2,3])
        Traceback (most recent call last):
        ...
        ValueError: Arguments to CRT_list should be lists

    The list of moduli must have the same length as the list of elements::

        sage: CRT_list([1,2,3],[2,3,5])
        23
        sage: CRT_list([1,2,3],[2,3])
        Traceback (most recent call last):
        ...
        ValueError: Arguments to CRT_list should be lists of the same length
        sage: CRT_list([1,2,3],[2,3,5,7])
        Traceback (most recent call last):
        ...
        ValueError: Arguments to CRT_list should be lists of the same length

    TESTS::

        sage: CRT([32r,2r,2r],[60r,90r,150r])
        452

    """
    if not isinstance(v,list) or not isinstance(moduli,list):
        raise ValueError, "Arguments to CRT_list should be lists"
    if len(v) != len(moduli):
        raise ValueError, "Arguments to CRT_list should be lists of the same length"
    if len(v) == 0:
        return ZZ(0)
    if len(v) == 1:
        return moduli[0].parent()(v[0])
    x = v[0]
    m = moduli[0]
    for i in range(1,len(v)):
        x = CRT(x,v[i],m,moduli[i])
        m = lcm(m,moduli[i])
    return x%m

def CRT_basis(moduli):
    r"""
    Returns a CRT basis for the given moduli.

    INPUT:

    - ``moduli`` - list of pairwise coprime moduli `m` which admit an
       extended Euclidean algorithm

    OUTPUT:

    - a list of elements `a_i` of the same length as `m` such that
      `a_i` is congruent to 1 modulo `m_i` and to 0 modulo `m_j` for
      `j\not=i`.

    .. note::

       The pairwise coprimality of the input is not checked.

    EXAMPLES::

        sage: a1 = ZZ(mod(42,5))
        sage: a2 = ZZ(mod(42,13))
        sage: c1,c2 = CRT_basis([5,13])
        sage: mod(a1*c1+a2*c2,5*13)
        42

    A polynomial example::

        sage: x=polygen(QQ)
        sage: mods = [x,x^2+1,2*x-3]
        sage: b = CRT_basis(mods)
        sage: b
        [-2/3*x^3 + x^2 - 2/3*x + 1, 6/13*x^3 - x^2 + 6/13*x, 8/39*x^3 + 8/39*x]
        sage: [[bi % mj for mj in mods] for bi in b]
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    """
    n = len(moduli)
    if n == 0:
        return []
    M = prod(moduli)
    return [((xgcd(m,M//m)[2])*(M//m))%M for m in moduli]

def CRT_vectors(X, moduli):
    r"""
    Vector form of the Chinese Remainder Theorem: given a list of integer
    vectors `v_i` and a list of coprime moduli `m_i`, find a vector `w` such
    that `w = v_i \pmod m_i` for all `i`. This is more efficient than applying
    :func:`CRT` to each entry.

    INPUT:

    -  ``X`` - list or tuple, consisting of lists/tuples/vectors/etc of
       integers of the same length
    -  ``moduli`` - list of len(X) moduli

    OUTPUT:

    -  ``list`` - application of CRT componentwise.

    EXAMPLES::

        sage: CRT_vectors([[3,5,7],[3,5,11]], [2,3])
        [3, 5, 5]

        sage: CRT_vectors([vector(ZZ, [2,3,1]), Sequence([1,7,8],ZZ)], [8,9])
        [10, 43, 17]
    """
    # First find the CRT basis:
    if len(X) == 0 or len(X[0]) == 0:
        return []
    n = len(X)
    if n != len(moduli):
        raise ValueError, "number of moduli must equal length of X"
    a = CRT_basis(moduli)
    modulus = misc.prod(moduli)
    return [sum([a[i]*X[i][j] for i in range(n)]) % modulus for j in range(len(X[0]))]

def binomial(x, m, **kwds):
    r"""
    Return the binomial coefficient

    .. math::

                \binom{x}{m} = x (x-1) \cdots (x-m+1) / m!


    which is defined for `m \in \ZZ` and any
    `x`. We extend this definition to include cases when
    `x-m` is an integer but `m` is not by

    .. math::

        \binom{x}{m}= \binom{x}{x-m}

    If `m < 0`, return `0`.

    INPUT:

    -  ``x``, ``m`` - numbers or symbolic expressions. Either ``m``
       or ``x-m`` must be an integer.

    OUTPUT: number or symbolic expression (if input is symbolic)

    EXAMPLES::

        sage: binomial(5,2)
        10
        sage: binomial(2,0)
        1
        sage: binomial(1/2, 0)
        1
        sage: binomial(3,-1)
        0
        sage: binomial(20,10)
        184756
        sage: binomial(-2, 5)
        -6
        sage: binomial(-5, -2)
        0
        sage: binomial(RealField()('2.5'), 2)
        1.87500000000000
        sage: n=var('n'); binomial(n,2)
        1/2*(n - 1)*n
        sage: n=var('n'); binomial(n,n)
        1
        sage: n=var('n'); binomial(n,n-1)
        n
        sage: binomial(2^100, 2^100)
        1

        sage: k, i = var('k,i')
        sage: binomial(k,i)
        binomial(k, i)

    If `x \in \ZZ`, there is an optional 'algorithm' parameter, which
    can be 'mpir' (faster for small values) or 'pari' (faster for
    large values)::

        sage: a = binomial(100, 45, algorithm='mpir')
        sage: b = binomial(100, 45, algorithm='pari')
        sage: a == b
        True

    TESTS:

    We test that certain binomials are very fast (this should be
    instant) -- see :trac:`3309`::

        sage: a = binomial(RR(1140000.78), 23310000)

    We test conversion of arguments to Integers -- see :trac:`6870`::

        sage: binomial(1/2,1/1)
        1/2
        sage: binomial(10^20+1/1,10^20)
        100000000000000000001
        sage: binomial(SR(10**7),10**7)
        1
        sage: binomial(3/2,SR(1/1))
        3/2

    Some floating point cases -- see :trac:`7562`, :trac:`9633`, and
    :trac:`12448`::

        sage: binomial(1.,3)
        0.000000000000000
        sage: binomial(-2.,3)
        -4.00000000000000
        sage: binomial(0.5r, 5)
        0.02734375
        sage: a = binomial(float(1001), float(1)); a
        1001.0
        sage: type(a)
        <type 'float'>
        sage: binomial(float(1000), 1001)
        0.0

    Test symbolic and uni/multivariate polynomials::

        sage: K.<x> = ZZ[]
        sage: binomial(x,3)
        1/6*x^3 - 1/2*x^2 + 1/3*x
        sage: binomial(x,3).parent()
        Univariate Polynomial Ring in x over Rational Field
        sage: K.<x,y> = Integers(7)[]
        sage: binomial(y,3)
        -y^3 + 3*y^2 - 2*y
        sage: binomial(y,3).parent()
        Multivariate Polynomial Ring in x, y over Ring of integers modulo 7
        sage: n = var('n')
        sage: binomial(n,2)
        1/2*(n - 1)*n
    """
    if isinstance(m,sage.symbolic.expression.Expression):
        try:
            # For performance reasons, we avoid to try to coerce
            # to Integer in the symbolic case (see #6870)
            m=m.pyobject()
            m=ZZ(m)
        except TypeError:
            pass
    else:
        try:
            m=ZZ(m)
        except TypeError:
            pass
    if not isinstance(m,integer.Integer):
        try:
            m = ZZ(x-m)
        except TypeError:
            try:
                return x.binomial(m)
            except AttributeError:
                pass
            raise TypeError('Either m or x-m must be an integer')
    if isinstance(x, (float, sage.rings.real_mpfr.RealNumber,
                      sage.rings.real_mpfr.RealLiteral)):
        P = parent(x)
        # Try to convert x to an integer if it *is* an integer but of type
        # float. Otherwise, pari fails below in the case m > x.
        try:
            x = ZZ(x)
        except TypeError:
            pass
        if m < 0:
            return P(0)
        return P(pari(x).binomial(m))
    if isinstance(x,sage.symbolic.expression.Expression):
        try:
            x=x.pyobject()
            x=ZZ(x)
        except TypeError:
            pass
    else:
        try:
            x=ZZ(x)
        except TypeError:
            pass
    if isinstance(x,integer.Integer):
        if m < 0 or (x >= 0 and m > x):
            return ZZ.zero()

        s = sys.maxint
        if m > s:
            m = x - m
            if m > s:
                raise ValueError("binomial not implemented for "
                        "m > " + str(s) + ".\nThis is probably OK, "
                        "since the answer would have "
                        "billions of digits.")

        return x.binomial(m, **kwds)
    try:
        P = x.parent()
    except AttributeError:
        P = type(x)
    if m < 0:
        return P(0)
    return misc.prod([x-i for i in xrange(m)])/factorial(m)

def multinomial(*ks):
    r"""
    Return the multinomial coefficient

    INPUT:

    - An arbitrary number of integer arguments `k_1,\dots,k_n`
    - A list of integers `[k_1,\dots,k_n]`

    OUTPUT:

    Returns the integer:

    .. math::

           \binom{k_1 + \cdots + k_n}{k_1, \cdots, k_n}
           =\frac{\left(\sum_{i=1}^n k_i\right)!}{\prod_{i=1}^n k_i!}
           = \prod_{i=1}^n \binom{\sum_{j=1}^i k_j}{k_i}

    EXAMPLES::

        sage: multinomial(0, 0, 2, 1, 0, 0)
        3
        sage: multinomial([0, 0, 2, 1, 0, 0])
        3
        sage: multinomial(3, 2)
        10
        sage: multinomial(2^30, 2, 1)
        618970023101454657175683075
        sage: multinomial([2^30, 2, 1])
        618970023101454657175683075

    AUTHORS:

    - Gabriel Ebner
    """
    if isinstance(ks[0],list):
        if len(ks) >1:
            raise ValueError, "multinomial takes only one list argument"
        ks=ks[0]

    s, c = 0, 1
    for k in ks:
        s += k
        c *= binomial(s, k)
    return c

def binomial_coefficients(n):
    r"""
    Return a dictionary containing pairs
    `\{(k_1,k_2) : C_{k,n}\}` where `C_{k_n}` are
    binomial coefficients and `n = k_1 + k_2`.

    INPUT:


    -  ``n`` - an integer


    OUTPUT: dict

    EXAMPLES::

        sage: sorted(binomial_coefficients(3).items())
        [((0, 3), 1), ((1, 2), 3), ((2, 1), 3), ((3, 0), 1)]

    Notice the coefficients above are the same as below::

        sage: R.<x,y> = QQ[]
        sage: (x+y)^3
        x^3 + 3*x^2*y + 3*x*y^2 + y^3

    AUTHORS:

    - Fredrik Johansson
    """
    d = {(0, n):1, (n, 0):1}
    a = 1
    for k in xrange(1, n//2+1):
        a = (a * (n-k+1))//k
        d[k, n-k] = d[n-k, k] = a
    return d

def multinomial_coefficients(m, n):
    r"""
    Return a dictionary containing pairs
    `\{(k_1, k_2, ..., k_m) : C_{k, n}\}` where
    `C_{k, n}` are multinomial coefficients such that
    `n = k_1 + k_2 + ...+ k_m`.

    INPUT:

    -  ``m`` - integer
    -  ``n`` - integer

    OUTPUT: dict

    EXAMPLES::

        sage: sorted(multinomial_coefficients(2, 5).items())
        [((0, 5), 1), ((1, 4), 5), ((2, 3), 10), ((3, 2), 10), ((4, 1), 5), ((5, 0), 1)]

    Notice that these are the coefficients of `(x+y)^5`::

        sage: R.<x,y> = QQ[]
        sage: (x+y)^5
        x^5 + 5*x^4*y + 10*x^3*y^2 + 10*x^2*y^3 + 5*x*y^4 + y^5

    ::

        sage: sorted(multinomial_coefficients(3, 2).items())
        [((0, 0, 2), 1), ((0, 1, 1), 2), ((0, 2, 0), 1), ((1, 0, 1), 2), ((1, 1, 0), 2), ((2, 0, 0), 1)]

    ALGORITHM: The algorithm we implement for computing the multinomial
    coefficients is based on the following result:

    ..math::

        \binom{n}{k_1, \cdots, k_m} =
        \frac{k_1+1}{n-k_1}\sum_{i=2}^m \binom{n}{k_1+1, \cdots, k_i-1, \cdots}

    e.g.::

        sage: k = (2, 4, 1, 0, 2, 6, 0, 0, 3, 5, 7, 1) # random value
        sage: n = sum(k)
        sage: s = 0
        sage: for i in range(1, len(k)):
        ...       ki = list(k)
        ...       ki[0] += 1
        ...       ki[i] -= 1
        ...       s += multinomial(n, *ki)
        sage: multinomial(n, *k) == (k[0] + 1) / (n - k[0]) * s
        True

    TESTS::

        sage: multinomial_coefficients(0, 0)
        {(): 1}
        sage: multinomial_coefficients(0, 3)
        {}

    """
    if not m:
        if n:
            return {}
        else:
            return {(): 1}
    if m == 2:
        return binomial_coefficients(n)
    t = [n] + [0] * (m - 1)
    r = {tuple(t): 1}
    if n:
        j = 0 # j will be the leftmost nonzero position
    else:
        j = m
    # enumerate tuples in co-lex order
    while j < m - 1:
        # compute next tuple
        tj = t[j]
        if j:
            t[j] = 0
            t[0] = tj
        if tj > 1:
            t[j + 1] += 1
            j = 0
            start = 1
            v = 0
        else:
            j += 1
            start = j + 1
            v = r[tuple(t)]
            t[j] += 1
        # compute the value
        # NB: the initialization of v was done above
        for k in xrange(start, m):
            if t[k]:
                t[k] -= 1
                v += r[tuple(t)]
                t[k] += 1
        t[0] -= 1
        r[tuple(t)] = (v * tj) // (n - t[0])
    return r

# (since trac 14496) gaussian_binomial = sage.combinat.q_analogues.q_binomial

def kronecker_symbol(x,y):
    """
    The Kronecker symbol `(x|y)`.

    INPUT:

    -  ``x`` - integer

    -  ``y`` - integer

    EXAMPLES::

        sage: kronecker_symbol(13,21)
        -1
        sage: kronecker_symbol(101,4)
        1

    IMPLEMENTATION: Using GMP.
    """
    x = QQ(x).numerator() * QQ(x).denominator()
    return ZZ(x.kronecker(y))

def kronecker(x,y):
    r"""
    Synonym for :func:`kronecker_symbol`.

    The Kronecker symbol `(x|y)`.

    INPUT:

    -  ``x`` - integer

    -  ``y`` - integer

    OUTPUT:

    - an integer

    EXAMPLES::

        sage: kronecker(3,5)
        -1
        sage: kronecker(3,15)
        0
        sage: kronecker(2,15)
        1
        sage: kronecker(-2,15)
        -1
        sage: kronecker(2/3,5)
        1
    """
    return kronecker_symbol(x,y)

def legendre_symbol(x,p):
    r"""
    The Legendre symbol `(x|p)`, for `p` prime.

    .. note::

       The :func:`kronecker_symbol` command extends the Legendre
       symbol to composite moduli and `p=2`.

    INPUT:


    -  ``x`` - integer

    -  ``p`` - an odd prime number


    EXAMPLES::

        sage: legendre_symbol(2,3)
        -1
        sage: legendre_symbol(1,3)
        1
        sage: legendre_symbol(1,2)
        Traceback (most recent call last):
        ...
        ValueError: p must be odd
        sage: legendre_symbol(2,15)
        Traceback (most recent call last):
        ...
        ValueError: p must be a prime
        sage: kronecker_symbol(2,15)
        1
        sage: legendre_symbol(2/3,7)
        -1
    """
    x = QQ(x).numerator() * QQ(x).denominator()
    p = ZZ(p)
    if not p.is_prime():
        raise ValueError, "p must be a prime"
    if p == 2:
        raise ValueError, "p must be odd"
    return x.kronecker(p)

def jacobi_symbol(a,b):
    r"""
    The Jacobi symbol of integers a and b, where b is odd.

    .. note::

       The :func:`kronecker_symbol` command extends the Jacobi
       symbol to all integers b.

    If

    `b = p_1^{e_1} * ... * p_r^{e_r}`

    then

    `(a|b) = (a|p_1)^{e_1} ... (a|p_r)^{e_r}`

    where `(a|p_j)` are Legendre Symbols.



    INPUT:

    -  ``a`` - an integer

    -  ``b`` - an odd integer

    EXAMPLES::

        sage: jacobi_symbol(10,777)
        -1
        sage: jacobi_symbol(10,5)
        0
        sage: jacobi_symbol(10,2)
        Traceback (most recent call last):
        ...
        ValueError: second input must be odd, 2 is not odd
    """

    if b%2==0:
        raise ValueError, "second input must be odd, %s is not odd"%b

    return kronecker_symbol(a,b)


def primitive_root(n, check=True):
    """
    Return a positive integer that generates the multiplicative group
    of integers modulo `n`, if one exists; otherwise, raise a
    ``ValueError``.

    A primitive root exists if `n=4` or `n=p^k` or `n=2p^k`, where `p`
    is an odd prime and `k` is a nonnegative number.

    INPUT:

    - ``n`` -- a non-zero integer
    - ``check`` -- bool (default: True); if False, then `n` is assumed
      to be a positive integer possessing a primitive root, and behavior
      is undefined otherwise.

    OUTPUT:

    A primitive root of `n`. If `n` is prime, this is the smallest
    primitive root.

    EXAMPLES::

        sage: primitive_root(23)
        5
        sage: primitive_root(-46)
        5
        sage: primitive_root(25)
        2
        sage: print [primitive_root(p) for p in primes(100)]
        [1, 2, 2, 3, 2, 2, 3, 2, 5, 2, 3, 2, 6, 3, 5, 2, 2, 2, 2, 7, 5, 3, 2, 3, 5]
        sage: primitive_root(8)
        Traceback (most recent call last):
        ...
        ValueError: no primitive root

    .. NOTE::

        It takes extra work to check if `n` has a primitive root; to
        avoid this, use ``check=False``, which may slightly speed things
        up (but could also result in undefined behavior).  For example,
        the second call below is an order of magnitude faster than the
        first:

    ::

        sage: n = 10^50 + 151   # a prime
        sage: primitive_root(n)
        11
        sage: primitive_root(n, check=False)
        11

    TESTS:

    Various special cases::

        sage: primitive_root(-1)
        0
        sage: primitive_root(0)
        Traceback (most recent call last):
        ...
        ValueError: no primitive root
        sage: primitive_root(1)
        0
        sage: primitive_root(2)
        1
        sage: primitive_root(4)
        3

    We test that various numbers without primitive roots give
    an error - see Trac 10836::

        sage: primitive_root(15)
        Traceback (most recent call last):
        ...
        ValueError: no primitive root
        sage: primitive_root(16)
        Traceback (most recent call last):
        ...
        ValueError: no primitive root
        sage: primitive_root(1729)
        Traceback (most recent call last):
        ...
        ValueError: no primitive root
        sage: primitive_root(4*7^8)
        Traceback (most recent call last):
        ...
        ValueError: no primitive root
    """
    if not check:
        return ZZ(pari(n).znprimroot())
    n = ZZ(n).abs()
    if n == 4:
        return ZZ(3)
    if n%2: # n odd
        if n.is_prime_power():
            return ZZ(pari(n).znprimroot())
    else:   # n even
        m = n // 2
        if m%2 and m.is_prime_power():
            return ZZ(pari(n).znprimroot())
    raise ValueError, "no primitive root"

def nth_prime(n):
    """

    Return the n-th prime number (1-indexed, so that 2 is the 1st prime.)

    INPUT:

    - ``n`` -- a positive integer

    OUTPUT:

    -  the n-th prime number

    EXAMPLES::

        sage: nth_prime(3)
        5
        sage: nth_prime(10)
        29

    ::

        sage: nth_prime(0)
        Traceback (most recent call last):
        ...
        ValueError: nth prime meaningless for non-positive n (=0)

    TESTS::

        sage: all(prime_pi(nth_prime(j)) == j for j in range(1, 1000, 10))
        True

    """
    return ZZ(pari.nth_prime(n))

def quadratic_residues(n):
    r"""
    Return a sorted list of all squares modulo the integer `n`
    in the range `0\leq x < |n|`.

    EXAMPLES::

        sage: quadratic_residues(11)
        [0, 1, 3, 4, 5, 9]
        sage: quadratic_residues(1)
        [0]
        sage: quadratic_residues(2)
        [0, 1]
        sage: quadratic_residues(8)
        [0, 1, 4]
        sage: quadratic_residues(-10)
        [0, 1, 4, 5, 6, 9]
        sage: v = quadratic_residues(1000); len(v);
        159
    """
    n = abs(int(n))
    X = list(set([ZZ((a*a)%n) for a in range(n/2+1)]))
    X.sort()
    return X

## This much slower than above, for obvious reasons.
## def quadratic_residues2(p):
##     return [x for x in range(p-1) if kronecker_symbol(x,p)==1]

class Moebius:
    r"""
    Returns the value of the Moebius function of abs(n), where n is an
    integer.

    DEFINITION: `\mu(n)` is 0 if `n` is not square
    free, and otherwise equals `(-1)^r`, where `n` has
    `r` distinct prime factors.

    For simplicity, if `n=0` we define `\mu(n) = 0`.

    IMPLEMENTATION: Factors or - for integers - uses the PARI C
    library.

    INPUT:


    -  ``n`` - anything that can be factored.


    OUTPUT: 0, 1, or -1

    EXAMPLES::

        sage: moebius(-5)
        -1
        sage: moebius(9)
        0
        sage: moebius(12)
        0
        sage: moebius(-35)
        1
        sage: moebius(-1)
        1
        sage: moebius(7)
        -1

    ::

        sage: moebius(0)   # potentially nonstandard!
        0

    The moebius function even makes sense for non-integer inputs.

    ::

        sage: x = GF(7)['x'].0
        sage: moebius(x+2)
        -1
    """
    def __call__(self, n):
        """
        EXAMPLES::

            sage: Moebius().__call__(7)
            -1
        """
        if isinstance(n, (int, long)):
            n = ZZ(n)
        elif not isinstance(n, integer.Integer):
            # Use a generic algorithm.
            if n < 0:
                n = -n
            F = factor(n)
            for _, e in F:
                if e >= 2:
                    return 0
            return (-1)**len(F)

        # Use fast PARI algorithm
        if n == 0:
            return integer.Integer(0)
        return ZZ(pari(n).moebius())


    def __repr__(self):
        """
        Returns a description of this function.

        EXAMPLES::

            sage: q = Moebius()
            sage: q.__repr__()
            'The Moebius function'
        """
        return "The Moebius function"

    def plot(self, xmin=0, xmax=50, pointsize=30, rgbcolor=(0,0,1), join=True,
             **kwds):
        """
        Plot the Moebius function.

        INPUT:


        -  ``xmin`` - default: 0

        -  ``xmax`` - default: 50

        -  ``pointsize`` - default: 30

        -  ``rgbcolor`` - default: (0,0,1)

        -  ``join`` - default: True; whether to join the points
           (very helpful in seeing their order).

        -  ``**kwds`` - passed on

        EXAMPLES::

            sage: p = Moebius().plot()
            sage: p.ymax()
            1.0
        """
        values = self.range(xmin, xmax + 1)
        v = [(n,values[n-xmin]) for n in range(xmin,xmax + 1)]
        from sage.plot.all import list_plot
        P = list_plot(v, pointsize=pointsize, rgbcolor=rgbcolor, **kwds)
        if join:
            P += list_plot(v, plotjoined=True, rgbcolor=(0.7,0.7,0.7), **kwds)
        return P

    def range(self, start, stop=None, step=None):
        """
        Return the Moebius function evaluated at the given range of values,
        i.e., the image of the list range(start, stop, step) under the
        Mobius function.

        This is much faster than directly computing all these values with a
        list comprehension.

        EXAMPLES::

            sage: v = moebius.range(-10,10); v
            [1, 0, 0, -1, 1, -1, 0, -1, -1, 1, 0, 1, -1, -1, 0, -1, 1, -1, 0, 0]
            sage: v == [moebius(n) for n in range(-10,10)]
            True
            sage: v = moebius.range(-1000, 2000, 4)
            sage: v == [moebius(n) for n in range(-1000,2000, 4)]
            True
        """
        if stop is None:
            start, stop = 1, int(start)
        else:
            start = int(start)
            stop = int(stop)
        if step is None:
            step = 1
        else:
            step = int(step)

        Z = integer.Integer

        if start <= 0 and 0 < stop and start % step == 0:
            return self.range(start, 0, step) + [Z(0)] +\
                   self.range(step, stop, step)

        if step == 1:
            v = pari('vector(%s, i, moebius(i-1+%s))'%(
                stop-start, start))
        else:
            n = len(range(start, stop, step)) # stupid
            v = pari('vector(%s, i, moebius(%s*(i-1) + %s))'%(
                n, step, start))
        return [Z(x) for x in v]

moebius = Moebius()

def farey(v, lim):
    """
    Return the Farey sequence associated to the floating point number
    v.

    INPUT:


    -  ``v`` - float (automatically converted to a float)

    -  ``lim`` - maximum denominator.


    OUTPUT: Results are (numerator, denominator); (1, 0) is "infinity".

    EXAMPLES::

        sage: farey(2.0, 100)
        (2, 1)
        sage: farey(2.0, 1000)
        (2, 1)
        sage: farey(2.1, 1000)
        (21, 10)
        sage: farey(2.1, 100000)
        (21, 10)
        sage: farey(pi, 100000)
        (312689, 99532)

    AUTHORS:

    - Scott David Daniels: Python Cookbook, 2nd Ed., Recipe 18.13
    """
    v = float(v)
    if v < 0:
        n, d = farey(-v, lim)
        return -n, d
    z = lim - lim    # Get a "0 of the right type" for denominator
    lower, upper = (z, z+1), (z+1, z)
    while True:
        mediant = (lower[0] + upper[0]), (lower[1] + upper[1])
        if v * mediant[1] > mediant[0]:
            if lim < mediant[1]:
                return upper
            lower = mediant
        elif v * mediant[1] == mediant[0]:
            if lim >= mediant[1]:
                return mediant
            if lower[1] < upper[1]:
                return lower
            return upper
        else:
            if lim < mediant[1]:
                return lower
            upper = mediant


## def convergents_pnqn(x):
##     """
##     Return the pairs (pn,qn) that are the numerators and denominators
##     of the partial convergents of the continued fraction of x.  We
##     include (0,1) and (1,0) at the beginning of the list (these are
##     the -2 and -1 th convergents).
##     """
##     v = pari(x).contfrac()
##     w = [(0,1), (1,0)]
##     for n in range(len(v)):
##         pn = w[n+1][0]*v[n] + w[n][0]
##         qn = w[n+1][1]*v[n] + w[n][1]
##         w.append(int(pn), int(qn))
##     return w

def continued_fraction_list(x, partial_convergents=False, bits=None, nterms=None):
    r"""
    Returns the continued fraction of x as a list.

    The continued fraction expansion of `x` are the coefficients `a_i` in

    .. math::

        x = a_1 + 1/(a_2+1/(...) ... )

    with `a_1` integer and `a_2`, `...` positive integers.

    .. note::

       This may be slow for real number input, since it's implemented in pure
       Python. For rational number input the PARI C library is used.

    .. SEEALSO::

         :func:`Hirzebruch_Jung_continued_fraction_list` for
         Hirzebruch-Jung continued fractions.

    INPUT:

    - ``x`` -- exact rational or floating-point number. The number to
      compute the continued fraction of.

    - ``partial_convergents`` -- boolean. Whether to return the partial convergents.

    - ``bits`` -- integer. the precision of the real interval field
      that is used internally.

    - ``nterms`` -- integer. The upper bound on the number of terms in
      the continued fraction expansion to return.

    OUTPUT:

    A lits of integers, the coefficients in the continued fraction
    expansion of ``x``. If ``partial_convergents=True`` is passed, a
    pair containing the coefficient list and the partial convergents
    list is returned.

    EXAMPLES::

        sage: continued_fraction_list(45/17)
        [2, 1, 1, 1, 5]
        sage: continued_fraction_list(e, bits=20)
        [2, 1, 2, 1, 1, 4, 1, 1]
        sage: continued_fraction_list(e, bits=30)
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8]
        sage: continued_fraction_list(sqrt(2))
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: continued_fraction_list(sqrt(4/19))
        [0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1]
        sage: continued_fraction_list(RR(pi), partial_convergents=True)
        ([3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3],
         [(3, 1),
          (22, 7),
          (333, 106),
          (355, 113),
          (103993, 33102),
          (104348, 33215),
          (208341, 66317),
          (312689, 99532),
          (833719, 265381),
          (1146408, 364913),
          (4272943, 1360120),
          (5419351, 1725033),
          (80143857, 25510582),
          (245850922, 78256779)])
        sage: continued_fraction_list(e)
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1]
        sage: continued_fraction_list(RR(e))
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1]
        sage: continued_fraction_list(RealField(200)(e))
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1,
         14, 1, 1, 16, 1, 1, 18, 1, 1, 20, 1, 1, 22, 1, 1, 24, 1, 1,
         26, 1, 1, 28, 1, 1, 30, 1, 1, 32, 1, 1, 34, 1, 1, 36, 1, 1, 38, 1, 1]

    TESTS::

        sage: continued_fraction_list(1 + 10^-10, nterms=3)
        [1, 10000000000]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, bits=10, nterms=3)
        [1, 100000000000000000000, 2688]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, bits=10, nterms=5)
        [1, 100000000000000000000, 2688, 8, 1]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, bits=1000, nterms=5)
        [1, 100000000000000000000, 2688, 8, 1]
    """
    if isinstance(x, (integer.Integer, int, long)):
        if partial_convergents:
            return [x], [(x,1)]
        else:
            return [x]

    if isinstance(x, sage.rings.rational.Rational):
        if bits is not None and nterms is None:
            x = RealIntervalField(bits)(x)
        else:
            # PARI is faster than the pure Python below, but doesn't give us the convergents.
            v = pari(x).contfrac().python()
            if nterms is not None:
                v = v[:nterms]
            if partial_convergents:
                w = [(0,1), (1,0)]
                for a in v:
                    pn = a*w[-1][0] + w[-2][0]
                    qn = a*w[-1][1] + w[-2][1]
                    w.append((pn, qn))
                return v, w[2:]
            else:
                return v

    # Work in interval field, increasing precision as needed.
    if bits is None:
        try:
            bits = x.prec()
        except AttributeError:
            bits = 53
    RIF = RealIntervalField(bits)
    v = []
    w = [(0,1), (1,0)]
    orig, x = x, RIF(x)

    while True:
        try:
            a = x.unique_floor()
        except ValueError:
            # Either we're done or we need more precision.
            if nterms is None:
                break
            else:
                RIF = RIF.to_prec(2*RIF.prec())
                x = RIF(orig)
                for a in v: x = ~(x-a)
                continue
        if partial_convergents:
            pn = a*w[-1][0] + w[-2][0]
            qn = a*w[-1][1] + w[-2][1]
            w.append((pn, qn))
        v.append(a)
        if x == a or nterms is not None and len(v) >= nterms:
            break
        x = ~(x-a)

    if partial_convergents:
        return v, w[2:]
    else:
        return v


def Hirzebruch_Jung_continued_fraction_list(x, bits=None, nterms=None):
    r"""
    Return the Hirzebruch-Jung continued fraction of ``x`` as a list.

    The Hirzebruch-Jung continued fraction of `x` is similar to the
    ordinary continued fraction expansion, but with minus signs. That
    is, the coefficients `a_i` in

    .. math::

        x = a_1 - 1/(a_2-1/(...) ... )

    with `a_1` integer and `a_2`, `...` positive integers.

    .. SEEALSO::

         :func:`continued_fraction_list` for ordinary continued fractions.

    INPUT:

    - ``x`` -- exact rational or something that can be numerically
      evaluated. The number to compute the continued fraction of.

    - ``bits`` -- integer (default: the precision of ``x``). the
      precision of the real interval field that is used
      internally. This is only used if ``x`` is not an exact fraction.

    - ``nterms`` -- integer (default: None). The upper bound on the
      number of terms in the continued fraction expansion to return.

    OUTPUT:

    A lits of integers, the coefficients in the Hirzebruch-Jung continued
    fraction expansion of ``x``.

    EXAMPLES::

        sage: Hirzebruch_Jung_continued_fraction_list(17/11)
        [2, 3, 2, 2, 2, 2]
        sage: Hirzebruch_Jung_continued_fraction_list(45/17)
        [3, 3, 6]
        sage: Hirzebruch_Jung_continued_fraction_list(e, bits=20)
        [3, 4, 3, 2, 2, 2, 3, 7]
        sage: Hirzebruch_Jung_continued_fraction_list(e, bits=30)
        [3, 4, 3, 2, 2, 2, 3, 8, 3, 2, 2, 2, 2, 2, 2, 2, 3]
        sage: Hirzebruch_Jung_continued_fraction_list(sqrt(2), bits=100)
        [2, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4,
         2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 2]
        sage: Hirzebruch_Jung_continued_fraction_list(sqrt(4/19))
        [1, 2, 7, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 7,
         2, 2, 2, 7, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: Hirzebruch_Jung_continued_fraction_list(pi)
        [4, 2, 2, 2, 2, 2, 2, 17, 294, 3, 4, 5, 16, 2, 2]
        sage: Hirzebruch_Jung_continued_fraction_list(e)
        [3, 4, 3, 2, 2, 2, 3, 8, 3, 2, 2, 2, 2, 2, 2, 2,
         3, 12, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 10]
        sage: Hirzebruch_Jung_continued_fraction_list(e, nterms=20)
        [3, 4, 3, 2, 2, 2, 3, 8, 3, 2, 2, 2, 2, 2, 2, 2, 3, 12, 3, 2]
        sage: len(_) == 20
        True

    TESTS::

        sage: Hirzebruch_Jung_continued_fraction_list(1 - 10^-10, nterms=3)
        [1, 10000000000]
        sage: Hirzebruch_Jung_continued_fraction_list(1 - 10^-10 - e^-100, bits=100, nterms=5)
        [1, 10000000000]
        sage: Hirzebruch_Jung_continued_fraction_list(1 - 10^-20 - e^-100, bits=1000, nterms=5)
        [1, 100000000000000000000, 2689, 2, 2]
   """
    if not isinstance(x, sage.rings.rational.Rational):
        try:
            x = QQ(x)
        except TypeError:
            # Numerically evaluate x
            if bits is None:
                try:
                    bits = x.prec()
                except AttributeError:
                    bits = 53
            x = QQ(x.n(bits))
    v = []
    while True:
        div, mod = divmod(x.numerator(), x.denominator())
        if mod == 0:
            v.append(div)
            break
        v.append(div+1)
        if nterms is not None and len(v) >= nterms:
            break
        x = 1/(div+1-x)
    return v


def convergent(v, n):
    r"""
    Return the n-th continued fraction convergent of the continued
    fraction defined by the sequence of integers v. We assume
    `n \geq 0`.

    INPUT:


    -  ``v`` - list of integers

    -  ``n`` - integer


    OUTPUT: a rational number

    If the continued fraction integers are

    .. math::

       v = [a_0, a_1, a_2, \ldots, a_k]


    then ``convergent(v,2)`` is the rational number

    .. math::

       a_0 + 1/a_1

    and ``convergent(v,k)`` is the rational number

    .. math::

       a1 + 1/(a2+1/(...) ... )

    represented by the continued fraction.

    EXAMPLES::

        sage: convergent([2, 1, 2, 1, 1, 4, 1, 1], 7)
        193/71
    """
    if hasattr(v, 'convergent'):
        return v.convergent(n)
    i = int(n)
    x = QQ(v[i])
    i -= 1
    while i >= 0:
        x = QQ(v[i]) + 1/x
        i -= 1
    return x


def convergents(v):
    """
    Return all the partial convergents of a continued fraction defined
    by the sequence of integers v.

    If v is not a list, compute the continued fraction of v and return
    its convergents (this is potentially much faster than calling
    continued_fraction first, since continued fractions are
    implemented using PARI and there is overhead moving the answer back
    from PARI).

    INPUT:


    -  ``v`` - list of integers or a rational number


    OUTPUT:


    -  ``list`` - of partial convergents, as rational
       numbers


    EXAMPLES::

        sage: convergents([2, 1, 2, 1, 1, 4, 1, 1])
        [2, 3, 8/3, 11/4, 19/7, 87/32, 106/39, 193/71]
    """
    if hasattr(v, 'convergents'):
        return v.convergents()
    if not isinstance(v, list):
        v = pari(v).contfrac()
    w = [(0,1), (1,0)]
    for n in range(len(v)):
        pn = w[n+1][0]*v[n] + w[n][0]
        qn = w[n+1][1]*v[n] + w[n][1]
        w.append((pn, qn))
    return [QQ(x) for x in w[2:]]


## def continuant(v, n=None):
##     """
##     Naive implementation with recursion.
##
##     See Graham, Knuth and Patashnik: Concrete Mathematics: 6.7 Continuants
##     """
##     m = len(v)
##     if n == None or m < n:
##         n = m
##     if n == 0:
##         return 1
##     if n == 1:
##         return v[0]
##     return continuant(v[:n-1])*v[n-1] + continuant(v[:n-2])

def continuant(v, n=None):
    r"""
    Function returns the continuant of the sequence `v` (list
    or tuple).

    Definition: see Graham, Knuth and Patashnik, *Concrete Mathematics*,
    section 6.7: Continuants. The continuant is defined by

    - `K_0() = 1`
    - `K_1(x_1) = x_1`
    - `K_n(x_1, \cdots, x_n) = K_{n-1}(x_n, \cdots x_{n-1})x_n + K_{n-2}(x_1,  \cdots, x_{n-2})`

    If ``n = None`` or ``n > len(v)`` the default
    ``n = len(v)`` is used.

    INPUT:

    -  ``v`` - list or tuple of elements of a ring
    -  ``n`` - optional integer

    OUTPUT: element of ring (integer, polynomial, etcetera).

    EXAMPLES::

        sage: continuant([1,2,3])
        10
        sage: p = continuant([2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10])
        sage: q = continuant([1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10])
        sage: p/q
        517656/190435
        sage: convergent([2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10],14)
        517656/190435
        sage: x = PolynomialRing(RationalField(),'x',5).gens()
        sage: continuant(x)
        x0*x1*x2*x3*x4 + x0*x1*x2 + x0*x1*x4 + x0*x3*x4 + x2*x3*x4 + x0 + x2 + x4
        sage: continuant(x, 3)
        x0*x1*x2 + x0 + x2
        sage: continuant(x,2)
        x0*x1 + 1

    We verify the identity

    .. math::

        K_n(z,z,\cdots,z) = \sum_{k=0}^n \binom{n-k}{k} z^{n-2k}

    for `n = 6` using polynomial arithmetic::

        sage: z = QQ['z'].0
        sage: continuant((z,z,z,z,z,z,z,z,z,z,z,z,z,z,z),6)
        z^6 + 5*z^4 + 6*z^2 + 1

        sage: continuant(9)
        Traceback (most recent call last):
        ...
        TypeError: object of type 'sage.rings.integer.Integer' has no len()

    AUTHORS:

    - Jaap Spies (2007-02-06)
    """
    m = len(v)
    if n == None or m < n:
        n = m
    if n == 0:
        return 1
    if n == 1:
        return v[0]
    a, b = 1, v[0]
    for k in range(1,n):
        a, b = b, a + b*v[k]
    return b

def number_of_divisors(n):
    """
    Return the number of divisors of the integer n.

    INPUT:

    - ``n`` - a nonzero integer

    OUTPUT:

    - an integer, the number of divisors of n

    EXAMPLES::

        sage: number_of_divisors(100)
        9
        sage: number_of_divisors(-720)
        30
    """
    m = ZZ(n)
    if m.is_zero():
        raise ValueError, "input must be nonzero"
    return ZZ(pari(m).numdiv())



def hilbert_symbol(a, b, p, algorithm="pari"):
    """
    Returns 1 if `ax^2 + by^2` `p`-adically represents
    a nonzero square, otherwise returns `-1`. If either a or b
    is 0, returns 0.

    INPUT:


    -  ``a, b`` - integers

    -  ``p`` - integer; either prime or -1 (which
       represents the archimedean place)

    -  ``algorithm`` - string

       -  ``'pari'`` - (default) use the PARI C library

       -  ``'direct'`` - use a Python implementation

       -  ``'all'`` - use both PARI and direct and check that
          the results agree, then return the common answer


    OUTPUT: integer (0, -1, or 1)

    EXAMPLES::

        sage: hilbert_symbol (-1, -1, -1, algorithm='all')
        -1
        sage: hilbert_symbol (2,3, 5, algorithm='all')
        1
        sage: hilbert_symbol (4, 3, 5, algorithm='all')
        1
        sage: hilbert_symbol (0, 3, 5, algorithm='all')
        0
        sage: hilbert_symbol (-1, -1, 2, algorithm='all')
        -1
        sage: hilbert_symbol (1, -1, 2, algorithm='all')
        1
        sage: hilbert_symbol (3, -1, 2, algorithm='all')
        -1

        sage: hilbert_symbol(QQ(-1)/QQ(4), -1, 2) == -1
        True
        sage: hilbert_symbol(QQ(-1)/QQ(4), -1, 3) == 1
        True

    AUTHORS:

    - William Stein and David Kohel (2006-01-05)
    """
    p = ZZ(p)
    if p != -1 and not p.is_prime():
        raise ValueError, "p must be prime or -1"
    a = QQ(a).numerator() * QQ(a).denominator()
    b = QQ(b).numerator() * QQ(b).denominator()

    if algorithm == "pari":
        return ZZ(pari(a).hilbert(b,p))

    elif algorithm == 'direct':
        if a == 0 or b == 0:
            return ZZ(0)

        p = ZZ(p)
        one = ZZ(1)

        if p != -1:
            p_sqr = p**2
            while a%p_sqr == 0: a //= p_sqr
            while b%p_sqr == 0: b //= p_sqr

        if p != 2 and True in ( kronecker(x,p) == 1 for x in (a,b,a+b) ):
            return one
        if a%p == 0:
            if b%p == 0:
                return hilbert_symbol(p,-(b//p),p)*hilbert_symbol(a//p,b,p)
            elif p == 2 and (b%4) == 3:
                if kronecker(a+b,p) == -1:
                    return -one
            elif kronecker(b,p) == -1:
                return -one
        elif b%p == 0:
            if p == 2 and (a%4) == 3:
                if kronecker(a+b,p) == -1:
                    return -one
            elif kronecker(a,p) == -1:
                return -one
        elif p == 2 and (a%4) == 3 and (b%4) == 3:
            return -one
        return one
    elif algorithm == 'all':
        ans_pari = hilbert_symbol(a,b,p,algorithm='pari')
        ans_direct = hilbert_symbol(a,b,p,algorithm='direct')
        if ans_pari != ans_direct:
            raise RuntimeError, "There is a bug in hilbert_symbol; two ways of computing the Hilbert symbol (%s,%s)_%s disagree"%(a,b,p)
        return ans_pari
    else:
        raise ValueError, "Algorithm %s not defined"%algorithm


def hilbert_conductor(a, b):
    """
    This is the product of all (finite) primes where the Hilbert symbol is -1.
    What is the same, this is the (reduced) discriminant of the quaternion
    algebra `(a,b)` over `\QQ`.

    INPUT:

    - ``a``, ``b`` -- integers

    OUTPUT:

    - squarefree positive integer

    EXAMPLES::

        sage: hilbert_conductor(-1, -1)
        2
        sage: hilbert_conductor(-1, -11)
        11
        sage: hilbert_conductor(-2, -5)
        5
        sage: hilbert_conductor(-3, -17)
        17

    AUTHOR:

    - Gonzalo Tornaria (2009-03-02)
    """
    a, b = ZZ(a), ZZ(b)
    d = ZZ(1)
    for p in union(union( [2], prime_divisors(a) ), prime_divisors(b) ):
        if hilbert_symbol(a, b, p) == -1:
            d *= p
    return d

def hilbert_conductor_inverse(d):
    """
    Finds a pair of integers `(a,b)` such that ``hilbert_conductor(a,b) == d``.
    The quaternion algebra `(a,b)` over `\QQ` will then have (reduced)
    discriminant `d`.

    INPUT:

    - ``d`` -- square-free positive integer

    OUTPUT: pair of integers

    EXAMPLES::

        sage: hilbert_conductor_inverse(2)
        (-1, -1)
        sage: hilbert_conductor_inverse(3)
        (-1, -3)
        sage: hilbert_conductor_inverse(6)
        (-1, 3)
        sage: hilbert_conductor_inverse(30)
        (-3, -10)
        sage: hilbert_conductor_inverse(4)
        Traceback (most recent call last):
        ...
        ValueError: d needs to be squarefree
        sage: hilbert_conductor_inverse(-1)
        Traceback (most recent call last):
        ...
        ValueError: d needs to be positive

    AUTHOR:

    - Gonzalo Tornaria (2009-03-02)

    TESTS::

        sage: for i in xrange(100):
        ...     d = ZZ.random_element(2**32).squarefree_part()
        ...     if hilbert_conductor(*hilbert_conductor_inverse(d)) != d:
        ...         print "hilbert_conductor_inverse failed for d =", d
    """
    Z = ZZ
    d = Z(d)
    if d <= 0:
        raise ValueError, "d needs to be positive"
    if d == 1:
        return (Z(-1), Z(1))
    if d == 2:
        return (Z(-1), Z(-1))
    if d.is_prime():
        if d%4 == 3:
            return (Z(-1), -d)
        if d%8 == 5:
            return (Z(-2), -d)
        q = 3
        while q%4 != 3 or kronecker_symbol(d,q) != -1:
            q = next_prime(q)
        return (Z(-q), -d)
    else:
        mo = moebius(d)
        if mo == 0:
            raise ValueError, "d needs to be squarefree"
        if d % 2 == 0 and mo*d % 16 != 2:
            dd = mo * d / 2
        else:
            dd = mo * d
        q = 1
        while hilbert_conductor(-q, dd) != d:
            q+=1;
        if dd%q == 0:
            dd /= q
        return (Z(-q), Z(dd))


##############################################################################
##  falling and rising factorials
##  By Jaap Spies
##
##       Copyright (C) 2006 Jaap Spies <j.spies@hccnet.nl>
##      Copyright (C) 2006 William Stein <wstein@gmail.com>
##
## Distributed under the terms of the GNU General Public License (GPL)
##                  http://www.gnu.org/licenses/
##############################################################################


def falling_factorial(x, a):
    r"""
    Returns the falling factorial `(x)_a`.

    The notation in the literature is a mess: often `(x)_a`,
    but there are many other notations: GKP: Concrete Mathematics uses
    `x^{\underline{a}}`.

    Definition: for integer `a \ge 0` we have
    `x(x-1) \cdots (x-a+1)`. In all other cases we use the
    GAMMA-function: `\frac {\Gamma(x+1)} {\Gamma(x-a+1)}`.

    INPUT:

    -  ``x`` - element of a ring

    -  ``a`` - a non-negative integer or

    OR

    -  ``x and a`` - any numbers

    OUTPUT: the falling factorial

    EXAMPLES::

        sage: falling_factorial(10, 3)
        720
        sage: falling_factorial(10, RR('3.0'))
        720.000000000000
        sage: falling_factorial(10, RR('3.3'))
        1310.11633396601
        sage: falling_factorial(10, 10)
        3628800
        sage: factorial(10)
        3628800
        sage: a = falling_factorial(1+I, I); a
        gamma(I + 2)
        sage: CC(a)
        0.652965496420167 + 0.343065839816545*I
        sage: falling_factorial(1+I, 4)
        4*I + 2
        sage: falling_factorial(I, 4)
        -10

    ::

        sage: M = MatrixSpace(ZZ, 4, 4)
        sage: A = M([1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
        sage: falling_factorial(A, 2) # A(A - I)
        [  1   0  10  10]
        [  1   0  10  10]
        [ 20   0 101 100]
        [  2   0  11  10]

    ::

        sage: x = ZZ['x'].0
        sage: falling_factorial(x, 4)
        x^4 - 6*x^3 + 11*x^2 - 6*x

    AUTHORS:

    - Jaap Spies (2006-03-05)
    """
    if isinstance(a, (integer.Integer, int, long)) and a >= 0:
        return misc.prod([(x - i) for i in range(a)])
    from sage.functions.all import gamma
    return gamma(x+1) / gamma(x-a+1)

def rising_factorial(x, a):
    r"""
    Returns the rising factorial `(x)^a`.

    The notation in the literature is a mess: often `(x)^a`,
    but there are many other notations: GKP: Concrete Mathematics uses
    `x^{\overline{a}}`.

    The rising factorial is also known as the Pochhammer symbol, see
    Maple and Mathematica.

    Definition: for integer `a \ge 0` we have
    `x(x+1) \cdots (x+a-1)`. In all other cases we use the
    GAMMA-function: `\frac {\Gamma(x+a)} {\Gamma(x)}`.

    INPUT:


    -  ``x`` - element of a ring

    -  ``a`` - a non-negative integer or

    -  ``x and a`` - any numbers


    OUTPUT: the rising factorial

    EXAMPLES::

        sage: rising_factorial(10,3)
        1320

    ::

        sage: rising_factorial(10,RR('3.0'))
        1320.00000000000

    ::

        sage: rising_factorial(10,RR('3.3'))
        2826.38895824964

    ::

        sage: a = rising_factorial(1+I, I); a
        gamma(2*I + 1)/gamma(I + 1)
        sage: CC(a)
        0.266816390637832 + 0.122783354006372*I

    ::

        sage: a = rising_factorial(I, 4); a
        -10

    See falling_factorial(I, 4).

    ::

        sage: x = polygen(ZZ)
        sage: rising_factorial(x, 4)
        x^4 + 6*x^3 + 11*x^2 + 6*x

    AUTHORS:

    - Jaap Spies (2006-03-05)
    """
    if isinstance(a, (integer.Integer, int, long)) and a >= 0:
        return misc.prod([(x + i) for i in range(a)])
    from sage.functions.all import gamma
    return gamma(x+a) / gamma(x)


def integer_ceil(x):
    """
    Return the ceiling of x.

    EXAMPLES::

        sage: integer_ceil(5.4)
        6
    """
    try:
        return sage.rings.all.Integer(x.ceil())
    except AttributeError:
        try:
            return sage.rings.all.Integer(int(math.ceil(float(x))))
        except TypeError:
            pass
    raise NotImplementedError, "computation of floor of %s not implemented"%repr(x)

def integer_floor(x):
    r"""
    Return the largest integer `\leq x`.

    INPUT:


    -  ``x`` - an object that has a floor method or is
       coercible to int


    OUTPUT: an Integer

    EXAMPLES::

        sage: integer_floor(5.4)
        5
        sage: integer_floor(float(5.4))
        5
        sage: integer_floor(-5/2)
        -3
        sage: integer_floor(RDF(-5/2))
        -3
    """
    try:
        return ZZ(x.floor())
    except AttributeError:
        try:
            return ZZ(int(math.floor(float(x))))
        except TypeError:
            pass
    raise NotImplementedError, "computation of floor of %s not implemented"%x



def two_squares(n, algorithm='gap'):
    """
    Write the integer n as a sum of two integer squares if possible;
    otherwise raise a ValueError.

    EXAMPLES::

        sage: two_squares(389)
        (10, 17)
        sage: two_squares(7)
        Traceback (most recent call last):
        ...
        ValueError: 7 is not a sum of two squares
        sage: a,b = two_squares(2009); a,b
        (28, 35)
        sage: a^2 + b^2
        2009

    TODO: Create an implementation using PARI's continued fraction
    implementation.
    """
    from sage.rings.all import Integer
    n = Integer(n)

    if algorithm == 'gap':
        import sage.interfaces.gap as gap
        a = gap.gap.eval('TwoSquares(%s)'%n)
        if a == 'fail':
            raise ValueError, "%s is not a sum of two squares"%n
        x, y = eval(a)
        return Integer(x), Integer(y)
    else:
        raise RuntimeError, "unknown algorithm '%s'"%algorithm

def _brute_force_four_squares(n):
    """
    Brute force search for decomposition into a sum of four squares,
    for cases that the main algorithm fails to handle.

    INPUT:

    - ``n`` - a positive integer

    OUTPUT:

    - a list of four numbers whose squares sum to n

    EXAMPLES::

        sage: from sage.rings.arith import _brute_force_four_squares
        sage: _brute_force_four_squares(567)
        [1, 1, 6, 23]
    """
    from math import sqrt
    for i in range(0,int(sqrt(n))+1):
        for j in range(i,int(sqrt(n-i**2))+1):
            for k in range(j, int(sqrt(n-i**2-j**2))+1):
                rem = n-i**2-j**2-k**2
                if rem >= 0:
                    l = int(sqrt(rem))
                    if rem-l**2==0:
                        return [i,j,k,l]

def four_squares(n):
    """
    Computes the decomposition into the sum of four squares,
    using an algorithm described by Peter Schorn at:
    http://www.schorn.ch/howto.html.

    INPUT:

    - ``n`` - an integer

    OUTPUT:

    - a list of four numbers whose squares sum to n

    EXAMPLES::

        sage: four_squares(3)
        [0, 1, 1, 1]
        sage: four_squares(130)
        [0, 0, 3, 11]
        sage: four_squares(1101011011004)
        [2, 1049178, 2370, 15196]
        sage: sum([i-sum([q^2 for q in four_squares(i)]) for i in range(2,10000)]) # long time
        0
    """
    from sage.rings.finite_rings.integer_mod import mod
    from math import sqrt
    from sage.rings.arith import _brute_force_four_squares
    try:
        ts = two_squares(n)
        return [0,0,ts[0],ts[1]]
    except ValueError:
        pass
    m = n
    v = 0
    while mod(m,4) == 0:
        v = v +1
        m = m // 4
    if mod(m,8) == 7:
        d = 1
        m = m - 1
    else:
        d = 0
    if mod(m,8)==3:
        x = int(sqrt(m))
        if mod(x,2) == 0:
            x = x - 1
        p = (m-x**2) // 2
        while not is_prime(p):
            x = x - 2
            p = (m-x**2) // 2
            if x < 0:
            # fall back to brute force
                m = m + d
                return [2**v*q for q in _brute_force_four_squares(m)]
        y,z = two_squares(p)
        return [2**v*q for q in [d,x,y+z,abs(y-z)]]
    x = int(sqrt(m))
    p = m - x**2
    if p == 1:
        return[2**v*q for q in [d,0,x,1]]
    while not is_prime(p):
        x = x - 1
        p = m - x**2
        if x < 0:
            # fall back to brute force
            m = m + d
            return [2**v*q for q in _brute_force_four_squares(m)]
    y,z = two_squares(p)
    return [2**v*q for q in [d,x,y,z]]


def subfactorial(n):
    r"""
    Subfactorial or rencontres numbers, or derangements: number of
    permutations of `n` elements with no fixed points.

    INPUT:


    -  ``n`` - non negative integer


    OUTPUT:


    -  ``integer`` - function value


    EXAMPLES::

        sage: subfactorial(0)
        1
        sage: subfactorial(1)
        0
        sage: subfactorial(8)
        14833

    AUTHORS:

    - Jaap Spies (2007-01-23)
    """
    return factorial(n)*sum(((-1)**k)/factorial(k) for k in range(n+1))


def is_power_of_two(n):
    r"""
    This function returns True if and only if `n` is a power of
    2

    INPUT:


    -  ``n`` - integer


    OUTPUT:


    -  ``True`` - if n is a power of 2

    -  ``False`` - if not


    EXAMPLES::

        sage: is_power_of_two(1024)
        True

    ::

        sage: is_power_of_two(1)
        True

    ::

        sage: is_power_of_two(24)
        False

    ::

        sage: is_power_of_two(0)
        False

    ::

        sage: is_power_of_two(-4)
        False

    AUTHORS:

    - Jaap Spies (2006-12-09)
    """
    # modification of is2pow(n) from the Programming Guide
    while n > 0 and n%2 == 0:
        n = n >> 1
    return n == 1

def differences(lis, n=1):
    """
    Returns the `n` successive differences of the elements in
    `lis`.

    EXAMPLES::

        sage: differences(prime_range(50))
        [1, 2, 2, 4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4]
        sage: differences([i^2 for i in range(1,11)])
        [3, 5, 7, 9, 11, 13, 15, 17, 19]
        sage: differences([i^3 + 3*i for i in range(1,21)])
        [10, 22, 40, 64, 94, 130, 172, 220, 274, 334, 400, 472, 550, 634, 724, 820, 922, 1030, 1144]
        sage: differences([i^3 - i^2 for i in range(1,21)], 2)
        [10, 16, 22, 28, 34, 40, 46, 52, 58, 64, 70, 76, 82, 88, 94, 100, 106, 112]
        sage: differences([p - i^2 for i, p in enumerate(prime_range(50))], 3)
        [-1, 2, -4, 4, -4, 4, 0, -6, 8, -6, 0, 4]

    AUTHORS:

    - Timothy Clemans (2008-03-09)
    """
    n = ZZ(n)
    if n < 1:
        raise ValueError, 'n must be greater than 0'
    lis = [lis[i + 1] - num for i, num in enumerate(lis[:-1])]
    if n == 1:
        return lis
    return differences(lis, n - 1)

def _cmp_complex_for_display(a, b):
    r"""
    Compare two complex numbers in a "pretty" (but mathematically
    meaningless) fashion, for display only.

    Real numbers (with a zero imaginary part) come before complex numbers,
    and are sorted.  Complex numbers are sorted by their real part
    unless their real parts are quite close, in which case they are
    sorted by their imaginary part.

    EXAMPLES::

        sage: import sage.rings.arith
        sage: cmp_c = sage.rings.arith._cmp_complex_for_display
        sage: teeny = 3e-11
        sage: cmp_c(CC(5), CC(3, 3))
        -1
        sage: cmp_c(CC(3), CC(5, 5))
        -1
        sage: cmp_c(CC(5), CC(3))
        1
        sage: cmp_c(CC(teeny, -1), CC(-teeny, 1))
        -1
        sage: cmp_c(CC(teeny, 1), CC(-teeny, -1))
        1
        sage: cmp_c(CC(0, 1), CC(1, 0.5))
        -1
        sage: cmp_c(CC(3+teeny, -1), CC(3-teeny, 1))
        -1
        sage: CIF200 = ComplexIntervalField(200)
        sage: cmp_c(CIF200(teeny, -1), CIF200(-teeny, 1))
        -1
        sage: cmp_c(CIF200(teeny, 1), CIF200(-teeny, -1))
        1
        sage: cmp_c(CIF200(0, 1), CIF200(1, 0.5))
        -1
        sage: cmp_c(CIF200(3+teeny, -1), CIF200(3-teeny, 1))
        -1
    """
    ar = a.real(); br = b.real()
    ai = a.imag(); bi = b.imag()
    epsilon = ar.parent()(1e-10)
    if ai:
        if bi:
            if abs(br) < epsilon:
                if abs(ar) < epsilon:
                    return cmp(ai, bi)
                return cmp(ar, 0)
            if abs((ar - br) / br) < epsilon:
                return cmp(ai, bi)
            return cmp(ar, br)
        else:
            return 1
    else:
        if bi:
            return -1
        else:
            return cmp(ar, br)

def sort_complex_numbers_for_display(nums):
    r"""
    Given a list of complex numbers (or a list of tuples, where the
    first element of each tuple is a complex number), we sort the list
    in a "pretty" order.  First come the real numbers (with zero
    imaginary part), then the complex numbers sorted according to
    their real part.  If two complex numbers have a real part which is
    sufficiently close, then they are sorted according to their
    imaginary part.

    This is not a useful function mathematically (not least because
    there's no principled way to determine whether the real components
    should be treated as equal or not).  It is called by various
    polynomial root-finders; its purpose is to make doctest printing
    more reproducible.

    We deliberately choose a cumbersome name for this function to
    discourage use, since it is mathematically meaningless.

    EXAMPLES::

        sage: import sage.rings.arith
        sage: sort_c = sort_complex_numbers_for_display
        sage: nums = [CDF(i) for i in range(3)]
        sage: for i in range(3):
        ...       nums.append(CDF(i + RDF.random_element(-3e-11, 3e-11),
        ...                       RDF.random_element()))
        ...       nums.append(CDF(i + RDF.random_element(-3e-11, 3e-11),
        ...                       RDF.random_element()))
        sage: shuffle(nums)
        sage: sort_c(nums)
        [0.0, 1.0, 2.0, -2.862406201e-11 - 0.708874026302*I, 2.2108362707e-11 - 0.436810529675*I, 1.00000000001 - 0.758765473764*I, 0.999999999976 - 0.723896589334*I, 1.99999999999 - 0.456080101207*I, 1.99999999999 + 0.609083628313*I]
    """
    if len(nums) == 0:
        return nums

    if isinstance(nums[0], tuple):
        return sorted(nums, cmp=_cmp_complex_for_display, key=lambda t: t[0])
    else:
        return sorted(nums, cmp=_cmp_complex_for_display)

def fundamental_discriminant(D):
    r"""
    Return the discriminant of the quadratic extension
    `K=Q(\sqrt{D})`, i.e. an integer d congruent to either 0 or
    1, mod 4, and such that, at most, the only square dividing it is
    4.

    INPUT:

    - ``D`` - an integer

    OUTPUT:

    - an integer, the fundamental discriminant

    EXAMPLES::

        sage: fundamental_discriminant(102)
        408
        sage: fundamental_discriminant(720)
        5
        sage: fundamental_discriminant(2)
        8

    """
    from sage.rings.all import Integer
    D = Integer(D)
    D = D.squarefree_part()
    if D%4 == 1:
        return D
    return 4*D

def squarefree_divisors(x):
    """
    Iterator over the squarefree divisors (up to units) of the element x.

    Depends on the output of the prime_divisors function.

    INPUT::

        x -- an element of any ring for which the prime_divisors
             function works.

    EXAMPLES::

        sage: list(squarefree_divisors(7))
        [1, 7]
        sage: list(squarefree_divisors(6))
        [1, 2, 3, 6]
        sage: list(squarefree_divisors(12))
        [1, 2, 3, 6]

    """
    p_list = prime_divisors(x)
    from sage.misc.misc import powerset
    for a in powerset(p_list):
        yield prod(a,1)

def dedekind_sum(p, q, algorithm='default'):
    r"""
    Return the Dedekind sum `s(p,q)` defined for integers `p`, `q` as

    .. MATH::

        s(p,q) = \sum_{i=0}^{q-1} \left(\!\left(\frac{i}{q}\right)\!\right)
                                  \left(\!\left(\frac{pi}{q}\right)\!\right)

    where

    .. MATH::

        ((x))=\begin{cases}
            x-\lfloor x \rfloor - \frac{1}{2} &\mbox{if }
                x \in \QQ \setminus \ZZ \\
            0 & \mbox{if } x \in \ZZ.
            \end{cases}

    .. WARNING::

        Caution is required as the Dedekind sum sometimes depends on the
        algorithm or is left undefined when `p` and `q` are not coprime.

    INPUT:

    -  ``p``, ``q`` -- integers
    -  ``algorithm`` -- must be one of the following

       -  ``'default'`` - (default) use FLINT
       -  ``'flint'`` - use FLINT
       -  ``'pari'`` - use PARI (gives different results if `p` and `q`
          are not coprime)

    OUTPUT: a rational number

    EXAMPLES:

    Several small values::

        sage: for q in range(10): print [dedekind_sum(p,q) for p in range(q+1)]
        [0]
        [0, 0]
        [0, 0, 0]
        [0, 1/18, -1/18, 0]
        [0, 1/8, 0, -1/8, 0]
        [0, 1/5, 0, 0, -1/5, 0]
        [0, 5/18, 1/18, 0, -1/18, -5/18, 0]
        [0, 5/14, 1/14, -1/14, 1/14, -1/14, -5/14, 0]
        [0, 7/16, 1/8, 1/16, 0, -1/16, -1/8, -7/16, 0]
        [0, 14/27, 4/27, 1/18, -4/27, 4/27, -1/18, -4/27, -14/27, 0]

    Check relations for restricted arguments::

        sage: q = 23; dedekind_sum(1, q); (q-1)*(q-2)/(12*q)
        77/46
        77/46
        sage: p, q = 100, 723    # must be coprime
        sage: dedekind_sum(p, q) + dedekind_sum(q, p)
        31583/86760
        sage: -1/4 + (p/q + q/p + 1/(p*q))/12
        31583/86760

    We check that evaluation works with large input::

        sage: dedekind_sum(3^54 - 1, 2^93 + 1)
        459340694971839990630374299870/29710560942849126597578981379
        sage: dedekind_sum(3^54 - 1, 2^93 + 1, algorithm='pari')
        459340694971839990630374299870/29710560942849126597578981379

    Pari uses a different definition if the inputs are not coprime::

        sage: dedekind_sum(5, 7, algorithm='default')
        -1/14
        sage: dedekind_sum(5, 7, algorithm='flint')
        -1/14
        sage: dedekind_sum(5, 7, algorithm='pari')
        -1/14
        sage: dedekind_sum(6, 8, algorithm='default')
        -1/8
        sage: dedekind_sum(6, 8, algorithm='flint')
        -1/8
        sage: dedekind_sum(6, 8, algorithm='pari')
        -1/24

    REFERENCES:

    .. [Apostol] T. Apostol, Modular functions and Dirichlet series
       in number theory, Springer, 1997 (2nd ed), section 3.7--3.9.

    - :wikipedia:`Dedekind\_sum`
    """
    if algorithm == 'default' or algorithm == 'flint':
        return flint_arith.dedekind_sum(p, q)

    if algorithm == 'pari':
        import sage.interfaces.gp
        x = sage.interfaces.gp.gp('sumdedekind(%s,%s)' % (p, q))
        return sage.rings.rational.Rational(x)

    raise ValueError('unknown algorithm')

