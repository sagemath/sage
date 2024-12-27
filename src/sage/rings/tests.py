"""
Tests for rings

TESTS::

    sage: K.<x> = FractionField(QQ['x'])
    sage: V.<z> = K[]
    sage: x+z
    z + x
"""

import sage.misc.prandom as random

from sage.misc.random_testing import random_testing


def prime_finite_field():
    """
    Create a random prime finite field with cardinality at most 10^20.

    OUTPUT: a prime finite field

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.prime_finite_field(); K                              # needs sage.rings.finite_rings
        Finite Field of size ...
        sage: K.cardinality().is_prime()                                                # needs sage.rings.finite_rings
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.finite_field_constructor import GF
    return GF(ZZ.random_element(x=2, y=10**20 - 12).next_prime())


def finite_field():
    """
    Create a random finite field with degree at most 20 and prime at most 10^6.

    OUTPUT: a finite field

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.finite_field(); K                                    # needs sage.rings.finite_rings
        Finite Field...of size ...
        sage: K.cardinality().is_prime_power()                                          # needs sage.rings.finite_rings
        True
        sage: while K.cardinality().is_prime():                                         # needs sage.rings.finite_rings
        ....:     K = sage.rings.tests.finite_field()
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.finite_field_constructor import GF
    p = ZZ.random_element(x=2, y=10**6 - 18).next_prime()
    d = ZZ.random_element(x=1, y=20)
    return GF(p**d, 'a')


def small_finite_field():
    """
    Create a random finite field with cardinality at most 2^16.

    OUTPUT: a finite field

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: import sage.rings.tests
        sage: K = sage.rings.tests.small_finite_field(); K
        Finite Field...of size ...
        sage: q = K.cardinality()
        sage: q.is_prime_power()
        True
        sage: q <= 2^16
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.finite_field_constructor import GF
    while True:
        q = ZZ.random_element(x=2, y=2**16)
        if q.is_prime_power():
            return GF(q, 'a')


def integer_mod_ring():
    """
    Return a random ring of integers modulo n with n at most 50000.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: R = sage.rings.tests.integer_mod_ring(); R
        Ring of integers modulo ...
        sage: R.cardinality() <= 50000
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    n = ZZ.random_element(x=2, y=50000)
    return IntegerModRing(n)


def padic_field():
    """
    Return a random `p`-adic field modulo n with p at most 10000
    and precision between 10 and 100.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: sage.rings.tests.padic_field()                                            # needs sage.rings.padics
        ...-adic Field with capped relative precision ...
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.padics.factory import Qp
    prec = ZZ.random_element(x=10, y=100)
    p = ZZ.random_element(x=2, y=10**4 - 30).next_prime()
    return Qp(p, prec)


def quadratic_number_field():
    """
    Return a quadratic extension of QQ.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.quadratic_number_field(); K                          # needs sage.rings.number_field
        Number Field in a with defining polynomial x^2 ... with a = ...
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.number_field.number_field import QuadraticField
    while True:
        d = ZZ.random_element(x=-10**5, y=10**5)
        if not d.is_square():
            return QuadraticField(d, 'a')


def absolute_number_field(maxdeg=10):
    """
    Return an absolute extension of QQ of degree at most 10.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.absolute_number_field(); K                           # needs sage.rings.number_field
        Number Field in a with defining polynomial ...
        sage: K.degree() <= 10                                                          # needs sage.rings.number_field
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.number_field.number_field import NumberField
    R = ZZ['x']
    while True:
        f = R.random_element(degree=ZZ.random_element(x=1, y=maxdeg),
                             x=-100, y=100)
        if f.degree() <= 0:
            continue
        f = f + R.gen()**(f.degree() + 1)  # make monic
        if f.is_irreducible():
            return NumberField(f, 'a')


def relative_number_field(n=2, maxdeg=2):
    """
    Return a tower of at most n extensions each of degree at most maxdeg.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: import sage.rings.tests
        sage: K = sage.rings.tests.relative_number_field(3); K
        Number Field in aaa with defining polynomial x^2 ... over its base field
        sage: K.relative_degree()
        2
        sage: L = K.base_ring()
        sage: L.relative_degree()
        2
        sage: M = L.base_ring()
        sage: M.relative_degree()
        2
        sage: M.base_ring() is QQ
        True

    TESTS:

    Check that :issue:`32117` is fixed::

        sage: set_random_seed(3030)
        sage: from sage.rings.tests import relative_number_field
        sage: _ = relative_number_field(3)                                              # needs sage.rings.number_field
    """
    from sage.rings.integer_ring import ZZ
    K = absolute_number_field(maxdeg)
    n -= 1
    var = 'aa'
    R = ZZ['x']
    R1 = K['x']
    while n >= 1:
        while True:
            f = R.random_element(degree=ZZ.random_element(x=1, y=maxdeg),
                                 x=-100, y=100)
            if f.degree() <= 0:
                continue
            f = f * f.denominator()  # bug trac #4781
            f = f + R.gen()**maxdeg  # make monic
            if R1(f).is_irreducible():
                break
        K = K.extension(f, var)
        R1 = K['x']
        var += 'a'
        n -= 1
    return K


def rings0():
    """
    Return a list of pairs (f, desc), where f is a function that when
    called creates a random ring of a certain representative type
    described by desc.

    RINGS:

    - ZZ
    - QQ
    - ZZ/nZZ
    - GF(p)
    - GF(q)
    - p-adic fields
    - quadratic number fields
    - absolute number fields
    - relative number fields

    EXAMPLES::

        sage: import sage.rings.tests
        sage: type(sage.rings.tests.rings0())
        <... 'list'>
    """
    from sage.rings.integer_ring import IntegerRing
    from sage.rings.rational_field import RationalField

    v = [(IntegerRing, 'ring of integers'),
         (RationalField, 'field of rational numbers'),
         (integer_mod_ring, 'integers modulo n for n at most 50000')]
    try:
        v += [(prime_finite_field, 'a prime finite field with cardinality at most 10^20'),
              (finite_field, 'finite field with degree at most 20 and prime at most 10^6'),
              (small_finite_field, 'finite field with cardinality at most 2^16')]
    except ImportError:
        pass

    try:
        v += [(padic_field, 'a p-adic field')]
    except ImportError:
        pass

    try:
        v += [(quadratic_number_field, 'a quadratic number field'),
              (absolute_number_field, 'an absolute number field of degree at most 10'),
              (relative_number_field, 'a tower of at most 2 extensions each of degree at most 2')]
    except ImportError:
        pass

    return v


def rings1():
    """
    Return an iterator over random rings.

    Return a list of pairs (f, desc), where f is a function that
    outputs a random ring that takes a ring and possibly
    some other data as constructor.

    RINGS:

    - polynomial ring in one variable over a rings0() ring.
    - polynomial ring over a rings1() ring.
    - multivariate polynomials
    - power series rings in one variable over a rings0() ring.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: type(sage.rings.tests.rings0())
        <... 'list'>
    """
    v = rings0()
    X = random_rings(level=0)
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.power_series_ring import PowerSeriesRing
    from sage.rings.integer_ring import ZZ

    v = [(lambda: PolynomialRing(next(X), names='x'),
          'univariate polynomial ring over level 0 ring'),
         (lambda: PowerSeriesRing(next(X), names='x'),
          'univariate power series ring over level 0 ring')]

    try:
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    except ImportError:
        pass
    else:
        v += [(lambda: LaurentPolynomialRing(next(X), names='x'),
               'univariate Laurent polynomial ring over level 0 ring')]

    v += [(lambda: PolynomialRing(next(X), abs(ZZ.random_element(x=2, y=10)),
                                  names='x'),
           'multivariate polynomial ring in between 2 and 10 variables over a level 0 ring')]

    return v


MAX_LEVEL = 99999


def random_rings(level=MAX_LEVEL):
    """
    Return an iterator over random rings up to the given "level" of complexity.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: type(sage.rings.tests.random_rings())
        <... 'generator'>
    """
    v = rings0()
    if level >= 1:
        v += rings1()
    while True:
        yield random.choice(v)[0]()


@random_testing
def check_random_elements(level=MAX_LEVEL, trials=1):
    """
    Create random elements of random rings until a crash occurs, in
    which case an exception is raised.  Defaults to running a single
    trial, but more can be specified.  To run tests in an infinite
    loop, you could use::

        while True: check_random_elements(trials=100, print_seed=True)

    INPUT:

    - ``level`` -- (default: ``MAX_LEVEL``) controls the types of rings to use
    - ``trials`` -- a positive integer (default: 1); the number of trials to run
    - ``seed`` -- the random seed to use; if not specified, uses a truly random seed
    - ``print_seed`` -- if ``True`` (default: ``False``), prints the random seed chosen

    EXAMPLES::

        sage: import sage.rings.tests
        sage: sage.rings.tests.check_random_elements(trials=2, seed=0)                   # needs sage.rings.number_field
        survived 0 tests
        Rational Field
        -1/2
        ----
        survived 1 tests
        Number Field in a with defining polynomial x^2 - 61891 with a = 248.7790184079036?
        -12
        ----

        sage: # needs sage.rings.finite_rings sage.rings.number_field sage.rings.padics
        sage: sage.rings.tests.check_random_elements(trials=10)
        survived 0 tests...
        sage: sage.rings.tests.check_random_elements(trials=1000)  # long time (5 seconds)
        survived 0 tests...
    """
    r = random_rings(level)
    i = 0
    for R in r:
        print("survived %s tests" % i)
        i += 1
        print(R)
        print(R.random_element())
        print("----")
        if i >= trials:
            return


@random_testing
def check_random_arith(level=MAX_LEVEL, trials=1):
    """
    Create random elements of random rings and do some arithmetic with them.

    Repeats until a crash occurs, in which case an exception is
    raised.  Defaults to running a single trial, but more can be
    specified.  To run tests in an infinite loop, you could use::

        while True: check_random_arith(trials=100, print_seed=True)

    INPUT:

    - ``level`` -- (default: ``MAX_LEVEL``) controls the types of rings to use
    - ``trials`` -- positive integer (default: 1); the number of trials to run
    - ``seed`` -- the random seed to use; if not specified, uses a truly random seed
    - ``print_seed`` -- if ``True`` (default: ``False``), prints the random seed chosen

    EXAMPLES::

        sage: # needs sage.rings.finite_rings sage.rings.number_field sage.rings.padics
        sage: import sage.rings.tests
        sage: sage.rings.tests.check_random_arith(trials=2, seed=0)
        survived 0 tests
        Rational Field
        -1/2 -1/95
        49/95
        survived 1 tests
        Number Field in a with defining polynomial x^2 - 15083 with a = 122.81286577553673?
        a -2*a - 1
        2*a - 30164
        sage: sage.rings.tests.check_random_arith(trials=10)
        survived 0 tests...
        sage: sage.rings.tests.check_random_arith(trials=1000)   # long time (5 seconds?)
        survived 0 tests...
    """
    i = 0
    for x in random_rings(level):
        print("survived %s tests" % i)
        i += 1
        print(x)
        a = x.random_element()
        b = x.random_element()
        print(a, b)
        print(a * b + a - b + 1)
        if i >= trials:
            return


@random_testing
def check_karatsuba_multiplication(base_ring, maxdeg1, maxdeg2,
                                  ref_mul=lambda f, g: f._mul_generic(g),
                                  base_ring_random_elt_args=[],
                                  numtests=10, verbose=False):
    """
    Test univariate Karatsuba multiplication against other multiplication algorithms.

    EXAMPLES:

    First check that random tests are reproducible::

        sage: from sage.rings.tests import check_karatsuba_multiplication
        sage: check_karatsuba_multiplication(ZZ, 6, 5, verbose=True, seed=42)
        check_karatsuba_multiplication: ring=Univariate Polynomial Ring in x over Integer Ring, threshold=2
          (x^6 + 4*x^5 + 4*x^4 - 3*x^3 - x^2 - x)*(2*x^4 + 3*x^3 - 20*x^2 - 2*x + 1)
          (4*x^5 + 16*x^2 + x - 41)*(x^2 + x - 1)
          (8*x^2 + 2*x + 1)*(3)
          (-4*x - 1)*(-8*x^2 - x)
          (-x^6 - x^3 - x^2 + x + 1)*(2*x^3 - x + 3)
          (-x^2 + x + 1)*(x^4 + x^3 - x^2 - x + 76)
          (4*x^3 + x^2 + 6)*(-x^2 - 5*x)
          (x + 4)*(-x + 5)
          (-2*x)*(3*x^2 - x)
          (x^6 + 21*x^5 + x^4 + 4*x^3 - x^2)*(14*x^4 + x^3 + 2*x^2 - 12*x)

    Test Karatsuba multiplication of polynomials of small degree over some common rings::

        sage: rings = [QQ]
        sage: rings += [ZZ[I], ZZ[I, sqrt(2)]]                                          # needs sage.rings.number_field sage.symbolic
        sage: rings += [GF(49, 'a')]                                                    # needs sage.rings.finite_rings
        sage: rings += [MatrixSpace(GF(17), 3)]                                         # needs sage.modules
        sage: for C in rings:                                                           # needs sage.modules
        ....:     check_karatsuba_multiplication(C, 10, 10)

    Zero-tests over ``QQbar`` are currently very slow, so we test only very small examples::

        sage: check_karatsuba_multiplication(QQbar, 3, 3, numtests=2)    # long time, needs sage.rings.number_field

    Larger degrees (over ``ZZ``, using FLINT)::

        sage: check_karatsuba_multiplication(ZZ, 1000, 1000,
        ....:                               ref_mul=lambda f,g: f*g,
        ....:                               base_ring_random_elt_args=[1000])

    Some more aggressive tests::

        sage: testrings = [ZZ[I, sqrt(2)], ZZ[I, sqrt(2), sqrt(3)]]     # long time
        sage: for C in testrings:                                       # long time
        ....:     check_karatsuba_multiplication(C, 100, 100)
        sage: check_karatsuba_multiplication(ZZ, 10000, 10000,           # long time
        ....:                               ref_mul=lambda f,g: f*g,
        ....:                               base_ring_random_elt_args=[100000])
    """
    from sage.misc.prandom import randint
    from sage.misc.sage_input import sage_input
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    threshold = randint(0, min(maxdeg1, maxdeg2))
    R = PolynomialRing(base_ring, 'x')
    if verbose:
        print(f"check_karatsuba_multiplication: ring={R}, threshold={threshold}")
    for _ in range(numtests):
        f = R.random_element(randint(0, maxdeg1), False, *base_ring_random_elt_args)
        g = R.random_element(randint(0, maxdeg2), False, *base_ring_random_elt_args)
        if verbose:
            print("  ({})*({})".format(f, g))
        if ref_mul(f, g) - f._mul_karatsuba(g, threshold) != 0:
            msg = "Multiplication failed for elements defined by\n"
            msg += f"{sage_input(f)}\n"
            msg += "and\n"
            msg += f"{sage_input(g)}"
            raise ValueError(msg)
    return
