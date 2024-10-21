# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.libs.pari sage.modules
r"""
Routines for computing special values of `L`-functions

- :func:`gamma__exact` -- exact values of the `\Gamma` function at integers and half-integers
- :func:`zeta__exact` -- exact values of the Riemann `\zeta` function at critical values
- :func:`quadratic_L_function__exact` -- exact values of the Dirichlet `L`-functions of quadratic characters at critical values
- :func:`quadratic_L_function__numerical` -- numerical values of the Dirichlet `L`-functions of quadratic characters in the domain of convergence
"""

import sage.rings.abc

from sage.arith.misc import (bernoulli,
                             factorial,
                             fundamental_discriminant,
                             kronecker as kronecker_symbol)
from sage.misc.functional import denominator
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ


# ---------------- The Gamma Function  ------------------

def gamma__exact(n):
    r"""
    Evaluates the exact value of the `\Gamma` function at an integer or
    half-integer argument.

    EXAMPLES::

        sage: gamma__exact(4)
        6
        sage: gamma__exact(3)
        2
        sage: gamma__exact(2)
        1
        sage: gamma__exact(1)
        1

        sage: # needs sage.symbolic
        sage: gamma__exact(1/2)
        sqrt(pi)
        sage: gamma__exact(3/2)
        1/2*sqrt(pi)
        sage: gamma__exact(5/2)
        3/4*sqrt(pi)
        sage: gamma__exact(7/2)
        15/8*sqrt(pi)

        sage: # needs sage.symbolic
        sage: gamma__exact(-1/2)
        -2*sqrt(pi)
        sage: gamma__exact(-3/2)
        4/3*sqrt(pi)
        sage: gamma__exact(-5/2)
        -8/15*sqrt(pi)
        sage: gamma__exact(-7/2)
        16/105*sqrt(pi)

    TESTS::

        sage: gamma__exact(1/3)
        Traceback (most recent call last):
        ...
        TypeError: you must give an integer or half-integer argument
    """
    n = QQ(n)

    if denominator(n) == 1:
        if n <= 0:
            return infinity
        return factorial(n - 1)

    if denominator(n) == 2:
        # now n = 1/2 + an integer
        from sage.misc.functional import sqrt
        from sage.symbolic.constants import pi

        ans = QQ.one()
        while n != QQ((1, 2)):
            if n < 0:
                ans /= n
                n += 1
            else:
                n += -1
                ans *= n

        ans *= sqrt(pi)
        return ans

    raise TypeError("you must give an integer or half-integer argument")


# ------------- The Riemann Zeta Function  --------------

def zeta__exact(n):
    r"""
    Return the exact value of the Riemann Zeta function.

    The argument must be a critical value, namely either positive even
    or negative odd.

    See for example [Iwa1972]_, p13, Special value of `\zeta(2k)`

    EXAMPLES:

    Let us test the accuracy for negative special values::

        sage: RR = RealField(100)
        sage: for i in range(1,10):                                                     # needs sage.symbolic
        ....:     print("zeta({}): {}".format(1 - 2*i,
        ....:                                 RR(zeta__exact(1-2*i)) - zeta(RR(1-2*i))))
        zeta(-1): 0.00000000000000000000000000000
        zeta(-3): 0.00000000000000000000000000000
        zeta(-5): 0.00000000000000000000000000000
        zeta(-7): 0.00000000000000000000000000000
        zeta(-9): 0.00000000000000000000000000000
        zeta(-11): 0.00000000000000000000000000000
        zeta(-13): 0.00000000000000000000000000000
        zeta(-15): 0.00000000000000000000000000000
        zeta(-17): 0.00000000000000000000000000000

    Let us test the accuracy for positive special values::

        sage: all(abs(RR(zeta__exact(2*i)) - zeta(RR(2*i))) < 10**(-28)                 # needs sage.symbolic
        ....:     for i in range(1,10))
        True

    TESTS::

        sage: zeta__exact(4)                                                            # needs sage.symbolic
        1/90*pi^4
        sage: zeta__exact(-3)
        1/120
        sage: zeta__exact(0)
        -1/2
        sage: zeta__exact(5)
        Traceback (most recent call last):
        ...
        TypeError: n must be a critical value (i.e. even > 0 or odd < 0)

    REFERENCES:

    - [Iwa1972]_
    - [IR1990]_
    - [Was1997]_
    """
    if n < 0:
        return bernoulli(1-n)/(n-1)
    elif n > 1:
        if (n % 2 == 0):
            from sage.symbolic.constants import pi

            return ZZ(-1)**(n//2 + 1) * ZZ(2)**(n-1) * pi**n * bernoulli(n) / factorial(n)
        else:
            raise TypeError("n must be a critical value (i.e. even > 0 or odd < 0)")
    elif n == 1:
        return infinity
    elif n == 0:
        return QQ((-1, 2))


# ---------- Dirichlet L-functions with quadratic characters ----------

def QuadraticBernoulliNumber(k, d):
    r"""
    Compute `k`-th Bernoulli number for the primitive
    quadratic character associated to `\chi(x) = \left(\frac{d}{x}\right)`.

    EXAMPLES:

    Let us create a list of some odd negative fundamental discriminants::

        sage: test_set = [d for d in srange(-163, -3, 4)
        ....:             if d.is_fundamental_discriminant()]

    In general, we have `B_{1, \chi_d} = -2 h/w` for odd negative fundamental
    discriminants::

        sage: all(QuadraticBernoulliNumber(1, d)
        ....:       == -len(BinaryQF_reduced_representatives(d))
        ....:     for d in test_set)
        True

    REFERENCES:

    - [Iwa1972]_, pp 7-16.
    """
    from sage.combinat.combinat import bernoulli_polynomial

    # Ensure the character is primitive
    d1 = fundamental_discriminant(d)
    f = abs(d1)

    # Make the (usual) k-th Bernoulli polynomial
    x = PolynomialRing(QQ, 'x').gen()
    bp = bernoulli_polynomial(x, k)

    # Make the k-th quadratic Bernoulli number
    total = sum([kronecker_symbol(d1, i) * bp(i / f) for i in range(f)])
    total *= f**(k - 1)

    return total


def quadratic_L_function__exact(n, d):
    r"""
    Return the exact value of a quadratic twist of the Riemann Zeta function
    by `\chi_d(x) = \left(\frac{d}{x}\right)`.

    The input `n` must be a critical value.

    EXAMPLES::

        sage: quadratic_L_function__exact(1, -4)                                        # needs sage.symbolic
        1/4*pi
        sage: quadratic_L_function__exact(-4, -4)
        5/2
        sage: quadratic_L_function__exact(2, 1)                                         # needs sage.symbolic
        1/6*pi^2

    TESTS::

        sage: quadratic_L_function__exact(2, -4)
        Traceback (most recent call last):
        ...
        TypeError: n must be a critical value (i.e. odd > 0 or even <= 0)

    REFERENCES:

    - [Iwa1972]_, pp 16-17, Special values of `L(1-n, \chi)` and `L(n, \chi)`
    - [IR1990]_
    - [Was1997]_
    """
    if n <= 0:
        return QuadraticBernoulliNumber(1-n, d)/(n-1)
    elif n >= 1:
        # Compute the kind of critical values (p10)
        if kronecker_symbol(fundamental_discriminant(d), -1) == 1:
            delta = 0
        else:
            delta = 1

        # Compute the positive special values (p17)
        if ((n - delta) % 2 == 0):
            from sage.misc.functional import sqrt
            from sage.symbolic.constants import I, pi
            from sage.symbolic.ring import SR

            f = abs(fundamental_discriminant(d))
            if delta == 0:
                GS = sqrt(f)
            else:
                GS = I * sqrt(f)
            ans = SR(ZZ(-1)**(1+(n-delta)/2))
            ans *= (2*pi/f)**n
            ans *= GS     # Evaluate the Gauss sum here! =0
            ans *= QQ.one()/(2 * I**delta)
            ans *= QuadraticBernoulliNumber(n, d)/factorial(n)
            return ans
        else:
            if delta == 0:
                raise TypeError("n must be a critical value (i.e. even > 0 or odd < 0)")
            if delta == 1:
                raise TypeError("n must be a critical value (i.e. odd > 0 or even <= 0)")


def quadratic_L_function__numerical(n, d, num_terms=1000):
    """
    Evaluate the Dirichlet `L`-function (for quadratic character) numerically
    (in a very naive way).

    EXAMPLES:

    First, let us test several values for a given character::

        sage: RR = RealField(100)
        sage: for i in range(5):                                                        # needs sage.symbolic
        ....:     print("L({}, (-4/.)): {}".format(1+2*i,
        ....:             RR(quadratic_L_function__exact(1+2*i, -4))
        ....:                - quadratic_L_function__numerical(RR(1+2*i), -4, 10000)))
        L(1, (-4/.)): 0.000049999999500000024999996962707
        L(3, (-4/.)): 4.99999970000003...e-13
        L(5, (-4/.)): 4.99999922759382...e-21
        L(7, (-4/.)): ...e-29
        L(9, (-4/.)): ...e-29

    This procedure fails for negative special values, as the Dirichlet
    series does not converge here::

        sage: quadratic_L_function__numerical(-3, -4, 10000)
        Traceback (most recent call last):
        ...
        ValueError: the Dirichlet series does not converge here

    Test for several characters that the result agrees with the exact
    value, to a given accuracy ::

        sage: for d in range(-20,0):            # long time (2s on sage.math 2014), needs sage.symbolic
        ....:     if abs(RR(quadratic_L_function__numerical(1, d, 10000)
        ....:                - quadratic_L_function__exact(1, d))) > 0.001:
        ....:         print("We have a problem at d = {}: exact = {}, numerical = {}".format(d,
        ....:                   RR(quadratic_L_function__exact(1, d)),
        ....:                   RR(quadratic_L_function__numerical(1, d))))
    """
    # Set the correct precision if it is given (for n).
    if isinstance(n.parent(), sage.rings.abc.RealField):
        R = n.parent()
    else:
        from sage.rings.real_mpfr import RealField
        R = RealField()

    if n < 0:
        raise ValueError('the Dirichlet series does not converge here')

    d1 = fundamental_discriminant(d)
    ans = R.zero()
    for i in range(1, num_terms):
        ans += R(kronecker_symbol(d1, i) / R(i)**n)
    return ans
