"""
Routines for computing special values of L-functions

- :func:`gamma__exact` -- Exact values of the `\Gamma` function at integers and half-integers
- :func:`zeta__exact` -- Exact values of the Riemann `\zeta` function at critical values
- :func:`quadratic_L_function__exact` -- Exact values of the Dirichlet L-functions of quadratic characters at critical values
- :func:`quadratic_L_function__numerical` -- Numerical values of the Dirichlet L-functions of quadratic characters in the domain of convergence
"""

from sage.combinat.combinat import bernoulli_polynomial
from sage.misc.functional import denominator
from sage.rings.all import RealField
from sage.rings.arith import kronecker_symbol, bernoulli, factorial, fundamental_discriminant
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import is_RealField
from sage.symbolic.constants import pi
from sage.symbolic.pynac import I

# ---------------- The Gamma Function  ------------------

def gamma__exact(n):
    """
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

        sage: gamma__exact(1/2)
        sqrt(pi)
        sage: gamma__exact(3/2)
        1/2*sqrt(pi)
        sage: gamma__exact(5/2)
        3/4*sqrt(pi)
        sage: gamma__exact(7/2)
        15/8*sqrt(pi)

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
    from sage.all import sqrt
    # SANITY CHECK
    if (not n in QQ) or (denominator(n) > 2):
        raise TypeError("you must give an integer or half-integer argument")

    if denominator(n) == 1:
        if n <= 0:
            return infinity
        if n > 0:
            return factorial(n-1)
    else:
        ans = QQ(1)
        while (n != QQ(1)/2):
            if n < 0:
                ans *= QQ(1)/n
                n += 1
            elif n > 0:
                n += -1
                ans *= n

        ans *= sqrt(pi)
        return ans

# ------------- The Riemann Zeta Function  --------------

def zeta__exact(n):
    r"""
    Returns the exact value of the Riemann Zeta function

    The argument must be a critical value, namely either positive even
    or negative odd.

    See for example [Iwasawa]_, p13, Special value of `\zeta(2k)`

    EXAMPLES:

    Let us test the accuracy for negative special values::

        sage: RR = RealField(100)
        sage: for i in range(1,10):
        ...       print "zeta(" + str(1-2*i) + "): ", RR(zeta__exact(1-2*i)) - zeta(RR(1-2*i))
        zeta(-1):  0.00000000000000000000000000000
        zeta(-3):  0.00000000000000000000000000000
        zeta(-5):  0.00000000000000000000000000000
        zeta(-7):  0.00000000000000000000000000000
        zeta(-9):  0.00000000000000000000000000000
        zeta(-11):  0.00000000000000000000000000000
        zeta(-13):  0.00000000000000000000000000000
        zeta(-15):  0.00000000000000000000000000000
        zeta(-17):  0.00000000000000000000000000000

    Let us test the accuracy for positive special values::

        sage: all(abs(RR(zeta__exact(2*i))-zeta(RR(2*i))) < 10**(-28) for i in range(1,10))
        True

    TESTS::

        sage: zeta__exact(5)
        Traceback (most recent call last):
        ...
        TypeError: n must be a critical value (i.e. even > 0 or odd < 0)

    REFERENCES:

    .. [Iwasawa] Iwasawa, *Lectures on p-adic L-functions*
    .. [IreRos] Ireland and Rosen, *A Classical Introduction to Modern Number Theory*
    .. [WashCyc] Washington, *Cyclotomic Fields*
    """
    if n < 0:
        return bernoulli(1-n)/(n-1)
    elif n > 1:
        if (n % 2 == 0):
            return ZZ(-1)**(n/2 + 1) * ZZ(2)**(n-1) * pi**n * bernoulli(n) / factorial(n)
        else:
            raise TypeError("n must be a critical value (i.e. even > 0 or odd < 0)")
    elif n==1:
        return infinity
    elif n==0:
        return -1/2

# ---------- Dirichlet L-functions with quadratic characters ----------

def QuadraticBernoulliNumber(k, d):
    r"""
    Compute `k`-th Bernoulli number for the primitive
    quadratic character associated to `\chi(x) = \left(\frac{d}{x}\right)`.

    EXAMPLES:

    Let us create a list of some odd negative fundamental discriminants::

        sage: test_set = [d for d in range(-163, -3, 4) if is_fundamental_discriminant(d)]

    In general, we have `B_{1, \chi_d} = -2 h/w` for odd negative fundamental
    discriminants::

        sage: all(QuadraticBernoulliNumber(1, d) == -len(BinaryQF_reduced_representatives(d)) for d in test_set)
        True

    REFERENCES:

    - [Iwasawa]_, pp 7-16.
    """
    # Ensure the character is primitive
    d1 = fundamental_discriminant(d)
    f = abs(d1)

    # Make the (usual) k-th Bernoulli polynomial
    x =  PolynomialRing(QQ, 'x').gen()
    bp = bernoulli_polynomial(x, k)

    # Make the k-th quadratic Bernoulli number
    total = sum([kronecker_symbol(d1, i) * bp(i/f)  for i in range(f)])
    total *= (f ** (k-1))

    return total

def quadratic_L_function__exact(n, d):
    r"""
    Returns the exact value of a quadratic twist of the Riemann Zeta function
    by `\chi_d(x) = \left(\frac{d}{x}\right)`.

    The input `n` must be a critical value.

    EXAMPLES::

        sage: quadratic_L_function__exact(1, -4)
        1/4*pi
        sage: quadratic_L_function__exact(-4, -4)
        5/2

    TESTS::

        sage: quadratic_L_function__exact(2, -4)
        Traceback (most recent call last):
        ...
        TypeError: n must be a critical value (i.e. odd > 0 or even <= 0)

    REFERENCES:

    - [Iwasawa]_, pp 16-17, Special values of `L(1-n, \chi)` and `L(n, \chi)`
    - [IreRos]_
    - [WashCyc]_
    """
    from sage.all import SR, sqrt
    if n <= 0:
        return QuadraticBernoulliNumber(1-n,d)/(n-1)
    elif n >= 1:
        # Compute the kind of critical values (p10)
        if kronecker_symbol(fundamental_discriminant(d), -1) == 1:
            delta = 0
        else:
            delta = 1

        # Compute the positive special values (p17)
        if ((n - delta) % 2 == 0):
            f = abs(fundamental_discriminant(d))
            if delta == 0:
                GS = sqrt(f)
            else:
                GS = I * sqrt(f)
            ans = SR(ZZ(-1)**(1+(n-delta)/2))
            ans *= (2*pi/f)**n
            ans *= GS     # Evaluate the Gauss sum here! =0
            ans *= 1/(2 * I**delta)
            ans *= QuadraticBernoulliNumber(n,d)/factorial(n)
            return ans
        else:
            if delta == 0:
                raise TypeError("n must be a critical value (i.e. even > 0 or odd < 0)")
            if delta == 1:
                raise TypeError("n must be a critical value (i.e. odd > 0 or even <= 0)")

def quadratic_L_function__numerical(n, d, num_terms=1000):
    """
    Evaluate the Dirichlet L-function (for quadratic character) numerically
    (in a very naive way).

    EXAMPLES:

    First, let us test several values for a given character::

        sage: RR = RealField(100)
        sage: for i in range(5):
        ...       print "L(" + str(1+2*i) + ", (-4/.)): ", RR(quadratic_L_function__exact(1+2*i, -4)) - quadratic_L_function__numerical(RR(1+2*i),-4, 10000)
        L(1, (-4/.)):  0.000049999999500000024999996962707
        L(3, (-4/.)):  4.99999970000003...e-13
        L(5, (-4/.)):  4.99999922759382...e-21
        L(7, (-4/.)):  ...e-29
        L(9, (-4/.)):  ...e-29

    This procedure fails for negative special values, as the Dirichlet
    series does not converge here::

        sage: quadratic_L_function__numerical(-3,-4, 10000)
        Traceback (most recent call last):
        ...
        ValueError: the Dirichlet series does not converge here

    Test for several characters that the result agrees with the exact
    value, to a given accuracy ::

        sage: for d in range(-20,0):  # long time (2s on sage.math 2014)
        ....:     if abs(RR(quadratic_L_function__numerical(1, d, 10000) - quadratic_L_function__exact(1, d))) > 0.001:
        ....:         print "Oops!  We have a problem at d = ", d, "    exact = ", RR(quadratic_L_function__exact(1, d)), "    numerical = ", RR(quadratic_L_function__numerical(1, d))
    """
    # Set the correct precision if it is given (for n).
    if is_RealField(n.parent()):
        R = n.parent()
    else:
        R = RealField()

    if n < 0:
        raise ValueError('the Dirichlet series does not converge here')

    d1 = fundamental_discriminant(d)
    ans = R(0)
    for i in range(1,num_terms):
        ans += R(kronecker_symbol(d1,i) / R(i)**n)
    return ans
