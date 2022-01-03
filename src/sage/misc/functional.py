# -*- coding: utf-8 -*-
"""
Functional notation

These are functions so that you can write foo(x) instead of x.foo()
in certain common cases.

AUTHORS:

- William Stein: Initial version

- David Joyner (2005-12-20): More Examples
"""
# ****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import builtins
import math

from sage.rings.complex_double import CDF
from sage.rings.real_double import RDF, RealDoubleElement
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.misc.superseded import deprecation

##############################################################################
# There are many functions on elements of a ring, which mathematicians
# usually write f(x), e.g., it is weird to write x.log() and natural
# to write log(x).  The functions below allow for the more familiar syntax.
##############################################################################


def additive_order(x):
    """
    Return the additive order of ``x``.

    EXAMPLES::

        sage: additive_order(5)
        +Infinity
        sage: additive_order(Mod(5,11))
        11
        sage: additive_order(Mod(4,12))
        3
    """
    return x.additive_order()


def base_ring(x):
    """
    Return the base ring over which ``x`` is defined.

    EXAMPLES::

        sage: R = PolynomialRing(GF(7), 'x')
        sage: base_ring(R)
        Finite Field of size 7
    """
    return x.base_ring()


def base_field(x):
    """
    Return the base field over which ``x`` is defined.

    EXAMPLES::

        sage: R = PolynomialRing(GF(7), 'x')
        sage: base_ring(R)
        Finite Field of size 7
        sage: base_field(R)
        Finite Field of size 7

    This catches base rings which are fields as well, but does
    not implement a ``base_field`` method for objects which do
    not have one::

        sage: R.base_field()
        Traceback (most recent call last):
        ...
        AttributeError: 'PolynomialRing_dense_mod_p_with_category' object has no attribute 'base_field'
    """
    try:
        return x.base_field()
    except AttributeError:
        y = x.base_ring()
        if y.is_field():
            return y
        else:
            raise AttributeError("The base ring of %s is not a field" % x)


def basis(x):
    """
    Return the fixed basis of ``x``.

    EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: S = V.subspace([[1,2,0],[2,2,-1]])
        sage: basis(S)
        [
        (1, 0, -1),
        (0, 1, 1/2)
        ]
    """
    return x.basis()


def category(x):
    """
    Return the category of ``x``.

    EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: category(V)
        Category of finite dimensional vector spaces with basis over
         (number fields and quotient fields and metric spaces)
    """
    try:
        return x.category()
    except AttributeError:
        from sage.categories.objects import Objects
        return Objects()


def characteristic_polynomial(x, var='x'):
    """
    Return the characteristic polynomial of ``x`` in the given variable.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: charpoly(A)
        x^3 - 15*x^2 - 18*x
        sage: charpoly(A, 't')
        t^3 - 15*t^2 - 18*t

        sage: k.<alpha> = GF(7^10); k
        Finite Field in alpha of size 7^10
        sage: alpha.charpoly('T')
        T^10 + T^6 + T^5 + 4*T^4 + T^3 + 2*T^2 + 3*T + 3
        sage: characteristic_polynomial(alpha, 'T')
        T^10 + T^6 + T^5 + 4*T^4 + T^3 + 2*T^2 + 3*T + 3

    Ensure the variable name of the polynomial does not conflict with
    variables used within the matrix, and that non-integral powers of
    variables do not confuse the computation (:trac:`14403`)::

        sage: y = var('y')
        sage: a = matrix([[x,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        sage: characteristic_polynomial(a).list()
        [x, -3*x - 1, 3*x + 3, -x - 3, 1]
        sage: b = matrix([[y^(1/2),0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        sage: charpoly(b).list()
        [sqrt(y), -3*sqrt(y) - 1, 3*sqrt(y) + 3, -sqrt(y) - 3, 1]
    """
    try:
        return x.charpoly(var)
    except AttributeError:
        raise NotImplementedError("computation of charpoly of M (={}) not implemented".format(x))


charpoly = characteristic_polynomial


def coerce(P, x):
    """
    Coerce ``x`` to type ``P`` if possible.

    EXAMPLES::

        sage: type(5)
        <class 'sage.rings.integer.Integer'>
        sage: type(coerce(QQ,5))
        <class 'sage.rings.rational.Rational'>
    """
    try:
        return P._coerce_(x)
    except AttributeError:
        return P(x)


def cyclotomic_polynomial(n, var='x'):
    """
    Return the `n^{th}` cyclotomic polynomial.

    EXAMPLES::

        sage: cyclotomic_polynomial(3)
        x^2 + x + 1
        sage: cyclotomic_polynomial(4)
        x^2 + 1
        sage: cyclotomic_polynomial(9)
        x^6 + x^3 + 1
        sage: cyclotomic_polynomial(10)
        x^4 - x^3 + x^2 - x + 1
        sage: cyclotomic_polynomial(11)
        x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
    """
    return ZZ[var].cyclotomic_polynomial(n)


def decomposition(x):
    """
    Return the decomposition of ``x``.

    EXAMPLES::

        sage: M = matrix([[2, 3], [3, 4]])
        sage: M.decomposition()
        [
        (Ambient free module of rank 2 over the principal ideal domain Integer Ring, True)
        ]

        sage: G.<a,b> = DirichletGroup(20)
        sage: c = a*b
        sage: d = c.decomposition(); d
        [Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1,
        Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4]
        sage: d[0].parent()
        Group of Dirichlet characters modulo 4 with values in Cyclotomic Field of order 4 and degree 2
    """
    return x.decomposition()


def denominator(x):
    """
    Return the denominator of ``x``.

    EXAMPLES::

        sage: denominator(17/11111)
        11111
        sage: R.<x> = PolynomialRing(QQ)
        sage: F = FractionField(R)
        sage: r = (x+1)/(x-1)
        sage: denominator(r)
        x - 1
    """
    if isinstance(x, int):
        return 1
    return x.denominator()


def det(x):
    """
    Return the determinant of ``x``.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: det(A)
        0
    """
    return x.det()


def dimension(x):
    """
    Return the dimension of ``x``.

    EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: S = V.subspace([[1,2,0],[2,2,-1]])
        sage: dimension(S)
        2
    """
    return x.dimension()


dim = dimension


def discriminant(x):
    """
    Return the discriminant of ``x``.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: S = R.quotient(x^29 - 17*x - 1, 'alpha')
        sage: K = S.number_field()
        sage: discriminant(K)
        -15975100446626038280218213241591829458737190477345113376757479850566957249523
    """
    return x.discriminant()


disc = discriminant


def eta(x):
    r"""
    Return the value of the `\eta` function at ``x``, which must be
    in the upper half plane.

    The `\eta` function is

    .. MATH::

        \eta(z) = e^{\pi i z / 12} \prod_{n=1}^{\infty}(1-e^{2\pi inz})

    EXAMPLES::

        sage: eta(1+I)
        0.7420487758365647 + 0.1988313702299107*I
    """
    try:
        return x.eta()
    except AttributeError:
        return CDF(x).eta()


def fcp(x, var='x'):
    """
    Return the factorization of the characteristic polynomial of ``x``.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: fcp(A, 'x')
        x * (x^2 - 15*x - 18)
    """
    try:
        return x.fcp(var)
    except AttributeError:
        return charpoly(x, var).factor()


def gen(x):
    """
    Return the generator of ``x``.

    EXAMPLES::

        sage: R.<x> = QQ[]; R
        Univariate Polynomial Ring in x over Rational Field
        sage: gen(R)
        x
        sage: gen(GF(7))
        1
        sage: A = AbelianGroup(1, [23])
        sage: gen(A)
        f
    """
    return x.gen()


def gens(x):
    """
    Return the generators of ``x``.

    EXAMPLES::

        sage: R.<x,y> = SR[]
        sage: R
        Multivariate Polynomial Ring in x, y over Symbolic Ring
        sage: gens(R)
        (x, y)
        sage: A = AbelianGroup(5, [5,5,7,8,9])
        sage: gens(A)
        (f0, f1, f2, f3, f4)
    """
    return x.gens()


def hecke_operator(x, n):
    r"""
    Return the `n`-th Hecke operator `T_n` acting on ``x``.

    EXAMPLES::

        sage: M = ModularSymbols(1,12)
        sage: hecke_operator(M,5)
        Hecke operator T_5 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
    """
    return x.hecke_operator(n)


def image(x):
    """
    Return the image of ``x``.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: image(A)
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1]
        [ 0  1  2]
    """
    return x.image()


def symbolic_sum(expression, *args, **kwds):
    r"""
    Return the symbolic sum `\sum_{v = a}^b expression` with respect
    to the variable `v` with endpoints `a` and `b`.

    INPUT:

    - ``expression`` - a symbolic expression

    - ``v`` - a variable or variable name

    - ``a`` - lower endpoint of the sum

    - ``b`` - upper endpoint of the sum

    - ``algorithm`` - (default: ``'maxima'``)  one of

      - ``'maxima'`` - use Maxima (the default)

      - ``'maple'`` - (optional) use Maple

      - ``'mathematica'`` - (optional) use Mathematica

      - ``'giac'`` - (optional) use Giac

      - ``'sympy'`` - use SymPy

    EXAMPLES::

        sage: k, n = var('k,n')
        sage: sum(k, k, 1, n).factor()
        1/2*(n + 1)*n

    ::

        sage: sum(1/k^4, k, 1, oo)
        1/90*pi^4

    ::

        sage: sum(1/k^5, k, 1, oo)
        zeta(5)

    .. WARNING::

        This function only works with symbolic expressions. To sum any
        other objects like list elements or function return values,
        please use python summation, see
        http://docs.python.org/library/functions.html#sum

        In particular, this does not work::

            sage: n = var('n')
            sage: mylist = [1,2,3,4,5]
            sage: sum(mylist[n], n, 0, 3)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert n to an integer

        Use python ``sum()`` instead::

            sage: sum(mylist[n] for n in range(4))
            10

        Also, only a limited number of functions are recognized in symbolic sums::

            sage: sum(valuation(n,2),n,1,5)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert n to an integer

        Again, use python ``sum()``::

            sage: sum(valuation(n+1,2) for n in range(5))
            3

        (now back to the Sage ``sum`` examples)

    A well known binomial identity::

        sage: sum(binomial(n,k), k, 0, n)
        2^n

    The binomial theorem::

        sage: x, y = var('x, y')
        sage: sum(binomial(n,k) * x^k * y^(n-k), k, 0, n)
        (x + y)^n

    ::

        sage: sum(k * binomial(n, k), k, 1, n)
        2^(n - 1)*n

    ::

        sage: sum((-1)^k*binomial(n,k), k, 0, n)
        0

    ::

        sage: sum(2^(-k)/(k*(k+1)), k, 1, oo)
        -log(2) + 1

    Another binomial identity (:trac:`7952`)::

        sage: t,k,i = var('t,k,i')
        sage: sum(binomial(i+t,t),i,0,k)
        binomial(k + t + 1, t + 1)

    Summing a hypergeometric term::

        sage: sum(binomial(n, k) * factorial(k) / factorial(n+1+k), k, 0, n)
        1/2*sqrt(pi)/factorial(n + 1/2)

    We check a well known identity::

        sage: bool(sum(k^3, k, 1, n) == sum(k, k, 1, n)^2)
        True

    A geometric sum::

        sage: a, q = var('a, q')
        sage: sum(a*q^k, k, 0, n)
        (a*q^(n + 1) - a)/(q - 1)

    The geometric series::

        sage: assume(abs(q) < 1)
        sage: sum(a*q^k, k, 0, oo)
        -a/(q - 1)

    A divergent geometric series.  Don't forget
    to forget your assumptions::

        sage: forget()
        sage: assume(q > 1)
        sage: sum(a*q^k, k, 0, oo)
        Traceback (most recent call last):
        ...
        ValueError: Sum is divergent.

    This summation only Mathematica can perform::

        sage: sum(1/(1+k^2), k, -oo, oo, algorithm = 'mathematica')     # optional - mathematica
        pi*coth(pi)

    Use Maple as a backend for summation::

        sage: sum(binomial(n,k)*x^k, k, 0, n, algorithm = 'maple')      # optional - maple
        (x + 1)^n

    Python ints should work as limits of summation (:trac:`9393`)::

        sage: sum(x, x, 1r, 5r)
        15

    .. note::

       #. Sage can currently only understand a subset of the output of Maxima, Maple and
          Mathematica, so even if the chosen backend can perform the summation the
          result might not be convertible into a Sage expression.

    """
    if hasattr(expression, 'sum'):
        return expression.sum(*args, **kwds)
    elif len(args) <= 1:
        return sum(expression, *args)
    else:
        from sage.symbolic.ring import SR
        return SR(expression).sum(*args, **kwds)


def symbolic_prod(expression, *args, **kwds):
    r"""
    Return the symbolic product `\prod_{v = a}^b expression` with respect
    to the variable `v` with endpoints `a` and `b`.

    INPUT:

    - ``expression`` - a symbolic expression

    - ``v`` - a variable or variable name

    - ``a`` - lower endpoint of the product

    - ``b`` - upper endpoint of the prduct

    - ``algorithm`` - (default: ``'maxima'``)  one of

      - ``'maxima'`` - use Maxima (the default)

      - ``'giac'`` - (optional) use Giac

      - ``'sympy'`` - use SymPy

    - ``hold`` - (default: ``False``) if ``True`` don't evaluate

    EXAMPLES::

        sage: i, k, n = var('i,k,n')
        sage: product(k,k,1,n)
        factorial(n)
        sage: product(x + i*(i+1)/2, i, 1, 4)
        x^4 + 20*x^3 + 127*x^2 + 288*x + 180
        sage: product(i^2, i, 1, 7)
        25401600
        sage: f = function('f')
        sage: product(f(i), i, 1, 7)
        f(7)*f(6)*f(5)*f(4)*f(3)*f(2)*f(1)
        sage: product(f(i), i, 1, n)
        product(f(i), i, 1, n)
        sage: assume(k>0)
        sage: product(integrate (x^k, x, 0, 1), k, 1, n)
        1/factorial(n + 1)
        sage: product(f(i), i, 1, n).log().log_expand()
        sum(log(f(i)), i, 1, n)

    """
    from .misc_c import prod as c_prod
    if hasattr(expression, 'prod'):
        return expression.prod(*args, **kwds)
    elif len(args) <= 1:
        return c_prod(expression, *args)
    else:
        from sage.symbolic.ring import SR
        return SR(expression).prod(*args, **kwds)


def integral(x, *args, **kwds):
    """
    Return an indefinite or definite integral of an object ``x``.

    First call ``x.integral()`` and if that fails make an object and
    integrate it using Maxima, maple, etc, as specified by algorithm.

    For symbolic expression calls
    :func:`sage.calculus.calculus.integral` - see this function for
    available options.

    EXAMPLES::

        sage: f = cyclotomic_polynomial(10)
        sage: integral(f)
        1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x

    ::

        sage: integral(sin(x),x)
        -cos(x)

    ::

        sage: y = var('y')
        sage: integral(sin(x),y)
        y*sin(x)

    ::

        sage: integral(sin(x), x, 0, pi/2)
        1
        sage: sin(x).integral(x, 0,pi/2)
        1
        sage: integral(exp(-x), (x, 1, oo))
        e^(-1)

    Numerical approximation::

        sage: h = integral(tan(x)/x, (x, 1, pi/3))
        ...
        sage: h
        integrate(tan(x)/x, x, 1, 1/3*pi)
        sage: h.n()
        0.07571599101...

    Specific algorithm can be used for integration::

        sage: integral(sin(x)^2, x, algorithm='maxima')
        1/2*x - 1/4*sin(2*x)
        sage: integral(sin(x)^2, x, algorithm='sympy')
        -1/2*cos(x)*sin(x) + 1/2*x

    TESTS:

    A symbolic integral from :trac:`11445` that was incorrect in
    earlier versions of Maxima::

        sage: f = abs(x - 1) + abs(x + 1) - 2*abs(x)
        sage: integrate(f, (x, -Infinity, Infinity))
        2

    Another symbolic integral, from :trac:`11238`, that used to return
    zero incorrectly; with Maxima 5.26.0 one gets
    ``1/2*sqrt(pi)*e^(1/4)``, whereas with 5.29.1, and even more so
    with 5.33.0, the expression is less pleasant, but still has the
    same value.  Unfortunately, the computation takes a very long time
    with the default settings, so we temporarily use the Maxima
    setting ``domain: real``::

        sage: sage.calculus.calculus.maxima('domain: real')
        real
        sage: f = exp(-x) * sinh(sqrt(x))
        sage: t = integrate(f, x, 0, Infinity); t            # long time
        1/4*sqrt(pi)*(erf(1) - 1)*e^(1/4) - 1/4*(sqrt(pi)*(erf(1) - 1) - sqrt(pi) + 2*e^(-1) - 2)*e^(1/4) + 1/4*sqrt(pi)*e^(1/4) - 1/2*e^(1/4) + 1/2*e^(-3/4)
        sage: t.canonicalize_radical()  # long time
        1/2*sqrt(pi)*e^(1/4)
        sage: sage.calculus.calculus.maxima('domain: complex')
        complex

    An integral which used to return -1 before maxima 5.28. See :trac:`12842`::

        sage: f = e^(-2*x)/sqrt(1-e^(-2*x))
        sage: integrate(f, x, 0, infinity)
        1

    This integral would cause a stack overflow in earlier versions of
    Maxima, crashing sage. See :trac:`12377`. We don't care about the
    result here, just that the computation completes successfully::

        sage: y = (x^2)*exp(x) / (1 + exp(x))^2
        sage: _ = integrate(y, x, -1000, 1000)

    When SymPy cannot solve an integral it gives it back, so we must
    be able to convert SymPy's ``Integral`` (:trac:`14723`)::

        sage: x, y, z = var('x,y,z')
        sage: f = function('f')
        sage: integrate(f(x), x, algorithm='sympy')
        integrate(f(x), x)
        sage: integrate(f(x), x, 0, 1,algorithm='sympy')
        integrate(f(x), x, 0, 1)
        sage: integrate(integrate(integrate(f(x,y,z), x, algorithm='sympy'), y, algorithm='sympy'), z, algorithm='sympy')
        integrate(integrate(integrate(f(x, y, z), x), y), z)
        sage: integrate(sin(x)*tan(x)/(1-cos(x)), x, algorithm='sympy')
        -integrate(sin(x)*tan(x)/(cos(x) - 1), x)
        sage: _ = var('a,b,x')
        sage: integrate(sin(x)*tan(x)/(1-cos(x)), x, a, b, algorithm='sympy')
        -integrate(sin(x)*tan(x)/(cos(x) - 1), x, a, b)
        sage: import sympy
        sage: x, y, z = sympy.symbols('x y z')
        sage: f = sympy.Function('f')
        sage: SR(sympy.Integral(f(x,y,z), x, y, z))
        integrate(integrate(integrate(f(x, y, z), x), y), z)

    Ensure that the following integral containing a signum term from
    :trac:`11590` can be integrated::

        sage: x = SR.symbol('x', domain='real')
        sage: result = integrate(x * sgn(x^2 - 1/4), x, -1, 0)
        ...
        sage: result
        -1/4

    """
    if hasattr(x, 'integral'):
        return x.integral(*args, **kwds)
    else:
        from sage.symbolic.ring import SR
        return SR(x).integral(*args, **kwds)


integrate = integral


def integral_closure(x):
    """
    Return the integral closure of ``x``.

    EXAMPLES::

        sage: integral_closure(QQ)
        Rational Field
        sage: K.<a> = QuadraticField(5)
        sage: O2 = K.order(2*a); O2
        Order in Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
        sage: integral_closure(O2)
        Maximal Order in Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
    """
    return x.integral_closure()


def interval(a, b):
    r"""
    Integers between `a` and `b` *inclusive* (`a` and `b` integers).

    EXAMPLES::

        sage: I = interval(1,3)
        sage: 2 in I
        True
        sage: 1 in I
        True
        sage: 4 in I
        False
    """
    return list(range(a, b + 1))


def xinterval(a, b):
    r"""
    Iterator over the integers between `a` and `b`, *inclusive*.

    EXAMPLES::

        sage: I = xinterval(2,5); I
        range(2, 6)
        sage: 5 in I
        True
        sage: 6 in I
        False
    """
    return range(a, b + 1)


def is_commutative(x):
    """
    Return whether or not ``x`` is commutative.

    EXAMPLES::

        sage: R = PolynomialRing(QQ, 'x')
        sage: is_commutative(R)
        doctest:...DeprecationWarning: use X.is_commutative() or X in Rings().Commutative()
        See https://trac.sagemath.org/32347 for details.
        True
    """
    deprecation(32347, "use X.is_commutative() or X in Rings().Commutative()")
    return x.is_commutative()


def is_even(x):
    """
    Return whether or not an integer ``x`` is even, e.g., divisible by 2.

    EXAMPLES::

        sage: is_even(-1)
        False
        sage: is_even(4)
        True
        sage: is_even(-2)
        True
    """
    try:
        return x.is_even()
    except AttributeError:
        return x % 2 == 0


def is_integrally_closed(x):
    """
    Return whether ``x`` is integrally closed.

    EXAMPLES::

        sage: is_integrally_closed(QQ)
        doctest:...DeprecationWarning: use X.is_integrally_closed()
        See https://trac.sagemath.org/32347 for details.
        True
        sage: K.<a> = NumberField(x^2 + 189*x + 394)
        sage: R = K.order(2*a)
        sage: is_integrally_closed(R)
        False
    """
    deprecation(32347, "use X.is_integrally_closed()")
    return x.is_integrally_closed()


def is_field(x, proof=True):
    """
    Return whether or not ``x`` is a field.

    Alternatively, one can use ``x in Fields()``.

    EXAMPLES::

        sage: R = PolynomialRing(QQ, 'x')
        sage: F = FractionField(R)
        sage: is_field(F)
        doctest:...DeprecationWarning: use X.is_field() or X in Fields()
        See https://trac.sagemath.org/32347 for details.
        True
    """
    deprecation(32347, "use X.is_field() or X in Fields()")
    return x.is_field(proof=proof)


def is_odd(x):
    """
    Return whether or not ``x`` is odd.

    This is by definition the complement of :func:`is_even`.

    EXAMPLES::

        sage: is_odd(-2)
        False
        sage: is_odd(-3)
        True
        sage: is_odd(0)
        False
        sage: is_odd(1)
        True
    """
    return not is_even(x)


def kernel(x):
    """
    Return the left kernel of ``x``.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,2)
        sage: A = M([1,2,3,4,5,6])
        sage: kernel(A)
        Vector space of degree 3 and dimension 1 over Rational Field
        Basis matrix:
        [ 1 -2  1]
        sage: kernel(A.transpose())
        Vector space of degree 2 and dimension 0 over Rational Field
        Basis matrix:
        []

    Here are two corner cases::

        sage: M = MatrixSpace(QQ,0,3)
        sage: A = M([])
        sage: kernel(A)
        Vector space of degree 0 and dimension 0 over Rational Field
        Basis matrix:
        []
        sage: kernel(A.transpose()).basis()
        [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
    """
    return x.kernel()


def krull_dimension(x):
    """
    Return the Krull dimension of ``x``.

    EXAMPLES::

        sage: krull_dimension(QQ)
        0
        sage: krull_dimension(ZZ)
        1
        sage: krull_dimension(ZZ[sqrt(5)])
        1
        sage: U.<x,y,z> = PolynomialRing(ZZ,3); U
        Multivariate Polynomial Ring in x, y, z over Integer Ring
        sage: U.krull_dimension()
        4
    """
    return x.krull_dimension()


def lift(x):
    """
    Lift an object of a quotient ring `R/I` to `R`.

    EXAMPLES:

    We lift an integer modulo `3`::

        sage: Mod(2,3).lift()
        2

    We lift an element of a quotient polynomial ring::

        sage: R.<x> = QQ['x']
        sage: S.<xmod> = R.quo(x^2 + 1)
        sage: lift(xmod-7)
        x - 7
    """
    try:
        return x.lift()
    except AttributeError:
        raise ArithmeticError("no lift defined.")


def log(*args, **kwds):
    """
    Return the logarithm of the first argument to the base of
    the second argument which if missing defaults to ``e``.

    It calls the ``log`` method of the first argument when computing
    the logarithm, thus allowing the use of logarithm on any object
    containing a ``log`` method. In other words, ``log`` works
    on more than just real numbers.

    .. NOTE::

        In Magma, the order of arguments is reversed from in Sage,
        i.e., the base is given first. We use the opposite ordering, so
        the base can be viewed as an optional second argument.

    EXAMPLES::

        sage: log(e^2)
        2

    To change the base of the logarithm, add a second parameter::

        sage: log(1000,10)
        3

    The synonym ``ln`` can only take one argument::

        sage: ln(RDF(10))
        2.302585092994046
        sage: ln(2.718)
        0.999896315728952
        sage: ln(2.0)
        0.693147180559945
        sage: ln(float(-1))
        3.141592653589793j
        sage: ln(complex(-1))
        3.141592653589793j

    You can use
    :class:`RDF<sage.rings.real_double.RealDoubleField_class>`,
    :class:`~sage.rings.real_mpfr.RealField` or ``n`` to get a
    numerical real approximation::

        sage: log(1024, 2)
        10
        sage: RDF(log(1024, 2))
        10.0
        sage: log(10, 4)
        1/2*log(10)/log(2)
        sage: RDF(log(10, 4))
        1.6609640474436813
        sage: log(10, 2)
        log(10)/log(2)
        sage: n(log(10, 2))
        3.32192809488736
        sage: log(10, e)
        log(10)
        sage: n(log(10, e))
        2.30258509299405

    The log function works for negative numbers, complex
    numbers, and symbolic numbers too, picking the branch
    with angle between `-\\pi` and `\\pi`::

        sage: log(-1+0*I)
        I*pi
        sage: log(CC(-1))
        3.14159265358979*I
        sage: log(-1.0)
        3.14159265358979*I

    Small integer powers are factored out immediately::

        sage: log(4)
        2*log(2)
        sage: log(1000000000)
        9*log(10)
        sage: log(8) - 3*log(2)
        0
        sage: bool(log(8) == 3*log(2))
        True

    The ``hold`` parameter can be used to prevent automatic evaluation::

        sage: log(-1,hold=True)
        log(-1)
        sage: log(-1)
        I*pi
        sage: I.log(hold=True)
        log(I)
        sage: I.log(hold=True).simplify()
        1/2*I*pi

    For input zero, the following behavior occurs::

        sage: log(0)
        -Infinity
        sage: log(CC(0))
        -infinity
        sage: log(0.0)
        -infinity

    The log function also works in finite fields as long as the
    argument lies in the multiplicative group generated by the base::

        sage: F = GF(13); g = F.multiplicative_generator(); g
        2
        sage: a = F(8)
        sage: log(a,g); g^log(a,g)
        3
        8
        sage: log(a,3)
        Traceback (most recent call last):
        ...
        ValueError: no logarithm of 8 found to base 3 modulo 13
        sage: log(F(9), 3)
        2

    The log function also works for p-adics (see documentation for
    p-adics for more information)::

        sage: R = Zp(5); R
        5-adic Ring with capped relative precision 20
        sage: a = R(16); a
        1 + 3*5 + O(5^20)
        sage: log(a)
        3*5 + 3*5^2 + 3*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 2*5^8 + 5^9 +
        5^11 + 2*5^12 + 5^13 + 3*5^15 + 2*5^16 + 4*5^17 + 3*5^18 +
        3*5^19 + O(5^20)


    TESTS:

    Check if :trac:`10136` is fixed::

        sage: ln(x).operator() is ln
        True
        sage: log(x).operator() is ln
        True

        sage: log(1000, 10)
        3
        sage: log(3,-1)
        -I*log(3)/pi
        sage: log(int(8),2)
        3
        sage: log(8,int(2))
        3
        sage: log(8,2)
        3
        sage: log(1/8,2)
        -3
        sage: log(1/8,1/2)
        3
        sage: log(8,1/2)
        -3

        sage: log(1000, 10, base=5)
        Traceback (most recent call last):
        ...
        TypeError: log takes at most 2 arguments (3 given)

    Check if :trac:`29164` is fixed::

        sage: log(0, 2)
        -Infinity
    """
    base = kwds.pop('base', None)
    if base:
        args = args + (base,)
    if not args:
        raise TypeError("log takes at least 1 arguments (0 given)")
    if len(args) == 1:
        from sage.functions.log import ln
        return ln(args[0], **kwds)
    if len(args) > 2:
        raise TypeError("log takes at most 2 arguments (%s given)" % (len(args) + 1 - (base is not None)))
    try:
        return args[0].log(args[1])
    except ValueError as ex:
        if ex.args[0].startswith("no logarithm"):
            raise
    except (AttributeError, TypeError):
        pass
    from sage.functions.log import logb
    return logb(args[0], args[1])


def minimal_polynomial(x, var='x'):
    """
    Return the minimal polynomial of ``x``.

    EXAMPLES::

        sage: a = matrix(ZZ, 2, [1..4])
        sage: minpoly(a)
        x^2 - 5*x - 2
        sage: minpoly(a,'t')
        t^2 - 5*t - 2
        sage: minimal_polynomial(a)
        x^2 - 5*x - 2
        sage: minimal_polynomial(a,'theta')
        theta^2 - 5*theta - 2
    """
    try:
        return x.minpoly(var=var)
    except AttributeError:
        return x.minimal_polynomial(var=var)


minpoly = minimal_polynomial


def multiplicative_order(x):
    r"""
    Return the multiplicative order of ``x``, if ``x`` is a unit, or
    raise ``ArithmeticError`` otherwise.

    EXAMPLES::

        sage: a = mod(5,11)
        sage: multiplicative_order(a)
        5
        sage: multiplicative_order(mod(2,11))
        10
        sage: multiplicative_order(mod(2,12))
        Traceback (most recent call last):
        ...
        ArithmeticError: multiplicative order of 2 not defined since it is not a unit modulo 12
    """
    return x.multiplicative_order()


def ngens(x):
    """
    Return the number of generators of ``x``.

    EXAMPLES::

        sage: R.<x,y> = SR[]; R
        Multivariate Polynomial Ring in x, y over Symbolic Ring
        sage: ngens(R)
        2
        sage: A = AbelianGroup(5, [5,5,7,8,9])
        sage: ngens(A)
        5
        sage: ngens(ZZ)
        1
    """
    return x.ngens()


def norm(x):
    r"""
    Return the norm of ``x``.

    For matrices and vectors, this returns the L2-norm. The L2-norm of a
    vector `\textbf{v} = (v_1, v_2, \dots, v_n)`, also called the Euclidean
    norm, is defined as

    .. MATH::

        |\textbf{v}|
        =
        \sqrt{\sum_{i=1}^n |v_i|^2}

    where `|v_i|` is the complex modulus of `v_i`. The Euclidean norm is often
    used for determining the distance between two points in two- or
    three-dimensional space.

    For complex numbers, the function returns the field norm. If
    `c = a + bi` is a complex number, then the norm of `c` is defined as the
    product of `c` and its complex conjugate:

    .. MATH::

        \text{norm}(c)
        =
        \text{norm}(a + bi)
        =
        c \cdot \overline{c}
        =
        a^2 + b^2.

    The norm of a complex number is different from its absolute value.
    The absolute value of a complex number is defined to be the square
    root of its norm. A typical use of the complex norm is in the
    integral domain `\ZZ[i]` of Gaussian integers, where the norm of
    each Gaussian integer `c = a + bi` is defined as its complex norm.

    For vector fields on a pseudo-Riemannian manifold `(M,g)`, the function
    returns the norm with respect to the metric `g`:

    .. MATH::

        |v| = \sqrt{g(v,v)}

    .. SEEALSO::

        - :meth:`sage.matrix.matrix2.Matrix.norm`

        - :meth:`sage.modules.free_module_element.FreeModuleElement.norm`

        - :meth:`sage.rings.complex_double.ComplexDoubleElement.norm`

        - :meth:`sage.rings.complex_mpfr.ComplexNumber.norm`

        - :meth:`sage.symbolic.expression.Expression.norm`

        - :meth:`sage.manifolds.differentiable.vectorfield.VectorField.norm`

    EXAMPLES:

    The norm of vectors::

        sage: z = 1 + 2*I
        sage: norm(vector([z]))
        sqrt(5)
        sage: v = vector([-1,2,3])
        sage: norm(v)
        sqrt(14)
        sage: _ = var("a b c d", domain='real')
        sage: v = vector([a, b, c, d])
        sage: norm(v)
        sqrt(a^2 + b^2 + c^2 + d^2)

    The norm of matrices::

        sage: z = 1 + 2*I
        sage: norm(matrix([[z]]))
        2.23606797749979
        sage: M = matrix(ZZ, [[1,2,4,3], [-1,0,3,-10]])
        sage: norm(M)  # abs tol 1e-14
        10.690331129154467
        sage: norm(CDF(z))
        5.0
        sage: norm(CC(z))
        5.00000000000000

    The norm of complex numbers::

        sage: z = 2 - 3*I
        sage: norm(z)
        13
        sage: a = randint(-10^10, 100^10)
        sage: b = randint(-10^10, 100^10)
        sage: z = a + b*I
        sage: bool(norm(z) == a^2 + b^2)
        True

    The complex norm of symbolic expressions::

        sage: a, b, c = var("a, b, c")
        sage: assume((a, 'real'), (b, 'real'), (c, 'real'))
        sage: z = a + b*I
        sage: bool(norm(z).simplify() == a^2 + b^2)
        True
        sage: norm(a + b).simplify()
        a^2 + 2*a*b + b^2
        sage: v = vector([a, b, c])
        sage: bool(norm(v).simplify() == sqrt(a^2 + b^2 + c^2))
        True
        sage: forget()
    """
    return x.norm()


def numerator(x):
    """
    Return the numerator of ``x``.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: F = FractionField(R)
        sage: r = (x+1)/(x-1)
        sage: numerator(r)
        x + 1
        sage: numerator(17/11111)
        17
    """
    if isinstance(x, int):
        return x
    return x.numerator()


def numerical_approx(x, prec=None, digits=None, algorithm=None):
    r"""
    Return a numerical approximation of ``self`` with ``prec`` bits
    (or decimal ``digits``) of precision.

    No guarantee is made about the accuracy of the result.

    .. NOTE::

        Lower case :func:`n` is an alias for :func:`numerical_approx`
        and may be used as a method.

    INPUT:

    - ``prec`` -- precision in bits

    - ``digits`` -- precision in decimal digits (only used if
      ``prec`` is not given)

    - ``algorithm`` -- which algorithm to use to compute this
      approximation (the accepted algorithms depend on the object)

    If neither ``prec`` nor ``digits`` is given, the default
    precision is 53 bits (roughly 16 digits).

    EXAMPLES::

        sage: numerical_approx(pi, 10)
        3.1
        sage: numerical_approx(pi, digits=10)
        3.141592654
        sage: numerical_approx(pi^2 + e, digits=20)
        12.587886229548403854
        sage: n(pi^2 + e)
        12.5878862295484
        sage: N(pi^2 + e)
        12.5878862295484
        sage: n(pi^2 + e, digits=50)
        12.587886229548403854194778471228813633070946500941
        sage: a = CC(-5).n(prec=40)
        sage: b = ComplexField(40)(-5)
        sage: a == b
        True
        sage: parent(a) is parent(b)
        True
        sage: numerical_approx(9)
        9.00000000000000

    You can also usually use method notation::

        sage: (pi^2 + e).n()
        12.5878862295484
        sage: (pi^2 + e).numerical_approx()
        12.5878862295484

    Vectors and matrices may also have their entries approximated::

        sage: v = vector(RDF, [1,2,3])
        sage: v.n()
        (1.00000000000000, 2.00000000000000, 3.00000000000000)

        sage: v = vector(CDF, [1,2,3])
        sage: v.n()
        (1.00000000000000, 2.00000000000000, 3.00000000000000)
        sage: _.parent()
        Vector space of dimension 3 over Complex Field with 53 bits of precision
        sage: v.n(prec=20)
        (1.0000, 2.0000, 3.0000)

        sage: u = vector(QQ, [1/2, 1/3, 1/4])
        sage: n(u, prec=15)
        (0.5000, 0.3333, 0.2500)
        sage: n(u, digits=5)
        (0.50000, 0.33333, 0.25000)

        sage: v = vector(QQ, [1/2, 0, 0, 1/3, 0, 0, 0, 1/4], sparse=True)
        sage: u = v.numerical_approx(digits=4)
        sage: u.is_sparse()
        True
        sage: u
        (0.5000, 0.0000, 0.0000, 0.3333, 0.0000, 0.0000, 0.0000, 0.2500)

        sage: A = matrix(QQ, 2, 3, range(6))
        sage: A.n()
        [0.000000000000000  1.00000000000000  2.00000000000000]
        [ 3.00000000000000  4.00000000000000  5.00000000000000]

        sage: B = matrix(Integers(12), 3, 8, srange(24))
        sage: N(B, digits=2)
        [0.00  1.0  2.0  3.0  4.0  5.0  6.0  7.0]
        [ 8.0  9.0  10.  11. 0.00  1.0  2.0  3.0]
        [ 4.0  5.0  6.0  7.0  8.0  9.0  10.  11.]

    Internally, numerical approximations of real numbers are stored in base-2.
    Therefore, numbers which look the same in their decimal expansion might be
    different::

        sage: x=N(pi, digits=3); x
        3.14
        sage: y=N(3.14, digits=3); y
        3.14
        sage: x==y
        False
        sage: x.str(base=2)
        '11.001001000100'
        sage: y.str(base=2)
        '11.001000111101'

    Increasing the precision of a floating point number is not allowed::

        sage: CC(-5).n(prec=100)
        Traceback (most recent call last):
        ...
        TypeError: cannot approximate to a precision of 100 bits, use at most 53 bits
        sage: n(1.3r, digits=20)
        Traceback (most recent call last):
        ...
        TypeError: cannot approximate to a precision of 70 bits, use at most 53 bits
        sage: RealField(24).pi().n()
        Traceback (most recent call last):
        ...
        TypeError: cannot approximate to a precision of 53 bits, use at most 24 bits

    As an exceptional case, ``digits=1`` usually leads to 2 digits (one
    significant) in the decimal output (see :trac:`11647`)::

        sage: N(pi, digits=1)
        3.2
        sage: N(pi, digits=2)
        3.1
        sage: N(100*pi, digits=1)
        320.
        sage: N(100*pi, digits=2)
        310.

    In the following example, ``pi`` and ``3`` are both approximated to two
    bits of precision and then subtracted, which kills two bits of precision::

        sage: N(pi, prec=2)
        3.0
        sage: N(3, prec=2)
        3.0
        sage: N(pi - 3, prec=2)
        0.00

    TESTS::

        sage: numerical_approx(I)
        1.00000000000000*I
        sage: x = QQ['x'].gen()
        sage: F.<k> = NumberField(x^2+2, embedding=sqrt(CC(2))*CC.0)
        sage: numerical_approx(k)
        1.41421356237309*I

        sage: type(numerical_approx(CC(1/2)))
        <class 'sage.rings.complex_mpfr.ComplexNumber'>

    The following tests :trac:`10761`, in which ``n()`` would break when
    called on complex-valued algebraic numbers.  ::

        sage: E = matrix(3, [3,1,6,5,2,9,7,3,13]).eigenvalues(); E
        [18.16815365088822?, -0.08407682544410650? - 0.2190261484802906?*I, -0.08407682544410650? + 0.2190261484802906?*I]
        sage: E[1].parent()
        Algebraic Field
        sage: [a.n() for a in E]
        [18.1681536508882, -0.0840768254441065 - 0.219026148480291*I, -0.0840768254441065 + 0.219026148480291*I]

    Make sure we've rounded up log(10,2) enough to guarantee
    sufficient precision (:trac:`10164`)::

        sage: ks = 4*10**5, 10**6
        sage: check_str_length = lambda k: len(str(numerical_approx(1+10**-k,digits=k+1)))-1 >= k+1
        sage: check_precision = lambda k: numerical_approx(1+10**-k,digits=k+1)-1 > 0
        sage: all(check_str_length(k) and check_precision(k) for k in ks)
        True

    Testing we have sufficient precision for the golden ratio (:trac:`12163`), note
    that the decimal point adds 1 to the string length::

        sage: len(str(n(golden_ratio, digits=5000)))
        5001
        sage: len(str(n(golden_ratio, digits=5000000)))  # long time (4s on sage.math, 2012)
        5000001

    Check that :trac:`14778` is fixed::

        sage: n(0, algorithm='foo')
        0.000000000000000
    """
    if prec is None:
        from sage.arith.numerical_approx import digits_to_bits
        prec = digits_to_bits(digits)
    try:
        n = x.numerical_approx
    except AttributeError:
        from sage.arith.numerical_approx import numerical_approx_generic
        return numerical_approx_generic(x, prec)
    else:
        return n(prec, algorithm=algorithm)


n = numerical_approx
N = numerical_approx


def objgens(x):
    """
    EXAMPLES::

        sage: R, x = objgens(PolynomialRing(QQ,3, 'x'))
        sage: R
        Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        sage: x
        (x0, x1, x2)
    """
    return x.objgens()


def objgen(x):
    """
    EXAMPLES::

        sage: R, x = objgen(FractionField(QQ['x']))
        sage: R
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: x
        x
    """
    return x.objgen()


def order(x):
    """
    Return the order of ``x``.

    If ``x`` is a ring or module element, this is
    the additive order of ``x``.

    EXAMPLES::

        sage: C = CyclicPermutationGroup(10)
        sage: order(C)
        10
        sage: F = GF(7)
        sage: order(F)
        7
    """
    return x.order()


def rank(x):
    """
    Return the rank of ``x``.

    EXAMPLES:

    We compute the rank of a matrix::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: rank(A)
        2

    We compute the rank of an elliptic curve::

        sage: E = EllipticCurve([0,0,1,-1,0])
        sage: rank(E)
        1
    """
    return x.rank()


def regulator(x):
    """
    Return the regulator of ``x``.

    EXAMPLES::

        sage: regulator(NumberField(x^2-2, 'a'))
        0.881373587019543
        sage: regulator(EllipticCurve('11a'))
        1.00000000000000
    """
    return x.regulator()


def round(x, ndigits=0):
    """
    round(number[, ndigits]) - double-precision real number

    Round a number to a given precision in decimal digits (default 0
    digits). If no precision is specified this just calls the element's
    .round() method.

    EXAMPLES::

        sage: round(sqrt(2),2)
        1.41
        sage: q = round(sqrt(2),5); q
        1.41421
        sage: type(q)
        <class 'sage.rings.real_double...RealDoubleElement...'>
        sage: q = round(sqrt(2)); q
        1
        sage: type(q)
        <class 'sage.rings.integer.Integer'>
        sage: round(pi)
        3
        sage: b = 5.4999999999999999
        sage: round(b)
        5

    This example addresses :trac:`23502`::

        sage: n = round(6); type(n)
        <class 'sage.rings.integer.Integer'>

    Since we use floating-point with a limited range, some roundings can't
    be performed::

        sage: round(sqrt(Integer('1'*1000)),2)
        +infinity

    IMPLEMENTATION: If ndigits is specified, it calls Python's builtin
    round function, and converts the result to a real double field
    element. Otherwise, it tries the argument's .round() method; if
    that fails, it reverts to the builtin round function, converted to
    a real double field element.

    .. NOTE::

        This is currently slower than the builtin round function, since it does
        more work - i.e., allocating an RDF element and initializing it. To
        access the builtin version do ``import builtins; builtins.round``.
    """
    try:
        if ndigits:
            x = float(x)
            return RealDoubleElement(builtins.round(x, ndigits))
        else:
            try:
                return x.round()
            except AttributeError:
                return RealDoubleElement(builtins.round(x, 0))
    except ArithmeticError:
        if not isinstance(x, RealDoubleElement):
            return round(RDF(x), ndigits)
        else:
            raise


def quotient(x, y, *args, **kwds):
    """
    Return the quotient object x/y, e.g., a quotient of numbers or of a
    polynomial ring x by the ideal generated by y, etc.

    EXAMPLES::

        sage: quotient(5,6)
        5/6
        sage: quotient(5.,6.)
        0.833333333333333
        sage: R.<x> = ZZ[]; R
        Univariate Polynomial Ring in x over Integer Ring
        sage: I = Ideal(R, x^2+1)
        sage: quotient(R, I)
        Univariate Quotient Polynomial Ring in xbar over Integer Ring with modulus x^2 + 1
    """
    try:
        return x.quotient(y, *args, **kwds)
    except AttributeError:
        return x / y


quo = quotient


def isqrt(x):
    """
    Return an integer square root, i.e., the floor of a square root.

    EXAMPLES::

        sage: isqrt(10)
        3
        sage: isqrt(10r)
        3
    """
    try:
        return x.isqrt()
    except AttributeError:
        from sage.functions.all import floor
        n = Integer(floor(x))
        return n.isqrt()


def squarefree_part(x):
    """
    Return the square free part of ``x``, i.e., a divisor
    `z` such that `x = z y^2`, for a perfect square
    `y^2`.

    EXAMPLES::

        sage: squarefree_part(100)
        1
        sage: squarefree_part(12)
        3
        sage: squarefree_part(10)
        10
        sage: squarefree_part(216r) # see #8976
        6

    ::

        sage: x = QQ['x'].0
        sage: S = squarefree_part(-9*x*(x-6)^7*(x-3)^2); S
        -9*x^2 + 54*x
        sage: S.factor()
        (-9) * (x - 6) * x

    ::

        sage: f = (x^3 + x + 1)^3*(x-1); f
        x^10 - x^9 + 3*x^8 + 3*x^5 - 2*x^4 - x^3 - 2*x - 1
        sage: g = squarefree_part(f); g
        x^4 - x^3 + x^2 - 1
        sage: g.factor()
        (x - 1) * (x^3 + x + 1)
    """
    try:
        return x.squarefree_part()
    except AttributeError:
        pass
    from sage.arith.all import factor
    from sage.structure.all import parent
    F = factor(x)
    n = parent(x)(1)
    for p, e in F:
        if e % 2:
            n *= p
    return n * F.unit()


def _do_sqrt(x, prec=None, extend=True, all=False):
    r"""
    Used internally to compute the square root of x.

    INPUT:

    -  ``x`` - a number

    -  ``prec`` - None (default) or a positive integer
       (bits of precision) If not None, then compute the square root
       numerically to prec bits of precision.

    -  ``extend`` - bool (default: True); this is a place
       holder, and is always ignored since in the symbolic ring everything
       has a square root.

    -  ``extend`` - bool (default: True); whether to extend
       the base ring to find roots. The extend parameter is ignored if
       prec is a positive integer.

    -  ``all`` - bool (default: False); whether to return
       a list of all the square roots of x.


    EXAMPLES::

        sage: from sage.misc.functional import _do_sqrt
        sage: _do_sqrt(3)
        sqrt(3)
        sage: _do_sqrt(3,prec=10)
        1.7
        sage: _do_sqrt(3,prec=100)
        1.7320508075688772935274463415
        sage: _do_sqrt(3,all=True)
        [sqrt(3), -sqrt(3)]

    Note that the extend parameter is ignored in the symbolic ring::

        sage: _do_sqrt(3,extend=False)
        sqrt(3)
    """
    if prec:
        if x >= 0:
            from sage.rings.real_mpfr import RealField
            return RealField(prec)(x).sqrt(all=all)
        else:
            from sage.rings.complex_mpfr import ComplexField
            return ComplexField(prec)(x).sqrt(all=all)
    if x == -1:
        from sage.symbolic.expression import I
        z = I
    else:
        from sage.symbolic.ring import SR
        z = SR(x).sqrt()

    if all:
        if z:
            return [z, -z]
        else:
            return [z]
    return z


def sqrt(x, *args, **kwds):
    r"""
    INPUT:

    -  ``x`` - a number

    -  ``prec`` - integer (default: None): if None, returns
       an exact square root; otherwise returns a numerical square root if
       necessary, to the given bits of precision.

    -  ``extend`` - bool (default: True); this is a place
       holder, and is always ignored or passed to the sqrt function for x,
       since in the symbolic ring everything has a square root.

    -  ``all`` - bool (default: False); if True, return all
       square roots of self, instead of just one.

    EXAMPLES::

        sage: sqrt(-1)
        I
        sage: sqrt(2)
        sqrt(2)
        sage: sqrt(2)^2
        2
        sage: sqrt(4)
        2
        sage: sqrt(4,all=True)
        [2, -2]
        sage: sqrt(x^2)
        sqrt(x^2)

    For a non-symbolic square root, there are a few options.
    The best is to numerically approximate afterward::

        sage: sqrt(2).n()
        1.41421356237310
        sage: sqrt(2).n(prec=100)
        1.4142135623730950488016887242

    Or one can input a numerical type.

        sage: sqrt(2.)
        1.41421356237310
        sage: sqrt(2.000000000000000000000000)
        1.41421356237309504880169
        sage: sqrt(4.0)
        2.00000000000000

    To prevent automatic evaluation, one can use the ``hold`` parameter
    after coercing to the symbolic ring::

        sage: sqrt(SR(4),hold=True)
        sqrt(4)
        sage: sqrt(4,hold=True)
        Traceback (most recent call last):
        ...
        TypeError: ..._do_sqrt() got an unexpected keyword argument 'hold'

    This illustrates that the bug reported in :trac:`6171` has been fixed::

        sage: a = 1.1
        sage: a.sqrt(prec=100)  # this is supposed to fail
        Traceback (most recent call last):
        ...
        TypeError: ...sqrt() got an unexpected keyword argument 'prec'
        sage: sqrt(a, prec=100)
        1.0488088481701515469914535137
        sage: sqrt(4.00, prec=250)
        2.0000000000000000000000000000000000000000000000000000000000000000000000000

    One can use numpy input as well::

        sage: import numpy
        sage: a = numpy.arange(2,5)
        sage: sqrt(a)
        array([1.41421356, 1.73205081, 2.        ])
    """
    if isinstance(x, float):
        return math.sqrt(x)
    elif type(x).__module__ == 'numpy':
        from numpy import sqrt
        return sqrt(x)
    try:
        return x.sqrt(*args, **kwds)
    # The following includes TypeError to catch cases where sqrt
    # is called with a "prec" keyword, for example, but the sqrt
    # method for x doesn't accept such a keyword.
    except (AttributeError, TypeError):
        pass
    return _do_sqrt(x, *args, **kwds)


def transpose(x):
    """
    Return the transpose of ``x``.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: transpose(A)
        [1 4 7]
        [2 5 8]
        [3 6 9]
    """
    return x.transpose()
