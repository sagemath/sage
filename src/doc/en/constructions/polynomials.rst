***********
Polynomials
***********

.. index::
   pair: polynomial; powers

.. _section-polynomialpower:

Polynomial powers
=================

How do I compute modular polynomial powers in Sage?

To compute :math:`x^{2006} \pmod {x^3 + 7}` in
:math:`GF(97)[x]`, we create the quotient ring
:math:`GF(97)[x]/(x^3+7)`, and compute :math:`x^{2006}` in it.
As a matter of Sage notation, we must distinguish between the
:math:`x` in :math:`GF(97)[x]` and the corresponding element
(which we denote by :math:`a`) in the quotient ring
:math:`GF(97)[x]/(x^3+7)`.

::

    sage: R = PolynomialRing(GF(97),'x')
    sage: x = R.gen()
    sage: S = R.quotient(x^3 + 7, 'a')
    sage: a = S.gen()
    sage: S
    Univariate Quotient Polynomial Ring in a over
    Finite Field of size 97 with modulus x^3 + 7
    sage: a^2006
    4*a^2

Another approach to this:

::

    sage: R = PolynomialRing(GF(97),'x')
    sage: x = R.gen()
    sage: S = R.quotient(x^3 + 7, 'a')
    sage: a = S.gen()
    sage: a^20062006
    80*a
    sage: libgap.eval("R:= PolynomialRing( GF(97))")
    GF(97)[x_1]
    sage: libgap.eval("i:= IndeterminatesOfPolynomialRing(R)")
    [ x_1 ]
    sage: libgap.eval("x:= i[1]"); libgap.eval("f:= x;")
    x_1
    x_1
    sage: libgap.eval("PowerMod( R, x, 20062006, x^3+7 );")
    Z(97)^41*x_1
    sage: libgap.eval("PowerMod( R, x, 20062006, x^3+7 );")
    Z(97)^41*x_1
    sage: libgap.eval("PowerMod( R, x, 2006200620062006, x^3+7 );")
    Z(97)^4*x_1^2
    sage: a^2006200620062006
    43*a^2
    sage: libgap.eval("PowerMod( R, x, 2006200620062006, x^3+7 );")
    Z(97)^4*x_1^2
    sage: libgap.eval("Int(Z(97)^4)")
    43

.. index::
   pair: polynomial; factorization

.. _section-factor:

Factorization
=============

You can factor a polynomial using Sage.

Using Sage to factor a univariate polynomial is a matter of
applying the method ``factor`` to the PolynomialRingElement object f.
In fact, this method actually calls Pari, so the computation is
fairly fast.

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f = (x^3 - 1)^2-(x^2-1)^2
    sage: f.factor()
    (x - 1)^2 * x^2 * (x^2 + 2*x + 2)

Using the Singular interface, Sage also factors multivariate
polynomials.

::

    sage: x, y = PolynomialRing(RationalField(), 2, ['x','y']).gens()
    sage: f =  (9*y^6 - 9*x^2*y^5 - 18*x^3*y^4 - 9*x^5*y^4 + 9*x^6*y^2 + 9*x^7*y^3
    ....:     + 18*x^8*y^2 - 9*x^11)
    sage: f.factor()
    (9) * (-x^5 + y^2) * (x^6 - 2*x^3*y^2 - x^2*y^3 + y^4)

.. index::
   pair: polynomial; gcd

Polynomial GCD's
================

This example illustrates single variable polynomial GCD's:

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f = 3*x^3 + x
    sage: g = 9*x*(x+1)
    sage: f.gcd(g)
    x

This example illustrates multivariate polynomial GCD's:

::

    sage: R.<x,y,z> = PolynomialRing(RationalField(), order='lex')
    sage: f = 3*x^2*(x+y)
    sage: g = 9*x*(y^2 - x^2)
    sage: f.gcd(g)
    x^2 + x*y

Here's another way to do this:

::

    sage: R2 = singular.ring(0, '(x,y,z)', 'lp')
    sage: a = singular.new('3x2*(x+y)')
    sage: b = singular.new('9x*(y2-x2)')
    sage: g = a.gcd(b)
    sage: g
    x^2+x*y

This example illustrates univariate polynomial GCD's via the GAP
interface.

::

    sage: R = libgap.PolynomialRing(GF(2)); R
    GF(2)[x_1]
    sage: i = R.IndeterminatesOfPolynomialRing(); i
    [ x_1 ]
    sage: x_1 = i[0]
    sage: f = (x_1^3 - x_1 + 1)*(x_1 + x_1^2); f
    x_1^5+x_1^4+x_1^3+x_1
    sage: g = (x_1^3 - x_1 + 1)*(x_1 + 1); g
    x_1^4+x_1^3+x_1^2+Z(2)^0
    sage: f.Gcd(g)
    x_1^4+x_1^3+x_1^2+Z(2)^0

We can, of course, do the same computation in , which uses the NTL
library (which does huge polynomial gcd's over finite fields very
quickly).

::

    sage: x = PolynomialRing(GF(2), 'x').gen()
    sage: f = (x^3 - x + 1)*(x + x^2); f
    x^5 + x^4 + x^3 + x
    sage: g = (x^3 - x + 1)*(x + 1)
    sage: f.gcd(g)
    x^4 + x^3 + x^2 + 1

.. index::
   pair: polynomial; roots

.. _section-roots:

Roots of polynomials
====================

Sage can compute roots of a univariant polynomial.

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f = x^3 - 1
    sage: f.roots()
    [(1, 1)]
    sage: f = (x^3 - 1)^2
    sage: f.roots()
    [(1, 2)]
    sage: x = PolynomialRing(CyclotomicField(3), 'x').gen()
    sage: f = x^3 - 1
    sage: f.roots()
    [(1, 1), (zeta3, 1), (-zeta3 - 1, 1)]

The first of the pair is the root, the second of the pair is its
multiplicity.

There are some situations where GAP does find the roots of a
univariate polynomial but GAP does not do this generally. (The
roots must generate either a finite field or a subfield of a
cyclotomic field.) However, there is a GAP package called
``RadiRoot``, which must be installed into 's installation of GAP,
which does help to do this for polynomials with rational
coefficients (``radiroot`` itself requires other packages to be
installed; please see its webpage for more details). The ``Factors``
command actually has an option which allows you to increase the
groundfield so that a factorization actually returns the roots.
Please see the examples given in section 64.10
"Polynomial Factorization" of the GAP Reference Manual for more
details.

.. index::
   pair: polynomial; evaluation

.. _section-evaluate:

Evaluation of multivariate functions
====================================

You can evaluate polynomials in Sage as usual by
substituting in points:

::

    sage: x = PolynomialRing(RationalField(), 3, 'x').gens()
    sage: f = x[0] + x[1] - 2*x[1]*x[2]
    sage: f
    -2*x1*x2 + x0 + x1
    sage: f(1,2,0)
    3
    sage: f(1,2,5)
    -17

This also will work with rational functions:

.. link

::

    sage: h = f /(x[1] + x[2])
    sage: h
    (-2*x1*x2 + x0 + x1)/(x1 + x2)
    sage: h(1,2,3)
    -9/5

.. index::
   pair: polynomial; symbolic manipulation

Sage also performs symbolic manipulation:

::

    sage: var('x,y,z')
    (x, y, z)
    sage: f = (x + 3*y + x^2*y)^3; f
    (x^2*y + x + 3*y)^3
    sage: f(x=1,y=2,z=3)
    729
    sage: f.expand()
    x^6*y^3 + 3*x^5*y^2 + 9*x^4*y^3 + 3*x^4*y + 18*x^3*y^2 +
    27*x^2*y^3 +
    x^3 + 9*x^2*y + 27*x*y^2 + 27*y^3
    sage: f(x = 5/z)
    (3*y + 25*y/z^2 + 5/z)^3
    sage: g = f.subs(x = 5/z); g
    (3*y + 25*y/z^2 + 5/z)^3
    sage: h = g.rational_simplify(); h
    (27*y^3*z^6 + 135*y^2*z^5 + 225*(3*y^3 + y)*z^4 + 125*(18*y^2 + 1)*z^3 +
    15625*y^3 + 9375*y^2*z + 1875*(3*y^3 + y)*z^2)/z^6

Roots of multivariate polynomials
=================================

Sage (using the interface to Singular) can solve multivariate polynomial
equations in some situations (they assume that the solutions form a
zero-dimensional variety) using Gröbner bases. Here is a simple
example:

::

    sage: R = PolynomialRing(QQ, 2, 'ab', order='lp')
    sage: a,b = R.gens()
    sage: I = (a^2-b^2-3, a-2*b)*R
    sage: B = I.groebner_basis(); B
    [a - 2*b, b^2 - 1]

So :math:`b=\pm 1` and :math:`a=2b`.

.. index:
   pair: polynomial; Groebner basis of ideal

.. _section-groebner:

Gröbner bases
=============

This computation uses Singular behind the scenes to
compute the Gröbner basis.

::

    sage: R = PolynomialRing(QQ, 4, 'abcd', order='lp')
    sage: a,b,c,d = R.gens()
    sage: I = (a+b+c+d, a*b+a*d+b*c+c*d, a*b*c+a*b*d+a*c*d+b*c*d, a*b*c*d-1)*R; I
    Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d + a*c*d + b*c*d,
    a*b*c*d - 1) of Multivariate Polynomial Ring in a, b, c, d over Rational Field
    sage: B = I.groebner_basis(); B
    [a + b + c + d,
     b^2 + 2*b*d + d^2,
     b*c - b*d + c^2*d^4 + c*d - 2*d^2,
     b*d^4 - b + d^5 - d,
     c^3*d^2 + c^2*d^3 - c - d,
     c^2*d^6 - c^2*d^2 - d^4 + 1]

You can work with multiple rings without having to switch back and
forth like in Singular. For example,

::

    sage: a,b,c = QQ['a,b,c'].gens()
    sage: X,Y = GF(7)['X,Y'].gens()
    sage: I = ideal(a, b^2, b^3+c^3)
    sage: J = ideal(X^10 + Y^10)

    sage: I.minimal_associated_primes ()
    [Ideal (c, b, a) of Multivariate Polynomial Ring in a, b, c over Rational Field]

    sage: J.minimal_associated_primes ()     # slightly random output
    [Ideal (Y^4 + 3*X*Y^3 + 4*X^2*Y^2 + 4*X^3*Y + X^4) of Multivariate Polynomial
    Ring in X, Y over Finite Field of size 7,
     Ideal (Y^4 + 4*X*Y^3 + 4*X^2*Y^2 + 3*X^3*Y + X^4) of Multivariate Polynomial
    Ring in X, Y over Finite Field of size 7,
     Ideal (Y^2 + X^2) of Multivariate Polynomial Ring in X, Y over Finite Field
    of size 7]

All the real work is done by Singular.

Sage also includes ``gfan`` which provides other fast algorithms for
computing Gröbner bases. See the section on "Gröbner fans" in the
Reference Manual for more details.
