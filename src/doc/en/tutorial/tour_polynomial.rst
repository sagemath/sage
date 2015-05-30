.. _section-poly:

Polynomials
===========

In this section we illustrate how to create and use
polynomials in Sage.


.. _section-univariate:

Univariate Polynomials
----------------------

There are three ways to create polynomial rings.

::

    sage: R = PolynomialRing(QQ, 't')
    sage: R
    Univariate Polynomial Ring in t over Rational Field

This creates a polynomial ring and tells Sage to use (the string)
't' as the indeterminate when printing to the screen. However, this
does not define the symbol ``t`` for use in Sage, so you cannot use it to
enter a polynomial (such as :math:`t^2+1`) belonging to ``R``.

An alternate way is

.. link

::

    sage: S = QQ['t']
    sage: S == R
    True

This has the same issue regarding ``t``.

A third very convenient way is

::

    sage: R.<t> = PolynomialRing(QQ)

or

::

    sage: R.<t> = QQ['t']

or even

::

    sage: R.<t> = QQ[]

This has the additional side effect that it defines the variable
``t`` to be the indeterminate of the polynomial ring, so you can
easily construct elements of ``R``, as follows. (Note that the third
way is very similar to the constructor notation in Magma, and just
as in Magma it can be used for a wide range of objects.)

.. link

::

    sage: poly = (t+1) * (t+2); poly
    t^2 + 3*t + 2
    sage: poly in R
    True

Whatever method you use to define a polynomial ring, you can
recover the indeterminate as the :math:`0^{th}` generator:

::

    sage: R = PolynomialRing(QQ, 't')
    sage: t = R.0
    sage: t in R
    True

Note that a similar construction works with the complex numbers:
the complex numbers can be viewed as being generated over the real
numbers by the symbol ``i``; thus we have the following:

::

    sage: CC
    Complex Field with 53 bits of precision
    sage: CC.0  # 0th generator of CC
    1.00000000000000*I

For polynomial rings, you can obtain both the ring and its
generator, or just the generator, during ring creation as follows:

::

    sage: R, t = QQ['t'].objgen()
    sage: t    = QQ['t'].gen()
    sage: R, t = objgen(QQ['t'])
    sage: t    = gen(QQ['t'])

Finally we do some arithmetic in :math:`\QQ[t]`.

::

    sage: R, t = QQ['t'].objgen()
    sage: f = 2*t^7 + 3*t^2 - 15/19
    sage: f^2
    4*t^14 + 12*t^9 - 60/19*t^7 + 9*t^4 - 90/19*t^2 + 225/361
    sage: cyclo = R.cyclotomic_polynomial(7); cyclo
    t^6 + t^5 + t^4 + t^3 + t^2 + t + 1
    sage: g = 7 * cyclo * t^5 * (t^5 + 10*t + 2)
    sage: g
    7*t^16 + 7*t^15 + 7*t^14 + 7*t^13 + 77*t^12 + 91*t^11 + 91*t^10 + 84*t^9
           + 84*t^8 + 84*t^7 + 84*t^6 + 14*t^5
    sage: F = factor(g); F
    (7) * t^5 * (t^5 + 10*t + 2) * (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)
    sage: F.unit()
    7
    sage: list(F)
    [(t, 5), (t^5 + 10*t + 2, 1), (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1, 1)]

Notice that the factorization correctly takes into account and
records the unit part.

If you were to use, e.g., the ``R.cyclotomic_polynomial`` function a
lot for some research project, in addition to citing Sage you should
make an attempt to find out what component of Sage is being used to
actually compute the cyclotomic polynomial and cite that as well.
In this case, if you type ``R.cyclotomic_polynomial??`` to see the
source code, you'll quickly see a line ``f = pari.polcyclo(n)`` which
means that PARI is being used for computation of the cyclotomic
polynomial. Cite PARI in your work as well.

Dividing two polynomials constructs an element of the fraction
field (which Sage creates automatically).

::

    sage: x = QQ['x'].0
    sage: f = x^3 + 1; g = x^2 - 17
    sage: h = f/g;  h
    (x^3 + 1)/(x^2 - 17)
    sage: h.parent()
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

Using Laurent series, one can compute series expansions in the
fraction field of ``QQ[x]``:

::

    sage: R.<x> = LaurentSeriesRing(QQ); R
    Laurent Series Ring in x over Rational Field
    sage: 1/(1-x) + O(x^10)
    1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)

If we name the variable differently, we obtain a different
univariate polynomial ring.

::

    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<y> = PolynomialRing(QQ)
    sage: x == y
    False
    sage: R == S
    False
    sage: R(y)
    x
    sage: R(y^2 - 17)
    x^2 - 17

The ring is determined by the variable. Note that making another
ring with variable called ``x`` does not return a different ring.

::

    sage: R = PolynomialRing(QQ, "x")
    sage: T = PolynomialRing(QQ, "x")
    sage: R == T
    True
    sage: R is T
    True
    sage: R.0 == T.0
    True

Sage also has support for power series and Laurent series rings
over any base ring. In the following example, we create an element
of :math:`\GF{7}[[T]]` and divide to create an element of
:math:`\GF{7}((T))`.

::

    sage: R.<T> = PowerSeriesRing(GF(7)); R
    Power Series Ring in T over Finite Field of size 7
    sage: f = T  + 3*T^2 + T^3 + O(T^4)
    sage: f^3
    T^3 + 2*T^4 + 2*T^5 + O(T^6)
    sage: 1/f
    T^-1 + 4 + T + O(T^2)
    sage: parent(1/f)
    Laurent Series Ring in T over Finite Field of size 7

You can also create power series rings using a double-brackets
shorthand:

::

    sage: GF(7)[['T']]
    Power Series Ring in T over Finite Field of size 7

Multivariate Polynomials
------------------------

To work with polynomials of several variables, we declare the
polynomial ring and variables first.

::

    sage: R = PolynomialRing(GF(5),3,"z") # here, 3 = number of variables
    sage: R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Just as for defining univariate polynomial rings, there are
alternative ways:

::

    sage: GF(5)['z0, z1, z2']
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5
    sage: R.<z0,z1,z2> = GF(5)[]; R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Also, if you want the variable names to be single letters then you
can use the following shorthand:

::

    sage: PolynomialRing(GF(5), 3, 'xyz')
    Multivariate Polynomial Ring in x, y, z over Finite Field of size 5

Next let's do some arithmetic.

::

    sage: z = GF(5)['z0, z1, z2'].gens()
    sage: z
    (z0, z1, z2)
    sage: (z[0]+z[1]+z[2])^2
    z0^2 + 2*z0*z1 + z1^2 + 2*z0*z2 + 2*z1*z2 + z2^2

You can also use more mathematical notation to construct a
polynomial ring.

::

    sage: R = GF(5)['x,y,z']
    sage: x,y,z = R.gens()
    sage: QQ['x']
    Univariate Polynomial Ring in x over Rational Field
    sage: QQ['x,y'].gens()
    (x, y)
    sage: QQ['x'].objgens()
    (Univariate Polynomial Ring in x over Rational Field, (x,))

Multivariate polynomials are implemented in Sage using Python
dictionaries and the "distributive representation" of a polynomial.
Sage makes some use of Singular [Si]_, e.g., for computation of
gcd's and Gröbner basis of ideals.

::

    sage: R, (x, y) = PolynomialRing(RationalField(), 2, 'xy').objgens()
    sage: f = (x^3 + 2*y^2*x)^2
    sage: g = x^2*y^2
    sage: f.gcd(g)
    x^2

Next we create the ideal :math:`(f,g)` generated by :math:`f`
and :math:`g`, by simply multiplying ``(f,g)`` by ``R`` (we could
also write ``ideal([f,g])`` or ``ideal(f,g)``).

.. link

::

    sage: I = (f, g)*R; I
    Ideal (x^6 + 4*x^4*y^2 + 4*x^2*y^4, x^2*y^2) of Multivariate Polynomial
    Ring in x, y over Rational Field
    sage: B = I.groebner_basis(); B
    [x^6, x^2*y^2]
    sage: x^2 in I
    False

Incidentally, the Gröbner basis above is not a list but an
immutable sequence. This means that it has a universe, parent, and
cannot be changed (which is good because changing the basis would
break other routines that use the Gröbner basis).

.. link

::

    sage: B.universe()
    Multivariate Polynomial Ring in x, y over Rational Field
    sage: B[1] = x
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Some (read: not as much as we would like) commutative algebra is
available in Sage, implemented via Singular. For example, we can
compute the primary decomposition and associated primes of
:math:`I`:

.. link

::

    sage: I.primary_decomposition()
    [Ideal (x^2) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y^2, x^6) of Multivariate Polynomial Ring in x, y over Rational Field]
    sage: I.associated_primes()
    [Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field]
