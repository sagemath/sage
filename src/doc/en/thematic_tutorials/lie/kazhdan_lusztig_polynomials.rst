---------------------------
Kazhdan-Lusztig Polynomials
---------------------------

Sage can compute ordinary Kazhdan-Lusztig polynomials for Weyl groups
or affine Weyl groups (and potentially other Coxeter groups).

You must create a Weyl group ``W`` and a ring containing an
indeterminate ``q``. The ring may be a univariate polynomial ring or a
univariate Laurent polynomial ring. Then you may calculate
Kazhdan-Lusztig polynomials as follows::

    sage: W = WeylGroup("A3", prefix="s")
    sage: [s1,s2,s3] = W.simple_reflections()
    sage: P.<q> = LaurentPolynomialRing(QQ)
    sage: KL = KazhdanLusztigPolynomial(W,q)
    sage: KL.R(s2, s2*s1*s3*s2)
    -1 + 3*q - 3*q^2 + q^3
    sage: KL.P(s2, s2*s1*s3*s2)
    1 + q

Thus we have the Kazhdan-Lusztig `R` and `P` polynomials.

Known algorithms for computing Kazhdan-Lusztig polynomials are highly
recursive, and caching of intermediate results is necessary for the
programs not to be prohibitively slow. Therefore intermediate results
are cached. This has the effect that as you run the program for any
given ``KazhdanLusztigPolynomial`` class, the calculations will be
slow at first but progressively faster as more polynomials are
computed.

You may see the results of the intermediate calculations by creating
the class with the option ``trace="true"``.

Since the parent of ``q`` must be a univariate ring, if you want to
work with other indeterminates, *first* create a univariate polynomial
or Laurent polynomial ring, and the Kazhdan-Lusztig class. *Then*
create a ring containing ``q`` and the other variables::

    sage: W = WeylGroup("B3", prefix="s")
    sage: [s1,s2,s3] = W.simple_reflections()
    sage: P.<q> = PolynomialRing(QQ)
    sage: KL = KazhdanLusztigPolynomial(W,q)
    sage: P1.<x,y> = PolynomialRing(P)
    sage: x*KL.P(s1*s3,s1*s3*s2*s1*s3)
    (q + 1)*x
