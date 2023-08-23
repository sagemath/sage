# -*- coding: utf-8 -*-
r"""
Field of Algebraic Numbers

AUTHOR:

- Carl Witty (2007-01-27): initial version
- Carl Witty (2007-10-29): massive rewrite to support complex as well as real numbers

This is an implementation of the algebraic numbers (the complex
numbers which are the zero of a polynomial in `\ZZ[x]`; in other
words, the algebraic closure of `\QQ`, with an embedding into `\CC`).
All computations are exact. We also include an implementation of the
algebraic reals (the intersection of the algebraic numbers with
`\RR`). The field of algebraic numbers `\QQbar` is available with
abbreviation ``QQbar``; the field of algebraic reals has abbreviation
``AA``.

As with many other implementations of the algebraic numbers, we try
hard to avoid computing a number field and working in the number
field; instead, we use floating-point interval arithmetic whenever
possible (basically whenever we need to prove non-equalities), and
resort to symbolic computation only as needed (basically to prove
equalities).

Algebraic numbers exist in one of the following forms:

- a rational number

- the sum, difference, product, or quotient of algebraic numbers

- the negation, inverse, absolute value, norm, real part,
  imaginary part, or complex conjugate of an algebraic number

- a particular root of a polynomial, given as a polynomial with
  algebraic coefficients together with an isolating interval (given as
  a ``RealIntervalFieldElement``) which encloses exactly one root, and
  the multiplicity of the root

- a polynomial in one generator, where the generator is an algebraic
  number given as the root of an irreducible polynomial with integral
  coefficients and the polynomial is given as a
  ``NumberFieldElement``.

An algebraic number can be coerced into ``ComplexIntervalField`` (or
``RealIntervalField``, for algebraic reals); every algebraic number has a
cached interval of the highest precision yet calculated.

In most cases, computations that need to compare two algebraic numbers
compute them with 128-bit precision intervals; if this does not suffice to
prove that the numbers are different, then we fall back on exact
computation.

Note that division involves an implicit comparison of the divisor against
zero, and may thus trigger exact computation.

Also, using an algebraic number in the leading coefficient of
a polynomial also involves an implicit comparison against zero, which
again may trigger exact computation.

Note that we work fairly hard to avoid computing new number fields;
to help, we keep a lattice of already-computed number fields and
their inclusions.

EXAMPLES::

    sage: sqrt(AA(2)) > 0
    True
    sage: (sqrt(5 + 2*sqrt(QQbar(6))) - sqrt(QQbar(3)))^2 == 2
    True
    sage: AA((sqrt(5 + 2*sqrt(6)) - sqrt(3))^2) == 2
    True

For a monic cubic polynomial `x^3 + bx^2 + cx + d` with roots `s1`,
`s2`, `s3`, the discriminant is defined as
`(s1-s2)^2(s1-s3)^2(s2-s3)^2` and can be computed as `b^2c^2 - 4b^3d -
4c^3 + 18bcd - 27d^2`. We can test that these definitions do give the
same result::

    sage: def disc1(b, c, d):
    ....:     return b^2*c^2 - 4*b^3*d - 4*c^3 + 18*b*c*d - 27*d^2
    sage: def disc2(s1, s2, s3):
    ....:     return ((s1-s2)*(s1-s3)*(s2-s3))^2
    sage: x = polygen(AA)
    sage: p = x*(x-2)*(x-4)
    sage: cp = AA.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = AA.polynomial_root(cp, RIF(-1, 1))
    sage: s2 = AA.polynomial_root(cp, RIF(1, 3))
    sage: s3 = AA.polynomial_root(cp, RIF(3, 5))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True
    sage: p = p + 1
    sage: cp = AA.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = AA.polynomial_root(cp, RIF(-1, 1))
    sage: s2 = AA.polynomial_root(cp, RIF(1, 3))
    sage: s3 = AA.polynomial_root(cp, RIF(3, 5))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True
    sage: p = (x-sqrt(AA(2)))*(x-AA(2).nth_root(3))*(x-sqrt(AA(3)))
    sage: cp = AA.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = AA.polynomial_root(cp, RIF(1.4, 1.5))
    sage: s2 = AA.polynomial_root(cp, RIF(1.7, 1.8))
    sage: s3 = AA.polynomial_root(cp, RIF(1.2, 1.3))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True

We can convert from symbolic expressions::

    sage: QQbar(sqrt(-5))
    2.236067977499790?*I
    sage: AA(sqrt(2) + sqrt(3))
    3.146264369941973?
    sage: QQbar(I)
    I
    sage: QQbar(I * golden_ratio)
    1.618033988749895?*I
    sage: AA(golden_ratio)^2 - AA(golden_ratio)
    1
    sage: QQbar((-8)^(1/3))
    1.000000000000000? + 1.732050807568878?*I
    sage: AA((-8)^(1/3))
    -2
    sage: QQbar((-4)^(1/4))
    1 + 1*I
    sage: AA((-4)^(1/4))
    Traceback (most recent call last):
    ...
    ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real

The coercion, however, goes in the other direction, since not all
symbolic expressions are algebraic numbers::

    sage: QQbar(sqrt(2)) + sqrt(3)
    sqrt(3) + 1.414213562373095?
    sage: QQbar(sqrt(2) + QQbar(sqrt(3)))
    3.146264369941973?

Note the different behavior in taking roots: for ``AA`` we prefer real
roots if they exist, but for ``QQbar`` we take the principal root::

    sage: AA(-1)^(1/3)
    -1
    sage: QQbar(-1)^(1/3)
    0.500000000000000? + 0.866025403784439?*I

However, implicit coercion from `\QQ[I]` is only allowed when it is equipped
with a complex embedding::

    sage: i.parent()
    Number Field in I with defining polynomial x^2 + 1 with I = 1*I
    sage: QQbar(1) + i
    I + 1

    sage: K.<im> = QuadraticField(-1, embedding=None)
    sage: QQbar(1) + im
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +: 'Algebraic Field' and
    'Number Field in im with defining polynomial x^2 + 1'

However, we can explicitly coerce from the abstract number field `\QQ[I]`.
(Technically, this is not quite kosher, since we do not know whether the field
generator is supposed to map to `+I` or `-I`. We assume that for any quadratic
field with polynomial `x^2+1`, the generator maps to `+I`.)::

    sage: pythag = QQbar(3/5 + 4*im/5); pythag
    4/5*I + 3/5
    sage: pythag.abs() == 1
    True

We can implicitly coerce from algebraic reals to algebraic numbers::

    sage: a = QQbar(1); a, a.parent()
    (1, Algebraic Field)
    sage: b = AA(1); b, b.parent()
    (1, Algebraic Real Field)
    sage: c = a + b; c, c.parent()
    (2, Algebraic Field)

Some computation with radicals::

    sage: phi = (1 + sqrt(AA(5))) / 2
    sage: phi^2 == phi + 1
    True
    sage: tau = (1 - sqrt(AA(5))) / 2
    sage: tau^2 == tau + 1
    True
    sage: phi + tau == 1
    True
    sage: tau < 0
    True

    sage: rt23 = sqrt(AA(2/3))
    sage: rt35 = sqrt(AA(3/5))
    sage: rt25 = sqrt(AA(2/5))
    sage: rt23 * rt35 == rt25
    True

The Sage rings ``AA`` and ``QQbar`` can decide equalities between radical
expressions (over the reals and complex numbers respectively)::

    sage: a = AA((2/(3*sqrt(3)) + 10/27)^(1/3) - 2/(9*(2/(3*sqrt(3)) + 10/27)^(1/3)) + 1/3)
    sage: a
    1.000000000000000?
    sage: a == 1
    True

Algebraic numbers which are known to be rational print as rationals;
otherwise they print as intervals (with 53-bit precision)::

    sage: AA(2)/3
    2/3
    sage: QQbar(5/7)
    5/7
    sage: QQbar(1/3 - 1/4*I)
    -1/4*I + 1/3
    sage: two = QQbar(4).nth_root(4)^2; two
    2.000000000000000?
    sage: two == 2; two
    True
    2
    sage: phi
    1.618033988749895?

We can find the real and imaginary parts of an algebraic number (exactly)::

    sage: r = QQbar.polynomial_root(x^5 - x - 1, CIF(RIF(0.1, 0.2), RIF(1.0, 1.1))); r
    0.1812324444698754? + 1.083954101317711?*I
    sage: r.real()
    0.1812324444698754?
    sage: r.imag()
    1.083954101317711?
    sage: r.minpoly()
    x^5 - x - 1
    sage: r.real().minpoly()
    x^10 + 3/16*x^6 + 11/32*x^5 - 1/64*x^2 + 1/128*x - 1/1024
    sage: r.imag().minpoly()  # long time (10s on sage.math, 2013)
    x^20 - 5/8*x^16 - 95/256*x^12 - 625/1024*x^10 - 5/512*x^8 - 1875/8192*x^6 + 25/4096*x^4 - 625/32768*x^2 + 2869/1048576

We can find the absolute value and norm of an algebraic number exactly.
(Note that we define the norm as the product of a number and its
complex conjugate; this is the algebraic definition of norm, if we
view ``QQbar`` as ``AA[I]``.)::

    sage: R.<x> = QQ[]
    sage: r = (x^3 + 8).roots(QQbar, multiplicities=False)[2]; r
    1.000000000000000? + 1.732050807568878?*I
    sage: r.abs() == 2
    True
    sage: r.norm() == 4
    True
    sage: (r+QQbar(I)).norm().minpoly()
    x^2 - 10*x + 13
    sage: r = AA.polynomial_root(x^2 - x - 1, RIF(-1, 0)); r
    -0.618033988749895?
    sage: r.abs().minpoly()
    x^2 + x - 1

We can compute the multiplicative order of an algebraic number::

    sage: QQbar(-1/2 + I*sqrt(3)/2).multiplicative_order()
    3
    sage: QQbar(-sqrt(3)/2 + I/2).multiplicative_order()
    12
    sage: (QQbar.zeta(23)**5).multiplicative_order()
    23

The paper "ARPREC: An Arbitrary Precision Computation Package" by
Bailey, Yozo, Li and Thompson discusses this result. Evidently it is
difficult to find, but we can easily verify it. ::

    sage: alpha = QQbar.polynomial_root(x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1, RIF(1, 1.2))
    sage: lhs = alpha^630 - 1
    sage: rhs_num = (alpha^315 - 1) * (alpha^210 - 1) * (alpha^126 - 1)^2 * (alpha^90 - 1) * (alpha^3 - 1)^3 * (alpha^2 - 1)^5 * (alpha - 1)^3
    sage: rhs_den = (alpha^35 - 1) * (alpha^15 - 1)^2 * (alpha^14 - 1)^2 * (alpha^5 - 1)^6 * alpha^68
    sage: rhs = rhs_num / rhs_den
    sage: lhs
    2.642040335819351?e44
    sage: rhs
    2.642040335819351?e44
    sage: lhs - rhs
    0.?e29
    sage: lhs == rhs
    True
    sage: lhs - rhs
    0
    sage: lhs._exact_value()
    -10648699402510886229334132989629606002223831*a^9 + 23174560249100286133718183712802529035435800*a^8 - 27259790692625442252605558473646959458901265*a^7 + 21416469499004652376912957054411004410158065*a^6 - 14543082864016871805545108986578337637140321*a^5 + 6458050008796664339372667222902512216589785*a^4 + 3052219053800078449122081871454923124998263*a^3 - 14238966128623353681821644902045640915516176*a^2 + 16749022728952328254673732618939204392161001*a - 9052854758155114957837247156588012516273410 where a^10 - a^9 + a^7 - a^6 + a^5 - a^4 + a^3 - a + 1 = 0 and a in -1.176280818259918?

Given an algebraic number, we can produce a string that will reproduce
that algebraic number if you type the string into Sage. We can see
that until exact computation is triggered, an algebraic number keeps
track of the computation steps used to produce that number::

    sage: rt2 = AA(sqrt(2))
    sage: rt3 = AA(sqrt(3))
    sage: n = (rt2 + rt3)^5; n
    308.3018001722975?
    sage: sage_input(n)
    R.<x> = AA[]
    v1 = AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))) + AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))
    v2 = v1*v1
    v2*v2*v1

But once exact computation is triggered, the computation tree is discarded,
and we get a way to produce the number directly::

    sage: n == 109*rt2 + 89*rt3
    True
    sage: sage_input(n)
    R.<y> = QQ[]
    v = AA.polynomial_root(AA.common_polynomial(y^4 - 4*y^2 + 1), RIF(-RR(1.9318516525781366), -RR(1.9318516525781364)))
    -109*v^3 + 89*v^2 + 327*v - 178

We can also see that some computations (basically, those which are
easy to perform exactly) are performed directly, instead of storing
the computation tree::

    sage: z3_3 = QQbar.zeta(3) * 3
    sage: z4_4 = QQbar.zeta(4) * 4
    sage: z5_5 = QQbar.zeta(5) * 5
    sage: sage_input(z3_3 * z4_4 * z5_5)
    R.<y> = QQ[]
    3*QQbar.polynomial_root(AA.common_polynomial(y^2 + y + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871))))*QQbar(4*I)*(5*QQbar.polynomial_root(AA.common_polynomial(y^4 + y^3 + y^2 + y + 1), CIF(RIF(RR(0.3090169943749474), RR(0.30901699437494745)), RIF(RR(0.95105651629515353), RR(0.95105651629515364)))))

Note that the ``verify=True`` argument to ``sage_input`` will always trigger
exact computation, so running ``sage_input`` twice in a row on the same number
will actually give different answers. In the following, running ``sage_input``
on ``n`` will also trigger exact computation on ``rt2``, as you can see by the
fact that the third output is different than the first::

    sage: rt2 = AA(sqrt(2))
    sage: n = rt2^2
    sage: sage_input(n, verify=True)
    # Verified
    R.<x> = AA[]
    v = AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
    v*v
    sage: sage_input(n, verify=True)
    # Verified
    AA(2)
    sage: n = rt2^2
    sage: sage_input(n, verify=True)
    # Verified
    AA(2)

Just for fun, let's try ``sage_input`` on a very complicated expression. The
output of this example changed with the rewriting of polynomial multiplication
algorithms in :trac:`10255`::

    sage: rt2 = sqrt(AA(2))
    sage: rt3 = sqrt(QQbar(3))
    sage: x = polygen(QQbar)
    sage: nrt3 = AA.polynomial_root((x-rt2)*(x+rt3), RIF(-2, -1))
    sage: one = AA.polynomial_root((x-rt2)*(x-rt3)*(x-nrt3)*(x-1-rt3-nrt3), RIF(0.9, 1.1))
    sage: one
    1.000000000000000?
    sage: sage_input(one, verify=True)
    # Verified
    R1.<x> = QQbar[]
    R2.<y> = QQ[]
    v = AA.polynomial_root(AA.common_polynomial(y^4 - 4*y^2 + 1), RIF(-RR(1.9318516525781366), -RR(1.9318516525781364)))
    AA.polynomial_root(AA.common_polynomial(x^4 + QQbar(v^3 - 3*v - 1)*x^3 + QQbar(-v^3 + 3*v - 3)*x^2 + QQbar(-3*v^3 + 9*v + 3)*x + QQbar(3*v^3 - 9*v)), RIF(RR(0.99999999999999989), RR(1.0000000000000002)))
    sage: one
    1

We can pickle and unpickle algebraic fields (and they are globally unique)::

    sage: loads(dumps(AlgebraicField())) is AlgebraicField()
    True
    sage: loads(dumps(AlgebraicRealField())) is AlgebraicRealField()
    True

We can pickle and unpickle algebraic numbers::

    sage: loads(dumps(QQbar(10))) == QQbar(10)
    True
    sage: loads(dumps(QQbar(5/2))) == QQbar(5/2)
    True
    sage: loads(dumps(QQbar.zeta(5))) == QQbar.zeta(5)
    True

    sage: t = QQbar(sqrt(2)); type(t._descr)
    <class 'sage.rings.qqbar.ANRoot'>
    sage: loads(dumps(t)) == QQbar(sqrt(2))
    True

    sage: t.exactify(); type(t._descr)
    <class 'sage.rings.qqbar.ANExtensionElement'>
    sage: loads(dumps(t)) == QQbar(sqrt(2))
    True

    sage: t = ~QQbar(sqrt(2)); type(t._descr)
    <class 'sage.rings.qqbar.ANUnaryExpr'>
    sage: loads(dumps(t)) == 1/QQbar(sqrt(2))
    True

    sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3)); type(t._descr)
    <class 'sage.rings.qqbar.ANBinaryExpr'>
    sage: loads(dumps(t)) == QQbar(sqrt(2)) + QQbar(sqrt(3))
    True

We can convert elements of ``QQbar`` and ``AA`` into the following
types: ``float``, ``complex``, ``RDF``, ``CDF``, ``RR``, ``CC``,
``RIF``, ``CIF``, ``ZZ``, and ``QQ``, with a few exceptions. (For the
arbitrary-precision types, ``RR``, ``CC``, ``RIF``, and ``CIF``, it
can convert into a field of arbitrary precision.)

Converting from ``QQbar`` to a real type (``float``, ``RDF``, ``RR``,
``RIF``, ``ZZ``, or ``QQ``) succeeds only if the ``QQbar`` is actually
real (has an imaginary component of exactly zero). Converting from
either ``AA`` or ``QQbar`` to ``ZZ`` or ``QQ`` succeeds only if the
number actually is an integer or rational. If conversion fails, a
ValueError will be raised.

Here are examples of all of these conversions::

    sage: all_vals = [AA(42), AA(22/7), AA(golden_ratio), QQbar(-13), QQbar(89/55), QQbar(-sqrt(7)), QQbar.zeta(5)]
    sage: def convert_test_all(ty):
    ....:     def convert_test(v):
    ....:         try:
    ....:             return ty(v)
    ....:         except (TypeError, ValueError):
    ....:             return None
    ....:     return [convert_test(_) for _ in all_vals]
    sage: convert_test_all(float)
    [42.0, 3.1428571428571432, 1.618033988749895, -13.0, 1.6181818181818182, -2.6457513110645907, None]
    sage: convert_test_all(complex)
    [(42+0j), (3.1428571428571432+0j), (1.618033988749895+0j), (-13+0j), (1.6181818181818182+0j), (-2.6457513110645907+0j), (0.30901699437494745+0.9510565162951536j)]
    sage: convert_test_all(RDF)
    [42.0, 3.1428571428571432, 1.618033988749895, -13.0, 1.6181818181818182, -2.6457513110645907, None]
    sage: convert_test_all(CDF)
    [42.0, 3.1428571428571432, 1.618033988749895, -13.0, 1.6181818181818182, -2.6457513110645907, 0.30901699437494745 + 0.9510565162951536*I]
    sage: convert_test_all(RR)
    [42.0000000000000, 3.14285714285714, 1.61803398874989, -13.0000000000000, 1.61818181818182, -2.64575131106459, None]
    sage: convert_test_all(CC)
    [42.0000000000000, 3.14285714285714, 1.61803398874989, -13.0000000000000, 1.61818181818182, -2.64575131106459, 0.309016994374947 + 0.951056516295154*I]
    sage: convert_test_all(RIF)
    [42, 3.142857142857143?, 1.618033988749895?, -13, 1.618181818181819?, -2.645751311064591?, None]
    sage: convert_test_all(CIF)
    [42, 3.142857142857143?, 1.618033988749895?, -13, 1.618181818181819?, -2.645751311064591?, 0.3090169943749474? + 0.9510565162951536?*I]
    sage: convert_test_all(ZZ)
    [42, None, None, -13, None, None, None]
    sage: convert_test_all(QQ)
    [42, 22/7, None, -13, 89/55, None, None]

Compute the exact coordinates of a 34-gon (the formulas used are from
Weisstein, Eric W. "Trigonometry Angles--Pi/17." and can be found at
http://mathworld.wolfram.com/TrigonometryAnglesPi17.html)::

    sage: rt17 = AA(17).sqrt()
    sage: rt2 = AA(2).sqrt()
    sage: eps = (17 + rt17).sqrt()
    sage: epss = (17 - rt17).sqrt()
    sage: delta = rt17 - 1
    sage: alpha = (34 + 6*rt17 + rt2*delta*epss - 8*rt2*eps).sqrt()
    sage: beta = 2*(17 + 3*rt17 - 2*rt2*eps - rt2*epss).sqrt()
    sage: x = rt2*(15 + rt17 + rt2*(alpha + epss)).sqrt()/8
    sage: y = rt2*(epss**2 - rt2*(alpha + epss)).sqrt()/8

    sage: cx, cy = 1, 0
    sage: for i in range(34):
    ....:    cx, cy = x*cx-y*cy, x*cy+y*cx
    sage: cx
    1.000000000000000?
    sage: cy
    0.?e-15

    sage: ax = polygen(AA)
    sage: x2 = AA.polynomial_root(256*ax**8 - 128*ax**7 - 448*ax**6 + 192*ax**5 + 240*ax**4 - 80*ax**3 - 40*ax**2 + 8*ax + 1, RIF(0.9829, 0.983))
    sage: y2 = (1-x2**2).sqrt()
    sage: x - x2
    0.?e-18
    sage: y - y2
    0.?e-17

Ideally, in the above example we should be able to test ``x == x2`` and ``y ==
y2`` but this is currently infinitely long.

TESTS:

Verify that :trac:`10981` is fixed::

    sage: x = AA['x'].gen()
    sage: P = 1/(1+x^4)
    sage: P.partial_fraction_decomposition()
    (0, [(-0.3535533905932738?*x + 1/2)/(x^2 - 1.414213562373095?*x + 1), (0.3535533905932738?*x + 1/2)/(x^2 + 1.414213562373095?*x + 1)])

Check that :trac:`22202` is fixed::

    sage: R1.<x> = AA[]; R2.<s> = QQbar[]
    sage: v = QQbar.polynomial_root(x^2 - x + 1, CIF(0.5, RIF(-0.87, -0.85)))
    sage: a = QQbar.polynomial_root((-4*v + 2)*s + (v - 1/2), CIF(RIF(0.24, 0.26), RIF(0)))
    sage: QQ(a)
    1/4

This example from :trac:`17896` should run in reasonable time, see also
:trac:`15600`::

    sage: x,y = polygens(QQ,"x,y")
    sage: p1 = x^5 + 6*x^4 - 42*x^3 - 142*x^2 + 467*x + 422
    sage: p2 = p1(x=(x-1)^2)
    sage: p3 = p2(x=x*y).resultant(p2,x).univariate_polynomial()
    sage: p4, = [f[0] for f in p3.factor() if f[0].degree() == 80]
    sage: ival = CIF((0.77, 0.78), (-0.08, -0.07))
    sage: z, = [r for r in p4.roots(QQbar, False) if r in ival]
    sage: z.exactify()

Check that :trac:`17895` is fixed. We check that ``polynomial_root``
runs under 5 seconds (used to take ~40sec)::

    sage: x,y = polygens(QQ,"x,y")
    sage: p1 = x^5 + 6*x^4 - 42*x^3 - 142*x^2 + 467*x + 422
    sage: p2 = p1(x=(x-1)^2)
    sage: p3 = p2(x=x*y).resultant(p2,x).univariate_polynomial()
    sage: p4, = [f[0] for f in p3.factor() if f[0].degree() == 80]
    sage: ival = CIF((0.77, 0.78), (-0.08, -0.07))
    sage: alarm(5.0)
    sage: z2 = QQbar.polynomial_root(p4, ival)
    sage: cancel_alarm()

Check that :trac:`28530` is fixed::

    sage: x = polygen(QQ)
    sage: K.<a> = NumberField(x^2 - x - 6256320, embedding=-2500.763730596996)
    sage: y = polygen(K)
    sage: lc = (253699680440307500000000000000000000000000*y^13 +
    ....: (-82964409970750000000000000000000000*a - 253907049983029389625000000000000000000000)*y^12 -
    ....: 1269011504560040442911087500000000000000000*y^11 +
    ....: (414989843657644100408750000000000000*a + 1270048771674262724340059170625000000000000)*y^10 +
    ....: 2539049473271641600616704837811000000000000*y^9 +
    ....: (-830315359762894607374452813100000000*a - 2541124846513368955687837282617343450000000)*y^8 -
    ....: 2540076196857756969319550626460768394380000*y^7 +
    ....: (830651117050319162421733536395645998*a + 2542152409324824242066023749434989311552001)*y^6 +
    ....: 1270551589939213408433336739488536788760000*y^5 +
    ....: (-415493479588780904355108633491291996*a - 1271590115891445566303772333517948273104002)*y^4 -
    ....: 254213042233365096819403450838768394380000*y^3 +
    ....: (83132288614462248899077910195645998*a + 254420831388756945210526696075302411552001)*y^2)
    sage: lc = lc.change_ring(QQbar)
    sage: lc.roots(CIF)
    [(-1.000505492239?, 2),
     (-1.000000000000?, 2),
     (-0.999999999662605?, 1),
     (0, 2),
     (1.000000000000?, 2),
     (1.000505492239?, 2),
     (0.999999587? + 0.?e-11*I, 1),
     (0.999999999? + 0.?e-11*I, 1)]
"""

import itertools
import operator

import sage.rings.ring
import sage.rings.abc
import sage.rings.number_field.number_field_base
from sage.misc.fast_methods import Singleton
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_string import lazy_string
from sage.misc.misc import increase_recursion_limit
from sage.structure.coerce import parent_is_numerical, parent_is_real_numerical
from sage.structure.sage_object import SageObject
from sage.structure.richcmp import (richcmp, richcmp_method,
                                    rich_to_bool, richcmp_not_equal,
                                    op_EQ, op_NE, op_GT)
from sage.rings.real_arb import RealBallField
from sage.rings.real_mpfr import RR
from sage.rings.real_mpfi import RealIntervalField, RIF, is_RealIntervalFieldElement, RealIntervalField_class
from sage.rings.cc import CC
from sage.rings.cif import CIF
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.complex_interval import is_ComplexIntervalFieldElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import NumberField, GaussianField, CyclotomicField
from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_gaussian
from sage.arith.misc import factor
from . import infinity
from sage.categories.action import Action

from sage.structure.global_options import GlobalOptions


class AlgebraicField_common(sage.rings.abc.AlgebraicField_common):
    r"""
    Common base class for the classes :class:`~AlgebraicRealField` and
    :class:`~AlgebraicField`.

    TESTS::

        sage: AA.is_finite()
        False
        sage: QQbar.is_finite()
        False
    """

    class options(GlobalOptions):
        NAME = 'AlgebraicField'
        display_format = dict(default='decimal',
                              values=dict(decimal='Always display a decimal approximation',
                                          radical='Display using radicals (if possible)'))

    def default_interval_prec(self):
        r"""
        Return the default interval precision used for root isolation.

        EXAMPLES::

            sage: AA.default_interval_prec()
            64
        """
        return 64

    def characteristic(self):
        r"""
        Return the characteristic of this field.

        Since this class is only used
        for fields of characteristic 0, this always returns 0.

        EXAMPLES::

            sage: AA.characteristic()
            0
        """
        return sage.rings.integer.Integer(0)

    def order(self):
        r"""
        Return the cardinality of ``self``.

        Since this class is only used for
        fields of characteristic 0, always returns Infinity.

        EXAMPLES::

            sage: QQbar.order()
            +Infinity
        """
        return infinity.infinity

    def common_polynomial(self, poly):
        """
        Given a polynomial with algebraic coefficients, return a
        wrapper that caches high-precision calculations and
        factorizations. This wrapper can be passed to :meth:`polynomial_root`
        in place of the polynomial.

        Using :meth:`common_polynomial` makes no semantic difference, but will
        improve efficiency if you are dealing with multiple roots
        of a single polynomial.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: p = AA.common_polynomial(x^2 - x - 1)
            sage: phi = AA.polynomial_root(p, RIF(1, 2))
            sage: tau = AA.polynomial_root(p, RIF(-1, 0))
            sage: phi + tau == 1
            True
            sage: phi * tau == -1
            True

            sage: x = polygen(SR)                                                       # optional - sage.symbolic
            sage: p = (x - sqrt(-5)) * (x - sqrt(3)); p                                 # optional - sage.symbolic
            x^2 + (-sqrt(3) - sqrt(-5))*x + sqrt(3)*sqrt(-5)
            sage: p = QQbar.common_polynomial(p)                                        # optional - sage.symbolic
            sage: a = QQbar.polynomial_root(p, CIF(RIF(-0.1, 0.1), RIF(2, 3))); a       # optional - sage.symbolic
            0.?e-18 + 2.236067977499790?*I
            sage: b = QQbar.polynomial_root(p, RIF(1, 2)); b                            # optional - sage.symbolic
            1.732050807568878?

        These "common polynomials" can be shared between real and
        complex roots::

             sage: p = AA.common_polynomial(x^3 - x - 1)
             sage: r1 = AA.polynomial_root(p, RIF(1.3, 1.4)); r1
             1.324717957244746?
             sage: r2 = QQbar.polynomial_root(p, CIF(RIF(-0.7, -0.6), RIF(0.5, 0.6))); r2
             -0.6623589786223730? + 0.5622795120623013?*I
        """
        return AlgebraicPolynomialTracker(poly)

    def _get_action_(self, G, op, self_on_left):
        """
        EXAMPLES::

            sage: coercion_model.get_action(QQbar, QQ, operator.pow)
            Right Rational Powering by Rational Field on Algebraic Field
            sage: print(coercion_model.get_action(QQ, QQbar, operator.pow))
            None
            sage: print(coercion_model.get_action(QQbar, QQ, operator.mul))
            None
            sage: coercion_model.get_action(QQbar, ZZ, operator.pow)
            Right Integer Powering by Integer Ring on Algebraic Field
        """
        if self_on_left and G is QQ and op is operator.pow:
            return AlgebraicNumberPowQQAction(G, self)

    def _factor_multivariate_polynomial(self, f, proof=True):
        r"""
        Factor the multivariate polynomial ``f``.

        INPUT:

        - ``f`` -- a multivariate polynomial defined over the algebraic field
          or the real algebraic field

        OUTPUT:

        - A factorization of ``f`` over the algebraic or real algebraic numbers
          into a unit and monic irreducible factors

        ALGORITHM:

        For rings over `\QQ`, uses Singular's ``absfact`` library.

        For rings over number fields, we reduce to the `\QQ` case by factoring
        the norm of the polynomial.

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict.factor`.

        REFERENCES:

        - [GCL1992]_ Section 8.8
        - :trac:`25390` and https://github.com/sagemath/sage/files/10659152/qqbar.pdf.gz

        .. TODO::

            Investigate whether performance can be improved by testing
            if the multivariate polynomial is actually univariate, and
            using the univariate code if so.

        TESTS::

            sage: R.<x,y> = QQbar[]
            sage: A.<u,v> = AA[]

            sage: L = QQbar._factor_multivariate_polynomial(x^2-y^2)
            sage: L
            (x - y) * (x + y)
            sage: L = QQbar._factor_multivariate_polynomial(x^2+y^2)
            sage: L
            (x + (-1*I)*y) * (x + 1*I*y)
            sage: L.value()
            x^2 + y^2

            sage: L = AA._factor_multivariate_polynomial(u^2-v^2)
            sage: L
            (u - v) * (u + v)
            sage: L = AA._factor_multivariate_polynomial(u^2+v^2)
            sage: L
            u^2 + v^2

        The test from Singular's ``absfact`` documentation::

            sage: p = (-7*x^2 + 2*x*y^2 + 6*x + y^4 + 14*y^2 + 47)*(5*x^2+y^2)^3*(x-y)^4
            sage: F = QQbar._factor_multivariate_polynomial(p)
            sage: F
            (125) * (x + (-0.4472135954999580?*I)*y)^3
            * (x + 0.4472135954999580?*I*y)^3 * (x - y)^4
            * (y^2 + (-1.828427124746191?)*x + 5.585786437626905?)
            * (y^2 + 3.828427124746190?*x + 8.414213562373095?)
            sage: F.value() == p
            True

            sage: p = (-7*u^2 + 2*u*v^2 + 6*u + v^4 + 14*v^2 + 47)*(5*u^2+v^2)^3*(u-v)^4
            sage: F = AA._factor_multivariate_polynomial(p)
            sage: F
            (125) * (u - v)^4
            * (v^2 - 1.828427124746191?*u + 5.585786437626905?)
            * (v^2 + 3.828427124746190?*u + 8.414213562373095?)
            * (u^2 + 1/5*v^2)^3
            sage: F.value() == p
            True

        A test requiring us to further extend a number field that was
        used to specify the polynomial::

            sage: p = x^2 + QQbar(sqrt(2))*y^2
            sage: F = QQbar._factor_multivariate_polynomial(p)
            sage: F
            (x + (-1.189207115002722?*I)*y) * (x + 1.189207115002722?*I*y)
            sage: F.value() == p
            True

            sage: p = u^2 + AA(sqrt(2))*v^2
            sage: F = AA._factor_multivariate_polynomial(p)
            sage: F
            u^2 + 1.414213562373095?*v^2
            sage: F.value() == p
            True

        A test requiring a number field different from the number field
        used to specify the polynomial::

            sage: p = QQbar(sqrt(2))*(x^2+y^2)
            sage: F = QQbar._factor_multivariate_polynomial(p)
            sage: F
            (1.414213562373095?) * (x + (-1*I)*y) * (x + 1*I*y)
            sage: F.value() == p
            True

            sage: p = AA(sqrt(2))*(u^2+v^2)
            sage: F = AA._factor_multivariate_polynomial(p)
            sage: F
            (1.414213562373095?) * (u^2 + v^2)
            sage: F.value() == p
            True

        A test where a factor introduces a number field that was already
        used to specify the polynomial::

            sage: p = QQbar(sqrt(2))*(x^2-2*y^2)^2
            sage: F = QQbar._factor_multivariate_polynomial(p)
            sage: F
            (1.414213562373095?)
            * (x + (-1.414213562373095?)*y)^2 * (x + 1.414213562373095?*y)^2
            sage: F.value() == p
            True

            sage: p = AA(sqrt(2))*(u^2-2*v^2)^2
            sage: F = AA._factor_multivariate_polynomial(p)
            sage: F
            (1.414213562373095?)
            * (u - 1.414213562373095?*v)^2 * (u + 1.414213562373095?*v)^2
            sage: F.value() == p
            True

        A test where two factors produce the same factor in the norm::

            sage: p = (x^2+QQbar(sqrt(2))*y^2)*(x^4-2*y^4)
            sage: F = QQbar._factor_multivariate_polynomial(p)
            sage: F
            (x + (-1.189207115002722?)*y) * (x + 1.189207115002722?*y)
            * (x + (-1.189207115002722?*I)*y)^2 * (x + 1.189207115002722?*I*y)^2
            sage: F.value() == p
            True

            sage: p = (u^2+AA(sqrt(2))*v^2)*(u^4-2*v^4)
            sage: F = AA._factor_multivariate_polynomial(p)
            sage: F
            (u - 1.189207115002722?*v) * (u + 1.189207115002722?*v)
            * (u^2 + 1.414213562373095?*v^2)^2
            sage: F.value() == p
            True

        A test where the number field that expresses the result is a subfield
        of the number field that expressed the polynomial::

            sage: p = (x^2+QQbar(2)^(1/2)*y^2)*(x+QQbar(2)^(1/8)*y)
            sage: F = QQbar._factor_multivariate_polynomial(p)
            sage: F
            (x + (-1.189207115002722?*I)*y) * (x + 1.189207115002722?*I*y)
            * (x + 1.090507732665258?*y)
            sage: F.value() == p
            True

        A test where the polynomial variable names conflict with the
        number field generator::

            sage: S.<a,b> = QQbar[]
            sage: p = a^2 + QQbar(sqrt(2))*b^2
            sage: F = QQbar._factor_multivariate_polynomial(p)
            sage: F
            (a + (-1.189207115002722?*I)*b) * (a + 1.189207115002722?*I*b)

        A test that led to :trac:`26898`::

            sage: R.<x> = QQ[]
            sage: minpoly = 4*x^7 + 27
            sage: NF.<b> = NumberField(minpoly)
            sage: for hom in NF.embeddings(QQbar):
            ....:    factor_f = (x - b).map_coefficients(hom)
            ....:    assert(minpoly % factor_f == 0)

        Test :trac:`29076`::

            sage: AA['x','y'](1).factor()   # indirect doctest
            1

        """
        from sage.structure.factorization import Factorization
        from sage.interfaces.singular import singular

        if f.degree() == 0:
            return Factorization([], f.lc())

        singular.lib('absfact.lib')

        orig_elems = f.coefficients()
        numfield, new_elems, morphism = number_field_elements_from_algebraics(orig_elems, same_field=True)

        elem_dict = dict(zip(orig_elems, new_elems))

        numfield_f = f.map_coefficients(elem_dict.__getitem__, new_base_ring=numfield)

        if numfield is not QQ:
            # We now want to compute the norm of the polynomial, i.e,
            #
            #    norm_f = prod([numfield_f.map_coefficients(h)
            #                   for h in numfield.embeddings(QQbar)])
            #
            # As nbruin pointed out during the review of Issue #25390,
            # this can be accomplished more efficiently using the resultant
            # of the polynomial with the number field's minimal polynomial.
            #
            # We use two auxiliary polynomial rings:
            #
            #   norm_ring - for the polynomial's norm (same var names,
            #               but over QQ instead of over the number field)
            #   flat_ring - all the polynomial and number field variables
            #               in a single flat polynomial ring, with the
            #               variables renamed to ensure they don't conflict
            #
            # norm_ring has the same number of generators as the polynomial
            # being factored, and flat_ring has one more.

            polynomial_gens = numfield_f.parent().gens()

            norm_ring = PolynomialRing(QQ, polynomial_gens)
            flat_ring = PolynomialRing(QQ, len(polynomial_gens) + 1, 'x')
            nf_gen = flat_ring.gen(0)

            numfield_polynomial_flat = numfield.polynomial()(nf_gen)

            polynomial_flat = sum(flat_ring({(0,)+tuple(k):1})
                                  * v.polynomial()(nf_gen)
                                  for k,v in numfield_f.dict().items())

            norm_flat = polynomial_flat.resultant(numfield_polynomial_flat, nf_gen)
            norm_f = norm_flat((0,)+norm_ring.gens())
        else:
            norm_f = numfield_f

        R = norm_f._singular_().absFactorize()

        singular.setring(R)
        L = singular('absolute_factors')

        # We're going to do some polynomial operations below that
        # require Singular to change to a different base ring, which
        # will make L "disappear".  Convert its contents now.

        factors = []

        for i in range(2, len(L[1])+1):
            factor = L[1][i].sage()
            #multiplicity = L[2][i].sage()
            minpoly = L[3][i].sage()
            factors.append((factor, minpoly))

        factorization = []

        for factor, minpoly in factors:

            # minpoly is in a multivariate polynomial ring
            # over a univariate fraction field

            assert minpoly.degree() == 0
            minpoly = minpoly.constant_coefficient()
            assert minpoly.denominator() == 1
            minpoly = minpoly.numerator()

            NF = NumberField(minpoly, minpoly.parent().gen(0))

            # Singular returns factor coefficients in a fraction field
            # of a univariate ring, which is actually a number field.

            def NF_elem_map(e):
                return NF(e.numerator()) / NF(e.denominator())

            factor_NF = factor.map_coefficients(NF_elem_map, new_base_ring=NF)
            factor_NF /= factor_NF.lc()

            # We now have a number field and a factor in that number field such
            # that the factor and all of its conjugates multiply together to
            # form a factor of the original polynomial's norm.  Each of those
            # conjugate factors may (or may not) be a factor of the original
            # polynomial, so check each one.  In the real case (AA), we expect
            # the embeddings to come in conjugate pairs, and we want to combine
            # each pair together so that the factor maps to a real polynomial.

            for hom in NF.embeddings(QQbar):
                factor_f = factor_NF.map_coefficients(hom)
                if f.base_ring() is AA:
                    target = hom(NF.gen(0))
                    if target.imag() < 0:
                        continue
                    elif target.imag() > 0:
                        conjugate_hom = NF.hom([target.conjugate()], QQbar)
                        factor_f *= factor_NF.map_coefficients(conjugate_hom)
                    factor_f = factor_f.change_ring(AA)
                for i in itertools.count(1):
                    if f % factor_f**i != 0:
                        multiplicity = i-1
                        break
                if multiplicity > 0:
                    factorization.append((factor_f, multiplicity))

        # What we'd now like to do is
        #     return Factorization(factorization, unit=QQ(L[1][1].sage()))
        # but Singular seems to have a bug and doesn't
        # always compute the unit correctly

        trial = Factorization(factorization).value()

        return Factorization(factorization, unit=f.lc() / trial.lc())


class AlgebraicRealField(Singleton, AlgebraicField_common, sage.rings.abc.AlgebraicRealField):
    r"""
    The field of algebraic reals.

    TESTS::

        sage: AA == loads(dumps(AA))
        True
    """

    def __new__(cls):
        r"""
        This method is there to ensure that pickles created before this class
        was made a :class:`~sage.misc.fast_methods.Singleton` still load.

        TESTS::

            sage: s = loads(b'x\x9cmQ\xcbR\x141\x14\xad\x11A\x083\xe2\x03T|'
            ....: b'\x82l`\xd3\xff\xe0\x86\x8de/\xba*\xcb\xa9[\xe9\xf4'
            ....: b'\xa5;e:=\'I+,\xa6J\x17B\xf9\xd7f\x08\xe2s\x95\xa4\xee9\xf7<'
            ....: b'\xf2\xe5\x8e\x0e\xaa\xe5"D?\xea8z.\x9a\x0b\xa7z\xa3I[\x15'
            ....: b'\x82\xf8\xf3\x85\xc9\xb1<xg[\xae\xbd2\xbabeO\r\xdb\x86>\x9b'
            ....: b'\xd8\x91V\x91\xdb\xc1_\xe0f\xa57\xae\r\x05P+/\xfe\xe5\x08'
            ....: b'\xaci\xa2z46\x1aG$Z\x8e*F/p\xf7oC\xa33\x18\x99</<\x07v\tf'
            ....: b'\x06\'F\xe7\xb9\x195\x0b\xacg\xc2\x8d\xbc\xe1P\x9c\xad\x04'
            ....: b'\x828\xcd\x076N\x96W\xb8WaSN\x17\xca\xa7\r9\r\xb6.+\x88Kl'
            ....: b'\x97e\xb7\x16+LO\xbeb\xb6\xc4\xfdc)\x88\xfb\x9a\x9b&\x05'
            ....: b'\xc0N)wI\x0f\xee\x13\xfbH=\xc7nh(U\xc2xP\xca\r\xd2\x8d'
            ....: b'\x8a\n\x0fK\xb9\xf5+\xfe\xa3n3MV\x98\x80\xc7rr\xfe\r\xbbr'
            ....: b'\x9bZv\xecU\x1c|\xc0\xde\x12O\xe4:\xd5*0\x9ev3\xb9C\x0b'
            ....: b'\xa3?Z\xa6\xa4\x11R6<{?I\xa2l\xb9\xbf6;\xb8\\\xc6\xe0\xb1'
            ....: b'\x9f\xb3\xf6&\xe8\xe2,\xb3R\x13\xf9\xf2\xe1\xda\x9c\xc0s'
            ....: b'\xb9\xf7?.\xe1E7\xeb\xa6W\x15^&\x80q&\x1aeo\x93Y\x13"^\xcd'
            ....: b'\xf1Z\xee\xdf\x92W\x18Z\xa4\xa6(\xd7\x867\xdf\x93\xad\x9fL'
            ....: b'\xa5W\xff\x90\x89\x07s\x1c\xfe6\xd2\x03{\xcdy\xf4v\x8e\xa3'
            ....: b'\xb1.~\x000\xc2\xe0\xa1')
            sage: s is AA
            True
        """
        try:
            return AA
        except BaseException:
            return AlgebraicField_common.__new__(cls)

    def __init__(self):
        r"""
        Standard initialization function.

        EXAMPLES:

        This function calls functions in superclasses which set the category, so we check that. ::

            sage: QQbar.category() # indirect doctest
            Category of infinite fields

        Coercions::

            sage: AA.has_coerce_map_from(ZZ)
            True
            sage: AA.has_coerce_map_from(int)
            True
        """
        from sage.categories.fields import Fields
        AlgebraicField_common.__init__(self, self, ('x',), normalize=False, category=Fields().Infinite())
        self._populate_coercion_lists_([ZZ, QQ])

    def _element_constructor_(self, x):
        r"""
        Construct an element of the field of algebraic real numbers from ``x``.

        EXAMPLES::

            sage: QQbar(sqrt(2)) in AA  # indirect doctest
            True
            sage: QQbar(I) in AA
            False
            sage: AA in AA
            False

        The following should both return ``True`` (this is a bug). ::

            sage: sqrt(2) in AA  # known bug
            False
            sage: K.<z> = CyclotomicField(5); z + 1/z in AA  # known bug
            False
        """
        if isinstance(x, AlgebraicReal):
            return x
        elif isinstance(x, AlgebraicNumber):
            if x.imag().is_zero():
                return x.real()
            else:
                raise ValueError("Cannot coerce algebraic number with non-zero imaginary part to algebraic real")
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(AA)
        return AlgebraicReal(x)

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: AA._repr_()
            'Algebraic Real Field'
        """
        return "Algebraic Real Field"

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: AA._repr_option('element_is_atomic')
            True
        """
        if key == 'element_is_atomic':
            return True
        return super()._repr_option(key)

    # Is there a standard representation for this?
    def _latex_(self):
        r"""
        Latex representation of ``self``.

        EXAMPLES::

            sage: AA._latex_()
            '\\mathbf{A}'
        """
        return "\\mathbf{A}"

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(AA, verify=True)
            # Verified
            AA
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: AA._sage_input_(SageInputBuilder(), False)
            {atomic:AA}
        """
        return sib.name('AA')

    def _coerce_map_from_(self, from_par):
        r"""
        Set up the coercion model.

        TESTS::

            sage: K.<a> = QuadraticField(7, embedding=AA(7).sqrt()); AA.has_coerce_map_from(K)
            True
            sage: a in AA
            True
            sage: a + AA(3)
            5.645751311064590?
            sage: AA.has_coerce_map_from(SR)                                            # optional - sage.symbolic
            False

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x^3 - 2, 'a', embedding=2.**(1/3))
            sage: AA.has_coerce_map_from(K)
            True
            sage: K.<s> = QuadraticField(3, embedding=-2.)
            sage: s + AA(1)
            -0.732050807568878?
            sage: K.<s> = QuadraticField(3, embedding=2.)
            sage: s + AA(1)
            2.732050807568878?
            sage: K.<s> = QuadraticField(-5)
            sage: AA.has_coerce_map_from(K)
            False
        """
        if isinstance(from_par, sage.rings.number_field.number_field_base.NumberField):
            emb = from_par.coerce_embedding()
            return emb is not None and parent_is_real_numerical(emb.codomain())

    def completion(self, p, prec, extras={}):
        r"""
        Return the completion of ``self`` at the place `p`.

        Only implemented for `p = \infty` at present.

        INPUT:

        - ``p`` -- either a prime (not implemented at present) or ``Infinity``
        - ``prec`` -- precision of approximate field to return
        - ``extras`` -- (optional) a dict of extra keyword arguments
          for the ``RealField`` constructor

        EXAMPLES::

            sage: AA.completion(infinity, 500)
            Real Field with 500 bits of precision
            sage: AA.completion(infinity, prec=53, extras={'type':'RDF'})
            Real Double Field
            sage: AA.completion(infinity, 53) is RR
            True
            sage: AA.completion(7, 10)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if p == infinity.Infinity:
            from sage.rings.real_field import create_RealField
            return create_RealField(prec, **extras)
        else:
            raise NotImplementedError

    def algebraic_closure(self):
        r"""
        Return the algebraic closure of this field, which is the field
        `\overline{\QQ}` of algebraic numbers.

        EXAMPLES::

            sage: AA.algebraic_closure()
            Algebraic Field
        """
        return QQbar

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=False):
        r"""
        Attempt to construct a homomorphism from ``self`` to ``codomain`` sending the
        generators to ``im_gens``.

        Since this field is not finitely generated,
        this cannot be implemented in a mathematically sensible way, and we
        just test that there exists a canonical coercion.

        EXAMPLES::

            sage: AA._is_valid_homomorphism_(QQbar, [QQbar(1)])
            True
            sage: AA._is_valid_homomorphism_(QQ, [QQ(1)])
            False
        """
        try:
            return im_gens[0] == codomain.coerce(self.gen(0))
        except TypeError:
            return False

    def gens(self):
        r"""
        Return a set of generators for this field.

        As this field is not
        finitely generated, we opt for just returning 1.

        EXAMPLES::

            sage: AA.gens()
            (1,)
        """
        return (self(1), )

    def gen(self, n=0):
        r"""
        Return the `n`-th element of the tuple returned by :meth:`gens`.

        EXAMPLES::

            sage: AA.gen(0)
            1
            sage: AA.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: n must be 0
        """
        if n == 0:
            return self(1)
        else:
            raise IndexError("n must be 0")

    def ngens(self):
        r"""
        Return the size of the tuple returned by :meth:`gens`.

        EXAMPLES::

            sage: AA.ngens()
            1
        """
        return 1

    def zeta(self, n=2):
        r"""
        Return an `n`-th root of unity in this field. This will raise a
        ``ValueError`` if `n \ne \{1, 2\}` since no such root exists.

        INPUT:

        - ``n`` (integer) -- default 2

        EXAMPLES::

            sage: AA.zeta(1)
            1
            sage: AA.zeta(2)
            -1
            sage: AA.zeta()
            -1
            sage: AA.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in algebraic reals

        Some silly inputs::

            sage: AA.zeta(Mod(-5, 7))
            -1
            sage: AA.zeta(0)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in algebraic reals
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        else:
            raise ValueError("no n-th root of unity in algebraic reals")

    def polynomial_root(self, poly, interval, multiplicity=1):
        r"""
        Given a polynomial with algebraic coefficients and an interval
        enclosing exactly one root of the polynomial, constructs
        an algebraic real representation of that root.

        The polynomial need not be irreducible, or even squarefree; but
        if the given root is a multiple root, its multiplicity must be
        specified. (IMPORTANT NOTE: Currently, multiplicity-`k` roots
        are handled by taking the `(k-1)`-st derivative of the polynomial.
        This means that the interval must enclose exactly one root
        of this derivative.)

        The conditions on the arguments (that the interval encloses exactly
        one root, and that multiple roots match the given multiplicity)
        are not checked; if they are not satisfied, an error may be
        thrown (possibly later, when the algebraic number is used),
        or wrong answers may result.

        Note that if you are constructing multiple roots of a single
        polynomial, it is better to use ``AA.common_polynomial`` (or
        ``QQbar.common_polynomial``; the two are equivalent) to get a
        shared polynomial.

        EXAMPLES::

            sage: x = polygen(AA)
            sage: phi = AA.polynomial_root(x^2 - x - 1, RIF(1, 2)); phi
            1.618033988749895?
            sage: p = (x-1)^7 * (x-2)
            sage: r = AA.polynomial_root(p, RIF(9/10, 11/10), multiplicity=7)
            sage: r; r == 1
            1.000000000000000?
            True
            sage: p = (x-phi)*(x-sqrt(AA(2)))
            sage: r = AA.polynomial_root(p, RIF(1, 3/2))
            sage: r; r == sqrt(AA(2))
            1.414213562373095?
            True

        We allow complex polynomials, as long as the particular root
        in question is real. ::

            sage: K.<im> = QQ[I]
            sage: x = polygen(K)
            sage: p = (im + 1) * (x^3 - 2); p
            (I + 1)*x^3 - 2*I - 2
            sage: r = AA.polynomial_root(p, RIF(1, 2)); r^3
            2.000000000000000?
        """
        if not is_RealIntervalFieldElement(interval):
            raise ValueError("interval argument of .polynomial_root on algebraic real field must be real")

        return AlgebraicReal(ANRoot(poly, interval, multiplicity))

    def random_element(self, poly_degree=2, *args, **kwds):
        r"""
        Return a random algebraic real number.

        INPUT:

        - ``poly_degree`` - default: 2 - degree of the random
          polynomial over the integers of which the returned algebraic
          real number is a (real part of a) root. This is not
          necessarily the degree of the minimal polynomial of the
          number. Increase this parameter to achieve a greater
          diversity of algebraic numbers, at a cost of greater
          computation time. You can also vary the distribution of the
          coefficients but that will not vary the degree of the
          extension containing the element.

        - ``args``, ``kwds`` - arguments and keywords passed to the random
          number generator for elements of ``ZZ``, the integers. See
          :meth:`~sage.rings.integer_ring.IntegerRing_class.random_element` for
          details, or see example below.

        OUTPUT:

        An element of ``AA``, the field of algebraic real numbers (see
        :mod:`sage.rings.qqbar`).

        ALGORITHM:

        We pass all arguments to :meth:`AlgebraicField.random_element`, and
        then take the real part of the result.

        EXAMPLES::

            sage: a = AA.random_element()
            sage: a in AA
            True

        ::

            sage: b = AA.random_element(poly_degree=5)
            sage: b in AA
            True

        Parameters for the distribution of the integer coefficients of
        the polynomials can be passed on to the random element method
        for integers. For example, we can rule out zero as a
        coefficient (and therefore as a root) by requesting
        coefficients between ``1`` and ``10``::

            sage: z = [AA.random_element(x=1, y=10) for _ in range(5)]
            sage: AA(0) in z
            False

        TESTS::

            sage: AA.random_element('junk')
            Traceback (most recent call last):
            ...
            TypeError: polynomial degree must be an integer, not junk
            sage: AA.random_element(poly_degree=0)
            Traceback (most recent call last):
            ...
            ValueError: polynomial degree must be greater than zero, not 0

        Random vectors already have a 'degree' keyword, so
        we cannot use that for the polynomial's degree::

            sage: v = random_vector(AA, degree=2, poly_degree=3)
            sage: v in AA^2
            True
        """
        return QQbar.random_element(poly_degree, *args, **kwds).real()

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the real algebraic field

        OUTPUT:

        - A factorization of ``f`` over the real algebraic numbers into a unit
          and monic irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

        TESTS::

            sage: R.<x> = AA[]
            sage: AA._factor_univariate_polynomial(x)
            x
            sage: AA._factor_univariate_polynomial(2*x)
            (2) * x
            sage: AA._factor_univariate_polynomial((x^2 + 1)^2)
            (x^2 + 1)^2
            sage: AA._factor_univariate_polynomial(x^8 + 1)
            (x^2 - 1.847759065022574?*x + 1.000000000000000?) * (x^2 - 0.7653668647301795?*x + 1.000000000000000?) * (x^2 + 0.7653668647301795?*x + 1.000000000000000?) * (x^2 + 1.847759065022574?*x + 1.000000000000000?)
            sage: AA._factor_univariate_polynomial(R(3))
            3
            sage: AA._factor_univariate_polynomial(12*x^2 - 4)
            (12) * (x - 0.5773502691896258?) * (x + 0.5773502691896258?)
            sage: AA._factor_univariate_polynomial(12*x^2 + 4)
            (12) * (x^2 + 0.3333333333333334?)
            sage: AA._factor_univariate_polynomial(EllipticCurve('11a1').change_ring(AA).division_polynomial(5))
            (5) * (x - 16.00000000000000?) * (x - 5.000000000000000?) * (x - 1.959674775249769?) * (x + 2.959674775249769?) * (x^2 - 2.854101966249685?*x + 15.47213595499958?) * (x^2 + 1.909830056250526?*x + 1.660606461254312?) * (x^2 + 3.854101966249685?*x + 6.527864045000421?) * (x^2 + 13.09016994374948?*x + 93.33939353874569?)

        """
        rr = f.roots()
        cr = [(r, e) for r, e in f.roots(QQbar) if r.imag() > 0]

        from sage.structure.factorization import Factorization
        return Factorization(
            [(f.parent()([-r, 1]), e) for r, e in rr] +
            [(f.parent()([r.norm(), -2 * r.real(), 1]), e) for r, e in cr],
            unit=f.leading_coefficient())


def is_AlgebraicRealField(F):
    r"""
    Check whether ``F`` is an :class:`~AlgebraicRealField` instance. For internal use.

    This function is deprecated. Use :func:`isinstance` with
    :class:`~sage.rings.abc.AlgebraicRealField` instead.

    EXAMPLES::

        sage: from sage.rings.qqbar import is_AlgebraicRealField
        sage: [is_AlgebraicRealField(x) for x in [AA, QQbar, None, 0, "spam"]]
        doctest:warning...
        DeprecationWarning: is_AlgebraicRealField is deprecated;
        use isinstance(..., sage.rings.abc.AlgebraicRealField instead
        See https://github.com/sagemath/sage/issues/32660 for details.
        [True, False, False, False, False]
    """
    from sage.misc.superseded import deprecation
    deprecation(32660, 'is_AlgebraicRealField is deprecated; use isinstance(..., sage.rings.abc.AlgebraicRealField instead')
    return isinstance(F, AlgebraicRealField)


# Create the globally unique AlgebraicRealField object.
AA = AlgebraicRealField()


class AlgebraicField(Singleton, AlgebraicField_common, sage.rings.abc.AlgebraicField):
    """
    The field of all algebraic complex numbers.
    """

    def __new__(cls):
        r"""
        This method is there to ensure that pickles created before this class
        was made a :class:`~sage.misc.fast_methods.Singleton` still load.

        TESTS::

            sage: s = loads(b'x\x9c}RMo\x131\x10U(-\xad\x9b\x92\x16ZJh\x80~'
            ....: b'\x00MZX~\x03\x97J\x08\xb1\x87H>F\x96\xd7;\xdd\xb1\xd8x3\xb6'
            ....: b'\x17\xe8!\x12\x1c\xda\xaa\xff\x9aI\xb7\x04\x8a*N\xb65\xef'
            ....: b'\xcd\xbc\xf7\xc6?\xee\x99\xa0\x0bHB\xf4\xb5\x89\xb5'
            ....: b'\x87$?szl\x8d2\xa5\x0eA\xdc~Q\xab/{\x1f\xca\x022\xaf\xad9'
            ....: b'\xb1P\xe6\xea\x9b\x8d\xa8\x8c\x8ePT\xfe\x8cn\xday\xeb\x8a'
            ....: b'\x90\x10e\xda\x8b\xdbxA\x0bF\xa9\xac\xb6e\xb4N)Q@\xd41zA'
            ....: b'\xf7\xff\x15R;K5(\x0f\x13\x0f\x01\x1c\xc3l\xe5D\xed<\xe4'
            ....: b'\xb5\x01A\x8b\r\xe1f\xb4\x85\x90\x9c\xce\x06\x04q\xd2\x1c'
            ....: b'\xb44\x98^\xd2\x83!-\xcb\xf6D{\xee\xd0\xb8\xa0\x95\x8b!\x89'
            ....: b'\x0bZMS\\\x88Cj\x0f~\xd2\xda\x94\x1e\xf6\xa5P0\xce \xcfY<uR'
            ....: b'\xb9\xa9L\xe5\xbe\x82\x8fj\x0c\x11\xab\\q\x14@\xeb\xa9\\R&'
            ....: b'\xd7Q\xd3F*W\xfeX\x7f\x84\xcb\\\x99a\x02=\x96\xad\x8f\xe7'
            ....: b'\xb4)WU\x01\x0e\xbc\x8e\x95\x0f\xb45\xa5\'rQe:\x00m#G\xb9;'
            ....: b'\x8ff\x08\xba\xbc+\xce\xa7\xff\x89s\xce\x11\xd4E\xf6\xf3'
            ....: b'\x8c\xfdt\xd9\xcf\x0e\xfb\xe9M\xe9y\x1f;)\xae\xa7\xb8'
            ....: b'\x91"KC\x96\xf4\xfd\x9c^ \xabx\x89\xdb\xd8\x93\x1d5\xb1'
            ....: b'\xe6K\t\x8a-\x06\x8e\x96v?\xb5\xd83\x940\xbe\xce\xaar'
            ....: b'\xcd.*O{\x8d\x8c\xb1\r&9mX\xbc\x88\xe6\xf2\xf9:\x1bA\xfbr'
            ....: b'\xeb.\xae\xa2\x03\xec\xe1\xce\xe5\x90^1\xc0:\x1b\xad.\xe7'
            ....: b'\xc1\x966Dz=\xa27\xb2;\'\xcf0j\xc2\x8bR\xcd\xd6\xe8\xf0'
            ....: b'\x8ae\xfdfj3\xfb\x06\r\xb1?\xa2\xc1_%S\x817\xd0\x94'
            ....: b'\x8eFt\\g\xc8\x96p\x0f\xf7\xf1\x00\xd7\xb0\xcd\x1a\xde"'
            ....: b'\x0f{\x87\x87W\xc8\xdc\x04\x19\xf5\xbe\xce\x92_p\'\x13\xc5')
            sage: s is QQbar
            True
        """
        try:
            return QQbar
        except BaseException:
            return AlgebraicField_common.__new__(cls)

    def __init__(self):
        r"""
        Standard init function.

        We test by setting the category::

            sage: QQbar.category() # indirect doctest
            Category of infinite fields
            sage: QQbar.base_ring()
            Algebraic Real Field

        TESTS::

            sage: QQbar._repr_option('element_is_atomic')
            False

            sage: QQbar.has_coerce_map_from(ZZ)
            True
            sage: QQbar.has_coerce_map_from(int)
            True
        """
        from sage.categories.fields import Fields
        AlgebraicField_common.__init__(self, AA, ('I',), normalize=False, category=Fields().Infinite())
        self._populate_coercion_lists_([ZZ, QQ])

    def _element_constructor_(self, x):
        """
        Try to construct an element of the field of algebraic numbers from `x`.

        EXAMPLES::

            sage: sqrt(2) in QQbar # indirect doctest
            True
            sage: 22/7 in QQbar
            True
            sage: pi in QQbar
            False
        """
        if isinstance(x, AlgebraicNumber):
            return x
        elif isinstance(x, AlgebraicReal):
            return AlgebraicNumber(x._descr)
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(QQbar)
        return AlgebraicNumber(x)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: QQbar._repr_()
            'Algebraic Field'
        """
        return "Algebraic Field"

    def _latex_(self):
        r"""
        Latex representation of ``self``.

        EXAMPLES::

            sage: QQbar._latex_()
            '\\overline{\\QQ}'
        """
        return "\\overline{\\QQ}"

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(QQbar, verify=True)
            # Verified
            QQbar
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: QQbar._sage_input_(SageInputBuilder(), False)
            {atomic:QQbar}
        """
        return sib.name('QQbar')

    def _coerce_map_from_(self, from_par):
        r"""
        Set up the coercion model.

        TESTS::

            sage: QQbar.has_coerce_map_from(AA)
            True
            sage: QQbar.has_coerce_map_from(CC)
            False
            sage: QQbar.has_coerce_map_from(SR)                                         # optional - sage.symbolic
            False

            sage: i + QQbar(2)
            I + 2
            sage: K.<ii> = QuadraticField(-1, embedding=ComplexField(13)(0,-1))
            sage: ii + QQbar(2)
            -I + 2

            sage: L.<a> = QuadraticField(-1, embedding=Zp(5).teichmuller(2))
            sage: QQbar.has_coerce_map_from(L)
            False
        """
        if from_par is AA:
            return True
        if isinstance(from_par, sage.rings.number_field.number_field_base.NumberField):
            emb = from_par.coerce_embedding()
            return emb is not None and parent_is_numerical(emb.codomain())

    def completion(self, p, prec, extras={}):
        r"""
        Return the completion of ``self`` at the place `p`.

        Only implemented for `p = \infty` at present.

        INPUT:

        - ``p`` -- either a prime (not implemented at present) or ``Infinity``
        - ``prec`` -- precision of approximate field to return
        - ``extras`` -- (optional) a dict of extra keyword arguments
          for the ``RealField`` constructor

        EXAMPLES::

            sage: QQbar.completion(infinity, 500)
            Complex Field with 500 bits of precision
            sage: QQbar.completion(infinity, prec=53, extras={'type':'RDF'})
            Complex Double Field
            sage: QQbar.completion(infinity, 53) is CC
            True
            sage: QQbar.completion(3, 20)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if p == infinity.Infinity:
            from sage.rings.real_field import create_RealField
            return create_RealField(prec, **extras).complex_field()
        else:
            raise NotImplementedError

    def algebraic_closure(self):
        """
        Return the algebraic closure of this field.

        As this field is already algebraically closed, just returns ``self``.

        EXAMPLES::

            sage: QQbar.algebraic_closure()
            Algebraic Field
        """
        return self

    def construction(self):
        """
        Return a functor that constructs ``self`` (used by the coercion machinery).

        EXAMPLES::

            sage: QQbar.construction()
            (AlgebraicClosureFunctor, Rational Field)
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        from sage.rings.rational_field import QQ
        return (AlgebraicClosureFunctor(), QQ)

    def gens(self):
        r"""
        Return a set of generators for this field.

        As this field is not
        finitely generated over its prime field, we opt for just returning I.

        EXAMPLES::

            sage: QQbar.gens()
            (I,)
        """
        return (QQbar_I,)

    def gen(self, n=0):
        r"""
        Return the `n`-th element of the tuple returned by :meth:`gens`.

        EXAMPLES::

            sage: QQbar.gen(0)
            I
            sage: QQbar.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: n must be 0
        """
        if n == 0:
            return QQbar_I
        else:
            raise IndexError("n must be 0")

    def ngens(self):
        r"""
        Return the size of the tuple returned by :meth:`gens`.

        EXAMPLES::

            sage: QQbar.ngens()
            1
        """
        return 1

    @cached_method
    def zeta(self, n=4):
        r"""
        Return a primitive `n`'th root of unity, specifically `\exp(2*\pi*i/n)`.

        INPUT:

        - ``n`` (integer) -- default 4

        EXAMPLES::

            sage: QQbar.zeta(1)
            1
            sage: QQbar.zeta(2)
            -1
            sage: QQbar.zeta(3)
            -0.500000000000000? + 0.866025403784439?*I
            sage: QQbar.zeta(4)
            I
            sage: QQbar.zeta()
            I
            sage: QQbar.zeta(5)
            0.3090169943749474? + 0.9510565162951536?*I
            sage: QQbar.zeta(3000)
            0.999997806755380? + 0.002094393571219374?*I
        """
        if n == 1:
            return self.one()
        elif n == 2:
            return -self.one()
        elif n == 4:
            return self.gen()
        else:
            nf = CyclotomicField(n)
            p = nf.polynomial()
            root = ANRoot(p, ComplexIntervalField(64).zeta(n))
            gen = AlgebraicGenerator(nf, root)
            return AlgebraicNumber(ANExtensionElement(gen, nf.gen()))

    def polynomial_root(self, poly, interval, multiplicity=1):
        r"""
        Given a polynomial with algebraic coefficients and an interval
        enclosing exactly one root of the polynomial, constructs
        an algebraic real representation of that root.

        The polynomial need not be irreducible, or even squarefree; but
        if the given root is a multiple root, its multiplicity must be
        specified. (IMPORTANT NOTE: Currently, multiplicity-`k` roots
        are handled by taking the `(k-1)`-st derivative of the polynomial.
        This means that the interval must enclose exactly one root
        of this derivative.)

        The conditions on the arguments (that the interval encloses exactly
        one root, and that multiple roots match the given multiplicity)
        are not checked; if they are not satisfied, an error may be
        thrown (possibly later, when the algebraic number is used),
        or wrong answers may result.

        Note that if you are constructing multiple roots of a single
        polynomial, it is better to use ``QQbar.common_polynomial``
        to get a shared polynomial.

        EXAMPLES::

            sage: x = polygen(QQbar)
            sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(0, 2)); phi
            1.618033988749895?
            sage: p = (x-1)^7 * (x-2)
            sage: r = QQbar.polynomial_root(p, RIF(9/10, 11/10), multiplicity=7)
            sage: r; r == 1
            1
            True
            sage: p = (x-phi)*(x-sqrt(QQbar(2)))
            sage: r = QQbar.polynomial_root(p, RIF(1, 3/2))
            sage: r; r == sqrt(QQbar(2))
            1.414213562373095?
            True
        """
        return AlgebraicNumber(ANRoot(poly, interval, multiplicity))

    def random_element(self, poly_degree=2, *args, **kwds):
        r"""
        Return a random algebraic number.

        INPUT:

        - ``poly_degree`` - default: 2 - degree of the random polynomial over
          the integers of which the returned algebraic number is a root. This
          is not necessarily the degree of the minimal polynomial of the
          number. Increase this parameter to achieve a greater diversity of
          algebraic numbers, at a cost of greater computation time. You can
          also vary the distribution of the coefficients but that will not vary
          the degree of the extension containing the element.

        - ``args``, ``kwds`` - arguments and keywords passed to the random
          number generator for elements of ``ZZ``, the integers. See
          :meth:`~sage.rings.integer_ring.IntegerRing_class.random_element` for
          details, or see example below.

        OUTPUT:

        An element of ``QQbar``, the field of algebraic numbers (see
        :mod:`sage.rings.qqbar`).

        ALGORITHM:

        A polynomial with degree between 1 and ``poly_degree``,
        with random integer coefficients is created. A root of this
        polynomial is chosen at random. The default degree is
        2 and the integer coefficients come from a distribution
        heavily weighted towards `0, \pm 1, \pm 2`.

        EXAMPLES::

            sage: a = QQbar.random_element()
            sage: a                         # random
            0.2626138748742799? + 0.8769062830975992?*I
            sage: a in QQbar
            True

            sage: b = QQbar.random_element(poly_degree=20)
            sage: b                         # random
            -0.8642649077479498? - 0.5995098147478391?*I
            sage: b in QQbar
            True

        Parameters for the distribution of the integer coefficients
        of the polynomials can be passed on to the random element method
        for integers. For example, current default behavior of this method
        returns zero about 15% of the time; if we do not include zero as a
        possible coefficient, there will never be a zero constant term, and
        thus never a zero root. ::

            sage: z = [QQbar.random_element(x=1, y=10) for _ in range(20)]
            sage: QQbar(0) in z
            False

        If you just want real algebraic numbers you can filter them out.
        Using an odd degree for the polynomials will ensure some degree of
        success. ::

            sage: r = []
            sage: while len(r) < 3:
            ....:   x = QQbar.random_element(poly_degree=3)
            ....:   if x in AA:
            ....:     r.append(x)
            sage: (len(r) == 3) and all(z in AA for z in r)
            True

        TESTS::

            sage: QQbar.random_element('junk')
            Traceback (most recent call last):
            ...
            TypeError: polynomial degree must be an integer, not junk
            sage: QQbar.random_element(poly_degree=0)
            Traceback (most recent call last):
            ...
            ValueError: polynomial degree must be greater than zero, not 0

        Random vectors already have a 'degree' keyword, so
        we cannot use that for the polynomial's degree. ::

            sage: v = random_vector(QQbar, degree=2, poly_degree=3)
            sage: v                                 # random
            (0.4694381338921299?, -0.500000000000000? + 0.866025403784439?*I)
        """
        from sage.rings.integer_ring import ZZ
        import sage.misc.prandom
        try:
            poly_degree = ZZ(poly_degree)
        except TypeError:
            msg = "polynomial degree must be an integer, not {0}"
            raise TypeError(msg.format(poly_degree))
        if poly_degree < 1:
            msg = "polynomial degree must be greater than zero, not {0}"
            raise ValueError(msg.format(poly_degree))
        R = PolynomialRing(ZZ, 'x')
        p = R.random_element(degree=poly_degree, *args, **kwds)
        # degree zero polynomials have no roots
        # totally zero poly has degree -1
        # add a random leading term
        if p.degree() < 1:
            g = R.gen(0)
            m = sage.misc.prandom.randint(1, poly_degree)
            p = p + g**m
        roots = p.roots(ring=QQbar, multiplicities=False)

        # p will have at least one root; pick one at random
        # could we instead just compute one root "randomly"?
        m = sage.misc.prandom.randint(0, len(roots) - 1)
        return roots[m]

    def _is_irreducible_univariate_polynomial(self, f):
        r"""
        Return whether ``f`` is irreducible.

        INPUT:

        - ``f`` -- a non-constant univariate polynomial defined over the
          algebraic field

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.is_irreducible`.

        EXAMPLES::

            sage: R.<x> = QQbar[]
            sage: (x^2).is_irreducible() # indirect doctest
            False

        Note that this method does not handle constant polynomials::

            sage: QQbar._is_irreducible_univariate_polynomial(R(1))
            Traceback (most recent call last):
            ...
            ValueError: polynomial must not be constant
            sage: R(1).is_irreducible()
            False
        """
        if f.degree() < 1:
            # this case is handled by the caller (PolynomialElement.is_irreducible())
            raise ValueError("polynomial must not be constant")

        return f.degree() == 1

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the algebraic field

        OUTPUT:

        - A factorization of ``f`` over the algebraic numbers into a unit and
          monic irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

        TESTS::

            sage: R.<x> = QQbar[]
            sage: QQbar._factor_univariate_polynomial(x)
            x
            sage: QQbar._factor_univariate_polynomial(2*x)
            (2) * x
            sage: QQbar._factor_univariate_polynomial((x^2 + 1)^2)
            (x - I)^2 * (x + I)^2
            sage: QQbar._factor_univariate_polynomial(x^8 - 1)
            (x - 1) * (x - 0.7071067811865475? - 0.7071067811865475?*I) * (x - 0.7071067811865475? + 0.7071067811865475?*I) * (x - I) * (x + I) * (x + 0.7071067811865475? - 0.7071067811865475?*I) * (x + 0.7071067811865475? + 0.7071067811865475?*I) * (x + 1)
            sage: QQbar._factor_univariate_polynomial(12*x^2 - 4)
            (12) * (x - 0.5773502691896258?) * (x + 0.5773502691896258?)
            sage: QQbar._factor_univariate_polynomial(R(-1))
            -1
            sage: QQbar._factor_univariate_polynomial(EllipticCurve('11a1').change_ring(QQbar).division_polynomial(5))
            (5) * (x - 16) * (x - 5) * (x - 1.959674775249769?) * (x - 1.427050983124843? - 3.665468789467727?*I) * (x - 1.427050983124843? + 3.665468789467727?*I) * (x + 0.9549150281252629? - 0.8652998037182486?*I) * (x + 0.9549150281252629? + 0.8652998037182486?*I) * (x + 1.927050983124843? - 1.677599044300515?*I) * (x + 1.927050983124843? + 1.677599044300515?*I) * (x + 2.959674775249769?) * (x + 6.545084971874737? - 7.106423590645660?*I) * (x + 6.545084971874737? + 7.106423590645660?*I)

        """
        from sage.structure.factorization import Factorization
        return Factorization([(f.parent()([-r, 1]), e) for r, e in f.roots()],
                             unit=f.leading_coefficient())


def is_AlgebraicField(F):
    r"""
    Check whether ``F`` is an :class:`~AlgebraicField` instance.

    This function is deprecated. Use :func:`isinstance` with
    :class:`~sage.rings.abc.AlgebraicField` instead.

    EXAMPLES::

        sage: from sage.rings.qqbar import is_AlgebraicField
        sage: [is_AlgebraicField(x) for x in [AA, QQbar, None, 0, "spam"]]
        doctest:warning...
        DeprecationWarning: is_AlgebraicField is deprecated;
        use isinstance(..., sage.rings.abc.AlgebraicField instead
        See https://github.com/sagemath/sage/issues/32660 for details.
        [False, True, False, False, False]
    """
    from sage.misc.superseded import deprecation
    deprecation(32660, 'is_AlgebraicField is deprecated; use isinstance(..., sage.rings.abc.AlgebraicField instead')
    return isinstance(F, AlgebraicField)


# Create the globally unique AlgebraicField object.
QQbar = AlgebraicField()


def is_AlgebraicField_common(F):
    r"""
    Check whether ``F`` is an :class:`~AlgebraicField_common` instance.

    This function is deprecated. Use :func:`isinstance` with
    :class:`~sage.rings.abc.AlgebraicField_common` instead.

    EXAMPLES::

        sage: from sage.rings.qqbar import is_AlgebraicField_common
        sage: [is_AlgebraicField_common(x) for x in [AA, QQbar, None, 0, "spam"]]
        doctest:warning...
        DeprecationWarning: is_AlgebraicField_common is deprecated;
        use isinstance(..., sage.rings.abc.AlgebraicField_common) instead
        See https://github.com/sagemath/sage/issues/32610 for details.
        [True, True, False, False, False]
    """
    from sage.misc.superseded import deprecation
    deprecation(32610, 'is_AlgebraicField_common is deprecated; use isinstance(..., sage.rings.abc.AlgebraicField_common) instead')
    return isinstance(F, AlgebraicField_common)


def prec_seq():
    r"""
    Return a generator object which iterates over an infinite increasing
    sequence of precisions to be tried in various numerical computations.

    Currently just returns powers of 2 starting at 64.

    EXAMPLES::

        sage: g = sage.rings.qqbar.prec_seq()
        sage: [next(g), next(g), next(g)]
        [64, 128, 256]
    """
    # XXX Should do some testing to see where the efficiency breaks are
    # in MPFR. We could also test variants like "bits = bits + bits // 2"
    # (I think this is what MPFR uses internally).
    bits = 64
    while True:
        yield bits
        bits = bits * 2


_short_prec_seq = (64, 128, None)


def short_prec_seq():
    r"""
    Return a sequence of precisions to try in cases when an infinite-precision
    computation is possible: returns a couple of small powers of 2 and then
    ``None``.

    EXAMPLES::

        sage: from sage.rings.qqbar import short_prec_seq
        sage: short_prec_seq()
        (64, 128, None)
    """
    return _short_prec_seq


def tail_prec_seq():
    r"""
    A generator over precisions larger than those in :func:`~short_prec_seq`.

    EXAMPLES::

        sage: from sage.rings.qqbar import tail_prec_seq
        sage: g = tail_prec_seq()
        sage: [next(g), next(g), next(g)]
        [256, 512, 1024]
    """
    bits = 256
    while True:
        yield bits
        bits = bits * 2


def rational_exact_root(r, d):
    r"""
    Check whether the rational `r` is an exact `d`'th power.

    If so, this returns the `d`'th root of `r`; otherwise, this returns ``None``.

    EXAMPLES::

        sage: from sage.rings.qqbar import rational_exact_root
        sage: rational_exact_root(16/81, 4)
        2/3
        sage: rational_exact_root(8/81, 3) is None
        True
    """
    num = r.numerator()
    den = r.denominator()

    (num_rt, num_exact) = num.nth_root(d, truncate_mode=1)
    if not num_exact:
        return None
    (den_rt, den_exact) = den.nth_root(d, truncate_mode=1)
    if not den_exact:
        return None
    return (num_rt / den_rt)


def clear_denominators(poly):
    r"""
    Take a monic polynomial and rescale the variable to get a monic
    polynomial with "integral" coefficients.

    This works on any univariate
    polynomial whose base ring has a ``denominator()`` method that returns
    integers; for example, the base ring might be `\QQ` or a number
    field.

    Returns the scale factor and the new polynomial.

    (Inspired by :pari:`primitive_pol_to_monic` .)

    We assume that coefficient denominators are "small"; the algorithm
    factors the denominators, to give the smallest possible scale factor.

    EXAMPLES::

        sage: from sage.rings.qqbar import clear_denominators

        sage: _.<x> = QQ['x']
        sage: clear_denominators(x + 3/2)
        (2, x + 3)
        sage: clear_denominators(x^2 + x/2 + 1/4)
        (2, x^2 + x + 1)

    TESTS::

        sage: R.<y> = QQ[]
        sage: coefficients_as_integer_ratios = [
        ....:     (-2774600080567517563395913264491323241652779066919616441429094563840,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (-24216324060414384566983400245979288839929814383090701293489050615808,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (325579773864372490083706670433410006284520887405882567940047555526656,
        ....:      180143564412823751787837643000565238870031999674661705128541236331583133845246252269971),
        ....:     (-86736048492777879473586471630941922517134071457946320753641122078523392,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (-2338058278498910195688689352766977573607428722429118859280880481590329344,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (105830270645785996318880019945503938356315302592627229453391693256551317504,
        ....:      1381100660498315430373421929671000164670245330839073072652149478542137359480221267403111),
        ....:     (1110926147990548796149597141538460730252912439930561079348611699181798425600,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (-89705438380888704653335165590083767769953879654958783855317882966200828559360,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (1151092895747371986483047191334923516591005329489629755485810229546333821625856,
        ....:      1381100660498315430373421929671000164670245330839073072652149478542137359480221267403111),
        ....:     (24725641793859400310483886670136079788266826658111372723121573233077840328938576,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (-31051495080139473677925068000403254349133134904365702868216464107777210775457136,
        ....:      153455628944257270041491325519000018296693925648785896961349942060237484386691251933679),
        ....:     (9431591461895130351865642769482226964622378075329823505708119342634182162193000560,
        ....:      4143301981494946291120265789013000494010735992517219217956448435626412078440663802209333),
        ....:     (1721694880863483428337378731387732043714427651970488363462560317808769716807148992,
        ....:      153455628944257270041491325519000018296693925648785896961349942060237484386691251933679),
        ....:     (255327752077837584624694974814916395144764296822788813014081161094149724325120096,
        ....:      27080405107810106477910233915117650287651869232138687699061754481218379597651397400061),
        ....:     (238105337335596176836773151768694069523377650990453522899627157538495252117232992338,
        ....:      27080405107810106477910233915117650287651869232138687699061754481218379597651397400061),
        ....:     (1255826892296350234297164500548658984205287902407560187136301197703464130999349114638,
        ....:      14336685057075938723599535602121108975815695475838128781856222960645024492874269211797),
        ....:     (1, 1)]
        sage: p = R(coefficients_as_integer_ratios)
        sage: a = QQbar.polynomial_root(
        ....:     AA.common_polynomial(p),
        ....:     CIF(RIF(-RR(0.036151142425748496), -RR(0.036151142425748489)),
        ....:     RIF(-RR(0.011298617187916445), -RR(0.011298617187916443))))
        sage: a.exactify()
        sage: a
        -0.03615114242574849? - 0.011298617187916444?*I
    """

    d = poly.denominator()
    if d == 1:
        return d, poly
    deg = poly.degree()
    denoms = [c.denominator() for c in poly]
    if all(d.nbits() < 128 for d in denoms):
        # Factor the polynomial denominators.
        factors = {}
        for i, d in enumerate(denoms):
            df = factor(d)
            for f, e in df:
                oe = 0
                if f in factors:
                    oe = factors[f]
                min_e = (e + (deg - i) - 1) // (deg - i)
                factors[f] = max(oe, min_e)
        change = 1
        for f, e in factors.items():
            change = change * f**e
    else:
        # Factoring would be too slow.
        change = poly.monic().denominator()
    poly = poly * (change**deg)
    poly = poly(poly.parent().gen() / change)
    return change, poly


def do_polred(poly, threshold=32):
    r"""
    Find a polynomial of reasonably small discriminant that generates
    the same number field as ``poly``, using the PARI ``polredbest``
    function.

    INPUT:

    - ``poly`` - a monic irreducible polynomial with integer coefficients
    - ``threshold`` - an integer used to decide whether to run ``polredbest``

    OUTPUT:

    A triple (``elt_fwd``, ``elt_back``, ``new_poly``), where:

    - ``new_poly`` is the new polynomial defining the same number field,
    - ``elt_fwd`` is a polynomial expression for a root of the new
      polynomial in terms of a root of the original polynomial,
    - ``elt_back`` is a polynomial expression for a root of the original
      polynomial in terms of a root of the new polynomial.

    EXAMPLES::

        sage: from sage.rings.qqbar import do_polred
        sage: R.<x> = QQ['x']
        sage: oldpol = x^2 - 5
        sage: fwd, back, newpol = do_polred(oldpol)
        sage: newpol
        x^2 - x - 1
        sage: Kold.<a> = NumberField(oldpol)
        sage: Knew.<b> = NumberField(newpol)
        sage: newpol(fwd(a))
        0
        sage: oldpol(back(b))
        0
        sage: do_polred(x^2 - x - 11)
        (1/3*x + 1/3, 3*x - 1, x^2 - x - 1)
        sage: do_polred(x^3 + 123456)
        (-1/4*x, -4*x, x^3 - 1929)

    This shows that :trac:`13054` has been fixed::

        sage: do_polred(x^4 - 4294967296*x^2 + 54265257667816538374400)
        (1/4*x, 4*x, x^4 - 268435456*x^2 + 211973662764908353025)
    """
    parent = poly.parent()
    bitsize = ZZ(poly[0].numerator().nbits() + poly[0].denominator().nbits())
    # time(polredbest) ≈ b²d⁵
    cost = 2 * bitsize.nbits() + 5 * poly.degree().nbits()
    if cost > threshold:
        return parent.gen(), parent.gen(), poly
    new_poly, elt_back = poly.numerator().__pari__().polredbest(flag=1)
    elt_fwd = elt_back.modreverse()
    return parent(elt_fwd.lift()), parent(elt_back.lift()), parent(new_poly)


def isolating_interval(intv_fn, pol):
    """
    ``intv_fn`` is a function that takes a precision and returns an
    interval of that precision containing some particular root of pol.
    (It must return better approximations as the precision increases.)
    pol is an irreducible polynomial with rational coefficients.

    Returns an interval containing at most one root of pol.

    EXAMPLES::

        sage: from sage.rings.qqbar import isolating_interval

        sage: _.<x> = QQ['x']
        sage: isolating_interval(lambda prec: sqrt(RealIntervalField(prec)(2)), x^2 - 2)
        1.4142135623730950488?

    And an example that requires more precision::

        sage: delta = 10^(-70)
        sage: p = (x - 1) * (x - 1 - delta) * (x - 1 + delta)
        sage: isolating_interval(lambda prec: RealIntervalField(prec)(1 + delta), p)
        1.000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000?

    The function also works with complex intervals and complex roots::

        sage: p = x^2 - x + 13/36
        sage: isolating_interval(lambda prec: ComplexIntervalField(prec)(1/2, 1/3), p)
        0.500000000000000000000? + 0.3333333333333333334?*I
    """
    dpol = pol.derivative()

    for prec in prec_seq():
        intv = intv_fn(prec)

        # We need to verify that pol has exactly one root in the
        # interval intv. We know (because it is a precondition of
        # calling this function) that it has at least one root in the
        # interval, so we only need to verify that it has at most one
        # root (that the interval is sufficiently narrow).

        # We do this by computing the derivative of the polynomial
        # over the interval. If the derivative is bounded away from zero,
        # then we know there can be at most one root.

        if not dpol(intv).contains_zero():
            return intv


def find_zero_result(fn, l):
    """
    ``l`` is a list of some sort. ``fn`` is a function which maps an element of
    ``l`` and a precision into an interval (either real or complex) of that
    precision, such that for sufficient precision, exactly one element of ``l``
    results in an interval containing 0. Returns that one element of ``l``.

    EXAMPLES::

        sage: from sage.rings.qqbar import find_zero_result
        sage: _.<x> = QQ['x']
        sage: delta = 10^(-70)
        sage: p1 = x - 1
        sage: p2 = x - 1 - delta
        sage: p3 = x - 1 + delta
        sage: p2 == find_zero_result(lambda p, prec: p(RealIntervalField(prec)(1 + delta)), [p1, p2, p3])
        True
    """
    for prec in prec_seq():
        result = None
        ambig = False
        for v in l:
            intv = fn(v, prec)
            if intv.contains_zero():
                if result is not None:
                    ambig = True
                    break
                result = v
        if ambig:
            continue
        if result is None:
            raise ValueError('find_zero_result could not find any zeroes')
        return result


def conjugate_expand(v):
    r"""
    If the interval ``v`` (which may be real or complex) includes some
    purely real numbers, return ``v'`` containing ``v`` such that
    ``v' == v'.conjugate()``. Otherwise return ``v`` unchanged. (Note that if
    ``v' == v'.conjugate()``, and ``v'`` includes one non-real root of a real
    polynomial, then ``v'`` also includes the conjugate of that root.
    Also note that the diameter of the return value is at most twice
    the diameter of the input.)

    EXAMPLES::

        sage: from sage.rings.qqbar import conjugate_expand
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(1, 2))).str(style='brackets')
        '[0.0000000000000000 .. 1.0000000000000000] + [1.0000000000000000 .. 2.0000000000000000]*I'
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(0, 1))).str(style='brackets')
        '[0.0000000000000000 .. 1.0000000000000000] + [-1.0000000000000000 .. 1.0000000000000000]*I'
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(-2, 1))).str(style='brackets')
        '[0.0000000000000000 .. 1.0000000000000000] + [-2.0000000000000000 .. 2.0000000000000000]*I'
        sage: conjugate_expand(RIF(1, 2)).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000]'
    """
    if is_RealIntervalFieldElement(v):
        return v
    im = v.imag()
    if not im.contains_zero():
        return v
    re = v.real()
    fld = ComplexIntervalField(v.prec())
    return fld(re, im.union(-im))


def conjugate_shrink(v):
    r"""
    If the interval ``v`` includes some purely real numbers, return
    a real interval containing only those real numbers. Otherwise
    return ``v`` unchanged.

    If ``v`` includes exactly one root of a real polynomial, and ``v`` was
    returned by ``conjugate_expand()``, then ``conjugate_shrink(v)`` still
    includes that root, and is a ``RealIntervalFieldElement`` iff the root
    in question is real.

    EXAMPLES::

        sage: from sage.rings.qqbar import conjugate_shrink
        sage: conjugate_shrink(RIF(3, 4)).str(style='brackets')
        '[3.0000000000000000 .. 4.0000000000000000]'
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(1, 2))).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000] + [1.0000000000000000 .. 2.0000000000000000]*I'
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(0, 1))).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000]'
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(-1, 2))).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000]'
    """
    if is_RealIntervalFieldElement(v):
        return v
    im = v.imag()
    if im.contains_zero():
        return v.real()
    return v

def number_field_elements_from_algebraics(numbers, minimal=False, same_field=False, embedded=False, prec=53):
    r"""
    Given a sequence of elements of either ``AA`` or ``QQbar``
    (or a mixture), computes a number field containing all of these
    elements, these elements as members of that number field, and a
    homomorphism from the number field back to ``AA`` or
    ``QQbar``.

    INPUT:

    - ``numbers`` -- a number or list of numbers.

    - ``minimal`` -- Boolean (default: ``False``). Whether to minimize the
      degree of the extension.

    - ``same_field`` -- Boolean (default: ``False``). See below.

    - ``embedded`` -- Boolean (default: ``False``). Whether to make the
      NumberField embedded.

    - ``prec`` -- integer (default: ``53``). The number of bit of precision
      to guarantee finding real roots.

    OUTPUT:

    A tuple with the NumberField, the numbers inside the NumberField,
    and a homomorphism from the number field back to ``AA`` or ``QQbar``.

    This may not return the smallest such number field, unless
    ``minimal=True`` is specified.

    If ``same_field=True`` is specified, and all of the elements are
    from the same field (either ``AA`` or ``QQbar``), the generated
    homomorphism will map back to that field.  Otherwise, if all specified
    elements are real, the homomorphism might map back to ``AA``
    (and will, if ``minimal=True`` is specified), even if the
    elements were in ``QQbar``.

    Also, a single number can be passed, rather than a sequence; and
    any values which are not elements of ``AA`` or ``QQbar``
    will automatically be coerced to ``QQbar``

    This function may be useful for efficiency reasons: doing exact
    computations in the corresponding number field will be faster
    than doing exact computations directly in ``AA`` or ``QQbar``.

    EXAMPLES:

    We can use this to compute the splitting field of a polynomial.
    (Unfortunately this takes an unreasonably long time for non-toy
    examples.)::

        sage: x = polygen(QQ)
        sage: p = x^3 + x^2 + x + 17
        sage: rts = p.roots(ring=QQbar, multiplicities=False)
        sage: splitting = number_field_elements_from_algebraics(rts)[0]; splitting
        Number Field in a with defining polynomial y^6 - 40*y^4 - 22*y^3 + 873*y^2 + 1386*y + 594
        sage: p.roots(ring=splitting)
        [(361/29286*a^5 - 19/3254*a^4 - 14359/29286*a^3 + 401/29286*a^2 + 18183/1627*a + 15930/1627, 1), (49/117144*a^5 - 179/39048*a^4 - 3247/117144*a^3 + 22553/117144*a^2 + 1744/4881*a - 17195/6508, 1), (-1493/117144*a^5 + 407/39048*a^4 + 60683/117144*a^3 - 24157/117144*a^2 - 56293/4881*a - 53033/6508, 1)]
        sage: rt2 = AA(sqrt(2)); rt2
        1.414213562373095?
        sage: rt3 = AA(sqrt(3)); rt3
        1.732050807568878?
        sage: rt3a = QQbar(sqrt(3)); rt3a
        1.732050807568878?
        sage: qqI = QQbar.zeta(4); qqI
        I
        sage: z3 = QQbar.zeta(3); z3
        -0.500000000000000? + 0.866025403784439?*I
        sage: rt2b = rt3 + rt2 - rt3; rt2b
        1.414213562373095?
        sage: rt2c = z3 + rt2 - z3; rt2c
        1.414213562373095? + 0.?e-19*I

        sage: number_field_elements_from_algebraics(rt2)
        (Number Field in a with defining polynomial y^2 - 2, a, Ring morphism:
           From: Number Field in a with defining polynomial y^2 - 2
           To:   Algebraic Real Field
           Defn: a |--> 1.414213562373095?)

        sage: number_field_elements_from_algebraics((rt2,rt3))
        (Number Field in a with defining polynomial y^4 - 4*y^2 + 1, [-a^3 + 3*a, a^2 - 2], Ring morphism:
            From: Number Field in a with defining polynomial y^4 - 4*y^2 + 1
            To:   Algebraic Real Field
            Defn: a |--> -1.931851652578137?)

    ``rt3a`` is a real number in ``QQbar``.  Ordinarily, we'd get a homomorphism
    to ``AA`` (because all elements are real), but if we specify ``same_field=True``,
    we'll get a homomorphism back to ``QQbar``::

        sage: number_field_elements_from_algebraics(rt3a)
        (Number Field in a with defining polynomial y^2 - 3, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 3
            To:   Algebraic Real Field
            Defn: a |--> 1.732050807568878?)

        sage: number_field_elements_from_algebraics(rt3a, same_field=True)
        (Number Field in a with defining polynomial y^2 - 3, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 3
            To:   Algebraic Field
            Defn: a |--> 1.732050807568878?)

    We've created ``rt2b`` in such a way that \sage does not initially know
    that it's in a degree-2 extension of `\QQ`::

        sage: number_field_elements_from_algebraics(rt2b)
        (Number Field in a with defining polynomial y^4 - 4*y^2 + 1, -a^3 + 3*a, Ring morphism:
            From: Number Field in a with defining polynomial y^4 - 4*y^2 + 1
            To:   Algebraic Real Field
            Defn: a |--> -1.931851652578137?)

    We can specify ``minimal=True`` if we want the smallest number field::

        sage: number_field_elements_from_algebraics(rt2b, minimal=True)
        (Number Field in a with defining polynomial y^2 - 2, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 2
            To:   Algebraic Real Field
            Defn: a |--> 1.414213562373095?)

    Things work fine with rational numbers, too::

        sage: number_field_elements_from_algebraics((QQbar(1/2), AA(17)))
        (Rational Field, [1/2, 17], Ring morphism:
            From: Rational Field
            To:   Algebraic Real Field
            Defn: 1 |--> 1)

    Or we can just pass in symbolic expressions, as long as they can be
    coerced into ``QQbar``::

        sage: number_field_elements_from_algebraics((sqrt(7), sqrt(9), sqrt(11)))
        (Number Field in a with defining polynomial y^4 - 9*y^2 + 1, [-a^3 + 8*a, 3, -a^3 + 10*a], Ring morphism:
            From: Number Field in a with defining polynomial y^4 - 9*y^2 + 1
            To:   Algebraic Real Field
            Defn: a |--> 0.3354367396454047?)

    Here we see an example of doing some computations with number field
    elements, and then mapping them back into ``QQbar``::

        sage: (fld,nums,hom) = number_field_elements_from_algebraics((rt2, rt3, qqI, z3))
        sage: fld,nums,hom  # random
        (Number Field in a with defining polynomial y^8 - y^4 + 1, [-a^5 + a^3 + a, a^6 - 2*a^2, a^6, -a^4], Ring morphism:
          From: Number Field in a with defining polynomial y^8 - y^4 + 1
          To:   Algebraic Field
          Defn: a |--> -0.2588190451025208? - 0.9659258262890683?*I)
        sage: (nfrt2, nfrt3, nfI, nfz3) = nums
        sage: hom(nfrt2)
        1.414213562373095? + 0.?e-18*I
        sage: nfrt2^2
        2
        sage: nfrt3^2
        3
        sage: nfz3 + nfz3^2
        -1
        sage: nfI^2
        -1
        sage: sum = nfrt2 + nfrt3 + nfI + nfz3; sum
        a^5 + a^4 - a^3 + 2*a^2 - a - 1
        sage: hom(sum)
        2.646264369941973? + 1.866025403784439?*I
        sage: hom(sum) == rt2 + rt3 + qqI + z3
        True
        sage: [hom(n) for n in nums] == [rt2, rt3, qqI, z3]
        True

    It is also possible to have an embedded Number Field::

        sage: x = polygen(ZZ)
        sage: my_num = AA.polynomial_root(x^3-2, RIF(0,3))
        sage: res = number_field_elements_from_algebraics(my_num,embedded=True)
        sage: res[0].gen_embedding()
        1.259921049894873?
        sage: res[2]
        Ring morphism:
          From: Number Field in a with defining polynomial y^3 - 2 with a = 1.259921049894873?
          To:   Algebraic Real Field
          Defn: a |--> 1.259921049894873?

    ::

        sage: nf,nums,hom = number_field_elements_from_algebraics([2^(1/3),3^(1/5)],embedded=True)
        sage: nf
        Number Field in a with defining polynomial y^15 - 9*y^10 + 21*y^5 - 3 with a = 0.6866813218928813?
        sage: nums
        [a^10 - 5*a^5 + 2, -a^8 + 4*a^3]
        sage: hom
        Ring morphism:
          From: Number Field in a with defining polynomial y^15 - 9*y^10 + 21*y^5 - 3 with a = 0.6866813218928813?
          To:   Algebraic Real Field
          Defn: a |--> 0.6866813218928813?

    Complex embeddings are possible as well::

        sage: elems = [sqrt(5), 2^(1/3)+sqrt(3)*I, 3/4]
        sage: nf, nums, hom = number_field_elements_from_algebraics(elems, embedded=True)
        sage: nf # random (polynomial and root not unique)
        Number Field in a with defining polynomial y^24 - 6*y^23 ...- 9*y^2 + 1
          with a = 0.2598679? + 0.0572892?*I
        sage: nf.is_isomorphic(NumberField(x^24 - 9*x^22 + 135*x^20 - 720*x^18 + 1821*x^16 - 3015*x^14 + 3974*x^12 - 3015*x^10 + 1821*x^8 - 720*x^6 + 135*x^4 - 9*x^2 + 1, 'a'))
        True
        sage: list(map(QQbar, nums)) == elems == list(map(hom, nums))
        True

    TESTS::

        sage: number_field_elements_from_algebraics(rt3)
        (Number Field in a with defining polynomial y^2 - 3, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 3
            To:   Algebraic Real Field
            Defn: a |--> 1.732050807568878?)
        sage: number_field_elements_from_algebraics((rt2,qqI))
        (Number Field in a with defining polynomial y^4 + 1, [-a^3 + a, a^2], Ring morphism:
            From: Number Field in a with defining polynomial y^4 + 1
            To:   Algebraic Field
            Defn: a |--> 0.7071067811865475? + 0.7071067811865475?*I)

    Note that for the first example, where \sage does not realize that
    the number is real, we get a homomorphism to ``QQbar``::

        sage: number_field_elements_from_algebraics(rt2c)   # random
        (Number Field in a with defining polynomial y^4 + 2*y^2 + 4, 1/2*a^3, Ring morphism:
            From: Number Field in a with defining polynomial y^4 + 2*y^2 + 4
            To:   Algebraic Field
            Defn: a |--> -0.7071067811865475? - 1.224744871391589?*I)

    But with ``minimal=True``, we get a homomorphism to ``AA``::

        sage: number_field_elements_from_algebraics(rt2c, minimal=True)
        (Number Field in a with defining polynomial y^2 - 2, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 2
            To:   Algebraic Real Field
            Defn: a |--> 1.414213562373095?)

    If we specify both ``minimal=True`` and ``same_field=True``, we get a second
    degree extension (minimal) that maps back to ``QQbar``::

        sage: number_field_elements_from_algebraics(rt2c, minimal=True, same_field=True)
        (Number Field in a with defining polynomial y^2 - 2, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 2
            To:   Algebraic Field
            Defn: a |--> 1.414213562373095?)

    Tests trivial cases::

        sage: number_field_elements_from_algebraics([], embedded=True)
        (Rational Field, [], Ring morphism:
           From: Rational Field
           To:   Algebraic Real Field
           Defn: 1 |--> 1)
        sage: number_field_elements_from_algebraics([1], embedded=True)
        (Rational Field, [1], Ring morphism:
           From: Rational Field
           To:   Algebraic Real Field
           Defn: 1 |--> 1)

    Tests more complicated combinations::

        sage: UCF = UniversalCyclotomicField()
        sage: E = UCF.gen(5)
        sage: L.<b> = NumberField(x^2-189*x+16, embedding=200)
        sage: x = polygen(ZZ)
        sage: my_nums = [-52*E - 136*E^2 - 136*E^3 - 52*E^4, \
                         L.gen()._algebraic_(AA), \
                         sqrt(2), AA.polynomial_root(x^3-3, RIF(0,3)), 11/9, 1]
        sage: res = number_field_elements_from_algebraics(my_nums, embedded=True)
        sage: res[0]
        Number Field in a with defining polynomial y^24 - 107010*y^22 - 24*y^21 + ... + 250678447193040618624307096815048024318853254384 with a = 93.32530798172420?
    """
    gen = qq_generator

    # Keep track of whether we were given a single value or a list.
    single_number = False
    try:
        len(numbers)
    except TypeError:
        numbers = [numbers]
        single_number = True

    if any(isinstance(nb, AlgebraicNumber) for nb in numbers):
        algebraic_field = QQbar
    else:
        algebraic_field = AA

    def mk_algebraic(x):
        if isinstance(x, AlgebraicNumber_base):
            return x
        return QQbar(x)

    # Try to cast into AA
    try:
        aa_numbers = [AA(_) for _ in numbers]
        numbers = aa_numbers
        real_case = True
    except (ValueError, TypeError):
        real_case = False
    # Make the numbers algebraic
    numbers = [mk_algebraic(_) for _ in numbers]

    # Make the numbers have a real exact underlying field
    real_numbers = []
    for v in numbers:
        if v._exact_field().is_complex() and real_case:
            # the number comes from a complex algebraic number field
            embedded_rt = v.interval_fast(RealIntervalField(prec))
            root = ANRoot(v.minpoly(), embedded_rt)
            real_nf = NumberField(v.minpoly(),'a')
            new_ef = AlgebraicGenerator(real_nf, root)
            real_numbers += [new_ef.root_as_algebraic()]
        else:
            real_numbers += [v]
    numbers = real_numbers

    # Get the union of the exact fields
    for v in numbers:
        if minimal:
            v.simplify()
        gen = gen.union(v._exact_field())

    fld = gen._field
    nums = [gen(v._exact_value()) for v in numbers]

    exact_generator = gen.root_as_algebraic()
    hom = fld.hom([exact_generator])

    if fld is not QQ and embedded:
        # creates the embedded field
        embedded_field = NumberField(fld.defining_polynomial(),fld.variable_name(),embedding=exact_generator)

        # embeds the numbers
        inter_hom = fld.hom([embedded_field.gen(0)])
        nums = [inter_hom(n) for n in nums]

        # get the field and homomorphism
        hom = embedded_field.hom([exact_generator])
        fld = embedded_field

    if single_number:
        nums = nums[0]

    if same_field:
        hom = fld.hom([exact_generator], codomain=algebraic_field)

    return (fld, nums, hom)


# Cache some commonly-used polynomial rings
QQx = QQ['x']
QQx_x = QQx.gen()
QQy = QQ['y']
QQy_y = QQy.gen()
QQxy = QQ['x', 'y']
QQxy_x = QQxy.gen(0)
QQxy_y = QQxy.gen(1)


def cmp_elements_with_same_minpoly(a, b, p):
    r"""
    Compare the algebraic elements ``a`` and ``b`` knowing that they have the
    same minimal polynomial ``p``.

    This is an helper function for comparison of algebraic elements (i.e. the
    methods :meth:`AlgebraicNumber._richcmp_` and
    :meth:`AlgebraicReal._richcmp_`).

    INPUT:

    - ``a`` and ``b`` -- elements of the algebraic or the real algebraic field
      with same minimal polynomial

    - ``p`` -- the minimal polynomial

    OUTPUT:

    `-1`, `0`, `1`, `None` depending on whether `a < b`, `a = b` or `a > b` or
    the function did not succeed with the given precision of `a` and `b`.

    EXAMPLES::

        sage: from sage.rings.qqbar import cmp_elements_with_same_minpoly
        sage: x = polygen(ZZ)
        sage: p = x^2 - 2
        sage: a = AA.polynomial_root(p, RIF(1,2))
        sage: b = AA.polynomial_root(p, RIF(-2,-1))
        sage: cmp_elements_with_same_minpoly(a, b, p)
        1
        sage: cmp_elements_with_same_minpoly(-a, b, p)
        0
    """
    ar = a._value.real()
    br = b._value.real()
    if not ar.overlaps(br):
        # NOTE: do not try to use "ar < br" here as it will coerce to a common
        # precision which is to be avoided. See
        # https://github.com/sagemath/sage/issues/29220
        return 1 if ar._richcmp_(br, op_GT) else -1

    ai = a._value.imag()
    bi = b._value.imag()

    if a.parent() is AA or b.parent() is AA:
        ring = AA
    else:
        ring = QQbar
    roots = p.roots(ring, False)

    real = ar.union(br)
    imag = ai.union(bi)
    oroots = [r for r in roots if r._value.real().overlaps(real)
             and r._value.imag().overlaps(imag)]
    if not oroots:
        raise RuntimeError('a = {}\nb = {}\np = {}'.format(a, b, p))
    if len(oroots) == 1:
        # There is a single root matching both descriptors
        # so they both must be that root and therefore equal.
        return 0

    # test whether we have a conjugated pair (in which situation
    # real part are equal)
    imag = ai.abs().union(bi.abs())
    oroots = [r for r in roots if r._value.real().overlaps(real)
             and r._value.imag().abs().overlaps(imag)]
    if (len(oroots) == 2 and
           not oroots[0]._value.imag().contains_zero()):
        # There is a complex conjugate pair of roots matching both
        # descriptors, so compare by imaginary value.
        while ai.contains_zero():
            a._more_precision()
            ai = a._value.imag()
        while bi.contains_zero():
            b._more_precision()
            bi = b._value.imag()
        if ai.overlaps(bi):
            return 0

        # NOTE: do not try to use "ai < bi" here as it will coerce to a common
        # precision which is to be avoided. See
        # https://github.com/sagemath/sage/issues/29220
        return 1 if ai._richcmp_(bi, op_GT) else -1

    # not able to determine equality
    return None


class AlgebraicGeneratorRelation(SageObject):
    """
    A simple class for maintaining relations in the lattice of algebraic
    extensions.
    """
    def __init__(self, child1, child1_poly, child2, child2_poly, parent):
        r"""
        EXAMPLES::

            sage: from sage.rings.qqbar import AlgebraicGeneratorRelation
            sage: AlgebraicGeneratorRelation(None, None, None, None, None)
            <sage.rings.qqbar.AlgebraicGeneratorRelation object at ...>
        """
        self.child1 = child1
        self.child1_poly = child1_poly
        self.child2 = child2
        self.child2_poly = child2_poly
        self.parent = parent


algebraic_generator_counter = 0


@richcmp_method
class AlgebraicGenerator(SageObject):
    r"""
    An ``AlgebraicGenerator`` represents both an algebraic number `\alpha` and
    the number field `\QQ[\alpha]`. There is a single ``AlgebraicGenerator``
    representing `\QQ` (with `\alpha=0`).

    The ``AlgebraicGenerator`` class is private, and should not be used
    directly.
    """

    def __init__(self, field, root):
        """
        Construct an ``AlgebraicGenerator`` object.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ, 'y')
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen
            Number Field in a with defining polynomial y^2 - y - 1 with a in 1.618033988749895?
            sage: gen.field()
            Number Field in a with defining polynomial y^2 - y - 1
            sage: gen.is_trivial()
            False
            sage: gen.union(qq_generator) is gen
            True
            sage: qq_generator.union(gen) is gen
            True
            sage: nf = NumberField(y^2 + 1, name='a', check=False)
            sage: root = ANRoot(x^2 + 1, CIF(0, 1))
            sage: x = AlgebraicGenerator(nf, root); x
            Number Field in a with defining polynomial y^2 + 1 with a in 1*I
        """
        self._field = field
        self._pari_field = None
        self._trivial = (field is QQ)
        self._root = root
        self._root_as_algebraic = (QQbar if root.is_complex() else AA)(root)
        self._unions = {}
        self._cyclotomic = False
        global algebraic_generator_counter
        self._index = algebraic_generator_counter
        algebraic_generator_counter += 1

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3))
            sage: t.exactify()
            sage: type(t._descr._generator)
            <class 'sage.rings.qqbar.AlgebraicGenerator'>
            sage: loads(dumps(t)) == t
            True
        """
        return (AlgebraicGenerator, (self._field, self._root))

    def __hash__(self):
        r"""
        Return a hash value for self. This will depend on the order that
        commands get executed at load time, so we do not test the value that is
        returned, just that it does not raise an error.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: hash(gen) # random
        """
        return self._index

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` with another ``AlgebraicGenerator`` object.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen > qq_generator
            True
        """
        return richcmp(self._index, other._index, op)

    def is_complex(self):
        r"""
        Return ``True`` if this is a generator for a non-real number field.

        EXAMPLES::

            sage: z7 = QQbar.zeta(7)
            sage: g = z7._descr._generator
            sage: g.is_complex()
            True

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator
            sage: y = polygen(QQ, 'y')
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen.is_complex()
            False
        """
        return self._root.is_complex()

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator._repr_()
            'Trivial generator'

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ)
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen._repr_()
            'Number Field in a with defining polynomial x^2 - x - 1 with a in 1.618033988749895?'
        """
        if self._trivial:
            return 'Trivial generator'
        else:
            return '%s with a in %s' % (self._field,
                                        self._root._interval_fast(53))

    def root_as_algebraic(self):
        r"""
        Return the root attached to self as an algebraic number.

        EXAMPLES::

            sage: t = sage.rings.qqbar.qq_generator.root_as_algebraic(); t
            1
            sage: t.parent()
            Algebraic Real Field
        """
        return self._root_as_algebraic

    def is_trivial(self):
        """
        Return ``True`` iff this is the trivial generator (alpha == 1), which
        does not actually extend the rationals.

        EXAMPLES::

            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator.is_trivial()
            True
        """
        return self._trivial

    def field(self):
        r"""
        Return the number field attached to self.

        EXAMPLES::

            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator.field()
            Rational Field
        """
        return self._field

    def pari_field(self):
        r"""
        Return the PARI field attached to this generator.

        EXAMPLES::


            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator.pari_field()
            Traceback (most recent call last):
            ...
            ValueError: No PARI field attached to trivial generator

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ)
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen.pari_field()
            [[y^2 - y - 1, [2, 0], ...]
        """
        if self.is_trivial():
            raise ValueError("No PARI field attached to trivial generator")
        if self._pari_field is None:
            pari_pol = self._field.pari_polynomial("y")
            self._pari_field = pari_pol.nfinit(1)
        return self._pari_field

    def conjugate(self):
        r"""
        If this generator is for the algebraic number `\alpha`, return a
        generator for the complex conjugate of `\alpha`.

        EXAMPLES::

            sage: from sage.rings.qqbar import AlgebraicGenerator
            sage: x = polygen(QQ); f = x^4 + x + 17
            sage: nf = NumberField(f,name='a')
            sage: b = f.roots(QQbar)[0][0]
            sage: root = b._descr
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen.conjugate()
            Number Field in a with defining polynomial x^4 + x + 17 with a in -1.436449997483091? + 1.374535713065812?*I
        """
        try:
            return self._conjugate
        except AttributeError:
            if not self.is_complex():
                self._conjugate = self
            else:
                conj = AlgebraicGenerator(self._field, self._root.conjugate(None))
                self._conjugate = conj
                conj._conjugate = self
            if self._cyclotomic:
                conj_rel = QQx_x ** (self._cyclotomic_order - 1)
                rel = AlgebraicGeneratorRelation(self, conj_rel, conj, conj_rel, self)
                self._unions[conj] = rel
                conj._unions[self] = rel
            return self._conjugate

    def _interval_fast(self, prec):
        """
        Return an interval containing this generator, to the specified
        precision.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ, 'y')
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen._interval_fast(128)
            1.61803398874989484820458683436563811773?
        """
        return self._root._interval_fast(prec)

    def union(self, other):
        r""" Given generators ``alpha`` and ``beta``,
        ``alpha.union(beta)`` gives a generator for the number field
        `\QQ[\alpha][\beta]`.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in 1.414213562373095?
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in 1.732050807568878?
            sage: gen2.union(qq_generator) is gen2
            True
            sage: qq_generator.union(gen3) is gen3
            True
            sage: gen2.union(gen3)
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in -1.931851652578137?
        """
        if self._trivial:
            return other
        if other._trivial:
            return self
        if self is other:
            return self
        if other in self._unions:
            return self._unions[other].parent
        if self._field.polynomial().degree() < other._field.polynomial().degree():
            self, other = other, self
        elif other._cyclotomic:
            self, other = other, self

        op = other._field.polynomial()
        op = QQx(op)
        # pari_nf = self._field.pari_nf()
        pari_nf = self.pari_field()
        factors = list(pari_nf.nffactor(op).lift())[0]
        x, y = QQxy.gens()
        factors_sage = [QQxy(p) for p in factors]

        def find_fn(p, prec):
            ifield = RealIntervalField(prec)
            if_poly = ifield['x', 'y']
            ip = if_poly(p)
            return ip(other._root._interval_fast(prec), self._root._interval_fast(prec))
        my_factor = find_zero_result(find_fn, factors_sage)

        if my_factor.degree(x) == 1 and my_factor.coefficient(x) == 1:
            value = (-my_factor + x).univariate_polynomial(QQy)
            rel = AlgebraicGeneratorRelation(self, QQy_y,
                                             other, value,
                                             self)
            self._unions[other] = rel
            other._unions[self] = rel
            return rel.parent

        # XXX need more caching here
        newpol, self_pol, k = pari_nf.rnfequation(my_factor, 1)
        k = int(k)

        newpol_sage = QQx(newpol)
        newpol_sage_y = QQy(newpol_sage)

        red_elt, red_back, red_pol = do_polred(newpol_sage_y)

        red_back_x = QQx(red_back)

        new_nf = NumberField(red_pol, name='a', check=False)

        self_pol_sage = QQx(self_pol.lift())

        def intv_fn(prec):
            return conjugate_expand(red_elt(self._root._interval_fast(prec) * k + other._root._interval_fast(prec)))
        new_intv = conjugate_shrink(isolating_interval(intv_fn, red_pol))

        new_gen = AlgebraicGenerator(new_nf, ANRoot(QQx(red_pol), new_intv))
        rel = AlgebraicGeneratorRelation(self, self_pol_sage(red_back_x),
                                         other,
                                         (QQx_x - k * self_pol_sage)(red_back_x),
                                         new_gen)
        self._unions[other] = rel
        other._unions[self] = rel
        return new_gen

    def super_poly(self, super, checked=None):
        r"""
        Given a generator ``gen`` and another generator ``super``, where ``super``
        is the result of a tree of ``union()`` operations where one of the
        leaves is ``gen``, ``gen.super_poly(super)`` returns a polynomial
        expressing the value of ``gen`` in terms of the value of ``super``
        (except that if ``gen`` is ``qq_generator``, ``super_poly()`` always
        returns None.)

        EXAMPLES::

            sage: from sage.rings.qqbar import AlgebraicGenerator, ANRoot, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in 1.414213562373095?
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in 1.732050807568878?
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in -1.931851652578137?
            sage: qq_generator.super_poly(gen2) is None
            True
            sage: gen2.super_poly(gen2_3)
            -a^3 + 3*a
            sage: gen3.super_poly(gen2_3)
            a^2 - 2

        """
        if checked is None:
            checked = {}
        checked[self] = True
        if super is self:
            return self._field.gen()
        for u in self._unions.values():
            if u.parent in checked:
                continue
            poly = u.parent.super_poly(super, checked)
            if poly is None:
                continue
            if self is u.child1:
                return u.child1_poly(poly)
            assert(self is u.child2)
            return u.child2_poly(poly)
        return None

    def __call__(self, elt):
        """
        Takes an algebraic number which is represented as either a
        rational or a number field element, and which is in a subfield
        of the field generated by this generator. Lifts the number
        into the field of this generator, and returns either a
        ``Rational`` or a ``NumberFieldElement`` depending on whether
        this is the trivial generator.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, ANExtensionElement, ANRational
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in 1.414213562373095?
            sage: sqrt2 = ANExtensionElement(gen2, nf2.gen())
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in 1.732050807568878?
            sage: sqrt3 = ANExtensionElement(gen3, nf3.gen())
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in -1.931851652578137?
            sage: gen2_3(sqrt2)
            -a^3 + 3*a
            sage: gen2_3(ANRational(1/7))
            1/7
            sage: gen2_3(sqrt3)
            a^2 - 2
        """
        if self._trivial:
            return elt._value
        if isinstance(elt, ANRational):
            return self._field(elt._value)
        if elt.generator() is self:
            return elt.field_element_value()
        gen = elt.generator()
        sp = gen.super_poly(self)
        assert(not(sp is None))
        return self._field(elt.field_element_value().polynomial()(sp))


# dictionary for multimethod dispatch
_binop_algo = {}


class ANDescr(SageObject):
    r"""
    An ``AlgebraicNumber`` or ``AlgebraicReal`` is a wrapper around an
    ``ANDescr`` object. ``ANDescr`` is an abstract base class, which should
    never be directly instantiated; its concrete subclasses are ``ANRational``,
    ``ANBinaryExpr``, ``ANUnaryExpr``, ``ANRoot``, and ``ANExtensionElement``.
    ``ANDescr`` and all of its subclasses are for internal use, and should not
    be used directly.
    """
    def is_simple(self):
        r"""
        Check whether this descriptor represents a value with the same
        algebraic degree as the number field associated with the descriptor.

        This returns ``True`` if self is an ``ANRational``, or a minimal
        ``ANExtensionElement``.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRational
            sage: ANRational(1/2).is_simple()
            True
            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt2b = rt3 + rt2 - rt3
            sage: rt2.exactify()
            sage: rt2._descr.is_simple()
            True
            sage: rt2b.exactify()
            sage: rt2b._descr.is_simple()
            False
            sage: rt2b.simplify()
            sage: rt2b._descr.is_simple()
            True
        """
        return False

    # Unitary operators: the second argument "n" is an AlgebraicNumber_base
    # wrapper around self.

    def neg(self, n):
        r"""
        Negation of self.

        EXAMPLES::

            sage: a = QQbar(sqrt(2))
            sage: b = a._descr
            sage: b.neg(a)
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        return ANUnaryExpr(n, '-')

    def invert(self, n):
        r"""
        1/self.

        EXAMPLES::

            sage: a = QQbar(sqrt(2))
            sage: b = a._descr
            sage: b.invert(a)
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        return ANUnaryExpr(n, '~')

    def abs(self, n):
        r"""
        Absolute value of self.

        EXAMPLES::

            sage: a = QQbar(sqrt(2))
            sage: b = a._descr
            sage: b.abs(a)
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        return ANUnaryExpr(n, 'abs')

    def real(self, n):
        r"""
        Real part of self.

        EXAMPLES::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.real(a)
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'real')
        else:
            return self

    def imag(self, n):
        r"""
        Imaginary part of self.

        EXAMPLES::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.imag(a)
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'imag')
        else:
            return ANRational(0)

    def conjugate(self, n):
        r"""
        Complex conjugate of self.

        EXAMPLES::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.conjugate(a)
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'conjugate')
        else:
            return self

    def norm(self, n):
        r"""
        Field norm of self from `\overline{\QQ}` to its real subfield
        `\mathbf{A}`, i.e.~the square of the usual complex absolute value.

        EXAMPLES::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.norm(a)
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'norm')
        else:
            return (n * n)._descr


class AlgebraicNumber_base(sage.structure.element.FieldElement):
    r"""
    This is the common base class for algebraic numbers (complex
    numbers which are the zero of a polynomial in `\ZZ[x]`) and algebraic
    reals (algebraic numbers which happen to be real).

    ``AlgebraicNumber`` objects can be created using ``QQbar`` (==
    ``AlgebraicNumberField()``), and ``AlgebraicReal`` objects can be created
    using ``AA`` (== ``AlgebraicRealField()``). They can be created either by
    coercing a rational or a symbolic expression, or by using the
    ``QQbar.polynomial_root()`` or ``AA.polynomial_root()`` method to
    construct a particular root of a polynomial with algebraic
    coefficients. Also, ``AlgebraicNumber`` and ``AlgebraicReal`` are closed
    under addition, subtraction, multiplication, division (except by
    0), and rational powers (including roots), except that for a
    negative ``AlgebraicReal``, taking a power with an even denominator returns
    an ``AlgebraicNumber`` instead of an ``AlgebraicReal``.

    ``AlgebraicNumber``   and   ``AlgebraicReal``   objects   can   be
    approximated  to  any desired  precision. They  can be  compared
    exactly; if the two numbers are very close, or are equal, this may
    require exact computation, which can be extremely slow.

    As long as exact computation is not triggered, computation with
    algebraic numbers should not be too much slower than computation with
    intervals. As mentioned above, exact computation is triggered
    when comparing two algebraic numbers which are very close together.
    This can be an explicit comparison in user code, but the following
    list of actions (not necessarily complete) can also trigger exact
    computation:

    - Dividing by an algebraic number which is very close to 0.

    - Using an algebraic number which is very close to 0 as the leading
      coefficient in a polynomial.

    - Taking a root of an algebraic number which is very close to 0.

    The exact definition of "very close" is subject to change; currently,
    we compute our best approximation of the two numbers using 128-bit
    arithmetic, and see if that's sufficient to decide the comparison.
    Note that comparing two algebraic numbers which are actually equal will
    always trigger exact computation, unless they are actually the same object.

    EXAMPLES::

        sage: sqrt(QQbar(2))
        1.414213562373095?
        sage: sqrt(QQbar(2))^2 == 2
        True
        sage: x = polygen(QQbar)
        sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(1, 2))
        sage: phi
        1.618033988749895?
        sage: phi^2 == phi+1
        True
        sage: AA(sqrt(65537))
        256.0019531175495?
    """

    def __init__(self, parent, x):
        r"""
        Initialize an algebraic number. The argument must be either
        a rational number, a Gaussian rational, or a subclass of ``ANDescr``.

        EXAMPLES::

            sage: AlgebraicReal(22/7)
            22/7
        """
        sage.structure.element.FieldElement.__init__(self, parent)
        if isinstance(x, (int, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            self._descr = ANRational(x)
        elif isinstance(x, ANDescr):
            self._descr = x
        elif parent is QQbar and isinstance(x, NumberFieldElement_gaussian):
            if x.parent()._standard_embedding:
                self._descr = ANExtensionElement(QQbar_I_generator, QQbar_I_nf(x.list()))
            else:
                self._descr = ANExtensionElement(QQbar_I_generator, QQbar_I_nf([x[0], -x[1]]))
        else:
            raise TypeError("Illegal initializer for algebraic number")

        prec = 64
        self._value = self._descr._interval_fast(prec)
        while self._value.is_NaN():
            prec = 2 * prec
            self._value = self._descr._interval_fast(prec)

    def _repr_(self):
        """
        Return the print representation of this number.

        :class:`AlgebraicField_common`'s option display_format controls
        whether irrational numbers will always be printed using a decimal
        approximation (display_format = 'decimal'), or whether an attempt
        will be made to print them as radicals (display_format = 'radical')

        EXAMPLES::

            sage: AA(22/7) # indirect doctest
            22/7
            sage: QQbar(1/3 + 2/7*I)
            2/7*I + 1/3
            sage: QQbar.zeta(4) + 5
            I + 5
            sage: QQbar.zeta(4)
            I
            sage: 3*QQbar.zeta(4)
            3*I
            sage: QQbar.zeta(17)
            0.9324722294043558? + 0.3612416661871530?*I
            sage: AA(19).sqrt()
            4.358898943540674?
            sage: AA.options.display_format = 'radical'
            sage: AA(19).sqrt()
            sqrt(19)
            sage: QQbar.zeta(6)
            1/2*I*sqrt(3) + 1/2
            sage: QQbar.zeta(17)
            0.9324722294043558? + 0.3612416661871530?*I
            sage: AA.options.display_format = 'decimal'
        """
        if isinstance(self._descr, ANRational):
            return repr(self._descr)
        if isinstance(self._descr, ANExtensionElement) and self._descr._generator is QQbar_I_generator:
            return repr(self._descr._value)
        if self.parent().options.display_format == 'radical':
            radical = self.radical_expression()
            if radical is not self:
                return repr(radical)
        if self.parent() is QQbar:
            return repr(CIF(self._value))
        else:
            return repr(RIF(self._value))

    def _latex_(self):
        r"""
        Return the latex representation of this number.

        EXAMPLES::

            sage: latex(AA(22/7))
            \frac{22}{7}
            sage: latex(QQbar(1/3 + 2/7*I))
            \frac{2}{7} i + \frac{1}{3}
            sage: latex(QQbar.zeta(4) + 5)
            i + 5
            sage: latex(QQbar.zeta(4))
            i
            sage: latex(3*QQbar.zeta(4))
            3 i
            sage: latex(QQbar.zeta(17))
            0.9324722294043558? + 0.3612416661871530? \sqrt{-1}
            sage: latex(AA(19).sqrt())
            4.358898943540674?
            sage: AA.options.display_format = 'radical'
            sage: latex(AA(19).sqrt())
            \sqrt{19}
            sage: latex(QQbar.zeta(6))
            \frac{1}{2} i \, \sqrt{3} + \frac{1}{2}
            sage: latex(QQbar.zeta(17))
            0.9324722294043558? + 0.3612416661871530? \sqrt{-1}
            sage: AA.options.display_format = 'decimal'
        """
        from sage.misc.latex import latex
        if isinstance(self._descr, ANRational):
            return latex(self._descr._value)
        if isinstance(self._descr, ANExtensionElement) and self._descr._generator is QQbar_I_generator:
            return latex(self._descr._value)
        if self.parent().options.display_format == 'radical':
            radical = self.radical_expression()
            if radical is not self:
                return latex(radical)
        return repr(self).replace('*I', r' \sqrt{-1}')

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES:

        These examples are mostly copied from the doctests of
        the ``handle_sage_input`` functions; see those for more examples::

            sage: sage_input(QQbar(3))
            QQbar(3)
            sage: sage_input(AA(22/7))
            AA(22/7)
            sage: sage_input(22/7*QQbar.zeta(4))
            QQbar(22/7*I)
            sage: sage_input(QQbar.zeta(5)^3)
            R.<y> = QQ[]
            QQbar.polynomial_root(AA.common_polynomial(y^4 + y^3 + y^2 + y + 1), CIF(RIF(RR(0.3090169943749474), RR(0.30901699437494745)), RIF(RR(0.95105651629515353), RR(0.95105651629515364))))^3
            sage: sage_input((AA(3)^(1/2))^(1/3))
            R.<x> = AA[]
            AA.polynomial_root(AA.common_polynomial(x^3 - AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))), RIF(RR(1.2009369551760025), RR(1.2009369551760027)))
            sage: sage_input(QQbar(3+4*I))
            QQbar(3 + 4*I)
            sage: sage_input(-sqrt(AA(2)))
            R.<x> = AA[]
            -AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2 + sqrt(AA(2)))
            R.<x> = AA[]
            2 + AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))

        And a nice big example::

            sage: K.<x> = QQ[]
            sage: p = K(-12*x^3 + 1/2*x^2 - 1/95*x - 1/2)
            sage: rts = p.roots(ring=QQbar, multiplicities=False); rts
            [-0.3325236940280402?, 0.1870951803473535? - 0.3004991638609601?*I, 0.1870951803473535? + 0.3004991638609601?*I]
            sage: sage_input(rts, verify=True)
            # Verified
            R.<y> = QQ[]
            cp = AA.common_polynomial(-12*y^3 + 1/2*y^2 - 1/95*y - 1/2)
            [QQbar.polynomial_root(cp, CIF(RIF(-RR(0.33252369402804022), -RR(0.33252369402804016)), RIF(RR(0)))), QQbar.polynomial_root(cp, CIF(RIF(RR(0.18709518034735342), RR(0.18709518034735345)), RIF(-RR(0.30049916386096009), -RR(0.30049916386096004)))), QQbar.polynomial_root(cp, CIF(RIF(RR(0.18709518034735342), RR(0.18709518034735345)), RIF(RR(0.30049916386096004), RR(0.30049916386096009))))]

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sqrt(QQbar(7))._sage_input_(sib, False)
            {call: {getattr: {atomic:QQbar}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:QQbar}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:7}})}, {call: {atomic:CIF}({call: {atomic:RIF}({call: {atomic:RR}({atomic:2.6457513110645903})}, {call: {atomic:RR}({atomic:2.6457513110645907})})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:0})})})})}
        """
        (v, complicated) = \
            self._descr.handle_sage_input(sib, coerce, self.parent() is QQbar)
        if complicated or True:
            sib.id_cache(self, v, 'v')
        return v

    def _mul_(self, other):
        """
        TESTS::

            sage: AA(sqrt(2)) * AA(sqrt(8)) # indirect doctest
            4.000000000000000?
        """
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk, ok](self, other, operator.mul))

    def _div_(self, other):
        """
        TESTS::

            sage: AA(sqrt(2)) / AA(sqrt(8)) # indirect doctest
            0.500000000000000?

            sage: z = QQbar(I).real()
            sage: 1 / z
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in algebraic field
        """
        if not other:
            raise ZeroDivisionError("division by zero in algebraic field")
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk, ok](self, other, operator.truediv))

    def __invert__(self):
        """
        TESTS::

            sage: ~AA(sqrt(~2))
            1.414213562373095?

            sage: z = QQbar(I).real()
            sage: a = ~z
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in algebraic field
        """
        if not self:
            raise ZeroDivisionError("division by zero in algebraic field")
        return type(self)(self._descr.invert(self))

    def _add_(self, other):
        """
        TESTS::

            sage: x = polygen(ZZ)
            sage: rt1, rt2 = (x^2 - x - 1).roots(ring=AA, multiplicities=False)
            sage: rt1 + rt2 # indirect doctest
            1.000000000000000?
        """
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk, ok](self, other, operator.add))

    def _sub_(self, other):
        """
        TESTS::

            sage: AA(golden_ratio) * 2 - AA(5).sqrt() # indirect doctest
            1.000000000000000?
        """
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk, ok](self, other, operator.sub))

    def _neg_(self):
        """
        TESTS::

            sage: -QQbar(I) # indirect doctest
            -I
        """
        return type(self)(self._descr.neg(self))

    def __abs__(self):
        """
        TESTS::

            sage: abs(AA(sqrt(2) - sqrt(3)))
            0.3178372451957823?
            sage: abs(QQbar(3+4*I))
            5
            sage: v = QQbar.zeta(3) + 1
            sage: v.exactify()
            sage: v.abs().minpoly()
            x - 1
        """
        return AlgebraicReal(self._descr.abs(self))

    def __hash__(self):
        """
        Compute a hash code for this number (equal algebraic numbers will
        have the same hash code, different algebraic numbers are likely
        to have different hash codes).

        This may trigger exact computation, but that is very unlikely.

        TESTS:

        The hash code is stable, even when the representation changes::

            sage: two = QQbar(4).nth_root(4)^2
            sage: two
            2.000000000000000?
            sage: h1 = hash(two)
            sage: two == 2
            True
            sage: two
            2
            sage: h2 = hash(two)
            sage: h1 == h2
            True

            sage: h1 = hash(QQbar.zeta(6))
            sage: h2 = hash(QQbar(1/2 + I*sqrt(3)/2))
            sage: h1 == h2
            True

        Unfortunately, the hash code for algebraic numbers which are close
        enough to each other are the same. (This is inevitable, if
        equal algebraic reals give the same hash code and hashing does
        not always trigger exact computation.)::

            sage: h1 = hash(QQbar(0))
            sage: h2 = hash(QQbar(1/2^100))
            sage: hash(h1) == hash(h2)
            True

        """

        # The only way I can think of to hash algebraic numbers without
        # always triggering exact computation is to use interval_exact().
        # However, interval_exact() always triggers exact computation
        # if the number is exactly representable in floating point, which
        # is presumably not too unlikely (algebraic reals like 0, 1/2,
        # 1, or 2 are presumably not uncommon).

        # So I modify the algebraic real by adding 1/123456789 to it before
        # calling interval_exact(). Then, exact computation will be triggered
        # by algebraic reals which are sufficiently close to
        # (some floating point number minus 1/123456789). Hopefully,
        # -1/123456789 comes up in algebraic real computations far less
        # often than 0 does. Algebraic numbers have a similar offset added,
        # with an additional complex component of 1/987654321*I.

        # All of this effort to avoid exact computation is probably wasted,
        # anyway... in almost all uses of hash codes, if the hash codes
        # match, the next step is to compare for equality; and comparing
        # for equality often requires exact computation. (If a==b,
        # then checking a==b requires exact computation unless (a is b).)

        if self.parent() is AA:
            return hash((self + AA_hash_offset).interval_exact(RIF))
        else:
            return hash((self + QQbar_hash_offset).interval_exact(CIF))

    def __bool__(self):
        """
        Check whether ``self`` is nonzero.

        This is fast if interval arithmetic proves it and in many other cases.
        Though, it might be slow in very particular cases where the number is
        actually zero or very close to zero.

        EXAMPLES::

            sage: bool(QQbar.zeta(2) + 1)
            False
            sage: bool(QQbar.zeta(7) / (2^500))
            True

            sage: bool(QQbar(I).imag())
            True
            sage: bool(QQbar(I).real())
            False

        The following is very fast, even though the number is really small::

            sage: a1 = QQbar(2).sqrt() - 16616132878186749607/11749380235262596085
            sage: a2 = QQbar(2).sqrt() - 16616132878186749607/11749380235262596085
            sage: bool(a1 + a2)
            True
            sage: bool(a1 - a2)
            False

            sage: a = QQbar(2).sqrt() - 16616132878186749607/11749380235262596085
            sage: b = QQbar(2).sqrt() - 6882627592338442563/4866752642924153522
            sage: c = QQbar(3).sqrt() - 142437039878091970439/82236063316189858921
            sage: d = (59/2)**(1000/7)
            sage: e = (a + b + c) * (a + b - c) * (a - b) * (a - b - c) / d
            sage: bool(e)
            True
            sage: bool(e.abs() < 2**-500)
            True

        An identity between roots of unity::

            sage: z3 = QQbar.zeta(3)
            sage: z4 = QQbar.zeta(4)
            sage: z5 = QQbar.zeta(5)
            sage: p1 = (z3 + z4 + z5)**2
            sage: p2 = (z3 - z4 - z5)**2
            sage: p3 = (z3 - z4 + z5)**2
            sage: p4 = (z3 + z4 - z5)**2
            sage: bool(p1 - p2 + p3 - p4 - 8 * QQbar.zeta(15)**8)
            False

        Test some non-trivial zeros::

            sage: x = polygen(ZZ)
            sage: a = (AA(2).sqrt() + AA(3).sqrt() + AA(5).sqrt())^2
            sage: b = 10 + 2*max((x^4 - 62*x^2 - 240*x - 239).roots(AA, False))
            sage: bool(a - b)
            False

            sage: d = sum(AA(k)**(1/k) for k in [2..100])
            sage: bool(d * (a - b))
            False
            sage: bool((a - b) * d)
            False
            sage: bool(d * (a - b) * d)
            False
            sage: bool((a - b) / d)
            False

            sage: d = sum(QQbar(-k)**(1/k) for k in [2..100])
            sage: bool(d * (a - b))
            False
            sage: bool((a - b) * d)
            False
            sage: bool(d * (a - b) * d)
            False
            sage: bool((a - b) / d)
            False
        """
        # case 0: trivial tests
        if not self._value.contains_zero():
            return True
        elif self._value.is_zero():
            if not isinstance(self._descr, ANRational):
                self._set_descr(ANRational(QQ.zero()))
            return False

        # case 1: cheap tests
        sd = self._descr
        if isinstance(sd, ANExtensionElement):
            # The ANExtensionElement returns an ANRational
            # instead, if the number is zero.
            return True
        elif isinstance(sd, ANRational):
            return bool(sd._value)
        elif isinstance(sd, ANUnaryExpr) and sd._op != 'real' and sd._op != 'imag':
            ans = bool(sd._arg)
            if not ans:
                self._set_descr(ANRational(QQ.zero()))
            return ans
        elif isinstance(sd, ANBinaryExpr) and sd._op is operator.mul:
            ans = bool(sd._left) and bool(sd._right)
            if not ans:
                self._set_descr(ANRational(QQ.zero()))
            return ans
        elif isinstance(sd, ANBinaryExpr) and sd._op is operator.truediv:
            ans = bool(sd._left)
            if not ans:
                self._set_descr(ANRational(QQ.zero()))
            return ans

        # case 2: try more precision
        if self._value.prec() < 128:
            self._more_precision()
            if not self._value.contains_zero():
                return True

        # case 3: try with minpoly in case of x+y or x-y
        if isinstance(sd, ANBinaryExpr):
            op = sd._op
            left = sd._left
            right = sd._right if op is operator.sub else -sd._right

            lp = left.minpoly()
            rp = right.minpoly()
            if lp != rp:
                return True

            c = cmp_elements_with_same_minpoly(left, right, left.minpoly())
            if c is not None:
                if c == 0:
                    self._set_descr(ANRational(QQ.zero()))
                return bool(c)

        # Sigh...
        self.exactify()
        return bool(self)

    def is_square(self):
        """
        Return whether or not this number is square.

        OUTPUT:

        (boolean)
        ``True`` in all cases for elements of ``QQbar``;
        ``True`` for non-negative elements of ``AA``;
        otherwise ``False``

        EXAMPLES::

            sage: AA(2).is_square()
            True
            sage: AA(-2).is_square()
            False
            sage: QQbar(-2).is_square()
            True
            sage: QQbar(I).is_square()
            True
        """
        if self.parent() is AA:
            return bool(self >= 0)
        else:
            return True

    def is_integer(self):
        """
        Return ``True`` if this number is a integer.

        EXAMPLES::

            sage: QQbar(2).is_integer()
            True
            sage: QQbar(1/2).is_integer()
            False
        """
        return self in ZZ

    def sqrt(self, all=False, extend=True):
        """
        Return the square root(s) of this number.

        INPUT:

        - ``extend`` - bool (default: True); ignored if self is in QQbar, or
          positive in AA. If self is negative in AA, do the following: if True,
          return a square root of self in QQbar, otherwise raise a ValueError.

        - ``all`` - bool (default: False); if True, return a list of all square
          roots. If False, return just one square root, or raise an ValueError
          if self is a negative element of AA and extend=False.

        OUTPUT:

        Either the principal square root of self, or a list of its
        square roots (with the principal one first).

        EXAMPLES::

            sage: AA(2).sqrt()
            1.414213562373095?

            sage: QQbar(I).sqrt()
            0.7071067811865475? + 0.7071067811865475?*I
            sage: QQbar(I).sqrt(all=True)
            [0.7071067811865475? + 0.7071067811865475?*I, -0.7071067811865475? - 0.7071067811865475?*I]

            sage: a = QQbar(0)
            sage: a.sqrt()
            0
            sage: a.sqrt(all=True)
            [0]

            sage: a = AA(0)
            sage: a.sqrt()
            0
            sage: a.sqrt(all=True)
            [0]

        This second example just shows that the program does not care where 0
        is defined, it gives the same answer regardless. After all, how many
        ways can you square-root zero?

        ::

            sage: AA(-2).sqrt()
            1.414213562373095?*I

            sage: AA(-2).sqrt(all=True)
            [1.414213562373095?*I, -1.414213562373095?*I]

            sage: AA(-2).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: -2 is not a square in AA, being negative. Use extend = True for a square root in QQbar.
        """
        # deal with 0 first:

        if self.is_zero():
            if all:
                return [self]
            else:
                return self

        # raise an error if appropriate:

        if self.parent() is AA and self < 0 and not extend:
            if not all:
                raise ValueError(lazy_string("%s is not a square in AA, being negative. Use extend = True for a square root in QQbar.", self))
            else:
                return []

        root = self ** ~ZZ(2)

        if all:
            return [root, -root]
        else:
            return root

    def nth_root(self, n, all=False):
        r"""
        Return the ``n``-th root of this number.

        INPUT:

        -  ``all`` - bool (default: ``False``). If ``True``, return a list of
           all `n`-th roots as complex algebraic numbers.

        .. WARNING::

            Note that for odd `n`, all=`False` and negative real numbers,
            ``AlgebraicReal`` and ``AlgebraicNumber`` values give different
            answers: ``AlgebraicReal`` values prefer real results, and
            ``AlgebraicNumber`` values return the principal root.

        EXAMPLES::

            sage: AA(-8).nth_root(3)
            -2
            sage: QQbar(-8).nth_root(3)
            1.000000000000000? + 1.732050807568878?*I
            sage: QQbar.zeta(12).nth_root(15)
            0.9993908270190957? + 0.03489949670250097?*I

        You can get all ``n``-th roots of algebraic numbers::

            sage: AA(-8).nth_root(3, all=True)
            [1.000000000000000? + 1.732050807568878?*I,
            -2.000000000000000? + 0.?e-18*I,
            1.000000000000000? - 1.732050807568878?*I]

            sage: QQbar(1+I).nth_root(4, all=True)
            [1.069553932363986? + 0.2127475047267431?*I,
             -0.2127475047267431? + 1.069553932363986?*I,
             -1.069553932363986? - 0.2127475047267431?*I,
             0.2127475047267431? - 1.069553932363986?*I]

        TESTS::

            sage: AA(-8).nth_root(3, all=True)[1]
            -2.000000000000000? + 0.?e-18*I
            sage: _.parent()
            Algebraic Field

            sage: AA(-2).nth_root(5, all=True) == QQbar(-2).nth_root(5, all=True)   # long time
            True
        """
        if not all:
            return self ** ~ZZ(n)
        else:
            root = QQbar(self) ** ~ZZ(n)
            zlist = [root]
            zeta = QQbar.zeta(n)
            for k in range(1, n):
                root *= zeta
                zlist.append(root)
            return zlist

    def as_number_field_element(self, minimal=False, embedded=False, prec=53):
        r"""
        Return a number field containing this value, a representation of
        this value as an element of that number field, and a homomorphism
        from the number field back to ``AA`` or ``QQbar``.

        INPUT:

        - ``minimal`` -- Boolean (default: ``False``). Whether to minimize the
          degree of the extension.

        - ``embedded`` -- Boolean (default: ``False``). Whether to make the
          NumberField embedded.

        - ``prec`` -- integer (default: ``53``). The number of bit of precision
          to guarantee finding real roots.

        This may not return the smallest such number field, unless
        ``minimal=True`` is specified.

        To compute a single number field containing multiple algebraic
        numbers, use the function
        ``number_field_elements_from_algebraics`` instead.

        EXAMPLES::

            sage: QQbar(sqrt(8)).as_number_field_element()
            (Number Field in a with defining polynomial y^2 - 2, 2*a, Ring morphism:
                From: Number Field in a with defining polynomial y^2 - 2
                To:   Algebraic Real Field
                Defn: a |--> 1.414213562373095?)

            sage: x = polygen(ZZ)
            sage: p = x^3 + x^2 + x + 17
            sage: (rt,) = p.roots(ring=AA, multiplicities=False); rt
            -2.804642726932742?

            sage: (nf, elt, hom) = rt.as_number_field_element()
            sage: nf, elt, hom
            (Number Field in a with defining polynomial y^3 - 2*y^2 - 31*y - 50,
             a^2 - 5*a - 19,
             Ring morphism:
               From: Number Field in a with defining polynomial y^3 - 2*y^2 - 31*y - 50
               To:   Algebraic Real Field
               Defn: a |--> 7.237653139801104?)
            sage: elt == rt
            False
            sage: AA(elt)
            Traceback (most recent call last):
            ...
            ValueError: need a real or complex embedding to convert a non rational element of a number field into an algebraic number
            sage: hom(elt) == rt
            True

        Creating an element of an embedded number field::

            sage: (nf, elt, hom) = rt.as_number_field_element(embedded=True)
            sage: nf.coerce_embedding()
            Generic morphism:
              From: Number Field in a with defining polynomial y^3 - 2*y^2 - 31*y - 50 with a = 7.237653139801104?
              To:   Algebraic Real Field
              Defn: a -> 7.237653139801104?
            sage: elt
            a^2 - 5*a - 19
            sage: elt.parent() == nf
            True
            sage: hom(elt).parent()
            Algebraic Real Field
            sage: hom(elt) == rt
            True
            sage: elt == rt
            True
            sage: AA(elt)
            -2.804642726932742?
            sage: RR(elt)
            -2.80464272693274

        A complex algebraic number as an element of an embedded number field::

            sage: num = QQbar(sqrt(2) + 3^(1/3)*I)
            sage: nf, elt, hom = num.as_number_field_element(embedded=True)
            sage: hom(elt).parent() is QQbar
            True
            sage: nf.coerce_embedding() is not None
            True
            sage: QQbar(elt) == num == hom(elt)
            True

        We see an example where we do not get the minimal number field unless
        we specify ``minimal=True``::

            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt3b = rt2 + rt3 - rt2
            sage: rt3b.as_number_field_element()
            (Number Field in a with defining polynomial y^4 - 4*y^2 + 1, a^2 - 2, Ring morphism:
                From: Number Field in a with defining polynomial y^4 - 4*y^2 + 1
                To:   Algebraic Real Field
                Defn: a |--> -1.931851652578137?)
            sage: rt3b.as_number_field_element(minimal=True)
            (Number Field in a with defining polynomial y^2 - 3, a, Ring morphism:
               From: Number Field in a with defining polynomial y^2 - 3
               To:   Algebraic Real Field
               Defn: a |--> 1.732050807568878?)
        """
        return number_field_elements_from_algebraics(self, minimal=minimal, embedded=embedded, prec=prec)

    def exactify(self):
        """
        Compute an exact representation for this number.

        EXAMPLES::

            sage: two = QQbar(4).nth_root(4)^2
            sage: two
            2.000000000000000?
            sage: two.exactify()
            sage: two
            2
        """
        od = self._descr
        if isinstance(od, (ANRational, ANExtensionElement)):
            return
        self._set_descr(self._descr.exactify())

    def _set_descr(self, new_descr):
        """
        Set ``self._descr`` to ``new_descr``, and update
        ``self._value`` accordingly.

        EXAMPLES::

            sage: c = QQbar(-1)**(1/3) - QQbar(3)**(1/2)/2*QQbar.gen()
            sage: c._value
            0.5000000000000000000? + 0.?e-19*I
            sage: c.exactify()   # indirect doctest
            sage: c._value
            0.500000000000000000000?
        """
        self._descr = new_descr
        new_val = self._descr._interval_fast(self.parent().default_interval_prec())
        if is_RealIntervalFieldElement(new_val) and is_ComplexIntervalFieldElement(self._value):
            self._value = self._value.real().intersection(new_val)
        elif is_RealIntervalFieldElement(self._value) and is_ComplexIntervalFieldElement(new_val):
            self._value = self._value.intersection(new_val.real())
        else:
            self._value = self._value.intersection(new_val)

    def simplify(self):
        """
        Compute an exact representation for this number, in the
        smallest possible number field.

        EXAMPLES::

            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt2b = rt3 + rt2 - rt3
            sage: rt2b.exactify()
            sage: rt2b._exact_value()
            a^3 - 3*a where a^4 - 4*a^2 + 1 = 0 and a in -0.5176380902050415?
            sage: rt2b.simplify()
            sage: rt2b._exact_value()
            a where a^2 - 2 = 0 and a in 1.414213562373095?
        """
        self.exactify()
        od = self._descr
        if od.is_simple():
            return
        self._set_descr(od.simplify(self))

    def _exact_field(self):
        """
        Return a generator for a number field that includes this number
        (not necessarily the smallest such number field).

        EXAMPLES::

            sage: QQbar(2)._exact_field()
            Trivial generator
            sage: (sqrt(QQbar(2)) + sqrt(QQbar(19)))._exact_field()
            Number Field in a with defining polynomial y^4 - 20*y^2 + 81 with a in -3.789313782671036?
            sage: (QQbar(7)^(3/5))._exact_field()
            Number Field in a with defining polynomial y^5 - 2*y^4 - 18*y^3 + 38*y^2 + 82*y - 181 with a in 2.554256611698490?
        """
        sd = self._descr
        if isinstance(sd, (ANRational, ANExtensionElement)):
            return sd.generator()
        self.exactify()
        return self._exact_field()

    def _exact_value(self):
        r"""
        Return an ``ANRational`` or an ``ANExtensionElement`` representing this
        value.

        EXAMPLES::

            sage: QQbar(2)._exact_value()
            2
            sage: (sqrt(QQbar(2)) + sqrt(QQbar(19)))._exact_value()
            -1/9*a^3 + a^2 + 11/9*a - 10 where a^4 - 20*a^2 + 81 = 0 and a in -3.789313782671036?
            sage: (QQbar(7)^(3/5))._exact_value()
            2*a^4 + 2*a^3 - 34*a^2 - 17*a + 150 where a^5 - 2*a^4 - 18*a^3 + 38*a^2 + 82*a - 181 = 0 and a in 2.554256611698490?
        """
        sd = self._descr
        if isinstance(sd, (ANRational, ANExtensionElement)):
            return sd
        self.exactify()
        return self._descr

    def _more_precision(self):
        """
        Recompute the interval bounding this number with higher-precision
        interval arithmetic.

        EXAMPLES::

            sage: rt2 = sqrt(QQbar(2))
            sage: rt2._value
            1.4142135623730950488?
            sage: rt2._more_precision()
            sage: rt2._value
            1.41421356237309504880168872420969807857?
            sage: rt2._more_precision()
            sage: rt2._value
            1.41421356237309504880168872420969807856967187537694807317667973799073247846211?
        """
        prec = self._value.prec()
        self._value = self._descr._interval_fast(prec * 2)

    def minpoly(self):
        """
        Compute the minimal polynomial of this algebraic number.
        The minimal polynomial is the monic polynomial of least degree
        having this number as a root; it is unique.

        EXAMPLES::

            sage: QQbar(4).sqrt().minpoly()
            x - 2
            sage: ((QQbar(2).nth_root(4))^2).minpoly()
            x^2 - 2
            sage: v = sqrt(QQbar(2)) + sqrt(QQbar(3)); v
            3.146264369941973?
            sage: p = v.minpoly(); p
            x^4 - 10*x^2 + 1
            sage: p(RR(v.real()))
            1.31006316905768e-14
        """
        try:
            return self._minimal_polynomial
        except AttributeError:
            self.exactify()
            self._minimal_polynomial = self._descr.minpoly()
            return self._minimal_polynomial

    def degree(self):
        """
        Return the degree of this algebraic number (the degree of its
        minimal polynomial, or equivalently, the degree of the smallest
        algebraic extension of the rationals containing this number).

        EXAMPLES::

            sage: QQbar(5/3).degree()
            1
            sage: sqrt(QQbar(2)).degree()
            2
            sage: QQbar(17).nth_root(5).degree()
            5
            sage: sqrt(3+sqrt(QQbar(8))).degree()
            2
        """
        return self.minpoly().degree()

    def interval_fast(self, field):
        r"""
        Given a :class:`RealIntervalField` or
        :class:`ComplexIntervalField`, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field. (More precision may be used
        in the computation.)  The returned interval may be arbitrarily
        imprecise, if this number is the result of a sufficiently long
        computation chain.

        EXAMPLES::

            sage: x = AA(2).sqrt()
            sage: x.interval_fast(RIF)
            1.414213562373095?
            sage: x.interval_fast(RealIntervalField(200))
            1.414213562373095048801688724209698078569671875376948073176680?
            sage: x = QQbar(I).sqrt()
            sage: x.interval_fast(CIF)
            0.7071067811865475? + 0.7071067811865475?*I
            sage: x.interval_fast(RIF)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert complex interval 0.7071067811865475244? + 0.7071067811865475244?*I to real interval
        """
        while self._value.prec() < field.prec():
            self._more_precision()
        return field(self._value)

    def interval_diameter(self, diam):
        """
        Compute an interval representation of self with ``diameter()`` at
        most ``diam``. The precision of the returned value is unpredictable.

        EXAMPLES::

            sage: AA(2).sqrt().interval_diameter(1e-10)
            1.4142135623730950488?
            sage: AA(2).sqrt().interval_diameter(1e-30)
            1.41421356237309504880168872420969807857?
            sage: QQbar(2).sqrt().interval_diameter(1e-10)
            1.4142135623730950488?
            sage: QQbar(2).sqrt().interval_diameter(1e-30)
            1.41421356237309504880168872420969807857?
        """
        if diam <= 0:
            raise ValueError('diameter must be positive in interval_diameter')

        while self._value.diameter() > diam:
            self._more_precision()

        return self._value

    def interval(self, field):
        r"""
        Given an interval (or ball) field (real or complex, as appropriate) of
        precision `p`, compute an interval representation of self with
        ``diameter()`` at most `2^{-p}`; then round that representation into
        the given field. Here ``diameter()`` is relative diameter for
        intervals not containing 0, and absolute diameter for
        intervals that do contain 0; thus, if the returned interval
        does not contain 0, it has at least `p-1` good bits.

        EXAMPLES::

            sage: RIF64 = RealIntervalField(64)
            sage: x = AA(2).sqrt()
            sage: y = x*x
            sage: y = 1000 * y - 999 * y
            sage: y.interval_fast(RIF64)
            2.000000000000000?
            sage: y.interval(RIF64)
            2.000000000000000000?
            sage: CIF64 = ComplexIntervalField(64)
            sage: x = QQbar.zeta(11)
            sage: x.interval_fast(CIF64)
            0.8412535328311811689? + 0.5406408174555975821?*I
            sage: x.interval(CIF64)
            0.8412535328311811689? + 0.5406408174555975822?*I
            sage: x.interval(CBF) # abs tol 1e-16
            [0.8412535328311812 +/- 3.12e-17] + [0.5406408174555976 +/- 1.79e-17]*I

        The following implicitly use this method::

            sage: RIF(AA(5).sqrt())
            2.236067977499790?
            sage: AA(-5).sqrt().interval(RIF)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2.236067977499790?*I to real interval

        TESTS:

        Check that :trac:`20209` is fixed::

            sage: RIF(QQbar(2).sqrt())
            1.414213562373095?
            sage: RIF(QQbar.gen() + QQbar(2).sqrt() - QQbar.gen())
            1.414213562373095?
            sage: RIF((QQbar.gen() + QQbar(2).sqrt() - QQbar.gen()).sqrt())
            1.189207115002722?

            sage: RealIntervalField(129)(QQbar(3).sqrt())
            1.73205080756887729352744634150587236695?
            sage: RIF(QQbar.gen())
            Traceback (most recent call last):
            ...
            TypeError: unable to convert I to real interval
        """
        target = RR(1.0) >> field.prec()
        val = self.interval_diameter(target)
        if (isinstance(field, (RealIntervalField_class, RealBallField))
                and is_ComplexIntervalFieldElement(val)):
            if val.imag().is_zero():
                return field(val.real())
            elif self.imag().is_zero():
                return field(self.real())
            else:
                raise TypeError(lazy_string("unable to convert %s to real interval", self))
        else:
            return field(val)

    _arb_ = _acb_ = _complex_mpfi_ = _real_mpfi_ = interval

    def radical_expression(self):
        r"""
        Attempt to obtain a symbolic expression using radicals. If no
        exact symbolic expression can be found, the algebraic number
        will be returned without modification.

        EXAMPLES::

            sage: AA(1/sqrt(5)).radical_expression()                                    # optional - sage.symbolic
            sqrt(1/5)
            sage: AA(sqrt(5 + sqrt(5))).radical_expression()                            # optional - sage.symbolic
            sqrt(sqrt(5) + 5)
            sage: QQbar.zeta(5).radical_expression()                                    # optional - sage.symbolic
            1/4*sqrt(5) + 1/2*sqrt(-1/2*sqrt(5) - 5/2) - 1/4
            sage: x = polygen(QQ, 'x')
            sage: a = (x^7 - x - 1).roots(AA, False)[0]
            sage: a.radical_expression()                                                # optional - sage.symbolic
            1.112775684278706?
            sage: a.radical_expression().parent() == SR                                 # optional - sage.symbolic
            False
            sage: a = sorted((x^7-x-1).roots(QQbar, False), key=imag)[0]
            sage: a.radical_expression()                                                # optional - sage.symbolic
            -0.3636235193291805? - 0.9525611952610331?*I
            sage: QQbar.zeta(5).imag().radical_expression()                             # optional - sage.symbolic
            1/2*sqrt(1/2*sqrt(5) + 5/2)
            sage: AA(5/3).radical_expression()                                          # optional - sage.symbolic
            5/3
            sage: AA(5/3).radical_expression().parent() == SR                           # optional - sage.symbolic
            True
            sage: QQbar(0).radical_expression()                                         # optional - sage.symbolic
            0

        TESTS:

        In this example we find the correct answer despite the fact that
        multiple roots overlap with the current value. As a consequence,
        the precision of the evaluation will have to be increased.

        ::

            sage: a = AA(sqrt(2) + 10^25)                                               # optional - sage.symbolic
            sage: p = a.minpoly()                                                       # optional - sage.symbolic
            sage: v = a._value                                                          # optional - sage.symbolic
            sage: f = ComplexIntervalField(v.prec())                                    # optional - sage.symbolic
            sage: var('x')                                                              # optional - sage.symbolic
            x
            sage: [f(b.rhs()).overlaps(f(v)) for b in SR(p).solve(x)]                   # optional - sage.symbolic
            [True, True]
            sage: a.radical_expression()                                                # optional - sage.symbolic
            sqrt(2) + 10000000000000000000000000
        """
        from sage.symbolic.ring import SR  # Lazy to avoid cyclic dependency

        # Adapted from NumberFieldElement._symbolic_()
        poly = self.minpoly()
        if is_ComplexIntervalFieldElement(self._value):
            interval_field = self._value.parent()
        else:
            interval_field = ComplexIntervalField(self._value.prec())
        roots = poly.roots(SR, multiplicities=False)
        if len(roots) != poly.degree():
            return self
        while True:
            candidates = []
            for root in roots:
                if interval_field(root).overlaps(interval_field(self._value)):
                    candidates.append(root)
            if len(candidates) == 1:
                return candidates[0]
            roots = candidates
            interval_field = interval_field.to_prec(interval_field.prec() * 2)

    def _maxima_init_(self, I=None):
        r"""
        EXAMPLES::

            sage: maxima(AA(7))
            7
            sage: maxima(QQbar(sqrt(5/2)))                                              # optional - sage.symbolic
            sqrt(10)/2
            sage: maxima(AA(-sqrt(5)))                                                  # optional - sage.symbolic
            -sqrt(5)
            sage: maxima(QQbar(sqrt(-2)))                                               # optional - sage.symbolic
            sqrt(2)*%i
            sage: maxima(AA(2+sqrt(5)))                                                 # optional - sage.symbolic
            sqrt(5)+2
            sage: maxima(QQ[x](x^7 - x - 1).roots(AA, False)[0])                        # optional - sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot find radical expression
        """
        try:
            return self._rational_()._maxima_init_()
        except ValueError:
            pass
        rad = self.radical_expression()
        if isinstance(rad.parent(), sage.rings.abc.SymbolicRing):
            return rad._maxima_init_()
        raise NotImplementedError('cannot find radical expression')


class AlgebraicNumber(AlgebraicNumber_base):
    r"""
    The class for algebraic numbers (complex numbers which are the roots
    of a polynomial with integer coefficients). Much of its functionality
    is inherited from :class:`AlgebraicNumber_base`.

    .. automethod:: _richcmp_
    """
    def __init__(self, x):
        r"""
        Initialize this AlgebraicNumber object.

        EXAMPLES::

            sage: t = QQbar.zeta(5)
            sage: type(t)
            <class 'sage.rings.qqbar.AlgebraicNumber'>
        """
        AlgebraicNumber_base.__init__(self, QQbar, x)

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = QQbar.zeta(5)
            sage: loads(dumps(t)) == t
            True
        """
        return (AlgebraicNumber, (self._descr, ))

    def _richcmp_(self, other, op):
        r"""
        Compare two algebraic numbers, lexicographically. (That is,
        first compare the real components; if the real components are
        equal, compare the imaginary components.)

        EXAMPLES::

            sage: x = QQbar.zeta(3); x
            -0.500000000000000? + 0.866025403784439?*I
            sage: QQbar(-1) < x
            True
            sage: QQbar(-1/2) < x
            True
            sage: QQbar(0) > x
            True

        One problem with this lexicographic ordering is the fact that if
        two algebraic numbers have the same real component, that real
        component has to be compared for exact equality, which can be
        a costly operation.  For the special case where both numbers
        have the same minimal polynomial, that cost can be avoided,
        though (see :trac:`16964`)::

            sage: x = polygen(ZZ)
            sage: p = 69721504*x^8 + 251777664*x^6 + 329532012*x^4 + 184429548*x^2 + 37344321
            sage: sorted(p.roots(QQbar,False))
            [-0.0221204634374361? - 1.090991904211621?*I,
             -0.0221204634374361? + 1.090991904211621?*I,
             -0.8088604911480535?*I,
             0.?e-215 - 0.7598602580415435?*I,
             0.?e-229 + 0.7598602580415435?*I,
             0.8088604911480535?*I,
             0.0221204634374361? - 1.090991904211621?*I,
             0.0221204634374361? + 1.090991904211621?*I]

        It also works for comparison of conjugate roots even in a degenerate
        situation where many roots have the same real part. In the following
        example, the polynomial ``p2`` is irreducible and all its roots have
        real part equal to `1`::

            sage: p1 = x^8 + 74*x^7 + 2300*x^6 + 38928*x^5 + \
            ....: 388193*x^4 + 2295312*x^3 + 7613898*x^2 + \
            ....: 12066806*x + 5477001
            sage: p2 = p1((x-1)^2)
            sage: sum(1 for r in p2.roots(CC,False) if abs(r.real() - 1) < 0.0001)
            16
            sage: r1 = QQbar.polynomial_root(p2, CIF(1, (-4.1,-4.0)))
            sage: r2 = QQbar.polynomial_root(p2, CIF(1, (4.0, 4.1)))
            sage: all([r1<r2, r1==r1, r2==r2, r2>r1])
            True

        Though, comparing roots which are not equal or conjugate is much
        slower because the algorithm needs to check the equality of the real
        parts::

            sage: sorted(p2.roots(QQbar,False))   # long time - 3 secs
            [1.000000000000000? - 4.016778562562223?*I,
             1.000000000000000? - 3.850538755978243?*I,
             1.000000000000000? - 3.390564396412898?*I,
             ...
             1.000000000000000? + 3.390564396412898?*I,
             1.000000000000000? + 3.850538755978243?*I,
             1.000000000000000? + 4.016778562562223?*I]

        TESTS::

            sage: QQbar.zeta(6) == QQbar(1/2 + I*sqrt(3)/2)                             # optional - sage.symbolic
            True
            sage: QQbar(I) == QQbar(I * (2^100+1)/(2^100))
            False
            sage: QQbar(2) == 2
            True
            sage: QQbar(2) == GF(7)(2)                                                  # optional - sage.rings.finite_rings
            False
            sage: GF(7)(2) in QQbar                                                     # optional - sage.rings.finite_rings
            False

            sage: QQbar.zeta(6) != QQbar(1/2 + I*sqrt(3)/2)                             # optional - sage.symbolic
            False
            sage: QQbar(I) != QQbar(I * (2^100+1)/(2^100))
            True
            sage: QQbar(2) != 2
            False
            sage: QQbar(2) != GF(7)(2)                                                  # optional - sage.rings.finite_rings
            True

            sage: QQbar.zeta(3).real() == -1/2
            True

        Check that :trac:`26593` is fixed (the test here has to be repeated
        twice)::

            sage: pi = x^7 - 2*x^6 + x^3 - 2*x^2 + 2*x - 1
            sage: b = pi.roots(ring=QQbar)[3][0]
            sage: pi = b.minpoly()
            sage: K = NumberField(pi, 'b', embedding=b)
            sage: pi = x^7 - 2*x^6 + x^3 - 2*x^2 + 2*x - 1
            sage: b = pi.roots(ring=QQbar)[3][0]
            sage: pi = b.minpoly()
            sage: K = NumberField(pi, 'b', embedding=b)

        Check that :trac:`29220` is fixed::

            sage: a = AA(2**(1/2) - 2**(1/3))
            sage: b = 808620184/5240825825
            sage: a < b
            True
            sage: a < b
            True

            sage: a = AA(2^(1/3))
            sage: r = 3085094589/2448641198
            sage: a < r
            False
            sage: a > r
            True
            sage: a < r
            False
            sage: a > r
            True
        """
        if self is other:
            return rich_to_bool(op, 0)

        # case 0: rationals
        sd = self._descr
        od = other._descr
        if isinstance(sd, ANRational) and isinstance(od, ANRational):
            return richcmp(sd._value, od._value, op)

        # case 1: real parts are clearly distinct
        ri1 = self._value.real()
        ri2 = other._value.real()
        if not ri1.overlaps(ri2):
            # NOTE: do not call richcmp here as self._value and other._value
            # might have different precisions. See
            # https://github.com/sagemath/sage/issues/29220
            return ri1._richcmp_(ri2, op)

        if op == op_EQ or op == op_NE:
            # some cheap and quite common tests where we can decide
            # equality or difference
            if not self._value.imag().overlaps(other._value.imag()):
                return op == op_NE
            if isinstance(sd, ANRational) and not sd._value:
                return bool(other) == (op == op_NE)
            elif isinstance(od, ANRational) and not od._value:
                return bool(self) == (op == op_NE)
            elif (isinstance(sd, ANExtensionElement) and
                  isinstance(od, ANExtensionElement) and
                  sd._generator is od._generator):
                return sd._value == od._value if op == op_EQ else sd._value != od._value

        # case 2: possibly equal or conjugate values
        # (this case happen a lot when sorting the roots of a real polynomial)
        ci1 = self._value.imag().abs()
        ci2 = other._value.imag().abs()
        if ci1.overlaps(ci2) and self.minpoly() == other.minpoly():
            c = cmp_elements_with_same_minpoly(self, other, self.minpoly())
            if c is not None:
                return rich_to_bool(op, c)

        # case 3: try hard to compare real parts and imaginary parts
        srp = self.real()
        orp = other.real()
        if srp != orp:
            return richcmp_not_equal(srp, orp, op)
        return richcmp(self.imag(), other.imag(), op)

    def _mpfr_(self, field):
        r"""
        Given a ``RealField``, compute a good approximation to self in
        that field. Works only if the imaginary component of self is
        exactly zero; otherwise it raises a ``ValueError``.

        EXAMPLES::

            sage: QQbar(sqrt(2))._mpfr_(RR)
            1.41421356237309
            sage: QQbar(-22/7)._mpfr_(RR)
            -3.14285714285714
            sage: QQbar.zeta(3)._mpfr_(RR)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
        """
        return AA(self)._mpfr_(field)

    def __float__(self):
        r"""
        Compute a good float approximation to self. Works only if the
        imaginary component of self is exactly zero; otherwise it
        raises a ``ValueError``.

        EXAMPLES::

            sage: QQbar(sqrt(2)).__float__()
            1.414213562373095
            sage: float(QQbar(-22/7))
            -3.1428571428571432
            sage: float(QQbar.zeta(3))
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
        """
        return AA(self).__float__()

    def __complex__(self):
        r"""
        Compute a good complex approximation to self.

        EXAMPLES::

            sage: QQbar(sqrt(2)).__complex__()
            (1.414213562373095+0j)
            sage: complex(QQbar.zeta(3))
            (-0.5+0.8660254037844386j)
        """
        return CC(self).__complex__()

    def _complex_double_(self, cdf):
        r"""
        Compute a good approximation to self in CDF.

        EXAMPLES::

            sage: QQbar(sqrt(-5))._complex_double_(CDF)
            2.23606797749979*I
            sage: CDF(QQbar.zeta(12))
            0.8660254037844386 + 0.5*I
        """
        return cdf(CC(self))

    def _interval_fast(self, prec):
        r"""
        Shortcut for :meth:`AlgebraicNumber_base.interval_fast` which uses the complex interval field.

        EXAMPLES::

            sage: QQbar(sqrt(-5))._interval_fast(100)
            2.236067977499789696409173...?*I
        """
        return self.interval_fast(ComplexIntervalField(prec))

    def _integer_(self, ZZ=None):
        """
        Return self as an Integer.

        EXAMPLES::

            sage: QQbar(0)._integer_()
            0
            sage: QQbar(0)._integer_().parent()
            Integer Ring
            sage: QQbar.zeta(6)._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
            sage: QQbar(sqrt(17))._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce non-integral Algebraic Real 4.123105625617660? to Integer
            sage: QQbar(sqrt(16))._integer_()
            4
            sage: v = QQbar(1 + I*sqrt(3))^5 + QQbar(16*sqrt(3)*I); v
            16.00000000000000? + 0.?e-17*I
            sage: v._integer_()
            16
        """
        return AA(self)._integer_(ZZ)

    def _rational_(self):
        """
        Return self as a Rational.

        EXAMPLES::

            sage: QQbar(-22/7)._rational_()
            -22/7
            sage: QQbar(3)._rational_().parent()
            Rational Field
            sage: (QQbar.zeta(7)^3)._rational_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
            sage: QQbar(sqrt(2))._rational_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce irrational Algebraic Real 1.414213562373095? to Rational
            sage: v1 = QQbar(1/3 + I*sqrt(5))^7
            sage: v2 = QQbar((100336/729*golden_ratio - 50168/729)*I)
            sage: v = v1 + v2; v
            -259.6909007773206? + 0.?e-15*I
            sage: v._rational_()
            -567944/2187
        """
        return AA(self)._rational_()

    def real(self):
        r"""
        Return the real part of ``self``.

        EXAMPLES::

            sage: QQbar.zeta(5).real()
            0.3090169943749474?
        """
        return AlgebraicReal(self._descr.real(self))

    def imag(self):
        r"""
        Return the imaginary part of ``self``.

        EXAMPLES::

            sage: QQbar.zeta(7).imag()
            0.7818314824680299?
        """
        return AlgebraicReal(self._descr.imag(self))

    def conjugate(self):
        """
        Return the complex conjugate of ``self``.

        EXAMPLES::

            sage: QQbar(3 + 4*I).conjugate()
            3 - 4*I
            sage: QQbar.zeta(7).conjugate()
            0.6234898018587335? - 0.7818314824680299?*I
            sage: QQbar.zeta(7) + QQbar.zeta(7).conjugate()
            1.246979603717467? + 0.?e-18*I
        """
        return AlgebraicNumber(self._descr.conjugate(self))

    def norm(self):
        r"""
        Return ``self * self.conjugate()``.

        This is the algebraic definition of norm, if we view ``QQbar``
        as ``AA[I]``.

        EXAMPLES::

            sage: QQbar(3 + 4*I).norm()
            25
            sage: type(QQbar(I).norm())
            <class 'sage.rings.qqbar.AlgebraicReal'>
            sage: QQbar.zeta(1007).norm()
            1.000000000000000?
        """
        return AlgebraicReal(self._descr.norm(self))

    def interval_exact(self, field):
        r"""
        Given a :class:`ComplexIntervalField`, compute the best possible
        approximation of this number in that field. Note that if
        either the real or imaginary parts of this number are
        sufficiently close to some floating-point number (and, in
        particular, if either is exactly representable in floating-point),
        then this will trigger exact computation, which may be very slow.

        EXAMPLES::

            sage: a = QQbar(I).sqrt(); a
            0.7071067811865475? + 0.7071067811865475?*I
            sage: a.interval_exact(CIF)
            0.7071067811865475? + 0.7071067811865475?*I
            sage: b = QQbar((1+I)*sqrt(2)/2)                                            # optional - sage.symbolic
            sage: (a - b).interval(CIF)                                                 # optional - sage.symbolic
            0.?e-19 + 0.?e-18*I
            sage: (a - b).interval_exact(CIF)                                           # optional - sage.symbolic
            0
        """
        if not isinstance(field, sage.rings.abc.ComplexIntervalField):
            raise ValueError("AlgebraicNumber interval_exact requires a ComplexIntervalField")
        rfld = field._real_field()
        re = self.real().interval_exact(rfld)
        im = self.imag().interval_exact(rfld)
        return field(re, im)

    def _complex_mpfr_field_(self, field):
        r"""
        Compute an approximation to self in the given field, which must
        be a complex field.

        EXAMPLES::

            sage: a = QQbar(1 + I).sqrt()
            sage: t = a._complex_mpfr_field_(ComplexField(100)); t
            1.0986841134678099660398011952 + 0.45508986056222734130435775782*I
            sage: parent(t)
            Complex Field with 100 bits of precision
        """
        return self.complex_number(field)

    def complex_number(self, field):
        r"""
        Given the complex field ``field``, compute an accurate approximation of
        this element in that field.

        The approximation will be off by at most two ulp's in each component,
        except for components which are very close to zero, which will have an
        absolute error at most `2^{-prec+1}` where ``prec`` is the precision of
        the field.

        EXAMPLES::

            sage: a = QQbar.zeta(5)
            sage: a.complex_number(CC)
            0.309016994374947 + 0.951056516295154*I

            sage: b = QQbar(2).sqrt() + QQbar(3).sqrt() * QQbar.gen()
            sage: b.complex_number(ComplexField(128))
            1.4142135623730950488016887242096980786 + 1.7320508075688772935274463415058723669*I
        """
        v = self.interval(ComplexIntervalField(field.prec()))
        return field(v)

    def complex_exact(self, field):
        r"""
        Given a :class:`ComplexField`, return the best possible approximation of
        this number in that field. Note that if either component is
        sufficiently close to the halfway point between two floating-point
        numbers in the corresponding :class:`RealField`, then this will trigger
        exact computation, which may be very slow.

        EXAMPLES::

            sage: a = QQbar.zeta(9) + QQbar(I) + QQbar.zeta(9).conjugate(); a
            1.532088886237957? + 1.000000000000000?*I
            sage: a.complex_exact(CIF)
            1.532088886237957? + 1*I
        """
        rfld = field._real_field()
        re = self.real().real_exact(rfld)
        im = self.imag().real_exact(rfld)
        return field(re, im)

    def multiplicative_order(self):
        r"""
        Compute the multiplicative order of this algebraic number.

        That is, find the smallest positive integer `n` such
        that `x^n = 1`. If there is no such `n`, returns ``+Infinity``.

        We first check that ``abs(x)`` is very close to 1. If so, we compute
        `x` exactly and examine its argument.

        EXAMPLES::

            sage: QQbar(-sqrt(3)/2 - I/2).multiplicative_order()                        # optional - sage.symbolic
            12
            sage: QQbar(1).multiplicative_order()
            1
            sage: QQbar(-I).multiplicative_order()
            4
            sage: QQbar(707/1000 + 707/1000*I).multiplicative_order()
            +Infinity
            sage: QQbar(3/5 + 4/5*I).multiplicative_order()
            +Infinity
        """
        if not (1 in CIF(self).norm()):
            return infinity.infinity
        if self.norm() != 1:
            return infinity.infinity
        d = self.minpoly().is_cyclotomic(True)
        return d if d else infinity.infinity

    def rational_argument(self):
        r"""
        Return the argument of ``self``, divided by `2\pi`, as long as this
        result is rational. Otherwise returns ``None``. Always triggers
        exact computation.

        EXAMPLES::

            sage: QQbar((1+I)*(sqrt(2)+sqrt(5))).rational_argument()                    # optional - sage.symbolic
            1/8
            sage: QQbar(-1 + I*sqrt(3)).rational_argument()                             # optional - sage.symbolic
            1/3
            sage: QQbar(-1 - I*sqrt(3)).rational_argument()                             # optional - sage.symbolic
            -1/3
            sage: QQbar(3+4*I).rational_argument() is None
            True
            sage: (QQbar(2)**(1/5) * QQbar.zeta(7)**2).rational_argument()  # long time
            2/7
            sage: (QQbar.zeta(73)**5).rational_argument()
            5/73
            sage: (QQbar.zeta(3)^65536).rational_argument()
            1/3
        """
        # This always triggers exact computation. An alternate method
        # could almost always avoid exact computation when the result
        # is None: if we can compute an upper bound on the degree of
        # this algebraic number without exact computation, we can use
        # the method of ANExtensionElement.rational_argument().

        # Even a very loose upper bound would suffice; for instance,
        # an upper bound of 2^100, when the true degree was 8, would
        # still be efficient.

        self.exactify()
        return self._descr.rational_argument(self)

    def _pow_(self, other):
        """
        Powering for ``QQbar(1)``.

        EXAMPLES::

            sage: QQbar(1) ^ QQbar(sqrt(2))                                             # optional - sage.symbolic
            1
            sage: 1 ^ QQbar(sqrt(2))                                                    # optional - sage.symbolic
            1
            sage: QQbar(2) ^ QQbar(2)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for ^: 'Algebraic Field' and 'Algebraic Field'
            sage: AA(1) ^ AA(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for ^: 'Algebraic Real Field' and 'Algebraic Real Field'
        """
        # For some crazy unspecified reason, we must allow this if the
        # base is QQbar(1). See Issue #22120 and #24490.
        if self == 1:
            return self
        raise TypeError("unsupported operand parent(s) for ^: '{0}' and '{0}'".format(self.parent()))


class AlgebraicReal(AlgebraicNumber_base):
    r"""
    A real algebraic number.

    .. automethod:: _richcmp_
    """
    def __init__(self, x):
        """
        Create an algebraic real from x, possibly taking the real part of x.

        TESTS:

        Both of the following examples, from :trac:`11728`, trigger
        taking the real part below. This is necessary because
        sometimes a very small (e.g., 1e-17) complex part appears in a
        complex interval used to create an AlgebraicReal.::

            sage: a = QQbar((-1)^(1/4)); b = AA(a^3-a); t = b.as_number_field_element()
            sage: b*1
            -1.414213562373095?
        """
        AlgebraicNumber_base.__init__(self, AA, x)
        self._ensure_real()

    def _ensure_real(self):
        """
        This is used internally by some methods to check if
        self._value is a complex interval, and if so, take the real
        part.

        EXAMPLES::

            sage: a = QQbar((-1)^(1/4)); b = AA(a^3-a); b._value
            -1.4142135623730950488?
            sage: b._value = a._value; b._value
            0.7071067811865475244? + 0.7071067811865475244?*I
            sage: b._ensure_real()
            sage: b._value
            0.7071067811865475244?
            sage: type(b._value)
            <class 'sage.rings.real_mpfi.RealIntervalFieldElement'>
        """
        if is_ComplexIntervalFieldElement(self._value):
            self._value = self._value.real()

    def _more_precision(self):
        """
        Recompute the interval bounding this number with higher-precision
        interval arithmetic.

        EXAMPLES::

            sage: a = QQbar(sqrt(2))                                                    # optional - sage.symbolic
            sage: a._more_precision()                                                   # optional - sage.symbolic

        TESTS:

        We have to ensure after doing this that self._value is still
        real which is not the case without calling _ensure_real (see
        :trac:`11728`)::

            sage: x = polygen(ZZ, 'x')
            sage: P = AA['x'](1 + x^4); a1,a2 = P.factor()[0][0], P.factor()[1][0]; a1*a2
            x^4 + 1.000000000000000?
            sage: a1,a2
            (x^2 - 1.414213562373095?*x + 1, x^2 + 1.414213562373095?*x + 1)
            sage: a1*a2
            x^4 + 1
        """
        AlgebraicNumber_base._more_precision(self)
        self._ensure_real()

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = AA(sqrt(2))                                                       # optional - sage.symbolic
            sage: loads(dumps(t)) == t                                                  # optional - sage.symbolic
            True
        """
        return (AlgebraicReal, (self._descr, ))

    def _richcmp_(self, other, op):
        """
        Compare two algebraic reals.

        EXAMPLES::

            sage: AA(2).sqrt() < AA(3).sqrt()
            True
            sage: ((5+AA(5).sqrt())/2).sqrt() == 2*QQbar.zeta(5).imag()
            True
            sage: AA(3).sqrt() + AA(2).sqrt() < 3
            False

        TESTS::

            sage: AA(golden_ratio) < AA(sqrt(5))                                        # optional - sage.symbolic
            True
            sage: AA(golden_ratio) == AA((sqrt(5)+1)/2)                                 # optional - sage.symbolic
            True
            sage: AA(7) >= AA(50/7)
            False

        Check for trivial equality with identical elements::

            sage: x1 = AA(2^(1/50))                                                     # optional - sage.symbolic
            sage: x2 = AA(2^(1/50))                                                     # optional - sage.symbolic
            sage: y = x1 - x2                                                           # optional - sage.symbolic
            sage: y == y                                                                # optional - sage.symbolic
            True
            sage: y >= y                                                                # optional - sage.symbolic
            True
            sage: y < y                                                                 # optional - sage.symbolic
            False

            sage: z = x1 - x2                                                           # optional - sage.symbolic
            sage: z == 0                                                                # optional - sage.symbolic
            True

            sage: a = x1 - x2                                                           # optional - sage.symbolic
            sage: b = x1 - x2                                                           # optional - sage.symbolic
            sage: a == b                                                                # optional - sage.symbolic
            True
        """
        if self is other:
            return rich_to_bool(op, 0)

        # note: we can assume that self is not other here
        sd = self._descr
        od = other._descr

        # case 0: rationals
        if type(sd) is ANRational and type(od) is ANRational:
            return richcmp(sd._value, od._value, op)

        # case 1: real parts are clearly distinct
        if not self._value.overlaps(other._value):
            # NOTE: do not call richcmp here as self._value and other._value
            # might have different precisions. See
            # https://github.com/sagemath/sage/issues/29220
            return self._value._richcmp_(other._value, op)

        if op == op_EQ or op == op_NE:
            # some cheap and quite common tests where we can decide equality or difference
            if type(sd) is ANRational and not sd._value:
                return bool(other) == (op == op_NE)
            elif type(od) is ANRational and not od._value:
                return bool(self) == (op == op_NE)
            elif (type(sd) is ANExtensionElement and
                  type(od) is ANExtensionElement and
                  sd._generator is od._generator):
                return sd._value == od._value if op == op_EQ else sd._value != od._value
            else:
                # Only compare the minimal polynomials if they have been computed
                #   as otherwise it calls exactify().
                try:
                    if self._minimal_polynomial != other._minimal_polynomial:
                        return op == op_NE
                except AttributeError:
                    pass

        # case 2: possibly equal values
        # (this case happen a lot when sorting the roots of a real polynomial)
        # Only compare the minimal polynomials if they have been computed
        #   as otherwise it calls exactify().
        try:
            if self._minimal_polynomial != other._minimal_polynomial:
                c = cmp_elements_with_same_minpoly(self, other, self.minpoly())
                if c is not None:
                    return rich_to_bool(op, c)
        except AttributeError:
            pass

        if self._value.prec() < 128:
            self._more_precision()
        if other._value.prec() < 128:
            other._more_precision()
        if not self._value.overlaps(other._value):
            # NOTE: do not call richcmp here as self._value and other._value
            # might have different precisions. See
            # https://github.com/sagemath/sage/issues/29220
            return self._value._richcmp_(other._value, op)

        return rich_to_bool(op, (self - other).sign())

    def _integer_(self, Z=None):
        """
        Return self as an Integer.

        EXAMPLES::

            sage: AA(42)._integer_()
            42
            sage: AA(42)._integer_().parent()
            Integer Ring
            sage: AA(golden_ratio)._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce non-integral Algebraic Real 1.618033988749895? to Integer
            sage: (AA(golden_ratio)^10 + AA(1-golden_ratio)^10)._integer_()
            123
            sage: AA(-22/7)._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce non-integral Algebraic Real -22/7 to Integer
        """
        if self._value.lower().ceiling() > self._value.upper().floor():
            # The value is known to be non-integral.
            raise ValueError(lazy_string("Cannot coerce non-integral Algebraic Real %s to Integer", self))

        self.exactify()
        if not isinstance(self._descr, ANRational):
            raise ValueError(lazy_string("Cannot coerce irrational Algebraic Real %s to Integer", self))

        return ZZ(self._descr._value)

    def _floor_ceil(self, method):
        r"""
        Helper method used by :meth:`floor()`, :meth:`ceil()`,
        :meth:`round()`, and :meth:`trunc()`.

        TESTS::

            sage: x = polygen(QQ)
            sage: a = AA.polynomial_root(x^5 - (1-2^(-80)), RIF((0,2)))
            sage: b = AA.polynomial_root(x^5 - (1+2^(-80)), RIF((0,2)))
            sage: two = (a+b)^5 - 5*(a^4*b+a*b^4) - 10*(a^3*b^2+a^2*b^3)
            sage: one_half = 1/two
            sage: [[z.floor(), z.ceil(), z.round(), z.trunc()] # indirect doctest
            ....:  for z in [a, -a, b, -b, 6*(a+two),
            ....:            AA(0), AA(1), AA(-1), AA(1/2), AA(-1/2)]]
            [[0, 1, 1, 0], [-1, 0, -1, 0], [1, 2, 1, 1], [-2, -1, -1, -1],
            [17, 18, 18, 17], [0, 0, 0, 0], [1, 1, 1, 1], [-1, -1, -1, -1],
            [0, 1, 1, 0], [-1, 0, -1, 0]]
            sage: [[z.floor(), z.ceil(), z.trunc()] for z in [two, a*b]] # long time
            [[2, 2, 2], [0, 1, 0]]
            sage: [one_half.round(), (-one_half).round()] # long time
            [1, -1]
        """
        for i in itertools.count():
            candidate = method(self._value.lower())
            if candidate == method(self._value.upper()):
                return candidate
            self._more_precision()
            # field elements are irrational by construction
            if i == 2 and not isinstance(self._descr, ANExtensionElement):
                try:
                    return method(self._rational_())
                except (ValueError, TypeError):
                    pass

    def floor(self):
        r"""
        Return the largest integer not greater than ``self``.

        EXAMPLES::

            sage: AA(sqrt(2)).floor()                                                   # optional - sage.symbolic
            1
            sage: AA(-sqrt(2)).floor()                                                  # optional - sage.symbolic
            -2
            sage: AA(42).floor()
            42

        TESTS:

        Check that :trac:`15501` is fixed::

            sage: a = QQbar((-1)^(1/4)).real()                                          # optional - sage.symbolic
            sage: (floor(a-a) + a).parent()                                             # optional - sage.symbolic
            Algebraic Real Field
        """
        return self._floor_ceil(lambda x: x.floor())

    def ceil(self):
        r"""
        Return the smallest integer not smaller than ``self``.

        EXAMPLES::

            sage: AA(sqrt(2)).ceil()                                                    # optional - sage.symbolic
            2
            sage: AA(-sqrt(2)).ceil()                                                   # optional - sage.symbolic
            -1
            sage: AA(42).ceil()
            42
        """
        return self._floor_ceil(lambda x: x.ceil())

    def round(self):
        r"""
        Round ``self`` to the nearest integer.

        EXAMPLES::

            sage: AA(sqrt(2)).round()                                                   # optional - sage.symbolic
            1
            sage: AA(1/2).round()
            1
            sage: AA(-1/2).round()
            -1
        """
        return self._floor_ceil(lambda x: x.round())

    def trunc(self):
        r"""
        Round ``self`` to the nearest integer toward zero.

        EXAMPLES::

            sage: AA(sqrt(2)).trunc()                                                   # optional - sage.symbolic
            1
            sage: AA(-sqrt(2)).trunc()                                                  # optional - sage.symbolic
            -1
            sage: AA(1).trunc()
            1
            sage: AA(-1).trunc()
            -1
        """
        return self._floor_ceil(lambda x: x.trunc())

    def _rational_(self):
        """
        Return self as a Rational.

        EXAMPLES::

            sage: AA(42)._rational_().parent()
            Rational Field
            sage: AA(-22/7)._rational_()
            -22/7
            sage: AA(sqrt(7))._rational_()                                              # optional - sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce irrational Algebraic Real 2.645751311064591? to Rational
            sage: v = AA(1/2 + sqrt(2))^3 - AA(11/4*sqrt(2)); v                         # optional - sage.symbolic
            3.125000000000000?
            sage: v._rational_()
            25/8
        """
        self.exactify()
        if not isinstance(self._descr, ANRational):
            raise ValueError(lazy_string("Cannot coerce irrational Algebraic Real %s to Rational", self))

        return QQ(self._descr._value)

    def real(self):
        """
        Return the real part of this algebraic real.

        It always returns ``self``.

        EXAMPLES::

            sage: a = AA(sqrt(2) + sqrt(3))                                             # optional - sage.symbolic
            sage: a.real()                                                              # optional - sage.symbolic
            3.146264369941973?
            sage: a.real() is a                                                         # optional - sage.symbolic
            True
        """
        return self

    def imag(self):
        """
        Return the imaginary part of this algebraic real.

        It always returns 0.

        EXAMPLES::

            sage: a = AA(sqrt(2) + sqrt(3))                                             # optional - sage.symbolic
            sage: a.imag()                                                              # optional - sage.symbolic
            0
            sage: parent(a.imag())                                                      # optional - sage.symbolic
            Algebraic Real Field
        """
        return AA_0

    def conjugate(self):
        """
        Return the complex conjugate of ``self``, i.e. returns itself.

        EXAMPLES::

            sage: a = AA(sqrt(2) + sqrt(3))                                             # optional - sage.symbolic
            sage: a.conjugate()                                                         # optional - sage.symbolic
            3.146264369941973?
            sage: a.conjugate() is a                                                    # optional - sage.symbolic
            True
        """
        return self

    def multiplicative_order(self):
        r"""
        Compute the multiplicative order of this real algebraic number.

        That is, find the smallest positive integer `n` such
        that `x^n = 1`. If there is no such `n`, returns ``+Infinity``.

        We first check that ``abs(x)`` is very close to 1. If so, we compute
        `x` exactly and compare it to `1` and `-1`.

        EXAMPLES::

            sage: AA(1).multiplicative_order()
            1
            sage: AA(-1).multiplicative_order()
            2
            sage: AA(5).sqrt().multiplicative_order()
            +Infinity
        """
        if not (1 in RIF(self).abs()):
            return infinity.infinity
        if self == 1:
            return 1
        elif self == -1:
            return 2
        else:
            return infinity.infinity

    def sign(self):
        """
        Compute the sign of this algebraic number (return `-1` if negative,
        `0` if zero, or `1` if positive).

        This computes an interval enclosing this number using 128-bit interval
        arithmetic; if this interval includes 0, then fall back to
        exact computation (which can be very slow).

        EXAMPLES::

            sage: AA(-5).nth_root(7).sign()
            -1
            sage: (AA(2).sqrt() - AA(2).sqrt()).sign()
            0

            sage: a = AA(2).sqrt() + AA(3).sqrt() - 58114382797550084497/18470915334626475921
            sage: a.sign()
            1
            sage: b = AA(2).sqrt() + AA(3).sqrt() - 2602510228533039296408/827174681630786895911
            sage: b.sign()
            -1

            sage: c = AA(5)**(1/3) - 1437624125539676934786/840727688792155114277
            sage: c.sign()
            1

            sage: (((a+b)*(a+c)*(b+c))**9 / (a*b*c)).sign()
            1
            sage: (a-b).sign()
            1
            sage: (b-a).sign()
            -1
            sage: (a*b).sign()
            -1
            sage: ((a*b).abs() + a).sign()
            1
            sage: (a*b - b*a).sign()
            0

            sage: a = AA(sqrt(1/2))                                                     # optional - sage.symbolic
            sage: b = AA(-sqrt(1/2))                                                    # optional - sage.symbolic
            sage: (a + b).sign()                                                        # optional - sage.symbolic
            0

        TESTS:

        We avoid calling :meth:`exactify()` for trivial differences. The
        following example will take a long time (more than 5 seconds)
        when calling ``y.exactify()``::

            sage: x1 = AA(2^(1/50))                                                     # optional - sage.symbolic
            sage: x2 = AA(2^(1/50))                                                     # optional - sage.symbolic
            sage: y = x1 - x2                                                           # optional - sage.symbolic
            sage: y.sign()                                                              # optional - sage.symbolic
            0

        Simplify to rationals for binary operations when computing the sign::

            sage: a = AA(2^(1/60))                                                      # optional - sage.symbolic
            sage: b = a - (a + 1)                                                       # optional - sage.symbolic
            sage: (b + 1).sign()                                                        # optional - sage.symbolic
            0
        """
        if not self._value.contains_zero():
            return self._value.unique_sign()

        sd = self._descr
        if isinstance(self._descr, ANRational):
            return sd._value.sign()
        elif isinstance(self._descr, ANExtensionElement):
            # All field elements are irrational by construction
            # (the ANExtensionElement constructor will return an ANRational
            # instead, if the number is actually rational).
            # An irrational number must eventually be different from 0
            while self._value.contains_zero():
                self._more_precision()
            return self._value.unique_sign()
        elif type(sd) is ANBinaryExpr:
            ls = sd._left.sign()
            rs = sd._right.sign()
            if sd._op is operator.mul or sd._op is operator.truediv:
                return ls * rs
            elif sd._op is operator.add:
                if ls == rs:
                    return ls
            else:
                if ls == -rs:
                    return ls
                elif not ls:
                    self._set_descr((-sd._right)._descr)
                    return -rs
                elif not rs:
                    self._set_descr(sd._left._descr)
                    return ls
                elif sd._left is sd._right:
                    self._set_descr(ANRational(QQ.zero()))
                    return 0
        elif type(sd) is ANUnaryExpr:
            if sd._op == 'abs':
                c = 1 if bool(sd._arg) else 0
                if not c:
                    self._set_descr(ANRational(QQ.zero()))
                return c
            elif sd._op == '-':
                return -(sd._arg.sign())
            elif sd._op == '~':
                return sd._arg.sign()

        if self._value.prec() < 128:
            # OK, we'll try adding precision one more time
            self._more_precision()
            if not self._value.contains_zero():
                return self._value.unique_sign()

        if type(sd) is ANBinaryExpr:
            # We will now exactify both sides and do another sign comparison.
            # We try to avoid making ourself exact if possible.
            # It will only reach this block if the operation is addition or subtraction.
            sd._left.exactify()
            sd._right.exactify()

            # Rationals
            if type(sd._left._descr) is ANRational and type(sd._right._descr) is ANRational:
                ret = sd._op(sd._left._descr._value, sd._right._descr._value)
                if ret == 0:
                    self._set_descr(ANRational(QQ.zero()))
                    return 0
                return ret.sign()

            if sd._left.minpoly() == sd._right.minpoly():
                # Negating the element does not change the minpoly
                right = sd._right if sd._op is operator.sub else -sd._right
                c = cmp_elements_with_same_minpoly(sd._left, right, sd._left.minpoly())
                if c == 0:
                    self._set_descr(ANRational(QQ.zero()))
                    return 0
                elif c is not None:
                    return c

            ret = sd._op(sd._left._value, sd._right._value)
            if not ret.contains_zero():
                return ret.unique_sign()
            if not ret:  # Known to be exactly 0
                self._set_descr(ANRational(QQ.zero()))
                return 0

        # Sigh...
        self.exactify()
        return self.sign()

    def _interval_fast(self, prec):
        r"""
        Compute an approximation to this ``AlgebraicReal`` object in a real interval field of precision prec.

        EXAMPLES::

            sage: t = AA(sqrt(7))                                                       # optional - sage.symbolic
            sage: t._interval_fast(100)                                                 # optional - sage.symbolic
            2.64575131106459059050161575364?
        """
        return self.interval_fast(RealIntervalField(prec))

    def interval_exact(self, field):
        """
        Given a :class:`RealIntervalField`, compute the best possible
        approximation of this number in that field. Note that if this
        number is sufficiently close to some floating-point number
        (and, in particular, if this number is exactly representable in
        floating-point), then this will trigger exact computation, which
        may be very slow.

        EXAMPLES::

            sage: x = AA(2).sqrt()
            sage: y = x*x
            sage: x.interval(RIF)
            1.414213562373095?
            sage: x.interval_exact(RIF)
            1.414213562373095?
            sage: y.interval(RIF)
            2.000000000000000?
            sage: y.interval_exact(RIF)
            2
            sage: z = 1 + AA(2).sqrt() / 2^200
            sage: z.interval(RIF)
            1.000000000000001?
            sage: z.interval_exact(RIF)
            1.000000000000001?

        TESTS:

        Check that :trac:`26898` is fixed.  This calculation triggers the 40 bits
        of extra precision below, and the point is not that the length of the list
        is seven, but that the code runs in a reasonable time::

            sage: R.<x> = QQbar[]
            sage: roots = (x^7 + 27/4).roots()
            sage: from sage.rings.qqbar import QQbar_hash_offset
            sage: len([(r[0] + QQbar_hash_offset).interval_exact(CIF) for r in roots])
            7

        """
        for extra in (0, 40):
            target = RR(1.0) >> (field.prec() + extra)
            # p==precise; pr==precise rounded
            pval = self.interval_diameter(target)
            pbot = pval.lower()
            ptop = pval.upper()
            val = field(pval)
            bot = val.lower()
            top = val.upper()
            prbot = pbot.parent()(bot)
            prtop = ptop.parent()(top)
            if bot == top or (bot.nextabove() == top and
                              prbot < pbot and ptop < prtop):
                return val

        # Even 40 extra bits of precision are not enough to prove that
        # self is not an exactly representable float.
        self.exactify()
        while True:
            # p==precise; pr==precise rounded
            pval = self._value
            pbot = pval.lower()
            ptop = pval.upper()
            val = field(pval)
            bot = val.lower()
            top = val.upper()
            prbot = pbot.parent()(bot)
            prtop = ptop.parent()(top)
            if bot == top or (bot.nextabove() == top and
                              prbot < pbot and ptop < prtop):
                return val

            self._more_precision()

    def real_number(self, field):
        """
        Given a :class:`RealField`, compute a good approximation to ``self`` in
        that field. The approximation will be off by at most two
        ulp's, except for numbers which are very close to 0, which
        will have an absolute error at most
        ``2**(-(field.prec()-1))``. Also, the rounding mode of the
        field is respected.

        EXAMPLES::

            sage: x = AA(2).sqrt()^2
            sage: x.real_number(RR)
            2.00000000000000
            sage: x.real_number(RealField(53, rnd='RNDD'))
            1.99999999999999
            sage: x.real_number(RealField(53, rnd='RNDU'))
            2.00000000000001
            sage: x.real_number(RealField(53, rnd='RNDZ'))
            1.99999999999999
            sage: (-x).real_number(RR)
            -2.00000000000000
            sage: (-x).real_number(RealField(53, rnd='RNDD'))
            -2.00000000000001
            sage: (-x).real_number(RealField(53, rnd='RNDU'))
            -1.99999999999999
            sage: (-x).real_number(RealField(53, rnd='RNDZ'))
            -1.99999999999999
            sage: (x-2).real_number(RR)
            5.42101086242752e-20
            sage: (x-2).real_number(RealField(53, rnd='RNDD'))
            -1.08420217248551e-19
            sage: (x-2).real_number(RealField(53, rnd='RNDU'))
            2.16840434497101e-19
            sage: (x-2).real_number(RealField(53, rnd='RNDZ'))
            0.000000000000000
            sage: y = AA(2).sqrt()
            sage: y.real_number(RR)
            1.41421356237309
            sage: y.real_number(RealField(53, rnd='RNDD'))
            1.41421356237309
            sage: y.real_number(RealField(53, rnd='RNDU'))
            1.41421356237310
            sage: y.real_number(RealField(53, rnd='RNDZ'))
            1.41421356237309
        """
        v = self.interval(RealIntervalField(field.prec()))
        return field(v)

    _mpfr_ = real_number

    def __float__(self):
        r"""
        Compute a good float approximation to self.

        EXAMPLES::

            sage: AA(golden_ratio).__float__()
            1.618033988749895
            sage: float(AA(sqrt(11)))
            3.3166247903554
        """
        return float(RR(self))

    def _complex_mpfr_field_(self, field):
        r"""
        Compute an approximation to this ``AlgebraicReal`` in the given field,
        which may be an interval field (in which case ``self.interval()`` is
        called) or any other real number field (in which case
        ``self.real_number()`` is called.

        Note that the field ``field`` should be a *complex* field (whose
        ``_real_field()`` method will be called to obtain a real subfield.)

        EXAMPLES::

            sage: AA(golden_ratio)._complex_mpfr_field_(ComplexIntervalField(100))
            1.618033988749894848204586834365?
            sage: AA(golden_ratio)._complex_mpfr_field_(ComplexField(100))
            1.6180339887498948482045868344
        """
        if isinstance(field, sage.rings.abc.ComplexIntervalField):
            return field(self.interval(field._real_field()))
        else:
            return field(self.real_number(field._real_field()))

    def real_exact(self, field):
        r"""
        Given a :class:`RealField`, compute the best possible approximation of
        this number in that field. Note that if this number is
        sufficiently close to the halfway point between two
        floating-point numbers in the field (for the default
        round-to-nearest mode) or if the number is sufficiently close
        to a floating-point number in the field (for directed rounding
        modes), then this will trigger exact computation, which may be
        very slow.

        The rounding mode of the field is respected.

        EXAMPLES::

            sage: x = AA(2).sqrt()^2
            sage: x.real_exact(RR)
            2.00000000000000
            sage: x.real_exact(RealField(53, rnd='RNDD'))
            2.00000000000000
            sage: x.real_exact(RealField(53, rnd='RNDU'))
            2.00000000000000
            sage: x.real_exact(RealField(53, rnd='RNDZ'))
            2.00000000000000
            sage: (-x).real_exact(RR)
            -2.00000000000000
            sage: (-x).real_exact(RealField(53, rnd='RNDD'))
            -2.00000000000000
            sage: (-x).real_exact(RealField(53, rnd='RNDU'))
            -2.00000000000000
            sage: (-x).real_exact(RealField(53, rnd='RNDZ'))
            -2.00000000000000
            sage: y = (x-2).real_exact(RR).abs()
            sage: y == 0.0 or y == -0.0 # the sign of 0.0 is not significant in MPFI
            True
            sage: y = (x-2).real_exact(RealField(53, rnd='RNDD'))
            sage: y == 0.0 or y == -0.0 # same as above
            True
            sage: y = (x-2).real_exact(RealField(53, rnd='RNDU'))
            sage: y == 0.0 or y == -0.0 # idem
            True
            sage: y = (x-2).real_exact(RealField(53, rnd='RNDZ'))
            sage: y == 0.0 or y == -0.0 # ibidem
            True
            sage: y = AA(2).sqrt()
            sage: y.real_exact(RR)
            1.41421356237310
            sage: y.real_exact(RealField(53, rnd='RNDD'))
            1.41421356237309
            sage: y.real_exact(RealField(53, rnd='RNDU'))
            1.41421356237310
            sage: y.real_exact(RealField(53, rnd='RNDZ'))
            1.41421356237309
        """
        for extra in (0, 40):
            target = RR(1.0) >> (field.prec() + extra)
            val = self.interval_diameter(target)
            fbot = field(val.lower())
            ftop = field(val.upper())
            if fbot == ftop:
                return ftop

        # Even 40 extra bits of precision are not enough to determine the
        # answer.
        rifp1 = RealIntervalField(field.prec() + 1)
        rifp2 = RealIntervalField(field.prec() + 2)

        val = self.interval_exact(rifp1)

        # Call the largest floating-point number <= self 'x'. Then
        # val may be [x .. x], [x .. x + 1/2 ulp],
        # [x + 1/2 ulp .. x + 1/2 ulp], or [x + 1/2 ulp .. x + 1 ulp];
        # in the second and fourth cases, the true value is not equal
        # to either of the interval endpoints.

        mid = rifp2(val).center()

        # Now mid may be x, x + 1/4 ulp, x + 1/2 ulp, or x + 3/4 ulp; in
        # the first and third cases, mid is the exact, true value of self;
        # in the second and fourth cases, self is close to mid, and is
        # neither x, x + 1/2 ulp, nor x + 1 ulp.

        # In all of these cases, in all rounding modes, the rounded value
        # of mid is the same as the rounded value of self.

        return field(mid)


class AlgebraicNumberPowQQAction(Action):
    """
    Implement powering of an algebraic number (an element of ``QQbar``
    or ``AA``) by a rational.

    This is always a right action.

    INPUT:

    - ``G`` -- must be ``QQ``

    - ``S`` -- the parent on which to act, either ``AA`` or ``QQbar``.

    .. NOTE::

        To compute ``x ^ (a/b)``, we take the `b`'th root of `x`; then
        we take that to the `a`'th power. If `x` is a negative algebraic
        real and `b` is odd, take the real `b`'th root; otherwise take
        the principal `b`'th root.

    EXAMPLES:

    In ``QQbar``::

        sage: QQbar(2)^(1/2)
        1.414213562373095?
        sage: QQbar(8)^(2/3)
        4
        sage: QQbar(8)^(2/3) == 4
        True
        sage: x = polygen(QQbar)
        sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(1, 2))
        sage: tau = QQbar.polynomial_root(x^2 - x - 1, RIF(-1, 0))
        sage: rt5 = QQbar(5)^(1/2)
        sage: phi^10 / rt5
        55.00363612324742?
        sage: tau^10 / rt5
        0.003636123247413266?
        sage: (phi^10 - tau^10) / rt5
        55.00000000000000?
        sage: (phi^10 - tau^10) / rt5 == fibonacci(10)
        True
        sage: (phi^50 - tau^50) / rt5 == fibonacci(50)
        True
        sage: QQbar(-8)^(1/3)
        1.000000000000000? + 1.732050807568878?*I
        sage: (QQbar(-8)^(1/3))^3
        -8
        sage: QQbar(32)^(1/5)
        2
        sage: a = QQbar.zeta(7)^(1/3); a
        0.9555728057861407? + 0.2947551744109043?*I
        sage: a == QQbar.zeta(21)
        True
        sage: QQbar.zeta(7)^6
        0.6234898018587335? - 0.7818314824680299?*I
        sage: (QQbar.zeta(7)^6)^(1/3) * QQbar.zeta(21)
        1.000000000000000? + 0.?e-17*I

    In ``AA``::

        sage: AA(2)^(1/2)
        1.414213562373095?
        sage: AA(8)^(2/3)
        4
        sage: AA(8)^(2/3) == 4
        True
        sage: x = polygen(AA)
        sage: phi = AA.polynomial_root(x^2 - x - 1, RIF(0, 2))
        sage: tau = AA.polynomial_root(x^2 - x - 1, RIF(-2, 0))
        sage: rt5 = AA(5)^(1/2)
        sage: phi^10 / rt5
        55.00363612324742?
        sage: tau^10 / rt5
        0.003636123247413266?
        sage: (phi^10 - tau^10) / rt5
        55.00000000000000?
        sage: (phi^10 - tau^10) / rt5 == fibonacci(10)
        True
        sage: (phi^50 - tau^50) / rt5 == fibonacci(50)
        True

    TESTS::

        sage: AA(-8)^(1/3)
        -2
        sage: AA(-8)^(2/3)
        4
        sage: AA(32)^(3/5)
        8
        sage: AA(-16)^(1/2)
        4*I
        sage: AA(-16)^(1/4)
        1.414213562373095? + 1.414213562373095?*I
        sage: AA(-16)^(1/4)/QQbar.zeta(8)
        2

    We check that :trac:`7859` is fixed::

        sage: (AA(2)^(1/2)-AA(2)^(1/2))^(1/2)
        0
    """
    def __init__(self, G, S):
        """
        EXAMPLES::

            sage: from sage.rings.qqbar import AlgebraicNumberPowQQAction
            sage: act = AlgebraicNumberPowQQAction(QQ, AA); act
            Right Rational Powering by Rational Field on Algebraic Real Field
            sage: act(AA(-2), 1/3)
            -1.259921049894873?

        ::

            sage: act = AlgebraicNumberPowQQAction(QQ, QQbar); act
            Right Rational Powering by Rational Field on Algebraic Field
            sage: act(QQbar(-2), 1/3)
            0.6299605249474365? + 1.091123635971722?*I
        """
        Action.__init__(self, G, S, False, operator.pow)

    def _act_(self, e, x):
        r"""
        Return the power ``x ^ e``.

        INPUT:

        - ``x`` -- an algebraic number

        - ``e`` -- a rational number
        """
        if not x:
            return x

        n = e.numerator()
        d = e.denominator()
        if d == 1:
            return x._pow_int(n)

        # Parent of the result
        S = self.codomain()
        if S is AA and d % 2 == 0 and x.sign() < 0:
            S = QQbar

        # First, check for exact roots.
        if isinstance(x._descr, ANRational):
            rt = rational_exact_root(abs(x._descr._value), d)
            if rt is not None:
                if x._descr._value < 0:
                    if S is AA:
                        return AlgebraicReal(ANRational((-rt)**n))
                    else:
                        z = QQbar.zeta(2 * d)._pow_int(n)
                        return z * AlgebraicNumber(ANRational(rt**n))
                return S(ANRational(rt**n))

        if S is AA:
            # Result lies in AA
            pow_n = x._pow_int(n)
            poly = AAPoly.gen()**d - pow_n
            range = pow_n.interval_fast(RIF)
            if d % 2 == 0:
                result_min = 0
            else:
                result_min = min(range.lower(), -1)
            result_max = max(range.upper(), 1)
            return AlgebraicReal(ANRoot(poly, RIF(result_min, result_max)))

        # Result lies in QQbar

        # Determine whether arg(x) equals pi.
        argument_is_pi = False
        for prec in short_prec_seq():
            if prec is None:
                # We know that x.real() < 0, since x._value
                # crosses the negative real line and x._value
                # is known to be non-zero.
                isgn = x.imag().sign()
                val = x._value
                argument = val.argument()
                if isgn == 0:
                    argument = argument.parent().pi()
                    argument_is_pi = True
                elif isgn > 0:
                    if argument < 0:
                        argument = argument + 2 * argument.parent().pi()
                else:
                    if argument > 0:
                        argument = argument - 2 * argument.parent().pi()
            else:
                val = x._interval_fast(prec)
                if is_RealIntervalFieldElement(val) or not val.crosses_log_branch_cut():
                    argument = val.argument()
                    if val.imag().is_zero() and val.real() < 0:
                        argument_is_pi = True
                    break

        target_abs = abs(val) ** e
        target_arg = argument * e

        for prec in tail_prec_seq():
            if target_abs.relative_diameter() < RR_1_10 and (target_arg * d).absolute_diameter() < RR_1_10:
                break

            val = x._interval_fast(prec)

            target_abs = abs(val) ** e
            argument = val.argument()
            if argument_is_pi:
                argument = argument.parent().pi()
            target_arg = argument * e

        pow_n = x**n
        poly = QQbarPoly.gen()**d - pow_n

        prec = target_abs.prec()
        if argument_is_pi and d == 2:
            target_real = 0
        else:
            target_real = target_arg.cos() * target_abs
        target = ComplexIntervalField(prec)(target_real,
                                            target_arg.sin() * target_abs)

        return AlgebraicNumber(ANRoot(poly, target))

    def _repr_name_(self):
        return "Rational Powering"


class ANRational(ANDescr):
    r"""
    The subclass of :class:`ANDescr` that represents an arbitrary
    rational. This class is private, and should not be used directly.
    """

    def __init__(self, x):
        """
        TESTS::

            sage: polygen(QQbar) / int(3)
            1/3*x
        """
        if isinstance(x, (sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            self._value = x
        elif isinstance(x, int):
            self._value = ZZ(x)
        else:
            raise TypeError("Illegal initializer for algebraic number rational")

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = AA(5/2); type(t._descr)
            <class 'sage.rings.qqbar.ANRational'>
            sage: loads(dumps(t)) == t
            True
        """
        return (ANRational, (self._value, ))

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: QQbar(2/3)._repr_()
            '2/3'
        """
        return repr(self._value)

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always
        False, for rationals).

        EXAMPLES::

            sage: sage_input(QQbar(22/7), verify=True)
            # Verified
            QQbar(22/7)
            sage: sage_input(-AA(3)/5, verify=True)
            # Verified
            AA(-3/5)
            sage: sage_input(vector(AA, (0, 1/2, 1/3)), verify=True)
            # Verified
            vector(AA, [0, 1/2, 1/3])
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: rat = ANRational(9/10)
            sage: rat.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({binop:/ {atomic:9} {atomic:10}})}, False)
        """
        v = sib(self._value, True)
        if not coerce:
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)
        return (v, False)

    def _interval_fast(self, prec):
        r"""
        Return an approximation to self in a real interval field of precision prec.

        EXAMPLES::

            sage: QQbar(355/113)._descr._interval_fast(30)
            3.14159292?
        """
        return RealIntervalField(prec)(self._value)

    def generator(self):
        r"""
        Return an :class:`AlgebraicGenerator` object associated to this
        element. Returns the trivial generator, since ``self`` is rational.

        EXAMPLES::

            sage: QQbar(0)._descr.generator()
            Trivial generator
        """
        return qq_generator

    def is_complex(self):
        r"""
        Return ``False``, since rational numbers are real

        EXAMPLES::

            sage: QQbar(1/7)._descr.is_complex()
            False
        """
        return False

    def exactify(self):
        r"""
        Calculate ``self`` exactly. Since ``self`` is a rational number, return ``self``.

        EXAMPLES::

            sage: a = QQbar(1/3)._descr
            sage: a.exactify() is a
            True
        """
        return self

    def is_simple(self):
        """
        Checks whether this descriptor represents a value with the same
        algebraic degree as the number field associated with the descriptor.

        This is always true for rational numbers.

        EXAMPLES::

            sage: AA(1/2)._descr.is_simple()
            True
        """
        return True

    def minpoly(self):
        r"""
        Return the min poly of ``self`` over `\QQ`.

        EXAMPLES::

            sage: QQbar(7)._descr.minpoly()
            x - 7
        """
        return QQx_x - self._value

    def neg(self, n):
        r"""
        Negation of ``self``.

        EXAMPLES::

            sage: a = QQbar(3)
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANRational'>
            sage: b.neg(a)
            -3
        """
        return ANRational(-self._value)

    def invert(self, n):
        r"""
        1/``self``.

        EXAMPLES::

            sage: a = QQbar(3)
            sage: b = a._descr
            sage: b.invert(a)
            1/3
        """
        return ANRational(~self._value)

    def abs(self, n):
        r"""
        Absolute value of ``self``.

        EXAMPLES::

            sage: a = QQbar(3)
            sage: b = a._descr
            sage: b.abs(a)
            3
        """
        return ANRational(abs(self._value))

    def rational_argument(self, n):
        r"""
        Return the argument of self divided by `2 \pi`, or ``None`` if this
        element is 0.

        EXAMPLES::

            sage: QQbar(3)._descr.rational_argument(None)
            0
            sage: QQbar(-3)._descr.rational_argument(None)
            1/2
            sage: QQbar(0)._descr.rational_argument(None) is None
            True
        """
        if self._value > 0:
            return QQ.zero()
        if self._value < 0:
            return QQ((1, 2))
        return None

    def angle(self):
        r"""
        Return a rational number `q \in (-1/2, 1/2]` such that ``self`` is a rational multiple of
        `e^{2\pi i q}`. Always returns 0, since this element is rational.

        EXAMPLES::

            sage: QQbar(3)._descr.angle()
            0
            sage: QQbar(-3)._descr.angle()
            0
            sage: QQbar(0)._descr.angle()
            0
        """
        return QQ_0

    def scale(self):
        r"""
        Return a rational number `r` such that ``self`` is equal to `r e^{2 \pi
        i q}` for some `q \in (-1/2, 1/2]`.  In other words, just return ``self``
        as a rational number.

        EXAMPLES::

            sage: QQbar(-3)._descr.scale()
            -3
        """
        return self._value


def is_AlgebraicReal(x):
    r"""
    Test if ``x`` is an instance of :class:`~AlgebraicReal`. For internal use.

    EXAMPLES::

        sage: from sage.rings.qqbar import is_AlgebraicReal
        sage: is_AlgebraicReal(AA(sqrt(2)))
        True
        sage: is_AlgebraicReal(QQbar(sqrt(2)))
        False
        sage: is_AlgebraicReal("spam")
        False
    """
    return isinstance(x, AlgebraicReal)


def is_AlgebraicNumber(x):
    r"""
    Test if ``x`` is an instance of :class:`~AlgebraicNumber`. For internal use.

    EXAMPLES::

        sage: from sage.rings.qqbar import is_AlgebraicNumber
        sage: is_AlgebraicNumber(AA(sqrt(2)))
        False
        sage: is_AlgebraicNumber(QQbar(sqrt(2)))
        True
        sage: is_AlgebraicNumber("spam")
        False
    """
    return isinstance(x, AlgebraicNumber)


QQbarPoly = PolynomialRing(QQbar, 'x')
AAPoly = PolynomialRing(AA, 'x')


class AlgebraicPolynomialTracker(SageObject):
    r"""
    Keeps track of a polynomial used for algebraic numbers.

    If multiple algebraic numbers are created as roots of a single
    polynomial, this allows the polynomial and information about
    the polynomial to be shared. This reduces work if the polynomial
    must be recomputed at higher precision, or if it must be factored.

    This class is private, and should only be constructed by
    ``AA.common_polynomial()`` or ``QQbar.common_polynomial()``, and should
    only be used as an argument to ``AA.polynomial_root()`` or
    ``QQbar.polynomial_root()``. (It does not matter whether you create
    the common polynomial with ``AA.common_polynomial()`` or
    ``QQbar.common_polynomial()``.)

    EXAMPLES::

        sage: x = polygen(QQbar)
        sage: P = QQbar.common_polynomial(x^2 - x - 1)
        sage: P
        x^2 - x - 1
        sage: QQbar.polynomial_root(P, RIF(1, 2))
        1.618033988749895?
    """

    def __init__(self, poly):
        r"""
        Initialize this AlgebraicPolynomialTracker object.

        EXAMPLES::

            sage: x = polygen(QQbar)
            sage: P = QQbar.common_polynomial(x^2 - x - 1)
            sage: type(P) # indirect doctest
            <class 'sage.rings.qqbar.AlgebraicPolynomialTracker'>
        """
        if not isinstance(poly, Polynomial):
            raise ValueError("Trying to create AlgebraicPolynomialTracker on non-Polynomial")
        B = poly.base_ring()

        if B is ZZ or B is QQ:
            poly = QQy(poly)
            complex = False
        elif isinstance(B, AlgebraicField_common):
            complex = isinstance(poly.base_ring(), AlgebraicField)
        else:
            try:
                poly = poly.change_ring(AA)
                complex = False
            except (TypeError, ValueError):
                poly = poly.change_ring(QQbar)
                complex = True
        self._poly = poly
        self._complex = complex
        self._exact = False
        self._roots_cache = {}

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: type(v._descr._poly)
            <class 'sage.rings.qqbar.AlgebraicPolynomialTracker'>
            sage: loads(dumps(v)) == v
            True
        """
        return (AlgebraicPolynomialTracker, (self._poly, ))

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: sage_input(AA.common_polynomial(x^3 - 7))
            R.<y> = QQ[]
            AA.common_polynomial(y^3 - 7)
            sage: x = polygen(AA)
            sage: p = sqrt(AA(2)) * x^2 - sqrt(AA(3))
            sage: cp = AA.common_polynomial(p)
            sage: sage_input((cp, cp))
            R.<x> = AA[]
            cp = AA.common_polynomial(AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))*x^2 - AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774))))
            (cp, cp)
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: cp._sage_input_(sib, False)
            {call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:* {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.4142135623730949})}, {call: {atomic:RR}({atomic:1.4142135623730951})})})} {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}}} {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:3}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.7320508075688772})}, {call: {atomic:RR}({atomic:1.7320508075688774})})})}})}
        """
        # XXX It would be nicer to skip the "AA.common_polynomial()"
        # wrapper if the polynomial is not actually shared. But
        # sage_input.py is not quite that generic.
        v = sib.name('AA').common_polynomial(self._poly)
        sib.id_cache(self, v, 'cp')
        return v

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: AA.common_polynomial(x^3 - 7)._repr_()
            'y^3 - 7'
        """
        return repr(self._poly)

    def poly(self):
        r"""
        Return the underlying polynomial of ``self``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^3 - 7
            sage: g = AA.common_polynomial(f)
            sage: g.poly()
            y^3 - 7
        """
        return self._poly

    def is_complex(self):
        r"""
        Return ``True`` if the coefficients of this polynomial are non-real.

        EXAMPLES::

            sage: x = polygen(QQ); f = x^3 - 7
            sage: g = AA.common_polynomial(f)
            sage: g.is_complex()
            False
            sage: QQbar.common_polynomial(x^3 - QQbar(I)).is_complex()
            True
        """
        return self._complex

    def complex_roots(self, prec, multiplicity):
        """
        Find the roots of ``self`` in the complex field to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: cp = AA.common_polynomial(x^4 - 2)

        Note that the precision is not guaranteed to find the tightest
        possible interval since :meth:`complex_roots` depends on the
        underlying BLAS implementation. ::

            sage: cp.complex_roots(30, 1)
            [-1.18920711500272...?,
             1.189207115002721?,
             -1.189207115002721?*I,
             1.189207115002721?*I]
        """
        if multiplicity in self._roots_cache:
            roots = self._roots_cache[multiplicity]
            if roots[0] >= prec:
                return roots[1]

        p = self._poly
        for i in range(multiplicity - 1):
            p = p.derivative()

        from sage.rings.polynomial.complex_roots import complex_roots
        roots_mult = complex_roots(p, min_prec=prec)
        roots = [rt for (rt, mult) in roots_mult if mult == 1]
        self._roots_cache[multiplicity] = (prec, roots)
        return roots

    def exactify(self):
        """
        Compute a common field that holds all of the algebraic coefficients
        of this polynomial, then factor the polynomial over that field.
        Store the factors for later use (ignoring multiplicity).

        EXAMPLES::

            sage: x = polygen(AA)
            sage: p = sqrt(AA(2)) * x^2 - sqrt(AA(3))
            sage: cp = AA.common_polynomial(p)
            sage: cp._exact
            False
            sage: cp.exactify()
            sage: cp._exact
            True
        """
        if self._exact:
            return

        self._exact = True

        if self._poly.base_ring() is QQ:
            self._factors = [fac_exp[0] for fac_exp in self._poly.factor()]
            self._gen = qq_generator
        else:
            gen = qq_generator
            for c in self._poly.list():
                c.exactify()
                gen = gen.union(c._exact_field())

            self._gen = gen

            coeffs = [gen(c._exact_value()) for c in self._poly.list()]

            if gen.is_trivial():
                self._poly = QQy(self._poly)
                self._factors = [fac_exp[0] for fac_exp in self._poly.factor()]
            else:
                fld = gen.field()
                fld_poly = fld['y']

                fp = fld_poly(coeffs)

                self._factors = [fac_exp[0] for fac_exp in fp.factor()]

    def factors(self):
        r"""
        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = QQbar.common_polynomial(x^4 + 4)
            sage: f.factors()
            [y^2 - 2*y + 2, y^2 + 2*y + 2]
        """
        self.exactify()
        return self._factors

    def generator(self):
        r"""
        Return an :class:`AlgebraicGenerator` for a number field containing all
        the coefficients of ``self``.

        EXAMPLES::

            sage: x = polygen(AA)
            sage: p = sqrt(AA(2)) * x^2 - sqrt(AA(3))
            sage: cp = AA.common_polynomial(p)
            sage: cp.generator()
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1
             with a in -0.5176380902050415?
        """
        self.exactify()
        return self._gen


class ANRoot(ANDescr):
    """
    The subclass of :class:`ANDescr` that represents a particular
    root of a polynomial with algebraic coefficients.
    This class is private, and should not be used directly.
    """
    def __init__(self, poly, interval, multiplicity=1):
        r"""
        Initialize this ``ANRoot`` object.

        EXAMPLES::

            sage: x = polygen(QQ); f = (x^3 + x + 1).roots(AA,multiplicities=False)[0]._descr
            sage: type(f) # indirect doctest
            <class 'sage.rings.qqbar.ANRoot'>
        """
        if not isinstance(poly, AlgebraicPolynomialTracker):
            poly = AlgebraicPolynomialTracker(poly)
        self._poly = poly
        self._multiplicity = multiplicity
        self._complex = is_ComplexIntervalFieldElement(interval)
        self._complex_poly = poly.is_complex()
        self._interval = self.refine_interval(interval, 64)

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: type(v._descr)
            <class 'sage.rings.qqbar.ANRoot'>
            sage: loads(dumps(v)) == v
            True
        """
        return (ANRoot, (self._poly, self._interval, self._multiplicity))

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: v._descr._repr_()
            'Root 1.618033988749894849? of y^2 - y - 1'
        """
        return 'Root %s of %s' % (self._interval, self._poly)

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always ``True``
        for :class:`ANRoot`).

        EXAMPLES::

            sage: sage_input((AA(3)^(1/2))^(1/3), verify=True)
            # Verified
            R.<x> = AA[]
            AA.polynomial_root(AA.common_polynomial(x^3 - AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))), RIF(RR(1.2009369551760025), RR(1.2009369551760027)))

        These two examples are too big to verify quickly. (Verification
        would create a field of degree 28.)::

            sage: sage_input((sqrt(AA(3))^(5/7))^(9/4))
            R.<x> = AA[]
            v1 = AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))
            v2 = v1*v1
            v3 = AA.polynomial_root(AA.common_polynomial(x^7 - v2*v2*v1), RIF(RR(1.4804728524798112), RR(1.4804728524798114)))
            v4 = v3*v3
            v5 = v4*v4
            AA.polynomial_root(AA.common_polynomial(x^4 - v5*v5*v3), RIF(RR(2.4176921938267877), RR(2.4176921938267881)))
            sage: sage_input((sqrt(QQbar(-7))^(5/7))^(9/4))
            R.<x> = QQbar[]
            v1 = QQbar.polynomial_root(AA.common_polynomial(x^2 + 7), CIF(RIF(RR(0)), RIF(RR(2.6457513110645903), RR(2.6457513110645907))))
            v2 = v1*v1
            v3 = QQbar.polynomial_root(AA.common_polynomial(x^7 - v2*v2*v1), CIF(RIF(RR(0.8693488875796217), RR(0.86934888757962181)), RIF(RR(1.8052215661454434), RR(1.8052215661454436))))
            v4 = v3*v3
            v5 = v4*v4
            QQbar.polynomial_root(AA.common_polynomial(x^4 - v5*v5*v3), CIF(RIF(-RR(3.8954086044650791), -RR(3.8954086044650786)), RIF(RR(2.7639398015408925), RR(2.7639398015408929))))
            sage: x = polygen(QQ)
            sage: sage_input(AA.polynomial_root(x^2-x-1, RIF(1, 2)), verify=True)
            # Verified
            R.<y> = QQ[]
            AA.polynomial_root(AA.common_polynomial(y^2 - y - 1), RIF(RR(1.6180339887498947), RR(1.6180339887498949)))
            sage: sage_input(QQbar.polynomial_root(x^3-5, CIF(RIF(-3, 0), RIF(0, 3))), verify=True)
            # Verified
            R.<y> = QQ[]
            QQbar.polynomial_root(AA.common_polynomial(y^3 - 5), CIF(RIF(-RR(0.85498797333834853), -RR(0.85498797333834842)), RIF(RR(1.4808826096823642), RR(1.4808826096823644))))
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: rt = ANRoot(x^3 - 2, RIF(0, 4))
            sage: rt.handle_sage_input(sib, False, True)
            ({call: {getattr: {atomic:QQbar}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:y {constr_parent: {subscr: {atomic:QQ}[{atomic:'y'}]} with gens: ('y',)}} {atomic:3}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.259921049894873})}, {call: {atomic:RR}({atomic:1.2599210498948732})})})},
             True)
        """
        parent = sib.name('QQbar' if is_qqbar else 'AA')
        poly = sib(self._poly)
        intv = self._interval
        # Check whether a 53-bit interval actually isolates the root.
        # If so, use it, because 53-bit intervals print prettier.
        if is_ComplexIntervalFieldElement(intv):
            loose_intv = CIF(intv)
        else:
            loose_intv = RIF(intv)
        # If the derivative of the polynomial is bounded away from 0
        # over this interval, then it definitely isolates a root.
        if self._poly._poly.derivative()(loose_intv) != 0:
            good_intv = loose_intv
        else:
            good_intv = intv
        return (parent.polynomial_root(poly, sib(good_intv)), True)

    def is_complex(self):
        r"""
        Whether this is a root in `\overline{\QQ}` (rather than `\mathbf{A}`).
        Note that this may return True even if the root is actually real, as
        the second example shows; it does *not* trigger exact computation to
        see if the root is real.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.is_complex()
            False
            sage: (x^2 - x - 1).roots(ring=QQbar, multiplicities=False)[1]._descr.is_complex()
            True
        """
        return self._complex

    def conjugate(self, n):
        r"""
        Complex conjugate of this :class:`ANRoot` object.

        EXAMPLES::

            sage: a = (x^2 + 23).roots(ring=QQbar, multiplicities=False)[0]             # optional - sage.symbolic
            sage: b = a._descr                                                          # optional - sage.symbolic
            sage: type(b)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANRoot'>
            sage: c = b.conjugate(a); c                                                 # optional - sage.symbolic
            <sage.rings.qqbar.ANUnaryExpr object at ...>
            sage: c.exactify()                                                          # optional - sage.symbolic
            -2*a + 1 where a^2 - a + 6 = 0 and a in 0.50000000000000000? - 2.397915761656360?*I
        """
        if not self._complex:
            return self
        if not self._complex_poly:
            return ANRoot(self._poly, self._interval.conjugate(), self._multiplicity)

        return ANUnaryExpr(n, 'conjugate')

    def refine_interval(self, interval, prec):
        r"""
        Takes an interval which is assumed to enclose exactly one root
        of the polynomial (or, with multiplicity=`k`, exactly one root
        of the `k-1`-st derivative); and a precision, in bits.

        Tries to find a narrow interval enclosing the root using
        interval arithmetic of the given precision. (No particular
        number of resulting bits of precision is guaranteed.)

        Uses a combination of Newton's method (adapted for interval
        arithmetic) and bisection. The algorithm will converge very
        quickly if started with a sufficiently narrow interval.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(AA)
            sage: rt2 = ANRoot(x^2 - 2, RIF(0, 2))
            sage: rt2.refine_interval(RIF(0, 2), 75)
            1.4142135623730950488017?
        """
        if self._complex or self._complex_poly:
            v = self._complex_refine_interval(interval, prec)
            if self._complex:
                return v
            else:
                return v.real()
        else:
            return self._real_refine_interval(interval, prec)

    def _real_refine_interval(self, interval, prec):
        r"""
        Does the calculation for ``refine_interval``.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(AA)
            sage: rt2 = ANRoot(x^2 - 2, RIF(0, 2))
            sage: rt2.refine_interval(RIF(0, 2), 75) # indirect doctest
            1.4142135623730950488017?
        """
        # Do not throw away bits in the original interval; doing so might
        # invalidate it (include an extra root)

        field = RealIntervalField(max(prec, interval.prec()))
        interval = field(interval)
        if interval.is_exact():
            return interval

        p = self._poly.poly()

        dp = p.derivative()
        for i in range(self._multiplicity - 1):
            p = dp
            dp = p.derivative()

        zero = field.zero()

        poly_ring = field['x']

        # interval_p = poly_ring(p)
        if p.base_ring() is QQ:
            interval_p = p.change_ring(field)
        else:
            coeffs = [c._interval_fast(prec) for c in p.list()]
            interval_p = poly_ring(coeffs)

        # This special case is important: this is the only way we could
        # refine "infinitely deep" (we could get an interval of diameter
        # about 2^{-2^31}, and then hit floating-point underflow); avoiding
        # this case here means we do not have to worry about iterating too
        # many times later
        if interval_p[0].is_zero() and interval.contains_zero():
            return zero

        # interval_dp = poly_ring(dp)
        if p.base_ring() is QQ:
            interval_dp = p.derivative().change_ring(field)
        else:
            dcoeffs = [c.interval_fast(field) for c in dp.list()]
            interval_dp = poly_ring(dcoeffs)

        l = interval.lower()
        pl = interval_p(field(l))
        u = interval.upper()
        pu = interval_p(field(u))

        if pl.contains_zero():
            if pu.contains_zero():
                return interval
            else:
                su = pu.unique_sign()
                sl = -su
        elif pu.contains_zero():
            sl = pl.unique_sign()
            su = -sl
        else:
            sl = pl.unique_sign()
            su = pu.unique_sign()
            if sl == su:
                # Oops...
                raise ValueError("Refining interval that does not bound unique root!")

        while True:
            assert l == interval.lower()
            assert u == interval.upper()
            assert pl.contains_zero() or \
                   pu.contains_zero() or \
                   pl.unique_sign() != pu.unique_sign()

            # Use a simple algorithm:
            # Try an interval Newton-Raphson step. If this does not add at
            # least one bit of information, or if it fails (because the
            # slope is not bounded away from zero), then try bisection.
            # If this fails because the value at the midpoint is not
            # bounded away from zero, then also try the 1/4 and 3/4 points.
            # If all of these fail, then return the current interval.

            slope = interval_dp(interval)
            diam = interval.diameter()
            newton_lower = True

            if not slope.contains_zero():
                newton_lower = not newton_lower

                if newton_lower:
                    interval = interval.intersection(field(l) - pl/slope)
                else:
                    interval = interval.intersection(field(u) - pu/slope)
                new_diam = interval.diameter()

                if new_diam == 0:
                    # Wow, we managed to find the answer exactly.
                    # (I think this can only happen with a linear polynomial,
                    # in which case we shouldn't have been in this
                    # function at all; but oh well.)
                    return interval

                if (new_diam << 1) <= diam:
                    # We got at least one bit
                    l = interval.lower()
                    u = interval.upper()
                    pl = interval_p(field(l))
                    pu = interval_p(field(u))
                    continue

            # bisection
            for i,j in [(2,2),(3,1),(1,3)]:
                c = (i*l + j*u) / 4
                pc = interval_p(field(c))

                if c <= l or c >= u:
                    # not enough precision
                    return interval

                if pc.contains_zero():
                    continue
                elif pc.unique_sign() != sl:
                    interval = field(l, c)
                    u = c
                    pu = pc
                    break
                else:
                    interval = field(c, u)
                    l = c
                    pl = pc
                    break
            else:
                # bisection failed
                return interval

    def _complex_refine_interval(self, interval, prec):
        r"""
        Takes an interval which is assumed to enclose exactly one root
        of the polynomial (or, with multiplicity=`k`, exactly one root
        of the `k-1`-st derivative); and a precision, in bits.

        Tries to find a narrow interval enclosing the root using
        interval arithmetic of the given precision. (No particular
        number of resulting bits of precision is guaranteed.)

        Uses Newton's method (adapted for interval arithmetic). The
        algorithm will converge very quickly if started with a
        sufficiently narrow interval. If Newton's method fails, then
        we falls back on computing all the roots of the polynomial
        numerically, and select the appropriate root.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: intv = CIF(RIF(0, 1), RIF(0.1, 1))
            sage: rt = ANRoot(x^5 - 1, intv)
            sage: new_intv = rt.refine_interval(intv, 53); new_intv # indirect doctest
            0.3090169943749474241? + 0.951056516295153573?*I
            sage: rt.refine_interval(new_intv, 70)
            0.30901699437494742411? + 0.95105651629515357212?*I
        """
        # Don't throw away bits in the original interval; doing so might
        # invalidate it (include an extra root)

        field = ComplexIntervalField(max(prec, interval.prec()))
        interval = field(interval)
        if interval.is_exact():
            return interval

        p = self._poly.poly()
        dp = p.derivative()
        for i in range(self._multiplicity - 1):
            p = dp
            dp = p.derivative()

        zero = field.zero()

        poly_ring = field['x']

        # interval_p = poly_ring(p)
        if p.base_ring() is QQ:
            interval_p = p.change_ring(field)
        else:
            coeffs = [c.interval_fast(field) for c in p.list()]
            interval_p = poly_ring(coeffs)

        # This special case is important: this is the only way we could
        # refine "infinitely deep" (we could get an interval of diameter
        # about 2^{-2^31}, and then hit floating-point underflow); avoiding
        # this case here means we do not have to worry about iterating too
        # many times later
        if interval_p[0].is_zero() and zero in interval:
            return zero

        # interval_dp = poly_ring(dp)
        if p.base_ring() is QQ:
            interval_dp = p.derivative().change_ring(field)
        else:
            dcoeffs = [c.interval_fast(field) for c in dp.list()]
            interval_dp = poly_ring(dcoeffs)

        while True:
            center = field(interval.center())
            val = interval_p(center)

            slope = interval_dp(interval)

            diam = interval.diameter()

            if zero in slope:
                # Give up and fall back on root isolation.
                return self._complex_isolate_interval(interval, prec)

            if not (zero in slope):
                new_range = center - val / slope
                interval = interval.intersection(new_range)

                new_diam = interval.diameter()

                if new_diam == 0:
                    # Wow; we nailed it exactly. (This may happen
                    # whenever the root is exactly equal to some
                    # floating-point number, and cannot happen
                    # if the root is not equal to a floating-point
                    # number.)  We just return the perfect answer.
                    return interval

                if new_diam == diam:
                    # We're not getting any better. There are two
                    # possible reasons for this. Either we have
                    # refined as much as possible given the imprecision
                    # of our interval polynomial, and we have the best
                    # answer we're going to get at this precision;
                    # or we started with a poor approximation to the
                    # root, resulting in a broad range of possible
                    # slopes in this interval, and Newton-Raphson refining
                    # is not going to help.

                    # I do not have a formal proof, but I believe the
                    # following test differentiates between these two
                    # behaviors. (If I'm wrong, we might get bad behavior
                    # like infinite loops, but we still won't actually
                    # return a wrong answer.)

                    if val.contains_zero():
                        # OK, center must be a good approximation
                        # to the root (in the current precision, anyway).
                        # And the expression "center - val / slope"
                        # above means that we have a pretty good interval,
                        # even if slope is a poor estimate.
                        return interval

                    # The center of the current interval is known
                    # not to be a root. This should let us divide
                    # the interval in half, and improve on our previous
                    # estimates. I can only think of two reasons why
                    # it might not:
                    # 1) the "center" of the interval is actually
                    # on one of the edges of the interval (because the
                    # interval is only one ulp wide), or
                    # 2) the slope estimate is so bad that
                    # "center - val / slope" doesn't give us information.

                    # With complex intervals implemented as
                    # rectangular regions of the complex plane, it's
                    # possible that "val / slope" includes zero even
                    # if both "val" and "slope" are bounded away from
                    # zero, if the diameter of the (interval) argument
                    # of val or slope is large enough.

                    # So we test the diameter of the argument of
                    # slope; if it's small, we decide that we must
                    # have a good interval, but if it's big, we decide
                    # that we probably can't make progress with
                    # Newton-Raphson.

                    # I think the relevant measure of "small" is
                    # whether the diameter is less than pi/2; in that
                    # case, no matter the value of "val" (as long as
                    # "val" is fairly precise), "val / slope" should
                    # be bounded away from zero. But we compare
                    # against 1 instead, in the hopes that this might
                    # be slightly faster.

                    if slope.argument().absolute_diameter() < 1:
                        return interval

                    # And now it's time to give up.
                    return self._complex_isolate_interval(interval, prec)

    def _complex_isolate_interval(self, interval, prec):
        """
        Find a precise approximation to the unique root in interval,
        by finding a precise approximation to all roots of the
        polynomial, and checking which one is in interval. Slow but sure.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: intv = CIF(RIF(0, 1), RIF(0.1, 1))
            sage: rt = ANRoot(x^5 - 1, intv)
            sage: rt._complex_isolate_interval(intv, 53)
            0.3090169943749474241? + 0.951056516295153573?*I
        """
        rts = self._poly.complex_roots(prec, self._multiplicity)

        # Find all the roots that overlap interval.
        our_root = [rt for rt in rts if rt.overlaps(interval)]

        if len(our_root) == 1:
            return our_root[0]

        if not our_root:
            raise ValueError("Complex root interval does not include any roots")

        # We have more than one root that overlap the current interval.
        # Technically, this might not be an error; perhaps the actual
        # root is just outside our interval, even though the (presumably
        # tight) interval containing that root touches our interval.

        # But it seems far more likely that the provided interval is
        # just too big.

        raise ValueError("Complex root interval probably includes multiple roots")

    def exactify(self):
        """
        Return either an :class:`ANRational` or an
        :class:`ANExtensionElement` with the same value as this number.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: two = ANRoot((x-2)*(x-sqrt(QQbar(2))), RIF(1.9, 2.1))
            sage: two.exactify()
            2
            sage: strange = ANRoot(x^2 + sqrt(QQbar(3))*x - sqrt(QQbar(2)), RIF(-0, 1))
            sage: strange.exactify()
            a where a^8 - 6*a^6 + 5*a^4 - 12*a^2 + 4 = 0 and a in 0.6051012265139511?

        TESTS:

        Verify that :trac:`12727` is fixed::

            sage: m = sqrt(sin(pi/5)); a = QQbar(m); b = AA(m)                          # optional - sage.symbolic
            sage: a.minpoly()                                                           # optional - sage.symbolic
            x^8 - 5/4*x^4 + 5/16
            sage: b.minpoly()                                                           # optional - sage.symbolic
            x^8 - 5/4*x^4 + 5/16
        """
        gen = self._poly.generator()

        if gen.is_trivial():
            qpf = self._poly.factors()

            def find_fn(factor, prec):
                return factor(self._interval_fast(prec))
            my_factor = find_zero_result(find_fn, qpf)

            # Factoring always returns monic polynomials over the rationals
            assert(my_factor.is_monic())

            if my_factor.degree() == 1:
                return ANRational(-my_factor[0])

            den, my_factor = clear_denominators(my_factor)

            red_elt, red_back, red_pol = do_polred(my_factor)

            field = NumberField(red_pol, 'a', check=False)

            def intv_fn(rif):
                return conjugate_expand(red_elt(self._interval_fast(rif) * den))
            new_intv = conjugate_shrink(isolating_interval(intv_fn, red_pol))
            root = ANRoot(QQx(red_pol), new_intv)
            new_gen = AlgebraicGenerator(field, root)

            return ANExtensionElement(new_gen, red_back(field.gen()) / den)
        else:
            fpf = self._poly.factors()

            def find_fn(factor, prec):
                # XXX
                ifield = (ComplexIntervalField if self.is_complex() else RealIntervalField)(prec)
                if_poly = ifield['x']
                gen_val = gen._interval_fast(prec)
                self_val = self._interval_fast(prec)
                v = [c.polynomial()(gen_val) for c in factor]
                # This try/except can be triggered if ifield is Real
                # but the entries in v have some imaginary part that
                # is only known to be 0 to very low precision, e.g.,
                # as in Issue #12727.  In such cases, we instead create
                # the polynomial over the appropriate complex interval
                # field, which is mathematically safe, unlike taking
                # real parts would be.
                try:
                    ip = if_poly(v)
                except TypeError:
                    if_poly = ComplexIntervalField(prec)['x']
                    ip = if_poly(v)
                return ip(self_val)
            my_factor = find_zero_result(find_fn, fpf)

            assert(my_factor.is_monic())

            if my_factor.degree() == 1:
                return ANExtensionElement(gen, -my_factor[0])

            # rnfequation needs a monic polynomial with integral coefficients.
            # We achieve this with a change of variables.

            den, my_factor = clear_denominators(my_factor)

            pari_nf = gen.pari_field()

            x, y = QQxy.gens()
            my_factor = QQxy['z']([c.polynomial()(y) for c in my_factor])(x)

            # XXX much duplicate code with AlgebraicGenerator.union()

            # XXX need more caching here
            newpol, self_pol, k = pari_nf.rnfequation(my_factor, 1)
            k = int(k)

            newpol_sage = QQx(newpol)
            newpol_sage_y = QQy(newpol_sage)

            red_elt, red_back, red_pol = do_polred(newpol_sage_y)

            new_nf = NumberField(red_pol, name='a', check=False)

            self_pol_sage = QQx(self_pol.lift())

            def intv_fn(prec):
                return conjugate_expand(red_elt(gen._interval_fast(prec) * k + self._interval_fast(prec) * den))
            new_intv = conjugate_shrink(isolating_interval(intv_fn, red_pol))

            root = ANRoot(QQx(red_pol), new_intv)
            new_gen = AlgebraicGenerator(new_nf, root)
            red_back_a = red_back(new_nf.gen())
            new_poly = ((QQx_x - k * self_pol_sage)(red_back_a) / den)
            return ANExtensionElement(new_gen, new_poly)

    def _more_precision(self):
        """
        Recompute the interval enclosing this ``ANRoot`` object at higher
        precision.

        EXAMPLES::

            sage: x = polygen(QQ); y = (x^3 + x + 1).roots(AA,multiplicities=False)[0]
            sage: z = y._descr
            sage: z._interval.prec()
            64
            sage: z._more_precision()
            sage: z._interval.prec()
            128
        """
        prec = self._interval.prec()
        self._interval = self.refine_interval(self._interval, prec * 2)

    def _interval_fast(self, prec):
        """
        Given a RealIntervalField, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field. (More precision may be used
        in the computation.)

        EXAMPLES::

            sage: x = polygen(QQ); y = (x^3 + x + 1).roots(AA,multiplicities=False)[0]._descr
            sage: y._interval_fast(128)
            -0.68232780382801932736948373971104825689?

        Check that :trac:`15493` is fixed::

            sage: y._interval_fast(20).parent() is RealIntervalField(20)
            True
        """
        if prec == self._interval.prec():
            return self._interval
        if prec < self._interval.prec():
            return self._interval.parent().to_prec(prec)(self._interval)
        self._more_precision()
        return self._interval_fast(prec)


class ANExtensionElement(ANDescr):
    r"""
    The subclass of :class:`ANDescr` that represents a number field
    element in terms of a specific generator. Consists of a polynomial
    with rational coefficients in terms of the generator, and the
    generator itself, an :class:`AlgebraicGenerator`.
    """

    def __new__(self, generator, value):
        if value.is_rational():
            return ANRational(value._rational_())
        else:
            return ANDescr.__new__(self)

    def __init__(self, generator, value):
        self._generator = generator
        self._value = value
        self._exactly_real = not generator.is_complex()

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: v.exactify()
            sage: type(v._descr)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: loads(dumps(v)) == v
            True
        """
        return (ANExtensionElement, (self._generator, self._value))

    def _repr_(self):
        fgen = self._generator._field.gen()
        sgen = str(fgen)
        return '%s where %s = 0 and %s in %s' % (self._value,
                                                 self._generator.field().polynomial()._repr(name=sgen),
                                                 sgen,
                                                 self._generator._interval_fast(53))

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always ``True``
        for :class:`ANExtensionElement`).

        EXAMPLES::

            sage: I = QQbar(I)
            sage: sage_input(3+4*I, verify=True)
            # Verified
            QQbar(3 + 4*I)
            sage: v = QQbar.zeta(3) + QQbar.zeta(5)
            sage: v - v == 0
            True
            sage: sage_input(vector(QQbar, (4-3*I, QQbar.zeta(7))), verify=True)
            # Verified
            R.<y> = QQ[]
            vector(QQbar, [4 - 3*I, QQbar.polynomial_root(AA.common_polynomial(y^6 + y^5 + y^4 + y^3 + y^2 + y + 1), CIF(RIF(RR(0.62348980185873348), RR(0.62348980185873359)), RIF(RR(0.7818314824680298), RR(0.78183148246802991))))])
            sage: sage_input(v, verify=True)
            # Verified
            R.<y> = QQ[]
            v = QQbar.polynomial_root(AA.common_polynomial(y^8 - y^7 + y^5 - y^4 + y^3 - y + 1), CIF(RIF(RR(0.91354545764260087), RR(0.91354545764260098)), RIF(RR(0.40673664307580015), RR(0.40673664307580021))))
            v^5 + v^3
            sage: v = QQbar(sqrt(AA(2)))
            sage: v.exactify()
            sage: sage_input(v, verify=True)
            # Verified
            R.<y> = QQ[]
            QQbar(AA.polynomial_root(AA.common_polynomial(y^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))))
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: extel = ANExtensionElement(QQbar_I_generator, QQbar_I_generator.field().gen() + 1)
            sage: extel.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({binop:+ {atomic:1} {atomic:I}})}, True)
        """
        if self._generator is QQbar_I_generator:
            assert(is_qqbar)
            re, im = self._value.list()
            im_part = sib.prod([sib(im, True), sib.name('I')], simplify=True)
            v = sib.sum([sib(re, True), im_part], simplify=True)
            if coerce != 2:
                v = sib.name('QQbar')(v)
                return (v, True)
            return (v, False)

        result_is_qqbar = self._generator.is_complex()

        rt = sib(self._generator.root_as_algebraic())
        # For the best fidelity, we really ought to somehow ensure
        # that rt is exactified, but sage_input doesn't support that
        # nicely. Skip it for now.
        # The following is copied with slight mods from polynomial_element.pyx
        coeffs = [sib(c, True) for c in self._value.list()]
        terms = []
        for i in range(len(coeffs) - 1, -1, -1):
            if i > 0:
                if i > 1:
                    rt_pow = rt**sib.int(i)
                else:
                    rt_pow = rt
                terms.append(sib.prod((coeffs[i], rt_pow), simplify=True))
            else:
                terms.append(coeffs[i])
        v = sib.sum(terms, simplify=True)
        if result_is_qqbar != is_qqbar:
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)
        return (v, True)

    def is_complex(self):
        r"""
        Return ``True`` if the number field that defines this element is not real.

        This does not imply that the element itself is definitely non-real, as
        in the example below.

        EXAMPLES::

            sage: rt2 = QQbar(sqrt(2))                                                  # optional - sage.symbolic
            sage: rtm3 = QQbar(sqrt(-3))                                                # optional - sage.symbolic
            sage: x = rtm3 + rt2 - rtm3                                                 # optional - sage.symbolic
            sage: x.exactify()                                                          # optional - sage.symbolic
            sage: y = x._descr                                                          # optional - sage.symbolic
            sage: type(y)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: y.is_complex()                                                        # optional - sage.symbolic
            True
            sage: x.imag() == 0                                                         # optional - sage.symbolic
            True
        """
        return not self._exactly_real

    def is_simple(self):
        r"""
        Check whether this descriptor represents a value with the same
        algebraic degree as the number field associated with the descriptor.

        For :class:`ANExtensionElement` elements, we check this by
        comparing the degree of the minimal polynomial to the degree
        of the field.

        EXAMPLES::

            sage: rt2 = AA(sqrt(2))                                                     # optional - sage.symbolic
            sage: rt3 = AA(sqrt(3))                                                     # optional - sage.symbolic
            sage: rt2b = rt3 + rt2 - rt3                                                # optional - sage.symbolic
            sage: rt2.exactify()                                                        # optional - sage.symbolic
            sage: rt2._descr                                                            # optional - sage.symbolic
            a where a^2 - 2 = 0 and a in 1.414213562373095?
            sage: rt2._descr.is_simple()                                                # optional - sage.symbolic
            True

            sage: rt2b.exactify()                                                       # optional - sage.symbolic
            sage: rt2b._descr                                                           # optional - sage.symbolic
            a^3 - 3*a where a^4 - 4*a^2 + 1 = 0 and a in -0.5176380902050415?
            sage: rt2b._descr.is_simple()                                               # optional - sage.symbolic
            False
        """
        try:
            return self._is_simple
        except AttributeError:
            self._is_simple = (self.minpoly().degree() == self.generator().field().degree())
            return self._is_simple

    def generator(self):
        r"""
        Return the :class:`~AlgebraicGenerator` object corresponding to ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: v.generator()
            Number Field in a with defining polynomial y^2 - y - 1 with a in 1.618033988749895?
        """
        return self._generator

    def exactify(self):
        r"""
        Return an exact representation of ``self``.

        Since ``self`` is already exact, just return ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: type(v)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: v.exactify() is v
            True
        """
        return self

    def field_element_value(self):
        r"""
        Return the underlying number field element.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: v.field_element_value()
            a
        """
        return self._value

    def minpoly(self):
        """
        Compute the minimal polynomial of this algebraic number.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: type(v)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: v.minpoly()
            x^2 - x - 1
        """
        try:
            return self._minpoly
        except AttributeError:
            self._minpoly = self._value.minpoly()
            return self._minpoly

    def simplify(self, n):
        """
        Compute an exact representation for this descriptor, in the
        smallest possible number field.

        INPUT:

        - ``n`` -- The element of ``AA`` or ``QQbar`` corresponding
          to this descriptor.

        EXAMPLES::

            sage: rt2 = AA(sqrt(2))                                                     # optional - sage.symbolic
            sage: rt3 = AA(sqrt(3))                                                     # optional - sage.symbolic
            sage: rt2b = rt3 + rt2 - rt3                                                # optional - sage.symbolic
            sage: rt2b.exactify()                                                       # optional - sage.symbolic
            sage: rt2b._descr                                                           # optional - sage.symbolic
            a^3 - 3*a where a^4 - 4*a^2 + 1 = 0 and a in -0.5176380902050415?
            sage: rt2b._descr.simplify(rt2b)                                            # optional - sage.symbolic
            a where a^2 - 2 = 0 and a in 1.414213562373095?
        """

        if self.is_simple():
            return self

        # This is very inefficient...
        # for instance, the .exactify() call will try to factor poly,
        # even though we know that poly is irreducible
        poly = self.minpoly()
        intv = isolating_interval(n._interval_fast, poly)
        new_v = QQbar.polynomial_root(poly, intv)
        new_v.exactify()
        return new_v._descr

    def _interval_fast(self, prec):
        gen_val = self._generator._interval_fast(prec)
        v = self._value.polynomial()(gen_val)
        if self._exactly_real and is_ComplexIntervalFieldElement(v):
            return v.real()
        return v

    # for these three functions the argument n is not used (but it is there
    # anyway for compatibility)

    def neg(self, n):
        r"""
        Negation of ``self``.

        EXAMPLES::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))                                 # optional - sage.symbolic
            sage: a.exactify()                                                          # optional - sage.symbolic
            sage: b = a._descr                                                          # optional - sage.symbolic
            sage: type(b)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: c = b.neg(None); c  # random (not uniquely represented)               # optional - sage.symbolic
            -1/3*a^3 + 1/3*a^2 - a - 1 where a^4 - 2*a^3 + a^2 + 6*a + 3 = 0
             and a in 1.724744871391589? + 1.573132184970987?*I
            sage: (c.generator() == b.generator()                                       # optional - sage.symbolic
            ....:  and c.field_element_value() + b.field_element_value() == 0)
            True

        The parameter is ignored::

            sage: (b.neg("random").generator() == c.generator()                         # optional - sage.symbolic
            ....:  and b.neg("random").field_element_value() == c.field_element_value())
            True
        """
        return ANExtensionElement(self._generator, -self._value)

    def invert(self, n):
        r"""
        Reciprocal of ``self``.

        EXAMPLES::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))                                 # optional - sage.symbolic
            sage: a.exactify()                                                          # optional - sage.symbolic
            sage: b = a._descr                                                          # optional - sage.symbolic
            sage: type(b)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: c = b.invert(None); c  # random (not uniquely represented)            # optional - sage.symbolic
            -7/3*a^3 + 19/3*a^2 - 7*a - 9 where a^4 - 2*a^3 + a^2 + 6*a + 3 = 0
             and a in 1.724744871391589? + 1.573132184970987?*I
            sage: (c.generator() == b.generator()                                       # optional - sage.symbolic
            ....:  and c.field_element_value() * b.field_element_value() == 1)
            True

        The parameter is ignored::

            sage: (b.invert("random").generator() == c.generator()                      # optional - sage.symbolic
            ....:  and b.invert("random").field_element_value() == c.field_element_value())
            True
        """
        return ANExtensionElement(self._generator, ~self._value)

    def conjugate(self, n):
        r"""
        Complex conjugate of ``self``.

        EXAMPLES::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))                                 # optional - sage.symbolic
            sage: a.exactify()                                                          # optional - sage.symbolic
            sage: b = a._descr                                                          # optional - sage.symbolic
            sage: type(b)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: c = b.conjugate(None); c  # random (not uniquely represented)         # optional - sage.symbolic
            1/3*a^3 - 1/3*a^2 + a + 1 where a^4 - 2*a^3 + a^2 + 6*a + 3 = 0
             and a in 1.724744871391589? - 1.573132184970987?*I

        Internally, complex conjugation is implemented by taking the
        same abstract field element but conjugating the complex embedding of
        the field::

            sage: c.generator() == b.generator().conjugate()                            # optional - sage.symbolic
            True
            sage: c.field_element_value() == b.field_element_value()                    # optional - sage.symbolic
            True

        The parameter is ignored::

            sage: (b.conjugate("random").generator() == c.generator()                   # optional - sage.symbolic
            ....:  and b.conjugate("random").field_element_value() == c.field_element_value())
            True

        """
        if self._exactly_real:
            return self
        else:
            return ANExtensionElement(self._generator.conjugate(), self._value)

    # The rest of these unary operations do actually use n, which is an
    # AlgebraicNumber pointing to self.

    def norm(self, n):
        r"""
        Norm of ``self`` (square of complex absolute value)

        EXAMPLES::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))                                 # optional - sage.symbolic
            sage: a.exactify()                                                          # optional - sage.symbolic
            sage: b = a._descr                                                          # optional - sage.symbolic
            sage: type(b)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.norm(a)                                                             # optional - sage.symbolic
            <sage.rings.qqbar.ANUnaryExpr object at ...>
        """
        if self._exactly_real:
            return (n * n)._descr
        elif self._generator is QQbar_I_generator:
            return ANRational(self._value.norm())
        else:
            return ANUnaryExpr(n, 'norm')

    def abs(self, n):
        r"""
        Return the absolute value of ``self`` (square root of the norm).

        EXAMPLES::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))                                 # optional - sage.symbolic
            sage: a.exactify()                                                          # optional - sage.symbolic
            sage: b = a._descr                                                          # optional - sage.symbolic
            sage: type(b)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.abs(a)                                                              # optional - sage.symbolic
            Root 3.146264369941972342? of x^2 - 9.89897948556636?
        """
        return AlgebraicReal(self.norm(n)).sqrt()._descr

    def rational_argument(self, n):
        r"""
        If the argument of ``self`` is `2\pi` times some rational number in `[1/2,
        -1/2)`, return that rational; otherwise, return ``None``.

        EXAMPLES::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(3))                                  # optional - sage.symbolic
            sage: a.exactify()                                                          # optional - sage.symbolic
            sage: b = a._descr                                                          # optional - sage.symbolic
            sage: type(b)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.rational_argument(a) is None                                        # optional - sage.symbolic
            True
            sage: x = polygen(QQ)
            sage: a = (x^4 + 1).roots(QQbar, multiplicities=False)[0]
            sage: a.exactify()
            sage: b = a._descr
            sage: b.rational_argument(a)
            -3/8
        """
        # If the argument of self is 2*pi times some rational number a/b,
        # then self/abs(self) is a root of the b'th cyclotomic polynomial.
        # This implies that the algebraic degree of self is at least
        # phi(b). Working backward, we know that the algebraic degree
        # of self is at most the degree of the generator, so that gives
        # an upper bound on phi(b). According to
        # http://mathworld.wolfram.com/TotientFunction.html,
        # phi(b) >= sqrt(b) for b > 6; this gives us an upper bound on b.
        # We then check to see if this is possible; if so, we test
        # to see if it actually holds.

        if self._exactly_real:
            if n > 0:
                return QQ.zero()
            else:
                return QQ((1, 2))

        gen_degree = self._generator._field.degree()
        if gen_degree <= 2:
            max_b = 6
        else:
            max_b = gen_degree * gen_degree
        rat_arg_fl = n._interval_fast(128).argument() / RealIntervalField(128).pi() / 2
        rat_arg = rat_arg_fl.simplest_rational()
        if rat_arg.denominator() > max_b:
            return None
        n_exp = n ** rat_arg.denominator()
        if n_exp.real() > AA.zero() and n_exp.imag().is_zero():
            return rat_arg
        # Strictly speaking, we need to look for the second-simplest
        # rational in rat_arg_fl and make sure its denominator is > max_b.
        # For now, we just punt.
        raise NotImplementedError


class ANUnaryExpr(ANDescr):
    def __init__(self, arg, op):
        r"""
        Initialize this ANUnaryExpr.

        EXAMPLES::

            sage: t = ~QQbar(sqrt(2)); type(t._descr)  # indirect doctest               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        self._arg = arg
        self._op = op
        self._complex = True

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = ~QQbar(sqrt(2)); type(t._descr)                                   # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANUnaryExpr'>
            sage: loads(dumps(t)) == 1/QQbar(sqrt(2))                                   # optional - sage.symbolic
            True
        """
        return (ANUnaryExpr, (self._arg, self._op))

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always
        ``True`` for :class:`ANUnaryExpr`).

        EXAMPLES::

            sage: sage_input(-sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            -AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(~sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            ~AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(sqrt(QQbar(-3)).conjugate(), verify=True)
            # Verified
            R.<x> = QQbar[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 + 3), CIF(RIF(RR(0)), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))).conjugate()
            sage: sage_input(QQbar.zeta(3).real(), verify=True)
            # Verified
            R.<y> = QQ[]
            QQbar.polynomial_root(AA.common_polynomial(y^2 + y + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871)))).real()
            sage: sage_input(QQbar.zeta(3).imag(), verify=True)
            # Verified
            R.<y> = QQ[]
            QQbar.polynomial_root(AA.common_polynomial(y^2 + y + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871)))).imag()
            sage: sage_input(abs(sqrt(QQbar(-3))), verify=True)
            # Verified
            R.<x> = QQbar[]
            abs(QQbar.polynomial_root(AA.common_polynomial(x^2 + 3), CIF(RIF(RR(0)), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))))
            sage: sage_input(sqrt(QQbar(-3)).norm(), verify=True)
            # Verified
            R.<x> = QQbar[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 + 3), CIF(RIF(RR(0)), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))).norm()
            sage: sage_input(QQbar(QQbar.zeta(3).real()), verify=True)
            # Verified
            R.<y> = QQ[]
            QQbar(QQbar.polynomial_root(AA.common_polynomial(y^2 + y + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871)))).real())
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: unexp = ANUnaryExpr(sqrt(AA(2)), '~')
            sage: unexp.handle_sage_input(sib, False, False)
            ({unop:~ {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.4142135623730949})}, {call: {atomic:RR}({atomic:1.4142135623730951})})})}},
             True)
            sage: unexp.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({unop:~ {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.4142135623730949})}, {call: {atomic:RR}({atomic:1.4142135623730951})})})}})},
             True)
        """
        arg_is_qqbar = self._arg.parent() is QQbar
        v = sib(self._arg)
        op = self._op
        if op == '-':
            v = -v
        elif op == '~':
            v = ~v
        elif op == 'conjugate':
            v = v.conjugate()
        elif op == 'real':
            v = v.real()
        elif op == 'imag':
            v = v.imag()
        elif op == 'abs':
            v = abs(v)
        elif op == 'norm':
            v = v.norm()
        else:
            raise NotImplementedError

        result_is_qqbar = arg_is_qqbar
        if op in ('real', 'imag', 'abs', 'norm'):
            result_is_qqbar = False
        if result_is_qqbar != is_qqbar:
            # The following version is not safe with respect to caching;
            # with the current sage_input.py, anything that gets entered
            # into the cache must be safe at all coercion levels.
            #             if is_qqbar and not coerce:
            #                 v = sib.name('QQbar')(v)
            #             if not is_qqbar and coerce != 2:
            #                 v = sib.name('AA')(v)
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)

        return (v, True)

    def is_complex(self):
        r"""
        Return whether or not this element is complex. Note that this is a data
        type check, and triggers no computations -- if it returns ``False``, the
        element might still be real, it just doesn't know it yet.

        EXAMPLES::

            sage: t = AA(sqrt(2))                                                       # optional - sage.symbolic
            sage: s = (-t)._descr                                                       # optional - sage.symbolic
            sage: s                                                                     # optional - sage.symbolic
            <sage.rings.qqbar.ANUnaryExpr object at ...>
            sage: s.is_complex()                                                        # optional - sage.symbolic
            False
            sage: QQbar(-sqrt(2))._descr.is_complex()                                   # optional - sage.symbolic
            True
        """
        return self._complex

    def _interval_fast(self, prec):
        r"""
        Calculate an approximation to this ``ANUnaryExpr`` object in an interval field of precision ``prec``.

        EXAMPLES::

            sage: t = AA(sqrt(2))                                                       # optional - sage.symbolic
            sage: s = (-t)._descr                                                       # optional - sage.symbolic
            sage: s                                                                     # optional - sage.symbolic
            <sage.rings.qqbar.ANUnaryExpr object at ...>
            sage: s._interval_fast(150)                                                 # optional - sage.symbolic
            -1.414213562373095048801688724209698078569671876?
        """
        op = self._op

        v = self._arg._interval_fast(prec)

        if not is_ComplexIntervalFieldElement(v):
            self._complex = False

        if op == '-':
            return -v

        if op == '~':
            return ~v

        if op == 'conjugate':
            if is_ComplexIntervalFieldElement(v):
                return v.conjugate()
            else:
                return v

        self._complex = False

        if op == 'real':
            if is_ComplexIntervalFieldElement(v):
                return v.real()
            else:
                return v

        if op == 'imag':
            if is_ComplexIntervalFieldElement(v):
                return v.imag()
            else:
                return RealIntervalField(prec)(0)

        if op == 'abs':
            return abs(v)

        if op == 'norm':
            if is_ComplexIntervalFieldElement(v):
                return v.norm()
            else:
                return v.square()

        raise NotImplementedError

    def exactify(self):
        r"""
        Trigger exact computation of ``self``.

        EXAMPLES::

            sage: v = (-QQbar(sqrt(2)))._descr                                          # optional - sage.symbolic
            sage: type(v)                                                               # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANUnaryExpr'>
            sage: v.exactify()                                                          # optional - sage.symbolic
            -a where a^2 - 2 = 0 and a in 1.414213562373095?
        """
        op = self._op
        arg = self._arg

        if op == '-':
            arg.exactify()
            return arg._descr.neg(None)

        if op == '~':
            arg.exactify()
            return arg._descr.invert(None)

        if op == 'real':
            arg.exactify()
            rv = (arg + arg.conjugate()) / 2
            rv.exactify()
            rvd = rv._descr
            rvd._exactly_real = True
            return rvd

        if op == 'imag':
            arg.exactify()
            iv = QQbar_I * (arg.conjugate() - arg) / 2
            iv.exactify()
            ivd = iv._descr
            ivd._exactly_real = True
            return ivd

        if op == 'abs':
            arg.exactify()
            if arg.parent() is AA:
                if arg.sign() > 0:
                    return arg._descr
                else:
                    return arg._descr.neg(None)

            v = (arg * arg.conjugate()).sqrt()
            v.exactify()
            vd = v._descr
            vd._exactly_real = True
            return vd

        if op == 'norm':
            arg.exactify()
            v = arg * arg.conjugate()
            v.exactify()
            vd = v._descr
            vd._exactly_real = True
            return vd

        if op == 'conjugate':
            arg.exactify()
            return arg._descr.conjugate(None)


class ANBinaryExpr(ANDescr):
    def __init__(self, left, right, op):
        r"""
        Initialize this ANBinaryExpr.

        EXAMPLES::

            sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3)); type(t._descr) # indirect doctest    # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANBinaryExpr'>
        """
        self._left = left
        self._right = right
        self._op = op
        self._complex = True

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3)); type(t._descr)                   # optional - sage.symbolic
            <class 'sage.rings.qqbar.ANBinaryExpr'>
            sage: loads(dumps(t)) == QQbar(sqrt(2)) + QQbar(sqrt(3))                    # optional - sage.symbolic
            True
        """
        return (ANBinaryExpr, (self._left, self._right, self._op))

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always
        ``True`` for :class:`ANBinaryExpr`).

        EXAMPLES::

            sage: sage_input(2 + sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2 + AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(sqrt(AA(2)) + 2, verify=True)
            # Verified
            R.<x> = AA[]
            AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))) + 2
            sage: sage_input(2 - sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2 - AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2 / sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2/AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2 + (-1*sqrt(AA(2))), verify=True)
            # Verified
            R.<x> = AA[]
            2 - AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2*sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2*AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: rt2 = sqrt(AA(2))
            sage: one = rt2/rt2
            sage: n = one+3
            sage: sage_input(n)
            R.<x> = AA[]
            v = AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            v/v + 3
            sage: one == 1
            True
            sage: sage_input(n)
            1 + AA(3)
            sage: rt3 = QQbar(sqrt(3))                                                  # optional - sage.symbolic
            sage: one = rt3/rt3                                                         # optional - sage.symbolic
            sage: n = sqrt(AA(2)) + one
            sage: one == 1                                                              # optional - sage.symbolic
            True
            sage: sage_input(n)
            R.<x> = AA[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))) + 1
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: binexp = ANBinaryExpr(AA(3), AA(5), operator.mul)
            sage: binexp.handle_sage_input(sib, False, False)
            ({binop:* {atomic:3} {call: {atomic:AA}({atomic:5})}}, True)
            sage: binexp.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({binop:* {atomic:3} {call: {atomic:AA}({atomic:5})}})}, True)
        """
        arg1 = self._left
        arg2 = self._right
        op = self._op

        # We want 2+QQbar.zeta(3) and QQbar.zeta(3)+2, not
        # QQbar(2)+QQbar.zeta(3). So we want to pass coerced=True to
        # an argument if it is rational (but if both arguments are
        # rational, we only want to set it for one of them).

        arg1_coerced = False
        arg2_coerced = False

        if isinstance(arg1._descr, ANRational):
            arg1_coerced = True
        elif isinstance(arg2._descr, ANRational):
            arg2_coerced = True

        arg1_is_qqbar = arg1.parent() is QQbar
        arg2_is_qqbar = arg2.parent() is QQbar

        result_is_qqbar = \
            (arg1_is_qqbar and not arg1_coerced) or \
            (arg2_is_qqbar and not arg2_coerced)

        v1 = sib(arg1, arg1_coerced)
        v2 = sib(arg2, arg2_coerced)

        if op is operator.add:
            v = sib.sum([v1, v2], simplify=True)
        elif op is operator.sub:
            v = sib.sum([v1, -v2], simplify=True)
        elif op is operator.mul:
            v = sib.prod([v1, v2], simplify=True)
        elif op is operator.truediv:
            v = v1 / v2
        else:
            raise RuntimeError("op is {}".format(op))

        if result_is_qqbar != is_qqbar:
            # The following version is not safe with respect to caching;
            # with the current sage_input.py, anything that gets entered
            # into the cache must be safe at all coercion levels.
            #             if is_qqbar and not coerce:
            #                 v = sib.name('QQbar')(v)
            #             if not is_qqbar and coerce != 2:
            #                 v = sib.name('AA')(v)
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)

        return (v, True)

    def is_complex(self):
        r"""
        Whether this element is complex. Does not trigger exact computation, so
        may return ``True`` even if the element is real.

        EXAMPLES::

            sage: x = (QQbar(sqrt(-2)) / QQbar(sqrt(-5)))._descr                        # optional - sage.symbolic
            sage: x.is_complex()                                                        # optional - sage.symbolic
            True
        """
        return self._complex

    def _interval_fast(self, prec):
        r"""
        Calculate an approximation to self in an interval field of precision prec.

        EXAMPLES::

            sage: x = (QQbar(sqrt(-2)) / QQbar(sqrt(-5)))._descr                        # optional - sage.symbolic
            sage: y= x._interval_fast(64); y                                            # optional - sage.symbolic
            0.632455532033675867?
            sage: y.parent()                                                            # optional - sage.symbolic
            Complex Interval Field with 64 bits of precision
        """
        op = self._op

        lv = self._left._interval_fast(prec)
        rv = self._right._interval_fast(prec)

        if not (is_ComplexIntervalFieldElement(lv) or is_ComplexIntervalFieldElement(rv)):
            self._complex = False

        return op(lv, rv)

    def exactify(self):
        """
        TESTS::

            sage: rt2c = QQbar.zeta(3) + AA(sqrt(2)) - QQbar.zeta(3)                    # optional - sage.symbolic
            sage: rt2c.exactify()                                                       # optional - sage.symbolic

        We check to make sure that this method still works even. We
        do this by increasing the recursion level at each step and
        decrease it before we return.
        We lower the recursion limit for this test to allow
        a test in reasonable time::

            sage: import sys
            sage: old_recursion_limit = sys.getrecursionlimit()
            sage: sys.setrecursionlimit(1000)

            sage: sys.getrecursionlimit()
            1000
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a=s([3,2]).expand(8)(flatten([[QQbar.zeta(3)^d for d in range(3)], [QQbar.zeta(5)^d for d in range(5)]]))
            sage: a.exactify(); a # long time
            0
            sage: sys.getrecursionlimit()
            1000

            sage: sys.setrecursionlimit(old_recursion_limit)
        """
        import sys
        with increase_recursion_limit(10):
            left = self._left
            right = self._right
            left.exactify()
            right.exactify()
            gen = left._exact_field().union(right._exact_field())
            left_value = gen(left._exact_value())
            right_value = gen(right._exact_value())

            value = self._op(left_value, right_value)

            if gen.is_trivial():
                return ANRational(value)
            else:
                return ANExtensionElement(gen, value)


# These are the functions used to add, subtract, multiply, and divide
# algebraic numbers. Basically, we try to compute exactly if both
# arguments are already known to be in the same number field. Otherwise
# we fall back to floating-point computation to be backed up by exact
# symbolic computation only as required.


def an_binop_rational(a, b, op):
    r"""
    Used to add, subtract, multiply or divide algebraic numbers.

    Used when both are actually rational.

    EXAMPLES::

        sage: from sage.rings.qqbar import an_binop_rational
        sage: f = an_binop_rational(QQbar(2), QQbar(3/7), operator.add)
        sage: f
        17/7
        sage: type(f)
        <class 'sage.rings.qqbar.ANRational'>

        sage: f = an_binop_rational(QQbar(2), QQbar(3/7), operator.mul)
        sage: f
        6/7
        sage: type(f)
        <class 'sage.rings.qqbar.ANRational'>
    """
    return ANRational(op(a._descr._value, b._descr._value))


def an_binop_expr(a, b, op):
    r"""
    Add, subtract, multiply or divide algebraic numbers represented as
    binary expressions.

    INPUT:

    - ``a``, ``b`` -- two elements

    - ``op`` -- an operator

    EXAMPLES::

        sage: a = QQbar(sqrt(2)) + QQbar(sqrt(3))                                       # optional - sage.symbolic
        sage: b = QQbar(sqrt(3)) + QQbar(sqrt(5))                                       # optional - sage.symbolic
        sage: type(a._descr); type(b._descr)                                            # optional - sage.symbolic
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: from sage.rings.qqbar import an_binop_expr
        sage: x = an_binop_expr(a, b, operator.add); x                                  # optional - sage.symbolic
        <sage.rings.qqbar.ANBinaryExpr object at ...>
        sage: x.exactify()                                                              # optional - sage.symbolic
        6/7*a^7 - 2/7*a^6 - 71/7*a^5 + 26/7*a^4 + 125/7*a^3 - 72/7*a^2 - 43/7*a + 47/7
        where a^8 - 12*a^6 + 23*a^4 - 12*a^2 + 1 = 0 and a in -0.3199179336182997?

        sage: a = QQbar(sqrt(2)) + QQbar(sqrt(3))                                       # optional - sage.symbolic
        sage: b = QQbar(sqrt(3)) + QQbar(sqrt(5))                                       # optional - sage.symbolic
        sage: type(a._descr)                                                            # optional - sage.symbolic
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: x = an_binop_expr(a, b, operator.mul); x                                  # optional - sage.symbolic
        <sage.rings.qqbar.ANBinaryExpr object at ...>
        sage: x.exactify()                                                              # optional - sage.symbolic
        2*a^7 - a^6 - 24*a^5 + 12*a^4 + 46*a^3 - 22*a^2 - 22*a + 9
        where a^8 - 12*a^6 + 23*a^4 - 12*a^2 + 1 = 0 and a in -0.3199179336182997?
    """
    return ANBinaryExpr(a, b, op)


def an_binop_element(a, b, op):
    r"""
    Add, subtract, multiply or divide two elements represented as elements of
    number fields.

    EXAMPLES::

        sage: sqrt2 = QQbar(2).sqrt()
        sage: sqrt3 = QQbar(3).sqrt()
        sage: sqrt5 = QQbar(5).sqrt()

        sage: a = sqrt2 + sqrt3; a.exactify()
        sage: b = sqrt3 + sqrt5; b.exactify()
        sage: type(a._descr)
        <class 'sage.rings.qqbar.ANExtensionElement'>
        sage: from sage.rings.qqbar import an_binop_element
        sage: an_binop_element(a, b, operator.add)
        <sage.rings.qqbar.ANBinaryExpr object at ...>
        sage: an_binop_element(a, b, operator.sub)
        <sage.rings.qqbar.ANBinaryExpr object at ...>
        sage: an_binop_element(a, b, operator.mul)
        <sage.rings.qqbar.ANBinaryExpr object at ...>
        sage: an_binop_element(a, b, operator.truediv)
        <sage.rings.qqbar.ANBinaryExpr object at ...>

    The code tries to use existing unions of number fields::

        sage: sqrt17 = QQbar(17).sqrt()
        sage: sqrt19 = QQbar(19).sqrt()
        sage: a = sqrt17 + sqrt19
        sage: b = sqrt17 * sqrt19 - sqrt17 + sqrt19 * (sqrt17 + 2)
        sage: b, type(b._descr)
        (40.53909377268655?, <class 'sage.rings.qqbar.ANBinaryExpr'>)
        sage: a.exactify()
        sage: b = sqrt17 * sqrt19 - sqrt17 + sqrt19 * (sqrt17 + 2)
        sage: b, type(b._descr)
        (40.53909377268655?, <class 'sage.rings.qqbar.ANExtensionElement'>)
    """
    ad = a._descr
    bd = b._descr
    adg = ad.generator()
    bdg = bd.generator()
    if adg == qq_generator or adg == bdg:
        return ANExtensionElement(bdg, op(ad._value, bd._value))

    if bdg == qq_generator:
        return ANExtensionElement(adg, op(ad._value, bd._value))

    if adg in bdg._unions or bdg in adg._unions:
        p = bdg._unions[adg] if adg in bdg._unions else adg._unions[bdg]
        p = p.parent
        adg2 = adg.super_poly(p)
        bdg2 = bdg.super_poly(p)
        av = ad._value.polynomial()(adg2)
        bv = bd._value.polynomial()(bdg2)
        return ANExtensionElement(p, op(av, bv))

    adg2 = adg.super_poly(bdg)
    if adg2 is not None:
        av = ad._value.polynomial()(adg2)
        return ANExtensionElement(bdg, op(av, bd._value))

    bdg2 = bdg.super_poly(adg)
    if bdg2 is not None:
        bv = bd._value.polynomial()(bdg2)
        return ANExtensionElement(adg, op(ad._value, bv))

    return ANBinaryExpr(a, b, op)


# instanciation of the multimethod dispatch
_binop_algo[ANRational, ANRational] = an_binop_rational
_binop_algo[ANRational, ANExtensionElement] = \
_binop_algo[ANExtensionElement, ANRational] = \
_binop_algo[ANExtensionElement, ANExtensionElement] = an_binop_element
for t1 in (ANRational, ANRoot, ANExtensionElement, ANUnaryExpr, ANBinaryExpr):
    for t2 in (ANUnaryExpr, ANBinaryExpr, ANRoot):
        _binop_algo[t1, t2] = _binop_algo[t2, t1] = an_binop_expr

qq_generator = AlgebraicGenerator(QQ, ANRoot(AAPoly.gen() - 1, RIF(1)))


def _init_qqbar():
    """
    This code indirectly uses a huge amount of sage, despite the fact
    that qqbar is imported rather early on in the sage loading. This function
    is called at the end of sage.all.

    EXAMPLES::

        sage: sage.rings.qqbar.QQbar_I_generator # indirect doctest
        Number Field in I with defining polynomial x^2 + 1 with I = 1*I with a in 1*I
    """
    global ZZX_x, AA_0, QQbar_I, AA_hash_offset, QQbar_hash_offset, QQbar_I_generator, QQbar_I_nf
    global QQ_0, QQ_1, QQ_1_2, QQ_1_4, RR_1_10

    RR_1_10 = RR(1) / 10
    QQ_0 = QQ.zero()
    QQ_1 = QQ.one()
    QQ_1_2 = QQ((1, 2))
    QQ_1_4 = QQ((1, 4))

    AA_0 = AA.zero()

    QQbar_I_nf = GaussianField()
    QQbar_I_generator = AlgebraicGenerator(QQbar_I_nf, ANRoot(AAPoly.gen()**2 + 1, CIF(0, 1)))
    QQbar_I = AlgebraicNumber(ANExtensionElement(QQbar_I_generator, QQbar_I_nf.gen()))

    AA_hash_offset = AA(~ZZ(123456789))

    QQbar_hash_offset = AlgebraicNumber(ANExtensionElement(QQbar_I_generator, ~ZZ(123456789) + QQbar_I_nf.gen() / ZZ(987654321)))

    ZZX_x = ZZ['x'].gen()


# This is used in the _algebraic_ method of the golden_ratio constant,
# in sage/symbolic/constants.py
AA_golden_ratio = None


def get_AA_golden_ratio():
    r"""
    Return the golden ratio as an element of the algebraic real field. Used by
    :meth:`sage.symbolic.constants.golden_ratio._algebraic_`.

    EXAMPLES::

        sage: AA(golden_ratio)  # indirect doctest                                      # optional - sage.symbolic
        1.618033988749895?
    """
    global AA_golden_ratio
    if AA_golden_ratio is None:
        AA_golden_ratio_nf = NumberField(ZZX_x**2 - ZZX_x - 1, 'phi')
        AA_golden_ratio_generator = AlgebraicGenerator(AA_golden_ratio_nf, ANRoot(AAPoly.gen()**2 - AAPoly.gen() - 1, RIF(1.618, 1.6181)))
        AA_golden_ratio = AlgebraicReal(ANExtensionElement(AA_golden_ratio_generator, AA_golden_ratio_nf.gen()))
    return AA_golden_ratio


# Support Python's numbers abstract base class
import numbers
numbers.Real.register(AlgebraicReal)
numbers.Complex.register(AlgebraicNumber)
