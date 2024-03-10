.. index:: elliptic curves

***************
Elliptic curves
***************

Conductor
=========
How do you compute the conductor of an elliptic curve (over
:math:`\QQ`) in Sage?

Once you define an elliptic curve :math:`E` in Sage, using the
``EllipticCurve`` command, the conductor is one of several "methods"
associated to :math:`E`. Here is an example of the syntax
(borrowed from section 2.4 "Modular forms" in the tutorial):

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E.conductor()
    10351

:math:`j`-invariant
=====================
How do you compute the :math:`j`-invariant of an elliptic curve
in Sage?

Other methods associated to the ``EllipticCurve`` class are
``j_invariant``, ``discriminant``, and ``weierstrass_model``. Here is
an example of their syntax.

::

    sage: E = EllipticCurve([0, -1, 1, -10, -20])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
    sage: E.j_invariant()
    -122023936/161051
    sage: E.short_weierstrass_model()
    Elliptic Curve defined by y^2  = x^3 - 13392*x - 1080432 over Rational Field
    sage: E.discriminant()
    -161051
    sage: E = EllipticCurve(GF(5),[0, -1, 1, -10, -20])
    sage: E.short_weierstrass_model()
    Elliptic Curve defined by y^2  = x^3 + 3*x + 3 over Finite Field of size 5
    sage: E.j_invariant()
    4

.. index:: elliptic curves

The :math:`GF(q)`-rational points on E
========================================

How do you compute the number of points of an elliptic curve over a
finite field?

Given an elliptic curve defined over :math:`\mathbb{F} = GF(q)`, Sage
can compute its set of :math:`\mathbb{F}`-rational points

::

    sage: E = EllipticCurve(GF(5),[0, -1, 1, -10, -20])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + 4*x^2 over Finite Field of size 5
    sage: E.points()
    [(0 : 1 : 0), (0 : 0 : 1), (0 : 4 : 1), (1 : 0 : 1), (1 : 4 : 1)]
    sage: E.cardinality()
    5
    sage: G = E.abelian_group()
    sage: G
    Additive abelian group isomorphic to Z/5 embedded in Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + 4*x^2 over Finite Field of size 5
    sage: G.permutation_group()
    Permutation Group with generators [(1,2,3,4,5)]

.. index::
   pair: modular form; elliptic curve

Modular form associated to an elliptic curve over :math:`\QQ`
========================================================================

Let :math:`E` be a "nice" elliptic curve whose equation has
integer coefficients, let :math:`N` be the conductor of
:math:`E` and, for each :math:`n`, let :math:`a_n` be the
number appearing in the Hasse-Weil :math:`L`-function of
:math:`E`. The Taniyama-Shimura conjecture (proven by Wiles)
states that there exists a modular form of weight two and level
:math:`N` which is an eigenform under the Hecke operators and has
a Fourier series :math:`\sum_{n = 0}^\infty a_n q^n`. Sage can
compute the sequence :math:`a_n` associated to :math:`E`. Here
is an example.

::

    sage: E = EllipticCurve([0, -1, 1, -10, -20])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
    sage: E.conductor()
    11
    sage: E.anlist(20)
    [0, 1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2]
    sage: E.analytic_rank()
    0
