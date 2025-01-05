Some More Advanced Mathematics
==============================

Algebraic Geometry
------------------

You can define arbitrary algebraic varieties in Sage, but sometimes
nontrivial functionality is limited to rings over :math:`\QQ` or
finite fields. For example, we compute the union of two affine
plane curves, then recover the curves as the irreducible components
of the union.

::

    sage: x, y = AffineSpace(2, QQ, 'xy').gens()
    sage: C2 = Curve(x^2 + y^2 - 1)
    sage: C3 = Curve(x^3 + y^3 - 1)
    sage: D = C2 + C3
    sage: D
    Affine Plane Curve over Rational Field defined by
       x^5 + x^3*y^2 + x^2*y^3 + y^5 - x^3 - y^3 - x^2 - y^2 + 1
    sage: D.irreducible_components()
    [Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x^2 + y^2 - 1,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x^3 + y^3 - 1]

We can also find all points of intersection of the two curves by
intersecting them and computing the irreducible components.

.. link

::

    sage: V = C2.intersection(C3)
    sage: V.irreducible_components()
    [Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
    y,
    x - 1,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
    y - 1,
    x,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
    x + y + 2,
    2*y^2 + 4*y + 3]
<<<<<<< HEAD
    [
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y - 1,
      x,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y,
      x - 1,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x + y + 2,
      2*y^2 + 4*y + 3
    ]
||||||| c9dd1e85118
    [
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y,
      x - 1,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y - 1,
      x,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x + y + 2,
      2*y^2 + 4*y + 3
    ]
=======
    [Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       y,
       x - 1,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       y - 1,
       x,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x + y + 2,
       2*y^2 + 4*y + 3]
>>>>>>> develop

Thus, e.g., :math:`(1,0)` and :math:`(0,1)` are on both curves
(visibly clear), as are certain (quadratic) points whose
:math:`y` coordinates satisfy :math:`2y^2 + 4y + 3=0`.

Sage can compute the toric ideal of the twisted cubic in projective 3
space:

::

    sage: R.<a,b,c,d> = PolynomialRing(QQ, 4)
    sage: I = ideal(b^2-a*c, c^2-b*d, a*d-b*c)
    sage: F = I.groebner_fan(); F
    Groebner fan of the ideal:
    Ideal (b^2 - a*c, c^2 - b*d, -b*c + a*d) of Multivariate Polynomial Ring
    in a, b, c, d over Rational Field
    sage: F.reduced_groebner_bases ()
    [[-c^2 + b*d, -b*c + a*d, -b^2 + a*c],
     [-b*c + a*d, -c^2 + b*d, b^2 - a*c],
     [-c^3 + a*d^2, -c^2 + b*d, b*c - a*d, b^2 - a*c],
     [-c^2 + b*d, b^2 - a*c, b*c - a*d, c^3 - a*d^2],
     [-b*c + a*d, -b^2 + a*c, c^2 - b*d],
     [-b^3 + a^2*d, -b^2 + a*c, c^2 - b*d, b*c - a*d],
     [-b^2 + a*c, c^2 - b*d, b*c - a*d, b^3 - a^2*d],
     [c^2 - b*d, b*c - a*d, b^2 - a*c]]
    sage: F.polyhedralfan()
    Polyhedral fan in 4 dimensions of dimension 4

Elliptic Curves
---------------

Elliptic curve functionality includes most of the elliptic curve
functionality of PARI, access to the data in Cremona's online
tables (this requires an optional database package), the
functionality of mwrank, i.e., 2-descents with computation of the
full Mordell-Weil group, the SEA algorithm, computation of all
isogenies, much new code for curves over :math:`\QQ`, and some of Denis
Simon's algebraic descent software.

The command ``EllipticCurve`` for creating an elliptic curve has many
forms:


-  EllipticCurve([:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6`]):
   Returns the elliptic curve

   .. math::  y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6,


   where the :math:`a_i`'s are coerced into the parent of
   :math:`a_1`. If all the :math:`a_i` have parent :math:`\ZZ`, they are
   coerced into :math:`\QQ`.

-  EllipticCurve([:math:`a_4`, :math:`a_6`]): Same as above, but
   :math:`a_1=a_2=a_3=0`.

-  EllipticCurve(label): Returns the elliptic curve over from the
   Cremona database with the given (new!) Cremona label. The label is
   a string, such as ``"11a"`` or ``"37b2"``. The letter must be lower
   case (to distinguish it from the old labeling).

-  EllipticCurve(j): Returns an elliptic curve with
   :math:`j`-invariant :math:`j`.

-  EllipticCurve(R,
   [:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6`]):
   Create the elliptic curve over a ring :math:`R` with given
   :math:`a_i`'s as above.


We illustrate each of these constructors:

::

    sage: EllipticCurve([0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve([GF(5)(0),0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    sage: EllipticCurve([1,2])
    Elliptic Curve defined by y^2  = x^3 + x + 2 over Rational Field

    sage: EllipticCurve('37a')
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve_from_j(1)
    Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455 over Rational Field

    sage: EllipticCurve(GF(5), [0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

The pair :math:`(0,0)` is a point on the elliptic curve
:math:`E` defined by :math:`y^2 +
y = x^3 - x`. To create this
point in Sage type ``E([0,0])``. Sage can add points on such an
elliptic curve (recall elliptic curves support an additive group
structure where the point at infinity is the zero element and three
co-linear points on the curve add to zero):

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
    sage: P = E([0,0])
    sage: P + P
    (1 : 0 : 1)
    sage: 10*P
    (161/16 : -2065/64 : 1)
    sage: 20*P
    (683916417/264517696 : -18784454671297/4302115807744 : 1)
    sage: E.conductor()
    37

The elliptic curves over the complex numbers are parameterized by
the :math:`j`-invariant. Sage computes :math:`j`-invariant as
follows:

::

    sage: E = EllipticCurve([0,0,0,-4,2]); E
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: E.conductor()
    2368
    sage: E.j_invariant()
    110592/37

If we make a curve with the same :math:`j`-invariant as that of
:math:`E`, it need not be isomorphic to :math:`E`. In the
following example, the curves are not isomorphic because their
conductors are different.

::

    sage: F = EllipticCurve_from_j(110592/37)
    sage: F.conductor()
    37

However, the twist of :math:`F` by 2 gives an isomorphic curve.

.. link

::

    sage: G = F.quadratic_twist(2); G
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: G.conductor()
    2368
    sage: G.j_invariant()
    110592/37

We can compute the coefficients :math:`a_n` of the
:math:`L`-series or modular form
:math:`\sum_{n=0}^\infty
a_nq^n` attached to the elliptic curve.
This computation uses the PARI C-library:

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E.anlist(30)
    [0, 1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0, -4,
     3, 10, 2, 0, -1, 4, -9, -2, 6, -12]
    sage: v = E.anlist(10000)

It only takes a second to compute all :math:`a_n` for
:math:`n\leq 10^5`:

.. skip

::

    sage: %time v = E.anlist(100000)
    CPU times: user 0.98 s, sys: 0.06 s, total: 1.04 s
    Wall time: 1.06

Elliptic curves can be constructed using their Cremona labels. This
pre-loads the elliptic curve with information about its rank,
Tamagawa numbers, regulator, etc.

::

    sage: E = EllipticCurve("37b2")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational
    Field
    sage: E = EllipticCurve("389a")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x  over Rational Field
    sage: E.rank()
    2
    sage: E = EllipticCurve("5077a")
    sage: E.rank()
    3

We can also access the Cremona database directly.

::

    sage: db = sage.databases.cremona.CremonaDatabase()
    sage: db.curves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1], 'b1': [[0, 1, 1, -23, -50], 0, 3]}
    sage: db.allcurves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1],
     'b1': [[0, 1, 1, -23, -50], 0, 3],
     'b2': [[0, 1, 1, -1873, -31833], 0, 1],
     'b3': [[0, 1, 1, -3, 1], 0, 3]}

The objects returned from the database are not of type
``EllipticCurve``. They are elements of a database and have a couple
of fields, and that's it. There is a small version of Cremona's
database, which is distributed by default with Sage, and contains
limited information about elliptic curves of conductor
:math:`\leq 10000`. There is also a large optional version, which
contains extensive data about all curves of conductor up to
:math:`120000` (as of October 2005). There is also a huge (2GB)
optional database package for Sage that contains the hundreds of
millions of elliptic curves in the Stein-Watkins database.

Dirichlet Characters
--------------------

A *Dirichlet character* is the extension of a homomorphism
:math:`(\ZZ/N\ZZ)^* \to R^*`, for some ring :math:`R`, to the map
:math:`\ZZ \to R` obtained by sending those integers :math:`x`
with :math:`\gcd(N,x)>1` to 0.

::

    sage: G = DirichletGroup(12)
    sage: G.list()
    [Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1,
    Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1]
    sage: G.gens()
    (Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1)
    sage: len(G)
    4

Having created the group, we next create an element and compute
with it.

.. link

::

    sage: G = DirichletGroup(21)
    sage: chi = G.1; chi
    Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6
    sage: chi.values()
    [0, 1, zeta6 - 1, 0, -zeta6, -zeta6 + 1, 0, 0, 1, 0, zeta6, -zeta6, 0, -1,
     0, 0, zeta6 - 1, zeta6, 0, -zeta6 + 1, -1]
    sage: chi.conductor()
    7
    sage: chi.modulus()
    21
    sage: chi.order()
    6
    sage: chi(19)
    -zeta6 + 1
    sage: chi(40)
    -zeta6 + 1

It is also possible to compute the action of the Galois group
:math:`\text{Gal}(\QQ(\zeta_N)/\QQ)` on these characters, as well
as the direct product decomposition corresponding to the
factorization of the modulus.

.. link

::

    sage: chi.galois_orbit()
    [Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> -zeta6 + 1,
     Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6]

    sage: go = G.galois_orbits()
    sage: [len(orbit) for orbit in go]
    [1, 2, 2, 1, 1, 2, 2, 1]

    sage: G.decomposition()
    [Group of Dirichlet characters modulo 3 with values in Cyclotomic Field of order 6 and degree 2,
     Group of Dirichlet characters modulo 7 with values in Cyclotomic Field of order 6 and degree 2]

Next, we construct the group of Dirichlet characters mod 20, but
with values in :math:`\QQ(i)`:

::

    sage: K.<i> = NumberField(x^2+1)
    sage: G = DirichletGroup(20,K)
    sage: G
    Group of Dirichlet characters modulo 20 with values in Number Field in i with defining polynomial x^2 + 1


We next compute several invariants of ``G``:

.. link

::

    sage: G.gens()
    (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1,
    Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> i)

    sage: G.unit_gens()
    (11, 17)
    sage: G.zeta()
    i
    sage: G.zeta_order()
    4

In this example we create a Dirichlet character with values in a
number field. We explicitly specify the choice of root of unity by
the third argument to ``DirichletGroup`` below.

::

    sage: x = polygen(QQ, 'x')
    sage: K = NumberField(x^4 + 1, 'a'); a = K.0
    sage: b = K.gen(); a == b
    True
    sage: K
    Number Field in a with defining polynomial x^4 + 1
    sage: G = DirichletGroup(5, K, a); G
    Group of Dirichlet characters modulo 5 with values in the group of order 8 generated by a in Number Field in a with defining polynomial x^4 + 1
    sage: chi = G.0; chi
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> a^2
    sage: [(chi^i)(2) for i in range(4)]
    [1, a^2, -1, -a^2]

Here ``NumberField(x^4 + 1, 'a')`` tells Sage to use the symbol "a" in
printing what ``K`` is (a Number Field in a with defining polynomial
:math:`x^4 + 1`). The name "a" is undeclared at this point. Once
``a = K.0`` (or equivalently ``a = K.gen()``) is evaluated, the symbol
"a" represents a root of the generating polynomial
:math:`x^4+1`.

Modular Forms
-------------

Sage can do some computations related to modular forms, including
dimensions, computing spaces of modular symbols, Hecke operators,
and decompositions.

There are several functions available for computing dimensions of
spaces of modular forms. For example,

::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(11),2)
    1
    sage: dimension_cusp_forms(Gamma0(1),12)
    1
    sage: dimension_cusp_forms(Gamma1(389),2)
    6112

Next we illustrate computation of Hecke operators on a space of
modular symbols of level :math:`1` and weight :math:`12`.

::

    sage: M = ModularSymbols(1,12)
    sage: M.basis()
    ([X^8*Y^2,(0,0)], [X^9*Y,(0,0)], [X^10,(0,0)])
    sage: t2 = M.T(2)
    sage: t2
    Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1)
    of weight 12 with sign 0 over Rational Field
    sage: t2.matrix()
    [ -24    0    0]
    [   0  -24    0]
    [4860    0 2049]
    sage: f = t2.charpoly('x'); f
    x^3 - 2001*x^2 - 97776*x - 1180224
    sage: factor(f)
    (x - 2049) * (x + 24)^2
    sage: M.T(11).charpoly('x').factor()
    (x - 285311670612) * (x - 534612)^2

We can also create spaces for :math:`\Gamma_0(N)` and
:math:`\Gamma_1(N)`.

::

    sage: ModularSymbols(11,2)
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
     0 over Rational Field
    sage: ModularSymbols(Gamma1(11),2)
    Modular Symbols space of dimension 11 for Gamma_1(11) of weight 2 with
    sign 0 over Rational Field

Let's compute some characteristic polynomials and
:math:`q`-expansions.

::

    sage: M = ModularSymbols(Gamma1(11),2)
    sage: M.T(2).charpoly('x')
    x^11 - 8*x^10 + 20*x^9 + 10*x^8 - 145*x^7 + 229*x^6 + 58*x^5 - 360*x^4
         + 70*x^3 - 515*x^2 + 1804*x - 1452
    sage: M.T(2).charpoly('x').factor()
    (x - 3) * (x + 2)^2 * (x^4 - 7*x^3 + 19*x^2 - 23*x + 11)
            * (x^4 - 2*x^3 + 4*x^2 + 2*x + 11)
    sage: S = M.cuspidal_submodule()
    sage: S.T(2).matrix()
    [-2  0]
    [ 0 -2]
    sage: S.q_expansion_basis(10)
    [q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)]

We can even compute spaces of modular symbols with character.

::

    sage: G = DirichletGroup(13)
    sage: e = G.0^2
    sage: M = ModularSymbols(e,2); M
    Modular Symbols space of dimension 4 and level 13, weight 2, character
    [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
    sage: M.T(2).charpoly('x').factor()
    (x - zeta6 - 2) * (x - 2*zeta6 - 1) * (x + zeta6 + 1)^2
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 2 of Modular Symbols space of
    dimension 4 and level 13, weight 2, character [zeta6], sign 0, over
    Cyclotomic Field of order 6 and degree 2
    sage: S.T(2).charpoly('x').factor()
    (x + zeta6 + 1)^2
    sage: S.q_expansion_basis(10)
    [q + (-zeta6 - 1)*q^2 + (2*zeta6 - 2)*q^3 + zeta6*q^4 + (-2*zeta6 + 1)*q^5 + (-2*zeta6 + 4)*q^6 + (2*zeta6 - 1)*q^8 - zeta6*q^9 + O(q^10)]

Here is another example of how Sage can compute the action of Hecke
operators on a space of modular forms.

::

    sage: T = ModularForms(Gamma0(11),2)
    sage: T
    Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of
    weight 2 over Rational Field
    sage: T.degree()
    2
    sage: T.level()
    11
    sage: T.group()
    Congruence Subgroup Gamma0(11)
    sage: T.dimension()
    2
    sage: T.cuspidal_subspace()
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
    Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: T.eisenstein_subspace()
    Eisenstein subspace of dimension 1 of Modular Forms space of dimension 2
    for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: M = ModularSymbols(11); M
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
    0 over Rational Field
    sage: M.weight()
    2
    sage: M.basis()
    ((1,0), (1,8), (1,9))
    sage: M.sign()
    0

Let :math:`T_p` denote the usual Hecke operators (:math:`p`
prime). How do the Hecke operators :math:`T_2`, :math:`T_3`,
:math:`T_5` act on the space of modular symbols?

.. link

::

    sage: M.T(2).matrix()
    [ 3  0 -1]
    [ 0 -2  0]
    [ 0  0 -2]
    sage: M.T(3).matrix()
    [ 4  0 -1]
    [ 0 -1  0]
    [ 0  0 -1]
    sage: M.T(5).matrix()
    [ 6  0 -1]
    [ 0  1  0]
    [ 0  0  1]
