Etwas weiter fortgeschrittene Mathematik
========================================

Algebraische Geometrie
----------------------

Sie können in Sage beliebige algebraische Varietäten definieren, aber
manchmal ist die nichttriviale Funktionalität auf Ringe über
:math:`\QQ` oder endliche Körper beschränkt. Zum Beispiel können wir
die Vereinigung zweier affiner, planarer Kurven berechnen, und dann
die Kurven als irreduzible Komponenten der Vereinigung zurück erhalten.

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

Wir können auch alle Punkte im Schnitt der beiden Kurven finden, indem
wir diese schneiden und dann die irreduziblen Komponenten berechnen.

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

Also sind zum Beispiel :math:`(1,0)` und :math:`(0,1)` auf beiden
Kurven (wie man sofort sieht), genauso wie bestimmte (quadratischen)
Punkte, deren :math:`y` Koordinaten :math:`2y^2 + 4y + 3=0` erfüllen.

Sage kann auch das torische Ideal der gedrehten Kubik im
dreidimensionalen projektiven Raum berechnen:


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

Elliptische Kurven
------------------

Die Funktionalität elliptischer Kurven beinhaltet die meisten von
PARIs Funktionen zu elliptischen Kurven, den Zugriff auf die Daten von
Cremonas Online-Tabellen (dies benötigt ein optionales
Datenbankpaket), die Funktionen von mwrank, d.h.
2-Abstiege mit der Berechnung der vollen Mordell-Weil-Gruppe, der SEA
Algorithmus, Berechnung aller Isogenien, viel neuem Code für Kurven
über :math:`\QQ` und Teile von Denis Simons "algebraic descent" Software.

Der Befehl ``EllipticCurve`` zum Erzeugen von Elliptischen Kurven hat
viele Formen:


-  EllipticCurve([:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6`]):
   Gibt die elliptische Kurve

   .. math::  y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6,


   zurück, wobei die :math:`a_i`'s umgewandelt werden zum Typ von
   :math:`a_1`. Falls alle :math:`a_i` den Typ :math:`\ZZ` haben,
   werden sie zu :math:`\QQ` umgewandelt.

-  EllipticCurve([:math:`a_4`, :math:`a_6`]): Das Gleiche wie oben,
   jedoch ist :math:`a_1=a_2=a_3=0`.

-  EllipticCurve(label): Gibt die elliptische Kurve aus der Datenbank
   von Cremona mit dem angegebenen (neuen!) Cremona-Label zurück. Das
   Label ist ein String, wie z.B. ``"11a"`` oder ``"37b2"``. Die
   Buchstaben müssen klein geschrieben sein (um sie von dem alten
   Label unterscheiden zu können).

-  EllipticCurve(j): Gibt die elliptische Kurve mit
   :math:`j`-Invariante :math:`j` zurück.

-  EllipticCurve(R,
   [:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6`]):
   Erzeugt die elliptische Kurve über dem Ring :math:`R` mit
   vorgegebenen :math:`a_i`'s wie oben.


Wir veranschaulichen jede dieser Konstruktionen:

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

Das Paar :math:`(0,0)` ist ein Punkt auf der elliptischen Kurve
:math:`E` definiert durch :math:`y^2 + y = x^3 - x`. Um diesen Punkt
in Sage zu erzeugen, geben Sie ``E([0,0])`` ein. Sage kann auf einer
solchen elliptischen Kurve Punkte addieren (erinnern Sie sich:
elliptische Kurven haben eine additive Gruppenstruktur, wobei der unendlich
ferne Punkt das Nullelement ist, und drei kollineare Punkte auf
der Kurve sich zu Null addieren):

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

Die elliptischen Kurven über den komplexen Zahlen sind durch die
:math:`j`-Invariante parametrisiert. Sage berechnet
:math:`j`-Invarianten wie folgt:

::

    sage: E = EllipticCurve([0,0,0,-4,2]); E
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: E.conductor()
    2368
    sage: E.j_invariant()
    110592/37

Wenn wir eine Kurve mit der gleichen :math:`j`-Invarianten wie
:math:`E` erstellen, muss diese nicht isomorph zu :math:`E` sein. Im
folgenden Beispiel sind die Kurven nicht isomorph, da ihre Führer
unterschiedlich sind.

::

    sage: F = EllipticCurve_from_j(110592/37)
    sage: F.conductor()
    37

Jedoch ergibt der Twist von :math:`F` mit 2 eine isomorphe Kurve.

.. link

::

    sage: G = F.quadratic_twist(2); G
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: G.conductor()
    2368
    sage: G.j_invariant()
    110592/37

Wir können die Koeffizienten :math:`a_n` der zur elliptischen Kurve
gehörenden :math:`L`-Reihe oder der Modulform
:math:`\sum_{n=0}^\infty a_nq^n` berechnen.
Die Berechnung benutzt die PARI C-Bibliothek:

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E.anlist(30)
    [0, 1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0, -4,
     3, 10, 2, 0, -1, 4, -9, -2, 6, -12]
    sage: v = E.anlist(10000)

Alle Koeffizienten :math:`a_n` bis zu :math:`n\leq 10^5` zu berechnen
dauert nur eine Sekunde:

.. skip

::

    sage: %time v = E.anlist(100000)
    CPU times: user 0.98 s, sys: 0.06 s, total: 1.04 s
    Wall time: 1.06

Elliptische Kurven können mit Hilfe ihres Cremona-Labels konstruiert
werden. Dies lädt die Kurve zusammen mit Informationen über ihren Rank, mit
Tamagawa Zahlen, Regulatoren, usw..

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

Wir können auch direkt auf die Cremona-Datenbank zugreifen.

::

    sage: db = sage.databases.cremona.CremonaDatabase()
    sage: db.curves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1], 'b1': [[0, 1, 1, -23, -50], 0, 3]}
    sage: db.allcurves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1],
     'b1': [[0, 1, 1, -23, -50], 0, 3],
     'b2': [[0, 1, 1, -1873, -31833], 0, 1],
     'b3': [[0, 1, 1, -3, 1], 0, 3]}

Die Objekte, die aus der Datenbank zurückgegeben werden, sind nicht
vom Typ ``EllipticCurve``. Sie sind Elemente einer Datenbank und haben
ein paar Komponenten, und das war's. Es gibt eine kleine Version von
Cremonas Datenbank, die standardmäßig zu Sage gehört und beschränkte
Information zu elliptischen Kurven mit Führer :math:`\leq 10000`
enthält. Es gibt auch eine große optionale Version, welche ausgiebige
Daten zu allen elliptischen Kurven mit Führer bis zu :math:`120000`
enthält (Stand Oktober 2005). Es gibt auch ein riesiges (2GB großes)
optionales Datenbank-Paket für Sage, das in der Stein-Watkins
Datenbank hunderte Millionen von elliptischen Kurven enthält.

Dirichlet-Charaktere
--------------------

Ein *Dirichlet Charakter* ist die Erweiterung eines Homomorphismus
:math:`(\ZZ/N\ZZ)^* \to R^*`, für einen Ring :math:`R`, zu der Abbildung
:math:`\ZZ \to R`, welche erhalten wird, wenn man diese ganzen Zahlen :math:`x`
mit :math:`\gcd(N,x)>1` nach :math:`0` schickt.

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

Nachdem wir dies Gruppe erzeugt haben, erstellen wir als nächstes ein
Element und rechnen damit.

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

Es ist auch möglich die Operation der Galoisgruppe
:math:`\text{Gal}(\QQ(\zeta_N)/\QQ)` auf diesen Charakteren zu
berechnen, sowie die Zerlegung in direkte Produkte, die der
Faktorisierung des Moduls entsprechen.


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

Als nächstes konstruieren wir die Gruppe der Dirichlet-Charaktere
mod 20, jedoch mit Werten in :math:`\QQ(i)`:

::

    sage: K.<i> = NumberField(x^2+1)
    sage: G = DirichletGroup(20,K)
    sage: G
    Group of Dirichlet characters modulo 20 with values in Number Field in i with defining polynomial x^2 + 1


Nun berechnen wir mehrere Invarianten von ``G``:

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

In diesem Beispiel erzeugen wir einen Dirichlet-Charakter mit Werten
in einem Zahlenfeld. Wir geben die Wahl der Einheitswurzel im dritten
Argument von ``DirichletGroup`` an.

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

Hier teilt ``NumberField(x^4 + 1, 'a')`` Sage mit, dass es das Symbol "a"
beim Ausgeben dessen was ``K`` ist (ein Zahlenfeld mit definierendem Polynom
:math:`x^4 + 1`) benutzen soll. Der Name "a" ist zu diesem Zeitpunkt
nicht deklariert. Sobald ``a = K.0`` (oder äquivalent ``a = K.gen()``)
evaluiert wurde, repräsentiert das Symbol "a" eine Wurzel des
erzeugenden Polynoms :math:`x^4+1`.


Modulformen
-------------

Sage kann einige Berechnungen im Zusammenhang mit Modulformen
durchführen, einschließlich Dimensionsberechnungen, das Berechnen von Räumen von
Symbolen von Modulformen, Hecke-Operatoren, und Dekompositionen.

Es stehen mehrere Funktionen für das Berechnen von Dimensionen von
Räumen von Modulformen zur Verfügung. Zum Beispiel,

::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(11),2)
    1
    sage: dimension_cusp_forms(Gamma0(1),12)
    1
    sage: dimension_cusp_forms(Gamma1(389),2)
    6112

Als nächstes illustrieren wir die Berechnung von Hecke-Operatoren auf
einem Raum von Modulformen von Level :math:`1` und Gewicht :math:`12`.

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

Wir können auch Räume für :math:`\Gamma_0(N)` und :math:`\Gamma_1(N)`
erzeugen.

::

    sage: ModularSymbols(11,2)
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
     0 over Rational Field
    sage: ModularSymbols(Gamma1(11),2)
    Modular Symbols space of dimension 11 for Gamma_1(11) of weight 2 with
    sign 0 over Rational Field

Nun berechnen wir ein paar charakteristische Polynome und
:math:`q`-Entwicklungen.

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

Wir können sogar Räume von Modulsymbolen mit Charakteren berechnen.

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

Hier ist ein weiteres Beispiel davon wie Sage mit den Operationen von
Hecke-Operatoren auf dem Raum von Modulformen rechnen kann.

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

Sei :math:`T_p` die Bezeichnung der gewöhnlichen Hecke-Operatoren (:math:`p`
prim). Wie operieren die Hecke-Operatoren :math:`T_2`, :math:`T_3`,
:math:`T_5` auf dem Raum der Modulsymbole?

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
