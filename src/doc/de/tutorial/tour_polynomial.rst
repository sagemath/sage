.. _section-poly:

Polynome
========

In diesem Abschnitt erläutern wir, wie man in Sage Polynome erzeugt
und benutzt.


.. _section-univariate:

Polynome in einer Unbestimmten
------------------------------

Es gibt drei Möglichkeiten Polynomringe zu erzeugen.

::

    sage: R = PolynomialRing(QQ, 't')
    sage: R
    Univariate Polynomial Ring in t over Rational Field

Dies erzeugt einen Polynomring und teilt Sage mit (den String) 't' als
Unbestimmte bei Ausgaben auf dem Bildschirm zu verwenden.
Jedoch definiert dies nicht das Symbol ``t`` zur Verwendung in Sage,
Sie können es also nicht verwenden um ein Polynom (wie
z.B. :math:`t^2+1`) einzugeben, welches zu ``R`` gehört.

Eine alternative Möglichkeit ist:

.. link

::

    sage: S = QQ['t']
    sage: S == R
    True

Dies verhält sich bezüglich ``t`` gleich.

Eine dritte sehr bequeme Möglichkeit ist:

::

    sage: R.<t> = PolynomialRing(QQ)

oder

::

    sage: R.<t> = QQ['t']

oder sogar nur

::

    sage: R.<t> = QQ[]

Dies hat den zusätzlichen Nebeneffekt, dass die Variable ``t`` als
Unbestimmte des Polynomrings definiert wird, Sie können daher nun wie
folgt auf einfache Weise Elemente von ``R`` definieren. (Beachten Sie,
dass die dritte Möglichkeit sehr ähnlich zu der Konstruktor-Notation
in Magma ist und, genau wie in Magma, kann Sie dazu verwendet werden
eine Vielzahl von Objekten zu definieren.)

.. link

::

    sage: poly = (t+1) * (t+2); poly
    t^2 + 3*t + 2
    sage: poly in R
    True

Unabhängig davon wie Sie den Polynomring definieren, können Sie die
Unbestimmte als den :math:`0^{th}` Erzeuger zurückerhalten:

::

    sage: R = PolynomialRing(QQ, 't')
    sage: t = R.0
    sage: t in R
    True

Beachten Sie, dass Sie bei den komplexen Zahlen eine ähnliche
Konstruktion verwenden können: Sie können die komplexen Zahlen
ansehen, als wären sie von dem Symbol ``i`` über den reellen Zahlen
erzeugt; wir erhalten also Folgendes:

::

    sage: CC
    Complex Field with 53 bits of precision
    sage: CC.0  # 0th generator of CC
    1.00000000000000*I

Beim Erzeugen von Polynomringen kann man sowohl den Ring, als auch den
Erzeuger, oder nur den Erzeuger wie folgt erhalten:

::

    sage: R, t = QQ['t'].objgen()
    sage: t    = QQ['t'].gen()
    sage: R, t = objgen(QQ['t'])
    sage: t    = gen(QQ['t'])

Schließlich treiben wir etwas Arithmetik in :math:`\QQ[t]`.

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

Beachten Sie, dass die Faktorisierung die Einheit korrekt in Betracht
zieht und ausgibt.

Falls Sie zum Beispiel die ``R.cyclotomic_polynomial``-Funktion in
einem Forschungsprojekt viel verwenden würden, sollten Sie neben Sage
zu zitieren, auch versuchen herauszufinden welche Komponente von Sage verwendet
wird um das zyklotomische Polynom zu berechnen, und diese ebenso angeben.
In diesen Fall sehen Sie im Quellcode der Funktion, welchen Sie mit
``R.cyclotomic_polynomial??`` erhalten schnell die Zeile ``f =
pari.polcyclo(n)``, was bedeutet, dass PARI verwendet wird um das
zyklotomische Polynom zu berechnen. Zitieren Sie PARI ebenso in Ihrer Arbeit.

Wenn Sie zwei Polynome teilen, erzeugen Sie ein Element des Quotientenkörpers
(den Sage automatisch erzeugt).

::

    sage: x = QQ['x'].0
    sage: f = x^3 + 1; g = x^2 - 17
    sage: h = f/g;  h
    (x^3 + 1)/(x^2 - 17)
    sage: h.parent()
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

Mit Hilfe von Laurentreihen können Sie die Reihenentwicklung im
Quotientenkörper von ``QQ[x]`` berechnen:

::

    sage: R.<x> = LaurentSeriesRing(QQ); R
    Laurent Series Ring in x over Rational Field
    sage: 1/(1-x) + O(x^10)
    1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)

Wenn wir die Variablen unterschiedlich benennen, erhalten wir einen
unterschiedlichen Polynomring.

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

Der Ring wird durch die Variable bestimmt. Beachten Sie, dass das
Erzeugen eines weiteren Rings mit einer ``x`` genannten Variablen
keinen unterschiedlichen Ring zurück gibt.

::

    sage: R = PolynomialRing(QQ, "x")
    sage: T = PolynomialRing(QQ, "x")
    sage: R == T
    True
    sage: R is T
    True
    sage: R.0 == T.0
    True

Sage unterstützt auch Ringe von Potenz- und Laurentreihen über
beliebigen Ringen. Im folgenden Beispiel erzeugen wir ein Element aus
:math:`\GF{7}[[T]]` und teilen es um ein Element aus :math:`\GF{7}((T))`
zu erhalten.

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

Sie können einen Potenzreihenring auch mit der Kurzschreibweise,
doppelter eckiger Klammern erzeugen:

::

    sage: GF(7)[['T']]
    Power Series Ring in T over Finite Field of size 7

Polynome in mehreren Unbestimmten
---------------------------------

Um mit Polynomringen in mehreren Variablen zu arbeiten, deklarieren
wir zunächst den Ring und die Variablen.

::

    sage: R = PolynomialRing(GF(5),3,"z") # here, 3 = number of variables
    sage: R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Genau wie bei dem Definieren von Polynomringen in einer Variablen,
gibt es mehrere Möglichkeiten:

::

    sage: GF(5)['z0, z1, z2']
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5
    sage: R.<z0,z1,z2> = GF(5)[]; R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Falls die Variablennamen nur einen Buchstaben lang sein sollen, können
Sie auch die folgende Kurzschreibweise verwenden:

::

    sage: PolynomialRing(GF(5), 3, 'xyz')
    Multivariate Polynomial Ring in x, y, z over Finite Field of size 5

Als nächstes treiben wir wieder etwas Arithmetik.

::

    sage: z = GF(5)['z0, z1, z2'].gens()
    sage: z
    (z0, z1, z2)
    sage: (z[0]+z[1]+z[2])^2
    z0^2 + 2*z0*z1 + z1^2 + 2*z0*z2 + 2*z1*z2 + z2^2

Sie können auch eine mathematisch etwas weiter verbreitete
Schreibweise verwenden um den Polynomring zu definieren.

::

    sage: R = GF(5)['x,y,z']
    sage: x,y,z = R.gens()
    sage: QQ['x']
    Univariate Polynomial Ring in x over Rational Field
    sage: QQ['x,y'].gens()
    (x, y)
    sage: QQ['x'].objgens()
    (Univariate Polynomial Ring in x over Rational Field, (x,))

Polynomringe in mehreren Variablen sind in Sage mit Hilfe von
Python-Dictionaries und der "distributiven Darstellung" eines Polynoms
implementiert. Sage benutzt Singular [Si]_, zum Beispiel bei der
Berechnung von ggTs und Gröbnerbasen von Idealen.

::

    sage: R, (x, y) = PolynomialRing(RationalField(), 2, 'xy').objgens()
    sage: f = (x^3 + 2*y^2*x)^2
    sage: g = x^2*y^2
    sage: f.gcd(g)
    x^2

Als nächstes erstellen wir das Ideal :math:`(f,g)` welches von
:math:`f` und :math:`g` erzeugt wird, indem wir einfach ``(f,g)`` mit
``R`` multiplizieren (wir könnten auch ``ideal([f,g])`` oder
``ideal(f,g)``) schreiben.

.. link

::

    sage: I = (f, g)*R; I
    Ideal (x^6 + 4*x^4*y^2 + 4*x^2*y^4, x^2*y^2) of Multivariate Polynomial
    Ring in x, y over Rational Field
    sage: B = I.groebner_basis(); B
    [x^6, x^2*y^2]
    sage: x^2 in I
    False

Übrigens ist die obige Gröbnerbasis keine Liste, sondern eine
unveränderliche Folge. Das bedeutet das sie die Attribute "universe" und
"parent" besitzt und nicht verändert werden kann. (Dies ist nützlich,
da nach dem Ändern der Basis andere Routinen, welche die Gröbnerbasis
verwenden, nicht mehr funktionieren könnten)

.. link

::

    sage: B.universe()
    Multivariate Polynomial Ring in x, y over Rational Field
    sage: B[1] = x
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Etwas (damit meinen wir: nicht so viel wie wir gerne hätten)
kommutative Algebra ist in Sage, mit Hilfe von Singular implementiert,
vorhanden. Zum Beispiel können wir die Zerlegung in Primideale und die
assoziierten Primideale von :math:`I` berechnen.

.. link

::

    sage: I.primary_decomposition()
    [Ideal (x^2) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y^2, x^6) of Multivariate Polynomial Ring in x, y over Rational Field]
    sage: I.associated_primes()
    [Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field]
