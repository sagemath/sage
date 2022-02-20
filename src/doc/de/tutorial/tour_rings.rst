.. _section-rings:

Wichtige Ringe
==============

Wenn wir Matrizen, Vektoren oder Polynome definieren ist es manchmal
nützlich, und manchmal notwendig, den "Ring" anzugeben, über dem diese
definiert sind. Ein *Ring* ist ein mathematisches Objekt, für das es
die wohldefinierten Operationen Addition und Multiplikation gibt; Falls
Sie davon noch nie gehört haben, brauchen Sie wahrscheinlich nur die
folgenden vier häufig verwendeten Ringe zu kennen.

* die ganzen Zahlen `\{..., -1, 0, 1, 2, ...\}`, welche ``ZZ`` in Sage
  genannt werden.
* die rationalen Zahlen -- d.h  Brüche oder Quotienten von ganzen
  Zahlen, welche ``QQ`` in Sage genannt werden.
* die reellen Zahlen, welche ``RR`` in Sage genannt werden.
* die komplexen Zahlen, welche ``CC`` in Sage genannt werden.

Sie müssen diese Unterschiede kennen, da das gleiche Polynom zum
Beispiel, unterschiedlich, abhängig von dem Ring über dem es definiert
ist, behandelt werden kann. Zum Beispiel hat das Polynom `x^2-2` die
beiden Nullstellen `\pm \sqrt{2}`.  Diese Nullstellen sind nicht
rational, wenn Sie also mit Polynomen über rationalen Koeffizienten
arbeiten, lässt sich das Polynom nicht faktorisieren. Mit reellen
Koeffizienten lässt es sich faktorisieren. Deshalb müssen Sie den Ring
angeben, um sicher zu gehen, dass Sie die Information
erhalten, die Sie erwarten. Die folgenden beiden Befehle definieren
jeweils die Mengen der Polynome mit rationalen und reellen
Koeffizienten. Diese Mengen werden "ratpoly" und "realpoly" genannt,
aber das ist hier nicht wichtig; beachten Sie jedoch, dass die Strings
".<t>" und ".<z>" die *Variablen* benennen, die in beiden Fällen benutzt
werden. ::

    sage: ratpoly.<t> = PolynomialRing(QQ)
    sage: realpoly.<z> = PolynomialRing(RR)

Jetzt verdeutlichen wir die Behauptung über das Faktorisieren von `x^2-2`:

.. link

::

    sage: factor(t^2-2)
    t^2 - 2
    sage: factor(z^2-2)
    (z - 1.41421356237310) * (z + 1.41421356237310)

Ähnliche Kommentare treffen auf Matrizen zu: Die zeilenreduzierte Form
eine Matrix kann vom Ring abhängen, über dem sie definiert ist,
genauso wie ihre Eigenwerte und Eigenvektoren. Um mehr über das
Konstruieren von Polynomen zu erfahren, lesen Sie :ref:`section-poly`,
und für mehr über Matrizen, lesen Sie :ref:`section-linalg`.

Das Symbol ``I`` steht für die Quadratwurzel von :math:`-1`; ``i`` ist
ein Synonym für ``I``. Natürlich ist dies keine rationale Zahl::

    sage: i  # Wurzel von -1
    I
    sage: i in QQ
    False

Beachten Sie: Der obige Code kann möglicherweise nicht wie erwartet
funktionieren. Zum Beispiel wenn der Variablen ``i`` ein
unterschiedlicher Wert, etwa wenn diese als Schleifenvariable
verwendet wurde, zugewiesen wurde. Falls dies der Fall ist, tippen Sie ::

    sage: reset('i')

um den ursprünglichen komplexen Wert der Variable ``i`` zu erhalten.

Es ist noch eine Feinheit beim Definieren von komplexen Zahlen zu
beachten: Wie oben erwähnt wurde, stellt das Symbol ``i`` eine
Quadratwurzel von `-1` dar, es ist jedoch eine *formale* Quadratwurzel
von `-1`, jedoch eine algebraische Zahl.  Das Aufrufen von ``CC(i)``
oder ``CC.0`` oder ``CC.gen(0)``, gibt die *komplexe* Quadratwurzel
von `-1` zurück. ::

    sage: i = CC(i)       # komplexe Gleitkommazahl
    sage: i == CC.0
    True
    sage: a, b = 4/3, 2/3
    sage: z = a + b*i
    sage: z
    1.33333333333333 + 0.666666666666667*I
    sage: z.imag()        # Imaginärteil
    0.666666666666667
    sage: z.real() == a   # automatische Umwandlung vor dem Vergleich
    True
    sage: a + b
    2
    sage: 2*b == a
    True
    sage: parent(2/3)
    Rational Field
    sage: parent(4/2)
    Rational Field
    sage: 2/3 + 0.1       # automatische Umwandlung vor der Addition
    0.766666666666667
    sage: 0.1 + 2/3       # Umwandlungsregeln sind symmetrisch in Sage
    0.766666666666667

Hier sind weitere Beispiele von Ringen in Sage. Wie oben angemerkt,
kann auf den Ring der rationalen Zahlen mit ``QQ`` zugegriffen werden,
ebenso wie mit ``RationalField()`` (ein  *Körper* (engl. *field*) ist
ein Ring in dem die Multiplikation kommutativ ist, und in dem jedes von
Null verschiedene Element in dem Ring einen Kehrwehrt besitzt. Die
rationalen Zahlen bilden also auch einen Körper, die ganzen Zahlen
jedoch nicht)::


    sage: RationalField()
    Rational Field
    sage: QQ
    Rational Field
    sage: 1/2 in QQ
    True

Die Dezimalzahl ``1.2`` wird als rationale Zahl in ``QQ`` gesehen:
Dezimalzahlen, die auch rational sind, können in rationale Zahlen
"umgewandelt" (engl. "coerced") werden. Die Zahlen `\pi` und `\sqrt{2}`
sind jedoch nicht rational::

    sage: 1.2 in QQ
    True
    sage: pi in QQ
    False
    sage: pi in RR
    True
    sage: sqrt(2) in QQ
    False
    sage: sqrt(2) in CC
    True

Für die Verwendung in der höheren Mathematik kennt Sage noch weitere
Ringe, wie z.B. endliche Körper, `p`-adische Zahlen, den Ring der
algebraischen Zahlen, Polynomringe und Matrizenringe. Hier sind
Konstruktionen einiger von ihnen::

    sage: GF(3)
    Finite Field of size 3
    sage: GF(27, 'a')  # Sie müssen den Names des Generators angeben \
    ....:              # wenn es sich um keinen Primkörper handelt
    Finite Field in a of size 3^3
    sage: Zp(5)
    5-adic Ring with capped relative precision 20
    sage: sqrt(3) in QQbar # algebraischer Abschluss von QQ
    True
