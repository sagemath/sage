**************
Programmierung
**************

.. _section-loadattach:

Sage-Dateien Laden und Anhängen
===============================

Als nächstes zeigen wir wie man Programme, die einer separaten Datei
geschrieben wurden in Sage lädt. Erstellen Sie eine Datei, welche Sie
``beispiel.sage`` nennen mit folgendem Inhalt:

.. skip

::

    print("Hello World")
    print(2^3)

Sie können ``beispiel.sage`` einlesen und ausführen, indem Sie den
``load``-Befehl verwenden.

.. skip

::

    sage: load("beispiel.sage")
    Hello World
    8

Sie können auch eine Sage-Datei an eine laufende Sitzung anhängen,
indem Sie den ``attach``-Befehl verwenden:

.. skip

::

    sage: attach("beispiel.sage")
    Hello World
    8

Wenn Sie nun ``beispiel.sage`` verändern und eine Leerzeile in Sage eingeben
(d.h. ``return`` drücken) wird der Inhalt von ``beispiel.sage``
automatisch in Sage neu geladen.

Insbesondere lädt der ``attach``-Befehl eine Datei jedesmal, wenn
diese verändert wird automatisch neu, was beim Debuggen von Code
nützlich sein kann, wobei der ``load``-Befehl eine Datei nur einmal lädt.

Wenn Sage die Datei ``beispiel.sage`` lädt, wird sie zu Python-Code konvertiert,
welcher dann vom Python-Interpreter ausgeführt wird. Diese
Konvertierung ist geringfügig; sie besteht hautsächlich daraus
Integer-Literale mit ``Integer()`` und Fließkomma-Literale mit
``RealNumber()`` zu versehen, ``^`` durch ``**`` zu ersetzen und
z.B. ``R.2`` durch ``R.gen(2)`` auszutauschen. Die konvertierte
Version von ``beispiel.sage`` befindet sich im gleichen Verzeichnis wie
``beispiel.sage`` und ist ``beispiel.sage.py`` genannt. Diese Datei
enthält den folgenden Code:

::

    print("Hello World")
    print(Integer(2)**Integer(3))

Integer-Literale wurden mit Integer() versehen und das ``^`` wurde
durch ein ``**`` ersetzt.
(In Python bedeutet ``^`` "exklusives oder" und ``**`` bedeutet
"Exponentiation".)

Dieses "preparsing" ist in ``sage/misc/interpreter.py`` implementiert.)

Sie können mehrzeiligen eingerückten Code in Sage einfügen, solange
Zeilenenden neue Blöcke kenntlich machen (dies ist in Dateien nicht notwendig).
Jedoch beseht die beste Möglichkeit solchen Code in Sage einzufügen darin,
diesen in einer Datei zu speichern und ``attach`` wie oben beschrieben zu verwenden.


.. _section-compile:

Kompilierten Code erzeugen
==========================

Geschwindigkeit ist bei mathematischen Berechnungen äußerst
wichtig. Python ist zwar eine komfortable Programmiersprache mit sehr
hohen Abstraktionsniveau, jedoch können bestimmte Berechnungen
mehrere Größenordnungen schneller als in Python sein, wenn sie in
einer kompilierten Sprache mit statischen Datentypen implementiert
wurden. Manche Teile von Sage würden zu langsam sein, wenn sie
komplett in Python geschrieben wären. Um dies zu berücksichtigen
unterstützt Sage eine kompilierte "Version" von Python, welche Cython
([Cyt]_ und [Pyr]_) genannt wird. Cython ist gleichzeitig sowohl zu Python,
als auch zu C ähnlich. Die meisten von Pythons Konstruktionen,
einschließlich "list comprehensions", bedingte Ausdrücke und Code wie
``+=`` sind erlaubt; Sie können auch Code importieren, den Sie in
anderen Python-Modulen geschrieben haben. Darüberhinaus können Sie
beliebige C Variablen definieren und beliebe C-Bibliothekaufrufe
direkt ausführen. Der daraus entstehende Code wird nach C konvertiert
und mithilfe eines C-Compilers kompiliert.

Um eigenen kompilierten Sagecode zu erstellen, geben Sie der Datei eine
``.spyx`` Endung (anstelle von ``.sage``). Falls Sie mit der
Kommandozeile arbeiten, können Sie kompilierten Code genau wie
interpretierten Code anhängen und laden. (Im Moment wird das Anhängen
von Cythoncode vom Notebook aus nicht unterstützt).
Die tatsächliche Kompilierung wird "hinter den Kulissen" durchgeführt
ohne dass Sie explizit etwas tun müssen. Die komplierte "shared object library"
wird unter ``$HOME/.sage/temp/hostname/pid/spyx`` gespeichert. Diese Dateien
werden gelöscht wenn Sie Sage beenden.

Auf spyx-Dateien wird kein "preparsing" angewendet, d.h. ``1/3`` wird
in einer spyx-Datei zu 0 ausgewertet, anstelle der rationalen Zahl
:math:`1/3`.
Wenn ``foo`` eine Funktion in der Sage-Bibliothek ist die Sie
verwenden möchten, müssen Sie ``sage.all`` importieren und
``sage.all.foo`` benutzen.

::

    import sage.all
    def foo(n):
        return sage.all.factorial(n)

Auf C-Funktionen in separaten Dateien zugreifen
-----------------------------------------------

Es ist auch nicht schwer auf C-Funktions zuzugreifen  welche in
separaten \*.c Dateien definiert sind. Hier ist ein Beispiel. Erzeugen
Sie die Dateien ``test.c`` und ``test.spyx`` in dem gleichen
Verzeichnis mit den Inhalten:

Der reine C-Code: ``test.c``

::

    int add_one(int n) {
      return n + 1;
    }

Der Cython-Code: ``test.spyx``:

::

    cdef extern from "test.c":
        int add_one(int n)

    def test(n):
        return add_one(n)

Dann funktioniert das Folgende:

.. skip

::

    sage: attach("test.spyx")
    Compiling (...)/test.spyx...
    sage: test(10)
    11

Wenn die zusätzliche Bibliothek ``foo`` gebraucht wird um den C-Code,
der aus einer Cython-Datei generiert wurde zu kompilieren, fügen Sie
die Zeile ``clib foo`` zu dem Cython-Quellcode hinzu. Auf ähnliche
Weise kann eine zusätzliche Datei ``bar`` zu der Kompilierung mit der
Deklaration ``cfile bar`` hinzugefügt werden.

.. _section-standalone:

eigenständige Python/Sage Skripte
=================================

Das folgende eigenständige Sageskript faktorisiert ganze Zahlen,
Polynome, usw.:

::

    #!/usr/bin/env sage

    import sys

    if len(sys.argv) != 2:
        print("Usage: %s <n>" % sys.argv[0])
        print("Outputs the prime factorization of n.")
        sys.exit(1)

    print(factor(sage_eval(sys.argv[1])))

Um dieses Skript benutzen zu können muss ``SAGE_ROOT`` in ihrer
PATH-Umgebungsvariable enthalten sein. Falls das das obige Skript
``factor`` genannt wurde, ist hier ein beispielhafter Aufruf:

::

    bash $ ./factor 2006
    2 * 17 * 59

Datentypen
==========

Jedes Objekt hat in Sage einen wohldefinierten Datentyp. Python
besitzt eine Vielzahl von standardmäßiger elementarer Datentypen und die
Sage-Bibliothek fügt noch viele weitere hinzu. Zu Pythons
standardmäßigen Datentypen gehören Strings, Listen, Tupel, Ganzzahlen und
Gleitkommazahlen, wie hier zu sehen ist:

::

    sage: s = "sage"; type(s)
    <... 'str'>
    sage: s = 'sage'; type(s)      # Sie können einfache oder doppelte Anführungszeichen verwenden
    <... 'str'>
    sage: s = [1,2,3,4]; type(s)
    <... 'list'>
    sage: s = (1,2,3,4); type(s)
    <... 'tuple'>
    sage: s = int(2006); type(s)
    <... 'int'>
    sage: s = float(2006); type(s)
    <... 'float'>

Hierzu fügt Sage noch viele weitere hinzu. Zum Beispiel Vektorräume:

::

    sage: V = VectorSpace(QQ, 1000000); V
    Vector space of dimension 1000000 over Rational Field
    sage: type(V)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category'>

Nur bestimmte Funktionen können auf ``V`` aufgerufen werden. In
anderen mathematischen Softwaresystemem würde dies mit der
"Funktionalen"-Notation ``foo(V,...)`` geschehen. In Sage sind
bestimmte Funktionen an den Typ (oder der Klasse) von ``V`` angehängt,
und diese werden unter Benutzung einer objektorientierten Syntax,
wie in Java oder C++ aufgerufen. Zum Beispiel ``V.foo(...)``. Dies
hilft dabei eine Überfüllung des globalen Namensraums mit tausenden
von Funktionen zu vermeiden. Das bedeutet auch, dass viele
verschiedene Funktionen mit unterschiedlichen Funktionsweisen foo
genannt werden können, ohne dass der Typ des Arguments überprüft (oder
Case-Anweisungen ausgeführt) werden muss, um zu entscheiden welche
aufgerufen werden soll. Weiterhin ist die Funktion auch dann noch
verfügbar, wenn ihr Name zu einem anderen Zweck verwendet wurde. (Zum
Beispiel wenn Sie etwas ``zeta`` nennen und dann den Wert der
Riemannschen Zeta-Funktion bei 0.5 berechnen wollen, können Sie
immernoch ``s=.5; s.zeta()`` benutzen).

::

    sage: zeta = -1
    sage: s=.5; s.zeta()
    -1.46035450880959

In manchen sehr oft auftretenden Fällen wird auch die gewöhnliche
funktionale Notation unterstützt, da dies bequem ist und manche
mathematische Ausdrücke in objektorientierter Notation verwirrend
aussehen könnten. Hier sind einige Beispiele:

::

    sage: n = 2; n.sqrt()
    sqrt(2)
    sage: sqrt(2)
    sqrt(2)
    sage: V = VectorSpace(QQ,2)
    sage: V.basis()
        [(1, 0), (0, 1)]
    sage: basis(V)
        [(1, 0), (0, 1)]
    sage: M = MatrixSpace(GF(7), 2); M
    Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
    sage: A = M([1,2,3,4]); A
    [1 2]
    [3 4]
    sage: A.charpoly('x')
    x^2 + 2*x + 5
    sage: charpoly(A, 'x')
    x^2 + 2*x + 5

Um alle Member-Funktionen von :math:`A` anzuzeigen, können Sie die
Tab-Vervollständigung benutzen. Tippen Sie einfach ``A.``, dann die
``[tab]``-Taste auf Ihrer Tastatur, wie es in
:ref:`section-tabcompletion` beschrieben ist.

Listen, Tupel, und Folgen
=========================

Der Listen-Datentyp speichert Elemente eines beliebigen Typs. Wie in
C, C++, usw. (jedoch anders als in vielen gewöhnlichen
Computer-Algebra-Systemen), die Elemente der Liste werden bei
:math:`0` beginnend indiziert:

::

    sage: v = [2, 3, 5, 'x', SymmetricGroup(3)]; v
    [2, 3, 5, 'x', Symmetric group of order 3! as a permutation group]
    sage: type(v)
    <... 'list'>
    sage: v[0]
    2
    sage: v[2]
    5

(Wenn man auf ein Listenelement zugreift ist es OK wenn der Index
kein Python int ist!)
Mit einem Sage-Integer (oder Rational, oder mit allem anderen mit einer ``__index__`` Methode)
funktioniert es genauso.

::

    sage: v = [1,2,3]
    sage: v[2]
    3
    sage: n = 2      # Sage Integer
    sage: v[n]       # Perfectly OK!
    3
    sage: v[int(n)]  # Also OK.
    3

Die ``range``-Funktion erzeugt eine Liste von Python int's (nicht
Sage-Integers):

::

    sage: list(range(1, 15))
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

Dies ist nützlich wenn man List-Comprehensions verwendet um Listen zu
konstruieren:

::

    sage: L = [factor(n) for n in range(1, 15)]
    sage: L
    [1, 2, 3, 2^2, 5, 2 * 3, 7, 2^3, 3^2, 2 * 5, 11, 2^2 * 3, 13, 2 * 7]
    sage: L[12]
    13
    sage: type(L[12])
     <class 'sage.structure.factorization_integer.IntegerFactorization'>
    sage: [factor(n) for n in range(1, 15) if is_odd(n)]
    [1, 3, 5, 7, 3^2, 11, 13]

Um mehr darüber zu erfahren wie man Listen mit Hilfe von
List-Comprehensions erzeugt, lesen Sie [PyT]_.

List-Slicing ist eine wunderbare Eigenschaft. Wenn ``L`` eine Liste
ist, dann gibt ``L[m:n]`` die Teilliste von ``L`` zurück, die erhalten
wird wenn man mit dem :math:`m^{ten}` Element beginnt und bei dem
:math:`(n-1)^{ten}` Element aufhört, wie unten gezeigt wird.

::

    sage: L = [factor(n) for n in range(1, 20)]
    sage: L[4:9]
    [5, 2 * 3, 7, 2^3, 3^2]
    sage: L[:4]
    [1, 2, 3, 2^2]
    sage: L[14:4]
    []
    sage: L[14:]
    [3 * 5, 2^4, 17, 2 * 3^2, 19]

Tupel sind ähnlich wie Listen, außer dass sie unveränderbar sind, was
bedeutet dass sie, sobald sie erzeugt wurden, nicht mehr verändert werden
können.

::

    sage: v = (1,2,3,4); v
    (1, 2, 3, 4)
    sage: type(v)
    <... 'tuple'>
    sage: v[1] = 5
    Traceback (most recent call last):
    ...
    TypeError: 'tuple' object does not support item assignment

Folgen sind ein dritter an Listen angelehnter Sage-Datentyp. Anders
als Listen und Tupel, sind Folgen kein gewöhnlicher Python-Datentyp.
Standardmäßig sind Folgen veränderbar, mit der
``Sequence``-Klassenmethode ``set_immutable`` können sie auf unveränderbar
gestellt werden, wie das folgende Beispiel zeigt. Alle Elemente einer
Folge haben einen gemeinsamen Obertyp, der das Folgenuniversum genannt wird.

::

    sage: v = Sequence([1,2,3,4/5])
    sage: v
    [1, 2, 3, 4/5]
    sage: type(v)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: type(v[1])
    <class 'sage.rings.rational.Rational'>
    sage: v.universe()
    Rational Field
    sage: v.is_immutable()
    False
    sage: v.set_immutable()
    sage: v[0] = 3
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Folgen sind von Listen abgeleitet und können überall dort verwendet werden, wo auch
Listen benutzt werden können.

::

    sage: v = Sequence([1,2,3,4/5])
    sage: isinstance(v, list)
    True
    sage: list(v)
    [1, 2, 3, 4/5]
    sage: type(list(v))
    <... 'list'>

Ein weiteres Beispiel von unveränderbaren Folgen sind Basen von
Vektorräumen. Es ist wichtig, dass sie nicht verändert werden können.

::

    sage: V = QQ^3; B = V.basis(); B
    [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    sage: type(B)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: B[0] = B[1]
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.
    sage: B.universe()
    Vector space of dimension 3 over Rational Field

Dictionaries
============

Ein Dictionary (manchmal auch assoziativer Array genannt) ist eine
Abbildung von 'hashbaren' Objekten (z.B. Strings, Zahlen und Tupel;
Lesen Sie die Python documentation
http://docs.python.org/3/tutorial/datastructures.html und
https://docs.python.org/3/library/stdtypes.html#typesmapping für weitere Details) zu
beliebigen Objekten.

::

    sage: d = {1:5, 'sage':17, ZZ:GF(7)}
    sage: type(d)
    <... 'dict'>
    sage: list(d.keys())
    [1, 'sage', Integer Ring]
    sage: d['sage']
    17
    sage: d[ZZ]
    Finite Field of size 7
    sage: d[1]
    5

Der dritte "key" zeigt, dass Indizes eines Dictionaries kompliziert,
also beispielsweise der Ring der ganzen Zahlen, sein können.

Sie können das obige Dictionary auch in eine Liste mit den gleichen
Daten umwandeln:

.. link

::

    sage: list(d.items())
    [(1, 5), ('sage', 17), (Integer Ring, Finite Field of size 7)]

Eine häufig vorkommende Ausdrucksweise ist über einem Paar in einem
Dictionary zu iterieren:

::

    sage: d = {2:4, 3:9, 4:16}
    sage: [a*b for a, b in d.items()]
    [8, 27, 64]

Ein Dictionary ist ungeordnet, wie die letzte Ausgabe verdeutlicht.

Mengen
======

Python hat einen standardmäßigen Mengen-Datentyp. Sein Hauptmerkmal
ist, neben weiteren typischen Mengenoperationen, dass das Nachschlagen
ob ein Element zu der Menge gehört oder nicht, sehr schnell geht.

::

    sage: X = set([1,19,'a']);   Y = set([1,1,1, 2/3])
    sage: X   # random sort order
    {1, 19, 'a'}
    sage: X == set(['a', 1, 1, 19])
    True
    sage: Y
    {2/3, 1}
    sage: 'a' in X
    True
    sage: 'a' in Y
    False
    sage: X.intersection(Y)
    {1}

Sage besitzt auch einen eigenen Mengen-Datentyp, welcher (manchmal)
mit Hilfe des standardmäßigen Python-Mengen-Datentyps implementiert
ist, jedoch darüberhinaus manche Sage-spezifischen Funktionen
aufweist. Sie können eine Sage-Menge erzeugen indem Sie ``Set(...)``
verwenden. Zum Beispiel,

::

    sage: X = Set([1,19,'a']);   Y = Set([1,1,1, 2/3])
    sage: X   # random sort order
    {'a', 1, 19}
    sage: X == Set(['a', 1, 1, 19])
    True
    sage: Y
    {1, 2/3}
    sage: X.intersection(Y)
    {1}
    sage: print(latex(Y))
    \left\{1, \frac{2}{3}\right\}
    sage: Set(ZZ)
    Set of elements of Integer Ring

Iteratoren
==========

Iteratoren sind seit Version 2.2 ein Teil von Python und erweisen sich
in mathematischen Anwendungen als besonders nützlich. Wir geben hier
ein paar Beispiele an; Lesen Sie [PyT]_ um weitere Details zu
erfahren. Wir erstellen einen Iterator über die Quadrate der
nichtnegativen ganzen Zahlen bis :math:`10000000`.

::

    sage: v = (n^2 for n in range(10000000))
    sage: next(v)
    0
    sage: next(v)
    1
    sage: next(v)
    4

Nun erzeugen wir einen Iterator über den Primzahlen der Form :math:`4p+1`
wobei auch :math:`p` prim ist und schauen uns die ersten Werte an.

::

    sage: w = (4*p + 1 for p in Primes() if is_prime(4*p+1))
    sage: w         # in the next line, 0xb0853d6c is a random 0x number
    <generator object at 0xb0853d6c>
    sage: next(w)
    13
    sage: next(w)
    29
    sage: next(w)
    53

Bestimmte Ringe, z. B. endliche Körper und die ganzen Zahlen, haben
zugehörige Iteratoren:


::

    sage: [x for x in GF(7)]
    [0, 1, 2, 3, 4, 5, 6]
    sage: W = ((x,y) for x in ZZ for y in ZZ)
    sage: next(W)
    (0, 0)
    sage: next(W)
    (0, 1)
    sage: next(W)
    (0, -1)

Schleifen, Funktionen, Kontrollstrukturen und Vergleiche
========================================================

Wir haben schon ein paar Beispiele gesehen in denen die
``for``-Schleife üblicherweise Verwendung findet. In Python hat eine
``for``-Schleife eine eingerückte Struktur, wie hier:

::

    >>> for i in range(5):
    ...     print(i)
    ...
    0
    1
    2
    3
    4

Beachten Sie den Doppelpunkt am Ende der for-Anweisung (dort befindet
sich kein "do" oder "od" wie in GAP oder Maple) und die Einrückung
vor dem Schleifenrumpf, dem ``print(i)``. Diese Einrückung ist
wichtig. In Sage wird die Einrückung automatisch hinzugefügt wenn Sie
nach einem ":" die ``enter``-Taste drücken, wie etwa im Folgenden
Beispiel.

::

    sage: for i in range(5):
    ....:     print(i)  # now hit enter twice
    ....:
    0
    1
    2
    3
    4


Das Symbol ``=`` wird bei Zuweisungen verwendet.
Das Symbol ``==`` wird verwendet um Gleichheit zu testen:

::

    sage: for i in range(15):
    ....:     if gcd(i,15) == 1:
    ....:         print(i)
    1
    2
    4
    7
    8
    11
    13
    14

Behalten Sie im Gedächtnis, dass die Block-Struktur von ``if``,
``for`` und ``while`` Ausdrücken durch die Einrückung bestimmt
wird:

::

    sage: def legendre(a,p):
    ....:     is_sqr_modp=-1
    ....:     for i in range(p):
    ....:         if a % p == i^2 % p:
    ....:             is_sqr_modp=1
    ....:     return is_sqr_modp

    sage: legendre(2,7)
    1
    sage: legendre(3,7)
    -1

Natürlich ist dies keine effiziente Implementierung des
Legendre-Symbols! Dies soll nur bestimmte Aspekte won Python/Sage
verdeutlichen. Die Funktion {kronecker}, welche zu Sage gehört,
berechnet das Legendre-Symbol effizient mittels eines Aufrufs von
PARIs C-Bibliothek.

Schließlich merken wir an, dass Vergleiche wie ``==``, ``!=``, ``<=``,
``>=``, ``>``, ``<`` von zwei Zahlen automatisch beide Zahlen in den
gleichen Typ konvertieren, falls dies möglich ist:

::

    sage: 2 < 3.1; 3.1 <= 1
    True
    False
    sage: 2/3 < 3/2;   3/2 < 3/1
    True
    True

Nutzen Sie bool für symbolische Ungleichungen:

::

    sage: x < x + 1
    x < x + 1
    sage: bool(x < x + 1)
    True

Beim Vergleichen von Objekten unterschiedlichen Typs versucht Sage in
den meisten Fällen eine kanonische Umwandlung beider Objekte in einen
gemeinsamen Typ zu finden. Falls erfolgreich wird der Vergleich auf den
umgewandelten Objekten durchgeführt; Falls nicht erfolgreich werden
die Objekte als ungleich angesehen. Um zu Testen, ob zwei Variablen
auf das gleiche Objekt zeigen, verwenden Sie ``is``. Zum Beispiel:

::

    sage: 1 is 2/2
    False
    sage: 1 is 1
    False
    sage: 1 == 2/2
    True

In den folgenden zwei Zeilen ist der erste Gleichheitstest ``False``,
da es keinen kanonischen Morphismus :math:`\QQ\ \to \GF{5}` gibt,
also gibt es keine kanonische Möglichkeit die  :math:`1` in :math:`\GF{5}`
mit der :math:`1 \in \QQ` zu vergleichen. Im Gegensatz dazu gibt es
eine kanonische Abbildung :math:`\ZZ \to \GF{5}`, also ist der zweite
Gleichheitstest ``True``. Beachten Sie auch, dass die Reihenfolge
keine Rolle spielt.

::

    sage: GF(5)(1) == QQ(1); QQ(1) == GF(5)(1)
    False
    False
    sage: GF(5)(1) == ZZ(1); ZZ(1) == GF(5)(1)
    True
    True
    sage: ZZ(1) == QQ(1)
    True

WARNUNG: Vergleiche in Sage sind restriktiver als in Magma, welches
die :math:`1 \in \GF{5}` gleich der :math:`1 \in \QQ` festlegt.

::

    sage: magma('GF(5)!1 eq Rationals()!1')            # optional - magma
    true

Profiling
=========

Autor des Abschnitts: Martin Albrecht (malb@informatik.uni-bremen.de)

    "Premature optimization is the root of all evil." - Donald Knuth


Manchmal ist es nützlich nach Engstellen im Code zu suchen, um zu
verstehen welche Abschnitte die meiste Berechnungszeit beanspruchen;
dies kann ein guter Hinweis darauf sein, welche Teile optimiert werden
sollten. Python, und daher auch Sage, stellen mehrere "Profiling" -- so
wird dieser Prozess genannt -- Optionen zur Verfügung.

Am einfachsten zu Benutzen ist das ``prun``-Kommando in der
interaktiven Shell. Es gibt eine Zusammenfassung zurück, die
beschreibt welche Funktionen wie viel Berechnungszeit veranschlagt haben.
Um die (zu diesem Zeitpunkt langsame) Matrixmultiplikation über
endlichen Körpern zu Profilieren, geben Sie z.B. folgendes ein:

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])

.. skip

::

    sage: %prun B = A*A
           32893 function calls in 1.100 CPU seconds

    Ordered by: internal time

    ncalls tottime percall cumtime percall filename:lineno(function)
     12127  0.160   0.000   0.160  0.000 :0(isinstance)
      2000  0.150   0.000   0.280  0.000 matrix.py:2235(__getitem__)
      1000  0.120   0.000   0.370  0.000 finite_field_element.py:392(__mul__)
      1903  0.120   0.000   0.200  0.000 finite_field_element.py:47(__init__)
      1900  0.090   0.000   0.220  0.000 finite_field_element.py:376(__compat)
       900  0.080   0.000   0.260  0.000 finite_field_element.py:380(__add__)
         1  0.070   0.070   1.100  1.100 matrix.py:864(__mul__)
      2105  0.070   0.000   0.070  0.000 matrix.py:282(ncols)
      ...

Hier ist ``ncalls`` die Anzahl der Aufrufe, ``tottime`` ist die
Gesamtzeit, die für die Funktion verwendet wurde (ausgenommen der
Zeit, die für Unterfunktionsaufrufe verwendet wurde), ``percall`` ist
der Quotient von  ``tottime`` geteilt durch ``ncalls``. ``cumtime``
ist die Gesamtzeit, die für diese Funktion und alle
Unterfunktionsaufrufe (d.h., vom Aufruf bis zum Ende) verwendet
wurde, ``percall`` ist der Quotient von ``cumtime`` geteilt durch die
Zeit elementarer Funktionsaufrufe, und ``filename:lineno(function)``
stellt die entsprechenden Daten jeder Funktion zur Verfügung. Die
Daumenregel ist hier: Je höher die Funktion in dieser Liste steht,
desto teurer ist sie. Also ist sie interessanter für Optimierungen.

Wie sonst auch stellt ``prun?`` Details zur Benutzung des Profilers
und zum Verstehen seines Outputs zur Verfügung.

Die Profilierungsdaten können auch in ein Objekt geschrieben werden um
eine weitere Untersuchung zu ermöglichen:

.. skip

::

    sage: %prun -r A*A
    sage: stats = _
    sage: stats?

Beachten Sie: das Eingeben von ``stats = prun -r A\*A`` erzeugt eine
Syntaxfehlermeldung, da prun ein IPython-Shell-Kommando ist und keine
reguläre Funktion.

Um eine schöne graphische Repräsentation der Profilerdaten zu
erhalten, können Sie den "hotshot-Profiler", ein kleines Skript
genannt ``hotshot2cachetree`` und das Programm ``kcachegrind`` (nur
für Unix) benutzen. Hier ist das gleiche Beispiel mit dem "hotshot-Profiler":

.. skip

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])
    sage: import hotshot
    sage: filename = "pythongrind.prof"
    sage: prof = hotshot.Profile(filename, lineevents=1)

.. skip

::

    sage: prof.run("A*A")
    <hotshot.Profile instance at 0x414c11ec>
    sage: prof.close()

Dies führt zu einer Datei ``pythongrind.prof`` in aktuellen
Datenverzeichnis. Diese kann nun zur Visualisierung in das
cachegrind-Format konvertiert werden.

Tippen Sie in einer System-Shell:

.. skip

::

    hotshot2calltree -o cachegrind.out.42 pythongrind.prof

Die Ausgabedatei ``cachegrind.out.42`` kann nun mit ``kcachegrind``
untersucht werden. Bitte beachten Sie, dass die Namenskonvention
``cachegrind.out.XX`` erhalten bleiben muss.
