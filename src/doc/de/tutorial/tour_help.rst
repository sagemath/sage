.. _chapter-help:

Hilfe
=====

Sage hat eine umfassende eingebaute Dokumentation, auf die zugegriffen
werden kann, indem der Name der Funktion oder Konstanten (zum
Beispiel) gefolgt von einem Fragezeichen eingegeben wird:

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring:

        The tangent function

        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring:

        The natural logarithm of the real number 2.

        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 32-bit
            0.69314718055994530941723212145817656807   # 64-bit
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <... 'function'>
    Definition:  sudoku(A)
    Docstring:

        Solve the 9x9 Sudoku puzzle defined by the matrix A.

        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

Sage stellt auch eine 'Tab-Vervollständigung' zur Verfügung: Schreiben Sie die
ersten Buchstaben einer Funktion und drücken Sie danach die Tabulatortaste.
Wenn Sie zum Beispiel ``ta`` gefolgt von ``TAB`` eingeben, wird Sage
``tachyon, tan, tanh, taylor`` ausgeben. Dies ist eine gute
Möglichkeit den Namen von Funktionen und anderen Strukturen in Sage herauszufinden.


.. _section-functions:

Funktionen, Einrückungen, und Zählen
====================================

Um in Sage eine neue Funktion zu definieren, können Sie den ``def``
Befehl und einen Doppelpunkt nach der Liste der Variablennamen
benutzen. Zum Beispiel:

::

    sage: def is_even(n):
    ....:     return n%2 == 0
    sage: is_even(2)
    True
    sage: is_even(3)
    False

Anmerkung: Abhängig von der Version des Tutorials, das Sie gerade lesen,
sehen Sie vielleicht drei Punkte ``....:`` in der zweiten Zeile dieses
Beispiels. Tippen Sie diese nicht; sie sind nur da um zu
verdeutlichen, dass der Code eingerückt ist. Wann immer dies der Fall
ist, drücken Sie [Return/Enter] einmal am Ende des Blocks um eine
Leerzeile einzufügen und die Funktionsdefinition zu beenden.

Sie bestimmen den Typ ihrer Eingabeargumente nicht. Sie können mehrere
Argumente festlegen, jedes davon kann einen optionalen Standardwert haben.
Zum Beispiel wird in der Funktion unterhalb standardmäßig der Wert
``divisor=2`` benutzt, falls ``divisor`` nicht angegeben wurde.

::

    sage: def is_divisible_by(number, divisor=2):
    ....:     return number%divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False

Sie können auch ein oder mehrere Eingabeargumente explizit angeben
wenn Sie die Funktion aufrufen; wenn Sie die Eingaben explizit
angeben, können Sie dies in beliebiger Reihenfolge tun:

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True

In Python werden Codeblöcke nicht mit geschweiften Klammern oder
"begin-" und "end-Blöcken" kenntlich gemacht. Stattdessen werden
Codeblöcke durch Einrückungen bestimmt, welche exakt zusammenpassen
müssen. Zum Beispiel ist das Folgende ein Syntaxfehler, da die
``return`` Anweisung nicht genauso weit eingerückt ist wie die anderen
Zeilen zuvor.

.. skip

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:    return v
    Syntax Error:
           return v

Wenn Sie die Einrückung korrigieren, funktioniert die Funktion:

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:     return v
    sage: even(10)
    [4, 6, 8]

Semikola sind an den Zeilenenden nicht notwendig; sie können jedoch
mehrere Anweisungen, mit Semikola getrennt, in eine Zeile schreiben:

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

Falls Sie möchten, dass sich eine einzelne Codezeile über mehrere
Zeilen erstreckt, können Sie einen terminierenden Backslash verwenden:

::

    sage: 2 + \
    ....:    3
    5

In Sage können Sie zählen indem Sie über einen Zahlenbereich
iterieren. Zum Beispiel ist nächste Zeile unterhalb gleichwertig zu
``for(i=0; i<3; i++)`` in C++ oder Java:

::

    sage: for i in range(3):
    ....:     print(i)
    0
    1
    2

Die nächste Zeile unterhalb ist gleichwertig zu ``for(i=2;i<5;i++)``.

::

    sage: for i in range(2,5):
    ....:     print(i)
    2
    3
    4

Das dritte Argument bestimmt die Schrittweite, also ist das Folgende
gleichwertig zu
``for(i=1;i<6;i+=2)``.

::

    sage: for i in range(1,6,2):
    ....:     print(i)
    1
    3
    5

Oft will man eine schöne Tabelle erstellen, um die mit Sage
berechneten Zahlen auszugeben. Eine einfache Möglichkeit dies zu tun
ist String-Formatierung zu verwenden. Unten erstellen wir drei Spalten,
jede genau 6 Zeichen breit, und erzeugen somit eine Tabelle mit
Quadrat- und Kubikzahlen.

::

    sage: for i in range(5):
    ....:     print('%6s %6s %6s' % (i, i^2, i^3))
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64

Die elementarste Datenstruktur in Sage ist die Liste. Sie ist -- wie
der Name schon sagt -- nichts anderes als eine Liste beliebiger
Objekte. Hier ist ein Beispiel::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

Listenindizierung beginnt, wie in vielen Programmiersprachen, bei 0.

.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)

Benutzen Sie ``len(v)`` um die Länge von ``v`` zu erhalten, benutzen
Sie ``v.append(obj)`` um ein neues Objekt an das Ende von ``v``
anzuhängen, und benutzen Sie ``del v[i]`` um den :math:`i^{ten}`
Eintrag von ``v`` zu löschen:

.. link

::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]

Eine weitere wichtige Datenstruktur ist das Dictionary (oder
assoziatives Array). Dies funktioniert wie eine Liste, außer dass
es mit fast jedem Objekt indiziert werden kann (die Indizes müssen
jedoch unveränderbar sein):

::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi

Sie können auch neue Datentypen definieren, indem Sie Klassen
verwenden. Mathematische Objekte mit Klassen zusammenzufassen ist eine
mächtige Technik, die dabei helfen kann Sage-Programme zu vereinfachen
und zu organisieren. Unten definieren wir eine Klasse, welche die Liste
der geraden Zahlen bis *n* darstellt;
Sie wird von dem Standard-Typ ``list`` abgeleitet.

::

    sage: class Evens(list):
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:         list.__init__(self, range(2, n+1, 2))
    ....:     def __repr__(self):
    ....:         return "Even positive numbers up to n."

Die ``__init__`` Methode wird aufgerufen um das Objekt zu
initialisieren, wenn es erzeugt wird; die ``__repr__`` Method gibt
einen Objekt-String aus. Wir rufen die Listen-Konstruktor-Methode in
der zweite Zeile der ``__init__`` Methode. Ein Objekt der Klasse
``Evens`` erzeugen wir wie folgt:

.. link

::

    sage: e = Evens(10)
    sage: e
    Even positive numbers up to n.

Beachten Sie, dass die Ausgabe von ``e`` die ``__repr__`` Methode
verwendet, die wir definiert haben. Um die eigentliche Liste
sehen zu können, benutzen wir die ``list``-Funktion:

.. link

::

    sage: list(e)
    [2, 4, 6, 8, 10]

Wir können auch das ``n`` Attribut verwenden oder ``e`` wie eine Liste
behandeln.

.. link

::

    sage: e.n
    10
    sage: e[2]
    6
