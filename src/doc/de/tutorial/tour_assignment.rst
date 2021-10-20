Zuweisung, Gleichheit und Arithmetik
====================================

Bis auf wenige Ausnahmen benutzt Sage die Programmiersprache Python,
deshalb werden Ihnen die meisten einführenden Bücher über Python dabei
helfen Sage zu lernen.

Sage benutzt ``=`` für Zuweisungen und ``==``, ``<=``, ``>=``,
``<`` und ``>`` für Vergleiche.

::

    sage: a = 5
    sage: a
    5
    sage: 2 == 2
    True
    sage: 2 == 3
    False
    sage: 2 < 3
    True
    sage: a == 5
    True

Sage unterstützt alle grundlegenden mathematischen Operationen:

::

    sage: 2**3    #  ** bedeutet hoch
    8
    sage: 2^3     #  ^ ist ein Synonym für ** (anders als in Python)
    8
    sage: 10 % 3  #  für ganzzahlige Argumente bedeutet % mod, d.h. Restbildung
    1
    sage: 10/4
    5/2
    sage: 10//4   #  für ganzzahlige Argumente gibt // den \
    ....:         #  ganzzahligen Quotienten zurück
    2
    sage: 4 * (10 // 4) + 10 % 4 == 10
    True
    sage: 3^2*4 + 2%5
    38

Die Berechnung eines Ausdrucks wie ``3^2*4 + 2%5`` hängt von der
Reihenfolge ab, in der die Operationen ausgeführt werden. Dies wird in
der "Operatorrangfolge-Tabelle" in :ref:`section-precedence` festgelegt.

Sage stellt auch viele bekannte mathematische Funktionen zur
Verfügung; hier sind nur ein paar Beispiele

::

    sage: sqrt(3.4)
    1.84390889145858
    sage: sin(5.135)
    -0.912021158525540
    sage: sin(pi/3)
    1/2*sqrt(3)

Wie das letzte Beispiel zeigt, geben manche mathematische Ausdrücke
'exakte' Werte anstelle von numerischen Approximationen zurück. Um
eine numerische Approximation zu bekommen, können Sie entweder die
Funktion ``n`` oder die Methode ``n`` verwenden (beide haben auch
den längeren Namen, ``numerical_approx``, und die Funktion ``N`` ist
die gleiche wie ``n``). Diese nehmen auch die optionalen Argumente
``prec``, welches die gewünschte Anzahl von Bit an Genauigkeit ist und
``digits``, welches die gewünschte Anzahl Dezimalstellen an Genauigkeit
ist, entgegen; der Standardwert ist 53 Bit Genauigkeit.

::

    sage: exp(2)
    e^2
    sage: n(exp(2))
    7.38905609893065
    sage: sqrt(pi).numerical_approx()
    1.77245385090552
    sage: sin(10).n(digits=5)
    -0.54402
    sage: N(sin(10),digits=10)
    -0.5440211109
    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

Python ist dynamisch typisiert, also ist dem Wert, auf den jede Variable
weist, ein Typ zugeordnet; jedoch darf eine Variable Werte eines
beliebigen Python-Typs innerhalb eines Sichtbarkeitsbereich aufnehmen.

::

    sage: a = 5   # a ist eine ganze Zahl
    sage: type(a)
    <class 'sage.rings.integer.Integer'>
    sage: a = 5/3  # jetzt ist a eine rationale Zahl
    sage: type(a)
    <class 'sage.rings.rational.Rational'>
    sage: a = 'hello'  # jetzt ist a ein String
    sage: type(a)
    <... 'str'>

Die Programmiersprache C, welche statisch typisiert ist, unterscheidet
sich hierzu stark; eine Variable, die dazu deklariert ist eine Ganzzahl (int)
aufzunehmen, kann in ihrem Sichtbarkeitsbereich auch nur ganze Zahlen aufnehmen.
