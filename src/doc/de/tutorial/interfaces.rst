**************
Schnittstellen
**************

Ein zentraler Aspekt von Sage ist, dass es Berechnungen mit Objekten
vieler verschiedener Computer Algebra Systeme unter einem Dach durch eine
einheitliche Schnittstelle und Programmiersprache vereinigt.

Die ``console`` und ``interact`` Methoden einer Schnittstelle unterstützen
viele verschiedene Dinge. Zum Beispiel, anhand von GAP:


#. ``gap.console()``: Öffnet die GAP Konsole und übergibt GAP die
   Kontrolle. Hier ist Sage nichts weiter als ein praktischer
   Programmstarter, ähnlich einer Linux-Bash-Konsole.

#. ``gap.interact()``: Ist eine einfache Art mit einer GAP Instanz
   zu interagieren, die Sage Objekte enthalten kann. Sie können Sage
   Objekte in diese GAP Sitzung importieren (sogar von der interaktiven
   Schnittstelle aus), usw.

.. index: PARI; GP

GP/PARI
=======

PARI ist ein kompaktes, sehr ausgereiftes, stark optimiertes C-Programm,
dessen primärer Fokus Zahlentheorie ist. Es gibt zwei sehr verschiedene
Schnittstellen, die Sie in Sage nutzen können:


-  ``gp`` -- Der **PARI** Interpreter und

-  ``pari`` -- Die **PARI** C-Bibliothek.


Die folgenden Zeilen zum Beispiel sind zwei Wege, genau das gleiche zu
tun. Sie sehen identisch aus, aber die Ausgabe ist verschieden, und was
hinter den Kulissen passiert ist völlig unterschiedlich.

::

    sage: gp('znprimroot(10007)')
    Mod(5, 10007)
    sage: pari('znprimroot(10007)')
    Mod(5, 10007)

Im ersten Fall wird eine separate Kopie des GP-Interpreters als Server
gestartet, die Zeichenkette ``'znprimroot(10007)'`` übergeben,
von GP ausgewertet und das Ergebnis wird einer Variable in GP zugewiesen
(was Platz im Speicher des GP-Unterprozesses benötigt, der nicht wieder
freigegeben wird). Dann wird der Wert der Variablen erst angezeigt.
Im zweiten Fall wird kein separates Programm gestartet, stattdessen
wird die Zeichenkette ``'znprimroot(10007)'`` von einer bestimmten
PARI-C-Bibliotheksfunktion ausgewertet. Das Ergebnis wird im Speicher
von Python gehalten, welcher freigegeben wird wenn die Variable nicht
mehr referenziert wird. Die Objekte haben außerdem verschiedene Typen:

::

    sage: type(gp('znprimroot(10007)'))
    <class 'sage.interfaces.gp.GpElement'>
    sage: type(pari('znprimroot(10007)'))
    <class 'cypari2.gen.Gen'>

Welche Variante sollten Sie also nutzen? Das kommt darauf an was
Sie tun. Die GP-Schnittstelle kann alles was ein normales
GP/PARI-Konsolenprogramm könnte, da es das Programm
startet. Genauergesagt könnten Sie komplizierte PARI-Programme laden
und laufen lassen. Im Gegensatz dazu ist die PARI-Schnittstelle
(mittels C-Bibliothek) wesentlich restriktiver. Zuerst einmal sind
nicht alle Unterfunktionen implementiert. Außerdem
wird relativ viel Quellcode nicht in der PARI-Schnittstelle funktionieren,
z.B. numerisches Integrieren. Abgesehen davon ist die PARI-Schnittstelle
wesentlich schneller und robuster als die von GP.

(Wenn der GP-Schnittstelle der Speicher ausgeht beim Auswerten einer
Zeile, wird sie automatisch und unbemerkt den Speicherbereich
verdoppeln und das Auswerten erneut versuchen. Dadurch wird Ihre
Berechnung nicht abstürzen, falls Sie den benötigen Speicher falsch
eingeschätzt haben. Das ist eine  hilfreiche Erweiterung, die der
gewöhnliche GP-Interpreter nicht bietet. Die PARI-C-Bibliothek
hingegen kopiert jedes neue Objekt sofort vom PARI-Stack, daher wird
der Stapel nicht größer. Allerdings muss jedes Objekt kleiner als 100
MB sein, da ansonsten der Stapel "überläuft", wenn das Objekt erstellt
wird. Dieses zusätzliche Kopieren erzeugt allerdings ein wenig
Leistungseinbußen.)

Zusammengefasst nutzt Sage also die PARI-C-Bibliothek um Funktionalitäten
eines GP/PARI-Interpreters bereitzustellen, allerdings mit einer anderen
komplizierten Speicherverwaltung und der Programmiersprache Python.

Zuerst erstellen wir eine PARI-Liste aus einer Python-Liste.

::

    sage: v = pari([1,2,3,4,5])
    sage: v
    [1, 2, 3, 4, 5]
    sage: type(v)
    <class 'cypari2.gen.Gen'>

Jedes PARI-Objekt ist vom Typ ``Gen``. Den PARI Typ des vorliegenden
Objekts können Sie mit der ``type`` Unterfunktion herausfinden.

.. link

::

    sage: v.type()
    't_VEC'

Um eine elliptische Kurve in PARI zu erstellen geben Sie
``ellinit([1,2,3,4,5])`` ein. Bei Sage ist es ähnlich, nur
dass ``ellinit`` eine Methode ist, die von jedem PARI-Objekt
aus aufgerufen werden kann, z.B. unser ``t_VEC`` `v`.

.. link

::

    sage: e = v.ellinit()
    sage: e.type()
    't_VEC'
    sage: pari(e)[:13]
    [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351]

Jetzt haben wir eine elliptische Kurve als Objekt und können einige
Dinge mit ihr berechnen.

.. link

::

    sage: e.elltors()
    [1, [], []]
    sage: e.ellglobalred()
    [10351, [1, -1, 0, -1], 1, [11, 1; 941, 1], [[1, 5, 0, 1], [1, 5, 0, 1]]]
    sage: f = e.ellchangecurve([1,-1,0,-1])
    sage: f[:5]
    [1, -1, 0, 4, 3]

.. index: GAP

.. _section-gap:

GAP
===

Sage enthält ausserdem GAP für diskrete Mathematik, insbesondere
Gruppentheorie.

Hier ist ein Beispiel für GAP's ``IdGroup``-Funktion, die die optionale kleine
Gruppen Datenbank benötigt, die separat installiert werden muss, wie unten beschrieben.

::

    sage: G = gap('Group((1,2,3)(4,5), (3,4))')
    sage: G
    Group( [ (1,2,3)(4,5), (3,4) ] )
    sage: G.Center()
    Group( () )
    sage: G.IdGroup()
    [ 120, 34 ]
    sage: G.Order()
    120

Wir können die gleiche Berechnung in Sage durchführen ohne vorher explizit
die GAP-Schnittstelle aufzurufen:

::

    sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
    sage: G.center()
    Subgroup generated by [()] of (Permutation Group with generators [(3,4), (1,2,3)(4,5)])
    sage: G.group_id()
    [120, 34]
    sage: n = G.order(); n
    120

Nach Installation eine optionale Sage-Pakete mit folgendem Befehl
sind weitere GAP-Funktionen verfügbar::

    sage -i gap_packages


Singular
========

Singular bietet eine sehr gute, ausgereifte Bibliothek für Gröbnerbasen,
größte gemeinsame Teiler von mehrdimensionalen Polynomen, Basen von
Riemann-Roch Räumen einer planaren Kurve und Faktorisierungen unter anderem.
Wir zeigen hier die Faktorisierung mehrdimensionaler Polynome mit
Sages Singular-Schnittstelle (ohne die ``....:``):

::

    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: R1
    polynomial ring, over a field, global ordering
    // coefficients: QQ...
    // number of vars : 2
    //        block   1 : ordering dp
    //                  : names    x y
    //        block   2 : ordering C
    sage: f = singular('9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 +'
    ....:     '9*x^6*y^4 + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 -'
    ....:     '9*x^12*y^3 - 18*x^13*y^2 + 9*x^16')

Wir haben also das Polynom :math:`f` definiert, nun geben wir es aus und faktorisieren es.

.. link

::

    sage: f
    9*x^16-18*x^13*y^2-9*x^12*y^3+9*x^10*y^4-18*x^11*y^2+36*x^8*y^4+18*x^7*y^5-18*x^5*y^6+9*x^6*y^4-18*x^3*y^6-9*x^2*y^7+9*y^8
    sage: f.parent()
    Singular
    sage: F = f.factorize(); F
    [1]:
       _[1]=9
       _[2]=x^6-2*x^3*y^2-x^2*y^3+y^4
       _[3]=-x^5+y^2
    [2]:
       1,1,2
    sage: F[1][2]
    x^6-2*x^3*y^2-x^2*y^3+y^4

Genau wie im GAP Beispiel in :ref:`section-gap`, können wir diese Faktorisierung
berechnen ohne explizit die Singular-Schnittstelle zu nutzen.
(Dennoch nutzt Sage im Hintergrund die Singular-Schnittstelle für die Berechnung.)
Bitte geben Sie ein ohne ``....:``:

::

    sage: x, y = QQ['x, y'].gens()
    sage: f = (9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 + 9*x^6*y^4
    ....:     + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 - 9*x^12*y^3
    ....:     - 18*x^13*y^2 + 9*x^16)
    sage: factor(f)
    (9) * (-x^5 + y^2)^2 * (x^6 - 2*x^3*y^2 - x^2*y^3 + y^4)

.. _section-maxima:

Maxima
======

Das in Lisp implementierte Maxima ist ein Teil von Sage. Hingegen wird
das gnuplot-Paket (welches Maxima standardmäßig zum plotten nutzt) als
optionales Sage-Paket angeboten. Neben anderen Dingen rechnet Maxima
mit Symbolen. Maxima integriert und differenziert Funktionen
symbolisch, löst gewöhnliche Differentialgleichungen ersten Grades
sowie viele lineare Differentialgleichungen zweiten Grades und besitzt
eine Methode zur Laplace Transformation linearer
Differentialgleichungen von beliebigem Grad. Maxima kennt eine große
Zahl spezieller Funktionen, plottet mittels gnuplot und hat Methoden,
um Polynomgleichungen oder Matrizen zu lösen oder zu verändern
(z.B. Zeilenelimination oder Eigenwerte und Eigenvektoren berechnen).

Wir zeigen die Sage/Maxima Schnittstelle, indem wir die Matrix konstruieren,
deren :math:`i,j` Eintrag gerade :math:`i/j` ist, für :math:`i,j=1,\ldots,4`.

::

    sage: f = maxima.eval('ij_entry[i,j] := i/j')
    sage: A = maxima('genmatrix(ij_entry,4,4)'); A
    matrix([1,1/2,1/3,1/4],[2,1,2/3,1/2],[3,3/2,1,3/4],[4,2,4/3,1])
    sage: A.determinant()
    0
    sage: A.echelon()
    matrix([1,1/2,1/3,1/4],[0,0,0,0],[0,0,0,0],[0,0,0,0])
    sage: A.eigenvalues()
    [[0,4],[3,1]]
    sage: A.eigenvectors().sage()
    [[[0, 4], [3, 1]], [[[1, 0, 0, -4], [0, 1, 0, -2], [0, 0, 1, -4/3]], [[1, 2, 3, 4]]]]

Hier ein anderes Beispiel:

::

    sage: A = maxima("matrix ([1, 0, 0], [1, -1, 0], [1, 3, -2])")
    sage: eigA = A.eigenvectors()
    sage: V = VectorSpace(QQ,3)
    sage: eigA
    [[[-2,-1,1],[1,1,1]],[[[0,0,1]],[[0,1,3]],[[1,1/2,5/6]]]]
    sage: v1 = V(sage_eval(repr(eigA[1][0][0]))); lambda1 = eigA[0][0][0]
    sage: v2 = V(sage_eval(repr(eigA[1][1][0]))); lambda2 = eigA[0][0][1]
    sage: v3 = V(sage_eval(repr(eigA[1][2][0]))); lambda3 = eigA[0][0][2]

    sage: M = MatrixSpace(QQ,3,3)
    sage: AA = M([[1,0,0],[1, - 1,0],[1,3, - 2]])
    sage: b1 = v1.base_ring()
    sage: AA*v1 == b1(lambda1)*v1
    True
    sage: b2 = v2.base_ring()
    sage: AA*v2 == b2(lambda2)*v2
    True
    sage: b3 = v3.base_ring()
    sage: AA*v3 == b3(lambda3)*v3
    True

Zuletzt noch ein Beispiel wie man Sage zum Plotten mittels
``openmath`` nutzt. Einige von ihnen wurden (verändert) aus dem Maxima
Benutzerhandbuch entnommen.

Ein 2D-Plot verschiedener Funktionen (ohne ``....:`` eingeben):

::

    sage: maxima.plot2d('[cos(7*x),cos(23*x)^4,sin(13*x)^3]','[x,0,1]', # not tested
    ....:     '[plot_format,openmath]')

Ein "live" 3D-Plot, den man mit der Maus bewegen kann:

::

    sage: maxima.plot3d ("2^(-u^2 + v^2)", "[u, -3, 3]", "[v, -2, 2]", # not tested
    ....:     '[plot_format, openmath]')
    sage: maxima.plot3d("atan(-x^2 + y^3/4)", "[x, -4, 4]", "[y, -4, 4]", # not tested
    ....:     "[grid, 50, 50]",'[plot_format, openmath]')

Der nächste Plot ist das berühmte Möbiusband:

::

    sage: maxima.plot3d("[cos(x)*(3 + y*cos(x/2)), sin(x)*(3 + y*cos(x/2)), y*sin(x/2)]", # not tested
    ....:     "[x, -4, 4]", "[y, -4, 4]",
    ....:     '[plot_format, openmath]')

Und der letzte ist die berühmte Kleinsche Flasche:

::

    sage: _ = maxima("expr_1: 5*cos(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0) - 10.0")
    sage: _ = maxima("expr_2: -5*sin(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0)")
    sage: _ = maxima("expr_3: 5*(-sin(x/2)*cos(y) + cos(x/2)*sin(2*y))")
    sage: maxima.plot3d ("[expr_1, expr_2, expr_3]", "[x, -%pi, %pi]", # not tested
    ....:     "[y, -%pi, %pi]", "['grid, 40, 40]",
    ....:     '[plot_format, openmath]')

