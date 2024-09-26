.. _chapter-interactive_shell:

*****************************
Die interaktive Kommandozeile
*****************************
In den meisten Teilen dieses Tutorials gehen wir davon aus, dass Sie
Sage mit dem ``sage``-Befehl starten. Dieser startet eine angepasste
Version der IPython Kommandozeile und lädt Funktionen und Klassen,
sodass sie in der Kommandozeile genutzt werden können. Weitere
Anpassungen können Sie in der Datei ``$SAGE_ROOT/ipythonrc``
vornehmen. Nach dem Start von Sage sehen Sie etwa folgendes:

.. skip

::

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘


    sage:

Um Sage zu beenden drücken Sie Strg-D oder geben Sie
``quit`` oder ``exit`` ein.

.. skip

::

    sage: quit
    Exiting Sage (CPU time 0m0.00s, Wall time 0m0.89s)

Unter "wall time" finden Sie die vergangene Echtzeit (der Uhr an Ihrer
Wand). Diese ist nötig, da die CPU Zeit Unterprozesse wie GAP oder
Singular nicht berücksichtigt.

(Vermeiden Sie es den Sage Prozess mit ``kill -9`` in der Konsole zu
beenden, da so möglicherweise Unterprozesse wie z.B. Maple-Prozesse
nicht beendet oder temporäre Dateien in ``$HOME/.sage/tmp`` nicht
gelöscht würden.)


Ihre Sage Sitzung
=================

Unter einer Sitzung verstehen wir die Ein- und Ausgaben von Sage vom
Starten bis zum Beenden. Sage speichert alle Eingaben mittels IPython. Wenn
Sie die interaktive Kommandozeile nutzen (im Gegensatz zur
Browser-Oberfläche "Notebook"), so können Sie jederzeit mittels
``%hist`` eine Liste aller bisher getätigten Eingaben sehen. Sie
können auch ``?`` eingeben, um mehr über IPython zu
erfahren. Z.B. "IPython unterstützt Zeilennummerierung ... sowie Ein-
und Ausgabezwischenspeicherung.  Alle Eingaben werden gespeichert und
können in Variablen abgerufen werden (neben der normalen
Pfeiltasten-Navigation). Die folgenden globalen Variablen existieren
immer (also bitte überschreiben Sie sie nicht!)":

::

      _:  letzte Eingabe (interaktive Kommandozeile und Browser-Oberfläche)
      __: vorletzte Eingabe (nur in der Kommandozeile)
      _oh: Liste aller Eingaben (nur in der Kommandozeile)

Hier ein Beispiel:

.. skip

::

    sage: factor(100)
     _1 = 2^2 * 5^2
    sage: kronecker_symbol(3,5)
     _2 = -1
    sage: %hist   # funktioniert nur in der Kommandozeile, nicht im Browser.
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    sage: _oh
     _4 = {1: 2^2 * 5^2, 2: -1}
    sage: _i1
     _5 = 'factor(ZZ(100))\n'
    sage: eval(_i1)
     _6 = 2^2 * 5^2
    sage: %hist
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    4: _oh
    5: _i1
    6: eval(_i1)
    7: %hist

Wir lassen die Zeilennummerierung im restlichen Tutorial sowie in der
weiteren Sage-Dokumentation weg. Sie können auch eine Liste von
Eingaben einer Sitzung in einem Makro für diese Sitzung speichern.

.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    1: E = EllipticCurve([1,2,3,4,5])
    2: M = ModularSymbols(37)
    3: %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).


.. skip

::

    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E = 5
    sage: M = None
    sage: em
    Executing Macro...
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field

Während Sie die interaktive Kommandozeile nutzen, können Sie jeden
UNIX-Kommandozeilenbefehl in Sage ausführen, indem Sie ihm ein
Ausrufezeichen ``!`` voranstellen. Zum Beispiel gibt

.. skip

::

    sage: !ls
    auto  example.sage glossary.tex  t  tmp  tut.log  tut.tex

den Inhalt des aktuellen Verzeichnisses aus.

In der ``PATH``-Variablen steht das Sage "bin" Verzeichnis vorne. Wenn
Sie also ``gp``, ``gap``, ``singular``, ``maxima``, usw. eingeben,
starten Sie die in Sage enthaltenden Versionen.

.. skip

::

    sage: !gp
    Reading GPRC: /etc/gprc ...Done.

                               GP/PARI CALCULATOR Version 2.2.11 (alpha)
                      i686 running linux (ix86/GMP-4.1.4 kernel) 32-bit version
    ...
    sage: !singular
                         SINGULAR                             /  Development
     A Computer Algebra System for Polynomial Computations   /   version 3-1-0
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   Mar 2009
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \

Ein- und Ausgaben loggen
========================

Die Sage Sitzung loggen bzw. speichern ist nicht das Gleiche (siehe
:ref:`section-save`). Um Eingaben (und optional auch Ausgaben) zu
loggen nutzen Sie den Befehl ``logstart``. Geben Sie ``logstart?`` ein
um weitere Informationen zu erhalten. Sie können diesen Befehl nutzen
um alle Eingaben und Ausgaben zu loggen, und diese sogar wiederholen
in einer zukünftigen Sitzung (indem Sie einfach die Log-Datei laden).

.. skip

::

    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: logstart setup
    Activating auto-logging. Current session state plus future input saved.
    Filename       : setup
    Mode           : backup
    Output logging : False
    Timestamping   : False
    State          : active
    sage: E = EllipticCurve([1,2,3,4,5]).minimal_model()
    sage: F = QQ^3
    sage: x,y = QQ['x,y'].gens()
    sage: G = E.gens()
    sage:
    Exiting Sage (CPU time 0m0.61s, Wall time 0m50.39s).
    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: load("setup")
    Loading log file <setup> one line at a time...
    Finished replaying log file <setup>
    sage: E
    Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 + 4*x + 3 over Rational
    Field
    sage: x*y
    x*y
    sage: G
    [(2 : 3 : 1)]

Wenn Sie Sage in der Linux KDE Konsole ``konsole`` verwenden, können
Sie Ihre Sitzung wie folgt speichern: Nachdem Sie Sage in ``konsole``
gestartet haben, wählen Sie "Einstellungen", dann "Verlauf...", dann
"auf unbegrenzt" setzen. Wenn Sie soweit sind Ihre Sitzung zu
speichern, wählen Sie "Bearbeiten" und dann "Verlauf speichern
unter..."  und geben einen Namen ein, um den Text ihrer Sitzung
auf dem Computer zu speichern. Nach dem Speichern der Datei können Sie
jene in einem Editor wie GNU Emacs öffnen und ausdrucken.


Einfügen ignoriert Eingabeaufforderungen
========================================

Stellen Sie sich vor, Sie lesen eine Sitzung von Sage oder Python
Berechnungen und  wollen sie in Sage kopieren, aber überall sind noch
die störenden ``>>>`` oder ``sage:``
Eingabeaufforderungen. Tatsächlich können Sie einfach die gewünschte
Stelle mit Eingabeaufforderungen in Sage einfügen. Der Sage Parser
wird standardmäßig die führenden ``>>>`` oder ``sage:``
Eingabeaufforderungen entfernen bevor er es an Python weitergibt. Zum
Beispiel:

.. skip

::

    sage: 2^10
    1024
    sage: sage: sage: 2^10
    1024
    sage: >>> 2^10
    1024

Befehle zur Zeitmessung
=======================

Wenn Sie den ``%time`` Befehl vor eine Eingabe schreiben wird die
Zeit, die der Aufruf benötigt, ausgegeben nachdem er gelaufen ist.
Zum Beispiel können wir die Laufzeit einer bestimmten Potenzierung auf
verschiedene Arten vergleichen. Die unten genannte Laufzeit wird unter
Umständen weit von der Laufzeit auf Ihrem Computer oder sogar zwischen
verschiedenen SageMath Versionen abweichen. Zuerst natives Python:

.. skip

::

    sage: %time a = int(1938)^int(99484)
    CPU times: user 0.66 s, sys: 0.00 s, total: 0.66 s
    Wall time: 0.66

Das bedeutet insgesamt 0,66 Sekunden wurden benötigt und die
vergangene "Wall time", also die vergangene Echtzeit (auf Ihrer
Wanduhr), betrug auch 0,66 Sekunden. Wenn auf Ihrem Computer viele
andere Programme gleichzeitig laufen kann die "Wall time"
wesentlich größer als die CPU Zeit sein.

Als nächstes messen wir die Laufzeit der Potenzierung unter Verwendung
des nativen Sage Ganzzahl-Typs, der (in Cython implementiert ist und)
die GMP Bibliothek nutzt:

.. skip

::

    sage: %time a = 1938^99484
    CPU times: user 0.04 s, sys: 0.00 s, total: 0.04 s
    Wall time: 0.04

Unter Verwendung der PARI C-Bibliothek:

.. skip

::

    sage: %time a = pari(1938)^pari(99484)
    CPU times: user 0.05 s, sys: 0.00 s, total: 0.05 s
    Wall time: 0.05

GMP ist also ein bisschen besser (wie erwartet, da die für Sage
verwendete PARI Version GMP für Ganzzahlarithmetik nutzt).
Sie können ebenso Befehlsblöcke messen, indem Sie ``cputime`` wie
unten verwenden:

::

    sage: t = cputime()
    sage: a = int(1938)^int(99484)
    sage: b = 1938^99484
    sage: c = pari(1938)^pari(99484)
    sage: cputime(t)                       # random output
    0.64

.. skip

::

    sage: cputime?
    ...
        Return the time in CPU second since Sage started, or with optional
        argument t, return the time since time t.
        INPUT:
            t -- (optional) float, time in CPU seconds
        OUTPUT:
            float -- time in CPU seconds

Der ``walltime`` Befehl entspricht ``cputime``, nur misst dieser die Echtzeit.

Wir können die oben genannte Potenz auch in einigen der Computer
Algebra Systeme, die Sage mitbringt berechnen. In jedem Fall führen wir
einen trivialen Befehl aus, um den entsprechenden Server dieses
Programms zu starten. Sollte es erhebliche Unterschiede zwischen
Echtzeit und CPU-Zeit geben, deutet dies auf ein Leistungsproblem hin,
dem man nachgehen sollte.

.. skip

::

    sage: time 1938^99484;
    CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
    Wall time: 0.01
    sage: gp(0)
    0
    sage: time g = gp('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: maxima(0)
    0
    sage: time g = maxima('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.30
    sage: kash(0)
    0
    sage: time g = kash('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: mathematica(0)
            0
    sage: time g = mathematica('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.03
    sage: maple(0)
    0
    sage: time g = maple('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.11
    sage: gap(0)
    0
    sage: time g = gap.eval('1938^99484;;')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 1.02

Achten Sie darauf, dass GAP und Maxima am langsamsten in diesem Test
sind (er lief auf dem Computer ``sage.math.washington.edu``). Aufgrund
des Pexpect-Schnittstellen-Overheads ist es aber vielleicht unfair
diese mit Sage zu vergleichen, welches am schnellsten war.

Fehlerbehandlung
================

Wenn irgendetwas schief geht, werden Sie normalerweise eine
Python-Fehlermeldung sehen. Python macht sogar einen Vorschlag, was den
Fehler ausgelöst hat. Oft sehen Sie den Namen der Fehlermeldung,
z.B. :class:`NameError` oder :class:`ValueError` (vgl. Python Library Reference
[PyLR]_ für eine komplette Liste der Fehlermeldungen). Zum Beispiel:

.. skip

::

    sage: 3_2
    ------------------------------------------------------------
       File "<console>", line 1
         ZZ(3)_2
               ^
    SyntaxError: invalid ...

    sage: EllipticCurve([0,infinity])
    ------------------------------------------------------------
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce Infinity (<class 'sage...Infinity'>) to Rational

Der interaktive Debugger ist manchmal hilfreich um zu verstehen was
schiefgelaufen ist. Sie können ihn ein- oder ausschalten indem Sie
``%pdb`` eingeben (standardmäßig ist er ausgeschaltet). Die
Eingabeaufforderung ``ipdb>`` erscheint wenn eine Fehlermeldung
geworfen wird und der Debugger eingeschaltet ist. Im Debugger können
Sie den Status jeder lokalen Variable ausgeben oder im Ausführungstack
hoch- und runterspringen.
Zum Beispiel:

.. skip

::

    sage: %pdb
    Automatic pdb calling has been turned ON
    sage: EllipticCurve([1,infinity])
    ---------------------------------------------------------------------------
    <class 'exceptions.TypeError'>             Traceback (most recent call last)
    ...

    ipdb>

Tippen Sie ``?`` in der ``ipdb>``-Eingabeaufforderung  um eine Liste
der Befehle des Debuggers zu erhalten.

::

    ipdb> ?

    Documented commands (type help <topic>):
    ========================================
    EOF    break  commands   debug    h       l     pdef   quit    tbreak
    a      bt     condition  disable  help    list  pdoc   r       u
    alias  c      cont       down     ignore  n     pinfo  return  unalias
    args   cl     continue   enable   j       next  pp     s       up
    b      clear  d          exit     jump    p     q      step    w
    whatis where

    Miscellaneous help topics:
    ==========================
    exec  pdb

    Undocumented commands:
    ======================
    retval  rv

Drücken Sie Strg-D oder geben Sie ``quit`` ein um zu Sage zurückzukehren.

.. _section-tabcompletion:

Rückwärtssuche und Tab-Vervollständigung
========================================

Definieren Sie zuerst einen dreidimensionalen Vektorraum
:math:`V=\QQ^3` wie folgt:

::

    sage: V = VectorSpace(QQ,3)
    sage: V
    Vector space of dimension 3 over Rational Field

Sie können auch die folgende verkürzte Schreibweise verwenden:

::

    sage: V = QQ^3

Schreiben Sie den Anfang eines Befehls und drücken Sie dann ``Strg-p``
(oder drücken Sie einfach die Pfeil-nach-oben-Taste) um zur vorher
eingegebenen Zeile zu gelangen, die ebenfalls so beginnt. Das
funktioniert auch nach einem kompletten Sage-Neustart noch. Sie können
den Verlauf auch mit ``Strg-r`` rückwärts durchsuchen.  Diese
Funktionalität wird vom ``readline``-Paket bereitgestellt, welches in
nahezu jeder Linux-Distribution verfügbar ist.

Es ist sehr einfach alle Unterfunktionen für :math:`V` mittels
Tab-Vervollständigung  aufzulisten, indem Sie erst ``V.`` eingeben,
und dann die ``[Tabulator Taste]`` drücken:

.. skip

::

    sage: V.[tab key]
    V._VectorSpace_generic__base_field
    ...
    V.ambient_space
    V.base_field
    V.base_ring
    V.basis
    V.coordinates
    ...
    V.zero_vector

Wenn Sie die ersten paar Buchstaben einer Funktion tippen und dann die
``[Tabulator Taste]`` drücken, bekommen Sie nur die Funktionen, die so
beginnen angezeigt.

.. skip

::

    sage: V.i[tab key]
    V.is_ambient  V.is_dense    V.is_full     V.is_sparse

Wenn sie wissen wollen, was eine bestimmte Funktion tut, z.B. die
"coordinates"-Funktion, so geben Sie ``V.coordinates?`` ein um die
Hilfe, und ``V.coordinates??`` um den Quelltext der Funktion zu
sehen.


Integriertes Hilfesystem
========================

Sage hat ein integriertes Hilfesystem. Hängen Sie an einen beliebigen
Funktionsnamen ein ``?`` an, um die Dokumentation dazu aufzurufen.

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <bound method FreeModule_ambient_field.coordinates of Vector
    space of dimension 3 over Rational Field>
    Namespace:      Interactive
    File:           /home/was/s/local/lib/python2.4/site-packages/sage/modules/f
    ree_module.py
    Definition:     V.coordinates(self, v)
    Docstring:
        Write v in terms of the basis for self.

        Returns a list c such that if B is the basis for self, then

                sum c_i B_i = v.

        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: M = FreeModule(IntegerRing(), 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]

Wie Sie sehen, beinhaltet die Ausgabe den Typ des Objekts, den
Dateinamen in welcher die Funktion definiert ist und eine Beschreibung der Funktionalität
mit Beispielen, die Sie direkt in Ihre aktuelle Sitzung einfügen können.
Fast alle dieser Beispiele werden regelmäßig automatisch getestet um sicherzustellen, dass sie
genau wie beschrieben funktionieren.

Eine andere Funktionalität, die sehr eng in Verbindung mit Open-Source-Gedanken steht ist,
dass Sie sich zu jeder Funktion den Quelltext anzeigen lassen
können. Sei ``f`` eine Sage oder Python Funktion, dann können Sie mit
``f??`` den Quellcode, der ``f`` definiert anzeigen. Zum Beispiel:

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates??
    Type:           instancemethod
    ...
    Source:
    def coordinates(self, v):
            """
            Write $v$ in terms of the basis for self.
            ...
            """
            return self.coordinate_vector(v).list()

Das zeigt uns, dass die ``coordinates``-Funktion nichts anderes tut,
als ``coordinates_vector``-Funktion aufruft und das Ergebnis in eine
Liste umwandelt. Aber was tut die ``coordinates``-Funktion?

.. skip

::

    sage: V = QQ^3
    sage: V.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            ...
            return self.ambient_vector_space()(v)

Die ``coordinate_vector``-Funktion steckt ihre Eingabe in den
umgebenden Raum, was zur Folge hat, dass der Koeffizientenvektor von
:math:`v` zur Basis des Vektorraums :math:`V` ausgerechnet wird.
Der Raum :math:`V` ist schon der umgebende, nämlich gerade
:math:`\QQ^3`. Es gibt auch eine ``coordinate_vector``-Funktion für
Unterräume, und sie funktioniert anders.
Wir definieren einen Unterraum und schauen uns das an:

.. skip

::

    sage: V = QQ^3; W = V.span_of_basis([V.0, V.1])
    sage: W.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            """
             ...
            """
            # First find the coordinates of v wrt echelon basis.
            w = self.echelon_coordinate_vector(v)
            # Next use transformation matrix from echelon basis to
            # user basis.
            T = self.echelon_to_user_matrix()
            return T.linear_combination_of_rows(w)

(Wenn Sie der Meinung sind, dass diese Implementation ineffizient ist,
helfen Sie uns bitte unsere Lineare Algebra zu optimieren.)

Sie können auch ``help(command_name)`` oder ``help(class)`` eingeben
um eine manpage-artige Hilfe zu bekommen.


.. skip

::

    sage: help(VectorSpace)
    Help on function VectorSpace in module sage.modules.free_module:

    VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args)
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x
         over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
    --More--

Wenn Sie ``q`` drücken um das Hilfesystem zu verlassen, kommen Sie genau
dahin zurück, wo Sie Ihre Sitzung verlassen haben. Die ``help`` Anzeige
bleibt nicht in Ihrer Sitzung zurück im Gegensatz zu ``funktion?``.
Es ist besonders hilfreich ``help(modul_name)`` zu nutzen. Zum Beispiel sind
Vektorräume in ``sage.modules.free_module`` definiert. Geben Sie also
``help(sage.modules.free_module)`` ein, um die Dokumentation des
ganzen Moduls zu sehen. Wenn Sie sich Die Dokumentation mit ``help``
ansehen, können Sie mit ``/`` vorwärts und mit ``?`` rückwärts suchen.

Speichern und Laden von individuellen Objekten
==============================================

Angenommen Sie berechnen eine Matrix oder schlimmer, einen
komplizierten Modulsymbolraum, und Sie wollen ihn für später
speichern. Was können Sie tun? Es gibt mehrere Möglichkeiten für
Computer Algebra Systeme solche individuellen Objekte zu speichern.


#. **speichern Ihres Spiels:** Unterstützt nur das Speichern und Laden kompletter
   Sitzungen (z.B. GAP, Magma).

#. **Einheitliche Ein-/Ausgabe:** Bringen Sie jedes Objekt in eine Form, die
         Sie wieder einlesen können in (GP/PARI).

#. **Eval**: Machen Sie beliebigen Code auswertbar im Interpreter (z.B. Sigular, PARI).


Da Sage Python nutzt, braucht es einen anderen Ansatz, nämlich dass
jedes Objekt serialisiert werden kann. Das heißt es in eine Zeichenkette
umzuwandeln, die man wieder einlesen kann. Das ist im Prinzip ähnlich zum
einheitlichen Ein-/Ausgabe Ansatz von PARI, abgesehen von der zu komplizierten
Darstellung auf dem Bildschirm. Außerdem ist das Laden und Speichern (meistens)
vollautomatisch und benötigt nicht einmal speziellen Programmieraufwand; es ist
einfach ein Merkmal, das von Grund auf in Python war.

Fast alle Objekte x in Sage können in komprimierter Form gespeichert werden
via ``save(x, Dateiname)`` (oder in vielen Fällen ``x.save(Dateiname)``).
Um das Objekt wieder zu laden, nutzen Sie ``load(Dateiname)``.

.. skip

::

    sage: A = MatrixSpace(QQ,3)(range(9))^2
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]
    sage: save(A, 'A')

Sie sollten Sage nun schließen und neu starten. Dann können Sie ``A`` wieder laden:

.. skip

::

    sage: A = load('A')
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]

Sie können das selbe mit komplizierteren Objekten, wie etwa elliptischen
Kurven machen. Alle Daten über das Objekt sind zwischengespeichert und
werden mit dem Objekt gespeichert. Zum Beispiel:

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: v = E.anlist(100000)              # dauert etwas länger
    sage: save(E, 'E')
    sage: quit

Die gespeicherte Version von ``E`` braucht 153 Kilobyte, da die ersten
100000 :math:`a_n` mitgespeichert werden.

.. skip

::

    ~/tmp$ ls -l E.sobj
    -rw-r--r--  1 was was 153500 2006-01-28 19:23 E.sobj
    ~/tmp$ sage [...]
    sage: E = load('E')
    sage: v = E.anlist(100000)              # sofort!

(In Python wird das Laden und Speichern mittels des ``cPickle``
Moduls umgesetzt. Genauer: Ein Sage Objekt ``x`` kann mit
``cPickle.dumps(x, 2)`` gespeichert werden.  Beachten Sie die ``2``!)

Sage kann allerdings keine individuellen Objekte anderer Computer Algebra Systeme
wie GAP, Singular, Maxima, usw. laden und speichern. Sie sind mit "invalid" gekennzeichnet nach dem Laden.
In GAP werden viele Objekte in einer Form dargestellt, die man wiederherstellen kann,
viele andere allerdings nicht. Deshalb ist das Wiederherstellen aus ihren Druckdarstellungen
nicht erlaubt.

.. skip

::

    sage: a = gap(2)
    sage: a.save('a')
    sage: load('a')
    Traceback (most recent call last):
    ...
    ValueError: The session in which this object was defined is no longer
    running.

GP/PARI Objekte können hingegen gespeichert und geladen werden, da
ihre Druckdarstellung ausreicht um sie wiederherzustellen.

.. skip

::

    sage: a = gp(2)
    sage: a.save('a')
    sage: load('a')
    2

Gespeicherte Objekte können auch auf Computern mit anderen Architekturen
oder Betriebssystemen wieder geladen werden. Zum Beispiel können Sie
eine riesige Matrix auf einem 32 Bit Mac OS X speichern und später auf
einem 64 Bit Linux System laden, dort die Stufenform herstellen und dann
wieder zurückladen. Außerdem können Sie in den meisten Fällen auch Objekte
laden, die mit anderen SageMath Versionen gespeichert wurden, solange der Quelltext
des Objekts nicht zu verschieden ist. Alle Attribute eines Objekts werden zusammen
mit seiner Klasse (aber nicht dem Quellcode) gespeichert. Sollte diese Klasse
in einer neueren SageMath Version nicht mehr existieren, kann das Objekt in dieser
neueren SageMath Version nicht mehr geladen werden. Aber Sie könnten es in der alten
SageMath Version laden, die Objekt Dictionaries mit ``x.__dict__`` laden und das Objekt
zusammen mit diesem in der neuen SageMath Version laden.

Als Text speichern
------------------

Sie können die ASCII Text Darstellung eines Objekts in eine Klartextdatei
schreiben, indem Sie die Datei einfach mit Schreibzugriff öffnen und die
Textdarstellung des Objekts hineinkopieren. (Sie können auch viele andere
Objekte auf diese Art speichern.) Wenn Sie alle Objekte hineinkopiert haben,
schließen Sie die Datei einfach.

.. skip

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = (x+y)^7
    sage: o = open('file.txt','w')
    sage: o.write(str(f))
    sage: o.close()

.. _section-save:

Speichern und Laden kompletter Sitzungen
========================================

Sage hat eine sehr flexible Unterstützung für das Speichern und Laden
kompletter Sitzungen.

Der Befehl ``save_session(sitzungsname)`` speichert alle Variablen,
die Sie während dieser Sitzung definiert haben als ein Dictionary
``sessionname``. (Im seltenen Fall, dass eine Variable nicht gespeichert
werden kann, fehlt sie anschließend einfach im Dictionary.)
Die erzeugte Datei ist eine ``.sobj``-Datei und kann genau wie jedes andere
Objekt geladen werden. Wenn Sie Objekte aus einer Sitzung laden, werden Sie
diese in einem Dictionary finden. Dessen Schlüssel sind die Variablen und
dessen Werte sind die Objekte.

Sie können den ``load_session(sitzungsname)`` Befehl nutzen um die Variablen
aus ``sitzungsname`` in die aktuelle Sitzung zu laden. Beachten Sie, dass
dieses Vorgehen nicht die Variablen der aktuellen Sitzung löscht, vielmehr
werden beide Sitzungen vereinigt.

Starten wir also zunächst Sage und definieren einige Variablen.

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: M = ModularSymbols(37)
    sage: a = 389
    sage: t = M.T(2003).matrix(); t.charpoly().factor()
     _4 = (x - 2004) * (x - 12)^2 * (x + 54)^2

Als nächstes speichern wir unsere Sitzung, was jede der Variablen
in eine Datei speichert. Dann sehen wir uns die Datei, die etwa
3 Kilobyte groß ist an.

.. skip

::

    sage: save_session('misc')
    Saving a
    Saving M
    Saving t
    Saving E
    sage: quit
    was@form:~/tmp$ ls -l misc.sobj
    -rw-r--r--  1 was was 2979 2006-01-28 19:47 misc.sobj

Zuletzt starten wir Sage neu, definieren uns eine extra Variable, und laden
unsere gespeicherte Sitzung.

.. skip

::

    sage: b = 19
    sage: load_session('misc')
    Loading a
    Loading M
    Loading E
    Loading t

Jede der gespeicherten Variablen ist wieder verfügbar und die
Variable ``b`` wurde nicht überschrieben.

.. skip

::

    sage: M
    Full Modular Symbols space for Gamma_0(37) of weight 2 with sign 0
    and dimension 5 over Rational Field
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational
    Field
    sage: b
    19
    sage: a
    389

