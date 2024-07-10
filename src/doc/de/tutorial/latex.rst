****************************
Sage, LaTeX und ihre Freunde
****************************

Sage und der TeX-Dialekt LaTeX haben eine sehr synergetische
Beziehung. Dieses Kapitel hat das Ziel die Vielfalt an Interaktionen,
von den einfachsten bis hin zu den ungewöhnlichen und fast schon
magischen, vorzustellen. (Sie sollten also nicht gleich das ganze
Kapitel im ersten Durchgang durch das Tutorial lesen.)

Überblick
=========

Es ist wahrscheinlich am einfachsten die verschiedenen
Einsatzmöglichkeiten von LaTeX zu verstehen, wenn man sich die drei
grundsätzlichen Methoden in Sage ansieht.

    #. Jedes Objekt in Sage muss eine LaTeX Darstellung haben.
       Sie können diese Darstellung erreichen, indem Sie im Notebook
       oder der Kommandozeile ``latex(foo)`` ausführen, wobei ``foo``
       ein Objekt in  Sage ist. Die Ausgabe ist eine Zeichenkette, die
       eine recht genaue Darstellung  im mathematischen Modus von TeX
       bietet (z.B. zwischen jeweils zwei Dollarzeichen). Einige
       Beispiele hierfür folgen unten.

       So kann Sage effektiv genutzt werden um Teile eines
       LaTeX-Dokuments zu erstellen:
       Erstellen oder berechnen Sie ein Objekt in Sage, drucken Sie es
       mit dem ``latex()``-Befehl  aus und fügen Sie es in Ihr Dokument ein.

    #. Die Notebook Schnittstelle ist konfiguriert
       `MathJax <http://www.mathjax.org>`_
       zu nutzen um mathematische Ausdrücke im Browser darzustellen.
       MathJax ist eine Kollektion aus JavaScript-Routinen und
       zugehörigen Schriftarten. Es ist also nichts zusätzlich
       einzustellen um mathematische Ausdrücke in Ihrem Browser
       anzuzeigen, wenn Sie das Sage-Notebook nutzen.

       MathJax wurde entwickelt um einen großen, aber nicht vollständigen
       Teil von TeX darstellen zu können. Es gibt keine Unterstützung
       für Dinge, wie komplizierte Tabellen, Kapiteleinteilung oder
       Dokument Management, da es für genaues Darstellen von TeX
       Ausdrücken konzipiert wurde. Die nahtlose Darstellung von
       mathematischen Ausdrücken im Sage Notebook wird durch
       Konvertierung der ``latex()``-Darstellung in MathJax
       gewährleistet.

       Da MathJax seine eigenen skalierbaren Schriftarten nutzt, ist es
       anderen Methoden überlegen, die auf Konvertierung in kleine
       Bilder beruhen.

    #. Sollte in der Sage Kommandozeile oder im Notebook mehr
       LaTeX-Code vorkommen als MathJax verarbeiten kann, kann eine
       systemweite Installation von LaTeX aushelfen. Sage beinhaltet
       fast alles, das Sie brauchen um Sage weiter zu entwickeln und
       zu nutzen. Eine Ausnahme hierzu ist TeX selbst. In diesen
       Situationen müssen also TeX und verschiedene Konverter
       installiert sein, um alle Möglichkeiten nutzen zu können.

Hier führen wir einige grundlegenden Funktionen von ``latex()`` vor. ::

    sage: var('z')
    z
    sage: latex(z^12)
    z^{12}
    sage: latex(integrate(z^4, z))
    \frac{1}{5} \, z^{5}
    sage: latex('a string')
    \text{\texttt{a{ }string}}
    sage: latex(QQ)
    \Bold{Q}
    sage: latex(matrix(QQ, 2, 3, [[2,4,6],[-1,-1,-1]]))
    \left(\begin{array}{rrr}
    2 & 4 & 6 \\
    -1 & -1 & -1
    \end{array}\right)

Grundlegende MathJax Funktionen gibt es im Notebook weitgehend automatisch,
aber wir können es teilweise mit Hilfe der ``MathJax`` Klasse demonstrieren.
Die ``eval`` Funktion dieser Klasse konvertiert ein Sage-Objekt in
seine LaTeX-Darstellung und dann in HTML mit der CSS ``math`` Klasse,
die dann MathJax verwendet. ::

    sage: from sage.misc.html import MathJax
    sage: mj = MathJax()
    sage: var('z')
    z
    sage: mj(z^12)
    <html>\[z^{12}\]</html>
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>
    sage: mj(ZZ['x'])
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}[x]\]</html>
    sage: mj(integrate(z^4, z))
    <html>\[\frac{1}{5} \, z^{5}\]</html>

Grundlegende Nutzung
====================

Wie schon im Überblick angekündigt, ist der einfachste Weg Sage's
LaTeX-Unterstützung zu nutzen die ``latex()`` Funktion um eine
legitime LaTeX-Darstellung eines mathematischen Objekts zu erhalten.
Diese Zeichenketten können dann in unabhängigen LaTeX-Dokumenten
genutzt werden. Das funktioniert im Notebook genauso wie in der
Sage-Kommandozeile.

Das andere Extrem ist der ``view()``-Befehl. In der Sage-Kommandozeile
wird der Befehl ``view()`` die LaTeX-Darstellung von ``foo`` in ein
einfaches  LaTeX Dokument packen, und dann dieses mit der systemweiten
TeX-Installation aufrufen. Zuletzt wird das passende Programm zum
Anzeigen der Ausgabe von TeX aufgerufen. Welche Version von TeX
genutzt wird, und damit auch wie die Ausgabe aussieht und welches
Anzeigeprogramm aufgerufen wird, kann angepasst werden (siehe
:ref:`sec-custom-processing`).

Im Notebook schafft der ``view(foo)`` Befehl die nötige Kombination
von HTML und CSS sodass MathJax die LaTeX Darstellung im Arbeitsblatt
anzeigt. Für den Anwender erstellt er einfach eine schön formatierte
Ausgabe, die sich von der normalen ASCII Ausgabe aus Sage
unterscheidet. Nicht jedes mathematische Objekt in Sage hat eine
LaTeX-Darstellung, die die eingeschränkten Möglichkeiten von MathJax
unterstützt. In diesen Fällen kann die MathJax Darstellung umgangen
werden, und stattdessen die systemweite TeX-Installation aufgerufen
werden. Dessen Ausgabe kann dann als Bild im Arbeitsblatt angezeigt
werden. Die Einstellungen und Auswirkungen dieses Prozesses wird im
Kapitel :ref:`sec-custom-generation` dargestellt.

Der interne ``pretty_print()`` Befehl zeigt die Konvertierung von Sage
Objekten in HTML Code der MathJax nutzt im Notebook.  ::

    sage: pretty_print(x^12)
    x^12
    sage: pretty_print(integrate(sin(x), x))
    -cos(x)

Das Notebook hat zwei weitere Möglichkeiten TeX zu nutzen. Die erste
ist der "Typeset"-Knopf über der ersten Zelle eines Arbeitsblatts,
rechts von den vier Drop-Down-Boxen. Ist er ausgewählt werden die
Ausgaben aller folgenden Berechnungen von MathJax
interpretiert. Beachten Sie, dass dieser Befehl nicht rückwirkend ist
-- alle vorher berechneten Zellen werden nicht neu berechnet. Im
Grunde ist der "Typeset"-Knopf gleichzusetzen mit dem Aufruf des
``view()``-Befehls in allen Zellen.

Die zweite Möglichkeit im Notebook ist das Eingeben von TeX
Kommentaren in einem Arbeitsblatt. Wenn der Cursor zwischen zwei
Zellen steht, und der erscheinende blaue Balken mit gedrückter Shift
Taste geklickt wird, wird ein  kleiner Texteditor TinyMCE
geöffnet. Dieser erlaubt die Eingabe von HTML und CSS formatiertem
Text mit einem WYSIWYG-Editor. Es ist also möglich den so formatierten
Text als Kommentar in einem  Arbeitsblatt unterzubringen. Text den Sie
hier zwischen ``$...$`` oder ``$$...$$`` eingeben wird ebenfalls von
MathJax in einer "inline" bzw. "display math" Umgebung gesetzt.

.. _sec-custom-generation:


Anpassen der LaTeX-Generierung
==============================

Es gibt verschiedene Arten den vom ``latex()``-Befehl generierten
LaTeX-Code anzupassen. Im Notebook und der Sage Kommandozeile gibt es
ein vordefiniertes Objekt Namens ``latex``, das verschiedene Methoden
hat, die Sie sich auflisten lassen können indem Sie ``latex.``
eingeben und die Tab Taste drücken (beachten Sie den Punkt).

Ein gutes Beispiel ist die ``latex.matrix_delimiters`` Methode. Es
kann benutzt werden um die Darstellung der Matrizen zu beeinflussen --
runde Klammern, eckige Klammern, geschwungene Klammern oder senkrechte
Striche. Sie müssen sich nicht für eine Darstellung entscheiden, Sie
können verschiedene miteinander kombinieren, wie Sie es
wünschen. Beachten Sie dass die in LaTeX benötigten  Backslashes einen
zusätzlichen Slash benötigen damit sie in Python korrekt erkannt
werden. ::

    sage: A = matrix(ZZ, 2, 2, range(4))
    sage: latex(A)
    \left(\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right)
    sage: latex.matrix_delimiters(left='[', right=']')
    sage: latex(A)
    \left[\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right]
    sage: latex.matrix_delimiters(left='\\{', right='\\}')
    sage: latex(A)
    \left\{\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right\}

Die ``latex.vector_delimiters`` Methode funktioniert ähnlich.

Die Darstellung von Ringen und Körpern (ganze, rationale, reelle
Zahlen, etc.) kann mit der ``latex.blackboard_bold`` Methode verändert
werden. Diese Mengen werden in standardmäßig in fett gedruckt,
alternativ können sie auch mit Doppelstrichen geschrieben
werden. Hierfür wird das  ``\Bold{}``-Makro genutzt, das in Sage
integriert ist. ::

    sage: latex(QQ)
    \Bold{Q}
    sage: from sage.misc.html import MathJax
    sage: mj=MathJax()
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>
    sage: latex.blackboard_bold(True)
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbb{#1}}\Bold{Q}\]</html>
    sage: latex.blackboard_bold(False)

Dank der Erweiterbarkeit von TeX können Sie selbst Makros und Pakete
einbinden. Individuelle Makros können hinzugefügt werden, die dann von
MathJax als TeX-Schnipsel  interpretiert werden. ::

    sage: latex.extra_macros()
    ''
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: latex.extra_macros()
    '\\newcommand{\\foo}{bar}'
    sage: var('x y')
    (x, y)
    sage: latex(x+y)
    x + y
    sage: from sage.misc.html import MathJax
    sage: mj=MathJax()
    sage: mj(x+y)
    <html>\[\newcommand{\foo}{bar}x + y\]</html>

Zusätzliche Makros, die so hinzugefügt wurden, werden auch vom
systemweiten TeX genutzt, wenn MathJax an seine Grenzen gestoßen ist.
Der Befehl ``latex_extra_preamble`` kann genutzt werden um eine
Präambel eines kompletten LaTeX Dokuments zu erzeugen, das folgende
Beispiel zeigt wie. Beachten Sie wiederrum die doppelten Backslashes
in den Python Zeichenketten. ::


    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: from sage.misc.latex import latex_extra_preamble
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    \newcommand{\foo}{bar}

Für größere oder kompliziertere LaTeX-Ausdrücke können mit
``latex.add_to_preamble`` Pakete (oder ähnliches) zur LaTeX-Präambel
hinzugefügt werden. Der zweite Befehl
``latex.add_package_to_preamble_if_available``  prüft hingegen erst ob
das Paket vorhanden ist, bevor es eingebunden wird.

Hier fügen wir das geometry-Paket zur Präambel hinzu, um die
Seitenränder einzustellen. Achten Sie wieder auf die doppelten
Backslashes in Python. ::


    sage: from sage.misc.latex import latex_extra_preamble
    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: latex.add_to_preamble('\\usepackage{geometry}')
    sage: latex.add_to_preamble('\\geometry{letterpaper,total={8in,10in}}')
    sage: latex.extra_preamble()
    '\\usepackage{geometry}\\geometry{letterpaper,total={8in,10in}}'
    sage: print(latex_extra_preamble())
    \usepackage{geometry}\geometry{letterpaper,total={8in,10in}}
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}

Ein bestimmtes Paket, dessen Existenz nicht sicher ist, wird wie folgt
eingebunden. ::

    sage: latex.extra_preamble('')
    sage: latex.extra_preamble()
    ''
    sage: latex.add_to_preamble('\\usepackage{foo-bar-unchecked}')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'
    sage: latex.add_package_to_preamble_if_available('foo-bar-checked')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'

.. _sec-custom-processing:

Anpassen der LaTeX-Verarbeitung
===============================

Es ist möglich zu entscheiden welche Variante von TeX für einen
systemweiten Aufruf genutzt werden soll, und somit auch wie die
Ausgabe aussehen soll. Ebenso ist es möglich zu beeinflussen, ob das
Notebook MathJax oder die systemweite LaTeX Installation nutzt.

Der Befehl ``latex.engine()`` entscheidet, ob die systemweiten
Anwendungen ``latex``, ``pdflatex`` oder ``xelatex`` genutzt werden
für kompliziertere LaTeX-Ausdrücke. Wenn ``view()`` in der Sage
Kommandozeile aufgerufen wird, und ``latex`` als Prozessor eingestellt
ist, wird eine .dvi Datei erzeugt, die dann mit einem dvi
Anzeigeprogramm (wie xdvi) angezeigt wird. Im Gegensatz hierzu wird
bei Aufruf von ``view()`` mit dem Prozessor ``pdflatex`` eine .PDF
Datei erzeugt, die mit dem Standard-PDF-Programm angezeigt
wird. (acrobat, okular, evince, etc.).

Im Notebook kann es nötig sein, dem System die Entscheidung
abzunehmen, ob MathJax für einige TeX-Schnipsel, oder das systemweite
LaTeX für kompliziertere Ausdrücke genutzt werden soll. Es gibt eine
Liste von Befehlen, die wenn einer von ihnen in einem Stück LaTeX-Code
erkannt wird, die Ausgabe von LaTeX (oder welcher Prozessor auch immer
durch ``latex.engine()`` gesetzt ist) statt von MathJax erstellen
lässt. Diese Liste wird verwaltet durch die Befehle
``latex.add_to_mathjax_avoid_list`` und
``latex.mathjax_avoid_list``. ::

    sage: # not tested
    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []
    sage: latex.mathjax_avoid_list(['foo', 'bar'])
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar']
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar', 'tikzpicture']
    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []

Nehmen wir an ein LaTeX-Ausdruck wurde im Notebook durch ``view()``
oder während aktiviertem "Typeset" Knopf erzeugt. Und dann wird
festgestellt, dass er die externe LaTeX-Installation benötigt, weil
er in der ``mathjax_avoid_list`` steht. Der Ausdruck wird nun vom
ausgewählten (durch ``latex.engine()``) Prozessor erzeugt, und statt
der Anzeige in einem externen Programm (was in der Kommandozeile
passieren  würde) wird Sage versuchen das Ergebnis in einem einzigen,
leicht beschnittenen Bild in der Ausgabezelle darzustellen.

Wie diese Umwandlung abläuft hängt von einigen Faktoren ab,
hauptsächlich vom verwendeten LaTeX-Prozessor und davon welche
Konvertierungswerkzeuge auf dem System vorhanden sind. Vier nützliche
Konverter, die alle Eventualitäten abdecken sind ``dvips``,
``ps2pdf``, ``dvipng`` und aus dem ``ImageMagick`` Paket,
``convert``. Das Ziel ist die Erzeugung einer .png Datei, die später
wieder im Arbeitsblatt eingebunden werden kann. Wenn ein
LaTeX-Ausdruck erfolgreich von ``latex`` in eine .dvi Datei verwandelt
wird, dann sollte dvipng die Umwandlung vornehmen. Wenn der LaTeX
Ausdruck und der gewählte LaTeX-Prozessor eine .dvi Datei mit
Erweiterungen erstellt, die dvipng nicht unterstützt, so wird dvips
eine PostScript-Datei erzeugen. So eine PostScript-Datei, oder eine
.pdf Datei aus dem Prozessor ``pdflatex``, wird dann von ``convert``
in eine .png Datei gewandelt. Das Vorhandensein von zweier solcher
Konverter kann mit Hilfe der ``have_dvipng()`` und ``have_convert()``
Routinen überprüft werden.

Diese Umwandlungen werden automatisch ausgeführt, wenn Sie die nötigen
Konverter installiert haben; falls nicht wird Ihnen eine Fehlermeldung
angezeigt, die Ihnen sagt was fehlt und wo Sie es herunterladen können.

Für ein konkretes Beispiel wie komplizierte LaTeX-Ausdrücke
verarbeitet werden können, sehen Sie sich das Beispiel des
``tkz-graph`` Pakets zum Erstellen von hochwertigen kombinatorischen
Graphen im nächsten Abschnitt (:ref:`sec-tkz-graph`) an. Für weitere
Beispiele gibt es einige vorgepackte Testfälle. Um diese zu nutzen,
müssen Sie das ``sage.misc.latex.latex_examples`` Objekt
importieren. Dieses ist eine Instanz der
``sage.misc.latex.LatexExamples`` Klasse, wie unten beschrieben. Diese
Klasse enthält momentan Beispiele von kommutativen Diagrammen,
kombinatorischen Graphen, Knotentheorie und Beispiele für Graphen mit
pstricks. Es werden damit die folgenden Pakete getestet: xy,
tkz-graph, xypic, pstricks.  Nach dem Import können Sie mittels
Tab-Vervollständigung von ``latex_examples`` die vorgepackten
Beispiele sehen. Bei Aufruf vom jedem Beispiel erhalten Sie eine
Erklärung was nötig ist, damit das Beispiel korrekt dargestellt
wird. Um die Darstellung tatsächlich zu sehen müssen Sie ``view()``
benutzen (sofern die Präambel, der LaTeX-Prozessor, etc richtig
eingestellt sind).
::

    sage: from sage.misc.latex import latex_examples
    sage: latex_examples.diagram()
    LaTeX example for testing display of a commutative diagram produced
    by xypic.
    <BLANKLINE>
    To use, try to view this object -- it will not work.  Now try
    'latex.add_to_preamble("\\usepackage[matrix,arrow,curve,cmtip]{xy}")',
    and try viewing again. You should get a picture (a part of the diagram arising
    from a filtered chain complex).

.. _sec-tkz-graph:

Ein Beispiel: Kombinatorische Graphen mit tkz-graph
===================================================

Hochwertige Darstellungen von kombinatorischen Graphen (fortan nur
noch "Graphen") sind mit Hilfe des ``tkz-graph`` Pakets möglich.
Dieses Paket wurde ausbauend auf das ``tikz`` front-end der ``pgf``
Bibliothek entwickelt. Es müssen also all diese Komponenten Teil der
systemweiten TeX-Installation sein, und es ist möglich, dass sie nicht
in ihrer neusten Version in der TeX-Implementation vorliegen. Es ist
also unter Umständen nötig oder ratsam diese Teile separat in Ihrem
persönlichen texmf Baum zu installieren. Das Erstellen, Anpassen und
Warten einer systemweiten oder persönlichen TeX-Installation würde
allerdings den Rahmen dieses Dokuments sprengen. Es sollte aber
einfach sein Anleitungen hierzu zu finden. Die nötigen Dateien sind
unter :ref:`sec-system-wide-tex` aufgeführt.

Um also zu beginnen, müssen wir sicher sein, dass die relevanten
Pakete eingefügt werden, indem wir sie in die Präambel des
LaTeX-Dokuments hinzufügen. Die Bilder der Graphen werden nicht
korrekt formatiert sein, wenn eine .dvi Datei als Zwischenergebnis
erzeugt wird. Es ist also ratsam, den LaTeX-Prozessor auf
``pdflatex`` zu stellen. Nun sollte ein Befehl wie
``view(graphs.CompleteGraph(4))`` in der Sage-Kommandozeile
erfolgreich eine .pdf Datei mit einem  Bild vom kompletten `K_4`
Graphen erzeugen.

Um das Gleiche im Notebook zu erstellen, müssen Sie MathJax
für die Verarbeitung von LaTeX-Code ausschalten, indem Sie
die "mathjax avoid list" benutzen. Graphen werden in einer
``tikzpicture`` Umgebung eingebunden, das ist also eine gute Wahl
für die Zeichenkette für die Ausschlussliste. Jetzt sollte
``view(graphs.CompleteGraph(4))`` in einem Arbeitsblatt
eine .pdf Datei mit pdflatex erstellen, mit dem
``convert`` Werkzeug eine .png Grafik erstellen und in die Ausgabezelle
des Arbeitsblatts einfügen.
Die folgenden Befehle veranschaulichen die Schritte einen Graphen
mittels LaTeX in einem Notebook darzustellen. ::

    sage: from sage.graphs.graph_latex import setup_latex_preamble
    sage: setup_latex_preamble()
    sage: latex.extra_preamble() # random - depends on system's TeX installation
    '\\usepackage{tikz}\n\\usepackage{tkz-graph}\n\\usepackage{tkz-berge}\n'
    sage: latex.engine('pdflatex')
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')  # not tested
    sage: latex.mathjax_avoid_list()                      # not tested
    ['tikz', 'tikzpicture']

Beachten Sie, dass es eine Vielzahl von Optionen gibt, die die
Darstellung des Graphen in LaTeX mit ``tkz-graph`` beeinflussen. Auch
das wiederrum ist nicht Ziel dieses Abschnitts. Sehen Sie sich hierfür
den Abschnitt "LaTeX-Optionen für Graphen" aus dem Handbuch für
weitere Anleitungen und Details an.

.. _sec-system-wide-tex:

Eine vollfunktionsfähige TeX-Installation
=========================================
Viele der erweiterten Integrationsmöglichkeiten von
TeX in Sage benötigen eine systemweite Installation von TeX.
Viele Linuxdistributionen bieten bereits TeX-Pakete basierend auf
TeX-live, für OSX gibt es TeXshop und für Windows MikTeX.
Das ``convert`` Werkzeug ist Teil der
`ImageMagick <http://www.imagemagick.org/>`_ Suite (welche ein
Paket oder zumindest ein simpler Download sein sollte). Die drei
Programme ``dvipng``, ``ps2pdf``, und ``dvips`` sind wahrscheinlich
bereits Teil Ihrer TeX Distribution.  Die ersten beiden sollten
auch von http://sourceforge.net/projects/dvipng/ als Teil von
`Ghostscript <http://www.ghostscript.com/>`_ bezogen werden können.

Um kombinatorische Graphen darstellen zu können,
benötigt man eine aktuelle Version der PGF-Bibliothek,
die Datei ``tkz-graph.sty`` (verfügbar
unter https://www.ctan.org/pkg/tkz-graph), die Datei
``tkz-arith.sty`` und vielleicht auch ``tkz-berge.sty``, beide verfügbar
auf der Seite https://www.ctan.org/pkg/tkz-berge.

Externe Programme
=================

Es sind drei Programme verfügbar um TeX weiter in Sage zu integrieren.
Das erste ist sagetex. Eine kurze Beschreibung von sagetex wäre: Es ist
eine Sammlung von TeX-Makros, die es einem LaTeX-Dokument erlauben
Anweisungen einzubinden, mit denen Sage genutzt wird um verschiedene
Objekte zu berechnen und/oder mittels eingebauter ``latex()``-Funktion darzustellen.
Als Zwischenschritt zum Kompilieren eines LaTeX-Dokuments werden also
alle Berechnungs- oder LaTeX-Formatierungseigenschaften von Sage automatisch genutzt.
Als Beispiel hierfür kann in einer mathematischen Betrachtung die korrekte Reihenfolge
von Fragen und Antworten beibehalten werden, indem sagetex dazu genutzt wird Sage die einen
aus den anderen berechnen zu lassen. Siehe hierfür auch :ref:`sec-sagetex`
