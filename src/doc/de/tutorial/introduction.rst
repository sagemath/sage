**********
Einleitung
**********

Um dieses Tutorial vollständig durchzuarbeiten sollten 3 bis 4 Stunden
genügen. Sie können es im HTML oder PDF Format
lesen oder im Sage Notebook. Dazu klicken Sie zuerst auf ``Help`` und
``Tutorial``, um innerhalb von Sage interaktiv mit dem
Tutorial zu arbeiten.

Obwohl große Teile von Sage mithilfe von Python implementiert sind,
ist kein tieferes Verständnis von Python notwendig um dieses Tutorial
lesen zu können. Sie werden Python zu einem gewissen Zeitpunkt lernen
wollen (Python kann sehr viel Spass bereiten) und es gibt viele
ausgezeichnete freie Quellen: das Python-Anfängerhandbuch [PyB]_
listet viele Optionen auf.
Wenn Sie nur kurz etwas in Sage ausprobieren möchten, ist dieses
Tutorial der richtige Ort um damit anzufangen. Zum Beispiel:

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223

    sage: A = matrix(4,4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]

    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)

    sage: m = matrix(ZZ,2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]

    sage: E = EllipticCurve([1,2,3,4,5]);
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1

    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    36/(20*sqrt(73) + 36*I*sqrt(3) + 27)
    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k,30)      # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

.. _installation:

Installation
============

Falls Sie Sage auf Ihrem Computer nicht installiert haben und nur ein
paar Befehle ausführen möchten, können Sie es online unter
http://sagecell.sagemath.org benutzen.

Schauen Sie sich den `Sage Installation Guide <http://doc.sagemath.org/html/en/installation/index.html>`_ an, um Anleitungen
zur Installation von Sage auf Ihrem Computer zu erhalten.
Hier geben wir nur ein paar Kommentare ab.


#. Die herunterladbare Sage-Datei wurde nach der *batteries included*
   Philosophie zusammengestellt. In anderen Worten, obwohl Sage
   Python, IPython, PARI, GAP, Singular, Maxima, NTL, GMP und so
   weiter benutzt, müssen Sie diese Programme  nicht separat
   installieren, da diese in der Sage-Distribution enthalten
   sind. Jedoch müssen Sie, um bestimmte Sage Zusatzfunktionen, zum
   Beispiel Macaulay oder KASH, nutzen zu können, die relevanten
   Programme auf ihrem Computer schon installiert haben.

#. Die vorkompilierte Binärversion von Sage (zu finden auf der
   Sage-Webseite) ist vielleicht einfacher und
   schneller zu installieren als die Quellcode-Version. Sie müssen
   die Datei nur entpacken und das Kommando ``sage`` ausführen.

#. Falls Sie das SageTeX-Paket benutzen möchten (mit welchem Sie
   die Ergebnisse von Sage Berechnungen in eine LaTeX-Datei
   einbauen können), müssen Sie SageTeX Ihrer TeX-Distribution bekannt
   machen. Um dies zu tun, lesen Sie den Abschnitt `Make SageTeX known
   to TeX <http://doc.sagemath.org/html/en/installation/sagetex.html>`_ im
   Sage Installation Guide
   (`Dieser Link <../../en/installation/index.html>`_ sollte Sie zu
   eine lokalen Kopie des Installation Guides führen). Es ist ziemlich
   einfach; Sie müssen  nur eine Umgebungsvariable setzen oder eine
   einzige Datei in ein Verzeichnis kopieren, welches TeX durchsucht.

   Die Dokumentation für SageTeX befindet sich in
   ``$SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/``, wobei
   "``$SAGE_ROOT``" auf das Verzeichnis zeigt, in welches Sie Sage
   installiert haben, zum Beispiel ``/opt/sage-9.6``.

Wie man Sage benutzen kann
==========================

Sie können Sage auf verschiedene Weise benutzen.

-  **graphisches Notebook-Interface:** rufen Sie ``sage -n jupyter`` auf; lesen Sie
   `Jupyter documentation on-line <https://jupyter-notebook.readthedocs.io/en/latest/notebook.html>`_,

-  **interaktive Kommandozeile:** lesen Sie :ref:`chapter-interactive_shell`,

-  **Programme:** Indem Sie interpretierte und kompilierte Programme in
   Sage schreiben (lesen Sie :ref:`section-loadattach` und :ref:`section-compile`), und

-  **Skripte:** indem Sie eigenständige Pythonskripte schreiben, welche
   die Sage-Bibliothek benutzen (lesen Sie :ref:`section-standalone`).


Langfristige Ziele von Sage
=============================

-  **nützlich**: Sages Zielgruppen sind Mathematikstudenten (von der
   Schule bis zur Universität), Lehrer und forschende
   Mathematiker. Das Ziel ist es, Software bereitzustellen, die benutzt
   werden kann, um mathematische Konstruktionen in der Algebra,
   Geometrie, Zahlentheorie, Analysis, Numerik, usw. zu erforschen und
   mit ihnen zu experimentieren. Sage hilft dabei, einfacher mit
   mathematischen Objekten experimentieren zu können.

-  **effizient:** Schnell sein. Sage benutzt hochoptimierte
   ausgereifte Software wie GMP, PARI, GAP und NTL, und ist somit bei
   vielen Aufgaben sehr schnell.

-  **frei und Open-Source:** Der Quellcode muss frei verfügbar und
   lesbar sein, damit Benutzer verstehen können, was das System gerade
   macht, und es einfacher erweitern zu können. Genauso wie
   Mathematiker ein tieferes Verständnis eines Theorems erlangen,
   indem sie den Beweis sorgfältig lesen oder zumindest überfliegen,
   sollten Leute, die Berechnungen durchführen, verstehen, wie die
   Berechnungen zustande kommen, indem sie den dokumentierten
   Quellcode lesen. Falls Sie Sage verwenden, um Berechnungen für ein
   Paper durchzuführen, welches Sie veröffentlichen, können Sie
   sicher sein, dass Ihre Leser immer freien Zugang zu Sage und
   seinem Quellcode haben und Sie dürfen sogar Ihre SageMath Version
   archivieren und weiterverteilen.

-  **einfach zu kompilieren:** Sage sollte für GNU/Linux, Mac OS X und
   Windowsbenutzer einfach aus dem Quellcode kompiliert werden können.

-  **kooperativ** Sage stellt robuste Schnittstelle zu vielen anderen
   Computeralgebrasystemen, einschließlich PARI, GAP, Singular, Maxima,
   KASH, Magma, Maple und Mathematica zur Verfügung. Sage ist dazu
   gedacht, bestehende Mathematik-Software zu vereinheitlichen und zu erweitern.

-  **gut dokumentiert:** Es gibt ein Tutorial, einen Programmierguide,
   ein Referenzhandbuch und Howtos mit zahlreichen Beispielen und
   Erläuterungen der dahinterstehenden Mathematik.

-  **erweiterbar:** Es ist möglich, neue Datentypen zu definieren oder
   von eingebauten Typen abzuleiten und Code vieler verschiedener Sprachen zu benutzen.

-  **benutzerfreundlich**: Es sollte einfach sein zu verstehen, welche
   Funktionalität für ein bestimmtes Objekt zur Verfügung gestellt
   wird und die Dokumentation und den Quellcode zu betrachten. Weiterhin sollte ein
   hochwertiger Benutzersupport erreicht werden.

