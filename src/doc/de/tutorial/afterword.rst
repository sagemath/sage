********
Nachwort
********

Warum Python?
=============

Vorteile von Python
-------------------

Sage ist hautsächlich in der Programmiersprache Python implementiert (siehe [Py]_).
Jedoch ist Code, bei dem Geschwindigkeit ausschlaggebend ist, in einer
kompilierten Sprache implementiert. Python hat folgende Vorteile:


-  **Speichern von Objekten** wird in Python gut unterstützt. Für das
   Speichern von (nahezu) beliebigen Objekten auf Festplatten oder in
   Datenbanken sind in Python weitgehende Hilfsmittel vorhanden.

-  Exzellente Unterstütztung für die **Dokumentation** von Funktionen
   und Paketen im Quellcode, einschließlich der automatischen
   Erstellung der Dokumentation und automatisches Testen aller
   Beispiele. Die Beispiele werden regelmäßig automatisch getestet und
   es wird garantiert, dass sie wie angegeben funktionieren.

-  **Speicherverwaltung**: Python besitzt nun einen gut durchdachten
   und robusten Speicherverwalter und einen Speicherbereiniger, der
   zirkuläre Referenzen korrekt behandelt und lokale Variablen in
   Dateien berücksichtigt.

-  Python besitzt mittlerweile **viele Pakete**, die für Sagenutzer
   sehr reizvoll sein könnten: numerische Analysis und lineare
   Algebra, 2D und 3D Visualisierungen, Vernetzungen (für verteilte
   Berechnungen und Server, z.B. mithilfe von twisted),
   Datenbankunterstützung, usw.

-  **Portabilität:** Python kann auf den meisten Systemen
   unkompliziert, innerhalb von Minuten aus dem Quellcode kompiliert
   werden.

-  **Fehlerbehandlung:** Python besitzt ein ausgeklügeltes und wohl
   durchdachtes System für die Behandlung von Ausnahmebedingungen,
   mit dem Programme sinnvoll weiterarbeiten können, sogar wenn bei
   ihrem Aufruf Fehler auftreten.

-  **Debugger:** Python beinhaltet einen Debugger. Folglich kann der
   Benutzer, falls der Code aus irgendeinem Grund fehlschlägt, auf
   eine ausgiebige Stack-Ablaufverfolgung zugreifen, den Zustand
   aller relevanter Variablen betrachten, und sich auf dem Stack nach
   oben oder unten bewegen.

-  **Profiler:** Es gibt einen Python-Profiler, welcher Code
   ausführt und einen Bericht erstellt, in dem detailliert
   aufgestellt wurde wie oft und wie lange jede Funkion aufgerufen
   wurde.

-  **Eine Sprache:** Anstatt eine **neue Sprache** für mathematische
   Software zu schreiben, wie es für Magma, Maple, Mathematica, Matlab,
   GP/PARI, GAP, Macaulay 2, Simath, usw. gemacht wurde, benutzen wir
   die Programmiersprache Python, eine beliebte Programmiersprache, die
   von hunderten begabten Softwareingenieuren rege weiterentwickelt und
   optimiert wird. Python ist eine bedeutete Open-Source
   Erfolgsgeschichte mit einem ausgereiften Entwicklungsprozess. (siehe [PyDev]_).


.. _section-mathannoy:

Der Pre-Parser: Unterschiede zwischen Sage und Python
-----------------------------------------------------

Aus mathematischer Sicht kann Python in verschiedener Weise verwirrend
sein, also verhält sich Sage an manchen Stellen anders als Python.

-  **Notation für Exponentiation:** ``**`` versus ``^``. In Python
   bedeutet  ``^`` "xor", und nicht Exponentiation, also gilt in
   Python:

   ::

       >>> 2^8
       10
       >>> 3^2
       1
       >>> 3**2
       9

   Diese Benutzung von ``^`` kann merkwürdig erscheinen und sie ist
   ineffizient für mathematische Anwender, da die
   "Exklusives-Oder"-Funktion nur selten verwendet wird.
   Um dies zu beheben parst Sage alle Kommandozeilen bevor es diese zu
   Python weitergibt und ersetzt jedes Auftreten ``^``, das in keinem
   String vorkommt mit ``**``:

   ::

       sage: 2^8
       256
       sage: 3^2
       9
       sage: "3^2"
       '3^2'

-  **Integerdivision:** Der Pythonaudruck ``2/3`` verhält sich nicht
   so, wie es Mathematiker erwarten würden. In Python 2 wird, falls ``m`` und
   ``n`` Integer sind, auch ``m/n`` als Integer behandelt, es ist
   nämlich der Quotient von ``m`` geteilt durch ``n``. Daher ist
   ``2/3=0``.  Es wurde in der Pythoncommunity darüber geredet, ob in
   Python die Division geändert werden sollte, so dass ``2/3`` die
   Gleitkommazahl ``0.6666...`` zurückgibt und ``2//3`` das Ergebnis
   ``0`` hat.

   Wir berücksichtigen dies im Sage-Interpreter indem wir
   Integer-Literale mit  ``Integer( )`` versehen und die Division als
   Konstruktor für rationale Zahlen behandeln. Zum Beispiel:

   ::

       sage: 2/3
       2/3
       sage: (2/3).parent()
       Rational Field
       sage: 2//3
       0

-  **Große ganze Zahlen:** Python besitzt von Hause aus Unterstützung
   für beliebig große ganze Zahlen zusätzlich zu C-ints. Diese sind
   bedeutend langsamer als die von GMP zur Verfügung gestellten und sie
   haben die Eigenschaft, dass die mit einem ``L`` am Ende ausgegeben
   werden um sie von ints unterscheiden zu können (und dies wird sich
   in naher Zeit nicht ändern). Sage implementiert beliebig große
   Integers mit Hilfe der GMP C-Bibliothek, und diese werden ohne
   ``L`` ausgegeben.


Anstatt den Python-Interpreter zu verändern (wie es mache Leute für
interne Projekte getan haben), benutzen wir die Sprache Python
unverändert und haben einen Prä-Parser geschrieben, so dass sich
die Kommandozeilen-IPython-Version so verhält, wie es Mathematiker
erwarten würden. Dies bedeutet, dass bereits existierender Python-Code
in Sage so verwendet werden kann wie er ist. Man muss jedoch immernoch
die standardmäßigen Python-Regeln beachten, wenn man Pakete schreibt,
die in Sage importiert werden können.

(Um eine Python-Bibliothek zu installieren, die Sie zum Beispiel im
Internet gefunden haben, folgen Sie den Anweisungen, aber verwenden
sie ``sage -python`` anstelle von ``python``.  Oft bedeutet dies, dass
``sage -python setup.py install`` eingegeben werden muss.)


Ich möchte einen Beitrag zu Sage leisten. Wie kann ich dies tun?
================================================================

Falls Sie für Sage einen Beitrag leisten möchten, wird Ihre Hilfe hoch
geschätzt! Sie kann von wesentlichen Code-Beiträge bis zum Hinzufügen
zur Sage-Dokumention oder zum Berichten von Fehlern reichen.


Schauen Sie sich die Sage-Webseite an um Informationen für Entwickler
zu erhalten; neben anderen Dingen können Sie eine lange Liste nach
Priorität und Kategorie geordneter, zu Sage gehörender Projekte finden.
Auch der `Sage Developer's Guide <http://doc.sagemath.org/html/en/developer/>`_
beinhaltet hilfreiche Informationen, und Sie können der ``sage-devel``
Google-Group beitreten.


Wie zitiere ich Sage?
=====================

Falls Sie ein Paper schreiben, das Sage verwendet, zitieren Sie bitte
die Berechnungen die Sie mithilfe von Sage durchgeführt haben, indem
Sie

::

    [Sage] SageMath, the Sage Mathematics Software System (Version 8.7),
           The Sage Developers, 2019, https://www.sagemath.org.

in Ihrem Literaturverzeichnis hinzufügen. (Ersetzen Sie hierbei 8.7 mit der von
Ihnen benutzten Version von Sage.) Versuchen Sie bitte weiterhin
festzustellen welche Komponenten von Sage in Ihrer Berechnung
verwendet wurden, z.B. PARI?, GAP?, Singular? Maxima? und zitieren Sie
diese Systeme ebenso. Falls Sie nicht sicher sind welche Software Ihre
Berechnung verwendet, können Sie dies gerne in der ``sage-devel``
Google-Gruppe fragen. Lesen Sie :ref:`section-univariate` um weitere
Information darüber zu erhalten.

------------

Falls Sie gerade das Tutorial vollständig durchgelesen haben, und noch
wissen wie lange Sie hierfür gebraucht haben, lassen Sie und dies bitte
in der ``sage-devel`` Google-Gruppe wissen.

Viel Spass mit Sage!
