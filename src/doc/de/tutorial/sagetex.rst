.. _sec-sagetex:

**************
SageTeX nutzen
**************

Das SageTeX Paket ermöglicht es Ihnen die Ergebnisse von Sage Berechnungen
direkt in ein LaTeX-Dokument zu setzen. Es wird standardmäßig mit Sage
installiert. Um es zu nutzen müssen Sie es lediglich in Ihrem lokalen
TeX-System "installieren", wobei "installieren" hier eine einzige Datei
kopieren bedeutet. Siehe hierfür auch :ref:`installation` in diesem
Tutorial und den Abschnitt "Make SageTeX known to TeX" des `Sage installation guide
<http://doc.sagemath.org/html/en/installation/index.html>`_ (`dieser Link
<../../en/installation/index.html>`_ sollte Sie zu einer lokalen Kopie der
Installationsanleitung führen) um weitere Informationen zu erhalten.

Hier stellen wir ein sehr kurzes Beispiel vor wie man SageTeX nutzt.
Die komplette Dokumentation finden Sie unter :file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`,
wobei ``SAGE_ROOT`` das Installationsverzeichnis von Sage ist. Dieses Verzeichnis
enthält die Dokumentation, eine Beispieldatei und einige nützliche Python Skripte.

Um zu sehen wie SageTeX funktioniert, folgen Sie den Anweisungen zur Installation von
SageTeX (in :ref:`installation`) und kopieren Sie den folgenden Text in eine Datei namens -
sagen wir ``st_example.tex``:

.. warning::

  Der folgende Text wird mehrere Fehler bezüglich unbekannten Kontrollsequenzen
  anzeigen, wenn Sie ihn in der "live" Hilfe ansehen. Nutzen Sie stattdessen
  die statische Version um den korrekten Text zu sehen.

.. code-block:: latex

    \documentclass{article}
    \usepackage{sagetex}

    \begin{document}

    Wenn Sie Sage\TeX nutzen, können Sie Sage nutzen um Dinge auszurechen und
    sie direkt in ein \LaTeX{} Dokument zu setzen. Zum Beispiel gibt es
    $\sage{number_of_partitions(1269)}$ ganzzahlige Partitionen von $1269$.
    Sie müssen die Zahl nicht selbst ausrechnen, oder aus einem anderen
    Programm herauskopieren.

    Hier ein wenig Sage Code:

    \begin{sageblock}
        f(x) = exp(x) * sin(2*x)
    \end{sageblock}

    Die zweite Ableitung von $f$ ist

    \[
      \frac{\mathrm{d}^{2}}{\mathrm{d}x^{2}} \sage{f(x)} =
      \sage{diff(f, x, 2)(x)}.
    \]

    Hier ein Plot von $f$ von $-1$ bis $1$:

    \sageplot{plot(f, -1, 1)}

    \end{document}

Lassen Sie LaTeX ganz normal über ``st_example.tex`` laufen. Beachten Sie dabei, dass LaTeX
sich über einige Dinge beschwert, z.B.::

    Package sagetex Warning: Graphics file
    sage-plots-for-st_example.tex/plot-0.eps on page 1 does not exist. Plot
    command is on input line 25.

    Package sagetex Warning: There were undefined Sage formulas and/or
    plots. Run Sage on st_example.sage, and then run LaTeX on
    st_example.tex again.

Beachten Sie, dass zusätzlich zu den Dateien, die LaTeX normalerweise produziert
noch eine Datei ``st_example.sage`` erscheint. Das ist das Sage Skript, das
erstellt wurde als Sie LaTeX mit ``st_example.tex`` aufgerufen haben. Wie Ihnen die
Warnmeldung mitteilte sollten Sie Sage über die Datei ``st_example.sage`` laufen lassen,
also tun Sie das bitte. Ihnen wird gesagt werden, dass Sie LaTeX erneut über die Datei
``st_example.tex`` laufen lassen sollen; bevor Sie dies tun beachten Sie, dass eine neue
Datei namens ``st_example.sout`` von Sage erstellt wurde. Diese Datei enthält die Ergebnisse
von Sages Berechnungen in einem Format, das LaTeX nutzen kann um es in Ihren Text einzufügen.
Ein neues Verzeichnis mit einer .eps Datei Ihres Plots wurde ebenfalls erstellt.
Lassen Sie LaTeX nun erneut laufen, und Sie werden sehen, dass alles was Sage berechnet und
geplottet hat nun in Ihrem Dokument erscheint.

Die verschiednenen verwendeten Makros sollten einfach zu verstehen sein.
Eine ``sageblock`` Umgebung setzt Ihren Code unverändert und führt ihn auch
aus wenn Sie Sage laufen lassen. Wenn Sie etwa ``\sage{foo}`` schreiben, wird
das Ergebnis des Aufrufs ``latex(foo)`` (in Sage) in Ihrem Dokument erscheinen.
Plot-Befehle sind etwas komplizierter, aber in Ihrer einfachsten Form fügt
``\sageplot{foo}`` das Bild ein, das Sie erhalten wenn Sie ``foo.save('filename.eps')``
in Sage aufrufen würden.

Grundsätzlich gilt:

    - lassen Sie LaTeX über Ihre .tex Datei laufen;
    - lassen Sie Sage über die neu generierte .sage Datei laufen;
    - lassen Sie LaTeX erneut laufen.

Sie können das Aufrufen von Sage weglassen, wenn Sie keine Änderung
an den Sage Befehlen in Ihrem Dokument vorgenommen haben.

Es gibt noch viel mehr über SageTeX zu sagen, aber da sowohl Sage alsauch
LaTeX komplexe und mächtige Werkzeuge sind, sollten Sie die Dokumentation
über SageTeX in :file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex` lesen.
