s.. _section-functions-issues:

Häufige Probleme mit Funktionen
===============================

Das Definieren von Funktionen kann in mancher Hinsicht verwirrend sein
(z.B. beim Ableiten oder Plotten). In diesem Abschnitt versuchen wir
die relevanten Probleme anzusprechen.

Nun erläutern wir verschiedene Möglichkeiten Dinge zu definieren, die das
Recht haben könnten "Funktionen" genannt zu werden:

1. Definition einer Python-Funktion wie in :ref:`section-functions`
beschrieben. Diese Funktionen können geplottet, aber nicht abgeleitet
oder integriert werden.

::

       sage: def f(z): return z^2
       sage: type(f)
       <... 'function'>
       sage: f(3)
       9
       sage: plot(f, 0, 2)
       Graphics object consisting of 1 graphics primitive

Beachten Sie die Syntax in der letzten Zeile. Falls Sie stattdessen
``plot(f(z), 0, 2)`` verwenden, erhalten Sie einen Fehler, da ``z``
eine Dummy-Variable in der Definition von ``f`` ist, und außerhalb
dieser nicht definiert ist. In der Tat gibt sogar nur ``f(z)`` einen
Fehler zurück. Das Folgende funktioniert in diesem Fall, obwohl es im
Allgemeinen Probleme verursachen kann und deshalb vermieden werden
sollte. (Beachten Sie unten den 4. Punkt)

.. link

::

       sage: var('z')   # z wird als Variable definiert
       z
       sage: f(z)
       z^2
       sage: plot(f(z), 0, 2)
       Graphics object consisting of 1 graphics primitive

Nun ist `f(z)`` ein symbolischer Ausdruck. Dies ist unser nächster Stichpunkt
unserer Aufzählung.

2. Definition eines  "aufrufbaren symbolischen Ausdrucks".  Diese
können geplottet, differenziert und integriert werden.

::

       sage: g(x) = x^2
       sage: g        # g bildet x auf x^2 ab
       x |--> x^2
       sage: g(3)
       9
       sage: Dg = g.derivative(); Dg
       x |--> 2*x
       sage: Dg(3)
       6
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(g, 0, 2)
       Graphics object consisting of 1 graphics primitive

Beachten Sie, dass während ``g`` ein aufrufbarer symbolischer Ausdruck
ist, ``g(x)`` ein verwandtes aber unterschiedliches Objekt ist,
welches auch geplottet, differenziert, usw. werden kann - wenn auch
mit einigen Problemen: Lesen Sie den 5. Stichpunkt unterhalb, um eine
Erläuterung zu erhalten.

.. link

::

       sage: g(x)
       x^2
       sage: type(g(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: g(x).derivative()
       2*x
       sage: plot(g(x), 0, 2)
       Graphics object consisting of 1 graphics primitive

3. Benutzung einer vordefinierten 'trigonometrischen Sage-Funktion'.
Diese können mit ein wenig Hilfestellung differenziert und integriert
werden.

::

       sage: type(sin)
       <class 'sage.functions.trig.Function_sin'>
       sage: plot(sin, 0, 2)
       Graphics object consisting of 1 graphics primitive
       sage: type(sin(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(sin(x), 0, 2)
       Graphics object consisting of 1 graphics primitive

Alleinestehend kann ``sin`` nicht differenziert werden, zumindest nicht
um ``cos`` zu erhalten.

::

       sage: f = sin
       sage: f.derivative()
       Traceback (most recent call last):
       ...
       AttributeError: ...

``f = sin(x)`` anstelle von ``sin`` zu benutzen funktioniert, aber
es ist wohl noch besser ``f(x) = sin(x)`` zu benutzen, um einen
aufrufbaren symbolischen Ausdruck zu definieren.

::

       sage: S(x) = sin(x)
       sage: S.derivative()
       x |--> cos(x)

Hier sind ein paar häufige Probleme mit Erklärungen:

\4. Versehentliche Auswertung.

::

       sage: def h(x):
       ....:     if x<2:
       ....:         return 0
       ....:     else:
       ....:         return x-2

Das Problem: ``plot(h(x), 0, 4)`` zeichnet die Linie `y=x-2` und nicht
die mehrzeilige Funktion, welche durch ``h`` definiert wird.  Der
Grund? In dem Befehl ``plot(h(x), 0, 4)`` wird zuerst ``h(x)``
ausgewertet: Das bedeutet, dass ``x`` in die Funktion ``h`` eingesetzt
wird, was wiederum bedeutet, dass ``x<2`` ausgewertet wird.

.. link

::

       sage: type(x<2)
       <class 'sage.symbolic.expression.Expression'>

Wenn eine symbolische Gleichung ausgewertet wird, wie in der
Definition von ``h``, wird falls sie nicht offensichtlicherweise wahr
ist, False zurück gegeben. Also wird ``h(x)`` zu ``x-2`` ausgewertet
und dies ist die Funktion, die geplottet wird.

Die Lösung: verwenden Sie nicht ``plot(h(x), 0, 4)``; benutzen Sie stattdessen:

.. link

::

       sage: plot(h, 0, 4)
       Graphics object consisting of 1 graphics primitive

\5. Versehentliches Erzeugen einer Konstanten anstelle von einer Funktion.

::

       sage: f = x
       sage: g = f.derivative()
       sage: g
       1

Das Problem: ``g(3)``, zum Beispiel, gibt folgenden Fehler zurück:
"ValueError: the number of arguments must be less than or equal to 0."

.. link

::

       sage: type(f)
       <class 'sage.symbolic.expression.Expression'>
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>

``g`` ist keine Funktion, es ist eine Konstante, hat also keine
zugehörigen Variablen, und man kann in sie nichts einsetzen.

Die Lösung: Es gibt mehrere Möglichkeiten.

- Definieren Sie ``f`` anfangs als symbolischen Ausdruck.

::

         sage: f(x) = x        # statt 'f = x'
         sage: g = f.derivative()
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- Oder mit der ursprünglichen Definition von ``f``, definieren Sie
  ``g`` als symbolischen Ausdruck.

::

         sage: f = x
         sage: g(x) = f.derivative()  # statt 'g = f.derivative()'
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- Oder mit den ursprünglichen Definitionen von ``f`` and ``g``, geben
  Sie die Variable an, in diese Sie den Wert einsetzen.

::

         sage: f = x
         sage: g = f.derivative()
         sage: g
         1
         sage: g(x=3)    # statt 'g(3)'
         1

Schließlich ist hier noch eine Möglichkeit den Unterschied zwischen der
Ableitung von ``f = x`` und der von ``f(x) = x`` zu erkennen:

::

       sage: f(x) = x
       sage: g = f.derivative()
       sage: g.variables()  # Die in g präsenten Variablen
       ()
       sage: g.arguments()  # Die Argumente die in g gesteckt werden können
       (x,)
       sage: f = x
       sage: h = f.derivative()
       sage: h.variables()
       ()
       sage: h.arguments()
       ()

Wie dieses Beispiel verdeutlichen sollte, nimmt ``h`` keine Argumente
an, und deshalb gibt ``h(3)`` einen Fehler zurück.
