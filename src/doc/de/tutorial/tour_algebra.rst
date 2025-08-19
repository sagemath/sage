Elementare Algebra und Analysis
===============================

Sage kann viele zur elementaren Algebra und Analysis gehörende
Probleme lösen. Zum Beispiel: Lösungen von Gleichungen finden,
Differentiation, Integration, und Laplace-Transformationen
berechnen. Lesen Sie die `Sage Constructions
<http://doc.sagemath.org/html/en/constructions/>`_ Dokumentation um
weitere Beispiele zu finden.

Lösen von Gleichungen
---------------------

Gleichungen exakt lösen
~~~~~~~~~~~~~~~~~~~~~~~

Die ``solve`` Funktion löst Gleichungen. Legen Sie zunächst Variablen
an, bevor Sie diese benutzen; Die Argumente von ``solve`` sind eine
Gleichung (oder ein System von Gleichungen) zusammen mit den
Variablen, nach welchen Sie auflösen möchten:

::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

Sie können eine Gleichung nach einer Variablen, in Abhängigkeit von den
anderen, auflösen:

::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]

Sie können auch nach mehreren Variablen auflösen:

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]

Das folgende Beispiel, in dem Sage benutzt wird um ein System von
nichtlinearen Gleichungen zu lösen, stammt von Jason Grout. Zunächst
lösen wir das System symbolisch:

::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3], [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]

Um eine numerische Approximation der Lösungen zu erhalten können Sie
stattdessen wie folgt vorgehen:

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]

(Die Funktion ``n`` gibt eine numerische Approximation zurück, ihr
Argument ist die Anzahl der Bits an Genauigkeit.)

Gleichungen numerisch lösen
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Oftmals kann ``solve`` keine exakte Lösung der angegebenen Gleichung
bzw. Gleichungen finden. Wenn dies passiert können Sie ``find_root``
verwenden um eine numerische Approximation zu finden. Beispielsweise
gibt ``solve`` bei folgender Gleichung nichts brauchbares zurück::

    sage: theta = var('theta')
    sage: solve(cos(theta)==sin(theta), theta)
    [sin(theta) == cos(theta)]

Wir können jedoch ``find_root`` verwenden um eine Lösung der obigen
Gleichung im Bereich :math:`0 < \phi < \pi/2` zu finden::

    sage: phi = var('phi')
    sage: find_root(cos(phi)==sin(phi),0,pi/2)
    0.785398163397448...

Differentiation, Integration, etc.
----------------------------------

Sage weiß wie man viele Funktionen differenziert und integriert. Zum
Beispiel können Sie folgendes eingeben um :math:`\sin(u)` nach
:math:`u` abzuleiten:

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

Um die vierte Ableitung :math:`\sin(x^2)` zu berechnen:

::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)

Um die partiellen Ableitungen von :math:`x^2+17y^2` nach `x`
beziehungsweise `y` zu berechnen:

::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)
    34*y

Wir machen weiter mit Integralen, sowohl bestimmt als auch
unbestimmt. Die Berechnung von :math:`\int x\sin(x^2)\, dx` und
:math:`\int_0^1 \frac{x}{x^2+1}\, dx`:

::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

Die Partialbruchzerlegung von :math:`\frac{1}{x^2-1}`:

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)

.. _section-systems:

Lösen von Differentialgleichungen
---------------------------------

Sie können Sage verwenden um gewöhnliche Differentialgleichungen zu
berechnen. Die Gleichung :math:`x'+x-1=0` berechnen Sie wie folgt:

::

    sage: t = var('t')    # definiere die Variable t
    sage: x = function('x')(t)   # definiere x als Funktion dieser Variablen
    sage: DE = diff(x, t) + x - 1
    sage: desolve(DE, [x,t])
    (_C + e^t)*e^(-t)

Dies benutzt Sages Schnittstelle zu Maxima [Max]_, daher kann sich die
Ausgabe ein wenig von anderen Ausgaben in Sage unterscheiden. In
diesem Fall wird mitgeteilt, dass :math:`x(t) = e^{-t}(e^{t}+c)`
die allgemeine Lösung der Differentialgleichung ist.

Sie können auch Laplace-Transformationen berechnen:
Die Laplace-Transformation von :math:`t^2e^t -\sin(t)` wird wie folgt
berechnet:

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3

Hier ist ein komplizierteres Beispiel. Die Verschiebung des
Gleichgewichts einer verkoppelten Feder, die an der linken Wand
befestigt ist,

::

    |------\/\/\/\/\---|Masse1|----\/\/\/\/\/----|Masse2|
             Feder1                  Feder2

wird durch dieses System der Differentialgleichungen zweiter Ordnung
modelliert,

.. math::

    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0

    m_2 x_2''+ k_2 (x_2-x_1) = 0,


wobei :math:`m_{i}` die Masse des Objekts *i*, :math:`x_{i}` die
Verschiebung des Gleichgewichts der Masse *i* und :math:`k_{i}` die
Federkonstante der Feder *i* ist.

**Beispiel:** Benutzen Sie Sage um das obige Problem mit folgenden
Werten zu lösen:
:math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`,
:math:`k_{2}=2`, :math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`,
:math:`x_{2}(0)=3`, :math:`x_{2}'(0)=0`.

Lösung: Berechnen Sie die Laplace-Transformierte der ersten Gleichung
(mit der Notation :math:`x=x_{1}`, :math:`y=x_{2}`):

::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

Das ist schwierig zu lesen, es besagt jedoch, dass

.. math:: -2x'(0) + 2s^2\cdot X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0


(wobei die Laplace-Transformierte der Funktion mit kleinem
Anfangsbuchstaben :math:`x(t)` die Funktion mit großem
Anfangsbuchstaben :math:`X(s)` ist). Berechnen Sie die
Laplace-Transformierte der zweiten Gleichung:

::

    sage: de2 = maxima("diff(y(t),t, 2) + 2*y(t) - 2*x(t)")
    sage: lde2 = de2.laplace("t","s"); lde2.sage()
    s^2*laplace(y(t), t, s) - s*y(0) - 2*laplace(x(t), t, s) + 2*laplace(y(t), t, s) - D[0](y)(0)

Dies besagt

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.


Setzen Sie die Anfangsbedingungen für :math:`x(0)`, :math:`x'(0)`,
:math:`y(0)` und :math:`y'(0)` ein, und lösen die beiden Gleichungen,
die Sie so erhalten:

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s]
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
      Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

Berechnen Sie jetzt die inverse Laplace-Transformierte um die Antwort
zu erhalten:

::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)

Also ist die Lösung:

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).


Die kann folgenderweise parametrisiert geplottet werden:

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),
    ....:     (t, 0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

Die einzelnen Komponenten können so geplottet werden:

::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), (t,0, 2*pi), rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), (t,0, 2*pi), rgbcolor=hue(0.6))
    sage: show(p1 + p2)

Um mehr über das Plotten zu erfahren lesen Sie :ref:`section-plot`. Lesen
Sie Abschnitt 5.5 von [NagleEtAl2004]_ um weitere Informationen über
Differentialgleichungen zu erhalten.


Das Euler-Verfahren zur Lösung von Systemen von Differentialgleichungen
-----------------------------------------------------------------------

Im nächsten Beispiel illustrieren wir das Euler-Verfahren für ODEs erster
und zweiter Ordnung. Wir rufen zunächst die grundlegende Idee für
Differentialgleichungen erster Ordnung in Erinnerung. Sei ein
Anfangswertproblem der Form

.. math::

    y'=f(x,y), \quad y(a)=c,

gegeben. Wir möchten eine Approximation des Wertes der Lösung bei
:math:`x=b` mit :math:`b>a` finden.

Machen Sie sich anhand der Definition der Ableitung klar, dass

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},


wobei :math:`h>0` vorgegeben und klein ist. Zusammen mit der
Differentialgleichung gibt dies :math:`f(x,y(x))\approx
\frac{y(x+h)-y(x)}{h}`. Jetzt lösen wir nach :math:`y(x+h)` auf:

.. math::   y(x+h) \approx y(x) + h\cdot f(x,y(x)).


Wenn wir :math:`h\cdot f(x,y(x))` den "Korrekturterm", :math:`y(x)`
den "alten Wert von `y`" und :math:`y(x+h)` den "neuen Wert von `y`"
nennen, kann diese Approximation neu ausgedrückt werden als:

.. math::   y_{new} \approx y_{old} + h\cdot f(x,y_{old}).


Wenn wir das Intervall von `a` bis `b` in `n` Teilintervalle
aufteilen, so dass :math:`h=\frac{b-a}{n}` gilt, können wir die
Information in folgender Tabelle festhalten.

============== =======================   =====================
:math:`x`      :math:`y`                 :math:`h\cdot f(x,y)`
============== =======================   =====================
:math:`a`      :math:`c`                 :math:`h\cdot f(a,c)`
:math:`a+h`    :math:`c+h\cdot f(a,c)`         ...
:math:`a+2h`   ...
...
:math:`b=a+nh` ???                             ...
============== =======================   =====================


Unser Ziel ist zeilenweise alle leeren Einträge der Tabelle
auszufüllen, bis wir den Eintrag ??? erreichen, welcher die
Approximation des Euler-Verfahrens für :math:`y(b)` ist.

Die Idee für Systeme von ODEs ist ähnlich.

**Beispiel:** Approximiere :math:`z(t)`, mit 4 Schritten der
 Eulermethode numerisch bei :math:`t=1` , wobei :math:`z''+tz'+z=0`,
 :math:`z(0)=1` und :math:`z'(0)=0` ist.

Wir müssen die ODE zweiter Ordnung auf ein System von zwei
Differentialgleichungen erster Ordnung reduzieren (wobei :math:`x=z`,
:math:`y=z'`) und das Euler-Verfahren anwenden:

::

    sage: t,x,y = PolynomialRing(RealField(10),3,"txy").gens()
    sage: f = y; g = -x - y * t
    sage: eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
          t                x            h*f(t,x,y)                y       h*g(t,x,y)
          0                1                  0.00                0           -0.25
        1/4              1.0                -0.062            -0.25           -0.23
        1/2             0.94                 -0.12            -0.48           -0.17
        3/4             0.82                 -0.16            -0.66          -0.081
          1             0.65                 -0.18            -0.74           0.022

Also ist :math:`z(1)\approx 0.75`.

Wir können auch die Punkte :math:`(x,y)` plotten um ein ungefähres
Bild der Kurve zu erhalten. Die Funktion ``eulers_method_2x2_plot``
macht dies; um sie zu benutzen, müssen wir die Funktionen  `f` und `g`
definieren, welche ein Argument mit drei Koordinaten (`t`, `x`, `y`)
erwarten.

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

Zu diesem Zeitpunkt enthält ``P`` die beiden Plots ``P[0]`` (der Plot
von `x` nach `t`) und ``P[1]`` (der Plot von `y` nach `t`). Wir können
beide wie folgt anzeigen:

.. link

::

    sage: show(P[0] + P[1])

(Um mehr über das Plotten zu erfahren, lesen Sie :ref:`section-plot`.)

Spezielle Funktionen
--------------------

Mehrere orthogonale Polynome und spezielle Funktionen sind
implementiert, wobei sowohl PARI [GP]_ als auch Maxima [Max]_
verwendet wird. Sie sind in den dazugehörigen Abschnitten ("Orthogonal polynomials"
beziehungsweise "Special functions") des Sage Referenzhandbuchs dokumentiert.

::

    sage: x = polygen(QQ, 'x')
    sage: chebyshev_U(2,x)
    4*x^2 - 1
    sage: bessel_I(1,1).n(250)
    0.56515910399248502720769602760986330732889962162109200948029448947925564096
    sage: bessel_I(1,1).n()
    0.565159103992485
    sage: bessel_I(2,1.1).n()
    0.167089499251049

Zum jetzigen Zeitpunkt, enthält Sage nur Wrapper-Funktionen für
numerische Berechnungen. Um symbolisch zu rechen, rufen Sie die
Maxima-Schnittstelle bitte, wie im folgenden Beispiel, direkt auf

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'
