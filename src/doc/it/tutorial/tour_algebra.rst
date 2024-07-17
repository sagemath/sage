Algebra di base e Analisi
=========================

Sage sa svolgere diversi calcoli legati all'algebra di base
ed all'analisi: per esempio, risoluzione di equazioni,
calcolo differenziale ed integrale e trasformate di Laplace.
Si veda la documentazione per le "Costruzioni di Sage" per
ulteriori esempi.

Risoluzione di equazioni
------------------------

La funzione ``solve`` risolve le equazioni. Per usarla,
bisogna anzitutto specificare alcune variabili; pertanto
gli argomenti di ``solve`` sono un'equazione (od un sistema
di equazioni), insieme con le variabili rispetto alle quali
risolvere:

::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

Si possono risolvere le equazioni rispetto ad una variabile in funzione
delle altre:

::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]

Si può anche risolvere rispetto a diverse variabili:

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]

Il seguente esempio dell'uso di Sage per risolvere un sistema di
equazioni non lineari è stato fornito da Jason Grout: per prima cosa,
si risolve il sistema simbolicamente:

::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3],
    [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]

Per una soluzione numerica, si può invece usare:

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]

(La funzione ``n`` scrive un'approssimazione numerica, e
l'argomento è il numero di bit di precisione.)

Differenziazione, Integrazione, etc.
------------------------------------

Sage è in grado di differenziare ed integrare molte funzioni. Per
esempio, per differenziare :math:`\sin(u)` rispetto a :math:`u`,
si procede come nelle righe seguenti:

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

Per calcolare la derivata quarta di :math:`\sin(x^2)`:

::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)

Per calcolare le derivate parziali di :math:`x^2+17y^2`
rispetto a *x* e *y*, rispettivamente:

::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)                                
    34*y

Passiamo agli integrali, sia indefiniti che definiti. Per calcolare
:math:`\int x\sin(x^2)\, dx` e
:math:`\int_0^1 \frac{x}{x^2+1}\, dx`

::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

Per calcolare la decomposizione in frazioni parziali di
:math:`\frac{1}{x^2-1}`:

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)

.. _section-systems:

Risoluzione di Equazioni Differenziali
--------------------------------------

Si può usare Sage per studiare le equazioni differenziali ordinarie.
Per risolvere l'equazione :math:`x'+x-1=0`:

::

    sage: t = var('t')    # definisce una variabile t
    sage: x = function('x')(t)   # definisce x come funzione di quella variabile
    sage: DE = diff(x,t) + x - 1
    sage: desolve(DE, [x,t])
    (_C + e^t)*e^(-t)

Questo metodo utilizza l'interfaccia di Sage per Maxima [Max]_, e così il suo
output può essere leggermente diverso dagli altri output di Sage. In questo caso,
risulta che la soluzione generale dell'equazione differenziale è
:math:`x(t) = e^{-t}(e^{t}+c)`.

Si può anche calcolare la trasformata di Laplace; la trasformata di Laplace di
:math:`t^2e^t -\sin(t)` è calcolata come segue:

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3

Il successivo è un esempio più articolato. Lo scostamento dall'equilibrio
(rispettivamente) per due molle accoppiate fissate ad un muro a sinistra

.. CODE-BLOCK:: text

    |------\/\/\/\/\---|massa1|----\/\/\/\/\/----|massa2|
             molla1                  molla2

è modellizzato dal sistema di equazioni differenziali del secondo ordine

.. math::

    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0
    m_2 x_2''+ k_2 (x_2-x_1) = 0,


dove :math:`m_{i}` è la massa dell'oggetto *i*, :math:`x_{i}` è
lo scostamento dall'equilibrio della massa *i*, e :math:`k_{i}`
è la costante elastica della molla *i*.

**Esempio:** Usare Sage per risolvere il problema precedente con
:math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`,
:math:`k_{2}=2`, :math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`,
:math:`x_{2}(0)=3`, :math:`x_{2}'(0)=0`.

Soluzione: Calcolare la trasformata di Laplace della prima equazione (con
la notazione :math:`x=x_{1}`, :math:`y=x_{2}`:

::

    sage: de1 = maxima("2*diff(x(t),t, 2) + 6*x(t) - 2*y(t)")
    sage: lde1 = de1.laplace("t","s"); lde1.sage()
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

Questo è di difficile lettura, ma dice che

.. math:: -2x'(0) + 2s^2*X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0


(dove la trasformata di Laplace di una funzione in minuscolo come
:math:`x(t)` è la funzione in maiuscolo :math:`X(s)`). Calcolare la
trasformata di Laplace della seconda equazione:

::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

che significa

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.


Imporre le condizioni iniziali per :math:`x(0)`, :math:`x'(0)`,
:math:`y(0)`, e :math:`y'(0)`, e risolvere le due equazioni
risultanti:

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s] 
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
    Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

Ora si calcola la trasformata inversa di Laplace per ottenere la risposta:

::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)

Pertanto, la soluzione è

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).


Essa può essere disegnata in forma parametrica usando

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),
    ....: (0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

Le singole componenti possono essere tracciate usando::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), 0, 2*pi, rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), 0, 2*pi, rgbcolor=hue(0.6))
    sage: show(p1 + p2)

BIBLIOGRAFIA: Nagle, Saff, Snider, Fundamentals of Differential
Equations, 6th ed, Addison-Wesley, 2004. (si veda § 5.5).

Metodo di Eulero per i sistemi di equazioni differenziali
---------------------------------------------------------

Nel prossimo esempio, si illustrerà il metodo di Eulero per le ODE
di primo e secondo ordine. Per prima cosa ricordiamo l'idea di base per
le equazioni di primo ordine. Dato un problema di Cauchy della forma

.. math::
    y'=f(x,y)
    y(a)=c 


si vuole trovare il valore approssimato della soluzione a
:math:`x=b` con :math:`b>a`.

Ricordando dalla definizione di derivata che

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},


dove :math:`h>0` è dato e piccolo. Questo e la DE insieme danno
give :math:`f(x,y(x))\approx
\frac{y(x+h)-y(x)}{h}`. Ora si risolve
per :math:`y(x+h)`:

.. math::   y(x+h) \approx y(x) + h*f(x,y(x)).


Se chiamiamo :math:`h f(x,y(x))` il "termine di correzione" (per mancanza
di un termine migliore), :math:`y(x)` il "vecchio valore di *y*", e
:math:`y(x+h)` il "nuovo valore di *y*", allora questa
approssimazione può essere espressa come

.. math::   y_{new} \approx y_{old} + h*f(x,y_{old}).


Se si spezza l'intervallo da *a* a *b* in *n* intervalli, dimodoché
:math:`h=\frac{b-a}{n}`, allora si possono registrare le informazioni per
questo metodo in una tabella.

============== ==================   ================
:math:`x`      :math:`y`            :math:`hf(x,y)`
============== ==================   ================
:math:`a`      :math:`c`            :math:`hf(a,c)`
:math:`a+h`    :math:`c+hf(a,c)`    ...
:math:`a+2h`   ...                   
...
:math:`b=a+nh` ???                  ...
============== ==================   ================  


L'obiettivo è riempire tutti gli spazi vuoti della tavella, una riga alla
volta, finché si arriva al valore ???, che è il
metodo di approssimazione di Eulero per :math:`y(b)`.

L'idea per sistemi di ODE è simile.

**Esempio:** Si approssimi numericamente :math:`z(t)` a :math:`t=1` usando 4
passi del metodo di Eulero, dove :math:`z''+tz'+z=0`,
:math:`z(0)=1`, :math:`z'(0)=0`.

Si deve ridurre l'ODE di secondo ordine ad un sistema di due equazioni del primo
ordine (usando :math:`x=z`, :math:`y=z'`) ed applicare il metodo di
Eulero:

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

Pertanto, :math:`z(1)\approx 0.75`.

Si possono anche tracciare i punti :math:`(x,y)` per ottenere un grafico
approssimato della curva. La funzione ``eulers_method_2x2_plot`` svolge
questa funzione; per usarla, bisogna definire le funzioni *f* e
*g* che prendono on argomento con tre coordinate: (*t*, *x*,
*y*).

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

A questo punto, ``P`` ha in memoria due grafici: ``P[0]``, il grafico di *x*
vs. *t*, e ``P[1]``, il grafico di *y* vs. *t*. Si possono tracciare entrambi
come mostrato qui in seguito:

.. link

::

    sage: show(P[0] + P[1])

Funzioni speciali
-----------------

Sono implementati diversi polinomi ortogonali e funzioni
speciali, usando sia PARI [GAP]_ che Maxima [Max]_. Essi
sono documentati nelle sezioni apposite ("Polinomi ortogonali"
e "Funzioni speciali", rispettivamente) del manuale di Sage.

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

A questo punto, Sage ha soltanto incorporato queste funzioni per l'uso numerico.
Per l'uso simbolico, si usi direttamente l'interfaccia di Maxima, come
nell'esempio seguente::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'

.. [GAP] (en) The GAP Group, GAP - Groups, Algorithms, and
   Programming, Version 4.11; 2021, https://www.gap-system.org

.. [Max] (en) Maxima, Version 5.45; 2021, http://maxima.sf.net/
