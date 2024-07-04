Algèbre de base et calcul infinitésimal
=======================================

Sage peut accomplir divers calculs d'algèbre et d'analyse de base : par
exemple, trouver les solutions d'équations, dériver, intégrer, calculer
des transformées de Laplace. Voir la documentation
`Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_
pour plus d'exemples.

Résolution d'équations
----------------------

La fonction ``solve`` résout des équations. Pour l'utiliser, il convient
de spécifier d'abord les variables. Les arguments de ``solve`` sont
alors une équation (ou un système d'équations) suivie des variables à
résoudre.


::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

On peut résoudre une équation en une variable en fonction des autres :

::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]

On peut également résoudre un système à plusieurs variables :

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]

L'exemple suivant, qui utilise Sage pour la résolution d'un système
d'équations non-linéaires, a été proposé par Jason Grout. D'abord, on
résout le système de façon symbolique :

::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3], [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]

Pour une résolution numérique, on peut utiliser à la place :

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]

(La fonction ``n`` affiche une approximation numérique ; son argument
indique le nombre de bits de précision.)

Dérivation, intégration, etc.
-----------------------------

Sage est capable de dériver et d'intégrer de nombreuses fonctions. Par
exemple, pour dériver :math:`\sin(u)` par rapport à :math:`u`, on
procède comme suit :

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

Pour calculer la dérivée quatrième de  :math:`\sin(x^2)`:

::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)

Pour calculer la dérivée partielle de  :math:`x^2+17y^2` par rapport à
`x` et `y` respectivement :

::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)
    34*y

Passons aux primitives et intégrales. Pour calculer
:math:`\int x\sin(x^2)\, dx` et
:math:`\int_0^1 \frac{x}{x^2+1}\, dx`

::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

Pour calculer la décomposition en éléments simples de
:math:`\frac{1}{x^2-1}`:

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)

.. _section-systems:

Résolution des équations différentielles
----------------------------------------

On peut utiliser Sage pour étudier les équations différentielles
ordinaires. Pour résoudre l'équation :math:`x'+x-1=0` :

::

    sage: t = var('t')    # on définit une variable t
    sage: function('x')(t)   # on déclare x fonction de cette variable
    x(t)
    sage: DE = lambda y: diff(y,t) + y - 1
    sage: desolve(DE(x(t)), [x(t),t])
    (_C + e^t)*e^(-t)

Ceci utilise l'interface de Sage vers Maxima [Max]_, aussi il se peut
que la sortie diffère un peu des sorties habituelles de Sage. Dans notre
cas, le résultat indique que la solution générale à l'équation
différentielle est :math:`x(t) = e^{-t}(e^{t}+C)`.

Il est aussi possible de calculer des transformées de Laplace. La
transformée de Laplace de :math:`t^2e^t -\sin(t)` s'obtient comme suit :

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3

Voici un exemple plus élaboré. L'élongation à partir du point
d'équilibre de ressorts couplés attachés à gauche à un mur

::

    |------\/\/\/\/\---|masse1|----\/\/\/\/\/----|masse2|
            ressort1                ressort2

est modélisée par le système d'équations différentielles d'ordre 2

.. math::
    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0
    m_2 x_2''+ k_2 (x_2-x_1) = 0,


où :math:`m_{i}` est la masse de l'objet *i*, :math:`x_{i}` est
l'élongation à partir du point d'équilibre de la masse  *i*, et
:math:`k_{i}` est la constante de raideur du ressort *i*.

**Exemple :** Utiliser Sage pour résoudre le problème ci-dessus
avec :math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`, :math:`k_{2}=2`,
:math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`, :math:`x_{2}(0)=3`,
:math:`x_{2}'(0)=0`.

Solution : Considérons la transformée de Laplace de la première équation
(avec les notations :math:`x=x_{1}`, :math:`y=x_{2}`):

::

    sage: de1 = maxima("2*diff(x(t),t, 2) + 6*x(t) - 2*y(t)")
    sage: lde1 = de1.laplace("t","s"); lde1.sage()
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

La réponse n'est pas très lisible, mais elle signifie que

.. math:: -2x'(0) + 2s^2\cdot X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0

(où la transformée de Laplace d'une fonction notée par une lettre
minuscule telle que :math:`x(t)` est désignée par la majuscule
correspondante  :math:`X(s)`). Considérons la transformée de Laplace de
la seconde équation :

::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

Ceci signifie

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.


Injectons les conditions initiales pour  :math:`x(0)`, :math:`x'(0)`,
:math:`y(0)` et :math:`y'(0)` et résolvons les deux équations qui en
résultent :

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s]
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
      Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

À présent, prenons la transformée de Laplace inverse pour obtenir la réponse :

::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)

Par conséquent, la solution est

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).


On peut en tracer le graphe paramétrique en utilisant

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),
    ....:     (t, 0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

Les coordonnées individuelles peuvent être tracées en utilisant

::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), (t, 0, 2*pi), rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), (t, 0, 2*pi), rgbcolor=hue(0.6))
    sage: show(p1 + p2)

Les fonctions de tracé de graphes sont décrites dans la section
:ref:`section-plot` de ce tutoriel. On pourra aussi consulter
[NagleEtAl2004]_, §5.5 pour plus d'informations sur les équations
différentielles.

Méthode d'Euler pour les systèmes d'équations différentielles
-------------------------------------------------------------

Dans l'exemple suivant, nous illustrons la méthode d'Euler pour des
équations différentielles ordinaires d'ordre un et deux. Rappelons
d'abord le principe de la méthode pour les équations du premier ordre.
Etant donné un problème donné avec une valeur initiale sous la forme

.. math::
    y'=f(x,y), \quad y(a)=c,


nous cherchons une valeur approchée de la solution au point
:math:`x=b` avec :math:`b>a`.

Rappelons que par définition de la dérivée

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},


où :math:`h>0` est fixé et petit. Ceci, combiné à l'équation
différentielle, donne
:math:`f(x,y(x))\approx
\frac{y(x+h)-y(x)}{h}`. Aussi :math:`y(x+h)` s'écrit:

.. math::   y(x+h) \approx y(x) + h\cdot f(x,y(x)).


Si nous notons :math:`h\cdot f(x,y(x))` le « terme de correction » (faute
d'un terme plus approprié), et si nous appelons :math:`y(x)`
« l'ancienne valeur de `y` » et :math:`y(x+h)` la « nouvelle valeur de
`y` », cette approximation se réécrit

.. math::   y_{nouveau} \approx y_{ancien} + h\cdot f(x,y_{ancien}).


Divisions l'intervalle entre  `a` et `b` en `n` pas, si bien que
:math:`h=\frac{b-a}{n}`. Nous pouvons alors remplir un tableau avec les
informations utilisées dans la méthode.

============== =======================   =====================
:math:`x`      :math:`y`                 :math:`h\cdot f(x,y)`
============== =======================   =====================
:math:`a`      :math:`c`                 :math:`h\cdot f(a,c)`
:math:`a+h`    :math:`c+h\cdot f(a,c)`         ...
:math:`a+2h`   ...
...
:math:`b=a+nh` ???                             ...
============== =======================   =====================


Le but est est de remplir tous les trous du tableau, ligne après ligne,
jusqu'à atteindre le coefficient « ??? », qui est l'approximation de
:math:`y(b)` au sens de la méthode d'Euler.

L'idée est la même pour les systèmes d'équations différentielles.

**Exemple:** Rechercher une approximation numérique de :math:`z(t)` en
:math:`t=1` en utilisant 4 étapes de la méthode d'Euler, où
:math:`z''+tz'+z=0`, :math:`z(0)=1`, :math:`z'(0)=0`.

Il nous faut réduire l'équation différentielle d'ordre 2 à un système de deux équations différentielles d'ordre 1 (en posant :math:`x=z`,
:math:`y=z'`) et appliquer la méthode d'Euler :

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

On en déduit :math:`z(1)\approx 0.65`.

On peut également tracer le graphe des points :math:`(x,y)` pour obtenir
une image approchée de la courbe. La fonction ``eulers_method_2x2_plot``
réalise cela ; pour l'utiliser, il faut définir les fonctions  `f` et
`g` qui prennent un argument à trois coordonnées : (`t`, `x`, `y`).

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

Arrivé à ce point, ``P`` conserve en mémoire deux graphiques : ``P[0]``,
le graphe de  `x` en fonction de `t`, et ``P[1]``, le graphique de `y`
par rapport à `t`. On peut tracer les deux graphiques simultanément par
:

.. link

::

    sage: show(P[0] + P[1])

(Pour plus d'information sur le tracé de graphiques, voir :ref:`section-plot`.)

Fonctions spéciales
-------------------

Plusieurs familles de polynômes orthogonaux et fonctions spéciales sont
implémentées via PARI [GAP]_ et Maxima [Max]_. Ces fonctions sont
documentées dans les sections correspondantes (*Orthogonal polynomials*
et *Special functions*, respectively) du manuel de référence de Sage
(*Sage reference manual*).

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

Pour l'instant, ces fonctions n'ont été adaptées à Sage que pour une
utilisation numérique. Pour faire du calcul formel, il faut utiliser
l'interface Maxima directement, comme le présente l'exemple suivant :

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'
