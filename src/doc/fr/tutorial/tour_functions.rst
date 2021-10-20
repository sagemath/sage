.. _section-functions-issues:

Problèmes fréquents concernant les fonctions
============================================

La définition de fonctions, par exemple pour calculer leurs dérivées ou
tracer leurs courbes représentatives, donne lieu à un certain nombre de
confusions. Le but de cette section est de clarifier quelques points à
l'origine de ces confusions.

Il y a plusieurs façons de définir un objet que l'on peut légitimement
appeler « fonction ».

1. Définir une fonction Python, comme expliqué dans la section :ref:`section-functions`. Les fonctions Python peuvent être utilisées
pour tracer des courbes, mais pas dérivées ou intégrées symboliquement::

    sage: def f(z): return z^2
    sage: type(f)
    <... 'function'>
    sage: f(3)
    9
    sage: plot(f, 0, 2)
    Graphics object consisting of 1 graphics primitive

Remarquez la syntaxe de la dernière ligne. Écrire plutôt ``plot(f(z), 0, 2)``
provoquerait une erreur : en effet, le ``z`` qui apparaît dans
la définition de ``f`` est une variable muette qui n'a pas de sens
en dehors de la définition. Un simple ``f(z)`` déclenche la même erreur.
En l'occurrence, faire de ``z`` une variable symbolique comme dans
l'exemple ci-dessous fonctionne, mais cette façon de faire soulève
d'autres problèmes (voir le point 4 ci-dessous), et il vaut mieux
s'abstenir de l'utiliser.

.. link

::

    sage: var('z')   # on définit z comme variable symbolique
    z
    sage: f(z)
    z^2
    sage: plot(f(z), 0, 2)
    Graphics object consisting of 1 graphics primitive

L'appel de fonction ``f(z)`` renvoie ici l'expression symbolique
``z^2``, qui est alors utilisée par la fonction ``plot``.

2. Définir une expression symbolique fonctionnelle (« appelable »). Une
telle expression représente une fonction dont on peut tracer le graphe,
et que l'on peut aussi dériver ou intégrer symboliquement ::

    sage: g(x) = x^2
    sage: g        # g envoie x sur x^2
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

Notez que, si ``g`` est une expression symbolique fonctionnelle
(``x |--> x^2``), l'objet ``g(x)`` (``x^2``) est d'une nature un
peu différente. Les expressions comme ``g(x)`` peuvent aussi être
tracées, dérivées, intégrées, etc., avec cependant quelques difficultés
illustrées dans le point 5 ci-dessous.

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

3. Utiliser une fonction usuelle prédéfinie de Sage. Celles-ci peuvent
servir à tracer des courbes, et, indirectement, être dérivées ou intégrées ::

    sage: type(sin)
    <class 'sage.functions.trig.Function_sin'>
    sage: plot(sin, 0, 2)
    Graphics object consisting of 1 graphics primitive
    sage: type(sin(x))
    <class 'sage.symbolic.expression.Expression'>
    sage: plot(sin(x), 0, 2)
    Graphics object consisting of 1 graphics primitive

Il n'est pas possible de dériver la fonction ``sin`` tout court pour
obtenir ``cos`` ::

    sage: f = sin
    sage: f.derivative()
    Traceback (most recent call last):
    ...
    AttributeError: ...

Une possibilité est de remplacer ``f = sin`` par ``f = sin(x)``, mais il
est généralement préférable de définir une expression symbolique
fonctionnelle ``f(x) = sin(x)`` ::

    sage: S(x) = sin(x)
    sage: S.derivative()
    x |--> cos(x)

Examinons maintenant quelques problèmes fréquents.

\4. Évaluation accidentelle ::

    sage: def h(x):
    ....:     if x < 2:
    ....:         return 0
    ....:     else:
    ....:         return x-2

Problème : ``plot(h(x), 0, 4)`` trace la droite `y = x - 2`, et non pas la
fonction affine par morceaux définie par ``h``. Pourquoi ? Lors de l'exécution,
``plot(h(x), 0, 4)`` évalue d'abord ``h(x)`` : la fonction
Python ``h`` est appelée avec le paramètre ``x``, et la condition ``x < 2``
est donc évaluée.

.. link

::

    sage: type(x < 2)
    <class 'sage.symbolic.expression.Expression'>

Or, l'évaluation d'une inégalité symbolique renvoie ``False`` quand la
condition n'est pas clairement vraie. Ainsi, ``h(x)`` s'évalue en
``x - 2``, et c'est cette expression-là qui est finalement tracée.

Solution : Il ne faut pas utiliser ``plot(h(x), 0, 4)``, mais plutôt

.. link

::

    sage: def h(x):
    ....:     if x < 2:
    ....:         return 0
    ....:     else:
    ....:         return x-2
    sage: plot(h, 0, 4)
    Graphics object consisting of 1 graphics primitive

\5. Constante plutôt que fonction ::

    sage: f = x
    sage: g = f.derivative()
    sage: g
    1

Problème : ``g(3)`` déclenche une erreur avec le message « ValueError:
the number of arguments must be less than or equal to 0 ».

.. link

::

    sage: type(f)
    <class 'sage.symbolic.expression.Expression'>
    sage: type(g)
    <class 'sage.symbolic.expression.Expression'>

En effet, ``g`` n'est pas une fonction, mais une constante, sans
variable en laquelle on peut l'évaluer.

Solution : il y a plusieurs possibilités.

- Définir ``f`` comme une expression symbolique fonctionnelle ::

    sage: f(x) = x        # au lieu de 'f = x'
    sage: g = f.derivative()
    sage: g
    x |--> 1
    sage: g(3)
    1
    sage: type(g)
    <class 'sage.symbolic.expression.Expression'>

- Ou, sans changer la définition de ``f``, définir ``g`` comme une
  expression symbolique fonctionnelle ::

    sage: f = x
    sage: g(x) = f.derivative()  # au lieu de 'g = f.derivative()'
    sage: g
    x |--> 1
    sage: g(3)
    1
    sage: type(g)
    <class 'sage.symbolic.expression.Expression'>

- Ou encore, avec ``f`` et ``g`` définies comme dans l'exemple de
  départ, donner explicitement la variable à remplacer par sa valeur ::

    sage: f = x
    sage: g = f.derivative()
    sage: g
    1
    sage: g(x=3)    # au lieu de  'g(3)'
    1

Nous terminons en mettant encore une fois en évidence la différence entre
les dérivées des expressions ``f`` définies par ``f = x`` et par ``f(x)
= x`` ::

    sage: f(x) = x
    sage: g = f.derivative()
    sage: g.variables()  # variables apparaissant dans g
    ()
    sage: g.arguments()  # paramètres auxquels on peut donner une valeur dans g
    (x,)
    sage: f = x
    sage: h = f.derivative()
    sage: h.variables()
    ()
    sage: h.arguments()
    ()

Comme l'illustre cet exemple, ``h`` n'accepte pas de paramètres. C'est
pour cela que ``h(3)`` déclenche une erreur.
