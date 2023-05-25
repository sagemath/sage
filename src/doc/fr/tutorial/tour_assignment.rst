Affectation, égalité et arithmétique
====================================

A quelques exceptions mineures près, Sage utilise le langage de
programmation Python, si bien que la plupart des ouvrages d'introduction
à Python vous aideront à apprendre Sage.

Sage utilise ``=`` pour les affectations. Il utilise ``==``, ``<=``,
``>=``, ``<`` et ``>`` pour les comparaisons.

::

    sage: a = 5
    sage: a
    5
    sage: 2 == 2
    True
    sage: 2 == 3
    False
    sage: 2 < 3
    True
    sage: a == 5
    True

Sage fournit toutes les opérations mathématiques de base :

::

    sage: 2**3    #  ** désigne l'exponentiation
    8
    sage: 2^3     #  ^ est un synonyme de ** (contrairement à Python)
    8
    sage: 10 % 3  #  pour des arguments entiers, % signifie mod, i.e., le reste dans la division euclidienne
    1
    sage: 10/4
    5/2
    sage: 10//4   #  pour des arguments entiers, // renvoie le quotient dans la division euclidienne
    2
    sage: 4 * (10 // 4) + 10 % 4 == 10
    True
    sage: 3^2*4 + 2%5
    38

Le calcul d'une expression telle que ``3^2*4 + 2%5``  dépend de l'ordre
dans lequel les opérations sont effectuées ; ceci est expliqué dans
l'annexe *Priorité des opérateurs arithmétiques binaires*
:ref:`section-precedence`.

Sage fournit également un grand nombre de fonctions mathématiques
usuelles ; en voici quelques exemples choisis :

::

    sage: sqrt(3.4)
    1.84390889145858
    sage: sin(5.135)
    -0.912021158525540
    sage: sin(pi/3)
    1/2*sqrt(3)

Comme le montre le dernier de ces exemples, certaines expressions
mathématiques renvoient des valeurs 'exactes' plutôt que des
approximations numériques. Pour obtenir une approximation numérique, on
utilise au choix la fonction ``n`` ou la méthode ``n`` (chacun de ces
noms possède le nom plus long ``numerical_approx``, la fonction ``N``
est identique à ``n``). Celles-ci acceptent, en argument optionnel,
``prec``, qui indique le nombre de bits de précisions requis, et
``digits``, qui indique le nombre de décimales demandées ; par défaut,
il y a 53 bits de précision.

::

    sage: exp(2)
    e^2
    sage: n(exp(2))
    7.38905609893065
    sage: sqrt(pi).numerical_approx()
    1.77245385090552
    sage: sin(10).n(digits=5)
    -0.54402
    sage: N(sin(10),digits=10)
    -0.5440211109
    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

Python est doté d'un typage dynamique. Ainsi la valeur à laquelle fait
référence une variable est toujours associée à un type donné, mais une
variable donnée peut contenir des valeurs de plusieurs types Python au
sein d'une même portée :


::

    sage: a = 5   # a est un entier
    sage: type(a)
    <class 'sage.rings.integer.Integer'>
    sage: a = 5/3  # a est maintenant un rationnel...
    sage: type(a)
    <class 'sage.rings.rational.Rational'>
    sage: a = 'hello'  # ...et maintenant une chaîne
    sage: type(a)
    <... 'str'>

Le langage de programmation C, qui est statiquement typé, est bien
différent : une fois déclarée de type int, une variable ne peut contenir
que des int au sein de sa portée.
