.. linkall

**********
Interfaces
**********

L'un des aspects essentiels de Sage est qu'il permet d'effectuer des calculs
utilisant des objets issus de nombreux systèmes de calcul formel de
façon unifiée, avec une interface commune et un langage de programmation
sain.

Les méthodes ``console`` et ``interact`` d'une interface avec un
programme externe font des choses tout-à-fait différentes. Prenons
l'exemple de GAP :

#. ``gap.console()`` : Cette commande ouvre la console GAP. Cela
   transfère le contrôle à GAP ; Sage n'est dans ce cas qu'un moyen
   commode de lancer des programmes, un peu comme le shell sous Unix.

#. ``gap.interact()`` : Cette commande permet d'interagir avec une
   instance de GAP en cours d'exécution, et éventuellement « remplie
   d'objets Sage ». Il est possible d'importer des objets Sage dans la
   session GAP (y compris depuis l'interface interactive), etc.

.. index: PARI; GP

GP/PARI
=======

PARI est un programme C compact, mature, fortement optimisé et
spécialisé en théorie des nombres. Il possède deux
interfaces très différentes utilisables depuis Sage :

-  ``gp`` -- l'interpréteur "**G** o **P** ARI", et

-  ``pari`` -- la bibliothèque C PARI.


Ainsi, les deux commandes suivantes font le même calcul de deux façons
différentes. Les deux sorties ont l'air identiques, mais elles ne le
sont pas en réalité, et ce qui ce passe en coulisses est radicalement
différent.

::

    sage: gp('znprimroot(10007)')
    Mod(5, 10007)
    sage: pari('znprimroot(10007)')
    Mod(5, 10007)

Dans le premier exemple, on démarre une instance de l'interpréteur GP et
on lui envoie la chaîne ``'znprimroot(10007)'``. Il l'évalue, et affecte
le résultat à une variable GP (ce qui occupe un espace qui ne sera pas
libéré dans la mémoire du processus fils GP). La valeur de la variable
est ensuite affichée. Dans le second cas, nul programme séparé n'est
démarré, la chaîne ``'znprimroot(10007)'`` est évaluée par une certaine
fonction de la bibliothèque C PARI. Le résultat est stocké sur le tas de
l'interpréteur Python, et la zone de mémoire utilisée est libérée
lorsque son contenu n'est plus utilisé. Les objets renvoyés par ces deux
commandes sont de types différents :

::

    sage: type(gp('znprimroot(10007)'))
    <class 'sage.interfaces.gp.GpElement'>
    sage: type(pari('znprimroot(10007)'))
    <class 'cypari2.gen.Gen'>

Alors, laquelle des interfaces utiliser ? Tout dépend de ce que vous
cherchez à faire. L'interface GP permet de faire absolument tout ce que
vous pourriez faire avec la ligne de commande GP/PARI habituelle,
puisqu'elle fait appel à celle-ci. En particulier, vous pouvez
l'utiliser pour charger et exécuter des scripts PARI compliqués.
L'interface PARI (via la bibliothèque C) est nettement plus restrictive.
Tout d'abord, toutes les méthodes ne sont pas implémentées.
Deuxièmement, beaucoup de code utilisant par exemple l'intégration
numérique ne fonctionne pas via l'interface PARI. Ceci dit, l'interface
PARI est souvent considérablement plus rapide et robuste que l'interface GP.

(Si l'interface GP manque de mémoire pour évaluer une ligne d'entrée
donnée, elle double silencieusement la taille de la pile et réessaie
d'évaluer la ligne. Ainsi votre calcul ne plantera pas même si vous
n'avez pas évalué convenablement l'espace qu'il nécessite. C'est une
caractéristique commode que l'interpréteur GP habituel ne semble pas
fournir. L'interface PARI, quant à elle, déplace immédiatement les
objets créés hors de la pile de PARI, de sorte que celle-ci ne
grossit pas. Cependant, la taille de chaque objet est limitée à 100 Mo,
sous peine que la pile ne déborde à la création de l'objet. Par
ailleurs, cette copie supplémentaire a un léger impact sur les
performances.)

En résumé, Sage fait appel à la bibliothèque C PARI pour fournir des
fonctionnalités similaires à celle de l'interpréteur PARI/GP, mais
depuis le langage Python, et avec un gestionnaire de mémoire plus
perfectionné.

Commençons par créer une liste PARI à partir d'une liste Python.

::

    sage: v = pari([1,2,3,4,5])
    sage: v
    [1, 2, 3, 4, 5]
    sage: type(v)
    <class 'cypari2.gen.Gen'>

En Sage, les objets PARI sont de type ``Gen``. Le type PARI de
l'objet sous-jacent est donné par la méthode ``type``.

::

    sage: v.type()
    't_VEC'

Pour créer une courbe elliptique en PARI, on utiliserait
``ellinit([1,2,3,4,5])``. La syntaxe Sage est semblable, à ceci près que
``ellinit`` devient une méthode qui peut être appelée sur n'importe quel
objet PARI, par exemle notre ``t_VEC`` `v`.

::

    sage: e = v.ellinit()
    sage: e.type()
    't_VEC'
    sage: pari(e)[:13]
    [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351]

À présent que nous disposons d'une courbe elliptique, faisons quelques
calculs avec.

::

    sage: e.elltors()
    [1, [], []]
    sage: e.ellglobalred()
    [10351, [1, -1, 0, -1], 1, [11, 1; 941, 1], [[1, 5, 0, 1], [1, 5, 0, 1]]]
    sage: f = e.ellchangecurve([1,-1,0,-1])
    sage: f[:5]
    [1, -1, 0, 4, 3]

.. index: GAP

.. _section-gap:

GAP
===

Pour les mathématiques discrètes effectives et principalement la théorie
des groupes, Sage utilise GAP.

Voici un exemple d'utilisation de la fonction GAP ``IdGroup``, qui
nécessite une base de données optionnelle de groupes de petit ordre, à
installer séparément comme décrit plus bas.

::

    sage: G = gap('Group((1,2,3)(4,5), (3,4))')
    sage: G
    Group( [ (1,2,3)(4,5), (3,4) ] )
    sage: G.Center()
    Group( () )
    sage: G.IdGroup()
    [ 120, 34 ]
    sage: G.Order()
    120

On peut faire le même calcul en SAGE sans invoquer explicitement
l'interface GAP comme suit :

::

    sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
    sage: G.center()
    Subgroup generated by [()] of (Permutation Group with generators [(3,4), (1,2,3)(4,5)])
    sage: G.group_id()
    [120, 34]
    sage: n = G.order(); n
    120

Pour utiliser certaines fonctionnalités de GAP,
vous devez installer un paquet Sage optionnel.
Cela peut être fait avec la commande::

    sage -i gap_packages


Singular
========

Singular fournit une bibliothèque consistante et mature qui permet, entre
autres, de calculer des pgcd de polynômes de plusieurs variables, des
factorisations, des bases de Gröbner ou encore des bases d'espaces de
Riemann-Roch de courbes planes. Considérons la factorisation de
polynômes de plusieurs variables à l'aide de l'interface à Singular
fournie par Sage (n'entrez pas les ``....:``) :

::

    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: R1
    polynomial ring, over a field, global ordering
    //   coefficients: QQ
    //   number of vars : 2
    //        block   1 : ordering dp
    //                  : names    x y
    //        block   2 : ordering C
    sage: f = singular('9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 +'
    ....:     '9*x^6*y^4 + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 -'
    ....:     '9*x^12*y^3 - 18*x^13*y^2 + 9*x^16')

Maintenant que nous avons défini :math:`f`, affichons-le puis
factorisons-le.

::

    sage: f
    9*x^16-18*x^13*y^2-9*x^12*y^3+9*x^10*y^4-18*x^11*y^2+36*x^8*y^4+18*x^7*y^5-18*x^5*y^6+9*x^6*y^4-18*x^3*y^6-9*x^2*y^7+9*y^8
    sage: f.parent()
    Singular
    sage: F = f.factorize(); F
    [1]:
       _[1]=9
       _[2]=x^6-2*x^3*y^2-x^2*y^3+y^4
       _[3]=-x^5+y^2
    [2]:
       1,1,2
    sage: F[1][2]
    x^6-2*x^3*y^2-x^2*y^3+y^4

Comme avec GAP dans la section :ref:`section-gap`, nous pouvons aussi
calculer la factorisation sans utiliser explicitement l'interface
Singular (Sage y fera tout de même appel en coulisses pour le calcul).

::

    sage: x, y = QQ['x, y'].gens()
    sage: f = (9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 + 9*x^6*y^4
    ....:     + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 - 9*x^12*y^3
    ....:     - 18*x^13*y^2 + 9*x^16)
    sage: factor(f)
    (9) * (-x^5 + y^2)^2 * (x^6 - 2*x^3*y^2 - x^2*y^3 + y^4)

.. _section-maxima:

Maxima
======

Le système de calcul formel Maxima est fourni avec Sage accompagné d'une
implémentation du langage Lisp. Le logiciel gnuplot (que Maxima utilise
par défaut pour tracer des graphiques) est disponible comme paquet
optionnel. Maxima fournit notamment des routines de calcul sur des
expressions formelles. Il permet de calculer des dérivées, primitives et
intégrales, de résoudre des équations différentielles d'ordre 1 et
souvent d'ordre 2, et de résoudre par transformée de Laplace les
équations différentielles linéaires d'ordre quelconque. Maxima dispose
aussi d'un grand nombre de fonctions spéciales, permet de tracer des
graphes de fonctions via gnuplot, et de manipuler des matrices
(réduction en lignes, valeurs propres, vecteurs propres...) ou encore
des équations polynomiales.

Utilisons par exemple l'interface Sage/Maxima pour construire
la matrice dont le coefficient d'indice :math:`i,j` vaut :math:`i/j`,
pour :math:`i,j=1,\ldots,4`.

::

    sage: f = maxima.eval('ij_entry[i,j] := i/j')
    sage: A = maxima('genmatrix(ij_entry,4,4)'); A
    matrix([1,1/2,1/3,1/4],[2,1,2/3,1/2],[3,3/2,1,3/4],[4,2,4/3,1])
    sage: A.determinant()
    0
    sage: A.echelon()
    matrix([1,1/2,1/3,1/4],[0,0,0,0],[0,0,0,0],[0,0,0,0])
    sage: A.eigenvalues()
    [[0,4],[3,1]]
    sage: A.eigenvectors().sage()
    [[[0, 4], [3, 1]], [[[1, 0, 0, -4], [0, 1, 0, -2], [0, 0, 1, -4/3]], [[1, 2, 3, 4]]]]

Un deuxième exemple :

::

    sage: A = maxima("matrix ([1, 0, 0], [1, -1, 0], [1, 3, -2])")
    sage: eigA = A.eigenvectors()
    sage: V = VectorSpace(QQ,3)
    sage: eigA
    [[[-2,-1,1],[1,1,1]],[[[0,0,1]],[[0,1,3]],[[1,1/2,5/6]]]]
    sage: v1 = V(sage_eval(repr(eigA[1][0][0]))); lambda1 = eigA[0][0][0]
    sage: v2 = V(sage_eval(repr(eigA[1][1][0]))); lambda2 = eigA[0][0][1]
    sage: v3 = V(sage_eval(repr(eigA[1][2][0]))); lambda3 = eigA[0][0][2]

    sage: M = MatrixSpace(QQ,3,3)
    sage: AA = M([[1,0,0],[1, - 1,0],[1,3, - 2]])
    sage: b1 = v1.base_ring()
    sage: AA*v1 == b1(lambda1)*v1
    True
    sage: b2 = v2.base_ring()
    sage: AA*v2 == b2(lambda2)*v2
    True
    sage: b3 = v3.base_ring()
    sage: AA*v3 == b3(lambda3)*v3
    True

Voici enfin quelques exemples de tracés de graphiques avec ``openmath``
depuis Sage. Un grand nombre de ces exemples sont des adaptations de
ceux du manuel de référence de Maxima.

Tracé en 2D de plusieurs fonctions (n'entrez pas les ``....:``) :

::

    sage: maxima.plot2d('[cos(7*x),cos(23*x)^4,sin(13*x)^3]','[x,0,1]', # not tested
    ....:     '[plot_format,openmath]')

Un graphique 3D interactif, que vous pouvez déplacer à la souris
(n'entrez pas les ``....:``) :

::

    sage: maxima.plot3d ("2^(-u^2 + v^2)", "[u, -3, 3]", "[v, -2, 2]", # not tested
    ....:     '[plot_format, openmath]')
    sage: maxima.plot3d("atan(-x^2 + y^3/4)", "[x, -4, 4]", "[y, -4, 4]", # not tested
    ....:     "[grid, 50, 50]",'[plot_format, openmath]')

Le célèbre ruban de Möbius (n'entrez pas les ``....:``) :

::

    sage: maxima.plot3d("[cos(x)*(3 + y*cos(x/2)), sin(x)*(3 + y*cos(x/2)), y*sin(x/2)]", # not tested
    ....:     "[x, -4, 4]", "[y, -4, 4]",
    ....:     '[plot_format, openmath]')

Et la fameuse bouteille de Klein (n'entrez pas les ``....:``):

::

    sage: _ = maxima("expr_1: 5*cos(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0) - 10.0")
    sage: _ = maxima("expr_2: -5*sin(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0)")
    sage: _ = maxima("expr_3: 5*(-sin(x/2)*cos(y) + cos(x/2)*sin(2*y))")
    sage: maxima.plot3d ("[expr_1, expr_2, expr_3]", "[x, -%pi, %pi]", # not tested
    ....:     "[y, -%pi, %pi]", "['grid, 40, 40]",
    ....:     '[plot_format, openmath]')

