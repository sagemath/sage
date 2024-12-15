*************
Programmation
*************

.. _section-loadattach:

Charger et attacher des fichiers Sage
=====================================

Nous décrivons maintenant la manière de charger dans Sage des programmes
écrits dans des fichiers séparés. Créons un fichier appelé
``example.sage`` avec le contenu suivant :

.. skip

::

    print("Hello World")
    print(2^3)

Nous pouvons lire et exécuter le contenu du fichier ``example.sage``
en utilisant la commande ``load``.

.. skip

::

    sage: load("example.sage")
    Hello World
    8

Nous pouvons aussi attacher un fichier Sage à la session en cours avec
la commande ``attach`` :

.. skip

::

    sage: attach("example.sage")
    Hello World
    8

L'effet de cette commande est que si nous modifions maintenant
``example.sage`` et entrons une ligne vierge (i.e. appuyons sur
``entrée``), le contenu de ``example.sage`` sera automatiquement
rechargé dans Sage.

Avec ``attach``, le fichier est rechargé automatiquement
dans Sage à chaque modification, ce qui est pratique pour déboguer du
code ; tandis qu'avec ``load`` il n'est chargé qu'une fois.

Lorsque Sage lit ``example.sage``, il le convertit en un programme
Python qui est ensuite exécuté par l'interpréteur Python. Cette
conversion est minimale, elle se résume essentiellement à encapsuler les
littéraux entiers dans des ``Integer()`` et les littéraux flottants dans
des ``RealNumber()``, à remplacer les ``^`` par des ``**``, et à
remplacer par exemple ``R.2`` par ``R.gen(2)``. La version convertie en
Python de ``example.sage`` est placée dans le même répertoire sous le
nom ``example.sage.py``. Elle contient le code suivant :

::

    print("Hello World")
    print(Integer(2)**Integer(3))

On voit que les littéraux entiers ont été encapsulés et que le ``^`` a
été remplacé par ``**`` (en effet, en Python, ``^`` représente le ou
exclusif et ``**`` l'exponentiation).

(Ce prétraitement est implémenté dans le fichier
``sage/misc/interpreter.py``.)

On peut coller dans Sage des blocs de code de plusieurs lignes avec des
indentations pourvu que les blocs soient délimités par des
retours à la ligne (cela n'est pas nécessaire quand le code est dans
un fichier). Cependant, la meilleure façon d'entrer ce genre de code
dans Sage est de l'enregistrer dans un fichier et d'utiliser la commande
``attach`` comme décrit ci-dessus.


.. _section-compile:

Écrire des programmes compilés
==============================

Dans les calculs mathématiques sur ordinateur, la vitesse a une
importance cruciale. Or, si Python est un langage commode et de très haut
niveau, certaines opérations peuvent être plus rapides de plusieurs
ordres de grandeur si elles sont implémentées sur des types statiques
dans un langage compilé. Certaines parties de Sage auraient été trop
lentes si elles avaient été écrites entièrement en Python. Pour pallier
ce problème, Sage supporte une sorte de « version compilée » de Python
appelée Cython (voir [Cyt]_ et [Pyr]_). Cython ressemble à la fois
à Python et à C. La plupart des constructions Python, dont la définition
de listes par compréhension, les expressions conditionnelles, les
constructions comme ``+=`` sont autorisées en Cython. Vous pouvez aussi
importer du code depuis d'autres modules Python. En plus de cela,
vous pouvez déclarer des variables C arbitraires, et faire directement
des appels arbitraires à des bibliothèques C. Le code Cython est
converti en C et compilé avec le compilateur C.

Pour créer votre propre code Sage compilé, donnez à votre fichier source
l'exension ``.spyx`` (à la place de ``.sage``). Avec l'interface en
ligne de commande, vous pouvez charger ou attacher des fichiers de code
compilé exactement comme les fichiers interprétés. Pour l'instant,
le *notebook* ne permet pas d'attacher des fichiers compilés. La
compilation proprement dite a lieu « en coulisse », sans que vous ayez à
la déclencher explicitement.  La bibliothèque d'objets partagés compilés se trouve
dans ``$HOME/.sage/temp/hostname/pid/spyx``. Ces fichiers sont supprimés
lorsque vous quittez Sage.

Attention, le prétraitement des fichiers Sage mentionné plus haut N'EST
PAS appliqué aux fichiers spyx, ainsi, dans un fichier spyx, ``1/3``
vaut 0 et non le nombre rationnel :math:`1/3`. Pour appeler une fonction
``foo`` de la bibliothèque Sage depuis un fichier spyx, importez
``sage.all`` et appelez ``sage.all.foo``.

::

    import sage.all
    def foo(n):
        return sage.all.factorial(n)

Appeler des fonctions définies dans des fichiers C séparés
----------------------------------------------------------

Il n'est pas difficile non plus d'accéder à des fonctions écrites en C,
dans des fichiers \*.c séparés. Créez dans un même répertoire deux
fichiers ``test.c`` et ``test.spyx`` avec les contenus suivants :

Le code C pur : ``test.c``

::

    int add_one(int n) {
      return n + 1;
    }

Le code Cython : ``test.spyx``:

::

    cdef extern from "test.c":
        int add_one(int n)

    def test(n):
        return add_one(n)

Vous pouvez alors faire :

.. skip

::

    sage: attach("test.spyx")
    Compiling (...)/test.spyx...
    sage: test(10)
    11

Si la compilation du code C généré à partir d'un fichier Cython
nécessite une bibliothèque supplémentaire ``foo``, ajoutez au source
Cython la ligne ``clib foo``. De même, il est possible d'ajouter un
fichier C supplémentaire ``bar`` aux fichiers à compiler avec la
déclaration ``cfile bar``.

.. _section-standalone:

Scripts Python/Sage autonomes
=============================

Le script autonome suivant, écrit en Sage, permet de factoriser des
entiers, des polynômes, etc. :

::

    #!/usr/bin/env sage

    import sys

    if len(sys.argv) != 2:
        print("Usage: %s <n>" % sys.argv[0])
        print("Outputs the prime factorization of n.")
        sys.exit(1)

    print(factor(sage_eval(sys.argv[1])))

Pour utiliser ce script, votre répertoire ``SAGE_ROOT`` doit apparaître
dans la variable d'environnement PATH. Supposons que le script ci-dessus
soit appelé ``factor``, il peut alors être utilisé comme dans l'exemple
suivant :

::

    bash $ ./factor 2006
    2 * 17 * 59

Types de données
================

Chaque objet Sage a un type bien défini. Python dispose d'une vaste
gamme de types intégrés et la bibliothèque Sage en fournit de nombreux
autres. Parmi les types intégrés de Python, citons les chaînes, les
listes, les n-uplets, les entiers et les flottants :

::

    sage: s = "sage"; type(s)
    <... 'str'>
    sage: s = 'sage'; type(s)      # guillemets simples ou doubles
    <... 'str'>
    sage: s = [1,2,3,4]; type(s)
    <... 'list'>
    sage: s = (1,2,3,4); type(s)
    <... 'tuple'>
    sage: s = int(2006); type(s)
    <... 'int'>
    sage: s = float(2006); type(s)
    <... 'float'>

Sage ajoute de nombreux autres types. Par exemple, les espaces
vectoriels :

::

    sage: V = VectorSpace(QQ, 1000000); V
    Vector space of dimension 1000000 over Rational Field
    sage: type(V)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category'>

Seules certaines fonctions peuvent être appelées sur ``V``. Dans
d'autres logiciels mathématiques, cela se fait en notation
« fonctionnelle », en écrivant ``foo(V,...)``. En Sage, certaines
fonctions sont attachés au type (ou classe) de l'objet et appelées avec
une syntaxe « orientée objet » comme en Java ou en C++, par exemple
``V.foo(...)``. Cela évite de polluer l'espace de noms global avec des
dizaines de milliers de fonctions, et cela permet d'avoir plusieurs
fonctions appelées ``foo``, avec des comportements différents, sans devoir
se reposer sur le type des arguments (ni sur des instructions case) pour
décider laquelle appeler. De plus, une fonction dont vous réutilisez le
nom demeure disponible : par exemple, si vous appelez quelque chose
``zeta`` et si ensuite vous voulez calculer la valeur de la fonction
zêta de Riemann au point 0.5, vous pouvez encore écrire ``s=.5;
s.zeta()``.

::

    sage: zeta = -1
    sage: s=.5; s.zeta()
    -1.46035450880959

La notation fonctionnelle usuelle est aussi acceptée dans certains cas
courants, par commodité et parce que certaines expressions mathématiques
ne sont pas claires en notation orientée objet. Voici quelques exemples.

::

    sage: n = 2; n.sqrt()
    sqrt(2)
    sage: sqrt(2)
    sqrt(2)
    sage: V = VectorSpace(QQ,2)
    sage: V.basis()
        [(1, 0), (0, 1)]
    sage: basis(V)
        [(1, 0), (0, 1)]
    sage: M = MatrixSpace(GF(7), 2); M
    Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
    sage: A = M([1,2,3,4]); A
    [1 2]
    [3 4]
    sage: A.charpoly('x')
    x^2 + 2*x + 5
    sage: charpoly(A, 'x')
    x^2 + 2*x + 5

Pour obtenir la liste de toutes les fonctions membres de :math:`A`,
utilisez la complétion de ligne de commande : tapez ``A.``, puis appuyez
sur la touche ``[tab]`` de votre clavier, comme expliqué dans la section
:ref:`section-tabcompletion`.


Listes, n-uplets et séquences
=============================

Une liste stocke des éléments qui peuvent être de type arbitraire. Comme
en C, en C++ etc. (mais au contraire de ce qu'il se passe dans la
plupart des systèmes de calcul formel usuels) les éléments de la liste
sont indexés à partir de :math:`0` :

::

    sage: v = [2, 3, 5, 'x', SymmetricGroup(3)]; v
    [2, 3, 5, 'x', Symmetric group of order 3! as a permutation group]
    sage: type(v)
    <... 'list'>
    sage: v[0]
    2
    sage: v[2]
    5

Lors d'un accès à une liste, l'index n'a pas besoin d'être un entier
Python. Un entier (Integer) Sage (ou un Rational, ou n'importe quoi
d'autre qui a une méthode ``__index__``) fait aussi l'affaire.

::

    sage: v = [1,2,3]
    sage: v[2]
    3
    sage: n = 2      # Integer (entier Sage)
    sage: v[n]       # ça marche !
    3
    sage: v[int(n)]  # Ok aussi
    3

La fonction ``range`` crée une liste d'entiers Python (et non d'entiers
Sage) :

::

    sage: list(range(1, 15))
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

Cela est utile pour construire des listes par compréhension :

::

    sage: L = [factor(n) for n in range(1, 15)]
    sage: L
    [1, 2, 3, 2^2, 5, 2 * 3, 7, 2^3, 3^2, 2 * 5, 11, 2^2 * 3, 13, 2 * 7]
    sage: L[12]
    13
    sage: type(L[12])
    <class 'sage.structure.factorization_integer.IntegerFactorization'>
    sage: [factor(n) for n in range(1, 15) if is_odd(n)]
    [1, 3, 5, 7, 3^2, 11, 13]

Pour plus d'information sur les compréhensions, voir [PyT]_.

Une fonctionnalité merveilleuse est l'extraction de tranches d'une
liste. Si ``L`` est une liste, ``L[m:n]`` renvoie la sous-liste de ``L``
formée des éléments d'indices :math:`m` à :math:`n-1` inclus :

::

    sage: L = [factor(n) for n in range(1, 20)]
    sage: L[4:9]
    [5, 2 * 3, 7, 2^3, 3^2]
    sage: L[:4]
    [1, 2, 3, 2^2]
    sage: L[14:4]
    []
    sage: L[14:]
    [3 * 5, 2^4, 17, 2 * 3^2, 19]

Les n-uplets ressemblent aux listes, à ceci près qu'ils sont non
mutables, ce qui signifie qu'ils ne peuvent plus être modifiés
une fois créés.

::

    sage: v = (1,2,3,4); v
    (1, 2, 3, 4)
    sage: type(v)
    <... 'tuple'>
    sage: v[1] = 5
    Traceback (most recent call last):
    ...
    TypeError: 'tuple' object does not support item assignment

Les séquences sont un troisième type Sage analogue aux listes.
Contrairement aux listes et aux n-uplets, il ne s'agit pas d'un type
interne de Python. Par défaut, les séquences sont mutables, mais on
peut interdire leur modification en utilisant la méthode
``set_immutable`` de la classe ``Sequence``, comme dans l'exemple
suivant. Tous les éléments d'une séquence ont un parent commun, appelé
l'univers de la séquence.

::

    sage: v = Sequence([1,2,3,4/5])
    sage: v
    [1, 2, 3, 4/5]
    sage: type(v)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: type(v[1])
    <class 'sage.rings.rational.Rational'>
    sage: v.universe()
    Rational Field
    sage: v.is_immutable()
    False
    sage: v.set_immutable()
    sage: v[0] = 3
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Les séquences sont des objets dérivés des listes, et peuvent être
utilisées partout où les listes peuvent l'être :

::

    sage: v = Sequence([1,2,3,4/5])
    sage: isinstance(v, list)
    True
    sage: list(v)
    [1, 2, 3, 4/5]
    sage: type(list(v))
    <... 'list'>

Autre exemple : les bases d'espaces vectoriels sont des séquences non
mutables, car il ne faut pas les modifier.

::

    sage: V = QQ^3; B = V.basis(); B
    [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    sage: type(B)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: B[0] = B[1]
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.
    sage: B.universe()
    Vector space of dimension 3 over Rational Field

Dictionnaires
=============

Un dictionnaire (parfois appelé un tableau associatif) est une
correspondance entre des objets « hachables » (par exemple des chaînes, des
nombres, ou des n-uplets de tels objets, voir
http://docs.python.org/tut/node7.html et
http://docs.python.org/lib/typesmapping.html dans la documentation de
Python pour plus de détails) vers des objets arbitraires.

::

    sage: d = {1:5, 'sage':17, ZZ:GF(7)}
    sage: type(d)
    <... 'dict'>
    sage: list(d.keys())
    [1, 'sage', Integer Ring]
    sage: d['sage']
    17
    sage: d[ZZ]
    Finite Field of size 7
    sage: d[1]
    5

La troisième clé utilisée ci-dessus, l'anneau des entiers relatifs,
montre que les indices d'un dictionnaire peuvent être des objets
compliqués.

Un dictionnaire peut être transformé en une liste de couples clé-objet
contenant les mêmes données :

.. link

::

    sage: list(d.items())
    [(1, 5), ('sage', 17), (Integer Ring, Finite Field of size 7)]

Le parcours itératifs des paires d'un dictionnaire est un idiome de
programmation fréquent :

::

    sage: d = {2:4, 3:9, 4:16}
    sage: [a*b for a, b in d.items()]
    [8, 27, 64]

Comme le montre la dernière sortie ci-dessus, un dictionnaire stocke ses
éléments sans ordre particulier.

Ensembles
=========

Python dispose d'un type ensemble intégré. Sa principale caractéristique
est qu'il est possible de tester très rapidement si un élément
appartient ou non à un ensemble. Le type ensemble fournit les opérations
ensemblistes usuelles.

::

    sage: X = set([1,19,'a']);   Y = set([1,1,1, 2/3])
    sage: X   # random sort order
    {1, 19, 'a'}
    sage: X == set(['a', 1, 1, 19])
    True
    sage: Y
    {2/3, 1}
    sage: 'a' in X
    True
    sage: 'a' in Y
    False
    sage: X.intersection(Y)
    {1}

Sage a son propre type ensemble, qui est (dans certains cas) implémenté
au-dessus du type Python, mais offre quelques fonctionnalités
supplémentaires utiles à Sage. Pour créer un ensemble Sage, on utilise
``Set(...)``. Par exemple,

::

    sage: X = Set([1,19,'a']);   Y = Set([1,1,1, 2/3])
    sage: X   # random sort order
    {'a', 1, 19}
    sage: X == Set(['a', 1, 1, 19])
    True
    sage: Y
    {1, 2/3}
    sage: X.intersection(Y)
    {1}
    sage: print(latex(Y))
    \left\{1, \frac{2}{3}\right\}
    sage: Set(ZZ)
    Set of elements of Integer Ring

Itérateurs
==========

Les itérateurs sont un ajout récent à Python, particulièrement utile
dans les applications mathématiques. Voici quelques exemples, consultez
[PyT]_ pour plus de détails. Fabriquons un itérateur sur les carrés
d'entiers positifs jusqu'à :math:`10000000`.

::

    sage: v = (n^2 for n in range(10000000))
    sage: next(v)
    0
    sage: next(v)
    1
    sage: next(v)
    4

Nous créons maintenant un itérateur sur les nombres premiers de la forme
:math:`4p+1` où :math:`p` est lui aussi premier, et nous examinons les
quelques premières valeurs qu'il prend.

::

    sage: w = (4*p + 1 for p in Primes() if is_prime(4*p+1))
    sage: w
    <generator object <genexpr> at 0x...>
    sage: next(w)
    13
    sage: next(w)
    29
    sage: next(w)
    53

Certains anneaux, par exemple les corps finis et les entiers, disposent
d'itérateurs associés :

::

    sage: [x for x in GF(7)]
    [0, 1, 2, 3, 4, 5, 6]
    sage: W = ((x,y) for x in ZZ for y in ZZ)
    sage: next(W)
    (0, 0)
    sage: next(W)
    (0, 1)
    sage: next(W)
    (0, -1)

Boucles, fonctions, structures de contrôle et comparaisons
==========================================================

Nous avons déjà vu quelques exemples courants d'utilisation des boucles
``for``. En Python, les boucles ``for`` ont la structure suivante, avec
une indentation :

::

    >>> for i in range(5):
    ...     print(i)
    ...
    0
    1
    2
    3
    4

Notez bien les deux points à la fin de l'instruction for (il n'y a pas de
« do » ou « od » comme en Maple ou en GAP) ainsi que l'indentation du
corps de la boucle, formé de l'unique instruction ``print(i)``. Cette
indentation est significative, c'est elle qui délimite le corps de la
boucle. Depuis la ligne de commande Sage, les lignes suivantes sont
automatiquement indentées quand vous appuyez sur ``entrée`` après un
signe « : », comme illustré ci-dessous.

::

    sage: for i in range(5):
    ....:     print(i)  # appuyez deux fois sur entrée ici
    0
    1
    2
    3
    4


Le signe ``=`` représente l'affectation.
L'opérateur ``==`` est le test d'égalité.

::

    sage: for i in range(15):
    ....:     if gcd(i,15) == 1:
    ....:         print(i)
    1
    2
    4
    7
    8
    11
    13
    14

Retenez bien que l'indentation détermine la structure en blocs des
instructions ``if``, ``for`` et ``while`` :

::

    sage: def legendre(a,p):
    ....:     is_sqr_modp=-1
    ....:     for i in range(p):
    ....:         if a % p == i^2 % p:
    ....:             is_sqr_modp=1
    ....:     return is_sqr_modp

    sage: legendre(2,7)
    1
    sage: legendre(3,7)
    -1

Naturellement, l'exemple précédent n'est pas une implémentation efficace du
symbole de Legendre ! Il est simplement destiné à illustrer différents
aspects de la programmation Python/Sage. La fonction {kronecker} fournie
avec Sage calcule le symbole de Legendre efficacement, en appelant la
bibliothèque C de PARI.

Remarquons aussi que les opérateurs de comparaison numériques comme ``==``,
``!=``, ``<=``, ``>=``, ``>``, ``<`` convertissent automatiquement leurs
deux membres en des nombres du même type lorsque c'est possible :

::

    sage: 2 < 3.1; 3.1 <= 1
    True
    False
    sage: 2/3 < 3/2;   3/2 < 3/1
    True
    True

Pour évaluer des inégalités symboliques, utilisez ``bool`` :

::

    sage: x < x + 1
    x < x + 1
    sage: bool(x < x + 1)
    True

Lorsque l'on cherche à comparer des objets de types différents, Sage essaie le
plus souvent de trouver une coercition canonique des deux objets dans un même
parent (voir la section :ref:`section-coercion` pour plus de détails). Si cela
réussit, la comparaison est faite entre les objets convertis ; sinon, les
objets sont simplement considérés comme différents. Pour tester si deux
variables font référence au même objet, on utilise l'opérateur ``is``. Ainsi,
l'entier Python (``int``) ``1`` est unique, mais pas l'entier Sage ``1`` :

::

    sage: 1 is 2/2
    False
    sage: 1 is 1
    False
    sage: 1 == 2/2
    True

Dans les deux lignes suivantes, la première égalité est fausse parce
qu'il n'y a pas de morphisme canonique :math:`\QQ\to
\GF{5}`, et donc pas de manière canonique de comparer l'élément
:math:`1` de :math:`\GF{5}` à :math:`1 \in \QQ`. En
revanche, il y a une projection canonique :math:`\ZZ \to
\GF{5}`, de sorte que la deuxième comparaison renvoie « vrai ».
Remarquez aussi que l'ordre des membres de l'égalité n'a pas
d'importance.

::

    sage: GF(5)(1) == QQ(1); QQ(1) == GF(5)(1)
    False
    False
    sage: GF(5)(1) == ZZ(1); ZZ(1) == GF(5)(1)
    True
    True
    sage: ZZ(1) == QQ(1)
    True

ATTENTION : La comparaison est plus restrictive en Sage qu'en Magma, qui
considère :math:`1 \in \GF{5}` comme égal à :math:`1 \in \QQ`.


::

    sage: magma('GF(5)!1 eq Rationals()!1')  # optional - magma
    true

Profilage (profiling)
=====================

Auteur de la section : Martin Albrecht (malb@informatik.uni-bremen.de)

    "Premature optimization is the root of all evil." - Donald Knuth
    (« L'optimisation prématurée est la source de tous les maux. »)


Il est parfois utile de rechercher dans un programme les goulets
d'étranglements qui représentent la plus grande partie du temps de
calcul : cela peut donner une idée des parties à optimiser. Cette
opération s'appelle profiler le code. Python, et donc Sage, offrent un
certain nombre de possibilités pour ce faire.

La plus simple consiste à utiliser la commande ``prun`` du shell
interactif. Elle renvoie un rapport qui résume les temps d'exécution des
fonctions les plus coûteuses. Pour profiler, par exemple, le produit de
matrices à coefficients dans un corps fini (qui, dans Sage 1.0, est
lent), on entre :

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])

.. skip

::

    sage: %prun B = A*A
           32893 function calls in 1.100 CPU seconds

    Ordered by: internal time

    ncalls tottime percall cumtime percall filename:lineno(function)
     12127  0.160   0.000   0.160  0.000 :0(isinstance)
      2000  0.150   0.000   0.280  0.000 matrix.py:2235(__getitem__)
      1000  0.120   0.000   0.370  0.000 finite_field_element.py:392(__mul__)
      1903  0.120   0.000   0.200  0.000 finite_field_element.py:47(__init__)
      1900  0.090   0.000   0.220  0.000 finite_field_element.py:376(__compat)
       900  0.080   0.000   0.260  0.000 finite_field_element.py:380(__add__)
         1  0.070   0.070   1.100  1.100 matrix.py:864(__mul__)
      2105  0.070   0.000   0.070  0.000 matrix.py:282(ncols)
      ...

Ici, ``ncalls`` désigne le nombre d'appels, ``tottime`` le temps total
passé dans une fonction (sans compter celui pris par les autres
fonctions appelées par la fonction en question), ``percall`` est le
rapport ``tottime`` divisé par ``ncalls``. ``cumtime`` donne le temps
total passé dans la fonction en comptant les appels qu'elle effectue, la
deuxième colonne ``percall`` est le quotient de ``cumtime`` par le
nombre d'appels primitifs, et ``filename:lineno(function)`` donne pour
chaque fonction le nom de fichier et le numéro de la ligne où elle est
définie. En règle générale, plus haut la fonction apparaît dans ce tableau,
plus elle est coûteuse — et donc intéressante à optimiser.

Comme d'habitude, ``prun?`` donne plus d'informations sur l'utilisation
du profileur et la signification de sa sortie.

Il est possible d'écrire les données de profilage dans un objet pour les
étudier de plus près :

.. skip

::

    sage: %prun -r A*A
    sage: stats = _
    sage: stats?

Remarque : entrer ``stats = prun -r A\*A`` à la place des deux premières
lignes ci-dessus provoque une erreur de syntaxe, car prun n'est pas une
fonction normale mais une commande du shell IPython.

Pour obtenir une jolie représentation graphique des données de
profilage, vous pouvez utiliser le profileur hotshot, un petit script
appelé ``hotshot2cachetree`` et (sous Unix uniquement) le programme
``kcachegrind``. Voici le même exemple que ci-dessus avec le profileur
hotshot :

.. skip

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])
    sage: import hotshot
    sage: filename = "pythongrind.prof"
    sage: prof = hotshot.Profile(filename, lineevents=1)

.. skip

::

    sage: prof.run("A*A")
    <hotshot.Profile instance at 0x414c11ec>
    sage: prof.close()

À ce stade le résultat est dans un fichier ``pythongrind.prof`` dans le
répertoire de travail courant. Convertissons-le au format cachegrind
pour le visualiser.

Dans le shell du système d'exploitation, tapez

.. skip

::

    hotshot2calltree -o cachegrind.out.42 pythongrind.prof

Le fichier ``cachegrind.out.42`` peut maintenant être examiné avec
``kcachegrind``. Notez qu'il est important de respecter la convention de
nommage ``cachegrind.out.XX``.

