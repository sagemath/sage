.. _chapter-help:

Obtenir de l'aide
=================

Sage est doté d'une importante documentation intégrée, accessible en
tapant (par exemple) le nom d'une fonction ou d'une constante suivi d'un
point d'interrogation :

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring:

        The tangent function

        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring:

        The natural logarithm of the real number 2.

        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 32-bit
            0.69314718055994530941723212145817656807   # 64-bit
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <... 'function'>
    Definition:  sudoku(A)
    Docstring:

        Solve the 9x9 Sudoku puzzle defined by the matrix A.

        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

Sage dispose aussi de la complétion de ligne de commande, accessible en
tapant les quelques premières lettres du nom d'une fonction puis en
appuyant sur la touche tabulation. Ainsi, si vous tapez ``ta`` suivi de
``TAB``, Sage affichera ``tachyon, tan, tanh, taylor``. C'est une façon
commode de voir quels noms de fonctions et d'autres structures sont
disponibles en Sage.


.. _section-functions:

Fonctions, indentation et itération
===================================

Les définitions de fonctions en Sage sont introduites par la commande
``def``, et la liste des noms des paramètres est suivie de deux points,
comme dans :

::

    sage: def is_even(n):
    ....:     return n%2 == 0
    sage: is_even(2)
    True
    sage: is_even(3)
    False

Remarque : suivant la version du *notebook* que vous utilisez, il est
possible que vous voyez trois points ``....:`` au début de la deuxième
ligne de l'exemple. Ne les entrez pas, ils servent uniquement à signaler
que le code est indenté.

Les types des paramètres ne sont pas spécifiés dans la définition de la
fonction. Il peut y avoir plusieurs paramètres, chacun accompagné
optionnellement d'une valeur par défaut. Par exemple, si la valeur de
``divisor`` n'est pas donnée lors d'un appel à la fonction ci-dessous,
la valeur par défaut ``divisor=2`` est utilisée.

::

    sage: def is_divisible_by(number, divisor=2):
    ....:     return number%divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False

Il est possible de spécifier un ou plusieurs des paramètres par leur nom
lors de l'appel de la fonction ; dans ce cas, les paramètres nommés
peuvent apparaître dans n'importe quel ordre :

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True

En Python, contrairement à de nombreux autres langages, les blocs de
code ne sont pas délimités par des accolades ou des mots-clés de début
et de fin de bloc. Au lieu de cela, la structure des blocs est donnée
par l'indentation, qui doit être la même dans tout le bloc. Par exemple,
le code suivant déclenche une erreur de syntaxe parce que l'instruction
``return`` n'est pas au même niveau d'indentation que les lignes
précédentes.

.. skip

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:    return v
    Syntax Error:
           return v

Une fois l'indentation corrigée, l'exemple fonctionne :

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:     return v
    sage: even(10)
    [4, 6, 8]

Il n'y a pas besoin de placer des points-virgules en fin de ligne ; une
instruction est en général terminée par un passage à la ligne. En
revanche, il est possible de placer plusieurs instructions sur la même
ligne en les séparant par des points-virgules :

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

Pour continuer une instruction sur la ligne suivante, placez une barre
oblique inverse en fin de ligne :

::

    sage: 2 + \
    ....:    3
    5

Pour compter en Sage, utilisez une boucle dont la variable d'itération
parcourt une séquence d'entiers. Par exemple, la première ligne
ci-dessous a exactement le même effet que ``for(i=0; i<3; i++)`` en C++
ou en Java :

::

    sage: for i in range(3):
    ....:     print(i)
    0
    1
    2

La première ligne ci-dessous correspond à ``for(i=2;i<5;i++)``.

::

    sage: for i in range(2,5):
    ....:     print(i)
    2
    3
    4

Le troisième paramètre contrôle le pas de l'itération. Ainsi, ce qui
suit est équivalent à ``for(i=1;i<6;i+=2)``.

::

    sage: for i in range(1,6,2):
    ....:     print(i)
    1
    3
    5

Vous souhaiterez peut-être regrouper dans un joli tableau les résultats
numériques que vous aurez calculés avec Sage. Une façon de faire commode
utilise les chaînes de format. Ici, nous affichons une table des carrés
et des cubes en trois colonnes, chacune d'une largeur de six caractères.

::

    sage: for i in range(5):
    ....:     print('%6s %6s %6s' % (i, i^2, i^3))
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64

La structure de données de base de Sage est la liste, qui est — comme
son nom l'indique — une liste d'objets arbitraires. Voici un exemple
de liste::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

Comme dans de nombreux langages de programmation, les listes sont
indexées à partir de 0.

.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)

La fonction ``len(v)`` donne la longueur de ``v``....:``v.append(obj)``
ajoute un nouvel objet à la fin de ``v`` ; et ``del v[i]`` supprime
l'élément d'indice ``i`` de ``v``.

.. link

::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]

Une autre structure de données importante est le dictionnaire (ou
tableau associatif). Un dictionnaire fonctionne comme une liste, à ceci
près que les indices peuvent être presque n'importe quels objets (les
objets mutables sont interdits) :

::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi

Vous pouvez définir de nouveaux types de données en utilisant les
classes. Encapsuler les objets mathématiques dans des classes représente
une technique puissante qui peut vous aider à simplifier et organiser
vos programmes Sage. Dans l'exemple suivant, nous définissons une classe
qui représente la liste des entiers impairs strictement positifs jusqu'à
*n*. Cette classe dérive du type interne ``list``.

::

    sage: class Evens(list):
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:         list.__init__(self, range(2, n+1, 2))
    ....:     def __repr__(self):
    ....:         return "Even positive numbers up to n."

La méthode ``__init__`` est appelée à la création de l'objet pour
l'initialiser ; la méthode ``__repr__`` affiche l'objet. À la seconde
ligne de la méthode ``__init__``, nous appelons le constructeur de la
classe ``list``. Pour créer un objet de classe ``Evens``, nous procédons
ensuite comme suit :

.. link

::

    sage: e = Evens(10)
    sage: e
    Even positive numbers up to n.

Notez que ``e`` s'affiche en utilisant la méthode ``__repr__`` que nous avons
définie plus haut. Pour voir la liste de nombres sous-jacente, on utilise
la fonction ``list`` :

.. link

::

    sage: list(e)
    [2, 4, 6, 8, 10]

Il est aussi possible d'accéder à l'attribut ``n``, ou encore d'utiliser
``e`` en tant que liste.

.. link

::

    sage: e.n
    10
    sage: e[2]
    6
