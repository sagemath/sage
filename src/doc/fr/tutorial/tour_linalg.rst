.. _section-linalg:

Algèbre linéaire
================

Sage fournit les constructions standards d'algèbre linéaire, par exemple
le polynôme caractéristique, la forme échelonnée, la trace, diverses
décompositions, etc. d'une matrice.

La création de matrices et la multiplication matricielle sont très
faciles et naturelles :

::

    sage: A = Matrix([[1,2,3],[3,2,1],[1,1,1]])
    sage: w = vector([1,1,-4])
    sage: w*A
    (0, 0, 0)
    sage: A*w
    (-9, 1, -2)
    sage: kernel(A)
    Free module of degree 3 and rank 1 over Integer Ring
    Echelon basis matrix:
    [ 1  1 -4]

Notez bien qu'avec Sage, le noyau d'une matrice  :math:`A` est le « noyau
à gauche », c'est-à-dire l'espace des vecteurs :math:`w` tels que
:math:`wA=0`.

La résolution d'équations matricielles est facile et se fait avec la
méthode ``solve_right``. L'évaluation de ``A.solve_right(Y)`` renvoie
une matrice (ou un vecteur)  :math:`X` tel que :math:`AX=Y`:

.. link

::

    sage: Y = vector([0, -4, -1])
    sage: X = A.solve_right(Y)
    sage: X
    (-2, 1, 0)
    sage: A * X   # vérifions la réponse...
    (0, -4, -1)

S'il n'y a aucune solution, Sage renvoie une erreur :

.. skip

::

    sage: A.solve_right(w)
    Traceback (most recent call last):
    ...
    ValueError: matrix equation has no solutions

De même, il faut utiliser ``A.solve_left(Y)`` pour résoudre en :math:`X`
l'équation :math:`XA=Y`.

Sage sait aussi calculer les valeurs propres et vecteurs propres::

    sage: A = matrix([[0, 4], [-1, 0]])
    sage: A.eigenvalues ()
    [-2*I, 2*I]
    sage: B = matrix([[1, 3], [3, 1]])
    sage: B.eigenvectors_left()
    [(4, [
    (1, 1)
    ], 1), (-2, [
    (1, -1)
    ], 1)]

(La sortie de ``eigenvectors_left`` est une liste de triplets (valeur propre,
vecteur propre, multiplicité).) Sur ``QQ`` et ``RR``, on peut aussi utiliser
Maxima (voir la section :ref:`section-maxima` ci-dessous).

Comme signalé en :ref:`section-rings`, l'anneau sur lequel une matrice est
définie a une influence sur les propriétés de la matrice. Dans l'exemple
suivant, le premier argument de la commande ``matrix`` indique à Sage s'il faut
traiter la matrice comme une matrice d'entier (``ZZ``), de rationnels (``QQ``)
ou de réels (``RR``)::

    sage: AZ = matrix(ZZ, [[2,0], [0,1]])
    sage: AQ = matrix(QQ, [[2,0], [0,1]])
    sage: AR = matrix(RR, [[2,0], [0,1]])
    sage: AZ.echelon_form()
    [2 0]
    [0 1]
    sage: AQ.echelon_form()
    [1 0]
    [0 1]
    sage: AR.echelon_form()
    [ 1.00000000000000 0.000000000000000]
    [0.000000000000000  1.00000000000000]

Pour le calcul de valeurs propres et vecteurs propres sur les nombres à virgule
flottante réels ou complexes, la matrice doit être respectivement à
coefficients dans ``RDF`` (*Real Double Field*, nombres réels à précision
machine) ou ``CDF`` (*Complex Double Field*). Lorsque l'on définit une matrice
avec des coefficients flottants sans spécifier explicitement l'anneau de base,
ce ne sont pas ``RDF`` ou ``CDF`` qui sont utilisés par défaut, mais ``RR`` et
``CC``, sur lesquels ces calculs ne sont pas implémentés dans tous les cas::

    sage: ARDF = matrix(RDF, [[1.2, 2], [2, 3]])
    sage: ARDF.eigenvalues() # rel tol 8e-16
    [-0.09317121994613098, 4.293171219946131]
    sage: ACDF = matrix(CDF, [[1.2, I], [2, 3]])
    sage: ACDF.eigenvectors_right()  # rel tol 3e-15
    [(0.8818456983293743 - 0.8209140653434135*I, [(0.7505608183809549, -0.616145932704589 + 0.2387941530333261*I)], 1),
    (3.3181543016706256 + 0.8209140653434133*I, [(0.14559469829270957 + 0.3756690858502104*I, 0.9152458258662108)], 1)]

Espaces de matrices
-------------------

Créons l'espace :math:`\text{Mat}_{3\times 3}(\QQ)`:

::

    sage: M = MatrixSpace(QQ,3)
    sage: M
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field

(Pour indiquer l'espace des matrices 3 par 4, il faudrait utiliser
``MatrixSpace(QQ,3,4)``. Si le nombre de colonnes est omis, il est égal
par défaut au nombre de lignes. Ainsi ``MatrixSpace(QQ,3)`` est un
synonyme de ``MatrixSpace(QQ,3,3)``). L'espace des matrices est muni
de sa base canonique:

.. link

::

    sage: B = M.basis()
    sage: len(B)
    9
    sage: B[0,1]
    [0 1 0]
    [0 0 0]
    [0 0 0]

Nous créons une matrice comme un élément de ``M``.

.. link

::

    sage: A = M(range(9)); A
    [0 1 2]
    [3 4 5]
    [6 7 8]

Puis, nous calculons sa forme échelonnée en ligne et son noyau.

.. link

::

    sage: A.echelon_form()
    [ 1  0 -1]
    [ 0  1  2]
    [ 0  0  0]
    sage: A.kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [ 1 -2  1]

Puis nous illustrons les possibilités de calcul de matrices définies sur
des corps finis :

::

    sage: M = MatrixSpace(GF(2),4,8)
    sage: A = M([1,1,0,0, 1,1,1,1, 0,1,0,0, 1,0,1,1,
    ....:        0,0,1,0, 1,1,0,1, 0,0,1,1, 1,1,1,0])
    sage: A
    [1 1 0 0 1 1 1 1]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 1 1 1 1 1 0]
    sage: rows = A.rows()
    sage: A.columns()
    [(1, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1), (0, 0, 0, 1),
     (1, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0)]
    sage: rows
    [(1, 1, 0, 0, 1, 1, 1, 1), (0, 1, 0, 0, 1, 0, 1, 1),
     (0, 0, 1, 0, 1, 1, 0, 1), (0, 0, 1, 1, 1, 1, 1, 0)]

Nous créons le sous-espace engendré sur `\GF{2}` par les
vecteurs lignes ci-dessus.

.. link

::

    sage: V = VectorSpace(GF(2),8)
    sage: S = V.subspace(rows)
    sage: S
    Vector space of degree 8 and dimension 4 over Finite Field of size 2
    Basis matrix:
    [1 0 0 0 0 1 0 0]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 0 1 0 0 1 1]
    sage: A.echelon_form()
    [1 0 0 0 0 1 0 0]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 0 1 0 0 1 1]

La base de `S` utilisée par Sage est obtenue à partir des lignes
non nulles de la matrice des générateurs de `S` réduite sous forme
échelonnée en lignes.

Algèbre linéaire creuse
-----------------------

Sage permet de travailler avec des matrices creuses sur des anneaux
principaux.

::

    sage: M = MatrixSpace(QQ, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()

L'algorithme multi-modulaire présent dans Sage fonctionne bien pour les
matrices carrées (mais moins pour les autres) :

::

    sage: M = MatrixSpace(QQ, 50, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()
    sage: M = MatrixSpace(GF(2), 20, 40, sparse=True)
    sage: A = M.random_element()
    sage: E = A.echelon_form()

Notez que Python distingue les majuscules des minuscules :

::

    sage: M = MatrixSpace(QQ, 10,10, Sparse=True)
    Traceback (most recent call last):
    ...
    TypeError: ...__init__() got an unexpected keyword argument 'Sparse'...
