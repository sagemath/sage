Quelques mathématiques plus avancées
====================================

Géométrie algébrique
--------------------

Il est possible de définir des variétés algébriques arbitraires avec Sage,
mais les fonctionnalités non triviales sont parfois limitées aux anneaux
sur  :math:`\QQ` ou sur les corps finis. Calculons par exemple
la réunion de deux courbes planes affines, puis récupérons chaque courbe
en tant que composante irréductible de la réunion.

::

    sage: x, y = AffineSpace(2, QQ, 'xy').gens()
    sage: C2 = Curve(x^2 + y^2 - 1)
    sage: C3 = Curve(x^3 + y^3 - 1)
    sage: D = C2 + C3
    sage: D
    Affine Plane Curve over Rational Field defined by
       x^5 + x^3*y^2 + x^2*y^3 + y^5 - x^3 - y^3 - x^2 - y^2 + 1
    sage: D.irreducible_components()
    [Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x^2 + y^2 - 1,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x^3 + y^3 - 1]

Nous pouvons également trouver tous les points d'intersection des deux
courbes en les intersectant et en calculant les composantes
irréductibles.

.. link

::

    sage: V = C2.intersection(C3)
    sage: V.irreducible_components()
    [Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       y,
       x - 1,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       y - 1,
       x,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x + y + 2,
       2*y^2 + 4*y + 3]

Ainsi, par exemple,  :math:`(1,0)` et :math:`(0,1)` appartiennent aux
deux courbes (ce dont on pouvait directement s'apercevoir) ; il en va de même des
points (quadratiques), dont la coordonnée en :math:`y` satisfait à
l'équation :math:`2y^2 + 4y + 3=0`.

Sage peut calculer l'idéal torique de la cubique gauche dans l'espace
projectif de dimension 3.

::

    sage: R.<a,b,c,d> = PolynomialRing(QQ, 4)
    sage: I = ideal(b^2-a*c, c^2-b*d, a*d-b*c)
    sage: F = I.groebner_fan(); F
    Groebner fan of the ideal:
    Ideal (b^2 - a*c, c^2 - b*d, -b*c + a*d) of Multivariate Polynomial Ring
    in a, b, c, d over Rational Field
    sage: F.reduced_groebner_bases ()
    [[-c^2 + b*d, -b*c + a*d, -b^2 + a*c],
     [-b*c + a*d, -c^2 + b*d, b^2 - a*c],
     [-c^3 + a*d^2, -c^2 + b*d, b*c - a*d, b^2 - a*c],
     [-c^2 + b*d, b^2 - a*c, b*c - a*d, c^3 - a*d^2],
     [-b*c + a*d, -b^2 + a*c, c^2 - b*d],
     [-b^3 + a^2*d, -b^2 + a*c, c^2 - b*d, b*c - a*d],
     [-b^2 + a*c, c^2 - b*d, b*c - a*d, b^3 - a^2*d],
     [c^2 - b*d, b*c - a*d, b^2 - a*c]]
    sage: F.polyhedralfan()
    Polyhedral fan in 4 dimensions of dimension 4

Courbes elliptiques
-------------------

Les fonctionnalités relatives aux courbes elliptiques comprennent la
plupart des fonctionnalités de PARI, l'accès aux données des tables en
ligne de Cremona (ceci requiert le chargement d'une base de donnée
optionnelle), les fonctionnalités de mwrank, c'est-à-dire la 2-descente
avec calcul du groupe de Mordell-Weil complet, l'algorithme SEA, le
calcul de toutes les isogénies, beaucoup de nouveau code pour les
courbes sur :math:`\QQ` et une partie du code de descente
algébrique de Denis Simon.

La commande ``EllipticCurve`` permet de créer une courbe elliptique avec
beaucoup de souplesse :


-  EllipticCurve([:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6`]) : renvoie la courbe elliptique

   .. math::  y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6,

   où les :math:`a_i`'s sont convertis par coercition dans le parent
   de :math:`a_1`. Si tous les :math:`a_i` ont pour parent
   :math:`\ZZ`, ils sont convertis par coercition dans
   :math:`\QQ`.

-  EllipticCurve([:math:`a_4`, :math:`a_6`]) : idem
   avec :math:`a_1=a_2=a_3=0`.

-  EllipticCurve(label) : Renvoie la courbe elliptique sur  :math:`\QQ`  de la
   base de données de Cremona selon son nom dans la (nouvelle !)
   nomenclature de Cremona. Les courbes sont étiquetées par une chaîne de
   caractère telle que ``"11a"`` ou ``"37b2"``. La lettre doit être en
   minuscule (pour faire la différence avec l'ancienne nomenclature).

-  EllipticCurve(j) : renvoie une courbe elliptique de
   :math:`j`-invariant :math:`j`.

-  EllipticCurve(R, [:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`,  :math:`a_6`]) : Crée la courbe elliptique sur l'anneau :math:`R` donnée
   par les coefficients :math:`a_i` comme ci-dessus.


Illustrons chacune de ces constructions :

::

    sage: EllipticCurve([0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve([GF(5)(0),0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    sage: EllipticCurve([1,2])
    Elliptic Curve defined by y^2  = x^3 + x + 2 over Rational Field

    sage: EllipticCurve('37a')
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve_from_j(1)
    Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455 over Rational Field

    sage: EllipticCurve(GF(5), [0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

Le couple :math:`(0,0)` est un point de la courbe elliptique :math:`E`
définie par :math:`y^2 + y = x^3 - x`. Pour créer ce point avec Sage, il
convient de taper ``E([0,0])``. Sage peut additionner des points sur une
telle courbe elliptique (rappelons qu'une courbe elliptique possède
une structure de groupe additif où le point à l'infini représente
l'élément neutre et où trois points alignés de la courbe sont de somme
nulle) :

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
    sage: P = E([0,0])
    sage: P + P
    (1 : 0 : 1)
    sage: 10*P
    (161/16 : -2065/64 : 1)
    sage: 20*P
    (683916417/264517696 : -18784454671297/4302115807744 : 1)
    sage: E.conductor()
    37

Les courbes elliptiques sur les nombres complexes sont paramétrées par
leur   :math:`j`-invariant. Sage calcule le :math:`j`-invariant comme
suit :

::

    sage: E = EllipticCurve([0,0,0,-4,2]); E
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: E.conductor()
    2368
    sage: E.j_invariant()
    110592/37

Si l'on fabrique une courbe avec le même :math:`j`-invariant que celui
de :math:`E`, elle n'est pas nécessairement isomorphe à :math:`E`. Dans
l'exemple suivant, les courbes ne sont pas isomorphes parce que leur
conducteur est différent.

::

    sage: F = EllipticCurve_from_j(110592/37)
    sage: F.conductor()
    37

Toutefois, le twist de :math:`F` par 2 donne une courbe isomorphe.

.. link

::

    sage: G = F.quadratic_twist(2); G
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: G.conductor()
    2368
    sage: G.j_invariant()
    110592/37

On peut calculer les coefficients :math:`a_n` de la série-:math:`L` ou
forme modulaire :math:`\sum_{n=0}^\infty a_nq^n` attachée à une courbe
elliptique.  Le calcul s'effectue en utilisant la bibliothèque PARI
écrite en C :

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E.anlist(30)
    [0, 1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0, -4,
     3, 10, 2, 0, -1, 4, -9, -2, 6, -12]
    sage: v = E.anlist(10000)

Il faut à peine quelques secondes pour calculer tous les coefficients
:math:`a_n` pour :math:`n\leq 10^5`:

.. skip

::

    sage: %time v = E.anlist(100000)
    CPU times: user 0.98 s, sys: 0.06 s, total: 1.04 s
    Wall time: 1.06

Les courbes elliptiques peuvent être construites en utilisant leur nom
dans la nomenclature de Cremona. Ceci charge par avance la courbe
elliptique avec les informations la concernant, telles que son rang, son
nombre de Tamagawa, son régulateur, etc.

::

    sage: E = EllipticCurve("37b2")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational
    Field
    sage: E = EllipticCurve("389a")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x  over Rational Field
    sage: E.rank()
    2
    sage: E = EllipticCurve("5077a")
    sage: E.rank()
    3

On peut aussi accéder à la base de données de Cremona directement.

::

    sage: db = sage.databases.cremona.CremonaDatabase()
    sage: db.curves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1], 'b1': [[0, 1, 1, -23, -50], 0, 3]}
    sage: db.allcurves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1],
     'b1': [[0, 1, 1, -23, -50], 0, 3],
     'b2': [[0, 1, 1, -1873, -31833], 0, 1],
     'b3': [[0, 1, 1, -3, 1], 0, 3]}

Les objets extraits de la base de données ne sont pas de type
``EllipticCurve``, mais de simples entrées de base de données formées de
quelques champs. Par défaut, Sage est distribué avec une version réduite
de la base de données de Cremona qui ne contient que des informations
limitées sur les courbes elliptiques de conducteur :math:`\leq 10000`.
Il existe également en option une version plus complète qui contient des
données étendues portant sur toute les courbes de conducteur jusqu'à
:math:`120000` (à la date d'octobre 2005). Une autre - énorme (2GB) -
base de données optionnelle, fournie dans un package séparé, contient
des centaines de millions de courbes elliptiques de la bases de donnée de
Stein-Watkins.

Caractères de Dirichlet
-----------------------

Un *caractère de Dirichlet* est une extension d'un homomorphisme
:math:`(\ZZ/N\ZZ)^* \to R^*`, pour un certain anneau
:math:`R`, à l'application :math:`\ZZ \to R` obtenue en envoyant
les entiers :math:`x` tels que  :math:`\gcd(N,x)>1` vers 0.

::

    sage: G = DirichletGroup(12)
    sage: G.list()
    [Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1,
    Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1]
    sage: G.gens()
    (Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1)
    sage: len(G)
    4

Une fois le groupe créé, on crée aussitôt un élément et on calcule avec lui.

.. link

::

    sage: G = DirichletGroup(21)
    sage: chi = G.1; chi
    Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6
    sage: chi.values()
    [0, 1, zeta6 - 1, 0, -zeta6, -zeta6 + 1, 0, 0, 1, 0, zeta6, -zeta6, 0, -1,
     0, 0, zeta6 - 1, zeta6, 0, -zeta6 + 1, -1]
    sage: chi.conductor()
    7
    sage: chi.modulus()
    21
    sage: chi.order()
    6
    sage: chi(19)
    -zeta6 + 1
    sage: chi(40)
    -zeta6 + 1

Il est possible aussi de calculer l'action d'un groupe de Galois
:math:`\text{Gal}(\QQ(\zeta_N)/\QQ)` sur l'un de ces
caractères, de même qu'une décomposition en produit direct correspondant
à la factorisation du module.

.. link

::

    sage: chi.galois_orbit()
    [Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> -zeta6 + 1,
     Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6]

    sage: go = G.galois_orbits()
    sage: [len(orbit) for orbit in go]
    [1, 2, 2, 1, 1, 2, 2, 1]

    sage: G.decomposition()
    [Group of Dirichlet characters modulo 3 with values in Cyclotomic Field of order 6 and degree 2,
     Group of Dirichlet characters modulo 7 with values in Cyclotomic Field of order 6 and degree 2]

Construisons à present le groupe de caractères de Dirichlet modulo 20,
mais à valeur dans  :math:`\QQ(i)`:

::

    sage: K.<i> = NumberField(x^2+1)
    sage: G = DirichletGroup(20,K)
    sage: G
    Group of Dirichlet characters modulo 20 with values in Number Field in i with defining polynomial x^2 + 1

Nous calculons ensuite différents invariants de ``G``:

.. link

::

    sage: G.gens()
    (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1,
    Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> i)

    sage: G.unit_gens()
    (11, 17)
    sage: G.zeta()
    i
    sage: G.zeta_order()
    4

Dans cet exemple, nous créons un caractère de Dirichlet à valeurs dans
un corps de nombres. Nous spécifions ci-dessous explicitement le choix
de la racine de l'unité par le troisième argument de la fonction
``DirichletGroup``.

::

    sage: x = polygen(QQ, 'x')
    sage: K = NumberField(x^4 + 1, 'a'); a = K.0
    sage: b = K.gen(); a == b
    True
    sage: K
    Number Field in a with defining polynomial x^4 + 1
    sage: G = DirichletGroup(5, K, a); G
    Group of Dirichlet characters modulo 5 with values in the group of order 8 generated by a in Number Field in a with defining polynomial x^4 + 1
    sage: chi = G.0; chi
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> a^2
    sage: [(chi^i)(2) for i in range(4)]
    [1, a^2, -1, -a^2]

Ici, ``NumberField(x^4 + 1, 'a')`` indique à Sage d'utiliser le symbole
"a" dans l'affichage de ce qu'est ``K`` (un corps de nombre en "a"
défini par le polynôme :math:`x^4 + 1`). Le nom "a" n'est pas déclaré à
ce point. Une fois que  ``a = K.0`` (ou de manière équivalente ``a =
K.gen()``) est évalué, le symbole "a" représente une racine du polynôme
générateur :math:`x^4+1`.

Formes modulaires
-----------------

Sage peut accomplir des calculs relatifs aux formes modulaires,
notamment des calculs de dimension, d'espace de symboles modulaires, d'opérateurs de Hecke et de décomposition.

Il y a plusieurs fonctions disponibles pour calculer la dimension
d'espaces de formes modulaires. Par exemple,

::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(11),2)
    1
    sage: dimension_cusp_forms(Gamma0(1),12)
    1
    sage: dimension_cusp_forms(Gamma1(389),2)
    6112

Nous illustrons ci-dessous le calcul des opérateurs de Hecke sur un
espace de symboles modulaires de niveau :math:`1` et de poids
:math:`12`.

::

    sage: M = ModularSymbols(1,12)
    sage: M.basis()
    ([X^8*Y^2,(0,0)], [X^9*Y,(0,0)], [X^10,(0,0)])
    sage: t2 = M.T(2)
    sage: t2
    Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1)
    of weight 12 with sign 0 over Rational Field
    sage: t2.matrix()
    [ -24    0    0]
    [   0  -24    0]
    [4860    0 2049]
    sage: f = t2.charpoly('x'); f
    x^3 - 2001*x^2 - 97776*x - 1180224
    sage: factor(f)
    (x - 2049) * (x + 24)^2
    sage: M.T(11).charpoly('x').factor()
    (x - 285311670612) * (x - 534612)^2

Nous pouvons aussi créer des espaces pour :math:`\Gamma_0(N)` et
:math:`\Gamma_1(N)`.

::

    sage: ModularSymbols(11,2)
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
     0 over Rational Field
    sage: ModularSymbols(Gamma1(11),2)
    Modular Symbols space of dimension 11 for Gamma_1(11) of weight 2 with
    sign 0 over Rational Field

Calculons quelques polynômes caractéristiques et développements en série de Fourier.

::

    sage: M = ModularSymbols(Gamma1(11),2)
    sage: M.T(2).charpoly('x')
    x^11 - 8*x^10 + 20*x^9 + 10*x^8 - 145*x^7 + 229*x^6 + 58*x^5 - 360*x^4
         + 70*x^3 - 515*x^2 + 1804*x - 1452
    sage: M.T(2).charpoly('x').factor()
    (x - 3) * (x + 2)^2 * (x^4 - 7*x^3 + 19*x^2 - 23*x + 11)
            * (x^4 - 2*x^3 + 4*x^2 + 2*x + 11)
    sage: S = M.cuspidal_submodule()
    sage: S.T(2).matrix()
    [-2  0]
    [ 0 -2]
    sage: S.q_expansion_basis(10)
    [q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)]

On peut même calculer des espaces de formes modulaires avec caractères.

::

    sage: G = DirichletGroup(13)
    sage: e = G.0^2
    sage: M = ModularSymbols(e,2); M
    Modular Symbols space of dimension 4 and level 13, weight 2, character
    [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
    sage: M.T(2).charpoly('x').factor()
    (x - zeta6 - 2) * (x - 2*zeta6 - 1) * (x + zeta6 + 1)^2
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 2 of Modular Symbols space of
    dimension 4 and level 13, weight 2, character [zeta6], sign 0, over
    Cyclotomic Field of order 6 and degree 2
    sage: S.T(2).charpoly('x').factor()
    (x + zeta6 + 1)^2
    sage: S.q_expansion_basis(10)
    [q + (-zeta6 - 1)*q^2 + (2*zeta6 - 2)*q^3 + zeta6*q^4 + (-2*zeta6 + 1)*q^5 + (-2*zeta6 + 4)*q^6 + (2*zeta6 - 1)*q^8 - zeta6*q^9 + O(q^10)]

Voici un autre exemple montrant comment Sage peut calculer l'action d'un
opérateur de Hecke sur un espace de formes modulaires.

::

    sage: T = ModularForms(Gamma0(11),2)
    sage: T
    Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of
    weight 2 over Rational Field
    sage: T.degree()
    2
    sage: T.level()
    11
    sage: T.group()
    Congruence Subgroup Gamma0(11)
    sage: T.dimension()
    2
    sage: T.cuspidal_subspace()
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
    Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: T.eisenstein_subspace()
    Eisenstein subspace of dimension 1 of Modular Forms space of dimension 2
    for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: M = ModularSymbols(11); M
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
    0 over Rational Field
    sage: M.weight()
    2
    sage: M.basis()
    ((1,0), (1,8), (1,9))
    sage: M.sign()
    0

Notons :math:`T_p` les opérateurs de Hecke usuels  (:math:`p` premier).
Comment agissent les opérateurs de Hecke :math:`T_2`, :math:`T_3`,
:math:`T_5`  sur l'espace des symboles modulaires ?

.. link

::

    sage: M.T(2).matrix()
    [ 3  0 -1]
    [ 0 -2  0]
    [ 0  0 -2]
    sage: M.T(3).matrix()
    [ 4  0 -1]
    [ 0 -1  0]
    [ 0  0 -1]
    sage: M.T(5).matrix()
    [ 6  0 -1]
    [ 0  1  0]
    [ 0  0  1]
