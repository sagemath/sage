.. _section-poly:

Polynômes
=========

Dans cette partie, nous expliquons comment créer et utiliser des
polynômes avec Sage.


.. _section-univariate:

Polynômes univariés
--------------------

Il existe trois façons de créer des anneaux de polynômes.

::

    sage: R = PolynomialRing(QQ, 't')
    sage: R
    Univariate Polynomial Ring in t over Rational Field

Ceci crée un anneau de polynômes et indique à Sage d'utiliser (la chaîne
de caractère) 't' comme indéterminée lors de l'affichage à l'écran.
Toutefois, ceci ne définit pas le symbole  ``t`` pour son utilisation
dans Sage. Aussi, il n'est pas possible de l'utiliser pour saisir un
polynôme  (comme :math:`t^2+1`) qui appartient à ``R``.

Une deuxième manière de procéder est

.. link

::

    sage: S = QQ['t']
    sage: S == R
    True

Ceci a les mêmes effets en ce qui concerne ``t``.

Une troisième manière de procéder, très pratique, consiste à entrer

::

    sage: R.<t> = PolynomialRing(QQ)

ou

::

    sage: R.<t> = QQ['t']

ou même

::

    sage: R.<t> = QQ[]

L'effet secondaire de ces dernières instructions est de définir la
variable ``t`` comme l'indéterminée de l'anneau de polynômes. Ceci
permet de construire très aisément des éléments de ``R``, comme décrit
ci-après. (Noter que cette troisième manière est très semblable à la
notation par constructeur de Magma et que, de même que dans Magma, ceci
peut servir pour une très large classe d'objets.)

.. link

::

    sage: poly = (t+1) * (t+2); poly
    t^2 + 3*t + 2
    sage: poly in R
    True

Quelle que soit la méthode utilisée pour définir l'anneau de polynômes, on
récupère l'indéterminée comme le :math:`0`-ième générateur :

::

    sage: R = PolynomialRing(QQ, 't')
    sage: t = R.0
    sage: t in R
    True

Notez que les nombres complexes peuvent être construits de façon
similaire : les nombres complexes peuvent être vus comme engendrés sur
les réels par le symbole ``i``. Aussi, on dispose de :

::

    sage: CC
    Complex Field with 53 bits of precision
    sage: CC.0  # 0ième générateur CC
    1.00000000000000*I

Pour un anneau de polynômes, on peut obtenir à la fois l'anneau et son
générateur ou juste le générateur au moment de la création de l'anneau
comme suit :

::

    sage: R, t = QQ['t'].objgen()
    sage: t    = QQ['t'].gen()
    sage: R, t = objgen(QQ['t'])
    sage: t    = gen(QQ['t'])

Finalement, on peut faire de l'arithmétique dans  :math:`\QQ[t]`.

::

    sage: R, t = QQ['t'].objgen()
    sage: f = 2*t^7 + 3*t^2 - 15/19
    sage: f^2
    4*t^14 + 12*t^9 - 60/19*t^7 + 9*t^4 - 90/19*t^2 + 225/361
    sage: cyclo = R.cyclotomic_polynomial(7); cyclo
    t^6 + t^5 + t^4 + t^3 + t^2 + t + 1
    sage: g = 7 * cyclo * t^5 * (t^5 + 10*t + 2)
    sage: g
    7*t^16 + 7*t^15 + 7*t^14 + 7*t^13 + 77*t^12 + 91*t^11 + 91*t^10 + 84*t^9
           + 84*t^8 + 84*t^7 + 84*t^6 + 14*t^5
    sage: F = factor(g); F
    (7) * t^5 * (t^5 + 10*t + 2) * (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)
    sage: F.unit()
    7
    sage: list(F)
    [(t, 5), (t^5 + 10*t + 2, 1), (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1, 1)]

On remarquera que la factorisation prend correctement en compte le
coefficient dominant, et ne l'oublie pas dans le résultat.

S'il arrive que vous utilisiez intensivement, par exemple, la fonction
``R.cyclotomic_polynomial``  dans un projet de recherche
quelconque, en plus de citer Sage, vous devriez chercher à quel
composant Sage fait appel pour calculer en réalité ce polynôme
cyclotomique et citer ce composant. Dans ce cas particulier, en tapant
``R.cyclotomic_polynomial??`` pour voir le code source, vous verriez
rapidement une ligne telle que ``f = pari.polcyclo(n)`` ce qui signifie
que PARI est utilisé pour le calcul du polynôme cyclotomique. Pensez à
citer PARI dans votre travail.

La division d'un polynôme par un autre produit un élément du corps des
fractions, que Sage crée automatiquement.

::

    sage: x = QQ['x'].0
    sage: f = x^3 + 1; g = x^2 - 17
    sage: h = f/g;  h
    (x^3 + 1)/(x^2 - 17)
    sage: h.parent()
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

En utilisant des séries de Laurent, on peut calculer des développements
en série dans le corps des fractions de ``QQ[x]``:

::

    sage: R.<x> = LaurentSeriesRing(QQ); R
    Laurent Series Ring in x over Rational Field
    sage: 1/(1-x) + O(x^10)
    1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)

Si l'on nomme les variables différemment, on obtient un anneau de
polynômes univariés différent.

::

    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<y> = PolynomialRing(QQ)
    sage: x == y
    False
    sage: R == S
    False
    sage: R(y)
    x
    sage: R(y^2 - 17)
    x^2 - 17

L'anneau est déterminé par sa variable. Notez que créer un autre anneau
avec la même variable ``x`` ne renvoie pas de nouvel anneau.

::

    sage: R = PolynomialRing(QQ, "x")
    sage: T = PolynomialRing(QQ, "x")
    sage: R == T
    True
    sage: R is T
    True
    sage: R.0 == T.0
    True

Sage permet aussi de travailler dans des anneaux de séries formelles et
de séries de Laurent sur un anneau de base quelconque. Dans l'exemple
suivant, nous créons un élément de  :math:`\GF{7}[[T]]` et
effectuons une division pour obtenir un élément de
:math:`\GF{7}((T))`.

::

    sage: R.<T> = PowerSeriesRing(GF(7)); R
    Power Series Ring in T over Finite Field of size 7
    sage: f = T  + 3*T^2 + T^3 + O(T^4)
    sage: f^3
    T^3 + 2*T^4 + 2*T^5 + O(T^6)
    sage: 1/f
    T^-1 + 4 + T + O(T^2)
    sage: parent(1/f)
    Laurent Series Ring in T over Finite Field of size 7

On peut aussi créer des anneaux de séries formelles en utilisant des
doubles crochets :

::

    sage: GF(7)[['T']]
    Power Series Ring in T over Finite Field of size 7

Polynômes multivariés
---------------------

Pour travailler avec des polynômes à plusieurs variables, on commence
par déclarer l'anneau des polynômes et les variables, de l'une des deux
manières suivantes.


::

    sage: R = PolynomialRing(GF(5),3,"z") # here, 3 = number of variables
    sage: R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

De même que pour les polynômes à une seule variable, les variantes
suivantes sont autorisées :

::

    sage: GF(5)['z0, z1, z2']
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5
    sage: R.<z0,z1,z2> = GF(5)[]; R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Si l'on désire de simples lettres comme noms de variables, on
peut utiliser les raccourcis suivants :

::

    sage: PolynomialRing(GF(5), 3, 'xyz')
    Multivariate Polynomial Ring in x, y, z over Finite Field of size 5

A présent, passons aux questions arithmétiques.

::

    sage: z = GF(5)['z0, z1, z2'].gens()
    sage: z
    (z0, z1, z2)
    sage: (z[0]+z[1]+z[2])^2
    z0^2 + 2*z0*z1 + z1^2 + 2*z0*z2 + 2*z1*z2 + z2^2

On peut aussi utiliser des notations plus mathématiques pour construire
un anneau de polynômes.

::

    sage: R = GF(5)['x,y,z']
    sage: x,y,z = R.gens()
    sage: QQ['x']
    Univariate Polynomial Ring in x over Rational Field
    sage: QQ['x,y'].gens()
    (x, y)
    sage: QQ['x'].objgens()
    (Univariate Polynomial Ring in x over Rational Field, (x,))

Sous Sage, les polynômes multivariés sont implémentés en représentation
« distributive » (par opposition à récursive), à l'aide de dictionnaires
Python. Sage a souvent recours à Singular [Si]_, par exemple, pour le
calcul de pgcd ou de bases de Gröbner d'idéaux.

::

    sage: R, (x, y) = PolynomialRing(RationalField(), 2, 'xy').objgens()
    sage: f = (x^3 + 2*y^2*x)^2
    sage: g = x^2*y^2
    sage: f.gcd(g)
    x^2

Créons ensuite l'idéal :math:`(f,g)` engendré par  :math:`f` et
:math:`g`, en multipliant simplement ``(f,g)`` par ``R`` (nous pourrions
aussi bien écrire ``ideal([f,g])`` ou ``ideal(f,g)``).

.. link

::

    sage: I = (f, g)*R; I
    Ideal (x^6 + 4*x^4*y^2 + 4*x^2*y^4, x^2*y^2) of Multivariate Polynomial
    Ring in x, y over Rational Field
    sage: B = I.groebner_basis(); B
    [x^6, x^2*y^2]
    sage: x^2 in I
    False

En passant, la base de Gröbner ci-dessus n'est pas une liste mais
une suite non mutable. Ceci signifie qu'elle possède un univers, un
parent, et qu'elle ne peut pas être modifiée (ce qui est une bonne chose
puisque changer la base perturberait d'autres routines qui utilisent la
base de Gröbner).

.. link

::

    sage: B.universe()
    Multivariate Polynomial Ring in x, y over Rational Field
    sage: B[1] = x
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Un peu (comprenez : pas assez à notre goût) d'algèbre commutative est
disponible en Sage. Ces routines font appel à Singular. Par exemple, il
est possible de calculer la décomposition en facteurs premiers et les
idéaux premiers associés de :math:`I`:

.. link

::

    sage: I.primary_decomposition()
    [Ideal (x^2) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y^2, x^6) of Multivariate Polynomial Ring in x, y over Rational Field]
    sage: I.associated_primes()
    [Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field]
