.. -*- coding: utf-8 -*-

.. _section-coercion:

=================================
Parents, conversions, coercitions
=================================

Cette section peut paraître plus technique que celles qui précèdent, mais nous
pensons qu'il est important de comprendre ce que sont les parents et les
coercitions pour utiliser comme il faut les structures algébriques fournies par
Sage.

Nous allons voir ici ce que ces notions signifient, mais pas comment les mettre
en œuvre pour implémenter une nouvelle structure algébrique. Un tutoriel
thématique couvrant ce point est disponible `ici <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_.

Éléments
--------

Une première approximation en Python de la notion mathématique d'anneau
pourrait consister à définir une classe pour les éléments ``X`` de l'anneau
concerné, de fournir les méthodes « double-underscore » nécessaires pour donner
un sens aux opérations de l'anneau, par exemple  ``__add__``, ``__sub__`` et
``__mul__``, et naturellement de s'assurer qu'elles respectent les axiomes de
la structure d'anneau.

Python étant un language (dynamiquement) fortement typé, on pourrait s'attendre
à devoir implémenter une classe pour chaque anneau. Après tout, Python définit
bien un type ``<int>`` pour les entiers, un type ``<float>`` pour les réels, et
ainsi de suite. Mais cette approche ne peut pas fonctionner : il y a une
infinité d'anneaux différents, et l'on ne peut pas implémenter une infinité de
classes !

Une autre idée est de créer une hiérarchie de classes destinées à implémenter
les éléments des structures algébriques usuelles : éléments de groupes,
d'anneaux, d'algèbres à division, d'anneaux commutatifs, de corps, d'algèbres,
etc.

Mais cela signifie que des éléments d'anneaux franchement différents peuvent
avoir le même type.

::

    sage: P.<x,y> = GF(3)[]
    sage: Q.<a,b> = GF(4,'z')[]
    sage: type(x)==type(a)
    True

On pourrait aussi vouloir avoir des classes Python différentes pour fournir
plusieurs implémentations d'une même structure mathématique (matrices denses
contre matrices creuses par exemple).

::

    sage: P.<a> = PolynomialRing(ZZ)
    sage: Q.<b> = PolynomialRing(ZZ, sparse=True)
    sage: R.<c> = PolynomialRing(ZZ, implementation='NTL')
    sage: type(a); type(b); type(c)
    <class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
    <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_integral_domain_with_category.element_class'>
    <class 'sage.rings.polynomial.polynomial_integer_dense_ntl.Polynomial_integer_dense_ntl'>

Deux problèmes se posent alors. D'une part, si deux éléments sont instances de
la même classe, on s'attend à ce que leur méthode ``__add__`` soit capable de
les additionner, alors que ce n'est pas ce que l'on souhaite si les éléments
appartiennent en fait à des anneaux différents. D'autre part, si l'on a deux
éléments qui appartiennent à des implémentations différentes d'un même anneau,
on veut pouvoir les ajouter, et ce n'est pas immédiats s'ils ne sont pas
instances de la même classe.

La solution à ces difficultés est fournie par le mécanisme de *coercition*
décrit ci-dessous.

Mais avant tout, il est essentiel que chaque élément « sache » de quoi il est
élément. Cette information est donnée par la méthode ``parent()``.

.. link

::

    sage: a.parent(); b.parent(); c.parent()
    Univariate Polynomial Ring in a over Integer Ring
    Sparse Univariate Polynomial Ring in b over Integer Ring
    Univariate Polynomial Ring in c over Integer Ring (using NTL)


Parents et catégories
---------------------

En plus d'une hiérarchie de classes destinée à implémenter les éléments de
structures algébriques, Sage fournit une hiérarchie similaire pour les
structures elles-mêmes. Ces structures s'appellent en Sage des *parents*, et
leurs classes dérivent d'une même classe de base. Celle-ci a des sous-classes
« ensemble », « anneau », « corps », et ainsi de suite, dont la hiérarchie
correspond à peu près à celle des concepts mathématiques qu'elles décrivent :

::

    sage: isinstance(QQ,Field)
    True
    sage: isinstance(QQ, Ring)
    True
    sage: isinstance(ZZ,Field)
    False
    sage: isinstance(ZZ, Ring)
    True

Or en algèbre, on regroupe les objets qui partagent le même genre de structure
algébrique en ce que l'on appelle des *catégories*. Il y a donc un parallèle
approximatif entre la hiérarchie des classes de Sage et la hiérarchie des
catégories. Mais cette correspondance n'est pas parfaite, et Sage implémente
par ailleurs les catégories en tant que telles :

::

    sage: Rings()
    Category of rings
    sage: ZZ.category()
    Join of Category of Dedekind domains
        and Category of euclidean domains
        and Category of noetherian rings
        and Category of infinite enumerated sets
        and Category of metric spaces
    sage: ZZ.category().is_subcategory(Rings())
    True
    sage: ZZ in Rings()
    True
    sage: ZZ in Fields()
    False
    sage: QQ in Fields()
    True

Tandis que la hiérarchie des classes est déterminée avant tout par des
considérations de programmation, l'infrastructure des catégories cherche plutôt
à respecter la structure mathématique. Elle permet de munir les objets d'une
catégorie de méthodes et de tests génériques, qui ne dépendent pas de
l'implémentation particulière d'un objet donné de la catégorie.

Les parents en tant qu'objets Python doivent être uniques. Ainsi, lorsqu'un
anneau de polynômes sur un anneau donné et avec une liste donnée de générateurs
est construit, il est conservé en cache et réutilisé par la suite :

::

    sage: RR['x','y'] is RR['x','y']
    True


Types et parents
----------------

Le type ``RingElement`` ne correspond pas parfaitement à la notion
mathématique d'élément d'anneau. Par exemple, bien que les matrices carrées
appartiennent à un anneau, elles ne sont pas de type ``RingElement`` :

::

    sage: M = Matrix(ZZ,2,2); M
    [0 0]
    [0 0]
    sage: isinstance(M, RingElement)
    False

Si les *parents* sont censés être uniques, des *éléments* égaux d'un parent ne
sont pas nécessairement identiques. Le comportement de Sage diffère ici de
celui de Python pour certains entiers (pas tous) :

::

    sage: int(1) is int(1) # Python int
    True
    sage: int(-15) is int(-15)
    False
    sage: 1 is 1           # Sage Integer
    False

Il faut bien comprendre que les éléments d'anneaux différents ne se distinguent
généralement pas par leur type, mais par leur parent :

::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: type(a) is type(b)
    True
    sage: parent(a)
    Finite Field of size 2
    sage: parent(b)
    Finite Field of size 5

Ainsi, **le parent d'un élément est plus important que son type** du point de
vue algébrique.

Conversion et coercition
-------------------------

Il est parfois possible de convertir un élément d'un certain parent en élément
d'un autre parent. Une telle conversion peut être explicite ou implicite. Les
conversions implicites sont appelées *coercitions*.

Le lecteur aura peut-être rencontré les notions de *conversion de type* et de
*coercition de type* dans le contexte du langage C par exemple. En Sage, il
existe aussi des notions de conversion et de coercition, mais elles
s'appliquent aux *parents* et non aux types. Attention donc à ne pas confondre
les conversions en Sage avec les conversions de type du C !

Nous nous limitons ici à une brève présentation, et renvoyons le lecteur à la
section du manuel de référence consacrée aux coercitions ainsi qu'au
`tutoriel <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_
spécifique pour plus de détails.

On peut adopter deux positions extrêmes sur les opérations arithmétiques entre
éléments d'anneaux *différents* :

* les anneaux différents sont des mondes indépendants, et l'addition ou la
  multiplication entre éléments d'anneaux différents n'ont aucun sens ; même
  ``1 + 1/2`` n'a pas de sens puisque le premier terme est un entier et le
  second un rationnel ;

ou

* si un élément ``r1`` d'un anneau ``R1`` peut, d'une manière ou d'une autre,
  s'interpréter comme élément d'un autre anneau ``R2``, alors toutes les
  opérations arithmétiques entre ``r1`` et un élément quelconque de ``R2`` sont
  permises. En particulier, les éléments neutres de la multiplication dans les
  corps et anneaux doivent tous être égaux entre eux.

Sage adopte un compromis. Si ``P1`` et ``P2`` sont des parents et si ``p1`` est
un élément de ``P1``, l'utilisateur peut demander explicitement comment ``P1``
s'interprète dans ``P2``. Cela n'a pas forcément de sens dans tous les cas, et
l'interprétation peut n'être définie que pour certains éléments de ``P1`` ;
c'est à l'utilisateur de s'assurer que la conversion a un sens. Cela s'appelle
une **conversion** :

::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: GF(5)(a) == b
    True
    sage: GF(2)(b) == a
    True

Cependant, une conversion *implicite* (c'est-à-dire automatique) n'est
possible que si elle peut se faire *systématiquement* et de manière
*cohérente*. Il faut ici absolument faire preuve de rigueur.

Une telle conversion implicite s'appelle une **coercition**. Si une coercition
est définie entre deux parents, elle doit coïncider avec la conversion. De
plus, les coercitions doivent obéir aux deux conditions suivantes :

#. Une coercition de ``P1`` dans ``P2`` doit être un morphisme (par exemple
   un morphisme d'anneaux). Elle doit être définie pour *tous* les éléments de
   ``P1``, et préserver la structure algébrique de celui-ci.
#. Le choix des applications de coercition doit être fait de manière cohérente.
   Si ``P3`` est un troisième parent, la composée de la coercition choisie
   de ``P1`` dans ``P2`` et de celle de ``P2`` dans ``P3`` doit être la
   coercition de ``P1`` dans ``P3``. En particulier, s'il existe des
   coercitions de ``P1`` dans ``P2`` et de ``P2`` dans ``P1``, leur composée
   doit être l'identité sur ``P1``.

Ainsi, bien qu'il soit possible de convertir tout élément de ``GF(2)`` en un
élément de ``GF(5)``, la conversion ne peut être une coercition, puisque il
n'existe pas de morphisme d'anneaux de ``GF(2)`` dans ``GF(5)``.

Le second point — la cohérence des choix — est un peu plus compliqué à
expliquer. Illustrons-le sur l'exemple des anneaux de polynômes multivariés.
Dans les applications, il s'avère utile que les coercitions respectent les noms
des variables. Nous avons donc :

::

    sage: R1.<x,y> = ZZ[]
    sage: R2 = ZZ['y','x']
    sage: R2.has_coerce_map_from(R1)
    True
    sage: R2(x)
    x
    sage: R2(y)
    y
    sage: R2.coerce(y)
    y

En l'absence d'un morphisme d'anneau qui préserve les noms de variable, la
coercition entre anneaux de polynômes multivariés n'est pas définie. Il peut
tout de même exister une conversion qui envoie les variables d'un anneau sur
celle de l'autre en fonction de leur position dans la liste des générateurs :

.. link

::

    sage: R3 = ZZ['z','x']
    sage: R3.has_coerce_map_from(R1)
    False
    sage: R3(x)
    z
    sage: R3(y)
    x
    sage: R3.coerce(y)
    Traceback (most recent call last):
    ...
    TypeError: no canonical coercion
    from Multivariate Polynomial Ring in x, y over Integer Ring
    to Multivariate Polynomial Ring in z, x over Integer Ring

Mais une telle conversion ne répond pas aux critères pour être une coercition :
en effet, en composant l'application de ``ZZ['x','y']`` dans ``ZZ['y','x']``
avec celle qui préserve les positions de ``ZZ['y','x']`` dans ``ZZ['a','b']``,
nous obtiendrions une application qui ne préserve ni les noms ni les positions,
ce qui viole la règle de cohérence.

Lorsqu'une coercition est définie, elle est souvent utilisée pour comparer des
éléments d'anneaux différents ou pour effectuer des opérations arithmétiques.
Cela est commode, mais il faut être prudent en étendant la relation d'égalité
``==`` au-delà des frontières d'un parent donné. Par exemple, si ``==`` est
bien censé être une relation d'équivalence entre éléments d'*un* anneau, il
n'en va pas forcément de même quand on compare des éléments d'anneaux
différents. Ainsi, les éléments ``1`` de ``ZZ`` et d'un corps fini sont
considérés comme égaux, puisqu'il existe une coercition canonique des entiers
dans tout corps fini. En revanche, il n'y a en général pas de coercition entre
deux corps finis quelconques. On a donc

.. link

::

    sage: GF(5)(1) == 1
    True
    sage: 1 == GF(2)(1)
    True
    sage: GF(5)(1) == GF(2)(1)
    False
    sage: GF(5)(1) != GF(2)(1)
    True

De même, on a

.. link

::

    sage: R3(R1.1) == R3.1
    True
    sage: R1.1 == R3.1
    False
    sage: R1.1 != R3.1
    True

Une autre conséquence de la condition de cohérence est que les coercitions ne
sont possibles que des anneaux exacts (comme les rationnels ``QQ``) vers les
anneaux inexacts (comme les réels à précision donnée ``RR``), jamais l'inverse.
En effet, pour qu'une conversion de ``RR`` dans ``QQ`` puisse être une
coercition, il faudrait que la composée de la coercition de ``QQ`` dans ``RR``
et de cette conversion soit l'identité sur ``QQ``, ce qui n'est pas possible
puisque des rationnels distincts peuvent très bien être envoyés sur le même
élément de ``RR`` :

::

    sage: RR(1/10^200+1/10^100) == RR(1/10^100)
    True
    sage: 1/10^200+1/10^100 == 1/10^100
    False

Lorsque l'on compare des éléments de deux parents ``P1`` et ``P2``, il peut
arriver qu'il n'existe pas de coercition entre ``P1`` et ``P2``, mais qu'il y
ait un choix canonique de parent ``P3`` tel que ``P1`` et ``P2`` admettent tous
deux des coercitions dans ``P3``. Dans ce cas aussi, la coercition a lieu. Un
exemple typique de ce mécanisme est l'addition d'un rationnel et d'un polynôme
à coefficients entiers, qui produit un polynôme à coefficients rationnels :

::

    sage: P1.<x> = ZZ[]
    sage: p = 2*x+3
    sage: q = 1/2
    sage: parent(p)
    Univariate Polynomial Ring in x over Integer Ring
    sage: parent(p+q)
    Univariate Polynomial Ring in x over Rational Field

Notons qu'en principe, on aurait très bien pu choisir pour ``P3`` le corps des
fractions de ``ZZ['x']``. Cependant, Sage tente de choisir un parent commun
*canonique* aussi naturel que possible (ici ``QQ['x']``). Afin que cela
fonctionne de façon fiable, Sage ne se contente *pas* de prendre n'importe
lequel lorsque plusieurs candidats semblent aussi naturels les uns que les
autres. La manière dont le choix est fait est décrite dans le `tutoriel
<http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_
spécifique déjà mentionné.

Dans l'exemple suivant, il n'y a pas de coercition vers un parent commun :

::

    sage: R.<x> = QQ[]
    sage: S.<y> = QQ[]
    sage: x+y
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +: 'Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

En effet, Sage refuse de choisir entre les candidats  ``QQ['x']['y']``,
``QQ['y']['x']``, ``QQ['x','y']`` et ``QQ['y','x']``, car ces quatre structures
deux à deux distinctes semblent toutes des parents communs naturels, et aucun
choix canonique ne s'impose.
