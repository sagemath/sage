.. _section-rings:

***************
Anneaux de base
***************

Nous illustrons la prise en main de quelques anneaux de base avec Sage.
Par exemple, ``RationalField()`` ou ``QQ`` désigneront dans ce qui
suit au corps des nombres rationnels :

::

    sage: RationalField()
    Rational Field
    sage: QQ
    Rational Field
    sage: 1/2 in QQ
    True

Le nombre décimal ``1.2`` est considéré comme un élément de ``QQ``,
puisqu'il existe une application de coercition entre les réels et les
rationnels :

::

    sage: 1.2 in QQ
    True

Néanmoins, il n'y a pas d'application de coercition entre le corps fini
à 3 éléments et les rationnels :

::

    sage: c = GF(3)(1)   # c est l'élément 1 du corps fini à 3 éléments
    sage: c in QQ
    False

De même, bien entendu, la constante symbolique :math:`\pi` n'appartient
pas aux rationnels :

::

    sage: pi in QQ
    False

Le symbole ``I`` représente la racine carrée de :math:`-1`; ``i`` est
synonyme de ``I``. Bien entendu, ``I`` n'appartient pas aux rationnels :

::

    sage: i  # i^2 = -1
    I
    sage: i in QQ
    False

À ce propos, d'autres anneaux sont prédéfinis en Sage : l'anneau des
entiers relatifs ``ZZ``, celui des nombres réels ``RR`` et celui des
nombres complexes ``CC``. Les anneaux de polynômes sont décrits dans
:ref:`section-poly`.

Passons maintenant à quelques éléments d'arithmétique.

::

    sage: a, b = 4/3, 2/3
    sage: a + b
    2
    sage: 2*b == a
    True
    sage: parent(2/3)
    Rational Field
    sage: parent(4/2)
    Rational Field
    sage: 2/3 + 0.1       # coercition automatique avant addition
    0.766666666666667
    sage: 0.1 + 2/3       # les règles de coercition sont symétriques en Sage
    0.766666666666667

Il y a une subtilité dans la définition des nombres complexes. Comme
mentionné ci-dessus, le symbole  ``i`` représente une racine carrée de
:math:`-1`, mais il s'agit d'une racine carrée *formelle* de :math:`-1`,
comme élément d'un corps de nombres quadratique. L'appel ``CC(i)`` renvoie
la racine carrée de :math:`-1` comme nombre complexe en virgule flottante.

.. link

::

    sage: i = CC(i)       # nombre complexe en virgule flottante
    sage: z = a + b*i
    sage: z
    1.33333333333333 + 0.666666666666667*I
    sage: z.imag()        # partie imaginaire
    0.666666666666667
    sage: z.real() == a   # coercition automatique avant comparaison
    True
    sage: QQ(11.1)
    111/10
