.. sage-doctest: needs sage.plot sage.symbolic

.. _section-plot:

Graphiques
==========

Sage peut produire des graphiques en deux ou trois dimensions.

Graphiques en deux dimensions
-----------------------------

En deux dimensions, Sage est capable de tracer des cercles, des droites,
des polygones, des graphes de fonctions en coordonnées cartésiennes, des
graphes en coordonnées polaires, des lignes de niveau et des
représentations de champs de vecteurs. Nous présentons quelques exemples
de ces objets ici. Pour plus d'exemples de graphiques avec Sage, on
consultera :ref:`section-systems`, :ref:`section-maxima` et aussi la
documentation
`Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_

La commande suivante produit un cercle jaune de rayon 1 centré à l'origine :

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0))
    Graphics object consisting of 1 graphics primitive

Il est également possible de produire un disque plein :

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0), fill=True)
    Graphics object consisting of 1 graphics primitive

Il est aussi possible de créer un cercle en l'affectant à une variable ;
ceci ne provoque pas son affichage.

::

    sage: c = circle((0,0), 1, rgbcolor=(1,1,0))

Pour l'afficher, on utilise ``c.show()``  ou  ``show(c)``, comme suit :

.. link

::

    sage: c.show()

Alternativement, l'évaluation de ``c.save('filename.png')`` enregistre
le graphique dans le fichier spécifié.

Toutefois, ces « cercles » ressemblent plus à des ellipses qu'à des
cercles, puisque les axes possèdent des échelles différentes. On peut
arranger ceci :

.. link

::

    sage: c.show(aspect_ratio=1)

La commande ``show(c, aspect_ratio=1)`` produit le même résultat. On
peut enregistrer l'image avec cette option par la commande
``c.save('filename.png', aspect_ratio=1)``.

Il est très facile de tracer le graphique de fonctions de base :

::

    sage: plot(cos, (-5,5))
    Graphics object consisting of 1 graphics primitive

En spécifiant un nom de variable, on peut aussi créer des graphes
paramétriques :

::

    sage: x = var('x')
    sage: parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    Graphics object consisting of 1 graphics primitive

Différents graphiques peuvent se combiner sur une même image :

::

    sage: x = var('x')
    sage: p1 = parametric_plot((cos(x),sin(x)),(x,0,2*pi),rgbcolor=hue(0.2))
    sage: p2 = parametric_plot((cos(x),sin(x)^2),(x,0,2*pi),rgbcolor=hue(0.4))
    sage: p3 = parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    sage: show(p1+p2+p3, axes=false)

Une manière commode de tracer des formes pleines est de préparer une
liste de points  (``L`` dans l'exemple ci-dessous) puis d'utiliser la
commande ``polygon`` pour tracer la forme pleine dont le bord est formé
par ces points. Par example, voici un deltoïde vert :

::

    sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),
    ....:     2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
    sage: polygon(L, rgbcolor=(1/8,3/4,1/2))
    Graphics object consisting of 1 graphics primitive

Pour visualiser le graphique en masquant les axes, tapez ``show(p,
axes=false)``.

On peut ajouter un texte à un graphique :

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),
    ....:     6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t)

En cours d'analyse, les professeurs font souvent le dessin suivant
au tableau : non pas une mais plusieurs branches de la fonction arcsin, autrement dit, le graphe d'équation :math:`y=\sin(x)`
pour  :math:`x` entre :math:`-2\pi` et :math:`2\pi`, renversé par
symétrie par rapport à la première bissectrice des axes. La commande
Sage suivante réalise cela :

::

    sage: v = [(sin(x),x) for x in srange(-2*float(pi),2*float(pi),0.1)]
    sage: line(v)
    Graphics object consisting of 1 graphics primitive

Comme les valeurs prises par la fonction tangente ne sont pas bornées,
pour utiliser la même astuce pour représenter la fonction arctangente,
il faut préciser les bornes de la coordonnée *x* :

::

    sage: v = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
    sage: show(line(v), xmin=-20, xmax=20)

Sage sait aussi tracer des graphiques en coordonnées polaires, des
lignes de niveau et (pour certains types de fonctions) des champs de
vecteurs.  Voici un exemple de lignes de niveau :

::

    sage: f = lambda x,y: cos(x*y)
    sage: contour_plot(f, (-4, 4), (-4, 4))
    Graphics object consisting of 1 graphics primitive

Graphiques en trois dimensions
------------------------------

Sage produit des graphes en trois dimensions en utilisant le package
open source appelé [ThreeJS]_. En voici quelques exemples :

Le parapluie de Whitney tracé en jaune
http://en.wikipedia.org/wiki/Whitney_umbrella:

::

    sage: u, v = var('u,v')
    sage: fx = u*v
    sage: fy = u
    sage: fz = v^2
    sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1),
    ....:   frame=False, color="yellow")
    Graphics3d Object


Une fois évaluée la commande ``parametric_plot3d``, qui affiche le
graphique,  il est possible de cliquer et de le tirer pour
faire pivoter la figure.

Le bonnet croisé (cf. http://en.wikipedia.org/wiki/Cross-cap ou
http://www.mathcurve.com/surfaces/bonnetcroise/bonnetcroise.shtml) :

::

    sage: u, v = var('u,v')
    sage: fx = (1+cos(v))*cos(u)
    sage: fy = (1+cos(v))*sin(u)
    sage: fz = -tanh((2/3)*(u-pi))*sin(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object

Un tore tordu :

::

    sage: u, v = var('u,v')
    sage: fx = (3+sin(v)+cos(u))*cos(2*v)
    sage: fy = (3+sin(v)+cos(u))*sin(2*v)
    sage: fz = sin(u)+2*cos(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object


