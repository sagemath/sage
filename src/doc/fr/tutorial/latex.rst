*********************************
Sage, LaTeX et compagnie
*********************************

AUTEUR:  Rob Beezer (2010-05-23)

Sage et le dialecte LaTeX de TeX entretiennent une forte synergie. L'objet de
cette section est de présenter leurs différentes interactions. Nous commençons
par les plus basiques, avant de passer à des fonctionnalités plus obscures.
(Une partie de cette section peut donc être sautée en première lecture.)

Vue d'ensemble
==============

Le plus simple pour comprendre comment Sage peut faire appel à LaTeX est
peut-être de passer en revue les trois principales techniques.

    #. Dans Sage, chaque « objet » doit disposer d'une représentation LaTeX. La
       représentation d'un objet ``foo`` est accessible par la commande
       ``latex(foo)``, utilisable dans le bloc-notes ou en ligne de commande.
       Celle-ci renvoie une chaîne de caractères qui, interprétée par TeX en
       mode mathématique (entre signes dollar simples par exemple) devrait
       produire une représentation raisonnable de l'objet ``foo``. Nous verrons
       quelques exemples plus bas.

       On peut ainsi utiliser Sage pour préparer des fragments de document
       LaTeX : il suffit de définir ou obtenir par un calcul un objet Sage et
       de copier-coller le résultat de ``latex()`` appliqué à l'objet dans le
       document en question.

    #. Le bloc-notes de Sage fait par défaut appel à
       `MathJax <http://www.mathjax.org>`_
       pour afficher proprement les formules mathématiques dans le navigateur
       web. MathJax est un moteur de rendu mathématique open source écrit en
       JavaScript compatible avec tous les navigateurs récents. MathJax
       interprète un sous-ensemble conséquent de TeX, mais pas la totalité.
       Orienté avant tout vers le rendu correct de petits fragments de TeX, il
       ne gère par exemple ni les tableaux compliqués ni le découpage en
       sections du document. Le rendu automatique des formules mathématiques
       dans le bloc-notes fonctionne en convertissant la sortie de ``latex()``
       mentionnée ci-dessus en une portion de document HTML que MathJax sait
       interpréter.

       MathJax utilise ses propres polices de caractères vectorielles, et
       fournit ainsi un rendu de meilleure qualité que les méthodes d'affichage
       d'équations ou d'autres fragments de TeX qui passent par des images
       bitmap statiques.

    #. Il est possible de faire appel à une installation extérieure de LaTeX
       depuis la ligne de commande de Sage, ou depuis le bloc-notes pour
       interpréter du code plus compliqué que ce que MathJax sait traiter. La
       distribution Sage inclut pratiquement tout le nécessaire pour compiler
       et utiliser le logiciel Sage, à la notable exception de TeX. Il faut
       donc installer par ailleurs TeX et quelques utilitaires de conversion
       associés pour l'utiliser depuis Sage de cette manière.

Voici quelques exemples d'utilisation élémentaire de la fonction ``latex()``. ::

    sage: var('z')
    z
    sage: latex(z^12)
    z^{12}
    sage: latex(integrate(z^4, z))
    \frac{1}{5} \, z^{5}
    sage: latex('a string')
    \text{\texttt{a{ }string}}
    sage: latex(QQ)
    \Bold{Q}
    sage: latex(matrix(QQ, 2, 3, [[2,4,6],[-1,-1,-1]]))
    \left(\begin{array}{rrr}
    2 & 4 & 6 \\
    -1 & -1 & -1
    \end{array}\right)

L'utilisation de base de MathJax dans le bloc-notes est largement automatique.
Nous pouvons tout de même en voir quelques exemples en employant la classe
``MathJax``. La méthode ``eval`` de celle-ci convertit un objet Sage en sa
représentation LaTeX, puis emballe le résultat dans du code HTML qui fait
possède la classe CSS "math", laquelle indique de faire appel à MathJax. ::

    sage: from sage.misc.html import MathJax
    sage: mj = MathJax()
    sage: var('z')
    z
    sage: mj(z^12)
    <html>\[z^{12}\]</html>
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>
    sage: mj(ZZ['x'])
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}[x]\]</html>
    sage: mj(integrate(z^4, z))
    <html>\[\frac{1}{5} \, z^{5}\]</html>

Utilisation de base
===================

Comme indiqué dans la vue d'ensemble, la manière la plus simple d'exploiter le
support LaTeX de Sage consiste à appeler la fonction ``latex()`` pour produire
du code LaTeX représentant un objet mathématique. Les chaînes obtenues peuvent
ensuite être incorporées dans des documents LaTeX indépendants. Cela fonctionne
de la même manière dans le bloc-notes et en ligne de commande.

L'autre extrême est la commande ``view()``, qui fait tout le nécessaire pour
afficher le rendu correspondant au code LaTeX.

En ligne de commande, ``view(foo)`` produit la représentation LaTeX de ``foo``,
la place dans un document LaTeX simple, et compile ce document en utilisant
l'installation de TeX du système. Elle appelle ensuite un autre programme
externe pour afficher le résultat de la compilation. La version de TeX (et donc
le format de sortie) ainsi que la visionneuse à utiliser sont configurables,
voir :ref:`sec-custom-processing`.

Dans le bloc-notes, ``view(foo)`` produit une combinaison de HTML et CSS qui
indique à MathJax de s'occuper du rendu de la représentation LaTeX.
L'utilisateur voit une version joliment formatée de la sortie à la place de la
sortie ASCII par défaut de Sage. Certains objets ont cependant des
représentations LaTeX trop compliquées pour être affichés par MathJax. Lorsque
c'est le cas, il est possible de contourner l'interprétation par MathJax,
d'appeler l'installation LaTeX du système, et de convertir le document produit
en image pour l'afficher dans le bloc-note. La section
:ref:`sec-custom-generation` ci-dessous explique comment configurer et
contrôler ce processus.

La commande interne ``pretty_print()`` permet de convertir un objet Sage en code
HTML utilisant MathJax. C'est le code qui sera ensuite utilisé dans le
bloc-notes ::

    sage: pretty_print(x^12)
    x^12
    sage: pretty_print(integrate(sin(x), x))
    -cos(x)

Le bloc-notes dispose de deux autres fonctionnalités pour appeler LaTeX.
Premièrement, lorsque la case « Typeset » (juste au-dessus de la première
cellule d'une feuille de travail, à droite des quatre listes déroulantes) est
cochée, le résultat de l'évaluation d'une cellule est automatiquement
interprété par MathJax et affiché sous forme de formule plutôt que de texte
brut. Les sorties déjà affichées ne sont pas modifiées tant que l'on ne
ré-évalue pas les cellules correspondantes. Cocher la case « Typeset » revient
essentiellement à appeler ``view()`` sur le résultat de chaque cellule.

Deuxièmement, le bloc-notes permet d'annoter une feuille de travail en
saisissant du TeX. Un clic en tenant la touche Maj enfoncée sur la barre bleue
qui apparaît lorsque l'on place le curseur de la souris entre deux cellules
ouvre un mini-traitement de texte appelé TinyMCE. Cela permet de saisir du
texte pour commenter la feuille de travail, et de le mettre en forme avec un
éditeur WYSIWIG de HTML et CSS. Mais le texte placé entre signes dollar simples
ou doubles est interprété par MathJax, respectivement comme formule composée en
ligne ou hors texte.

.. _sec-custom-generation:

Personnaliser le code LaTeX produit
===================================

Les méthodes de l'objet prédéfini ``latex`` permettent de personnaliser le code
LaTeX produit par la commande ``latex()`` de différentes manières. Cela
s'applique dans le bloc-notes comme en ligne de commande. On obtient la liste
des méthodes en saisissant ``latex.`` (noter la présence du point) puis en
appuyant sur la touche tabulation.

Un bon exemple est la méthode ``latex.matrix_delimiters``, qui sert à modifier
les symboles entourant les matrices -- parenthèses, crochets, accolades
ou barres verticales par exemple. Les délimiteurs gauche et droit sont
donnés par des chaînes LaTeX. Ils n'ont pas besoin de se correspondre. Notons
comment les contre-obliques qui doivent être interprétées par TeX sont
protégées par une seconde contre-oblique dans l'exemple ci-dessous. ::

    sage: A = matrix(ZZ, 2, 2, range(4))
    sage: latex(A)
    \left(\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right)
    sage: latex.matrix_delimiters(left='[', right=']')
    sage: latex(A)
    \left[\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right]
    sage: latex.matrix_delimiters(left='\\{', right='\\}')
    sage: latex(A)
    \left\{\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right\}

La méthode ``latex.vector_delimiters`` fonctionne de manière analogue.

Les anneaux et corps usuels (entiers, rationnels, réels, etc.) sont par défaut
composés en gras. La méthode ``latex.blackboard_bold`` permet de changer pour
des lettres ajourées. Elle ne change pas la sortie de la commande ``latex()``
mais la définition de la macro TeX ``\Bold{}`` fournie par Sage. ::

    sage: latex(QQ)
    \Bold{Q}
    sage: from sage.misc.html import MathJax
    sage: mj=MathJax()
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>
    sage: latex.blackboard_bold(True)
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbb{#1}}\Bold{Q}\]</html>
    sage: latex.blackboard_bold(False)

On peut aussi définir de nouvelles macros TeX ou charger des packages
supplémentaires. L'exemple suivant montre comment ajouter des macros qui seront
utilisées à chaque fois que MathJax interprète un fragment de TeX dans le
bloc-notes. ::

    sage: latex.extra_macros()
    ''
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: latex.extra_macros()
    '\\newcommand{\\foo}{bar}'
    sage: var('x y')
    (x, y)
    sage: latex(x+y)
    x + y
    sage: from sage.misc.html import MathJax
    sage: mj=MathJax()
    sage: mj(x+y)
    <html>\[\newcommand{\foo}{bar}x + y\]</html>

Ces macros supplémentaires sont disponibles aussi quand Sage appelle TeX pour
compiler un fragment de document trop gros pour MathJax. C'est la fonction
``latex_extra_preamble``, appelée pour préparer le préambule du document LaTeX,
qui les définit, comme l'illustre l'exemple suivant. Notons à nouveau le
dédoublement des ``\`` dans les chaînes Python. ::

    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: from sage.misc.latex import latex_extra_preamble
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    \newcommand{\foo}{bar}

On peut aussi charger des packages LaTeX, ou ajouter n'importe quelle autre
commande au préambule, grâce à la méthode ``latex.add_package_to_preamble``. Sa
variante plus spécialisée ``latex.add_package_to_preamble_if_available``
vérifie qu'un package donné est disponible avant de l'ajouter au préambule si
c'est bien le cas.

Dans l'exemple suivant, nous ajoutons au préambule la commande qui charge le
package ``geometry``, puis nous l'utilisons pour régler la taille de la zone de
texte (et donc indirectement les marges) du document TeX. Une fois encore, les
contre-obliques sont dédoublées. ::

    sage: from sage.misc.latex import latex_extra_preamble
    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: latex.add_to_preamble('\\usepackage{geometry}')
    sage: latex.add_to_preamble('\\geometry{letterpaper,total={8in,10in}}')
    sage: latex.extra_preamble()
    '\\usepackage{geometry}\\geometry{letterpaper,total={8in,10in}}'
    sage: print(latex_extra_preamble())
    \usepackage{geometry}\geometry{letterpaper,total={8in,10in}}
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}

Voici enfin comment ajouter un package en vérifiant sa disponibilité, et ce
qu'il se passe quand le package n'existe pas. ::

    sage: latex.extra_preamble('')
    sage: latex.extra_preamble()
    ''
    sage: latex.add_to_preamble('\\usepackage{foo-bar-unchecked}')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'
    sage: latex.add_package_to_preamble_if_available('foo-bar-checked')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'

.. _sec-custom-processing:

Personnaliser le traitement du code par LaTeX
=============================================

En plus de modifier LaTeX produit par Sage, on peut choisir la variante de
TeX appelée pour le traiter, et donc la nature du document produit. De même, il
est possible de contrôler dans quelles circonstances le bloc-notes utilisera
MathJax (c'est-à-dire quels fragments de code TeX seront jugés suffisamment
simples) et quand il choisira de se rabattre sur l'installation de TeX du système.

La méthode ``latex.engine()`` permet de choisir lequel des moteurs TeX
``latex``, ``pdflatex`` et ``xelatex`` doit servir à compiler les expressions
LaTeX complexes. Lorsque l'on appelle ``view`` en ligne de commande, si le
moteur actif est ``latex``, celui-ci produit un fichier dvi, puis Sage fait
appel à une visionneuse dvi (par exemple xdvi) pour afficher le résultat. Si en
revanche le moteur est ``pdflatex``, il produit par défaut un fichier PDF, que
Sage affiche grâce à la visionneuse PDF du système (Adobe Reader, Okular,
evince...).

Dans le bloc-notes, la première étape est de décider s'il faut utiliser MathJax
ou LaTeX pour interpréter un fragment de TeX donné. La décision se fonde sur
une liste de chaînes « interdites » dont la présence dans le fragment indique
d'appeler latex (ou plus généralement le moteur choisi via ``latex.engine()``)
au lieu MathJax. Les méthodes ``latex.add_to_mathjax_avoid_list`` et
``latex.mathjax_avoid_list`` permettent de gérer le contenu de cette liste. ::

    sage: # not tested
    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []
    sage: latex.mathjax_avoid_list(['foo', 'bar'])
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar']
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar', 'tikzpicture']
    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []

Supposons maintenant que, dans le bloc-notes, un appel à ``view()`` ou
l'évaluation d'une cellule lorsque la case "Typeset" est cochée produise un
résultat dont le mécanisme décrit ci-dessus détermine qu'il doit être passé au
moteur LaTeX externe. Comme en ligne de commande, l'exécutable spécifié par
``latex.engine()`` traite alors le document. Cependant, au lieu d'appeler une
visionneuse externe pour afficher le document produit, Sage tente de recadrer
le document en rognant les zones blanches, et de le convertir en une image qui
est ensuite insérée dans le bloc-notes comme sortie associée à la cellule.

Plusieurs facteurs influencent la conversion, principalement le moteur TeX
choisi et la palette d'utilitaires de conversion disponibles sur le système.
Les convertisseurs suivants couvrent à eux quatre toutes les situations :
``dvips``, ``ps2pdf``, ``dvipng`` et ``convert`` (de la collection ImageMagick).
Dans tous les cas, il s'agit d'arriver à produire une image PNG à insérer dans
la feuille de travail. Lorsque le moteur LaTeX produit un fichier dvi, le
programme dvipng suffit en général à effectuer la conversion. Il peut cependant
arriver que le fichier dvi contienne des instructions spécifiques à un pilote
(commande TeX ``\special``) que dvipng est incapable d'interpréter, auquel cas
on utilise ``dvips`` pour créer un fichier PostScript. Celui-ci, de même que
le fichier PDF produit par le moteur le cas échéant, est ensuite converti en
image png avec ``convert``. Les commandes ``have_dvipng()`` et
``have_convert()`` permettent de tester la présence sur le système des
utilitaires en question.

Toutes ces conversions sont automatiques lorsque les outils nécessaires sont
installés. Dans le cas contraire, un message d'erreur indique ce qu'il manque
et où le télécharger.

La section suivante (:ref:`sec-tkz-graph`) présente un exemple concret de
traitement d'expressions LaTeX complexes, en l'occurrence pour obtenir un rendu
de qualité de graphes grâce au package LaTeX ``tkz-graph``. Les tests inclus
dans Sage contiennent d'autres exemples. On y accède en important l'objet
prédéfini ``sage.misc.latex.latex_examples``, instance de la classe
``sage.misc.latex.LatexExamples``, comme illustré ci-dessous. Les exemples
fournis actuellement couvrent les types d'objets suivants : diagrammes
commutatifs (utilisant le package `xy`), graphes combinatoires (`tkz-graph`),
nœuds (`xypic`), schémas pstricks (`pstricks`). Pour obtenir la liste des
exemples, utilisez la complétion de ligne de commande après avoir importé
``latex_examples``. Chaque exemple affiche quand on l'appelle des instructions
sur la configuration nécessaire pour qu'il fonctionne correctement. Une fois le
préambule, le moteur LaTeX etc. configurés comme indiqué, il suffit d'appeler
la commande ``view()`` pour visualiser l'exemple. ::

    sage: from sage.misc.latex import latex_examples
    sage: latex_examples.diagram()
    LaTeX example for testing display of a commutative diagram produced
    by xypic.
    <BLANKLINE>
    To use, try to view this object -- it will not work.  Now try
    'latex.add_to_preamble("\\usepackage[matrix,arrow,curve,cmtip]{xy}")',
    and try viewing again. You should get a picture (a part of the diagram arising
    from a filtered chain complex).

.. _sec-tkz-graph:

Example : rendu de graphes avec tkz-graph
=========================================

Le package ``tkz-graph`` permet de produire des dessins de graphes
(combinatoires) de qualité. Il repose sur TikZ, lui-même une interface pour la
bibliothèque TeX pgf : pgf, TikZ et tkz-graph doivent donc tous être présents
dans l'installation TeX du système pour que cet exemple fonctionne. Les
versions fournies par certaines distributions TeX sont parfois trop anciennes,
et il peut donc être souhaitable de les installer manuellement dans son arbre
texmf personnel. On consultera la documentation de la distribution TeX pour la
procédure à suivre, qui dépasse le cadre de ce document. La section
:ref:`sec-system-wide-tex` donne la liste des fichiers nécessaires.

Il nous faut tout d'abord nous assurer que les packages requis sont inclus dans
le document LaTeX, en les ajoutant au préambule. Le rendu des graphes n'est pas
correct quand on passe par le format dvi, aussi il est préférable de
sélectionner ``pdflatex`` comme moteur TeX. Après ces réglages, une instruction
du genre ``view(graphs.CompleteGraph(4))`` saisie dans l'interface en ligne de
commande doit produire un fichier PDF contenant un dessin du graphe complet
`K_4`.

Pour que la même chose fonctionne dans le bloc-notes, il faut de plus
désactiver l'interprétation du code LaTeX produisant le graphe par MathJax, à
l'aide de la liste de motifs exclus. Le nom de l'environnement ``tikzpicture``,
dans lequel sont placés les graphes, est un bon choix de chaîne à exclure. Une
fois cela fait, la commande ``view(graphs.CompleteGraph(4))`` dans une feuille
de travail du bloc-notes appelle pdflatex pour produire un fichier PDF, puis
``convert`` pour en extraire une image PNG à placer dans la zone de sortie de
la feuille de travail. Les commandes suivantes reprennent l'ensemble des
étapes de configuration. ::

    sage: from sage.graphs.graph_latex import setup_latex_preamble
    sage: setup_latex_preamble()
    sage: latex.extra_preamble() # random - depends on system's TeX installation
    '\\usepackage{tikz}\n\\usepackage{tkz-graph}\n\\usepackage{tkz-berge}\n'
    sage: latex.engine('pdflatex')
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')  # not tested
    sage: latex.mathjax_avoid_list()                      # not tested
    ['tikz', 'tikzpicture']

La mise en forme du graphe est faite en traitant des commandes ``tkz-graph``
qui le décrivent avec ``pdflatex``. Diverses options pour influencer ce rendu,
qui sortent du cadre de cette section, sont décrites dans la section intitulée
"LaTeX Options for Graphs" du manuel de référence de Sage.

.. _sec-system-wide-tex:

Une installation TeX pleinement opérationnelle
==============================================

Beaucoup de fonctionnalités avancées de l'intégration Sage-TeX nécessitent
qu'une installation extérieure de TeX soit disponible sur le système. Les
distributions Linux en fournissent généralement, sous forme de paquets basés
sur TeX Live ; sous OS X, on peut installer TeXshop ; et sous Windows, MikTeX.
L'utilitaire ``convert`` fait partie de la boîte à outils `ImageMagick
<http://www.imagemagick.org/>`_ (probablement disponible dans l'archive de
paquets de votre système ou facile à télécharger et installer). Les programmes
``dvipng``, ``ps2pdf``, and ``dvips`` sont parfois inclus dans les
installations de TeX, et les deux premiers sont par ailleurs disponibles
respectivement à l'adresse http://sourceforge.net/projects/dvipng/ et dans
`Ghostscript <http://www.ghostscript.com/>`_.

Le rendu des graphes nécessite une version suffisamment récente de PGF, ainsi
que les fichiers ``tkz-graph.sty``, disponible sur le site web
https://www.ctan.org/pkg/tkz-graph, ``tkz-arith.sty`` et suivant les
cas ``tkz-berge.sty``, disponibles sur le site web
https://www.ctan.org/pkg/tkz-berge.

Programmes externes
===================

Trois programmes séparés contribuent encore à l'intégration TeX-Sage.

Le premier, sagetex, est (pour simplifier) une collection de macros TeX qui
permettent d'introduire dans un document LaTeX des instructions qui seront
interprétées par Sage pour effectuer des calculs et/ou mettre en forme des
objets mathématiques avec la commande ``latex()`` de Sage. Il est donc possible
de faire faire des calculs à Sage et de produire les sorties LaTeX associées
comme étape intermédiaire de la compilation d'un document LaTeX. Par exemple,
on peut imaginer de maintenir la correspondance entre questions et réponses
dans un sujet d'examen en utilisant Sage pour calculer les unes à partir des
autres. Sagetex est décrit plus en détail en section :ref:`sec-sagetex` de ce
document.
