.. _chapter-interactive_shell:

********************************
La ligne de commande interactive
********************************

Dans la plus grande partie de ce tutoriel, nous supposons que vous avez
lancé l'interpréteur Sage avec la commande ``sage``. Cela démarre une
version adaptée du shell (interpréteur de commandes) IPython et importe
un grand nombre de fonctions et de classes qui sont ainsi prêtes à
l'emploi depuis l'invite de commande. D'autres personnalisations sont
possibles en éditant le fichier ``$SAGE_ROOT/ipythonrc``. Au démarrage,
le shell Sage affiche un message de ce genre :

.. skip

::

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘


    sage:

Pour quitter Sage, tapez Ctrl-D ou tapez ``quit`` ou ``exit``.

.. skip

::

    sage: quit
    Exiting Sage (CPU time 0m0.00s, Wall time 0m0.89s)

L'indication *wall time* donne le temps écoulé à votre montre (ou
l'horloge suspendue au mur) pendant l'exécution. C'est une donnée
pertinente car le temps processeur (*CPU time*) ne tient pas compte du
temps utilisé par les sous-processus comme GAP et Singular.

(Il vaut mieux éviter de tuer un processus Sage depuis un terminal avec
``kill -9``, car il est possible que Sage ne tue pas ses processus
enfants, par exemple des processus Maple qu'il aurait lancés, ou encore
qu'il ne nettoie pas les fichiers temporaires de
``$HOME/.sage/tmp``.)

Votre session Sage
==================

Une session est la suite des entrées et sorties qui interviennent entre
le moment où vous démarrez Sage et celui où vous le quittez. Sage
enregistre un journal de toutes les entrées via IPython. Si vous
utilisez le shell interactif (par opposition à l'interface *jupyter*),
vous pouvez taper ``%history`` (ou ``%hist``) à n'importe
quel moment pour obtenir la
liste de toutes les lignes de commandes entrées depuis le début de la
session. Tapez ``?`` à l'invite de commande Sage pour plus
d'informations sur IPython. Par exemple : « IPython fournit des invites
de commande numérotées [...] avec un cache des entrées-sorties. Toutes
les entrées sont sauvegardées et peuvent être rappelées comme variables
(en plus de la navigation habituelle dans l'historique avec les flèches
du clavier). Les variables GLOBALES suivantes existent toujours (ne les
écrasez pas !) » :

::

      _ : dernière entrée
      __ : avant-dernière entrée
      _oh : liste de toutes les entrées précédentes

Voici un exemple :

.. skip

::

    sage: factor(100)
     _1 = 2^2 * 5^2
    sage: kronecker_symbol(3,5)
     _2 = -1
    sage: %hist   #Fonctionne depuis le shell mais pas depuis le bloc-note.
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    sage: _oh
     _4 = {1: 2^2 * 5^2, 2: -1}
    sage: _i1
     _5 = 'factor(ZZ(100))\n'
    sage: eval(_i1)
     _6 = 2^2 * 5^2
    sage: %hist
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    4: _oh
    5: _i1
    6: eval(_i1)
    7: %hist

Dans la suite de ce tutoriel et le reste de la documentation de Sage,
nous omettrons la numérotation des sorties.

Il est possible de créer (pour la durée d'une session) une macro qui
rappelle une liste de plusieurs lignes d'entrée.

.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    1: E = EllipticCurve([1,2,3,4,5])
    2: M = ModularSymbols(37)
    3: %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).


.. skip

::

    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E = 5
    sage: M = None
    sage: em
    Executing Macro...
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field

Depuis le shell interactif Sage, il est possible d'exécuter une commande
Unix en la faisant précéder d'un point d'exclamation ``!``. Par exemple,

.. skip

::

    sage: !ls
    auto  example.sage glossary.tex  t  tmp  tut.log  tut.tex

renvoie la liste des fichiers du répertoire courant.

Dans ce contexte, le ``PATH`` commence par le répertoire des binaires
de Sage, de sorte que les commandes ``gp``, ``gap``, ``singular``,
``maxima``, etc. appellent les versions inclues dans Sage.

.. skip

::

    sage: !gp
    Reading GPRC: /etc/gprc ...Done.

                               GP/PARI CALCULATOR Version 2.2.11 (alpha)
                      i686 running linux (ix86/GMP-4.1.4 kernel) 32-bit version
    ...
    sage: !singular
                         SINGULAR                             /  Development
     A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \

Journal des entrées-sorties
===========================

Enregistrer le journal d'une session Sage n'est pas la même chose que
sauvegarder la session (voir :ref:`section-save` pour cette possibilité).
Pour tenir un journal des entrées (et optionnellement des sorties) de
Sage, utilisez la commande ``logstart``. Tapez ``logsatart?`` pour plus
d'informations. Cette commande permet d'enregistrer toutes les entrées que
vous tapez, toutes les sorties, et de rejouer ces entrées dans une
session future (en rechargeant le fichier journal).

.. skip

::

    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: logstart setup
    Activating auto-logging. Current session state plus future input saved.
    Filename       : setup
    Mode           : backup
    Output logging : False
    Timestamping   : False
    State          : active
    sage: E = EllipticCurve([1,2,3,4,5]).minimal_model()
    sage: F = QQ^3
    sage: x,y = QQ['x,y'].gens()
    sage: G = E.gens()
    sage:
    Exiting Sage (CPU time 0m0.61s, Wall time 0m50.39s).
    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: load("setup")
    Loading log file <setup> one line at a time...
    Finished replaying log file <setup>
    sage: E
    Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 + 4*x + 3 over Rational
    Field
    sage: x*y
    x*y
    sage: G
    [(2 : 3 : 1)]

Si vous utilisez le terminal Konsole de KDE, vous pouvez aussi sauver votre
session comme suit : après avoir lancé Sage dans la ``konsole``, ouvrez
le menu « Configuration » et choisissez « Historique... » puis comme
nombre de lignes « Illimité ». Ensuite, lorsque vous souhaitez
enregistrer l'état de votre session, sélectionnez « Enregistrer
l'historique sous... » dans le menu « Édition » et entrez le nom d'un
fichier où enregistrer le texte de votre session. Une fois le fichier
sauvegardé, vous pouvez par exemple l'ouvrir dans un éditeur comme xemacs et
l'imprimer.

Coller du texte ignore les invites
==================================

Imaginons que vous lisiez une session Sage ou Python et que vous vouliez
copier-coller les calculs dans Sage. Le problème est qu'il y a des
invites ``>>>`` ou ``sage:`` en plus des entrées. En fait, vous pouvez
tout à fait copier un exemple complet, invites comprises : par défaut,
l'analyseur syntaxique de Sage supprime les ``>>>`` et ``sage:`` en
début de ligne avant de passer la ligne à Python. Par exemple, les
lignes suivantes sont interprétées correctement :

.. skip

::

    sage: 2^10
    1024
    sage: sage: sage: 2^10
    1024
    sage: >>> 2^10
    1024

Mesure du temps d'exécution d'une commande
==========================================

Si une ligne d'entrée commence par ``%time``, le temps d'exécution de la
commande correspondante est affiché après la sortie. Nous pouvons par
exemple comparer le temps que prend le calcul d'une certaine puissance
entière par diverses méthodes. Les temps de calcul ci-dessous seront
sans doute très différents suivant l'ordinateur, voire la version de
Sage utilisés. Premièrement, en pur Python :

.. skip

::

    sage: %time a = int(1938)^int(99484)
    CPU times: user 0.66 s, sys: 0.00 s, total: 0.66 s
    Wall time: 0.66

Le calcul a pris 0.66 seconde, pendant un intervalle de *wall time* (le
temps de votre montre) lui aussi de 0.66 seconde. Si d'autres programmes
qui s'exécutent en même temps que Sage chargent l'ordinateur avec de gros
calculs, le *wall time* peut être nettement plus important que le temps
processeur.

Chronométrons maintenant le calcul de la même puissance avec le type
``Integer`` de Sage, qui est implémenté (en Cython) en utilisant la
bibliothèque GMP :

.. skip

::

    sage: %time a = 1938^99484
    CPU times: user 0.04 s, sys: 0.00 s, total: 0.04 s
    Wall time: 0.04

Avec l'interface à la bibliothèque C PARI :

.. skip

::

    sage: %time a = pari(1938)^pari(99484)
    CPU times: user 0.05 s, sys: 0.00 s, total: 0.05 s
    Wall time: 0.05

GMP est plus rapide, mais de peu (ce n'est pas une surprise, car la
version de PARI incluse dans Sage utilise GMP pour l'arithmétique
entière).

Il est aussi possible de chronométrer tout un bloc de commandes avec la
commande ``cputime``, comme dans l'exemple suivant :

::

    sage: t = cputime()
    sage: a = int(1938)^int(99484)
    sage: b = 1938^99484
    sage: c = pari(1938)^pari(99484)
    sage: cputime(t)                       #random
    0.64

.. skip

::

    sage: cputime?
    ...
        Return the time in CPU second since Sage started, or with optional
        argument t, return the time since time t.
        INPUT:
            t -- (optional) float, time in CPU seconds
        OUTPUT:
            float -- time in CPU seconds

La commande ``walltime`` fonctionne comme ``cputime``, à ceci près
qu'elle mesure le temps total écoulé « à la montre ».

Nous pouvons aussi faire faire le calcul de puissance ci-dessus à chacun
des systèmes de calcul formel inclus dans Sage. Dans chaque cas, nous
commençons par lancer une commande triviale dans le système en question,
de façon à démarrer son serveur. La mesure la plus pertinente est le
*wall time*. Cependant, si la différence entre celui-ci et le temps
processeur est importante, cela peut indiquer un problème de performance
qui mérite d'être examiné.

.. skip

::

    sage: time 1938^99484;
    CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
    Wall time: 0.01
    sage: gp(0)
    0
    sage: time g = gp('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: maxima(0)
    0
    sage: time g = maxima('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.30
    sage: kash(0)
    0
    sage: time g = kash('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: mathematica(0)
            0
    sage: time g = mathematica('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.03
    sage: maple(0)
    0
    sage: time g = maple('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.11
    sage: libgap(0)
    0
    sage: time g = libgap.eval('1938^99484;')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 1.02

Nous voyons que GAP et Maxima sont les plus lents sur ce test (lancé sur
la machine ``sage.math.washington.edu``). Mais en raison du surcoût de
l'interface pexpect, la comparaison avec Sage, qui est le plus rapide,
n'est pas vraiment équitable.

Trucs et astuces IPython
========================

Comme signalé plus haut, Sage utilise l'interpréteur de commandes IPython, et
met donc à votre disposition toutes les commandes et fonctionnalités de
celui-ci. Vous voudrez peut-être consulter la `documentation complète de IPython
<https://ipython.readthedocs.io/en/stable/>`_. Voici en attendant quelques
astuces utiles -- qui reposent sur ce que IPython appelle des « commandes
magiques » :

- Lorsque l'on souhaite saisir un morceau de code complexe, on peut utiliser
  ``%edit`` (ou ``%ed``, ou ``ed``) pour ouvrir un éditeur de texte.
  Assurez-vous que la variable d'environnement :envvar:`EDITOR` est réglée à
  votre éditeur favori au démarrage de Sage (en plaçant si nécessaire quelque
  chose du genre ``export EDITOR=/usr/bin/emacs`` ou encore  ``export
  EDITOR=/usr/bin/vim`` dans un fichier de configuration convenable, par
  exemple ``.profile``). La commande ``%edit`` à l'invite de Sage ouvrira
  l'éditeur sélectionné. Vous pouvez alors par exemple saisir une définition de
  fonction::

    def some_function(n):
        return n**2 + 3*n + 2

  puis enregistrer le fichier et quitter l'éditeur. La fonction
  ``some_function`` est désormais disponible dans votre session Sage, et vous
  pouvez la modifier en saisissant ``edit some_function`` à l'invite de
  commande.

- Si vous souhaitez reprendre une version modifiée du résultat d'un calcul dans
  une nouvelle commande, tapez ``%rep`` après avoir fait le calcul. Cela
  récupère le texte du résultat et le place sur la ligne de commande, prêt à
  être modifié. ::

    sage: f(x) = cos(x)
    sage: f(x).derivative(x)
    -sin(x)

  Ainsi, après les commandes ci-dessus, la commande ``%rep`` fournit un nouvel
  invite de commande pré-rempli avec le texte ``-sin(x)`` et le curseur en fin
  de ligne.

Pour plus d'information, entrez la commande ``%quickref`` pour un résumé des
possibilités de IPython. Au moment où cette documentation est écrite
(avril 2011), Sage emploie IPython 0.9.1. La `documentation des commandes
magiques
<https://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions>`_
est disponible en ligne, et divers aspects un peu plus avancés de leur
fonctionnement sont décrits `ici  <https://ipython.org/ipython-doc/stable/interactive/reference.html#magic-command-system>`_.

Erreurs et exceptions
=====================

Quand quelque chose ne marche pas, cela se manifeste habituellement par
une « exception » Python. Python essaie de plus de donner une idée de ce
qui a pu déclencher l'exception. Bien souvent, il affiche le nom de
l'exception (par exemple :class:`NameError` ou :class:`ValueError`, voir
le manuel de référence de la bibliothèque de Python [PyLR]_ pour une liste
complète). Par exemple :

.. skip

::

    sage: 3_2
    ------------------------------------------------------------
       File "<console>", line 1
         ZZ(3)_2
               ^
    SyntaxError: invalid ...

    sage: EllipticCurve([0,infinity])
    ------------------------------------------------------------
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce Infinity (<class 'sage...Infinity'>) to Rational

Le débogueur interactif est parfois utile pour comprendre ce qu'il s'est
passé. Il s'active ou se désactive avec ``%pdb`` (et est désactivé par
défaut). L'invite ``ipdb>>`` du débogueur apparaît si une exception a
lieu alors que celui-ci est actif. Le débogueur permet d'afficher l'état
de n'importe quelle variable locale et de monter ou descendre dans la
pile d'exécution. Par exemple :

.. skip

::

    sage: %pdb
    Automatic pdb calling has been turned ON
    sage: EllipticCurve([1,infinity])
    ---------------------------------------------------------------------------
    <class 'exceptions.TypeError'>             Traceback (most recent call last)
    ...

    ipdb>

Pour obtenir une liste des commandes disponibles dans le débogueur,
tapez ``?`` à l'invite ``ipdb>`` :

::

    ipdb> ?

    Documented commands (type help <topic>):
    ========================================
    EOF    break  commands   debug    h       l     pdef   quit    tbreak
    a      bt     condition  disable  help    list  pdoc   r       u
    alias  c      cont       down     ignore  n     pinfo  return  unalias
    args   cl     continue   enable   j       next  pp     s       up
    b      clear  d          exit     jump    p     q      step    w
    whatis where

    Miscellaneous help topics:
    ==========================
    exec  pdb

    Undocumented commands:
    ======================
    retval  rv

Tapez Ctrl-D ou ``quit`` pour revenir à Sage.

.. _section-tabcompletion:

Recherche en arrière et complétion de ligne de commande
=======================================================

Commençons par créer l'espace vectoriel de dimension trois
:math:`V=\QQ^3` comme suit :

::

    sage: V = VectorSpace(QQ,3)
    sage: V
    Vector space of dimension 3 over Rational Field

Nous pouvons aussi utiliser la variante plus concise :

::

    sage: V = QQ^3

Tapez ensuite le début d'une commande, puis ``Ctrl-p`` (ou flèche vers
le haut) pour passer en revue les lignes qui commencent par les mêmes
lettres parmi celles que vous avez entrées jusque-là. Cela fonctionne
même si vous avez quitté et relancé Sage entre-temps. Vous pouvez aussi
rechercher une portion de commande en remontant dans l'historique avec
``Ctrl-r``. Toutes ces fonctionnalités reposent sur la bibliothèque
``readline``, qui existe pour la plupart des variantes de Linux.

La complétion de ligne de commande permet d'obtenir facilement la liste
des fonctions membres de :math:`V` : tapez simplement ``V.`` puis
appuyez sur la touche tabulation.

.. skip

::

    sage: V.[tab key]
    V._VectorSpace_generic__base_field
    ...
    V.ambient_space
    V.base_field
    V.base_ring
    V.basis
    V.coordinates
    ...
    V.zero_vector

Si vous tapez les quelques premières lettres d'un nom de fonction avant
d'appuyer sur ``tab``, vous n'obtiendrez que les fonctions qui
commencent par ces quelques lettres :

.. skip

::

    sage: V.i[tab key]
    V.is_ambient  V.is_dense    V.is_full     V.is_sparse

Si vous cherchez à savoir ce que fait une fonction, par exemple la
fonction coordinates, ``V.coordinates?`` affiche un message d'aide et
``V.coordinates??`` le code source de la fonction, comme expliqué dans
la section suivante.


Aide en ligne
=============

Sage dispose d'un système d'aide intégré. Pour obtenir la documentation
d'une fonction, tapez son nom suivi d'un point d'interrogation.

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <bound method FreeModule_ambient_field.coordinates of Vector
    space of dimension 3 over Rational Field>
    Namespace:      Interactive
    File:           /home/was/s/local/lib/python2.4/site-packages/sage/modules/f
    ree_module.py
    Definition:     V.coordinates(self, v)
    Docstring:
        Write v in terms of the basis for self.

        Returns a list c such that if B is the basis for self, then

                sum c_i B_i = v.

        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: M = FreeModule(IntegerRing(), 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]

Comme nous pouvons le voir ci-dessus, la sortie indique le type de
l'objet, le nom du fichier où il est défini, et donne une description
de l'effet de la fonction, avec des exemples que vous pouvez copier dans
votre session Sage. Pratiquement tous ces exemples sont automatiquement
testés régulièrement pour s'assurer qu'ils se comportent exactement
comme indiqué.

Une autre fonctionnalité, nettement dans l'esprit du caractère ouvert de
Sage, est que lorsque ``f`` est une fonction Python, taper ``f??``
affiche son code source. Par exemple,

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates??
    Type:           instancemethod
    ...
    Source:
    def coordinates(self, v):
            """
            Write $v$ in terms of the basis for self.
            ...
            """
            return self.coordinate_vector(v).list()

Nous voyons que la fonction ``coordinates`` ne fait qu'appeler
``coordinate_vector`` et transformer le résultat en une liste. Mais alors,
que fait la fonction ``coordinate_vector`` ?

.. skip

::

    sage: V = QQ^3
    sage: V.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            ...
            return self.ambient_vector_space()(v)

La fonction ``coordinate_vector`` convertit son entrée en un élément de
l'espace ambiant, ce qui a pour effet de calculer le vecteur des
coefficients de  :math:`v` dans :math:`V`. L'espace :math:`V` est déjà
« l'espace ambiant » puisque c'est simplement :math:`\QQ^3`. Il y
a aussi une fonction ``coordinate_vector`` différente pour les
sous-espaces. Créons un sous-espace et examinons-là :


.. skip

::

    sage: V = QQ^3; W = V.span_of_basis([V.0, V.1])
    sage: W.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            """
             ...
            """
            # First find the coordinates of v wrt echelon basis.
            w = self.echelon_coordinate_vector(v)
            # Next use transformation matrix from echelon basis to
            # user basis.
            T = self.echelon_to_user_matrix()
            return T.linear_combination_of_rows(w)

(Si vous pensez que cette implémentation est inefficace, venez nous
aider à optimiser l'algèbre linéaire !)

Vous pouvez aussi taper ``help(commande)`` ou ``help(classe)`` pour
appeler une sorte de page de manuel relative à une commande ou une
classe.

.. skip

::

    sage: help(VectorSpace)
    Help on function VectorSpace in module sage.modules.free_module:

    VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args)
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x
         over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
    --More--

Pour quitter la page d'aide, appuyez sur ``q``. Votre session revient à
l'écran comme elle était : contrairement à la sortie de ``fonction?``,
celle de ``help`` n'encombre pas votre session. Une possibilité
particulièrement utile est de consulter l'aide d'un module entier avec
``help(nom_du_module``. Par exemple, les espaces vectoriels sont définis
dans  ``sage.modules.free_module``, et on accède à la documentation de
ce module en tapant ``help(sage.modules.free_module)``. Lorsque vous
lisez une page de documentation avec la commande ``help``, vous pouvez
faire des recherches en avant en tapant ``/`` et en arrière en tapant
``?``.


Enregistrer et charger des objets individuellement
==================================================

Imaginons que nous calculions une matrice, ou pire, un espace compliqué
de symboles modulaires, et que nous souhaitions les sauvegarder pour
un usage futur. Les systèmes de calcul formel ont différentes approches
pour permettre cela.


#. **Sauver la partie :** il n'est possible de sauver que la session
   entière (p.ex. GAP, Magma).

#. **Format d'entrée/sortie unifié :** chaque objet est
   affiché sous une forme qui peut être relue (GP/PARI).

#. **Eval :** permettre d'évaluer facilement du code arbitraire dans
   l'interpréteur (p.ex. Singular, PARI).


Utilisant Python, Sage adopte une approche différente, à savoir que tous
les objets peuvent être sérialisés, i.e. transformés en chaînes de
caractères à partir desquelles ils peuvent être reconstruits. C'est une
méthode semblable dans l'esprit à l'unification des entrées et sorties
de PARI, avec l'avantage que l'affichage normal des objets n'a pas
besoin d'être trop compliqué. En outre, cette fonction de sauvegarde et
de relecture des objets ne nécessite (dans la plupart des cas) aucune
programmation supplémentaire : il s'agit simplement une fonctionnalité de
Python fournie par le langage depuis la base.

Quasiment n'importe quel objet Sage ``x`` peut être enregistré sur le
disque, dans un format compressé, avec ``save(x, nom_de_fichier)`` (ou
dans bien des cas ``x.save(nom_de_fichier)``). Pour recharger les
objets, on utilise ``load(nom_de_fichier)``.


.. skip

::

    sage: A = MatrixSpace(QQ,3)(range(9))^2
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]
    sage: save(A, 'A')

Quittez puis redémarrez maintenant Sage. Vous pouvez récupérer ``A`` :

.. skip

::

    sage: A = load('A')
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]

Vous pouvez faire de même avec des objets plus compliqués, par exemple
des courbes elliptiques. Toute l'information en cache sur l'objet est
stockée avec celui-ci :

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: v = E.anlist(100000)              # prend un moment
    sage: save(E, 'E')
    sage: quit

Ainsi, la version sauvegardée de ``E`` prend 153 kilo-octets car elle
contient les 100000 premiers :math:`a_n`.

.. skip

::

    ~/tmp$ ls -l E.sobj
    -rw-r--r--  1 was was 153500 2006-01-28 19:23 E.sobj
    ~/tmp$ sage [...]
    sage: E = load('E')
    sage: v = E.anlist(100000)              # instantané !

(En Python, les sauvegardes et rechargements s'effectuent à l'aide du
module ``pickle``. En particulier, on peut sauver un objet Sage ``x``
par la commande ``pickle.dumps(x, 2)``.  Attention au ``2`` !)

Sage n'est pas capable de sauvegarder les objets créés dans d'autres systèmes
de calcul formel comme GAP, Singular, Maxima etc. : au rechargement, ils
sont dans un état marqué « invalide ». Concernant GAP, un certain
nombre d'objets sont affichés sous une forme qui permet de les
reconstruire, mais d'autres non, aussi la reconstruction d'objets GAP
à partir de leur affichage est intentionnellement interdite.

.. skip

::

    sage: a = libgap(2)
    sage: a.save('a')
    sage: load('a')
    Traceback (most recent call last):
    ...
    ValueError: The session in which this object was defined is no longer
    running.

Les objets GP/PARI, en revanche, peuvent être sauvegardés et rechargés,
puisque la forme imprimée d'un objet suffit à reconstruire celui-ci.

.. skip

::

    sage: a = gp(2)
    sage: a.save('a')
    sage: load('a')
    2

Un objet sauvegardé peut être rechargé y compris sur un ordinateur doté
d'une architecture ou d'un système d'exploitation différent. Ainsi, il
est possible de sauvegarder une immense matrice sur un OS-X 32 bits, la
recharger sur un Linux 64 bits, l'y mettre en forme échelon et rapatrier
le résultat. Bien souvent, un objet peut même être rechargé avec une
version de Sage différente de celle utilisée pour le sauver, pourvu que
le code qui gère cet objet n'ait pas trop changé d'une version sur
l'autre. Sauver un objet enregistre tous ses attributs ainsi que la
classe à laquelle il appartient (mais pas son code source). Si cette
classe n'existe plus dans une version ultérieure de Sage, l'objet ne
peut pas y être rechargé. Mais il demeure possible de le charger dans
l'ancienne version pour récupérer son dictionnaire (avec
``x.__dict__``), sauver celui-ci, et le recharger dans la nouvelle
version.

Enregistrer un objet comme texte
--------------------------------

Une autre possibilité consiste à sauvegarder la représentation texte ASCII
dans un fichier texte brut, ce qui se fait simplement en ouvrant le
fichier en écriture et en y écrivant la représentation de l'objet (il
est tout à fait possible d'écrire plusieurs objets). Une fois l'écriture
terminée, nous refermons le fichier.

.. skip

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = (x+y)^7
    sage: o = open('file.txt','w')
    sage: o.write(str(f))
    sage: o.close()

.. _section-save:

Enregister et recharger des sessions entières
=============================================

Sage dispose de fonctions très souples de sauvegarde et relecture de
sessions entières.

La commande ``save_session(nom_de_session)`` enregistre toutes les
variables définies dans la session courante sous forme de dictionnaire
dans le fichier ``nom_de_session.sobj``. (Les éventuelles variables qui
ne supportent pas la sauvegarde sont ignorées.) Le fichier ``.sobj`` obtenu
peut être rechargé comme n'importe quel objet sauvegardé ; on obtient en
le rechargeant un dictionnaire dont les clés sont les noms de variables
et les valeurs les objets correspondants.

La commande ``reload_session(nom_de_session)`` charge toutes les
variables sauvées dans ``nom_de_session``. Cela n'efface pas les
variables déjà définies dans la session courante : les deux sessions
sont fusionnées.

Commençons par démarrer Sage et par définir quelques variables.

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: M = ModularSymbols(37)
    sage: a = 389
    sage: t = M.T(2003).matrix(); t.charpoly().factor()
     _4 = (x - 2004) * (x - 12)^2 * (x + 54)^2

Nous sauvons maintenant notre session, ce qui a pour effet d'enregistrer
dans un même fichier toutes les variables ci-dessus. Nous pouvons
constater que le fichier fait environ 3 ko.

.. skip

::

    sage: save_session('misc')
    Saving a
    Saving M
    Saving t
    Saving E
    sage: quit
    was@form:~/tmp$ ls -l misc.sobj
    -rw-r--r--  1 was was 2979 2006-01-28 19:47 misc.sobj

Enfin, nous redémarrons Sage, nous définissons une nouvelle variable, et
nous rechargeons la session précédente.

.. skip

::

    sage: b = 19
    sage: load_session('misc')
    Loading a
    Loading M
    Loading E
    Loading t

Toutes les variables sauvegardées sont à nouveau disponibles. En outre,
la variable ``b`` n'a pas été écrasée.

.. skip

::

    sage: M
    Full Modular Symbols space for Gamma_0(37) of weight 2 with sign 0
    and dimension 5 over Rational Field
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational
    Field
    sage: b
    19
    sage: a
    389

