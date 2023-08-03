********
Postface
********

Pourquoi Python ?
=================

Les avantages de Python
-----------------------

Le langage d'implémentation de la base de Sage est le langage Python (voir
[Py]_), même si le code qui doit s'exécuter rapidement est écrit
dans un langage compilé. Python présente plusieurs avantages :

-  L'**enregistrement d'objets** est très facile en Python. Il existe en
   Python un vaste support pour enregistrer (presque) n'importe quel
   objet dans des fichiers sur le disque ou dans une base de données.

-  Python fournit d'excellents outils pour la  **documentation** des
   fonctions et des packages du code source, ce qui comprend
   l'extraction automatique de documentation et le test
   automatique de tous les exemples. Ces tests automatiques sont
   exécutés régulièrement de façon automatique, ce qui garantit que
   les exemples donnés dans la documentation fonctionnent comme
   indiqué.

-  **Gestion de la mémoire :** Python possède désormais un gestionnaire
   de mémoire bien pensé et robuste ainsi qu'un ramasse-miettes (*garbage
   collector*) qui traite correctement les références circulaires et
   tient compte des variables locales dans les fichiers.

-  **Énormément de packages** d'ores et déjà disponibles pour Python
   pourraient se révéler d'un grand intérêt pour les utilisateurs de
   Sage : analyse numérique et algèbre linéaire, visualisation 2D et 3D,
   réseau (pour le calcul distribué, la mise en place de serveurs - par
   exemple twisted), bases de données, etc.

-  **Portabilité :** Python se compile sans difficulté et en quelques
   minutes sur la plupart des plates-formes.

-  **Gestion des exception :** Python possède un système sophistiqué et
   bien pensé de gestion des exceptions, grâce auquel les programmes
   peuvent rétablir leur fonctionnement normal même si des erreurs
   surviennent dans le code qu'ils appellent.

-  **Débogueur :** Python comprend un débogueur. Ainsi, quand un
   programme échoue pour une raison quelconque, l'utilisateur peut
   consulter la trace complète de la pile d'exécution, inspecter l'état de
   toutes les variables pertinentes et se déplacer dans la pile.

-  **Profileur :** Il existe un profileur Python, qui exécute le code et
   renvoie un rapport qui détaille combien de fois et pendant combien de
   temps chaque fonction a été appelée.

-  **Un langage :** Au lieu d'écrire un **nouveau langage** pour les
   mathématiques comme cela a été fait pour Magma, Maple, Mathematica,
   Matlab, GP/PARI, GAP, Macaulay 2, Simath, etc., nous utilisons le
   langage Python, qui est un langage de programmation répandu,
   activement développé et optimisé par des centaines de développeurs
   qualifiés. Python, avec son processus de développement
   éprouvé, fait partie des success stories majeures de
   l'open source (see [PyDev]_).


.. _section-mathannoy:

Le préprocesseur Sage et les différences entre Sage et Python
-------------------------------------------------------------

Certains aspects mathématiques de Python peuvent induire des confusions.
Aussi, Sage se comporte différemment de Python à plusieurs égards.


-  **Notation de l'exponentiation :** ``**`` au lieu de ``^``. En Python,
   ``^`` désigne le "xor" (ou exclusif bit à bit) et non
   l'exponentiation. Ainsi en Python, on a

   ::

       >>> 2^8
       10
       >>> 3^2
       1
       >>> 3**2
       9

   Cette utilisation de ``^`` peut paraître étrange et surtout
   inefficace pour une utilisation purement mathématique puisque le ou
   exclusif n'est que rarement utilisé. Par commodité, Sage prétraite
   chaque ligne de commande avant de la transmettre
   à Python, en remplaçant par exemple les apparitions de ``^``
   (en dehors des chaînes de caractères) par des ``**``:

   ::

       sage: 2^8
       256
       sage: 3^2
       9
       sage: "3^2"
       '3^2'

    Le ou exclusif bit à bit est quant à lui noté ``^^``, et l'opération en
    place ``^^=`` fonctionne comme on s'y attend :

    ::

        sage: 3^^2
        1
        sage: a = 2
        sage: a ^^= 8
        sage: a
        10

- **Division entière :** L'expression Python ``2/3`` ne se comporte
   pas de la manière à laquelle s'attendraient des mathématiciens. En
   Python 3, si ``m`` et ``n`` sont de type ``int``, alors ``m/n`` est
   de type ``float``, c'est le quotient réel de ``m`` par ``n``. Par
   exemple, ``2/3`` renvoie ``0.6666...``. Pour obtenir le quotient
   entier, il faut utiliser ``2//3`` qui renvoie ``0``.

   Dans l'interpréteur Sage, nous réglons cela en encapsulant
   automatiquement les entiers litéraux par ``Integer( )`` et en faisant
   de la division un constructeur pour les nombres rationnels. Par
   exemple :

   ::

       sage: 2/3
       2/3
       sage: (2/3).parent()
       Rational Field
       sage: 2//3
       0

-  **Entiers longs :** Python possède nativement un support pour les entiers de
   précision arbitraire, en plus des int du langage C. Les entiers longs
   Python sont significativement plus lents que ceux que GMP fournit et
   sont marqués à l'affichage par un ``L`` qui les distingue des int (il
   est pas prévu de changer cela à court terme). Sage implémente les
   entiers en précision arbitraire en utilisant la bibliothèque C GMP.
   Les entiers longs GMP utilisés par Sage s'affichent sans le ``L``.

Plutôt que de modifier l'interpréteur Python (comme l'ont fait certaines
personnes pour leurs projets internes), nous utilisons le langage Python
exactement comme il est et rajoutons un pré-parseur pour IPython de sorte
que la ligne de commande de IPython se comporte comme l'attend un
mathématicien. Ceci signifie que tout code Python existant peut être
utilisé sous Sage. Toutefois, il faut toujours respecter les règles
standards de Python lorsque l'on écrit des packages à importer dans
Sage.

(Pour installer une bibliothèque Python, trouvée sur Internet par
exemple, suivez les instructions mais exécutez  ``sage -python`` au lieu
de ``python``.  La plupart du temps, ceci signifie concrètement qu'il
faut taper ``sage -python setup.py install``.)

Comment puis-je contribuer ?
============================

Si vous souhaitez contribuer au développement de Sage, votre aide sera grandement
appréciée ! Cela peut aller de contributions substantielles en code au
signalement de bogues en passant par l'enrichissement de la documentation.

Parcourez la page web de Sage pour y trouver les informations pour les
développeurs. Entre autres choses, vous trouverez une longue liste de
projets en lien avec Sage rangés par priorité et catégorie. Le Guide du
développeur Sage (`Sage Developer's Guide
<http://doc.sagemath.org/html/en/developer/>`_) contient également des
informations utiles. Vous pouvez aussi faire un tour sur le groupe
Google ``sage-devel``.

Comment citer Sage ?
====================

Si vous écrivez un article qui utilise Sage, merci d'y préciser les
calculs faits avec Sage en citant

::

    [Sage] SageMath, the Sage Mathematics Software System (Version 8.7),
           The Sage Developers, 2019, https://www.sagemath.org.

dans votre bibliographie (en remplaçant 8.7 par la version de Sage que
vous avez utilisée). De plus, pensez à rechercher les composants de Sage
que vous avez utilisés pour vos calculs, par exemple PARI, GAP, Singular,
Maxima et citez également ces systèmes. Si vous vous demandez quel
logiciel votre calcul utilise, n'hésitez pas à poser la question sur le
groupe Google ``sage-devel``. Voir :ref:`section-univariate` pour une
discussion plus approfondie de ce point.

------------

Si vous venez de lire d'une traite ce tutoriel et que vous avez une idée du temps qu'il vous a fallu pour le parcourir, merci de nous le faire savoir sur le groupe Google ``sage-devel``.

Amusez-vous bien avec Sage !
