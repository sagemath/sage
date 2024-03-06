.. _sec-sagetex:

****************
Utiliser SageTeX
****************

Le paquet SageTeX permet d'inclure dans un document LaTeX les résultats
de calculs effectués avec Sage. Il est fourni avec Sage, mais pour
l'utiliser, vous aurez besoin de l'ajouter à votre installation TeX.
Cette opération se résume à copier un fichier ; voyez la section
:ref:`installation` du présent tutoriel ainsi que "Make SageTeX known to
TeX" dans le guide d'installation de Sage (`Sage installation guide
<http://doc.sagemath.org/html/en/installation/>`_, `ce lien
<../../en/installation/index.html>`_ devrait conduire à une copie
locale) pour plus de détails.

Voici un bref exemple d'utilisation de SageTeX. La documentation
complète se trouve dans
:file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`, où ``SAGE_ROOT``
désigne le répertoire racine de votre installation Sage. Elle est
accompagnée d'un fichier exemple et de scripts Python potentiellement
utiles.

Pour essayer SageTeX, suivez les instructions d'installation puis copiez
le texte suivant dans un fichier ``st_example.tex``:

.. warning::

    Attention, tel qu'il est affiché dans la version interactive de
    l'aide de Sage, l'exemple suivant provoque des erreurs lors de la
    compilation par LaTeX. Consultez la version statique de la
    documentation pour voir le texte correct.


.. code-block:: latex

    \documentclass{article}
    \usepackage[T1]{fontenc}
    \usepackage[frenchb]{babel}
    \usepackage{sagetex}

    \begin{document}

    Sage\TeX{} permet de faire des calculs avec Sage en pla\,cant
    automatiquement les r\'esultats dans un document \LaTeX.
    Par exemple, l'entier $1269$ admet
    $\sage{number_of_partitions(1269)}$ partitions. Le r\'esultat de ce
    calcul est affich\'e sans que vous ayez \`a le calculer
    vous-m\^eme ou m\^eme \`a le copier-coller dans votre document.

    Un peu de code Sage :

    \begin{sageblock}
        f(x) = exp(x) * sin(2*x)
    \end{sageblock}

    La d\'eriv\'ee seconde de $f$ est

    \[
      \frac{\mathrm{d}^{2}}{\mathrm{d}x^{2}} \sage{f(x)} =
      \sage{diff(f, x, 2)(x)}.
    \]

    Voici enfin le graphe de $f$ sur $[-1,1]$:

    \sageplot{plot(f, -1, 1)}

    \end{document}

Lancez LaTeX pour compiler ``st_example.tex`` comme à l'accoutumée.
LaTeX va afficher un certain nombre de messages d'avertissement, dont::

    Package sagetex Warning: Graphics file
    sage-plots-for-st_example.tex/plot-0.pdf on page 1 does not exist.
    Plot command is on input line 30.

    Package sagetex Warning: There were undefined Sage formulas and/or
    plots. Run Sage on st_example.sage, and then run LaTeX on
    st_example.tex again.

En plus des fichiers habituellement produits par LaTeX, le répertoire où
vous travaillez contient maintenant un fichier ``st_example.sage``. Il
s'agit d'un script Sage produit par la compilation du fichier LaTeX.
L'avertissement précédent vous demande de lancer Sage sur
``st_example.sage``, faites-le. Un nouveau message vous demande de
relancer LaTeX sur ``st_example.tex``, mais avant de le faire, observez
qu'un nouveau fichier ``st_example.sout`` a été créé. C'est ce fichier
qui contient le résultat des calculs effectués par Sage, dans un format
que LaTeX est capable de lire pour insérer les résultats dans votre
document. Un nouveau répertoire contenant un fichier EPS avec votre
graphique est également apparu. Après une nouvelle exécution de LaTeX,
les résultats des commandes Sage sont présents dans votre document
compilé.

Les macros LaTeX utilisées dans l'exemple ci-dessus ne sont guère
compliquées à comprendre. Un environnement ``sageblock`` compose le code
qu'il contient verbatim, et le fait exécuter par Sage. Avec
``\sage{toto}``, le résultat placé dans votre document est celui de la
commande Sage ``latex(toto)``. Les commandes de tracé de graphiques sont
un peu plus compliquées. L'exemple ci-dessus utilise la forme la plus
simple, ``\sageplot{toto}``, qui insère dans le document l'image obtenue
sous Sage par ``toto.save('fichier.eps')``.

Pour utiliser SageTeX, la procédure est donc la suivante :

    - lancez LaTeX sur votre fichier .tex ;
    - lancez Sage sur le fichier .sage produit pas LaTeX ;
    - lancez LaTeX une seconde fois.

(Il n'est pas nécessaire de lancer Sage si vous recompilez un document
LaTeX sans avoir modifié les commandes Sage qu'il contient depuis la
compilation précédente.)

SageTeX offre bien d'autres possibilités. Puisque Sage
comme LaTeX sont des outils complexes et puissants, le mieux est sans
doute de consulter la documentation complète de SageTeX, qui se trouve
dans :file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`.
