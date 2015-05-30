******
Annexe
******

.. _section-precedence:

Priorité des opérateurs arithmétiques binaires
==============================================

Combien font ``3^2*4 + 2%5`` ?  Le résultat (38) est déterminé par le
« tableau de priorité des opérateurs » suivant. Il est dérivé de celui
donné § 5.14 du manuel de référence de Python (*Python Language
Reference Manual*, de G. Rossum et F. Drake.) Les opérations sont
données par priorités croissantes.


==========================  =================
Opérateur                   Description
==========================  =================
or                          ou booléen
and                         et booléen
not                         négation booléenne
in, not in                  appartenance
is, is not                  test d'identité
>, <=, >, >=, ==, !=        comparaisons
+, -                        addition, soustraction
\*, /, %                    multiplication, division, reste
\*\*, ^                     exponentiation
==========================  =================

Ainsi, pour calculer ``3^2*4 + 2%5``, Sage « met les parenthèses » comme
suit : ``((3^2)*4) + (2%5)``. Il calcule donc d'abord ``3^2``, ce qui
fait ``9``, puis ``(3^2)*4`` et ``2%5``, et enfin ajoute les valeurs de
ces deux dernières expressions.
