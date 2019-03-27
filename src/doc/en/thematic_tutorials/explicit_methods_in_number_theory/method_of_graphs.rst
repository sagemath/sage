Method of Graphs
================

The Mestre Method of Graphs is an intriguing
algorithm for computing the action of Hecke operators on yet
another module :math:`X` that is isomorphic to
:math:`M_2(\Gamma_0(N))`. The implementation in Sage
unfortunately only works when :math:`N` is prime; in contrast, my
implementation in Magma works when :math:`N=pM` and
:math:`S_2(\Gamma_0(M))=0`.

The matrices of Hecke operators on :math:`X` are vastly sparser
than on any basis of :math:`M_2(\Gamma_0(N))` that you are
likely to use.

::

    sage: X = SupersingularModule(389); X
    Module of supersingular points on X_0(1)/F_389 over Integer Ring
    sage: t2_mestre = X.T(2).matrix()
    sage: p_mestre = t2_mestre.charpoly()
    sage: factor(p_mestre)
    (x - 3) * (x + 2) * (x^2 - 2) * (x^3 - 4*x - 2) * ...

    sage: t2 = ModularSymbols(389,sign=1).hecke_matrix(2)
    sage: p = t2.charpoly()
    sage: factor(p)
    (x - 3) * (x + 2) * (x^2 - 2) * (x^3 - 4*x - 2) * ...

    sage: p_mestre == p
    True

The method of graphs is also used in computer science to construct
expander graphs with good properties. And it is important in my
algorithm for computing Tamagawa numbers of purely toric modular
abelian varieties. This algorithm is not implemented in Sage yet,
since it is only interesting in the case of non-prime level, as it
turns out.

REFERENCE:

- Jean-François Mestre, *La méthode des graphes. Exemples et
  applications*, Proceedings of the international conference on class
  numbers and fundamental units of algebraic number fields (Katata,
  1986), 217-242
