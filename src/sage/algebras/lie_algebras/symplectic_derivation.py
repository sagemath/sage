# -*- coding: utf-8 -*-
r"""
Symplectic Derivation Lie Algebras

AUTHORS:

- Travis Scrimshaw (2020-10): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2020 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.lie_algebras import LieAlgebras
from sage.sets.family import Family
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.indexed_generators import IndexedGenerators
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import InfinitelyGeneratedLieAlgebra
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.combinat.partition import _Partitions, Partitions


class SymplecticDerivationLieAlgebra(InfinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    The symplectic derivation Lie algebra.

    Fix a `g \geq 4` and let `R` be a commutative ring. Let `H = R^{2g}`
    be equipped with a symplectic form `\mu` with the basis
    `a_1, \ldots, a_g, b_1, \ldots, b_g` such that

    .. MATH::

        \mu(a_i, a_j) = \mu(b_i, b_j) = 0,
        \qquad\qquad
        \mu(a_i, b_j) = -\mu(b_j, a_i) = \delta_{ij},

    for all `i, j`. The *symplectic derivation Lie algebra* is the Lie
    algebra

    .. MATH::

        \mathfrak{c}_g := \bigoplus_{w \geq 0} S^{w+2} H

    with the Lie bracket on basis elements

    .. MATH::

        [x_1 \cdots x_{m+2}, y_1 \cdots y_{n+2}] =
        \sum_{i,j} \mu(x_i, y_j) x_1 \cdots \widehat{x}_i \cdots x_{m+2}
        \cdot y_1 \cdots \widehat{y}_j \cdots y_{n+2},

    where `\widehat{z}` denotes that factor is missing. When `R = \QQ`, this
    corresponds to the classical Poisson bracket on `C^{\infty}(\RR^{2g})`
    restricted to polynomials with coefficients in `\QQ`.

    EXAMPLES::

        sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
        sage: elts = L.some_elements()
        sage: list(elts)
        [a1*a2, b1*b3, a1*a1*a2, b3*b4,
         a1*a4*b3, a1*a2 - 1/2*a1*a2*a2*a5 + a1*a1*a2*b1*b4]
        sage: [[elts[i].bracket(elts[j]) for i in range(len(elts))]
        ....:  for j in range(len(elts))]
        [[0, -a2*b3, 0, 0, 0, -a1*a1*a2*a2*b4],
         [a2*b3, 0, 2*a1*a2*b3, 0, a4*b3*b3, a2*b3 - 1/2*a2*a2*a5*b3 + 2*a1*a2*b1*b3*b4],
         [0, -2*a1*a2*b3, 0, 0, 0, -2*a1*a1*a1*a2*a2*b4],
         [0, 0, 0, 0, a1*b3*b3, 0],
         [0, -a4*b3*b3, 0, -a1*b3*b3, 0, -a1*a1*a1*a2*b1*b3 - a1*a1*a2*a4*b3*b4],
         [a1*a1*a2*a2*b4, -a2*b3 + 1/2*a2*a2*a5*b3 - 2*a1*a2*b1*b3*b4, 2*a1*a1*a1*a2*a2*b4,
          0, a1*a1*a1*a2*b1*b3 + a1*a1*a2*a4*b3*b4, 0]]
        sage: x = L.monomial(Partition([8,8,6,6,4,2,2,1,1,1])); x
        a1*a1*a1*a2*a2*a4*b1*b1*b3*b3
        sage: [L[x, elt] for elt in elts]
        [-2*a1*a1*a1*a2*a2*a2*a4*b1*b3*b3,
         3*a1*a1*a2*a2*a4*b1*b1*b3*b3*b3,
         -4*a1*a1*a1*a1*a2*a2*a2*a4*b1*b3*b3,
         a1*a1*a1*a2*a2*b1*b1*b3*b3*b3,
         -2*a1*a1*a1*a2*a2*a4*a4*b1*b3*b3*b3,
         -2*a1*a1*a1*a2*a2*a2*a4*b1*b3*b3 + a1*a1*a1*a2*a2*a2*a2*a4*a5*b1*b3*b3
          + a1*a1*a1*a1*a1*a2*a2*a2*b1*b1*b1*b3*b3 - a1*a1*a1*a1*a2*a2*a2*a4*b1*b1*b3*b3*b4]

    REFERENCES:

    - [Harako2020]_
    """
    def __init__(self, R, g):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: TestSuite(L).run()
        """
        if g < 4:
            raise ValueError("g must be at least 4")
        cat = LieAlgebras(R).WithBasis().Graded()
        self._g = g
        d = Family(NonNegativeIntegers(), lambda n: Partitions(n, min_length=2, max_part=2*g))
        indices = DisjointUnionEnumeratedSets(d)
        InfinitelyGeneratedLieAlgebra.__init__(self, R, index_set=indices, category=cat)
        IndexedGenerators.__init__(self, indices, sorting_key=self._basis_key)

    def _basis_key(self, x):
        r"""
        Return the key used to compare two basis element indices.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L._basis_key( [7, 5, 2, 1] )
            (4, [7, 5, 2, 1])

            sage: x = L.an_element(); x
            a1*a2 - 1/2*a1*a2*a2*a5 + a1*a1*a2*b1*b4
            sage: sorted(map(L._basis_key, x.support()))
            [(2, [2, 1]), (4, [5, 2, 2, 1]), (5, [9, 6, 2, 1, 1])]
        """
        return (len(x), x)

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L._repr_term([7, 5, 2, 1])
            'a1*a2*a5*b2'
        """
        g = self._g
        def label(i):
            return "a{}".format(i) if i <= g else "b{}".format(i-g)
        return "*".join(label(i) for i in reversed(m))

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the term indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L._latex_term([7, 5, 2, 1])
            'a_{1} a_{2} a_{5} b_{2}'
        """
        g = self._g
        def label(i):
            return "a_{{{}}}".format(i) if i <= g else "b_{{{}}}".format(i-g)
        return " ".join(label(i) for i in reversed(m))

    def _unicode_art_term(self, m):
        r"""
        Return a unicode art representation of the term indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L._unicode_art_term([7, 5, 2, 1])
            a₁·a₂·a₅·b₂
        """
        from sage.typeset.unicode_art import unicode_art, unicode_subscript
        g = self._g
        def label(i):
            return "a{}".format(unicode_subscript(i)) if i <= g else "b{}".format(unicode_subscript(i-g))
        return unicode_art("·".join(label(i) for i in reversed(m)))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.SymplecticDerivation(QQ, 5)
            Symplectic derivation Lie algebra of rank 5 over Rational Field
        """
        return "Symplectic derivation Lie algebra of rank {} over {}".format(self._g, self.base_ring())

    def degree_on_basis(self, x):
        r"""
        Return the degree of the basis element indexed by ``x``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L.degree_on_basis([5,2,1])
            1
            sage: L.degree_on_basis([1,1])
            0
            sage: elt = L.monomial(Partition([5,5,2,1])) + 3*L.monomial(Partition([3,3,2,1]))
            sage: elt.degree()
            2
        """
        return len(x) - 2

    def bracket_on_basis(self, x, y):
        r"""
        Return the bracket of basis elements indexed by ``x`` and ``y``,
        where ``i < j``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L.bracket_on_basis([5,2,1], [5,1,1])
            0
            sage: L.bracket_on_basis([6,1], [3,1,1])
            -2*a1*a1*a3
            sage: L.bracket_on_basis([9,2,1], [4,1,1])
            -a1*a1*a1*a2
            sage: L.bracket_on_basis([5,5,2], [6,1,1])
            0
            sage: L.bracket_on_basis([5,5,5], [10,3])
            3*a3*a5*a5
            sage: L.bracket_on_basis([10,10,10], [5,3])
            -3*a3*b5*b5
        """
        g = self._g
        ret = {}
        one = self.base_ring().one()
        for i,xi in enumerate(x):
            for j,yj in enumerate(y):
                # The symplectic form will be 0
                if (xi <= g and yj <= g) or (xi > g and yj > g):
                    continue
                if xi <= g and yj > g:
                    if xi != yj - g:
                        continue
                    m = _Partitions(sorted(x[:i] + x[i+1:] + y[:j] + y[j+1:], reverse=True))
                    if m in ret:
                        ret[m] += one
                    else:
                        ret[m] = one
                else: # if ci > g and yj <= g:
                    if xi - g != yj:
                        continue
                    m = _Partitions(sorted(x[:i] + x[i+1:] + y[:j] + y[j+1:], reverse=True))
                    if m in ret:
                        ret[m] -= one
                    else:
                        ret[m] = -one
        return self._from_dict(ret, remove_zeros=True)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L.an_element()
            a1*a2 - 1/2*a1*a2*a2*a5 + a1*a1*a2*b1*b4
        """
        d = self.monomial
        return (
                 d( _Partitions([2,1]) )
                 - self.base_ring().an_element() * d( _Partitions([5,2,2,1]) )
                 + d( _Partitions([2*self._g-1, self._g+1, 2, 1, 1]) )
                )

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.SymplecticDerivation(QQ, 5)
            sage: L.some_elements()
            [a1*a2, b1*b3, a1*a1*a2, b3*b4, a1*a4*b3,
             a1*a2 - 1/2*a1*a2*a2*a5 + a1*a1*a2*b1*b4]
        """
        d = self.monomial
        g = self._g
        return [d( _Partitions([2,1]) ), d( _Partitions([g+3,g+1]) ), d( _Partitions([2,1,1])),
                d( _Partitions([2*g-1,2*g-2]) ), d( _Partitions([2*g-2,g-1,1]) ),
                self.an_element()]

    class Element(LieAlgebraElement):
        pass

