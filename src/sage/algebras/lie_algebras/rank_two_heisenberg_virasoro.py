r"""
Rank Two Heisenberg-Virasoro Algebras

AUTHORS:

- Travis Scrimshaw (2018-08): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.lie_algebras import LieAlgebras
from sage.rings.integer_ring import ZZ
from sage.sets.family import Family
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.indexed_generators import IndexedGenerators
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import (InfinitelyGeneratedLieAlgebra)


class RankTwoHeisenbergVirasoro(InfinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    The rank 2 Heisenberg-Virasoro algebra.

    The *rank 2 Heisenberg-Virasoro* (Lie) algebra is the Lie algebra
    `L` spaned by the elements

    .. MATH::

        \{t^{\alpha}, E(\alpha) \mid \alpha \in \ZZ^2 \setminus \{(0,0)\} \}
        \cup \{K_1, K_2, K_3, K_4\},

    which satisfy the relations

    .. MATH::

        \begin{aligned}
        \mbox{ } [t^{\alpha}, t^{\beta}] & = [K_i, L] = 0,
        \\
        [t^{\alpha}, E(\beta)] & =
          \det\begin{pmatrix} \beta \\ \alpha \end{pmatrix} t^{\alpha+\beta}
          + \delta_{\alpha,-\beta} (\alpha_1 K_1 + \alpha_2 K_2),
        \\
        [E(\alpha), E(\beta)] & =
          \det\begin{pmatrix} \beta \\ \alpha \end{pmatrix} E(\alpha+\beta)
          + \delta_{\alpha,-\beta} (\alpha_1 K_3 + \alpha_2 K_4),
        \end{aligned}

    where `\alpha = (\alpha_1, \alpha_2)` and `\delta_{xy}` is the
    Kronecker delta.

    EXAMPLES::

        sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
        sage: K1,K2,K3,K4 = L.K()
        sage: E2m1 = L.E(2,-1)
        sage: Em21 = L.E(-2,1)
        sage: t2m1 = L.t(2,-1)
        sage: t53 = L.t(5,3)

        sage: Em21.bracket(t2m1)
        -2*K1 + K2
        sage: t53.bracket(E2m1)
        11*t(7, 2)
        sage: E2m1.bracket(Em21)
        2*K3 - K4
        sage: E2m1.bracket(t2m1)
        0

        sage: all(x.bracket(y) == 0 for x in [K1,K2,K3,K4] for y in [E2m1, Em21, t2m1])
        True

    REFERENCES:

    - [LT2018]_
    """
    def __init__(self, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: TestSuite(L).run()
        """
        cat = LieAlgebras(R).WithBasis()
        self._KI = FiniteEnumeratedSet([1,2,3,4])
        self._V = ZZ**2
        d = {'K': self._KI, 'E': self._V, 't': self._V}
        indices = DisjointUnionEnumeratedSets(d, keepkey=True, facade=True)
        InfinitelyGeneratedLieAlgebra.__init__(self, R, index_set=indices, category=cat)
        IndexedGenerators.__init__(self, indices, sorting_key=self._basis_key)

    def _basis_key(self, x):
        r"""
        Return the key used to compare two basis element indices.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L._basis_key( ('K',2) )
            (0, 2)
            sage: L._basis_key( ('t',(-1,2)) )
            (1, (-1, 2))
            sage: L._basis_key( ('E',(-1,2)) )
            (2, (-1, 2))

            sage: x = L.an_element(); x
            K3 - 1/2*t(-1, 3) + E(1, -3) + E(2, 2)
            sage: sorted(map(L._basis_key, x.support()))
            [(0, 3), (1, (-1, 3)), (2, (1, -3)), (2, (2, 2))]
        """
        if x[0] == 'K':
            return (0, x[1])
        if x[0] == 't':
            return (1, tuple(x[1]))
        return (2, tuple(x[1])) # x[0] == 'E'

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L._repr_term(('K', 2))
            'K2'
            sage: L._repr_term(('t', (2,-4)))
            't(2, -4)'
            sage: L._repr_term(('E', (2,-4)))
            'E(2, -4)'
        """
        return "{}{}".format(*m)

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the term indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L._latex_term(('K', 2))
            'K_2'
            sage: L._latex_term(('t', (2,-4)))
            't^{(2,-4)}'
            sage: L._latex_term(('E', (2,-4)))
            'E(2,-4)'
        """
        if m[0] == 'K':
            return 'K_{}'.format(m[1])
        if m[0] == 't':
            return 't^{{({},{})}}'.format(m[1][0], m[1][1])
        return 'E({},{})'.format(m[1][0], m[1][1])

    def _unicode_art_term(self, m):
        r"""
        Return a unicode art representation of the term indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L._unicode_art_term(('K', 2))
            K₂
            sage: L._unicode_art_term(('t', (2,-4)))
            t⁽²˴⁻⁴⁾
            sage: L._unicode_art_term(('E', (2,-4)))
            E(2,-4)
        """
        from sage.typeset.unicode_art import unicode_art, unicode_subscript, unicode_superscript
        if m[0] == 'K':
            return unicode_art('K' + unicode_subscript(m[1]))
        if m[0] == 't':
            return unicode_art('t⁽{}˴{}⁾'.format(unicode_superscript(m[1][0]),
                                                unicode_superscript(m[1][1])))
        return unicode_art('E({},{})'.format(m[1][0], m[1][1]))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            The Rank 2 Heisenberg-Virasoro algebra over Rational Field
        """
        return "The Rank 2 Heisenberg-Virasoro algebra over {}".format(self.base_ring())

    @cached_method
    def K(self, i=None):
        r"""
        Return the basis element `K_i` of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L.K(1)
            K1
            sage: list(L.K())
            [K1, K2, K3, K4]
        """
        if i is None:
            return Family(self._KI, self.K)
        return self.monomial( ('K', i) )

    def t(self, a, b):
        r"""
        Return the basis element `t^{(a,b)}` of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L.t(1,-2)
            t(1, -2)
        """
        if a == b == 0:
            raise ValueError("no t(0, 0) element")
        return self.monomial( ('t', self._v(a,b)) )

    def E(self, a, b):
        r"""
        Return the basis element `E(a,b)` of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L.E(1,-2)
            E(1, -2)
        """
        if a == b == 0:
            raise ValueError("no E(0, 0) element")
        return self.monomial( ('E', self._v(a,b)) )

    def _v(self, a, b):
        r"""
        Construct an immutable vector ``(a, b)``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: v = L._v(1,2)
            sage: v
            (1, 2)
            sage: v.parent()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: hash(v) == hash(v)
            True
        """
        ret = self._V((a,b))
        ret.set_immutable()
        return ret

    def bracket_on_basis(self, i, j):
        r"""
        Return the bracket of basis elements indexed by ``i`` and ``j``,
        where ``i < j``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: v = L._v
            sage: L.bracket_on_basis(('K',2), ('t', v(3,-1)))
            0
            sage: L.bracket_on_basis(('K', 4), ('E', v(3,-1)))
            0
            sage: L.bracket_on_basis(('t', v(3,-1)), ('t', v(4,3)))
            0
            sage: L.bracket_on_basis(('t', v(3,-1)), ('E', v(4,3)))
            -13*t(7, 2)
            sage: L.bracket_on_basis(('t', v(2,2)),  ('E', v(1,1)))
            0
            sage: L.bracket_on_basis(('t', v(3,-1)), ('E', v(-3,1)))
            3*K1 - K2
            sage: L.bracket_on_basis(('E', v(3,-1)), ('E', v(4,3)))
            -13*E(7, 2)
            sage: L.bracket_on_basis(('E', v(2,2)),  ('E', v(1,1)))
            0
            sage: L.bracket_on_basis(('E', v(3,-1)), ('E', v(-3,1)))
            3*K3 - K4
        """
        if i[0] == 'K' or j[0] == 'K':
            return self.zero()

        if i[0] == 't':
            if j[0] == 't':
                return self.zero()
            k = ('t', i[1] + j[1])
            if not k[1]: # == 0
                d = {('K',1): i[1][0], ('K',2): i[1][1]}     # Kronecker delta summand
            else:
                k[1].set_immutable()
                d = {k: j[1][0]*i[1][1] - j[1][1]*i[1][0]}   # determinant summand
            return self._from_dict(d)

        # else i[0] == 'E'
        k = ('E', i[1] + j[1])
        if not k[1]: # == 0
            d = {('K',3): i[1][0], ('K',4): i[1][1]}      # Kronecker delta summand
        else:
            k[1].set_immutable()
            d = {k: (j[1][0]*i[1][1] - j[1][1]*i[1][0])}  # determinant summand
        return self._from_dict(d)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L.an_element()
            K3 - 1/2*t(-1, 3) + E(1, -3) + E(2, 2)
        """
        d = self.monomial
        v = self._v
        return (
                 d( ('E',v(1,-3)) )
                 - self.base_ring().an_element() * d( ('t',v(-1,3)) )
                 + d( ('E',v(2,2)) )
                 + d( ('K',3) )
                )

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.RankTwoHeisenbergVirasoro(QQ)
            sage: L.some_elements()
            [E(1, 1), E(-2, -2), E(0, 1),
             t(1, 1), t(4, -1), t(2, 3),
             K2, K4,
             K3 - 1/2*t(-1, 3) + E(1, -3) + E(2, 2)]
        """
        d = self.monomial
        v = self._v
        return [d( ('E',v(1,1)) ), d( ('E',v(-2,-2)) ), d( ('E',v(0,1)) ),
                d( ('t',v(1,1)) ), d( ('t',v(4,-1)) ), d( ('t',v(2,3)) ),
                d( ('K',2) ), d( ('K',4) ), self.an_element()]

    class Element(LieAlgebraElement):
        pass
