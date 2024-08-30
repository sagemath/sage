"""
Virasoro Algebra and Related Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.modules import Modules
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.sets.family import Family
from sage.sets.set import Set
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.indexed_generators import IndexedGenerators
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import (InfinitelyGeneratedLieAlgebra,
                                                    FinitelyGeneratedLieAlgebra)
from sage.combinat.free_module import CombinatorialFreeModule


class LieAlgebraRegularVectorFields(InfinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    The Lie algebra of regular vector fields on `\CC^{\times}`.

    This is the Lie algebra with basis `\{d_i\}_{i \in \ZZ}` and subject
    to the relations

    .. MATH::

        [d_i, d_j] = (i - j) d_{i+j}.

    This is also known as the Witt (Lie) algebra.

    .. NOTE::

        This differs from some conventions (e.g., [Ka1990]_), where
        we have `d'_i \mapsto -d_i`.

    REFERENCES:

    - :wikipedia:`Witt_algebra`

    .. SEEALSO::

        :class:`WittLieAlgebra_charp`
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: TestSuite(L).run()
        """
        cat = LieAlgebras(R).WithBasis().Graded()
        InfinitelyGeneratedLieAlgebra.__init__(self, R, index_set=ZZ, category=cat)
        IndexedGenerators.__init__(self, ZZ, prefix='d', bracket='[')

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.regular_vector_fields(QQ)
            The Lie algebra of regular vector fields over Rational Field
        """
        return "The Lie algebra of regular vector fields over {}".format(self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.regular_vector_fields(QQ)
            sage: latex(g)
            \mathcal{W}_{\Bold{Q}}
        """
        from sage.misc.latex import latex
        return r"\mathcal{{W}}_{{{}}}".format(latex(self.base_ring()))

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: L.lie_algebra_generators()
            Lazy family (generator map(i))_{i in Integer Ring}
        """
        return Family(self._indices, self.monomial, name='generator map')

    def bracket_on_basis(self, i, j):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``
        where ``x < y``.

        (This particular implementation actually does not require
        ``x < y``.)

        EXAMPLES::

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: L.bracket_on_basis(2, -2)
            4*d[0]
            sage: L.bracket_on_basis(2, 4)
            -2*d[6]
            sage: L.bracket_on_basis(4, 4)
            0
        """
        return self.term(i + j, i - j)

    def degree_on_basis(self, i):
        r"""
        Return the degree of the basis element indexed by ``i``,
        which is ``i``.

        EXAMPLES::

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: L.degree_on_basis(2)
            2
        """
        return i

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: L.an_element()
            d[-1] + d[0] - 3*d[1]
        """
        return self.monomial(0) - 3*self.monomial(1) + self.monomial(-1)

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.regular_vector_fields(QQ)
            sage: L.some_elements()
            [d[0], d[2], d[-2], d[-1] + d[0] - 3*d[1]]
        """
        return [self.monomial(0), self.monomial(2), self.monomial(-2), self.an_element()]

    class Element(LieAlgebraElement):
        pass


class WittLieAlgebra_charp(FinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    The `p`-Witt Lie algebra over a ring `R` in which
    `p \cdot 1_R = 0`.

    Let `R` be a ring and `p` be a positive integer such that
    `p \cdot 1_R = 0`. The `p`-Witt Lie algebra over `R` is
    the Lie algebra with basis `\{d_0, d_1, \ldots, d_{p-1}\}`
    and subject to the relations

    .. MATH::

        [d_i, d_j] = (i - j) d_{i+j},

    where the `i+j` on the right hand side is identified with its
    remainder modulo `p`.

    .. SEEALSO::

        :class:`LieAlgebraRegularVectorFields`
    """
    def __init__(self, R, p):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.pwitt(GF(5), 5); L
            The 5-Witt Lie algebra over Finite Field of size 5
            sage: TestSuite(L).run()

        We skip the grading test as we need to be able to echelonize a
        matrix over the base ring as part of the test::

            sage: L = lie_algebras.pwitt(Zmod(6), 6)
            sage: TestSuite(L).run(skip='_test_grading')
        """
        if R(p) != 0:
            raise ValueError("{} is not 0 in {}".format(p, R))
        cat = LieAlgebras(R).FiniteDimensional().WithBasis().Graded()
        FinitelyGeneratedLieAlgebra.__init__(self, R, index_set=list(range(p)),
                                             category=cat)
        IndexedGenerators.__init__(self, list(range(p)), prefix='d',
                                   bracket='[')
        self._p = p

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.pwitt(Zmod(5), 5)
            The 5-Witt Lie algebra over Ring of integers modulo 5
            sage: lie_algebras.pwitt(Zmod(5), 15)
            The 15-Witt Lie algebra over Ring of integers modulo 5
        """
        return "The {}-Witt Lie algebra over {}".format(self._p, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.pwitt(GF(3), 15)
            sage: latex(g)
            \mathcal{W}(15)_{\Bold{F}_{3}}
        """
        from sage.misc.latex import latex
        return r"\mathcal{{W}}({})_{{{}}}".format(latex(self._p), latex(self.base_ring()))

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.pwitt(Zmod(5), 5)
            sage: L.lie_algebra_generators()
            Finite family {0: d[0], 1: d[1], 2: d[2], 3: d[3], 4: d[4]}
        """
        return Family(self._indices, self.monomial, name='generator map')

    def bracket_on_basis(self, i, j):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``
        where ``x < y``.

        (This particular implementation actually does not require
        ``x < y``.)

        EXAMPLES::

            sage: L = lie_algebras.pwitt(Zmod(5), 5)
            sage: L.bracket_on_basis(2, 3)
            4*d[0]
            sage: L.bracket_on_basis(3, 2)
            d[0]
            sage: L.bracket_on_basis(2, 2)
            0
            sage: L.bracket_on_basis(1, 3)
            3*d[4]
        """
        return self.term((i + j) % self._p, i - j)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.pwitt(Zmod(5), 5)
            sage: L.an_element()
            d[0] + 2*d[1] + d[4]
        """
        return self.monomial(0) - 3*self.monomial(1 % self._p) + self.monomial((-1) % self._p)

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.pwitt(Zmod(5), 5)
            sage: L.some_elements()
            [d[0], d[2], d[3], d[0] + 2*d[1] + d[4]]
        """
        return [self.monomial(0), self.monomial(2 % self._p),
                self.monomial((-2) % self._p),
                self.an_element()]

    def degree_on_basis(self, i):
        r"""
        Return the degree of the basis element indexed by ``i``,
        which is ``i`` mod `p`.

        EXAMPLES::

            sage: L = lie_algebras.pwitt(Zmod(5), 5)
            sage: L.degree_on_basis(7)
            2
            sage: L.degree_on_basis(2).parent()
            Ring of integers modulo 5
        """
        return IntegerModRing(self._p)(i)

    class Element(LieAlgebraElement):
        pass


def _basis_key(x):
    """
    Helper function that generates a key for the basis elements
    of the Virasoro algebra.

    EXAMPLES::

        sage: from sage.algebras.lie_algebras.virasoro import _basis_key
        sage: _basis_key('c')
        +Infinity
        sage: _basis_key(2)
        2
    """
    if x == 'c':
        from sage.rings.infinity import infinity
        return infinity
    return x


class VirasoroAlgebra(InfinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    The Virasoro algebra.

    This is the Lie algebra with basis `\{d_i\}_{i \in \ZZ} \cup \{c\}`
    and subject to the relations

    .. MATH::

        [d_i, d_j] = (i - j) d_{i+j} + \frac{1}{12}(i^3 - i) \delta_{i,-j} c

    and

    .. MATH::

        [d_i, c] = 0.

    (Here, it is assumed that the base ring `R` has `2` invertible.)

    This is the universal central extension `\widetilde{\mathfrak{d}}` of
    the Lie algebra `\mathfrak{d}` of
    :class:`regular vector fields <LieAlgebraRegularVectorFields>`
    on `\CC^{\times}`.

    EXAMPLES::

        sage: d = lie_algebras.VirasoroAlgebra(QQ)

    REFERENCES:

    - :wikipedia:`Virasoro_algebra`
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: TestSuite(d).run()
        """
        cat = LieAlgebras(R).WithBasis().Graded()
        InfinitelyGeneratedLieAlgebra.__init__(self, R, index_set=ZZ, category=cat)
        IndexedGenerators.__init__(self, ZZ, prefix='d', bracket='[',
                                   sorting_key=_basis_key)

    def _basis_key(self, m):
        """
        Return a key for sorting for the index ``m``.

        TESTS::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d._basis_key(3)
            3
            sage: d._basis_key('c')
            +Infinity
            sage: d._basis_key(4) < d._basis_key('c')
            True
        """
        return _basis_key(m)

    def _repr_term(self, m):
        """
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d._repr_term('c')
            'c'
            sage: d._repr_term(2)
            'd[2]'
        """
        if isinstance(m, str):
            return m
        return IndexedGenerators._repr_generator(self, m)

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the term indexed by ``m``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d._latex_term('c')
            'c'
            sage: d._latex_term(2)
            'd_{2}'
            sage: d._latex_term(-13)
            'd_{-13}'
        """
        if isinstance(m, str):
            return m
        return IndexedGenerators._latex_generator(self, m)

    def _unicode_art_term(self, m):
        r"""
        Return a unicode art representation of the term indexed by ``m``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d._unicode_art_term('c')
            c
            sage: d._unicode_art_term(2)
            d₂
            sage: d._unicode_art_term(-13)
            d₋₁₃
        """
        from sage.typeset.unicode_art import unicode_art, unicode_subscript
        if isinstance(m, str):
            return unicode_art(m)
        return unicode_art('d' + unicode_subscript(m))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.VirasoroAlgebra(QQ)
            The Virasoro algebra over Rational Field
        """
        return "The Virasoro algebra over {}".format(self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.VirasoroAlgebra(QQ)
            sage: latex(g)
            \mathcal{V}_{\Bold{Q}}
        """
        from sage.misc.latex import latex
        return r"\mathcal{{V}}_{{{}}}".format(latex(self.base_ring()))

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d.lie_algebra_generators()
            Lazy family (generator map(i))_{i in Integer Ring}
        """
        return Family(self._indices, self.monomial, name='generator map')

    @cached_method
    def basis(self):
        """
        Return a basis of ``self``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: B = d.basis(); B
            Lazy family (basis map(i))_{i in Disjoint union of
                                        Family ({'c'}, Integer Ring)}
            sage: B['c']
            c
            sage: B[3]
            d[3]
            sage: B[-15]
            d[-15]
        """
        I = DisjointUnionEnumeratedSets([Set(['c']), ZZ])
        return Family(I, self.monomial, name='basis map')

    def d(self, i):
        """
        Return the element `d_i` in ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: L.d(2)
            d[2]
        """
        return self.monomial(i)

    def c(self):
        """
        The central element `c` in ``self``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d.c()
            c
        """
        return self.monomial('c')

    def bracket_on_basis(self, i, j):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``
        where ``x < y``.

        (This particular implementation actually does not require
        ``x < y``.)

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d.bracket_on_basis('c', 2)
            0
            sage: d.bracket_on_basis(2, -2)
            4*d[0] + 1/2*c
        """
        if i == 'c' or j == 'c':
            return self.zero()
        ret = self._from_dict({i + j: i-j})
        R = self.base_ring()
        if i == -j:
            ret += R(i ** 3 - i) / R(12) * self.c()
        return ret

    def degree_on_basis(self, i):
        r"""
        Return the degree of the basis element indexed by ``i``,
        which is ``i`` and `0` for ``'c'``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d.degree_on_basis(2)
            2
            sage: d.c().degree()
            0
            sage: (d.c() + d.basis()[0]).is_homogeneous()
            True
        """
        if i == 'c':
            return ZZ.zero()
        return i

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d.an_element()
            d[-1] + d[0] - 1/2*d[1] + c
        """
        d = self.monomial
        return d(0) - self.base_ring().an_element()*d(1) + d(-1) + d('c')

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: d = lie_algebras.VirasoroAlgebra(QQ)
            sage: d.some_elements()
            [d[0], d[2], d[-2], c, d[-1] + d[0] - 1/2*d[1] + c]
        """
        d = self.monomial
        return [d(0), d(2), d(-2), d('c'), self.an_element()]

    def chargeless_representation(self, a, b):
        """
        Return the chargeless representation of ``self`` with
        parameters ``a`` and ``b``.

        .. SEEALSO::

            :class:`~sage.algebras.lie_algebras.virasoro.ChargelessRepresentation`

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: L.chargeless_representation(3, 2)
            Chargeless representation (3, 2) of
             The Virasoro algebra over Rational Field
        """
        return ChargelessRepresentation(self, a, b)

    def verma_module(self, c, h):
        """
        Return the Verma module with central charge ``c`` and
        conformal (or highest) weight ``h``.

        .. SEEALSO::

            :class:`~sage.algebras.lie_algebras.virasoro.VermaModule`

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: L.verma_module(3, 2)
            Verma module with charge 3 and conformal weight 2 of
             The Virasoro algebra over Rational Field
        """
        return VermaModule(self, c, h)

    class Element(LieAlgebraElement):
        pass

#####################################################################
## Representations


class ChargelessRepresentation(CombinatorialFreeModule):
    r"""
    A chargeless representation of the Virasoro algebra.

    Let `L` be the Virasoro algebra over the field `F` of
    characteristic `0`. For `\alpha, \beta \in R`, we denote `V_{a,b}`
    as the `(a, b)`-*chargeless representation* of `L`, which is the
    `F`-span of `\{v_k \mid k \in \ZZ\}` with `L` action

    .. MATH::

        \begin{aligned}
        d_n \cdot v_k & = (a n + b - k) v_{n+k},
        \\ c \cdot v_k & = 0,
        \end{aligned}

    This comes from the action of `d_n = -t^{n+1} \frac{d}{dt}` on
    `F[t, t^{-1}]` (recall that `L` is the central extension of the
    :class:`algebra of derivations <LieAlgebraRegularVectorFields>`
    of `F[t, t^{-1}]`), where

    .. MATH::

        V_{a,b} = F[t, t^{-1}] t^{a-b} (dt)^{-a}

    and `v_k = t^{a-b+k} (dz)^{-a}`.

    The chargeless representations are either irreducible or
    contains exactly two simple subquotients, one of which is the
    trivial representation and the other is `F[t, t^{-1}] / F`.
    The non-trivial simple subquotients are called the
    *intermediate series*.

    The module `V_{a,b}` is irreducible if and only if
    `a \neq 0, -1` or `b \notin \ZZ`. When `a = 0` and `b \in \ZZ`,
    then there exists a subrepresentation isomorphic to the trivial
    representation. If `a = -1` and `b \in \ZZ`, then there exists
    a subrepresentation `V` such that `V_{a,b} / V` is isomorphic
    to `K \frac{dt}{t}` and `V` is irreducible.

    In characteristic `p`, the non-trivial simple subquotient
    is isomorphic to `F[t, t^{-1}] / F[t^p, t^{-p}]`. For
    `p \neq 2,3`, then the action is given as above.

    EXAMPLES:

    We first construct the irreducible `V_{1/2, 3/4}` and do some
    basic computations::

        sage: L = lie_algebras.VirasoroAlgebra(QQ)
        sage: M = L.chargeless_representation(1/2, 3/4)
        sage: d = L.basis()
        sage: v = M.basis()
        sage: d[3] * v[2]
        1/4*v[5]
        sage: d[3] * v[-1]
        13/4*v[2]
        sage: (d[3] - d[-2]) * (v[-1] + 1/2*v[0] - v[4])
        -3/4*v[-3] + 1/8*v[-2] - v[2] + 9/8*v[3] + 7/4*v[7]

    We construct the reducible `V_{0,2}` and the trivial
    subrepresentation given by the span of `v_2`. We verify
    this for `\{d_i \mid -10 \leq i < 10\}`::

        sage: M = L.chargeless_representation(0, 2)
        sage: v = M.basis()
        sage: all(d[i] * v[2] == M.zero() for i in range(-10, 10))
        True

    REFERENCES:

    - [Mat1992]_
    - [IK2010]_
    """
    def __init__(self, V, a, b):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.chargeless_representation(1/2, 3/4)
            sage: TestSuite(M).run()
        """
        self._a = a
        self._b = b
        self._V = V
        R = V.base_ring()
        if R.characteristic() in [2, 3]:
            raise NotImplementedError("not implemented for characteristic 2,3")
        cat = Modules(R).WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, R, ZZ, prefix='v', category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: L.chargeless_representation(1/2, 3/4)
            Chargeless representation (1/2, 3/4) of
             The Virasoro algebra over Rational Field
        """
        return "Chargeless representation ({}, {}) of {}".format(
                    self._a, self._b, self._V)

    def parameters(self):
        """
        Return the parameters `(a, b)` of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.chargeless_representation(1/2, 3/4)
            sage: M.parameters()
            (1/2, 3/4)
        """
        return (self._a, self._b)

    def virasoro_algebra(self):
        """
        Return the Virasoro algebra ``self`` is a representation of.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.chargeless_representation(1/2, 3/4)
            sage: M.virasoro_algebra() is L
            True
        """
        return self._V

    def degree_on_basis(self, i):
        r"""
        Return the degree of the basis element indexed by ``i``,
        which is `i`.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.chargeless_representation(1/2, 3/4)
            sage: M.degree_on_basis(-3)
            -3
        """
        return i

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: L = lie_algebras.VirasoroAlgebra(QQ)
                sage: d = L.basis()
                sage: M = L.chargeless_representation(1/2, 3/4)
                sage: x = d[-5] * M.an_element() + M.basis()[10]; x
                -9/4*v[-6] - 7/4*v[-5] - 33/4*v[-4] + v[10]
                sage: d[2] * x
                -279/16*v[-4] - 189/16*v[-3] - 759/16*v[-2] - 33/4*v[12]

                sage: v = M.basis()
                sage: all(d[i]*(d[j]*v[k]) - d[j]*(d[i]*v[k]) == d[i].bracket(d[j])*v[k]
                ....:     for i in range(-5, 5) for j in range(-5, 5) for k in range(-5, 5))
                True
            """
            P = self.parent()
            # We implement only a left action
            if not self_on_left and scalar in P._V:
                scalar = P._V(scalar)
                return P.sum_of_terms((n+k, (P._a * n + P._b - k) * cv * cm)
                                      for n,cv in scalar.monomial_coefficients(copy=False).items() if n != 'c'
                                      for k,cm in self.monomial_coefficients(copy=False).items())
            return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)

        _rmul_ = _lmul_ = _acted_upon_


class VermaModule(CombinatorialFreeModule):
    r"""
    A Verma module of the Virasoro algebra.

    The Virasoro algebra admits a triangular decomposition

    .. MATH::

        V_- \oplus R d_0 \oplus R \hat{c} \oplus V_+,

    where `V_-` (resp. `V_+`) is the span of `\{d_i \mid i < 0\}`
    (resp. `\{d_i \mid i > 0\}`). We can construct the *Verma module*
    `M_{c,h}` as the induced representation of the `R d_0 \oplus
    R \hat{c} \oplus V_+` representation `R_{c,H} = Rv`, where

    .. MATH::

        V_+ v = 0, \qquad \hat{c} v = c v, \qquad d_0 v = h v.

    Therefore, we have a basis of `M_{c,h}`

    .. MATH::

        \{ L_{i_1} \cdots L_{i_k} v \mid i_1 \leq \cdots \leq i_k < 0 \}.

    Moreover, the Verma modules are the free objects in the category of
    highest weight representations of `V` and are indecomposable.
    The Verma module `M_{c,h}` is irreducible for generic values of `c`
    and `h` and when it is reducible, the quotient by the maximal
    submodule is the unique irreducible highest weight representation
    `V_{c,h}`.

    EXAMPLES:

    We construct a Verma module and do some basic computations::

        sage: L = lie_algebras.VirasoroAlgebra(QQ)
        sage: M = L.verma_module(3, 0)
        sage: d = L.basis()
        sage: v = M.highest_weight_vector()
        sage: d[3] * v
        0
        sage: d[-3] * v
        d[-3]*v
        sage: d[-1] * (d[-3] * v)
        2*d[-4]*v + d[-3]*d[-1]*v
        sage: d[2] * (d[-1] * (d[-3] * v))
        12*d[-2]*v + 5*d[-1]*d[-1]*v

    We verify that `d_{-1} v` is a singular vector for
    `\{d_i \mid 1 \leq i < 20\}`::

        sage: w = M.basis()[-1]; w
        d[-1]*v
        sage: all(d[i] * w == M.zero() for i in range(1,20))
        True

    We also verify a singular vector for `V_{-2,1}`::

        sage: M = L.verma_module(-2, 1)
        sage: B = M.basis()
        sage: w = B[-1,-1] - 2 * B[-2]
        sage: d = L.basis()
        sage: all(d[i] * w == M.zero() for i in range(1,20))
        True

    REFERENCES:

    - :wikipedia:`Virasoro_algebra#Representation_theory`
    """
    @staticmethod
    def __classcall_private__(cls, V, c, h):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(3, 1/2)
            sage: M2 = L.verma_module(int(3), 1/2)
            sage: M is M2
            True
        """
        R = V.base_ring()
        return super().__classcall__(cls, V, R(c), R(h))

    @staticmethod
    def _partition_to_neg_tuple(x):
        """
        Helper function to convert a partition to an increasing
        sequence of negative numbers.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.virasoro import VermaModule
            sage: VermaModule._partition_to_neg_tuple([3,2,2,1])
            (-3, -2, -2, -1)
        """
        # The entries of the partition are likely ints, but we need to
        #   make sure they are Integers.
        return tuple([ZZ(-i) for i in x])

    def __init__(self, V, c, h):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(3, 1/2)
            sage: TestSuite(M).run()
        """
        self._c = c
        self._h = h
        self._V = V
        from sage.combinat.partition import _Partitions
        indices = _Partitions.map(VermaModule._partition_to_neg_tuple)
        R = V.base_ring()
        cat = Modules(R).WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, R, indices, prefix='v', category=cat)

    def _repr_term(self, k):
        """
        Return a string representation for the term indexed by ``k``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(1, -2)
            sage: M._repr_term((-3,-2,-2,-1))
            'd[-3]*d[-2]*d[-2]*d[-1]*v'
        """
        if not k:
            return 'v'
        d = self._V.basis()
        return '*'.join(repr(d[i]) for i in k) + '*v'

    def _latex_term(self, k):
        """
        Return a latex representation for the term indexed by ``k``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(1, -2)
            sage: M._latex_term((-3,-2,-2,-1))
            'd_{-3} d_{-2} d_{-2} d_{-1} v'
        """
        if not k:
            return 'v'
        d = self._V.basis()
        from sage.misc.latex import latex
        return ' '.join(latex(d[i]) for i in k) + ' v'

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(3, 0)
            sage: M
            Verma module with charge 3 and conformal weight 0 of
             The Virasoro algebra over Rational Field
        """
        return "Verma module with charge {} and conformal weight {} of {}".format(
                    self._c, self._h, self._V)

    def _monomial(self, index):
        """
        TESTS::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(3, 0)
            sage: v = M.basis()
            sage: v[-3]  # indirect doctest
            d[-3]*v
            sage: v[-3,-2,-2]  # indirect doctest
            d[-3]*d[-2]*d[-2]*v
        """
        if index in ZZ:
            if index >= 0:
                raise ValueError("sequence must have nonpositive entries")
            index = (index,)
        return super()._monomial(index)

    def central_charge(self):
        """
        Return the central charge of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(3, 0)
            sage: M.central_charge()
            3
        """
        return self._c

    def conformal_weight(self):
        """
        Return the conformal weight of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(3, 0)
            sage: M.conformal_weight()
            3
        """
        return self._c

    def virasoro_algebra(self):
        """
        Return the Virasoro algebra ``self`` is a representation of.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(1/2, 3/4)
            sage: M.virasoro_algebra() is L
            True
        """
        return self._V

    @cached_method
    def highest_weight_vector(self):
        """
        Return the highest weight vector of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(-2/7, 3)
            sage: M.highest_weight_vector()
            v
        """
        return self.monomial(())

    def _d_action_on_basis(self, n, k):
        """
        Return the action of `d_n` on `v_k`.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(-2/7, 3)
            sage: M._d_action_on_basis(-3, ())
            d[-3]*v
            sage: M._d_action_on_basis(0, ())
            3*v
            sage: M._d_action_on_basis('c', ())
            -2/7*v
            sage: M._d_action_on_basis('c', (-4,-2,-2,-1))
            -2/7*d[-4]*d[-2]*d[-2]*d[-1]*v
            sage: M._d_action_on_basis(3, (-4,-2,-2,-1))
            7*d[-5]*d[-1]*v + 60*d[-4]*d[-2]*v + 15*d[-4]*d[-1]*d[-1]*v
             + 14*d[-3]*d[-2]*d[-1]*v + 7*d[-2]*d[-2]*d[-1]*d[-1]*v
            sage: M._d_action_on_basis(-1, (-4,-2,-2,-1))
            d[-9]*d[-1]*v + d[-5]*d[-4]*d[-1]*v + 3*d[-5]*d[-2]*d[-2]*d[-1]*v
             + 2*d[-4]*d[-3]*d[-2]*d[-1]*v + d[-4]*d[-2]*d[-2]*d[-1]*d[-1]*v
        """
        # c acts my multiplication by self._c on all elements
        if n == 'c':
            return self.term(k, self._c)

        # when k corresponds to the highest weight vector
        if not k:
            if n > 0:
                return self.zero()
            if n == 0:
                return self.term(k, self._h)
            return self.monomial((n,))

        # The basis are eigenvectors for d_0
        if n == 0:
            return self.term(k, self._h - sum(k))

        # We keep things in order
        if n <= k[0]:
            return self.monomial((n,) + k)

        # [L_n, L_m] v = L_n L_m v - L_m L_n v
        # L_n L_m v = L_m L_n v + [L_n, L_m] v
        d = self._V.basis()
        m = k[0]
        k = k[1:]
        # We need to explicitly call the action as this method is
        #   used in discovering the action
        return (self._d_action_on_basis(n, k)._acted_upon_(d[m], False)
                + self.monomial(k)._acted_upon_(d[n].bracket(d[m]), False))

    def degree_on_basis(self, d):
        r"""
        Return the degree of the basis element indexed by ``d``, which
        is the sum of the entries of ``d``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: M = L.verma_module(-2/7, 3)
            sage: M.degree_on_basis((-3,-3,-1))
            -7
        """
        return sum(d)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: L = lie_algebras.VirasoroAlgebra(QQ)
                sage: d = L.basis()
                sage: M = L.verma_module(1/2, 3/4)
                sage: x = d[-5] * M.an_element() + M.basis()[-10]; x
                d[-10]*v + 2*d[-5]*v + 3*d[-5]*d[-2]*v + 2*d[-5]*d[-1]*v
                sage: d[2] * x
                12*d[-8]*v + 39/4*d[-5]*v + 14*d[-3]*v + 21*d[-3]*d[-2]*v
                 + 14*d[-3]*d[-1]*v
                sage: v = M.highest_weight_vector()
                sage: d[2] * (d[-2] * v)
                13/4*v

                sage: it = iter(M.basis())
                sage: B = [next(it) for _ in range(10)]
                sage: all(d[i]*(d[j]*v) - d[j]*(d[i]*v) == d[i].bracket(d[j])*v
                ....:     for i in range(-5, 5) for j in range(-5, 5) for v in B)
                True
            """
            P = self.parent()
            # We implement only a left action
            if not self_on_left:
                S = scalar.parent()
                R = P.base_ring()
                if S is R or scalar in R:
                    scalar = R(scalar)
                    return P._from_dict({k: scalar*c for k,c in self._monomial_coefficients.items()})
                elif S is P._V or scalar in P._V:
                    scalar = P._V(scalar)
                    return P.linear_combination((P._d_action_on_basis(n, k), cv * cm)
                                                for n,cv in scalar.monomial_coefficients(copy=False).items()
                                                for k,cm in self._monomial_coefficients.items())
            return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)

        _rmul_ = _lmul_ = _acted_upon_
