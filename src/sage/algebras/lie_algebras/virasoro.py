"""
Virasoro Algebra and Related Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.lie_algebras import LieAlgebras
from sage.rings.all import ZZ
from sage.sets.family import Family
from sage.sets.set import Set
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.indexed_generators import IndexedGenerators
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import (InfinitelyGeneratedLieAlgebra,
                                                    FinitelyGeneratedLieAlgebra)

class LieAlgebraRegularVectorFields(InfinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    The Lie algebra of regular vector fields on `\CC^{\times}`.

    This is the Lie algebra with basis `\{d_i\}_{i \in \ZZ}` and subject
    to the relations

    .. MATH::

        [d_i, d_j] = (j - i) d_{i+j}.

    This is also known as the Witt (Lie) algebra.

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
        cat = LieAlgebras(R).WithBasis()
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
            -4*d[0]
            sage: L.bracket_on_basis(2, 4)
            2*d[6]
            sage: L.bracket_on_basis(4, 4)
            0
        """
        return self.term(i + j, j - i)

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

        [d_i, d_j] = (j - i) d_{i+j},

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
            sage: L = lie_algebras.pwitt(Zmod(6), 6)
            sage: TestSuite(L).run()  # not tested -- universal envelope doesn't work
            sage: L._test_jacobi_identity()
        """
        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        FinitelyGeneratedLieAlgebra.__init__(self, R, index_set=range(p), category=cat)
        IndexedGenerators.__init__(self, range(p), prefix='d', bracket='[')
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
            d[0]
            sage: L.bracket_on_basis(3, 2)
            4*d[0]
            sage: L.bracket_on_basis(2, 2)
            0
            sage: L.bracket_on_basis(1, 3)
            2*d[4]
        """
        return self.term((i + j) % self._p, j - i)

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

        [d_i, d_j] = (j - i) d_{i+j} + \frac{1}{12}(j^3 - j) \delta_{i,-j} c

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
        cat = LieAlgebras(R).WithBasis()
        InfinitelyGeneratedLieAlgebra.__init__(self, R, index_set=ZZ, category=cat)
        IndexedGenerators.__init__(self, ZZ, prefix='d', bracket='[',
                                   sorting_key=_basis_key)

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

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.VirasoroAlgebra(QQ)
            The Virasoro algebra over Rational Field
        """
        return "The Virasoro algebra over {}".format(self.base_ring())

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
            -4*d[0] - 1/2*c
        """
        if i == 'c' or j == 'c':
            return self.zero()
        ret = self._from_dict({i + j: j-i})
        R = self.base_ring()
        if i == -j:
            ret += R(j ** 3 - j) / R(12) * self.c()
        return ret

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

    class Element(LieAlgebraElement):
        pass

