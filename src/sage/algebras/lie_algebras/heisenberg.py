"""
Heisenberg Algebras

AUTHORS:

- Travis Scrimshaw (2013-08-13): Initial version
"""

# ***************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.indexed_generators import IndexedGenerators

from sage.algebras.lie_algebras.lie_algebra import (LieAlgebraFromAssociative,
                                                    LieAlgebraWithGenerators)
from sage.algebras.lie_algebras.lie_algebra_element import (LieAlgebraElement,
                                                            LieAlgebraMatrixWrapper)
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.cartesian_product import cartesian_product
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer import Integer
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.sets.set import Set


class HeisenbergAlgebra_abstract(IndexedGenerators):
    """
    The common methods for the (non-matrix) Heisenberg algebras.
    """
    def __init__(self, I):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo) # indirect doctest
        """
        IndexedGenerators.__init__(self, I, prefix='', bracket=False,
                                   latex_bracket=False, string_quotes=False)

    def p(self, i):
        """
        The generator `p_i` of the Heisenberg algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.p(2)
            p2
        """
        return self.element_class(self, {'p%i' % i: self.base_ring().one()})

    def q(self, i):
        """
        The generator `q_i` of the Heisenberg algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.q(2)
            q2
        """
        return self.element_class(self, {'q%i' % i: self.base_ring().one()})

    def z(self):
        """
        Return the basis element `z` of the Heisenberg algebra.

        The element `z` spans the center of the Heisenberg algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.z()
            z
        """
        return self.element_class(self, {'z': self.base_ring().one()})

    def bracket_on_basis(self, x, y):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``
        where ``x < y``.

        The basis of a Heisenberg algebra is ordered in such a way that
        the `p_i` come first, the `q_i` come next, and the `z` comes last.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: p1 = ('p', 1)
            sage: q1 = ('q', 1)
            sage: H.bracket_on_basis(p1, q1)
            z
        """
        if y == 'z': # No need to test for x == 'z' since x < y is assumed.
            return self.zero()
        if x[0] == 'p' and y[0] == 'q' and x[1] == y[1]:
            return self.z()
        return self.zero()

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: H._repr_term('p1')
            'p1'
            sage: H._repr_term('z')
            'z'
        """
        return m

    def _ascii_art_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: H._ascii_art_term('p1')
            p1
            sage: H._ascii_art_term('z')
            z
            sage: ascii_art(sum(i * b for i, b in enumerate(H.basis())))
            p2 + 2*p3 + 3*q1 + 4*q2 + 5*q3 + 6*z
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(m)

    def _latex_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 10)
            sage: H._latex_term('p1')
            'p_{1}'
            sage: H._latex_term('z')
            'z'
            sage: latex(H.p(10))
            p_{10}
        """
        if len(m) == 1:
            return m
        return "%s_{%s}" % (m[0], m[1:]) # else it is of length at least 2

    def _unicode_art_term(self, m):
        r"""
        Return a unicode art representation of the term indexed by ``m``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 10)
            sage: H._unicode_art_term('p1')
            p₁
            sage: H._unicode_art_term('z')
            z
            sage: unicode_art(H.p(10))
            p₁₀
        """
        from sage.typeset.unicode_art import unicode_art, unicode_subscript
        if len(m) == 1:
            return unicode_art(m)
        return unicode_art(str(m[0]) + unicode_subscript(m[1:])) # else it is of length at least 2

    def step(self):
        r"""
        Return the nilpotency step of ``self``.

        EXAMPLES::

            sage: h = lie_algebras.Heisenberg(ZZ, 10)
            sage: h.step()
            2

            sage: h = lie_algebras.Heisenberg(ZZ, oo)
            sage: h.step()
            2
        """
        return Integer(2)

    class Element(LieAlgebraElement):
        pass


class HeisenbergAlgebra_fd:
    """
    Common methods for finite-dimensional Heisenberg algebras.
    """
    def __init__(self, n):
        """
        Initialize ``self``.

        INPUT:

        - ``n`` -- the rank

        TESTS::

            sage: H = lie_algebras.Heisenberg(QQ, 3) # indirect doctest
        """
        self._n = n

    def n(self):
        """
        Return the rank of the Heisenberg algebra ``self``.

        This is the ``n`` such that ``self`` is the `n`-th Heisenberg
        algebra. The dimension of this Heisenberg algebra is then
        `2n + 1`.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: H.n()
            3
            sage: H = lie_algebras.Heisenberg(QQ, 3, representation='matrix')
            sage: H.n()
            3
        """
        return self._n

    @cached_method
    def gens(self) -> tuple:
        """
        Return the Lie algebra generators of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 2)
            sage: H.gens()
            (p1, p2, q1, q2)
            sage: H = lie_algebras.Heisenberg(QQ, 0)
            sage: H.gens()
            (z,)
        """
        return tuple(self.lie_algebra_generators())

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 2)
            sage: H.gen(0)
            p1
            sage: H.gen(3)
            q2
        """
        return self.gens()[i]

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the Lie algebra generators of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: H.lie_algebra_generators()
            Finite family {'p1': p1, 'q1': q1}
            sage: H = lie_algebras.Heisenberg(QQ, 0)
            sage: H.lie_algebra_generators()
            Finite family {'z': z}
        """
        if self._n == 0:
            return Family(['z'], lambda i: self.z())
        k = ['p%s' % i for i in range(1, self._n+1)]
        k += ['q%s' % i for i in range(1, self._n+1)]
        d = {}
        for i in range(1, self._n+1):
            d['p%s' % i] = self.p(i)
            d['q%s' % i] = self.q(i)
        return Family(k, lambda i: d[i])

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: H.basis()
            Finite family {'p1': p1, 'q1': q1, 'z': z}
        """
        d = {}
        for i in range(1, self._n+1):
            d['p%s' % i] = self.p(i)
            d['q%s' % i] = self.q(i)
        d['z'] = self.z()
        return Family(self._indices, lambda i: d[i])

    def _coerce_map_from_(self, H):
        """
        Return the coercion map from ``H`` to ``self`` if one exists,
        otherwise return ``None``.

        EXAMPLES::

            sage: HB = lie_algebras.Heisenberg(QQ, 3)
            sage: HM = lie_algebras.Heisenberg(QQ, 3, representation='matrix')
            sage: HB.has_coerce_map_from(HM)
            True
            sage: HM.has_coerce_map_from(HB)
            True
            sage: HB(HM.p(2))
            p2
            sage: HM(-HB.q(3)) == -HM.q(3)
            True
            sage: HB(HM.z())
            z
            sage: HM(HB.z()) == HM.z()
            True
            sage: HQ = lie_algebras.Heisenberg(QQ, 2)
            sage: HB.has_coerce_map_from(HQ)
            True
            sage: HB(HQ.p(2))
            p2
            sage: HZ = lie_algebras.Heisenberg(ZZ, 2)
            sage: HB.has_coerce_map_from(HZ)
            True
            sage: HB(HZ.p(2))
            p2
            sage: HZ = lie_algebras.Heisenberg(ZZ, 2, representation='matrix')
            sage: HB.has_coerce_map_from(HZ)
            True
            sage: HB(HZ.p(2))
            p2
        """
        if isinstance(H, HeisenbergAlgebra_fd):
            if H._n <= self._n and self.base_ring().has_coerce_map_from(H.base_ring()):
                return H.module_morphism(lambda i: self.basis()[i], codomain=self)
            return None # Otherwise no coercion
        return super()._coerce_map_from_(H)


class HeisenbergAlgebra(HeisenbergAlgebra_fd, HeisenbergAlgebra_abstract,
                        LieAlgebraWithGenerators):
    r"""
    A Heisenberg algebra defined using structure coefficients.

    The `n`-th Heisenberg algebra (where `n` is a nonnegative
    integer or infinity) is the Lie algebra with basis
    `\{p_i\}_{1 \leq i \leq n} \cup \{q_i\}_{1 \leq i \leq n} \cup \{z\}`
    with the following relations:

    .. MATH::

        [p_i, q_j] = \delta_{ij} z, \quad [p_i, z] = [q_i, z] = [p_i, p_j]
        = [q_i, q_j] = 0.

    This Lie algebra is also known as the Heisenberg algebra of rank `n`.

    .. NOTE::

        The relations `[p_i, q_j] = \delta_{ij} z`, `[p_i, z] = 0`, and
        `[q_i, z] = 0` are known as canonical commutation relations. See
        :wikipedia:`Canonical_commutation_relations`.

    .. WARNING::

        The `n` in the above definition is called the "rank" of the
        Heisenberg algebra; it is not, however, a rank in any of the usual
        meanings that this word has in the theory of Lie algebras.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the rank of the Heisenberg algebra

    REFERENCES:

    - :wikipedia:`Heisenberg_algebra`

    EXAMPLES::

        sage: L = lie_algebras.Heisenberg(QQ, 2)
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 2)
            sage: TestSuite(L).run()
            sage: L = lie_algebras.Heisenberg(QQ, 0)  # not tested -- :issue:`18224`
            sage: TestSuite(L).run()
        """
        HeisenbergAlgebra_fd.__init__(self, n)
        names = tuple(['p%s' % i for i in range(1,n+1)]
                      + ['q%s' % i for i in range(1,n+1)]
                      + ['z'])
        LieAlgebraWithGenerators.__init__(self, R, names=names, index_set=names,
            category=LieAlgebras(R).Nilpotent().FiniteDimensional().WithBasis())
        HeisenbergAlgebra_abstract.__init__(self, names)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.Heisenberg(QQ, 3)
            Heisenberg algebra of rank 3 over Rational Field
        """
        return "Heisenberg algebra of rank {} over {}".format(self._n, self.base_ring())


class InfiniteHeisenbergAlgebra(HeisenbergAlgebra_abstract, LieAlgebraWithGenerators):
    r"""
    The infinite Heisenberg algebra.

    This is the Heisenberg algebra on an infinite number of generators. In
    other words, this is the Heisenberg algebra of rank `\infty`. See
    :class:`HeisenbergAlgebra` for more information.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: TestSuite(L).run()
            sage: L.p(1).bracket(L.q(1)) == L.z()
            True
            sage: L.q(1).bracket(L.p(1)) == -L.z()
            True
        """
        S = cartesian_product([PositiveIntegers(), ['p','q']])
        cat = LieAlgebras(R).Nilpotent().WithBasis()
        LieAlgebraWithGenerators.__init__(self, R, index_set=S, category=cat)
        HeisenbergAlgebra_abstract.__init__(self, S)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.Heisenberg(QQ, oo)
            Infinite Heisenberg algebra over Rational Field
        """
        return "Infinite Heisenberg algebra over {}".format(self.base_ring())

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L._an_element_()
            p2 + q2 - 1/2*q3 + z
        """
        c = self.base_ring().an_element()
        return self.p(2) + self.q(2) - c * self.q(3) + self.z()

    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.lie_algebra_generators()
            Lazy family (generator map(i))_{i in The Cartesian product of
                                            (Positive integers, {'p', 'q'})}
        """
        return Family(self._indices, lambda x: self.monomial(x[1] + str(x[0])),
                      name='generator map')

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.basis()
            Lazy family (basis map(i))_{i in Disjoint union of Family ({'z'},
             The Cartesian product of (Positive integers, {'p', 'q'}))}
            sage: L.basis()['z']
            z
            sage: L.basis()[(12, 'p')]
            p12
        """
        S = cartesian_product([PositiveIntegers(), ['p', 'q']])
        I = DisjointUnionEnumeratedSets([Set(['z']), S])

        def basis_elt(x):
            if isinstance(x, str):
                return self.monomial(x)
            return self.monomial(x[1] + str(x[0]))
        return Family(I, basis_elt, name="basis map")

    def _from_fd_on_basis(self, i):
        """
        Return the monomial in ``self`` corresponding to the
        basis element indexed by ``i``, where ``i`` is a basis index for
        a *finite-dimensional* Heisenberg algebra.

        This is used for coercion.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, oo)
            sage: H._from_fd_on_basis('p2')
            p2
            sage: H._from_fd_on_basis('q3')
            q3
            sage: H._from_fd_on_basis('z')
            z
        """
        if i == 'z':
            return self.z()
        if i[0] == 'p':
            return self.p(Integer(i[1:]))
        return self.q(Integer(i[1:]))

    def _coerce_map_from_(self, H):
        """
        Return the coercion map from ``H`` to ``self`` if one exists,
        otherwise return ``None``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, oo)
            sage: HZ = lie_algebras.Heisenberg(ZZ, oo)
            sage: phi = H.coerce_map_from(HZ)
            sage: phi(HZ.p(3)) == H.p(3)
            True
            sage: phi(HZ.p(3)).leading_coefficient().parent()
            Rational Field
            sage: HF = lie_algebras.Heisenberg(QQ, 3, representation='matrix')
            sage: H.has_coerce_map_from(HF)
            True
            sage: H(HF.p(2))
            p2
            sage: H(HF.z())
            z
            sage: HF = lie_algebras.Heisenberg(QQ, 3)
            sage: H.has_coerce_map_from(HF)
            True
            sage: H(HF.p(2))
            p2
            sage: H(HF.z())
            z
        """
        if isinstance(H, HeisenbergAlgebra_fd):
            if self.base_ring().has_coerce_map_from(H.base_ring()):
                return H.module_morphism(self._from_fd_on_basis, codomain=self)
            return None  # Otherwise no coercion
        if isinstance(H, InfiniteHeisenbergAlgebra):
            if self.base_ring().has_coerce_map_from(H.base_ring()):
                return lambda C,x: self._from_dict(x._monomial_coefficients, coerce=True)
            return None  # Otherwise no coercion
        return super()._coerce_map_from_(H)


#######################################################
# Finite rank Heisenberg algebra using matrices

class HeisenbergAlgebra_matrix(HeisenbergAlgebra_fd, LieAlgebraFromAssociative):
    r"""
    A Heisenberg algebra represented using matrices.

    The `n`-th Heisenberg algebra over `R` is a Lie algebra which is
    defined as the Lie algebra of the `(n+2) \times (n+2)`-matrices:

    .. MATH::

        \begin{bmatrix}
        0 & p^T & k \\
        0 & 0_n & q \\
        0 & 0 & 0
        \end{bmatrix}

    where `p, q \in R^n` and `0_n` in the `n \times n` zero matrix. It has
    a basis consisting of

    .. MATH::

        \begin{aligned}
        p_i & = \begin{bmatrix}
        0 & e_i^T & 0 \\
        0 & 0_n & 0 \\
        0 & 0 & 0
        \end{bmatrix} \qquad \text{for } 1 \leq i \leq n ,
        \\ q_i & = \begin{bmatrix}
        0 & 0 & 0 \\
        0 & 0_n & e_i \\
        0 & 0 & 0
        \end{bmatrix} \qquad \text{for } 1 \leq i \leq n ,
        \\ z & = \begin{bmatrix}
        0 & 0 & 1 \\
        0 & 0_n & 0 \\
        0 & 0 & 0
        \end{bmatrix},
        \end{aligned}

    where `\{e_i\}` is the standard basis of `R^n`. In other words, it has
    the basis `(p_1, p_2, \ldots, p_n, q_1, q_2, \ldots, q_n, z)`, where
    `p_i = E_{1, i+1}`, `q_i = E_{i+1, n+2}` and `z = E_{1, n+2}` are
    elementary matrices.

    This Lie algebra is isomorphic to the `n`-th Heisenberg algebra
    constructed in :class:`HeisenbergAlgebra`; the bases correspond to
    each other.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the nonnegative integer `n`

    EXAMPLES::

        sage: L = lie_algebras.Heisenberg(QQ, 1, representation='matrix')
        sage: p = L.p(1)
        sage: q = L.q(1)
        sage: z = L.bracket(p, q); z
        [0 0 1]
        [0 0 0]
        [0 0 0]
        sage: z == L.z()
        True
        sage: L.dimension()
        3

        sage: L = lie_algebras.Heisenberg(QQ, 2, representation='matrix')
        sage: sorted(dict(L.basis()).items())
        [(
              [0 1 0 0]
              [0 0 0 0]
              [0 0 0 0]
        'p1', [0 0 0 0]
        ),
         (
              [0 0 1 0]
              [0 0 0 0]
              [0 0 0 0]
        'p2', [0 0 0 0]
        ),
         (
              [0 0 0 0]
              [0 0 0 1]
              [0 0 0 0]
        'q1', [0 0 0 0]
        ),
         (
              [0 0 0 0]
              [0 0 0 0]
              [0 0 0 1]
        'q2', [0 0 0 0]
        ),
         (
             [0 0 0 1]
             [0 0 0 0]
             [0 0 0 0]
        'z', [0 0 0 0]
        )]

        sage: L = lie_algebras.Heisenberg(QQ, 0, representation='matrix')
        sage: sorted(dict(L.basis()).items())
        [(
             [0 1]
        'z', [0 0]
        )]
        sage: L.gens()
        (
        [0 1]
        [0 0]
        )
        sage: L.lie_algebra_generators()
        Finite family {'z': [0 1]
        [0 0]}
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 2, representation='matrix')
            sage: TestSuite(L).run()
        """
        HeisenbergAlgebra_fd.__init__(self, n)
        MS = MatrixSpace(R, n+2, sparse=True)
        one = R.one()
        p = tuple(MS({(0,i): one}) for i in range(1, n+1))
        q = tuple(MS({(i,n+1): one}) for i in range(1, n+1))
        z = (MS({(0,n+1): one}),)
        names = tuple('p%s' % i for i in range(1,n+1))
        names = names + tuple('q%s' % i for i in range(1,n+1)) + ('z',)
        cat = LieAlgebras(R).Nilpotent().FiniteDimensional().WithBasis()
        LieAlgebraFromAssociative.__init__(self, MS, p + q + z, names=names,
                                           index_set=names, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.Heisenberg(QQ, 3, representation='matrix')
            Heisenberg algebra of rank 3 over Rational Field
        """
        return "Heisenberg algebra of rank {} over {}".format(self._n, self.base_ring())

    def p(self, i):
        r"""
        Return the generator `p_i` of the Heisenberg algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation='matrix')
            sage: L.p(1)
            [0 1 0]
            [0 0 0]
            [0 0 0]
        """
        return self._gens['p%s' % i]

    def q(self, i):
        r"""
        Return the generator `q_i` of the Heisenberg algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation='matrix')
            sage: L.q(1)
            [0 0 0]
            [0 0 1]
            [0 0 0]
        """
        return self._gens['q%s' % i]

    def z(self):
        """
        Return the basis element `z` of the Heisenberg algebra.

        The element `z` spans the center of the Heisenberg algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation='matrix')
            sage: L.z()
            [0 0 1]
            [0 0 0]
            [0 0 0]
        """
        return self._gens['z']

    def step(self):
        r"""
        Return the nilpotency step of ``self``.

        EXAMPLES::

            sage: h = lie_algebras.Heisenberg(ZZ, 2, representation='matrix')
            sage: h.step()
            2
        """
        return Integer(2)

    class Element(LieAlgebraMatrixWrapper, LieAlgebraFromAssociative.Element):
        def monomial_coefficients(self, copy=True):
            """
            Return a dictionary whose keys are indices of basis elements in
            the support of ``self`` and whose values are the corresponding
            coefficients.

            INPUT:

            - ``copy`` -- ignored

            EXAMPLES::

                sage: L = lie_algebras.Heisenberg(QQ, 3, representation='matrix')
                sage: elt = L(Matrix(QQ, [[0, 1, 3, 0, 3], [0, 0, 0, 0, 0], [0, 0, 0, 0, -3],
                ....:                     [0, 0, 0, 0, 7], [0, 0, 0, 0, 0]]))
                sage: elt
                [ 0  1  3  0  3]
                [ 0  0  0  0  0]
                [ 0  0  0  0 -3]
                [ 0  0  0  0  7]
                [ 0  0  0  0  0]
                sage: sorted(elt.monomial_coefficients().items())
                [('p1', 1), ('p2', 3), ('q2', -3), ('q3', 7), ('z', 3)]
            """
            d = {}
            n = self.parent()._n
            for i, mon in enumerate(self.parent().basis().keys()):
                if i < n:
                    entry = self[0, i+1]
                elif i < 2 * n:
                    entry = self[i-n+1, n+1]
                else:
                    entry = self[0, n+1]
                if entry:
                    d[mon] = entry
            return d
