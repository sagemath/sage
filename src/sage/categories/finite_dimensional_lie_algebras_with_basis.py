# sage_setup: distribution = sagemath-categories
r"""
Finite Dimensional Lie Algebras With Basis

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2013-2024 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.subobjects import SubobjectsCategory
from sage.sets.family import Family


def _ce_complex_key(self, M, d, s, n):
    """
    The key for caching the Chevalley-Eilenberg complex.

    TESTS::

        sage: from sage.categories.finite_dimensional_lie_algebras_with_basis import _ce_complex_key
        sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
        sage: _ce_complex_key(L, None, False, True, 5)
        (None, False, True)
        sage: f = ({x: Matrix([[1,0],[0,0]]), y: Matrix([[0,1],[0,0]])})
        sage: _ce_complex_key(L, f, False, True, 5)
        (Representation of Lie algebra on 2 generators (x, y) over Rational Field defined by:
                [1 0]
         x |--> [0 0]
                [0 1]
         y |--> [0 0],
         False,
         True)
    """
    if isinstance(M, dict):
        M = self.representation(M)
    return (M, d, s)


class FiniteDimensionalLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    Category of finite dimensional Lie algebras with a basis.

    .. TODO::

        Many of these tests should use non-abelian Lie algebras and need to
        be added after :issue:`16820`.
    """
    _base_category_class_and_axiom = (LieAlgebras.FiniteDimensional, "WithBasis")

    def example(self, n=3):
        """
        Return an example of a finite dimensional Lie algebra with basis as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis()
            sage: C.example()                                                           # needs sage.modules
            An example of a finite dimensional Lie algebra with basis:
             the 3-dimensional abelian Lie algebra over Rational Field

        Other dimensions can be specified as an optional argument::

            sage: C.example(5)                                                          # needs sage.modules
            An example of a finite dimensional Lie algebra with basis:
             the 5-dimensional abelian Lie algebra over Rational Field
        """
        from sage.categories.examples.finite_dimensional_lie_algebras_with_basis import Example
        return Example(self.base_ring(), n)

    Nilpotent = LazyImport('sage.categories.finite_dimensional_nilpotent_lie_algebras_with_basis',
                           'FiniteDimensionalNilpotentLieAlgebrasWithBasis')

    class ParentMethods:
        @cached_method
        def _construct_UEA(self):
            r"""
            Construct the universal enveloping algebra of ``self``.

            EXAMPLES::

                sage: # needs sage.combinat sage.libs.singular sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: UEA = L._construct_UEA(); UEA
                Noncommutative Multivariate Polynomial Ring in b0, b1, b2
                 over Rational Field, nc-relations: {}
                sage: UEA.relations(add_commutative=True)
                {b1*b0: b0*b1, b2*b0: b0*b2, b2*b1: b1*b2}

            ::

                sage: # needs sage.combinat sage.libs.singular sage.modules
                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1},
                ....:                             ('y','z'): {'x':1},
                ....:                             ('z','x'):{'y':1}})
                sage: UEA = L._construct_UEA(); UEA
                Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
                 nc-relations: {...}
                sage: sorted(UEA.relations().items(), key=str)
                [(y*x, x*y - z), (z*x, x*z + y), (z*y, y*z - x)]

            Singular's ``nc_algebra`` does not work over `\ZZ/6\ZZ`,
            so we fallback to the PBW basis in this case::

                sage: L = lie_algebras.pwitt(Zmod(6), 6)                                # needs sage.combinat sage.modules
                sage: L._construct_UEA()                                                # needs sage.combinat sage.libs.singular sage.modules
                Universal enveloping algebra of
                 The 6-Witt Lie algebra over Ring of integers modulo 6
                 in the Poincare-Birkhoff-Witt basis

            Corner case for the trivial (0-dimensional) Lie algebra::

                sage: L.<a,b,c> = LieAlgebra(QQ, abelian=True)
                sage: I = L.product_space(L)
                sage: I._construct_UEA()
                Free Algebra on 0 generators () over Rational Field
            """
            from sage.algebras.free_algebra import FreeAlgebra

            # Create the UEA relations
            # We need to get names for the basis elements, not just the generators
            I = self._basis_ordering
            if not I:  # trivial Lie algebra
                return FreeAlgebra(self.base_ring(), [])
            try:
                names = [str(x) for x in I]

                def names_map(x):
                    return x
                F = FreeAlgebra(self.base_ring(), names)
            except ValueError:
                names = ['b{}'.format(i) for i in range(self.dimension())]
                self._UEA_names_map = {g: names[i] for i,g in enumerate(I)}
                names_map = self._UEA_names_map.__getitem__
                F = FreeAlgebra(self.base_ring(), names)
            # ``F`` is the free algebra over the basis of ``self``. The
            # universal enveloping algebra of ``self`` will be constructed
            # as a quotient of ``F``.
            d = F.gens_dict()
            rels = {}
            S = self.structure_coefficients(True)
            # Construct the map from indices to names of the UEA

            def get_var(g):
                return d[names_map(g)]
            # The function ``get_var`` sends an element of the basis of
            # ``self`` to the corresponding element of ``F``.
            for k in S.keys():
                g0 = get_var(k[0])
                g1 = get_var(k[1])
                if g0 < g1:
                    rels[g1*g0] = g0*g1 - F.sum(val*get_var(g) for g, val in S[k])
                else:
                    rels[g0*g1] = g1*g0 + F.sum(val*get_var(g) for g, val in S[k])
            try:
                return F.g_algebra(rels)
            except RuntimeError:
                # Something went wrong with the computation, so fallback to
                #   the generic PBW basis implementation
                return self.pbw_basis()

        @lazy_attribute
        def _basis_ordering(self):
            """
            Return the indices of the basis of ``self`` as a tuple in
            a fixed order.

            Override this attribute to get a specific ordering of the basis.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()     # needs sage.modules
                sage: L._basis_ordering                                                 # needs sage.modules
                (0, 1, 2)
            """
            return tuple(self.basis().keys())

        @lazy_attribute
        def _basis_key_inverse(self):
            """
            A dictionary for keys to their appropriate index given by
            ``self._basis_ordering``.

            EXAMPLES::

                sage: # needs sage.combinat sage.groups sage.modules
                sage: G = SymmetricGroup(3)
                sage: S = GroupAlgebra(G, QQ)
                sage: L = LieAlgebra(associative=S)
                sage: [L._basis_key_inverse[k] for k in L._basis_ordering]
                [0, 1, 2, 3, 4, 5]
            """
            return {k: i for i,k in enumerate(self._basis_ordering)}

        def _basis_key(self, x):
            """
            Return a key for sorting for the index ``x``.

            TESTS::

                sage: L = lie_algebras.three_dimensional_by_rank(QQ, 3,                 # needs sage.groups sage.modules
                ....:                                            names=['E','F','H'])
                sage: PBW = L.pbw_basis()                                               # needs sage.groups sage.modules
                sage: PBW._basis_key('E') < PBW._basis_key('H')                         # needs sage.groups sage.modules
                True

            ::

                sage: L = lie_algebras.sl(QQ, 2)                                        # needs sage.groups sage.modules
                sage: def neg_key(x):
                ....:     return -L.basis().keys().index(x)
                sage: PBW = L.pbw_basis(basis_key=neg_key)                              # needs sage.groups sage.modules
                sage: prod(PBW.gens())  # indirect doctest                              # needs sage.groups sage.modules
                PBW[-alpha[1]]*PBW[alphacheck[1]]*PBW[alpha[1]]
                 - 4*PBW[-alpha[1]]*PBW[alpha[1]] + PBW[alphacheck[1]]^2
                 - 2*PBW[alphacheck[1]]

            Check that :issue:`23266` is fixed::

                sage: # needs sage.groups sage.modules
                sage: sl2 = lie_algebras.sl(QQ, 2, 'matrix')
                sage: sl2.indices()
                {'e1', 'f1', 'h1'}
                sage: type(sl2.basis().keys())
                <class 'list'>
                sage: Usl2 = sl2.pbw_basis()
                sage: Usl2._basis_key(2)
                2
                sage: Usl2._basis_key(3)
                Traceback (most recent call last):
                ...
                KeyError: 3
            """
            return self._basis_key_inverse[x]

        def _dense_free_module(self, R=None):
            """
            Return a dense free module associated to ``self`` over ``R``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()     # needs sage.modules
                sage: L._dense_free_module()                                            # needs sage.modules
                Vector space of dimension 3 over Rational Field
            """
            if R is None:
                R = self.base_ring()
            from sage.modules.free_module import FreeModule
            return FreeModule(R, self.dimension())

        module = _dense_free_module

        def from_vector(self, v, order=None):
            """
            Return the element of ``self`` corresponding to the
            vector ``v`` in ``self.module()``.

            Implement this if you implement :meth:`module`; see the
            documentation of
            :meth:`sage.categories.lie_algebras.LieAlgebras.module`
            for how this is to be done.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()     # needs sage.modules
                sage: u = L.from_vector(vector(QQ, (1, 0, 0))); u                       # needs sage.modules
                (1, 0, 0)
                sage: parent(u) is L                                                    # needs sage.modules
                True
            """
            if order is None:
                order = self._basis_ordering
            B = self.basis()
            return self.sum(v[i] * B[k] for i,k in enumerate(order) if v[i] != 0)

        def killing_matrix(self, x, y):
            r"""
            Return the Killing matrix of ``x`` and ``y``, where ``x``
            and ``y`` are two elements of ``self``.

            The Killing matrix is defined as the matrix corresponding
            to the action of
            `\operatorname{ad}_x \circ \operatorname{ad}_y` in the
            basis of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()     # needs sage.modules
                sage: a, b, c = L.lie_algebra_generators()                              # needs sage.modules
                sage: L.killing_matrix(a, b)                                            # needs sage.modules
                [0 0 0]
                [0 0 0]
                [0 0 0]

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})                    # needs sage.combinat sage.modules
                sage: L.killing_matrix(y, x)                                            # needs sage.combinat sage.modules
                [ 0 -1]
                [ 0  0]
            """
            return x.adjoint_matrix() * y.adjoint_matrix()

        def killing_form(self, x, y):
            r"""
            Return the Killing form on ``x`` and ``y``, where ``x``
            and ``y`` are two elements of ``self``.

            The Killing form is defined as

            .. MATH::

                \langle x \mid y \rangle
                = \operatorname{tr}\left( \operatorname{ad}_x
                \circ \operatorname{ad}_y \right).

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()     # needs sage.modules
                sage: a, b, c = L.lie_algebra_generators()                              # needs sage.modules
                sage: L.killing_form(a, b)                                              # needs sage.modules
                0
            """
            return self.killing_matrix(x, y).trace()

        @cached_method
        def killing_form_matrix(self):
            """
            Return the matrix of the Killing form of ``self``.

            The rows and the columns of this matrix are indexed by the
            elements of the basis of ``self`` (in the order provided by
            :meth:`basis`).

            EXAMPLES::

                sage: # needs sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.killing_form_matrix()
                [0 0 0]
                [0 0 0]
                [0 0 0]
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example(0)
                sage: m = L.killing_form_matrix(); m
                []
                sage: parent(m)
                Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            """
            from sage.matrix.constructor import matrix

            B = self.basis()
            m = matrix(self.base_ring(),
                       [[self.killing_form(x, y) for x in B] for y in B])
            m.set_immutable()
            return m

        @cached_method
        def structure_coefficients(self, include_zeros=False):
            r"""
            Return the structure coefficients of ``self``.

            INPUT:

            - ``include_zeros`` -- boolean (default: ``False``); if ``True``,
              then include the `[x, y] = 0` pairs in the output

            OUTPUT:

            A dictionary whose keys are pairs of basis indices `(i, j)`
            with `i < j`, and whose values are the corresponding
            *elements* `[b_i, b_j]` in the Lie algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()     # needs sage.modules
                sage: L.structure_coefficients()                                        # needs sage.modules
                Finite family {}
                sage: L.structure_coefficients(True)                                    # needs sage.modules
                Finite family {(0, 1): (0, 0, 0), (0, 2): (0, 0, 0), (1, 2): (0, 0, 0)}

            ::

                sage: # needs sage.combinat sage.groups sage.modules
                sage: G = SymmetricGroup(3)
                sage: S = GroupAlgebra(G, QQ)
                sage: L = LieAlgebra(associative=S)
                sage: L.structure_coefficients()
                Finite family {((2,3), (1,2)): (1,2,3) - (1,3,2),
                               ((2,3), (1,3)): -(1,2,3) + (1,3,2),
                               ((1,2,3), (2,3)): -(1,2) + (1,3),
                               ((1,2,3), (1,2)): (2,3) - (1,3),
                               ((1,2,3), (1,3)): -(2,3) + (1,2),
                               ((1,3,2), (2,3)): (1,2) - (1,3),
                               ((1,3,2), (1,2)): -(2,3) + (1,3),
                               ((1,3,2), (1,3)): (2,3) - (1,2),
                               ((1,3), (1,2)): -(1,2,3) + (1,3,2)}
            """
            d = {}
            B = self.basis()
            K = list(B.keys())
            zero = self.zero()
            for i, x in enumerate(K):
                for y in K[i + 1:]:
                    bx = B[x]
                    by = B[y]
                    val = self.bracket(bx, by)
                    if not include_zeros and val == zero:
                        continue
                    if self._basis_key(x) > self._basis_key(y):
                        d[y,x] = -val
                    else:
                        d[x,y] = val
            return Family(d)

        def centralizer_basis(self, S):
            """
            Return a basis of the centralizer of ``S`` in ``self``.

            INPUT:

            - ``S`` -- a subalgebra of ``self`` or a list of elements that
              represent generators for a subalgebra

            .. SEEALSO::

                :meth:`centralizer`

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a, b, c = L.lie_algebra_generators()
                sage: L.centralizer_basis([a + b, 2*a + c])
                [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

                sage: # needs sage.combinat sage.modules
                sage: H = lie_algebras.Heisenberg(QQ, 2)
                sage: H.centralizer_basis(H)
                [z]

                sage: # needs sage.combinat sage.groups sage.modules
                sage: D = DescentAlgebra(QQ, 4).D()
                sage: L = LieAlgebra(associative=D)
                sage: L.centralizer_basis(L)
                [D{},
                 D{1} + D{1, 2} + D{2, 3} + D{3},
                 D{1, 2, 3} + D{1, 3} + D{2}]
                sage: D.center_basis()
                (D{},
                 D{1} + D{1, 2} + D{2, 3} + D{3},
                 D{1, 2, 3} + D{1, 3} + D{2})

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.centralizer_basis([a, c])
                [a, b, c]
                sage: L.centralizer_basis([a, e])
                [c]
            """
            from sage.matrix.constructor import matrix

            #from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
            #if isinstance(S, LieSubalgebra) or S is self:
            if S is self:
                from sage.matrix.special import identity_matrix
                m = identity_matrix(self.base_ring(), self.dimension())
            elif isinstance(S, (list, tuple)):
                m = matrix([v.to_vector() for v in self.echelon_form(S)])
            else:
                m = self.subalgebra(S).basis_matrix()

            S = self.structure_coefficients()
            sc = {}
            for k in S.keys():
                v = S[k].to_vector()
                sc[k] = v
                sc[k[1],k[0]] = -v
            X = self.basis().keys()
            d = len(X)
            c_mat = matrix(self.base_ring(),
                           [[sum(m[i,j] * sc[x,xp][k] for j,xp in enumerate(X)
                                 if (x, xp) in sc)
                             for x in X]
                            for i in range(m.nrows()) for k in range(d)])
            C = c_mat.right_kernel().basis_matrix()
            return [self.from_vector(c) for c in C]

        def centralizer(self, S):
            """
            Return the centralizer of ``S`` in ``self``.

            INPUT:

            - ``S`` -- a subalgebra of ``self`` or a list of elements that
              represent generators for a subalgebra

            .. SEEALSO::

                :meth:`centralizer_basis`

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a, b, c = L.lie_algebra_generators()
                sage: S = L.centralizer([a + b, 2*a + c]); S
                An example of a finite dimensional Lie algebra with basis:
                 the 3-dimensional abelian Lie algebra over Rational Field
                sage: S.basis_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            return self.ideal(self.centralizer_basis(S))

        @cached_method
        def center(self):
            """
            Return the center of ``self``.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: Z = L.center(); Z
                An example of a finite dimensional Lie algebra with basis: the
                 3-dimensional abelian Lie algebra over Rational Field
                sage: Z.basis_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            return self.centralizer(self)

        def normalizer_basis(self, S):
            r"""
            Return a basis of the normalizer of ``S`` in ``self``.

            INPUT:

            - ``S`` -- a subalgebra of ``self`` or a list of elements that
              represent generators for a subalgebra

            .. SEEALSO::

                :meth:`normalizer`

            EXAMPLES::

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.normalizer_basis([a, e])
                [b, c]

                sage: S = L.subalgebra([a, e])
                sage: L.normalizer_basis(S)
                [a, b, c, e]

            When the subalgebra is the ambient Lie algebra, we return the
            basis of the ambient Lie algebra::

                sage: L.normalizer_basis(L)
                Finite family {'a': a, 'b': b, 'c': c, 'd': d, 'e': e}
                sage: L.normalizer_basis([a, b, c, a, d + e, a + e])
                Finite family {'a': a, 'b': b, 'c': c, 'd': d, 'e': e}
            """
            from sage.matrix.constructor import matrix

            if S is self:
                return self.basis()
            if isinstance(S, (list, tuple)):
                m = matrix([v.to_vector() for v in self.echelon_form(S)])
            else:
                m = self.subalgebra(S).basis_matrix()

            if m.nrows() == self.dimension():
                return self.basis()

            S = self.structure_coefficients()
            sc = {}
            for k in S.keys():
                v = S[k].to_vector()
                sc[k] = v
                sc[k[1], k[0]] = -v
            X = self.basis().keys()
            d = len(X)
            ret = []
            t = m.nrows()
            c_mat = matrix(self.base_ring(),
                           [[sum(m[i,j] * sc[x,xp][k] for j, xp in enumerate(X)
                                 if (x, xp) in sc)
                             for x in X]
                            + [0]*(i*t) + [-m[j,k] for j in range(t)] + [0]*((t-i-1)*t)
                            for i in range(t) for k in range(d)])
            C = c_mat.right_kernel().basis_matrix()
            return [self.from_vector(c[:d]) for c in C]

        def normalizer(self, S):
            r"""
            Return the normalizer of ``S`` in ``self``.

            INPUT:

            - ``S`` -- a subalgebra of ``self`` or a list of elements that
              represent generators for a subalgebra

            .. SEEALSO::

                :meth:`normalizer_basis`

            EXAMPLES::

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.normalizer([a, e])
                Subalgebra generated by (b, c) of Lie algebra on
                 5 generators (a, b, c, d, e) over Rational Field
                sage: L.normalizer([a, c, e])
                Subalgebra generated by (b, c, d) of Lie algebra on
                 5 generators (a, b, c, d, e) over Rational Field
            """
            return self.subalgebra(self.normalizer_basis(S))

        @cached_method
        def derivations_basis(self):
            r"""
            Return a basis for the Lie algebra of derivations
            of ``self`` as matrices.

            A derivation `D` of an algebra is an endomorphism of `A`
            such that

            .. MATH::

                D([a, b]) = [D(a), b] + [a, D(b)]

            for all `a, b \in A`. The set of all derivations
            form a Lie algebra.

            EXAMPLES:

            We construct the derivations of the Heisenberg Lie algebra::

                sage: # needs sage.combinat sage.modules
                sage: H = lie_algebras.Heisenberg(QQ, 1)
                sage: H.derivations_basis()
                (
                [1 0 0]  [0 1 0]  [0 0 0]  [0 0 0]  [0 0 0]  [0 0 0]
                [0 0 0]  [0 0 0]  [1 0 0]  [0 1 0]  [0 0 0]  [0 0 0]
                [0 0 1], [0 0 0], [0 0 0], [0 0 1], [1 0 0], [0 1 0]
                )

            We construct the derivations of `\mathfrak{sl}_2`::

                sage: # needs sage.combinat sage.modules
                sage: sl2 = lie_algebras.sl(QQ, 2)
                sage: sl2.derivations_basis()
                (
                [ 1  0  0]  [   0    1    0]  [ 0  0  0]
                [ 0  0  0]  [   0    0 -1/2]  [ 1  0  0]
                [ 0  0 -1], [   0    0    0], [ 0 -2  0]
                )

            We verify these are derivations::

                sage: # needs sage.combinat sage.modules
                sage: D = [sl2.module_morphism(matrix=M, codomain=sl2)
                ....:      for M in sl2.derivations_basis()]
                sage: all(d(a.bracket(b)) == d(a).bracket(b) + a.bracket(d(b))
                ....:     for a in sl2.basis() for b in sl2.basis() for d in D)
                True

            REFERENCES:

            :wikipedia:`Derivation_(differential_algebra)`
            """
            from sage.matrix.constructor import matrix

            R = self.base_ring()
            B = self.basis()
            keys = list(B.keys())
            scoeffs = {(j,y,i): c for y in keys for i in keys
                       for j,c in self.bracket(B[y], B[i])
                      }
            zero = R.zero()
            data = {}
            N = len(keys)
            for ii,i in enumerate(keys):
                for ij,j in enumerate(keys[ii+1:]):
                    ijp = ij + ii + 1
                    for il,l in enumerate(keys):
                        row = ii + N * il + N**2 * ij
                        for ik,k in enumerate(keys):
                            data[row,ik+N*il] = (data.get((row,ik+N*il), zero)
                                                 + scoeffs.get((k, i, j), zero))
                            data[row,ii+N*ik] = (data.get((row,ii+N*ik), zero)
                                                 - scoeffs.get((l, k, j), zero))
                            data[row,ijp+N*ik] = (data.get((row,ijp+N*ik), zero)
                                                  - scoeffs.get((l, i, k), zero))
            mat = matrix(R, data, sparse=True)
            return tuple([matrix(R, N, N, list(b)) for b in mat.right_kernel().basis()])

        @cached_method
        def inner_derivations_basis(self):
            r"""
            Return a basis for the Lie algebra of inner derivations
            of ``self`` as matrices.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: H = lie_algebras.Heisenberg(QQ, 1)
                sage: H.inner_derivations_basis()
                (
                [0 0 0]  [0 0 0]
                [0 0 0]  [0 0 0]
                [1 0 0], [0 1 0]
                )
            """
            from sage.matrix.constructor import matrix

            R = self.base_ring()
            IDer = matrix(R, [b.adjoint_matrix().list() for b in self.basis()])
            N = self.dimension()
            return tuple([matrix(R, N, N, list(b))
                          for b in IDer.row_module().basis()])

        @cached_method
        def nilradical_basis(self):
            r"""
            Return a basis of the nilradical of ``self``.

            .. SEEALSO::

                :meth:`nilradical`

            EXAMPLES::

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.nilradical_basis()
                (a, b, c)
                sage: L.is_nilpotent()
                False

                sage: sl3 = LieAlgebra(QQ, cartan_type=['A',2])
                sage: sl3.nilradical_basis()
                ()

                sage: scoeffs = {('a','e'): {'a':1}, ('b','e'): {'a':1,'b':1},
                ....:            ('c','d'): {'a':1}, ('c','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.nilradical_basis()
                (a, b, c, d)
                sage: L.is_solvable()
                True
                sage: L.is_nilpotent()
                False

                sage: K1 = L.quotient([a])
                sage: K1.nilradical_basis()
                (b, c, d)

                sage: SL = L.subalgebra([a,b,c,d]); SL
                Subalgebra generated by (a, b, c, d) of
                 Lie algebra on 5 generators (a, b, c, d, e) over Rational Field
                sage: SL.nilradical_basis()
                (a, b, c, d)

                sage: scoeffs = {('x','z'): {'x':1, 'y':1}, ('y','z'): {'y':1}}
                sage: L.<x,y,z> = LieAlgebra(GF(3), scoeffs)
                sage: L.nilradical_basis()
                (x, y)

            We check against the generic algorithm::

                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','z'): {'x':1,'y':1}, ('y','z'): {'y':1}})
                sage: L.nilradical_basis()
                (x, y)

                sage: dim = L.dimension()
                sage: MS = MatrixSpace(L.base_ring(), dim)
                sage: gens = [b.adjoint_matrix() for b in L.basis()]
                sage: A = MS.subalgebra(gens)
                sage: RB = A.radical_basis()
                sage: mat = matrix(L.base_ring(),
                ....:              [g._vector_() for g in gens]
                ....:              + [A.lift(r)._vector_() for r in RB])
                sage: tuple([L.from_vector(w) for v in mat.right_kernel().basis()
                ....:               if (w := v[:dim])])
                (x, y)

            A positive characteristic example::

                sage: scoeffs = {('x','z'): {'x':1,'y':1}, ('y','z'): {'y':1}}
                sage: L.<x,y,z> = LieAlgebra(GF(3), scoeffs)
                sage: L.nilradical_basis()
                (x, y)
            """
            if self.base_ring().characteristic() == 0:
                L = self.solvable_radical()
                if not L.dimension():
                    return ()
                P = L.ideal(list(L.product_space(L).basis()))

                I = P.derived_subalgebra()
                if I.dimension():
                    Q = L.quotient(I)
                    ret = [Q.lift(b) for b in Q.nilradical_basis()]
                    ret.extend(b.value for b in I.basis())
                    return tuple([self(L.lift(b)) for b in ret])

                H = L.hypercenter()
                if H.dimension():
                    # the hypercenter is everything, so we don't need to compute the quotient
                    if L.dimension() == H.dimension():
                        return tuple([self(H.lift(b)) for b in H.basis()])
                    Q = L.quotient(H)
                    ret = [Q.lift(b) for b in Q.nilradical_basis()]
                    ret.extend(b.value for b in H.basis())
                    return tuple([self(L.lift(b)) for b in ret])

                from sage.matrix.constructor import matrix
                s = P.dimension()
                QP = L.quotient(P)
                MP = P.module()
                for b in QP.basis():
                    yi = QP.lift(b)
                    adj = matrix([MP.coordinate_vector(yi.bracket(P.lift(b)).to_vector())
                                  for b in P.basis()]).transpose()
                    if adj.rank() < s:
                        J = L.ideal([yi.bracket(p) for p in P.basis()])
                        QJ = L.quotient(J)
                        M = L.ideal([QJ.lift(b) for b in QJ.nilradical_basis()]
                                    + list(J.basis()))
                        return tuple([self(L.lift(b.value)) for b in M.nilradical_basis()])

                    f = adj.minimal_polynomial()
                    if not f.is_squarefree():
                        poly_ring = f.parent()
                        g = poly_ring(f / f.gcd(f.diff()))
                        phi = P.module_morphism(codomain=P, matrix=g(adj))
                        I = L.ideal([phi(p) for p in P.basis()])
                        QI = L.quotient(I)
                        M = L.ideal([QI.lift(b) for b in QI.nilradical_basis()]
                                    + list(I.basis()))
                        return tuple([self(L.lift(b.value)) for b in M.nilradical_basis()])
                return tuple([self(L.lift(b.value)) for b in P.basis()])

            # positive characteristic
            if self.is_nilpotent():
                return tuple(self.basis())

            from sage.matrix.matrix_space import MatrixSpace
            from sage.matrix.constructor import matrix
            dim = self.dimension()
            MS = MatrixSpace(self.base_ring(), dim)
            gens = [b.adjoint_matrix() for b in self.basis()]
            A = MS.subalgebra(gens)
            RB = A.radical_basis()
            mat = matrix(self.base_ring(),
                         [g._vector_() for g in gens]
                         + [A.lift(r)._vector_() for r in RB])
            return tuple([self.from_vector(w) for v in mat.left_kernel().basis()
                          if (w := v[:dim])])

        def nilradical(self):
            r"""
            Return the nilradical of ``self``.

            The *nilradical* of a Lie algebra `L` is the largest
            nilpotent ideal of `L`.

            .. SEEALSO::

                :meth:`nilradical_basis`

            EXAMPLES::

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.solvable_radical()
                Ideal (a, b, c, d, e) of
                 Lie algebra on 5 generators (a, b, c, d, e) over Rational Field
            """
            return self.ideal(self.nilradical_basis())

        @cached_method
        def solvable_radical_basis(self):
            r"""
            Return a basis of the solvable radical of ``self``.

            .. SEEALSO::

                :meth:`solvable_radical`

            EXAMPLES::

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.solvable_radical_basis()
                (a, b, c, d, e)
                sage: L.is_solvable()
                True

                sage: sl3 = LieAlgebra(QQ, cartan_type=['A',2])
                sage: sl3.solvable_radical_basis()
                ()

                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','z'): {'x':1,'y':1}, ('y','z'): {'y':1}})
                sage: S = L.subalgebra([x, y])
                sage: S.solvable_radical_basis()
                (x, y)
                sage: S.is_solvable()
                True

            Positive characteristic examples::

                sage: scoeffs = {('x','z'): {'x':1,'y':1}, ('y','z'): {'y':1}}
                sage: L.<x,y,z> = LieAlgebra(GF(3), scoeffs)
                sage: L.solvable_radical_basis()
                (x, y, z)
                sage: sl3 = LieAlgebra(GF(3), cartan_type=['A',2])
                sage: sl3.solvable_radical_basis()
                (2*h1 + h2,)
            """
            if self.base_ring().characteristic() == 0:
                P = self.derived_subalgebra()  # same ambient space as self
                if not P.dimension():
                    return tuple(self.basis())
                Bad = [b.adjoint_matrix() for b in self.basis()]
                Pad = [self(p).adjoint_matrix() for p in P.basis()]
                from sage.matrix.constructor import matrix
                mat = matrix(self.base_ring(),
                             [[(B * P).trace() for B in Bad] for P in Pad])
                return tuple([self.from_vector(c) for c in mat.right_kernel().basis_matrix()])

            # positive characteristic
            if not self.nilradical_basis():
                return ()
            dim = self.dimension()
            R = self.nilradical()
            while True:
                Q = self.quotient(R)
                RQ = Q.nilradical()
                if not RQ.dimension():  # we did not add anything
                    return tuple([self(b) for b in R.basis()])
                new_gens = [Q.lift(b.value) for b in RQ.basis()]
                R = self.ideal(list(R.basis()) + new_gens)
                if R.dimension() == dim:
                    return tuple(self.basis())

        def solvable_radical(self):
            r"""
            Return the solvable radical of ``self``.

            The *solvable radical* of a Lie algebra `L` is the largest
            solvable ideal of `L`.

            .. SEEALSO::

                :meth:`solvable_radical_basis`

            EXAMPLES::

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.solvable_radical()
                Ideal (a, b, c, d, e) of Lie algebra on 5 generators (a, b, c, d, e) over Rational Field
            """
            return self.ideal(self.solvable_radical_basis())

        def subalgebra(self, *gens, **kwds):
            r"""
            Return the subalgebra of ``self`` generated by ``gens``.

            INPUT:

            - ``gens`` -- list of generators of the subalgebra
            - ``category`` -- (optional) a subcategory of subobjects of finite
              dimensional Lie algebras with basis

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: H = lie_algebras.Heisenberg(QQ, 2)
                sage: p1,p2,q1,q2,z = H.basis()
                sage: S = H.subalgebra([p1, q1])
                sage: S.basis().list()
                [p1, q1, z]
                sage: S.basis_matrix()
                [1 0 0 0 0]
                [0 0 1 0 0]
                [0 0 0 0 1]

            Passing an extra category to a subalgebra::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, 3, step=2)
                sage: x,y,z = L.homogeneous_component_basis(1)
                sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis()
                sage: C = C.Subobjects().Graded().Stratified()
                sage: S = L.subalgebra([x, y], category=C)
                sage: S.homogeneous_component_basis(2).list()
                [X_12]
            """
            from sage.algebras.lie_algebras.subalgebra import LieSubalgebra_finite_dimensional_with_basis
            if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
                gens = gens[0]
            category = kwds.pop('category', None)
            return LieSubalgebra_finite_dimensional_with_basis(
                self, gens, category=category, **kwds)

        def ideal(self, *gens, **kwds):
            r"""
            Return the ideal of ``self`` generated by ``gens``.

            INPUT:

            - ``gens`` -- list of generators of the ideal
            - ``category`` -- (optional) a subcategory of subobjects of finite
              dimensional Lie algebras with basis

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: H = lie_algebras.Heisenberg(QQ, 2)
                sage: p1,p2,q1,q2,z = H.basis()
                sage: I = H.ideal([p1 - p2, q1 - q2])
                sage: I.basis().list()
                [-p1 + p2, -q1 + q2, z]
                sage: I.reduce(p1 + p2 + q1 + q2 + z)
                2*p1 + 2*q1

            Passing an extra category to an ideal::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y,z> = LieAlgebra(QQ, abelian=True)
                sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis()
                sage: C = C.Subobjects().Graded().Stratified()
                sage: I = L.ideal(x, y, category=C)
                sage: I.homogeneous_component_basis(1).list()
                [x, y]
            """
            from sage.algebras.lie_algebras.subalgebra import LieSubalgebra_finite_dimensional_with_basis
            if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
                gens = gens[0]
            category = kwds.pop('category', None)
            return LieSubalgebra_finite_dimensional_with_basis(
                self, gens, ideal=True, category=category, **kwds)

        @cached_method
        def is_ideal(self, A):
            """
            Return if ``self`` is an ideal of ``A``.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a, b, c = L.lie_algebra_generators()
                sage: I = L.ideal([2*a - c, b + c])
                sage: I.is_ideal(L)
                True

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})                     # needs sage.combinat sage.modules
                sage: L.is_ideal(L)                                                     # needs sage.combinat sage.modules
                True

                sage: F = LieAlgebra(QQ, 'F', representation='polynomial')              # needs sage.combinat sage.modules
                sage: L.is_ideal(F)                                                     # needs sage.combinat sage.modules
                Traceback (most recent call last):
                ...
                NotImplementedError: A must be a finite dimensional Lie algebra
                 with basis
            """
            if A == self:
                return True
            if A not in LieAlgebras(self.base_ring()).FiniteDimensional().WithBasis():
                raise NotImplementedError("A must be a finite dimensional"
                                          " Lie algebra with basis")

            from sage.matrix.constructor import matrix

            B = self.basis()
            AB = A.basis()
            try:
                b_mat = matrix(A.base_ring(), [A.bracket(b, ab).to_vector()
                                               for b in B for ab in AB])
            except (ValueError, TypeError):
                return False
            return b_mat.row_space().is_submodule(self.module())

        def quotient(self, I, names=None, category=None):
            r"""
            Return the quotient of ``self`` by the ideal ``I``.

            A quotient Lie algebra.

            INPUT:

            - ``I`` -- an ideal or a list of generators of the ideal
            - ``names`` -- (optional) string or list of strings;
              names for the basis elements of the quotient. If ``names`` is a
              string, the basis will be named ``names_1``,...,``names_n``.

            EXAMPLES:

            The Engel Lie algebra as a quotient of the free nilpotent Lie algebra
            of step 3 with 2 generators::

                sage: # needs sage.combinat sage.modules
                sage: L.<X,Y,Z,W,U> = LieAlgebra(QQ, 2, step=3)
                sage: E = L.quotient(U); E
                Lie algebra quotient L/I of dimension 4 over Rational Field where
                 L: Free Nilpotent Lie algebra on 5 generators (X, Y, Z, W, U)
                    over Rational Field
                 I: Ideal (U)
                sage: E.basis().list()
                [X, Y, Z, W]
                sage: E(X).bracket(E(Y))
                Z
                sage: Y.bracket(Z)
                -U
                sage: E(Y).bracket(E(Z))
                0
                sage: E(U)
                0

            Quotients when the base ring is not a field are not implemented::

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.Heisenberg(ZZ, 1)
                sage: L.quotient(L.an_element())
                Traceback (most recent call last):
                ...
                NotImplementedError: quotients over non-fields not implemented
            """
            from sage.algebras.lie_algebras.quotient import LieQuotient_finite_dimensional_with_basis
            return LieQuotient_finite_dimensional_with_basis(I, ambient=self,
                                                             names=names,
                                                             category=category)

        def product_space(self, L, submodule=False):
            r"""
            Return the product space ``[self, L]``.

            INPUT:

            - ``L`` -- a Lie subalgebra of ``self``
            - ``submodule`` -- boolean (default: ``False``); if ``True``, then
              the result is forced to be a submodule of ``self``

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: X = L.subalgebra([a, b + c])
                sage: L.product_space(X)
                An example of a finite dimensional Lie algebra with basis:
                 the 0-dimensional abelian Lie algebra over Rational Field
                  with basis matrix: []
                sage: Y = L.subalgebra([a, 2*b - c])
                sage: X.product_space(Y)
                An example of a finite dimensional Lie algebra with basis:
                 the 0-dimensional abelian Lie algebra over Rational Field
                  with basis matrix: []

            ::

                sage: # needs sage.combinat sage.modules
                sage: H = lie_algebras.Heisenberg(ZZ, 4)
                sage: Hp = H.product_space(H, submodule=True).basis()
                sage: [H.from_vector(v) for v in Hp]
                [z]

            ::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: Lp = L.product_space(L)   # not implemented
                sage: Lp                        # not implemented
                Subalgebra generated of
                 Lie algebra on 2 generators (x, y) over Rational Field
                 with basis: (x,)
                sage: Lp.product_space(L)       # not implemented
                Subalgebra generated of
                 Lie algebra on 2 generators (x, y) over Rational Field
                 with basis: (x,)
                sage: L.product_space(Lp)       # not implemented
                Subalgebra generated of
                 Lie algebra on 2 generators (x, y) over Rational Field
                 with basis: (x,)
                sage: Lp.product_space(Lp)      # not implemented
                Subalgebra generated of
                 Lie algebra on 2 generators (x, y) over Rational Field
                 with basis: ()
            """
            from sage.matrix.constructor import matrix

            # Make sure we lift everything to the ambient space
            if self in LieAlgebras(self.base_ring()).Subobjects():
                A = self.ambient()
            elif L in LieAlgebras(L.base_ring()).Subobjects():
                A = L.ambient()
            else:
                A = self

            if L not in self.category():
                # L might be a submodule of A.module()
                LB = [self.from_vector(b) for b in L.basis()]
            else:
                LB = L.basis()

            B = self.basis()
            b_mat = matrix(A.base_ring(), [A.bracket(b, lb).to_vector()
                                           for b in B for lb in LB])
            if submodule is True or not (self.is_ideal(A) and L.is_ideal(A)):
                return b_mat.row_space()
            # We echelonize the matrix here
            # TODO: Do we want to?
            b_mat.echelonize()
            r = b_mat.rank()
            gens = [A.from_vector(row) for row in b_mat.rows()[:r]]
            return A.ideal(gens)

        @cached_method
        def derived_subalgebra(self):
            """
            Return the derived subalgebra of ``self``.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.derived_subalgebra()
                An example of a finite dimensional Lie algebra with basis:
                 the 0-dimensional abelian Lie algebra over Rational Field
                 with basis matrix:
                []

            If ``self`` is semisimple, then the derived subalgebra is ``self``::

                sage: # needs sage.combinat sage.modules
                sage: sl3 = LieAlgebra(QQ, cartan_type=['A', 2])
                sage: sl3.derived_subalgebra()
                Lie algebra of ['A', 2] in the Chevalley basis
                sage: sl3 is sl3.derived_subalgebra()
                True
            """
            if self.is_semisimple():
                return self
            else:
                return self.product_space(self)

        @cached_method
        def derived_series(self):
            r"""
            Return the derived series `(\mathfrak{g}^{(i)})_i` of ``self``
            where the rightmost
            `\mathfrak{g}^{(k)} = \mathfrak{g}^{(k+1)} = \cdots`.

            We define the derived series of a Lie algebra `\mathfrak{g}`
            recursively by `\mathfrak{g}^{(0)} := \mathfrak{g}` and

            .. MATH::

                \mathfrak{g}^{(k+1)} =
                [\mathfrak{g}^{(k)}, \mathfrak{g}^{(k)}]

            and recall that
            `\mathfrak{g}^{(k)} \supseteq \mathfrak{g}^{(k+1)}`.
            Alternatively we can express this as

            .. MATH::

                \mathfrak{g} \supseteq [\mathfrak{g}, \mathfrak{g}] \supseteq
                \bigl[ [\mathfrak{g}, \mathfrak{g}], [\mathfrak{g},
                \mathfrak{g}] \bigr] \supseteq
                \biggl[ \bigl[ [\mathfrak{g}, \mathfrak{g}], [\mathfrak{g},
                \mathfrak{g}] \bigr], \bigl[ [\mathfrak{g}, \mathfrak{g}],
                [\mathfrak{g}, \mathfrak{g}] \bigr] \biggr] \supseteq \cdots.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.derived_series()
                (An example of a finite dimensional Lie algebra with basis:
                    the 3-dimensional abelian Lie algebra over Rational Field,
                 An example of a finite dimensional Lie algebra with basis:
                    the 0-dimensional abelian Lie algebra over Rational Field
                    with basis matrix: [])

            ::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.derived_series()
                (Lie algebra on 2 generators (x, y) over Rational Field,
                 Ideal (x) of Lie algebra on 2 generators (x, y) over Rational Field,
                 Ideal () of Lie algebra on 2 generators (x, y) over Rational Field)

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.derived_series()
                (Lie algebra on 5 generators (a, b, c, d, e) over Rational Field,
                 Ideal (a, b, c) of Lie algebra on 5 generators (a, b, c, d, e) over Rational Field,
                 Ideal () of Lie algebra on 5 generators (a, b, c, d, e) over Rational Field)
            """
            L = [self]
            while L[-1].dimension() > 0:
                p = L[-1].derived_subalgebra()
                if L[-1].dimension() == p.dimension():
                    break
                L.append(p)
            return tuple(L)

        @cached_method
        def lower_central_series(self, submodule=False):
            r"""
            Return the lower central series `(\mathfrak{g}_{i})_i`
            of ``self`` where the rightmost
            `\mathfrak{g}_k = \mathfrak{g}_{k+1} = \cdots`.

            INPUT:

            - ``submodule`` -- boolean (default: ``False``); if ``True``, then
              the result is given as submodules of ``self``

            We define the lower central series of a Lie algebra `\mathfrak{g}`
            recursively by `\mathfrak{g}_0 := \mathfrak{g}` and

            .. MATH::

                \mathfrak{g}_{k+1} = [\mathfrak{g}, \mathfrak{g}_{k}]

            and recall that `\mathfrak{g}_{k} \supseteq \mathfrak{g}_{k+1}`.
            Alternatively we can express this as

            .. MATH::

                \mathfrak{g} \supseteq [\mathfrak{g}, \mathfrak{g}] \supseteq
                \bigl[ [\mathfrak{g}, \mathfrak{g}], \mathfrak{g} \bigr]
                \supseteq \Bigl[\bigl[ [\mathfrak{g}, \mathfrak{g}],
                \mathfrak{g} \bigr], \mathfrak{g}\Bigr] \supseteq \cdots.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.derived_series()
                (An example of a finite dimensional Lie algebra with basis:
                  the 3-dimensional abelian Lie algebra over Rational Field,
                 An example of a finite dimensional Lie algebra with basis:
                  the 0-dimensional abelian Lie algebra over Rational Field
                   with basis matrix: [])

            The lower central series as submodules::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.lower_central_series(submodule=True)
                (Sparse vector space of dimension 2 over Rational Field,
                 Vector space of degree 2 and dimension 1 over Rational Field
                  Basis matrix: [1 0])

            ::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.lower_central_series()
                (Lie algebra on 2 generators (x, y) over Rational Field,
                 Ideal (x) of Lie algebra on 2 generators (x, y) over Rational Field)

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.lower_central_series()
                (Lie algebra on 5 generators (a, b, c, d, e) over Rational Field,
                 Ideal (a, b, c) of Lie algebra on 5 generators (a, b, c, d, e) over Rational Field,
                 Ideal (a, b) of Lie algebra on 5 generators (a, b, c, d, e) over Rational Field)
            """
            if submodule:
                L = [self.module()]
            else:
                L = [self]
            while L[-1].dimension() > 0:
                s = self.product_space(L[-1], submodule=submodule)
                if L[-1].dimension() == s.dimension():
                    break
                L.append(s)
            return tuple(L)

        @cached_method
        def upper_central_series(self):
            r"""
            Return the upper central series `(Z_i(\mathfrak{g}))_i`
            of ``self`` where the rightmost
            `Z_k(\mathfrak{g}) = Z_{k+1}(\mathfrak{g}) = \cdots`.

            The *upper central series* of a Lie algebra `\mathfrak{g}` is
            defined recursively by `Z_0(\mathfrak{g}) := Z(\mathfrak{g})` and

            .. MATH::

                Z_{k+1}(\mathfrak{g}) / Z_k(\mathfrak{g})
                = Z(\mathfrak{g} / Z_k(\mathfrak{g}),

            and recall that `Z(\mathfrak{g})` is the :meth:`center`
            of `\mathfrak{g}`.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.upper_central_series()
                [An example of a finite dimensional Lie algebra with basis:
                 the 3-dimensional abelian Lie algebra over Rational Field]

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.upper_central_series()
                [Ideal () of Lie algebra on 2 generators (x, y) over Rational Field]

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.upper_central_series()
                [Ideal (c) of Lie algebra on 5 generators (a, b, c, d, e) over Rational Field]

                sage: L = lie_algebras.Heisenberg(QQ, 3)
                sage: L.upper_central_series()
                [Ideal (z) of Heisenberg algebra of rank 3 over Rational Field,
                 Heisenberg algebra of rank 3 over Rational Field]
            """
            I = self.center()
            if I.dimension() == 0:
                return [I]
            ret = [I]
            dim = self.dimension()
            while True:
                Q = self.quotient(I)
                Z = Q.center()
                if not Z.dimension():  # we did not add anything
                    return ret
                new_gens = [Q.lift(b.value) for b in Z.basis()]
                I = self.ideal(list(I.basis()) + new_gens)
                if I.dimension() == dim:
                    ret.append(self)
                    return ret
                ret.append(I)

        def hypercenter(self):
            r"""
            Return the hypercenter of ``self``.

            EXAMPLES::

                sage: SGA3 = SymmetricGroup(3).algebra(QQ)
                sage: L = LieAlgebra(associative=SGA3)
                sage: L.hypercenter()
                Ideal ((), (1,2,3) + (1,3,2), (2,3) + (1,2) + (1,3)) of
                 Lie algebra of Symmetric group algebra of order 3
                 over Rational Field

                sage: L = lie_algebras.Heisenberg(QQ, 3)
                sage: L.hypercenter()
                Heisenberg algebra of rank 3 over Rational Field
            """
            return self.upper_central_series()[-1]

        def is_abelian(self):
            """
            Return if ``self`` is an abelian Lie algebra.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_abelian()
                True

            ::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.is_abelian()
                False
            """
            return len(self.structure_coefficients()) == 0
            # TODO: boolean handling of empty family
            #return not self.structure_coefficients()

        def is_solvable(self):
            r"""
            Return if ``self`` is a solvable Lie algebra.

            A Lie algebra is solvable if the derived series eventually
            becomes `0`.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_solvable()
                True

            ::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.is_solvable()           # not implemented
                False
            """
            return not self.derived_series()[-1].dimension()

        def is_nilpotent(self):
            r"""
            Return if ``self`` is a nilpotent Lie algebra.

            A Lie algebra is nilpotent if the lower central series eventually
            becomes `0`.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_nilpotent()
                True
            """
            return not self.lower_central_series()[-1].dimension()

        def is_semisimple(self):
            """
            Return if ``self`` if a semisimple Lie algebra.

            A Lie algebra is semisimple if the solvable radical is zero. In
            characteristic 0, this is equivalent to saying the Killing form
            is non-degenerate.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_semisimple()
                False

            Positive characteristic examples::

                sage: L.<x,y,z> = LieAlgebra(GF(3), {('x','z'): {'x':1, 'y':1}, ('y','z'): {'y':1}})
                sage: L.is_semisimple()
                False

                sage: sp4 = LieAlgebra(GF(3), cartan_type=['C',2])
                sage: sp4.killing_form_matrix().det()
                0
                sage: sp4.solvable_radical_basis()  # long time
                ()
                sage: sp4.is_semisimple()  # long time
                True
            """
            if self.base_ring().characteristic() == 0:
                return not self.killing_form_matrix().is_singular()
            if not self.killing_form_matrix().is_singular():
                return True
            return not self.solvable_radical_basis()

        @cached_method(key=_ce_complex_key)
        def chevalley_eilenberg_complex(self, M=None, dual=False, sparse=True, ncpus=None):
            r"""
            Return the Chevalley-Eilenberg complex of ``self``.

            Let `\mathfrak{g}` be a Lie algebra and `M` be a right
            `\mathfrak{g}`-module. The *Chevalley-Eilenberg complex*
            is the chain complex on

            .. MATH::

                C_{\bullet}(\mathfrak{g}, M) =
                M \otimes \bigwedge\nolimits^{\bullet} \mathfrak{g},

            where the differential is given by

            .. MATH::

                d(m \otimes g_1 \wedge \cdots \wedge g_p) =
                \sum_{i=1}^p (-1)^{i+1}
                  (m g_i) \otimes g_1 \wedge \cdots \wedge
                  \hat{g}_i \wedge \cdots \wedge g_p +
                \sum_{1 \leq i < j \leq p} (-1)^{i+j}
                  m \otimes [g_i, g_j] \wedge
                  g_1 \wedge \cdots \wedge \hat{g}_i
                  \wedge \cdots \wedge \hat{g}_j
                  \wedge \cdots \wedge g_p.

            INPUT:

            - ``M`` -- (default: the trivial 1-dimensional module)
              one of the following:

              * a module `M` with an action of ``self``
              * a dictionary whose keys are basis elements and values
                are matrices representing a Lie algebra homomorphism
                defining the representation

            - ``dual`` -- boolean (default: ``False``); if ``True``, causes
              the dual of the complex to be computed
            - ``sparse`` -- boolean (default: ``True``); whether to use sparse
              or dense matrices
            - ``ncpus`` -- (optional) how many cpus to use

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.sl(ZZ, 2)
                sage: C = L.chevalley_eilenberg_complex(); C
                Chain complex with at most 4 nonzero terms over Integer Ring
                sage: ascii_art(C)
                                          [-2  0  0]       [0]
                                          [ 0  1  0]       [0]
                            [0 0 0]       [ 0  0 -2]       [0]
                 0 <-- C_0 <-------- C_1 <----------- C_2 <---- C_3 <-- 0

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['C',2])
                sage: C = L.chevalley_eilenberg_complex()  # long time
                sage: [C.free_module_rank(i) for i in range(11)]  # long time
                [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1]

                sage: # needs sage.combinat sage.modules
                sage: g = lie_algebras.sl(QQ, 2)
                sage: E, F, H = g.basis()
                sage: n = g.subalgebra([F, H])
                sage: ascii_art(n.chevalley_eilenberg_complex())
                                        [ 0]
                            [0 0]       [-2]
                 0 <-- C_0 <------ C_1 <----- C_2 <-- 0

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
                sage: f = ({x: Matrix([[1,0],[0,0]]), y: Matrix([[0,1],[0,0]])})
                sage: C = L.chevalley_eilenberg_complex(f); C
                Chain complex with at most 3 nonzero terms over Rational Field
                sage: ascii_art(C)
                                            [ 0 -1]
                                            [ 2  0]
                            [1 0 0 1]       [ 0  0]
                            [0 0 0 0]       [ 0  1]
                 0 <-- C_0 <---------- C_1 <-------- C_2 <-- 0

                sage: ascii_art(L.chevalley_eilenberg_complex(f, sparse=False))
                                            [ 0 -1]
                                            [ 2  0]
                            [1 0 0 1]       [ 0  0]
                            [0 0 0 0]       [ 0  1]
                 0 <-- C_0 <---------- C_1 <-------- C_2 <-- 0

            REFERENCES:

            - :wikipedia:`Lie_algebra_cohomology#Chevalley-Eilenberg_complex`
            - [Wei1994]_ Chapter 7
            """
            if dual:
                return self.chevalley_eilenberg_complex(M, dual=False,
                                                        sparse=sparse,
                                                        ncpus=ncpus).dual()

            import itertools
            from itertools import combinations, product
            from sage.arith.misc import binomial
            from sage.matrix.matrix_space import MatrixSpace
            from sage.algebras.lie_algebras.representation import Representation_abstract
            R = self.base_ring()
            zero = R.zero()
            mone = -R.one()

            # Make sure we specify the ordering of the basis
            LB = self.basis()
            LK = list(LB.keys())
            LB = [LB[k] for k in LK]
            LI = list(range(len(LK)))
            Lmod = self.module()
            ambient = Lmod.is_ambient()

            if M is not None:
                if not isinstance(M, Representation_abstract):
                    M = self.representation(M)

                MB = M.basis()
                MK = list(MB.keys())
                MB = [MB[k] for k in MK]
                MI = list(range(len(MK)))

            def sgn(k, X):
                """
                Insert a new entry ``k`` into a strictly increasing
                list ``X`` in such a way that the resulting list is
                still strictly increasing.
                The return value is the pair ``(s, Y)``, where ``Y``
                is the resulting list (as tuple) and ``s`` is the
                Koszul sign incurred by the insertion (with the
                understanding that ``k`` originally stood to the
                left of the list).
                If ``k`` is already in ``X``, then the return value
                is ``(zero, None)``.
                """
                Y = list(X)
                for i in range(len(X)-1, -1, -1):
                    val = X[i]
                    if val == k:
                        return zero, None
                    if k > val:
                        Y.insert(i+1, k)
                        return mone**(i+1), tuple(Y)
                Y.insert(0, k)
                return R.one(), tuple(Y)

            from sage.parallel.decorate import parallel
            from sage.matrix.constructor import matrix

            @parallel(ncpus=ncpus)
            def compute_diff(k):
                """
                Build the ``k``-th differential (in parallel).
                """
                # The indices for the exterior algebra
                ext_ind = {tuple(X): i for i, X in enumerate(combinations(LI, k-1))}

                # Compute the part independent of the module first ("part 2" of the computation)
                if sparse:
                    p2_data = {}
                    row = 0
                else:
                    p2_data = []
                if not sparse:
                    zv = [zero] * len(ext_ind)
                for X in combinations(LI, k):
                    if not sparse:
                        ret = list(zv)
                    for i in range(k):
                        Y = list(X)
                        Y.pop(i)
                        # This is where we would do the action on
                        #   the coefficients module
                        #ret[indices[tuple(Y)]] += mone**i * zero
                        for j in range(i+1,k):
                            # We shift j by 1 because we already removed
                            #   an earlier element from X.
                            Z = tuple(Y[:j-1] + Y[j:])
                            elt = mone**(i+j+1) * LB[X[i]].bracket(LB[X[j]])
                            if not elt:
                                continue
                            if ambient:
                                vec = elt.to_vector()
                            else:
                                vec = Lmod.coordinate_vector(elt.to_vector())
                            for key, coeff in vec.iteritems():
                                if not coeff:
                                    continue
                                s, A = sgn(key, Z)
                                if A is None:
                                    continue
                                if sparse:
                                    coords = (row, ext_ind[A])
                                    if coords in p2_data:
                                        p2_data[coords] += s * coeff
                                    else:
                                        p2_data[coords] = s * coeff
                                else:
                                    ret[ext_ind[A]] += s * coeff
                    if sparse:
                        row += 1
                    else:
                        p2_data.append(ret)

                nrows = binomial(len(LI), k)
                ncols = binomial(len(LI), k-1)
                MS = MatrixSpace(R, nrows, ncols, sparse=sparse)
                if M is None:
                    p2 = MS(p2_data).transpose()
                    p2.set_immutable()
                    return p2
                p2 = matrix.identity(len(MI)).tensor_product(MS(p2_data)).transpose()

                ten_ind = {tuple(Y): i for i, Y in enumerate(product(MI, ext_ind))}

                # Now compute the part from the module ("part 1")
                if sparse:
                    p1_data = {}
                    row = 0
                else:
                    p1_data = []
                if not sparse:
                    zv = [zero] * len(ten_ind)
                for v, X in product(MI, combinations(LI, k)):
                    if not sparse:
                        ret = list(zv)
                    for i in range(k):
                        # We do mone**i because we are 0-based
                        elt = mone**i * LB[X[i]] * MB[v]
                        if not elt:
                            continue
                        Y = X[:i] + X[i+1:]
                        for j in MI:
                            coeff = elt[MK[j]]
                            if not coeff:
                                continue
                            if sparse:
                                coords = (row, ten_ind[j, Y])
                                if coords in p1_data:
                                    p1_data[coords] += coeff
                                else:
                                    p1_data[coords] = coeff
                            else:
                                ret[ten_ind[j, Y]] += coeff
                    if sparse:
                        row += 1
                    else:
                        p1_data.append(ret)

                nrows = len(MI) * binomial(len(LI), k)
                ncols = len(ten_ind)
                MS = MatrixSpace(R, nrows, ncols, sparse=sparse)
                ret = MS(p1_data).transpose() + p2
                ret.set_immutable()
                return ret

            from sage.homology.chain_complex import ChainComplex
            ind = list(range(1, len(LI) + 1))
            chain_data = {X[0][0]: M for X, M in compute_diff(ind)}
            C = ChainComplex(chain_data, degree_of_differential=-1)
            return C

        def homology(self, deg=None, M=None, sparse=True, ncpus=None):
            r"""
            Return the Lie algebra homology of ``self``.

            The Lie algebra homology is the homology of the
            Chevalley-Eilenberg chain complex.

            INPUT:

            - ``deg`` -- the degree of the homology (optional)
            - ``M`` -- (default: the trivial module) a right module
              of ``self``
            - ``sparse`` -- boolean (default: ``True``); whether to use sparse
              matrices for the Chevalley-Eilenberg chain complex
            - ``ncpus`` -- (optional) how many cpus to use when
              computing the Chevalley-Eilenberg chain complex

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.cross_product(QQ)
                sage: L.homology()
                {0: Vector space of dimension 1 over Rational Field,
                 1: Vector space of dimension 0 over Rational Field,
                 2: Vector space of dimension 0 over Rational Field,
                 3: Vector space of dimension 1 over Rational Field}

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.pwitt(GF(5), 5)
                sage: L.homology()
                {0: Vector space of dimension 1 over Finite Field of size 5,
                 1: Vector space of dimension 0 over Finite Field of size 5,
                 2: Vector space of dimension 1 over Finite Field of size 5,
                 3: Vector space of dimension 1 over Finite Field of size 5,
                 4: Vector space of dimension 0 over Finite Field of size 5,
                 5: Vector space of dimension 1 over Finite Field of size 5}

                sage: # needs sage.combinat sage.modules
                sage: d = {('x', 'y'): {'y': 2}}
                sage: L.<x,y> = LieAlgebra(ZZ, d)
                sage: L.homology()
                {0: Z, 1: Z x C2, 2: 0}

            .. SEEALSO::

                :meth:`chevalley_eilenberg_complex`
            """
            C = self.chevalley_eilenberg_complex(M=M, sparse=sparse,
                                                 ncpus=ncpus)
            return C.homology(deg=deg)

        def cohomology(self, deg=None, M=None, sparse=True, ncpus=None):
            r"""
            Return the Lie algebra cohomology of ``self``.

            The Lie algebra cohomology is the cohomology of the
            Chevalley-Eilenberg cochain complex (which is the dual
            of the Chevalley-Eilenberg chain complex).

            Let `\mathfrak{g}` be a Lie algebra and `M` a left
            `\mathfrak{g}`-module. It is known that `H^0(\mathfrak{g}; M)`
            is the subspace of `\mathfrak{g}`-invariants of `M`:

            .. MATH::

                H^0(\mathfrak{g}; M) = M^{\mathfrak{g}}
                = \{ m \in M \mid g m = 0
                    \text{ for all } g \in \mathfrak{g} \}.

            Additionally, `H^1(\mathfrak{g}; M)` is the space of
            derivations `\mathfrak{g} \to M`
            modulo the space of inner derivations, and
            `H^2(\mathfrak{g}; M)` is the space of equivalence classes
            of Lie algebra extensions of `\mathfrak{g}` by `M`.

            INPUT:

            - ``deg`` -- the degree of the homology (optional)
            - ``M`` -- (default: the trivial module) a right module
              of ``self``
            - ``sparse`` -- boolean (default: ``True``); whether to use sparse
              matrices for the Chevalley-Eilenberg chain complex
            - ``ncpus`` -- (optional) how many cpus to use when
              computing the Chevalley-Eilenberg chain complex

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.so(QQ, 4)
                sage: L.cohomology()
                {0: Vector space of dimension 1 over Rational Field,
                 1: Vector space of dimension 0 over Rational Field,
                 2: Vector space of dimension 0 over Rational Field,
                 3: Vector space of dimension 2 over Rational Field,
                 4: Vector space of dimension 0 over Rational Field,
                 5: Vector space of dimension 0 over Rational Field,
                 6: Vector space of dimension 1 over Rational Field}

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.Heisenberg(QQ, 2)
                sage: L.cohomology()
                {0: Vector space of dimension 1 over Rational Field,
                 1: Vector space of dimension 4 over Rational Field,
                 2: Vector space of dimension 5 over Rational Field,
                 3: Vector space of dimension 5 over Rational Field,
                 4: Vector space of dimension 4 over Rational Field,
                 5: Vector space of dimension 1 over Rational Field}

                sage: # needs sage.combinat sage.modules
                sage: d = {('x', 'y'): {'y': 2}}
                sage: L.<x,y> = LieAlgebra(ZZ, d)
                sage: L.cohomology()
                {0: Z, 1: Z, 2: C2}

            .. SEEALSO::

                :meth:`chevalley_eilenberg_complex`

            REFERENCES:

            - :wikipedia:`Lie_algebra_cohomology`
            """
            C = self.chevalley_eilenberg_complex(M=M, dual=True, sparse=sparse,
                                                 ncpus=ncpus)
            return C.homology(deg=deg)

        def as_finite_dimensional_algebra(self):
            """
            Return ``self`` as a :class:`FiniteDimensionalAlgebra`.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.cross_product(QQ)
                sage: x, y, z = L.basis()
                sage: F = L.as_finite_dimensional_algebra()
                sage: X, Y, Z = F.basis()
                sage: x.bracket(y)
                Z
                sage: X * Y
                Z
            """
            from sage.matrix.constructor import matrix

            K = self._basis_ordering
            mats = []
            R = self.base_ring()
            S = dict(self.structure_coefficients())
            V = self._dense_free_module()
            zero_vec = V.zero()
            for k in K:
                M = []
                for kp in K:
                    if (k, kp) in S:
                        M.append( -S[k,kp].to_vector() )
                    elif (kp, k) in S:
                        M.append( S[kp,k].to_vector() )
                    else:
                        M.append( zero_vec )
                mats.append(matrix(R, M))
            from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra import FiniteDimensionalAlgebra
            return FiniteDimensionalAlgebra(R, mats, names=self._names)

        def morphism(self, on_generators, codomain=None, base_map=None, check=True):
            r"""
            Return a Lie algebra morphism defined by images of a Lie
            generating subset of ``self``.

            INPUT:

            - ``on_generators`` -- dictionary ``{X: Y}`` of the images `Y`
              in ``codomain`` of elements `X` of ``domain``
            - ``codomain`` -- a Lie algebra (optional); this is inferred
              from the values of ``on_generators`` if not given
            - ``base_map`` -- a homomorphism from the base ring to something
              coercing into the codomain
            - ``check`` -- boolean (default: ``True``); if ``False`` the
              values  on the Lie brackets implied by ``on_generators`` will
              not be checked for contradictory values

            .. NOTE::

                The keys of ``on_generators`` need to generate ``domain``
                as a Lie algebra.

            .. SEEALSO::

                :class:`sage.algebras.lie_algebras.morphism.LieAlgebraMorphism_from_generators`

            EXAMPLES:

            A quotient type Lie algebra morphism ::

                sage: # needs sage.combinat sage.modules
                sage: L.<X,Y,Z,W> = LieAlgebra(QQ, {('X','Y'): {'Z': 1},
                ....:                               ('X','Z'): {'W': 1}})
                sage: K.<A,B> = LieAlgebra(QQ, abelian=True)
                sage: L.morphism({X: A, Y: B})
                Lie algebra morphism:
                  From: Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
                  To:   Abelian Lie algebra on 2 generators (A, B) over Rational Field
                  Defn: X |--> A
                        Y |--> B
                        Z |--> 0
                        W |--> 0

            The reverse map `A \mapsto X`, `B \mapsto Y` does not define a Lie
            algebra morphism, since `[A,B] = 0`, but `[X,Y] \neq 0`::

                sage: # needs sage.combinat sage.modules
                sage: K.morphism({A:X, B: Y})
                Traceback (most recent call last):
                ...
                ValueError: this does not define a Lie algebra morphism;
                 contradictory values for brackets of length 2

            However, it is still possible to create a morphism that acts nontrivially
            on the coefficients, even though it's not a Lie algebra morphism
            (since it isn't linear)::

                sage: # needs sage.combinat sage.modules sage.rings.number_fields
                sage: R.<x> = ZZ[]
                sage: K.<i> = NumberField(x^2 + 1)
                sage: cc = K.hom([-i])
                sage: L.<X,Y,Z,W> = LieAlgebra(K, {('X','Y'): {'Z': 1},
                ....:                              ('X','Z'): {'W': 1}})
                sage: M.<A,B> = LieAlgebra(K, abelian=True)
                sage: phi = L.morphism({X: A, Y: B}, base_map=cc)
                sage: phi(X)
                A
                sage: phi(i*X)
                -i*A
            """
            from sage.algebras.lie_algebras.morphism import LieAlgebraMorphism_from_generators
            return LieAlgebraMorphism_from_generators(on_generators, domain=self,
                                                      codomain=codomain, base_map=base_map, check=check)

        @cached_method
        def universal_polynomials(self):
            r"""
            Return the family of universal polynomials of ``self``.

            The *universal polynomials* of a Lie algebra `L` with
            basis `\{e_i\}_{i \in I}` and structure coefficients
            `[e_i, e_j] = \tau_{ij}^a e_a` is given by

            .. MATH::

                P_{aij} = \sum_{u \in I} \tau_{ij}^u X_{au}
                - \sum_{s,t \in I} \tau_{st}^a X_{si} X_{tj},

            where `a,i,j \in I`.

            REFERENCES:

            - [AM2020]_

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.universal_polynomials()
                Finite family {('x', 'x', 'y'): X01*X10 - X00*X11 + X00,
                               ('y', 'x', 'y'): X10}

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['A',1])
                sage: list(L.universal_polynomials())
                [-2*X01*X10 + 2*X00*X11 - 2*X00,
                 -2*X02*X10 + 2*X00*X12 + X01,
                 -2*X02*X11 + 2*X01*X12 - 2*X02,
                 X01*X20 - X00*X21 - 2*X10,
                 X02*X20 - X00*X22 + X11,
                 X02*X21 - X01*X22 - 2*X12,
                 -2*X11*X20 + 2*X10*X21 - 2*X20,
                 -2*X12*X20 + 2*X10*X22 + X21,
                 -2*X12*X21 + 2*X11*X22 - 2*X22]

                sage: # long time, needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['B', 2])
                sage: al = RootSystem(['B', 2]).root_lattice().simple_roots()
                sage: k = list(L.basis().keys())[0]
                sage: UP = L.universal_polynomials()
                sage: len(UP)
                450
                sage: UP[al[2], al[1], -al[1]]
                X0_7*X4_1 - X0_1*X4_7 - 2*X0_7*X5_1 + 2*X0_1*X5_7 + X2_7*X7_1
                 - X2_1*X7_7 - X3_7*X8_1 + X3_1*X8_7 + X0_4
            """
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            I = self.basis().keys()
            n = len(I)
            s_coeffs = self.structure_coefficients(True)
            zero = self.base_ring().zero()

            def sc(i, j):
                if i == j:
                    return zero
                if i > j:
                    return -s_coeffs[I[j], I[i]]
                return s_coeffs[I[i], I[j]]
            d = {}
            keys = []
            if n >= 10:
                vs = 'X{}_{}'
            else:
                vs = 'X{}{}'
            R = PolynomialRing(self.base_ring(), ','.join(vs.format(i,j)
                                                          for i in range(n)
                                                          for j in range(n)))
            X = [[R.gen(i+n*j) for i in range(n)] for j in range(n)]
            for a in range(n):
                for i in range(n):
                    for j in range(i+1, n):
                        k = (I[a], I[i], I[j])
                        keys.append(k)
                        if i != j:
                            s = sc(i, j)
                            d[k] = (R.sum(s[I[u]] * X[a][u] for u in range(n))
                                    - R.sum(sc(s,t)[I[a]] * X[s][i] * X[t][j]
                                            for s in range(n) for t in range(n) if s != t))
                        else:
                            d[k] = -R.sum(sc(s,t)[I[a]] * X[s][i] * X[t][j]
                                          for s in range(n) for t in range(n) if s != t)
            return Family(keys, d.__getitem__)

        @cached_method
        def universal_commutative_algebra(self):
            r"""
            Return the universal commutative algebra associated to ``self``.

            Let `I` be the index set of the basis of ``self``. Let
            `\mathcal{P} = \{P_{a,i,j}\}_{a,i,j \in I}` denote the
            universal polynomials of a Lie algebra `L`. The *universal
            commutative algebra* associated to `L` is the quotient
            ring `R[X_{ij}]_{i,j \in I} / (\mathcal{P})`.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: A = L.universal_commutative_algebra()
                sage: a, b, c, d = A.gens()
                sage: a, b, c, d
                (X00bar, X01bar, 0, X11bar)
                sage: a*d - a
                0
            """
            P = list(self.universal_polynomials())
            R = P[0].parent()
            return R.quotient(P)

        def casimir_element(self, order=2, UEA=None, force_generic=False, basis=False):
            r"""
            Return a Casimir element of order ``order`` in the universal
            enveloping algebra of ``self``.

            A *Casimir element* of order `k` is a distinguished basis element
            for the center of `U(\mathfrak{g})` of homogeneous degree `k`
            (that is, it is an element of `U_k \setminus U_{k-1}`, where
            `\{U_i\}_{i=0}^{\infty}` is the natural filtration of
            `U(\mathfrak{g})`). When `\mathfrak{g}` is a simple Lie algebra,
            then this spans `Z(U(\mathfrak{g}))_k`.

            INPUT:

            - ``order`` -- (default: ``2``) the order of the Casimir element
            - ``UEA`` -- (optional) the universal enveloping algebra
              implementation to return the result in
            - ``force_generic`` -- boolean (default: ``False``); if ``True``
              for the quadratic order, then this uses the default algorithm
              (otherwise this is ignored)
            - ``basis`` -- boolean (default: ``False``); if ``True``, this
              returns a basis of all Casimir elements of order ``order`` as a
              list

            ALGORITHM:

            For the quadratic order (i.e., ``order=2``), then this uses
            `K^{ij}`, the inverse of the Killing form matrix, to compute
            `C_{(2)} = \sum_{i,j} K^{ij} X_i \cdots X_j`, where `\{X_1, \ldots,
            X_n\}` is a basis for `\mathfrak{g}`. Otherwise this solves the
            system of equations

            .. MATH::

                f_{aj}^b \kappa^{jc\cdots d} + f_{aj}^c \kappa^{cj\cdots d}
                \cdots + f_{aj}^d \kappa^{bc \cdots j}

            for the symmetric tensor `\kappa^{i_1 \cdots i_k}`, where `k` is
            the ``order``. This system comes from `[X_i, C_{(k)}] = 0` with

            .. MATH::

                C_{(k)} = \sum_{i_1, \ldots, i_k}^n
                \kappa^{i_1 \cdots i_k} X_{i_1} \cdots X_{i_k}.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['A', 1])
                sage: C = L.casimir_element(); C
                1/8*b1^2 + 1/2*b0*b2 - 1/4*b1
                sage: U = L.universal_enveloping_algebra()
                sage: all(g * C == C * g for g in U.gens())
                True
                sage: U = L.pbw_basis()
                sage: C = L.casimir_element(UEA=U); C
                1/2*PBW[alpha[1]]*PBW[-alpha[1]] + 1/8*PBW[alphacheck[1]]^2
                 - 1/4*PBW[alphacheck[1]]
                sage: all(g * C == C * g for g in U.algebra_generators())
                True

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['B', 2])
                sage: U = L.pbw_basis()
                sage: C = L.casimir_element(UEA=U)
                sage: all(g * C == C * g for g in U.algebra_generators())
                True

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['C', 3])
                sage: U = L.pbw_basis()
                sage: C = L.casimir_element(UEA=U)
                sage: all(g * C == C * g for g in U.algebra_generators())
                True

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['A', 1])
                sage: C4 = L.casimir_element(order=4, UEA=L.pbw_basis()); C4
                4*PBW[alpha[1]]^2*PBW[-alpha[1]]^2
                 + 2*PBW[alpha[1]]*PBW[alphacheck[1]]^2*PBW[-alpha[1]]
                 + 1/4*PBW[alphacheck[1]]^4 - PBW[alphacheck[1]]^3
                 - 4*PBW[alpha[1]]*PBW[-alpha[1]] + 2*PBW[alphacheck[1]]
                sage: all(g * C4 == C4 * g for g in L.pbw_basis().algebra_generators())
                True

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.Heisenberg(QQ, 2)
                sage: L.casimir_element()
                0

                sage: # needs sage.combinat sage.modules
                sage: g = LieAlgebra(QQ, cartan_type=['D',2])
                sage: U = g.pbw_basis()
                sage: U.casimir_element(2, basis=True)
                [2*PBW[alpha[2]]*PBW[-alpha[2]] + 1/2*PBW[alphacheck[2]]^2 - PBW[alphacheck[2]],
                 2*PBW[alpha[1]]*PBW[-alpha[1]] + 1/2*PBW[alphacheck[1]]^2 - PBW[alphacheck[1]]]

            TESTS::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, cartan_type=['A', 1])
                sage: L.casimir_element(1)
                Traceback (most recent call last):
                ...
                ValueError: invalid order
                sage: 4 * L.casimir_element() == L.casimir_element(force_generic=True)
                True

            .. TODO::

                Use the symmetry of the tensor to reduce the number of
                equations and/or variables to solve.
            """
            if order < 2:
                raise ValueError("invalid order")

            if UEA is None:
                UEA = self.universal_enveloping_algebra()

            B = self.basis()

            if order == 2 and not force_generic and not basis:
                # Special case for the quadratic using the Killing form
                try:
                    K = self.killing_form_matrix().inverse()
                    return UEA.sum(K[i, j] * UEA(x) * UEA(y) for i, x in enumerate(B)
                                   for j, y in enumerate(B) if K[i, j])
                except (ValueError, TypeError, ZeroDivisionError):
                    # fall back to finding solutions to the system of equations
                    pass

            keys = self.get_order()
            dim = len(keys)
            s_coeffs = dict(self.structure_coefficients())
            for k in list(s_coeffs.keys()):
                s_coeffs[k[1], k[0]] = -s_coeffs[k]

            # setup the equations
            from sage.matrix.constructor import matrix
            from itertools import product
            eqns = matrix.zero(self.base_ring(), dim**(order+1), dim**order, sparse=True)
            for ii, p in enumerate(product(range(dim), repeat=order+1)):
                i = p[0]
                a = keys[i]
                for j, b in enumerate(keys):
                    if (a, b) not in s_coeffs:
                        continue
                    sc_val = s_coeffs[a, b]
                    for k in range(order):
                        c = keys[p[k+1]]
                        if not sc_val[c]:
                            continue
                        pp = list(p[1:])
                        pp[k] = j
                        jj = sum(dim**m * pp[m] for m in range(order))
                        eqns[ii, jj] += sc_val[c]

            ker = eqns.right_kernel()
            if ker.dimension() == 0:
                return self.zero()

            del eqns  # no need to hold onto the matrix

            def to_prod(vec, index):
                coeff = vec[index]
                p = [0] * order
                base = dim ** (order-1)
                for i in range(order):
                    p[i] = index // base
                    index %= base
                    base //= dim
                p.reverse()
                return coeff * UEA.prod(UEA(B[keys[i]]) for i in p)

            tens = ker.basis()

            if not basis:
                vec = tens[0]
                return UEA.sum(to_prod(vec, index) for index in vec.support())

            return [UEA.sum(to_prod(vec, index) for index in vec.support())
                    for vec in tens]

        def faithful_representation(self, algorithm=None):
            r"""
            Return a faithful representation of ``self``.

            By Ado's and Iwasawa's theorems, every finite dimensional
            Lie algebra has a faithful finite dimensional representation.

            INPUT:

            - ``algorithm`` -- one of the following depending on the
              classification of the Lie algebra:

              Nilpotent:

              * ``'regular'`` -- use the universal enveloping algebra quotient
                :class:`~sage.algebras.lie_algebras.representation.FaithfulRepresentationNilpotentPBW`
              * ``'minimal'`` -- construct the minimal representation (for
                precise details, see the documentation of
                :class:`~sage.algebras.lie_algebras.representation.FaithfulRepresentationNilpotentPBW`)

              Solvable:

              * Not implemented

              Semisimple:

              * Not implemented

              General case

              * ``'generic'`` -- generic algorithm (only implemented currently
                for positive characteristic)

            Note that the algorithm for any more generic cases can be used
            in the specialized cases. For instance, using ``'generic'`` for
            any Lie algebra (e.g., even if nilpotent) will use the generic
            implementation.

            EXAMPLES::

                sage: H2 = lie_algebras.Heisenberg(QQ, 2)
                sage: H2.is_nilpotent()
                True
                sage: F = H2.faithful_representation(); F
                Faithful 16 dimensional representation of
                 Heisenberg algebra of rank 2 over Rational Field
                sage: M = H2.faithful_representation(algorithm='minimal'); M
                Minimal faithful representation of
                 Heisenberg algebra of rank 2 over Rational Field
                sage: M.dimension()
                4
                sage: H2.faithful_representation(algorithm='invalid')
                Traceback (most recent call last):
                ...
                ValueError: invalid algorithm 'invalid'

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: L.is_nilpotent()
                False
                sage: L.is_solvable()
                True
                sage: L.faithful_representation()
                Traceback (most recent call last):
                ...
                NotImplementedError: only implemented for nilpotent Lie algebras

                sage: sl3 = LieAlgebra(QQ, cartan_type=['A', 2])
                sage: sl3.is_semisimple()
                True
                sage: sl3.faithful_representation()
                Traceback (most recent call last):
                ...
                NotImplementedError: only implemented for nilpotent Lie algebras
            """
            if self.is_nilpotent():
                if algorithm is None:
                    algorithm = "regular"
                if algorithm == "regular":
                    from sage.algebras.lie_algebras.representation import FaithfulRepresentationNilpotentPBW
                    return FaithfulRepresentationNilpotentPBW(self, minimal=False)
                if algorithm == "minimal":
                    from sage.algebras.lie_algebras.representation import FaithfulRepresentationNilpotentPBW
                    return FaithfulRepresentationNilpotentPBW(self, minimal=True)
            if algorithm is None or algorithm == "generic":
                if self.base_ring().characteristic() > 0:
                    from sage.algebras.lie_algebras.representation import FaithfulRepresentationPBWPosChar
                    return FaithfulRepresentationPBWPosChar(self)
                raise NotImplementedError("only implemented for nilpotent Lie algebras")
            raise ValueError("invalid algorithm '{}'".format(algorithm))

    class ElementMethods:
        def adjoint_matrix(self, sparse=False): # In #11111 (more or less) by using matrix of a morphism
            """
            Return the matrix of the adjoint action of ``self``.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.an_element().adjoint_matrix()
                [0 0 0]
                [0 0 0]
                [0 0 0]
                sage: L.an_element().adjoint_matrix(sparse=True).is_sparse()
                True

            ::

                sage: # needs sage.combinat sage.modules
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: x.adjoint_matrix()
                [0 1]
                [0 0]
                sage: y.adjoint_matrix()
                [-1  0]
                [ 0  0]

            We verify that this forms a representation::

                sage: # needs sage.combinat sage.modules
                sage: sl3 = lie_algebras.sl(QQ, 3)
                sage: e1, e2 = sl3.e(1), sl3.e(2)
                sage: e12 = e1.bracket(e2)
                sage: E1, E2 = e1.adjoint_matrix(), e2.adjoint_matrix()
                sage: E1 * E2 - E2 * E1 == e12.adjoint_matrix()
                True

            TESTS::

                sage: scoeffs = {('a','d'): {'a':1}, ('a','e'): {'b':-1},
                ....:            ('b','d'): {'b':1}, ('b','e'): {'a':1},
                ....:            ('d','e'): {'c':1}}
                sage: L.<a,b,c,d,e> = LieAlgebra(QQ, scoeffs)
                sage: S = L.solvable_radical()
                sage: elt = S.derived_subalgebra().an_element()
                sage: elt.adjoint_matrix()
                [0 0 0]
                [0 0 0]
                [0 0 0]
            """
            from sage.matrix.constructor import matrix

            P = self.parent()
            basis = P.basis()
            return matrix(self.base_ring(),
                          [P.bracket(self, b).to_vector(sparse=sparse) for b in basis],
                          sparse=sparse).transpose()

        def to_vector(self, sparse=False, order=None):
            r"""
            Return the vector in ``g.module()`` corresponding to the
            element ``self`` of ``g`` (where ``g`` is the parent of
            ``self``).

            Implement this if you implement ``g.module()``.
            See :meth:`sage.categories.lie_algebras.LieAlgebras.module`
            for how this is to be done.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.an_element().to_vector()
                (0, 0, 0)
                sage: L.an_element().to_vector(sparse=True)
                (0, 0, 0)

                sage: # needs sage.combinat sage.groupssage.modules
                sage: D = DescentAlgebra(QQ, 4).D()
                sage: L = LieAlgebra(associative=D)
                sage: L.an_element().to_vector()
                (1, 1, 1, 1, 1, 1, 1, 1)

            TESTS:

            Check that the error raised agrees with the one
            from ``monomial_coefficients()`` (see :issue:`25007`)::

                sage: # needs sage.combinat sage.modules
                sage: L = lie_algebras.sp(QQ, 4, representation='matrix')
                sage: x = L.an_element()
                sage: x.monomial_coefficients()
                Traceback (most recent call last):
                ...
                NotImplementedError: the basis is not defined
                sage: x.to_vector()
                Traceback (most recent call last):
                ...
                NotImplementedError: the basis is not defined
            """
            mc = self.monomial_coefficients(copy=False)
            if sparse:
                from sage.modules.free_module import FreeModule
                M = FreeModule(self.parent().base_ring(), self.dimension(), sparse=True)
                if order is None:
                    order = {b: i for i,b in enumerate(self.parent()._basis_ordering)}
                return M({order[k]: c for k, c in mc.items()})
            else:
                M = self.parent().module()
                B = M.basis()
                if order is None:
                    order = self.parent()._basis_ordering
                return M.sum(mc[k] * B[i] for i, k in enumerate(order) if k in mc)

        _vector_ = to_vector

    class Subobjects(SubobjectsCategory):
        """
        A category for subalgebras of a finite dimensional Lie algebra
        with basis.
        """
        class ParentMethods:
            @abstract_method
            def ambient(self):
                """
                Return the ambient Lie algebra of ``self``.

                EXAMPLES::

                    sage: # needs sage.combinat sage.modules
                    sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis()
                    sage: L = C.example()
                    sage: a, b, c = L.lie_algebra_generators()
                    sage: S = L.subalgebra([2*a + b, b + c])
                    sage: S.ambient() == L
                    True
                """

            @abstract_method
            def basis_matrix(self):
                """
                Return the basis matrix of ``self``.

                EXAMPLES::

                    sage: # needs sage.combinat sage.modules
                    sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis()
                    sage: L = C.example()
                    sage: a, b, c = L.lie_algebra_generators()
                    sage: S = L.subalgebra([2*a + b, b + c])
                    sage: S.basis_matrix()
                    [   1    0 -1/2]
                    [   0    1    1]
                """

            def reduce(self, X):
                r"""
                Reduce an element of the ambient Lie algebra modulo the
                ideal ``self``.

                INPUT:

                - ``X`` -- an element of the ambient Lie algebra

                OUTPUT:

                An element `Y` of the ambient Lie algebra that is contained
                in a fixed complementary submodule `V` to ``self`` such that
                `X = Y` mod ``self``.

                When the base ring of ``self`` is a field, the complementary
                submodule `V` is spanned by the elements of the basis that
                are not the leading supports of the basis of ``self``.

                EXAMPLES:

                An example reduction in a 6 dimensional Lie algebra::

                    sage: sc = {('a','b'): {'d': 1}, ('a','c'): {'e': 1},
                    ....:       ('b','c'): {'f': 1}}
                    sage: L.<a,b,c,d,e,f> = LieAlgebra(QQ, sc)
                    sage: I = L.ideal(c)
                    sage: I.reduce(a + b + c + d + e + f)
                    a + b + d

                The reduction of an element is zero if and only if the
                element belongs to the subalgebra::

                    sage: I.reduce(c + e)
                    0
                    sage: c + e in I
                    True

                Over non-fields, the complementary submodule may not be spanned
                by a subset of the basis of the ambient Lie algebra::

                    sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
                    sage: I = L.ideal(Y)
                    sage: I.basis()
                    Family (Y, 3*Z)
                    sage: I.reduce(3*Z)
                    0
                    sage: I.reduce(Y + 14*Z)
                    2*Z
                """
                R = self.base_ring()
                from sage.categories.fields import Fields
                is_field = R in Fields()
                for Y in self.basis():
                    Y = self.lift(Y)
                    k, c = Y.leading_item(key=self._order)

                    if is_field:
                        X -= (X[k] / c) * Y
                    else:
                        try:
                            q, _ = X[k].quo_rem(c)
                            X -= q * Y
                        except AttributeError:
                            break

                return X
