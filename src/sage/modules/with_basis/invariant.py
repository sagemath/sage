# sage.doctest: needs sage.groups
r"""
Invariant modules
"""

# ****************************************************************************
#       Copyright (C) 2021 Trevor K. Karn <karnx018 at umn.edu>
#                     2021 Travis Scrimshaw
#                     2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator
from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.modules.with_basis.representation import Representation
from sage.categories.finitely_generated_semigroups import FinitelyGeneratedSemigroups
from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModulesWithBasis
from sage.sets.family import Family
from sage.matrix.constructor import Matrix
from sage.libs.gap.libgap import libgap


class FiniteDimensionalInvariantModule(SubmoduleWithBasis):
    r"""
    The invariant submodule under a semigroup action.

    When a semigroup `S` acts on a module `M`, the invariant module is the
    set of elements `m \in M` such that `s \cdot m = m` for all `s \in S`:

    .. MATH::

        M^S := \{m \in M : s \cdot m = m,\, \forall s \in S \}.

    INPUT:

    - ``M`` -- a module in the category of
      :class:`~sage.categories.finite_dimensional_modules_with_basis.FiniteDimensionalModulesWithBasis`

    - ``S`` -- a semigroup in the category of
      :class:`~sage.categories.finitely_generated_semigroups.FinitelyGeneratedSemigroups`

    - ``action`` -- (default: ``operator.mul``) the action of ``S`` on ``M``

    - ``side`` -- (default: ``'left'``) the side on which ``S`` acts

    EXAMPLES:

    First, we create the invariant defined by the cyclic group action on the
    free module with basis `\{1,2,3\}`::

        sage: G = CyclicPermutationGroup(3)
        sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
        sage: action = lambda g, m: M.monomial(g(m))  # cyclically permute coordinates

    In order to give the module an action of ``G``, we create a
    :class:`~sage.modules.with_basis.representation.Representation`::

        sage: from sage.modules.with_basis.representation import Representation
        sage: R = Representation(G, M, action)
        sage: I = R.invariant_module()

    Then we can lift the basis from the invariant to the original module::

        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]

    The we could also have the action be a right-action, instead of the
    default left-action::

        sage: def rt_action(g, m): return M.monomial(g(m))  # cyclically permute coordinates
        sage: R = Representation(G, M, rt_action, side='right')  # same as last but on right
        sage: g = G.an_element(); g
        (1,2,3)
        sage: r = R.an_element(); r
        2*M[1] + 2*M[2] + 3*M[3]
        sage: R.side()
        'right'

    So now we can see that multiplication with ``g`` on the right sends
    ``M[1]`` to ``M[2]`` and so on::

        sage: r * g
        3*M[1] + 2*M[2] + 2*M[3]
        sage: I = R.invariant_module()
        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]

    Now we will take the regular representation of the symmetric group on
    three elements to be the module, and compute its invariant submodule::

        sage: G = SymmetricGroup(3)
        sage: R = G.regular_representation(QQ)
        sage: I = R.invariant_module()
        sage: [I.lift(b).to_vector() for b in I.basis()]
        [(1, 1, 1, 1, 1, 1)]

    We can also check the scalar multiplication by elements of the base ring
    (for this example, the rational field)::

        sage: [I.lift(3*b).to_vector() for b in I.basis()]
        [(3, 3, 3, 3, 3, 3)]

    A more subtle example is the invariant submodule of a skew-commutative
    module, for example the exterior algebra `E[x_0,x_1,x_2]` generated
    by three elements::

        sage: G = CyclicPermutationGroup(3)
        sage: M = algebras.Exterior(QQ, 'x', 3)
        sage: def cyclic_ext_action(g, m):
        ....:     # cyclically permute generators
        ....:     return M.prod([M.monomial(FrozenBitset([g(j+1)-1])) for j in m])

    If you care about being able to exploit the algebra structure of the
    exterior algebra (i.e. if you want to multiply elements together), you
    should make sure the representation knows it is also an algebra with
    the semigroup action being by algebra endomorphisms::

        sage: cat = Algebras(QQ).WithBasis().FiniteDimensional()
        sage: R = Representation(G, M, cyclic_ext_action, category=cat)
        sage: I = R.invariant_module()

    We can express the basis in the ambient algebra (`E[x_0,x_1,x_2]`)::

        sage: [I.lift(b) for b in I.basis()]
        [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]

    or we can express the basis intrinsicallly to the invariant ``I``::

        sage: B = I.basis()
        sage: m = 3*B[0] + 2*B[1] + 7*B[3]

    This lifts to the exterior algebra::

        sage: I.lift(m)
        3 + 2*x0 + 7*x0*x1*x2 + 2*x1 + 2*x2

    We can also check using the invariant element ``m`` that arithmetic works::

        sage: m^2
        9*B[0] + 12*B[1] + 42*B[3]
        sage: m+m
        6*B[0] + 4*B[1] + 14*B[3]

    To see the actual elements expressed in the exterior algebra, we lift them
    again::

        sage: I.lift(m+m)
        6 + 4*x0 + 14*x0*x1*x2 + 4*x1 + 4*x2
        sage: 7*m
        21*B[0] + 14*B[1] + 49*B[3]
        sage: I.lift(7*m)
        21 + 14*x0 + 49*x0*x1*x2 + 14*x1 + 14*x2

    The classic example of an invariant module is the module of symmetric
    functions, which is the invariant module of polynomials whose variables
    are acted upon by permutation. We can create a module isomorphic to the
    homogeneous component of a a polynomial ring in `n` variable of a fixed
    degree `d` by looking at  weak  compositions of `d` of length `n`, which
    we consider as the exponent vector. For example, `x^2yz \in \QQ[x,y,z]`
    would have the exponent vector `(2,1,1)`. The vector `(2,1,1)` is a
    weak composition of `4`, with length `3`, and so we can think of it as
    being in the degree-`4` homogeneous component of a polynomial ring
    in three variables::

        sage: C = IntegerVectors(4, length=3, min_part=0)  # representing degree-4 monomials
        sage: M = CombinatorialFreeModule(QQ, C)  # isomorphic to deg-4 homog. polynomials
        sage: G = SymmetricGroup(3)
        sage: def perm_action(g,x): return M.monomial(C(g(list(x))))
        sage: perm_action(G((1,2,3)), C([4,3,2]))
        B[[3, 2, 4]]
        sage: R = Representation(G, M, perm_action)
        sage: I = R.invariant_module()
        sage: [I.lift(b) for b in I.basis()]
        [B[[0, 0, 4]] + B[[0, 4, 0]] + B[[4, 0, 0]],
         B[[0, 1, 3]] + B[[0, 3, 1]] + B[[1, 0, 3]]
         + B[[1, 3, 0]] + B[[3, 0, 1]] + B[[3, 1, 0]],
         B[[0, 2, 2]] + B[[2, 0, 2]] + B[[2, 2, 0]],
         B[[1, 1, 2]] + B[[1, 2, 1]] + B[[2, 1, 1]]]

    These are the monomial symmetric functions, which are a well-known
    basis for the symmetric functions. For comparison::

        sage: Sym = SymmetricFunctions(QQ)
        sage: m = Sym.monomial()
        sage: [m[mu].expand(3) for mu in Partitions(4)]
        [x0^4 + x1^4 + x2^4,
         x0^3*x1 + x0*x1^3 + x0^3*x2 + x1^3*x2 + x0*x2^3 + x1*x2^3,
         x0^2*x1^2 + x0^2*x2^2 + x1^2*x2^2,
         x0^2*x1*x2 + x0*x1^2*x2 + x0*x1*x2^2,
         0]

    .. NOTE::

        The current implementation works when `S` is a finitely-generated
        semigroup, and when `M` is a finite-dimensional free module with
        a distinguished basis.

    .. TODO::

        Extend this to have multiple actions, including actions on both sides.

    .. TODO::

        Extend when `M` does not have a basis and `S` is a permutation
        group using:

        - :arxiv:`0812.3082`
        - https://www.dmtcs.org/pdfpapers/dmAA0123.pdf
    """
    def __init__(self, M, S, action=operator.mul, side='left', *args, **kwargs):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: R = G.regular_representation()
            sage: I = R.invariant_module()
            sage: TestSuite(I).run()

        TESTS::

            sage: G = GroupExp()(QQ) # a group that is not finitely generated
            sage: M = CombinatorialFreeModule(QQ, [1,2,3])
            sage: def on_basis(g,m): return M.monomial(m)  # trivial rep'n
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, M, on_basis)
            sage: R.invariant_module()
            Traceback (most recent call last):
            ...
            ValueError: Multiplicative form of Rational Field is not finitely generated
        """
        if S not in FinitelyGeneratedSemigroups():
            raise ValueError(f"{S} is not finitely generated")
        if M not in FiniteDimensionalModulesWithBasis:
            raise ValueError(f"{M} is not a finite dimensional module with a distinguished basis")

        if side == "left":
            def _invariant_map(g, x):
                return action(g, x) - x
        elif side == "right":
            def _invariant_map(g, x):
                return action(x, g) - x
        else:
            raise ValueError("side must either be 'left' or 'right'")

        self._side = side
        self._action = action
        self._semigroup = S

        category = kwargs.pop("category", M.category().Subobjects())

        # Give the intersection of kernels of the map `s*x-x` to determine when
        # `s*x = x` for all generators `s` of `S`
        basis = M.annihilator_basis(S.gens(), action=_invariant_map, side='left')

        super().__init__(Family(basis),
                         support_order=M._compute_support_order(basis),
                         ambient=M,
                         unitriangular=False,
                         category=category,
                         *args, **kwargs)

    def construction(self):
        r"""
        Return the functorial construction of ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: R = G.regular_representation(); R
            Left Regular Representation of Cyclic group of order 3 as a permutation group over Integer Ring
            sage: I = R.invariant_module()
            sage: I.construction()
            (EquivariantSubobjectConstructionFunctor,
            Left Regular Representation of Cyclic group of order 3 as a permutation group over Integer Ring)
        """
        from sage.categories.pushout import EquivariantSubobjectConstructionFunctor
        return (EquivariantSubobjectConstructionFunctor(self._semigroup,
                                                        self._action,
                                                        self._side),
                self.ambient())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: R = G.trivial_representation()
            sage: R.invariant_module()
            (Cyclic group of order 3 as a permutation group)-invariant submodule of
             Trivial representation of Cyclic group of order 3 as a permutation group over Integer Ring

            sage: G = CyclicPermutationGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: action = lambda g, m: M.monomial(g(m))  # cyclically permute coordinates
            sage: M.invariant_module(G, action_on_basis=action)
            (Cyclic group of order 3 as a permutation group)-invariant submodule of
             Free module generated by {1, 2, 3} over Rational Field
        """
        M = self._ambient
        if isinstance(self._ambient, Representation):
            M = M._module
        return f"({self._semigroup})-invariant submodule of {M}"

    def _latex_(self):
        r"""
        Return a latex representaion of ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: R = G.algebra(QQ)
            sage: latex(R.invariant_module(G))
            \left( \Bold{Q}[\langle (1,2,3) \rangle] \right)^{\langle (1,2,3) \rangle}
        """
        M = self._ambient
        if isinstance(self._ambient, Representation):
            M = M._module
        from sage.misc.latex import latex
        return "\\left( {} \\right)^{{{}}}".format(latex(M), latex(self._semigroup))

    def _test_invariant(self, **options):
        """
        Check (on some elements) that ``self`` is invariant.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: def action(g, x): return M.monomial(g(x))
            sage: I = M.invariant_module(G, action_on_basis=action)
            sage: I._test_invariant()

            sage: G = SymmetricGroup(10)
            sage: M = CombinatorialFreeModule(QQ, list(range(1,11)), prefix='M')
            sage: def action(g, x): return M.monomial(g(x))
            sage: I = M.invariant_module(G, action_on_basis=action)
            sage: I._test_invariant(max_runs=10)
        """
        tester = self._tester(**options)
        X = tester.some_elements()
        L = []
        max_len = tester._max_runs

        # FIXME: This is max_len * dim number of runs!!!
        for i, x in enumerate(self._semigroup):
            L.append(x)
            if i >= max_len:
                break

        for x in L:
            for elt in X:
                lifted = self.lift(elt)
                if self._side == 'left':
                    tester.assertEqual(self._action(x, lifted), lifted)
                else:
                    tester.assertEqual(self._action(lifted, x), lifted)

    def semigroup(self):
        r"""
        Return the semigroup `S` whose action ``self`` is invariant under.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: def action(g,x): return M.monomial(g(x))
            sage: I = M.invariant_module(G, action_on_basis=action)
            sage: I.semigroup()
            Symmetric group of order 3! as a permutation group
        """
        return self._semigroup

    semigroup_representation = SubmoduleWithBasis.ambient

    class Element(SubmoduleWithBasis.Element):
        def _mul_(self, other):
            r"""
            Multiply ``self`` and ``other``.

            EXAMPLES:

            In general, there is not a well defined multiplication between
            two elements of a given module, but there is a multiplication
            with scalars::

                sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M');
                sage: G = CyclicPermutationGroup(3); G.rename('G')
                sage: g = G.an_element(); g
                (1,2,3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: R = Representation(G,M,lambda g,x:M.monomial(g(x))); R.rename('R')
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [M[1] + M[2] + M[3]]
                sage: v = B[0]
                sage: v*v
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'R' and 'R'
                sage: (1/2) * v
                1/2*B[0]
                sage: v * (1/2)
                1/2*B[0]
                sage: R.rename()  # reset name

            Sometimes, the module is also a ring. To ensure the multiplication
            works as desired, we should be sure to pass the correct category to
            the :class:`~sage.modules.with_basis.representation.Representation`.
            In the following example, we use the exterior algebra over `\QQ`
            with three generators, which is in the category of finite
            dimensional `\QQ`-algebras with a basis::

                sage: G = CyclicPermutationGroup(3); G.rename('G')
                sage: M = algebras.Exterior(QQ, 'x', 3)
                sage: def on_basis(g,m): return M.prod([M.monomial(FrozenBitset([g(j+1)-1])) for j in m])  # cyclically permute generators
                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional(), side='right')
                sage: I = R.invariant_module(); I.rename('I')
                sage: B = I.basis()
                sage: v = B[0] + 2*B[1]; I.lift(v)
                1 + 2*x0 + 2*x1 + 2*x2
                sage: w = B[2]; I.lift(w)
                x0*x1 - x0*x2 + x1*x2
                sage: v * w
                B[2] + 6*B[3]
                sage: I.lift(v*w)
                x0*x1 + 6*x0*x1*x2 - x0*x2 + x1*x2
                sage: w * v
                B[2] + 6*B[3]
                sage: (1/2) * v
                1/2*B[0] + B[1]
                sage: w * (1/2)
                1/2*B[2]
                sage: g = G((1,3,2))
                sage: v * g
                B[0] + 2*B[1]
                sage: w * g
                B[2]
                sage: g * v
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'G' and 'I'
                sage: I.rename()  # reset name

                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
                sage: I = R.invariant_module(); I.rename('I')
                sage: B = I.basis()
                sage: v = B[0] + 2*B[1]; I.lift(v)
                1 + 2*x0 + 2*x1 + 2*x2
                sage: w = B[2]; I.lift(w)
                x0*x1 - x0*x2 + x1*x2
                sage: v * w
                B[2] + 6*B[3]
                sage: I.lift(v*w)
                x0*x1 + 6*x0*x1*x2 - x0*x2 + x1*x2
                sage: w * v
                B[2] + 6*B[3]
                sage: (1/2) * v
                1/2*B[0] + B[1]
                sage: w * (1/2)
                1/2*B[2]
                sage: g = G((1,3,2))
                sage: v * v
                B[0] + 4*B[1]
                sage: g * w
                B[2]
                sage: v * g
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'I' and 'G'
                sage: G.rename(); I.rename()  # reset names
            """
            P = self.parent()
            return P.retract(P.lift(self) * P.lift(other))

        def _acted_upon_(self, scalar, self_on_left=False):
            """
            EXAMPLES::

                sage: G = CyclicPermutationGroup(3)
                sage: g = G.an_element(); g
                (1,2,3)
                sage: M = CombinatorialFreeModule(QQ, [1,2,3])
                sage: from sage.modules.with_basis.representation import Representation
                sage: R = Representation(G, M, lambda g,x: M.monomial(g(x)))
                sage: I = R.invariant_module()
                sage: [b._acted_upon_(G((1,3,2))) for b in I.basis()]
                [B[0]]
                sage: v = I.an_element(); v
                2*B[0]
                sage: g * v
                2*B[0]
                sage: [g * v for g in G.list()]
                [2*B[0], 2*B[0], 2*B[0]]


                sage: E = algebras.Exterior(QQ, 'x', 3)
                sage: def on_basis(g,m): return E.prod([E.monomial(FrozenBitset([g(j+1)-1])) for j in m])  # cyclically permute generators
                sage: R = Representation(G, E, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]
                sage: [[g*b for g in G] for b in B]
                [[B[0], B[0], B[0]],
                 [B[1], B[1], B[1]],
                 [B[2], B[2], B[2]],
                 [B[3], B[3], B[3]]]
                sage: 3 * I.basis()[0]
                3*B[0]
                sage: 3*B[0] + B[1]*2
                3*B[0] + 2*B[1]

                sage: R = G.regular_representation(QQ)
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [() + (1,2,3) + (1,3,2)]
                sage: B[0]._acted_upon_(G((1,3,2)))
                B[0]
                sage: B[0]._acted_upon_(G((1,3,2)), self_on_left=True) is None
                True

                sage: R = G.regular_representation(QQ, side='right')
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [() + (1,2,3) + (1,3,2)]
                sage: g = G((1,3,2))
                sage: B[0]._acted_upon_(g, self_on_left=True)
                B[0]
                sage: B[0]._acted_upon_(g, self_on_left=False) is None
                True

                sage: R = Representation(G, M, lambda g,x: M.monomial(g(x)), side='right')
                sage: I = R.invariant_module()
                sage: v = I.an_element(); v
                2*B[0]
                sage: v * g
                2*B[0]
                sage: [v * g for g in G.list()]
                [2*B[0], 2*B[0], 2*B[0]]
                sage: [b._acted_upon_(G((1,3,2)), self_on_left=True) for b in I.basis()]
                [B[0]]

                sage: def on_basis(g,m): return E.prod([E.monomial(FrozenBitset([g(j+1)-1])) for j in m])  # cyclically permute generators
                sage: R = Representation(G, E, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional(), side='right')
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]
                sage: [[b * g for g in G] for b in B]
                [[B[0], B[0], B[0]],
                 [B[1], B[1], B[1]],
                 [B[2], B[2], B[2]],
                 [B[3], B[3], B[3]]]
                sage: 3 * B[0] + B[1] * 2
                3*B[0] + 2*B[1]
            """
            if scalar in self.parent()._semigroup and self_on_left == (self.parent()._side == 'right'):
                return self
            return super()._acted_upon_(scalar, self_on_left)


class FiniteDimensionalTwistedInvariantModule(SubmoduleWithBasis):
    r"""
    Construct the `\chi`-twisted invariant submodule of `M`.

    When a group `G` acts on a module `M`, the `\chi`-*twisted invariant
    submodule* of `M` is the isotypic component of the representation `M`
    corresponding to the irreducible character `\chi`.

    For more information, see [Sta1979]_.

    INPUT:

    - ``M`` -- a module in the category of
      :class:`~sage.categories.finite_dimensional_modules_with_basis.FiniteDimensionalModulesWithBasis`
      and whose base ring contains all the values passed to ``chi`` and `1/|G|`

    - ``G`` -- a finitely generated group

    - ``chi`` -- list/tuple of the character values of the irreducible representation
      onto which you want to project. The order of values of `chi` must
      agree with the order of ``G.conjugacy_classes()``

    - ``action`` -- (default: ``operator.mul``) the action of ``G`` on ``M``

    - ``side`` -- (default: ``'left'``) the side on which ``G`` acts

    .. WARNING::

        The current implementation does not check if ``chi`` is irreducible.
        Passing character values of non-irreducible representations may lead
        to mathematically incorrect results.

    EXAMPLES:

    Suppose that the symmetric group `S_3` acts on a four dimensional
    vector space by permuting the first three coordinates only::

        sage: M = CombinatorialFreeModule(QQ, [1,2,3,4], prefix='M')
        sage: G = SymmetricGroup(3)
        sage: action = lambda g,x: M.term(g(x))

    The trivial representation corresponds to the usual invariant module,
    so trying to create the twisted invariant module when there is no twist
    returns a :class:`~sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule`::

        sage: chi = ClassFunction(G, (1,1,1))
        sage: T = M.twisted_invariant_module(G, chi, action_on_basis=action)
        sage: type(T)
        <class 'sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule_with_category'>

    In this case, there are two copies of the trivial representation, one
    coming from the first three coordinates and the other coming from the
    fact that `S_3` does not touch the fourth coordinate::

        sage: T.basis()
        Finite family {0: B[0], 1: B[1]}
        sage: [T.lift(b) for b in T.basis()]
        [M[1] + M[2] + M[3], M[4]]

    The character values of the standard representation are `2,0,-1`::

        sage: chi = ClassFunction(G, [2,0,-1])
        sage: T = M.twisted_invariant_module(G, chi, action_on_basis=action)
        sage: type(T)
        <class 'sage.modules.with_basis.invariant.FiniteDimensionalTwistedInvariantModule_with_category'>
        sage: T.basis()
        Finite family {0: B[0], 1: B[1]}
        sage: [T.lift(b) for b in T.basis()]
        [M[1] - M[3], M[2] - M[3]]

    The permutation representation is the direct sum of the standard
    representation with the trivial representation, and the action on the
    basis element ``B[4]`` is itself a copy of the trivial representation,
    so the sign representation does not appear in the decomposition::

        sage: T = M.twisted_invariant_module(G, [1,-1,1], action_on_basis=action)
        sage: T.basis()
        Finite family {}

    We can also get two copies of the standard representation by looking at
    two copies of the permutation representation, found by reduction modulo
    three on the indices of a six-dimensional module::

        sage: M = CombinatorialFreeModule(QQ, [0,1,2,3,4,5], prefix='M')
        sage: action = lambda g,x: M.term(g(x%3 + 1)-1 + (x>=3)*3)
        sage: T = M.twisted_invariant_module(G, [2,0,-1], action_on_basis=action)
        sage: T.basis()
        Finite family {0: B[0], 1: B[1], 2: B[2], 3: B[3]}
        sage: [T.lift(b) for b in T.basis()]
        [M[0] - M[2], M[1] - M[2], M[3] - M[5], M[4] - M[5]]

        sage: T = M.twisted_invariant_module(G, [1,1,1], action_on_basis=action)
        sage: T.basis()
        Finite family {0: B[0], 1: B[1]}
        sage: [T.lift(b) for b in T.basis()]
        [M[0] + M[1] + M[2], M[3] + M[4] + M[5]]

    There are still no copies of the sign representation::

        sage: T = M.twisted_invariant_module(G, [1,-1,1], action_on_basis=action)
        sage: T.basis()
        Finite family {}

    The trivial representation also contains no copies of the sign
    representation::

        sage: R = G.trivial_representation(QQ)
        sage: T = R.twisted_invariant_module([1,-1,1])
        sage: T.basis()
        Finite family {}

    The regular representation contains two copies of the standard
    representation and one copy each of the trivial and the sign::

        sage: R = G.regular_representation(QQ)
        sage: std = R.twisted_invariant_module([2,0,-1])
        sage: std.basis()
        Finite family {0: B[0], 1: B[1], 2: B[2], 3: B[3]}
        sage: [std.lift(b) for b in std.basis()]
        [() - (1,2,3), -(1,2,3) + (1,3,2), (2,3) - (1,2), -(1,2) + (1,3)]

        sage: triv = R.twisted_invariant_module([1,1,1])
        sage: triv.basis()
        Finite family {0: B[0]}
        sage: [triv.lift(b) for b in triv.basis()]
        [() + (2,3) + (1,2) + (1,2,3) + (1,3,2) + (1,3)]

        sage: sgn = R.twisted_invariant_module([1,-1,1])
        sage: sgn.basis()
        Finite family {0: B[0]}
        sage: [sgn.lift(b) for b in sgn.basis()]
        [() - (2,3) - (1,2) + (1,2,3) + (1,3,2) - (1,3)]

    For the next example, we construct a twisted invariant by the character
    for the 2 dimensional representation of `S_3` on the natural action on
    the exterior algebra. While `S_3` acts by automorphisms, the twisted
    invariants do not form an algebra in this case::

        sage: G = SymmetricGroup(3); G.rename('S3')
        sage: E = algebras.Exterior(QQ, 'x', 3); E.rename('E')
        sage: def action(g,m): return E.prod([E.monomial(FrozenBitset([g(j+1)-1])) for j in m])
        sage: from sage.modules.with_basis.representation import Representation
        sage: EA = Representation(G, E, action, category=Algebras(QQ).WithBasis().FiniteDimensional())
        sage: T = EA.twisted_invariant_module([2,0,-1])
        sage: t = T.an_element(); t
        2*B[0] + 2*B[1] + 3*B[2]

    We can still get meaningful information about the product
    by taking the product in the ambient space::

        sage: T.lift(t) * T.lift(t)
        -36*x0*x1*x2

    We can see this does not lie in this twisted invariant algebra::

        sage: T.retract(T.lift(t) * T.lift(t))
        Traceback (most recent call last):
        ...
        ValueError: -36*x0*x1*x2 is not in the image

        sage: [T.lift(b) for b in T.basis()]
        [x0 - x2, x1 - x2, x0*x1 - x1*x2, x0*x2 + x1*x2]

    It happens to be in the trivial isotypic component (equivalently in
    the usual invariant algebra) but Sage does not know this.

    ::

        sage: G.rename(); E.rename()  # reset the names

    .. TODO::

        - Replace ``G`` by ``S`` in :class:`~sage.categories.finitely_generated_semigroups.FinitelyGeneratedSemigroups`
        - Allow for ``chi`` to be a :class:`~sage.modules.with_basis.representation.Representation`
        - Add check for irreducibility of ``chi``
    """

    @staticmethod
    def __classcall_private__(cls, M, G, chi,
                              action=operator.mul, side='left', **kwargs):
        r"""
        TESTS:

        Check that it works for lists::

            sage: M = CombinatorialFreeModule(QQ, [1,2,3])
            sage: G = SymmetricGroup(3)
            sage: def action(g,x): return M.term(g(x))
            sage: T = M.twisted_invariant_module(G, [2,0,-1], action_on_basis=action)

        Check that it works for tuples::

            sage: T2 = M.twisted_invariant_module(G, (2,0,-1), action_on_basis=action)
            sage: T is T2
            True

        Check that it works for class functions::

            sage: chi = ClassFunction(G, [2,0,-1])
            sage: T3 = M.twisted_invariant_module(G, chi, action_on_basis=action)
            sage: T is T3
            True

        Check that it works when the character values are not an instance of
        :class:`~sage.rings.integer.Integer`::

            sage: chi = [QQ(2), QQ(0), QQ(-1)]
            sage: T4 = M.twisted_invariant_module(G, chi, action_on_basis=action)
            sage: T is T4
            True

        Check that the trivial character returns an instance of
        :class:`~sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule`::

            sage: chi = [1, 1, 1] # check for list
            sage: T = M.twisted_invariant_module(G, chi, action_on_basis=action)
            sage: type(T)
            <class 'sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule_with_category'>

            sage: chi = (1, 1, 1) # check for tuple
            sage: T = M.twisted_invariant_module(G, chi, action_on_basis=action)
            sage: type(T)
            <class 'sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule_with_category'>

            sage: chi = ClassFunction(G, [1,1,1]) # check for class function
            sage: T = M.twisted_invariant_module(G, chi, action_on_basis=action)
            sage: type(T)
            <class 'sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule_with_category'>

        Check the :exc:`ValueError`::

            sage: T = M.twisted_invariant_module(G, "ichigo", action_on_basis=action)
            Traceback (most recent call last):
            ...
            ValueError: chi must be a list/tuple or a class function of the group G
        """
        from sage.groups.class_function import ClassFunction, ClassFunction_libgap

        if isinstance(chi, (list, tuple)):
            chi = ClassFunction(G, libgap(chi))
        elif not isinstance(chi, ClassFunction_libgap):
            raise ValueError("chi must be a list/tuple or a class function of the group G")

        try:
            is_trivial = all(chi(next(iter(conj))) == 1 for conj in G.conjugacy_classes())
        except AttributeError:  # to handle ReflectionGroups
            is_trivial = all(chi(G(next(iter(conj)))) == 1 for conj in G.conjugacy_classes())

        if is_trivial:
            action_on_basis = kwargs.pop('action_on_basis', None)
            if action_on_basis is not None:
                return M.invariant_module(G, action_on_basis=action_on_basis)
            return M.invariant_module(G, action=action)

        return super().__classcall__(cls, M, G, chi, action=operator.mul,
                                     side='left', **kwargs)

    def __init__(self, M, G, chi, action=operator.mul, side='left', **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES:

        As a first example we will consider the permutation representation
        of `S_3`::

            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M');
            sage: G = SymmetricGroup(3); G.conjugacy_classes()
            [Conjugacy class of cycle type [1, 1, 1] in Symmetric group of order 3! as a permutation group,
             Conjugacy class of cycle type [2, 1] in Symmetric group of order 3! as a permutation group,
             Conjugacy class of cycle type [3] in Symmetric group of order 3! as a permutation group]
            sage: from sage.groups.class_function import ClassFunction
            sage: chi = ClassFunction(G, [2,0,-1]) # the standard representation character values
            sage: def action(g,x): return M.term(g(x))
            sage: import __main__
            sage: __main__.action = action
            sage: T = M.twisted_invariant_module(G, chi, action_on_basis=action)
            sage: TestSuite(T).run()

        We know that the permutation representation decomposes as a direct
        sum of one copy of the standard representation which is two-dimensional
        and one copy of the trivial representation::

            sage: T.basis()
            Finite family {0: B[0], 1: B[1]}
            sage: [T.lift(b) for b in T.basis()]
            [M[1] - M[3], M[2] - M[3]]
        """

        if G not in FinitelyGeneratedSemigroups():
            raise ValueError(f"{G} is not finitely generated")
        if M not in FiniteDimensionalModulesWithBasis:
            raise ValueError(f"{M} is not a finite dimensional module with a distinguished basis")

        self._chi = chi
        self._group = G
        self._action = action
        self._side = side

        # define a private action to deal
        # with sidedness issues in the action.
        if side == 'left':
            self.__sided_action__ = action
        elif side == 'right':
            # flip the sides since the second argument
            # to action should be the group element
            def __sided_action__(g, x):
                return action(x, g)
            self.__sided_action__ = __sided_action__
        else:
            raise ValueError("side must either be 'left' or 'right'")

        proj_matrix = Matrix(M.dimension())  # initialize the zero-matrix
        for g in self._group:
            proj_matrix += self._chi(g)*Matrix((self.__sided_action__(g, b)).to_vector() for b in M.basis())

        n = self._chi(self._group.identity())  # chi(1) is the dimension
        g = self._group.order()

        self._projection_matrix = (n/g)*proj_matrix

        self._project_ambient = M.module_morphism(matrix=self._projection_matrix,
                                                  codomain=M)

        category = kwargs.pop("category", M.category().Subobjects())

        # Give the kernel of the map `\pi(x)-x` to determine when `x` lies
        # within the isotypic component of `R`.

        def proj_difference(g, x):
            return self._project_ambient(x) - x

        basis = M.annihilator_basis(M.basis(),
                                    action=proj_difference,
                                    side='left')

        super().__init__(Family(basis),
                         support_order=M._compute_support_order(basis),
                         ambient=M,
                         unitriangular=False,
                         category=category,
                         **kwargs)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: action = lambda g, m: M.monomial(g(m))  # cyclically permute coordinates
            sage: M.twisted_invariant_module(G, [2,0,-1], action_on_basis=action)
            Twist of (Cyclic group of order 3 as a permutation group)-invariant submodule of
             Free module generated by {1, 2, 3} over Rational Field by character [2, 0, -1]
        """
        M = self._ambient
        if isinstance(self._ambient, Representation):
            M = M._module
        return f"Twist of ({self._group})-invariant submodule of {M} by character {self._chi.values()}"

    def project(self, x):
        r"""
        Project ``x`` in the ambient module onto ``self``.

        EXAMPLES:

        The standard representation is the orthogonal complement
        of the trivial representation inside of the permutation
        representation, so the basis for the trivial representation
        projects to `0`::

            sage: M = CombinatorialFreeModule(QQ, [1,2,3]); M.rename('M')
            sage: B = M.basis()
            sage: G = SymmetricGroup(3); G.rename('S3')
            sage: def action(g,x): return M.term(g(x))
            sage: T = M.twisted_invariant_module(G, [2,0,-1], action_on_basis=action)
            sage: m = B[1] + B[2] + B[3]
            sage: parent(m)
            M
            sage: t = T.project(m); t
            0
            sage: parent(t)
            Twist of (S3)-invariant submodule of M by character [2, 0, -1]

            sage: G.rename(); M.rename()  # reset names
        """
        return self.retract(self.project_ambient(x))

    def project_ambient(self, x):
        r"""
        Project ``x`` in the ambient representation onto the submodule of the
        ambient representation to which ``self`` is isomorphic as a module.

        .. NOTE::

            The image of ``self.project_ambient`` is not in ``self`` but
            rather is in ``self.ambient()``.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(QQ, [1,2,3]); M.rename('M')
            sage: B = M.basis()
            sage: G = SymmetricGroup(3); G.rename('S3')
            sage: def action(g,x): return M.term(g(x))
            sage: T = M.twisted_invariant_module(G, [2,0,-1], action_on_basis=action)

        To compare with ``self.project``, we can inspect the parents.
        The image of ``self.project`` is in ``self``, while the image
        of ``self.project_ambient`` is in ``self._ambient``::

            sage: t = T.project(B[1] + B[2] + B[3]); t
            0
            sage: parent(t)
            Twist of (S3)-invariant submodule of M by character [2, 0, -1]
            sage: s = T.project_ambient(B[1] + B[2] + B[3]); s
            0
            sage: parent(s)
            Representation of S3 indexed by {1, 2, 3} over Rational Field

        Note that because of the construction of ``T``, ``self._ambient``
        is an instance of
        :class:`~sage.modules.with_basis.representation.Representation`,
        but you still may pass elements of ``M``, which is an instance of
        :class:`~sage.combinat.free_module.CombinatorialFreeModule`,
        because the underlying ``Representation`` is built off of ``M``
        and we can cannonically construct elements of the ``Representation``
        from elements of ``M``.

        ::

            sage: G.rename(); M.rename()  # reset names
        """
        if (isinstance(self._ambient, Representation)
                and x.parent() is self._ambient._module):
            x = self._ambient._element_constructor_(x)
        return self._project_ambient(x)

    def projection_matrix(self):
        r"""
        Return the matrix defining the projection map from
        the ambient representation onto ``self``.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(QQ, [1,2,3])
            sage: def action(g,x): return(M.term(g(x)))
            sage: G = SymmetricGroup(3)

        If the matrix `A` has columns form a basis for
        the subspace onto which we are trying to project,
        then we can find the projection matrix via the
        formula `P = A (A^T A)^{-1} A^T`. Recall that the
        standard representation twisted invariant has basis
        ``(B[1] - B[3], B[2] - B[3])``, hence::

            sage: A = Matrix([[1,0],[0,1],[-1,-1]])
            sage: P = A*(A.transpose()*A).inverse()*A.transpose()
            sage: T = M.twisted_invariant_module(G, [2,0,-1], action_on_basis=action)
            sage: P == T.projection_matrix()
            True

        Moreover, since there is no component of the sign
        representation in this representation, the projection
        matrix is just the zero matrix::

            sage: T = M.twisted_invariant_module(G, [1,-1,1], action_on_basis=action)
            sage: T.projection_matrix()
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        return self._projection_matrix

    class Element(SubmoduleWithBasis.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: G = SymmetricGroup(3)
                sage: R = G.regular_representation(QQ)
                sage: T = R.twisted_invariant_module([2,0,-1])
                sage: t = T.an_element()
                sage: 5 * t
                10*B[0] + 10*B[1] + 15*B[2]
                sage: t * -2/3
                -4/3*B[0] - 4/3*B[1] - 2*B[2]
                sage: [g * t for g in G]
                [2*B[0] + 2*B[1] + 3*B[2],
                 -4*B[0] + 2*B[1] - 3*B[3],
                 2*B[0] - 4*B[1] - 3*B[2] + 3*B[3],
                 3*B[0] + 2*B[2] + 2*B[3],
                 -3*B[1] - 4*B[2] + 2*B[3],
                 -3*B[0] + 3*B[1] + 2*B[2] - 4*B[3]]
            """
            P = self.parent()
            if scalar in P._group and self_on_left == (P._side == 'right'):
                return P.retract(scalar * P.lift(self))
            return super()._acted_upon_(scalar, self_on_left)
