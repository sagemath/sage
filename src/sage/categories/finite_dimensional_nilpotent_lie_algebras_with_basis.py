# sage_setup: distribution = sagemath-categories
r"""
Finite Dimensional Nilpotent Lie Algebras With Basis

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras


class FiniteDimensionalNilpotentLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    r"""
    Category of finite dimensional nilpotent Lie algebras with basis.

    TESTS::

        sage: C1 = LieAlgebras(QQ).FiniteDimensional().WithBasis().Nilpotent()
        sage: C2 = LieAlgebras(QQ).FiniteDimensional().Nilpotent().WithBasis()
        sage: C3 = LieAlgebras(QQ).Nilpotent().FiniteDimensional().WithBasis()
        sage: C4 = LieAlgebras(QQ).Nilpotent().WithBasis().FiniteDimensional()
        sage: C5 = LieAlgebras(QQ).WithBasis().Nilpotent().FiniteDimensional()
        sage: C6 = LieAlgebras(QQ).WithBasis().FiniteDimensional().Nilpotent()
        sage: C1 is C2
        True
        sage: C2 is C3
        True
        sage: C3 is C4
        True
        sage: C4 is C5
        True
        sage: C5 is C6
        True
        sage: TestSuite(C1).run()
    """
    _base_category_class_and_axiom = (LieAlgebras.FiniteDimensional.WithBasis, "Nilpotent")

    class ParentMethods:

        def _test_nilpotency(self, **options):
            r"""
            Test that ``self`` is nilpotent and has the correct step.

            INPUT:

            - ``options`` -- any keyword arguments accepted by
              :meth:`_tester`

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
                sage: L._test_nilpotency()
                sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}},
                ....:                nilpotent=True, step=3)
                sage: L._test_nilpotency()
                Traceback (most recent call last):
                ...
                AssertionError: claimed nilpotency step 3
                does not match the actual nilpotency step 2
                sage: L = LieAlgebra(QQ, {('X','Y'): {'X': 1}}, nilpotent=True)
                sage: L._test_nilpotency()
                Traceback (most recent call last):
                ...
                AssertionError: final term of lower central series is nonzero

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)

            lcs = self.lower_central_series(submodule=True)
            tester.assertEqual(lcs[-1].dimension(), 0,
                msg="final term of lower central series is nonzero")

            step = self.step()
            tester.assertEqual(len(lcs) - 1, step,
                msg="claimed nilpotency step %d does not match the "
                "actual nilpotency step %d" % (step, len(lcs) - 1))

        def lie_group(self, name='G', **kwds):
            r"""
            Return the Lie group associated to ``self``.

            INPUT:

            - ``name`` -- string (default: ``'G'``);
              the name (symbol) given to the Lie group

            EXAMPLES:

            We define the Heisenberg group::

                sage: L = lie_algebras.Heisenberg(QQ, 1)                                # needs sage.combinat sage.modules
                sage: G = L.lie_group('G'); G                                           # needs sage.combinat sage.modules sage.symbolic
                Lie group G of Heisenberg algebra of rank 1 over Rational Field

            We test multiplying elements of the group::

                sage: # needs sage.combinat sage.modules sage.symbolic
                sage: p, q, z = L.basis()
                sage: g = G.exp(p); g
                exp(p1)
                sage: h = G.exp(q); h
                exp(q1)
                sage: g * h
                exp(p1 + q1 + 1/2*z)

            We extend an element of the Lie algebra to a left-invariant
            vector field::

                sage: X = G.left_invariant_extension(2*p + 3*q, name='X'); X            # needs sage.combinat sage.modules sage.symbolic
                Vector field X on the Lie group G of
                 Heisenberg algebra of rank 1 over Rational Field
                sage: X.at(G.one()).display()                                           # needs sage.combinat sage.modules sage.symbolic
                X = 2 ∂/∂x_0 + 3 ∂/∂x_1
                sage: X.display()                                                       # needs sage.combinat sage.modules sage.symbolic
                X = 2 ∂/∂x_0 + 3 ∂/∂x_1 + (3/2*x_0 - x_1) ∂/∂x_2

            .. SEEALSO::

                :class:`~sage.groups.lie_gps.nilpotent_lie_group.NilpotentLieGroup`
            """
            from sage.groups.lie_gps.nilpotent_lie_group import NilpotentLieGroup
            return NilpotentLieGroup(self, name, **kwds)

        def step(self):
            r"""
            Return the nilpotency step of ``self``.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
                sage: L.step()
                2
                sage: sc = {('X','Y'): {'Z': 1}, ('X','Z'): {'W': 1}}
                sage: LieAlgebra(QQ, sc, nilpotent=True).step()
                3
            """
            if not hasattr(self, '_step'):
                self._step = len(self.lower_central_series(submodule=True)) - 1
            return self._step

        def is_nilpotent(self):
            r"""
            Return ``True`` since ``self`` is nilpotent.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, {('x','y'): {'z': 1}}, nilpotent=True)         # needs sage.combinat sage.modules
                sage: L.is_nilpotent()                                                  # needs sage.combinat sage.modules
                True
            """
            return True
