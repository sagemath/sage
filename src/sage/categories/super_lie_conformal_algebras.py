# sage_setup: distribution = sagemath-categories
"""
Super Lie Conformal Algebras

AUTHORS:

- Reimundo Heluani (2019-10-05): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.misc.abstract_method import abstract_method
from sage.categories.lambda_bracket_algebras import LambdaBracketAlgebras


class SuperLieConformalAlgebras(SuperModulesCategory):
    r"""
    The category of super Lie conformal algebras.

    EXAMPLES::

        sage: LieConformalAlgebras(AA).Super()                                          # needs sage.rings.number_field
        Category of super Lie conformal algebras over Algebraic Real Field

    Notice that we can force to have a *purely even* super Lie
    conformal algebra::

        sage: bosondict = {('a','a'): {1:{('K',0):1}}}
        sage: R = LieConformalAlgebra(QQ, bosondict, names=('a',),                      # needs sage.combinat sage.modules
        ....:                         central_elements=('K',), super=True)
        sage: [g.is_even_odd() for g in R.gens()]                                       # needs sage.combinat sage.modules
        [0, 0]
    """
    def extra_super_categories(self):
        """
        The extra super categories of ``self``.

        EXAMPLES::

            sage: LieConformalAlgebras(QQ).Super().super_categories()
            [Category of super modules over Rational Field,
             Category of Lambda bracket algebras over Rational Field]
        """
        return [LambdaBracketAlgebras(self.base_ring())]

    def example(self):
        """
        An example parent in this category.

        EXAMPLES::

            sage: LieConformalAlgebras(QQ).Super().example()                            # needs sage.combinat sage.modules
            The Neveu-Schwarz super Lie conformal algebra over Rational Field
        """
        from sage.algebras.lie_conformal_algebras.neveu_schwarz_lie_conformal_algebra\
                                      import NeveuSchwarzLieConformalAlgebra
        return NeveuSchwarzLieConformalAlgebra(self.base_ring())

    class ParentMethods:

        def _test_jacobi(self, **options):
            """
            Test the Jacobi axiom of this super Lie conformal algebra.

            INPUT:

            - ``options`` -- any keyword arguments acceptde by :meth:`_tester`

            EXAMPLES:

            By default, this method tests only the elements returned by
            ``self.some_elements()``::

                sage: V = lie_conformal_algebras.Affine(QQ, 'B2')                       # needs sage.combinat sage.modules
                sage: V._test_jacobi()          # long time (6 seconds)                 # needs sage.combinat sage.modules

            It works for super Lie conformal algebras too::

                sage: V = lie_conformal_algebras.NeveuSchwarz(QQ)                       # needs sage.combinat sage.modules
                sage: V._test_jacobi()                                                  # needs sage.combinat sage.modules

            We can use specific elements by passing the ``elements``
            keyword argument::

                sage: V = lie_conformal_algebras.Affine(QQ, 'A1',                       # needs sage.combinat sage.modules
                ....:                                   names=('e', 'h', 'f'))
                sage: V.inject_variables()                                              # needs sage.combinat sage.modules
                Defining e, h, f, K
                sage: V._test_jacobi(elements=(e, 2*f + h, 3*h))                        # needs sage.combinat sage.modules

            TESTS::

                sage: wrongdict = {('a', 'a'): {0: {('b', 0): 1}},
                ....:              ('b', 'a'): {0: {('a', 0): 1}}}
                sage: V = LieConformalAlgebra(QQ, wrongdict,                            # needs sage.combinat sage.modules
                ....:                         names=('a', 'b'), parity=(1, 0))
                sage: V._test_jacobi()                                                  # needs sage.combinat sage.modules
                Traceback (most recent call last):
                ...
                AssertionError: {(0, 0): -3*a} != {} - {(0, 0): -3*a} + {}
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            # Try our best to avoid non-homogeneous elements
            elements = []
            for s in S:
                try:
                    s.is_even_odd()
                except ValueError:
                    try:
                        elements.extend([s.even_component(), s.odd_component()])
                    except (AttributeError, ValueError):
                        pass
                    continue
                elements.append(s)
            S = elements
            from sage.misc.misc import some_tuples
            from sage.arith.misc import binomial
            pz = tester._instance.zero()
            for x,y,z in some_tuples(S, 3, tester._max_runs):
                if x.is_zero() or y.is_zero():
                    sgn = 1
                elif x.is_even_odd() * y.is_even_odd():
                    sgn = -1
                else:
                    sgn = 1
                brxy = x.bracket(y)
                brxz = x.bracket(z)
                bryz = y.bracket(z)
                br1 = {k: x.bracket(v) for k,v in bryz.items()}
                br2 = {k: v.bracket(z) for k,v in brxy.items()}
                br3 = {k: y.bracket(v) for k,v in brxz.items()}
                jac1 = {(j,k): v for k in br1 for j,v in br1[k].items()}
                jac3 = {(k,j): v for k in br3 for j,v in br3[k].items()}
                jac2 = {}
                for k,br in br2.items():
                    for j,v in br.items():
                        for r in range(j+1):
                            jac2[(k+r, j-r)] = (jac2.get((k+r, j-r), pz)
                                                + binomial(k+r, r)*v)
                for k,v in jac2.items():
                    jac1[k] = jac1.get(k, pz) - v
                for k,v in jac3.items():
                    jac1[k] = jac1.get(k, pz) - sgn*v
                jacobiator = {k: v for k,v in jac1.items() if v}
                tester.assertDictEqual(jacobiator, {})

    class ElementMethods:

        @abstract_method
        def is_even_odd(self):
            """
            Return ``0`` if this element is *even* and ``1`` if it is
            *odd*.

            EXAMPLES::

                sage: R = lie_conformal_algebras.NeveuSchwarz(QQ)                       # needs sage.combinat sage.modules
                sage: R.inject_variables()                                              # needs sage.combinat sage.modules
                Defining L, G, C
                sage: G.is_even_odd()                                                   # needs sage.combinat sage.modules
                1
            """

    class Graded(GradedModulesCategory):
        """
        The category of H-graded super Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(AA).Super().Graded()                             # needs sage.rings.number_field
            Category of H-graded super Lie conformal algebras over Algebraic Real Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of this category.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).Graded()                              # needs sage.rings.number_field
                Category of H-graded Lie conformal algebras over Algebraic Field
            """
            return "H-graded {}".format(self.base_category()._repr_object_names())
