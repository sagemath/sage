# sage_setup: distribution = sagemath-categories
r"""
Filtered Algebras
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.categories.filtered_modules import FilteredModulesCategory


class FilteredAlgebras(FilteredModulesCategory):
    r"""
    The category of filtered algebras.

    An algebra `A` over a commutative ring `R` is *filtered* if
    `A` is endowed with a structure of a filtered `R`-module
    (whose underlying `R`-module structure is identical with
    that of the `R`-algebra `A`) such that the indexing set `I`
    (typically `I = \NN`) is also an additive abelian monoid,
    the unity `1` of `A` belongs to `F_0`, and we have
    `F_i \cdot F_j \subseteq F_{i+j}` for all `i, j \in I`.

    EXAMPLES::

        sage: Algebras(ZZ).Filtered()
        Category of filtered algebras over Integer Ring
        sage: Algebras(ZZ).Filtered().super_categories()
        [Category of algebras over Integer Ring,
         Category of filtered modules over Integer Ring]

    TESTS::

        sage: TestSuite(Algebras(ZZ).Filtered()).run()

    REFERENCES:

    - :wikipedia:`Filtered_algebra`
    """
    class ParentMethods:
        @abstract_method(optional=True)
        def graded_algebra(self):
            """
            Return the associated graded algebra to ``self``.

            .. TODO::

                Implement a version of the associated graded algebra
                which does not require ``self`` to have a
                distinguished basis.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: A.graded_algebra()
                Graded Algebra of An example of a filtered algebra with basis:
                 the universal enveloping algebra of
                 Lie algebra of RR^3 with cross product over Integer Ring
            """
