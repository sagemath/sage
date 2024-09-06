# sage_setup: distribution = sagemath-categories
r"""
Lie Groups
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

#from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.groups import Groups
from sage.categories.manifolds import Manifolds


class LieGroups(Category_over_base_ring):
    r"""
    The category of Lie groups.

    A Lie group is a topological group with a smooth manifold structure.

    EXAMPLES::

        sage: from sage.categories.lie_groups import LieGroups
        sage: C = LieGroups(QQ); C
        Category of Lie groups over Rational Field

    TESTS::

        sage: TestSuite(C).run(skip='_test_category_over_bases')
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups(QQ).super_categories()
            [Category of topological groups,
             Category of smooth manifolds over Rational Field]
        """
        return [Groups().Topological(), Manifolds(self.base()).Smooth()]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of Lie groups defines no new
        structure: a morphism of topological spaces and of smooth
        manifolds is a morphism as Lie groups.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups(QQ).additional_structure()
        """
        return None

    # Because Lie is a name that deserves to be capitalized
    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups(QQ) # indirect doctest
            Category of Lie groups over Rational Field
        """
        return "Lie groups over {}".format(self.base_ring())
