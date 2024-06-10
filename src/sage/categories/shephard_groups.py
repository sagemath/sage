# sage_setup: distribution = sagemath-categories
r"""
Shephard Groups
"""
#*****************************************************************************
#  Copyright (C) 2016 Travis Scrimshaw <tscrim at ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.generalized_coxeter_groups import GeneralizedCoxeterGroups
from sage.misc.cachefunc import cached_method


class ShephardGroups(Category_singleton):
    r"""
    The category of Shephard groups.

    EXAMPLES::

        sage: from sage.categories.shephard_groups import ShephardGroups
        sage: C = ShephardGroups(); C
        Category of shephard groups

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.shephard_groups import ShephardGroups
            sage: ShephardGroups().super_categories()
            [Category of finite generalized Coxeter groups]
        """
        return [GeneralizedCoxeterGroups().Finite()]
