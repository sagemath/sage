r"""
Finitely enumerated families
"""

#*****************************************************************************
#       Copyright (C) 2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.enumerated_families import EnumeratedFamilies


class FinitelyEnumeratedFamilies(Category_singleton):
    r"""
    The category of finitely enumerated families

    An *finitely enumerated family* is an enumerated family whose keys
    form a finite enumerated set; this implies that the family is a finite set.

    EXAMPLES::

        sage: from sage.categories.enumerated_families import EnumeratedFamilies
        sage: from sage.categories.finitely_enumerated_families import FinitelyEnumeratedFamilies
        sage: FinitelyEnumeratedFamilies()
        Category of finitely enumerated families

        sage: FinitelyEnumeratedFamilies() == EnumeratedFamilies() & FiniteSets()
        False
        sage: F = Family(ZZ, lambda x: x % 7, is_injective=False, category=FiniteSets())
        sage: F in EnumeratedFamilies() & FiniteSets()
        True
        sage: F in FinitelyEnumeratedFamilies()
        False
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.finitely_enumerated_families import FinitelyEnumeratedFamilies
            sage: FinitelyEnumeratedFamilies().super_categories()
            [Category of enumerated families, Category of finite enumerated sets]
        """
        return [EnumeratedFamilies(), FiniteEnumeratedSets()]
