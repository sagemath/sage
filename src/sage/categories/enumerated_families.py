r"""
Enumerated families
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
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.families import Families


class EnumeratedFamilies(Category_singleton):
    r"""
    The category of enumerated families

    An *enumerated family* is a family whose keys are an enumerated set
    that induces an enumeration of the family.

    - ``F[key]``: the value (element) indexed by ``key``.

    The standard methods for a family ``F`` are:

    - ``F.keys()``, ``F.values()``, ``F.items()``: methods similar to those of
      :class:`dict` (or :class:`collections.abc.Mapping`).

    The map is not necessarily injective.

    - ``iter(F)``: returns an iterator for the elements (values) of ``F`` as a set
      in the same order as ``iter(F.values())``, but with duplicates removed.

    EXAMPLES::

        sage: from sage.categories.enumerated_families import EnumeratedFamilies
        sage: from sage.categories.families import Families
        sage: EnumeratedFamilies()
        Category of enumerated families

        sage: EnumeratedFamilies() == Families() & EnumeratedSets()
        False

        sage: Families() & EnumeratedSets()  # BUG
        Category of enumerated families
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.enumerated_families import EnumeratedFamilies
            sage: EnumeratedFamilies().super_categories()
            [Category of families, Category of enumerated sets]
        """
        return [Families(), EnumeratedSets()]
