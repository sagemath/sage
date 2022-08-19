r"""
Families
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
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_singleton
from sage.categories.sets_cat import Sets


class Families(Category_singleton):
    r"""
    The category of families

    A *family* is a set together with a map from a set of "keys" onto it.

    Morphisms of :class:`Families` preserve this map. This is the additional
    structure compared to :class:`Sets`. Hence, equality of families takes
    this map into account.

    The standard methods for a family ``F`` are:

    - ``F.keys()``, ``F.values()``, ``F.items()``: methods similar to those of
      :class:`dict` (or :class:`collections.abc.Mapping`).

    - ``F[key]``: the value (element) indexed by ``key``.

    A family is not necessarily enumerated, and the map is not necessarily injective.
    However, if it is iterable, then:

    - ``iter(F)``: returns an iterator for the elements (values) of ``F`` as a set,
      i.e., no value appears twice.

    Stronger guarantees are given by the category :class:`EnumeratedFamilies`.

    EXAMPLES::

        sage: from sage.categories.families import Families
        sage: Families().is_full_subcategory(Sets())
        False

    TESTS::

        sage: C = Families()
        sage: TestSuite(C).run()
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.families import Families
            sage: Families().super_categories()
            [Category of sets]
        """
        return [Sets().Facade()]

    class ParentMethods:

        #def __getitem__(self, i):

        @abstract_method
        def keys(self):
            """
            Return the keys of the family.

            This may or may not be an iterable.

            EXAMPLES::

                sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
                sage: sorted(f.keys())
                [3, 4, 7]
            """

        @abstract_method(optional=True)
        def values(self):
            """
            Return the elements (values) of the family.

            If :meth:`keys` returns an iterable, then :meth:`values` will
            return an iterable parallel to that. When the family is not injective,
            values will appear multiple times in the iteration.

            EXAMPLES::

                sage: f = Family(["c", "a", "b"], lambda x: x + x)
                sage: sorted(f.values())
                ['aa', 'bb', 'cc']
            """

        def items(self):
            """
            Return the key-value pairs of the family.

            This may or may not be an iterable.

            A key can only appear once, but if the function is not injective, values will
            appear multiple times.

            EXAMPLES::

                sage: f = Family(["a", "ab", "bc", "def"], len)
                sage: sorted(f.items())
                [('a', 1), ('ab', 2), ('bc', 2), ('def', 3)]
            """
            return zip(self.keys(), self.values())

        @cached_method
        def as_set(self):
            """
            Return the elements (values) of this family as a set.

            EXAMPLES::

                sage: f = Family({1: 'a', 2: 'b', 3: 'c'})
                sage: g = Family({1: 'b', 2: 'c', 3: 'a'})
                sage: f == g
                False
                sage: f.as_set() == g.as_set()
                True

            This is the same as calling :func:`~sage.sets.set.Set` on ``self``::

                sage: f.as_set()
                {...}
                sage: Set(f)
                {...}
            """
            return Set(self.values())
