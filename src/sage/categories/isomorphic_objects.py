# sage_setup: distribution = sagemath-objects
"""
Isomorphic Objects Functorial Construction

AUTHORS:

 - Nicolas M. Thiery (2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category import Category
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory


class IsomorphicObjectsCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "IsomorphicObjects"

    @classmethod
    def default_super_categories(cls, category):
        """
        Return the default super categories of ``category.IsomorphicObjects()``.

        Mathematical meaning: if `A` is the image of `B` by an
        isomorphism in the category `C`, then `A` is both a subobject
        of `B` and a quotient of `B` in the category `C`.

        INPUT:

        - ``cls`` -- the class ``IsomorphicObjectsCategory``
        - ``category`` -- a category `Cat`

        OUTPUT: a (join) category

        In practice, this returns ``category.Subobjects()`` and
        ``category.Quotients()``, joined together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories() <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category`` and
        ``cat.IsomorphicObjects()`` for each ``cat`` in the super
        categories of ``category``).

        EXAMPLES:

        Consider ``category=Groups()``, which has ``cat=Monoids()`` as
        super category. Then, the image of a group `G'` by a group
        isomorphism is simultaneously a subgroup of `G`, a subquotient
        of `G`, a group by itself, and the image of `G` by a monoid
        isomorphism::

            sage: Groups().IsomorphicObjects().super_categories()
            [Category of groups,
             Category of subquotients of monoids,
             Category of quotients of semigroups,
             Category of isomorphic objects of sets]

        Mind the last item above: there is indeed currently nothing
        implemented about isomorphic objects of monoids.

        This resulted from the following call::

            sage: sage.categories.isomorphic_objects.IsomorphicObjectsCategory.default_super_categories(Groups())
            Join of Category of groups and
            Category of subquotients of monoids and
            Category of quotients of semigroups and
            Category of isomorphic objects of sets
        """
        return Category.join([category.Subobjects(), category.Quotients(),
                              super().default_super_categories(category)])
