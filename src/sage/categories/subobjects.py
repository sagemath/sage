"""
Subobjects Functorial Construction

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

class SubobjectsCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "Subobjects"

    @classmethod
    def default_super_categories(cls, category):
        """
        Returns the default super categories of ``category.Subobjects()``

        Mathematical meaning: if `A` is a subobject of `B` in the
        category `C`, then `A` is also a subquotient of `B` in the
        category `C`.

        INPUT:

         - ``cls`` -- the class ``SubobjectsCategory``
         - ``category`` -- a category `Cat`

        OUTPUT: a (join) category

        In practice, this returns ``category.Subquotients()``, joined
        together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories() <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category`` and ``cat.Subobjects()`` for
        each ``cat`` in the super categories of ``category``).

        EXAMPLES:

        Consider ``category=Groups()``, which has ``cat=Monoids()`` as
        super category. Then, a subgroup of a group `G` is
        simultaneously a subquotient of `G`, a group by itself, and a
        submonoid of `G`::

            sage: Groups().Subobjects().super_categories()
            [Category of groups, Category of subquotients of monoids, Category of subobjects of sets]

        Mind the last item above: there is indeed currently nothing
        implemented about submonoids.

        This resulted from the following call::

            sage: sage.categories.subobjects.SubobjectsCategory.default_super_categories(Groups())
            Join of Category of groups and Category of subquotients of monoids and Category of subobjects of sets
        """
        return Category.join([category.Subquotients(), super(SubobjectsCategory, cls).default_super_categories(category)])
