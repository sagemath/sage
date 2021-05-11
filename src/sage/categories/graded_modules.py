r"""
Graded modules
"""
# ****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory


class GradedModulesCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = GradedAlgebras(QQ)
            sage: C
            Category of graded algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of filtered algebras over Rational Field,
             Category of graded vector spaces over Rational Field]

            sage: AlgebrasWithBasis(QQ).Graded().base_ring()
            Rational Field
            sage: GradedHopfAlgebrasWithBasis(QQ).base_ring()
            Rational Field

        TESTS::

            sage: GradedModules(ZZ)
            Category of graded modules over Integer Ring
            sage: Modules(ZZ).Graded()
            Category of graded modules over Integer Ring
            sage: GradedModules(ZZ) is Modules(ZZ).Graded()
            True
        """
        super(GradedModulesCategory, self).__init__(base_category, base_category.base_ring())

    _functor_category = "Graded"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Graded()  # indirect doctest
            Category of graded algebras with basis over Rational Field
        """
        return "graded {}".format(self.base_category()._repr_object_names())

    @classmethod
    def default_super_categories(cls, category, *args):
        r"""
        Return the default super categories of ``category.Graded()``.

        Mathematical meaning: every graded object (module, algebra,
        etc.) is a filtered object with the (implicit) filtration
        defined by `F_i = \bigoplus_{j \leq i} G_j`.

        INPUT:

        - ``cls`` -- the class ``GradedModulesCategory``
        - ``category`` -- a category

        OUTPUT: a (join) category

        In practice, this returns ``category.Filtered()``, joined
        together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories() <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category.Filtered()`` and ``cat`` for
        each ``cat`` in the super categories of ``category``).

        EXAMPLES:

        Consider ``category=Algebras()``, which has ``cat=Modules()``
        as super category. Then, a grading of an algebra `G`
        is also a filtration of `G`::

            sage: Algebras(QQ).Graded().super_categories()
            [Category of filtered algebras over Rational Field,
             Category of graded vector spaces over Rational Field]

        This resulted from the following call::

            sage: sage.categories.graded_modules.GradedModulesCategory.default_super_categories(Algebras(QQ))
            Join of Category of filtered algebras over Rational Field
             and Category of graded vector spaces over Rational Field
        """
        cat = super(GradedModulesCategory, cls).default_super_categories(category, *args)
        return Category.join([category.Filtered(), cat])

class GradedModules(GradedModulesCategory):
    r"""
    The category of graded modules.

    We consider every graded module `M = \bigoplus_i M_i` as a
    filtered module under the (natural) filtration given by

    .. MATH::

        F_i = \bigoplus_{j < i} M_j.

    EXAMPLES::

        sage: GradedModules(ZZ)
        Category of graded modules over Integer Ring
        sage: GradedModules(ZZ).super_categories()
        [Category of filtered modules over Integer Ring]

    The category of graded modules defines the graded structure which
    shall be preserved by morphisms::

        sage: Modules(ZZ).Graded().additional_structure()
        Category of graded modules over Integer Ring

    TESTS::

        sage: TestSuite(GradedModules(ZZ)).run()
    """
    class ParentMethods:
        pass

    class ElementMethods:
        pass

