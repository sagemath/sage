# sage_setup: distribution = sagemath-categories
r"""
Additive groups
"""
#*****************************************************************************
#  Copyright (C) 2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_singleton, CategoryWithAxiom
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.additive_monoids import AdditiveMonoids
from sage.cpython.getattr import raw_getattr
Groups = LazyImport('sage.categories.groups', 'Groups', at_startup=True)


class AdditiveGroups(CategoryWithAxiom_singleton):
    r"""
    The category of additive groups.

    An *additive group* is a set with an internal binary operation `+` which
    is associative, admits a zero, and where every element can be negated.

    EXAMPLES::

        sage: from sage.categories.additive_groups import AdditiveGroups
        sage: from sage.categories.additive_monoids import AdditiveMonoids
        sage: AdditiveGroups()
        Category of additive groups
        sage: AdditiveGroups().super_categories()
        [Category of additive inverse additive unital additive magmas,
         Category of additive monoids]
        sage: AdditiveGroups().all_super_categories()
        [Category of additive groups,
         Category of additive inverse additive unital additive magmas,
         Category of additive monoids,
         Category of additive unital additive magmas,
         Category of additive semigroups,
         Category of additive magmas,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

        sage: AdditiveGroups().axioms()
        frozenset({'AdditiveAssociative', 'AdditiveInverse', 'AdditiveUnital'})
        sage: AdditiveGroups() is AdditiveMonoids().AdditiveInverse()
        True

    TESTS::

        sage: C = AdditiveGroups()
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (AdditiveMonoids, "AdditiveInverse")

    class Algebras(AlgebrasCategory):
        class ParentMethods:
            group = raw_getattr(Groups.Algebras.ParentMethods, "group")

    class Finite(CategoryWithAxiom):
        class Algebras(AlgebrasCategory):
            extra_super_categories = raw_getattr(Groups.Finite.Algebras, "extra_super_categories")

            class ParentMethods:
                __init_extra__ = raw_getattr(Groups.Finite.Algebras.ParentMethods, "__init_extra__")

    AdditiveCommutative = LazyImport('sage.categories.commutative_additive_groups', 'CommutativeAdditiveGroups', at_startup=True)
