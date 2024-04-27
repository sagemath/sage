r"""
Objects
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.homsets import HomsetsCategory

#############################################################
# Generic category (default when requesting category of
# an object using misc.functional.category
#############################################################

class Objects(Category_singleton):
    """
    The category of all objects
    the basic category

    EXAMPLES::

        sage: Objects()
        Category of objects
        sage: Objects().super_categories()
        []

    TESTS::

        sage: TestSuite(Objects()).run()
    """

    def additional_structure(self):
        """
        Return ``None``

        Indeed, by convention, the category of objects defines no
        additional structure.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: Objects().additional_structure()
        """
        return None

    def super_categories(self):
        """
        EXAMPLES::

            sage: Objects().super_categories()
            []
        """
        return []

    def __contains__(self, x):
        """
        Anything is in the category of objects.

        EXAMPLES::

            sage: int(1) in Objects()
            True
            sage: ZZ     in Objects()
            True
            sage: 2/3    in Objects()
            True
        """
        return True

    class SubcategoryMethods:
        @cached_method
        def Homsets(self):
            r"""
            Return the category of homsets between objects of this category.

            EXAMPLES::

                sage: Sets().Homsets()
                Category of homsets of sets

                sage: Rings().Homsets()
                Category of homsets of unital magmas and additive unital additive magmas

            .. NOTE:: Background

                Information, code, documentation, and tests about the
                category of homsets of a category ``Cs`` should go in
                the nested class ``Cs.Homsets``. They will then be
                made available to homsets of any subcategory of
                ``Cs``.

                Assume, for example, that homsets of ``Cs`` are ``Cs``
                themselves. This information can be implemented in the
                method ``Cs.Homsets.extra_super_categories`` to make
                ``Cs.Homsets()`` a subcategory of ``Cs()``.

                Methods about the homsets themselves should go in the
                nested class ``Cs.Homsets.ParentMethods``.

                Methods about the morphisms can go in the nested class
                ``Cs.Homsets.ElementMethods``. However it's generally
                preferable to put them in the nested class
                ``Cs.MorphimMethods``; indeed they will then apply to
                morphisms of all subcategories of ``Cs``, and not only
                full subcategories.


            .. SEEALSO::

                :class:`~.covariant_functorial_construction.FunctorialConstruction`

            .. TODO::

                - Design a mechanism to specify that an axiom is
                  compatible with taking subsets. Examples:
                  ``Finite``, ``Associative``, ``Commutative`` (when
                  meaningful), but not ``Infinite`` nor ``Unital``.

                - Design a mechanism to specify that, when `B` is a
                  subcategory of `A`, a `B`-homset is a subset of the
                  corresponding `A` homset. And use it to recover all
                  the relevant axioms from homsets in super categories.

                - For instances of redundant code due to this missing
                  feature, see:

                  - :meth:`AdditiveMonoids.Homsets.extra_super_categories`
                  - :meth:`HomsetsCategory.extra_super_categories`
                    (slightly different nature)
                  - plus plenty of spots where this is not implemented.
            """
            return HomsetsCategory.category_of(self)

        @cached_method
        def Endsets(self):
            r"""
            Return the category of endsets between objects of this category.

            EXAMPLES::

                sage: Sets().Endsets()
                Category of endsets of sets

                sage: Rings().Endsets()
                Category of endsets of unital magmas and additive unital additive magmas

            .. SEEALSO::

                - :meth:`Homsets`
            """
            return self.Homsets()._with_axiom("Endset")

    class ParentMethods:
        """
        Methods for all category objects
        """
