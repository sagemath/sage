# sage_setup: distribution = sagemath-objects
"""
Realizations Covariant Functorial Construction

.. SEEALSO::

    - :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`
      for an introduction to *realizations* and *with realizations*.
    - :mod:`sage.categories.covariant_functorial_construction`
      for an introduction to covariant functorial constructions.
    - :mod:`sage.categories.examples.with_realizations` for an example.
"""
#*****************************************************************************
#  Copyright (C) 2010-2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.bindable_class import BindableClass
from sage.categories.category import Category
from sage.categories.category_types import Category_over_base
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory


class RealizationsCategory(RegressiveCovariantConstructionCategory):
    """
    An abstract base class for all categories of realizations category.

    Realization are implemented as
    :class:`~sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory`.
    See there for the documentation of how the various bindings such
    as ``Sets().Realizations()`` and ``P.Realizations()``, where ``P``
    is a parent, work.

    .. SEEALSO:: :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`

    TESTS::

        sage: Sets().Realizations
        <bound method Realizations of Category of sets>
        sage: Sets().Realizations()
        Category of realizations of sets
        sage: Sets().Realizations().super_categories()
        [Category of sets]
        sage: Groups().Realizations().super_categories()
        [Category of groups, Category of realizations of unital magmas]
    """

    _functor_category = "Realizations"


def Realizations(self):
    """
    Return the category of realizations of the parent ``self`` or of objects
    of the category ``self``

    INPUT:

    - ``self`` -- a parent or a concrete category

    .. NOTE:: this *function* is actually inserted as a *method* in the class
       :class:`~sage.categories.category.Category` (see
       :meth:`~sage.categories.category.Category.Realizations`). It is defined
       here for code locality reasons.

    EXAMPLES:

    The category of realizations of some algebra::

        sage: Algebras(QQ).Realizations()
        Join of Category of algebras over Rational Field
            and Category of realizations of unital magmas

    The category of realizations of a given algebra::

        sage: A = Sets().WithRealizations().example(); A                                # needs sage.modules
        The subset algebra of {1, 2, 3} over Rational Field
        sage: A.Realizations()                                                          # needs sage.modules
        Category of realizations of
         The subset algebra of {1, 2, 3} over Rational Field

        sage: C = GradedHopfAlgebrasWithBasis(QQ).Realizations(); C
        Join of Category of graded Hopf algebras with basis over Rational Field
            and Category of realizations of Hopf algebras over Rational Field
        sage: C.super_categories()
        [Category of graded Hopf algebras with basis over Rational Field,
         Category of realizations of Hopf algebras over Rational Field]

        sage: TestSuite(C).run()

    .. SEEALSO::

        - :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`
        - :class:`ClasscallMetaclass`

    .. TODO::

        Add an optional argument to allow for::

            sage: Realizations(A, category=Blahs()) # todo: not implemented
    """
    if isinstance(self, Category):
        return RealizationsCategory.category_of(self)
    else:
        return getattr(self.__class__, "Realizations")(self)


Category.Realizations = Realizations


class Category_realization_of_parent(Category_over_base, BindableClass):
    """
    An abstract base class for categories of all realizations of a given parent.

    INPUT:

    - ``parent_with_realization`` -- a parent

    .. SEEALSO:: :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`

    EXAMPLES::

        sage: A = Sets().WithRealizations().example(); A                                # needs sage.modules
        The subset algebra of {1, 2, 3} over Rational Field

    The role of this base class is to implement some technical goodies, like
    the binding ``A.Realizations()`` when a subclass ``Realizations`` is
    implemented as a nested class in ``A``
    (see the :mod:`code of the example <sage.categories.examples.with_realizations.SubsetAlgebra>`)::

        sage: C = A.Realizations(); C                                                   # needs sage.modules
        Category of realizations of
         The subset algebra of {1, 2, 3} over Rational Field

    as well as the name for that category.
    """
    def __init__(self, parent_with_realization):
        """
        TESTS::

            sage: # needs sage.combinat sage.modules
            sage: from sage.categories.realizations import Category_realization_of_parent
            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: C = A.Realizations(); C
            Category of realizations of
             The subset algebra of {1, 2, 3} over Rational Field
            sage: isinstance(C, Category_realization_of_parent)
            True
            sage: C.parent_with_realization
            The subset algebra of {1, 2, 3} over Rational Field
            sage: TestSuite(C).run(skip=["_test_category_over_bases"])

        .. TODO::

            Fix the failing test by making ``C`` a singleton
            category. This will require some fiddling with the
            assertion in :meth:`Category_singleton.__classcall__`
        """
        Category_over_base.__init__(self, parent_with_realization)
        self.parent_with_realization = parent_with_realization

    def _get_name(self):
        """
        Return a human readable string specifying which kind of bases this category is for.

        It is obtained by splitting and lower casing the last part of
        the class name.

        EXAMPLES::

            sage: from sage.categories.realizations import Category_realization_of_parent
            sage: class MultiplicativeBasesOnPrimitiveElements(Category_realization_of_parent):
            ....:     def super_categories(self): return [Objects()]
            sage: Sym = SymmetricFunctions(QQ); Sym.rename('Sym')                       # needs sage.combinat sage.modules
            sage: MultiplicativeBasesOnPrimitiveElements(Sym)._get_name()               # needs sage.combinat sage.modules
            'multiplicative bases on primitive elements'
        """
        import re
        return re.sub(".[A-Z]", lambda s: s.group()[0]+" "+s.group()[1], self.__class__.__base__.__name__.split(".")[-1]).lower()

    def _repr_object_names(self):
        """
        Return the name of the objects of this category.

        .. SEEALSO:: :meth:`Category._repr_object_names`

        EXAMPLES::

            sage: from sage.categories.realizations import Category_realization_of_parent
            sage: class MultiplicativeBasesOnPrimitiveElements(Category_realization_of_parent):
            ....:     def super_categories(self): return [Objects()]
            sage: Sym = SymmetricFunctions(QQ); Sym.rename('Sym')                       # needs sage.combinat sage.modules
            sage: C = MultiplicativeBasesOnPrimitiveElements(Sym); C                    # needs sage.combinat sage.modules
            Category of multiplicative bases on primitive elements of Sym
            sage: C._repr_object_names()                                                # needs sage.combinat sage.modules
            'multiplicative bases on primitive elements of Sym'
        """
        return "{} of {}".format(self._get_name(), self.base())
