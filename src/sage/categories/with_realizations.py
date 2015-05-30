"""
With Realizations Covariant Functorial Construction

.. SEEALSO::

    - :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`
      for an introduction to *realizations* and *with realizations*.
    - :mod:`sage.categories.covariant_functorial_construction`
      for an introduction to covariant functorial constructions.
"""
#*****************************************************************************
#  Copyright (C) 2010-2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category import Category
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory

def WithRealizations(self):
    """
    Returns the category of parents in ``self`` endowed with multiple realizations

    INPUT:

    - ``self`` -- a category

    .. SEEALSO::

        - the documentation and code
          (:mod:`sage.categories.examples.with_realizations`) of
          ``Sets().WithRealizations().example()`` for more on how to use and
          implement a parent with several realizations.

        - :mod:`sage.categories.realizations`


    .. NOTE:: this *function* is actually inserted as a *method* in the class
       :class:`~sage.categories.category.Category` (see
       :meth:`~sage.categories.category.Category.WithRealizations`). It is defined
       here for code locality reasons.

    EXAMPLES::

        sage: Sets().WithRealizations()
        Category of sets with realizations

    .. RUBRIC:: Parent with realizations

    Let us now explain the concept of realizations. A *parent with
    realizations* is a facade parent (see :class:`Sets.Facade`)
    admitting multiple concrete realizations where its elements are
    represented. Consider for example an algebra `A` which admits
    several natural bases::

        sage: A = Sets().WithRealizations().example(); A
        The subset algebra of {1, 2, 3} over Rational Field

    For each such basis `B` one implements a parent `P_B` which
    realizes `A` with its elements represented by expanding them on
    the basis `B`::

        sage: A.F()
        The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
        sage: A.Out()
        The subset algebra of {1, 2, 3} over Rational Field in the Out basis
        sage: A.In()
        The subset algebra of {1, 2, 3} over Rational Field in the In basis

        sage: A.an_element()
        F[{}] + 2*F[{1}] + 3*F[{2}] + F[{1, 2}]

    If `B` and `B'` are two bases, then the change of basis from `B`
    to `B'` is implemented by a canonical coercion between `P_B` and
    `P_{B'}`::

        sage: F = A.F(); In = A.In(); Out = A.Out()
        sage: i = In.an_element(); i
        In[{}] + 2*In[{1}] + 3*In[{2}] + In[{1, 2}]
        sage: F(i)
        7*F[{}] + 3*F[{1}] + 4*F[{2}] + F[{1, 2}]
        sage: F.coerce_map_from(Out)
        Generic morphism:
          From: The subset algebra of {1, 2, 3} over Rational Field in the Out basis
          To:   The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis

    allowing for mixed arithmetic::

        sage: (1 + Out.from_set(1)) * In.from_set(2,3)
        Out[{}] + 2*Out[{1}] + 2*Out[{2}] + 2*Out[{3}] + 2*Out[{1, 2}] + 2*Out[{1, 3}] + 4*Out[{2, 3}] + 4*Out[{1, 2, 3}]

    In our example, there are three realizations::

        sage: A.realizations()
        [The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis,
         The subset algebra of {1, 2, 3} over Rational Field in the In basis,
         The subset algebra of {1, 2, 3} over Rational Field in the Out basis]

    The set of all realizations of `A`, together with the coercion morphisms
    is a category (whose class inherits from
    :class:`~sage.categories.realizations.Category_realization_of_parent`)::

        sage: A.Realizations()
        Category of realizations of The subset algebra of {1, 2, 3} over Rational Field

    The various parent realizing `A` belong to this category::

        sage: A.F() in A.Realizations()
        True

    `A` itself is in the category of algebras with realizations::

        sage: A in Algebras(QQ).WithRealizations()
        True

    The (mostly technical) ``WithRealizations`` categories are the
    analogs of the ``*WithSeveralBases`` categories in
    MuPAD-Combinat. They provide support tools for handling the
    different realizations and the morphisms between them.

    Typically, ``FiniteDimensionalVectorSpaces(QQ).WithRealizations()``
    will eventually be in charge, whenever a coercion `\phi: A\mapsto B` is
    registered, to register `\phi^{-1}` as coercion `B \mapsto A`
    if there is none defined yet. To achieve this,
    ``FiniteDimensionalVectorSpaces`` would provide a nested class
    ``WithRealizations`` implementing the appropriate logic.

    ``WithRealizations`` is a :mod:`regressive covariant functorial
    construction <sage.categories.covariant_functorial_construction>`.
    On our example, this simply means that `A` is automatically in the
    category of rings with realizations (covariance)::

        sage: A in Rings().WithRealizations()
        True

    and in the category of algebras (regressiveness)::

        sage: A in Algebras(QQ)
        True

    .. NOTE::

        For ``C`` a category, ``C.WithRealizations()`` in fact calls
        ``sage.categories.with_realizations.Realizations(C)``. The
        later is responsible for building the hierarchy of the
        categories with realizations in parallel to that of their base
        categories, optimizing away those categories that do not
        provide a ``WithRealizations`` nested class. See
        :mod:`sage.categories.covariant_functorial_construction` for
        the technical details.

    .. NOTE::

        Design question: currently ``WithRealizations`` is a
        regressive construction. That is ``self.WithRealizations()``
        is a subcategory of ``self`` by default::

            sage: Algebras(QQ).WithRealizations().super_categories()
            [Category of algebras over Rational Field,
             Category of monoids with realizations,
             Category of additive unital additive magmas with realizations]

        Is this always desirable? For example,
        ``AlgebrasWithBasis(QQ).WithRealizations()`` should certainly
        be a subcategory of ``Algebras(QQ)``, but not of
        ``AlgebrasWithBasis(QQ)``. This is because
        ``AlgebrasWithBasis(QQ)`` is specifying something about the
        concrete realization.

    TESTS::

        sage: Semigroups().WithRealizations()
        Join of Category of semigroups and Category of sets with realizations
        sage: C = GradedHopfAlgebrasWithBasis(QQ).WithRealizations(); C
        Category of graded hopf algebras with basis over Rational Field with realizations
        sage: C.super_categories()
        [Join of Category of hopf algebras over Rational Field
             and Category of graded algebras over Rational Field]
        sage: TestSuite(Semigroups().WithRealizations()).run()

    """
    return WithRealizationsCategory.category_of(self)

Category.WithRealizations = WithRealizations

class WithRealizationsCategory(RegressiveCovariantConstructionCategory):
    """
    An abstract base class for all categories of parents with multiple
    realizations.

    .. SEEALSO:: :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`

    The role of this base class is to implement some technical goodies, such
    as the name for that category.
    """

    _functor_category = "WithRealizations"

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: C = GradedHopfAlgebrasWithBasis(QQ).WithRealizations(); C #indirect doctest
            Category of graded hopf algebras with basis over Rational Field with realizations
        """
        s = repr(self.base_category())
        return s+" with realizations"


