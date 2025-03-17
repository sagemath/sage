# sage_setup: distribution = sagemath-objects
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
    r"""
    Return the category of parents in ``self`` endowed with multiple realizations.

    INPUT:

    - ``self`` -- a category

    .. SEEALSO::

        - The documentation and code
          (:mod:`sage.categories.examples.with_realizations`) of
          ``Sets().WithRealizations().example()`` for more on how to use and
          implement a parent with several realizations.

        - Various use cases:

          - :class:`SymmetricFunctions`
          - :class:`QuasiSymmetricFunctions`
          - :class:`NonCommutativeSymmetricFunctions`
          - :class:`SymmetricFunctionsNonCommutingVariables`
          - :class:`DescentAlgebra`
          - :class:`algebras.Moebius`
          - :class:`IwahoriHeckeAlgebra`
          - :class:`ExtendedAffineWeylGroup`

        - The `Implementing Algebraic Structures
          <../../../../../thematic_tutorials/tutorial-implementing-algebraic-structures>`_
          thematic tutorial.

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

        sage: A = Sets().WithRealizations().example(); A                                # needs sage.modules
        The subset algebra of {1, 2, 3} over Rational Field

    For each such basis `B` one implements a parent `P_B` which
    realizes `A` with its elements represented by expanding them on
    the basis `B`::

        sage: # needs sage.modules
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

        sage: # needs sage.modules
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

        sage: (1 + Out.from_set(1)) * In.from_set(2,3)                                  # needs sage.modules
        Out[{}] + 2*Out[{1}] + 2*Out[{2}] + 2*Out[{3}] + 2*Out[{1, 2}]
        + 2*Out[{1, 3}] + 4*Out[{2, 3}] + 4*Out[{1, 2, 3}]

    In our example, there are three realizations::

        sage: A.realizations()                                                          # needs sage.modules
        [The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis,
         The subset algebra of {1, 2, 3} over Rational Field in the In basis,
         The subset algebra of {1, 2, 3} over Rational Field in the Out basis]

    Instead of manually defining the shorthands ``F``, ``In``, and
    ``Out``, as above one can just do::

        sage: A.inject_shorthands()                                                     # needs sage.modules
        Defining F as shorthand for
         The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
        Defining In as shorthand for
         The subset algebra of {1, 2, 3} over Rational Field in the In basis
        Defining Out as shorthand for
         The subset algebra of {1, 2, 3} over Rational Field in the Out basis

    .. RUBRIC:: Rationale

    Besides some goodies described below, the role of `A` is threefold:

    - To provide, as illustrated above, a single entry point for the
      algebra as a whole: documentation, access to its properties and
      different realizations, etc.

    - To provide a natural location for the initialization of the
      bases and the coercions between, and other methods that are
      common to all bases.

    - To let other objects refer to `A` while allowing elements to be
      represented in any of the realizations.

    We now illustrate this second point by defining the polynomial
    ring with coefficients in `A`::

        sage: P = A['x']; P                                                             # needs sage.modules
        Univariate Polynomial Ring in x over
         The subset algebra of {1, 2, 3} over Rational Field
        sage: x = P.gen()                                                               # needs sage.modules

    In the following examples, the coefficients turn out to be all
    represented in the `F` basis::

        sage: P.one()                                                                   # needs sage.modules
        F[{}]
        sage: (P.an_element() + 1)^2                                                    # needs sage.modules
        F[{}]*x^2 + 2*F[{}]*x + F[{}]

    However we can create a polynomial with mixed coefficients, and
    compute with it::

        sage: p = P([1, In[{1}], Out[{2}] ]); p                                         # needs sage.modules
        Out[{2}]*x^2 + In[{1}]*x + F[{}]
        sage: p^2                                                                       # needs sage.modules
        Out[{2}]*x^4
        + (-8*In[{}] + 4*In[{1}] + 8*In[{2}] + 4*In[{3}]
           - 4*In[{1, 2}] - 2*In[{1, 3}] - 4*In[{2, 3}] + 2*In[{1, 2, 3}])*x^3
        + (F[{}] + 3*F[{1}] + 2*F[{2}] - 2*F[{1, 2}] - 2*F[{2, 3}] + 2*F[{1, 2, 3}])*x^2
        + (2*F[{}] + 2*F[{1}])*x
        + F[{}]

    Note how each coefficient involves a single basis which need not
    be that of the other coefficients. Which basis is used depends on
    how coercion happened during mixed arithmetic and needs not be
    deterministic.

    One can easily coerce all coefficient to a given basis with::

        sage: p.map_coefficients(In)                                                    # needs sage.modules
        (-4*In[{}] + 2*In[{1}] + 4*In[{2}] + 2*In[{3}]
         - 2*In[{1, 2}] - In[{1, 3}] - 2*In[{2, 3}] + In[{1, 2, 3}])*x^2
        + In[{1}]*x + In[{}]

    Alas, the natural notation for constructing such polynomials does
    not yet work::

        sage: In[{1}] * x                                                               # needs sage.modules
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for *:
        'The subset algebra of {1, 2, 3} over Rational Field in the In basis'
        and 'Univariate Polynomial Ring in x over
        The subset algebra of {1, 2, 3} over Rational Field'

    .. RUBRIC:: The category of realizations of `A`

    The set of all realizations of `A`, together with the coercion morphisms
    is a category (whose class inherits from
    :class:`~sage.categories.realizations.Category_realization_of_parent`)::

        sage: A.Realizations()                                                          # needs sage.modules
        Category of realizations of
         The subset algebra of {1, 2, 3} over Rational Field

    The various parent realizing `A` belong to this category::

        sage: A.F() in A.Realizations()                                                 # needs sage.modules
        True

    `A` itself is in the category of algebras with realizations::

        sage: A in Algebras(QQ).WithRealizations()                                      # needs sage.modules
        True

    The (mostly technical) ``WithRealizations`` categories are the
    analogs of the ``*WithSeveralBases`` categories in
    MuPAD-Combinat. They provide support tools for handling the
    different realizations and the morphisms between them.

    Typically, ``VectorSpaces(QQ).FiniteDimensional().WithRealizations()``
    will eventually be in charge, whenever a coercion `\phi: A\mapsto B` is
    registered, to register `\phi^{-1}` as coercion `B \mapsto A`
    if there is none defined yet. To achieve this,
    ``FiniteDimensionalVectorSpaces`` would provide a nested class
    ``WithRealizations`` implementing the appropriate logic.

    ``WithRealizations`` is a :mod:`regressive covariant functorial
    construction <sage.categories.covariant_functorial_construction>`.
    On our example, this simply means that `A` is automatically in the
    category of rings with realizations (covariance)::

        sage: A in Rings().WithRealizations()                                           # needs sage.modules
        True

    and in the category of algebras (regressiveness)::

        sage: A in Algebras(QQ)                                                         # needs sage.modules
        True

    .. NOTE::

        For ``C`` a category, ``C.WithRealizations()`` in fact calls
        ``sage.categories.with_realizations.WithRealizations(C)``. The
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
        Category of graded Hopf algebras with basis over Rational Field with realizations
        sage: C.super_categories()
        [Join of Category of Hopf algebras over Rational Field
             and Category of graded algebras over Rational Field
             and Category of graded coalgebras over Rational Field]
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

            sage: C = GradedHopfAlgebrasWithBasis(QQ).WithRealizations(); C  # indirect doctest
            Category of graded Hopf algebras with basis over Rational Field with realizations
        """
        s = repr(self.base_category())
        return s + " with realizations"
