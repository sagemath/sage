"""
Algebra Functorial Construction

AUTHORS:

- Nicolas M. Thiery (2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.covariant_functorial_construction import CovariantFunctorialConstruction, CovariantConstructionCategory
from sage.categories.category_types import Category_over_base_ring

class AlgebraFunctor(CovariantFunctorialConstruction):
    """
    A singleton class for the algebra functor.
    """
    _functor_name = "algebra"
    _functor_category = "Algebras"

    def __init__(self, base_ring):
        """
        EXAMPLES::

            sage: from sage.categories.algebra_functor import AlgebraFunctor
            sage: F = AlgebraFunctor(QQ); F
            The algebra functorial construction
            sage: TestSuite(F).run()
        """
        from sage.categories.rings import Rings
        assert base_ring in Rings()
        self._base_ring = base_ring

    def base_ring(self):
        """
        Return the base ring for this functor.

        EXAMPLES::

            sage: from sage.categories.algebra_functor import AlgebraFunctor
            sage: AlgebraFunctor(QQ).base_ring()
            Rational Field
        """
        return self._base_ring

class AlgebrasCategory(CovariantConstructionCategory, Category_over_base_ring):
    """
    An abstract base class for categories of monoid algebras,
    groups algebras, and the like.

    .. SEEALSO::

        - :meth:`Sets.ParentMethods.algebra`
        - :meth:`Sets.SubcategoryMethods.Algebras`
        - :class:`~sage.categories.covariant_functorial_construction.CovariantFunctorialConstruction`

    INPUT:

     - ``base_ring`` -- a ring

    EXAMPLES::

        sage: C = Monoids().Algebras(QQ); C
        Category of monoid algebras over Rational Field
        sage: C = Groups().Algebras(QQ); C
        Category of group algebras over Rational Field

        sage: C._short_name()
        'Algebras'
        sage: latex(C) # todo: improve that
        \mathbf{Algebras}(\mathbf{Groups})
    """

    _functor_category = "Algebras"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Semigroups().Algebras(QQ) # indirect doctest
            Category of semigroup algebras over Rational Field
        """
        return "{} algebras over {}".format(self.base_category()._repr_object_names()[:-1],
                                            self.base_ring())

