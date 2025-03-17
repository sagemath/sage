# sage_setup: distribution = sagemath-categories
r"""
Division rings
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.rings import Rings


class DivisionRings(CategoryWithAxiom):
    """
    The category of division rings.

    A division ring (or skew field) is a not necessarily commutative
    ring where all nonzero elements have multiplicative inverses

    EXAMPLES::

      sage: DivisionRings()
      Category of division rings
      sage: DivisionRings().super_categories()
      [Category of domains]

    TESTS::

        sage: TestSuite(DivisionRings()).run()
    """

    # This information could be deduced from the name. It's set here
    # just to get ``division rings`` for the name of the objects of
    # this category and not ``no zero divisors division rings``. See
    # :meth:`Category_with_axiom._repr_object_names` and the ``named``
    # option of :meth:`Category_with_axiom._without_axioms
    _base_category_class_and_axiom = (Rings, "Division")

    def extra_super_categories(self):
        r"""
        Return the :class:`Domains` category.

        This method specifies that a division ring has no zero
        divisors, i.e. is a domain.

        .. SEEALSO::

            The :ref:`axioms-deduction-rules` section in the
            documentation of axioms

        EXAMPLES::

            sage: DivisionRings().extra_super_categories()
            (Category of domains,)
            sage: "NoZeroDivisors" in DivisionRings().axioms()
            True
        """
        return (Rings().NoZeroDivisors(),)

    Commutative = LazyImport('sage.categories.fields', 'Fields', at_startup=True)

    def Finite_extra_super_categories(self):
        r"""
        Return extraneous super categories for ``DivisionRings().Finite()``.

        EXAMPLES:

        Any field is a division ring::

            sage: Fields().is_subcategory(DivisionRings())
            True

        This methods specifies that, by Weddeburn theorem, the
        reciprocal holds in the finite case: a finite division ring is
        commutative and thus a field::

            sage: DivisionRings().Finite_extra_super_categories()
            (Category of commutative magmas,)
            sage: DivisionRings().Finite()
            Category of finite enumerated fields

        .. WARNING::

            This is not implemented in
            ``DivisionRings.Finite.extra_super_categories`` because
            the categories of finite division rings and of finite
            fields coincide. See the section
            :ref:`axioms-deduction-rules` in the documentation of
            axioms.

        TESTS::

            sage: DivisionRings().Finite() is Fields().Finite()
            True

        This works also for subcategories::

            sage: class Foo(Category):
            ....:     def super_categories(self): return [DivisionRings()]
            sage: Foo().Finite().is_subcategory(Fields())
            True
        """
        from sage.categories.magmas import Magmas
        return (Magmas().Commutative(),)

    class ParentMethods:
        pass

    class ElementMethods:
        pass
