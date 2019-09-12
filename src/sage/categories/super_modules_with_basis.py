r"""
Super modules with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory

class SuperModulesWithBasis(SuperModulesCategory):
    """
    The category of super modules with a distinguished basis.

    An `R`-*super module with a distinguished basis* is an
    `R`-super module equipped with an `R`-module basis whose elements are
    homogeneous.

    EXAMPLES::

        sage: C = GradedModulesWithBasis(QQ); C
        Category of graded vector spaces with basis over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of filtered vector spaces with basis over Rational Field,
         Category of graded modules with basis over Rational Field,
         Category of graded vector spaces over Rational Field]
        sage: C is ModulesWithBasis(QQ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:
        def _even_odd_on_basis(self, m):
            """
            Return the parity of the basis element indexed by ``m``.

            OUTPUT:

            ``0`` if ``m`` is for an even element or ``1`` if ``m``
            is for an odd element.

            .. NOTE::

                The default implementation assumes that the even/odd is
                determined by the parity of :meth:`degree`.

                Overwrite this method if the even/odd behavior is desired
                to be independent.

            EXAMPLES::

                sage: Q = QuadraticForm(QQ, 2, [1,2,3])
                sage: C.<x,y> = CliffordAlgebra(Q)
                sage: C._even_odd_on_basis((0,))
                1
                sage: C._even_odd_on_basis((0,1))
                0
            """
            return self.degree_on_basis(m) % 2

    class ElementMethods:
        def is_super_homogeneous(self):
            r"""
            Return whether this element is homogeneous, in the sense
            of a super module (i.e., is even or odd).

            EXAMPLES::

                sage: Q = QuadraticForm(QQ, 2, [1,2,3])
                sage: C.<x,y> = CliffordAlgebra(Q)
                sage: a = x + y
                sage: a.is_super_homogeneous()
                True
                sage: a = x*y + 4
                sage: a.is_super_homogeneous()
                True
                sage: a = x*y + x - 3*y + 4
                sage: a.is_super_homogeneous()
                False

            The exterior algebra has a `\ZZ` grading, which induces the
            `\ZZ / 2\ZZ` grading. However the definition of homogeneous
            elements differs because of the different gradings::

                sage: E.<x,y> = ExteriorAlgebra(QQ)
                sage: a = x*y + 4
                sage: a.is_super_homogeneous()
                True
                sage: a.is_homogeneous()
                False
            """
            even_odd = self.parent()._even_odd_on_basis
            degree = None
            for m in self.support():
                if degree is None:
                    degree = even_odd(m)
                else:
                    if degree != even_odd(m):
                        return False
            return True

        def is_even_odd(self):
            """
            Return ``0`` if ``self`` is an even element and ``1`` if
            ``self`` is an odd element.

            EXAMPLES::

                sage: Q = QuadraticForm(QQ, 2, [1,2,3])
                sage: C.<x,y> = CliffordAlgebra(Q)
                sage: a = x + y
                sage: a.is_even_odd()
                1
                sage: a = x*y + 4
                sage: a.is_even_odd()
                0
                sage: a = x + 4
                sage: a.is_even_odd()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

                sage: E.<x,y> = ExteriorAlgebra(QQ)
                sage: (x*y).is_even_odd()
                0
            """
            if not self.support():
                raise ValueError("the zero element does not have a well-defined degree")
            if not self.is_super_homogeneous():
                raise ValueError("element is not homogeneous")
            return self.parent()._even_odd_on_basis(self.leading_support())

        def even_component(self):
            """
            Return the even component of ``self``.

            EXAMPLES::

                sage: Q = QuadraticForm(QQ, 2, [1,2,3])
                sage: C.<x,y> = CliffordAlgebra(Q)
                sage: a = x*y + x - 3*y + 4
                sage: a.even_component()
                x*y + 4

            TESTS:

            Check that this really return ``A.zero()`` and not a plain ``0``::

                sage: a = x + y
                sage: a.even_component().parent() is C
                True
            """
            even_odd = self.parent()._even_odd_on_basis
            return self.parent().sum_of_terms((i, c)
                                              for (i, c) in self
                                              if even_odd(i) == 0)

        def odd_component(self):
            """
            Return the odd component of ``self``.

            EXAMPLES::

                sage: Q = QuadraticForm(QQ, 2, [1,2,3])
                sage: C.<x,y> = CliffordAlgebra(Q)
                sage: a = x*y + x - 3*y + 4
                sage: a.odd_component()
                x - 3*y

            TESTS:

            Check that this really return ``A.zero()`` and not a plain ``0``::

                sage: a = x*y
                sage: a.odd_component().parent() is C
                True
            """
            even_odd = self.parent()._even_odd_on_basis
            return self.parent().sum_of_terms((i, c)
                                              for (i, c) in self
                                              if even_odd(i) == 1)

