r"""
Algebra ideals
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from .category_types import Category_ideal
from .algebra_modules import AlgebraModules


class AlgebraIdeals(Category_ideal):
    """
    The category of two-sided ideals in a fixed algebra `A`.

    EXAMPLES::

        sage: AlgebraIdeals(QQ['a'])
        Category of algebra ideals in Univariate Polynomial Ring in a over Rational Field

    .. TODO::

        - Add support for non commutative rings (this is currently not
          supported by the subcategory :class:`AlgebraModules`).
        - Make ``AlgebraIdeals(R)``, return ``CommutativeAlgebraIdeals(R)``
          when ``R`` is commutative.
        - If useful, implement ``AlgebraLeftIdeals`` and
          ``AlgebraRightIdeals`` of which ``AlgebraIdeals``
          would be a subcategory.
    """
    def __init__(self, A):
        """
        EXAMPLES::

            sage: AlgebraIdeals(QQ['a'])
            Category of algebra ideals in Univariate Polynomial Ring in a over Rational Field
            sage: AlgebraIdeals(QQ)
            Traceback (most recent call last):
            ...
            TypeError: A (=Rational Field) must be an algebra

        TESTS::

            sage: TestSuite(AlgebraIdeals(QQ['a'])).run()
        """
        from sage.algebras.algebra import is_Algebra
        if not is_Algebra(A): # A not in Algebras() ?
            raise TypeError("A (=%s) must be an algebra"%A)
        Category_ideal.__init__(self, A)

    def algebra(self):
        """
        EXAMPLES::

            sage: AlgebraIdeals(QQ['x']).algebra()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self.ambient()

    def super_categories(self):
        """
        The category of algebra modules should be a super category of this category.

        However, since algebra modules are currently only available over commutative rings,
        we have to omit it if our ring is non-commutative.

        EXAMPLES::

            sage: AlgebraIdeals(QQ['x']).super_categories()
            [Category of algebra modules over Univariate Polynomial Ring in x over Rational Field]
            sage: C = AlgebraIdeals(FreeAlgebra(QQ,2,'a,b'))
            sage: C.super_categories()
            []

        """
        R = self.algebra()
        try:
            if R.is_commutative():
                return [AlgebraModules(R)]
        except (AttributeError,NotImplementedError):
            pass
        return []
