"""
Graded Lie Conformal Algebras

AUTHORS:

- Reimundo Heluani (2019-10-05): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.misc.cachefunc import cached_method

class GradedLieConformalAlgebrasCategory(GradedModulesCategory):
    @cached_method
    def Super(self, base_ring=None):
        r"""
        Return the super-analogue category of ``self``.

        INPUT:

        - ``base_ring`` -- this is ignored

        EXAMPLES::

            sage: C = LieConformalAlgebras(QQbar)
            sage: C.Graded().Super() is C.Super().Graded()
            True
            sage: Cp = C.WithBasis()
            sage: Cp.Graded().Super() is Cp.Super().Graded()
            True
        """
        return self.base_category().Super(base_ring).Graded()

    def _repr_object_names(self):
        """
        The names of the objects of ``self``.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).Graded()
            Category of H-graded Lie conformal algebras over Algebraic Field

            sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated().Graded()
            Category of H-graded finitely generated Lie conformal algebras with basis over Algebraic Field
        """
        return "H-graded {}".format(self.base_category()._repr_object_names())

class GradedLieConformalAlgebras(GradedLieConformalAlgebrasCategory):
    """
    The category of graded Lie conformal algebras.

    EXAMPLES::

        sage: C = LieConformalAlgebras(QQbar).Graded(); C
        Category of H-graded Lie conformal algebras over Algebraic Field

        sage: CS = LieConformalAlgebras(QQ).Graded().Super(); CS
        Category of H-graded super Lie conformal algebras over Rational Field
        sage: CS is LieConformalAlgebras(QQ).Super().Graded()
        True
    """
