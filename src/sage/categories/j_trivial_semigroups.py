# sage_setup: distribution = sagemath-categories
r"""
J-trivial semigroups
"""
#*****************************************************************************
#  Copyright (C) 2016 Nicolas M. Thiéry <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.semigroups import Semigroups


class JTrivialSemigroups(CategoryWithAxiom):
    def extra_super_categories(self):
        """
        Implement the fact that a `J`-trivial semigroup is `L` and `R`-trivial.

        EXAMPLES::

            sage: Semigroups().JTrivial().extra_super_categories()
            [Category of l trivial semigroups, Category of r trivial semigroups]
        """
        return [Semigroups().LTrivial(), Semigroups().RTrivial()]
