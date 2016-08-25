# -*- coding: utf-8 -*-
r"""
Aperiodic semigroups
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

class AperiodicSemigroups(CategoryWithAxiom):
    def extra_super_categories(self):
        r"""
        Implement the fact that an aperiodic semigroup is `H`-trivial.

        EXAMPLES::

            sage: Semigroups().Aperiodic().extra_super_categories()
            [Category of h trivial semigroups]
        """
        return [Semigroups().HTrivial()]
