# -*- coding: utf-8 -*-
r"""
The category of ChowSchemes

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""

# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category import Category
from sage.categories.category_types import Category_over_base
from sage.categories.sets_cat import Sets


class ChowSchemes(Category):
    r"""
    Construct a category of ChowSchemes.
    """

    @staticmethod
    def __classcall_private__(cls, X=None):
        if X is not None:
            from .scheme import is_chowScheme
            if not is_chowScheme(X):
                X = ChowSchemes()(X)
            return ChowSchemes_over_base(X)
        else:
            return super(ChowSchemes, cls).__classcall__(cls)

    def super_categories(self):
        return [Sets()]

    def _call_(self, x):
        from .scheme import is_chowScheme
        if is_chowScheme(x):
            return x
        from .morphism import is_chowSchemeMorphism
        if is_chowSchemeMorphism(x):
            return x
        else:
            m = "Can't create an object or morphism in %s from %s" % (self, x)
            raise TypeError(m)



#############################################################
# Schemes over a given base scheme.
#############################################################

# noinspection PyAbstractClass
class ChowSchemes_over_base(Category_over_base):
    r"""
    The category of schemes over a given base scheme.
    """

    def base_chowscheme(self):
        return self.base()

    def super_categories(self):
        return [ChowSchemes()]

    def _repr_object_names(self):
        return "chowschemes over %s" % self.base_chowscheme()
