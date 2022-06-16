r"""
Drinfeld modules
"""
#*****************************************************************************
#  Copyright (C) 2022      Xavier Caruso <xavier.caruso@normalesup.org>
#                          Antoine Leudière <antoine.leudiere@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import CategoryWithParameters

# from sage.misc.cachefunc import cached_method
# from sage.categories.basic import Fields

class DrinfeldModules(CategoryWithParameters):
    r"""
    The category of Drinfeld modules.

    EXAMPLES:

    We create the category of function fields::

        sage: C = FunctionFields()
        sage: C
        Category of function fields

    TESTS::

        sage: TestSuite(FunctionFields()).run()
    """

    def __init__(self, domain, codomain):
        r"""
        """
        self._domain = domain
        self._codomain = codomain

    def _call_(self, phi_X):
        r"""
        Constructs an object in this category from the data in ``x``,
        or throws a TypeError.
        """
        return FiniteDrinfeldModule(self._domain, self._codomain(x), characteristic)

    def super_categories(self):
        return []

    def _repr_(self):
        return f'Category of Drinfeld modules:\n'\
                f'  Domain: {self._domain}\n' \
                f'  Codomain: {self._codomain}'

    def _make_named_class_key(self, name):
        return self._domain.category()

    class ParentMethods:
        pass

    class ElementMethods:
        pass
