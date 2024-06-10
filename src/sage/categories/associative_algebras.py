# sage_setup: distribution = sagemath-categories
r"""
Associative algebras
"""
# ****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.magmatic_algebras import MagmaticAlgebras


class AssociativeAlgebras(CategoryWithAxiom_over_base_ring):
    r"""
    The category of associative algebras over a given base ring.

    An associative algebra over a ring `R` is a module over `R` which
    is also a not necessarily unital ring.

    .. WARNING::

        Until :issue:`15043` is implemented, :class:`Algebras` is the
        category of associative unital algebras; thus, unlike the name
        suggests, :class:`AssociativeAlgebras` is not a subcategory of
        :class:`Algebras` but of
        :class:`~.magmatic_algebras.MagmaticAlgebras`.

    EXAMPLES::

        sage: from sage.categories.associative_algebras import AssociativeAlgebras
        sage: C = AssociativeAlgebras(ZZ); C
        Category of associative algebras over Integer Ring

    TESTS::

        sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
        sage: C is MagmaticAlgebras(ZZ).Associative()
        True
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (MagmaticAlgebras, "Associative")

    Unital = LazyImport('sage.categories.algebras', 'Algebras', at_startup=True)
