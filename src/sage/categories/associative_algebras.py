r"""
Associative algebras
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.magmas import Magmas
from sage.categories.magmatic_algebras import MagmaticAlgebras

class AssociativeAlgebras(CategoryWithAxiom_over_base_ring):
    r"""
    The category of associative algebras over a given base ring.

    An associative algebra over a ring `R` is a module over `R` which
    is also a not necessarily unital ring.

    .. WARNING::

        Until :trac:`15043` is implemented, :class:`Algebras` is the
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

    class ElementMethods:
        """
        An abstract class for elements of an associative algebra.

        .. NOTE::

            ``Magmas.Element.__mul__`` is preferable to
            ``Modules.Element.__mul__`` since the later does not
            handle products of two elements of ``self``.

        TESTS::

            sage: A = AlgebrasWithBasis(QQ).example()
            sage: a = A.an_element()
            sage: a
            2*B[word: ] + 2*B[word: a] + 3*B[word: b]
            sage: a.__mul__(a)
            4*B[word: ] + 8*B[word: a] + 4*B[word: aa] + 6*B[word: ab] + 12*B[word: b] + 6*B[word: ba] + 9*B[word: bb]
        """
        __mul__ = Magmas.ElementMethods.__mul__.im_func

#        __imul__ = __mul__


    Unital = LazyImport('sage.categories.algebras', 'Algebras', at_startup=True)
