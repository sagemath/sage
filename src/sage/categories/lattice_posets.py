# sage_setup: distribution = sagemath-categories
r"""
Lattice posets
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category import Category
from sage.categories.posets import Posets


class LatticePosets(Category):
    r"""
    The category of lattices, i.e. partially ordered sets in which any
    two elements have a unique supremum (the elements' least upper
    bound; called their *join*) and a unique infimum (greatest lower bound;
    called their *meet*).

    EXAMPLES::

        sage: LatticePosets()
        Category of lattice posets
        sage: LatticePosets().super_categories()
        [Category of posets]
        sage: LatticePosets().example()
        NotImplemented

    .. SEEALSO:: :class:`~sage.categories.posets.Posets`, :class:`FiniteLatticePosets`, :func:`LatticePoset`

    TESTS::

        sage: C = LatticePosets()
        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        r"""
        Return a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: LatticePosets().super_categories()
            [Category of posets]
        """
        return [Posets()]

    Finite = LazyImport('sage.categories.finite_lattice_posets', 'FiniteLatticePosets')

    class ParentMethods:

        @abstract_method
        def meet(self, x, y):
            """
            Return the meet of `x` and `y` in this lattice.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: D = LatticePoset((divisors(30), attrcall("divides")))             # needs sage.graphs sage.modules
                sage: D.meet( D(6), D(15) )                                             # needs sage.graphs sage.modules
                3
            """

        @abstract_method
        def join(self, x, y):
            """
            Return the join of `x` and `y` in this lattice.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: D = LatticePoset((divisors(60), attrcall("divides")))             # needs sage.graphs sage.modules
                sage: D.join( D(6), D(10) )                                             # needs sage.graphs sage.modules
                30
            """
