"""
Examples of manifolds
"""
# ****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.manifolds import Manifolds


class Plane(UniqueRepresentation, Parent):
    r"""
    An example of a manifold: the `n`-dimensional plane.

    This class illustrates a minimal implementation of a manifold.

    EXAMPLES::

        sage: from sage.categories.manifolds import Manifolds
        sage: M = Manifolds(QQ).example(); M
        An example of a Rational Field manifold: the 3-dimensional plane

        sage: M.category()
        Category of manifolds over Rational Field

    We conclude by running systematic tests on this manifold::

        sage: TestSuite(M).run()
    """

    def __init__(self, n=3, base_ring=None):
        r"""
        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: M = Manifolds(QQ).example(6); M
            An example of a Rational Field manifold: the 6-dimensional plane

        TESTS::

            sage: TestSuite(M).run()
        """
        self._n = n
        Parent.__init__(self, base=base_ring, category=Manifolds(base_ring))

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds(QQ).example()
            An example of a Rational Field manifold: the 3-dimensional plane
        """
        return "An example of a {} manifold: the {}-dimensional plane".format(
                self.base_ring(), self._n)

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: M = Manifolds(QQ).example()
            sage: M.dimension()
            3
        """
        return self._n

    def an_element(self):
        r"""
        Return an element of the manifold, as per
        :meth:`Sets.ParentMethods.an_element`.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: M = Manifolds(QQ).example()
            sage: M.an_element()
            (0, 0, 0)
        """
        zero = self.base_ring().zero()
        return self(tuple([zero]*self._n))

    Element = ElementWrapper

Example = Plane
