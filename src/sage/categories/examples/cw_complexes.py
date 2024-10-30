# sage_setup: distribution = sagemath-categories
"""
Examples of CW complexes
"""
# ****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.categories.cw_complexes import CWComplexes
from sage.sets.family import Family


class Surface(UniqueRepresentation, Parent):
    r"""
    An example of a CW complex: a (2-dimensional) surface.

    This class illustrates a minimal implementation of a CW complex.

    EXAMPLES::

        sage: from sage.categories.cw_complexes import CWComplexes
        sage: X = CWComplexes().example(); X
        An example of a CW complex: the surface given by the boundary map (1, 2, 1, 2)

        sage: X.category()
        Category of finite finite dimensional CW complexes

    We conclude by running systematic tests on this manifold::

        sage: TestSuite(X).run()
    """
    def __init__(self, bdy=(1, 2, 1, 2)):
        r"""
        EXAMPLES::

            sage: from sage.categories.cw_complexes import CWComplexes
            sage: X = CWComplexes().example((1, 2)); X
            An example of a CW complex: the surface given by the boundary map (1, 2)

        TESTS::

            sage: TestSuite(X).run()
        """
        self._bdy = bdy
        self._edges = frozenset(bdy)
        Parent.__init__(self, category=CWComplexes().Finite())

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.categories.cw_complexes import CWComplexes
            sage: CWComplexes().example()
            An example of a CW complex: the surface given by the boundary map (1, 2, 1, 2)
        """
        return "An example of a CW complex: the surface given by the boundary map {}".format(self._bdy)

    def cells(self):
        """
        Return the cells of ``self``.

        EXAMPLES::

            sage: from sage.categories.cw_complexes import CWComplexes
            sage: X = CWComplexes().example()
            sage: C = X.cells()
            sage: sorted((d, C[d]) for d in C.keys())
            [(0, (0-cell v,)),
             (1, (0-cell e1, 0-cell e2)),
             (2, (2-cell f,))]
        """
        d = {0: (self.element_class(self, 0, 'v'),)}
        d[1] = tuple([self.element_class(self, 0, 'e'+str(e)) for e in self._edges])
        d[2] = (self.an_element(),)
        return Family(d)

    def an_element(self):
        r"""
        Return an element of the CW complex, as per
        :meth:`Sets.ParentMethods.an_element`.

        EXAMPLES::

            sage: from sage.categories.cw_complexes import CWComplexes
            sage: X = CWComplexes().example()
            sage: X.an_element()
            2-cell f
        """
        return self.element_class(self, 2, 'f')

    class Element(Element):
        """
        A cell in a CW complex.
        """
        def __init__(self, parent, dim, name):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: from sage.categories.cw_complexes import CWComplexes
                sage: X = CWComplexes().example()
                sage: f = X.an_element()
                sage: TestSuite(f).run()
            """
            Element.__init__(self, parent)
            self._dim = dim
            self._name = name

        def _repr_(self):
            """
            Return a string representation of ``self``.

            EXAMPLES::

                sage: from sage.categories.cw_complexes import CWComplexes
                sage: X = CWComplexes().example()
                sage: X.an_element()
                2-cell f
            """
            return "{}-cell {}".format(self._dim, self._name)

        def __eq__(self, other):
            """
            Check equality.

            EXAMPLES::

                sage: from sage.categories.cw_complexes import CWComplexes
                sage: X = CWComplexes().example()
                sage: f = X.an_element()
                sage: f == X(2, 'f')
                True
                sage: e1 = X(1, 'e1')
                sage: e1 == f
                False
            """
            return (isinstance(other, Surface.Element)
                    and self.parent() is other.parent()
                    and self._dim == other._dim
                    and self._name == other._name)

        def dimension(self):
            """
            Return the dimension of ``self``.

            EXAMPLES::

                sage: from sage.categories.cw_complexes import CWComplexes
                sage: X = CWComplexes().example()
                sage: f = X.an_element()
                sage: f.dimension()
                2
            """
            return self._dim


Example = Surface
