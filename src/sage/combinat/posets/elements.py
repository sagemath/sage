# sage_setup: distribution = sagemath-graphs
# sage.doctest: needs sage.modules
r"""
Elements of posets, lattices, semilattices, etc.
"""
# ****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.element import Element
from sage.structure.element import have_same_parent


class PosetElement(Element):

    def __init__(self, poset, element, vertex) -> None:
        r"""
        Establish the parent-child relationship between ``poset``
        and ``element``, where ``element`` is associated to the
        vertex ``vertex`` of the Hasse diagram of the poset.

        INPUT:

        - ``poset`` -- a poset object

        - ``element`` -- any object

        - ``vertex`` -- a vertex of the Hasse diagram of the poset

        TESTS::

            sage: from sage.combinat.posets.elements import PosetElement
            sage: P = Poset([[1,2],[4],[3],[4],[]], facade = False)
            sage: e = P(0)
            sage: e.parent() is P
            True
            sage: TestSuite(e).run()
        """
        Element.__init__(self, poset)
        if isinstance(element, self.parent().element_class):
            self.element = element.element
        else:
            self.element = element
        self.vertex = vertex

    def __hash__(self):
        r"""
        TESTS::

            sage: P = Poset([[1,2],[4],[3],[4],[]], facade = False)
            sage: e = P(0)
            sage: hash(e)
            0
        """
        return hash(self.element)

    def _repr_(self):
        """
        TESTS::

            sage: Poset([[1,2],[4],[3],[4],[]], facade = False)(0)._repr_()
            '0'
        """
        return "%s" % str(self.element)

    def _latex_(self):
        r"""
        Return the latex code of the poset element.

        EXAMPLES::

            sage: m = matrix(2, [1,2,3,4])
            sage: m.set_immutable()
            sage: P = Poset(([m],[]), facade=False)
            sage: [e] = P
            sage: type(e)
            <class 'sage.combinat.posets.posets.FinitePoset_with_category.element_class'>
            sage: latex(e)                 #indirect doctest
            \left(\begin{array}{rr}
            1 & 2 \\
            3 & 4
            \end{array}\right)
        """
        from sage.misc.latex import latex
        return latex(self.element)

    def __eq__(self, other):
        """
        TESTS::

            sage: P = Poset([["a","b"],["d"],["c"],["d"],[]], facade = False)
            sage: Q = Poset([["a","b"],["d"],["c"],[],[]], facade = False)
            sage: P(0).__eq__(P(4))
            False
            sage: from sage.combinat.posets.elements import PosetElement
            sage: PosetElement(P,0,"c") == PosetElement(P,0,"c")
            True
            sage: PosetElement(P,0,"c") == PosetElement(Q,0,"c")
            False
            sage: PosetElement(P,0,"b") == PosetElement(P,0,"c")
            False

        .. warning:: as an optimization, this only compares the parent
           and vertex, using the invariant that, in a proper poset
           element, ``self.element == other.element`` if and only
           ``self.vertex == other.vertex``::

            sage: PosetElement(P,1,"c") == PosetElement(P,0,"c")
            True

        Test that :issue:`12351` is fixed::

            sage: P(0) == int(0)
            False
        """
        # This should instead exploit unique representation, using
        # self is other, or best inherit __eq__ from there. But there
        # are issues around pickling and rich comparison functions.
        return have_same_parent(self, other) \
            and self.vertex == other.vertex

    def __ne__(self, other):
        r"""
        TESTS::

            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: P = Poset([["a","b"],["d"],["c"],["d"],[]])
            sage: P(0).__ne__(P(4))
            True
            sage: from sage.combinat.posets.elements import PosetElement
            sage: PosetElement(P,0,"c") != PosetElement(P,0,"c")
            False
            sage: PosetElement(P,0,"b") != PosetElement(P,0,"c")
            True

        For this one, see comment in :meth:`__eq__`::

            sage: PosetElement(P,1,"c") != PosetElement(P,0,"c")
            False
        """
        return not self == other

    def _cmp(self, other):
        """
        TESTS::

            sage: P = Poset([[1,2],[4],[3],[4],[]], facade = False)
            sage: P(0)._cmp(P(4))
            -1
            sage: P(4)._cmp(P(0))
            1
            sage: P(0)._cmp(P(0))
            0
            sage: P(1)._cmp(P(2))
        """
        return self.parent().compare_elements(self, other)

    def __lt__(self, other):
        """
        TESTS::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag, facade = False)
            sage: P(0) < P(1)
            False
            sage: P(4) < P(1)
            False
            sage: P(0) < P(0)
            False
        """
        return self._cmp(other) == -1 or False

    def __le__(self, other):
        """
        TESTS::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag, facade = False)
            sage: P(1) <= P(0)
            False
            sage: P(0) <= P(1)
            False
            sage: P(0) <= P(3)
            True
            sage: P(0) <= P(0)
            True
        """
        return self == other or self._cmp(other) == -1 or False

    def __gt__(self, other):
        """
        TESTS::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag)
            sage: P(0).__gt__(P(5))
            False
            sage: P(5).__gt__(P(0))
            True
            sage: P(0).__gt__(P(0))
            False
        """
        return self._cmp(other) == 1 or False

    def __ge__(self, other):
        """
        TESTS::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag)
            sage: P(0).__ge__(P(5))
            False
            sage: P(5).__ge__(P(0))
            True
            sage: P(0).__ge__(P(0))
            True
        """
        return self == other or self._cmp(other) == 1 or False


class MeetSemilatticeElement(PosetElement):
    def __mul__(self, other):
        r"""
        Return the meet of ``self`` and ``other`` in the lattice.

        EXAMPLES::

            sage: D = posets.DiamondPoset(5, facade=False)
            sage: D(1) * D(2)
            0
            sage: D(1) * D(1)
            1
            sage: D(1) * D(0)
            0
            sage: D(1) * D(4)
            1
        """
        return self.parent().meet(self, other)


class JoinSemilatticeElement(PosetElement):
    def __add__(self, other):
        r"""
        Return the join of ``self`` and ``other`` in the lattice.

        EXAMPLES::

            sage: D = posets.DiamondPoset(5,facade=False)
            sage: D(1) + D(2)
            4
            sage: D(1) + D(1)
            1
            sage: D(1) + D(4)
            4
            sage: D(1) + D(0)
            1
        """
        return self.parent().join(self, other)


class LatticePosetElement(MeetSemilatticeElement, JoinSemilatticeElement):
    pass
