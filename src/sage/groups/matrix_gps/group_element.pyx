"""
Matrix Group Elements

EXAMPLES::

    sage: F = GF(3); MS = MatrixSpace(F, 2, 2)
    sage: gens = [MS([[1,0], [0,1]]), MS([[1,1], [0,1]])]
    sage: G = MatrixGroup(gens); G
    Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1] )
    sage: g = G([[1,1], [0,1]])
    sage: h = G([[1,2], [0,1]])
    sage: g*h
    [1 0]
    [0 1]

You cannot add two matrices, since this is not a group operation.
You can coerce matrices back to the matrix space and add them
there::

    sage: g + h
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +:
    'Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )' and
    'Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )'

    sage: g.matrix() + h.matrix()
    [2 0]
    [0 2]

Similarly, you cannot multiply group elements by scalars but you can
do it with the underlying matrices::

    sage: 2*g
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for *: 'Integer Ring'
    and 'Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1] )'

AUTHORS:

- David Joyner (2006-05): initial version David Joyner

- David Joyner (2006-05): various modifications to address William
  Stein's TODO's.

- William Stein (2006-12-09): many revisions.

- Volker Braun (2013-1) port to new Parent, libGAP.

- Travis Scrimshaw (2016-01): reworks class hierarchy in order
  to cythonize
"""

#*****************************************************************************
#       Copyright (C) 2006      David Joyner and William Stein <wstein@gmail.com>
#                     2013      Volker Braun <vbraun.name@gmail.com>
#                     2016      Travis Scrimshaw <tscrimsh at umn.edu>
#                     2016-2018 Jeroen Demeyer
#                     2023      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport richcmp


try:
    from sage.groups.matrix_gps.group_element_gap import MatrixGroupElement_gap
except ImportError:
    MatrixGroupElement_gap = ()


cpdef is_MatrixGroupElement(x):
    """
    Test whether ``x`` is a matrix group element.

    INPUT:

    - ``x`` -- anything

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.groups.matrix_gps.group_element import is_MatrixGroupElement
        sage: is_MatrixGroupElement('helloooo')
        False

        sage: G = GL(2,3)
        sage: is_MatrixGroupElement(G.an_element())
        True
    """
    return isinstance(x, (MatrixGroupElement_generic, MatrixGroupElement_gap))

###################################################################
#
# Matrix group elements implemented in Sage
#
###################################################################

cdef class MatrixGroupElement_generic(MultiplicativeGroupElement):
    """
    Element of a matrix group over a generic ring.

    The group elements are implemented as Sage matrices.

    INPUT:

    - ``M`` -- a matrix

    - ``parent`` -- the parent

    - ``check`` -- boolean (default: ``True``); if ``True``, then
      do some type checking

    - ``convert`` -- boolean (default: ``True``); if ``True``, then
      convert ``M`` to the right matrix space

    EXAMPLES::

        sage: W = CoxeterGroup(['A',3], base_ring=ZZ)                                   # needs sage.graphs
        sage: g = W.an_element(); g                                                     # needs sage.graphs
        [ 0  0 -1]
        [ 1  0 -1]
        [ 0  1 -1]
    """
    def __init__(self, parent, M, check=True, convert=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)                               # needs sage.graphs
            sage: g = W.an_element()                                                    # needs sage.graphs
            sage: TestSuite(g).run()                                                    # needs sage.graphs
        """
        if convert:
            M = parent.matrix_space()(M)
        if check:
            if not isinstance(M, Matrix):
                raise TypeError('M must be a matrix')
            if M.parent() is not parent.matrix_space():
                raise TypeError('M must be a in the matrix space of the group')
            parent._check_matrix(M)
        super().__init__(parent)
        if M.is_immutable():
            self._matrix = M
        else:
            self._matrix = M.__copy__()
            self._matrix.set_immutable()

    def __hash__(self):
        r"""
        TESTS::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)                               # needs sage.graphs
            sage: g = W.an_element()                                                    # needs sage.graphs
            sage: hash(g)                                                               # needs sage.graphs
            660522311176098153  # 64-bit
            -606138007          # 32-bit
        """
        return hash(self._matrix)

    def __reduce__(self):
        """
        Implement pickling.

        TESTS::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)                               # needs sage.graphs
            sage: g = W.an_element()                                                    # needs sage.graphs
            sage: loads(g.dumps()) == g                                                 # needs sage.graphs
            True
        """
        return (_unpickle_generic_element, (self.parent(), self._matrix,))

    def _repr_(self):
        """
        Return string representation of this matrix.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)                               # needs sage.graphs
            sage: W.an_element()                                                        # needs sage.graphs
            [ 0  0 -1]
            [ 1  0 -1]
            [ 0  1 -1]
        """
        return str(self._matrix)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)                               # needs sage.graphs
            sage: g = W.an_element()                                                    # needs sage.graphs
            sage: latex(g)                                                              # needs sage.graphs
            \left(\begin{array}{rrr}
            0 & 0 & -1 \\
            1 & 0 & -1 \\
            0 & 1 & -1
            \end{array}\right)
        """
        return self._matrix._latex_()

    cpdef _act_on_(self, x, bint self_on_left):
        """
        EXAMPLES::

            sage: # needs sage.combinat sage.libs.gap
            sage: W = CoxeterGroup(['A',4], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: g * vector([1,1,1,1])
            (0, 1, 1, 1)
            sage: v = vector([3,2,1,-1])
            sage: g = W.gen(1)
            sage: v * g == v * g.matrix()   # indirect doctest
            True
        """
        if not is_MatrixGroupElement(x) and x not in self.parent().base_ring():
            try:
                if self_on_left:
                    return self._matrix * x
                else:
                    return x * self._matrix
            except TypeError:
                return None

    cpdef _richcmp_(self, other, int op):
        """
        EXAMPLES::

            sage: # needs sage.combinat sage.libs.gap
            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: TestSuite(g).run()
            sage: h = W.gen(0) * W.gen(1) * W.gen(2)
            sage: g == h
            True
            sage: a = W.gen(0)
            sage: a == g
            False
            sage: a != g
            True
        """
        cdef MatrixGroupElement_generic x = <MatrixGroupElement_generic>self
        cdef MatrixGroupElement_generic y = <MatrixGroupElement_generic>other
        return richcmp(x._matrix, y._matrix, op)

    cpdef list list(self):
        """
        Return list representation of this matrix.

        EXAMPLES::

            sage: # needs sage.combinat sage.libs.gap
            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: g
            [-1  1  0]
            [ 0  1  0]
            [ 0  0  1]
            sage: g.list()
            [[-1, 1, 0], [0, 1, 0], [0, 0, 1]]
        """
        return [r.list() for r in self._matrix.rows()]

    def matrix(self):
        """
        Obtain the usual matrix (as an element of a matrix space)
        associated to this matrix group element.

        One reason to compute the associated matrix is that matrices
        support a huge range of functionality.

        EXAMPLES::

            sage: # needs sage.combinat sage.libs.gap
            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: g.matrix()
            [-1  1  0]
            [ 0  1  0]
            [ 0  0  1]
            sage: parent(g.matrix())
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

        Matrices have extra functionality that matrix group elements
        do not have::

            sage: g.matrix().charpoly('t')                                              # needs sage.combinat sage.libs.gap
            t^3 - t^2 - t + 1
        """
        return self._matrix

    def _matrix_(self, base=None):
        """
        Method used by the :func:`matrix` constructor.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], base_ring=ZZ)                              # needs sage.graphs
            sage: g = W.gen(0)                                                          # needs sage.graphs
            sage: matrix(RDF, g)                                                        # needs sage.graphs
            [-1.0  1.0  0.0]
            [ 0.0  1.0  0.0]
            [ 0.0  0.0  1.0]
        """
        return self.matrix()

    cpdef _mul_(self, other):
        """
        Return the product of ``self`` and`` other``, which must
        have identical parents.

        EXAMPLES::

            sage: # needs sage.combinat sage.libs.gap
            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.gen(0)
            sage: h = W.an_element()
            sage: g * h
            [ 1  0  0]
            [ 1  0 -1]
            [ 0  1 -1]
        """
        cdef Parent parent = self.parent()
        cdef MatrixGroupElement_generic y = <MatrixGroupElement_generic>other
        cdef Matrix M = self._matrix * y._matrix
        # Make it immutable so the constructor doesn't make a copy
        M.set_immutable()
        return parent.element_class(parent, M, check=False, convert=False)

    def is_one(self):
        """
        Return whether ``self`` is the identity of the group.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: W = CoxeterGroup(['A',3])
            sage: g = W.gen(0)
            sage: g.is_one()
            False
            sage: W.an_element().is_one()
            False
            sage: W.one().is_one()
            True
        """
        return self._matrix.is_one()

    def __invert__(self):
        """
        Return the inverse group element.

        OUTPUT: a matrix group element

        EXAMPLES::

            sage: # needs sage.combinat sage.libs.gap
            sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
            sage: g = W.an_element()
            sage: ~g
            [-1  1  0]
            [-1  0  1]
            [-1  0  0]
            sage: g * ~g == W.one()
            True
            sage: ~g * g == W.one()
            True

            sage: # needs sage.combinat sage.libs.gap sage.rings.number_field
            sage: W = CoxeterGroup(['B',3])
            sage: W.base_ring()
            Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?
            sage: g = W.an_element()
            sage: ~g
            [-1  1  0]
            [-1  0  a]
            [-a  0  1]
        """
        cdef Parent parent = self.parent()
        cdef Matrix M = self._matrix
        # We have a special method for dense matrices over ZZ
        if M.base_ring() is ZZ and M.is_dense():
            M = M.inverse_of_unit()
        else:
            M = ~M
            if M.base_ring() is not parent.base_ring():
                M = M.change_ring(parent.base_ring())
        # Make it immutable so the constructor doesn't make a copy
        M.set_immutable()
        return parent.element_class(parent, M, check=False, convert=False)

    inverse = __invert__


def _unpickle_generic_element(G, mat):
    """
    Unpickle the element in ``G`` given by ``mat``.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: m1 = matrix(SR, [[1,2], [3,4]])
        sage: m2 = matrix(SR, [[1,3], [-1,0]])
        sage: G = MatrixGroup(m1, m2)
        sage: m = G.an_element()
        sage: from sage.groups.matrix_gps.group_element import _unpickle_generic_element
        sage: _unpickle_generic_element(G, m.matrix()) == m
        True
    """
    return G.element_class(G, mat, False, False)
