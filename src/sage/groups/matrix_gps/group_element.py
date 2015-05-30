"""
Matrix Group Elements

EXAMPLES::

    sage: F = GF(3); MS = MatrixSpace(F,2,2)
    sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
    sage: G = MatrixGroup(gens); G
    Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )
    sage: g = G([[1,1],[0,1]])
    sage: h = G([[1,2],[0,1]])
    sage: g*h
    [1 0]
    [0 1]

You cannot add two matrices, since this is not a group operation.
You can coerce matrices back to the matrix space and add them
there::

    sage: g + h
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand type(s) for +:
    'FinitelyGeneratedMatrixGroup_gap_with_category.element_class' and
    'FinitelyGeneratedMatrixGroup_gap_with_category.element_class'

    sage: g.matrix() + h.matrix()
    [2 0]
    [0 2]

Similarly, you cannot multiply group elements by scalars but you can
do it with the underlying matrices::

    sage: 2*g
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and 'Matrix group over Finite Field of size 3 with 2 generators (
    [1 0]  [1 1]
    [0 1], [0 1]
    )'

AUTHORS:

- David Joyner (2006-05): initial version David Joyner

- David Joyner (2006-05): various modifications to address William
  Stein's TODO's.

- William Stein (2006-12-09): many revisions.

- Volker Braun (2013-1) port to new Parent, libGAP.
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.interfaces.gap import gap
from sage.structure.element import MultiplicativeGroupElement
from sage.matrix.matrix import Matrix, is_Matrix
from sage.structure.factorization import Factorization
from sage.structure.element import have_same_parent
from sage.libs.gap.element import GapElement, GapElement_List
from sage.misc.cachefunc import cached_method
from sage.groups.libgap_wrapper import ElementLibGAP
from sage.groups.libgap_mixin import GroupElementMixinLibGAP


def is_MatrixGroupElement(x):
    """
    Test whether ``x`` is a matrix group element

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.group_element import is_MatrixGroupElement
        sage: is_MatrixGroupElement('helloooo')
        False

        sage: G = GL(2,3)
        sage: is_MatrixGroupElement(G.an_element())
        True
    """
    return isinstance(x, MatrixGroupElement_base)


class MatrixGroupElement_base(MultiplicativeGroupElement):
    """
    Base class for elements of matrix groups.

    You should use one of the two subclasses:

    * :class:`MatrixGroupElement_sage` implements the group
      multiplication using Sage matrices.

    * :class:`MatrixGroupElement_gap` implements the group
      multiplication using libGAP matrices.

    The base class only assumes that derived classes implement
    :meth:`matrix`.

    EXAMPLES::

        sage: F = GF(3); MS = MatrixSpace(F,2,2)
               sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
        sage: G = MatrixGroup(gens)
        sage: g = G.random_element()
        sage: type(g)
        <class 'sage.groups.matrix_gps.group_element.FinitelyGeneratedMatrixGroup_gap_with_category.element_class'>
    """

    def _repr_(self):
        """
        Return string representation of this matrix.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: g  # indirect doctest
            [1 1]
            [0 1]
            sage: g._repr_()
            '[1 1]\n[0 1]'
        """
        return str(self.matrix())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: print g._latex_()
            \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right)

        Type ``view(g._latex_())`` to see the object in an
        xdvi window (assuming you have latex and xdvi installed).
        """
        return self.matrix()._latex_()

    def _act_on_(self, x, self_on_left):
        """
        EXAMPLES::

            sage: G = GL(4,7)
            sage: G.0 * vector([1,2,3,4])
            (3, 2, 3, 4)
            sage: v = vector(GF(7), [3,2,1,-1])
            sage: g = G.1
            sage: v * g == v * g.matrix()   # indirect doctest
            True
        """
        if not is_MatrixGroupElement(x) and x not in self.parent().base_ring():
            try:
                if self_on_left:
                    return self.matrix() * x
                else:
                    return x * self.matrix()
            except TypeError:
                return None

    def _cmp_(self, other):
        """
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = G([1,1, 0,1])
            sage: g == h
            True
            sage: g == G.one()
            False
        """
        return cmp(self.matrix(), other.matrix())

    __cmp__ = _cmp_

    def list(self):
        """
        Return list representation of this matrix.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G.0
            sage: g
            [1 0]
            [0 1]
            sage: g.list()
            [[1, 0], [0, 1]]
        """
        return map(list, self.matrix().rows())


###################################################################
#
# Matrix groups elements implemented with Sage matrices
#
###################################################################

class MatrixGroupElement_generic(MatrixGroupElement_base):

    def __init__(self, parent, M, check=True, convert=True):
        r"""
        Element of a matrix group over a generic ring.

        The group elements are implemented as Sage matrices.

        INPUT:

        - ``M`` -- a matrix.

        - ``parent`` -- the parent.

        - ``check`` -- bool (default: ``True``). If true does some
           type checking.

        - ``convert`` -- bool (default: ``True``). If true convert
          ``M`` to the right matrix space.

        TESTS::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G.random_element()
            sage: TestSuite(g).run()
        """
        if convert:
            M = parent.matrix_space()(M)
        if check:
            if not is_Matrix(M):
                raise TypeError('M must be a matrix')
            if M.parent() is not parent.matrix_space():
                raise TypeError('M must be a in the matrix space of the group')
            parent._check_matrix(M)
        super(MatrixGroupElement_generic, self).__init__(parent)
        if M.is_immutable():
            self._matrix = M
        else:
            self._matrix = M.__copy__()
            self._matrix.set_immutable()

    def matrix(self):
        """
        Obtain the usual matrix (as an element of a matrix space)
        associated to this matrix group element.

        One reason to compute the associated matrix is that matrices
        support a huge range of functionality.

        EXAMPLES::

            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])])
            sage: g = G.0
            sage: g.matrix()
            [1 1]
            [0 1]
            sage: parent(g.matrix())
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7

        Matrices have extra functionality that matrix group elements
        do not have::

            sage: g.matrix().charpoly('t')
            t^2 + 5*t + 1
        """
        return self._matrix

    def _mul_(self,other):
        """
        Return the product of self and other, which must have identical
        parents.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = G([1,1, 0,1])
            sage: g*h         # indirect doctest
            [1 2]
            [0 1]
        """
        parent = self.parent()
        return parent.element_class(parent, self._matrix * other._matrix,
                                    check=False, convert=False)

    def __invert__(self):
        """
        Return the inverse group element

        OUTPUT:

        A matrix group element.

        EXAMPLES::

           sage: G = GL(2,3)
           sage: g = G([1,2,1,0]); g
           [1 2]
           [1 0]
           sage: g.__invert__()
           [0 1]
           [2 1]
           sage: g * (~g)
           [1 0]
           [0 1]
        """
        parent = self.parent()
        return parent.element_class(parent, ~self._matrix, check=False, convert=False)

    inverse = __invert__

###################################################################
#
# Matrix group elements implemented in GAP
#
###################################################################

class MatrixGroupElement_gap(GroupElementMixinLibGAP, MatrixGroupElement_base, ElementLibGAP):

    def __init__(self, parent, M, check=True, convert=True):
        r"""
        Element of a matrix group over a generic ring.

        The group elements are implemented as Sage matrices.

        INPUT:

        - ``M`` -- a matrix.

        - ``parent`` -- the parent.

        - ``check`` -- bool (default: ``True``). If true does some
          type checking.

        - ``convert`` -- bool (default: ``True``). If true convert
          ``M`` to the right matrix space.

        TESTS::

            sage: MS = MatrixSpace(GF(3),2,2)
            sage: G = MatrixGroup(MS([[1,0],[0,1]]), MS([[1,1],[0,1]]))
            sage: G.gen(0)
            [1 0]
            [0 1]
            sage: g = G.random_element()
            sage: TestSuite(g).run()
        """
        if isinstance(M, GapElement):
            ElementLibGAP.__init__(self, parent, M)
            return
        if convert:
            M = parent.matrix_space()(M)
        from sage.libs.gap.libgap import libgap
        M_gap = libgap(M)
        if check:
            if not is_Matrix(M):
                raise TypeError('M must be a matrix')
            if M.parent() is not parent.matrix_space():
                raise TypeError('M must be a in the matrix space of the group')
            parent._check_matrix(M, M_gap)
        ElementLibGAP.__init__(self, parent, M_gap)

    def __reduce__(self):
        """
        Implement pickling.

        TESTS::

            sage: MS = MatrixSpace(GF(3), 2, 2)
            sage: G = MatrixGroup(MS([[1,0],[0,1]]), MS([[1,1],[0,1]]))
            sage: loads(G.gen(0).dumps())
            [1 0]
            [0 1]
        """
        return (self.parent(), (self.matrix(),))

    @cached_method
    def matrix(self):
        """
        Obtain the usual matrix (as an element of a matrix space)
        associated to this matrix group element.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.gen(0).matrix()
            [1 0]
            [0 1]
            sage: _.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 3
        """
        g = self.gap()
        m = g.matrix(self.base_ring())
        m.set_immutable()
        return m

