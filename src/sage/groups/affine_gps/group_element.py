# sage_setup: distribution = sagemath-modules
r"""
Elements of Affine Groups

The class in this module is used to represent the elements of
:func:`~sage.groups.affine_gps.affine_group.AffineGroup` and its
subgroups.

EXAMPLES::

    sage: F = AffineGroup(3, QQ)
    sage: F([1,2,3,4,5,6,7,8,0], [10,11,12])
          [1 2 3]     [10]
    x |-> [4 5 6] x + [11]
          [7 8 0]     [12]

    sage: G = AffineGroup(2, ZZ)
    sage: g = G([[1,1],[0,1]], [1,0])
    sage: h = G([[1,2],[0,1]], [0,1])
    sage: g*h
          [1 3]     [2]
    x |-> [0 1] x + [1]
    sage: h*g
          [1 3]     [1]
    x |-> [0 1] x + [1]
    sage: g*h != h*g
    True

AUTHORS:

- Volker Braun
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Matrix
from sage.misc.cachefunc import cached_method
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.richcmp import richcmp, richcmp_not_equal


class AffineGroupElement(MultiplicativeGroupElement):
    r"""
    An affine group element.

    INPUT:

    - ``A`` -- an invertible matrix, or something defining a
      matrix if ``convert==True``

    - ``b`` -- a vector, or something defining a vector if
      ``convert==True`` (default: ``0``, defining the zero
      vector)

    - ``parent`` -- the parent affine group

    - ``convert`` -- boolean (default: ``True``); whether to convert
      ``A`` into the correct matrix space and ``b`` into the
      correct vector space

    - ``check`` -- boolean (default: ``True``); whether to do some
      checks or just accept the input as valid

    As a special case, ``A`` can be a matrix obtained from
    :meth:`matrix`, that is, one row and one column larger. In
    that case, the group element defining that matrix is
    reconstructed.

    OUTPUT: the affine group element `x \mapsto Ax + b`

    EXAMPLES::

        sage: G = AffineGroup(2, GF(3))

        sage: # needs sage.libs.gap
        sage: g = G.random_element()
        sage: type(g)
        <class 'sage.groups.affine_gps.affine_group.AffineGroup_with_category.element_class'>
        sage: G(g.matrix()) == g
        True

        sage: G(2)
              [2 0]     [0]
        x |-> [0 2] x + [0]

    Conversion from a matrix and a matrix group element::

        sage: M = Matrix(4, 4, [0, 0, -1, 1, 0, -1, 0, 1, -1, 0, 0, 1, 0, 0, 0, 1])
        sage: A = AffineGroup(3, ZZ)
        sage: A(M)
              [ 0  0 -1]     [1]
        x |-> [ 0 -1  0] x + [1]
              [-1  0  0]     [1]
        sage: G = MatrixGroup([M])
        sage: A(G.0)
              [ 0  0 -1]     [1]
        x |-> [ 0 -1  0] x + [1]
              [-1  0  0]     [1]
    """
    def __init__(self, parent, A, b=0, convert=True, check=True):
        r"""
        Create element of an affine group.

        TESTS::

            sage: # needs sage.libs.gap
            sage: G = AffineGroup(4, GF(5))
            sage: g = G.random_element()
            sage: TestSuite(g).run()
        """
        try:
            A = A.matrix()
        except AttributeError:
            pass
        if isinstance(A, Matrix) and A.nrows() == A.ncols() == parent.degree()+1:
            g = A
            d = parent.degree()
            A = g.submatrix(0, 0, d, d)
            b = [ g[i,d] for i in range(d) ]
            convert = True
        if convert:
            A = parent.matrix_space()(A)
            b = parent.vector_space()(b)
        if check:
            # Note: the coercion framework expects that we raise TypeError for invalid input
            if not isinstance(A, Matrix):
                raise TypeError('A must be a matrix')
            if not (A.parent() is parent.matrix_space()):
                raise TypeError('A must be an element of ' + str(parent.matrix_space()))
            if not (b.parent() is parent.vector_space()):
                raise TypeError('b must be an element of ' + str(parent.vector_space()))
            parent._element_constructor_check(A, b)
        super().__init__(parent)
        self._A = A
        self._b = b

    def A(self):
        """
        Return the general linear part of an affine group element.

        OUTPUT: the matrix `A` of the affine group element `Ax + b`

        EXAMPLES::

            sage: G = AffineGroup(3, QQ)
            sage: g = G([1,2,3,4,5,6,7,8,0], [10,11,12])
            sage: g.A()
            [1 2 3]
            [4 5 6]
            [7 8 0]
        """
        return self._A

    def b(self):
        """
        Return the translation part of an affine group element.

        OUTPUT: the vector `b` of the affine group element `Ax + b`

        EXAMPLES::

            sage: G = AffineGroup(3, QQ)
            sage: g = G([1,2,3,4,5,6,7,8,0], [10,11,12])
            sage: g.b()
            (10, 11, 12)
        """
        return self._b

    @cached_method
    def matrix(self):
        """
        Return the standard matrix representation of ``self``.

        .. SEEALSO::

            - :meth:`AffineGroup.linear_space()`

        EXAMPLES::

            sage: G = AffineGroup(3, GF(7))
            sage: g = G([1,2,3,4,5,6,7,8,0], [10,11,12])
            sage: g
                  [1 2 3]     [3]
            x |-> [4 5 6] x + [4]
                  [0 1 0]     [5]
            sage: g.matrix()
            [1 2 3|3]
            [4 5 6|4]
            [0 1 0|5]
            [-----+-]
            [0 0 0|1]
            sage: parent(g.matrix())
            Full MatrixSpace of 4 by 4 dense matrices over Finite Field of size 7
            sage: g.matrix() == matrix(g)
            True

        Composition of affine group elements equals multiplication of
        the matrices::

            sage: # needs sage.libs.gap
            sage: g1 = G.random_element()
            sage: g2 = G.random_element()
            sage: g1.matrix() * g2.matrix() == (g1*g2).matrix()
            True
        """
        A = self._A
        b = self._b
        parent = self.parent()
        d = parent.degree()
        from sage.matrix.constructor import matrix, zero_matrix, block_matrix
        zero = zero_matrix(parent.base_ring(), 1, d)
        one = matrix(parent.base_ring(), [[1]])
        m = block_matrix(2,2, [A, b.column(), zero, one])
        m.set_immutable()
        return m

    _matrix_ = matrix

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = AffineGroup(2, QQ)
            sage: g = G([[1, 1], [0, 1]], [3,4])
            sage: g
                  [1 1]     [3]
            x |-> [0 1] x + [4]
        """
        A = str(self._A)
        b = str(self._b.column())
        deg = self.parent().degree()
        indices = range(deg)
        s = []
        for Ai, bi, i in zip(A.splitlines(), b.splitlines(), indices):
            if i == deg//2:
                s.append('x |-> '+Ai+' x + '+bi)
            else:
                s.append('      '+Ai+'     '+bi)
        return '\n'.join(s)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: G = AffineGroup(2, QQ)
            sage: g = G([[1, 1], [0, 1]], [3,4])
            sage: latex(g)
            \vec{x}\mapsto \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right)\vec{x} + \left(\begin{array}{r}
            3 \\
            4
            \end{array}\right)
            sage: g._latex_()
            '\\vec{x}\\mapsto \\left(\\begin{array}{rr}\n1 & 1 \\\\\n0 &
            1\n\\end{array}\\right)\\vec{x} + \\left(\\begin{array}{r}\n3
            \\\\\n4\n\\end{array}\\right)'
        """
        return r'\vec{x}\mapsto '+self.A()._latex_()+r'\vec{x} + '+self.b().column()._latex_()

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: G2 = AffineGroup(2, QQ)
            sage: g2 = G2([[1, 1], [0, 1]], [3,4])
            sage: ascii_art(g2)
            x |-> [1 1] x + [3]
                  [0 1]     [4]

            sage: G3 = AffineGroup(3, QQ)
            sage: g3 = G3([[1,1,-1], [0,1,2], [0,10,2]], [3,4,5/2])
            sage: ascii_art(g3)
                  [ 1  1 -1]     [  3]
            x |-> [ 0  1  2] x + [  4]
                  [ 0 10  2]     [5/2]
        """
        from sage.typeset.ascii_art import ascii_art
        deg = self.parent().degree()
        A = ascii_art(self._A, baseline=deg//2)
        b = ascii_art(self._b.column(), baseline=deg//2)
        return ascii_art("x |-> ") + A + ascii_art(" x + ") + b

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: G2 = AffineGroup(2, QQ)
            sage: g2 = G2([[1, 1], [0, 1]], [3,4])
            sage: unicode_art(g2)
            x ↦ ⎛1 1⎞ x + ⎛3⎞
                ⎝0 1⎠     ⎝4⎠

            sage: G3 = AffineGroup(3, QQ)
            sage: g3 = G3([[1,1,-1], [0,1,2], [0,10,2]], [3,4,5/2])
            sage: unicode_art(g3)
                ⎛ 1  1 -1⎞     ⎛  3⎞
            x ↦ ⎜ 0  1  2⎟ x + ⎜  4⎟
                ⎝ 0 10  2⎠     ⎝5/2⎠
        """
        from sage.typeset.unicode_art import unicode_art
        deg = self.parent().degree()
        A = unicode_art(self._A, baseline=deg//2)
        b = unicode_art(self._b.column(), baseline=deg//2)
        return unicode_art("x ↦ ") + A + unicode_art(" x + ") + b

    def _mul_(self, other):
        """
        Return the composition of ``self`` and ``other``.

        INPUT:

        - ``other`` -- another element of the same affine group

        OUTPUT: the product of the affine group elements ``self`` and
        ``other`` defined by the composition of the two affine transformations

        EXAMPLES::

            sage: G = AffineGroup(2, GF(3))
            sage: g = G([1,1, 0,1], [0,1])
            sage: h = G([1,1, 0,1], [1,2])
            sage: g*h
                  [1 2]     [0]
            x |-> [0 1] x + [0]
            sage: g.matrix() * h.matrix() == (g*h).matrix()
            True
        """
        parent = self.parent()
        A = self._A * other._A
        b = self._b + self._A * other._b
        return parent.element_class(parent, A, b, check=False)

    def __call__(self, v):
        """
        Apply the affine transformation to ``v``.

        INPUT:

        - ``v`` -- a polynomial, a multivariate polynomial, a polyhedron, a
          vector, or anything that can be converted into a vector

        OUTPUT: the image of ``v`` under the affine group element

        EXAMPLES::

            sage: G = AffineGroup(2, QQ)
            sage: g = G([0,1,-1,0],[2,3]);  g
                  [ 0  1]     [2]
            x |-> [-1  0] x + [3]
            sage: v = vector([4,5])
            sage: g(v)
            (7, -1)

            sage: R.<x,y> = QQ[]
            sage: g(x), g(y)
            (y + 2, -x + 3)
            sage: p = x^2 + 2*x*y + y + 1
            sage: g(p)
            -2*x*y + y^2 - 5*x + 10*y + 20

        The action on polynomials is such that it intertwines with
        evaluation. That is::

            sage: p(*g(v)) == g(p)(*v)
            True

        Test that the univariate polynomial ring is covered::

            sage: H = AffineGroup(1, QQ)
            sage: h = H([2],[3]);  h
            x |-> [2] x + [3]
            sage: R.<z> = QQ[]
            sage: h(z+1)
            3*z + 2

        The action on a polyhedron is defined (see :issue:`30327`)::

            sage: F = AffineGroup(3, QQ)
            sage: M = matrix(3, [-1, -2, 0, 0, 0, 1, -2, 1, -1])
            sage: v = vector(QQ,(1,2,3))
            sage: f = F(M, v)
            sage: cube = polytopes.cube()                                               # needs sage.geometry.polyhedron
            sage: f(cube)                                                               # needs sage.geometry.polyhedron
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices
        """
        parent = self.parent()

        # start with the most probable case, i.e., v is in the vector space
        if v in parent.vector_space():
            return self._A*v + self._b

        from sage.rings.polynomial.polynomial_element import Polynomial
        if isinstance(v, Polynomial) and parent.degree() == 1:
            ring = v.parent()
            return ring([self._A[0,0], self._b[0]])

        from sage.rings.polynomial.multi_polynomial import MPolynomial
        if isinstance(v, MPolynomial) and parent.degree() == v.parent().ngens():
            ring = v.parent()
            from sage.modules.free_module_element import vector
            image_coords = self._A * vector(ring, ring.gens()) + self._b
            return v(*image_coords)

        import sage.geometry.abc
        if isinstance(v, sage.geometry.abc.Polyhedron):
            return self._A*v + self._b

        # otherwise, coerce v into the vector space
        v = parent.vector_space()(v)
        return self._A*v + self._b

    def _act_on_(self, x, self_on_left):
        """
        Define the multiplicative action of the affine group elements.

        EXAMPLES::

            sage: G = AffineGroup(2, GF(3))
            sage: g = G([1,2,3,4], [5,6])
            sage: g
                  [1 2]     [2]
            x |-> [0 1] x + [0]
            sage: v = vector(GF(3), [1,-1]); v
            (1, 2)
            sage: g*v
            (1, 2)
            sage: g*v == g.A() * v + g.b()
            True
        """
        if self_on_left:
            return self(x)

    def __invert__(self):
        """
        Return the inverse group element.

        OUTPUT: another affine group element

        EXAMPLES::

            sage: G = AffineGroup(2, GF(3))
            sage: g = G([1,2,3,4], [5,6])
            sage: g
                  [1 2]     [2]
            x |-> [0 1] x + [0]
            sage: ~g
                  [1 1]     [1]
            x |-> [0 1] x + [0]
            sage: g * g.inverse()   # indirect doctest
                  [1 0]     [0]
            x |-> [0 1] x + [0]
            sage: g * g.inverse() == g.inverse() * g == G(1)
            True
        """
        parent = self.parent()
        A = parent.matrix_space()(~self._A)
        b = -A * self.b()
        return parent.element_class(parent, A, b, check=False)

    def _richcmp_(self, other, op):
        """
        Compare ``self`` with ``other``.

        OUTPUT: boolean

        EXAMPLES::

            sage: F = AffineGroup(3, QQ)
            sage: g = F([1,2,3,4,5,6,7,8,0], [10,11,12])
            sage: h = F([1,2,3,4,5,6,7,8,0], [10,11,0])
            sage: g == h
            False
            sage: g == g
            True
        """
        lx = self._A
        rx = other._A
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        return richcmp(self._b, other._b, op)

    def list(self):
        """
        Return list representation of ``self``.

        EXAMPLES::

            sage: F = AffineGroup(3, QQ)
            sage: g = F([1,2,3,4,5,6,7,8,0], [10,11,12])
            sage: g
                  [1 2 3]     [10]
            x |-> [4 5 6] x + [11]
                  [7 8 0]     [12]
            sage: g.matrix()
            [ 1  2  3|10]
            [ 4  5  6|11]
            [ 7  8  0|12]
            [--------+--]
            [ 0  0  0| 1]
            sage: g.list()
            [[1, 2, 3, 10], [4, 5, 6, 11], [7, 8, 0, 12], [0, 0, 0, 1]]
        """
        return [r.list() for r in self.matrix().rows()]
