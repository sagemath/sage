"""
Affine Subspaces of a Vector Space

An affine subspace of a vector space is a translation of a linear
subspace. The affine subspaces here are only used internally in
hyperplane arrangements. You should not use them for interactive work
or return them to the user.

EXAMPLES::

    sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
    sage: a = AffineSubspace([1,0,0,0], QQ^4)
    sage: a.dimension()
    4
    sage: a.point()
    (1, 0, 0, 0)
    sage: a.linear_part()
    Vector space of dimension 4 over Rational Field
    sage: a
    Affine space p + W where:
      p = (1, 0, 0, 0)
      W = Vector space of dimension 4 over Rational Field
    sage: b = AffineSubspace((1,0,0,0), matrix(QQ, [[1,2,3,4]]).right_kernel())
    sage: c = AffineSubspace((0,2,0,0), matrix(QQ, [[0,0,1,2]]).right_kernel())
    sage: b.intersection(c)
    Affine space p + W where:
      p = (-3, 2, 0, 0)
      W = Vector space of degree 4 and dimension 2 over Rational Field
    Basis matrix:
    [  1   0  -1 1/2]
    [  0   1  -2   1]
    sage: b < a
    True
    sage: c < b
    False
    sage: A = AffineSubspace([8,38,21,250], VectorSpace(GF(19),4))
    sage: A
    Affine space p + W where:
       p = (8, 0, 2, 3)
       W = Vector space of dimension 4 over Finite Field of size 19

TESTS::

    sage: A = AffineSubspace([2], VectorSpace(QQ, 1))
    sage: A.point()
    (2)
    sage: A.linear_part()
    Vector space of dimension 1 over Rational Field
    sage: A.linear_part().basis_matrix()
    [1]
    sage: A = AffineSubspace([], VectorSpace(QQ, 0))
    sage: A.point()
    ()
    sage: A.linear_part()
    Vector space of dimension 0 over Rational Field
    sage: A.linear_part().basis_matrix()
    []
"""

#*****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.matrix.constructor import vector


class AffineSubspace(SageObject):
    """
    An affine subspace.

    INPUT:

    - ``p`` -- list/tuple/iterable representing a point on the
      affine space

    - ``V`` -- vector subspace

    OUTPUT:

    Affine subspace parallel to ``V`` and passing through ``p``.

    EXAMPLES::

        sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
        sage: a = AffineSubspace([1,0,0,0], VectorSpace(QQ,4))
        sage: a
        Affine space p + W where:
          p = (1, 0, 0, 0)
          W = Vector space of dimension 4 over Rational Field
    """
    def __init__(self, p, V):
        r"""
        Construct an :class:`AffineSubspace`.

        TESTS::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: a = AffineSubspace([1,0,0,0], VectorSpace(QQ,4))
            sage: TestSuite(a).run()
            sage: AffineSubspace(0, VectorSpace(QQ,4)).point()
            (0, 0, 0, 0)
        """
        R = V.base_ring()
        from sage.categories.all import Fields
        if R not in Fields():
            R = R.fraction_field()
            V = V.change_ring(R)
        self._base_ring = R
        self._linear_part = V
        p = V.ambient_vector_space()(p)
        p.set_immutable()
        self._point = p

    def __hash__(self):
        """
        Return a hash value.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: a = AffineSubspace([1,0,0,0], VectorSpace(QQ,4))
            sage: a.__hash__()    # random output
            -3713096828371451969
        """
        # note that the point is not canonically chosen, but the linear part is
        return hash(self._linear_part)
    
    def _repr_(self):
        r"""
        String representation for an :class:`AffineSubspace`.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
            sage: a
            Affine space p + W where:
              p = (1, 0, 0, 0)
              W = Vector space of dimension 4 over Rational Field
        """
        return "Affine space p + W where:\n  p = "+str(self._point)+"\n  W = "+str(self._linear_part)

    def __eq__(self, other):
        r"""
        Test whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- an :class:`AffineSubspace`

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: a = AffineSubspace([1,0,0], matrix([[1,0,0]]).right_kernel())
            sage: b = AffineSubspace([2,0,0], matrix([[1,0,0]]).right_kernel())
            sage: c = AffineSubspace([1,1,0], matrix([[1,0,0]]).right_kernel())
            sage: a == b
            False
            sage: a == c
            True
        """
        V = self._linear_part
        W = other._linear_part
        return V == W and self._point - other._point in V

    def __ne__(self, other):
        r"""
        Test whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- an :class:`AffineSubspace`

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: a = AffineSubspace([1,0,0],matrix([[1,0,0]]).right_kernel())
            sage: b = AffineSubspace([2,0,0],matrix([[1,0,0]]).right_kernel())
            sage: a == b
            False
            sage: a != b
            True
            sage: a != a
            False
        """
        return not self == other

    def __le__(self, other):
        r"""
        Test whether ``self`` is an affine subspace of ``other``.

        INPUT:

        - ``other`` -- an :class:`AffineSubspace`

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: V = VectorSpace(QQ, 3)
            sage: W1 = V.subspace([[1,0,0],[0,1,0]])
            sage: W2 = V.subspace([[1,0,0]])
            sage: a = AffineSubspace([1,2,3], W1)
            sage: b = AffineSubspace([1,2,3], W2)
            sage: a <= b
            False
            sage: a <= a
            True
            sage: b <= a
            True
        """
        V = self._linear_part
        W = other._linear_part
        return V.is_subspace(W) and self._point-other._point in W

    def __lt__(self, other):
        r"""
        Test whether ``self`` is a proper affine subspace of ``other``.

        INPUT:

        - ``other`` -- an :class:`AffineSubspace`

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: V = VectorSpace(QQ, 3)
            sage: W1 = V.subspace([[1,0,0], [0,1,0]])
            sage: W2 = V.subspace([[1,0,0]])
            sage: a = AffineSubspace([1,2,3], W1)
            sage: b = AffineSubspace([1,2,3], W2)
            sage: a < b
            False
            sage: a < a
            False
            sage: b < a
            True
        """
        if self._linear_part == other._linear_part:
            return False
        return self <= other

    def __contains__(self, q):
        r"""
        Test whether the point ``q`` is in the affine space.

        INPUT:

        - ``q`` -- point as a list/tuple/iterable

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: a = AffineSubspace([1,0,0], matrix([[1,0,0]]).right_kernel())
            sage: (1,1,0) in a
            True
            sage: (0,0,0) in a
            False
        """
        q = vector(self._base_ring, q)
        return self._point - q in self._linear_part

    def linear_part(self):
        r"""
        Return the linear part of the affine space.

        OUTPUT:

        A vector subspace of the ambient space.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: A = AffineSubspace([2,3,1], matrix(QQ, [[1,2,3]]).right_kernel())
            sage: A.linear_part()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 -1/3]
            [   0    1 -2/3]
            sage: A.linear_part().ambient_vector_space()
            Vector space of dimension 3 over Rational Field
        """
        return self._linear_part

    def point(self):
        r"""
        Return a point ``p`` in the affine space.

        OUTPUT:

        A point of the affine space as a vector in the ambient space.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: A = AffineSubspace([2,3,1], VectorSpace(QQ,3))
            sage: A.point()
            (2, 3, 1)
        """
        return self._point

    def dimension(self):
        r"""
        Return the dimension of the affine space.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
            sage: a.dimension()
            4
        """
        return self.linear_part().dimension()

    def intersection(self, other):
        r"""
        Return the intersection of ``self`` with ``other``.

        INPUT:

        - ``other`` -- an :class:`AffineSubspace`

        OUTPUT:

        A new affine subspace, (or ``None`` if the intersection is
        empty).

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
            sage: V = VectorSpace(QQ,3)
            sage: U = V.subspace([(1,0,0), (0,1,0)])
            sage: W = V.subspace([(0,1,0), (0,0,1)])
            sage: A = AffineSubspace((0,0,0), U)
            sage: B = AffineSubspace((1,1,1), W)
            sage: A.intersection(B)
            Affine space p + W where:
              p = (1, 1, 0)
              W = Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            sage: C = AffineSubspace((0,0,1), U)
            sage: A.intersection(C)
            sage: C = AffineSubspace((7,8,9), U.complement())
            sage: A.intersection(C)
            Affine space p + W where:
              p = (7, 8, 0)
              W = Vector space of degree 3 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: A.intersection(C).intersection(B)

            sage: D = AffineSubspace([1,2,3], VectorSpace(GF(5),3))
            sage: E = AffineSubspace([3,4,5], VectorSpace(GF(5),3))
            sage: D.intersection(E)
            Affine space p + W where:
              p = (3, 4, 0)
              W = Vector space of dimension 3 over Finite Field of size 5
        """
        if self.linear_part().ambient_vector_space() != \
           other.linear_part().ambient_vector_space():
            raise ValueError('incompatible ambient vector spaces')
        m = self.linear_part().matrix()
        n = other.linear_part().matrix()
        p = self.point()
        q = other.point()
        M = m.stack(n)
        v = q - p
        try:
            t = M.solve_left(v)
        except ValueError:
            return None  # empty intersection
        new_p = p + t[:m.nrows()]*m
        new_V = self.linear_part().intersection(other._linear_part)
        return AffineSubspace(new_p, new_V)

