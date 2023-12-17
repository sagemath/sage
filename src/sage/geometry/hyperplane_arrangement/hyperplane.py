r"""
Hyperplanes

.. NOTE::

    If you want to learn about Sage's hyperplane arrangements then you
    should start with
    :mod:`sage.geometry.hyperplane_arrangement.arrangement`. This
    module is used to represent the individual hyperplanes, but you
    should never construct the classes from this module directly (but
    only via the
    :class:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangements`.

A linear expression, for example, `3x+3y-5z-7` stands for the
hyperplane with the equation `x+3y-5z=7`. To create it in Sage, you
first have to create a
:class:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangements`
object to define the variables `x`, `y`, `z`::

    sage: H.<x,y,z> = HyperplaneArrangements(QQ)
    sage: h = 3*x + 2*y - 5*z - 7;  h
    Hyperplane 3*x + 2*y - 5*z - 7
    sage: h.coefficients()
    [-7, 3, 2, -5]
    sage: h.normal()
    (3, 2, -5)
    sage: h.constant_term()
    -7
    sage: h.change_ring(GF(3))
    Hyperplane 0*x + 2*y + z + 2
    sage: h.point()
    (21/38, 7/19, -35/38)
    sage: h.linear_part()
    Vector space of degree 3 and dimension 2 over Rational Field
    Basis matrix:
    [  1   0 3/5]
    [  0   1 2/5]

Another syntax to create hyperplanes is to specify coefficients and a
constant term::

    sage: V = H.ambient_space();  V
    3-dimensional linear space over Rational Field with coordinates x, y, z
    sage: h in V
    True
    sage: V([3, 2, -5], -7)
    Hyperplane 3*x + 2*y - 5*z - 7

Or constant term and coefficients together in one list/tuple/iterable::

    sage: V([-7, 3, 2, -5])
    Hyperplane 3*x + 2*y - 5*z - 7
    sage: v = vector([-7, 3, 2, -5]);  v
    (-7, 3, 2, -5)
    sage: V(v)
    Hyperplane 3*x + 2*y - 5*z - 7

Note that the constant term comes first, which matches the notation
for Sage's :func:`~sage.geometry.polyhedron.constructor.Polyhedron` ::

    sage: Polyhedron(ieqs=[(4,1,2,3)]).Hrepresentation()
    (An inequality (1, 2, 3) x + 4 >= 0,)

The difference between hyperplanes as implemented in this module and
hyperplane arrangements is that:

* hyperplane arrangements contain multiple hyperplanes (of course),

* linear expressions are a module over the base ring, and these module
  structure is inherited by the hyperplanes.

The latter means that you can add and multiply by a scalar::

    sage: h = 3*x + 2*y - 5*z - 7;  h
    Hyperplane 3*x + 2*y - 5*z - 7
    sage: -h
    Hyperplane -3*x - 2*y + 5*z + 7
    sage: h + x
    Hyperplane 4*x + 2*y - 5*z - 7
    sage: h + 7
    Hyperplane 3*x + 2*y - 5*z + 0
    sage: 3*h
    Hyperplane 9*x + 6*y - 15*z - 21
    sage: h * RDF(3)
    Hyperplane 9.0*x + 6.0*y - 15.0*z - 21.0

Which you can't do with hyperplane arrangements::

    sage: arrangement = H(h, x, y, x+y-1);  arrangement
    Arrangement <y | x | x + y - 1 | 3*x + 2*y - 5*z - 7>
    sage: arrangement + x
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +:
    'Hyperplane arrangements in 3-dimensional linear space
         over Rational Field with coordinates x, y, z' and
    'Hyperplane arrangements in 3-dimensional linear space
         over Rational Field with coordinates x, y, z'
"""

# *****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.geometry.linear_expression import LinearExpression, LinearExpressionModule


class Hyperplane(LinearExpression):
    """
    A hyperplane.

    You should always use :class:`AmbientVectorSpace` to construct
    instances of this class.

    INPUT:

    - ``parent`` -- the parent :class:`AmbientVectorSpace`

    - ``coefficients`` -- a vector of coefficients of the linear variables

    - ``constant`` -- the constant term for the linear expression

    EXAMPLES::

        sage: H.<x,y> = HyperplaneArrangements(QQ)
        sage: x+y-1
        Hyperplane x + y - 1

        sage: ambient = H.ambient_space()
        sage: ambient._element_constructor_(x+y-1)
        Hyperplane x + y - 1

    For technical reasons, we must allow the degenerate cases of
    an empty space and of a full space::

        sage: 0*x
        Hyperplane 0*x + 0*y + 0
        sage: 0*x + 1
        Hyperplane 0*x + 0*y + 1
        sage: x + 0 == x + ambient(0)    # because coercion requires them
        True
    """
    def __init__(self, parent, coefficients, constant):
        """
        Initialize ``self``.

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: x.change_ring(RR)
            Hyperplane 1.00000000000000*x + 0.000000000000000*y + 0.000000000000000
            sage: TestSuite(x+y-1).run()
        """
        super().__init__(parent, coefficients, constant)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: H.<x> = HyperplaneArrangements(QQ)
            sage: x._repr_()
            'Hyperplane x + 0'
        """
        return 'Hyperplane {0}'.format(self._repr_linear())

    def _latex_(self):
        r"""
        Return a LaTeX representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: H.<x> = HyperplaneArrangements(QQ)
            sage: V = H.ambient_space()
            sage: V([2, -3])._latex_()
            '$-3x = -2$'

            sage: H.<x, y, z> = HyperplaneArrangements(QQ)
            sage: V = H.ambient_space()
            sage: V([-5, 1, 3, 0])._latex_()
            '$x + 3y = 5$'
            sage: V([4, 1, 0, -1])._latex_()
            '$x - z = -4$'
        """
        linear = self._repr_linear(include_zero=False, include_constant=False, multiplication='')
        s = '{0} = {1}'.format(linear, -self.b())
        return '${0}$'.format(s)

    def normal(self):
        """
        Return the normal vector.

        OUTPUT:

        A vector over the base ring.

        EXAMPLES::

            sage: H.<x, y, z> = HyperplaneArrangements(QQ)
            sage: x.normal()
            (1, 0, 0)
            sage: x.A(), x.b()
            ((1, 0, 0), 0)
            sage: (x + 2*y + 3*z + 4).normal()
            (1, 2, 3)
        """
        return self.A()

    def _normal_pivot(self):
        """
        Return the index of the largest entry of the normal vector.

        OUTPUT:

        An integer. The index of the largest entry.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: V = H.ambient_space()
            sage: (x + 3/2*y - 2*z)._normal_pivot()
            2

            sage: H.<x,y,z> = HyperplaneArrangements(GF(5))
            sage: V = H.ambient_space()
            sage: (x + 3*y - 4*z)._normal_pivot()
            1
        """
        try:
            values = [abs(x) for x in self.A()]
        except ArithmeticError:
            from sage.rings.real_double import RDF
            values = [abs(RDF(x)) for x in self.A()]
        max_pos = 0
        max_value = values[max_pos]
        for i in range(1, len(values)):
            if values[i] > max_value:
                max_pos = i
                max_value = values[i]
        return max_pos

    def __contains__(self, q):
        r"""
        Test whether the point ``q`` is in the hyperplane.

        INPUT:

        - ``q`` -- point (as a vector, list, or tuple)

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + y + z - 1
            sage: (1/3, 1/3, 1/3) in h
            True
            sage: (0,0,0) in h
            False
        """
        V = self.parent().ambient_vector_space()
        q = V(q)
        return self.A() * q + self._const == 0

    @cached_method
    def polyhedron(self, **kwds):
        """
        Return the hyperplane as a polyhedron.

        OUTPUT:

        A :func:`~sage.geometry.polyhedron.constructor.Polyhedron` instance.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + 2*y + 3*z - 4
            sage: P = h.polyhedron();  P
            A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex and 2 lines
            sage: P.Hrepresentation()
            (An equation (1, 2, 3) x - 4 == 0,)
            sage: P.Vrepresentation()
            (A line in the direction (0, 3, -2),
             A line in the direction (3, 0, -1),
             A vertex at (0, 0, 4/3))
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        R = kwds.pop('base_ring', None)
        if R is None:
            R = self.parent().base_ring()
        return Polyhedron(eqns=[self.coefficients()], base_ring=R, **kwds)

    @cached_method
    def linear_part(self):
        r"""
        The linear part of the affine space.

        OUTPUT:

        Vector subspace of the ambient vector space, parallel to the
        hyperplane.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + 2*y + 3*z - 1
            sage: h.linear_part()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 -1/3]
            [   0    1 -2/3]
        """
        AA = self.parent().ambient_module()
        from sage.matrix.constructor import matrix
        return matrix(AA.base_ring(), [self.A()]).right_kernel()

    def linear_part_projection(self, point):
        """
        Orthogonal projection onto the linear part.

        INPUT:

        - ``point`` -- vector of the ambient space, or anything that
          can be converted into one; not necessarily on the
          hyperplane

        OUTPUT:

        Coordinate vector of the projection of ``point`` with respect
        to the basis of :meth:`linear_part`. In particular, the length
        of this vector is one less than the ambient space
        dimension.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + 2*y + 3*z - 4
            sage: h.linear_part()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 -1/3]
            [   0    1 -2/3]
            sage: p1 = h.linear_part_projection(0);  p1
            (0, 0)
            sage: p2 = h.linear_part_projection([3,4,5]);  p2
            (8/7, 2/7)
            sage: h.linear_part().basis()
            [
            (1, 0, -1/3),
            (0, 1, -2/3)
            ]
            sage: p3 = h.linear_part_projection([1,1,1]);  p3
            (4/7, 1/7)
        """
        point = self.orthogonal_projection(point) - self.point()
        return self.linear_part().coordinate_vector(point)

    @cached_method
    def point(self):
        """
        Return the point closest to the origin.

        OUTPUT:

        A vector of the ambient vector space. The closest point to the
        origin in the `L^2`-norm.

        In finite characteristic a random point will be returned if
        the norm of the hyperplane normal vector is zero.

        EXAMPLES::


            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + 2*y + 3*z - 4
            sage: h.point()
            (2/7, 4/7, 6/7)
            sage: h.point() in h
            True

            sage: # needs sage.rings.finite_rings
            sage: H.<x,y,z> = HyperplaneArrangements(GF(3))
            sage: h = 2*x + y + z + 1
            sage: h.point()
            (1, 0, 0)
            sage: h.point().base_ring()
            Finite Field of size 3

            sage: H.<x,y,z> = HyperplaneArrangements(GF(3))
            sage: h = x + y + z + 1
            sage: h.point()
            (2, 0, 0)
        """
        P = self.parent()
        AA = P.ambient_module()
        R = P.base_ring()
        norm2 = sum(x**2 for x in self.A())
        if norm2 == 0:
            from sage.matrix.constructor import matrix, vector
            solution = matrix(R, self.A()).solve_right(vector(R, [-self.b()]))
        else:
            solution = [-x * self.b() / norm2 for x in self.A()]
        return AA(solution)

    def dimension(self):
        r"""
        The dimension of the hyperplane.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + y + z - 1
            sage: h.dimension()
            2
        """
        return self.linear_part().dimension()

    def intersection(self, other):
        r"""
        The intersection of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a hyperplane, a polyhedron, or something that
          defines a polyhedron

        OUTPUT:

        A polyhedron.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + y + z - 1
            sage: h.intersection(x - y)
            A 1-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex and 1 line
            sage: h.intersection(polytopes.cube())
            A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices
        """
        from sage.geometry.polyhedron.base import is_Polyhedron
        from sage.geometry.polyhedron.constructor import Polyhedron
        if not is_Polyhedron(other):
            try:
                other = other.polyhedron()
            except AttributeError:
                other = Polyhedron(other)
        return self.polyhedron().intersection(other)

    def orthogonal_projection(self, point):
        """
        Return the orthogonal projection of a point.

        INPUT:

        - ``point`` -- vector of the ambient space, or anything that
          can be converted into one; not necessarily on the
          hyperplane

        OUTPUT:

        A vector in the ambient vector space that lies on the
        hyperplane.

        In finite characteristic, a ``ValueError`` is raised if the
        the norm of the hyperplane normal is zero.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: h = x + 2*y + 3*z - 4
            sage: p1 = h.orthogonal_projection(0);  p1
            (2/7, 4/7, 6/7)
            sage: p1 in h
            True
            sage: p2 = h.orthogonal_projection([3,4,5]);  p2
            (10/7, 6/7, 2/7)
            sage: p1 in h
            True
            sage: p3 = h.orthogonal_projection([1,1,1]);  p3
            (6/7, 5/7, 4/7)
            sage: p3 in h
            True
        """
        P = self.parent()
        norm2 = sum(x**2 for x in self.A())
        if norm2 == 0:
            raise ValueError('norm of hyperplane normal is zero')
        point = P.ambient_vector_space()(point)
        n = self.normal()
        return point - n * (self.b() + point*n) / norm2

    def primitive(self, signed=True):
        """
        Return hyperplane defined by primitive equation.

        INPUT:

        - ``signed`` -- boolean (optional, default: ``True``); whether
          to preserve the overall sign

        OUTPUT:

        Hyperplane whose linear expression has common factors and
        denominators cleared. That is, the same hyperplane (with the
        same sign) but defined by a rescaled equation. Note that
        different linear expressions must define different hyperplanes
        as comparison is used in caching.

        If ``signed``, the overall rescaling is by a positive constant
        only.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = -1/3*x + 1/2*y - 1;  h
            Hyperplane -1/3*x + 1/2*y - 1
            sage: h.primitive()
            Hyperplane -2*x + 3*y - 6
            sage: h == h.primitive()
            False
            sage: (4*x + 8).primitive()
            Hyperplane x + 0*y + 2

            sage: (4*x - y - 8).primitive(signed=True)   # default
            Hyperplane 4*x - y - 8
            sage: (4*x - y - 8).primitive(signed=False)
            Hyperplane -4*x + y + 8

        TESTS:

        Check that :trac:`30078` is fixed::

            sage: # needs sage.rings.number_field
            sage: R.<sqrt2> = QuadraticField(2)
            sage: H.<x,y> = HyperplaneArrangements(base_ring=R)
            sage: B = H([1,1,0], [2,2,0], [sqrt2,sqrt2,0])
            sage: B
            Arrangement <x + 1>

        Check that :trac:`30749` is fixed::

            sage: # needs sage.rings.number_field
            sage: tau = (1+AA(5).sqrt()) / 2
            sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1]]
            sage: ncn += [[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1]]
            sage: ncn += [[1,1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
            sage: H = HyperplaneArrangements(AA,names='xyz')
            sage: A = H([[0]+v for v in ncn])
            sage: A.n_regions()
            60
        """
        from sage.rings.rational_field import QQ
        base_ring = self.parent().base_ring()
        coeffs = self.coefficients()
        # first check if the linear expression even defines a hyperplane
        if self.is_zero():
            raise ValueError('linear expression must be non-constant to define a hyperplane')
        # for scalar adjustment over the base ring QQ,
        # get rid of the denominators and use gcd
        if base_ring is QQ:
            from sage.arith.functions import lcm
            from sage.arith.misc import gcd
            d = lcm(x.denominator() for x in coeffs)
            n = gcd(x.numerator() for x in coeffs)
            adjustment = d/n
        # over other base rings, rescale the coefficients so that
        # the first nonzero of the normal vector is one or negative one
        else:
            for x in coeffs[1:]:
                if not x.is_zero():
                    adjustment = x.inverse_of_unit()
                    if x < 0:  # avoid accidental sign reversal
                        adjustment = -adjustment
                    break
        # if ``signed`` is not set, adjust the sign
        if not signed:
            for x in coeffs:
                if not x.is_zero():
                    if x < 0:
                        adjustment = -adjustment
                    break
        # return the rescaled hyperplane
        return self.parent(adjustment * self)

    @cached_method
    def _affine_subspace(self):
        """
        Return the hyperplane as affine subspace.

        OUTPUT:

        The hyperplane as a
        :class:`~sage.geometry.hyperplane_arrangement.affine_subspace.AffineSubspace`.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = -1/3*x + 1/2*y - 1;  h
            Hyperplane -1/3*x + 1/2*y - 1
            sage: h._affine_subspace()
            Affine space p + W where:
               p = (-12/13, 18/13)
               W = Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [  1 2/3]
        """
        from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
        return AffineSubspace(self.point(), self.linear_part())

    def plot(self, **kwds):
        """
        Plot the hyperplane.

        OUTPUT:

        A graphics object.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: (x + y - 2).plot()                                                    # needs sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        from sage.geometry.hyperplane_arrangement.plot import plot_hyperplane
        return plot_hyperplane(self, **kwds)

    def __or__(self, other):
        """
        Construct hyperplane arrangement from bitwise or.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: x | y + 1
            Arrangement <y + 1 | x>
            sage: x | [(0,1), 1]
            Arrangement <y + 1 | x>

        TESTS::

            sage: (x | y).parent() is L
            True
        """
        from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
        parent = self.parent()
        arrangement = HyperplaneArrangements(parent.base_ring(), names=parent._names)
        return arrangement(self, other)

    def to_symmetric_space(self):
        """
        Return ``self`` considered as an element in the corresponding
        symmetric space.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: h = -1/3*x + 1/2*y
            sage: h.to_symmetric_space()
            -1/3*x + 1/2*y

            sage: hp = -1/3*x + 1/2*y - 1
            sage: hp.to_symmetric_space()
            Traceback (most recent call last):
            ...
            ValueError: the hyperplane must pass through the origin
        """
        coeff = self.coefficients()
        if coeff[0] != 0:
            raise ValueError("the hyperplane must pass through the origin")
        S = self.parent().symmetric_space()
        G = S.gens()
        # We skip the first coefficient since it corresponds to the constant term
        return S.sum(G[i]*c for i, c in enumerate(coeff[1:]))


class AmbientVectorSpace(LinearExpressionModule):
    """
    The ambient space for hyperplanes.

    This class is the parent for the :class:`Hyperplane` instances.

    TESTS::

        sage: from sage.geometry.hyperplane_arrangement.hyperplane import AmbientVectorSpace
        sage: V = AmbientVectorSpace(QQ, ('x', 'y'))
        sage: V.change_ring(QQ) is V
        True
    """

    Element = Hyperplane

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.geometry.hyperplane_arrangement.hyperplane import AmbientVectorSpace
            sage: AmbientVectorSpace(QQ, ('x', 'y'))
            2-dimensional linear space over Rational Field with coordinates x, y
        """
        return '{0}-dimensional linear space over {3} with coordinate{1} {2}'.format(
            self.dimension(),
            's' if self.ngens() > 1 else '',
            ', '.join(self._names),
            self.base_ring())

    def dimension(self):
        """
        Return the ambient space dimension.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: M.<x,y> = HyperplaneArrangements(QQ)
            sage: x.parent().dimension()
            2
            sage: x.parent() is M.ambient_space()
            True
            sage: x.dimension()
            1
        """
        return self.ngens()

    def change_ring(self, base_ring):
        """
        Return a ambient vector space with a changed base ring.

        INPUT:

        - ``base_ring`` -- a ring; the new base ring

        OUTPUT:

        A new :class:`AmbientVectorSpace`.

        EXAMPLES::

            sage: M.<y> = HyperplaneArrangements(QQ)
            sage: V = M.ambient_space()
            sage: V.change_ring(RR)
            1-dimensional linear space over Real Field with 53 bits of precision with coordinate y

        TESTS::

            sage: V.change_ring(QQ) is V
            True
        """
        return AmbientVectorSpace(base_ring, self._names)

    def symmetric_space(self):
        r"""
        Construct the symmetric space of ``self``.

        Consider a hyperplane arrangement `A` in the vector space
        `V = k^n`, for some field `k`. The symmetric space is the
        symmetric algebra `S(V^*)` as the polynomial ring
        `k[x_1, x_2, \ldots, x_n]` where `(x_1, x_2, \ldots, x_n)` is
        a basis for `V`.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H.ambient_space()
            sage: A.symmetric_space()
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(self.base_ring(), self.variable_names())
