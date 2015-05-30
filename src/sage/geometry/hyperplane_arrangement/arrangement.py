r"""
Hyperplane Arrangements

Before talking about hyperplane arrangements, let us start with
individual hyperplanes. This package uses certain linear expressions
to represent hyperplanes, that is, a linear expression `3x + 3y - 5z - 7`
stands for the hyperplane with the equation `x + 3y - 5z = 7`. To create it
in Sage, you first have to create a :class:`HyperplaneArrangements`
object to define the variables `x`, `y`, `z`::

    sage: H.<x,y,z> = HyperplaneArrangements(QQ)
    sage: h = 3*x + 2*y - 5*z - 7;  h
    Hyperplane 3*x + 2*y - 5*z - 7
    sage: h.normal()
    (3, 2, -5)
    sage: h.constant_term()
    -7

The individual hyperplanes behave like the linear expression with
regard to addition and scalar multiplication, which is why you can do
linear combinations of the coordinates::

    sage: -2*h
    Hyperplane -6*x - 4*y + 10*z + 14
    sage: x, y, z
    (Hyperplane x + 0*y + 0*z + 0, Hyperplane 0*x + y + 0*z + 0, Hyperplane 0*x + 0*y + z + 0)

See :mod:`sage.geometry.hyperplane_arrangement.hyperplane` for more
functionality of the individual hyperplanes.

Arrangements
------------

There are several ways to create hyperplane arrangements:

Notation (i): by passing individual hyperplanes to the
:class:`HyperplaneArrangements` object::

    sage: H.<x,y> = HyperplaneArrangements(QQ)
    sage: box = x | y | x-1 | y-1;  box
    Arrangement <y - 1 | y | x - 1 | x>
    sage: box == H(x, y, x-1, y-1)    # alternative syntax 
    True
 
Notation (ii): by passing anything that defines a hyperplane, for
example a coefficient vector and constant term::

    sage: H = HyperplaneArrangements(QQ, ('x', 'y'))
    sage: triangle = H([(1, 0), 0], [(0, 1), 0], [(1,1), -1]);  triangle
    Arrangement <y | x | x + y - 1>

    sage: H.inject_variables()
    Defining x, y
    sage: triangle == x | y | x+y-1
    True

The default base field is `\QQ`, the rational numbers.  Finite fields are also
supported::

    sage: H.<x,y,z> = HyperplaneArrangements(GF(5))
    sage: a = H([(1,2,3), 4], [(5,6,7), 8]);  a
    Arrangement <y + 2*z + 3 | x + 2*y + 3*z + 4>

Notation (iii): a list or tuple of hyperplanes::

    sage: H.<x,y,z> = HyperplaneArrangements(GF(5))
    sage: k = [x+i for i in range(4)];  k
    [Hyperplane x + 0*y + 0*z + 0, Hyperplane x + 0*y + 0*z + 1, 
     Hyperplane x + 0*y + 0*z + 2, Hyperplane x + 0*y + 0*z + 3]
    sage: H(k)
    Arrangement <x | x + 1 | x + 2 | x + 3>

Notation (iv): using the library of arrangements::

    sage: hyperplane_arrangements.braid(4)
    Arrangement of 6 hyperplanes of dimension 4 and rank 3
    sage: hyperplane_arrangements.semiorder(3)
    Arrangement of 6 hyperplanes of dimension 3 and rank 2
    sage: hyperplane_arrangements.graphical(graphs.PetersenGraph())
    Arrangement of 15 hyperplanes of dimension 10 and rank 9
    sage: hyperplane_arrangements.Ish(5)
    Arrangement of 20 hyperplanes of dimension 5 and rank 4

Notation (v): from the bounding hyperplanes of a polyhedron::

    sage: a = polytopes.cube().hyperplane_arrangement();  a
    Arrangement of 6 hyperplanes of dimension 3 and rank 3
    sage: a.n_regions()
    27

New arrangements from old::

    sage: a = hyperplane_arrangements.braid(3)
    sage: b = a.add_hyperplane([4, 1, 2, 3])     
    sage: b
    Arrangement <t1 - t2 | t0 - t1 | t0 - t2 | t0 + 2*t1 + 3*t2 + 4>
    sage: c = b.deletion([4, 1, 2, 3])
    sage: a == c
    True

    sage: a = hyperplane_arrangements.braid(3)
    sage: b = a.union(hyperplane_arrangements.semiorder(3))
    sage: b == a | hyperplane_arrangements.semiorder(3)    # alternate syntax
    True
    sage: b == hyperplane_arrangements.Catalan(3)
    True

    sage: a
    Arrangement <t1 - t2 | t0 - t1 | t0 - t2>
    sage: a = hyperplane_arrangements.coordinate(4)
    sage: h = a.hyperplanes()[0]
    sage: b = a.restriction(h)
    sage: b == hyperplane_arrangements.coordinate(3)
    True

A hyperplane arrangement is *essential* is the normals to its
hyperplane span the ambient space.  Otherwise, it is *inessential*.
The essentialization is formed by intersecting the hyperplanes by this
normal space (actually, it is a bit more complicated over finite
fields)::

    sage: a = hyperplane_arrangements.braid(4);  a
    Arrangement of 6 hyperplanes of dimension 4 and rank 3
    sage: a.is_essential()
    False
    sage: a.rank() < a.dimension()  # double-check
    True
    sage: a.essentialization()
    Arrangement of 6 hyperplanes of dimension 3 and rank 3

The connected components of the complement of the hyperplanes of an arrangement
in `\RR^n` are called the *regions* of the arrangement::

    sage: a = hyperplane_arrangements.semiorder(3)
    sage: b = a.essentialization();   b
    Arrangement of 6 hyperplanes of dimension 2 and rank 2
    sage: b.n_regions()
    19
    sage: b.regions()
    (A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays)
    sage: b.bounded_regions()
    (A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, 
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices)
    sage: b.n_bounded_regions()
    7
    sage: a.unbounded_regions()
    (A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line)

The distance between regions is defined as the number of hyperplanes
separating them. For example::

    sage: r1 = b.regions()[0]
    sage: r2 = b.regions()[1]  
    sage: b.distance_between_regions(r1, r2)
    1
    sage: [hyp for hyp in b if b.is_separating_hyperplane(r1, r2, hyp)]
    [Hyperplane 2*t1 + t2 + 1]
    sage: b.distance_enumerator(r1)  # generating function for distances from r1
    6*x^3 + 6*x^2 + 6*x + 1

.. NOTE::

    *bounded region* really mean *relatively bounded* here.  A region is
    relatively bounded if its intersection with space spanned by the normals
    of the hyperplanes in the arrangement is bounded.

The intersection poset of a hyperplane arrangement is the collection
of all nonempty intersections of hyperplanes in the arrangement,
ordered by reverse inclusion.  It includes the ambient space of the
arrangement (as the intersection over the empty set)::

    sage: a = hyperplane_arrangements.braid(3)
    sage: p = a.intersection_poset()
    sage: p.is_ranked()
    True
    sage: p.order_polytope()
    A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 10 vertices

The characteristic polynomial is a basic invariant of a hyperplane
arrangement. It is defined as

.. MATH::

    \chi(x) := \sum_{w\in P} \mu(w) x^{dim(w)}

where the sum is `P` is the
:meth:`~HyperplaneArrangementElement.intersection_poset` of the
arrangement and `\mu` is the Moebius function of `P`::

    sage: a = hyperplane_arrangements.semiorder(5)
    sage: a.characteristic_polynomial()               # long time (about a second on Core i7)
    x^5 - 20*x^4 + 180*x^3 - 790*x^2 + 1380*x
    sage: a.poincare_polynomial()                     # long time
    1380*x^4 + 790*x^3 + 180*x^2 + 20*x + 1
    sage: a.n_regions()                               # long time
    2371
    sage: charpoly = a.characteristic_polynomial()    # long time
    sage: charpoly(-1)                                # long time
    -2371
    sage: a.n_bounded_regions()                       # long time
    751
    sage: charpoly(1)                                 # long time
    751

For finer invariants derived from the intersection poset, see
:meth:`~HyperplaneArrangementElement.whitney_number` and
:meth:`~HyperplaneArrangementElement.doubly_indexed_whitney_number`.

Miscellaneous methods (see documentation for an explanation)::

    sage: a = hyperplane_arrangements.semiorder(3)
    sage: a.has_good_reduction(5)
    True
    sage: b = a.change_ring(GF(5))
    sage: pa = a.intersection_poset()
    sage: pb = b.intersection_poset()
    sage: pa.is_isomorphic(pb)
    True
    sage: a.face_vector()
    (0, 12, 30, 19)
    sage: a.face_vector()
    (0, 12, 30, 19)
    sage: a.is_central()
    False
    sage: a.is_linear()
    False
    sage: a.sign_vector((1,1,1))
    (-1, 1, -1, 1, -1, 1)
    sage: a.varchenko_matrix()
    [          1          h2       h2*h4       h2*h3    h2*h3*h4 h2*h3*h4*h5]
    [         h2           1          h4          h3       h3*h4    h3*h4*h5]
    [      h2*h4          h4           1       h3*h4          h3       h3*h5]
    [      h2*h3          h3       h3*h4           1          h4       h4*h5]
    [   h2*h3*h4       h3*h4          h3          h4           1          h5]
    [h2*h3*h4*h5    h3*h4*h5       h3*h5       h4*h5          h5           1]

There are extensive methods for visualizing hyperplane arrangements in
low dimensions.  See :meth:`~HyperplaneArrangementElement.plot` for
details.

TESTS::

    sage: H.<x,y> = HyperplaneArrangements(QQ)
    sage: h = H([(1, 106), 106266], [(83, 101), 157866], [(111, 110), 186150], [(453, 221), 532686], 
    ....:       [(407, 237), 516882], [(55, 32), 75620], [(221, 114), 289346], [(452, 115), 474217], 
    ....:       [(406, 131), 453521], [(28, 9), 32446], [(287, 19), 271774], [(241, 35), 244022],
    ....:       [(231, 1), 210984], [(185, 17), 181508], [(23, -8), 16609])
    sage: h.n_regions()
    85

    sage: H()
    Empty hyperplane arrangement of dimension 2

    sage: Zero = HyperplaneArrangements(QQ)
    sage: Zero
    Hyperplane arrangements in 0-dimensional linear space over Rational Field with coordinate 
    sage: Zero()
    Empty hyperplane arrangement of dimension 0
    sage: Zero.an_element()
    Empty hyperplane arrangement of dimension 0

AUTHORS:

- David Perkinson (2013-06): initial version

- Qiaoyu Yang (2013-07)

- Kuai Yu (2013-07)

- Volker Braun (2013-10): Better Sage integration, major code refactoring.

This module implements hyperplane arrangements defined over the
rationals or over finite fields.  The original motivation was to make
a companion to Richard Stanley's notes [RS]_ on hyperplane
arrangements.

REFERENCES:

..  [RS] 
    Stanley, Richard: *Hyperplane Arrangements*,
    Geometric Combinatorics (E. Miller, V. Reiner, and B. Sturmfels, eds.),
    IAS/Park City Mathematics Series, vol. 13, American Mathematical Society,
    Providence, RI, 2007, pp. 389-496.
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

# Possible extensions for hyperplane_arrangement.py:
# - the big face lattice
# - Orlik-Solomon algebras
# - create ties with the Sage matroid methods
# - hyperplane arrangements over other fields

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import QQ, ZZ
from sage.misc.cachefunc import cached_method
from sage.misc.misc import uniq
from sage.matrix.constructor import matrix, vector

from sage.geometry.hyperplane_arrangement.hyperplane import AmbientVectorSpace, Hyperplane



class HyperplaneArrangementElement(Element):
    """
    An element in a hyperplane arrangement.

    .. WARNING::

        You should never create
        :class:`HyperplaneArrangementElement` instances directly,
        always use the parent.
    """
    def __init__(self, parent, hyperplanes, check=True):
        """
        Construct a hyperplane arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`HyperplaneArrangements`
       
        - ``hyperplanes`` -- a tuple of hyperplanes

        - ``check`` -- boolean (optional; default ``True``); whether
          to check input

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement <y | x>
            sage: TestSuite(elt).run()
        """
        super(HyperplaneArrangementElement, self).__init__(parent)
        self._hyperplanes = hyperplanes
        if check:
            if not isinstance(hyperplanes, tuple):
                raise ValueError("the hyperplanes must be given as a tuple")
            if not all(isinstance(h, Hyperplane) for h in hyperplanes):
                raise ValueError("not all elements are hyperplanes")
            if not all(h.parent() is self.parent().ambient_space() for h in hyperplanes):
                raise ValueError("not all hyperplanes are in the ambient space")

    def _first_ngens(self, n):
        """
        Workaround to support the construction with names.

        INPUT/OUTPUT:

        See :meth:`HyperplaneArrangements._first_ngens`.

        EXAMPLES::

            sage: a.<x,y,z> = hyperplane_arrangements.braid(3)   # indirect doctest
            sage: (x, y) == a._first_ngens(2)
            True
        """
        return self.parent()._first_ngens(n)

    def __getitem__(self, i):
        """
        Return the `i`-th hyperplane.

        INPUT:
        
        - ``i`` -- integer

        OUTPUT:

        The `i`-th hyperplane.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = x|y;  h
            Arrangement <y | x>
            sage: h[0]
            Hyperplane 0*x + y + 0
            sage: h[1] 
            Hyperplane x + 0*y + 0
        """
        return self._hyperplanes[i]
        
    def n_hyperplanes(self):
        r"""
        Return the number of hyperplanes in the arrangement.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,1,0], [2,3,-1], [4,5,3])
            sage: A.n_hyperplanes()
            3
            sage: len(A)    # equivalent
            3
        """
        return len(self._hyperplanes)

    __len__ = n_hyperplanes

    def hyperplanes(self):
        r"""
        Return the number of hyperplanes in the arrangement.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,1,0], [2,3,-1], [4,5,3])
            sage: A.hyperplanes()
            (Hyperplane x + 0*y + 1, Hyperplane 3*x - y + 2, Hyperplane 5*x + 3*y + 4)

        Note that the hyperplanes can be indexed as if they were a list::

            sage: A[0]
            Hyperplane x + 0*y + 1
        """
        return self._hyperplanes

    def _repr_(self):
        r"""
        String representation for a hyperplane arrangement.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: H(x, y, x-1, y-1)
            Arrangement <y - 1 | y | x - 1 | x>
            sage: x | y | x - 1 | y - 1 | x + y | x - y
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
            sage: H()
            Empty hyperplane arrangement of dimension 2
        """
        if len(self) == 0:
            return 'Empty hyperplane arrangement of dimension {0}'.format(self.dimension())
        elif len(self) < 5:
            hyperplanes = ' | '.join(h._repr_linear(include_zero=False) for h in self._hyperplanes)
            return 'Arrangement <{0}>'.format(hyperplanes)
        return 'Arrangement of {0} hyperplanes of dimension {1} and rank {2}'.format(
            len(self), self.dimension(), self.rank())

    def dimension(self):
        """
        Return the ambient space dimension of the arrangement.

        OUTPUT:

        An integer.
        
        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: (x | x-1 | x+1).dimension()
            2
            sage: H(x).dimension()
            2
        """
        return self.parent().ngens()

    def rank(self):
        """
        Return the rank.

        OUTPUT:

        The dimension of the span of the normals to the
        hyperplanes in the arrangement.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H([[0, 1, 2, 3],[-3, 4, 5, 6]])
            sage: A.dimension()
            3
            sage: A.rank()
            2

            sage: B = hyperplane_arrangements.braid(3)
            sage: B.hyperplanes()
            (Hyperplane 0*t0 + t1 - t2 + 0, 
             Hyperplane t0 - t1 + 0*t2 + 0, 
             Hyperplane t0 + 0*t1 - t2 + 0)
            sage: B.dimension()
            3
            sage: B.rank()
            2

            sage: p = polytopes.simplex(5, project=True)
            sage: H = p.hyperplane_arrangement()
            sage: H.rank()
            5
        """
        R = self.parent().base_ring()
        normals = [h.normal() for h in self]
        return matrix(R, normals).rank()

    def __cmp__(self, other):
        """
        Compare two hyperplane arrangements.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: H(x) == H(y)
            False
            sage: H(x) == 0
            False
        """
        assert (isinstance(self, type(other))) and (self.parent() is other.parent()) # guaranteed by framework
        return cmp(self._hyperplanes, other._hyperplanes)

    def union(self, other):
        r"""
        The union of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a hyperplane arrangement or something that can
          be converted into a hyperplane arrangement

        OUTPUT:

        A new hyperplane arrangement.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,2,3], [0,1,1], [0,1,-1], [1,-1,0], [1,1,0])
            sage: B = H([1,1,1], [1,-1,1], [1,0,-1])
            sage: A.union(B)
            Arrangement of 8 hyperplanes of dimension 2 and rank 2
            sage: A | B   # syntactic sugar
            Arrangement of 8 hyperplanes of dimension 2 and rank 2

        A single hyperplane is coerced into a hyperplane arrangement
        if necessary::

            sage: A.union(x+y-1)
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
            sage: A.add_hyperplane(x+y-1)    # alias
            Arrangement of 6 hyperplanes of dimension 2 and rank 2

            sage: P.<x,y> = HyperplaneArrangements(RR)
            sage: C = P(2*x + 4*y + 5)
            sage: C.union(A)
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
        """
        P = self.parent()
        other = P(other)
        hyperplanes = self._hyperplanes + other._hyperplanes
        return P(*hyperplanes)

    add_hyperplane = union

    __or__ = union

    def plot(self, **kwds):
        """
        Plot the hyperplane arrangement.

        OUTPUT:

        A graphics object.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: L(x, y, x+y-2).plot()
            Graphics object consisting of 3 graphics primitives
        """
        from sage.geometry.hyperplane_arrangement.plot import plot
        return plot(self, **kwds)

    def cone(self, variable='t'):
        r"""
        Return the cone over the hyperplane arrangement.  

        INPUT:

        - ``variable`` -- string; the name of the additional variable

        OUTPUT:

        A new yperplane arrangement. Its equations consist of
        `[0, -d, a_1, \ldots, a_n]` for each `[d, a_1, \ldots, a_n]` in the
        original arrangement and the equation `[0, 1, 0, \ldots, 0]`.

        EXAMPLES::

            sage: a.<x,y,z> = hyperplane_arrangements.semiorder(3)
            sage: b = a.cone()
            sage: a.characteristic_polynomial().factor()
            x * (x^2 - 6*x + 12)
            sage: b.characteristic_polynomial().factor()
            (x - 1) * x * (x^2 - 6*x + 12)
            sage: a.hyperplanes()
            (Hyperplane 0*x + y - z - 1, 
             Hyperplane 0*x + y - z + 1, 
             Hyperplane x - y + 0*z - 1, 
             Hyperplane x - y + 0*z + 1, 
             Hyperplane x + 0*y - z - 1, 
             Hyperplane x + 0*y - z + 1)
            sage: b.hyperplanes()
            (Hyperplane -t + 0*x + y - z + 0, 
             Hyperplane -t + x - y + 0*z + 0, 
             Hyperplane -t + x + 0*y - z + 0, 
             Hyperplane t + 0*x + 0*y + 0*z + 0, 
             Hyperplane t + 0*x + y - z + 0, 
             Hyperplane t + x - y + 0*z + 0, 
             Hyperplane t + x + 0*y - z + 0)
        """
        hyperplanes = []
        for h in self.hyperplanes():
            new_h = [0] + [h.b()] + list(h.A())
            hyperplanes.append(new_h)
        hyperplanes.append([0, 1] + [0] * self.dimension())
        P = self.parent()
        names = (variable,) + P._names
        H = HyperplaneArrangements(self.parent().base_ring(), names=names)
        return H(*hyperplanes)
        
    @cached_method
    def intersection_poset(self):
        r"""
        Return the intersection poset of the hyperplane arrangement.

        OUTPUT:

        The poset of non-empty intersections of hyperplanes.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a.intersection_poset()
            Finite poset containing 4 elements

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.intersection_poset()
            Finite poset containing 19 elements
        """
        K = self.base_ring()
        from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
        from sage.modules.all import VectorSpace
        whole_space = AffineSubspace(0, VectorSpace(K, self.dimension()))
        L = [[whole_space]]
        active = True
        codim = 0
        while active:
            active = False
            new_level = []
            for T in L[codim]:
                for H in self:
                    I = H._affine_subspace().intersection(T)
                    if I is not None and I != T and I not in new_level:
                        new_level.append(I)
                        active = True
            if active:
                L.append(new_level)
            codim += 1
        from sage.misc.flatten import flatten
        L = flatten(L)
        t = {}
        for i in range(len(L)):
            t[i] = L[i]
        cmp_fn = lambda p, q: t[q] < t[p]
        from sage.combinat.posets.posets import Poset
        return Poset((t, cmp_fn))

    def _slow_characteristic_polynomial(self):
        """
        Return the characteristic polynomial of the hyperplane arrangement.

        This is the slow computation directly from the definition. For
        educational use only.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a._slow_characteristic_polynomial()
            x^2 - 2*x + 1
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        x = polygen(QQ, 'x')
        P = self.intersection_poset()
        n = self.dimension()
        return sum([P.mobius_function(0, p) * x**(n - P.rank(p)) for p in P])
        
    @cached_method
    def characteristic_polynomial(self):
        r"""
        Return the characteristic polynomial of the hyperplane arrangement.

        OUTPUT:

        The characteristic polynomial in `\QQ[x]`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a.characteristic_polynomial()
            x^2 - 2*x + 1
        
        TESTS::
          
            sage: H.<s,t,u,v> = HyperplaneArrangements(QQ)
            sage: m = matrix([(0, -1, 0, 1, -1), (0, -1, 1, -1, 0), (0, -1, 1, 0, -1),
            ....:   (0, 1, 0, 0, 0), (0, 1, 0, 1, -1), (0, 1, 1, -1, 0), (0, 1, 1, 0, -1)])
            sage: R.<x> = QQ[]
            sage: expected_charpoly = (x - 1) * x * (x^2 - 6*x + 12)
            sage: for s in SymmetricGroup(4):   # long time (about a second on a Core i7)
            ....:     m_perm = [m.column(i) for i in [0, s(1), s(2), s(3), s(4)]]
            ....:     m_perm = matrix(m_perm).transpose()
            ....:     charpoly = H(m_perm.rows()).characteristic_polynomial()
            ....:     assert charpoly == expected_charpoly
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        x = polygen(QQ, 'x')
        if self.rank() == 1:
            return x**(self.dimension() - 1) * (x - len(self))

        H = self[0]
        R = self.restriction(H)
        charpoly_R = R.characteristic_polynomial()
        D = self.deletion(H)
        charpoly_D = D.characteristic_polynomial()
        return charpoly_D - charpoly_R

    @cached_method
    def poincare_polynomial(self):
        r"""
        Return the Poincare polynomial of the hyperplane arrangement.

        OUTPUT:

        The Poincare polynomial in `\QQ[x]`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a.poincare_polynomial()
            x^2 + 2*x + 1
        """
        charpoly = self.characteristic_polynomial()
        R = charpoly.parent()
        x = R.gen(0)
        poincare = (-x)**self.dimension() * charpoly(-QQ(1)/x)
        return R(poincare)

    def deletion(self, hyperplanes):
        r"""
        Return the hyperplane arrangement obtained by removing ``h``.

        INPUT:

        - ``h`` -- a hyperplane or hyperplane arrangement

        OUTPUT:

        A new hyperplane arrangement with the given hyperplane(s)
        ``h`` removed.

        .. SEEALSO::

            :meth:`restriction`

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([0,1,0], [1,0,1], [-1,0,1], [0,1,-1], [0,1,1]);  A
            Arrangement of 5 hyperplanes of dimension 2 and rank 2
            sage: A.deletion(x)
            Arrangement <y - 1 | y + 1 | x - y | x + y>
            sage: h = H([0,1,0], [0,1,1])
            sage: A.deletion(h)
            Arrangement <y - 1 | y + 1 | x - y>

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([0,1,0], [1,0,1], [-1,0,1], [0,1,-1], [0,1,1])
            sage: h = H([0,4,0])
            sage: A.deletion(h)
            Arrangement <y - 1 | y + 1 | x - y | x + y>
            sage: l = H([1,2,3])
            sage: A.deletion(l)
            Traceback (most recent call last):
            ...
            ValueError: hyperplane is not in the arrangement
        """
        parent = self.parent()
        hyperplanes = parent(hyperplanes)
        planes = list(self)
        for hyperplane in hyperplanes:
            try:
                planes.remove(hyperplane)
            except ValueError:
                raise ValueError('hyperplane is not in the arrangement')
        return parent(*planes)

    def restriction(self, hyperplane):
        r"""
        Return the restriction to a hyperplane.

        INPUT:

        - ``hyperplane`` -- a hyperplane of the hyperplane arrangement

        OUTPUT:

        The restriction of the hyperplane arrangement to the given
        ``hyperplane``.

        EXAMPLES::

            sage: A.<u,x,y,z> = hyperplane_arrangements.braid(4);  A
            Arrangement of 6 hyperplanes of dimension 4 and rank 3
            sage: H = A[0];  H
            Hyperplane 0*u + 0*x + y - z + 0
            sage: R = A.restriction(H);  R
            Arrangement <x - z | u - x | u - z>
            sage: D = A.deletion(H);  D
            Arrangement of 5 hyperplanes of dimension 4 and rank 3
            sage: ca = A.characteristic_polynomial()
            sage: cr = R.characteristic_polynomial()
            sage: cd = D.characteristic_polynomial()
            sage: ca
            x^4 - 6*x^3 + 11*x^2 - 6*x
            sage: cd - cr
            x^4 - 6*x^3 + 11*x^2 - 6*x

        .. SEEALSO::

            :meth:`deletion`
        """
        parent = self.parent()
        hyperplane = parent(hyperplane)[0]
        if hyperplane not in self.hyperplanes():
            raise ValueError('hyperplane not in arrangement')
        pivot = hyperplane._normal_pivot()
        hyperplanes = []
        for h in self:
            rescale = h.A()[pivot] / hyperplane.A()[pivot]
            h = h - rescale * hyperplane
            A = list(h.A())
            A_pivot = A.pop(pivot)
            assert A_pivot == 0
            if all(a == 0 for a in A):
                continue
            b = h.b()
            hyperplanes.append([A, b])
        names = list(parent._names)
        names.pop(pivot)
        H = HyperplaneArrangements(parent.base_ring(), names=tuple(names))
        return H(*hyperplanes, signed=False)

    def change_ring(self, base_ring):
        """
        Return hyperplane arrangement over the new base ring.
        
        INPUT:

        - ``base_ring`` -- the new base ring; must be a field for
          hyperplane arrangements

        OUTPUT:

        The hyperplane arrangement obtained by changing the base
        field, as a new hyperplane arrangement.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,1), 0], [(2,3), -1])
            sage: A.change_ring(FiniteField(2))
            Arrangement <y + 1 | x + y>
        """
        parent = self.parent().change_ring(base_ring)
        return parent(self)

    @cached_method
    def n_regions(self):
        r"""
        The number of regions of the hyperplane arrangement.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.n_regions()
            19

        TESTS::
         
            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,1), 0], [(2,3), -1], [(4,5), 3])
            sage: B = A.change_ring(FiniteField(7))
            sage: B.n_regions()
            Traceback (most recent call last):
            ...
            TypeError: base field must have characteristic zero
        """
        if self.base_ring().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        charpoly = self.characteristic_polynomial()
        return (-1)**self.dimension() * charpoly(-1)

    @cached_method
    def n_bounded_regions(self):
        r"""
        Return the number of (relatively) bounded regions.

        OUTPUT:

        An integer. The number of relatively bounded regions of the
        hyperplane arrangement.

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.n_bounded_regions()
            7

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,1),0], [(2,3),-1], [(4,5),3])
            sage: B = A.change_ring(FiniteField(7))
            sage: B.n_bounded_regions()
            Traceback (most recent call last):
            ...
            TypeError: base field must have characteristic zero
        """
        if self.base_ring().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        charpoly = self.characteristic_polynomial()
        return (-1)**self.rank() * charpoly(1)

    def has_good_reduction(self, p):
        r"""
        Return whether the hyperplane arrangement has good reduction mod `p`.

        Let `A` be a hyperplane arrangement with equations defined
        over the integers, and let `B` be the hyperplane arrangement
        defined by reducing these equations modulo a prime `p`.  Then
        `A` has good reduction modulo `p` if the intersection posets
        of `A` and `B` are isomorphic.

        INPUT:

        - ``p`` -- prime number

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a.has_good_reduction(5)
            True
            sage: a.has_good_reduction(3)
            False
            sage: b = a.change_ring(GF(3))
            sage: a.characteristic_polynomial()
            x^3 - 6*x^2 + 12*x
            sage: b.characteristic_polynomial()  # not equal to that for a
            x^3 - 6*x^2 + 10*x
        """
        if self.base_ring() != QQ:
            raise TypeError('arrangement must be defined over QQ')
        if not p.is_prime():
            raise TypeError('must reduce modulo a prime number')
        from sage.rings.all import GF
        a = self.change_ring(GF(p))
        p = self.intersection_poset()
        q = a.intersection_poset()
        return p.is_isomorphic(q)

    def is_linear(self):
        r"""
        Test whether all hyperplanes pass through the origin.

        OUTPUT:

        A boolean. Whether all the hyperplanes pass through the origin.

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a.is_linear()
            False
            sage: b = hyperplane_arrangements.braid(3)
            sage: b.is_linear()
            True
        
            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: c = H(x+1, y+1)
            sage: c.is_linear()
            False
            sage: c.is_central()
            True
        """
        return all(hyperplane.b() == 0 for hyperplane in self)

    def is_essential(self):
        r"""
        Test whether the hyperplane arrangement is essential.

        A hyperplane arrangement is essential if the span of the normals
        of its hyperplanes spans the ambient space.

        .. SEEALSO::

            :meth:`essentialization`

        OUTPUT:

        A boolean indicating whether the hyperplane arrangement is essential.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: H(x, x+1).is_essential()
            False
            sage: H(x, y).is_essential()
            True
        """
        return self.rank() == self.dimension()

    @cached_method
    def is_central(self):
        r"""
        Test whether the intersection of all the hyperplanes is nonempty.

        OUTPUT:

        A boolean whether the hyperplane arrangement is such that the
        intersection of all the hyperplanes in the arrangement is
        nonempty.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)
            sage: a.is_central()
            True
        """
        R = self.base_ring()
        m = matrix(R, [h.normal() for h in self])
        b = vector(R, [h.b() for h in self])
        try:
            m.solve_right(b)
            return True
        except ValueError:
            return False

    @cached_method
    def essentialization(self):
        r"""
        Return the essentialization of the hyperplane arrangement.

        The essentialization of a hyperplane arrangement whose base field
        has characteristic 0 is obtained by intersecting the hyperplanes by
        the space spanned by their normal vectors.

        OUTPUT:

        The essentialization as a new hyperplane arrangement.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(3)
            sage: a.is_essential()
            False
            sage: a.essentialization()
            Arrangement <t1 - t2 | t1 + 2*t2 | 2*t1 + t2>

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: B = H([(1,0),1], [(1,0),-1])
            sage: B.is_essential()
            False
            sage: B.essentialization()
            Arrangement <-x + 1 | x + 1>
            sage: B.essentialization().parent()
            Hyperplane arrangements in 1-dimensional linear space over 
            Rational Field with coordinate x
          
            sage: H.<x,y> = HyperplaneArrangements(GF(2))
            sage: C = H([(1,1),1], [(1,1),0])
            sage: C.essentialization()
            Arrangement <y | y + 1>

            sage: h = hyperplane_arrangements.semiorder(4)
            sage: h.essentialization()
            Arrangement of 12 hyperplanes of dimension 3 and rank 3

        TESTS::

            sage: b = hyperplane_arrangements.coordinate(2)
            sage: b.is_essential()
            True
            sage: b.essentialization() is b
            True
        """
        def echelon_col_iter(row_iter):
            """helper to iterat over the echelon pivot column indices"""
            for row in row_iter:
                if row == 0:
                    raise StopIteration
                for pivot in range(self.dimension()):
                    if row[pivot] != 0:
                        break
                assert row[pivot] == 1
                yield pivot, row

        if self.is_essential():
            return self
        parent = self.parent()
        H = parent.ambient_space()
        R = parent.base_ring()
        hyperplanes = self.hyperplanes()
        normals = matrix(R, [h.normal() for h in self]).transpose()
        # find a (any) complement to the normals
        if R.characteristic() == 0:
            complement_basis = normals.kernel().echelonized_basis()
        else:
            # we don't necessarily have an orthogonal complement, pick any complement
            complement_basis = []
            for pivot, row in echelon_col_iter(normals.echelon_form().rows()):
                v = [0] * self.dimension()
                v[pivot] = 1
                complement_basis.append(vector(R, v))
        # reduce the hyperplane equations
        echelon_pivots = []   # the column indices where N has 1's from the echelonization
        for pivot, row in echelon_col_iter(complement_basis):
            assert row[pivot] == 1
            echelon_pivots.append(pivot)
            hyperplanes = [h - h.A()[pivot] * H(row, 0) for h in hyperplanes]
        # eliminate the pivot'ed coodinates
        restricted = []
        for h in hyperplanes:
            A = h.A()
            if A == 0:
                continue
            A = [A[i] for i in range(self.dimension()) if i not in echelon_pivots]
            b = h.b()
            restricted.append([A, b])
        names = tuple(name for i, name in enumerate(parent._names) if i not in echelon_pivots)
        # Construct the result
        restricted_parent = HyperplaneArrangements(R, names=names)
        return restricted_parent(*restricted, signed=False)

    def sign_vector(self, p):
        r"""
        Indicates on which side of each hyperplane the given
        point `p` lies.

        The base field must have characteristic zero.

        INPUT:

        - ``p`` -- point as a list/tuple/iterable

        OUTPUT:

        A vector whose entries are in `[-1, 0, +1]`.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,0), 0], [(0,1), 1]);  A
            Arrangement <y + 1 | x>
            sage: A.sign_vector([2, -2])
            (-1, 1)
            sage: A.sign_vector((-1, -1))
            (0, -1)

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(GF(3))
            sage: A = H(x, y)
            sage: A.sign_vector([1, 2])
            Traceback (most recent call last):
            ...
            ValueError: characteristic must be zero
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('characteristic must be zero')
        from sage.functions.generalized import sign
        values = [hyperplane(p) for hyperplane in self]
        signs = vector(ZZ, map(sign, values))
        signs.set_immutable()
        return signs
    
    def face_vector(self):
        r"""
        Return the face vector.

        OUTPUT:

        A vector of integers.

        The `d`-th entry is the number of faces of dimension `d`.  A
        *face* is is the intersection of a region with a hyperplane of
        the arrangehment.

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.face_vector()
            (0, 6, 21, 16)
        """
        m = self.whitney_data()[0]
        v = list(sum(m.transpose().apply_map(abs)))
        v.reverse()
        v = vector(ZZ, [0]*(self.dimension() - self.rank()) + v)
        v.set_immutable()
        return v

    @cached_method
    def _parallel_hyperplanes(self):
        """
        Return the hyperplanes grouped into parallel sets.

        OUTPUT:

        A tuple with one entry per set of parallel hyperplanes. Each
        entry is a tuple of triples, one for each parallel hyperplane
        in the parallel set. The triple consists of the hyperplane,
        the normal vector `A`, and the constant `b` of the hyperplane
        equation `Ax+b`. The normalization is such that `A` is the
        same for each hyperplane of the parallel set, and the order is
        in increasing order of the `b` values.

        In other words, each parallel set of hyperplanes is also
        ordered by the order with which a common normal passes through
        them.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = (x + 2*y | 2*x + 4*y + 1 | -x/4 - y/2 + 1);  h
            Arrangement <-x - 2*y + 4 | x + 2*y | 2*x + 4*y + 1>
            sage: h._parallel_hyperplanes()[0]
            ((Hyperplane -x - 2*y + 4, (1, 2), -4),
             (Hyperplane x + 2*y + 0, (1, 2), 0),
             (Hyperplane 2*x + 4*y + 1, (1, 2), 1/2))

           sage: hyperplane_arrangements.Shi(3)._parallel_hyperplanes()
           (((Hyperplane 0*t0 + t1 - t2 - 1, (0, 1, -1), -1),
             (Hyperplane 0*t0 + t1 - t2 + 0, (0, 1, -1), 0)),
            ((Hyperplane t0 - t1 + 0*t2 - 1, (1, -1, 0), -1),
             (Hyperplane t0 - t1 + 0*t2 + 0, (1, -1, 0), 0)),
            ((Hyperplane t0 + 0*t1 - t2 - 1, (1, 0, -1), -1),
             (Hyperplane t0 + 0*t1 - t2 + 0, (1, 0, -1), 0)))
        """
        V = self.parent().ambient_space()
        parallels = dict()
        for hyperplane in self:
            through_origin = V([list(hyperplane.A()), 0]).primitive(signed=False)
            parallel_planes = parallels.get(through_origin, [])
            A = through_origin.A()
            b = hyperplane.b() * (A / hyperplane.A())
            parallel_planes.append([b, (hyperplane, A, b)])
            parallels[through_origin] = parallel_planes
        parallels = [tuple(tuple(hyperplane[1]
                           for hyperplane in sorted(parallels[key])))
                     for key in parallels.keys()]
        return tuple(sorted(parallels))

    def vertices(self, exclude_sandwiched=False):
        """
        Return the vertices.

        The vertices are the zero-dimensional faces, see
        :meth:`face_vector`.
    
        INPUT:
        
        - ``exclude_sandwiched`` -- boolean (default:
          ``False``). Whether to exclude hyperplanes that are
          sandwiched between parallel hyperplanes. Useful if you only
          need the convex hull.

        OUTPUT:

        The vertices in a sorted tuple. Each vertex is returned as a
        vector in the ambient vector space.
        
        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3).essentialization()
            sage: A.dimension()
            2
            sage: A.face_vector()
            (6, 21, 16)
            sage: A.vertices()
            ((-2/3, 1/3), (-1/3, -1/3), (0, -1), (0, 0), (1/3, -2/3), (2/3, -1/3))
            sage: point2d(A.vertices(), size=20) + A.plot()
            Graphics object consisting of 7 graphics primitives
        
            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: chessboard = []
            sage: N = 8
            sage: for x0, y0 in CartesianProduct(range(N+1), range(N+1)):
            ....:     chessboard.extend([x-x0, y-y0])
            sage: chessboard = H(chessboard)
            sage: len(chessboard.vertices())
            81
            sage: chessboard.vertices(exclude_sandwiched=True)
            ((0, 0), (0, 8), (8, 0), (8, 8))
        """
        from sage.matroids.all import Matroid
        from sage.combinat.cartesian_product import CartesianProduct
        R = self.parent().base_ring()
        parallels = self._parallel_hyperplanes()
        A_list = [parallel[0][1] for parallel in parallels]
        b_list_list = [[-hyperplane[2] for hyperplane in parallel]
                       for parallel in parallels]
        if exclude_sandwiched:
            def skip(b_list):
                if len(b_list) == 1:
                    return b_list
                return [b_list[0], b_list[-1]]
            b_list_list = map(skip, b_list_list)
        M = Matroid(groundset=range(len(parallels)), matrix=matrix(A_list).transpose())
        d = self.dimension()
        # vertices are solutions v * lhs = rhs
        lhs = matrix(R, d, d)
        rhs = vector(R, d)
        vertices = set()
        for indices in M.independent_r_sets(d):
            for row, i in enumerate(indices):
                lhs[row] = A_list[i]
            b_list = [b_list_list[i] for i in indices]
            for b in CartesianProduct(*b_list):
                for i in range(d):
                    rhs[i] = b[i]
                vertex = lhs.solve_right(rhs)
                vertex.set_immutable()
                vertices.add(vertex)
        return tuple(sorted(vertices))

    def _make_region(self, hyperplanes):
        """
        Helper method to construct a region.

        INPUT:

        - ``hyperplanes`` -- a list/tuple/iterable of hyperplanes

        OUTPUT:

        The polyhedron constructed from taking the linear expressions
        as inequalities.
        
        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = H(x)
            sage: h._make_region([x, 1-x, y, 1-y])
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        """
        ieqs = [h.coefficients() for h in hyperplanes]
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(ieqs=ieqs, ambient_dim=self.dimension(), 
                          base_ring=self.parent().base_ring())

    @cached_method
    def regions(self):
        r"""
        Return the regions of the hyperplane arrangement.

        The base field must have characteristic zero.

        OUTPUT:

        A tuple containing the regions as polyhedra.

        The regions are the connected components of the complement of
        the union of the hyperplanes as a subset of `\RR^n`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)
            sage: a.regions()
            (A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex, 1 ray, 1 line, 
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex, 1 ray, 1 line)

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H(x, y+1)
            sage: A.regions()
            (A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays)
        
            sage: chessboard = []
            sage: N = 8
            sage: for x0, y0 in CartesianProduct(range(N+1), range(N+1)):
            ....:     chessboard.extend([x-x0, y-y0])
            sage: chessboard = H(chessboard)
            sage: len(chessboard.bounded_regions())   # long time, 359 ms on a Core i7
            64
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('base field must have characteristic zero')
        from sage.geometry.polyhedron.constructor import Polyhedron
        R = self.base_ring()
        dim = self.dimension()
        universe = Polyhedron(eqns=[[0] + [0] * dim], base_ring=R)
        regions = [universe]
        for hyperplane in self:
            ieq = vector(R, hyperplane.coefficients())
            pos_half = Polyhedron(ieqs=[ ieq], base_ring=R)
            neg_half = Polyhedron(ieqs=[-ieq], base_ring=R)
            subdivided = []
            for region in regions:
                for half_space in pos_half, neg_half:
                    part = region.intersection(half_space)
                    if part.dim() == dim:
                        subdivided.append(part)
            regions = subdivided
        return tuple(regions)

    def region_containing_point(self, p):
        r"""
        The region in the hyperplane arrangement containing a given point.

        The base field must have characteristic zero.

        INPUT:

        - ``p`` -- point

        OUTPUT:

        A polyhedron. A ``ValueError`` is raised if the point is not
        interior to a region, that is, sits on a hyperplane.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,0), 0], [(0,1), 1], [(0,1), -1], [(1,-1), 0], [(1,1), 0])
            sage: A.region_containing_point([1,2])
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 2 rays

        TESTS::

            sage: A = H([(1,1),0], [(2,3),-1], [(4,5),3])
            sage: B = A.change_ring(FiniteField(7))
            sage: B.region_containing_point((1,2))
            Traceback (most recent call last):
            ...
            ValueError: base field must have characteristic zero

            sage: A = H([(1,1),0], [(2,3),-1], [(4,5),3])
            sage: A.region_containing_point((1,-1))
            Traceback (most recent call last):
            ...
            ValueError: point sits on a hyperplane
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('base field must have characteristic zero')
        sign_vector = self.sign_vector(p)
        ieqs = []
        for i, hyperplane in enumerate(self):
            sign = sign_vector[i]
            if sign == 1:
                ieqs.append(hyperplane)
            elif sign == -1:
                ieqs.append(-hyperplane)
            else:
                assert sign == 0
                raise ValueError('point sits on a hyperplane')
        return self._make_region(ieqs)

    @cached_method
    def _bounded_region_indices(self):
        r"""
        Return the relatively bounded regions.

        OUTPUT:

        Tuple of integers. The positions of the relatively bounded
        regions in :meth:`regions`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a._bounded_region_indices()
            (2, 7, 8, 9, 10, 11, 16)
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        normal = Polyhedron(vertices=[[0]*self.dimension()], 
                            lines=[hyperplane.normal() for hyperplane in self])
        if normal.dim() == 0:
            transverse = lambda poly: poly
        else:
            transverse = lambda poly: poly.intersection(normal)
        return tuple(i for i, region in enumerate(self.regions())
                     if transverse(region).is_compact())

    def bounded_regions(self):
        r"""
        Return the relatively bounded regions of the arrangement.

        A region is relatively bounded if its intersection with the space
        spanned by the normals to the hyperplanes is bounded. This is the
        same as being bounded in the case that the hyperplane arrangement
        is essential. It is assumed that the arrangement is defined over
        the rationals.

        OUTPUT:

        Tuple of polyhedra. The relatively bounded regions of the
        arrangement.

        .. SEEALSO::

            :meth:`unbounded_regions`

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.bounded_regions()
            (A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line)
            sage: A.bounded_regions()[0].is_compact()    # the regions are only *relatively* bounded
            False
            sage: A.is_essential()
            False
        """
        return tuple(self.regions()[i] for i in self._bounded_region_indices())

    def unbounded_regions(self):
        r"""
        Return the relatively bounded regions of the arrangement.

        OUTPUT:

        Tuple of polyhedra. The regions of the arrangement that are not
        relatively bounded.  It is assumed that the arrangement is
        defined over the rationals.

        .. SEEALSO::

            :meth:`bounded_regions`

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: B = A.essentialization()
            sage: B.n_regions() - B.n_bounded_regions()
            12
            sage: B.unbounded_regions()
            (A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays, 
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays)
        """
        s = set(range(self.n_regions())).difference(set(self._bounded_region_indices()))
        return tuple(self.regions()[i] for i in s)

    @cached_method
    def whitney_data(self):
        r"""
        Return the Whitney numbers.

        .. SEEALSO::

            :meth:`whitney_number`,
            :meth:`doubly_indexed_whitney_number`

        OUTPUT:

        A pair of integer matrices. The two matrices are the
        doubly-indexed Whitney numbers of the first or second kind,
        respectively. The `i,j`-th entry is the `i,j`-th
        doubly-indexed Whitney number.

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.whitney_data()
            (
            [  1  -6   9]  [ 1  6  6]
            [  0   6 -15]  [ 0  6 15]
            [  0   0   6], [ 0  0  6]
            )
        """
        p = self.intersection_poset()
        r = p.rank_function()
        top = r(p.maximal_elements()[0])
        from sage.matrix.constructor import zero_matrix
        m1 = zero_matrix(ZZ, top+1, top+1)
        m2 = zero_matrix(ZZ, top+1, top+1)
        for i, j in p.relations_iterator():
            m1[r(i), r(j)] += p.mobius_function(i, j)
            m2[r(i), r(j)] += 1
        m1.set_immutable()
        m2.set_immutable()
        return (m1, m2)

    def doubly_indexed_whitney_number(self, i, j, kind=1):
        r"""
        Return the `i,j`-th  doubly-indexed Whitney number.

        If ``kind=1``, this number is obtained by adding the Moebius function
        values `mu(x,y)` over all `x, y` in the intersection poset with
        `\mathrm{rank}(x) = i` and `\mathrm{rank}(y) = j`.

        If `kind=2`, this number is the number of elements `x,y` in the
        intersection poset such that `x \leq y` with ranks `i` and `j`,
        respectively.

        INPUT:

        - ``i``, ``j`` -- integers

        - ``kind`` -- (default: 1) 1 or 2

        OUTPUT:

        Integer. The `(i,j)`-th entry of the ``kind`` Whitney number.

        .. SEEALSO::

            :meth:`whitney_number`,
            :meth:`whitney_data`

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.doubly_indexed_whitney_number(0, 2)
            9
            sage: A.whitney_number(2)
            9
            sage: A.doubly_indexed_whitney_number(1, 2)
            -15

        REFERENCES:

        ..  [GZ] Greene; Zaslavsky
            "On the Interpretation of Whitney Numbers Through Arrangements of
            Hyperplanes, Zonotopes, Non-Radon Partitions, and Orientations of
            Graphs"
            Transactions of the American Mathematical Society, Vol. 280, No. 1.
            (Nov., 1983), pp. 97-126.
        """
        if 0 <= i and j <= self.dimension():
            if kind == 1:
                return self.whitney_data()[0][i, j]
            elif kind == 2:
                return self.whitney_data()[1][i, j]
        raise ValueError('argument out of range')

    def whitney_number(self, k, kind=1):
        r"""
        Return the ``k``-th Whitney number.

        If ``kind=1``, this number is obtained by summing the Moebius function
        values `mu(0, x)` over all `x` in the intersection poset with
        `\mathrm{rank}(x) = k`.

        If ``kind=2``, this number is the number of elements `x, y` in the
        intersection poset such that `x \leq y` with ranks `i` and `j`,
        respectively. 

        See [GZ]_ for more details.

        INPUT:

        - ``k`` -- integer

        - ``kind`` -- 1 or 2 (default: 1)

        OUTPUT:

        Integer. The ``k``-th Whitney number.

        .. SEEALSO::

            :meth:`doubly_indexed_whitney_number`
            :meth:`whitney_data`

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.whitney_number(0)
            1
            sage: A.whitney_number(1)
            -6
            sage: A.whitney_number(2)
            9
            sage: A.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x
            sage: A.whitney_number(1,kind=2)
            6
            sage: p = A.intersection_poset()
            sage: r = p.rank_function()
            sage: len([i for i in p if r(i) == 1])
            6
        """
        if k >= 0 and k <= self.dimension():
            if kind == 1:
                return self.whitney_data()[0][0, k]
            elif kind == 2:
                return self.whitney_data()[1][0, k]
        raise ValueError('argument out of range')

    def is_separating_hyperplane(self, region1, region2, hyperplane):
        r"""
        Test whether the ``hyperplane`` separates the given regions.

        INPUT:

        - ``region1``, ``region2`` -- polyhedra or list/tuple/iterable
          of coordinates which are regions of the arrangement or an interior
          point of a region

        - ``hyperplane`` -- a hyperplane

        OUTPUT:

        A boolean. Whether the hyperplane ``hyperplane`` separate the given
        regions.

        EXAMPLES::

            sage: A.<x,y> = hyperplane_arrangements.coordinate(2)
            sage: A.is_separating_hyperplane([1,1], [2,1], y)
            False
            sage: A.is_separating_hyperplane([1,1], [-1,1], x)
            True
            sage: r = A.region_containing_point([1,1])
            sage: s = A.region_containing_point([-1,1])
            sage: A.is_separating_hyperplane(r, s, x)
            True
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('requires characteristic zero')
        try:
            p1 = region1.representative_point()
        except AttributeError:
            p1 = list(region1)
        try:
            p2 = region2.representative_point()
        except AttributeError:
            p2 = list(region2)
        from sage.functions.generalized import sign
        s = sign(hyperplane(p1)) * sign(hyperplane(p2))
        if s < 0:
            return True
        if s > 0:
            return False
        raise ValueError('point lies on hyperplane')

    def distance_between_regions(self, region1, region2):
        r"""
        Return the number of hyperplanes separating the two regions.

        INPUT:

        - ``region1``, ``region2`` -- regions of the arrangement or
          representative points of regions

        OUTPUT:

        An integer. The number of hyperplanes separating the two regions.

        EXAMPLES::

            sage: c = hyperplane_arrangements.coordinate(2)
            sage: r = c.region_containing_point([-1, -1])
            sage: s = c.region_containing_point([1, 1])
            sage: c.distance_between_regions(r, s)
            2
            sage: c.distance_between_regions(s, s)
            0
        """
        count = sum(1 for hyperplane in self 
                    if self.is_separating_hyperplane(region1, region2, hyperplane))
        return ZZ(count)

    def distance_enumerator(self, base_region):
        r"""
        Return the generating function for the number of hyperplanes
        at given distance.

        INPUT:

        - ``base_region`` -- region of arrangement or point in region

        OUTPUT:

        A polynomial `f(x)` for which the coefficient of `x^i` is the
        number of hyperplanes of distance `i` from ``base_region``,
        i.e., the number of hyperplanes separated by `i` hyperplanes
        from ``base_region``.

        EXAMPLES::

            sage: c = hyperplane_arrangements.coordinate(3)
            sage: c.distance_enumerator(c.region_containing_point([1,1,1]))
            x^3 + 3*x^2 + 3*x + 1
        """
        d = [self.distance_between_regions(r,base_region) for r in self.regions()]
        d = [d.count(i) for i in range(max(d)+1)]
        from sage.rings.polynomial.polynomial_ring import polygen
        x = polygen(QQ, 'x')
        return sum([d[i]*x**i for i in range(len(d))])

    @cached_method
    def varchenko_matrix(self, names='h'):
        r"""
        Return the Varchenko matrix of the arrangement.

        Let `H_1, \ldots, H_s` and `R_1, \ldots, R_t` denote the hyperplanes
        and regions, respectively, of the arrangement.  Let `S =
        \QQ[h_1, \ldots, h_s]`, a polynomial ring with indeterminate `h_i`
        corresponding to hyperplane `H_i`.  The Varchenko matrix is
        the `t \times t` matrix with `i,j`-th entry the product of
        those `h_k` such that `H_k` separates `R_i` and `R_j`.

        INPUT:

        - ``names`` -- string or list/tuple/iterable of strings. The
          variable names for the polynomial ring `S`.

        OUTPUT:

        The Varchenko matrix.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(3)
            sage: v = a.varchenko_matrix();  v
            [    1    h2    h1]
            [   h2     1 h1*h2]
            [   h1 h1*h2     1]
            sage: factor(det(v))
            (h2 - 1) * (h2 + 1) * (h1 - 1) * (h1 + 1)
        """
        from sage.rings.all import PolynomialRing
        from sage.matrix.constructor import identity_matrix
        from sage.misc.all import prod
        k = len(self)
        R = PolynomialRing(QQ, names, k)
        h = R.gens()
        region = self.regions()
        v = identity_matrix(R, k, k)
        for i in range(k):
            for j in range(i+1, k):
                t = prod(h[p] for p in range(k) if
                         self.is_separating_hyperplane(region[i], region[j], self[p]))
                v[i,j] = v[j,i] = t
        v.set_immutable()
        return v



class HyperplaneArrangements(Parent, UniqueRepresentation):
    """
    Hyperplane arrangements.

    For more information on hyperplane arrangements, see
    :mod:`sage.geometry.hyperplane_arrangements.arrangement`.

    INPUT:

    - ``base_ring`` -- ring; the base ring

    - ``names`` -- tuple of strings; the variable names

    EXAMPLES::
    
        sage: H.<x,y> = HyperplaneArrangements(QQ)
        sage: x
        Hyperplane x + 0*y + 0
        sage: x + y
        Hyperplane x + y + 0
        sage: H(x, y, x-1, y-1)
        Arrangement <y - 1 | y | x - 1 | x>
    """
    Element = HyperplaneArrangementElement

    def __init__(self, base_ring, names=tuple()):
        """
        Initialize ``self``.

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: K = HyperplaneArrangements(QQ, names=('x', 'y'))
            sage: H is K
            True
            sage: type(K)
            <class 'sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangements_with_category'>
            sage: K.change_ring(RR).gen(0)
            Hyperplane 1.00000000000000*x + 0.000000000000000*y + 0.000000000000000

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: TestSuite(H).run()
            sage: K = HyperplaneArrangements(QQ)
            sage: TestSuite(K).run()
        """
        from sage.categories.all import Fields, Sets
        if not base_ring in Fields: 
            raise ValueError('base ring must be a field')
        super(HyperplaneArrangements, self).__init__(category=Sets())
        self._base_ring = base_ring
        self._names = names

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        The base ring of the hyperplane arrangement.

        EXAMPLES::

            sage: L.<x,y> = HyperplaneArrangements(QQ)
            sage: L.base_ring()
            Rational Field
        """
        return self._base_ring

    def change_ring(self, base_ring):
        """
        Return hyperplane arrangements over a different base ring.

        INPUT:

        - ``base_ring`` -- a ring; the new base ring.

        OUTPUT:
        
        A new :class:`HyperplaneArrangements` instance over the new
        base ring.

        EXAMPLES::

            sage: L.<x,y> = HyperplaneArrangements(QQ)
            sage: L.gen(0)
            Hyperplane x + 0*y + 0
            sage: L.change_ring(RR).gen(0)
            Hyperplane 1.00000000000000*x + 0.000000000000000*y + 0.000000000000000      

        TESTS::

            sage: L.change_ring(QQ) is L
            True
        """
        return HyperplaneArrangements(base_ring, names=self._names)

    @cached_method
    def ambient_space(self):
        """
        Return the ambient space.

        The ambient space is the parent of hyperplanes. That is, new
        hyperplanes are always constructed internally from the ambient
        space instance.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: L.ambient_space()([(1,0), 0])
            Hyperplane x + 0*y + 0
            sage: L.ambient_space()([(1,0), 0]) == x
            True
        """
        return AmbientVectorSpace(self.base_ring(), self._names)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 2-dimensional linear space over Rational Field with coordinates x, y
        """
        return 'Hyperplane arrangements in {0}'.format(self.ambient_space())

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        INPUT:

        - ``*args`` -- positional arguments, each defining a
          hyperplane; alternatively, a single polytope or a single
          hyperplane arrangement

        - ``signed`` -- boolean (optional, default: ``True``); whether to
          preserve signs of hyperplane equations

        - ``warn_duplicates`` -- boolean (optional, default: ``False``);
          whether to issue a warning if duplicate hyperplanes were
          passed -- note that duplicate hyperplanes are always removed,
          whether or not there is a warning shown

        - ``check`` -- boolean (optional, default: ``True``); whether to
          perform argument checking.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: L._element_constructor_(x, y)
            Arrangement <y | x>
            sage: L._element_constructor_([x, y])
            Arrangement <y | x>
            sage: L._element_constructor_([0, 1, 0], [0, 0, 1])
            Arrangement <y | x>
            sage: L._element_constructor_([[0, 1, 0], [0, 0, 1]])
            Arrangement <y | x>

            sage: L._element_constructor(polytopes.hypercube(2))
            Arrangement <-x + 1 | -y + 1 | y + 1 | x + 1>

            sage: L(x, x, warn_duplicates=True)
            doctest:...: UserWarning: Input contained 2 hyperplanes, but only 1 are distinct.
            Arrangement <x>
            sage: L(-x, x + y - 1, signed=False)
            Arrangement <-x - y + 1 | x>

        TESTS::

            sage: L()
            Empty hyperplane arrangement of dimension 2
            sage: L(0)        # zero is equivalent to no argument, Trac #8648
            Empty hyperplane arrangement of dimension 2
            sage: L(0*x)      # degenerate hyperplane is NOT allowed
            Traceback (most recent call last):
            ...
            ValueError: linear expression must be non-constant to define a hyperplane
            sage: L(0*x, y)   # ditto
            Traceback (most recent call last):
            ...
            ValueError: linear expression must be non-constant to define a hyperplane
       """
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, HyperplaneArrangementElement) and args[0].parent() is self:
                # optimization if argument is already a hyperplane arrangement
                return arg
            if arg == 0 and not isinstance(arg, Hyperplane):
                # zero = neutral element under addition = the empty hyperplane arrangement
                args = []
        # process keyword arguments
        not_char2 = (self.base_ring().characteristic() != 2)
        signed = kwds.pop('signed', not_char2)
        warn_duplicates = kwds.pop('warn_duplicates', False)
        check = kwds.pop('check', True)
        if len(kwds) > 0:
            raise ValueError('unknown keyword argument')
        # process positional arguments
        AA = self.ambient_space()
        try:
            hyperplanes = map(AA, args)
        except (TypeError, ValueError, AttributeError):
            if len(args) > 1:
                raise
            arg = args[0]
            if hasattr(arg, 'Hrepresentation'):
                hyperplanes = [AA(h) for h in arg.Hrepresentation()]
            else:
                hyperplanes = map(AA, arg)
        hyperplanes = [h.primitive(signed) for h in hyperplanes]
        n = len(hyperplanes)
        hyperplanes = tuple(uniq(hyperplanes))
        if warn_duplicates and n != len(hyperplanes):
            from warnings import warn
            warn('Input contained {0} hyperplanes, but only {1} are distinct.'.format(n, len(hyperplanes)))
        # argument checking (optional but recommended)
        if check:
            if signed and not not_char2:
                raise ValueError('cannot be signed in characteristic 2')
            hyperplane_set = set(hyperplanes)
            for h in hyperplanes:
                if h.A() == 0:
                    raise ValueError('linear expression must be non-constant to define a hyperplane')
                if not_char2 and -h in hyperplane_set:
                    raise ValueError('arrangement cannot simultaneouly have h and -h as hyperplane')
        return self.element_class(self, hyperplanes)

    @cached_method
    def ngens(self):
        """
        Return the number of linear variables.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 3-dimensional linear space over Rational Field with coordinates x, y, z
            sage: L.ngens()
            3
        """
        return len(self._names)

    @cached_method
    def gens(self):
        """
        Return the coordinate hyperplanes.
        
        OUTPUT:

        A tuple of linear expressions, one for each linear variable.

        EXAMPLES::

            sage: L = HyperplaneArrangements(QQ, ('x', 'y', 'z'))
            sage: L.gens()
            (Hyperplane x + 0*y + 0*z + 0, 
             Hyperplane 0*x + y + 0*z + 0, 
             Hyperplane 0*x + 0*y + z + 0)
        """
        return self.ambient_space().gens()

    def gen(self, i):
        """
        Return the `i`-th coordinate hyperplane.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        A linear expression.

        EXAMPLES::

            sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 3-dimensional linear space over Rational Field with coordinates x, y, z
            sage: L.gen(0)
            Hyperplane x + 0*y + 0*z + 0
        """
        return self.gens()[i]

    def _coerce_map_from_(self, P):
        """
        Return whether there is a coercion.
   
        TESTS::
        
            sage: L.<x> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
            sage: M.<y> = HyperplaneArrangements(RR);  M
            Hyperplane arrangements in 1-dimensional linear space over Real Field with 53 bits of precision with coordinate y

            sage: L.coerce_map_from(ZZ)
            Conversion map:
              From: Integer Ring
              To:   Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
            sage: M.coerce_map_from(L)
            Conversion map:
              From: Hyperplane arrangements in 1-dimensional linear space over
                    Rational Field with coordinate x
              To:   Hyperplane arrangements in 1-dimensional linear space over
                    Real Field with 53 bits of precision with coordinate y
            sage: L.coerce_map_from(M)
        """
        if self.ambient_space().has_coerce_map_from(P):
            return True
        if isinstance(P, HyperplaneArrangements):
            return self.base_ring().has_coerce_map_from(P.base_ring())
        return super(HyperplaneArrangements, self)._coerce_map_from_(P)

