# sage_setup: distribution = sagemath-polyhedra
r"""
Hyperplane Arrangements

Before talking about hyperplane arrangements, let us start with
individual hyperplanes. This package uses certain linear expressions
to represent hyperplanes, that is, a linear expression `3x + 3y - 5z - 7`
stands for the hyperplane with the equation `3x + 3y - 5z = 7`. To create it
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
    (Hyperplane x + 0*y + 0*z + 0,
     Hyperplane 0*x + y + 0*z + 0,
     Hyperplane 0*x + 0*y + z + 0)

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

Number fields are also possible::

    sage: # needs sage.rings.number_field
    sage: x = polygen(QQ, 'x')
    sage: NF.<a> = NumberField(x**4 - 5*x**2 + 5, embedding=1.90)
    sage: H.<y,z> = HyperplaneArrangements(NF)
    sage: A = H([[(-a**3 + 3*a, -a**2 + 4), 1], [(a**3 - 4*a, -1), 1],
    ....:        [(0, 2*a**2 - 6), 1], [(-a**3 + 4*a, -1), 1],
    ....:        [(a**3 - 3*a, -a**2 + 4), 1]])
    sage: A
    Arrangement of 5 hyperplanes of dimension 2 and rank 2
    sage: A.base_ring()
    Number Field in a with defining polynomial x^4 - 5*x^2 + 5
     with a = 1.902113032590308?

Notation (iii): a list or tuple of hyperplanes::

    sage: H.<x,y,z> = HyperplaneArrangements(GF(5))
    sage: k = [x+i for i in range(4)];  k
    [Hyperplane x + 0*y + 0*z + 0, Hyperplane x + 0*y + 0*z + 1,
     Hyperplane x + 0*y + 0*z + 2, Hyperplane x + 0*y + 0*z + 3]
    sage: H(k)
    Arrangement <x | x + 1 | x + 2 | x + 3>

Notation (iv): using the library of arrangements::

    sage: hyperplane_arrangements.braid(4)                                              # needs sage.graphs
    Arrangement of 6 hyperplanes of dimension 4 and rank 3
    sage: hyperplane_arrangements.semiorder(3)
    Arrangement of 6 hyperplanes of dimension 3 and rank 2
    sage: hyperplane_arrangements.graphical(graphs.PetersenGraph())                     # needs sage.graphs
    Arrangement of 15 hyperplanes of dimension 10 and rank 9
    sage: hyperplane_arrangements.Ish(5)
    Arrangement of 20 hyperplanes of dimension 5 and rank 4

Notation (v): from the bounding hyperplanes of a polyhedron::

    sage: a = polytopes.cube().hyperplane_arrangement();  a
    Arrangement of 6 hyperplanes of dimension 3 and rank 3
    sage: a.n_regions()
    27

New arrangements from old::

    sage: # needs sage.graphs
    sage: a = hyperplane_arrangements.braid(3)
    sage: b = a.add_hyperplane([4, 1, 2, 3])
    sage: b
    Arrangement <t1 - t2 | t0 - t1 | t0 - t2 | t0 + 2*t1 + 3*t2 + 4>
    sage: c = b.deletion([4, 1, 2, 3])
    sage: a == c
    True

    sage: # needs sage.combinat sage.graphs
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

Properties of Arrangements
--------------------------

A hyperplane arrangement is *essential* if the normals to its
hyperplanes span the ambient space.  Otherwise, it is *inessential*.
The essentialization is formed by intersecting the hyperplanes by this
normal space (actually, it is a bit more complicated over finite
fields)::

    sage: # needs sage.graphs
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

    sage: # needs sage.combinat
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

    sage: # needs sage.graphs
    sage: a = hyperplane_arrangements.braid(3)
    sage: p = a.intersection_poset()
    sage: p.is_ranked()
    True
    sage: p.order_polytope()
    A 5-dimensional polyhedron in ZZ^5 defined as the convex hull of 10 vertices

The characteristic polynomial is a basic invariant of a hyperplane
arrangement. It is defined as

.. MATH::

    \chi(x) := \sum_{w\in P} \mu(w) x^{dim(w)}

where `P` is the
:meth:`~HyperplaneArrangementElement.intersection_poset` of the
arrangement and `\mu` is the Möbius function of `P`::

    sage: # long time
    sage: a = hyperplane_arrangements.semiorder(5)
    sage: a.characteristic_polynomial()               # about a second on Core i7
    x^5 - 20*x^4 + 180*x^3 - 790*x^2 + 1380*x
    sage: a.poincare_polynomial()
    1380*x^4 + 790*x^3 + 180*x^2 + 20*x + 1
    sage: a.n_regions()
    2371
    sage: charpoly = a.characteristic_polynomial()
    sage: charpoly(-1)
    -2371
    sage: a.n_bounded_regions()
    751
    sage: charpoly(1)
    751

For finer invariants derived from the intersection poset, see
:meth:`~HyperplaneArrangementElement.whitney_number` and
:meth:`~HyperplaneArrangementElement.doubly_indexed_whitney_number`.

Miscellaneous methods (see documentation for an explanation)::

    sage: a = hyperplane_arrangements.semiorder(3)
    sage: a.has_good_reduction(5)                                                       # needs sage.rings.finite_rings
    True
    sage: b = a.change_ring(GF(5))
    sage: pa = a.intersection_poset()                                                   # needs sage.graphs
    sage: pb = b.intersection_poset()                                                   # needs sage.rings.finite_rings
    sage: pa.is_isomorphic(pb)                                                          # needs sage.graphs sage.rings.finite_rings
    True
    sage: a.face_vector()                                                               # needs sage.graphs
    (0, 12, 30, 19)
    sage: a.face_vector()                                                               # needs sage.graphs
    (0, 12, 30, 19)
    sage: a.is_central()
    False
    sage: a.is_linear()
    False
    sage: a.sign_vector((1,1,1))
    (-1, 1, -1, 1, -1, 1)
    sage: a.varchenko_matrix()[:6, :6]
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
a companion to Richard Stanley's notes [Sta2007]_ on hyperplane
arrangements.
"""

# *****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Possible extensions for hyperplane_arrangement.py:
# - the big face lattice
# - create ties with the Sage matroid methods
# - hyperplane arrangements over other fields

from sage.geometry.hyperplane_arrangement.hyperplane import AmbientVectorSpace, Hyperplane
from sage.matrix.constructor import matrix, vector
from sage.misc.cachefunc import cached_method
from sage.modules.free_module import VectorSpace
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation


class HyperplaneArrangementElement(Element):
    """
    A hyperplane arrangement.

    .. WARNING::

        You should never create
        :class:`HyperplaneArrangementElement` instances directly,
        always use the parent.
    """
    def __init__(self, parent, hyperplanes, check=True, backend=None):
        """
        Construct a hyperplane arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`HyperplaneArrangements`

        - ``hyperplanes`` -- tuple of hyperplanes

        - ``check`` -- boolean (default: ``True``); whether to check input

        - ``backend`` -- string (optional); the backend to
          use for the related polyhedral objects

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement <y | x>
            sage: TestSuite(elt).run()

        It is possible to specify a backend for polyhedral computations::

            sage: # needs sage.rings.number_field
            sage: R.<sqrt5> = QuadraticField(5)
            sage: H = HyperplaneArrangements(R, names='xyz')
            sage: x, y, z = H.gens()
            sage: A = H(sqrt5*x + 2*y + 3*z, backend='normaliz')
            sage: A.backend()
            'normaliz'
            sage: A.regions()[0].backend()                              # optional - pynormaliz
            'normaliz'
        """
        super().__init__(parent)
        self._hyperplanes = hyperplanes
        self._backend = backend
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

        INPUT/OUTPUT: see :meth:`HyperplaneArrangements._first_ngens`

        EXAMPLES::

            sage: a.<x,y,z> = hyperplane_arrangements.braid(3)   # indirect doctest     # needs sage.graphs
            sage: (x, y) == a._first_ngens(2)                                           # needs sage.graphs
            True
        """
        return self.parent()._first_ngens(n)

    def __getitem__(self, i):
        """
        Return the `i`-th hyperplane.

        INPUT:

        - ``i`` -- integer

        OUTPUT: the `i`-th hyperplane

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

    def __hash__(self):
        r"""
        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = x|y; h
            Arrangement <y | x>
            sage: len_dict = {h: len(h)}
        """
        return hash(self.hyperplanes())

    def n_hyperplanes(self):
        r"""
        Return the number of hyperplanes in the arrangement.

        OUTPUT: integer

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
        Return the hyperplanes in the arrangement as a tuple.

        OUTPUT: a tuple

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

        OUTPUT: string

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

        OUTPUT: integer

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

            sage: # needs sage.graphs
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

    def backend(self):
        """
        Return the backend used for polyhedral objects.

        OUTPUT: string giving the backend or ``None`` if none is specified

        EXAMPLES:

        By default, no backend is specified::

           sage: H = HyperplaneArrangements(QQ)
           sage: A = H()
           sage: A.backend()

        Otherwise, one may specify a polyhedral backend::

           sage: A = H(backend='ppl')
           sage: A.backend()
           'ppl'
           sage: A = H(backend='normaliz')
           sage: A.backend()
           'normaliz'
        """
        return self._backend

    def _richcmp_(self, other, op):
        """
        Compare two hyperplane arrangements.

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: H(x) == H(y)
            False

        TESTS::

            sage: H(x) == 0
            False
        """
        return richcmp(self._hyperplanes, other._hyperplanes, op)

    def union(self, other):
        r"""
        The union of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a hyperplane arrangement or something that can
          be converted into a hyperplane arrangement

        OUTPUT: a new hyperplane arrangement

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,2,3], [0,1,1], [0,1,-1], [1,-1,0], [1,1,0])
            sage: B = H([1,1,1], [1,-1,1], [1,0,-1])
            sage: C = A.union(B); C
            Arrangement of 8 hyperplanes of dimension 2 and rank 2
            sage: C == A | B   # syntactic sugar
            True

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
        other_h = P(other)
        hyperplanes = self._hyperplanes + other_h._hyperplanes
        result = P(*hyperplanes, backend=self._backend)
        return result

    add_hyperplane = union

    __or__ = union

    def plot(self, **kwds):
        """
        Plot the hyperplane arrangement.

        OUTPUT: a graphics object

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: L(x, y, x+y-2).plot()                                                 # needs sage.plot
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

        A new hyperplane arrangement `L`.
        Its equations consist of `[0, -d, a_1, \ldots, a_n]` for each
        `[d, a_1, \ldots, a_n]` in the original arrangement and the
        equation `[0, 1, 0, \ldots, 0]` (maybe not in this order).

        .. WARNING::

            While there is an almost-one-to-one correspondence between the
            hyperplanes of ``self`` and those of ``self.cone()``, there is
            no guarantee that the order in which they appear in
            ``self.hyperplanes()`` will match the order in which their
            counterparts in ``self.cone()`` will appear in
            ``self.cone().hyperplanes()``! This warning does not apply
            to ordered hyperplane arrangements.

        EXAMPLES::

            sage: # needs sage.combinat
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
        H = type(P).__base__(P.base_ring(), names=names)
        return H(*hyperplanes, backend=self._backend)

    @cached_method
    def intersection_poset(self, element_label='int'):
        r"""
        Return the intersection poset of the hyperplane arrangement.

        INPUT:

        - ``element_label`` -- (default: ``'int'``) specify how an
          intersection should be represented; must be one of the following:

          * ``'subspace'`` -- as a subspace
          * ``'subset'`` -- as a subset of the defining hyperplanes
          * ``'int'`` -- as an integer

        OUTPUT:

        The poset of non-empty intersections of hyperplanes, with intersections
        represented by integers, subsets of integers or subspaces (see the
        examples for more details).

        EXAMPLES:

        By default, the elements of the poset are the integers from `0` through
        the cardinality of the poset *minus one*. The element labelled `0`
        always corresponds to the ambient vector space, and the hyperplanes
        themselves are labelled `1, 2, \ldots, n`, where `n` is the number
        of hyperplanes of the arrangement. ::

            sage: A = hyperplane_arrangements.coordinate(2)
            sage: L = A.intersection_poset(); L                                         # needs sage.combinat
            Finite poset containing 4 elements
            sage: sorted(L)                                                             # needs sage.combinat
            [0, 1, 2, 3]
            sage: L.level_sets()                                                        # needs sage.combinat
            [[0], [1, 2], [3]]

        ::

            sage: # needs sage.combinat
            sage: A = hyperplane_arrangements.semiorder(3)
            sage: L = A.intersection_poset(); L
            Finite poset containing 19 elements
            sage: sorted(L)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
            sage: [sorted(level_set) for level_set in L.level_sets()]
            [[0], [1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]]

        By passing the argument ``element_label="subset"``, each element of the
        intesection poset is labelled by the set of indices of the hyperplanes
        whose intersection is said element. The index of a hyperplane is its
        index in ``self.hyperplanes()``. ::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: L = A.intersection_poset(element_label='subset')                      # needs sage.combinat
            sage: [sorted(level, key=sorted) for level in L.level_sets()]               # needs sage.combinat
            [[{}],
             [{0}, {1}, {2}, {3}, {4}, {5}],
             [{0, 2}, {0, 3}, {0, 4}, {0, 5}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 4}, {2, 5}, {3, 4}, {3, 5}]]

        ::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H((y, y-1, y+1, x-y, x+y))
            sage: L = A.intersection_poset(element_label='subset')                      # needs sage.combinat
            sage: sorted(L, key=sorted)                                                 # needs sage.combinat
            [{}, {0}, {0, 3}, {0, 4}, {1}, {1, 3, 4}, {2}, {2, 3}, {2, 4}, {3}, {4}]

        One can instead use affine subspaces as elements,
        which is what is used to compute the poset in the first place::

            sage: A = hyperplane_arrangements.coordinate(2)
            sage: L = A.intersection_poset(element_label='subspace'); L                 # needs sage.combinat
            Finite poset containing 4 elements
            sage: sorted(L, key=lambda S: (S.dimension(),                               # needs sage.combinat
            ....:                          S.linear_part().basis_matrix()))
            [Affine space p + W where:
               p = (0, 0)
               W = Vector space of degree 2 and dimension 0 over Rational Field
                   Basis matrix: [],
             Affine space p + W where:
               p = (0, 0)
               W = Vector space of degree 2 and dimension 1 over Rational Field
                   Basis matrix: [0 1],
             Affine space p + W where:
               p = (0, 0)
               W = Vector space of degree 2 and dimension 1 over Rational Field
                   Basis matrix: [1 0],
             Affine space p + W where:
               p = (0, 0)
               W = Vector space of dimension 2 over Rational Field]
        """
        if element_label == "int":
            def update(mapping, val, I0):
                mapping[val] = len(mapping)
        elif element_label == "subset":
            from sage.sets.set import Set

            def update(mapping, val, I0):
                mapping[val] = Set(val)
        elif element_label == "subspace":
            def update(mapping, val, I0):
                mapping[val] = I0
        else:
            raise ValueError("invalid element label type")

        K = self.base_ring()
        from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace

        whole_space = AffineSubspace(0, VectorSpace(K, self.dimension()))
        hyperplanes = [H._affine_subspace() for H in self.hyperplanes()]

        mapping = {}
        update(mapping, frozenset(), whole_space)
        for i, H in enumerate(hyperplanes):
            update(mapping, frozenset([i]), H)

        hasse = {frozenset(): [frozenset([i]) for i in range(len(hyperplanes))]}
        cur_level = [(frozenset([i]), H) for i, H in enumerate(hyperplanes)]

        while cur_level:
            new_level = {}
            for label, T in cur_level:
                edges = []
                for i, H in enumerate(hyperplanes):
                    I0 = H.intersection(T)
                    if I0 is not None and I0 != T:
                        try:
                            target = new_level[I0]
                        except KeyError:
                            target = set(label)
                            new_level[I0] = target
                        target.add(i)
                        edges.append(target)
                hasse[label] = edges
            for label, T in cur_level:
                # Freeze them in place now
                hasse[label] = [frozenset(X) for X in hasse[label]]
            cur_level = [(frozenset(X), T) for T, X in new_level.items()]
            for label, T in cur_level:
                update(mapping, label, T)

        from sage.combinat.posets.posets import Poset
        return Poset({mapping[i]: [mapping[j] for j in val] for i, val in hasse.items()})

    def _slow_characteristic_polynomial(self):
        """
        Return the characteristic polynomial of the hyperplane arrangement.

        This is the slow computation directly from the definition. For
        educational use only.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a._slow_characteristic_polynomial()                                   # needs sage.combinat
            x^2 - 2*x + 1
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        x = polygen(QQ, 'x')
        P = self.intersection_poset()
        n = self.dimension()
        return sum([P.moebius_function(0, p) * x**(n - P.rank(p)) for p in P])

    @cached_method
    def characteristic_polynomial(self):
        r"""
        Return the characteristic polynomial of the hyperplane arrangement.

        OUTPUT: the characteristic polynomial in `\QQ[x]`

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

        Check the corner case of the empty arrangement::

            sage: E = H()
            sage: E.characteristic_polynomial()
            1
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        x = polygen(QQ, 'x')
        if self.rank() == 1:
            return x**(self.dimension() - 1) * (x - len(self))
        if self.rank() == 0:
            return x ** 0

        H = self[0]
        R = self.restriction(H)
        charpoly_R = R.characteristic_polynomial()
        D = self.deletion(H)
        charpoly_D = D.characteristic_polynomial()
        return charpoly_D - charpoly_R

    @cached_method
    def poincare_polynomial(self):
        r"""
        Return the Poincaré polynomial of the hyperplane arrangement.

        OUTPUT: the Poincaré polynomial in `\QQ[x]`

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

    @cached_method
    def cocharacteristic_polynomial(self):
        r"""
        Return the cocharacteristic polynomial of ``self``.

        The cocharacteristic polynomial of a hyperplane arrangement `A`
        is defined by

        .. MATH::

            \Psi_A(z) := \sum_{X \in L} |\mu(B,X)| z^{\dim X},

        where `L` is the intersection poset of `A`, `B` is the minimal
        element of `L` (here, the `0` dimensional subspace), and
        `\mu` is the Möbius function of `L`.

        OUTPUT: the cocharacteristic polynomial in `\ZZ[z]`

        EXAMPLES::

            sage: A = hyperplane_arrangements.coordinate(2)
            sage: A.cocharacteristic_polynomial()                                       # needs sage.graphs
            z^2 + 2*z + 1
            sage: B = hyperplane_arrangements.braid(3)
            sage: B.cocharacteristic_polynomial()                                       # needs sage.graphs
            2*z^3 + 3*z^2 + z

        TESTS::

            sage: I = hyperplane_arrangements.Ish(2)
            sage: I.is_central()
            False
            sage: I.cocharacteristic_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: only defined for central hyperplane arrangements
        """
        if not self.is_central():
            raise ValueError("only defined for central hyperplane arrangements")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(ZZ, 'z')
        z = R.gen()
        L = self.intersection_poset(element_label='subspace').dual()
        B = L.minimal_elements()[0]
        return R.sum(abs(L.moebius_function(B, X)) * z**X.dimension()
                     for X in L)

    @cached_method
    def primitive_eulerian_polynomial(self):
        r"""
        Return the primitive Eulerian polynomial of ``self``.

        The primitive Eulerian polynomial of a hyperplane arrangement `A`
        is defined [BHS2023]_ by

        .. MATH::

            P_A(z) := \sum_{X \in L} |\mu(B,X)| (z - 1)^{\mathrm{codim} X},

        where `L` is the intersection poset of `A`, `B` is the minimal
        element of `L` (here, the `0` dimensional subspace), and
        `\mu` is the Möbius function of `L`.

        OUTPUT: the primitive Eulerian polynomial in `\ZZ[z]`

        EXAMPLES::

            sage: A = hyperplane_arrangements.coordinate(2)
            sage: A.primitive_eulerian_polynomial()                                     # needs sage.graphs
            z^2
            sage: B = hyperplane_arrangements.braid(3)
            sage: B.primitive_eulerian_polynomial()                                     # needs sage.graphs
            z^2 + z

            sage: H = hyperplane_arrangements.Shi(['B',2]).cone()
            sage: H.is_simplicial()
            False
            sage: H.primitive_eulerian_polynomial()                                     # needs sage.graphs
            z^3 + 11*z^2 + 4*z

            sage: H = hyperplane_arrangements.graphical(graphs.CycleGraph(4))
            sage: H.primitive_eulerian_polynomial()                                     # needs sage.graphs
            z^3 + 3*z^2 - z

        We verify Example 2.4 in [BHS2023]_ for `k = 2,3,4,5`::

            sage: R.<x,y> = HyperplaneArrangements(QQ)
            sage: for k in range(2,6):                                                  # needs sage.graphs
            ....:     H = R([x+j*y for j in range(k)])
            ....:     H.primitive_eulerian_polynomial()
            z^2
            z^2 + z
            z^2 + 2*z
            z^2 + 3*z

        We verify Equation (4) in [BHS2023]_ on some examples::

            sage: # needs sage.graphs
            sage: R.<x> = ZZ[]
            sage: Arr = [hyperplane_arrangements.braid(n) for n in range(2,6)]
            sage: all(R(A.cocharacteristic_polynomial()(1/(x-1)) * (x-1)^A.dimension())
            ....:     == R(A.primitive_eulerian_polynomial()) for A in Arr)
            True

        We compute types `H_3` and `F_4` in Table 1 of [BHS2023]_::

            sage: # needs sage.libs.gap
            sage: W = CoxeterGroup(['H',3], implementation='matrix')
            sage: A = HyperplaneArrangements(W.base_ring(), tuple(f'x{s}' for s in range(W.rank())))
            sage: H = A([[0] + list(r) for r in W.positive_roots()])
            sage: H.is_simplicial()                                                     # needs sage.graphs
            True
            sage: H.primitive_eulerian_polynomial()
            z^3 + 28*z^2 + 16*z

            sage: W = CoxeterGroup(['F',4], implementation='permutation')
            sage: A = HyperplaneArrangements(QQ, tuple(f'x{s}' for s in range(W.rank())))
            sage: H = A([[0] + list(r) for r in W.positive_roots()])
            sage: H.primitive_eulerian_polynomial()     # long time                     # needs sage.graphs
            z^4 + 116*z^3 + 220*z^2 + 48*z

        We verify Proposition 2.5 in [BHS2023]_ on the braid arrangement
        `B_k` for `k = 2,3,4,5`::

            sage: B = [hyperplane_arrangements.braid(k) for k in range(2,6)]
            sage: all(H.is_simplicial() for H in B)
            True
            sage: all(c > 0 for H in B                                                  # needs sage.graphs
            ....:     for c in H.primitive_eulerian_polynomial().coefficients())
            True

        We verify Example 9.4 in [BHS2023]_ showing a hyperplane arrangement
        whose primitive Eulerian polynomial does not have real roots (in
        general, the graphical arrangement of a cycle graph corresponds
        to the arrangements in Example 9.4)::

            sage: # needs sage.graphs
            sage: H = hyperplane_arrangements.graphical(graphs.CycleGraph(5))
            sage: pep = H.primitive_eulerian_polynomial(); pep
            z^4 + 6*z^3 - 4*z^2 + z
            sage: pep.roots(QQbar)
            [(-6.626418492719221?, 1),
             (0, 1),
             (0.3132092463596102? - 0.2298065541510677?*I, 1),
             (0.3132092463596102? + 0.2298065541510677?*I, 1)]
            sage: pep.roots(AA)
            [(-6.626418492719221?, 1), (0, 1)]

        TESTS::

            sage: I = hyperplane_arrangements.Ish(2)
            sage: I.is_central()
            False
            sage: I.primitive_eulerian_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: only defined for central hyperplane arrangements
        """
        if not self.is_central():
            raise ValueError("only defined for central hyperplane arrangements")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(ZZ, 'z')
        z = R.gen()
        L = self.intersection_poset(element_label='subspace').dual()
        B = L.minimal_elements()[0]
        n = self.dimension()
        return R.sum(abs(L.moebius_function(B, X)) * (z - 1)**(n-X.dimension())
                     for X in L)

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

        Checks that deletion preserves the backend::

            sage: H = HyperplaneArrangements(QQ, names='xyz')
            sage: x,y,z = H.gens()
            sage: h1,h2 = [1*x+2*y+3*z, 3*x+2*y+1*z]
            sage: A = H(h1,h2,backend='normaliz')
            sage: A.deletion(h2).backend()
            'normaliz'
        """
        parent = self.parent()
        hyperplanes = parent(hyperplanes)
        planes = list(self)
        for hyperplane in hyperplanes:
            try:
                planes.remove(hyperplane)
            except ValueError:
                raise ValueError('hyperplane is not in the arrangement')
        return parent(*planes, backend=self._backend)

    def restriction(self, hyperplane, repetitions=False):
        r"""
        Return the restriction to a hyperplane.

        INPUT:

        - ``hyperplane`` -- a hyperplane of the hyperplane arrangement

        - ``repetitions`` -- boolean (default: ``False``); eliminate
          repetitions for ordered arrangements

        OUTPUT:

        The restriction `\mathcal{A}_H` of the
        hyperplane arrangement `\mathcal{A}` to the given ``hyperplane`` `H`.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: A.<u,x,y,z> = hyperplane_arrangements.braid(4);  A
            Arrangement of 6 hyperplanes of dimension 4 and rank 3
            sage: H = A[0];  H
            Hyperplane 0*u + 0*x + y - z + 0
            sage: R = A.restriction(H); R
            Arrangement <x - z | u - x | u - z>
            sage: A.add_hyperplane(z).restriction(z)
            Arrangement of 6 hyperplanes of dimension 3 and rank 3
            sage: A.add_hyperplane(u).restriction(u)
            Arrangement of 6 hyperplanes of dimension 3 and rank 3
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

        TESTS:

        Checks that restriction preserves the backend::

            sage: H = HyperplaneArrangements(QQ, names='xyz')
            sage: x,y,z = H.gens()
            sage: h1,h2 = [1*x+2*y+3*z, 3*x+2*y+1*z]
            sage: A = H(h1, h2, backend='normaliz')
            sage: A.restriction(h2).backend()
            'normaliz'
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
        from sage.geometry.hyperplane_arrangement.ordered_arrangement import OrderedHyperplaneArrangements
        if isinstance(parent, OrderedHyperplaneArrangements):
            H = OrderedHyperplaneArrangements(parent.base_ring(), names=tuple(names))
            if not repetitions:
                L = list(hyperplanes)
                hyperplanes = ()
                for h in L:
                    if h not in hyperplanes:
                        hyperplanes += (h,)
        else:
            H = HyperplaneArrangements(parent.base_ring(), names=tuple(names))
        result = H(*hyperplanes, signed=False, backend=self._backend)
        return result

    def change_ring(self, base_ring):
        """
        Return hyperplane arrangement over the new base ring.

        INPUT:

        - ``base_ring`` -- the new base ring; must be a field for
          hyperplane arrangements

        OUTPUT:

        The hyperplane arrangement obtained by changing the base
        field, as a new hyperplane arrangement.

        .. WARNING::

            While there is often a one-to-one correspondence between the
            hyperplanes of ``self`` and those of
            ``self.change_ring(base_ring)``, there is
            no guarantee that the order in which they appear in
            ``self.hyperplanes()`` will match the order in which their
            counterparts in ``self.cone()`` will appear in
            ``self.change_ring(base_ring).hyperplanes()``!

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,1), 0], [(2,3), -1])
            sage: A.change_ring(FiniteField(2))
            Arrangement <y + 1 | x + y>

        TESTS:

        Checks that changing the ring preserves the backend::

            sage: H = HyperplaneArrangements(QQ, names='xyz')
            sage: x,y,z = H.gens()
            sage: h1, h2 = [1*x+2*y+3*z, 3*x+2*y+1*z]
            sage: A = H(h1, h2, backend='normaliz')
            sage: A.change_ring(RDF).backend()
            'normaliz'
        """
        parent = self.parent().change_ring(base_ring)
        return parent(self, backend=self._backend)

    @cached_method
    def n_regions(self):
        r"""
        The number of regions of the hyperplane arrangement.

        OUTPUT: integer

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

        Check that :issue:`30749` is fixed::

            sage: # needs sage.rings.number_field
            sage: R.<y> = QQ[]
            sage: v1 = AA.polynomial_root(AA.common_polynomial(y^2 - 3),
            ....:                         RIF(RR(1.7320508075688772), RR(1.7320508075688774)))
            sage: v2 = QQbar.polynomial_root(AA.common_polynomial(y^4 - y^2 + 1),
            ....:                            CIF(RIF(RR(0.8660254037844386), RR(0.86602540378443871)),
            ....:                                RIF(-RR(0.50000000000000011), -RR(0.49999999999999994))))
            sage: my_vectors = (vector(AA, [-v1, -1, 1]), vector(AA, [0, 2, 1]), vector(AA, [v1, -1, 1]),
            ....:               vector(AA, [1, 0, 0]), vector(AA, [1/2, AA(-1/2*v2^3 + v2),0]),
            ....:               vector(AA, [-1/2, AA(-1/2*v2^3 + v2), 0]))
            sage: H = HyperplaneArrangements(AA, names='xyz')
            sage: x,y,z = H.gens()
            sage: A = H(backend='normaliz')                                     # optional - pynormaliz
            sage: for v in my_vectors:                                          # optional - pynormaliz
            ....:     a, b, c = v
            ....:     A = A.add_hyperplane(a*x + b*y + c*z)
            sage: A.n_regions()                                                 # optional - pynormaliz
            24
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

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.combinat
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
        from sage.rings.finite_rings.finite_field_constructor import GF
        a = self.change_ring(GF(p))
        p = self.intersection_poset()
        q = a.intersection_poset()
        return p.is_isomorphic(q)

    def is_linear(self):
        r"""
        Test whether all hyperplanes pass through the origin.

        OUTPUT: boolean

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a.is_linear()
            False
            sage: b = hyperplane_arrangements.braid(3)                                  # needs sage.graphs
            sage: b.is_linear()                                                         # needs sage.graphs
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

        OUTPUT: boolean

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: H(x, x+1).is_essential()
            False
            sage: H(x, y).is_essential()
            True
        """
        return self.rank() == self.dimension()

    @cached_method
    def is_central(self, certificate=False):
        r"""
        Test whether the intersection of all the hyperplanes is nonempty.

        A hyperplane arrangement is central if the intersection of all the
        hyperplanes in the arrangement is nonempty.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return the center as a polyhedron (possibly empty) as part
          of the output

        OUTPUT: if ``certificate`` is ``True``, returns a tuple containing:

        1. A boolean
        2. The polyhedron defined to be the intersection of all the hyperplanes

        If ``certificate`` is ``False``, returns a boolean.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)                                  # needs sage.graphs
            sage: a.is_central()                                                        # needs sage.graphs
            True

        The Catalan arrangement in dimension 3 is not central::

            sage: b = hyperplane_arrangements.Catalan(3)
            sage: b.is_central(certificate=True)
            (False, The empty polyhedron in QQ^3)

        The empty arrangement in dimension 5 is central::

            sage: H = HyperplaneArrangements(QQ, names=tuple(['x'+str(i) for i in range(7)]))
            sage: c = H()
            sage: c.is_central(certificate=True)
            (True, A 7-dimensional polyhedron in QQ^7 defined
                   as the convex hull of 1 vertex and 7 lines)
        """
        R = self.base_ring()
        # If there are no hyperplanes in the arrangement,
        # the center is the entire ambient space
        if self.n_hyperplanes() == 0:
            if certificate:
                from sage.geometry.polyhedron.parent import Polyhedra
                pp = Polyhedra(R, self.dimension(), backend=self._backend)
                return (True, pp.universe())
            else:
                return True
        # The center is the set of points contained in all hyperplanes,
        # expressible as the solution set of m*x=b with m and b as follows:
        m = matrix(R, [h.normal() for h in self])
        b = vector(R, [h.b() for h in self])
        try:
            x = m.solve_right(b)
        except ValueError:
            # The solution set is empty, therefore the center is empty
            if certificate:
                from sage.geometry.polyhedron.parent import Polyhedra
                pp = Polyhedra(R, self.dimension(), backend=self._backend)
                return (False, pp.empty())
            else:
                return False
        # The center is the kernel of m translated by x.
        if certificate:
            Ker = m.right_kernel()
            from sage.geometry.polyhedron.constructor import Polyhedron
            return (True, Polyhedron(base_ring=R, vertices=[x],
                                     lines=Ker.basis(),
                                     backend=self._backend))
        else:
            return True

    def center(self):
        r"""
        Return the center of the hyperplane arrangement.

        The polyhedron defined to be the set of all points in the
        ambient space of the arrangement that lie on all of the
        hyperplanes.

        OUTPUT: a polyhedron

        EXAMPLES:

        The empty hyperplane arrangement has the entire ambient space as its
        center::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H()
            sage: A.center()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines

        The Shi arrangement in dimension 3 has an empty center::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.center()
            The empty polyhedron in QQ^3

        The Braid arrangement in dimension 3 has a center that is neither
        empty nor full-dimensional::

            sage: A = hyperplane_arrangements.braid(3)                                  # needs sage.combinat
            sage: A.center()                                                            # needs sage.combinat
            A 1-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex and 1 line
        """
        return self.is_central(certificate=True)[1]

    @cached_method
    def is_simplicial(self):
        r"""
        Test whether the arrangement is simplicial.

        A region is simplicial if the normal vectors of its bounding hyperplanes
        are linearly independent. A hyperplane arrangement is said to be
        simplicial if every region is simplicial.

        OUTPUT: boolean; whether the hyperplane arrangement is simplicial

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H([[0,1,1,1], [0,1,2,3]])
            sage: A.is_simplicial()
            True
            sage: A = H([[0,1,1,1], [0,1,2,3], [0,1,3,2]])
            sage: A.is_simplicial()
            True
            sage: A = H([[0,1,1,1], [0,1,2,3], [0,1,3,2], [0,2,1,3]])
            sage: A.is_simplicial()
            False
            sage: hyperplane_arrangements.braid(3).is_simplicial()                      # needs sage.graphs
            True
        """
        # if the arr is not essential, grab the essential version and check there.
        if not self.is_essential():
            return self.essentialization().is_simplicial()

        # Check that the number of facets for each region is equal to rank
        rank = self.rank()
        return all(R.n_facets() == rank for R in self.regions())

    @cached_method
    def essentialization(self):
        r"""
        Return the essentialization of the hyperplane arrangement.

        The essentialization of a hyperplane arrangement whose base field
        has characteristic 0 is obtained by intersecting the hyperplanes by
        the space spanned by their normal vectors.

        OUTPUT:

        The essentialization `\mathcal{A}'` of `\mathcal{A}` as a
        new hyperplane arrangement.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(3)                                  # needs sage.graphs
            sage: a.is_essential()                                                      # needs sage.graphs
            False
            sage: a.essentialization()                                                  # needs sage.graphs
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
                    return
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
        echelon_pivots = []   # the column indices where N has 1s from the echelonization
        for pivot, row in echelon_col_iter(complement_basis):
            assert row[pivot] == 1
            echelon_pivots.append(pivot)
            hyperplanes = [h - h.A()[pivot] * H(row, 0) for h in hyperplanes]
        # eliminate the pivot'ed coordinates
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
        return restricted_parent(*restricted, signed=False, backend=self._backend)

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
        signs = vector(ZZ, [sign(_) for _ in values])
        signs.set_immutable()
        return signs

    def face_vector(self):
        r"""
        Return the face vector.

        OUTPUT: a vector of integers

        The `d`-th entry is the number of faces of dimension `d`.  A
        *face* is the intersection of a region with a hyperplane of
        the arrangement.

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.face_vector()                                                       # needs sage.combinat
            (0, 6, 21, 16)
        """
        m = self.whitney_data()[0]
        v = list(sum(m.transpose().apply_map(abs)))
        v.reverse()
        v = vector(ZZ, [0]*(self.dimension() - self.rank()) + v)
        v.set_immutable()
        return v

    @cached_method
    def _parallel_hyperplanes(self) -> tuple:
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
        parallels = {}
        for hyperplane in self:
            through_origin = V([list(hyperplane.A()), 0]).primitive(signed=False)
            parallel_planes = parallels.get(through_origin, [])
            A = through_origin.A()
            b = hyperplane.b() * (A / hyperplane.A())
            parallel_planes.append([b, (hyperplane, A, b)])
            parallels[through_origin] = parallel_planes
        parallels = sorted(tuple(hyperplane[1] for hyperplane in sorted(value))
                           for key, value in parallels.items())
        return tuple(parallels)

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

            sage: # needs sage.combinat
            sage: A = hyperplane_arrangements.Shi(3).essentialization()
            sage: A.dimension()
            2
            sage: A.face_vector()
            (6, 21, 16)
            sage: A.vertices()
            ((-2/3, 1/3), (-1/3, -1/3), (0, -1), (0, 0), (1/3, -2/3), (2/3, -1/3))
            sage: point2d(A.vertices(), size=20) + A.plot()                             # needs sage.plot
            Graphics object consisting of 7 graphics primitives

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: chessboard = []
            sage: N = 8
            sage: for x0 in range(N + 1):
            ....:     for y0 in range(N + 1):
            ....:         chessboard.extend([x-x0, y-y0])
            sage: chessboard = H(chessboard)
            sage: len(chessboard.vertices())
            81
            sage: chessboard.vertices(exclude_sandwiched=True)
            ((0, 0), (0, 8), (8, 0), (8, 8))
        """
        import itertools
        from sage.matroids.constructor import Matroid
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
            b_list_list = [skip(_) for _ in b_list_list]
        M = Matroid(groundset=range(len(parallels)), matrix=matrix(A_list).transpose())
        d = self.dimension()
        # vertices are solutions v * lhs = rhs
        lhs = matrix(R, d, d)
        rhs = vector(R, d)
        vertices = set()
        for indices in M.independent_sets(d):
            for row, i in enumerate(indices):
                lhs[row] = A_list[i]
            b_list = [b_list_list[i] for i in indices]
            for b in itertools.product(*b_list):
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

        - ``hyperplanes`` -- list/tuple/iterable of hyperplanes

        OUTPUT:

        The polyhedron constructed from taking the linear expressions
        as inequalities.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = H(x)
            sage: h._make_region([x, 1-x, y, 1-y])
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        TESTS:

        Checks that it creates the regions with the appropriate backend::

            sage: h = H(x,backend='normaliz')
            sage: h._make_region([x, 1-x, y, 1-y]).backend()            # optional - pynormaliz
            'normaliz'
        """
        ieqs = [h.dense_coefficient_list() for h in hyperplanes]
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(ieqs=ieqs, ambient_dim=self.dimension(),
                          base_ring=self.parent().base_ring(),
                          backend=self._backend)

    @cached_method
    def regions(self):
        r"""
        Return the regions of the hyperplane arrangement.

        The base field must have characteristic zero.

        OUTPUT: a tuple containing the regions as polyhedra

        The regions are the connected components of the complement of
        the union of the hyperplanes as a subset of `\RR^n`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)                                  # needs sage.graphs
            sage: a.regions()                                                           # needs sage.graphs
            (A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex, 1 ray, 1 line,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex, 1 ray, 1 line)

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H(x, y+1)
            sage: A.regions()
            (A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays)

            sage: chessboard = []
            sage: N = 8
            sage: for x0 in range(N + 1):
            ....:     for y0 in range(N + 1):
            ....:         chessboard.extend([x-x0, y-y0])
            sage: chessboard = H(chessboard)
            sage: len(chessboard.bounded_regions())   # long time, 359 ms on a Core i7
            64

        Example 6 of [KP2020]_::

            sage: from itertools import product
            sage: def zero_one(d):
            ....:     for x in product([0,1], repeat=d):
            ....:         if any(x):
            ....:             yield [0] + list(x)

            sage: K.<x,y> = HyperplaneArrangements(QQ)
            sage: A = K(*zero_one(2))
            sage: len(A.regions())
            6
            sage: K.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = K(*zero_one(3))
            sage: len(A.regions())
            32
            sage: K.<x,y,z,w> = HyperplaneArrangements(QQ)
            sage: A = K(*zero_one(4))
            sage: len(A.regions())
            370
            sage: K.<x,y,z,w,r> = HyperplaneArrangements(QQ)
            sage: A = K(*zero_one(5))
            sage: len(A.regions())            # not tested (~25s)
            11292

        It is possible to specify the backend::

            sage: # needs sage.rings.number_field
            sage: K.<q> = CyclotomicField(9)
            sage: L.<r9> = NumberField((q + q**(-1)).minpoly(),
            ....:                      embedding=AA(q + q**-1))
            sage: norms = [[1, 1/3*(-2*r9**2-r9+1), 0],
            ....:          [1, -r9**2 - r9, 0],
            ....:          [1, -r9**2 + 1, 0],
            ....:          [1, -r9**2, 0],
            ....:          [1, r9**2 - 4, -r9**2+3]]
            sage: H.<x,y,z> = HyperplaneArrangements(L)
            sage: A = H(backend='normaliz')
            sage: for v in norms:
            ....:     a,b,c = v
            ....:     A = A.add_hyperplane(a*x + b*y + c*z)
            sage: R = A.regions()                                       # optional - pynormaliz
            sage: R[0].backend()                                        # optional - pynormaliz
            'normaliz'

        TESTS::

            sage: K.<x,y,z,w,r> = HyperplaneArrangements(QQ)
            sage: A = K()
            sage: A.regions()
            (A 5-dimensional polyhedron in QQ^5
                 defined as the convex hull of 1 vertex and 5 lines,)
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('base field must have characteristic zero')
        from sage.geometry.polyhedron.constructor import Polyhedron
        R = self.base_ring()
        dim = self.dimension()
        be = self._backend
        universe = Polyhedron(eqns=[[0] + [0] * dim],
                              base_ring=R,
                              backend=be)
        regions = [universe]
        if self.is_linear() and self.n_hyperplanes():
            # We only take the positive half w.r. to the first hyperplane.
            # We fix this by appending all negative regions in the end.
            regions = None

        for hyperplane in self:
            ieq = vector(R, hyperplane.dense_coefficient_list())
            pos_half = Polyhedron(ieqs=[ieq], base_ring=R, backend=be)
            neg_half = Polyhedron(ieqs=[-ieq], base_ring=R, backend=be)
            if not regions:
                # See comment above.
                regions = [pos_half]
                continue
            subdivided = []
            for region in regions:
                # For each region we determine, if the hyperplane splits it.
                splits = False

                # Determine if all vertices lie on one side of the hyperplane.
                # If so, we determine on which side.
                valuations = tuple(ieq[0] + ieq[1:]*v[:] for v in region.vertices())
                direction = 0
                if any(x > 0 for x in valuations):
                    direction = 1
                if any(x < 0 for x in valuations):
                    if direction:
                        splits = True
                    else:
                        direction = -1

                if not splits:
                    # All vertices lie in one closed halfspace of the hyperplane.
                    region_lines = region.lines()
                    if direction == 0:
                        # In this case all vertices lie on the hyperplane and we must
                        # check if rays are contained in one closed halfspace given by the hyperplane.
                        valuations = tuple(ieq[1:]*ray[:] for ray in region.rays())
                        if region_lines:
                            valuations += tuple(ieq[1:]*line[:] for line in region_lines)
                            valuations += tuple(-ieq[1:]*line[:] for line in region_lines)
                        if any(x > 0 for x in valuations) and any(x < 0 for x in valuations):
                            splits = True
                    else:
                        # In this case, at least one of the vertices is not on the hyperplane.
                        # So we check if any ray or line pokes the hyperplane.
                        if (any(ieq[1:]*r[:]*direction < 0 for r in region.rays()) or
                                any(ieq[1:]*ll[:] != 0 for ll in region_lines)):
                            splits = True

                if splits:
                    subdivided.append(region.intersection(pos_half))
                    subdivided.append(region.intersection(neg_half))
                else:
                    subdivided.append(region)
            regions = subdivided

        if self.is_linear() and self.n_hyperplanes():
            # We have treated so far only the positive half space w.r. to the first hyperplane.
            return tuple(regions) + tuple(-x for x in regions)
        else:
            return tuple(regions)

    @cached_method
    def poset_of_regions(self, B=None, numbered_labels=True):
        r"""
        Return the poset of regions for a central hyperplane arrangement.

        The poset of regions is a partial order on the set of regions
        where the regions are ordered by `R\leq R'` if and only if
        `S(R) \subseteq S(R')` where `S(R)` is the set of hyperplanes which
        separate the region `R` from the base region `B`.

        INPUT:

        - ``B`` -- a region (optional); if ``None``, then
          an arbitrary region is chosen as the base region

        - ``numbered_labels`` -- boolean (default: ``True``); if ``True``,
          then the elements of the poset are numbered. Else they are labelled
          with the regions themselves.

        OUTPUT: a Poset object containing the poset of regions

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H([[0,1,1,1], [0,1,2,3]])
            sage: A.poset_of_regions()                                                  # needs sage.combinat
            Finite poset containing 4 elements

            sage: # needs sage.combinat sage.graphs
            sage: A = hyperplane_arrangements.braid(3)
            sage: A.poset_of_regions()
            Finite poset containing 6 elements
            sage: A.poset_of_regions(numbered_labels=False)
            Finite poset containing 6 elements
            sage: A = hyperplane_arrangements.braid(4)
            sage: A.poset_of_regions()
            Finite poset containing 24 elements

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H([[0,1,1,1], [0,1,2,3], [0,1,3,2], [0,2,1,3]])
            sage: R = A.regions()
            sage: base_region = R[3]
            sage: A.poset_of_regions(B=base_region)                                     # needs sage.combinat
            Finite poset containing 14 elements
        """
        from sage.combinat.posets.posets import Poset

        # We use RX to keep track of indexes and R to keep track of which regions
        # we've already hit. This poset is graded, so we can go one set at a time
        RX = self.regions()
        R = set(RX)
        if B in R:
            R.discard(B)
        else:
            B = R.pop()

        # Will record the edges in our poset
        edges = []

        # Start with rank=0 for the poset
        nextTest = [B]

        # While we have objects in our set R
        while R:
            # Transfer the "next step" to the "current step"
            curTest = list(nextTest)
            nextTest = set([])
            # we want to test each region that we haven't hit yet
            for r in R:
                # Since it's graded, it suffices to look at the regions of the previous rank
                for b in curTest:
                    if self.distance_between_regions(b, r) == 1:
                        nextTest.add(r)
                        if numbered_labels:
                            edges.append([RX.index(b), RX.index(r)])
                        else:
                            edges.append([b, r])
            for x in nextTest:
                R.discard(x)

        if numbered_labels:
            return Poset([range(len(RX)), edges])
        else:
            return Poset([RX, edges])

    @cached_method
    def closed_faces(self, labelled=True):
        r"""
        Return the closed faces of the hyperplane arrangement ``self``
        (provided that ``self`` is defined over a totally ordered field).

        Let `\mathcal{A}` be a hyperplane arrangement in the vector
        space `K^n`, whose hyperplanes are the zero sets of the
        affine-linear functions `u_1, u_2, \ldots, u_N`. (We consider
        these functions `u_1, u_2, \ldots, u_N`, and not just the
        hyperplanes, as given. We also assume the field `K` to be
        totally ordered.) For any point `x \in K^n`, we define the
        *sign vector* of `x` to be the vector
        `(v_1, v_2, \ldots, v_N) \in \{-1, 0, 1\}^N` such that (for each
        `i`) the number `v_i` is the sign of `u_i(x)`. For any
        `v \in \{-1, 0, 1\}^N`, we let `F_v` be the set of all `x \in K^n`
        which have sign vector `v`. The nonempty ones among all these
        subsets `F_v` are called the *open faces* of `\mathcal{A}`. They
        form a partition of the set `K^n`.

        Furthermore, for any
        `v = (v_1, v_2, \ldots, v_N) \in \{-1, 0, 1\}^N`, we let `G_v` be
        the set of all `x \in K^n` such that, for every `i`, the sign of
        `u_i(x)` is either `0` or `v_i`.
        Then, `G_v` is a polyhedron. The nonempty ones among all these
        polyhedra `G_v` are called the *closed faces* of `\mathcal{A}`.
        While several sign vectors `v` can lead to one and the same
        closed face `G_v`, we can assign to every closed face a canonical
        choice of a sign vector: Namely, if `G` is a closed face of
        `\mathcal{A}`, then the *sign vector* of `G` is defined to be the
        vector `(v_1, v_2, \ldots, v_N) \in \{-1, 0, 1\}^N` where `x` is
        any point in the relative interior of `G` and where, for each `i`,
        the number `v_i` is the sign of `u_i(x)`. (This does not depend on
        the choice of `x`.)

        There is a one-to-one correspondence between the closed faces and
        the open faces of `\mathcal{A}`. It sends a closed face `G` to
        the open face `F_v`, where `v` is the sign vector of `G`; this
        `F_v` is also the relative interior of `G_v`. The inverse map
        sends any open face `O` to the closure of `O`.

        INPUT:

        - ``labelled`` -- boolean (default: ``True``); if ``True``, then
          this method returns not the faces itself but rather pairs
          `(v, F)` where `F` is a closed face and `v` is its sign vector
          (here, the order and the orientation of the
          `u_1, u_2, \ldots, u_N` is as given by ``self.hyperplanes()``).

        OUTPUT:

        A tuple containing the closed faces as polyhedra, or (if
        ``labelled`` is set to ``True``) the pairs of sign vectors and
        corresponding closed faces.

        .. TODO::

            Should the output rather be a dictionary where the keys are
            the sign vectors and the values are the faces?

        EXAMPLES::

            sage: # needs sage.graphs
            sage: a = hyperplane_arrangements.braid(2)
            sage: a.hyperplanes()
            (Hyperplane t0 - t1 + 0,)
            sage: a.closed_faces()
            (((0,),  A 1-dimensional polyhedron in QQ^2 defined
                     as the convex hull of 1 vertex and 1 line),
             ((1,),  A 2-dimensional polyhedron in QQ^2 defined
                     as the convex hull of 1 vertex, 1 ray, 1 line),
             ((-1,), A 2-dimensional polyhedron in QQ^2 defined
                     as the convex hull of 1 vertex, 1 ray, 1 line))
            sage: a.closed_faces(labelled=False)
            (A 1-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 1 line,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex, 1 ray, 1 line,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex, 1 ray, 1 line)
            sage: [(v, F, F.representative_point()) for v, F in a.closed_faces()]
            [((0,),  A 1-dimensional polyhedron in QQ^2 defined
                     as the convex hull of 1 vertex and 1 line,      (0, 0)),
             ((1,),  A 2-dimensional polyhedron in QQ^2 defined
                     as the convex hull of 1 vertex, 1 ray, 1 line,  (0, -1)),
             ((-1,), A 2-dimensional polyhedron in QQ^2 defined
                     as the convex hull of 1 vertex, 1 ray, 1 line,  (-1, 0))]

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: a = H(x, y+1)
            sage: a.hyperplanes()
            (Hyperplane 0*x + y + 1, Hyperplane x + 0*y + 0)
            sage: [(v, F, F.representative_point()) for v, F in a.closed_faces()]
            [((0, 0),   A 0-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex,             (0, -1)),
             ((0, 1),   A 1-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 1 ray,   (1, -1)),
             ((0, -1),  A 1-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 1 ray,   (-1, -1)),
             ((1, 0),   A 1-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 1 ray,   (0, 0)),
             ((1, 1),   A 2-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 2 rays,  (1, 0)),
             ((1, -1),  A 2-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 2 rays,  (-1, 0)),
             ((-1, 0),  A 1-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 1 ray,   (0, -2)),
             ((-1, 1),  A 2-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 2 rays,  (1, -2)),
             ((-1, -1), A 2-dimensional polyhedron in QQ^2 defined
                        as the convex hull of 1 vertex and 2 rays,  (-1, -2))]

            sage: # needs sage.graphs
            sage: a = hyperplane_arrangements.braid(3)
            sage: a.hyperplanes()
            (Hyperplane 0*t0 + t1 - t2 + 0,
             Hyperplane t0 - t1 + 0*t2 + 0,
             Hyperplane t0 + 0*t1 - t2 + 0)
            sage: [(v, F, F.representative_point()) for v, F in a.closed_faces()]
            [((0, 0, 0),    A 1-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex and 1 line,      (0, 0, 0)),
             ((0, 1, 1),    A 2-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 1 ray, 1 line,  (0, -1, -1)),
             ((0, -1, -1),  A 2-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 1 ray, 1 line,  (-1, 0, 0)),
             ((1, 0, 1),    A 2-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 1 ray, 1 line,  (1, 1, 0)),
             ((1, 1, 1),    A 3-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 2 rays, 1 line, (0, -1, -2)),
             ((1, -1, 0),   A 2-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 1 ray, 1 line,  (-1, 0, -1)),
             ((1, -1, 1),   A 3-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 2 rays, 1 line, (1, 2, 0)),
             ((1, -1, -1),  A 3-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 2 rays, 1 line, (-2, 0, -1)),
             ((-1, 0, -1),  A 2-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 1 ray, 1 line,  (0, 0, 1)),
             ((-1, 1, 0),   A 2-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 1 ray, 1 line,  (1, 0, 1)),
             ((-1, 1, 1),   A 3-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 2 rays, 1 line, (0, -2, -1)),
             ((-1, 1, -1),  A 3-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 2 rays, 1 line, (1, 0, 2)),
             ((-1, -1, -1), A 3-dimensional polyhedron in QQ^3 defined
                            as the convex hull of 1 vertex, 2 rays, 1 line, (-1, 0, 1))]

        Let us check that the number of closed faces with a given
        dimension computed using ``self.closed_faces()`` equals the one
        computed using :meth:`face_vector`::

            sage: def test_number(a):
            ....:     Qx = PolynomialRing(QQ, 'x'); x = Qx.gen()
            ....:     RHS = Qx.sum(vi * x ** i for i, vi in enumerate(a.face_vector()))
            ....:     LHS = Qx.sum(x ** F[1].dim() for F in a.closed_faces())
            ....:     return LHS == RHS
            sage: a = hyperplane_arrangements.Catalan(2)
            sage: test_number(a)                                                        # needs sage.combinat
            True
            sage: a = hyperplane_arrangements.Shi(3)
            sage: test_number(a)                # long time                             # needs sage.combinat
            True

        TESTS:

        An empty border case::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: a = H()
            sage: a.closed_faces()
            (((),
              A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines),)
        """
        R = self.base_ring()
        if R.characteristic() != 0:
            raise ValueError('base field must have characteristic zero')
        from sage.geometry.polyhedron.constructor import Polyhedron
        dim = self.dimension()
        hypes = self.hyperplanes()
        be = self._backend
        universe = Polyhedron(eqns=[[0] + [0] * dim], base_ring=R, backend=be)
        faces = [((), universe)]
        for k, hyperplane in enumerate(hypes):
            # Loop invariant:
            # ``faces == Hk.closed_faces()``, where ``Hk`` is the
            # hyperplane arrangement given by the first ``k`` hyperplanes
            # in the list ``hypes`` (that is, by ``hypes[:k]``).
            ieq = vector(R, hyperplane.dense_coefficient_list())
            zero_half = Polyhedron(eqns=[ieq], base_ring=R, backend=be)
            # ``zero_half`` is the hyperplane ``hyperplane`` itself
            # (viewed as a polyhedron).
            pos_half = Polyhedron(ieqs=[ieq], base_ring=R, backend=be)
            neg_half = Polyhedron(ieqs=[-ieq], base_ring=R, backend=be)
            subdivided = []
            for signs, face in faces:
                # So ``face`` is a face of the hyperplane arrangement
                # given by the first ``k`` hyperplanes in the list
                # ``hypes``, and ``signs`` is the corresponding
                # (length-``k``) sign vector.
                face_dim = face.dim()
                # Adding the intersection of ``face`` with ``hyperplane``:
                zero_part = face.intersection(zero_half)
                zero_part_dim = zero_part.dim()
                if zero_part_dim == face_dim:
                    # If the intersection of ``face`` with ``hyperplane``
                    # has the same dimension as ``face``, then this
                    # intersection *is* ``face``, so we can continue
                    # (without adding the other two intersections, since
                    # those are empty):
                    subdivided.append((signs + (0,), face))
                    continue
                # If we are here, then ``face`` is not contained in
                # ``hyperplane``.
                if zero_part_dim >= 0:
                    # Do not append ``zero_part`` yet! It might be
                    # redundant (in the sense that some of its defining
                    # inequalities are always equalities on it). Check for
                    # this:
                    zero_part_point = zero_part.representative_point()
                    for l, testhype in enumerate(hypes[:k]):
                        if signs[l] != 0:
                            h = testhype.dense_coefficient_list()
                            testval = R.sum(h[i+1] * gi for i, gi in enumerate(zero_part_point)) + h[0]
                            if testval == 0:
                                break
                    else:
                        # Now we know ``zero_part`` is not redundant.
                        subdivided.append((signs + (0,), zero_part))
                # Adding the intersection of ``face`` with the positive
                # halfspace:
                pos_part = face.intersection(pos_half)
                pos_part_dim = pos_part.dim()
                if pos_part_dim == face_dim:
                    # If this condition is not satisfied, then
                    # ``pos_part`` is either ``zero_part`` or the empty
                    # set; in either case we need not add it. Conversely,
                    # if it is satisfied, then ``pos_part`` is not yet in
                    # ``subdivided``, nor is it redundant.
                    subdivided.append((signs + (1,), pos_part))
                neg_part = face.intersection(neg_half)
                neg_part_dim = neg_part.dim()
                if neg_part_dim == face_dim:
                    # If this condition is not satisfied, then
                    # ``neg_part`` is either ``zero_part`` or the empty
                    # set; in either case we need not add it. Conversely,
                    # if it is satisfied, then ``neg_part`` is not yet in
                    # ``subdivided``, nor is it redundant.
                    subdivided.append((signs + (-1,), neg_part))
            faces = subdivided
        if labelled:
            return tuple(faces)
            # Or, if we want a dictionary:
            # return {F[0]: F[1] for F in faces}
        return tuple(x[1] for x in faces)

    def face_product(self, F, G, normalize=True):
        r"""
        Return the product `FG` in the face semigroup of ``self``, where
        `F` and `G` are two closed faces of ``self``.

        The face semigroup of a hyperplane arrangement `\mathcal{A}` is
        defined as follows: As a set, it is the set of all open faces
        of ``self`` (see :meth:`closed_faces`). Its product is defined by
        the following rule: If `F` and `G` are two open faces of
        `\mathcal{A}`, then `FG` is an open face of `\mathcal{A}`, and
        for every hyperplane `H \in \mathcal{A}`, the open face `FG` lies
        on the same side of `H` as `F` unless `F \subseteq H`, in which
        case `FG` lies on the same side of `H` as `G`. Alternatively,
        `FG` can be defined as follows: If `f` and `g` are two points in
        `F` and `G`, respectively, then `FG` is the face that contains
        the point `(f + \varepsilon g) / (1 + \varepsilon)` for any
        sufficiently small positive `\varepsilon`.

        In our implementation, the face semigroup consists of closed faces
        rather than open faces (thanks to the 1-to-1 correspondence
        between open faces and closed faces, this is not really a
        different semigroup); these closed faces are given as polyhedra.

        The face semigroup of a hyperplane arrangement is always a
        left-regular band (i.e., a semigroup satisfying the identities
        `x^2 = x` and `xyx = xy`). When the arrangement is central, then
        this semigroup is a monoid. See [Br2000]_ (Appendix A in
        particular) for further properties.

        INPUT:

        - ``F``, ``G`` -- two faces of ``self`` (as polyhedra)

        - ``normalize`` -- boolean (default: ``True``); if ``True``, then
          this method returns the precise instance of `FG` in the list
          returned by ``self.closed_faces()``, rather than creating a new
          instance

        EXAMPLES::

            sage: # needs sage.graphs
            sage: a = hyperplane_arrangements.braid(3)
            sage: a.hyperplanes()
            (Hyperplane 0*t0 + t1 - t2 + 0,
             Hyperplane t0 - t1 + 0*t2 + 0,
             Hyperplane t0 + 0*t1 - t2 + 0)
            sage: faces = {F0: F1 for F0, F1 in a.closed_faces()}
            sage: xGyEz = faces[(0, 1, 1)]   # closed face x >= y = z
            sage: xGyEz.representative_point()
            (0, -1, -1)
            sage: xGyEz = faces[(0, 1, 1)]   # closed face x >= y = z
            sage: xGyEz.representative_point()
            (0, -1, -1)
            sage: yGxGz = faces[(1, -1, 1)]  # closed face y >= x >= z
            sage: xGyGz = faces[(1, 1, 1)]   # closed face x >= y >= z
            sage: a.face_product(xGyEz, yGxGz) == xGyGz
            True
            sage: a.face_product(yGxGz, xGyEz) == yGxGz
            True
            sage: xEzGy = faces[(-1, 1, 0)]  # closed face x = z >= y
            sage: xGzGy = faces[(-1, 1, 1)]  # closed face x >= z >= y
            sage: a.face_product(xEzGy, yGxGz) == xGzGy
            True
        """
        f = F.representative_point()
        g = G.representative_point()
        n = len(f)
        R = self.base_ring()
        from sage.geometry.polyhedron.constructor import Polyhedron
        eqns = [[0] + [0] * n]
        ieqs = []
        signs = []
        for hyperplane in self.hyperplanes():
            # Decide which side of ``hyperplane`` our face ``FG`` will be
            # on.
            H = hyperplane.dense_coefficient_list()
            ieq = vector(R, H)
            x = R.sum(H[i+1] * fi for i, fi in enumerate(f)) + H[0]
            if x < 0:
                side = -1
            elif x > 0:
                side = 1
            else:
                x = R.sum(H[i+1] * gi for i, gi in enumerate(g)) + H[0]
                if x < 0:
                    side = -1
                elif x > 0:
                    side = 1
                else:
                    side = 0
            signs.append(side)
            if side == 0:
                eqns.append(ieq)
            elif side == -1:
                ieqs.append(-ieq)
            else:
                ieqs.append(ieq)
        face = Polyhedron(eqns=eqns, ieqs=ieqs, base_ring=R, backend=self._backend)
        if not normalize:
            return face
        # Look for ``I`` in ``self.closed_faces()``:
        for I in self.closed_faces():
            if I[0] == tuple(signs):
                return I[1]

    def face_semigroup_algebra(self, field=None, names='e'):
        r"""
        Return the face semigroup algebra of ``self``.

        This is the semigroup algebra of the face semigroup of ``self``
        (see :meth:`face_product` for the definition of the semigroup).

        Due to limitations of the current Sage codebase (e.g., semigroup
        algebras do not profit from the functionality of the
        :class:`FiniteDimensionalAlgebra` class), this is implemented not
        as a semigroup algebra, but as a
        :class:`FiniteDimensionalAlgebra`. The closed faces of ``self``
        (in the order in which the :meth:`closed_faces` method outputs
        them) are identified with the vectors `(0, 0, \ldots, 0, 1, 0, 0,
        \ldots, 0)` (with the `1` moving from left to right).

        INPUT:

        - ``field`` -- a field (default: `\QQ`), to be used as the
          base ring for the algebra (can also be a commutative ring, but
          then certain representation-theoretical methods might misbehave)

        - ``names`` -- (default: ``'e'``) string; names for the basis
          elements of the algebra

        .. TODO::

            Also implement it as an actual semigroup algebra?

        EXAMPLES::

            sage: # needs sage.graphs
            sage: a = hyperplane_arrangements.braid(3)
            sage: [(i, F[0]) for i, F in enumerate(a.closed_faces())]
            [(0, (0, 0, 0)),
             (1, (0, 1, 1)),
             (2, (0, -1, -1)),
             (3, (1, 0, 1)),
             (4, (1, 1, 1)),
             (5, (1, -1, 0)),
             (6, (1, -1, 1)),
             (7, (1, -1, -1)),
             (8, (-1, 0, -1)),
             (9, (-1, 1, 0)),
             (10, (-1, 1, 1)),
             (11, (-1, 1, -1)),
             (12, (-1, -1, -1))]
            sage: U = a.face_semigroup_algebra(); U
            Finite-dimensional algebra of degree 13 over Rational Field
            sage: e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12 = U.basis()
            sage: e0 * e1
            e1
            sage: e0 * e5
            e5
            sage: e5 * e0
            e5
            sage: e3 * e2
            e6
            sage: e7 * e12
            e7
            sage: e3 * e12
            e6
            sage: e4 * e8
            e4
            sage: e8 * e4
            e11
            sage: e8 * e1
            e11
            sage: e5 * e12
            e7
            sage: (e3 + 2*e4) * (e1 - e7)
            e4 - e6

            sage: U3 = a.face_semigroup_algebra(field=GF(3)); U3                        # needs sage.graphs sage.rings.finite_rings
            Finite-dimensional algebra of degree 13 over Finite Field of size 3

        TESTS:

        The ``names`` keyword works::

            sage: # needs sage.graphs
            sage: a = hyperplane_arrangements.braid(3)
            sage: U = a.face_semigroup_algebra(names='x'); U
            Finite-dimensional algebra of degree 13 over Rational Field
            sage: e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12 = U.basis()
            sage: e0 * e1
            x1
        """
        if field is None:
            from sage.rings.rational_field import QQ
            field = QQ
        zero = field.zero()
        one = field.one()
        from sage.matrix.matrix_space import MatrixSpace
        Fs = [F0 for F0, F1 in self.closed_faces()]
        # ``Fs`` is the list of the sign vectors of all closed faces of
        # ``self``.
        Fdict = {v: i for i, v in enumerate(Fs)}
        # ``Fdict`` is a dictionary whose keys are the sign vectors of the
        # closed faces of ``self``, and whose values are their positions
        # in the list ``Fs``.
        N = len(Fs)
        # Some hackery to generate a matrix quickly and without
        # unnecessary sanitization/ducktyping:
        MS = MatrixSpace(field, N, N)
        table = []
        for j, sj in enumerate(Fs):
            matrix_j = []
            for i, si in enumerate(Fs):
                row_i = [zero] * N
                sk = [sil if sil != 0 else sj[l]
                      for l, sil in enumerate(si)]
                k = Fdict[tuple(sk)]
                row_i[k] = one
                matrix_j += row_i
            table.append(MS(matrix_j, coerce=False))
        from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra import FiniteDimensionalAlgebra as FDA
        return FDA(field, table, names=names, assume_associative=True)

    def region_containing_point(self, p):
        r"""
        The region in the hyperplane arrangement containing a given point.

        The base field must have characteristic zero.

        INPUT:

        - ``p`` -- point

        OUTPUT:

        A polyhedron. A :exc:`ValueError` is raised if the point is not
        interior to a region, that is, sits on a hyperplane.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,0), 0], [(0,1), 1], [(0,1), -1], [(1,-1), 0], [(1,1), 0])
            sage: A.region_containing_point([1,2])
            A 2-dimensional polyhedron in QQ^2 defined
            as the convex hull of 2 vertices and 2 rays

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
                            lines=[hyperplane.normal() for hyperplane in self],
                            backend=self._backend)
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

            sage: # needs sage.combinat
            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.bounded_regions()
            (A 3-dimensional polyhedron in QQ^3 defined
                 as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined
                 as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined
                 as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined
                 as the convex hull of 6 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined
                 as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined
                 as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined
                 as the convex hull of 3 vertices and 1 line)
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

            sage: # needs sage.combinat
            sage: A = hyperplane_arrangements.semiorder(3)
            sage: B = A.essentialization()
            sage: B.n_regions() - B.n_bounded_regions()
            12
            sage: B.unbounded_regions()
            (A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined
                as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined
                as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined
                 as the convex hull of 1 vertex and 2 rays)
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

            sage: # needs sage.combinat
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
            m1[r(i), r(j)] += p.moebius_function(i, j)
            m2[r(i), r(j)] += 1
        m1.set_immutable()
        m2.set_immutable()
        return (m1, m2)

    def doubly_indexed_whitney_number(self, i, j, kind=1):
        r"""
        Return the `i,j`-th  doubly-indexed Whitney number.

        If ``kind=1``, this number is obtained by adding the Möbius function
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

            sage: # needs sage.combinat
            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.doubly_indexed_whitney_number(0, 2)
            9
            sage: A.whitney_number(2)
            9
            sage: A.doubly_indexed_whitney_number(1, 2)
            -15

        REFERENCES:

        - [GZ1983]_
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

        If ``kind=1``, this number is obtained by summing the Möbius function
        values `mu(0, x)` over all `x` in the intersection poset with
        `\mathrm{rank}(x) = k`.

        If ``kind=2``, this number is the number of elements `x, y` in the
        intersection poset such that `x \leq y` with ranks `i` and `j`,
        respectively.

        See [GZ1983]_ for more details.

        INPUT:

        - ``k`` -- integer

        - ``kind`` -- 1 or 2 (default: 1)

        OUTPUT:

        Integer. The ``k``-th Whitney number.

        .. SEEALSO::

            :meth:`doubly_indexed_whitney_number`
            :meth:`whitney_data`

        EXAMPLES::

            sage: # needs sage.combinat
            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.whitney_number(0)
            1
            sage: A.whitney_number(1)
            -6
            sage: A.whitney_number(2)
            9
            sage: A.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x
            sage: A.whitney_number(1, kind=2)
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

        OUTPUT: boolean; whether the hyperplane ``hyperplane`` separate the
        given regions

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

        OUTPUT: integer; the number of hyperplanes separating the two regions

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
        d = [self.distance_between_regions(r, base_region) for r in self.regions()]
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
          variable names for the polynomial ring `S`

        OUTPUT: the Varchenko matrix

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(3)
            sage: v = a.varchenko_matrix();  v
            [       1       h2       h1    h1*h2 h0*h1*h2    h0*h1    h0*h2       h0]
            [      h2        1    h1*h2       h1    h0*h1 h0*h1*h2       h0    h0*h2]
            [      h1    h1*h2        1       h2    h0*h2       h0 h0*h1*h2    h0*h1]
            [   h1*h2       h1       h2        1       h0    h0*h2    h0*h1 h0*h1*h2]
            [h0*h1*h2    h0*h1    h0*h2       h0        1       h2       h1    h1*h2]
            [   h0*h1 h0*h1*h2       h0    h0*h2       h2        1    h1*h2       h1]
            [   h0*h2       h0 h0*h1*h2    h0*h1       h1    h1*h2        1       h2]
            [      h0    h0*h2    h0*h1 h0*h1*h2    h1*h2       h1       h2        1]
            sage: factor(det(v))
            (h2 - 1)^4 * (h2 + 1)^4 * (h1 - 1)^4 * (h1 + 1)^4 * (h0 - 1)^4 * (h0 + 1)^4

        TESTS:

        Verify that :issue:`36490` is fixed::

            sage: hyperplane_arrangements.coordinate(1).varchenko_matrix()
            [1 h]
            [h 1]
        """
        from sage.matrix.constructor import identity_matrix
        from sage.misc.misc_c import prod
        k = len(self)
        R = PolynomialRing(QQ, names, k)
        h = R.gens()
        region = self.regions()
        n = len(region)
        v = identity_matrix(R, n, n)
        for i in range(n):
            for j in range(i + 1, n):
                t = prod(h[p] for p in range(k) if
                         self.is_separating_hyperplane(region[i], region[j], self[p]))
                v[i, j] = v[j, i] = t
        v.set_immutable()
        return v

    @cached_method
    def matroid(self):
        r"""
        Return the matroid associated to ``self``.

        Let `A` denote a central hyperplane arrangement and `n_H` the
        normal vector of some hyperplane `H \in A`. We define a matroid
        `M_A` as the linear matroid spanned by `\{ n_H | H \in A \}`.
        The matroid `M_A` is such that the lattice of flats of `M` is
        isomorphic to the intersection lattice of `A`
        (Proposition 3.6 in [Sta2007]_).

        EXAMPLES::

            sage: P.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = P(x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z)
            sage: M = A.matroid(); M
            Linear matroid of rank 3 on 7 elements represented over the Rational Field

        We check the lattice of flats is isomorphic to the
        intersection lattice::

            sage: f = sum([list(M.flats(i)) for i in range(M.rank() + 1)], [])
            sage: PF = Poset([f, lambda x, y: x < y])                                   # needs sage.combinat
            sage: PF.is_isomorphic(A.intersection_poset())                              # needs sage.combinat
            True
        """
        if not self.is_central():
            raise ValueError("the hyperplane arrangement must be central")
        norms = [p.normal() for p in self]
        from sage.matroids.constructor import Matroid
        return Matroid(matrix=matrix(norms).transpose())

    def orlik_solomon_algebra(self, base_ring=None, ordering=None, **kwds):
        """
        Return the Orlik-Solomon algebra of ``self``.

        INPUT:

        - ``base_ring`` -- (default: the base field of ``self``) the ring
          over which the Orlik-Solomon algebra will be defined
        - ``ordering`` -- (optional) an ordering of the ground set

        EXAMPLES::

            sage: P.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = P(x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z)
            sage: A.orlik_solomon_algebra()
            Orlik-Solomon algebra of Linear matroid of rank 3 on 7 elements
             represented over the Rational Field
            sage: A.orlik_solomon_algebra(base_ring=ZZ)
            Orlik-Solomon algebra of Linear matroid of rank 3 on 7 elements
             represented over the Rational Field
        """
        if base_ring is None:
            base_ring = self.base_ring()
        return self.matroid().orlik_solomon_algebra(base_ring, ordering, **kwds)

    def orlik_terao_algebra(self, base_ring=None, ordering=None, **kwds):
        """
        Return the Orlik-Terao algebra of ``self``.

        INPUT:

        - ``base_ring`` -- (default: the base field of ``self``) the ring
          over which the Orlik-Terao algebra will be defined
        - ``ordering`` -- (optional) an ordering of the ground set

        EXAMPLES::

            sage: P.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = P(x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z)
            sage: A.orlik_terao_algebra()
            Orlik-Terao algebra of Linear matroid of rank 3 on 7 elements
             represented over the Rational Field over Rational Field
            sage: A.orlik_terao_algebra(base_ring=QQ['t'])
            Orlik-Terao algebra of Linear matroid of rank 3 on 7 elements
             represented over the Rational Field
             over Univariate Polynomial Ring in t over Rational Field
        """
        if base_ring is None:
            base_ring = self.base_ring()
        return self.matroid().orlik_terao_algebra(base_ring, ordering, **kwds)

    @cached_method
    def minimal_generated_number(self):
        r"""
        Return the minimum `k` such that ``self`` is `k`-generated.

        Let `A` be a central hyperplane arrangement. Let `W_k` denote
        the solution space of the linear system corresponding to the
        linear dependencies among the hyperplanes of `A` of length at
        most `k`. We say `A` is `k`-*generated* if
        `\dim W_k = \operatorname{rank} A`.

        Equivalently this says all dependencies forming the Orlik-Terao
        ideal are generated by at most `k` hyperplanes.

        EXAMPLES:

        We construct Example 2.2 from [Yuz1993]_::

            sage: P.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = P(x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z, 3*x+5*z, 3*x+4*y+5*z)
            sage: B = P(x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z, x+3*z, x+2*y+3*z)
            sage: A.minimal_generated_number()
            3
            sage: B.minimal_generated_number()
            4

        TESTS:

        Check that :issue:`26705` is fixed::

            sage: # needs sage.combinat sage.groups
            sage: w = WeylGroup(['A', 4]).from_reduced_word([3, 4, 2, 1])
            sage: I = w.inversion_arrangement()
            sage: I
            Arrangement <a4 | a1 | a1 + a2 | a1 + a2 + a3 + a4>
            sage: I.minimal_generated_number()
            0
            sage: I.is_formal()
            True
        """
        V = VectorSpace(self.base_ring(), self.dimension())
        W = VectorSpace(self.base_ring(), self.n_hyperplanes())
        r = self.rank()
        M = self.matroid()
        if len(M.groundset()) == r:  # there are no circuits
            return ZZ.zero()
        norms = M.representation().columns()
        circuits = M.circuits()
        for i in range(2, self.n_hyperplanes()):
            sol = []
            for d in circuits:
                if len(d) > i:
                    continue
                d = list(d)
                dep = V.linear_dependence([norms[j] for j in d])
                w = W.zero().list()
                for j, k in enumerate(d):
                    w[k] = dep[0][j]
                sol.append(w)
            mat = matrix(sol)
            if mat.right_kernel().dimension() == r:
                return i
        return self.n_hyperplanes()

    def is_formal(self):
        """
        Return if ``self`` is formal.

        A hyperplane arrangement is *formal* if it is 3-generated [Yuz1993]_,
        where `k`-generated is defined in :meth:`minimal_generated_number`.

        EXAMPLES::

            sage: P.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = P(x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z, 3*x+5*z, 3*x+4*y+5*z)
            sage: B = P(x, y, z, x+y+z, 2*x+y+z, 2*x+3*y+z, 2*x+3*y+4*z, x+3*z, x+2*y+3*z)
            sage: A.is_formal()
            True
            sage: B.is_formal()
            False
        """
        return self.minimal_generated_number() <= 3

    def defining_polynomial(self):
        r"""
        Return the defining polynomial of ``A``.

        Let `A = (H_i)_i` be a hyperplane arrangement in a vector space `V`
        corresponding to the null spaces of `\alpha_{H_i} \in V^*`. Then
        the *defining polynomial* of `A` is given by

        .. MATH::

            Q(A) = \prod_i \alpha_{H_i} \in S(V^*).

        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H([2*x + y - z, -x - 2*y + z])
            sage: p = A.defining_polynomial(); p
            -2*x^2 - 5*x*y - 2*y^2 + 3*x*z + 3*y*z - z^2
            sage: p.factor()
            (-1) * (x + 2*y - z) * (2*x + y - z)
        """
        S = self.parent().ambient_space().symmetric_space()
        return S.prod(H.to_symmetric_space() for H in self)

    @cached_method
    def derivation_module_free_chain(self):
        r"""
        Return a free chain for the derivation module if one
        exists, otherwise return ``None``.

        .. SEEALSO::

            :meth:`is_free`

        EXAMPLES::

            sage: # needs sage.combinat sage.groups
            sage: W = WeylGroup(['A',3], prefix='s')
            sage: A = W.long_element().inversion_arrangement()
            sage: for M in A.derivation_module_free_chain(): print("%s\n"%M)
            [ 1  0  0]
            [ 0  1  0]
            [ 0  0 a3]
            <BLANKLINE>
            [ 1  0  0]
            [ 0  0  1]
            [ 0 a2  0]
            <BLANKLINE>
            [  1   0   0]
            [  0  -1  -1]
            [  0  a2 -a3]
            <BLANKLINE>
            [ 0  1  0]
            [ 0  0  1]
            [a1  0  0]
            <BLANKLINE>
            [ 1  0 -1]
            [a3 -1  0]
            [a1  0 a2]
            <BLANKLINE>
            [       1        0        0]
            [      a3       -1       -1]
            [       0       a1 -a2 - a3]
            <BLANKLINE>
        """
        if not self.is_central():
            raise NotImplementedError("only implemented for central arrangements")
        from sage.geometry.hyperplane_arrangement.check_freeness import construct_free_chain
        return construct_free_chain(self)

    @cached_method(key=lambda self, a: None)
    def is_free(self, algorithm='singular'):
        r"""
        Return if ``self`` is free.

        A hyperplane arrangement `A` is free if the module
        of derivations `\operatorname{Der}(A)` is a free `S`-module,
        where `S` is the corresponding symmetric space.

        INPUT:

        - ``algorithm`` -- (default: ``'singular'``) can be one of
          the following:

          * ``'singular'`` -- use Singular's minimal free resolution
          * ``'BC'`` -- use the algorithm given by Barakat and Cuntz
            in [BC2012]_ (much slower than using Singular)

        ALGORITHM:

        .. RUBRIC:: singular

        Check that the minimal free resolution has length at most 2
        by using Singular.

        .. RUBRIC:: BC

        This implementation follows [BC2012]_ by constructing a chain
        of free modules

        .. MATH::

            D(A) = D(A_n) < D(A_{n-1}) < \cdots < D(A_1) < D(A_0)

        corresponding to some ordering of the arrangements `A_0 \subset
        A_1 \subset \cdots \subset A_{n-1} \subset A_n = A`. Such a
        chain is found by using a backtracking algorithm.

        EXAMPLES:

        For type `A` arrangements, chordality is equivalent to freeness.
        We verify that in type `A_3`::

            sage: W = WeylGroup(['A', 3], prefix='s')                                   # needs sage.combinat sage.groups
            sage: for x in W:                                                           # needs sage.combinat sage.groups
            ....:    A = x.inversion_arrangement()
            ....:    assert A.matroid().is_chordal() == A.is_free()

        TESTS:

        We check that the algorithms agree::

            sage: W = WeylGroup(['B', 3], prefix='s')                                   # needs sage.combinat sage.groups
            sage: for x in W:                   # long time                             # needs sage.combinat sage.groups
            ....:    A = x.inversion_arrangement()
            ....:    assert (A.is_free(algorithm='BC')
            ....:            == A.is_free(algorithm='singular'))
        """
        if not self.is_central():
            raise NotImplementedError("only implemented for central arrangements")
        if algorithm == "singular":
            # TODO: Implement this using libSingular
            mres = self.defining_polynomial().jacobian_ideal()._singular_().mres(0)
            return len(mres) <= 2
        elif algorithm == "BC":
            return self.derivation_module_free_chain() is not None
        else:
            raise ValueError("invalid algorithm")

    def derivation_module_basis(self, algorithm='singular'):
        """
        Return a basis for the derivation module of ``self`` if
        one exists, otherwise return ``None``.

        .. SEEALSO::

            :meth:`derivation_module_free_chain`, :meth:`is_free`

        INPUT:

        - ``algorithm`` -- (default: ``'singular'``) can be one of
          the following:

          * ``'singular'`` -- use Singular's minimal free resolution
          * ``'BC'`` -- use the algorithm given by Barakat and Cuntz
            in [BC2012]_ (much slower than using Singular)

        OUTPUT:

        A basis for the derivation module (over `S`, the
        :meth:`symmetric space
        <sage.geometry.hyperplane_arrangement.hyperplane.AmbientVectorSpace.symmetric_space>`)
        as vectors of a free module over `S`.

        ALGORITHM:

        .. RUBRIC:: Singular

        This gets the reduced syzygy module of the Jacobian ideal of
        the defining polynomial `f` of ``self``. It then checks Saito's
        criterion that the determinant of the basis matrix is a scalar
        multiple of `f`. If the basis matrix is not square or it fails
        Saito's criterion, then we check if the arrangement is free.
        If it is free, then we fall back to the Barakat-Cuntz algorithm.

        .. RUBRIC:: BC

        Return the product of the derivation module free chain matrices.
        See Section 6 of [BC2012]_.

        EXAMPLES::

            sage: # needs sage.combinat sage.groups
            sage: W = WeylGroup(['A', 2], prefix='s')
            sage: A = W.long_element().inversion_arrangement()
            sage: A.derivation_module_basis()
            [(a1, a2), (0, a1*a2 + a2^2)]

        TESTS:

        We check the algorithms produce a basis with the same exponents::

            sage: W = WeylGroup(['A', 2], prefix='s')                                   # needs sage.combinat sage.groups
            sage: def exponents(B):
            ....:     return sorted([max(x.degree() for x in b) for b in B])
            sage: for x in W:                   # long time                             # needs sage.combinat sage.groups
            ....:     A = x.inversion_arrangement()
            ....:     B = A.derivation_module_basis(algorithm='singular')
            ....:     Bp = A.derivation_module_basis(algorithm='BC')
            ....:     if B is None:
            ....:         assert Bp is None
            ....:     else:
            ....:         assert exponents(B) == exponents(Bp)
        """
        alg = algorithm  # prevent possible changes to a global variable
        if alg == "singular":
            # import sage.libs.singular.function_factory
            # syz = sage.libs.singular.function_factory.ff.syz
            f = self.defining_polynomial()
            I = f + f.jacobian_ideal()
            IS = I._singular_()
            ISS = IS.syz()
            MSTD = ISS.mstd()
            basis = MSTD[2]._sage_().transpose().submatrix(0, 1)
            try:
                det = basis.det()
                # Check using Saito's criterion
                if det / f in f.parent().base_ring() and not det.is_zero():
                    return basis.rows()
            except ValueError:  # Non-square matrix or det = 0
                pass
            # Check if it is free
            if not self.is_free(algorithm=alg):
                return None
            # The syzygy module did not give a basis, but since it is free,
            #    fallback to the Barakat-Cuntz method
            alg = "BC"
        if alg == "BC":
            C = self.derivation_module_free_chain()
            if C is not None:
                if not C:  # C is an empty list
                    S = self.parent().ambient_space().symmetric_space()
                    return matrix.identity(S, self.dimension()).rows()
                from sage.misc.misc_c import prod
                return prod(reversed(C)).rows()
            return None
        else:
            raise ValueError("invalid algorithm")


class HyperplaneArrangements(Parent, UniqueRepresentation):
    """
    Hyperplane arrangements.

    For more information on hyperplane arrangements, see
    :mod:`sage.geometry.hyperplane_arrangement.arrangement`.

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
        from sage.categories.sets_cat import Sets
        from sage.rings.ring import _Fields
        if base_ring not in _Fields:
            raise ValueError('base ring must be a field')
        super().__init__(category=Sets())
        self._base_ring = base_ring
        self._names = names

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT: the base ring of the hyperplane arrangement

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

        - ``base_ring`` -- a ring; the new base ring

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

        OUTPUT: string

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

        - ``signed`` -- boolean (default: ``True``); whether to
          preserve signs of hyperplane equations

        - ``warn_duplicates`` -- boolean (default: ``False``);
          whether to issue a warning if duplicate hyperplanes were
          passed -- note that duplicate hyperplanes are always removed,
          whether or not there is a warning shown

        - ``check`` -- boolean (default: ``True``); whether to
          perform argument checking

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

            sage: L._element_constructor_(polytopes.hypercube(2))
            Arrangement <-x + 1 | -y + 1 | y + 1 | x + 1>

            sage: L(x, x, warn_duplicates=True)
            doctest:...: UserWarning: Input contained 2 hyperplanes, but only 1 are distinct.
            Arrangement <x>
            sage: L(-x, x + y - 1, signed=False)
            Arrangement <-x - y + 1 | x>

        TESTS::

            sage: L()
            Empty hyperplane arrangement of dimension 2
            sage: L(0)        # zero is equivalent to no argument, Issue #8648
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
        backend = kwds.pop('backend', None)
        if len(kwds) > 0:
            raise ValueError('unknown keyword argument')
        # process positional arguments
        AA = self.ambient_space()
        try:
            hyperplanes = [AA(_) for _ in args]
        except (TypeError, ValueError, AttributeError):
            if len(args) > 1:
                raise
            arg = args[0]
            if hasattr(arg, 'Hrepresentation'):
                hyperplanes = [AA(h) for h in arg.Hrepresentation()]
            else:
                hyperplanes = [AA(_) for _ in arg]
        hyperplanes = [h.primitive(signed) for h in hyperplanes]
        n = len(hyperplanes)
        hyperplanes = set(hyperplanes)
        if warn_duplicates and n != len(hyperplanes):
            from warnings import warn
            warn('Input contained {0} hyperplanes, but only {1} are distinct.'.format(n, len(hyperplanes)))
        # argument checking (optional but recommended)
        if check:
            if signed and not not_char2:
                raise ValueError('cannot be signed in characteristic 2')
            for h in hyperplanes:
                if h.A() == 0:
                    raise ValueError('linear expression must be non-constant to define a hyperplane')
                if not_char2 and -h in hyperplanes:
                    raise ValueError('arrangement cannot simultaneously have h and -h as hyperplane')
        return self.element_class(self, tuple(sorted(hyperplanes)), backend=backend)

    @cached_method
    def ngens(self):
        """
        Return the number of linear variables.

        OUTPUT: integer

        EXAMPLES::

            sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 3-dimensional linear space
             over Rational Field with coordinates x, y, z
            sage: L.ngens()
            3
        """
        return len(self._names)

    @cached_method
    def gens(self) -> tuple:
        """
        Return the coordinate hyperplanes.

        OUTPUT: a tuple of linear expressions, one for each linear variable

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

        OUTPUT: a linear expression

        EXAMPLES::

            sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in
             3-dimensional linear space over Rational Field with coordinates x, y, z
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
            Coercion map:
              From: Integer Ring
              To:   Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
            sage: M.coerce_map_from(L)
            Coercion map:
              From: Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
              To:   Hyperplane arrangements in 1-dimensional linear space over Real Field with 53 bits of precision with coordinate y
            sage: L.coerce_map_from(M)
        """
        if self.ambient_space().has_coerce_map_from(P):
            return True
        if isinstance(P, HyperplaneArrangements):
            return self.base_ring().has_coerce_map_from(P.base_ring())
        return super()._coerce_map_from_(P)
