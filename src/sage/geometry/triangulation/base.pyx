# distutils: sources = sage/geometry/triangulation/functions.cc sage/geometry/triangulation/data.cc sage/geometry/triangulation/triangulations.cc
# distutils: depends = sage/geometry/triangulation/functions.h sage/geometry/triangulation/data.h sage/geometry/triangulation/triangulations.h
# distutils: language = c++

r"""
Base classes for triangulations

We provide (fast) cython implementations here.

AUTHORS:

- Volker Braun (2010-09-14): initial version.
"""


########################################################################
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

from sage.misc.fast_methods cimport hash_by_id
from sage.structure.sage_object cimport SageObject
from sage.structure.parent cimport Parent
from sage.categories.sets_cat import Sets
from sage.matrix.constructor import matrix

from sage.geometry.triangulation.functions cimport binomial
from sage.geometry.triangulation.triangulations cimport \
    triangulations_ptr, init_triangulations, next_triangulation, delete_triangulations


########################################################################
cdef class Point(SageObject):
    r"""
    A point of a point configuration.

    Note that the coordinates of the points of a point configuration
    are somewhat arbitrary. What counts are the abstract linear
    relations between the points, for example encoded by the
    :meth:`~sage.geometry.triangulation.point_configuration.PointConfiguration.circuits`.

    .. WARNING::

        You should not create :class:`Point` objects manually. The
        constructor of :class:`PointConfiguration_base` takes care of
        this for you.

    INPUT:

    - ``point_configuration`` -- :class:`PointConfiguration_base`; the
      point configuration to which the point belongs

    - ``i`` -- integer; the index of the point in the point
      configuration

    - ``projective`` -- the projective coordinates of the point

    - ``affine`` -- the affine coordinates of the point

    - ``reduced`` -- the reduced (with linearities removed)
      coordinates of the point

    EXAMPLES::

        sage: pc = PointConfiguration([(0,0)])
        sage: from sage.geometry.triangulation.base import Point
        sage: Point(pc, 123, (0,0,1), (0,0), ())
        P(0, 0)
    """

    cdef int _index
    cdef tuple _projective, _affine, _reduced_affine
    cdef object _point_configuration
    cdef object _reduced_affine_vector, _reduced_projective_vector

    def __init__(self, point_configuration, i, projective, affine, reduced):
        r"""
        Construct a :class:`Point`.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0)])
            sage: from sage.geometry.triangulation.base import Point
            sage: Point(pc, 123, (0,0,1), (0,0), ())   # indirect doctest
            P(0, 0)
        """
        self._index = i
        self._projective = tuple(projective)
        self._affine = tuple(affine)
        self._reduced_affine = tuple(reduced)
        self._point_configuration = point_configuration
        V = point_configuration.reduced_affine_vector_space()
        self._reduced_affine_vector = V(self._reduced_affine)
        P = point_configuration.reduced_projective_vector_space()
        self._reduced_projective_vector = P(self.reduced_projective())

    def __hash__(self):
        r"""
        Hash value for a point in a point configuration.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,1]])
            sage: hash(p[0]) # random
            35822008390213632
        """
        return hash(self._point_configuration) ^ (<long>self._index)

    cpdef point_configuration(self):
        r"""
        Return the point configuration to which the point belongs.

        OUTPUT: a :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`

        EXAMPLES::

            sage: pc = PointConfiguration([ (0,0), (1,0), (0,1) ])
            sage: p = pc.point(0)
            sage: p is pc.point(0)
            True
            sage: p.point_configuration() is pc
            True
        """
        return self._point_configuration

    def __iter__(self):
        r"""
        Iterate through the affine ambient space coordinates of the point.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0)])
            sage: from sage.geometry.triangulation.base import Point
            sage: p = Point(pc, 123, (1,2,1), (3,4), ())
            sage: list(p)  # indirect doctest
            [3, 4]
        """
        return iter(self._affine)

    def __len__(self):
        r"""
        Return the affine ambient space dimension.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0)])
            sage: from sage.geometry.triangulation.base import Point
            sage: p = Point(pc, 123, (1,2,1), (3,4), ())
            sage: len(p)
            2
            sage: p.__len__()
            2
        """
        return len(self._affine)

    cpdef index(self):
        """
        Return the index of the point in the point configuration.

        EXAMPLES::

            sage: pc = PointConfiguration([[0, 1], [0, 0], [1, 0]])
            sage: p = pc.point(2); p
            P(1, 0)
            sage: p.index()
            2
        """
        return self._index

    cpdef projective(self):
        r"""
        Return the projective coordinates of the point in the ambient space.

        OUTPUT: a tuple containing the coordinates

        EXAMPLES::

            sage: pc = PointConfiguration([[10, 0, 1], [10, 0, 0], [10, 2, 3]])
            sage: p = pc.point(2); p
            P(10, 2, 3)
            sage: p.affine()
            (10, 2, 3)
            sage: p.projective()
            (10, 2, 3, 1)
            sage: p.reduced_affine()
            (2, 2)
            sage: p.reduced_projective()
            (2, 2, 1)
            sage: p.reduced_affine_vector()
            (2, 2)
        """
        return self._projective

    cpdef affine(self):
        r"""
        Return the affine coordinates of the point in the ambient space.

        OUTPUT: a tuple containing the coordinates

        EXAMPLES::

            sage: pc = PointConfiguration([[10, 0, 1], [10, 0, 0], [10, 2, 3]])
            sage: p = pc.point(2); p
            P(10, 2, 3)
            sage: p.affine()
            (10, 2, 3)
            sage: p.projective()
            (10, 2, 3, 1)
            sage: p.reduced_affine()
            (2, 2)
            sage: p.reduced_projective()
            (2, 2, 1)
            sage: p.reduced_affine_vector()
            (2, 2)
        """
        return self._affine

    cpdef reduced_affine(self):
        r"""
        Return the affine coordinates of the point on the hyperplane
        spanned by the point configuration.

        OUTPUT: a tuple containing the coordinates

        EXAMPLES::

            sage: pc = PointConfiguration([[10, 0, 1], [10, 0, 0], [10, 2, 3]])
            sage: p = pc.point(2); p
            P(10, 2, 3)
            sage: p.affine()
            (10, 2, 3)
            sage: p.projective()
            (10, 2, 3, 1)
            sage: p.reduced_affine()
            (2, 2)
            sage: p.reduced_projective()
            (2, 2, 1)
            sage: p.reduced_affine_vector()
            (2, 2)
        """
        return self._reduced_affine

    cpdef reduced_projective(self):
        r"""
        Return the projective coordinates of the point on the hyperplane
        spanned by the point configuration.

        OUTPUT: a tuple containing the coordinates

        EXAMPLES::

            sage: pc = PointConfiguration([[10, 0, 1], [10, 0, 0], [10, 2, 3]])
            sage: p = pc.point(2); p
            P(10, 2, 3)
            sage: p.affine()
            (10, 2, 3)
            sage: p.projective()
            (10, 2, 3, 1)
            sage: p.reduced_affine()
            (2, 2)
            sage: p.reduced_projective()
            (2, 2, 1)
            sage: p.reduced_affine_vector()
            (2, 2)
        """
        return tuple(self._reduced_affine)+(1,)

    cpdef reduced_affine_vector(self):
        """
        Return the affine coordinates of the point on the hyperplane
        spanned by the point configuration.

        OUTPUT: a tuple containing the coordinates

        EXAMPLES::

            sage: pc = PointConfiguration([[10, 0, 1], [10, 0, 0], [10, 2, 3]])
            sage: p = pc.point(2); p
            P(10, 2, 3)
            sage: p.affine()
            (10, 2, 3)
            sage: p.projective()
            (10, 2, 3, 1)
            sage: p.reduced_affine()
            (2, 2)
            sage: p.reduced_projective()
            (2, 2, 1)
            sage: p.reduced_affine_vector()
            (2, 2)
        """
        return self._reduced_affine_vector

    cpdef reduced_projective_vector(self):
        """
        Return the affine coordinates of the point on the hyperplane
        spanned by the point configuration.

        OUTPUT: a tuple containing the coordinates

        EXAMPLES::

            sage: pc = PointConfiguration([[10, 0, 1], [10, 0, 0], [10, 2, 3]])
            sage: p = pc.point(2); p
            P(10, 2, 3)
            sage: p.affine()
            (10, 2, 3)
            sage: p.projective()
            (10, 2, 3, 1)
            sage: p.reduced_affine()
            (2, 2)
            sage: p.reduced_projective()
            (2, 2, 1)
            sage: p.reduced_affine_vector()
            (2, 2)
            sage: type(p.reduced_affine_vector())
            <class 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        """
        return self._reduced_projective_vector

    cpdef _repr_(self):
        """
        Return a string representation of the point.

        OUTPUT: string

        EXAMPLES::

            sage: pc = PointConfiguration([[10, 0, 1], [10, 0, 0]])
            sage: from sage.geometry.triangulation.base import Point
            sage: p = Point(pc, 123, (0,0,1), (0,0), (0,))
            sage: p._repr_()
            'P(0, 0)'
        """
        return 'P' + str(self._affine)


########################################################################
cdef class PointConfiguration_base(Parent):
    r"""
    The cython abstract base class for
    :class:`~sage.geometry.triangulation.PointConfiguration`.

    .. WARNING::

        You should not instantiate this base class, but only its
        derived class
        :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`.
    """

    def __init__(self, points, defined_affine):
        r"""
        Construct a :class:`PointConfiguration_base`.

        INPUT:

        - ``points`` -- tuple of tuples of projective coordinates
          with ``1`` as the final coordinate

        - ``defined_affine`` -- boolean; whether the point
          configuration is defined as a configuration of affine (as
          opposed to projective) points

        TESTS::

            sage: from sage.geometry.triangulation.base import PointConfiguration_base
            sage: PointConfiguration_base(((1,2,1),(2,3,1),(3,4,1)), False)
            <sage.geometry.triangulation.base.PointConfiguration_base object at ...>
        """
        Parent.__init__(self, category = Sets())
        self._init_points(points)
        self._is_affine = defined_affine

    cdef tuple _pts
    cdef int _ambient_dim
    cdef int _dim
    cdef object _base_ring
    cdef bint _is_affine
    cdef object _reduced_affine_vector_space, _reduced_projective_vector_space

    cdef _init_points(self, tuple projective_points):
        """
        Internal method to determine coordinates of points.

        EXAMPLES::

            sage: p = PointConfiguration([[0, 1], [0, 0], [1, 0]])   # indirect doctest
            sage: p.points()
            (P(0, 1), P(0, 0), P(1, 0))

        Special cases::

            sage: PointConfiguration([])
            The pointless empty configuration
            sage: PointConfiguration([(1,2,3)])
            A point configuration in affine 3-space over Integer Ring
            consisting of 1 point. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
        """
        if not projective_points:
            self._ambient_dim = 0
            self._dim = -1
            self._pts = tuple()
            return

        # We now are sure that projective_points is not empty
        n = len(projective_points)
        self._ambient_dim = len(projective_points[0]) - 1
        assert all(len(p) == self._ambient_dim+1 for p in projective_points), \
            'The given point coordinates must all have the same length.'
        assert len(set(projective_points)) == len(projective_points), \
            'Not all points are pairwise distinct.'

        proj = matrix(projective_points).transpose()
        self._base_ring = proj.base_ring()

        if all(x == 1 for x in proj.row(self.ambient_dim())):
            aff = proj.submatrix(0, 0, nrows=self.ambient_dim())
        else:
            raise NotImplementedError  # TODO

        if n > 1:
            # shift first point to origin
            red = matrix([aff.column(i)-aff.column(0) for i in range(n)]).transpose()
            # pick linearly independent rows
            red = matrix([red.row(i) for i in red.pivot_rows()])
        else:
            red = matrix(0, 1)
        self._dim = red.nrows()

        from sage.modules.free_module import VectorSpace
        self._reduced_affine_vector_space = VectorSpace(self._base_ring.fraction_field(), self._dim)
        self._reduced_projective_vector_space = VectorSpace(self._base_ring.fraction_field(), self._dim+1)
        self._pts = tuple([Point(self, i, proj.column(i),
                                 aff.column(i), red.column(i))
                           for i in range(n)])

    def __hash__(self):
        r"""
        Hash function.

        TESTS::

            sage: p = PointConfiguration([[0,0],[0,1]])
            sage: hash(p) # random
            8746748042501
        """
        return hash_by_id(<void *> self)

    cpdef reduced_affine_vector_space(self):
        """
        Return the vector space that contains the affine points.

        OUTPUT: a vector space over the fraction field of :meth:`base_ring`

        EXAMPLES::

            sage: p = PointConfiguration([[0,0,0], [1,2,3]])
            sage: p.base_ring()
            Integer Ring
            sage: p.reduced_affine_vector_space()
            Vector space of dimension 1 over Rational Field
            sage: p.reduced_projective_vector_space()
            Vector space of dimension 2 over Rational Field
        """
        return self._reduced_affine_vector_space

    cpdef reduced_projective_vector_space(self):
        """
        Return the vector space that is spanned by the homogeneous
        coordinates.

        OUTPUT: a vector space over the fraction field of :meth:`base_ring`

        EXAMPLES::

            sage: p = PointConfiguration([[0,0,0], [1,2,3]])
            sage: p.base_ring()
            Integer Ring
            sage: p.reduced_affine_vector_space()
            Vector space of dimension 1 over Rational Field
            sage: p.reduced_projective_vector_space()
            Vector space of dimension 2 over Rational Field
        """
        return self._reduced_projective_vector_space

    cpdef ambient_dim(self):
        """
        Return the dimension of the ambient space of the point
        configuration.

        See also :meth:`dimension`

        EXAMPLES::

            sage: p = PointConfiguration([[0,0,0]])
            sage: p.ambient_dim()
            3
            sage: p.dim()
            0
        """
        return self._ambient_dim

    cpdef dim(self):
        """
        Return the actual dimension of the point configuration.

        See also :meth:`ambient_dim`

        EXAMPLES::

            sage: p = PointConfiguration([[0,0,0]])
            sage: p.ambient_dim()
            3
            sage: p.dim()
            0
        """
        return self._dim

    cpdef base_ring(self):
        r"""
        Return the base ring, that is, the ring containing the
        coordinates of the points.

        OUTPUT: a ring

        EXAMPLES::

            sage: p = PointConfiguration([(0,0)])
            sage: p.base_ring()
            Integer Ring

            sage: p = PointConfiguration([(1/2,3)])
            sage: p.base_ring()
            Rational Field

            sage: p = PointConfiguration([(0.2, 5)])
            sage: p.base_ring()
            Real Field with 53 bits of precision
        """
        return self._base_ring

    cpdef bint is_affine(self) noexcept:
        """
        Return whether the configuration is defined by affine points.

        OUTPUT: boolean; if true, the homogeneous coordinates all have `1` as
        their last entry

        EXAMPLES::

            sage: p = PointConfiguration([(0.2, 5), (3, 0.1)])
            sage: p.is_affine()
            True

            sage: p = PointConfiguration([(0.2, 5, 1), (3, 0.1, 1)], projective=True)
            sage: p.is_affine()
            False
        """
        return self._is_affine

    def _assert_is_affine(self):
        """
        Raise a :exc:`ValueError` if the point configuration is not
        defined by affine points.

        EXAMPLES::

            sage: p = PointConfiguration([(0.2, 5), (3, 0.1)])
            sage: p._assert_is_affine()
            sage: p = PointConfiguration([(0.2, 5, 1), (3, 0.1, 1)], projective=True)
            sage: p._assert_is_affine()
            Traceback (most recent call last):
            ...
            ValueError: The point configuration contains projective points.
        """
        if not self.is_affine():
            raise ValueError('The point configuration contains projective points.')

    def __getitem__(self, i):
        """
        Return the ``i``-th point.

        Same as :meth:`point`.

        INPUT:

        - ``i`` -- integer

        OUTPUT: the ``i``-th point of the point configuration

        EXAMPLES::

            sage: p = PointConfiguration([[1,0], [2,3], [3,2]])
            sage: [ p[i] for i in range(0,p.n_points()) ]
            [P(1, 0), P(2, 3), P(3, 2)]
            sage: list(p)
            [P(1, 0), P(2, 3), P(3, 2)]
            sage: list(p.points())
            [P(1, 0), P(2, 3), P(3, 2)]
            sage: [ p.point(i) for i in range(0,p.n_points()) ]
            [P(1, 0), P(2, 3), P(3, 2)]
        """
        return self._pts[i]

    cpdef n_points(self):
        """
        Return the number of points.

        Same as ``len(self)``.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p
            A point configuration in affine 2-space over Integer Ring
            consisting of 5 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
            sage: len(p)
            5
            sage: p.n_points()
            5
        """
        return len(self._pts)

    cpdef points(self):
        """
        Return a list of the points.

        OUTPUT:

        A list of the points. See also the :meth:`__iter__`
        method, which returns the corresponding generator.

        EXAMPLES::

            sage: pconfig = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: list(pconfig)
            [P(0, 0), P(0, 1), P(1, 0), P(1, 1), P(-1, -1)]
            sage: [ p for p in pconfig.points() ]
            [P(0, 0), P(0, 1), P(1, 0), P(1, 1), P(-1, -1)]
            sage: pconfig.point(0)
            P(0, 0)
            sage: pconfig.point(1)
            P(0, 1)
            sage: pconfig.point( pconfig.n_points()-1 )
            P(-1, -1)
        """
        return self._pts

    def point(self, i):
        """
        Return the `i`-th point of the configuration.

        Same as :meth:`__getitem__`

        INPUT:

        - ``i`` -- integer

        OUTPUT: a point of the point configuration

        EXAMPLES::

            sage: pconfig = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: list(pconfig)
            [P(0, 0), P(0, 1), P(1, 0), P(1, 1), P(-1, -1)]
            sage: [ p for p in pconfig.points() ]
            [P(0, 0), P(0, 1), P(1, 0), P(1, 1), P(-1, -1)]
            sage: pconfig.point(0)
            P(0, 0)
            sage: pconfig[0]
            P(0, 0)
            sage: pconfig.point(1)
            P(0, 1)
            sage: pconfig.point( pconfig.n_points()-1 )
            P(-1, -1)
        """
        return self._pts[i]

    def __len__(self):
        """
        Return the number of points.

        Same as :meth:`n_points`

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p
            A point configuration in affine 2-space over Integer Ring
            consisting of 5 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
            sage: len(p)
            5
            sage: p.n_points()
            5
        """
        return len(self._pts)

    cpdef simplex_to_int(self, simplex):
        r"""
        Return an integer that uniquely identifies the given simplex.

        See also the inverse method :meth:`int_to_simplex`.

        The enumeration is compatible with [PUNTOS]_.

        INPUT:

        - ``simplex`` -- iterable, for example a list. The elements
          are the vertex indices of the simplex

        OUTPUT: integer that uniquely specifies the simplex

        EXAMPLES::

            sage: U=matrix([
            ....:    [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ....:    [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ....:    [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ....:    [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ....:    [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ....: ])
            sage: pc = PointConfiguration(U.columns())
            sage: pc.simplex_to_int([1,3,4,7,10,13])
            1678
            sage: pc.int_to_simplex(1678)
            (1, 3, 4, 7, 10, 13)
        """
        cdef int s = 1
        cdef int k = 1
        cdef int n = self.n_points()
        cdef int d = len(simplex)
        assert d == self.dim()+1
        cdef int i, j
        for i in range(1, d+1):
            l = simplex[i-1]+1
            for j in range(k, l):
                s += binomial(n-j, d-i)
            k = l+1
        return s

    cpdef int_to_simplex(self, int s):
        r"""
        Reverse the enumeration of possible simplices in
        :meth:`simplex_to_int`.

        The enumeration is compatible with [PUNTOS]_.

        INPUT:

        - ``s`` -- integer that uniquely specifies a simplex

        OUTPUT:

        An ordered tuple consisting of the indices of the vertices of
        the simplex.

        EXAMPLES::

            sage: U=matrix([
            ....:    [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ....:    [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ....:    [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ....:    [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ....:    [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ....: ])
            sage: pc = PointConfiguration(U.columns())
            sage: pc.simplex_to_int([1,3,4,7,10,13])
            1678
            sage: pc.int_to_simplex(1678)
            (1, 3, 4, 7, 10, 13)
        """
        simplex = []
        cdef int l = 0
        cdef int n = self.n_points()
        cdef int d = self.dim() + 1
        cdef int k, b
        for k in range(1, d):
            l += 1
            j = 1
            b = binomial(n - l, d - k)
            while s > b > 0:
                j += 1
                l += 1
                s -= b
                b = binomial(n - l, d - k)
            simplex.append(l - 1)
        simplex.append(s + l - 1)
        assert len(simplex) == d
        return tuple(simplex)


########################################################################
cdef class ConnectedTriangulationsIterator(SageObject):
    r"""
    A Python shim for the C++-class 'triangulations'.

    INPUT:

    - ``point_configuration`` -- a
      :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`

    - ``seed`` -- a regular triangulation or ``None`` (default). In
      the latter case, a suitable triangulation is generated
      automatically. Otherwise, you can explicitly specify the seed
      triangulation as

        * A
          :class:`~sage.geometry.triangulation.element.Triangulation`
          object, or

        * an iterable of iterables specifying the vertices of the simplices, or

        * an iterable of integers, which are then considered the
          enumerated simplices (see
          :meth:`~PointConfiguration_base.simplex_to_int`.

    - ``star`` -- either ``None`` (default) or an integer. If an
      integer is passed, all returned triangulations will be star with
      respect to the

    - ``fine`` -- boolean (default: ``False``); whether to return only
      fine triangulations, that is, simplicial decompositions that
      make use of all the points of the configuration.

    OUTPUT:

    An iterator. The generated values are tuples of
    integers, which encode simplices of the triangulation. The output
    is a suitable input to
    :class:`~sage.geometry.triangulation.element.Triangulation`.

    EXAMPLES::

        sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
        sage: from sage.geometry.triangulation.base import ConnectedTriangulationsIterator
        sage: ci = ConnectedTriangulationsIterator(p)
        sage: next(ci)
        (9, 10)
        sage: next(ci)
        (2, 3, 4, 5)
        sage: next(ci)
        (7, 8)
        sage: next(ci)
        (1, 3, 5, 7)
        sage: next(ci)
        Traceback (most recent call last):
        ...
        StopIteration

    You can reconstruct the triangulation from the compressed output via::

        sage: from sage.geometry.triangulation.element import Triangulation
        sage: Triangulation((2, 3, 4, 5), p)
        (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)

    How to use the restrictions::

        sage: ci = ConnectedTriangulationsIterator(p, fine=True)
        sage: list(ci)
        [(2, 3, 4, 5), (1, 3, 5, 7)]
        sage: ci = ConnectedTriangulationsIterator(p, star=1)
        sage: list(ci)
        [(7, 8)]
        sage: ci = ConnectedTriangulationsIterator(p, star=1, fine=True)
        sage: list(ci)
        []
    """

    cdef triangulations_ptr _tp

    def __cinit__(self):
        """
        The Cython constructor.

        TESTS::

            sage: from sage.geometry.triangulation.base import ConnectedTriangulationsIterator
            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: ConnectedTriangulationsIterator(p, fine=True)   # indirect doctest
            <sage.geometry.triangulation.base.ConnectedTriangulationsIterator object at ...>
        """
        self._tp = NULL

    def __init__(self, point_configuration, seed=None, star=None, fine=False):
        r"""
        The Python constructor.

        See :class:`ConnectedTriangulationsIterator` for a description
        of the arguments.

        TESTS::

            sage: p = PointConfiguration([[0,4],[2,3],[3,2],[4,0],[3,-2],[2,-3],[0,-4],[-2,-3],[-3,-2],[-4,0],[-3,2],[-2,3]])
            sage: from sage.geometry.triangulation.base import ConnectedTriangulationsIterator
            sage: ci = ConnectedTriangulationsIterator(p)
            sage: len(list(ci))  # long time (26s on sage.math, 2012)
            16796
            sage: ci = ConnectedTriangulationsIterator(p, star=3)
            sage: len(list(ci))  # long time (26s on sage.math, 2012)
            1
        """
        if star is None:
            star = -1
        if seed is None:
            seed = point_configuration.lexicographic_triangulation().enumerate_simplices()
        try:
            enumerated_simplices_seed = seed.enumerated_simplices()
        except AttributeError:
            enumerated_simplices_seed = tuple([int(t) for t in seed])
        assert self._tp == NULL
        self._tp = init_triangulations(point_configuration.n_points(),
                                       point_configuration.dim()+1,
                                       star, fine,
                                       enumerated_simplices_seed,
                                       point_configuration.bistellar_flips())

    def __dealloc__(self):
        r"""
        The Cython destructor.
        """
        delete_triangulations(self._tp)

    def __iter__(self):
        r"""
        The iterator interface: Start iterating.

        TESTS::

            sage: from sage.geometry.triangulation.base import ConnectedTriangulationsIterator
            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: ci = ConnectedTriangulationsIterator(p, fine=True)
            sage: ci.__iter__()
            <sage.geometry.triangulation.base.ConnectedTriangulationsIterator object at ...>
            sage: ci.__iter__() is ci
            True
        """
        return self

    def __next__(self):
        r"""
        The iterator interface: Next iteration.

        EXAMPLES::

            sage: from sage.geometry.triangulation.base import ConnectedTriangulationsIterator
            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: ci = ConnectedTriangulationsIterator(p)
            sage: ci.__next__()
            (9, 10)
        """
        t = next_triangulation(self._tp)
        if not t:
            raise StopIteration
        return t
