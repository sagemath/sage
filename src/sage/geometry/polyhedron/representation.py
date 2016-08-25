"""
H(yperplane) and V(ertex) representation objects for polyhedra
"""

#*****************************************************************************
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.sage_object import SageObject
from sage.structure.element import is_Vector
from sage.rings.all import QQ, ZZ, RDF
from sage.modules.free_module_element import vector



#########################################################################
#                      PolyhedronRepresentation
#                       /                     \
#                      /                       \
#              Hrepresentation            Vrepresentation
#                   /     \                 /     |   \
#                  /       \               /      |    \
#           Inequality  Equation        Vertex   Ray   Line

class PolyhedronRepresentation(SageObject):
    """
    The internal base class for all representation objects of
    ``Polyhedron`` (vertices/rays/lines and inequalites/equations)

    .. note::

        You should not (and cannot) instantiate it yourself. You can
        only obtain them from a Polyhedron() class.

    TESTS::

            sage: import sage.geometry.polyhedron.representation as P
            sage: P.PolyhedronRepresentation()
            <class 'sage.geometry.polyhedron.representation.PolyhedronRepresentation'>
    """

    # Numeric values for the output of the type() method
    INEQUALITY = 0
    EQUATION = 1
    VERTEX = 2
    RAY = 3
    LINE = 4

    def __len__(self):
        """
        Returns the length of the representation data.

        TESTS::

            sage: p = Polyhedron(vertices=[[1,2,3]])
            sage: v = p.Vrepresentation(0)
            sage: v.__len__()
            3
        """
        return self._vector.degree()

    def __getitem__(self, i):
        """
        Supports indexing.

        TESTS::

            sage: p = Polyhedron(vertices=[[1,2,3]])
            sage: v = p.Vrepresentation(0)
            sage: v.__getitem__(1)
            2
        """
        return self._vector[i]

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.geometry.polyhedron.representation import Hrepresentation
            sage: pr = Hrepresentation(Polyhedron(vertices = [[1,2,3]]).parent())
            sage: hash(pr)
            1647257843           # 32-bit
            4686581268940269811  # 64-bit
        """
        # TODO: ideally the argument self._vector of self should be immutable.
        # So that we could change the line below by hash(self._vector). The
        # mutability is kept because this argument might be reused (see e.g.
        # Hrepresentation._set_data below).
        return hash(tuple(self._vector))

    def __cmp__(self, other):
        """
        Compare two representation objects

        They are equal if and only if they define the same
        vertex/ray/line or inequality/equation in the ambient space,
        regardless of the polyhedron that they belong to.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        One of `-1`, `0`, `+1`.  ``True`` if and only if ``other`` represents the same
        H-representation object.

        EXAMPLES::

            sage: triangle = Polyhedron([(0,0), (1,0), (0,1)])
            sage: ieq = next(triangle.inequality_generator());  ieq
            An inequality (1, 0) x + 0 >= 0
            sage: ieq == copy(ieq)
            True
            sage: cmp(ieq, copy(ieq))
            0

            sage: cmp(ieq, 'a string')
            -1

            sage: square = Polyhedron([(0,0), (1,0), (0,1), (1,1)], base_ring=QQ)
            sage: cmp(square.Vrepresentation(0), triangle.Vrepresentation(0))
            0

            sage: ieq = square.Hrepresentation(0);  ieq.vector()
            (0, 1, 0)
            sage: abs(cmp(ieq, Polyhedron([(0,1,0)]).Vrepresentation(0)))
            1
        """
        if not isinstance(other, PolyhedronRepresentation):
            return -1
        return cmp(type(self), type(other)) or cmp(self._vector, other._vector)

    def vector(self, base_ring=None):
        """
        Returns the vector representation of the H/V-representation object.

        INPUT:

        - ``base_ring`` -- the base ring of the vector.

        OUTPUT:

        For a V-representation object, a vector of length
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.ambient_dim`. For
        a H-representation object, a vector of length
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.ambient_dim`
        + 1.

        EXAMPLES::

            sage: s = polytopes.cuboctahedron()
            sage: v = next(s.vertex_generator())
            sage: v
            A vertex at (-1, -1, 0)
            sage: v.vector()
            (-1, -1, 0)
            sage: v()
            (-1, -1, 0)
            sage: type(v())
            <type 'sage.modules.vector_integer_dense.Vector_integer_dense'>

       Conversion to a different base ring can be forced with the optional argument::

            sage: v.vector(RDF)
            (-1.0, -1.0, 0.0)
            sage: vector(RDF, v)
            (-1.0, -1.0, 0.0)
        """
        if (base_ring is None) or (base_ring is self._base_ring):
            return self._vector
        else:
            return vector(base_ring, self._vector)

    _vector_ = vector

    def polyhedron(self):
        """
        Returns the underlying polyhedron.

        TESTS::

            sage: p = Polyhedron(vertices=[[1,2,3]])
            sage: v = p.Vrepresentation(0)
            sage: v.polyhedron()
            A 0-dimensional polyhedron in ZZ^3 defined as the convex hull of 1 vertex
        """
        return self._polyhedron

    def __call__(self):
        """
        Returns the vector representation of the representation
        object. Shorthand for the vector() method.

        TESTS::

            sage: p = Polyhedron(vertices=[[1,2,3]])
            sage: v = p.Vrepresentation(0)
            sage: v.__call__()
            (1, 2, 3)
        """
        return self._vector

    def index(self):
        """
        Returns an arbitrary but fixed number according to the internal
        storage order.

        NOTES:

        H-representation and V-representation objects are enumerated
        independently. That is, amongst all vertices/rays/lines there
        will be one with ``index()==0``, and amongs all
        inequalities/equations there will be one with ``index()==0``,
        unless the polyhedron is empty or spans the whole space.

        EXAMPLES::

            sage: s = Polyhedron(vertices=[[1],[-1]])
            sage: first_vertex = next(s.vertex_generator())
            sage: first_vertex.index()
            0
            sage: first_vertex == s.Vrepresentation(0)
            True
        """
        return self._index

    def __add__(self, coordinate_list):
        """
        Return the coordinates concatenated with ``coordinate_list``.

        INPUT:

        - ``coordinate_list`` -- a list.

        OUTPUT:

        The coordinates of ``self`` concatenated with ``coordinate_list``.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: v = p.Vrepresentation(0); v
            A vertex at (1, 0)
            sage: v + [4,5]
            [1, 0, 4, 5]
        """
        if not isinstance(coordinate_list, list):
            raise TypeError('Can only concatenate with a list of coordinates')
        return list(self) + coordinate_list

    def __radd__(self, coordinate_list):
        """
        Return ``coordinate_list`` concatenated with the coordinates.

        INPUT:

        - ``coordinate_list`` -- a list.

        OUTPUT:

        ``coordinate_list`` concatenated with the coordinates of ``self``.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: v = p.Vrepresentation(0); v
            A vertex at (1, 0)
            sage: [4,5] + v
            [4, 5, 1, 0]
        """
        if not isinstance(coordinate_list, list):
            raise TypeError('Can only concatenate with a list of coordinates')
        return coordinate_list + list(self)

    def count(self, i):
        """
        Count the number of occurrences of ``i`` in the coordinates.

        INPUT:

        - ``i`` -- Anything.

        OUTPUT:

        Integer. The number of occurrences of ``i`` in the coordinates.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1,1,2,1)])
            sage: v = p.Vrepresentation(0); v
            A vertex at (0, 1, 1, 2, 1)
            sage: v.count(1)
            3
        """
        return sum([1 for j in self if i==j])


class Hrepresentation(PolyhedronRepresentation):
    """
    The internal base class for H-representation objects of
    a polyhedron. Inherits from ``PolyhedronRepresentation``.
    """

    def __init__(self, polyhedron_parent):
        """
        Initializes the PolyhedronRepresentation object.

        TESTS::

            sage: from sage.geometry.polyhedron.representation import Hrepresentation
            sage: pr = Hrepresentation(Polyhedron(vertices = [[1,2,3]]).parent())
            sage: tuple(pr)
            (0, 0, 0, 0)
            sage: TestSuite(pr).run(skip='_test_pickling')
        """
        self._polyhedron_parent = polyhedron_parent
        self._base_ring = polyhedron_parent.base_ring()
        self._vector = polyhedron_parent.Hrepresentation_space()(0)
        self._A = polyhedron_parent.ambient_space()(0)
        self._b = polyhedron_parent.base_ring()(0)
        self._index = 0

    def _set_data(self, polyhedron, data):
        """
        Initialization function.

        The H/V-representation objects are kept in a pool, and this
        function is used to reassign new values to already existing
        (but unused) objects. You must not call this function on
        objects that are in normal use.

        INPUT:

        - ``polyhedron`` -- the new polyhedron.

        - ``data`` -- the H-representation data.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0) # indirect doctest
            sage: TestSuite(pH).run(skip='_test_pickling')
        """
        assert polyhedron.parent() is self._polyhedron_parent
        if len(data) != self._vector.degree():
            raise ValueError('H-representation data requires a list of length ambient_dim+1')

        self._vector[:] = data
        self._A[:] = data[1:]
        self._b = self._base_ring(data[0])

        self._index = len(polyhedron._Hrepresentation)
        polyhedron._Hrepresentation.append(self)
        self._polyhedron = polyhedron

    def is_H(self):
        """
        Returns True if the object is part of a H-representation
        (inequality or equation).

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_H()
            True
        """
        return True

    def is_inequality(self):
        """
        Returns True if the object is an inequality of the H-representation.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_inequality()
            True
        """
        return False

    def is_equation(self):
        """
        Returns True if the object is an equation of the H-representation.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]], eqns = [[1,1,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_equation()
            True
        """
        return False

    def A(self):
        r"""
        Returns the coefficient vector `A` in `A\vec{x}+b`.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(2)
            sage: pH.A()
            (1, 0)
        """
        return self._A

    def b(self):
        r"""
        Returns the constant `b` in `A\vec{x}+b`.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(2)
            sage: pH.b()
            0
        """
        return self._b

    def neighbors(self):
        """
        Iterate over the adjacent facets (i.e. inequalities/equations)

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],
            ...                          [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: a = list(pH.neighbors())
            sage: a[0]
            An inequality (0, -1, 0) x + 1 >= 0
            sage: list(a[0])
            [1, 0, -1, 0]
        """
        adjacency_matrix = self.polyhedron().facet_adjacency_matrix()
        for x in self.polyhedron().Hrep_generator():
            if adjacency_matrix[self.index(), x.index()] == 1:
                yield x

    def adjacent(self):
        """
        Alias for neighbors().

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,0,2],[0,0,1,0,],[0,10,0,0],
            ...       [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: a = list(pH.neighbors())
            sage: b = list(pH.adjacent())
            sage: a==b
            True
        """
        return self.neighbors()

    def is_incident(self, Vobj):
        """
        Returns whether the incidence matrix element (Vobj,self) == 1

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],
            ...       [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_incident(p.Vrepresentation(1))
            True
            sage: pH.is_incident(p.Vrepresentation(5))
            False
        """
        return self.polyhedron().incidence_matrix()[Vobj.index(), self.index()] == 1

    def __mul__(self, Vobj):
        """
        Shorthand for ``self.eval(x)``

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],
            ...        [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH*p.Vrepresentation(5)
            1
        """
        return self.eval(Vobj)

    def eval(self, Vobj):
        r"""
        Evaluates the left hand side `A\vec{x}+b` on the given
        vertex/ray/line.

        NOTES:

          * Evaluating on a vertex returns `A\vec{x}+b`
          * Evaluating on a ray returns `A\vec{r}`. Only the sign or
            whether it is zero is meaningful.
          * Evaluating on a line returns `A\vec{l}`. Only whether it
            is zero or not is meaningful.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[-1,-1]])
            sage: ineq = next(triangle.inequality_generator())
            sage: ineq
            An inequality (2, -1) x + 1 >= 0
            sage: [ ineq.eval(v) for v in triangle.vertex_generator() ]
            [0, 0, 3]
            sage: [ ineq * v for v in triangle.vertex_generator() ]
            [0, 0, 3]

        If you pass a vector, it is assumed to be the coordinate vector of a point::

            sage: ineq.eval( vector(ZZ, [3,2]) )
            5
        """
        if is_Vector(Vobj):
            return self.A() * Vobj + self.b()
        return Vobj.evaluated_on(self)

    def incident(self):
        """
        Returns a generator for the incident H-representation objects,
        that is, the vertices/rays/lines satisfying the (in)equality.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[-1,-1]])
            sage: ineq = next(triangle.inequality_generator())
            sage: ineq
            An inequality (2, -1) x + 1 >= 0
            sage: [ v for v in ineq.incident()]
            [A vertex at (-1, -1), A vertex at (0, 1)]
            sage: p = Polyhedron(vertices=[[0,0,0],[0,1,0],[0,0,1]], rays=[[1,-1,-1]])
            sage: ineq = p.Hrepresentation(2)
            sage: ineq
            An inequality (1, 0, 1) x + 0 >= 0
            sage: [ x for x in ineq.incident() ]
            [A vertex at (0, 0, 0),
             A vertex at (0, 1, 0),
             A ray in the direction (1, -1, -1)]
        """
        incidence_matrix = self.polyhedron().incidence_matrix()
        for V in self.polyhedron().Vrep_generator():
            if incidence_matrix[V.index(), self.index()] == 1:
                yield V


class Inequality(Hrepresentation):
    """
    A linear inequality (supporting hyperplane) of the
    polyhedron. Inherits from ``Hrepresentation``.
    """

    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: repr_obj = next(p.inequality_generator())
            sage: repr_obj.type()
            0
            sage: repr_obj.type() == repr_obj.INEQUALITY
            True
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.INEQUALITY


    def is_inequality(self):
        """
        Returns True since this is, by construction, an inequality.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = next(p.inequality_generator())
            sage: a.is_inequality()
            True
        """
        return True

    def _repr_(self):
        """
        The string representation of the inequality.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = next(p.inequality_generator())
            sage: a._repr_()
            'An inequality (-1, 1, 0) x + 0 >= 0'
            sage: Polyhedron(ieqs=[(1,-1),(-1,2)]).Hrepresentation()
            (An inequality (-1) x + 1 >= 0, An inequality (2) x - 1 >= 0)
            sage: Polyhedron(eqns=[(1,0)]).Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(eqns=[(-1,0)]).Hrepresentation()
            (An equation -1 == 0,)

        TESTS:

        Test that :trac:`21105` has been fixed::

            sage: K.<cbrt2> = NumberField(x^3 - 2, 'a', embedding=1.26)
            sage: P = Polyhedron(vertices=[(1,1,cbrt2),(cbrt2,1,1)])
            sage: P.inequalities()
            (An inequality (-cbrt2^2 - cbrt2 - 1, 0, 0) x + cbrt2^2 + cbrt2 + 2 >= 0,
             An inequality (cbrt2^2 + cbrt2 + 1, 0, 0) x - cbrt2^2 + cbrt2 + 1 >= 0)
        """
        s = 'An inequality '
        have_A = not self.A().is_zero()
        if have_A:
            s += repr(self.A()) + ' x '
        if self.b()>=0:
            if have_A:
                s += '+'
        else:
            s += '-'
        if have_A:
            s += ' '
        s += repr(abs(self.b())) + ' >= 0'
        return s

    def contains(self, Vobj):
        """
        Tests whether the halfspace (including its boundary) defined
        by the inequality contains the given vertex/ray/line.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: i1 = next(p.inequality_generator())
            sage: [i1.contains(q) for q in p.vertex_generator()]
            [True, True, True, True, True, True]
            sage: p2 = 3*polytopes.hypercube(3)
            sage: [i1.contains(q) for q in p2.vertex_generator()]
            [True, False, False, False, True, True, True, False]
        """
        try:
            if Vobj.is_vector(): # assume we were passed a point
                return self.polyhedron()._is_nonneg( self.eval(Vobj) )
        except AttributeError:
            pass

        if Vobj.is_line():
            return self.polyhedron()._is_zero( self.eval(Vobj) )
        else:
            return self.polyhedron()._is_nonneg( self.eval(Vobj) )

    def interior_contains(self, Vobj):
        """
        Tests whether the interior of the halfspace (excluding its
        boundary) defined by the inequality contains the given
        vertex/ray/line.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: i1 = next(p.inequality_generator())
            sage: [i1.interior_contains(q) for q in p.vertex_generator()]
            [False, True, True, False, False, True]
            sage: p2 = 3*polytopes.hypercube(3)
            sage: [i1.interior_contains(q) for q in p2.vertex_generator()]
            [True, False, False, False, True, True, True, False]

        If you pass a vector, it is assumed to be the coordinate vector of a point::

            sage: P = Polyhedron(vertices=[[1,1],[1,-1],[-1,1],[-1,-1]])
            sage: p = vector(ZZ, [1,0] )
            sage: [ ieq.interior_contains(p) for ieq in P.inequality_generator() ]
            [True, True, False, True]
        """
        try:
            if Vobj.is_vector(): # assume we were passed a point
                return self.polyhedron()._is_positive( self.eval(Vobj) )
        except AttributeError:
            pass

        if Vobj.is_line():
            return self.polyhedron()._is_zero( self.eval(Vobj) )
        elif Vobj.is_vertex():
            return self.polyhedron()._is_positive( self.eval(Vobj) )
        else: # Vobj.is_ray()
            return self.polyhedron()._is_nonneg( self.eval(Vobj) )


class Equation(Hrepresentation):
    """
    A linear equation of the polyhedron. That is, the polyhedron is
    strictly smaller-dimensional than the ambient space, and contained
    in this hyperplane. Inherits from ``Hrepresentation``.
    """

    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: repr_obj = next(p.equation_generator())
            sage: repr_obj.type()
            1
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            True
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.EQUATION


    def is_equation(self):
        """
        Tests if this object is an equation.  By construction, it must be.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = next(p.equation_generator())
            sage: a.is_equation()
            True
        """
        return True

    def _repr_(self):
        """
        A string representation of this object.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = next(p.equation_generator())
            sage: a._repr_()
            'An equation (0, 0, 1) x + 0 == 0'
            sage: Polyhedron().Hrepresentation(0)
            An equation -1 == 0
        """
        s = 'An equation '
        have_A = not self.A().is_zero()
        if have_A:
            s += repr(self.A()) + ' x '
        if self.b()>=0:
            if have_A:
                s += '+'
        else:
            s += '-'
        if have_A:
            s += ' '
        s += repr(abs(self.b())) + ' == 0'
        return s

    def contains(self, Vobj):
        """
        Tests whether the hyperplane defined by the equation contains
        the given vertex/ray/line.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: v = next(p.vertex_generator())
            sage: v
            A vertex at (0, 0, 0)
            sage: a = next(p.equation_generator())
            sage: a
            An equation (0, 0, 1) x + 0 == 0
            sage: a.contains(v)
            True
        """
        return self.polyhedron()._is_zero( self.eval(Vobj) )

    def interior_contains(self, Vobj):
        """
        Tests whether the interior of the halfspace (excluding its
        boundary) defined by the inequality contains the given
        vertex/ray/line.

        NOTE:

        Returns False for any equation.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: v = next(p.vertex_generator())
            sage: v
            A vertex at (0, 0, 0)
            sage: a = next(p.equation_generator())
            sage: a
            An equation (0, 0, 1) x + 0 == 0
            sage: a.interior_contains(v)
            False
        """
        return False


class Vrepresentation(PolyhedronRepresentation):
    """
    The base class for V-representation objects of a
    polyhedron. Inherits from ``PolyhedronRepresentation``.
    """

    def __init__(self, polyhedron_parent):
        """
        Initializes the PolyhedronRepresentation object.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = next(p.inequality_generator())
            sage: a
            An inequality (-1, 1, 0) x + 0 >= 0
            sage: TestSuite(a).run(skip='_test_pickling')
        """
        self._polyhedron_parent = polyhedron_parent
        self._base_ring = polyhedron_parent.base_ring()
        self._vector = polyhedron_parent.Vrepresentation_space()(0)
        self._index = 0

    def _set_data(self, polyhedron, data):
        """
        Initialization function.

        The H/V-representation objects are kept in a pool, and this
        function is used to reassign new values to already existing
        (but unused) objects. You must not call this function on
        objects that are in normal use.

        INPUT:

        - ``polyhedron`` -- the new polyhedron.

        - ``data`` -- the V-representation data.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pV = p.Vrepresentation(0)   # indirect doctest
            sage: TestSuite(pV).run(skip='_test_pickling')
        """
        assert polyhedron.parent() is self._polyhedron_parent
        if len(data) != self._vector.degree():
            raise ValueError('V-representation data requires a list of length ambient_dim')

        self._vector[:] = data

        self._index = len(polyhedron._Vrepresentation)
        polyhedron._Vrepresentation.append(self)
        self._polyhedron = polyhedron

    def is_V(self):
        """
        Returns True if the object is part of a V-representation
        (a vertex, ray, or line).

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,3]])
            sage: v = next(p.vertex_generator())
            sage: v.is_V()
            True
        """
        return True

    def is_vertex(self):
        """
        Returns True if the object is a vertex of the V-representation.
        This method is over-ridden by the corresponding method in the
        derived class Vertex.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,3]])
            sage: v = next(p.vertex_generator())
            sage: v.is_vertex()
            True
            sage: p = Polyhedron(ieqs = [[1, 0, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0]])
            sage: r1 = next(p.ray_generator())
            sage: r1.is_vertex()
            False
        """
        return False

    def is_ray(self):
        """
        Returns True if the object is a ray of the V-representation.
        This method is over-ridden by the corresponding method in the
        derived class Ray.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0]])
            sage: r1 = next(p.ray_generator())
            sage: r1.is_ray()
            True
            sage: v1 = next(p.vertex_generator())
            sage: v1
            A vertex at (-1, -1, 0, -1)
            sage: v1.is_ray()
            False
        """
        return False

    def is_line(self):
        """
        Returns True if the object is a line of the V-representation.
        This method is over-ridden by the corresponding method in the
        derived class Line.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0]])
            sage: line1 = next(p.line_generator())
            sage: line1.is_line()
            True
            sage: v1 = next(p.vertex_generator())
            sage: v1.is_line()
            False
        """
        return False

    def neighbors(self):
        """
        Returns a generator for the adjacent vertices/rays/lines.

        EXAMPLES::

             sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,4]])
             sage: v = next(p.vertex_generator())
             sage: next(v.neighbors())
             A vertex at (0, 3)
        """
        adjacency_matrix = self.polyhedron().vertex_adjacency_matrix()
        for x in self.polyhedron().Vrep_generator():
            if adjacency_matrix[self.index(), x.index()] == 1:
                yield x

    def adjacent(self):
        """
        Alias for neighbors().

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,4]])
            sage: v = next(p.vertex_generator())
            sage: a = next(v.neighbors())
            sage: b = next(v.adjacent())
            sage: a==b
            True
        """
        return self.neighbors()

    def is_incident(self, Hobj):
        """
        Returns whether the incidence matrix element (self,Hobj) == 1

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: h1 = next(p.inequality_generator())
            sage: h1
            An inequality (0, 0, -1) x + 1 >= 0
            sage: v1 = next(p.vertex_generator())
            sage: v1
            A vertex at (-1, -1, -1)
            sage: v1.is_incident(h1)
            False
        """
        return self.polyhedron().incidence_matrix()[self.index(), Hobj.index()] == 1

    def __mul__(self, Hobj):
        """
        Shorthand for self.evaluated_on(Hobj)

        TESTS::

            sage: p = polytopes.hypercube(3)
            sage: h1 = next(p.inequality_generator())
            sage: v1 = next(p.vertex_generator())
            sage: v1.__mul__(h1)
            2
        """
        return self.evaluated_on(Hobj)

    def incident(self):
        """
        Returns a generator for the equations/inequalities that are satisfied on the given
        vertex/ray/line.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[-1,-1]])
            sage: ineq = next(triangle.inequality_generator())
            sage: ineq
            An inequality (2, -1) x + 1 >= 0
            sage: [ v for v in ineq.incident()]
            [A vertex at (-1, -1), A vertex at (0, 1)]
            sage: p = Polyhedron(vertices=[[0,0,0],[0,1,0],[0,0,1]], rays=[[1,-1,-1]])
            sage: ineq = p.Hrepresentation(2)
            sage: ineq
            An inequality (1, 0, 1) x + 0 >= 0
            sage: [ x for x in ineq.incident() ]
            [A vertex at (0, 0, 0),
             A vertex at (0, 1, 0),
             A ray in the direction (1, -1, -1)]
        """
        incidence_matrix = self.polyhedron().incidence_matrix()
        for H in self.polyhedron().Hrep_generator():
            if incidence_matrix[self.index(), H.index()] == 1:
                yield H


class Vertex(Vrepresentation):
    """
    A vertex of the polyhedron. Inherits from ``Vrepresentation``.
    """

    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: repr_obj = next(p.vertex_generator())
            sage: repr_obj.type()
            2
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            True
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.VERTEX

    def is_vertex(self):
        """
        Tests if this object is a vertex.  By construction it always is.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = next(p.vertex_generator())
            sage: a.is_vertex()
            True
        """
        return True

    def _repr_(self):
        """
        Returns a string representation of the vertex.

        OUTPUT:

        String.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: v = next(p.vertex_generator())
            sage: v.__repr__()
            'A vertex at (1, 0)'
        """
        return 'A vertex at ' + repr(self.vector());

    def homogeneous_vector(self, base_ring=None):
        """
        Return homogeneous coordinates for this vertex.

        Since a vertex is given by an affine point, this is the vector
        with a 1 appended.

        INPUT:

        - ``base_ring`` -- the base ring of the vector.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(2,0)], rays=[(1,0)], lines=[(3,2)])
            sage: P.vertices()[0].homogeneous_vector()
            (2, 0, 1)
            sage: P.vertices()[0].homogeneous_vector(RDF)
            (2.0, 0.0, 1.0)
        """
        v = list(self._vector) + [1]
        return vector(base_ring or self._base_ring, v)

    def evaluated_on(self, Hobj):
        r"""
        Returns `A\vec{x}+b`

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: v = next(p.vertex_generator())
            sage: h = next(p.inequality_generator())
            sage: v
            A vertex at (-1, -1, -1)
            sage: h
            An inequality (0, 0, -1) x + 1 >= 0
            sage: v.evaluated_on(h)
            2
        """
        return Hobj.A() * self.vector() + Hobj.b()

    def is_integral(self):
        r"""
        Return whether the coordinates of the vertex are all integral.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(1/2,3,5), (0,0,0), (2,3,7)])
            sage: [ v.is_integral() for v in p.vertex_generator() ]
            [True, False, True]
        """
        return (self._base_ring is ZZ) or all(x in ZZ for x in self)


class Ray(Vrepresentation):
    """
    A ray of the polyhedron. Inherits from ``Vrepresentation``.
    """

    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: repr_obj = next(p.ray_generator())
            sage: repr_obj.type()
            3
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            True
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.RAY

    def is_ray(self):
        """
        Tests if this object is a ray.  Always True by construction.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = next(p.ray_generator())
            sage: a.is_ray()
            True
        """
        return True

    def _repr_(self):
        """
        A string representation of the ray.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = next(p.ray_generator())
            sage: a._repr_()
            'A ray in the direction (0, 1)'
        """
        return 'A ray in the direction ' + repr(self.vector());

    def homogeneous_vector(self, base_ring=None):
        """
        Return homogeneous coordinates for this ray.

        Since a ray is given by a direction, this is the vector with a
        0 appended.

        INPUT:

        - ``base_ring`` -- the base ring of the vector.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(2,0)], rays=[(1,0)], lines=[(3,2)])
            sage: P.rays()[0].homogeneous_vector()
            (1, 0, 0)
            sage: P.rays()[0].homogeneous_vector(RDF)
            (1.0, 0.0, 0.0)
        """
        v = list(self._vector) + [0]
        return vector(base_ring or self._base_ring, v)

    def evaluated_on(self, Hobj):
        r"""
        Returns `A\vec{r}`

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = next(p.ray_generator())
            sage: h = next(p.inequality_generator())
            sage: a.evaluated_on(h)
            0
        """
        return Hobj.A() * self.vector()


class Line(Vrepresentation):
    r"""
    A line (Minkowski summand `\simeq\RR`) of the
    polyhedron. Inherits from ``Vrepresentation``.
    """

    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: repr_obj = next(p.line_generator())
            sage: repr_obj.type()
            4
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            True
        """
        return self.LINE

    def is_line(self):
        """
        Tests if the object is a line.  By construction it must be.

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = next(p.line_generator())
            sage: a.is_line()
            True
        """
        return True

    def _repr_(self):
        """
        A string representation of the line.

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = next(p.line_generator())
            sage: a.__repr__()
            'A line in the direction (0, 1, 0)'
        """
        return 'A line in the direction ' + repr(self.vector());

    def homogeneous_vector(self, base_ring=None):
        """
        Return homogeneous coordinates for this line.

        Since a line is given by a direction, this is the vector with a
        0 appended.

        INPUT:

        - ``base_ring`` -- the base ring of the vector.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(2,0)], rays=[(1,0)], lines=[(3,2)])
            sage: P.lines()[0].homogeneous_vector()
            (3, 2, 0)
            sage: P.lines()[0].homogeneous_vector(RDF)
            (3.0, 2.0, 0.0)
        """
        v = list(self._vector) + [0]
        return vector(base_ring or self._base_ring, v)

    def evaluated_on(self, Hobj):
        r"""
        Returns `A\vec{\ell}`

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = next(p.line_generator())
            sage: h = next(p.inequality_generator())
            sage: a.evaluated_on(h)
            0
        """
        return Hobj.A() * self.vector()
