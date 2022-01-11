"""
H(yperplane) and V(ertex) representation objects for polyhedra
"""

# ****************************************************************************
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.sage_object import SageObject
from sage.structure.element import is_Vector
from sage.structure.richcmp import richcmp_method, richcmp
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from copy import copy


#########################################################################
#                      PolyhedronRepresentation
#                       /                     \
#                      /                       \
#              Hrepresentation            Vrepresentation
#                   /     \                 /     |   \
#                  /       \               /      |    \
#           Inequality  Equation        Vertex   Ray   Line


@richcmp_method
class PolyhedronRepresentation(SageObject):
    """
    The internal base class for all representation objects of
    ``Polyhedron`` (vertices/rays/lines and inequalities/equations)

    .. note::

        You should not (and cannot) instantiate it yourself. You can
        only obtain them from a Polyhedron() class.

    TESTS::

        sage: import sage.geometry.polyhedron.representation as P
        sage: P.PolyhedronRepresentation()
        <sage.geometry.polyhedron.representation.PolyhedronRepresentation object at ...>
    """

    # Numeric values for the output of the type() method
    INEQUALITY = 0
    EQUATION = 1
    VERTEX = 2
    RAY = 3
    LINE = 4

    def __len__(self):
        """
        Return the length of the representation data.

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
            sage: hash(pr) == hash(tuple([0,0,0,0]))
            True
        """
        # TODO: ideally the argument self._vector of self should be immutable.
        # So that we could change the line below by hash(self._vector). The
        # mutability is kept because this argument might be reused (see e.g.
        # Hrepresentation._set_data below).
        return hash(tuple(self._vector))

    def __richcmp__(self, other, op):
        """
        Compare two representation objects

        This method defines a linear order on the H/V-representation objects.
        The order is first determined by the types of the objects,
        such that inequality < equation < vertex < ray < line.
        Then, representation objects with the same type are ordered
        lexicographically according to their canonical vectors.

        Thus, two representation objects are equal if and only if they define
        the same vertex/ray/line or inequality/equation in the ambient space,
        regardless of the polyhedron that they belong to.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        boolean

        EXAMPLES::

            sage: triangle = Polyhedron([(0,0), (1,0), (0,1)])
            sage: ieq = next(triangle.inequality_generator());  ieq
            An inequality (1, 0) x + 0 >= 0
            sage: ieq == copy(ieq)
            True

            sage: square = Polyhedron([(0,0), (1,0), (0,1), (1,1)], base_ring=QQ)
            sage: square.Vrepresentation(0) == triangle.Vrepresentation(0)
            True

            sage: ieq = square.Hrepresentation(0); ieq.vector()
            (0, 1, 0)
            sage: ieq != Polyhedron([(0,1,0)]).Vrepresentation(0)
            True

            sage: H = Polyhedron(vertices=[(4,0)], rays=[(1,1)], lines=[(-1,1)])
            sage: H.vertices()[0] < H.rays()[0] < H.lines()[0]
            True

        TESTS:

        Check :trac:`30954`::

            sage: P = (1/2)*polytopes.cube()
            sage: Q = (1/2)*polytopes.cube(backend='field')
            sage: for p in P.inequalities():
            ....:     assert p in Q.inequalities()
        """
        if not isinstance(other, PolyhedronRepresentation):
            return NotImplemented
        return richcmp((self.type(), self._vector*self._comparison_scalar()),
                (other.type(), other._vector*other._comparison_scalar()), op)

    def _comparison_scalar(self):
        r"""
        Return a number ``a`` such that ``a*self._vector`` is canonical.

        Except for vertices, ``self._vector`` is only unique up to a positive scalar.

        This is overwritten for the vertex class.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0],[1,5]], rays=[[3,4]])
            sage: P.Vrepresentation()
            (A vertex at (0, 0), A vertex at (1, 5), A ray in the direction (3, 4))
            sage: P.Vrepresentation()[0]._comparison_scalar()
            1
            sage: P.Vrepresentation()[1]._comparison_scalar()
            1
            sage: P.Vrepresentation()[2]._comparison_scalar()
            1/4
            sage: P.Hrepresentation()
            (An inequality (5, -1) x + 0 >= 0,
             An inequality (-4, 3) x + 0 >= 0,
             An inequality (4, -3) x + 11 >= 0)
            sage: P.Hrepresentation()[0]._comparison_scalar()
            1
            sage: P.Hrepresentation()[1]._comparison_scalar()
            1/3
            sage: P.Hrepresentation()[2]._comparison_scalar()
            1/3

        ::

            sage: P = Polyhedron(vertices=[[1,3]], lines=[[-1,3]])
            sage: P.Vrepresentation()
            (A line in the direction (1, -3), A vertex at (2, 0))
            sage: P.Vrepresentation()[0]._comparison_scalar()
            -1/3
            sage: P.Vrepresentation()[1]._comparison_scalar()
            1
        """
        if self.type() == self.VERTEX:
            return 1

        lcf = self._vector.leading_coefficient()
        if self.type() == self.EQUATION or self.type() == self.LINE:
            return 1/lcf
        else:
            return 1/lcf.abs()

    def vector(self, base_ring=None):
        """
        Return the vector representation of the H/V-representation object.

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
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

       Conversion to a different base ring can be forced with the optional argument::

            sage: v.vector(RDF)
            (-1.0, -1.0, 0.0)
            sage: vector(RDF, v)
            (-1.0, -1.0, 0.0)

        TESTS:

        Checks that :trac:`27709` is fixed::

            sage: C = polytopes.cube()
            sage: C.vertices()[0].vector()[0] = 3
            sage: C.vertices()[0]
            A vertex at (1, -1, -1)
        """
        if (base_ring is None) or (base_ring is self._base_ring):
            return copy(self._vector)
        else:
            return vector(base_ring, self._vector)

    _vector_ = vector

    def polyhedron(self):
        """
        Return the underlying polyhedron.

        TESTS::

            sage: p = Polyhedron(vertices=[[1,2,3]])
            sage: v = p.Vrepresentation(0)
            sage: v.polyhedron()
            A 0-dimensional polyhedron in ZZ^3 defined as the convex hull of 1 vertex
        """
        return self._polyhedron

    def __call__(self):
        """
        Return the vector representation of the representation
        object. Shorthand for the vector() method.

        TESTS::

            sage: p = Polyhedron(vertices=[[1,2,3]])
            sage: v = p.Vrepresentation(0)
            sage: v.__call__()
            (1, 2, 3)
        """
        return copy(self._vector)

    def index(self):
        """
        Return an arbitrary but fixed number according to the internal
        storage order.

        .. NOTE::

            H-representation and V-representation objects are enumerated
            independently. That is, amongst all vertices/rays/lines there
            will be one with ``index()==0``, and amongst all
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
        return sum([1 for j in self if i == j])


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
        if polyhedron.is_mutable():
            polyhedron._add_dependent_object(self)

    def is_H(self):
        """
        Return True if the object is part of a H-representation
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
        Return True if the object is an inequality of the H-representation.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_inequality()
            True
        """
        return False

    def is_equation(self):
        """
        Return True if the object is an equation of the H-representation.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]], eqns = [[1,1,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_equation()
            True
        """
        return False

    def A(self):
        r"""
        Return the coefficient vector `A` in `A\vec{x}+b`.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(2)
            sage: pH.A()
            (1, 0)

        TESTS:

        Checks that :trac:`27709` is fixed::

            sage: C = polytopes.cube()
            sage: C.inequalities()[0].A()[2] = 5
            sage: C.inequalities()[0]
            An inequality (-1, 0, 0) x + 1 >= 0
        """
        return copy(self._A)

    def b(self):
        r"""
        Return the constant `b` in `A\vec{x}+b`.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(2)
            sage: pH.b()
            0
        """
        return self._b

    def neighbors(self):
        """
        Iterate over the adjacent facets (i.e. inequalities).

        Only defined for inequalities.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],
            ....:                        [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: a = list(pH.neighbors())
            sage: a[0]
            An inequality (0, -1, 0) x + 1 >= 0
            sage: list(a[0])
            [1, 0, -1, 0]

        TESTS:

        Checking that :trac:`28463` is fixed::

            sage: P = polytopes.simplex()
            sage: F1 = P.Hrepresentation()[1]
            sage: list(F1.neighbors())
            [An inequality (0, 1, 0, 0) x + 0 >= 0,
             An inequality (0, 0, 1, 0) x + 0 >= 0,
             An inequality (0, 0, 0, 1) x + 0 >= 0]

        Does not work for equalities::

            sage: F0 = P.Hrepresentation()[0]
            sage: list(F0.neighbors())
            Traceback (most recent call last):
            ...
            TypeError: must be inequality
        """
        # The adjacency matrix does not include equations.
        n_eqs = self.polyhedron().n_equations()
        if not self.is_inequality():
            raise TypeError("must be inequality")

        adjacency_matrix = self.polyhedron().facet_adjacency_matrix()
        for x in self.polyhedron().Hrep_generator():
            if not x.is_equation():
                if adjacency_matrix[self.index()-n_eqs, x.index()-n_eqs] == 1:
                    yield x

    def adjacent(self):
        """
        Alias for neighbors().

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,0,2],[0,0,1,0,],[0,10,0,0],
            ....:     [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: a = list(pH.neighbors())
            sage: b = list(pH.adjacent())
            sage: a==b
            True
        """
        return self.neighbors()

    def is_incident(self, Vobj):
        """
        Return whether the incidence matrix element (Vobj,self) == 1

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],
            ....:     [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
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
            ....:      [1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH*p.Vrepresentation(5)
            1
        """
        return self.eval(Vobj)

    def eval(self, Vobj):
        r"""
        Evaluate the left hand side `A\vec{x}+b` on the given
        vertex/ray/line.

        .. NOTE:

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
        Return a generator for the incident H-representation objects,
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

    def repr_pretty(self, **kwds):
        r"""
        Return a pretty representation of this equality/inequality.

        INPUT:

        - ``prefix`` -- a string

        - ``indices`` -- a tuple or other iterable

        - ``latex`` -- a boolean

        OUTPUT:

        A string

        EXAMPLES::

            sage: P = Polyhedron(ieqs=[(0, 1, 0, 0), (1, 2, 1, 0)],
            ....:                eqns=[(1, -1, -1, 1)])
            sage: for h in P.Hrepresentation():
            ....:     print(h.repr_pretty())
            x0 + x1 - x2 == 1
            x0 >= 0
            2*x0 + x1 >= -1
        """
        return repr_pretty(self.vector(), self.type(), **kwds)

    def _latex_(self):
        r"""
        Return a LaTeX-representation of this equality/inequality.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P = Polyhedron(ieqs=[(0, 1, 0, 0), (1, 2, 1, 0)],
            ....:                eqns=[(1, -1, -1, 1)])
            sage: for h in P.Hrepresentation():
            ....:     print(latex(h))
            x_{0} + x_{1} - x_{2} = 1
            x_{0} \geq 0
            2 x_{0} + x_{1} \geq -1
        """
        return self.repr_pretty(latex=True)


class Inequality(Hrepresentation):
    """
    A linear inequality (supporting hyperplane) of the
    polyhedron. Inherits from ``Hrepresentation``.
    """

    def type(self):
        r"""
        Return the type (equation/inequality/vertex/ray/line) as an
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
        Return True since this is, by construction, an inequality.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = next(p.inequality_generator())
            sage: a.is_inequality()
            True
        """
        return True

    def is_facet_defining_inequality(self, other):
        r"""
        Check if ``self`` defines a facet of ``other``.

        INPUT:

        - ``other`` -- a polyhedron

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.slack_matrix`
            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0,0],[0,1,0]], rays=[[1,0,0]])
            sage: P.inequalities()
            (An inequality (1, 0, 0) x + 0 >= 0,
             An inequality (0, 1, 0) x + 0 >= 0,
             An inequality (0, -1, 0) x + 1 >= 0)
            sage: Q = Polyhedron(ieqs=[[0,1,0,0]])
            sage: Q.inequalities()[0].is_facet_defining_inequality(P)
            True
            sage: Q = Polyhedron(ieqs=[[0,2,0,3]])
            sage: Q.inequalities()[0].is_facet_defining_inequality(P)
            True
            sage: Q = Polyhedron(ieqs=[[0,AA(2).sqrt(),0,3]])                   # optional - sage.rings.number_field
            sage: Q.inequalities()[0].is_facet_defining_inequality(P)           # optional - sage.rings.number_field
            True
            sage: Q = Polyhedron(ieqs=[[1,1,0,0]])
            sage: Q.inequalities()[0].is_facet_defining_inequality(P)
            False

        ::

            sage: P = Polyhedron(vertices=[[0,0,0],[0,1,0]], lines=[[1,0,0]])
            sage: P.inequalities()
            (An inequality (0, 1, 0) x + 0 >= 0, An inequality (0, -1, 0) x + 1 >= 0)
            sage: Q = Polyhedron(ieqs=[[0,1,0,0]])
            sage: Q.inequalities()[0].is_facet_defining_inequality(P)
            False
            sage: Q = Polyhedron(ieqs=[[0,-1,0,0]])
            sage: Q.inequalities()[0].is_facet_defining_inequality(P)
            False
            sage: Q = Polyhedron(ieqs=[[0,0,1,3]])
            sage: Q.inequalities()[0].is_facet_defining_inequality(P)
            True

        TESTS::

            sage: p1 = Polyhedron(backend='normaliz', base_ring=QQ, vertices=[  # optional - pynormaliz
            ....:     (2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3),
            ....:     (1, 1, 1, 9/10, 4/5, 7/10, 3/5, 0, 0),
            ....:     (1, 1, 1, 1, 4/5, 3/5, 1/2, 1/10, 0),
            ....:     (1, 1, 1, 1, 9/10, 1/2, 2/5, 1/5, 0),
            ....:     (1, 1, 1, 1, 1, 2/5, 3/10, 1/5, 1/10)])
            sage: p2 = Polyhedron(backend='ppl', base_ring=QQ, vertices=[
            ....:     (2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3),
            ....:     (1, 1, 1, 9/10, 4/5, 7/10, 3/5, 0, 0),
            ....:     (1, 1, 1, 1, 4/5, 3/5, 1/2, 1/10, 0),
            ....:     (1, 1, 1, 1, 9/10, 1/2, 2/5, 1/5, 0),
            ....:     (1, 1, 1, 1, 1, 2/5, 3/10, 1/5, 1/10)])
            sage: p2 == p1                                                      # optional - pynormaliz
            True
            sage: for ieq in p1.inequalities():                                 # optional - pynormaliz
            ....:     assert ieq.is_facet_defining_inequality(p2)
            sage: for ieq in p2.inequalities():                                 # optional - pynormaliz
            ....:     assert ieq.is_facet_defining_inequality(p1)
        """
        from sage.geometry.polyhedron.base import Polyhedron_base
        if not isinstance(other, Polyhedron_base):
            raise ValueError("other must be a polyhedron")

        if not other.n_Vrepresentation():
            # An empty polytope does not have facets.
            return False

        # We evaluate ``self`` on the Vrepresentation of other.

        from sage.matrix.constructor import matrix
        Vrep_matrix = matrix(other.base_ring(), other.Vrepresentation())

        # Getting homogeneous coordinates of the Vrepresentation.
        hom_helper = matrix(other.base_ring(), [1 if v.is_vertex() else 0 for v in other.Vrepresentation()])
        hom_Vrep = hom_helper.stack(Vrep_matrix.transpose())

        self_matrix = matrix(self.vector())

        cross_slack_matrix = self_matrix * hom_Vrep

        # First of all ``self`` should not evaluate negative on anything.
        # If it has the same incidences as an inequality of ``other``,
        # all ``Vrepresentatives`` lie on the same (closed) side.
        if not any(x > 0 for x in cross_slack_matrix):
            return False

        # Also it should evaluate ``0`` on all lines.
        if any(self.A() * line.vector() for line in other.lines()):
            return False

        incidences = cross_slack_matrix.zero_pattern_matrix(ZZ)

        # See if ``self`` has the same incidences as an inequality of ``other``.
        # If this is the case, then the above check suffices to guarantee that all
        # entries of ``cross_slack_matrix`` are non-negative.
        return incidences.row(0) in other.incidence_matrix().columns()

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

            sage: K.<cbrt2> = NumberField(x^3 - 2, 'a', embedding=1.26)         # optional - sage.rings.number_field
            sage: P = Polyhedron(vertices=[(1,1,cbrt2),(cbrt2,1,1)])            # optional - sage.rings.number_field
            sage: P.inequalities()                                              # optional - sage.rings.number_field
            (An inequality (-cbrt2^2 - cbrt2 - 1, 0, 0) x + cbrt2^2 + cbrt2 + 2 >= 0,
             An inequality (cbrt2^2 + cbrt2 + 1, 0, 0) x - cbrt2^2 + cbrt2 + 1 >= 0)
        """
        s = 'An inequality '
        have_A = not self.A().is_zero()
        if have_A:
            s += repr(self.A()) + ' x '
        if self.b() >= 0:
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
            [True, True, False, True, False, True, False, False]
        """
        try:
            if Vobj.is_vector():  # assume we were passed a point
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
            [True, True, False, True, False, True, False, False]

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

    def outer_normal(self):
        r"""
        Return the outer normal vector of ``self``.

        OUTPUT:

        The normal vector directed away from the interior of the polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = next(p.inequality_generator())
            sage: a.outer_normal()
            (1, -1, 0)
        """
        return -self.A()


class Equation(Hrepresentation):
    """
    A linear equation of the polyhedron. That is, the polyhedron is
    strictly smaller-dimensional than the ambient space, and contained
    in this hyperplane. Inherits from ``Hrepresentation``.
    """

    def type(self):
        r"""
        Return the type (equation/inequality/vertex/ray/line) as an
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

        .. NOTE::

            Return False for any equation.

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
        data = list(data)
        if len(data) != self._vector.degree():
            raise ValueError('V-representation data requires a list of length ambient_dim')

        self._vector[:] = data

        self._index = len(polyhedron._Vrepresentation)
        polyhedron._Vrepresentation.append(self)
        self._polyhedron = polyhedron
        if polyhedron.is_mutable():
            polyhedron._add_dependent_object(self)

    def is_V(self):
        """
        Return True if the object is part of a V-representation
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
        Return True if the object is a vertex of the V-representation.
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
        Return True if the object is a ray of the V-representation.
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
        Return True if the object is a line of the V-representation.
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
        Return a generator for the adjacent vertices/rays/lines.

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
        Return whether the incidence matrix element (self,Hobj) == 1

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: h1 = next(p.inequality_generator())
            sage: h1
            An inequality (-1, 0, 0) x + 1 >= 0
            sage: v1 = next(p.vertex_generator())
            sage: v1
            A vertex at (1, -1, -1)
            sage: v1.is_incident(h1)
            True
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
            0
        """
        return self.evaluated_on(Hobj)

    def incident(self):
        """
        Return a generator for the equations/inequalities that are satisfied on the given
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
        Return the type (equation/inequality/vertex/ray/line) as an
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
        Return a string representation of the vertex.

        OUTPUT:

        String.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: v = next(p.vertex_generator())
            sage: v.__repr__()
            'A vertex at (1, 0)'
        """
        return 'A vertex at ' + repr(self.vector())

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
        Return `A\vec{x}+b`

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: v = next(p.vertex_generator())
            sage: h = next(p.inequality_generator())
            sage: v
            A vertex at (1, -1, -1)
            sage: h
            An inequality (-1, 0, 0) x + 1 >= 0
            sage: v.evaluated_on(h)
            0
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
        Return the type (equation/inequality/vertex/ray/line) as an
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
        return 'A ray in the direction ' + repr(self.vector())

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
        Return `A\vec{r}`

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
        Return the type (equation/inequality/vertex/ray/line) as an
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
        return 'A line in the direction ' + repr(self.vector())

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
        Return `A\vec{\ell}`

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = next(p.line_generator())
            sage: h = next(p.inequality_generator())
            sage: a.evaluated_on(h)
            0
        """
        return Hobj.A() * self.vector()


def repr_pretty(coefficients, type, prefix='x', indices=None,
                latex=False, style='>=', split=False):
    r"""
    Return a pretty representation of equation/inequality represented
    by the coefficients.

    INPUT:

    - ``coefficients`` -- a tuple or other iterable

    - ``type`` -- either ``0`` (``PolyhedronRepresentation.INEQUALITY``)
      or ``1`` (``PolyhedronRepresentation.EQUATION``)

    - ``prefix`` -- a string

    - ``indices`` -- a tuple or other iterable

    - ``latex`` -- a boolean

    - ``split`` -- a boolean; (Default: ``False``). If set to ``True``,
                   the output is split into a 3-tuple containing the left-hand side,
                   the relation, and the right-hand side of the object.

    - ``style`` -- either ``"positive"`` (making all coefficients positive), or
                   ``"<="`` or ``">="``.

    OUTPUT:

    A string or 3-tuple of strings (depending on ``split``).

    EXAMPLES::

        sage: from sage.geometry.polyhedron.representation import repr_pretty
        sage: from sage.geometry.polyhedron.representation import PolyhedronRepresentation
        sage: print(repr_pretty((0, 1, 0, 0), PolyhedronRepresentation.INEQUALITY))
        x0 >= 0
        sage: print(repr_pretty((1, 2, 1, 0), PolyhedronRepresentation.INEQUALITY))
        2*x0 + x1 >= -1
        sage: print(repr_pretty((1, -1, -1, 1), PolyhedronRepresentation.EQUATION))
        -x0 - x1 + x2 == -1
    """
    from sage.misc.repr import repr_lincomb

    coeffs = list(coefficients)
    if indices is None:
        indices = range(len(coeffs)-1)
    vars = [1]
    if latex:
        vars += ['x_{{{}}}'.format(i) for i in indices]
    else:
        vars += ['x{}'.format(i) for i in indices]
    if type == PolyhedronRepresentation.EQUATION:
        rel = '=' if latex else '=='
    elif type == PolyhedronRepresentation.INEQUALITY:
        if style == '<=':
            rel = r'\leq' if latex else '<='
        else:
            rel = r'\geq' if latex else '>='
    else:
        raise NotImplementedError(
            'no pretty printing available: wrong type {}'.format(type))

    rvars = range(len(vars))

    if style == 'positive':
        pos_part = [max(c, 0) for c in coeffs]
        neg_part = [pos_part[i] - coeffs[i] for i in rvars]
        assert all(coeffs[i] == pos_part[i] - neg_part[i] for i in rvars)
        left_part = repr_lincomb([[vars[i], pos_part[i]] for i in rvars], is_latex=latex, strip_one=True)
        right_part = repr_lincomb([[vars[i], neg_part[i]] for i in rvars], is_latex=latex, strip_one=True)
    elif style == '>=':
        left_part = repr_lincomb([[vars[i], coeffs[i]] for i in rvars[1:]], is_latex=latex)
        right_part = repr_lincomb([[vars[0], -coeffs[0]]], is_latex=latex, strip_one=True)
    elif style == '<=':
        left_part = repr_lincomb([[vars[i], -coeffs[i]] for i in rvars[1:]], is_latex=latex)
        right_part = repr_lincomb([[vars[0], coeffs[0]]], is_latex=latex, strip_one=True)
    else:
        raise NotImplementedError('no pretty printing available: wrong style {}'.format(style))

    if not split:
        return '{} {} {}'.format(left_part, rel, right_part)
    else:
        return (str(left_part), rel, str(right_part))
