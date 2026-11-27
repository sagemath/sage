"""
The PPL (Parma Polyhedra Library) backend for polyhedral computations
"""

from sage.structure.element import Element
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.arith.functions import LCM_list
from sage.misc.functional import denominator
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from .base_mutable import Polyhedron_mutable
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ
from .representation import VERTEX, RAY, LINE, INEQUALITY, EQUATION

from sage.misc.lazy_import import lazy_import
from sage.features import PythonModule
lazy_import('ppl', ['C_Polyhedron', 'Generator_System', 'Constraint_System',
                    'Linear_Expression', 'line', 'ray', 'point'],
                    feature=PythonModule("ppl", spkg='pplpy', type='standard'))


#########################################################################
class Polyhedron_ppl(Polyhedron_mutable):
    """
    Polyhedra with ppl.

    INPUT:

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[], backend='ppl')
        sage: TestSuite(p).run()
    """

    _backend_object_name = "ppl_polyhedron"
    _is_mutable = True

    def __init__(self, parent, Vrep, Hrep, ppl_polyhedron=None, mutable=False, **kwds):
        """
        Initialize the polyhedron.

        See :class:`Polyhedron_ppl` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron()
            sage: TestSuite(p).run()
            sage: p = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)])
            sage: TestSuite(p).run()
            sage: q = polytopes.cube()
            sage: p = q.parent().element_class(q.parent(), None, None, q._ppl_polyhedron)
            sage: TestSuite(p).run()
        """
        # This is important. For some reason the element constructor copies the list sometimes.
        self._dependent_objects = []
        if ppl_polyhedron:
            if Hrep is not None or Vrep is not None:
                raise ValueError("only one of Vrep, Hrep, or ppl_polyhedron can be different from None")
            Element.__init__(self, parent=parent)
            minimize = bool('minimize' in kwds and kwds['minimize'])
            self._init_from_ppl_polyhedron(ppl_polyhedron, minimize)
        else:
            Polyhedron_mutable.__init__(self, parent, Vrep, Hrep, **kwds)
        if not mutable:
            self.set_immutable()

    def _init_from_Vrepresentation(self, vertices, rays, lines, minimize=True, verbose=False):
        """
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of point. Each point can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``rays`` -- list of rays. Each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``lines`` -- list of lines. Each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Vrepresentation(p, [], [], [])
        """
        gs = self._convert_generators_to_ppl(vertices, rays, lines)
        if gs.empty():
            ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'empty')
        else:
            ppl_polyhedron = C_Polyhedron(gs)
        self._init_from_ppl_polyhedron(ppl_polyhedron, minimize)

    def _init_from_Hrepresentation(self, ieqs, eqns, minimize=True, verbose=False):
        """
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Hrepresentation(p, [], [])
        """
        cs = self._convert_constraints_to_ppl(ieqs, eqns)
        if cs.empty():
            ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'universe')
        else:
            ppl_polyhedron = C_Polyhedron(cs)
        self._init_from_ppl_polyhedron(ppl_polyhedron, minimize)

    def _init_from_ppl_polyhedron(self, ppl_polyhedron, minimize=True):
        """
        Create the V-/Hrepresentation objects from the ppl polyhedron.

        TESTS::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Hrepresentation(p, [], [])  # indirect doctest
        """
        self._ppl_polyhedron = ppl_polyhedron

    def set_immutable(self):
        r"""
        Make this polyhedron immutable. This operation cannot be undone.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: p.is_mutable()
            True
            sage: hasattr(p, "_Vrepresentation")
            False
            sage: p.set_immutable()
            sage: hasattr(p, "_Vrepresentation")
            True

        TESTS:

        Check that :issue:`33666` is fixed::

            sage: cube = polytopes.cube()
            sage: parent = cube.parent()
            sage: smaller_cube_ZZ = parent._element_constructor_(1/2 * cube, mutable=True)
            sage: smaller_cube_ZZ.set_immutable()
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: smaller_cube_ZZ.is_immutable()
            False
            sage: smaller_cube_ZZ.set_immutable()
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: smaller_cube_ZZ.is_immutable()
            False
            sage: smaller_cube_QQ = smaller_cube_ZZ.base_extend(QQ)
            sage: smaller_cube_QQ.set_immutable()
            sage: smaller_cube_QQ.is_immutable()
            True
        """
        if not hasattr(self, '_Vrepresentation'):
            try:
                self._init_Vrepresentation_from_ppl(True)
            except TypeError as e:
                # Apparently the polyhedron is (no longer) integral.
                self._clear_cache()
                raise e
        if not hasattr(self, '_Hrepresentation'):
            self._init_Hrepresentation_from_ppl(True)
        self._is_mutable = False

    def Vrepresentation(self, index=None):
        """
        Return the objects of the V-representation. Each entry is
        either a vertex, a ray, or a line.

        See :mod:`sage.geometry.polyhedron.constructor` for a
        definition of vertex/ray/line.

        INPUT:

        - ``index`` -- either an integer or ``None``

        OUTPUT:

        The optional argument is an index running from ``0`` to
        ``self.n_Vrepresentation()-1``. If present, the
        V-representation object at the given index will be
        returned. Without an argument, returns the list of all
        V-representation objects.

        EXAMPLES::

            sage: p = polytopes.cube()
            sage: p.Vrepresentation(0)
            A vertex at (1, -1, -1)

        ::

            sage: P = p.parent()
            sage: p = P._element_constructor_(p, mutable=True)
            sage: p.Vrepresentation(0)
            A vertex at (-1, -1, -1)
            sage: p._clear_cache()
            sage: p.Vrepresentation(0)
            A vertex at (-1, -1, -1)
            sage: TestSuite(p).run()

        TESTS:

        Check that :issue:`33666` is fixed::

            sage: cube = polytopes.cube()
            sage: parent = cube.parent()
            sage: smaller_cube_ZZ = parent._element_constructor_(1/2 * cube, mutable=True)
            sage: smaller_cube_ZZ.Hrepresentation()
            (An inequality (0, 0, -2) x + 1 >= 0,
            An inequality (0, -2, 0) x + 1 >= 0,
            An inequality (-2, 0, 0) x + 1 >= 0,
            An inequality (2, 0, 0) x + 1 >= 0,
            An inequality (0, 0, 2) x + 1 >= 0,
            An inequality (0, 2, 0) x + 1 >= 0)
            sage: smaller_cube_ZZ.Vrepresentation()
            Traceback (most recent call last):
            ...
            TypeError: the polyhedron is not integral; do a base extension ``self.base_extend(QQ)``
            sage: smaller_cube_ZZ.Vrepresentation()
            Traceback (most recent call last):
            ...
            TypeError: the polyhedron is not integral; do a base extension ``self.base_extend(QQ)``
            sage: smaller_cube_QQ = smaller_cube_ZZ.base_extend(QQ)
            sage: smaller_cube_QQ.Hrepresentation()
            (An inequality (0, 0, -2) x + 1 >= 0,
            An inequality (0, -2, 0) x + 1 >= 0,
            An inequality (-2, 0, 0) x + 1 >= 0,
            An inequality (2, 0, 0) x + 1 >= 0,
            An inequality (0, 0, 2) x + 1 >= 0,
            An inequality (0, 2, 0) x + 1 >= 0)
        """
        if not hasattr(self, '_Vrepresentation'):
            try:
                self._init_Vrepresentation_from_ppl(True)
            except TypeError:
                # Apparently the polyhedron is (no longer) integral.
                self._clear_cache()
                raise TypeError("the polyhedron is not integral; do a base extension ``self.base_extend(QQ)``")
        if index is None:
            return self._Vrepresentation
        return self._Vrepresentation[index]

    def _init_Vrepresentation_from_ppl(self, minimize):
        """
        Create the Vrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest
            ....:                backend='ppl')
            sage: p.Hrepresentation()
            (An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0)
            sage: p._ppl_polyhedron.minimized_constraints()
            Constraint_System {x0+4*x1-2>=0, x0-12*x1+6>=0, -5*x0+12*x1+10>=0}
            sage: p.Vrepresentation()
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
            sage: p._ppl_polyhedron.minimized_generators()
            Generator_System {point(0/2, 1/2), point(2/1, 0/1), point(24/6, 5/6)}

        PPL may choose a non-integral point generator even when its
        class modulo the lineality has an integral representative.  In
        the ``ZZ`` backend, we replace it by such a representative when
        one exists (:issue:`42142`)::

            sage: P = Polyhedron(eqns=[(-1, 1, 1, 1, -2)], base_ring=ZZ)
            sage: P.Vrepresentation()
            (A line in the direction (0, 0, 2, 1),
             A line in the direction (2, 0, 0, 1),
             A line in the direction (0, 2, 0, 1),
             A vertex at (1, 0, 0, 0))

        Point generators are treated independently.  Here PPL produces
        two non-integral point generators, both of which have integral
        representatives modulo the lineality::

            sage: P = Polyhedron(ieqs=[(-1, 1, -2), (3, -1, 2)],
            ....:                base_ring=ZZ, backend='ppl')
            sage: P._ppl_polyhedron.minimized_generators()
            Generator_System {line(2, 1), point(0/2, -3/2), point(0/2, -1/2)}
            sage: P.n_vertices()
            2
            sage: all(v.vector() in ZZ^2 for v in P.vertices())
            True

        If no integral representative exists, construction fails::

            sage: Polyhedron(eqns=[(-1, 2)], base_ring=ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: Polyhedron(eqns=[(-1, 2, 2)], base_ring=ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        Merely containing an integral point is not sufficient.  Every
        point-generator class modulo the lineality must have an integral
        representative.  The following strip contains ``(1, 0)`` in its
        interior, but neither boundary class has an integral representative::

            sage: P = Polyhedron(ieqs=[(-1, 2, 0), (3, -2, 0)],
            ....:                base_ring=QQ, backend='ppl')
            sage: P.interior_contains((1, 0))
            True
            sage: Polyhedron(ieqs=[(-1, 2, 0), (3, -2, 0)],
            ....:            base_ring=ZZ, backend='ppl')
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
        """
        if not self._is_mutable:
            raise TypeError("Vrepresentation of mutable polyhedra cannot be recomputed")
        gs = self._ppl_polyhedron.minimized_generators()
        gs = tuple(gs)
        parent = self.parent()
        lines = [[Integer(mpz) for mpz in g.coefficients()]
                 for g in gs if g.is_line()]
        self._Vrepresentation = []
        for g in gs:
            coefficients = [Integer(mpz) for mpz in g.coefficients()]
            if g.is_point():
                d = Integer(g.divisor())
                if d.is_one():
                    parent._make_Vertex(self, coefficients)
                else:
                    point = [x/d for x in coefficients]
                    if parent.base_ring() is ZZ and lines:
                        point = self._integral_representative_mod_lines(point, lines)
                    parent._make_Vertex(self, point)
            elif g.is_ray():
                parent._make_Ray(self, coefficients)
            elif g.is_line():
                parent._make_Line(self, coefficients)
            else:
                assert False
        self._Vrepresentation = tuple(self._Vrepresentation)

    @staticmethod
    def _integral_representative_mod_lines(point, lines):
        r"""
        Return an integral point equivalent to ``point`` modulo ``lines``.

        A minimized PPL generator system represents a polyhedron as
        ``conv(points) + cone(rays) + span(lines)``.  Replacing any point
        generator ``p`` independently by a point in ``p + span(lines)``
        leaves the represented polyhedron unchanged.

        Modulo the lineality, the point generators are the vertices of the
        quotient polyhedron.  Consequently, an integral V-representation
        exists only if every one of these classes has an integral
        representative; an integral point elsewhere in the polyhedron is
        not sufficient.
        """
        line_matrix = matrix(ZZ, lines)
        equations = line_matrix.right_kernel_matrix()
        if not equations.nrows():
            return [ZZ.zero()] * len(point)

        point = vector(point)
        rhs = equations.change_ring(point.base_ring()) * point
        if any(denominator(x) != 1 for x in rhs):
            return point

        hermite_transpose, transformation_transpose = equations.transpose().hermite_form(
            transformation=True)
        hermite = hermite_transpose.transpose()
        transformation = transformation_transpose.transpose()
        rhs = vector(ZZ, rhs)
        coordinates = vector(ZZ, equations.ncols())
        for i in range(hermite.nrows()):
            residual = rhs[i] - sum(hermite[i, j] * coordinates[j] for j in range(i))
            if not residual:
                continue
            d = hermite[i, i]
            if not d or residual % d:
                return point
            coordinates[i] = residual // d
        return transformation * coordinates

    def _init_Hrepresentation_from_ppl(self, minimize):
        """
        Create the Hrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest
            ....:                backend='ppl')
            sage: p.Hrepresentation()
            (An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0)
            sage: p._ppl_polyhedron.minimized_constraints()
            Constraint_System {x0+4*x1-2>=0, x0-12*x1+6>=0, -5*x0+12*x1+10>=0}
            sage: p.Vrepresentation()
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
            sage: p._ppl_polyhedron.minimized_generators()
            Generator_System {point(0/2, 1/2), point(2/1, 0/1), point(24/6, 5/6)}
        """
        if not self._is_mutable:
            raise TypeError("Hrepresentation of mutable polyhedra cannot be recomputed")
        self._Hrepresentation = []
        cs = self._ppl_polyhedron.minimized_constraints()
        parent = self.parent()
        for c in cs:
            if c.is_inequality():
                parent._make_Inequality(self, (c.inhomogeneous_term(),) + c.coefficients())
            elif c.is_equality():
                parent._make_Equation(self, (c.inhomogeneous_term(),) + c.coefficients())
        self._Hrepresentation = tuple(self._Hrepresentation)

    def Hrepresentation(self, index=None):
        """
        Return the objects of the H-representation. Each entry is
        either an inequality or a equation.

        INPUT:

        - ``index`` -- either an integer or ``None``

        OUTPUT:

        The optional argument is an index running from ``0`` to
        ``self.n_Hrepresentation()-1``. If present, the
        H-representation object at the given index will be
        returned. Without an argument, returns the list of all
        H-representation objects.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p.Hrepresentation(0)
            An inequality (-1, 0, 0) x + 1 >= 0
            sage: p.Hrepresentation(0) == p.Hrepresentation()[0]
            True

        ::

            sage: P = p.parent()
            sage: p = P._element_constructor_(p, mutable=True)
            sage: p.Hrepresentation(0)
            An inequality (0, 0, -1) x + 1 >= 0
            sage: p._clear_cache()
            sage: p.Hrepresentation(0)
            An inequality (0, 0, -1) x + 1 >= 0
            sage: TestSuite(p).run()
        """
        if not hasattr(self, '_Hrepresentation'):
            self._init_Hrepresentation_from_ppl(True)
        if index is None:
            return self._Hrepresentation
        return self._Hrepresentation[index]

    def _init_empty_polyhedron(self):
        """
        Initialize an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='ppl'); empty
            The empty polyhedron in ZZ^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [], backend='ppl')
            The empty polyhedron in ZZ^0
            sage: Polyhedron(backend='ppl')._init_empty_polyhedron()
        """
        super()._init_empty_polyhedron()
        self._ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'empty')

    @staticmethod
    def _convert_generator_to_ppl(v, typ):
        r"""
        Convert a generator to ``ppl``.

        INPUT:

        - ``v`` -- a vertex, ray, or line

        - ``typ`` -- integer according to
          `:sage:`~sage.geometry.polyhedron.representation.LINE` etc.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.representation import VERTEX, RAY, LINE
            sage: P = Polyhedron()
            sage: P._convert_generator_to_ppl([1, 1/2, 3], VERTEX)
            point(2/2, 1/2, 6/2)
            sage: P._convert_generator_to_ppl([1, 1/2, 3], RAY)
            ray(2, 1, 6)
            sage: P._convert_generator_to_ppl([1, 1/2, 3], LINE)
            line(2, 1, 6)
        """
        if typ == VERTEX:
            ob = point
        elif typ == RAY:
            ob = ray
        else:
            ob = line

        d = LCM_list([denominator(v_i) for v_i in v])
        if d.is_one():
            return ob(Linear_Expression(v, 0))
        dv = [ d*v_i for v_i in v ]
        if typ == VERTEX:
            return ob(Linear_Expression(dv, 0), d)
        return ob(Linear_Expression(dv, 0))

    @staticmethod
    def _convert_generators_to_ppl(vertices, rays, lines):
        r"""
        Convert generators to a ``ppl`` generator system.

        INPUT:

        - ``vertices`` -- iterable of vertices or ``None``

        - ``rays`` -- iterable of rays or ``None``

        - ``lines`` -- iterable of lines or ``None``

        EXAMPLES::

            sage: P = Polyhedron()
            sage: P._convert_generators_to_ppl([[1, 1/2, 3]], [[0, 1, 3/2]], [[0, 0, 1]])
            Generator_System {point(2/2, 1/2, 6/2), ray(0, 2, 3), line(0, 0, 1)}
        """
        gs = Generator_System()
        if vertices is None:
            vertices = []
        for v in vertices:
            gs.insert(Polyhedron_ppl._convert_generator_to_ppl(v, VERTEX))
        if rays is None:
            rays = []
        for r in rays:
            gs.insert(Polyhedron_ppl._convert_generator_to_ppl(r, RAY))
        if lines is None:
            lines = []
        for l in lines:
            gs.insert(Polyhedron_ppl._convert_generator_to_ppl(l, LINE))
        return gs

    @staticmethod
    def _convert_constraint_to_ppl(c, typ):
        r"""
        Convert a constraint to ``ppl``.

        INPUT:

        - ``c`` -- an inequality or equation

        - ``typ`` -- integer according to
          `:sage:`~sage.geometry.polyhedron.representation.INEQUALITY` etc.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.representation import INEQUALITY, EQUATION
            sage: P = Polyhedron()
            sage: P._convert_constraint_to_ppl([1, 1/2, 3], INEQUALITY)
            x0+6*x1+2>=0
            sage: P._convert_constraint_to_ppl([1, 1/2, 3], EQUATION)
            x0+6*x1+2==0
        """
        d = LCM_list([denominator(c_i) for c_i in c])
        dc = [ ZZ(d*c_i) for c_i in c ]
        b = dc[0]
        A = dc[1:]
        if typ == INEQUALITY:
            return Linear_Expression(A, b) >= 0
        return Linear_Expression(A, b) == 0

    @staticmethod
    def _convert_constraints_to_ppl(ieqs, eqns):
        r"""
        Convert constraints to a ``ppl`` constraint system.

        INPUT:

        - ``ieqs`` -- iterable of inequalities or ``None``

        - ``eqns`` -- iterable of equations or ``None``

        EXAMPLES::

            sage: P = Polyhedron()
            sage: P._convert_constraints_to_ppl([[1, 1/2, 3]], None)
            Constraint_System {x0+6*x1+2>=0}
        """
        cs = Constraint_System()
        if ieqs is None:
            ieqs = []
        for ieq in ieqs:
            cs.insert(Polyhedron_ppl._convert_constraint_to_ppl(ieq, INEQUALITY))
        if eqns is None:
            eqns = []
        for eqn in eqns:
            cs.insert(Polyhedron_ppl._convert_constraint_to_ppl(eqn, EQUATION))
        return cs


#########################################################################
class Polyhedron_QQ_ppl(Polyhedron_ppl, Polyhedron_QQ):
    r"""
    Polyhedra over `\QQ` with ppl.

    INPUT:

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[],
        ....:                backend='ppl', base_ring=QQ)
        sage: TestSuite(p).run()
    """
    pass


#########################################################################
class Polyhedron_ZZ_ppl(Polyhedron_ppl, Polyhedron_ZZ):
    r"""
    Polyhedra over `\ZZ` with ppl.

    INPUT:

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[],
        ....:                backend='ppl', base_ring=ZZ)
        sage: TestSuite(p).run()
    """
    pass
