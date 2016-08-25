"""
The PPL (Parma Polyhedra Library) backend for polyhedral computations
"""
from __future__ import absolute_import

from sage.rings.all import ZZ, QQ
from sage.rings.integer import LCM_list
from sage.misc.functional import denominator
from sage.matrix.constructor import matrix
from sage.libs.ppl import (
    C_Polyhedron, Constraint_System, Generator_System,
    Variable, Linear_Expression,
    line, ray, point )

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ


#########################################################################
class Polyhedron_ppl(Polyhedron_base):
    """
    Polyhedra with ppl

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[], backend='ppl')
        sage: TestSuite(p).run()
    """

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

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Vrepresentation(p, [], [], [])
        """
        gs = Generator_System()
        if vertices is None: vertices = []
        for v in vertices:
            d = LCM_list([denominator(v_i) for v_i in v])
            if d.is_one():
                gs.insert(point(Linear_Expression(v, 0)))
            else:
                dv = [ d*v_i for v_i in v ]
                gs.insert(point(Linear_Expression(dv, 0), d))
        if rays is None: rays = []
        for r in rays:
            d = LCM_list([denominator(r_i) for r_i in r])
            if d.is_one():
                gs.insert(ray(Linear_Expression(r, 0)))
            else:
                dr = [ d*r_i for r_i in r ]
                gs.insert(ray(Linear_Expression(dr, 0)))
        if lines is None: lines = []
        for l in lines:
            d = LCM_list([denominator(l_i) for l_i in l])
            if d.is_one():
                gs.insert(line(Linear_Expression(l, 0)))
            else:
                dl = [ d*l_i for l_i in l ]
                gs.insert(line(Linear_Expression(dl, 0)))
        if gs.empty():
            self._ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'empty')
        else:
            self._ppl_polyhedron = C_Polyhedron(gs)
        self._init_Vrepresentation_from_ppl(minimize)
        self._init_Hrepresentation_from_ppl(minimize)

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

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Hrepresentation(p, [], [])
        """
        cs = Constraint_System()
        if ieqs is None: ieqs = []
        for ieq in ieqs:
            d = LCM_list([denominator(ieq_i) for ieq_i in ieq])
            dieq = [ ZZ(d*ieq_i) for ieq_i in ieq ]
            b = dieq[0]
            A = dieq[1:]
            cs.insert(Linear_Expression(A, b) >= 0)
        if eqns is None: eqns = []
        for eqn in eqns:
            d = LCM_list([denominator(eqn_i) for eqn_i in eqn])
            deqn = [ ZZ(d*eqn_i) for eqn_i in eqn ]
            b = deqn[0]
            A = deqn[1:]
            cs.insert(Linear_Expression(A, b) == 0)
        if cs.empty():
            self._ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'universe')
        else:
            self._ppl_polyhedron = C_Polyhedron(cs)
        self._init_Vrepresentation_from_ppl(minimize)
        self._init_Hrepresentation_from_ppl(minimize)

    def _init_Vrepresentation_from_ppl(self, minimize):
        """
        Create the Vrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ...                  backend='ppl')  # indirect doctest
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
        self._Vrepresentation = []
        gs = self._ppl_polyhedron.minimized_generators()
        parent = self.parent()
        for g in gs:
            if g.is_point():
                d = g.divisor()
                if d.is_one():
                    parent._make_Vertex(self, g.coefficients())
                else:
                    parent._make_Vertex(self, [x/d for x in g.coefficients()])
            elif g.is_ray():
                parent._make_Ray(self, g.coefficients())
            elif g.is_line():
                parent._make_Line(self, g.coefficients())
            else:
                assert False
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_from_ppl(self, minimize):
        """
        Create the Hrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ...                  backend='ppl')  # indirect doctest
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
        self._Hrepresentation = []
        cs = self._ppl_polyhedron.minimized_constraints()
        parent = self.parent()
        for c in cs:
            if c.is_inequality():
                parent._make_Inequality(self, (c.inhomogeneous_term(),) + c.coefficients())
            elif c.is_equality():
                parent._make_Equation(self, (c.inhomogeneous_term(),) + c.coefficients())
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron.

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
        super(Polyhedron_ppl, self)._init_empty_polyhedron()
        self._ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'empty')




#########################################################################
class Polyhedron_QQ_ppl(Polyhedron_ppl, Polyhedron_QQ):
    """
    Polyhedra over `\QQ` with ppl

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[],
        ...                  backend='ppl', base_ring=QQ)
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    pass


#########################################################################
class Polyhedron_ZZ_ppl(Polyhedron_ppl, Polyhedron_ZZ):
    """
    Polyhedra over `\ZZ` with ppl

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[])
        ...                  backend='ppl', base_ring=ZZ)
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    pass
