r"""
Base class for polyhedra
"""

# ****************************************************************************
#       Copyright (C) 2008-2012 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011-2015 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2012-2018 Frederic Chapoton
#       Copyright (C) 2013      Andrey Novoseltsev
#       Copyright (C) 2014-2017 Moritz Firsching
#       Copyright (C) 2014-2019 Thierry Monteil
#       Copyright (C) 2015      Nathann Cohen
#       Copyright (C) 2015-2017 Jeroen Demeyer
#       Copyright (C) 2015-2017 Vincent Delecroix
#       Copyright (C) 2015-2018 Dima Pasechnik
#       Copyright (C) 2015-2020 Jean-Philippe Labbe <labbe at math.huji.ac.il>
#       Copyright (C) 2015-2021 Matthias Koeppe
#       Copyright (C) 2016-2019 Daniel Krenn
#       Copyright (C) 2017      Marcelo Forets
#       Copyright (C) 2017-2018 Mark Bell
#       Copyright (C) 2019      Julian Ritter
#       Copyright (C) 2019-2020 Laith Rastanawi
#       Copyright (C) 2019-2020 Sophia Elia
#       Copyright (C) 2019-2021 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools

from sage.structure.element import Element, coerce_binop, is_Vector, is_Matrix
from sage.structure.richcmp import rich_to_bool, op_NE
from sage.cpython.string import bytes_to_str

from sage.misc.cachefunc import cached_method
from sage.misc.all import prod
from sage.misc.randstate import current_randstate
from sage.misc.superseded import deprecated_function_alias

from sage.rings.all import QQ, ZZ, AA
from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.modules.vector_space_morphism import linear_transformation
from sage.matrix.constructor import matrix
from sage.functions.other import sqrt, floor, ceil
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.graphs.graph import Graph
from sage.geometry.convex_set import ConvexSet_closed, AffineHullProjectionData

from .constructor import Polyhedron
from sage.geometry.relative_interior import RelativeInterior
from sage.categories.sets_cat import EmptySetError

#########################################################################
# Notes if you want to implement your own backend:
#
#  * derive from Polyhedron_base
#
#  * you must implement _init_from_Vrepresentation and
#    _init_from_Hrepresentation
#
#  * You might want to override _init_empty_polyhedron
#
#  * You may implement _init_from_Vrepresentation_and_Hrepresentation
#
#  * You can of course also override any other method for which you
#    have a faster implementation.
#########################################################################


#########################################################################
def is_Polyhedron(X):
    """
    Test whether ``X`` is a Polyhedron.

    INPUT:

    - ``X`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = polytopes.hypercube(2)
        sage: from sage.geometry.polyhedron.base import is_Polyhedron
        sage: is_Polyhedron(p)
        True
        sage: is_Polyhedron(123456)
        False
    """
    return isinstance(X, Polyhedron_base)


#########################################################################
class Polyhedron_base(Element, ConvexSet_closed):
    """
    Base class for Polyhedron objects

    INPUT:

    - ``parent`` -- the parent, an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedra`.

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``. The
      V-representation of the polyhedron. If ``None``, the polyhedron
      is determined by the H-representation.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``. The
      H-representation of the polyhedron. If ``None``, the polyhedron
      is determined by the V-representation.

    - ``Vrep_minimal`` (optional) -- see below

    - ``Hrep_minimal`` (optional) -- see below

    - ``pref_rep`` -- string (default: ``None``);
       one of``Vrep`` or ``Hrep`` to pick this in case the backend
       cannot initialize from complete double description

    - ``mutable`` -- ignored

    If both ``Vrep`` and ``Hrep`` are provided, then
    ``Vrep_minimal`` and ``Hrep_minimal`` must be set to ``True``.

    TESTS::

        sage: p = Polyhedron()
        sage: TestSuite(p).run()

    ::

        sage: p = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)], base_ring=ZZ)
        sage: TestSuite(p).run()

    ::

        sage: p=polytopes.flow_polytope(digraphs.DeBruijn(3,2))
        sage: TestSuite(p).run()

    ::

        sage: TestSuite(Polyhedron([[]])).run()
        sage: TestSuite(Polyhedron([[0]])).run()
        sage: TestSuite(Polyhedron([[1]])).run()

    ::

        sage: P = polytopes.permutahedron(3) * Polyhedron(rays=[[0,0,1],[0,1,1],[1,2,3]])
        sage: TestSuite(P).run()

    ::

        sage: P = polytopes.permutahedron(3)*Polyhedron(rays=[[0,0,1],[0,1,1]], lines=[[1,0,0]])
        sage: TestSuite(P).run()

    ::

        sage: M = random_matrix(ZZ, 5, 5, distribution='uniform')
        sage: while True:
        ....:     M = random_matrix(ZZ, 5, 5, distribution='uniform')
        ....:     if M.rank() != 5:
        ....:         break
        ....:
        sage: P = Polyhedron(M)
        sage: TestSuite(P).run()
    """

    def __init__(self, parent, Vrep, Hrep, Vrep_minimal=None, Hrep_minimal=None, pref_rep=None, mutable=False, **kwds):
        """
        Initializes the polyhedron.

        See :class:`Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron()    # indirect doctests

            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: from sage.geometry.polyhedron.parent import Polyhedra_field
            sage: parent = Polyhedra_field(AA, 1, 'field')
            sage: Vrep = [[[0], [1/2], [1]], [], []]
            sage: Hrep = [[[0, 1], [1, -1]], []]
            sage: p = Polyhedron_field(parent, Vrep, Hrep,
            ....:                      Vrep_minimal=False, Hrep_minimal=True)
            Traceback (most recent call last):
            ...
            ValueError: if both Vrep and Hrep are provided, they must be minimal...

        Illustration of ``pref_rep``.
        Note that ``ppl`` doesn't support precomputed data::

            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_QQ_ppl
            sage: from sage.geometry.polyhedron.parent import Polyhedra_QQ_ppl
            sage: parent = Polyhedra_QQ_ppl(QQ, 1, 'ppl')
            sage: p = Polyhedron_QQ_ppl(parent, Vrep, 'nonsense',
            ....:                       Vrep_minimal=True, Hrep_minimal=True, pref_rep='Vrep')
            sage: p = Polyhedron_QQ_ppl(parent, 'nonsense', Hrep,
            ....:                       Vrep_minimal=True, Hrep_minimal=True, pref_rep='Hrep')
            sage: p = Polyhedron_QQ_ppl(parent, 'nonsense', Hrep,
            ....:                       Vrep_minimal=True, Hrep_minimal=True, pref_rep='Vrepresentation')
            Traceback (most recent call last):
            ...
            ValueError: ``pref_rep`` must be one of ``(None, 'Vrep', 'Hrep')``

        If the backend supports precomputed data, ``pref_rep`` is ignored::

            sage: p = Polyhedron_field(parent, Vrep, 'nonsense',  # py3
            ....:                      Vrep_minimal=True, Hrep_minimal=True, pref_rep='Vrep')
            Traceback (most recent call last):
            ...
            TypeError: _init_Hrepresentation() takes 3 positional arguments but 9 were given
            sage: p = Polyhedron_field(parent, Vrep, 'nonsense',  # py2
            ....:                      Vrep_minimal=True, Hrep_minimal=True, pref_rep='Vrep')
            Traceback (most recent call last):
            ...
            TypeError: _init_Hrepresentation() takes exactly 3 arguments (9 given)

        The empty polyhedron is detected when the Vrepresentation is given with generator;
        see :trac:`29899`::

            sage: from sage.geometry.polyhedron.backend_cdd import Polyhedron_QQ_cdd
            sage: from sage.geometry.polyhedron.parent import Polyhedra_QQ_cdd
            sage: parent = Polyhedra_QQ_cdd(QQ, 0, 'cdd')
            sage: p = Polyhedron_QQ_cdd(parent, [iter([]), iter([]), iter([])], None)
        """
        Element.__init__(self, parent=parent)
        if Vrep is not None and Hrep is not None:
            if not (Vrep_minimal is True and Hrep_minimal is True):
                raise ValueError("if both Vrep and Hrep are provided, they must be minimal"
                                 " and Vrep_minimal and Hrep_minimal must both be True")
            if hasattr(self, "_init_from_Vrepresentation_and_Hrepresentation"):
                self._init_from_Vrepresentation_and_Hrepresentation(Vrep, Hrep)
                return
            else:
                if pref_rep is None:
                    # Initialize from Hrepresentation if this seems simpler.
                    Vrep = [tuple(Vrep[0]), tuple(Vrep[1]), Vrep[2]]
                    Hrep = [tuple(Hrep[0]), Hrep[1]]
                    if len(Hrep[0]) < len(Vrep[0]) + len(Vrep[1]):
                        pref_rep = 'Hrep'
                    else:
                        pref_rep = 'Vrep'
                if pref_rep == 'Vrep':
                    Hrep = None
                elif pref_rep == 'Hrep':
                    Vrep = None
                else:
                    raise ValueError("``pref_rep`` must be one of ``(None, 'Vrep', 'Hrep')``")
        if Vrep is not None:
            vertices, rays, lines = Vrep

            # We build tuples out of generators now to detect the empty polyhedron.

            # The damage is limited:
            # The backend will have to obtain all elements from the generator anyway.
            # The generators are mainly for saving time with initializing from
            # Vrepresentation and Hrepresentation.
            # If we dispose of one of them (see above), it is wasteful to have generated it.

            # E.g. the dilate will be set up with new Vrepresentation and Hrepresentation
            # regardless of the backend along with the argument ``pref_rep``.
            # As we only use generators, there is no penalty to this approach
            # (and the method ``dilation`` does not have to distinguish by backend).

            if not isinstance(vertices, (tuple, list)):
                vertices = tuple(vertices)
            if not isinstance(rays, (tuple, list)):
                rays = tuple(rays)
            if not isinstance(lines, (tuple, list)):
                lines = tuple(lines)

            if vertices or rays or lines:
                self._init_from_Vrepresentation(vertices, rays, lines, **kwds)
            else:
                self._init_empty_polyhedron()
        elif Hrep is not None:
            ieqs, eqns = Hrep
            self._init_from_Hrepresentation(ieqs, eqns, **kwds)
        else:
            self._init_empty_polyhedron()

    def __hash__(self):
        r"""
        TESTS::

            sage: K.<a> = QuadraticField(2)
            sage: p = Polyhedron(vertices=[(0,1,a),(3,a,5)],
            ....:                rays=[(a,2,3), (0,0,1)],
            ....:                base_ring=K)
            sage: q = Polyhedron(vertices=[(3,a,5),(0,1,a)],
            ....:                rays=[(0,0,1), (a,2,3)],
            ....:                base_ring=K)
            sage: hash(p) == hash(q)
            True
        """
        # TODO: find something better *but* fast
        return hash((self.dim(),
                     self.ambient_dim(),
                     self.n_Hrepresentation(),
                     self.n_Vrepresentation(),
                     self.n_equations(),
                     self.n_facets(),
                     self.n_inequalities(),
                     self.n_lines(),
                     self.n_rays(),
                     self.n_vertices()))

    def _sage_input_(self, sib, coerced):
        """
        Return Sage command to reconstruct ``self``.

        See :mod:`sage.misc.sage_input` for details.

        .. TODO::

            Add the option ``preparse`` to the method.

        EXAMPLES::

            sage: P = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]], backend='ppl')
            sage: sage_input(P)
            Polyhedron(backend='ppl', base_ring=QQ, rays=[(QQ(1), QQ(1))], vertices=[(QQ(0), QQ(1)), (QQ(1), QQ(0))])
            sage: P = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]], backend='normaliz') # optional - pynormaliz
            sage: sage_input(P)                                                                    # optional - pynormaliz
            Polyhedron(backend='normaliz', base_ring=QQ, rays=[(QQ(1), QQ(1))], vertices=[(QQ(0), QQ(1)), (QQ(1), QQ(0))])
            sage: P = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]], backend='polymake') # optional - polymake
            sage: sage_input(P)                                                                    # optional - polymake
            Polyhedron(backend='polymake', base_ring=QQ, rays=[(QQ(1), QQ(1))], vertices=[(QQ(1), QQ(0)), (QQ(0), QQ(1))])
       """
        kwds = dict()
        kwds['base_ring'] = sib(self.base_ring())
        kwds['backend'] = sib(self.backend())
        if self.n_vertices() > 0:
            kwds['vertices'] = [sib(tuple(v)) for v in self.vertices()]
        if self.n_rays() > 0:
            kwds['rays'] = [sib(tuple(r)) for r in self.rays()]
        if self.n_lines() > 0:
            kwds['lines'] = [sib(tuple(l)) for l in self.lines()]
        return sib.name('Polyhedron')(**kwds)

    def _init_from_Vrepresentation(self, vertices, rays, lines, **kwds):
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

        EXAMPLES::

            sage: p = Polyhedron()
            sage: from sage.geometry.polyhedron.base import Polyhedron_base
            sage: Polyhedron_base._init_from_Vrepresentation(p, [], [], [])
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this method
        """
        raise NotImplementedError('a derived class must implement this method')

    def _init_from_Hrepresentation(self, ieqs, eqns, **kwds):
        """
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        EXAMPLES::

            sage: p = Polyhedron()
            sage: from sage.geometry.polyhedron.base import Polyhedron_base
            sage: Polyhedron_base._init_from_Hrepresentation(p, [], [])
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this method
        """
        raise NotImplementedError('a derived class must implement this method')

    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [])
            The empty polyhedron in ZZ^0
            sage: Polyhedron(vertices = [])._init_empty_polyhedron()
            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: Polyhedra(QQ,7)()
            A 0-dimensional polyhedron in QQ^7 defined as the convex hull of 1 vertex
        """
        self._Vrepresentation = []
        self._Hrepresentation = []
        self.parent()._make_Equation(self, [-1] + [0] * self.ambient_dim())
        self._Vrepresentation = tuple(self._Vrepresentation)
        self._Hrepresentation = tuple(self._Hrepresentation)

        V_matrix = matrix(ZZ, 0, 0, 0)
        V_matrix.set_immutable()
        self.vertex_adjacency_matrix.set_cache(V_matrix)

        H_matrix = matrix(ZZ, 1, 1, 0)
        H_matrix.set_immutable()
        self.facet_adjacency_matrix.set_cache(H_matrix)

    def _facet_adjacency_matrix(self):
        """
        Compute the facet adjacency matrix in case it has not been
        computed during initialization.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)])
            sage: p._facet_adjacency_matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]

        Checks that :trac:`22455` is fixed::

            sage: s = polytopes.simplex(2)
            sage: s._facet_adjacency_matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]

        """
        # TODO: This implementation computes the whole face lattice,
        # which is much more information than necessary.
        M = matrix(ZZ, self.n_facets(), self.n_facets(), 0)
        codim = self.ambient_dim()-self.dim()

        def set_adjacent(h1, h2):
            if h1 is h2:
                return
            i = h1.index() - codim
            j = h2.index() - codim
            M[i, j] = 1
            M[j, i] = 1

        for face in self.faces(self.dim()-2):
            Hrep = face.ambient_Hrepresentation()
            assert(len(Hrep) == codim+2)
            set_adjacent(Hrep[-2], Hrep[-1])
        M.set_immutable()
        return M

    def _vertex_adjacency_matrix(self):
        """
        Compute the vertex adjacency matrix in case it has not been
        computed during initialization.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)])
            sage: p._vertex_adjacency_matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        # TODO: This implementation computes the whole face lattice,
        # which is much more information than necessary.
        M = matrix(ZZ, self.n_Vrepresentation(), self.n_Vrepresentation(), 0)

        def set_adjacent(v1, v2):
            if v1 is v2:
                return
            i = v1.index()
            j = v2.index()
            M[i, j] = 1
            M[j, i] = 1

        face_lattice = self.face_lattice()
        for face in face_lattice:
            Vrep = face.ambient_Vrepresentation()
            if len(Vrep) == 2:
                set_adjacent(Vrep[0], Vrep[1])
        M.set_immutable()
        return M

    def _delete(self):
        """
        Delete this polyhedron.

        This speeds up creation of new polyhedra by reusing
        objects. After recycling a polyhedron object, it is not in a
        consistent state any more and neither the polyhedron nor its
        H/V-representation objects may be used any more.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.parent.Polyhedra_base.recycle`

        EXAMPLES::

            sage: p = Polyhedron([(0,0),(1,0),(0,1)])
            sage: p._delete()

            sage: vertices = [(0,0,0,0),(1,0,0,0),(0,1,0,0),(1,1,0,0),(0,0,1,0),(0,0,0,1)]
            sage: def loop_polyhedra():
            ....:     for i in range(100):
            ....:         p = Polyhedron(vertices)

            sage: timeit('loop_polyhedra()')                   # not tested - random
            5 loops, best of 3: 79.5 ms per loop

            sage: def loop_polyhedra_with_recycling():
            ....:     for i in range(100):
            ....:         p = Polyhedron(vertices)
            ....:         p._delete()

            sage: timeit('loop_polyhedra_with_recycling()')    # not tested - random
            5 loops, best of 3: 57.3 ms per loop
        """
        self.parent().recycle(self)

    def _test_basic_properties(self, tester=None, **options):
        """
        Run some basic tests to see, that some general assertion on polyhedra hold.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_basic_properties()
        """
        if tester is None:
            tester = self._tester(**options)

        tester.assertEqual(self.n_vertices() + self.n_rays() + self.n_lines(), self.n_Vrepresentation())
        tester.assertEqual(self.n_inequalities() + self.n_equations(), self.n_Hrepresentation())
        if self.n_vertices():
            # Depending on the backend, this does not hold for the empty polyhedron.
            tester.assertEqual(self.dim() + self.n_equations(), self.ambient_dim())

        tester.assertTrue(all(len(v[::]) == self.ambient_dim() for v in self.Vrep_generator()))
        tester.assertTrue(all(len(h[::]) == self.ambient_dim() + 1 for h in self.Hrep_generator()))

        if self.n_vertices() + self.n_rays() < 40:
            tester.assertEqual(self, Polyhedron(vertices=self.vertices(), rays=self.rays(), lines=self.lines(), ambient_dim=self.ambient_dim()))
        if self.n_inequalities() < 40:
            tester.assertEqual(self, Polyhedron(ieqs=self.inequalities(), eqns=self.equations(), ambient_dim=self.ambient_dim()))

    def base_extend(self, base_ring, backend=None):
        """
        Return a new polyhedron over a larger base ring.

        This method can also be used to change the backend.

        INPUT:

        - ``base_ring`` -- the new base ring

        - ``backend`` -- the new backend, see
          :func:`~sage.geometry.polyhedron.constructor.Polyhedron`.
          If ``None`` (the default), attempt to keep the same backend.
          Otherwise, use the same defaulting behavior
          as described there.

        OUTPUT:

        The same polyhedron, but over a larger base ring and possibly with a changed backend.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)], base_ring=ZZ);  P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.base_extend(QQ)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.base_extend(QQ) == P
            True

        TESTS:

        Test that :trac:`22575` is fixed::

            sage: Q = P.base_extend(ZZ, backend='field')
            sage: Q.backend()
            'field'

        """
        new_parent = self.parent().base_extend(base_ring, backend)
        return new_parent(self, copy=True)

    def change_ring(self, base_ring, backend=None):
        """
        Return the polyhedron obtained by coercing the entries of the
        vertices/lines/rays of this polyhedron into the given ring.

        This method can also be used to change the backend.

        INPUT:

        - ``base_ring`` -- the new base ring

        - ``backend`` -- the new backend or ``None`` (default), see
          :func:`~sage.geometry.polyhedron.constructor.Polyhedron`.
          If ``None`` (the default), attempt to keep the same backend.
          Otherwise, use the same defaulting behavior
          as described there.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)], base_ring=QQ); P
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.change_ring(ZZ)
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.change_ring(ZZ) == P
            True

            sage: P = Polyhedron(vertices=[(-1.3,0), (0,2.3)], base_ring=RDF); P.vertices()
            (A vertex at (-1.3, 0.0), A vertex at (0.0, 2.3))
            sage: P.change_ring(QQ).vertices()
            (A vertex at (-13/10, 0), A vertex at (0, 23/10))
            sage: P == P.change_ring(QQ)
            True
            sage: P.change_ring(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: cannot change the base ring to the Integer Ring

            sage: P = polytopes.regular_polygon(3); P
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 3 vertices
            sage: P.vertices()
            (A vertex at (0.?e-16, 1.000000000000000?),
             A vertex at (0.866025403784439?, -0.500000000000000?),
             A vertex at (-0.866025403784439?, -0.500000000000000?))
            sage: P.change_ring(QQ)
            Traceback (most recent call last):
            ...
            TypeError: cannot change the base ring to the Rational Field

        .. WARNING::

            The base ring ``RDF`` should be used with care. As it is
            not an exact ring, certain computations may break or
            silently produce wrong results, for example changing the
            base ring from an exact ring into ``RDF`` may cause a
            loss of data::

                sage: P = Polyhedron([[2/3,0],[6666666666666667/10^16,0]], base_ring=AA); P
                A 1-dimensional polyhedron in AA^2 defined as the convex hull of 2 vertices
                sage: Q = P.change_ring(RDF); Q
                A 0-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex
                sage: P.n_vertices() == Q.n_vertices()
                False
       """

        from sage.categories.all import Rings

        if base_ring not in Rings:
            raise ValueError("invalid base ring")

        try:
            vertices = [[base_ring(x) for x in vertex] for vertex in self.vertices_list()]
            rays = [[base_ring(x) for x in ray] for ray in self.rays_list()]
            lines = [[base_ring(x) for x in line] for line in self.lines_list()]

        except (TypeError, ValueError):
            raise TypeError("cannot change the base ring to the {0}".format(base_ring))

        new_parent = self.parent().change_ring(base_ring, backend)
        return new_parent([vertices, rays, lines], None)

    def _richcmp_(self, other, op):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- a polyhedron

        OUTPUT:

        If ``other`` is a polyhedron, then the comparison
        operator "less or equal than" means "is contained in", and
        "less than" means "is strictly contained in".

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P >= Q
            True
            sage: Q <= P
            True
            sage: P == P
            True

       The polytope ``Q`` is strictly contained in ``P``::

            sage: P > Q
            True
            sage: P < Q
            False
            sage: P == Q
            False

        Test that we have fixed a problem revealed in :trac:`31701`,
        where neither of the two polyhedra contains the other::

            sage: P = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: Q = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: Q < P
            False
            sage: P > Q
            False
         """
        if self.Vrepresentation() is None or other.Vrepresentation() is None:
            raise RuntimeError('some V representation is missing')
            # make sure deleted polyhedra are not used in cache

        if self.ambient_dim() != other.ambient_dim():
            return op == op_NE

        c0 = self._is_subpolyhedron(other)
        c1 = other._is_subpolyhedron(self)
        if c0 and c1:
            return rich_to_bool(op, 0)
        elif c0:
            return rich_to_bool(op, -1)
        elif c1:
            return rich_to_bool(op, 1)
        else:
            return op == op_NE

    @coerce_binop
    def _is_subpolyhedron(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`

        OUTPUT:

        Boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P._is_subpolyhedron(Q)
            False
            sage: Q._is_subpolyhedron(P)
            True
        """
        return all(other_H.contains(self_V)
                   for other_H in other.Hrepresentation()
                   for self_V in self.Vrepresentation())

    def is_mutable(self):
        r"""
        Return True if the polyhedron is mutable, i.e. it can be modified in place.

        EXAMPLES::

            sage: p = polytopes.cube(backend='field')
            sage: p.is_mutable()
            False
        """
        return False

    def is_immutable(self):
        r"""
        Return True if the polyhedron is immutable, i.e. it cannot be modified in place.

        EXAMPLES::

            sage: p = polytopes.cube(backend='field')
            sage: p.is_immutable()
            True
        """
        return True

    @cached_method
    def vertex_facet_graph(self, labels=True):
        r"""
        Return the vertex-facet graph.

        This function constructs a directed bipartite graph.
        The nodes of the graph correspond to the vertices of the polyhedron
        and the facets of the polyhedron. There is an directed edge
        from a vertex to a face if and only if the vertex is incident to the face.

        INPUT:

        - ``labels`` -- boolean (default: ``True``); decide how the nodes
          of the graph are labelled. Either with the original vertices/facets
          of the Polyhedron or with integers.

        OUTPUT:

        - a bipartite DiGraph. If ``labels`` is ``True``, then the nodes
          of the graph will actually be the vertices and facets of ``self``,
          otherwise they will be integers.

        .. SEEALSO::

            :meth:`combinatorial_automorphism_group`,
            :meth:`is_combinatorially_isomorphic`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: G = P.vertex_facet_graph(); G
            Digraph on 14 vertices
            sage: G.vertices(key = lambda v: str(v))
            [A vertex at (-1, -1, -1),
             A vertex at (-1, -1, 1),
             A vertex at (-1, 1, -1),
             A vertex at (-1, 1, 1),
             A vertex at (1, -1, -1),
             A vertex at (1, -1, 1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             An inequality (-1, 0, 0) x + 1 >= 0,
             An inequality (0, -1, 0) x + 1 >= 0,
             An inequality (0, 0, -1) x + 1 >= 0,
             An inequality (0, 0, 1) x + 1 >= 0,
             An inequality (0, 1, 0) x + 1 >= 0,
             An inequality (1, 0, 0) x + 1 >= 0]
            sage: G.automorphism_group().is_isomorphic(P.hasse_diagram().automorphism_group())
            True
            sage: O = polytopes.octahedron(); O
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: O.vertex_facet_graph()
            Digraph on 14 vertices
            sage: H = O.vertex_facet_graph()
            sage: G.is_isomorphic(H)
            False
            sage: G2 = copy(G)
            sage: G2.reverse_edges(G2.edges())
            sage: G2.is_isomorphic(H)
            True

        TESTS:

        Check that :trac:`28828` is fixed::

            sage: G._immutable
            True

        Check that :trac:`29188` is fixed::

            sage: P = polytopes.cube()
            sage: P.vertex_facet_graph().is_isomorphic(P.vertex_facet_graph(False))
            True
        """
        return self.combinatorial_polyhedron().vertex_facet_graph(names=labels)

    def plot(self,
             point=None, line=None, polygon=None,  # None means unspecified by the user
             wireframe='blue', fill='green',
             position=None,
             orthonormal=True,  # whether to use orthonormal projections
             **kwds):
        """
        Return a graphical representation.

        INPUT:

        - ``point``, ``line``, ``polygon`` -- Parameters to pass to
          point (0d), line (1d), and polygon (2d) plot commands.
          Allowed values are:

          * A Python dictionary to be passed as keywords to the plot
            commands.

          * A string or triple of numbers: The color. This is
            equivalent to passing the dictionary ``{'color':...}``.

          * ``False``: Switches off the drawing of the corresponding
            graphics object

        - ``wireframe``, ``fill`` -- Similar to ``point``, ``line``,
          and ``polygon``, but ``fill`` is used for the graphics
          objects in the dimension of the polytope (or of dimension 2
          for higher dimensional polytopes) and ``wireframe`` is used
          for all lower-dimensional graphics objects
          (default: 'green' for ``fill`` and 'blue' for ``wireframe``)

        - ``position`` -- positive number; the position to take the projection
          point in Schlegel diagrams.

        - ``orthonormal`` -- Boolean (default: True); whether to use
          orthonormal projections.

        - ``**kwds`` -- optional keyword parameters that are passed to
          all graphics objects.

        OUTPUT:

        A (multipart) graphics object.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: point = Polyhedron([[1,1]])
            sage: line = Polyhedron([[1,1],[2,1]])
            sage: cube = polytopes.hypercube(3)
            sage: hypercube = polytopes.hypercube(4)

        By default, the wireframe is rendered in blue and the fill in green::

            sage: square.plot()
            Graphics object consisting of 6 graphics primitives
            sage: point.plot()
            Graphics object consisting of 1 graphics primitive
            sage: line.plot()
            Graphics object consisting of 2 graphics primitives
            sage: cube.plot()
            Graphics3d Object
            sage: hypercube.plot()
            Graphics3d Object

        Draw the lines in red and nothing else::

            sage: square.plot(point=False, line='red', polygon=False)
            Graphics object consisting of 4 graphics primitives
            sage: point.plot(point=False, line='red', polygon=False)
            Graphics object consisting of 0 graphics primitives
            sage: line.plot(point=False, line='red', polygon=False)
            Graphics object consisting of 1 graphics primitive
            sage: cube.plot(point=False, line='red', polygon=False)
            Graphics3d Object
            sage: hypercube.plot(point=False, line='red', polygon=False)
            Graphics3d Object

        Draw points in red, no lines, and a blue polygon::

            sage: square.plot(point={'color':'red'}, line=False, polygon=(0,0,1))
            Graphics object consisting of 2 graphics primitives
            sage: point.plot(point={'color':'red'}, line=False, polygon=(0,0,1))
            Graphics object consisting of 1 graphics primitive
            sage: line.plot(point={'color':'red'}, line=False, polygon=(0,0,1))
            Graphics object consisting of 1 graphics primitive
            sage: cube.plot(point={'color':'red'}, line=False, polygon=(0,0,1))
            Graphics3d Object
            sage: hypercube.plot(point={'color':'red'}, line=False, polygon=(0,0,1))
            Graphics3d Object

        If we instead use the ``fill`` and ``wireframe`` options, the
        coloring depends on the dimension of the object::

            sage: square.plot(fill='green', wireframe='red')
            Graphics object consisting of 6 graphics primitives
            sage: point.plot(fill='green', wireframe='red')
            Graphics object consisting of 1 graphics primitive
            sage: line.plot(fill='green', wireframe='red')
            Graphics object consisting of 2 graphics primitives
            sage: cube.plot(fill='green', wireframe='red')
            Graphics3d Object
            sage: hypercube.plot(fill='green', wireframe='red')
            Graphics3d Object

        It is possible to draw polyhedra up to dimension 4, no matter what the
        ambient dimension is::

            sage: hcube = polytopes.hypercube(5)
            sage: facet = hcube.facets()[0].as_polyhedron();facet
            A 4-dimensional polyhedron in ZZ^5 defined as the convex hull of 16 vertices
            sage: facet.plot()
            Graphics3d Object

        TESTS::

            sage: for p in square.plot():
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            blue Point set defined by 4 point(s)
            blue Line defined by 2 points
            blue Line defined by 2 points
            blue Line defined by 2 points
            blue Line defined by 2 points
            green Polygon defined by 4 points

            sage: for p in line.plot():
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            blue Point set defined by 2 point(s)
            green Line defined by 2 points

            sage: for p in point.plot():
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            green Point set defined by 1 point(s)

        Draw the lines in red and nothing else::

            sage: for p in square.plot(point=False, line='red', polygon=False):
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Line defined by 2 points
            red Line defined by 2 points
            red Line defined by 2 points
            red Line defined by 2 points

        Draw vertices in red, no lines, and a blue polygon::

            sage: for p in square.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 4 point(s)
            (0, 0, 1) Polygon defined by 4 points

            sage: for p in line.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 2 point(s)

            sage: for p in point.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 1 point(s)

        Draw in red without wireframe::

            sage: for p in square.plot(wireframe=False, fill="red"):
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Polygon defined by 4 points

            sage: for p in line.plot(wireframe=False, fill="red"):
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Line defined by 2 points

            sage: for p in point.plot(wireframe=False, fill="red"):
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 1 point(s)

        We try to draw the polytope in 2 or 3 dimensions::

            sage: type(Polyhedron(ieqs=[(1,)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(polytopes.hypercube(1).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(polytopes.hypercube(2).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(polytopes.hypercube(3).plot())
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>

        In 4d a projection to 3d is used::

            sage: type(polytopes.hypercube(4).plot())
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
            sage: type(polytopes.hypercube(5).plot())
            Traceback (most recent call last):
            ...
            NotImplementedError: plotting of 5-dimensional polyhedra not implemented

        If the polyhedron is not full-dimensional, the :meth:`affine_hull_projection` is used if necessary::

            sage: type(Polyhedron([(0,), (1,)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(Polyhedron([(0,0), (1,1)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(Polyhedron([(0,0,0), (1,1,1)]).plot())
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
            sage: type(Polyhedron([(0,0,0,0), (1,1,1,1)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(Polyhedron([(0,0,0,0,0), (1,1,1,1,1)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(Polyhedron([(0,0,0,0), (1,1,1,1), (1,0,0,0)]).plot())
            <class 'sage.plot.graphics.Graphics'>

        TESTS:

        Check that :trac:`30015` is fixed::

            sage: fcube = polytopes.hypercube(4)
            sage: tfcube = fcube.face_truncation(fcube.faces(0)[0])
            sage: sp = tfcube.schlegel_projection()
            sage: for face in tfcube.faces(2):
            ....:     vertices = face.ambient_Vrepresentation()
            ....:     indices = [sp.coord_index_of(vector(x)) for x in vertices]
            ....:     projected_vertices = [sp.transformed_coords[i] for i in indices]
            ....:     assert Polyhedron(projected_vertices).dim() == 2
        """
        def merge_options(*opts):
            merged = dict()
            for i in range(len(opts)):
                opt = opts[i]
                if opt is None:
                    continue
                elif opt is False:
                    return False
                elif isinstance(opt, (str, list, tuple)):
                    merged['color'] = opt
                else:
                    merged.update(opt)
            return merged

        d = min(self.dim(), 2)
        opts = [wireframe] * d + [fill] + [False] * (2-d)
        # The point/line/polygon options take precedence over wireframe/fill
        opts = [merge_options(opt1, opt2, kwds)
                for opt1, opt2 in zip(opts, [point, line, polygon])]

        def project(polyhedron, ortho):
            if polyhedron.ambient_dim() <= 3:
                return polyhedron.projection()
            elif polyhedron.dim() <= 3:
                if ortho:
                    return polyhedron.affine_hull_projection(orthonormal=True, extend=True).projection()
                else:
                    return polyhedron.affine_hull_projection().projection()
            elif polyhedron.dimension() == 4:
                # For 4d-polyhedron, we can use schlegel projections:
                return polyhedron.schlegel_projection(position=position)
            else:
                return polyhedron.projection()

        projection = project(self, orthonormal)
        try:
            plot_method = projection.plot
        except AttributeError:
            raise NotImplementedError('plotting of {0}-dimensional polyhedra not implemented'
                                          .format(self.ambient_dim()))
        return plot_method(*opts)

    def show(self, **kwds):
        """
        Display graphics immediately

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        INPUT:

        - ``kwds`` -- optional keyword arguments. See :meth:`plot` for
          the description of available options.

        OUTPUT:

        This method does not return anything. Use :meth:`plot` if you
        want to generate a graphics object that can be saved or
        further transformed.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: square.show(point='red')
        """
        self.plot(**kwds).show()

    def tikz(self, view=[0, 0, 1], angle=0, scale=1,
             edge_color='blue!95!black', facet_color='blue!95!black',
             opacity=0.8, vertex_color='green', axis=False):
        r"""
        Return a string ``tikz_pic`` consisting of a tikz picture of ``self``
        according to a projection ``view`` and an angle ``angle``
        obtained via the threejs viewer.

        INPUT:

        - ``view`` - list (default: [0,0,1]) representing the rotation axis (see note below).
        - ``angle`` - integer (default: 0) angle of rotation in degree from 0 to 360 (see note
          below).
        - ``scale`` - integer (default: 1) specifying the scaling of the tikz picture.
        - ``edge_color`` - string (default: 'blue!95!black') representing colors which tikz
          recognize.
        - ``facet_color`` - string (default: 'blue!95!black') representing colors which tikz
          recognize.
        - ``vertex_color`` - string (default: 'green') representing colors which tikz
          recognize.
        - ``opacity`` - real number (default: 0.8) between 0 and 1 giving the opacity of
          the front facets.
        - ``axis`` - Boolean (default: False) draw the axes at the origin or not.

        OUTPUT:

        - LatexExpr -- containing the TikZ picture.

        .. NOTE::

            This is a wrapper of a method of the projection object
            `self.projection()`. See :meth:`~sage.geometry.polyhedron.plot.Projection.tikz`
            for more detail.

            The inputs ``view`` and ``angle`` can be obtained by visualizing it
            using ``.show(aspect_ratio=1)``. This will open an interactive view
            in your default browser, where you can rotate the polytope. Once
            the desired view angle is found, click on the information icon in
            the lower right-hand corner and select *Get Viewpoint*. This will
            copy a string of the form '[x,y,z],angle' to your local clipboard.
            Go back to Sage and type ``Img = P.tikz([x,y,z],angle)``.

            The inputs ``view`` and ``angle`` can also be obtained from the
            viewer Jmol::

                1) Right click on the image
                2) Select ``Console``
                3) Select the tab ``State``
                4) Scroll to the line ``moveto``

            It reads something like::

                moveto 0.0 {x y z angle} Scale

            The ``view`` is then [x,y,z] and ``angle`` is angle.
            The following number is the scale.

            Jmol performs a rotation of ``angle`` degrees along the
            vector [x,y,z] and show the result from the z-axis.


        EXAMPLES::

            sage: co = polytopes.cuboctahedron()
            sage: Img = co.tikz([0,0,1], 0)
            sage: print('\n'.join(Img.splitlines()[:9]))
            \begin{tikzpicture}%
                [x={(1.000000cm, 0.000000cm)},
                y={(0.000000cm, 1.000000cm)},
                z={(0.000000cm, 0.000000cm)},
                scale=1.000000,
                back/.style={loosely dotted, thin},
                edge/.style={color=blue!95!black, thick},
                facet/.style={fill=blue!95!black,fill opacity=0.800000},
                vertex/.style={inner sep=1pt,circle,draw=green!25!black,fill=green!75!black,thick}]
            sage: print('\n'.join(Img.splitlines()[12:21]))
            %% with the command: ._tikz_3d_in_3d and parameters:
            %% view = [0, 0, 1]
            %% angle = 0
            %% scale = 1
            %% edge_color = blue!95!black
            %% facet_color = blue!95!black
            %% opacity = 0.8
            %% vertex_color = green
            %% axis = False
            sage: print('\n'.join(Img.splitlines()[22:26]))
            %% Coordinate of the vertices:
            %%
            \coordinate (-1.00000, -1.00000, 0.00000) at (-1.00000, -1.00000, 0.00000);
            \coordinate (-1.00000, 0.00000, -1.00000) at (-1.00000, 0.00000, -1.00000);
        """
        return self.projection().tikz(view, angle, scale,
                                      edge_color, facet_color,
                                      opacity, vertex_color, axis)

    def _repr_(self):
        """
        Return a description of the polyhedron.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1]])
            sage: poly_test._repr_()
            'A 2-dimensional polyhedron in ZZ^4 defined as the convex hull of 3 vertices'
            sage: grammar_test = Polyhedron(vertices = [[1,1,1,1,1,1]])
            sage: grammar_test._repr_()
            'A 0-dimensional polyhedron in ZZ^6 defined as the convex hull of 1 vertex'
        """
        desc = ''
        if self.n_vertices() == 0:
            desc += 'The empty polyhedron'
        else:
            desc += 'A ' + repr(self.dim()) + '-dimensional polyhedron'
        desc += ' in '
        desc += self.parent()._repr_ambient_module()

        if self.n_vertices() > 0:
            desc += ' defined as the convex hull of '
            desc += repr(self.n_vertices())
            if self.n_vertices() == 1:
                desc += ' vertex'
            else:
                desc += ' vertices'

            if self.n_rays() > 0:
                if self.n_lines() > 0:
                    desc += ", "
                else:
                    desc += " and "
                desc += repr(self.n_rays())
                if self.n_rays() == 1:
                    desc += ' ray'
                else:
                    desc += ' rays'

            if self.n_lines() > 0:
                if self.n_rays() > 0:
                    desc += ", "
                else:
                    desc += " and "
                desc += repr(self.n_lines())
                if self.n_lines() == 1:
                    desc += ' line'
                else:
                    desc += ' lines'

        return desc

    def _rich_repr_(self, display_manager, **kwds):
        r"""
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: polytopes.hypercube(2)._rich_repr_(dm)
            OutputPlainText container

        The ``supplemental_plot`` preference lets us control whether
        this object is shown as text or picture+text::

            sage: dm.preferences.supplemental_plot
            'never'
            sage: del dm.preferences.supplemental_plot
            sage: polytopes.hypercube(3)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices (use the .plot() method to plot)
            sage: dm.preferences.supplemental_plot = 'never'
        """
        prefs = display_manager.preferences
        is_small = (self.ambient_dim() <= 2)
        can_plot = (prefs.supplemental_plot != 'never')
        plot_graph = can_plot and (prefs.supplemental_plot == 'always' or is_small)
        # Under certain circumstances we display the plot as graphics
        if plot_graph:
            plot_kwds = dict(kwds)
            plot_kwds.setdefault('title', repr(self))
            output = self.plot(**plot_kwds)._rich_repr_(display_manager)
            if output is not None:
                return output
        # create text for non-graphical output
        if can_plot:
            text = '{0} (use the .plot() method to plot)'.format(repr(self))
        else:
            text = repr(self)
        # latex() produces huge tikz environment, override
        tp = display_manager.types
        if (prefs.text == 'latex' and tp.OutputLatex in display_manager.supported_output()):
            return tp.OutputLatex(r'\text{{{0}}}'.format(text))
        return tp.OutputPlainText(text)

    def cdd_Hrepresentation(self):
        r"""
        Write the inequalities/equations data of the polyhedron in
        cdd's H-representation format.

        .. SEEALSO::

            :meth:`write_cdd_Hrepresentation` -- export the polyhedron as a
            H-representation to a file.

        OUTPUT: a string

        EXAMPLES::

            sage: p = polytopes.hypercube(2)
            sage: print(p.cdd_Hrepresentation())
            H-representation
            begin
             4 3 rational
             1 -1 0
             1 0 -1
             1 1 0
             1 0 1
            end
            <BLANKLINE>

            sage: triangle = Polyhedron(vertices = [[1,0],[0,1],[1,1]],base_ring=AA)
            sage: triangle.base_ring()
            Algebraic Real Field
            sage: triangle.cdd_Hrepresentation()
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be ZZ, QQ, or RDF
        """
        from .cdd_file_format import cdd_Hrepresentation
        try:
            cdd_type = self._cdd_type
        except AttributeError:
            if self.base_ring() is ZZ or self.base_ring() is QQ:
                cdd_type = 'rational'
            elif self.base_ring() is RDF:
                cdd_type = 'real'
            else:
                raise TypeError('the base ring must be ZZ, QQ, or RDF')
        return cdd_Hrepresentation(cdd_type,
                                   list(self.inequality_generator()),
                                   list(self.equation_generator()))

    def write_cdd_Hrepresentation(self, filename):
        r"""
        Export the polyhedron as a H-representation to a file.

        INPUT:

        - ``filename`` -- the output file.

        .. SEEALSO::

            :meth:`cdd_Hrepresentation` -- return the H-representation of the
            polyhedron as a string.

        EXAMPLES::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename(ext='.ext')
            sage: polytopes.cube().write_cdd_Hrepresentation(filename)
        """
        with open(filename, 'w') as f:
            f.write(self.cdd_Hrepresentation())

    def cdd_Vrepresentation(self):
        r"""
        Write the vertices/rays/lines data of the polyhedron in cdd's
        V-representation format.

        .. SEEALSO::

            :meth:`write_cdd_Vrepresentation` -- export the polyhedron as a
            V-representation to a file.

        OUTPUT: a string

        EXAMPLES::

            sage: q = Polyhedron(vertices = [[1,1],[0,0],[1,0],[0,1]])
            sage: print(q.cdd_Vrepresentation())
            V-representation
            begin
             4 3 rational
             1 0 0
             1 0 1
             1 1 0
             1 1 1
            end
        """
        from .cdd_file_format import cdd_Vrepresentation
        try:
            cdd_type = self._cdd_type
        except AttributeError:
            if self.base_ring() is ZZ or self.base_ring() is QQ:
                cdd_type = 'rational'
            elif self.base_ring() is RDF:
                cdd_type = 'real'
            else:
                raise TypeError('the base ring must be ZZ, QQ, or RDF')
        return cdd_Vrepresentation(cdd_type,
                                   list(self.vertex_generator()),
                                   list(self.ray_generator()),
                                   list(self.line_generator()))

    def write_cdd_Vrepresentation(self, filename):
        r"""
        Export the polyhedron as a V-representation to a file.

        INPUT:

        - ``filename`` -- the output file.

        .. SEEALSO::

            :meth:`cdd_Vrepresentation` -- return the V-representation of the
            polyhedron as a string.

        EXAMPLES::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename(ext='.ext')
            sage: polytopes.cube().write_cdd_Vrepresentation(filename)
        """
        with open(filename, 'w') as f:
            f.write(self.cdd_Vrepresentation())

    @cached_method
    def n_equations(self):
        """
        Return the number of equations. The representation will
        always be minimal, so the number of equations is the
        codimension of the polyhedron in the ambient space.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0,0],[0,1,0],[0,0,1]])
            sage: p.n_equations()
            1
        """
        return len(self.equations())

    @cached_method
    def n_inequalities(self):
        """
        Return the number of inequalities. The representation will
        always be minimal, so the number of inequalities is the
        number of facets of the polyhedron in the ambient space.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0,0],[0,1,0],[0,0,1]])
            sage: p.n_inequalities()
            3

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)])
            sage: p.n_facets()
            8
        """
        return len(self.inequalities())

    n_facets = n_inequalities

    @cached_method
    def n_vertices(self):
        """
        Return the number of vertices. The representation will
        always be minimal.

        .. WARNING::

            If the polyhedron has lines, return the number of vertices in
            the ``Vrepresentation``. As the represented polyhedron has
            no 0-dimensional faces (i.e. vertices), ``n_vertices`` corresponds
            to the number of `k`-faces, where `k` is the number of lines::

                sage: P = Polyhedron(rays=[[1,0,0]],lines=[[0,1,0]])
                sage: P.n_vertices()
                1
                sage: P.faces(0)
                ()
                sage: P.f_vector()
                (1, 0, 1, 1)

                sage: P = Polyhedron(rays=[[1,0,0]],lines=[[0,1,0],[0,1,1]])
                sage: P.n_vertices()
                1
                sage: P.f_vector()
                (1, 0, 0, 1, 1)

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]])
            sage: p.n_vertices()
            2
        """
        return len(self.vertices())

    @cached_method
    def n_rays(self):
        """
        Return the number of rays. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]])
            sage: p.n_rays()
            1
        """
        return len(self.rays())

    @cached_method
    def n_lines(self):
        """
        Return the number of lines. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]])
            sage: p.n_lines()
            1
        """
        return len(self.lines())

    def to_linear_program(self, solver=None, return_variable=False, base_ring=None):
        r"""
        Return a linear optimization problem over the polyhedron in the form of
        a :class:`MixedIntegerLinearProgram`.

        INPUT:

        - ``solver`` -- select a solver (MIP backend). See the documentation
          of for :class:`MixedIntegerLinearProgram`. Set to ``None`` by default.

        - ``return_variable`` -- (default: ``False``) If ``True``, return a tuple
          ``(p, x)``, where ``p`` is the :class:`MixedIntegerLinearProgram` object
          and ``x`` is the vector-valued MIP variable in this problem, indexed
          from 0.  If ``False``, only return ``p``.

        - ``base_ring`` -- select a field over which the linear program should be
          set up.  Use ``RDF`` to request a fast inexact (floating point) solver
          even if ``self`` is exact.

        Note that the :class:`MixedIntegerLinearProgram` object will have the
        null function as an objective to be maximized.

        .. SEEALSO::

            :meth:`~MixedIntegerLinearProgram.polyhedron` -- return the
            polyhedron associated with a :class:`MixedIntegerLinearProgram`
            object.

        EXAMPLES:

        Exact rational linear program::

            sage: p = polytopes.cube()
            sage: p.to_linear_program()
            Linear Program (no objective, 3 variables, 6 constraints)
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(2*x[0] + 1*x[1] + 39*x[2])
            sage: lp.solve()
            42
            sage: lp.get_values(x[0], x[1], x[2])
            [1, 1, 1]

        Floating-point linear program::

            sage: lp, x = p.to_linear_program(return_variable=True, base_ring=RDF)
            sage: lp.set_objective(2*x[0] + 1*x[1] + 39*x[2])
            sage: lp.solve()
            42.0

        Irrational algebraic linear program over an embedded number field::

            sage: p=polytopes.icosahedron()
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve()
            1/4*sqrt5 + 3/4

        Same example with floating point::

            sage: lp, x = p.to_linear_program(return_variable=True, base_ring=RDF)
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve() # tol 1e-5
            1.3090169943749475

        Same example with a specific floating point solver::

            sage: lp, x = p.to_linear_program(return_variable=True, solver='GLPK')
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve() # tol 1e-8
            1.3090169943749475

        Irrational algebraic linear program over `AA`::

            sage: p=polytopes.icosahedron(base_ring=AA)
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve()  # long time
            1.309016994374948?

        TESTS::

            sage: p=polytopes.flow_polytope(digraphs.DeBruijn(3,2)); p
            A 19-dimensional polyhedron in QQ^27 defined as the convex hull of 1 vertex and 148 rays
            sage: p.to_linear_program().polyhedron() == p
            True
            sage: p=polytopes.icosahedron()
            sage: p.to_linear_program(solver='PPL')
            Traceback (most recent call last):
            ...
            TypeError: The PPL backend only supports rational data.

        Test that equations are handled correctly (:trac:`24154`)::

            sage: p = Polyhedron(vertices=[[19]])
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(x[0])
            sage: lp.solve()
            19
        """
        if base_ring is None:
            base_ring = self.base_ring()
        base_ring = base_ring.fraction_field()
        from sage.numerical.mip import MixedIntegerLinearProgram
        p = MixedIntegerLinearProgram(solver=solver, base_ring=base_ring)
        x = p.new_variable(real=True, nonnegative=False)

        for ineqn in self.inequalities_list():
            b = -ineqn.pop(0)
            p.add_constraint(p.sum([x[i]*ineqn[i] for i in range(len(ineqn))]) >= b)

        for eqn in self.equations_list():
            b = -eqn.pop(0)
            p.add_constraint(p.sum([x[i]*eqn[i] for i in range(len(eqn))]) == b)

        if return_variable:
            return p, x
        else:
            return p

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

            sage: p = polytopes.hypercube(3, backend='field')
            sage: p.Hrepresentation(0)
            An inequality (-1, 0, 0) x + 1 >= 0
            sage: p.Hrepresentation(0) == p.Hrepresentation()[0]
            True
        """
        if index is None:
            return self._Hrepresentation
        else:
            return self._Hrepresentation[index]

    def Hrepresentation_str(self, separator='\n', latex=False, style='>=', align=None, **kwds):
        r"""
        Return a human-readable string representation of the Hrepresentation of this
        polyhedron.

        INPUT:

        - ``separator`` -- a string. Default is ``"\n"``.

        - ``latex`` -- a boolean. Default is ``False``.

        - ``style`` -- either ``"positive"`` (making all coefficients positive)
                       or ``"<="``, or ``">="``. Default is ``">="``.

        - ``align`` -- a boolean or ``None''. Default is ``None`` in which case
                       ``align`` is ``True`` if ``separator`` is the newline character.
                       If set, then the lines of the output string are aligned
                       by the comparison symbol by padding blanks.

        Keyword parameters of
        :meth:`~sage.geometry.polyhedron.representation.Hrepresentation.repr_pretty`
        are passed on:

        - ``prefix`` -- a string

        - ``indices`` -- a tuple or other iterable

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P = polytopes.permutahedron(3)
            sage: print(P.Hrepresentation_str())
            x0 + x1 + x2 ==  6
                 x0 + x1 >=  3
                -x0 - x1 >= -5
                      x1 >=  1
                     -x0 >= -3
                      x0 >=  1
                     -x1 >= -3

            sage: print(P.Hrepresentation_str(style='<='))
            -x0 - x1 - x2 == -6
                 -x0 - x1 <= -3
                  x0 + x1 <=  5
                      -x1 <= -1
                       x0 <=  3
                      -x0 <= -1
                       x1 <=  3

            sage: print(P.Hrepresentation_str(style='positive'))
            x0 + x1 + x2 == 6
                 x0 + x1 >= 3
                       5 >= x0 + x1
                      x1 >= 1
                       3 >= x0
                      x0 >= 1
                       3 >= x1

            sage: print(P.Hrepresentation_str(latex=True))
            \begin{array}{rcl}
            x_{0} + x_{1} + x_{2} & =    &  6 \\
                    x_{0} + x_{1} & \geq &  3 \\
                   -x_{0} - x_{1} & \geq & -5 \\
                            x_{1} & \geq &  1 \\
                           -x_{0} & \geq & -3 \\
                            x_{0} & \geq &  1 \\
                           -x_{1} & \geq & -3
            \end{array}

            sage: print(P.Hrepresentation_str(align=False))
            x0 + x1 + x2 == 6
            x0 + x1 >= 3
            -x0 - x1 >= -5
            x1 >= 1
            -x0 >= -3
            x0 >= 1
            -x1 >= -3

            sage: c = polytopes.cube()
            sage: c.Hrepresentation_str(separator=', ', style='positive')
            '1 >= x0, 1 >= x1, 1 >= x2, x0 + 1 >= 0, x2 + 1 >= 0, x1 + 1 >= 0'
        """
        pretty_hs = [h.repr_pretty(split=True, latex=latex, style=style, **kwds) for h in self.Hrepresentation()]
        shift = any(pretty_h[2].startswith('-') for pretty_h in pretty_hs)

        if align is None:
            align = separator == "\n"
        if align:
            lengths = [(len(s[0]), len(s[1]), len(s[2])) for s in pretty_hs]
            from operator import itemgetter
            length_left = max(lengths, key=itemgetter(0))[0]
            length_middle = max(lengths, key=itemgetter(1))[1]
            length_right = max(lengths, key=itemgetter(2))[2]
            if shift:
                length_right += 1
            if latex:
                h_line = "{:>" + "{}".format(length_left) + "} & {:" + \
                         "{}".format(length_middle) + "} & {:" + \
                         "{}".format(length_right) + "}\\\\"
            else:
                h_line = "{:>" + "{}".format(length_left) \
                         + "} {:" + "{}".format(length_middle) \
                         + "} {:" + "{}".format(length_right) + "}"
        elif latex:
            h_line = "{} & {} & {}\\\\"
        else:
            h_line = "{} {} {}"

        def pad_non_minus(s):
            if align and shift and not s.startswith('-'):
                return ' ' + s
            else:
                return s
        h_list = [h_line.format(pretty_h[0], pretty_h[1], pad_non_minus(pretty_h[2]))
                  for pretty_h in pretty_hs]
        pretty_print = separator.join(h_list)

        if not latex:
            return pretty_print
        else:
            # below we remove the 2 unnecessary backslashes at the end of pretty_print
            return "\\begin{array}{rcl}\n" + pretty_print[:-2] + "\n\\end{array}"

    def Hrep_generator(self):
        """
        Return an iterator over the objects of the H-representation
        (inequalities or equations).

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: next(p.Hrep_generator())
            An inequality (-1, 0, 0) x + 1 >= 0
        """
        for H in self.Hrepresentation():
            yield H

    @cached_method
    def n_Hrepresentation(self):
        """
        Return the number of objects that make up the
        H-representation of the polyhedron.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: p.n_Hrepresentation()
            16
            sage: p.n_Hrepresentation() == p.n_inequalities() + p.n_equations()
            True
        """
        return len(self.Hrepresentation())

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

            sage: p = polytopes.simplex(4, project=True)
            sage: p.Vrepresentation(0)
            A vertex at (0.7071067812, 0.4082482905, 0.2886751346, 0.2236067977)
            sage: p.Vrepresentation(0) == p.Vrepresentation() [0]
            True
        """
        if index is None:
            return self._Vrepresentation
        else:
            return self._Vrepresentation[index]

    @cached_method
    def n_Vrepresentation(self):
        """
        Return the number of objects that make up the
        V-representation of the polyhedron.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.simplex(4)
            sage: p.n_Vrepresentation()
            5
            sage: p.n_Vrepresentation() == p.n_vertices() + p.n_rays() + p.n_lines()
            True
        """
        return len(self.Vrepresentation())

    def Vrep_generator(self):
        """
        Return an iterator over the objects of the V-representation
        (vertices, rays, and lines).

        EXAMPLES::

            sage: p = polytopes.cyclic_polytope(3,4)
            sage: vg = p.Vrep_generator()
            sage: next(vg)
            A vertex at (0, 0, 0)
            sage: next(vg)
            A vertex at (1, 1, 1)
        """
        for V in self.Vrepresentation():
            yield V

    def inequality_generator(self):
        """
        Return  a generator for the defining inequalities of the
        polyhedron.

        OUTPUT:

        A generator of the inequality Hrepresentation objects.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.inequality_generator(): print(v)
            An inequality (1, 1) x - 1 >= 0
            An inequality (0, -1) x + 1 >= 0
            An inequality (-1, 0) x + 1 >= 0
            sage: [ v for v in triangle.inequality_generator() ]
            [An inequality (1, 1) x - 1 >= 0,
             An inequality (0, -1) x + 1 >= 0,
             An inequality (-1, 0) x + 1 >= 0]
            sage: [ [v.A(), v.b()] for v in triangle.inequality_generator() ]
            [[(1, 1), -1], [(0, -1), 1], [(-1, 0), 1]]
        """
        for H in self.Hrepresentation():
            if H.is_inequality():
                yield H

    @cached_method
    def inequalities(self):
        """
        Return all inequalities.

        OUTPUT:

        A tuple of inequalities.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[1,0,0],[2,2,2]])
            sage: p.inequalities()[0:3]
            (An inequality (1, 0, 0) x + 0 >= 0,
             An inequality (0, 1, 0) x + 0 >= 0,
             An inequality (0, 0, 1) x + 0 >= 0)
            sage: p3 = Polyhedron(vertices = Permutations([1,2,3,4]))
            sage: ieqs = p3.inequalities()
            sage: ieqs[0]
            An inequality (0, 1, 1, 1) x - 6 >= 0
            sage: list(_)
            [-6, 0, 1, 1, 1]
        """
        return tuple(self.inequality_generator())

    def inequalities_list(self):
        """
        Return a list of inequalities as coefficient lists.

        .. NOTE::

            It is recommended to use :meth:`inequalities` or
            :meth:`inequality_generator` instead to iterate over the
            list of :class:`Inequality` objects.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[1,0,0],[2,2,2]])
            sage: p.inequalities_list()[0:3]
            [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
            sage: p3 = Polyhedron(vertices = Permutations([1,2,3,4]))
            sage: ieqs = p3.inequalities_list()
            sage: ieqs[0]
            [-6, 0, 1, 1, 1]
            sage: ieqs[-1]
            [-3, 0, 1, 0, 1]
            sage: ieqs == [list(x) for x in p3.inequality_generator()]
            True
        """
        return [list(x) for x in self.inequality_generator()]

    def equation_generator(self):
        """
        Return a generator for the linear equations satisfied by the
        polyhedron.

        EXAMPLES::

            sage: p = polytopes.regular_polygon(8,base_ring=RDF)
            sage: p3 = Polyhedron(vertices = [x+[0] for x in p.vertices()], base_ring=RDF)
            sage: next(p3.equation_generator())
            An equation (0.0, 0.0, 1.0) x + 0.0 == 0
        """
        for H in self.Hrepresentation():
            if H.is_equation():
                yield H

    @cached_method
    def equations(self):
        """
        Return all linear constraints of the polyhedron.

        OUTPUT:

        A tuple of equations.

        EXAMPLES::

            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.equations()
            (An equation (1, 1, 1, 1) x - 10 == 0,)
        """
        return tuple(self.equation_generator())

    def equations_list(self):
        """
        Return the linear constraints of the polyhedron. As with
        inequalities, each constraint is given as [b -a1 -a2 ... an]
        where for variables x1, x2,..., xn, the polyhedron satisfies
        the equation b = a1*x1 + a2*x2 + ... + an*xn.

        .. NOTE::

            It is recommended to use :meth:`equations` or
            :meth:`equation_generator()` instead to iterate over the
            list of
            :class:`~sage.geometry.polyhedron.representation.Equation`
            objects.

        EXAMPLES::

            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.equations_list()
            [[-10, 1, 1, 1, 1]]
        """
        return [list(eq) for eq in self.equation_generator()]

    def vertices_list(self):
        """
        Return a list of vertices of the polyhedron.

        .. NOTE::

            It is recommended to use :meth:`vertex_generator` instead to
            iterate over the list of :class:`Vertex` objects.

        .. WARNING::

            If the polyhedron has lines, return the vertices
            of the ``Vrepresentation``. However, the represented polyhedron
            has no 0-dimensional faces (i.e. vertices)::

                sage: P = Polyhedron(rays=[[1,0,0]],lines=[[0,1,0]])
                sage: P.vertices_list()
                [[0, 0, 0]]
                sage: P.faces(0)
                ()

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: triangle.vertices_list()
            [[0, 1], [1, 0], [1, 1]]
            sage: a_simplex = Polyhedron(ieqs = [
            ....:          [0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]
            ....:      ], eqns = [[1,-1,-1,-1,-1]])
            sage: a_simplex.vertices_list()
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
            sage: a_simplex.vertices_list() == [list(v) for v in a_simplex.vertex_generator()]
            True
        """
        return [list(x) for x in self.vertex_generator()]

    def vertex_generator(self):
        """
        Return a generator for the vertices of the polyhedron.

        .. WARNING::

            If the polyhedron has lines, return a generator for the vertices
            of the ``Vrepresentation``. However, the represented polyhedron
            has no 0-dimensional faces (i.e. vertices)::

                sage: P = Polyhedron(rays=[[1,0,0]],lines=[[0,1,0]])
                sage: list(P.vertex_generator())
                [A vertex at (0, 0, 0)]
                sage: P.faces(0)
                ()

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.vertex_generator(): print(v)
            A vertex at (0, 1)
            A vertex at (1, 0)
            A vertex at (1, 1)
            sage: v_gen = triangle.vertex_generator()
            sage: next(v_gen)   # the first vertex
            A vertex at (0, 1)
            sage: next(v_gen)   # the second vertex
            A vertex at (1, 0)
            sage: next(v_gen)   # the third vertex
            A vertex at (1, 1)
            sage: try: next(v_gen)   # there are only three vertices
            ....: except StopIteration: print("STOP")
            STOP
            sage: type(v_gen)
            <... 'generator'>
            sage: [ v for v in triangle.vertex_generator() ]
            [A vertex at (0, 1), A vertex at (1, 0), A vertex at (1, 1)]
        """
        for V in self.Vrepresentation():
            if V.is_vertex():
                yield V

    @cached_method
    def vertices(self):
        """
        Return all vertices of the polyhedron.

        OUTPUT:

        A tuple of vertices.

        .. WARNING::

            If the polyhedron has lines, return the vertices
            of the ``Vrepresentation``. However, the represented polyhedron
            has no 0-dimensional faces (i.e. vertices)::

                sage: P = Polyhedron(rays=[[1,0,0]],lines=[[0,1,0]])
                sage: P.vertices()
                (A vertex at (0, 0, 0),)
                sage: P.faces(0)
                ()

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: triangle.vertices()
            (A vertex at (0, 1), A vertex at (1, 0), A vertex at (1, 1))
            sage: a_simplex = Polyhedron(ieqs = [
            ....:          [0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]
            ....:      ], eqns = [[1,-1,-1,-1,-1]])
            sage: a_simplex.vertices()
            (A vertex at (1, 0, 0, 0), A vertex at (0, 1, 0, 0),
             A vertex at (0, 0, 1, 0), A vertex at (0, 0, 0, 1))
        """
        return tuple(self.vertex_generator())

    @cached_method
    def vertices_matrix(self, base_ring=None):
        """
        Return the coordinates of the vertices as the columns of a matrix.

        INPUT:

        - ``base_ring`` -- A ring or ``None`` (default). The base ring
          of the returned matrix. If not specified, the base ring of
          the polyhedron is used.

        OUTPUT:

        A matrix over ``base_ring`` whose columns are the coordinates
        of the vertices. A ``TypeError`` is raised if the coordinates
        cannot be converted to ``base_ring``.

        .. WARNING::

            If the polyhedron has lines, return the coordinates of the vertices
            of the ``Vrepresentation``. However, the represented polyhedron
            has no 0-dimensional faces (i.e. vertices)::

                sage: P = Polyhedron(rays=[[1,0,0]],lines=[[0,1,0]])
                sage: P.vertices_matrix()
                [0]
                [0]
                [0]
                sage: P.faces(0)
                ()

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: triangle.vertices_matrix()
            [0 1 1]
            [1 0 1]
            sage: (triangle/2).vertices_matrix()
            [  0 1/2 1/2]
            [1/2   0 1/2]
            sage: (triangle/2).vertices_matrix(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        TESTS:

        Check that :trac:`28828` is fixed::

                sage: P.vertices_matrix().is_immutable()
                True
        """
        if base_ring is None:
            base_ring = self.base_ring()
        m = matrix(base_ring, self.ambient_dim(), self.n_vertices())
        for i, v in enumerate(self.vertices()):
            for j in range(self.ambient_dim()):
                m[j, i] = v[j]
        m.set_immutable()
        return m

    def an_affine_basis(self):
        """
        Return points in ``self`` that are a basis for the affine span of the polytope.

        This implementation of the method :meth:`ConvexSet_base.an_affine_basis`
        for polytopes guarantees the following:

        - All points are vertices.

        - The basis is obtained by considering a maximal chain of faces
          in the face lattice and picking for each cover relation
          one vertex that is in the difference. Thus this method
          is independent of the concrete realization of the polytope.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.an_affine_basis()
            [A vertex at (-1, -1, -1),
             A vertex at (1, -1, -1),
             A vertex at (1, -1, 1),
             A vertex at (1, 1, -1)]

            sage: P = polytopes.permutahedron(5)
            sage: P.an_affine_basis()
            [A vertex at (1, 2, 3, 5, 4),
             A vertex at (2, 1, 3, 5, 4),
             A vertex at (1, 3, 2, 5, 4),
             A vertex at (4, 1, 3, 5, 2),
             A vertex at (4, 2, 5, 3, 1)]

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.an_affine_basis()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is not implemented for unbounded polyhedra
        """
        if not self.is_compact():
            raise NotImplementedError("this function is not implemented for unbounded polyhedra")

        chain = self.a_maximal_chain()[1:]  # we exclude the empty face
        chain_indices = [face.ambient_V_indices() for face in chain]
        basis_indices = []

        # We use in the following that elements in ``chain_indices`` are sorted lists
        # of V-indices.
        # Thus for each two faces we can easily find the first vertex that differs.
        for dim, face in enumerate(chain_indices):
            if dim == 0:
                # Append the vertex.
                basis_indices.append(face[0])
                continue

            prev_face = chain_indices[dim-1]
            for i in range(len(prev_face)):
                if prev_face[i] != face[i]:
                    # We found a vertex that ``face`` has, but its facet does not.
                    basis_indices.append(face[i])
                    break
            else:  # no break
                # ``prev_face`` contains all the same vertices as ``face`` until now.
                # But ``face`` is guaranteed to contain one more vertex (at least).
                basis_indices.append(face[len(prev_face)])

        return [self.Vrepresentation()[i] for i in basis_indices]

    def _test_an_affine_basis(self, tester=None, **options):
        """
        Run tests on the method :meth:`.an_affine_basis`

        TESTS::

            sage: polytopes.cross_polytope(3)._test_an_affine_basis()
        """
        if tester is None:
            tester = self._tester(**options)
        if self.is_compact():
            b = self.an_affine_basis()
            m = matrix([1] + list(v) for v in b)
            tester.assertEqual(m.rank(), self.dim() + 1)
            for v in b:
                tester.assertIn(v, self.vertices())

    def ray_generator(self):
        """
        Return a generator for the rays of the polyhedron.

        EXAMPLES::

            sage: pi = Polyhedron(ieqs = [[1,1,0],[1,0,1]])
            sage: pir = pi.ray_generator()
            sage: [x.vector() for x in pir]
            [(1, 0), (0, 1)]
        """
        for V in self.Vrepresentation():
            if V.is_ray():
                yield V

    @cached_method
    def rays(self):
        """
        Return a list of rays of the polyhedron.

        OUTPUT:

        A tuple of rays.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: p.rays()
            (A ray in the direction (1, 0, 0),
             A ray in the direction (0, 1, 0),
             A ray in the direction (0, 0, 1))
        """
        return tuple(self.ray_generator())

    def rays_list(self):
        """
        Return a list of rays as coefficient lists.

        .. NOTE::

            It is recommended to use :meth:`rays` or
            :meth:`ray_generator` instead to iterate over the list of
            :class:`Ray` objects.

        OUTPUT:

        A list of rays as lists of coordinates.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: p.rays_list()
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            sage: p.rays_list() == [list(r) for r in p.ray_generator()]
            True
        """
        return [list(x) for x in self.ray_generator()]

    def line_generator(self):
        """
        Return a generator for the lines of the polyhedron.

        EXAMPLES::

            sage: pr = Polyhedron(rays = [[1,0],[-1,0],[0,1]], vertices = [[-1,-1]])
            sage: next(pr.line_generator()).vector()
            (1, 0)
        """
        for V in self.Vrepresentation():
            if V.is_line():
                yield V

    @cached_method
    def lines(self):
        """
        Return all lines of the polyhedron.

        OUTPUT:

        A tuple of lines.

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.lines()
            (A line in the direction (1, 0),)
        """
        return tuple(self.line_generator())

    def lines_list(self):
        """
        Return a list of lines of the polyhedron.  The line data is given
        as a list of coordinates rather than as a Hrepresentation object.

        .. NOTE::

            It is recommended to use :meth:`line_generator` instead to
            iterate over the list of :class:`Line` objects.

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.lines_list()
            [[1, 0]]
            sage: p.lines_list() == [list(x) for x in p.line_generator()]
            True
        """
        return [list(x) for x in self.line_generator()]

    def bounded_edges(self):
        """
        Return the bounded edges (excluding rays and lines).

        OUTPUT:

        A generator for pairs of vertices, one pair per edge.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
            sage: [ e for e in p.bounded_edges() ]
            [(A vertex at (0, 1), A vertex at (1, 0))]
            sage: for e in p.bounded_edges(): print(e)
            (A vertex at (0, 1), A vertex at (1, 0))
        """
        obj = self.Vrepresentation()
        for i in range(len(obj)):
            if not obj[i].is_vertex():
                continue
            for j in range(i+1, len(obj)):
                if not obj[j].is_vertex():
                    continue
                if self.vertex_adjacency_matrix()[i, j] == 0:
                    continue
                yield (obj[i], obj[j])

    def Vrepresentation_space(self):
        r"""
        Return the ambient free module.

        OUTPUT:

        A free module over the base ring of dimension :meth:`ambient_dim`.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.Vrepresentation_space()
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
            sage: poly_test.ambient_space() is poly_test.Vrepresentation_space()
            True
        """
        return self.parent().Vrepresentation_space()

    ambient_space = Vrepresentation_space

    def ambient_vector_space(self, base_field=None):
        r"""
        Return the ambient vector space.

        It is the ambient free module (:meth:`Vrepresentation_space`) tensored
        with a field.

        INPUT:

        - ``base_field`` -- (default: the fraction field of the base ring) a field.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_vector_space()
            Vector space of dimension 4 over Rational Field
            sage: poly_test.ambient_vector_space() is poly_test.ambient()
            True

            sage: poly_test.ambient_vector_space(AA)
            Vector space of dimension 4 over Algebraic Real Field
            sage: poly_test.ambient_vector_space(RR)
            Vector space of dimension 4 over Real Field with 53 bits of precision
            sage: poly_test.ambient_vector_space(SR)
            Vector space of dimension 4 over Symbolic Ring
        """
        return self.Vrepresentation_space().vector_space(base_field=base_field)

    ambient = ambient_vector_space

    def Hrepresentation_space(self):
        r"""
        Return the linear space containing the H-representation vectors.

        OUTPUT:

        A free module over the base ring of dimension :meth:`ambient_dim` + 1.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.Hrepresentation_space()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
        """
        return self.parent().Hrepresentation_space()

    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_dim()
            4
        """
        return self.parent().ambient_dim()

    def dim(self):
        """
        Return the dimension of the polyhedron.

        OUTPUT:

        -1 if the polyhedron is empty, otherwise a non-negative integer.

        EXAMPLES::

            sage: simplex = Polyhedron(vertices = [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])
            sage: simplex.dim()
            3
            sage: simplex.ambient_dim()
            4

        The empty set is a special case (:trac:`12193`)::

            sage: P1=Polyhedron(vertices=[[1,0,0],[0,1,0],[0,0,1]])
            sage: P2=Polyhedron(vertices=[[2,0,0],[0,2,0],[0,0,2]])
            sage: P12 = P1.intersection(P2)
            sage: P12
            The empty polyhedron in ZZ^3
            sage: P12.dim()
            -1
        """
        if self.n_Vrepresentation() == 0:
            return -1   # the empty set
        else:
            return self.ambient_dim() - self.n_equations()

    dimension = dim

    def is_empty(self):
        """
        Test whether the polyhedron is the empty polyhedron

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,0,0],[0,1,0],[0,0,1]]);  P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: P.is_empty(), P.is_universe()
            (False, False)

            sage: Q = Polyhedron(vertices=());  Q
            The empty polyhedron in ZZ^0
            sage: Q.is_empty(), Q.is_universe()
            (True, False)

            sage: R = Polyhedron(lines=[(1,0),(0,1)]);  R
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
            sage: R.is_empty(), R.is_universe()
            (False, True)
        """
        return self.n_Vrepresentation() == 0

    def is_universe(self):
        """
        Test whether the polyhedron is the whole ambient space

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,0,0],[0,1,0],[0,0,1]]);  P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: P.is_empty(), P.is_universe()
            (False, False)

            sage: Q = Polyhedron(vertices=());  Q
            The empty polyhedron in ZZ^0
            sage: Q.is_empty(), Q.is_universe()
            (True, False)

            sage: R = Polyhedron(lines=[(1,0),(0,1)]);  R
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
            sage: R.is_empty(), R.is_universe()
            (False, True)
        """
        return self.n_Hrepresentation() == 0

    @cached_method
    def vertex_adjacency_matrix(self):
        """
        Return the binary matrix of vertex adjacencies.

        EXAMPLES::

            sage: polytopes.simplex(4).vertex_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]

        The rows and columns of the vertex adjacency matrix correspond
        to the :meth:`Vrepresentation` objects: vertices, rays, and
        lines. The `(i,j)` matrix entry equals `1` if the `i`-th and
        `j`-th V-representation object are adjacent.

        Two vertices are adjacent if they are the endpoints of an
        edge, that is, a one-dimensional face. For unbounded polyhedra
        this clearly needs to be generalized and we define two
        V-representation objects (see
        :mod:`sage.geometry.polyhedron.constructor`) to be adjacent if
        they together generate a one-face. There are three possible
        combinations:

        * Two vertices can bound a finite-length edge.

        * A vertex and a ray can generate a half-infinite edge
          starting at the vertex and with the direction given by the
          ray.

        * A vertex and a line can generate an infinite edge. The
          position of the vertex on the line is arbitrary in this
          case, only its transverse position matters. The direction of
          the edge is given by the line generator.

        For example, take the half-plane::

            sage: half_plane = Polyhedron(ieqs=[(0,1,0)])
            sage: half_plane.Hrepresentation()
            (An inequality (1, 0) x + 0 >= 0,)

        Its (non-unique) V-representation consists of a vertex, a ray,
        and a line. The only edge is spanned by the vertex and the
        line generator, so they are adjacent::

            sage: half_plane.Vrepresentation()
            (A line in the direction (0, 1), A ray in the direction (1, 0), A vertex at (0, 0))
            sage: half_plane.vertex_adjacency_matrix()
            [0 0 1]
            [0 0 0]
            [1 0 0]

        In one dimension higher, that is for a half-space in 3
        dimensions, there is no one-dimensional face. Hence nothing is
        adjacent::

            sage: Polyhedron(ieqs=[(0,1,0,0)]).vertex_adjacency_matrix()
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

        EXAMPLES:

        In a bounded polygon, every vertex has precisely two adjacent ones::

            sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (3, 0), (4, 1)])
            sage: for v in P.Vrep_generator():
            ....:     print("{} {}".format(P.adjacency_matrix().row(v.index()), v))
            (0, 1, 0, 1) A vertex at (0, 1)
            (1, 0, 1, 0) A vertex at (1, 0)
            (0, 1, 0, 1) A vertex at (3, 0)
            (1, 0, 1, 0) A vertex at (4, 1)

        If the V-representation of the polygon contains vertices and
        one ray, then each V-representation object is adjacent to two
        V-representation objects::

            sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (3, 0), (4, 1)],
            ....:                rays=[(0,1)])
            sage: for v in P.Vrep_generator():
            ....:       print("{} {}".format(P.adjacency_matrix().row(v.index()), v))
            (0, 1, 0, 0, 1) A ray in the direction (0, 1)
            (1, 0, 1, 0, 0) A vertex at (0, 1)
            (0, 1, 0, 1, 0) A vertex at (1, 0)
            (0, 0, 1, 0, 1) A vertex at (3, 0)
            (1, 0, 0, 1, 0) A vertex at (4, 1)

        If the V-representation of the polygon contains vertices and
        two distinct rays, then each vertex is adjacent to two
        V-representation objects (which can now be vertices or
        rays). The two rays are not adjacent to each other::

            sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (3, 0), (4, 1)],
            ....:                rays=[(0,1), (1,1)])
            sage: for v in P.Vrep_generator():
            ....:     print("{} {}".format(P.adjacency_matrix().row(v.index()), v))
            (0, 1, 0, 0, 0) A ray in the direction (0, 1)
            (1, 0, 1, 0, 0) A vertex at (0, 1)
            (0, 1, 0, 0, 1) A vertex at (1, 0)
            (0, 0, 0, 0, 1) A ray in the direction (1, 1)
            (0, 0, 1, 1, 0) A vertex at (3, 0)

        The vertex adjacency matrix has base ring integers. This way one can express various
        counting questions::

            sage: P = polytopes.cube()
            sage: Q = P.stack(P.faces(2)[0])
            sage: M = Q.vertex_adjacency_matrix()
            sage: sum(M)
            (4, 4, 3, 3, 4, 4, 4, 3, 3)
            sage: G = Q.vertex_graph()
            sage: G.degree()
            [4, 4, 3, 3, 4, 4, 4, 3, 3]

        TESTS:

        Check that :trac:`28828` is fixed::

                sage: P.adjacency_matrix().is_immutable()
                True
        """
        return self._vertex_adjacency_matrix()

    adjacency_matrix = vertex_adjacency_matrix

    def boundary_complex(self):
        """
        Return the simplicial complex given by the boundary faces of ``self``,
        if it is simplicial.

        OUTPUT:

        A (spherical) simplicial complex

        EXAMPLES:

        The boundary complex of the octahedron::

            sage: oc = polytopes.octahedron()
            sage: sc_oc = oc.boundary_complex()
            sage: fl_oc = oc.face_lattice()
            sage: fl_sc = sc_oc.face_poset()
            sage: [len(x) for x in fl_oc.level_sets()]
            [1, 6, 12, 8, 1]
            sage: [len(x) for x in fl_sc.level_sets()]
            [6, 12, 8]
            sage: sc_oc.euler_characteristic()
            2
            sage: sc_oc.homology()
            {0: 0, 1: 0, 2: Z}

        The polyhedron should be simplicial::

            sage: c = polytopes.cube()
            sage: c.boundary_complex()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is only implemented for simplicial polytopes

        TESTS::

            sage: p = Polyhedron(rays=[[1,1]])
            sage: p.boundary_complex()
            Traceback (most recent call last):
            ...
            ValueError: self should be compact
        """
        from sage.topology.simplicial_complex import SimplicialComplex
        if not self.is_compact():
            raise ValueError("self should be compact")

        if self.is_simplicial():
            inc_mat_cols = self.incidence_matrix().columns()
            ineq_indices = [inc_mat_cols[i].nonzero_positions()
                            for i in range(self.n_Hrepresentation())
                            if self.Hrepresentation()[i].is_inequality()]
            return SimplicialComplex(ineq_indices, maximality_check=False)
        else:
            raise NotImplementedError("this function is only implemented for simplicial polytopes")

    @cached_method
    def facet_adjacency_matrix(self):
        """
        Return the adjacency matrix for the facets and hyperplanes.

        EXAMPLES::

            sage: s4 = polytopes.simplex(4, project=True)
            sage: s4.facet_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]

        The facet adjacency matrix has base ring integers. This way one can express various
        counting questions::

            sage: P = polytopes.cube()
            sage: Q = P.stack(P.faces(2)[0])
            sage: M = Q.facet_adjacency_matrix()
            sage: sum(M)
            (4, 4, 4, 4, 3, 3, 3, 3, 4)

        TESTS:

        Check that :trac:`28828` is fixed::

                sage: s4.facet_adjacency_matrix().is_immutable()
                True
        """
        return self._facet_adjacency_matrix()

    @cached_method
    def incidence_matrix(self):
        """
        Return the incidence matrix.

        .. NOTE::

            The columns correspond to inequalities/equations in the
            order :meth:`Hrepresentation`, the rows correspond to
            vertices/rays/lines in the order
            :meth:`Vrepresentation`.

        .. SEEALSO::

            :meth:`slack_matrix`.

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: p.incidence_matrix()
            [0 0 1 1 0 1 0 0 0 0 1 0 0 0]
            [0 0 0 1 0 0 1 0 1 0 1 0 0 0]
            [0 0 1 1 1 0 0 1 0 0 0 0 0 0]
            [1 0 0 1 1 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 1 1 1 0 0 0]
            [0 0 1 0 0 1 0 1 0 0 0 1 0 0]
            [1 0 0 0 0 0 1 0 1 0 0 0 1 0]
            [1 0 0 0 1 0 0 1 0 0 0 0 0 1]
            [0 1 0 0 0 1 0 0 0 1 0 1 0 0]
            [0 1 0 0 0 0 0 0 1 1 0 0 1 0]
            [0 1 0 0 0 0 0 1 0 0 0 1 0 1]
            [1 1 0 0 0 0 0 0 0 0 0 0 1 1]
            sage: v = p.Vrepresentation(0)
            sage: v
            A vertex at (-1, -1, 0)
            sage: h = p.Hrepresentation(2)
            sage: h
            An inequality (1, 1, -1) x + 2 >= 0
            sage: h.eval(v)        # evaluation (1, 1, -1) * (-1/2, -1/2, 0) + 1
            0
            sage: h*v              # same as h.eval(v)
            0
            sage: p.incidence_matrix() [0,2]   # this entry is (v,h)
            1
            sage: h.contains(v)
            True
            sage: p.incidence_matrix() [2,0]   # note: not symmetric
            0

        The incidence matrix depends on the ambient dimension::

            sage: simplex = polytopes.simplex(); simplex
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
            sage: simplex.incidence_matrix()
            [1 1 1 1 0]
            [1 1 1 0 1]
            [1 1 0 1 1]
            [1 0 1 1 1]
            sage: simplex = simplex.affine_hull_projection(); simplex
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: simplex.incidence_matrix()
            [1 1 1 0]
            [1 1 0 1]
            [1 0 1 1]
            [0 1 1 1]

        An incidence matrix does not determine a unique
        polyhedron::

            sage: P = Polyhedron(vertices=[[0,1],[1,1],[1,0]])
            sage: P.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 1]

            sage: Q = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,1]])
            sage: Q.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 1]


        An example of two polyhedra with isomorphic face lattices
        but different incidence matrices::

            sage: Q.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 1]

            sage: R = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,3/2], [3/2,1]])
            sage: R.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 0]
            [0 0 1]

        The incidence matrix has base ring integers. This way one can express various
        counting questions::

            sage: P = polytopes.twenty_four_cell()
            sage: M = P.incidence_matrix()
            sage: sum(sum(x) for x in M) == P.flag_f_vector(0,3)
            True

        TESTS:

        Check that :trac:`28828` is fixed::

            sage: R.incidence_matrix().is_immutable()
            True

        Test that this method works for inexact base ring
        (``cdd`` sets the cache already)::

            sage: P = polytopes.dodecahedron(exact=False)
            sage: M = P.incidence_matrix.cache
            sage: P.incidence_matrix.clear_cache()
            sage: M == P.incidence_matrix()
            True
        """
        if self.base_ring() in (ZZ, QQ):
            # Much faster for integers or rationals.
            incidence_matrix = self.slack_matrix().zero_pattern_matrix(ZZ)
            incidence_matrix.set_immutable()
            return incidence_matrix

        incidence_matrix = matrix(ZZ, self.n_Vrepresentation(),
                                  self.n_Hrepresentation(), 0)

        Vvectors_vertices = tuple((v.vector(), v.index())
                                  for v in self.Vrep_generator()
                                  if v.is_vertex())
        Vvectors_rays_lines = tuple((v.vector(), v.index())
                                    for v in self.Vrep_generator()
                                    if not v.is_vertex())

        # Determine ``is_zero`` to save lots of time.
        if self.base_ring().is_exact():
            def is_zero(x):
                return not x
        else:
            is_zero = self._is_zero

        for H in self.Hrep_generator():
            Hconst = H.b()
            Hvec = H.A()
            Hindex = H.index()
            for Vvec, Vindex in Vvectors_vertices:
                if is_zero(Hvec*Vvec + Hconst):
                    incidence_matrix[Vindex, Hindex] = 1

            # A ray or line is considered incident with a hyperplane,
            # if it is orthogonal to the normal vector of the hyperplane.
            for Vvec, Vindex in Vvectors_rays_lines:
                if is_zero(Hvec*Vvec):
                    incidence_matrix[Vindex, Hindex] = 1

        incidence_matrix.set_immutable()
        return incidence_matrix

    @cached_method
    def slack_matrix(self):
        r"""
        Return the slack matrix.

        The entries correspond to the evaluation of the Hrepresentation
        elements on the  Vrepresentation elements.

        .. NOTE::

            The columns correspond to inequalities/equations in the
            order :meth:`Hrepresentation`, the rows correspond to
            vertices/rays/lines in the order
            :meth:`Vrepresentation`.

        .. SEEALSO::

            :meth:`incidence_matrix`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.slack_matrix()
            [0 2 2 2 0 0]
            [0 0 2 2 0 2]
            [0 0 0 2 2 2]
            [0 2 0 2 2 0]
            [2 2 0 0 2 0]
            [2 2 2 0 0 0]
            [2 0 2 0 0 2]
            [2 0 0 0 2 2]

            sage: P = polytopes.cube(intervals='zero_one')
            sage: P.slack_matrix()
            [0 1 1 1 0 0]
            [0 0 1 1 0 1]
            [0 0 0 1 1 1]
            [0 1 0 1 1 0]
            [1 1 0 0 1 0]
            [1 1 1 0 0 0]
            [1 0 1 0 0 1]
            [1 0 0 0 1 1]

            sage: P = polytopes.dodecahedron().faces(2)[0].as_polyhedron()
            sage: P.slack_matrix()
            [1/2*sqrt5 - 1/2               0               0               1 1/2*sqrt5 - 1/2               0]
            [              0               0 1/2*sqrt5 - 1/2 1/2*sqrt5 - 1/2               1               0]
            [              0 1/2*sqrt5 - 1/2               1               0 1/2*sqrt5 - 1/2               0]
            [              1 1/2*sqrt5 - 1/2               0 1/2*sqrt5 - 1/2               0               0]
            [1/2*sqrt5 - 1/2               1 1/2*sqrt5 - 1/2               0               0               0]

            sage: P = Polyhedron(rays=[[1, 0], [0, 1]])
            sage: P.slack_matrix()
            [0 0]
            [0 1]
            [1 0]

        TESTS::

            sage: Polyhedron().slack_matrix()
            []
            sage: Polyhedron(base_ring=QuadraticField(2)).slack_matrix().base_ring()
            Number Field in a with defining polynomial x^2 - 2 with a = 1.41...
        """
        if not self.n_Vrepresentation() or not self.n_Hrepresentation():
            slack_matrix = matrix(self.base_ring(), self.n_Vrepresentation(),
                                  self.n_Hrepresentation(), 0)
        else:
            Vrep_matrix = matrix(self.base_ring(), self.Vrepresentation())
            Hrep_matrix = matrix(self.base_ring(), self.Hrepresentation())

            # Getting homogeneous coordinates of the Vrepresentation.
            hom_helper = matrix(self.base_ring(), [1 if v.is_vertex() else 0 for v in self.Vrepresentation()])
            hom_Vrep = hom_helper.stack(Vrep_matrix.transpose())

            slack_matrix = (Hrep_matrix * hom_Vrep).transpose()

        slack_matrix.set_immutable()
        return slack_matrix

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        The ring over which the polyhedron is defined. Must be a
        sub-ring of the reals to define a polyhedron, in particular
        comparison must be defined. Popular choices are

        * ``ZZ`` (the ring of integers, lattice polytope),

        * ``QQ`` (exact arithmetic using gmp),

        * ``RDF`` (double precision floating-point arithmetic), or

        * ``AA`` (real algebraic field).

        EXAMPLES::

            sage: triangle = Polyhedron(vertices = [[1,0],[0,1],[1,1]])
            sage: triangle.base_ring() == ZZ
            True
        """
        return self.parent().base_ring()

    def backend(self):
        """
        Return the backend used.

        OUTPUT:

        The name of the backend used for computations. It will be one of
        the following backends:

         * ``ppl`` the Parma Polyhedra Library

         * ``cdd`` CDD

         * ``normaliz`` normaliz

         * ``polymake`` polymake

         * ``field`` a generic Sage implementation

        EXAMPLES::

            sage: triangle = Polyhedron(vertices = [[1, 0], [0, 1], [1, 1]])
            sage: triangle.backend()
            'ppl'
            sage: D = polytopes.dodecahedron()
            sage: D.backend()
            'field'
            sage: P = Polyhedron([[1.23]])
            sage: P.backend()
            'cdd'
        """
        return self.parent().backend()

    @cached_method
    def center(self):
        """
        Return the average of the vertices.

        .. SEEALSO::

            :meth:`representative_point`.

        OUTPUT:

        The center of the polyhedron. All rays and lines are
        ignored. Raises a ``ZeroDivisionError`` for the empty
        polytope.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p = p + vector([1,0,0])
            sage: p.center()
            (1, 0, 0)
        """
        if self.dim() == 0:
            return self.vertices()[0].vector()
        else:
            vertex_sum = vector(self.base_ring(), [0]*self.ambient_dim())
            for v in self.vertex_generator():
                vertex_sum += v.vector()
            vertex_sum.set_immutable()
            return vertex_sum / self.n_vertices()

    @cached_method(do_pickle=True)
    def centroid(self, engine='auto', **kwds):
        r"""
        Return the center of the mass of the polytope.

        The mass is taken with respect to the induced Lebesgue measure,
        see :meth:`volume`.

        If the polyhedron is not compact, a ``NotImplementedError`` is
        raised.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal',
          'TOPCOM', or 'normaliz'.  The 'internal' and 'TOPCOM' instruct
          this package to always use its own triangulation algorithms
          or TOPCOM's algorithms, respectively. By default ('auto'),
          TOPCOM is used if it is available and internal routines otherwise.

        - ``**kwds`` -- keyword arguments that are passed to the
          triangulation engine (see :meth:`triangulate`).

        OUTPUT: The centroid as vector.

        ALGORITHM:

        We triangulate the polytope and find the barycenter of the simplices.
        We add the individual barycenters weighted by the fraction of the total
        mass.

        EXAMPLES::

            sage: P = polytopes.hypercube(2).pyramid()
            sage: P.centroid()
            (1/4, 0, 0)

            sage: P = polytopes.associahedron(['A',2])
            sage: P.centroid()
            (2/21, 2/21)

            sage: P = polytopes.permutahedron(4, backend='normaliz')  # optional - pynormaliz
            sage: P.centroid()                                        # optional - pynormaliz
            (5/2, 5/2, 5/2, 5/2)

        The method is not implemented for unbounded polyhedra::

            sage: P = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: P.centroid()
            Traceback (most recent call last):
            ...
            NotImplementedError: the polyhedron is not compact

        The centroid of an empty polyhedron is not defined::

            sage: Polyhedron().centroid()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        TESTS::

            sage: Polyhedron(vertices=[[0,1]]).centroid()
            (0, 1)
        """
        if not self.is_compact():
            raise NotImplementedError("the polyhedron is not compact")
        if self.n_vertices() == self.dim() + 1:
            # The centroid of a simplex is its center.
            return self.center()

        triangulation = self.triangulate(engine=engine, **kwds)

        if self.ambient_dim() == self.dim():
            pc = triangulation.point_configuration()
        else:
            from sage.geometry.triangulation.point_configuration import PointConfiguration
            A, b = self.affine_hull_projection(as_affine_map=True, orthogonal=True, orthonormal=True, extend=True)
            pc = PointConfiguration((A(v.vector()) for v in self.Vrep_generator()))

        barycenters = [sum(self.Vrepresentation(i).vector() for i in simplex)/(self.dim() + 1) for simplex in triangulation]
        volumes = [pc.volume(simplex) for simplex in triangulation]

        centroid = sum(volumes[i]*barycenters[i] for i in range(len(volumes)))/sum(volumes)
        if self.ambient_dim() != self.dim():
            # By the affine hull projection, the centroid has base ring ``AA``,
            # we try return the centroid in a reasonable ring.
            try:
                return centroid.change_ring(self.base_ring().fraction_field())
            except ValueError:
                pass
        return centroid

    @cached_method
    def representative_point(self):
        """
        Return a "generic" point.

        .. SEEALSO::

            :meth:`center`.

        OUTPUT:

        A point as a coordinate vector. The point is chosen to be
        interior as far as possible. If the polyhedron is not
        full-dimensional, the point is in the relative interior. If
        the polyhedron is zero-dimensional, its single point is
        returned.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(3,2)], rays=[(1,-1)])
            sage: p.representative_point()
            (4, 1)
            sage: p.center()
            (3, 2)

            sage: Polyhedron(vertices=[(3,2)]).representative_point()
            (3, 2)
        """
        accumulator = vector(self.base_ring(), [0]*self.ambient_dim())
        for v in self.vertex_generator():
            accumulator += v.vector()
        accumulator /= self.n_vertices()
        for r in self.ray_generator():
            accumulator += r.vector()
        accumulator.set_immutable()
        return accumulator

    def _some_elements_(self):
        r"""
        Generate some points of ``self``.

        If ``self`` is empty, no points are generated; no exception will be raised.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: P.an_element()            # indirect doctest
            (1/4, 1/4, 1/4, 1/4)
            sage: P.some_elements()         # indirect doctest
            [(1/4, 1/4, 1/4, 1/4),
             (0, 0, 0, 1),
             (0, 0, 1/2, 1/2),
             (0, 1/2, 1/4, 1/4),
             (1/2, 1/4, 1/8, 1/8)]
        """
        if self.is_empty():
            return
        yield self.representative_point()
        vertex_iter = iter(self.vertex_generator())
        try:
            p = next(vertex_iter).vector()
            yield vector(p, immutable=True)
            for i in range(4):
                p = (p + next(vertex_iter).vector()) / 2
                yield vector(p, immutable=True)
        except StopIteration:
            pass

    def a_maximal_chain(self):
        r"""
        Return a maximal chain of the face lattice in increasing order.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.a_maximal_chain()
            [A -1-dimensional face of a Polyhedron in ZZ^3,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices]
            sage: P = polytopes.cube()
            sage: chain = P.a_maximal_chain(); chain
            [A -1-dimensional face of a Polyhedron in ZZ^3,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices]
            sage: [face.ambient_V_indices() for face in chain]
            [(), (5,), (0, 5), (0, 3, 4, 5), (0, 1, 2, 3, 4, 5, 6, 7)]

        TESTS::

        Check output for the empty polyhedron::

            sage: P = Polyhedron()
            sage: P.a_maximal_chain()
            [A -1-dimensional face of a Polyhedron in ZZ^0]
        """
        comb_chain = self.combinatorial_polyhedron().a_maximal_chain()

        from sage.geometry.polyhedron.face import combinatorial_face_to_polyhedral_face
        empty_face = self.faces(-1)[0]
        universe = self.faces(self.dim())[0]

        if self.dim() == -1:
            return [empty_face]

        return [empty_face] + \
               [combinatorial_face_to_polyhedral_face(self, face)
                for face in comb_chain] + \
               [universe]

    @cached_method
    def radius_square(self):
        """
        Return the square of the maximal distance from the
        :meth:`center` to a vertex. All rays and lines are ignored.

        OUTPUT:

        The square of the radius, which is in :meth:`base_ring`.

        EXAMPLES::

            sage: p = polytopes.permutahedron(4, project = False)
            sage: p.radius_square()
            5
        """
        vertices = [v.vector() - self.center() for v in self.vertex_generator()]
        return max(v.dot_product(v) for v in vertices)

    def radius(self):
        """
        Return the maximal distance from the center to a vertex. All
        rays and lines are ignored.

        OUTPUT:

        The radius for a rational polyhedron is, in general, not
        rational.  use :meth:`radius_square` if you need a rational
        distance measure.

        EXAMPLES::

            sage: p = polytopes.hypercube(4)
            sage: p.radius()
            2
        """
        return sqrt(self.radius_square())

    def is_inscribed(self, certificate=False):
        """
        This function tests whether the vertices of the polyhedron are
        inscribed on a sphere.

        The polyhedron is expected to be compact and full-dimensional.
        A full-dimensional compact polytope is inscribed if there exists
        a point in space which is equidistant to all its vertices.

        ALGORITHM:

        The function first computes the circumsphere of a full-dimensional
        simplex with vertices of ``self``. It is found by lifting the points on a
        paraboloid to find the hyperplane on which the circumsphere is lifted.
        Then, it checks if all other vertices are equidistant to the
        circumcenter of that simplex.

        INPUT:

        - ``certificate`` -- (default: ``False``) boolean; specifies whether to
          return the circumcenter, if found.

        OUTPUT:

        If ``certificate`` is true, returns a tuple containing:

        1. Boolean.
        2. The circumcenter of the polytope or None.

        If ``certificate`` is false:

        - a Boolean.

        EXAMPLES::

            sage: q = Polyhedron(vertices = [[1,1,1,1],[-1,-1,1,1],[1,-1,-1,1],
            ....:                            [-1,1,-1,1],[1,1,1,-1],[-1,-1,1,-1],
            ....:                            [1,-1,-1,-1],[-1,1,-1,-1],[0,0,10/13,-24/13],
            ....:                            [0,0,-10/13,-24/13]])
            sage: q.is_inscribed(certificate=True)
            (True, (0, 0, 0, 0))

            sage: cube = polytopes.cube()
            sage: cube.is_inscribed()
            True

            sage: translated_cube = Polyhedron(vertices=[v.vector() + vector([1,2,3])
            ....:                                        for v in cube.vertices()])
            sage: translated_cube.is_inscribed(certificate=True)
            (True, (1, 2, 3))

            sage: truncated_cube = cube.face_truncation(cube.faces(0)[0])
            sage: truncated_cube.is_inscribed()
            False

        The method is not implemented for non-full-dimensional polytope or
        unbounded polyhedra::

            sage: square = Polyhedron(vertices=[[1,0,0],[0,1,0],[1,1,0],[0,0,0]])
            sage: square.is_inscribed()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for full-dimensional polyhedra only

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.is_inscribed()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is not implemented for unbounded polyhedra

        TESTS:

        We check that :trac:`28464` is fixed::

            sage: P = Polyhedron(vertices=[(-130658298093891402635075/416049251842505144482473,
            ....: 177469511761879509172000/1248147755527515433447419,
            ....: 485550543257132133136169/2496295511055030866894838,
            ....: 2010744967797898733758669/2496295511055030866894838),
            ....: (-146945725603929909850/706333405676769433081,
            ....: -84939725782618445000/706333405676769433081,
            ....: 560600045283000988081/1412666811353538866162,
            ....: 969778382942371268081/1412666811353538866162),
            ....: (-46275018824497300/140422338198040641,
            ....: -5747688262110000/46807446066013547, 1939357556329/7033601552658,
            ....: 1939357556329/7033601552658), (-17300/59929, -10000/59929, 39929/119858,
            ....: 39929/119858), (-4700/32209, -10000/32209, 12209/64418, 12209/64418),
            ....: (QQ(0), QQ(0), QQ(0), QQ(1)), (QQ(0), QQ(0), 1/2, 1/2), (300/10027,
            ....: -10000/30081, 10081/60162, 10081/60162), (112393975400/1900567733649,
            ....: 117311600000/633522577883, 43678681/95197362, 43678681/95197362),
            ....: (6109749955400/133380598418321, 37106807920000/133380598418321,
            ....: 2677964249/6680888498, 2677964249/6680888498),
            ....: (29197890764005600/402876806828660641,
            ....: -2150510776960000/402876806828660641,
            ....: 398575785274740641/805753613657321282,
            ....: 398575785274740641/805753613657321282),
            ....: (5576946899441759759983005325/110078073300232813237456943251,
            ....: -29071211718677797926570478000/110078073300232813237456943251,
            ....: 59439312069347378584317232001/220156146600465626474913886502,
            ....: 181346577228466312205473034501/220156146600465626474913886502),
            ....: (150040732779124914266530235300/6774574358246204311268446913881,
            ....: -2813827375989039189507000218000/6774574358246204311268446913881,
            ....: 1260217414021285074925933133881/13549148716492408622536893827762,
            ....: 3232518047094242684574253773881/13549148716492408622536893827762),
            ....: (3816349407976279597850158016285000/88842127448735433741180809504357161,
            ....: 27965821247423216557301387453968000/88842127448735433741180809504357161,
            ....: 68546256000224819256028677086357161/177684254897470867482361619008714322,
            ....: 86062257922545755787315412690197161/177684254897470867482361619008714322)])
            sage: P.is_inscribed()
            True

            sage: P = Polyhedron(vertices=[[0, -1, 0, 0],
            ....:                          [0, 0, -1, 0],
            ....:                          [0, 0, 0, -1],
            ....:                          [0, 0, +1, 0],
            ....:                          [0, 0, 0, +1],
            ....:                          [+1, 0, 0, 0]])
            sage: P.is_inscribed()
            True

        We check that :trac:`29125` is fixed::

            sage: P = Polyhedron(vertices=[[-2,-1], [-2,1], [0,-1], [0,1]], backend='field')
            sage: P.is_inscribed()
            True
            sage: V = P.Vrepresentation()
            sage: H = P.Hrepresentation()
            sage: parent = P.parent()
            sage: for V1 in Permutations(V):
            ....:     P1 = parent._element_constructor_(
            ....:         [V1, [], []], [H, []], Vrep_minimal=True, Hrep_minimal=True)
            ....:     assert P1.is_inscribed()
        """

        if not self.is_compact():
            raise NotImplementedError("this function is not implemented for unbounded polyhedra")

        if not self.is_full_dimensional():
            raise NotImplementedError("this function is implemented for full-dimensional polyhedra only")

        dimension = self.dimension()
        vertices = self.vertices()

        # We obtain vertices that are an affine basis of the affine hull.
        affine_basis = self.an_affine_basis()
        raw_data = []
        for vertex in affine_basis:
            vertex_vector = vertex.vector()
            raw_data += [[sum(i**2 for i in vertex_vector)] +
                         [i for i in vertex_vector] + [1]]
        matrix_data = matrix(raw_data)

        # The determinant "a" should not be zero because
        # the vertices in ``affine_basis`` are an affine basis.
        a = matrix_data.matrix_from_columns(range(1, dimension+2)).determinant()

        minors = [(-1)**(i)*matrix_data.matrix_from_columns([j for j in range(dimension+2) if j != i]).determinant()
                  for i in range(1, dimension+1)]
        c = (-1)**(dimension+1)*matrix_data.matrix_from_columns(range(dimension+1)).determinant()

        circumcenter = vector([minors[i]/(2*a) for i in range(dimension)])
        squared_circumradius = (sum(m**2 for m in minors) - 4 * a * c) / (4*a**2)

        # Checking if the circumcenter has the correct sign
        if not all(sum(i**2 for i in v.vector() - circumcenter) == squared_circumradius
                   for v in vertices if v in affine_basis):
            circumcenter = - circumcenter

        is_inscribed = all(sum(i**2 for i in v.vector() - circumcenter) == squared_circumradius
                           for v in vertices if v not in affine_basis)

        if certificate:
            if is_inscribed:
                return (True, circumcenter)
            else:
                return (False, None)
        else:
            return is_inscribed

    def is_compact(self):
        """
        Test for boundedness of the polytope.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: p.is_compact()
            True
            sage: p = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,-1,0,0]])
            sage: p.is_compact()
            False
        """
        return self.n_rays() == 0 and self.n_lines() == 0

    @cached_method
    def combinatorial_polyhedron(self):
        r"""
        Return the combinatorial type of ``self``.

        See :class:`sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

        EXAMPLES::

            sage: polytopes.cube().combinatorial_polyhedron()
            A 3-dimensional combinatorial polyhedron with 6 facets

            sage: polytopes.cyclic_polytope(4,10).combinatorial_polyhedron()
            A 4-dimensional combinatorial polyhedron with 35 facets

            sage: Polyhedron(rays=[[0,1], [1,0]]).combinatorial_polyhedron()
            A 2-dimensional combinatorial polyhedron with 2 facets
        """
        from sage.geometry.polyhedron.combinatorial_polyhedron.base import CombinatorialPolyhedron
        return CombinatorialPolyhedron(self)

    def _test_combinatorial_polyhedron(self, tester=None, **options):
        """
        Run test suite of combinatorial polyhedron.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_combinatorial_polyhedron()
        """
        from sage.misc.sage_unittest import TestSuite

        tester = self._tester(tester=tester, **options)
        tester.info("\n  Running the test suite of self.combinatorial_polyhedron()")
        TestSuite(self.combinatorial_polyhedron()).run(verbose=tester._verbose,
                                                       prefix=tester._prefix+"  ")
        tester.info(tester._prefix+" ", newline = False)

    def simplicity(self):
        r"""
        Return the largest integer `k` such that the polytope is `k`-simple.

        A polytope `P` is `k`-simple, if every `(d-1-k)`-face
        is contained in exactly `k+1` facets of `P` for `1 \leq k \leq d-1`.
        Equivalently it is `k`-simple if the polar/dual polytope is `k`-simplicial.
        If `self` is a simplex, it returns its dimension.

        EXAMPLES::

            sage: polytopes.hypersimplex(4,2).simplicity()
            1
            sage: polytopes.hypersimplex(5,2).simplicity()
            2
            sage: polytopes.hypersimplex(6,2).simplicity()
            3
            sage: polytopes.simplex(3).simplicity()
            3
            sage: polytopes.simplex(1).simplicity()
            1

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.simplicity()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only
        """
        if not(self.is_compact()):
            raise NotImplementedError("this function is implemented for polytopes only")
        return self.combinatorial_polyhedron().simplicity()

    def is_simple(self):
        """
        Test for simplicity of a polytope.

        See :wikipedia:`Simple_polytope`

        EXAMPLES::

            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p.is_simple()
            True
            sage: p = Polyhedron([[0,0,0],[4,4,0],[4,0,0],[0,4,0],[2,2,2]])
            sage: p.is_simple()
            False
        """
        if not self.is_compact():
            return False
        return self.combinatorial_polyhedron().is_simple()

    def simpliciality(self):
        r"""
        Return the largest integer `k` such that the polytope is `k`-simplicial.

        A polytope is `k`-simplicial, if every `k`-face is a simplex.
        If `self` is a simplex, returns its dimension.

        EXAMPLES::

            sage: polytopes.cyclic_polytope(10,4).simpliciality()
            3
            sage: polytopes.hypersimplex(5,2).simpliciality()
            2
            sage: polytopes.cross_polytope(4).simpliciality()
            3
            sage: polytopes.simplex(3).simpliciality()
            3
            sage: polytopes.simplex(1).simpliciality()
            1

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.simpliciality()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only
        """
        if not(self.is_compact()):
            raise NotImplementedError("this function is implemented for polytopes only")
        return self.combinatorial_polyhedron().simpliciality()

    def is_simplicial(self):
        """
        Tests if the polytope is simplicial

        A polytope is simplicial if every facet is a simplex.

        See :wikipedia:`Simplicial_polytope`

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p.is_simplicial()
            False
            sage: q = polytopes.simplex(5, project=True)
            sage: q.is_simplicial()
            True
            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p.is_simplicial()
            True
            sage: q = Polyhedron([[1,1,1],[-1,1,1],[1,-1,1],[-1,-1,1],[1,1,-1]])
            sage: q.is_simplicial()
            False
            sage: P = polytopes.simplex(); P
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
            sage: P.is_simplicial()
            True

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.is_simplicial()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only
        """
        if not(self.is_compact()):
            raise NotImplementedError("this function is implemented for polytopes only")
        return self.combinatorial_polyhedron().is_simplicial()

    def is_pyramid(self, certificate=False):
        """
        Test whether the polytope is a pyramid over one of its facets.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return a vertex of the polytope which is the apex of a pyramid,
          if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. The apex of the pyramid or ``None``.

        If ``certificate`` is ``False`` returns a boolean.

        EXAMPLES::

            sage: P = polytopes.simplex(3)
            sage: P.is_pyramid()
            True
            sage: P.is_pyramid(certificate=True)
            (True, A vertex at (1, 0, 0, 0))
            sage: egyptian_pyramid = polytopes.regular_polygon(4).pyramid()
            sage: egyptian_pyramid.is_pyramid()
            True
            sage: Q = polytopes.octahedron()
            sage: Q.is_pyramid()
            False

        For the `0`-dimensional polyhedron, the output is ``True``,
        but it cannot be constructed as a pyramid over the empty polyhedron::

            sage: P = Polyhedron([[0]])
            sage: P.is_pyramid()
            True
            sage: Polyhedron().pyramid()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        if not self.is_compact():
            raise ValueError("polyhedron has to be compact")

        return self.combinatorial_polyhedron().is_pyramid(certificate)

    def is_bipyramid(self, certificate=False):
        r"""
        Test whether the polytope is combinatorially equivalent to a
        bipyramid over some polytope.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return two vertices of the polytope which are the apices of a
          bipyramid, if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. ``None`` or a tuple containing:
            a. The first apex.
            b. The second apex.

        If ``certificate`` is ``False`` returns a boolean.

        EXAMPLES::

            sage: P = polytopes.octahedron()
            sage: P.is_bipyramid()
            True
            sage: P.is_bipyramid(certificate=True)
            (True, [A vertex at (-1, 0, 0), A vertex at (1, 0, 0)])
            sage: Q = polytopes.cyclic_polytope(3,7)
            sage: Q.is_bipyramid()
            False
            sage: R = Q.bipyramid()
            sage: R.is_bipyramid(certificate=True)
            (True, [A vertex at (-1, 3, 13, 63), A vertex at (1, 3, 13, 63)])

        TESTS::

            sage: P = polytopes.permutahedron(4).bipyramid()
            sage: P.is_bipyramid()
            True

            sage: P = polytopes.cube()
            sage: P.is_bipyramid()
            False

            sage: P = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,1]])
            sage: P.is_bipyramid()
            Traceback (most recent call last):
            ...
            ValueError: polyhedron has to be compact

        ALGORITHM:

        Assume all faces of a polyhedron to be given as lists of vertices.

        A polytope is a bipyramid with apexes `v`, `w` if and only if for each
        proper face `v \in F` there exists a face `G` with
        `G \setminus \{w\} = F \setminus \{v\}`
        and vice versa (for each proper face
        `w \in F` there exists ...).

        To check this property it suffices to check for all facets of the polyhedron.
        """
        if not self.is_compact():
            raise ValueError("polyhedron has to be compact")

        from sage.misc.functional import is_odd
        n_verts = self.n_vertices()
        n_facets = self.n_facets()
        if is_odd(n_facets):
            if certificate:
                return (False, None)
            return False

        IM = self.incidence_matrix()
        if self.n_equations():
            # Remove equations from the incidence matrix,
            # such that this is the vertex-facet incidences matrix.
            I1 = IM.transpose()
            I2 = I1[[i for i in range(self.n_Hrepresentation())
                     if not self.Hrepresentation()[i].is_equation()]]
            IM = I2.transpose()

        facets_incidences = [set(column.nonzero_positions()) for column in IM.columns()]
        verts_incidences = dict()
        for i in range(n_verts):
            v_i = set(IM.row(i).nonzero_positions())
            if len(v_i) == n_facets/2:
                verts_incidences[i] = v_i

        # Find two vertices ``vert1`` and ``vert2`` such that one of them
        # lies on exactly half of the facets, and the other one lies on
        # exactly the other half.
        from itertools import combinations
        for index1, index2 in combinations(verts_incidences, 2):
            vert1_incidences = verts_incidences[index1]
            vert2_incidences = verts_incidences[index2]
            vert1and2 = vert1_incidences.union(vert2_incidences)
            if len(vert1and2) == n_facets:
                # We have found two candidates for apexes.
                # Remove from each facet ``index1`` resp. ``index2``.
                test_facets = set(frozenset(facet_inc.difference({index1, index2}))
                                  for facet_inc in facets_incidences)
                if len(test_facets) == n_facets/2:
                    # For each `F` containing `index1` there is
                    # `G` containing `index2` such that
                    # `F \setminus \{index1\} =  G \setminus \{index2\}
                    # and vice versa.
                    if certificate:
                        V = self.vertices()
                        return (True, [V[index1], V[index2]])
                    return True

        if certificate:
            return (False, None)
        return False

    def is_prism(self, certificate=False):
        """
        Test whether the polytope is combinatorially equivalent to a prism of
        some polytope.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return two facets of the polytope which are the bases of a prism,
          if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. ``None`` or a tuple containing:
            a. List of the vertices of the first base facet.
            b. List of the vertices of the second base facet.

        If ``certificate`` is ``False`` returns a boolean.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.is_prism()
            True
            sage: P.is_prism(certificate=True)
            (True,
             [[A vertex at (1, -1, -1),
               A vertex at (1, 1, -1),
               A vertex at (1, 1, 1),
               A vertex at (1, -1, 1)],
              [A vertex at (-1, -1, 1),
               A vertex at (-1, -1, -1),
               A vertex at (-1, 1, -1),
               A vertex at (-1, 1, 1)]])
            sage: Q = polytopes.cyclic_polytope(3,8)
            sage: Q.is_prism()
            False
            sage: R = Q.prism()
            sage: R.is_prism(certificate=True)
            (True,
             [[A vertex at (1, 6, 36, 216),
               A vertex at (1, 0, 0, 0),
               A vertex at (1, 7, 49, 343),
               A vertex at (1, 5, 25, 125),
               A vertex at (1, 1, 1, 1),
               A vertex at (1, 2, 4, 8),
               A vertex at (1, 4, 16, 64),
               A vertex at (1, 3, 9, 27)],
              [A vertex at (0, 3, 9, 27),
               A vertex at (0, 6, 36, 216),
               A vertex at (0, 0, 0, 0),
               A vertex at (0, 7, 49, 343),
               A vertex at (0, 5, 25, 125),
               A vertex at (0, 1, 1, 1),
               A vertex at (0, 2, 4, 8),
               A vertex at (0, 4, 16, 64)]])

        TESTS::

            sage: P = polytopes.cross_polytope(5)
            sage: P.is_prism()
            False

            sage: P = polytopes.permutahedron(4).prism()
            sage: P.is_prism()
            True

            sage: P = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,1]])
            sage: P.is_prism()
            Traceback (most recent call last):
            ...
            NotImplementedError: polyhedron has to be compact

        ALGORITHM:

        See :meth:`Polyhedron_base.is_bipyramid`.
        """
        if not self.is_compact():
            raise NotImplementedError("polyhedron has to be compact")

        from sage.misc.functional import is_odd
        n_verts = self.n_vertices()
        n_facets = self.n_facets()
        if is_odd(n_verts):
            if certificate:
                return (False, None)
            return False

        IM = self.incidence_matrix()
        if self.n_equations():
            # Remove equations from the incidence matrix,
            # such that this is the vertex-facet incidences matrix.
            I1 = IM.transpose()
            I2 = I1[[i for i in range(self.n_Hrepresentation())
                     if not self.Hrepresentation()[i].is_equation()]]
            IM = I2.transpose()

        verts_incidences = [set(row.nonzero_positions()) for row in IM.rows()]
        facets_incidences = dict()
        for j in range(n_facets):
            F_j = set(IM.column(j).nonzero_positions())
            if len(F_j) == n_verts/2:
                facets_incidences[j] = F_j

        # Find two vertices ``facet1`` and ``facet2`` such that one of them
        # contains exactly half of the vertices, and the other one contains
        # exactly the other half.
        from itertools import combinations
        for index1, index2 in combinations(facets_incidences, 2):
            facet1_incidences = facets_incidences[index1]
            facet2_incidences = facets_incidences[index2]
            facet1and2 = facet1_incidences.union(facet2_incidences)
            if len(facet1and2) == n_verts:
                # We have found two candidates for base faces.
                # Remove from each vertex ``index1`` resp. ``index2``.
                test_verts = set(frozenset(vert_inc.difference({index1, index2}))
                                 for vert_inc in verts_incidences)
                if len(test_verts) == n_verts/2:
                    # For each vertex containing `index1` there is
                    # another one contained in `index2`
                    # and vice versa.
                    # Other than `index1` and `index2` both are contained in
                    # exactly the same facets.
                    if certificate:
                        V = self.vertices()
                        facet1_vertices = [V[i] for i in facet1_incidences]
                        facet2_vertices = [V[i] for i in facet2_incidences]
                        return (True, [facet1_vertices, facet2_vertices])
                    return True

        if certificate:
            return (False, None)
        return False

    def hyperplane_arrangement(self):
        """
        Return the hyperplane arrangement defined by the equations and
        inequalities.

        OUTPUT:

        A :class:`hyperplane arrangement
        <sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement>`
        consisting of the hyperplanes defined by the
        :meth:`Hrepresentation`.
        If the polytope is full-dimensional, this is the hyperplane
        arrangement spanned by the facets of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.hypercube(2)
            sage: p.hyperplane_arrangement()
            Arrangement <-t0 + 1 | -t1 + 1 | t1 + 1 | t0 + 1>
        """
        names = tuple('t' + str(i) for i in range(self.ambient_dim()))
        from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
        field = self.base_ring().fraction_field()
        H = HyperplaneArrangements(field, names)
        return H(self)

    @cached_method
    def gale_transform(self):
        """
        Return the Gale transform of a polytope as described in the
        reference below.

        OUTPUT:

        A list of vectors, the Gale transform.  The dimension is the
        dimension of the affine dependencies of the vertices of the
        polytope.

        EXAMPLES:

        This is from the reference, for a triangular prism::

            sage: p = Polyhedron(vertices = [[0,0],[0,1],[1,0]])
            sage: p2 = p.prism()
            sage: p2.gale_transform()
            ((-1, 0), (0, -1), (1, 1), (-1, -1), (1, 0), (0, 1))

        REFERENCES:

            Lectures in Geometric Combinatorics, R.R.Thomas, 2006, AMS Press.

        .. SEEALSO::

            :func`~sage.geometry.polyhedron.library.gale_transform_to_polyhedron`.

        TESTS::

            sage: P = Polyhedron(rays=[[1,0,0]])
            sage: P.gale_transform()
            Traceback (most recent call last):
            ...
            ValueError: not a polytope

        Check that :trac:`29073` is fixed::

            sage: P = polytopes.icosahedron(exact=False)
            sage: sum(P.gale_transform()).norm() < 1e-15
            True
        """
        if not self.is_compact():
            raise ValueError('not a polytope')

        A = matrix(self.n_vertices(),
                   [[1]+x for x in self.vertex_generator()])
        A = A.transpose()
        A_ker = A.right_kernel_matrix(basis='computed')
        return tuple(A_ker.columns())

    def _test_gale_transform(self, tester=None, **options):
        """
        Run tests on the method :meth:`.gale_transform` and its inverse
        :meth:`~sage.geometry.polyhedron.library.gale_transform_to_polytope`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_gale_transform()
        """
        if tester is None:
            tester = self._tester(**options)

        if not self.is_compact():
            with tester.assertRaises(ValueError):
                self.gale_transform()
            return

        # Check :trac:`29073`.
        if not self.base_ring().is_exact() and self.ambient_dim() > 0:
            g = self.gale_transform()
            tester.assertTrue(sum(g).norm() < 1e-10 or sum(g).norm()/matrix(g).norm() < 1e-13)
            return

        # Prevent very long doctests.
        if self.n_vertices() + self.n_rays() > 50 or self.n_facets() > 50:
            return

        if not self.is_empty():
            # ``gale_transform_to_polytope`` needs at least one vertex to work.
            from sage.geometry.polyhedron.library import gale_transform_to_polytope
            g = self.gale_transform()
            P = gale_transform_to_polytope(g, base_ring=self.base_ring(), backend=self.backend())

            tester.assertTrue(self.is_combinatorially_isomorphic(P))

    @cached_method
    def normal_fan(self, direction='inner'):
        r"""
        Return the normal fan of a compact full-dimensional rational polyhedron.

        This returns the inner normal fan of ``self``. For the outer normal fan,
        use ``direction='outer'``.

        INPUT:

        - ``direction`` -- either ``'inner'`` (default) or ``'outer'``; if
          set to ``'inner'``, use the inner normal vectors to span the cones of
          the fan, if set to ``'outer'``, use the outer normal vectors.

        OUTPUT:

        A complete fan of the ambient space as a
        :class:`~sage.geometry.fan.RationalPolyhedralFan`.

        .. SEEALSO::

            :meth:`face_fan`.

        EXAMPLES::

            sage: S = Polyhedron(vertices = [[0, 0], [1, 0], [0, 1]])
            sage: S.normal_fan()
            Rational polyhedral fan in 2-d lattice N

            sage: C = polytopes.hypercube(4)
            sage: NF = C.normal_fan(); NF
            Rational polyhedral fan in 4-d lattice N

        Currently, it is only possible to get the normal fan of a bounded rational polytope::

            sage: P = Polyhedron(rays = [[1, 0], [0, 1]])
            sage: P.normal_fan()
            Traceback (most recent call last):
            ...
            NotImplementedError: the normal fan is only supported for polytopes (compact polyhedra).

            sage: Q = Polyhedron(vertices = [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: Q.normal_fan()
            Traceback (most recent call last):
            ...
            ValueError: the normal fan is only defined for full-dimensional polytopes

            sage: R = Polyhedron(vertices = [[0, 0], [AA(sqrt(2)), 0], [0, AA(sqrt(2))]])
            sage: R.normal_fan()
            Traceback (most recent call last):
            ...
            NotImplementedError: normal fan handles only polytopes over the rationals

            sage: P = Polyhedron(vertices=[[0,0],[2,0],[0,2],[2,1],[1,2]])
            sage: P.normal_fan(direction=None)
            Traceback (most recent call last):
            ...
            TypeError: the direction should be 'inner' or 'outer'

            sage: inner_nf = P.normal_fan()
            sage: inner_nf.rays()
            N( 1,  0),
            N( 0, -1),
            N( 0,  1),
            N(-1,  0),
            N(-1, -1)
            in 2-d lattice N

            sage: outer_nf = P.normal_fan(direction='outer')
            sage: outer_nf.rays()
            N( 1,  0),
            N( 1,  1),
            N( 0,  1),
            N(-1,  0),
            N( 0, -1)
            in 2-d lattice N

        REFERENCES:

        For more information, see Chapter 7 of [Zie2007]_.
        """
        from sage.geometry.fan import NormalFan

        if not QQ.has_coerce_map_from(self.base_ring()):
            raise NotImplementedError('normal fan handles only polytopes over the rationals')
        if direction == 'inner':
            return NormalFan(self)
        elif direction == 'outer':
            return NormalFan(-self)
        else:
            raise TypeError("the direction should be 'inner' or 'outer'")

    @cached_method
    def face_fan(self):
        r"""
        Return the face fan of a compact rational polyhedron.

        OUTPUT:

        A fan of the ambient space as a
        :class:`~sage.geometry.fan.RationalPolyhedralFan`.

        .. SEEALSO::

            :meth:`normal_fan`.

        EXAMPLES::

            sage: T = polytopes.cuboctahedron()
            sage: T.face_fan()
            Rational polyhedral fan in 3-d lattice M

        The polytope should contain the origin in the interior::

            sage: P = Polyhedron(vertices = [[1/2, 1], [1, 1/2]])
            sage: P.face_fan()
            Traceback (most recent call last):
            ...
            ValueError: face fans are defined only for polytopes containing the origin as an interior point!

            sage: Q = Polyhedron(vertices = [[-1, 1/2], [1, -1/2]])
            sage: Q.contains([0,0])
            True
            sage: FF = Q.face_fan(); FF
            Rational polyhedral fan in 2-d lattice M

        The polytope has to have rational coordinates::

            sage: S = polytopes.dodecahedron()
            sage: S.face_fan()
            Traceback (most recent call last):
            ...
            NotImplementedError: face fan handles only polytopes over the rationals

        REFERENCES:

        For more information, see Chapter 7 of [Zie2007]_.
        """
        from sage.geometry.fan import FaceFan

        if not QQ.has_coerce_map_from(self.base_ring()):
            raise NotImplementedError('face fan handles only polytopes over the rationals')

        return FaceFan(self)

    def _triangulate_normaliz(self):
        r"""
        Gives a triangulation of the polyhedron using normaliz

        OUTPUT:

        A tuple of pairs ``(simplex,simplex_volume)`` used in the
        triangulation.

        .. NOTE::

            This function depends on Normaliz (i.e. the ``pynormaliz`` optional
            package). See the Normaliz documentation for further details.

        TESTS::

            sage: K = Polyhedron(vertices=[[1,1]], rays=[[1,0],[1,2]])
            sage: K._triangulate_normaliz()
            Traceback (most recent call last):
            ...
            TypeError: the polyhedron's backend should be 'normaliz'
        """
        raise TypeError("the polyhedron's backend should be 'normaliz'")

    def triangulate(self, engine='auto', connected=True, fine=False, regular=None, star=None):
        r"""
        Return a triangulation of the polytope.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal',
          'TOPCOM', or 'normaliz'.  The 'internal' and 'TOPCOM' instruct
          this package to always use its own triangulation algorithms
          or TOPCOM's algorithms, respectively. By default ('auto'),
          TOPCOM is used if it is available and internal routines otherwise.

        The remaining keyword parameters are passed through to the
        :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`
        constructor:

        - ``connected`` -- boolean (default: ``True``). Whether the
          triangulations should be connected to the regular
          triangulations via bistellar flips. These are much easier to
          compute than all triangulations.

        - ``fine`` -- boolean (default: ``False``). Whether the
          triangulations must be fine, that is, make use of all points
          of the configuration.

        - ``regular`` -- boolean or ``None`` (default:
          ``None``). Whether the triangulations must be regular. A
          regular triangulation is one that is induced by a
          piecewise-linear convex support function. In other words,
          the shadows of the faces of a polyhedron in one higher
          dimension.

          * ``True``: Only regular triangulations.

          * ``False``: Only non-regular triangulations.

          * ``None`` (default): Both kinds of triangulation.

        - ``star`` -- either ``None`` (default) or a point. Whether
          the triangulations must be star. A triangulation is star if
          all maximal simplices contain a common point. The central
          point can be specified by its index (an integer) in the
          given points or by its coordinates (anything iterable.)

        OUTPUT:

        A triangulation of the convex hull of the vertices as a
        :class:`~sage.geometry.triangulation.point_configuration.Triangulation`. The
        indices in the triangulation correspond to the
        :meth:`Vrepresentation` objects.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: triangulation = cube.triangulate(
            ....:    engine='internal') # to make doctest independent of TOPCOM
            sage: triangulation
            (<0,1,2,7>, <0,1,5,7>, <0,2,3,7>, <0,3,4,7>, <0,4,5,7>, <1,5,6,7>)
            sage: simplex_indices = triangulation[0]; simplex_indices
            (0, 1, 2, 7)
            sage: simplex_vertices = [ cube.Vrepresentation(i) for i in simplex_indices ]
            sage: simplex_vertices
            [A vertex at (1, -1, -1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             A vertex at (-1, 1, 1)]
            sage: Polyhedron(simplex_vertices)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices

        It is possible to use ``'normaliz'`` as an engine. For this, the
        polyhedron should have the backend set to normaliz::

            sage: P = Polyhedron(vertices=[[0,0,1],[1,0,1],[0,1,1],[1,1,1]],backend='normaliz')  # optional - pynormaliz
            sage: P.triangulate(engine='normaliz')  # optional - pynormaliz
            (<0,1,2>, <1,2,3>)

            sage: P = Polyhedron(vertices=[[0,0,1],[1,0,1],[0,1,1],[1,1,1]])
            sage: P.triangulate(engine='normaliz')
            Traceback (most recent call last):
            ...
            TypeError: the polyhedron's backend should be 'normaliz'

        The normaliz engine can triangulate pointed cones::

            sage: C1 = Polyhedron(rays=[[0,0,1],[1,0,1],[0,1,1],[1,1,1]],backend='normaliz')  # optional - pynormaliz
            sage: C1.triangulate(engine='normaliz')  # optional - pynormaliz
            (<0,1,2>, <1,2,3>)
            sage: C2 = Polyhedron(rays=[[1,0,1],[0,0,1],[0,1,1],[1,1,10/9]],backend='normaliz')  # optional - pynormaliz
            sage: C2.triangulate(engine='normaliz')  # optional - pynormaliz
            (<0,1,2>, <1,2,3>)

        They can also be affine cones::

            sage: K = Polyhedron(vertices=[[1,1,1]],rays=[[1,0,0],[0,1,0],[1,1,-1],[1,1,1]], backend='normaliz')  # optional - pynormaliz
            sage: K.triangulate(engine='normaliz')  # optional - pynormaliz
            (<0,1,2>, <0,1,3>)
        """
        if self.lines():
            raise NotImplementedError('triangulation of polyhedra with lines is not supported')
        if len(self.vertices_list()) >= 2 and self.rays_list():
            raise NotImplementedError('triangulation of non-compact polyhedra that are not cones is not supported')
        if not self.is_compact() and engine != 'normaliz':
            raise NotImplementedError("triangulation of pointed polyhedra requires 'normaliz'")
        from sage.geometry.triangulation.point_configuration import PointConfiguration
        if self.is_compact():
            pc = PointConfiguration((v.vector() for v in self.vertex_generator()),
                                    connected=connected, fine=fine, regular=regular, star=star)
            # If the engine is not normaliz, we pass directly to the
            # PointConfiguration module.
            if engine != 'normaliz':
                pc.set_engine(engine)
                return pc.triangulate()
            else:
                return pc(self._triangulate_normaliz())
        else:  # From above, we have a pointed cone and the engine is normaliz
            try:
                pc = PointConfiguration((v.vector() for v in self.ray_generator()),
                                        connected=connected, fine=fine, regular=regular, star=star)
                return pc(self._triangulate_normaliz())
            except AssertionError:
                # PointConfiguration is not adapted to inhomogeneous cones
                # This is a hack. TODO: Implement the necessary things in
                # PointConfiguration to accept such cases.
                c = self.representative_point()
                normed_v = ((1/(r.vector()*c))*r.vector() for r in self.ray_generator())
                pc = PointConfiguration(normed_v, connected=connected, fine=fine, regular=regular, star=star)
                return pc(self._triangulate_normaliz())

    @coerce_binop
    def minkowski_sum(self, other):
        r"""
        Return the Minkowski sum.

        Minkowski addition of two subsets of a vector space is defined
        as

        .. MATH::

            X \oplus Y =
            \cup_{y\in Y} (X+y) =
            \cup_{x\in X, y\in Y} (x+y)

        See :meth:`minkowski_difference` for a partial inverse operation.

        INPUT:

        - ``other`` -- a :class:`Polyhedron_base`

        OUTPUT:

        The Minkowski sum of ``self`` and ``other``

        EXAMPLES::

            sage: X = polytopes.hypercube(3)
            sage: Y = Polyhedron(vertices=[(0,0,0), (0,0,1/2), (0,1/2,0), (1/2,0,0)])
            sage: X+Y
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 13 vertices

            sage: four_cube = polytopes.hypercube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: four_cube + four_simplex
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 36 vertices
            sage: four_cube.minkowski_sum(four_simplex) == four_cube + four_simplex
            True

            sage: poly_spam = Polyhedron([[3,4,5,2],[1,0,0,1],[0,0,0,0],[0,4,3,2],[-3,-3,-3,-3]], base_ring=ZZ)
            sage: poly_eggs = Polyhedron([[5,4,5,4],[-4,5,-4,5],[4,-5,4,-5],[0,0,0,0]], base_ring=QQ)
            sage: poly_spam + poly_spam + poly_eggs
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 12 vertices
        """
        new_vertices = []
        for v1 in self.vertex_generator():
            for v2 in other.vertex_generator():
                new_vertices.append(list(v1() + v2()))
        if new_vertices != []:
            new_rays = self.rays() + other.rays()
            new_lines = self.lines() + other.lines()
            return self.parent().element_class(self.parent(), [new_vertices, new_rays, new_lines], None)
        else:
            return self.parent().element_class(self.parent(), None, None)

    _add_ = minkowski_sum

    @coerce_binop
    def minkowski_difference(self, other):
        r"""
        Return the Minkowski difference.

        Minkowski subtraction can equivalently be defined via
        Minkowski addition (see :meth:`minkowski_sum`) or as
        set-theoretic intersection via

        .. MATH::

            X \ominus Y =
            (X^c \oplus Y)^c =
            \cap_{y\in Y} (X-y)

        where superscript-"c" means the complement in the ambient
        vector space. The Minkowski difference of convex sets is
        convex, and the difference of polyhedra is again a
        polyhedron. We only consider the case of polyhedra in the
        following. Note that it is not quite the inverse of
        addition. In fact:

        * `(X+Y)-Y = X` for any polyhedra `X`, `Y`.

        * `(X-Y)+Y \subseteq X`

        * `(X-Y)+Y = X` if and only if Y is a Minkowski summand of X.

        INPUT:

        - ``other`` -- a :class:`Polyhedron_base`

        OUTPUT:

        The Minkowski difference of ``self`` and ``other``. Also known
        as Minkowski subtraction of ``other`` from ``self``.

        EXAMPLES::

            sage: X = polytopes.hypercube(3)
            sage: Y = Polyhedron(vertices=[(0,0,0), (0,0,1), (0,1,0), (1,0,0)]) / 2
            sage: (X+Y)-Y == X
            True
            sage: (X-Y)+Y < X
            True

        The polyhedra need not be full-dimensional::

            sage: X2 = Polyhedron(vertices=[(-1,-1,0),(1,-1,0),(-1,1,0),(1,1,0)])
            sage: Y2 = Polyhedron(vertices=[(0,0,0), (0,1,0), (1,0,0)]) / 2
            sage: (X2+Y2)-Y2 == X2
            True
            sage: (X2-Y2)+Y2 < X2
            True

        Minus sign is really an alias for :meth:`minkowski_difference`
        ::

            sage: four_cube = polytopes.hypercube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: four_cube - four_simplex
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 16 vertices
            sage: four_cube.minkowski_difference(four_simplex) == four_cube - four_simplex
            True

        Coercion of the base ring works::

            sage: poly_spam = Polyhedron([[3,4,5,2],[1,0,0,1],[0,0,0,0],[0,4,3,2],[-3,-3,-3,-3]], base_ring=ZZ)
            sage: poly_eggs = Polyhedron([[5,4,5,4],[-4,5,-4,5],[4,-5,4,-5],[0,0,0,0]], base_ring=QQ) / 100
            sage: poly_spam - poly_eggs
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices

        TESTS::

            sage: X = polytopes.hypercube(2)
            sage: Y = Polyhedron(vertices=[(1,1)])
            sage: (X-Y).Vrepresentation()
            (A vertex at (0, -2), A vertex at (0, 0), A vertex at (-2, 0), A vertex at (-2, -2))

            sage: Y = Polyhedron(vertices=[(1,1), (0,0)])
            sage: (X-Y).Vrepresentation()
            (A vertex at (0, -1), A vertex at (0, 0), A vertex at (-1, 0), A vertex at (-1, -1))

            sage: X = X + Y   # now Y is a Minkowski summand of X
            sage: (X+Y)-Y == X
            True
            sage: (X-Y)+Y == X
            True

        Testing that :trac:`28506` is fixed::

            sage: Q = Polyhedron([[1,0],[0,1]])
            sage: S = Polyhedron([[0,0],[1,2]])
            sage: S.minkowski_difference(Q)
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices
        """
        if other.is_empty():
            return self.parent().universe()   # empty intersection = everything
        if not other.is_compact():
            raise NotImplementedError('only subtracting compact polyhedra is implemented')
        new_eqns = []
        for eq in self.equations():
            values = [ eq.A() * v.vector() for v in other.vertices() ]
            eq = list(eq)
            eq[0] += min(values)   # shift constant term
            new_eqns.append(eq)
        P = self.parent()
        new_ieqs = []
        for ieq in self.inequalities():
            values = [ ieq.A() * v.vector() for v in other.vertices() ]
            ieq = list(ieq)
            ieq[0] += min(values)   # shift constant term
            new_ieqs.append(ieq)

        # Some vertices might need fractions.
        P = self.parent().change_ring(self.base_ring().fraction_field())
        return P.element_class(P, None, [new_ieqs, new_eqns])

    def __sub__(self, other):
        r"""
        Implement minus binary operation

        Polyhedra are not a ring with respect to dilatation and
        Minkowski sum, for example `X\oplus(-1)*Y \not= X\ominus Y`.

        INPUT:

        - ``other`` -- a translation vector or a polyhedron

        OUTPUT:

        Either translation by the negative of the given vector or
        Minkowski subtraction by the given polyhedron.

        EXAMPLES::

            sage: X = polytopes.hypercube(2)
            sage: v = vector([1,1])
            sage: (X - v/2).Vrepresentation()
            (A vertex at (-3/2, -3/2), A vertex at (-3/2, 1/2),
             A vertex at (1/2, -3/2), A vertex at (1/2, 1/2))
            sage: (X-v)+v == X
            True

            sage: Y = Polyhedron(vertices=[(1/2,0),(0,1/2)])
            sage: (X-Y).Vrepresentation()
            (A vertex at (1/2, -1), A vertex at (1/2, 1/2),
             A vertex at (-1, 1/2), A vertex at (-1, -1))
            sage: (X+Y)-Y == X
            True
        """
        if is_Polyhedron(other):
            return self.minkowski_difference(other)
        return self + (-other)

    def is_minkowski_summand(self, Y):
        r"""
        Test whether ``Y`` is a Minkowski summand.

        See :meth:`minkowski_sum`.

        OUTPUT:

        Boolean. Whether there exists another polyhedron `Z` such that
        ``self`` can be written as `Y\oplus Z`.

        EXAMPLES::

            sage: A = polytopes.hypercube(2)
            sage: B = Polyhedron(vertices=[(0,1), (1/2,1)])
            sage: C = Polyhedron(vertices=[(1,1)])
            sage: A.is_minkowski_summand(B)
            True
            sage: A.is_minkowski_summand(C)
            True
            sage: B.is_minkowski_summand(C)
            True
            sage: B.is_minkowski_summand(A)
            False
            sage: C.is_minkowski_summand(A)
            False
            sage: C.is_minkowski_summand(B)
            False
        """
        return self.minkowski_difference(Y).minkowski_sum(Y) == self

    def translation(self, displacement):
        """
        Return the translated polyhedron.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector

        OUTPUT:

        The translated polyhedron.

        EXAMPLES::

            sage: P = Polyhedron([[0,0],[1,0],[0,1]], base_ring=ZZ)
            sage: P.translation([2,1])
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: P.translation( vector(QQ,[2,1]) )
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices

        TESTS::

            sage: P = Polyhedron([[0,0],[1,0],[0,1]], base_ring=ZZ, backend='field')
            sage: P.translation([2,1]).backend()
            'field'

        Check that precomputed data is set up correctly::

            sage: P = polytopes.permutahedron(4)*Polyhedron(lines=[[1]])
            sage: Q = P.change_ring(P.base_ring(), backend='field')
            sage: P + vector([1,2,3,4,5]) == Q + vector([1,2,3,4,5])
            True
            sage: P + vector([1,2,3,4,5/2]) == Q + vector([1,2,3,4,5/2])
            True
        """
        Vrep, Hrep, parent = self._translation_double_description(displacement)

        pref_rep = 'Vrep' if self.n_vertices() + self.n_rays() <= self.n_inequalities() else 'Hrep'

        return parent.element_class(parent, Vrep, Hrep,
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _translation_double_description(self, displacement):
        r"""
        Return the input parameters for the translation.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector

        OUTPUT: Tuple of consisting of new Vrepresentation, Hrepresentation and parent.

        .. SEEALSO::

            :meth:`translation`

        EXAMPLES::

            sage: P = Polyhedron([[0,0],[1,0],[0,1]], base_ring=ZZ)
            sage: Vrep, Hrep, parent = P._translation_double_description([2,1])
            sage: [tuple(x) for x in Vrep], [tuple(x) for x in Hrep], parent
            ([((2, 1), (2, 2), (3, 1)), (), ()],
             [((-2, 1, 0), (-1, 0, 1), (4, -1, -1)), ()],
             Polyhedra in ZZ^2)
        """
        displacement = vector(displacement)
        new_vertices = (x.vector()+displacement for x in self.vertex_generator())
        new_rays = self.rays()
        new_lines = self.lines()
        parent = self.parent().base_extend(displacement)

        # Replace a hyperplane of the form A*x + b >= 0 by
        # A(x-displacement) + b >= 0 <=> Ax + b - A*displacement >= 0.
        # Likewise for equations.
        def get_new(x):
            y = x.vector().change_ring(parent.base_ring())
            y[0] -= x.A()*displacement
            return y

        new_ieqs = (get_new(x) for x in self.inequality_generator())
        new_eqns = (get_new(x) for x in self.equation_generator())
        return [new_vertices, new_rays, new_lines], [new_ieqs, new_eqns], parent

    def product(self, other):
        """
        Return the Cartesian product.

        INPUT:

        - ``other`` -- a :class:`Polyhedron_base`

        OUTPUT:

        The Cartesian product of ``self`` and ``other`` with a
        suitable base ring to encompass the two.

        EXAMPLES::

            sage: P1 = Polyhedron([[0],[1]], base_ring=ZZ)
            sage: P2 = Polyhedron([[0],[1]], base_ring=QQ)
            sage: P1.product(P2)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        The Cartesian product is the product in the semiring of polyhedra::

            sage: P1 * P1
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: P1 * P2
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: P2 * P2
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: 2 * P1
            A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 2 vertices
            sage: P1 * 2.0
            A 1-dimensional polyhedron in RDF^1 defined as the convex hull of 2 vertices

        An alias is :meth:`cartesian_product`::

            sage: P1.cartesian_product(P2) == P1.product(P2)
            True

        TESTS:

        Check that :trac:`15253` is fixed::

            sage: polytopes.hypercube(1) * polytopes.hypercube(2)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        """
        try:
            new_ring = self.parent()._coerce_base_ring(other)
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                             + " and " + str(other.parent()))

        from itertools import chain

        new_vertices = (tuple(x) + tuple(y)
                        for x in self.vertex_generator() for y in other.vertex_generator())

        self_zero  = tuple(0 for _ in range( self.ambient_dim()))
        other_zero = tuple(0 for _ in range(other.ambient_dim()))

        rays = chain((tuple(r) + other_zero for r in  self.ray_generator()),
                     (self_zero + tuple(r)  for r in other.ray_generator()))

        lines = chain((tuple(l) + other_zero for l in  self.line_generator()),
                      (self_zero + tuple(l)  for l in other.line_generator()))

        if self.n_vertices() == 0 or other.n_vertices() == 0:
            # In this case we obtain the empty polyhedron.
            # There is not vertex to attach the rays or lines to.
            # By our convention, in this case the polyhedron shall also not have rays or lines.
            rays = ()
            lines = ()

        ieqs = chain((tuple(i) + other_zero               for i in  self.inequality_generator()),
                     ((i.b(),) + self_zero + tuple(i.A()) for i in other.inequality_generator()))

        eqns = chain((tuple(e) + other_zero               for e in  self.equation_generator()),
                     ((e.b(),) + self_zero + tuple(e.A()) for e in other.equation_generator()))

        pref_rep = 'Vrep' if self.n_vertices() + self.n_rays() + other.n_vertices() + other.n_rays() \
                             <= self.n_inequalities() + other.n_inequalities() else 'Hrep'

        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim())
        return parent.element_class(parent, [new_vertices, rays, lines],
                                    [ieqs, eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    _mul_ = product

    cartesian_product = product

    def _test_product(self, tester=None, **options):
        """
        Run tests on the method :meth:`.product`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_product()
        """
        from sage.geometry.polyhedron.library import polytopes
        if tester is None:
            tester = self._tester(**options)

        if self.n_vertices() + self.n_rays() < 40 and self.n_facets() < 40:
            # Check that the product preserves the backend, where possible.
            P = polytopes.simplex(backend="cdd")
            tester.assertEqual((self*P).backend(), self.backend())
            Q = polytopes.simplex(backend="ppl")
            tester.assertEqual((self*Q).backend(), self.backend())

            # And that it changes the backend correctly where necessary.
            if self.base_ring() is not AA and AA.has_coerce_map_from(self.base_ring()):
                R = self*polytopes.regular_polygon(5, exact=True)
                assert R
            if RDF.has_coerce_map_from(self.base_ring()):
                R = self*polytopes.regular_polygon(5, exact=False)
                assert R

        if self.base_ring() in (ZZ, QQ):
            # Check that the double description is set up correctly.
            self_field = self.base_extend(self.base_ring(), backend='field')
            P = polytopes.permutahedron(4, backend='field').base_extend(QQ)
            Q = Polyhedron(rays=[[1,0,0,0],[0,1,1,0]], lines=[[0,1,0,1]], backend='field')
            (self_field * P)._test_basic_properties(tester)
            (self_field * Q)._test_basic_properties(tester)

    def join(self, other):
        """
        Return the join of ``self`` and ``other``.

        The join of two polyhedra is obtained by first placing the two objects in
        two non-intersecting affine subspaces `V`, and `W` whose affine hull is
        the whole ambient space, and finally by taking the convex hull of their
        union. The dimension of the join is the sum of the dimensions of the
        two polyhedron plus 1.

        INPUT:

        - ``other`` -- a polyhedron

        EXAMPLES::

            sage: P1 = Polyhedron([[0],[1]], base_ring=ZZ)
            sage: P2 = Polyhedron([[0],[1]], base_ring=QQ)
            sage: P1.join(P2)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: P1.join(P1)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: P2.join(P2)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices

        An unbounded example::

            sage: R1 = Polyhedron(rays=[[1]])
            sage: R1.join(R1)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices and 2 rays

        TESTS::

            sage: C = polytopes.hypercube(5)
            sage: S = Polyhedron([[1]])
            sage: C.join(S).is_combinatorially_isomorphic(C.pyramid())
            True

            sage: P = polytopes.simplex(backend='cdd')
            sage: Q = polytopes.simplex(backend='ppl')
            sage: P.join(Q).backend()
            'cdd'
            sage: Q.join(P).backend()
            'ppl'
        """
        try:
            new_ring = self.parent()._coerce_base_ring(other)
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                     + " and " + str(other.parent()))

        dim_self = self.ambient_dim()
        dim_other = other.ambient_dim()

        new_vertices = [list(x)+[0]*dim_other+[0] for x in self.vertex_generator()] + \
                       [[0]*dim_self+list(x)+[1] for x in other.vertex_generator()]
        new_rays = []
        new_rays.extend( [ r+[0]*dim_other+[0]
                           for r in self.ray_generator() ] )
        new_rays.extend( [ [0]*dim_self+r+[1]
                           for r in other.ray_generator() ] )
        new_lines = []
        new_lines.extend( [ l+[0]*dim_other+[0]
                            for l in self.line_generator() ] )
        new_lines.extend( [ [0]*dim_self+l+[1]
                            for l in other.line_generator() ] )

        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim() + 1)
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    def subdirect_sum(self, other):
        """
        Return the subdirect sum of ``self`` and ``other``.

        The subdirect sum of two polyhedron is a projection of the join of the
        two polytopes. It is obtained by placing the two objects in orthogonal subspaces
        intersecting at the origin.

        INPUT:

        - ``other`` -- a :class:`Polyhedron_base`

        EXAMPLES::

            sage: P1 = Polyhedron([[1],[2]], base_ring=ZZ)
            sage: P2 = Polyhedron([[3],[4]], base_ring=QQ)
            sage: sds = P1.subdirect_sum(P2);sds
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4
            vertices
            sage: sds.vertices()
            (A vertex at (0, 3),
             A vertex at (0, 4),
             A vertex at (1, 0),
             A vertex at (2, 0))

        .. SEEALSO::

            :meth:`join`
            :meth:`direct_sum`

        TESTS::

            sage: P = polytopes.simplex(backend='cdd')
            sage: Q = polytopes.simplex(backend='ppl')
            sage: P.subdirect_sum(Q).backend()
            'cdd'
            sage: Q.subdirect_sum(P).backend()
            'ppl'
        """
        try:
            new_ring = self.parent()._coerce_base_ring(other)
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                     + " and " + str(other.parent()))

        dim_self = self.ambient_dim()
        dim_other = other.ambient_dim()

        new_vertices = [list(x)+[0]*dim_other for x in self.vertex_generator()] + \
                       [[0]*dim_self+list(x) for x in other.vertex_generator()]
        new_rays = []
        new_rays.extend( [ r+[0]*dim_other
                           for r in self.ray_generator() ] )
        new_rays.extend( [ [0]*dim_self+r
                           for r in other.ray_generator() ] )
        new_lines = []
        new_lines.extend( [ l+[0]*dim_other
                            for l in self.line_generator() ] )
        new_lines.extend( [ [0]*dim_self+l
                            for l in other.line_generator() ] )

        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim())
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    def direct_sum(self, other):
        """
        Return the direct sum of ``self`` and ``other``.

        The direct sum of two polyhedron is the subdirect sum of the two, when
        they have the origin in their interior. To avoid checking if the origin
        is contained in both, we place the affine subspace containing ``other``
        at the center of ``self``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron_base`

        EXAMPLES::

            sage: P1 = Polyhedron([[1],[2]], base_ring=ZZ)
            sage: P2 = Polyhedron([[3],[4]], base_ring=QQ)
            sage: ds = P1.direct_sum(P2);ds
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: ds.vertices()
            (A vertex at (1, 0),
             A vertex at (2, 0),
             A vertex at (3/2, -1/2),
             A vertex at (3/2, 1/2))

        .. SEEALSO::

            :meth:`join`
            :meth:`subdirect_sum`

        TESTS:

        Check that the backend is preserved::

            sage: P = polytopes.simplex(backend='cdd')
            sage: Q = polytopes.simplex(backend='ppl')
            sage: P.direct_sum(Q).backend()
            'cdd'
            sage: Q.direct_sum(P).backend()
            'ppl'

        Check that :trac:`28506` is fixed::

            sage: s2 = polytopes.simplex(2)
            sage: s3 = polytopes.simplex(3)
            sage: s2.direct_sum(s3)
            A 5-dimensional polyhedron in QQ^7 defined as the convex hull of 7 vertices
        """
        try:
            # Some vertices might need fractions.
            new_ring = self.parent()._coerce_base_ring(other).fraction_field()
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                     + " and " + str(other.parent()))

        dim_self = self.ambient_dim()
        dim_other = other.ambient_dim()

        new_vertices = [list(x) + [0]*dim_other for x in self.vertex_generator()] + \
                       [list(self.center()) + list(x.vector() - other.center()) for x in other.vertex_generator()]
        new_rays = []
        new_rays.extend( [ r + [0]*dim_other
                           for r in self.ray_generator() ] )
        new_rays.extend( [ [0]*dim_self + r
                           for r in other.ray_generator() ] )
        new_lines = []
        new_lines.extend( [ l + [0]*dim_other
                            for l in self.line_generator() ] )
        new_lines.extend( [ [0]*dim_self + l
                            for l in other.line_generator() ] )

        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim())
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    def dilation(self, scalar):
        """
        Return the dilated (uniformly stretched) polyhedron.

        INPUT:

        - ``scalar`` -- A scalar, not necessarily in :meth:`base_ring`

        OUTPUT:

        The polyhedron dilated by that scalar, possibly coerced to a
        bigger base ring.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,6)])
            sage: next(p.vertex_generator())
            A vertex at (2, 4, 8)
            sage: p2 = p.dilation(2)
            sage: next(p2.vertex_generator())
            A vertex at (4, 8, 16)
            sage: p.dilation(2) == p * 2
            True

        TESTS:

        Dilation of empty polyhedra works, see :trac:`14987`::

            sage: p = Polyhedron(ambient_dim=2); p
            The empty polyhedron in ZZ^2
            sage: p.dilation(3)
            The empty polyhedron in ZZ^2

            sage: p = Polyhedron(vertices=[(1,1)], rays=[(1,0)], lines=[(0,1)])
            sage: (-p).rays()
            (A ray in the direction (-1, 0),)
            sage: (-p).lines()
            (A line in the direction (0, 1),)

            sage: (0*p).rays()
            ()
            sage: (0*p).lines()
            ()
        """
        parent = self.parent().base_extend(scalar)

        if scalar == 0:
            new_vertices = tuple(self.ambient_space().zero() for v in self.vertex_generator())
            new_rays = []
            new_lines = []
            return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

        one = parent.base_ring().one()
        sign = one if scalar > 0 else -one

        make_new_Hrep = lambda h: tuple(scalar*sign*x if i == 0 else sign*x
                                        for i, x in enumerate(h._vector))

        new_vertices = (tuple(scalar*x for x in v._vector) for v in self.vertex_generator())
        new_rays = (tuple(sign*x for x in r._vector) for r in self.ray_generator())
        new_lines = self.line_generator()
        new_inequalities = map(make_new_Hrep, self.inequality_generator())
        new_equations = map(make_new_Hrep, self.equation_generator())

        pref_rep = 'Vrep' if self.n_vertices() + self.n_rays() <= self.n_inequalities() else 'Hrep'

        return parent.element_class(parent, [new_vertices, new_rays, new_lines],
                                    [new_inequalities, new_equations],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _test_dilation(self, tester=None, **options):
        """
        Run tests on the method :meth:`.dilation`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_dilation()
        """
        if tester is None:
            tester = self._tester(**options)

        # Testing that the backend is preserved.
        tester.assertEqual(self.dilation(2*self.base_ring().gen()).backend(), self.backend())
        tester.assertEqual(self.dilation(ZZ(3)).backend(), self.backend())

        if self.n_vertices() + self.n_rays() > 40:
            # Avoid long time computations.
            return

        # Testing that the double description is set up correctly.
        if self.base_ring().is_exact():
            if self.base_ring() in (QQ, ZZ):
                p = self.base_extend(self.base_ring(), backend='field')
                (ZZ(2)*p)._test_basic_properties(tester)
                (ZZ(2)/2*p)._test_basic_properties(tester)
                (ZZ(-3)*p)._test_basic_properties(tester)
                (ZZ(-1)/2*p)._test_basic_properties(tester)
        else:
            tester.assertIsInstance(ZZ(1)/3*self, Polyhedron_base)

        if self.n_vertices() > 20 or self.base_ring() is AA:
            # Avoid long time computations.
            return

        # Some sanity check on the volume (only run for relatively small instances).
        if self.dim() > -1 and self.is_compact() and self.base_ring().is_exact():
            tester.assertEqual(self.dilation(3).volume(measure='induced'), self.volume(measure='induced')*3**self.dim())

        # Testing coercion with algebraic numbers.
        from sage.rings.number_field.number_field import QuadraticField
        K1 = QuadraticField(2, embedding=AA(2).sqrt())
        sqrt2 = K1.gen()
        K2 = QuadraticField(3, embedding=AA(3).sqrt())
        sqrt3 = K2.gen()

        if self.base_ring() in (QQ,ZZ,AA,RDF):
            tester.assertIsInstance(sqrt2*self, Polyhedron_base)
            tester.assertIsInstance(sqrt3*self, Polyhedron_base)
        elif hasattr(self.base_ring(), "composite_fields"):
            for scalar, K in ((sqrt2, K1), (sqrt3, K2)):
                new_ring = None
                try:
                    new_ring = self.base_ring().composite_fields()[0]
                except:
                    # This isn't about testing composite fields.
                    pass
                if new_ring:
                    p = self.change_ring(new_ring)
                    tester.assertIsInstance(scalar*p, Polyhedron_base)

    def linear_transformation(self, linear_transf, new_base_ring=None):
        """
        Return the linear transformation of ``self``.

        INPUT:

        - ``linear_transf`` -- a matrix, not necessarily in :meth:`base_ring`
        - ``new_base_ring`` -- ring (optional); specify the new base ring;
          may avoid coercion failure

        OUTPUT:

        The polyhedron transformed by that matrix, possibly coerced to a
        bigger base ring.

        EXAMPLES::

            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: proj_mat=matrix([[0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0]])
            sage: b3_proj = proj_mat * b3; b3_proj
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices

            sage: square = polytopes.regular_polygon(4)
            sage: square.vertices_list()
            [[0, -1], [1, 0], [-1, 0], [0, 1]]
            sage: transf = matrix([[1,1],[0,1]])
            sage: sheared = transf * square
            sage: sheared.vertices_list()
            [[-1, -1], [1, 0], [-1, 0], [1, 1]]
            sage: sheared == square.linear_transformation(transf)
            True

        Specifying the new base ring may avoid coercion failure::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: L.<sqrt3> = QuadraticField(3)
            sage: P = polytopes.cube()*sqrt2
            sage: M = matrix([[sqrt3, 0, 0], [0, sqrt3, 0], [0, 0, 1]])
            sage: P.linear_transformation(M, new_base_ring=K.composite_fields(L)[0])
            A 3-dimensional polyhedron in (Number Field in sqrt2sqrt3 with defining polynomial x^4 - 10*x^2 + 1 with sqrt2sqrt3 = 0.3178372451957823?)^3 defined as the convex hull of 8 vertices

        Linear transformation without specified new base ring fails in this case::

            sage: M*P
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Full MatrixSpace of 3 by 3 dense matrices over Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?' and 'Full MatrixSpace of 3 by 8 dense matrices over Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?'

        TESTS:

        Linear transformation respects backend::

            sage: P = polytopes.simplex(backend='field')
            sage: t = matrix([[1,1,1,1],[0,1,1,1],[0,0,1,1],[0,0,0,1]])
            sage: P.linear_transformation(t).backend()
            'field'

        Check that coercion works::

            sage: (1.0 * proj_mat) * b3
            A 3-dimensional polyhedron in RDF^4 defined as the convex hull of 5 vertices
            sage: (1/1 * proj_mat) * b3
            A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices
            sage: (AA(2).sqrt() * proj_mat) * b3
            A 3-dimensional polyhedron in AA^4 defined as the convex hull of 5 vertices

        Check that zero-matrices act correctly::

            sage: Matrix([]) * b3
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: Matrix([[0 for _ in range(9)]]) * b3
            A 0-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex
            sage: Matrix([[0 for _ in range(9)] for _ in range(4)]) * b3
            A 0-dimensional polyhedron in ZZ^4 defined as the convex hull of 1 vertex
            sage: Matrix([[0 for _ in range(8)]]) * b3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Full MatrixSpace of 1 by 8 dense matrices over Integer Ring' and 'Full MatrixSpace of 9 by 6 dense matrices over Integer Ring'
            sage: Matrix(ZZ, []) * b3
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: Matrix(ZZ, [[],[]]) * b3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Full MatrixSpace of 2 by 0 dense matrices over Integer Ring' and 'Full MatrixSpace of 9 by 6 dense matrices over Integer Ring'

        Check that the precomputed double description is correct::

            sage: P = polytopes.permutahedron(4)
            sage: Q = P.change_ring(QQ, backend='field')
            sage: P.affine_hull_projection() == Q.affine_hull_projection()
            True

            sage: M = matrix([[1, 2, 3, 4], [2, 3, 4, 5], [0, 0, 5, 1], [0, 2, 0, 3]])
            sage: M*P == M*Q
            True

            sage: M = matrix([[1, 2, 3, 4], [2, 3, 4, 5], [0, 0, 5, 1], [0, 2, 0, 3], [0, 1, 0, -3]])
            sage: M*P == M*Q
            True
        """
        is_injective = False
        if linear_transf.nrows() != 0:
            if new_base_ring:
                R = new_base_ring
            else:
                R = self.base_ring()

            # Multiplying a matrix with a vector is slow.
            # So we multiply the entire vertex matrix etc.
            # Still we create generators, as possibly the Vrepresentation will be discarded later on.
            if self.n_vertices():
                new_vertices = ( v for v in ((linear_transf*self.vertices_matrix(R)).transpose()) )
            else:
                new_vertices = ()
            if self.n_rays():
                new_rays = ( r for r in matrix(R, self.rays())*linear_transf.transpose() )
            else:
                new_rays = ()
            if self.n_lines():
                new_lines = ( l for l in matrix(R, self.lines())*linear_transf.transpose() )
            else:
                new_lines = ()

            if self.is_compact() and self.n_vertices() and self.n_inequalities():
                homogeneous_basis = matrix(R, ( [1] + list(v) for v in self.an_affine_basis() )).transpose()

                # To convert first to a list and then to a matrix seems to be necessary to obtain a meaningful error,
                # in case the number of columns doesn't match the dimension.
                new_homogeneous_basis = matrix(list( [1] + list(linear_transf*vector(R, v)) for v in self.an_affine_basis()) ).transpose()

                if self.dim() + 1 == new_homogeneous_basis.rank():
                    # The transformation is injective on the polytope.
                    is_injective = True

                    # Let V be the homogeneous vertex matrix (each vertex a column)
                    # and M the linear transformation.
                    # Then M*V is the new homogeneous vertex matrix.

                    # Let H be the inequalities matrix (each inequality a row).
                    # If we find N such that N*M*V = V than the new inequalities are
                    # given by H*N.

                    # Note that such N must exist, as our map is injective on the polytope.
                    # It is uniquely defined by considering a basis of the homogeneous vertices.
                    N = new_homogeneous_basis.solve_left(homogeneous_basis)
                    new_inequalities = ( h for h in matrix(R, self.inequalities())*N )

                    # The equations are the left kernel matrix of the homogeneous vertices
                    # or equivalently a basis thereof.
                    new_equations = (new_homogeneous_basis.transpose()).right_kernel_matrix()

        else:
            new_vertices = [[] for v in self.vertex_generator() ]
            new_rays = []
            new_lines = []

        new_dim = linear_transf.nrows()
        par = self.parent()

        if new_base_ring:
            new_parent = par.change_ring(new_base_ring, ambient_dim=new_dim)
        else:
            new_parent = par.base_extend(linear_transf.base_ring(), ambient_dim=new_dim)

        if is_injective:
            # Set up with both Vrepresentation and Hrepresentation.
            pref_rep = 'Vrep' if self.n_vertices() <= self.n_inequalities() else 'Hrep'

            return new_parent.element_class(new_parent, [new_vertices, new_rays, new_lines],
                                            [new_inequalities, new_equations],
                                            Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

        return new_parent.element_class(new_parent, [tuple(new_vertices), tuple(new_rays), tuple(new_lines)], None)

    def _test_linear_transformation(self, tester=None, **options):
        """
        Run some tests on linear transformation.

        TESTS::

            sage: Polyhedron(rays=[(0,1)])._test_linear_transformation()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.n_vertices() > 200 or self.n_facets() > 200:
            # Avoid very long doctests.
            return

        # Check that :trac:`30146` is fixed.
        from sage.matrix.special import identity_matrix
        tester.assertEqual(self, self.linear_transformation(identity_matrix(self.ambient_dim())))

    def _acted_upon_(self, actor, self_on_left):
        """
        Implement the action by scalars, vectors, matrices or other polyhedra.

        INPUT:

        - ``actor`` -- one of the following:
          - a scalar, not necessarily in :meth:`base_ring`,
          - a :class:`Polyhedron`,
          - a :class:`sage.modules.free_module_element.vector`,
          - a :class:`sage.matrix.constructor.matrix`,
        - ``self_on_right`` -- must be ``False`` for actor a matrix;
          ignored otherwise

        OUTPUT:

        - Dilation for a scalar
        - Product for a polyhedron
        - Translation for a vector
        - Linear transformation for a matrix

        EXAMPLES:

        ``actor`` is a scalar::

             sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,6)])
             sage: p._acted_upon_(2, True) == p.dilation(2)
             True
             sage: p*2 == p.dilation(2)
             True

        ``actor`` is a polyhedron::

             sage: p*p == p.product(p)
             True

        ``actor`` is a vector::

             sage: p + vector(ZZ,[1,2,3]) == p.translation([1,2,3])
             True

        ``actor`` is a matrix::

             sage: matrix(ZZ,[[1,2,3]]) * p
             A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 2 vertices

        A matrix must act from the left::

             sage: p * matrix(ZZ, [[1,2,3]]*3)
             Traceback (most recent call last):
             ...
             ValueError: matrices should act on the left
        """
        if is_Polyhedron(actor):
            return self.product(actor)
        elif is_Vector(actor):
            return self.translation(actor)
        elif is_Matrix(actor):
            if self_on_left:
                raise ValueError("matrices should act on the left")
            else:
                return self.linear_transformation(actor)
        else:
            return self.dilation(actor)

    def __neg__(self):
        """
        Negation of a polytope is defined as inverting the coordinates.

        EXAMPLES::

            sage: t = polytopes.simplex(3,project=False);  t.vertices()
            (A vertex at (0, 0, 0, 1), A vertex at (0, 0, 1, 0),
             A vertex at (0, 1, 0, 0), A vertex at (1, 0, 0, 0))
            sage: neg_ = -t
            sage: neg_.vertices()
            (A vertex at (-1, 0, 0, 0), A vertex at (0, -1, 0, 0),
             A vertex at (0, 0, -1, 0), A vertex at (0, 0, 0, -1))

        TESTS::

            sage: p = Polyhedron(ieqs=[[1,1,0]])
            sage: p.rays()
            (A ray in the direction (1, 0),)
            sage: pneg = p.__neg__()
            sage: pneg.rays()
            (A ray in the direction (-1, 0),)
        """
        return self.dilation(-1)

    def __truediv__(self, scalar):
        """
        Divide by a scalar factor.

        See :meth:`dilation` for details.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,4)])
            sage: (p/5).Vrepresentation()
            (A vertex at (2/5, 4/5, 8/5), A vertex at (3/5, 9/5, 27/5))
            sage: (p/int(5)).Vrepresentation()
            (A vertex at (0.4, 0.8, 1.6), A vertex at (0.6, 1.8, 5.4))
        """
        return self.dilation(1/scalar)

    @coerce_binop
    def convex_hull(self, other):
        """
        Return the convex hull of the set-theoretic union of the two
        polyhedra.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`

        OUTPUT:

        The convex hull.

        EXAMPLES::

            sage: a_simplex = polytopes.simplex(3, project=True)
            sage: verts = a_simplex.vertices()
            sage: verts = [[x[0]*3/5+x[1]*4/5, -x[0]*4/5+x[1]*3/5, x[2]] for x in verts]
            sage: another_simplex = Polyhedron(vertices = verts)
            sage: simplex_union = a_simplex.convex_hull(another_simplex)
            sage: simplex_union.n_vertices()
            7
        """
        hull_vertices = self.vertices() + other.vertices()
        hull_rays = self.rays() + other.rays()
        hull_lines = self.lines() + other.lines()
        return self.parent().element_class(self.parent(), [hull_vertices, hull_rays, hull_lines], None)

    @coerce_binop
    def intersection(self, other):
        r"""
        Return the intersection of one polyhedron with another.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`

        OUTPUT:

        The intersection.

        Note that the intersection of two `\ZZ`-polyhedra might not be
        a `\ZZ`-polyhedron. In this case, a `\QQ`-polyhedron is
        returned.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: oct = polytopes.cross_polytope(3)
            sage: cube.intersection(oct*2)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 12 vertices

        As a shorthand, one may use::

            sage: cube & oct*2
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 12 vertices

        The intersection of two `\ZZ`-polyhedra is not necessarily a `\ZZ`-polyhedron::

            sage: P = Polyhedron([(0,0),(1,1)], base_ring=ZZ)
            sage: P.intersection(P)
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: Q = Polyhedron([(0,1),(1,0)], base_ring=ZZ)
            sage: P.intersection(Q)
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex
            sage: _.Vrepresentation()
            (A vertex at (1/2, 1/2),)

        TESTS:

        Check that :trac:`19012` is fixed::

            sage: K.<a> = QuadraticField(5)
            sage: P = Polyhedron([[0,0],[0,a],[1,1]])
            sage: Q = Polyhedron(ieqs=[[-1,a,1]])
            sage: P.intersection(Q)
            A 2-dimensional polyhedron in (Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?)^2 defined as the convex hull of 4 vertices
        """
        new_ieqs = self.inequalities() + other.inequalities()
        new_eqns = self.equations() + other.equations()
        parent = self.parent()
        try:
            intersection = parent.element_class(parent, None, [new_ieqs, new_eqns])

            # Force calculation of the vertices.
            _ = intersection.n_vertices()
            return intersection
        except TypeError as msg:
            if self.base_ring() is ZZ:
                parent = parent.base_extend(QQ)
                return parent.element_class(parent, None, [new_ieqs, new_eqns])
            else:
                raise TypeError(msg)

    __and__ = intersection

    def truncation(self, cut_frac=None):
        r"""
        Return a new polyhedron formed from two points on each edge
        between two vertices.

        INPUT:

        - ``cut_frac`` -- integer, how deeply to cut into the edge.
          Default is `\frac{1}{3}`.

        OUTPUT:

        A Polyhedron object, truncated as described above.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: trunc_cube = cube.truncation()
            sage: trunc_cube.n_vertices()
            24
            sage: trunc_cube.n_inequalities()
            14

        TESTS::

            sage: polytopes.simplex(backend='field').truncation().backend()
            'field'
        """
        if cut_frac is None:
            cut_frac = ZZ.one() / 3

        new_vertices = []
        for e in self.bounded_edges():
            new_vertices.append((1 - cut_frac) * e[0]() + cut_frac * e[1]())
            new_vertices.append(cut_frac * e[0]() + (1 - cut_frac) * e[1]())

        new_vertices = [list(v) for v in new_vertices]
        new_rays = self.rays()
        new_lines = self.lines()

        parent = self.parent().base_extend(cut_frac)
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    def face_truncation(self, face, linear_coefficients=None, cut_frac=None):
        r"""
        Return a new polyhedron formed by truncating a face by an hyperplane.

        By default, the normal vector of the hyperplane used to truncate the
        polyhedron is obtained by taking the barycenter vector of the cone
        corresponding to the truncated face in the normal fan of the
        polyhedron. It is possible to change the direction using the option
        ``linear_coefficients``.

        To determine how deep the truncation is done, the method uses the
        parameter ``cut_frac``. By default it is equal to `\frac{1}{3}`. Once
        the normal vector of the cutting hyperplane is chosen, the vertices of
        polyhedron are evaluated according to the corresponding linear
        function. The parameter `\frac{1}{3}` means that the cutting
        hyperplane is placed `\frac{1}{3}` of the way from the vertices of the
        truncated face to the next evaluated vertex.

        INPUT:

        - ``face`` -- a PolyhedronFace
        - ``linear_coefficients`` -- tuple of integer. Specifies the coefficient
          of the normal vector of the cutting hyperplane used to truncate the
          face.
          The default direction is determined using the normal fan of the
          polyhedron.
        - ``cut_frac`` -- number between 0 and 1. Determines where the
           hyperplane cuts the polyhedron. A value close to 0 cuts very close
           to the face, whereas a value close to 1 cuts very close to the next
           vertex (according to the normal vector of the cutting hyperplane).
           Default is `\frac{1}{3}`.

        OUTPUT:

        A Polyhedron object, truncated as described above.

        EXAMPLES::

            sage: Cube = polytopes.hypercube(3)
            sage: vertex_trunc1 = Cube.face_truncation(Cube.faces(0)[0])
            sage: vertex_trunc1.f_vector()
            (1, 10, 15, 7, 1)
            sage: tuple(f.ambient_V_indices() for f in vertex_trunc1.faces(2))
            ((4, 5, 6, 7, 9),
             (0, 3, 4, 8, 9),
             (0, 1, 6, 7, 8),
             (7, 8, 9),
             (2, 3, 4, 5),
             (1, 2, 5, 6),
             (0, 1, 2, 3))
            sage: vertex_trunc1.vertices()
            (A vertex at (1, -1, -1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             A vertex at (1, -1, 1),
             A vertex at (-1, -1, 1),
             A vertex at (-1, 1, 1),
             A vertex at (-1, 1, -1),
             A vertex at (-1, -1/3, -1),
             A vertex at (-1/3, -1, -1),
             A vertex at (-1, -1, -1/3))
            sage: vertex_trunc2 = Cube.face_truncation(Cube.faces(0)[0],cut_frac=1/2)
            sage: vertex_trunc2.f_vector()
            (1, 10, 15, 7, 1)
            sage: tuple(f.ambient_V_indices() for f in vertex_trunc2.faces(2))
            ((4, 5, 6, 7, 9),
             (0, 3, 4, 8, 9),
             (0, 1, 6, 7, 8),
             (7, 8, 9),
             (2, 3, 4, 5),
             (1, 2, 5, 6),
             (0, 1, 2, 3))
            sage: vertex_trunc2.vertices()
            (A vertex at (1, -1, -1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             A vertex at (1, -1, 1),
             A vertex at (-1, -1, 1),
             A vertex at (-1, 1, 1),
             A vertex at (-1, 1, -1),
             A vertex at (-1, 0, -1),
             A vertex at (0, -1, -1),
             A vertex at (-1, -1, 0))
            sage: vertex_trunc3 = Cube.face_truncation(Cube.faces(0)[0],cut_frac=0.3)
            sage: vertex_trunc3.vertices()
            (A vertex at (-1.0, -1.0, 1.0),
             A vertex at (-1.0, 1.0, -1.0),
             A vertex at (-1.0, 1.0, 1.0),
             A vertex at (1.0, 1.0, -1.0),
             A vertex at (1.0, 1.0, 1.0),
             A vertex at (1.0, -1.0, 1.0),
             A vertex at (1.0, -1.0, -1.0),
             A vertex at (-0.4, -1.0, -1.0),
             A vertex at (-1.0, -0.4, -1.0),
             A vertex at (-1.0, -1.0, -0.4))
            sage: edge_trunc = Cube.face_truncation(Cube.faces(1)[11])
            sage: edge_trunc.f_vector()
            (1, 10, 15, 7, 1)
            sage: tuple(f.ambient_V_indices() for f in edge_trunc.faces(2))
            ((0, 5, 6, 7),
             (1, 4, 5, 6, 8),
             (6, 7, 8, 9),
             (0, 2, 3, 7, 9),
             (1, 2, 8, 9),
             (0, 3, 4, 5),
             (1, 2, 3, 4))
             sage: face_trunc = Cube.face_truncation(Cube.faces(2)[2])
             sage: face_trunc.vertices()
             (A vertex at (1, -1, -1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (-1/3, -1, 1),
              A vertex at (-1/3, 1, 1),
              A vertex at (-1/3, 1, -1),
              A vertex at (-1/3, -1, -1))
             sage: face_trunc.face_lattice().is_isomorphic(Cube.face_lattice())
             True

        TESTS:

        Testing that the backend is preserved::

            sage: Cube = polytopes.cube(backend='field')
            sage: face_trunc = Cube.face_truncation(Cube.faces(2)[0])
            sage: face_trunc.backend()
            'field'

        Testing that :trac:`28506` is fixed::

            sage: P = polytopes.twenty_four_cell()
            sage: P = P.dilation(6)
            sage: P = P.change_ring(ZZ)
            sage: P.face_truncation(P.faces(2)[0], cut_frac=1)
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 27 vertices
        """
        if cut_frac is None:
            cut_frac = ZZ.one() / 3

        face_vertices = face.vertices()

        normal_vectors = []

        for facet in self.Hrepresentation():
            if all(facet.contains(x) and not facet.interior_contains(x)
                   for x in face_vertices):
                # The facet contains the face
                normal_vectors.append(facet.A())

        if linear_coefficients is not None:
            normal_vector = sum(linear_coefficients[i]*normal_vectors[i]
                                for i in range(len(normal_vectors)))
        else:
            normal_vector = sum(normal_vectors)

        B = - normal_vector * (face_vertices[0].vector())

        linear_evaluation = set(-normal_vector * (v.vector()) for v in self.vertices())

        if B == max(linear_evaluation):
            C = max(linear_evaluation.difference(set([B])))
        else:
            C = min(linear_evaluation.difference(set([B])))

        cut_height = (1 - cut_frac) * B + cut_frac * C
        ineq_vector = tuple([cut_height]) + tuple(normal_vector)

        new_ieqs = self.inequalities_list() + [ineq_vector]
        new_eqns = self.equations_list()

        # Some vertices might need fractions.
        parent = self.parent().base_extend(cut_frac/1)
        return parent.element_class(parent, None, [new_ieqs, new_eqns])

    def stack(self, face, position=None):
        r"""
        Return a new polyhedron formed by stacking onto a ``face``. Stacking a
        face adds a new vertex located slightly outside of the designated face.

        INPUT:

        - ``face`` -- a PolyhedronFace

        - ``position`` -- a positive number. Determines a relative distance
          from the barycenter of ``face``. A value close to 0 will place the
          new vertex close to the face and a large value further away. Default
          is `1`. If the given value is too large, an error is returned.

        OUTPUT:

        A Polyhedron object

        EXAMPLES::

            sage: cube = polytopes.cube()
            sage: square_face = cube.facets()[2]
            sage: stacked_square = cube.stack(square_face)
            sage: stacked_square.f_vector()
            (1, 9, 16, 9, 1)

            sage: edge_face = cube.faces(1)[3]
            sage: stacked_edge = cube.stack(edge_face)
            sage: stacked_edge.f_vector()
            (1, 9, 17, 10, 1)

            sage: cube.stack(cube.faces(0)[0])
            Traceback (most recent call last):
            ...
            ValueError: cannot stack onto a vertex

            sage: stacked_square_half = cube.stack(square_face,position=1/2)
            sage: stacked_square_half.f_vector()
            (1, 9, 16, 9, 1)
            sage: stacked_square_large = cube.stack(square_face,position=10)

            sage: hexaprism = polytopes.regular_polygon(6).prism()
            sage: hexaprism.f_vector()
            (1, 12, 18, 8, 1)
            sage: square_face = hexaprism.faces(2)[2]
            sage: stacked_hexaprism = hexaprism.stack(square_face)
            sage: stacked_hexaprism.f_vector()
            (1, 13, 22, 11, 1)

            sage: hexaprism.stack(square_face,position=4)
            Traceback (most recent call last):
            ...
            ValueError: the chosen position is too large

            sage: s = polytopes.simplex(7)
            sage: f = s.faces(3)[69]
            sage: sf = s.stack(f); sf
            A 7-dimensional polyhedron in QQ^8 defined as the convex hull of 9 vertices
            sage: sf.vertices()
            (A vertex at (-4, -4, -4, -4, 17/4, 17/4, 17/4, 17/4),
             A vertex at (0, 0, 0, 0, 0, 0, 0, 1),
             A vertex at (0, 0, 0, 0, 0, 0, 1, 0),
             A vertex at (0, 0, 0, 0, 0, 1, 0, 0),
             A vertex at (0, 0, 0, 0, 1, 0, 0, 0),
             A vertex at (0, 0, 0, 1, 0, 0, 0, 0),
             A vertex at (0, 0, 1, 0, 0, 0, 0, 0),
             A vertex at (0, 1, 0, 0, 0, 0, 0, 0),
             A vertex at (1, 0, 0, 0, 0, 0, 0, 0))

        It is possible to stack on unbounded faces::

            sage: Q = Polyhedron(vertices=[[0,1],[1,0]],rays=[[1,1]])
            sage: E = Q.faces(1)
            sage: Q.stack(E[0],1/2).Vrepresentation()
            (A vertex at (0, 1),
             A vertex at (1, 0),
             A ray in the direction (1, 1),
             A vertex at (2, 0))
            sage: Q.stack(E[1],1/2).Vrepresentation()
            (A vertex at (0, 1),
             A vertex at (0, 2),
             A vertex at (1, 0),
             A ray in the direction (1, 1))
            sage: Q.stack(E[2],1/2).Vrepresentation()
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (1, 0),
             A ray in the direction (1, 1))

        Stacking requires a proper face::

            sage: Q.stack(Q.faces(2)[0])
            Traceback (most recent call last):
            ...
            ValueError: can only stack on proper face

        TESTS:

        Checking that the backend is preserved::

            sage: Cube = polytopes.cube(backend='field')
            sage: stack = Cube.stack(Cube.faces(2)[0])
            sage: stack.backend()
            'field'

        Taking the stacking vertex too far with the parameter ``position``
        may result in a failure to produce the desired
        (combinatorial type of) polytope.
        The interval of permitted values is always open.
        This is the smallest unpermitted value::

            sage: P = polytopes.octahedron()
            sage: P.stack(P.faces(2)[0], position=4)
            Traceback (most recent call last):
            ...
            ValueError: the chosen position is too large

        Testing that :trac:`29057` is fixed::

            sage: P = polytopes.cross_polytope(4)
            sage: P.stack(P.faces(3)[0])
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 9 vertices
        """
        from sage.geometry.polyhedron.face import PolyhedronFace
        if not isinstance(face, PolyhedronFace):
            raise TypeError("{} should be a PolyhedronFace of {}".format(face, self))
        elif face.dim() == 0:
            raise ValueError("cannot stack onto a vertex")
        elif face.dim() == -1 or face.dim() == self.dim():
            raise ValueError("can only stack on proper face")
        if position is None:
            position = 1

        barycenter = ZZ.one()*sum([v.vector() for v in face.vertices()]) / len(face.vertices())
        locus_polyhedron = face.stacking_locus()
        repr_point = locus_polyhedron.representative_point()
        new_vertex = (1-position)*barycenter + position*repr_point
        if not locus_polyhedron.relative_interior_contains(new_vertex):
            raise ValueError("the chosen position is too large")

        parent = self.parent().base_extend(new_vertex)
        return parent.element_class(parent, [self.vertices() + (new_vertex,), self.rays(), self.lines()], None)

    def wedge(self, face, width=1):
        r"""
        Return the wedge over a ``face`` of the polytope ``self``.

        The wedge over a face `F` of a polytope `P` with width `w \not= 0`
        is defined as:

        .. MATH::

            (P \times \mathbb{R}) \cap \{a^\top x + |w x_{d+1}| \leq b\}

        where `\{x | a^\top x = b\}` is a supporting hyperplane defining `F`.

        INPUT:

        - ``face`` -- a PolyhedronFace of ``self``, the face which we take
          the wedge over
        - ``width`` -- a nonzero number (default: ``1``);
          specifies how wide the wedge will be

        OUTPUT:

        A (bounded) polyhedron

        EXAMPLES::

            sage: P_4 = polytopes.regular_polygon(4)
            sage: W1 = P_4.wedge(P_4.faces(1)[0]); W1
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 6 vertices
            sage: triangular_prism = polytopes.regular_polygon(3).prism()
            sage: W1.is_combinatorially_isomorphic(triangular_prism)
            True

            sage: Q = polytopes.hypersimplex(4,2)
            sage: W2 = Q.wedge(Q.faces(2)[7]); W2
            A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 9 vertices
            sage: W2.vertices()
            (A vertex at (1, 1, 0, 0, 1),
             A vertex at (1, 1, 0, 0, -1),
             A vertex at (1, 0, 1, 0, 1),
             A vertex at (1, 0, 1, 0, -1),
             A vertex at (1, 0, 0, 1, 1),
             A vertex at (1, 0, 0, 1, -1),
             A vertex at (0, 0, 1, 1, 0),
             A vertex at (0, 1, 1, 0, 0),
             A vertex at (0, 1, 0, 1, 0))

            sage: W3 = Q.wedge(Q.faces(1)[11]); W3
            A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 10 vertices
            sage: W3.vertices()
            (A vertex at (1, 1, 0, 0, -2),
             A vertex at (1, 1, 0, 0, 2),
             A vertex at (1, 0, 1, 0, -2),
             A vertex at (1, 0, 1, 0, 2),
             A vertex at (1, 0, 0, 1, 1),
             A vertex at (1, 0, 0, 1, -1),
             A vertex at (0, 1, 0, 1, 0),
             A vertex at (0, 1, 1, 0, 1),
             A vertex at (0, 0, 1, 1, 0),
             A vertex at (0, 1, 1, 0, -1))

            sage: C_3_7 = polytopes.cyclic_polytope(3,7)
            sage: P_6 = polytopes.regular_polygon(6)
            sage: W4 = P_6.wedge(P_6.faces(1)[0])
            sage: W4.is_combinatorially_isomorphic(C_3_7.polar())
            True

        REFERENCES:

        For more information, see Chapter 15 of [HoDaCG17]_.

        TESTS:

        The backend should be preserved as long as the value of width permits.
        The base_ring will change to the field of fractions of the current
        base_ring, unless width forces a different ring. ::

            sage: P = polytopes.cyclic_polytope(3,7, base_ring=ZZ, backend='field')
            sage: W1 = P.wedge(P.faces(2)[0]); W1.base_ring(); W1.backend()
            Rational Field
            'field'
            sage: W2 = P.wedge(P.faces(2)[0], width=5/2); W2.base_ring(); W2.backend()
            Rational Field
            'field'
            sage: W2 = P.wedge(P.faces(2)[9], width=4/2); W2.base_ring(); W2.backend()
            Rational Field
            'field'
            sage: W2.vertices()
            (A vertex at (3, 9, 27, -1/2),
             A vertex at (4, 16, 64, -2),
             A vertex at (6, 36, 216, -10),
             A vertex at (5, 25, 125, -5),
             A vertex at (2, 4, 8, 0),
             A vertex at (1, 1, 1, 0),
             A vertex at (0, 0, 0, 0),
             A vertex at (3, 9, 27, 1/2),
             A vertex at (4, 16, 64, 2),
             A vertex at (6, 36, 216, 10),
             A vertex at (5, 25, 125, 5))
            sage: W2 = P.wedge(P.faces(2)[2], width=1.0); W2.base_ring(); W2.backend()
            Real Double Field
            'cdd'
        """
        width = width*ZZ.one()

        if not self.is_compact():
            raise ValueError("polyhedron 'self' must be a polytope")

        if width == 0:
            raise ValueError("the width should be nonzero")

        from sage.geometry.polyhedron.face import PolyhedronFace
        if not isinstance(face, PolyhedronFace):
            raise TypeError("{} should be a PolyhedronFace of {}".format(face, self))

        F_Hrep = vector([0]*(self.ambient_dim()+1))
        for facet in face.ambient_Hrepresentation():
            if facet.is_inequality():
                F_Hrep = F_Hrep + facet.vector()
        F_Hrep = list(F_Hrep)

        # Preserve the backend, if value of ``width`` permits.
        backend = None
        from .parent import does_backend_handle_base_ring
        if does_backend_handle_base_ring(width.base_ring().fraction_field(), self.backend()):
            backend = self.backend()

        L = Polyhedron(lines=[[1]])
        Q = self.product(L)
        ieqs = [F_Hrep + [width], F_Hrep + [-width]]
        H = Polyhedron(ieqs=ieqs, backend=backend)
        return Q.intersection(H)

    def lawrence_extension(self, v):
        """
        Return the Lawrence extension of ``self`` on the point ``v``.

        Let `P` be a polytope and `v` be a vertex of `P` or a point outside
        `P`. The Lawrence extension of `P` on `v` is the convex hull of
        `(v,1),(v,2)` and `(u,0)` for all vertices `u` in `P` other than `v`
        if `v` is a vertex.

        INPUT:
            - ``v`` -- a vertex of ``self`` or a point outside it

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.lawrence_extension(P.vertices()[0])
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 9 vertices
            sage: P.lawrence_extension([-1,-1,-1])
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 9 vertices

        REFERENCES:

            For more information, see Section 6.6 of [Zie2007]_.
        """
        if not self.is_compact():
            raise NotImplementedError("self must be a polytope")

        V = self.vertices_list()
        v = list(v)

        if self.contains(v) and (v not in V):
            raise ValueError("{} must not be a vertex or outside self".format(v))

        lambda_V = [u + [0] for u in V if u != v] + [v+[1]] + [v+[2]]
        parent = self.parent().base_extend(vector(v), ambient_dim=self.ambient_dim() + 1)
        return parent.element_class(parent, [lambda_V, [], []], None)

    def lawrence_polytope(self):
        r"""
        Return the Lawrence polytope of ``self``.

        Let `P` be a `d`-polytope in `\RR^r` with `n` vertices. The Lawrence
        polytope of `P` is the polytope whose vertices are the columns of the
        following `(r+n)`-by-`2n` matrix.

        .. MATH::

            \begin{pmatrix}
             V      &   V    \\
             I_n    &   2I_n
            \end{pmatrix},

        where `V` is the `r`-by-`n` vertices matrix of `P`.

        EXAMPLES::

            sage: P = polytopes.octahedron()
            sage: L = P.lawrence_polytope(); L
            A 9-dimensional polyhedron in ZZ^9 defined as the convex hull of 12 vertices
            sage: V = P.vertices_list()
            sage: i = 0
            sage: for v in V:
            ....:     v = v + i*[0]
            ....:     P = P.lawrence_extension(v)
            ....:     i = i + 1
            sage: P == L
            True

        REFERENCES:

            For more information, see Section 6.6 of [Zie2007]_.
        """
        from sage.matrix.constructor import block_matrix

        if not self.is_compact():
            raise NotImplementedError("self must be a polytope")

        V = self.vertices_matrix().transpose()
        n = self.n_vertices()
        I_n = matrix.identity(n)
        lambda_V = block_matrix([[V, I_n], [V, 2*I_n]])
        parent = self.parent().change_ring(self.base_ring(), ambient_dim=self.ambient_dim() + n)
        return parent.element_class(parent, [lambda_V, [], []], None)

    def is_lawrence_polytope(self):
        """
        Return ``True`` if ``self`` is a Lawrence polytope.

        A polytope is called a Lawrence polytope if it has a centrally
        symmetric (normalized) Gale diagram.

        EXAMPLES::

            sage: P = polytopes.hypersimplex(5,2)
            sage: L = P.lawrence_polytope()
            sage: L.is_lattice_polytope()
            True
            sage: egyptian_pyramid = polytopes.regular_polygon(4).pyramid()
            sage: egyptian_pyramid.is_lawrence_polytope()
            True
            sage: polytopes.octahedron().is_lawrence_polytope()
            False

        REFERENCES:

            For more information, see [BaSt1990]_.
        """
        if not self.is_compact():
            raise NotImplementedError("self must be a polytope")

        return self.combinatorial_polyhedron().is_lawrence_polytope()

    def _test_lawrence(self, tester=None, **options):
        """
        Run tests on the methods related to lawrence extensions.

        TESTS:

        Check that :trac:`28725` is fixed::

            sage: polytopes.regular_polygon(3)._test_lawrence()

        Check that :trac:`30293` is fixed::

            sage: polytopes.cube()._test_lawrence()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.backend() == 'normaliz' and not self.base_ring() in (ZZ, QQ):
            # Speeds up the doctest for significantly.
            self = self.change_ring(self._normaliz_field)

        if not self.is_compact():
            with tester.assertRaises(NotImplementedError):
                self.lawrence_polytope()
            with tester.assertRaises(NotImplementedError):
                self.lawrence_extension(self.vertices()[0])
            return

        if self.n_vertices() > 1:
            # ``v`` must be a vertex or outside ``self``.
            with tester.assertRaises(ValueError):
                self.lawrence_extension(self.center())

        if self.n_vertices() >= 40 or self.n_facets() > 40:
            # Avoid very long tests.
            return

        if self.n_vertices():
            from sage.misc.prandom import randint
            v = self.vertices()[randint(0, self.n_vertices()-1)].vector()

            # A lawrence extension with a vertex.
            P = self.lawrence_extension(v)
            tester.assertEqual(self.dim() + 1, P.dim())
            tester.assertEqual(self.n_vertices() + 1, P.n_vertices())
            tester.assertEqual(self.backend(), P.backend())

            if self.n_vertices() > 1:
                # A lawrence extension with a point outside of the polyhedron.
                Q = self.lawrence_extension(2*v - self.center())
                tester.assertEqual(self.dim() + 1, Q.dim())
                tester.assertEqual(self.n_vertices() + 2, Q.n_vertices())
                tester.assertEqual(self.backend(), Q.backend())  # Any backend should handle the fraction field.

                import warnings

                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    try:
                        # Implicitly checks :trac:`30328`.
                        R = self.lawrence_extension(2.0*v - self.center())
                        tester.assertEqual(self.dim() + 1, R.dim())
                        tester.assertEqual(self.n_vertices() + 2, R.n_vertices())

                        tester.assertTrue(Q.is_combinatorially_isomorphic(R))
                    except UserWarning:
                        # Data is numerically complicated.
                        pass
                    except ValueError as err:
                        if "Numerical inconsistency" not in err.args[0]:
                            raise err

        if self.n_vertices() >= 12 or (self.base_ring() not in (ZZ, QQ) and self.backend() == 'field'):
            # Avoid very long tests.
            return

        P = self.lawrence_polytope()
        tester.assertEqual(self.dim() + self.n_vertices(), P.dim())
        tester.assertEqual(self.n_vertices()*2, P.n_vertices())
        tester.assertEqual(self.backend(), P.backend())
        tester.assertTrue(P.is_lawrence_polytope())

        # Construct the lawrence polytope iteratively by lawrence extensions.
        V = self.vertices_list()
        Q = self
        i = 0
        for v in V:
            v = v + i*[0]
            Q = Q.lawrence_extension(v)
            i = i + 1
        tester.assertEqual(P, Q)

    def barycentric_subdivision(self, subdivision_frac=None):
        r"""
        Return the barycentric subdivision of a compact polyhedron.

        DEFINITION:

        The barycentric subdivision of a compact polyhedron is a standard way
        to triangulate its faces in such a way that maximal faces correspond to
        flags of faces of the starting polyhedron (i.e. a maximal chain in the
        face lattice of the polyhedron). As a simplicial complex, this is known
        as the order complex of the face lattice of the polyhedron.

        REFERENCE:

        See :wikipedia:`Barycentric_subdivision`
        Section 6.6, Handbook of Convex Geometry, Volume A, edited by P.M. Gruber and J.M.
        Wills. 1993, North-Holland Publishing Co..

        INPUT:

        - ``subdivision_frac`` -- number. Gives the proportion how far the new
          vertices are pulled out of the polytope. Default is `\frac{1}{3}` and
          the value should be smaller than `\frac{1}{2}`. The subdivision is
          computed on the polar polyhedron.

        OUTPUT:

        A Polyhedron object, subdivided as described above.

        EXAMPLES::

            sage: P = polytopes.hypercube(3)
            sage: P.barycentric_subdivision()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull
            of 26 vertices
            sage: P = Polyhedron(vertices=[[0,0,0],[0,1,0],[1,0,0],[0,0,1]])
            sage: P.barycentric_subdivision()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull
            of 14 vertices
            sage: P = Polyhedron(vertices=[[0,1,0],[0,0,1],[1,0,0]])
            sage: P.barycentric_subdivision()
            A 2-dimensional polyhedron in QQ^3 defined as the convex hull
            of 6 vertices
            sage: P = polytopes.regular_polygon(4, base_ring=QQ)
            sage: P.barycentric_subdivision()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 8
            vertices

        TESTS::

            sage: P.barycentric_subdivision(1/2)
            Traceback (most recent call last):
            ...
            ValueError: the subdivision fraction should be between 0 and 1/2
            sage: P = Polyhedron(ieqs=[[1,0,1],[0,1,0],[1,0,0],[0,0,1]])
            sage: P.barycentric_subdivision()
            Traceback (most recent call last):
            ...
            ValueError: the polytope has to be compact
            sage: P = Polyhedron(vertices=[[0,0,0],[0,1,0],[1,0,0],[0,0,1]], backend='field')
            sage: P.barycentric_subdivision()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 14 vertices

            sage: polytopes.simplex(backend='field').barycentric_subdivision().backend()
            'field'
            sage: polytopes.cube(backend='cdd').barycentric_subdivision().backend()
            'cdd'
        """
        if subdivision_frac is None:
            subdivision_frac = ZZ.one() / 3

        if not self.is_compact():
            raise ValueError("the polytope has to be compact")
        if not (0 < subdivision_frac < ZZ.one() / 2):
            raise ValueError("the subdivision fraction should be "
                             "between 0 and 1/2")

        barycenter = self.center()
        parent = self.parent().base_extend(subdivision_frac)

        start_polar = (self - barycenter).polar(in_affine_span=True)
        polar = (self - barycenter).polar(in_affine_span=True)

        for i in range(self.dimension() - 1):

            new_ineq = []
            subdivided_faces = list(start_polar.faces(i))
            Hrep = polar.Hrepresentation()

            for face in subdivided_faces:

                face_vertices = face.vertices()
                normal_vectors = []

                for facet in Hrep:
                    if all(facet.contains(v) and not facet.interior_contains(v)
                           for v in face_vertices):
                        # The facet contains the face
                        normal_vectors.append(facet.A())

                normal_vector = sum(normal_vectors)
                B = - normal_vector * (face_vertices[0].vector())
                linear_evaluation = set([-normal_vector * (v.vector())
                                         for v in polar.vertices()])

                if B == max(linear_evaluation):
                    C = max(linear_evaluation.difference(set([B])))
                else:
                    C = min(linear_evaluation.difference(set([B])))

                ineq_vector = [(1 - subdivision_frac) * B + subdivision_frac * C] + list(normal_vector)
                new_ineq += [ineq_vector]

            new_ieqs = polar.inequalities_list() + new_ineq
            new_eqns = polar.equations_list()

            polar = parent.element_class(parent, None, [new_ieqs, new_eqns])

        return (polar.polar(in_affine_span=True)) + barycenter

    def face_lattice(self):
        """
        Return the face-lattice poset.

        OUTPUT:

        A :class:`~sage.combinat.posets.posets.FinitePoset`. Elements
        are given as
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`.

        In the case of a full-dimensional polytope, the faces are
        pairs (vertices, inequalities) of the spanning vertices and
        corresponding saturated inequalities. In general, a face is
        defined by a pair (V-rep. objects, H-rep. objects). The
        V-representation objects span the face, and the corresponding
        H-representation objects are those inequalities and equations
        that are saturated on the face.

        The bottom-most element of the face lattice is the "empty
        face". It contains no V-representation object. All
        H-representation objects are incident.

        The top-most element is the "full face". It is spanned by all
        V-representation objects. The incident H-representation
        objects are all equations and no inequalities.

        In the case of a full-dimensional polytope, the "empty face"
        and the "full face" are the empty set (no vertices, all
        inequalities) and the full polytope (all vertices, no
        inequalities), respectively.

        ALGORITHM:

        See :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator`.

        .. NOTE::

            The face lattice is not cached, as long as this creates a memory leak, see :trac:`28982`.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: fl = square.face_lattice();fl
            Finite lattice containing 10 elements
            sage: list(f.ambient_V_indices() for f in fl)
            [(), (0,), (1,), (0, 1), (2,), (1, 2), (3,), (0, 3), (2, 3), (0, 1, 2, 3)]
            sage: poset_element = fl[5]
            sage: a_face = poset_element
            sage: a_face
            A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: a_face.ambient_V_indices()
            (1, 2)
            sage: set(a_face.ambient_Vrepresentation()) == \
            ....: set([square.Vrepresentation(1), square.Vrepresentation(2)])
            True
            sage: a_face.ambient_Vrepresentation()
            (A vertex at (1, 1), A vertex at (-1, 1))
            sage: a_face.ambient_Hrepresentation()
            (An inequality (0, -1) x + 1 >= 0,)

        A more complicated example::

            sage: c5_10 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5] for i in range(1,11)])
            sage: c5_10_fl = c5_10.face_lattice()
            sage: [len(x) for x in c5_10_fl.level_sets()]
            [1, 10, 45, 100, 105, 42, 1]

        Note that if the polyhedron contains lines then there is a
        dimension gap between the empty face and the first non-empty
        face in the face lattice::

            sage: line = Polyhedron(vertices=[(0,)], lines=[(1,)])
            sage: [ fl.dim() for fl in line.face_lattice() ]
            [-1, 1]

        TESTS::

            sage: c5_20 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5]
            ....:     for i in range(1,21)])
            sage: c5_20_fl = c5_20.face_lattice() # long time
            sage: [len(x) for x in c5_20_fl.level_sets()] # long time
            [1, 20, 190, 580, 680, 272, 1]
            sage: polytopes.hypercube(2).face_lattice().plot()
            Graphics object consisting of 27 graphics primitives
            sage: level_sets = polytopes.cross_polytope(2).face_lattice().level_sets()
            sage: level_sets[0][0].ambient_V_indices(), level_sets[-1][0].ambient_V_indices()
            ((), (0, 1, 2, 3))

        Various degenerate polyhedra::

            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0]]).face_lattice().level_sets()]
            [[()], [(0,), (1,), (2,)], [(0, 1), (0, 2), (1, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(1,0,0),(0,1,0)], rays=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(1,), (2,)], [(0, 1), (0, 2), (1, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(rays=[(1,0,0),(0,1,0)], vertices=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(0,)], [(0, 1), (0, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(rays=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0,)], [(0, 1), (0, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(1,),(0,)]).face_lattice().level_sets()]
            [[()], [(0,), (1,)], [(0, 1)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(1,0,0),(0,1,0)], lines=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(0, 1), (0, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0,0)], vertices=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(0, 1)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0)], rays=[(0,1)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0, 1)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(0,)], lines=[(1,)]).face_lattice().level_sets()]
            [[()], [(0, 1)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0, 1)]]

        Test that computing the face lattice does not lead to a memory leak::

            sage: import gc
            sage: _ = gc.collect()
            sage: P = polytopes.cube()
            sage: a = P.face_lattice()
            sage: n = get_memory_usage()
            sage: P = polytopes.cube()
            sage: a = P.face_lattice()
            sage: _ = gc.collect()
            sage: n == get_memory_usage()
            True
        """
        from sage.combinat.posets.lattices import FiniteLatticePoset
        return FiniteLatticePoset(self.hasse_diagram())

    @cached_method
    def hasse_diagram(self):
        r"""
        Return the Hasse diagram of the face lattice of ``self``.

        This is the Hasse diagram of the poset of the faces of ``self``.

        OUTPUT: a directed graph

        EXAMPLES::

            sage: P = polytopes.regular_polygon(4).pyramid()
            sage: D = P.hasse_diagram(); D
            Digraph on 20 vertices
            sage: D.degree_polynomial()
            x^5 + x^4*y + x*y^4 + y^5 + 4*x^3*y + 8*x^2*y^2 + 4*x*y^3

        Faces of an mutable polyhedron are not hashable. Hence those are not suitable as
        vertices of the hasse diagram. Use the combinatorial polyhedron instead::

            sage: P = polytopes.regular_polygon(4).pyramid()
            sage: parent = P.parent()
            sage: parent = parent.change_ring(QQ, backend='ppl')
            sage: Q = parent._element_constructor_(P, mutable=True)
            sage: Q.hasse_diagram()
            Traceback (most recent call last):
            ...
            TypeError: mutable polyhedra are unhashable
            sage: C = Q.combinatorial_polyhedron()
            sage: D = C.hasse_diagram()
            sage: set(D.vertices()) == set(range(20))
            True
            sage: def index_to_combinatorial_face(n):
            ....:     return C.face_by_face_lattice_index(n)
            sage: D.relabel(index_to_combinatorial_face, inplace=True)
            sage: D.vertices()
            [A -1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 3-dimensional face of a 3-dimensional combinatorial polyhedron]
            sage: D.degree_polynomial()
            x^5 + x^4*y + x*y^4 + y^5 + 4*x^3*y + 8*x^2*y^2 + 4*x*y^3
        """

        from sage.geometry.polyhedron.face import combinatorial_face_to_polyhedral_face
        C = self.combinatorial_polyhedron()
        D = C.hasse_diagram()

        def index_to_polyhedron_face(n):
            return combinatorial_face_to_polyhedral_face(
                    self, C.face_by_face_lattice_index(n))

        return D.relabel(index_to_polyhedron_face, inplace=False, immutable=True)

    def face_generator(self, face_dimension=None, dual=None):
        r"""
        Return an iterator over the faces of given dimension.

        If dimension is not specified return an iterator over all faces.

        INPUT:

        - ``face_dimension`` -- integer (default ``None``),
          yield only faces of this dimension if specified
        - ``dual`` -- boolean (default ``None``);
          if ``True``, generate the faces using the vertices;
          if ``False``, generate the faces using the facets;
          if ``None``, pick automatically

        OUTPUT:

        A :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator_geom`.
        This class iterates over faces as
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`. See
        :mod:`~sage.geometry.polyhedron.face` for details. The order
        is random but fixed.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: it = P.face_generator()
            sage: it
            Iterator over the faces of a 3-dimensional polyhedron in ZZ^3
            sage: list(it)
            [A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices,
             A -1-dimensional face of a Polyhedron in ZZ^3,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices]

             sage: P = polytopes.hypercube(4)
             sage: list(P.face_generator(2))[:4]
             [A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices,
              A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices,
              A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices,
              A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices]

        If a polytope has more facets than vertices, the dual mode is chosen::

            sage: P = polytopes.cross_polytope(3)
            sage: list(P.face_generator())
            [A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 6 vertices,
             A -1-dimensional face of a Polyhedron in ZZ^3,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices]

        The face iterator can also be slightly modified.
        In non-dual mode we can skip subfaces of the current (proper) face::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(dual=False)
            sage: _ = next(it), next(it)
            sage: face = next(it)
            sage: face.ambient_H_indices()
            (5,)
            sage: it.ignore_subfaces()
            sage: face = next(it)
            sage: face.ambient_H_indices()
            (4,)
            sage: it.ignore_subfaces()
            sage: [face.ambient_H_indices() for face in it]
            [(3,),
             (2,),
             (1,),
             (0,),
             (2, 3),
             (1, 3),
             (1, 2, 3),
             (1, 2),
             (0, 2),
             (0, 1, 2),
             (0, 1)]

        In dual mode we can skip supfaces of the current (proper) face::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(dual=True)
            sage: _ = next(it), next(it)
            sage: face = next(it)
            sage: face.ambient_V_indices()
            (7,)
            sage: it.ignore_supfaces()
            sage: next(it)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: face = next(it)
            sage: face.ambient_V_indices()
            (5,)
            sage: it.ignore_supfaces()
            sage: [face.ambient_V_indices() for face in it]
            [(4,),
             (3,),
             (2,),
             (1,),
             (0,),
             (1, 6),
             (3, 4),
             (2, 3),
             (0, 3),
             (0, 1, 2, 3),
             (1, 2),
             (0, 1)]

        In non-dual mode, we cannot skip supfaces::

            sage: it = P.face_generator(dual=False)
            sage: _ = next(it), next(it)
            sage: next(it)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: it.ignore_supfaces()
            Traceback (most recent call last):
            ...
            ValueError: only possible when in dual mode

        In dual mode, we cannot skip subfaces::

            sage: it = P.face_generator(dual=True)
            sage: _ = next(it), next(it)
            sage: next(it)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: it.ignore_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: only possible when not in dual mode

        We can only skip sub-/supfaces of proper faces::

            sage: it = P.face_generator(dual=False)
            sage: next(it)
            A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: it.ignore_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: iterator not set to a face yet

        .. SEEALSO::

            :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator_geom`.

        ALGORITHM:

        See :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`.

        TESTS::

            sage: P = polytopes.simplex()
            sage: list(P.face_generator(-2))
            []
            sage: list(P.face_generator(-1))
            [A -1-dimensional face of a Polyhedron in ZZ^4]
            sage: list(P.face_generator(3))
            [A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices]

            sage: list(Polyhedron().face_generator())
            [A -1-dimensional face of a Polyhedron in ZZ^0]

        Check that :trac:`29155` is fixed::

            sage: P = polytopes.permutahedron(3)
            sage: [f] = P.face_generator(2)
            sage: f.ambient_Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0,)
        """
        from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator import FaceIterator_geom
        return FaceIterator_geom(self, output_dimension=face_dimension, dual=dual)

    def faces(self, face_dimension):
        """
        Return the faces of given dimension

        INPUT:

        - ``face_dimension`` -- integer. The dimension of the faces
          whose representation will be returned.

        OUTPUT:

        A tuple of
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`. See
        :mod:`~sage.geometry.polyhedron.face` for details. The order
        is random but fixed.

        .. SEEALSO::

            :meth:`face_generator`,
            :meth:`facet`.

        EXAMPLES:

        Here we find the vertex and face indices of the eight three-dimensional
        facets of the four-dimensional hypercube::

            sage: p = polytopes.hypercube(4)
            sage: list(f.ambient_V_indices() for f in p.faces(3))
            [(0, 5, 6, 7, 8, 9, 14, 15),
             (1, 4, 5, 6, 10, 13, 14, 15),
             (1, 2, 6, 7, 8, 10, 11, 15),
             (8, 9, 10, 11, 12, 13, 14, 15),
             (0, 3, 4, 5, 9, 12, 13, 14),
             (0, 2, 3, 7, 8, 9, 11, 12),
             (1, 2, 3, 4, 10, 11, 12, 13),
             (0, 1, 2, 3, 4, 5, 6, 7)]

            sage: face = p.faces(3)[3]
            sage: face.ambient_Hrepresentation()
            (An inequality (1, 0, 0, 0) x + 1 >= 0,)
            sage: face.vertices()
            (A vertex at (-1, -1, 1, -1),
             A vertex at (-1, -1, 1, 1),
             A vertex at (-1, 1, -1, -1),
             A vertex at (-1, 1, 1, -1),
             A vertex at (-1, 1, 1, 1),
             A vertex at (-1, 1, -1, 1),
             A vertex at (-1, -1, -1, 1),
             A vertex at (-1, -1, -1, -1))

        You can use the
        :meth:`~sage.geometry.polyhedron.representation.PolyhedronRepresentation.index`
        method to enumerate vertices and inequalities::

            sage: def get_idx(rep): return rep.index()
            sage: [get_idx(_) for _ in face.ambient_Hrepresentation()]
            [4]
            sage: [get_idx(_) for _ in face.ambient_Vrepresentation()]
            [8, 9, 10, 11, 12, 13, 14, 15]

            sage: [ ([get_idx(_) for _ in face.ambient_Vrepresentation()],
            ....:    [get_idx(_) for _ in face.ambient_Hrepresentation()])
            ....:   for face in p.faces(3) ]
            [([0, 5, 6, 7, 8, 9, 14, 15], [7]),
             ([1, 4, 5, 6, 10, 13, 14, 15], [6]),
             ([1, 2, 6, 7, 8, 10, 11, 15], [5]),
             ([8, 9, 10, 11, 12, 13, 14, 15], [4]),
             ([0, 3, 4, 5, 9, 12, 13, 14], [3]),
             ([0, 2, 3, 7, 8, 9, 11, 12], [2]),
             ([1, 2, 3, 4, 10, 11, 12, 13], [1]),
             ([0, 1, 2, 3, 4, 5, 6, 7], [0])]

        TESTS::

            sage: pr = Polyhedron(rays = [[1,0,0],[-1,0,0],[0,1,0]], vertices = [[-1,-1,-1]], lines=[(0,0,1)])
            sage: pr.faces(4)
            ()
            sage: pr.faces(3)[0].ambient_V_indices()
            (0, 1, 2, 3)
            sage: pr.facets()[0].ambient_V_indices()
            (0, 1, 2)
            sage: pr.faces(1)
            ()
            sage: pr.faces(0)
            ()
            sage: pr.faces(-1)
            (A -1-dimensional face of a Polyhedron in QQ^3,)
        """
        return tuple(self.face_generator(face_dimension))

    def facets(self):
        r"""
        Return the facets of the polyhedron.

        Facets are the maximal nontrivial faces of polyhedra.
        The empty face and the polyhedron itself are trivial.

        A facet of a `d`-dimensional polyhedron is a face of dimension
        `d-1`. For `d \neq 0` the converse is true as well.

        OUTPUT:

        A tuple of
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`. See
        :mod:`~sage.geometry.polyhedron.face` for details. The order
        is random but fixed.

        .. SEEALSO:: :meth:`facets`

        EXAMPLES:

        Here we find the eight three-dimensional facets of the
        four-dimensional hypercube::

            sage: p = polytopes.hypercube(4)
            sage: p.facets()
            (A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices)

        This is the same result as explicitly finding the
        three-dimensional faces::

            sage: dim = p.dimension()
            sage: p.faces(dim-1)
            (A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices)

        The ``0``-dimensional polyhedron does not have facets::

            sage: P = Polyhedron([[0]])
            sage: P.facets()
            ()
        """
        if self.dimension() == 0:
            return ()
        return self.faces(self.dimension()-1)

    def join_of_Vrep(self, *Vrepresentatives):
        r"""
        Return the smallest face that contains ``Vrepresentatives``.

        INPUT:

        - ``Vrepresentatives`` -- vertices/rays/lines of ``self`` or indices of such

        OUTPUT: a :class:`~sage.geometry.polyhedron.face.PolyhedronFace`

        .. NOTE::

            In the case of unbounded polyhedra, the join of rays etc. may not be well-defined.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: P.join_of_Vrep(1)
            A 0-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 1 vertex
            sage: P.join_of_Vrep()
            A -1-dimensional face of a Polyhedron in ZZ^5
            sage: P.join_of_Vrep(0,12,13).ambient_V_indices()
            (0, 12, 13, 68)

        The input is flexible::

            sage: P.join_of_Vrep(2, P.vertices()[3], P.Vrepresentation(4))
            A 2-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 6 vertices

        ::

            sage: P = polytopes.cube()
            sage: a, b = P.faces(0)[:2]
            sage: P.join_of_Vrep(a, b)
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices

        In the case of an unbounded polyhedron, the join may not be well-defined::

            sage: P = Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,1]])
            sage: P.join_of_Vrep(0)
            A 0-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex
            sage: P.join_of_Vrep(0,1)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 2 vertices
            sage: P.join_of_Vrep(0,2)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: P.join_of_Vrep(1,2)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: P.join_of_Vrep(2)
            Traceback (most recent call last):
            ...
            ValueError: the join is not well-defined

        The ``Vrepresentatives`` must be of ``self``::

            sage: P = polytopes.cube(backend='ppl')
            sage: Q = polytopes.cube(backend='field')
            sage: v = P.vertices()[0]
            sage: P.join_of_Vrep(v)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: Q.join_of_Vrep(v)
            Traceback (most recent call last):
            ...
            ValueError: not a Vrepresentative of ``self``
            sage: f = P.faces(0)[0]
            sage: P.join_of_Vrep(v)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: Q.join_of_Vrep(v)
            Traceback (most recent call last):
            ...
            ValueError: not a Vrepresentative of ``self``

        TESTS:

        ``least_common_superface_of_Vrep`` is an alias::

            sage: P.least_common_superface_of_Vrep(v)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: P.least_common_superface_of_Vrep == P.join_of_Vrep
            True

        Error message for invalid input::

            sage: P.join_of_Vrep('foo')
            Traceback (most recent call last):
            ...
            ValueError: foo is not a Vrepresentative
        """
        from sage.geometry.polyhedron.representation import Vrepresentation
        from sage.geometry.polyhedron.face import PolyhedronFace

        new_indices = [0]*len(Vrepresentatives)
        for i, v in enumerate(Vrepresentatives):
            if isinstance(v, PolyhedronFace) and v.dim() == 0:
                if v.polyhedron() is not self:
                    raise ValueError("not a Vrepresentative of ``self``")
                new_indices[i] = v.ambient_V_indices()[0]
            elif v in ZZ:
                new_indices[i] = v
            elif isinstance(v, Vrepresentation):
                if v.polyhedron() is not self:
                    raise ValueError("not a Vrepresentative of ``self``")
                new_indices[i] = v.index()
            else:
                raise ValueError("{} is not a Vrepresentative".format(v))

        return self.face_generator().join_of_Vrep(*new_indices)

    least_common_superface_of_Vrep = join_of_Vrep

    def meet_of_Hrep(self, *Hrepresentatives):
        r"""
        Return the largest face that is contained in ``Hrepresentatives``.

        INPUT:

        - ``Hrepresentatives`` -- facets or indices of Hrepresentatives;
          the indices are assumed to be the indices of the Hrepresentation

        OUTPUT: a :class:`~sage.geometry.polyhedron.face.PolyhedronFace`

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: P.meet_of_Hrep()
            A 4-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 120 vertices
            sage: P.meet_of_Hrep(1)
            A 3-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 24 vertices
            sage: P.meet_of_Hrep(4)
            A 3-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 12 vertices
            sage: P.meet_of_Hrep(1,3,7)
            A 1-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 2 vertices
            sage: P.meet_of_Hrep(1,3,7).ambient_H_indices()
            (0, 1, 3, 7)

        The indices are the indices of the Hrepresentation.
        ``0`` corresponds to an equation and is ignored::

            sage: P.meet_of_Hrep(0)
            A 4-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 120 vertices

        The input is flexible::

            sage: P.meet_of_Hrep(P.facets()[-1], P.inequalities()[2], 7)
            A 1-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 2 vertices

        The ``Hrepresentatives`` must belong to ``self``::

            sage: P = polytopes.cube(backend='ppl')
            sage: Q = polytopes.cube(backend='field')
            sage: f = P.facets()[0]
            sage: P.meet_of_Hrep(f)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: Q.meet_of_Hrep(f)
            Traceback (most recent call last):
            ...
            ValueError: not a facet of ``self``
            sage: f = P.inequalities()[0]
            sage: P.meet_of_Hrep(f)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: Q.meet_of_Hrep(f)
            Traceback (most recent call last):
            ...
            ValueError: not a facet of ``self``

        TESTS:

        Equations are not considered by the combinatorial polyhedron.
        We check that the index corresponds to the Hrepresentation index::

            sage: P = polytopes.permutahedron(3, backend='field')
            sage: P.Hrepresentation()
            (An inequality (0, 0, 1) x - 1 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (0, 1, 1) x - 3 >= 0,
             An inequality (1, 0, 0) x - 1 >= 0,
             An inequality (1, 0, 1) x - 3 >= 0,
             An inequality (1, 1, 0) x - 3 >= 0,
             An equation (1, 1, 1) x - 6 == 0)
            sage: P.meet_of_Hrep(0).ambient_Hrepresentation()
            (An inequality (0, 0, 1) x - 1 >= 0, An equation (1, 1, 1) x - 6 == 0)

            sage: P = polytopes.permutahedron(3, backend='ppl')
            sage: P.Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0,
             An inequality (1, 1, 0) x - 3 >= 0,
             An inequality (-1, -1, 0) x + 5 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (-1, 0, 0) x + 3 >= 0,
             An inequality (1, 0, 0) x - 1 >= 0,
             An inequality (0, -1, 0) x + 3 >= 0)
            sage: P.meet_of_Hrep(1).ambient_Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0, An inequality (1, 1, 0) x - 3 >= 0)

        ``greatest_common_subface_of_Hrep`` is an alias::

            sage: P.greatest_common_subface_of_Hrep(1).ambient_Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0, An inequality (1, 1, 0) x - 3 >= 0)
            sage: P.greatest_common_subface_of_Hrep == P.meet_of_Hrep
            True

        Error message for invalid input::

            sage: P.meet_of_Hrep('foo')
            Traceback (most recent call last):
            ...
            ValueError: foo is not a Hrepresentative
        """
        from sage.geometry.polyhedron.representation import Inequality, Equation
        from sage.geometry.polyhedron.face import PolyhedronFace

        # Equations are ignored by combinatorial polyhedron for indexing.
        offset = 0
        if self.n_equations() and self.Hrepresentation(0).is_equation():
            offset = self.n_equations()

        new_indices = []
        for i, facet in enumerate(Hrepresentatives):
            if isinstance(facet, PolyhedronFace) and facet.dim() + 1 == self.dim():
                if facet.polyhedron() is not self:
                    raise ValueError("not a facet of ``self``")
                H_indices = facet.ambient_H_indices()
                facet = H_indices[0] if H_indices[0] >= offset else H_indices[-1]

            if facet in ZZ and facet >= offset:
                # Note that ``CombinatorialPolyhedron`` ignores indices of equations
                # and has equations last.
                new_indices.append(facet - offset)
            elif isinstance(facet, Inequality):
                if facet.polyhedron() is not self:
                    raise ValueError("not a facet of ``self``")
                new_indices.append(facet.index() - offset)
            elif isinstance(facet, Equation):
                # Ignore equations.
                continue
            elif facet in ZZ and 0 <= facet < offset:
                # Ignore equations.
                continue
            else:
                raise ValueError("{} is not a Hrepresentative".format(facet))

        return self.face_generator().meet_of_Hrep(*new_indices)

    greatest_common_subface_of_Hrep = meet_of_Hrep

    def _test_combinatorial_face_as_combinatorial_polyhedron(self, tester=None, **options):
        """
        Run tests on obtaining the combinatorial face as combinatorial polyhedron.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_combinatorial_face_as_combinatorial_polyhedron()
        """
        if not self.is_compact():
            return
        if self.dim() > 7 or self.n_vertices() > 100 or self.n_facets() > 100:
            # Avoid very long tests.
            return
        if self.dim() < 1:
            # Avoid trivial cases.
            return

        from sage.misc.prandom import random

        if tester is None:
            tester = self._tester(**options)

        it = self.face_generator()
        _ = next(it), next(it)  # get rid of non_proper faces
        C1 = self.combinatorial_polyhedron()
        it1 = C1.face_iter()
        C2 = C1.dual()
        it2 = C2.face_iter(dual=not it1.dual)

        for f in it:
            f1 = next(it1)
            f2 = next(it2)
            if random() < 0.95:
                # Only test a random 5 percent of the faces.
                continue

            P = f.as_polyhedron()
            D1 = f1.as_combinatorial_polyhedron()
            D2 = f2.as_combinatorial_polyhedron(quotient=True).dual()
            D1._test_bitsets(tester, **options)
            D2._test_bitsets(tester, **options)
            tester.assertTrue(P.combinatorial_polyhedron().vertex_facet_graph().is_isomorphic(D1.vertex_facet_graph()))
            tester.assertTrue(P.combinatorial_polyhedron().vertex_facet_graph().is_isomorphic(D2.vertex_facet_graph()))

    @cached_method(do_pickle=True)
    def f_vector(self, num_threads=None, parallelization_depth=None):
        r"""
        Return the f-vector.

        INPUT:

        - ``num_threads`` -- integer (optional); specify the number of threads;
          otherwise determined by :func:`~sage.parallel.ncpus.ncpus`

        - ``parallelization_depth`` -- integer (optional); specify
          how deep in the lattice the parallelization is done

        OUTPUT:

        Return a vector whose `i`-th entry is the number of
        `i-2`-dimensional faces of the polytope.

        .. NOTE::

            The ``vertices`` as given by :meth:`Polyhedron_base.vertices`
            do not need to correspond to `0`-dimensional faces. If a polyhedron
            contains `k` lines they correspond to `k`-dimensional faces.
            See example below

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1, 2, 3], [1, 3, 2],
            ....:     [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1], [0, 0, 0]])
            sage: p.f_vector()
            (1, 7, 12, 7, 1)

            sage: polytopes.cyclic_polytope(4,10).f_vector()
            (1, 10, 45, 70, 35, 1)

            sage: polytopes.hypercube(5).f_vector()
            (1, 32, 80, 80, 40, 10, 1)

        Polyhedra with lines do not have `0`-faces::

            sage: Polyhedron(ieqs=[[1,-1,0,0],[1,1,0,0]]).f_vector()
            (1, 0, 0, 2, 1)

        However, the method :meth:`Polyhedron_base.vertices` returns
        two points that belong to the ``Vrepresentation``::

            sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
            sage: P.vertices()
            (A vertex at (1, 0), A vertex at (-1, 0))
            sage: P.f_vector()
            (1, 0, 2, 1)

        TESTS:

        Check that :trac:`28828` is fixed::

            sage: P.f_vector().is_immutable()
            True

        The cache of the f-vector is being pickled::

            sage: P = polytopes.cube()
            sage: P.f_vector()
            (1, 8, 12, 6, 1)
            sage: Q = loads(dumps(P))
            sage: Q.f_vector.is_in_cache()
            True
        """
        return self.combinatorial_polyhedron().f_vector(num_threads, parallelization_depth)

    def flag_f_vector(self, *args):
        r"""
        Return the flag f-vector.

        For each `-1 < i_0 < \dots < i_n < d` the flag f-vector
        counts the number of flags `F_0 \subset \dots \subset F_n`
        with `F_j` of dimension `i_j` for each `0 \leq j \leq n`,
        where `d` is the dimension of the polyhedron.

        INPUT:

        - ``args`` -- integers (optional); specify an entry of the
          flag-f-vector; must be an increasing sequence of integers

        OUTPUT:

        - a dictionary, if no arguments were given

        - an Integer, if arguments were given

        EXAMPLES:

        Obtain the entire flag-f-vector::

            sage: P = polytopes.twenty_four_cell()
            sage: P.flag_f_vector()
                {(-1,): 1,
                 (0,): 24,
                 (0, 1): 192,
                 (0, 1, 2): 576,
                 (0, 1, 2, 3): 1152,
                 (0, 1, 3): 576,
                 (0, 2): 288,
                 (0, 2, 3): 576,
                 (0, 3): 144,
                 (1,): 96,
                 (1, 2): 288,
                 (1, 2, 3): 576,
                 (1, 3): 288,
                 (2,): 96,
                 (2, 3): 192,
                 (3,): 24,
                 (4,): 1}

        Specify an entry::

            sage: P.flag_f_vector(0,3)
            144
            sage: P.flag_f_vector(2)
            96

        Leading ``-1`` and trailing entry of dimension are allowed::

            sage: P.flag_f_vector(-1,0,3)
            144
            sage: P.flag_f_vector(-1,0,3,4)
            144

        One can get the number of trivial faces::

            sage: P.flag_f_vector(-1)
            1
            sage: P.flag_f_vector(4)
            1

        Polyhedra with lines, have ``0`` entries accordingly::

            sage: P = (Polyhedron(lines=[[1]]) * polytopes.cross_polytope(3))
            sage: P.flag_f_vector()
            {(-1,): 1,
             (0, 1): 0,
             (0, 1, 2): 0,
             (0, 1, 3): 0,
             (0, 2): 0,
             (0, 2, 3): 0,
             (0, 3): 0,
             (0,): 0,
             (1, 2): 24,
             (1, 2, 3): 48,
             (1, 3): 24,
             (1,): 6,
             (2, 3): 24,
             (2,): 12,
             (3,): 8,
             4: 1}

        If the arguments are not stricly increasing or out of range, a key error is raised::

            sage: P.flag_f_vector(-1,0,3,6)
            Traceback (most recent call last):
            ...
            KeyError: (0, 3, 6)
            sage: P.flag_f_vector(-1,3,0)
            Traceback (most recent call last):
            ...
            KeyError: (3, 0)
        """
        flag = self._flag_f_vector()
        if len(args) == 0:
            return flag
        elif len(args) == 1:
            return flag[(args[0],)]
        else:
            dim = self.dimension()
            if args[0] == -1:
                args = args[1:]
            if args[-1] == dim:
                args = args[:-1]
            return flag[tuple(args)]

    @cached_method(do_pickle=True)
    def _flag_f_vector(self):
        r"""
        Return the flag-f-vector.

        See :meth:`flag_f_vector`.

        TESTS::

            sage: polytopes.hypercube(4)._flag_f_vector()
            {(-1,): 1,
            (0,): 16,
            (0, 1): 64,
            (0, 1, 2): 192,
            (0, 1, 2, 3): 384,
            (0, 1, 3): 192,
            (0, 2): 96,
            (0, 2, 3): 192,
            (0, 3): 64,
            (1,): 32,
            (1, 2): 96,
            (1, 2, 3): 192,
            (1, 3): 96,
            (2,): 24,
            (2, 3): 48,
            (3,): 8,
            (4,): 1}
        """
        return self.combinatorial_polyhedron()._flag_f_vector()

    def vertex_graph(self):
        """
        Return a graph in which the vertices correspond to vertices
        of the polyhedron, and edges to edges.

        ..NOTE::

            The graph of a polyhedron with lines has no vertices,
            as the polyhedron has no vertices (`0`-faces).

            The method :meth:`Polyhedron_base:vertices` returns
            the defining points in this case.

        EXAMPLES::

            sage: g3 = polytopes.hypercube(3).vertex_graph(); g3
            Graph on 8 vertices
            sage: g3.automorphism_group().cardinality()
            48
            sage: s4 = polytopes.simplex(4).vertex_graph(); s4
            Graph on 5 vertices
            sage: s4.is_eulerian()
            True

        The graph of an unbounded polyhedron
        is the graph of the bounded complex::

            sage: open_triangle = Polyhedron(vertices=[[1,0], [0,1]],
            ....:                            rays    =[[1,1]])
            sage: open_triangle.vertex_graph()
            Graph on 2 vertices

        The graph of a polyhedron with lines has no vertices::

            sage: line = Polyhedron(lines=[[0,1]])
            sage: line.vertex_graph()
            Graph on 0 vertices

        TESTS:

        Check for a line segment (:trac:`30545`)::

            sage: polytopes.simplex(1).graph().edges()
            [(A vertex at (0, 1), A vertex at (1, 0), None)]
        """
        return self.combinatorial_polyhedron().vertex_graph()

    graph = vertex_graph

    def vertex_digraph(self, f, increasing=True):
        r"""
        Return the directed graph of the polyhedron according to a linear form.

        The underlying undirected graph is the graph of vertices and edges.

        INPUT:

        - ``f`` -- a linear form. The linear form can be provided as:

            - a vector space morphism with one-dimensional codomain, (see
              :meth:`sage.modules.vector_space_morphism.linear_transformation`
              and
              :class:`sage.modules.vector_space_morphism.VectorSpaceMorphism`)
            - a vector ; in this case the linear form is obtained by duality
              using the dot product: ``f(v) = v.dot_product(f)``.

        - ``increasing`` -- boolean (default ``True``) whether to orient
          edges in the increasing or decreasing direction.

        By default, an edge is oriented from `v` to `w` if
        `f(v) \leq f(w)`.

        If `f(v)=f(w)`, then two opposite edges are created.

        EXAMPLES::

            sage: penta = Polyhedron([[0,0],[1,0],[0,1],[1,2],[3,2]])
            sage: G = penta.vertex_digraph(vector([1,1])); G
            Digraph on 5 vertices
            sage: G.sinks()
            [A vertex at (3, 2)]

            sage: A = matrix(ZZ, [[1], [-1]])
            sage: f = linear_transformation(A)
            sage: G = penta.vertex_digraph(f) ; G
            Digraph on 5 vertices
            sage: G.is_directed_acyclic()
            False

        .. SEEALSO::

            :meth:`vertex_graph`
        """
        from sage.modules.vector_space_morphism import VectorSpaceMorphism
        if isinstance(f, VectorSpaceMorphism):
            if f.codomain().dimension() == 1:
                orientation_check = lambda v: f(v) >= 0
            else:
                raise TypeError('the linear map f must have '
                                'one-dimensional codomain')
        else:
            try:
                if f.is_vector():
                    orientation_check = lambda v: v.dot_product(f) >= 0
                else:
                    raise TypeError('f must be a linear map or a vector')
            except AttributeError:
                raise TypeError('f must be a linear map or a vector')
        if not increasing:
            f = -f
        from sage.graphs.digraph import DiGraph
        dg = DiGraph()
        for j in range(self.n_vertices()):
            vj = self.Vrepresentation(j)
            for vi in vj.neighbors():
                if orientation_check(vj.vector() - vi.vector()):
                    dg.add_edge(vi, vj)
        return dg

    def polar(self, in_affine_span=False):
        """
        Return the polar (dual) polytope.

        The original vertices are translated so that their barycenter
        is at the origin, and then the vertices are used as the
        coefficients in the polar inequalities.

        The polytope must be full-dimensional, unless ``in_affine_span`` is ``True``.
        If ``in_affine_span`` is ``True``, then the operation will be performed in the
        linear/affine span of the polyhedron (after translation).

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,1],[0,1,0],[1,0,0],[0,0,0],[1,1,1]], base_ring=QQ)
            sage: p
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 5 vertices
            sage: p.polar()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices

            sage: cube = polytopes.hypercube(3)
            sage: octahedron = polytopes.cross_polytope(3)
            sage: cube_dual = cube.polar()
            sage: octahedron == cube_dual
            True

        ``in_affine_span`` somewhat ignores equations, performing the polar in the
        spanned subspace (after translating barycenter to origin)::

            sage: P = polytopes.simplex(3, base_ring=QQ)
            sage: P.polar(in_affine_span=True)
            A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 4 vertices

        Embedding the polytope in a higher dimension, commutes with polar in this case::

            sage: point = Polyhedron([[0]])
            sage: P = polytopes.cube().change_ring(QQ)
            sage: (P*point).polar(in_affine_span=True) == P.polar()*point
            True

        TESTS::

            Check that :trac:`25081` is fixed::

            sage: C = polytopes.hypercube(4,backend='cdd')
            sage: C.polar().backend()
            'cdd'

        Check that :trac:`28850` is fixed::

            sage: P = polytopes.simplex(3, base_ring=QQ)
            sage: P.polar()
            Traceback (most recent call last):
            ...
            ValueError: not full-dimensional; try with 'in_affine_span=True'

        Check that the double description is set up correctly::

            sage: P = Polyhedron([[1,0],[0,1],[-1,-1]], backend='field')
            sage: Q = P.change_ring(QQ, backend='ppl')
            sage: P.polar() == Q.polar()
            True

            sage: P = polytopes.simplex(4, backend='field')
            sage: Q = P.change_ring(QQ, backend='ppl')
            sage: P.polar(in_affine_span=True) == Q.polar(in_affine_span=True)
            True

        Check that it works, even when the equations are not orthogonal to each other::

            sage: P = polytopes.cube()*Polyhedron([[0,0,0]])
            sage: P = P.change_ring(QQ)

            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: from sage.geometry.polyhedron.parent import Polyhedra_field
            sage: parent = Polyhedra_field(QQ, 6, 'field')
            sage: equations = [[0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, -1, 1, -1], [0, 0, 0, 0, 1, -1, -1]]
            sage: Q = Polyhedron_field(parent, [P.vertices(), [], []], [P.inequalities(), equations],
            ....:                      Vrep_minimal=True, Hrep_minimal=True)
            sage: Q == P
            True
            sage: Q.polar(in_affine_span=True) == P.polar(in_affine_span=True)
            True
        """
        if not self.is_compact():
            raise ValueError("not a polytope")
        if not in_affine_span and not self.dim() == self.ambient_dim():
            raise ValueError("not full-dimensional; try with 'in_affine_span=True'")

        t_Vrep, t_Hrep, parent = self._translation_double_description(-self.center())
        t_verts = t_Vrep[0]
        t_ieqs = t_Hrep[0]
        t_eqns = t_Hrep[1]

        new_ieqs = ((1,) + tuple(-v) for v in t_verts)
        if self.n_vertices() == 1:
            new_verts = self.vertices()
        elif not self.n_equations():
            new_verts = ((-h/h[0])[1:] for h in t_ieqs)
        else:
            # Transform the equations such that the normals are pairwise orthogonal.
            t_eqns = list(t_eqns)
            for i, h in enumerate(t_eqns):
                for h1 in t_eqns[:i]:
                    a = h[1:]*h1[1:]
                    if a:
                        b = h1[1:]*h1[1:]
                        t_eqns[i] = b*h - a*h1

            def move_vertex_to_subspace(vertex):
                for h in t_eqns:
                    offset = vertex*h[1:]+h[0]
                    vertex = vertex-h[1:]*offset/(h[1:]*h[1:])
                return vertex

            new_verts = (move_vertex_to_subspace((-h/h[0])[1:]) for h in t_ieqs)

        pref_rep = 'Hrep' if self.n_vertices() <= self.n_inequalities() else 'Vrep'

        return parent.element_class(parent, [new_verts, [], []],
                                    [new_ieqs, t_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def is_self_dual(self):
        r"""
        Return whether the polytope is self-dual.

        A polytope is self-dual if its face lattice is isomorphic to the face
        lattice of its dual polytope.

        EXAMPLES::

            sage: polytopes.simplex().is_self_dual()
            True
            sage: polytopes.twenty_four_cell().is_self_dual()
            True
            sage: polytopes.cube().is_self_dual()
            False
            sage: polytopes.hypersimplex(5,2).is_self_dual()
            False
            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]]).is_self_dual()
            Traceback (most recent call last):
            ...
            ValueError: polyhedron has to be compact

        """
        if not self.is_compact():
            raise ValueError("polyhedron has to be compact")

        n = self.n_vertices()
        m = self.n_facets()
        if n != m:
            return False

        G1 = self.vertex_facet_graph()
        G2 = G1.reverse()
        return G1.is_isomorphic(G2)

    def pyramid(self):
        """
        Return a polyhedron that is a pyramid over the original.

        EXAMPLES::

            sage: square = polytopes.hypercube(2);  square
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: egyptian_pyramid = square.pyramid();  egyptian_pyramid
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 5 vertices
            sage: egyptian_pyramid.n_vertices()
            5
            sage: for v in egyptian_pyramid.vertex_generator(): print(v)
            A vertex at (0, -1, -1)
            A vertex at (0, -1, 1)
            A vertex at (0, 1, -1)
            A vertex at (0, 1, 1)
            A vertex at (1, 0, 0)

        TESTS::

            sage: polytopes.simplex(backend='cdd').pyramid().backend()
            'cdd'
        """
        assert self.is_compact(), "Not a polytope."
        c = self.center()

        from itertools import chain
        new_verts = chain(([0] + x for x in self.Vrep_generator()),
                          [[1] + list(c)])
        new_ieqs = chain(([i.b()] + [-c*i.A() - i.b()] + list(i.A()) for i in self.inequalities()),
                         [[0, 1] + [0]*self.ambient_dim()])
        new_eqns = ([e.b()] + [0] + list(e.A()) for e in self.equations())

        pref_rep = 'Hrep' if self.n_vertices() > self.n_inequalities() else 'Vrep'
        parent = self.parent().base_extend(self.center().parent(), ambient_dim=self.ambient_dim()+1)

        if self.n_vertices() == 1:
            # Fix the polyhedron with one vertex.
            return parent.element_class(parent, [new_verts, [], []], None)

        return parent.element_class(parent, [new_verts, [], []],
                                    [new_ieqs, new_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _test_pyramid(self, tester=None, **options):
        """
        Run tests on the methods related to pyramids.

        TESTS:

            sage: polytopes.regular_polygon(4)._test_pyramid()
        """
        if tester is None:
            tester = self._tester(**options)

        def check_pyramid_certificate(P, cert):
            others = set(v for v in P.vertices() if not v == cert)
            if len(others):
                tester.assertTrue(any(set(f.ambient_Vrepresentation()) == others for f in P.facets()))

        if self.is_compact():
            b, cert = self.is_pyramid(certificate=True)
            if b:
                check_pyramid_certificate(self, cert)

            if 1 < self.n_vertices() < 50 and self.n_facets() < 50:
                pyr = self.pyramid()
                b, cert = pyr.is_pyramid(certificate=True)
                tester.assertTrue(b)
                check_pyramid_certificate(pyr, cert)
        else:
            with tester.assertRaises(AssertionError):
                pyr = self.pyramid()

        if self.is_compact() and 1 < self.n_vertices() < 50 and self.n_facets() < 50:
            # Check the pyramid of the polar.
            self_fraction_field = self.base_extend(QQ)

            polar = self_fraction_field.polar(in_affine_span=True)
            pyr_polar = polar.pyramid()
            b, cert = pyr_polar.is_pyramid(certificate=True)
            tester.assertTrue(b)
            check_pyramid_certificate(pyr_polar, cert)

            pyr = self_fraction_field.pyramid()
            polar_pyr = pyr.polar(in_affine_span=True)
            b, cert = polar_pyr.is_pyramid(certificate=True)
            tester.assertTrue(b)
            check_pyramid_certificate(polar_pyr, cert)

            tester.assertTrue(pyr_polar.is_combinatorially_isomorphic(pyr_polar))

            # Basic properties of the pyramid.

            # Check that the prism preserves the backend.
            tester.assertEqual(pyr.backend(), self.backend())

            tester.assertEqual(1 + self.n_vertices(), pyr.n_vertices())
            tester.assertEqual(self.n_equations(), pyr.n_equations())
            tester.assertEqual(1 + self.n_inequalities(), pyr.n_inequalities())

            if self.n_vertices() < 15 and self.n_facets() < 15:
                pyr._test_basic_properties()

    def bipyramid(self):
        """
        Return a polyhedron that is a bipyramid over the original.

        EXAMPLES::

            sage: octahedron = polytopes.cross_polytope(3)
            sage: cross_poly_4d = octahedron.bipyramid()
            sage: cross_poly_4d.n_vertices()
            8
            sage: q = [list(v) for v in cross_poly_4d.vertex_generator()]
            sage: q
            [[-1, 0, 0, 0],
             [0, -1, 0, 0],
             [0, 0, -1, 0],
             [0, 0, 0, -1],
             [0, 0, 0, 1],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [1, 0, 0, 0]]

        Now check that bipyramids of cross-polytopes are cross-polytopes::

            sage: q2 = [list(v) for v in polytopes.cross_polytope(4).vertex_generator()]
            sage: [v in q2 for v in q]
            [True, True, True, True, True, True, True, True]

        TESTS::

            sage: polytopes.simplex(backend='cdd').bipyramid().backend()
            'cdd'
        """
        c = self.center()
        from itertools import chain
        new_verts = chain(([0] + list(x) for x in self.vertex_generator()),
                          [[1] + list(c), [-1] + list(c)])
        new_rays =  ([0] + r for r in self.rays())
        new_lines = ([0] + l for l in self.lines())
        new_ieqs = chain(([i.b()] + [ c*i.A() + i.b()] + list(i.A()) for i in self.inequalities()),
                         ([i.b()] + [-c*i.A() - i.b()] + list(i.A()) for i in self.inequalities()))
        new_eqns = ([e.b()] + [0] + list(e.A()) for e in self.equations())

        pref_rep = 'Hrep' if 2 + (self.n_vertices() + self.n_rays()) >= 2*self.n_inequalities() else 'Vrep'
        parent = self.parent().base_extend(self.center().parent(), ambient_dim=self.ambient_dim()+1)

        if c not in self.relative_interior():
            # Fix polyhedra with non-proper center.
            return parent.element_class(parent, [new_verts, new_rays, new_lines], None)

        return parent.element_class(parent, [new_verts, new_rays, new_lines],
                                    [new_ieqs, new_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _test_bipyramid(self, tester=None, **options):
        """
        Run tests on the method :meth:`.bipyramid`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_bipyramid()
        """
        if tester is None:
            tester = self._tester(**options)

        if (self.n_vertices() + self.n_rays() >= 40
                or self.n_facets() >= 40
                or self.n_vertices() <= 1):
            return

        bipyramid = self.bipyramid()

        # Check that the double description is set up correctly.
        if self.base_ring().is_exact() and self.n_vertices() + self.n_rays() < 15 and self.n_facets() < 15:
            bipyramid._test_basic_properties(tester)

        # Check that the bipyramid preserves the backend.
        tester.assertEqual(bipyramid.backend(), self.backend())

        if self.center() not in self.relative_interior():
            # In this case (unbounded) the bipyramid behaves a bit different.
            return

        tester.assertEqual(2 + self.n_vertices(), bipyramid.n_vertices())
        tester.assertEqual(self.n_rays(), bipyramid.n_rays())
        tester.assertEqual(self.n_lines(), bipyramid.n_lines())
        tester.assertEqual(self.n_equations(), bipyramid.n_equations())
        tester.assertEqual(2*self.n_inequalities(), bipyramid.n_inequalities())

        if not self.is_compact():
            # ``is_bipyramid`` is only implemented for compact polyhedra.
            return

        b, cert = bipyramid.is_bipyramid(certificate=True)
        tester.assertTrue(b)

        if not self.is_bipyramid() and self.base_ring().is_exact():
            # In this case the certificate is unique.

            R = self.base_ring()
            a = (R(1),) + tuple(self.center())
            b = (R(-1),) + tuple(self.center())
            c, d = [tuple(v) for v in cert]
            tester.assertEqual(sorted([a, b]), sorted([c, d]))

    def prism(self):
        """
        Return a prism of the original polyhedron.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: cube = square.prism()
            sage: cube
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: hypercube = cube.prism()
            sage: hypercube.n_vertices()
            16

        TESTS::

            sage: polytopes.simplex(backend='cdd').prism().backend()
            'cdd'
        """
        from itertools import chain
        new_verts = chain(([0] + v for v in self.vertices()),
                          ([1] + v for v in self.vertices()))
        new_rays =  ([0] + r for r in self.rays())
        new_lines = ([0] + l for l in self.lines())
        new_eqns = ([e.b()] + [0] + list(e[1:]) for e in self.equations())
        new_ieqs = chain(([i.b()] + [0] + list(i[1:]) for i in self.inequalities()),
                         [[0, 1] + [0]*self.ambient_dim(), [1, -1] + [0]*self.ambient_dim()])

        pref_rep = 'Hrep' if 2*(self.n_vertices() + self.n_rays()) >= self.n_inequalities() + 2 else 'Vrep'
        parent = self.parent().change_ring(self.base_ring(), ambient_dim=self.ambient_dim()+1)

        if not self.vertices():
            # Fix the empty polyhedron.
            return parent.element_class(parent, [[], [], []], None)

        return parent.element_class(parent, [new_verts, new_rays, new_lines],
                                    [new_ieqs, new_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _test_prism(self, tester=None, **options):
        """
        Run tests on the method :meth:`.prism`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_prism()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.n_vertices() + self.n_rays() < 40 and self.n_facets() < 40:
            prism = self.prism()

            # Check that the double description is set up correctly.
            if self.base_ring().is_exact() and self.n_vertices() + self.n_rays() < 15 and self.n_facets() < 15:
                prism._test_basic_properties(tester)

            # Check that the prism preserves the backend.
            tester.assertEqual(prism.backend(), self.backend())

            tester.assertEqual(2*self.n_vertices(), prism.n_vertices())
            tester.assertEqual(self.n_rays(), prism.n_rays())
            tester.assertEqual(self.n_lines(), prism.n_lines())
            tester.assertEqual(self.n_equations(), prism.n_equations())
            if self.is_empty():
                return

            tester.assertEqual(2 + self.n_inequalities(), prism.n_inequalities())

            if not self.is_compact():
                # ``is_prism`` only implemented for compact polyhedra.
                return

            b, cert = prism.is_prism(certificate=True)
            tester.assertTrue(b)

            if not self.is_prism() and self.base_ring().is_exact():
                # In this case the certificate is unique.

                R = self.base_ring()
                cert_set = set(frozenset(tuple(v) for v in f) for f in cert)
                expected_cert = set(frozenset((i,) + tuple(v)
                                              for v in self.vertices())
                                    for i in (R(0), R(1)))
                tester.assertEqual(cert_set, expected_cert)

    def one_point_suspension(self, vertex):
        """
        Return the one-point suspension of ``self`` by splitting the vertex
        ``vertex``.

        The resulting polyhedron has one more vertex and its dimension
        increases by one.

        INPUT:

        - ``vertex`` -- a Vertex of ``self``

        EXAMPLES::

            sage: cube = polytopes.cube()
            sage: v = cube.vertices()[0]
            sage: ops_cube = cube.one_point_suspension(v)
            sage: ops_cube.f_vector()
            (1, 9, 24, 24, 9, 1)

            sage: pentagon  = polytopes.regular_polygon(5)
            sage: v = pentagon.vertices()[0]
            sage: ops_pentagon = pentagon.one_point_suspension(v)
            sage: ops_pentagon.f_vector()
            (1, 6, 12, 8, 1)

        It works with a polyhedral face as well::

            sage: vv = cube.faces(0)[1]
            sage: ops_cube2 = cube.one_point_suspension(vv)
            sage: ops_cube == ops_cube2
            True

        .. SEEALSO::

            :meth:`face_split`

        TESTS::

            sage: e = cube.faces(1)[0]
            sage: cube.one_point_suspension(e)
            Traceback (most recent call last):
            ...
            TypeError: the vertex A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices should be a Vertex or PolyhedronFace of dimension 0
        """
        from sage.geometry.polyhedron.representation import Vertex
        from sage.geometry.polyhedron.face import PolyhedronFace
        if isinstance(vertex, Vertex):
            return self.face_split(vertex)
        elif isinstance(vertex, PolyhedronFace) and vertex.dim() == 0:
            return self.face_split(vertex)
        else:
            raise TypeError("the vertex {} should be a Vertex or PolyhedronFace of dimension 0".format(vertex))

    def face_split(self, face):
        """
        Return the face splitting of the face ``face``.

        Splitting a face correspond to the bipyramid (see :meth:`bipyramid`)
        of ``self`` where the two new vertices are placed above and below
        the center of ``face`` instead of the center of the whole polyhedron.
        The two new vertices are placed in the new dimension at height `-1` and
        `1`.

        INPUT:

        - ``face`` -- a PolyhedronFace or a Vertex

        EXAMPLES::

            sage: pentagon  = polytopes.regular_polygon(5)
            sage: f = pentagon.faces(1)[0]
            sage: fsplit_pentagon = pentagon.face_split(f)
            sage: fsplit_pentagon.f_vector()
            (1, 7, 14, 9, 1)

        TESTS:

        Check that :trac:`28668` is fixed::

            sage: P = polytopes.octahedron()
            sage: P.face_split(P.faces(2)[0])
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 8 vertices

        .. SEEALSO::

            :meth:`one_point_suspension`
        """
        from sage.geometry.polyhedron.representation import Vertex
        from sage.geometry.polyhedron.face import PolyhedronFace
        if isinstance(face, Vertex):
            new_vertices = [list(x) + [0] for x in self.vertex_generator()] + \
                           [list(face) + [x] for x in [-1, 1]]  # Splitting the vertex
        elif isinstance(face, PolyhedronFace):
            new_vertices = [list(x) + [0] for x in self.vertex_generator()] + \
                           [list(face.as_polyhedron().center()) + [x] for x in [-1, 1]]  # Splitting the face
        else:
            raise TypeError("the face {} should be a Vertex or PolyhedronFace".format(face))

        new_rays = []
        new_rays.extend( [ r + [0] for r in self.ray_generator() ] )

        new_lines = []
        new_lines.extend( [ l + [0] for l in self.line_generator() ] )

        parent = self.parent().change_ring(self.base_ring().fraction_field(), ambient_dim=self.ambient_dim()+1)
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    def projection(self, projection=None):
        """
        Return a projection object.

        INPUT:

        - ``proj`` -- a projection function

        OUTPUT:

        The identity projection. This is useful for plotting
        polyhedra.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.schlegel_projection` for a more interesting projection.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: proj = p.projection()
            sage: proj
            The projection of a polyhedron into 3 dimensions
        """
        from .plot import Projection
        if projection is not None:
            self.projection = Projection(self, projection)
        else:
            self.projection = Projection(self)
        return self.projection

    def render_solid(self, **kwds):
        """
        Return a solid rendering of a 2- or 3-d polytope.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p_solid = p.render_solid(opacity = .7)
            sage: type(p_solid)
            <type 'sage.plot.plot3d.index_face_set.IndexFaceSet'>
        """
        proj = self.projection()
        if self.ambient_dim() == 3:
            return proj.render_solid_3d(**kwds)
        if self.ambient_dim() == 2:
            return proj.render_fill_2d(**kwds)
        raise ValueError("render_solid is only defined for 2 and 3 dimensional polyhedra")

    def render_wireframe(self, **kwds):
        """
        For polytopes in 2 or 3 dimensions, return the edges
        as a list of lines.

        EXAMPLES::

            sage: p = Polyhedron([[1,2,],[1,1],[0,0]])
            sage: p_wireframe = p.render_wireframe()
            sage: p_wireframe._objects
            [Line defined by 2 points, Line defined by 2 points, Line defined by 2 points]
        """
        proj = self.projection()
        if self.ambient_dim() == 3:
            return proj.render_wireframe_3d(**kwds)
        if self.ambient_dim() == 2:
            return proj.render_outline_2d(**kwds)
        raise ValueError("render_wireframe is only defined for 2 and 3 dimensional polyhedra")

    def schlegel_projection(self, facet=None, position=None):
        """
        Return the Schlegel projection.

        * The facet is orthonormally transformed into its affine hull.

        * The position specifies a point coming out of the barycenter of the
          facet from which the other vertices will be projected into the facet.

        INPUT:

        - ``facet`` -- a PolyhedronFace. The facet into which the Schlegel
          diagram is created. The default is the first facet.

        - ``position`` -- a positive number. Determines a relative distance
          from the barycenter of ``facet``. A value close to 0 will place the
          projection point close to the facet and a large value further away.
          Default is `1`. If the given value is too large, an error is returned.

        OUTPUT:

        A :class:`~sage.geometry.polyhedron.plot.Projection` object.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: sch_proj = p.schlegel_projection()
            sage: schlegel_edge_indices = sch_proj.lines
            sage: schlegel_edges = [sch_proj.coordinates_of(x) for x in schlegel_edge_indices]
            sage: len([x for x in schlegel_edges if x[0][0] > 0])
            8

        The Schlegel projection preserves the convexity of facets, see :trac:`30015`::

            sage: fcube = polytopes.hypercube(4)
            sage: tfcube = fcube.face_truncation(fcube.faces(0)[0])
            sage: tfcube.facets()[-1]
            A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 8 vertices
            sage: sp = tfcube.schlegel_projection(tfcube.facets()[-1])
            sage: sp.plot()
            Graphics3d Object

        The same truncated cube but see inside the tetrahedral facet::

            sage: tfcube.facets()[4]
            A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 4 vertices
            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4])
            sage: sp.plot()
            Graphics3d Object

        A different values of ``position`` changes the projection::

            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4],1/2)
            sage: sp.plot()
            Graphics3d Object
            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4],4)
            sage: sp.plot()
            Graphics3d Object

        A value which is too large give a projection point that sees more than
        one facet resulting in a error::

            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4],5)
            Traceback (most recent call last):
            ...
            ValueError: the chosen position is too large
        """
        proj = self.projection()
        return proj.schlegel(facet, position)

    def _volume_lrs(self, verbose=False):
        """
        Computes the volume of a polytope using lrs.

        OUTPUT:

        The volume, cast to RDF (although lrs seems to output a
        rational value this must be an approximation in some cases).

        EXAMPLES::

            sage: polytopes.hypercube(3)._volume_lrs() # optional - lrslib
            8.0
            sage: (polytopes.hypercube(3)*2)._volume_lrs() # optional - lrslib
            64.0
            sage: polytopes.twenty_four_cell()._volume_lrs() # optional - lrslib
            2.0

        REFERENCES:

        - David Avis's lrs program.
        """
        from sage.features.lrs import Lrs
        Lrs().require()

        from sage.misc.temporary_file import tmp_filename
        from subprocess import Popen, PIPE

        in_str = self.cdd_Vrepresentation()
        in_str += 'volume'
        in_filename = tmp_filename()
        in_file = open(in_filename, 'w')
        in_file.write(in_str)
        in_file.close()
        if verbose:
            print(in_str)

        lrs_procs = Popen(['lrs', in_filename],
                          stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = lrs_procs.communicate()
        ans = bytes_to_str(ans)
        err = bytes_to_str(err)
        if verbose:
            print(ans)
        # FIXME: check err

        for a_line in ans.splitlines():
            if 'Volume=' in a_line:
                volume = a_line.split('Volume=')[1]
                volume = RDF(QQ(volume))
                return volume

        raise ValueError("lrs did not return a volume")

    def _volume_latte(self, verbose=False, algorithm='triangulate', **kwargs):
        """
        Computes the volume of a polytope using LattE integrale.

        INPUT:

        - ``arg`` -- a cdd or LattE description string

        - ``algorithm`` -- (default: 'triangulate') the integration method. Use 'triangulate' for
          polytope triangulation or 'cone-decompose' for tangent cone decomposition method.

        - ``raw_output`` -- if ``True`` then return directly the output string from LattE.

        - ``verbose`` -- if ``True`` then return directly verbose output from LattE.

        - For all other options, consult the LattE manual.

        OUTPUT:

        A rational value, or a string if ``raw_output`` if set to ``True``.

        .. NOTE::

            This function depends on LattE (i.e., the ``latte_int`` optional
            package). See the LattE documentation for further details.

        EXAMPLES::

            sage: polytopes.hypercube(3)._volume_latte() # optional - latte_int
            8
            sage: (polytopes.hypercube(3)*2)._volume_latte() # optional - latte_int
            64
            sage: polytopes.twenty_four_cell()._volume_latte() # optional - latte_int
            2
            sage: polytopes.cuboctahedron()._volume_latte() # optional - latte_int
            20/3

        TESTS:

        Testing triangulate algorithm::

            sage: polytopes.cuboctahedron()._volume_latte(algorithm='triangulate') # optional - latte_int
            20/3

        Testing cone decomposition algorithm::

            sage: polytopes.cuboctahedron()._volume_latte(algorithm='cone-decompose') # optional - latte_int
            20/3

        Testing raw output::

            sage: polytopes.cuboctahedron()._volume_latte(raw_output=True) # optional - latte_int
            '20/3'

        Testing inexact rings::

            sage: P = Polyhedron(vertices=[[0,0],[1,0],[0,1]],base_ring=RDF)
            sage: P.volume(engine='latte')
            Traceback (most recent call last):
            ...
            ValueError: LattE integrale cannot be applied over inexact rings
        """
        from sage.interfaces.latte import integrate
        if self.base_ring() == RDF:
            raise ValueError("LattE integrale cannot be applied over inexact rings")
        else:
            return integrate(self.cdd_Hrepresentation(), algorithm=algorithm, cdd=True, verbose=verbose, **kwargs)

    def _volume_normaliz(self, measure='induced'):
        r"""
        Computes the volume of a polytope using normaliz.

        INPUT:

        - ``measure`` -- (default: 'induced') the measure to take. 'induced'
          correspond to ``EuclideanVolume`` in normaliz and 'induced_lattice'
          correspond to ``Volume`` in normaliz

        OUTPUT:

        A float value (when ``measure`` is 'induced') or a rational number
        (when ``measure`` is 'induced_lattice')

        .. NOTE::

            This function depends on Normaliz (i.e., the ``pynormaliz`` optional
            package). See the Normaliz documentation for further details.

        TESTS::

            sage: P = Polyhedron(vertices=[[0,0],[1,0],[0,1],[1,1]])
            sage: P._volume_normaliz()
            Traceback (most recent call last):
            ...
            TypeError: the backend should be normaliz
        """
        raise TypeError("the backend should be normaliz")

    @cached_method(do_pickle=True)
    def volume(self, measure='ambient', engine='auto', **kwds):
        """
        Return the volume of the polytope.

        INPUT:

        - ``measure`` -- string. The measure to use. Allowed values are:

          * ``ambient`` (default): Lebesgue measure of ambient space (volume)
          * ``induced``: Lebesgue measure of the affine hull (relative volume)
          * ``induced_rational``: Scaling of the Lebesgue measure for rational
            polytopes, such that the unit hypercube has volume 1
          * ``induced_lattice``: Scaling of the Lebesgue measure, such that the
            volume of the hypercube is factorial(n)

        - ``engine`` -- string. The backend to use. Allowed values are:

          * ``'auto'`` (default): choose engine according to measure
          * ``'internal'``: see :meth:`triangulate`
          * ``'TOPCOM'``: see :meth:`triangulate`
          * ``'lrs'``: use David Avis's lrs program (optional)
          * ``'latte'``: use LattE integrale program (optional)
          * ``'normaliz'``: use Normaliz program (optional)

        - ``**kwds`` -- keyword arguments that are passed to the
          triangulation engine

        OUTPUT:

        The volume of the polytope

        EXAMPLES::

            sage: polytopes.hypercube(3).volume()
            8
            sage: (polytopes.hypercube(3)*2).volume()
            64
            sage: polytopes.twenty_four_cell().volume()
            2

        Volume of the same polytopes, using the optional package lrslib
        (which requires a rational polytope).  For mysterious historical
        reasons, Sage casts lrs's exact answer to a float::

            sage: I3 = polytopes.hypercube(3)
            sage: I3.volume(engine='lrs') # optional - lrslib
            8.0
            sage: C24 = polytopes.twenty_four_cell()
            sage: C24.volume(engine='lrs') # optional - lrslib
            2.0

        If the base ring is exact, the answer is exact::

            sage: P5 = polytopes.regular_polygon(5)
            sage: P5.volume()
            2.377641290737884?

            sage: polytopes.icosahedron().volume()
            5/12*sqrt5 + 5/4
            sage: numerical_approx(_) # abs tol 1e9
            2.18169499062491

        When considering lower-dimensional polytopes, we can ask for the
        ambient (full-dimensional), the induced measure (of the affine
        hull) or, in the case of lattice polytopes, for the induced rational measure.
        This is controlled by the parameter `measure`. Different engines
        may have different ideas on the definition of volume of a
        lower-dimensional object::

            sage: P = Polyhedron([[0, 0], [1, 1]])
            sage: P.volume()
            0
            sage: P.volume(measure='induced')
            1.414213562373095?
            sage: P.volume(measure='induced_rational') # optional -- latte_int
            1

            sage: S = polytopes.regular_polygon(6); S
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 6 vertices
            sage: edge = S.faces(1)[4].as_polyhedron()
            sage: edge.vertices()
            (A vertex at (0.866025403784439?, 1/2), A vertex at (0, 1))
            sage: edge.volume()
            0
            sage: edge.volume(measure='induced')
            1

            sage: P = Polyhedron(backend='normaliz',vertices=[[1,0,0],[0,0,1],[-1,1,1],[-1,2,0]]) # optional - pynormaliz
            sage: P.volume()  # optional - pynormaliz
            0
            sage: P.volume(measure='induced')  # optional - pynormaliz
            2.598076211353316?
            sage: P.volume(measure='induced',engine='normaliz')  # optional - pynormaliz
            2.598076211353316
            sage: P.volume(measure='induced_rational')  # optional - pynormaliz, latte_int
            3/2
            sage: P.volume(measure='induced_rational',engine='normaliz')  # optional - pynormaliz
            3/2
            sage: P.volume(measure='induced_lattice')  # optional - pynormaliz
            3

        The same polytope without normaliz backend::

            sage: P = Polyhedron(vertices=[[1,0,0],[0,0,1],[-1,1,1],[-1,2,0]])
            sage: P.volume(measure='induced_lattice',engine='latte')  # optional - latte_int
            3

            sage: Dexact = polytopes.dodecahedron()
            sage: v = Dexact.faces(2)[0].as_polyhedron().volume(measure='induced', engine='internal'); v
            1.53406271079097?
            sage: v = Dexact.faces(2)[4].as_polyhedron().volume(measure='induced', engine='internal'); v
            1.53406271079097?
            sage: RDF(v)    # abs tol 1e-9
            1.53406271079044

            sage: Dinexact = polytopes.dodecahedron(exact=False)
            sage: w = Dinexact.faces(2)[2].as_polyhedron().volume(measure='induced', engine='internal'); RDF(w) # abs tol 1e-9
            1.5340627082974878

            sage: [polytopes.simplex(d).volume(measure='induced') for d in range(1,5)] == [sqrt(d+1)/factorial(d) for d in range(1,5)]
            True

            sage: I = Polyhedron([[-3, 0], [0, 9]])
            sage: I.volume(measure='induced')
            9.48683298050514?
            sage: I.volume(measure='induced_rational') # optional -- latte_int
            3

            sage: T = Polyhedron([[3, 0, 0], [0, 4, 0], [0, 0, 5]])
            sage: T.volume(measure='induced')
            13.86542462386205?
            sage: T.volume(measure='induced_rational') # optional -- latte_int
            1/2

            sage: Q = Polyhedron(vertices=[(0, 0, 1, 1), (0, 1, 1, 0), (1, 1, 0, 0)])
            sage: Q.volume(measure='induced')
            1
            sage: Q.volume(measure='induced_rational') # optional -- latte_int
            1/2

        The volume of a full-dimensional unbounded polyhedron is infinity::

            sage: P = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]])
            sage: P.volume()
            +Infinity

        The volume of a non full-dimensional unbounded polyhedron depends on the measure used::

            sage: P = Polyhedron(ieqs = [[1,1,1],[-1,-1,-1],[3,1,0]]); P
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: P.volume()
            0
            sage: P.volume(measure='induced')
            +Infinity
            sage: P.volume(measure='ambient')
            0
            sage: P.volume(measure='induced_rational')  # optional - pynormaliz
            +Infinity
            sage: P.volume(measure='induced_rational',engine='latte')  # optional - latte_int
            +Infinity

        The volume in `0`-dimensional space is taken by counting measure::

            sage: P = Polyhedron(vertices=[[]]); P
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: P.volume()
            1
            sage: P = Polyhedron(vertices=[]); P
            The empty polyhedron in ZZ^0
            sage: P.volume()
            0

        TESTS:

        The cache of the volume is being pickled::

            sage: P = polytopes.cube()
            sage: P.volume()
            8
            sage: Q = loads(dumps(P))
            sage: Q.volume.is_in_cache()
            True
        """
        from sage.features import FeatureNotPresentError, PythonModule
        if measure == 'induced_rational' and engine not in ['auto', 'latte', 'normaliz']:
            raise RuntimeError("the induced rational measure can only be computed with the engine set to `auto`, `latte`, or `normaliz`")
        if measure == 'induced_lattice' and engine not in ['auto', 'latte', 'normaliz']:
            raise RuntimeError("the induced lattice measure can only be computed with the engine set to `auto`, `latte`, or `normaliz`")
        if engine == 'auto' and measure == 'induced_rational':
            # Enforce a default choice, change if a better engine is found.
            from sage.features.latte import Latte
            try:
                Latte().require()
                engine = 'latte'
            except FeatureNotPresentError:
                try:
                    PythonModule("PyNormaliz", spkg="pynormaliz").require()
                    engine = 'normaliz'
                except FeatureNotPresentError:
                    raise RuntimeError("the induced rational measure can only be computed with the optional packages `latte_int`, or `pynormaliz`")

        if engine == 'auto' and measure == 'induced_lattice':
            # Enforce a default choice, change if a better engine is found.
            try:
                PythonModule("PyNormaliz", spkg="pynormaliz").require()
                engine = 'normaliz'
            except FeatureNotPresentError:
                try:
                    from sage.features.latte import Latte
                    Latte().require()
                    engine = 'latte'
                except FeatureNotPresentError:
                    raise RuntimeError("the induced rational measure can only be computed with the optional packages `latte_int`, or `pynormaliz`")

        if engine == 'auto' and measure == 'ambient' and self.backend() == 'normaliz':
            engine = 'normaliz'

        if measure == 'ambient':
            if self.dim() < self.ambient_dim():
                return self.base_ring().zero()
            elif self.dim() == 0:
                return 1
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'lrs':
                return self._volume_lrs(**kwds)
            elif engine == 'latte':
                return self._volume_latte(**kwds)
            elif engine == 'normaliz':
                return self._volume_normaliz(measure='ambient')

            triangulation = self.triangulate(engine=engine, **kwds)
            pc = triangulation.point_configuration()
            return sum([pc.volume(simplex) for simplex in triangulation]) / ZZ(self.dim()).factorial()
        elif measure == 'induced':
            # if polyhedron is actually full-dimensional, return volume with ambient measure
            if self.dim() == self.ambient_dim():
                return self.volume(measure='ambient', engine=engine, **kwds)
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'normaliz':
                return self._volume_normaliz(measure='euclidean')
            # use an orthogonal transformation, which preserves volume up to a factor provided by the transformation matrix
            affine_hull_data = self.affine_hull_projection(orthogonal=True, as_polyhedron=True, as_affine_map=True)
            A = affine_hull_data.projection_linear_map.matrix()
            Adet = (A.transpose() * A).det()
            scaled_volume = affine_hull_data.image.volume(measure='ambient', engine=engine, **kwds)
            if Adet.is_square():
                sqrt_Adet = Adet.sqrt()
            else:
                sqrt_Adet = AA(Adet).sqrt()
                scaled_volume = AA(scaled_volume)
            return scaled_volume / sqrt_Adet
        elif measure == 'induced_rational':
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'latte':
                return self._volume_latte(**kwds)
            else:  # engine is 'normaliz'
                return self._volume_normaliz(measure='induced_lattice') / ZZ(self.dim()).factorial()
        elif measure == 'induced_lattice':
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'latte':
                return self._volume_latte(**kwds) * ZZ(self.dim()).factorial()
            else:  # engine is 'normaliz'
                return self._volume_normaliz(measure='induced_lattice')
        else:
            raise TypeError("the measure should be `ambient`, `induced`, `induced_rational`, or `induced_lattice`")

    def integrate(self, function, measure='ambient', **kwds):
        r"""
        Return the integral of ``function`` over this polytope.

        INPUT:

        - ``self`` -- Polyhedron

        - ``function`` -- a multivariate polynomial or
          a valid LattE description string for polynomials

        - ``measure`` -- string, the measure to use

          Allowed values are:

          * ``ambient`` (default): Lebesgue measure of ambient space,
          * ``induced``: Lebesgue measure of the affine hull,
          * ``induced_nonnormalized``: Lebesgue measure of the affine hull
            without the normalization by `\sqrt{\det(A^\top A)}` (with
            `A` being the affine transformation matrix; see :meth:`affine_hull`).

        - ``**kwds`` -- additional keyword arguments that
          are passed to the engine

        OUTPUT:

        The integral of the polynomial over the polytope

        .. NOTE::

            The polytope triangulation algorithm is used. This function depends
            on LattE (i.e., the ``latte_int`` optional package).

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: x, y, z = polygens(QQ, 'x, y, z')
            sage: P.integrate(x^2*y^2*z^2)    # optional - latte_int
            8/27

        If the polyhedron has floating point coordinates, an inexact result can
        be obtained if we transform to rational coordinates::

            sage: P = 1.4142*polytopes.cube()
            sage: P_QQ = Polyhedron(vertices=[[QQ(vi) for vi in v] for v in P.vertex_generator()])
            sage: RDF(P_QQ.integrate(x^2*y^2*z^2))                  # optional - latte_int
            6.703841212195228

        Integral over a non full-dimensional polytope::

            sage: x, y = polygens(QQ, 'x, y')
            sage: P = Polyhedron(vertices=[[0,0],[1,1]])
            sage: P.integrate(x*y)    # optional - latte_int
            0
            sage: ixy = P.integrate(x*y, measure='induced'); ixy    # optional - latte_int
            0.4714045207910317?
            sage: ixy.parent()                                      # optional - latte_int
            Algebraic Real Field

        Convert to a symbolic expression::

            sage: ixy.radical_expression()                          # optional - latte_int
            1/3*sqrt(2)

        Another non full-dimensional polytope integration::

            sage: R.<x, y, z> = QQ[]
            sage: P = polytopes.simplex(2)
            sage: V = AA(P.volume(measure='induced')); V.radical_expression()
            1/2*sqrt(3)
            sage: P.integrate(R(1), measure='induced') == V                      # optional - latte_int
            True

        Computing the mass center::

            sage: (P.integrate(x, measure='induced') / V).radical_expression()   # optional - latte_int
            1/3
            sage: (P.integrate(y, measure='induced') / V).radical_expression()   # optional - latte_int
            1/3
            sage: (P.integrate(z, measure='induced') / V).radical_expression()   # optional - latte_int
            1/3

        TESTS:

        Testing a three-dimensional integral::

            sage: P = polytopes.octahedron()
            sage: x, y, z = polygens(QQ, 'x, y, z')
            sage: P.integrate(2*x^2*y^4*z^6+z^2)    # optional - latte_int
            630632/4729725

        Testing a polytope with non-rational vertices::

            sage: P = polytopes.icosahedron()
            sage: P.integrate(x^2*y^2*z^2)    # optional - latte_int
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be ZZ, QQ, or RDF

        Testing a univariate polynomial::

            sage: P = Polyhedron(vertices=[[0],[1]])
            sage: x = polygen(QQ, 'x')
            sage: P.integrate(x)    # optional - latte_int
            1/2

        Testing a polytope with floating point coordinates::

            sage: P = Polyhedron(vertices = [[0, 0], [1, 0], [1.1, 1.1], [0, 1]])
            sage: P.integrate('[[1,[2,2]]]')    # optional - latte_int
            Traceback (most recent call last):
            ...
            TypeError: LattE integrale cannot be applied over inexact rings

        Integration of zero-polynomial::

            sage: R.<x, y, z> = QQ[]
            sage: P = polytopes.simplex(2)
            sage: P.integrate(R(0))
            0
            sage: P.integrate('[]')  # with LattE description string
            0

        ::

            sage: R.<x, y, z> = QQ[]
            sage: P = Polyhedron(vertices=[(0, 0, 1), (0, 1, 0)])
            sage: P.integrate(x^2)
            0
        """
        if function == 0 or function == '[]':
            return self.base_ring().zero()

        if not self.is_compact():
            raise NotImplementedError(
                'integration over non-compact polyhedra not allowed')

        if measure == 'ambient':
            if not self.is_full_dimensional():
                return self.base_ring().zero()

            return self._integrate_latte_(function, **kwds)

        elif measure == 'induced' or measure == 'induced_nonnormalized':
            # if polyhedron is actually full-dimensional,
            # return with ambient measure
            if self.is_full_dimensional():
                return self.integrate(function, measure='ambient', **kwds)

            if isinstance(function, str):
                raise NotImplementedError(
                    'LattE description strings for polynomials not allowed '
                    'when using measure="induced"')

            # use an orthogonal transformation
            affine_hull_data = self.affine_hull_projection(orthogonal=True, return_all_data=True)
            polyhedron = affine_hull_data.image
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(affine_hull_data.section_linear_map.base_ring(), 'x', self.dim())
            coordinate_images = affine_hull_data.section_linear_map.matrix().transpose() * vector(R.gens()) + affine_hull_data.section_translation

            hom = function.parent().hom(coordinate_images)
            function_in_affine_hull = hom(function)

            I = polyhedron.integrate(function_in_affine_hull,
                                     measure='ambient', **kwds)
            if measure == 'induced_nonnormalized':
                return I
            else:
                A = affine_hull_data.projection_linear_map.matrix()
                Adet = (A.transpose() * A).det()
                try:
                    Adet = AA.coerce(Adet)
                except TypeError:
                    pass
                return I / sqrt(Adet)

        else:
            raise ValueError('unknown measure "{}"'.format(measure))

    def _integrate_latte_(self, polynomial, **kwds):
        r"""
        Return the integral of a polynomial over this polytope by calling LattE.

        INPUT:

        - ``polynomial`` -- a multivariate polynomial or
          a valid LattE description string for polynomials

        - ``**kwds`` -- additional keyword arguments that are passed
          to the engine

        OUTPUT:

        The integral of the polynomial over the polytope.

        .. NOTE::

            The polytope triangulation algorithm is used. This function depends
            on LattE (i.e., the ``latte_int`` optional package).

        TESTS::

            sage: P = polytopes.cube()
            sage: x, y, z = polygens(QQ, 'x, y, z')
            sage: P._integrate_latte_(x^2 + y^2*z^2)    # optional - latte_int
            32/9

        ::

            sage: R = PolynomialRing(QQ, '', 0)
            sage: Polyhedron(vertices=[()]).integrate(R(42))
            42
        """
        from sage.interfaces.latte import integrate

        if self.base_ring() == RDF:
            raise TypeError("LattE integrale cannot be applied over inexact rings")
        if self.dimension() == 0:
            vertices = self.vertices()
            assert len(self.vertices()) == 1
            vertex = tuple(vertices[0])
            return polynomial(vertex)
        return integrate(self.cdd_Hrepresentation(),
                         polynomial,
                         cdd=True, **kwds)

    def contains(self, point):
        """
        Test whether the polyhedron contains the given ``point``.

        .. SEEALSO::

            :meth:`interior_contains`, :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point (an iterable)

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,1],[1,-1],[0,0]])
            sage: P.contains( [1,0] )
            True
            sage: P.contains( P.center() )  # true for any convex set
            True

        As a shorthand, one may use the usual ``in`` operator::

            sage: P.center() in P
            True
            sage: [-1,-1] in P
            False

        The point need not have coordinates in the same field as the
        polyhedron::

            sage: ray = Polyhedron(vertices=[(0,0)], rays=[(1,0)], base_ring=QQ)
            sage: ray.contains([sqrt(2)/3,0])        # irrational coordinates are ok
            True
            sage: a = var('a')
            sage: ray.contains([a,0])                # a might be negative!
            False
            sage: assume(a>0)
            sage: ray.contains([a,0])
            True
            sage: ray.contains(['hello', 'kitty'])   # no common ring for coordinates
            False

        The empty polyhedron needs extra care, see :trac:`10238`::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.contains([])
            False
            sage: empty.contains([0])               # not a point in QQ^0
            False
            sage: full = Polyhedron(vertices=[()]); full
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: full.contains([])
            True
            sage: full.contains([0])
            False

        TESTS:

        Passing non-iterable objects does not cause an exception, see :trac:`32013`::

            sage: None in Polyhedron(vertices=[(0,0)], rays=[(1,0)], base_ring=QQ)
            False
        """
        try:
            p = vector(point)
        except TypeError:  # point not iterable or no common ring for elements
            try:
                l = len(point)
            except TypeError:
                return False
            if l > 0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p) != self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.contains(p):
                return False
        return True

    __contains__ = contains

    @cached_method
    def interior(self):
        """
        The interior of ``self``.

        OUTPUT:

        - either an empty polyhedron or an instance of
          :class:`~sage.geometry.relative_interior.RelativeInterior`

        EXAMPLES:

        If the polyhedron is full-dimensional, the result is the
        same as that of :meth:`relative_interior`::

            sage: P_full = Polyhedron(vertices=[[0,0],[1,1],[1,-1]])
            sage: P_full.interior()
            Relative interior of
             a 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices

        If the polyhedron is of strictly smaller dimension than the
        ambient space, its interior is empty::

            sage: P_lower = Polyhedron(vertices=[[0,1], [0,-1]])
            sage: P_lower.interior()
            The empty polyhedron in ZZ^2

        TESTS::

            sage: Empty = Polyhedron(ambient_dim=2); Empty
            The empty polyhedron in ZZ^2
            sage: Empty.interior() is Empty
            True
        """
        if self.is_open():
            return self
        if not self.is_full_dimensional():
            return self.parent().element_class(self.parent(), None, None)
        return self.relative_interior()

    def interior_contains(self, point):
        """
        Test whether the interior of the polyhedron contains the
        given ``point``.

        .. SEEALSO::

            :meth:`contains`, :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0],[1,1],[1,-1]])
            sage: P.contains( [1,0] )
            True
            sage: P.interior_contains( [1,0] )
            False

        If the polyhedron is of strictly smaller dimension than the
        ambient space, its interior is empty::

            sage: P = Polyhedron(vertices=[[0,1],[0,-1]])
            sage: P.contains( [0,0] )
            True
            sage: P.interior_contains( [0,0] )
            False

        The empty polyhedron needs extra care, see :trac:`10238`::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError:  # point not iterable or no common ring for elements
            try:
                l = len(point)
            except TypeError:
                return False
            if l > 0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p) != self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.interior_contains(p):
                return False
        return True

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is relatively open.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)]); P
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: P.is_relatively_open()
            False

            sage: P0 = Polyhedron(vertices=[[1, 2]]); P0
            A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: P0.is_relatively_open()
            True

            sage: Empty = Polyhedron(ambient_dim=2); Empty
            The empty polyhedron in ZZ^2
            sage: Empty.is_relatively_open()
            True

            sage: Line = Polyhedron(vertices=[(1, 1)], lines=[(1, 0)]); Line
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 line
            sage: Line.is_relatively_open()
            True

        """
        return not self.inequalities()

    @cached_method
    def relative_interior(self):
        """
        Return the relative interior of ``self``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)])
            sage: ri_P = P.relative_interior(); ri_P
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: (0, 0) in ri_P
            True
            sage: (1, 0) in ri_P
            False

            sage: P0 = Polyhedron(vertices=[[1, 2]])
            sage: P0.relative_interior() is P0
            True

            sage: Empty = Polyhedron(ambient_dim=2)
            sage: Empty.relative_interior() is Empty
            True

            sage: Line = Polyhedron(vertices=[(1, 1)], lines=[(1, 0)])
            sage: Line.relative_interior() is Line
            True
        """
        if self.is_relatively_open():
            return self
        return RelativeInterior(self)

    def relative_interior_contains(self, point):
        """
        Test whether the relative interior of the polyhedron
        contains the given ``point``.

        .. SEEALSO::

            :meth:`contains`, :meth:`interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point

        OUTPUT:

        ``True`` or ``False``

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)])
            sage: P.contains( (0,0) )
            True
            sage: P.interior_contains( (0,0) )
            False
            sage: P.relative_interior_contains( (0,0) )
            True
            sage: P.relative_interior_contains( (1,0) )
            False

        The empty polyhedron needs extra care, see :trac:`10238`::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.relative_interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError:  # point not iterable or no common ring for elements
            try:
                l = len(point)
            except TypeError:
                return False
            if l > 0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p) != self.ambient_dim():
            return False

        for eq in self.equation_generator():
            if not eq.contains(p):
                return False

        for ine in self.inequality_generator():
            if not ine.interior_contains(p):
                return False

        return True

    def is_simplex(self):
        r"""
        Return whether the polyhedron is a simplex.

        A simplex is a bounded polyhedron with `d+1` vertices, where
        `d` is the dimension.

        EXAMPLES::

            sage: Polyhedron([(0,0,0), (1,0,0), (0,1,0)]).is_simplex()
            True
            sage: polytopes.simplex(3).is_simplex()
            True
            sage: polytopes.hypercube(3).is_simplex()
            False
        """
        return self.is_compact() and (self.dim()+1 == self.n_vertices())

    def neighborliness(self):
        r"""
        Return the largest ``k``, such that the polyhedron is ``k``-neighborly.

        A polyhedron is `k`-neighborly if every set of `n` vertices forms a face
        for `n` up to `k`.

        In case of the `d`-dimensional simplex, it returns `d + 1`.

        .. SEEALSO::

            :meth:`is_neighborly`

        EXAMPLES::

            sage: cube = polytopes.cube()
            sage: cube.neighborliness()
            1
            sage: P = Polyhedron(); P
            The empty polyhedron in ZZ^0
            sage: P.neighborliness()
            0
            sage: P = Polyhedron([[0]]); P
            A 0-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex
            sage: P.neighborliness()
            1
            sage: S = polytopes.simplex(5); S
            A 5-dimensional polyhedron in ZZ^6 defined as the convex hull of 6 vertices
            sage: S.neighborliness()
            6
            sage: C = polytopes.cyclic_polytope(7,10); C
            A 7-dimensional polyhedron in QQ^7 defined as the convex hull of 10 vertices
            sage: C.neighborliness()
            3
            sage: C = polytopes.cyclic_polytope(6,11); C
            A 6-dimensional polyhedron in QQ^6 defined as the convex hull of 11 vertices
            sage: C.neighborliness()
            3
            sage: [polytopes.cyclic_polytope(5,n).neighborliness() for n in range(6,10)]
            [6, 2, 2, 2]

        """
        return self.combinatorial_polyhedron().neighborliness()

    def is_neighborly(self, k=None):
        r"""
        Return whether the polyhedron is neighborly.

        If the input ``k`` is provided, then return whether the polyhedron is ``k``-neighborly

        A polyhedron is neighborly if every set of `n` vertices forms a face
        for `n` up to floor of half the dimension of the polyhedron.
        It is `k`-neighborly if this is true for `n` up to `k`.

        INPUT:

        - ``k`` -- the dimension up to which to check if every set of ``k``
          vertices forms a face. If no ``k`` is provided, check up to floor
          of half the dimension of the polyhedron.

        OUTPUT:

        - ``True`` if every set of up to ``k`` vertices forms a face,
        - ``False`` otherwise

        .. SEEALSO::

            :meth:`neighborliness`

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: cube.is_neighborly()
            True
            sage: cube = polytopes.hypercube(4)
            sage: cube.is_neighborly()
            False

        Cyclic polytopes are neighborly::

            sage: all(polytopes.cyclic_polytope(i, i + 1 + j).is_neighborly() for i in range(5) for j in range(3))
            True

        The neighborliness of a polyhedron equals floor of dimension half
        (or larger in case of a simplex) if and only if the polyhedron
        is neighborly::

            sage: testpolys = [polytopes.cube(), polytopes.cyclic_polytope(6, 9), polytopes.simplex(6)]
            sage: [(P.neighborliness()>=floor(P.dim()/2)) == P.is_neighborly() for P in  testpolys]
            [True, True, True]

        """
        return self.combinatorial_polyhedron().is_neighborly()

    @cached_method
    def is_lattice_polytope(self):
        r"""
        Return whether the polyhedron is a lattice polytope.

        OUTPUT:

        ``True`` if the polyhedron is compact and has only integral
        vertices, ``False`` otherwise.

        EXAMPLES::

            sage: polytopes.cross_polytope(3).is_lattice_polytope()
            True
            sage: polytopes.regular_polygon(5).is_lattice_polytope()
            False
        """
        if not self.is_compact():
            return False
        if self.base_ring() is ZZ:
            return True
        return all(v.is_integral() for v in self.vertex_generator())

    @cached_method
    def lattice_polytope(self, envelope=False):
        r"""
        Return an encompassing lattice polytope.

        INPUT:

        - ``envelope`` -- boolean (default: ``False``). If the
          polyhedron has non-integral vertices, this option decides
          whether to return a strictly larger lattice polytope or
          raise a ``ValueError``. This option has no effect if the
          polyhedron has already integral vertices.

        OUTPUT:

        A :class:`LatticePolytope
        <sage.geometry.lattice_polytope.LatticePolytopeClass>`. If the
        polyhedron is compact and has integral vertices, the lattice
        polytope equals the polyhedron. If the polyhedron is compact
        but has at least one non-integral vertex, a strictly larger
        lattice polytope is returned.

        If the polyhedron is not compact, a ``NotImplementedError`` is
        raised.

        If the polyhedron is not integral and ``envelope=False``, a
        ``ValueError`` is raised.

        ALGORITHM:

        For each non-integral vertex, a bounding box of integral
        points is added and the convex hull of these integral points
        is returned.

        EXAMPLES:

        First, a polyhedron with integral vertices::

            sage: P = Polyhedron( vertices = [(1, 0), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope(); lp
            2-d reflexive polytope #3 in 2-d lattice M
            sage: lp.vertices()
            M(-1,  0),
            M( 0, -1),
            M( 0,  1),
            M( 1,  0)
            in 2-d lattice M

        Here is a polyhedron with non-integral vertices::

            sage: P = Polyhedron( vertices = [(1/2, 1/2), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope()
            Traceback (most recent call last):
            ...
            ValueError: Some vertices are not integral. You probably want
            to add the argument "envelope=True" to compute an enveloping
            lattice polytope.
            sage: lp = P.lattice_polytope(True); lp
            2-d reflexive polytope #5 in 2-d lattice M
            sage: lp.vertices()
            M(-1,  0),
            M( 0, -1),
            M( 1,  1),
            M( 0,  1),
            M( 1,  0)
            in 2-d lattice M
        """
        if not self.is_compact():
            raise NotImplementedError('only compact lattice polytopes are allowed')

        try:
            vertices = self.vertices_matrix(ZZ).columns()
        except TypeError:
            if not envelope:
                raise ValueError('Some vertices are not integral. '
                    'You probably want to add the argument '
                    '"envelope=True" to compute an enveloping lattice polytope.')
            vertices = []
            for v in self.vertex_generator():
                vbox = [ set([floor(x), ceil(x)]) for x in v ]
                vertices.extend( itertools.product(*vbox) )

        # construct the (enveloping) lattice polytope
        from sage.geometry.lattice_polytope import LatticePolytope
        return LatticePolytope(vertices)

    def _integral_points_PALP(self):
        r"""
        Return the integral points in the polyhedron using PALP.

        This method is for testing purposes and will eventually be removed.

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])._integral_points_PALP()
            [M(-1, -1), M(0, 1), M(1, 0), M(1, 1), M(0, 0)]
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)]).lattice_polytope(True).points()
            M(-1, -1),
            M(-1,  0),
            M( 0, -1),
            M( 1,  1),
            M( 0,  1),
            M( 1,  0),
            M( 0,  0)
            in 2-d lattice M
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)])._integral_points_PALP()
            [M(1, 1), M(0, 1), M(1, 0), M(0, 0)]
        """
        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')
        lp = self.lattice_polytope(True)
        # remove cached values to get accurate timings
        try:
            del lp._points
            del lp._npoints
        except AttributeError:
            pass
        if self.is_lattice_polytope():
            return list(lp.points())
        return [p for p in lp.points() if self.contains(p)]

    @cached_method(do_pickle=True)
    def h_star_vector(self):
        r"""
        Return the `h^*`-vector of the lattice polytope.

        The `h^*`-vector records the coefficients of the polynomial in the
        numerator of the Ehrhart series of a lattice polytope.

        INPUT:

        - ``self`` -- A lattice polytope.

        OUTPUT:

        A list whose entries give the `h^*`-vector.

        .. NOTE:

            The backend of ``self`` should be ``'normaliz'``.
            This function depends on Normaliz (i.e. the ``'pynormaliz'`` optional
            package). See the Normaliz documentation for further details.

        EXAMPLES:

        The `h^*`-vector of a unimodular simplex S (a simplex with
        volume = `\frac{1}{dim(S)!}`) is always 1. Here we test this on
        simplices up to dimension 3::

            sage: s1 = polytopes.simplex(1,backend='normaliz')              # optional - pynormaliz
            sage: s2 = polytopes.simplex(2,backend='normaliz')              # optional - pynormaliz
            sage: s3 = polytopes.simplex(3,backend='normaliz')              # optional - pynormaliz
            sage: [s1.h_star_vector(),s2.h_star_vector(),s3.h_star_vector()]  # optional - pynormaliz
            [[1], [1], [1]]

        For a less trivial example, we compute the `h^*`-vector of the
        `0/1`-cube, which has the Eulerian numbers `(3,i)` for `i \in [0,2]`
        as an `h^*`-vector::

            sage: cube = polytopes.cube(intervals='zero_one', backend='normaliz') # optional - pynormaliz
            sage: cube.h_star_vector()   # optional - pynormaliz
            [1, 4, 1]
            sage: from sage.combinat.combinat import eulerian_number
            sage: [eulerian_number(3,i) for i in range(3)]
            [1, 4, 1]

        TESTS::

            sage: s3 = polytopes.simplex(3)
            sage: s3.h_star_vector()
            Traceback (most recent call last):
            ...
            TypeError: The backend of self must be normaliz

            sage: t = Polyhedron(vertices=[[0],[1/2]])
            sage: t.h_star_vector()
            Traceback (most recent call last):
            ...
            TypeError: The h_star vector is only defined for lattice polytopes

            sage: t2 = Polyhedron(vertices=[[AA(sqrt(2))],[1/2]])
            sage: t2.h_star_vector()
            Traceback (most recent call last):
            ...
            TypeError: The h_star vector is only defined for lattice polytopes

        Check that the `h^*`-vector is pickled::

            sage: new_cube = loads(dumps(cube))         # optional - pynormaliz
            sage: new_cube.h_star_vector.is_in_cache()  # optional - pynormaliz
            True
        """
        if self.is_empty():
            return 0
        if not self.is_lattice_polytope():
            raise TypeError('The h_star vector is only defined for lattice polytopes')
        if not self.backend() == 'normaliz':
            raise TypeError('The backend of self must be normaliz')
        return self._h_star_vector_normaliz()

    def _h_star_vector_normaliz(self):
        r"""
        Return the `h^*`-vector of a lattice polytope with backend = 'normaliz'.

        INPUT:

        - ``self`` -- A lattice polytope.

        OUTPUT:

        The `h^*`-vector as a list.

        .. NOTE:

        The backend of ``self`` should be ``'normaliz'``.

        TESTS::

            sage: s3 = polytopes.simplex(3)
            sage: s3._h_star_vector_normaliz()
            Traceback (most recent call last):
            ...
            TypeError: the backend should be normaliz
        """
        raise TypeError("the backend should be normaliz")

    @cached_method
    def bounding_box(self, integral=False, integral_hull=False):
        r"""
        Return the coordinates of a rectangular box containing the non-empty polytope.

        INPUT:

        - ``integral`` -- Boolean (default: ``False``). Whether to
          only allow integral coordinates in the bounding box.

        - ``integral_hull`` -- Boolean (default: ``False``). If ``True``, return a
          box containing the integral points of the polytope, or ``None, None`` if it
          is known that the polytope has no integral points.

        OUTPUT:

        A pair of tuples ``(box_min, box_max)`` where ``box_min`` are
        the coordinates of a point bounding the coordinates of the
        polytope from below and ``box_max`` bounds the coordinates
        from above.

        EXAMPLES::

            sage: Polyhedron([ (1/3,2/3), (2/3, 1/3) ]).bounding_box()
            ((1/3, 1/3), (2/3, 2/3))
            sage: Polyhedron([ (1/3,2/3), (2/3, 1/3) ]).bounding_box(integral=True)
            ((0, 0), (1, 1))
            sage: Polyhedron([ (1/3,2/3), (2/3, 1/3) ]).bounding_box(integral_hull=True)
            (None, None)
            sage: Polyhedron([ (1/3,2/3), (3/3, 4/3) ]).bounding_box(integral_hull=True)
            ((1, 1), (1, 1))
            sage: polytopes.buckyball(exact=False).bounding_box()
            ((-0.8090169944, -0.8090169944, -0.8090169944), (0.8090169944, 0.8090169944, 0.8090169944))

        TESTS::

            sage: Polyhedron().bounding_box()
            Traceback (most recent call last):
            ...
            ValueError: empty polytope is not allowed
        """
        box_min = []
        box_max = []
        if not self.is_compact():
            raise ValueError("only polytopes (compact polyhedra) are allowed")
        if self.n_vertices() == 0:
            raise ValueError("empty polytope is not allowed")
        for i in range(self.ambient_dim()):
            coords = [v[i] for v in self.vertex_generator()]
            max_coord = max(coords)
            min_coord = min(coords)
            if integral_hull:
                a = ceil(min_coord)
                b = floor(max_coord)
                if a > b:
                    return None, None
                box_max.append(b)
                box_min.append(a)
            elif integral:
                box_max.append(ceil(max_coord))
                box_min.append(floor(min_coord))
            else:
                box_max.append(max_coord)
                box_min.append(min_coord)
        return (tuple(box_min), tuple(box_max))

    def integral_points_count(self, **kwds):
        r"""
        Return the number of integral points in the polyhedron.

        This generic version of this method simply calls ``integral_points``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.integral_points_count()
            27

        We shrink the polyhedron a little bit::

            sage: Q = P*(8/9)
            sage: Q.integral_points_count()
            1

        Same for a polyhedron whose coordinates are not rationals.  Note that
        the answer is an integer even though there are no guarantees for
        exactness::

            sage: Q = P*RDF(8/9)
            sage: Q.integral_points_count()
            1

        Unbounded polyhedra (with or without lattice points) are not supported::

            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]])
            sage: P.integral_points_count()
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
            sage: P = Polyhedron(vertices=[[1, 1]], rays=[[1, 1]])
            sage: P.integral_points_count()
            Traceback (most recent call last):
            ...
            NotImplementedError: ...

        """
        return len(self.integral_points())

    def integral_points(self, threshold=100000):
        r"""
        Return the integral points in the polyhedron.

        Uses either the naive algorithm (iterate over a rectangular
        bounding box) or triangulation + Smith form.

        INPUT:

        - ``threshold`` -- integer (default: 100000). Use the naive
          algorithm as long as the bounding box is smaller than this.

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)]).integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = Polyhedron([(1,2,3), (2,3,7), (-2,-3,-11)])
            sage: simplex.integral_points()
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)])
            sage: simplex.integral_points()
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = Polyhedron([(2,3,7)])
            sage: point.integral_points()
            ((2, 3, 7),)

            sage: empty = Polyhedron()
            sage: empty.integral_points()
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = Polyhedron(v); simplex
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
            sage: len(simplex.integral_points())
            49

        A case where rounding in the right direction goes a long way::

            sage: P = 1/10*polytopes.hypercube(14, backend='field')
            sage: P.integral_points()
            ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),)

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ....:      (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v)
            sage: pts1 = P.integral_points()                     # Sage's own code
            sage: all(P.contains(p) for p in pts1)
            True
            sage: pts2 = LatticePolytope(v).points()          # PALP
            sage: for p in pts1: p.set_immutable()
            sage: set(pts1) == set(pts2)
            True

            sage: timeit('Polyhedron(v).integral_points()')   # not tested - random
            625 loops, best of 3: 1.41 ms per loop
            sage: timeit('LatticePolytope(v).points()')       # not tested - random
            25 loops, best of 3: 17.2 ms per loop

        TESTS:

        Test some trivial cases (see :trac:`17937`)::

            sage: P = Polyhedron(ambient_dim=1)  # empty polyhedron in 1 dimension
            sage: P.integral_points()
            ()
            sage: P = Polyhedron(ambient_dim=0)  # empty polyhedron in 0 dimensions
            sage: P.integral_points()
            ()
            sage: P = Polyhedron([[3]])  # single point in 1 dimension
            sage: P.integral_points()
            ((3),)
            sage: P = Polyhedron([[1/2]])  # single non-integral point in 1 dimension
            sage: P.integral_points()
            ()
            sage: P = Polyhedron([[]])  # single point in 0 dimensions
            sage: P.integral_points()
            ((),)

        Test unbounded polyhedron::

            sage: P = Polyhedron(rays=[[1,0,0]])
            sage: P.integral_points()
            Traceback (most recent call last):
            ...
            ValueError: can only enumerate points in a compact polyhedron
        """
        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')
        # Trivial cases: polyhedron with 0 or 1 vertices
        if self.n_vertices() == 0:
            return ()
        if self.n_vertices() == 1:
            v = self.vertices_list()[0]
            try:
                return (vector(ZZ, v),)
            except TypeError:  # vertex not integral
                return ()

        # for small bounding boxes, it is faster to naively iterate over the points of the box
        box_min, box_max = self.bounding_box(integral_hull=True)
        if box_min is None:
            return ()
        box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))
        if not self.is_lattice_polytope() or \
                (self.is_simplex() and box_points < 1000) or \
                box_points < threshold:
            from sage.geometry.integral_points import rectangular_box_points
            return rectangular_box_points(list(box_min), list(box_max), self)

        # for more complicate polytopes, triangulate & use smith normal form
        from sage.geometry.integral_points import simplex_points
        if self.is_simplex():
            return simplex_points(self.Vrepresentation())
        triangulation = self.triangulate()
        points = set()
        for simplex in triangulation:
            triang_vertices = [self.Vrepresentation(i) for i in simplex]
            new_points = simplex_points(triang_vertices)
            for p in new_points:
                p.set_immutable()
            points.update(new_points)
        # assert all(self.contains(p) for p in points)   # slow
        return tuple(points)

    def get_integral_point(self, index, **kwds):
        r"""
        Return the ``index``-th integral point in this polyhedron.

        This is equivalent to ``sorted(self.integral_points())[index]``.
        However, so long as self.integral_points_count() does not need to
        enumerate all integral points, neither does this method. Hence it can
        be significantly faster. If the polyhedron is not compact, a
        ``ValueError`` is raised.

        INPUT:

        - ``index`` -- integer. The index of the integral point to be found. If
          this is not in [0, ``self.integral_point_count()``), an ``IndexError``
          is raised.

        - ``**kwds`` -- optional keyword parameters that are passed to
          :meth:`self.integral_points_count`.

        ALGORITHM:

        The function computes each of the components of the requested point in
        turn. To compute x_i, the ith component, it bisects the upper and lower
        bounds on x_i given by the bounding box. At each bisection, it uses
        :meth:`integral_points_count` to determine on which side of the
        bisecting hyperplane the requested point lies.

        .. SEEALSO::

            :meth:`integral_points_count`.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])
            sage: P.get_integral_point(1)
            (0, 0)
            sage: P.get_integral_point(4)
            (1, 1)
            sage: sorted(P.integral_points())
            [(-1, -1), (0, 0), (0, 1), (1, 0), (1, 1)]
            sage: P.get_integral_point(5)
            Traceback (most recent call last):
            ...
            IndexError: ...

            sage: Q = Polyhedron([(1,3), (2, 7), (9, 77)])
            sage: [Q.get_integral_point(i) for i in range(Q.integral_points_count())] == sorted(Q.integral_points())
            True
            sage: Q.get_integral_point(0, explicit_enumeration_threshold=0, triangulation='cddlib')  # optional - latte_int
            (1, 3)
            sage: Q.get_integral_point(0, explicit_enumeration_threshold=0, triangulation='cddlib', foo=True)  # optional - latte_int
            Traceback (most recent call last):
            ...
            RuntimeError: ...

            sage: R = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]])
            sage: R.get_integral_point(0)
            Traceback (most recent call last):
            ...
            ValueError: ...
        """

        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')

        if not 0 <= index < self.integral_points_count(**kwds):
            raise IndexError('polytope index out of range')

        D = self.ambient_dim()
        lower_bounds, upper_bounds = self.bounding_box()
        coordinate = []
        P = self
        S = self.parent()
        for i in range(D):  # Now compute x_i, the ith component of coordinate.
            lower, upper = ceil(lower_bounds[i]), floor(upper_bounds[i]) + 1  # So lower <= x_i < upper.
            while lower < upper-1:
                guess = (lower + upper) // 2  # > lower.
                # Build new polyhedron by intersecting P with the halfspace {x_i < guess}.
                P_lt_guess = P.intersection(S(None, ([[guess-1] + [0] * i + [-1] + [0] * (D - i - 1)], [])))
                # Avoid computing P_geq_guess = P.intersection({x_i >= guess}) right now, it might not be needed.
                P_lt_guess_count = P_lt_guess.integral_points_count(**kwds)
                if P_lt_guess_count > index:  # Move upper down to guess.
                    upper = guess
                    index -= 0
                    P = P_lt_guess
                else:  # P_lt_guess_count <= index:  # Move lower up to guess.
                    lower = guess
                    index -= P_lt_guess_count
                    P_geq_guess = P.intersection(S(None, ([[-guess] + [0] * i + [1] + [0] * (D - i - 1)], [])))
                    P = P_geq_guess
            coordinate.append(lower)  # Record the new component that we have found.
        point = vector(ZZ, coordinate)
        point.set_immutable()
        return point

    def random_integral_point(self, **kwds):
        r"""
        Return an integral point in this polyhedron chosen uniformly at random.

        INPUT:

        - ``**kwds`` -- optional keyword parameters that are passed to
          :meth:`self.get_integral_point`.

        OUTPUT:

        The integral point in the polyhedron chosen uniformly at random. If the
        polyhedron is not compact, a ``ValueError`` is raised. If the
        polyhedron does not contain any integral points, an ``EmptySetError`` is
        raised.

        .. SEEALSO::

            :meth:`get_integral_point`.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])
            sage: P.random_integral_point()  # random
            (0, 0)
            sage: P.random_integral_point() in P.integral_points()
            True
            sage: P.random_integral_point(explicit_enumeration_threshold=0, triangulation='cddlib')  # random, optional - latte_int
            (1, 1)
            sage: P.random_integral_point(explicit_enumeration_threshold=0, triangulation='cddlib', foo=7)  # optional - latte_int
            Traceback (most recent call last):
            ...
            RuntimeError: ...

            sage: Q = Polyhedron(vertices=[(2, 1/3)], rays=[(1, 2)])
            sage: Q.random_integral_point()
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: R = Polyhedron(vertices=[(1/2, 0), (1, 1/2), (0, 1/2)])
            sage: R.random_integral_point()
            Traceback (most recent call last):
            ...
            EmptySetError: ...
        """

        if not self.is_compact():
            raise ValueError('can only sample integral points in a compact polyhedron')

        count = self.integral_points_count()
        if count == 0:
            raise EmptySetError('polyhedron does not contain any integral points')

        return self.get_integral_point(current_randstate().python_random().randint(0, count-1), **kwds)

    @cached_method
    def combinatorial_automorphism_group(self, vertex_graph_only=False):
        """
        Computes the combinatorial automorphism group.

        If ``vertex_graph_only`` is ``True``,  the automorphism group
        of the vertex-edge graph of the polyhedron is returned. Otherwise
        the automorphism group of the vertex-facet graph, which is
        isomorphic to the automorphism group of the face lattice is returned.

        INPUT:

        - ``vertex_graph_only`` -- boolean (default: ``False``); whether
          to return the automorphism group of the vertex edges graph or
          of the lattice

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic_with_category'>`
        that is isomorphic to the combinatorial automorphism group is
        returned.

        - if ``vertex_graph_only`` is ``True``:
          The automorphism group of the vertex-edge graph of the polyhedron

        - if ``vertex_graph_only`` is ``False`` (default):
          The automorphism group of the vertex-facet graph of the polyhedron,
          see :meth:`vertex_facet_graph`. This group is isomorphic to the
          automorphism group of the face lattice of the polyhedron.

        NOTE:

            Depending on ``vertex_graph_only``, this method returns groups
            that are not necessarily isomorphic, see the examples below.

        .. SEEALSO::

            :meth:`is_combinatorially_isomorphic`,
            :meth:`graph`,
            :meth:`vertex_facet_graph`.

        EXAMPLES::

            sage: quadrangle = Polyhedron(vertices=[(0,0),(1,0),(0,1),(2,3)])
            sage: quadrangle.combinatorial_automorphism_group().is_isomorphic(groups.permutation.Dihedral(4))
            True
            sage: quadrangle.restricted_automorphism_group()
            Permutation Group with generators [()]

        Permutations of the vertex graph only exchange vertices with vertices::

            sage: P = Polyhedron(vertices=[(1,0), (1,1)], rays=[(1,0)])
            sage: P.combinatorial_automorphism_group(vertex_graph_only=True)
            Permutation Group with generators [(A vertex at (1,0),A vertex at (1,1))]

        This shows an example of two polytopes whose vertex-edge graphs are isomorphic,
        but their face_lattices are not isomorphic::

            sage: Q=Polyhedron([[-123984206864/2768850730773, -101701330976/922950243591, -64154618668/2768850730773, -2748446474675/2768850730773],
            ....: [-11083969050/98314591817, -4717557075/98314591817, -32618537490/98314591817, -91960210208/98314591817],
            ....: [-9690950/554883199, -73651220/554883199, 1823050/554883199, -549885101/554883199], [-5174928/72012097, 5436288/72012097, -37977984/72012097, 60721345/72012097],
            ....: [-19184/902877, 26136/300959, -21472/902877, 899005/902877], [53511524/1167061933, 88410344/1167061933, 621795064/1167061933, 982203941/1167061933],
            ....: [4674489456/83665171433, -4026061312/83665171433, 28596876672/83665171433, -78383796375/83665171433], [857794884940/98972360190089, -10910202223200/98972360190089, 2974263671400/98972360190089, -98320463346111/98972360190089]])
            sage: C = polytopes.cyclic_polytope(4,8)
            sage: C.is_combinatorially_isomorphic(Q)
            False
            sage: C.combinatorial_automorphism_group(vertex_graph_only=True).is_isomorphic(Q.combinatorial_automorphism_group(vertex_graph_only=True))
            True
            sage: C.combinatorial_automorphism_group(vertex_graph_only=False).is_isomorphic(Q.combinatorial_automorphism_group(vertex_graph_only=False))
            False

        The automorphism group of the face lattice is isomorphic to the combinatorial automorphism group::

            sage: CG = C.hasse_diagram().automorphism_group()
            sage: C.combinatorial_automorphism_group().is_isomorphic(CG)
            True
            sage: QG = Q.hasse_diagram().automorphism_group()
            sage: Q.combinatorial_automorphism_group().is_isomorphic(QG)
            True

        """
        if vertex_graph_only:
            G = self.graph()
        else:
            G = self.vertex_facet_graph()
        return G.automorphism_group(edge_labels=True)

    @cached_method
    def restricted_automorphism_group(self, output="abstract"):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the affine group `AGL(d,\RR) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional polyhedron. The affine group
        acts in the usual way `\vec{x}\mapsto A\vec{x}+b` on the
        ambient space.

        The restricted automorphism group is the subgroup of the linear
        automorphism group generated by permutations of the generators
        of the same type. That is, vertices can only be permuted with
        vertices, ray generators with ray generators, and line
        generators with line generators.

        For example, take the first quadrant

        .. MATH::

            Q = \Big\{ (x,y) \Big| x\geq 0,\; y\geq0 \Big\}
            \subset \QQ^2

        Then the linear automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            a & 0 \\ 0 & b
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & c \\ d & 0
            \end{pmatrix}
            :~
            a, b, c, d \in \QQ_{>0}
            \right\}
            \subset
            GL(2,\QQ)
            \subset
            E(d)

        Note that there are no translations that map the quadrant `Q`
        to itself, so the linear automorphism group is contained in
        the general linear group (the subgroup of transformations
        preserving the origin). The restricted automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            1 & 0 \\ 0 & 1
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & 1 \\ 1 & 0
            \end{pmatrix}
            \right\}
            \simeq \ZZ_2

        INPUT:

        - ``output`` -- how the group should be represented:

          - ``"abstract"`` (default) -- return an abstract permutation
            group without further meaning.

          - ``"permutation"`` -- return a permutation group on the
            indices of the polyhedron generators. For example, the
            permutation ``(0,1)`` would correspond to swapping
            ``self.Vrepresentation(0)`` and ``self.Vrepresentation(1)``.

          - ``"matrix"`` -- return a matrix group representing affine
            transformations. When acting on affine vectors, you should
            append a `1` to every vector. If the polyhedron is not full
            dimensional, the returned matrices act as the identity on
            the orthogonal complement of the affine space spanned by
            the polyhedron.

          - ``"matrixlist"`` -- like ``matrix``, but return the list of
            elements of the matrix group. Useful for fields without a
            good implementation of matrix groups or to avoid the
            overhead of creating the group.

        OUTPUT:

        - For ``output="abstract"`` and ``output="permutation"``:
          a :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`.

        - For ``output="matrix"``: a :class:`MatrixGroup`.

        - For ``output="matrixlist"``: a list of matrices.

        REFERENCES:

        - [BSS2009]_

        EXAMPLES:

        A cross-polytope example::

            sage: P = polytopes.cross_polytope(3)
            sage: P.restricted_automorphism_group() == PermutationGroup([[(3,4)], [(2,3),(4,5)],[(2,5)],[(1,2),(5,6)],[(1,6)]])
            True
            sage: P.restricted_automorphism_group(output="permutation") == PermutationGroup([[(2,3)],[(1,2),(3,4)],[(1,4)],[(0,1),(4,5)],[(0,5)]])
            True
            sage: mgens = [[[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,1]], [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], [[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]]

        We test groups for equality in a fool-proof way; they can have different generators, etc::

            sage: poly_g = P.restricted_automorphism_group(output="matrix")
            sage: matrix_g = MatrixGroup([matrix(QQ,t) for t in mgens])
            sage: all(t.matrix() in poly_g for t in matrix_g.gens())
            True
            sage: all(t.matrix() in matrix_g for t in poly_g.gens())
            True

        24-cell example::

            sage: P24 = polytopes.twenty_four_cell()
            sage: AutP24 = P24.restricted_automorphism_group()
            sage: PermutationGroup([
            ....:     '(1,20,2,24,5,23)(3,18,10,19,4,14)(6,21,11,22,7,15)(8,12,16,17,13,9)',
            ....:     '(1,21,8,24,4,17)(2,11,6,15,9,13)(3,20)(5,22)(10,16,12,23,14,19)'
            ....: ]).is_isomorphic(AutP24)
            True
            sage: AutP24.order()
            1152

        Here is the quadrant example mentioned in the beginning::

            sage: P = Polyhedron(rays=[(1,0),(0,1)])
            sage: P.Vrepresentation()
            (A vertex at (0, 0), A ray in the direction (0, 1), A ray in the direction (1, 0))
            sage: P.restricted_automorphism_group(output="permutation")
            Permutation Group with generators [(1,2)]

        Also, the polyhedron need not be full-dimensional::

            sage: P = Polyhedron(vertices=[(1,2,3,4,5),(7,8,9,10,11)])
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(1,2)]
            sage: G = P.restricted_automorphism_group(output="matrixlist")
            sage: G
            (
            [1 0 0 0 0 0]  [ -87/55  -82/55    -2/5   38/55   98/55   12/11]
            [0 1 0 0 0 0]  [-142/55  -27/55    -2/5   38/55   98/55   12/11]
            [0 0 1 0 0 0]  [-142/55  -82/55     3/5   38/55   98/55   12/11]
            [0 0 0 1 0 0]  [-142/55  -82/55    -2/5   93/55   98/55   12/11]
            [0 0 0 0 1 0]  [-142/55  -82/55    -2/5   38/55  153/55   12/11]
            [0 0 0 0 0 1], [      0       0       0       0       0       1]
            )
            sage: g = AffineGroup(5, QQ)(G[1])
            sage: g
                  [ -87/55  -82/55    -2/5   38/55   98/55]     [12/11]
                  [-142/55  -27/55    -2/5   38/55   98/55]     [12/11]
            x |-> [-142/55  -82/55     3/5   38/55   98/55] x + [12/11]
                  [-142/55  -82/55    -2/5   93/55   98/55]     [12/11]
                  [-142/55  -82/55    -2/5   38/55  153/55]     [12/11]
            sage: g^2
                  [1 0 0 0 0]     [0]
                  [0 1 0 0 0]     [0]
            x |-> [0 0 1 0 0] x + [0]
                  [0 0 0 1 0]     [0]
                  [0 0 0 0 1]     [0]
            sage: g(list(P.vertices()[0]))
            (7, 8, 9, 10, 11)
            sage: g(list(P.vertices()[1]))
            (1, 2, 3, 4, 5)

        Affine transformations do not change the restricted automorphism
        group. For example, any non-degenerate triangle has the
        dihedral group with 6 elements, `D_6`, as its automorphism
        group::

            sage: initial_points = [vector([1,0]), vector([0,1]), vector([-2,-1])]
            sage: points = initial_points
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[0] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - 2*initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]

        The ``output="matrixlist"`` can be used over fields without a
        complete implementation of matrix groups::

            sage: P = polytopes.dodecahedron(); P
            A 3-dimensional polyhedron in (Number Field in sqrt5 with defining polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^3 defined as the convex hull of 20 vertices
            sage: G = P.restricted_automorphism_group(output="matrixlist")
            sage: len(G)
            120

        Floating-point computations are supported with a simple fuzzy
        zero implementation::

            sage: P = Polyhedron(vertices=[(1/3,0,0,1),(0,1/4,0,1),(0,0,1/5,1)], base_ring=RDF)
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: len(P.restricted_automorphism_group(output="matrixlist"))
            6

        TESTS::

            sage: P = Polyhedron(vertices=[(1,0), (1,1)], rays=[(1,0)])
            sage: P.restricted_automorphism_group(output="permutation")
            Permutation Group with generators [(1,2)]
            sage: P.restricted_automorphism_group(output="matrix")
            Matrix group over Rational Field with 1 generators (
            [ 1  0  0]
            [ 0 -1  1]
            [ 0  0  1]
            )
            sage: P.restricted_automorphism_group(output="foobar")
            Traceback (most recent call last):
            ...
            ValueError: unknown output 'foobar', valid values are ('abstract', 'permutation', 'matrix', 'matrixlist')

        Check that :trac:`28828` is fixed::

            sage: P.restricted_automorphism_group(output="matrixlist")[0].is_immutable()
            True
        """
        # The algorithm works as follows:
        #
        # Let V be the matrix where every column is a homogeneous
        # coordinate of a V-representation object (vertex, ray, line).
        # Let us assume that V has full rank, that the polyhedron is
        # full dimensional.
        #
        # Let Q = V Vt and C = Vt Q^-1 V. The rows and columns of C
        # can be thought of as being indexed by the V-rep objects of the
        # polytope.
        #
        # It turns out that we can identify the restricted automorphism
        # group with the automorphism group of the edge-colored graph
        # on the V-rep objects with colors determined by the symmetric
        # matrix C.
        #
        # An automorphism of this graph is equivalent to a permutation
        # matrix P such that C = Pt C P. If we now define
        # A = V P Vt Q^-1, then one can check that V P = A V.
        # In other words: permuting the generators is the same as
        # applying the affine transformation A on the generators.
        #
        # If the given polyhedron is not fully-dimensional,
        # then Q will be not invertible. In this case, we use a
        # pseudoinverse Q+ instead of Q^-1. The formula for A acting on
        # the space spanned by V then simplifies to A = V P V+ where V+
        # denotes the pseudoinverse of V, which also equals V+ = Vt Q+.
        #
        # If we are asked to return the (group of) transformation
        # matrices to the user, we also require that those
        # transformations act as the identity on the orthogonal
        # complement of the space spanned by V. This complement is the
        # space spanned by the columns of W = 1 - V V+. One can check
        # that B = (V P V+) + W is the correct matrix: it acts the same
        # as A on V and it satisfies B W = W.

        outputs = ("abstract", "permutation", "matrix", "matrixlist")
        if output not in outputs:
            raise ValueError("unknown output {!r}, valid values are {}".format(output, outputs))

        # For backwards compatibility, we treat "abstract" as
        # "permutation", but where we add 1 to the indices of the
        # permutations.
        index0 = 0
        if output == "abstract":
            index0 = 1
            output = "permutation"

        if self.base_ring().is_exact():
            def rational_approximation(c):
                return c
        else:
            c_list = []

            def rational_approximation(c):
                # Implementation detail: Return unique integer if two
                # c-values are the same up to machine precision. But
                # you can think of it as a uniquely-chosen rational
                # approximation.
                for i, x in enumerate(c_list):
                    if self._is_zero(x - c):
                        return i
                c_list.append(c)
                return len(c_list) - 1

        if self.is_compact():
            def edge_label(i, j, c_ij):
                return c_ij
        else:
            # In the non-compact case, we also label the edges by the
            # type of the V-representation object. This ensures that
            # vertices, rays, and lines are only permuted amongst
            # themselves.
            def edge_label(i, j, c_ij):
                return (self.Vrepresentation(i).type(), c_ij, self.Vrepresentation(j).type())

        # Homogeneous coordinates for the V-representation objects.
        # Mathematically, V is a matrix. For efficiency however, we
        # represent it as a list of column vectors.
        V = [v.homogeneous_vector() for v in self.Vrepresentation()]

        # Pseudoinverse of V Vt
        Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()

        # Construct the graph.
        G = Graph()
        for i in range(len(V)):
            for j in range(i+1, len(V)):
                c_ij = rational_approximation(V[i] * Qplus * V[j])
                G.add_edge(index0+i, index0+j, edge_label(i, j, c_ij))

        permgroup = G.automorphism_group(edge_labels=True)
        if output == "permutation":
            return permgroup
        elif output == "matrix":
            permgroup = permgroup.gens()

        # Compute V+ = Vt Q+ as list of row vectors
        Vplus = list(matrix(V) * Qplus)  # matrix(V) is Vt

        # Compute W = 1 - V V+
        W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

        # Convert the permutation group to a matrix group.
        # If P is a permutation, then we return the matrix
        # B = (V P V+) + W.
        #
        # If output == "matrix", we loop over the generators of the group.
        # Otherwise, we loop over all elements.
        matrices = []
        for perm in permgroup:
            A = sum(V[perm(i)].column() * Vplus[i].row() for i in range(len(V)))
            matrices.append(A + W)

        for mat in matrices:
            mat.set_immutable()

        if output == "matrixlist":
            return tuple(matrices)
        else:
            return MatrixGroup(matrices)

    def is_combinatorially_isomorphic(self, other, algorithm='bipartite_graph'):
        r"""
        Return whether the polyhedron is combinatorially isomorphic to another polyhedron.

        We only consider bounded polyhedra. By definition, they are
        combinatorially isomorphic if their face lattices are isomorphic.

        INPUT:

        - ``other`` -- a polyhedron object
        - ``algorithm`` (default = ``bipartite_graph``) -- the algorithm to use.
          The other possible value is ``face_lattice``.

        OUTPUT:

        - ``True`` if the two polyhedra are combinatorially isomorphic
        - ``False`` otherwise

        .. SEEALSO::

            :meth:`combinatorial_automorphism_group`,
            :meth:`vertex_facet_graph`.

        REFERENCES:

        For the equivalence of the two algorithms see [KK1995]_, p. 877-878

        EXAMPLES:

        The square is combinatorially isomorphic to the 2-dimensional cube::

            sage: polytopes.hypercube(2).is_combinatorially_isomorphic(polytopes.regular_polygon(4))
            True

        All the faces of the 3-dimensional permutahedron are either
        combinatorially isomorphic to a square or a hexagon::

            sage: H = polytopes.regular_polygon(6)
            sage: S = polytopes.hypercube(2)
            sage: P = polytopes.permutahedron(4)
            sage: all(F.as_polyhedron().is_combinatorially_isomorphic(S) or F.as_polyhedron().is_combinatorially_isomorphic(H) for F in P.faces(2))
            True

        Checking that a regular simplex intersected with its reflection
        through the origin is combinatorially isomorphic to the intersection
        of a cube with a hyperplane perpendicular to its long diagonal::

            sage: def simplex_intersection(k):
            ....:   S1 = Polyhedron([vector(v)-vector(polytopes.simplex(k).center()) for v in polytopes.simplex(k).vertices_list()])
            ....:   S2 = Polyhedron([-vector(v) for v in S1.vertices_list()])
            ....:   return S1.intersection(S2)
            sage: def cube_intersection(k):
            ....:    C = polytopes.hypercube(k+1)
            ....:    H = Polyhedron(eqns=[[0]+[1 for i in range(k+1)]])
            ....:    return C.intersection(H)
            sage: [simplex_intersection(k).is_combinatorially_isomorphic(cube_intersection(k)) for k in range(2,5)]
            [True, True, True]
            sage: simplex_intersection(2).is_combinatorially_isomorphic(polytopes.regular_polygon(6))
            True
            sage: simplex_intersection(3).is_combinatorially_isomorphic(polytopes.octahedron())
            True

        Two polytopes with the same `f`-vector, but different combinatorial types::

            sage: P = Polyhedron([[-605520/1525633, -605520/1525633, -1261500/1525633, -52200/1525633, 11833/1525633],\
             [-720/1769, -600/1769, 1500/1769, 0, -31/1769], [-216/749, 240/749, -240/749, -432/749, 461/749], \
             [-50/181, 50/181, 60/181, -100/181, -119/181], [-32/51, -16/51, -4/51, 12/17, 1/17],\
             [1, 0, 0, 0, 0], [16/129, 128/129, 0, 0, 1/129], [64/267, -128/267, 24/89, -128/267, 57/89],\
             [1200/3953, -1200/3953, -1440/3953, -360/3953, -3247/3953], [1512/5597, 1512/5597, 588/5597, 4704/5597, 2069/5597]])
            sage: C = polytopes.cyclic_polytope(5,10)
            sage: C.f_vector() == P.f_vector(); C.f_vector()
            True
            (1, 10, 45, 100, 105, 42, 1)
            sage: C.is_combinatorially_isomorphic(P)
            False

            sage: S = polytopes.simplex(3)
            sage: S = S.face_truncation(S.faces(0)[3])
            sage: S = S.face_truncation(S.faces(0)[4])
            sage: S = S.face_truncation(S.faces(0)[5])
            sage: T = polytopes.simplex(3)
            sage: T = T.face_truncation(T.faces(0)[3])
            sage: T = T.face_truncation(T.faces(0)[4])
            sage: T = T.face_truncation(T.faces(0)[4])
            sage: T.is_combinatorially_isomorphic(S)
            False
            sage: T.f_vector(), S.f_vector()
            ((1, 10, 15, 7, 1), (1, 10, 15, 7, 1))

            sage: C = polytopes.hypercube(5)
            sage: C.is_combinatorially_isomorphic(C)
            True
            sage: C.is_combinatorially_isomorphic(C, algorithm='magic')
            Traceback (most recent call last):
            ...
            AssertionError: `algorithm` must be 'bipartite graph' or 'face_lattice'

            sage: G = Graph()
            sage: C.is_combinatorially_isomorphic(G)
            Traceback (most recent call last):
            ...
            AssertionError: input `other` must be a polyhedron

            sage: H = Polyhedron(eqns=[[0,1,1,1,1]]); H
            A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex and 3 lines
            sage: C.is_combinatorially_isomorphic(H)
            Traceback (most recent call last):
            ...
            AssertionError: polyhedron `other` must be bounded

        """
        assert isinstance(other, Polyhedron_base), "input `other` must be a polyhedron"
        assert self.is_compact(), "polyhedron `self` must be bounded"
        assert other.is_compact(), "polyhedron `other` must be bounded"
        assert algorithm in ['bipartite_graph', 'face_lattice'], "`algorithm` must be 'bipartite graph' or 'face_lattice'"

        # For speed, we check if the polyhedra have the same number of facets and vertices.
        # This is faster than building the bipartite graphs first and
        # then check that they won't be isomorphic.
        if self.n_vertices() != other.n_vertices() or self.n_facets() != other.n_facets():
            return False

        if algorithm == 'bipartite_graph':
            G_self = self.vertex_facet_graph(False)
            G_other = other.vertex_facet_graph(False)

            return G_self.is_isomorphic(G_other)
        else:
            return self.face_lattice().is_isomorphic(other.face_lattice())

    def _test_is_combinatorially_isomorphic(self, tester=None, **options):
        """
        Run tests on the method :meth:`.is_combinatorially_isomorphic`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_is_combinatorially_isomorphic()
        """
        if tester is None:
            tester = self._tester(**options)

        if not self.is_compact():
            with tester.assertRaises(AssertionError):
                self.is_combinatorially_isomorphic(self)
            return

        if self.n_vertices() > 200 or self.n_facets() > 200:
            # Avoid very long doctests.
            return

        tester.assertTrue(self.is_combinatorially_isomorphic(ZZ(4)*self))
        if self.n_vertices():
            tester.assertTrue(self.is_combinatorially_isomorphic(self + self.center()))

        if self.n_vertices() < 20 and self.n_facets() < 20 and self.is_immutable():
            tester.assertTrue(self.is_combinatorially_isomorphic(ZZ(4)*self, algorithm='face_lattice'))
            if self.n_vertices():
                tester.assertTrue(self.is_combinatorially_isomorphic(self + self.center(), algorithm='face_lattice'))

    def affine_hull(self, *args, **kwds):
        r"""
        Return the affine hull of ``self`` as a polyhedron.

        EXAMPLES::

            sage: half_plane_in_space = Polyhedron(ieqs=[(0,1,0,0)], eqns=[(0,0,0,1)])
            sage: half_plane_in_space.affine_hull().Hrepresentation()
            (An equation (0, 0, 1) x + 0 == 0,)

            sage: polytopes.cube().affine_hull().is_universe()
            True
        """
        if args or kwds:
            raise TypeError("the method 'affine_hull' does not take any parameters; perhaps you meant 'affine_hull_projection'")
        if not self.inequalities():
            return self
        self_as_face = self.faces(self.dimension())[0]
        return self_as_face.affine_tangent_cone()

    @cached_method
    def _affine_hull_projection(self, *,
                                as_convex_set=True, as_affine_map=True, as_section_map=True,
                                orthogonal=False, orthonormal=False,
                                extend=False, minimal=False):
        r"""
        Return ``self`` projected into its affine hull.

        INPUT:

        See :meth:`affine_hull_projection`.

        OUTPUT:

        An instance of :class:`~sage.geometry.convex_set.AffineHullProjectionData`.
        See :meth:`affine_hull_projection` for details.

        TESTS:

        Check that :trac:`23355` is fixed::

            sage: P = Polyhedron([[7]]); P
            A 0-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex
            sage: P.affine_hull_projection()
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: P.affine_hull_projection(orthonormal='True')
            A 0-dimensional polyhedron in QQ^0 defined as the convex hull of 1 vertex
            sage: P.affine_hull_projection(orthogonal='True')
            A 0-dimensional polyhedron in QQ^0 defined as the convex hull of 1 vertex

        Check that :trac:`24047` is fixed::

            sage: P1 = Polyhedron(vertices=([[-1, 1], [0, -1], [0, 0], [-1, -1]]))
            sage: P2 = Polyhedron(vertices=[[1, 1], [1, -1], [0, -1], [0, 0]])
            sage: P = P1.intersection(P2)
            sage: A, b = P.affine_hull_projection(as_affine_map=True, orthonormal=True, extend=True)

            sage: Polyhedron([(2,3,4)]).affine_hull_projection()
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex

        Check that backend is preserved::

            sage: polytopes.simplex(backend='field').affine_hull_projection().backend()
            'field'

            sage: P = Polyhedron(vertices=[[0,0], [1,0]], backend='field')
            sage: P.affine_hull_projection(orthogonal=True, orthonormal=True, extend=True).backend()
            'field'

        Check that :trac:`29116` is fixed::

            sage: V =[
            ....:    [1, 0, -1, 0, 0],
            ....:    [1, 0, 0, -1, 0],
            ....:    [1, 0, 0, 0, -1],
            ....:    [1, 0, 0, +1, 0],
            ....:    [1, 0, 0, 0, +1],
            ....:    [1, +1, 0, 0, 0]
            ....:     ]
            sage: P = Polyhedron(V)
            sage: P.affine_hull_projection()
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 6 vertices
            sage: P.affine_hull_projection(orthonormal=True)
            Traceback (most recent call last):
            ...
            ValueError: the base ring needs to be extended; try with "extend=True"
            sage: P.affine_hull_projection(orthonormal=True, extend=True)
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 6 vertices
        """
        result = AffineHullProjectionData()

        if self.is_empty():
            raise ValueError('affine hull projection of an empty polyhedron is undefined')

        # handle trivial full-dimensional case
        if self.ambient_dim() == self.dim():
            if as_convex_set:
                result.image = self
            if as_affine_map:
                identity = linear_transformation(matrix(self.base_ring(),
                                                        self.dim(),
                                                        self.dim(),
                                                        self.base_ring().one()))
                result.projection_linear_map = result.section_linear_map = identity
                result.projection_translation = result.section_translation = self.ambient_space().zero()
        elif orthogonal or orthonormal:
            # see TODO
            if not self.is_compact():
                raise NotImplementedError('"orthogonal=True" and "orthonormal=True" work only for compact polyhedra')
            affine_basis = self.an_affine_basis()
            v0 = affine_basis[0].vector()
            # We implicitly translate the first vertex of the affine basis to zero.
            vi = tuple(v.vector() - v0 for v in affine_basis[1:])
            M = matrix(self.base_ring(), self.dim(), self.ambient_dim(), vi)

            # Switch base_ring to AA if necessary,
            # since gram_schmidt needs to be able to take square roots.
            # Pick orthonormal basis and transform all vertices accordingly
            # if the orthonormal transform makes it necessary, change base ring.
            try:
                A, G = M.gram_schmidt(orthonormal=orthonormal)
            except TypeError:
                if not extend:
                    raise ValueError('the base ring needs to be extended; try with "extend=True"')
                M = matrix(AA, M)
                A = M.gram_schmidt(orthonormal=orthonormal)[0]
                if minimal:
                    from sage.rings.qqbar import number_field_elements_from_algebraics
                    new_ring = number_field_elements_from_algebraics(A.list(), embedded=True, minimal=True)[0]
                    A = A.change_ring(new_ring)
            L = linear_transformation(A, side='right')
            ambient_translation = -vector(A.base_ring(), affine_basis[0])
            image_translation = A * ambient_translation
            # Note the order. We compute ``A*self`` and then translate the image.
            # ``A*self`` uses the incidence matrix and we avoid recomputation.
            # Also, if the new base ring is ``AA``, we want to avoid computing the incidence matrix in that ring.
            # ``convert=True`` takes care of the case, where there might be no coercion (``AA`` and quadratic field).
            if as_convex_set:
                result.image = self.linear_transformation(A, new_base_ring=A.base_ring()) + image_translation
            if as_affine_map:
                result.projection_linear_map = L
                result.projection_translation = image_translation
            if as_section_map:
                L_dagger = linear_transformation(A.transpose() * (A * A.transpose()).inverse(), side='right')
                result.section_linear_map = L_dagger
                result.section_translation = v0.change_ring(A.base_ring())
        else:
            # translate one vertex to the origin
            v0 = self.vertices()[0].vector()
            gens = []
            for v in self.vertices()[1:]:
                gens.append(v.vector() - v0)
            for r in self.rays():
                gens.append(r.vector())
            for l in self.lines():
                gens.append(l.vector())

            # Pick subset of coordinates to coordinatize the affine span
            M = matrix(gens)
            pivots = M.pivots()

            A = matrix(self.base_ring(), len(pivots), self.ambient_dim(),
                       [[1 if j == i else 0 for j in range(self.ambient_dim())] for i in pivots])
            if as_affine_map:
                image_translation = vector(self.base_ring(), self.dim())
                L = linear_transformation(A, side='right')
                result.projection_linear_map = L
                result.projection_translation = image_translation
            if as_convex_set:
                result.image = A*self
            if as_section_map:
                if self.dim():
                    B = M.transpose()/(A*M.transpose())
                else:
                    B = matrix(self.ambient_dim(), 0)
                L_section = linear_transformation(B, side='right')
                result.section_linear_map = L_section
                result.section_translation = v0 - L_section(L(v0) + image_translation)

        return result

    def affine_hull_projection(self,
                               as_polyhedron=None, as_affine_map=False,
                               orthogonal=False, orthonormal=False,
                               extend=False, minimal=False,
                               return_all_data=False,
                               *, as_convex_set=None):
        r"""Return the polyhedron projected into its affine hull.

        Each polyhedron is contained in some smallest affine subspace
        (possibly the entire ambient space) -- its affine hull.  We
        provide an affine linear map that projects the ambient space of
        the polyhedron to the standard Euclidean space of dimension of
        the polyhedron, which restricts to a bijection from the affine
        hull.

        The projection map is not unique; some parameters control the
        choice of the map.  Other parameters control the output of the
        function.

        INPUT:

        - ``as_polyhedron`` (or ``as_convex_set``) -- (boolean or the default
          ``None``) and

        - ``as_affine_map`` -- (boolean, default ``False``) control the output

          The default ``as_polyhedron=None`` translates to
          ``as_polyhedron=not as_affine_map``,
          therefore to ``as_polyhedron=True`` if nothing is specified.

          If exactly one of either ``as_polyhedron`` or ``as_affine_map`` is
          set, then either a polyhedron or the affine transformation
          is returned. The affine transformation
          sends the embedded polytope to a fulldimensional one.
          It is given as a pair ``(A, b)``, where A is a linear transformation
          and `b` is a vector, and the affine transformation sends ``v`` to
          ``A(v)+b``.

          If both ``as_polyhedron`` and ``as_affine_map`` are set, then
          both are returned, encapsulated in an instance of
          :class:`~sage.geometry.convex_set.AffineHullProjectionData`.

        - ``return_all_data`` -- (boolean, default ``False``)

          If set, then ``as_polyhedron`` and ``as_affine_map`` will set
          (possibly overridden) and additional (internal) data concerning
          the transformation is returned. Everything is encapsulated
          in an instance of
          :class:`~sage.geometry.convex_set.AffineHullProjectionData` in
          this case.

        - ``orthogonal`` -- boolean (default: ``False``); if ``True``,
          provide an orthogonal transformation.

        - ``orthonormal`` -- boolean (default: ``False``); if ``True``,
          provide an orthonormal transformation. If the base ring does not
          provide the necessary square roots, the extend parameter
          needs to be set to ``True``.

        - ``extend`` -- boolean (default: ``False``); if ``True``,
          allow base ring to be extended if necessary. This becomes
          relevant when requiring an orthonormal transformation.

        - ``minimal`` -- boolean (default: ``False``); if ``True``,
          when doing an extension, it computes the minimal base ring of the
          extension, otherwise the base ring is ``AA``.

        OUTPUT:

        A full-dimensional polyhedron or an affine transformation,
        depending on the parameters ``as_polyhedron`` and ``as_affine_map``,
        or an instance of :class:`~sage.geometry.convex_set.AffineHullProjectionData`
        containing all data (parameter ``return_all_data``).

        If the output is an instance of
        :class:`~sage.geometry.convex_set.AffineHullProjectionData`, the
        following fields may be set:

        - ``image`` -- the projection of the original polyhedron

        - ``projection_map`` -- the affine map as a pair whose first component
          is a linear transformation and its second component a shift;
          see above.

        - ``section_map`` -- an affine map as a pair whose first component
          is a linear transformation and its second component a shift.
          It maps the codomain of ``affine_map`` to the affine hull of
          ``self``.  It is a right inverse of ``projection_map``.

        Note that all of these data are compatible.

         .. TODO::

            - make the parameters ``orthogonal`` and ``orthonormal`` work
              with unbounded polyhedra.

        EXAMPLES::

            sage: triangle = Polyhedron([(1,0,0), (0,1,0), (0,0,1)]);  triangle
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: triangle.affine_hull_projection()
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices

            sage: half3d = Polyhedron(vertices=[(3,2,1)], rays=[(1,0,0)])
            sage: half3d.affine_hull_projection().Vrepresentation()
            (A ray in the direction (1), A vertex at (3))

        The resulting affine hulls depend on the parameter ``orthogonal`` and ``orthonormal``::

            sage: L = Polyhedron([[1,0],[0,1]]); L
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: A = L.affine_hull_projection(); A
            A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (0), A vertex at (1))
            sage: A = L.affine_hull_projection(orthogonal=True); A
            A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (0), A vertex at (2))
            sage: A = L.affine_hull_projection(orthonormal=True)
            Traceback (most recent call last):
            ...
            ValueError: the base ring needs to be extended; try with "extend=True"
            sage: A = L.affine_hull_projection(orthonormal=True, extend=True); A
            A 1-dimensional polyhedron in AA^1 defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (1.414213562373095?), A vertex at (0.?e-18))

        More generally::

            sage: S = polytopes.simplex(); S
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
            sage: S.vertices()
            (A vertex at (0, 0, 0, 1),
             A vertex at (0, 0, 1, 0),
             A vertex at (0, 1, 0, 0),
             A vertex at (1, 0, 0, 0))
            sage: A = S.affine_hull_projection(); A
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: A.vertices()
            (A vertex at (0, 0, 0),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 0),
             A vertex at (1, 0, 0))
            sage: A = S.affine_hull_projection(orthogonal=True); A
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: A.vertices()
            (A vertex at (0, 0, 0),
             A vertex at (2, 0, 0),
             A vertex at (1, 3/2, 0),
             A vertex at (1, 1/2, 4/3))
            sage: A = S.affine_hull_projection(orthonormal=True, extend=True); A
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 4 vertices
            sage: A.vertices()
            (A vertex at (0.7071067811865475?, 0.4082482904638630?, 1.154700538379252?),
             A vertex at (0.7071067811865475?, 1.224744871391589?, 0.?e-18),
             A vertex at (1.414213562373095?, 0.?e-18, 0.?e-18),
             A vertex at (0.?e-18, 0.?e-18, 0.?e-18))

        With the parameter ``minimal`` one can get a minimal base ring::

            sage: s = polytopes.simplex(3)
            sage: s_AA = s.affine_hull_projection(orthonormal=True, extend=True)
            sage: s_AA.base_ring()
            Algebraic Real Field
            sage: s_full = s.affine_hull_projection(orthonormal=True, extend=True, minimal=True)
            sage: s_full.base_ring()
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a = 0.5176380902050415?

        More examples with the ``orthonormal`` parameter::

            sage: P = polytopes.permutahedron(3); P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: set([F.as_polyhedron().affine_hull_projection(orthonormal=True, extend=True).volume() for F in P.affine_hull_projection().faces(1)]) == {1, sqrt(AA(2))}
            True
            sage: set([F.as_polyhedron().affine_hull_projection(orthonormal=True, extend=True).volume() for F in P.affine_hull_projection(orthonormal=True, extend=True).faces(1)]) == {sqrt(AA(2))}
            True
            sage: D = polytopes.dodecahedron()
            sage: F = D.faces(2)[0].as_polyhedron()
            sage: F.affine_hull_projection(orthogonal=True)
            A 2-dimensional polyhedron in (Number Field in sqrt5 with defining polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^2 defined as the convex hull of 5 vertices
            sage: F.affine_hull_projection(orthonormal=True, extend=True)
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 5 vertices
            sage: K.<sqrt2> = QuadraticField(2)
            sage: P = Polyhedron([2*[K.zero()],2*[sqrt2]])
            sage: K.<sqrt2> = QuadraticField(2)
            sage: P = Polyhedron([2*[K.zero()],2*[sqrt2]]); P
            A 1-dimensional polyhedron in (Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?)^2 defined as the convex hull of 2 vertices
            sage: P.vertices()
            (A vertex at (0, 0), A vertex at (sqrt2, sqrt2))
            sage: A = P.affine_hull_projection(orthonormal=True); A
            A 1-dimensional polyhedron in (Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?)^1 defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (0), A vertex at (2))
            sage: K.<sqrt3> = QuadraticField(3)
            sage: P = Polyhedron([2*[K.zero()],2*[sqrt3]]); P
            A 1-dimensional polyhedron in (Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?)^2 defined as the convex hull of 2 vertices
            sage: P.vertices()
            (A vertex at (0, 0), A vertex at (sqrt3, sqrt3))
            sage: A = P.affine_hull_projection(orthonormal=True)
            Traceback (most recent call last):
            ...
            ValueError: the base ring needs to be extended; try with "extend=True"
            sage: A = P.affine_hull_projection(orthonormal=True, extend=True); A
            A 1-dimensional polyhedron in AA^1 defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (0), A vertex at (2.449489742783178?))
            sage: sqrt(6).n()
            2.44948974278318

        The affine hull is combinatorially equivalent to the input::

            sage: P.is_combinatorially_isomorphic(P.affine_hull_projection())
            True
            sage: P.is_combinatorially_isomorphic(P.affine_hull_projection(orthogonal=True))
            True
            sage: P.is_combinatorially_isomorphic(P.affine_hull_projection(orthonormal=True, extend=True))
            True

        The ``orthonormal=True`` parameter preserves volumes;
        it provides an isometric copy of the polyhedron::

            sage: Pentagon = polytopes.dodecahedron().faces(2)[0].as_polyhedron()
            sage: P = Pentagon.affine_hull_projection(orthonormal=True, extend=True)
            sage: _, c= P.is_inscribed(certificate=True)
            sage: c
            (0.4721359549995794?, 0.6498393924658126?)
            sage: circumradius = (c-vector(P.vertices()[0])).norm()
            sage: p = polytopes.regular_polygon(5)
            sage: p.volume()
            2.377641290737884?
            sage: P.volume()
            1.53406271079097?
            sage: p.volume()*circumradius^2
            1.534062710790965?
            sage: P.volume() == p.volume()*circumradius^2
            True

        One can also use ``orthogonal`` parameter to calculate volumes;
        in this case we don't need to switch base rings. One has to divide
        by the square root of the determinant of the linear part of the
        affine transformation times its transpose::

            sage: Pentagon = polytopes.dodecahedron().faces(2)[0].as_polyhedron()
            sage: Pnormal = Pentagon.affine_hull_projection(orthonormal=True, extend=True)
            sage: Pgonal = Pentagon.affine_hull_projection(orthogonal=True)
            sage: A, b = Pentagon.affine_hull_projection(orthogonal=True, as_affine_map=True)
            sage: Adet = (A.matrix().transpose()*A.matrix()).det()
            sage: Pnormal.volume()
            1.53406271079097?
            sage: Pgonal.volume()/Adet.sqrt(extend=True)
            -80*(55*sqrt(5) - 123)/sqrt(-6368*sqrt(5) + 14240)
            sage: Pgonal.volume()/AA(Adet).sqrt().n(digits=20)
            1.5340627107909646813
            sage: AA(Pgonal.volume()^2) == (Pnormal.volume()^2)*AA(Adet)
            True

        Another example with ``as_affine_map=True``::

            sage: P = polytopes.permutahedron(4)
            sage: A, b = P.affine_hull_projection(orthonormal=True, as_affine_map=True, extend=True)
            sage: Q = P.affine_hull_projection(orthonormal=True, extend=True)
            sage: Q.center()
            (0.7071067811865475?, 1.224744871391589?, 1.732050807568878?)
            sage: A(P.center()) + b == Q.center()
            True

        For unbounded, non full-dimensional polyhedra, the ``orthogonal=True`` and ``orthonormal=True``
        is not implemented::

            sage: P = Polyhedron(ieqs=[[0, 1, 0], [0, 0, 1], [0, 0, -1]]); P
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: P.is_compact()
            False
            sage: P.is_full_dimensional()
            False
            sage: P.affine_hull_projection(orthogonal=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: "orthogonal=True" and "orthonormal=True" work only for compact polyhedra
            sage: P.affine_hull_projection(orthonormal=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: "orthogonal=True" and "orthonormal=True" work only for compact polyhedra

        Setting ``as_affine_map`` to ``True``
        without ``orthogonal`` or ``orthonormal`` set to ``True``::

            sage: S = polytopes.simplex()
            sage: S.affine_hull_projection(as_affine_map=True)
            (Vector space morphism represented by the matrix:
             [1 0 0]
             [0 1 0]
             [0 0 1]
             [0 0 0]
             Domain: Vector space of dimension 4 over Rational Field
             Codomain: Vector space of dimension 3 over Rational Field,
             (0, 0, 0))

        If the polyhedron is full-dimensional, it is returned::

            sage: polytopes.cube().affine_hull_projection()
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: polytopes.cube().affine_hull_projection(as_affine_map=True)
            (Vector space morphism represented by the matrix:
             [1 0 0]
             [0 1 0]
             [0 0 1]
             Domain: Vector space of dimension 3 over Rational Field
             Codomain: Vector space of dimension 3 over Rational Field,
             (0, 0, 0))

        Return polyhedron and affine map::

            sage: S = polytopes.simplex(2)
            sage: data = S.affine_hull_projection(orthogonal=True,
            ....:                                 as_polyhedron=True,
            ....:                                 as_affine_map=True); data
            AffineHullProjectionData(image=A 2-dimensional polyhedron in QQ^2
                    defined as the convex hull of 3 vertices,
                projection_linear_map=Vector space morphism represented by the matrix:
                    [  -1 -1/2]
                    [   1 -1/2]
                    [   0    1]
                    Domain: Vector space of dimension 3 over Rational Field
                    Codomain: Vector space of dimension 2 over Rational Field,
                projection_translation=(1, 1/2),
                section_linear_map=None,
                section_translation=None)

        Return all data::

            sage: data = S.affine_hull_projection(orthogonal=True, return_all_data=True); data
            AffineHullProjectionData(image=A 2-dimensional polyhedron in QQ^2
                    defined as the convex hull of 3 vertices,
                projection_linear_map=Vector space morphism represented by the matrix:
                    [  -1 -1/2]
                    [   1 -1/2]
                    [   0    1]
                    Domain: Vector space of dimension 3 over Rational Field
                    Codomain: Vector space of dimension 2 over Rational Field,
                projection_translation=(1, 1/2),
                section_linear_map=Vector space morphism represented by the matrix:
                    [-1/2  1/2    0]
                    [-1/3 -1/3  2/3]
                    Domain: Vector space of dimension 2 over Rational Field
                    Codomain: Vector space of dimension 3 over Rational Field, section_translation=(1, 0, 0))

        The section map is a right inverse of the projection map::

            sage: data.image.linear_transformation(data.section_linear_map.matrix().transpose()) + data.section_translation == S
            True

        Same without ``orthogonal=True``::

            sage: data = S.affine_hull_projection(return_all_data=True); data
            AffineHullProjectionData(image=A 2-dimensional polyhedron in ZZ^2
                    defined as the convex hull of 3 vertices,
                projection_linear_map=Vector space morphism represented by the matrix:
                    [1 0]
                    [0 1]
                    [0 0]
                    Domain: Vector space of dimension 3 over Rational Field
                    Codomain: Vector space of dimension 2 over Rational Field, projection_translation=(0, 0),
                section_linear_map=Vector space morphism represented by the matrix:
                    [ 1  0 -1]
                    [ 0  1 -1]
                    Domain: Vector space of dimension 2 over Rational Field
                    Codomain: Vector space of dimension 3 over Rational Field, section_translation=(0, 0, 1))
            sage: data.image.linear_transformation(data.section_linear_map.matrix().transpose()) + data.section_translation == S
            True

        ::

            sage: P0 = Polyhedron(
            ....:     ieqs=[(0, -1, 0, 1, 1, 1), (0, 1, 1, 0, -1, -1), (0, -1, 1, 1, 0, 0),
            ....:           (0, 1, 0, 0, 0, 0), (0, 0, 1, 1, -1, -1), (0, 0, 0, 0, 0, 1),
            ....:           (0, 0, 0, 0, 1, 0), (0, 0, 0, 1, 0, -1), (0, 0, 1, 0, 0, 0)])
            sage: P = P0.intersection(Polyhedron(eqns=[(-1, 1, 1, 1, 1, 1)]))
            sage: P.dim()
            4
            sage: P.affine_hull_projection(orthogonal=True, as_affine_map=True)[0]
            Vector space morphism represented by the matrix:
            [    0     0     0   1/3]
            [ -2/3  -1/6     0 -1/12]
            [  1/3  -1/6   1/2 -1/12]
            [    0   1/2     0 -1/12]
            [  1/3  -1/6  -1/2 -1/12]
            Domain: Vector space of dimension 5 over Rational Field
            Codomain: Vector space of dimension 4 over Rational Field
        """
        if as_polyhedron is not None:
            as_convex_set = as_polyhedron
        return super().affine_hull_projection(
            as_convex_set=as_convex_set, as_affine_map=as_affine_map,
            orthogonal=orthogonal, orthonormal=orthonormal,
            extend=extend, minimal=minimal,
            return_all_data=return_all_data)

    def _test_affine_hull_projection(self, tester=None, verbose=False, **options):
        """
        Run tests on the method :meth:`.affine_hull_projection`.

        TESTS::

            sage: D = polytopes.dodecahedron()
            sage: D.facets()[0].as_polyhedron()._test_affine_hull_projection()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.is_empty():
            # Undefined, nothing to test
            return

        if self.n_vertices() > 30 or self.n_facets() > 30 or self.dim() > 6:
            # Avoid very long doctests.
            return

        data_sets = [None]*4
        data_sets[0] = self.affine_hull_projection(return_all_data=True)
        if self.is_compact():
            data_sets[1] = self.affine_hull_projection(return_all_data=True,
                                                       orthogonal=True,
                                                       extend=True)
            data_sets[2] = self.affine_hull_projection(return_all_data=True,
                                                       orthonormal=True,
                                                       extend=True)
            data_sets[3] = self.affine_hull_projection(return_all_data=True,
                                                       orthonormal=True,
                                                       extend=True,
                                                       minimal=True)
        else:
            data_sets = data_sets[:1]

        for i, data in enumerate(data_sets):
            if verbose:
                print("Running test number {}".format(i))
            M = data.projection_linear_map.matrix().transpose()
            tester.assertEqual(self.linear_transformation(M, new_base_ring=M.base_ring())
                               + data.projection_translation,
                               data.image)

            M = data.section_linear_map.matrix().transpose()
            if M.base_ring() is AA:
                self_extend = self.change_ring(AA)
            else:
                self_extend = self
            tester.assertEqual(data.image.linear_transformation(M)
                               + data.section_translation,
                               self_extend)
            if i == 0:
                tester.assertEqual(data.image.base_ring(), self.base_ring())
            else:
                # Test whether the map is orthogonal.
                M = data.projection_linear_map.matrix()
                tester.assertTrue((M.transpose() * M).is_diagonal())
                if i > 1:
                    # Test whether the map is orthonormal.
                    tester.assertTrue((M.transpose() * M).is_one())
            if i == 3:
                # Test that the extension is indeed minimal.
                if self.base_ring() is not AA:
                    tester.assertIsNot(data.image.base_ring(), AA)

    def affine_hull_manifold(self, name=None, latex_name=None, start_index=0, ambient_space=None,
                             ambient_chart=None, names=None, **kwds):
        r"""
        Return the affine hull of ``self`` as a manifold.

        If ``self`` is full-dimensional, it is just the ambient Euclidean space.
        Otherwise, it is a Riemannian submanifold of the ambient Euclidean space.

        INPUT:

        - ``ambient_space`` -- a :class:`~sage.manifolds.differentiable.examples.euclidean.EuclideanSpace`
          of the ambient dimension (default: the manifold of ``ambient_chart``, if provided;
          otherwise, a new instance of ``EuclideanSpace``).

        - ``ambient_chart`` -- a chart on ``ambient_space``.

        - ``names`` -- names for the coordinates on the affine hull.

        - optional arguments accepted by :meth:`~sage.geometry.polyhedron.base.affine_hull_projection`.

        The default chart is determined by the optional arguments of
        :meth:`~sage.geometry.polyhedron.base.affine_hull_projection`.

        EXAMPLES::

            sage: triangle = Polyhedron([(1,0,0), (0,1,0), (0,0,1)]);  triangle
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: A = triangle.affine_hull_manifold(name='A'); A
            2-dimensional Riemannian submanifold A embedded in the Euclidean space E^3
            sage: A.embedding().display()
            A → E^3
               (x0, x1) ↦ (x, y, z) = (t0 + x0, t0 + x1, t0 - x0 - x1 + 1)
            sage: A.embedding().inverse().display()
            E^3 → A
               (x, y, z) ↦ (x0, x1) = (x, y)
            sage: A.adapted_chart()
            [Chart (E^3, (x0_E3, x1_E3, t0_E3))]
            sage: A.normal().display()
            n = 1/3*sqrt(3) e_x + 1/3*sqrt(3) e_y + 1/3*sqrt(3) e_z
            sage: A.induced_metric()       # Need to call this before volume_form
            Riemannian metric gamma on the 2-dimensional Riemannian submanifold A embedded in the Euclidean space E^3
            sage: A.volume_form()
            2-form eps_gamma on the 2-dimensional Riemannian submanifold A embedded in the Euclidean space E^3

        Orthogonal version::

            sage: A = triangle.affine_hull_manifold(name='A', orthogonal=True); A
            2-dimensional Riemannian submanifold A embedded in the Euclidean space E^3
            sage: A.embedding().display()
            A → E^3
               (x0, x1) ↦ (x, y, z) = (t0 - 1/2*x0 - 1/3*x1 + 1, t0 + 1/2*x0 - 1/3*x1, t0 + 2/3*x1)
            sage: A.embedding().inverse().display()
            E^3 → A
               (x, y, z) ↦ (x0, x1) = (-x + y + 1, -1/2*x - 1/2*y + z + 1/2)

        Arrangement of affine hull of facets::

            sage: D = polytopes.dodecahedron()
            sage: E3 = EuclideanSpace(3)
            sage: submanifolds = [
            ....:     F.as_polyhedron().affine_hull_manifold(name=f'F{i}', orthogonal=True, ambient_space=E3)
            ....:     for i, F in enumerate(D.facets())]
            sage: sum(FM.plot({}, srange(-2, 2, 0.1), srange(-2, 2, 0.1), opacity=0.2)  # not tested
            ....:     for FM in submanifolds) + D.plot()
            Graphics3d Object

        Full-dimensional case::

            sage: cube = polytopes.cube(); cube
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: cube.affine_hull_manifold()
            Euclidean space E^3

        """
        if ambient_space is None:
            if ambient_chart is not None:
                ambient_space = ambient_chart.manifold()
            else:
                from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
                ambient_space = EuclideanSpace(self.ambient_dim(), start_index=start_index)
        if ambient_space.dimension() != self.ambient_dim():
            raise ValueError('ambient_space and ambient_chart must match the ambient dimension')

        if self.is_full_dimensional():
            return ambient_space

        if ambient_chart is None:
            ambient_chart = ambient_space.default_chart()
        CE = ambient_chart

        from sage.manifolds.manifold import Manifold
        if name is None:
            name, latex_name = self._affine_hull_name_latex_name()
        H = Manifold(self.dim(), name, ambient=ambient_space, structure="Riemannian",
                     latex_name=latex_name, start_index=start_index)
        if names is None:
            names = tuple(f'x{i}' for i in range(self.dim()))
        CH = H.chart(names=names)

        data = self.affine_hull_projection(return_all_data=True, **kwds)
        projection_matrix = data.projection_linear_map.matrix().transpose()
        projection_translation_vector = data.projection_translation
        section_matrix = data.section_linear_map.matrix().transpose()
        section_translation_vector = data.section_translation

        from sage.symbolic.ring import SR
        # We use the slacks of the (linear independent) equations as the foliation parameters
        foliation_parameters = vector(SR.var(f't{i}') for i in range(self.ambient_dim() - self.dim()))
        normal_matrix = matrix(equation.A() for equation in self.equation_generator()).transpose()
        slack_matrix = normal_matrix.pseudoinverse()

        phi = H.diff_map(ambient_space, {(CH, CE):
                                         (section_matrix * vector(CH._xx) + section_translation_vector
                                          + normal_matrix * foliation_parameters).list()})
        phi_inv = ambient_space.diff_map(H, {(CE, CH):
                                             (projection_matrix * vector(CE._xx) + projection_translation_vector).list()})

        foliation_scalar_fields = {parameter:
                                   ambient_space.scalar_field({CE: slack_matrix.row(i) * (vector(CE._xx) - section_translation_vector)})
                                   for i, parameter in enumerate(foliation_parameters)}

        H.set_embedding(phi, inverse=phi_inv,
                        var=list(foliation_parameters), t_inverse=foliation_scalar_fields)
        return H

    def _affine_hull_name_latex_name(self, name=None, latex_name=None):
        r"""
        Return the default name of the affine hull.

        EXAMPLES::

            sage: polytopes.cube()._affine_hull_name_latex_name('C', r'\square')
            ('aff_C', '\\mathop{\\mathrm{aff}}(\\square)')

            sage: Polyhedron(vertices=[[0, 1], [1, 0]])._affine_hull_name_latex_name()
            ('aff_P', '\\mathop{\\mathrm{aff}}(P)')
        """

        if name is None:
            name = 'P'
        if latex_name is None:
            latex_name = name
        operator = 'aff'
        aff_name = f'{operator}_{name}'
        aff_latex_name = r'\mathop{\mathrm{' + operator + '}}(' + latex_name + ')'
        return aff_name, aff_latex_name

    def _polymake_init_(self):
        """
        Return a polymake "Polytope" object corresponding to ``self``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: PP = polymake(P)         # optional - polymake
            sage: PP.N_VERTICES            # optional - polymake
            8

        Lower-dimensional polyhedron::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1]])
            sage: PP = polymake(P)         # optional - polymake
            sage: PP.COMBINATORIAL_DIM     # optional - polymake
            1
            sage: PP.AFFINE_HULL           # optional - polymake
            -1 1 1

        Empty polyhedron::

            sage: P = Polyhedron(ambient_dim=2, vertices=[])
            sage: PP = polymake(P)         # optional - polymake
            sage: PP.COMBINATORIAL_DIM     # optional - polymake
            -1

        Pointed unbounded polyhedron::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1]], rays=[[1, 0]])
            sage: PP = polymake(P)         # optional - polymake
            sage: PP.VERTICES              # optional - polymake
            1 0 1
            1 1 0
            0 1 0
            sage: PP.FACETS                # optional - polymake
            1 0 -1
            -1 1 1
            0 0 1

        Non-pointed polyhedron::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1]], lines=[[1, 0]])
            sage: PP = polymake(P)         # optional - polymake
            sage: PP.VERTICES              # optional - polymake
            1 0 1
            1 0 0
            sage: PP.FACETS                # optional - polymake
            1 0 -1
            0 0 1
            sage: PP.LINEALITY_SPACE       # optional - polymake
            0 1 0

        Algebraic polyhedron::

            sage: P = polytopes.dodecahedron(); P
            A 3-dimensional polyhedron in (Number Field in sqrt5 with defining polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^3 defined as the convex hull of 20 vertices
            sage: print("There may be a recompilation warning"); PP = polymake(P); PP # optional - polymake
            There may be a recompilation warning...
            Polytope<QuadraticExtension<Rational>>[...]
            sage: sorted(PP.VERTICES[:], key=repr)[0]  # optional - polymake
            1 -1+1r5 -4+2r5 0

        Floating-point polyhedron::

            sage: P = polytopes.dodecahedron(exact=False); P
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 20 vertices
            sage: print("There may be a recompilation warning"); PP = polymake(P); PP # optional - polymake
            There may be a recompilation warning...
            Polytope<Float>[...]
            sage: sorted(PP.VERTICES[:], key=repr)[0] # optional - polymake
            1 -0.472135955 0 -1.236067978

        """
        from sage.interfaces.polymake import polymake
        polymake_field = polymake(self.base_ring().fraction_field())
        polymake_class = "Polytope<{}>".format(polymake_field)
        if self.is_empty():
            # Polymake 3.1 cannot enter an empty polyhedron using
            # FACETS and AFFINE_HULL.  Use corresponding input properties instead.
            # https://forum.polymake.org/viewtopic.php?f=8&t=545
            return polymake.new_object(polymake_class,
                                       INEQUALITIES=self.inequalities_list(),
                                       EQUATIONS=self.equations_list())
        else:
            return polymake.new_object(polymake_class,
                                       FACETS=self.inequalities_list(),
                                       AFFINE_HULL=self.equations_list(),
                                       VERTICES=   [ [1] + v for v in self.vertices_list() ] \
                                                 + [ [0] + r for r in self.rays_list() ],
                                       LINEALITY_SPACE=[ [0] + l for l in self.lines_list() ])
