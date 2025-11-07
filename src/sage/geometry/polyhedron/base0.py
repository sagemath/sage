r"""
Base class for polyhedra: Initialization and access to Vrepresentation and Hrepresentation
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
#       Copyright (C) 2019-2021 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.structure.element import Element
import sage.geometry.abc


class Polyhedron_base0(Element, sage.geometry.abc.Polyhedron):
    """
    Initialization and basic access for polyhedra.

    See :class:`sage.geometry.polyhedron.base.Polyhedron_base`.

    TESTS::

        sage: from sage.geometry.polyhedron.base0 import Polyhedron_base0
        sage: P = Polyhedron(rays=[[1, 0, 0]], lines=[[0, 1, 0]])
        sage: Polyhedron_base0.Vrepresentation(P)
        (A line in the direction (0, 1, 0),
        A vertex at (0, 0, 0),
        A ray in the direction (1, 0, 0))
        sage: Polyhedron_base0.vertices.f(P)
        (A vertex at (0, 0, 0),)
        sage: Polyhedron_base0.rays.f(P)
        (A ray in the direction (1, 0, 0),)
        sage: Polyhedron_base0.lines.f(P)
        (A line in the direction (0, 1, 0),)
        sage: Polyhedron_base0.Hrepresentation(P)
        (An equation (0, 0, 1) x + 0 == 0, An inequality (1, 0, 0) x + 0 >= 0)
        sage: Polyhedron_base0.inequalities.f(P)
        (An inequality (1, 0, 0) x + 0 >= 0,)
        sage: Polyhedron_base0.equations.f(P)
        (An equation (0, 0, 1) x + 0 == 0,)
        sage: Polyhedron_base0.base_ring(P)
        Integer Ring
        sage: Polyhedron_base0.backend(P)
        'ppl'
        sage: Polyhedron_base0.change_ring(P, ZZ, backend='field').backend()
        'field'
        sage: Polyhedron_base0.base_extend(P, QQ)
        A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 1 ray, 1 line
    """
    def __init__(self, parent, Vrep, Hrep, Vrep_minimal=None, Hrep_minimal=None, pref_rep=None, mutable=False, **kwds):
        """
        Initialize the polyhedron.

        See :class:`sage.geometry.polyhedron.base.Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron()    # indirect doctests

            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: from sage.geometry.polyhedron.parent import Polyhedra_field
            sage: parent = Polyhedra_field(AA, 1, 'field')                              # needs sage.rings.number_field
            sage: Vrep = [[[0], [1/2], [1]], [], []]
            sage: Hrep = [[[0, 1], [1, -1]], []]
            sage: p = Polyhedron_field(parent, Vrep, Hrep,                              # needs sage.rings.number_field
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

            sage: p = Polyhedron_field(parent, Vrep, 'nonsense',
            ....:                      Vrep_minimal=True, Hrep_minimal=True, pref_rep='Vrep')
            Traceback (most recent call last):
            ...
            TypeError: ..._init_Hrepresentation() takes 3 positional arguments but 9 were given

        The empty polyhedron is detected when the Vrepresentation is given with generator;
        see :issue:`29899`::

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

    @abstract_method
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
            TypeError: 'AbstractMethod' object is not callable
        """

    @abstract_method
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
            TypeError: 'AbstractMethod' object is not callable
        """

    def _init_empty_polyhedron(self):
        """
        Initialize an empty polyhedron.

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
            sage: P = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]], backend='polymake') # optional - jupymake
            sage: sage_input(P)                                                                    # optional - jupymake
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

    def base_extend(self, base_ring, backend=None):
        """
        Return a new polyhedron over a larger base ring.

        This method can also be used to change the backend.

        INPUT:

        - ``base_ring`` -- the new base ring

        - ``backend`` -- the new backend, see
          :func:`~sage.geometry.polyhedron.constructor.Polyhedron`
          If ``None`` (the default), attempt to keep the same backend.
          Otherwise, use the same defaulting behavior
          as described there.

        OUTPUT: the same polyhedron, but over a larger base ring and possibly with a changed backend

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)], base_ring=ZZ);  P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.base_extend(QQ)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.base_extend(QQ) == P
            True

        TESTS:

        Test that :issue:`22575` is fixed::

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

            sage: P = polytopes.regular_polygon(3); P                                   # needs sage.rings.number_field
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 3 vertices
            sage: P.vertices()                                                          # needs sage.rings.number_field
            (A vertex at (0.?e-16, 1.000000000000000?),
             A vertex at (0.866025403784439?, -0.500000000000000?),
             A vertex at (-0.866025403784439?, -0.500000000000000?))
            sage: P.change_ring(QQ)                                                     # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: cannot change the base ring to the Rational Field

        .. WARNING::

            The base ring ``RDF`` should be used with care. As it is
            not an exact ring, certain computations may break or
            silently produce wrong results, for example changing the
            base ring from an exact ring into ``RDF`` may cause a
            loss of data::

                sage: P = Polyhedron([[2/3,0],[6666666666666667/10^16,0]], base_ring=AA); P         # needs sage.rings.number_field
                A 1-dimensional polyhedron in AA^2 defined as the convex hull of 2 vertices
                sage: Q = P.change_ring(RDF); Q                                         # needs sage.rings.number_field
                A 0-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex
                sage: P.n_vertices() == Q.n_vertices()                                  # needs sage.rings.number_field
                False
        """
        from sage.categories.rings import Rings

        if base_ring not in Rings():
            raise ValueError("invalid base ring")

        try:
            vertices = [[base_ring(x) for x in vertex] for vertex in self.vertices_list()]
            rays = [[base_ring(x) for x in ray] for ray in self.rays_list()]
            lines = [[base_ring(x) for x in line] for line in self.lines_list()]

        except (TypeError, ValueError):
            raise TypeError("cannot change the base ring to the {0}".format(base_ring))

        new_parent = self.parent().change_ring(base_ring, backend)
        return new_parent([vertices, rays, lines], None)

    def is_mutable(self):
        r"""
        Return ``True`` if the polyhedron is mutable, i.e. it can be modified in place.

        EXAMPLES::

            sage: p = polytopes.cube(backend='field')
            sage: p.is_mutable()
            False
        """
        return False

    def is_immutable(self):
        r"""
        Return ``True`` if the polyhedron is immutable, i.e. it cannot be modified in place.

        EXAMPLES::

            sage: p = polytopes.cube(backend='field')
            sage: p.is_immutable()
            True
        """
        return True

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

    def is_compact(self):
        """
        Test for boundedness of the polytope.

        EXAMPLES::

            sage: p = polytopes.icosahedron()                                           # needs sage.groups sage.rings.number_field
            sage: p.is_compact()                                                        # needs sage.groups sage.rings.number_field
            True
            sage: p = Polyhedron(ieqs=[[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,-1,0,0]])
            sage: p.is_compact()
            False
        """
        return self.n_rays() == 0 and self.n_lines() == 0

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

        - ``separator`` -- string (default: ``'\n'``)

        - ``latex`` -- boolean (default: ``False``)

        - ``style`` -- either ``'positive'`` (making all coefficients positive)
                       or ``'<='``, or ``'>='``; default is ``'>='``

        - ``align`` -- boolean or ``None''; default is ``None`` in which case
                       ``align`` is ``True`` if ``separator`` is the newline character.
                       If set, then the lines of the output string are aligned
                       by the comparison symbol by padding blanks.

        Keyword parameters of
        :meth:`~sage.geometry.polyhedron.representation.Hrepresentation.repr_pretty`
        are passed on:

        - ``prefix`` -- string

        - ``indices`` -- tuple or other iterable

        OUTPUT: string

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
            '1 >= x0, 1 >= x1, 1 >= x2, 1 + x0 >= 0, 1 + x2 >= 0, 1 + x1 >= 0'
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
        yield from self.Hrepresentation()

    @cached_method
    def n_Hrepresentation(self):
        """
        Return the number of objects that make up the
        H-representation of the polyhedron.

        OUTPUT: integer

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

        OUTPUT: integer

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
        yield from self.Vrepresentation()

    def inequality_generator(self):
        """
        Return  a generator for the defining inequalities of the
        polyhedron.

        OUTPUT: a generator of the inequality Hrepresentation objects

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

        OUTPUT: a tuple of inequalities

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[1,0,0],[2,2,2]])
            sage: p.inequalities()[0:3]
            (An inequality (1, 0, 0) x + 0 >= 0,
             An inequality (0, 1, 0) x + 0 >= 0,
             An inequality (0, 0, 1) x + 0 >= 0)

            sage: # needs sage.combinat
            sage: p3 = Polyhedron(vertices=Permutations([1, 2, 3, 4]))
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
            list of :class:`~sage.geometry.polyhedron.representation.Inequality` objects.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[1,0,0],[2,2,2]])
            sage: p.inequalities_list()[0:3]
            [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

            sage: # needs sage.combinat
            sage: p3 = Polyhedron(vertices=Permutations([1, 2, 3, 4]))
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

        OUTPUT: a tuple of equations

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
            iterate over the list of :class:`~sage.geometry.polyhedron.representation.Vertex` objects.

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

        OUTPUT: a tuple of vertices

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

        - ``base_ring`` -- a ring or ``None`` (default); the base ring
          of the returned matrix. If not specified, the base ring of
          the polyhedron is used.

        OUTPUT:

        A matrix over ``base_ring`` whose columns are the coordinates
        of the vertices. A :exc:`TypeError` is raised if the coordinates
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

        Check that :issue:`28828` is fixed::

                sage: P.vertices_matrix().is_immutable()
                True
        """
        from sage.matrix.constructor import matrix

        if base_ring is None:
            base_ring = self.base_ring()
        m = matrix(base_ring, self.ambient_dim(), self.n_vertices())
        for i, v in enumerate(self.vertices()):
            for j in range(self.ambient_dim()):
                m[j, i] = v[j]
        m.set_immutable()
        return m

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

        OUTPUT: a tuple of rays

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
            :class:`~sage.geometry.polyhedron.representation.Ray` objects.

        OUTPUT: list of rays as lists of coordinates

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

        OUTPUT: a tuple of lines

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
            iterate over the list of :class:`~sage.geometry.polyhedron.representation.Line` objects.

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.lines_list()
            [[1, 0]]
            sage: p.lines_list() == [list(x) for x in p.line_generator()]
            True
        """
        return [list(x) for x in self.line_generator()]

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

            sage: triangle = Polyhedron(vertices=[[1, 0], [0, 1], [1, 1]])
            sage: triangle.backend()
            'ppl'
            sage: D = polytopes.dodecahedron()                                          # needs sage.groups sage.rings.number_field
            sage: D.backend()                                                           # needs sage.groups sage.rings.number_field
            'field'
            sage: P = Polyhedron([[1.23]])
            sage: P.backend()
            'cdd'
        """
        return self.parent().backend()

    def cdd_Hrepresentation(self):
        r"""
        Write the inequalities/equations data of the polyhedron in
        cdd's H-representation format.

        .. SEEALSO::

            :meth:`write_cdd_Hrepresentation` -- export the polyhedron as a
            H-representation to a file.

        OUTPUT: string

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

            sage: triangle = Polyhedron(vertices=[[1,0], [0,1], [1,1]], base_ring=AA)   # needs sage.rings.number_field
            sage: triangle.base_ring()                                                  # needs sage.rings.number_field
            Algebraic Real Field
            sage: triangle.cdd_Hrepresentation()                                        # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be ZZ, QQ, or RDF
        """
        from .cdd_file_format import cdd_Hrepresentation
        try:
            cdd_type = self._cdd_type
        except AttributeError:
            from sage.rings.integer_ring import ZZ
            from sage.rings.rational_field import QQ
            if self.base_ring() is ZZ or self.base_ring() is QQ:
                cdd_type = 'rational'
            elif isinstance(self.base_ring(), sage.rings.abc.RealDoubleField):
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

        - ``filename`` -- the output file

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

        OUTPUT: string

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
            from sage.rings.integer_ring import ZZ
            from sage.rings.rational_field import QQ
            if self.base_ring() is ZZ or self.base_ring() is QQ:
                cdd_type = 'rational'
            elif isinstance(self.base_ring(), sage.rings.abc.RealDoubleField):
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

        - ``filename`` -- the output file

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
