r"""
Base class for polyhedra
"""


########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.structure.element import Element, coerce_binop, is_Vector

from sage.misc.all import cached_method, prod
from sage.misc.package import is_package_installed

from sage.rings.all import Integer, QQ, ZZ
from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.functions.other import sqrt, floor, ceil

from sage.graphs.graph import Graph

from sage.combinat.cartesian_product import CartesianProduct

from constructor import Polyhedron


#########################################################################
# Notes if you want to implement your own backend:
#
#  * derive from Polyhedron_base
#
#  * you must implement _init_from_Vrepresentation and
#    _init_from_Vrepresentationa
#
#  * You might want to override _init_empty_polyhedron,
#    _init_facet_adjacency_matrix, _init_vertex_adjacency_matrix, and
#    _make_polyhedron_face.
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
class Polyhedron_base(Element):
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

    Only one of ``Vrep`` or ``Hrep`` can be different from ``None``.

    TESTS::

        sage: p = Polyhedron()
        sage: TestSuite(p).run()
    """

    def __init__(self, parent, Vrep, Hrep, **kwds):
        """
        Initializes the polyhedron.

        See :class:`Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron()    # indirect doctests
        """
        Element.__init__(self, parent=parent)
        if Vrep is not None:
            vertices, rays, lines = Vrep
            self._init_from_Vrepresentation(vertices, rays, lines, **kwds)
        elif Hrep is not None:
            ieqs, eqns = Hrep
            self._init_from_Hrepresentation(ieqs, eqns, **kwds)
        else:
            self._init_empty_polyhedron()

    def _sage_input_(self, sib, coerced):
        """
        Return Sage command to reconstruct ``self``.

        See :mod:`sage.misc.sage_input` for details.

        EXAMPLES::

            sage: P = Polyhedron([(1,0), (0,1)], rays=[(1,1)])
            sage: sage_input(P)
            Polyhedron(base_ring=ZZ, rays=[(1, 1)], vertices=[(0, 1), (1, 0)])
       """
        kwds = dict()
        kwds['base_ring'] = sib(self.base_ring())
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
            NotImplementedError: A derived class must implement this method.
        """
        raise NotImplementedError('A derived class must implement this method.')

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
            NotImplementedError: A derived class must implement this method.
        """
        raise NotImplementedError('A derived class must implement this method.')

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
        self.parent()._make_Equation(self, [-1] + [0]*self.ambient_dim());
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
        """
        # TODO: This implementation computes the whole face lattice,
        # which is much more information than necessary.
        M = matrix(ZZ, self.n_Hrepresentation(), self.n_Hrepresentation(), 0)
        def set_adjacent(h1,h2):
            if h1 is h2:
                return
            i = h1.index()
            j = h2.index()
            M[i,j]=1
            M[j,i]=1

        face_lattice = self.face_lattice()
        for face in face_lattice:
            Hrep = face.ambient_Hrepresentation()
            if len(Hrep) == 2:
                set_adjacent(Hrep[0], Hrep[1])
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
        def set_adjacent(v1,v2):
            if v1 is v2:
                return
            i = v1.index()
            j = v2.index()
            M[i,j]=1
            M[j,i]=1

        face_lattice = self.face_lattice()
        for face in face_lattice:
            Vrep = face.ambient_Vrepresentation()
            if len(Vrep) == 2:
                set_adjacent(Vrep[0], Vrep[1])
        return M

    def delete(self):
        """
        Delete this polyhedron.

        This speeds up creation of new polyhedra by reusing
        objects. After recycling a polyhedron object, it is not in a
        consistent state any more and neither the polyhedron nor its
        H/V-representation objects may be used any more.

        .. seealso:: :meth:`~sage.geometry.polyhedron.parent.Polyhedra_base.recycle`

        EXAMPLES::

            sage: p = Polyhedron([(0,0),(1,0),(0,1)])
            sage: p.delete()

            sage: vertices = [(0,0,0,0),(1,0,0,0),(0,1,0,0),(1,1,0,0),(0,0,1,0),(0,0,0,1)]
            sage: def loop_polyhedra():
            ...       for i in range(0,100):
            ...           p = Polyhedron(vertices)

            sage: timeit('loop_polyhedra()')                   # not tested - random
            5 loops, best of 3: 79.5 ms per loop

            sage: def loop_polyhedra_with_recycling():
            ...       for i in range(0,100):
            ...           p = Polyhedron(vertices)
            ...           p.delete()

            sage: timeit('loop_polyhedra_with_recycling()')    # not tested - random
            5 loops, best of 3: 57.3 ms per loop
        """
        self.parent().recycle(self)

    def base_extend(self, base_ring, backend=None):
        """
        Return a new polyhedron over a larger field.

        INPUT:

        - ``base_ring`` -- the new base ring.

        - ``backend`` -- the new backend, see
          :func:`~sage.geometry.polyhedron.constructor.Polyhedron`.

        OUTPUT:

        The same polyhedron, but over a larger base ring.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)], base_ring=ZZ);  P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.base_extend(QQ)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 1 ray
            sage: P.base_extend(QQ) == P
            True
        """
        new_parent = self.parent().base_extend(base_ring, backend)
        return new_parent(self)

    def __cmp__(self, other):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        `-1, 0, +1` depending on how ``self`` and ``other``
        compare. If ``other`` is a polyhedron, then the comparison
        operator "less or equal than" means "is contained in", and
        "less than" means "is strictly contained in".

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: cmp(P,Q)
            1
            sage: cmp(Q,P)
            -1
            sage: cmp(P,P)
            0
            sage: abs(cmp(P, 'anything'))
            1

       The polytope ``Q`` is contained in ``P``::

            sage: P > Q
            True
            sage: P < Q
            False
            sage: P == Q
            False

        TESTS::

            sage: abs(cmp(P, 'string'))
            1
         """
        if not isinstance(other, Polyhedron_base):
            return -1
        if self._Vrepresentation is None or other._Vrepresentation is None:
            return -1   # make sure deleted polyhedra are not used in cache
        c = cmp(self.ambient_dim(), other.ambient_dim())
        if c != 0: return c
        c0 = self._is_subpolyhedron(other)
        c1 = other._is_subpolyhedron(self)
        if c0 and c1:
            return 0
        if c0:
            return -1
        else:
            return +1

    @coerce_binop
    def _is_subpolyhedron(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        sub-polyhdedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P._is_subpolyhedron(Q)
            False
            sage: Q._is_subpolyhedron(P)
            True
        """
        return all( other_H.contains(self_V)
                    for other_H, self_V in
                    CartesianProduct(other.Hrepresentation(), self.Vrepresentation()) )

    def plot(self,
             point=None, line=None, polygon=None, # None means unspecified by the user
             wireframe='blue', fill='green',
             projection_direction=None,
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

        - ``projection_direction`` -- coordinate list/tuple/iterable
          or ``None`` (default). The direction to use for the
          :meth:`schlegel_projection`` of the polytope. If not
          specified, no projection is used in dimensions `< 4` and
          parallel projection is used in dimension `4`.

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

        TESTS::

            sage: for p in square.plot():
            ...       print p.options()['rgbcolor'], p
            blue Point set defined by 4 point(s)
            blue Line defined by 2 points
            blue Line defined by 2 points
            blue Line defined by 2 points
            blue Line defined by 2 points
            green Polygon defined by 4 points

            sage: for p in line.plot():
            ...       print p.options()['rgbcolor'], p
            blue Point set defined by 2 point(s)
            green Line defined by 2 points

            sage: for p in point.plot():
            ...       print p.options()['rgbcolor'], p
            green Point set defined by 1 point(s)

        Draw the lines in red and nothing else::

            sage: for p in square.plot(point=False, line='red', polygon=False):
            ...       print p.options()['rgbcolor'], p
            red Line defined by 2 points
            red Line defined by 2 points
            red Line defined by 2 points
            red Line defined by 2 points

        Draw vertices in red, no lines, and a blue polygon::

            sage: for p in square.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):
            ...       print p.options()['rgbcolor'], p
            red Point set defined by 4 point(s)
            (0, 0, 1) Polygon defined by 4 points

            sage: for p in line.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):
            ...       print p.options()['rgbcolor'], p
            red Point set defined by 2 point(s)

            sage: for p in point.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):
            ...       print p.options()['rgbcolor'], p
            red Point set defined by 1 point(s)

        Draw in red without wireframe::

            sage: for p in square.plot(wireframe=False, fill="red"):
            ...       print p.options()['rgbcolor'], p
            red Polygon defined by 4 points

            sage: for p in line.plot(wireframe=False, fill="red"):
            ...       print p.options()['rgbcolor'], p
            red Line defined by 2 points

            sage: for p in point.plot(wireframe=False, fill="red"):
            ...       print p.options()['rgbcolor'], p
            red Point set defined by 1 point(s)

        The ``projection_direction`` option::

            sage: line3d = Polyhedron([(-1,-1,-1), (1,1,1)])
            sage: print(line3d.plot(projection_direction=[2,3,4]).description())
            Line defined by 2 points:           [(-0.00..., 0.126...), (0.131..., -1.93...)]
            Point set defined by 2 point(s):    [(-0.00..., 0.126...), (0.131..., -1.93...)]

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

        If the polyhedron is not full-dimensional, the :meth:`affine_hull` is used if necessary::

            sage: type(Polyhedron([(0,), (1,)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(Polyhedron([(0,0), (1,1)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(Polyhedron([(0,0,0), (1,1,1)]).plot())
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
            sage: type(Polyhedron([(0,0,0,0), (1,1,1,1)]).plot())
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
            sage: type(Polyhedron([(0,0,0,0,0), (1,1,1,1,1)]).plot())
            <class 'sage.plot.graphics.Graphics'>
        """
        def merge_options(*opts):
            merged = dict()
            for i in range(len(opts)):
                opt = opts[i]
                if opt is None:
                    continue
                elif opt is False:
                    return False
                elif isinstance(opt, (basestring, list, tuple)):
                    merged['color'] = opt
                else:
                    merged.update(opt)
            return merged

        d = min(self.dim(), 2)
        opts = [wireframe] * d + [fill] + [False] * (2-d)
        # The point/line/polygon options take precedence over wireframe/fill
        opts = [merge_options(opt1, opt2, kwds)
                for opt1, opt2 in zip(opts, [point, line, polygon])]

        def project(polyhedron):
            if projection_direction is not None:
                return polyhedron.schlegel_projection(projection_direction)
            elif polyhedron.ambient_dim() == 4:
                # There is no 4-d screen, we must project down to 3d
                return polyhedron.schlegel_projection()
            else:
                return polyhedron.projection()

        projection = project(self)
        try:
            plot_method = projection.plot
        except AttributeError:
            projection = project(self.affine_hull())
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
        if self.n_vertices()==0:
            desc += 'The empty polyhedron'
        else:
            desc += 'A ' + repr(self.dim()) + '-dimensional polyhedron'
        desc += ' in '
        desc += self.parent()._repr_ambient_module()

        if self.n_vertices()>0:
            desc += ' defined as the convex hull of '
            desc += repr(self.n_vertices())
            if self.n_vertices()==1: desc += ' vertex'
            else:                    desc += ' vertices'

            if self.n_rays()>0:
                if self.n_lines()>0: desc += ", "
                else:                desc += " and "
                desc += repr(self.n_rays())
                if self.n_rays()==1: desc += ' ray'
                else:                desc += ' rays'

            if self.n_lines()>0:
                if self.n_rays()>0: desc += ", "
                else:               desc += " and "
                desc += repr(self.n_lines())
                if self.n_lines()==1: desc +=' line'
                else:                 desc +=' lines'

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
        """
        Write the inequalities/equations data of the polyhedron in
        cdd's H-representation format.

        OUTPUT:

        A string. If you save the output to filename.ine then you can
        run the stand-alone cdd via ``cddr+ filename.ine``

        EXAMPLES::

            sage: p = polytopes.hypercube(2)
            sage: print p.cdd_Hrepresentation()
            H-representation
            begin
             4 3 rational
             1 1 0
             1 0 1
             1 -1 0
             1 0 -1
            end
        """
        from cdd_file_format import cdd_Hrepresentation
        try:
            cdd_type = self._cdd_type
        except AttributeError:
            if self.base_ring() is ZZ or self.base_ring() is QQ:
                cdd_type = 'rational'
            elif self.base_ring() is RDF:
                cdd_type = 'real'
            else:
                raise TypeError('The base ring must be ZZ, QQ, or RDF')
        return cdd_Hrepresentation(cdd_type,
                                   list(self.inequality_generator()),
                                   list(self.equation_generator()) )

    def cdd_Vrepresentation(self):
        """
        Write the vertices/rays/lines data of the polyhedron in cdd's
        V-representation format.

        OUTPUT:

        A string. If you save the output to filename.ext then you can
        run the stand-alone cdd via ``cddr+ filename.ext``

        EXAMPLES::

            sage: q = Polyhedron(vertices = [[1,1],[0,0],[1,0],[0,1]])
            sage: print q.cdd_Vrepresentation()
            V-representation
            begin
             4 3 rational
             1 0 0
             1 0 1
             1 1 0
             1 1 1
            end
        """
        from cdd_file_format import cdd_Vrepresentation
        try:
            cdd_type = self._cdd_type
        except AttributeError:
            if self.base_ring() is ZZ or self.base_ring() is QQ:
                cdd_type = 'rational'
            elif self.base_ring() is RDF:
                cdd_type = 'real'
            else:
                raise TypeError('The base ring must be ZZ, QQ, or RDF')
        return cdd_Vrepresentation(cdd_type,
                                   list(self.vertex_generator()),
                                   list(self.ray_generator()),
                                   list(self.line_generator()) )

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

    def Hrepresentation(self, index=None):
        """
        Return the objects of the H-representaton. Each entry is
        either an inequality or a equation.

        INPUT:

        - ``index`` -- either an integer or ``None``.

        OUTPUT:

        The optional argument is an index running from ``0`` to
        ``self.n_Hrepresentation()-1``. If present, the
        H-representation object at the given index will be
        returned. Without an argument, returns the list of all
        H-representation objects.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p.Hrepresentation(0)
            An inequality (0, 0, -1) x + 1 >= 0
            sage: p.Hrepresentation(0) == p.Hrepresentation() [0]
            True
        """
        if index is None:
            return self._Hrepresentation
        else:
            return self._Hrepresentation[index]


    def Hrep_generator(self):
        """
        Return an iterator over the objects of the H-representation
        (inequalities or equations).

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: next(p.Hrep_generator())
            An inequality (0, 0, -1) x + 1 >= 0
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

        - ``index`` -- either an integer or ``None``.

        OUTPUT:

        The optional argument is an index running from ``0`` to
        `self.n_Vrepresentation()-1``. If present, the
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
        Returns an iterator over the objects of the V-representation
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

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: triangle.vertices_list()
            [[0, 1], [1, 0], [1, 1]]
            sage: a_simplex = Polyhedron(ieqs = [
            ...            [0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]
            ...        ], eqns = [[1,-1,-1,-1,-1]])
            sage: a_simplex.vertices_list()
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
            sage: a_simplex.vertices_list() == [list(v) for v in a_simplex.vertex_generator()]
            True
        """
        return [list(x) for x in self.vertex_generator()]

    def vertex_generator(self):
        """
        Return a generator for the vertices of the polyhedron.

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
            ... except StopIteration: print "STOP"
            STOP
            sage: type(v_gen)
            <type 'generator'>
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

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: triangle.vertices()
            (A vertex at (0, 1), A vertex at (1, 0), A vertex at (1, 1))
            sage: a_simplex = Polyhedron(ieqs = [
            ...            [0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]
            ...        ], eqns = [[1,-1,-1,-1,-1]])
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
        """
        if base_ring is None:
            base_ring = self.base_ring()
        m = matrix(base_ring, self.ambient_dim(), self.n_vertices())
        for i,v in enumerate(self.vertices()):
            for j in range(0,self.ambient_dim()):
                m[j,i] = v[j]
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
            sage: for e in p.bounded_edges(): print e
            (A vertex at (0, 1), A vertex at (1, 0))
        """
        obj = self.Vrepresentation()
        for i in range(len(obj)):
            if not obj[i].is_vertex(): continue
            for j in range(i+1,len(obj)):
                if not obj[j].is_vertex(): continue
                if self.vertex_adjacency_matrix()[i,j] == 0: continue
                yield (obj[i], obj[j])

    def Vrepresentation_space(self):
        r"""
        Return the ambient vector space.

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

        The empty set is a special case (Trac #12193)::

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
            ...      print P.adjacency_matrix().row(v.index()), v
            (0, 1, 0, 1) A vertex at (0, 1)
            (1, 0, 1, 0) A vertex at (1, 0)
            (0, 1, 0, 1) A vertex at (3, 0)
            (1, 0, 1, 0) A vertex at (4, 1)

        If the V-representation of the polygon contains vertices and
        one ray, then each V-representation object is adjacent to two
        V-representation objects::

            sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (3, 0), (4, 1)],
            ...                  rays=[(0,1)])
            sage: for v in P.Vrep_generator():
            ...       print P.adjacency_matrix().row(v.index()), v
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
            ...                  rays=[(0,1), (1,1)])
            sage: for v in P.Vrep_generator():
            ...       print P.adjacency_matrix().row(v.index()), v
            (0, 1, 0, 0, 0) A ray in the direction (0, 1)
            (1, 0, 1, 0, 0) A vertex at (0, 1)
            (0, 1, 0, 0, 1) A vertex at (1, 0)
            (0, 0, 0, 0, 1) A ray in the direction (1, 1)
            (0, 0, 1, 1, 0) A vertex at (3, 0)
        """
        return self._vertex_adjacency_matrix()

    adjacency_matrix = vertex_adjacency_matrix

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
            :meth:`Vrepresentation`

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
        """
        incidence_matrix = matrix(ZZ, self.n_Vrepresentation(),
                                  self.n_Hrepresentation(), 0)
        for V in self.Vrep_generator():
            for H in self.Hrep_generator():
                if self._is_zero(H*V):
                    incidence_matrix[V.index(),H.index()] = 1
        return incidence_matrix

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        Either ``QQ`` (exact arithmetic using gmp, default) or ``RDF``
        (double precision floating-point arithmetic)

        EXAMPLES::

            sage: triangle = Polyhedron(vertices = [[1,0],[0,1],[1,1]])
            sage: triangle.base_ring() == ZZ
            True
        """
        return self.parent().base_ring()

    field = base_ring

    @cached_method
    def center(self):
        """
        Return the average of the vertices.

        See also :meth:`interior_point`.

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
        vertex_sum = vector(self.base_ring(), [0]*self.ambient_dim())
        for v in self.vertex_generator():
            vertex_sum += v.vector()
        vertex_sum.set_immutable()
        return vertex_sum / self.n_vertices()

    @cached_method
    def representative_point(self):
        """
        Return a "generic" point.

        See also :meth:`center`.

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


    @cached_method
    def radius_square(self):
        """
        Return the square of the maximal distance from the
        :meth:`center` to a vertex. All rays and lines are ignored.

        OUTPUT:

        The square of the radius, which is in :meth:`field`.

        EXAMPLES::

            sage: p = polytopes.permutahedron(4, project = False)
            sage: p.radius_square()
            5
        """
        vertices = [ v.vector() - self.center() for v in self.vertex_generator() ]
        return max( v.dot_product(v) for v in vertices )


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
        return self.n_rays()==0 and self.n_lines()==0

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
        if not self.is_compact(): return False

        for v in self.vertex_generator():
            adj = [a for a in v.neighbors()]
            if len(adj) != self.dim():
                return False
        return True

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

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.is_simplicial()
            Traceback (most recent call last):
            ...
            NotImplementedError: This function is implemented for polytopes only.
        """
        if not(self.is_compact()):
            raise NotImplementedError("This function is implemented for polytopes only.")
        d = self.dim()
        return all(len([vertex for vertex in face.incident()]) == d
                   for face in self.Hrepresentation())

    def hyperplane_arrangement(self):
        """
        Return the hyperplane arrangement defined by the equations and
        inequalities.

        OUTPUT:

        A :class:`hyperplane arrangement
        <sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement>`
        consisting of the hyperplanes defined by the
        :meth:`~sage.geometric.hyperplane_arragement.arrangement.HyperplaneArrangementElement.Hrepresentation`. 
        If the polytope is full-dimensional, this is the hyperplane
        arrangement spanned by the facets of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.hypercube(2)
            sage: p.hyperplane_arrangement()
            Arrangement <-t0 + 1 | -t1 + 1 | t1 + 1 | t0 + 1>
        """
        names = tuple('t'+str(i) for i in range(self.ambient_dim()))
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
            [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]

        REFERENCES:

            Lectures in Geometric Combinatorics, R.R.Thomas, 2006, AMS Press.
        """
        if not self.is_compact(): raise ValueError('Not a polytope.')

        A = matrix(self.n_vertices(),
                   [ [1]+x for x in self.vertex_generator()])
        A = A.transpose()
        A_ker = A.right_kernel()
        return A_ker.basis_matrix().transpose().rows()

    def triangulate(self, engine='auto', connected=True, fine=False, regular=None, star=None):
        r"""
        Returns a triangulation of the polytope.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal', or
          'TOPCOM'.  The latter two instruct this package to always
          use its own triangulation algorithms or TOPCOM's algorithms,
          respectively. By default ('auto'), TOPCOM is used if it is
          available and internal routines otherwise.

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
            ...      engine='internal') # to make doctest independent of TOPCOM
            sage: triangulation
            (<0,1,2,7>, <0,1,4,7>, <0,2,4,7>, <1,2,3,7>, <1,4,5,7>, <2,4,6,7>)
            sage: simplex_indices = triangulation[0]; simplex_indices
            (0, 1, 2, 7)
            sage: simplex_vertices = [ cube.Vrepresentation(i) for i in simplex_indices ]
            sage: simplex_vertices
            [A vertex at (-1, -1, -1), A vertex at (-1, -1, 1),
             A vertex at (-1, 1, -1), A vertex at (1, 1, 1)]
            sage: Polyhedron(simplex_vertices)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
        """
        if not self.is_compact():
            raise NotImplementedError('I can only triangulate compact polytopes.')
        from sage.geometry.triangulation.point_configuration import PointConfiguration
        pc = PointConfiguration((v.vector() for v in self.vertex_generator()),
                                connected=connected, fine=fine, regular=regular, star=star)
        pc.set_engine(engine)
        return pc.triangulate()

    @coerce_binop
    def Minkowski_sum(self, other):
        """
        Return the Minkowski sum.

        Minkowski addition of two subsets of a vector space is defined
        as

        .. math::

            X \oplus Y =
            \cup_{y\in Y} (X+y) =
            \cup_{x\in X, y\in Y} (x+y)

        See :meth:`Minkowski_difference` for a partial inverse operation.

        INPUT:

        - ``other`` -- a :class:`Polyhedron_base`.

        OUTPUT:

        The Minkowski sum of ``self`` and ``other``.

        EXAMPLES::

            sage: X = polytopes.hypercube(3)
            sage: Y = Polyhedron(vertices=[(0,0,0), (0,0,1/2), (0,1/2,0), (1/2,0,0)])
            sage: X+Y
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 13 vertices

            sage: four_cube = polytopes.hypercube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: four_cube + four_simplex
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 36 vertices
            sage: four_cube.Minkowski_sum(four_simplex) == four_cube + four_simplex
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

    _add_ = Minkowski_sum

    @coerce_binop
    def Minkowski_difference(self, other):
        """
        Return the Minkowski difference.

        Minkowski subtraction can equivalently be defined via
        Minkowski addition (see :meth:`Minkowski_sum`) or as
        set-theoretic intersection via

        .. math::

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

        - ``other`` -- a :class:`Polyhedron_base`.

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

        Minus sign is really an alias for :meth:`Minkowski_difference`
        ::

            sage: four_cube = polytopes.hypercube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: four_cube - four_simplex
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 16 vertices
            sage: four_cube.Minkowski_difference(four_simplex) == four_cube - four_simplex
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
        P = self.parent()
        return P.element_class(P, None, [new_ieqs, new_eqns])

    def __sub__(self, other):
        """
        Implement minus binary operation

        Polyhedra are not a ring with respect to dilatation and
        Minkowski sum, for example `X\oplus(-1)*Y \not= X\ominus Y`.

        INPUT:

        - ``other`` -- a translation vector or a polyhedron.

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
            return self.Minkowski_difference(other)
        return self + (-other)

    def is_Minkowski_summand(self, Y):
        """
        Test whether ``Y`` is a Minkowski summand.

        See :meth:`Minkowski_sum`.

        OUTPUT:

        Boolean. Whether there exists another polyhedron `Z` such that
        ``self`` can be written as `Y\oplus Z`.

        EXAMPLES::

            sage: A = polytopes.hypercube(2)
            sage: B = Polyhedron(vertices=[(0,1), (1/2,1)])
            sage: C = Polyhedron(vertices=[(1,1)])
            sage: A.is_Minkowski_summand(B)
            True
            sage: A.is_Minkowski_summand(C)
            True
            sage: B.is_Minkowski_summand(C)
            True
            sage: B.is_Minkowski_summand(A)
            False
            sage: C.is_Minkowski_summand(A)
            False
            sage: C.is_Minkowski_summand(B)
            False
        """
        return self.Minkowski_difference(Y).Minkowski_sum(Y) == self

    def translation(self, displacement):
        """
        Return the translated polyhedron.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector.

        OUTPUT:

        The translated polyhedron.

        EXAMPLES::

            sage: P = Polyhedron([[0,0],[1,0],[0,1]], base_ring=ZZ)
            sage: P.translation([2,1])
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: P.translation( vector(QQ,[2,1]) )
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
        """
        displacement = vector(displacement)
        new_vertices = [x.vector()+displacement for x in self.vertex_generator()]
        new_rays = self.rays()
        new_lines = self.lines()
        new_ring = self.parent()._coerce_base_ring(displacement)
        return Polyhedron(vertices=new_vertices, rays=new_rays, lines=new_lines, base_ring=new_ring)

    @coerce_binop
    def product(self, other):
        """
        Return the cartesian product.

        INPUT:

        - ``other`` -- a :class:`Polyhedron_base`.

        OUTPUT:

        The cartesian product of ``self`` and ``other`` with a
        suitable base ring to encompass the two.

        EXAMPLES::

            sage: P1 = Polyhedron([[0],[1]], base_ring=ZZ)
            sage: P2 = Polyhedron([[0],[1]], base_ring=QQ)
            sage: P1.product(P2)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        The cartesian product is the product in the semiring of polyhedra::

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
        """
        new_vertices = [ list(x)+list(y)
                         for x in self.vertex_generator() for y in other.vertex_generator()]
        new_rays = []
        new_rays.extend( [ r+[0]*other.ambient_dim()
                           for r in self.ray_generator() ] )
        new_rays.extend( [ [0]*self.ambient_dim()+r
                           for r in other.ray_generator() ] )
        new_lines = []
        new_lines.extend( [ l+[0]*other.ambient_dim()
                            for l in self.line_generator() ] )
        new_lines.extend( [ [0]*self.ambient_dim()+l
                            for l in other.line_generator() ] )
        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines,
                          base_ring=self.parent()._coerce_base_ring(other))

    _mul_ = product

    def dilation(self, scalar):
        """
        Return the dilated (uniformly stretched) polyhedron.

        INPUT:

        - ``scalar`` -- A scalar, not necessarily in :meth:`base_ring`.

        OUTPUT:

        The polyhedron dilated by that scalar, possibly coerced to a
        bigger field.

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

        Dilation of empty polyhedrons works, see :trac:`14987`::

            sage: p = Polyhedron(ambient_dim=2); p
            The empty polyhedron in ZZ^2
            sage: p.dilation(3)
            The empty polyhedron in ZZ^2

        TESTS::

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
        if scalar > 0:
            new_vertices = [ list(scalar*v.vector()) for v in self.vertex_generator() ]
            new_rays = self.rays()
            new_lines = self.lines()
        elif scalar < 0:
            new_vertices = [ list(scalar*v.vector()) for v in self.vertex_generator() ]
            new_rays = [ list(-r.vector()) for r in self.ray_generator()]
            new_lines = self.lines()
        else:
            new_vertices = [ self.ambient_space().zero() for v in self.vertex_generator() ]
            new_rays = []
            new_lines = []
        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines,
                          base_ring=self.parent()._coerce_base_ring(scalar),
                          ambient_dim=self.ambient_dim())

    def _acted_upon_(self, actor, self_on_left):
        """
        Implement the multiplicative action by scalars or other polyhedra.

        INPUT:

        - ``actor`` -- A scalar, not necessarily in :meth:`base_ring`,
          or a :class:`Polyhedron`.

        OUTPUT:

        Multiplication by another polyhedron returns the product
        polytope. Multiplication by a scalar returns the polytope
        dilated by that scalar, possibly coerced to the bigger field.

        EXAMPLES::

             sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,6)])
             sage: p._acted_upon_(2, True) == p.dilation(2)
             True
             sage: p*2 == p.dilation(2)
             True
             sage: p*p == p.product(p)
             True
             sage: p + vector(ZZ,[1,2,3]) == p.translation([1,2,3])
             True
        """
        if is_Polyhedron(actor):
            return self.product(actor)
        if is_Vector(actor):
            return self.translation(actor)
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

    def __div__(self, scalar):
        """
        Divide by a scalar factor.

        See :meth:`dilation` for details.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,4)])
            sage: (p/5).Vrepresentation()
            (A vertex at (2/5, 4/5, 8/5), A vertex at (3/5, 9/5, 27/5))
        """
        return self.dilation(1/scalar)

    @coerce_binop
    def convex_hull(self, other):
        """
        Return the convex hull of the set-theoretic union of the two
        polyhedra.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

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
        """
        Return the intersection of one polyhedron with another.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

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
        """
        new_ieqs = self.inequalities() + other.inequalities()
        new_eqns = self.equations() + other.equations()
        parent = self.parent()
        try:
            return parent.element_class(parent, None, [new_ieqs, new_eqns])
        except TypeError as msg:
            if self.base_ring() is ZZ:
                parent = parent.base_extend(QQ)
                return parent.element_class(parent, None, [new_ieqs, new_eqns])
            else:
                raise TypeError(msg)

    __and__ = intersection

    def edge_truncation(self, cut_frac = Integer(1)/3):
        r"""
        Return a new polyhedron formed from two points on each edge
        between two vertices.

        INPUT:

        - ``cut_frac`` -- integer. how deeply to cut into the edge.
            Default is `\frac{1}{3}`.

        OUTPUT:

        A Polyhedron object, truncated as described above.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: trunc_cube = cube.edge_truncation()
            sage: trunc_cube.n_vertices()
            24
            sage: trunc_cube.n_inequalities()
            14
        """
        new_vertices = []
        for e in self.bounded_edges():
            new_vertices.append((1-cut_frac)*e[0]() + cut_frac *e[1]())
            new_vertices.append(cut_frac *e[0]() + (1-cut_frac)*e[1]())

        new_vertices = [list(v) for v in new_vertices]
        new_rays =  self.rays()
        new_lines = self.lines()

        return Polyhedron(vertices=new_vertices, rays=new_rays,
                          lines=new_lines,
                          base_ring=self.parent()._coerce_base_ring(cut_frac))

    def _make_polyhedron_face(self, Vindices, Hindices):
        """
        Construct a face of the polyhedron.

        INPUT:

        - ``Vindices`` -- a tuple of integers. The indices of the
          V-represenation objects that span the face.

        - ``Hindices`` -- a tuple of integers. The indices of the
          H-representation objects that hold as equalities on the
          face.

        OUTPUT:

        A new :class:`~sage.geometry.polyhedron.face.PolyhedronFace` instance. It is not checked
        whether the input data actually defines a face.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: square._make_polyhedron_face((0,2), (1,))
            <0,2>
        """
        from sage.geometry.polyhedron.face import PolyhedronFace
        return PolyhedronFace(self, Vindices, Hindices)

    @cached_method
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

        For a full-dimensional polytope, the basic algorithm is
        described in
        :func:`~sage.geometry.hasse_diagram.Hasse_diagram_from_incidences`.
        There are three generalizations of [KP2002]_ necessary to deal
        with more general polytopes, corresponding to the extra
        H/V-representation objects:

        * Lines are removed before calling
          :func:`Hasse_diagram_from_incidences`, and then added back
          to each face V-representation except for the "empty face".

        * Equations are removed before calling
          :func:`Hasse_diagram_from_incidences`, and then added back
          to each face H-representation.

        * Rays: Consider the half line as an example. The
          V-representation objects are a point and a ray, which we can
          think of as a point at infinity. However, the point at
          infinity has no inequality associated to it, so there is
          only one H-representation object alltogether. The face
          lattice does not contain the "face at infinity". This means
          that in :func:`Hasse_diagram_from_incidences`, one needs to
          drop faces with V-representations that have no matching
          H-representation. In addition, one needs to ensure that
          every non-empty face contains at least one vertex.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: square.face_lattice()
            Finite poset containing 10 elements with distinguished linear extension
            sage: list(_)
            [<>, <0>, <1>, <2>, <3>, <0,1>, <0,2>, <2,3>, <1,3>, <0,1,2,3>]
            sage: poset_element = _[6]
            sage: a_face = poset_element
            sage: a_face
            <0,2>
            sage: a_face.dim()
            1
            sage: set(a_face.ambient_Vrepresentation()) == \
            ...   set([square.Vrepresentation(0), square.Vrepresentation(2)])
            True
            sage: a_face.ambient_Vrepresentation()
            (A vertex at (-1, -1), A vertex at (1, -1))
            sage: a_face.ambient_Hrepresentation()
            (An inequality (0, 1) x + 1 >= 0,)

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
            ...       for i in range(1,21)])
            sage: c5_20_fl = c5_20.face_lattice() # long time
            sage: [len(x) for x in c5_20_fl.level_sets()] # long time
            [1, 20, 190, 580, 680, 272, 1]
            sage: polytopes.hypercube(2).face_lattice().plot()
            Graphics object consisting of 27 graphics primitives
            sage: level_sets = polytopes.cross_polytope(2).face_lattice().level_sets()
            sage: print level_sets[0], level_sets[-1]
            [<>] [<0,1,2,3>]

        Various degenerate polyhedra::

            sage: Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0]]).face_lattice().level_sets()
            [[<>], [<0>, <1>, <2>], [<0,1>, <0,2>, <1,2>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(1,0,0),(0,1,0)], rays=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<1>, <2>], [<0,1>, <0,2>, <1,2>], [<0,1,2>]]
            sage: Polyhedron(rays=[(1,0,0),(0,1,0)], vertices=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0>], [<0,1>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(rays=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0>], [<0,1>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(1,),(0,)]).face_lattice().level_sets()
            [[<>], [<0>, <1>], [<0,1>]]
            sage: Polyhedron(vertices=[(1,0,0),(0,1,0)], lines=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0,1>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(lines=[(1,0,0)], vertices=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0,1>]]
            sage: Polyhedron(lines=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0,1,2>]]
            sage: Polyhedron(lines=[(1,0)], rays=[(0,1)], vertices=[(0,0)])\
            ...       .face_lattice().level_sets()
            [[<>], [<0,1>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(0,)], lines=[(1,)]).face_lattice().level_sets()
            [[<>], [<0,1>]]
            sage: Polyhedron(lines=[(1,0)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0,1>]]

        REFERENCES:

        ..  [KP2002]

            Volker Kaibel and Marc E. Pfetsch, "Computing the Face
            Lattice of a Polytope from its Vertex-Facet Incidences",
            Computational Geometry: Theory and Applications, Volume
            23, Issue 3 (November 2002), 281-290.  Available at
            http://portal.acm.org/citation.cfm?id=763203 and free of
            charge at http://arxiv.org/abs/math/0106043
        """
        coatom_to_Hindex = [ h.index() for h in self.inequality_generator() ]
        Hindex_to_coatom = [None] * self.n_Hrepresentation()
        for i in range(0,len(coatom_to_Hindex)):
            Hindex_to_coatom[ coatom_to_Hindex[i] ] = i

        atom_to_Vindex = [ v.index() for v in self.Vrep_generator() if not v.is_line() ]
        Vindex_to_atom = [None] * self.n_Vrepresentation()
        for i in range(0,len(atom_to_Vindex)):
                        Vindex_to_atom[ atom_to_Vindex[i] ] = i

        atoms_incidences   = [ tuple([ Hindex_to_coatom[h.index()]
                                       for h in v.incident() if h.is_inequality() ])
                               for v in self.Vrepresentation() if not v.is_line() ]

        coatoms_incidences = [ tuple([ Vindex_to_atom[v.index()]
                                       for v in h.incident() if not v.is_line() ])
                               for h in self.Hrepresentation() if h.is_inequality() ]

        atoms_vertices = [ Vindex_to_atom[v.index()] for v in self.vertex_generator() ]
        equations = [ e.index() for e in self.equation_generator() ]
        lines     = [ l.index() for l in self.line_generator() ]

        def face_constructor(atoms,coatoms):
            if len(atoms)==0:
                Vindices = ()
            else:
                Vindices = tuple(sorted([   atom_to_Vindex[i] for i in   atoms ]+lines))
            Hindices = tuple(sorted([ coatom_to_Hindex[i] for i in coatoms ]+equations))
            return self._make_polyhedron_face(Vindices, Hindices)

        from sage.geometry.hasse_diagram import Hasse_diagram_from_incidences
        return Hasse_diagram_from_incidences\
            (atoms_incidences, coatoms_incidences,
             face_constructor=face_constructor, required_atoms=atoms_vertices)

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
        random but fixed.

        EXAMPLES:

        Here we find the vertex and face indices of the eight three-dimensional
        facets of the four-dimensional hypercube::

            sage: p = polytopes.hypercube(4)
            sage: p.faces(3)
            (<0,1,2,3,4,5,6,7>, <0,1,2,3,8,9,10,11>, <0,1,4,5,8,9,12,13>,
             <0,2,4,6,8,10,12,14>, <2,3,6,7,10,11,14,15>, <8,9,10,11,12,13,14,15>,
             <4,5,6,7,12,13,14,15>, <1,3,5,7,9,11,13,15>)

            sage: face = p.faces(3)[0]
            sage: face.ambient_Hrepresentation()
            (An inequality (1, 0, 0, 0) x + 1 >= 0,)
            sage: face.vertices()
            (A vertex at (-1, -1, -1, -1), A vertex at (-1, -1, -1, 1),
             A vertex at (-1, -1, 1, -1), A vertex at (-1, -1, 1, 1),
             A vertex at (-1, 1, -1, -1), A vertex at (-1, 1, -1, 1),
             A vertex at (-1, 1, 1, -1), A vertex at (-1, 1, 1, 1))

        You can use the
        :meth:`~sage.geometry.polyhedron.representation.PolyhedronRepresentation.index`
        method to enumerate vertices and inequalities::

            sage: def get_idx(rep): return rep.index()
            sage: map(get_idx, face.ambient_Hrepresentation())
            [4]
            sage: map(get_idx, face.ambient_Vrepresentation())
            [0, 1, 2, 3, 4, 5, 6, 7]

            sage: [ (map(get_idx, face.ambient_Vrepresentation()), map(get_idx, face.ambient_Hrepresentation()))
            ...     for face in p.faces(3) ]
            [([0, 1, 2, 3, 4, 5, 6, 7], [4]),
             ([0, 1, 2, 3, 8, 9, 10, 11], [5]),
             ([0, 1, 4, 5, 8, 9, 12, 13], [6]),
             ([0, 2, 4, 6, 8, 10, 12, 14], [7]),
             ([2, 3, 6, 7, 10, 11, 14, 15], [2]),
             ([8, 9, 10, 11, 12, 13, 14, 15], [0]),
             ([4, 5, 6, 7, 12, 13, 14, 15], [1]),
             ([1, 3, 5, 7, 9, 11, 13, 15], [3])]

        TESTS::

            sage: pr = Polyhedron(rays = [[1,0,0],[-1,0,0],[0,1,0]], vertices = [[-1,-1,-1]], lines=[(0,0,1)])
            sage: pr.faces(4)
            ()
            sage: pr.faces(3)
            (<0,1,2,3>,)
            sage: pr.faces(2)
            (<0,1,2>,)
            sage: pr.faces(1)
            ()
            sage: pr.faces(0)
            ()
            sage: pr.faces(-1)
            ()
        """
        fl = self.face_lattice().level_sets()
        codim = self.dim() - face_dimension
        index = len(fl) - 1 - codim
        if index>=len(fl) or index<1:
            return tuple()
        return tuple(fl[index])

    @cached_method
    def f_vector(self):
        r"""
        Return the f-vector.

        OUTPUT:

        Returns a vector whose ``i``-th entry is the number of
        ``i``-dimensional faces of the polytope.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1, 2, 3], [1, 3, 2],
            ...       [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1], [0, 0, 0]])
            sage: p.f_vector()
            (1, 7, 12, 7, 1)
        """
        return vector(ZZ, [len(x) for x in self.face_lattice().level_sets()])

    def vertex_graph(self):
        """
        Return a graph in which the vertices correspond to vertices
        of the polyhedron, and edges to edges.

        EXAMPLES::

            sage: g3 = polytopes.hypercube(3).vertex_graph(); g3
            Graph on 8 vertices
            sage: g3.automorphism_group().cardinality()
            48
            sage: s4 = polytopes.simplex(4).vertex_graph(); s4
            Graph on 5 vertices
            sage: s4.is_eulerian()
            True
        """
        return Graph(self.vertex_adjacency_matrix(), loops=False)

    graph = vertex_graph

    def vertex_digraph(self, f, increasing=True):
        """
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
                raise TypeError('The linear map f must have '
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

    def polar(self):
        """
        Return the polar (dual) polytope.

        The original vertices are translated so that their barycenter
        is at the origin, and then the vertices are used as the
        coefficients in the polar inequalities.

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
        """
        assert self.is_compact(), "Not a polytope."

        verts = [list(self.center() - v.vector()) for v in self.vertex_generator()]
        base_ring = self.parent()._coerce_base_ring(self.center().parent())
        return Polyhedron(ieqs=[[1] + list(v) for v in verts], base_ring=base_ring)

    def pyramid(self):
        """
        Returns a polyhedron that is a pyramid over the original.

        EXAMPLES::

            sage: square = polytopes.hypercube(2);  square
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: egyptian_pyramid = square.pyramid();  egyptian_pyramid
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 5 vertices
            sage: egyptian_pyramid.n_vertices()
            5
            sage: for v in egyptian_pyramid.vertex_generator(): print v
            A vertex at (0, -1, -1)
            A vertex at (0, -1, 1)
            A vertex at (0, 1, -1)
            A vertex at (0, 1, 1)
            A vertex at (1, 0, 0)
        """
        new_verts = \
            [[0] + x for x in self.Vrep_generator()] + \
            [[1] + list(self.center())]
        return Polyhedron(vertices=new_verts)

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
        """
        new_verts = \
            [[ 0] + list(x) for x in self.vertex_generator()] + \
            [[ 1] + list(self.center())] + \
            [[-1] + list(self.center())]
        new_rays = [[0] + r for r in self.rays()]
        new_lines = [[0] + list(l) for l in self.lines()]
        return Polyhedron(vertices=new_verts, rays=new_rays, lines=new_lines)

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
        """
        new_verts = []
        new_verts.extend( [ [0] + v for v in self.vertices()] )
        new_verts.extend( [ [1] + v for v in self.vertices()] )
        new_rays =        [ [0] + r for r in self.rays()]
        new_lines =       [ [0] + l for l in self.lines()]
        return Polyhedron(vertices=new_verts, rays=new_rays, lines=new_lines,
                          base_ring=self.base_ring())

    def projection(self):
        """
        Return a projection object.

        See also
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.schlegel_projection`
        for a more interesting projection.

        OUTPUT:

        The identity projection. This is useful for plotting
        polyhedra.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: proj = p.projection()
            sage: proj
            The projection of a polyhedron into 3 dimensions
        """
        from plot import Projection
        self.projection = Projection(self)
        return self.projection

    def render_solid(self, **kwds):
        """
        Return a solid rendering of a 2- or 3-d polytope.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p_solid = p.render_solid(opacity = .7)
            sage: type(p_solid)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        proj = self.projection()
        if self.ambient_dim()==3:
            return proj.render_solid_3d(**kwds)
        if self.ambient_dim()==2:
            return proj.render_fill_2d(**kwds)
        raise ValueError("render_solid is only defined for 2 and 3 dimensional polyhedra.")

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
        if self.ambient_dim()==3:
            return proj.render_wireframe_3d(**kwds)
        if self.ambient_dim()==2:
            return proj.render_outline_2d(**kwds)
        raise ValueError("render_wireframe is only defined for 2 and 3 dimensional polyhedra.")

    def schlegel_projection(self, projection_dir=None, height=1.1):
        """
        Return the Schlegel projection.

        * The polyhedron is translated such that its
          :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.center`
          is at the origin.

        * The vertices are then normalized to the unit sphere

        * The normalized points are stereographically projected from a
          point slightly outside of the sphere.

        INPUT:

        - ``projection_direction`` -- coordinate list/tuple/iterable
          or ``None`` (default). The direction of the Schlegel
          projection. For a full-dimensional polyhedron, the default
          is the first facet normal; Otherwise, the vector consisting
          of the first n primes is chosen.

        - ``height`` -- float (default: `1.1`). How far outside of the
          unit sphere the focal point is.

        OUTPUT:

        A :class:`~sage.geometry.polyhedron.plot.Projection` object.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: sch_proj = p.schlegel_projection()
            sage: schlegel_edge_indices = sch_proj.lines
            sage: schlegel_edges = [sch_proj.coordinates_of(x) for x in schlegel_edge_indices]
            sage: len([x for x in schlegel_edges if x[0][0] > 0])
            4
        """
        proj = self.projection()
        if projection_dir is None:
            vertices = self.vertices()
            facet = self.Hrepresentation(0)
            f0 = [ v.index() for v in facet.incident() ]
            projection_dir = [sum([vertices[f0[i]][j]/len(f0) for i in range(len(f0))])
                              for j in range(self.ambient_dim())]
        return proj.schlegel(projection_direction=projection_dir, height=height)

    def _volume_lrs(self, verbose=False):
        """
        Computes the volume of a polytope using lrs.

        OUTPUT:

        The volume, cast to RDF (although lrs seems to output a
        rational value this must be an approximation in some cases).

        EXAMPLES::

            sage: polytopes.hypercube(3)._volume_lrs() #optional - lrs
            8.0
            sage: (polytopes.hypercube(3)*2)._volume_lrs() #optional - lrs
            64.0
            sage: polytopes.twenty_four_cell()._volume_lrs() #optional - lrs
            2.0

        REFERENCES:

             David Avis's lrs program.
        """
        if is_package_installed('lrs') != True:
            print 'You must install the optional lrs package ' \
                  'for this function to work'
            raise NotImplementedError

        from sage.misc.temporary_file import tmp_filename
        from subprocess import Popen, PIPE
        in_str = self.cdd_Vrepresentation()
        in_str += 'volume'
        in_filename = tmp_filename()
        in_file = open(in_filename, 'w')
        in_file.write(in_str)
        in_file.close()
        if verbose: print in_str

        lrs_procs = Popen(['lrs',in_filename],
                          stdin = PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = lrs_procs.communicate()
        if verbose:
            print ans
        # FIXME: check err

        for a_line in ans.splitlines():
            if 'Volume=' in a_line:
                volume = a_line.split('Volume=')[1]
                volume = RDF(QQ(volume))
                return volume

        raise ValueError("lrs did not return a volume")

    @cached_method
    def volume(self, engine='auto', **kwds):
        """
        Return the volume of the polytope.

        INPUT:

        - ``engine`` -- string. The backend to use. Allowed values are:

          * ``'auto'`` (default): see :meth:`triangulate`.
          * ``'internal'``: see :meth:`triangulate`.
          * ``'TOPCOM'``: see :meth:`triangulate`.
          * ``'lrs'``: use David Avis's lrs program (optional).

        - ``**kwds`` -- keyword arguments that are passed to the
          triangulation engine.

        OUTPUT:

        The volume of the polytope.

        EXAMPLES::

            sage: polytopes.hypercube(3).volume()
            8
            sage: (polytopes.hypercube(3)*2).volume()
            64
            sage: polytopes.twenty_four_cell().volume()
            2

        If the base ring is exact, the answer is exact::

            sage: P5 = polytopes.regular_polygon(5)
            sage: P5.volume()
            2.377641290737884?

            sage: polytopes.icosahedron().volume()
            5/12*sqrt5 + 5/4
            sage: numerical_approx(_)
            2.18169499062491

        Volume of the same polytopes, using the optional package lrs::

            sage: P5 = polytopes.regular_polygon(5, base_ring=RDF)
            sage: P5.volume(engine='lrs') #optional - lrs
            2.37764129...
        """
        if engine=='lrs':
            return self._volume_lrs(**kwds)
        dim = self.dim()
        if dim < self.ambient_dim():
            return self.base_ring().zero()
        triangulation = self.triangulate(engine=engine, **kwds)
        pc = triangulation.point_configuration()
        return sum([ pc.volume(simplex) for simplex in triangulation ]) / ZZ(dim).factorial()

    def contains(self, point):
        """
        Test whether the polyhedron contains the given ``point``.

        See also :meth:`interior_contains` and
        :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point (an iterable).

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

        The empty polyhedron needs extra care, see trac #10238::

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
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p)!=self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.contains(p):
                return False
        return True

    __contains__ = contains

    def interior_contains(self, point):
        """
        Test whether the interior of the polyhedron contains the
        given ``point``.

        See also :meth:`contains` and
        :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point.

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

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p)!=self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.interior_contains(p):
                return False
        return True

    def relative_interior_contains(self, point):
        """
        Test whether the relative interior of the polyhedron
        contains the given ``point``.

        See also :meth:`contains` and :meth:`interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point.

        OUTPUT:

        ``True`` or ``False``.

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

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.relative_interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p)!=self.ambient_dim():
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

        EXAMPLES::

            sage: Polyhedron([(0,0,0), (1,0,0), (0,1,0)]).is_simplex()
            True
            sage: polytopes.simplex(3).is_simplex()
            True
            sage: polytopes.hypercube(3).is_simplex()
            False
        """
        return self.is_compact() and (self.dim()+1==self.n_vertices())

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
            sage: lp.vertices_pc()
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
            sage: lp.vertices_pc()
            M(-1,  0),
            M( 0, -1),
            M( 0,  1),
            M( 1,  0),
            M( 1,  1)
            in 2-d lattice M
        """
        if not self.is_compact():
            raise NotImplementedError('Only compact lattice polytopes are allowed.')

        try:
            vertices = self.vertices_matrix(ZZ).columns()
        except TypeError:
            if not envelope:
                raise ValueError('Some vertices are not integral. '
                    'You probably want to add the argument '
                    '"envelope=True" to compute an enveloping lattice polytope.')
            vertices = []
            for v in self.vertex_generator():
                vbox = [ set([floor(x),ceil(x)]) for x in v ]
                vertices.extend( CartesianProduct(*vbox) )

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
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)]).lattice_polytope(True).points_pc()
            M( 0, -1),
            M(-1,  0),
            M(-1, -1),
            M( 0,  1),
            M( 1,  0),
            M( 1,  1),
            M( 0,  0)
            in 2-d lattice M
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)])._integral_points_PALP()
            [M(0, 1), M(1, 0), M(1, 1), M(0, 0)]
        """
        if not self.is_compact():
            raise ValueError('Can only enumerate points in a compact polyhedron.')
        lp = self.lattice_polytope(True)
        # remove cached values to get accurate timings
        try:
            del lp._points
            del lp._npoints
        except AttributeError:
            pass
        if self.is_lattice_polytope():
            return list(lp.points_pc())
        return [p for p in lp.points_pc() if self.contains(p)]

    @cached_method
    def bounding_box(self, integral=False):
        r"""
        Return the coordinates of a rectangular box containing the non-empty polytope.

        INPUT:

        - ``integral`` -- Boolean (default: ``False``). Whether to
          only allow integral coordinates in the bounding box.

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
            sage: polytopes.buckyball(exact=False).bounding_box()
            ((-0.8090169944, -0.8090169944, -0.8090169944), (0.8090169944, 0.8090169944, 0.8090169944))
        """
        box_min = []
        box_max = []
        if self.n_vertices==0:
            raise ValueError('Empty polytope is not allowed')
        if not self.is_compact():
            raise ValueError('Only polytopes (compact polyhedra) are allowed.')
        for i in range(0,self.ambient_dim()):
            coords = [ v[i] for v in self.vertex_generator() ]
            max_coord = max(coords)
            min_coord = min(coords)
            if integral:
                box_max.append(ceil(max_coord))
                box_min.append(floor(min_coord))
            else:
                box_max.append(max_coord)
                box_min.append(min_coord)
        return (tuple(box_min), tuple(box_max))

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

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ...        (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v)
            sage: pts1 = P.integral_points()                     # Sage's own code
            sage: all(P.contains(p) for p in pts1)
            True
            sage: pts2 = LatticePolytope(v).points_pc()          # PALP
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
        """
        if not self.is_compact():
            raise ValueError('Can only enumerate points in a compact polyhedron.')
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
        box_min, box_max = self.bounding_box(integral=True)
        box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))
        if  not self.is_lattice_polytope() or \
                (self.is_simplex() and box_points<1000) or \
                box_points<threshold:
            from sage.geometry.integral_points import rectangular_box_points
            return rectangular_box_points(box_min, box_max, self)

        # for more complicate polytopes, triangulate & use smith normal form
        from sage.geometry.integral_points import simplex_points
        if self.is_simplex():
            return simplex_points(self.Vrepresentation())
        triangulation = self.triangulate()
        points = set()
        for simplex in triangulation:
            triang_vertices = [ self.Vrepresentation(i) for i in simplex ]
            new_points = simplex_points(triang_vertices)
            for p in new_points:
                p.set_immutable()
            points.update(new_points)
        # assert all(self.contains(p) for p in points)   # slow
        return tuple(points)

    @cached_method
    def combinatorial_automorphism_group(self):
        """
        Computes the combinatorial automorphism group of the vertex
        graph of the polyhedron.

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the combinatorial automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while ``self.Vrepresentation()`` is indexed by
        nonnegative integers. The indexing of the permutation group is
        chosen to be shifted by ``+1``. That is, ``i`` in the
        permutation group corresponds to the V-representation object
        ``self.Vrepresentation(i-1)``.

        EXAMPLES::

            sage: quadrangle = Polyhedron(vertices=[(0,0),(1,0),(0,1),(2,3)])
            sage: quadrangle.combinatorial_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)(3,4)]
            sage: quadrangle.restricted_automorphism_group()
            Permutation Group with generators [()]

        Permutations can only exchange vertices with vertices, rays
        with rays, and lines with lines::

            sage: P = Polyhedron(vertices=[(1,0,0), (1,1,0)], rays=[(1,0,0)], lines=[(0,0,1)])
            sage: P.combinatorial_automorphism_group()
            Permutation Group with generators [(3,4)]
        """
        G = Graph()
        for edge in self.vertex_graph().edges():
            i = edge[0]
            j = edge[1]
            G.add_edge(i+1, j+1, (self.Vrepresentation(i).type(), self.Vrepresentation(j).type()) )

        group = G.automorphism_group(edge_labels=True)
        self._combinatorial_automorphism_group = group
        return group

    def _affine_coordinates(self, Vrep_object):
        r"""
        Return affine coordinates for a V-representation object.

        INPUT:

        - ``v`` -- a V-representation object or any iterable
          containing ``self.ambient_dim()`` coordinates. The
          coordinates must specify a point in the affine plane
          containing the polyhedron, or the output will be invalid. No
          checks on the input are performed.

        OUTPUT:

        A ``self.dim()``-dimensional coordinate vector. It contains
        the coordinates of ``v`` in an arbitrary but fixed basis for
        the affine span of the polyhedron.

        EXAMPLES::

            sage: P = Polyhedron(rays=[(1,0,0),(0,1,0)])
            sage: P._affine_coordinates( (-1,-2,-3) )
            (-1, -2)
            sage: [ P._affine_coordinates(v) for v in P.Vrep_generator() ]
            [(0, 0), (0, 1), (1, 0)]
        """
        if '_affine_coordinates_pivots' not in self.__dict__:
            v_list = [ vector(v) for v in self.Vrepresentation() ]
            if len(v_list)>0:
                origin = v_list[0]
                v_list = [ v - origin for v in v_list ]
            coordinates = matrix(v_list)
            self._affine_coordinates_pivots = coordinates.pivots()

        v = list(Vrep_object)
        if len(v) != self.ambient_dim():
            raise ValueError('Incorrect dimension: '+str(v))

        return vector(self.base_ring(), [ v[i] for i in self._affine_coordinates_pivots ])

    @cached_method
    def restricted_automorphism_group(self):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the Euclidean group `E(d) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional polyhedron. The Euclidean group
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
        the subgroup of rotations of the whole Euclidean group. The
        restricted automorphism group is

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

        OUTPUT:

        A :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the restricted automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while ``self.Vrepresentation()`` is indexed by
        nonnegative integers. The indexing of the permutation group is
        chosen to be shifted by ``+1``. That is, ``i`` in the
        permutation group corresponds to the V-representation object
        ``self.Vrepresentation(i-1)``.

        REFERENCES:

        ..  [BSS]
            David Bremner, Mathieu Dutour Sikiric, Achill Schuermann:
            Polyhedral representation conversion up to symmetries.
            http://arxiv.org/abs/math/0702239

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3)
            sage: AutP = P.restricted_automorphism_group();  AutP
            Permutation Group with generators [(3,4), (2,3)(4,5), (2,5), (1,2)(5,6), (1,6)]
            sage: P24 = polytopes.twenty_four_cell()
            sage: AutP24 = P24.restricted_automorphism_group()
            sage: PermutationGroup([
            ...     '(3,6)(4,7)(10,11)(14,15)(18,21)(19,22)',
            ...     '(2,3)(7,8)(11,12)(13,14)(17,18)(22,23)',
            ...     '(2,5)(3,10)(6,11)(8,17)(9,13)(12,16)(14,19)(15,22)(20,23)',
            ...     '(2,10)(3,5)(6,12)(7,18)(9,14)(11,16)(13,19)(15,23)(20,22)',
            ...     '(2,11)(3,12)(4,21)(5,6)(9,15)(10,16)(13,22)(14,23)(19,20)',
            ...     '(1,2)(3,4)(6,7)(8,9)(12,13)(16,17)(18,19)(21,22)(23,24)',
            ...     '(1,24)(2,13)(3,14)(5,9)(6,15)(10,19)(11,22)(12,23)(16,20)'
            ...   ]) == AutP24
            True

        Here is the quadrant example mentioned in the beginning::

            sage: P = Polyhedron(rays=[(1,0),(0,1)])
            sage: P.Vrepresentation()
            (A vertex at (0, 0), A ray in the direction (0, 1), A ray in the direction (1, 0))
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(2,3)]

        Also, the polyhedron need not be full-dimensional::

            sage: P = Polyhedron(vertices=[(1,2,3,4,5),(7,8,9,10,11)])
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(1,2)]

        Translations do not change the restricted automorphism
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

        Floating-point computations are supported with a simple fuzzy
        zero implementation::

            sage: P = Polyhedron(vertices=[(1.0/3.0,0,0),(0,1.0/3.0,0),(0,0,1.0/3.0)], base_ring=RDF)
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]

        TESTS::

            sage: p = Polyhedron(vertices=[(1,0), (1,1)], rays=[(1,0)])
            sage: p.restricted_automorphism_group()
            Permutation Group with generators [(2,3)]
        """
        if self.base_ring() is ZZ or self.base_ring() is QQ:
            def rational_approximation(c):
                return c

        elif self.base_ring() is RDF:
            c_list = []
            def rational_approximation(c):
                # Implementation detail: Return unique integer if two
                # c-values are the same up to machine precision. But
                # you can think of it as a uniquely-chosen rational
                # approximation.
                for i,x in enumerate(c_list):
                    if self._is_zero(x-c):
                        return i
                c_list.append(c)
                return len(c_list)-1

        # The algorithm identifies the restricted automorphism group
        # with the automorphism group of a edge-colored graph. The
        # nodes of the graph are the V-representation objects. If all
        # V-representation objects are vertices, the edges are
        # labelled by numbers (to be computed below). Roughly
        # speaking, the edge label is the inner product of the
        # coordinate vectors with some orthogonalization thrown in
        # [BSS].
        def edge_label_compact(i,j,c_ij):
            return c_ij

        # In the non-compact case we also label the edges by the type
        # of the V-representation object. This ensures that vertices,
        # rays, and lines are only permuted amongst themselves.
        def edge_label_noncompact(i,j,c_ij):
            return (self.Vrepresentation(i).type(), c_ij, self.Vrepresentation(j).type())

        if self.is_compact():
            edge_label = edge_label_compact
        else:
            edge_label = edge_label_noncompact

        # good coordinates for the V-representation objects
        v_list = []
        for v in self.Vrepresentation():
            v_coords = list(self._affine_coordinates(v))
            if v.is_vertex():
                v_coords = [1]+v_coords
            else:
                v_coords = [0]+v_coords
            v_list.append(vector(v_coords))

        # Finally, construct the graph
        Qinv = sum( v.column() * v.row() for v in v_list ).inverse()
        G = Graph()
        for i in range(0,len(v_list)):
            for j in range(i+1,len(v_list)):
                v_i = v_list[i]
                v_j = v_list[j]
                c_ij = rational_approximation( v_i * Qinv * v_j )
                G.add_edge(i+1,j+1, edge_label(i,j,c_ij))

        group = G.automorphism_group(edge_labels=True)
        self._restricted_automorphism_group = group
        return group

    def is_full_dimensional(self):
        """
        Return whether the polyhedron is full dimensional.

        OUTPUT:

        Boolean. Whether the polyhedron is not contained in any strict
        affine subspace.

        EXAMPLES::

            sage: polytopes.hypercube(3).is_full_dimensional()
            True
            sage: Polyhedron(vertices=[(1,2,3)], rays=[(1,0,0)]).is_full_dimensional()
            False
        """
        return self.dim() == self.ambient_dim()

    def affine_hull(self):
        """
        Return the affine hull.

        Each polyhedron is contained in some smallest affine subspace
        (possibly the entire ambient space). The affine hull is the
        same polyhedron but thought of as a full-dimensional
        polyhedron in this subspace.

        OUTPUT:

        A full-dimensional polyhedron.

        EXAMPLES::

            sage: triangle = Polyhedron([(1,0,0), (0,1,0), (0,0,1)]);  triangle
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: triangle.affine_hull()
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices

            sage: half3d = Polyhedron(vertices=[(3,2,1)], rays=[(1,0,0)])
            sage: half3d.affine_hull().Vrepresentation()
            (A ray in the direction (1), A vertex at (3))

        TESTS::

            sage: Polyhedron([(2,3,4)]).affine_hull()
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
        """
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
        pivots = matrix(gens, base_ring=self.base_ring()).pivots()
        def pivot(indexed):
            return [indexed[i] for i in pivots]

        vertices = map(pivot, self.vertices())
        rays = map(pivot, self.rays())
        lines = map(pivot, self.lines())
        return Polyhedron(vertices=vertices, rays=rays, lines=lines, base_ring=self.base_ring())
