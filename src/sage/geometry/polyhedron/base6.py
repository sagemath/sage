r"""
Base class for polyhedra: Methods for plotting and affine hull projection
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
from sage.modules.vector_space_morphism import linear_transformation
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.geometry.convex_set import AffineHullProjectionData
from .base5 import Polyhedron_base5


class Polyhedron_base6(Polyhedron_base5):
    r"""
    Methods related to plotting including affine hull projection.

    TESTS::

        sage: from sage.geometry.polyhedron.base6 import Polyhedron_base6
        sage: P = polytopes.cube()
        sage: Polyhedron_base6.plot(P)                                                  # needs sage.plot
        Graphics3d Object
        sage: print(Polyhedron_base6.tikz(P, output_type='TikzPicture'))                # needs sage.plot
        \RequirePackage{luatex85}
        \documentclass[tikz]{standalone}
        \begin{document}
        \begin{tikzpicture}%
            [x={(1.000000cm, 0.000000cm)},
            y={(-0.000000cm, 1.000000cm)},
            z={(0.000000cm, -0.000000cm)},
            scale=1.000000,
            back/.style={loosely dotted, thin},
            edge/.style={color=blue!95!black, thick},
            facet/.style={fill=blue!95!black,fill opacity=0.800000},
            vertex/.style={inner sep=1pt,circle,draw=green!25!black,fill=green!75!black,thick}]
        %
        %
        %% This TikZ-picture was produced with Sagemath version ...
        %% with the command: ._tikz_3d_in_3d and parameters:
        %% view = [0, 0, 1]
        %% angle = 0
        %% scale = 1
        %% edge_color = blue!95!black
        %% facet_color = blue!95!black
        %% opacity = 0.8
        %% vertex_color = green
        %% axis = False
        %%
        %% Coordinate of the vertices:
        %%
        \coordinate (1.00000, -1.00000, -1.00000) at (1.00000, -1.00000, -1.00000);
        \coordinate (1.00000, 1.00000, -1.00000) at (1.00000, 1.00000, -1.00000);
        \coordinate (1.00000, 1.00000, 1.00000) at (1.00000, 1.00000, 1.00000);
        \coordinate (1.00000, -1.00000, 1.00000) at (1.00000, -1.00000, 1.00000);
        \coordinate (-1.00000, -1.00000, 1.00000) at (-1.00000, -1.00000, 1.00000);
        \coordinate (-1.00000, -1.00000, -1.00000) at (-1.00000, -1.00000, -1.00000);
        \coordinate (-1.00000, 1.00000, -1.00000) at (-1.00000, 1.00000, -1.00000);
        \coordinate (-1.00000, 1.00000, 1.00000) at (-1.00000, 1.00000, 1.00000);
        %%
        %%
        %% Drawing edges in the back
        %%
        \draw[edge,back] (1.00000, -1.00000, -1.00000) -- (1.00000, 1.00000, -1.00000);
        \draw[edge,back] (1.00000, -1.00000, -1.00000) -- (1.00000, -1.00000, 1.00000);
        \draw[edge,back] (1.00000, -1.00000, -1.00000) -- (-1.00000, -1.00000, -1.00000);
        \draw[edge,back] (1.00000, 1.00000, -1.00000) -- (1.00000, 1.00000, 1.00000);
        \draw[edge,back] (1.00000, 1.00000, -1.00000) -- (-1.00000, 1.00000, -1.00000);
        \draw[edge,back] (-1.00000, -1.00000, 1.00000) -- (-1.00000, -1.00000, -1.00000);
        \draw[edge,back] (-1.00000, -1.00000, -1.00000) -- (-1.00000, 1.00000, -1.00000);
        \draw[edge,back] (-1.00000, 1.00000, -1.00000) -- (-1.00000, 1.00000, 1.00000);
        %%
        %%
        %% Drawing vertices in the back
        %%
        \node[vertex] at (1.00000, -1.00000, -1.00000)     {};
        \node[vertex] at (1.00000, 1.00000, -1.00000)     {};
        \node[vertex] at (-1.00000, 1.00000, -1.00000)     {};
        \node[vertex] at (-1.00000, -1.00000, -1.00000)     {};
        %%
        %%
        %% Drawing the facets
        %%
        \fill[facet] (-1.00000, 1.00000, 1.00000) -- (1.00000, 1.00000, 1.00000) -- (1.00000, -1.00000, 1.00000) -- (-1.00000, -1.00000, 1.00000) -- cycle {};
        %%
        %%
        %% Drawing edges in the front
        %%
        \draw[edge] (1.00000, 1.00000, 1.00000) -- (1.00000, -1.00000, 1.00000);
        \draw[edge] (1.00000, 1.00000, 1.00000) -- (-1.00000, 1.00000, 1.00000);
        \draw[edge] (1.00000, -1.00000, 1.00000) -- (-1.00000, -1.00000, 1.00000);
        \draw[edge] (-1.00000, -1.00000, 1.00000) -- (-1.00000, 1.00000, 1.00000);
        %%
        %%
        %% Drawing the vertices in the front
        %%
        \node[vertex] at (1.00000, 1.00000, 1.00000)     {};
        \node[vertex] at (1.00000, -1.00000, 1.00000)     {};
        \node[vertex] at (-1.00000, -1.00000, 1.00000)     {};
        \node[vertex] at (-1.00000, 1.00000, 1.00000)     {};
        %%
        %%
        \end{tikzpicture}
        \end{document}

        sage: Q = polytopes.hypercube(4)
        sage: Polyhedron_base6.show(Q)                                                  # needs sage.plot
        sage: Polyhedron_base6.schlegel_projection(Q)
        The projection of a polyhedron into 3 dimensions

        sage: R = polytopes.simplex(5)
        sage: Polyhedron_base6.affine_hull(R)
        A 5-dimensional polyhedron in ZZ^6 defined as the convex hull of 1 vertex and 5 lines
        sage: Polyhedron_base6.affine_hull_projection(R)
        A 5-dimensional polyhedron in ZZ^5 defined as the convex hull of 6 vertices
    """
    def plot(self,
             point=None, line=None, polygon=None,  # None means unspecified by the user
             wireframe='blue', fill='green',
             position=None,
             orthonormal=True,  # whether to use orthonormal projections
             **kwds):
        r"""
        Return a graphical representation.

        INPUT:

        - ``point``, ``line``, ``polygon`` -- parameters to pass to
          point (0d), line (1d), and polygon (2d) plot commands.
          Allowed values are:

          * A Python dictionary to be passed as keywords to the plot
            commands.

          * A string or triple of numbers: The color. This is
            equivalent to passing the dictionary ``{'color':...}``.

          * ``False``: Switches off the drawing of the corresponding
            graphics object

        - ``wireframe``, ``fill`` -- similar to ``point``, ``line``,
          and ``polygon``, but ``fill`` is used for the graphics
          objects in the dimension of the polytope (or of dimension 2
          for higher dimensional polytopes) and ``wireframe`` is used
          for all lower-dimensional graphics objects
          (default: 'green' for ``fill`` and 'blue' for ``wireframe``)

        - ``position`` -- positive number; the position to take the projection
          point in Schlegel diagrams

        - ``orthonormal`` -- boolean (default: ``True``); whether to use
          orthonormal projections

        - ``**kwds`` -- optional keyword parameters that are passed to
          all graphics objects

        OUTPUT:

        A (multipart) graphics object.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: point = Polyhedron([[1,1]])
            sage: line = Polyhedron([[1,1],[2,1]])
            sage: cube = polytopes.hypercube(3)
            sage: hypercube = polytopes.hypercube(4)

        By default, the wireframe is rendered in blue and the fill in green::

            sage: # needs sage.plot
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

            sage: # needs sage.plot
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

            sage: # needs sage.plot
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

            sage: # needs sage.plot
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
            sage: facet = hcube.facets()[0].as_polyhedron(); facet
            A 4-dimensional polyhedron in ZZ^5 defined as the convex hull of 16 vertices
            sage: facet.plot()                                                          # needs sage.plot
            Graphics3d Object

        For a 3d plot, we may draw the polygons with rainbow colors, using any of the following ways::

            sage: cube.plot(polygon='rainbow')                                          # needs sage.plot
            Graphics3d Object
            sage: cube.plot(polygon={'color':'rainbow'})                                # needs sage.plot
            Graphics3d Object
            sage: cube.plot(fill='rainbow')                                             # needs sage.plot
            Graphics3d Object

        For a 3d plot, the size of a point, the thickness of a line and the width of an arrow
        are controlled by the respective parameters::

            sage: prism = Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0]], rays=[[0,0,1]])
            sage: prism.plot(size=20, thickness=30, width=1)                            # needs sage.plot
            Graphics3d Object
            sage: prism.plot(point={'size':20, 'color':'black'},                        # needs sage.plot
            ....:            line={'thickness':30, 'width':1, 'color':'black'},
            ....:            polygon='rainbow')
            Graphics3d Object

        TESTS::

            sage: for p in square.plot():                                               # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            blue Point set defined by 4 point(s)
            blue Line defined by 2 points
            blue Line defined by 2 points
            blue Line defined by 2 points
            blue Line defined by 2 points
            green Polygon defined by 4 points

            sage: for p in line.plot():                                                 # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            blue Point set defined by 2 point(s)
            green Line defined by 2 points

            sage: for p in point.plot():                                                # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            green Point set defined by 1 point(s)

        Draw the lines in red and nothing else::

            sage: for p in square.plot(point=False, line='red', polygon=False):         # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Line defined by 2 points
            red Line defined by 2 points
            red Line defined by 2 points
            red Line defined by 2 points

        Draw vertices in red, no lines, and a blue polygon::

            sage: for p in square.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):         # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 4 point(s)
            (0, 0, 1) Polygon defined by 4 points

            sage: for p in line.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):           # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 2 point(s)

            sage: for p in point.plot(point={'color':'red'}, line=False, polygon=(0,0,1)):          # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 1 point(s)

        Draw in red without wireframe::

            sage: for p in square.plot(wireframe=False, fill='red'):                    # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Polygon defined by 4 points

            sage: for p in line.plot(wireframe=False, fill='red'):                      # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Line defined by 2 points

            sage: for p in point.plot(wireframe=False, fill='red'):                     # needs sage.plot
            ....:     print("{} {}".format(p.options()['rgbcolor'], p))
            red Point set defined by 1 point(s)

        We try to draw the polytope in 2 or 3 dimensions::

            sage: # needs sage.plot
            sage: type(Polyhedron(ieqs=[(1,)]).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(polytopes.hypercube(1).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(polytopes.hypercube(2).plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: type(polytopes.hypercube(3).plot())
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>

        In 4d a projection to 3d is used::

            sage: type(polytopes.hypercube(4).plot())                                   # needs sage.plot
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
            sage: type(polytopes.hypercube(5).plot())
            Traceback (most recent call last):
            ...
            NotImplementedError: plotting of 5-dimensional polyhedra not implemented

        If the polyhedron is not full-dimensional, the :meth:`affine_hull_projection` is used if necessary::

            sage: # needs sage.plot
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

        Check that :issue:`30015` is fixed::

            sage: fcube = polytopes.hypercube(4)
            sage: tfcube = fcube.face_truncation(fcube.faces(0)[0])
            sage: sp = tfcube.schlegel_projection()
            sage: for face in tfcube.faces(2):
            ....:     vertices = face.ambient_Vrepresentation()
            ....:     indices = [sp.coord_index_of(vector(x)) for x in vertices]
            ....:     projected_vertices = [sp.transformed_coords[i] for i in indices]
            ....:     assert Polyhedron(projected_vertices).dim() == 2

        Check that :issue:`31802` is fixed::

            sage: # needs sage.plot
            sage: halfspace = Polyhedron(rays=[(0, 0, 1)], lines=[(1, 0, 0), (0, 1, 0)])
            sage: len(halfspace.projection().arrows)
            5
            sage: halfspace.plot(fill=(0, 1, 0))
            Graphics3d Object
            sage: fullspace = Polyhedron(lines=[(1, 0, 0), (0, 1, 0), (0, 0, 1)])
            sage: len(fullspace.projection().arrows)
            6
            sage: fullspace.plot(color=(1, 0, 0), alpha=0.5)
            Graphics3d Object
            sage: cone = Polyhedron(rays=[(1, 0, 0), (0, 1, 0), (0, 0, 1)])
            sage: cone.plot(fill='rainbow', alpha=0.6)
            Graphics3d Object
            sage: p = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0)], rays=[(-1, 1, 0), (1, 1, 0), (0, 0, 1)])
            sage: p.plot(fill='mediumspringgreen', point='red', size=30, width=2)
            Graphics3d Object

            sage: cylinder = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0), (0, 1, 0)], lines=[(0, 0, 1)])
            sage: cylinder.plot(fill='red')  # check it is not all black                # needs sage.plot
            Graphics3d Object
            sage: quarter = Polyhedron(rays=[(-1, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
            sage: quarter.plot(fill='rainbow')  # check it is not all black nor with too many colors                    # needs sage.plot
            Graphics3d Object
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
        r"""
        Display graphics immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        INPUT:

        - ``kwds`` -- optional keyword arguments; see :meth:`plot` for
          the description of available options

        OUTPUT:

        This method does not return anything. Use :meth:`plot` if you
        want to generate a graphics object that can be saved or
        further transformed.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: square.show(point='red')                                              # needs sage.plot
        """
        self.plot(**kwds).show()

    def tikz(self, view=[0, 0, 1], angle=0, scale=1,
             edge_color='blue!95!black', facet_color='blue!95!black',
             opacity=0.8, vertex_color='green', axis=False,
             output_type=None):
        r"""
        Return a tikz picture of ``self`` as a string or as a
        :class:`~sage.misc.latex_standalone.TikzPicture`
        according to a projection ``view`` and an angle ``angle``
        obtained via the threejs viewer. ``self`` must be bounded.

        INPUT:

        - ``view`` -- list (default: [0,0,1]) representing the rotation axis (see note below)
        - ``angle`` -- integer (default: 0); angle of rotation in degree from 0 to 360 (see note
          below)
        - ``scale`` -- integer (default: 1); the scaling of the tikz picture
        - ``edge_color`` -- string (default: ``'blue!95!black'``); representing colors which tikz
          recognizes
        - ``facet_color`` -- string (default: ``'blue!95!black'``); representing colors which tikz
          recognizes
        - ``vertex_color`` -- string (default: ``'green'``); representing colors which tikz
          recognizes
        - ``opacity`` -- real number (default: 0.8) between 0 and 1 giving the opacity of
          the front facets
        - ``axis`` -- boolean (default: ``False``); draw the axes at the origin or not
        - ``output_type`` -- string (default: ``None``); valid values
          are ``None`` (deprecated), ``'LatexExpr'`` and ``'TikzPicture'``,
          whether to return a LatexExpr object (which inherits from Python
          str) or a ``TikzPicture`` object from module
          :mod:`sage.misc.latex_standalone`

        OUTPUT: LatexExpr object or TikzPicture object

        .. NOTE::

            This is a wrapper of a method of the projection object
            ``self.projection()``. See :meth:`~sage.geometry.polyhedron.plot.Projection.tikz`
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

            sage: # needs sage.plot
            sage: co = polytopes.cuboctahedron()
            sage: Img = co.tikz([0, 0, 1], 0, output_type='TikzPicture')
            sage: Img
            \documentclass[tikz]{standalone}
            \begin{document}
            \begin{tikzpicture}%
                    [x={(1.000000cm, 0.000000cm)},
                    y={(0.000000cm, 1.000000cm)},
                    z={(0.000000cm, 0.000000cm)},
                    scale=1.000000,
            ...
            Use print to see the full content.
            ...
            \node[vertex] at (1.00000, 0.00000, 1.00000)     {};
            \node[vertex] at (1.00000, 1.00000, 0.00000)     {};
            %%
            %%
            \end{tikzpicture}
            \end{document}
            sage: print('\n'.join(Img.content().splitlines()[12:21]))
            %% with the command: ._tikz_3d_in_3d and parameters:
            %% view = [0, 0, 1]
            %% angle = 0
            %% scale = 1
            %% edge_color = blue!95!black
            %% facet_color = blue!95!black
            %% opacity = 0.8
            %% vertex_color = green
            %% axis = False
            sage: print('\n'.join(Img.content().splitlines()[22:26]))
            %% Coordinate of the vertices:
            %%
            \coordinate (-1.00000, -1.00000, 0.00000) at (-1.00000, -1.00000, 0.00000);
            \coordinate (-1.00000, 0.00000, -1.00000) at (-1.00000, 0.00000, -1.00000);

        When output type is a :class:`sage.misc.latex_standalone.TikzPicture`::

            sage: # needs sage.plot
            sage: co = polytopes.cuboctahedron()
            sage: t = co.tikz([674, 108, -731], 112, output_type='TikzPicture'); t
            \documentclass[tikz]{standalone}
            \begin{document}
            \begin{tikzpicture}%
                    [x={(0.249656cm, -0.577639cm)},
                    y={(0.777700cm, -0.358578cm)},
                    z={(-0.576936cm, -0.733318cm)},
                    scale=1.000000,
            ...
            Use print to see the full content.
            ...
            \node[vertex] at (1.00000, 0.00000, 1.00000)     {};
            \node[vertex] at (1.00000, 1.00000, 0.00000)     {};
            %%
            %%
            \end{tikzpicture}
            \end{document}
            sage: path_to_file = t.pdf()        # not tested
        """
        return self.projection().tikz(view, angle, scale,
                                      edge_color, facet_color,
                                      opacity, vertex_color, axis,
                                      output_type=output_type)

    def _rich_repr_(self, display_manager, **kwds):
        r"""
        Rich Output Magic Method.

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

    @cached_method
    def gale_transform(self):
        r"""
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

            :func:`~sage.geometry.polyhedron.library.gale_transform_to_polyhedron`.

        TESTS::

            sage: P = Polyhedron(rays=[[1,0,0]])
            sage: P.gale_transform()
            Traceback (most recent call last):
            ...
            ValueError: not a polytope

        Check that :issue:`29073` is fixed::

            sage: P = polytopes.icosahedron(exact=False)                                # needs sage.groups
            sage: sum(P.gale_transform()).norm() < 1e-15                                # needs sage.groups
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
        r"""
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

        # Check :issue:`29073`.
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

            try:
                import sage.graphs.graph
                assert sage.graphs.graph  # to muffle pyflakes
            except ImportError:
                pass
            else:
                tester.assertTrue(self.is_combinatorially_isomorphic(P))

    def projection(self, projection=None):
        r"""
        Return a projection object.

        INPUT:

        - ``proj`` -- a projection function

        OUTPUT:

        The identity projection. This is useful for plotting
        polyhedra.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base6.Polyhedron_base6.schlegel_projection` for a more interesting projection.

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
        r"""
        Return a solid rendering of a 2- or 3-d polytope.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p_solid = p.render_solid(opacity=.7)                                  # needs sage.plot
            sage: type(p_solid)                                                         # needs sage.plot
            <class 'sage.plot.plot3d.index_face_set.IndexFaceSet'>
        """
        proj = self.projection()
        if self.ambient_dim() == 3:
            return proj.render_solid_3d(**kwds)
        if self.ambient_dim() == 2:
            return proj.render_fill_2d(**kwds)
        raise ValueError("render_solid is only defined for 2 and 3 dimensional polyhedra")

    def render_wireframe(self, **kwds):
        r"""
        For polytopes in 2 or 3 dimensions, return the edges
        as a list of lines.

        EXAMPLES::

            sage: p = Polyhedron([[1,2,],[1,1],[0,0]])
            sage: p_wireframe = p.render_wireframe()                                    # needs sage.plot
            sage: p_wireframe._objects                                                  # needs sage.plot
            [Line defined by 2 points, Line defined by 2 points, Line defined by 2 points]
        """
        proj = self.projection()
        if self.ambient_dim() == 3:
            return proj.render_wireframe_3d(**kwds)
        if self.ambient_dim() == 2:
            return proj.render_outline_2d(**kwds)
        raise ValueError("render_wireframe is only defined for 2 and 3 dimensional polyhedra")

    def schlegel_projection(self, facet=None, position=None):
        r"""
        Return the Schlegel projection.

        * The facet is orthonormally transformed into its affine hull.

        * The position specifies a point coming out of the barycenter of the
          facet from which the other vertices will be projected into the facet.

        INPUT:

        - ``facet`` -- a :class:`~sage.geometry.polyhedron.face.PolyhedronFace`
          The facet into which the Schlegel diagram is created. The default is the first facet

        - ``position`` -- a positive number. Determines a relative distance
          from the barycenter of ``facet``. A value close to 0 will place the
          projection point close to the facet and a large value further away.
          Default is `1`. If the given value is too large, an error is returned.

        OUTPUT: a :class:`~sage.geometry.polyhedron.plot.Projection` object

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: sch_proj = p.schlegel_projection()
            sage: schlegel_edge_indices = sch_proj.lines
            sage: schlegel_edges = [sch_proj.coordinates_of(x) for x in schlegel_edge_indices]
            sage: len([x for x in schlegel_edges if x[0][0] > 0])
            8

        The Schlegel projection preserves the convexity of facets, see :issue:`30015`::

            sage: fcube = polytopes.hypercube(4)
            sage: tfcube = fcube.face_truncation(fcube.faces(0)[0])
            sage: tfcube.facets()[-1]
            A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 8 vertices
            sage: sp = tfcube.schlegel_projection(tfcube.facets()[-1])
            sage: sp.plot()                                                             # needs sage.plot
            Graphics3d Object

        The same truncated cube but see inside the tetrahedral facet::

            sage: tfcube.facets()[4]
            A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 4 vertices
            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4])                   # needs sage.symbolic
            sage: sp.plot()                                                             # needs sage.plot sage.symbolic
            Graphics3d Object

        A different values of ``position`` changes the projection::

            sage: # needs sage.symbolic
            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4], 1/2)
            sage: sp.plot()                                                             # needs sage.plot
            Graphics3d Object
            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4], 4)
            sage: sp.plot()                                                             # needs sage.plot
            Graphics3d Object

        A value which is too large give a projection point that sees more than
        one facet resulting in a error::

            sage: sp = tfcube.schlegel_projection(tfcube.facets()[4], 5)
            Traceback (most recent call last):
            ...
            ValueError: the chosen position is too large
        """
        proj = self.projection()
        return proj.schlegel(facet, position)

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

        Check that :issue:`23355` is fixed::

            sage: P = Polyhedron([[7]]); P
            A 0-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex
            sage: P.affine_hull_projection()
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: P.affine_hull_projection(orthonormal='True')
            A 0-dimensional polyhedron in QQ^0 defined as the convex hull of 1 vertex
            sage: P.affine_hull_projection(orthogonal='True')
            A 0-dimensional polyhedron in QQ^0 defined as the convex hull of 1 vertex

        Check that :issue:`24047` is fixed::

            sage: P1 = Polyhedron(vertices=[[-1, 1], [0, -1], [0, 0], [-1, -1]])
            sage: P2 = Polyhedron(vertices=[[1, 1], [1, -1], [0, -1], [0, 0]])
            sage: P = P1.intersection(P2)
            sage: A, b = P.affine_hull_projection(as_affine_map=True,
            ....:                                 orthonormal=True, extend=True)

            sage: Polyhedron([(2,3,4)]).affine_hull_projection()
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex

        Check that backend is preserved::

            sage: polytopes.simplex(backend='field').affine_hull_projection().backend()
            'field'

            sage: P = Polyhedron(vertices=[[0,0], [1,0]], backend='field')
            sage: P.affine_hull_projection(orthogonal=True, orthonormal=True,
            ....:                          extend=True).backend()
            'field'

        Check that :issue:`29116` is fixed::

            sage: V = [[1, 0, -1, 0, 0],
            ....:      [1, 0, 0, -1, 0],
            ....:      [1, 0, 0, 0, -1],
            ....:      [1, 0, 0, +1, 0],
            ....:      [1, 0, 0, 0, +1],
            ....:      [1, +1, 0, 0, 0]]
            sage: P = Polyhedron(V)
            sage: P.affine_hull_projection()
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 6 vertices
            sage: P.affine_hull_projection(orthonormal=True)                            # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: the base ring needs to be extended; try with "extend=True"
            sage: P.affine_hull_projection(orthonormal=True, extend=True)               # needs sage.rings.number_field
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
                from sage.rings.qqbar import AA
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
        r"""
        Return the polyhedron projected into its affine hull.

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

        - ``as_polyhedron``, ``as_convex_set`` -- boolean or the default
          ``None``; one of the two to be set

        - ``as_affine_map`` -- boolean (default: ``False``); control the output

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

        - ``return_all_data`` -- boolean (default: ``False``)

          If set, then ``as_polyhedron`` and ``as_affine_map`` will set
          (possibly overridden) and additional (internal) data concerning
          the transformation is returned. Everything is encapsulated
          in an instance of
          :class:`~sage.geometry.convex_set.AffineHullProjectionData` in
          this case.

        - ``orthogonal`` -- boolean (default: ``False``); if ``True``,
          provide an orthogonal transformation

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

            sage: L = Polyhedron([[1,0], [0,1]]); L
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: A = L.affine_hull_projection(); A
            A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (0), A vertex at (1))
            sage: A = L.affine_hull_projection(orthogonal=True); A
            A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (0), A vertex at (2))
            sage: A = L.affine_hull_projection(orthonormal=True)                        # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            ValueError: the base ring needs to be extended; try with "extend=True"
            sage: A = L.affine_hull_projection(orthonormal=True, extend=True); A        # needs sage.rings.number_field
            A 1-dimensional polyhedron in AA^1 defined as the convex hull of 2 vertices
            sage: A.vertices()                                                          # needs sage.rings.number_field
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
            sage: A = S.affine_hull_projection(orthonormal=True, extend=True); A        # needs sage.rings.number_field
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 4 vertices
            sage: A.vertices()                                                          # needs sage.rings.number_field
            (A vertex at (0.7071067811865475?, 0.4082482904638630?, 1.154700538379252?),
             A vertex at (0.7071067811865475?, 1.224744871391589?, 0.?e-18),
             A vertex at (1.414213562373095?, 0.?e-18, 0.?e-18),
             A vertex at (0.?e-18, 0.?e-18, 0.?e-18))

        With the parameter ``minimal`` one can get a minimal base ring::

            sage: # needs sage.rings.number_field
            sage: s = polytopes.simplex(3)
            sage: s_AA = s.affine_hull_projection(orthonormal=True, extend=True)
            sage: s_AA.base_ring()
            Algebraic Real Field
            sage: s_full = s.affine_hull_projection(orthonormal=True, extend=True,
            ....:                                   minimal=True)
            sage: s_full.base_ring()
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1
             with a = 0.5176380902050415?

        More examples with the ``orthonormal`` parameter::

            sage: P = polytopes.permutahedron(3); P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: set([F.as_polyhedron().affine_hull_projection(                        # needs sage.combinat sage.rings.number_field
            ....:          orthonormal=True, extend=True).volume()
            ....:     for F in P.affine_hull_projection().faces(1)]) == {1, sqrt(AA(2))}
            True
            sage: set([F.as_polyhedron().affine_hull_projection(                        # needs sage.combinat sage.rings.number_field
            ....:          orthonormal=True, extend=True).volume()
            ....:     for F in P.affine_hull_projection(
            ....:             orthonormal=True, extend=True).faces(1)]) == {sqrt(AA(2))}
            True

            sage: # needs sage.rings.number_field
            sage: D = polytopes.dodecahedron()
            sage: F = D.faces(2)[0].as_polyhedron()
            sage: F.affine_hull_projection(orthogonal=True)
            A 2-dimensional polyhedron in
             (Number Field in sqrt5 with defining polynomial x^2 - 5
              with sqrt5 = 2.236067977499790?)^2
             defined as the convex hull of 5 vertices
            sage: F.affine_hull_projection(orthonormal=True, extend=True)
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 5 vertices

            sage: # needs sage.rings.number_field
            sage: K.<sqrt2> = QuadraticField(2)
            sage: P = Polyhedron([2*[K.zero()],2*[sqrt2]]); P
            A 1-dimensional polyhedron in
             (Number Field in sqrt2 with defining polynomial x^2 - 2
              with sqrt2 = 1.414213562373095?)^2
             defined as the convex hull of 2 vertices
            sage: P.vertices()
            (A vertex at (0, 0), A vertex at (sqrt2, sqrt2))
            sage: A = P.affine_hull_projection(orthonormal=True); A
            A 1-dimensional polyhedron in
             (Number Field in sqrt2 with defining polynomial x^2 - 2
              with sqrt2 = 1.414213562373095?)^1
             defined as the convex hull of 2 vertices
            sage: A.vertices()
            (A vertex at (0), A vertex at (2))

            sage: # needs sage.rings.number_field
            sage: K.<sqrt3> = QuadraticField(3)
            sage: P = Polyhedron([2*[K.zero()], 2*[sqrt3]]); P
            A 1-dimensional polyhedron in
             (Number Field in sqrt3 with defining polynomial x^2 - 3
              with sqrt3 = 1.732050807568878?)^2
             defined as the convex hull of 2 vertices
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

            sage: P.is_combinatorially_isomorphic(P.affine_hull_projection())           # needs sage.rings.number_field
            True
            sage: P.is_combinatorially_isomorphic(P.affine_hull_projection(             # needs sage.rings.number_field
            ....:     orthogonal=True))
            True
            sage: P.is_combinatorially_isomorphic(P.affine_hull_projection(             # needs sage.rings.number_field
            ....:     orthonormal=True, extend=True))
            True

        The ``orthonormal=True`` parameter preserves volumes;
        it provides an isometric copy of the polyhedron::

            sage: # needs sage.rings.number_field
            sage: Pentagon = polytopes.dodecahedron().faces(2)[0].as_polyhedron()
            sage: P = Pentagon.affine_hull_projection(orthonormal=True, extend=True)
            sage: _, c= P.is_inscribed(certificate=True)
            sage: c
            (0.4721359549995794?, 0.6498393924658126?)
            sage: circumradius = (c - vector(P.vertices()[0])).norm()
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

            sage: # needs sage.rings.number_field
            sage: Pentagon = polytopes.dodecahedron().faces(2)[0].as_polyhedron()
            sage: Pnormal = Pentagon.affine_hull_projection(orthonormal=True,
            ....:                                           extend=True)
            sage: Pgonal  = Pentagon.affine_hull_projection(orthogonal=True)
            sage: A, b    = Pentagon.affine_hull_projection(orthogonal=True,
            ....:                                           as_affine_map=True)
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

            sage: # needs sage.combinat sage.rings.number_field
            sage: P = polytopes.permutahedron(4)
            sage: Q    = P.affine_hull_projection(orthonormal=True, extend=True)
            sage: A, b = P.affine_hull_projection(orthonormal=True, extend=True,
            ....:                                 as_affine_map=True)
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
            NotImplementedError: "orthogonal=True" and "orthonormal=True"
            work only for compact polyhedra
            sage: P.affine_hull_projection(orthonormal=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: "orthogonal=True" and "orthonormal=True"
            work only for compact polyhedra

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
                    Domain:   Vector space of dimension 3 over Rational Field
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
                    Domain:   Vector space of dimension 3 over Rational Field
                    Codomain: Vector space of dimension 2 over Rational Field,
                projection_translation=(1, 1/2),
                section_linear_map=Vector space morphism represented by the matrix:
                    [-1/2  1/2    0]
                    [-1/3 -1/3  2/3]
                    Domain:   Vector space of dimension 2 over Rational Field
                    Codomain: Vector space of dimension 3 over Rational Field,
                section_translation=(1, 0, 0))

        The section map is a right inverse of the projection map::

            sage: mat = data.section_linear_map.matrix().transpose()
            sage: data.image.linear_transformation(mat) + data.section_translation == S
            True

        Same without ``orthogonal=True``::

            sage: data = S.affine_hull_projection(return_all_data=True); data
            AffineHullProjectionData(image=A 2-dimensional polyhedron in ZZ^2
                                           defined as the convex hull of 3 vertices,
                projection_linear_map=Vector space morphism represented by the matrix:
                    [1 0]
                    [0 1]
                    [0 0]
                    Domain:   Vector space of dimension 3 over Rational Field
                    Codomain: Vector space of dimension 2 over Rational Field,
                projection_translation=(0, 0),
                section_linear_map=Vector space morphism represented by the matrix:
                    [ 1  0 -1]
                    [ 0  1 -1]
                    Domain:   Vector space of dimension 2 over Rational Field
                    Codomain: Vector space of dimension 3 over Rational Field,
                section_translation=(0, 0, 1))
            sage: mat = data.section_linear_map.matrix().transpose()
            sage: data.image.linear_transformation(mat) + data.section_translation == S
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
        r"""
        Run tests on the method :meth:`.affine_hull_projection`.

        TESTS::

            sage: D = polytopes.dodecahedron()                                          # needs sage.rings.number_field
            sage: D.facets()[0].as_polyhedron()._test_affine_hull_projection()          # needs sage.rings.number_field
        """
        if tester is None:
            tester = self._tester(**options)

        if self.is_empty():
            # Undefined, nothing to test
            return

        if self.n_vertices() > 30 or self.n_facets() > 30 or self.dim() > 6:
            # Avoid very long doctests.
            return

        try:
            from sage.rings.qqbar import AA
        except ImportError:
            AA = None

        data_sets = []
        data_sets.append(self.affine_hull_projection(return_all_data=True))
        if self.is_compact():
            data_sets.append(self.affine_hull_projection(return_all_data=True,
                                                         orthogonal=True,
                                                         extend=True))
            if AA is not None:
                try:
                    data_sets.append(self.affine_hull_projection(return_all_data=True,
                                                                 orthonormal=True,
                                                                 extend=True))
                    data_sets.append(self.affine_hull_projection(return_all_data=True,
                                                                 orthonormal=True,
                                                                 extend=True,
                                                                 minimal=True))
                except ModuleNotFoundError:
                    pass

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

        - ``ambient_chart`` -- a chart on ``ambient_space``

        - ``names`` -- names for the coordinates on the affine hull

        - optional arguments accepted by :meth:`affine_hull_projection`

        The default chart is determined by the optional arguments of
        :meth:`affine_hull_projection`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: triangle = Polyhedron([(1, 0, 0), (0, 1, 0), (0, 0, 1)]);  triangle
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
            Riemannian metric gamma on the
             2-dimensional Riemannian submanifold A embedded in the Euclidean space E^3
            sage: A.volume_form()
            2-form eps_gamma on the
             2-dimensional Riemannian submanifold A embedded in the Euclidean space E^3

        Orthogonal version::

            sage: A = triangle.affine_hull_manifold(name='A', orthogonal=True); A       # needs sage.symbolic
            2-dimensional Riemannian submanifold A embedded in the Euclidean space E^3
            sage: A.embedding().display()                                               # needs sage.symbolic
            A → E^3
               (x0, x1) ↦ (x, y, z) = (t0 - 1/2*x0 - 1/3*x1 + 1,
                                       t0 + 1/2*x0 - 1/3*x1, t0 + 2/3*x1)
            sage: A.embedding().inverse().display()                                     # needs sage.symbolic
            E^3 → A
               (x, y, z) ↦ (x0, x1) = (-x + y + 1, -1/2*x - 1/2*y + z + 1/2)

        Arrangement of affine hull of facets::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: D = polytopes.dodecahedron()
            sage: E3 = EuclideanSpace(3)
            sage: submanifolds = [              # long time
            ....:     F.as_polyhedron().affine_hull_manifold(name=f'F{i}',
            ....:                                            orthogonal=True, ambient_space=E3)
            ....:     for i, F in enumerate(D.facets())]
            sage: sum(FM.plot({},               # long time, not tested                 # needs sage.plot
            ....:             srange(-2, 2, 0.1), srange(-2, 2, 0.1),
            ....:             opacity=0.2)
            ....:     for FM in submanifolds) + D.plot()
            Graphics3d Object

        Full-dimensional case::

            sage: cube = polytopes.cube(); cube
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: cube.affine_hull_manifold()                                           # needs sage.symbolic
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
        H = Manifold(self.dim(), name, ambient=ambient_space, structure='Riemannian',
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
