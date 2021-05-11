r"""
Introduction

Sage has a wide support for 3D graphics, from basic shapes to implicit and
parametric plots.

The following graphics functions are supported:


-  :func:`~plot3d` - plot a 3d function

-  :func:`~sage.plot.plot3d.parametric_plot3d.parametric_plot3d` - a parametric three-dimensional space curve or surface

-  :func:`~sage.plot.plot3d.revolution_plot3d.revolution_plot3d` - a plot of a revolved curve

-  :func:`~sage.plot.plot3d.plot_field3d.plot_vector_field3d` - a plot of a 3d vector field

-  :func:`~sage.plot.plot3d.implicit_plot3d.implicit_plot3d` - a plot of an isosurface of a function

-  :func:`~sage.plot.plot3d.list_plot3d.list_plot3d`- a 3-dimensional plot of a surface defined by a list of points in 3-dimensional space

-  :func:`~sage.plot.plot3d.list_plot3d.list_plot3d_matrix` - a 3-dimensional plot of a surface defined by a matrix defining points in 3-dimensional space

-  :func:`~sage.plot.plot3d.list_plot3d.list_plot3d_array_of_arrays`- A 3-dimensional plot of a surface defined by a list of lists defining points in 3-dimensional space

-  :func:`~sage.plot.plot3d.list_plot3d.list_plot3d_tuples` - a 3-dimensional plot of a surface defined by a list of points in 3-dimensional space

The following classes for basic shapes are supported:


-  :class:`~sage.plot.plot3d.shapes.Box` - a box given its three magnitudes

-  :class:`~sage.plot.plot3d.shapes.Cone` - a cone, with base in the xy-plane pointing up the z-axis

-  :class:`~sage.plot.plot3d.shapes.Cylinder` - a cylinder, with base in the xy-plane pointing up the z-axis

-  :class:`~sage.plot.plot3d.shapes2.Line` - a 3d line joining a sequence of points

-  :class:`~sage.plot.plot3d.shapes.Sphere` - a sphere centered at the origin

-  :class:`~sage.plot.plot3d.shapes.Text` - a text label attached to a point in 3d space

-  :class:`~sage.plot.plot3d.shapes.Torus` - a 3d torus

-  :class:`~sage.plot.plot3d.shapes2.Point` - a position in 3d, represented by a sphere of fixed size


The following plotting functions for basic shapes are supported


-  :func:`~sage.plot.plot3d.shapes.ColorCube` - a cube with given size and sides with given colors

-  :func:`~sage.plot.plot3d.shapes.LineSegment` - a line segment, which is drawn as a cylinder from start to end with given radius

-  :func:`~sage.plot.plot3d.shapes2.line3d` - a 3d line joining a sequence of points

-  :func:`~sage.plot.plot3d.shapes.arrow3d` - a 3d arrow

-  :func:`~sage.plot.plot3d.shapes2.point3d` - a point or list of points in 3d space

-  :func:`~sage.plot.plot3d.shapes2.bezier3d` - a 3d bezier path

-  :func:`~sage.plot.plot3d.shapes2.frame3d` - a frame in 3d

-  :func:`~sage.plot.plot3d.shapes2.frame_labels` - labels for a given frame in 3d

-  :func:`~sage.plot.plot3d.shapes2.polygon3d` - draw a polygon in 3d

-  :func:`~sage.plot.plot3d.shapes2.polygons3d` - draw the union of several polygons in 3d

-  :func:`~sage.plot.plot3d.shapes2.ruler` - draw a ruler in 3d, with major and minor ticks

-  :func:`~sage.plot.plot3d.shapes2.ruler_frame` - draw a frame made of 3d rulers, with major and minor ticks

-  :func:`~sage.plot.plot3d.shapes2.sphere` - plot of a sphere given center and radius

-  :func:`~sage.plot.plot3d.shapes2.text3d` - 3d text

Sage also supports platonic solids with the following functions:


-  :func:`~sage.plot.plot3d.platonic.tetrahedron`

-  :func:`~sage.plot.plot3d.platonic.cube`

-  :func:`~sage.plot.plot3d.platonic.octahedron`

-  :func:`~sage.plot.plot3d.platonic.dodecahedron`

-  :func:`~sage.plot.plot3d.platonic.icosahedron`

Different viewers are supported: a web-based interactive viewer using the
Three.js JavaScript library (the default), Jmol, and the Tachyon ray tracer.
The viewer is invoked by adding the keyword argument
``viewer='threejs'`` (respectively ``'jmol'`` or ``'tachyon'``)
to the command ``show()`` on any three-dimensional graphic.


-  :class:`~sage.plot.plot3d.tachyon.Tachyon` - create a scene the can be rendered using the Tachyon ray tracer

-  :class:`~sage.plot.plot3d.tachyon.Axis_aligned_box` - box with axis-aligned edges with the given min and max coordinates

-  :class:`~sage.plot.plot3d.tachyon.Cylinder` - an infinite cylinder

-  :class:`~sage.plot.plot3d.tachyon.FCylinder` - a finite cylinder

-  :class:`~sage.plot.plot3d.tachyon.FractalLandscape`- axis-aligned fractal landscape

-  :class:`~sage.plot.plot3d.tachyon.Light` - represents lighting objects

-  :class:`~sage.plot.plot3d.tachyon.ParametricPlot` - parametric plot routines

-  :class:`~sage.plot.plot3d.tachyon.Plane` - an infinite plane

-  :class:`~sage.plot.plot3d.tachyon.Ring` - an annulus of zero thickness

-  :class:`~sage.plot.plot3d.tachyon.Sphere`- a sphere

-  :class:`~sage.plot.plot3d.tachyon.TachyonSmoothTriangle` - a triangle along with a normal vector, which is used for smoothing

-  :class:`~sage.plot.plot3d.tachyon.TachyonTriangle` - basic triangle class

-  :class:`~sage.plot.plot3d.tachyon.TachyonTriangleFactory` - class to produce triangles of various rendering types

-  :class:`~sage.plot.plot3d.tachyon.Texfunc` - creates a texture function

-  :class:`~sage.plot.plot3d.tachyon.Texture` - stores texture information

-  :func:`~sage.plot.plot3d.tachyon.tostr` - converts vector information to a space-separated string

"""
