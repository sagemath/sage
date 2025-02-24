.. sage-doctest: needs sage.plot sage.symbolic

.. _section-plot:

Plotting
========

Sage can produce two-dimensional and three-dimensional plots.

Two-dimensional Plots
---------------------

In two dimensions, Sage can draw circles, lines, and polygons;
plots of functions in rectangular coordinates; and also polar
plots, contour plots and vector field plots. We present examples of
some of these here. For more examples of plotting with Sage, see
:ref:`section-systems` and :ref:`section-maxima`, and also the
`Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_
documentation.

This command produces a yellow circle of radius 1, centered at the
origin:

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0))
    Graphics object consisting of 1 graphics primitive

You can also produce a filled circle:

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0), fill=True)
    Graphics object consisting of 1 graphics primitive

You can also create a circle by assigning it to a variable; this
does not plot it:

::

    sage: c = circle((0,0), 1, rgbcolor=(1,1,0))

To plot it, use ``c.show()`` or ``show(c)``, as follows:

.. link

::

    sage: c.show()

Alternatively, evaluating ``c.save('filename.png')`` will save the
plot to the given file.

Now, these 'circles' look more like ellipses because the axes are
scaled differently.  You can fix this:

.. link

::

    sage: c.show(aspect_ratio=1)

The command ``show(c, aspect_ratio=1)`` accomplishes the same
thing, or you can save the picture using
``c.save('filename.png', aspect_ratio=1)``.

It's easy to plot basic functions:

::

    sage: plot(cos, (-5,5))
    Graphics object consisting of 1 graphics primitive

Once you specify a variable name, you can create parametric plots
also:

::

    sage: x = var('x')
    sage: parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    Graphics object consisting of 1 graphics primitive

It's important to notice that the axes of the plots will only
intersect if the origin is in the viewing range of the graph, and
that with sufficiently large values scientific notation may be used:

::

    sage: plot(x^2,(x,300,500))
    Graphics object consisting of 1 graphics primitive

You can combine several plots by adding them:

::

    sage: x = var('x')
    sage: p1 = parametric_plot((cos(x),sin(x)),(x,0,2*pi),rgbcolor=hue(0.2))
    sage: p2 = parametric_plot((cos(x),sin(x)^2),(x,0,2*pi),rgbcolor=hue(0.4))
    sage: p3 = parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    sage: show(p1+p2+p3, axes=false)

A good way to produce filled-in shapes is to produce a list of
points (``L`` in the example below) and then use the ``polygon``
command to plot the shape with boundary formed by those points. For
example, here is a green deltoid:

::

    sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),
    ....:     2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))
    sage: p
    Graphics object consisting of 1 graphics primitive

Type ``show(p, axes=false)`` to see this without any axes.

You can add text to a plot:

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),
    ....:     6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t)

Calculus teachers draw the following plot frequently on the board:
not just one branch of arcsin but rather several of them: i.e., the
plot of :math:`y=\sin(x)` for :math:`x` between :math:`-2\pi`
and :math:`2\pi`, flipped about the 45 degree line. The following
Sage commands construct this:

::

    sage: v = [(sin(x),x) for x in srange(-2*float(pi),2*float(pi),0.1)]
    sage: line(v)
    Graphics object consisting of 1 graphics primitive

Since the tangent function has a larger range than sine, if you use
the same trick to plot the inverse tangent, you should change the
minimum and maximum coordinates for the *x*-axis:

::

    sage: v = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
    sage: show(line(v), xmin=-20, xmax=20)

Sage also computes polar plots, contour plots and vector field plots
(for special types of functions). Here is an example of a contour
plot:

::

    sage: f = lambda x,y: cos(x*y)
    sage: contour_plot(f, (-4, 4), (-4, 4))
    Graphics object consisting of 1 graphics primitive

Three-Dimensional Plots
-----------------------

Sage can also be used to create three-dimensional plots. In both
the notebook and the REPL, these plots will be displayed by default
using the open source package [ThreeJS]_, which supports interactively
rotating and zooming the figure with the mouse.

Use ``plot3d`` to graph a function of the form `f(x, y) = z`:

::

    sage: x, y = var('x,y')
    sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
    Graphics3d Object

Alternatively, you can use ``parametric_plot3d`` to graph a
parametric surface where each of `x, y, z` is determined by
a function of one or two variables (the parameters, typically
`u` and `v`). The previous plot can be expressed parametrically
as follows:

::

    sage: u, v = var('u, v')
    sage: f_x(u, v) = u
    sage: f_y(u, v) = v
    sage: f_z(u, v) = u^2 + v^2
    sage: parametric_plot3d([f_x, f_y, f_z], (u, -2, 2), (v, -2, 2))
    Graphics3d Object

The third way to plot a 3D surface in Sage is ``implicit_plot3d``,
which graphs a contour of a function like `f(x, y, z) = 0` (this
defines a set of points). We graph a sphere using the classical
formula:

::

    sage: x, y, z = var('x, y, z')
    sage: implicit_plot3d(x^2 + y^2 + z^2 - 4, (x,-2, 2), (y,-2, 2), (z,-2, 2))
    Graphics3d Object

Here are some more examples:

`Yellow Whitney's umbrella <http://en.wikipedia.org/wiki/Whitney_umbrella>`__:

::

    sage: u, v = var('u,v')
    sage: fx = u*v
    sage: fy = u
    sage: fz = v^2
    sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1),
    ....:   frame=False, color="yellow")
    Graphics3d Object

`Cross cap <http://en.wikipedia.org/wiki/Cross-cap>`__:

::

    sage: u, v = var('u,v')
    sage: fx = (1+cos(v))*cos(u)
    sage: fy = (1+cos(v))*sin(u)
    sage: fz = -tanh((2/3)*(u-pi))*sin(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object

Twisted torus:

::

    sage: u, v = var('u,v')
    sage: fx = (3+sin(v)+cos(u))*cos(2*v)
    sage: fy = (3+sin(v)+cos(u))*sin(2*v)
    sage: fz = sin(u)+2*cos(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object

Lemniscate:

::

    sage: x, y, z = var('x,y,z')
    sage: f(x, y, z) = 4*x^2 * (x^2 + y^2 + z^2 + z) + y^2 * (y^2 + z^2 - 1)
    sage: implicit_plot3d(f, (x, -0.5, 0.5), (y, -1, 1), (z, -1, 1))
    Graphics3d Object
