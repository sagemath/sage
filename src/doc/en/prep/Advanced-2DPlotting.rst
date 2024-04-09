.. -*- coding: utf-8 -*-

.. linkall

.. _prep-advanced-2dplotting:

Tutorial for Advanced 2d Plotting
=================================

This `Sage <http://www.sagemath.org>`_ document is one of the tutorials
developed for the MAA PREP Workshop "Sage: Using Open\-Source
Mathematics Software with Undergraduates" (funding provided by NSF DUE
0817071).  It is licensed under the Creative Commons
Attribution\-ShareAlike 3.0 license (`CC BY\-SA
<http://creativecommons.org/licenses/by-sa/3.0/>`_).

Thanks to Sage's integration of projects like `matplotlib
<http://matplotlib.sourceforge.net/>`_, Sage has comprehensive
two-dimensional plotting capabilities.  This worksheet consists of the
following sections:

- :ref:`Cartesian`

- :ref:`Parametric`

- :ref:`Polar`

- :ref:`Data`

- :ref:`Contour`

- :ref:`Vector`

- :ref:`Complex`

- :ref:`Region`

- :ref:`Builtin`

- :ref:`Saving`

This tutorial assumes that one is familiar with the basics of Sage,
such as evaluating a cell by clicking the "evaluate" link, or by
pressing :kbd:`Shift` + :kbd:`Enter` (hold down :kbd:`Shift` while
pressing the :kbd:`Enter` key).

.. fixme - if log plots are in by the time this makes it in, put them in!!!

Graphing Functions and Plotting Curves
--------------------------------------

.. _Cartesian:

Cartesian Plots
~~~~~~~~~~~~~~~

A simple quadratic is easy.

::

    sage: plot(x^2, (x,-2,2))
    Graphics object consisting of 1 graphics primitive

You can combine "plot objects" by adding them.

::

    sage: regular = plot(x^2, (x,-2,2), color= 'purple')
    sage: skinny = plot(4*x^2, (x,-2,2), color = 'green')
    sage: regular + skinny
    Graphics object consisting of 2 graphics primitives

**Problem** : Plot a green :math:`y=\sin(x)` together with a red
:math:`y=2\,\cos(x)`.  (Hint: you can use ``pi`` as part of your range.)

Boundaries of a plot can be specified, in addition to the overall size.

::

    sage: plot(1+e^(-x^2), xmin=-2, xmax=2, ymin=0, ymax=2.5, figsize=10)
    Graphics object consisting of 1 graphics primitive

**Problem** : Plot :math:`y=5+3\,\sin(4x)` with suitable boundaries.

You can add lots of extra information.

::

    sage: exponential = plot(1+e^(-x^2), xmin=-2, xmax=2, ymin=0, ymax=2.5)
    sage: max_line = plot(2, xmin=-2, xmax=2, linestyle='-.', color = 'red')
    sage: min_line = plot(1, xmin=-2, xmax=2, linestyle=':', color = 'red')
    sage: exponential + max_line + min_line
    Graphics object consisting of 3 graphics primitives

You can fill regions with transparent color, and thicken the curve.
This example uses several options to fine\-tune our graphic.

::

    sage: exponential = plot(1+e^(-x^2), xmin=-2, xmax=2, ymin=0, ymax=2.5, fill=0.5, fillcolor='grey', fillalpha=0.3)
    sage: min_line = plot(1, xmin=-2, xmax=2, linestyle='-', thickness= 6, color = 'red')
    sage: exponential + min_line
    Graphics object consisting of 3 graphics primitives

::

    sage: sum([plot(x^n,(x,0,1),color=rainbow(5)[n]) for n in [0..4]])
    Graphics object consisting of 5 graphics primitives

**Problem** : Create a plot showing the cross-section area for the
following solid of revolution problem: Consider the area bounded by
:math:`y=x^2-3x+6` and the line :math:`y=4`.  Find the volume created by
rotating this area around the line :math:`y=1`.

.. _Parametric:

Parametric Plots
~~~~~~~~~~~~~~~~

A parametric plot needs a list of two functions of the parameter; in
Sage, we use *square* brackets to delimit the list.  Notice also that we
must declare ``t`` as a variable first.  Because the graphic is slightly
wider than it is tall, we use the ``aspect_ratio`` option (such options
are called *keywords* ) to ensure the axes are correct for how we want
to view this object.

::

    sage: t = var('t')
    sage: parametric_plot([cos(t) + 3 * cos(t/9), sin(t) - 3 * sin(t/9)], (t, 0, 18*pi), fill = True, aspect_ratio=1)
    Graphics object consisting of 2 graphics primitives

**Problem** : These parametric equations will create a hypocycloid.

.. MATH::

    x(t)=17\cos(t)+3\cos(17t/3)

.. MATH::

    y(t)=17\sin(t)-3\sin(17t/3)

Create this as a parametric plot.

Sage automatically plots a 2d or 3d plot, and a curve or a surface,
depending on how many variables and coordinates you specify.

::

    sage: t = var('t')
    sage: parametric_plot((t^2,sin(t)), (t,0,pi))
    Graphics object consisting of 1 graphics primitive

::

    sage: parametric_plot((t^2,sin(t),cos(t)), (t,0,pi))
    Graphics3d Object

::

    sage: r = var('r')
    sage: parametric_plot((t^2,sin(r*t),cos(r*t)), (t,0,pi),(r,-1,1))
    Graphics3d Object

.. _Polar:

Polar Plots
~~~~~~~~~~~

Sage can also do polar plots.

::

    sage: polar_plot(2 + 2*cos(x), (x, 0, 2*pi), color=hue(0.5), thickness=4)
    Graphics object consisting of 1 graphics primitive

Although they aren't essential, many of these examples try to
demonstrate things like coloring, fills, and shading to give you a sense
of the possibilities.

More than one polar curve can be specified in a list (square brackets).
Notice the automatic graded shading of the fill color.

::

    sage: t = var('t')
    sage: polar_plot([cos(4*t) + 1.5,  0.5 * cos(4*t) + 2.5], (t, 0, 2*pi),
    ....:            color='black', thickness=2, fill=True, fillcolor='orange')
    Graphics object consisting of 4 graphics primitives

Problem: Create a plot for the following problem. Find the area that is
inside the circle :math:`r=2`, but outside the cardioid
:math:`2+2\cos(\theta)`.

Interactive Demonstration
~~~~~~~~~~~~~~~~~~~~~~~~~

It may be of interest to see all these things put together in a very
nice pedagogical graphic.  Even though this is fairly advanced, and so
you may want to skip the code, it is not as difficult as you might think
to put together.

.. skip

::

    sage: html('<h2>Sine and unit circle (by Jurgis Pralgauskis)</h2> inspired by <a href="http://www.youtube.com/watch?v=Ohp6Okk_tww&feature=related">this video</a>' )
    sage: # http://doc.sagemath.org/html/en/reference/sage/plot/plot.html
    sage: radius = 100 # scale for radius of "unit" circle
    sage: graph_params = dict(xmin = -2*radius,    xmax = 360,
    ....:                    ymin = -(radius+30), ymax = radius+30,
    ....:                    aspect_ratio=1,
    ....:                    axes = False
    ....:                    )
    sage: def sine_and_unit_circle( angle=30, instant_show = True, show_pi=True ):
    ....:     ccenter_x, ccenter_y = -radius, 0  # center of circle on real coords
    ....:     sine_x = angle # the big magic to sync both graphs :)
    ....:     current_y = circle_y = sine_y = radius * sin(angle*pi/180)
    ....:     circle_x = ccenter_x + radius * cos(angle*pi/180)
    ....:     graph = Graphics()
    ....:     # we'll put unit circle and sine function on the same graph
    ....:     # so there will be some coordinate mangling ;)
    ....:     # CIRCLE
    ....:     unit_circle = circle((ccenter_x, ccenter_y), radius, color="#ccc")
    ....:     # SINE
    ....:     x = var('x')
    ....:     sine = plot( radius * sin(x*pi/180) , (x, 0, 360), color="#ccc" )
    ....:     graph += unit_circle + sine
    ....:     # CIRCLE axis
    ....:     # x axis
    ....:     graph +=  arrow( [-2*radius, 0], [0, 0], color = "#666" )
    ....:     graph += text("$(1, 0)$",  [-16, 16],  color = "#666")
    ....:     # circle y axis
    ....:     graph +=  arrow( [ccenter_x,-radius], [ccenter_x, radius], color = "#666" )
    ....:     graph += text("$(0, 1)$",  [ccenter_x, radius+15],  color = "#666")
    ....:     # circle center
    ....:     graph += text("$(0, 0)$",  [ccenter_x, 0],  color = "#666")
    ....:     # SINE x axis
    ....:     graph +=  arrow( [0,0], [360, 0], color = "#000" )
    ....:     # let's set tics
    ....:     # or http://aghitza.org/posts/tweak_labels_and_ticks_in_2d_plots_using_matplotlib/
    ....:     # or wayt for https://github.com/sagemath/sage/issues/1431
    ....:     # ['$-\pi/3$', '$2\pi/3$', '$5\pi/3$']
    ....:     for x in range(0, 361, 30):
    ....:         graph += point( [x, 0] )
    ....:         angle_label = ".  $%3d^{\circ}$ " % x
    ....:         if show_pi: angle_label += " $(%s\pi) $"% x/180
    ....:         graph += text(angle_label,  [x, 0], rotation=-90,
    ....:         vertical_alignment='top', fontsize=8, color="#000" )
    ....:     # CURRENT VALUES
    ....:     # SINE -- y
    ....:     graph +=  arrow( [sine_x,0], [sine_x, sine_y], width=1, arrowsize=3)
    ....:     graph +=  arrow( [circle_x,0], [circle_x, circle_y], width=1, arrowsize=3)
    ....:     graph +=  line(([circle_x, current_y], [sine_x, current_y]), rgbcolor="#0F0", linestyle = "--", alpha=0.5)
    ....:     # LABEL on sine
    ....:     graph += text("$(%d^{\circ}, %3.2f)$"%(sine_x, float(current_y)/radius),  [sine_x+30, current_y],  color = "#0A0")
    ....:     # ANGLE -- x
    ....:     # on sine
    ....:     graph += arrow( [0,0], [sine_x, 0], width=1, arrowsize=1, color='red')
    ....:     # on circle
    ....:     graph += disk( (ccenter_x, ccenter_y), float(radius)/4, (0, angle*pi/180), color='red', fill=False, thickness=1)
    ....:     graph +=  arrow( [ccenter_x, ccenter_y], [circle_x, circle_y],
    ....:                  rgbcolor="#cccccc", width=1, arrowsize=1)
    ....:     if instant_show:
    ....:         show (graph,  **graph_params)
    ....:     return graph
    sage: #####################
    sage: # make Interaction
    sage: ######################
    sage: @interact
    sage: def _( angle = slider([0..360], default=30, step_size=5,
    ....:          label="Pasirinkite kampą:    ", display_value=True) ):
    ....:     sine_and_unit_circle(angle, show_pi = False)

Plotting Data
-------------

.. _Data:

Plotting Data Points
~~~~~~~~~~~~~~~~~~~~

Sometimes one wishes to simply plot data.  Here, we demonstrate several
ways of plotting points and data via the simple approximation to the
Fibonacci numbers given by

.. MATH::

    F_n=\frac{1}{\sqrt{5}}\left(\frac{1+\sqrt{5}}{2}\right)^n\; ,

which is quite good after about :math:`n=5`.

First, we notice that the Fibonacci numbers are built in.

::

    sage: fibonacci_sequence(6)
    <generator object fibonacci_sequence at ...>

::

    sage: list(fibonacci_sequence(6))
    [0, 1, 1, 2, 3, 5]

The ``enumerate`` command is useful for taking a list and coordinating
it with the counting numbers.

::

    sage: list(enumerate(fibonacci_sequence(6)))
    [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 5)]

So we just define the numbers and coordinate pairs we are about to plot.

::

    sage: fibonacci = list(enumerate(fibonacci_sequence(6)))
    sage: f(n)=(1/sqrt(5))*((1+sqrt(5))/2)^n
    sage: asymptotic = [(i, f(i)) for i in range(6)]
    sage: fibonacci
    [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 5)]
    sage: asymptotic
    [(0, 1/5*sqrt(5)), (1, 1/10*sqrt(5)*(sqrt(5) + 1)), (2, 1/20*sqrt(5)*(sqrt(5) + 1)^2), (3, 1/40*sqrt(5)*(sqrt(5) + 1)^3), (4, 1/80*sqrt(5)*(sqrt(5) + 1)^4), (5, 1/160*sqrt(5)*(sqrt(5) + 1)^5)]

Now we can plot not just the two sets of points, but also use several of
the documented options for plotting points. Those coming from other
systems may prefer ``list_plot``.

::

    sage: fib_plot=list_plot(fibonacci, color='red', pointsize=30)
    sage: asy_plot = list_plot(asymptotic, marker='D',color='black',thickness=2,plotjoined=True)
    sage: show(fib_plot+asy_plot, aspect_ratio=1)

Other options include ``line``, ``points``, and ``scatter_plot``.
Having the choice of markers for different data is particularly helpful
for generating publishable graphics.

::

    sage: fib_plot=scatter_plot(fibonacci, facecolor='red', marker='o',markersize=40)
    sage: asy_plot = line(asymptotic, marker='D',color='black',thickness=2)
    sage: show(fib_plot+asy_plot, aspect_ratio=1)

Contour\-type Plots
-------------------

.. _Contour:

Contour Plots
~~~~~~~~~~~~~

Contour plotting can be very useful when trying to get a handle on
multivariable functions, as well as modeling.  The basic syntax is
essentially the same as for 3D plotting \- simply an extension of the 2D
plotting syntax.

::

    sage: f(x,y)=y^2+1-x^3-x
    sage: contour_plot(f, (x,-pi,pi), (y,-pi,pi))
    Graphics object consisting of 1 graphics primitive

We can change colors, specify contours, label curves, and many other
things.  When there are many levels, the ``colorbar`` keyword becomes
quite useful for keeping track of them.  Notice that, as opposed to many
other options, it can only be ``True`` or ``False`` (corresponding to
whether it appears or does not appear).

::

    sage: contour_plot(f, (x,-pi,pi), (y,-pi,pi),colorbar=True,labels=True)
    Graphics object consisting of 1 graphics primitive

This example is fairly self\-explanatory, but demonstrates the power of
formatting, labeling, and the wide variety of built\-in color gradations
(colormaps or ``cmap``).  The strange\-looking construction
corresponding to ``label_fmt`` is a Sage/Python data type called a
*dictionary* , and turns out to be useful for more advanced Sage use; it
consists of pairs connected by a colon, all inside curly braces.

::

    sage: contour_plot(f, (x,-pi,pi), (y,-pi,pi), contours=[-4,0,4], fill=False,
    ....:     cmap='cool', labels=True, label_inline=True, 
    ....:     label_fmt={-4:"low", 0:"medium", 4: "hi"}, label_colors='black')
    Graphics object consisting of 1 graphics primitive

Implicit plots are a special type of contour plot (they just plot the
zero contour).

::

    sage: f(x,y)
    -x^3 + y^2 - x + 1

::

    sage: implicit_plot(f(x,y)==0,(x,-pi,pi),(y,-pi,pi))
    Graphics object consisting of 1 graphics primitive

A density plot is like a contour plot, but without discrete levels.

::

    sage: density_plot(f, (x, -2, 2), (y, -2, 2))
    Graphics object consisting of 1 graphics primitive

Sometimes contour plots can be a little misleading (which makes for a
*great* classroom discussion about the problems of ignorantly relying on
technology). Here we combine a density plot and contour plot to show
even better what is happening with the function.

::

    sage: density_plot(f,(x,-2,2),(y,-2,2))+contour_plot(f,(x,-2,2),(y,-2,2),fill=False,labels=True,label_inline=True,cmap='jet')
    Graphics object consisting of 2 graphics primitives

It can be worth getting familiar with the various options for different
plots, especially if you will be doing a lot of them in a given
worksheet or pedagogical situation.

Here are the options for contour plots.

- They are given as an "attribute" \- no parentheses \- of the
  ``contour_plot`` object.

- They are given as a dictionary (see :ref:`the programming tutorial
  <Advanced>`).

::

    sage: contour_plot.options
    {'aspect_ratio': 1,
     'axes': False,
     'colorbar': False,
     'contours': None,
     'fill': True,
     'frame': True,
     'labels': False,
     'legend_label': None,
     'linestyles': None,
     'linewidths': None,
     'plot_points': 100,
     'region': None}

Let's change it so that all future contour plots don't have the fill.
That's how some of us might use them in a class.  We'll also check that
the change happened.

::

    sage: contour_plot.options["fill"]=False
    sage: contour_plot.options
    {'aspect_ratio': 1,
     'axes': False,
     'colorbar': False,
     'contours': None,
     'fill': False,
     'frame': True,
     'labels': False,
     'legend_label': None,
     'linestyles': None,
     'linewidths': None,
     'plot_points': 100,
     'region': None}

And it works!

::

    sage: contour_plot(f,(x,-2,2),(y,-2,2))
    Graphics object consisting of 1 graphics primitive

We can always access the default options, of course, to remind us.

::

    sage: contour_plot.defaults()
    {'aspect_ratio': 1,
     'axes': False,
     'colorbar': False,
     'contours': None,
     'fill': True,
     'frame': True,
     'labels': False,
     'legend_label': None,
     'linestyles': None,
     'linewidths': None,
     'plot_points': 100,
     'region': None}

.. _Vector:

Vector fields
~~~~~~~~~~~~~

The syntax for vector fields is very similar to other multivariate
constructions.  Notice that the arrows are scaled appropriately, and
colored by length in the 3D case.

::

    sage: var('x,y')
    (x, y)
    sage: plot_vector_field((-y+x,y*x),(x,-3,3),(y,-3,3))
    Graphics object consisting of 1 graphics primitive

::

    sage: var('x,y,z')
    (x, y, z)
    sage: plot_vector_field3d((-y,-z,x), (x,-3,3),(y,-3,3),(z,-3,3))
    Graphics3d Object

3d vector field plots are ideally viewed with 3d glasses (right\-click
on the plot and select "Style" and "Stereographic")

.. _Complex:

Complex Plots
~~~~~~~~~~~~~

We can plot functions of complex variables, where the magnitude is
indicated by the brightness (black is zero magnitude) and the argument
is indicated by the hue (red is a positive real number).

::

    sage: f(z) = exp(z) #z^5 + z - 1 + 1/z
    sage: complex_plot(f, (-5,5),(-5,5))
    Graphics object consisting of 1 graphics primitive

.. _Region:

Region plots
~~~~~~~~~~~~

These plot where an expression is true, and are useful for plotting inequalities.

::

    sage: region_plot(cos(x^2+y^2) <= 0, (x, -3, 3), (y, -3, 3),aspect_ratio=1)
    Graphics object consisting of 1 graphics primitive

We can get fancier options as well.

::

    sage: region_plot(sin(x)*sin(y) >= 1/4, (x,-10,10), (y,-10,10), incol='yellow', bordercol='black', borderstyle='dashed', plot_points=250,aspect_ratio=1)
    Graphics object consisting of 2 graphics primitives

Remember, what command would give full information about the syntax,
options, and examples?

Miscellaneous Plot Information
------------------------------

.. _Builtin:

Builtin Graphics Objects
~~~~~~~~~~~~~~~~~~~~~~~~

Sage includes a variety of built\-in graphics objects.  These are
particularly useful for adding to one's plot certain objects which are
difficult to describe with equations, but which are basic geometric
objects nonetheless.  In this section we will try to demonstrate the
syntax of some of the most useful of them; for most of them the
contextual (remember, append ``?``) help will give more details.

Points
######

To make one point, a coordinate pair suffices.

::

    sage: point((3,5))
    Graphics object consisting of 1 graphics primitive

It doesn't matter how multiple point are generated; they must go
in as input via a list (square brackets).  Here, we demonstrate the
hard (but naive) and easy (but a little more sophisticated) way to
do this.

::

    sage: f(x)=x^2
    sage: points([(0,f(0)), (1,f(1)), (2,f(2)), (3,f(3)), (4,f(4))])
    Graphics object consisting of 1 graphics primitive

::

    sage: points([(x,f(x)) for x in range(5)])
    Graphics object consisting of 1 graphics primitive

Sage tries to tell how many dimensions you are working in automatically.

::

    sage: f(x,y)=x^2-y^2
    sage: points([(x,y,f(x,y)) for x in range(5) for y in range(5)])
    Graphics3d Object

Lines
#####

The syntax for lines is the same as that for points, but you get...
well, you get connecting lines too!

::

    sage: f(x)=x^2
    sage: line([(x,f(x)) for x in range(5)])
    Graphics object consisting of 1 graphics primitive

Balls
#####

Sage has disks and spheres of various types available.  Generally the
center and radius are all that is needed, but other options are
possible.

::

    sage: circle((0,1),1,aspect_ratio=1)
    Graphics object consisting of 1 graphics primitive

::

    sage: disk((0,0), 1, (pi, 3*pi/2), color='yellow',aspect_ratio=1)
    Graphics object consisting of 1 graphics primitive

There are also ellipses and various arcs; see the `full plot
documentation <http://doc.sagemath.org/html/en/reference/plotting/index.html>`_.

Arrows
######

::

    sage: arrow((0,0), (1,1))
    Graphics object consisting of 1 graphics primitive

Polygons
########

Polygons will try to complete themselves and fill in the interior;
otherwise the syntax is fairly self\-evident.

::

    sage: polygon([[0,0],[1,1],[1,2]])
    Graphics object consisting of 1 graphics primitive

Text
####

In 2d, one can typeset mathematics using the ``text`` command.  This can
be used to fine-tune certain types of labels.  Unfortunately, in 3D the
text is just text.

::

    sage: text(r'$\int_0^2 x^2\, dx$', (0.5,2))+plot(x^2,(x,0,2),fill=True)
    Graphics object consisting of 3 graphics primitives

.. _Saving:

Saving Plots
~~~~~~~~~~~~

We can save 2d plots to many different formats.  Sage can determine the
format based on the filename for the image.

::

    sage: p=plot(x^2,(x,-1,1))
    sage: p
    Graphics object consisting of 1 graphics primitive

For testing purposes, we use the Sage standard temporary filename;
however, you could use any string for a name that you wanted, like
``"my_plot.png"``.

::

    sage: name = tmp_filename() # this is a string
    sage: png_savename = name+'.png'
    sage: p.save(png_savename)

In the notebook, these are usually ready for downloading in little links
by the cells.

::

    sage: pdf_savename = name+'.pdf'
    sage: p.save(pdf_savename)

Notably, we can export in formats ready for inclusion in web pages.

::

    sage: svg_savename = name+'.svg'
    sage: p.save(svg_savename)

