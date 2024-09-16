.. sage-doctest: needs sage.plot sage.symbolic

.. linkall

.. _prep-symbolics-and-basic-plotting:

Tutorial for Symbolics and Plotting
===================================

This `Sage <http://www.sagemath.org>`_ document is one of the tutorials
developed for the MAA PREP Workshop "Sage: Using Open\-Source
Mathematics Software with Undergraduates" (funding provided by NSF DUE
0817071).  It is licensed under the Creative Commons
Attribution\-ShareAlike 3.0 license (`CC BY\-SA
<http://creativecommons.org/licenses/by-sa/3.0/>`_).

This tutorial has the following sections:

- :ref:`Symbolic`

- :ref:`2DPlotting`

- :ref:`3DPlotting`

It assumes that one is familiar with the absolute basics of functions
and evaluation in Sage.  We provide a (very) brief refresher.

#. Make sure the syntax below for defining a function and getting a
   value makes sense.

#. Then evaluate the cell by clicking the "evaluate" link, or by
   pressing :kbd:`Shift` + :kbd:`Enter` (hold down :kbd:`Shift`
   while pressing the :kbd:`Enter` key).

::

    sage: f(x)=x^3+1
    sage: f(2)
    9

.. _Symbolic:

Symbolic Expressions
--------------------

In the first tutorial, we defined *functions* using notation similar to
that one would use in (say) a calculus course.

There is a useful variant on this \- defining *expressions* involving
variables.  This will give us the opportunity to point out several
important, and sometimes subtle, things.

In the cell below, we define an expression :math:`FV` which is the
future value of an investment of \$100, compounded continuously.  We
then substitute in values for :math:`r` and :math:`t` which calculate
the future value for :math:`t=5` years and :math:`r=5\%` nominal
interest.

::

    sage: var('r,t')
    (r, t)
    sage: FV=100*e^(r*t)

::

    sage: FV(r=.05,t=5)
    128.402541668774

The previous cells point out several things to remember when working
with symbolic expressions.  Some are fairly standard.

- An asterisk (``*``) signifies multiplication.  This should be how you
  always do multiplication.

 - Although it is possible to allow implicit multiplication, this can
   easily lead to ambiguity.

- We can access the most important constants; for instance, :math:`e`
  stands for the constant :math:`2.71828...`.  Likewise, ``pi`` (or
  :math:`\pi`) and :math:`I` (think complex numbers) are also defined.

 - Of course, if you redefine :math:`e` to be something else, all bets
   are off!

However, two others may be unfamiliar, especially if you have not used
much mathematical software before.

- You must tell Sage what the variables are before using them in a
  symbolic expression.

 - We did that above by typing ``var('r,t')``.

 - This is automatically done with the :math:`f(x)` notation, but
   without that it is necessary, so that Sage knows one intends ``t`` (for
   instance) is a symbolic variable and not a number or something else.

- If you then wish to substitute some values into the expression, you
  must explicitly tell Sage which variables are being assigned which
  values.

 - For instance, above we used ``FV(r=.05,t=5)`` to indicate the precise
   values of :math:`r` and :math:`t`.

Notice that when we define a function, we don't need to specify which
variable has which value.  In the function defined below, we have
already specified an order.

::

    sage: FV2(r,t)=100*e^(r*t)
    sage: FV2(.05,5)
    128.402541668774

In this case it is clear that :math:`r` is first and :math:`t` is second.

But with ``FV=100*e^(r*t)``, there is no particular reason :math:`r` or
:math:`t` should be first.

::

    sage: FV(r=.05,t=5); FV(t=5,r=.05)
    128.402541668774
    128.402541668774

Also remember that when we don't use function notation, we'll need to
define our variables.

One of the great things we can do with expressions is manipulate them.
Let's make a typical expression.

::

    sage: z = (x+1)^3

In the cells below, you'll notice something new: the character ``#``.
In Sage (and in `Python <http://www.python.org>`_ ), anything on a
single line after the number/pound sign (the `octothorp
<http://en.wikipedia.org/wiki/Number_sign>`_ ) is ignored.  We say that
``#`` is a comment character.  We use it below to mention alternative
ways to do the same thing.

::

    sage: expand(z) # or z.expand()
    x^3 + 3*x^2 + 3*x + 1

::

    sage: y = expand(z)
    sage: y.factor() # or factor(y)
    (x + 1)^3

In the previous cell, we *assigned* the expression which is the
expansion of :math:`z` to the variable :math:`y` with the first line.
After that, anything we want to do to the expansion of :math:`z` can be
done by doing it to :math:`y`.

There are more commands like this as well.  Notice that :math:`z` will
no longer be :math:`(x+1)^3` after this cell is evaluated, since we've
assigned :math:`z` to a (much more complex) expression.

::

    sage: z = ((x - 1)^(3/2) - (x + 1)*sqrt(x - 1))/sqrt((x - 1)*(x + 1))
    sage: z.simplify_full()
    -2*sqrt(x - 1)/sqrt(x^2 - 1)

This is a good place for a few reminders of basic help.

- You can see various methods for simplifying an expression by using tab
  completion.  Put your cursor at the end of the next cell (after the
  ``simplify``) and press tab to see lots of different methods.

- Also remember that you can use the question mark (e.g.,
  ``z.simplify_rational?``) to get help about a particular method.

::

    sage: z.simplify
    <built-in method simplify of sage.symbolic.expression.Expression object at ...>

Finally, recall that you can get nicely typeset versions of the output
in several ways.

- One option is to click the 'Typeset' button at the top.

- Another \- which does not require scrolling! \- is to use the ``show``
  command.

.. skip

::

    sage: show(z.simplify_rational())

.. MATH::

    -\frac{2 \, \sqrt{x - 1}}{\sqrt{x^{2} - 1}}

Another Sage command that is useful in this context is ``solve``.

Here, we solve the simple equation :math:`x^2=-1`.

::

    sage: solve(x^2==-1,x) # solve x^2==-1 for x
    [x == -I, x == I]

- In the ``solve`` command, one types an equals sign in the equation as
  *two* equal signs.

 - This is because the single equals sign means assignment to a
   variable, as we've done above, so Sage (along with Python) uses the
   double equals sign for symbolic equality.

- We also include the variable we'd like to solve for after the comma.

It's also possible to solve more than one expression simultaneously.

::

    sage: solve([x^2==1,x^3==1],x)
    [[x == 1]]

.. _2DPlotting:

Basic 2D Plotting
-----------------

One of the other basic uses of mathematics software is easy plotting.
Here, we include a brief introduction to the sorts of plotting which
will prepare us to use Sage in calculus.  (There will be a separate
tutorial for more advanced plotting techniques.)

Recall that we can generate a plot using fairly simple syntax.  Here, we
define a function of :math:`x` and plot it between :math:`-1` and
:math:`1`.

::

    sage: f(x)=x^3+1
    sage: plot(f,(x,-1,1))
    Graphics object consisting of 1 graphics primitive

We can give the plot a name, so that if we want to do something with the
plot later, we don't have to type out the entire plot command.
Remember, this is called *assigning* the plot to the name/variable.

In the next cell, we give the plot the name :math:`P`.

::

    sage: P=plot(f,(x,-1,1))

One plot is nice, but one might want to superimpose plots as well.  For
instance, the tangent line to :math:`f` at :math:`x=0` is just the line
:math:`y=1`, and we might want to show this together with the plot.

So let's plot this line in a different color, and with a different style
for the line, but over the same interval.

::

    sage: Q=plot(1,(x,-1,1),color="red", linestyle="--")
    sage: Q
    Graphics object consisting of 1 graphics primitive

Because we put :math:`Q` in a line by itself at the end, it shows.  We
were able to use just one cell to define :math:`Q` and show it, by
putting each command in a separate line in the same input cell.

Now to show the plots superimposed on each other, we simply add them.

::

    sage: P+Q
    Graphics object consisting of 2 graphics primitives

Suppose we wanted to view a detail of this.

- We could create another plot with different endpoints.

- Another way is to keep the currently created plots, but to set the
  viewing window using the ``show`` command, as below.

::

    sage: (P+Q).show(xmin=-.1,xmax=.1,ymin=.99,ymax=1.01)

Since the axes no longer cross in the frame of reference, Sage shows a
short gap between the horizontal and vertical axes.

There are many options one can pass in for various purposes.

- Some require quotes around the values.

 - Such as the ``color`` option when we made a ``"red"`` line.

- Some do not.

 - Such as the ``xmin`` in the previous plot, where the minimum x value
   was just :math:`-.1`

Usually (though not always) quotes are required for option values which
are words or strings of characters, and not required for numerical
values.

Two of the most useful of these options help in labeling graphs.

- The ``axes_labels`` option labels the axes.

 - As with the word processor, we can use dollar signs (like in LaTeX)
   to make the labels typeset nicely.

 - Here we need both quotes and brackets for proper syntax, since there
   are two axes to label and the labels are not actually numbers.

::

    sage: plot(f,(x,-1,1),axes_labels=['$x$','$y$'],legend_label='$f(x)$',show_legend=True)
    Graphics object consisting of 1 graphics primitive

- The ``legend_label`` option is especially useful with multiple plots.

 - LaTeX notation works here too.

 - In the graphic above, we needed to explicitly ask to show the label.
   With multiple graphs this should not be necessary.

::

    sage: P1 = plot(f,(x,-1,1),axes_labels=['$x$','$y$'],legend_label='$f(x)$')
    sage: P2 = plot(sin,(x,-1,1),axes_labels=['$x$','$y$'],legend_label=r'$\sin(x)$',color='red')
    sage: P1+P2
    Graphics object consisting of 2 graphics primitives

One additional useful note is that plots of functions with vertical
asymptotes may need their vertical viewing range set manually; otherwise
the asymptote may really go to infinity!

::

    sage: plot(1/x^2,(x,-10,10),ymax=10)
    Graphics object consisting of 1 graphics primitive

Remember, you can use the command ``plot?`` to find out about most of
the options demonstrated above.

Below, you can experiment with several of the plotting options.

- Just evaluate the cell and play with the sliders, buttons, color
  picker, etc., to change the plot options.

- You can access low\-level options like the initial number of plotted
  points, or high\-level ones like whether axes are shown or not.

- This uses a feature of Sage called "interacts", which is a very
  powerful way to engage students in exploring a problem.

.. skip

::

    sage: x = var('x')
    sage: @interact
    sage: def plot_example(f=sin(x^2),r=range_slider(-5,5,step_size=1/4,default=(-3,3)),
    ....:                  color=color_selector(widget='colorpicker'),
    ....:                  thickness=(3,(1..10)),
    ....:                  adaptive_recursion=(5,(0..10)), adaptive_tolerance=(0.01,(0.001,1)),
    ....:                  plot_points=(20,(1..100)),
    ....:                  linestyle=['-','--','-.',':'],
    ....:                  gridlines=False, fill=False,
    ....:                  frame=False, axes=True
    ....:                  ):
    ....:     show(plot(f, (x,r[0],r[1]), color=color, thickness=thickness,
    ....:                  adaptive_recursion=adaptive_recursion,
    ....:                  adaptive_tolerance=adaptive_tolerance, plot_points=plot_points,
    ....:                  linestyle=linestyle, fill=fill if fill else None),
    ....:                  gridlines=gridlines, frame=frame, axes=axes)

.. _3DPlotting:

Basic 3D Plotting
-----------------

Plotting a 3D plot is similar to plotting a 2D plot, but we need to
specify ranges for two variables instead of one.

::

    sage: g(x,y)=sin(x^2+y^2)
    sage: plot3d(g,(x,-5,5),(y,-5,5))
    Graphics3d Object

There is a lot you can do with the 3D plots.

- Try rotating the plot above by clicking and dragging the mouse inside
  of the plot.

- Also, right\-click (Control\-click if you have only one mouse button)
  just to the right of the plot to see other options in a menu.

- If you have a wheel on your mouse or a multi\-touch trackpad, you can
  scroll to zoom.

- You can also right\-click to see other options, such as

 - spinning the plot,

 - changing various colors,

 - and even making the plot suitable for viewing through 3D glasses
   (under the "style", then "stereographic" submenus),

When using the ``plot3d`` command, the first variable range specified is
plotted along the usual "x" axis, while the second range specified is
plotted along the usual "y" axis.

The plot above is somewhat crude because the function is not sampled
enough times \- this is fairly rapidly changing function, after all.  We
can make the plot smoother by telling Sage to sample the function using
a grid of 300 by 300 points.  Sage then samples the function at 90,000
points!

::

    sage: plot3d(g,(x,-5,5),(y,-5,5),plot_points=300)
    Graphics3d Object

As with 2D plots, we can superimpose 3D plots by adding them together.

Note that in this one, we do not define the functions, but only use
expressions (see the first set of topics in this tutorial), so it is
wisest to define the variables ahead of time.

::

    sage: var('x,y')
    (x, y)
    sage: b = 2.2
    sage: P=plot3d(sin(x^2-y^2),(x,-b,b),(y,-b,b), opacity=.7)
    sage: Q=plot3d(0, (x,-b,b), (y,-b,b), color='red')
    sage: P+Q
    Graphics3d Object

As usual, only the last command shows up in the notebook, though clearly
all are evaluated.  This also demonstrates that many of the same options
work for 3D plots as for 2D plots.

We close this tutorial with a cool plot that we define *implicitly* as a
3D contour plot.

::

    sage: var('x,y,z')
    (x, y, z)
    sage: T = golden_ratio
    sage: p = 2 - (cos(x + T*y) + cos(x - T*y) + cos(y + T*z) + cos(y - T*z) + cos(z - T*x) + cos(z + T*x))
    sage: r = 4.78
    sage: implicit_plot3d(p, (x, -r, r), (y, -r, r), (z, -r, r), plot_points=50, color='yellow')
    Graphics3d Object

The next tutorial will use all that you have learned about Sage basics,
symbolics, and plotting in a specific mathematical venue \- the
:ref:`calculus sequence <prep-calculus>`!

