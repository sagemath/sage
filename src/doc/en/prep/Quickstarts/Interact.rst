.. -*- coding: utf-8 -*-

.. linkall

.. _prep-quickstart-interact:

Sage Interact Quickstart
========================

This `Sage <http://www.sagemath.org/>`_ quickstart tutorial was
developed for the MAA PREP Workshop "Sage: Using Open\-Source
Mathematics Software with Undergraduates" (funding provided by NSF DUE
0817071).

Invaluable resources are the Sage wiki
`http://wiki.sagemath.org/interact <http://wiki.sagemath.org/interact>`_
(type "sage interact" into Google), `UTMOST Sage Cell Repository <http://utmost-sage-cell.org/interacts>`_
(a collection of contributed interacts).

Start with just one command
---------------------------

How would one create an interactive cell?  First, let's focus on a new
thing to do!  Perhaps we just want a graph plotter that has some
options.

So let's start by getting the commands for what you want the output to
look like.  Here we just want a simple plot.

::

    sage: plot(x^2,(x,-3,3))
    Graphics object consisting of 1 graphics primitive

Then abstract out the parts you want to change.  We'll be letting the
user change the function, so let's make that a variable ``f``.

::

    sage: f=x^3
    sage: plot(f,(x,-3,3))
    Graphics object consisting of 1 graphics primitive

This was important because it allowed you to step back and think about
what you would really be doing.

Now for the technical part.  We make this a ``def`` function \- see the
:doc:`programming tutorial <../Programming>`.

::

    sage: def myplot(f=x^2):
    ....:     show(plot(f,(x,-3,3)))

.. note:
   The ``show`` or ``print`` is needed since the output is not
   automatically printed from within a function. Note also that we give
   the variable a default value of ``x^2``.  This is what ``f`` is if
   the user does not specify a value for ``f``.

Let's test the ``def`` function ``myplot`` by just calling it.

::

    sage: myplot()

If we call it with a different value for ``f``, we should get a
different plot.

::

    sage: myplot(x^3)

So far, we've only defined a new function, so this was review.  To make
a "control" to allow the user to interactively enter the function, we just preface the function with
``@interact``.

.. skip

::

    sage: @interact
    sage: def myplot(f=x^2):
    ....:     show(plot(f,(x,-3,3)))

.. note::
   Technically what ``@interact`` does is wrap the function, so the
   above is equivalent to::

       def myplot(..): ...
       myplot=interact(myplot)

Note that we can still call our function, even when we've used
``@interact``.  This is often useful in debugging it.

::

    sage: myplot(x^4)

Adding Complexity
-----------------

We can go ahead and replace other parts of the expression with
variables.  Note that ``_`` is the function name now. That is a just
convention for throw\-away names that we don't care about.

.. skip

::

    sage: @interact
    sage: def _(f=x^2, a=-3, b=3):
    ....:     show(plot(f,(x,a,b)))

If we pass ``('label', default_value)`` in for a control, then the
control gets the label when printed. Here, we've put in some text for
all three of them.  Remember that the text must be in quotes!  Otherwise
Sage will think that you are referring (for example) to some variable
called "lower", which it will think you forgot to define.

.. skip

::

    sage: @interact
    sage: def _(f=('$f$', x^2), a=('lower', -3), b=('upper', 3)):
    ....:     show(plot(f,(x,a,b)))

We can specify the type of control explicitly, along with options.
See :ref:`below <Control>` for more detail on the possibilities.

.. skip

::

    sage: @interact
    sage: def _(f=input_box(x^2, width=20, label="$f$")):
    ....:     show(plot(f,(x,-3,3)))

..
   Currently doesn't work.

   Here's another type of control: a color picker.

   .. skip

   ::

       sage: @interact
       sage: def _(f=input_box(x^2, width=20), color=color_selector()):
       ....:     show(plot(f,(x,-3,3), color=color))

Here we demonstrate a bunch of options.  Notice the new controls:

- ``range_slider``, which passes in  *two* values, ``zoom[0]`` and ``zoom[1]``

- ``True``/``False`` gets converted to checkboxes for the end user

.. skip

::

    sage: @interact
    sage: def _(f=input_box(x^2,width=20),
    ....: color=color_selector(widget='colorpicker', label=""),
    ....: axes=True,
    ....: fill=True,
    ....: zoom=range_slider(-3,3,default=(-3,3))):
    ....:     show(plot(f,(x,zoom[0], zoom[1]), color=color, axes=axes,fill=fill))

There is also one button type to :ref:`disable automatic updates <NoUpdate>`.

The previous interact was a bit ugly, because all of the controls were
stacked on top of each other. We can control the layout of the widget
controls in a grid (at the top, bottom, left, or right) using the
``layout`` parameter.

.. skip

::

    sage: @interact(layout=dict(top=[['f', 'color']],
    ....: left=[['axes'],['fill']],
    ....: bottom=[['zoom']]))
    sage: def _(f=input_box(x^2,width=20),
    ....: color=color_selector(widget='colorpicker', label=""),
    ....: axes=True,
    ....: fill=True,
    ....: zoom=range_slider(-3,3, default=(-3,3))):
    ....:     show(plot(f,(x,zoom[0], zoom[1]), color=color, axes=axes,fill=fill))

.. _Control:

Control Types
-------------

There are many potential types of widgets one might want to use for
interactive control.  Sage has all of the following:

- boxes

- sliders

- range sliders

- checkboxes

- selectors (dropdown lists or buttons)

- grid of boxes

- color selectors

- plain text

We illustrate some more of these below.

.. skip

::

    sage: @interact
    sage: def _(frame=checkbox(True, label='Use frame')):
    ....:     show(plot(sin(x), (x,-5,5)), frame=frame)

.. skip

::

    sage: var('x,y')
    sage: colormaps=sage.plot.colors.colormaps.keys()
    sage: @interact
    sage: def _(cmap=selector(colormaps)):
    ....:     contour_plot(x^2-y^2,(x,-2,2),(y,-2,2),cmap=cmap).show()

.. skip

::

    sage: var('x,y')
    sage: colormaps=sage.plot.colors.colormaps.keys()
    sage: @interact
    sage: def _(cmap=selector(['RdBu', 'jet', 'gray','gray_r'],buttons=True),
    sage: type=['density','contour']):
    ....:     if type=='contour':
    ....:         contour_plot(x^2-y^2,(x,-2,2),(y,-2,2),cmap=cmap, aspect_ratio=1).show()
    ....:     else:
    ....:         density_plot(x^2-y^2,(x,-2,2),(y,-2,2),cmap=cmap, frame=True,axes=False,aspect_ratio=1).show()

By default, ranges are sliders that divide the range into 50 steps.

.. skip

::

    sage: @interact
    sage: def _(n=(1,20)):
    ....:     print(factorial(n))

You can set the step size to get, for example, just integer values.

.. skip

::

    sage: @interact
    sage: def _(n=slider(1,20, step_size=1)):
    ....:     print(factorial(n))

Or you can explicitly specify the slider values.

.. skip

::

    sage: @interact
    sage: def _(n=slider([1..20])):
    ....:     print(factorial(n))

And the slider values don't even have to be numbers!

.. skip

::

    sage: @interact
    sage: def _(fun=('function', slider([sin,cos,tan,sec,csc,cot]))):
    ....:     print(fun(4.39293))

Matrices are automatically converted to a grid of input boxes.

.. skip

::

    sage: @interact
    sage: def _(m=('matrix', identity_matrix(2))):
    ....:     print(m.eigenvalues())

Here's how to get vectors from a grid of boxes.

.. skip

::

    sage: @interact
    sage: def _(v=('vector', input_grid(1, 3, default=[[1,2,3]], to_value=lambda x: vector(flatten(x))))):
    ....:     print(v.norm())

.. _NoUpdate:

The option not to update
------------------------

As a final problem, what happens when the controls get so complicated
that it would counterproductive to see the interact update for each of
the changes one wants to make?  Think changing the endpoints and order
of integration for a triple integral, for instance, or the example below
where a whole matrix might be changed.

In this situation, where we don't want any updates until we specifically
say so, we can use the ``auto_update=False`` option.  This will create a
button to enable the user to update as soon as he or she is ready.

.. skip

::

    sage: @interact
    sage: def _(m=('matrix', identity_matrix(2)), auto_update=False):
    ....:     print(m.eigenvalues())

