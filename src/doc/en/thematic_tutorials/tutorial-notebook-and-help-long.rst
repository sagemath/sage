.. _tutorial-notebook-and-help-long:

==============================================================================
Tutorial: Using the Sage notebook, navigating the help system, first exercises
==============================================================================

.. linkall

This worksheet is based on William Stein's `JPL09__intro_to_sage.sws
<http://modular.math.washington.edu/talks/20090701-sage_graphics_tutorial/JPL09___intro_to_sage.sws>`_
worksheet and the `Sage days 20.5_demo <http://wiki.sagemath.org/days20.5>`_
worksheet and aims to be an interactive introduction to Sage through exercises.
You will learn how to use the notebook and call the help.

Making this help page into a worksheet
======================================

If you are browsing this document as a static web page, you can see all the
examples; however you need to copy-paste them one by one to experiment with them.
Use the ``Upload worksheet`` button of the notebook and copy-paste the URL of
this page to obtain an editable copy in your notebook.

If you are browsing this document as part of Sage's live documentation, you can
play with the examples directly here; however your changes will be lost when
you close this page. Use ``Copy worksheet`` from the ``File...`` menu at the
top of this page to get an editable copy in your notebook.

Both in the live tutorial and in the notebook, you can clear all output by
selecting ``Delete All Output`` from the ``Action...`` menu next to the
``File...`` menu at the top of the worksheet.

Entering, Editing and Evaluating Input
======================================

To *evaluate code* in the Sage Notebook, type the code into an input cell and
press ``shift-enter`` or click the ``evaluate`` link. Try it now with a simple
expression (e.g., `2+3`). The first time you evaluate a cell takes longer
than subsequent times since a new Sage process is started::

    sage: 2 + 3
    5

    sage: # edit here

    sage: # edit here

To create *new input cells*, click the blue line that appears between
cells when you move your mouse around. Try it now::

    sage: 1 + 1
    2

    sage: # edit here

You can *go back* and edit any cell by clicking in it (or using the
arrow keys on your keyboard to move up or down). Go back and change
your `2+3` above to `3+3` and re-evaluate it. An empty cell can be
*deleted* with backspace.

You can also *edit this text* right here by double clicking on it,
which will bring up the TinyMCE Javascript text editor. You can even
put embedded mathematics like this $\sin(x) - y^3$ by using dollar signs
just like in TeX or LaTeX.

Help systems
============

There are various ways of getting help in Sage.

- navigate through the documentation (there is a link ``Help`` at the top right
  of the worksheet),
- ``tab`` completion,
- contextual help.

We detail below the latter two methods through examples.

Completion and contextual documentation
=======================================

Start typing something and press the ``tab`` key. The interface tries to
complete it with a command name. If there is more than one completion, then
they are all presented to you. Remember that Sage is case sensitive, i.e. it
differentiates upper case from lower case. Hence the ``tab`` completion of
``klein`` won't show you the ``KleinFourGroup`` command that builds the group
`\ZZ/2 \times \ZZ/2` as a permutation group. Try it on the next cells:

.. skip

::

    sage: klein<tab>

    sage: Klein<tab>

To see documentation and examples for a command, type a question mark ``?`` at
the end of the command name and press the ``tab`` key as in:

.. skip

::

    sage: KleinFourGroup?<tab>

::

    sage: # edit here

.. TOPIC:: Exercise A

    What is the largest prime factor of `600851475143`?

    .. skip

    ::

        sage: factor?<tab>

    ::

        sage: # edit here

In the above manipulations we have not stored any data for
later use. This can be done in Sage with the ``=`` symbol as in::

    sage: a = 3
    sage: b = 2
    sage: print a+b
    5

This can be understood as Sage evaluating the expression to the right
of the ``=`` sign and creating the appropriate object, and then
associating that object with a label, given by the left-hand side (see
the foreword of :ref:`tutorial-objects-and-classes` for
details). Multiple assignments can be done at once::

    sage: a,b = 2,3
    sage: print a,b
    2 3

This allows us to swap the values of two variables directly::

    sage: a,b = 2,3
    sage: a,b = b,a
    sage: print a,b
    3 2

We can also assign a common value to several variables simultaneously::

    sage: c = d = 1
    sage: c, d
    (1, 1)
    sage: d = 2
    sage: c, d
    (1, 2)

Note that when we use the word *variable* in the computer-science sense we
mean "a label attached to some data stored by Sage". Once an object is
created, some *methods* apply to it. This means *functions* but instead of
writing **f(my_object)** you write **my_object.f()**::

    sage: p = 17
    sage: p.is_prime()
    True

See :ref:`tutorial-objects-and-classes` for details.
To know all methods of an object you can once more use tab-completion. Write the
name of the object followed by a dot and then press ``tab``:

.. skip

::

    sage: a.<tab>

    sage: # edit here

.. TOPIC:: Exercise B

    Create the permutation 51324 and assign it to the variable ``p``.

    .. skip

    ::

        sage: Permutation?<tab>

    ::

        sage: # edit here


    What is the ``inverse`` of ``p``?

    .. skip

    ::

        sage: p.inv<tab>

        sage: # edit here

    Does ``p`` have the ``pattern`` 123? What about 1234? And 312? (even if you don't
    know what a pattern is, you should be able to find a command that does this).

    .. skip

    ::

        sage: p.pat<tab>

        sage: # edit here

Some linear algebra
===================

.. TOPIC:: Exercise C

    Use the :func:`matrix` command to create the following matrix.

    .. MATH::

        M = \left(\begin{array}{rrrr}
        10 & 4 & 1 & 1 \\
        4 & 6 & 5 & 1 \\
        1 & 5 & 6 & 4 \\
        1 & 1 & 4 & 10
        \end{array}\right)

    .. skip

    ::

        sage: matrix?<tab>

    ::

        sage: # edit here

    Then, using methods of the matrix,

    1. Compute the determinant of the matrix.
    2. Compute the echelon form of the matrix.
    3. Compute the eigenvalues of the matrix.
    4. Compute the kernel of the matrix.
    5. Compute the LLL decomposition of the matrix (and lookup the
       documentation for what LLL is if needed!)

    ::

        sage: # edit here

        sage: # edit here

    Now that you know how to access the different methods of matrices,

    6. Create the vector `v = (1,-1,-1,1)`.
    7. Compute the two products: `M\cdot v` and `v\cdot M`. What mathematically
       borderline operation is Sage doing implicitly?

    .. skip

    ::

        sage: vector?<tab>

    ::

        sage: # edit here

.. NOTE::

    Vectors in Sage are row vectors. A method such as ``eigenspaces`` might not
    return what you expect, so it is best to specify ``eigenspaces_left`` or
    ``eigenspaces_right`` instead. Same thing for kernel (``left_kernel`` or
    ``right_kernel``), and so on.


Some Plotting
=============

The :func:`plot` command allows you to draw plots of functions. Recall
that you can access the documentation by pressing the ``tab`` key
after writing ``plot?`` in a cell:

.. skip

::

    sage: plot?<tab>

::

    sage: # edit here

Here is a simple example::

    sage: var('x')   # make sure x is a symbolic variable
    x
    sage: plot(sin(x^2), (x,0,10))
    Graphics object consisting of 1 graphics primitive

Here is a more complicated plot. Try to change every single input to the plot
command in some way, evaluating to see what happens::

    sage: P = plot(sin(x^2), (x,-2,2), rgbcolor=(0.8,0,0.2), thickness=3, linestyle='--', fill='axis')
    sage: show(P, gridlines=True)

Above we used the :func:`show` command to show a plot after it was created. You can
also use ``P.show`` instead::

    sage: P.show(gridlines=True)

Try putting the cursor right after ``P.show(`` and pressing tab to get a list of
the options for how you can change the values of the given inputs.

.. skip

::

    sage: P.show(

Plotting multiple functions at once is as easy as adding them together::

    sage: P1 = plot(sin(x), (x,0,2*pi))
    sage: P2 = plot(cos(x), (x,0,2*pi), rgbcolor='red')
    sage: P1 + P2
    Graphics object consisting of 2 graphics primitives

Symbolic Expressions
====================

Here is an example of a symbolic function::

    sage: f(x) = x^4 - 8*x^2 - 3*x + 2
    sage: f(x)
    x^4 - 8*x^2 - 3*x + 2

    sage: f(-3)
    20

This is an example of a function in the *mathematical* variable `x`. When Sage
starts, it defines the symbol `x` to be a mathematical variable. If you want
to use other symbols for variables, you must define them first::

    sage: x^2
    x^2
    sage: u + v
    Traceback (most recent call last):
    ...
    NameError: name 'u' is not defined

    sage: var('u v')
    (u, v)
    sage: u + v
    u + v

Still, it is possible to define symbolic functions without first
defining their variables::

    sage: f(w) = w^2
    sage: f(3)
    9

In this case those variables are defined implicitly::

    sage: w
    w

.. TOPIC:: Exercise D

    Define the symbolic function `f(x) = x \sin(x^2)`. Plot `f` on the
    domain `[-3,3]` and color it red. Use the :func:`find_root` method to
    numerically approximate the root of `f` on the interval `[1,2]`::

        sage: # edit here

    Compute the tangent line to `f` at `x=1`::

        sage: # edit here

    Plot `f` and the tangent line to `f` at `x=1` in one image::

        sage: # edit here

.. TOPIC:: Exercise E (Advanced)

     Solve the following equation for `y`:

    .. MATH::

        y = 1 + x y^2

    There are two solutions, take the one for which `\lim_{x\to0}y(x)=1`.
    (Don't forget to create the variables `x` and `y`!).

    ::

        sage: # edit here

    Expand `y` as a truncated Taylor series around `0` and containing
    `n=10` terms.

    ::

        sage: # edit here

    Do you recognize the coefficients of the Taylor series expansion? You might
    want to use the `On-Line Encyclopedia of Integer Sequences
    <http://oeis.org>`_, or better yet, Sage's class :class:`OEIS` which
    queries the encyclopedia:

    .. skip

    ::


        sage: oeis?<tab>

    ::

        sage: # edit here

Congratulations for completing your first Sage tutorial!
