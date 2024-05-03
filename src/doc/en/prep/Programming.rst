.. -*- coding: utf-8 -*-

.. linkall

.. _prep-programming:

Sage Introductory Programming Tutorial
======================================

This `Sage <http://www.sagemath.org>`_ document is one of the tutorials
developed for the MAA PREP Workshop "Sage: Using Open\-Source
Mathematics Software with Undergraduates" (funding provided by NSF DUE
0817071). It is licensed under the Creative Commons
Attribution\-ShareAlike 3.0 license (`CC BY\-SA
<http://creativecommons.org/licenses/by-sa/3.0/>`_).

This tutorial will cover the following topics (and various others
throughout):

- :ref:`Methods`

- :ref:`Lists`

- :ref:`Defs`

- :ref:`Gotchas`

- :ref:`Advanced`

We will motivate our examples using basic matrices and situations one
might want to handle in everyday use of matrices in the classroom.

.. _Methods:

Methods and Dot Notation
-------------------------

Making a new matrix is not too hard in Sage.

::

    sage: A = matrix([[1,2],[3,4]])

As we can see, this gives the matrix with rows given in the two
bracketed sets.

::

    sage: A
    [1 2]
    [3 4]

Some commands are available right off the bat, like derivatives are.
This is the determinant.

::

    sage: det(A)
    -2

But some things are not available this way \- for instance a
row\-reduced echelon form.  We can 'tab' after this to make sure.

.. skip

::

    sage: r

So, as we've already seen in previous tutorials, *many* of the commands
in Sage are "methods" of objects.

That is, we access them by typing:

- the name of the mathematical object,

- a dot/period,

- the name of the method, and

- parentheses (possibly with an argument).

This is a huge advantage, once you get familiar with it, because it
allows you to do *only* the things that are possible, and *all* such
things.

First, let's do the determinant again.

::

    sage: A.det()
    -2

Then we do the row\-reduced echelon form.

::

    sage: A.rref()
    [1 0]
    [0 1]

It is very important to keep in the parentheses.

.. note::
   Things that would be legal without them would be called 'attributes',
   but Sage prefers stylistically to hide them, since math is made of functions
   and not elements of sets. Or so a category\-theorist would say.

   ::

       sage: A.det # Won't work
       <built-in method det of sage.matrix.matrix_integer_dense.Matrix_integer_dense object at ...>

This is so useful because we can use the 'tab' key, remember!

.. skip

::

    sage: A.

Sometimes you will have surprises.  Subtle changes in an object can
affect what commands are available, or what their outcomes are.

::

    sage: A.echelon_form()
    [1 0]
    [0 2]

This is because our original matrix had only integer coefficients, and
you can't make the last entry one via elementary operations unless you
multiply by a rational number!

::

    sage: B = A.change_ring(QQ); B.echelon_form()
    [1 0]
    [0 1]

Another question is whether one needs an argument.  Remember, it's easy
to just read the documentation!

Below, let's see whether we need an argument to get a column.

.. skip

::

    sage: A.column?

It looks like we do.  Let's input ``1``.

::

    sage: A.column(1)
    (2, 4)

Notice that this gives the SECOND column!

::

    sage: A
    [1 2]
    [3 4]

What is that about?

.. _Lists:

Lists, Loops, and Set Builders
------------------------------

.. rubric:: (Especially List Comprehensions!)

In the previous example, we saw that the ``1`` choice for the column of
a matrix gives the *second* column.

::

    sage: matrix([[1,2],[3,4]]).column(1)
    (2, 4)

You might have thought that the would give the first column, but Sage
(along with the Python programming language) begins numbering of
anything that is like a sequence at zero.  We've mentioned this once
before, but it's very important to remember.

To reinforce this, let's formally introduce a fundamental object we've
seen once or twice before, called a *list*.

You should think of a list as an ordered set, where the elements of the
set can be pretty much anything \- including other lists.

::

    sage: my_list=[2,'Grover',[3,2,1]]; my_list
    [2, 'Grover', [3, 2, 1]]

You can access any elements of such a list quite easily using square
brackets.  Just remember that the counting starts at zero.

::

    sage: my_list[0]; my_list[2]
    2
    [3, 2, 1]

There are lots of advanced things one can do with lists.

::

    sage: my_list[0:2]
    [2, 'Grover']

However, our main reason for introducing this is more practical, as
we'll now see.

- One of the best uses of the computer in the classroom is to quickly
  show tedious things.

- One of the most tedious things to do by hand in linear algebra is
  taking powers of matrices.

- Here we make the first four powers of our matrix 'by hand'.

::

    sage: A = matrix([[1,2],[3,4]])
    sage: A^0; A^1; A^2; A^3; A^4
    [1 0]
    [0 1]
    [1 2]
    [3 4]
    [ 7 10]
    [15 22]
    [ 37  54]
    [ 81 118]
    [199 290]
    [435 634]

This is not terrible, but it's not exactly nice either, particularly if
you might want to do something *with* these new matrices.

Instead, we can do what is known as a *loop* construction.  See the
notation below; it's at least vaguely mathematical.

::

    sage: for i in [0,1,2,3,4]:
    ....:     A^i
    [1 0]
    [0 1]
    [1 2]
    [3 4]
    [ 7 10]
    [15 22]
    [ 37  54]
    [ 81 118]
    [199 290]
    [435 634]

What did we do?

- For each :math:`i` in the set :math:`\{0,1,2,3,4\}`, return :math:`A^i`.

Yeah, that makes sense.  The square brackets created a list, and the
powers of the original matrix come in the same order as the list.

(The colon in the first line and the indentation in the second line are
**extremely** important; they are the basic syntactical structure of
Python.)

For the curious: this is better, but still not perfect.  It would be
best to find a quicker way to write the possible values for
:math:`i`.  There are two ways to do this in Sage.

::

    sage: for i in [0..4]:
    ....:     det(A^i)
    1
    -2
    4
    -8
    16

::

    sage: for i in range(5):
    ....:     det(A^i)
    1
    -2
    4
    -8
    16

These ways of constructing lists are very useful \- and demonstrate
that, like many Sage/Python things, that counting begins at zero and
ends at one less than the "end" in things like ``range``.

Below, we show that one can get step sizes other than one as well.

::

    sage: list(range(3, 23, 2))
    [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    sage: [3,5..21]
    [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]

.. note::
   It is also important to emphasize that the ``range`` command does
   *not* include its last value!  For a quick quiz, confirm this in the
   examples above.

This all works well.  However, after a short time this will seem tedious
as well (you may have to trust us on this).  It turns out that there is
a very powerful way to create such lists in a way that very strongly
resembles the so\-called set builder notation, called a *list
comprehension* .

We start with a relatively easy example:

.. MATH::

    \{n^2\mid n\in\ZZ, 3 \leq n \leq 12\}

Who hasn't written something like this at some point in a course?

This is a natural for the list comprehension, and can be very powerful
when used in Sage.

::

    sage: [n^2 for n in [3..12]]
    [9, 16, 25, 36, 49, 64, 81, 100, 121, 144]

That's it.  This sort of turns the loop around.

- The notation is easiest if you think of it mathematically; "The set of
  :math:`n^2`, for (all) :math:`n` in the range between 3 and 13."

This is phenomenally useful.  Here is a nice plotting example::

    sage: plot([x^n for n in [2..6]],(x,0,1))
    Graphics object consisting of 5 graphics primitives

Now we apply it to the example we were doing in the first place.  Notice
we now have a nice concise description of all determinants of these
matrices, without the syntax of colon and indentation::

    sage: [det(A^i) for i in [0..4]]
    [1, -2, 4, -8, 16]

.. _Tables:

Tables
~~~~~~

Finally, getting away from strictly programming, here is a useful tip.

Some of you may be familiar with a way to take such data and put it in
tabular form from other programs. The ``table`` command does this
for us::

    sage: table( [ (i,det(A^i)) for i in [0..4] ] )
      0   1
      1   -2
      2   4
      3   -8
      4   16

Notice that each element of *this* list is two items in parentheses (a
so\-called *tuple*).

Even better, we can put a header line on it to make it really clear what
we are doing, by adding lists. We've seen keywords like ``header=True``
when doing some of our plotting and limits. What do you think will
happen if you put dollar signs around the labels in the header? ::

    sage: table( [('i', 'det(A^i)')] + [ (i,det(A^i)) for i in [0..4] ], header_row=True)
        i   det(A^i)
      ├───┼──────────┤
        0   1
        1   -2
        2   4
        3   -8
        4   16


.. _Defs:

Defining Functions
------------------

.. rubric:: Or, Extending Sage

It is often the case that Sage can do something, but doesn't have a
simple command for it.  For instance, you might want to take a matrix
and output the square of that matrix minus the original matrix.

::

    sage: A = matrix([[1,2],[3,4]])
    sage: A^2-A
    [ 6  8]
    [12 18]

How might one do this for other matrices?  Of course, you could just
always do :math:`A^2-A` again and again.  But this would be tedious and
hard to follow, as with so many things that motivate a little
programming.  Here is how we solve this problem.

::

    sage: def square_and_subtract(mymatrix):
    ....:     return mymatrix^2-mymatrix

The ``def`` command has created a new function called
``square_and_subtract``.  It should even be available using
tab\-completion.

Here are things to note about its construction:

- The input is inside the parentheses.

- The indentation and colon are crucial, as above.

- There will usually be a return value, given by ``return``.  This is
  what Sage will give below the input cell.

::

    sage: square_and_subtract(A)
    [ 6  8]
    [12 18]

::

    sage: square_and_subtract(matrix([[1.5,0],[0,2]]))
    [0.750000000000000 0.000000000000000]
    [0.000000000000000  2.00000000000000]

We can get a documentation string available by putting it in triple
quotes ``"""``.

.. skip

::

    sage: def square_and_subtract(mymatrix):
    ....:     """
    ....:     Return `A^2-A`
    ....:     """
    ....:     return mymatrix^2-mymatrix

.. skip

::

    sage: square_and_subtract?

Pretty cool!  And potentially quite helpful to students - and you -
especially if the function is complicated.  The :math:`A` typesets
properly because we put it in backticks (see above).

.. tip:
   For the *real* experts, one can use "raw strings" to include
   backslashes (say, for LaTeX) in these documentation strings, like
   ``r"""\frac{a}{b}"""``.

A very careful reader *may* have noticed that there is nothing that
requires the input ``mymatrix`` to be a matrix.  Sage will just try to
square whatever you give it and subtract the original thing.

::

    sage: square_and_subtract(sqrt(5))
    -sqrt(5) + 5

This is a typical thing to watch out for; just because you define
something doesn't mean it's useful (though in this case it was).

Try to define a function which inputs a matrix and returns the
determinant of the cube of the matrix.  (There are a few ways to do
this, of course!)

.. _Gotchas:

Gotchas from names and copies
-----------------------------

.. rubric:: Or, What's in a Name

Before we finish the tutorial, we want to point out a few
programming\-related things that often trip people up.

The first 'gotcha' is that it's possible to clobber constants!

::

    sage: i
    4

Can you figure out why ``i=4``?  Look carefully above to see when this
happened.

- This gives a valuable lesson; *any time* you use a name there is
  potential for renaming.

This may seem quite bad, but could be quite logical to do \- for
instance, if you are only dealing with real matrices. It is definitely
is something a Sage user needs to know, though.

Luckily, it's possible to restore symbolic constants.

::

    sage: reset('i')
    sage: i; i^2
    I
    -1

::

    sage: type(e)
    <class 'sage.symbolic.expression.E'>

::

    sage: type(pi)
    <class 'sage.symbolic.expression.Expression'>

Variables are another thing to keep in mind.  As mentioned briefly in
earlier tutorials, in order to maintain maximum flexibility while not
allowing things to happen which shouldn't, only ``x`` is predefined,
nothing else.

::

    sage: type(x)
    <class 'sage.symbolic.expression.Expression'>

::

    sage: type(y)
    Traceback (most recent call last):
    ...
    NameError: name 'y' is not defined

.. warning::
   There *is* a way to get around this, but it unleashes a horde of
   potential misuse. See the cells below if you are interested in this.

   ::

       sage: automatic_names(True) # not tested
       sage: trig_expand((2*x + 4*y + sin(2*theta))^2) # not tested
       4*(sin(theta)*cos(theta) + x + 2*y)^2

   This only works in the notebook.  Now we'll turn it off.

   ::

       sage: automatic_names(False) # not tested

Another related issue is that a few names are "reserved" by Python/Sage,
and which aren't allowed as variable names.

It's not surprising that 'for' is not allowed, but neither is 'lambda'
(:math:`\lambda`)!  People often request a workaround for that.

::

    sage: var('lambda')
    Traceback (most recent call last):
    ...
    ValueError: The name "lambda" is not a valid Python identifier.

There are lots of ways to get around this.  One popular, though
annoying, way is this.

::

    sage: var('lambda_')
    lambda_

::

    sage: lambda_^2-1
    lambda_^2 - 1

Still, in this one case, showing the expression still shows the Greek letter.

.. skip

::

    sage: show(lambda_^2-1)

.. MATH::

    \lambda^{2} - 1

Finally, there is another thing that can happen if you rename things too
loosely.

::

    sage: A = matrix(QQ,[[1,2],[3,4]])
    sage: B = A
    sage: C = copy(A)

This actually has just made B and A refer to the same matrix.  B isn't
like A, it *is* A. The ``copy`` command gets around this (though not
always).

::

    sage: A[0,0]=987

.. skip

::

    sage: show([A,B,C])

.. MATH::

    \left[\left(\begin{array}{rr}
    987 & 2 \\
    3 & 4
    \end{array}\right), \left(\begin{array}{rr}
    987 & 2 \\
    3 & 4
    \end{array}\right), \left(\begin{array}{rr}
    1 & 2 \\
    3 & 4
    \end{array}\right)\right]

This is very subtle if you've never programmed before.  Suffice it to
say that it is safest to let each ``=`` sign stand for one thing, and to
avoid redundant equals (unlike students at times).

.. _Advanced:

Appendix: Advanced Introductory Topics
---------------------------------------

There are several things which are useful to know about, but which are
not always introduced immediately in programming.  We give a few
examples here, but they are mainly here to make sure you have seen them
so that they are not completely surprising when they come up again.

We saw the "block" structure of Python earlier, with the indentation.
This gives the opportunity to introduce conditional statements and
comparisons.  Here, we just give an example for those who have seen
conditionals ("if" clauses) before.

::

    sage: B = matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0]])
    sage: for i in range(5): # all integers from 0 to 4, remember
    ....:     if B^i==0: # We ask if the power is the zero matrix
    ....:         print(i)
    4

We use the double equals sign to test for equality, because ``=``
assigns something to a variable name.  Notice again that colons and
indentation are the primary way that Sage/Python indicate syntax, just
as commas and spaces do in English.

Another useful concept is that of a *dictionary* .  This can be thought
of as a mathematical mapping from "keys" to "values".  The order is
*not* important and *not* guaranteed.  A dictionary is delimited by
curly brackets and correspondence is indicated by colons.

Again, we will just give a small example to illustrate the idea.

What if one wants to specify a matrix using just the nonzero entries?  A
dictionary is a great way to do this.

This one puts 3 as an entry in the :math:`(2,3)` spot, for example
(remember, this is the *third* row and *fourth* column, since we start
with zero).

::

    sage: D = {(2,3):3, (4,5):6, (6,0):-3}
    sage: C = matrix(D)
    sage: C
    [ 0  0  0  0  0  0]
    [ 0  0  0  0  0  0]
    [ 0  0  0  3  0  0]
    [ 0  0  0  0  0  0]
    [ 0  0  0  0  0  6]
    [ 0  0  0  0  0  0]
    [-3  0  0  0  0  0]

That was a lot easier than inputting the whole matrix!

Finally, although Sage tries to anticipate what you want, sometimes it
does matter how you define a given element in Sage.

- We saw this above with matrices over the rationals versus integers,
  for instance.

Here's an example with straight\-up numbers.

::

    sage: a = 2
    sage: b = 2/1
    sage: c = 2.0
    sage: d = 2 + 0*I
    sage: e = 2.0 + 0.0*I

We will not go in great depth about this, either, but it is worth
knowing about.  Notice that each of these types of numbers has or does
not have :math:`I=\sqrt{-1}`, decimal points, or division.

::

    sage: parent(a)
    Integer Ring
    sage: parent(b)
    Rational Field
    sage: parent(c)
    Real Field with 53 bits of precision
    sage: parent(d)
    Number Field in I with defining polynomial x^2 + 1 with I = 1*I
    sage: parent(e)
    Complex Field with 53 bits of precision

