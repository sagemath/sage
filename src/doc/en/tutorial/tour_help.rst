.. _chapter-help:

Getting Help
============

Sage has extensive built-in documentation, accessible by typing the
name of a function or a constant (for example), followed by a
question mark:

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring:

        The tangent function

        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring:

        The natural logarithm of the real number 2.

        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 32-bit
            0.69314718055994530941723212145817656807   # 64-bit
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <... 'function'>
    Definition:  sudoku(A)
    Docstring:

        Solve the 9x9 Sudoku puzzle defined by the matrix A.

        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

Sage also provides 'Tab completion': type the first few letters of a
function and then hit the :kbd:`Tab` key. For example, if you type
``ta`` followed by :kbd:`Tab`, Sage will print ``tachyon, tan, tanh,
taylor``. This provides a good way to find the names of functions and
other structures in Sage.


.. _section-functions:

Functions, Indentation, and Counting
====================================

To define a new function in Sage, use the ``def`` command and a colon
after the list of variable names. For example:

::

    sage: def is_even(n):
    ....:     return n%2 == 0
    sage: is_even(2)
    True
    sage: is_even(3)
    False

Note: Depending on which version of the tutorial you are viewing, you
may see three dots ``....:`` on the second line of this example.  Do
not type them; they are just to emphasize that the code is indented.
Whenever this is the case, press [Return/Enter] once at the end of the block to
insert a blank line and conclude the function definition.

You do not specify the types of any of the input arguments. You can
specify multiple inputs, each of which may have an optional default
value. For example, the function below defaults to ``divisor=2`` if
``divisor`` is not specified.

::

    sage: def is_divisible_by(number, divisor=2):
    ....:     return number%divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False

You can also explicitly specify one or either of the inputs when
calling the function; if you specify the inputs explicitly, you can
give them in any order:

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True

In Python, blocks of code are not indicated by curly braces or
begin and end blocks as in many other languages. Instead, blocks of
code are indicated by indentation, which must match up exactly. For
example, the following is a syntax error because the ``return``
statement is not indented the same amount as the other lines above
it.

.. skip

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:    return v
    Syntax Error:
           return v

If you fix the indentation, the function works:

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:     return v
    sage: even(10)
    [4, 6, 8]

Semicolons are not needed at the ends of lines; a line is in most
cases ended by a newline. However, you can put multiple statements
on one line, separated by semicolons:

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

If you would like a single line of code to span multiple lines, use
a terminating backslash:

::

    sage: 2 + \
    ....:    3
    5

In Sage, you count by iterating over a range of integers. For example,
the first line below is exactly like ``for(i=0; i<3; i++)`` in C++ or
Java:

::

    sage: for i in range(3):
    ....:     print(i)
    0
    1
    2

The first line below is like ``for(i=2;i<5;i++)``.

::

    sage: for i in range(2,5):
    ....:     print(i)
    2
    3
    4

The third argument controls the step, so the following is like
``for(i=1;i<6;i+=2)``.

::

    sage: for i in range(1,6,2):
    ....:     print(i)
    1
    3
    5

Often you will want to create a nice table to display numbers you
have computed using Sage. One easy way to do this is to use string
formatting. Below, we create three columns each of width exactly 6
and make a table of squares and cubes.

::

    sage: for i in range(5):
    ....:     print('%6s %6s %6s' % (i, i^2, i^3))
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64

The most basic data structure in Sage is the list, which is -- as
the name suggests -- just a list of arbitrary objects. For example,
using ``range``, the following command creates a list::

    sage: list(range(2,10))
    [2, 3, 4, 5, 6, 7, 8, 9]

Here is a more complicated list:

::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

List indexing is 0-based, as in many programming languages.

.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)

Use ``len(v)`` to get the length of ``v``, use ``v.append(obj)`` to
append a new object to the end of ``v``, and use ``del v[i]`` to delete
the :math:`i^{th}` entry of ``v``:

.. link

::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]

Another important data structure is the dictionary (or associative
array). This works like a list, except that it can be indexed with
almost any object (the indices must be immutable):

::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi

You can also define new data types using classes. Encapsulating
mathematical objects with classes is a powerful technique that can
help to simplify and organize your Sage programs. Below, we define a
class that represents the list of even positive integers up to *n*;
it derives from the builtin type ``list``.

::

    sage: class Evens(list):
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:         list.__init__(self, range(2, n+1, 2))
    ....:     def __repr__(self):
    ....:         return "Even positive numbers up to n."

The ``__init__`` method is called to initialize the object when
it is created; the ``__repr__`` method prints the object out. We
call the list constructor method in the second line of the
``__init__`` method. We create an object of class ``Evens`` as
follows:

.. link

::

    sage: e = Evens(10)
    sage: e
    Even positive numbers up to n.

Note that ``e`` prints using the ``__repr__`` method that we
defined. To see the underlying list of numbers, use the ``list``
function:

.. link

::

    sage: list(e)
    [2, 4, 6, 8, 10]

We can also access the ``n`` attribute or treat ``e`` like a list.

.. link

::

    sage: e.n
    10
    sage: e[2]
    6
