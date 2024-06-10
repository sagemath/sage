
Assignment, Equality, and Arithmetic
====================================

With some minor exceptions, Sage uses the Python programming language,
so most introductory books on Python will help you to learn Sage.

Sage uses ``=`` for assignment. It uses ``==``, ``<=``, ``>=``, ``<`` and ``>`` for
comparison:

::

    sage: a = 5
    sage: a
    5
    sage: 2 == 2
    True
    sage: 2 == 3
    False
    sage: 2 < 3
    True
    sage: a == 5
    True

Sage provides all of the basic mathematical operations:

::

    sage: 2**3    #  ** means exponent
    8
    sage: 2^3     #  ^ is a synonym for ** (unlike in Python)
    8
    sage: 10 % 3  #  for integer arguments, % means mod, i.e., remainder
    1
    sage: 10/4
    5/2
    sage: 10//4   #  for integer arguments, // returns the integer quotient
    2
    sage: 4 * (10 // 4) + 10 % 4 == 10
    True
    sage: 3^2*4 + 2%5
    38

The computation of an expression like ``3^2*4 + 2%5`` depends on
the order in which the operations are applied; this is specified in
the "operator precedence table" in :ref:`section-precedence`.

Sage also provides many familiar mathematical functions; here are
just a few examples:

::

    sage: sqrt(3.4)
    1.84390889145858
    sage: sin(5.135)
    -0.912021158525540
    sage: sin(pi/3)
    1/2*sqrt(3)

As the last example shows, some mathematical expressions return
'exact' values, rather than numerical approximations. To get a
numerical approximation, use either the function ``N`` or the method
``n`` (and both of these have a longer name, ``numerical_approx``, and
the function ``N`` is the same as ``n``)). These take optional
arguments ``prec``, which is the requested number of bits of
precision, and ``digits``, which is the requested number of decimal
digits of precision; the default is 53 bits of precision.

::

    sage: exp(2)
    e^2
    sage: n(exp(2))
    7.38905609893065
    sage: sqrt(pi).numerical_approx()
    1.77245385090552
    sage: sin(10).n(digits=5)
    -0.54402
    sage: N(sin(10),digits=10)
    -0.5440211109
    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

Python is dynamically typed, so the value referred to by each
variable has a type associated with it, but a given variable may
hold values of any Python type within a given scope:

::

    sage: a = 5   # a is an integer
    sage: type(a)
    <class 'sage.rings.integer.Integer'>
    sage: a = 5/3  # now a is a rational number
    sage: type(a)
    <class 'sage.rings.rational.Rational'>
    sage: a = 'hello'  # now a is a string
    sage: type(a)
    <... 'str'>

The C programming language, which is statically typed, is much
different; a variable declared to hold an int can only hold an int
in its scope.
