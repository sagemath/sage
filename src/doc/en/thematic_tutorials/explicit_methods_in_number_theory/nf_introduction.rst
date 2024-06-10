.. section-nf:

Introduction to Number Fields
=============================

In Sage, we can create the number field
:math:`\QQ(\sqrt[3]{2})` as follows.

::

    sage: K.<alpha> = NumberField(x^3 - 2)

The above creates *two* Sage objects, :math:`K` and
:math:`\alpha`. Here :math:`K` "is" (isomorphic to) the number
field :math:`\QQ(\sqrt[3]{2})`, as we confirm below:

.. link

::

    sage: K
    Number Field in alpha with defining polynomial x^3 - 2

and :math:`\alpha` is a root of :math:`x^3 - 2`, so
:math:`\alpha` is an abstract choice of :math:`\sqrt[3]{2}` (no
specific embedding of the number field :math:`K` into
:math:`\CC` is chosen by default in Sage-3.1.2):

.. link

::

    sage: alpha^3
    2
    sage: (alpha+1)^3
    3*alpha^2 + 3*alpha + 3


The variable :math:`x`
----------------------

Note that we did *not* define :math:`x` above before using it.
You could "break" the above example by redefining :math:`x` to be
something funny::

    sage: x = 1
    sage: K.<alpha> = NumberField(x^3 - 2)
    Traceback (most recent call last):
    ...
    TypeError: polynomial (=-1) must be a polynomial.

The *Traceback* above indicates that there was an error.
Potentially lots of detailed information about the error (a
"traceback") may be given after the word ``Traceback``
and before the last line, which contains the actual error
messages.


.. note::

   *Important*: whenever you use Sage and get a big error, look at the
   last line for the actual error, and only look at the rest if you are
   feeling adventurous.  In the notebook, the part indicated
   by ``...`` above is not displayed; to see it, click just to the left of
   the word *Traceback* and the traceback will appear.

If you redefine :math:`x` as above, but need to define a number
field using the indeterminate :math:`x`, you have several
options. You can reset :math:`x` to its default value at the
start of Sage, you can redefine :math:`x` to be a symbolic
variable, or you can define :math:`x` to be a polynomial
indeterminate (a polygen)::

    sage: reset('x')
    sage: x
    x
    sage: x = 1
    sage: x = var('x')
    sage: x
    x
    sage: x = 1
    sage: x = polygen(QQ, 'x')
    sage: x
    x
    sage: x = 1
    sage: R.<x> = PolynomialRing(QQ)
    sage: x
    x


Using tab completion to get the methods of an object
----------------------------------------------------

One you have created a number field :math:`K`, type ``K.[tab key]`` to
see a list of functions. Type, e.g., ``K.minkowski_embedding?[tab
key]`` to see help on the ``minkowski_embedding`` command. To see
source code, type ``K.minkowski_embedding??[tab key]``.

.. skip

::

    sage: K.<alpha> = NumberField(x^3 - 2)
    sage: K.[tab key]


.. _section-symbolic:

Symbolic Expressions
--------------------

Another natural way for us to create certain number fields is to
create a symbolic expression and adjoin it to the rational
numbers. Unlike Pari and Magma (and like Mathematica and Maple), Sage
also supports manipulation of symbolic expressions and solving
equations, without defining abstract structures such as a number
fields. For example, we can define a variable :math:`a=\sqrt{2}` as an
abstract symbolic object by simply typing ``a = sqrt(2)``. When we
type ``parent(a)`` below, Sage tells us the mathematical object that
it views :math:`a` as being an element of; in this case, it's the ring
of all symbolic expressions.

::

    sage: a = sqrt(2)
    sage: parent(a)
    Symbolic Ring


sqrt(2) in Pari and Magma
-------------------------

In particular, typing ``sqrt(2)`` does *not* numerically extract an
approximation to :math:`\sqrt{2}`, like it would in Pari or Magma. We
illustrate this below by calling Pari (via the gp interpreter) and
Magma directly from within Sage. After we evaluate the following two
input lines, copies of GP/Pari and Magma are running, and there is a
persistent connection between Sage and those sessions.

::

    sage: gp('sqrt(2)')
    1.414213562373095048801688724...
    sage: magma('Sqrt(2)')                  # optional - magma
    1.414213562373095048801688724...

You probably noticed a pause when evaluated the second line as
Magma started up. Also, note the ``# optional``
comment, which indicates that the line won't work if you don't have
Magma installed.


Numerically evaluating sqrt(2)
------------------------------

Incidentally, if you want to numerically evaluate :math:`\sqrt{2}` in
Sage, just give the optional ``prec`` argument to the ``sqrt``
function, which takes the required number of *bits* (binary digits)
of precision.

::

    sage: sqrt(2, prec=100)
    1.4142135623730950488016887242

It's important to note in computations like this that there is not an
*a priori* guarantee that ``prec`` bits of the *answer* are all
correct. Instead, what happens is that Sage creates the number
:math:`2` as a floating point number with :math:`100` bits of
accuracy, then asks Paul Zimmerman's MPFR C library to compute the
square root of that approximate number.


Arithmetic with sqrt(2)
-----------------------

We return now to our symbolic expression :math:`a = \sqrt{2}`. If
you ask to square :math:`a+1` you simply get the formal square.
To expand out this formal square, we use the expand command.

::

    sage: a = sqrt(2)
    sage: (a+1)^2
    (sqrt(2) + 1)^2
    sage: expand((a+1)^2)
    2*sqrt(2) + 3


Adjoining a symbolic expression
-------------------------------

Given any symbolic expression for which Sage can compute its
minimal polynomial, you can construct the number field obtained by
adjoining that expression to :math:`\QQ`. The notation is
quite simple - just type ``QQ[a]`` where ``a`` is the symbolic expression.

::

    sage: a = sqrt(2)
    sage: K.<b> = QQ[a]
    sage: K
    Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?
    sage: b
    sqrt2
    sage: (b+1)^2
    2*sqrt2 + 3
    sage: QQ[a/3 + 5]
    Number Field in a with defining polynomial x^2 - 10*x + 223/9 with a = 5.471404520791032?


Coercion: QQ[a] versus QQ(a)
----------------------------

You can't create the number field :math:`\QQ(a)` in Sage by
typing ``QQ(a)``, which has a *very different* meaning in Sage. It
means "try to create a rational number from :math:`a`." Thus ``QQ(a)``
in Sage is the analogue of ``QQ!a`` in Magma (Pari has no notion of
rings such as ``QQ``).

::

    sage: a = sqrt(2)
    sage: QQ(a)
    Traceback (most recent call last):
    ...
    TypeError: unable to convert sqrt(2) to a rational

In general, if :math:`X` is a ring, or vector space or other "parent
structure" in Sage, and :math:`a` is an element, type ``X(a)`` to
make an element of :math:`X` from :math:`a`. For example, if :math:`X`
is the finite field of order :math:`7`, and :math:`a=2/5` is a
rational number, then ``X(a)`` is the finite field element :math:`6`
(as a quick exercise, check that this is mathematically the correct
interpretation).

::

    sage: X = GF(7); a = 2/5
    sage: X(a)
    6


Solving a cubic equation
------------------------

As a slightly less trivial illustration of symbolic manipulation,
consider the cubic equation

.. cubic:

.. math::

  x^3 + \sqrt{2} x + 5 = 0.


In Sage, we can create this equation, and find an exact symbolic
solution.

::

    sage: x = var('x')
    sage: eqn =  x^3 + sqrt(2)*x + 5 == 0
    sage: a = solve(eqn, x)[0].rhs()

The first line above makes sure that the symbolic variable :math:`x`
is defined, the second creates the equation ``eqn``, and the third
line solves ``eqn`` for :math:`x`, extracts the first solution (there
are three), and takes the right hand side of that solution and assigns
it to the variable ``a``.


Viewing complicated symbolic expressions
----------------------------------------

To see the solution nicely typeset, use the ``pretty_print``
command

.. link

::

    sage: pretty_print(a)
    -1/2*(I*sqrt(3) + 1)*(1/6*sqrt(8/3*sqrt(2) + 225) - 5/2)^(1/3) + 1/6*sqrt(2)*(-I*sqrt(3) + 1)/(1/6*sqrt(8/3*sqrt(2) + 225) - 5/2)^(1/3)

.. math::

    -\frac{1}{2} \, {(i \, \sqrt{3} + 1)} {(\frac{1}{18} \, \sqrt{8 \, \sqrt{2} + 675} \sqrt{3} - \frac{5}{2})}^{\left(\frac{1}{3}\right)} + \frac{1}{6} \, \frac{{(-i \, \sqrt{3} + 1)} \sqrt{2}}{{(\frac{1}{18} \, \sqrt{8 \, \sqrt{2} + 675} \sqrt{3} - \frac{5}{2})}^{\left(\frac{1}{3}\right)}}


You can also see the latex needed to paste :math:`a` into a paper
by typing ``latex(a)``. The ``latex``
command works on most Sage objects.

.. link

::

    sage: latex(a)
    -\frac{1}{2} \, {\left(i \, \sqrt{3} + 1\right)} ...


Adjoining a root of the cubic
-----------------------------

Next, we construct the number field obtained by adjoining the solution
``a`` to :math:`\QQ`. Notice that the minimal polynomial of the
root is :math:`x^6 + 10x^3 - 2x^2 + 25`.

.. warning::

   The following tests are currently broken until :issue:`5338` is
   fixed.

.. skip

::

    sage: K.<b> = QQ[a]
    sage: K
    Number Field in a with defining
    polynomial x^6 + 10*x^3 - 2*x^2 + 25
    sage: a.minpoly()
    x^6 + 10*x^3 - 2*x^2 + 25
    sage: b.minpoly()
    x^6 + 10*x^3 - 2*x^2 + 25

We can now compute interesting invariants of the number field
:math:`K`

.. skip

::

    sage: K.class_number()
    5
    sage: K.galois_group().order()
    72


