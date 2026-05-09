.. sage-doctest: needs sage.plot sage.symbolic

.. linkall

.. _prep-quickstart-differential-equations:

Sage Quickstart for Differential Equations
==========================================

This `Sage <http://www.sagemath.org>`_ quickstart tutorial was developed
for the MAA PREP Workshop "Sage: Using Open\-Source Mathematics Software
with Undergraduates" (funding provided by NSF DUE 0817071).  It is
licensed under the Creative Commons Attribution\-ShareAlike 3.0 license
(`CC BY\-SA <http://creativecommons.org/licenses/by-sa/3.0/>`_).

Solving differential equations is a combination of exact and numerical
methods, and hence a great place to explore with the computer.  We have
already seen one example of this in the calculus tutorial, which is
worth reviewing.

Basic Symbolic Techniques
-------------------------

::

    sage: y = function('y')(x)
    sage: de = diff(y,x) + y -2
    sage: h = desolve(de, y)

Forgetting about plotting for the moment, notice that there are three
things one needs to solve a differential equation symbolically:

- an abstract function ``y``;

- a differential equation, which here we put in a separate line;

- the actual **d**\ ifferential **e**\ quation **solve** command (bold for
  the acronym ``desolve`` ).

.. skip

::

    sage: show(expand(h))

.. MATH::

    c e^{\left(-x\right)} + 2

Since we did not specify any initial conditions, Sage (from Maxima) puts
in a parameter.  If we want to put in an initial condition, we use
``ics`` (for **i**\ nitial **c**\ ondition\ **s**). For example, we
set ``ics=[0,3]`` to specify that when `x=0`, `y=3`.

::

    sage: h = desolve(de, y, ics=[0,3]); h
    (2*e^x + 1)*e^(-x)

And of course we have already noted that we can plot all this with a
slope field.

::

    sage: y = var('y') # Needed so we can plot
    sage: Plot1=plot_slope_field(2-y,(x,0,3),(y,0,5))
    sage: Plot2=plot(h,x,0,3)
    sage: Plot1+Plot2
    Graphics object consisting of 2 graphics primitives

.. note::
   Regarding symbolic functions versus symbolic variables:

   - If you wanted to make ``y`` an abstract function again instead of a
     variable, you'd have to do that separately.  A differential equation
     requires a ``function`` but plotting requires a ``var``.

   - Another option is to let ``z`` be the name of the vertical axis
     variable.

   - Either way something will have to give, since in common speaking
     about these things we treat ``y`` as both a variable and a function,
     which is much trickier to accomplish with a computer.

There are many other differential equation facilities in Sage.  We can't
cover all the variants of this in a quickstart, but the documentation is
good for symbolic solvers.

.. skip

::

    sage: desolvers?

For instance, Maxima can do systems, as well as use Laplace transforms,
and we include versions of these wrapped for ease of use.

In all differential equation solving routines, it is important to pay
attention to the syntax!  In the following example, we have placed the
differential equation in the body of the command, and had to specify
that ``f`` was the **d**\ ependent **var**\ iable (``dvar``), as well as
give initial conditions :math:`f(0)=1` and :math:`f'(0)=2`, which
gives the last list in the example.

::

    sage: f=function('f')(x)
    sage: desolve_laplace(diff(f,x,2) == 2*diff(f,x)-f, dvar = f, ics = [0,1,2])
    x*e^x + e^x

::

    sage: g(x)=x*e^x+e^x
    sage: derivative(g,x,2)-2*derivative(g,x)+g
    x |--> 0

Numerical and Power Series Methods
----------------------------------

There are also numerical methods.

For instance, one of the options above was ``desolve_rk4``.  This is a
fourth\-order Runge\-Kutta method, and returns appropriate (numerical)
output.  Here, we *must* give the dependent variable *and* initial
conditions.

::

    sage: y = function('y')(x)
    sage: de = diff(y,x) + y -2
    sage: h = desolve_rk4(de, y, step=.05, ics=[0,3])

It can be fun to compare this with the original, symbolic solution.  We
use the ``points`` command from the advanced plotting tutorial.

::

    sage: h1 = desolve(de, y, ics=[0,3])
    sage: plot(h1,(x,0,5),color='red')+points(h)
    Graphics object consisting of 2 graphics primitives

The primary use of numerical routines from here is pedagogical in nature.

For more advanced numerical routines, we primarily use the GNU
scientific library.  Using this is a little more sophisticated, but
gives a wealth of options.

.. skip

::

    sage: ode_solver?

We can even do power series solutions.  In order to do this, we must
first define a special *power series ring* , including the precision.

::

    sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)
    sage: a = -1 + 0*t
    sage: b = 2 + 0*t
    sage: h=a.solve_linear_de(b=b,f0=3,prec=10)

This power series solution is pretty good for a while!

::

    sage: h = h.polynomial()
    sage: plot(h,-2,5)+plot(2+e^-x,(x,-2,5),color='red',linestyle=':',thickness=3)
    Graphics object consisting of 2 graphics primitives

This was just an introduction; there are a lot of resources for
differential equations using Sage elsewhere, including a book by David
Joyner, who wrote much of the original code wrapping Maxima for Sage to
do just this.

