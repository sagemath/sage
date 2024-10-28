Basic Algebra and Calculus
==========================

Sage can perform various computations related to basic algebra and
calculus: for example, finding solutions to equations,
differentiation, integration, and Laplace transforms. See the
`Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_
documentation for more examples.

In all these examples, it is important to note that the variables in the
functions are defined to be ``var(...)``. As an example:

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

If you get a :class:`NameError`, check to see if you misspelled something,
or forgot to define a variable with ``var(...)``.


Solving Equations
-----------------

Solving Equations Exactly
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``solve`` function solves equations. To use it, first specify
some variables; then the arguments to ``solve`` are an equation (or a
system of equations), together with the variables for which to
solve:

::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

You can solve equations for one variable in terms of others:

::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]

You can also solve for several variables:

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]

The following example of using Sage to solve a system of non-linear
equations was provided by Jason Grout: first, we solve the system
symbolically:

::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3], [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]

For numerical approximations of the solutions, you can instead use:

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]

(The function ``n`` prints a numerical approximation, and the
argument is the number of bits of precision.)

Solving Equations Numerically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often times, ``solve`` will not be able to find an exact solution to
the equation or equations specified.  When it fails, you can use
``find_root`` to find a numerical solution.  For example, solve does
not return anything interesting for the following equation::

    sage: theta = var('theta')
    sage: solve(cos(theta)==sin(theta), theta)
    [sin(theta) == cos(theta)]

On the other hand, we can use ``find_root`` to find a solution to the
above equation in the range :math:`0 < \phi < \pi/2`::

    sage: phi = var('phi')
    sage: find_root(cos(phi)==sin(phi),0,pi/2)
    0.785398163397448...

Differentiation, Integration, etc.
----------------------------------

Sage knows how to differentiate and integrate many functions. For
example, to differentiate :math:`\sin(u)` with respect to :math:`u`,
do the following:

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

To compute the fourth derivative of :math:`\sin(x^2)`:

::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)

To compute the partial derivatives of :math:`x^2+17y^2` with
respect to `x` and `y`, respectively:

::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)
    34*y

We move on to integrals, both indefinite and definite. To compute
:math:`\int x\sin(x^2)\, dx` and
:math:`\int_0^1 \frac{x}{x^2+1}\, dx`

::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

To compute the partial fraction decomposition of
:math:`\frac{1}{x^2-1}`:

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)

.. _section-systems:

Solving Differential Equations
------------------------------

You can use Sage to investigate ordinary differential equations. To
solve the equation :math:`x'+x-1=0`:

::

    sage: t = var('t')    # define a variable t
    sage: x = function('x')(t)   # define x to be a function of that variable
    sage: DE = diff(x, t) + x - 1
    sage: desolve(DE, [x,t])
    (_C + e^t)*e^(-t)

This uses Sage's interface to Maxima [Max]_, and so its output may be
a bit different from other Sage output. In this case, this says
that the general solution to the differential equation is
:math:`x(t) = e^{-t}(e^{t}+c)`.

You can compute Laplace transforms also; the Laplace transform of
:math:`t^2e^t -\sin(t)` is computed as follows:

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3

Here is a more involved example. The displacement from equilibrium
(respectively) for a coupled spring attached to a wall on the left

.. CODE-BLOCK:: text

    |------\/\/\/\/\---|mass1|----\/\/\/\/\/----|mass2|
             spring1               spring2

is modeled by the system of 2nd order differential equations

.. math::

    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0

    m_2 x_2''+ k_2 (x_2-x_1) = 0,


where :math:`m_{i}` is the mass of object *i*, :math:`x_{i}` is
the displacement from equilibrium of mass *i*, and :math:`k_{i}`
is the spring constant for spring *i*.

**Example:** Use Sage to solve the above problem with
:math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`,
:math:`k_{2}=2`, :math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`,
:math:`x_{2}(0)=3`, :math:`x_{2}'(0)=0`.

Solution: Take the Laplace transform of the first equation (with
the notation :math:`x=x_{1}`, :math:`y=x_{2}`):

::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

This is hard to read, but it says that

.. math:: -2x'(0) + 2s^2 \cdot X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0


(where the Laplace transform of a lower case function like
:math:`x(t)` is the upper case function :math:`X(s)`). Take the
Laplace transform of the second equation:

::

    sage: de2 = maxima("diff(y(t),t, 2) + 2*y(t) - 2*x(t)")
    sage: lde2 = de2.laplace("t","s"); lde2.sage()
    s^2*laplace(y(t), t, s) - s*y(0) - 2*laplace(x(t), t, s) + 2*laplace(y(t), t, s) - D[0](y)(0)

This says

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.


Plug in the initial conditions for :math:`x(0)`, :math:`x'(0)`,
:math:`y(0)`, and :math:`y'(0)`, and solve the resulting two
equations:

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s]
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
      Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

Now take inverse Laplace transforms to get the answer:

::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)

Therefore, the solution is

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).


This can be plotted parametrically using

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),
    ....:     (t, 0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

The individual components can be plotted using

::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), (t,0, 2*pi), rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), (t,0, 2*pi), rgbcolor=hue(0.6))
    sage: show(p1 + p2)

For more on plotting, see :ref:`section-plot`. See section 5.5 of
[NagleEtAl2004]_ for further information on differential equations.


Euler's Method for Systems of Differential Equations
----------------------------------------------------

In the next example, we will illustrate Euler's method for first
and second order ODEs. We first recall the basic idea for first
order equations. Given an initial value problem of the form

.. math::

    y'=f(x,y), \quad y(a)=c,

we want to find the approximate value of the solution at
:math:`x=b` with :math:`b>a`.

Recall from the definition of the derivative that

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},


where :math:`h>0` is given and small. This and the DE together
give :math:`f(x,y(x))\approx
\frac{y(x+h)-y(x)}{h}`. Now solve
for :math:`y(x+h)`:

.. math::   y(x+h) \approx y(x) + h\cdot f(x,y(x)).


If we call :math:`h \cdot f(x,y(x))` the "correction term" (for lack of
anything better), call :math:`y(x)` the "old value of `y`", and
call :math:`y(x+h)` the "new value of `y`", then this
approximation can be re-expressed as

.. math::   y_{new} \approx y_{old} + h\cdot f(x,y_{old}).


If we break the interval from `a` to `b` into `n` steps, so that
:math:`h=\frac{b-a}{n}`, then we can record the information for
this method in a table.

============== =======================   =====================
:math:`x`      :math:`y`                 :math:`h\cdot f(x,y)`
============== =======================   =====================
:math:`a`      :math:`c`                 :math:`h\cdot f(a,c)`
:math:`a+h`    :math:`c+h\cdot f(a,c)`         ...
:math:`a+2h`   ...
...
:math:`b=a+nh` ???                             ...
============== =======================   =====================


The goal is to fill out all the blanks of the table, one row at a
time, until we reach the ??? entry, which is the
Euler's method approximation for  :math:`y(b)`.

The idea for systems of ODEs is similar.

**Example:** Numerically approximate :math:`z(t)` at :math:`t=1` using 4
steps of Euler's method, where :math:`z''+tz'+z=0`,
:math:`z(0)=1`, :math:`z'(0)=0`.

We must reduce the 2nd order ODE down to a system of two first
order DEs (using :math:`x=z`, :math:`y=z'`) and apply Euler's
method:

::

    sage: t,x,y = PolynomialRing(RealField(10),3,"txy").gens()
    sage: f = y; g = -x - y * t
    sage: eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
          t                x            h*f(t,x,y)                y       h*g(t,x,y)
          0                1                  0.00                0           -0.25
        1/4              1.0                -0.062            -0.25           -0.23
        1/2             0.94                 -0.12            -0.48           -0.17
        3/4             0.82                 -0.16            -0.66          -0.081
          1             0.65                 -0.18            -0.74           0.022

Therefore, :math:`z(1)\approx 0.65`.

We can also plot the points :math:`(x,y)` to get an approximate
picture of the curve. The function ``eulers_method_2x2_plot`` will
do this; in order to use it, we need to define functions `f` and
`g` which takes one argument with three coordinates: (`t`, `x`,
`y`).

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

At this point, ``P`` is storing two plots: ``P[0]``, the plot of `x`
vs. `t`, and ``P[1]``, the plot of `y` vs. `t`. We can plot both of
these as follows:

.. link

::

    sage: show(P[0] + P[1])

(For more on plotting, see :ref:`section-plot`.)

Special functions
-----------------

Several orthogonal polynomials and special functions are
implemented, using both PARI [GAP]_ and Maxima [Max]_. These are
documented in the appropriate sections ("Orthogonal polynomials"
and "Special functions", respectively) of the Sage reference
manual.

::

    sage: x = polygen(QQ, 'x')
    sage: chebyshev_U(2,x)
    4*x^2 - 1
    sage: bessel_I(1,1).n(250)
    0.56515910399248502720769602760986330732889962162109200948029448947925564096
    sage: bessel_I(1,1).n()
    0.565159103992485
    sage: bessel_I(2,1.1).n()
    0.167089499251049

At this point, Sage has only wrapped these functions for numerical use.
For symbolic use, please use the Maxima interface directly, as in
the following example:

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'


Vector calculus
---------------

See the
`Vector Calculus Tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/vector_calculus.html>`__.
