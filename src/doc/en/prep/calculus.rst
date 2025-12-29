.. sage-doctest: needs sage.plot sage.symbolic

.. linkall

.. _prep-calculus:

Tutorial for Calculus
=====================

This `Sage <http://www.sagemath.org>`_ document is one of the tutorials
developed for the MAA PREP Workshop "Sage: Using Open\-Source
Mathematics Software with Undergraduates" (funding provided by NSF DUE
0817071).  It is licensed under the Creative Commons
Attribution\-ShareAlike 3.0 license (`CC BY\-SA
<http://creativecommons.org/licenses/by-sa/3.0/>`_).

This tutorial has the following sections.  The first three are based on
the topics encountered in a typical three\-semester calculus sequence in
the United States; the final section is a checkpoint of sorts.

- :ref:`Calc1`

- :ref:`Calc2`

- :ref:`Calc3`

- :ref:`Exam`

The tutorial assumes that one is familiar with the basics of Sage, such
as outlined in the previous tutorials.

For a refresher, make sure the syntax below for defining a function
and getting a value makes sense; then evaluate the cell by clicking
the "evaluate" link, or by pressing :kbd:`Shift` + :kbd:`Enter` (hold
down :kbd:`Shift` while pressing the :kbd:`Enter` key).

::

    sage: f(x)=x^3+1
    sage: f(2)
    9

We'll use this function several times in this tutorial, so it's
important that it is evaluated!

.. _Calc1:

Calculus 1
----------

The calculus is amazing not so much because it solves problems of
tangency or area \- individual solutions to these problems have been
known before.  It is amazing because it gives a remarkably comprehensive
set of rules for symbolic manipulation for solving such problems in
great generality!

Thus, the most typical use of a computer system like Sage in this
context is simply to help check (or make less tedious) basic symbolic
manipulation of things like derivatives.  Let's first show what our
function :math:`f` is again, then do things with it.

.. skip

::

    sage: show(f)

.. MATH::

    x \ {\mapsto}\ x^{3} + 1

Something most (but not all) curricula for first\-semester calculus do
is spend a fair amount of time with limits.

What is the limit of this function as :math:`x` approaches :math:`1`?

::

    sage: lim(f,x=1)
    x |--> 2

Sage gives the answer we expect.  The syntax for limits is pretty
straightforward, though it may differ slightly from that of other
systems.

Since a limit may exist even if the function is not defined, Sage uses
the syntax ``x=1`` to show the input value approached in all situations.

::

    sage: lim((x^2-1)/(x-1),x=1)
    2

Next, we'll return to the original function, and check the two
directional limits.

The syntax uses the extra  *keyword*  ``dir``.

- Notice that we use ``dir='right'``, not ``dir=right``.

- It's also okay to use ``dir='+'``, like we use ``dir='-'`` below.

 - This is basically because otherwise we might have done something
   like let ``right=55`` earlier, so instead Sage \- and Python \-
   requires us to put ``'right'`` and ``'-'`` in quotes to make it
   clear they are not variables.

::

    sage: lim(f,x=1,dir='-'); lim(f,x=1,dir='right'); f(1)
    x |--> 2
    x |--> 2
    2

By comparing with :math:`f(1)`, we see that :math:`f(x)` is continuous
at :math:`x=1`.

This cell also reminds us of something else:

- By separating several commands with semicolons, we can tell Sage to
  evaluate all of them while still seeing all their outputs.

What we spend the most time on in Calculus 1 is derivatives, and Sage is
fully\-featured here.

For example, here are three ways to get the basic, single-variable
derivative of :math:`f(x)=x^3+1`.

::

    sage: diff(f,x); derivative(f,x); f.derivative(x)
    x |--> 3*x^2
    x |--> 3*x^2
    x |--> 3*x^2

Naturally, Sage knows all of the derivatives you want.

::

    sage: derivative(sinh(x^2+sqrt(x-1)),x)
    1/2*(4*x + 1/sqrt(x - 1))*cosh(x^2 + sqrt(x - 1))

And maybe even knows those you do not want.  In this case, we put the
computation inside ``show()`` since the output is so long.

.. skip

::

    sage: show(derivative(sinh(x^2+sqrt(x-1)),x,3))

.. MATH::

    \frac{1}{8} \, {\left(4 \, x + \frac{1}{\sqrt{x - 1}}\right)}^{3} \cosh\left(\sqrt{x - 1} + x^{2}\right) - \frac{3}{8} \, {\left(\frac{1}{{\left(x - 1\right)}^{\frac{3}{2}}} - 8\right)} {\left(4 \, x + \frac{1}{\sqrt{x - 1}}\right)} \sinh\left(\sqrt{x - 1} + x^{2}\right) + \frac{3 \, \cosh\left(\sqrt{x - 1} + x^{2}\right)}{8 \, {\left(x - 1\right)}^{\left(\frac{5}{2}\right)}}

A common question is why Sage might not check automatically if there is
some "simpler" version.  But simplifying an expression, or indeed, even
defining what a 'simple' expression is, turns out to be a very hard
technical and mathematical problem in general.  Computers will not solve
every problem!

As a brief interlude, let's consider an application of our ability to do
some basic differential calculus.

First, as a reminder of plotting from the previous tutorial, try to plot
the function :math:`f(x)`, together with its tangent line at
:math:`x=1`, in the empty cells below.  (If you can, do it without
looking at the previous tutorial for a reminder.)

Did you get it?

Of course, in general we might want to see tangent lines at lots of
different points \- for instance, for a class demonstration.

So in the following cell, there are several auxiliary elements:

- We define the plot ``P`` of the original function.

- There is a parameter ``c=1/3``, which is the :math:`x`-value where we
  want a tangent line.

- We let ``fprime`` be the derivative function simply by declaring it
  equal to the derivative.

- We let ``L`` be the tangent line defined by the point-slope formula at
  :math:`x=c`.

- We make ``Q`` be the plot of this line.

Finally, we plot everything together in the last line by adding
:math:`P+Q`.

::

    sage: P=plot(f,(x,-1,1))
    sage: c=1/3
    sage: fprime=derivative(f,x)
    sage: L(x)=fprime(c)*(x-c)+f(c)
    sage: Q=plot(L,(x,-1,1), color="red", linestyle="--")
    sage: P+Q
    Graphics object consisting of 2 graphics primitives

You may want to experiment by

- Changing ``c`` to some other value, or

- Changing the function ``f``, or

- Changing the colors, or

- Changing something else (like the ``linestyle`` used for the tangent line).

Ideally, it would be *extremely* easy to change that parameter
:math:`c`.  In the cell below, we show our second example of a Sage
"interact" (or "Sagelet").

In this one, dragging a slider will show the tangent line moving.

- Future tutorials will explain this process in more detail; here it's
  just as an example.

- However, the reader will note how very similar the code for the two
  cells is.

.. skip

::

    sage: %auto
    sage: f(x)=x^3+1
    sage: @interact
    sage: def _(c=(1/3,(-1,1))):
    ....:     P=plot(f,(x,-1,1))
    ....:     fprime=derivative(f,x)
    ....:     L(x)=fprime(c)*(x-c)+f(c)
    ....:     Q=plot(L,(x,-1,1),color="red", linestyle="--")
    ....:     show(P+Q+point((c,f(c)), pointsize=40, color='red'),ymin=0,ymax=2)

A very sharp\-eyed reader will also have noticed that the previous cell
had ``%auto`` at the very top, and that it was not necessary to evaluate
the cell to use it.

- The command ``%auto`` allows us to have a cell, especially an
  interactive one, all loaded up as soon as we start \- particularly
  convenient for a classroom situation.

- Such instructions are called *percent directives* .  Most are
  documented in the notebook help one can access at the top of any
  worksheet.

A final topic in Calculus 1 usually is basic integration.  The syntax
for indefinite integration is similar to that for differentiation.

::

    sage: integral(cos(x),x)
    sin(x)

We do not get the whole indefinite integral, just a convenient
antiderivative.

- (If you were to get a different answer 'by hand', remember that being
  an antiderivative means the answer is correct *up to a constant* \- and
  deciding whether this is the case with two expressions has the same
  problems as the issue of "simplification" above.)

Definite integration has similar syntax to plotting.

::

    sage: integral(cos(x),(x,0,pi/2))
    1

.. _Calc2:

Calculus 2
----------

Second\-semester calculus is typically more challenging.

- One reason for that is that the computational problems are not so
  straightforward as computing derivatives and basic integrals.

- Another reason is that the second semester is usually where the harder
  versions of problems from the first semester show up.

Nonetheless, Sage can handle this as well.

Sage includes a large number of indefinite integrals (via Maxima),
though not all the ones you will find in a comprehensive table.

::

    sage: h(x)=sec(x)
    sage: h.integrate(x)
    x |--> log(sec(x) + tan(x))

Since I defined ``h`` as a function, the answer I get is also a
function.  If I just want an expression as the answer, I can do the
following.

::

    sage: integrate(sec(x),x)
    log(sec(x) + tan(x))

Here is another (longer) example.  Do you remember what command would
help it look nicer in the browser?

::

    sage: integrate(1/(1+x^5),x)
    1/5*sqrt(5)*(sqrt(5) + 1)*arctan((4*x + sqrt(5) - 1)/sqrt(2*sqrt(5) + 10))/sqrt(2*sqrt(5) + 10)
    + 1/5*sqrt(5)*(sqrt(5) - 1)*arctan((4*x - sqrt(5) - 1)/sqrt(-2*sqrt(5) + 10))/sqrt(-2*sqrt(5) + 10)
    - 1/10*(sqrt(5) + 3)*log(2*x^2 - x*(sqrt(5) + 1) + 2)/(sqrt(5) + 1)
    - 1/10*(sqrt(5) - 3)*log(2*x^2 + x*(sqrt(5) - 1) + 2)/(sqrt(5) - 1)
    + 1/5*log(x + 1)

Some integrals are a little tricky, of course. Sage tries hard to integrate using Maxima, Sympy, and (optionally, if installed) Giac. The following can only be integrated by Giac::

    sage: # needs sage.libs.giac
    sage: integral(1/(1+x^10),x)
    ...1/20*(sqrt(5) + 1)*arctan((4*x + sqrt(-2*sqrt(5) + 10))/(sqrt(5) + 1))
    + 1/20*(sqrt(5) + 1)*arctan((4*x - sqrt(-2*sqrt(5) + 10))/(sqrt(5) + 1))
    + 1/20*(sqrt(5) - 1)*arctan((4*x + sqrt(2*sqrt(5) + 10))/(sqrt(5) - 1))
    + 1/20*(sqrt(5) - 1)*arctan((4*x - sqrt(2*sqrt(5) + 10))/(sqrt(5) - 1))
    + 1/40*sqrt(2*sqrt(5) + 10)*log(x^2 + 1/2*x*sqrt(2*sqrt(5) + 10) + 1)
    - 1/40*sqrt(2*sqrt(5) + 10)*log(x^2 - 1/2*x*sqrt(2*sqrt(5) + 10) + 1)
    + 1/40*sqrt(-2*sqrt(5) + 10)*log(x^2 + 1/2*x*sqrt(-2*sqrt(5) + 10) + 1)
    - 1/40*sqrt(-2*sqrt(5) + 10)*log(x^2 - 1/2*x*sqrt(-2*sqrt(5) + 10) + 1)
    + 1/5*arctan(x)

If you ask for Maxima specifically, the result can be partial::

    sage: integral(1/(1+x^10),x, algorithm='maxima')
    1/5*arctan(x)
    - 1/5*integrate((x^6 - 2*x^4 + 3*x^2 - 4)/(x^8 - x^6 + x^4 - x^2 + 1), x)

If no antiderivative is found, the result is just the input::

    sage: result = integral(sinh(x^2+sqrt(x-1)),x)  # long time (15s on sage.math, 2012)
    ...
    sage: result  # long time
    integrate(sinh(x^2 + sqrt(x - 1)), x)

This last one stumps other systems too.

However, if there is a special function which helps compute the
integral, Sage will look for it.  In the following case there is no
elementary antiderivative, but the ``erf`` function helps us out.

::

    sage: integral(e^(-x^2),x)
    1/2*sqrt(pi)*erf(x)

Do not forget, if this function is unfamiliar to you (as it might be to
students trying this integral), Sage's contextual help system comes to
the rescue.

.. skip

::

    sage: erf?

There are several ways to do definite integrals in Sage.

The most obvious one is simply turning

.. MATH::

    \int f(x)dx

into

.. MATH::

    \int_a^b f(x)dx\; ,

as indicated in the Calculus I section.

::

    sage: integral(cos(x),(x,0,pi/2))
    1

The preferred syntax puts the variable and endpoints together in parentheses.

Just like with derivatives, we can visualize this integral using some of
the plotting options from the plotting tutorial.

::

    sage: plot(cos(x),(x,0,pi/2),fill=True,ticks=[[0,pi/4,pi/2],None],tick_formatter=pi)
    Graphics object consisting of 2 graphics primitives

It is possible to be completely symbolic in doing integration.  If you
do this, you'll have to make sure you define anything that's a symbolic
variable \- which includes constants, naturally.

::

    sage: var('a,b')
    (a, b)
    sage: integral(cos(x),(x,a,b))
    -sin(a) + sin(b)

On the numerical side, sometimes the answer one gets from the
Fundamental Theorem of Calculus is not entirely helpful.  Recall that
:math:`h` is the secant function.

::

    sage: integral(h,(x,0,pi/7))
    1/2*log(sin(1/7*pi) + 1) - 1/2*log(-sin(1/7*pi) + 1)

Here, just a number might be more helpful.  Sage has several ways of
numerical evaluating integrals.

- Doing a definite integral symbolically, then approximating it numerically

- The  ``numerical_integral``  function

- The . ``nintegrate``  method

The first one, using the n or N function for numerical approximation,
was also mentioned in the introductory tutorial.

::

    sage: N(integral(h,(x,0,pi/8)))
    0.403199719161511

The second function, ``numerical_integral``, uses a powerful numerical
program (the GNU Scientific Library).

- Unfortunately, the syntax for this function is not yet consistent with
  the rest of Sage.

- Helpfully, the output has two elements \- the answer you desire, and
  its error tolerance.

::

    sage: numerical_integral(h,0,pi/8)
    (0.4031997191615114, 4.476416117355069e-15)

To access just the number, one asks for the 'zeroth' element of this
sequence of items.  This is done with the following bracket notation.

::

    sage: numerical_integral(h,0,pi/8)[0]
    0.4031997191615114

Notice that we began counting at zero.  This is fairly typical in
computer programs (though certainly not universal).

To aid readability (more important than one might think), we often
assign the numerical integral to a variable, and then take the zeroth
element of that.

::

    sage: ni = numerical_integral(h,0,pi/8)
    sage: ni[0]
    0.4031997191615114

Finally, the ``.nintegrate()`` method from Maxima gives even more extra
information.

- Notice again the period/dot needed to use this.

- It is only possible to use ``h(x)``; doing ``h.nintegrate()`` raises an error.

::

    sage: h(x).nintegrate(x,0,pi/8)
    (0.4031997191615114, 4.47641611735507e-15, 21, 0)

Second\-semester calculus usually also covers various topics in
summation.  Sage can sum many abstract series; the notation is similar
to plotting and integration.

::

    sage: var('n') # Do not forget to declare your variables
    n
    sage: sum((1/3)^n,n,0,oo)
    3/2

This is the geometric series, of course.

The next one is the famous result that a row of Pascal's triangle is a
power of 2 \-

.. MATH::

    \binom{n}{0}+\binom{n}{1}+\binom{n}{2}+\cdots+\binom{n}{n-1}+\binom{n}{n}=2^n\; ,

which has many pleasing combinatorial interpretations.

::

    sage: k = var('k') # We already declared n, so now we just need k
    sage: sum(binomial(n,k), k, 0, n)
    2^n

Do you remember what to do to see how we typed the nice sum in the text
above? That's right, we can double\-click the text area/cell to see
this.

Sage also can compute Taylor polynomials.

Taylor expansions depend on a lot of things.  Whenever there are several
inputs, keeping syntax straight is important.  Here we have as inputs:

- the function,

- the variable,

- the point around which we are expanding the function, and

- the degree.

In the next cell, we call :math:`g(x)` the Taylor polynomial in question.

::

    sage: g(x)=taylor(log(x),x,1,6); g(x)
    -1/6*(x - 1)^6 + 1/5*(x - 1)^5 - 1/4*(x - 1)^4 + 1/3*(x - 1)^3 - 1/2*(x - 1)^2 + x - 1

Notice how close the approximation is to the function on this interval!

::

    sage: plot(g,(x,0,2))+plot(log(x),(x,0,2),color='red')
    Graphics object consisting of 2 graphics primitives

.. _Calc3:

Calculus 3
----------

We have already seen three\-dimensional plotting, so it is not
surprising that Sage has support for a variety of multivariable calculus
problems.

.. warning::
   We will often need to define all variables other than :math:`x`.

::

    sage: var('y')
    y
    sage: f(x,y)=3*sin(x)-2*cos(y)-x*y

Above, we have defined a typical function of two variables.

Below, we use the separating semicolons to demonstrate several things
one might do with such a function, including:

- The gradient vector of all :math:`\frac{\partial f}{\partial x_i}`

- The Hessian of all possible second derivatives

- A double partial derivative of :math:`f` with respect to :math:`x`,
  then :math:`y` (that is, :math:`\frac{\partial f}{\partial y\partial
  x}`)

::

    sage: f.gradient(); f.hessian(); f.diff(x,y)
    (x, y) |--> (-y + 3*cos(x), -x + 2*sin(y))
    [(x, y) |--> -3*sin(x)        (x, y) |--> -1]
    [       (x, y) |--> -1  (x, y) |--> 2*cos(y)]
    (x, y) |--> -1

In an effort to make the syntax simpler, the gradient and Hessian are
also available by asking for a total derivative.  We also ask for nicer
output again.

.. skip

::

    sage: show(f.diff()); show(f.diff(2))

.. MATH::

    \left( x, y \right) \ {\mapsto} \ \left(-y + 3 \, \cos\left(x\right),\,-x + 2 \, \sin\left(y\right)\right)

.. MATH::

    \left(\begin{array}{rr}
    \left( x, y \right) \ {\mapsto} \ -3 \, \sin\left(x\right) & \left( x, y \right) \ {\mapsto} \ -1 \\
    \left( x, y \right) \ {\mapsto} \ -1 & \left( x, y \right) \ {\mapsto} \ 2 \, \cos\left(y\right)
    \end{array}\right)

If we take the determinant of the Hessian, we get something useful for
evaluating (the two-dimensional) critical points of :math:`f`.

.. skip

::

    sage: show(f.diff(2).det())

.. MATH::

    \left( x, y \right) \ {\mapsto} \ -6 \, \sin\left(x\right) \cos\left(y\right) - 1

These ideas are particularly helpful if one wants to plot a vector field.

The following example is of the gradient.  The vector plotted in the
cell below is the unit vector in the direction :math:`(1,2)`.

::

    sage: P=plot_vector_field(f.diff(), (x,-3,3), (y,-3,3))
    sage: u=vector([1,2])
    sage: Q=plot(u/u.norm())
    sage: P+Q
    Graphics object consisting of 2 graphics primitives

Rather than actually figure out the unit vector in that direction, it's
easier to let Sage compute it by dividing the vector by its norm.

The directional derivative itself (in that direction, at the origin) can
also be computed in this way.

::

    sage: (f.diff()*u/u.norm())(0,0)
    3/5*sqrt(5)

Another useful type of plot in these situations is a contour plot.

Notice that the one below uses several options.  Try to correlate the
options with features of the graphic.

::

    sage: y = var('y')
    sage: contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),\
    ....:    contours=[-8,-4,0,4,8], colorbar=True, labels=True, label_colors='red')
    Graphics object consisting of 1 graphics primitive

In this one, we have used options to:

- Explicitly list the contours we want to show,

- Label these contours,

- Place a color bar on the side to show the different levels.

(Incidentally, the ``True`` and ``False`` valued options are some of the
few non\-numerical ones that do *not* need quotes.)

This is another good time to remind us we must explicitly ask for
:math:`y` to be a variable here, as will be the case a few more times.

As you gain experience in Sage, we will slowly explain less and less of
the syntax of commands in these tutorials.  You can think of places
where not everything is explained as a mini\-quiz.

For example, the next example shows how one currently does a multiple
integral.  What have we done here?

::

    sage: integrate(integrate(f,(x,0,pi)),(y,0,pi))
    6*pi - 1/4*pi^4

Answer: notice that ``integrate(f,(x,0,pi))`` has been itself placed as
the function inside ``integrate(...,(y,0,pi))``.

We could use a 3D plot to help visualize this; these were already
mentioned in the symbolics and plotting tutorial.

::

    sage: plot3d(f,(x,0,pi),(y,0,pi),color='red')+plot3d(0,(x,0,pi),(y,0,pi))
    Graphics3d Object

In addition to multivariate calculus, Calculus 3 often covers parametric
calculus of a single variable.  Sage can do arbitrary parametric plots,
with fairly natural syntax.

This plot shows the tangent line to the most basic Lissajous curve at
:math:`t=1`.  The commands should be strongly reminiscent of the ones at
the beginning of this tutorial.

::

    sage: t = var('t')
    sage: my_curve(t)=(sin(t), sin(2*t))
    sage: PP=parametric_plot( my_curve, (t, 0, 2*pi), color="purple" )
    sage: my_prime=my_curve.diff(t)
    sage: L=my_prime(1)*t+my_curve(1) # tangent line at t=1
    sage: parametric_plot(L, (t,-2,2))+PP
    Graphics object consisting of 2 graphics primitives

.. tip::
  - After a while, you'll find that giving things names other than ``f``
    and ``g`` becomes quite helpful in distinguishing things from each
    other.  Use descriptive names!  We have tried to do so here.
  - If you are adventurous, try turning this into an interactive cell
    along the lines of the single variable example earlier!

.. _Exam:

'Exam'
------

Before moving out of the calculus world, it is good to have a sort of
miniature exam.

In the cell below, we have plotted and given:

- A slope field for a differential equation,

- A solution to an initial value problem,

- And a symbolic formula for that solution.

We assume you have never seen several of the commands before.  Can you
nonetheless figure out which commands are doing each piece, and what
their syntax is?  How would you look for help to find out more?

.. skip

::

    sage: y = var('y')
    sage: Plot1=plot_slope_field(2-y,(x,0,3),(y,0,20))
    sage: y = function('y',x) # declare y to be a function of x
    sage: h = desolve(diff(y,x) + y - 2, y, ics=[0,7])
    sage: Plot2=plot(h,0,3)
    sage: show(expand(h)); show(Plot1+Plot2)

.. MATH::

    5 \, e^{\left(-x\right)} + 2

Ready to see the answers?  Do not peek until you have really tried it.

In this cell we do the following:

- Make sure that ``y`` is indeed a variable for the first plot.

- Create a slope field for the DE for appropriate inputs of :math:`x`
  and :math:`y`, and give the plot the name ``Plot1``.

- Use the formalism of the function command to get ready for the DE.

 - Notice we have here once again used ``#`` to indicate a comment.

 - In this case, in order to use common terminology, we now have told
   Sage ``y`` is no longer a variable, but instead a function (abstract) of
   the variable ``x``.

- Use the differential equation solving command, with **I**nitial
  **C**ondition **S** of 2 and 2.

- Plot the solution and give it the name ``Plot2``.

- Show a simplification of the symbolic version of the solution (which
  we did not know ahead of time!) as well as the sum of the two graphs \-
  the solution against the slope field.

As you gain experience, you will see how to glean what *you* are looking
for from examples in the documentation like this \- which is one of the
real goals of these tutorials.

Congratulations!  You are now armed with the basics of deploying Sage in
the calculus sequence.
