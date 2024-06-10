SciPy
=====
Again I recommend this
http://www.scipy.org/Wiki/Documentation?action=AttachFile&do=get&target=scipy_tutorial.pdf.
There are many useful SciPy modules, in particular scipy.optimize,
scipy.stats, scipy.linalg, scipy.linsolve, scipy.sparse,
scipy.integrate, scipy.fftpack, scipy.signal, scipy.special. Most
of these have relatively good documentation and often you can
figure out what things do from the names of functions. I recommend
exploring them. For example if you do

::

    sage: import scipy
    sage: from scipy import optimize

Then

.. skip

::

    sage: optimize.[tab]

will show a list of available functions. You should see a bunch of
routines for finding minimum of functions. In particular if you do

.. skip

::

    sage: optimize.fmin_cg?

you find it is a routine that uses the conjugate gradient algorithm
to find the minima of a function.

.. skip

::

    sage: scipy.special.[tab]

will show all the special functions that SciPy has. Spending a
little bit of time looking around is a good way to familiarize
yourself with SciPy. One thing that is sort of annoying, is that
often if you do scipy.:math:`\langle` tab :math:`\rangle`. You
won't see a module that is importable. For example
scipy.:math:`\langle` tab :math:`\rangle` will not show a
signal module but

::

    sage: from scipy import signal

and then

.. skip

::

    signal.[tab]

will show you a large number of functions for signal processing and
filter design. All the modules I listed above can be imported even
if you can't see them initially.

scipy.integrate
---------------

This module has routines related to numerically solving ODE's and
numerical integration. Lets give an example of using an ODE solver.
Suppose you want to solve the ode

    :math:`x''(t) + ux'(t)(x(t)^2-1)+x(t)=0`


which as a system reads

    :math:`x'=y`


    :math:`y'=-x+\mu y(1-x^2).`


The module we want to use is odeint in scipy.integrate. We can
solve this ode, computing the value of :math:`(y,y')`, at 1000
points between :math:`0`, and :math:`100` using the following
code.

::

    sage: import numpy
    sage: from scipy import integrate
    sage: def f_1(y,t):
    ....:    return[y[1],-y[0]-10*y[1]*(y[0]**2-1)]
    sage: def j_1(y,t):
    ....:    return [ [0, 1.0],[-2.0*10*y[0]*y[1]-1.0,-10*(y[0]*y[0]-1.0)] ]
    sage: x = numpy.arange(0,100,.1)
    sage: y = integrate.odeint(f_1,[1,0],x,Dfun=j_1)

We could plot the solution if we wanted by doing

.. link

::

    sage: pts = [(x[i],y[i][0]) for i in range(len(x))]
    sage: point2d(pts).show()

Optimization
------------

The Optimization module has routines related to finding roots,
least squares fitting, and minimization. :math:`\langle` To be
Written :math:`\rangle`
