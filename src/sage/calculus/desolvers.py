r"""
Solving ordinary differential equations

This file contains functions useful for solving differential equations
which occur commonly in a 1st semester differential equations
course. For another numerical solver see the :meth:`ode_solver` function
and the optional package Octave.

Solutions from the Maxima package can contain the three constants
``_C``, ``_K1``, and ``_K2`` where the underscore is used to distinguish
them from symbolic variables that the user might have used. You can
substitute values for them, and make them into accessible usable
symbolic variables, for example with ``var("_C")``.

Commands:

- :func:`desolve` -- compute the "general solution" to a 1st or 2nd order
  ODE via Maxima

- :func:`desolve_laplace` -- solve an ODE using Laplace transforms via
  Maxima. Initial conditions are optional

- :func:`desolve_rk4` -- solve numerically an IVP for one first order
  equation, return list of points or plot

- :func:`desolve_system_rk4` -- solve numerically an IVP for a system of first
  order equations, return list of points

- :func:`desolve_odeint` -- solve numerically a system of first-order ordinary
  differential equations using :func:`~scipy:scipy.integrate.odeint` from
  the module :mod:`scipy:scipy.integrate`.

- :func:`desolve_system` -- solve a system of 1st order ODEs of any size using
  Maxima. Initial conditions are optional

- :func:`eulers_method` -- approximate solution to a 1st order DE,
  presented as a table

- :func:`eulers_method_2x2` -- approximate solution to a 1st order system
  of DEs, presented as a table

- :func:`eulers_method_2x2_plot` -- plot the sequence of points obtained
  from Euler's method

The following functions require the optional package ``tides``:

- :func:`desolve_mintides` -- numerical solution of a system of 1st order ODEs via
  the Taylor series integrator method implemented in TIDES

- :func:`desolve_tides_mpfr` -- arbitrary precision Taylor series integrator implemented in TIDES

AUTHORS:

- David Joyner (3-2006) - Initial version of functions

- Marshall Hampton (7-2007) - Creation of Python module and testing

- Robert Bradshaw (10-2008) - Some interface cleanup.

- Robert Marik (10-2009) - Some bugfixes and enhancements

- Miguel Marco (06-2014) - Tides desolvers
"""

##########################################################################
#  Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>, Marshall Hampton,
#  Robert Marik <marik@mendelu.cz>
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  https://www.gnu.org/licenses/
##########################################################################

import shutil
import os

from sage.interfaces.maxima import Maxima
from sage.misc.functional import N
from sage.rings.real_mpfr import RealField
from sage.structure.element import Expression

from .functional import diff


maxima = Maxima()


def fricas_desolve(de, dvar, ics, ivar):
    r"""
    Solve an ODE using FriCAS.

    EXAMPLES::

        sage: x = var('x')
        sage: y = function('y')(x)
        sage: desolve(diff(y,x) + y - 1, y, algorithm='fricas')             # optional - fricas
        _C0*e^(-x) + 1

        sage: desolve(diff(y, x) + y == y^3*sin(x), y, algorithm='fricas')  # optional - fricas
        -1/5*(2*cos(x)*y(x)^2 + 4*sin(x)*y(x)^2 - 5)*e^(-2*x)/y(x)^2

    TESTS::

        sage: from sage.calculus.desolvers import fricas_desolve
        sage: Y = fricas_desolve(diff(y,x) + y - 1, y, [42,1783], x)  # optional - fricas
        sage: Y.subs(x=42)                                            # optional - fricas
        1783
    """
    from sage.interfaces.fricas import fricas
    from sage.symbolic.ring import SR
    if ics is None:
        y = fricas(de).solve(dvar.operator(), ivar).sage()
    else:
        eq = fricas.equation(ivar, ics[0])
        y = fricas(de).solve(dvar.operator(), eq, ics[1:]).sage()

    if isinstance(y, dict):
        basis = y["basis"]
        particular = y["particular"]
        return particular + sum(SR.var("_C" + str(i)) * v
                                for i, v in enumerate(basis))
    return y


def fricas_desolve_system(des, dvars, ics, ivar):
    r"""
    Solve a system of first order ODEs using FriCAS.

    EXAMPLES::

        sage: t = var('t')
        sage: x = function('x')(t)
        sage: y = function('y')(t)
        sage: de1 = diff(x,t) + y - 1 == 0
        sage: de2 = diff(y,t) - x + 1 == 0
        sage: desolve_system([de1, de2], [x, y], algorithm='fricas')          # optional - fricas
        [x(t) == _C0*cos(t) + cos(t)^2 + _C1*sin(t) + sin(t)^2,
         y(t) == -_C1*cos(t) + _C0*sin(t) + 1]

        sage: desolve_system([de1, de2], [x,y], [0,1,2], algorithm='fricas')  # optional - fricas
        [x(t) == cos(t)^2 + sin(t)^2 - sin(t), y(t) == cos(t) + 1]

    TESTS::

        sage: from sage.calculus.desolvers import fricas_desolve_system
        sage: t = var('t')
        sage: x = function('x')(t)
        sage: fricas_desolve_system([diff(x,t) + 1 == 0], [x], None, t)   # optional - fricas
        [x(t) == _C0 - t]

        sage: y = function('y')(t)
        sage: de1 = diff(x,t) + y - 1 == 0
        sage: de2 = diff(y,t) - x + 1 == 0
        sage: sol = fricas_desolve_system([de1,de2], [x,y], [0,1,-1], t)  # optional - fricas
        sage: sol                                                         # optional - fricas
        [x(t) == cos(t)^2 + sin(t)^2 + 2*sin(t), y(t) == -2*cos(t) + 1]
    """
    from sage.interfaces.fricas import fricas
    from sage.symbolic.ring import SR
    from sage.symbolic.relation import solve
    ops = [dvar.operator() for dvar in dvars]
    y = fricas(des).solve(ops, ivar).sage()
    basis = y["basis"]
    particular = y["particular"]
    pars = [SR.var("_C" + str(i)) for i in range(len(basis))]
    solv = particular + sum(p * v for p, v in zip(pars, basis))

    if ics is None:
        sols = solv
    else:
        ics0 = ics[0]
        eqs = [val == sol.subs({ivar: ics0}) for val, sol in zip(ics[1:], solv)]
        pars_values = solve(eqs, pars, solution_dict=True)
        sols = [sol.subs(pars_values[0]) for sol in solv]

    return [dvar == sol for dvar, sol in zip(dvars, sols)]


def desolve(de, dvar, ics=None, ivar=None, show_method=False, contrib_ode=False,
            algorithm='maxima'):
    r"""
    Solve a 1st or 2nd order linear ODE, including IVP and BVP.

    INPUT:

    - ``de`` -- an expression or equation representing the ODE

    - ``dvar`` -- the dependent variable (hereafter called `y`)

    - ``ics`` -- (optional) the initial or boundary conditions

      - for a first-order equation, specify the initial `x` and `y`

      - for a second-order equation, specify the initial `x`, `y`,
        and `dy/dx`, i.e. write `[x_0, y(x_0), y'(x_0)]`

      - for a second-order boundary solution, specify initial and
        final `x` and `y` boundary conditions, i.e. write `[x_0, y(x_0), x_1, y(x_1)]`.

      - gives an error if the solution is not SymbolicEquation (as happens for
        example for a Clairaut equation)

    - ``ivar`` -- (optional) the independent variable (hereafter called
      `x`), which must be specified if there is more than one
      independent variable in the equation

    - ``show_method`` -- (optional) if ``True``, then Sage returns pair
      ``[solution, method]``, where method is the string describing
      the method which has been used to get a solution (Maxima uses the
      following order for first order equations: linear, separable,
      exact (including exact with integrating factor), homogeneous,
      bernoulli, generalized homogeneous) - use carefully in class,
      see below the example of an equation which is separable but
      this property is not recognized by Maxima and the equation is solved
      as exact.

    - ``contrib_ode`` -- (optional) if ``True``, ``desolve`` allows to solve
      Clairaut, Lagrange, Riccati and some other equations. This may take
      a long time and is thus turned off by default.  Initial conditions
      can be used only if the result is one SymbolicEquation (does not
      contain a singular solution, for example).

    - ``algorithm`` -- (default: ``'maxima'``) one of

      * ``'maxima'`` -- use maxima
      * ``'fricas'`` -- use FriCAS (the optional fricas spkg has to be installed)

    OUTPUT:

    In most cases return a SymbolicEquation which defines the solution
    implicitly.  If the result is in the form `y(x)=\ldots` (happens for
    linear eqs.), return the right-hand side only. The possible
    constant solutions of separable ODEs are omitted.

    .. NOTE::

        Use ``desolve? <tab>`` if the output in the Sage notebook is truncated.

    EXAMPLES::

        sage: x = var('x')
        sage: y = function('y')(x)
        sage: desolve(diff(y,x) + y - 1, y)
        (_C + e^x)*e^(-x)

    ::

        sage: f = desolve(diff(y,x) + y - 1, y, ics=[10,2]); f
        (e^10 + e^x)*e^(-x)

    ::

        sage: plot(f)
        Graphics object consisting of 1 graphics primitive

    We can also solve second-order differential equations::

        sage: x = var('x')
        sage: y = function('y')(x)
        sage: de = diff(y,x,2) - y == x
        sage: desolve(de, y)
        _K2*e^(-x) + _K1*e^x - x

    ::

        sage: f = desolve(de, y, [10,2,1]); f
        -x + 7*e^(x - 10) + 5*e^(-x + 10)

    ::

        sage: f(x=10)
        2

    ::

        sage: diff(f,x)(x=10)
        1

    ::

        sage: de = diff(y,x,2) + y == 0
        sage: desolve(de, y)
        _K2*cos(x) + _K1*sin(x)

    ::

        sage: desolve(de, y, [0,1,pi/2,4])
        cos(x) + 4*sin(x)

    ::

        sage: desolve(y*diff(y,x)+sin(x)==0,y)
        -1/2*y(x)^2 == _C - cos(x)

    Clairaut equation: general and singular solutions::

        sage: desolve(diff(y,x)^2+x*diff(y,x)-y==0,y,contrib_ode=True,show_method=True)
        [[y(x) == _C^2 + _C*x, y(x) == -1/4*x^2], 'clairau...']

    For equations involving more variables we specify an independent variable::

        sage: a,b,c,n=var('a b c n')
        sage: desolve(x^2*diff(y,x)==a+b*x^n+c*x^2*y^2,y,ivar=x,contrib_ode=True)
        [[y(x) == 0, (b*x^(n - 2) + a/x^2)*c^2*u == 0]]

    ::

        sage: desolve(x^2*diff(y,x)==a+b*x^n+c*x^2*y^2,y,ivar=x,contrib_ode=True,show_method=True)
        [[[y(x) == 0, (b*x^(n - 2) + a/x^2)*c^2*u == 0]], 'riccati']


    Higher order equations, not involving independent variable::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y).expand()
        1/6*y(x)^3 + _K1*y(x) == _K2 + x

    ::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,1,3]).expand()
        1/6*y(x)^3 - 5/3*y(x) == x - 3/2

    ::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,1,3],show_method=True)
        [1/6*y(x)^3 - 5/3*y(x) == x - 3/2, 'freeofx']

    Separable equations - Sage returns solution in implicit form::

        sage: desolve(diff(y,x)*sin(y) == cos(x),y)
        -cos(y(x)) == _C + sin(x)

    ::

        sage: desolve(diff(y,x)*sin(y) == cos(x),y,show_method=True)
        [-cos(y(x)) == _C + sin(x), 'separable']

    ::

        sage: desolve(diff(y,x)*sin(y) == cos(x),y,[pi/2,1])
        -cos(y(x)) == -cos(1) + sin(x) - 1

    Linear equation - Sage returns the expression on the right hand side only::

        sage: desolve(diff(y,x)+(y) == cos(x),y)
        1/2*((cos(x) + sin(x))*e^x + 2*_C)*e^(-x)

    ::

        sage: desolve(diff(y,x)+(y) == cos(x),y,show_method=True)
        [1/2*((cos(x) + sin(x))*e^x + 2*_C)*e^(-x), 'linear']

    ::

        sage: desolve(diff(y,x)+(y) == cos(x),y,[0,1])
        1/2*(cos(x)*e^x + e^x*sin(x) + 1)*e^(-x)

    This ODE with separated variables is solved as
    exact. Explanation - factor does not split `e^{x-y}` in Maxima
    into `e^{x}e^{y}`::

        sage: desolve(diff(y,x)==exp(x-y),y,show_method=True)
        [-e^x + e^y(x) == _C, 'exact']

    You can solve Bessel equations, also using initial
    conditions, but you cannot put (sometimes desired) the initial
    condition at `x=0`, since this point is a singular point of the
    equation. Anyway, if the solution should be bounded at `x=0`, then
    ``_K2=0``. ::

        sage: desolve(x^2*diff(y,x,x)+x*diff(y,x)+(x^2-4)*y==0,y)
        _K1*bessel_J(2, x) + _K2*bessel_Y(2, x)

    Example of difficult ODE producing an error::

        sage: desolve(sqrt(y)*diff(y,x)+e^(y)+cos(x)-sin(x+y)==0,y) # not tested
        Traceback (click to the left for traceback)
        ...
        NotImplementedError, "Maxima was unable to solve this ODE. Consider to set option contrib_ode to True."

    Another difficult ODE with error - moreover, it takes a long time::

        sage: desolve(sqrt(y)*diff(y,x)+e^(y)+cos(x)-sin(x+y)==0,y,contrib_ode=True) # not tested

    Some more types of ODEs::

        sage: desolve(x*diff(y,x)^2-(1+x*y)*diff(y,x)+y==0,y,contrib_ode=True,show_method=True)
        [[y(x) == _C + log(x), y(x) == _C*e^x], 'factor']

    ::

        sage: desolve(diff(y,x)==(x+y)^2,y,contrib_ode=True,show_method=True)
        [[[x == _C - arctan(sqrt(t)), y(x) == -x - sqrt(t)], [x == _C + arctan(sqrt(t)), y(x) == -x + sqrt(t)]], 'lagrange']

    These two examples produce an error (as expected, Maxima 5.18 cannot
    solve equations from initial conditions). Maxima 5.18
    returns false answer in this case! ::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,2]).expand() # not tested
        Traceback (click to the left for traceback)
        ...
        NotImplementedError, "Maxima was unable to solve this ODE. Consider to set option contrib_ode to True."

    ::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,2],show_method=True) # not tested
        Traceback (click to the left for traceback)
        ...
        NotImplementedError, "Maxima was unable to solve this ODE. Consider to set option contrib_ode to True."

    Second order linear ODE::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y)
        (_K2*x + _K1)*e^(-x) + 1/2*sin(x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,show_method=True)
        [(_K2*x + _K1)*e^(-x) + 1/2*sin(x), 'variationofparameters']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,1])
        1/2*(7*x + 6)*e^(-x) + 1/2*sin(x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,1],show_method=True)
        [1/2*(7*x + 6)*e^(-x) + 1/2*sin(x), 'variationofparameters']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,pi/2,2])
        3*(x*(e^(1/2*pi) - 2)/pi + 1)*e^(-x) + 1/2*sin(x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,pi/2,2],show_method=True)
        [3*(x*(e^(1/2*pi) - 2)/pi + 1)*e^(-x) + 1/2*sin(x), 'variationofparameters']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y)
        (_K2*x + _K1)*e^(-x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,show_method=True)
        [(_K2*x + _K1)*e^(-x), 'constcoeff']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,1])
        (4*x + 3)*e^(-x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,1],show_method=True)
        [(4*x + 3)*e^(-x), 'constcoeff']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,pi/2,2])
        (2*x*(2*e^(1/2*pi) - 3)/pi + 3)*e^(-x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,pi/2,2],show_method=True)
        [(2*x*(2*e^(1/2*pi) - 3)/pi + 3)*e^(-x), 'constcoeff']

    Using ``algorithm='fricas'`` we can invoke the differential
    equation solver from FriCAS.  For example, it can solve higher
    order linear equations::

        sage: de = x^3*diff(y, x, 3) + x^2*diff(y, x, 2) - 2*x*diff(y, x) + 2*y - 2*x^4
        sage: Y = desolve(de, y, algorithm='fricas'); Y               # optional - fricas
        (2*x^3 - 3*x^2 + 1)*_C0/x + (x^3 - 1)*_C1/x
         + (x^3 - 3*x^2 - 1)*_C2/x + 1/15*(x^5 - 10*x^3 + 20*x^2 + 4)/x

    The initial conditions are then interpreted as `[x_0, y(x_0),
    y'(x_0), \ldots, y^(n)(x_0)]`::

        sage: Y = desolve(de, y, ics=[1,3,7], algorithm='fricas'); Y  # optional - fricas
        1/15*(x^5 + 15*x^3 + 50*x^2 - 21)/x

    FriCAS can also solve some non-linear equations::

        sage: de = diff(y, x) == y / (x+y*log(y))
        sage: Y = desolve(de, y, algorithm='fricas'); Y               # optional - fricas
        1/2*(log(y(x))^2*y(x) - 2*x)/y(x)

    TESTS:

    :issue:`9961` fixed (allow assumptions on the dependent variable in desolve)::

        sage: y=function('y')(x); assume(x>0); assume(y>0)
        sage: sage.calculus.calculus.maxima('domain:real')  # needed since Maxima 5.26.0 to get the answer as below
        real
        sage: desolve(x*diff(y,x)-x*sqrt(y^2+x^2)-y == 0, y, contrib_ode=True)
        [x - arcsinh(y(x)/x) == _C]

    :issue:`10682` updated Maxima to 5.26, and it started to show a different
    solution in the complex domain for the ODE above::

        sage: forget()
        sage: sage.calculus.calculus.maxima('domain:complex')  # back to the default complex domain
        complex
        sage: assume(x>0)
        sage: assume(y>0)
        sage: desolve(x*diff(y,x)-x*sqrt(y^2+x^2)-y == 0, y, contrib_ode=True)
        [x - arcsinh(y(x)^2/(x*sqrt(y(x)^2))) - arcsinh(y(x)/x) + 1/2*log(4*(x^2 + 2*y(x)^2 + 2*sqrt(x^2*y(x)^2 + y(x)^4))/x^2) == _C]

    :issue:`6479` fixed::

        sage: x = var('x')
        sage: y = function('y')(x)
        sage: desolve( diff(y,x,x) == 0, y, [0,0,1])
        x

    ::

        sage: desolve( diff(y,x,x) == 0, y, [0,1,1])
        x + 1

    :issue:`9835` fixed::

        sage: x = var('x')
        sage: y = function('y')(x)
        sage: desolve(diff(y,x,2)+y*(1-y^2)==0,y,[0,-1,1,1])
        Traceback (most recent call last):
        ...
        NotImplementedError: Unable to use initial condition for this equation (freeofx).

    :issue:`8931` fixed::

        sage: x=var('x'); f=function('f')(x); k=var('k'); assume(k>0)
        sage: desolve(diff(f,x,2)/f==k,f,ivar=x)
        _K1*e^(sqrt(k)*x) + _K2*e^(-sqrt(k)*x)

    :issue:`15775` fixed::

        sage: forget()
        sage: y = function('y')(x)
        sage: desolve(diff(y, x) == sqrt(abs(y)), dvar=y, ivar=x)
        integrate(1/sqrt(abs(y(x))), y(x)) == _C + x

    AUTHORS:

    - David Joyner (1-2006)
    - Robert Bradshaw (10-2008)
    - Robert Marik (10-2009)
    """
    if isinstance(de, Expression) and de.is_relational():
        de = de.lhs() - de.rhs()
    if isinstance(dvar, Expression) and dvar.is_symbol():
        raise ValueError("You have to declare dependent variable as a function evaluated at the independent variable, eg. y=function('y')(x)")
    # for backwards compatibility
    if isinstance(dvar, list):
        dvar, ivar = dvar
    elif ivar is None:
        ivars = de.variables()
        ivars = [t for t in ivars if t is not dvar]
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = ivars[0]

    if algorithm == "fricas":
        return fricas_desolve(de, dvar, ics, ivar)
    elif algorithm != "maxima":
        raise ValueError("unknown algorithm %s" % algorithm)

    de00 = de._maxima_()
    P = de00.parent()
    dvar_str = P(dvar.operator()).str()
    ivar_str = P(ivar).str()
    de00 = de00.str()

    def sanitize_var(exprs):
        return exprs.replace("'" + dvar_str + "(" + ivar_str + ")", dvar_str)
    de0 = sanitize_var(de00)
    ode_solver = "ode2"
    cmd = "(TEMP:%s(%s,%s,%s), if TEMP=false then TEMP else substitute(%s=%s(%s),TEMP))" % (ode_solver, de0, dvar_str, ivar_str, dvar_str, dvar_str, ivar_str)
    # we produce string like this
    # ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y(x),x)
    soln = P(cmd)

    if str(soln).strip() == 'false':
        if contrib_ode:
            ode_solver = "contrib_ode"
            P("load('contrib_ode)")
            cmd = "(TEMP:%s(%s,%s,%s), if TEMP=false then TEMP else substitute(%s=%s(%s),TEMP))" % (ode_solver, de0, dvar_str, ivar_str, dvar_str, dvar_str, ivar_str)
            # we produce string like this
            # (TEMP:contrib_ode(x*('diff(y,x,1))^2-(x*y+1)*'diff(y,x,1)+y,y,x), if TEMP=false then TEMP else substitute(y=y(x),TEMP))
            soln = P(cmd)
            if str(soln).strip() == 'false':
                raise NotImplementedError("Maxima was unable to solve this ODE.")
        else:
            raise NotImplementedError("Maxima was unable to solve this ODE. Consider to set option contrib_ode to True.")

    if show_method:
        maxima_method = P("method")

    if ics is not None:
        if not (isinstance(soln.sage(), Expression) and soln.sage().is_relational()):
            if not show_method:
                maxima_method = P("method")
            raise NotImplementedError("Unable to use initial condition for this equation (%s)." % (str(maxima_method).strip()))
        if len(ics) == 2:
            tempic = (ivar == ics[0])._maxima_().str()
            tempic = tempic + "," + (dvar == ics[1])._maxima_().str()
            cmd = "(TEMP:ic1(%s(%s,%s,%s),%s),substitute(%s=%s(%s),TEMP))" % (ode_solver, de00, dvar_str, ivar_str, tempic, dvar_str, dvar_str, ivar_str)
            cmd = sanitize_var(cmd)
            # we produce string like this
            # (TEMP:ic2(ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y,x),x=0,y=3,'diff(y,x)=1),substitute(y=y(x),TEMP))
            soln = P(cmd)
        if len(ics) == 3:
            # fixed ic2 command from Maxima - we have to ensure that %k1, %k2 do not depend on variables, should be removed when fixed in Maxima
            P("ic2_sage(soln,xa,ya,dya):=block([programmode:true,backsubst:true,singsolve:true,temp,%k2,%k1,TEMP_k], \
                noteqn(xa), noteqn(ya), noteqn(dya), boundtest('%k1,%k1), boundtest('%k2,%k2), \
                temp: lhs(soln) - rhs(soln), \
                TEMP_k:solve([subst([xa,ya],soln), subst([dya,xa], lhs(dya)=-subst(0,lhs(dya),diff(temp,lhs(xa)))/diff(temp,lhs(ya)))],[%k1,%k2]), \
                if not freeof(lhs(ya),TEMP_k) or not freeof(lhs(xa),TEMP_k) then return (false), \
                temp: maplist(lambda([zz], subst(zz,soln)), TEMP_k), \
                if length(temp)=1 then return(first(temp)) else return(temp))")
            tempic = P(ivar == ics[0]).str()
            tempic += "," + P(dvar == ics[1]).str()
            tempic += ",'diff(" + dvar_str + "," + ivar_str + ")=" + P(ics[2]).str()
            cmd = "(TEMP:ic2_sage(%s(%s,%s,%s),%s),substitute(%s=%s(%s),TEMP))" % (ode_solver, de00, dvar_str, ivar_str, tempic, dvar_str, dvar_str, ivar_str)
            cmd = sanitize_var(cmd)
            # we produce string like this
            # (TEMP:ic2(ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y,x),x=0,y=3,'diff(y,x)=1),substitute(y=y(x),TEMP))
            soln = P(cmd)
            if str(soln).strip() == 'false':
                raise NotImplementedError("Maxima was unable to solve this IVP. Remove the initial condition to get the general solution.")
        if len(ics) == 4:
            # fixed bc2 command from Maxima - we have to ensure that %k1, %k2 do not depend on variables, should be removed when fixed in Maxima
            P("bc2_sage(soln,xa,ya,xb,yb):=block([programmode:true,backsubst:true,singsolve:true,temp,%k1,%k2,TEMP_k], \
                noteqn(xa), noteqn(ya), noteqn(xb), noteqn(yb), boundtest('%k1,%k1), boundtest('%k2,%k2), \
                TEMP_k:solve([subst([xa,ya],soln), subst([xb,yb],soln)], [%k1,%k2]), \
                if not freeof(lhs(ya),TEMP_k) or not freeof(lhs(xa),TEMP_k) then return (false), \
                temp: maplist(lambda([zz], subst(zz,soln)),TEMP_k), \
                if length(temp)=1 then return(first(temp)) else return(temp))")
            cmd = "bc2_sage(%s(%s,%s,%s),%s,%s=%s,%s,%s=%s)" % (ode_solver, de00, dvar_str, ivar_str, P(ivar == ics[0]).str(), dvar_str, P(ics[1]).str(), P(ivar == ics[2]).str(), dvar_str, P(ics[3]).str())
            cmd = "(TEMP:%s,substitute(%s=%s(%s),TEMP))" % (cmd, dvar_str, dvar_str, ivar_str)
            cmd = sanitize_var(cmd)
            # we produce string like this
            # (TEMP:bc2(ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y,x),x=0,y=3,x=%pi/2,y=2),substitute(y=y(x),TEMP))
            soln = P(cmd)
            if str(soln).strip() == 'false':
                raise NotImplementedError("Maxima was unable to solve this BVP. Remove the initial condition to get the general solution.")

    soln = soln.sage()
    if isinstance(soln, Expression) and soln.is_relational() and soln.lhs() == dvar:
        # Remark: Here we do not check that the right hand side does not depend on dvar.
        # This probably will not happen for solutions obtained via ode2, anyway.
        soln = soln.rhs()
    if show_method:
        return [soln, maxima_method.str()]
    return soln


def desolve_laplace(de, dvar, ics=None, ivar=None):
    """
    Solve an ODE using Laplace transforms. Initial conditions are optional.

    INPUT:

    - ``de`` -- a lambda expression representing the ODE (e.g. ``de =
      diff(y,x,2) == diff(y,x)+sin(x)``)

    - ``dvar`` -- the dependent variable (e.g. ``y``)

    - ``ivar`` -- (optional) the independent variable (hereafter called
      `x`), which must be specified if there is more than one
      independent variable in the equation.

    - ``ics`` -- list of numbers representing initial conditions, (e.g.
      ``f(0)=1``, ``f'(0)=2`` corresponds to ``ics = [0,1,2]``)

    OUTPUT: solution of the ODE as symbolic expression

    EXAMPLES::

        sage: u=function('u')(x)
        sage: eq = diff(u,x) - exp(-x) - u == 0
        sage: desolve_laplace(eq,u)
        1/2*(2*u(0) + 1)*e^x - 1/2*e^(-x)

    We can use initial conditions::

        sage: desolve_laplace(eq,u,ics=[0,3])
        -1/2*e^(-x) + 7/2*e^x

    The initial conditions do not persist in the system (as they persisted
    in previous versions)::

        sage: desolve_laplace(eq,u)
        1/2*(2*u(0) + 1)*e^x - 1/2*e^(-x)

    ::

        sage: f=function('f')(x)
        sage: eq = diff(f,x) + f == 0
        sage: desolve_laplace(eq,f,[0,1])
        e^(-x)

    ::

        sage: x = var('x')
        sage: f = function('f')(x)
        sage: de = diff(f,x,x) - 2*diff(f,x) + f
        sage: desolve_laplace(de,f)
        -x*e^x*f(0) + x*e^x*D[0](f)(0) + e^x*f(0)

    ::

        sage: desolve_laplace(de,f,ics=[0,1,2])
        x*e^x + e^x

    TESTS:

    Check that :issue:`4839` is fixed::

        sage: t = var('t')
        sage: x = function('x')(t)
        sage: soln = desolve_laplace(diff(x,t)+x==1, x, ics=[0,2])
        sage: soln
        e^(-t) + 1

    ::

        sage: soln(t=3)
        e^(-3) + 1

    AUTHORS:

    - David Joyner (1-2006,8-2007)

    - Robert Marik (10-2009)
    """
    # This is the original code from David Joyner (inputs and outputs strings)
    # maxima("de:"+de._repr_()+"=0;")
    # if ics is not None:
    #     d = len(ics)
    #     for i in range(0,d-1):
    #         ic = "atvalue(diff("+vars[1]+"("+vars[0]+"),"+str(vars[0])+","+str(i)+"),"+str(vars[0])+"="+str(ics[0])+","+str(ics[1+i])+")"
    #         maxima(ic)
    #
    # cmd = "desolve("+de._repr_()+","+vars[1]+"("+vars[0]+"));"
    # return maxima(cmd).rhs()._maxima_init_()

    # verbatim copy from desolve - begin
    if isinstance(de, Expression) and de.is_relational():
        de = de.lhs() - de.rhs()
    if isinstance(dvar, Expression) and dvar.is_symbol():
        raise ValueError("You have to declare dependent variable as a function evaluated at the independent variable, eg. y=function('y')(x)")
    # for backwards compatibility
    if isinstance(dvar, list):
        dvar, ivar = dvar
    elif ivar is None:
        ivars = de.variables()
        ivars = [t for t in ivars if t is not dvar]
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = ivars[0]
    # verbatim copy from desolve - end

    dvar_str = str(dvar)

    def sanitize_var(exprs):
        # 'y(x) -> y(x)
        return exprs.replace("'" + dvar_str, dvar_str)

    de0 = de._maxima_()
    P = de0.parent()
    i = dvar_str.find('(')
    dvar_str = dvar_str[:i + 1] + '_SAGE_VAR_' + dvar_str[i + 1:]
    cmd = sanitize_var("desolve(" + de0.str() + "," + dvar_str + ")")
    soln = P(cmd).rhs()
    if str(soln).strip() == 'false':
        raise NotImplementedError("Maxima was unable to solve this ODE.")
    soln = soln.sage()
    if ics is not None:
        d = len(ics)
        for i in range(d - 1):
            soln = eval('soln.substitute(diff(dvar,ivar,i)(' + str(ivar) + '=ics[0])==ics[i+1])')
    return soln


def desolve_system(des, vars, ics=None, ivar=None, algorithm='maxima'):
    r"""
    Solve a system of any size of 1st order ODEs. Initial conditions
    are optional.

    One dimensional systems are passed to :meth:`desolve_laplace`.

    INPUT:

    - ``des`` -- list of ODEs

    - ``vars`` -- list of dependent variables

    - ``ics`` -- (optional) list of initial values for ``ivar`` and ``vars``;
      if ``ics`` is defined, it should provide initial conditions for each
      variable, otherwise an exception would be raised

    - ``ivar`` -- (optional) the independent variable, which must be
      specified if there is more than one independent variable in the
      equation

    - ``algorithm`` -- (default: ``'maxima'``) one of

      * ``'maxima'`` -- use maxima
      * ``'fricas'`` -- use FriCAS (the optional fricas spkg has to be installed)

    EXAMPLES::

        sage: t = var('t')
        sage: x = function('x')(t)
        sage: y = function('y')(t)
        sage: de1 = diff(x,t) + y - 1 == 0
        sage: de2 = diff(y,t) - x + 1 == 0
        sage: desolve_system([de1, de2], [x,y])
        [x(t) == (x(0) - 1)*cos(t) - (y(0) - 1)*sin(t) + 1,
         y(t) == (y(0) - 1)*cos(t) + (x(0) - 1)*sin(t) + 1]

    The same system solved using FriCAS::

        sage: desolve_system([de1, de2], [x,y], algorithm='fricas')  # optional - fricas
        [x(t) == _C0*cos(t) + cos(t)^2 + _C1*sin(t) + sin(t)^2,
         y(t) == -_C1*cos(t) + _C0*sin(t) + 1]

    Now we give some initial conditions::

        sage: sol = desolve_system([de1, de2], [x,y], ics=[0,1,2]); sol
        [x(t) == -sin(t) + 1, y(t) == cos(t) + 1]

    ::

        sage: solnx, solny = sol[0].rhs(), sol[1].rhs()
        sage: plot([solnx,solny],(0,1))  # not tested
        sage: parametric_plot((solnx,solny),(0,1))  # not tested

    TESTS:

    Check that :issue:`9823` is fixed::

        sage: t = var('t')
        sage: x = function('x')(t)
        sage: de1 = diff(x,t) + 1 == 0
        sage: desolve_system([de1], [x])
        -t + x(0)

    Check that :issue:`16568` is fixed::

        sage: t = var('t')
        sage: x = function('x')(t)
        sage: y = function('y')(t)
        sage: de1 = diff(x,t) + y - 1 == 0
        sage: de2 = diff(y,t) - x + 1 == 0
        sage: des = [de1,de2]
        sage: ics = [0,1,-1]
        sage: vars = [x,y]
        sage: sol = desolve_system(des, vars, ics); sol
        [x(t) == 2*sin(t) + 1, y(t) == -2*cos(t) + 1]

    ::

        sage: solx, soly = sol[0].rhs(), sol[1].rhs()
        sage: RR(solx(t=3))
        1.28224001611973

    ::

        sage: P1 = plot([solx,soly], (0,1))                                             # needs sage.plot
        sage: P2 = parametric_plot((solx,soly), (0,1))                                  # needs sage.plot

    Now type ``show(P1)``, ``show(P2)`` to view these plots.

    Check that :issue:`9824` is fixed::

        sage: t = var('t')
        sage: epsilon = var('epsilon')
        sage: x1 = function('x1')(t)
        sage: x2 = function('x2')(t)
        sage: de1 = diff(x1,t) == epsilon
        sage: de2 = diff(x2,t) == -2
        sage: desolve_system([de1, de2], [x1, x2], ivar=t)
        [x1(t) == epsilon*t + x1(0), x2(t) == -2*t + x2(0)]
        sage: desolve_system([de1, de2], [x1, x2], ics=[1,1], ivar=t)
        Traceback (most recent call last):
        ...
        ValueError: Initial conditions aren't complete: number of vars is different
        from number of dependent variables. Got ics = [1, 1], vars = [x1(t), x2(t)]

    Check that :issue:`9825` is fixed::

        sage: t = var('t')
        sage: x1, x2=function("x1, x2")
        sage: de1=x1(t).diff(t)==-3*(x2(t)-1)
        sage: de2=x2(t).diff(t)==1
        sage: Sol=desolve_system([de1, de2],[x1(t),x2(t)],ivar=t) ; Sol
        [x1(t) == -3/2*t^2 - 3*t*x2(0) + 3*t + x1(0), x2(t) == t + x2(0)]


    AUTHORS:

    - Robert Bradshaw (10-2008)
    - Sergey Bykov (10-2014)
    """
    if ics is not None:
        if len(ics) != (len(vars) + 1):
            raise ValueError("Initial conditions aren't complete: number of vars is different from number of dependent variables. Got ics = {0}, vars = {1}".format(ics, vars))

    if len(des) == 1 and algorithm == "maxima":
        return desolve_laplace(des[0], vars[0], ics=ics, ivar=ivar)
    ivars = set([])
    for i, de in enumerate(des):
        if not (isinstance(de, Expression) and de.is_relational()):
            des[i] = de == 0
        ivars = ivars.union(set(de.variables()))
    if ivar is None:
        ivars = ivars - set(vars)
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = list(ivars)[0]

    if algorithm == "fricas":
        return fricas_desolve_system(des, vars, ics, ivar)
    elif algorithm != "maxima":
        raise ValueError("unknown algorithm %s" % algorithm)

    dvars = [v._maxima_() for v in vars]
    if ics is not None:
        ivar_ic = ics[0]
        for dvar, ic in zip(dvars, ics[1:]):
            dvar.atvalue(ivar == ivar_ic, ic)
    soln = dvars[0].parent().desolve(des, dvars)
    if str(soln).strip() == 'false':
        raise NotImplementedError("Maxima was unable to solve this system.")
    soln = list(soln)
    for i, sol in enumerate(soln):
        soln[i] = sol.sage()
    if ics is not None:
        ivar_ic = ics[0]
        for dvar, ic in zip(dvars, ics[:1]):
            dvar.atvalue(ivar == ivar_ic, dvar)
    return soln


def eulers_method(f, x0, y0, h, x1, algorithm='table'):
    r"""
    This implements Euler's method for finding numerically the
    solution of the 1st order ODE `y' = f(x,y)`, `y(a)=c`. The ``x``
    column of the table increments from `x_0` to `x_1` by `h` (so
    `(x_1-x_0)/h` must be an integer). In the ``y`` column, the new
    `y`-value equals the old `y`-value plus the corresponding entry in the
    last column.

    .. NOTE::

        This function is for pedagogical purposes only.

    EXAMPLES::

        sage: from sage.calculus.desolvers import eulers_method
        sage: x,y = PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1)
             x                    y                  h*f(x,y)
             0                    1                   -2
           1/2                   -1                 -7/4
             1                -11/4                -11/8

    ::

        sage: x,y = PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1,algorithm='none')
        [[0, 1], [1/2, -1], [1, -11/4], [3/2, -33/8]]

    ::

        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y = PolynomialRing(RR,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1,algorithm="None")
        [[0, 1], [1/2, -1.0], [1, -2.7], [3/2, -4.0]]

    ::

        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y=PolynomialRing(RR,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1)
             x                    y                  h*f(x,y)
             0                    1                 -2.0
           1/2                 -1.0                 -1.7
             1                 -2.7                 -1.3

    ::

        sage: x,y=PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,1,1,1/3,2)
                 x                    y                  h*f(x,y)
                 1                    1                  1/3
               4/3                  4/3                    1
               5/3                  7/3                 17/9
                 2                 38/9                83/27

    ::

        sage: eulers_method(5*x+y-5,0,1,1/2,1,algorithm='none')
        [[0, 1], [1/2, -1], [1, -11/4], [3/2, -33/8]]

    ::

        sage: pts = eulers_method(5*x+y-5,0,1,1/2,1,algorithm='none')
        sage: P1 = list_plot(pts)                                                       # needs sage.plot
        sage: P2 = line(pts)                                                            # needs sage.plot
        sage: (P1 + P2).show()                                                          # needs sage.plot

    AUTHORS:

    - David Joyner
    """
    if algorithm == "table":
        print("%10s %20s %25s" % ("x", "y", "h*f(x,y)"))
    n = int((1.0) * (x1 - x0) / h)
    x00 = x0
    y00 = y0
    soln = [[x00, y00]]
    for i in range(n + 1):
        if algorithm == "table":
            print("%10r %20r %20r" % (x00, y00, h * f(x00, y00)))
        y00 = y00 + h * f(x00, y00)
        x00 = x00 + h
        soln.append([x00, y00])
    if algorithm != "table":
        return soln


def eulers_method_2x2(f, g, t0, x0, y0, h, t1, algorithm='table'):
    r"""
    This implements Euler's method for finding numerically the
    solution of the 1st order system of two ODEs

    .. MATH::

        \begin{aligned}
        x' &= f(t, x, y), x(t_0)=x_0 \\
        y' &= g(t, x, y), y(t_0)=y_0.
        \end{aligned}

    The ``t`` column of the table increments from `t_0` to `t_1` by `h`
    (so `\frac{t_1-t_0}{h}` must be an integer). In the ``x`` column,
    the new `x`-value equals the old `x`-value plus the corresponding
    entry in the next (third) column.  In the ``y`` column, the new
    `y`-value equals the old `y`-value plus the corresponding entry in the
    next (last) column.

    .. NOTE::

        This function is for pedagogical purposes only.

    EXAMPLES::

        sage: from sage.calculus.desolvers import eulers_method_2x2
        sage: t, x, y = PolynomialRing(QQ,3,"txy").gens()
        sage: f = x+y+t; g = x-y
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1,algorithm='none')
        [[0, 0, 0], [1/3, 0, 0], [2/3, 1/9, 0], [1, 10/27, 1/27], [4/3, 68/81, 4/27]]

    ::

        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    0                         0                    0                    0
           1/3                    0                       1/9                    0                    0
           2/3                  1/9                      7/27                    0                 1/27
             1                10/27                     38/81                 1/27                  1/9

    ::

        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: t,x,y=PolynomialRing(RR,3,"txy").gens()
        sage: f = x+y+t; g = x-y
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    0                      0.00                    0                 0.00
           1/3                 0.00                      0.13                 0.00                 0.00
           2/3                 0.13                      0.29                 0.00                0.043
             1                 0.41                      0.57                0.043                 0.15

    To numerically approximate `y(1)`, where `(1+t^2)y''+y'-y=0`,
    `y(0)=1`, `y'(0)=-1`, using 4 steps of Euler's method, first
    convert to a system: `y_1' = y_2`, `y_1(0)=1`; `y_2' =
    \frac{y_1-y_2}{1+t^2}`, `y_2(0)=-1`.::

         sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
         sage: t, x, y=PolynomialRing(RR,3,"txy").gens()
         sage: f = y; g = (x-y)/(1+t^2)
         sage: eulers_method_2x2(f,g, 0, 1, -1, 1/4, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    1                     -0.25                   -1                 0.50
           1/4                 0.75                     -0.12                -0.50                 0.29
           1/2                 0.63                    -0.054                -0.21                 0.19
           3/4                 0.63                   -0.0078               -0.031                 0.11
             1                 0.63                     0.020                0.079                0.071

    To numerically approximate `y(1)`, where `y''+ty'+y=0`, `y(0)=1`, `y'(0)=0`::

        sage: t,x,y=PolynomialRing(RR,3,"txy").gens()
        sage: f = y; g = -x-y*t
        sage: eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    1                      0.00                    0                -0.25
           1/4                  1.0                    -0.062                -0.25                -0.23
           1/2                 0.94                     -0.11                -0.46                -0.17
           3/4                 0.88                     -0.15                -0.62                -0.10
             1                 0.75                     -0.17                -0.68               -0.015

    AUTHORS:

    - David Joyner
    """
    if algorithm == "table":
        print("%10s %20s %25s %20s %20s" % ("t", "x", "h*f(t,x,y)", "y", "h*g(t,x,y)"))
    n = int((1.0) * (t1 - t0) / h)
    t00 = t0
    x00 = x0
    y00 = y0
    soln = [[t00, x00, y00]]
    for i in range(n + 1):
        if algorithm == "table":
            print("%10r %20r %25r %20r %20r" % (t00, x00, h * f(t00, x00, y00), y00, h * g(t00, x00, y00)))
        x01 = x00 + h * f(t00, x00, y00)
        y00 = y00 + h * g(t00, x00, y00)
        x00 = x01
        t00 = t00 + h
        soln.append([t00, x00, y00])
    if algorithm != "table":
        return soln


def eulers_method_2x2_plot(f, g, t0, x0, y0, h, t1):
    r"""
    Plot solution of ODE.

    This plots the solution in the rectangle with sides ``(xrange[0],xrange[1])`` and
    ``(yrange[0],yrange[1])``, and plots using Euler's method the
    numerical solution of the 1st order ODEs `x' = f(t,x,y)`,
    `x(a)=x_0`, `y' = g(t,x,y)`, `y(a) = y_0`.

    .. NOTE::

        This function is for pedagogical purposes only.

    EXAMPLES:

    The following example plots the solution to
    `\theta''+\sin(\theta)=0`, `\theta(0)=\frac 34`, `\theta'(0) =
    0`.  Type ``P[0].show()`` to plot the solution,
    ``(P[0]+P[1]).show()`` to plot `(t,\theta(t))` and
    `(t,\theta'(t))`::

        sage: from sage.calculus.desolvers import eulers_method_2x2_plot
        sage: f = lambda z : z[2]; g = lambda z : -sin(z[1])
        sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)                 # needs sage.plot
    """
    from sage.plot.line import line

    n = int((1.0) * (t1 - t0) / h)
    t00 = t0
    x00 = x0
    y00 = y0
    soln = [[t00, x00, y00]]
    for i in range(n + 1):
        x01 = x00 + h * f([t00, x00, y00])
        y00 = y00 + h * g([t00, x00, y00])
        x00 = x01
        t00 = t00 + h
        soln.append([t00, x00, y00])
    Q1 = line([[x[0], x[1]] for x in soln], rgbcolor=(.25, .125, .75))
    Q2 = line([[x[0], x[2]] for x in soln], rgbcolor=(.5, .125, .25))
    return [Q1, Q2]


def desolve_rk4_determine_bounds(ics, end_points=None):
    """
    Used to determine bounds for numerical integration.

    - If ``end_points`` is None, the interval for integration is from ``ics[0]``
      to ``ics[0]+10``

    - If ``end_points`` is ``a`` or ``[a]``, the interval for integration is from ``min(ics[0],a)``
      to ``max(ics[0],a)``

    - If ``end_points`` is ``[a,b]``, the interval for integration is from ``min(ics[0],a)``
      to ``max(ics[0],b)``

    EXAMPLES::

        sage: from sage.calculus.desolvers import desolve_rk4_determine_bounds
        sage: desolve_rk4_determine_bounds([0,2],1)
        (0, 1)

    ::

        sage: desolve_rk4_determine_bounds([0,2])
        (0, 10)

    ::

        sage: desolve_rk4_determine_bounds([0,2],[-2])
        (-2, 0)

    ::

        sage: desolve_rk4_determine_bounds([0,2],[-2,4])
        (-2, 4)
    """
    if end_points is None:
        return ics[0], ics[0] + 10
    if not isinstance(end_points, list):
        end_points = [end_points]
    if len(end_points) == 1:
        return min(ics[0], end_points[0]), max(ics[0], end_points[0])
    else:
        return min(ics[0], end_points[0]), max(ics[0], end_points[1])


def desolve_rk4(de, dvar, ics=None, ivar=None, end_points=None, step=0.1, output='list', **kwds):
    """
    Solve numerically one first-order ordinary differential
    equation.

    INPUT:

    Input is similar to ``desolve`` command. The differential equation can be
    written in a form close to the ``plot_slope_field`` or ``desolve`` command.

    - Variant 1 (function in two variables)

      - ``de`` -- right hand side, i.e. the function `f(x,y)` from ODE `y'=f(x,y)`

      - ``dvar`` -- dependent variable (symbolic variable declared by var)

    - Variant 2 (symbolic equation)

      - ``de`` -- equation, including term with ``diff(y,x)``

      - ``dvar`` -- dependent variable (declared as function of independent variable)

    - Other parameters

      - ``ivar`` -- should be specified, if there are more variables or if the equation is autonomous

      - ``ics`` -- initial conditions in the form ``[x0,y0]``

      - ``end_points`` -- the end points of the interval

        - if ``end_points`` is a or [a], we integrate between ``min(ics[0],a)`` and ``max(ics[0],a)``
        - if ``end_points`` is None, we use ``end_points=ics[0]+10``

        - if end_points is [a,b] we integrate between ``min(ics[0], a)`` and ``max(ics[0], b)``

      - ``step`` -- (default: 0.1) the length of the step (positive number)

      - ``output`` -- (default: ``'list'``) one of ``'list'``,
        ``'plot'``, ``'slope_field'`` (graph of the solution with slope field)

    OUTPUT:

    Return a list of points, or plot produced by ``list_plot``,
    optionally with slope field.

    .. SEEALSO::

        :func:`ode_solver`.

    EXAMPLES::

        sage: from sage.calculus.desolvers import desolve_rk4

    Variant 2 for input - more common in numerics::

        sage: x,y = var('x,y')
        sage: desolve_rk4(x*y*(2-y),y,ics=[0,1],end_points=1,step=0.5)
        [[0, 1], [0.5, 1.12419127424558], [1.0, 1.46159016228882...]]

    Variant 1 for input - we can pass ODE in the form used by
    desolve function In this example we integrate backwards, since
    ``end_points < ics[0]``::

        sage: y = function('y')(x)
        sage: desolve_rk4(diff(y,x)+y*(y-1) == x-2,y,ics=[1,1],step=0.5, end_points=0)
        [[0.0, 8.904257108962112], [0.5, 1.90932794536153...], [1, 1]]

    Here we show how to plot simple pictures. For more advanced
    applications use list_plot instead. To see the resulting picture
    use ``show(P)`` in Sage notebook. ::

        sage: x,y = var('x,y')
        sage: P=desolve_rk4(y*(2-y),y,ics=[0,.1],ivar=x,output='slope_field',end_points=[-4,6],thickness=3)

    ALGORITHM:

    `4`-th order Runge-Kutta method. Wrapper for command ``rk`` in
    Maxima's dynamics package.  Perhaps could be faster by using
    fast_float instead.

    AUTHORS:

    - Robert Marik (10-2009)
    """
    if ics is None:
        raise ValueError("No initial conditions, specify with ics=[x0,y0].")

    if output not in ['list', 'plot', 'slope_field']:
        raise ValueError("Option output should be 'list', 'plot' or 'slope_field'.")

    if ivar is None:
        ivars = de.variables()
        ivars = [t for t in ivars if t != dvar]
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = ivars[0]

    step = abs(step)

    def desolve_rk4_inner(de, dvar):
        de0 = de._maxima_()
        maxima("load('dynamics)")
        lower_bound, upper_bound = desolve_rk4_determine_bounds(ics, end_points)
        sol_1, sol_2 = [], []
        if lower_bound < ics[0]:
            cmd = "rk(%s,%s,%s,[%s,%s,%s,%s])\
            " % (de0.str(), '_SAGE_VAR_' + str(dvar), str(ics[1]), '_SAGE_VAR_' + str(ivar), str(ics[0]), lower_bound, -step)
            sol_1 = maxima(cmd).sage()
            sol_1.pop(0)
            sol_1.reverse()
        if upper_bound > ics[0]:
            cmd = "rk(%s,%s,%s,[%s,%s,%s,%s])\
            " % (de0.str(), '_SAGE_VAR_' + str(dvar), str(ics[1]), '_SAGE_VAR_' + str(ivar), str(ics[0]), upper_bound, step)
            sol_2 = maxima(cmd).sage()
            sol_2.pop(0)
        sol = sol_1
        sol.extend([[ics[0], ics[1]]])
        sol.extend(sol_2)

        if output == 'list':
            return sol
        from sage.plot.plot import list_plot
        from sage.plot.plot_field import plot_slope_field
        R = list_plot(sol, plotjoined=True, **kwds)
        if output == 'plot':
            return R
        if output == 'slope_field':
            XMIN = sol[0][0]
            YMIN = sol[0][1]
            XMAX = XMIN
            YMAX = YMIN
            for s, t in sol:
                XMAX = max(s, XMAX)
                XMIN = min(s, XMIN)
                YMAX = max(t, YMAX)
                YMIN = min(t, YMIN)
            return plot_slope_field(de, (ivar, XMIN, XMAX), (dvar, YMIN, YMAX)) + R

    if not (isinstance(dvar, Expression) and dvar.is_symbol()):
        from sage.symbolic.ring import SR
        from sage.calculus.all import diff
        from sage.symbolic.relation import solve
        if isinstance(de, Expression) and de.is_relational():
            de = de.lhs() - de.rhs()
        # consider to add warning if the solution is not unique
        de = solve(de, diff(dvar, ivar), solution_dict=True)
        if len(de) != 1:
            raise NotImplementedError("Sorry, cannot find explicit formula for right-hand side of the ODE.")
        with SR.temp_var() as dummy_dvar:
            return desolve_rk4_inner(de[0][diff(dvar, ivar)].subs({dvar: dummy_dvar}), dummy_dvar)
    else:
        return desolve_rk4_inner(de, dvar)


def desolve_system_rk4(des, vars, ics=None, ivar=None, end_points=None, step=0.1):
    r"""
    Solve numerically a system of first-order ordinary differential
    equations using the `4`-th order Runge-Kutta method. Wrapper for
    Maxima command ``rk``.

    INPUT:

    Input is similar to ``desolve_system`` and ``desolve_rk4`` commands

    - ``des`` -- right hand sides of the system

    - ``vars`` -- dependent variables

    - ``ivar`` -- (optional) should be specified, if there are more variables or
      if the equation is autonomous and the independent variable is
      missing

    - ``ics`` -- initial conditions in the form ``[x0,y01,y02,y03,....]``

    - ``end_points`` -- the end points of the interval

      - if ``end_points`` is a or [a], we integrate on between ``min(ics[0], a)`` and ``max(ics[0], a)``
      - if ``end_points`` is None, we use ``end_points=ics[0]+10``

      - if ``end_points`` is [a,b] we integrate on between ``min(ics[0], a)`` and ``max(ics[0], b)``

    - ``step`` -- (default: 0.1) the length of the step

    OUTPUT: a list of points

    .. SEEALSO::

        :func:`ode_solver`.

    EXAMPLES::

        sage: from sage.calculus.desolvers import desolve_system_rk4

    Lotka Volterra system::

        sage: from sage.calculus.desolvers import desolve_system_rk4
        sage: x,y,t = var('x y t')
        sage: P = desolve_system_rk4([x*(1-y),-y*(1-x)], [x,y], ics=[0,0.5,2],
        ....:                        ivar=t, end_points=20)
        sage: Q = [[i,j] for i,j,k in P]
        sage: LP = list_plot(Q)                                                         # needs sage.plot

        sage: Q = [[j,k] for i,j,k in P]
        sage: LP = list_plot(Q)                                                         # needs sage.plot

    ALGORITHM:

    `4`-th order Runge-Kutta method. Wrapper for command ``rk`` in Maxima's
    dynamics package.  Perhaps could be faster by using ``fast_float``
    instead.

    AUTHOR:

    - Robert Marik (10-2009)
    """

    if ics is None:
        raise ValueError("No initial conditions, specify with ics=[x0,y01,y02,...].")

    ivars = set([])

    for de in des:
        ivars = ivars.union(set(de.variables()))
    if ivar is None:
        ivars = ivars - set(vars)
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = list(ivars)[0]

    dess = [de._maxima_().str() for de in des]
    desstr = "[" + ",".join(dess) + "]"
    varss = [varsi._maxima_().str() for varsi in vars]
    varstr = "[" + ",".join(varss) + "]"
    x0 = ics[0]
    icss = [ics[i]._maxima_().str() for i in range(1, len(ics))]
    icstr = "[" + ",".join(icss) + "]"
    step = abs(step)

    maxima("load('dynamics)")
    lower_bound, upper_bound = desolve_rk4_determine_bounds(ics, end_points)
    sol_1, sol_2 = [], []
    if lower_bound < ics[0]:
        cmd = "rk(%s,%s,%s,[%s,%s,%s,%s])\
        " % (desstr, varstr, icstr, '_SAGE_VAR_' + str(ivar), str(x0), lower_bound, -step)
        sol_1 = maxima(cmd).sage()
        sol_1.pop(0)
        sol_1.reverse()
    if upper_bound > ics[0]:
        cmd = "rk(%s,%s,%s,[%s,%s,%s,%s])\
        " % (desstr, varstr, icstr, '_SAGE_VAR_' + str(ivar), str(x0), upper_bound, step)
        sol_2 = maxima(cmd).sage()
        sol_2.pop(0)
    sol = sol_1
    sol.append(ics)
    sol.extend(sol_2)

    return sol


def desolve_odeint(des, ics, times, dvars, ivar=None, compute_jac=False, args=(),
                   rtol=None, atol=None, tcrit=None, h0=0.0, hmax=0.0, hmin=0.0, ixpr=0,
                   mxstep=0, mxhnil=0, mxordn=12, mxords=5, printmessg=0):
    r"""
    Solve numerically a system of first-order ordinary differential equations
    using :func:`scipy:scipy.integrate.odeint`.

    INPUT:

    - ``des`` -- right hand sides of the system

    - ``ics`` -- initial conditions

    - ``times`` -- a sequence of time points in which the solution must be found

    - ``dvars`` -- dependent variables. ATTENTION: the order must be the same as
      in ``des``, that means: ``d(dvars[i])/dt=des[i]``

    - ``ivar`` -- independent variable, optional

    - ``compute_jac`` -- boolean (default: ``False``); if ``True``, the
      Jacobian of ``des`` is computed and used during the integration of stiff
      systems

    Other Parameters (taken from the documentation of the
    :func:`~scipy:scipy.integrate.odeint` function from
    :mod:`scipy:scipy.integrate`):

    - ``rtol``, ``atol`` : float
      The input parameters ``rtol`` and ``atol`` determine the error
      control performed by the solver.  The solver will control the
      vector, `e`, of estimated local errors in `y`, according to an
      inequality of the form:

        max-norm of (e / ewt) <= 1

      where ewt is a vector of positive error weights computed as:

        ewt = rtol * abs(y) + atol

      ``rtol`` and ``atol`` can be either vectors the same length as `y` or scalars.

    - ``tcrit`` : array
      Vector of critical points (e.g. singularities) where integration
      care should be taken.

    - ``h0`` : float, (0: solver-determined)
      The step size to be attempted on the first step.

    - ``hmax`` : float, (0: solver-determined)
      The maximum absolute step size allowed.

    - ``hmin`` : float, (0: solver-determined)
      The minimum absolute step size allowed.

    - ``ixpr`` : boolean.
      Whether to generate extra printing at method switches.

    - ``mxstep`` : integer, (0: solver-determined)
      Maximum number of (internally defined) steps allowed for each
      integration point in t.

    - ``mxhnil`` : integer, (0: solver-determined)
      Maximum number of messages printed.

    - ``mxordn`` : integer, (0: solver-determined)
      Maximum order to be allowed for the nonstiff (Adams) method.

    - ``mxords`` : integer, (0: solver-determined)
      Maximum order to be allowed for the stiff (BDF) method.

    OUTPUT: a list with the solution of the system at each time in ``times``

    EXAMPLES:

    Lotka Volterra Equations::

        sage: from sage.calculus.desolvers import desolve_odeint
        sage: x,y = var('x,y')
        sage: f = [x*(1-y), -y*(1-x)]
        sage: sol = desolve_odeint(f, [0.5,2], srange(0,10,0.1), [x,y])                 # needs scipy
        sage: p = line(zip(sol[:,0],sol[:,1]))                                          # needs scipy sage.plot
        sage: p.show()                                                                  # needs scipy sage.plot

    Lorenz Equations::

        sage: x,y,z = var('x,y,z')
        sage: # Next we define the parameters
        sage: sigma = 10
        sage: rho = 28
        sage: beta = 8/3
        sage: # The Lorenz equations
        sage: lorenz = [sigma*(y-x),x*(rho-z)-y,x*y-beta*z]
        sage: # Time and initial conditions
        sage: times = srange(0,50.05,0.05)
        sage: ics = [0,1,1]
        sage: sol = desolve_odeint(lorenz, ics, times, [x,y,z],                         # needs scipy
        ....:                      rtol=1e-13, atol=1e-14)

    One-dimensional stiff system::

        sage: y = var('y')
        sage: epsilon = 0.01
        sage: f = y^2*(1-y)
        sage: ic = epsilon
        sage: t = srange(0,2/epsilon,1)
        sage: sol = desolve_odeint(f, ic, t, y,                                         # needs scipy
        ....:                      rtol=1e-9, atol=1e-10, compute_jac=True)
        sage: p = points(zip(t,sol[:,0]))                                               # needs scipy sage.plot
        sage: p.show()                                                                  # needs scipy sage.plot

    Another stiff system with some optional parameters with no
    default value::

        sage: y1,y2,y3 = var('y1,y2,y3')
        sage: f1 = 77.27*(y2+y1*(1-8.375*1e-6*y1-y2))
        sage: f2 = 1/77.27*(y3-(1+y1)*y2)
        sage: f3 = 0.16*(y1-y3)
        sage: f = [f1,f2,f3]
        sage: ci = [0.2,0.4,0.7]
        sage: t = srange(0,10,0.01)
        sage: v = [y1,y2,y3]
        sage: sol = desolve_odeint(f, ci, t, v, rtol=1e-3, atol=1e-4,                   # needs scipy
        ....:                      h0=0.1, hmax=1, hmin=1e-4, mxstep=1000, mxords=17)

    AUTHOR:

    - Oriol Castejon (05-2010)
    """

    from scipy.integrate import odeint
    from sage.ext.fast_eval import fast_float
    from sage.calculus.functions import jacobian

    def desolve_odeint_inner(ivar):
        # one-dimensional systems:
        if len(dvars) == 1:
            assert len(des) == 1
            dvar = dvars[0]
            de = des[0]
            func = fast_float(de, dvar, ivar)
            if not compute_jac:
                Dfun = None
            else:
                J = diff(de, dvar)
                J = fast_float(J, dvar, ivar)

                def Dfun(y, t):
                    return [J(y.item(), t)]

        # n-dimensional systems:
        else:
            variabs = dvars[:]
            variabs.append(ivar)
            desc = [fast_float(de, *variabs) for de in des]

            def func(y, t):
                v = list(y[:])
                v.append(t)
                return [dec(*v) for dec in desc]

            if not compute_jac:
                Dfun = None
            else:
                J = jacobian(des, dvars)
                J = [list(v) for v in J]
                J = fast_float(J, *variabs)

                def Dfun(y, t):
                    v = list(y[:])
                    v.append(t)
                    return [[element(*v) for element in row] for row in J]

        return odeint(func, ics, times, args=args, Dfun=Dfun, rtol=rtol, atol=atol,
                      tcrit=tcrit, h0=h0, hmax=hmax, hmin=hmin, ixpr=ixpr, mxstep=mxstep,
                      mxhnil=mxhnil, mxordn=mxordn, mxords=mxords, printmessg=printmessg)

    if isinstance(dvars, Expression) and dvars.is_symbol():
        dvars = [dvars]

    if not isinstance(des, (list, tuple)):
        des = [des]

    if ivar is None:
        all_vars = set().union(*[de.variables() for de in des])
        ivars = all_vars - set(dvars)

        if len(ivars) == 1:
            return desolve_odeint_inner(next(iter(ivars)))
        elif not ivars:
            from sage.symbolic.ring import SR
            with SR.temp_var() as ivar:
                return desolve_odeint_inner(ivar)
        else:
            raise ValueError("Unable to determine independent variable, please specify.")
    return desolve_odeint_inner(ivar)


def desolve_mintides(f, ics, initial, final, delta, tolrel=1e-16, tolabs=1e-16):
    r"""
    Solve numerically a system of first order differential equations using the
    taylor series integrator implemented in mintides.

    INPUT:

    - ``f`` -- symbolic function; its first argument will be the independent
      variable, . Its output should be de derivatives of the dependent variables.

    - ``ics`` -- list or tuple with the initial conditions

    - ``initial`` -- the starting value for the independent variable

    - ``final`` -- the final value for the independent value

    - ``delta`` -- the size of the steps in the output

    - ``tolrel`` -- the relative tolerance for the method

    - ``tolabs`` -- the absolute tolerance for the method

    OUTPUT: list with the positions of the IVP

    EXAMPLES:

    We integrate a periodic orbit of the Kepler problem along 50 periods::

        sage: var('t,x,y,X,Y')
        (t, x, y, X, Y)
        sage: f(t,x,y,X,Y)=[X, Y, -x/(x^2+y^2)^(3/2), -y/(x^2+y^2)^(3/2)]
        sage: ics = [0.8, 0, 0, 1.22474487139159]
        sage: t = 100*pi
        sage: sol = desolve_mintides(f, ics, 0, t, t, 1e-12, 1e-12) # optional -tides
        sage: sol # optional -tides # abs tol 1e-5
        [[0.000000000000000,
        0.800000000000000,
        0.000000000000000,
        0.000000000000000,
        1.22474487139159],
        [314.159265358979,
        0.800000000028622,
        -5.91973525754241e-9,
        7.56887091890590e-9,
        1.22474487136329]]


    ALGORITHM:

    Uses TIDES.

    REFERENCES:

    - A. Abad, R. Barrio, F. Blesa, M. Rodriguez. Algorithm 924. *ACM
      Transactions on Mathematical Software* , *39* (1), 1-28.

    - A. Abad, R. Barrio, F. Blesa, M. Rodriguez.
      `TIDES tutorial: Integrating ODEs by using the Taylor Series Method.
      <http://www.unizar.es/acz/05Publicaciones/Monografias/MonografiasPublicadas/Monografia36/IndMonogr36.htm>`_
    """
    import subprocess
    if subprocess.call('command -v gcc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
        raise RuntimeError('Unable to run because gcc cannot be found')
    from sage.interfaces.tides import genfiles_mintides
    from sage.misc.temporary_file import tmp_dir
    tempdir = tmp_dir()
    intfile = os.path.join(tempdir, 'integrator.c')
    drfile = os.path.join(tempdir, 'driver.c')
    fileoutput = os.path.join(tempdir, 'output')
    runmefile = os.path.join(tempdir, 'runme')
    genfiles_mintides(intfile, drfile, f, [N(_) for _ in ics],
                      N(initial), N(final), N(delta), N(tolrel),
                      N(tolabs), fileoutput)
    subprocess.check_call('gcc -o ' + runmefile + ' ' + os.path.join(tempdir, '*.c ') +
                          os.path.join('$SAGE_LOCAL', 'lib', 'libTIDES.a') + ' $LDFLAGS '
                          + os.path.join('-L$SAGE_LOCAL', 'lib ') + ' -lm  -O2 ' +
                          os.path.join('-I$SAGE_LOCAL', 'include '),
                          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.check_call(os.path.join(tempdir, 'runme'), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open(fileoutput) as outfile:
        res = outfile.readlines()
    for i in range(len(res)):
        res[i] = [RealField()(_) for _ in res[i].split(' ') if len(_) > 2]
    shutil.rmtree(tempdir)
    return res


def desolve_tides_mpfr(f, ics, initial, final, delta, tolrel=1e-16, tolabs=1e-16, digits=50):
    r"""
    Solve numerically a system of first order differential equations using the
    taylor series integrator in arbitrary precision implemented in tides.

    INPUT:

    - ``f`` -- symbolic function; its first argument will be the independent
      variable. Its output should be de derivatives of the dependent variables.

    - ``ics`` -- list or tuple with the initial conditions

    - ``initial`` -- the starting value for the independent variable

    - ``final`` -- the final value for the independent value

    - ``delta`` -- the size of the steps in the output

    - ``tolrel`` -- the relative tolerance for the method

    - ``tolabs`` -- the absolute tolerance for the method

    - ``digits`` -- the digits of precision used in the computation

    OUTPUT: list with the positions of the IVP

    EXAMPLES:

    We integrate the Lorenz equations with Saltzman values for the parameters
    along 10 periodic orbits with 100 digits of precision::

        sage: var('t,x,y,z')
        (t, x, y, z)
        sage: s = 10
        sage: r = 28
        sage: b = 8/3
        sage: f(t,x,y,z)= [s*(y-x),x*(r-z)-y,x*y-b*z]
        sage: x0 = -13.7636106821342005250144010543616538641008648540923684535378642921202827747268115852940239346395038284
        sage: y0 = -19.5787519424517955388380414460095588661142400534276438649791334295426354746147526415973165506704676171
        sage: z0 = 27
        sage: T = 15.586522107161747275678702092126960705284805489972439358895215783190198756258880854355851082660142374
        sage: sol = desolve_tides_mpfr(f, [x0, y0, z0], 0, T, T, 1e-100, 1e-100, 100) # optional - tides
        sage: sol # optional -tides # abs tol 1e-50
        [[0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
        -13.7636106821342005250144010543616538641008648540923684535378642921202827747268115852940239346395038,
        -19.5787519424517955388380414460095588661142400534276438649791334295426354746147526415973165506704676,
        27.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000],
        [15.5865221071617472756787020921269607052848054899724393588952157831901987562588808543558510826601424,
        -13.7636106821342005250144010543616538641008648540923684535378642921202827747268115852940239346315658,
        -19.5787519424517955388380414460095588661142400534276438649791334295426354746147526415973165506778440,
        26.9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999636628]]


    ALGORITHM:

    Uses TIDES.


    .. WARNING::

        This requires the package tides.

    REFERENCES:

    - [ABBR2011]_
    - [ABBR2012]_
    """
    import subprocess
    if subprocess.call('command -v gcc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
        raise RuntimeError('Unable to run because gcc cannot be found')
    from sage.interfaces.tides import genfiles_mpfr
    from sage.functions.other import ceil
    from sage.functions.log import log
    from sage.misc.temporary_file import tmp_dir
    tempdir = tmp_dir()
    intfile = os.path.join(tempdir, 'integrator.c')
    drfile = os.path.join(tempdir, 'driver.c')
    fileoutput = os.path.join(tempdir, 'output')
    runmefile = os.path.join(tempdir, 'runme')
    genfiles_mpfr(intfile, drfile, f, ics, initial, final, delta, [], [],
                  digits, tolrel, tolabs, fileoutput)
    subprocess.check_call('gcc -o ' + runmefile + ' ' + os.path.join(tempdir, '*.c ') +
                          os.path.join('$SAGE_LOCAL', 'lib', 'libTIDES.a') + ' $LDFLAGS '
                          + os.path.join('-L$SAGE_LOCAL', 'lib ') + '-lmpfr -lgmp -lm  -O2 -w ' +
                          os.path.join('-I$SAGE_LOCAL', 'include '),
                          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.check_call(os.path.join(tempdir, 'runme'), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open(fileoutput) as outfile:
        res = outfile.readlines()
    for i in range(len(res)):
        res[i] = [RealField(ceil(digits * log(10, 2)))(piece)
                  for piece in res[i].split(' ') if len(piece) > 2]
    shutil.rmtree(tempdir)
    return res
