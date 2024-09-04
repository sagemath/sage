r"""
Solving ODE numerically by GSL

AUTHORS:

- Joshua Kantor (2004-2006)

- Robert Marik (2010 - fixed docstrings)
"""

# ****************************************************************************
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.misc.sageinspect import sage_getargspec
from sage.libs.gsl.all cimport *
import sage.calculus.interpolation


cdef class PyFunctionWrapper:
    cdef object the_function
    cdef object the_jacobian
    cdef object the_parameters
    cdef int y_n

    cdef set_yn(self, x):
        self.y_n = x

cdef class ode_system:
    cdef int c_j(self, double t, double *y, double *dfdy, double *dfdt) noexcept:
        return 0

    cdef int c_f(self, double t, double* y, double* dydt) noexcept:
        return 0

cdef int c_jac_compiled(double t, const double *y, double *dfdy, double *dfdt, void *params) noexcept:
    cdef int status
    cdef ode_system wrapper
    wrapper = <ode_system> params
    status = wrapper.c_j(t, y, dfdy, dfdt)    # Could add parameters
    return status

cdef int c_f_compiled(double t, const double *y, double *dydt, void *params) noexcept:
    cdef int status
    cdef ode_system wrapper
    wrapper = <ode_system> params
    status = wrapper.c_f(t, y, dydt)  # Could add parameters
    return status

cdef int c_jac(double t, const double *y, double *dfdy, double *dfdt, void *params) noexcept:
    cdef int i
    cdef int j
    cdef int y_n
    cdef PyFunctionWrapper wrapper
    wrapper = <PyFunctionWrapper > params
    y_n = wrapper.y_n
    y_list = [y[i] for i in range(y_n)]
    try:
        if len(wrapper.the_parameters) == 0:
            jac_list = wrapper.the_jacobian(t, y_list)
        else:
            jac_list = wrapper.the_jacobian(t, y_list, wrapper.the_parameters)
        for i in range(y_n):
            for j in range(y_n):
                dfdy[i * y_n + j] = jac_list[i][j]

        for i in range(y_n):
            dfdt[i] = jac_list[y_n][i]

        return GSL_SUCCESS
    except Exception:
        return -1

cdef int c_f(double t, const double *y, double *dydt, void *params) noexcept:
    cdef int i
    cdef int y_n

    cdef PyFunctionWrapper wrapper
    wrapper = <PyFunctionWrapper> params
    y_n = wrapper.y_n
    y_list = []
    for i in range(y_n):
        y_list.append(y[i])
    try:
        if len(wrapper.the_parameters) != 0:
            dydt_list = wrapper.the_function(t, y_list, wrapper.the_parameters)
        else:
            dydt_list = wrapper.the_function(t, y_list)
        for i in range(y_n):
            dydt[i] = dydt_list[i]
        return GSL_SUCCESS
    except Exception:
        return -1


class ode_solver():
    r"""
    :meth:`ode_solver` is a class that wraps the GSL library's ode solver routines.

    To use it, instantiate the class::

        sage: T = ode_solver()

    To solve a system of the form `dy_i/dt=f_i(t,y)`, you must
    supply a vector or tuple/list valued function ``f`` representing
    `f_i`.  The functions `f` and the jacobian should have the
    form ``foo(t,y)`` or ``foo(t,y,params)``.  ``params`` which is
    optional allows for your function to depend on one or a tuple of
    parameters.  Note if you use it, ``params`` must be a tuple even
    if it only has one component.  For example if you wanted to solve
    `y''+y=0`, you would need to write it as a first order system::

        y_0' = y_1
        y_1' = -y_0

    In code::

        sage: f = lambda t, y: [y[1], -y[0]]
        sage: T.function = f

    For some algorithms, the jacobian must be supplied as well, the
    form of this should be a function returning a list of lists of the
    form ``[ [df_1/dy_1,...,df_1/dy_n], ...,
    [df_n/dy_1,...,df_n,dy_n], [df_1/dt,...,df_n/dt] ]``.

    There are examples below, if your jacobian was the function
    ``my_jacobian`` you would do::

        sage: T.jacobian = my_jacobian     # not tested, since it doesn't make sense to test this

    There are a variety of algorithms available for different types of systems. Possible algorithms are

    - ``'rkf45'`` -- Runge-Kutta-Fehlberg (4,5)

    - ``'rk2'`` -- embedded Runge-Kutta (2,3)

    - ``'rk4'`` -- `4`-th order classical Runge-Kutta

    - ``'rk8pd'`` -- Runge-Kutta Prince-Dormand (8,9)

    - ``'rk2imp'`` -- implicit 2nd order Runge-Kutta at gaussian points

    - ``'rk4imp'`` -- implicit `4`-th order Runge-Kutta at gaussian points

    - ``'bsimp'`` -- implicit Burlisch-Stoer (requires jacobian)

    - ``'gear1'`` -- M=1 implicit gear

    - ``'gear2'`` -- M=2 implicit gear

    The default algorithm is ``'rkf45'``. If you instead wanted to use
    ``'bsimp'`` you would do::

        sage: T.algorithm = "bsimp"

    The user should supply initial conditions in ``y_0``. For example if
    your initial conditions are `y_0=1, y_1=1`, do::

        sage: T.y_0 = [1,1]

    The actual solver is invoked by the method :meth:`ode_solve`.  It
    has arguments ``t_span``, ``y_0``, ``num_points``, ``params``.
    ``y_0`` must be supplied either as an argument or above by
    assignment.  Params which are optional and only necessary if your
    system uses ``params`` can be supplied to ``ode_solve`` or by
    assignment.

    ``t_span`` is the time interval on which to solve the ode.  There
    are two ways to specify ``t_span``:

    * If ``num_points`` is not specified, then the sequence ``t_span``
      is used as the time points for the solution.  Note that the
      first element ``t_span[0]`` is the initial time, where the
      initial condition ``y_0`` is the specified solution, and
      subsequent elements are the ones where the solution is computed.

    * If ``num_points`` is specified and ``t_span`` is a sequence with
      just 2 elements, then these are the starting and ending times,
      and the solution will be computed at ``num_points`` equally
      spaced points between ``t_span[0]`` and ``t_span[1]``.  The
      initial condition is also included in the output so that
      ``num_points + 1`` total points are returned.  E.g. if ``t_span
      = [0.0, 1.0]`` and ``num_points = 10``, then solution is
      returned at the 11 time points ``[0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
      0.6, 0.7, 0.8, 0.9, 1.0]``.

    (Note that if ``num_points`` is specified and ``t_span`` is not
    length 2 then ``t_span`` are used as the time points and
    ``num_points`` is ignored.)

    Error is estimated via the expression ``D_i =
    error_abs*s_i+error_rel*(a|y_i|+a_dydt*h*|y_i'|)``.  The user can
    specify

    - ``error_abs`` (1e-10 by default),
    - ``error_rel`` (1e-10 by default),
    - ``a`` (1 by default),
    - ``a_dydt`` (0 by default) and
    - ``s_i`` (as ``scaling_abs`` which should be a tuple and is 1 in all
      components by default).

    If you specify one of ``a`` or ``a_dydt``
    you must specify the other.  You may specify ``a`` and ``a_dydt``
    without ``scaling_abs`` (which will be taken =1 be default).
    ``h`` is the initial step size, which is 1e-2 by default.

    ``ode_solve`` solves the solution as a list of tuples of the form,
    ``[ (t_0,[y_1,...,y_n]),(t_1,[y_1,...,y_n]),...,(t_n,[y_1,...,y_n])]``.

    This data is stored in the variable solutions::

        sage: T.solution               # not tested

    EXAMPLES:

    Consider solving the Van der Pol oscillator `x''(t) +
    ux'(t)(x(t)^2-1)+x(t)=0` between `t=0` and `t= 100`.  As a first
    order system it is `x'=y`, `y'=-x+uy(1-x^2)`. Let us take `u=10`
    and use initial conditions `(x,y)=(1,0)` and use the Runge-Kutta
    Prince-Dormand algorithm. ::

        sage: def f_1(t, y, params):
        ....:    return [y[1], -y[0] - params[0]*y[1]*(y[0]**2-1.0)]

        sage: def j_1(t, y, params):
        ....:    return [[0.0, 1.0],
        ....:            [-2.0*params[0]*y[0]*y[1] - 1.0, -params[0]*(y[0]*y[0]-1.0)],
        ....:            [0.0, 0.0]]

        sage: T = ode_solver()
        sage: T.algorithm = "rk8pd"
        sage: T.function = f_1
        sage: T.jacobian = j_1
        sage: T.ode_solve(y_0=[1,0], t_span=[0,100], params=[10.0], num_points=1000)
        sage: import tempfile
        sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:                     # needs sage.plot
        ....:     T.plot_solution(filename=f.name)

    The solver line is equivalent to::

        sage: T.ode_solve(y_0=[1,0], t_span=[x/10.0 for x in range(1000)], params=[10.0])

    Let's try a system::

        y_0'=y_1*y_2
        y_1'=-y_0*y_2
        y_2'=-.51*y_0*y_1

    We will not use the jacobian this time and will change the
    error tolerances. ::

        sage: g_1 = lambda t,y: [y[1]*y[2], -y[0]*y[2], -0.51*y[0]*y[1]]
        sage: T.function = g_1
        sage: T.y_0 = [0,1,1]
        sage: T.scale_abs = [1e-4, 1e-4, 1e-5]
        sage: T.error_rel = 1e-4
        sage: T.ode_solve(t_span=[0,12], num_points=100)

    By default ``T.plot_solution()`` plots the `y_0`; to plot general `y_i`, use::

        sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:                     # needs sage.plot
        ....:     T.plot_solution(i=0, filename=f.name)
        ....:     T.plot_solution(i=1, filename=f.name)
        ....:     T.plot_solution(i=2, filename=f.name)

    The method interpolate_solution will return a spline interpolation
    through the points found by the solver. By default, `y_0` is
    interpolated.  You can interpolate `y_i` through the keyword
    argument ``i``. ::

        sage: f = T.interpolate_solution()
        sage: plot(f,0,12).show()                                                       # needs sage.plot
        sage: f = T.interpolate_solution(i=1)
        sage: plot(f,0,12).show()                                                       # needs sage.plot
        sage: f = T.interpolate_solution(i=2)
        sage: plot(f,0,12).show()                                                       # needs sage.plot
        sage: f = T.interpolate_solution()
        sage: from math import pi
        sage: f(pi)
        0.5379...

    The solver attributes may also be set up using arguments to
    ode_solver.  The previous example can be rewritten as::

        sage: T = ode_solver(g_1, y_0=[0,1,1], scale_abs=[1e-4,1e-4,1e-5],
        ....:                error_rel=1e-4, algorithm='rk8pd')
        sage: T.ode_solve(t_span=[0,12], num_points=100)
        sage: f = T.interpolate_solution()
        sage: f(pi)
        0.5379...

    Unfortunately because Python functions are used, this solver
    is slow on systems that require many function evaluations.  It
    is possible to pass a compiled function by deriving from the
    class :class:`ode_system` and overloading ``c_f`` and ``c_j`` with C
    functions that specify the system. The following will work in the
    notebook:

    .. code-block:: cython

          %cython
          cimport sage.calculus.ode
          import sage.calculus.ode
          from sage.libs.gsl.all cimport *

          cdef class van_der_pol(sage.calculus.ode.ode_system):
              cdef int c_f(self,double t, double *y,double *dydt):
                  dydt[0]=y[1]
                  dydt[1]=-y[0]-1000*y[1]*(y[0]*y[0]-1)
                  return GSL_SUCCESS
              cdef int c_j(self, double t,double *y,double *dfdy,double *dfdt):
                  dfdy[0]=0
                  dfdy[1]=1.0
                  dfdy[2]=-2.0*1000*y[0]*y[1]-1.0
                  dfdy[3]=-1000*(y[0]*y[0]-1.0)
                  dfdt[0]=0
                  dfdt[1]=0
                  return GSL_SUCCESS

    After executing the above block of code you can do the
    following (WARNING: the following is *not* automatically
    doctested)::

        sage: # not tested
        sage: T = ode_solver()
        sage: T.algorithm = "bsimp"
        sage: vander = van_der_pol()
        sage: T.function = vander
        sage: T.ode_solve(y_0=[1, 0], t_span=[0, 2000],
        ....:             num_points=1000)
        sage: from tempfile import NamedTemporaryFile
        sage: with NamedTemporaryFile(suffix='.png') as f:
        ....:     T.plot_solution(i=0, filename=f.name)
    """
    def __init__(self, function=None, jacobian=None, h=1e-2, error_abs=1e-10, error_rel=1e-10, a=False, a_dydt=False, scale_abs=False, algorithm='rkf45', y_0=None, t_span=None, params=[]):
        self.function = function
        self.jacobian = jacobian
        self.h = h
        self.error_abs = error_abs
        self.error_rel = error_rel
        self.a = a
        self.a_dydt = a_dydt
        self.scale_abs = scale_abs
        self.algorithm = algorithm
        self.y_0 = y_0
        self.t_span = t_span
        self.params = params
        self.solution = []

    def __setattr__(self, name, value):
        if hasattr(self, 'solution'):
            object.__setattr__(self, 'solution', [])
        object.__setattr__(self, name, value)

    def interpolate_solution(self, i=0):
        pts = [(t, y[i]) for t, y in self.solution]
        return sage.calculus.interpolation.spline(pts)

    def plot_solution(self, i=0, filename=None, interpolate=False, **kwds):
        r"""
        Plot a one dimensional projection of the solution.

        INPUT:

        - ``i`` -- nonnegative integer; composant of the projection

        - ``filename`` -- string or ``None``; whether to plot the picture or
          save it in a file

        - ``interpolate`` -- whether to interpolate between the points of the
          discretized solution

        - additional keywords are passed to the graphics primitive

        EXAMPLES::

            sage: T = ode_solver()
            sage: T.function = lambda t,y: [cos(y[0]) * sin(t)]
            sage: T.jacobian = lambda t,y: [[-sin(y[0]) * sin(t)]]
            sage: T.ode_solve(y_0=[1],t_span=[0,20],num_points=1000)
            sage: T.plot_solution()                                                     # needs sage.plot

        And with some options::

            sage: T.plot_solution(color='red', axes_labels=["t", "x(t)"])               # needs sage.plot
        """
        if interpolate:
            from sage.plot.line import line2d
            pts = self.interpolate_solution(i)
            G = line2d(pts, **kwds)
        else:
            pts = [(t, y[i]) for t, y in self.solution]
            from sage.plot.point import point2d
            G = point2d([(t, y[i]) for t, y in self.solution], **kwds)
        if filename is None:
            G.show()
        else:
            G.save(filename=filename)

    def ode_solve(self, t_span=False, y_0=False, num_points=False, params=[]):
        cdef double h  # step size
        h = self.h
        cdef int i
        cdef int j
        cdef int type
        cdef int dim
        cdef PyFunctionWrapper wrapper  # struct to pass information into GSL C function
        self.params = params

        if t_span:
            self.t_span = t_span
        if y_0:
            self.y_0 = y_0

        dim = len(self.y_0)
        type = isinstance(self.function, ode_system)
        if type == 0:
            wrapper = PyFunctionWrapper()
            if self.function is not None:
                wrapper.the_function = self.function
            else:
                raise ValueError("ODE system not yet defined")
            if self.jacobian is None:
                wrapper.the_jacobian = None
            else:
                wrapper.the_jacobian = self.jacobian
            if self.params == [] and len(sage_getargspec(wrapper.the_function)[0]) == 2:
                wrapper.the_parameters = []
            elif self.params == [] and len(sage_getargspec(wrapper.the_function)[0]) > 2:
                raise ValueError("ODE system has a parameter but no parameters specified")
            elif self.params != []:
                wrapper.the_parameters = self.params
            wrapper.y_n = dim

        cdef double t
        cdef double t_end
        cdef double *y
        cdef double * scale_abs_array
        scale_abs_array=NULL

        y = <double*> sig_malloc(sizeof(double) * dim)
        if y == NULL:
            raise MemoryError("error allocating memory")
        result = []
        v = [0] * dim
        cdef gsl_odeiv_step_type * T

        for i in range(dim):  # copy initial conditions into C array
            y[i] = self.y_0[i]

        if self.algorithm == "rkf45":
            T = gsl_odeiv_step_rkf45
        elif self.algorithm == "rk2":
            T = gsl_odeiv_step_rk2
        elif self.algorithm == "rk4":
            T = gsl_odeiv_step_rk4
        elif self.algorithm == "rkck":
            T = gsl_odeiv_step_rkck
        elif self.algorithm == "rk8pd":
            T = gsl_odeiv_step_rk8pd
        elif self.algorithm == "rk2imp":
            T = gsl_odeiv_step_rk2imp
        elif self.algorithm == "rk4imp":
            T = gsl_odeiv_step_rk4imp
        elif self.algorithm == "bsimp":
            T = gsl_odeiv_step_bsimp
            if not type and self.jacobian is None:
                raise TypeError("The jacobian must be provided for the implicit Burlisch-Stoer method")
        elif self.algorithm == "gear1":
            T = gsl_odeiv_step_gear1
        elif self.algorithm == "gear2":
            T = gsl_odeiv_step_gear2
        else:
            raise TypeError("algorithm not valid")

        cdef gsl_odeiv_step * s
        s = gsl_odeiv_step_alloc(T, dim)
        if s == NULL:
            sig_free(y)
            raise MemoryError("error setting up solver")

        cdef gsl_odeiv_control * c

        if not self.a and not self.a_dydt:
            c = gsl_odeiv_control_y_new(self.error_abs, self.error_rel)
        elif self.a and self.a_dydt:
            if not self.scale_abs:
                c = gsl_odeiv_control_standard_new(self.error_abs, self.error_rel, self.a, self.a_dydt)
            elif hasattr(self.scale_abs, '__len__'):
                if len(self.scale_abs) == dim:
                    scale_abs_array = <double *> sig_malloc(dim * sizeof(double))
                    for i in range(dim):
                        scale_abs_array[i] = self.scale_abs[i]
                    c = gsl_odeiv_control_scaled_new(self.error_abs, self.error_rel, self.a, self.a_dydt, scale_abs_array, dim)

        if c == NULL:
            gsl_odeiv_control_free(c)
            gsl_odeiv_step_free(s)
            sig_free(y)
            sig_free(scale_abs_array)
            raise MemoryError("error setting up solver")

        cdef gsl_odeiv_evolve * e
        e = gsl_odeiv_evolve_alloc(dim)

        if e == NULL:
            gsl_odeiv_control_free(c)
            gsl_odeiv_step_free(s)
            sig_free(y)
            sig_free(scale_abs_array)
            raise MemoryError("error setting up solver")

        cdef gsl_odeiv_system sys
        if type:               # The user has passed a class with a compiled function, use that for the system
            sys.function = c_f_compiled
            sys.jacobian = c_jac_compiled
#         (<ode_system>self.function).the_parameters = self.params
            sys.params = <void *> self.function
        else:                  # The user passed a python function.
            sys.function = c_f
            sys.jacobian = c_jac
            sys.params = <void *> wrapper
        sys.dimension = dim

        cdef int status
        import copy
        cdef int n

        if len(self.t_span) == 2 and num_points:
            try:
                n = num_points
            except TypeError:
                gsl_odeiv_evolve_free(e)
                gsl_odeiv_control_free(c)
                gsl_odeiv_step_free(s)
                sig_free(y)
                sig_free(scale_abs_array)
                raise TypeError("numpoints must be integer")
            result.append((self.t_span[0], self.y_0))
            delta = (self.t_span[1] - self.t_span[0]) / (1.0 * num_points)
            t = self.t_span[0]
            t_end = self.t_span[0] + delta
            for i in range(1, n + 1):
                while t < t_end:
                    try:
                        sig_on()
                        status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t_end, &h, y)
                        sig_off()
                        if status != GSL_SUCCESS:
                            raise RuntimeError
                    except RuntimeError:
                        gsl_odeiv_evolve_free(e)
                        gsl_odeiv_control_free(c)
                        gsl_odeiv_step_free(s)
                        sig_free(y)
                        sig_free(scale_abs_array)
                        raise ValueError("error solving")

                for j in range(dim):
                    v[j] = <double> y[j]
                result.append((t, copy.copy(v)))
                t = t_end
                t_end = t + delta

        else:
            n = len(self.t_span)
            result.append((self.t_span[0], self.y_0))
            t = self.t_span[0]
            for i in range(1, n):
                t_end = self.t_span[i]
                while t < t_end:
                    try:
                        sig_on()
                        status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t_end, &h, y)
                        sig_off()
                        if status != GSL_SUCCESS:
                            raise RuntimeError
                    except RuntimeError:
                        gsl_odeiv_evolve_free(e)
                        gsl_odeiv_control_free(c)
                        gsl_odeiv_step_free(s)
                        sig_free(y)
                        sig_free(scale_abs_array)
                        raise ValueError("error solving")

                for j in range(dim):
                    v[j] = <double> y[j]
                result.append((t, copy.copy(v)))

                t = self.t_span[i]

        gsl_odeiv_evolve_free(e)
        gsl_odeiv_control_free(c)
        gsl_odeiv_step_free(s)
        sig_free(y)
        sig_free(scale_abs_array)
        self.solution = result
