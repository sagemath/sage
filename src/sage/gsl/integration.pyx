"""
Numerical Integration

AUTHORS:

- Josh Kantor (2007-02): first version

- William Stein (2007-02): rewrite of docs, conventions, etc.

- Robert Bradshaw (2008-08): fast float integration

- Jeroen Demeyer (2011-11-23): :trac:`12047`: return 0 when the
  integration interval is a point; reformat documentation and add to
  the reference manual.
"""

#*****************************************************************************
#       Copyright (C) 2004,2005,2006,2007 Joshua Kantor <kantor.jm@gmail.com>
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

include "cysignals/signals.pxi"
from sage.libs.gsl.all cimport *
from sage.ext.fast_eval cimport FastDoubleFunc


cdef class PyFunctionWrapper:
   cdef object the_function
   cdef object the_parameters

cdef class compiled_integrand:
   cdef int c_f(self,double t):  #void *params):
      return 0

cdef double c_f(double t,void *params):
   cdef double value
   cdef PyFunctionWrapper wrapper
   wrapper = <PyFunctionWrapper> params
   try:
      if len(wrapper.the_parameters)!=0:
         value=wrapper.the_function(t,wrapper.the_parameters)
      else:
         value=wrapper.the_function(t)
   except Exception as msg:
      print(msg)
      return 0

   return value


cdef double c_ff(double t, void *params):
    return (<FastDoubleFunc>params)._call_c(&t)


def numerical_integral(func, a, b=None,
                       algorithm='qag',
                       max_points=87, params=[], eps_abs=1e-6,
                       eps_rel=1e-6, rule=6):
   r"""
    Returns the numerical integral of the function on the interval
    from a to b and an error bound.

    INPUT:

    - ``a``, ``b`` -- The interval of integration, specified as two
      numbers or as a tuple/list with the first element the lower bound
      and the second element the upper bound.  Use ``+Infinity`` and
      ``-Infinity`` for plus or minus infinity.
    - ``algorithm`` -- valid choices are:

      * 'qag' -- for an adaptive integration
      * 'qags' -- for an adaptive integration with (integrable) singularities
      * 'qng' -- for a non-adaptive Gauss-Kronrod (samples at a maximum of 87pts)

    - ``max_points`` -- sets the maximum number of sample points
    - ``params`` -- used to pass parameters to your function
    - ``eps_abs``, ``eps_rel`` -- absolute and relative error tolerances
    - ``rule`` -- This controls the Gauss-Kronrod rule used in the adaptive integration:

      * rule=1 -- 15 point rule
      * rule=2 -- 21 point rule
      * rule=3 -- 31 point rule
      * rule=4 -- 41 point rule
      * rule=5 -- 51 point rule
      * rule=6 -- 61 point rule

      Higher key values are more accurate for smooth functions but lower
      key values deal better with discontinuities.

    OUTPUT:

    A tuple whose first component is the answer and whose second
    component is an error estimate.

    REMARK:

    There is also a method ``nintegral`` on symbolic expressions
    that implements numerical integration using Maxima.  It is potentially
    very useful for symbolic expressions.

    EXAMPLES:

    To integrate the function $x^2$ from 0 to 1, we do ::

        sage: numerical_integral(x^2, 0, 1, max_points=100)
        (0.3333333333333333, 3.700743415417188e-15)

    To integrate the function $\sin(x)^3 + \sin(x)$ we do ::

        sage: numerical_integral(sin(x)^3 + sin(x),  0, pi)
        (3.333333333333333, 3.700743415417188e-14)

    The input can be any callable::

        sage: numerical_integral(lambda x: sin(x)^3 + sin(x),  0, pi)
        (3.333333333333333, 3.700743415417188e-14)

    We check this with a symbolic integration::

        sage: (sin(x)^3+sin(x)).integral(x,0,pi)
        10/3

    If we want to change the error tolerances and gauss rule used::

        sage: f = x^2
        sage: numerical_integral(f, 0, 1, max_points=200, eps_abs=1e-7, eps_rel=1e-7, rule=4)
        (0.3333333333333333, 3.700743415417188e-15)

    For a Python function with parameters::

        sage: f(x,a) = 1/(a+x^2)
        sage: [numerical_integral(f, 1, 2, max_points=100, params=[n]) for n in range(10)]  # random output (architecture and os dependent)
        [(0.49999999999998657, 5.5511151231256336e-15),
         (0.32175055439664557, 3.5721487367706477e-15),
         (0.24030098317249229, 2.6678768435816325e-15),
         (0.19253082576711697, 2.1375215571674764e-15),
         (0.16087527719832367, 1.7860743683853337e-15),
         (0.13827545676349412, 1.5351659583939151e-15),
         (0.12129975935702741, 1.3466978571966261e-15),
         (0.10806674191683065, 1.1997818507228991e-15),
         (0.09745444625548845, 1.0819617008493815e-15),
         (0.088750683050217577, 9.8533051773561173e-16)]
        sage: y = var('y')
        sage: numerical_integral(x*y, 0, 1)
        Traceback (most recent call last):
        ...
        ValueError: The function to be integrated depends on 2 variables (x, y),
        and so cannot be integrated in one dimension. Please fix additional
        variables with the 'params' argument

    Note the parameters are always a tuple even if they have one component.

    It is possible to integrate on infinite intervals as well by using
    +Infinity or -Infinity in the interval argument. For example::

        sage: f = exp(-x)
        sage: numerical_integral(f, 0, +Infinity)  # random output
        (0.99999999999957279, 1.8429811298996553e-07)

    Note the coercion to the real field RR, which prevents underflow::

        sage: f = exp(-x**2)
        sage: numerical_integral(f, -Infinity, +Infinity)  # random output
        (1.7724538509060035, 3.4295192165889879e-08)

    One can integrate any real-valued callable function::

        sage: numerical_integral(lambda x: abs(zeta(x)), [1.1,1.5])  # random output
        (1.8488570602160455, 2.052643677492633e-14)

    We can also numerically integrate symbolic expressions using either this
    function (which uses GSL) or the native integration (which uses Maxima)::

        sage: exp(-1/x).nintegral(x, 1, 2)  # via maxima
        (0.50479221787318..., 5.60431942934407...e-15, 21, 0)
        sage: numerical_integral(exp(-1/x), 1, 2)
        (0.50479221787318..., 5.60431942934407...e-15)

    We can also integrate constant expressions::

        sage: numerical_integral(2, 1, 7)
        (12.0, 0.0)

    If the interval of integration is a point, then the result is
    always zero (this makes sense within the Lebesgue theory of
    integration), see :trac:`12047`::

        sage: numerical_integral(log, 0, 0)
        (0.0, 0.0)
        sage: numerical_integral(lambda x: sqrt(x), (-2.0, -2.0) )
        (0.0, 0.0)

    In the presence of integrable singularity, the default adaptive method might
    fail and it is advised to use ``'qags'``::

        sage: b = 1.81759643554688
        sage: F(x) = sqrt((-x + b)/((x - 1.0)*x))
        sage: numerical_integral(F, 1, b)
        (inf, nan)
        sage: numerical_integral(F, 1, b, algorithm='qags')    # abs tol 1e-10
        (1.1817104238446596, 3.387268288079781e-07)

    AUTHORS:

    - Josh Kantor
    - William Stein
    - Robert Bradshaw
    - Jeroen Demeyer

    ALGORITHM: Uses calls to the GSL (GNU Scientific Library) C library.

    TESTS:

    Make sure that constant Expressions, not merely uncallable arguments,
    can be integrated (:trac:`10088`), at least if we can coerce them
    to float::

        sage: f, g = x, x-1
        sage: numerical_integral(f-g, -2, 2)
        (4.0, 0.0)
        sage: numerical_integral(SR(2.5), 5, 20)
        (37.5, 0.0)
        sage: numerical_integral(SR(1+3j), 2, 3)
        Traceback (most recent call last):
        ...
        TypeError: unable to simplify to float approximation
   """

   import inspect
   cdef double abs_err # step size
   cdef double result
   cdef int i
   cdef int j
   cdef double _a, _b
   cdef PyFunctionWrapper wrapper  # struct to pass information into GSL C function

   if b is None or isinstance(a, (list, tuple)):
      b = a[1]
      a = a[0]

   # The integral over a point is always zero
   if a == b:
       return (0.0, 0.0)

   if not callable(func):
        # handle the constant case
        return (((<double>b - <double>a) * <double>func), 0.0)

   cdef gsl_function F
   cdef gsl_integration_workspace* W
   W=NULL

   if not isinstance(func, FastDoubleFunc):
        try:
            if hasattr(func, 'arguments'):
                vars = func.arguments()
            else:
                vars = func.variables()
            if len(vars) == 0:
               # handle the constant case
               return (((<double>b - <double>a) * <double>func), 0.0)
            if len(vars) != 1:
                if len(params) + 1 != len(vars):
                   raise ValueError(("The function to be integrated depends on "
                                     "{} variables {}, and so cannot be "
                                     "integrated in one dimension. Please fix "
                                     "additional variables with the 'params' "
                                     "argument").format(len(vars),tuple(vars)))

                to_sub = dict(zip(vars[1:], params))
                func = func.subs(to_sub)
            func = func._fast_float_(str(vars[0]))
        except (AttributeError):
            pass

   if isinstance(func, FastDoubleFunc):
        F.function = c_ff
        F.params = <void *>func

   elif not isinstance(func, compiled_integrand):
      wrapper = PyFunctionWrapper()
      if not func is None:
         wrapper.the_function = func
      else:
         raise ValueError("No integrand defined")
      try:
         if params==[] and len(inspect.getargspec(wrapper.the_function)[0])==1:
            wrapper.the_parameters=[]
         elif params==[] and len(inspect.getargspec(wrapper.the_function)[0])>1:
            raise ValueError("Integrand has parameters but no parameters specified")
         elif params!=[]:
            wrapper.the_parameters = params
      except TypeError:
         wrapper.the_function = eval("lambda x: func(x)", {'func':func})
         wrapper.the_parameters = []

      F.function=c_f
      F.params=<void *> wrapper


   cdef size_t n
   n=max_points

   gsl_set_error_handler_off()

   if algorithm=="qng":
      _a=a
      _b=b
      sig_on()
      gsl_integration_qng(&F,_a,_b,eps_abs,eps_rel,&result,&abs_err,&n)
      sig_off()

   elif algorithm=="qag":
      from sage.rings.infinity import Infinity
      if a is -Infinity and b is +Infinity:
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         sig_on()
         gsl_integration_qagi(&F,eps_abs,eps_rel,n,W,&result,&abs_err)
         sig_off()

      elif a is -Infinity:
         _b=b
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         sig_on()
         gsl_integration_qagil(&F,_b,eps_abs,eps_rel,n,W,&result,&abs_err)
         sig_off()

      elif b is +Infinity:
         _a=a
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         sig_on()
         gsl_integration_qagiu(&F,_a,eps_abs,eps_rel,n,W,&result,&abs_err)
         sig_off()

      else:
         _a=a
         _b=b
         W = <gsl_integration_workspace*> gsl_integration_workspace_alloc(n)
         sig_on()
         gsl_integration_qag(&F,_a,_b,eps_abs,eps_rel,n,rule,W,&result,&abs_err)
         sig_off()


   elif algorithm == "qags":

        W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
        sig_on()
        _a=a
        _b=b
        gsl_integration_qags(&F,_a,_b,eps_abs,eps_rel,n,W,&result,&abs_err)
        sig_off()

   else:
      raise TypeError("invalid integration algorithm")

   if W != NULL:
      gsl_integration_workspace_free(W)

   return result, abs_err
