# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.symbolic
"""
Numerical Integration

AUTHORS:

- Josh Kantor (2007-02): first version

- William Stein (2007-02): rewrite of docs, conventions, etc.

- Robert Bradshaw (2008-08): fast float integration

- Jeroen Demeyer (2011-11-23): :issue:`12047`: return 0 when the
  integration interval is a point; reformat documentation and add to
  the reference manual.
"""

# ****************************************************************************
#       Copyright (C) 2004,2005,2006,2007 Joshua Kantor <kantor.jm@gmail.com>
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2019 Vincent Klein <vinklein@gmail.com>
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator

from sage.rings.real_double import RDF
from sage.libs.gsl.all cimport *
from sage.misc.sageinspect import sage_getargspec
from sage.ext.interpreters.wrapper_rdf cimport Wrapper_rdf
from sage.ext.fast_callable import fast_callable


cdef class PyFunctionWrapper:
    cdef object the_function
    cdef object the_parameters
    cdef list lx

cdef class compiled_integrand:
    cdef int c_f(self, double t) noexcept:  # void *params):
        return 0

cdef double c_f(double t, void *params) noexcept:
    cdef double value
    cdef PyFunctionWrapper wrapper
    wrapper = <PyFunctionWrapper> params
    try:
        if len(wrapper.the_parameters) != 0:
            value = wrapper.the_function(t, wrapper.the_parameters)
        else:
            value = wrapper.the_function(t)
    except Exception as msg:
        print(msg)
        return 0

    return value


def numerical_integral(func, a, b=None,
                       algorithm='qag',
                       max_points=87, params=[], eps_abs=1e-6,
                       eps_rel=1e-6, rule=6):
    r"""
    Return the numerical integral of the function on the interval
    from a to b and an error bound.

    INPUT:

    - ``a``, ``b`` -- the interval of integration, specified as two
      numbers or as a tuple/list with the first element the lower bound
      and the second element the upper bound.  Use ``+Infinity`` and
      ``-Infinity`` for plus or minus infinity.
    - ``algorithm`` -- valid choices are:

      * 'qag' -- for an adaptive integration
      * 'qags' -- for an adaptive integration with (integrable) singularities
      * 'qng' -- for a non-adaptive Gauss-Kronrod (samples at a maximum of 87pts)

    - ``max_points`` -- sets the maximum number of sample points
    - ``params`` -- used to pass parameters to your function
    - ``eps_abs``, ``eps_rel`` -- sets the absolute and relative error
      tolerances which satisfies the relation ``|RESULT - I|  <= max(eps_abs,
      eps_rel * |I|)``, where ``I = \int_a^b f(x) d x``.
    - ``rule`` -- this controls the Gauss-Kronrod rule used in the adaptive integration:

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

    To integrate the function `x^2` from 0 to 1, we do ::

        sage: numerical_integral(x^2, 0, 1, max_points=100)
        (0.3333333333333333, 3.700743415417188e-15)

    To integrate the function `\sin(x)^3 + \sin(x)` we do ::

        sage: numerical_integral(sin(x)^3 + sin(x),  0, pi)
        (3.333333333333333, 3.700743415417188e-14)

    The input can be any callable::

        sage: numerical_integral(lambda x: sin(x)^3 + sin(x),  0, pi)
        (3.333333333333333, 3.700743415417188e-14)

    We check this with a symbolic integration::

        sage: (sin(x)^3+sin(x)).integral(x,0,pi)
        10/3

    If we want to change the error tolerances and Gauss rule used::

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
    integration), see :issue:`12047`::

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
    Documentation can be found in [GSL]_ chapter "Numerical Integration".

    TESTS:

    Make sure that constant Expressions, not merely uncallable arguments,
    can be integrated (:issue:`10088`), at least if we can coerce them
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

    Check for :issue:`15496`::

        sage: f = x^2/exp(-1/(x^2+1))/(x^2+1)
        sage: D = integrate(f,(x,-infinity,infinity),hold=True)
        sage: D.n()
        Traceback (most recent call last):
        ...
        ValueError: integral does not converge at -infinity

    Symbolic functions can be integrated as conveniently as symbolic
    expressions, as in :issue:`15219`::

        sage: h(x) = x
        sage: numerical_integral(h,0,1)[0] # abs tol 1e-8
        0.5
    """
    cdef double abs_err  # step size
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
    W = NULL

    if True:
        from sage.rings.infinity import Infinity
        try:
            if hasattr(func, 'arguments'):
                vars = func.arguments()
            else:
                vars = func.variables()
        except (AttributeError):
            pass
        else:
            if not vars:
                # handle the constant case
                return (((<double>b - <double>a) * <double>func), 0.0)
            if len(vars) != 1:
                if len(params) + 1 != len(vars):
                    raise ValueError(("The function to be integrated depends on "
                                      "{} variables {}, and so cannot be "
                                      "integrated in one dimension. Please fix "
                                      "additional variables with the 'params' "
                                      "argument").format(len(vars), tuple(vars)))

                to_sub = dict(zip(vars[1:], params))
                func = func.subs(to_sub)

            # sanity checks for integration up to infinity
            v = str(vars[0])
            if a is -Infinity:
                try:
                    ell = func.limit(**{v: -Infinity})
                except (AttributeError, ValueError):
                    pass
                else:
                    if ell.is_numeric() and not ell.is_zero():
                        raise ValueError('integral does not converge at -infinity')
            if b is Infinity:
                try:
                    ell = func.limit(**{v: Infinity})
                except (AttributeError, ValueError):
                    pass
                else:
                    if ell.is_numeric() and not ell.is_zero():
                        raise ValueError('integral does not converge at infinity')
            func = fast_callable(func, vars=[v], domain=float)

    if not isinstance(func, compiled_integrand):
        wrapper = PyFunctionWrapper()
        if func is not None:
            wrapper.the_function = func
        else:
            raise ValueError("No integrand defined")
        try:
            if not params and len(sage_getargspec(wrapper.the_function)[0]) == 1:
                wrapper.the_parameters = []
            elif not params and len(sage_getargspec(wrapper.the_function)[0]) > 1:
                raise ValueError("Integrand has parameters but no parameters specified")
            elif params:
                wrapper.the_parameters = params
        except TypeError:
            wrapper.the_function = eval("lambda x: func(x)", {'func': func})
            wrapper.the_parameters = []

        F.function = c_f
        F.params = <void *> wrapper

    cdef size_t n
    n = max_points

    gsl_set_error_handler_off()

    if algorithm == "qng":
        _a = a
        _b = b
        sig_on()
        gsl_integration_qng(&F, _a, _b, eps_abs, eps_rel, &result, &abs_err, &n)
        sig_off()

    elif algorithm == "qag":
        if a is -Infinity and b is +Infinity:
            W = <gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
            sig_on()
            gsl_integration_qagi(&F, eps_abs, eps_rel, n, W, &result, &abs_err)
            sig_off()

        elif a is -Infinity:
            _b = b
            W = <gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
            sig_on()
            gsl_integration_qagil(&F, _b, eps_abs, eps_rel, n, W, &result, &abs_err)
            sig_off()

        elif b is +Infinity:
            _a = a
            W = <gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
            sig_on()
            gsl_integration_qagiu(&F, _a, eps_abs, eps_rel, n, W, &result, &abs_err)
            sig_off()

        else:
            _a = a
            _b = b
            W = <gsl_integration_workspace*> gsl_integration_workspace_alloc(n)
            sig_on()
            gsl_integration_qag(&F,_a,_b,eps_abs,eps_rel,n,rule,W,&result,&abs_err)
            sig_off()

    elif algorithm == "qags":

        W = <gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
        sig_on()
        _a = a
        _b = b
        gsl_integration_qags(&F, _a, _b, eps_abs, eps_rel, n, W, &result, &abs_err)
        sig_off()

    else:
        raise TypeError("invalid integration algorithm")

    if W != NULL:
        gsl_integration_workspace_free(W)

    return result, abs_err


cdef double c_monte_carlo_f(double *t, size_t dim, void *params) noexcept:
    cdef double value
    cdef PyFunctionWrapper wrapper
    wrapper = <PyFunctionWrapper> params

    for i in range(dim):
        wrapper.lx[i] = t[i]

    try:
        if len(wrapper.the_parameters) != 0:
            value = wrapper.the_function(*wrapper.lx, *wrapper.the_parameters)
        else:
            value = wrapper.the_function(*wrapper.lx)
    except Exception as msg:
        print(msg)
        return 0

    return value


cdef double c_monte_carlo_ff(double *x, size_t dim, void *params) noexcept:
    cdef double result
    (<Wrapper_rdf> params).call_c(x, &result)
    return result


def monte_carlo_integral(func, xl, xu, size_t calls, algorithm='plain',
                         params=None):
    """
    Integrate ``func`` by Monte-Carlo method.

    Integrate ``func`` over the ``dim``-dimensional hypercubic region
    defined by the lower and upper limits in the arrays ``xl`` and
    ``xu``, each of size ``dim``.

    The integration uses a fixed number of function calls and obtains
    random sampling points using the default gsl's random number generator.

    ALGORITHM: Uses calls to the GSL (GNU Scientific Library) C library.
    Documentation can be found in [GSL]_ chapter "Monte Carlo Integration".

    INPUT:

    - ``func`` -- the function to integrate
    - ``params`` -- used to pass parameters to your function
    - ``xl`` -- list of lower limits
    - ``xu`` -- list of upper limits
    - ``calls`` -- number of functions calls used
    - ``algorithm`` -- valid choices are:

      * 'plain' -- The plain Monte Carlo algorithm samples points randomly
        from the integration region to estimate the integral and its error.
      * 'miser' -- The MISER algorithm of Press and Farrar is based on
        recursive stratified sampling
      * 'vegas' -- The VEGAS algorithm of Lepage is based on importance
        sampling.

    OUTPUT:

    A tuple whose first component is the approximated integral and whose second
    component is an error estimate.

    EXAMPLES::

        sage: x, y = SR.var('x,y')
        sage: monte_carlo_integral(x*y, [0,0], [2,2], 10000)   # abs tol 0.1
        (4.0, 0.0)
        sage: integral(integral(x*y, (x,0,2)), (y,0,2))
        4

    An example with a parameter::

        sage: x, y, z = SR.var('x,y,z')
        sage: monte_carlo_integral(x*y*z, [0,0], [2,2], 10000, params=[1.2])   # abs tol 0.1
        (4.8, 0.0)

    Integral of a constant::

        sage: monte_carlo_integral(3, [0,0], [2,2], 10000)   # abs tol 0.1
        (12, 0.0)

    Test different algorithms::

        sage: x, y, z = SR.var('x,y,z')
        sage: f(x,y,z) = exp(z) * cos(x + sin(y))
        sage: for algo in ['plain', 'miser', 'vegas']:  # abs tol 0.01
        ....:   monte_carlo_integral(f, [0,0,-1], [2,2,1], 10^6, algorithm=algo)
        (-1.06, 0.01)
        (-1.06, 0.01)
        (-1.06, 0.01)

    Tests with Python functions::

        sage: def f(u, v): return u * v
        sage: monte_carlo_integral(f, [0,0], [2,2], 10000)  # abs tol 0.1
        (4.0, 0.0)
        sage: monte_carlo_integral(lambda u,v: u*v, [0,0], [2,2], 10000)  # abs tol 0.1
        (4.0, 0.0)
        sage: def f(x1,x2,x3,x4): return x1*x2*x3*x4
        sage: monte_carlo_integral(f, [0,0], [2,2], 1000, params=[0.6,2])  # abs tol 0.2
        (4.8, 0.0)

    TESTS::

        sage: monte_carlo_integral(f, [0,0,0], [2,2], 10)
        Traceback (most recent call last):
        ...
        TypeError: xl and xu must be lists of floating point values of identical lengths
        sage: monte_carlo_integral(f, [0,0], [2,2], 1, algorithm='unicorn')
        Traceback (most recent call last):
        ...
        ValueError: 'unicorn' is an invalid value for algorithm
        sage: monte_carlo_integral(lambda x,y: y*x, [], [], 1)
        Traceback (most recent call last):
        ...
        NotImplementedError: 0 dimensional integration not available
        sage: monte_carlo_integral(x*y, [0,0,0], [2,2,2], 1)
        Traceback (most recent call last):
        ...
        ValueError: The function to be integrated depends on 2 variables (x, y),
        and so cannot be integrated in 3 dimensions. Please fix additional
        variables with the 'params' argument
        sage: def f(x,y): return x*y
        sage: monte_carlo_integral(f, [0,0,0], [2,2,2], 100)
        Traceback (most recent call last):
        ...
        ValueError: The function to be integrated depends on 2 variables ('x', 'y'),
        and so cannot be integrated in 3 dimensions. Please fix additional
        variables with the 'params' argument
        sage: monte_carlo_integral(x*y, [0], [2], 1)
        Traceback (most recent call last):
        ...
        ValueError: The function to be integrated depends on 2 variables (x, y),
        and so cannot be integrated in 1 dimensions. Please add more items in
        upper and lower limits
        sage: monte_carlo_integral(f, [0], [2], 100)
        Traceback (most recent call last):
        ...
        ValueError: The function to be integrated depends on 2 variables ('x', 'y'),
        and so cannot be integrated in 1 dimensions. Please add more items in
        upper and lower limits

    AUTHORS:

    - Vincent Delecroix
    - Vincent Klein
    """
    cdef double result
    cdef double abs_err
    cdef gsl_monte_function F
    cdef PyFunctionWrapper wrapper  # struct to pass information into GSL Monte C function
    cdef gsl_monte_plain_state* state_plain = NULL
    cdef gsl_monte_miser_state* state_miser = NULL
    cdef gsl_monte_vegas_state* state_vegas = NULL
    cdef gsl_rng_type *type_rng
    cdef gsl_rng *_rng
    cdef size_t dim
    cdef double *_xl
    cdef double *_xu
    cdef MemoryAllocator mem = MemoryAllocator()

    if not isinstance(xl, (tuple, list)) or \
       not isinstance(xu, (tuple, list)) or \
       len(xl) != len(xu):
        raise TypeError("xl and xu must be lists of floating point values of"
                        " identical lengths")

    if algorithm not in ('plain', 'miser', 'vegas'):
        raise ValueError("'{}' is an invalid value for algorithm".format(algorithm))

    dim = len(xl)
    if not dim:
        raise NotImplementedError("0 dimensional integration not available")

    if params is None:
        params = []

    # Initialize hypercubic region's lower and upper limits
    _xl = <double *> mem.calloc(dim, sizeof(double))
    _xu = <double *> mem.calloc(dim, sizeof(double))
    for i in range(dim):
        _xl[i] = <double> xl[i]
        _xu[i] = <double> xu[i]

    if not callable(func):
        # constant. Note that all Expression objects are callable.
        v = float(1)
        for i in range(dim):
            v *= _xu[i] - _xl[i]
        return (v * <double?> func, 0.0)

    elif not isinstance(func, Wrapper_rdf):
        # func is either an Expression or another callable.
        try:
            vars = func.arguments()
        except AttributeError:
            try:
                vars = func.variables()
            except AttributeError:
                vars = sage_getargspec(func)[0]

        target_dim = dim + len(params)
        if len(vars) < target_dim:
            raise ValueError(("The function to be integrated depends on "
                              "{} variables {}, and so cannot be "
                                 "integrated in {} dimensions. Please fix "
                              "additional variables with the 'params' "
                              "argument").format(len(vars), tuple(vars),
                                                 target_dim))
        elif len(vars) > target_dim:
            raise ValueError(("The function to be integrated depends on "
                              "{} variables {}, and so cannot be "
                              "integrated in {} dimensions. Please add "
                              "more items in upper and lower limits"
                             ).format(len(vars), tuple(vars), target_dim))

        from sage.structure.element import Expression
        if isinstance(func, Expression):
            if params:
                to_sub = dict(zip(vars[-len(params):], params))
                func = func.subs(to_sub)
                vars = vars[:dim]

            func = fast_callable(func, domain=RDF, vars=vars)

    if isinstance(func, Wrapper_rdf):
        F.dim = dim
        F.f = c_monte_carlo_ff
        F.params = <void *>func
    else:
        wrapper = PyFunctionWrapper()
        wrapper.the_function = func

        if not params and len(sage_getargspec(wrapper.the_function)[0]) == dim:
            wrapper.the_parameters = []
        elif not params and len(sage_getargspec(wrapper.the_function)[0]) > dim:
            raise ValueError("Integrand has parameters but no parameters specified")
        elif params:
            wrapper.the_parameters = params
        wrapper.lx = [None] * dim

        F.dim = dim
        F.f = c_monte_carlo_f
        F.params = <void *> wrapper

    try:
        # Initialize the random number generator
        gsl_rng_env_setup()
        type_rng = gsl_rng_default
        _rng = gsl_rng_alloc(type_rng)

        if algorithm == 'plain':
            state_plain = <gsl_monte_plain_state*> gsl_monte_plain_alloc(dim)
            sig_on()
            gsl_monte_plain_integrate(&F, _xl, _xu, dim, calls, _rng,
                                      state_plain, &result, &abs_err)
            sig_off()
        elif algorithm == 'miser':
            state_miser = <gsl_monte_miser_state*> gsl_monte_miser_alloc(dim)
            sig_on()
            gsl_monte_miser_integrate(&F, _xl, _xu, dim, calls, _rng,
                                      state_miser, &result, &abs_err)
            sig_off()
        elif algorithm == 'vegas':
            state_vegas = <gsl_monte_vegas_state*> gsl_monte_vegas_alloc(dim)
            sig_on()
            gsl_monte_vegas_integrate(&F, _xl, _xu, dim, calls, _rng,
                                      state_vegas, &result, &abs_err)
            sig_off()
    finally:
        gsl_rng_free(_rng)

        if state_plain != NULL:
            gsl_monte_plain_free(state_plain)
        elif state_miser != NULL:
            gsl_monte_miser_free(state_miser)
        elif state_vegas != NULL:
            gsl_monte_vegas_free(state_vegas)

    return result, abs_err
