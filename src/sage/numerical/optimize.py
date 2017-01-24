"""
Numerical Root Finding and Optimization

AUTHOR:

- William Stein (2007): initial version
- Nathann Cohen (2008) : Bin Packing


Functions and Methods
----------------------
"""
from six.moves import range

from sage.modules.free_module_element import vector
from sage.rings.real_double import RDF


def find_root(f, a, b, xtol=10e-13, rtol=4.5e-16, maxiter=100, full_output=False):
    """
    Numerically find a root of ``f`` on the closed interval `[a,b]`
    (or `[b,a]`) if possible, where ``f`` is a function in the one variable.
    Note: this function only works in fixed (machine) precision, it is not
    possible to get arbitrary precision approximations with it.

    INPUT:

    - ``f`` -- a function of one variable or symbolic equality

    - ``a``, ``b`` -- endpoints of the interval

    - ``xtol``, ``rtol`` -- the routine converges when a root is known
      to lie within ``xtol`` of the value return. Should be `\geq 0`.
      The routine modifies this to take into account the relative precision
      of doubles.

    - ``maxiter`` -- integer; if convergence is not achieved in
      ``maxiter`` iterations, an error is raised. Must be `\geq 0`.

    - ``full_output`` -- bool (default: ``False``), if ``True``, also return
      object that contains information about convergence.


    EXAMPLES:

    An example involving an algebraic polynomial function::

        sage: R.<x> = QQ[]
        sage: f = (x+17)*(x-3)*(x-1/8)^3
        sage: find_root(f, 0,4)
        2.999999999999995
        sage: find_root(f, 0,1)  # note -- precision of answer isn't very good on some machines.
        0.124999...
        sage: find_root(f, -20,-10)
        -17.0

    In Pomerance's book on primes he asserts that the famous Riemann
    Hypothesis is equivalent to the statement that the function `f(x)`
    defined below is positive for all `x \geq 2.01`::

        sage: def f(x):
        ...       return sqrt(x) * log(x) - abs(Li(x) - prime_pi(x))

    We find where `f` equals, i.e., what value that is slightly smaller
    than `2.01` that could have been used in the formulation of the Riemann
    Hypothesis::

        sage: find_root(f, 2, 4, rtol=0.0001)
        2.0082...

    This agrees with the plot::

        sage: plot(f,2,2.01)
        Graphics object consisting of 1 graphics primitive
    """
    try:
        return f.find_root(a=a,b=b,xtol=xtol,rtol=rtol,maxiter=maxiter,full_output=full_output)
    except AttributeError:
        pass
    a = float(a); b = float(b)
    if a > b:
        a, b = b, a
    left = f(a)
    right = f(b)
    if left > 0 and right > 0:
        # Refine further -- try to find a point where this
        # function is negative in the interval
        val, s = find_local_minimum(f, a, b)
        if val > 0:
            if val < rtol:
                if full_output:
                    return s, "No extra data"
                else:
                    return s
            raise RuntimeError("f appears to have no zero on the interval")
        # If we found such an s, then we just instead find
        # a root between left and s or s and right.
        a = s   # arbitrary choice -- maybe should try both and take one that works?

    elif left < 0 and right < 0:
        # Refine further
        val, s = find_local_maximum(f, a, b)
        if val < 0:
            if abs(val) < rtol:
                if full_output:
                    return s, "No extra data"
                else:
                    return s
            raise RuntimeError("f appears to have no zero on the interval")
        a = s

    import scipy.optimize
    return scipy.optimize.brentq(f, a, b,
                                 full_output=full_output, xtol=xtol, rtol=rtol, maxiter=maxiter)

def find_local_maximum(f, a, b, tol=1.48e-08, maxfun=500):
    """
    Numerically find a local maximum of the expression `f` on the interval
    `[a,b]` (or `[b,a]`) along with the point at which the maximum is attained.

    Note that this function only finds a *local* maximum, and not the
    global maximum on that interval -- see the examples with
    :func:`find_local_maximum`.

    See the documentation for :func:`find_local_maximum` for more
    details and possible workarounds for finding the global minimum on
    an interval.

    EXAMPLES::

        sage: f = lambda x: x*cos(x)
        sage: find_local_maximum(f, 0, 5)
        (0.561096338191..., 0.8603335890...)
        sage: find_local_maximum(f, 0, 5, tol=0.1, maxfun=10)
        (0.561090323458..., 0.857926501456...)
        sage: find_local_maximum(8*e^(-x)*sin(x) - 1, 0, 7)
        (1.579175535558..., 0.7853981...)
    """
    try:
        return f.find_local_maximum(a=a, b=b, tol=tol, maxfun=maxfun)
    except AttributeError:
        pass
    minval, x = find_local_minimum(lambda z: -f(z), a=a, b=b, tol=tol, maxfun=maxfun)
    return -minval, x

def find_local_minimum(f, a, b, tol=1.48e-08, maxfun=500):
    """
    Numerically find a local minimum of the expression ``f`` on the
    interval `[a,b]` (or `[b,a]`) and the point at which it attains that
    minimum.  Note that ``f`` must be a function of (at most) one
    variable.

    Note that this function only finds a *local* minimum, and not the
    global minimum on that interval -- see the examples below.

    INPUT:

    - ``f`` -- a function of at most one variable.

    - ``a``, ``b`` -- endpoints of interval on which to minimize self.

    - ``tol`` -- the convergence tolerance

    - ``maxfun`` -- maximum function evaluations


    OUTPUT:

    - ``minval`` -- (float) the minimum value that self takes on in the
      interval `[a,b]`

    - ``x`` -- (float) the point at which self takes on the minimum value


    EXAMPLES::

        sage: f = lambda x: x*cos(x)
        sage: find_local_minimum(f, 1, 5)
        (-3.28837139559..., 3.4256184695...)
        sage: find_local_minimum(f, 1, 5, tol=1e-3)
        (-3.28837136189098..., 3.42575079030572...)
        sage: find_local_minimum(f, 1, 5, tol=1e-2, maxfun=10)
        (-3.28837084598..., 3.4250840220...)
        sage: show(plot(f, 0, 20))
        sage: find_local_minimum(f, 1, 15)
        (-9.4772942594..., 9.5293344109...)

    Only local minima are found; if you enlarge the interval, the
    returned minimum may be *larger*! See :trac:`2607`.

    ::

        sage: f(x) = -x*sin(x^2)
        sage: find_local_minimum(f, -2.5, -1)
        (-2.182769784677722, -2.1945027498534686)

    Enlarging the interval returns a larger minimum::

        sage: find_local_minimum(f, -2.5, 2)
        (-1.3076194129914434, 1.3552111405712108)

    One work-around is to plot the function and grab the minimum from
    that, although the plotting code does not necessarily do careful
    numerics (observe the small number of decimal places that we
    actually test)::

        sage: plot(f, (x,-2.5, -1)).ymin()
        -2.1827...
        sage: plot(f, (x,-2.5, 2)).ymin()
        -2.1827...

    ALGORITHM:

    Uses `scipy.optimize.fminbound
    <http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html>`_
    which uses Brent's method.


    AUTHOR:

    - William Stein (2007-12-07)
    """
    try:
        return f.find_local_minimum(a=a, b=b, tol=tol, maxfun=maxfun)
    except AttributeError:
        pass
    a = float(a); b = float(b)
    import scipy.optimize
    xmin, fval, iter, funcalls = scipy.optimize.fminbound(f, a, b, full_output=1, xtol=tol, maxfun=maxfun)
    return fval, xmin

def minimize(func,x0,gradient=None,hessian=None,algorithm="default",**args):
    r"""
    This function is an interface to a variety of algorithms for computing
    the minimum of a function of several variables.


    INPUT:

    - ``func`` -- Either a symbolic function or a Python function whose
      argument is a tuple with `n` components

    - ``x0`` -- Initial point for finding minimum.

    - ``gradient`` -- Optional gradient function. This will be computed
      automatically for symbolic functions.  For Python functions, it allows
      the use of algorithms requiring derivatives.  It should accept a
      tuple of arguments and return a NumPy array containing the partial
      derivatives at that point.

    - ``hessian`` --  Optional hessian function. This will be computed
      automatically for symbolic functions. For Python functions, it allows
      the use of algorithms requiring derivatives. It should accept a tuple
      of arguments and return a NumPy array containing the second partial
      derivatives of the function.

    - ``algorithm`` -- String specifying algorithm to use. Options are
      ``'default'`` (for Python functions, the simplex method is the default)
      (for symbolic functions bfgs is the default):

       - ``'simplex'``

       - ``'powell'``

       - ``'bfgs'`` -- (Broyden-Fletcher-Goldfarb-Shanno) requires
         ``gradient``

       - ``'cg'`` -- (conjugate-gradient) requires gradient

       - ``'ncg'`` -- (newton-conjugate gradient) requires gradient and hessian


    EXAMPLES::

        sage: vars=var('x y z')
        sage: f=100*(y-x^2)^2+(1-x)^2+100*(z-y^2)^2+(1-y)^2
        sage: minimize(f,[.1,.3,.4],disp=0)
        (1.00..., 1.00..., 1.00...)

        sage: minimize(f,[.1,.3,.4],algorithm="ncg",disp=0)
        (0.9999999..., 0.999999..., 0.999999...)

    Same example with just Python functions::

        sage: def rosen(x): # The Rosenbrock function
        ...      return sum(100.0r*(x[1r:]-x[:-1r]**2.0r)**2.0r + (1r-x[:-1r])**2.0r)
        sage: minimize(rosen,[.1,.3,.4],disp=0)
        (1.00..., 1.00..., 1.00...)

    Same example with a pure Python function and a Python function to
    compute the gradient::

        sage: def rosen(x): # The Rosenbrock function
        ...      return sum(100.0r*(x[1r:]-x[:-1r]**2.0r)**2.0r + (1r-x[:-1r])**2.0r)
        sage: import numpy
        sage: from numpy import zeros
        sage: def rosen_der(x):
        ...      xm = x[1r:-1r]
        ...      xm_m1 = x[:-2r]
        ...      xm_p1 = x[2r:]
        ...      der = zeros(x.shape,dtype=float)
        ...      der[1r:-1r] = 200r*(xm-xm_m1**2r) - 400r*(xm_p1 - xm**2r)*xm - 2r*(1r-xm)
        ...      der[0] = -400r*x[0r]*(x[1r]-x[0r]**2r) - 2r*(1r-x[0])
        ...      der[-1] = 200r*(x[-1r]-x[-2r]**2r)
        ...      return der
        sage: minimize(rosen,[.1,.3,.4],gradient=rosen_der,algorithm="bfgs",disp=0)
        (1.00...,  1.00..., 1.00...)
    """
    from sage.symbolic.expression import Expression
    from sage.ext.fast_eval import fast_callable
    import scipy
    from scipy import optimize
    if isinstance(func, Expression):
        var_list=func.variables()
        var_names = [str(_) for _ in var_list]
        fast_f=fast_callable(func, vars=var_names, domain=float)
        f=lambda p: fast_f(*p)
        gradient_list=func.gradient()
        fast_gradient_functions=[fast_callable(gradient_list[i], vars=var_names, domain=float)  for i in range(len(gradient_list))]
        gradient=lambda p: scipy.array([ a(*p) for a in fast_gradient_functions])
    else:
        f=func

    if algorithm=="default":
        if gradient is None:
            min = optimize.fmin(f, [float(_) for _ in x0], **args)
        else:
            min= optimize.fmin_bfgs(f, [float(_) for _ in x0],fprime=gradient, **args)
    else:
        if algorithm=="simplex":
            min= optimize.fmin(f, [float(_) for _ in x0], **args)
        elif algorithm=="bfgs":
            min= optimize.fmin_bfgs(f, [float(_) for _ in x0], fprime=gradient, **args)
        elif algorithm=="cg":
            min= optimize.fmin_cg(f, [float(_) for _ in x0], fprime=gradient, **args)
        elif algorithm=="powell":
            min= optimize.fmin_powell(f, [float(_) for _ in x0], **args)
        elif algorithm=="ncg":
            if isinstance(func, Expression):
                hess=func.hessian()
                hess_fast= [ [fast_callable(a, vars=var_names, domain=float) for a in row] for row in hess]
                hessian=lambda p: [[a(*p) for a in row] for row in hess_fast]
                hessian_p=lambda p,v: scipy.dot(scipy.array(hessian(p)),v)
                min= optimize.fmin_ncg(f, [float(_) for _ in x0], fprime=gradient, fhess=hessian, fhess_p=hessian_p, **args)
    return vector(RDF,min)

def minimize_constrained(func,cons,x0,gradient=None,algorithm='default', **args):
    r"""
    Minimize a function with constraints.


    INPUT:

    - ``func`` -- Either a symbolic function, or a Python function whose
      argument is a tuple with n components

    - ``cons`` -- constraints. This should be either a function or list of
      functions that must be positive. Alternatively, the constraints can
      be specified as a list of intervals that define the region we are
      minimizing in. If the constraints are specified as functions, the
      functions should be functions of a tuple with `n` components
      (assuming `n` variables). If the constraints are specified as a list
      of intervals and there are no constraints for a given variable, that
      component can be (``None``, ``None``).

    - ``x0`` -- Initial point for finding minimum

    - ``algorithm`` -- Optional, specify the algorithm to use:

      - ``'default'``  -- default choices

      - ``'l-bfgs-b'`` -- only effective if you specify bound constraints.
        See [ZBN1997]_.

    - ``gradient`` -- Optional gradient function. This will be computed
      automatically for symbolic functions. This is only used when the
      constraints are specified as a list of intervals.


    EXAMPLES:

    Let us maximize `x + y - 50` subject to the following constraints:
    `50x + 24y \leq 2400`, `30x + 33y \leq 2100`, `x \geq 45`,
    and `y \geq 5`::

        sage: y = var('y')
        sage: f = lambda p: -p[0]-p[1]+50
        sage: c_1 = lambda p: p[0]-45
        sage: c_2 = lambda p: p[1]-5
        sage: c_3 = lambda p: -50*p[0]-24*p[1]+2400
        sage: c_4 = lambda p: -30*p[0]-33*p[1]+2100
        sage: a = minimize_constrained(f,[c_1,c_2,c_3,c_4],[2,3])
        sage: a
        (45.0, 6.25...)

    Let's find a minimum of `\sin(xy)`::

        sage: x,y = var('x y')
        sage: f = sin(x*y)
        sage: minimize_constrained(f, [(None,None),(4,10)],[5,5])
        (4.8..., 4.8...)

    Check, if L-BFGS-B finds the same minimum::

        sage: minimize_constrained(f, [(None,None),(4,10)],[5,5], algorithm='l-bfgs-b')
        (4.7..., 4.9...)

    Rosenbrock function, [http://en.wikipedia.org/wiki/Rosenbrock_function]::

        sage: from scipy.optimize import rosen, rosen_der
        sage: minimize_constrained(rosen, [(-50,-10),(5,10)],[1,1],gradient=rosen_der,algorithm='l-bfgs-b')
        (-10.0, 10.0)
        sage: minimize_constrained(rosen, [(-50,-10),(5,10)],[1,1],algorithm='l-bfgs-b')
        (-10.0, 10.0)
    """
    from sage.symbolic.expression import Expression
    import scipy
    from scipy import optimize
    function_type=type(lambda x,y: x+y)

    if isinstance(func, Expression):
        var_list=func.variables()
        var_names = [str(_) for _ in var_list]
        fast_f=func._fast_float_(*var_names)
        f=lambda p: fast_f(*p)
        gradient_list=func.gradient()
        fast_gradient_functions=[gradient_list[i]._fast_float_(*var_names)  for i in range(len(gradient_list))]
        gradient=lambda p: scipy.array([ a(*p) for a in fast_gradient_functions])
    else:
        f=func

    if isinstance(cons,list):
        if isinstance(cons[0],tuple) or isinstance(cons[0],list) or cons[0] is None:
            if gradient is not None:
                if algorithm=='l-bfgs-b':
                    min= optimize.fmin_l_bfgs_b(f,x0,gradient,bounds=cons, iprint=-1, **args)[0]
                else:
                    min= optimize.fmin_tnc(f,x0,gradient,bounds=cons,messages=0,**args)[0]
            else:
                if algorithm=='l-bfgs-b':
                    min= optimize.fmin_l_bfgs_b(f,x0,approx_grad=True,bounds=cons,iprint=-1, **args)[0]
                else:
                    min= optimize.fmin_tnc(f,x0,approx_grad=True,bounds=cons,messages=0,**args)[0]

        elif isinstance(cons[0],function_type):
            min= optimize.fmin_cobyla(f,x0,cons,iprint=0,**args)
    elif isinstance(cons, function_type):
        min= optimize.fmin_cobyla(f,x0,cons,iprint=0,**args)
    return vector(RDF,min)


def linear_program(c,G,h,A=None,b=None,solver=None):
    """
    Solves the dual linear programs:

    - Minimize  `c'x` subject to `Gx + s = h`, `Ax = b`, and `s \geq 0` where
      `'` denotes transpose.

    - Maximize  `-h'z - b'y` subject to `G'z + A'y + c = 0` and `z \geq 0`.


    INPUT:

    - ``c`` -- a vector

    - ``G`` -- a matrix

    - ``h`` -- a vector

    - ``A`` -- a matrix

    - ``b`` --- a vector

    - ``solver`` (optional) --- solver to use. If None, the cvxopt's lp-solver
                                is used. If it is 'glpk', then glpk's solver
                                is used.

    These can be over any field that can be turned into a floating point
    number.


    OUTPUT:

    A dictionary ``sol`` with keys ``x``, ``s``, ``y``, ``z`` corresponding
    to the variables above:

    - ``sol['x']`` -- the solution to the linear program

    - ``sol['s']`` -- the slack variables for the solution

    - ``sol['z']``, ``sol['y']`` -- solutions to the dual program


    EXAMPLES:

    First, we minimize `-4x_1 - 5x_2` subject to `2x_1 + x_2 \leq 3`,
    `x_1 +  2x_2 \leq 3`, `x_1 \geq 0`, and `x_2 \geq 0`::

        sage: c=vector(RDF,[-4,-5])
        sage: G=matrix(RDF,[[2,1],[1,2],[-1,0],[0,-1]])
        sage: h=vector(RDF,[3,3,0,0])
        sage: sol=linear_program(c,G,h)
        sage: sol['x']
        (0.999..., 1.000...)

    Next, we maximize `x+y-50` subject to `50x + 24y \leq 2400`,
    `30x + 33y \leq 2100`, `x \geq 45`, and `y \geq 5`::

        sage: v=vector([-1.0,-1.0,-1.0])
        sage: m=matrix([[50.0,24.0,0.0],[30.0,33.0,0.0],[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0],[0.0,0.0,-1.0]])
        sage: h=vector([2400.0,2100.0,-45.0,-5.0,1.0,-1.0])
        sage: sol=linear_program(v,m,h)
        sage: sol['x']
        (45.000000..., 6.2499999..., 1.00000000...)
        sage: sol=linear_program(v,m,h,solver='glpk')
        GLPK Simplex Optimizer...
        OPTIMAL LP SOLUTION FOUND
        sage: sol['x']
        (45.0..., 6.25..., 1.0...)
    """
    from cvxopt.base import matrix as m
    from cvxopt import solvers
    solvers.options['show_progress']=False
    if solver=='glpk':
        from cvxopt import glpk
        glpk.options['LPX_K_MSGLEV'] = 0
    c_=m(c.base_extend(RDF).numpy())
    G_=m(G.base_extend(RDF).numpy())
    h_=m(h.base_extend(RDF).numpy())
    if A is not None and b is not None:
        A_=m(A.base_extend(RDF).numpy())
        b_=m(b.base_extend(RDF).numpy())
        sol=solvers.lp(c_,G_,h_,A_,b_,solver=solver)
    else:
        sol=solvers.lp(c_,G_,h_,solver=solver)
    status=sol['status']
    if status != 'optimal':
       return  {'primal objective':None,'x':None,'s':None,'y':None,
              'z':None,'status':status}
    x=vector(RDF,list(sol['x']))
    s=vector(RDF,list(sol['s']))
    y=vector(RDF,list(sol['y']))
    z=vector(RDF,list(sol['z']))
    return  {'primal objective':sol['primal objective'],'x':x,'s':s,'y':y,
               'z':z,'status':status}


def find_fit(data, model, initial_guess = None, parameters = None, variables = None, solution_dict = False):
    r"""
    Finds numerical estimates for the parameters of the function model to
    give a best fit to data.


    INPUT:

    - ``data`` -- A two dimensional table of floating point numbers of the
      form `[[x_{1,1}, x_{1,2}, \ldots, x_{1,k}, f_1],
      [x_{2,1}, x_{2,2}, \ldots, x_{2,k}, f_2],
      \ldots,
      [x_{n,1}, x_{n,2}, \ldots, x_{n,k}, f_n]]` given as either a list of
      lists, matrix, or numpy array.

    - ``model`` -- Either a symbolic expression, symbolic function, or a
      Python function. ``model`` has to be a function of the variables
      `(x_1, x_2, \ldots, x_k)` and free parameters
      `(a_1, a_2, \ldots, a_l)`.

    - ``initial_guess`` -- (default: ``None``) Initial estimate for the
      parameters `(a_1, a_2, \ldots, a_l)`, given as either a list, tuple,
      vector or numpy array. If ``None``, the default estimate for each
      parameter is `1`.

    - ``parameters`` -- (default: ``None``) A list of the parameters
      `(a_1, a_2, \ldots, a_l)`. If model is a symbolic function it is
      ignored, and the free parameters of the symbolic function are used.

    - ``variables`` -- (default: ``None``) A list of the variables
      `(x_1, x_2, \ldots, x_k)`. If model is a symbolic function it is
      ignored, and the variables of the symbolic function are used.

    - ``solution_dict`` -- (default: ``False``) if ``True``, return the
      solution as a dictionary rather than an equation.


    EXAMPLES:

    First we create some data points of a sine function with some random
    perturbations::

        sage: data = [(i, 1.2 * sin(0.5*i-0.2) + 0.1 * normalvariate(0, 1)) for i in xsrange(0, 4*pi, 0.2)]
        sage: var('a, b, c, x')
        (a, b, c, x)

    We define a function with free parameters `a`, `b` and `c`::

        sage: model(x) = a * sin(b * x - c)

    We search for the parameters that give the best fit to the data::

        sage: find_fit(data, model)
        [a == 1.21..., b == 0.49..., c == 0.19...]

    We can also use a Python function for the model::

        sage: def f(x, a, b, c): return a * sin(b * x - c)
        sage: fit = find_fit(data, f, parameters = [a, b, c], variables = [x], solution_dict = True)
        sage: fit[a], fit[b], fit[c]
        (1.21..., 0.49..., 0.19...)

    We search for a formula for the `n`-th prime number::

        sage: dataprime = [(i, nth_prime(i)) for i in range(1, 5000, 100)]
        sage: find_fit(dataprime, a * x * log(b * x), parameters = [a, b], variables = [x])
        [a == 1.11..., b == 1.24...]


    ALGORITHM:

    Uses ``scipy.optimize.leastsq`` which in turn uses MINPACK's lmdif and
    lmder algorithms.
    """
    import numpy

    if not isinstance(data, numpy.ndarray):
        try:
            data = numpy.array(data, dtype = float)
        except (ValueError, TypeError):
            raise TypeError("data has to be a list of lists, a matrix, or a numpy array")
    elif data.dtype == object:
        raise ValueError("the entries of data have to be of type float")

    if data.ndim != 2:
        raise ValueError("data has to be a two dimensional table of floating point numbers")

    from sage.symbolic.expression import Expression

    if isinstance(model, Expression):
        if variables is None:
            variables = list(model.arguments())
        if parameters is None:
            parameters = list(model.variables())
            for v in variables:
                parameters.remove(v)

    if data.shape[1] != len(variables) + 1:
        raise ValueError("each row of data needs %d entries, only %d entries given" % (len(variables) + 1, data.shape[1]))

    if parameters is None or len(parameters) == 0 or \
       variables is None or len(variables) == 0:
        raise ValueError("no variables given")

    if initial_guess is None:
        initial_guess = len(parameters) * [1]

    if not isinstance(initial_guess, numpy.ndarray):
        try:
            initial_guess = numpy.array(initial_guess, dtype = float)
        except (ValueError, TypeError):
            raise TypeError("initial_guess has to be a list, tuple, or numpy array")
    elif initial_guess.dtype == object:
        raise ValueError("the entries of initial_guess have to be of type float")

    if len(initial_guess) != len(parameters):
        raise ValueError("length of initial_guess does not coincide with the number of parameters")

    if isinstance(model, Expression):
        var_list = variables + parameters
        var_names = [str(_) for _ in var_list]
        func = model._fast_float_(*var_names)
    else:
        func = model

    def function(x_data, params):
        result = numpy.zeros(len(x_data))
        for row in range(len(x_data)):
            fparams = numpy.hstack((x_data[row], params)).tolist()
            result[row] = func(*fparams)
        return result

    def error_function(params, x_data, y_data):
        result = numpy.zeros(len(x_data))
        for row in range(len(x_data)):
            fparams = x_data[row].tolist() + params.tolist()
            result[row] = func(*fparams)
        return result - y_data

    x_data = data[:, 0:len(variables)]
    y_data = data[:, -1]

    from scipy.optimize import leastsq
    estimated_params, d = leastsq(error_function, initial_guess, args = (x_data, y_data))

    if isinstance(estimated_params, float):
        estimated_params = [estimated_params]
    else:
        estimated_params = estimated_params.tolist()

    if solution_dict:
       dict = {}
       for item in zip(parameters, estimated_params):
           dict[item[0]] = item[1]
       return dict

    return [item[0] == item[1] for item in zip(parameters, estimated_params)]

def binpacking(items,maximum=1,k=None):
    r"""
    Solves the bin packing problem.

    The Bin Packing problem is the following :

    Given a list of items of weights `p_i` and a real value `K`, what is
    the least number of bins such that all the items can be put in the
    bins, while keeping sure that each bin contains a weight of at most `K` ?

    For more informations : http://en.wikipedia.org/wiki/Bin_packing_problem

    Two version of this problem are solved by this algorithm :
         * Is it possible to put the given items in `L` bins ?
         * What is the assignment of items using the
           least number of bins with the given list of items ?

    INPUT:

    - ``items`` -- A list of real values (the items' weight)

    - ``maximum``   -- The maximal size of a bin

    - ``k``     -- Number of bins

      - When set to an integer value, the function returns a partition
        of the items into `k` bins if possible, and raises an
        exception otherwise.

      - When set to ``None``, the function returns a partition of the items
        using the least number possible of bins.

    OUTPUT:

    A list of lists, each member corresponding to a box and containing
    the list of the weights inside it. If there is no solution, an
    exception is raised (this can only happen when ``k`` is specified
    or if ``maximum`` is less that the size of one item).

    EXAMPLES:

    Trying to find the minimum amount of boxes for 5 items of weights
    `1/5, 1/4, 2/3, 3/4, 5/7`::

        sage: from sage.numerical.optimize import binpacking
        sage: values = [1/5, 1/3, 2/3, 3/4, 5/7]
        sage: bins = binpacking(values)
        sage: len(bins)
        3

    Checking the bins are of correct size ::

        sage: all([ sum(b)<= 1 for b in bins ])
        True

    Checking every item is in a bin ::

        sage: b1, b2, b3 = bins
        sage: all([ (v in b1 or v in b2 or v in b3) for v in values ])
        True

    One way to use only three boxes (which is best possible) is to put
    `1/5 + 3/4` together in a box, `1/3+2/3` in another, and `5/7`
    by itself in the third one.

    Of course, we can also check that there is no solution using only two boxes ::

        sage: from sage.numerical.optimize import binpacking
        sage: binpacking([0.2,0.3,0.8,0.9], k=2)
        Traceback (most recent call last):
        ...
        ValueError: This problem has no solution !
    """

    if max(items) > maximum:
        raise ValueError("This problem has no solution !")

    if k is None:
        from sage.functions.other import ceil
        k=ceil(sum(items)/maximum)
        while True:
            from sage.numerical.mip import MIPSolverException
            try:
                return binpacking(items,k=k,maximum=maximum)
            except MIPSolverException:
                k = k + 1

    from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
    p=MixedIntegerLinearProgram()

    # Boolean variable indicating whether
    # the i th element belongs to box b
    box=p.new_variable(binary = True)

    # Each bin contains at most max
    for b in range(k):
        p.add_constraint(p.sum([items[i]*box[i,b] for i in range(len(items))]) <= maximum)

    # Each item is assigned exactly one bin
    for i in range(len(items)):
        p.add_constraint(p.sum([box[i,b] for b in range(k)]) == 1)

    p.set_objective(None)

    try:
        p.solve()
    except MIPSolverException:
        raise ValueError("This problem has no solution !")

    box=p.get_values(box)

    boxes=[[] for i in range(k)]

    for (i,b),value in box.iteritems():
        if value == 1:
            boxes[b].append(items[i])

    return boxes



