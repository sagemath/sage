# sage_setup: distribution = sagemath-categories
r"""
Fast Numerical Evaluation

For many applications such as numerical integration, differential
equation approximation, plotting a 3d surface, optimization problems,
monte-carlo simulations, etc., one wishes to pass around and evaluate
a single algebraic expression many, many times at various floating
point values. Doing this via recursive calls over a python
representation of the object (even if Maxima or other outside packages
are not involved) is extremely inefficient.

The solution implemented in this module, by Robert Bradshaw (2008-10),
has been superseded by :func:`~sage.ext.fast_callable.fast_callable`.
All that remains here is a compatible interface function :func:`fast_float`.

AUTHORS:

- Robert Bradshaw (2008-10): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.ext.fast_callable import fast_callable, Wrapper


def fast_float(f, *vars, expect_one_var=False):
    """
    Try to create a function that evaluates f quickly using
    floating-point numbers, if possible.

    On failure, returns the input unchanged.

    This is an alternative interface to :func:`sage.ext.fast_callable.fast_callable`.
    :issue:`32268` proposes to deprecate this function.

    INPUT:

    - ``f`` -- an expression
    - ``vars`` -- the names of the arguments
    - ``expect_one_var`` -- don't give deprecation warning if ``vars`` is
      omitted, as long as expression has only one var

    EXAMPLES::

        sage: from sage.ext.fast_eval import fast_float
        sage: x,y = var('x,y')                                                          # needs sage.symbolic
        sage: f = fast_float(sqrt(x^2+y^2), 'x', 'y')                                   # needs sage.symbolic
        sage: f(3,4)                                                                    # needs sage.symbolic
        5.0

    Specifying the argument names is essential, as fast_float objects
    only distinguish between arguments by order. ::

        sage: # needs sage.symbolic
        sage: f = fast_float(x-y, 'x','y')
        sage: f(1,2)
        -1.0
        sage: f = fast_float(x-y, 'y','x')
        sage: f(1,2)
        1.0
    """
    if isinstance(f, (tuple, list)):
        return tuple([fast_float(x, *vars, expect_one_var=expect_one_var) for x in f])

    cdef int i
    for i from 0 <= i < len(vars):
        if not isinstance(vars[i], str):
            v = str(vars[i])
            # inexact generators display as 1.00..0*x
            if '*' in v:
                v = v[v.index('*')+1:]
            vars = vars[:i] + (v,) + vars[i+1:]

    try:
        return fast_callable(f, vars=vars, domain=float,
                             expect_one_var=expect_one_var)
    except AttributeError:
        pass

    try:
        from sage.symbolic.ring import SR
        return fast_float(SR(f), *vars)
    except TypeError:
        pass

    if f is None:
        raise TypeError("no way to make fast_float from None")

    return f


def is_fast_float(x):
    return isinstance(x, Wrapper)
