# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

r"""
Utility functions for making derivative() behave uniformly across Sage.

To use these functions:

1. attach the following method to your class::

        def _derivative(self, var=None):
            [ should differentiate wrt the single variable var and return
                result; var==None means attempt to differentiate wrt a 'default'
                variable. ]

2. from sage.misc.derivative import multi_derivative

3. add the following method to your class::

        def derivative(self, *args):
            return multi_derivative(self, args)

Then your object will support the standard parameter format for derivative().
For example::

    F.derivative():
        diff wrt. default variable (calls F._derivative(None))

    F.derivative(3):
        diff three times wrt default variable (calls F._derivative(None) three times)

    F.derivative(x):
        diff wrt x (calls F._derivative(x))

    F.derivative(1, x, z, 3, y, 2, 2):
        diff once wrt default variable, then once wrt x, then three times wrt z,
        then twice wrt y, then twice wrt default variable.

    F.derivative([None, x, z, z, z, y, y, None, None]):
        identical to previous example

For the precise specification see documentation for derivative_parse().

AUTHORS:

- David Harvey (2008-02)
"""

from sage.rings.integer cimport Integer


def derivative_parse(args):
    r"""
    Translates a sequence consisting of 'variables' and iteration counts into
    a single sequence of variables.

    INPUT:

    - ``args`` -- any iterable, interpreted as a sequence of 'variables' and
      iteration counts. An iteration count is any integer type (python int
      or Sage Integer). Iteration counts must be nonnegative. Any object
      which is not an integer is assumed to be a variable.

    OUTPUT:

        A sequence, the 'expanded' version of the input, defined as follows.
        Read the input from left to right. If you encounter a variable V
        followed by an iteration count N, then output N copies of V. If V
        is not followed by an iteration count, output a single copy of V.
        If you encounter an iteration count N (not attached to a preceding
        variable), then output N copies of None.

        Special case: if input is empty, output [None] (i.e. "differentiate
        once with respect to the default variable").

        Special case: if the input is a 1-tuple containing a single list,
        then the return value is simply that list.

    EXAMPLES::

        sage: x = var("x")                                                              # needs sage.symbolic
        sage: y = var("y")                                                              # needs sage.symbolic
        sage: from sage.misc.derivative import derivative_parse

    Differentiate twice with respect to x, then once with respect to y,
    then once with respect to x::

        sage: derivative_parse([x, 2, y, x])                                            # needs sage.symbolic
        [x, x, y, x]

    Differentiate twice with respect to x, then twice with respect to
    the 'default variable'::

        sage: derivative_parse([x, 2, 2])                                               # needs sage.symbolic
        [x, x, None, None]

    Special case with empty input list::

        sage: derivative_parse([])
        [None]

        sage: derivative_parse([-1])
        Traceback (most recent call last):
        ...
        ValueError: derivative counts must be nonnegative

    Special case with single list argument provided::

        sage: derivative_parse(([x, y], ))                                              # needs sage.symbolic
        [x, y]

    If only the count is supplied::

        sage: derivative_parse([0])
        []
        sage: derivative_parse([1])
        [None]
        sage: derivative_parse([2])
        [None, None]
        sage: derivative_parse([int(2)])
        [None, None]

    Various other cases::

        sage: # needs sage.symbolic
        sage: derivative_parse([x])
        [x]
        sage: derivative_parse([x, x])
        [x, x]
        sage: derivative_parse([x, 2])
        [x, x]
        sage: derivative_parse([x, 0])
        []
        sage: derivative_parse([x, y, x, 2, 2, y])
        [x, y, x, x, None, None, y]
    """
    if not args:
        return [None]
    if len(args) == 1 and isinstance(args[0], list):
        return args[0]

    output = []
    cdef bint got_var = 0   # have a variable saved up from last loop?
    cdef int count, i

    for arg in args:
        if isinstance(arg, (int, Integer)):
            # process iteration count
            count = int(arg)
            if count < 0:
                raise ValueError("derivative counts must be nonnegative")
            if not got_var:
                var = None
            for i from 0 <= i < count:
                output.append(var)
            got_var = 0
        else:
            # process variable
            if got_var:
                output.append(var)
            got_var = 1
            var = arg

    if got_var:
        output.append(var)

    return output


def multi_derivative(F, args):
    r"""
    Call F._derivative(var) for a sequence of variables specified by args.

    INPUT:

    - ``F`` -- any object with a ``_derivative(var)`` method
    - ``args`` -- any tuple that can be processed by :func:`derivative_parse`

    EXAMPLES::

        sage: from sage.misc.derivative import multi_derivative
        sage: R.<x, y, z> = PolynomialRing(QQ)
        sage: f = x^3 * y^4 * z^5
        sage: multi_derivative(f, (x,))   # like f.derivative(x)
        3*x^2*y^4*z^5
        sage: multi_derivative(f, (x, y, x))     # like f.derivative(x, y, x)
        24*x*y^3*z^5
        sage: multi_derivative(f, ([x, y, x],))    # like f.derivative([x, y, x])
        24*x*y^3*z^5
        sage: multi_derivative(f, (x, 2))     # like f.derivative(x, 2)
        6*x*y^4*z^5

    ::

        sage: R.<x> = PolynomialRing(QQ)
        sage: f = x^4 + x^2 + 1
        sage: multi_derivative(f, [])   # like f.derivative()
        4*x^3 + 2*x
        sage: multi_derivative(f, [[]])   # like f.derivative([])
        x^4 + x^2 + 1
        sage: multi_derivative(f, [x])     # like f.derivative(x)
        4*x^3 + 2*x
    """
    if not args:
        # fast version where no arguments supplied
        return F._derivative()

    for arg in derivative_parse(args):
        F = F._derivative(arg)
    return F
