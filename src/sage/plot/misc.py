# sage.doctest: needs sage.symbolic
"""
Plotting utilities
"""

# ****************************************************************************
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
from sage.ext.fast_callable import FastCallableFloatWrapper
from collections.abc import Iterable


def setup_for_eval_on_grid(funcs,
                           ranges,
                           plot_points=None,
                           return_vars=False,
                           imaginary_tolerance=1e-8):
    r"""
    Calculate the necessary parameters to construct a list of points,
    and make the functions fast_callable.

    INPUT:

    - ``funcs`` -- a function, or a list, tuple, or vector of functions

    - ``ranges`` -- a list of ranges.  A range can be a 2-tuple of
      numbers specifying the minimum and maximum, or a 3-tuple giving
      the variable explicitly.

    - ``plot_points`` -- a tuple of integers specifying the number of
      plot points for each range.  If a single number is specified, it
      will be the value for all ranges.  This defaults to 2.

    - ``return_vars`` -- (default ``False``) If ``True``, return the variables,
      in order.

    - ``imaginary_tolerance`` -- (default: ``1e-8``); if an imaginary
      number arises (due, for example, to numerical issues), this
      tolerance specifies how large it has to be in magnitude before
      we raise an error. In other words, imaginary parts smaller than
      this are ignored in your plot points.

    OUTPUT:


    - ``fast_funcs`` - if only one function passed, then a fast
      callable function.  If funcs is a list or tuple, then a tuple
      of fast callable functions is returned.

    - ``range_specs`` - a list of range_specs: for each range, a
      tuple is returned of the form (range_min, range_max,
      range_step) such that ``srange(range_min, range_max,
      range_step, include_endpoint=True)`` gives the correct points
      for evaluation.

    EXAMPLES::

        sage: x,y,z=var('x,y,z')
        sage: f(x,y)=x+y-z
        sage: g(x,y)=x+y
        sage: h(y)=-y
        sage: sage.plot.misc.setup_for_eval_on_grid(f, [(0, 2),(1,3),(-4,1)], plot_points=5)
        (<sage...>, [(0.0, 2.0, 0.5), (1.0, 3.0, 0.5), (-4.0, 1.0, 1.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid([g,h], [(0, 2),(-1,1)], plot_points=5)
        ((<sage...>, <sage...>), [(0.0, 2.0, 0.5), (-1.0, 1.0, 0.5)])
        sage: sage.plot.misc.setup_for_eval_on_grid([sin,cos], [(-1,1)], plot_points=9)
        ((<sage...>, <sage...>), [(-1.0, 1.0, 0.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid([lambda x: x^2,cos], [(-1,1)], plot_points=9)
        ((<function <lambda> ...>, <sage...>), [(-1.0, 1.0, 0.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid([x+y], [(x,-1,1),(y,-2,2)])
        ((<sage...>,), [(-1.0, 1.0, 2.0), (-2.0, 2.0, 4.0)])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,-1,1),(y,-1,1)], plot_points=[4,9])
        (<sage...>, [(-1.0, 1.0, 0.6666666666666666), (-1.0, 1.0, 0.25)])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,-1,1),(y,-1,1)], plot_points=[4,9,10])
        Traceback (most recent call last):
        ...
        ValueError: plot_points must be either an integer or a list of integers, one for each range
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(1,-1),(y,-1,1)], plot_points=[4,9,10])
        Traceback (most recent call last):
        ...
        ValueError: Some variable ranges specify variables while others do not

    Beware typos: a comma which should be a period, for instance::

        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x, 1, 2), (y, 0,1, 0.2)], plot_points=[4,9,10])
        Traceback (most recent call last):
        ...
        ValueError: At least one variable range has more than 3 entries: each should either have 2 or 3 entries, with one of the forms (xmin, xmax) or (x, xmin, xmax)

        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(y,1,-1),(x,-1,1)], plot_points=5)
        (<sage...>, [(-1.0, 1.0, 0.5), (-1.0, 1.0, 0.5)])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,1,-1),(x,-1,1)], plot_points=5)
        Traceback (most recent call last):
        ...
        ValueError: range variables should be distinct, but there are duplicates
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,1,1),(y,-1,1)])
        Traceback (most recent call last):
        ...
        ValueError: plot start point and end point must be different
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(x,1,-1),(y,-1,1)], return_vars=True)
        (<sage...>, [(-1.0, 1.0, 2.0), (-1.0, 1.0, 2.0)], [x, y])
        sage: sage.plot.misc.setup_for_eval_on_grid(x+y, [(y,1,-1),(x,-1,1)], return_vars=True)
        (<sage...>, [(-1.0, 1.0, 2.0), (-1.0, 1.0, 2.0)], [y, x])

    TESTS:

    Ensure that we can plot expressions with intermediate complex
    terms as in :trac:`8450`::

        sage: x, y = SR.var('x y')
        sage: contour_plot(abs(x+i*y), (x,-1,1), (y,-1,1))
        Graphics object consisting of 1 graphics primitive
        sage: density_plot(abs(x+i*y), (x,-1,1), (y,-1,1))
        Graphics object consisting of 1 graphics primitive
        sage: plot3d(abs(x+i*y), (x,-1,1),(y,-1,1))
        Graphics3d Object
        sage: streamline_plot(abs(x+i*y), (x,-1,1),(y,-1,1))
        Graphics object consisting of 1 graphics primitive

    """
    if max(map(len, ranges)) > 3:
        raise ValueError("At least one variable range has more than 3 entries: each should either have 2 or 3 entries, with one of the forms (xmin, xmax) or (x, xmin, xmax)")
    if max(map(len, ranges)) != min(map(len, ranges)):
        raise ValueError("Some variable ranges specify variables while others do not")

    if len(ranges[0]) == 3:
        vars = [r[0] for r in ranges]
        ranges = [r[1:] for r in ranges]
        if len(set(vars)) < len(vars):
            raise ValueError("range variables should be distinct, but there are duplicates")
    else:
        vars, free_vars = unify_arguments(funcs)

    # check for invalid range (xmin > xmax or ymin > ymax) and swap
    if len(ranges) > 1:
        for i in range(len(ranges)):
            if ranges[i][-2] > ranges[i][-1]:
                ranges[i] = list(ranges[i])
                ranges[i][-1], ranges[i][-2] = ranges[i][-2], ranges[i][-1]
                ranges[i] = tuple(ranges[i])

    # pad the variables if we don't have enough
    nargs = len(ranges)
    if len(vars) < nargs:
        vars += ('_',)*(nargs-len(vars))

    ranges = [[float(z) for z in r] for r in ranges]

    if plot_points is None:
        plot_points = 2

    if not isinstance(plot_points, (list, tuple)):
        plot_points = [plot_points]*len(ranges)
    elif len(plot_points) != nargs:
        raise ValueError("plot_points must be either an integer or a list of integers, one for each range")

    plot_points = [int(p) if p >= 2 else 2 for p in plot_points]
    range_steps = [abs(range[1] - range[0])/(p-1) for range, p in zip(ranges, plot_points)]
    if min(range_steps) == float(0):
        raise ValueError("plot start point and end point must be different")

    eov = False  # eov = "expect one value"
    if nargs == 1:
        eov = True

    from sage.ext.fast_callable import fast_callable

    def try_make_fast(f):
        # If "f" supports fast_callable(), use it. We can't guarantee
        # that our arguments will actually support fast_callable()
        # because, for example, the user may already have done it
        # himself, and the result of fast_callable() can't be
        # fast-callabled again.
        from sage.rings.complex_double import CDF
        from sage.ext.interpreters.wrapper_cdf import Wrapper_cdf

        if hasattr(f, '_fast_callable_'):
            ff = fast_callable(f, vars=vars, expect_one_var=eov, domain=CDF)
            return FastCallablePlotWrapper(ff, imag_tol=imaginary_tolerance)
        elif isinstance(f, Wrapper_cdf):
            # Already a fast-callable, just wrap it. This can happen
            # if, for example, a symbolic expression is passed to a
            # higher-level plot() function that converts it to a
            # fast-callable with expr._plot_fast_callable() before
            # we ever see it.
            return FastCallablePlotWrapper(f, imag_tol=imaginary_tolerance)
        elif callable(f):
            # This will catch python functions, among other things. We don't
            # wrap these yet because we don't know what type they'll return.
            return f
        else:
            # Convert things like ZZ(0) into constant functions.
            from sage.symbolic.ring import SR
            ff = fast_callable(SR(f),
                               vars=vars,
                               expect_one_var=eov,
                               domain=CDF)
            return FastCallablePlotWrapper(ff, imag_tol=imaginary_tolerance)

    # Handle vectors, lists, tuples, etc.
    if isinstance(funcs, Iterable):
        funcs = tuple( try_make_fast(f) for f in funcs )
    else:
        funcs = try_make_fast(funcs)

    #TODO: raise an error if there is a function/method in funcs that takes more values than we have ranges

    if return_vars:
        return (funcs,
                [tuple(_range + [range_step])
                 for _range, range_step in zip(ranges, range_steps)],
                vars)
    else:
        return (funcs,
                [tuple(_range + [range_step])
                 for _range, range_step in zip(ranges, range_steps)])


def unify_arguments(funcs):
    """
    Return a tuple of variables of the functions, as well as the
    number of "free" variables (i.e., variables that defined in a
    callable function).

    INPUT:

    - ``funcs`` -- a list of functions; these can be symbolic
      expressions, polynomials, etc

    OUTPUT: functions, expected arguments

    - A tuple of variables in the functions

    - A tuple of variables that were "free" in the functions

    EXAMPLES::

        sage: x,y,z=var('x,y,z')
        sage: f(x,y)=x+y-z
        sage: g(x,y)=x+y
        sage: h(y)=-y
        sage: sage.plot.misc.unify_arguments((f,g,h))
        ((x, y, z), (z,))
        sage: sage.plot.misc.unify_arguments((g,h))
        ((x, y), ())
        sage: sage.plot.misc.unify_arguments((f,z))
        ((x, y, z), (z,))
        sage: sage.plot.misc.unify_arguments((h,z))
        ((y, z), (z,))
        sage: sage.plot.misc.unify_arguments((x+y,x-y))
        ((x, y), (x, y))
    """
    vars=set()
    free_variables=set()
    if not isinstance(funcs, (list, tuple)):
        funcs = [funcs]

    from sage.structure.element import Expression
    for f in funcs:
        if isinstance(f, Expression) and f.is_callable():
            f_args = set(f.arguments())
            vars.update(f_args)
        else:
            f_args = set()

        try:
            free_vars = set(f.variables()).difference(f_args)
            vars.update(free_vars)
            free_variables.update(free_vars)
        except AttributeError:
            # we probably have a constant
            pass
    return tuple(sorted(vars, key=str)), tuple(sorted(free_variables, key=str))


def _multiple_of_constant(n, pos, const):
    r"""
    Function for internal use in formatting ticks on axes with
    nice-looking multiples of various symbolic constants, such
    as `\pi` or `e`.  Should only be used via keyword argument
    `tick_formatter` in :meth:`plot.show`.  See documentation
    for the matplotlib.ticker module for more details.

    EXAMPLES:

    Here is the intended use::

        sage: plot(sin(x), (x,0,2*pi), ticks=pi/3, tick_formatter=pi)
        Graphics object consisting of 1 graphics primitive

    Here is an unintended use, which yields unexpected (and probably
    undesired) results::

        sage: plot(x^2, (x, -2, 2), tick_formatter=pi)
        Graphics object consisting of 1 graphics primitive

    We can also use more unusual constant choices::

        sage: plot(ln(x), (x,0,10), ticks=e, tick_formatter=e)
        Graphics object consisting of 1 graphics primitive
        sage: plot(x^2, (x,0,10), ticks=[sqrt(2),8], tick_formatter=sqrt(2))
        Graphics object consisting of 1 graphics primitive
    """
    from sage.misc.latex import latex
    from sage.rings.continued_fraction import continued_fraction
    from sage.rings.infinity import Infinity
    cf = continued_fraction(n/const)
    k = 1
    while cf.quotient(k) != Infinity and cf.denominator(k) < 12:
        k += 1
    return '$%s$'%latex(cf.convergent(k-1)*const)


def get_matplotlib_linestyle(linestyle, return_type):
    """
    Function which translates between matplotlib linestyle in short notation
    (i.e. '-', '--', ':', '-.') and long notation (i.e. 'solid', 'dashed',
    'dotted', 'dashdot' ).

    If linestyle is none of these allowed options, the function raises
    a ValueError.

    INPUT:

    - ``linestyle`` - The style of the line, which is one of
       - ``"-"`` or ``"solid"``
       - ``"--"`` or ``"dashed"``
       - ``"-."`` or ``"dash dot"``
       - ``":"`` or ``"dotted"``
       - ``"None"`` or ``" "`` or ``""`` (nothing)

       The linestyle can also be prefixed with a drawing style (e.g., ``"steps--"``)

       - ``"default"`` (connect the points with straight lines)
       - ``"steps"`` or ``"steps-pre"`` (step function; horizontal
         line is to the left of point)
       - ``"steps-mid"`` (step function; points are in the middle of
         horizontal lines)
       - ``"steps-post"`` (step function; horizontal line is to the
         right of point)

       If ``linestyle`` is ``None`` (of type NoneType), then we return it
       back unmodified.

    - ``return_type`` - The type of linestyle that should be output. This
      argument takes only two values - ``"long"`` or ``"short"``.

    EXAMPLES:

    Here is an example how to call this function::

        sage: from sage.plot.misc import get_matplotlib_linestyle
        sage: get_matplotlib_linestyle(':', return_type='short')
        ':'

        sage: get_matplotlib_linestyle(':', return_type='long')
        'dotted'

    TESTS:

    Make sure that if the input is already in the desired format, then it
    is unchanged::

        sage: get_matplotlib_linestyle(':', 'short')
        ':'

    Empty linestyles should be handled properly::

        sage: get_matplotlib_linestyle("", 'short')
        ''
        sage: get_matplotlib_linestyle("", 'long')
        'None'
        sage: get_matplotlib_linestyle(None, 'short') is None
        True

    Linestyles with ``"default"`` or ``"steps"`` in them should also be
    properly handled.  For instance, matplotlib understands only the short
    version when ``"steps"`` is used::

        sage: get_matplotlib_linestyle("default", "short")
        ''
        sage: get_matplotlib_linestyle("steps--", "short")
        'steps--'
        sage: get_matplotlib_linestyle("steps-predashed", "long")
        'steps-pre--'

    Finally, raise error on invalid linestyles::

        sage: get_matplotlib_linestyle("isthissage", "long")
        Traceback (most recent call last):
        ...
        ValueError: WARNING: Unrecognized linestyle 'isthissage'. Possible
        linestyle options are:
        {'solid', 'dashed', 'dotted', dashdot', 'None'}, respectively {'-',
        '--', ':', '-.', ''}

    """
    long_to_short_dict={'solid' : '-','dashed' : '--', 'dotted' : ':',
                        'dashdot':'-.'}
    short_to_long_dict={'-' : 'solid','--' : 'dashed', ':' : 'dotted',
                        '-.':'dashdot'}

    # We need this to take care of region plot. Essentially, if None is
    # passed, then we just return back the same thing.
    if linestyle is None:
        return None

    if linestyle.startswith("default"):
        return get_matplotlib_linestyle(linestyle.strip("default"), "short")
    elif linestyle.startswith("steps"):
        if linestyle.startswith("steps-mid"):
            return "steps-mid" + get_matplotlib_linestyle(
                                    linestyle.strip("steps-mid"), "short")
        elif linestyle.startswith("steps-post"):
            return "steps-post" + get_matplotlib_linestyle(
                                    linestyle.strip("steps-post"), "short")
        elif linestyle.startswith("steps-pre"):
            return "steps-pre" + get_matplotlib_linestyle(
                                    linestyle.strip("steps-pre"), "short")
        else:
            return "steps" + get_matplotlib_linestyle(
                                    linestyle.strip("steps"), "short")

    if return_type == 'short':
        if linestyle in short_to_long_dict.keys():
            return linestyle
        elif linestyle == "" or linestyle == " " or linestyle == "None":
            return ''
        elif linestyle in long_to_short_dict.keys():
            return long_to_short_dict[linestyle]
        else:
            raise ValueError("WARNING: Unrecognized linestyle '%s'. "
                             "Possible linestyle options are:\n{'solid', "
                             "'dashed', 'dotted', dashdot', 'None'}, "
                             "respectively {'-', '--', ':', '-.', ''}"%
                             (linestyle))

    elif return_type == 'long':
        if linestyle in long_to_short_dict.keys():
            return linestyle
        elif linestyle == "" or linestyle == " " or linestyle == "None":
            return "None"
        elif linestyle in short_to_long_dict.keys():
            return short_to_long_dict[linestyle]
        else:
            raise ValueError("WARNING: Unrecognized linestyle '%s'. "
                             "Possible linestyle options are:\n{'solid', "
                             "'dashed', 'dotted', dashdot', 'None'}, "
                             "respectively {'-', '--', ':', '-.', ''}"%
                             (linestyle))


class FastCallablePlotWrapper(FastCallableFloatWrapper):
    r"""
    A fast-callable wrapper for plotting that returns ``nan`` instead
    of raising an error whenever the imaginary tolerance is exceeded.

    A detailed rationale for this can be found in the superclass
    documentation.

    EXAMPLES:

    The ``float`` incarnation of "not a number" is returned instead
    of an error being thrown if the answer is complex::

        sage: from sage.plot.misc import FastCallablePlotWrapper
        sage: f = sqrt(x)
        sage: ff = fast_callable(f, vars=[x], domain=CDF)
        sage: fff = FastCallablePlotWrapper(ff, imag_tol=1e-8)
        sage: fff(1)
        1.0
        sage: fff(-1)
        nan
    """
    def __call__(self, *args):
        r"""
        Evaluate the underlying fast-callable and convert the result to
        ``float``.

        TESTS:

        Evaluation never fails and always returns a ``float``::

            sage: from sage.plot.misc import FastCallablePlotWrapper
            sage: f = x
            sage: ff = fast_callable(f, vars=[x], domain=CDF)
            sage: fff = FastCallablePlotWrapper(ff, imag_tol=0.1)
            sage: type(fff(CDF.random_element())) is float
            True

        """
        try:
            return super().__call__(*args)
        except ValueError:
            return float("nan")
