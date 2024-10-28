r"""
Symbolic minimum and maximum

Sage provides a symbolic maximum and minimum due to the fact that the
Python builtin :func:`max` and :func:`min` are not able to deal with variables
as users might expect. These functions wait to evaluate if there are variables.

Here you can see some differences::

   sage: max(x, x^2)                                                                    # needs sage.symbolic
   x
   sage: max_symbolic(x, x^2)                                                           # needs sage.symbolic
   max(x, x^2)
   sage: f(x) = max_symbolic(x, x^2); f(1/2)                                            # needs sage.symbolic
   1/2

This works as expected for more than two entries::

   sage: # needs sage.symbolic
   sage: max(3, 5, x)
   5
   sage: min(3, 5, x)
   3
   sage: max_symbolic(3, 5, x)
   max(x, 5)
   sage: min_symbolic(3, 5, x)
   min(x, 3)
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2010 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
###############################################################################

from builtins import max as builtin_max, min as builtin_min

from sage.misc.lazy_import import lazy_import
from sage.structure.element import Expression
from sage.symbolic.function import BuiltinFunction

lazy_import('sage.symbolic.ring', 'SR')


class MinMax_base(BuiltinFunction):
    def eval_helper(self, this_f, builtin_f, initial_val, args):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: max_symbolic(3, 5, x)  # indirect doctest
            max(x, 5)
            sage: max_symbolic([5.0r])   # indirect doctest
            5.0
            sage: min_symbolic(3, 5, x)
            min(x, 3)
            sage: min_symbolic([5.0r])   # indirect doctest
            5.0
        """
        # __call__ ensures that if args is a singleton, the element is iterable
        arg_is_iter = False
        if len(args) == 1:
            arg_is_iter = True
            args = args[0]

        symb_args = []
        res = initial_val
        num_non_symbolic_args = 0
        for x in args:
            if isinstance(x, Expression):
                symb_args.append(x)
            else:
                num_non_symbolic_args += 1
                if res is None:
                    res = x
                else:
                    res = builtin_f(res, x)

        # if no symbolic arguments, return the result
        if len(symb_args) == 0:
            return res

        # if all arguments were symbolic return
        if num_non_symbolic_args <= 1 and not arg_is_iter:
            return None

        if res is not None:
            symb_args.append(res)
        return this_f(*symb_args)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: max_symbolic(3, 5, x)                                                 # needs sage.symbolic
            max(x, 5)
            sage: max_symbolic(3, 5, x, hold=True)                                      # needs sage.symbolic
            max(3, 5, x)
            sage: max_symbolic([3, 5, x])                                               # needs sage.symbolic
            max(x, 5)

        ::

            sage: min_symbolic(3, 5, x)                                                 # needs sage.symbolic
            min(x, 3)
            sage: min_symbolic(3, 5, x, hold=True)                                      # needs sage.symbolic
            min(3, 5, x)
            sage: min_symbolic([3, 5, x])                                               # needs sage.symbolic
            min(x, 3)

        TESTS:

        We get an exception if no arguments are given::

            sage: max_symbolic()
            Traceback (most recent call last):
            ...
            ValueError: number of arguments must be > 0

        Check if a single argument which is not iterable works::

            sage: # needs sage.symbolic
            sage: max_symbolic(None)
            Traceback (most recent call last):
            ...
            TypeError: 'NoneType' object is not iterable
            sage: max_symbolic(5)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
            sage: max_symbolic(x)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.symbolic.expression.Expression' object is not iterable
            sage: min_symbolic(5)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
            sage: min_symbolic(x)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.symbolic.expression.Expression' object is not iterable
        """
        if len(args) == 0:
            raise ValueError("number of arguments must be > 0")
        if len(args) == 1:
            try:
                args = (SR._force_pyobject(iter(args[0])),)
            except TypeError:
                raise

        try:
            return BuiltinFunction.__call__(self, *args, **kwds)
        except ValueError:
            pass


class MaxSymbolic(MinMax_base):
    def __init__(self):
        r"""
        Symbolic `\max` function.

        The Python builtin :func:`max` function does not work as expected when symbolic
        expressions are given as arguments. This function delays evaluation
        until all symbolic arguments are substituted with values.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: max_symbolic(3, x)
            max(3, x)
            sage: max_symbolic(3, x).subs(x=5)
            5
            sage: max_symbolic(3, 5, x)
            max(x, 5)
            sage: max_symbolic([3, 5, x])
            max(x, 5)

        TESTS::

            sage: loads(dumps(max_symbolic(x, 5)))                                      # needs sage.symbolic
            max(x, 5)
            sage: latex(max_symbolic(x, 5))                                             # needs sage.symbolic
            \max\left(x, 5\right)
            sage: max_symbolic(x, 5)._sympy_()                                          # needs sympy sage.symbolic
            Max(5, x)
        """
        BuiltinFunction.__init__(self, 'max', nargs=0, latex_name=r"\max",
                                 conversions=dict(sympy='Max'))

    def _eval_(self, *args):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: t = max_symbolic(x, 5); t
            max(x, 5)
            sage: t.subs(x=3)  # indirect doctest
            5
            sage: max_symbolic(5, 3)
            5
            sage: u = max_symbolic(*(list(range(10)) + [x])); u
            max(x, 9)
            sage: u.subs(x=-1)
            9
            sage: u.subs(x=10)
            10
            sage: max_symbolic([0, x])
            max(x, 0)

        TESTS::

            sage: max_symbolic()
            Traceback (most recent call last):
            ...
            ValueError: number of arguments must be > 0
        """
        return self.eval_helper(max_symbolic, builtin_max, None, args)

    def _evalf_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: t = max_symbolic(sin(x), cos(x))
            sage: t.subs(x=1).n(200)
            0.84147098480789650665250232163029899962256306079837106567275
            sage: var('y')
            y
            sage: t = max_symbolic(sin(x), cos(x), y)
            sage: u = t.subs(x=1); u
            max(sin(1), cos(1), y)
            sage: u.n()
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression numerically

        We can usually integrate these expressions, but can't
        guarantee a symbolic answer in closed form::

            sage: f = max_symbolic(sin(x), cos(x))                                      # needs sage.symbolic
            sage: r = integral(f, x, 0, 1)                                              # needs sage.symbolic
            ...
            sage: r.n()  # abs tol 1e-8                                                 # needs sage.symbolic
            0.873911256504955
        """
        return max_symbolic(args)


max_symbolic = MaxSymbolic()


class MinSymbolic(MinMax_base):
    def __init__(self):
        r"""
        Symbolic `\min` function.

        The Python builtin :func:`min` function does not work as expected when symbolic
        expressions are given as arguments. This function delays evaluation
        until all symbolic arguments are substituted with values.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: min_symbolic(3, x)
            min(3, x)
            sage: min_symbolic(3, x).subs(x=5)
            3
            sage: min_symbolic(3, 5, x)
            min(x, 3)
            sage: min_symbolic([3, 5, x])
            min(x, 3)

        TESTS::

            sage: loads(dumps(min_symbolic(x, 5)))                                      # needs sage.symbolic
            min(x, 5)
            sage: latex(min_symbolic(x, 5))                                             # needs sage.symbolic
            \min\left(x, 5\right)
            sage: min_symbolic(x, 5)._sympy_()                                          # needs sympy sage.symbolic
            Min(5, x)
        """
        BuiltinFunction.__init__(self, 'min', nargs=0, latex_name=r"\min",
                                 conversions=dict(sympy='Min'))

    def _eval_(self, *args):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: t = min_symbolic(x, 5); t
            min(x, 5)
            sage: t.subs(x=3)  # indirect doctest
            3
            sage: min_symbolic(5, 3)
            3
            sage: u = min_symbolic(*(list(range(10)) + [x])); u
            min(x, 0)
            sage: u.subs(x=-1)
            -1
            sage: u.subs(x=10)
            0
            sage: min_symbolic([3, x])
            min(x, 3)

        TESTS::

            sage: min_symbolic()
            Traceback (most recent call last):
            ...
            ValueError: number of arguments must be > 0
        """
        return self.eval_helper(min_symbolic, builtin_min, float('inf'), args)

    def _evalf_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: t = min_symbolic(sin(x), cos(x))
            sage: t.subs(x=1).n(200)
            0.54030230586813971740093660744297660373231042061792222767010
            sage: var('y')
            y
            sage: t = min_symbolic(sin(x), cos(x), y)
            sage: u = t.subs(x=1); u
            min(sin(1), cos(1), y)
            sage: u.n()
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression numerically
        """
        return min_symbolic(args)


min_symbolic = MinSymbolic()
