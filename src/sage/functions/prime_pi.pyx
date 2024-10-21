# sage_setup: distribution = sagemath-categories
# sage.doctest: needs primecountpy
r"""
Counting primes

EXAMPLES::

    sage: z = sage.functions.prime_pi.PrimePi()
    sage: loads(dumps(z))
    prime_pi
    sage: loads(dumps(z)) == z
    True

AUTHORS:

- \R. Andrew Ohana (2009): initial version of efficient prime_pi

- William Stein (2009): fix plot method

- \R. Andrew Ohana (2011): complete rewrite, ~5x speedup

- Dima Pasechnik (2021): removed buggy cython code, replaced it with
  calls to primecount/primecountpy spkg
"""

# ****************************************************************************
#       Copyright (C) 2009,2011 R. Andrew Ohana <andrew.ohana@gmail.com>
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.rings.integer cimport Integer
from sage.symbolic.function cimport BuiltinFunction

_prime_pi = LazyImport('primecountpy.primecount', 'prime_pi', as_name='prime_pi')
_phi = LazyImport('primecountpy.primecount', 'phi', as_name='_phi')


cdef class PrimePi(BuiltinFunction):
    def __init__(self):
        r"""
        The prime counting function, which counts the number of primes less
        than or equal to a given value.

        INPUT:

        - ``x`` -- a real number
        - ``prime_bound`` -- (default: 0) a real number `< 2^32`; :func:`prime_pi`
          will make sure to use all the primes up to ``prime_bound`` (although,
          possibly more) in computing ``prime_pi``, this can potentially
          speedup the time of computation, at a cost to memory usage.

        OUTPUT: integer; the number of primes :math:`\leq` ``x``

        EXAMPLES:

        These examples test common inputs::

            sage: # needs sage.symbolic
            sage: prime_pi(7)
            4
            sage: prime_pi(100)
            25
            sage: prime_pi(1000)
            168
            sage: prime_pi(100000)
            9592
            sage: prime_pi(500509)
            41581

        The following test is to verify that :issue:`4670` has been essentially
        resolved::

            sage: prime_pi(10^10)                                                       # needs sage.symbolic
            455052511

        The :func:`prime_pi` function also has a special plotting method, so it
        plots quickly and perfectly as a step function::

            sage: P = plot(prime_pi, 50, 100)                                           # needs sage.plot sage.symbolic
        """
        super(PrimePi, self).__init__('prime_pi', latex_name=r"\pi",
                                      conversions={'mathematica': 'PrimePi',
                                                   'pari': 'primepi',
                                                   'sympy': 'primepi'})

    def __call__(self, *args, coerce=True, hold=False):
        r"""
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: prime_pi.__call__(756)
            133
            sage: prime_pi.__call__(6574, 577)
            850
            sage: f(x) = prime_pi.__call__(x^2); f(x)
            prime_pi(x^2)
            sage: f(5)
            9
            sage: prime_pi.__call__(1, 2, 3)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function prime_pi takes 1 or 2 arguments (3 given)
        """
        if len(args) > 2:
            raise TypeError(f"Symbolic function {self._name} takes 1 or 2"
                            f" arguments ({len(args)} given)")
        return super(PrimePi, self).__call__(args[0], coerce=coerce, hold=hold)

    def _eval_(self, x):
        r"""
        EXAMPLES::

            sage: prime_pi._eval_(7)
            4
            sage: prime_pi._eval_(100)
            25
            sage: prime_pi._eval_(1000)
            168
            sage: prime_pi._eval_(100000)
            9592
            sage: prime_pi._eval_(500509)
            41581
            sage: prime_pi._eval_(mod(30957, 9750979))
            3337

        Make sure we actually compute correct results for 64-bit entries::

            sage: for i in (32..42): prime_pi(2^i)      # long time (13s on sage.math, 2011)
            203280221
            393615806
            762939111
            1480206279
            2874398515
            5586502348
            10866266172
            21151907950
            41203088796
            80316571436
            156661034233

        This implementation uses 64-bit ints and does not support
        :math:`x \geq 2^63`::

            sage: prime_pi(2^63)
            Traceback (most recent call last):
            ...
            OverflowError: ...to convert...

        TESTS:

        Check that :issue:`24960` is fixed::

            sage: prime_pi(642763101936913)
            19439675999019
            sage: prime_pi(10.5)
            4
        """
        from sage.functions.other import floor
        try:
            z = Integer(x)
        except TypeError:
            try:
                z = Integer(floor(x))
            except TypeError:
                return None
        return _prime_pi(z)

    def plot(self, xmin=0, xmax=100, vertical_lines=True, **kwds):
        """
        Draw a plot of the prime counting function from ``xmin`` to ``xmax``.
        All additional arguments are passed on to the line command.

        WARNING: we draw the plot of ``prime_pi`` as a stairstep function with
        explicitly drawn vertical lines where the function jumps. Technically
        there should not be any vertical lines, but they make the graph look
        much better, so we include them. Use the option ``vertical_lines=False``
        to turn these off.

        EXAMPLES::

            sage: plot(prime_pi, 1, 100)                                                # needs sage.plot sage.symbolic
            Graphics object consisting of 1 graphics primitive
            sage: prime_pi.plot(1, 51, thickness=2, vertical_lines=False)               # needs sage.plot sage.symbolic
            Graphics object consisting of 16 graphics primitives
        """
        from sage.plot.step import plot_step_function
        if xmax < xmin:
            return plot_step_function([], **kwds)
        if xmax < 2:
            return plot_step_function([(xmin, 0), (xmax, 0)], **kwds)
        y = self(xmin)
        v = [(xmin, y)]
        from sage.rings.fast_arith import prime_range
        for p in prime_range(xmin+1, xmax+1, py_ints=True):
            y += 1
            v.append((p, y))
        v.append((xmax, y))
        return plot_step_function(v, vertical_lines=vertical_lines, **kwds)


prime_pi = PrimePi()


cpdef Integer legendre_phi(x, a):
    r"""
    Legendre's formula, also known as the partial sieve function, is a useful
    combinatorial function for computing the prime counting function (the
    ``prime_pi`` method in Sage). It counts the number of positive integers
    :math:`\leq` ``x`` that are not divisible by the first ``a`` primes.

    INPUT:

    - ``x`` -- a real number

    - ``a`` -- nonnegative integer

    OUTPUT: integer; the number of positive integers :math:`\leq` ``x`` that
    are not divisible by the first ``a`` primes

    EXAMPLES::

        sage: legendre_phi(100, 0)
        100
        sage: legendre_phi(29375, 1)
        14688
        sage: legendre_phi(91753, 5973)
        2893
        sage: legendre_phi(4215701455, 6450023226)
        1
    """
    if not isinstance(a, Integer):
        a = Integer(a)
    if a < Integer(0):
        raise ValueError("a (=%s) must be nonnegative" % a)
    y = Integer(x)

    # legendre_phi(x, a) = 0 when x <= 0
    if not y:
        return Integer(0)

    # legendre_phi(x, 0) = x
    if a == Integer(0):
        return Integer(y)

    # If a > prime_pi(2^32), we compute phi(x,a) = max(pi(x)-a+1,1)
    if a > Integer(203280221):
        ret = prime_pi(x)-a+Integer(1)
        if ret < Integer(1):
            return Integer(1)
        return ret

    # Deal with the general case
    return Integer(_phi(y, a))


partial_sieve_function = legendre_phi
