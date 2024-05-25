"""
Parallel Iterator built using Python's multiprocessing module
"""

################################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of (any version of) the GNU
#  General Public License (GPL). The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
################################################################################

from functools import partial
from sage.misc.fpickle import pickle_function, call_pickled_function
from . import ncpus


def pyprocessing(processes=0):
    """
    Return a parallel iterator using a given number of processes
    implemented using pyprocessing.

    INPUT:

    - ``processes`` -- integer (default: 0); if 0, set to the number
      of processors on the computer

    OUTPUT: a (partially evaluated) function

    EXAMPLES::

        sage: from sage.parallel.multiprocessing_sage import pyprocessing
        sage: p_iter = pyprocessing(4)
        sage: P = parallel(p_iter=p_iter)
        sage: def f(x): return x+x
        sage: v = list(P(f)(list(range(10)))); v.sort(); v
        [(((0,), {}), 0), (((1,), {}), 2), (((2,), {}), 4), (((3,), {}), 6), (((4,), {}), 8), (((5,), {}), 10), (((6,), {}), 12), (((7,), {}), 14), (((8,), {}), 16), (((9,), {}), 18)]
    """
    if processes == 0:
        processes = ncpus.ncpus()
    return partial(parallel_iter, processes)


def parallel_iter(processes, f, inputs):
    """
    Return a parallel iterator.

    INPUT:

    - ``processes`` -- integer
    - ``f`` -- function
    - ``inputs`` -- an iterable of pairs (args, kwds)

    OUTPUT: iterator over values of ``f`` at ``args, kwds`` in some random order

    EXAMPLES::

        sage: def f(x): return x+x
        sage: import sage.parallel.multiprocessing_sage
        sage: v = list(sage.parallel.multiprocessing_sage.parallel_iter(2, f, [((2,), {}), ((3,),{})]))
        sage: v.sort(); v
        [(((2,), {}), 4), (((3,), {}), 6)]
    """
    from multiprocessing import Pool

    if processes == 0:
        processes = ncpus.ncpus()
    p = Pool(processes)
    fp = pickle_function(f)

    result = p.imap_unordered(call_pickled_function,
                              [(fp, t) for t in inputs])
    yield from result
    p.close()
    p.join()
