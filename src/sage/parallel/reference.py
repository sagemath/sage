# sage_setup: distribution = sagemath-categories
"""
Reference Parallel Primitives

These are reference implementations of basic parallel
primitives. These are not actually parallel, but work the same way.
They are good for testing.
"""

from sage.misc.prandom import shuffle


def parallel_iter(f, inputs):
    """
    Reference parallel iterator implementation.

    INPUT:

    - ``f`` -- a Python function that can be pickled using
      the ``pickle_function`` command

    - ``inputs`` -- list of pickleable pairs (args, kwds), where
      args is a tuple and kwds is a dictionary

    OUTPUT: iterator over 2-tuples ``(inputs[i], f(inputs[i]))``, where the
    order may be completely random

    EXAMPLES::

        sage: def f(N,M=10): return N*M
        sage: inputs = [((2,3),{}),  (tuple(), {'M':5,'N':3}), ((2,),{})]
        sage: set_random_seed(0)
        sage: for a, val in sage.parallel.reference.parallel_iter(f, inputs):
        ....:     print((a, val))
        (((2,), {}), 20)
        (((), {'M': 5, 'N': 3}), 15)
        (((2, 3), {}), 6)
        sage: for a, val in sage.parallel.reference.parallel_iter(f, inputs):
        ....:     print((a, val))
        (((), {'M': 5, 'N': 3}), 15)
        (((2,), {}), 20)
        (((2, 3), {}), 6)
    """
    v = list(inputs)
    shuffle(v)
    for args, kwds in v:
        yield ((args, kwds), f(*args, **kwds))
