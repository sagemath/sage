r"""
Support for persistent functions in .sage files

Persistent functions are functions whose values are stored on disk
so they do not have to be recomputed.

The inputs to the function must be hashable (so lists are not
allowed). Though a hash is used, in the incredibly unlikely event
that a hash collision occurs, your function will not return an
incorrect result because of this (though the cache might not be
used either).

This is meant to be used from ``.sage`` files, not from
library ``.py`` files.

To use this disk caching mechanism, just put
``@func_persist`` right before your function
definition. For example,

::

    @func_persist
    def bern(n):
        "Return the n-th Bernoulli number, caching the result to disk."
        return bernoulli(n)

You can then use the function ``bern`` as usual, except
it will almost instantly return values that have already been
computed, even if you quit and restart.

The disk cache files are stored by default in the subdirectory
``func_persist`` of the current working directory,
with one file for each evaluation of the function.
"""
########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

import inspect
import os

from sage.misc import persist


class func_persist:
    r"""
    Put ``@func_persist`` right before your function
    definition to cache values it computes to disk.
    """
    def __init__(self, f, dir='func_persist'):
        self.__func = f
        self.__dir = dir
        os.makedirs(dir, exist_ok=True)
        self.__doc__ = '%s%s%s' % (
            f.__name__,
            inspect.signature(f),
            f.__doc__)

    def __call__(self, *args, **kwds):
        key = (tuple(args), tuple(kwds.items()))
        h = hash(key)
        name = '%s/%s_%s.sobj' % (self.__dir, self.__func.__name__, h)

        if os.path.exists(name):
            key2, val = persist.load(name)
            if key == key2:
                # We save and test equality of keys to avoid
                # the (extremely remote) possibility of a hash
                # collision.  Correctness is crucial in mathematics.
                return val

        val = self.__func(*args, **kwds)
        persist.save((key, val), name)
        return val
