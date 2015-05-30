"""
A weakref cache factory
"""

#*****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import weakref

class Cache:
    """
    A weakref cache for arbitrary objects via an arbitrary function.
    """
    def __init__(self, factory_function): #, canonical_params_function=None):
        """
        INPUT:
            factory_function -- a function that returns objects (which this
            class will weakref cache)
        """
        self.cache = {}
        self.factory = factory_function
        #self.canonical_params_function = canonical_params_function

    def key(self, *args, **kwds):
        """
        Return the key associated to the given args.  Note that
        the values must be hashable.
        """
        return (tuple(args), tuple(kwds.items()))

    def __call__(self, *args, **kwds):
        """
        Return the object from the cache defined by the arguments,
        or make it using the factoring function if it doesn't exist.
        """
        #if self.canonical_params_function is not None:
        #    args = self.canonical_params_function(*args, **kwds)
        #    kwds = {}

        key = self.key(*args, **kwds)
        try:
            x = self.cache[key]()
            if not (x is None):
                return x
        except KeyError:
            pass
        x = self.factory(*args, **kwds)
        self.cache[key] = weakref.ref(x)
        return x

    def has_object(self, *args, **kwds):
        """
        Returns true if cache has object defined by the
        given args.
        """
        key = self.key(args, kwds)
        try:
            return not (self.cache[key]() is None)
        except KeyError:
            return False

    def __setitem__(self, key, x):
        """
        Make a weakref to x with given key.
        """
        self.cache[key] = weakref.ref(x)

    def format_names(self, names, n=1):
        """
        Many objects have a names argument.  This function formats it,
        i.e., if it is a list it turns it into a tuple, and if it
        is a string and the number of names needed is bigger than 1,
        it makes each letter of the string an element of a tuple.
        """
        if isinstance(names, list) or (isinstance(names, str) and n > 1):
            return tuple(names)
        return names



