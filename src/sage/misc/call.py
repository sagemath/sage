# -*- coding: utf-8 -*-
"""
Attribute and method calling
"""

# ****************************************************************************
#       Copyright (C) 2008 Mike Hansen
#       Copyright (C) 2010, 2013 Nicolas M. Thiery
#       Copyright (C) 2018 Frédéric Chapoton
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

#############################################
# Operators
#############################################
class AttrCallObject(object):
    def __init__(self, name, args, kwds):
        """
        TESTS::

            sage: f = attrcall('core', 3); f
            *.core(3)
            sage: TestSuite(f).run()
        """
        self.name = name
        self.args = args
        self.kwds = kwds

    def __call__(self, x, *args):
        """
        Gets the ``self.name`` method from ``x``, calls it with
        ``self.args`` and ``args`` as positional parameters and
        ``self.kwds`` as keyword parameters, and returns the result.

        EXAMPLES::

            sage: core = attrcall('core', 3)
            sage: core(Partition([4,2]))
            [4, 2]

            sage: series = attrcall('series', x)
            sage: series(sin(x), 4)
            1*x + (-1/6)*x^3 + Order(x^4)
        """
        return getattr(x, self.name)(*(self.args + args), **self.kwds)

    def __repr__(self):
        """
        Return a string representation of this object.

        The star in the output represents the object passed into ``self``.

        EXAMPLES::

            sage: attrcall('core', 3)
            *.core(3)
            sage: attrcall('hooks', flatten=True)
            *.hooks(flatten=True)
            sage: attrcall('hooks', 3, flatten=True)
            *.hooks(3, flatten=True)
        """
        s = "*.%s(%s" % (self.name, ", ".join(map(repr, self.args)))
        if self.kwds:
            if self.args:
                s += ", "
            s += ", ".join("%s=%s" % keyvalue for keyvalue in self.kwds.items())
        s += ")"
        return s

    def __eq__(self, other):
        """
        Equality testing

        EXAMPLES::

            sage: attrcall('core', 3, flatten = True) == attrcall('core', 3, flatten = True)
            True
            sage: attrcall('core', 2) == attrcall('core', 3)
            False
            sage: attrcall('core', 2) == 1
            False
        """
        return self.__class__ == other.__class__ and self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Equality testing

        EXAMPLES::

            sage: attrcall('core', 3, flatten = True) != attrcall('core', 3, flatten = True)
            False
            sage: attrcall('core', 2) != attrcall('core', 3)
            True
            sage: attrcall('core', 2) != 1
            True
        """
        return not self == other

    def __hash__(self):
        """
        Hash value

        This method tries to ensure that, when two ``attrcall``
        objects are equal, they have the same hash value.

        .. warning:: dicts are not hashable, so we instead hash their
        items; however the order of those items might differ. The
        proper fix would be to use a frozen dict for ``kwds``, when
        frozen dicts will be available in Python.

        EXAMPLES::

            sage: x = attrcall('core', 3, flatten = True, blah = 1)
            sage: hash(x)       # random # indirect doctest
            210434060
            sage: type(hash(x))
            <class 'int'>
            sage: y = attrcall('core', 3, blah = 1, flatten = True)
            sage: hash(y) == hash(x)
            True
            sage: y = attrcall('core', 3, flatten = True, blah = 2)
            sage: hash(y) != hash(x)
            True
            sage: hash(attrcall('core', 2)) != hash(attrcall('core', 3))
            True
            sage: hash(attrcall('core', 2)) != hash(1)
            True

        Note: a missing ``__hash__`` method here used to break the
        unique representation of parents taking ``attrcall`` objects
        as input; see :trac:`8911`.
        """
        return hash((self.args, tuple(sorted(self.kwds.items()))))


def attrcall(name, *args, **kwds):
    """
    Return a callable which takes in an object, gets the method named
    name from that object, and calls it with the specified arguments
    and keywords.

    INPUT:

    -  ``name`` - a string of the name of the method you
       want to call

    -  ``args, kwds`` - arguments and keywords to be passed
       to the method

    EXAMPLES::

        sage: f = attrcall('core', 3); f
        *.core(3)
        sage: [f(p) for p in Partitions(5)]
        [[2], [1, 1], [1, 1], [3, 1, 1], [2], [2], [1, 1]]
    """
    return AttrCallObject(name, args, kwds)


def call_method(obj, name, *args, **kwds):
    """
    Call the method ``name`` on ``obj``.

    This has to exist somewhere in Python!!!

    .. SEEALSO:: :func:`operator.methodcaller` :func:`attrcal`

    EXAMPLES::

        sage: from sage.misc.call import call_method
        sage: call_method(1, "__add__", 2)
        3
    """
    return getattr(obj, name)(*args, **kwds)

from sage.misc.persist import register_unpickle_override
register_unpickle_override("sage.misc.misc", "call_method", call_method)
