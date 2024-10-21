# sage_setup: distribution = sagemath-categories
"""
Multiplex calls to one object to calls to many objects

AUTHORS:

- Martin Albrecht (2011): initial version
"""
# ****************************************************************************
#       Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
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


class MultiplexFunction:
    """
    A simple wrapper object for functions that are called on a list of
    objects.
    """
    def __init__(self, multiplexer, name):
        """
        EXAMPLES::

            sage: from sage.misc.object_multiplexer import Multiplex, MultiplexFunction
            sage: m = Multiplex(1,1/2)
            sage: f = MultiplexFunction(m,'str')
            sage: f
            <sage.misc.object_multiplexer.MultiplexFunction object at 0x...>
        """
        self.multiplexer = multiplexer
        self.name = name

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: from sage.misc.object_multiplexer import Multiplex, MultiplexFunction
            sage: m = Multiplex(1,1/2)
            sage: f = MultiplexFunction(m,'str')
            sage: f()
            ('1', '1/2')
        """
        l = [getattr(child, self.name)(*args, **kwds)
             for child in self.multiplexer.children]
        if all(e is None for e in l):
            return None
        return tuple(l)


class Multiplex:
    """
    Object for a list of children such that function calls on this
    new object implies that the same function is called on all
    children.
    """
    def __init__(self, *args):
        """
        EXAMPLES::

            sage: from sage.misc.object_multiplexer import Multiplex
            sage: m = Multiplex(1,1/2)
            sage: m.str()
            ('1', '1/2')
        """
        self.children = [arg for arg in args if arg is not None]

    def __getattr__(self, name):
        """
        EXAMPLES::

            sage: from sage.misc.object_multiplexer import Multiplex
            sage: m = Multiplex(1,1/2)
            sage: m.str
            <sage.misc.object_multiplexer.MultiplexFunction object at 0x...>
            sage: m.__bork__
            Traceback (most recent call last):
            ...
            AttributeError: 'Multiplex' has no attribute '__bork__'...
        """
        if name.startswith("__"):
            raise AttributeError("'Multiplex' has no attribute '%s'" % name)
        return MultiplexFunction(self, name)
