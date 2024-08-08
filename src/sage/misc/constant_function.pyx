# sage_setup: distribution = sagemath-objects
r"""
Constant functions
"""

#*****************************************************************************
#       Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.richcmp cimport richcmp
from sage.structure.sage_object cimport SageObject


cdef class ConstantFunction(SageObject):
    """
    A class for function objects implementing constant functions.

    EXAMPLES::

        sage: f = ConstantFunction(3)
        sage: f
        The constant function (...) -> 3
        sage: f()
        3
        sage: f(5)
        3

    Such a function could be implemented as a lambda expression, but
    this is not (currently) picklable::

        sage: g = lambda x: 3
        sage: g == loads(dumps(g))
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle ...: attribute lookup ... failed
        sage: f == loads(dumps(f))
        True

    Also, in the long run, the information that this function is
    constant could be used by some algorithms.

    .. TODO::

        - Should constant functions have unique representation?
        - Should the number of arguments be specified in the input?
        - Should this go into ``sage.categories.maps``?
          Then what should be the parent (e.g. for ``lambda x: True``)?

    TESTS:

    These tests do fail if we try to use ``UniqueRepresentation``::

        sage: f = ConstantFunction(True)
        sage: g = ConstantFunction(1)
        sage: f(), g()
        (True, 1)

    That's because ``1`` and ``True`` cannot be distinguished as keys
    in a dictionary (argl!)::

        sage: { 1: 'a', True: 'b' }
        {1: 'b'}
    """
    cdef object _value

    def __init__(self, value):
        """
        EXAMPLES::

            sage: ConstantFunction(1)()
            1
        """
        self._value = value

    def __reduce__(self):
        """
        TESTS::

            sage: loads(dumps(ConstantFunction(5))) == ConstantFunction(5) # indirect doctest
            True
        """
        return ConstantFunction, (self._value,)

    def _repr_(self):
        """
        EXAMPLES::

            sage: ConstantFunction(1)
            The constant function (...) -> 1
        """
        return "The constant function (...) -> %s"%self._value

    def __call__(self, *args):
        """
        EXAMPLES::

            sage: ConstantFunction(1)()
            1
            sage: ConstantFunction(1)(5,3)
            1
            sage: ConstantFunction(True)()
            True
        """
        return self._value

    def __richcmp__(self, other, op):
        """
        EXAMPLES::

            sage: ConstantFunction(1) == ConstantFunction(1)
            True
            sage: ConstantFunction(1) == ConstantFunction(3)
            False
            sage: ConstantFunction(1) == 1
            False
            sage: ConstantFunction(True) == ConstantFunction(1)  # argl!
            True
        """
        if not isinstance(other, ConstantFunction):
            return NotImplemented
        return richcmp((<ConstantFunction>self)._value,
                       (<ConstantFunction>other)._value, op)
