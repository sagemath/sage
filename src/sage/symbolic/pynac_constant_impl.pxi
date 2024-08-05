"""
Wrapper around Pynac's constants
"""

# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2009 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.cpython.string cimport str_to_bytes


cdef class PynacConstant:

    cdef GConstant* pointer
    cdef GConstant* _object
    cdef _name

    def __cinit__(self, name, texname, domain_string):
        """
        Create a constant in Pynac.

        EXAMPLES::

            sage: from sage.symbolic.expression import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f
            foo

        Note that this just creates a 'constant' object and not an
        expression.  If you want to work with this constant, you'll
        have to use the result of :meth:`expression`::

            sage: foo = f.expression()
            sage: foo + 2
            foo + 2
        """
        cdef unsigned domain
        if domain_string == 'complex':
            domain = domain_complex
        elif domain_string == 'real':
            domain = domain_real
        elif domain_string == 'positive':
            domain = domain_positive
        else:
            raise ValueError

        self._name = name

        # For the constants explicitly defined in constant.cpp in the
        # Pynac library, we use those symbols. Otherwise, we create a
        # new constant stored in *self._object
        if self._name == "pi":
            self.pointer = <GConstant *>&g_Pi
        elif self._name == "catalan":
            self.pointer = <GConstant *>&g_Catalan
        elif self._name == "euler_gamma":
            self.pointer = <GConstant *>&g_Euler
        elif self._name == "NaN":
            self.pointer = <GConstant *>&g_NaN
        else:
            self._object = new GConstant(str_to_bytes(name), ConstantEvalf,
                                         str_to_bytes(texname), domain)
            self.pointer = self._object

    def __dealloc__(self):
        del self._object

    def serial(self):
        """
        Return the underlying Pynac serial for this constant.

        EXAMPLES::

            sage: from sage.symbolic.expression import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f.serial()  #random
            15
        """
        return int(self.pointer.get_serial())

    def name(self):
        """
        Return the name of this constant.

        EXAMPLES::

            sage: from sage.symbolic.expression import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f.name()
            'foo'
        """
        return self._name

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real'); f
            foo
        """
        return self.name()

    def expression(self):
        """
        Return this constant as an Expression.

        EXAMPLES::

            sage: from sage.symbolic.expression import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f + 2
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +: '<class 'sage.symbolic.expression.PynacConstant'>' and 'Integer Ring'

            sage: foo = f.expression(); foo
            foo
            sage: foo + 2
            foo + 2
        """
        from sage.symbolic.ring import SR
        return new_Expression_from_GEx(SR, <GEx>(self.pointer[0]))
