"""
Context Managers for LibGAP

This module implements a context manager for global variables. This is
useful since the behavior of GAP is sometimes controlled by global
variables, which you might want to switch to a different value for a
computation. Here is an example how you are suppose to use it from
your code. First, let us set a dummy global variable for our example::

    sage: libgap.set_global('FooBar', 123)

Then, if you want to switch the value momentarily you can write::

    sage: with libgap.global_context('FooBar', 'test'):
    ....:     print(libgap.get_global('FooBar'))
    test

Afterward, the global variable reverts to the previous value::

    sage: print(libgap.get_global('FooBar'))
    123

The value is reset even if exceptions occur::

    sage: with libgap.global_context('FooBar', 'test'):
    ....:     print(libgap.get_global('FooBar'))
    ....:     raise ValueError(libgap.get_global('FooBar'))
    Traceback (most recent call last):
    ...
    ValueError: test
    sage: print(libgap.get_global('FooBar'))
    123
"""


###############################################################################
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################

from sage.libs.gap.libgap import libgap


class GlobalVariableContext():

    def __init__(self, variable, value):
        """
        Context manager for GAP global variables.

        It is recommended that you use the
        :meth:`sage.libs.gap.libgap.Gap.global_context` method and not
        construct objects of this class manually.

        INPUT:

        - ``variable`` -- string; the variable name

        - ``value`` -- anything that defines a GAP object

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: with libgap.global_context('FooBar', 2):
            ....:     print(libgap.get_global('FooBar'))
            2
            sage: libgap.get_global('FooBar')
            1
        """
        self._variable = variable
        self._new_value = value

    def __enter__(self):
        """
        Called when entering the with-block.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: with libgap.global_context('FooBar', 2):
            ....:     print(libgap.get_global('FooBar'))
            2
            sage: libgap.get_global('FooBar')
            1
        """
        self._old_value = libgap.get_global(self._variable)
        libgap.set_global(self._variable, self._new_value)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Called when exiting the with-block.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: with libgap.global_context('FooBar', 2):
            ....:     print(libgap.get_global('FooBar'))
            2
            sage: libgap.get_global('FooBar')
            1
        """
        libgap.set_global(self._variable, self._old_value)
        return False
