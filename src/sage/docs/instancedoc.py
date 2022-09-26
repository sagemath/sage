r"""
Dynamic documentation for instances of classes

The functionality in this module allows to define specific docstrings
of *instances* of a class, which are different from the class docstring.
A typical use case is given by cached methods: the documentation of a
cached method should not be the documentation of the class
:class:`CachedMethod`; it should be the documentation of the underlying
method.

In order to use this, define a class docstring as usual. Also define a
method ``def _instancedoc_(self)`` which should return the docstring of
the instance ``self``. Finally, add the decorator ``@instancedoc`` to
the class.

.. WARNING::

    Since the ``__doc__`` attribute is never inherited, the decorator
    ``@instancedoc`` must be added to all subclasses of the class
    defining ``_instancedoc_``. Doing it on the base class is not
    sufficient.

EXAMPLES::

    sage: from sage.misc.instancedoc import instancedoc
    sage: @instancedoc
    ....: class X():
    ....:     "Class docstring"
    ....:     def _instancedoc_(self):
    ....:         return "Instance docstring"
    sage: X.__doc__
    'Class docstring'
    sage: X().__doc__
    'Instance docstring'

For a Cython ``cdef class``, a decorator cannot be used. Instead, call
:func:`instancedoc` as a function after defining the class::

    sage: cython('''  # optional - sage.misc.cython
    ....: from sage.misc.instancedoc import instancedoc
    ....: cdef class Y:
    ....:     "Class docstring"
    ....:     def _instancedoc_(self):
    ....:         return "Instance docstring"
    ....: instancedoc(Y)
    ....: ''')
    sage: Y.__doc__
    'File:...\nClass docstring'
    sage: Y().__doc__
    'Instance docstring'

One can still add a custom ``__doc__`` attribute on a particular
instance::

    sage: obj = X()
    sage: obj.__doc__ = "Very special doc"
    sage: print(obj.__doc__)
    Very special doc

This normally does not work on extension types::

    sage: Y().__doc__ = "Very special doc"
    Traceback (most recent call last):
    ...
    AttributeError: attribute '__doc__' of 'Y' objects is not writable

This is an example involving a metaclass, where the instances are
classes. In this case, the ``_instancedoc_`` from the metaclass is only
used if the instance of the metaclass (the class) does not have a
docstring::

    sage: @instancedoc
    ....: class Meta(type):
    ....:     "Metaclass doc"
    ....:     def _instancedoc_(self):
    ....:         return "Docstring for {}".format(self)
    sage: class T(metaclass=Meta):
    ....:     pass
    sage: print(T.__doc__)
    Docstring for <class '__main__.T'>
    sage: class U(metaclass=Meta):
    ....:     "Special doc for U"
    sage: print(U.__doc__)
    Special doc for U

TESTS:

Check that inheritance works (after passing the subclass to
:func:`instancedoc`)::

    sage: @instancedoc
    ....: class A():
    ....:     "Class A docstring"
    ....:     def _instancedoc_(self):
    ....:         return "Instance docstring"
    sage: class B(A):
    ....:     "Class B docstring"
    sage: B.__doc__
    'Class B docstring'
    sage: B().__doc__  # Ideally, this would return the instance docstring
    'Class B docstring'
    sage: B = instancedoc(B)
    sage: B.__doc__
    'Class B docstring'
    sage: B().__doc__
    'Instance docstring'
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(33763, 'This module is deprecated. Use "sage.misc.instancedoc" instead.')

from sage.misc.instancedoc import instancedoc
