# sage_setup: distribution = sagemath-objects
"""
Metaclasses for Cython extension types

Cython does not support metaclasses, but this module can be used to
implement metaclasses for extension types.

.. WARNING::

    This module has many caveats and you can easily get segfaults if you
    make a mistake. It relies on undocumented Python and Cython
    behaviour, so things might break in future versions.

How to use
==========

To enable this metaclass mechanism, you need to put
``cimport sage.cpython.cython_metaclass`` in your module (in the ``.pxd``
file if you are using one).

In the extension type (a.k.a. ``cdef class``) for which you want to
define a metaclass, define a method ``__getmetaclass__`` with a single
unused argument, and turn off the Cython directive
``always_allow_keywords``. This method should return a type to be used as
metaclass:

.. code-block:: cython

    cimport cython
    cimport sage.cpython.cython_metaclass
    cdef class MyCustomType():
        @cython.always_allow_keywords(False)
        def __getmetaclass__(_):
            from foo import MyMetaclass
            return MyMetaclass

.. WARNING::

    ``__getmetaclass__`` must be defined as an ordinary method taking a
    single argument, but this argument should not be used in the
    method (it will be ``None``).

When a type ``cls`` is being constructed with metaclass ``meta``,
then ``meta.__init__(cls, None, None, None)`` is called from Cython.
In Python, this would be ``meta.__init__(cls, name, bases, dict)``.

.. WARNING::

    The ``__getmetaclass__`` method is called while the type is being
    created during the import of the module. Therefore,
    ``__getmetaclass__`` should not refer to any global objects,
    including the type being created or other types defined or imported
    in the module (unless you are very careful). Note that non-imported
    ``cdef`` functions are not Python objects, so those are safe to call.

    The same warning applies to the ``__init__`` method of the
    metaclass.

.. WARNING::

    The ``__new__`` method of the metaclass (including the ``__cinit__``
    method for Cython extension types) is never called if you're using
    this from Cython. In particular, the metaclass cannot have any
    attributes or virtual methods.

EXAMPLES::

    sage: cython(                                                                       # needs sage.misc.cython
    ....: '''
    ....: cimport cython
    ....: cimport sage.cpython.cython_metaclass
    ....: cdef class MyCustomType():
    ....:     @cython.always_allow_keywords(False)
    ....:     def __getmetaclass__(_):
    ....:         class MyMetaclass(type):
    ....:             def __init__(*args):
    ....:                 print("Calling MyMetaclass.__init__{}".format(args))
    ....:         return MyMetaclass
    ....:
    ....: cdef class MyDerivedType(MyCustomType):
    ....:     pass
    ....: ''')
    Calling MyMetaclass.__init__(<class '...MyCustomType'>, None, None, None)
    Calling MyMetaclass.__init__(<class '...MyDerivedType'>, None, None, None)
    sage: MyCustomType.__class__                                                        # needs sage.misc.cython
    <class '...MyMetaclass'>
    sage: class MyPythonType(MyDerivedType):                                            # needs sage.misc.cython
    ....:     pass
    Calling MyMetaclass.__init__(<class '...MyPythonType'>, 'MyPythonType', (<class '...MyDerivedType'>,), {...})

Implementation
==============

All this is implemented by defining

.. code-block:: c

    #define PyTypeReady(t)  Sage_PyType_Ready(t)

and then implementing the function ``Sage_PyType_Ready(t)`` which first
calls ``PyType_Ready(t)`` and then handles the metaclass stuff.

TESTS:

Check that a proper exception is raised if ``__getmetaclass__``
returns a non-type::

    sage: cython(                                                                       # needs sage.misc.cython
    ....: '''
    ....: cimport cython
    ....: cimport sage.cpython.cython_metaclass
    ....: cdef class MyCustomType():
    ....:     @cython.always_allow_keywords(False)
    ....:     def __getmetaclass__(_):
    ....:         return 2
    ....: ''')
    Traceback (most recent call last):
    ...
    TypeError: __getmetaclass__ did not return a type
"""

# ****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
