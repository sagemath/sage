# sage_setup: distribution = sagemath-objects
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

    sage: cython(                                                                       # needs sage.misc.cython
    ....: '''
    ....: from sage.misc.instancedoc import instancedoc
    ....: cdef class Y:
    ....:     "Class docstring"
    ....:     def _instancedoc_(self):
    ....:         return "Instance docstring"
    ....: instancedoc(Y)
    ....: ''')
    sage: Y.__doc__                                                                     # needs sage.misc.cython
    'File:...\nClass docstring'
    sage: Y().__doc__                                                                   # needs sage.misc.cython
    'Instance docstring'

One can still add a custom ``__doc__`` attribute on a particular
instance::

    sage: obj = X()
    sage: obj.__doc__ = "Very special doc"
    sage: print(obj.__doc__)
    Very special doc

This normally does not work on extension types::

    sage: Y().__doc__ = "Very special doc"                                              # needs sage.misc.cython
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

from cpython.object cimport PyObject, PyTypeObject

cdef extern from *:
    cdef int PyDict_SetItemString(PyObject*, const char*, object) except -1
    cdef void PyType_Modified(PyTypeObject*)

cdef inline PyTypeObject* TypeObject(cls) except NULL:
    if not isinstance(cls, type):
        raise TypeError(f"expected type, got {cls!r}")
    return <PyTypeObject*>cls


cdef class InstanceDocDescriptor:
    """
    Descriptor for dynamic documentation, to be installed as the
    ``__doc__`` attribute.

    INPUT:

    - ``classdoc`` -- string; class documentation

    - ``instancedoc`` -- (method) documentation for an instance

    - ``attr`` -- string (default: ``__doc__``); attribute name to use
      for custom docstring on the instance

    EXAMPLES::

        sage: from sage.misc.instancedoc import InstanceDocDescriptor
        sage: def instancedoc(self):
        ....:     return "Instance doc"
        sage: docattr = InstanceDocDescriptor("Class doc", instancedoc)
        sage: class Z():
        ....:     __doc__ = InstanceDocDescriptor("Class doc", instancedoc)
        sage: Z.__doc__
        'Class doc'
        sage: Z().__doc__
        'Instance doc'

    We can still override the ``__doc__`` attribute of the instance::

        sage: obj = Z()
        sage: obj.__doc__ = "Custom doc"
        sage: obj.__doc__
        'Custom doc'
        sage: del obj.__doc__
        sage: obj.__doc__
        'Instance doc'
    """
    cdef classdoc
    cdef instancedoc
    cdef attr

    def __init__(self, classdoc, instancedoc, attr='__doc__'):
        """
        TESTS::

            sage: from sage.misc.instancedoc import InstanceDocDescriptor
            sage: InstanceDocDescriptor(None, None)
            <sage.misc.instancedoc.InstanceDocDescriptor object at ...>
        """
        self.classdoc = classdoc
        self.instancedoc = instancedoc
        self.attr = intern(attr)

    def __get__(self, obj, typ):
        """
        TESTS::

            sage: from sage.misc.instancedoc import InstanceDocDescriptor
            sage: def instancedoc(self):
            ....:     return "Doc for {!r}".format(self)
            sage: descr = InstanceDocDescriptor("Class doc", instancedoc)
            sage: descr.__get__(None, object)
            'Class doc'
            sage: descr.__get__(42, type(42))
            'Doc for 42'
        """
        if obj is None:
            return self.classdoc

        # First, try the attribute self.attr (typically __doc__)
        # on the instance
        try:
            objdict = obj.__dict__
        except AttributeError:
            pass
        else:
            doc = objdict.get(self.attr)
            if doc is not None:
                return doc

        return self.instancedoc(obj)

    def __set__(self, obj, value):
        """
        TESTS::

            sage: from sage.misc.instancedoc import InstanceDocDescriptor
            sage: def instancedoc(self):
            ....:     return "Doc for {!r}".format(self)
            sage: descr = InstanceDocDescriptor("Class doc", instancedoc)
            sage: class X(): pass
            sage: obj = X()
            sage: descr.__set__(obj, "Custom doc")
            sage: obj.__doc__
            'Custom doc'

            sage: descr.__set__([], "Custom doc")
            Traceback (most recent call last):
            ...
            AttributeError: attribute '__doc__' of 'list' objects is not writable
            sage: descr.__set__(object, "Custom doc")
            Traceback (most recent call last):
            ...
            AttributeError: attribute '__doc__' of 'type' objects is not writable
        """
        try:
            obj.__dict__[self.attr] = value
        except (AttributeError, TypeError):
            raise AttributeError(f"attribute '{self.attr}' of '{type(obj).__name__}' objects is not writable")

    def __delete__(self, obj):
        """
        TESTS::

            sage: from sage.misc.instancedoc import InstanceDocDescriptor
            sage: def instancedoc(self):
            ....:     return "Doc for {!r}".format(self)
            sage: descr = InstanceDocDescriptor("Class doc", instancedoc)
            sage: class X(): pass
            sage: obj = X()
            sage: obj.__doc__ = "Custom doc"
            sage: descr.__delete__(obj)
            sage: print(obj.__doc__)
            None
            sage: descr.__delete__(obj)
            Traceback (most recent call last):
            ...
            AttributeError: 'X' object has no attribute '__doc__'...

            sage: descr.__delete__([])
            Traceback (most recent call last):
            ...
            AttributeError: attribute '__doc__' of 'list' objects is not writable
            sage: descr.__delete__(object)
            Traceback (most recent call last):
            ...
            AttributeError: attribute '__doc__' of 'type' objects is not writable
        """
        try:
            del obj.__dict__[self.attr]
        except (AttributeError, TypeError):
            raise AttributeError(f"attribute '{self.attr}' of '{type(obj).__name__}' objects is not writable")
        except KeyError:
            raise AttributeError(f"'{type(obj).__name__}' object has no attribute '{self.attr}'")


def instancedoc(cls):
    """
    Add support for ``_instancedoc_`` to the class ``cls``.

    Typically, this will be used as decorator.

    INPUT:

    - ``cls`` -- a new-style class

    OUTPUT: ``cls``

    .. WARNING::

        ``instancedoc`` mutates the given class. So you are *not* supposed
        to use it as ``newcls = instancedoc(cls)`` because that would
        mutate ``cls`` (and ``newcls`` would be the same object as ``cls``)

    TESTS:

    We get a useful error message if ``_instancedoc_`` is not defined::

        sage: from sage.misc.instancedoc import instancedoc
        sage: class X(): pass
        sage: instancedoc(X)
        Traceback (most recent call last):
        ...
        TypeError: instancedoc requires <class '__main__.X'> to have an '_instancedoc_' attribute

    This does not work on old-style classes or things which are not a
    class at all::

        sage: instancedoc(7)
        Traceback (most recent call last):
        ...
        TypeError: expected type, got 7

        sage: class OldStyle: pass
        sage: instancedoc(OldStyle)
        Traceback (most recent call last):
        ...
        TypeError: instancedoc requires <class '__main__.OldStyle'> to have an '_instancedoc_' attribute
    """
    cdef PyTypeObject* tp = TypeObject(cls)
    try:
        instdoc = cls._instancedoc_
    except AttributeError:
        raise TypeError(f"instancedoc requires {cls!r} to have an '_instancedoc_' attribute")
    docattr = InstanceDocDescriptor(cls.__doc__, instdoc)
    PyDict_SetItemString(tp.tp_dict, "__doc__", docattr)
    tp.tp_doc = NULL
    PyType_Modified(tp)
    return cls
