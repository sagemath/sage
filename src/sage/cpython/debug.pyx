"""
Various functions to debug Python internals
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport PyObject, Py_TYPE, descrgetfunc, descrsetfunc

cdef extern from "Python.h":
    # Helper to get a pointer to an object's __dict__ slot, if any
    PyObject** _PyObject_GetDictPtr(obj)

from .getattr cimport AttributeErrorMessage


def shortrepr(obj, max=50):
    """
    Return ``repr(obj)`` bounded to ``max`` characters. If the string
    is too long, it is truncated and ``~~~`` is added to the end.

    EXAMPLES::

        sage: from sage.cpython.debug import shortrepr
        sage: print(shortrepr("Hello world!"))
        'Hello world!'
        sage: print(shortrepr("Hello world!" * 4))
        'Hello world!Hello world!Hello world!Hello world!'
        sage: print(shortrepr("Hello world!" * 5))
        'Hello world!Hello world!Hello world!Hello worl~~~
    """
    r = repr(obj)
    if len(r) > max:
        r = r[:max-3] + "~~~"
    return r


cdef object _no_default = object()  # Unique object

def getattr_debug(obj, name, default=_no_default):
    r"""
    A re-implementation of ``getattr()`` with lots of debugging info.

    This will correctly use ``__getattr__`` if needed. On the other
    hand, it assumes a generic (not overridden) implementation of
    ``__getattribute__``. Note that Cython implements ``__getattr__``
    for a cdef class using ``__getattribute__``, so this will not
    detect a ``__getattr__`` in that case.

    INPUT:

    - ``obj`` -- the object whose attribute is requested

    - ``name`` -- (string) the name of the attribute

    - ``default`` -- default value to return if attribute was not found

    EXAMPLES::

        sage: _ = getattr_debug(list, "reverse")
        getattr_debug(obj=<type 'list'>, name='reverse'):
          type(obj) = <type 'type'>
          object has __dict__ slot (<type 'dict'>)
          did not find 'reverse' in MRO classes
          found 'reverse' in object __dict__
          returning <method 'reverse' of 'list' objects> (<type 'method_descriptor'>)
        sage: _ = getattr_debug([], "reverse")
        getattr_debug(obj=[], name='reverse'):
          type(obj) = <type 'list'>
          object does not have __dict__ slot
          found 'reverse' in dict of <type 'list'>
          got <method 'reverse' of 'list' objects> (<type 'method_descriptor'>)
          attribute is ordinary descriptor (has __get__)
          calling __get__()
          returning <built-in method reverse of list object at 0x... (<type 'builtin_function_or_method'>)
        sage: _ = getattr_debug([], "__doc__")
        getattr_debug(obj=[], name='__doc__'):
          type(obj) = <type 'list'>
          object does not have __dict__ slot
          found '__doc__' in dict of <type 'list'>
          got "list() -> new empty list\nlist(iterable) -> ne~~~ (<type 'str'>)
          returning "list() -> new empty list\nlist(iterable) -> ne~~~ (<type 'str'>)
        sage: _ = getattr_debug(gp(1), "log")
        getattr_debug(obj=1, name='log'):
          type(obj) = <class 'sage.interfaces.gp.GpElement'>
          object has __dict__ slot (<type 'dict'>)
          did not find 'log' in MRO classes
          object __dict__ does not have 'log'
          calling __getattr__()
          returning log (<class 'sage.interfaces.expect.FunctionElement'>)
        sage: from ipywidgets import IntSlider
        sage: _ = getattr_debug(IntSlider(), "value")
        getattr_debug(obj=IntSlider(value=0, min=0, max=100, step=1), name='value'):
          type(obj) = <class 'ipywidgets.widgets.widget_int.IntSlider'>
          object has __dict__ slot (<type 'dict'>)
          found 'value' in dict of <class 'ipywidgets.widgets.widget_int._Int'>
          got <traitlets.traitlets.CInt object at ... (<class 'traitlets.traitlets.CInt'>)
          attribute is data descriptor (has __get__ and __set__)
          ignoring __dict__ because we have a data descriptor
          calling __get__()
          returning 0 (<type 'int'>)
        sage: _ = getattr_debug(1, "foo")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'foo'
        sage: _ = getattr_debug(1, "foo", "xyz")
        getattr_debug(obj=1, name='foo'):
          type(obj) = <type 'sage.rings.integer.Integer'>
          object does not have __dict__ slot
          did not find 'foo' in MRO classes
          class does not have __getattr__
          attribute not found
          returning default 'xyz'
    """
    if default is not _no_default:
        try:
            return getattr_debug(obj, name)
        except AttributeError:
            print(f"  returning default {shortrepr(default)}")
            return default

    name = str(name)
    print(f"getattr_debug(obj={shortrepr(obj)}, name={name!r}):")
    print(f"  type(obj) = {type(obj)}")

    cdef object attr = None
    cdef object dct = None

    cdef PyObject** dictptr = _PyObject_GetDictPtr(obj)
    if dictptr is not NULL:
        if dictptr[0] is not NULL:
            dct = <object>(dictptr[0])
            print(f"  object has __dict__ slot ({type(dct)})")
        else:
            print(f"  object has uninitialized __dict__ slot")
    else:
        print(f"  object does not have __dict__ slot")

    cdef descrgetfunc get = NULL
    cdef descrsetfunc set = NULL

    # Look for name in dicts of types in MRO
    cdef bint attr_in_class = False
    for cls in type(obj).__mro__:
        if name in cls.__dict__:
            attr = cls.__dict__[name]
            print(f"  found {name!r} in dict of {cls}")
            print(f"  got {shortrepr(attr)} ({type(attr)})")

            get = Py_TYPE(attr).tp_descr_get
            set = Py_TYPE(attr).tp_descr_set
            if get is not NULL:
                if set is not NULL:
                    print(f"  attribute is data descriptor (has __get__ and __set__)")
                    if dct is not None:
                        print(f"  ignoring __dict__ because we have a data descriptor")
                    print(f"  calling __get__()")
                    attr = get(attr, obj, type(obj))
                    print(f"  returning {shortrepr(attr)} ({type(attr)})")
                    return attr
                else:
                    print(f"  attribute is ordinary descriptor (has __get__)")
            attr_in_class = True
            break

    if not attr_in_class:
        print(f"  did not find {name!r} in MRO classes")

    if dct is not None:
        if name in dct:
            attr = dct[name]
            print(f"  found {name!r} in object __dict__")
            print(f"  returning {shortrepr(attr)} ({type(attr)})")
            return attr
        else:
            print(f"  object __dict__ does not have {name!r}")

    if attr_in_class:
        if get is not NULL:
            print(f"  calling __get__()")
            attr = get(attr, obj, type(obj))
        print(f"  returning {shortrepr(attr)} ({type(attr)})")
        return attr

    try:
        tpgetattr = type(obj).__getattr__
    except AttributeError:
        print(f"  class does not have __getattr__")
    else:
        print(f"  calling __getattr__()")
        attr = tpgetattr(obj, name)
        print(f"  returning {shortrepr(attr)} ({type(attr)})")
        return attr

    print(f"  attribute not found")
    raise AttributeError(AttributeErrorMessage(obj, name))
