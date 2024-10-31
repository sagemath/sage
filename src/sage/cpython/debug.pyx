# sage_setup: distribution = sagemath-objects
"""
Various functions to debug Python internals
"""
# ****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport (PyObject, PyTypeObject, Py_TYPE,
                             descrgetfunc, descrsetfunc)

cdef extern from "Python.h":
    # Helper to get a pointer to an object's __dict__ slot, if any
    PyObject** _PyObject_GetDictPtr(obj)

cdef extern from "debugimpl.c":
    void _type_debug(PyTypeObject*)

from sage.cpython.getattr cimport AttributeErrorMessage


# Determine subtype_traverse, subtype_clear, subtype_dealloc functions
# for type_debug(). These are the default tp_traverse, tp_clear and
# tp_dealloc functions for heap types (= Python classes).
cdef:
    X = type("X", (), {})  # heap type
    void* subtype_traverse "subtype_traverse" = (<PyTypeObject*>X).tp_traverse
    void* subtype_clear "subtype_clear" = (<PyTypeObject*>X).tp_clear
    void* subtype_dealloc "subtype_dealloc" = (<PyTypeObject*>X).tp_dealloc


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

    - ``name`` -- string; the name of the attribute

    - ``default`` -- default value to return if attribute was not found

    EXAMPLES::

        sage: _ = getattr_debug(list, "reverse")  # not tested - broken in python 3.12
        getattr_debug(obj=<class 'list'>, name='reverse'):
          type(obj) = <class 'type'>
          object has __dict__ slot (<class 'dict'>)
          did not find 'reverse' in MRO classes
          found 'reverse' in object __dict__
          returning <method 'reverse' of 'list' objects> (<class 'method_descriptor'>)
        sage: _ = getattr_debug([], "reverse")
        getattr_debug(obj=[], name='reverse'):
          type(obj) = <class 'list'>
          object does not have __dict__ slot
          found 'reverse' in dict of <class 'list'>
          got <method 'reverse' of 'list' objects> (<class 'method_descriptor'>)
          attribute is ordinary descriptor (has __get__)
          calling __get__()
          returning <built-in method reverse of list object at 0x... (<class 'builtin_function_or_method'>)
        sage: _ = getattr_debug([], "__doc__")
        getattr_debug(obj=[], name='__doc__'):
          type(obj) = <class 'list'>
          object does not have __dict__ slot
          found '__doc__' in dict of <class 'list'>
          got ... 'str'>)
          returning ... 'str'>)
        sage: _ = getattr_debug(gp(1), "log")                                           # needs sage.libs.pari
        getattr_debug(obj=1, name='log'):
          type(obj) = <class 'sage.interfaces.gp.GpElement'>
          object has __dict__ slot (<class 'dict'>)
          did not find 'log' in MRO classes
          object __dict__ does not have 'log'
          calling __getattr__()
          returning log (<class 'sage.interfaces.expect.FunctionElement'>)
        sage: from ipywidgets import IntSlider
        sage: _ = getattr_debug(IntSlider(), "value")
        getattr_debug(obj=IntSlider(value=0), name='value'):
          type(obj) = <class 'ipywidgets.widgets.widget_int.IntSlider'>
          object has __dict__ slot (<class 'dict'>)
          found 'value' in dict of <class 'ipywidgets.widgets.widget_int._Int'>
          got <traitlets.traitlets.CInt object at ... (<class 'traitlets.traitlets.CInt'>)
          attribute is data descriptor (has __get__ and __set__)
          ignoring __dict__ because we have a data descriptor
          calling __get__()
          returning 0 (<class 'int'>)
        sage: _ = getattr_debug(1, "foo")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'foo'...
        sage: _ = getattr_debug(1, "foo", "xyz")
        getattr_debug(obj=1, name='foo'):
          type(obj) = <class 'sage.rings.integer.Integer'>
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
            print("  object has uninitialized __dict__ slot")
    else:
        print("  object does not have __dict__ slot")

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
                    print("  attribute is data descriptor (has __get__ and __set__)")
                    if dct is not None:
                        print("  ignoring __dict__ because we have a data descriptor")
                    print("  calling __get__()")
                    attr = get(attr, obj, type(obj))
                    print(f"  returning {shortrepr(attr)} ({type(attr)})")
                    return attr
                else:
                    print("  attribute is ordinary descriptor (has __get__)")
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
            print("  calling __get__()")
            attr = get(attr, obj, type(obj))
        print(f"  returning {shortrepr(attr)} ({type(attr)})")
        return attr

    try:
        tpgetattr = type(obj).__getattr__
    except AttributeError:
        print("  class does not have __getattr__")
    else:
        print("  calling __getattr__()")
        attr = tpgetattr(obj, name)
        print(f"  returning {shortrepr(attr)} ({type(attr)})")
        return attr

    print("  attribute not found")
    raise AttributeError(AttributeErrorMessage(obj, name))


def type_debug(cls):
    """
    Print all internals of the type ``cls``.

    EXAMPLES::

        sage: type_debug(object)  # random
        <class 'object'> (0x7fc57da7f040)
          ob_refcnt: 9739
          ob_type: <class 'type'>
          tp_name: object
          tp_basicsize: 16
          tp_itemsize: 0
          tp_dictoffset: 0
          tp_weaklistoffset: 0
          tp_base (__base__): NULL
          tp_bases (__bases__): tuple:
          tp_mro (__mro__): tuple:
            <class 'object'>
          tp_dict (__dict__): dict:
            '__setattr__': <slot wrapper '__setattr__' of 'object' objects>
            '__reduce_ex__': <method '__reduce_ex__' of 'object' objects>
            '__new__': <built-in method __new__ of type object at 0x7fc57da7f040>
            '__reduce__': <method '__reduce__' of 'object' objects>
            '__str__': <slot wrapper '__str__' of 'object' objects>
            '__format__': <method '__format__' of 'object' objects>
            '__getattribute__': <slot wrapper '__getattribute__' of 'object' objects>
            '__class__': <attribute '__class__' of 'object' objects>
            '__delattr__': <slot wrapper '__delattr__' of 'object' objects>
            '__subclasshook__': <method '__subclasshook__' of 'object' objects>
            '__repr__': <slot wrapper '__repr__' of 'object' objects>
            '__hash__': <slot wrapper '__hash__' of 'object' objects>
            '__sizeof__': <method '__sizeof__' of 'object' objects>
            '__doc__': 'The most base type'
            '__init__': <slot wrapper '__init__' of 'object' objects>
          tp_alloc: PyType_GenericAlloc
          tp_new (__new__): 0x7fc57d7594f0
          tp_init (__init__): 0x7fc57d758ee0
          tp_dealloc (__dealloc__): 0x7fc57d757010
          tp_free: PyObject_Del
          tp_repr (__repr__): 0x7fc57d75b990
          tp_print: NULL
          tp_hash (__hash__): _Py_HashPointer
          tp_call (__call__): NULL
          tp_str (__str__): 0x7fc57d757020
          tp_compare (__cmp__): NULL
          tp_richcompare (__richcmp__): NULL
          tp_getattr (__getattribute__): NULL
          tp_setattr (__setattribute__): NULL
          tp_getattro (__getattribute__): PyObject_GenericGetAttr
          tp_setattro (__setattribute__): PyObject_GenericSetAttr
          tp_iter (__iter__): NULL
          tp_iternext (__next__): NULL
          tp_descr_get (__get__): NULL
          tp_descr_set (__set__): NULL
          tp_cache: NULL
          tp_weaklist: NULL
          tp_traverse: NULL
          tp_clear: NULL
          tp_is_gc: NULL
          tp_as_number: NULL
          tp_as_sequence: NULL
          tp_as_mapping: NULL
          tp_as_buffer: NULL
          tp_flags:
            HAVE_GETCHARBUFFER
            HAVE_SEQUENCE_IN
            HAVE_INPLACEOPS
            HAVE_RICHCOMPARE
            HAVE_WEAKREFS
            HAVE_ITER
            HAVE_CLASS
            BASETYPE
            READY
            HAVE_INDEX
            HAVE_VERSION_TAG
            VALID_VERSION_TAG
          tp_version_tag: 2
        sage: type_debug(None)
        Traceback (most recent call last):
        ...
        TypeError: None is not a type
    """
    if not isinstance(cls, type):
        raise TypeError(f"{cls!r} is not a type")

    _type_debug(<PyTypeObject*>cls)
