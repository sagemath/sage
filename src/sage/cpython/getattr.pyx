# sage_setup: distribution = sagemath-objects
"""
Variants of getattr()
"""

from cpython.object cimport PyObject, PyTypeObject, Py_TYPE, descrgetfunc

from sage.cpython.string cimport bytes_to_str

cdef extern from "Python.h":
    r"""
    /* Coded in C because class internals need to be accessed. */
    static PyObject*
    instance_getattr(PyObject* obj, PyObject* name)
    {
        if (PyType_Check(obj)) {
            return _PyType_Lookup((PyTypeObject*)obj, name);
        }

        PyObject** dptr = _PyObject_GetDictPtr(obj);
        if (dptr == NULL) return NULL;
        PyObject* dict = *dptr;
        if (dict == NULL) return NULL;
        return PyDict_GetItem(dict, name);
    }
    """

    # Return the attribute "name" from "obj". This only looks up the
    # attribute in the instance "obj" and not in obj.__class__.
    # If "obj" is a class, this searches for the attribute in the base
    # classes.
    #
    # Return a borrowed reference or NULL if the attribute was not found.
    PyObject* instance_getattr(obj, name)

    int PyDescr_IsData(PyObject*)


cdef class AttributeErrorMessage:
    """
    Try to emulate the standard Python :exc:`AttributeError` message.

    .. NOTE::

        The typical fate of an attribute error is being caught. Hence,
        under normal circumstances, nobody will ever see the error
        message. The idea for this class is to provide an object that
        is fast to create and whose string representation is an attribute
        error's message. That string representation is only created if
        someone wants to see it.

    EXAMPLES::

        sage: 1.bla  #indirect doctest
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'bla'...
        sage: x = polygen(ZZ, 'x')
        sage: QQ[x].gen().bla                                                           # needs sage.libs.flint
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint' object has no attribute 'bla'...

    ::

        sage: from sage.cpython.getattr import AttributeErrorMessage
        sage: AttributeErrorMessage(int(1), 'bla')
        'int' object has no attribute 'bla'

    TESTS:

    The error message used for the :exc:`AttributeError` is a unique object
    and is changed inplace. This is for reasons of efficiency.
    Hence, if one really needs the error message as a string, then one should
    make a copy of its string representation before it changes. ::

        sage: try:
        ....:     1.__bla
        ....: except AttributeError as exc:
        ....:     ElementError = exc
        sage: ElementError
        AttributeError('sage.rings.integer.Integer' object has no attribute '__bla'...)
        sage: try:
        ....:     x.__bla
        ....: except AttributeError as exc:
        ....:     ElementError2 = exc
        sage: ElementError
        AttributeError('sage.rings.polynomial...' object has no attribute '__bla'...)
        sage: ElementError2.args[0] is ElementError.args[0]
        True
        sage: isinstance(ElementError.args[0], sage.cpython.getattr.AttributeErrorMessage)
        True

    AUTHOR:

    - Simon King (2011-05-21)
    """
    def __init__(self, obj=None, name=""):
        self.cls = type(obj)
        self.name = name

    def __repr__(self):
        cls = bytes_to_str((<PyTypeObject*>self.cls).tp_name, 'utf-8',
                           'replace')
        # Go directly through tp_name since __name__ can be overridden--this is
        # almost verbatim how CPython formats this message except we don't cut
        # off the class name after 50 characters, and non-strings are displayed
        # with their repr :)
        return f"'{cls}' object has no attribute {self.name!r}"


cdef AttributeErrorMessage dummy_error_message = AttributeErrorMessage()


cpdef raw_getattr(obj, name):
    """
    Like ``getattr(obj, name)`` but without invoking the binding
    behavior of descriptors under normal attribute access.
    This can be used to easily get unbound methods or other
    descriptors.

    This ignores ``__getattribute__`` hooks but it does support
    ``__getattr__``.

    .. NOTE::

        For Cython classes, ``__getattr__`` is actually implemented as
        ``__getattribute__``, which means that it is not supported by
        ``raw_getattr``.

    EXAMPLES::

        sage: class X:
        ....:     @property
        ....:     def prop(self):
        ....:         return 42
        ....:     def method(self):
        ....:         pass
        ....:     def __getattr__(self, name):
        ....:         return "magic " + name
        sage: raw_getattr(X, "prop")
        <property object at ...>
        sage: raw_getattr(X, "method")
        <function ...method at ...>
        sage: raw_getattr(X, "attr")
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'attr'...
        sage: x = X()
        sage: raw_getattr(x, "prop")
        <property object at ...>
        sage: raw_getattr(x, "method")
        <function ...method at ...>
        sage: raw_getattr(x, "attr")
        'magic attr'
        sage: x.__dict__["prop"] = 'no'
        sage: x.__dict__["method"] = 'yes'
        sage: x.__dict__["attr"] = 'ok'
        sage: raw_getattr(x, "prop")
        <property object at ...>
        sage: raw_getattr(x, "method")
        'yes'
        sage: raw_getattr(x, "attr")
        'ok'

    The same tests with an inherited new-style class::

        sage: class Y(X, object):
        ....:     pass
        sage: raw_getattr(Y, "prop")
        <property object at ...>
        sage: raw_getattr(Y, "method")
        <function ...method at ...>
        sage: raw_getattr(Y, "attr")
        Traceback (most recent call last):
        ...
        AttributeError: '...' object has no attribute 'attr'...
        sage: y = Y()
        sage: raw_getattr(y, "prop")
        <property object at ...>
        sage: raw_getattr(y, "method")
        <function ...method at ...>
        sage: raw_getattr(y, "attr")
        'magic attr'
        sage: y.__dict__["prop"] = 'no'
        sage: y.__dict__["method"] = 'yes'
        sage: y.__dict__["attr"] = 'ok'
        sage: raw_getattr(y, "prop")
        <property object at ...>
        sage: raw_getattr(y, "method")
        'yes'
        sage: raw_getattr(y, "attr")
        'ok'
    """
    cdef PyObject* class_attr = NULL

    try:
        cls = obj.__class__
    except AttributeError:
        # Old-style classes don't have a __class__ (because these do
        # not support metaclasses).
        cls = None
    else:
        class_attr = instance_getattr(cls, name)

        # We honor the order prescribed by the descriptor protocol:
        # a data descriptor overrides instance attributes. In other
        # cases, the instance attribute takes priority.
        if class_attr is not NULL and PyDescr_IsData(class_attr):
            return <object>class_attr

    instance_attr = instance_getattr(obj, name)
    if instance_attr is not NULL:
        return <object>instance_attr
    if class_attr is not NULL:
        return <object>class_attr

    if cls is not None:
        try:
            cls_getattr = cls.__getattr__
        except AttributeError:
            pass
        else:
            return cls_getattr(obj, name)

    dummy_error_message.cls = type(obj)
    dummy_error_message.name = name
    raise AttributeError(dummy_error_message)


cpdef getattr_from_other_class(self, cls, name):
    """
    Emulate ``getattr(self, name)``, as if ``self`` was an instance of
    ``cls``.

    INPUT:

    - ``self`` -- some object

    - ``cls`` -- a new-style class

    - ``name`` -- string

    If ``self`` is an instance of cls, raises an :exc:`AttributeError`, to
    avoid a double lookup. This function is intended to be called from
    __getattr__, and so should not be called if name is an attribute
    of ``self``.

    EXAMPLES::

        sage: from sage.cpython.getattr import getattr_from_other_class
        sage: class A():
        ....:      def inc(self):
        ....:          return self + 1
        ....:
        ....:      @staticmethod
        ....:      def greeting():
        ....:          print("Hello World!")
        ....:
        ....:      @lazy_attribute
        ....:      def lazy_attribute(self):
        ....:          return repr(self)
        sage: getattr_from_other_class(1, A, "inc")
        <bound method A.inc of 1>
        sage: getattr_from_other_class(1, A, "inc")()
        2

    Static methods work::

        sage: getattr_from_other_class(1, A, "greeting")()
        Hello World!

    Caveat: lazy attributes work with extension types only
    if they allow attribute assignment or have a public attribute
    ``_cached_methods`` of type ``<dict>``. This condition
    is satisfied, e.g., by any class that is derived from
    :class:`Parent`::

        sage: getattr_from_other_class(1, A, "lazy_attribute")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'lazy_attribute'...

    The integer ring is a parent, so, lazy attributes work::

        sage: getattr_from_other_class(ZZ, A, "lazy_attribute")
        'Integer Ring'
        sage: getattr_from_other_class(PolynomialRing(QQ, name='x', sparse=True).one(), A, "lazy_attribute")
        '1'
        sage: getattr_from_other_class(17, A, "lazy_attribute")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'lazy_attribute'...

    In general, descriptors are not yet well supported, because they
    often do not accept to be cheated with the type of their instance::

        sage: A.__weakref__.__get__(1)
        Traceback (most recent call last):
        ...
        TypeError: descriptor '__weakref__' for 'A' objects doesn't apply
        to ...'sage.rings.integer.Integer' object

    When this occurs, an :exc:`AttributeError` is raised::

        sage: getattr_from_other_class(1, A, "__weakref__")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__weakref__'...

    This was caught by :issue:`8296` for which we do a couple more tests::

        sage: "__weakref__" in dir(A)
        True
        sage: 1.__weakref__
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__weakref__'...

        sage: n = 1
        sage: ip = get_ipython()                 # not tested: only works in interactive shell
        sage: ip.magic_psearch('n.N')            # not tested: only works in interactive shell
        n.N
        sage: ip.magic_psearch('n.__weakref__')  # not tested: only works in interactive shell

    Caveat: When __call__ is not defined for instances, using
    ``A.__call__`` yields the method ``__call__`` of the class. We use
    a workaround but there is no guarantee for robustness.

        sage: getattr_from_other_class(1, A, "__call__")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__call__'...

    TESTS:

    Check that we do not pick up special attributes from the ``type``
    class, see :issue:`20686`::

        sage: getattr_from_other_class(1, type, "__name__")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__name__'...

    Non-strings as "name" are handled gracefully::

        sage: getattr_from_other_class(1, type, None)
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute None...
    """
    if not isinstance(cls, type):
        raise TypeError(f"{cls!r} is not a type")

    if isinstance(self, cls):
        dummy_error_message.cls = type(self)
        dummy_error_message.name = name
        raise AttributeError(dummy_error_message)
    cdef PyObject* attr = instance_getattr(cls, name)
    if attr is NULL:
        dummy_error_message.cls = type(self)
        dummy_error_message.name = name
        raise AttributeError(dummy_error_message)
    attribute = <object>attr
    # Check for a descriptor (__get__ in Python)
    cdef descrgetfunc getter = Py_TYPE(attribute).tp_descr_get
    if getter is NULL:
        # Not a descriptor
        return attribute
    # Conditionally defined lazy_attributes don't work well with fake subclasses
    # (a :exc:`TypeError` is raised if the lazy attribute is not defined).
    # For the moment, we ignore that when this occurs.
    # Other descriptors (including __weakref__) also break.
    try:
        return getter(attribute, self, cls)
    except TypeError:
        pass
    dummy_error_message.cls = type(self)
    dummy_error_message.name = name
    raise AttributeError(dummy_error_message)


def dir_with_other_class(self, *cls):
    r"""
    Emulates ``dir(self)``, as if ``self`` was also an instance ``cls``,
    right after ``caller_class`` in the method resolution order
    (``self.__class__.mro()``)

    EXAMPLES::

        sage: class A():
        ....:    a = 1
        ....:    b = 2
        ....:    c = 3
        sage: class B():
        ....:    b = 2
        ....:    c = 3
        ....:    d = 4
        sage: x = A()
        sage: x.c = 1; x.e = 1
        sage: from sage.cpython.getattr import dir_with_other_class
        sage: dir_with_other_class(x, B)
        [..., 'a', 'b', 'c', 'd', 'e']
        sage: class C():
        ....:    f = 6
        sage: dir_with_other_class(x, B, C)
        [..., 'a', 'b', 'c', 'd', 'e', 'f']

    Check that objects without dicts are well handled::

        sage: # needs sage.misc.cython
        sage: cython("cdef class A:\n    cdef public int a")
        sage: cython("cdef class B:\n    cdef public int b")
        sage: x = A()
        sage: x.a = 1
        sage: hasattr(x,'__dict__')
        False
        sage: dir_with_other_class(x, B)
        [..., 'a', 'b']

    TESTS:

    Check that :issue:`13043` is fixed::

        sage: len(dir(RIF))==len(set(dir(RIF)))                                         # needs sage.rings.real_interval_field
        True
    """
    ret = set()
    # This tries to emulate the standard dir function
    # Is there a better way to call dir on self, while ignoring this
    # __dir__? Using dir(super(A, self)) does not work since the
    # attributes coming from subclasses of A will be ignored
    ret.update(dir(self.__class__))
    if hasattr(self, "__dict__"):
        ret.update(list(self.__dict__))
    for c in cls:
        if not isinstance(self, c):
            ret.update(dir(c))
    return sorted(ret)
