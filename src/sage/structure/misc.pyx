"""
Miscellaneous utilities
"""

from cpython.object cimport PyObject, Py_TYPE, descrgetfunc

cdef extern from "Python.h":
    # Internal API to look for a name through the MRO.
    # This returns a borrowed reference, and doesn't set an exception!
    PyObject* _PyType_Lookup(type t, name)


def is_extension_type(cls):
    """
    INPUT:

    - cls: a class

    Tests whether cls is an extension type (int, list, cython compiled classes, ...)

    EXAMPLES::

        sage: from sage.structure.parent import is_extension_type
        sage: is_extension_type(int)
        True
        sage: is_extension_type(list)
        True
        sage: is_extension_type(ZZ.__class__)
        True
        sage: is_extension_type(QQ.__class__)
        False
    """
    # Robert B claims that this should be robust
    try:
        return cls.__dictoffset__ == 0
    except AttributeError:
        pass
    return False


cdef class AttributeErrorMessage:
    """
    Tries to emulate the standard Python ``AttributeError`` message.

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
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'bla'
        sage: QQ[x].gen().bla
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint' object has no attribute 'bla'
        sage: from sage.structure.misc import AttributeErrorMessage
        sage: AttributeErrorMessage(int(1),'bla')
        'int' object has no attribute 'bla'

    TESTS:

    By :trac:`14100`, the attribute errors raised on elements and parents are
    unique objects. The error message of this unique error object is changed
    inplace. This is for reasons of efficiency. ::

        sage: try:
        ....:     1.__bla
        ....: except AttributeError as ElementError:
        ....:     pass
        sage: ElementError
        AttributeError('sage.rings.integer.Integer' object has no attribute '__bla',)
        sage: try:
        ....:     x.__bla
        ....: except AttributeError as ElementError2:
        ....:     pass
        sage: ElementError2 is ElementError
        True
        sage: ElementError
        AttributeError('sage.symbolic.expression.Expression' object has no attribute '__bla',)
        sage: isinstance(ElementError.args[0], sage.structure.misc.AttributeErrorMessage)
        True

    Hence, if one really needs the error message as a string, then one should
    make a copy of its string representation before it changes. Attribute
    Errors of parents behave similarly::

        sage: try:
        ....:     QQ.__bla
        ....: except AttributeError as ParentError:
        ....:     pass
        sage: ParentError
        AttributeError('RationalField_with_category' object has no attribute '__bla',)
        sage: try:
        ....:     ZZ.__bla
        ....: except AttributeError as ParentError2:
        ....:     pass
        sage: ParentError2 is ParentError
        True
        sage: ParentError2
        AttributeError('sage.rings.integer_ring.IntegerRing_class' object has no attribute '__bla',)
        sage: ParentError2 is ElementError
        True

    AUTHOR:

    - Simon King (2011-05-21)
    """
    def __init__(self, P, name):
        """
        INPUT:

        - ``P`` -- any object

        - ``name`` -- a string

        TESTS::

            sage: from sage.structure.misc import AttributeErrorMessage
            sage: AttributeErrorMessage(int(1), 'bla')
            'int' object has no attribute 'bla'
        """
        self.cls = type(P)
        self.name = name

    def __repr__(self):
        """
        TESTS::

            sage: from sage.structure.misc import AttributeErrorMessage
            sage: AttributeErrorMessage(int(1), 'bla')
            'int' object has no attribute 'bla'
        """
        cdef int dictoff = 1
        try:
            dictoff = self.cls.__dictoffset__
        except AttributeError:
            pass
        if dictoff:
            cls = repr(self.cls.__name__)
        else:
            cls = repr(self.cls)[6:-1]
        return f"{cls} object has no attribute {self.name!r}"


cdef AttributeErrorMessage dummy_error_message = AttributeErrorMessage(None, '')
dummy_attribute_error = AttributeError(dummy_error_message)


cpdef getattr_from_other_class(self, cls, name):
    """
    Emulate ``getattr(self, name)``, as if ``self`` was an instance of
    ``cls``.

    INPUT:

    - ``self`` -- some object

    - ``cls`` -- a new-style class

    - ``name`` -- a string

    If self is an instance of cls, raises an ``AttributeError``, to
    avoid a double lookup. This function is intended to be called from
    __getattr__, and so should not be called if name is an attribute
    of self.

    EXAMPLES::

        sage: from sage.structure.misc import getattr_from_other_class
        sage: class A(object):
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
    ``__cached_methods`` of type ``<dict>``. This condition
    is satisfied, e.g., by any class that is derived from
    :class:`Parent`::

        sage: getattr_from_other_class(1, A, "lazy_attribute")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'lazy_attribute'

    The integer ring is a parent, so, lazy attributes work::

        sage: getattr_from_other_class(ZZ, A, "lazy_attribute")
        'Integer Ring'
        sage: getattr_from_other_class(PolynomialRing(QQ, name='x', sparse=True).one(), A, "lazy_attribute")
        '1'
        sage: getattr_from_other_class(PolynomialRing(QQ, name='x', implementation="FLINT").one(), A, "lazy_attribute")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'
        object has no attribute 'lazy_attribute'

    In general, descriptors are not yet well supported, because they
    often do not accept to be cheated with the type of their instance::

        sage: A.__weakref__.__get__(1)
        Traceback (most recent call last):
        ...
        TypeError: descriptor '__weakref__' for 'A' objects doesn't apply
        to 'sage.rings.integer.Integer' object

    When this occurs, an ``AttributeError`` is raised::

        sage: getattr_from_other_class(1, A, "__weakref__")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__weakref__'

    This was caught by :trac:`8296` for which we do a couple more tests::

        sage: "__weakref__" in dir(A)
        True
        sage: "__weakref__" in dir(1)
        True
        sage: 1.__weakref__
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__weakref__'

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
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__call__'

    TESTS:

    Check that we do not pick up special attributes from the ``type``
    class, see :trac:`20686`::

        sage: getattr_from_other_class(1, type, "__name__")
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute '__name__'

    This does not work with an old-style class::

        sage: class OldStyle:
        ....:     pass
        sage: getattr_from_other_class(1, OldStyle, "foo")
        Traceback (most recent call last):
        ...
        TypeError: <class __main__.OldStyle at ...> is not a type

    Non-strings as "name" are handled gracefully::

        sage: getattr_from_other_class(1, type, None)
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute None
    """
    if not isinstance(cls, type):
        raise TypeError(f"{cls!r} is not a type")
    if isinstance(self, cls):
        dummy_error_message.cls = type(self)
        dummy_error_message.name = name
        raise dummy_attribute_error
    cdef PyObject* attr = _PyType_Lookup(<type>cls, name)
    if attr is NULL:
        dummy_error_message.cls = type(self)
        dummy_error_message.name = name
        raise dummy_attribute_error
    attribute = <object>attr
    # Check for a descriptor (__get__ in Python)
    cdef descrgetfunc getter = Py_TYPE(attribute).tp_descr_get
    if getter is NULL:
        # Not a descriptor
        return attribute
    # Conditionally defined lazy_attributes don't work well with fake subclasses
    # (a TypeError is raised if the lazy attribute is not defined).
    # For the moment, we ignore that when this occurs.
    # Other descriptors (including __weakref__) also break.
    try:
        return getter(attribute, self, cls)
    except TypeError:
        pass
    dummy_error_message.cls = type(self)
    dummy_error_message.name = name
    raise dummy_attribute_error


def dir_with_other_class(self, cls):
    r"""
    Emulates ``dir(self)``, as if self was also an instance ``cls``,
    right after ``caller_class`` in the method resolution order
    (``self.__class__.mro()``)

    EXAMPLES::

        sage: class A(object):
        ....:    a = 1
        ....:    b = 2
        ....:    c = 3
        sage: class B(object):
        ....:    b = 2
        ....:    c = 3
        ....:    d = 4
        sage: x = A()
        sage: x.c = 1; x.e = 1
        sage: from sage.structure.misc import dir_with_other_class
        sage: dir_with_other_class(x, B)
        [..., 'a', 'b', 'c', 'd', 'e']

    Check that objects without dicts are well handled::

        sage: cython("cdef class A:\n    cdef public int a")
        sage: cython("cdef class B:\n    cdef public int b")
        sage: x = A()
        sage: x.a = 1
        sage: hasattr(x,'__dict__')
        False
        sage: dir_with_other_class(x, B)
        [..., 'a', 'b']

    TESTS:

    Check that :trac:`13043` is fixed::

        sage: len(dir(RIF))==len(set(dir(RIF)))
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
    if not isinstance(self, cls):
        ret.update(dir(cls))
    return sorted(ret)
