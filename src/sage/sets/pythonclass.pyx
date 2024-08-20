# sage_setup: distribution = sagemath-objects
"""
Set of all objects of a given Python class
"""

#*****************************************************************************
#       Copyright (C) 2018 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport Py_EQ, Py_NE
from sage.structure.richcmp cimport rich_to_bool
from sage.categories.sets_cat import Sets


cdef dict _type_set_cache = {}

cpdef Set_PythonType(typ):
    """
    Return the (unique) Parent that represents the set of Python objects
    of a specified type.

    EXAMPLES::

        sage: from sage.sets.pythonclass import Set_PythonType
        sage: Set_PythonType(list)
        Set of Python objects of class 'list'
        sage: Set_PythonType(list) is Set_PythonType(list)
        True
        sage: S = Set_PythonType(tuple)
        sage: S([1,2,3])
        (1, 2, 3)

    S is a parent which models the set of all lists::

        sage: S.category()
        Category of sets
    """
    try:
        return _type_set_cache[typ]
    except KeyError:
        _type_set_cache[typ] = theSet = Set_PythonType_class(typ)
        return theSet


cdef class Set_PythonType_class(Set_generic):
    r"""
    The set of Python objects of a given class.

    The elements of this set are not instances of
    :class:`~sage.structure.element.Element`; they are instances of
    the given class.

    INPUT:

    - ``typ`` -- a Python (new-style) class

    EXAMPLES::

        sage: from sage.sets.pythonclass import Set_PythonType
        sage: S = Set_PythonType(int); S
        Set of Python objects of class 'int'
        sage: int('1') in S
        True
        sage: Integer('1') in S
        False

        sage: Set_PythonType(2)
        Traceback (most recent call last):
        ...
        TypeError: must be initialized with a class, not 2
    """

    def __init__(self, typ):
        """
        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: Set_PythonType(float).category()
            Category of sets
        """
        if not isinstance(typ, type):
            raise TypeError(f"must be initialized with a class, not {typ!r}")
        super().__init__(category=Sets())
        self._type = <type>typ

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an instance of the class.

        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: S = Set_PythonType(complex)
            sage: S._element_constructor_(5)
            (5+0j)
            sage: S._element_constructor_(1, 5/2)
            (1+2.5j)
        """
        return self._type(*args, **kwds)

    def __reduce__(self):
        r"""
        Pickling support.

        TESTS::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: S = Set_PythonType(object)
            sage: loads(dumps(S))
            Set of Python objects of class 'object'
        """
        return type(self), (self._type,)

    def __call__(self, x):
        """
        Construct a new instance from ``x``. If ``x`` is already an
        instance of the correct class, directly return ``x`` itself.

        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: S = Set_PythonType(float)
            sage: S(5)
            5.0
            sage: S(9/3)
            3.0
            sage: S(1/3)
            0.333333333333333...
            sage: a = float(3); S(a) is a
            True
        """
        if isinstance(x, self._type):
            return x
        return self._type(x)

    def __hash__(self):
        """
        TESTS::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: S = Set_PythonType(int)
            sage: hash(S) == -hash(int)
            True
        """
        return -hash(self._type)

    def __richcmp__(self, other, int op):
        """
        Two Python class sets are considered the same if they contain
        the same class.

        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: S = Set_PythonType(int)
            sage: T = Set_PythonType(int)
            sage: U = type(S)(int)  # bypass caching
            sage: S is T
            True
            sage: S == T
            True
            sage: S is U
            False
            sage: S == U
            True
            sage: S == Set_PythonType(float)
            False
            sage: S == int
            False
        """
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        if self is other:
            return rich_to_bool(op, 0)
        if not isinstance(other, Set_PythonType_class):
            return rich_to_bool(op, 1)
        s = (<Set_PythonType_class>self)._type
        o = (<Set_PythonType_class>other)._type
        return rich_to_bool(op, s is not o)

    def __contains__(self, x):
        """
        Only things of the right class (or subclasses thereof) are
        considered to belong to the set.

        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: S = Set_PythonType(tuple)
            sage: (1,2,3) in S
            True
            sage: () in S
            True
            sage: [1,2] in S
            False
        """
        return isinstance(x, self._type)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: Set_PythonType(tuple)
            Set of Python objects of class 'tuple'
            sage: Set_PythonType(Integer)
            Set of Python objects of class 'Integer'
            sage: Set_PythonType(Parent)
            Set of Python objects of class 'Parent'
        """
        return f"Set of Python objects of class '{self._type.__name__}'"

    def object(self):
        """
        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: Set_PythonType(tuple).object()
            <... 'tuple'>
        """
        return self._type

    def cardinality(self):
        """
        EXAMPLES::

            sage: from sage.sets.pythonclass import Set_PythonType
            sage: S = Set_PythonType(bool)
            sage: S.cardinality()
            2
            sage: S = Set_PythonType(int)
            sage: S.cardinality()
            +Infinity
        """
        if self._type is bool:
            from sage.rings.integer import Integer
            return Integer(2)
        else:
            # Probably infinite
            import sage.rings.infinity
            return sage.rings.infinity.infinity
