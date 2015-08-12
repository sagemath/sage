"""
Metaclass for inheriting comparison functions

This module defines a metaclass :class:`InheritComparisonMetaclass` to
inherit comparison functions in Cython extension types. In Python 2,
the special methods ``__richcmp__``, ``__cmp__`` and ``__hash__`` are
only inherited as a whole: defining just 1 or 2 of these will prevent
the others from being inherited.

To solve this issue, you can use :class:`InheritComparisonMetaclass`
as a Cython "metaclass" (see :mod:`sage.misc.cython_metaclass` for the
general mechanism). If you do this for an extension type which defines
neither ``__richcmp__`` nor ``__cmp__``, then both these methods are
inherited from the base class (the MRO is not used).

In Sage, this is in particular used for
:class:`sage.structure.element.Element` to support comparisons using
the coercion framework.

None of this is relevant to Python classes, which inherit comparison
methods anyway.

AUTHOR:

- Jeroen Demeyer (2015-05-22): initial version, see :trac:`18329`
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport PyTypeObject
from sage.misc.classcall_metaclass cimport ClasscallMetaclass

cdef class InheritComparisonMetaclass(type):
    """
    If the type does not define ``__richcmp__`` nor ``__cmp__``,
    inherit both these methods from the base class. The difference with
    plain extension types is that comparison is inherited even if
    ``__hash__`` is defined.

    EXAMPLES::

        sage: cython('''
        ....: from sage.misc.inherit_comparison cimport InheritComparisonMetaclass
        ....:
        ....: cdef class Base(object):
        ....:     def __richcmp__(left, right, int op):
        ....:         print("Calling Base.__richcmp__")
        ....:         return left is right
        ....:
        ....: cdef class Derived(Base):
        ....:     def __hash__(self):
        ....:         return 1
        ....:
        ....: cdef class DerivedWithRichcmp(Base):
        ....:     def __getmetaclass__(_):
        ....:         from sage.misc.inherit_comparison import InheritComparisonMetaclass
        ....:         return InheritComparisonMetaclass
        ....:     def __hash__(self):
        ....:         return 1
        ....: ''')
        sage: a = Derived()
        sage: a == a
        True
        sage: b = DerivedWithRichcmp()
        sage: b == b
        Calling Base.__richcmp__
        True
    """
    def __init__(self, *args):
        cdef PyTypeObject* t = <PyTypeObject*>self
        cdef PyTypeObject* b = t.tp_base
        if b:
            if not t.tp_richcompare and not t.tp_compare:
                t.tp_richcompare = b.tp_richcompare
                t.tp_compare = b.tp_compare
        super(InheritComparisonMetaclass, self).__init__(*args)

class InheritComparisonClasscallMetaclass(InheritComparisonMetaclass, ClasscallMetaclass):
    pass
