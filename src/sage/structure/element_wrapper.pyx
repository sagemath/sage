# sage_setup: distribution = sagemath-objects
"""
Element Wrapper

Wrapping Sage or Python objects as Sage elements.

AUTHORS:

- Nicolas Thiery (2008-2010): Initial version
- Travis Scrimshaw (2013-05-04): Cythonized version
"""

# ****************************************************************************
#       Copyright (C) 2008-2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport Py_EQ, Py_NE, Py_LE, Py_GE, PyObject_RichCompare

from sage.structure.coerce cimport coercion_model
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.structure.unique_representation import UniqueRepresentation
from copy import copy


cdef class ElementWrapper(Element):
    r"""
    A class for wrapping Sage or Python objects as Sage elements.

    EXAMPLES::

        sage: from sage.structure.element_wrapper import DummyParent
        sage: parent = DummyParent("A parent")
        sage: o = ElementWrapper(parent, "bla"); o
        'bla'
        sage: isinstance(o, sage.structure.element.Element)
        True
        sage: o.parent()
        A parent
        sage: o.value
        'bla'

    Note that ``o`` is not *an instance of* ``str``, but rather
    *contains a* ``str``. Therefore, ``o`` does not inherit the string
    methods. On the other hand, it is provided with reasonable default
    implementations for equality testing, hashing, etc.

    The typical use case of ``ElementWrapper`` is for trivially
    constructing new element classes from pre-existing Sage or Python
    classes, with a containment relation. Here we construct the
    tropical monoid of integers endowed with ``min`` as
    multiplication. There, it is desirable *not* to inherit the
    ``factor`` method from ``Integer``::

        sage: class MinMonoid(Parent):
        ....:     def _repr_(self):
        ....:         return "The min monoid"
        ....:
        sage: M = MinMonoid()
        sage: class MinMonoidElement(ElementWrapper):
        ....:     wrapped_class = Integer
        ....:
        ....:     def __mul__(self, other):
        ....:         return MinMonoidElement(self.parent(), min(self.value, other.value))
        sage: x = MinMonoidElement(M, 5); x
        5
        sage: x.parent()
        The min monoid
        sage: x.value
        5
        sage: y = MinMonoidElement(M, 3)
        sage: x * y
        3

    This example was voluntarily kept to a bare minimum. See the
    examples in the categories (e.g. ``Semigroups().example()``) for
    several full featured applications.

    .. WARNING::

        Versions before :issue:`14519` had parent as the second argument and
        the value as the first.
    """
    def __init__(self, parent, value):
        """
        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: a = ElementWrapper(DummyParent("A parent"), 1)

        TESTS::

            sage: TestSuite(a).run(skip = "_test_category")

        .. NOTE::

            :class:`ElementWrapper` is not intended to be used directly,
            hence the failing category test.
        """
        Element.__init__(self, parent=parent)
        self.value = value

    def __getstate__(self):
        """
        Return a tuple describing the state of your object.

        This emulates :meth:`Element.__getstate__`, playing as if
        ``self.value`` was in the dictionary of ``self`` as it used to
        be before :issue:`14519`.

        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: P = DummyParent("A parent")
            sage: a = ElementWrapper(P, 1)
            sage: a.__getstate__()
            (A parent, {'value': 1})
            sage: class A(ElementWrapper):
            ....:     pass
            sage: a = A(P, 1)
            sage: a.x = 2
            sage: a.__getstate__() == (P, {'value': 1, 'x': 2})
            True
        """
        try:
            d = self.__dict__
        except AttributeError:
            d = {}
        else:
            d = d.copy()
        d['value'] = self.value
        return (self._parent, d)

    def __setstate__(self, state):
        r"""
        Initialize the state of the object from data saved in a pickle.

        This emulates :meth:`Element.__setstate__`, playing as if
        ``self.value`` was to be put in the dictionary of ``self`` as
        it used to be before :issue:`14519`.

        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: class A(ElementWrapper):
            ....:     pass
            sage: a = A(DummyParent("A parent"), 1)
            sage: a.__setstate__((DummyParent("Another parent"), {'value':0,'x':3}))
            sage: a.parent()
            Another parent
            sage: a.value
            0
            sage: a.x
            3

        TESTS::

            sage: a = A(DummyParent("A parent"), 1)
            sage: import __main__; __main__.A = A # Fake A being defined in a python module
            sage: a.x = 2
            sage: a == loads(dumps(a))
            True
            sage: a = ElementWrapper(DummyParent("A parent"), 1)
            sage: a == loads(dumps(a))
            True
        """
        # Make sure the first part of the state is the parent
        if not isinstance(state[0], Parent):
            state[0], state[1] = state[1], state[0]
        self._set_parent(state[0])
        d = state[1].copy()
        self.value = d.pop('value')
        if d:
            self.__dict__ = d

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: ElementWrapper(DummyParent("A parent"), 1)
            1
        """
        return repr(self.value)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: ElementWrapper(DummyParent("A parent"), 1)._latex_()
            1
            sage: ElementWrapper(DummyParent("A parent"), 3/5)._latex_()
            \frac{3}{5}
        """
        from sage.misc.latex import latex
        return latex(self.value)

    def _ascii_art_(self):
        r"""
        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: ElementWrapper(DummyParent("A parent"), 1)._ascii_art_()
            1
            sage: x = var('x')                                                          # needs sage.symbolic
            sage: ElementWrapper(DummyParent("A parent"), x^2 + x)._ascii_art_()        # needs sage.symbolic
             2
            x  + x
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(self.value)

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: ElementWrapper(DummyParent("A parent"), 1)._ascii_art_()
            1
            sage: x = var('x')                                                          # needs sage.symbolic
            sage: ElementWrapper(DummyParent("A parent"), x^2 + x)._unicode_art_()      # needs sage.symbolic
             2
            x  + x
        """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art(self.value)

    def __hash__(self):
        """
        Return the same hash as for the wrapped element.

        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: parent1 = DummyParent("A parent")
            sage: parent2 = DummyParent("Another parent")
            sage: hash(ElementWrapper(parent1, 1))
            1
            sage: hash(ElementWrapper(parent2, 1))
            1

        .. TODO::

            Should this take the parent and/or the class into account?
        """
        return hash(self.value)

    def __richcmp__(left, right, int op):
        """
        Return ``True`` if ``left`` compares with ``right`` based on ``op``.

        Default implementation of ``self == other``: two elements are
        equal if they have equal parents and equal values.

        Default implementation of ``self < other``: two elements are
        always incomparable.

        .. NOTE::

            Another option would be to not define ``__lt__``, but
            given the current implementation of SageObject, sorted(...)
            would break.

        TESTS:

        Check that elements of equal-but-not-identical parents compare
        properly (see :issue:`19488`)::

            sage: from sage.misc.test_nested_class import TestParent4
            sage: P = TestParent4()
            sage: Q = TestParent4()
            sage: P == Q
            True
            sage: P is Q
            False
            sage: x = P.an_element(); x
            '_an_element_'
            sage: y = Q.an_element(); y
            '_an_element_'
            sage: x == y
            True
        """
        if isinstance(right, ElementWrapper) and left.parent() == right.parent():
            return left._richcmp_(right, op)
        return coercion_model.richcmp(left, right, op)

    cpdef _richcmp_(left, right, int op):
        """
        Return ``True`` if ``left`` compares with ``right`` based on ``op``.

        TESTS:

        Testing equality::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: parent1 = DummyParent("A parent")
            sage: parent2 = DummyParent("Another parent")
            sage: parent1 == parent2
            False
            sage: l11 = ElementWrapper(parent1, 1)
            sage: l12 = ElementWrapper(parent1, 2)
            sage: l21 = ElementWrapper(parent2, 1)
            sage: l22 = ElementWrapper(parent2, 2)
            sage: l11 == l11
            True
            sage: l11 == l12
            False
            sage: l11 == l21
            False

        Testing inequality::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: parent1 = DummyParent("A parent")
            sage: parent2 = DummyParent("Another parent")
            sage: parent1 == parent2
            False
            sage: l11 = ElementWrapper(parent1, 1)
            sage: l12 = ElementWrapper(parent1, 2)
            sage: l21 = ElementWrapper(parent2, 1)
            sage: l22 = ElementWrapper(parent2, 2)
            sage: l11 != l11
            False
            sage: l11 != l12
            True
            sage: l11 != l21
            True

        Testing less than::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: parent = DummyParent("A parent")
            sage: x = ElementWrapper(parent, 1)
            sage: y = ElementWrapper(parent, 2)
            sage: x.__lt__(x), x.__lt__(y), y.__lt__(x)
            (False, False, False)
            sage: x < x, x < y, y < x
            (False, False, False)
            sage: sorted([x,y])
            [1, 2]
            sage: sorted([y,x])
            [2, 1]
        """
        cdef ElementWrapper self = left
        if op == Py_EQ or op == Py_LE or op == Py_GE:
            return self.value == (<ElementWrapper>right).value
        if op == Py_NE:
            return self.value != (<ElementWrapper>right).value
        return False

    cpdef bint _lt_by_value(self, other) noexcept:
        """
        Return whether ``self`` is strictly smaller than ``other``.

        With this implementation 'by value', they are always
        incomparable unless ``self`` and ``other`` have the same
        class, parent, and ``self.value < other.value``.

        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: class MyElement(ElementWrapper):
            ....:     __lt__ = ElementWrapper._lt_by_value
            ....:
            sage: parent1 = DummyParent("A parent")
            sage: parent2 = DummyParent("Another parent")
            sage: l11 = MyElement(parent1, 1)
            sage: l12 = MyElement(parent1, 2)
            sage: l21 = MyElement(parent2, 1)
            sage: l22 = MyElement(parent2, 2)
            sage: l11 < l11
            False
            sage: l11 < l12, l12 < l11   # values differ
            (True, False)
            sage: l11 < l21              # parents differ
            False
            sage: l11 < 1                # class differ
            False
        """
        return (self.__class__ is other.__class__
                and self._parent is other.parent()
                and self.value < (<ElementWrapper>other).value)

    def __copy__(self):
        """
        Copy ``self`` and in particular its ``value`` attribute.

        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: parent = DummyParent("A parent")
            sage: o1 = ElementWrapper(parent, [1]); o1
            [1]
            sage: o2 = copy(o1); o2
            [1]
            sage: o1 is o2, o1.value is o2.value
            (False, False)
            sage: o2.value[0] = 3; o2
            [3]
            sage: o1
            [1]
            sage: class bla(ElementWrapper): pass
            sage: o3 = bla(parent, [1])
            sage: o4 = copy(o3)
            sage: o3.value[0] = 3; o4
            [1]
            sage: o3.__class__
            <class '__main__.bla'>
            sage: o4.__class__
            <class '__main__.bla'>
        """
        # Note : copy(super()) does not work.
        res = super().__copy__()
        res.value = copy(self.value)
        return res


class DummyParent(UniqueRepresentation, Parent):
    """
    A class for creating dummy parents for testing :class:`ElementWrapper`
    """
    def __init__(self, name):
        """
        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: parent = DummyParent("A Parent")
            sage: skipped = ["_test_an_element", "_test_category",
            ....:            "_test_elements", "_test_elements_eq_reflexive",
            ....:            "_test_elements_eq_symmetric",
            ....:            "_test_elements_eq_transitive",
            ....:            "_test_elements_neq", "_test_some_elements"]
            sage: TestSuite(parent).run(skip=skipped)
        """
        self.name = name

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.structure.element_wrapper import DummyParent
            sage: DummyParent("A Parent") # indirect doctest
            A Parent
        """
        return self.name


class ElementWrapperTester(ElementWrapper):
    """
    Test class for the default :meth:`.__copy` method of subclasses of
    :class:`ElementWrapper`.

    TESTS::

        sage: from sage.structure.element_wrapper import ElementWrapperTester
        sage: x = ElementWrapperTester()
        sage: x.append(2); y = copy(x); y.append(42)
        sage: type(y)
        <class 'sage.structure.element_wrapper.ElementWrapperTester'>
        sage: x, y
        ([n=1, value=[2]], [n=2, value=[2, 42]])
        sage: x.append(21); x.append(7)
        sage: x, y
        ([n=3, value=[2, 21, 7]], [n=2, value=[2, 42]])
        sage: x.value, y.value
        ([2, 21, 7], [2, 42])
        sage: x.__dict__, y.__dict__
        ({'n': 3}, {'n': 2})
    """
    def __init__(self):
        """
        TESTS::

            sage: from sage.structure.element_wrapper import ElementWrapperTester
            sage: x = ElementWrapperTester(); x
            [n=0, value=[]]
        """
        from sage.categories.sets_cat import Sets
        super().__init__(Sets().example("facade"), [])
        self.n = 0

    def append(self, x):
        """
        TESTS::

            sage: from sage.structure.element_wrapper import ElementWrapperTester
            sage: x = ElementWrapperTester()
            sage: x.append(2); x
            [n=1, value=[2]]
        """
        self.n +=1
        self.value.append(x)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.structure.element_wrapper import ElementWrapperTester
            sage: x = ElementWrapperTester
            sage: x = ElementWrapperTester(); x
            [n=0, value=[]]
            sage: x.value = [2,32]; x # indirect doctest
            [n=0, value=[2, 32]]
        """
        return "[n=%s, value=%s]" % (self.n, self.value)


cdef class ElementWrapperCheckWrappedClass(ElementWrapper):
    """
    An :class:`element wrapper <ElementWrapper>` such that comparison
    operations are done against subclasses of ``wrapped_class``.
    """
    wrapped_class = object

    def __hash__(self):
        """
        Return the same hash as for the wrapped element.

        EXAMPLES::

            sage: A = cartesian_product([ZZ, ZZ])
            sage: e1 = A((6,9))
            sage: e2 = A((3,8))
            sage: e3 = A((6,9))
            sage: hash(e1) == hash(e2)
            False
            sage: hash(e1) == hash(e3)
            True
        """
        return hash(self.value)

    def __richcmp__(self, right, int op):
        """
        Return ``True`` if ``self`` compares with ``right`` based on ``op``.

        EXAMPLES::

            sage: A = cartesian_product([ZZ, ZZ])
            sage: elt = A((1,1))
            sage: (1, 1) == elt
            True
            sage: elt == (1, 1)
            True
            sage: A((1, 2)) > elt
            True
            sage: elts = [A((x,y)) for y in range(3) for x in range(2)]
            sage: elts
            [(0, 0), (1, 0), (0, 1), (1, 1), (0, 2), (1, 2)]
            sage: sorted(elts)
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

        ::

            sage: A = cartesian_product([ZZ, ZZ])
            sage: B = cartesian_product([GF(3), GF(5)])
            sage: A((3,5)) == B((0,0))
            True
        """
        if type(self) is type(right):
            # Both are instances of ElementWrapperCheckWrappedClass:
            # compare using wrapped element if the parents are the same
            other = <ElementWrapperCheckWrappedClass>right
            if self._parent is other._parent:
                return PyObject_RichCompare(self.value, other.value, op)
        elif not isinstance(right, Element):
            # Right is not an Element: compare using wrapped element
            return PyObject_RichCompare(self.value, right, op)
        elif self._parent is (<Element>right)._parent:
            # Different types but same parent? This should not happen
            raise TypeError(f"cannot compare {type(self).__name__} with {type(right).__name__} if parents are equal")

        # Different parents => use coercion model
        return coercion_model.richcmp(self, right, op)
