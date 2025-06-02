# cython: binding=True
# sage_setup: distribution = sagemath-objects
r"""
Cython-like rich comparisons in Python

With "rich comparisons", we mean the Python 3 comparisons which are
usually implemented in Python using methods like ``__eq__`` and
``__lt__``. Internally in Python, there is only one rich comparison
slot ``tp_richcompare``. The actual operator is passed as an integer
constant (defined in this module as
``op_LT``, ``op_LE``, ``op_EQ``, ``op_NE``, ``op_GT``, ``op_GE``).

Cython exposes rich comparisons in ``cdef`` classes as the
``__richcmp__`` special method. The Sage coercion model also supports
rich comparisons this way: for two instances ``x`` and ``y``
of :class:`~sage.structure.element.Element`, ``x._richcmp_(y, op)``
is called when the user does something like ``x <= y``
(possibly after coercion if ``x`` and ``y`` have different parents).

Various helper functions exist to make it easier to implement rich
comparison: the most important one is the :func:`richcmp` function.
This is analogous to the Python 2 function ``cmp()`` but implements
rich comparison, with the comparison operator (e.g. ``op_GE``) as
third argument. There is also :func:`richcmp_not_equal` which is like
:func:`richcmp` but it is optimized assuming that the compared objects
are not equal.

The functions :func:`rich_to_bool` and :func:`rich_to_bool_sgn` can be
used to convert results of ``cmp()`` (i.e. -1, 0 or 1) to a boolean
``True``/``False`` for rich comparisons.

AUTHORS:

- Jeroen Demeyer
"""

# ****************************************************************************
#       Copyright (C) 2017-2018 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport PyTypeObject
from sage.cpython.wrapperdescr cimport get_slotdef, wrapperbase, PyDescr_NewWrapper

cdef extern from *:
    void PyType_Modified(PyTypeObject* cls)


op_LT = Py_LT   # operator <
op_LE = Py_LE   # operator <=
op_EQ = Py_EQ   # operator ==
op_NE = Py_NE   # operator !=
op_GT = Py_GT   # operator >
op_GE = Py_GE   # operator >=


# Slotdefs for richcmp methods (the "bytes" below is arbitrary,
# any type which implements these methods would work)
cdef wrapperbase* richcmp_slotdef[6]
richcmp_slotdef[Py_EQ] = get_slotdef(bytes.__eq__)
richcmp_slotdef[Py_NE] = get_slotdef(bytes.__ne__)
richcmp_slotdef[Py_LT] = get_slotdef(bytes.__lt__)
richcmp_slotdef[Py_GT] = get_slotdef(bytes.__gt__)
richcmp_slotdef[Py_LE] = get_slotdef(bytes.__le__)
richcmp_slotdef[Py_GE] = get_slotdef(bytes.__ge__)


cpdef richcmp_item(x, y, int op):
    """
    This function is meant to implement lexicographic rich comparison
    of sequences (lists, vectors, polynomials, ...).
    The inputs ``x`` and ``y`` are corresponding items of such lists
    which should compared.

    INPUT:

    - ``x``, ``y`` -- arbitrary Python objects; typically, these are
      ``X[i]`` and ``Y[i]`` for sequences ``X`` and ``Y``

    - ``op`` -- comparison operator (one of ``op_LT``, ``op_LE``,
      ``op_EQ``, ``op_NE``, ``op_GT``, ``op_GE``)

    OUTPUT: assuming that ``x = X[i]`` and ``y = Y[i]``:

    - if the comparison ``X {op} Y`` (where ``op`` is the given
      operation) could not be decided yet (i.e. we should compare the
      next items in the list): return ``NotImplemented``

    - otherwise, if the comparison ``X {op} Y`` could be decided:
      return ``x {op} y``, which should then also be the result for
      ``X {op} Y``.

    .. NOTE::

        Since ``x {op} y`` cannot return ``NotImplemented``, the two
        cases above are mutually exclusive.

    The semantics of the comparison is different from Python lists or
    tuples in the case that the order is not total. Assume that ``A``
    and ``B`` are lists whose rich comparison is implemented using
    ``richcmp_item`` (as in EXAMPLES below). Then

    - ``A == B`` iff ``A[i] == B[i]`` for all indices `i`.

    - ``A != B`` iff ``A[i] != B[i]`` for some index `i`.

    - ``A < B`` iff ``A[i] < B[i]`` for some index `i` and
      for all `j < i`, ``A[j] <= B[j]``.

    - ``A <= B`` iff ``A < B`` or ``A[i] <= B[i]`` for all `i`.

    - ``A > B`` iff ``A[i] > B[i]`` for some index `i` and
      for all `j < i`, ``A[j] >= B[j]``.

    - ``A >= B`` iff ``A > B`` or ``A[i] >= B[i]`` for all `i`.

    See below for a detailed description of the exact semantics of
    ``richcmp_item`` in general.

    EXAMPLES::

        sage: from sage.structure.richcmp import *
        sage: @richcmp_method
        ....: class Listcmp(list):
        ....:     def __richcmp__(self, other, op):
        ....:         for i in range(len(self)):  # Assume equal lengths
        ....:             res = richcmp_item(self[i], other[i], op)
        ....:             if res is not NotImplemented:
        ....:                 return res
        ....:         return rich_to_bool(op, 0)  # Consider the lists to be equal
        sage: a = Listcmp([0, 1, 3])
        sage: b = Listcmp([0, 2, 1])
        sage: a == a
        True
        sage: a != a
        False
        sage: a < a
        False
        sage: a <= a
        True
        sage: a > a
        False
        sage: a >= a
        True
        sage: a == b, b == a
        (False, False)
        sage: a != b, b != a
        (True, True)
        sage: a < b, b > a
        (True, True)
        sage: a <= b, b >= a
        (True, True)
        sage: a > b, b < a
        (False, False)
        sage: a >= b, b <= a
        (False, False)

    The above tests used a list of integers, where the result of
    comparisons are the same as for Python lists.

    If we want to see the difference, we need more general entries in
    the list. The comparison rules are made to be consistent with
    setwise operations. If `A` and `B` are sets, we define ``A {op} B``
    to be true if ``a {op} B`` is true for every `a` in `A` and
    `b` in `B`. Interval comparisons are a special case of this. For
    lists of non-empty(!) sets, we want that ``[A1, A2] {op} [B1, B2]``
    is true if and only if ``[a1, a2] {op} [b1, b2]`` is true for all
    elements. We verify this::

        sage: @richcmp_method
        ....: class Setcmp(tuple):
        ....:     def __richcmp__(self, other, op):
        ....:         return all(richcmp(x, y, op) for x in self for y in other)
        sage: sym = {op_EQ: "==", op_NE: "!=", op_LT: "<", op_GT: ">", op_LE: "<=", op_GE: ">="}
        sage: for A1 in Set(range(4)).subsets():  # long time
        ....:     if not A1: continue
        ....:     for B1 in Set(range(4)).subsets():
        ....:         if not B1: continue
        ....:         for A2 in Set(range(4)).subsets():
        ....:             if not A2: continue
        ....:             for B2 in Set(range(3)).subsets():
        ....:                 if not B2: continue
        ....:                 L1 = Listcmp([Setcmp(A1), Setcmp(A2)])
        ....:                 L2 = Listcmp([Setcmp(B1), Setcmp(B2)])
        ....:                 for op in range(6):
        ....:                     reslist = richcmp(L1, L2, op)
        ....:                     reselt = all(richcmp([a1, a2], [b1, b2], op) for a1 in A1 for a2 in A2 for b1 in B1 for b2 in B2)
        ....:                     assert reslist is reselt

    EXACT SEMANTICS:

    Above, we only described how ``richcmp_item`` behaves when it is
    used to compare sequences. Here, we specify the exact semantics.
    First of all, recall that the result of ``richcmp_item(x, y, op)``
    is either ``NotImplemented`` or ``x {op} y``.

    - if ``op`` is ``==``: return ``NotImplemented`` if ``x == y``.
      If ``x == y`` is false, then return ``x == y``.

    - if ``op`` is ``!=``: return ``NotImplemented`` if not ``x != y``.
      If ``x != y`` is true, then return ``x != y``.

    - if ``op`` is ``<``: return ``NotImplemented`` if ``x == y``.
      If ``x < y`` or not ``x <= y``, return ``x < y``.
      Otherwise (if both ``x == y`` and ``x < y`` are false but
      ``x <= y`` is true), return ``NotImplemented``.

    - if ``op`` is ``<=``: return ``NotImplemented`` if ``x == y``.
      If ``x < y`` or not ``x <= y``, return ``x <= y``.
      Otherwise (if both ``x == y`` and ``x < y`` are false but
      ``x <= y`` is true), return ``NotImplemented``.

    - the ``>`` and ``>=`` operators are analogous to ``<`` and ``<=``.
    """
    if op == Py_NE:
        res = (x != y)
        if not res:
            return NotImplemented
        return res  # (x != y)  --> True

    # If x and y are equal, we cannot decide
    res = (x == y)
    if res:
        return NotImplemented

    if op == Py_EQ:
        return res  # not (x == y)  --> False

    # At this point, {op} is < or <= or > or >=. In the comments below,
    # we always refer to < and <= but > and >= are obviously analogous.

    # Compute x {op} y and convert to boolean (0 or 1)
    res = PyObject_RichCompare(x, y, op)
    cdef bint bres = res

    # true (1) for <= and >= and false (0) for < and >
    cdef bint op_is_not_strict = op & 1

    if bres != op_is_not_strict:
        # If we are asked to compute (x < y) and (x < y) is true,
        # return (x < y) which is true.
        # If we are asked to compute (x <= y) and (x <= y) is false,
        # return (x <= y) which is false.
        return res

    # Finally, check the inequality with the other strictness
    # (< becomes <= and vice versa; this corresponds to replacing op by
    # op ^ 1). This check is redundant in the typical case that (x <= y)
    # is equivalent to (x < y or x == y): since (x == y) returned false,
    # we expect that this PyObject_RichCompare() call returns the same
    # boolean result as the previous one.
    cdef bint xres = PyObject_RichCompare(x, y, op ^ 1)
    if xres == bres:  # As expected
        return res

    # OK, we are in a special case now. We checked that (x == y) is
    # false, that (x < y) is false but (x <= y) is true.
    # Since we want to give more importance to the < and <= results than
    # the == result, we treat this case as equality. Therefore, we
    # cannot decide.
    return NotImplemented


cdef slot_tp_richcompare(self, other, int op):
    """
    Function to put in the ``tp_richcompare`` slot.
    """
    return self.__richcmp__(other, op)


def richcmp_method(cls):
    """
    Class decorator to implement rich comparison using the special
    method ``__richcmp__`` (analogous to Cython) instead of the 6
    methods ``__eq__`` and friends.

    This changes the class in-place and returns the given class.

    EXAMPLES::

        sage: from sage.structure.richcmp import *
        sage: sym = {op_EQ: "==", op_NE: "!=", op_LT: "<", op_GT: ">", op_LE: "<=", op_GE: ">="}
        sage: @richcmp_method
        ....: class A(str):
        ....:     def __richcmp__(self, other, op):
        ....:         print("%s %s %s" % (self, sym[op], other))
        sage: A("left") < A("right")
        left < right
        sage: object() <= A("right")
        right >= <object object at ...>

    We can call this comparison with the usual Python special methods::

        sage: x = A("left"); y = A("right")
        sage: x.__eq__(y)
        left == right
        sage: A.__eq__(x, y)
        left == right

    Everything still works in subclasses::

        sage: class B(A):
        ....:     pass
        sage: x = B("left"); y = B("right")
        sage: x != y
        left != right
        sage: x.__ne__(y)
        left != right
        sage: B.__ne__(x, y)
        left != right

    We can override ``__richcmp__`` with standard Python rich
    comparison methods and conversely::

        sage: class C(A):
        ....:     def __ne__(self, other):
        ....:         return False
        sage: C("left") != C("right")
        False
        sage: C("left") == C("right")  # Calls __eq__ from class A
        left == right

        sage: class Base():
        ....:     def __eq__(self, other):
        ....:         return False
        sage: @richcmp_method
        ....: class Derived(Base):
        ....:     def __richcmp__(self, other, op):
        ....:         return True
        sage: Derived() == Derived()
        True

    TESTS::

        sage: richcmp_method(None)
        Traceback (most recent call last):
        ...
        TypeError: None is not a class

        sage: @richcmp_method
        ....: class X():
        ....:     def __eq__(self, other):
        ....:         pass
        ....:     def __richcmp__(self, other, op):
        ....:         pass
        Traceback (most recent call last):
        ...
        TypeError: class <class '__main__.X'> defines __eq__ which cannot be combined with @richcmp_method
    """
    if not isinstance(cls, type):
        raise TypeError(f"{cls!r} is not a class")
    cdef PyTypeObject* tp = <PyTypeObject*>cls
    tp.tp_richcompare = slot_tp_richcompare

    # Install slot wrappers in the class dict
    cdef dict D = <dict>(tp.tp_dict)
    cdef wrapperbase* slotdef
    for slotdef in richcmp_slotdef:
        name = <object>slotdef.name_strobj
        if name in D:
            raise TypeError("class %r defines %s which cannot be combined with @richcmp_method" % (cls, name))
        D[name] = PyDescr_NewWrapper(tp, slotdef, <void*>slot_tp_richcompare)

    PyType_Modified(tp)

    return cls


def richcmp_by_eq_and_lt(eq_attr, lt_attr):
    r"""
    Create a rich comparison method for a partial order, where the
    order is specified by methods called ``eq_attr`` and ``lt_attr``.

    INPUT when creating the method:

    - ``eq_attr`` -- attribute name for equality comparison

    - ``lt_attr`` -- attribute name for less-than comparison

    INPUT when calling the method:

    - ``self`` -- objects having methods ``eq_attr`` and ``lt_attr``

    - ``other`` -- arbitrary object. If it does have ``eq_attr`` and
      ``lt_attr`` methods, these are used for the comparison. Otherwise,
      the comparison is undefined.

    - ``op`` -- a rich comparison operation (e.g. ``op_EQ``)

    .. NOTE::

        For efficiency, identical objects (when ``self is other``)
        always compare equal.

    .. NOTE::

        The order is partial, so ``x <= y`` is implemented as
        ``x == y or x < y``. It is not required that this is the
        negation of ``y < x``.

    .. NOTE::

        This function is intended to be used as a method ``_richcmp_``
        in a class derived from :class:`sage.structure.element.Element`
        or a method ``__richcmp__`` in a class using
        :func:`richcmp_method`.

    EXAMPLES::

        sage: from sage.structure.richcmp import richcmp_by_eq_and_lt
        sage: from sage.structure.element import Element

        sage: class C(Element):
        ....:     def __init__(self, a, b):
        ....:         super().__init__(ZZ)
        ....:         self.a = a
        ....:         self.b = b
        ....:     _richcmp_ = richcmp_by_eq_and_lt("eq", "lt")
        ....:     def eq(self, other):
        ....:         return self.a == other.a and self.b == other.b
        ....:     def lt(self, other):
        ....:         return self.a < other.a and self.b < other.b

        sage: x = C(1,2); y = C(2,1); z = C(3,3)

        sage: x == x, x <= x, x == C(1,2), x <= C(1,2)  # indirect doctest
        (True, True, True, True)
        sage: y == z, y != z
        (False, True)

        sage: x < y, y < x, x > y, y > x, x <= y, y <= x, x >= y, y >= x
        (False, False, False, False, False, False, False, False)
        sage: y < z, z < y, y > z, z > y, y <= z, z <= y, y >= z, z >= y
        (True, False, False, True, True, False, False, True)
        sage: z < x, x < z, z > x, x > z, z <= x, x <= z, z >= x, x >= z
        (False, True, True, False, False, True, True, False)

    A simple example using ``richcmp_method``::

        sage: from sage.structure.richcmp import richcmp_method, richcmp_by_eq_and_lt
        sage: @richcmp_method
        ....: class C():
        ....:     __richcmp__ = richcmp_by_eq_and_lt("_eq", "_lt")
        ....:     def _eq(self, other):
        ....:         return True
        ....:     def _lt(self, other):
        ....:         return True
        sage: a = C(); b = C()
        sage: a == b
        True
        sage: a > b  # Calls b._lt(a)
        True
        sage: class X(): pass
        sage: x = X()
        sage: a == x  # Does not call a._eq(x) because x does not have _eq
        False
    """
    def richcmp(self, other, int op):
        if self is other:
            return rich_to_bool(op, 0)

        cdef bint equal_types = (type(self) is type(other))
        # Check whether "other" is of the right type
        if not equal_types:
            try:
                other_eq = getattr(other, eq_attr)
                other_lt = getattr(other, lt_attr)
            except AttributeError:
                return NotImplemented

        # Check equality first if needed
        if op != Py_LT and op != Py_GT:
            if equal_types:
                other_eq = getattr(other, eq_attr)
            if other_eq(self):
                return rich_to_bool(op, 0)
            if op == Py_EQ:
                return False
            if op == Py_NE:
                return True

        if op == Py_LT or op == Py_LE:
            self_lt = getattr(self, lt_attr)
            return self_lt(other)
        else:
            if equal_types:
                other_lt = getattr(other, lt_attr)
            return other_lt(self)

    return richcmp
