"""
GAP element wrapper

This document describes the individual wrappers for various GAP
elements. For general information about GAP, you should read the
:mod:`~sage.libs.gap.libgap` module documentation.
"""

# hi

# ****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport Py_EQ, Py_NE, Py_LE, Py_GE, Py_LT, Py_GT
from cysignals.signals cimport sig_on, sig_off

from sage.libs.gap.gap_includes cimport *
from sage.libs.gap.libgap import libgap
from sage.libs.gap.util cimport *
from sage.libs.gap.util import GAPError
from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.pylong cimport mpz_get_pylong
from sage.cpython.string cimport str_to_bytes, char_to_str
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF

from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement
from sage.combinat.permutation import Permutation
from sage.structure.coerce cimport coercion_model as cm

############################################################################
### helper functions to construct lists and records ########################
############################################################################

cdef Obj make_gap_list(sage_list) except NULL:
    """
    Convert Sage lists into Gap lists

    INPUT:

    - ``a`` -- list of :class:`GapElement`

    OUTPUT: list of the elements in ``a`` as a Gap ``Obj``
    """
    cdef Obj l
    cdef GapElement elem
    cdef int i
    try:
        GAP_Enter()
        l = GAP_NewPlist(0)

        for i, x in enumerate(sage_list):
            if not isinstance(x, GapElement):
                elem = <GapElement>libgap(x)
            else:
                elem = <GapElement>x

            GAP_AssList(l, i + 1, elem.value)
        return l
    finally:
        GAP_Leave()


cdef Obj make_gap_matrix(sage_list, gap_ring) except NULL:
    """
    Convert Sage lists into Gap matrices

    INPUT:

    - ``sage_list`` -- list of :class:`GapElement`

    - ``gap_ring`` -- the base ring

    If ``gap_ring`` is ``None``, nothing is made to make sure
    that all coefficients live in the same Gap ring. The resulting Gap list
    may not be recognized as a matrix by Gap.

    OUTPUT:

    The list of the elements in ``sage_list`` as a Gap ``Obj``.
    """
    cdef Obj l
    cdef GapElement elem
    cdef GapElement one
    cdef int i
    if gap_ring is not None:
        one = <GapElement>gap_ring.One()
    else:
        one = <GapElement>libgap(1)

    try:
        GAP_Enter()
        l = GAP_NewPlist(0)

        for i, x in enumerate(sage_list):
            if not isinstance(x, GapElement):
                elem = <GapElement>libgap(x)
                elem = elem * one
            else:
                elem = <GapElement>x

            GAP_AssList(l, i + 1, elem.value)
        return l
    finally:
        GAP_Leave()


cdef char *capture_stdout(Obj func, Obj obj) noexcept:
    """
    Call a single-argument GAP function ``func`` with the argument ``obj``
    and return the stdout from that function call.

    This can be used to capture the output of GAP functions that are used to
    print objects such as ``Print()`` and ``ViewObj()``.
    """
    cdef Obj s, stream, output_text_string
    cdef Obj l
    # The only way to get a string representation of an object that is truly
    # consistent with how it would be represented at the GAP REPL is to call
    # ViewObj on it.  Unfortunately, ViewObj *prints* to the output stream,
    # and there is no equivalent that simply returns the string that would be
    # printed.  The closest approximation would be DisplayString, but this
    # bypasses any type-specific overrides for ViewObj so for many objects
    # that does not give consistent results.
    # TODO: This is probably needlessly slow, but we might need better
    # support from GAP to improve this...
    try:
        GAP_Enter()
        s = NEW_STRING(0)
        output_text_string = GAP_ValueGlobalVariable("OutputTextString")
        stream = GAP_CallFunc2Args(output_text_string, s, GAP_True)

        l = GAP_NewPlist(1)
        GAP_AssList(l, 1, obj)

        CALL_WITH_STREAM = GAP_ValueGlobalVariable("CALL_WITH_STREAM")
        GAP_CallFunc3Args(CALL_WITH_STREAM, stream, func, l)
        return GAP_CSTR_STRING(s)
    finally:
        GAP_Leave()


cdef char *gap_element_repr(Obj obj) noexcept:
    """
    Implement ``repr()`` of ``GapElement``s using the ``ViewObj()`` function,
    which is by default closest to what you get when displaying an object in
    GAP on the command-line (i.e. when evaluating an expression that returns
    that object.
    """
    cdef Obj func
    try:
        GAP_Enter()
        func = GAP_ValueGlobalVariable("ViewObj")
        return capture_stdout(func, obj)
    finally:
        GAP_Leave()


cdef char *gap_element_str(Obj obj) noexcept:
    """
    Implement ``str()`` of ``GapElement``s using the ``Print()`` function.

    This mirrors somewhat how Python uses ``str()`` on an object when passing
    it to the ``print()`` function.  This is also how the GAP pexpect interface
    has traditionally repr'd objects; for the libgap interface we take a
    slightly different approach more closely mirroring Python's str/repr
    difference (though this does not map perfectly onto GAP).
    """
    cdef Obj func
    try:
        GAP_Enter()
        func = GAP_ValueGlobalVariable("Print")
        return capture_stdout(func, obj)
    finally:
        GAP_Leave()


cdef Obj make_gap_record(sage_dict) except NULL:
    """
    Convert Sage dictionary into Gap record

    INPUT:

    - ``a`` -- list of :class:`GapElement`

    OUTPUT: the list of the elements in ``a`` as a Gap ``Obj``

    TESTS::

        sage: libgap({'a': 1, 'b':123})   # indirect doctest
        rec( a := 1, b := 123 )
    """
    data = [ (str(key), libgap(value)) for key, value in sage_dict.iteritems() ]

    cdef Obj rec
    cdef GapElement val
    cdef UInt rnam

    try:
        GAP_Enter()
        rec = GAP_NewPrecord(len(data))
        for d in data:
            key, val = d
            rnam = RNamName(str_to_bytes(key))
            AssPRec(rec, rnam, val.value)
        return rec
    finally:
        GAP_Leave()


cdef Obj make_gap_integer(sage_int) except NULL:
    """
    Convert Sage integer into Gap integer

    INPUT:

    - ``sage_int`` -- Sage integer

    OUTPUT: the integer as a GAP ``Obj``

    TESTS::

        sage: libgap(1)   # indirect doctest
        1
    """
    return GAP_NewObjIntFromInt(<int>sage_int)


cdef Obj make_gap_string(sage_string) except NULL:
    """
    Convert a Python string to a Gap string

    INPUT:

    - ``sage_string`` -- a Python str

    OUTPUT: the string as a GAP ``Obj``

    TESTS::

        sage: libgap('string')   # indirect doctest
        "string"
    """
    cdef Obj result
    try:
        GAP_Enter()
        b = str_to_bytes(sage_string)
        result = GAP_MakeStringWithLen(b, len(b))
        return result
    finally:
        GAP_Leave()


############################################################################
### generic construction of GapElements ####################################
############################################################################

cdef GapElement make_any_gap_element(parent, Obj obj):
    """
    Return the GAP element wrapper of ``obj``.

    The most suitable subclass of GapElement is determined
    automatically. Use this function to wrap GAP objects unless you
    know exactly which type it is (then you can use the specialized
    ``make_GapElement_...``)

    TESTS::

        sage: x = libgap.eval("'c'");  x
        "c"
        sage: type(x)
        <class 'sage.libs.gap.element.GapElement_String'>

        sage: libgap.eval("['a', 'b', 'c']")   # gap strings are also lists of chars
        "abc"
        sage: t = libgap.UnorderedTuples('abc', 2);  t
        [ "aa", "ab", "ac", "bb", "bc", "cc" ]
        sage: t[1]
        "ab"
        sage: t[1].sage()
        'ab'
        sage: t.sage()
        ['aa', 'ab', 'ac', 'bb', 'bc', 'cc']

    Check that :issue:`18158` is fixed::

        sage: S = SymmetricGroup(5)
        sage: irr = libgap.Irr(S)[3]
        sage: irr[0]
        6
        sage: irr[1]
        0
    """
    cdef int num

    try:
        GAP_Enter()
        if obj is NULL:
            return make_GapElement(parent, obj)
        num = TNUM_OBJ(obj)
        if GAP_IsInt(obj):
            return make_GapElement_Integer(parent, obj)
        elif GAP_IsMacFloat(obj):
            return make_GapElement_Float(parent, obj)
        elif num == T_CYC:
            return make_GapElement_Cyclotomic(parent, obj)
        elif num == T_FFE:
            return make_GapElement_FiniteField(parent, obj)
        elif num == T_RAT:
            return make_GapElement_Rational(parent, obj)
        elif num == T_BOOL:
            return make_GapElement_Boolean(parent, obj)
        elif num == T_FUNCTION:
            return make_GapElement_Function(parent, obj)
        elif num == T_PERM2 or num == T_PERM4:
            return make_GapElement_Permutation(parent, obj)
        elif GAP_IsRecord(obj):
            return make_GapElement_Record(parent, obj)
        elif GAP_IsList(obj) and GAP_LenList(obj) == 0:
            # Empty lists are lists and not strings in Python
            return make_GapElement_List(parent, obj)
        elif IsStringConv(obj):
            # GAP strings are lists, too. Make sure this comes before non-empty make_GapElement_List
            return make_GapElement_String(parent, obj)
        elif GAP_IsList(obj):
            return make_GapElement_List(parent, obj)
        elif GAP_ValueOfChar(obj) != -1:
            ch = make_GapElement(parent, obj).IntChar().sage()
            return make_GapElement_String(parent, make_gap_string(chr(ch)))
        result = make_GapElement(parent, obj)
        if num == T_POSOBJ:
            if result.IsZmodnZObj():
                return make_GapElement_IntegerMod(parent, obj)
        if num == T_COMOBJ:
            if result.IsRing():
                return make_GapElement_Ring(parent, obj)
        return result
    finally:
        GAP_Leave()


############################################################################
### GapElement #############################################################
############################################################################

cdef GapElement make_GapElement(parent, Obj obj):
    r"""
    Turn a Gap C object (of type ``Obj``) into a Cython ``GapElement``.

    INPUT:

    - ``parent`` -- the parent of the new :class:`GapElement`

    - ``obj`` -- a GAP object

    OUTPUT:

    A :class:`GapElement_Function` instance, or one of its derived
    classes if it is a better fit for the GAP object.

    EXAMPLES::

        sage: libgap(0)
        0
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Integer'>

        sage: libgap.eval('')
        sage: libgap(None)
        Traceback (most recent call last):
        ...
        AttributeError: 'NoneType' object has no attribute '_gap_init_'...
    """
    cdef GapElement r = GapElement.__new__(GapElement)
    r._initialize(parent, obj)
    return r


cpdef _from_sage(elem):
    """
    Currently just used for unpickling; equivalent to calling ``libgap(elem)``
    to convert a Sage object to a `GapElement` where possible.
    """
    if isinstance(elem, str):
        return libgap.eval(elem)

    return libgap(elem)


cdef class GapElement(RingElement):
    r"""
    Wrapper for all Gap objects.

    .. NOTE::

        In order to create ``GapElements`` you should use the
        ``libgap`` instance (the parent of all Gap elements) to
        convert things into ``GapElement``. You must not create
        ``GapElement`` instances manually.

    EXAMPLES::

        sage: libgap(0)
        0

    If Gap finds an error while evaluating, a :exc:`GAPError`
    exception is raised::

        sage: libgap.eval('1/0')
        Traceback (most recent call last):
        ...
        GAPError: Error, Rational operations: <divisor> must not be zero

    Also, a :exc:`GAPError` is raised if the input is not a simple expression::

        sage: libgap.eval('1; 2; 3')
        Traceback (most recent call last):
        ...
        GAPError: can only evaluate a single statement
    """

    def __cinit__(self):
        """
        The Cython constructor.

        EXAMPLES::

            sage: libgap.eval('1')
            1
        """
        self.value = NULL
        self._compare_by_id = False

    def __init__(self):
        """
        The ``GapElement`` constructor.

        Users must use the ``libgap`` instance to construct instances
        of :class:`GapElement`. Cython programmers must use
        :func:`make_GapElement` factory function.

        TESTS::

            sage: from sage.libs.gap.element import GapElement
            sage: GapElement()
            Traceback (most recent call last):
            ...
            TypeError: this class cannot be instantiated from Python
        """
        raise TypeError('this class cannot be instantiated from Python')

    cdef _initialize(self, parent, Obj obj):
        r"""
        Initialize the GapElement.

        This Cython method is called from :func:`make_GapElement` to
        initialize the newly-constructed object. You must never call
        it manually.

        TESTS::

            sage: n_before = libgap.count_GAP_objects()
            sage: a = libgap.eval('123')
            sage: b = libgap.eval('456')
            sage: c = libgap.eval('CyclicGroup(3)')
            sage: d = libgap.eval('"a string"')
            sage: libgap.collect()
            sage: del c
            sage: libgap.collect()
            sage: n_after = libgap.count_GAP_objects()
            sage: n_after - n_before
            3
        """
        assert self.value is NULL
        self._parent = parent
        self.value = obj
        if obj is NULL:
            return
        reference_obj(obj)

    def __dealloc__(self):
        r"""
        The Cython destructor.

        TESTS::

            sage: pre_refcount = libgap.count_GAP_objects()
            sage: def f():
            ....:     local_variable = libgap.eval('"This is a new string"')
            sage: f()
            sage: f()
            sage: f()
            sage: post_refcount = libgap.count_GAP_objects()
            sage: post_refcount - pre_refcount
            0
        """
        if self.value is NULL:
            return
        dereference_obj(self.value)

    def __copy__(self):
        r"""
        TESTS::

            sage: a = libgap(1)
            sage: a.__copy__() is a
            True

            sage: a = libgap(1/3)
            sage: a.__copy__() is a
            True

            sage: a = libgap([1,2])
            sage: b = a.__copy__()
            sage: a is b
            False
            sage: a[0] = 3
            sage: a
            [ 3, 2 ]
            sage: b
            [ 1, 2 ]

            sage: a = libgap([[0,1],[2,3,4]])
            sage: b = a.__copy__()
            sage: b[0][1] = -2
            sage: b
            [ [ 0, -2 ], [ 2, 3, 4 ] ]
            sage: a
            [ [ 0, -2 ], [ 2, 3, 4 ] ]
        """
        if IS_MUTABLE_OBJ(self.value):
            return make_any_gap_element(self.parent(), SHALLOW_COPY_OBJ(self.value))
        else:
            return self

    cpdef GapElement deepcopy(self, bint mut):
        r"""
        Return a deepcopy of this Gap object.

        Note that this is the same thing as calling ``StructuralCopy`` but much
        faster.

        INPUT:

        - ``mut`` -- boolean; whether to return a mutable copy

        EXAMPLES::

            sage: a = libgap([[0,1],[2,3]])
            sage: b = a.deepcopy(1)
            sage: b[0,0] = 5
            sage: a
            [ [ 0, 1 ], [ 2, 3 ] ]
            sage: b
            [ [ 5, 1 ], [ 2, 3 ] ]

            sage: l = libgap([0,1])
            sage: l.deepcopy(0).IsMutable()
            false
            sage: l.deepcopy(1).IsMutable()
            true
        """
        if IS_MUTABLE_OBJ(self.value):
            return make_any_gap_element(self.parent(), CopyObj(self.value, mut))
        else:
            return self

    def __deepcopy__(self, memo):
        r"""
        TESTS::

            sage: a = libgap([[0,1],[2]])
            sage: b = deepcopy(a)
            sage: a[0,0] = -1
            sage: a
            [ [ -1, 1 ], [ 2 ] ]
            sage: b
            [ [ 0, 1 ], [ 2 ] ]
        """
        return self.deepcopy(0)

    def __reduce__(self):
        """
        Attempt to pickle GAP elements from libgap.

        This is inspired in part by
        ``sage.interfaces.interface.Interface._reduce``, though for a fallback
        we use ``str(self)`` instead of ``repr(self)``, since the former is
        equivalent in the libgap interface to the latter in the pexpect
        interface.

        TESTS:

        This workaround was motivated in particular by this example from the
        permutation groups implementation::

            sage: CC = libgap.eval('ConjugacyClass(SymmetricGroup([ 1 .. 5 ]), (1,2)(3,4))')
            sage: CC.sage()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot construct equivalent Sage object
            sage: libgap.eval(str(CC))
            (1,2)(3,4)^G
            sage: loads(dumps(CC))
            (1,2)(3,4)^G
        """

        if self.is_string():
            elem = repr(self.sage())
        try:
            elem = self.sage()
        except NotImplementedError:
            elem = str(self)

        return (_from_sage, (elem,))

    def __contains__(self, other):
        r"""
        TESTS::

            sage: libgap(1) in libgap.eval('Integers')
            True
            sage: 1 in libgap.eval('Integers')
            True

            sage: 3 in libgap([1,5,3,2])
            True
            sage: -5 in libgap([1,5,3,2])
            False

            sage: libgap.eval('Integers') in libgap(1)
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found! Error, no 1st choice method found for `in' on 2 arguments
        """
        GAP_IN = libgap.eval(r'\in')
        return GAP_IN(other, self).sage()

    cpdef _type_number(self):
        """
        Return the GAP internal type number.

        This is only useful for libgap development purposes.

        OUTPUT: integer

        EXAMPLES::

            sage: x = libgap(1)
            sage: x._type_number()
            (0, b'integer')
        """
        n = TNUM_OBJ(self.value)
        name = TNAM_OBJ(self.value)
        return (n, name)

    def __dir__(self):
        """
        Customize tab completion.

        EXAMPLES::

            sage: G = libgap.DihedralGroup(4)
            sage: 'GeneratorsOfMagmaWithInverses' in dir(G)
            True
            sage: 'GeneratorsOfGroup' in dir(G)    # known bug
            False
            sage: x = libgap(1)
            sage: len(dir(x)) > 100
            True
        """
        from sage.libs.gap.operations import OperationInspector
        ops = OperationInspector(self).op_names()
        return dir(self.__class__) + ops

    def __getattr__(self, name):
        r"""
        Return functionoid implementing the function ``name``.

        EXAMPLES::

            sage: lst = libgap([])
            sage: 'Add' in dir(lst)    # This is why tab-completion works
            True
            sage: lst.Add(1)    # this is the syntactic sugar
            sage: lst
            [ 1 ]

        The above is equivalent to the following calls::

            sage: lst = libgap.eval('[]')
            sage: libgap.eval('Add') (lst, 1)
            sage: lst
            [ 1 ]

        TESTS::

            sage: lst.Adddddd(1)
            Traceback (most recent call last):
            ...
            AttributeError: 'Adddddd' is not defined in GAP

            sage: libgap.eval('some_name := 1')
            1
            sage: lst.some_name
            Traceback (most recent call last):
            ...
            AttributeError: 'some_name' does not define a GAP function
        """
        if name in ('__dict__', '_getAttributeNames', '_SageObject__custom_name', 'keys'):
            raise AttributeError('Python special name, not a GAP function.')
        try:
            proxy = make_GapElement_MethodProxy\
                (self.parent(), gap_eval(name), self)
        except GAPError:
            raise AttributeError(f"'{name}' is not defined in GAP")
        if not proxy.is_function():
            raise AttributeError(f"'{name}' does not define a GAP function")
        return proxy

    def __str__(self):
        r"""
        Return a string representation of ``self`` for printing.

        EXAMPLES::

            sage: libgap(0)
            0
            sage: print(libgap.eval(''))
            None
            sage: print(libgap('a'))
            a
            sage: print(libgap.eval('SymmetricGroup(3)'))
            SymmetricGroup( [ 1 .. 3 ] )
            sage: libgap(0).__str__()
            '0'
        """
        if self.value == NULL:
            return 'NULL'

        s = char_to_str(gap_element_str(self.value))
        return s.strip()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: libgap(0)
            0
            sage: libgap.eval('')
            sage: libgap('a')
            "a"
            sage: libgap.eval('SymmetricGroup(3)')
            Sym( [ 1 .. 3 ] )
            sage: libgap(0)._repr_()
            '0'
        """
        if self.value == NULL:
            return 'NULL'

        s = char_to_str(gap_element_repr(self.value))
        return s.strip()

    cpdef _set_compare_by_id(self):
        """
        Set comparison to compare by ``id``.

        By default, GAP is used to compare GAP objects. However,
        this is not defined for all GAP objects. To have GAP play
        nice with ``UniqueRepresentation``, comparison must always
        work. This method allows one to override the comparison to
        sort by the (unique) Python ``id``.

        Obviously it is a bad idea to change the comparison of objects
        after you have inserted them into a set/dict. You also must
        not mix GAP objects with different sort methods in the same
        container.

        EXAMPLES::

            sage: F1 = libgap.FreeGroup(['a'])
            sage: F2 = libgap.FreeGroup(['a'])
            sage: F1 < F2
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `<' on 2 arguments

            sage: F1._set_compare_by_id()
            sage: F1 != F2
            Traceback (most recent call last):
            ...
            ValueError: comparison style must be the same for both operands

            sage: F1._set_compare_by_id()
            sage: F2._set_compare_by_id()
            sage: F1 != F2
            True
        """
        self._compare_by_id = True

    cpdef _assert_compare_by_id(self):
        """
        Ensure that comparison is by ``id``.

        See :meth:`_set_compare_by_id`.

        OUTPUT:

        This method returns nothing. A :exc:`ValueError` is raised if
        :meth:`_set_compare_by_id` has not been called on this libgap
        object.

        EXAMPLES::

            sage: x = libgap.FreeGroup(1)
            sage: x._assert_compare_by_id()
            Traceback (most recent call last):
            ...
            ValueError: this requires a GAP object whose comparison is by "id"

            sage: x._set_compare_by_id()
            sage: x._assert_compare_by_id()
        """
        if not self._compare_by_id:
            raise ValueError('this requires a GAP object whose comparison is by "id"')

    def __hash__(self):
        """
        Make hashable.

        EXAMPLES::

            sage: hash(libgap(123))   # random output
            163512108404620371
        """
        return hash(str(self))

    cpdef _richcmp_(self, other, int op):
        """
        Compare ``self`` with ``other``.

        Uses the GAP comparison by default, or the Python ``id`` if
        :meth:`_set_compare_by_id` was called.

        OUTPUT:

        Boolean, depending on the comparison of ``self`` and
        ``other``.  Raises a :exc:`ValueError` if GAP does not support
        comparison of ``self`` and ``other``, unless
        :meth:`_set_compare_by_id` was called on both ``self`` and
        ``other``.

        EXAMPLES::

            sage: a = libgap(123)
            sage: a == a
            True
            sage: b = libgap('string')
            sage: a._richcmp_(b, 0)
            1
            sage: (a < b) or (a > b)
            True
            sage: a._richcmp_(libgap(123), 2)
            True

        GAP does not have a comparison function for two ``FreeGroup``
        objects. LibGAP signals this by raising a :exc:`ValueError` ::

            sage: F1 = libgap.FreeGroup(['a'])
            sage: F2 = libgap.FreeGroup(['a'])
            sage: F1 < F2
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `<' on 2 arguments

            sage: F1._set_compare_by_id()
            sage: F1 < F2
            Traceback (most recent call last):
            ...
            ValueError: comparison style must be the same for both operands

            sage: F1._set_compare_by_id()
            sage: F2._set_compare_by_id()
            sage: F1 < F2 or F1 > F2
            True

        Check that :issue:`26388` is fixed::

            sage: 1 > libgap(1)
            False
            sage: libgap(1) > 1
            False
            sage: 1 >= libgap(1)
            True
            sage: libgap(1) >= 1
            True
        """
        if self._compare_by_id != (<GapElement>other)._compare_by_id:
            raise ValueError("comparison style must be the same for both operands")
        if op==Py_LT:
            return self._compare_less(other)
        elif op==Py_LE:
            return self._compare_equal(other) or self._compare_less(other)
        elif op == Py_EQ:
            return self._compare_equal(other)
        elif op == Py_GT:
            return not self._compare_less(other) and not self._compare_equal(other)
        elif op == Py_GE:
            return not self._compare_less(other)
        elif op == Py_NE:
            return not self._compare_equal(other)
        else:
            assert False  # unreachable

    cdef bint _compare_equal(self, Element other) except -2:
        """
        Compare ``self`` with ``other``.

        Helper for :meth:`_richcmp_`

        EXAMPLES::

            sage: libgap(1) == libgap(1)   # indirect doctest
            True
        """
        if self._compare_by_id:
            return id(self) == id(other)
        cdef GapElement c_other = <GapElement>other
        sig_on()
        try:
            GAP_Enter()
            return GAP_EQ(self.value, c_other.value)
        finally:
            GAP_Leave()
            sig_off()

    cdef bint _compare_less(self, Element other) except -2:
        """
        Compare ``self`` with ``other``.

        Helper for :meth:`_richcmp_`

        EXAMPLES::

            sage: libgap(1) < libgap(2)   # indirect doctest
            True
        """
        if self._compare_by_id:
            return id(self) < id(other)
        cdef GapElement c_other = <GapElement>other
        sig_on()
        try:
            GAP_Enter()
            return GAP_LT(self.value, c_other.value)
        finally:
            GAP_Leave()
            sig_off()

    cpdef _add_(self, right):
        r"""
        Add two GapElement objects.

        EXAMPLES::

            sage: g1 = libgap(1)
            sage: g2 = libgap(2)
            sage: g1._add_(g2)
            3
            sage: g1 + g2    # indirect doctest
            3

            sage: libgap(1) + libgap.CyclicGroup(2)
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `+' on 2 arguments
        """
        cdef Obj result
        try:
            sig_GAP_Enter()
            sig_on()
            result = GAP_SUM(self.value, (<GapElement>right).value)
            sig_off()
        finally:
            GAP_Leave()
        return make_any_gap_element(self.parent(), result)

    cpdef _sub_(self, right):
        r"""
        Subtract two GapElement objects.

        EXAMPLES::

            sage: g1 = libgap(1)
            sage: g2 = libgap(2)
            sage: g1._sub_(g2)
            -1
            sage: g1 - g2    # indirect doctest
            -1

            sage: libgap(1) - libgap.CyclicGroup(2)
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found! ...
        """
        cdef Obj result
        try:
            sig_GAP_Enter()
            sig_on()
            result = GAP_DIFF(self.value, (<GapElement>right).value)
            sig_off()
        finally:
            GAP_Leave()
        return make_any_gap_element(self.parent(), result)

    cpdef _mul_(self, right):
        r"""
        Multiply two GapElement objects.

        EXAMPLES::

            sage: g1 = libgap(3)
            sage: g2 = libgap(5)
            sage: g1._mul_(g2)
            15
            sage: g1 * g2    # indirect doctest
            15

            sage: libgap(1) * libgap.CyclicGroup(2)
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `*' on 2 arguments
        """
        cdef Obj result
        try:
            sig_GAP_Enter()
            sig_on()
            result = GAP_PROD(self.value, (<GapElement>right).value)
            sig_off()
        finally:
            GAP_Leave()
        return make_any_gap_element(self.parent(), result)

    cpdef _div_(self, right):
        r"""
        Divide two GapElement objects.

        EXAMPLES::

            sage: g1 = libgap(3)
            sage: g2 = libgap(5)
            sage: g1._div_(g2)
            3/5
            sage: g1 / g2    # indirect doctest
            3/5

            sage: libgap(1) / libgap.CyclicGroup(2)
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `/' on 2 arguments

            sage: libgap(1) / libgap(0)
            Traceback (most recent call last):
            ...
            GAPError: Error, Rational operations: <divisor> must not be zero
        """
        cdef Obj result
        try:
            sig_GAP_Enter()
            sig_on()
            result = GAP_QUO(self.value, (<GapElement>right).value)
            sig_off()
        finally:
            GAP_Leave()
        return make_any_gap_element(self.parent(), result)

    cpdef _mod_(self, right):
        r"""
        Modulus of two GapElement objects.

        EXAMPLES::

            sage: g1 = libgap(5)
            sage: g2 = libgap(2)
            sage: g1 % g2
            1

            sage: libgap(1) % libgap.CyclicGroup(2)
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `mod' on 2 arguments
        """
        cdef Obj result
        try:
            sig_GAP_Enter()
            sig_on()
            result = GAP_MOD(self.value, (<GapElement>right).value)
            sig_off()
        finally:
            GAP_Leave()
        return make_any_gap_element(self.parent(), result)

    cpdef _pow_(self, other):
        r"""
        Exponentiation of two GapElement objects.

        EXAMPLES::

            sage: r = libgap(5) ^ 2; r
            25
            sage: parent(r)
            C library interface to GAP
            sage: r = 5 ^ libgap(2); r
            25
            sage: parent(r)
            C library interface to GAP
            sage: g, = libgap.CyclicGroup(5).GeneratorsOfGroup()
            sage: g ^ 5
            <identity> of ...

        TESTS:

        Check that this can be interrupted gracefully::

            sage: a, b = libgap.GL(1000, 3).GeneratorsOfGroup(); g = a * b
            sage: alarm(0.5); g ^ (2 ^ 10000)
            Traceback (most recent call last):
            ...
            AlarmInterrupt

            sage: libgap.CyclicGroup(2) ^ 2
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `^' on 2 arguments

            sage: libgap(3) ^ Infinity
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found! Error, no 1st choice
            method found for `InverseMutable' on 1 arguments
        """
        try:
            sig_GAP_Enter()
            sig_on()
            result = GAP_POW(self.value, (<GapElement>other).value)
            sig_off()
        finally:
            GAP_Leave()
        return make_any_gap_element(self._parent, result)

    cpdef _pow_int(self, other):
        """
        TESTS::

            sage: libgap(5)._pow_int(int(2))
            25
        """
        return self._pow_(self._parent(other))

    def is_function(self):
        """
        Return whether the wrapped GAP object is a function.

        OUTPUT: boolean

        EXAMPLES::

            sage: a = libgap.eval("NormalSubgroups")
            sage: a.is_function()
            True
            sage: a = libgap(2/3)
            sage: a.is_function()
            False
        """
        return IS_FUNC(self.value)

    def is_list(self):
        r"""
        Return whether the wrapped GAP object is a GAP List.

        OUTPUT: boolean

        EXAMPLES::

            sage: libgap.eval('[1, 2,,,, 5]').is_list()
            True
            sage: libgap.eval('3/2').is_list()
            False
        """
        return GAP_IsList(self.value)

    def is_record(self):
        r"""
        Return whether the wrapped GAP object is a GAP record.

        OUTPUT: boolean

        EXAMPLES::

            sage: libgap.eval('[1, 2,,,, 5]').is_record()
            False
            sage: libgap.eval('rec(a:=1, b:=3)').is_record()
            True
        """
        return GAP_IsRecord(self.value)

    cpdef is_bool(self):
        r"""
        Return whether the wrapped GAP object is a GAP boolean.

        OUTPUT: boolean

        EXAMPLES::

            sage: libgap(True).is_bool()
            True
        """
        libgap = self.parent()
        cdef GapElement r_sage = libgap.IsBool(self)
        cdef Obj r_gap = r_sage.value
        return r_gap == GAP_True

    def is_string(self):
        r"""
        Return whether the wrapped GAP object is a GAP string.

        OUTPUT: boolean

        EXAMPLES::

            sage: libgap('this is a string').is_string()
            True
        """
        return IS_STRING(self.value)

    def is_permutation(self):
        r"""
        Return whether the wrapped GAP object is a GAP permutation.

        OUTPUT: boolean

        EXAMPLES::

            sage: perm = libgap.PermList( libgap([1,5,2,3,4]) );  perm
            (2,5,4,3)
            sage: perm.is_permutation()
            True
            sage: libgap('this is a string').is_permutation()
            False
        """
        return (TNUM_OBJ(self.value) == T_PERM2 or
                TNUM_OBJ(self.value) == T_PERM4)

    def sage(self):
        r"""
        Return the Sage equivalent of the :class:`GapElement`.

        EXAMPLES::

            sage: libgap(1).sage()
            1
            sage: type(_)
            <class 'sage.rings.integer.Integer'>

            sage: libgap(3/7).sage()
            3/7
            sage: type(_)
            <class 'sage.rings.rational.Rational'>

            sage: libgap.eval('5 + 7*E(3)').sage()
            7*zeta3 + 5

            sage: libgap(Infinity).sage()
            +Infinity
            sage: libgap(-Infinity).sage()
            -Infinity

            sage: libgap(True).sage()
            True
            sage: libgap(False).sage()
            False
            sage: type(_)
            <... 'bool'>

            sage: libgap('this is a string').sage()
            'this is a string'
            sage: type(_)
            <... 'str'>

            sage: x = libgap.Integers.Indeterminate("x")

            sage: p = x^2 - 2*x + 3
            sage: p.sage()
            x^2 - 2*x + 3
            sage: p.sage().parent()
            Univariate Polynomial Ring in x over Integer Ring

            sage: p = x^-2 + 3*x
            sage: p.sage()
            x^-2 + 3*x
            sage: p.sage().parent()
            Univariate Laurent Polynomial Ring in x over Integer Ring

            sage: p = (3 * x^2 + x) / (x^2 - 2)
            sage: p.sage()
            (3*x^2 + x)/(x^2 - 2)
            sage: p.sage().parent()
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring

            sage: G0 = libgap.SymplecticGroup(4,2)
            sage: P = G0.IsomorphismFpGroup().Range()
            sage: G = P.sage()
            sage: G.gap() == P
            True

            sage: F0 = libgap.FreeGroup(2)
            sage: F = F0.sage()
            sage: F.gap() is F0
            True

        TESTS:

        Check :issue:`30496`::

            sage: x = libgap.Integers.Indeterminate("x")
            sage: p = x^2 - 2*x
            sage: p.sage()
            x^2 - 2*x
        """
        if self.value is NULL:
            return None

        if self.IsInfinity():
            from sage.rings.infinity import Infinity
            return Infinity

        elif self.IsNegInfinity():
            from sage.rings.infinity import Infinity
            return -Infinity

        elif self.IsUnivariateRationalFunction():
            var = self.IndeterminateOfUnivariateRationalFunction().String()
            var = var.sage()
            num, den, val = self.CoefficientsOfUnivariateRationalFunction()
            num = num.sage()
            den = den.sage()
            val = val.sage()
            base_ring = cm.common_parent(*(num + den))

            if self.IsUnivariatePolynomial():
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                R = PolynomialRing(base_ring, var)
                x = R.gen()
                return x**val * R(num)

            elif self.IsLaurentPolynomial():
                from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
                R = LaurentPolynomialRing(base_ring, var)
                x = R.gen()
                return x**val * R(num)

            else:
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                R = PolynomialRing(base_ring, var)
                x = R.gen()
                return x**val * R(num) / R(den)

        elif self.IsList():
            # May be a list-like collection of some other type of GapElements
            # that we can convert
            return [item.sage() for item in self.AsList()]

        elif self.IsFreeGroup():
            from sage.groups.free_group import FreeGroup_class
            names = tuple(str(g) for g in self.GeneratorsOfGroup())
            return FreeGroup_class(names, gap_group=self)

        elif self.IsFpGroup():
            from sage.groups.free_group import FreeGroup
            from sage.groups.finitely_presented import FinitelyPresentedGroup
            # names = tuple(str(g).replace(".", "_") for g in self.FreeGroupOfFpGroup().GeneratorsOfGroup())
            F = self.FreeGroupOfFpGroup().sage()
            relations = tuple(F(rel.LetterRepAssocWord().sage())
                              for rel in self.RelatorsOfFpGroup())
            return FinitelyPresentedGroup(F, relations, libgap_fpgroup=self)

        raise NotImplementedError('cannot construct equivalent Sage object')


############################################################################
### GapElement_Integer #####################################################
############################################################################

cdef GapElement_Integer make_GapElement_Integer(parent, Obj obj):
    r"""
    Turn a Gap integer object into a GapElement_Integer Sage object.

    EXAMPLES::

        sage: libgap(123)
        123
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Integer'>
    """
    cdef GapElement_Integer r = GapElement_Integer.__new__(GapElement_Integer)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Integer(GapElement):
    r"""
    Derived class of GapElement for GAP integers.

    EXAMPLES::

        sage: i = libgap(123)
        sage: type(i)
        <class 'sage.libs.gap.element.GapElement_Integer'>
        sage: ZZ(i)
        123
    """

    def is_C_int(self):
        r"""
        Return whether the wrapped GAP object is a immediate GAP integer.

        An immediate integer is one that is stored as a C integer, and
        is subject to the usual size limits. Larger integers are
        stored in GAP as GMP integers.

        OUTPUT: boolean

        EXAMPLES::

            sage: n = libgap(1)
            sage: type(n)
            <class 'sage.libs.gap.element.GapElement_Integer'>
            sage: n.is_C_int()
            True
            sage: n.IsInt()
            true

            sage: N = libgap(2^130)
            sage: type(N)
            <class 'sage.libs.gap.element.GapElement_Integer'>
            sage: N.is_C_int()
            False
            sage: N.IsInt()
            true
        """
        return GAP_IsSmallInt(self.value)

    def _rational_(self):
        r"""
        EXAMPLES::

            sage: QQ(libgap(1))  # indirect doctest
            1
            sage: QQ(libgap(-2**200)) == -2**200
            True
        """
        return self.sage(ring=QQ)

    def sage(self, ring=None):
        r"""
        Return the Sage equivalent of the :class:`GapElement_Integer`.

        - ``ring`` -- integer ring or ``None`` (default); if not
          specified, the default Sage integer ring is used

        OUTPUT: a Sage integer

        EXAMPLES::

            sage: libgap([ 1, 3, 4 ]).sage()
            [1, 3, 4]
            sage: all( x in ZZ for x in _ )
            True

            sage: libgap(132).sage(ring=IntegerModRing(13))
            2
            sage: parent(_)
            Ring of integers modulo 13

        TESTS::

            sage: libgap(0).sage()
            0
            sage: large = libgap.eval('2^130');  large
            1361129467683753853853498429727072845824
            sage: large.sage()
            1361129467683753853853498429727072845824

            sage: huge = libgap.eval('10^9999');  huge     # gap abbreviates very long ints
            <integer 100...000 (10000 digits)>
            sage: huge.sage().ndigits()
            10000
        """
        cdef UInt* x
        cdef Int size
        cdef int c_sign
        cdef int c_size
        cdef mpz_t output
        if ring is None:
            ring = ZZ
        try:
            GAP_Enter()
            if self.is_C_int():
                return ring(GAP_ValueInt(self.value))
            else:
                # gap integers are stored as a mp_limb_t
                size = GAP_SizeInt(self.value) # count limbs and extract sign
                if size > 0:
                    c_sign = 1
                    c_size = size
                else: # Must have size < 0, or else self.value == 0 and self.is_C_int() == True
                    c_sign = -1
                    c_size = -size
                x = GAP_AddrInt(self.value) # pointer to limbs
                mpz_roinit_n(output, <mp_limb_t *>x, c_size)
                return ring(c_sign*mpz_get_pylong(output))
        finally:
            GAP_Leave()

    _integer_ = sage

    def __int__(self):
        r"""
        TESTS::

            sage: int(libgap(3))
            3
            sage: type(_)
            <class 'int'>

            sage: int(libgap(2)**128)
            340282366920938463463374607431768211456
            sage: type(_)
            <class 'int'>
        """
        return self.sage(ring=int)

    def __index__(self):
        r"""
        TESTS:

        Check that gap integers can be used as indices (:issue:`23878`)::

            sage: s = 'abcd'
            sage: s[libgap(1)]
            'b'
        """
        return int(self)


##########################################################################
### GapElement_Float #####################################################
##########################################################################

cdef GapElement_Float make_GapElement_Float(parent, Obj obj):
    r"""
    Turn a Gap macfloat object into a GapElement_Float Sage object.

    EXAMPLES::

        sage: libgap(123.5)
        123.5
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Float'>
    """
    cdef GapElement_Float r = GapElement_Float.__new__(GapElement_Float)
    r._initialize(parent, obj)
    return r

cdef class GapElement_Float(GapElement):
    r"""
    Derived class of GapElement for GAP floating point numbers.

    EXAMPLES::

        sage: i = libgap(123.5)
        sage: type(i)
        <class 'sage.libs.gap.element.GapElement_Float'>
        sage: RDF(i)
        123.5
        sage: float(i)
        123.5

    TESTS::

        sage: a = RDF.random_element()
        sage: libgap(a).sage() == a
        True
    """
    def sage(self, ring=None):
        r"""
        Return the Sage equivalent of the :class:`GapElement_Float`.

        - ``ring`` -- a floating point field or ``None`` (default); if not
          specified, the default Sage ``RDF`` is used

        OUTPUT: a Sage double precision floating point number

        EXAMPLES::

            sage: a = libgap.eval("Float(3.25)").sage()
            sage: a
            3.25
            sage: parent(a)
            Real Double Field
        """
        if ring is None:
            ring = RDF
        return ring(GAP_ValueMacFloat(self.value))

    def __float__(self):
        r"""
        TESTS::

            sage: float(libgap.eval("Float(3.5)"))
            3.5
        """
        return GAP_ValueMacFloat(self.value)


############################################################################
#  GapElement_IntegerMod  ##################################################
############################################################################

cdef GapElement_IntegerMod make_GapElement_IntegerMod(parent, Obj obj):
    r"""
    Turn a Gap integer object into a :class:`GapElement_IntegerMod` Sage object.

    EXAMPLES::

        sage: n = IntegerModRing(123)(13)
        sage: libgap(n)
        ZmodnZObj( 13, 123 )
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_IntegerMod'>
    """
    cdef GapElement_IntegerMod r = GapElement_IntegerMod.__new__(GapElement_IntegerMod)
    r._initialize(parent, obj)
    return r

cdef class GapElement_IntegerMod(GapElement):
    r"""
    Derived class of GapElement for GAP integers modulo an integer.

    EXAMPLES::

        sage: n = IntegerModRing(123)(13)
        sage: i = libgap(n)
        sage: type(i)
        <class 'sage.libs.gap.element.GapElement_IntegerMod'>
    """

    cpdef GapElement_Integer lift(self):
        """
        Return an integer lift.

        OUTPUT:

        A :class:`GapElement_Integer` that equals ``self`` in the
        integer mod ring.

        EXAMPLES::

            sage: n = libgap.eval('One(ZmodnZ(123)) * 13')
            sage: n.lift()
            13
            sage: type(_)
            <class 'sage.libs.gap.element.GapElement_Integer'>
        """
        return self.Int()

    def sage(self, ring=None):
        r"""
        Return the Sage equivalent of the :class:`GapElement_IntegerMod`.

        INPUT:

        - ``ring`` -- Sage integer mod ring or ``None`` (default); if
          not specified, a suitable integer mod ring is used
          automatically

        OUTPUT: a Sage integer modulo another integer

        EXAMPLES::

            sage: n = libgap.eval('One(ZmodnZ(123)) * 13')
            sage: n.sage()
            13
            sage: parent(_)
            Ring of integers modulo 123
        """
        if ring is None:
            # ring = self.DefaultRing().sage()
            characteristic = self.Characteristic().sage()
            ring = ZZ.quotient_ring(characteristic)
        return ring(self.lift())

    def _integer_(self, R):
        r"""
        TESTS::

            sage: n = libgap.ZmodnZObj(13, 123)
            sage: ZZ(n)
            13
            sage: ZZ(-n)
            110
            sage: ZZ(0*n)
            0
        """
        return self.lift()._integer_(R)

############################################################################
### GapElement_FiniteField #####################################################
############################################################################

cdef GapElement_FiniteField make_GapElement_FiniteField(parent, Obj obj):
    r"""
    Turn a GAP finite field object into a :class:`GapElement_FiniteField` Sage object.

    EXAMPLES::

        sage: libgap.eval('Z(5)^2')
        Z(5)^2
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_FiniteField'>
    """
    cdef GapElement_FiniteField r = GapElement_FiniteField.__new__(GapElement_FiniteField)
    r._initialize(parent, obj)
    return r


cdef class GapElement_FiniteField(GapElement):
    r"""
    Derived class of GapElement for GAP finite field elements.

    EXAMPLES::

        sage: libgap.eval('Z(5)^2')
        Z(5)^2
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_FiniteField'>
    """

    cpdef GapElement_Integer lift(self):
        """
        Return an integer lift.

        OUTPUT:

        The smallest positive :class:`GapElement_Integer` that equals
        ``self`` in the prime finite field.

        EXAMPLES::

            sage: n = libgap.eval('Z(5)^2')
            sage: n.lift()
            4
            sage: type(_)
            <class 'sage.libs.gap.element.GapElement_Integer'>

            sage: n = libgap.eval('Z(25)')
            sage: n.lift()
            Traceback (most recent call last):
            TypeError: not in prime subfield
        """
        degree = self.DegreeFFE().sage()
        if degree == 1:
            return self.IntFFE()
        else:
            raise TypeError('not in prime subfield')

    def sage(self, ring=None, var='a'):
        r"""
        Return the Sage equivalent of the :class:`GapElement_FiniteField`.

        INPUT:

        - ``ring`` -- a Sage finite field or ``None`` (default). The
          field to return ``self`` in. If not specified, a suitable
          finite field will be constructed.

        OUTPUT:

        A Sage finite field element. The isomorphism is chosen such
        that the Gap ``PrimitiveRoot()`` maps to the Sage
        :meth:`~sage.rings.finite_rings.finite_field_prime_modn.multiplicative_generator`.

        EXAMPLES::

            sage: n = libgap.eval('Z(25)^2')
            sage: n.sage()
            a + 3
            sage: parent(_)
            Finite Field in a of size 5^2

            sage: n.sage(ring=GF(5))
            Traceback (most recent call last):
            ...
            ValueError: the given ring is incompatible ...

        TESTS::

            sage: n = libgap.eval('Z(2^4)^2 + Z(2^4)^1 + Z(2^4)^0')
            sage: n
            Z(2^2)^2
            sage: n.sage()
            a + 1
            sage: parent(_)
            Finite Field in a of size 2^2
            sage: n.sage(ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: the given ring is incompatible ...
            sage: n.sage(ring=CC)
            Traceback (most recent call last):
            ...
            ValueError: the given ring is incompatible ...
            sage: n.sage(ring=GF(5))
            Traceback (most recent call last):
            ...
            ValueError: the given ring is incompatible ...
            sage: n.sage(ring=GF(2^3))
            Traceback (most recent call last):
            ...
            ValueError: the given ring is incompatible ...
            sage: n.sage(ring=GF(2^2, 'a'))
            a + 1
            sage: n.sage(ring=GF(2^4, 'a'))
            a^2 + a + 1
            sage: n.sage(ring=GF(2^8, 'a'))
            a^7 + a^6 + a^4 + a^2 + a + 1
            sage: (n^3).sage()
            1

        Check that :issue:`23153` is fixed::

            sage: n = libgap.eval('Z(2^4)^2 + Z(2^4)^1 + Z(2^4)^0')
            sage: n.sage(ring=GF(2^4, 'a'))
            a^2 + a + 1
        """
        deg = self.DegreeFFE().sage()
        char = self.Characteristic().sage()
        if ring is None:
            from sage.rings.finite_rings.finite_field_constructor import GF
            ring = GF(char**deg, name=var)
        elif not (ring.is_field() and ring.is_finite() and
                  ring.characteristic() == char and ring.degree() % deg == 0):
            raise ValueError(('the given ring is incompatible (must be a '
                              'finite field of characteristic {} and degree '
                              'divisible by {})').format(char, deg))

        if self.IsOne():
            return ring.one()
        if deg == 1 and char == ring.characteristic():
            return ring(self.lift().sage())
        else:
            gap_field = make_GapElement_Ring(self.parent(), gap_eval(ring._gap_init_()))
            exp = self.LogFFE(gap_field.PrimitiveRoot())
            return ring.multiplicative_generator() ** exp.sage()

    def __int__(self):
        r"""
        TESTS::

            sage: int(libgap.eval("Z(53)"))
            2
        """
        return int(self.Int())

    def _integer_(self, R):
        r"""
        TESTS::

            sage: ZZ(libgap.eval("Z(53)"))
            2
        """
        return R(self.Int())


############################################################################
### GapElement_Cyclotomic #####################################################
############################################################################

cdef GapElement_Cyclotomic make_GapElement_Cyclotomic(parent, Obj obj):
    r"""
    Turn a Gap cyclotomic object into a :class:`GapElement_Cyclotomic` Sage
    object.

    EXAMPLES::

        sage: libgap.eval('E(3)')
        E(3)
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Cyclotomic'>
    """
    cdef GapElement_Cyclotomic r = GapElement_Cyclotomic.__new__(GapElement_Cyclotomic)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Cyclotomic(GapElement):
    r"""
    Derived class of GapElement for GAP universal cyclotomics.

    EXAMPLES::

        sage: libgap.eval('E(3)')
        E(3)
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Cyclotomic'>
    """

    def sage(self, ring=None):
        r"""
        Return the Sage equivalent of the :class:`GapElement_Cyclotomic`.

        INPUT:

        - ``ring`` -- a Sage cyclotomic field or ``None``
          (default); if not specified, a suitable minimal cyclotomic
          field will be constructed

        OUTPUT: a Sage cyclotomic field element

        EXAMPLES::

            sage: n = libgap.eval('E(3)')
            sage: n.sage()
            zeta3
            sage: parent(_)
            Cyclotomic Field of order 3 and degree 2

            sage: n.sage(ring=CyclotomicField(6))
            zeta6 - 1

            sage: libgap.E(3).sage(ring=CyclotomicField(3))
            zeta3
            sage: libgap.E(3).sage(ring=CyclotomicField(6))
            zeta6 - 1

        TESTS:

        Check that :issue:`15204` is fixed::

            sage: libgap.E(3).sage(ring=UniversalCyclotomicField())
            E(3)
            sage: libgap.E(3).sage(ring=CC)
            -0.500000000000000 + 0.866025403784439*I
        """
        if ring is None:
            conductor = self.Conductor()
            from sage.rings.number_field.number_field import CyclotomicField
            ring = CyclotomicField(conductor.sage())
        else:
            try:
                conductor = ring._n()
            except AttributeError:
                from sage.rings.number_field.number_field import CyclotomicField
                conductor = self.Conductor()
                cf = CyclotomicField(conductor.sage())
                return ring(cf(self.CoeffsCyc(conductor).sage()))
        coeff = self.CoeffsCyc(conductor).sage()
        return ring(coeff)


############################################################################
### GapElement_Rational ####################################################
############################################################################

cdef GapElement_Rational make_GapElement_Rational(parent, Obj obj):
    r"""
    Turn a Gap Rational number (of type ``Obj``) into a Cython ``GapElement_Rational``.

    EXAMPLES::

        sage: libgap(123/456)
        41/152
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Rational'>
    """
    cdef GapElement_Rational r = GapElement_Rational.__new__(GapElement_Rational)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Rational(GapElement):
    r"""
    Derived class of GapElement for GAP rational numbers.

    EXAMPLES::

        sage: r = libgap(123/456)
        sage: type(r)
        <class 'sage.libs.gap.element.GapElement_Rational'>
    """
    def _rational_(self):
        r"""
        EXAMPLES::

            sage: r = libgap(-1/3)
            sage: QQ(r)  # indirect doctest
            -1/3
            sage: QQ(libgap(2**300 / 3**300)) == 2**300 / 3**300
            True
        """
        return self.sage(ring=QQ)

    def sage(self, ring=None):
        r"""
        Return the Sage equivalent of the :class:`GapElement`.

        INPUT:

        - ``ring`` -- the Sage rational ring or ``None`` (default); if
          not specified, the rational ring is used automatically

        OUTPUT: a Sage rational number

        EXAMPLES::

            sage: r = libgap(123/456);  r
            41/152
            sage: type(_)
            <class 'sage.libs.gap.element.GapElement_Rational'>
            sage: r.sage()
            41/152
            sage: type(_)
            <class 'sage.rings.rational.Rational'>
        """
        if ring is None:
            ring = ZZ
        libgap = self.parent()
        return libgap.NumeratorRat(self).sage(ring=ring) / libgap.DenominatorRat(self).sage(ring=ring)


############################################################################
### GapElement_Ring #####################################################
############################################################################

cdef GapElement_Ring make_GapElement_Ring(parent, Obj obj):
    r"""
    Turn a Gap integer object into a :class:`GapElement_Ring` Sage
    object.

    EXAMPLES::

        sage: libgap(GF(5))
        GF(5)
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Ring'>
    """
    cdef GapElement_Ring r = GapElement_Ring.__new__(GapElement_Ring)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Ring(GapElement):
    r"""
    Derived class of GapElement for GAP rings (parents of ring elements).

    EXAMPLES::

        sage: i = libgap(ZZ)
        sage: type(i)
        <class 'sage.libs.gap.element.GapElement_Ring'>
    """

    def ring_integer(self):
        """
        Construct the Sage integers.

        EXAMPLES::

            sage: libgap.eval('Integers').ring_integer()
            Integer Ring
        """
        return ZZ

    def ring_rational(self):
        """
        Construct the Sage rationals.

        EXAMPLES::

            sage: libgap.eval('Rationals').ring_rational()
            Rational Field
        """
        return ZZ.fraction_field()

    def ring_integer_mod(self):
        """
        Construct a Sage integer mod ring.

        EXAMPLES::

            sage: libgap.eval('ZmodnZ(15)').ring_integer_mod()
            Ring of integers modulo 15
        """
        characteristic = self.Characteristic().sage()
        return ZZ.quotient_ring(characteristic)

    def ring_finite_field(self, var='a'):
        """
        Construct an integer ring.

        EXAMPLES::

            sage: libgap.GF(3,2).ring_finite_field(var='A')
            Finite Field in A of size 3^2
        """
        size = self.Size().sage()
        from sage.rings.finite_rings.finite_field_constructor import GF
        return GF(size, name=var)

    def ring_cyclotomic(self):
        """
        Construct an integer ring.

        EXAMPLES::

            sage: libgap.CyclotomicField(6).ring_cyclotomic()
            Cyclotomic Field of order 3 and degree 2
        """
        conductor = self.Conductor()
        from sage.rings.number_field.number_field import CyclotomicField
        return CyclotomicField(conductor.sage())

    def ring_polynomial(self):
        """
        Construct a polynomial ring.

        EXAMPLES::

            sage: B = libgap(QQ['x'])
            sage: B.ring_polynomial()
            Univariate Polynomial Ring in x over Rational Field

            sage: B = libgap(ZZ['x','y'])
            sage: B.ring_polynomial()
            Multivariate Polynomial Ring in x, y over Integer Ring
        """
        base_ring = self.CoefficientsRing().sage()
        vars = [x.String().sage()
                for x in self.IndeterminatesOfPolynomialRing()]
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(base_ring, vars)

    def sage(self, **kwds):
        r"""
        Return the Sage equivalent of the :class:`GapElement_Ring`.

        INPUT:

        - ``**kwds`` -- keywords that are passed on to the ``ring_``
          method

        OUTPUT: a Sage ring

        EXAMPLES::

            sage: libgap.eval('Integers').sage()
            Integer Ring

            sage: libgap.eval('Rationals').sage()
            Rational Field

            sage: libgap.eval('ZmodnZ(15)').sage()
            Ring of integers modulo 15

            sage: libgap.GF(3,2).sage(var='A')
            Finite Field in A of size 3^2

            sage: libgap.CyclotomicField(6).sage()
            Cyclotomic Field of order 3 and degree 2

            sage: libgap(QQ['x','y']).sage()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        if self.IsField():
            if self.IsRationals():
                return self.ring_rational(**kwds)
            if self.IsCyclotomicField():
                return self.ring_cyclotomic(**kwds)
            if self.IsFinite():
                return self.ring_finite_field(**kwds)
        else:
            if self.IsIntegers():
                return self.ring_integer(**kwds)
            if self.IsFinite():
                return self.ring_integer_mod(**kwds)
            if self.IsPolynomialRing():
                return self.ring_polynomial(**kwds)
        raise NotImplementedError('cannot convert GAP ring to Sage')


############################################################################
### GapElement_Boolean #####################################################
############################################################################

cdef GapElement_Boolean make_GapElement_Boolean(parent, Obj obj):
    r"""
    Turn a Gap Boolean number (of type ``Obj``) into a Cython ``GapElement_Boolean``.

    EXAMPLES::

        sage: libgap(True)
        true
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Boolean'>
    """
    cdef GapElement_Boolean r = GapElement_Boolean.__new__(GapElement_Boolean)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Boolean(GapElement):
    r"""
    Derived class of GapElement for GAP boolean values.

    EXAMPLES::

        sage: b = libgap(True)
        sage: type(b)
        <class 'sage.libs.gap.element.GapElement_Boolean'>
    """

    def sage(self):
        r"""
        Return the Sage equivalent of the :class:`GapElement`.

        OUTPUT:

        A Python boolean if the values is either true or false. GAP
        booleans can have the third value ``Fail``, in which case a
        :exc:`ValueError` is raised.

        EXAMPLES::

            sage: b = libgap.eval('true');  b
            true
            sage: type(_)
            <class 'sage.libs.gap.element.GapElement_Boolean'>
            sage: b.sage()
            True
            sage: type(_)
            <... 'bool'>

            sage: libgap.eval('fail')
            fail
            sage: _.sage()
            Traceback (most recent call last):
            ...
            ValueError: the GAP boolean value "fail" cannot be represented in Sage
        """
        if self.value == GAP_True:
            return True
        if self.value == GAP_False:
            return False
        raise ValueError('the GAP boolean value "fail" cannot be represented in Sage')

    def __bool__(self):
        """
        Check that the boolean is "true".

        This is syntactic sugar for using libgap. See the examples below.

        OUTPUT: boolean

        EXAMPLES::

            sage: gap_bool = [libgap.eval('true'), libgap.eval('false'), libgap.eval('fail')]
            sage: for x in gap_bool:
            ....:     if x:     # this calls __bool__
            ....:         print("{} {}".format(x, type(x)))
            true <class 'sage.libs.gap.element.GapElement_Boolean'>

            sage: for x in gap_bool:
            ....:     if not x:     # this calls __bool__
            ....:         print("{} {}".format( x, type(x)))
            false <class 'sage.libs.gap.element.GapElement_Boolean'>
            fail <class 'sage.libs.gap.element.GapElement_Boolean'>
        """
        return self.value == GAP_True


############################################################################
### GapElement_String ####################################################
############################################################################

cdef GapElement_String make_GapElement_String(parent, Obj obj):
    r"""
    Turn a Gap String (of type ``Obj``) into a Cython ``GapElement_String``.

    EXAMPLES::

        sage: libgap('this is a string')
        "this is a string"
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_String'>
    """
    cdef GapElement_String r = GapElement_String.__new__(GapElement_String)
    r._initialize(parent, obj)
    return r


cdef class GapElement_String(GapElement):
    r"""
    Derived class of GapElement for GAP strings.

    EXAMPLES::

        sage: s = libgap('string')
        sage: type(s)
        <class 'sage.libs.gap.element.GapElement_String'>
        sage: s
        "string"
        sage: print(s)
        string
    """
    def __str__(self):
        r"""
        Convert this :class:`GapElement_String` to a Python string.

        OUTPUT: a Python string

        EXAMPLES::

            sage: s = libgap.eval(' "string" '); s
            "string"
            sage: type(_)
            <class 'sage.libs.gap.element.GapElement_String'>
            sage: str(s)
            'string'
            sage: s.sage()
            'string'
            sage: type(_)
            <class 'str'>
        """
        s = char_to_str(GAP_CSTR_STRING(self.value))
        return s

    sage = __str__

############################################################################
### GapElement_Function ####################################################
############################################################################

cdef GapElement_Function make_GapElement_Function(parent, Obj obj):
    r"""
    Turn a Gap C function object (of type ``Obj``) into a Cython ``GapElement_Function``.

    INPUT:

    - ``parent`` -- the parent of the new :class:`GapElement`

    - ``obj`` -- a GAP function object

    OUTPUT: a :class:`GapElement_Function` instance

    EXAMPLES::

        sage: libgap.CycleLength
        <Gap function "CycleLength">
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Function'>
    """
    cdef GapElement_Function r = GapElement_Function.__new__(GapElement_Function)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Function(GapElement):
    r"""
    Derived class of GapElement for GAP functions.

    EXAMPLES::

        sage: f = libgap.Cycles
        sage: type(f)
        <class 'sage.libs.gap.element.GapElement_Function'>
    """

    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: libgap.Orbits
            <Gap function "Orbits">
        """
        libgap = self.parent()
        name = libgap.NameFunction(self)
        s = '<Gap function "'+name.sage()+'">'
        return s

    def __call__(self, *args):
        """
        Call syntax for functions.

        INPUT:

        - ``*args`` -- arguments. Will be converted to `GapElement` if
          they are not already of this type

        OUTPUT:

        A :class:`GapElement` encapsulating the functions return
        value, or ``None`` if it does not return anything.

        EXAMPLES::

            sage: a = libgap.NormalSubgroups
            sage: b = libgap.SymmetricGroup(4)
            sage: libgap.collect()
            sage: a
            <Gap function "NormalSubgroups">
            sage: b
            Sym( [ 1 .. 4 ] )
            sage: [x.StructureDescription() for x in sorted(a(b))]
            ["1", "S4", "A4", "C2 x C2"]

            sage: libgap.eval("a := NormalSubgroups")
            <Gap function "NormalSubgroups">
            sage: libgap.eval("b := SymmetricGroup(4)")
            Sym( [ 1 .. 4 ] )
            sage: libgap.collect()
            sage: [x.StructureDescription() for x in sorted(libgap.eval('a') (libgap.eval('b')))]
            ["1", "S4", "A4", "C2 x C2"]

            sage: a = libgap.eval('a')
            sage: b = libgap.eval('b')
            sage: libgap.collect()
            sage: [x.StructureDescription() for x in sorted(a(b))]
            ["1", "S4", "A4", "C2 x C2"]

        Not every ``GapElement`` is callable::

            sage: f = libgap(3)
            sage: f()
            Traceback (most recent call last):
            ...
            TypeError: 'sage.libs.gap.element.GapElement_Integer' object is not callable

        We illustrate appending to a list which returns None::

            sage: a = libgap([]); a
            [  ]
            sage: a.Add(5); a
            [ 5 ]
            sage: a.Add(10); a
            [ 5, 10 ]

        TESTS::

            sage: s = libgap.Sum
            sage: s(libgap([1,2]))
            3
            sage: s(libgap(1), libgap(2))
            Traceback (most recent call last):
            ...
            GAPError: Error, no method found!
            Error, no 1st choice method found for `SumOp' on 2 arguments

            sage: for i in range(100):
            ....:     rnd = [ randint(-10,10) for i in range(randint(0,7)) ]
            ....:     # compute the sum in GAP
            ....:     _ = libgap.Sum(rnd)
            ....:     try:
            ....:         libgap.Sum(*rnd)
            ....:         print('This should have triggered a ValueError')
            ....:         print('because Sum needs a list as argument')
            ....:     except ValueError:
            ....:         pass

            sage: libgap_exec = libgap.eval("Exec")
            sage: libgap_exec('echo hello from the shell')
            hello from the shell
        """
        cdef Obj result = NULL
        cdef Obj arg_list
        cdef int n = len(args)
        cdef volatile Obj v2

        if n > 0 and n <= 3:
            libgap = self.parent()
            a = [x if isinstance(x, GapElement) else libgap(x) for x in args]

        try:
            sig_GAP_Enter()
            sig_on()
            if n == 0:
                result = GAP_CallFunc0Args(self.value)
            elif n == 1:
                result = GAP_CallFunc1Args(self.value,
                                           (<GapElement>a[0]).value)
            elif n == 2:
                result = GAP_CallFunc2Args(self.value,
                                           (<GapElement>a[0]).value,
                                           (<GapElement>a[1]).value)
            elif n == 3:
                v2 = (<GapElement>a[2]).value
                result = GAP_CallFunc3Args(self.value,
                                           (<GapElement>a[0]).value,
                                           (<GapElement>a[1]).value,
                                           v2)
            else:
                arg_list = make_gap_list(args)
                result = GAP_CallFuncList(self.value, arg_list)
            sig_off()
        finally:
            GAP_Leave()
        if result == NULL:
            # We called a procedure that does not return anything
            return None
        return make_any_gap_element(self.parent(), result)

    def _instancedoc_(self):
        r"""
        Return the help string.

        EXAMPLES::

            sage: f = libgap.CyclicGroup
            sage: 'constructs  the  cyclic  group' in f.__doc__
            True

        You would get the full help by typing ``f?`` in the command line.
        """
        libgap = self.parent()
        from sage.interfaces.gap import gap
        return gap.help(libgap.NameFunction(self).sage(), pager=False)


############################################################################
### GapElement_MethodProxy #################################################
############################################################################

cdef GapElement_MethodProxy make_GapElement_MethodProxy(parent, Obj function, GapElement base_object):
    r"""
    Turn a Gap C rec object (of type ``Obj``) into a Cython ``GapElement_Record``.

    This class implement syntactic sugar so that you can write
    ``gapelement.f()`` instead of ``libgap.f(gapelement)`` for any GAP
    function ``f``.

    INPUT:

    - ``parent`` -- the parent of the new :class:`GapElement`

    - ``obj`` -- a GAP function object

    - ``base_object`` -- the first argument to be inserted into the function

    OUTPUT: a :class:`GapElement_MethodProxy` instance

    EXAMPLES::

        sage: lst = libgap([])
        sage: type( lst.Add )
        <class 'sage.libs.gap.element.GapElement_MethodProxy'>
    """
    cdef GapElement_MethodProxy r = GapElement_MethodProxy.__new__(GapElement_MethodProxy)
    r._initialize(parent, function)
    r.first_argument = base_object
    return r


cdef class GapElement_MethodProxy(GapElement_Function):
    r"""
    Helper class returned by ``GapElement.__getattr__``.

    Derived class of GapElement for GAP functions. Like its parent,
    you can call instances to implement function call syntax. The only
    difference is that a fixed first argument is prepended to the
    argument list.

    EXAMPLES::

        sage: lst = libgap([])
        sage: lst.Add
        <Gap function "Add">
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_MethodProxy'>
        sage: lst.Add(1)
        sage: lst
        [ 1 ]
    """

    def __call__(self, *args):
        """
        Call syntax for methods.

        This method is analogous to
        :meth:`GapElement_Function.__call__`, except that it inserts a
        fixed :class:`GapElement` in the first slot of the function.

        INPUT:

        - ``*args`` -- arguments. Will be converted to `GapElement` if
          they are not already of this type

        OUTPUT:

        A :class:`GapElement` encapsulating the functions return
        value, or ``None`` if it does not return anything.

        EXAMPLES::

            sage: lst = libgap.eval('[1,,3]')
            sage: lst.Add.__call__(4)
            sage: lst.Add(5)
            sage: lst
            [ 1,, 3, 4, 5 ]
        """
        if len(args) > 0:
            return GapElement_Function.__call__(self, * ([self.first_argument] + list(args)))
        else:
            return GapElement_Function.__call__(self, self.first_argument)


############################################################################
### GapElement_List ########################################################
############################################################################

cdef GapElement_List make_GapElement_List(parent, Obj obj):
    r"""
    Turn a Gap C List object (of type ``Obj``) into a Cython ``GapElement_List``.

    EXAMPLES::

        sage: libgap([0, 2, 3])
        [ 0, 2, 3 ]
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_List'>
    """
    cdef GapElement_List r = GapElement_List.__new__(GapElement_List)
    r._initialize(parent, obj)
    return r


cdef class GapElement_List(GapElement):
    r"""
    Derived class of GapElement for GAP Lists.

    .. NOTE::

        Lists are indexed by `0..len(l)-1`, as expected from
        Python. This differs from the GAP convention where lists start
        at `1`.

    EXAMPLES::

        sage: lst = libgap.SymmetricGroup(3).List(); lst
        [ (), (1,3), (1,2,3), (2,3), (1,3,2), (1,2) ]
        sage: type(lst)
        <class 'sage.libs.gap.element.GapElement_List'>
        sage: len(lst)
        6
        sage: lst[3]
        (2,3)

    We can easily convert a Gap ``List`` object into a Python ``list``::

        sage: list(lst)
        [(), (1,3), (1,2,3), (2,3), (1,3,2), (1,2)]
        sage: type(_)
        <... 'list'>

    Range checking is performed::

        sage: lst[10]
        Traceback (most recent call last):
        ...
        IndexError: index out of range.
    """

    def __bool__(self):
        r"""
        Return ``True`` if the list is non-empty, as with Python ``list``s.

        EXAMPLES::

            sage: lst = libgap.eval('[1,,,4]')
            sage: bool(lst)
            True
            sage: lst = libgap.eval('[]')
            sage: bool(lst)
            False
        """
        return bool(len(self))

    def __len__(self):
        r"""
        Return the length of the list.

        OUTPUT: integer

        EXAMPLES::

            sage: lst = libgap.eval('[1,,,4]')   # a sparse list
            sage: len(lst)
            4
        """
        return GAP_LenList(self.value)

    def __getitem__(self, i):
        r"""
        Return the ``i``-th element of the list.

        As usual in Python, indexing starts at `0` and not at `1` (as
        in GAP). This can also be used with multi-indices.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        The ``i``-th element as a :class:`GapElement`.

        EXAMPLES::

            sage: lst = libgap.eval('["first",,,"last"]')   # a sparse list
            sage: lst[0]
            "first"

            sage: l = libgap.eval('[ [0, 1], [2, 3] ]')
            sage: l[0,0]
            0
            sage: l[0,1]
            1
            sage: l[1,0]
            2
            sage: l[0,2]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: l[2,0]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: l[0,0,0]
            Traceback (most recent call last):
            ...
            ValueError: too many indices
        """
        cdef int j
        cdef Obj obj = self.value

        if isinstance(i, tuple):
            for j in i:
                if not GAP_IsList(obj):
                    raise ValueError('too many indices')
                if j < 0 or j >= GAP_LenList(obj):
                    raise IndexError('index out of range')
                obj = ELM_LIST(obj, j+1)

        else:
            j = i
            if j < 0 or j >= GAP_LenList(obj):
                raise IndexError('index out of range.')
            obj = ELM_LIST(obj, j+1)

        return make_any_gap_element(self.parent(), obj)

    def __setitem__(self, i, elt):
        r"""
        Set the `i`-th item of this list.

        EXAMPLES::

            sage: l = libgap.eval('[0, 1]')
            sage: l
            [ 0, 1 ]
            sage: l[0] = 3
            sage: l
            [ 3, 1 ]

        Contrarily to Python lists, setting an element beyond the limit extends the list::

            sage: l[12] = -2
            sage: l
            [ 3, 1,,,,,,,,,,, -2 ]

        This function also handles multi-indices::

            sage: l = libgap.eval('[[[0,1],[2,3]],[[4,5], [6,7]]]')
            sage: l[0,1,0] = -18
            sage: l
            [ [ [ 0, 1 ], [ -18, 3 ] ], [ [ 4, 5 ], [ 6, 7 ] ] ]
            sage: l[0,0,0,0]
            Traceback (most recent call last):
            ...
            ValueError: too many indices

        Assignment to immutable objects gives error::

            sage: l = libgap([0,1])
            sage: u = l.deepcopy(0)
            sage: u[0] = 5
            Traceback (most recent call last):
            ...
            TypeError: immutable Gap object does not support item assignment

        TESTS::

            sage: m = libgap.eval('[[0,0],[0,0]]')
            sage: m[0,0] = 1
            sage: m[0,1] = 2
            sage: m[1,0] = 3
            sage: m[1,1] = 4
            sage: m
            [ [ 1, 2 ], [ 3, 4 ] ]
        """
        if not IS_MUTABLE_OBJ(self.value):
            raise TypeError('immutable Gap object does not support item assignment')

        cdef int j
        cdef Obj obj = self.value

        if isinstance(i, tuple):
            for j in i[:-1]:
                if not GAP_IsList(obj):
                    raise ValueError('too many indices')
                if j < 0 or j >= GAP_LenList(obj):
                    raise IndexError('index out of range')
                obj = ELM_LIST(obj, j+1)
            if not GAP_IsList(obj):
                raise ValueError('too many indices')
            j = i[-1]
        else:
            j = i

        if j < 0:
            raise IndexError('index out of range.')

        cdef GapElement celt
        if isinstance(elt, GapElement):
            celt = <GapElement> elt
        else:
            celt= self.parent()(elt)

        GAP_AssList(obj, j+1, celt.value)

    def sage(self, **kwds):
        r"""
        Return the Sage equivalent of the :class:`GapElement`.

        OUTPUT: a Python list

        EXAMPLES::

            sage: libgap([ 1, 3, 4 ]).sage()
            [1, 3, 4]
            sage: all( x in ZZ for x in _ )
            True
        """
        return [x.sage(**kwds) for x in self]

    def matrix(self, ring=None):
        """
        Return the list as a matrix.

        GAP does not have a special matrix data type, they are just
        lists of lists. This function converts a GAP list of lists to
        a Sage matrix.

        OUTPUT: a Sage matrix

        EXAMPLES::

            sage: F = libgap.GF(4)
            sage: a = F.PrimitiveElement()
            sage: m = libgap([[a,a^0],[0*a,a^2]]); m
            [ [ Z(2^2), Z(2)^0 ],
              [ 0*Z(2), Z(2^2)^2 ] ]
            sage: m.IsMatrix()
            true
            sage: matrix(m)
            [    a     1]
            [    0 a + 1]
            sage: matrix(GF(4,'B'), m)
            [    B     1]
            [    0 B + 1]

            sage: M = libgap.eval('SL(2,GF(5))').GeneratorsOfGroup()[1]
            sage: type(M)
            <class 'sage.libs.gap.element.GapElement_List'>
            sage: M[0][0]
            Z(5)^2
            sage: M.IsMatrix()
            true
            sage: M.matrix()
            [4 1]
            [4 0]
        """
        if not self.IsMatrix():
            raise ValueError('not a GAP matrix')
        entries = self.Flat()
        n = self.Length().sage()
        m = len(entries) // n
        if len(entries) % n != 0:
            raise ValueError('not a rectangular list of lists')
        from sage.matrix.matrix_space import MatrixSpace
        if ring is None:
            ring = entries.DefaultRing().sage()
        MS = MatrixSpace(ring, n, m)
        return MS([ring(x) for x in entries])

    _matrix_ = matrix

    def vector(self, ring=None):
        """
        Return the list as a vector.

        GAP does not have a special vector data type, they are just
        lists. This function converts a GAP list to a Sage vector.

        OUTPUT: a Sage vector

        EXAMPLES::

            sage: F = libgap.GF(4)
            sage: a = F.PrimitiveElement()
            sage: m = libgap([0*a, a, a^3, a^2]); m
            [ 0*Z(2), Z(2^2), Z(2)^0, Z(2^2)^2 ]
            sage: type(m)
            <class 'sage.libs.gap.element.GapElement_List'>
            sage: m[3]
            Z(2^2)^2
            sage: vector(m)
            (0, a, 1, a + 1)
            sage: vector(GF(4,'B'), m)
            (0, B, 1, B + 1)
        """
        if not self.IsVector():
            raise ValueError('not a GAP vector')
        from sage.modules.free_module_element import vector
        entries = self.Flat()
        n = self.Length().sage()
        if ring is None:
            ring = entries.DefaultRing().sage()
        return vector(ring, n, self.sage(ring=ring))

    _vector_ = vector


############################################################################
### GapElement_Permutation #################################################
############################################################################


cdef GapElement_Permutation make_GapElement_Permutation(parent, Obj obj):
    r"""
    Turn a Gap C permutation object (of type ``Obj``) into a Cython ``GapElement_Permutation``.

    EXAMPLES::

        sage: libgap.eval('(1,3,2)(4,5,8)')
        (1,3,2)(4,5,8)
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Permutation'>
    """
    cdef GapElement_Permutation r = GapElement_Permutation.__new__(GapElement_Permutation)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Permutation(GapElement):
    r"""
    Derived class of GapElement for GAP permutations.

    .. NOTE::

        Permutations in GAP act on the numbers starting with 1.

    EXAMPLES::

        sage: perm = libgap.eval('(1,5,2)(4,3,8)')
        sage: type(perm)
        <class 'sage.libs.gap.element.GapElement_Permutation'>
    """

    def sage(self, parent=None):
        r"""
        Return the Sage equivalent of the :class:`GapElement`.

        If the permutation group is given as parent, this method is
        *much* faster.

        EXAMPLES::

            sage: perm_gap = libgap.eval('(1,5,2)(4,3,8)');  perm_gap
            (1,5,2)(3,8,4)
            sage: perm_gap.sage()
            [5, 1, 8, 3, 2, 6, 7, 4]
            sage: type(_)
            <class 'sage.combinat.permutation.StandardPermutations_all_with_category.element_class'>
            sage: perm_gap.sage(PermutationGroup([(1,2),(1,2,3,4,5,6,7,8)]))
            (1,5,2)(3,8,4)
            sage: type(_)
            <class 'sage.groups.perm_gps.permgroup_element.PermutationGroupElement'>
        """
        cdef PermutationGroupElement one_c

        libgap = self.parent()
        lst = libgap.ListPerm(self)

        if parent is None:
            return Permutation(lst.sage(), check=False)
        else:
            return parent.one()._generate_new_GAP(lst)

############################################################################
### GapElement_Record ######################################################
############################################################################

cdef GapElement_Record make_GapElement_Record(parent, Obj obj):
    r"""
    Turn a Gap C rec object (of type ``Obj``) into a Cython ``GapElement_Record``.

    EXAMPLES::

        sage: libgap.eval('rec(a:=0, b:=2, c:=3)')
        rec( a := 0, b := 2, c := 3 )
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement_Record'>
    """
    cdef GapElement_Record r = GapElement_Record.__new__(GapElement_Record)
    r._initialize(parent, obj)
    return r


cdef class GapElement_Record(GapElement):
    r"""
    Derived class of GapElement for GAP records.

    EXAMPLES::

        sage: rec = libgap.eval('rec(a:=123, b:=456)')
        sage: type(rec)
        <class 'sage.libs.gap.element.GapElement_Record'>
        sage: len(rec)
        2
        sage: rec['a']
        123

    We can easily convert a Gap ``rec`` object into a Python ``dict``::

        sage: dict(rec)
        {'a': 123, 'b': 456}
        sage: type(_)
        <... 'dict'>

    Range checking is performed::

        sage: rec['no_such_element']
        Traceback (most recent call last):
        ...
        GAPError: Error, Record Element: '<rec>.no_such_element' must have an assigned value
    """

    def __len__(self):
        r"""
        Return the length of the record.

        OUTPUT: integer; the number of entries in the record

        EXAMPLES::

            sage: rec = libgap.eval('rec(a:=123, b:=456, S3:=SymmetricGroup(3))')
            sage: len(rec)
            3
        """
        return LEN_PREC(self.value)

    def __iter__(self):
        r"""
        Iterate over the elements of the record.

        OUTPUT: a :class:`GapElement_RecordIterator`

        EXAMPLES::

            sage: rec = libgap.eval('rec(a:=123, b:=456)')
            sage: iter = rec.__iter__()
            sage: type(iter)
            <class 'sage.libs.gap.element.GapElement_RecordIterator'>
            sage: sorted(rec)
            [('a', 123), ('b', 456)]
        """
        return GapElement_RecordIterator(self)

    cpdef UInt record_name_to_index(self, name) noexcept:
        r"""
        Convert string to GAP record index.

        INPUT:

        - ``py_name`` -- a python string

        OUTPUT:

        A ``UInt``, which is a GAP hash of the string. If this is the
        first time the string is encountered, a new integer is
        returned(!)

        EXAMPLES::

            sage: rec = libgap.eval('rec(first:=123, second:=456)')
            sage: rec.record_name_to_index('first')   # random output
            1812
            sage: rec.record_name_to_index('no_such_name') # random output
            3776
        """
        name = str_to_bytes(name)
        return RNamName(name)

    def __getitem__(self, name):
        r"""
        Return the ``name``-th element of the GAP record.

        INPUT:

        - ``name`` -- string

        OUTPUT: the record element labelled by ``name`` as a :class:`GapElement`

        EXAMPLES::

            sage: rec = libgap.eval('rec(first:=123, second:=456)')
            sage: rec['first']
            123
        """
        cdef UInt i = self.record_name_to_index(name)
        cdef Obj result
        sig_on()
        try:
            GAP_Enter()
            result = ELM_REC(self.value, i)
        finally:
            GAP_Leave()
            sig_off()
        return make_any_gap_element(self.parent(), result)

    def sage(self):
        r"""
        Return the Sage equivalent of the :class:`GapElement`.

        EXAMPLES::

            sage: libgap.eval('rec(a:=1, b:=2)').sage()
            {'a': 1, 'b': 2}
            sage: all( isinstance(key,str) and val in ZZ for key,val in _.items() )
            True

            sage: rec = libgap.eval('rec(a:=123, b:=456, Sym3:=SymmetricGroup(3))')
            sage: rec.sage()
            {'Sym3': NotImplementedError('cannot construct equivalent Sage object'...),
             'a': 123,
             'b': 456}
        """
        result = {}
        for key, val in self:
            try:
                val = val.sage()
            except Exception as ex:
                val = ex
            result[key] = val
        return result


cdef class GapElement_RecordIterator():
    r"""
    Iterator for :class:`GapElement_Record`.

    Since Cython does not support generators yet, we implement the
    older iterator specification with this auxiliary class.

    INPUT:

    - ``rec`` -- the :class:`GapElement_Record` to iterate over

    EXAMPLES::

        sage: rec = libgap.eval('rec(a:=123, b:=456)')
        sage: sorted(rec)
        [('a', 123), ('b', 456)]
        sage: dict(rec)
        {'a': 123, 'b': 456}
    """

    def __cinit__(self, rec):
        r"""
        The Cython constructor.

        INPUT:

        - ``rec`` -- the :class:`GapElement_Record` to iterate over

        EXAMPLES::

            sage: libgap.eval('rec(a:=123, b:=456)')
            rec( a := 123, b := 456 )
        """
        self.rec = rec
        self.i = 1

    def __next__(self):
        r"""
        Return the next element in the record.

        OUTPUT:

        A tuple ``(key, value)`` where ``key`` is a string and
        ``value`` is the corresponding :class:`GapElement`.

        EXAMPLES::

            sage: rec = libgap.eval('rec(a:=123, b:=456)')
            sage: iter = rec.__iter__()
            sage: a = iter.__next__()
            sage: b = next(iter)
            sage: sorted([a, b])
            [('a', 123), ('b', 456)]
        """
        cdef UInt i = self.i
        if i>len(self.rec):
            raise StopIteration
        # note the abs: negative values mean the rec keys are not sorted
        key_index = abs(GET_RNAM_PREC(self.rec.value, i))
        key = char_to_str(GAP_CSTR_STRING(NAME_RNAM(key_index)))
        cdef Obj result = GET_ELM_PREC(self.rec.value, i)
        val = make_any_gap_element(self.rec.parent(), result)
        self.i += 1
        return (key, val)


# Add support for _instancedoc_
from sage.misc.instancedoc import instancedoc
instancedoc(GapElement_Function)
instancedoc(GapElement_MethodProxy)
