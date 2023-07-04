#*****************************************************************************
#       Copyright (C) 2009 Carl Witty <Carl.Witty@gmail.com>
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

"""Implements different data storage types."""

from __future__ import print_function, absolute_import

from .utils import je, reindent_lines as ri


class StorageType(object):
    r"""
    A StorageType specifies the C types used to deal with values of a
    given type.

    We currently support three categories of types.

    First are the "simple" types.  These are types where: the
    representation is small, functions expect arguments to be passed
    by value, and the C/C++ assignment operator works.  This would
    include built-in C types (long, float, etc.) and small structs
    (like gsl_complex).

    Second is 'PyObject*'.  This is just like a simple type, except
    that we have to incref/decref at appropriate places.

    Third is "auto-reference" types.  This is how
    GMP/MPIR/MPFR/MPFI/FLINT types work.  For these types, functions
    expect arguments to be passed by reference, and the C assignment
    operator does not do what we want.  In addition, they take
    advantage of a quirk in C (where arrays are automatically
    converted to pointers) to automatically pass arguments by
    reference.

    Support for further categories would not be difficult to add (such
    as reference-counted types other than PyObject*, or
    pass-by-reference types that don't use the GMP auto-reference
    trick), if we ever run across a use for them.
    """

    def __init__(self):
        r"""
        Initialize an instance of StorageType.

        This sets several properties:

        class_member_declarations:
        A string giving variable declarations that must be members of any
        wrapper class using this type.

        class_member_initializations:
        A string initializing the class_member_declarations; will be
        inserted into the __init__ method of any wrapper class using this
        type.

        local_declarations:
        A string giving variable declarations that must be local variables
        in Cython methods using this storage type.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.class_member_declarations
            ''
            sage: ty_double.class_member_initializations
            ''
            sage: ty_double.local_declarations
            ''
            sage: ty_mpfr.class_member_declarations
            'cdef RealField_class domain\n'
            sage: ty_mpfr.class_member_initializations
            "self.domain = args['domain']\n"
            sage: ty_mpfr.local_declarations
            'cdef RealNumber rn\n'
        """
        self.class_member_declarations = ''
        self.class_member_initializations = ''
        self.local_declarations = ''

    def cheap_copies(self):
        r"""
        Returns True or False, depending on whether this StorageType
        supports cheap copies -- whether it is cheap to copy values of
        this type from one location to another.  This is true for
        primitive types, and for types like PyObject* (where you are only
        copying a pointer, and possibly changing some reference counts).
        It is false for types like mpz_t and mpfr_t, where copying values
        can involve arbitrarily much work (including memory allocation).

        The practical effect is that if cheap_copies is True,
        instructions with outputs of this type write the results into
        local variables, and the results are then copied to their
        final locations.  If cheap_copies is False, then the addresses
        of output locations are passed into the instruction and the
        instruction writes outputs directly in the final location.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.cheap_copies()
            True
            sage: ty_python.cheap_copies()
            True
            sage: ty_mpfr.cheap_copies()
            False
        """
        return False

    def python_refcounted(self):
        r"""
        Says whether this storage type is a Python type, so we need to
        use INCREF/DECREF.

        (If we needed to support any non-Python refcounted types, it
        might be better to make this object-oriented and have methods
        like "generate an incref" and "generate a decref".  But as
        long as we only support Python, this way is probably simpler.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.python_refcounted()
            False
            sage: ty_python.python_refcounted()
            True
        """
        return False

    def cython_decl_type(self):
        r"""
        Give the Cython type for a single value of this type (as a string).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.cython_decl_type()
            'double'
            sage: ty_python.cython_decl_type()
            'object'
            sage: ty_mpfr.cython_decl_type()
            'mpfr_t'
        """
        return self.c_decl_type()

    def cython_array_type(self):
        r"""
        Give the Cython type for referring to an array of values of
        this type (as a string).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.cython_array_type()
            'double*'
            sage: ty_python.cython_array_type()
            'PyObject**'
            sage: ty_mpfr.cython_array_type()
            'mpfr_t*'
        """
        return self.c_ptr_type()

    def needs_cython_init_clear(self):
        r"""
        Says whether values/arrays of this type need to be initialized
        before use and cleared before the underlying memory is freed.

        (We could remove this method, always call .cython_init() to
        generate initialization code, and just let .cython_init()
        generate empty code if no initialization is required; that would
        generate empty loops, which are ugly and potentially might not
        be optimized away.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.needs_cython_init_clear()
            False
            sage: ty_mpfr.needs_cython_init_clear()
            True
            sage: ty_python.needs_cython_init_clear()
            True
        """
        return False

    def c_decl_type(self):
        r"""
        Give the C type for a single value of this type (as a string).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.c_decl_type()
            'double'
            sage: ty_python.c_decl_type()
            'PyObject*'
            sage: ty_mpfr.c_decl_type()
            'mpfr_t'
        """
        raise NotImplementedError

    def c_ptr_type(self):
        r"""
        Give the C type for a pointer to this type (as a reference to
        either a single value or an array) (as a string).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.c_ptr_type()
            'double*'
            sage: ty_python.c_ptr_type()
            'PyObject**'
            sage: ty_mpfr.c_ptr_type()
            'mpfr_t*'
        """
        return self.c_decl_type() + '*'

    def c_reference_type(self):
        r"""
        Give the C type which should be used for passing a reference
        to a single value in a call. This is used as the type for the
        return value.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.c_reference_type()
            'double*'
            sage: ty_python.c_reference_type()
            'PyObject**'
        """
        return self.c_ptr_type()

    def c_local_type(self):
        r"""
        Give the C type used for a value of this type inside an
        instruction.  For assignable/cheap_copy types, this is the
        same as c_decl_type; for auto-reference types, this is the
        pointer type.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.c_local_type()
            'double'
            sage: ty_python.c_local_type()
            'PyObject*'
            sage: ty_mpfr.c_local_type()
            'mpfr_ptr'
        """
        raise NotImplementedError

    def assign_c_from_py(self, c, py):
        r"""
        Given a Cython variable/array reference/etc. of this storage type,
        and a Python expression, generate code to assign to the Cython
        variable from the Python expression.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.assign_c_from_py('foo', 'bar')
            'foo = bar'
            sage: ty_python.assign_c_from_py('foo[i]', 'bar[j]')
            'foo[i] = <PyObject *>bar[j]; Py_INCREF(foo[i])'
            sage: ty_mpfr.assign_c_from_py('foo', 'bar')
            'rn = self.domain(bar)\nmpfr_set(foo, rn.value, MPFR_RNDN)'
        """
        return je("{{ c }} = {{ py }}", c=c, py=py)

    def declare_chunk_class_members(self, name):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for a memory chunk with this storage type
        and the given name.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.declare_chunk_class_members('args')
            '    cdef int _n_args\n    cdef mpfr_t* _args\n'
        """
        return je(ri(0,
            """
            {# XXX Variables here (and everywhere, really) should actually be Py_ssize_t #}
                cdef int _n_{{ name }}
                cdef {{ myself.cython_array_type() }} _{{ name }}
            """), myself=self, name=name)

    def alloc_chunk_data(self, name, len):
        r"""
        Return a string allocating the memory for the class members for
        a memory chunk with this storage type and the given name.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: print(ty_mpfr.alloc_chunk_data('args', 'MY_LENGTH'))
                    self._n_args = MY_LENGTH
                    self._args = <mpfr_t*>check_allocarray(self._n_args, sizeof(mpfr_t))
                    for i in range(MY_LENGTH):
                        mpfr_init2(self._args[i], self.domain.prec())
            <BLANKLINE>
        """
        return je(ri(0,
            """
                    self._n_{{ name }} = {{ len }}
                    self._{{ name }} = <{{ myself.c_ptr_type() }}>check_allocarray(self._n_{{ name }}, sizeof({{ myself.c_decl_type() }}))
            {% if myself.needs_cython_init_clear() %}
                    for i in range({{ len }}):
                        {{ myself.cython_init('self._%s[i]' % name) }}
            {% endif %}
            """), myself=self, name=name, len=len)

    def dealloc_chunk_data(self, name):
        r"""
        Return a string to be put in the __dealloc__ method of a
        wrapper class using a memory chunk with this storage type, to
        deallocate the corresponding class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: print(ty_double.dealloc_chunk_data('args'))
                    if self._args:
                        sig_free(self._args)
            <BLANKLINE>
            sage: print(ty_mpfr.dealloc_chunk_data('constants'))
                    if self._constants:
                        for i in range(self._n_constants):
                            mpfr_clear(self._constants[i])
                        sig_free(self._constants)
            <BLANKLINE>
        """
        return je(ri(0, """
                    if self._{{ name }}:
            {%     if myself.needs_cython_init_clear() %}
                        for i in range(self._n_{{ name }}):
                            {{ myself.cython_clear('self._%s[i]' % name) }}
            {%     endif %}
                        sig_free(self._{{ name }})
            """), myself=self, name=name)


class StorageTypeAssignable(StorageType):
    r"""
    StorageTypeAssignable is a subtype of StorageType that deals with
    types with cheap copies, like primitive types and PyObject*.
    """

    def __init__(self, ty):
        r"""
        Initializes the property type (the C/Cython name for this type),
        as well as the properties described in the documentation for
        StorageType.__init__.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.class_member_declarations
            ''
            sage: ty_double.class_member_initializations
            ''
            sage: ty_double.local_declarations
            ''
            sage: ty_double.type
            'double'
            sage: ty_python.type
            'PyObject*'
        """
        StorageType.__init__(self)
        self.type = ty

    def cheap_copies(self):
        r"""
        Returns True or False, depending on whether this StorageType
        supports cheap copies -- whether it is cheap to copy values of
        this type from one location to another.  (See StorageType.cheap_copies
        for more on this property.)

        Since having cheap copies is essentially the definition of
        StorageTypeAssignable, this always returns True.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.cheap_copies()
            True
            sage: ty_python.cheap_copies()
            True
        """
        return True

    def c_decl_type(self):
        r"""
        Give the C type for a single value of this type (as a string).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.c_decl_type()
            'double'
            sage: ty_python.c_decl_type()
            'PyObject*'
        """
        return self.type

    def c_local_type(self):
        r"""
        Give the C type used for a value of this type inside an
        instruction.  For assignable/cheap_copy types, this is the
        same as c_decl_type; for auto-reference types, this is the
        pointer type.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_double.c_local_type()
            'double'
            sage: ty_python.c_local_type()
            'PyObject*'
        """
        return self.type


class StorageTypeSimple(StorageTypeAssignable):
    r"""
    StorageTypeSimple is a subtype of StorageTypeAssignable that deals
    with non-reference-counted types with cheap copies, like primitive
    types.  As of yet, it has no functionality differences from
    StorageTypeAssignable.
    """
    pass

ty_int = StorageTypeSimple('int')
ty_double = StorageTypeSimple('double')


class StorageTypeDoubleComplex(StorageTypeSimple):
    r"""
    This is specific to the complex double type. It behaves exactly
    like a StorageTypeSimple in C, but needs a little help to do
    conversions in Cython.

    This uses functions defined in CDFInterpreter, and is for use in
    that context.
    """
    def assign_c_from_py(self, c, py):
        """
        sage: from sage_setup.autogen.interpreters import ty_double_complex
        sage: ty_double_complex.assign_c_from_py('z_c', 'z_py')
        'z_c = CDE_to_dz(z_py)'
        """
        return je("{{ c }} = CDE_to_dz({{ py }})", c=c, py=py)

ty_double_complex = StorageTypeDoubleComplex('double_complex')


class StorageTypePython(StorageTypeAssignable):
    r"""
    StorageTypePython is a subtype of StorageTypeAssignable that deals
    with Python objects.

    Just allocating an array full of PyObject* leads to problems,
    because the Python garbage collector must be able to get to every
    Python object, and it wouldn't know how to get to these arrays.
    So we allocate the array as a Python list, but then we immediately
    pull the ob_item out of it and deal only with that from then on.

    We often leave these lists with NULL entries.  This is safe for
    the garbage collector and the deallocator, which is all we care
    about; but it would be unsafe to provide Python-level access to
    these lists.

    There is one special thing about StorageTypePython: memory that is
    used by the interpreter as scratch space (for example, the stack)
    must be cleared after each call (so we don't hold on to
    potentially-large objects and waste memory).  Since we have to do
    this anyway, the interpreter gains a tiny bit of speed by assuming
    that the scratch space is cleared on entry; for example, when
    pushing a value onto the stack, it doesn't bother to XDECREF the
    previous value because it's always NULL.
    """

    def __init__(self):
        r"""
        Initializes the properties described in the documentation
        for StorageTypeAssignable.__init__.  The type is always
        'PyObject*'.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.class_member_declarations
            ''
            sage: ty_python.class_member_initializations
            ''
            sage: ty_python.local_declarations
            ''
            sage: ty_python.type
            'PyObject*'
        """
        super(StorageTypePython, self).__init__('PyObject*')

    def python_refcounted(self):
        r"""
        Says whether this storage type is a Python type, so we need to
        use INCREF/DECREF.

        Returns True.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.python_refcounted()
            True
        """
        return True

    def cython_decl_type(self):
        r"""
        Give the Cython type for a single value of this type (as a string).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.cython_decl_type()
            'object'
        """
        return 'object'

    def declare_chunk_class_members(self, name):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for a memory chunk with this storage type
        and the given name.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.declare_chunk_class_members('args')
            '    cdef object _list_args\n    cdef int _n_args\n    cdef PyObject** _args\n'
        """
        return je(ri(4,
            """
            cdef object _list_{{ name }}
            cdef int _n_{{ name }}
            cdef {{ myself.cython_array_type() }} _{{ name }}
            """), myself=self, name=name)

    def alloc_chunk_data(self, name, len):
        r"""
        Return a string allocating the memory for the class members for
        a memory chunk with this storage type and the given name.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: print(ty_python.alloc_chunk_data('args', 'MY_LENGTH'))
                    self._n_args = MY_LENGTH
                    self._list_args = PyList_New(self._n_args)
                    self._args = (<PyListObject *>self._list_args).ob_item
            <BLANKLINE>
        """
        return je(ri(8,
            """
            self._n_{{ name }} = {{ len }}
            self._list_{{ name }} = PyList_New(self._n_{{ name }})
            self._{{ name }} = (<PyListObject *>self._list_{{ name }}).ob_item
            """), myself=self, name=name, len=len)

    def dealloc_chunk_data(self, name):
        r"""
        Return a string to be put in the __dealloc__ method of a
        wrapper class using a memory chunk with this storage type, to
        deallocate the corresponding class members.

        Our array was allocated as a Python list; this means we actually
        don't need to do anything to deallocate it.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.dealloc_chunk_data('args')
            ''
        """
        return ''

    def needs_cython_init_clear(self):
        r"""
        Says whether values/arrays of this type need to be initialized
        before use and cleared before the underlying memory is freed.

        Returns True.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.needs_cython_init_clear()
            True
        """
        return True

    def assign_c_from_py(self, c, py):
        r"""
        Given a Cython variable/array reference/etc. of this storage type,
        and a Python expression, generate code to assign to the Cython
        variable from the Python expression.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.assign_c_from_py('foo[i]', 'bar[j]')
            'foo[i] = <PyObject *>bar[j]; Py_INCREF(foo[i])'
        """
        return je("""{{ c }} = <PyObject *>{{ py }}; Py_INCREF({{ c }})""",
                  c=c, py=py)

    def cython_init(self, loc):
        r"""
        Generates code to initialize a variable (or array reference)
        holding a PyObject*.  Sets it to NULL.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.cython_init('foo[i]')
            'foo[i] = NULL'
        """
        return je("{{ loc }} = NULL", loc=loc)

    def cython_clear(self, loc):
        r"""
        Generates code to clear a variable (or array reference) holding
        a PyObject*.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.cython_clear('foo[i]')
            'Py_CLEAR(foo[i])'
        """
        return je("Py_CLEAR({{ loc }})", loc=loc)

ty_python = StorageTypePython()


class StorageTypeAutoReference(StorageType):
    r"""
    StorageTypeAutoReference is a subtype of StorageType that deals with
    types in the style of GMP/MPIR/MPFR/MPFI/FLINT, where copies are
    not cheap, functions expect arguments to be passed by reference,
    and the API takes advantage of the C quirk where arrays are
    automatically converted to pointers to automatically pass
    arguments by reference.
    """
    def __init__(self, decl_ty, ref_ty):
        r"""
        Initializes the properties decl_type and ref_type (the C type
        names used when declaring variables and function parameters,
        respectively), as well as the properties described in
        the documentation for StorageType.__init__.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.class_member_declarations
            'cdef RealField_class domain\n'
            sage: ty_mpfr.class_member_initializations
            "self.domain = args['domain']\n"
            sage: ty_mpfr.local_declarations
            'cdef RealNumber rn\n'
            sage: ty_mpfr.decl_type
            'mpfr_t'
            sage: ty_mpfr.ref_type
            'mpfr_ptr'
        """
        super(StorageTypeAutoReference, self).__init__()
        self.decl_type = decl_ty
        self.ref_type = ref_ty

    def c_decl_type(self):
        r"""
        Give the C type for a single value of this type (as a string).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.c_decl_type()
            'mpfr_t'
        """
        return self.decl_type

    def c_local_type(self):
        r"""
        Give the C type used for a value of this type inside an
        instruction.  For assignable/cheap_copy types, this is the
        same as c_decl_type; for auto-reference types, this is the
        pointer type.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.c_local_type()
            'mpfr_ptr'
        """
        return self.ref_type

    def c_reference_type(self):
        r"""
        Give the C type which should be used for passing a reference
        to a single value in a call. This is used as the type for the
        return value.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.c_reference_type()
            'mpfr_t'
        """
        return self.decl_type

    def needs_cython_init_clear(self):
        r"""
        Says whether values/arrays of this type need to be initialized
        before use and cleared before the underlying memory is freed.

        All known examples of auto-reference types do need a special
        initialization call, so this always returns True.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.needs_cython_init_clear()
            True
        """
        return True


class StorageTypeMPFR(StorageTypeAutoReference):
    r"""
    StorageTypeMPFR is a subtype of StorageTypeAutoReference that deals
    the MPFR's mpfr_t type.

    For any given program that we're interpreting, ty_mpfr can only
    refer to a single precision.  An interpreter that needs to use
    two precisions of mpfr_t in the same program should instantiate two
    separate instances of StorageTypeMPFR.  (Interpreters that need
    to handle arbitrarily many precisions in the same program are not
    handled at all.)
    """

    def __init__(self, id=''):
        r"""
        Initializes the id property, as well as the properties described
        in the documentation for StorageTypeAutoReference.__init__.

        The id property is used if you want to have an interpreter
        that handles two instances of StorageTypeMPFR (that is,
        handles mpfr_t variables at two different precisions
        simultaneously).  It's a string that's used to generate
        variable names that don't conflict.  (The id system has
        never actually been used, so bugs probably remain.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.class_member_declarations
            'cdef RealField_class domain\n'
            sage: ty_mpfr.class_member_initializations
            "self.domain = args['domain']\n"
            sage: ty_mpfr.local_declarations
            'cdef RealNumber rn\n'
            sage: ty_mpfr.decl_type
            'mpfr_t'
            sage: ty_mpfr.ref_type
            'mpfr_ptr'

        TESTS::

            sage: ty_mpfr2 = StorageTypeMPFR(id='_the_second')
            sage: ty_mpfr2.class_member_declarations
            'cdef RealField_class domain_the_second\n'
            sage: ty_mpfr2.class_member_initializations
            "self.domain_the_second = args['domain_the_second']\n"
            sage: ty_mpfr2.local_declarations
            'cdef RealNumber rn_the_second\n'
        """

        super(StorageTypeMPFR, self).__init__('mpfr_t', 'mpfr_ptr')
        self.id = id
        self.class_member_declarations = "cdef RealField_class domain%s\n" % self.id
        self.class_member_initializations = \
            "self.domain%s = args['domain%s']\n" % (self.id, self.id)
        self.local_declarations = "cdef RealNumber rn%s\n" % self.id

    def cython_init(self, loc):
        r"""
        Generates code to initialize an mpfr_t reference (a variable, an
        array reference, etc.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.cython_init('foo[i]')
            'mpfr_init2(foo[i], self.domain.prec())'
        """
        return je("mpfr_init2({{ loc }}, self.domain{{ myself.id }}.prec())",
                  myself=self, loc=loc)

    def cython_clear(self, loc):
        r"""
        Generates code to clear an mpfr_t reference (a variable, an
        array reference, etc.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.cython_clear('foo[i]')
            'mpfr_clear(foo[i])'
        """
        return 'mpfr_clear(%s)' % loc

    def assign_c_from_py(self, c, py):
        r"""
        Given a Cython variable/array reference/etc. of this storage type,
        and a Python expression, generate code to assign to the Cython
        variable from the Python expression.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpfr.assign_c_from_py('foo[i]', 'bar[j]')
            'rn = self.domain(bar[j])\nmpfr_set(foo[i], rn.value, MPFR_RNDN)'
        """
        return je(ri(0, """
            rn{{ myself.id }} = self.domain({{ py }})
            mpfr_set({{ c }}, rn.value, MPFR_RNDN)"""),
            myself=self, c=c, py=py)

ty_mpfr = StorageTypeMPFR()

class StorageTypeMPC(StorageTypeAutoReference):
    r"""
    StorageTypeMPC is a subtype of StorageTypeAutoReference that deals
    the MPC's mpc_t type.

    For any given program that we're interpreting, ty_mpc can only
    refer to a single precision.  An interpreter that needs to use
    two precisions of mpc_t in the same program should instantiate two
    separate instances of StorageTypeMPC.  (Interpreters that need
    to handle arbitrarily many precisions in the same program are not
    handled at all.)
    """

    def __init__(self, id=''):
        r"""
        Initializes the id property, as well as the properties described
        in the documentation for StorageTypeAutoReference.__init__.

        The id property is used if you want to have an interpreter
        that handles two instances of StorageTypeMPC (that is,
        handles mpC_t variables at two different precisions
        simultaneously).  It's a string that's used to generate
        variable names that don't conflict.  (The id system has
        never actually been used, so bugs probably remain.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpc.class_member_declarations
            'cdef object domain\ncdef ComplexNumber domain_element\n'
            sage: ty_mpc.class_member_initializations
            "self.domain = args['domain']\nself.domain_element = self.domain.zero()\n"
            sage: ty_mpc.local_declarations
            'cdef ComplexNumber cn\n'
            sage: ty_mpc.decl_type
            'mpc_t'
            sage: ty_mpc.ref_type
            'mpc_ptr'

        TESTS::

            sage: ty_mpfr2 = StorageTypeMPC(id='_the_second')
            sage: ty_mpfr2.class_member_declarations
            'cdef object domain_the_second\ncdef ComplexNumber domain_element_the_second\n'
            sage: ty_mpfr2.class_member_initializations
            "self.domain_the_second = args['domain_the_second']\nself.domain_element_the_second = self.domain.zero()\n"
            sage: ty_mpfr2.local_declarations
            'cdef ComplexNumber cn_the_second\n'

        """
        StorageTypeAutoReference.__init__(self, 'mpc_t', 'mpc_ptr')
        self.id = id
        self.class_member_declarations = "cdef object domain%s\ncdef ComplexNumber domain_element%s\n" % (self.id, self.id)
        self.class_member_initializations = \
            "self.domain%s = args['domain%s']\nself.domain_element%s = self.domain.zero()\n" % (self.id, self.id, self.id)
        self.local_declarations = "cdef ComplexNumber cn%s\n" % self.id

    def cython_init(self, loc):
        r"""
        Generates code to initialize an mpc_t reference (a variable, an
        array reference, etc.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpc.cython_init('foo[i]')
            'mpc_init2(foo[i], self.domain_element._prec)'
        """
        return je("mpc_init2({{ loc }}, self.domain_element{{ myself.id }}._prec)",
                  myself=self, loc=loc)

    def cython_clear(self, loc):
        r"""
        Generates code to clear an mpfr_t reference (a variable, an
        array reference, etc.)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpc.cython_clear('foo[i]')
            'mpc_clear(foo[i])'
        """
        return 'mpc_clear(%s)' % loc

    def assign_c_from_py(self, c, py):
        r"""
        Given a Cython variable/array reference/etc. of this storage type,
        and a Python expression, generate code to assign to the Cython
        variable from the Python expression.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_mpc.assign_c_from_py('foo[i]', 'bar[j]')
            'cn = self.domain(bar[j])\nmpc_set_fr_fr(foo[i], cn.__re, cn.__im, MPC_RNDNN)'
        """
        return je("""
cn{{ myself.id }} = self.domain({{ py }})
mpc_set_fr_fr({{ c }}, cn.__re, cn.__im, MPC_RNDNN)""", myself=self, c=c, py=py)

ty_mpc = StorageTypeMPC()
