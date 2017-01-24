r"""
Generate interpreters for fast_callable

AUTHORS:

- Carl Witty

This file is part of the Sage support for "planned" computations;
that is, computations that are separated into a planning stage and
a plan-execution stage.  Here, we generate fast interpreters for plan
executions.

There are at least two kinds of computations that are often planned in
this fashion.  First is arithmetic expression evaluation, where we
take an arbitrary user-specified arithmetic expression and compile it
into a bytecode form for fast interpretation.  Second is things like
FFTs and large multiplications, where large problems are split into
multiple smaller problems... we can do the logical "splitting" for a
given size only once, producing a plan which can be reused as often as
we want for different problems of the same size.  Currently only
arithmetic expression evaluation is implemented, but other kinds of
planned computations should be easy to add.

Typically, for arithmetic expressions, we want the storage of
intermediate results to be handled automatically (on a stack); for
FFTs/multiplications/etc., the planner will keep track of intermediate
results itself.

For arithmetic expression evaluation, we want to have lots of
interpreters (at least one, and possibly several, per
specially-handled type).  Also, for any given type, we have many
possible variants of instruction encoding, etc.; some of these could
be handled with conditional compilation, but some are more
complicated.  So we end up writing an interpreter generator.

We want to share as much code as possible across all of these
interpreters, while still maintaining the freedom to make drastic
changes in the interpretation strategy (which may change the
generated code, the calling convention for the interpreter, etc.)

To make this work, the interpreter back-end is divided into three
parts:

1. The interpreter itself, in C or C++.

2. The wrapper, which is a Cython object holding the
   constants, code, etc., and which actually calls the interpreter.

3. The code generator.

We generate parts 1 and 2.  The code generator is table-driven,
and we generate the tables for the code generator.

There are a lot of techniques for fast interpreters that we do not yet
use; hopefully at least some of these will eventually be implemented:

- using gcc's "labels as values" extension where available

- top-of-stack caching

- superinstructions and/or superoperators

- static stack caching

- context threading/subrouting threading

- selective inlining/dynamic superinstructions

- automatic replication

Interpreters may be stack-based or register-based.  Recent research
suggests that register-based interpreters are better, but the
researchers are investigating interpreters for entire programming
languages, rather than interpreters for expressions.  I suspect
that stack-based expression interpreters may be better.  However,
we'll implement both varieties and see what's best.

The relative costs of stack- and register-based interpreters will
depend on the costs of moving values.  For complicated types (like
mpz_t), a register-based interpreter will quite likely be better,
since it will avoid moving values.

We will NOT support any sort of storage of bytecode; instead, the
code must be re-generated from expression trees in every Sage run.
This means that we can trivially experiment with different styles of
interpreter, or even use quite different interpreters depending on
the architecture, without having to worry about forward and backward
compatibility.
"""

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

from __future__ import print_function

import os
import re
import six
from jinja2 import Environment
from jinja2.runtime import StrictUndefined
from collections import defaultdict
from distutils.extension import Extension

from os.path import getmtime

##############################
# This module is used during the Sage build process, so it should not
# use any other Sage modules.  (In particular, it MUST NOT use any
# Cython modules -- they won't be built yet!)
# Also, we have some trivial dependency tracking, where we don't
# rebuild the interpreters if this file hasn't changed; if
# interpreter configuration is split out into a separate file,
# that will have to be changed.
##############################


# We share a single jinja2 environment among all templating in this
# file.  We use trim_blocks=True (which means that we ignore white
# space after "%}" jinja2 command endings), and set undefined to
# complain if we use an undefined variable.
jinja_env = Environment(trim_blocks=True, undefined=StrictUndefined)

# Allow 'i' as a shorter alias for the built-in 'indent' filter.
jinja_env.filters['i'] = jinja_env.filters['indent']

autogen_warn = "Automatically generated by {}.  Do not edit!".format(__file__)

def indent_lines(n, text):
    r"""
    INPUTS:

    - n -- indentation amount
    - text -- text to indent

    Indents each line in text by n spaces.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import indent_lines
        sage: indent_lines(3, "foo")
        '   foo'
        sage: indent_lines(3, "foo\nbar")
        '   foo\n   bar'
        sage: indent_lines(3, "foo\nbar\n")
        '   foo\n   bar\n'
    """
    lines = text.splitlines(True)
    spaces = ' ' * n
    return ''.join(spaces + line for line in lines)

def je(template, **kwargs):
    r"""
    A convenience method for creating strings with Jinja templates.
    The name je stands for "Jinja evaluate".

    The first argument is the template string; remaining keyword
    arguments define Jinja variables.

    If the first character in the template string is a newline, it is
    removed (this feature is useful when using multi-line templates defined
    with triple-quoted strings -- the first line doesn't have to be on
    the same line as the quotes, which would screw up the indentation).

    (This is very inefficient, because it recompiles the Jinja
    template on each call; don't use it in situations where
    performance is important.)

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import je
        sage: je("{{ a }} > {{ b }} * {{ c }}", a='"a suffusion of yellow"', b=3, c=7)
        u'"a suffusion of yellow" > 3 * 7'
    """
    if len(template) > 0 and template[0] == '\n':
        template = template[1:]

    # It looks like Jinja2 automatically removes one trailing newline?
    if len(template) > 0 and template[-1] == '\n':
        template = template + '\n'

    tmpl = jinja_env.from_string(template)
    return tmpl.render(kwargs)

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
        primitive types, and for types like PyObject* (where you're only
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
            u'foo = bar'
            sage: ty_python.assign_c_from_py('foo[i]', 'bar[j]')
            u'foo[i] = <PyObject *>bar[j]; Py_INCREF(foo[i])'
            sage: ty_mpfr.assign_c_from_py('foo', 'bar')
            u'rn = self.domain(bar)\nmpfr_set(foo, rn.value, MPFR_RNDN)'
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
            u'    cdef int _n_args\n    cdef mpfr_t* _args\n'
        """
        return je("""
{# XXX Variables here (and everywhere, really) should actually be Py_ssize_t #}
    cdef int _n_{{ name }}
    cdef {{ myself.cython_array_type() }} _{{ name }}
""", myself=self, name=name)

    def alloc_chunk_data(self, name, len):
        r"""
        Return a string allocating the memory for the class members for
        a memory chunk with this storage type and the given name.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: print(ty_mpfr.alloc_chunk_data('args', 'MY_LENGTH'))
                    self._n_args = MY_LENGTH
                    self._args = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * MY_LENGTH)
                    if self._args == NULL: raise MemoryError
                    for i in range(MY_LENGTH):
                        mpfr_init2(self._args[i], self.domain.prec())
            <BLANKLINE>
        """
        return je("""
        self._n_{{ name }} = {{ len }}
        self._{{ name }} = <{{ myself.c_ptr_type() }}>sig_malloc(sizeof({{ myself.c_decl_type() }}) * {{ len }})
        if self._{{ name }} == NULL: raise MemoryError
{% if myself.needs_cython_init_clear() %}
        for i in range({{ len }}):
            {{ myself.cython_init('self._%s[i]' % name) }}
{% endif %}
""", myself=self, name=name, len=len)

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
        return je("""
        if self._{{ name }}:
{%     if myself.needs_cython_init_clear() %}
            for i in range(self._n_{{ name }}):
                {{ myself.cython_clear('self._%s[i]' % name) }}
{%     endif %}
            sig_free(self._{{ name }})
""", myself=self, name=name)

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
        u'z_c = CDE_to_dz(z_py)'
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
        StorageTypeAssignable.__init__(self, 'PyObject*')

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
            u'    cdef object _list_args\n    cdef int _n_args\n    cdef PyObject** _args\n'
        """
        return je("""
    cdef object _list_{{ name }}
    cdef int _n_{{ name }}
    cdef {{ myself.cython_array_type() }} _{{ name }}
""", myself=self, name=name)

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
        return je("""
        self._n_{{ name }} = {{ len }}
        self._list_{{ name }} = PyList_New(self._n_{{ name }})
        self._{{ name }} = (<PyListObject *>self._list_{{ name }}).ob_item
""", myself=self, name=name, len=len)

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
            u'foo[i] = <PyObject *>bar[j]; Py_INCREF(foo[i])'
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
            u'foo[i] = NULL'
        """
        return je("{{ loc }} = NULL", loc=loc)

    def cython_clear(self, loc):
        r"""
        Generates code to clear a variable (or array reference) holding
        a PyObject*.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: ty_python.cython_clear('foo[i]')
            u'Py_CLEAR(foo[i])'
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
        StorageType.__init__(self)
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
        StorageTypeAutoReference.__init__(self, 'mpfr_t', 'mpfr_ptr')
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
            u'mpfr_init2(foo[i], self.domain.prec())'
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
            u'rn = self.domain(bar[j])\nmpfr_set(foo[i], rn.value, MPFR_RNDN)'
        """
        return je("""
rn{{ myself.id }} = self.domain({{ py }})
mpfr_set({{ c }}, rn.value, MPFR_RNDN)""", myself=self, c=c, py=py)

ty_mpfr = StorageTypeMPFR()

class MemoryChunk(object):
    r"""
    Memory chunks control allocation, deallocation, initialization,
    etc.  of the vectors and objects in the interpreter.  Basically,
    there is one memory chunk per argument to the C interpreter.

    There are three "generic" varieties of memory chunk: "constants",
    "arguments", and "scratch".  These are named after their most
    common use, but they could be used for other things in some
    interpreters.

    All three kinds of chunks are allocated in the wrapper class.
    Constants are initialized when the wrapper is constructed;
    arguments are initialized in the __call__ method, from the
    caller's arguments.  "scratch" chunks are not initialized at all;
    they are used for scratch storage (often, but not necessarily, for
    a stack) in the interpreter.

    Interpreters which need memory chunks that don't fit into these
    categories can create new subclasses of MemoryChunk.
    """

    def __init__(self, name, storage_type):
        r"""
        Initialize an instance of MemoryChunk.

        This sets the properties "name" (the name of this memory chunk;
        used in generated variable names, etc.) and "storage_type",
        which is a StorageType object.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: mc.name
            'args'
            sage: mc.storage_type is ty_mpfr
            True
        """
        self.name = name
        self.storage_type = storage_type

    def __repr__(self):
        r"""
        Give a string representation of this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: mc
            {MC:args}
            sage: mc.__repr__()
            '{MC:args}'
        """
        return '{MC:%s}' % self.name

    def declare_class_members(self):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: mc.declare_class_members()
            u'    cdef int _n_args\n    cdef mpfr_t* _args\n'
        """
        return self.storage_type.declare_chunk_class_members(self.name)

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: print(mc.init_class_members())
                    count = args['args']
                    self._n_args = count
                    self._args = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * count)
                    if self._args == NULL: raise MemoryError
                    for i in range(count):
                        mpfr_init2(self._args[i], self.domain.prec())
            <BLANKLINE>
        """
        return ""

    def dealloc_class_members(self):
        r"""
        Return a string to be put in the __dealloc__ method of a wrapper
        class using this memory chunk, to deallocate the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: print(mc.dealloc_class_members())
                    if self._args:
                        for i in range(self._n_args):
                            mpfr_clear(self._args[i])
                        sig_free(self._args)
            <BLANKLINE>
        """
        return ""

    def declare_parameter(self):
        r"""
        Return the string to use to declare the interpreter parameter
        corresponding to this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: mc.declare_parameter()
            'mpfr_t* args'
        """
        return '%s %s' % (self.storage_type.c_ptr_type(), self.name)

    def declare_call_locals(self):
        r"""
        Return a string to put in the __call__ method of a wrapper
        class using this memory chunk, to allocate local variables.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.declare_call_locals()
            u'        cdef RealNumber retval = (self.domain)()\n'
        """
        return ""

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkConstants('constants', ty_mpfr)
            sage: mc.pass_argument()
            'self._constants'
        """
        raise NotImplementedError

    def pass_call_c_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter, for use in the call_c method.
        Almost always the same as pass_argument.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkConstants('constants', ty_mpfr)
            sage: mc.pass_call_c_argument()
            'self._constants'
        """
        return self.pass_argument()

    def needs_cleanup_on_error(self):
        r"""
        In an interpreter that can terminate prematurely (due to an
        exception from calling Python code, or divide by zero, or
        whatever) it will just return at the end of the current instruction,
        skipping the rest of the program.  Thus, it may still have
        values pushed on the stack, etc.

        This method returns True if this memory chunk is modified by the
        interpreter and needs some sort of cleanup when an error happens.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkConstants('constants', ty_mpfr)
            sage: mc.needs_cleanup_on_error()
            False
        """
        return False

    def is_stack(self):
        r"""
        Says whether this memory chunk is a stack.  This affects code
        generation for instructions using this memory chunk.

        It would be nicer to make this object-oriented somehow, so
        that the code generator called MemoryChunk methods instead of
        using::

            if ch.is_stack():
                ... hardcoded stack code
            else:
                ... hardcoded non-stack code

        but that hasn't been done yet.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkScratch('scratch', ty_mpfr)
            sage: mc.is_stack()
            False
            sage: mc = MemoryChunkScratch('stack', ty_mpfr, is_stack=True)
            sage: mc.is_stack()
            True
        """
        return False

    def is_python_refcounted_stack(self):
        r"""
        Says whether this memory chunk refers to a stack where the entries
        need to be INCREF/DECREF'ed.

        It would be nice to make this object-oriented, so that the
        code generator called MemoryChunk methods to do the potential
        INCREF/DECREF and didn't have to explicitly test
        is_python_refcounted_stack.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkScratch('args', ty_python)
            sage: mc.is_python_refcounted_stack()
            False
            sage: mc = MemoryChunkScratch('args', ty_python, is_stack=True)
            sage: mc.is_python_refcounted_stack()
            True
            sage: mc = MemoryChunkScratch('args', ty_mpfr, is_stack=True)
            sage: mc.is_python_refcounted_stack()
            False
        """
        return self.is_stack() and self.storage_type.python_refcounted()

class MemoryChunkLonglivedArray(MemoryChunk):
    r"""
    MemoryChunkLonglivedArray is a subtype of MemoryChunk that deals
    with memory chunks that are both 1) allocated as class members (rather
    than being allocated in __call__) and 2) are arrays.
    """

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_double)
            sage: print(mc.init_class_members())
                    count = args['args']
                    self._n_args = count
                    self._args = <double*>sig_malloc(sizeof(double) * count)
                    if self._args == NULL: raise MemoryError
            <BLANKLINE>
        """
        return je("""
        count = args['{{ myself.name }}']
{% print(myself.storage_type.alloc_chunk_data(myself.name, 'count')) %}
""", myself=self)

    def dealloc_class_members(self):
        r"""
        Return a string to be put in the __dealloc__ method of a wrapper
        class using this memory chunk, to deallocate the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: print(mc.dealloc_class_members())
                    if self._args:
                        for i in range(self._n_args):
                            mpfr_clear(self._args[i])
                        sig_free(self._args)
            <BLANKLINE>
        """
        return self.storage_type.dealloc_chunk_data(self.name)

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkConstants('constants', ty_mpfr)
            sage: mc.pass_argument()
            'self._constants'
        """
        return 'self._%s' % self.name

class MemoryChunkConstants(MemoryChunkLonglivedArray):
    r"""
    MemoryChunkConstants is a subtype of MemoryChunkLonglivedArray.

    MemoryChunkConstants chunks have their contents set in the
    wrapper's __init__ method (and not changed afterward).
    """

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkConstants('constants', ty_mpfr)
            sage: print(mc.init_class_members())
                    val = args['constants']
                    self._n_constants = len(val)
                    self._constants = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * len(val))
                    if self._constants == NULL: raise MemoryError
                    for i in range(len(val)):
                        mpfr_init2(self._constants[i], self.domain.prec())
                    for i in range(len(val)):
                        rn = self.domain(val[i])
                        mpfr_set(self._constants[i], rn.value, MPFR_RNDN)
            <BLANKLINE>
        """
        return je("""
        val = args['{{ myself.name }}']
{% print(myself.storage_type.alloc_chunk_data(myself.name, 'len(val)')) %}
        for i in range(len(val)):
            {{ myself.storage_type.assign_c_from_py('self._%s[i]' % myself.name, 'val[i]') | i(12) }}
""", myself=self)

class MemoryChunkArguments(MemoryChunkLonglivedArray):
    r"""
    MemoryChunkArguments is a subtype of MemoryChunkLonglivedArray,
    for dealing with arguments to the wrapper's ``__call__`` method.

    Currently the ``__call__`` method is declared to take a varargs
    `*args` argument tuple.  We assume that the MemoryChunk named `args`
    will deal with that tuple.
    """

    def setup_args(self):
        r"""
        Handle the arguments of __call__ -- copy them into a pre-allocated
        array, ready to pass to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: print(mc.setup_args())
            cdef mpfr_t* c_args = self._args
            cdef int i
            for i from 0 <= i < len(args):
                rn = self.domain(args[i])
                mpfr_set(self._args[i], rn.value, MPFR_RNDN)
            <BLANKLINE>
        """
        return je("""
cdef {{ myself.storage_type.c_ptr_type() }} c_args = self._args
cdef int i
for i from 0 <= i < len(args):
    {{ myself.storage_type.assign_c_from_py('self._args[i]', 'args[i]') | i(4) }}
""", myself=self)

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: mc.pass_argument()
            'c_args'
        """
        return 'c_args'

class MemoryChunkScratch(MemoryChunkLonglivedArray):
    r"""
    MemoryChunkScratch is a subtype of MemoryChunkLonglivedArray
    for dealing with memory chunks that are allocated in the wrapper,
    but only used in the interpreter -- stacks, scratch registers, etc.

    (Currently these are only used as stacks.)
    """

    def __init__(self, name, storage_type, is_stack=False):
        r"""
        Initialize an instance of MemoryChunkScratch.

        Initializes the _is_stack property, as well as
        the properties described in the documentation for
        MemoryChunk.__init__.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkScratch('stack', ty_double, is_stack=True)
            sage: mc.name
            'stack'
            sage: mc.storage_type is ty_double
            True
            sage: mc._is_stack
            True
        """
        MemoryChunkLonglivedArray.__init__(self, name, storage_type)
        self._is_stack = is_stack

    def is_stack(self):
        r"""
        Says whether this memory chunk is a stack.  This affects code
        generation for instructions using this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkScratch('stack', ty_mpfr, is_stack=True)
            sage: mc.is_stack()
            True
        """
        return self._is_stack

    def needs_cleanup_on_error(self):
        r"""
        In an interpreter that can terminate prematurely (due to an
        exception from calling Python code, or divide by zero, or
        whatever) it will just return at the end of the current instruction,
        skipping the rest of the program.  Thus, it may still have
        values pushed on the stack, etc.

        This method returns True if this memory chunk is modified by the
        interpreter and needs some sort of cleanup when an error happens.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkScratch('registers', ty_python)
            sage: mc.needs_cleanup_on_error()
            True
        """
        return self.storage_type.python_refcounted()

    def handle_cleanup(self):
        r"""
        Handle the cleanup if the interpreter exits with an error.

        For scratch/stack chunks that hold Python-refcounted values,
        we assume that they are filled with NULL on every entry to the
        interpreter.  If the interpreter exited with an error, it may
        have left values in the chunk, so we need to go through
        the chunk and Py_CLEAR it.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkScratch('registers', ty_python)
            sage: print(mc.handle_cleanup())
            for i in range(self._n_registers):
                Py_CLEAR(self._registers[i])
            <BLANKLINE>
        """
        # XXX This is a lot slower than it needs to be, because
        # we don't have a "cdef int i" in scope here.
        return je("""
for i in range(self._n_{{ myself.name }}):
    Py_CLEAR(self._{{ myself.name }}[i])
""", myself=self)

class MemoryChunkRRRetval(MemoryChunk):
    r"""
    A special-purpose memory chunk, for dealing with the return value
    of the RR-based interpreter.
    """

    def declare_class_members(self):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.declare_class_members()
            ''
        """
        return ""

    def declare_call_locals(self):
        r"""
        Return a string to put in the __call__ method of a wrapper
        class using this memory chunk, to allocate local variables.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.declare_call_locals()
            u'        cdef RealNumber retval = (self.domain)()\n'
        """
        return je("""
        cdef RealNumber {{ myself.name }} = (self.domain)()
""", myself=self)

    def declare_parameter(self):
        r"""
        Return the string to use to declare the interpreter parameter
        corresponding to this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.declare_parameter()
            'mpfr_t retval'
        """
        return '%s %s' % (self.storage_type.c_reference_type(), self.name)

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.pass_argument()
            u'retval.value'
        """
        return je("""{{ myself.name }}.value""", myself=self)

    def pass_call_c_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter, for use in the call_c method.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.pass_call_c_argument()
            'result'
        """
        return "result"

class MemoryChunkPythonArguments(MemoryChunk):
    r"""
    A special-purpose memory chunk, for the generic Python-object based
    interpreter.  Rather than copy the arguments into an array allocated
    in the wrapper, we use the PyTupleObject internals and pass the array
    that's inside the argument tuple.
    """

    def declare_class_members(self):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
        """
        return "    cdef int _n_%s\n" % self.name

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
            sage: mc.init_class_members()
            u"        count = args['args']\n        self._n_args = count\n"
        """
        return je("""
        count = args['{{ myself.name }}']
        self._n_args = count
""", myself=self)

    def setup_args(self):
        r"""
        Handle the arguments of __call__.  Nothing to do.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
            sage: mc.setup_args()
            ''
        """
        return ''

    def pass_argument(self):
        r"""
        Pass the innards of the argument tuple to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
            sage: mc.pass_argument()
            '(<PyTupleObject*>args).ob_item'
        """
        return "(<PyTupleObject*>args).ob_item"

class MemoryChunkElementArguments(MemoryChunkPythonArguments):
    r"""
    A special-purpose memory chunk, for the Python-object based
    interpreters that want to process (and perhaps modify) the data.

    We allocate a new list (via the map function) on every call to
    hold the modified arguments.  That's not strictly necessary --
    we could pre-allocate a list and map into it -- but this lets us
    use simpler code for a very-likely-negligible efficiency cost.
    (The Element interpreter is going to allocate lots of objects
    as it runs, anyway.)
    """

    def setup_args(self):
        r"""
        Handle the arguments of __call__.  Note: This hardcodes
        "self._domain".

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkElementArguments('args', ty_python)
            sage: mc.setup_args()
            'mapped_args = map(self._domain, args)\n'
        """
        return "mapped_args = map(self._domain, args)\n"

    def pass_argument(self):
        r"""
        Pass the innards of the argument tuple to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkElementArguments('args', ty_python)
            sage: mc.pass_argument()
            '(<PyListObject*>mapped_args).ob_item'
        """
        return "(<PyListObject*>mapped_args).ob_item"

class MemoryChunkPyConstant(MemoryChunk):
    r"""
    A special-purpose memory chunk, for holding a single Python constant
    and passing it to the interpreter as a PyObject*.
    """

    def __init__(self, name):
        r"""
        Initialize an instance of MemoryChunkPyConstant.

        Always uses the type ty_python.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.name
            'domain'
            sage: mc.storage_type is ty_python
            True
        """
        MemoryChunk.__init__(self, name, ty_python)

    def declare_class_members(self):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.declare_class_members()
            u'    cdef object _domain\n'
        """
        return je("""
    cdef object _{{ myself.name }}
""", myself=self)

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.init_class_members()
            u"        self._domain = args['domain']\n"
        """
        return je("""
        self._{{ myself.name }} = args['{{ myself.name }}']
""", myself=self)

    def declare_parameter(self):
        r"""
        Return the string to use to declare the interpreter parameter
        corresponding to this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.declare_parameter()
            'PyObject* domain'
        """
        return 'PyObject* %s' % self.name

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.pass_argument()
            '<PyObject*>self._domain'
        """
        return '<PyObject*>self._%s' % self.name

def params_gen(**chunks):
    r"""
    Instructions have a parameter specification that says where they get
    their inputs and where their outputs go.  Each parameter has
    the same form: it is a triple (chunk, addr, len).  The chunk says
    where the parameter is read from/written to.  The addr says which
    value in the chunk is used.  If the chunk is a stack chunk, then
    addr must be null; the value will be read from/written to the top
    of the stack.  Otherwise, addr must be an integer, or another chunk;
    if addr is another chunk, then the next value is read from that chunk
    to be the address.

    The len says how many values to read/write.  It can be either None
    (meaning to read/write only a single value), an integer, or
    another chunk; if it is a chunk, then the next value is read from that
    chunk to be the len.  Note that specifying len changes the types
    given to the instruction, so len==None is different than len==1 even
    though both mean to use a single value.

    These parameter specifications are cumbersome to write by hand, so
    there's also a simple string format for them.  This (curried)
    function parses the simple string format and produces parameter
    specifications.  The params_gen function takes keyword arguments
    mapping single-character names to memory chunks.  The string format
    uses these names.  The params_gen function returns another function,
    that takes two strings and returns a pair of lists of parameter
    specifications.

    Each string is the concatenation of arbitrarily many specifications.
    Each specification consists of an address and a length.  The
    address is either a single character naming a stack chunk,
    or a string of the form 'A[B]' where A names a non-stack chunk
    and B names the code chunk.  The length is either empty, or '@n'
    for a number n (meaning to use that many arguments), or '@C', where
    C is the code chunk.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: mc_stack = MemoryChunkScratch('stack', ty_double, is_stack=True)
        sage: mc_args = MemoryChunkArguments('args', ty_double)
        sage: mc_code = MemoryChunkConstants('code', ty_int)

        sage: pg = params_gen(D=mc_code, A=mc_args, S=mc_stack)
        sage: pg('S', '')
        ([({MC:stack}, None, None)], [])
        sage: pg('A[D]', '')
        ([({MC:args}, {MC:code}, None)], [])
        sage: pg('S@5', '')
        ([({MC:stack}, None, 5)], [])
        sage: pg('S@D', '')
        ([({MC:stack}, None, {MC:code})], [])
        sage: pg('A[D]@D', '')
        ([({MC:args}, {MC:code}, {MC:code})], [])
        sage: pg('SSS@D', 'A[D]S@D')
        ([({MC:stack}, None, None), ({MC:stack}, None, None), ({MC:stack}, None, {MC:code})], [({MC:args}, {MC:code}, None), ({MC:stack}, None, {MC:code})])
    """

    def make_params(s):
        p = []
        s = s.strip()
        while s:
            chunk_code = s[0]
            s = s[1:]
            chunk = chunks[chunk_code]
            addr = None
            ch_len = None
            # shouldn't hardcode 'code' here
            if chunk.is_stack() or chunk.name == 'code':
                pass
            else:
                m = re.match(r'\[(?:([0-9]+)|([a-zA-Z]))\]', s)
                if m.group(1):
                    addr = int(m.group(1))
                else:
                    ch = chunks[m.group(2)]
                    assert ch.storage_type is ty_int
                    addr = ch
                s = s[m.end():].strip()
            if len(s) and s[0] == '@':
                m = re.match(r'@(?:([0-9]+)|([a-zA-Z]))', s)
                if m.group(1):
                    ch_len = int(m.group(1))
                else:
                    ch = chunks[m.group(2)]
                    assert ch.storage_type is ty_int
                    ch_len = ch
                s = s[m.end():].strip()
            p.append((chunk, addr, ch_len))
        return p

    def params(s_ins, s_outs):
        ins = make_params(s_ins)
        outs = make_params(s_outs)
        return (ins, outs)

    return params

def string_of_addr(a):
    r"""
    An address or a length from a parameter specification may be
    either None, an integer, or a MemoryChunk.  If the address or
    length is an integer or a MemoryChunk, this function will convert
    it to a string giving an expression that will evaluate to the correct
    address or length.  (See the docstring for params_gen for more
    information on parameter specifications.)

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: mc_code = MemoryChunkConstants('code', ty_int)
        sage: string_of_addr(mc_code)
        '*code++'
        sage: string_of_addr(42r)
        '42'
    """
    if isinstance(a, six.integer_types):
        return str(a)
    assert(isinstance(a, MemoryChunk))
    return '*%s++' % a.name

class InstrSpec(object):
    r"""
    Each instruction in an interpreter is represented as an InstrSpec.
    This contains all the information that we need to generate code
    to interpret the instruction; it also is used to build the tables
    that fast_callable uses, so this is the nexus point between
    users of the interpreter (possibly pure Python) and the
    generated C interpreter.

    The underlying instructions are matched to the caller by name.
    For instance, fast_callable assumes that if the interpreter has an
    instruction named 'cos', then it will take a single argument,
    return a single result, and implement the cos() function.

    The print representation of an instruction (which will probably
    only be used when doctesting this file) consists of the name,
    a simplified stack effect, and the code (truncated if it's long).
    The stack effect has two parts, the input and the output, separated
    by '->'; the input shows what will be popped from the stack,
    the output what will be placed on the stack.  Each consists of
    a sequence of 'S' and '*' characters, where 'S' refers to a single
    argument and '*' refers to a variable number of arguments.

    The code for an instruction is a small snippet of C code.  It has
    available variables 'i0', 'i1', ..., 'o0', 'o1', ...; one variable
    for each input and output; its job is to assign values to the output
    variables, based on the values of the input variables.

    Normally, in an interpreter that uses doubles, each of the input
    and output variables will be a double.  If i0 actually represents
    a variable number of arguments, then it will be a pointer to
    double instead, and there will be another variable n_i0 giving
    the actual number of arguments.

    When instructions refer to auto-reference types, they actually
    get a pointer to the data in its original location; it is
    not copied into a local variable.  Mostly, this makes no difference,
    but there is one potential problem to be aware of.  It is possible
    for an output variable to point to the same object as an input
    variable; in fact, this usually will happen when you're working
    with the stack.  If the instruction maps to a single function call,
    then this is fine; the standard auto-reference implementations
    (GMP, MPFR, etc.) are careful to allow having the input and output
    be the same.  But if the instruction maps to multiple function
    calls, you may need to use a temporary variable.

    Here's an example of this issue.  Suppose you want to make an
    instruction that does ``out = a+b*c``.  You write code like this::

        out = b*c
        out = a+out

    But out will actually share the same storage as a; so the first line
    modifies a, and you actually end up computing 2*(b+c).  The fix
    is to only write to the output once, at the very end of your
    instruction.

    Instructions are also allowed to access memory chunks (other than
    the stack and code) directly.  They are available as C variables
    with the same name as the chunk.  This is useful if some type of
    memory chunk doesn't fit well with the params_gen interface.

    There are additional reference-counting rules that must be
    followed if your interpreter operates on Python objects; these
    rules are described in the docstring of the PythonInterpreter
    class.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: pg = RDFInterpreter().pg
        sage: InstrSpec('add', pg('SS','S'), code='o0 = i0+i1;')
        add: SS->S = 'o0 = i0+i1;'
    """

    def __init__(self, name, io, code=None, uses_error_handler=False, handles_own_decref=False):
        r"""
        Initialize an InstrSpec.

        INPUTS:
            name -- the name of the instruction
            io -- a pair of lists of parameter specifications for I/O of the
                  instruction
            code -- a string containing a snippet of C code to read
                    from the input variables and write to the output variables
            uses_error_handler -- True if the instruction calls Python
                                  and jumps to error: on a Python error
            handles_own_decref -- True if the instruction handles Python
                                  objects and includes its own
                                  reference-counting

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *

            sage: pg = RDFInterpreter().pg
            sage: InstrSpec('add', pg('SS','S'), code='o0 = i0+i1;')
            add: SS->S = 'o0 = i0+i1;'
            sage: instr = InstrSpec('py_call', pg('P[D]S@D', 'S'), code=('This is very complicated.  ' + 'blah ' * 30)); instr
            py_call: *->S = 'This is very compli... blah blah blah '
            sage: instr.name
            'py_call'
            sage: instr.inputs
            [({MC:py_constants}, {MC:code}, None), ({MC:stack}, None, {MC:code})]
            sage: instr.outputs
            [({MC:stack}, None, None)]
            sage: instr.code
            'This is very complicated.  blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah '
            sage: instr.parameters
            ['py_constants', 'n_inputs']
            sage: instr.n_inputs
            0
            sage: instr.n_outputs
            1
        """
        self.name = name
        self.inputs = io[0]
        self.outputs = io[1]
        self.uses_error_handler = uses_error_handler
        self.handles_own_decref = handles_own_decref
        if code is not None:
            self.code = code
        # XXX We assume that there is only one stack
        n_inputs = 0
        n_outputs = 0
        in_effect = ''
        out_effect = ''
        p = []
        for (ch, addr, len) in self.inputs:
            if ch.is_stack():
                if len is None:
                    n_inputs += 1
                    in_effect += 'S'
                elif isinstance(len, six.integer_types):
                    n_inputs += len
                    in_effect += 'S%d' % len
                else:
                    p.append('n_inputs')
                    in_effect += '*'
            else:
                p.append(ch.name)
        for (ch, addr, len) in self.outputs:
            if ch.is_stack():
                if len is None:
                    n_outputs += 1
                    out_effect += 'S'
                elif isinstance(len, six.integer_types):
                    n_outputs += len
                    out_effect += 'S%d' % len
                else:
                    p.append('n_outputs')
                    out_effect += '*'
            else:
                p.append(ch.name)
        self.parameters = p
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs
        self.in_effect = in_effect
        self.out_effect = out_effect

    def __repr__(self):
        r"""
        Produce a string representing a given instruction, consisting
        of its name, a brief stack specification, and its code
        (possibly abbreviated).

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: pg = RDFInterpreter().pg
            sage: InstrSpec('add', pg('SS','S'), code='o0 = i0+i1;')
            add: SS->S = 'o0 = i0+i1;'
        """
        rcode = repr(self.code)
        if len(rcode) > 40:
            rcode = rcode[:20] + '...' + rcode[-17:]
        return '%s: %s->%s = %s' % \
            (self.name, self.in_effect, self.out_effect, rcode)

# Now we have a series of helper functions that make it slightly easier
# to create instructions.

def instr_infix(name, io, op):
    r"""
    A helper function for creating instructions implemented by
    a single infix binary operator.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: pg = RDFInterpreter().pg
        sage: instr_infix('mul', pg('SS', 'S'), '*')
        mul: SS->S = 'o0 = i0 * i1;'
    """
    return InstrSpec(name, io, code='o0 = i0 %s i1;' % op)

def instr_funcall_2args(name, io, op):
    r"""
    A helper function for creating instructions implemented by
    a two-argument function call.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: pg = RDFInterpreter().pg
        sage: instr_funcall_2args('atan2', pg('SS', 'S'), 'atan2')
        atan2: SS->S = 'o0 = atan2(i0, i1);'
    """
    return InstrSpec(name, io, code='o0 = %s(i0, i1);' % op)

def instr_unary(name, io, op):
    r"""
    A helper function for creating instructions with one input
    and one output.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: pg = RDFInterpreter().pg
        sage: instr_unary('sin', pg('S','S'), 'sin(i0)')
        sin: S->S = 'o0 = sin(i0);'
        sage: instr_unary('neg', pg('S','S'), '-i0')
        neg: S->S = 'o0 = -i0;'
    """
    return InstrSpec(name, io, code='o0 = ' + op + ';')

def instr_funcall_2args_mpfr(name, io, op):
    r"""
    A helper function for creating MPFR instructions with two inputs
    and one output.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: pg = RRInterpreter().pg
        sage: instr_funcall_2args_mpfr('add', pg('SS','S'), 'mpfr_add')
        add: SS->S = 'mpfr_add(o0, i0, i1, MPFR_RNDN);'
    """
    return InstrSpec(name, io, code='%s(o0, i0, i1, MPFR_RNDN);' % op)

def instr_funcall_1arg_mpfr(name, io, op):
    r"""
    A helper function for creating MPFR instructions with one input
    and one output.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: pg = RRInterpreter().pg
        sage: instr_funcall_1arg_mpfr('exp', pg('S','S'), 'mpfr_exp')
        exp: S->S = 'mpfr_exp(o0, i0, MPFR_RNDN);'
    """
    return InstrSpec(name, io, code='%s(o0, i0, MPFR_RNDN);' % op)

class InterpreterSpec(object):
    r"""
    Each interpreter to be generated by this module is represented
    by an InterpreterSpec.
    """

    def __init__(self):
        r"""
        Initialize an InterpreterSpec.

        Initializes the following fields:

        - ``c_header`` -- a code snippet to go at the top of the C
           interpreter source file
        - ``pxd_header`` -- a code snippet to go at the top of the
           wrapper class .pxd file
        - ``pyx_header`` -- a code snippet to go at the top of the
          wrapper class source file
        - ``err_return`` -- a string indicating the value to be
          returned in case of a Python exception
        - ``mc_code`` -- a memory chunk to use for the interpreted code
        - ``extra_class_members`` -- Class members for the wrapper that
          don't correspond to memory chunks
        - ``extra_members_initialize`` -- Code to initialize
          extra_class_members

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: interp.c_header
            '#include <gsl/gsl_math.h>'
            sage: interp.pxd_header
            ''
            sage: interp.pyx_header
            'cimport sage.libs.gsl.math  # Add dependency on GSL'
            sage: interp.err_return
            '-1094648009105371'
            sage: interp.mc_code
            {MC:code}
            sage: interp = RRInterpreter()
            sage: interp.extra_class_members
            ''
            sage: interp.extra_members_initialize
            ''
        """
        self.c_header = ''
        self.pxd_header = ''
        self.pyx_header = ''
        self.err_return = 'NULL'
        self.mc_code = MemoryChunkConstants('code', ty_int)
        self.extra_class_members = ''
        self.extra_members_initialize = ''

    def _set_opcodes(self):
        r"""
        Assign opcodes to the instructions in this interpreter.

        Must be called at the end of __init__ by any subclass of
        InterpreterSpec.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: interp.instr_descs[5].opcode
            5
        """
        for i in range(len(self.instr_descs)):
            self.instr_descs[i].opcode = i


class StackInterpreter(InterpreterSpec):
    r"""
    A subclass of InterpreterSpec, specialized for stack-based
    interpreters.  (Currently all interpreters are stack-based.)
    """

    def __init__(self, type, mc_retval=None):
        r"""
        Initialize a StackInterpreter.

        INPUTS:
            type -- A StorageType; the basic type that this interpreter
                    operates on
            mc_retval -- default None; if not None, a special-purpose
                         MemoryChunk to use as a return value

        Initializes the fields described in the documentation for
        InterpreterSpec.__init__, as well as the following:

        mc_args, mc_constants, mc_stack -- MemoryChunk values
        return_type -- the type returned by the C interpreter (None for int,
                       where 1 means success and 0 means error)
        mc_retval -- None, or the MemoryChunk to use as a return value
        ipow_range -- the range of exponents supported by the ipow
                      instruction (default is False, meaning never use ipow)
        adjust_retval -- None, or a string naming a function to call
                         in the wrapper's __call__ to modify the return
                         value of the interpreter
        implement_call_c -- True if the wrapper should have a fast cdef call_c
                            method (that bypasses the Python call overhead)
                            (default True)

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: rdf = RDFInterpreter()
            sage: rr = RRInterpreter()
            sage: el = ElementInterpreter()
            sage: rdf.mc_args
            {MC:args}
            sage: rdf.mc_constants
            {MC:constants}
            sage: rdf.mc_stack
            {MC:stack}
            sage: rr.mc_retval
            {MC:retval}
            sage: rr.return_type is None
            True
            sage: rdf.return_type.type
            'double'
            sage: rdf.implement_call_c
            True
            sage: el.implement_call_c
            False
        """
        InterpreterSpec.__init__(self)
        self.mc_args = MemoryChunkArguments('args', type)
        self.mc_constants = MemoryChunkConstants('constants', type)
        self.mc_stack = MemoryChunkScratch('stack', type, is_stack=True)
        if isinstance(type, StorageTypeAssignable):
            self.return_type = type
        else:
            self.return_type = None
        self.mc_retval = mc_retval
        self.ipow_range = False
        self.adjust_retval = None
        self.implement_call_c = True

class RDFInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    machine-floating-point values (C doubles).  This is used for
    both domain=RDF and domain=float; currently the only difference
    between the two is the type of the value returned from the
    wrapper (they use the same wrapper and interpreter).
    """

    def __init__(self):
        r"""
        Initialize an RDFInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: interp.name
            'rdf'
            sage: interp.extra_class_members
            'cdef object _domain\n'
            sage: interp.extra_members_initialize
            "self._domain = args['domain']\n"
            sage: interp.adjust_retval
            'self._domain'
            sage: interp.mc_py_constants
            {MC:py_constants}
            sage: interp.chunks
            [{MC:args}, {MC:constants}, {MC:py_constants}, {MC:stack}, {MC:code}]
            sage: interp.pg('A[D]', 'S')
            ([({MC:args}, {MC:code}, None)], [({MC:stack}, None, None)])
            sage: instrs = dict([(ins.name, ins) for ins in interp.instr_descs])
            sage: instrs['add']
            add: SS->S = 'o0 = i0 + i1;'
            sage: instrs['py_call']
            py_call: *->S = '\nPyObject *py_arg...goto error;\n}\n'

        Make sure that pow behaves reasonably::

            sage: var('x,y')
            (x, y)
            sage: ff = fast_callable(x^y, vars=[x,y], domain=RDF)
            sage: ff(1.5, 3)
            3.375
            sage: ff(-2, 3)
            -8.0
            sage: ff(-2, 1/3)
            Traceback (most recent call last):
            ...
            ValueError: negative number to a fractional power not real
        """

        StackInterpreter.__init__(self, ty_double)
        self.name = 'rdf'
        self.mc_py_constants = MemoryChunkConstants('py_constants', ty_python)
        # This is a randomly chosen number.  Whenever this number is
        # returned, the wrapper has to check whether an exception actually
        # happened, so if an expression evaluates to this number execution
        # is slightly slower.  Hopefully that won't happen too often :)
        self.err_return = '-1094648009105371'
        self.chunks = [self.mc_args, self.mc_constants, self.mc_py_constants,
                       self.mc_stack,
                       self.mc_code]
        pg = params_gen(A=self.mc_args, C=self.mc_constants, D=self.mc_code,
                        S=self.mc_stack, P=self.mc_py_constants)
        self.pg = pg
        self.c_header = '#include <gsl/gsl_math.h>'
        self.pyx_header = 'cimport sage.libs.gsl.math  # Add dependency on GSL'
        instrs = [
            InstrSpec('load_arg', pg('A[D]', 'S'),
                       code='o0 = i0;'),
            InstrSpec('load_const', pg('C[D]', 'S'),
                       code='o0 = i0;'),
            InstrSpec('return', pg('S', ''),
                       code='return i0;'),
            InstrSpec('py_call', pg('P[D]S@D', 'S'),
                       uses_error_handler=True,
                       code="""
PyObject *py_args = PyTuple_New(n_i1);
if (py_args == NULL) goto error;
int i;
for (i = 0; i < n_i1; i++) {
  PyObject *arg = PyFloat_FromDouble(i1[i]);
  if (arg == NULL) {
    Py_DECREF(py_args);
    goto error;
  }
  PyTuple_SET_ITEM(py_args, i, arg);
}
PyObject *result = PyObject_CallObject(i0, py_args);
Py_DECREF(py_args);
if (result == NULL) goto error;
/* If result is not a float, then this will turn it into a float first. */
o0 = PyFloat_AsDouble(result);
Py_DECREF(result);
if (o0 == -1 && PyErr_Occurred()) {
  goto error;
}
"""),
            InstrSpec('pow', pg('SS', 'S'),
                       uses_error_handler=True,
                       code="""
/* See python's pow in floatobject.c */
if (i0 == 0) o0 = 1.0;
else {
    if (i0 < 0 && i1 != floor(i1)) {
        PyErr_SetString(PyExc_ValueError, "negative number to a fractional power not real");
        goto error;
    }
    o0 = pow(i0, i1);
}
""")
            ]
        for (name, op) in [('add', '+'), ('sub', '-'),
                           ('mul', '*'), ('div', '/')]:
            instrs.append(instr_infix(name, pg('SS', 'S'), op))
        instrs.append(instr_funcall_2args('ipow', pg('SD', 'S'), 'gsl_pow_int'))
        for (name, op) in [('neg', '-i0'), ('invert', '1/i0'),
                           ('abs', 'fabs(i0)')]:
            instrs.append(instr_unary(name, pg('S', 'S'), op))
        for name in ['sqrt', 'ceil', 'floor', 'sin', 'cos', 'tan',
                     'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh',
                     'asinh', 'acosh', 'atanh', 'exp', 'log']:
            instrs.append(instr_unary(name, pg('S',  'S'), "%s(i0)" % name))
        self.instr_descs = instrs
        self._set_opcodes()
        # supported for exponents that fit in an int
        self.ipow_range = (int(-2**31), int(2**31-1))
        self.extra_class_members = "cdef object _domain\n"
        self.extra_members_initialize = "self._domain = args['domain']\n"
        self.adjust_retval = 'self._domain'


class CDFInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    complex machine-floating-point values (C doubles).
    """

    def __init__(self):
        r"""
        Initialize a CDFInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = CDFInterpreter()
            sage: interp.name
            'cdf'
            sage: interp.mc_py_constants
            {MC:py_constants}
            sage: interp.chunks
            [{MC:args}, {MC:constants}, {MC:py_constants}, {MC:stack}, {MC:code}]
            sage: interp.pg('A[D]', 'S')
            ([({MC:args}, {MC:code}, None)], [({MC:stack}, None, None)])
            sage: instrs = dict([(ins.name, ins) for ins in interp.instr_descs])
            sage: instrs['add']
            add: SS->S = 'o0 = i0 + i1;'
            sage: instrs['sin']
            sin: S->S = 'o0 = csin(i0);'
            sage: instrs['py_call']
            py_call: *->S = '\nif (!cdf_py_call_...goto error;\n}\n'

        A test of integer powers::

            sage: f(x) = sum(x^k for k in [-20..20])
            sage: f(CDF(1+2j))  # rel tol 4e-16
            -10391778.999999996 + 3349659.499999962*I
            sage: ff = fast_callable(f, CDF)
            sage: ff(1 + 2j)  # rel tol 1e-14
            -10391779.000000004 + 3349659.49999997*I
            sage: ff.python_calls()
            []

            sage: f(x) = sum(x^k for k in [0..5])
            sage: ff = fast_callable(f, CDF)
            sage: ff(2)
            63.0
            sage: ff(2j)
            13.0 + 26.0*I
        """

        StackInterpreter.__init__(self, ty_double_complex)
        self.name = 'cdf'
        self.mc_py_constants = MemoryChunkConstants('py_constants', ty_python)
        # See comment for RDFInterpreter
        self.err_return = '-1094648119105371'
        self.adjust_retval = "dz_to_CDE"
        self.chunks = [self.mc_args, self.mc_constants, self.mc_py_constants,
                       self.mc_stack,
                       self.mc_code]
        pg = params_gen(A=self.mc_args, C=self.mc_constants, D=self.mc_code,
                        S=self.mc_stack, P=self.mc_py_constants)
        self.pg = pg
        self.c_header = """
#include <stdlib.h>
#include <complex.h>
#include "interpreters/wrapper_cdf.h"

/* On Solaris, we need to define _Imaginary_I when compiling with GCC,
 * otherwise the constant I doesn't work. The definition below is based
 * on glibc. */
#ifdef __GNUC__
#undef  _Imaginary_I
#define _Imaginary_I  (__extension__ 1.0iF)
#endif

typedef double complex double_complex;

static inline double complex csquareX(double complex z) {
    double complex res;
    __real__(res) = __real__(z) * __real__(z) - __imag__(z) * __imag__(z);
    __imag__(res) = 2 * __real__(z) * __imag__(z);
    return res;
}

static inline double complex cpow_int(double complex z, int exp) {
    if (exp < 0) return 1/cpow_int(z, -exp);
    switch (exp) {
        case 0: return 1;
        case 1: return z;
        case 2: return csquareX(z);
        case 3: return csquareX(z) * z;
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        {
            double complex z2 = csquareX(z);
            double complex z4 = csquareX(z2);
            if (exp == 4) return z4;
            if (exp == 5) return z4 * z;
            if (exp == 6) return z4 * z2;
            if (exp == 7) return z4 * z2 * z;
            if (exp == 8) return z4 * z4;
        }
    }
    if (cimag(z) == 0) return pow(creal(z), exp);
    if (creal(z) == 0) {
        double r = pow(cimag(z), exp);
        switch (exp % 4) {
            case 0:
                return r;
            case 1:
                return r * I;
            case 2:
                return -r;
            default /* case 3 */:
                return -r * I;
        }
    }
    return cpow(z, exp);
}
"""
        self.pxd_header = """
# This is to work around a header incompatibility with PARI using
# "I" as variable conflicting with the complex "I".
# If we cimport pari earlier, we avoid this problem.
cimport sage.libs.cypari2.types

# We need the type double_complex to work around
#   http://trac.cython.org/ticket/869
# so this is a bit hackish.
cdef extern from "complex.h":
    ctypedef double double_complex "double complex"
"""
        self.pyx_header = """
from sage.rings.complex_double cimport ComplexDoubleElement
import sage.rings.complex_double
cdef object CDF = sage.rings.complex_double.CDF

cdef extern from "solaris_fixes.h":
    pass

cdef extern from "complex.h":
    cdef double creal(double_complex)
    cdef double cimag(double_complex)
    cdef double_complex _Complex_I

cdef inline double_complex CDE_to_dz(zz):
    cdef ComplexDoubleElement z = <ComplexDoubleElement>(zz if isinstance(zz, ComplexDoubleElement) else CDF(zz))
    return z._complex.dat[0] + _Complex_I * z._complex.dat[1]

cdef inline ComplexDoubleElement dz_to_CDE(double_complex dz):
    cdef ComplexDoubleElement z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
    z._complex.dat[0] = creal(dz)
    z._complex.dat[1] = cimag(dz)
    return z

cdef public bint cdf_py_call_helper(object fn,
                                    int n_args,
                                    double_complex* args, double_complex* retval) except 0:
    py_args = []
    cdef int i
    for i from 0 <= i < n_args:
        py_args.append(dz_to_CDE(args[i]))
    py_result = fn(*py_args)
    cdef ComplexDoubleElement result
    if isinstance(py_result, ComplexDoubleElement):
        result = <ComplexDoubleElement>py_result
    else:
        result = CDF(py_result)
    retval[0] = CDE_to_dz(result)
    return 1

"""[1:]

        instrs = [
            InstrSpec('load_arg', pg('A[D]', 'S'),
                       code='o0 = i0;'),
            InstrSpec('load_const', pg('C[D]', 'S'),
                       code='o0 = i0;'),
            InstrSpec('return', pg('S', ''),
                       code='return i0;'),
            InstrSpec('py_call', pg('P[D]S@D', 'S'),
                       uses_error_handler=True,
                       code="""
if (!cdf_py_call_helper(i0, n_i1, i1, &o0)) {
  goto error;
}
""")
            ]
        for (name, op) in [('add', '+'), ('sub', '-'),
                           ('mul', '*'), ('div', '/')]:
            instrs.append(instr_infix(name, pg('SS', 'S'), op))
        instrs.append(instr_funcall_2args('pow', pg('SS', 'S'), 'cpow'))
        instrs.append(instr_funcall_2args('ipow', pg('SD', 'S'), 'cpow_int'))
        for (name, op) in [('neg', '-i0'), ('invert', '1/i0'),
                           ('abs', 'cabs(i0)')]:
            instrs.append(instr_unary(name, pg('S', 'S'), op))
        for name in ['sqrt', 'sin', 'cos', 'tan',
                     'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh',
                     'asinh', 'acosh', 'atanh', 'exp', 'log']:
            instrs.append(instr_unary(name, pg('S',  'S'), "c%s(i0)" % name))
        self.instr_descs = instrs
        self._set_opcodes()
        # supported for exponents that fit in an int
        self.ipow_range = (int(-2**31), int(2**31-1))


class RRInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    MPFR arbitrary-precision floating-point numbers.
    """

    def __init__(self):
        r"""
        Initialize an RDFInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RRInterpreter()
            sage: interp.name
            'rr'
            sage: interp.mc_py_constants
            {MC:py_constants}
            sage: interp.chunks
            [{MC:args}, {MC:retval}, {MC:constants}, {MC:py_constants}, {MC:stack}, {MC:code}, {MC:domain}]
            sage: interp.pg('A[D]', 'S')
            ([({MC:args}, {MC:code}, None)], [({MC:stack}, None, None)])
            sage: instrs = dict([(ins.name, ins) for ins in interp.instr_descs])
            sage: instrs['add']
            add: SS->S = 'mpfr_add(o0, i0, i1, MPFR_RNDN);'
            sage: instrs['py_call']
            py_call: *->S = '\nif (!rr_py_call_h...goto error;\n}\n'

        That py_call instruction is particularly interesting, and
        demonstrates a useful technique to let you use Cython code
        in an interpreter.  Let's look more closely::

            sage: print(instrs['py_call'].code)
            if (!rr_py_call_helper(domain, i0, n_i1, i1, o0)) {
              goto error;
            }

        This instruction makes use of the function ``rr_py_call_helper``,
        which is declared in ``wrapper_rr.h``::

            sage: print(interp.c_header)
            <BLANKLINE>
            #include <mpfr.h>
            #include "interpreters/wrapper_rr.h"
            <BLANKLINE>

        The function ``rr_py_call_helper`` is implemented in Cython::

            sage: print(interp.pyx_header)
            # distutils: libraries = mpfr gmp
            <BLANKLINE>
            cdef public bint rr_py_call_helper(object domain, object fn,
                                               int n_args,
                                               mpfr_t* args, mpfr_t retval) except 0:
                py_args = []
                cdef int i
                cdef RealNumber rn
                for i from 0 <= i < n_args:
                    rn = domain()
                    mpfr_set(rn.value, args[i], MPFR_RNDN)
                    py_args.append(rn)
                cdef RealNumber result = domain(fn(*py_args))
                mpfr_set(retval, result.value, MPFR_RNDN)
                return 1


        So instructions where you need to interact with Python can
        call back into Cython code fairly easily.
        """

        StackInterpreter.__init__(self, ty_mpfr, mc_retval= MemoryChunkRRRetval('retval', ty_mpfr))
        self.name = 'rr'
        self.err_return = '0'
        self.mc_py_constants = MemoryChunkConstants('py_constants', ty_python)
        self.mc_domain = MemoryChunkPyConstant('domain')
        self.chunks = [self.mc_args, self.mc_retval, self.mc_constants,
                       self.mc_py_constants,
                       self.mc_stack, self.mc_code, self.mc_domain]
        pg = params_gen(A=self.mc_args, C=self.mc_constants, D=self.mc_code,
                        S=self.mc_stack,
                        P=self.mc_py_constants)
        self.pg = pg
        self.c_header = '''
#include <mpfr.h>
#include "interpreters/wrapper_rr.h"
'''
        self.pxd_header = """
from sage.rings.real_mpfr cimport RealField_class, RealNumber
from sage.libs.mpfr cimport *
"""
        self.pyx_header = """# distutils: libraries = mpfr gmp

cdef public bint rr_py_call_helper(object domain, object fn,
                                   int n_args,
                                   mpfr_t* args, mpfr_t retval) except 0:
    py_args = []
    cdef int i
    cdef RealNumber rn
    for i from 0 <= i < n_args:
        rn = domain()
        mpfr_set(rn.value, args[i], MPFR_RNDN)
        py_args.append(rn)
    cdef RealNumber result = domain(fn(*py_args))
    mpfr_set(retval, result.value, MPFR_RNDN)
    return 1

"""
        instrs = [
            InstrSpec('load_arg', pg('A[D]', 'S'),
                       code='mpfr_set(o0, i0, MPFR_RNDN);'),
            InstrSpec('load_const', pg('C[D]', 'S'),
                       code='mpfr_set(o0, i0, MPFR_RNDN);'),
            InstrSpec('return', pg('S', ''),
                       code='mpfr_set(retval, i0, MPFR_RNDN);\nreturn 1;\n'),
            InstrSpec('py_call', pg('P[D]S@D', 'S'),
                       uses_error_handler=True,
                       code="""
if (!rr_py_call_helper(domain, i0, n_i1, i1, o0)) {
  goto error;
}
""")
            ]
        for (name, op) in [('add', 'mpfr_add'), ('sub', 'mpfr_sub'),
                           ('mul', 'mpfr_mul'), ('div', 'mpfr_div'),
                           ('pow', 'mpfr_pow')]:
            instrs.append(instr_funcall_2args_mpfr(name, pg('SS', 'S'), op))
        instrs.append(instr_funcall_2args_mpfr('ipow', pg('SD', 'S'), 'mpfr_pow_si'))
        for name in ['neg', 'abs',
                     'log', 'log2', 'log10',
                     'exp', 'exp2', 'exp10',
                     'cos', 'sin', 'tan',
                     'sec', 'csc', 'cot',
                     'acos', 'asin', 'atan',
                     'cosh', 'sinh', 'tanh',
                     'sech', 'csch', 'coth',
                     'acosh', 'asinh', 'atanh',
                     'log1p', 'expm1', 'eint',
                     'gamma', 'lngamma',
                     'zeta', 'erf', 'erfc',
                     'j0', 'j1', 'y0', 'y1']:
            instrs.append(instr_funcall_1arg_mpfr(name, pg('S', 'S'), 'mpfr_' + name))
        # mpfr_ui_div constructs a temporary mpfr_t and then calls mpfr_div;
        # it would probably be (slightly) faster to use a permanent copy
        # of "one" (on the other hand, the constructed temporary copy is
        # on the stack, so it's very likely to be in the cache).
        instrs.append(InstrSpec('invert', pg('S', 'S'),
                                 code='mpfr_ui_div(o0, 1, i0, MPFR_RNDN);'))
        self.instr_descs = instrs
        self._set_opcodes()
        # Supported for exponents that fit in a long, so we could use
        # a much wider range on a 64-bit machine.  On the other hand,
        # it's easier to write the code this way, and constant integer
        # exponents outside this range probably aren't very common anyway.
        self.ipow_range = (int(-2**31), int(2**31-1))

class PythonInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    Python objects.

    Let's discuss how the reference-counting works in Python-object
    based interpreters.

    There is a simple rule to remember: when executing the code
    snippets, the input variables contain borrowed references;
    you must fill in the output variables with references you own.

    As an optimization, an instruction may set .handles_own_decref; in
    that case, it must decref any input variables that came from the
    stack.  (Input variables that came from arguments/constants chunks
    must NOT be decref'ed!)  In addition, with .handles_own_decref, if
    any of your input variables are arbitrary-count, then you must
    NULL out these variables as you decref them.  (Use Py_CLEAR to do
    this, unless you understand the documentation of Py_CLEAR and why
    it's different than Py_XDECREF followed by assigning NULL.)

    Note that as a tiny optimization, the interpreter always assumes
    (and ensures) that empty parts of the stack contain NULL, so
    it doesn't bother to Py_XDECREF before it pushes onto the stack.
    """

    def __init__(self):
        r"""
        Initialize a PythonInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = PythonInterpreter()
            sage: interp.name
            'py'
            sage: interp.mc_args
            {MC:args}
            sage: interp.chunks
            [{MC:args}, {MC:constants}, {MC:stack}, {MC:code}]
            sage: instrs = dict([(ins.name, ins) for ins in interp.instr_descs])
            sage: instrs['add']
            add: SS->S = 'o0 = PyNumber_Add(i0, i1);'
            sage: instrs['py_call']
            py_call: *->S = '\nPyObject *py_args...CREF(py_args);\n'
        """

        StackInterpreter.__init__(self, ty_python)
        self.name = 'py'
        # StackInterpreter.__init__ gave us a MemoryChunkArguments.
        # Override with MemoryChunkPythonArguments.
        self.mc_args = MemoryChunkPythonArguments('args', ty_python)
        self.chunks = [self.mc_args, self.mc_constants, self.mc_stack,
                       self.mc_code]
        pg = params_gen(A=self.mc_args, C=self.mc_constants, D=self.mc_code,
                        S=self.mc_stack)
        self.pg = pg
        self.c_header = """
#define CHECK(x) (x != NULL)
"""
        instrs = [
            InstrSpec('load_arg', pg('A[D]', 'S'),
                       code='o0 = i0; Py_INCREF(o0);'),
            InstrSpec('load_const', pg('C[D]', 'S'),
                       code='o0 = i0; Py_INCREF(o0);'),
            InstrSpec('return', pg('S', ''),
                       code='return i0;',
                       handles_own_decref=True),
            InstrSpec('py_call', pg('C[D]S@D', 'S'),
                       handles_own_decref=True,
                       code="""
PyObject *py_args = PyTuple_New(n_i1);
if (py_args == NULL) goto error;
int i;
for (i = 0; i < n_i1; i++) {
  PyObject *arg = i1[i];
  PyTuple_SET_ITEM(py_args, i, arg);
  i1[i] = NULL;
}
o0 = PyObject_CallObject(i0, py_args);
Py_DECREF(py_args);
""")
            ]
        for (name, op) in [('add', 'PyNumber_Add'),
                           ('sub', 'PyNumber_Subtract'),
                           ('mul', 'PyNumber_Multiply'),
                           ('div', 'PyNumber_Divide')]:
            instrs.append(instr_funcall_2args(name, pg('SS', 'S'), op))
        instrs.append(InstrSpec('pow', pg('SS', 'S'),
                                code='o0 = PyNumber_Power(i0, i1, Py_None);'))
        instrs.append(InstrSpec('ipow', pg('SC[D]', 'S'),
                                code='o0 = PyNumber_Power(i0, i1, Py_None);'))
        for (name, op) in [('neg', 'PyNumber_Negative'),
                           ('invert', 'PyNumber_Invert'),
                           ('abs', 'PyNumber_Absolute')]:
            instrs.append(instr_unary(name, pg('S', 'S'), '%s(i0)'%op))
        self.instr_descs = instrs
        self._set_opcodes()
        # Always use ipow
        self.ipow_range = True
        # We don't yet support call_c for Python-object interpreters
        # (the default implementation doesn't work, because of
        # object vs. PyObject* confusion)
        self.implement_call_c = False

class ElementInterpreter(PythonInterpreter):
    r"""
    A subclass of PythonInterpreter, specifying an interpreter over
    Sage elements with a particular parent.

    This is very similar to the PythonInterpreter, but after every
    instruction, the result is checked to make sure it actually an
    element with the correct parent; if not, we attempt to convert it.

    Uses the same instructions (with the same implementation) as
    PythonInterpreter.
    """

    def __init__(self):
        r"""
        Initialize an ElementInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = ElementInterpreter()
            sage: interp.name
            'el'
            sage: interp.mc_args
            {MC:args}
            sage: interp.chunks
            [{MC:args}, {MC:constants}, {MC:stack}, {MC:domain}, {MC:code}]
            sage: instrs = dict([(ins.name, ins) for ins in interp.instr_descs])
            sage: instrs['add']
            add: SS->S = 'o0 = PyNumber_Add(i0, i1);'
            sage: instrs['py_call']
            py_call: *->S = '\nPyObject *py_args...CREF(py_args);\n'
        """

        PythonInterpreter.__init__(self)
        self.name = 'el'
        # PythonInterpreter.__init__ gave us a MemoryChunkPythonArguments.
        # Override with MemoryChunkElementArguments.
        self.mc_args = MemoryChunkElementArguments('args', ty_python)
        self.mc_domain_info = MemoryChunkPyConstant('domain')
        self.chunks = [self.mc_args, self.mc_constants, self.mc_stack,
                       self.mc_domain_info, self.mc_code]
        self.c_header = """
#include "interpreters/wrapper_el.h"

#define CHECK(x) do_check(&(x), domain)

static inline int do_check(PyObject **x, PyObject *domain) {
  if (*x == NULL) return 0;
  PyObject *new_x = el_check_element(*x, domain);
  Py_DECREF(*x);
  *x = new_x;
  if (*x == NULL) return 0;
  return 1;
}
"""
        self.pyx_header = """
from sage.structure.element cimport Element

cdef public object el_check_element(object v, parent):
    cdef Element v_el

    if isinstance(v, Element):
        v_el = <Element>v
        if v_el._parent is parent:
            return v_el

    return parent(v)

"""[1:]

class InterpreterGenerator(object):
    r"""
    This class takes an InterpreterSpec and generates the corresponding
    C interpreter and Cython wrapper.

    See the documentation for methods get_wrapper and get_interpreter
    for more information.
    """

    def __init__(self, spec):
        r"""
        Initialize an InterpreterGenerator.

        INPUT:

        - ``spec`` -- an InterpreterSpec

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: gen = InterpreterGenerator(interp)
            sage: gen._spec is interp
            True
            sage: gen.uses_error_handler
            False
        """
        self._spec = spec
        self.uses_error_handler = False

    def gen_code(self, instr_desc, write):
        r"""
        Generates code for a single instruction.

        INPUTS:
            instr_desc -- an InstrSpec
            write -- a Python callable

        This function calls its write parameter successively with
        strings; when these strings are concatenated, the result is
        the code for the given instruction.

        See the documentation for the get_interpreter method for more
        information.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: gen = InterpreterGenerator(interp)
            sage: from six.moves import cStringIO as StringIO
            sage: buff = StringIO()
            sage: instrs = dict([(ins.name, ins) for ins in interp.instr_descs])
            sage: gen.gen_code(instrs['div'], buff.write)
            sage: print(buff.getvalue())
                case 8: /* div */
                  {
                    double i1 = *--stack;
                    double i0 = *--stack;
                    double o0;
                    o0 = i0 / i1;
                    *stack++ = o0;
                  }
                  break;
            <BLANKLINE>
        """

        d = instr_desc
        w = write
        s = self._spec

        if d.uses_error_handler:
            self.uses_error_handler = True

        w(je("""
    case {{ d.opcode }}: /* {{ d.name }} */
      {
""", d=d))

        # If the inputs to an instruction come from the stack,
        # then we want to generate code for the inputs in reverse order:
        # for instance, the divide instruction, which takes inputs A and B
        # and generates A/B, needs to pop B off the stack first.
        # On the other hand, if the inputs come from the constant pool,
        # then we want to generate code for the inputs in normal order,
        # because the addresses in the code stream will be in that order.
        # We handle this by running through the inputs in two passes:
        # first a forward pass, where we handle non-stack inputs
        # (and lengths for stack inputs), and then a reverse pass,
        # where we handle stack inputs.
        for i in range(len(d.inputs)):
            (ch, addr, input_len) = d.inputs[i]
            chst = ch.storage_type
            if addr is not None:
                w("        int ai%d = %s;\n" % (i, string_of_addr(addr)))
            if input_len is not None:
                w("        int n_i%d = %s;\n" % (i, string_of_addr(input_len)))
            if not ch.is_stack():
                # Shouldn't hardcode 'code' here
                if ch.name == 'code':
                    w("        %s i%d = %s;\n" % (chst.c_local_type(), i, string_of_addr(ch)))
                elif input_len is not None:
                    w("        %s i%d = %s + ai%d;\n" %
                      (chst.c_ptr_type(), i, ch.name, i))
                else:
                    w("        %s i%d = %s[ai%d];\n" %
                      (chst.c_local_type(), i, ch.name, i))

        for i in reversed(range(len(d.inputs))):
            (ch, addr, input_len) = d.inputs[i]
            chst = ch.storage_type
            if ch.is_stack():
                if input_len is not None:
                    w("        %s -= n_i%d;\n" % (ch.name, i))
                    w("        %s i%d = %s;\n" % (chst.c_ptr_type(), i, ch.name))
                else:
                    w("        %s i%d = *--%s;\n" % (chst.c_local_type(), i, ch.name))
                    if ch.is_python_refcounted_stack():
                        w("        *%s = NULL;\n" % ch.name)

        for i in range(len(d.outputs)):
            (ch, addr, output_len) = d.outputs[i]
            chst = ch.storage_type
            if addr is not None:
                w("        int ao%d = %s;\n" % (i, string_of_addr(addr)))
            if output_len is not None:
                w("        int n_o%d = %s;\n" % (i, string_of_addr(output_len)))
                if ch.is_stack():
                    w("        %s o%d = %s;\n" %
                      (chst.c_ptr_type(), i, ch.name))
                    w("        %s += n_o%d;\n" % (ch.name, i))
                else:
                    w("        %s o%d = %s + ao%d;\n" %
                      (chst.c_ptr_type(), i, ch.name, i))
            else:
                if not chst.cheap_copies():
                    if ch.is_stack():
                        w("        %s o%d = *%s++;\n" %
                          (chst.c_local_type(), i, ch.name))
                    else:
                        w("        %s o%d = %s[ao%d];\n" %
                          (chst.c_local_type(), i, ch.name, i))
                else:
                    w("        %s o%d;\n" % (chst.c_local_type(), i))
        w(indent_lines(8, d.code.rstrip('\n') + '\n'))

        stack_offsets = defaultdict(int)
        for i in range(len(d.inputs)):
            (ch, addr, input_len) = d.inputs[i]
            chst = ch.storage_type
            if ch.is_python_refcounted_stack() and not d.handles_own_decref:
                if input_len is None:
                    w("        Py_DECREF(i%d);\n" % i)
                    stack_offsets[ch] += 1
                else:
                    w(je("""
        int {{ iter }};
        for ({{ iter }} = 0; {{ iter }} < n_i{{ i }}; {{ iter }}++) {
          Py_CLEAR(i{{ i }}[{{ iter }}]);
        }
""", iter='_interp_iter_%d' % i, i=i))

        for i in range(len(d.outputs)):
            ch = d.outputs[i][0]
            chst = ch.storage_type
            if chst.python_refcounted():
                # We don't yet support code chunks
                # that produce multiple Python values, because of
                # the way it complicates error handling.
                assert i == 0
                w("        if (!CHECK(o%d)) {\n" % i)
                w("          Py_XDECREF(o%d);\n" % i)
                w("          goto error;\n")
                w("        }\n")
                self.uses_error_handler = True
            if chst.cheap_copies():
                if ch.is_stack():
                    w("        *%s++ = o%d;\n" % (ch.name, i))
                else:
                    w("        %s[ao%d] = o%d;\n" % (ch.name, i, i))

        w(je("""
      }
      break;
"""))

    def func_header(self, cython=False):
        r"""
        Generates the function header for the declaration (in the Cython
        wrapper) or the definition (in the C interpreter) of the interpreter
        function.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = ElementInterpreter()
            sage: gen = InterpreterGenerator(interp)
            sage: print(gen.func_header())
            PyObject* interp_el(PyObject** args,
                    PyObject** constants,
                    PyObject** stack,
                    PyObject* domain,
                    int* code)
            sage: print(gen.func_header(cython=True))
            object interp_el(PyObject** args,
                    PyObject** constants,
                    PyObject** stack,
                    PyObject* domain,
                    int* code)
        """
        s = self._spec
        ret_ty = 'bint' if cython else 'int'
        if s.return_type:
            ret_ty = s.return_type.c_decl_type()
            if cython:
                ret_ty = s.return_type.cython_decl_type()
        return je("""{{ ret_ty }} interp_{{ s.name }}(
{%- for ch in s.chunks %}
{%    if not loop.first %},
        {% endif %}{{ ch.declare_parameter() }}
{%- endfor %})""", ret_ty=ret_ty, s=s)

    def write_interpreter(self, write):
        r"""
        Generate the code for the C interpreter.

        This function calls its write parameter successively with
        strings; when these strings are concatenated, the result is
        the code for the interpreter.

        See the documentation for the get_interpreter method for more
        information.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: gen = InterpreterGenerator(interp)
            sage: from six.moves import cStringIO as StringIO
            sage: buff = StringIO()
            sage: gen.write_interpreter(buff.write)
            sage: print(buff.getvalue())
            /* Automatically generated by ...
        """
        s = self._spec
        w = write
        w(je("""
/* {{ warn }} */
#include <Python.h>
{% print(s.c_header) %}

{{ myself.func_header() }} {
  while (1) {
    switch (*code++) {
""", s=s, myself=self, i=indent_lines, warn=autogen_warn))
        for instr_desc in s.instr_descs:
            self.gen_code(instr_desc, w)
        w(je("""
    }
  }
{% if myself.uses_error_handler %}
error:
  return {{ s.err_return }};
{% endif %}
}

""", s=s, i=indent_lines, myself=self))

    def write_wrapper(self, write):
        r"""
        Generate the code for the Cython wrapper.
        This function calls its write parameter successively with
        strings; when these strings are concatenated, the result is
        the code for the wrapper.

        See the documentation for the get_wrapper method for more
        information.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: gen = InterpreterGenerator(interp)
            sage: from six.moves import cStringIO as StringIO
            sage: buff = StringIO()
            sage: gen.write_wrapper(buff.write)
            sage: print(buff.getvalue())
            # Automatically generated by ...
        """
        s = self._spec
        w = write
        types = set()
        do_cleanup = False
        for ch in s.chunks:
            if ch.storage_type is not None:
                types.add(ch.storage_type)
            do_cleanup = do_cleanup or ch.needs_cleanup_on_error()
        for ch in s.chunks:
            if ch.name == 'args':
                arg_ch = ch

        the_call = je("""
        {% if s.return_type %}return {% endif -%}
{% if s.adjust_retval %}{{ s.adjust_retval }}({% endif %}
interp_{{ s.name }}({{ arg_ch.pass_argument() }}
{% for ch in s.chunks[1:] %}
            , {{ ch.pass_argument() }}
{% endfor %}
            ){% if s.adjust_retval %}){% endif %}

""", s=s, arg_ch=arg_ch)

        the_call_c = je("""
        {% if s.return_type %}result[0] = {% endif %}
interp_{{ s.name }}(args
{% for ch in s.chunks[1:] %}
            , {{ ch.pass_call_c_argument() }}
{% endfor %}
            )

""", s=s, arg_ch=arg_ch)

        w(je("""
# {{ warn }}
# distutils: sources = sage/ext/interpreters/interp_{{ s.name }}.c
{{ s.pyx_header }}

include "cysignals/memory.pxi"
from cpython.ref cimport PyObject
cdef extern from "Python.h":
    void Py_DECREF(PyObject *o)
    void Py_INCREF(PyObject *o)
    void Py_CLEAR(PyObject *o)

    object PyList_New(Py_ssize_t len)
    ctypedef struct PyListObject:
        PyObject **ob_item

    ctypedef struct PyTupleObject:
        PyObject **ob_item

from sage.ext.fast_callable cimport Wrapper

cdef extern:
    {{ myself.func_header(cython=true) -}}
{% if s.err_return != 'NULL' %}
 except? {{ s.err_return }}
{% endif %}

cdef class Wrapper_{{ s.name }}(Wrapper):
    # attributes are declared in corresponding .pxd file

    def __init__(self, args):
        Wrapper.__init__(self, args, metadata)
        cdef int i
        cdef int count
{% for ty in types %}
{% print(indent_lines(8, ty.local_declarations)) %}
{% print(indent_lines(8, ty.class_member_initializations)) %}
{% endfor %}
{% for ch in s.chunks %}
{% print(ch.init_class_members()) %}
{% endfor %}
{% print(indent_lines(8, s.extra_members_initialize)) %}

    def __dealloc__(self):
        cdef int i
{% for ch in s.chunks %}
{% print(ch.dealloc_class_members()) %}
{% endfor %}

    def __call__(self, *args):
        if self._n_args != len(args): raise ValueError
{% for ty in types %}
{% print(indent_lines(8, ty.local_declarations)) %}
{% endfor %}
{% print(indent_lines(8, arg_ch.setup_args())) %}
{% for ch in s.chunks %}
{% print(ch.declare_call_locals()) %}
{% endfor %}
{% if do_cleanup %}
        try:
{% print(indent_lines(4, the_call)) %}
        except BaseException:
{%   for ch in s.chunks %}
{%     if ch.needs_cleanup_on_error() %}
{%       print(indent_lines(12, ch.handle_cleanup())) %}
{%     endif %}
{%   endfor %}
            raise
{% else %}
{% print(the_call) %}
{% endif %}
{% if not s.return_type %}
        return retval
{% endif %}

{% if s.implement_call_c %}
    cdef bint call_c(self,
                     {{ arg_ch.storage_type.c_ptr_type() }} args,
                     {{ arg_ch.storage_type.c_reference_type() }} result) except 0:
{% if do_cleanup %}
        try:
{% print(indent_lines(4, the_call_c)) %}
        except BaseException:
{%   for ch in s.chunks %}
{%     if ch.needs_cleanup_on_error() %}
{%       print(indent_lines(12, ch.handle_cleanup())) %}
{%     endif %}
{%   endfor %}
            raise
{% else %}
{% print(the_call_c) %}
{% endif %}
        return 1
{% endif %}

from sage.ext.fast_callable import CompilerInstrSpec, InterpreterMetadata
metadata = InterpreterMetadata(by_opname={
{% for instr in s.instr_descs %}
  '{{ instr.name }}':
  (CompilerInstrSpec({{ instr.n_inputs }}, {{ instr.n_outputs }}, {{ instr.parameters }}), {{ instr.opcode }}),
{% endfor %}
 },
 by_opcode=[
{% for instr in s.instr_descs %}
  ('{{ instr.name }}',
   CompilerInstrSpec({{ instr.n_inputs }}, {{ instr.n_outputs }}, {{ instr.parameters }})),
{% endfor %}
 ],
 ipow_range={{ s.ipow_range }})
""", s=s, myself=self, types=types, arg_ch=arg_ch,
     indent_lines=indent_lines, the_call=the_call,
     the_call_c=the_call_c, do_cleanup=do_cleanup, warn=autogen_warn))

    def write_pxd(self, write):
        r"""
        Generate the pxd file for the Cython wrapper.
        This function calls its write parameter successively with
        strings; when these strings are concatenated, the result is
        the code for the pxd file.

        See the documentation for the get_pxd method for more
        information.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: interp = RDFInterpreter()
            sage: gen = InterpreterGenerator(interp)
            sage: from six.moves import cStringIO as StringIO
            sage: buff = StringIO()
            sage: gen.write_pxd(buff.write)
            sage: print(buff.getvalue())
            # Automatically generated by ...
        """
        s = self._spec
        w = write
        types = set()
        for ch in s.chunks:
            if ch.storage_type is not None:
                types.add(ch.storage_type)
        for ch in s.chunks:
            if ch.name == 'args':
                arg_ch = ch

        w(je("""
# {{ warn }}

from cpython cimport PyObject

from sage.ext.fast_callable cimport Wrapper
{% print(s.pxd_header) %}

cdef class Wrapper_{{ s.name }}(Wrapper):
{% for ty in types %}
{% print(indent_lines(4, ty.class_member_declarations)) %}
{% endfor %}
{% for ch in s.chunks %}
{% print(ch.declare_class_members()) %}
{% endfor %}
{% print(indent_lines(4, s.extra_class_members)) %}
{% if s.implement_call_c %}
    cdef bint call_c(self,
                     {{ arg_ch.storage_type.c_ptr_type() }} args,
                     {{ arg_ch.storage_type.c_reference_type() }} result) except 0
{% endif %}
""", s=s, myself=self, types=types, indent_lines=indent_lines,
     arg_ch=arg_ch, warn=autogen_warn))

    def get_interpreter(self):
        r"""
        Return the code for the C interpreter.

        EXAMPLES:

        First we get the InterpreterSpec for several interpreters::

            sage: from sage_setup.autogen.interpreters import *
            sage: rdf_spec = RDFInterpreter()
            sage: rr_spec = RRInterpreter()
            sage: el_spec = ElementInterpreter()

        Then we get the actual interpreter code::

            sage: rdf_interp = InterpreterGenerator(rdf_spec).get_interpreter()
            sage: rr_interp = InterpreterGenerator(rr_spec).get_interpreter()
            sage: el_interp = InterpreterGenerator(el_spec).get_interpreter()

        Now we can look through these interpreters.

        Each interpreter starts with a file header; this can be
        customized on a per-interpreter basis::

            sage: print(rr_interp)
            /* Automatically generated by ... */
            ...

        Next is the function header, with one argument per memory chunk
        in the interpreter spec::

            sage: print(el_interp)
            /* ... */ ...
            PyObject* interp_el(PyObject** args,
                    PyObject** constants,
                    PyObject** stack,
                    PyObject* domain,
                    int* code) {
            ...

        Currently, the interpreters have a very simple structure; just
        grab the next instruction and execute it, in a switch
        statement::

            sage: print(rdf_interp)
            /* ... */ ...
              while (1) {
                switch (*code++) {
            ...

        Then comes the code for each instruction.  Here is one of the
        simplest instructions::

            sage: print(rdf_interp)
            /* ... */ ...
                case 10: /* neg */
                  {
                    double i0 = *--stack;
                    double o0;
                    o0 = -i0;
                    *stack++ = o0;
                  }
                  break;
            ...

        We simply pull the top of the stack into a variable, negate it,
        and write the result back onto the stack.

        Let's look at the MPFR-based version of this instruction.
        This is an example of an interpreter with an auto-reference
        type::

            sage: print(rr_interp)
            /* ... */ ...
                case 10: /* neg */
                  {
                    mpfr_ptr i0 = *--stack;
                    mpfr_ptr o0 = *stack++;
                    mpfr_neg(o0, i0, MPFR_RNDN);
                  }
                  break;
            ...

        Here we see that the input and output variables are actually
        just pointers into the stack.  But due to the auto-reference
        trick, the actual code snippet, ``mpfr_net(o0, i0, MPFR_RNDN);``,
        is exactly the same as if i0 and o0 were declared as local
        mpfr_t variables.

        For completeness, let's look at this instruction in the
        Python-object element interpreter::

            sage: print(el_interp)
            /* ... */ ...
                case 10: /* neg */
                  {
                    PyObject* i0 = *--stack;
                    *stack = NULL;
                    PyObject* o0;
                    o0 = PyNumber_Negative(i0);
                    Py_DECREF(i0);
                    if (!CHECK(o0)) {
                      Py_XDECREF(o0);
                      goto error;
                    }
                    *stack++ = o0;
                  }
                  break;
            ...

        The original code snippet was only ``o0 = PyNumber_Negative(i0);``;
        all the rest is automatically generated.  For ElementInterpreter,
        the CHECK macro actually checks for an exception (makes sure that
        o0 is not NULL), tests if the o0 is an element with the correct
        parent, and if not converts it into the correct parent.  (That is,
        it can potentially modify the variable o0.)
        """
        from six.moves import cStringIO as StringIO
        buff = StringIO()
        self.write_interpreter(buff.write)
        return buff.getvalue()

    def get_wrapper(self):
        r"""
        Return the code for the Cython wrapper.

        EXAMPLES:

        First we get the InterpreterSpec for several interpreters::

            sage: from sage_setup.autogen.interpreters import *
            sage: rdf_spec = RDFInterpreter()
            sage: rr_spec = RRInterpreter()
            sage: el_spec = ElementInterpreter()

        Then we get the actual wrapper code::

            sage: rdf_wrapper = InterpreterGenerator(rdf_spec).get_wrapper()
            sage: rr_wrapper = InterpreterGenerator(rr_spec).get_wrapper()
            sage: el_wrapper = InterpreterGenerator(el_spec).get_wrapper()

        Now we can look through these wrappers.

        Each wrapper starts with a file header; this can be
        customized on a per-interpreter basis (some blank lines have been
        elided below)::

            sage: print(rdf_wrapper)
            # Automatically generated by ...
            include "cysignals/memory.pxi"
            from cpython.ref cimport PyObject
            cdef extern from "Python.h":
                void Py_DECREF(PyObject *o)
                void Py_INCREF(PyObject *o)
                void Py_CLEAR(PyObject *o)
            <BLANKLINE>
                object PyList_New(Py_ssize_t len)
                ctypedef struct PyListObject:
                    PyObject **ob_item
            <BLANKLINE>
                ctypedef struct PyTupleObject:
                    PyObject **ob_item
            <BLANKLINE>
            from sage.ext.fast_callable cimport Wrapper
            ...

        We need a way to propagate exceptions back to the wrapper,
        even though we only return a double from interp_rdf.  The
        ``except? -1094648009105371`` (that's a randomly chosen
        number) means that we will return that number if there's an
        exception, but the wrapper still has to check whether that's a
        legitimate return or an exception.  (Cython does this
        automatically.)

        Next comes the actual wrapper class.  The member declarations
        are in the corresponding pxd file; see the documentation for
        get_pxd to see them::

            sage: print(rdf_wrapper)
            # ...
            cdef class Wrapper_rdf(Wrapper):
                # attributes are declared in corresponding .pxd file
            ...

        Next is the __init__ method, which starts like this::

            sage: print(rdf_wrapper)
            # ...
                def __init__(self, args):
                    Wrapper.__init__(self, args, metadata)
                    cdef int i
                    cdef int count
            ...

        To make it possible to generate code for all expression
        interpreters with a single code generator, all wrappers
        have the same API.  The __init__ method takes a single
        argument (here called *args*), which is a dictionary holding
        all the information needed to initialize this wrapper.

        We call Wrapper.__init__, which saves a copy of this arguments
        object and of the interpreter metadata in the wrapper.  (This is
        only used for debugging.)

        Now we allocate memory for each memory chunk.  (We allocate
        the memory here, and reuse it on each call of the
        wrapper/interpreter.  This is for speed reasons; in a fast
        interpreter like RDFInterpreter, there are no memory allocations
        involved in a call of the wrapper, except for the ones that
        are required by the Python calling convention.  Eventually
        we will support alternate Cython-only entry points that do
        absolutely no memory allocation.)

        Basically the same code is repeated, with minor variations, for
        each memory chunk; for brevity, we'll only show the code
        for 'constants'::

            sage: print(rdf_wrapper)
            # ...
                    val = args['constants']
                    self._n_constants = len(val)
                    self._constants = <double*>sig_malloc(sizeof(double) * len(val))
                    if self._constants == NULL: raise MemoryError
                    for i in range(len(val)):
                        self._constants[i] = val[i]
            ...

        Recall that _n_constants is an int, and _constants is a
        double*.

        The RRInterpreter version is more complicated, because it has to
        call mpfr_init::

            sage: print(rr_wrapper)
            # ...
                    cdef RealNumber rn
            ...
                    val = args['constants']
                    self._n_constants = len(val)
                    self._constants = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * len(val))
                    if self._constants == NULL: raise MemoryError
                    for i in range(len(val)):
                        mpfr_init2(self._constants[i], self.domain.prec())
                    for i in range(len(val)):
                        rn = self.domain(val[i])
                        mpfr_set(self._constants[i], rn.value, MPFR_RNDN)
            ...

        And as described in the documentation for get_pxd, in
        Python-object based interpreters we actually allocate the
        memory as a Python list::

            sage: print(el_wrapper)
            # ...
                    val = args['constants']
                    self._n_constants = len(val)
                    self._list_constants = PyList_New(self._n_constants)
                    self._constants = (<PyListObject *>self._list_constants).ob_item
                    for i in range(len(val)):
                        self._constants[i] = <PyObject *>val[i]; Py_INCREF(self._constants[i])
            ...

        Of course, once we've allocated the memory, we eventually have
        to free it.  (Again, we'll only look at 'constants'.)::

            sage: print(rdf_wrapper)
            # ...
                def __dealloc__(self):
            ...
                    if self._constants:
                        sig_free(self._constants)
            ...

        The RRInterpreter code is more complicated again because it has
        to call mpfr_clear::

            sage: print(rr_wrapper)
            # ...
                def __dealloc__(self):
                    cdef int i
            ...
                    if self._constants:
                        for i in range(self._n_constants):
                            mpfr_clear(self._constants[i])
                        sig_free(self._constants)
            ...

        But the ElementInterpreter code is extremely simple --
        it doesn't have to do anything to deallocate constants!
        (Since the memory for constants is actually allocated as a
        Python list, and Cython knows how to deallocate Python lists.)

        Finally we get to the __call__ method.  We grab the arguments
        passed by the caller, stuff them in our pre-allocated
        argument array, and then call the C interpreter.

        We optionally adjust the return value of the interpreter
        (currently only the RDF/float interpreter performs this step;
        this is the only place where domain=RDF differs than
        domain=float)::

            sage: print(rdf_wrapper)
            # ...
                def __call__(self, *args):
                    if self._n_args != len(args): raise ValueError
                    cdef double* c_args = self._args
                    cdef int i
                    for i from 0 <= i < len(args):
                        self._args[i] = args[i]
                    return self._domain(interp_rdf(c_args
                        , self._constants
                        , self._py_constants
                        , self._stack
                        , self._code
                        ))
            ...

        In Python-object based interpreters, the call to the C
        interpreter has to be a little more complicated.  We don't
        want to hold on to Python objects from an old computation by
        leaving them referenced from the stack.  In normal operation,
        the C interpreter clears out the stack as it runs, leaving the
        stack totally clear when the interpreter finishes.  However,
        this doesn't happen if the C interpreter raises an exception.
        In that case, we have to clear out any remnants from the stack
        in the wrapper::

            sage: print(el_wrapper)
            # ...
                    try:
                        return interp_el((<PyListObject*>mapped_args).ob_item
                            , self._constants
                            , self._stack
                            , <PyObject*>self._domain
                            , self._code
                            )
                    except BaseException:
                        for i in range(self._n_stack):
                            Py_CLEAR(self._stack[i])
                        raise
            ...

        Finally, we define a cdef call_c method, for quickly calling
        this object from Cython.  (The method is omitted from
        Python-object based interpreters.)::

            sage: print(rdf_wrapper)
            # ...
                cdef bint call_c(self,
                                 double* args,
                                 double* result) except 0:
                    result[0] = interp_rdf(args
                        , self._constants
                        , self._py_constants
                        , self._stack
                        , self._code
                        )
                    return 1
            ...

        The method for the RR interpreter is slightly different, because
        the interpreter takes a pointer to a result location instead of
        returning the value::

            sage: print(rr_wrapper)
            # ...
                cdef bint call_c(self,
                                 mpfr_t* args,
                                 mpfr_t result) except 0:
                    interp_rr(args
                        , result
                        , self._constants
                        , self._py_constants
                        , self._stack
                        , self._code
                        , <PyObject*>self._domain
                        )
                    return 1
            ...

        That's it for the wrapper class.  The only thing remaining is
        the interpreter metadata.  This is the information necessary
        for the code generator to map instruction names to opcodes; it
        also gives information about stack usage, etc.  This is fully
        documented at InterpreterMetadata; for now, we'll just show
        what it looks like.

        Currently, there are three parts to the metadata; the first maps
        instruction names to instruction descriptions.  The second one
        maps opcodes to instruction descriptions.  Note that we don't
        use InstrSpec objects here; instead, we use CompilerInstrSpec
        objects, which are much simpler and contain only the information
        we'll need at runtime.  The third part says what range the
        ipow instruction is defined over.

        First the part that maps instruction names to
        (CompilerInstrSpec, opcode) pairs::

            sage: print(rdf_wrapper)
            # ...
            from sage.ext.fast_callable import CompilerInstrSpec, InterpreterMetadata
            metadata = InterpreterMetadata(by_opname={
            ...
              'return':
              (CompilerInstrSpec(1, 0, []), 2),
              'py_call':
              (CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs']), 3),
              'pow':
              (CompilerInstrSpec(2, 1, []), 4),
              'add':
              (CompilerInstrSpec(2, 1, []), 5),
            ...
             }, ...)

        There's also a table that maps opcodes to (instruction name,
        CompilerInstrSpec) pairs::

            sage: print(rdf_wrapper)
            # ...
            metadata = InterpreterMetadata(...,  by_opcode=[
            ...
              ('return',
               CompilerInstrSpec(1, 0, [])),
              ('py_call',
               CompilerInstrSpec(0, 1, ['py_constants', 'n_inputs'])),
              ('pow',
               CompilerInstrSpec(2, 1, [])),
              ('add',
               CompilerInstrSpec(2, 1, [])),
            ...
             ], ...)

        And then the ipow range::

            sage: print(rdf_wrapper)
            # ...
            metadata = InterpreterMetadata(...,
              ipow_range=(-2147483648, 2147483647))

        And that's it for the wrapper.
        """
        from six.moves import cStringIO as StringIO
        buff = StringIO()
        self.write_wrapper(buff.write)
        return buff.getvalue()

    def get_pxd(self):
        r"""
        Return the code for the Cython .pxd file.

        EXAMPLES:

        First we get the InterpreterSpec for several interpreters::

            sage: from sage_setup.autogen.interpreters import *
            sage: rdf_spec = RDFInterpreter()
            sage: rr_spec = RRInterpreter()
            sage: el_spec = ElementInterpreter()

        Then we get the corresponding .pxd::

            sage: rdf_pxd = InterpreterGenerator(rdf_spec).get_pxd()
            sage: rr_pxd = InterpreterGenerator(rr_spec).get_pxd()
            sage: el_pxd = InterpreterGenerator(el_spec).get_pxd()

        Now we can look through these pxd files.

        Each .pxd starts with a file header; this can be
        customized on a per-interpreter basis (some blank lines have been
        elided below)::

            sage: print(rdf_pxd)
            # Automatically generated by ...
            from cpython cimport PyObject
            from sage.ext.fast_callable cimport Wrapper
            ...
            sage: print(rr_pxd)
            # ...
            from sage.rings.real_mpfr cimport RealField_class, RealNumber
            from sage.libs.mpfr cimport *
            ...

        Next and last is the declaration of the wrapper class, which
        starts off with a list of member declarations::

            sage: print(rdf_pxd)
            # ...
            cdef class Wrapper_rdf(Wrapper):
                cdef int _n_args
                cdef double* _args
                cdef int _n_constants
                cdef double* _constants
                cdef object _list_py_constants
                cdef int _n_py_constants
                cdef PyObject** _py_constants
                cdef int _n_stack
                cdef double* _stack
                cdef int _n_code
                cdef int* _code
            ...

        Contrast the declaration of ``_stack`` here with the
        ElementInterpreter version.  To simplify our handling of
        reference counting and garbage collection, in a Python-object
        based interpreter, we allocate arrays as Python lists,
        and then pull the array out of the innards of the list::

            sage: print(el_pxd)
            # ...
                cdef object _list_stack
                cdef int _n_stack
                cdef PyObject** _stack
            ...

        Then, at the end of the wrapper class, we declare a cdef method
        for quickly calling the wrapper object from Cython.  (This method
        is omitted from Python-object based interpreters.)::

            sage: print(rdf_pxd)
            # ...
                cdef bint call_c(self,
                                 double* args,
                                 double* result) except 0
            sage: print(rr_pxd)
            # ...
                cdef bint call_c(self,
                                 mpfr_t* args,
                                 mpfr_t result) except 0

        """
        from six.moves import cStringIO as StringIO
        buff = StringIO()
        self.write_pxd(buff.write)
        return buff.getvalue()

def write_if_changed(fn, value):
    r"""
    Write value to the file named fn, if value is different than
    the current contents.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: def last_modification(fn): return os.stat(fn).st_mtime
        sage: fn = tmp_filename('gen_interp')
        sage: write_if_changed(fn, 'Hello, world')
        sage: t1 = last_modification(fn)
        sage: open(fn).read()
        'Hello, world'
        sage: sleep(2)            # long time
        sage: write_if_changed(fn, 'Goodbye, world')
        sage: t2 = last_modification(fn)
        sage: open(fn).read()
        'Goodbye, world'
        sage: sleep(2)            # long time
        sage: write_if_changed(fn, 'Goodbye, world')
        sage: t3 = last_modification(fn)
        sage: open(fn).read()
        'Goodbye, world'
        sage: t1 == t2            # long time
        False
        sage: t2 == t3
        True
    """
    old_value = None
    try:
        with open(fn) as file:
            old_value = file.read()
    except IOError:
        pass

    if value != old_value:
        # We try to remove the file, in case it exists.  This is to
        # automatically break hardlinks... see #5350 for motivation.
        try:
            os.remove(fn)
        except OSError:
            pass

        with open(fn, 'w') as file:
            file.write(value)

def build_interp(interp_spec, dir):
    r"""
    Given an InterpreterSpec, write the C interpreter and the Cython
    wrapper (generate a pyx and a pxd file).

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: testdir = tmp_dir()
        sage: rdf_interp = RDFInterpreter()
        sage: build_interp(rdf_interp, testdir)
        sage: open(testdir + '/interp_rdf.c').readline()
        '/* Automatically generated by ... */\n'
    """
    ig = InterpreterGenerator(interp_spec)
    interp_fn = '%s/interp_%s.c' % (dir, interp_spec.name)
    header_fn = '%s/interp_%s.h' % (dir, interp_spec.name)
    wrapper_fn = '%s/wrapper_%s.pyx' % (dir, interp_spec.name)
    pxd_fn = '%s/wrapper_%s.pxd' % (dir, interp_spec.name)
    interp = ig.get_interpreter()
    wrapper = ig.get_wrapper()
    pxd = ig.get_pxd()
    write_if_changed(interp_fn, interp)
    write_if_changed(wrapper_fn, wrapper)
    write_if_changed(pxd_fn, pxd)

def rebuild(dir, force=False):
    r"""
    Check whether the interpreter and wrapper sources have been written
    since the last time this module was changed.  If not, write them.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: testdir = tmp_dir()
        sage: rebuild(testdir)
        Building interpreters for fast_callable
        sage: open(testdir + '/wrapper_el.pyx').readline()
        '# Automatically generated by ...\n'
    """
    # This line will show up in "sage -b" (once per upgrade, not every time
    # you run it).
    print("Building interpreters for fast_callable")

    try:
        os.makedirs(dir)
    except OSError:
        pass

    # Although multiple files are generated by this function, since
    # they are all generated at once it suffices to make sure if just
    # one of the generated files is older than the generator sources
    src_file = __file__
    gen_file = os.path.join(dir, '__init__.py')

    if os.path.exists(gen_file):
        if getmtime(gen_file) > getmtime(src_file) and not force:
            return

    interp = RDFInterpreter()
    build_interp(interp, dir)

    interp = CDFInterpreter()
    build_interp(interp, dir)

    interp = RRInterpreter()
    build_interp(interp, dir)

    interp = PythonInterpreter()
    build_interp(interp, dir)

    interp = ElementInterpreter()
    build_interp(interp, dir)

    with open(os.path.join(dir, '__init__.py'), 'w') as f:
        f.write("# " + autogen_warn)
