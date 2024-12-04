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

"""General purpose MemoryChunk types and related utilities"""

from __future__ import print_function, absolute_import

from .utils import je, reindent_lines as ri


def string_of_addr(a):
    r"""
    An address or a length from a parameter specification may be
    either None, an integer, or a MemoryChunk.  If the address or
    length is an integer or a MemoryChunk, this function will convert
    it to a string giving an expression that will evaluate to the correct
    address or length.  (See the docstring for params_gen for more
    information on parameter specifications.)

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters.internal import *
        sage: mc_code = MemoryChunkConstants('code', ty_int)
        sage: string_of_addr(mc_code)
        '*code++'
        sage: string_of_addr(42r)
        '42'
    """
    if isinstance(a, int):
        return str(a)
    assert(isinstance(a, MemoryChunk))
    return '*%s++' % a.name


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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: mc.declare_class_members()
            '    cdef int _n_args\n    cdef mpfr_t* _args\n'
        """
        return self.storage_type.declare_chunk_class_members(self.name)

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: print(mc.init_class_members())
                    count = args['args']
                    self._n_args = count
                    self._args = <mpfr_t*>check_allocarray(self._n_args, sizeof(mpfr_t))
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.rr import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.declare_call_locals()
            '        cdef RealNumber retval = (self.domain)()\n'
        """
        return ""

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: mc = MemoryChunkArguments('args', ty_double)
            sage: print(mc.init_class_members())
                    count = args['args']
                    self._n_args = count
                    self._args = <double*>check_allocarray(self._n_args, sizeof(double))
            <BLANKLINE>
        """
        return je(ri(0, """
                    count = args['{{ myself.name }}']
            {% print(myself.storage_type.alloc_chunk_data(myself.name, 'count')) %}
            """), myself=self)

    def dealloc_class_members(self):
        r"""
        Return a string to be put in the __dealloc__ method of a wrapper
        class using this memory chunk, to deallocate the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: mc = MemoryChunkConstants('constants', ty_mpfr)
            sage: print(mc.init_class_members())
                    val = args['constants']
                    self._n_constants = len(val)
                    self._constants = <mpfr_t*>check_allocarray(self._n_constants, sizeof(mpfr_t))
                    for i in range(len(val)):
                        mpfr_init2(self._constants[i], self.domain.prec())
                    for i in range(len(val)):
                        rn = self.domain(val[i])
                        mpfr_set(self._constants[i], rn.value, MPFR_RNDN)
            <BLANKLINE>
        """
        return je(ri(0, """
                    val = args['{{ myself.name }}']
            {% print(myself.storage_type.alloc_chunk_data(myself.name, 'len(val)')) %}
                    for i in range(len(val)):
                        {{ myself.storage_type.assign_c_from_py('self._%s[i]' % myself.name, 'val[i]') | i(12) }}
            """), myself=self)


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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: mc = MemoryChunkArguments('args', ty_mpfr)
            sage: print(mc.setup_args())
            cdef mpfr_t* c_args = self._args
            cdef int i
            for i from 0 <= i < len(args):
                rn = self.domain(args[i])
                mpfr_set(self._args[i], rn.value, MPFR_RNDN)
            <BLANKLINE>
        """
        return je(ri(0, """
            cdef {{ myself.storage_type.c_ptr_type() }} c_args = self._args
            cdef int i
            for i from 0 <= i < len(args):
                {{ myself.storage_type.assign_c_from_py('self._args[i]', 'args[i]') | i(4) }}
            """), myself=self)

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: mc = MemoryChunkScratch('stack', ty_double, is_stack=True)
            sage: mc.name
            'stack'
            sage: mc.storage_type is ty_double
            True
            sage: mc._is_stack
            True
        """

        super(MemoryChunkScratch, self).__init__(name, storage_type)
        self._is_stack = is_stack

    def is_stack(self):
        r"""
        Says whether this memory chunk is a stack.  This affects code
        generation for instructions using this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: mc = MemoryChunkScratch('registers', ty_python)
            sage: print(mc.handle_cleanup())
            for i in range(self._n_registers):
                Py_CLEAR(self._registers[i])
            <BLANKLINE>
        """
        # XXX This is a lot slower than it needs to be, because
        # we don't have a "cdef int i" in scope here.
        return je(ri(0, """
            for i in range(self._n_{{ myself.name }}):
                Py_CLEAR(self._{{ myself.name }}[i])
            """), myself=self)
