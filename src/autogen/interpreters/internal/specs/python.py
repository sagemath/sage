# ****************************************************************************
#       Copyright (C) 2009 Carl Witty <Carl.Witty@gmail.com>
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from .base import StackInterpreter
from ..instructions import (params_gen, instr_funcall_2args, instr_unary,
                            InstrSpec)
from ..memory import MemoryChunk
from ..storage import ty_python
from ..utils import je, reindent_lines as ri


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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
        """
        return "    cdef int _n_%s\n" % self.name

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
            sage: mc.init_class_members()
            "        count = args['args']\n        self._n_args = count\n"
        """
        return je(ri(8,
            """
            count = args['{{ myself.name }}']
            self._n_args = count
            """), myself=self)

    def setup_args(self):
        r"""
        Handle the arguments of __call__.  Nothing to do.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
            sage: mc.setup_args()
            ''
        """
        return ''

    def pass_argument(self):
        r"""
        Pass the innards of the argument tuple to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPythonArguments('args', ty_python)
            sage: mc.pass_argument()
            '(<PyTupleObject*>args).ob_item'
        """
        return "(<PyTupleObject*>args).ob_item"


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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.name
            'domain'
            sage: mc.storage_type is ty_python
            True
        """
        super(MemoryChunkPyConstant, self).__init__(name, ty_python)

    def declare_class_members(self):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.declare_class_members()
            '    cdef object _domain\n'
        """
        return je(ri(4,
            """
            cdef object _{{ myself.name }}
            """), myself=self)

    def init_class_members(self):
        r"""
        Return a string to be put in the __init__ method of a wrapper
        class using this memory chunk, to initialize the corresponding
        class members.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.init_class_members()
            "        self._domain = args['domain']\n"
        """
        return je(ri(8,
            """
            self._{{ myself.name }} = args['{{ myself.name }}']
            """), myself=self)

    def declare_parameter(self):
        r"""
        Return the string to use to declare the interpreter parameter
        corresponding to this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
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

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
            sage: mc = MemoryChunkPyConstant('domain')
            sage: mc.pass_argument()
            '<PyObject*>self._domain'
        """
        return '<PyObject*>self._%s' % self.name


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

    name = 'py'

    def __init__(self):
        r"""
        Initialize a PythonInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.python import *
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

        super(PythonInterpreter, self).__init__(ty_python)
        # StackInterpreter.__init__ gave us a MemoryChunkArguments.
        # Override with MemoryChunkPythonArguments.
        self.mc_args = MemoryChunkPythonArguments('args', ty_python)
        self.chunks = [self.mc_args, self.mc_constants, self.mc_stack,
                       self.mc_code]
        self.c_header = ri(0,
            """
            #define CHECK(x) (x != NULL)
            """)

        self.pyx_header = ri(0,
            """\
            from cpython.number cimport PyNumber_TrueDivide
            """)

        pg = params_gen(A=self.mc_args, C=self.mc_constants, D=self.mc_code,
                        S=self.mc_stack)
        self.pg = pg

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
                       code=ri(0, """
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
                           """))
        ]

        binops = [
            ('add', 'PyNumber_Add'),
            ('sub', 'PyNumber_Subtract'),
            ('mul', 'PyNumber_Multiply'),
            ('div', 'PyNumber_TrueDivide'),
            ('floordiv', 'PyNumber_FloorDivide')
        ]

        for (name, op) in binops:
            instrs.append(instr_funcall_2args(name, pg('SS', 'S'), op))
        instrs.append(InstrSpec('pow', pg('SS', 'S'),
                                code='o0 = PyNumber_Power(i0, i1, Py_None);'))
        instrs.append(InstrSpec('ipow', pg('SC[D]', 'S'),
                                code='o0 = PyNumber_Power(i0, i1, Py_None);'))
        for (name, op) in [('neg', 'PyNumber_Negative'),
                           ('invert', 'PyNumber_Invert'),
                           ('abs', 'PyNumber_Absolute')]:
            instrs.append(instr_unary(name, pg('S', 'S'), '%s(i0)' % op))
        self.instr_descs = instrs
        self._set_opcodes()
        # Always use ipow
        self.ipow_range = True
        # We don't yet support call_c for Python-object interpreters
        # (the default implementation doesn't work, because of
        # object vs. PyObject* confusion)
        self.implement_call_c = False
