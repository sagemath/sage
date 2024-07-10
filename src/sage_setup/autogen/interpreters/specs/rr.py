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
from .python import MemoryChunkPyConstant
from ..instructions import (params_gen, instr_funcall_1arg_mpfr,
                            instr_funcall_2args_mpfr, InstrSpec)
from ..memory import MemoryChunk, MemoryChunkConstants
from ..storage import ty_mpfr, ty_python
from ..utils import je, reindent_lines as ri


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
            sage: from sage_setup.autogen.interpreters.specs.rr import *
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
            sage: from sage_setup.autogen.interpreters.specs.rr import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.declare_call_locals()
            '        cdef RealNumber retval = (self.domain)()\n'
        """
        return je(ri(8,
            """
            cdef RealNumber {{ myself.name }} = (self.domain)()
            """), myself=self)

    def declare_parameter(self):
        r"""
        Return the string to use to declare the interpreter parameter
        corresponding to this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.specs.rr import *
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
            sage: from sage_setup.autogen.interpreters.specs.rr import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.pass_argument()
            'retval.value'
        """
        return je("""{{ myself.name }}.value""", myself=self)

    def pass_call_c_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter, for use in the call_c method.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.specs.rr import *
            sage: mc = MemoryChunkRRRetval('retval', ty_mpfr)
            sage: mc.pass_call_c_argument()
            'result'
        """
        return "result"


class RRInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    MPFR arbitrary-precision floating-point numbers.
    """

    name = 'rr'

    def __init__(self):
        r"""
        Initialize an RDFInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.specs.rr import *
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
            #include "wrapper_rr.h"
            <BLANKLINE>

        The function ``rr_py_call_helper`` is implemented in Cython::

            sage: print(interp.pyx_header)
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

        mc_retval = MemoryChunkRRRetval('retval', ty_mpfr)
        super(RRInterpreter, self).__init__(ty_mpfr, mc_retval=mc_retval)
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
        self.c_header = ri(0,
            '''
            #include <mpfr.h>
            #include "wrapper_rr.h"
            ''')

        self.pxd_header = ri(0,
            """
            from sage.rings.real_mpfr cimport RealField_class, RealNumber
            from sage.libs.mpfr cimport *

            """)

        self.pyx_header = ri(0,
            """\
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

            """)

        instrs = [
            InstrSpec('load_arg', pg('A[D]', 'S'),
                       code='mpfr_set(o0, i0, MPFR_RNDN);'),
            InstrSpec('load_const', pg('C[D]', 'S'),
                       code='mpfr_set(o0, i0, MPFR_RNDN);'),
            InstrSpec('return', pg('S', ''),
                       code='mpfr_set(retval, i0, MPFR_RNDN);\nreturn 1;\n'),
            InstrSpec('py_call', pg('P[D]S@D', 'S'),
                       uses_error_handler=True,
                       code=ri(0,
                           """
                           if (!rr_py_call_helper(domain, i0, n_i1, i1, o0)) {
                             goto error;
                           }
                           """))
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
