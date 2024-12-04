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
from ..instructions import (params_gen, instr_funcall_1arg_mpc,
                            instr_funcall_2args_mpc, InstrSpec)
from ..memory import MemoryChunk, MemoryChunkConstants
from ..storage import ty_mpc, ty_python
from ..utils import je, reindent_lines as ri


class MemoryChunkCCRetval(MemoryChunk):
    r"""
    A special-purpose memory chunk, for dealing with the return value
    of the CC-based interpreter.
    """

    def declare_class_members(self):
        r"""
        Return a string giving the declarations of the class members
        in a wrapper class for this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.cc import *
            sage: mc = MemoryChunkCCRetval('retval', ty_mpc)
            sage: mc.declare_class_members()
            ''
        """
        return ""

    def declare_call_locals(self):
        r"""
        Return a string to put in the __call__ method of a wrapper
        class using this memory chunk, to allocate local variables.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.cc import *
            sage: mc = MemoryChunkCCRetval('retval', ty_mpc)
            sage: mc.declare_call_locals()
            '        cdef ComplexNumber retval = (self.domain_element._new())\n'
        """
        return je(ri(8,
            """
            cdef ComplexNumber {{ myself.name }} = (self.domain_element._new())
            """), myself=self)

    def declare_parameter(self):
        r"""
        Return the string to use to declare the interpreter parameter
        corresponding to this memory chunk.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.cc import *
            sage: mc = MemoryChunkCCRetval('retval', ty_mpc)
            sage: mc.declare_parameter()
            'mpc_t retval'
        """
        return '%s %s' % (self.storage_type.c_reference_type(), self.name)

    def pass_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.cc import *
            sage: mc = MemoryChunkCCRetval('retval', ty_mpc)
            sage: mc.pass_argument()
            '(<mpc_t>(retval.__re))'
        """
        return je("""(<mpc_t>({{ myself.name }}.__re))""", myself=self)

    def pass_call_c_argument(self):
        r"""
        Return the string to pass the argument corresponding to this
        memory chunk to the interpreter, for use in the call_c method.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.cc import *
            sage: mc = MemoryChunkCCRetval('retval', ty_mpc)
            sage: mc.pass_call_c_argument()
            'result'
        """
        return "result"

class CCInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    MPFR arbitrary-precision floating-point numbers.
    """

    name = 'cc'

    def __init__(self):
        r"""
        Initialize a CCInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.cc import *
            sage: interp = CCInterpreter()
            sage: interp.name
            'cc'
            sage: interp.mc_py_constants
            {MC:py_constants}
            sage: interp.chunks
            [{MC:args}, {MC:retval}, {MC:constants}, {MC:py_constants}, {MC:stack}, {MC:code}, {MC:domain}]
            sage: interp.pg('A[D]', 'S')
            ([({MC:args}, {MC:code}, None)], [({MC:stack}, None, None)])
            sage: instrs = dict([(ins.name, ins) for ins in interp.instr_descs])
            sage: instrs['add']
            add: SS->S = 'mpc_add(o0, i0, i1, MPC_RNDNN);'
            sage: instrs['py_call']
            py_call: *->S = '\n  if (!cc_py_call...goto error;\n}\n'

        That py_call instruction is particularly interesting, and
        demonstrates a useful technique to let you use Cython code
        in an interpreter.  Let's look more closely::

            sage: print(instrs['py_call'].code)
            <BLANKLINE>
            if (!cc_py_call_helper(domain, i0, n_i1, i1, o0)) {
              goto error;
            }
            <BLANKLINE>

        This instruction makes use of the function cc_py_call_helper,
        which is declared::

            sage: print(interp.c_header)
            <BLANKLINE>
            #include <mpc.h>
            <BLANKLINE>

        So instructions where you need to interact with Python can
        call back into Cython code fairly easily.
        """

        mc_retval = MemoryChunkCCRetval('retval', ty_mpc)
        super(CCInterpreter, self).__init__(ty_mpc, mc_retval=mc_retval)
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
            #include <mpc.h>
            ''')

        self.pxd_header = ri(0,
            """
            from sage.rings.real_mpfr cimport RealNumber
            from sage.libs.mpfr cimport *
            from sage.rings.complex_mpfr cimport ComplexNumber
            from sage.libs.mpc cimport *
            """)

        self.pyx_header = ri(0,
            """\
            # distutils: libraries = mpfr mpc gmp

            cdef public bint cc_py_call_helper(object domain, object fn,
                                               int n_args,
                                               mpc_t* args, mpc_t retval) except 0:
                py_args = []
                cdef int i
                cdef ComplexNumber ZERO=domain.zero()
                cdef ComplexNumber cn
                for i from 0 <= i < n_args:
                    cn = ZERO._new()
                    mpfr_set(cn.__re, mpc_realref(args[i]), MPFR_RNDN)
                    mpfr_set(cn.__im, mpc_imagref(args[i]), MPFR_RNDN)
                    py_args.append(cn)
                cdef ComplexNumber result = domain(fn(*py_args))
                mpc_set_fr_fr(retval, result.__re,result.__im, MPC_RNDNN)
                return 1
            """)

        instrs = [
            InstrSpec('load_arg', pg('A[D]', 'S'),
                       code='mpc_set(o0, i0, MPC_RNDNN);'),
            InstrSpec('load_const', pg('C[D]', 'S'),
                       code='mpc_set(o0, i0, MPC_RNDNN);'),
            InstrSpec('return', pg('S', ''),
                       code='mpc_set(retval, i0, MPC_RNDNN);\nreturn 1;\n'),
            InstrSpec('py_call', pg('P[D]S@D', 'S'),
                       uses_error_handler=True,
                       code="""
  if (!cc_py_call_helper(domain, i0, n_i1, i1, o0)) {
  goto error;
}
""")
            ]
        for (name, op) in [('add', 'mpc_add'), ('sub', 'mpc_sub'),
                           ('mul', 'mpc_mul'), ('div', 'mpc_div'),
                           ('pow', 'mpc_pow')]:
            instrs.append(instr_funcall_2args_mpc(name, pg('SS', 'S'), op))
        instrs.append(instr_funcall_2args_mpc('ipow', pg('SD', 'S'), 'mpc_pow_si'))
        for name in ['neg',
                     'log', 'log10',
                     'exp',
                     'cos', 'sin', 'tan',
                     'acos', 'asin', 'atan',
                     'cosh', 'sinh', 'tanh',
                     'acosh', 'asinh', 'atanh']:
            instrs.append(instr_funcall_1arg_mpc(name, pg('S', 'S'), 'mpc_' + name))
        # mpc_ui_div constructs a temporary mpc_t and then calls mpc_div;
        # it would probably be (slightly) faster to use a permanent copy
        # of "one" (on the other hand, the constructed temporary copy is
        # on the stack, so it's very likely to be in the cache).
        instrs.append(InstrSpec('invert', pg('S', 'S'),
                             code='mpc_ui_div(o0, 1, i0, MPC_RNDNN);'))
        self.instr_descs = instrs
        self._set_opcodes()
        # Supported for exponents that fit in a long, so we could use
        # a much wider range on a 64-bit machine.  On the other hand,
        # it's easier to write the code this way, and constant integer
        # exponents outside this range probably aren't very common anyway.
        self.ipow_range = (int(-2**31), int(2**31-1))
