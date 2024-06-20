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

"""Implements generic interpreter instructions and related utilities."""


import re

from .storage import ty_int


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
    given to the instruction, so len=None is different than len=1 even
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
        sage: from sage_setup.autogen.interpreters.specs.rdf import RDFInterpreter
        sage: pg = RDFInterpreter().pg
        sage: InstrSpec('add', pg('SS','S'), code='o0 = i0+i1;')
        add: SS->S = 'o0 = i0+i1;'
    """

    def __init__(self, name, io, code=None, uses_error_handler=False,
                 handles_own_decref=False):
        r"""
        Initialize an InstrSpec.

        INPUT:

        - name -- the name of the instruction
        - io -- a pair of lists of parameter specifications for I/O of the
          instruction
        - code -- a string containing a snippet of C code to read
          from the input variables and write to the output variables
        - uses_error_handler -- True if the instruction calls Python
          and jumps to error: on a Python error
        - handles_own_decref -- True if the instruction handles Python
          objects and includes its own
          reference-counting

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.specs.rdf import RDFInterpreter
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
                elif isinstance(len, int):
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
                elif isinstance(len, int):
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
            sage: from sage_setup.autogen.interpreters.specs.rdf import RDFInterpreter
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
        sage: from sage_setup.autogen.interpreters.specs.rdf import RDFInterpreter
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
        sage: from sage_setup.autogen.interpreters.specs.rdf import RDFInterpreter
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
        sage: from sage_setup.autogen.interpreters.specs.rdf import RDFInterpreter
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
        sage: from sage_setup.autogen.interpreters.specs.rr import RRInterpreter
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
        sage: from sage_setup.autogen.interpreters.specs.rr import RRInterpreter
        sage: pg = RRInterpreter().pg
        sage: instr_funcall_1arg_mpfr('exp', pg('S','S'), 'mpfr_exp')
        exp: S->S = 'mpfr_exp(o0, i0, MPFR_RNDN);'
    """
    return InstrSpec(name, io, code='%s(o0, i0, MPFR_RNDN);' % op)

def instr_funcall_2args_mpc(name, io, op):
    r"""
    A helper function for creating MPC instructions with two inputs
    and one output.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: from sage_setup.autogen.interpreters.specs.cc import CCInterpreter
        sage: pg = CCInterpreter().pg
        sage: instr_funcall_2args_mpc('add', pg('SS','S'), 'mpc_add')
        add: SS->S = 'mpc_add(o0, i0, i1, MPC_RNDNN);'
    """
    return InstrSpec(name, io, code='%s(o0, i0, i1, MPC_RNDNN);' % op)

def instr_funcall_1arg_mpc(name, io, op):
    r"""
    A helper function for creating MPC instructions with one input
    and one output.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters import *
        sage: from sage_setup.autogen.interpreters.specs.cc import CCInterpreter
        sage: pg = CCInterpreter().pg
        sage: instr_funcall_1arg_mpc('exp', pg('S','S'), 'mpc_exp')
        exp: S->S = 'mpc_exp(o0, i0, MPC_RNDNN);'
    """
    return InstrSpec(name, io, code='%s(o0, i0, MPC_RNDNN);' % op)
