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


from .base import StackInterpreter
from ..instructions import (params_gen, instr_infix, instr_funcall_2args,
                            instr_unary, InstrSpec)
from ..memory import MemoryChunkConstants
from ..storage import ty_double, ty_python
from ..utils import reindent_lines as ri


class RDFInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    machine-floating-point values (C doubles).  This is used for
    both domain=RDF and domain=float; currently the only difference
    between the two is the type of the value returned from the
    wrapper (they use the same wrapper and interpreter).
    """

    name = 'rdf'

    def __init__(self):
        r"""
        Initialize an RDFInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.specs.rdf import *
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

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=('x','y'))
            sage: x = etb.var('x')
            sage: y = etb.var('y')
            sage: ff = fast_callable(x^y, domain=RDF)
            sage: ff(1.5, 3)
            3.375
            sage: ff(-2, 3)
            -8.0
            sage: ff(-2, 1/3)
            Traceback (most recent call last):
            ...
            ValueError: negative number to a fractional power not real
        """

        super(RDFInterpreter, self).__init__(ty_double)
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
                       code=ri(0, """
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
                           """)),
            InstrSpec('pow', pg('SS', 'S'),
                       uses_error_handler=True,
                       code=ri(0, """
                           /* See python's pow in floatobject.c */
                           if (i0 == 0) o0 = 1.0;
                           else {
                             if (i0 < 0 && i1 != floor(i1)) {
                                 PyErr_SetString(PyExc_ValueError, "negative number to a fractional power not real");
                                 goto error;
                             }
                             o0 = pow(i0, i1);
                           }
                           """))
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
