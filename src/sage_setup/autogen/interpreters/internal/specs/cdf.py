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

from __future__ import print_function, absolute_import

from .base import StackInterpreter
from ..instructions import (params_gen, instr_infix, instr_funcall_2args,
                            instr_unary, InstrSpec)
from ..memory import MemoryChunkConstants
from ..storage import ty_double_complex, ty_python
from ..utils import reindent_lines as ri


class CDFInterpreter(StackInterpreter):
    r"""
    A subclass of StackInterpreter, specifying an interpreter over
    complex machine-floating-point values (C doubles).
    """

    name = 'cdf'

    def __init__(self):
        r"""
        Initialize a CDFInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters.internal import *
            sage: from sage_setup.autogen.interpreters.internal.specs.cdf import *
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

        super(CDFInterpreter, self).__init__(ty_double_complex)
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
        self.c_header = ri(0,"""
            #include <stdlib.h>
            #include <complex.h>

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
            """)

        self.pxd_header = ri(0, """
            # This is to work around a header incompatibility with PARI using
            # "I" as variable conflicting with the complex "I".
            # If we cimport pari earlier, we avoid this problem.
            cimport cypari2.types

            # We need the type double_complex to work around
            #   http://trac.cython.org/ticket/869
            # so this is a bit hackish.
            cdef extern from "complex.h":
                ctypedef double double_complex "double complex"
            """)

        self.pyx_header = ri(0, """
            from sage.libs.gsl.complex cimport *
            from sage.rings.complex_double cimport ComplexDoubleElement
            import sage.rings.complex_double
            cdef object CDF = sage.rings.complex_double.CDF

            cdef extern from "complex.h":
                cdef double creal(double_complex)
                cdef double cimag(double_complex)
                cdef double_complex _Complex_I

            cdef inline double_complex CDE_to_dz(zz) noexcept:
                cdef ComplexDoubleElement z = <ComplexDoubleElement>(zz if isinstance(zz, ComplexDoubleElement) else CDF(zz))
                return GSL_REAL(z._complex) + _Complex_I * GSL_IMAG(z._complex)

            cdef inline ComplexDoubleElement dz_to_CDE(double_complex dz):
                cdef ComplexDoubleElement z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
                GSL_SET_COMPLEX(&z._complex, creal(dz), cimag(dz))
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

            """[1:])

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
                           ('mul', '*'), ('div', '/'),
                           ('truediv', '/')]:
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
