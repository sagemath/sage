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
from .python import (MemoryChunkPyConstant, MemoryChunkPythonArguments,
                     PythonInterpreter)
from ..storage import ty_python
from ..utils import reindent_lines as ri


class MemoryChunkElementArguments(MemoryChunkPythonArguments):
    r"""
    A special-purpose memory chunk, for the Python-object based
    interpreters that want to process (and perhaps modify) the data.

    We allocate a new list on every call to hold the modified arguments.
    That's not strictly necessary -- we could pre-allocate a list and map into
    it -- but this lets us use simpler code for a very-likely-negligible
    efficiency cost.  (The Element interpreter is going to allocate lots of
    objects as it runs, anyway.)
    """

    def setup_args(self):
        r"""
        Handle the arguments of __call__.  Note: This hardcodes
        "self._domain".

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.internal.specs.element import *
            sage: mc = MemoryChunkElementArguments('args', ty_python)
            sage: mc.setup_args()
            'mapped_args = [self._domain(a) for a in args]\n'
        """
        return "mapped_args = [self._domain(a) for a in args]\n"

    def pass_argument(self):
        r"""
        Pass the innards of the argument tuple to the interpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.internal.specs.element import *
            sage: mc = MemoryChunkElementArguments('args', ty_python)
            sage: mc.pass_argument()
            '(<PyListObject*>mapped_args).ob_item'
        """
        return "(<PyListObject*>mapped_args).ob_item"


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

    name = 'el'

    def __init__(self):
        r"""
        Initialize an ElementInterpreter.

        EXAMPLES::

            sage: from sage_setup.autogen.interpreters import *
            sage: from sage_setup.autogen.interpreters.internal.specs.element import *
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

        super(ElementInterpreter, self).__init__()
        # PythonInterpreter.__init__ gave us a MemoryChunkPythonArguments.
        # Override with MemoryChunkElementArguments.
        self.mc_args = MemoryChunkElementArguments('args', ty_python)
        self.mc_domain_info = MemoryChunkPyConstant('domain')
        self.chunks = [self.mc_args, self.mc_constants, self.mc_stack,
                       self.mc_domain_info, self.mc_code]
        self.c_header = ri(0, """
            #include "wrapper_el.h"

            #define CHECK(x) do_check(&(x), domain)

            static inline int do_check(PyObject **x, PyObject *domain) {
              if (*x == NULL) return 0;
              PyObject *new_x = el_check_element(*x, domain);
              Py_DECREF(*x);
              *x = new_x;
              if (*x == NULL) return 0;
              return 1;
            }
            """)

        self.pyx_header += ri(0, """
            from sage.structure.element cimport Element

            cdef public object el_check_element(object v, parent):
                cdef Element v_el

                if isinstance(v, Element):
                    v_el = <Element>v
                    if v_el._parent is parent:
                        return v_el

                return parent(v)

            """[1:])
