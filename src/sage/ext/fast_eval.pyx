r"""
Fast Numerical Evaluation

For many applications such as numerical integration, differential
equation approximation, plotting a 3d surface, optimization problems,
monte-carlo simulations, etc., one wishes to pass around and evaluate
a single algebraic expression many, many times at various floating
point values. Doing this via recursive calls over a python
representation of the object (even if Maxima or other outside packages
are not involved) is extremely inefficient.

Up until now the solution has been to use lambda expressions, but this
is neither intuitive, Sage-like, nor efficient (compared to operating
on raw C doubles).  This module provides a representation of algebraic
expression in Reverse Polish Notation, and provides an efficient
interpreter on C double values as a callable python object. It does
what it can in C, and will call out to Python if necessary.

Essential to the understanding of this class is the distinction
between symbolic expressions and callable symbolic expressions (where
the latter binds argument names to argument positions). The
``*vars`` parameter passed around encapsulates this information.

See the function ``fast_float(f, *vars)`` to create a fast-callable
version of f.

.. NOTE::

    Sage temporarily has two implementations of this functionality ;
    one in this file, which will probably be deprecated soon, and one in
    fast_callable.pyx.  The following instructions are for the old
    implementation; you probably want to be looking at fast_callable.pyx
    instead.

To provide this interface for a class, implement ``fast_float_(self, *vars)``.  The basic building blocks are
provided by the functions ``fast_float_constant`` (returns a
constant function), ``fast_float_arg`` (selects the ``n``-th value
when called with ``\ge_n`` arguments), and ``fast_float_func`` which
wraps a callable Python function. These may be combined with the
standard Python arithmetic operators, and support many of the basic
math functions such ``sqrt``, ``exp``, and trig functions.

EXAMPLES::

    sage: from sage.ext.fast_eval import fast_float
    sage: f = fast_float(sqrt(x^7+1), 'x', old=True)
    sage: f(1)
    1.4142135623730951
    sage: f.op_list()
    ['load 0', 'push 7.0', 'pow', 'push 1.0', 'add', 'call sqrt(1)']

To interpret that last line, we load argument 0 (``x`` in this case) onto
the stack, push the constant 2.0 onto the stack, call the pow function
(which takes 2 arguments from the stack), push the constant 1.0, add the
top two arguments of the stack, and then call sqrt.

Here we take ``sin`` of the first argument and add it to ``f``::

    sage: from sage.ext.fast_eval import fast_float_arg
    sage: g = fast_float_arg(0).sin()
    sage: (f+g).op_list()
    ['load 0', 'push 7.0', 'pow', 'push 1.0', 'add', 'call sqrt(1)', 'load 0', 'call sin(1)', 'add']

TESTS:

This used to segfault because of an assumption that assigning None to a
variable would raise a TypeError::

    sage: from sage.ext.fast_eval import fast_float_arg, fast_float
    sage: fast_float_arg(0)+None
    Traceback (most recent call last):
    ...
    TypeError

AUTHORS:

- Robert Bradshaw (2008-10): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.memory cimport sig_malloc, sig_free

from sage.ext.fast_callable import fast_callable, Wrapper
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool

cimport cython
from cpython.ref cimport Py_INCREF
from cpython.object cimport PyObject_CallObject
from cpython.int cimport PyInt_AS_LONG
from cpython.tuple cimport PyTuple_New, PyTuple_SET_ITEM


cdef extern from "math.h":
    double sqrt(double)
    double pow(double, double)

    double ceil(double)
    double floor(double)

    double sin(double)
    double cos(double)
    double tan(double)

    double asin(double)
    double acos(double)
    double atan(double)
    double atan2(double, double)

    double sinh(double)
    double cosh(double)
    double tanh(double)

    double asinh(double)
    double acosh(double)
    double atanh(double)

    double exp(double)
    double log(double)
    double log10(double)
    double log2_ "log2"(double)


# This is only needed on Cygwin since log2 is a macro.
# If we don't do this the cygwin GCC gets very confused.
cdef inline double log2(double x):
    return log2_(x)

cdef extern from *:
    void* memcpy(void* dst, void* src, size_t len)

cdef inline int max(int a, int b):
    return a if a > b else b

cdef inline int min(int a, int b):
    return a if a < b else b

cdef enum:
# stack
    LOAD_ARG    # push input argument n onto the stack
    PUSH_CONST
    POP
    POP_N
    DUP

# basic arithmetic
    ADD
    SUB
    MUL
    DIV
    NEG
    ABS
    INVERT
    POW

# basic comparison
    LT
    LE
    EQ
    NE
    GT
    GE

# functional
    ONE_ARG_FUNC
    TWO_ARG_FUNC
    PY_FUNC


# These two dictionaries are for printable and machine independent representation.

op_names = {
    LOAD_ARG: 'load',
    PUSH_CONST: 'push',
    POP: 'pop',
    POP_N: 'popn',
    DUP: 'dup',

    ADD: 'add',
    SUB: 'sub',
    MUL: 'mul',
    DIV: 'div',
    NEG: 'neg',
    ABS: 'abs',
    INVERT: 'invert',
    POW: 'pow',


    LT: 'lt',
    LE: 'le',
    EQ: 'eq',
    NE: 'ne',
    GT: 'gt',
    GE: 'ge',


    ONE_ARG_FUNC: 'call',
    TWO_ARG_FUNC: 'call',
    PY_FUNC: 'py_call',
}

cfunc_names = {
    <size_t>&sqrt: 'sqrt',
    <size_t>&pow: 'pow',

    <size_t>&ceil: 'ceil',
    <size_t>&floor: 'floor',

    <size_t>&sin: 'sin',
    <size_t>&cos: 'cos',
    <size_t>&tan: 'tan',

    <size_t>&asin: 'asin',
    <size_t>&atan: 'atan',
    <size_t>&atan2: 'atan2',

    <size_t>&sinh: 'sinh',
    <size_t>&cosh: 'cosh',
    <size_t>&tanh: 'tanh',

    <size_t>&asinh: 'asinh',
    <size_t>&acosh: 'acosh',
    <size_t>&atanh: 'atanh',

    <size_t>&exp: 'exp',
    <size_t>&log: 'log',
    <size_t>&log2: 'log2',
    <size_t>&log10: 'log10',

}

cdef reverse_map(m):
    r = {}
    for key, value in m.iteritems():
        r[value] = key
    return r

# With all the functionality around the op struct, perhaps there should be
# a wrapper class, though we still wish to operate on pure structs for speed.

cdef op_to_string(fast_double_op op):
    s = op_names[op.type]
    if op.type in [LOAD_ARG, POP_N]:
        s += " %s" % op.params.n
    elif op.type == PUSH_CONST:
        s += " %s" % op.params.c
    elif op.type in [ONE_ARG_FUNC, TWO_ARG_FUNC]:
        try:
            cname = cfunc_names[<size_t>op.params.func]
        except KeyError:
            cname = "0x%x" % <size_t>op.params.func
        s += " %s(%s)" % (cname, 1 if op.type == ONE_ARG_FUNC else 2)
    elif op.type == PY_FUNC:
        n, func = <object>(op.params.func)
        s += " %s(%s)" % (func, n)
    return s

cdef op_to_tuple(fast_double_op op):
    s = op_names[op.type]
    if op.type in [LOAD_ARG, POP_N]:
        param = op.params.n
    elif op.type == PUSH_CONST:
        param = op.params.c
    elif op.type in [ONE_ARG_FUNC, TWO_ARG_FUNC]:
        param_count = 1 if op.type == ONE_ARG_FUNC else 2
        try:
            param = param_count, cfunc_names[<size_t>op.params.func]
        except KeyError:
            raise ValueError("Unknown C function: 0x%x"
                             % <size_t>op.params.func)
    elif op.type == PY_FUNC:
        param = <object>(op.params.func)
    else:
        param = None
    if param is None:
        return (s,)
    else:
        return s, param

def _unpickle_FastDoubleFunc(nargs, max_height, op_list):
    cdef FastDoubleFunc self = FastDoubleFunc.__new__(FastDoubleFunc)
    self.nops = len(op_list)
    self.nargs = nargs
    self.max_height = max_height
    self.ops = <fast_double_op *>sig_malloc(sizeof(fast_double_op) * self.nops)
    self.allocate_stack()
    cfunc_addresses = reverse_map(cfunc_names)
    op_enums = reverse_map(op_names)
    cdef size_t address
    cdef int i = 0, type
    for op in op_list:
        self.ops[i].type = type = op_enums[op[0]]
        if type in [LOAD_ARG, POP_N]:
            self.ops[i].params.n = op[1]
        elif type == PUSH_CONST:
            self.ops[i].params.c = op[1]
        elif type in [ONE_ARG_FUNC, TWO_ARG_FUNC]:
            param_count, cfunc = op[1]
            address = cfunc_addresses[cfunc]
            self.ops[i].params.func = <PyObject*>address
            self.ops[i].type = ['', ONE_ARG_FUNC, TWO_ARG_FUNC][param_count]
        elif type == PY_FUNC:
            if self.py_funcs is None:
                self.py_funcs = op[1]
            else:
                self.py_funcs = self.py_funcs + (op[1],)
            self.ops[i].params.func = <PyObject*>op[1]
        i += 1
    return self


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int process_op(fast_double_op op, double* stack, double* argv, int top) except -2:
    cdef int i, n
    cdef object arg
    cdef tuple py_args

    if op.type == LOAD_ARG:
        stack[top+1] = argv[op.params.n]
        return top+1

    elif op.type == PUSH_CONST:
        stack[top+1] = op.params.c
        return top+1

    elif op.type == POP:
        return top-1

    elif op.type == POP_N:
        return top-op.params.n

    elif op.type == DUP:
        stack[top+1] = stack[top]
        return top+1

    elif op.type == ADD:
        stack[top-1] += stack[top]
        return top-1

    elif op.type == SUB:
        stack[top-1] -= stack[top]
        return top-1

    elif op.type == MUL:
        stack[top-1] *= stack[top]
        return top-1

    elif op.type == DIV:
        stack[top-1] /= stack[top]
        return top-1

    elif op.type == NEG:
        stack[top] = -stack[top]
        return top

    elif op.type == ABS:
        if stack[top] < 0:
            stack[top] = -stack[top]
        return top

    elif op.type == INVERT:
        stack[top] = 1/stack[top]
        return top

    elif op.type == POW:
        if stack[top-1] < 0 and stack[top] != floor(stack[top]):
            raise ValueError("negative number to a fractional power not real")
        stack[top-1] = pow(stack[top-1], stack[top])
        return top-1

    elif op.type == LT:
        stack[top-1] = 1.0 if stack[top-1] < stack[top] else 0.0
        return top-1

    elif op.type == LE:
        stack[top-1] = 1.0 if stack[top-1] <= stack[top] else 0.0
        return top-1

    elif op.type == EQ:
        stack[top-1] = 1.0 if stack[top-1] == stack[top] else 0.0
        return top-1

    elif op.type == NE:
        stack[top-1] = 1.0 if stack[top-1] != stack[top] else 0.0
        return top-1

    elif op.type == GT:
        stack[top-1] = 1.0 if stack[top-1] > stack[top] else 0.0
        return top-1

    elif op.type == GE:
        stack[top-1] = 1.0 if stack[top-1] >= stack[top] else 0.0
        return top-1

    elif op.type == ONE_ARG_FUNC:
        stack[top] = (op.params.f)(stack[top])
        return top

    elif op.type == TWO_ARG_FUNC:
        stack[top-1] = (op.params.ff)(stack[top-1], stack[top])
        return top-1

    elif op.type == PY_FUNC:
        # We use a few direct C/API calls here because Cython itself
        # doesn't generate optimal code for this.
        n = PyInt_AS_LONG((<tuple>op.params.func)[0])
        top = top - n + 1
        py_args = PyTuple_New(n)
        for i in range(n):
            arg = stack[top+i]
            Py_INCREF(arg)  # PyTuple_SET_ITEM() steals a reference
            PyTuple_SET_ITEM(py_args, i, arg)
        stack[top] = PyObject_CallObject((<tuple>op.params.func)[1], py_args)
        return top

    raise RuntimeError("Bad op code %s" % op.type)


cdef class FastDoubleFunc:
    """
    This class is for fast evaluation of algebraic expressions over
    the real numbers (e.g. for plotting). It represents an expression
    as a stack-based series of operations.

    EXAMPLES::

        sage: from sage.ext.fast_eval import FastDoubleFunc
        sage: f = FastDoubleFunc('const', 1.5) # the constant function
        sage: f()
        1.5
        sage: g = FastDoubleFunc('arg', 0) # the first argument
        sage: g(5)
        5.0
        sage: h = f+g
        sage: h(17)
        18.5
        sage: h = h.sin()
        sage: h(pi/2-1.5)
        1.0
        sage: h.is_pure_c()
        True
        sage: list(h)
        ['push 1.5', 'load 0', 'add', 'call sin(1)']

    We can wrap Python functions too::

        sage: h = FastDoubleFunc('callable', lambda x,y: x*x*x - y, g, f)
        sage: h(10)
        998.5
        sage: h.is_pure_c()
        False
        sage: list(h)
        ['load 0', 'push 1.5', 'py_call <function <lambda> at 0x...>(2)']

    Here's a more complicated expression::

        sage: from sage.ext.fast_eval import fast_float_constant, fast_float_arg
        sage: a = fast_float_constant(1.5)
        sage: b = fast_float_constant(3.14)
        sage: c = fast_float_constant(7)
        sage: x = fast_float_arg(0)
        sage: y = fast_float_arg(1)
        sage: f = a*x^2 + b*x + c - y/sqrt(sin(y)^2+a)
        sage: f(2,3)
        16.846610528508116
        sage: f.max_height
        4
        sage: f.is_pure_c()
        True
        sage: list(f)
        ['push 1.5', 'load 0', 'dup', 'mul', 'mul', 'push 3.14', 'load 0', 'mul', 'add', 'push 7.0', 'add', 'load 1', 'load 1', 'call sin(1)', 'dup', 'mul', 'push 1.5', 'add', 'call sqrt(1)', 'div', 'sub']

    AUTHORS:

    - Robert Bradshaw
    """
    def __init__(self, type, param, *args):

        cdef FastDoubleFunc arg
        cdef int i

        if type == 'arg':
            self.nargs = param+1
            self.nops = 1
            self.max_height = 1
            self.ops = <fast_double_op *>sig_malloc(sizeof(fast_double_op))
            self.ops[0].type = LOAD_ARG
            self.ops[0].params.n = param

        elif type == 'const':
            self.nargs = 0
            self.nops = 1
            self.max_height = 1
            self.ops = <fast_double_op *>sig_malloc(sizeof(fast_double_op))
            self.ops[0].type = PUSH_CONST
            self.ops[0].params.c = param

        elif type == 'callable':
            py_func = len(args), param
            self.py_funcs = (py_func,) # just so it doesn't get garbage collected
            self.nops = 1
            self.nargs = 0
            for i from 0 <= i < len(args):
                a = args[i]
                if not isinstance(a, FastDoubleFunc):
                     a = FastDoubleFunc('const', a)
                     args = args[:i] + (a,) + args[i+1:]
                arg = a
                self.nops += arg.nops
                if arg.py_funcs is not None:
                    self.py_funcs += arg.py_funcs
                self.nargs = max(self.nargs, arg.nargs)
                self.max_height = max(self.max_height, arg.max_height+i)
            self.ops = <fast_double_op *>sig_malloc(sizeof(fast_double_op) * self.nops)
            if self.ops == NULL:
                raise MemoryError
            i = 0
            for arg in args:
                memcpy(self.ops + i, arg.ops, sizeof(fast_double_op) * arg.nops)
                i += arg.nops
            self.ops[self.nops-1].type = PY_FUNC
            self.ops[self.nops-1].params.func = <PyObject*>py_func

        else:
            raise ValueError("Unknown operation: %s" % type)

        self.allocate_stack()

    cdef int allocate_stack(FastDoubleFunc self) except -1:
        self.argv = <double*>sig_malloc(sizeof(double) * self.nargs)
        if self.argv == NULL:
            raise MemoryError
        self.stack = <double*>sig_malloc(sizeof(double) * self.max_height)
        if self.stack == NULL:
            raise MemoryError

    def __dealloc__(self):
        sig_free(self.ops)
        sig_free(self.stack)
        sig_free(self.argv)

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.ext.fast_eval import fast_float_arg, fast_float_func
            sage: f = fast_float_arg(0).sin() * 10 + fast_float_func(hash, fast_float_arg(1))
            sage: loads(dumps(f)) == f
            True
        """
        L = [op_to_tuple(self.ops[i]) for i from 0 <= i < self.nops]
        return _unpickle_FastDoubleFunc, (self.nargs, self.max_height, L)

    def __richcmp__(self, other, op):
        """
        Two functions are considered equal if they represent the same
        exact sequence of operations.

        TESTS::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: fast_float_arg(0) == fast_float_arg(0)
            True
            sage: fast_float_arg(0) == fast_float_arg(1)
            False
            sage: fast_float_arg(0) == fast_float_arg(0).sin()
            False
        """
        cdef int c, i
        cdef FastDoubleFunc left, right
        try:
            left = <FastDoubleFunc?>self
            right = <FastDoubleFunc?>other

            lx = left.nargs
            rx = right.nargs
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

            lx = left.nops
            rx = right.nops
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

            lx = left.max_height
            rx = right.max_height
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

            for i from 0 <= i < self.nops:
                lx = left.ops[i].type
                rx = right.ops[i].type
                if lx != rx:
                    return richcmp_not_equal(lx, rx, op)

            for i from 0 <= i < self.nops:
                lx = op_to_tuple(left.ops[i])
                rx = op_to_tuple(right.ops[i])
                if lx != rx:
                    return richcmp_not_equal(lx, rx, op)

            return rich_to_bool(op, 0)
        except TypeError:
            return NotImplemented

    def __call__(FastDoubleFunc self, *args):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(2)
            sage: f(0,1,2,3)
            2.0
            sage: f(10)
            Traceback (most recent call last):
            ...
            TypeError: Wrong number of arguments (need at least 3, got 1)
            sage: f('blah', 1, 2, 3) # py2
            Traceback (most recent call last):
            ...
            TypeError: a float is required
            sage: f('blah', 1, 2, 3) # py3
            Traceback (most recent call last):
            ...
            TypeError: must be real number, not str
        """
        if len(args) < self.nargs:
            raise TypeError("Wrong number of arguments (need at least %s, got %s)" % (self.nargs, len(args)))
        cdef int i = 0
        for i from 0 <= i < self.nargs:
            self.argv[i] = args[i]
        res = self._call_c(self.argv)
        return res

    cdef double _call_c(FastDoubleFunc self, double* argv) except? -2:
        # The caller must assure that argv has length at least self.nargs
        # The bulk of this function is in the (inlined) function process_op.
        cdef int i, top = -1
        for i from 0 <= i < self.nops:
            top = process_op(self.ops[i], self.stack, argv, top)
        cdef double res = self.stack[0]
        return res

    def _fast_float_(self, *vars):
        r"""
        Returns ``self`` if there are enough arguments, otherwise raises a ``TypeError``.

        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(1)
            sage: f._fast_float_('x','y') is f
            True
            sage: f._fast_float_('x') is f
            Traceback (most recent call last):
            ...
            TypeError: Needs at least 2 arguments (1 provided)
        """
        if self.nargs > len(vars):
            raise TypeError("Needs at least %s arguments (%s provided)" % (self.nargs, len(vars)))
        return self

    def op_list(self):
        """
        Returns a list of string representations of the
        operations that make up this expression.

        Python and C function calls may be only available by function
        pointer addresses.

        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_constant, fast_float_arg
            sage: a = fast_float_constant(17)
            sage: x = fast_float_arg(0)
            sage: a.op_list()
            ['push 17.0']
            sage: x.op_list()
            ['load 0']
            sage: (a*x).op_list()
            ['push 17.0', 'load 0', 'mul']
            sage: (a+a*x^2).sqrt().op_list()
            ['push 17.0', 'push 17.0', 'load 0', 'dup', 'mul', 'mul', 'add', 'call sqrt(1)']
        """
        cdef int i
        return [op_to_string(self.ops[i]) for i from 0 <= i < self.nops]

    def __iter__(self):
        """
        Returns the list of operations of ``self``.

        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0)*2 + 3
            sage: list(f)
            ['load 0', 'push 2.0', 'mul', 'push 3.0', 'add']
        """
        return iter(self.op_list())

    cpdef bint is_pure_c(self):
        """
        Returns ``True`` if this function can be evaluated without
        any python calls (at any level).

        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_constant, fast_float_arg, fast_float_func
            sage: fast_float_constant(2).is_pure_c()
            True
            sage: fast_float_arg(2).sqrt().sin().is_pure_c()
            True
            sage: fast_float_func(lambda _: 2).is_pure_c()
            False
        """
        cdef int i
        for i from 0 <= i < self.nops:
            if self.ops[i].type == PY_FUNC:
                return 0
        return 1

    def python_calls(self):
        """
        Returns a list of all python calls used by function.

        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_func, fast_float_arg
            sage: x = fast_float_arg(0)
            sage: f = fast_float_func(hash, sqrt(x))
            sage: f.op_list()
            ['load 0', 'call sqrt(1)', 'py_call <built-in function hash>(1)']
            sage: f.python_calls()
            [<built-in function hash>]
        """
        L = []
        cdef int i
        for i from 0 <= i < self.nops:
            if self.ops[i].type == PY_FUNC:
                L.append((<object>self.ops[i].params.func)[1])
        return L

    ###################################################################
    #   Basic Arithmetic
    ###################################################################

    def __add__(left, right):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0) + fast_float_arg(1)
            sage: f(3,4)
            7.0
        """
        return binop(left, right, ADD)

    def __sub__(left, right):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0) - fast_float_arg(2)
            sage: f(3,4,5)
            -2.0
        """
        return binop(left, right, SUB)

    def __mul__(left, right):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0) * 2
            sage: f(17)
            34.0
        """
        return binop(left, right, MUL)

    def __truediv__(left, right):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).__truediv__(7)
            sage: f(14)
            2.0
        """
        return binop(left, right, DIV)

    def __pow__(FastDoubleFunc left, right, dummy):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import FastDoubleFunc
            sage: f = FastDoubleFunc('arg', 0)^2
            sage: f(2)
            4.0
            sage: f = FastDoubleFunc('arg', 0)^4
            sage: f(2)
            16.0
            sage: f = FastDoubleFunc('arg', 0)^-3
            sage: f(2)
            0.125
            sage: f = FastDoubleFunc('arg', 0)^FastDoubleFunc('arg', 1)
            sage: f(5,3)
            125.0

        TESTS::

            sage: var('a,b')
            (a, b)
            sage: ff = (a^b)._fast_float_(a,b)
            sage: ff(2, 9)
            512.0
            sage: ff(-2, 9)
            -512.0
            sage: ff(-2, 9.1)
            Traceback (most recent call last):
            ...
            ValueError: negative number to a fractional power not real
        """
        if isinstance(right, FastDoubleFunc) and right.nargs == 0:
            right = float(right)
        if not isinstance(right, FastDoubleFunc):
            if right == int(float(right)):
                if right == 1:
                    return left
                elif right == 2:
                    return left.unop(DUP).unop(MUL)
                elif right == 3:
                    return left.unop(DUP).unop(DUP).unop(MUL).unop(MUL)
                elif right == 4:
                    return left.unop(DUP).unop(MUL).unop(DUP).unop(MUL)
                elif right < 0:
                    return (~left)**(-right)
            right = FastDoubleFunc('const', right)
        cdef FastDoubleFunc feval = binop(left, right, POW)
        return feval

    def __neg__(FastDoubleFunc self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = -fast_float_arg(0)
            sage: f(3.5)
            -3.5
        """
        return self.unop(NEG)

    def __abs__(FastDoubleFunc self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = abs(fast_float_arg(0))
            sage: f(-3)
            3.0
        """
        return self.unop(ABS)

    def __float__(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_constant, fast_float_arg
            sage: ff = fast_float_constant(17)
            sage: float(ff)
            17.0
            sage: ff = fast_float_constant(17) - fast_float_constant(2)^2
            sage: float(ff)
            13.0
            sage: ff = fast_float_constant(17) - fast_float_constant(2)^2 + fast_float_arg(1)
            sage: float(ff)
            Traceback (most recent call last):
            ...
            TypeError: Not a constant.
        """
        if self.nargs == 0:
            return self._call_c(NULL)
        else:
            raise TypeError("Not a constant.")

    def abs(FastDoubleFunc self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).abs()
            sage: f(3)
            3.0
        """
        return self.unop(ABS)

    def __invert__(FastDoubleFunc self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = ~fast_float_arg(0)
            sage: f(4)
            0.25
        """
        return self.unop(INVERT)

    def sqrt(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).sqrt()
            sage: f(4)
            2.0
        """
        return self.cfunc(&sqrt)

    ###################################################################
    #   Exponential and log
    ###################################################################

    def log(self, base=None):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).log()
            sage: f(2)
            0.693147180559945...
            sage: f = fast_float_arg(0).log(2)
            sage: f(2)
            1.0
            sage: f = fast_float_arg(0).log(3)
            sage: f(9)
            2.0...
        """
        if base is None:
            return self.cfunc(&log)
        elif base == 2:
            return self.cfunc(&log2)
        elif base == 10:
            return self.cfunc(&log10)
        else:
            try:
                base = fast_float_constant(log(float(base)))
            except TypeError as e:
                base = fast_float(base.log())
            return binop(self.cfunc(&log), base, DIV)

    def exp(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).exp()
            sage: f(1)
            2.718281828459045...
            sage: f(100)
            2.6881171418161356e+43
        """
        return self.cfunc(&exp)

    ###################################################################
    #   Rounding
    ###################################################################

    def ceil(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).ceil()
            sage: f(1.5)
            2.0
            sage: f(-1.5)
            -1.0
        """
        return self.cfunc(&ceil)

    def floor(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).floor()
            sage: f(11.5)
            11.0
            sage: f(-11.5)
            -12.0
        """
        return self.cfunc(&floor)

    ###################################################################
    #   Trigonometric
    ###################################################################

    def sin(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).sin()
            sage: f(pi/2)
            1.0
        """
        return self.cfunc(&sin)

    def cos(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).cos()
            sage: f(0)
            1.0
        """
        return self.cfunc(&cos)

    def tan(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).tan()
            sage: f(pi/3)
            1.73205080756887...
        """
        return self.cfunc(&tan)

    def csc(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).csc()
            sage: f(pi/2)
            1.0
        """
        return ~self.sin()

    def sec(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).sec()
            sage: f(pi)
            -1.0
        """
        return ~self.cos()

    def cot(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).cot()
            sage: f(pi/4)
            1.0...
        """
        return ~self.tan()

    def arcsin(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).arcsin()
            sage: f(0.5)
            0.523598775598298...
        """
        return self.cfunc(&asin)

    def arccos(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).arccos()
            sage: f(sqrt(3)/2)
            0.5235987755982989...
        """
        return self.cfunc(&acos)

    def arctan(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).arctan()
            sage: f(1)
            0.785398163397448...
        """
        return self.cfunc(&atan)

    ###################################################################
    #   Hyperbolic
    ###################################################################

    def sinh(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).sinh()
            sage: f(log(2))
            0.75
        """
        return self.cfunc(&sinh)

    def cosh(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).cosh()
            sage: f(log(2))
            1.25
        """
        return self.cfunc(&cosh)

    def tanh(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).tanh()
            sage: f(0)
            0.0
        """
        return self.cfunc(&tanh)

    def arcsinh(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).arcsinh()
            sage: f(sinh(5))
            5.0
        """
        return self.cfunc(&asinh)

    def arccosh(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).arccosh()
            sage: f(cosh(5))
            5.0
        """
        return self.cfunc(&acosh)

    def arctanh(self):
        """
        EXAMPLES::

            sage: from sage.ext.fast_eval import fast_float_arg
            sage: f = fast_float_arg(0).arctanh()
            sage: abs(f(tanh(0.5)) - 0.5) < 0.0000001
            True
        """
        return self.cfunc(&atanh)

    cdef FastDoubleFunc cfunc(FastDoubleFunc self, void* func):
        cdef FastDoubleFunc feval = self.unop(ONE_ARG_FUNC)
        feval.ops[feval.nops - 1].params.func = <PyObject*>func
        feval.allocate_stack()
        return feval

    ###################################################################
    #   Utility functions
    ###################################################################

    cdef FastDoubleFunc unop(FastDoubleFunc self, char type):
        cdef FastDoubleFunc feval = FastDoubleFunc.__new__(FastDoubleFunc)
        feval.nargs = self.nargs
        feval.nops = self.nops + 1
        feval.max_height = self.max_height
        if type == DUP:
            feval.max_height += 1
        feval.ops = <fast_double_op *>sig_malloc(sizeof(fast_double_op) * feval.nops)
        memcpy(feval.ops, self.ops, sizeof(fast_double_op) * self.nops)
        feval.ops[feval.nops - 1].type = type
        feval.py_funcs = self.py_funcs
        feval.allocate_stack()
        return feval

cdef FastDoubleFunc binop(_left, _right, char type):
    r"""
    Returns a function that calculates left and right on the stack, leaving
    their results on the top, and then calls operation ``type``.

    EXAMPLES::

        sage: from sage.ext.fast_eval import fast_float_arg
        sage: f = fast_float_arg(1)
        sage: g = fast_float_arg(2) * 11
        sage: f.op_list()
        ['load 1']
        sage: g.op_list()
        ['load 2', 'push 11.0', 'mul']
        sage: (f+g).op_list()
        ['load 1', 'load 2', 'push 11.0', 'mul', 'add']

    Correctly calculates the maximum stack heights and number of arguments::

        sage: f.max_height
        1
        sage: g.max_height
        2
        sage: (f+g).max_height
        3
        sage: (g+f).max_height
        2

        sage: f.nargs
        2
        sage: g.nargs
        3
        sage: (f+g).nargs
        3
    """
    cdef FastDoubleFunc left, right
    try:
        left = _left
    except TypeError:
        left = fast_float(_left)
    try:
        right = _right
    except TypeError:
        right = fast_float(_right)

    # In Cython assigning None does NOT raise a TypeError above.
    if left is None or right is None:
        raise TypeError

    cdef FastDoubleFunc feval = FastDoubleFunc.__new__(FastDoubleFunc)
    feval.nargs = max(left.nargs, right.nargs)
    feval.nops = left.nops + right.nops + 1
    feval.max_height = max(left.max_height, right.max_height+1)
    feval.ops = <fast_double_op *>sig_malloc(sizeof(fast_double_op) * feval.nops)
    memcpy(feval.ops, left.ops, sizeof(fast_double_op) * left.nops)
    memcpy(feval.ops + left.nops, right.ops, sizeof(fast_double_op) * right.nops)
    feval.ops[feval.nops - 1].type = type
    if left.py_funcs is None:
        feval.py_funcs = right.py_funcs
    elif right.py_funcs is None:
        feval.py_funcs = left.py_funcs
    else:
        feval.py_funcs = left.py_funcs + right.py_funcs
    feval.allocate_stack()
    return feval


def fast_float_constant(x):
    """
    Return a fast-to-evaluate constant function.

    EXAMPLES::

        sage: from sage.ext.fast_eval import fast_float_constant
        sage: f = fast_float_constant(-2.75)
        sage: f()
        -2.75

    This is all that goes on under the hood::

        sage: fast_float_constant(pi).op_list()
        ['push 3.1415926535...']
    """
    return FastDoubleFunc('const', x)

def fast_float_arg(n):
    """
    Return a fast-to-evaluate argument selector.

    INPUT:

    - ``n`` -- the (zero-indexed) argument to select

    EXAMPLES::

        sage: from sage.ext.fast_eval import fast_float_arg
        sage: f = fast_float_arg(0)
        sage: f(1,2)
        1.0
        sage: f = fast_float_arg(1)
        sage: f(1,2)
        2.0

    This is all that goes on under the hood::

        sage: fast_float_arg(10).op_list()
        ['load 10']
    """
    return FastDoubleFunc('arg', n)

def fast_float_func(f, *args):
    """
    Returns a wrapper around a python function.

    INPUT:

    - ``f`` -- a callable python object
    - ``args`` -- a list of FastDoubleFunc inputs

    EXAMPLES::

        sage: from sage.ext.fast_eval import fast_float_func, fast_float_arg
        sage: f = fast_float_arg(0)
        sage: g = fast_float_arg(1)
        sage: h = fast_float_func(lambda x,y: x-y, f, g)
        sage: h(5, 10)
        -5.0

    This is all that goes on under the hood::

        sage: h.op_list()
        ['load 0', 'load 1', 'py_call <function <lambda> at 0x...>(2)']
    """
    return FastDoubleFunc('callable', f, *args)


new_fast_float=True

def fast_float(f, *vars, old=None, expect_one_var=False):
    """
    Tries to create a function that evaluates f quickly using
    floating-point numbers, if possible.  There are two implementations
    of fast_float in Sage; by default we use the newer, which is
    slightly faster on most tests.

    On failure, returns the input unchanged.

    INPUT:

    - ``f``    -- an expression
    - ``vars`` -- the names of the arguments
    - ``old``  -- use the original algorithm for fast_float
    - ``expect_one_var`` -- don't give deprecation warning if ``vars`` is
      omitted, as long as expression has only one var

    EXAMPLES::

        sage: from sage.ext.fast_eval import fast_float
        sage: x,y = var('x,y')
        sage: f = fast_float(sqrt(x^2+y^2), 'x', 'y')
        sage: f(3,4)
        5.0

    Specifying the argument names is essential, as fast_float objects
    only distinguish between arguments by order. ::

        sage: f = fast_float(x-y, 'x','y')
        sage: f(1,2)
        -1.0
        sage: f = fast_float(x-y, 'y','x')
        sage: f(1,2)
        1.0
    """
    if old is None:
        old = not new_fast_float

    if isinstance(f, (tuple, list)):
        return tuple([fast_float(x, *vars, expect_one_var=expect_one_var) for x in f])

    cdef int i
    for i from 0 <= i < len(vars):
        if not isinstance(vars[i], str):
            v = str(vars[i])
            # inexact generators display as 1.00..0*x
            if '*' in v:
                v = v[v.index('*')+1:]
            vars = vars[:i] + (v,) + vars[i+1:]

    try:
        if old:
            return f._fast_float_(*vars)
        else:
            return fast_callable(f, vars=vars, domain=float,
                                 expect_one_var=expect_one_var)
    except AttributeError:
        pass

    try:
        return FastDoubleFunc('const', float(f))
    except TypeError:
        pass

    try:
        from sage.symbolic.ring import SR
        return fast_float(SR(f), *vars)
    except TypeError:
        pass

    if f is None:
        raise TypeError("no way to make fast_float from None")

    return f


def is_fast_float(x):
    return isinstance(x, FastDoubleFunc) or isinstance(x, Wrapper)
