"Operands"

# ****************************************************************************
#       Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object cimport SageObject


cdef inline int normalize_index(object arg, int nops, object err_msg) except -1:
    """
    Given an index ``arg`` and the number of operands ``nops`` return
    an integer between 0 and ``nops`` which will be used as the index to fetch
    the data from the underlying vector.

    TESTS::

        sage: from sage.symbolic.expression import normalize_index_for_doctests as normalize_index
        sage: normalize_index(-1, 4)
        3
        sage: normalize_index(1, 4)
        1
        sage: normalize_index(1.5, 4)
        Traceback (most recent call last):
        ...
        TypeError: some error message
        sage: normalize_index(-5, 4)
        Traceback (most recent call last):
        ...
        IndexError: operand index out of range, got -5, expect between -4 and 3
        sage: normalize_index(5, 4)
        Traceback (most recent call last):
        ...
        IndexError: operand index out of range, got 5, expect between -4 and 3
    """
    cdef int i
    i = arg
    if i != arg:
        raise TypeError(err_msg)
    if i < 0:
        i = nops + i
    if i >= nops or i < 0:
        raise IndexError("operand index out of range, got %s, expect between %s and %s" % (arg, -nops, nops-1))
    return i


def normalize_index_for_doctests(arg, nops):
    """
    Wrapper function to test ``normalize_index``.

    TESTS::

        sage: from sage.symbolic.expression import normalize_index_for_doctests
        sage: normalize_index_for_doctests(-1, 4)
        3
    """
    return normalize_index(arg, nops, "some error message")

cdef class OperandsWrapper(SageObject):
    """
    Operands wrapper for symbolic expressions.

    EXAMPLES::

        sage: x,y,z = var('x,y,z')
        sage: e = x + x*y + z^y + 3*y*z; e
        x*y + 3*y*z + x + z^y
        sage: e.op[1]
        3*y*z
        sage: e.op[1,1]
        z
        sage: e.op[-1]
        z^y
        sage: e.op[1:]
        [3*y*z, x, z^y]
        sage: e.op[:2]
        [x*y, 3*y*z]
        sage: e.op[-2:]
        [x, z^y]
        sage: e.op[:-2]
        [x*y, 3*y*z]
        sage: e.op[-5]
        Traceback (most recent call last):
        ...
        IndexError: operand index out of range, got -5, expect between -4 and 3
        sage: e.op[5]
        Traceback (most recent call last):
        ...
        IndexError: operand index out of range, got 5, expect between -4 and 3
        sage: e.op[1,1,0]
        Traceback (most recent call last):
        ...
        TypeError: expressions containing only a numeric coefficient, constant or symbol have no operands
        sage: e.op[:1.5]
        Traceback (most recent call last):
        ...
        TypeError: slice indices must be integers or None or have an __index__ method
        sage: e.op[:2:1.5]
        Traceback (most recent call last):
        ...
        ValueError: step value must be an integer
    """

    cdef Expression _expr

    def __getitem__(self, arg):
        """
        TESTS::

           sage: t = 1+x+x^2
           sage: t.op[1:]
           [x, 1]
        """
        cdef int ind, nops
        cdef int bind, eind, step
        cdef GEx cur_ex
        if isinstance(arg, slice):
            nops = self._expr._gobj.nops()

            slice_err_msg = "slice indices must be integers or None or have an __index__ method"
            if arg.start:
                bind = normalize_index(arg.start, nops, slice_err_msg)
            else:
                bind = 0
            if arg.stop:
                eind = normalize_index(arg.stop, nops, slice_err_msg)
            else:
                eind = nops
            if arg.step:
                step = arg.step
                if step != arg.step:
                    raise ValueError("step value must be an integer")
            else:
                step = 1
            return [new_Expression_from_GEx(self._expr._parent,
                self._expr._gobj.op(ind)) for ind in range(bind, eind, step)]

        ind_err_msg = "index should either be a slice object, an integer or a list of integers"
        if isinstance(arg, (list, tuple)):
            # handle nested index
            if len(arg) == 0: # or not all(lambda x: x in ZZ for t in args):
                raise TypeError(ind_err_msg)
            cur_ex = GEx(self._expr._gobj)
            for x in arg:
                nops = cur_ex.nops()
                if nops == 0:
                    raise TypeError("expressions containing only a numeric coefficient, constant or symbol have no operands")
                ind = normalize_index(x, nops, ind_err_msg)
                cur_ex = cur_ex.op(ind)
            return new_Expression_from_GEx(self._expr._parent, cur_ex)

        ind = normalize_index(arg, self._expr._gobj.nops(), ind_err_msg)
        return new_Expression_from_GEx(self._expr._parent,
                self._expr._gobj.op(ind))

    def _repr_(self):
        """
        TESTS::

            sage: (x^2).op
            Operands of x^2
        """
        return "Operands of %s" % (self._expr)

    def _latex_(self):
        r"""
        TESTS::

            sage: latex((x^2).op)
            \text{Operands wrapper for expression }x^{2}
        """
        return r"\text{Operands wrapper for expression }%s" % (self._expr._latex_())

    def __reduce__(self):
        """
        TESTS::

            sage: (x^2).op.__reduce__()
            (<built-in function restore_op_wrapper>, (x^2,))
            sage: loads(dumps((x^2).op))
            Operands of x^2
        """
        return restore_op_wrapper, (self._expr,)


def restore_op_wrapper(expr):
    """
    TESTS::

        sage: from sage.symbolic.expression import restore_op_wrapper
        sage: restore_op_wrapper(x^2)
        Operands of x^2
    """
    return expr.op
