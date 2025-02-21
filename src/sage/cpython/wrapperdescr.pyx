# sage_setup: distribution = sagemath-objects
"""
Slot wrappers

A slot wrapper is installed in the dict of an extension type to
access a special method implemented in C. For example,
``object.__init__`` or ``Integer.__lt__``. Note that slot wrappers
are always unbound (there is a bound variant called method-wrapper).

EXAMPLES::

    sage: int.__add__
    <slot wrapper '__add__' of 'int' objects>

Pure Python classes have normal methods, not slot wrappers::

    sage: class X():
    ....:     def __add__(self, other):
    ....:         return NotImplemented
    sage: X.__add__
    <function X.__add__ at ...>
"""

# ****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.cpython.string import bytes_to_str


def wrapperdescr_call(slotwrapper, self, *args, **kwds):
    """
    Call a slot wrapper without any type checks.

    The main reason to use this is to call arithmetic slots like
    ``__mul__`` without having to worry about whether to call
    ``T.__mul__(a, b)`` or ``T.__rmul__(b, a)``.

    INPUT:

    - ``slotwrapper`` -- a slot wrapper (for example ``int.__add__``)

    - ``self`` -- the first positional argument. Normally, this should
      be of the correct type (an ``int`` when calling ``int.__add__``).
      However, this check is skipped: you can pass an arbitrary object.

    - ``*args``, ``**kwds`` -- further arguments

    .. WARNING::

        Since this skips type checks, it can easily crash Python if
        used incorrectly.

    EXAMPLES::

        sage: from sage.cpython.wrapperdescr import wrapperdescr_call
        sage: wrapperdescr_call(Integer.__mul__, 6, 9)
        54
        sage: wrapperdescr_call(Integer.__mul__, 7/5, 9)
        63/5
        sage: from sage.structure.element import Element
        sage: wrapperdescr_call(Element.__mul__, 6, 9)
        54
        sage: wrapperdescr_call(Element.__mul__, 7/5, 9)
        63/5
        sage: from sage.numerical.mip import MixedIntegerLinearProgram                  # needs sage.numerical.mip
        sage: wrapperdescr_call(type.__call__,                                          # needs sage.numerical.mip
        ....:                   MixedIntegerLinearProgram, maximization=False)
        Mixed Integer Program (no objective, 0 variables, 0 constraints)

    TESTS::

        sage: wrapperdescr_call(Integer.__mul__, 1, 2, 3)
        Traceback (most recent call last):
        ...
        TypeError: expected 1 arg..., got 2
        sage: wrapperdescr_call(Integer.__mul__, 6, other=9)
        Traceback (most recent call last):
        ...
        TypeError: wrapper __mul__ slotdef doesn't take keyword arguments
    """
    return wrapperdescr_fastcall(slotwrapper, self, args, kwds)


cdef wrapperdescr_fastcall(wrapper_descriptor slotwrapper, self, args, kwds):
    # Cython implementation of wrapperdescr_call
    cdef wrapperbase* slotdef = slotwrapper.d_base

    cdef wrapperfunc_kwds wk
    if slotdef.flags & PyWrapperFlag_KEYWORDS:
        wk = <wrapperfunc_kwds>(slotdef.wrapper)
        return wk(self, args, slotwrapper.d_wrapped, kwds)

    if <PyObject*>kwds is not NULL and kwds:
        raise TypeError(f"wrapper {bytes_to_str(slotdef.name)} slotdef "
                        "doesn't take keyword arguments")

    return slotdef.wrapper(self, args, slotwrapper.d_wrapped)
