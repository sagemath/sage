r"""

EXAMPLES::

    sage: R = PolynomialRing(ZZ, 3, 'x')
    sage: R.ngens()
    3
    sage: R.gen(0)
    x0
    sage: R.gens()
    (x0, x1, x2)
    sage: R.variable_names()
    ('x0', 'x1', 'x2')

    sage: M = FreeModule(ZZ, 4)
    sage: M
    Ambient free module of rank 4 over the principal ideal domain Integer Ring
    sage: M.ngens()
    4
    sage: M.gen(0)
    (1, 0, 0, 0)
    sage: M.gens()
    ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))
"""

# ****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport sage.structure.category_object as category_object

cdef class localvars:
    r"""
    Context manager for safely temporarily changing the variables
    names of an object with generators.

    Objects with named generators are globally unique in Sage.
    Sometimes, though, it is very useful to be able to temporarily
    display the generators differently.   The new Python ``with``
    statement and the localvars context manager make this easy and
    safe (and fun!)

    Suppose X is any object with generators.  Write

    ::

        with localvars(X, names[, latex_names] [,normalize=False]):
             some code
             ...

    and the indented code will be run as if the names in X are changed
    to the new names.  If you give normalize=True, then the names are
    assumed to be a tuple of the correct number of strings.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: with localvars(R, 'z,w'):
        ....:     print(x^3 + y^3 - x*y)
        z^3 + w^3 - z*w

    .. NOTE::

       I wrote this because it was needed to print elements of the
       quotient of a ring R by an ideal I using the print function for
       elements of R.  See the code in
       ``quotient_ring_element.pyx``.

    AUTHOR:

    - William Stein (2006-10-31)
    """
    cdef object _obj
    cdef object _names
    cdef object _latex_names
    cdef object _orig

    def __init__(self, obj, names, latex_names=None, normalize=True):
        self._obj = obj
        if normalize:
            self._names = category_object.normalize_names(obj.ngens(), names)
            self._latex_names = latex_names
        else:
            self._names = names
            self._latex_names = latex_names

    def __enter__(self):
        self._orig = self._obj._temporarily_change_names(self._names, self._latex_names)

    def __exit__(self, type, value, traceback):
        self._obj._temporarily_change_names(self._orig[0], self._orig[1])
