r"""
Base class for old-style parent objects with a base ring
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport sage.structure.parent as parent

from .coerce_exceptions import CoercionException

cdef inline check_old_coerce(parent.Parent p):
    if p._element_constructor is not None:
        raise RuntimeError("%s still using old coercion framework" % p)


def is_ParentWithBase(x):
    """
    Return True if x is a parent object with base.
    """
    return isinstance(x, ParentWithBase)

cdef class ParentWithBase(Parent_old):
    """
    This class is being deprecated, see parent.Parent for the new model.
    """
    def __init__(self, base, *args, **kwds):
        Parent_old.__init__(self, *args, **kwds)
        self._base = base

    cdef _coerce_c_impl(self,x):
       check_old_coerce(self)
       if not self._base is self:
           return self._coerce_try(x,(self._base))
       else:
           raise TypeError("No canonical coercion found.")

    # Derived class *must* define base_extend.
    def base_extend(self, X):
        check_old_coerce(self)
        raise CoercionException("BUG: the base_extend method must be defined for '%s' (class '%s')" %
                                (self, type(self)))
