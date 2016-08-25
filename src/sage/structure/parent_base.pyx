r"""
Base class for old-style parent objects with a base ring
"""
###############################################################################
#   Sage: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "sage/ext/stdsage.pxi"

cimport parent

from coerce_exceptions import CoercionException

cdef inline check_old_coerce(parent.Parent p):
    if p._element_constructor is not None:
        raise RuntimeError("%s still using old coercion framework" % p)


# TODO: Unpickled parents with base sometimes have their base set to None.
# This causes a segfault in the module arithmetic architecture.
#
# sage: H = HomsetWithBase(QQ, RR, base=ZZ); H
# sage: H0 = loads(dumps(H))
# sage: H.base_ring(), H0.base_ring()
# (Integer Ring, None)
#
# Perhaps the code below would help (why was it commented out?).

## def make_parent_with_base_v0(_class, _dict, base, has_coerce_map_from):
##     """
##     This should work for any Python class deriving from this, as long
##     as it doesn't implement some screwy __new__() method.
##     """
##     new_object = _class.__new__(_class)
##     if base is None:
##         (<ParentWithBase>new_object)._base = new_object
##     else:
##         (<ParentWithBase>new_object)._base = base
##     (<ParentWithBase>new_object)._has_coerce_map_from = has_coerce_map_from
##     if not _dict is None:
##         new_object.__dict__ = _dict
##     return new_object

def is_ParentWithBase(x):
    """
    Return True if x is a parent object with base.
    """
    return isinstance(x, ParentWithBase)

cdef class ParentWithBase(Parent_old):
    """
    This class is being deprecated, see parent.Parent for the new model.
    """
    def __init__(self, base, coerce_from=[], actions=[], embeddings=[], category=None):
        # TODO: SymbolicExpressionRing has base RR, which makes this bad
#        print type(self), "base", base, coerce_from
#        if base != self and not base in coerce_from:
#            coerce_from.append(base)
        Parent_old.__init__(self, coerce_from=coerce_from, actions=actions, embeddings=embeddings, category=category)
        self._base = base

    cdef _coerce_c_impl(self,x):
       check_old_coerce(self)
       if not self._base is self:
           return self._coerce_try(x,(self._base))
       else:
           raise TypeError("No canonical coercion found.")

##     def x__reduce__(self):
##         if HAS_DICTIONARY(self):
##             _dict = self.__dict__
##         else:
##             _dict = None
##         if self._base is self:
##             return (make_parent_with_base_v0, (self.__class__, _dict, None, self._has_coerce_map_from))
##         else:
##             return (make_parent_with_base_v0, (self.__class__, _dict, self._base, self._has_coerce_map_from))

    # Derived class *must* define base_extend.
    def base_extend(self, X):
        check_old_coerce(self)
        raise CoercionException("BUG: the base_extend method must be defined for '%s' (class '%s')" %
                                (self, type(self)))
