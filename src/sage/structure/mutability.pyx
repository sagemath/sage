"""
Mutability Cython Implementation
"""

##########################################################################
#
#   Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
##########################################################################
from __future__ import print_function


cdef class Mutability:

    def __init__(self, is_immutable=False):
        self._is_immutable = is_immutable

    def _require_mutable(self):
        if self._is_immutable:
            raise ValueError("object is immutable; please change a copy instead.")

    cdef _require_mutable_cdef(self):
        if self._is_immutable:
            raise ValueError("object is immutable; please change a copy instead.")

    def set_immutable(self):
        """
        Make this object immutable, so it can never again be changed.

        EXAMPLES::

            sage: v = Sequence([1,2,3,4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.set_immutable()
            sage: v[3] = 7
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        self._is_immutable = 1

    def is_immutable(self):
        """
        Return True if this object is immutable (can not be changed)
        and False if it is not.

        To make this object immutable use self.set_immutable().

        EXAMPLE::

            sage: v = Sequence([1,2,3,4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.is_immutable()
            False
            sage: v.set_immutable()
            sage: v.is_immutable()
            True
        """
        self._is_immutable

    def is_mutable(self):
        return not self._is_immutable

    def __reduce__(self):
        return Mutability, (self._is_immutable, )

#################################################################################
## Method decorators for mutating methods resp. methods that assume immutability
from sage.misc.decorators import sage_wraps

def require_mutable(f):
    """
    A decorator that requires mutability for a method to be called.

    EXAMPLES::

        sage: from sage.structure.mutability import require_mutable, require_immutable
        sage: class A:
        ...    def __init__(self, val):
        ...        self._m = val
        ...    @require_mutable
        ...    def change(self, new_val):
        ...        'change self'
        ...        self._m = new_val
        ...    @require_immutable
        ...    def __hash__(self):
        ...        'implement hash'
        ...        return hash(self._m)
        sage: a = A(5)
        sage: a.change(6)
        sage: hash(a)
        Traceback (most recent call last):
        ...
        ValueError: <type 'instance'> instance is mutable, <function __hash__ at ...> must not be called
        sage: a._is_immutable = True
        sage: hash(a)
        6
        sage: a.change(7)   # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: <type 'instance'> instance is immutable, <function change at ...> must not be called
        sage: from sage.misc.sageinspect import sage_getdoc
        sage: print(sage_getdoc(a.change))
        change self

    AUTHORS:

    - Simon King <simon.king@uni-jena.de>
    """
    @sage_wraps(f)
    def new_f(self, *args,**kwds):
        if getattr(self,'_is_immutable',False):
            raise ValueError("%s instance is immutable, %s must not be called"%(type(self), repr(f)))
        return f(self, *args,**kwds)
    return new_f

def require_immutable(f):
    """
    A decorator that requires mutability for a method to be called.

    EXAMPLES::

        sage: from sage.structure.mutability import require_mutable, require_immutable
        sage: class A:
        ...    def __init__(self, val):
        ...        self._m = val
        ...    @require_mutable
        ...    def change(self, new_val):
        ...        'change self'
        ...        self._m = new_val
        ...    @require_immutable
        ...    def __hash__(self):
        ...        'implement hash'
        ...        return hash(self._m)
        sage: a = A(5)
        sage: a.change(6)
        sage: hash(a)   # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: <type 'instance'> instance is mutable, <function __hash__ at ...> must not be called
        sage: a._is_immutable = True
        sage: hash(a)
        6
        sage: a.change(7)
        Traceback (most recent call last):
        ...
        ValueError: <type 'instance'> instance is immutable, <function change at ...> must not be called
        sage: from sage.misc.sageinspect import sage_getdoc
        sage: print(sage_getdoc(a.__hash__))
        implement hash

    AUTHORS:

    - Simon King <simon.king@uni-jena.de>
    """
    @sage_wraps(f)
    def new_f(self, *args,**kwds):
        if not getattr(self,'_is_immutable',False):
            raise ValueError("%s instance is mutable, %s must not be called"%(type(self), repr(f)))
        return f(self, *args,**kwds)
    return new_f
