# sage_setup: distribution = sagemath-objects
"""
Mutability Cython Implementation
"""
##########################################################################
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
##########################################################################

from sage.misc.decorators import sage_wraps

cdef class Mutability:
    r"""
    Class to mix in mutability feature.

    EXAMPLES::

        sage: class A(SageObject, Mutability):
        ....:     def __init__(self, val):
        ....:         self._val = val
        ....:     def change(self, val):
        ....:         self._require_mutable()
        ....:         self._val = val
        ....:     def __hash__(self):
        ....:         self._require_immutable()
        ....:         return hash(self._val)
        sage: a = A(4)
        sage: a._val
        4
        sage: a.change(6); a._val
        6
        sage: hash(a)
        Traceback (most recent call last):
        ...
        ValueError: object is mutable; please make it immutable first
        sage: a.set_immutable()
        sage: a.change(4)
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead
        sage: hash(a)
        6
    """

    def __init__(self, is_immutable=False):
        r"""
        TESTS::

            sage: class A(SageObject, Mutability):
            ....:     def __init__(self, val):
            ....:         self._val = val
            ....:     def change(self, val):
            ....:         self._require_mutable()
            ....:         self._val = val
            ....:     def __hash__(self):
            ....:         self._require_immutable()
            ....:         return hash(self._val)
            sage: a = A(4)
            sage: TestSuite(a).run(skip ='_test_pickling')
        """
        self._is_immutable = is_immutable

    cpdef _require_mutable(self):
        r"""
        Whenever mutability is required, this method can be called.

        EXAMPLES::

            sage: class A(SageObject, Mutability):
            ....:     def __init__(self, val):
            ....:         self._val = val
            ....:     def change(self, val):
            ....:         self._require_mutable()
            ....:         self._val = val
            ....:     def __hash__(self):
            ....:         self._require_immutable()
            ....:         return hash(self._val)
            sage: a = A(4)
            sage: a.set_immutable()
            sage: a.change(6)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead
        """
        if self._is_immutable:
            raise ValueError("object is immutable; please change a copy instead")

    cpdef _require_immutable(self):
        r"""
        Whenever immutability is required, this method can be called.

        EXAMPLES::

            sage: class A(SageObject, Mutability):
            ....:     def __init__(self, val):
            ....:         self._val = val
            ....:     def change(self, val):
            ....:         self._require_mutable()
            ....:         self._val = val
            ....:     def __hash__(self):
            ....:         self._require_immutable()
            ....:         return hash(self._val)
            sage: a = A(4)
            sage: hash(a)
            Traceback (most recent call last):
            ...
            ValueError: object is mutable; please make it immutable first
        """
        if not self._is_immutable:
            raise ValueError("object is mutable; please make it immutable first")

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

    cpdef bint is_immutable(self) noexcept:
        """
        Return ``True`` if this object is immutable (cannot be changed)
        and ``False`` if it is not.

        To make this object immutable use self.set_immutable().

        EXAMPLES::

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
        return self._is_immutable

    cpdef bint is_mutable(self) noexcept:
        """
        Return ``True`` if this object is mutable (can be changed)
        and ``False`` if it is not.

        To make this object immutable use ``self.set_immutable()``.

        EXAMPLES::

            sage: v = Sequence([1,2,3,4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.is_mutable()
            True
            sage: v.set_immutable()
            sage: v.is_mutable()
            False
        """
        return not self._is_immutable

    def __getstate__(self):
        r"""
        Get the current state of ``self`` including the mutability status.

        TESTS::

            sage: class A(SageObject, Mutability):
            ....:     def __init__(self, val):
            ....:         self._val = val
            ....:     def change(self, val):
            ....:         self._require_mutable()
            ....:         self._val = val
            ....:     def __hash__(self):
            ....:         self._require_immutable()
            ....:         return hash(self._val)
            sage: a = A(4)
            sage: a.__dict__
            {'_val': 4}
            sage: a.__getstate__()
            {'_is_immutable': False, '_val': 4}
            sage: a.__reduce__()  # indirect doctest
            (<function _reconstructor at ...>,
             (<class '__main__.A'>,
              <class 'sage.structure.sage_object.SageObject'>,
              <sage.structure.sage_object.SageObject object at ...>),
             {'_is_immutable': False, '_val': 4})
        """
        state = getattr(self, '__dict__', {})
        state['_is_immutable'] = self._is_immutable
        return state

    def __setstate__(self, state):
        r"""
        Set the state of ``self`` from the dictionary ``state`` including the
        mutability status.

        TESTS::

            sage: class A(SageObject, Mutability):
            ....:     def __init__(self, val):
            ....:         self._val = val
            ....:     def change(self, val):
            ....:         self._require_mutable()
            ....:         self._val = val
            ....:     def __hash__(self):
            ....:         self._require_immutable()
            ....:         return hash(self._val)
            sage: a = A(4)
            sage: a.is_immutable()
            False
            sage: d = a.__getstate__(); d
            {'_is_immutable': False, '_val': 4}
            sage: d['_is_immutable'] = True
            sage: a.__setstate__(d)
            sage: a.is_immutable()
            True
            sage: a.__getstate__()
            {'_is_immutable': True, '_val': 4}
        """
        if hasattr(self, '__dict__'):
            self.__dict__ = state
        self._is_immutable = state['_is_immutable']


##########################################################################
# Method decorators for mutating methods resp. methods that assume immutability

def require_mutable(f):
    """
    A decorator that requires mutability for a method to be called.

    .. NOTE::

        Objects whose methods use this decorator should have an attribute
        ``_is_immutable``. Otherwise, the object is assumed to be mutable.

    EXAMPLES::

        sage: from sage.structure.mutability import require_mutable, require_immutable
        sage: class A():
        ....:     def __init__(self, val):
        ....:         self._m = val
        ....:     @require_mutable
        ....:     def change(self, new_val):
        ....:         'change self'
        ....:         self._m = new_val
        ....:     @require_immutable
        ....:     def __hash__(self):
        ....:         'implement hash'
        ....:         return hash(self._m)
        sage: a = A(5)
        sage: a.change(6)
        sage: hash(a)
        Traceback (most recent call last):
        ...
        ValueError: <class '__main__.A'> instance is mutable, <function ...__hash__ at ...> must not be called
        sage: a._is_immutable = True
        sage: hash(a)
        6
        sage: a.change(7)   # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: <class '__main__.A'> instance is immutable, <function ...change at ...> must not be called
        sage: from sage.misc.sageinspect import sage_getdoc
        sage: print(sage_getdoc(a.change))
        change self

    AUTHORS:

    - Simon King <simon.king@uni-jena.de>
    """
    @sage_wraps(f)
    def new_f(self, *args, **kwds):
        if getattr(self, '_is_immutable', False):
            raise ValueError("{} instance is immutable, {} must not be called".format(type(self), repr(f)))
        return f(self, *args, **kwds)
    return new_f


def require_immutable(f):
    """
    A decorator that requires immutability for a method to be called.

    .. NOTE::

        Objects whose methods use this decorator should have an attribute
        ``_is_immutable``. Otherwise, the object is assumed to be mutable.

    EXAMPLES::

        sage: from sage.structure.mutability import require_mutable, require_immutable
        sage: class A():
        ....:  def __init__(self, val):
        ....:      self._m = val
        ....:  @require_mutable
        ....:  def change(self, new_val):
        ....:      'change self'
        ....:      self._m = new_val
        ....:  @require_immutable
        ....:  def __hash__(self):
        ....:      'implement hash'
        ....:      return hash(self._m)
        sage: a = A(5)
        sage: a.change(6)
        sage: hash(a)   # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: <class '__main__.A'> instance is mutable, <function ...__hash__ at ...> must not be called
        sage: a._is_immutable = True
        sage: hash(a)
        6
        sage: a.change(7)
        Traceback (most recent call last):
        ...
        ValueError: <class '__main__.A'> instance is immutable, <function ...change at ...> must not be called
        sage: from sage.misc.sageinspect import sage_getdoc
        sage: print(sage_getdoc(a.__hash__))
        implement hash

    AUTHORS:

    - Simon King <simon.king@uni-jena.de>
    """
    @sage_wraps(f)
    def new_f(self, *args, **kwds):
        if not getattr(self,'_is_immutable',False):
            raise ValueError("{} instance is mutable, {} must not be called".format(type(self), repr(f)))
        return f(self, *args, **kwds)
    return new_f
