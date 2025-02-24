# sage_setup: distribution = sagemath-objects
r"""
Elements, Array and Lists With Clone Protocol

This module defines several classes which are subclasses of
:class:`Element<sage.structure.element.Element>` and which roughly implement
the "prototype" design pattern (see [Prototype_pattern]_, [GHJV1994]_). Those classes are
intended to be used to model *mathematical* objects, which are by essence
immutable. However, in many occasions, one wants to construct the
data-structure encoding of a new mathematical object by small modifications of
the data structure encoding some already built object. For the resulting
data-structure to correctly encode the mathematical object, some structural
invariants must hold. One problem is that, in many cases, during the
modification process, there is no possibility but to break the invariants.

For example, in a list modeling a permutation of `\{1,\dots,n\}`, the integers
must be distinct. A very common operation is to take a permutation to make a
copy with some small modifications, like exchanging two consecutive values in
the list or cycling some values. Though the result is clearly a permutations
there's no way to avoid breaking the permutations invariants at some point
during the modifications.

The main purpose of this module is to define the class

- :class:`ClonableElement` as an abstract super class,

and its subclasses:

- :class:`ClonableArray` for arrays (lists of fixed length) of objects;
- :class:`ClonableList` for (resizable) lists of objects;
- :class:`NormalizedClonableList` for lists of objects with a normalization method;
- :class:`ClonableIntArray` for arrays of int.

.. SEEALSO:: The following parents from :mod:`sage.structure.list_clone_demo`
    demonstrate how to use them:

    - ``IncreasingArrays()`` (see
      :class:`~sage.structure.list_clone_demo.IncreasingArray`
      and the parent class
      :class:`~sage.structure.list_clone_demo.IncreasingArrays`)
    - ``IncreasingLists()`` (see
      :class:`~sage.structure.list_clone_demo.IncreasingList`
      and the parent class
      :class:`~sage.structure.list_clone_demo.IncreasingLists`)
    - ``SortedLists()`` (see
      :class:`~sage.structure.list_clone_demo.SortedList`
      and the parent class
      :class:`~sage.structure.list_clone_demo.SortedLists`)
    - ``IncreasingIntArrays()`` (see
      :class:`~sage.structure.list_clone_demo.IncreasingIntArray`
      and the parent class
      :class:`~sage.structure.list_clone_demo.IncreasingIntArrays`)

EXAMPLES:

We now demonstrate how
:class:`~sage.structure.list_clone_demo.IncreasingArray` works, creating an
instance ``el`` through its parent ``IncreasingArrays()``::

    sage: from sage.structure.list_clone_demo import IncreasingArrays
    sage: P = IncreasingArrays()
    sage: P([1, 4 ,8])
    [1, 4, 8]

If one tries to create this way a list which in not increasing, an error is
raised::

    sage: IncreasingArrays()([5, 4 ,8])
    Traceback (most recent call last):
    ...
    ValueError: array is not increasing

Once created modifying ``el`` is forbidden::

    sage: el = P([1, 4 ,8])
    sage: el[1] = 3
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

However, you can modify a temporarily mutable clone::

    sage: with el.clone() as elc:
    ....:      elc[1] = 3
    sage: [el, elc]
    [[1, 4, 8], [1, 3, 8]]

We check that the original and the modified copy now are in a proper immutable
state::

    sage: el.is_immutable(), elc.is_immutable()
    (True, True)
    sage: elc[1] = 5
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

You can break the property that the list is increasing during the
modification::

    sage: with el.clone() as elc2:
    ....:     elc2[1] = 12
    ....:     print(elc2)
    ....:     elc2[2] = 25
    [1, 12, 8]
    sage: elc2
    [1, 12, 25]

But this property must be restored at the end of the ``with`` block; otherwise
an error is raised::

    sage: with elc2.clone() as el3:
    ....:     el3[1] = 100
    Traceback (most recent call last):
    ...
    ValueError: array is not increasing

Finally, as an alternative to the ``with`` syntax one can use::

    sage: el4 = copy(elc2)
    sage: el4[1] = 10
    sage: el4.set_immutable()
    sage: el4.check()

REFERENCES:

- [Prototype_pattern]_
- [GHJV1994]_

AUTHORS:

- Florent Hivert (2010-03): initial revision
"""
# ****************************************************************************
#       Copyright (C) 2009-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.list cimport *
from cpython.long cimport *
from cpython.ref cimport *

from cysignals.memory cimport check_reallocarray, sig_free

from sage.ext.stdsage cimport HAS_DICTIONARY
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport richcmp, rich_to_bool

############################################################################
#                          Basic clone elements                            #
############################################################################
cdef class ClonableElement(Element):
    r"""
    Abstract class for elements with clone protocol.

    This class is a subclass of :class:`Element<sage.structure.element.Element>`
    and implements the "prototype" design pattern (see [Prototype_pattern]_, [GHJV1994]_). The role
    of this class is:

    - to manage copy and mutability and hashing of elements
    - to ensure that at the end of a piece of code an object is restored in a
      meaningful mathematical state.

    A class ``C`` inheriting from :class:`ClonableElement` must implement
    the following two methods

    - ``obj.__copy__()`` -- returns a fresh copy of obj
    - ``obj.check()`` -- returns nothing, raise an exception if ``obj``
      doesn't satisfy the data structure invariants

    and ensure to call ``obj._require_mutable()`` at the beginning of any
    modifying method.

    Additionally, one can also implement

    - ``obj._hash_()`` -- return the hash value of ``obj``

    Then, given an instance ``obj`` of ``C``, the following sequences of
    instructions ensures that the invariants of ``new_obj`` are properly
    restored at the end::

       with obj.clone() as new_obj:
           ...
           # lot of invariant breaking modifications on new_obj
           ...
       # invariants are ensured to be fulfilled

    The following equivalent sequence of instructions can be used if speed is
    needed, in particular in Cython code::

       new_obj = obj.__copy__()
       ...
       # lot of invariant breaking modifications on new_obj
       ...
       new_obj.set_immutable()
       new_obj.check()
       # invariants are ensured to be fulfilled

    Finally, if the class implements the ``_hash_`` method, then
    :class:`ClonableElement` ensures that the hash value can only be
    computed on an immutable object. It furthermore caches it so that
    it is only computed once.

    .. warning:: for the hash caching mechanism to work correctly, the hash
       value cannot be 0.

    EXAMPLES:

    The following code shows a minimal example of usage of
    :class:`ClonableElement`. We implement a class or pairs `(x,y)`
    such that `x < y`::

        sage: from sage.structure.list_clone import ClonableElement
        sage: class IntPair(ClonableElement):
        ....:      def __init__(self, parent, x, y):
        ....:          ClonableElement.__init__(self, parent=parent)
        ....:          self._x = x
        ....:          self._y = y
        ....:          self.set_immutable()
        ....:          self.check()
        ....:      def _repr_(self):
        ....:          return "(x=%s, y=%s)"%(self._x, self._y)
        ....:      def check(self):
        ....:          if self._x >= self._y:
        ....:              raise ValueError("Incorrectly ordered pair")
        ....:      def get_x(self): return self._x
        ....:      def get_y(self): return self._y
        ....:      def set_x(self, v): self._require_mutable(); self._x = v
        ....:      def set_y(self, v): self._require_mutable(); self._y = v

    .. NOTE:: we don't need to define ``__copy__`` since it is properly
       inherited from :class:`Element<sage.structure.element.Element>`.

    We now demonstrate the behavior. Let's create an ``IntPair``::

        sage: myParent = Parent()
        sage: el = IntPair(myParent, 1, 3); el
        (x=1, y=3)
        sage: el.get_x()
        1

    Modifying it is forbidden::

        sage: el.set_x(4)
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead.

    However, you can modify a mutable copy::

        sage: with el.clone() as el1:
        ....:      el1.set_x(2)
        sage: [el, el1]
        [(x=1, y=3), (x=2, y=3)]

    We check that the original and the modified copy are in a proper immutable
    state::

        sage: el.is_immutable(), el1.is_immutable()
        (True, True)
        sage: el1.set_x(4)
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead.

    A modification which doesn't restore the invariant `x < y` at the end is
    illegal and raise an exception::

        sage: with el.clone() as elc2:
        ....:      elc2.set_x(5)
        Traceback (most recent call last):
        ...
        ValueError: Incorrectly ordered pair
    """
    def __cinit__(self):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3]) # indirect doctest
            sage: el.is_immutable()
            True
        """
        self._needs_check = True
        self._is_immutable = False
        self._hash = 0

    cpdef bint _require_mutable(self) except -2:
        """
        Check that ``self`` is mutable.

        Raise a :exc:`ValueError` if ``self`` is immutable.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el._require_mutable()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        if self._is_immutable:
            raise ValueError("object is immutable; please change a copy instead.")

    cpdef bint is_mutable(self) noexcept:
        """
        Return ``True`` if ``self`` is mutable (can be changed) and ``False``
        if it is not.

        To make this object immutable use ``self.set_immutable()``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el.is_mutable()
            False
            sage: copy(el).is_mutable()
            True
            sage: with el.clone() as el1:
            ....:      print([el.is_mutable(), el1.is_mutable()])
            [False, True]
        """
        return not self._is_immutable

    cpdef bint is_immutable(self) noexcept:
        """
        Return ``True`` if ``self`` is immutable (cannot be changed)
        and ``False`` if it is not.

        To make ``self`` immutable use ``self.set_immutable()``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el.is_immutable()
            True
            sage: copy(el).is_immutable()
            False
            sage: with el.clone() as el1:
            ....:      print([el.is_immutable(), el1.is_immutable()])
            [True, False]
        """
        return self._is_immutable

    cpdef set_immutable(self):
        """
        Makes ``self`` immutable, so it can never again be changed.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el1 = copy(el); el1.is_mutable()
            True
            sage: el1.set_immutable();  el1.is_mutable()
            False
            sage: el1[2] = 4
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        self._is_immutable = True

    cpdef _set_mutable(self):
        """
        Makes ``self`` mutable, so it can be changed.

        .. warning:: for internal use only. Casual users should make a copy
             using either the :meth:`__copy__` method or the :meth:`clone`
             protocol. Use only if you really know what you are doing. You
             should in particular make sure that you are the only owner of
             your object.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el._set_mutable(); el.is_mutable()
            True
            sage: hash(el)
            Traceback (most recent call last):
            ...
            ValueError: cannot hash a mutable object.
        """
        self._hash = 0
        self._is_immutable = False

    def __hash__(self):
        """
        Return the hash value of ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: hash(el)    # random
            -309690657
            sage: el1 = copy(el); hash(el1)
            Traceback (most recent call last):
            ...
            ValueError: cannot hash a mutable object.
        """
        if self._hash == 0:
            if not self._is_immutable:
                raise ValueError("cannot hash a mutable object.")
            else:
                self._hash = self._hash_()
        return self._hash

    cpdef ClonableElement clone(self, bint check=True):
        """
        Return a clone that is mutable copy of ``self``.

        INPUT:

        - ``check`` -- boolean indicating if ``self.check()`` must be called
          after modifications

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: with el.clone() as el1:
            ....:      el1[2] = 5
            sage: el1
            [1, 2, 5]
        """
        cdef ClonableElement res
        res = self.__copy__()
        res._needs_check = check
        return res

    def __enter__(self):
        """
        Implement the ``self`` guarding clone protocol.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el.clone().__enter__()
            [1, 2, 3]
        """
        self._require_mutable()
        return self

    def __exit__(self, typ, value, tracback):
        """
        Implement the ``self`` guarding clone protocol.

        .. NOTE:: The input argument are required by the ``with`` protocol but
           are ignored.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el1 = el.clone().__enter__()
            sage: el1.__exit__(None, None, None)
            False

        Lets try to make a broken list::

            sage: elc2 = el.clone().__enter__()
            sage: elc2[1] = 7
            sage: elc2.__exit__(None, None, None)
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing
        """
        self.set_immutable()
        if __debug__ and self._needs_check:
            self.check()
        # is there a way if check() fails to replace self by a invalid object ?
        return False


############################################################################
#       The most common case of clone object : list with constraints       #
############################################################################
cdef class ClonableArray(ClonableElement):
    """
    Array with clone protocol.

    The class of objects which are
    :class:`Element<sage.structure.element.Element>` behave as arrays
    (i.e. lists of fixed length) and implement the clone protocol. See
    :class:`ClonableElement` for details about clone protocol.

    INPUT:

    - ``parent`` -- a :class:`Parent<sage.structure.parent.Parent>`
    - ``lst`` -- list
    - ``check`` -- boolean specifying if the invariant must be checked
      using method :meth:`check`
    - ``immutable`` -- boolean (default: ``True``); whether the created element
      is immutable

    .. SEEALSO:: :class:`~sage.structure.list_clone_demo.IncreasingArray` for
                 an example of usage.

    EXAMPLES::

        sage: from sage.structure.list_clone_demo import IncreasingArrays
        sage: IA = IncreasingArrays()
        sage: ia1 = IA([1, 4, 6]); ia1
        [1, 4, 6]
        sage: with ia1.clone() as ia2:
        ....:      ia2[1] = 5
        sage: ia2
        [1, 5, 6]
        sage: with ia1.clone() as ia2:
        ....:      ia2[1] = 7
        Traceback (most recent call last):
        ...
        ValueError: array is not increasing
    """
    def __init__(self, Parent parent, lst, check=True, immutable=True):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: IncreasingArrays()([1,2,3])
            [1, 2, 3]

            sage: el = IncreasingArrays()([3,2,1])
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing

            sage: IncreasingArrays()(None)
            Traceback (most recent call last):
            ...
            TypeError: 'NoneType' object is not iterable

        You are not supposed to do the following (giving a wrong list and
        desactivating checks)::

            sage: broken = IncreasingArrays()([3,2,1], False)
        """
        self._parent = parent
        self._list = list(lst)
        self._is_immutable = immutable
        if check:
            self.check()

    def _repr_(self):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: IncreasingArrays()([1,2,3])
            [1, 2, 3]
        """
        return repr(self._list)

    def __bool__(self):
        """
        Test if ``self`` is not empty.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: bool(IncreasingArrays()([1,2,3]))
            True
            sage: bool(IncreasingArrays()([]))
            False
        """
        return bool(self._list)

    cpdef list _get_list(self):
        """
        Return the list embedded in ``self``.

        .. warning:: No copy is performed. As a consequence, modifying the
           returned list is not allowed.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el._get_list()
            [1, 2, 3]
        """
        return self._list

    cpdef _set_list(self, list lst):
        """
        Set the list embedded in ``self``.

        .. warning:: No copy is performed. Modifying the list after calling
           ``_set_list`` on it is not allowed.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el._set_list([1,4,5])
            sage: el
            [1, 4, 5]
        """
        self._list = lst

    def __len__(self):
        """
        Return the ``len`` of ``self``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: len(IncreasingArrays()([1,2,3]))
            3
        """
        return len(self._list)

    def __getitem__(self, key):
        """
        Return the ``key``-th element of ``self``.

        It also works with slice returning a python list in this case.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: IncreasingArrays()([1,2,3])[1]
            2
            sage: IncreasingArrays()([1,2,3])[7]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: IncreasingArrays()([1,2,3])[-1]
            3
            sage: IncreasingArrays()([1,2,3])[-1:]
            [3]
            sage: IncreasingArrays()([1,2,3])[:]
            [1, 2, 3]
            sage: type(IncreasingArrays()([1,2,3])[:])
            <... 'list'>
        """
        if isinstance(key, slice):
            self._list[key.start:key.stop:key.step]
        return self._list[key]

    def __setitem__(self, int key, value):
        """
        Set the ``i``-th element of ``self``.

        An exception is raised if ``self`` is immutable.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,4,10])
            sage: elc = copy(el)
            sage: elc[1] = 3; elc
            [1, 3, 4, 10]
            sage: el[1] = 3
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: elc[5] = 3
            Traceback (most recent call last):
            ...
            IndexError: list assignment index out of range
        """
        self._require_mutable()
        self._list[key] = value

    cpdef object _getitem(self, int key):
        """
        Same as :meth:`__getitem__`.

        This is much faster when used with Cython and ``key`` is known to be
        an ``int``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: IncreasingArrays()([1,2,3])._getitem(1)
            2
            sage: IncreasingArrays()([1,2,3])._getitem(5)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self._list[key]

    cpdef _setitem(self, int key, value):
        """
        Same as :meth:`__setitem__`.

        This is much faster when used with Cython and ``key`` is known to be
        an ``int``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,4])
            sage: elc = copy(el)
            sage: elc._setitem(1, 3); elc
            [1, 3, 4]
            sage: el._setitem(1, 3)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: elc._setitem(5, 5)
            Traceback (most recent call last):
            ...
            IndexError: list assignment index out of range
        """
        self._require_mutable()
        self._list[key] = value

    def __iter__(self):
        """
        Return an iterator for ``self``::

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,4])
            sage: list(iter(el))
            [1, 2, 4]

        As a consequence ``min``, ``max`` behave properly::

            sage: el = IncreasingArrays()([1,4,8])
            sage: min(el), max(el)
            (1, 8)

        TESTS::

            sage: list(iter(IncreasingArrays()([])))
            []
        """
        return iter(self._list)

    def __contains__(self, item):
        """
        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: c = IncreasingArrays()([1,2,4])
            sage: 1 in c
            True
            sage: 5 in c
            False
        """
        return item in self._list

    cpdef int index(self, x, start=None, stop=None) except -1:
        """
        Return the smallest ``k`` such that ``s[k] == x`` and ``i <= k < j``

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: c = IncreasingArrays()([1,2,4])
            sage: c.index(1)
            0
            sage: c.index(4)
            2
            sage: c.index(5)
            Traceback (most recent call last):
            ...
            ValueError: 5 is not in list
        """
        if start is None:
            return self._list.index(x)
        elif stop is None:
            return self._list.index(x, start)
        else:
            return self._list.index(x, start, stop)

    cpdef int count(self, key) except -1:
        """
        Return number of ``i``'s for which ``s[i] == key``

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: c = IncreasingArrays()([1,2,2,4])
            sage: c.count(1)
            1
            sage: c.count(2)
            2
            sage: c.count(3)
            0
        """
        return self._list.count(key)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: hash(el)    # random
            -309690657
            sage: el1 = copy(el); hash(el1)
            Traceback (most recent call last):
            ...
            ValueError: cannot hash a mutable object.
        """
        if self._hash == 0:
            if not self._is_immutable:
                raise ValueError("cannot hash a mutable object.")
            else:
                self._hash = self._hash_()
        return self._hash

    # See protocol in comment in sage/structure/element.pyx
    cpdef _richcmp_(left, right, int op):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,4])
            sage: elc = copy(el)
            sage: elc == el             # indirect doctest
            True

        ::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el1 = IncreasingArrays()([1,2,4])
            sage: el2 = IncreasingArrays()([1,2,3])
            sage: el1 == el1, el2 == el2, el1 == el2    # indirect doctest
            (True, True, False)
            sage: el1 <= el2, el1 >= el2, el2 <= el1    # indirect doctest
            (False, True, True)
        """
        cdef ClonableArray rgt = <ClonableArray>right
        return richcmp(left._list, rgt._list, op)

    cpdef ClonableArray __copy__(self):
        """
        Return a copy of ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,4])
            sage: elc = copy(el)
            sage: el[:] == elc[:]
            True
            sage: el is elc
            False

        We check that empty lists are correctly copied::

            sage: el = IncreasingArrays()([])
            sage: elc = copy(el)
            sage: el is elc
            False
            sage: bool(elc)
            False
            sage: elc.is_mutable()
            True

        We check that element with a ``__dict__`` are correctly copied::

            sage: IL = IncreasingArrays()
            sage: class myClass(IL.element_class): pass
            sage: el = myClass(IL, [])
            sage: el.toto = 2
            sage: elc = copy(el)
            sage: elc.toto
            2
        """
        cdef ClonableArray res
        cdef type t = type(self)
        res = t.__new__(t)
        res._parent = self._parent
        res._list = self._list[:]
        if HAS_DICTIONARY(self):
            res.__dict__ = self.__dict__.copy()
        return res

    cpdef check(self):
        """
        Check that ``self`` fulfill the invariants.

        This is an abstract method. Subclasses are supposed to overload
        ``check``.

        EXAMPLES::

            sage: from sage.structure.list_clone import ClonableArray
            sage: ClonableArray(Parent(), [1,2,3]) # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: this should never be called, please overload the check method
            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,4]) # indirect doctest
        """
        raise NotImplementedError("this should never be called, please overload the check method")

    cpdef long _hash_(self) except? -1:
        """
        Return the hash value of ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,3])
            sage: el._hash_()    # random
            -309711137
            sage: type(el._hash_()) == int
            True
        """
        return hash(tuple(self._list))

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: el = IncreasingArrays()([1,2,4])
            sage: loads(dumps(el))
            [1, 2, 4]
            sage: t = el.__reduce__(); t
            (<built-in function _make_array_clone>,
             (<class 'sage.structure.list_clone_demo.IncreasingArray'>,
              <sage.structure.list_clone_demo.IncreasingArrays_with_category object at ...>,
              [1, 2, 4],
              True,
              True,
              None))
            sage: t[0](*t[1])
            [1, 2, 4]
        """
        # Warning: don't pickle the hash value as it can change upon unpickling.
        if HAS_DICTIONARY(self):
            dic = self.__dict__
        else:
            dic = None
        return (_make_array_clone,
                (type(self), self._parent, self._list,
                 self._needs_check, self._is_immutable, dic))


##### Needed for unpickling #####
def _make_array_clone(clas, parent, list, needs_check, is_immutable, dic):
    """
    Helpler to unpickle :class:`list_clone` instances.

    TESTS::

        sage: from sage.structure.list_clone import _make_array_clone
        sage: from sage.structure.list_clone_demo import IncreasingArrays
        sage: ILs = IncreasingArrays()
        sage: el = _make_array_clone(ILs.element_class, ILs, [1,2,3], True, True, None)
        sage: el
        [1, 2, 3]
        sage: el == ILs([1,2,3])
        True

    We check that element with a ``__dict__`` are correctly pickled::

        sage: IL = IncreasingArrays()
        sage: class myClass(IL.element_class): pass
        sage: import __main__
        sage: __main__.myClass = myClass
        sage: el = myClass(IL, [])
        sage: el.toto = 2
        sage: elc = loads(dumps(el))
        sage: elc.toto
        2
    """
    cdef ClonableArray res
    res = <ClonableArray> clas.__new__(clas)
    res._parent = parent
    res._list = list
    res._needs_check = needs_check
    res._is_immutable = is_immutable
    if dic is not None:
        res.__dict__ = dic
    return res


############################################################################
###                      Clonable (Resizable) Lists                      ###
############################################################################
cdef class ClonableList(ClonableArray):
    """
    List with clone protocol.

    The class of objects which are
    :class:`Element<sage.structure.element.Element>` behave as lists and
    implement the clone protocol. See :class:`ClonableElement` for details
    about clone protocol.

    .. SEEALSO:: :class:`~sage.structure.list_clone_demo.IncreasingList` for
                 an example of usage.
    """
    cpdef append(self, el):
        """
        Appends ``el`` to ``self``.

        INPUT:

        - ``el`` -- any object

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: el = IncreasingLists()([1])
            sage: el.append(3)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with el.clone() as elc:
            ....:      elc.append(4)
            ....:      elc.append(6)
            sage: elc
            [1, 4, 6]
            sage: with el.clone() as elc:
            ....:      elc.append(4)
            ....:      elc.append(2)
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing
        """
        self._require_mutable()
        self._list.append(el)

    cpdef extend(self, it):
        """
        Extend ``self`` by the content of the iterable ``it``.

        INPUT:

        - ``it`` -- any iterable

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: el = IncreasingLists()([1, 4, 5, 8, 9])
            sage: el.extend((10,11))
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.

            sage: with el.clone() as elc:
            ....:      elc.extend((10,11))
            sage: elc
            [1, 4, 5, 8, 9, 10, 11]

            sage: el2 = IncreasingLists()([15, 16])
            sage: with el.clone() as elc:
            ....:      elc.extend(el2)
            sage: elc
            [1, 4, 5, 8, 9, 15, 16]

            sage: with el.clone() as elc:
            ....:      elc.extend((6,7))
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing
        """
        self._require_mutable()
        self._list.extend(it)

    cpdef insert(self, int index, el):
        """
        Inserts ``el`` in ``self`` at position ``index``.

        INPUT:

        - ``el`` -- any object
        - ``index`` -- any int

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: el = IncreasingLists()([1, 4, 5, 8, 9])
            sage: el.insert(3, 6)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with el.clone() as elc:
            ....:      elc.insert(3, 6)
            sage: elc
            [1, 4, 5, 6, 8, 9]
            sage: with el.clone() as elc:
            ....:      elc.insert(2, 6)
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing
        """
        self._require_mutable()
        self._list.insert(index, el)

    cpdef pop(self, int index=-1):
        """
        Remove ``self[index]`` from ``self`` and returns it.

        INPUT:

        - ``index`` -- integer (default: -1)

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: el = IncreasingLists()([1, 4, 5, 8, 9])
            sage: el.pop()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with el.clone() as elc:
            ....:      print(elc.pop())
            9
            sage: elc
            [1, 4, 5, 8]
            sage: with el.clone() as elc:
            ....:      print(elc.pop(2))
            5
            sage: elc
            [1, 4, 8, 9]
        """
        self._require_mutable()
        return self._list.pop(index)

    cpdef remove(self, el):
        """
        Remove the first occurrence of ``el`` from ``self``.

        INPUT:

        - ``el`` -- any object

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: el = IncreasingLists()([1, 4, 5, 8, 9])
            sage: el.remove(4)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with el.clone() as elc:
            ....:      elc.remove(4)
            sage: elc
            [1, 5, 8, 9]
            sage: with el.clone() as elc:
            ....:      elc.remove(10)
            Traceback (most recent call last):
            ...
            ValueError: list.remove(x): x not in list
        """
        self._require_mutable()
        return self._list.remove(el)

    def __setitem__(self, key, value):
        """
        Set the i-th element of ``self``.

        An exception is raised if ``self`` is immutable.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: el = IncreasingLists()([1,2,4,10,15,17])
            sage: el[1] = 3
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.

            sage: with el.clone() as elc:
            ....:      elc[3] = 7
            sage: elc
            [1, 2, 4, 7, 15, 17]

            sage: with el.clone(check=False) as elc:
            ....:      elc[1:3]  = [3,5,6,8]
            sage: elc
            [1, 3, 5, 6, 8, 10, 15, 17]
        """
        self._require_mutable()
        self._list[key] = value

    def __delitem__(self, key):
        """
        Remove the i-th element of ``self``.

        An exception is raised if ``self`` is immutable.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: el = IncreasingLists()([1, 4, 5, 8, 9])
            sage: del el[3]
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with el.clone() as elc:
            ....:      del elc[3]
            sage: elc
            [1, 4, 5, 9]
            sage: with el.clone() as elc:
            ....:      del elc[1:3]
            sage: elc
            [1, 8, 9]
        """
        self._require_mutable()
        del self._list[key]


############################################################################
###                     Clonable Arrays of int                            ##
############################################################################
cdef class ClonableIntArray(ClonableElement):
    """
    Array of integers with clone protocol.

    The class of objects which are
    :class:`Element<sage.structure.element.Element>` behave as list of int and
    implement the clone protocol. See :class:`ClonableElement` for details
    about clone protocol.


    INPUT:

    - ``parent`` -- a :class:`Parent<sage.structure.parent.Parent>`
    - ``lst`` -- list
    - ``check`` -- boolean specifying if the invariant must be checked
      using method :meth:`check`
    - ``immutable`` -- boolean (default: ``True``); whether the created element
      is immutable

    .. SEEALSO:: :class:`~sage.structure.list_clone_demo.IncreasingIntArray`
                 for an example of usage.
    """
    def __cinit__(self):
        self._len = -1
        self._list = NULL

    def __init__(self, Parent parent, lst, check=True, immutable=True):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: IncreasingIntArrays()([1,2,3])
            [1, 2, 3]
            sage: IncreasingIntArrays()((1,2,3))
            [1, 2, 3]

            sage: IncreasingIntArrays()(None)
            Traceback (most recent call last):
            ...
            TypeError: object of type 'NoneType' has no len()

            sage: el = IncreasingIntArrays()([3,2,1])
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing

            sage: el = IncreasingIntArrays()([1,2,4])
            sage: list(iter(el))
            [1, 2, 4]
            sage: list(iter(IncreasingIntArrays()([])))
            []

        You are not supposed to do the following (giving a wrong list and
        desactivating checks)::

            sage: broken = IncreasingIntArrays()([3,2,1], False)
        """
        cdef int i
        self._parent = parent

        if self._list is not NULL:
            raise ValueError("resizing is forbidden")
        self._alloc_(len(lst))
        for i in range(self._len):
            self._list[i] = lst[i]

        self._is_immutable = immutable
        if check:
            self.check()

    cpdef _alloc_(self, int size):
        """
        Allocate the array part of ``self`` for a given size.

        This can be used to initialize ``self`` without passing a list

        INPUT:

        - ``size`` -- integer

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([], check=False)
            sage: el._alloc_(3)
            sage: el._setitem(0, 1); el._setitem(1, 5); el._setitem(2, 8)
            sage: el
            [1, 5, 8]
            sage: copy(el)
            [1, 5, 8]

        TESTS::

            sage: el._alloc_(-1)
            Traceback (most recent call last):
            ...
            AssertionError: Negative size is forbidden
        """
        assert size >= 0, "Negative size is forbidden"
        self._is_immutable = False
        self._list = <int *>check_reallocarray(self._list, size, sizeof(int))
        self._len = size

    def __dealloc__(self):
        sig_free(self._list)
        self._len = -1
        self._list = NULL

    def _repr_(self):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: IncreasingIntArrays()([1,2,3])
            [1, 2, 3]
        """
        return '[' + ', '.join("%i" % self._list[i]
                               for i in range(self._len)) + ']'

    def __bool__(self):
        """
        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: bool(IncreasingIntArrays()([1,2,3]))
            True
            sage: bool(IncreasingIntArrays()([]))
            False
        """
        return self._len != 0

    def __len__(self):
        """
        Return the len of ``self``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: len(IncreasingIntArrays()([1,2,3]))
            3
        """
        return self._len

    def __iter__(self):
        """
        Iterate over the items of ``self``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: I = IncreasingIntArrays()(range(5))
            sage: I == list(range(5))
            False
            sage: list(I) == list(range(5))  # indirect doctest
            True
        """
        return iter(self.list())

    cpdef list list(self):
        """
        Convert ``self`` into a Python list.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: I = IncreasingIntArrays()(range(5))
            sage: I == list(range(5))
            False
            sage: I.list() == list(range(5))
            True
            sage: I = IncreasingIntArrays()(range(1000))
            sage: I.list() == list(range(1000))
            True
        """
        cdef int i
        cdef list L = <list> PyList_New(self._len)
        cdef object o
        for i in range(self._len):
            o = PyLong_FromLong(self._list[i])
            Py_INCREF(o)
            PyList_SET_ITEM(L, i, o)
        return L

    def __getitem__(self, key):
        """
        Return the i-th element of ``self``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,3])
            sage: el[1]
            2
            sage: el[1:2]
            [2]
            sage: el[4]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: el[-1]
            3
            sage: el[-1:]
            [3]
            sage: el[:]
            [1, 2, 3]
            sage: el[1:3]
            [2, 3]
            sage: type(el[:])
            <... 'list'>
            sage: list(el)
            [1, 2, 3]
            sage: it = iter(el); next(it), next(it)
            (1, 2)
        """
        cdef int start, stop, step, keyi
        cdef list res
        cdef slice keysl
        if isinstance(key, slice):
            keysl = <slice> key
            start, stop, step = keysl.indices(self._len)
            res = []
            for i in range(start, stop, step):
                res.append(self._getitem(i))
            return res
        keyi = <int> key
        if keyi < 0:
            keyi += self._len
        if 0 <= keyi < self._len:
            return self._list[keyi]
        else:
            raise IndexError("list index out of range")

    def __setitem__(self, int key, value):
        """
        Set the i-th element of ``self``.

        An exception is raised if ``self`` is immutable.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,4])
            sage: elc = copy(el)
            sage: elc[1] = 3; elc
            [1, 3, 4]
            sage: el[1] = 3
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        if 0 <= key < self._len:
            self._require_mutable()
            self._list[key] = value
        else:
            raise IndexError("list index out of range")

    cpdef object _getitem(self, int key):
        """
        Same as :meth:`__getitem__`.

        This is much faster when used with Cython and the index is known to be
        an int.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: IncreasingIntArrays()([1,2,3])._getitem(1)
            2
        """
        if 0 <= key < self._len:
            return self._list[key]
        else:
            raise IndexError("list index out of range")

    cpdef _setitem(self, int key, value):
        """
        Same as :meth:`__setitem__`.

        This is much faster when used with Cython and the index is known to be
        an int.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,4])
            sage: elc = copy(el)
            sage: elc._setitem(1, 3); elc
            [1, 3, 4]
            sage: el._setitem(1, 3)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        if 0 <= key < self._len:
            self._require_mutable()
            self._list[key] = value
        else:
            raise IndexError("list index out of range")

    def __contains__(self, int item):
        """
        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: c = IncreasingIntArrays()([1,2,4])
            sage: 1 in c
            True
            sage: 5 in c
            False
        """
        cdef int i
        for i in range(self._len):
            if item == self._list[i]:
                return True
        return False

    cpdef int index(self, int item) except -1:
        """
        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: c = IncreasingIntArrays()([1,2,4])
            sage: c.index(1)
            0
            sage: c.index(4)
            2
            sage: c.index(5)
            Traceback (most recent call last):
            ...
            ValueError: list.index(x): x not in list
        """
        cdef int i
        for i in range(self._len):
            if item == self._list[i]:
                return i
        raise ValueError("list.index(x): x not in list")

    # __hash__ is not properly inherited if comparison is changed
    # see <http://groups.google.com/group/cython-users/t/e89a9bd2ff20fd5a>
    def __hash__(self):
        """
        Return the hash value of ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,3])
            sage: hash(el)    # random
            -309690657
            sage: el1 = copy(el); hash(el1)
            Traceback (most recent call last):
            ...
            ValueError: cannot hash a mutable object.
        """
        if self._hash == 0:
            if not self._is_immutable:
                raise ValueError("cannot hash a mutable object.")
            else:
                self._hash = self._hash_()
        return self._hash

    # See protocol in comment in sage/structure/element.pyx
    cpdef _richcmp_(left, right, int op):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,4])
            sage: elc = copy(el)
            sage: elc == el             # indirect doctest
            True

        ::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el1 = IncreasingIntArrays()([1,2,4])
            sage: el2 = IncreasingIntArrays()([1,2,3])
            sage: el1 == el1, el2 == el2, el1 == el2    # indirect doctest
            (True, True, False)
            sage: el1 <= el2, el1 >= el2, el2 <= el1    # indirect doctest
            (False, True, True)
        """
        cdef int i, minlen, reslen
        cdef ClonableIntArray rgt = <ClonableIntArray>right
        if left is right:
            return rich_to_bool(op, 0)
        if left._list is NULL:
            if rgt._list is NULL:
                return rich_to_bool(op, 0)
            else:
                return rich_to_bool(op, -1)
        elif rgt._list is NULL:
            return rich_to_bool(op, 1)
        if left._len < rgt._len:
            minlen = left._len
            reslen = -1
        elif left._len > rgt._len:
            minlen = rgt._len
            reslen = +1
        else:
            minlen = rgt._len
            reslen = 0
        for i in range(minlen):
            if left._list[i] != rgt._list[i]:
                if left._list[i] < rgt._list[i]:
                    return rich_to_bool(op, -1)
                else:
                    return rich_to_bool(op, 1)
        return rich_to_bool(op, reslen)

    cpdef ClonableIntArray __copy__(self):
        """
        Return a copy of ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,4])
            sage: elc = copy(el)
            sage: el[:] == elc[:]
            True
            sage: el is elc
            False

        We check that void lists are correctly copied::

            sage: el = IncreasingIntArrays()([])
            sage: elc = copy(el)
            sage: el is elc
            False
            sage: bool(elc)
            True
            sage: elc.is_mutable()
            True

        We check that element with a ``__dict__`` are correctly copied::

            sage: IL = IncreasingIntArrays()
            sage: class myClass(IL.element_class): pass
            sage: el = myClass(IL, [])
            sage: el.toto = 2
            sage: elc = copy(el)
            sage: elc.toto
            2
        """
        cdef ClonableIntArray res
        cdef type t = type(self)
        res = t.__new__(t)
        res._parent = self._parent
        if self:
            res._alloc_(self._len)
            for i in range(self._len):
                res._list[i] = self._list[i]
        if HAS_DICTIONARY(self):
            res.__dict__ = self.__dict__.copy()
        return res

    cpdef check(self):
        """
        Check that ``self`` fulfill the invariants.

        This is an abstract method. Subclasses are supposed to overload
        ``check``.

        EXAMPLES::

            sage: from sage.structure.list_clone import ClonableArray
            sage: ClonableArray(Parent(), [1,2,3]) # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: this should never be called, please overload the check method
            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,4]) # indirect doctest
        """
        raise NotImplementedError("this should never be called, please overload the check method")

    cpdef long _hash_(self) except? -1:
        """
        Return the hash value of ``self``.

        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,3])
            sage: el._hash_()    # random
            -309711137
            sage: type(el._hash_()) == int
            True
        """
        cdef long hv
        if self._list == NULL:
            hv = hash(None)
        else:
            hv = hash(tuple(self))
        return hv

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: el = IncreasingIntArrays()([1,2,4])
            sage: loads(dumps(el))
            [1, 2, 4]
            sage: t = el.__reduce__(); t
            (<built-in function _make_int_array_clone>,
             (<class 'sage.structure.list_clone_demo.IncreasingIntArray'>,
              <sage.structure.list_clone_demo.IncreasingIntArrays_with_category object at ...>,
              [1, 2, 4],
              True,
              True,
              None))
            sage: t[0](*t[1])
            [1, 2, 4]
        """
        # Warning: don't pickle the hash value as it can change upon unpickling.
        if HAS_DICTIONARY(self):
            dic = self.__dict__
        else:
            dic = None
        return (_make_int_array_clone,
                (type(self), self._parent, self[:],
                 self._needs_check, self._is_immutable, dic))


##### Needed for unpickling #####
def _make_int_array_clone(clas, parent, lst, needs_check, is_immutable, dic):
    """
    Helpler to unpickle :class:`list_clone` instances.

    TESTS::

        sage: from sage.structure.list_clone import _make_int_array_clone
        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: ILs = IncreasingIntArrays()
        sage: el = _make_int_array_clone(ILs.element_class, ILs, [1,2,3], True, True, None)
        sage: el
        [1, 2, 3]
        sage: el == ILs([1,2,3])
        True

    We check that element with a ``__dict__`` are correctly pickled::

        sage: IL = IncreasingIntArrays()
        sage: class myClass(IL.element_class): pass
        sage: import __main__
        sage: __main__.myClass = myClass
        sage: el = myClass(IL, [])
        sage: el.toto = 2
        sage: elc = loads(dumps(el))
        sage: elc.toto
        2
    """
    cdef ClonableIntArray res
    res = <ClonableIntArray> clas.__new__(clas)
    ClonableIntArray.__init__(res, parent, lst, needs_check, is_immutable)
    if dic is not None:
        res.__dict__ = dic
    return res


cdef class NormalizedClonableList(ClonableList):
    """
    List with clone protocol and normal form.

    This is a subclass of :class:`ClonableList` which call a method
    :meth:`normalize` at creation and after any modification of its instance.

    .. SEEALSO:: :class:`~sage.structure.list_clone_demo.SortedList` for an
                 example of usage.

    EXAMPLES:

    We construct a sorted list through its parent::

        sage: from sage.structure.list_clone_demo import SortedLists
        sage: SL = SortedLists()
        sage: sl1 = SL([4,2,6,1]); sl1
        [1, 2, 4, 6]

    Normalization is also performed atfer modification::

        sage: with sl1.clone() as sl2:
        ....:      sl2[1] = 12
        sage: sl2
        [1, 4, 6, 12]
    """
    def __init__(self, Parent parent, lst, check=True, immutable=True):
        r"""
        TESTS::

            sage: from sage.structure.list_clone_demo import SortedList, SortedLists
            sage: SortedList(SortedLists(), [2,3,1])
            [1, 2, 3]
        """
        ClonableList.__init__(self, parent, lst, False, False)
        self.normalize()
        self._is_immutable = immutable
        if check:
            self.check()

    def __exit__(self, typ, value, tracback):
        r"""
        TESTS::

            sage: from sage.structure.list_clone_demo import SortedList, SortedLists
            sage: l = SortedList(SortedLists(), [2,3,1], immutable=False); l
            [1, 2, 3]
            sage: l[1] = 5; l
            [1, 5, 3]
            sage: l.__exit__(None, None, None)
            False
            sage: l
            [1, 3, 5]
        """
        self.normalize()
        return ClonableList.__exit__(self, typ, value, tracback)

    cpdef normalize(self):
        """
        Normalize ``self``.

        This is an abstract method. Subclasses are supposed to overload
        :meth:`normalize`. The call ``self.normalize()`` is supposed to

        - call ``self._require_mutable()`` to check that ``self`` is in a
          proper mutable state
        - modify ``self`` to put it in a normal form

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import SortedList, SortedLists
            sage: l = SortedList(SortedLists(), [2,3,2], False, False)
            sage: l
            [2, 2, 3]
            sage: l.check()
            Traceback (most recent call last):
            ...
            ValueError: list is not strictly increasing
        """
        raise NotImplementedError("This should never be called, please overload the normalize method")
