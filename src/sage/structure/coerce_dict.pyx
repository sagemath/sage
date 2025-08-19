# sage_setup: distribution = sagemath-objects
"""
Containers for storing coercion data

This module provides :class:`TripleDict` and :class:`MonoDict`. These are
structures similar to :class:`~weakref.WeakKeyDictionary` in Python's weakref
module, and are optimized for lookup speed. The keys for :class:`TripleDict`
consist of triples (k1,k2,k3) and are looked up by identity rather than
equality. The keys are stored by weakrefs if possible. If any one of the
components k1, k2, k3 gets garbage collected, then the entry is removed from
the :class:`TripleDict`.

Key components that do not allow for weakrefs are stored via a normal
refcounted reference. That means that any entry stored using a triple
(k1,k2,k3) so that none of the k1,k2,k3 allows a weak reference behaves
as an entry in a normal dictionary: Its existence in :class:`TripleDict`
prevents it from being garbage collected.

That container currently is used to store coercion and conversion maps between
two parents (:issue:`715`) and to store homsets of pairs of objects of a
category (:issue:`11521`). In both cases, it is essential that the parent
structures remain garbage collectable, it is essential that the data access is
faster than with a usual :class:`~weakref.WeakKeyDictionary`, and we enforce
the "unique parent condition" in Sage (parent structures should be identical
if they are equal).

:class:`MonoDict` behaves similarly, but it takes a single item as a key. It
is used for caching the parents which allow a coercion map into a fixed other
parent (:issue:`12313`).

By :issue:`14159`, :class:`MonoDict` and :class:`TripleDict` can be optionally
used with weak references on the values.

Note that this kind of dictionary is also used for caching actions and
coerce maps. In previous versions of Sage, the cache was by strong
references and resulted in a memory leak in the following example.
However, this leak was fixed by :issue:`715`, using weak references::

    sage: # needs sage.combinat sage.modules sage.rings.finite_rings
    sage: K.<t> = GF(2^55)
    sage: for i in range(20):
    ....:     a = K.random_element()
    ....:     E = EllipticCurve(j=a)
    ....:     P = E.random_point()
    ....:     Q = 2*P
    sage: L = [Partitions(n) for n in range(200)]  # purge strong cache in CachedRepresentation
    sage: import gc
    sage: n = gc.collect()
    sage: from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
    sage: LE = [x for x in gc.get_objects() if isinstance(x, EllipticCurve_finite_field)]
    sage: len(LE)
    1
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                     2012 Simon King <simon.king@uni-jena.de>
#                     2013 Nils Bruin <nbruin@sfu.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport cython
from cpython.object cimport *
from cpython.ref cimport Py_XINCREF, Py_XDECREF, Py_CLEAR
from cpython.tuple cimport PyTuple_New
from cpython.weakref cimport PyWeakref_GetObject, PyWeakref_GET_OBJECT
from cysignals.memory cimport check_calloc, sig_free

cdef extern from "Python.h":
    void PyTuple_SET_ITEM(object tuple, Py_ssize_t index, PyObject* item)

cdef extern from "../cpython/pyx_visit.h":
    void Py_VISIT3(PyObject*, visitproc, void*)

cdef type KeyedRef, ref
from weakref import KeyedRef, ref

cdef inline bint is_dead_keyedref(x) noexcept:
    """
    Check whether ``x`` is a ``KeyedRef`` which is dead.
    """
    if type(x) is not KeyedRef:
        return False
    return PyWeakref_GET_OBJECT(x) is <PyObject*>None


# Unique sentinel to indicate a deleted cell
cdef object dummy = object()
cdef PyObject* deleted_key = <PyObject*>dummy


cdef inline bint valid(PyObject* obj) noexcept:
    """
    Check whether ``obj`` points to a valid object
    """
    return obj is not NULL and obj is not deleted_key


@cython.freelist(256)
cdef class ObjectWrapper:
    """
    A simple fast wrapper around a Python object. This is like a
    1-element tuple except that it does not keep a reference to the
    wrapped object.
    """
    cdef PyObject* obj


cdef inline ObjectWrapper wrap(obj):
    """
    Wrap a given Python object in an :class:`ObjectWrapper`.
    """
    cdef ObjectWrapper w = <ObjectWrapper>(ObjectWrapper.__new__(ObjectWrapper))
    w.obj = <PyObject*>obj
    return w


cdef inline PyObject* unwrap(w) except? NULL:
    """
    Return the object wrapped by an :class:`ObjectWrapper`.
    """
    return (<ObjectWrapper?>w).obj


cdef extract_mono_cell(mono_cell* cell):
    """
    Take the refcounted components from a mono_cell, put them in a
    tuple and return it. The mono_cell itself is marked as "freed".
    The refcounts originally accounting for the presence in the
    mono_cell now account for the presence in the returned tuple,
    which steals those references.

    The returned result is only used to throw away: an advantage is
    that the containing tuple participates in CPython's trashcan,
    which prevents stack overflow on large dereffing cascades.

    A slight disadvantage is that this routine needs to allocate a
    tuple (mainly just to be thrown away)
    """
    assert valid(cell.key_id)
    t = PyTuple_New(2)
    PyTuple_SET_ITEM(t, 0, cell.key_weakref)
    PyTuple_SET_ITEM(t, 1, cell.value)
    cell.key_id = deleted_key
    cell.key_weakref = NULL
    cell.value = NULL
    return t


cdef extract_triple_cell(triple_cell* cell):
    # See extract_mono_cell for documentation
    assert valid(cell.key_id1)
    t = PyTuple_New(4)
    PyTuple_SET_ITEM(t, 0, cell.key_weakref1)
    PyTuple_SET_ITEM(t, 1, cell.key_weakref2)
    PyTuple_SET_ITEM(t, 2, cell.key_weakref3)
    PyTuple_SET_ITEM(t, 3, cell.value)
    cell.key_id1 = deleted_key
    cell.key_id2 = NULL
    cell.key_id3 = NULL
    cell.key_weakref1 = NULL
    cell.key_weakref2 = NULL
    cell.key_weakref3 = NULL
    cell.value = NULL
    return t


cdef class MonoDictEraser:
    """
    Erase items from a :class:`MonoDict` when a weak reference becomes
    invalid.

    This is of internal use only. Instances of this class will be passed as a
    callback function when creating a weak reference.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import MonoDict
        sage: class A: pass
        sage: a = A()
        sage: M = MonoDict()
        sage: M[a] = 1
        sage: len(M)
        1
        sage: del a
        sage: import gc
        sage: n = gc.collect()
        sage: len(M)    # indirect doctest
        0

    AUTHOR:

    - Simon King (2012-01)
    - Nils Bruin (2013-11)
    """
    cdef D

    def __init__(self, D):
        """
        INPUT:

        - ``D`` -- a :class:`MonoDict`

        EXAMPLES::

            sage: k = set([1])
            sage: D = sage.structure.coerce_dict.MonoDict([(k,1)])
            sage: len(D)
            1
            sage: del k
            sage: len(D) # indirect doctest
            0
        """
        self.D = ref(D)

    def __call__(self, r):
        """
        INPUT:

        - ``r`` -- a weak reference with key

        For internal use only.

        EXAMPLES::

            sage: k = set([1])
            sage: D = sage.structure.coerce_dict.MonoDict([(k,1)])
            sage: len(D)
            1
            sage: del k
            sage: len(D) # indirect doctest
            0
        """
        cdef MonoDict md = <MonoDict>PyWeakref_GetObject(self.D)
        if md is None or not md.mask:
            return
        cdef mono_cell* cursor = md.lookup(unwrap(r.key))
        cdef PyObject* r_ = <PyObject*>r
        if valid(cursor.key_id):
            if cursor.key_weakref is r_ or cursor.value is r_:
                L = extract_mono_cell(cursor)
                md.used -= 1
            else:
                raise AssertionError("MonoDictEraser: key match but no weakref match")


cdef class MonoDict:
    """
    This is a hashtable specifically designed for (read) speed in
    the coercion model.

    It differs from a python WeakKeyDictionary in the following important ways:

       - Comparison is done using the 'is' rather than '==' operator.
       - Only weak references to the keys are stored if at all possible.
         Keys that do not allow for weak references are stored with a normal
         refcounted reference.
       - The callback of the weak references is safe against recursion, see below.

    There are special cdef set/get methods for faster access.
    It is bare-bones in the sense that not all dictionary methods are
    implemented.

    IMPLEMENTATION:

    It is implemented as a hash table with open addressing, similar to python's
    dict.

    INPUT:

    - ``data`` -- (optional) iterable defining initial data, as dict or
      iterable of (key, value) pairs

    - ``weak_values`` -- boolean (default: ``False``); if it is
      ``True``, weak references to the values in this dictionary will be used,
      when possible

    EXAMPLES::

        sage: from sage.structure.coerce_dict import MonoDict
        sage: L = MonoDict()
        sage: a = 'a'; b = 'ab'; c = '-15'
        sage: L[a] = 1
        sage: L[b] = 2
        sage: L[c] = 3

    The key is expected to be a unique object. Hence, the item stored for ``c``
    cannot be obtained by providing another equal string::

        sage: L[a]
        1
        sage: L[b]
        2
        sage: L[c]
        3
        sage: L['-15']
        Traceback (most recent call last):
        ...
        KeyError: '-15'

    Not all features of Python dictionaries are available, but iteration over
    the dictionary items is possible::

        sage: sorted(L.items())
        [('-15', 3), ('a', 1), ('ab', 2)]
        sage: del L[c]
        sage: sorted(L.items())
        [('a', 1), ('ab', 2)]
        sage: len(L)
        2
        sage: for i in range(1000):
        ....:     L[i] = i
        sage: len(L)
        1002
        sage: L['a']
        1
        sage: L['c']
        Traceback (most recent call last):
        ...
        KeyError: 'c'

    TESTS:

    Here, we demonstrate the use of weak values::

        sage: M = MonoDict()
        sage: MW = MonoDict(weak_values=True)
        sage: class Foo: pass
        sage: a = Foo()
        sage: b = Foo()
        sage: k = 1
        sage: M[k] = a
        sage: MW[k] = b
        sage: M[k] is a
        True
        sage: MW[k] is b
        True
        sage: k in M
        True
        sage: k in MW
        True

    While ``M`` uses a strong reference to ``a``, ``MW`` uses a *weak*
    reference to ``b``, and after deleting ``b``, the corresponding item of
    ``MW`` will be removed during the next garbage collection::

        sage: import gc
        sage: del a,b
        sage: _ = gc.collect()
        sage: k in M
        True
        sage: k in MW
        False
        sage: len(MW)
        0
        sage: len(M)
        1

   Note that ``MW`` also accepts values that do not allow for weak references::

        sage: MW[k] = int(5)
        sage: MW[k]
        5

    The following demonstrates that :class:`MonoDict` is safer than
    :class:`~weakref.WeakKeyDictionary` against recursions created by nested
    callbacks; compare :issue:`15069` (the mechanism used now is different, though)::

        sage: M = MonoDict()
        sage: class A: pass
        sage: a = A()
        sage: prev = a
        sage: for i in range(1000):
        ....:     newA = A()
        ....:     M[prev] = newA
        ....:     prev = newA
        sage: len(M)
        1000
        sage: del a
        sage: len(M)
        0

    The corresponding example with a Python :class:`weakref.WeakKeyDictionary`
    would result in a too deep recursion during deletion of the dictionary
    items::

        sage: import weakref
        sage: M = weakref.WeakKeyDictionary()
        sage: a = A()
        sage: prev = a
        sage: for i in range(1000):
        ....:     newA = A()
        ....:     M[prev] = newA
        ....:     prev = newA
        sage: len(M)
        1000

    Check that also in the presence of circular references, :class:`MonoDict`
    gets properly collected::

        sage: import gc
        sage: def count_type(T):
        ....:     return len([c for c in gc.get_objects() if isinstance(c,T)])
        sage: gc.freeze()  # so that gc.collect() only deals with our trash
        sage: N = count_type(MonoDict)
        sage: for i in range(100):
        ....:     V = [MonoDict({"id":j+100*i}) for j in range(100)]
        ....:     n = len(V)
        ....:     for i in range(n): V[i][V[(i+1)%n]] = (i+1)%n
        ....:     del V
        ....:     _ = gc.collect()
        ....:     assert count_type(MonoDict) == N
        sage: count_type(MonoDict) == N
        True
        sage: gc.unfreeze()

    AUTHORS:

    - Simon King (2012-01)
    - Nils Bruin (2012-08)
    - Simon King (2013-02)
    - Nils Bruin (2013-11)
    """
    cdef mono_cell* lookup(self, PyObject* key) noexcept:
        """
        Return a pointer to where ``key`` should be stored in this
        :class:`MonoDict`.

        This routine is used for all cases where a (potential) spot for
        a key is looked up. The returned value is a pointer into the dictionary
        store that either contains an entry with the requested key or a free spot
        where an entry for that key should go.
        """
        assert valid(key)

        cdef size_t mask = self.mask
        cdef mono_cell* table = self.table
        cdef mono_cell* first_deleted = NULL

        # Use the memory location of the key as starting point for our
        # hash.
        cdef size_t h = <size_t>key

        # The size of a Python object is at least 2 * sizeof(size_t).
        # Therefore, we don't lose any information by dividing by that.
        # Instead, the lower order bits become more interesting.
        h //= 2 * sizeof(size_t)

        # Bring some higher-order bits in with this permutation.
        cdef size_t i = (h >> 8) ^ h

        cdef size_t perturb = h

        # The probing algorithm is heavily inspired by Python dicts.
        # There is always at least one NULL entry in the store, and the
        # probe sequence eventually covers the entire store (see Theorem
        # below), so the loop below does terminate. Since this loop does
        # not change any refcounts, we know that table will not change
        # during iteration.

        # Theorem: when iterating the function i -> 5*i + 1, every
        # element of Z/(2^n Z) is reached.
        # Proof: define f(x) = 4*x + 1. Then f(5*i + 1) = 5*f(i).
        # Therefore, the iteration is really a transformation of
        # i -> 5*i on the group of 1 mod 4 elements of (Z/2^(n+2) Z).
        # It is a well known fact that this is a cyclic group generated
        # by 5 (or any element which is 5 mod 8).
        cdef mono_cell* cursor
        while True:
            cursor = &(table[i & mask])
            perturb >>= 5
            if cursor.key_id is key:
                return cursor
            elif cursor.key_id is NULL:
                return first_deleted or cursor
            elif cursor.key_id is deleted_key:
                if first_deleted is NULL:
                    first_deleted = cursor
            i = (5*i + 1) + perturb

    cdef int resize(self) except -1:
        """
        Resize dictionary. That can also mean shrink! Size is always a power of 2.
        """
        cdef mono_cell* old_table = self.table
        cdef size_t old_mask = self.mask
        cdef size_t newsize = 8
        cdef size_t minsize = 2 * self.used
        cdef mono_cell* cursor
        cdef mono_cell* entry
        while newsize < minsize:
            newsize *= 2
        cdef mono_cell* table = <mono_cell*>check_calloc(newsize, sizeof(mono_cell))

        # We are done with memory activity. We can move the new (empty)
        # table into place:
        self.table = table
        self.mask = newsize - 1
        self.used = 0
        self.fill = 0

        # We now move all entries over. We are not changing any
        # refcounts here, so this is a very tight loop that doesn't need
        # to worry about tables changing.
        cdef size_t i
        for i in range(old_mask + 1):
            entry = &(old_table[i])
            if valid(entry.key_id):
                cursor = self.lookup(entry.key_id)
                assert cursor.key_id is NULL
                cursor[0] = entry[0]
                self.used += 1
                self.fill += 1
        sig_free(old_table)

    def __cinit__(self):
        """
        Setup basic data structure.

        TESTS::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: TripleDict.__new__(TripleDict)
            <sage.structure.coerce_dict.TripleDict object at ...>
        """
        cdef size_t newsize = 8
        # The order is important here: the object must be in a
        # consistent state even if exceptions are raised.
        self.eraser = MonoDictEraser(self)
        self.table = <mono_cell*>check_calloc(newsize, sizeof(mono_cell))
        self.mask = newsize - 1
        self.used = 0
        self.fill = 0

    def __init__(self, data=None, *, weak_values=False):
        """
        Create a special dict using singletons for keys.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: a = 'a'
            sage: L[a] = 1
            sage: L[a]
            1
            sage: L = MonoDict({a: 1})
            sage: L[a]
            1
            sage: L = MonoDict([(a, 1)])
            sage: L[a]
            1
        """
        self.weak_values = weak_values
        if data:
            try:
                data = data.items()
            except AttributeError:
                pass
            for k, v in data:
                self.set(k,v)

    def __dealloc__(self):
        MonoDict_clear(self)
        sig_free(self.table)

    def __len__(self):
        """
        The number of items in ``self``.
        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a] = 1
            sage: L[a] = -1 # re-assign
            sage: L[b] = 1
            sage: L[c] = None
            sage: len(L)
            3
        """
        return self.used

    def __contains__(self, k):
        """
        Test if the dictionary contains a given key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: a = 'a'; b = 'ab'; c = 15
            sage: L[a] = 1
            sage: L[b] = 2
            sage: L[c] = 3
            sage: c in L         # indirect doctest
            True

        The keys are compared by identity, not by equality. Hence, we have::

            sage: c == 15
            True
            sage: 15 in L
            False
        """
        cdef mono_cell* cursor = self.lookup(<PyObject*>k)
        if not valid(cursor.key_id):
            return False
        if not self.weak_values:
            return True
        value = <object>cursor.value
        return not is_dead_keyedref(value)

    def __getitem__(self, k):
        """
        Get the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: a = 'a'; b = 'b'; c = 15
            sage: L[a] = 1
            sage: L[b] = 2
            sage: L[c] = 3
            sage: L[c]                  # indirect doctest
            3

        Note that the keys are supposed to be unique::

            sage: c == 15
            True
            sage: c is 15
            False
            sage: L[15]
            Traceback (most recent call last):
            ...
            KeyError: 15
        """
        return self.get(k)

    cdef get(self, k):
        cdef mono_cell* cursor = self.lookup(<PyObject*>k)
        if not valid(cursor.key_id):
            raise KeyError(k)
        # We need to check that the value is a live reference.
        # Items with dead references (in the key or value) are deleted
        # from the MonoDict by the MonoDictEraser. However, if we are
        # in the middle of a deallocation, we may see a dead reference
        # for the value. This cannot happen for the key: we are passed a
        # strong reference to the key as argument of this function, so
        # we know that it's alive.
        value = <object>cursor.value
        if type(value) is KeyedRef:
            value = <object>PyWeakref_GET_OBJECT(value)
            if value is None:
                raise KeyError(k)
        return value

    def __setitem__(self, k, value):
        """
        Set the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: a = 'a'
            sage: L[a] = -1   # indirect doctest
            sage: L[a]
            -1
            sage: L[a] = 1
            sage: L[a]
            1
            sage: len(L)
            1
        """
        self.set(k, value)

    cdef int set(self, k, value) except -1:
        cdef mono_cell entry
        cdef PyObject* old_value
        cdef bint maybe_resize = False
        entry.key_id = <PyObject*>k
        if self.weak_values:
            wrap_k = wrap(k)
            try:
                value_store = KeyedRef(value, self.eraser, wrap_k)
                entry.value = <PyObject*>value_store
            except TypeError:
                entry.value = <PyObject*>value
        else:
            entry.value = <PyObject*>value
        Py_XINCREF(entry.value)
        cursor = self.lookup(<PyObject*>k)
        if not valid(cursor.key_id):
            self.used += 1
            if cursor.key_id is NULL:
                self.fill += 1
                maybe_resize = True
            if not self.weak_values:
                wrap_k = wrap(k)
            try:
                key_store = KeyedRef(k, self.eraser, wrap_k)
                entry.key_weakref = <PyObject*>key_store
            except TypeError:
                entry.key_weakref = <PyObject*>k
            Py_XINCREF(entry.key_weakref)

            # We are taking a bit of a gamble here: we're assuming the
            # dictionary has not been resized (otherwise cursor might
            # not be a valid location anymore). The only way in which
            # that could happen is if the allocation activity above
            # forced a GC that triggered code that *adds* entries to
            # this dictionary: the dictionary can only get reshaped if
            # self.fill increases (as happens below). Note that we're
            # holding a strong ref to the dict itself, so that's not
            # liable to disappear. For the truly paranoid: we could
            # detect a change by checking if self.table has changed
            # value.
            cursor[0] = entry

            if maybe_resize and 3*self.fill > 2*self.mask:
                self.resize()
        else:
            old_value = cursor.value
            cursor.value = entry.value
            Py_XDECREF(old_value)

    def __delitem__(self, k):
        """
        Delete the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: a = 15
            sage: L[a] = -1
            sage: len(L)
            1

        Note that the keys are unique, hence using a key that is equal but not
        identical to a results in an error::

            sage: del L[15]
            Traceback (most recent call last):
            ...
            KeyError: 15
            sage: a in L
            True
            sage: del L[a]
            sage: len(L)
            0
            sage: a in L
            False
        """
        cdef mono_cell* cursor = self.lookup(<PyObject*>k)
        if not valid(cursor.key_id):
            raise KeyError(k)
        L = extract_mono_cell(cursor)
        self.used -= 1

    def items(self):
        """
        Iterate over the ``(key, value)`` pairs of this :class:`MonoDict`.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: L[1] = None
            sage: L[2] = True
            sage: L.items()
            <...generator object at ...>
            sage: sorted(L.items())
            [(1, None), (2, True)]
        """
        # Iteration is tricky because the table could change from under
        # us. The following iterates properly if the dictionary does
        # not get resized, which is guaranteed if no NEW entries in the
        # dictionary are introduced. At least we make sure to get our
        # data fresh from "self" every iteration, so that at least we're
        # not reading random memory. If the dictionary changes, it's not
        # guaranteed you get to see any particular entry.
        cdef size_t i = 0
        while i <= self.mask:
            cursor = &(self.table[i])
            i += 1
            if valid(cursor.key_id):
                key = <object>(cursor.key_weakref)
                value = <object>(cursor.value)
                if type(key) is KeyedRef:
                    key = <object>PyWeakref_GET_OBJECT(key)
                    if key is None:
                        print("found defunct key")
                        continue
                if type(value) is KeyedRef:
                    value = <object>PyWeakref_GET_OBJECT(value)
                    if value is None:
                        print("found defunct value")
                        continue
                yield (key, value)

    def copy(self):
        """
        Return a copy of this :class:`MonoDict` as Python dict.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: L[1] = 42
            sage: L.copy()
            {1: 42}
        """
        return dict(self.items())

    def __reduce__(self):
        """
        Note that we don't expect equality as this class concerns itself with
        object identity rather than object equality.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict()
            sage: L[1] = True
            sage: loads(dumps(L)) == L
            False
            sage: list(loads(dumps(L)).items())
            [(1, True)]
        """
        return MonoDict, (self.copy(),)

# The Cython supplied tp_traverse and tp_clear do not take the
# dynamically allocated table into account, so we have to supply our
# own. The only additional link to follow (that Cython does pick up
# and we have to replicate here) is the "eraser" which in its closure
# stores a reference back to the dictionary itself (meaning that
# MonoDicts only disappear on cyclic GC).
cdef int MonoDict_traverse(MonoDict self, visitproc visit, void* arg) noexcept:
    if not self.mask:
        return 0
    Py_VISIT3(<PyObject*>self.eraser, visit, arg)
    cdef size_t i
    for i in range(self.mask + 1):
        cursor = &self.table[i]
        if valid(cursor.key_id):
            Py_VISIT3(cursor.key_weakref, visit, arg)
            Py_VISIT3(cursor.value, visit, arg)


cdef int MonoDict_clear(MonoDict self) noexcept:
    """
    We clear a monodict by taking first taking away the table before
    dereffing its contents. That shortcuts callbacks, so we deref the
    entries straight here. That means this code does not participate in
    Python's trashcan the way that deletion code based on
    extract_mono_cell does, so there is probably a way this code can be
    used to overflow the C stack. It would have to be a pretty devious
    example, though.
    """
    if not self.mask:
        return 0
    cdef size_t mask = self.mask
    self.mask = 0  # Setting mask to 0 immediately prevents recursion
    self.used = 0
    self.fill = 0
    # Set self.eraser to None safely
    cdef object eraser = self.eraser
    self.eraser = None
    for i in range(mask+1):
        cursor = &(self.table[i])
        if valid(cursor.key_id):
            cursor.key_id = deleted_key
            Py_CLEAR(cursor.key_weakref)
            Py_CLEAR(cursor.value)


(<PyTypeObject*>MonoDict).tp_traverse = <traverseproc>(&MonoDict_traverse)
(<PyTypeObject*>MonoDict).tp_clear = <inquiry>(&MonoDict_clear)


cdef class TripleDictEraser:
    """
    Erases items from a :class:`TripleDict` when a weak reference becomes
    invalid.

    This is of internal use only. Instances of this class will be passed as a
    callback function when creating a weak reference.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import TripleDict
        sage: class A: pass
        sage: a = A()
        sage: T = TripleDict()
        sage: T[a,ZZ,None] = 1
        sage: T[ZZ,a,1] = 2
        sage: T[a,a,ZZ] = 3
        sage: len(T)
        3
        sage: del a
        sage: import gc
        sage: n = gc.collect()
        sage: len(T) # indirect doctest
        0

    AUTHOR:

    - Simon King (2012-01)
    - Nils Bruin (2013-11)
    """
    cdef D

    def __init__(self, D):
        """
        INPUT:

        - ``D`` -- a :class:`TripleDict`. For internal use only.

        EXAMPLES::

            sage: D = sage.structure.coerce_dict.TripleDict()
            sage: k = set([1])
            sage: D[k,1,1] = 1
            sage: len(D)
            1
            sage: del k
            sage: len(D) # indirect doctest
            0
        """
        self.D = ref(D)

    def __call__(self, r):
        """
        INPUT:

        - ``r`` -- a weak reference with key

        For internal use only.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: class A: pass
            sage: a = A()
            sage: T = TripleDict()
            sage: T[a,ZZ,None] = 1
            sage: T[ZZ,a,1] = 2
            sage: T[a,a,ZZ] = 3
            sage: len(T)
            3
            sage: del a
            sage: import gc
            sage: n = gc.collect()
            sage: len(T)    # indirect doctest
            0
        """
        cdef TripleDict td = <TripleDict>PyWeakref_GetObject(self.D)
        if td is None or not td.mask:
            return
        k1, k2, k3 = r.key
        cdef triple_cell* cursor = td.lookup(unwrap(k1), unwrap(k2), unwrap(k3))
        cdef PyObject* r_ = <PyObject*>r
        if valid(cursor.key_id1):
            if (cursor.key_weakref1 is r_ or
                    cursor.key_weakref2 is r_ or
                    cursor.key_weakref3 is r_ or
                    cursor.value is r_):
                L = extract_triple_cell(cursor)
                td.used -= 1
            else:
                raise AssertionError("TripleDictEraser: key match but no weakref match")


cdef class TripleDict:
    """
    This is a hashtable specifically designed for (read) speed in
    the coercion model.

    It differs from a python dict in the following important ways:

       - All keys must be sequence of exactly three elements. All sequence
         types (tuple, list, etc.) map to the same item.

       - Any of the three key components that support weak-refs are stored
         via a weakref. If any of these components gets garbage collected
         then the entire entry is removed. In that sense, this structure
         behaves like a nested :class:`~weakref.WeakKeyDictionary`.

       - Comparison is done using the 'is' rather than '==' operator.

    There are special cdef set/get methods for faster access.
    It is bare-bones in the sense that not all dictionary methods are
    implemented.

    INPUT:

    - ``data`` -- (optional) iterable defining initial data, as dict or
      iterable of (key, value) pairs

    - ``weak_values`` -- boolean (default: ``False``); if it is
      ``True``, weak references to the values in this dictionary will be used,
      when possible

    IMPLEMENTATION:

    It is implemented as a hash table with open addressing, similar to python's
    dict.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import TripleDict
        sage: L = TripleDict()
        sage: a = 'a'; b = 'b'; c = 'c'
        sage: L[a,b,c] = 1
        sage: L[a,b,c]
        1
        sage: L[c,b,a] = -1
        sage: sorted(L.items())
        [(('a', 'b', 'c'), 1), (('c', 'b', 'a'), -1)]
        sage: del L[a,b,c]
        sage: list(L.items())
        [(('c', 'b', 'a'), -1)]
        sage: len(L)
        1
        sage: for i in range(1000):
        ....:     L[i,i,i] = i
        sage: len(L)
        1001
        sage: L = TripleDict(L)
        sage: L[c,b,a]
        -1
        sage: L[a,b,c]
        Traceback (most recent call last):
        ...
        KeyError: ('a', 'b', 'c')
        sage: L[a]
        Traceback (most recent call last):
        ...
        KeyError: 'a'
        sage: L[a] = 1
        Traceback (most recent call last):
        ...
        KeyError: 'a'

    TESTS:

    Here, we demonstrate the use of weak values::

        sage: class Foo: pass
        sage: T = TripleDict()
        sage: TW = TripleDict(weak_values=True)
        sage: a = Foo()
        sage: b = Foo()
        sage: k = 1
        sage: T[a,k,k]=1
        sage: T[k,a,k]=2
        sage: T[k,k,a]=3
        sage: T[k,k,k]=a
        sage: TW[b,k,k]=1
        sage: TW[k,b,k]=2
        sage: TW[k,k,b]=3
        sage: TW[k,k,k]=b
        sage: len(T)
        4
        sage: len(TW)
        4
        sage: (k,k,k) in T
        True
        sage: (k,k,k) in TW
        True
        sage: T[k,k,k] is a
        True
        sage: TW[k,k,k] is b
        True

    Now, ``T`` holds a strong reference to ``a``, namely in ``T[k,k,k]``. Hence,
    when we delete ``a``, *all* items of ``T`` survive::

        sage: import gc
        sage: del a
        sage: _ = gc.collect()
        sage: len(T)
        4

    Only when we remove the strong reference, the items become collectable::

        sage: del T[k,k,k]
        sage: _ = gc.collect()
        sage: len(T)
        0

    The situation is different for ``TW``, since it only holds *weak*
    references to ``a``. Therefore, all items become collectable after
    deleting ``a``::

        sage: del b
        sage: _ = gc.collect()
        sage: len(TW)
        0

    AUTHORS:

    - Robert Bradshaw, 2007-08

    - Simon King, 2012-01

    - Nils Bruin, 2012-08

    - Simon King, 2013-02

    - Nils Bruin, 2013-11
    """
    cdef triple_cell* lookup(self, PyObject* key1, PyObject* key2, PyObject* key3) noexcept:
        """
        Return a pointer to where ``(key1, key2, key3)`` should be
        stored in this :class:`MonoDict`.

        This routine is used for all cases where a (potential) spot for
        a key is looked up. The returned value is a pointer into the dictionary
        store that either contains an entry with the requested key or a free spot
        where an entry for that key should go.
        """
        cdef size_t mask = self.mask
        cdef triple_cell* table = self.table
        cdef triple_cell* first_deleted = NULL

        # A random linear combination of the memory locations of the keys
        cdef size_t C2 = 0x7de83cbb
        cdef size_t C3 = 0x32354bf3
        cdef size_t h = (<size_t>key1) + C2*(<size_t>key2) + C3*(<size_t>key3)

        # See MonoDict.lookup() for comments about the algorithm
        h //= 2 * sizeof(size_t)

        cdef size_t i = (h >> 8) ^ h
        cdef size_t perturb = h

        cdef triple_cell* cursor
        while True:
            cursor = &(table[i & mask])
            perturb >>= 5
            if cursor.key_id1 is key1:
                if cursor.key_id2 is key2 and cursor.key_id3 is key3:
                    return cursor
            elif cursor.key_id1 is NULL:
                return first_deleted or cursor
            elif cursor.key_id1 is deleted_key:
                if first_deleted is NULL:
                    first_deleted = cursor
            i = (5*i + 1) + perturb

    cdef int resize(self) except -1:
        cdef triple_cell* old_table = self.table
        cdef size_t old_mask = self.mask
        cdef size_t newsize = 8
        cdef size_t minsize = 2 * self.used
        cdef triple_cell* cursor
        cdef triple_cell* entry
        while newsize < minsize:
            newsize *= 2
        cdef triple_cell* table = <triple_cell*>check_calloc(newsize, sizeof(triple_cell))
        self.table = table
        self.mask = newsize - 1
        self.used = 0
        self.fill = 0
        cdef size_t i
        for i in range(old_mask + 1):
            entry = &(old_table[i])
            if valid(entry.key_id1):
                cursor = self.lookup(entry.key_id1, entry.key_id2, entry.key_id3)
                assert cursor.key_id1 is NULL
                cursor[0] = entry[0]
                self.used +=1
                self.fill +=1
        sig_free(old_table)

    def __cinit__(self):
        """
        Setup basic data structure.

        TESTS::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: MonoDict.__new__(MonoDict)
            <sage.structure.coerce_dict.MonoDict object at ...>
        """
        cdef size_t newsize = 8
        # The order is important here: the object must be in a
        # consistent state even if exceptions are raised.
        self.eraser = TripleDictEraser(self)
        self.table = <triple_cell*>check_calloc(newsize, sizeof(triple_cell))
        self.mask = newsize - 1
        self.used = 0
        self.fill = 0

    def __init__(self, data=None, *, weak_values=False):
        """
        Create a special dict using triples for keys.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c]
            1
            sage: key = ("x", "y", "z")
            sage: L = TripleDict([(key, 42)])
            sage: L[key]
            42
            sage: L = TripleDict({key: 42})
            sage: L[key]
            42
        """
        self.weak_values = weak_values
        if data:
            try:
                data = data.items()
            except AttributeError:
                pass
            for k, v in data:
                self[k] = v

    def __dealloc__(self):
        TripleDict_clear(self)
        sig_free(self.table)

    def __len__(self):
        """
        The number of items in ``self``.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c] = -1 # re-assign
            sage: L[a,c,b] = 1
            sage: L[a,a,None] = None
            sage: len(L)
            3
        """
        return self.used

    def __contains__(self, k):
        """
        Test if the dictionary contains a given key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: a = 'a'; b = 'ab'; c = 15
            sage: L[a,b,c] = 123
            sage: (a,b,c) in L         # indirect doctest
            True

        The keys are compared by identity, not by equality. Hence, we have::

            sage: c == 15
            True
            sage: (a, b, 15) in L
            False

        TESTS::

            sage: a in L
            Traceback (most recent call last):
            ...
            KeyError: 'a'
            sage: (a, b) in L
            Traceback (most recent call last):
            ...
            KeyError: ('a', 'ab')
        """
        try:
            k1, k2, k3 = k
        except (TypeError, ValueError):
            raise KeyError(k)
        cdef triple_cell* cursor = self.lookup(<PyObject*>k1, <PyObject*>k2, <PyObject*>k3)
        if not valid(cursor.key_id1):
            return False
        if not self.weak_values:
            return True
        value = <object>cursor.value
        return not is_dead_keyedref(value)

    def __getitem__(self, k):
        """
        Get the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c]
            1
        """
        try:
            k1, k2, k3 = k
        except (TypeError, ValueError):
            raise KeyError(k)
        return self.get(k1, k2, k3)

    cdef get(self, k1, k2, k3):
        cdef triple_cell* cursor = self.lookup(<PyObject*>k1, <PyObject*>k2, <PyObject*>k3)
        if not valid(cursor.key_id1):
            raise KeyError((k1, k2, k3))
        value = <object>cursor.value
        if type(value) is KeyedRef:
            value = <object>PyWeakref_GET_OBJECT(value)
            if value is None:
                raise KeyError((k1, k2, k3))
        return value

    def __setitem__(self, k, value):
        """
        Set the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = -1
            sage: L[a,b,c]
            -1
        """
        try:
            k1, k2, k3 = k
        except (TypeError, ValueError):
            raise KeyError(k)
        self.set(k1, k2, k3, value)

    cdef int set(self, k1, k2, k3, value) except -1:
        cdef triple_cell entry
        cdef PyObject* old_value
        cdef bint maybe_resize = False
        entry.key_id1 = <PyObject*>k1
        entry.key_id2 = <PyObject*>k2
        entry.key_id3 = <PyObject*>k3
        if self.weak_values:
            wrap_k = (wrap(k1), wrap(k2), wrap(k3))
            try:
                value_store = KeyedRef(value, self.eraser, wrap_k)
                entry.value = <PyObject*>value_store
            except TypeError:
                entry.value = <PyObject*>value
        else:
            entry.value = <PyObject*>value
        Py_XINCREF(entry.value)
        cursor = self.lookup(<PyObject*>k1, <PyObject*>k2, <PyObject*>k3)
        if not valid(cursor.key_id1):
            self.used += 1
            if cursor.key_id1 is NULL:
                self.fill += 1
                maybe_resize = True
            if not self.weak_values:
                wrap_k = (wrap(k1), wrap(k2), wrap(k3))
            try:
                key_store = KeyedRef(k1, self.eraser, wrap_k)
                entry.key_weakref1 = <PyObject*>key_store
            except TypeError:
                entry.key_weakref1 = <PyObject*>k1
            Py_XINCREF(entry.key_weakref1)
            try:
                key_store = KeyedRef(k2, self.eraser, wrap_k)
                entry.key_weakref2 = <PyObject*>key_store
            except TypeError:
                entry.key_weakref2 = <PyObject*>k2
            Py_XINCREF(entry.key_weakref2)
            try:
                key_store = KeyedRef(k3, self.eraser, wrap_k)
                entry.key_weakref3 = <PyObject*>key_store
            except TypeError:
                entry.key_weakref3 = <PyObject*>k3
            Py_XINCREF(entry.key_weakref3)

            # We are taking a bit of a gamble here: we're assuming the
            # dictionary has not been resized (otherwise cursor might
            # not be a valid location anymore). The only way in which
            # that could happen is if the allocation activity above
            # forced a GC that triggered code that *adds* entries to
            # this dictionary: the dictionary can only get reshaped if
            # self.fill increases (as happens below). Note that we're
            # holding a strong ref to the dict itself, so that's not
            # liable to disappear. For the truly paranoid: we could
            # detect a change by checking if self.table has changed
            # value.
            cursor[0] = entry

            if maybe_resize and 3*self.fill > 2*self.mask:
                self.resize()
        else:
            old_value = cursor.value
            cursor.value = entry.value
            Py_XDECREF(old_value)

    def __delitem__(self, k):
        """
        Delete the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = -1
            sage: (a,b,c) in L
            True
            sage: del L[a,b,c]
            sage: len(L)
            0
            sage: (a,b,c) in L
            False
        """
        try:
            k1, k2, k3 = k
        except (TypeError, ValueError):
            raise KeyError(k)
        cdef triple_cell* cursor = self.lookup(<PyObject*>k1, <PyObject*>k2, <PyObject*>k3)
        if not valid(cursor.key_id1):
            raise KeyError(k)
        L = extract_triple_cell(cursor)
        self.used -= 1

    def items(self):
        """
        Iterate over the ``(key, value)`` pairs of this :class:`TripleDict`.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: L[1,2,3] = None
            sage: L.items()
            <...generator object at ...>
            sage: list(L.items())
            [((1, 2, 3), None)]
        """
        cdef size_t i = 0
        while i <= self.mask:
            cursor = &(self.table[i])
            i += 1
            if valid(cursor.key_id1):
                key1 = <object>(cursor.key_weakref1)
                key2 = <object>(cursor.key_weakref2)
                key3 = <object>(cursor.key_weakref3)
                value = <object>(cursor.value)
                if type(key1) is KeyedRef:
                    key1 = <object>PyWeakref_GET_OBJECT(key1)
                    if key1 is None:
                        print("found defunct key1")
                        continue
                if type(key2) is KeyedRef:
                    key2 = <object>PyWeakref_GET_OBJECT(key2)
                    if key2 is None:
                        print("found defunct key2")
                        continue
                if type(key3) is KeyedRef:
                    key3 = <object>PyWeakref_GET_OBJECT(key3)
                    if key3 is None:
                        print("found defunct key3")
                        continue
                if type(value) is KeyedRef:
                    value = <object>PyWeakref_GET_OBJECT(value)
                    if value is None:
                        print("found defunct value")
                        continue
                yield ((key1, key2, key3), value)

    def copy(self):
        """
        Return a copy of this :class:`TripleDict` as Python dict.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: L[1,2,3] = 42
            sage: L.copy()
            {(1, 2, 3): 42}
        """
        return dict(self.items())

    def __reduce__(self):
        """
        Note that we don't expect equality as this class concerns itself with
        object identity rather than object equality.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict()
            sage: L[1,2,3] = True
            sage: loads(dumps(L)) == L
            False
            sage: list(loads(dumps(L)).items())
            [((1, 2, 3), True)]
        """
        return TripleDict, (self.copy(),)

# The Cython supplied tp_traverse and tp_clear do not take the
# dynamically allocated table into account, so we have to supply our
# own. The only additional link to follow (that Cython does pick up
# and we have to replicate here) is the "eraser" which in its closure
# stores a reference back to the dictionary itself (meaning that
# TripleDicts only disappear on cyclic GC).
cdef int TripleDict_traverse(TripleDict self, visitproc visit, void* arg) noexcept:
    if not self.mask:
        return 0
    Py_VISIT3(<PyObject*>self.eraser, visit, arg)
    cdef size_t i
    for i in range(self.mask + 1):
        cursor = &self.table[i]
        if valid(cursor.key_id1):
            Py_VISIT3(cursor.key_weakref1, visit, arg)
            Py_VISIT3(cursor.key_weakref2, visit, arg)
            Py_VISIT3(cursor.key_weakref3, visit, arg)
            Py_VISIT3(cursor.value, visit, arg)


cdef int TripleDict_clear(TripleDict self) noexcept:
    if not self.mask:
        return 0
    cdef size_t mask = self.mask
    self.mask = 0  # Setting mask to 0 immediately prevents recursion
    self.used = 0
    self.fill = 0
    # Set self.eraser to None safely
    cdef object eraser = self.eraser
    self.eraser = None
    for i in range(mask + 1):
        cursor = &(self.table[i])
        if valid(cursor.key_id1):
            cursor.key_id1 = deleted_key
            Py_CLEAR(cursor.key_weakref1)
            Py_CLEAR(cursor.key_weakref2)
            Py_CLEAR(cursor.key_weakref3)
            Py_CLEAR(cursor.value)


(<PyTypeObject*>TripleDict).tp_traverse = <traverseproc>(&TripleDict_traverse)
(<PyTypeObject*>TripleDict).tp_clear = <inquiry>(&TripleDict_clear)
