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
two parents (:trac:`715`) and to store homsets of pairs of objects of a
category (:trac:`11521`). In both cases, it is essential that the parent
structures remain garbage collectable, it is essential that the data access is
faster than with a usual :class:`~weakref.WeakKeyDictionary`, and we enforce
the "unique parent condition" in Sage (parent structures should be identical
if they are equal).

:class:`MonoDict` behaves similarly, but it takes a single item as a key. It
is used for caching the parents which allow a coercion map into a fixed other
parent (:trac:`12313`).

By :trac:`14159`, :class:`MonoDict` and :class:`TripleDict` can be optionally
used with weak references on the values.
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                     2012 Simon King <simon.king@uni-jena.de>
#                     2013 Nils Bruin <nbruin@sfu.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from cpython.object cimport *
from cpython.list cimport PyList_New
from cpython.mem cimport *
from cpython.weakref cimport PyWeakref_GetObject
from cpython.string cimport PyString_FromString
from cpython cimport Py_XINCREF, Py_XDECREF
from libc.string cimport memset
from weakref import KeyedRef, ref
#to avoid having to go and look in the module dictionary every time, we put
#a pointer to this class in here. No monkeypatching! or if you want to, perhaps
#cdef public this so that it can still be accessed.
#furthermore, by declaring this as "type", cython will compile
#isinstance(O,fixed_KeyedRef) to PyObject_CheckType(O, fixed_KeyedRef)
#which is considerable faster than PyObject_IsInstance(O, fixed_KeyedRef)
cdef type fixed_KeyedRef = KeyedRef
cdef type fixed_ref = ref

#we use a capsule to wrap our void* in the callback parameter.
#no assumptions on whether a void* fits in a Py_ssize_t anymore!
from cpython.pycapsule cimport PyCapsule_New, PyCapsule_GetPointer

cdef extern from "Python.h":
    PyObject* Py_None
    int PyList_SetItem(object list, Py_ssize_t index, PyObject * item) except -1

cdef extern from "pyx_visit.h":
    void Py_VISIT3(PyObject*, visitproc, void* arg)

#it's important that this is not an interned string: this object
#must be a unique sentinel. We could reuse the "dummy" sentinel
#that is defined in python's dictobject.c

cdef object dummy_object = PyString_FromString("dummy")
cdef PyObject* dummy = <PyObject*><void *>dummy_object

cdef struct mono_cell:
    void* key_id
    PyObject* key_weakref
    PyObject* value

cdef object extract_mono_cell(mono_cell *cell):
    #takes the refcounted components from a mono_cell
    #and puts them in a list and returns it.
    #The mono_cell itself is marked as "freed".
    #The refcounts originally accounting for the
    #presence in the mono_cell now account for the presence
    #in the returned list.
    #
    #the returned result is only used to throw away:
    #an advantage is that the containing list participates
    #in CPython's trashcan, which prevents stack overflow
    #on large dereffing cascades.
    #
    #a slight disadvantage is that this routine needs to
    #allocate a list (mainly just to be thrown away)
    if cell.key_id != NULL and cell.key_id != dummy :
        L=PyList_New(2)
        PyList_SetItem(L,0,cell.key_weakref)
        PyList_SetItem(L,1,cell.value)
        cell.key_id=dummy
        cell.key_weakref=NULL
        cell.value=NULL
        return L
    else:
        raise RuntimeError("unused mono_cell")

cdef struct triple_cell:
    void* key_id1
    void* key_id2
    void* key_id3
    PyObject* key_weakref1
    PyObject* key_weakref2
    PyObject* key_weakref3
    PyObject* value

cdef object extract_triple_cell(triple_cell *cell):
    #see extract_mono_cell for the rationale
    #behind this routine.
    if cell.key_id1 != NULL and cell.key_id1 != dummy :
        L=PyList_New(4)
        PyList_SetItem(L,0,cell.key_weakref1)
        PyList_SetItem(L,1,cell.key_weakref2)
        PyList_SetItem(L,2,cell.key_weakref3)
        PyList_SetItem(L,3,cell.value)
        cell.key_id1=dummy
        cell.key_id2=NULL
        cell.key_id3=NULL
        cell.key_weakref1=NULL
        cell.key_weakref2=NULL
        cell.key_weakref3=NULL
        cell.value=NULL
        return L
    else:
        raise RuntimeError("unused triple_cell")

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

        A :class:`MonoDict`.

        EXAMPLES::

            sage: k = set([1])
            sage: D = sage.structure.coerce_dict.MonoDict([(k,1)])
            sage: len(D)
            1
            sage: del k
            sage: len(D) # indirect doctest
            0
        """
        self.D = fixed_ref(D)

    def __call__(self, r):
        """
        INPUT:

        A weak reference with key.

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
        cdef MonoDict md = <object> PyWeakref_GetObject(self.D)
        if md is None:
            return
        if md.table == NULL:
            return
        cdef mono_cell* cursor = md.lookup(<PyObject *>PyCapsule_GetPointer(r.key,NULL))
        if (cursor.key_id != NULL and  cursor.key_id != dummy):
            if (cursor.key_weakref == <PyObject*><void*>r or cursor.value == <PyObject*><void*>r):
                L=extract_mono_cell(cursor)
                md.used -= 1
            else:
                #this should probably never happen
                raise RuntimeError("eraser: key match but no weakref match")

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

    If ki supports weak references then ri is a weak reference to ki with a
    callback to remove the entry from the dictionary if ki gets garbage
    collected. If ki is does not support weak references then ri is identical to ki.
    In the latter case the presence of the key in the dictionary prevents it from
    being garbage collected.

    INPUT:

    - ``size`` -- unused parameter, present for backward compatibility.
    - ``data`` -- optional iterable defining initial data.
    - ``threshold`` -- unused parameter, present for backward compatibility.
    - ``weak_values`` -- optional bool (default False). If it is true, weak references
      to the values in this dictionary will be used, when possible.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import MonoDict
        sage: L = MonoDict()
        sage: a = 'a'; b = 'ab'; c = '-15'
        sage: L[a] = 1
        sage: L[b] = 2
        sage: L[c] = 3

    The key is expected to be a unique object. Hence, the item stored for ``c``
    can not be obtained by providing another equal string::

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

        sage: # for some reason the following failed in "make ptest"
        sage: # on some installations, see #12313 for details
        sage: sorted(L.iteritems()) # random layout
        [('-15', 3), ('a', 1), ('ab', 2)]
        sage: # the following seems to be more consistent
        sage: set(L.iteritems())
        {('-15', 3), ('a', 1), ('ab', 2)}
        sage: del L[c]
        sage: sorted(L.iteritems())
        [('a', 1), ('ab', 2)]
        sage: len(L)
        2
        sage: for i in range(1000):
        ...       L[i] = i
        sage: len(L)
        1002
        sage: L['a']
        1
        sage: L['c']
        Traceback (most recent call last):
        ...
        KeyError: 'c'

    Note that this kind of dictionary is also used for caching actions
    and coerce maps. In previous versions of Sage, the cache was by
    strong references and resulted in a memory leak in the following
    example. However, this leak was fixed by :trac:`715`, using
    weak references::

        sage: K = GF(1<<55,'t')
        sage: for i in range(50):
        ...     a = K.random_element()
        ...     E = EllipticCurve(j=a)
        ...     P = E.random_point()
        ...     Q = 2*P
        sage: import gc
        sage: n = gc.collect()
        sage: from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
        sage: LE = [x for x in gc.get_objects() if isinstance(x, EllipticCurve_finite_field)]
        sage: len(LE)    # indirect doctest
        1

    TESTS:

    Here, we demonstrate the use of weak values.
    ::

        sage: M = MonoDict(13)
        sage: MW = MonoDict(13, weak_values=True)
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
    callbacks; compare :trac:`15069` (the mechanism used now is different, though)::

        sage: M = MonoDict(11)
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
        sage: del a
        Exception RuntimeError: 'maximum recursion depth exceeded while calling a Python object' in <function remove at ...> ignored
        sage: len(M)>0
        True

    Check that also in the presence of circular references, :class:`MonoDict`
    gets properly collected::

        sage: import gc
        sage: def count_type(T):
        ....:     return len([c for c in gc.get_objects() if isinstance(c,T)])
        sage: _=gc.collect()
        sage: N=count_type(MonoDict)
        sage: for i in range(100):
        ....:     V = [ MonoDict(11,{"id":j+100*i}) for j in range(100)]
        ....:     n= len(V)
        ....:     for i in range(n): V[i][V[(i+1)%n]]=(i+1)%n
        ....:     del V
        ....:     _=gc.collect()
        ....:     assert count_type(MonoDict) == N
        sage: count_type(MonoDict) == N
        True

    AUTHORS:

    - Simon King (2012-01)
    - Nils Bruin (2012-08)
    - Simon King (2013-02)
    - Nils Bruin (2013-11)
    """
    cdef mono_cell* lookup(self, PyObject* key):
        """
        This routine is used for all cases where a (potential) spot for
        a key is looked up. The returned value is a pointer into the dictionary
        store that either contains an entry with the requested key or a free spot
        where an entry for that key should go.
        """
        cdef size_t perturb
        cdef size_t mask = self.mask
        cdef mono_cell* table = self.table
        cdef mono_cell* cursor
        cdef mono_cell* first_freed = NULL

        #We seed our starting probe using the higher bits of the key as well.
        #Our hash is a memory location, so the bottom bits are likely 0.

        cdef size_t i = ((<size_t>key)>>8+(<size_t>key))
        if key == NULL or key == dummy:
            print("Request to look up invalid key")
        cursor = &(table[i & mask])
        # if the cell was never used, the entry wasn't there
        perturb = (<size_t>key)>>3

        #the probing algorithm is heavily inspired by python's own dict.
        #there is always at least one NULL entry in the store, and the probe
        #sequence eventually covers the entire store (see python's dictobject.c),
        #so the loop below does terminate. Since this loop executes only
        #straightforward C, we know the table will not change.

        while (cursor.key_id != key):
            if cursor.key_id == NULL:
                return first_freed if (first_freed != NULL) else cursor
            if first_freed == NULL and cursor.key_id == dummy:
                first_freed = cursor
            i = 5*i + perturb +1
            cursor = &(table[i & mask])
            perturb = perturb >> 5
        return cursor

    cdef int resize(self) except -1:
        """
        Resize dictionary. That can also mean shrink! Size is always a power of 2.
        """
        cdef mono_cell* old_table=self.table
        cdef size_t old_mask = self.mask
        cdef size_t newsize = 8
        cdef size_t minsize = 2*self.used
        cdef size_t i
        cdef mono_cell* cursor
        cdef mono_cell* entry
        while newsize < minsize:
            newsize = newsize<<1
        cdef mono_cell* table = <mono_cell*> PyMem_Malloc((newsize)*sizeof(mono_cell))
        if table == NULL:
            raise MemoryError()
        memset(table,0,(newsize)*sizeof(mono_cell))

        #we're done with memory activity. We can move the new (empty) table into place:

        self.table = table
        self.mask = newsize-1
        self.used = 0
        self.fill = 0

        #and move all entries over. We're not changing any refcounts here, so this is a very
        #tight loop that doesn't need to worry about tables changing.

        for i in range(old_mask+1):
            entry = &(old_table[i])
            if entry.key_id != NULL and entry.key_id != dummy:
                cursor=self.lookup(<PyObject*>(entry.key_id))
                assert cursor.key_id == NULL
                cursor[0]=entry[0]
                self.used +=1
                self.fill +=1
        PyMem_Free(old_table)
        return 0

    def __init__(self, size=11, data=None, threshold=0.7, weak_values=False):
        """
        Create a special dict using singletons for keys.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: a = 'a'
            sage: L[a] = 1
            sage: L[a]
            1
            sage: L = MonoDict({a: 1})
            sage: L[a]
            1
        """
        #if only one argument is supplied and it's iterable, use it for data rather than
        #for size. This way we're compatible with the old calling sequence (an integer is
        #just ignored) and we can also use the more usual construction.
        if data is None:
            try:
                data=size.iteritems()
            except AttributeError:
                try:
                    data=iter(size)
                except TypeError:
                    pass
        else:
            try:
                data=data.iteritems()
            except AttributeError:
                pass
        if self.table != NULL:
            raise RuntimeError("table already initialized. Called __init__ second time?")
        cdef minsize = 8
        cdef size_t newsize = 1<<3
        while newsize < minsize:
            newsize = newsize <<1
        self.mask = newsize - 1
        cdef mono_cell* table = <mono_cell*> PyMem_Malloc(newsize*sizeof(mono_cell))
        if table == NULL:
            raise MemoryError()
        memset(table,0,newsize*sizeof(mono_cell))
        self.table = table
        self.used = 0
        self.fill = 0
        self.eraser = MonoDictEraser(self)
        self.weak_values = weak_values
        if data:
            for k,v in data:
                self.set(k,v)

    def __dealloc__(self):
        """
        Ensure that the memory allocated by a MonoDict is properly freed.
        """
        #is this required or does this get called anyway if we set tp_clear properly?
        #at least it's safe, since MonoDict_clear checks if it has already run.
        MonoDict_clear(self)

    def __len__(self):
        """
        The number of items in self.
        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(37)
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
            sage: L = MonoDict(31)
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
        cdef mono_cell* cursor = self.lookup(<PyObject*><void*>k)
        if cursor.key_id == NULL or cursor.key_id == dummy:
            return False
        r = <object>cursor.key_weakref
        if isinstance(r, fixed_KeyedRef) and PyWeakref_GetObject(r) == Py_None:
            return False
        elif not(self.weak_values):
            return True
        else:
            value = <object>cursor.value
            return (not isinstance(value,fixed_KeyedRef)) or PyWeakref_GetObject(value) != Py_None

    def __getitem__(self, k):
        """
        Get the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: a = 'a'; b = 'b'; c = 15
            sage: L[a] = 1
            sage: L[b] = 2
            sage: L[c] = 3
            sage: L[c]                  # indirect doctest
            3

        Note that the keys are supposed to be unique::

            sage: c==15
            True
            sage: c is 15
            False
            sage: L[15]
            Traceback (most recent call last):
            ...
            KeyError: 15
        """
        return self.get(k)

    cdef get(self, object k):
        cdef mono_cell* cursor = self.lookup(<PyObject*><void *>k)
        if cursor.key_id == NULL or cursor.key_id == dummy:
            raise KeyError(k)
        r = <object>cursor.key_weakref
        if isinstance(r, fixed_KeyedRef) and PyWeakref_GetObject(r) == Py_None:
            raise KeyError(k)
        value = <object>cursor.value
        if self.weak_values and isinstance(value,fixed_KeyedRef):
            value = <object>PyWeakref_GetObject(value)
            if value is None:
                raise KeyError(k)
        return value

    def __setitem__(self, k, value):
        """
        Set the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
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
        self.set(k,value)

    cdef set(self,object k, object value):
        cdef mono_cell entry
        cdef PyObject* old_value = NULL
        cdef bint maybe_resize = False
        entry.key_id = <void*>k
        if self.weak_values:
            cap_k=PyCapsule_New(<void *>(k),NULL,NULL)
            try:
                value_store = fixed_KeyedRef(value,self.eraser,cap_k)
                entry.value = <PyObject*><void*>value_store
            except TypeError:
                entry.value = <PyObject*><void*>value
        else:
            entry.value = <PyObject*><void*>value
        Py_XINCREF(entry.value)
        cursor = self.lookup(<PyObject*><void*>k)
        if cursor.key_id == NULL or cursor.key_id == dummy:
            self.used += 1
            if cursor.key_id == NULL:
                self.fill += 1
                maybe_resize = True
            if not(self.weak_values):
                cap_k=PyCapsule_New(<void *>(k),NULL,NULL)
            try:
                key_store=fixed_KeyedRef(k,self.eraser,cap_k)
                entry.key_weakref = <PyObject*><void*>key_store
            except TypeError:
                entry.key_weakref = <PyObject*><void*>k
            Py_XINCREF(entry.key_weakref)

            #we're taking a bit of a gamble here: we're assuming the dictionary has
            #not been resized (otherwise cursor might not be a valid location
            #anymore). The only way in which that could happen is if the allocation
            #activity above forced a GC that triggered code that ADDS entries to this
            #dictionary: the dictionary can only get reshaped if self.fill increases.
            #(as happens below). Note that we're holding a strong ref to the dict
            #itself, so that's not liable to disappear.
            #for the truly paranoid: we could detect a change by checking if self.table
            #has changed value
            cursor[0] = entry

            #this is the one place where resize gets called:
            if maybe_resize and 3*self.fill > 2*self.mask: self.resize()
        else:
            old_value = cursor.value
            cursor.value = entry.value
            Py_XDECREF(old_value)

    def __delitem__(self, k):
        """
        Delete the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
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
        cdef mono_cell* cursor = self.lookup(<PyObject *><void *>k)
        if cursor.key_id == NULL or cursor.key_id==dummy:
            raise KeyError(k)
        L=extract_mono_cell(cursor)
        self.used -= 1

    def iteritems(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: L[1] = None
            sage: L[2] = True
            sage: list(sorted(L.iteritems()))
            [(1, None), (2, True)]
        """
        #iteration is tricky because the table could change from under us.
        #the following iterates properly if the dictionary does not
        #get "resize"d, which is guaranteed if no NEW entries in the
        #dictionary are introduced. At least we make sure to get our data fresh
        #from "self" every iteration, so that at least we're not reading random memory
        #(if the dictionary changes, it's not guaranteed you get to see any particular entry)
        cdef size_t i = 0
        while i <= self.mask:
            cursor = &(self.table[i])
            i += 1
            if cursor.key_id != NULL and cursor.key_id != dummy:
                key = <object>(cursor.key_weakref)
                value = <object>(cursor.value)
                if isinstance(key,fixed_KeyedRef):
                    key = <object>PyWeakref_GetObject(key)
                    if key is None:
                        print("found defunct key")
                        continue
                if self.weak_values and isinstance(value,fixed_KeyedRef):
                    value = <object>PyWeakref_GetObject(value)
                    if value is None:
                        print("found defunct value")
                        continue
                yield (key, value)

    def __reduce__(self):
        """
        Note that we don't expect equality as this class concerns itself with
        object identity rather than object equality.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: L[1] = True
            sage: loads(dumps(L)) == L
            False
            sage: list(loads(dumps(L)).iteritems())
            [(1, True)]
        """
        return MonoDict, (11, dict(self.iteritems()), 0.7)

#the cython supplied tp_traverse and tp_clear do not take the dynamically
#allocated table into account, so we have to supply our own.
#the only additional link to follow (that cython does pick up and we have
#to replicate here) is the "eraser" which in its closure stores a reference
#back to the dictionary itself (meaning that MonoDicts only disappear
#on cyclic GC)

cdef int MonoDict_traverse(MonoDict op, visitproc visit, void *arg):
    if op.table == NULL:
        return 0
    Py_VISIT3(<PyObject*>op.eraser, visit, arg)
    cdef size_t i
    for i in range(op.mask + 1):
        cursor = &op.table[i]
        if cursor.key_id != NULL and cursor.key_id != dummy:
            Py_VISIT3(cursor.key_weakref, visit, arg)
            Py_VISIT3(cursor.value, visit, arg)
    return 0


#we clear a monodict by taking first taking away the table before dereffing
#its contents. That shortcuts callbacks, so we deref the entries straight here.
#that means this code does not participate in Python's trashcan the way that
#deletion code based on extract_mono_cell does, so there is probably a way
#this code can be used to overflow the C stack. It would have to be a pretty
#devious example, though.
cdef int MonoDict_clear(MonoDict op):
    if op.table == NULL:
        return 0

    tmp=op.eraser
    cdef mono_cell* table=op.table
    cdef size_t mask=op.mask
    op.table=NULL

    op.eraser=None
    op.mask=0
    op.used=0
    op.fill=0

    #this deletion method incurs an extra refcount +-, but it seems very difficult in cython to do it
    #any other way and still be sure that the reference op.eraser is gone by the time the object gets
    #deallocated.
    del tmp
    for i in range(mask+1):
        cursor = &(table[i])
        if cursor.key_id != NULL and cursor.key_id != dummy:
            cursor.key_id = dummy
            Py_XDECREF(cursor.key_weakref)
            Py_XDECREF(cursor.value)
    PyMem_Free(table)
    return 0

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

        A :class:`TripleDict`. For internal use only.

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
        self.D = fixed_ref(D)

    def __call__(self, r):
        """
        INPUT:

        A weak reference with key.

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
        cdef TripleDict td = <object>PyWeakref_GetObject(self.D)
        if td is None:
            return
        if td.table == NULL:
            return

        k1,k2,k3 = r.key
        cdef triple_cell* cursor = td.lookup(<PyObject*>PyCapsule_GetPointer(k1,NULL),
                                             <PyObject*>PyCapsule_GetPointer(k2,NULL),
                                             <PyObject*>PyCapsule_GetPointer(k3,NULL))
        if (cursor.key_id1 != NULL and  cursor.key_id1 != dummy):
            if (cursor.key_weakref1 == <PyObject*><void*>r or
                    cursor.key_weakref2 == <PyObject*><void*>r or
                    cursor.key_weakref3 == <PyObject*><void*>r or
                    cursor.value == <PyObject*><void*>r):
                L=extract_triple_cell(cursor)
                td.used -= 1
            else:
                raise RuntimeError("eraser: key match but no weakref match")

cdef class TripleDict:
    """
    This is a hashtable specifically designed for (read) speed in
    the coercion model.

    It differs from a python dict in the following important ways:

       - All keys must be sequence of exactly three elements. All sequence
         types (tuple, list, etc.) map to the same item.
       - Comparison is done using the 'is' rather than '==' operator.

    There are special cdef set/get methods for faster access.
    It is bare-bones in the sense that not all dictionary methods are
    implemented.

    It is implemented as a list of lists (hereafter called buckets). The bucket is
    chosen according to a very simple hash based on the object pointer, and each
    bucket is of the form [id(k1), id(k2), id(k3), r1, r2, r3, value, id(k1),
    id(k2), id(k3), r1, r2, r3, value, ...], on which a linear search is performed.
    If a key component ki supports weak references then ri is a weak reference to
    ki; otherwise ri is identical to ki.

    INPUT:

    - ``size`` -- an integer, the initial number of buckets. To spread objects
      evenly, the size should ideally be a prime, and certainly not divisible
      by 2.
    - ``data`` -- optional iterable defining initial data.
    - ``threshold`` -- optional number, default `0.7`. It determines how frequently
      the dictionary will be resized (large threshold implies rare resizing).
    - ``weak_values`` -- optional bool (default False). If it is true, weak references
      to the values in this dictionary will be used, when possible.

    If any of the key components k1,k2,k3 (this can happen for a key component
    that supports weak references) gets garbage collected then the entire
    entry disappears. In that sense this structure behaves like a nested
    :class:`~weakref.WeakKeyDictionary`.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import TripleDict
        sage: L = TripleDict()
        sage: a = 'a'; b = 'b'; c = 'c'
        sage: L[a,b,c] = 1
        sage: L[a,b,c]
        1
        sage: L[c,b,a] = -1
        sage: list(L.iteritems())     # random order of output.
        [(('c', 'b', 'a'), -1), (('a', 'b', 'c'), 1)]
        sage: del L[a,b,c]
        sage: list(L.iteritems())
        [(('c', 'b', 'a'), -1)]
        sage: len(L)
        1
        sage: for i in range(1000):
        ...       L[i,i,i] = i
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

    Note that this kind of dictionary is also used for caching actions
    and coerce maps. In previous versions of Sage, the cache was by
    strong references and resulted in a memory leak in the following
    example. However, this leak was fixed by :trac:`715`, using weak
    references::

        sage: K = GF(1<<55,'t')
        sage: for i in range(50):
        ...     a = K.random_element()
        ...     E = EllipticCurve(j=a)
        ...     P = E.random_point()
        ...     Q = 2*P
        sage: import gc
        sage: n = gc.collect()
        sage: from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
        sage: LE = [x for x in gc.get_objects() if isinstance(x, EllipticCurve_finite_field)]
        sage: len(LE)    # indirect doctest
        1

    TESTS:

    Here, we demonstrate the use of weak values.
    ::

        sage: class Foo: pass
        sage: T = TripleDict(13)
        sage: TW = TripleDict(13, weak_values=True)
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

    .. NOTE::

        The index `h` corresponding to the key [k1, k2, k3] is computed as a
        value of unsigned type size_t as follows:

        .. MATH::

            h = id(k1) + 13*id(k2) xor 503 id(k3)

        The natural type for this quantity is Py_ssize_t, which is a signed
        quantity with the same length as size_t. Storing it in a signed way gives the most
        efficient storage into PyInt, while preserving sign information.

        In previous situations there were some problems with ending up with negative
        indices, which required casting to an unsigned type, i.e.,
        (<size_t> h)% N
        since C has a sign-preserving % operation This caused problems on 32 bits systems,
        see :trac:`715` for details. This is irrelevant for the current implementation.

    AUTHORS:

    - Robert Bradshaw, 2007-08

    - Simon King, 2012-01

    - Nils Bruin, 2012-08

    - Simon King, 2013-02

    - Nils Bruin, 2013-11
    """
    cdef triple_cell* lookup(self, PyObject* key1, PyObject* key2, PyObject* key3):
        #returns a pointer to where key should be stored in this dictionary.
        cdef size_t perturb
        cdef size_t mask = self.mask
        cdef triple_cell* table = self.table
        cdef triple_cell* cursor
        cdef triple_cell* first_freed = NULL
        cdef int j
        global summed_expected, samples
        #we device some hash that reasonably depends on bits of all keys
        #making sure it's not commutative and also involves some of the higher order
        #bits at an early stage.
        cdef size_t key = <size_t>key1 + 13*(<size_t>key2) ^ 503*(<size_t>key3)
        cdef size_t i = key>>8 + key
        cursor = &(table[i & mask])
        perturb = (key)>>3
        j=1
        while (cursor.key_id1 != key1 or cursor.key_id2 != key2 or cursor.key_id3 != key3):
            if cursor.key_id1 == NULL:
                return first_freed if (first_freed != NULL) else cursor
            if first_freed == NULL and cursor.key_id1 == dummy:
                first_freed = cursor
            i = 5*i + perturb +1
            cursor = &(table[i & mask])
            perturb = perturb >> 5
        return cursor

    cdef int resize(self) except -1:
        cdef triple_cell* old_table=self.table
        cdef size_t old_mask = self.mask
        cdef size_t newsize = 8
        cdef size_t minsize = 2*self.used
        cdef size_t i
        cdef triple_cell* cursor
        cdef triple_cell* entry
        while newsize < minsize:
            newsize = newsize<<1
        cdef triple_cell* table = <triple_cell*> PyMem_Malloc((newsize)*sizeof(triple_cell))
        if table == NULL:
            raise MemoryError()
        memset(table,0,(newsize)*sizeof(triple_cell))
        self.table = table
        self.mask = newsize-1
        self.used = 0
        self.fill = 0
        for i in range(old_mask+1):
            entry = &(old_table[i])
            if entry.key_id1 != NULL and entry.key_id1 != dummy:
                cursor=self.lookup(<PyObject*>(entry.key_id1),<PyObject*>(entry.key_id2),<PyObject*>(entry.key_id3))
                assert cursor.key_id1 == NULL
                cursor[0]=entry[0]
                self.used +=1
                self.fill +=1
        PyMem_Free(old_table)
        return 0

    def __init__(self, size=11, data=None, threshold=0.7, weak_values=False):
        """
        Create a special dict using triples for keys.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c]
            1
        """
        #if only one argument is supplied and it's iterable, use it for data rather than
        #for size
        if data is None:
            try:
                data=size.iteritems()
            except AttributeError:
                try:
                    data=iter(size)
                except TypeError:
                    pass
        else:
            try:
                data=data.iteritems()
            except AttributeError:
                pass
        if self.table != NULL:
            raise RuntimeError("table already initialized. Called __init__ second time?")
        cdef minsize = 8
        cdef size_t newsize = 1<<3
        while newsize < minsize:
            newsize = newsize <<1
        self.mask = newsize - 1
        cdef triple_cell* table = <triple_cell*> PyMem_Malloc(newsize*sizeof(triple_cell))
        if table == NULL:
            raise MemoryError()
        memset(table,0,newsize*sizeof(triple_cell))
        self.table = table
        self.used = 0
        self.fill = 0
        self.eraser = TripleDictEraser(self)
        self.weak_values = weak_values
        if data:
            for k,v in data:
                self[k]=v

    def __dealloc__(self):
        #is this required? (TripleDict_clear is safe to call multiple times)
        TripleDict_clear(self)

    def __len__(self):
        """
        The number of items in self.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(37)
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
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'ab'; c = 15
            sage: L[a,b,c] = 123
            sage: (a,b,c) in L         # indirect doctest
            True

        The keys are compared by identity, not by equality. Hence, we have::

            sage: c == 15
            True
            sage: (a,b,15) in L
            False
        """
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError(k)
        cdef triple_cell* cursor = self.lookup(<PyObject*><void*>k1,<PyObject*><void*>k2,<PyObject*><void*>k3)
        if cursor.key_id1 == NULL or cursor.key_id1 == dummy:
            return False
        r = <object>cursor.key_weakref1
        if isinstance(r, fixed_KeyedRef) and PyWeakref_GetObject(r) == Py_None:
            return False
        r = <object>cursor.key_weakref2
        if isinstance(r, fixed_KeyedRef) and PyWeakref_GetObject(r) == Py_None:
            return False
        r = <object>cursor.key_weakref3
        if isinstance(r, fixed_KeyedRef) and PyWeakref_GetObject(r) == Py_None:
            return False
        if not(self.weak_values):
            return True
        value = <object>cursor.value
        return (not isinstance(value,fixed_KeyedRef)) or PyWeakref_GetObject(value) != Py_None

    def __getitem__(self, k):
        """
        Get the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c]
            1
        """
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError(k)
        return self.get(k1, k2, k3)

    cdef get(self, object k1, object k2, object k3):
        cdef triple_cell* cursor = self.lookup(<PyObject*><void *>k1,<PyObject*><void *>k2,<PyObject*><void *>k3)
        if cursor.key_id1 == NULL or cursor.key_id1 == dummy:
            raise KeyError((k1, k2, k3))
        r1 = <object>cursor.key_weakref1
        r2 = <object>cursor.key_weakref2
        r3 = <object>cursor.key_weakref3
        if (isinstance(r1, fixed_KeyedRef) and PyWeakref_GetObject(r1) == Py_None) or \
                (isinstance(r2, fixed_KeyedRef) and PyWeakref_GetObject(r2) == Py_None) or \
                (isinstance(r3, fixed_KeyedRef) and PyWeakref_GetObject(r3) == Py_None):
            raise KeyError((k1, k2, k3))
        value = <object>cursor.value
        if self.weak_values and isinstance(value,fixed_KeyedRef):
            value = <object>PyWeakref_GetObject(value)
            if value is None:
                raise KeyError((k1, k2, k3))
        return value

    def __setitem__(self, k, value):
        """
        Set the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = -1
            sage: L[a,b,c]
            -1
        """
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError(k)
        self.set(k1, k2, k3, value)

    cdef set(self, object k1, object k2, object k3, value):
        cdef triple_cell entry
        cdef PyObject* old_value = NULL
        cdef bint maybe_resize = False
        entry.key_id1 = <void*>k1
        entry.key_id2 = <void*>k2
        entry.key_id3 = <void*>k3
        if self.weak_values:
            k = (PyCapsule_New(<void *>(k1),NULL,NULL),
                 PyCapsule_New(<void *>(k2),NULL,NULL),
                 PyCapsule_New(<void *>(k3),NULL,NULL))
            try:
                value_store = fixed_KeyedRef(value,self.eraser,k)
                entry.value = <PyObject*><void*>value_store
            except TypeError:
                entry.value = <PyObject*><void*>value
        else:
            entry.value = <PyObject*><void*>value
        Py_XINCREF(entry.value)
        cursor = self.lookup(<PyObject*><void*>k1,<PyObject*><void*>k2,<PyObject*><void*>k3)
        if cursor.key_id1 == NULL or cursor.key_id1 == dummy:
            self.used += 1
            if cursor.key_id1 == NULL:
                self.fill += 1
                maybe_resize = True
            if not(self.weak_values):
                k = (PyCapsule_New(<void *>(k1),NULL,NULL),
                     PyCapsule_New(<void *>(k2),NULL,NULL),
                     PyCapsule_New(<void *>(k3),NULL,NULL))
            try:
                key_store=fixed_KeyedRef(k1,self.eraser,k)
                entry.key_weakref1 = <PyObject*><void*>key_store
            except TypeError:
                entry.key_weakref1 = <PyObject*><void*>k1
            Py_XINCREF(entry.key_weakref1)
            try:
                key_store=fixed_KeyedRef(k2,self.eraser,k)
                entry.key_weakref2 = <PyObject*><void*>key_store
            except TypeError:
                entry.key_weakref2 = <PyObject*><void*>k2
            Py_XINCREF(entry.key_weakref2)
            try:
                key_store=fixed_KeyedRef(k3,self.eraser,k)
                entry.key_weakref3 = <PyObject*><void*>key_store
            except TypeError:
                entry.key_weakref3 = <PyObject*><void*>k3
            Py_XINCREF(entry.key_weakref3)

            #we're taking a bit of a gamble here: we're assuming the dictionary has
            #not been resized (otherwise cursor might not be a valid location
            #anymore). The only way in which that could happen is if the allocation
            #activity above forced a GC that triggered code that ADDS entries to this
            #dictionary: the dictionary can only get reshaped if self.fill increases.
            #(as happens below). Note that we're holding a strong ref to the dict
            #itself, so that's not liable to disappear.
            #for the truly paranoid: we could detect a change by checking if self.table
            #has changed value
            cursor[0] = entry

             #this is the only place where resize gets called:
            if maybe_resize and 3*self.fill > 2*self.mask: self.resize()
        else:
            old_value = cursor.value
            cursor.value = entry.value
            Py_XDECREF(old_value)

    def __delitem__(self, k):
        """
        Delete the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
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
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError(k)
        cdef triple_cell* cursor = self.lookup(<PyObject *><void *>k1,<PyObject *><void *>k2,<PyObject *><void *>k3)
        if cursor.key_id1 == NULL or cursor.key_id1==dummy:
            raise KeyError(k)
        L=extract_triple_cell(cursor)
        self.used -= 1

    def iteritems(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: L[1,2,3] = None
            sage: list(L.iteritems())
            [((1, 2, 3), None)]
        """

        cdef size_t i = 0
        while i <= self.mask:
            cursor = &(self.table[i])
            i += 1
            if cursor.key_id1 != NULL and cursor.key_id1 != dummy:
                key1 = <object>(cursor.key_weakref1)
                key2 = <object>(cursor.key_weakref2)
                key3 = <object>(cursor.key_weakref3)
                value = <object>(cursor.value)
                if isinstance(key1,fixed_KeyedRef):
                    key1 = <object>PyWeakref_GetObject(key1)
                    if key1 is None:
                        print("found defunct key1")
                        continue
                if isinstance(key2,fixed_KeyedRef):
                    key2 = <object>PyWeakref_GetObject(key2)
                    if key2 is None:
                        print("found defunct key2")
                        continue
                if isinstance(key3,fixed_KeyedRef):
                    key3 = <object>PyWeakref_GetObject(key3)
                    if key3 is None:
                        print("found defunct key3")
                        continue
                if self.weak_values and isinstance(value,fixed_KeyedRef):
                    value = <object>PyWeakref_GetObject(value)
                    if value is None:
                        print("found defunct value")
                        continue
                yield ((key1, key2, key3), value)

    def __reduce__(self):
        """
        Note that we don't expect equality as this class concerns itself with
        object identity rather than object equality.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: L[1,2,3] = True
            sage: loads(dumps(L)) == L
            False
            sage: list(loads(dumps(L)).iteritems())
            [((1, 2, 3), True)]
        """
        return TripleDict, (11, dict(self.iteritems()), 0.7)

#the cython supplied tp_traverse and tp_clear do not take the dynamically
#allocated table into account, so we have to supply our own.
#the only additional link to follow (that cython does pick up and we have
#to replicate here) is the "eraser" which in its closure stores a reference
#back to the dictionary itself (meaning that MonoDicts only disappear
#on cyclic GC)

cdef int TripleDict_traverse(TripleDict op, visitproc visit, void *arg):
    if op.table == NULL:
        return 0
    Py_VISIT3(<PyObject*>op.eraser, visit, arg)
    cdef size_t i
    for i in range(op.mask + 1):
        cursor = &op.table[i]
        if cursor.key_id1 != NULL and cursor.key_id1 != dummy:
            Py_VISIT3(cursor.key_weakref1, visit, arg)
            Py_VISIT3(cursor.key_weakref2, visit, arg)
            Py_VISIT3(cursor.key_weakref3, visit, arg)
            Py_VISIT3(cursor.value, visit, arg)
    return 0

cdef int TripleDict_clear(TripleDict op):
    if op.table == NULL:
        return 0

    tmp=op.eraser
    cdef triple_cell* table=op.table
    cdef size_t mask=op.mask
    op.table=NULL
    op.eraser=None
    op.mask=0
    op.used=0
    op.fill=0

    #refcount dance to ensure op.eraser==None when the actual object gets deallocated
    del tmp
    for i in range(mask+1):
        cursor = &(table[i])
        if cursor.key_id1 != NULL and cursor.key_id1 != dummy:
            cursor.key_id1 = dummy
            Py_XDECREF(cursor.key_weakref1)
            Py_XDECREF(cursor.key_weakref2)
            Py_XDECREF(cursor.key_weakref3)
            Py_XDECREF(cursor.value)
    PyMem_Free(table)
    return 0

(<PyTypeObject*>TripleDict).tp_traverse = <traverseproc>(&TripleDict_traverse)
(<PyTypeObject*>TripleDict).tp_clear = <inquiry>(&TripleDict_clear)
