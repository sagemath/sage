# sage_setup: distribution = sagemath-objects
"""
Delete item from PyDict by exact value and hash

Beware that the implementation of the routine here relies on implementation
details of CPython's dict that go beyond the published API.

AUTHORS:

- Nils Bruin (2017-05)
"""

# ****************************************************************************
#       Copyright (C) 2017 Nils Bruin <nbruin@sfu.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.list cimport PyList_New

cdef extern from "Python.h":
    ctypedef struct PyDictKeysObject

    ctypedef struct PyDictObject:
        Py_ssize_t ma_used
        PyDictKeysObject * ma_keys
        PyObject ** ma_values

    int PyList_SetItem(object list, Py_ssize_t index, PyObject * item) except -1

cdef extern from "dict_internal.h":
    Py_ssize_t DK_MASK(PyDictKeysObject *)
    PyDictKeyEntry * DK_ENTRIES(PyDictKeysObject *keys)

    Py_ssize_t dictkeys_get_index (PyDictKeysObject *keys, Py_ssize_t i)
    void dictkeys_set_index (PyDictKeysObject *keys, Py_ssize_t i, Py_ssize_t ix)

    Py_ssize_t DKIX_EMPTY, DKIX_DUMMY
    int PERTURB_SHIFT

    ctypedef struct PyDictKeyEntry:
        Py_hash_t me_hash
        PyObject * me_key
        PyObject * me_value

cdef int del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, Py_hash_t hash) except -1:
    """
    This is used in callbacks for the weak values of :class:`WeakValueDictionary`.

    INPUT:

    - ``PyDictObject *mp`` -- pointer to a dict
    - ``PyObject *value`` -- pointer to a value of the dictionary
    - ``Py_hash_t hash`` -- hash of the key by which the value is stored in the dict

    The hash bucket determined by the given hash is searched for the item
    containing the given value. If this item cannot be found, the function is
    silently returning. Otherwise, the item is removed from the dict.

    TESTS:

    The following is an indirect doctest, as discussed on :issue:`13394`.
    ::

        sage: from sage.misc.weak_dict import WeakValueDictionary
        sage: V = [set(range(n)) for n in range(5)]
        sage: D = WeakValueDictionary(enumerate(V))

    The line ``V[k] = None`` triggers execution of the callback functions of
    the dict values. However, the actual deletion is postponed till after the
    iteration over the dictionary has finished. Hence, when the callbacks are
    executed, the values which the callback belongs to has already been
    overridden by a new value. Therefore, the callback does not delete the
    item::

        sage: for k in D:    # indirect doctest
        ....:     V[k] = None
        ....:     D[k] = ZZ
        sage: len(D)
        5
        sage: D[1]
        Integer Ring

    TESTS:

    The following shows that the deletion of deeply nested structures does not
    result in an error, by :issue:`15506`::

        sage: class A: pass
        sage: a = A(); prev = a
        sage: M = WeakValueDictionary()
        sage: for i in range(10^3+10): newA = A(); M[newA] = prev; prev = newA
        sage: del a
    """
    keys = mp.ma_keys
    cdef size_t perturb
    cdef size_t mask = DK_MASK(keys)
    cdef PyDictKeyEntry *entries = DK_ENTRIES(keys)
    cdef PyDictKeyEntry *ep

    if mp.ma_values is not NULL:
        raise TypeError("del_dictitem_by_exact_value cannot be applied to a shared key dict")

    cdef size_t i = <size_t>hash & mask
    ix = dictkeys_get_index(keys, i)

    if ix == DKIX_EMPTY:
        # key not found
        return 0

    ep = &(entries[ix])
    perturb = hash
    while (ep.me_value != value or ep.me_hash != hash):
        perturb = perturb >> PERTURB_SHIFT
        i = mask & (i * 5 + perturb + 1)
        ix = dictkeys_get_index(keys, i)
        if ix == DKIX_EMPTY:
            # key not found
            return 0
        ep = &(entries[ix])

    T = PyList_New(2)
    PyList_SetItem(T, 0, ep.me_key)
    PyList_SetItem(T, 1, ep.me_value)
    ep.me_key = NULL
    ep.me_value = NULL
    mp.ma_used -= 1
    dictkeys_set_index(keys, i, DKIX_DUMMY)
    # We have transferred the to-be-deleted references to the list T
    # we now delete the list so that the actual decref happens through a
    # deallocation routine that uses the Python Trashcan macros to
    # avoid stack overflow in deleting deep structures.
    del T


def test_del_dictitem_by_exact_value(D, value, h):
    """
    This function helps testing some cdef function used to delete dictionary items.

    INPUT:

    - ``D`` -- a Python ``<dict>``
    - ``value`` -- an object that is value ``D``
    - ``h`` -- the hash of the key under which to find ``value`` in ``D``

    The underlying cdef function deletes an item from ``D`` that is in the
    hash bucket determined by ``h`` and whose value is identic with
    ``value``. Of course, this only makes sense if the pairs ``(h, value)``
    corresponding to items in ``D`` are pair-wise distinct.

    If a matching item cannot be found, the function does nothing and
    silently returns.

    TESTS:

    See :issue:`13394` for a discussion.
    ::

        sage: from sage.cpython.dict_del_by_value import test_del_dictitem_by_exact_value
        sage: B=1000
        sage: L=list(range(B))
        sage: D1=dict()
        sage: D2=dict()
        sage: for i in range(100000):        # long time
        ....:     ki=L[floor(random()*B)]
        ....:     vi=L[floor(random()*B)]
        ....:     D1[ki]=vi
        ....:     D2[ki]=vi
        ....:     ko=L[floor(random()*B)]
        ....:     if ko in D1:
        ....:         vo=D1[ko]
        ....:         del D1[ko]
        ....:         test_del_dictitem_by_exact_value(D2,vo,hash(ko))
        ....:     assert D1 == D2

    No action is taken if the item prescribed by key hash and value does not
    exist in the dictionary::

        sage: D = {1: ZZ}
        sage: test_del_dictitem_by_exact_value(D, ZZ, 2)
        sage: D
        {1: Integer Ring}
        sage: test_del_dictitem_by_exact_value(D, QQ, 1)
        sage: D
        {1: Integer Ring}
    """
    del_dictitem_by_exact_value(<PyDictObject *>D, <PyObject *>value, h)
