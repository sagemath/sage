from cpython.ref cimport Py_INCREF

cdef _hook_tp_functions_type(PyTypeObject *t, newfunc tp_new, destructor tp_dealloc, bint useGC):
    """
    Initialize the fast integer creation functions.
    """
    cdef long flag

    # By default every object created in Pyrex is garbage
    # collected. This means it may have references to other objects
    # the Garbage collector has to look out for. We remove this flag
    # as the only reference an Integer has is to the global Integer
    # ring. As this object is unique we don't need to garbage collect
    # it as we always have a module level reference to it. If another
    # attribute is added to the Integer class this flag removal so as
    # the alloc and free functions may not be used anymore.
    # This object will still be reference counted.
    if not useGC:
        flag = Py_TPFLAGS_HAVE_GC
        t.tp_flags = <long>(t.tp_flags & (~flag))

    # Finally replace the functions called when an Integer needs
    # to be constructed/destructed.
    t.tp_new = tp_new
    t.tp_dealloc = tp_dealloc


cdef hook_tp_functions_type(object tp, newfunc tp_new, destructor tp_dealloc, bint useGC):
    cdef PyTypeObject *t = <PyTypeObject *>tp
    _hook_tp_functions_type(t, tp_new, tp_dealloc, useGC)


cdef hook_tp_functions(object global_dummy, newfunc tp_new, destructor tp_dealloc, bint useGC):
    """
    Initialize the fast integer creation functions.
    """
    # Make sure this never, ever gets collected.
    # This is not necessary for cdef'ed variables as the global
    # dummy integer, as such objects do not get automatically collected.
    # In fact there is no obvious reason to prevent collection when Sage quits
    # and we are certain no further call to the allocation function will be
    # made; so this could be removed when the code is clean enough.
    Py_INCREF(global_dummy)

    cdef PyTypeObject* t = Py_TYPE(global_dummy)
    _hook_tp_functions_type(t, tp_new, tp_dealloc, useGC)
