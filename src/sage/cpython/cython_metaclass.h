/*****************************************************************************
*       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 2 of the License, or
* (at your option) any later version.
*                  http://www.gnu.org/licenses/
*****************************************************************************/

/* Compatibility for python 3.8, can be removed later */
#if PY_VERSION_HEX < 0x030900A4 && !defined(Py_SET_TYPE)
static inline void _Py_SET_TYPE(PyObject *ob, PyTypeObject *type)
{ ob->ob_type = type; }
#define Py_SET_TYPE(ob, type) _Py_SET_TYPE((PyObject*)(ob), type)
#endif

/* Tuple (None, None, None), initialized as needed */
static PyObject* NoneNoneNone;

/* All args flags of a PyMethod */
#define METH_ALLARGS (METH_VARARGS|METH_KEYWORDS|METH_NOARGS|METH_O)

/* Given an unbound method "desc" (this is not checked!) with only a
 * single "self" argument, call "desc(self)" without checking "self".
 * This can in particular be used to call any method as class or
 * static method. */
static CYTHON_INLINE PyObject* PyMethodDescr_CallSelf(PyMethodDescrObject* desc, PyObject* self)
{
    PyMethodDef* meth = desc->d_method;

    /* This must be a METH_NOARGS method */
    if (meth == NULL || (meth->ml_flags & METH_ALLARGS) != METH_NOARGS)
    {
        PyErr_SetString(PyExc_TypeError,
                "PyMethodDescr_CallSelf requires a method without arguments");
        return NULL;
    }

    return meth->ml_meth(self, NULL);
}

/*
 * This function calls PyType_Ready(t) and then calls
 * t.__getmetaclass__(None) (if that method exists) which should
 * return the metaclass for t. Then type(t) is set to this metaclass
 * and metaclass.__init__(t, None, None, None) is called.
 */
static CYTHON_INLINE int Sage_PyType_Ready(PyTypeObject* t)
{
    int r = PyType_Ready(t);
    if (r < 0)
        return r;

#if PY_VERSION_HEX >= 0x03050000
    // Cython 3 sets Py_TPFLAGS_HEAPTYPE before calling PyType_Ready,
    // and resets just after the call. We need to reset it earlier,
    // since otherwise the call to metaclass.__init__ below may have
    // illegal memory accesses.
    // See also:
    // https://github.com/cython/cython/issues/3603
    t->tp_flags &= ~Py_TPFLAGS_HEAPTYPE;
#endif

    /* Set or get metaclass (the type of t) */
    PyTypeObject* metaclass;

    PyObject* getmetaclass;
    getmetaclass = PyObject_GetAttrString((PyObject*)t, "__getmetaclass__");
    if (getmetaclass)
    {
        /* Call getmetaclass with self=None */
        metaclass = (PyTypeObject*)(PyMethodDescr_CallSelf((PyMethodDescrObject*)getmetaclass, Py_None));
        Py_DECREF(getmetaclass);
        if (!metaclass)
            return -1;

        if (!PyType_Check(metaclass))
        {
            PyErr_SetString(PyExc_TypeError,
                    "__getmetaclass__ did not return a type");
            return -1;
        }

        /* Now, set t.__class__ to metaclass */
        Py_SET_TYPE(t, metaclass);
        PyType_Modified(t);
    }
    else
    {
        /* No __getmetaclass__ method: read metaclass... */
        PyErr_Clear();
        metaclass = Py_TYPE(t);
    }

    /* Now call metaclass.__init__(t, None, None, None) unless
     * we would be calling type.__init__ */
    initproc init = metaclass->tp_init;
    if (init == NULL || init == PyType_Type.tp_init)
        return 0;

    /* Safety check: since we didn't call tp_new of metaclass,
     * we cannot safely call tp_init if the size of the structure
     * differs. */
    if (metaclass->tp_basicsize != PyType_Type.tp_basicsize)
    {
        PyErr_SetString(PyExc_TypeError,
                "metaclass is not compatible with 'type' (you cannot use cdef attributes in Cython metaclasses)");
        return -1;
    }

    /* Initialize a tuple (None, None, None) */
    if (!NoneNoneNone)
    {
        NoneNoneNone = PyTuple_Pack(3, Py_None, Py_None, Py_None);
        if (!NoneNoneNone) return -1;
    }
    return init((PyObject*)t, NoneNoneNone, NULL);
}


/* Use the above function in Cython code instead of the default
 * PyType_Ready() function */
#define PyType_Ready(t)  Sage_PyType_Ready(t)
