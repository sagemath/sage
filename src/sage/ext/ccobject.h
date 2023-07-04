#ifndef __SAGE_CCOBJECT_H__
#define __SAGE_CCOBJECT_H__

#ifdef __cplusplus

#include <iostream>
#include <sstream>


template <class T>
static CYTHON_INLINE int ccreadstr(T& x, PyObject* b)
{
    PyObject* converted = NULL;

    // Accept "str" input
    if (PyUnicode_Check(b))
    {
        converted = PyUnicode_EncodeFSDefault(b);
        if (!converted) {return -1;}
        b = converted;
    }

    char* buffer;
    Py_ssize_t length;

    if (PyBytes_AsStringAndSize(b, &buffer, &length) == -1)
        {Py_XDECREF(converted); return -1;}
    std::istringstream input(std::string(buffer, length));
    Py_XDECREF(converted);

    input >> x;

    return 0;
}


template <class T>
static CYTHON_INLINE PyObject* ccrepr(const T& x)
{
    std::ostringstream instore;
    instore << x;
    std::string instr = instore.str();
    return PyUnicode_DecodeFSDefaultAndSize(instr.c_str(), instr.size());
}


/* Arrays */
template <class T>
static inline T* Allocate_array(size_t n)
{
    return new T[n];
}

template <class T>
static inline void Delete_array(T* v)
{
    delete[] v;
}

#endif

#endif /* ifndef __SAGE_CCOBJECT_H__ */
