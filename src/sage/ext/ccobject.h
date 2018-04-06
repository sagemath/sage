/******************************************************************************
       Copyright (C) 2007 Joel B. Mohler <joel@kiwistrawberry.us>

   Distributed under the terms of the GNU General Public License (GPL)
   as published by the Free Software Foundation; either version 2 of
   the License, or (at your option) any later version.
                   http://www.gnu.org/licenses/

******************************************************************************/

/**
 * @file ccobject.h
 *
 * @author Joel B. Mohler <joel@kiwistrawberry.us>
 *
 * @brief These functions provide assistance for constructing and destructing
 *  C++ objects from Cython.
 *
 */

#ifndef __SAGE_CCOBJECT_H__
#define __SAGE_CCOBJECT_H__

#ifdef __cplusplus

#include <iostream>
#include <sstream>

/* Ok, we are in C++ mode.  Lets define some templated functions for construction
and destruction of allocated memory chunks. */

/* Allocate and Construct */
template <class T>
T* New(){
  return new T();
}

/* Construct */
template <class T>
T* Construct(void* mem){
  return new(mem) T();
}

/* Construct with one parameter */
template <class T, class P>
T* Construct_p(void* mem, const P &d){
  return new(mem) T(d);
}

/* Construct with two parameters */
template <class T, class P, class Q>
T* Construct_pp(void* mem, const P &d, const Q &e){
  return new(mem) T(d, e);
}

/* Construct with three parameters */
template <class T, class P, class Q, class R>
T* Construct_ppp(void* mem, const P &d, const Q &e, const R &f){
  return new(mem) T(d, e, f);
}

/* Destruct */
template <class T>
void Destruct(T* mem){
  mem->~T();
}

/* Deallocate and Destruct */
template <class T>
void Delete(T* mem){
  delete mem;
}


template <class T>
static CYTHON_INLINE int ccreadstr(T& x, PyObject* b)
{
    PyObject* converted = NULL;

#if PY_MAJOR_VERSION >= 3
    // Accept "str" input on Python 3
    if (PyUnicode_Check(b))
    {
        converted = PyUnicode_EncodeFSDefault(b);
        if (!converted) {return -1;}
        b = converted;
    }
#endif

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
#if PY_MAJOR_VERSION <= 2
    return PyString_FromStringAndSize(instr.c_str(), instr.size());
#else
    return PyUnicode_DecodeFSDefaultAndSize(instr.c_str(), instr.size());
#endif
}

/* Arrays */
template <class T>
static inline T* Allocate_array(size_t n){
  return new T[n];
}

template <class T>
static inline void Delete_array(T* v){
  delete[] v;
}

#endif

#endif /* ifndef __SAGE_CCOBJECT_H__ */
