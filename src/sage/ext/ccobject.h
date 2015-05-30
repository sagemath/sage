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

/* Construct with four parameters */
template <class T, class P, class Q, class R, class S>
T* Construct_pppp(void* mem, const P &d, const Q &e, const R &f, const S &g){
  return new(mem) T(d, e, f, g);
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
bool _equal(T lhs, T rhs)
{
    return lhs == rhs;
}

template <class T>
void _from_str(T* dest, const char* src){
  std::istringstream out(src);
  out >> *dest;
}

template <class T>
void _from_str_len(T* dest, const char* src, unsigned int len){
  std::istringstream out(std::string(src, len));
  out >> *dest;
}

template <class T>
PyObject* _to_PyString(const T *x)
{
  std::ostringstream instore;
  instore << (*x);
  std::string instr = instore.str();
  // using PyString_FromString truncates the output if whitespace is
  // encountered so we use Py_BuildValue and specify the length
  return Py_BuildValue("s#",instr.c_str(), instr.size());
}

#endif

#endif /* ifndef __SAGE_CCOBJECT_H__ */
