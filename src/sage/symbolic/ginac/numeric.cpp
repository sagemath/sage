/** @file numeric.cpp
 *
 *  This file contains the interface to the underlying numerics packages.
 *  Its most important design principle is to completely hide the inner
 *  working of that other package from the user of GiNaC.  It must either 
 *  provide implementation of arithmetic operators and numerical evaluation
 *  of special functions or implement the interface to the numericsd package. */

/*
 *  This is a modified version of code included with Ginac.  
 *  The modifications and modified version is:
 * 
 *      GiNaC-SAGE Copyright (C) 2008 William Stein
 *      Copyright (C) 2009 Burcin Erocal <burcin@erocal.org>
 *                (C) 2015-2017 Ralf Stephan <ralf@ark.in-berlin.de>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


/*  The original copyright:
 * 
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "gmp.h"
#include "flint/fmpz.h"
#include "flint/fmpz_factor.h"

#include <cmath>
#include <utility>

#include "numeric.h"
#include "operators.h"
#include "ex.h"
#include "mul.h"
#include "power.h"
#include "function.h"
#include "archive.h"
#include "tostring.h"
#include "utils.h"
#include "../../cpython/pycore_long.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-register"
#include "factory/factory.h"
#pragma clang diagnostic pop

//#define Logging_refctr
#if defined(Logging_refctr)
#undef Py_INCREF
#undef Py_DECREF

#define Py_INCREF(op) (                         \
    _Py_INC_REFTOTAL  _Py_REF_DEBUG_COMMA       \
    ((PyObject*)(op))->ob_refcnt++) ; \
std::cerr << "+ " << long(op) << ", " << Py_REFCNT(op) << ", " << Py_TYPE(op)->tp_name << std::endl; std::cerr.flush();

#define Py_DECREF(op)                                   \
    do {                                                \
std::cerr << "- " << long(op) << ", " << Py_REFCNT(op) << ", " << Py_TYPE(op)->tp_name << std::endl; std::cerr.flush(); \
        if (_Py_DEC_REFTOTAL  _Py_REF_DEBUG_COMMA       \
        --((PyObject*)(op))->ob_refcnt != 0)            \
            _Py_CHECK_REFCNT(op)                        \
        else                                            \
        _Py_Dealloc((PyObject *)(op));                  \
    } while (0); std::cerr << "Done" << std::endl; std::cerr.flush();
#endif

const GiNaC::numeric& to_numeric(GiNaC::ex& e)
{
        return GiNaC::ex_to<GiNaC::numeric>(e);
}

// Call the Python function f on *this as input and return the result
// as a PyObject*.
#define PY_RETURN(f)  PyObject *a = to_pyobject();		 \
  PyObject *ans = f(a);						 \
  Py_DECREF(a);                                                  \
  if (!ans) py_error("error calling function");			 \
  return ans; 

// Call the Python function f on *this and return the result
// as a PyObject*.
#define PY_RETURN2(f, b)  PyObject *aa = to_pyobject();	 \
  PyObject* bb = b.to_pyobject();				 \
  PyObject *ans = f(aa, bb);					 \
  if (!ans) py_error("error calling function");			 \
  Py_DECREF(aa); Py_DECREF(bb); return ans; 

// Call the Python function f on *this and b, and get back
// a 2-tuple (z,w).  Set c = w, where c should be a 
// reference in the caller, and return z.  This is used
// to return two inputs from a Python function call.  See
// its usage in code below. 
#define PY_RETURN3(f, b, c)						\
  PyObject *aa = to_pyobject();					\
  PyObject* bb = b.to_pyobject();					\
  PyObject *ans = f(aa, bb);						\
  if (!ans) py_error("error calling function");				\
  if (!PyTuple_CheckExact(ans) || PyTuple_GET_SIZE(ans) != 2) 		\
    py_error("error calling function -- return not a 2-tuple.");	\
  PyObject* z =  PyTuple_GET_ITEM(ans, 0); Py_INCREF(z);                \
  PyObject* w =  PyTuple_GET_ITEM(ans, 1); Py_INCREF(w);          \
  c = w;                                                         \
  Py_DECREF(aa); Py_DECREF(bb); Py_DECREF(ans);                   \
  return z; 

//#define DEBUG
//#define VERBOSE

#ifdef DEBUG
#define todo(s) std::cerr << "TODO: " << s << std::endl;
#define stub(s) { std::cerr << "Hit STUB: " << s << std::endl; throw std::runtime_error("stub"); }
#define fake(s) std::cerr << "fake: " << s << std::endl;
#define ASSERT(s, msg) if (!s) { std::cerr << "Failed assertion: " << msg << std::endl; }
#else
#define todo(s)
#define stub(s) { std::cerr << "** Hit STUB**: " << s << std::endl; throw std::runtime_error("stub"); }
#define fake(s)
#endif

#ifdef VERBOSE
#define verbose(s) std::cerr << s << std::endl;
#define verbose2(s,t) std::cerr << s << " " << t << std::endl;
#define verbose3(s,t,u) std::cerr << s << " " << t << ", " << u << std::endl;
#else
#define verbose(s)
#define verbose2(s,t)
#define verbose3(s,t,u)
#endif


//////////////////////////////////////////////////////////////
// Python Interface
//////////////////////////////////////////////////////////////

inline void py_error(const char* errmsg) {
        throw std::runtime_error((PyErr_Occurred() != nullptr) ? errmsg:
                        "pyerror() called but no error occurred!");
}

#define PyString_FromString PyUnicode_FromString

// The following variable gets changed to true once
// this library has been imported by the Python
// interpreter.  This is because the interpreter calls
// ginac_pyinit_I below, which sets this to true.
// Once this is done we can call all the py_* functions
// defined in Cython modules, which is often much faster
// than what we can do without those. 
static bool initialized = false;

static PyObject* pyfunc_Float = nullptr;

void ginac_pyinit_Float(PyObject* f) {
        Py_INCREF(f);
        pyfunc_Float = f;
}

void ginac_pyinit_I(PyObject* z) {
        initialized = true;
        Py_INCREF(z);
        GiNaC::I = z; // I is a global constant defined below.
}

PyObject* RR_get()
{
        static PyObject* ptr = nullptr;
        if (ptr == nullptr) {
                PyObject* m = PyImport_ImportModule("sage.rings.real_mpfr");
                if (m == nullptr)
                        py_error("Error importing sage.rings.real_mpfr");
                ptr = PyObject_GetAttrString(m, "RR");
                if (ptr == nullptr)
                        py_error("Error getting RR attribute");
                Py_INCREF(ptr);
        }
        return ptr;
}

PyObject* CC_get()
{
        static PyObject* ptr = nullptr;
        if (ptr)
                return ptr;
        PyObject* m = PyImport_ImportModule("sage.rings.cc");
        if (m == nullptr)
                py_error("Error importing sage.rings.cc");
        ptr = PyObject_GetAttrString(m, "CC");
        if (ptr == nullptr)
                py_error("Error getting CC attribute");
        Py_INCREF(ptr);
        return ptr;
}

static PyObject* pyfunc_Integer = nullptr;

void ginac_pyinit_Integer(PyObject* f) {
        Py_INCREF(f);
        pyfunc_Integer = f;
}

PyObject* Integer_pyclass()
{
        if (not initialized)
                throw std::runtime_error("can't happen");
        static PyObject *ptr = nullptr;
        if (ptr)
                return ptr;
        PyObject *m = PyImport_ImportModule("sage.rings.integer");
        if (m == nullptr)
                py_error("Error importing sage.rings.integer");
        ptr = PyObject_GetAttrString(m, "Integer");
        if (ptr == nullptr)
                py_error("Error getting Integer attribute");
        return ptr;
}

PyObject* Integer(const long int& x) {
        if (initialized)
                return GiNaC::py_funcs.py_integer_from_long(x);

        PyObject *Integer = Integer_pyclass();
        PyObject* ans = PyObject_CallFunction(Integer, const_cast<char*> ("l"), x);
        return ans;
}

int precision(const GiNaC::numeric& num, PyObject*& a_parent) {
        int prec = 0;
        PyObject* the_parent = a_parent;
        if (a_parent == nullptr) {
                PyObject* m = PyImport_ImportModule("sage.structure.element");
                if (m == nullptr)
                        py_error("Error importing element");
                PyObject* f = PyObject_GetAttrString(m, "parent");
                if (f == nullptr)
                        py_error("Error getting parent attribute");
                PyObject* obj = num.to_pyobject();
                the_parent = PyObject_CallFunctionObjArgs(f, obj, NULL);
                Py_DECREF(obj);
                Py_DECREF(f);
                Py_DECREF(m);
        }
        else if (PyDict_Check(a_parent)) {
                PyObject* pkey = PyString_FromString(const_cast<char*>("parent"));
                the_parent = PyDict_GetItem(a_parent, pkey);
                Py_DECREF(pkey);
        }
        PyObject *mprec = nullptr;
        if (the_parent != nullptr)
                mprec = PyObject_CallMethod(the_parent, const_cast<char*>("precision"), NULL);
        if (mprec == nullptr) {
                prec = 53;
                PyErr_Clear();
        }
        else {
                prec = PyLong_AsLong(mprec);
                Py_DECREF(mprec);
        }
        if (a_parent == nullptr) {
                a_parent = PyDict_New();
                PyDict_SetItemString(a_parent, "parent", the_parent);
        }
        return prec;
}

PyObject* CBF(int res) {
        PyObject* m = PyImport_ImportModule("sage.rings.complex_arb");
        if (m == nullptr)
                py_error("Error importing sage.rings.complex_arb");
        PyObject* f = PyObject_GetAttrString(m, "ComplexBallField");
        if (f == nullptr)
                py_error("Error getting ComplexBallField attribute");
        PyObject* list = PyTuple_New(1);
        if (list == nullptr)
                throw(std::runtime_error("GiNaC::CBF(): PyTuple_New returned NULL"));
        PyObject *iobj = Integer(res);
        int ret = PyTuple_SetItem(list, 0, iobj);
        if (ret != 0)
                throw(std::runtime_error("GiNaC::CBF(): PyTuple_SetItem unsuccessful"));
        PyObject *obj = PyObject_Call(f, list, NULL);
        if (obj == nullptr)
                throw(std::runtime_error("GiNaC::CBF(): PyObject_Call unsuccessful"));
        Py_DECREF(m);
        Py_DECREF(f);
        Py_DECREF(list);
        return obj;
}

// Convert a to field elt, return a.meth()
PyObject* CallBallMethod0Arg(PyObject* field, const char* meth, const GiNaC::numeric& a) {
        PyObject* list1 = PyTuple_New(1);
        if (list1 == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyTuple_New returned NULL"));
        PyObject *aobj = a.to_pyobject();
        int ret = PyTuple_SetItem(list1, 0, aobj);
        if (ret != 0)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyTuple_SetItem unsuccessful"));
        PyObject *aball = PyObject_Call(field, list1, NULL);
        if (aball == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyObject_Call unsuccessful"));
        PyObject* name = PyString_FromString(meth);
        PyObject* ret_ball = PyObject_CallMethodObjArgs(aball, name, NULL);
        if (ret_ball == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyObject_CallMethodObjArgs unsuccessful"));
        Py_DECREF(list1);
        Py_DECREF(aball);
        Py_DECREF(name);

        return ret_ball;
}

// Convert a to field elt, return a.meth(b)
PyObject* CallBallMethod1Arg(PyObject* field, const char* meth, const GiNaC::numeric& a, const GiNaC::numeric& b) {
        PyObject* list1 = PyTuple_New(1);
        if (list1 == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyTuple_New returned NULL"));
        PyObject *aobj = a.to_pyobject();
        int ret = PyTuple_SetItem(list1, 0, aobj);
        if (ret != 0)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyTuple_SetItem unsuccessful"));
        PyObject *aball = PyObject_Call(field, list1, NULL);
        if (aball == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyObject_Call unsuccessful"));
        PyObject* list2 = PyTuple_New(1);
        if (list2 == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyTuple_New returned NULL"));
        PyObject *bobj = b.to_pyobject();
        ret = PyTuple_SetItem(list2, 0, bobj);
        if (ret != 0)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyTuple_SetItem unsuccessful"));
        PyObject *bball = PyObject_Call(field, list2, NULL);
        if (bball == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyObject_Call unsuccessful"));

        PyObject* name = PyString_FromString(meth);
        PyObject* ret_ball = PyObject_CallMethodObjArgs(aball, name, bball, NULL);
        if (ret_ball == nullptr)
                throw(std::runtime_error("GiNaC::CallBallMethod1Arg(): PyObject_CallMethodObjArgs unsuccessful"));
        Py_DECREF(list1);
        Py_DECREF(list2);
        Py_DECREF(aball);
        Py_DECREF(bball);
        Py_DECREF(name);

        return ret_ball;
}

PyObject* CoerceBall(PyObject* ball, int prec) {
        PyObject* m = PyImport_ImportModule("sage.rings.complex_mpfr");
        if (m == nullptr)
                py_error("Error importing sage.rings.complex_mpfr");
        PyObject* f = PyObject_GetAttrString(m, "ComplexField");
        if (f == nullptr)
                py_error("Error getting ComplexField attribute");
        PyObject* list1 = PyTuple_New(1);
        if (list1 == nullptr)
                throw(std::runtime_error("GiNaC: PyTuple_New returned NULL"));
        PyObject *iobj = Integer(prec);
        int ret = PyTuple_SetItem(list1, 0, iobj);
        if (ret != 0)
                throw(std::runtime_error("GiNaC: PyTuple_SetItem unsuccessful"));
        PyObject *field = PyObject_CallObject(f, list1);
        if (field == nullptr)
                throw(std::runtime_error("GiNaC: PyObject_Call unsuccessful"));
        PyObject* list2 = PyTuple_New(1);
        if (list2 == nullptr)
                throw(std::runtime_error("GiNaC: PyTuple_New returned NULL"));
        ret = PyTuple_SetItem(list2, 0, ball);
        if (ret != 0)
                throw(std::runtime_error("GiNaC: PyTuple_SetItem unsuccessful"));
        PyObject *celt = PyObject_CallObject(field, list2);
        if (celt == nullptr)
                throw(std::runtime_error("GiNaC: PyObject_Call unsuccessful"));

        Py_INCREF(ball);
        Py_DECREF(list1);
        Py_DECREF(list2);
        Py_DECREF(field);
        Py_DECREF(f);
        Py_DECREF(m);

        PyObject* bret = PyObject_CallMethod(celt, const_cast<char*>("is_real"), NULL);
        if (PyObject_IsTrue(bret)) {
                PyObject* rret = PyObject_CallMethod(celt, const_cast<char*>("real"), NULL);
                Py_DECREF(bret);
                Py_DECREF(celt);
                return rret;
        }
        Py_DECREF(bret);
        return celt;
}

namespace GiNaC {

numeric I;

PyObject* common_parent(const numeric& x, const numeric& y)
{
        PyObject* m = PyImport_ImportModule("sage.structure.element");
        if (m == nullptr)
                py_error("Error importing coerce");
        PyObject* f = PyObject_GetAttrString(m, "coercion_model");
        if (f == nullptr)
                py_error("Error getting coercion_model attribute");

        PyObject* mname = PyString_FromString("common_parent");
        PyObject* xobj = x.to_pyobject();
        PyObject* yobj = y.to_pyobject();
        PyObject* ret = PyObject_CallMethodObjArgs(f, mname,
                                                   xobj, yobj, NULL);
        if (ret == nullptr)
                throw(std::runtime_error("GiNaC::common_parent: PyObject_CallMethodObjArgs unsuccessful"));
        Py_DECREF(xobj);
        Py_DECREF(yobj);
        Py_DECREF(m);
        Py_DECREF(f);
        Py_DECREF(mname);
        return ret;
}

const numeric numeric::arbfunc_0arg(const char* name, PyObject* parent) const
{
        int prec = precision(*this, parent);
        PyObject* field = CBF(prec+15);
        PyObject* ret = CallBallMethod0Arg(field, name, *this);
        Py_DECREF(field);

        numeric rnum(ret);
        return ex_to<numeric>(rnum.evalf(0, parent));
}

static PyObject* py_tuple_from_numvector(const std::vector<numeric>& vec)
{
        PyObject* list = PyTuple_New(vec.size());
        if (list == nullptr)
                throw(std::runtime_error("py_list_from_numvector(): PyList_New returned NULL"));
        int pos = 0;
        for (const numeric& num : vec) {
                PyObject *numobj = num.to_pyobject();
                int ret = PyTuple_SetItem(list, pos++, numobj);
                if (ret != 0)
                        throw(std::runtime_error("py_list_from_numvector(): PyList_Append unsuccessful"));
        }
        return list;
}

///////////////////////////////////////////////////////////////////////////////
// class numeric
///////////////////////////////////////////////////////////////////////////////

PyObject* ZERO = PyLong_FromLong(0); // todo: never freed
PyObject* ONE = PyLong_FromLong(1); // todo: never freed
PyObject* TWO = PyLong_FromLong(2); // todo: never freed

std::ostream& operator<<(std::ostream& os, const numeric& s) {
        switch (s.t) {
        case LONG:
                return os << s.v._long;
        case MPZ: {
                std::vector<char> cp(2 + mpz_sizeinbase(s.v._bigint, 10));
                mpz_get_str(&cp[0], 10, s.v._bigint);
                return os << &cp[0];
        }
        case MPQ: {
                size_t size = mpz_sizeinbase(mpq_numref(s.v._bigrat), 10)
                + mpz_sizeinbase(mpq_denref(s.v._bigrat), 10) + 5;
                std::vector<char> cp(size);
                mpq_get_str(&cp[0], 10, s.v._bigrat);
                return os << &cp[0];
        }
        case PYOBJECT:
                return os << *py_funcs.py_repr(s.v._pyobject, 0);
        default:
                stub("operator <<: type not yet handled");
        }
}

const numeric& numeric::operator=(const numeric& x) {
        switch (t) {
        case LONG:
                break;
        case MPZ:
                mpz_clear(v._bigint);
                break;
        case MPQ:
                mpq_clear(v._bigrat);
                break;
        case PYOBJECT:
                Py_DECREF(v._pyobject);
                break;
        }
        t = x.t;
        hash = x.hash;
        switch (x.t) {
        case LONG:
                v._long = x.v._long;
                break;
        case MPZ:
                mpz_init(v._bigint);
                mpz_set(v._bigint, x.v._bigint);
                break;
        case MPQ:
                mpq_init(v._bigrat);
                mpq_set(v._bigrat, x.v._bigrat);
                break;
        case PYOBJECT:
                v = x.v;
                Py_INCREF(v._pyobject);
                break;
        }
        return *this;
}

void numeric::swap(numeric& other)
{
        std::swap(t, other.t);
        std::swap(v, other.v);
        std::swap(hash, other.hash);
        std::swap(is_hashable, other.is_hashable);
}

int numeric::compare_same_type(const numeric& right) const {
        verbose("compare_same_type");
        if (t != right.t) {
                int ret;
                if (t == MPZ and right.t == MPQ) {
                        ret = -mpq_cmp_z(right.v._bigrat, v._bigint);
                }
                else if (t == MPQ and right.t == MPZ) {
                        ret = mpq_cmp_z(v._bigrat, right.v._bigint);
                }
                else if (t == LONG and right.t == MPZ) {
                        ret = -mpz_cmp_si(right.v._bigint, v._long);
                }
                else if (t == LONG and right.t == MPQ) {
                        ret = -mpq_cmp_si(right.v._bigrat, v._long, 1);
                }
                else if (t == MPZ and right.t == LONG) {
                        ret = mpz_cmp_si(v._bigint, right.v._long);
                }
                else if (t == MPQ and right.t == LONG) {
                        ret = mpq_cmp_si(v._bigrat, right.v._long, 1);
                }
                else {
                        numeric a, b;
                        coerce(a, b, *this, right);
                        return a.compare_same_type(b);
                }
                if (ret > 0)
                        ret = 1;
                else if (ret < 0)
                        ret = -1;
                return ret;
        }
        int ret;
        switch (t) {
        case LONG:
                if (v._long > right.v._long)
                        return 1;
                else if (v._long < right.v._long)
                        return -1;
                else
                        return 0;
        case MPZ:
                ret = mpz_cmp(v._bigint, right.v._bigint);
                if (ret > 0)
                        ret = 1;
                else if (ret < 0)
                        ret = -1;
                return ret;
        case MPQ:
                ret = mpq_cmp(v._bigrat, right.v._bigrat);
                if (ret > 0)
                        ret = 1;
                else if (ret < 0)
                        ret = -1;
                return ret;
        case PYOBJECT: {
                int result = PyObject_RichCompareBool(v._pyobject,
                right.v._pyobject, Py_LT);
                if (result == 1)
                        return -1;
                if (result == -1)
                        py_error("richcmp failed");
                else { // result == 0
                        result = PyObject_RichCompareBool(v._pyobject,
                        right.v._pyobject, Py_GT);
                        if (result == -1)
                                py_error("richcmp failed");
                }
                return result;
        }
        default:
                stub("invalid type: compare_same_type type not handled");
        }
}

#if PY_MAJOR_VERSION < 3 || defined(PYPY_VERSION)
#define hash_bits (8 * sizeof(void*))
#else
#define hash_bits _PyHASH_BITS
#endif
#define limb_bits (8 * sizeof(mp_limb_t))

/* By convention hashes of PyObjects must be identical
   with their Python hashes, this applies to our MPZ
   and MPQ objects too. This implementation is copied
   from sage.libs.gmp.pylong.mpz_pythonhash. */
static long _mpz_pythonhash_raw(mpz_t the_int)
{
    if (mpz_sgn(the_int) == 0)
        return 0;

    mp_limb_t modulus = ((((mp_limb_t)(1) << (hash_bits - 1)) - 1) * 2) + 1;
    mp_limb_t h=0, x, y;
    size_t n = mpz_size(the_int);
    unsigned int r;
    for (unsigned i=0; i<n; ++i) {
        x = mpz_getlimbn(the_int, i);
        if (limb_bits == hash_bits)
            y = x;
        else {
            r = (limb_bits * i) % hash_bits;
            y = (x << r) & modulus;
            y += (x >> (hash_bits - r)) & modulus;
            if (r > 2 * hash_bits - limb_bits)
                y += (x >> (2 * hash_bits - r));
            if (y > modulus)
                y -= modulus;
        }
        if (h < modulus - y)
            h = h + y;
        else
            h = h - (modulus - y);
    }
    if (limb_bits == hash_bits && h == 0)
        h = -1;

    if (mpz_sgn(the_int) < 0)
        return -h;
    return h;
}

static long _mpz_pythonhash(mpz_t the_int)
{
    long h = _mpz_pythonhash_raw(the_int);
    if (h == -1)
        return -2;
    return h;
}

static long _mpq_pythonhash(mpq_t the_rat)
{
    mpq_t rat;
    mpq_init(rat);
    mpq_set(rat, the_rat);
    long n = _mpz_pythonhash_raw(mpq_numref(rat));
    long d = _mpz_pythonhash_raw(mpq_denref(rat));
    if (d != 1L)
        n = n + (d-1) * 7461864723258187525;
    mpq_clear(rat);
    if (n == -1)
        return -2;
    return n;
}


// Initialize an mpz_t from a Python long integer
static void _mpz_set_pylong(mpz_t z, PyLongObject* l)
{
    Py_ssize_t pylong_size = _PyLong_DigitCount(l);

    mpz_import(z, pylong_size, -1, sizeof(digit), 0,
               8*sizeof(digit) - PyLong_SHIFT, ob_digit(l));

    if (_PyLong_IsNegative(l))
        mpz_neg(z, z);
}



///////////////////////////////////////////////////////////////////////////////
// class numeric
///////////////////////////////////////////////////////////////////////////////

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(numeric, basic,
        print_func<print_context>(&numeric::do_print).
        print_func<print_latex>(&numeric::do_print_latex).
        print_func<print_tree>(&numeric::do_print_tree).
        print_func<print_python_repr>(&numeric::do_print_python_repr))


numeric::numeric() : basic(&numeric::tinfo_static), t(LONG), v(0) {
        setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

// public

numeric::numeric(const numeric& other) : basic(&numeric::tinfo_static) {
        t = other.t;
        hash = other.hash;
        switch (t) {
        case LONG:
                v._long = other.v._long;
                return;
        case PYOBJECT:
                v = other.v;
                Py_INCREF(v._pyobject);
                return;
        case MPZ:
                mpz_init(v._bigint);
                mpz_set(v._bigint, other.v._bigint);
                return;
        case MPQ:
                mpq_init(v._bigrat);
                mpq_set(v._bigrat, other.v._bigrat);
                return;
        }
}

bool is_Sage_Integer(PyObject* o)
{
        PyObject *Integer = Integer_pyclass();
        int res = PyObject_IsInstance(o, Integer);
        if (res < 0)
                py_error("Error testing Integer attribute");
        return res != 0;
}

void set_from(Type& t, Value& v, long& hash, mpz_t bigint)
{
        if (mpz_fits_sint_p(bigint)) {
                t = LONG;
                v._long = mpz_get_si(bigint);
                hash = (v._long==-1) ? -2 : v._long;
        }
        else {
                t = MPZ;
                mpz_init_set(v._bigint, bigint);
                hash = _mpz_pythonhash(v._bigint);
        }
}

void set_from(Type& t, Value& v, long& hash, mpq_t bigrat)
{
        t = MPQ;
        mpq_init(v._bigrat);
        mpq_set(v._bigrat, bigrat);
        hash = _mpq_pythonhash(v._bigrat);
}

numeric::numeric(PyObject* o, bool force_py) : basic(&numeric::tinfo_static) {
        if (o == nullptr) py_error("Error");
        if (not force_py) {
                if (PyLong_Check(o)) {
                    t = MPZ;
                    mpz_init(v._bigint);
                    _mpz_set_pylong(v._bigint, reinterpret_cast<PyLongObject*>( o));
                    hash = _mpz_pythonhash(v._bigint);
                    setflag(status_flags::evaluated | status_flags::expanded);
                    Py_DECREF(o);
                    return;
                }
                if (initialized) {
                        if (is_Sage_Integer(o)) {
                                mpz_ptr mpz = py_funcs.py_mpz_from_integer(o);
                                set_from(t, v, hash, mpz);
                                Py_DECREF(o);
                                setflag(status_flags::evaluated
                                                | status_flags::expanded);
                                return;
                        }
                        if (py_funcs.py_is_Rational(o) != 0) {
                                mpq_ptr mpq = py_funcs.py_mpq_from_rational(o);
                                set_from(t, v, hash, mpq);
                                Py_DECREF(o);
                                setflag(status_flags::evaluated
                                                | status_flags::expanded);
                                return;
                        }
                }
        }

        t = PYOBJECT;
        hash = PyObject_Hash(o);
        if (hash == -1 && (PyErr_Occurred() != nullptr)) {
            // error is thrown on first hash request
            PyErr_Clear();
            is_hashable = false;
            }
        v._pyobject = o; // STEAL a reference
        setflag(status_flags::evaluated | status_flags::expanded);

}

numeric::numeric(mpz_t bigint) : basic(&numeric::tinfo_static)
{
        set_from(t, v, hash, bigint);
        mpz_clear(bigint);
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(mpq_t bigrat) : basic(&numeric::tinfo_static)
{
        set_from(t, v, hash, bigrat);
        mpq_clear(bigrat);
        setflag(status_flags::evaluated | status_flags::expanded);
}

/** Constructor for rational numerics a/b.
 *
 *  @exception overflow_error (division by zero) */
numeric::numeric(long num, long den) : basic(&numeric::tinfo_static) {
        if (den == 0)
                throw std::overflow_error("numeric::div(): division by zero");
        if ((num%den) == 0) {
                t = LONG;
                v._long = num/den;
                hash = (v._long==-1) ? -2 : v._long;
        }
        else
        {
                t = MPQ;
                mpq_init(v._bigrat);
                mpq_set_si(v._bigrat, num, den);
                mpq_canonicalize(v._bigrat);
                hash = _mpq_pythonhash(v._bigrat);
        }
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::numeric(double d) : basic(&numeric::tinfo_static) {
        t = PYOBJECT;
        if ((v._pyobject = PyFloat_FromDouble(d)) == nullptr)
                py_error("Error creating double");
        setflag(status_flags::evaluated | status_flags::expanded);
}

numeric::~numeric() {
        switch (t) {
        case LONG:
                return;
        case PYOBJECT:
                Py_DECREF(v._pyobject);
                return;
        case MPZ:
                mpz_clear(v._bigint);
                return;
        case MPQ:
                mpq_clear(v._bigrat);
                return;
        }
}

//////////
// archiving
//////////

numeric::numeric(const archive_node &n, lst &sym_lst) :
inherited(n, sym_lst) {
        // get type information
        unsigned int t_tmp;
        if (!n.find_unsigned(std::string("T"), t_tmp))
                throw std::runtime_error("archive error: cannot read type info");
        t = Type(t_tmp);
        std::string str;
        // read object to a string
        if (!n.find_string("S", str))
                throw (std::runtime_error("archive error: cannot read object data"));
        PyObject *arg;
        switch (t) {
        case LONG:
                v._long = std::stol(str);
                hash = (v._long == -1) ? -2 : v._long;
                return;
        case MPZ:
                mpz_init(v._bigint);
                mpz_set_str(v._bigint, str.c_str(), 10);
                hash = _mpz_pythonhash(v._bigint);
                return;
        case MPQ:
                mpq_init(v._bigrat);
                mpq_set_str(v._bigrat, str.c_str(), 10);
                hash = _mpq_pythonhash(v._bigrat);
                return;
        case PYOBJECT:
                // read pickled python object to a string
                if (!n.find_string("S", str))
                        throw(std::runtime_error("archive error: cannot read pyobject data"));
                arg = Py_BuildValue("s#", str.c_str(), str.size());
                // unpickle
                v._pyobject = py_funcs.py_loads(arg);
                Py_DECREF(arg);
                if (PyErr_Occurred() != nullptr) {
                        throw(std::runtime_error("archive error: caught exception in py_loads"));
                }
                hash = PyObject_Hash(v._pyobject);
                if (hash == -1 && (PyErr_Occurred() != nullptr)) {
                        PyErr_Clear();
                        is_hashable = false;
                }
                Py_INCREF(v._pyobject);
                return;
        default:
                stub("unarchiving numeric");
                return;
        }
        setflag(status_flags::evaluated | status_flags::expanded);
}

void numeric::archive(archive_node &n) const {
        // store type information
        n.add_unsigned("T", t);

        // create a string representation of this object
        std::string *tstr;
        switch (t) {
        case LONG:
                tstr = new std::string(std::to_string(v._long));
                break;
        case MPZ: {
                std::vector<char> cp(2 + mpz_sizeinbase(v._bigint, 10));
                mpz_get_str(&cp[0], 10, v._bigint);
                tstr = new std::string(&cp[0]);
                break;
        }
        case MPQ: {
                size_t size = mpz_sizeinbase(mpq_numref(v._bigrat), 10)
                + mpz_sizeinbase(mpq_denref(v._bigrat), 10) + 5;
                std::vector<char> cp(size);
                mpq_get_str(&cp[0], 10, v._bigrat);
                tstr = new std::string(&cp[0]);
                break;
        }
        case PYOBJECT:
                tstr = py_funcs.py_dumps(v._pyobject);
                if (PyErr_Occurred() != nullptr) {
                        throw(std::runtime_error("archive error: exception in py_dumps"));
                }
                break;
        default:
                stub("archive numeric");
        }

        n.add_string("S", *tstr);
        delete tstr;
        inherited::archive(n);
}

//////////
// functions overriding virtual functions from base classes
//////////

template<typename T1, typename T2>
static inline bool coerce(T1& dst, const T2& arg);

void numeric::print_numeric(const print_context & c, const char* /*unused*/,
                            const char* /*unused*/, const char* /*unused*/, const char* /*unused*/,
                            unsigned level, bool latex = false) const {
        std::string* out;
        PyObject* obj = to_pyobject();
        if (latex) {
                out = py_funcs.py_latex(obj, level);
        } else {
                out = py_funcs.py_repr(obj, level);
        }
        c.s << *out;
        Py_DECREF(obj);
        delete out;
}

void numeric::do_print(const print_context & c, unsigned level) const {
        print_numeric(c, "(", ")", "I", "*", level, false);
}

void numeric::do_print_latex(const print_latex & c, unsigned level) const {
        print_numeric(c, "{(", ")}", "i", " ", level, true);
}

void numeric::do_print_python_repr(const print_python_repr & c, unsigned level) const {
        c.s << class_name() << "('";
        print_numeric(c, "(", ")", "I", "*", level);
        c.s << "')";
}

static std::string dbgstring(const numeric& num)
{
        std::string ts;
        switch (num.t) {
        case LONG:
                ts = "LONG";
                break;
        case MPZ:
                ts = "MPZ";
                break;
        case MPQ:
                ts = "MPQ";
                break;
        case PYOBJECT:
                {
                ts = "PYOBJECT: ";
                PyObject *to = PyObject_Type(num.v._pyobject);
                if (to == nullptr)
                        ts.append("NULL");
                else {
                        PyObject *ro = PyObject_Repr(to);
                        if (ro == nullptr)
                                ts.append("NULL");
                        else {
                                ts.append(*py_funcs.py_repr(ro, 0));
                                Py_DECREF(ro);
                        }
                        Py_DECREF(to);
                }
                break;
        }
        default: stub("typestr()");
        }

        std::stringstream ss;
        ss << num << " (numeric)" << " @" << &num << std::hex << ", hash=0x"
                << num.hash << ", flags=0x" << num.flags << std::dec
                << ", type " << ts;
        return ss.str();
}

void numeric::do_print_tree(const print_tree & c, unsigned level) const
{
        c.s << std::string(level, ' ') << dbgstring(*this) << std::endl;
}

void numeric::dbgprint() const
{
        std::cerr << dbgstring(*this) << std::endl;
}

bool numeric::info(unsigned inf) const {
        switch (inf) {
        case info_flags::numeric:
        case info_flags::polynomial:
        case info_flags::rational_function:
        case info_flags::expanded:
                return true;
        case info_flags::real:
                return is_real();
        case info_flags::rational:
        case info_flags::rational_polynomial:
                return is_rational();
        case info_flags::crational:
        case info_flags::crational_polynomial:
                return is_crational();
        case info_flags::integer:
        case info_flags::integer_polynomial:
                return is_integer();
        case info_flags::cinteger:
        case info_flags::cinteger_polynomial:
                return is_cinteger();
        case info_flags::positive:
                return is_positive();
        case info_flags::negative:
                return is_negative();
        case info_flags::nonnegative:
                return is_zero() || is_positive();
        case info_flags::nonzero:
                return not is_zero();
        case info_flags::posint:
                return is_pos_integer();
        case info_flags::negint:
                return is_integer() && is_negative();
        case info_flags::nonnegint:
                return is_nonneg_integer();
        case info_flags::even:
                return is_even();
        case info_flags::odd:
                return is_odd();
        case info_flags::prime:
                return is_prime();
        case info_flags::algebraic:
                return !is_real();
        case info_flags::infinity:
                return false;
        case info_flags::inexact:
                return not is_exact();
        default:
                throw(std::runtime_error("numeric::info()"));
        }
        return false;
}

bool numeric::is_polynomial(const ex & var) const {
        return true;
}

numeric numeric::degree(const ex & s) const {
        // In sage deg (0) != 0 !!!
        return *_num0_p;
}

numeric numeric::ldegree(const ex & s) const {
        return *_num0_p;
}

ex numeric::coeff(const ex & s, const ex & n) const {
        return n.is_zero() ? * this : _ex0;
}

/** Disassemble real part and imaginary part to scan for the occurrence of a
 *  single number.  Also handles the imaginary unit.  It ignores the sign on
 *  both this and the argument, which may lead to what might appear as funny
 *  results:  (2+I).has(-2) -> true.  But this is consistent, since we also
 *  would like to have (-2+I).has(2) -> true and we want to think about the
 *  sign as a multiplicative factor. */
bool numeric::has(const ex &other, unsigned /*options*/) const {
        if (!is_exactly_a<numeric>(other))
                return false;
        const numeric &o = ex_to<numeric>(other);
        if (this->is_equal(o) || this->is_equal(-o))
                return true;
        if (o.imag().is_zero()) { // e.g. scan for 3 in -3*I
                if (!this->real().is_equal(*_num0_p))
                        if (this->real().is_equal(o) || this->real().is_equal(-o))
                                return true;
                if (!this->imag().is_equal(*_num0_p))
                        if (this->imag().is_equal(o) || this->imag().is_equal(-o))
                                return true;
                return false;
        } 
        if (o.is_equal(I)) // e.g scan for I in 42*I
                return !this->is_real();
        if (o.real().is_zero() // e.g. scan for 2*I in 2*I+1
            and not this->imag().is_equal(*_num0_p)
            and (this->imag().is_equal(o * I)
                 or this->imag().is_equal(-o * I)))
                return true;
        
        return false;
}

numeric numeric::conj() const {
        switch (t) {
        case LONG:
        case MPZ:
        case MPQ: return *this;
        case PYOBJECT: {
                PyObject *obj = PyObject_GetAttrString(v._pyobject,
                "conjugate");
                if (obj == nullptr)
                        return *this;
                obj = PyObject_CallObject(obj, NULL);
                if (obj == nullptr)
                        py_error("Error calling Python conjugate");
                return obj;
        }
        default:
                stub("invalid type: ::conjugate() type not handled");
        }
}

ex numeric::real_part() const {
        return real();
}

ex numeric::imag_part() const {
        return imag();
}

// protected

/** This method establishes a canonical order on all numbers.  For complex
 *  numbers this is not possible in a mathematically consistent way but we need
 *  to establish some order and it ought to be fast.  So we simply define it
 *  to be compatible with our method csgn.
 *
 *  @return csgn(*this-other) */
int numeric::compare_same_type(const basic &other) const {
        GINAC_ASSERT(is_exactly_a<numeric>(other));
        const numeric &o = static_cast<const numeric &> (other);

        return compare_same_type(o);
}

bool numeric::is_equal_same_type(const basic &other) const {
        GINAC_ASSERT(is_exactly_a<numeric>(other));
        const numeric &o = static_cast<const numeric &> (other);

        return this->is_equal(o);
}

long numeric::calchash() const {
        switch (t) {
        case LONG:
                return (v._long == -1) ? -2 : v._long;
        case MPZ:
        case MPQ:
        case PYOBJECT:
                if (is_hashable)
                        return hash;
                throw(std::runtime_error("Python object not hashable"));
        default:
                stub("invalid type: ::hash() type not handled");
        }
}


//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// public

/** Numerical addition method.  Adds argument to *this and returns result as
 *  a numeric object. */
const numeric numeric::add(const numeric &other) const {
        verbose("operator+");
        if (other.is_zero())
                return *this;
        if (is_zero())
                return other;
        if (t != other.t) {
                if (t == MPZ and other.t == MPQ) {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_z(bigrat, v._bigint);
                        mpq_add(bigrat, bigrat, other.v._bigrat);
                        return bigrat;
                }
                if (t == MPQ and other.t == MPZ) {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_z(bigrat, other.v._bigint);
                        mpq_add(bigrat, bigrat, v._bigrat);
                        return bigrat;
                }

                numeric a, b;
                coerce(a, b, *this, other);
                return a + b;
        }
        switch (t) {
        case LONG: {
                if ((v._long > 0
                    and v._long < std::numeric_limits<long>::max() / 2
                    and other.v._long < std::numeric_limits<long>::max() / 2)
                or (v._long < 0
                   and v._long > std::numeric_limits<long>::min() / 2
                   and other.v._long > std::numeric_limits<long>::min() / 2))
                        return v._long + other.v._long;
                mpz_t bigint;
                mpz_init_set_si(bigint, v._long);
                if (other.v._long < 0)
                        mpz_sub_ui(bigint, bigint, -other.v._long);
                else
                        mpz_add_ui(bigint, bigint, other.v._long);
                return bigint;
        }
        case MPZ:
                mpz_t bigint;
                mpz_init(bigint);
                mpz_add(bigint, v._bigint, other.v._bigint);
                return bigint;
        case MPQ:
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_add(bigrat, v._bigrat, other.v._bigrat);
                return bigrat;
        case PYOBJECT:
                return PyNumber_Add(v._pyobject, other.v._pyobject);
        default:
                stub("invalid type: operator+() type not handled");
        }
}

/** Numerical subtraction method.  Subtracts argument from *this and returns
 *  result as a numeric object. */
const numeric numeric::sub(const numeric &other) const {
        verbose("operator-");
        if (other.is_zero())
                return *this;
        if (is_zero())
                return other.negative();
        if (t != other.t) {
                if (t == MPZ and other.t == MPQ) {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_z(bigrat, v._bigint);
                        mpq_sub(bigrat, bigrat, other.v._bigrat);
                        return bigrat;
                }
                if (t == MPQ and other.t == MPZ) {
                        mpq_t bigrat, tmp;
                        mpq_init(bigrat);
                        mpq_init(tmp);
                        mpq_set(bigrat, v._bigrat);
                        mpq_set_z(tmp, other.v._bigint);
                        mpq_sub(bigrat, bigrat, tmp);
                        mpq_clear(tmp);
                        return bigrat;
                }

                numeric a, b;
                coerce(a, b, *this, other);
                return a - b;
        }
        switch (t) {
        case LONG: {
                if ((v._long > 0
                    and v._long < std::numeric_limits<long>::max() / 2
                    and -other.v._long < std::numeric_limits<long>::max() / 2)
                or (v._long < 0
                   and v._long > std::numeric_limits<long>::min() / 2
                   and -other.v._long > std::numeric_limits<long>::min() / 2))
                        return v._long - other.v._long;
                mpz_t bigint;
                mpz_init_set_si(bigint, v._long);
                if (other.v._long < 0)
                        mpz_add_ui(bigint, bigint, -other.v._long);
                else
                        mpz_sub_ui(bigint, bigint, other.v._long);
                return bigint;
        }
        case MPZ:
                mpz_t bigint;
                mpz_init(bigint);
                mpz_sub(bigint, v._bigint, other.v._bigint);
                return bigint;
        case MPQ:
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_sub(bigrat, v._bigrat, other.v._bigrat);
                return bigrat;
        case PYOBJECT:
                return PyNumber_Subtract(v._pyobject, other.v._pyobject);
        default:
                stub("invalid type: operator-() type not handled");
        }
}

#if defined __has_builtin
#  if __has_builtin (__builtin_smull_overflow)
#    define smull_overflow __builtin_smull_overflow
#  endif
#endif

#if !defined smull_overflow
static int smull_overflow(long a, long b, long *result) {
        static long lsqrt = std::lround(std::sqrt(std::numeric_limits<long>::max()));
        if (-lsqrt < a and a < lsqrt and -lsqrt < b and b < lsqrt) {
                // Will not overflow.
                *result = a * b;
                return 0;
        }
        // May overflow.
        return 1;
}
#endif

/** Numerical multiplication method.  Multiplies *this and argument and returns
 *  result as a numeric object. */
const numeric numeric::mul(const numeric &other) const {
        verbose("operator*");
        if ((is_zero() and t != PYOBJECT)
            or (other.is_zero() and other.t != PYOBJECT))
                return *_num0_p;
        if (other.is_one())
                return *this;
        if (is_one())
                return other;
        if (t != other.t) {
                numeric a, b;
                coerce(a, b, *this, other);
                return a * b;
        }
        switch (t) {
        case LONG: {
                long result;
                if (!smull_overflow(v._long, other.v._long, & result))
                        return result;
                // the multiplication overflowed, so use mpz
                mpz_t bigint;
                mpz_init_set_si(bigint, v._long);
                mpz_mul_si(bigint, bigint, other.v._long);
                return bigint;
        }
        case MPZ:
                mpz_t bigint;
                mpz_init(bigint);
                mpz_mul(bigint, v._bigint, other.v._bigint);
                return bigint;
        case MPQ:
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_mul(bigrat, v._bigrat, other.v._bigrat);
                return bigrat;
        case PYOBJECT:
                return PyNumber_Multiply(v._pyobject, other.v._pyobject);
        default:
                stub("invalid type: operator*() type not handled");
        }
}

/** Numerical division method.  Divides *this by argument and returns result as
 *  a numeric object.
 *
 *  @exception overflow_error (division by zero) */
const numeric numeric::div(const numeric &other) const {
        verbose("operator/");
        if (other.is_zero())
                throw std::overflow_error("numeric::div(): division by zero");
        if (is_zero())
                return *_num0_p;
        if (other.is_one())
                return *this;
        if (t != other.t) {
                numeric a, b;
                coerce(a, b, *this, other);
                return a / b;
        }
        switch (t) {
        case LONG: {
                if (v._long == 1 and other.v._long == 2)
                        return *_num1_2_p;
                if (other.v._long == -1)  // use multiplication to avoid possible overflow
                        return *this * -1;
                auto ld = std::div(v._long, other.v._long);
                if (ld.rem == 0)
                        return ld.quot;

                mpq_t bigrat, obigrat;
                mpq_init(bigrat);
                mpq_init(obigrat);
                mpq_set_si(bigrat, v._long, 1);
                mpq_set_si(obigrat, other.v._long, 1);
                mpq_div(bigrat, bigrat, obigrat);
                mpq_clear(obigrat);
                return bigrat;
        }
        case MPZ: {
                if (mpz_divisible_p(v._bigint, other.v._bigint)) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_divexact(bigint, v._bigint,
                        other.v._bigint);
                        return bigint;
                }
                mpq_t bigrat, obigrat;
                mpq_init(bigrat);
                mpq_init(obigrat);
                mpq_set_z(bigrat, v._bigint);
                mpq_set_z(obigrat, other.v._bigint);
                mpq_div(bigrat, bigrat, obigrat);
                mpq_clear(obigrat);
                return bigrat;
        }
        case MPQ: {
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_div(bigrat, v._bigrat, other.v._bigrat);
                return bigrat;
        }
        case PYOBJECT:
                if (PyLong_Check(v._pyobject)) {
                        PyObject *n = py_funcs.
                                py_integer_from_python_obj(v._pyobject);
                        PyObject *ans = PyNumber_TrueDivide(n,
                                        other.v._pyobject);
                        Py_DECREF(n);
                        return ans;
                }
                return PyNumber_TrueDivide(v._pyobject, other.v._pyobject);

        default:
                stub("invalid type: operator/() type not handled");
        }
}


// Compute `a^b` as an integer, where a is an integer. Assign to ``res``` if it is integral, or return ``false``.
// The nonnegative real root is taken for even denominators. To be used inside numeric::integer_rational_power,
// to handle the special case of integral ``a``.
bool integer_rational_power_of_mpz(
        numeric& res,
        const numeric& a,
        const numeric& b
) {
        if (a.t != MPZ)
                throw std::runtime_error("integer_rational_power_of_mpz: bad input");
        mpz_t z;
        mpz_init(z);
        mpz_set_ui(z, 0);
        int sgn = mpz_sgn(a.v._bigint);
        if (mpz_cmp_ui(a.v._bigint, 1) == 0
            or mpz_cmp_ui(mpq_numref(b.v._bigrat), 0) == 0)
                mpz_set_ui(z, 1);
        else if (sgn == 0) {
                res = *_num0_p;
                mpz_clear(z);
                return true;
        }
        else if (sgn < 0 and mpz_cmp_ui(mpq_denref(b.v._bigrat), 1)) {
                mpz_clear(z);
                return false;
        } else {
                if (not mpz_fits_ulong_p(mpq_numref(b.v._bigrat))
                    or not mpz_fits_ulong_p(mpq_denref(b.v._bigrat))) {
                        // too big to take roots/powers
                        mpz_clear(z);
                        return false;
                }
                if (mpz_cmp_ui(mpq_denref(b.v._bigrat), 2) == 0) {
                        if (mpz_perfect_square_p(a.v._bigint)) {
                                mpz_sqrt(z, a.v._bigint);
                        } else {
                                mpz_clear(z);
                                return false;
                        }
                }
                else {
                        bool exact = mpz_root(z, a.v._bigint,
                                              mpz_get_ui(mpq_denref(b.v._bigrat)));
                        if (not exact) {
                                mpz_clear(z);
                                return false;
                        }
                }
                mpz_pow_ui(z, z, mpz_get_ui(mpq_numref(b.v._bigrat)));
        }
        res = numeric(z);   // transfers ownership, no mpz_clear
        return true;
}


// Compute `a^b` as an integer, if it is integral, or return ``false``.
// The nonnegative real root is taken for even denominators.
bool numeric::integer_rational_power(numeric& res,
                const numeric& a, const numeric& b)
{
        if (b.t != MPQ)
                throw std::runtime_error("integer_rational_power: bad input");
        if (mpz_sgn(mpq_numref(b.v._bigrat)) < 0)
                throw std::runtime_error("integer_rational_power: bad input");
        if (a.t == LONG) {
                if (a.v._long == 1
                    or mpz_cmp_ui(mpq_numref(b.v._bigrat), 0) == 0) {
                        res = 1;
                        return true;
                }
                if (a.v._long == 0) {
                        res = 0;
                        return true;
                }
                if (a.v._long < 0
                    and mpz_cmp_ui(mpq_denref(b.v._bigrat), 1))
                        return false;
                if (not mpz_fits_ulong_p(mpq_numref(b.v._bigrat))
                    or not mpz_fits_ulong_p(mpq_denref(b.v._bigrat)))
                // too big to take roots/powers
                        return false;
                if (b.is_equal(*_num1_2_p)) {
                        long z = std::lround(std::sqrt(a.v._long));
                        if (a.v._long == z*z) {
                                res = numeric(z);
                                return true;
                        }
                        return false;
                }
                return integer_rational_power(res, a.to_bigint(), b);
        }
        // otherwise: a is integer
        return integer_rational_power_of_mpz(res, a, b);
}


// for a^b return c,d such that a^b = c*d^b
// only for MPZ/MPQ base and MPQ exponent
void rational_power_parts(const numeric& a_orig, const numeric& b_orig,
                numeric& c, numeric& d, bool& c_unit)
{
        if (b_orig.t == LONG or b_orig.t == MPZ) {
                c = a_orig.pow_intexp(b_orig);
                d = *_num1_p;
                c_unit = c.is_one();
                return;
        }

        if ((a_orig.t != LONG
             and a_orig.t != MPZ
             and a_orig.t != MPQ)
            or b_orig.t != MPQ) {
                d = a_orig;
                c = *_num1_p;
                c_unit = true;
                return;
        }
        bool b_negative = b_orig.is_negative();
        const numeric& a = b_negative? a_orig.inverse() : a_orig;
        const numeric& b = b_negative? b_orig.negative() : b_orig;
        if (a.t == MPQ) {
                numeric c1, c2, d1, d2;
                rational_power_parts(a.numer(), b, c1, d1, c_unit);
                rational_power_parts(a.denom(), b, c2, d2, c_unit);
                c = c1 / c2;
                if (b_negative)
                        d = d2 / d1;
                else
                        d = d1 / d2;
                c_unit = c.is_one();
                return;
        }

        // remaining a types are LONG and MPZ
        if (numeric::integer_rational_power(c, a, b)) {
                c_unit = c.is_one();
                d = *_num1_p;
                return;
        }
        numeric numer = b.numer();
        numeric denom = b.denom();
        if (denom.t == MPZ
            and not mpz_fits_slong_p(denom.v._bigint)) {
                c = *_num1_p;
                c_unit = true;
                d = a;
                return;
        }
        long denoml = denom.to_long();
        if (denoml > 1 and a.is_minus_one()) {
                c = *_num1_p;
                c_unit = true;
                d = *_num_1_p;
                return;
        }
        std::vector<std::pair<numeric, int>> factors;
        long max_prime_idx;
        static numeric maxnum = numeric(10).power(200);
        if (a.abs() < maxnum) {
                double ad;
                if (a.t == MPZ)
                        ad = mpz_get_d(a.v._bigint);
                else
                        ad = a.v._long;
                ad = std::abs(ad);
                double adrootlog = std::log(ad) / denoml;
                double adroot = std::exp(adrootlog);
                double mpid = adroot*1.25506/adrootlog;
                max_prime_idx = mpid>2000? 2000L : long(mpid) + 1;
        }
        else // can't convert a to double
                max_prime_idx = 2000L;
        a.factor(factors, max_prime_idx);
        c = *_num1_p;
        d = *_num1_p;
        for (auto p : factors) {
                c = c * p.first.pow_intexp(numer * numeric(p.second / 
                                        denoml));
                d = d * p.first.pow_intexp(numeric(p.second % denoml));
        }
        if (a.is_negative() and numer.is_odd())
                d = d.negative();
        if (b_negative)
                d = d.inverse();
        c_unit = c.is_one();
}

const numeric numeric::power(signed long exp_si) const
{
        PyObject *o, *r;
        if (exp_si == 0)
                return *_num1_p;
        if (exp_si == 1)
                return *this;
        switch (t) {
        case LONG:
                if (exp_si >= 0) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_set_si(bigint, v._long);
                        mpz_pow_ui(bigint, bigint, exp_si);
                        return numeric(bigint);
                } else {
                        mpz_t bigint;
                        mpz_init_set_si(bigint, v._long);
                        if (exp_si != -1)
                                mpz_pow_ui(bigint, bigint, -exp_si);
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_z(bigrat, bigint);
                        mpq_inv(bigrat, bigrat);
                        mpz_clear(bigint);
                        return numeric(bigrat);
                }
        case MPZ:
                if (exp_si >= 0) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_pow_ui(bigint, v._bigint, exp_si);
                        return numeric(bigint);
                } else {
                        mpz_t bigint;
                        mpz_init_set(bigint, v._bigint);
                        mpz_pow_ui(bigint, bigint, -exp_si);
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_z(bigrat, bigint);
                        mpq_inv(bigrat, bigrat);
                        mpz_clear(bigint);
                        return numeric(bigrat);
                }
        case MPQ:
                mpz_t bigint;
                mpq_t bigrat, obigrat;
                mpz_init(bigint);
                mpq_init(bigrat);
                mpq_init(obigrat);
                if (exp_si >= 0) {
                        mpz_pow_ui(bigint, mpq_numref(v._bigrat), exp_si);
                        mpq_set_z(bigrat, bigint);
                        mpz_pow_ui(bigint, mpq_denref(v._bigrat), exp_si);
                        mpq_set_z(obigrat, bigint);
                        mpq_div(bigrat, bigrat, obigrat);
                } else {
                        mpz_pow_ui(bigint, mpq_denref(v._bigrat), -exp_si);
                        mpq_set_z(bigrat, bigint);
                        mpz_pow_ui(bigint, mpq_numref(v._bigrat), -exp_si);
                        mpq_set_z(obigrat, bigint);
                        mpq_div(bigrat, bigrat, obigrat);
                }
                mpz_clear(bigint);
                mpq_clear(obigrat);
                return numeric(bigrat);
        case PYOBJECT:
                o = Integer(exp_si);
                r = PyNumber_Power(v._pyobject, o, Py_None);
                Py_DECREF(o);
                return numeric(r);
        default:
                stub("invalid type: pow_intexp numeric");
        }
}

const numeric numeric::pow_intexp(const numeric &exponent) const
{
        if (not exponent.is_integer())
                throw std::runtime_error("numeric::pow_intexp: exponent not integer");
	// Because Sage may have simplified a fraction to an integer,
	// we should check whether the type is MPQ as well as MPZ
        if (exponent.t == MPZ) {
                if (not mpz_fits_sint_p(exponent.v._bigint))
                        throw std::runtime_error("size of exponent exceeds signed long size");
                return power(mpz_get_si(exponent.v._bigint));
        }
	if (exponent.t == MPQ) {
		if (not mpz_fits_sint_p(mpq_numref(exponent.v._bigrat)))
			throw std::runtime_error("size of exponent exceeds signed long size");
		return power(mpz_get_si(mpq_numref(exponent.v._bigrat)));
        }
        return power(exponent.v._long);
}

/** Numerical exponentiation.
 * Raises *this to the power given as argument and returns the result as ex,
 * because it is possible that the result is converted to symbolic by Sage,
 * and we don't want symbolic as PyObject inside numeric. */
const ex numeric::power(const numeric &exponent) const {
        verbose("pow");
        numeric expo(exponent);

        // any PyObjects castable to long are casted
        if (exponent.t == PYOBJECT) {
                if (PyLong_Check(exponent.v._pyobject)) {
                        expo.t = MPZ;
                        _mpz_set_pylong(expo.v._bigint,
                        reinterpret_cast<PyLongObject *>(exponent.v._pyobject));
                }
        }

        // inexact PyObjects in base or exponent
        if (t == PYOBJECT and not is_exact()) {
                if (expo.t == PYOBJECT)
                        return numeric(PyNumber_Power(v._pyobject,
                                        exponent.v._pyobject,
                                        Py_None));
                else {
                        PyObject* obj = exponent.to_pyobject();
                        const numeric& ret = numeric(PyNumber_Power(v._pyobject,
                                        obj,
                                        Py_None));
                        Py_DECREF(obj);
                        return ret;
                }
        }
        if (expo.t == PYOBJECT and not expo.is_exact()) {
                if (t == PYOBJECT)
                        return numeric(PyNumber_Power(v._pyobject,
                                        exponent.v._pyobject,
                                        Py_None));
                else {
                        PyObject* obj = to_pyobject();
                        const numeric& ret = numeric(PyNumber_Power(obj,
                                        exponent.v._pyobject,
                                        Py_None));
                        Py_DECREF(obj);
                        return ret;
                }
        }

        // handle all integer exponents
        if (expo.t == LONG or expo.t == MPZ)
                return pow_intexp(expo);
        if (expo.t == MPQ and expo.is_integer())
                return power(exponent.to_long());

        using POW = class power;
        using MUL = class mul;
        // return any complex power unchanged
        if (not is_rational()
            or not exponent.is_rational())
                return (new POW(*this, expo))->
                        setflag(status_flags::dynallocated
                                        | status_flags::evaluated);
        // rational exponent
        if ((t == LONG or t == MPZ or t == MPQ) and expo.t == MPQ) {
                numeric c, d;
                bool c_unit;
                rational_power_parts(*this, expo, c, d, c_unit);
                if (d.is_one())
                        return c;
                if (d.is_minus_one()
                                and expo.denom().is_equal(*_num2_p)) {
                        switch (expo.numer().to_long() % 4) {
                        case 0: return c;
                        case 1: return I * c;
                        case 2: return -c;
                        case 3: return -I * c;
                        };
                }
                else if (not c_unit)
                        return d.power(expo) * c;
                if (exponent.is_negative()) {
                        long int_exp = -(expo.to_long());
                        numeric nexp = expo + numeric(int_exp);
                        ex p = (new POW(*this, nexp))->
                                setflag(status_flags::dynallocated
                                                | status_flags::evaluated);
                        return p * pow_intexp(int_exp).inverse();
                }
                if (expo < *_num1_p)
                        return (new POW(*this, expo))->
                                setflag(status_flags::dynallocated
                                                | status_flags::evaluated);
                // rational expo > 1
                const numeric m = expo.denom();
                numeric rem;
                numeric quo = expo.numer().iquo(m, rem);
                if (rem.is_negative()) {
                        rem += m;
                        quo -= *_num1_p;
                }
                if (quo.is_zero()) {
                        ex res = (new POW(d, expo))->
                                         setflag(status_flags::dynallocated
                                                | status_flags::evaluated);
                        if (c_unit)
                                return (new MUL(res, c))->
                                        setflag(status_flags::dynallocated
                                                | status_flags::evaluated);
                }
                if (rem.is_zero())
                        return this->power(quo);
                return (new MUL(this->power(rem / m), this->power(quo)))->
                                        setflag(status_flags::dynallocated
                                                | status_flags::evaluated);
        }

        // exact PyObjects in base or exponent
        if (t == PYOBJECT and exponent.t == PYOBJECT) 
                return numeric(PyNumber_Power(v._pyobject,
                                        exponent.v._pyobject,
                                        Py_None));
        if (t == PYOBJECT) {
                PyObject* obj = exponent.to_pyobject();
                const numeric& ret = numeric(PyNumber_Power(v._pyobject,
                                        obj,
                                        Py_None));
                Py_DECREF(obj);
                return ret;
        }
        if (expo.t == PYOBJECT) {
                PyObject* obj = to_pyobject();
                const numeric& ret = numeric(PyNumber_Power(to_pyobject(),
                                        obj,
                                        Py_None));
                Py_DECREF(obj);
                return ret;
        }

        throw std::runtime_error("numeric::power: can't happen");
}


/** Numerical addition method.  Adds argument to *this and returns result as
 *  a numeric object on the heap.  Use internally only for direct wrapping into
 *  an ex object, where the result would end up on the heap anyways. */
const numeric &numeric::add_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral element by pointer.  This hack
        // is supposed to keep the number of distinct numeric objects low.
        if (this == _num0_p)
                return other;
        if (&other == _num0_p)
                return *this;

        return static_cast<const numeric &> ((new numeric(*this + other))->
                setflag(status_flags::dynallocated));
}

/** Numerical subtraction method.  Subtracts argument from *this and returns
 *  result as a numeric object on the heap.  Use internally only for direct
 *  wrapping into an ex object, where the result would end up on the heap
 *  anyways. */
const numeric &numeric::sub_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral exponent (first by pointer).  This
        // hack is supposed to keep the number of distinct numeric objects low.
        if (&other == _num0_p || (other.is_zero()))
                return *this;

        return static_cast<const numeric &> ((new numeric(*this - other))->
                setflag(status_flags::dynallocated));
}

/** Numerical multiplication method.  Multiplies *this and argument and returns
 *  result as a numeric object on the heap.  Use internally only for direct
 *  wrapping into an ex object, where the result would end up on the heap
 *  anyways. */
const numeric &numeric::mul_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral element by pointer.  This hack
        // is supposed to keep the number of distinct numeric objects low.
        if (this == _num1_p)
                return other;
        if (&other == _num1_p)
                return *this;

        return static_cast<const numeric &> ((new numeric(*this * other))->
                setflag(status_flags::dynallocated));
}

/** Numerical division method.  Divides *this by argument and returns result as
 *  a numeric object on the heap.  Use internally only for direct wrapping
 *  into an ex object, where the result would end up on the heap
 *  anyways.
 *
 *  @exception overflow_error (division by zero) */
const numeric &numeric::div_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral element by pointer.  This hack
        // is supposed to keep the number of distinct numeric objects low.
        if (&other == _num1_p)
                return *this;
        if (other.is_zero())
                throw std::overflow_error("division by zero");
        return static_cast<const numeric &> ((new numeric(*this / other))->
                setflag(status_flags::dynallocated));
}

/** Numerical exponentiation.  Raises *this to the power given as argument and
 *  returns result as a numeric object on the heap.  Use internally only for
 *  direct wrapping into an ex object, where the result would end up on the
 *  heap anyways. */
// NOTE: Use only if sure to return a non-expression (see power()).
const numeric &numeric::power_dyn(const numeric &other) const {
        // Efficiency shortcut: trap the neutral exponent (first try by pointer, then
        // try harder, since calls to cln::expt() below may return amazing results for
        // floating point exponent 1.0).
        if (&other == _num1_p || (other == *_num1_p))
                return *this;

        return static_cast<const numeric &> ((new numeric(ex_to<numeric>(pow(*this, other))))->setflag(status_flags::dynallocated));
}

const numeric &numeric::operator=(int i) {
        return operator=(numeric(i));
}

const numeric &numeric::operator=(unsigned int i) {
        return operator=(numeric(i));
}

const numeric &numeric::operator=(long i) {
        return operator=(numeric(i));
}

const numeric &numeric::operator=(double d) {
        return operator=(numeric(d));
}

const numeric numeric::negative() const {
        verbose("operator-");
        switch (t) {
        case LONG:
                if (v._long != std::numeric_limits<long>::min())
                        return -v._long;
                else  // use multiplication to avoid negation overflow
                        return -1 * *this;
        case MPZ:
                mpz_t bigint;
                mpz_init_set(bigint, v._bigint);
                mpz_neg(bigint, bigint);
                return bigint;
        case MPQ:
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_set(bigrat, v._bigrat);
                mpq_neg(bigrat, bigrat);
                return bigrat;
        case PYOBJECT:
                return PyNumber_Negative(v._pyobject);
        default:
                stub("invalid type: operator-() type not handled");
        }
}

// binary arithmetic assignment operators with numeric

numeric & operator+=(numeric & lh, const numeric & rh)
{
        if (rh.is_zero())
                return lh;
        if (lh.is_zero()) {
                lh = rh;
                return lh;
        }
        if (lh.t != rh.t) {
                if (lh.t == MPZ and rh.t == MPQ) {
                        mpz_t bigint;
                        mpz_init_set(bigint, lh.v._bigint);
                        mpz_clear(lh.v._bigint);
                        lh.t = MPQ;
                        mpq_init(lh.v._bigrat);
                        mpq_set_z(lh.v._bigrat, bigint);
                        mpq_add(lh.v._bigrat, lh.v._bigrat, rh.v._bigrat);
                        lh.hash = _mpq_pythonhash(lh.v._bigrat);
                        mpz_clear(bigint);
                        return lh;
                }
                if (lh.t == MPQ and rh.t == MPZ) {
                        mpq_t tmp;
                        mpq_init(tmp);
                        mpq_set_z(tmp, rh.v._bigint);
                        mpq_add(lh.v._bigrat, lh.v._bigrat, tmp);
                        lh.hash = _mpq_pythonhash(lh.v._bigrat);
                        mpq_clear(tmp);
                        return lh;
                }

                numeric a, b;
                coerce(a, b, lh, rh);
                lh = a + b;
                return lh;
        }
        switch (lh.t) {
        case LONG: {
                if ((lh.v._long > 0
                    and lh.v._long < std::numeric_limits<long>::max() / 2
                    and rh.v._long < std::numeric_limits<long>::max() / 2)
                or (lh.v._long < 0
                   and lh.v._long > std::numeric_limits<long>::min() / 2
                   and rh.v._long > std::numeric_limits<long>::min() / 2)) {
                        lh.v._long += rh.v._long;
                        lh.hash = (lh.v._long == -1) ? -2 : lh.v._long;
                        return lh;
                }
                lh.t = MPZ;
                mpz_init_set_si(lh.v._bigint, lh.v._long);
                if (rh.v._long < 0)
                        mpz_sub_ui(lh.v._bigint,
                        lh.v._bigint, -rh.v._long);
                else
                        mpz_add_ui(lh.v._bigint, lh.v._bigint,
                        rh.v._long);
                lh.hash = _mpz_pythonhash(lh.v._bigint);
                return lh;
        }
        case MPZ:
                mpz_add(lh.v._bigint, lh.v._bigint, rh.v._bigint);
                lh.hash = _mpz_pythonhash(lh.v._bigint);
                return lh;
        case MPQ:
                mpq_add(lh.v._bigrat, lh.v._bigrat, rh.v._bigrat);
                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                return lh;
        case PYOBJECT: {
                PyObject *p = lh.v._pyobject;
                lh.v._pyobject = PyNumber_Add(p, rh.v._pyobject);
                if (lh.v._pyobject == nullptr) {
                        lh.v._pyobject = p;
                        py_error("numeric operator+=");
                }
                lh.hash = PyObject_Hash(lh.v._pyobject);
                Py_DECREF(p);
                return lh;
        }
        default:
                stub("invalid type: operator+=() type not handled");
        }
}

numeric & operator-=(numeric & lh, const numeric & rh)
{
        if (rh.is_zero())
                return lh;
        if (lh.is_zero()) {
                lh = rh.negative();
                return lh;
        }
        if (lh.t != rh.t) {
                if (lh.t == MPZ and rh.t == MPQ) {
                        mpz_t bigint;
                        mpz_init_set(bigint, lh.v._bigint);
                        mpz_clear(lh.v._bigint);
                        lh.t = MPQ;
                        mpq_init(lh.v._bigrat);
                        mpq_set_z(lh.v._bigrat, bigint);
                        mpq_sub(lh.v._bigrat, lh.v._bigrat, rh.v._bigrat);
                        lh.hash = _mpq_pythonhash(lh.v._bigrat);
                        mpz_clear(bigint);
                        return lh;
                }
                if (lh.t == MPQ and rh.t == MPZ) {
                        mpq_t tmp;
                        mpq_init(tmp);
                        mpq_set_z(tmp, rh.v._bigint);
                        mpq_sub(lh.v._bigrat, lh.v._bigrat, tmp);
                        lh.hash = _mpq_pythonhash(lh.v._bigrat);
                        mpq_clear(tmp);
                        return lh;
                }

                numeric a, b;
                coerce(a, b, lh, rh);
                lh = a - b;
                return lh;
        }
        switch (lh.t) {
        case LONG: {
                if ((lh.v._long > 0
                    and lh.v._long < std::numeric_limits<long>::max() / 2
                    and -rh.v._long < std::numeric_limits<long>::max() / 2)
                or (lh.v._long < 0
                   and lh.v._long > std::numeric_limits<long>::min() / 2
                   and -rh.v._long > std::numeric_limits<long>::min() / 2)) {
                        lh.v._long -= rh.v._long;
                        lh.hash = (lh.v._long == -1) ? -2 : lh.v._long;
                        return lh;
                }
                lh.t = MPZ;
                mpz_init_set_si(lh.v._bigint, lh.v._long);
                if (rh.v._long < 0)
                        mpz_add_ui(lh.v._bigint,
                        lh.v._bigint, -rh.v._long);
                else
                        mpz_sub_ui(lh.v._bigint, lh.v._bigint,
                        rh.v._long);
                lh.hash = _mpz_pythonhash(lh.v._bigint);
                return lh;
        }
        case MPZ:
                mpz_sub(lh.v._bigint, lh.v._bigint, rh.v._bigint);
                lh.hash = _mpz_pythonhash(lh.v._bigint);
                return lh;
        case MPQ:
                mpq_sub(lh.v._bigrat, lh.v._bigrat, rh.v._bigrat);
                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                return lh;
        case PYOBJECT: {
                PyObject *p = lh.v._pyobject;
                lh.v._pyobject = PyNumber_Subtract(p, rh.v._pyobject);
                if (lh.v._pyobject == nullptr) {
                        lh.v._pyobject = p;
                        py_error("numeric operator-=");
                }
                lh.hash = PyObject_Hash(lh.v._pyobject);
                Py_DECREF(p);
                return lh;
        }
        default:
                stub("invalid type: operator-() type not handled");
        }
}

numeric & operator*=(numeric & lh, const numeric & rh)
{
        if (rh.is_one())
                return lh;
        if (lh.is_one()) {
                lh = rh;
                return lh;
        }
        if ((lh.is_zero() and lh.t != PYOBJECT)
            or (rh.is_zero() and rh.t != PYOBJECT)) {
                lh = *_num0_p;
                return lh;
        }
        if (lh.t != rh.t) {
                if (lh.t == MPZ and rh.t == MPQ) {
                        mpq_t tmp;
                        mpq_init(tmp);
                        mpq_set_z(tmp, lh.v._bigint);
                        mpq_mul(tmp, tmp, rh.v._bigrat);
                        if (mpz_cmp_ui(mpq_denref(tmp),1) == 0) {
                                mpz_set(lh.v._bigint, mpq_numref(tmp));
                                lh.hash = _mpz_pythonhash(lh.v._bigint);
                                mpq_clear(tmp);
                                return lh;
                        }
                        mpz_clear(lh.v._bigint);
                        lh.t = MPQ;
                        mpq_init(lh.v._bigrat);
                        mpq_set(lh.v._bigrat, tmp);
                        lh.hash = _mpq_pythonhash(lh.v._bigrat);
                        mpq_clear(tmp);
                        return lh;
                }
                if (lh.t == MPQ and rh.t == MPZ) {
                        mpq_t tmp;
                        mpq_init(tmp);
                        mpq_set_z(tmp, rh.v._bigint);
                        mpq_mul(tmp, tmp, lh.v._bigrat);
                        if (mpz_cmp_ui(mpq_denref(tmp),1) != 0) {
                                mpq_set(lh.v._bigrat, tmp);
                                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                                mpq_clear(tmp);
                                return lh;
                        }
                        mpq_clear(lh.v._bigrat);
                        lh.t = MPZ;
                        mpz_init(lh.v._bigint);
                        mpz_set(lh.v._bigint, mpq_numref(tmp));
                        lh.hash = _mpz_pythonhash(lh.v._bigint);
                        mpq_clear(tmp);
                        return lh;
                }

                numeric a, b;
                coerce(a, b, lh, rh);
                lh = a * b;
                return lh;
        }
        switch (lh.t) {
        case LONG: {
                long result;
                if (!smull_overflow(lh.v._long, rh.v._long, & result)) {
                        lh.v._long = result;
                        lh.hash = (lh.v._long==-1) ? -2 : lh.v._long;
                        return lh;
                }
                // the multiplication overflowed, so use mpz
                lh.t = MPZ;
                mpz_init_set_si(lh.v._bigint, lh.v._long);
                mpz_mul_si(lh.v._bigint, lh.v._bigint, rh.v._long);
                lh.hash = _mpz_pythonhash(lh.v._bigint);
                return lh;
        }
        case MPZ:
                mpz_mul(lh.v._bigint, lh.v._bigint, rh.v._bigint);
                lh.hash = _mpz_pythonhash(lh.v._bigint);
                return lh;
        case MPQ:
                mpq_mul(lh.v._bigrat, lh.v._bigrat, rh.v._bigrat);
                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                return lh;
        case PYOBJECT: {
                PyObject *p = lh.v._pyobject;
                lh.v._pyobject = PyNumber_Multiply(p, rh.v._pyobject);
                if (lh.v._pyobject == nullptr) {
                        lh.v._pyobject = p;
                        py_error("numeric operator*=");
                }
                lh.hash = PyObject_Hash(lh.v._pyobject);
                Py_DECREF(p);
                return lh;
        }
        default:
                stub("invalid type: operator*=() type not handled");
        }
}

template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}
numeric & operator/=(numeric & lh, const numeric & rh)
{
        if (rh.is_zero())
                throw std::overflow_error("numeric::/=(): division by zero");
        if (rh.is_one())
                return lh;
        if (lh.t != rh.t) {
                if (lh.t == MPZ and rh.t == MPQ) {
                        mpq_t tmp;
                        mpq_init(tmp);
                        mpq_set_z(tmp, lh.v._bigint);
                        mpq_div(tmp, tmp, rh.v._bigrat);
                        if (mpz_cmp_ui(mpq_denref(tmp),1) == 0) {
                                mpz_set(lh.v._bigint, mpq_numref(tmp));
                                lh.hash = _mpz_pythonhash(lh.v._bigint);
                                mpq_clear(tmp);
                                return lh;
                        }
                        mpz_clear(lh.v._bigint);
                        lh.t = MPQ;
                        mpq_init(lh.v._bigrat);
                        mpq_set(lh.v._bigrat, tmp);
                        lh.hash = _mpq_pythonhash(lh.v._bigrat);
                        mpq_clear(tmp);
                        return lh;
                }
                if (lh.t == MPQ and rh.t == MPZ) {
                        mpq_t tmp;
                        mpq_init(tmp);
                        mpq_set_z(tmp, rh.v._bigint);
                        mpq_div(tmp, lh.v._bigrat, tmp);
                        if (mpz_cmp_ui(mpq_denref(tmp),1) != 0) {
                                mpq_set(lh.v._bigrat, tmp);
                                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                                mpq_clear(tmp);
                                return lh;
                        }
                        mpq_clear(lh.v._bigrat);
                        lh.t = MPZ;
                        mpz_init(lh.v._bigint);
                        mpz_set(lh.v._bigint, mpq_numref(tmp));
                        lh.hash = _mpz_pythonhash(lh.v._bigint);
                        mpq_clear(tmp);
                        return lh;
                }

                numeric a, b;
                coerce(a, b, lh, rh);
                lh = a / b;
                return lh;
        }
        switch (lh.t) {
        case LONG: {
                if (rh.v._long == -1) { // use multiplication to avoid possible overflow
                        lh *= -1;
                        return lh;
                }
                auto ld = std::div(lh.v._long, rh.v._long);
                if (ld.rem == 0) {
                        lh.v._long = ld.quot;
                        lh.hash = (ld.quot == -1) ? -2 : ld.quot;
                        return lh;
                }

                mpq_t obigrat;
                lh.t = MPQ;
                unsigned long l = std::labs(lh.v._long);
                unsigned long r = std::labs(rh.v._long);
                int sign = sgn(lh.v._long) * sgn(rh.v._long);
                mpq_init(obigrat);
                mpq_init(lh.v._bigrat);
                mpq_set_ui(lh.v._bigrat, l, 1);
                mpq_set_ui(obigrat, r, 1);
                mpq_div(lh.v._bigrat, lh.v._bigrat, obigrat);
                if (sign == -1)
                        mpq_neg(lh.v._bigrat, lh.v._bigrat);
                mpq_clear(obigrat);
                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                return lh;
        }
        case MPZ: {
                if (mpz_divisible_p(lh.v._bigint, rh.v._bigint)) {
                        mpz_divexact(lh.v._bigint,
                        lh.v._bigint,
                        rh.v._bigint);
                        lh.hash = _mpz_pythonhash(lh.v._bigint);
                        return lh;
                }
                mpq_t bigrat, obigrat;
                mpq_init(bigrat);
                mpq_init(obigrat);
                mpq_set_z(bigrat, lh.v._bigint);
                mpq_set_z(obigrat, rh.v._bigint);
                mpz_clear(lh.v._bigint);
                lh.t = MPQ;
                mpq_init(lh.v._bigrat);
                mpq_div(lh.v._bigrat, bigrat, obigrat);
                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                mpq_clear(bigrat);
                mpq_clear(obigrat);
                return lh;
        }
        case MPQ:
                mpq_div(lh.v._bigrat, lh.v._bigrat, rh.v._bigrat);
                lh.hash = _mpq_pythonhash(lh.v._bigrat);
                return lh;
        case PYOBJECT: {
                PyObject *p = lh.v._pyobject;
                {
                        if (PyLong_Check(p)) {
                                PyObject *n = py_funcs.py_integer_from_python_obj(p);
                                lh.v._pyobject = PyNumber_TrueDivide(n, rh.v._pyobject);
                                if (lh.v._pyobject == nullptr) {
                                        lh.v._pyobject = p;
                                        py_error("numeric operator/=");
                                }
                                lh.hash = (long)PyObject_Hash(lh.v._pyobject);
                                Py_DECREF(n);
                                Py_DECREF(p);
                                return lh;
                        }
                }

                lh.v._pyobject = PyNumber_TrueDivide(p, rh.v._pyobject);
                if (lh.v._pyobject == nullptr) {
                        lh.v._pyobject = p;
                        py_error("numeric operator/=");
                }
                lh.hash = PyObject_Hash(lh.v._pyobject);
                Py_DECREF(p);
                return lh;
        }
        default:
                stub("invalid type: operator/=() type not handled");
        }
}


/** Inverse of a number. */
const numeric numeric::inverse() const {
        if (is_zero())
                throw std::overflow_error("numeric::inverse(): division by zero");
        return numeric(1) / *this;
}

/** Return the step function of a numeric. The imaginary part of it is
 *  ignored because the step function is generally considered real but
 *  a numeric may develop a small imaginary part due to rounding errors.
 */
const numeric numeric::step() const {
        switch (t) {
        case LONG: return v._long > 0;
        case MPZ:
        case MPQ:
                if (is_positive())
                        return 1;
                else
                        return 0;
        case PYOBJECT:
                return py_funcs.py_step(v._pyobject);
        default:
                stub("invalid type: step() type not handled");
        }
}

/** Return the complex half-plane (left or right) in which the number lies.
 *  csgn(x)==0 for x==0, csgn(x)==1 for Re(x)>0 or Re(x)=0 and Im(x)>0,
 *  csgn(x)==-1 for Re(x)<0 or Re(x)=0 and Im(x)<0.
 *
 *  @see numeric::compare(const numeric &other) */
int numeric::csgn() const {
        verbose("csgn");
        switch (t) {
        case LONG:
                if (v._long == 0)
                        return 0;
                return v._long < 0 ? -1 : 1;
        case MPZ:
                return mpz_sgn(v._bigint);
        case MPQ:
                return mpq_sgn(v._bigrat);
        case PYOBJECT: {
                int result;
                if (is_real()) {
                        numeric z(ZERO);
                        Py_INCREF(ZERO);
                        return compare_same_type(z);
                }
                numeric re = real();
                numeric z(ZERO);
                Py_INCREF(ZERO);
                result = re.compare_same_type(z);
                if (result == 0) {
                        numeric im = imag();
                        result = im.compare_same_type(z);
                }
                return result;
        }
        default:
                stub("invalid type: csgn() type not handled");
        }
}

bool numeric::is_equal(const numeric &other) const {
        return *this == other;
}

/** True if object is zero. */
bool numeric::is_zero() const {
        verbose("is_zero");
        int a;
        switch (t) {
        case LONG:
                return v._long == 0;
        case MPZ:
                return mpz_cmp_si(v._bigint, 0) == 0;
        case MPQ:
                return mpq_cmp_si(v._bigrat, 0, 1) == 0;
        case PYOBJECT:
                a = PyObject_Not(v._pyobject);
                if (a == -1)
                        py_error("is_zero");
                return a == 1;
        default:
                std::cerr << "type = " << t << "\n";
                stub("invalid type: is_zero() type not handled");
        }
}

bool numeric::is_inexact_one() const {
        verbose("is_one");
        switch (t) {
        case LONG:
                return v._long == 1;
        case MPZ:
                return mpz_cmp_si(v._bigint, 1) == 0;
        case MPQ:
                return mpq_cmp_si(v._bigrat, 1, 1) == 0;
        case PYOBJECT:
                return is_equal(*_num1_p);
        default:
                std::cerr << "type = " << t << "\n";
                stub("invalid type: is_one() type not handled");
        }
}

bool numeric::is_one() const {
        verbose("is_one");
        switch (t) {
        case LONG:
                return v._long == 1;
        case MPZ:
                return mpz_cmp_si(v._bigint, 1) == 0;
        case MPQ:
                return mpq_cmp_si(v._bigrat, 1, 1) == 0;
        case PYOBJECT:
                return is_exact() and is_equal(*_num1_p);
        default:
                std::cerr << "type = " << t << "\n";
                stub("invalid type: is_one() type not handled");
        }
}

bool numeric::is_minus_one() const {
        verbose("is_one");
        switch (t) {
        case LONG:
                return v._long == -1;
        case MPZ:
                return mpz_cmp_si(v._bigint, -1) == 0;
        case MPQ:
                return mpq_cmp_si(v._bigrat, -1, 1) == 0;
        case PYOBJECT:
                return is_exact() and is_equal(*_num_1_p);
        default:
                std::cerr << "type = " << t << "\n";
                stub("invalid type: is_minus_one() type not handled");
        }
}

/** True if object is not complex and greater than zero. */
bool numeric::is_positive() const {
        verbose("is_positive");
        switch (t) {
        case LONG:
                return v._long > 0;
        case MPZ:
                return mpz_cmp_si(v._bigint, 0) > 0;
        case MPQ:
                return mpq_cmp_si(v._bigrat, 0, 1) > 0;
        case PYOBJECT:
                if (is_real()) {
                        int result;
                        result = PyObject_RichCompareBool(v._pyobject,
                                        ZERO, Py_GT);
                        if (result == 1)
                                return true;
                        if (result == -1)
                                PyErr_Clear();
                }
                return false;
        default:
                stub("invalid type: is_positive() type not handled");
        }
}

/** True if object is not complex and less than zero. */
bool numeric::is_negative() const {
        verbose("is_negative");
        switch (t) {
        case LONG:
                return v._long < 0;
        case MPZ:
                return mpz_cmp_si(v._bigint, 0) < 0;
        case MPQ:
                return mpq_cmp_si(v._bigrat, 0, 1) < 0;
        case PYOBJECT:
                if (is_real()) {
                        int result;
                        result = PyObject_RichCompareBool(v._pyobject,
                                        ZERO, Py_LT);
                        if (result == 1)
                                return true;
                        if (result == -1)
                                PyErr_Clear();
                }
                return false;
        default:
                stub("invalid type: is_negative() type not handled");
        }
}

/** True if object is a non-complex integer. */
bool numeric::is_integer() const {
        verbose2("is_integer", *this);

        bool ret;
        switch (t) {
        case LONG:
        case MPZ:
                return true;
        case MPQ: {
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_set(bigrat, v._bigrat);
                mpq_canonicalize(bigrat);
                ret = mpz_cmp_ui(mpq_denref(bigrat), 1) == 0;
                mpq_clear(bigrat);
                return ret;
        }
        case PYOBJECT:
                return py_funcs.py_is_integer(v._pyobject) != 0;
        default:
                stub("invalid type: is_integer() type not handled");
        }
}

/** True if object is an exact integer greater than zero. */
bool numeric::is_pos_integer() const {
        verbose("is_pos_integer");
        switch (t) {
        case LONG:
                return v._long > 0;
        case MPZ:
                return is_positive();
        case MPQ:
                return (is_integer() && is_positive());
        case PYOBJECT:
                return (is_integer() && is_positive());
        default:
                stub("invalid type: is_pos_integer() type not handled");
        }
}

/** True if object is an exact integer greater or equal zero. */
bool numeric::is_nonneg_integer() const {
        verbose("is_nonneg_integer");
        switch (t) {
        case LONG:
                return v._long >= 0;
        case MPZ:
                return is_positive() or is_zero();
        case MPQ:
                return (is_integer() and (is_positive() or is_zero()));
        case PYOBJECT:
                if (is_integer()) {
                        int result;
                        result = PyObject_RichCompareBool(v._pyobject,
                                        ZERO, Py_GE);
                        if (result == 1)
                                return true;
                        if (result == -1)
                                PyErr_Clear();
                }
                return false;
        default:
                stub("invalid type: is_nonneg_integer() type not handled");
        }
}

/** True if object is an exact even integer. */
bool numeric::is_even() const {
        verbose("is_even");

        if (!is_integer())
                return false;

        switch (t) {
        case LONG:
                return v._long % 2 == 0;
        case MPZ:
                return mpz_tstbit(v._bigint, 0) == 0;
        case MPQ:
                return is_integer()
                and mpz_tstbit(mpq_numref(v._bigrat), 0) == 0;
        case PYOBJECT:
                return py_funcs.py_is_even(v._pyobject) != 0;
        default:
                stub("invalid type: is_even() type not handled");
        }
}

/** True if object is an exact odd integer. */
bool numeric::is_odd() const {
        switch (t) {
        case LONG:
                return v._long & 1;
        case MPZ:
                return mpz_tstbit(v._bigint, 0) == 1;
        case MPQ:
                return is_integer()
                and mpz_tstbit(mpq_numref(v._bigrat), 0) == 1;
        case PYOBJECT:
                return !is_even();
        default:
                stub("invalid type: is_odd() type not handled");
        }
}

/** Probabilistic primality test.
 *
 *  @return  true if object is exact integer and prime. */
bool numeric::is_prime() const {
        verbose("is_prime");
        switch (t) {
        case LONG: {
                mpz_t bigint;
                mpz_init_set_si(bigint, v._long);
                bool ret = mpz_probab_prime_p(bigint, 25) > 0;
                mpz_clear(bigint);
                return ret;
        }
        case MPZ:
                return mpz_probab_prime_p(v._bigint, 25) > 0;
        case MPQ:
                return is_integer()
                        and mpz_probab_prime_p(mpq_numref(v._bigrat), 25) > 0;
        case PYOBJECT:
                return py_funcs.py_is_prime(v._pyobject) != 0;
        default:
                stub("invalid type: is_prime() type not handled");
        }
}

/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_rational() const {
        verbose("is_rational");
        switch (t) {
        case LONG:
        case MPZ:
        case MPQ:
                return true;
        case PYOBJECT:
                return false;
        default:
                stub("invalid type -- is_rational() type not handled");
        }
}

/** True if object is a real integer, rational or float (but not complex). */
bool numeric::is_real() const {
        verbose("is_real");
        switch (t) {
        case LONG:
        case MPZ:
        case MPQ:
                return true;
        case PYOBJECT:
                return py_funcs.py_is_real(v._pyobject) != 0;
        default:
                stub("invalid type -- is_real() type not handled");
        }
}

/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_exact() const {
        verbose("is_exact");
        switch (t) {
        case LONG:
        case MPZ:
        case MPQ:
                return true;
        case PYOBJECT:
                return py_funcs.py_is_exact(v._pyobject) != 0;
        default:
                stub("invalid type -- is_exact() type not handled");
        }
}

static std::map<long,std::pair<int,int>> small_powers;

static void fill_small_powers()
{
        static int lim[] = {30, 18, 15, 12, 11, 10, 10, 9, 9};
        for (size_t b = sizeof(lim)/sizeof(int) + 1; b >= 2; --b) {
                long p = b*b;
                int c = 2;
                while (c <= lim[b-2]) {
                        small_powers[p] = std::make_pair(int(b), c);
                        p *= b;
                        ++c;
                }
        }
}

bool numeric::is_small_power(std::pair<int,int>& p) const
{
        int i;
        switch (t) {
        case LONG:
                if (v._long < 2)
                        return false;
                i = v._long;
                break;
        case MPZ:
                if (not mpz_fits_sint_p(v._bigint))
                        return false;
                i = mpz_get_si(v._bigint);
                if (i < 2)
                        return false;
                break;
        case MPQ:
        case PYOBJECT:
                return false;
        default:
                stub("invalid type -- is_small_power() type not handled");
        }
        if (small_powers.empty())
                fill_small_powers();
        auto it = small_powers.find(i);
        if (it == small_powers.end())
                return false;
        p = it->second;
        return true;
} 

bool numeric::operator==(const numeric &right) const {
        verbose3("operator==", *this, right);

        if (this == &right)
                return true;
        if (t != right.t) {
                if (t == LONG and right.t == MPZ)
                        return mpz_cmp_si(right.v._bigint, v._long) ==0;
                if (right.t == LONG and t == MPZ)
                        return mpz_cmp_si(v._bigint, right.v._long) ==0;
                if (t == MPZ and right.t == MPQ) {
                        if (mpz_cmp_ui(mpq_denref(right.v._bigrat), 1) != 0)
                                return false;
                        return mpz_cmp(v._bigint, mpq_numref(right.v._bigrat)) ==0;
                }
                if (t == MPQ and right.t == MPZ) {
                        if (mpz_cmp_ui(mpq_denref(v._bigrat), 1) != 0)
                                return false;
                        return mpz_cmp(right.v._bigint, mpq_numref(v._bigrat)) ==0;
                }
                numeric a, b;
                coerce(a, b, *this, right);
                return a == b;
        }
        switch (t) {
        case LONG:
                return v._long == right.v._long;
        case MPZ:
                return mpz_cmp(v._bigint, right.v._bigint) == 0;
        case MPQ:
                return mpq_equal(v._bigrat, right.v._bigrat) != 0;
        case PYOBJECT:
                if (v._pyobject == right.v._pyobject)
                        return true;
                return py_funcs.py_is_equal(v._pyobject,
                                right.v._pyobject) != 0;
        default:
                stub("invalid type: operator== type not handled");
        }
}

bool numeric::operator!=(const numeric &right) const {
        verbose("operator!=");
        if (t != right.t) {
                if (t == LONG and right.t == MPZ)
                        return mpz_cmp_si(right.v._bigint, v._long) !=0;
                if (right.t == LONG and t == MPZ)
                        return mpz_cmp_si(v._bigint, right.v._long) !=0;
                if (t == MPZ and right.t == MPQ) {
                        if (mpz_cmp_ui(mpq_denref(right.v._bigrat), 1) != 0)
                                return true;
                        return mpz_cmp(v._bigint, mpq_numref(right.v._bigrat)) !=0;
                }
                if (t == MPQ and right.t == MPZ) {
                        if (mpz_cmp_ui(mpq_denref(v._bigrat), 1) != 0)
                                return true;
                        return mpz_cmp(right.v._bigint, mpq_numref(v._bigrat)) !=0;
                }
                numeric a, b;
                coerce(a, b, *this, right);
                return a != b;
        }
        switch (t) {
        case LONG:
                return v._long != right.v._long;
        case MPZ:
                return mpz_cmp(v._bigint, right.v._bigint) != 0;
        case MPQ:
                return mpq_equal(v._bigrat, right.v._bigrat) == 0;
        case PYOBJECT:
                return (py_funcs.py_is_equal(v._pyobject,
                                        right.v._pyobject) == 0);
        default:
                stub("invalid type: operator!= type not handled");
        }
}

/** True if object is element of the domain of integers extended by I, i.e. is
 *  of the form a+b*I, where a and b are integers. */
bool numeric::is_cinteger() const {
        verbose("is_crational");
        switch (t) {
        case LONG:
        case MPZ:
                return true;
        case MPQ:
                return is_integer();
        case PYOBJECT:
                return real().is_integer()
                and imag().is_integer();
        default:
                stub("invalid type -- is_cinteger() type not handled");
        }
}

/** True if object is an exact rational number, may even be complex
 *  (denominator may be unity). */
bool numeric::is_crational() const {
        verbose("is_crational");
        switch (t) {
        case LONG:
        case MPZ:
        case MPQ:
                return true;
        case PYOBJECT:
                return real().is_rational()
                and imag().is_rational();
        default:
                stub("invalid type -- is_crational() type not handled");
        }
}

/** Numerical comparison: less.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator<(const numeric &right) const {
        verbose("operator<");
        if (t == MPZ and right.t == LONG)
                return mpz_cmp_si(v._bigint, right.v._long) < 0;
        if (t == LONG and right.t == MPZ)
                return mpz_cmp_si(right.v._bigint, v._long) > 0;
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a < b;
        }

        switch (t) {
        case LONG:
                return v._long < right.v._long;
        case MPZ:
                return mpz_cmp(v._bigint, right.v._bigint) < 0;
        case MPQ:
                return mpq_cmp(v._bigrat, right.v._bigrat) < 0;
        case PYOBJECT: {
                int result;
                result = PyObject_RichCompareBool(v._pyobject,
                                right.v._pyobject, Py_LT);
                if (result == -1)
                        py_error("richcmp failed");

                return (result == 1);
        }
        default:
                stub("invalid type: operator< type not handled");
        }
}

/** Numerical comparison: less or equal.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator<=(const numeric &right) const {
        verbose("operator<=");
        if (t == MPZ and right.t == LONG)
                return mpz_cmp_si(v._bigint, right.v._long) <= 0;
        if (t == LONG and right.t == MPZ)
                return mpz_cmp_si(right.v._bigint, v._long) >= 0;
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a <= b;
        }
        switch (t) {
        case LONG:
                return v._long <= right.v._long;
        case MPZ:
                return mpz_cmp(v._bigint, right.v._bigint) <= 0;
        case MPQ:
                return mpq_cmp(v._bigrat, right.v._bigrat) <= 0;
        case PYOBJECT:
                int result;
                result = PyObject_RichCompareBool(v._pyobject,
                                right.v._pyobject, Py_LE);
                if (result == -1)
                        py_error("richcmp failed");

                return (result == 1);
        default:
                stub("invalid type: operator<= type not handled");
        }
}

/** Numerical comparison: greater.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator>(const numeric &right) const {
        verbose("operator>");
        if (t == MPZ and right.t == LONG)
                return mpz_cmp_si(v._bigint, right.v._long) > 0;
        if (t == LONG and right.t == MPZ)
                return mpz_cmp_si(right.v._bigint, v._long) < 0;
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a > b;
        }
        switch (t) {
        case LONG:
                return v._long > right.v._long;
        case MPZ:
                return mpz_cmp(v._bigint, right.v._bigint) > 0;
        case MPQ:
                return mpq_cmp(v._bigrat, right.v._bigrat) > 0;
        case PYOBJECT:
                int result;
                result = PyObject_RichCompareBool(v._pyobject,
                                right.v._pyobject, Py_GT);
                if (result == -1)
                        py_error("richcmp failed");

                return (result == 1);
        default:
                stub("invalid type: operator> type not handled");
        }
}

/** Numerical comparison: greater or equal.
 *
 *  @exception invalid_argument (complex inequality) */
bool numeric::operator>=(const numeric &right) const {
        verbose("operator>=");
        if (t == MPZ and right.t == LONG)
                return mpz_cmp_si(v._bigint, right.v._long) >= 0;
        if (t == LONG and right.t == MPZ)
                return mpz_cmp_si(right.v._bigint, v._long) <= 0;
        if (t != right.t) {
                numeric a, b;
                coerce(a, b, *this, right);
                return a >= b;
        }
        switch (t) {
        case LONG:
                return v._long >= right.v._long;
        case MPZ:
                return mpz_cmp(v._bigint, right.v._bigint) >= 0;
        case MPQ:
                return mpq_cmp(v._bigrat, right.v._bigrat) >= 0;
        case PYOBJECT:
                int result;
                result = PyObject_RichCompareBool(v._pyobject,
                                right.v._pyobject, Py_GE);
                if (result == -1)
                        py_error("richcmp failed");

                return (result == 1);
        default:
                stub("invalid type: operator!= type not handled");
        }
}

int numeric::to_int() const
{
        switch (t) {
        case LONG:
                if (v._long < std::numeric_limits<int>::max()
                and v._long > std::numeric_limits<int>::min())
                        return v._long;
                throw std::runtime_error("to_int");
        case MPZ:
                if (mpz_fits_sint_p(v._bigint))
                        return mpz_get_si(v._bigint);
                throw conversion_error();
        case MPQ: {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_fdiv_q(bigint, mpq_numref(v._bigrat),
                                mpq_denref(v._bigrat));
                if (mpz_fits_sint_p(bigint)) {
                        int n = mpz_get_si(bigint);
                        mpz_clear(bigint);
                        return n;
                }
                mpz_clear(bigint);
                throw conversion_error();
        }
        case PYOBJECT:
                return to_bigint().to_int();
        default:
                stub("invalid type: operator long int() type not handled");
        }
}
/** Converts numeric types to machine's long.  You should check with
 *  is_integer() if the number is really an integer before calling this method.
 *  You may also consider checking the range first. */
long numeric::to_long() const {
        verbose("operator long int");
        switch (t) {
        case LONG:
                return v._long;
        case MPZ:
                if (mpz_fits_slong_p(v._bigint))
                        return (long int)mpz_get_si(v._bigint);
                throw conversion_error();
        case MPQ: {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_fdiv_q(bigint, mpq_numref(v._bigrat),
                                mpq_denref(v._bigrat));
                if (mpz_fits_slong_p(v._bigint)) {
                        long n = mpz_get_si(bigint);
                        mpz_clear(bigint);
                        return n;
                }
                mpz_clear(bigint);
                throw conversion_error();
        }
        case PYOBJECT:
                return to_bigint().to_long();
        default:
                stub("invalid type: operator long int() type not handled");
        }
}

// Use this only if o was tested integer.
const numeric numeric::to_bigint() const {
        switch (t) {
        case LONG: {
                numeric nu;
                mpz_init_set_si(nu.v._bigint, v._long);
                nu.t = MPZ;
                nu.hash = _mpz_pythonhash(nu.v._bigint);
                return nu;
        }
        case MPZ: return *this;
        case MPQ:
                if (not denom().is_one())
                        throw std::runtime_error("not integer in numeric::to_mpz_num()");
                return numer();
        case PYOBJECT: {
                PyObject *Integer = Integer_pyclass();
                PyObject *ans = PyObject_CallFunctionObjArgs(Integer,
                                v._pyobject, NULL);
                return ans;
        }
        default:
                stub("invalid type: operator long int() type not handled");
        }
}

const mpz_t& numeric::as_mpz() const
{
        if (t != MPZ)
                throw std::runtime_error("mpz_t requested from non-mpz numeric");
        return v._bigint;
}

const mpq_t& numeric::as_mpq() const
{
        if (t != MPQ)
                throw std::runtime_error("mpq_t requested from non-mpq numeric");
        return v._bigrat;
}

void numeric::canonicalize()
{
        if (t == MPQ) {
                mpq_canonicalize(v._bigrat);
                if (mpz_cmp_ui(mpq_denref(v._bigrat), 1) == 0) {
                        mpz_t tmp;
                        mpz_init_set(tmp, mpq_numref(v._bigrat));
                        mpq_clear(v._bigrat);
                        set_from(t, v, hash, tmp);
                        mpz_clear(tmp);
                }
        }
}

CanonicalForm numeric::to_canonical() const
{
        if (t == LONG)
                return CanonicalForm(v._long);
        if (t == MPZ) {
                if (mpz_fits_sint_p(v._bigint))
                        return CanonicalForm(to_int());
                mpz_t bigint;
                mpz_init_set(bigint, v._bigint);
                return make_cf(mpz_ptr(const_cast<__mpz_struct*>(bigint)));
        }
        if (t == MPQ) {
                mpz_t num;
                mpz_init_set(num, mpz_ptr(mpq_numref(v._bigrat)));
                if (is_integer())
                        return make_cf(mpz_ptr(const_cast<__mpz_struct*>(num)));
                mpz_t den;
                mpz_init_set(den, mpz_ptr(mpq_denref(v._bigrat)));
                return make_cf(num, den, false);
        }

        throw (std::runtime_error("can't happen in numeric::to_canonical"));
}

/* Return the underlying Python object corresponding to this
   numeric.  If this numeric isn't implemented using a Python
   object, the corresponding Python object is constructed on
   the fly.

   Returns a NEW REFERENCE.
 */
PyObject* numeric::to_pyobject() const {
        // Returns a New Reference
        PyObject* o;
        switch (t) {
        case LONG: {
                mpz_t bigint;
                mpz_init_set_si(bigint, v._long);
                o = py_funcs.py_integer_from_mpz(bigint);
                mpz_clear(bigint);
                return o;
        }
        case MPZ: {
                mpz_t bigint;
                mpz_init_set(bigint, v._bigint);
                o = py_funcs.py_integer_from_mpz(bigint);
                mpz_clear(bigint);
                return o;
        }
        case MPQ: {
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_set(bigrat, v._bigrat);
                mpq_canonicalize(bigrat);
                o = py_funcs.py_rational_from_mpq(bigrat);
                mpq_clear(bigrat);
                return o;
        }
        case PYOBJECT:
                Py_INCREF(v._pyobject);
                return v._pyobject;

        default:
                std::cout << t << std::endl;
                stub("numeric::to_pyobject -- not able to do conversion to pyobject; everything else will be nonsense");
                return nullptr;
        }
}


/** Converts numeric types to machine's double. You should check with is_real()
 *  if the number is really not complex before calling this method. */
double numeric::to_double() const {
        GINAC_ASSERT(this->is_real());
        verbose("operator double");
        double d;
        switch (t) {
        case LONG:
                return v._long;
        case MPZ:
                return mpz_get_d(v._bigint);
        case MPQ:
                return mpq_get_d(v._bigrat);
        case PYOBJECT:
                d = PyFloat_AsDouble(v._pyobject);
                if (d == -1 && (PyErr_Occurred() != nullptr))
                        py_error("Error converting to a double.");
                return d;
        default:
                std::cerr << "type = " << t << std::endl;
                stub("invalid type: operator double() type not handled");
        }
}

/** Evaluation of numbers doesn't do anything at all. */
ex numeric::eval(int /*level*/) const {
        // Warning: if this is ever gonna do something, the ex constructors from all kinds
        // of numbers should be checking for status_flags::evaluated.
        return this->hold();
}

/** Cast numeric into a floating-point object.  For example exact numeric(1) is
 *  returned as a 1.0000000000000000000000 and so on according to how Digits is
 *  currently set.  In case the object already was a floating point number the
 *  precision is trimmed to match the currently set default.
 *
 *  @param level  ignored, only needed for overriding basic::evalf.
 *  @return  an ex-handle to a numeric. */
ex numeric::evalf(int /*level*/, PyObject* parent) const {
        PyObject *ans, *a = to_pyobject();
        if (parent == nullptr)
                parent = RR_get();
        if (not PyDict_CheckExact(parent)) {
                PyObject* dict = PyDict_New();
                if (dict == nullptr)
                        throw(std::runtime_error("PyDict_New returned NULL"));
                int r = PyDict_SetItemString(dict, "parent", parent);
                if (r<0)
                        throw(std::runtime_error("PyDict_SetItemString failed"));
                ans = py_funcs.py_float(a, dict);
                Py_DECREF(dict);
        }
        else
                ans = py_funcs.py_float(a, parent);
        Py_DECREF(a);
        if (ans == nullptr)
                throw (std::runtime_error("numeric::evalf(): error calling py_float()"));

        return ans;
}

const numeric numeric::try_py_method(const std::string& s) const
{
        PyObject *obj = to_pyobject();
        PyObject* ret = PyObject_CallMethod(obj,
                        const_cast<char*>(s.c_str()), NULL);
        Py_DECREF(obj);
        if (ret == nullptr) {
                PyErr_Clear();
                throw std::logic_error("");
        }
        
        return numeric(ret);
}

const numeric numeric::try_py_method(const std::string& s,
                const numeric& num2) const
{
        PyObject *obj1 = to_pyobject();
        PyObject *obj2 = num2.to_pyobject();
        PyObject *mstr = PyString_FromString(const_cast<char*>(s.c_str()));
        PyObject* ret = PyObject_CallMethodObjArgs(obj1, mstr, obj2, NULL);
        Py_DECREF(obj1);
        Py_DECREF(obj2);
        Py_DECREF(mstr);
        if (ret == nullptr) {
                PyErr_Clear();
                throw std::logic_error("");
        }
        return numeric(ret);
}

const numeric numeric::to_dict_parent(PyObject* obj) const
{
        PyObject *ret = nullptr;
        PyObject *the_parent = nullptr;
        PyObject *the_arg = to_pyobject();
        if (obj != nullptr and PyDict_Check(obj)) {
                PyObject* pkey = PyString_FromString(const_cast<char*>("parent"));
                the_parent = PyDict_GetItem(obj, pkey);
                Py_DECREF(pkey);
                if (the_parent != nullptr
                    and PyCallable_Check(the_parent)) {
                        ret = PyObject_CallFunctionObjArgs(the_parent, the_arg, nullptr);
                        Py_DECREF(the_arg);
                        if (ret == nullptr) {
                                PyErr_Clear();
                                throw std::logic_error("");
                        }
                        return numeric(ret);
                }
        }
        ret = PyObject_CallFunctionObjArgs(RR_get(), the_arg, NULL);
        if (ret == nullptr) {
                PyErr_Clear();
                ret = PyObject_CallFunctionObjArgs(CC_get(), the_arg, NULL);
                Py_DECREF(the_arg);
                if (ret == nullptr) {
                        PyErr_Clear();
                        throw std::logic_error("");
                }
        }
        else
                Py_DECREF(the_arg);
        return numeric(ret);
}

ex numeric::subs(const exmap & m, unsigned) const
{
        const numeric& im = imag();
        if (im.is_zero()) {
                if (is_zero() or is_one() or is_minus_one())
                        return *this;
                for (const auto & pair : m)
                        if (is_exactly_a<numeric>(pair.first)
                            and is_equal(ex_to<numeric>(pair.first)))
                                return pair.second;
                return *this;
        }

        const numeric& re = real();
        ex ret1 = re;
        ex ret2 = im;
        ex ret3 = I;
        bool changed1=false, changed2=false, changed3=false;
        for (const auto & pair : m) {
                if (not is_exactly_a<numeric>(pair.first))
                        continue;
                const numeric& p = ex_to<numeric>(pair.first);
                const numeric& pim = p.imag();
                const numeric& pre = p.real();
                if (pim.is_zero()) {
                        if (p.is_zero() or p.is_one() or p.is_minus_one())
                                continue;
                        if (re.is_equal(pre)) {
                                ret1 = pair.second;
                                changed1 = true;
                        }
                        if (im.is_equal(pre)) {
                                ret2 = pair.second;
                                changed2 = true;
                        }
                        continue;
                }
                if (pim.is_one() and pre.is_zero()) {
                        ret3 = pair.second;
                        changed3 = true;
                        continue;
                }
                if (re.is_equal(pre) and im.is_equal(pim))
                        return pair.second;
        }
        if (changed1 or changed2 or changed3)
                return ret1 + ret2*ret3;
        return *this;
}
        
///////////////////////////////////////////////////////////////////////

/** Real part of a number. */
const numeric numeric::real() const {
        switch (t) {
        case LONG:
        case MPZ:
        case MPQ:
                return *this;
        case PYOBJECT:
        {
                if (PyFloat_Check(v._pyobject))
                        return *this;
                if (PyComplex_Check(v._pyobject))
                        return PyComplex_RealAsDouble(v._pyobject);
                try {
                        return try_py_method("real");
                }
                catch (std::logic_error) {}
                try {
                        return try_py_method("real_part");
                }
                catch (std::logic_error) {}
                return *this;
        }
        default:
                stub("invalid type");
        }
}

/** Imaginary part of a number. */
const numeric numeric::imag() const {
        switch (t) {
        case LONG:
        case MPZ:
        case MPQ:
                return *_num0_p;
        case PYOBJECT:
        {
                if (PyFloat_Check(v._pyobject))
                        return *_num0_p;
                if (PyComplex_Check(v._pyobject))
                        return PyComplex_ImagAsDouble(v._pyobject);
                try {
                        return try_py_method("imag");
                }
                catch (std::logic_error) {}
                try {
                        return try_py_method("imag_part");
                }
                catch (std::logic_error) {}
                return *_num0_p;
        }
        default:
                stub("invalid type");
        }
}

/** Numerator.  Computes the numerator of rational numbers, rationalized
 *  numerator of complex if real and imaginary part are both rational numbers
 *  (i.e numer(4/3+5/6*I) == 8+5*I), the number carrying the sign in all other
 *  cases. */
const numeric numeric::numer() const {

        switch (t) {
        case LONG:
        case MPZ:
                return *this;
        case MPQ: {
                mpz_t bigint;
                mpz_init_set(bigint, mpq_numref(v._bigrat));
                return bigint;
        }
        case PYOBJECT: {
                PyObject *a;
                a = py_funcs.py_numer(v._pyobject);
                if (a == nullptr) py_error("numer");
                return a;
        }
        default:
                stub("invalid type -- numer() type not handled");
        }
}

/** Denominator.  Computes the denominator of rational numbers, common integer
 *  denominator of complex if real and imaginary part are both rational numbers
 *  (i.e denom(4/3+5/6*I) == 6), one in all other cases. */
const numeric numeric::denom() const {

        switch (t) {
        case LONG:
        case MPZ:
                return 1;
        case MPQ: {
                mpz_t bigint;
                mpz_init_set(bigint, mpq_denref(v._bigrat));
                return bigint;
        }
        case PYOBJECT: {
                PyObject *a;
                a = py_funcs.py_denom(v._pyobject);
                if (a == nullptr) py_error("denom");
                return a;
        }
        default:
                stub("invalid type -- denom() type not handled");
        }
}

const numeric numeric::floor() const {
        numeric d = denom();
        if (d.is_one())
                return *this;
        return numer().iquo(d);
}

const numeric numeric::frac() const {
        numeric d = denom();
        if (d.is_one())
                return 0;
        return *this - numer().iquo(d);
}

const numeric numeric::fibonacci() const {
    PY_RETURN(py_funcs.py_fibonacci);
}

const numeric numeric::exp(PyObject* parent) const {
        static numeric tentt20 = ex_to<numeric>(_num10_p->power(*_num20_p));
        // avoid Flint aborts
        if (real() > tentt20) {
                if (imag().is_zero())
                        return py_funcs.py_eval_infinity();
                else
                        return py_funcs.py_eval_unsigned_infinity();
        }
        if (real() < -tentt20)
                return ex_to<numeric>(_num0_p->evalf(0, parent));

        return arbfunc_0arg("exp", parent);
}

const numeric numeric::log(PyObject* parent) const {
        return arbfunc_0arg("log", parent);
}

// General log
const numeric numeric::log(const numeric &b, PyObject* parent) const {
        if (b.is_one()) {
                if (is_one())
		        throw (std::runtime_error("log(1,1) encountered"));
                else
                        return py_funcs.py_eval_unsigned_infinity();
        }
        if (b.is_zero())
                return *_num0_p;

        if ((t != LONG and t != MPZ and t != MPQ)
            or (b.t != LONG and b.t != MPZ and b.t != MPQ))
                return log(parent)/b.log(parent);
        
        bool israt;
        numeric ret = ratlog(b, israt);
        if (not israt)
                return log(parent)/b.log(parent);
        
        return ret;
}

// General log
// Handle special cases here that return MPZ/MPQ (or an infinity)
const numeric numeric::ratlog(const numeric &b, bool& israt) const {
        israt = true;
        if (b.is_one()) {
                if (is_one())
		        throw (std::runtime_error("log(1,1) encountered"));
                else
                        return py_funcs.py_eval_unsigned_infinity();
        }
        if (b.is_zero())
                return *_num0_p;

        if ((t != LONG and t != MPZ and t != MPQ)
            or (b.t != LONG and b.t != MPZ and b.t != MPQ)) {
                israt = false;
                return *_num0_p;
        }

        if (t == LONG and b.t == LONG) {
                if (b.v._long <= 0) {
                        israt = false;
                        return *_num0_p;
                }
                if (v._long == 0) {
                        return py_funcs.py_eval_neg_infinity();
                }
                int c = 0;
                std::ldiv_t ld;
                ld.quot = v._long;
                ld.rem = 0;
                do {
                        ld = std::div(ld.quot, b.v._long);
                        ++c;
                }
                while (ld.quot != 1 and ld.rem == 0);
                if (ld.quot != 1 or ld.rem != 0) {
                        israt = false;
                        return *_num0_p;
                }
                return c;
        }
        if (b.t == LONG)
                return ratlog(b.to_bigint(), israt);
        if (t == LONG)
                return to_bigint().ratlog(b, israt);
        if (t == MPZ and b.t == MPZ) {
                if (b > *_num0_p) {
                        mpz_t ret;
                        mpz_init(ret);
                        mp_bitcnt_t r = mpz_remove(ret, v._bigint, b.v._bigint);
                        if (mpz_cmp_ui(ret, 1) == 0) {
                                mpz_clear(ret);
                                return long(r);
                        }
                        mpz_clear(ret);
                }
                israt = false;
                return *_num0_p;
        }

        if (t == MPZ) {
                if (b.numer().is_one())
                        return -ratlog(b.denom(), israt);
                
                israt = false;
                return *_num0_p;
        }
        if (b.t == MPZ) {
                if (numer().is_one())
                        return -denom().ratlog(b, israt);
                
                israt = false;
                return *_num0_p;
        }

        // from here both MPQ
        numeric d = GiNaC::log(denom(), b.denom(), nullptr);
        if (numer().is_one() and b.numer().is_one())
                return d;
        numeric n = GiNaC::log(numer(), b.numer(), nullptr);
        if (n == d)
                return n;
        return log(nullptr)/b.log(nullptr);
}

const numeric numeric::sin() const {
        PY_RETURN(py_funcs.py_sin);
}

const numeric numeric::cos() const {
        PY_RETURN(py_funcs.py_cos);
}

const numeric numeric::tan() const {
        PY_RETURN(py_funcs.py_tan);
}

const numeric numeric::asin(PyObject* parent) const {
        return arbfunc_0arg("arcsin", parent);
}

const numeric numeric::acos(PyObject* parent) const {
        return arbfunc_0arg("arccos", parent);
}

const numeric numeric::atan(PyObject* parent) const {
        return arbfunc_0arg("arctan", parent);
}

const numeric numeric::atan(const numeric& y, PyObject* parent) const {
        PyObject *cparent = common_parent(*this, y);
        bool newdict = false;
        if (parent == nullptr) {
                parent = PyDict_New();
                PyDict_SetItemString(parent, "parent", cparent);
                newdict = true;
        }
        numeric rnum;
        if (not imag().is_zero()
            or not y.imag().is_zero()) {
                ex r = this->add(y.mul(I))
                        / (mul(*this) + y*y).power(*_num1_2_p);
                rnum = -ex_to<numeric>(r.evalf(0, parent)).log(parent) * I;
                Py_DECREF(cparent);
                return rnum;
        }
        if (y.is_zero()) {
                if (is_zero())
                        throw (std::runtime_error("atan2(): division by zero"));
                if (real().is_positive())
                        rnum = *_num0_p;
                else
                        rnum = ex_to<numeric>(Pi.evalf(0, parent));
        }
        else {
                if (is_zero())
                        rnum = ex_to<numeric>(Pi.evalf(0, parent)) / *_num2_p;
                else if (real().is_positive())
                        rnum = (y/(*this)).abs().atan(parent);
                else
                        rnum = (ex_to<numeric>(Pi.evalf(0, parent))
                                        - (y/(*this)).abs().atan(parent));
                if (not y.real().is_positive())
                        rnum = rnum.negative();
        }
        rnum = ex_to<numeric>(rnum.evalf(0, parent));
        Py_DECREF(cparent);
        if (newdict)
                Py_DECREF(parent);
        return rnum;
}

const numeric numeric::sinh(PyObject* parent) const {
        return (exp(parent) - negative().exp(parent)) / *_num2_p;
}

const numeric numeric::cosh(PyObject* parent) const {
        return (exp(parent) + negative().exp(parent)) / *_num2_p;
}

const numeric numeric::tanh(PyObject* parent) const {
        const numeric& e2x = exp(parent);
        const numeric& e2nx = negative().exp(parent);
        return (e2x - e2nx)/(e2x + e2nx);
}

const numeric numeric::asinh(PyObject* parent) const {
        return arbfunc_0arg("arcsinh", parent);
}

const numeric numeric::acosh(PyObject* parent) const {
        return arbfunc_0arg("arccosh", parent);
}

const numeric numeric::atanh(PyObject* parent) const {
        int prec = precision(*this, parent);
        PyObject* field = CBF(prec+15);
        PyObject* ret = CallBallMethod0Arg(field, const_cast<char*>("arctanh"), *this);
        Py_DECREF(field);

        numeric rnum(ret);
        if ((is_real() or imag().is_zero())
            and abs()<(*_num1_p))
                return ex_to<numeric>(rnum.real().evalf(0, parent));
        
        return ex_to<numeric>(rnum.evalf(0, parent));
}

const numeric numeric::Li2(const numeric &n, PyObject* parent) const {
        PyObject *cparent = common_parent(*this, n);
        if (parent == nullptr)
               parent = cparent;
        int prec = precision(*this, parent);
        PyObject* field = CBF(prec+15);
        PyObject* ret = CallBallMethod1Arg(field, const_cast<char*>("polylog"), *this, n);
        Py_DECREF(field);

        numeric rnum(ret), nret;
        if ((is_real() or imag().is_zero())
            and n.is_integer() and rnum.real()<(*_num1_p))
                nret = ex_to<numeric>(rnum.real().evalf(0, parent));
        else
                nret = ex_to<numeric>(rnum.evalf(0, parent));
        Py_DECREF(cparent);
        return nret;
}

const numeric numeric::lgamma(PyObject* parent) const {
        int prec = precision(*this, parent);
        PyObject* field = CBF(prec+15);
        PyObject* ret = CallBallMethod0Arg(field, const_cast<char*>("log_gamma"), *this);
        Py_DECREF(field);

        numeric rnum(ret);
        return ex_to<numeric>(rnum.evalf(0, parent));
}

const numeric numeric::gamma(PyObject* parent) const {
        return arbfunc_0arg("gamma", parent);
}

const numeric numeric::rgamma(PyObject* parent) const {
        return arbfunc_0arg("rgamma", parent);
}

const numeric numeric::psi(PyObject* parent) const {
        return arbfunc_0arg("psi", parent);
}

const numeric numeric::psi(const numeric& y) const {
        PY_RETURN2(py_funcs.py_psi2, y);
}

const numeric numeric::zeta() const {
        PY_RETURN(py_funcs.py_zeta);
}

const numeric numeric::stieltjes() const {
        PY_RETURN(py_funcs.py_stieltjes);
}

const numeric numeric::factorial() const {
        static long fac[] = {1, 1, 2, 6, 24, 120, 720,
                5040, 40320, 362880, 3628800, 39916800,
                479001600};
        if (is_integer()) {
                if (is_positive() and *this < 13)
                        return fac[to_long()];
                mpz_t bigint;
                mpz_init(bigint);
                mpz_fac_ui(bigint, to_long());
                return bigint;
        }
        PY_RETURN(py_funcs.py_factorial);
}

const numeric numeric::doublefactorial() const {
        PY_RETURN(py_funcs.py_doublefactorial);
}

const numeric numeric::binomial(const numeric &k) const {
        if ((t == LONG or t == MPZ)
            and k.is_integer()) {
                if (is_positive() and k.is_positive()
                    and *this < 13) {
                        static long fac[] = {1, 1, 2, 6, 24, 120, 720,
                                5040, 40320, 362880, 3628800, 39916800,
                                479001600};
                        long a = to_long();
                        long b = k.to_long();
                        if (b<=0 or b>12)
                                return *_num0_p;
                        return fac[a]/fac[b]/fac[a-b];
                }
                mpz_t bigint;
                mpz_init(bigint);
                if (t == MPZ) {
                        mpz_bin_ui(bigint, v._bigint, k.to_long());
                }
                else {
                        mpz_set_ui(bigint, v._long);
                        mpz_bin_ui(bigint, bigint, k.to_long());
                }
                return bigint;
        }
        numeric res;
        PyObject *nobj, *kobj;
        nobj = this->to_pyobject();
        kobj = k.to_pyobject();

        PyObject* m = PyImport_ImportModule("sage.arith.misc");
        if (m == nullptr)
                py_error("Error importing arith.misc");
        PyObject* binfunc = PyObject_GetAttrString(m, "binomial");
        if (binfunc == nullptr)
                py_error("Error getting binomial");

        PyObject* pyresult = PyObject_CallFunctionObjArgs(binfunc, nobj, kobj, NULL);
        Py_DECREF(kobj);
        Py_DECREF(nobj);
        Py_DECREF(m);
        Py_DECREF(binfunc);
        if (pyresult == nullptr) {
                throw(std::runtime_error("numeric::binomial(): python function binomial raised exception"));
        }
        if ( pyresult == Py_None ) {
                throw(std::runtime_error("numeric::binomial: python function binomial returned None"));
        }
        // convert output Expression to an ex
        ex eval_result = py_funcs.pyExpression_to_ex(pyresult);
        Py_DECREF(pyresult);
        if (PyErr_Occurred() != nullptr) {
                throw(std::runtime_error("numeric::binomial(): python function (Expression_to_ex) raised exception"));
        }
        return ex_to<numeric>(eval_result);
}

const numeric numeric::binomial(unsigned long n, unsigned long k)
{
        if (n < 13) {
                static long fac[] = {1, 1, 2, 6, 24, 120, 720,
                        5040, 40320, 362880, 3628800, 39916800,
                        479001600};
                if (k == 0)
                        return *_num1_p;
                if (k > n)
                        return *_num0_p;
                return fac[n]/fac[k]/fac[n-k];
        }
        mpz_t bigint;
        mpz_init(bigint);
        mpz_bin_uiui(bigint, n, k);
        return bigint;
}

const numeric numeric::bernoulli() const {
        PY_RETURN(py_funcs.py_bernoulli);
}

const numeric numeric::hypergeometric_pFq(const std::vector<numeric>& a, const std::vector<numeric>& b, PyObject *parent) const {
        PyObject *lista = py_tuple_from_numvector(a);
        PyObject *listb = py_tuple_from_numvector(b);
        PyObject *z = to_pyobject();

        // call opt.evalf_f with this tuple
        PyObject* m = PyImport_ImportModule("sage.functions.hypergeometric");
        if (m == nullptr)
                py_error("Error importing hypergeometric");
        PyObject* hypfunc = PyObject_GetAttrString(m, "hypergeometric");
        if (hypfunc == nullptr)
                py_error("Error getting hypergeometric attribute");

        if ((parent != nullptr) and PyDict_CheckExact(parent)) {
                Py_DECREF(z);
                z = ex_to<numeric>(this->evalf(0, parent)).to_pyobject();
        }
        PyObject* name = PyString_FromString(const_cast<char*>("_evalf_try_"));
        PyObject* pyresult = PyObject_CallMethodObjArgs(hypfunc, name, lista, listb, z, NULL);
        Py_DECREF(m);
        Py_DECREF(z);
        Py_DECREF(name);
        Py_DECREF(hypfunc);
        if (pyresult == nullptr) {
                throw(std::runtime_error("numeric::hypergeometric_pFq(): python function hypergeometric::_evalf_ raised exception"));
        }
        if ( pyresult == Py_None ) {
                throw(std::runtime_error("numeric::hypergeometric_pFq(): python function hypergeometric::_evalf_ returned None"));
        }
        // convert output Expression to an ex
        ex eval_result = py_funcs.pyExpression_to_ex(pyresult);
        Py_DECREF(pyresult);
        if (PyErr_Occurred() != nullptr) {
                throw(std::runtime_error("numeric::hypergeometric_pFq(): python function (Expression_to_ex) raised exception"));
        }
        return ex_to<numeric>(eval_result);
}

const numeric numeric::isqrt() const {
        PY_RETURN(py_funcs.py_isqrt);
}

const numeric numeric::sqrt() const {
        PY_RETURN(py_funcs.py_sqrt);
}

bool numeric::is_square() const
{
        if (is_negative())
                return false;
        if (is_zero() or is_one())
                return true;
        switch (t) {
        case LONG: {
                long rt = std::lround(std::sqrt(v._long));
                return rt * rt == v._long;
        }
        case MPZ:
                return mpz_perfect_square_p(v._bigint);
        case MPQ:
                return (mpz_perfect_square_p(mpq_numref(v._bigrat))
                and mpz_perfect_square_p(mpq_denref(v._bigrat)));
        case PYOBJECT:
        default:
                stub("invalid type: type not handled");
        }
        return false;
}

// Fast sqrt in case of square, symbolic sqrt else.
const ex numeric::sqrt_as_ex() const
{
        if (is_negative())
                return I * (-*this).sqrt_as_ex();
        if (is_zero())
                return _ex0;
        if (is_one())
                return _ex1;
        switch (t) {
        case LONG: {
                long rt = std::lround(std::sqrt(v._long));
                if (rt * rt == v._long)
                        return ex(numeric(rt));
                break;
        }
        case MPZ: {
                if (mpz_perfect_square_p(v._bigint)) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_sqrt(bigint, v._bigint);
                        return ex(numeric(bigint));
                }
                break;
        }
        case MPQ: {
                if (mpz_perfect_square_p(mpq_numref(v._bigrat))
                and mpz_perfect_square_p(mpq_denref(v._bigrat))) {
                        mpz_t bigint;
                        mpq_t bigrat, obigrat;
                        mpz_init(bigint);
                        mpq_init(bigrat);
                        mpq_init(obigrat);
                        mpz_sqrt(bigint, mpq_numref(v._bigrat));
                        mpq_set_z(bigrat, bigint);
                        mpz_sqrt(bigint, mpq_denref(v._bigrat));
                        mpq_set_z(obigrat, bigint);
                        mpq_div(bigrat, bigrat, obigrat);
                        mpz_clear(bigint);
                        mpq_clear(obigrat);
                        return ex(numeric(bigrat));
                }
                break;
        }
        case PYOBJECT:
                return sqrt();
        default:
                stub("invalid type: type not handled");
        }
        return (new GiNaC::power(ex(*this), _ex1_2))->setflag(status_flags::dynallocated | status_flags::evaluated);
}

const numeric numeric::abs() const {
        switch (t) {
        case LONG:
                if (*this >= 0)
                        return *this;
                return negative();
        case MPZ: {
                mpz_t bigint;
                mpz_init(bigint);
                mpz_abs(bigint, v._bigint);
                return bigint;
        }
        case MPQ: {
                mpq_t bigrat;
                mpq_init(bigrat);
                mpq_abs(bigrat, v._bigrat);
                return bigrat;
        }
        case PYOBJECT: {
                PyObject *ret = PyNumber_Absolute(v._pyobject);
                if (ret == NULL) {
                        PyErr_Clear();
                        return *this;
                }
                return ret;
        }
        default:
                stub("invalid type: type not handled");
        }
}

const numeric numeric::mod(const numeric &b) const {
        switch (t) {
        case LONG:
                if (b.t == LONG)
                        return v._long % b.v._long;
                if (b.t == MPZ)
                        return to_bigint().mod(b);
                throw std::runtime_error("unsupported type in numeric::md");
        case MPZ:
                if (b.t == LONG)
                        return mod(b.to_bigint());
                if (b.t == MPZ) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_mod(bigint, v._bigint, b.v._bigint);
                        return bigint;
                }
                throw std::runtime_error("unsupported type in numeric::md");
        case MPQ:
        case PYOBJECT: {
                PY_RETURN2(py_funcs.py_mod, b);
        }
        default:
                stub("invalid type: type not handled");
        }
}

const numeric numeric::_smod(const numeric &b) const {
        PY_RETURN2(py_funcs.py_smod, b);
}

const numeric numeric::irem(const numeric &b) const {
        switch (t) {
        case LONG:
                if (b.t == LONG)
                        return std::remainder(v._long, b.v._long);
                if (b.t == MPZ)
                        return to_bigint().irem(b);
                throw std::runtime_error("unsupported type in numeric::irem");
        case MPZ:
                if (b.t == LONG)
                        return irem(b.to_bigint());
                if (b.t == MPZ) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_fdiv_r(bigint, v._bigint, b.v._bigint);
                        return bigint;
                }
                throw std::runtime_error("unsupported type in numeric::irem");
        case MPQ:
        case PYOBJECT: {
                PY_RETURN2(py_funcs.py_irem, b);
        }
        default:
                stub("invalid type: type not handled");
        }
}

const numeric numeric::iquo(const numeric &b) const {
        switch (t) {
        case LONG:
                if (b.t == LONG)
                        return v._long / b.v._long;
                if (b.t == MPZ)
                        return to_bigint().iquo(b);
                throw std::runtime_error("unsupported type in numeric::iquo");
        case MPZ:
                if (b.t == LONG) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_fdiv_q_ui(bigint, v._bigint,
                        std::labs(b.v._long));
                        if (b.v._long < 0)
                                mpz_neg(bigint, bigint);
                        return bigint;
                }
                if (b.t == MPZ) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_fdiv_q(bigint, v._bigint, b.v._bigint);
                        return bigint;
                }
                throw std::runtime_error("unsupported type in numeric::iquo");
        case MPQ:
        case PYOBJECT:
        default:
                stub("invalid type: type not handled");
        }
        throw std::runtime_error("iquo: bad input");
}

const numeric numeric::iquo(const numeric &b, numeric& r) const {
        switch (t) {
        case LONG:
                if (b.t == LONG) {
                        auto ld = std::ldiv(v._long, b.v._long);
                        r = ld.rem;
                        return ld.quot;
                }
                if (b.t == MPZ)
                        return to_bigint().iquo(b, r);
                throw std::runtime_error("unsupported type in numeric::iquo");
        case MPZ:
                if (b.t == LONG) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        r = long(mpz_fdiv_q_ui(bigint, v._bigint,
                        std::labs(b.v._long)));
                        return bigint;
                }
                if (b.t == MPZ) {
                        mpz_t bigint, tmp;
                        mpz_init(bigint);
                        mpz_init(tmp);
                        mpz_fdiv_q(bigint, v._bigint, b.v._bigint);
                        mpz_mul(tmp, bigint, b.v._bigint);
                        mpz_sub(tmp, v._bigint, tmp);
                        r = numeric(tmp);
                        return bigint;
                }
                throw std::runtime_error("unsupported type in numeric::iquo");
        case MPQ:
        case PYOBJECT:
        default:
                stub("invalid type: type not handled");
        }
        throw std::runtime_error("iquo2: bad input");
}

const numeric numeric::gcd(const numeric &B) const
{
        if (is_zero())
                return B.abs();
        if (B.is_zero())
                return abs();
        numeric a,b;
        if (is_negative())
                a = negative();
        else
                a = *this;
        if (B.is_negative())
                b = B.negative();
        else
                b = B;
        if (a.is_one() or b.is_one())
                return *_num1_p;

        switch (a.t) {
        case LONG:
                switch (b.t) {
                case LONG: {
                        // C++17 will have a function for this
                        long al = a.v._long;
                        long bl = b.v._long;
                        long c;
                        while (al != 0) {
                                c = al;
                                al = bl % al;
                                bl = c;
                        }
                        return bl;
                }
                case MPZ:
                        return a.to_bigint().gcd(b);
                case MPQ:
                        return (a * b.denom()).gcd(b.numer()) / b.denom();
                case PYOBJECT:
                        return *_num1_p;
                default: stub("invalid type: type not handled");
                }
        case MPZ:
                switch (b.t) {
                case LONG: {
                        mpz_t bigint;
                        mpz_init(bigint);
                        long l = mpz_gcd_ui(bigint,
                        a.v._bigint,
                        std::labs(b.v._long));
                        mpz_clear(bigint);
                        return l;
                }
                case MPZ: {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_gcd(bigint, a.v._bigint, b.v._bigint);
                        return bigint;
                }
                case MPQ:
                        return (a * b.denom()).gcd(b.numer()) / b.denom();
                case PYOBJECT:
                        return *_num1_p;
                default:
                        stub("invalid type: type not handled");
                }
        case MPQ:
                if (b.t == PYOBJECT)
                        return *_num1_p;
                return (a.numer() * b.denom()).gcd(a.denom() * b.numer()) / (a.denom() * b.denom());
                stub("invalid type: type not handled");
        case PYOBJECT:
                        return *_num1_p;
        default:
                stub("invalid type: type not handled");
        }
}

const numeric numeric::lcm(const numeric &b) const
{
        if (is_zero() or b.is_zero())
                return *_num0_p;
        if (is_one())
               return b;
        if (b.is_one())
                return *this;
        switch (t) {
        case LONG:
                if (b.t == LONG) {
                        // C++17 will have a function for this
                        long a = v._long;
                        long bl = b.v._long;
                        long c;
                        while (a != 0) {
                                c = a;
                                a = bl % a;
                                bl = c;
                        }
                        return (v._long / bl) * b.v._long;
                }
                if (b.t == MPZ)
                        return to_bigint().lcm(b);
                throw std::runtime_error("unsupported type in numeric::lcm");
        case MPZ:
                if (b.t == LONG)
                        return lcm(b.to_bigint());
                if (b.t == MPZ) {
                        mpz_t bigint;
                        mpz_init(bigint);
                        mpz_lcm(bigint, v._bigint, b.v._bigint);
                        return bigint;
                }
                throw std::runtime_error("unsupported type in numeric::lcm");
        case MPQ:
        case PYOBJECT: {
                PY_RETURN2(py_funcs.py_lcm, b);
        }
        default:
                stub("invalid type: type not handled");
        }
}

// version of factor for factors fitting into long
void numeric::factorsmall(std::vector<std::pair<long, int>>& factors, long range) const
{
        if (is_one() or is_zero() or is_minus_one())
                return;
        switch (t) {
        case LONG:
                return to_bigint().factorsmall(factors, range);
        case MPZ: {
                fmpz_t f;
                fmpz_init(f);
                mpz_t bigint;
                mpz_init(bigint);
                mpz_abs(bigint, v._bigint);
                fmpz_set_mpz(f, bigint);
                fmpz_factor_t fs;
                fmpz_factor_init(fs);
                if (range == 0)
                        fmpz_factor(fs, f);
                else
                        fmpz_factor_trial_range(fs, f, 0, range);
                for (slong i = 0; i < fs->num; ++i) {
                        fmpz_get_mpz(bigint, fs->p + i);
                        factors.push_back(std::make_pair(mpz_get_si(bigint),
                        int(*(fs->exp + i))));
                }
                mpz_clear(bigint);
                fmpz_factor_clear(fs);
                fmpz_clear(f);
                return;
        }
        case MPQ:
                return to_bigint().factorsmall(factors);
        case PYOBJECT:
        default:
                stub("invalid type: type not handled");
        }
}

// may take long time if no range given
void numeric::factor(std::vector<std::pair<numeric, int>>& factors, long range) const
{
        if (is_one() or is_minus_one())
                return;
        switch (t) {
        case LONG: {
                std::vector<std::pair<long, int>> f;
                factorsmall(f, range);
                for (const auto& p : f)
                        factors.emplace_back(numeric(p.first), p.second);
                return;
        }
        case MPZ: {
                fmpz_t f;
                fmpz_init(f);
                mpz_t bigint, tmp;
                mpz_init(bigint);
                mpz_abs(bigint, v._bigint);
                fmpz_set_mpz(f, bigint);
                fmpz_factor_t fs;
                fmpz_factor_init(fs);
                if (range == 0)
                        fmpz_factor(fs, f);
                else
                        fmpz_factor_trial_range(fs, f, 0, range);
                for (slong i = 0; i < fs->num; ++i) {
                        mpz_init(tmp);
                        fmpz_get_mpz(tmp, fs->p + i);
                        if (range != 0)
                                for (int j = 0; j < int(*(fs->exp + i)); ++j) {
                                        mpz_divexact(bigint,
                                        bigint,
                                        tmp);
                                }
                        factors.emplace_back(numeric(tmp),
                                        int(*(fs->exp + i)));
                }
                fmpz_clear(f);
                fmpz_factor_clear(fs);
                if (range != 0 and mpz_cmp_ui(bigint, 1) != 0)
                        factors.push_back(std::make_pair(numeric(bigint), 1));
                else
                        mpz_clear(bigint);
                return;
        }
        case MPQ:
                to_bigint().factor(factors);
                return;
        case PYOBJECT:
        default:
                stub("invalid type: type not handled");
        }
}

static void setDivisors(long n, int i, std::set<int>& divisors,
        const std::vector<std::pair<long, int>>& factors)
// author of this little gem:
// http://stackoverflow.com/users/3864373/rocky-johnson
{
    for (size_t j = i; j<factors.size(); j++) {
        long x = factors[j].first * n;
        for (int k = 0; k<factors[j].second; k++) {
            divisors.insert(x);
            setDivisors(x, j + 1, divisors, factors);
            x *= factors[j].first;
        }
    }
}

void numeric::divisors(std::set<int>& divs) const
{
        divs.insert(1);
        if (is_one() or is_zero() or is_minus_one())
                return;
        switch (t) {
        case LONG:
        case MPZ: {
                std::vector<std::pair<long, int>> factors;
                factorsmall(factors);
                setDivisors(1, 0, divs, factors);
                return;
        }
        case MPQ:
                return to_bigint().divisors(divs);
        case PYOBJECT:
        default:
                stub("invalid type: type not handled");
        }
}

//////////
// global 
//////////

void coerce(numeric& new_left, numeric& new_right, const numeric& left, const numeric& right) {
        verbose("coerce");
        // Return a numeric 
        if (left.t == right.t) {
                new_left = left;
                new_right = right;
                return;
        }
        PyObject *o;
        switch (left.t) {
        case LONG:
                switch (right.t) {
                case MPZ:
                        if (mpz_fits_sint_p(right.v._bigint)) {
                                new_right = numeric(mpz_get_si(right.v._bigint));
                                new_left = left;
                        } else {
                                numeric n;
                                mpz_init(n.v._bigint);
                                n.t = MPZ;
                                mpz_set_si(n.v._bigint, left.v._long);
                                n.hash = (left.v._long == -1) ? -2 : left.v._long;
                                new_left = n;
                                new_right = right;
                        }
                        return;
                case MPQ: {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_si(bigrat, left.v._long, 1);
                        new_left = numeric(bigrat);
                        new_right = right;
                        return;
                }
                case PYOBJECT: {
                        mpz_t bigint;
                        mpz_init_set_si(bigint, left.v._long);
                        o = py_funcs.py_integer_from_mpz(bigint);
                        new_left = numeric(o, true);
                        new_right = right;
                        mpz_clear(bigint);
                        return;
                }
                default:
                        std::cerr << "type = " << right.t << "\n";
                        stub("** invalid coercion -- left MPZ**");
                }
        case MPZ:
                switch (right.t) {
                case LONG:
                        if (mpz_fits_sint_p(left.v._bigint)) {
                                new_left = numeric(mpz_get_si(left.v._bigint));
                                new_right = right;
                        } else {
                                numeric n;
                                mpz_init(n.v._bigint);
                                n.t = MPZ;
                                mpz_set_si(n.v._bigint, right.v._long);
                                n.hash = (right.v._long == -1) ? -2 : right.v._long;
                                new_right = n;
                                new_left = left;
                        }
                        return;
                case MPQ: {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_z(bigrat, left.v._bigint);
                        new_left = numeric(bigrat);
                        new_right = right;
                        return;
                }
                case PYOBJECT: {
                        mpz_t bigint;
                        mpz_init_set(bigint, left.v._bigint);
                        o = py_funcs.py_integer_from_mpz(bigint);
                        new_left = numeric(o, true);
                        new_right = right;
                        mpz_clear(bigint);
                        return;
                }
                default:
                        std::cerr << "type = " << right.t << "\n";
                        stub("** invalid coercion -- left MPZ**");
                }
        case MPQ:
                switch (right.t) {
                case LONG: {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_si(bigrat, right.v._long, 1);
                        new_right = numeric(bigrat);
                        new_left = left;
                        return;
                }
                case MPZ: {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set_z(bigrat, right.v._bigint);
                        new_left = left;
                        new_right = numeric(bigrat);
                        return;
                }
                case PYOBJECT: {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set(bigrat, left.v._bigrat);
                        o = py_funcs.py_rational_from_mpq(bigrat);
                        mpq_clear(bigrat);
                        new_left = numeric(o, true);
                        new_right = right;
                        return;
                }
                default:
                        std::cerr << "type = " << right.t << "\n";
                        stub("** invalid coercion -- left MPQ**");
                }
        case PYOBJECT:
                new_left = left;
                switch (right.t) {
                case LONG: {
                        mpz_t bigint;
                        mpz_init_set_si(bigint, right.v._long);
                        o = py_funcs.py_integer_from_mpz(bigint);
                        mpz_clear(bigint);
                        new_right = numeric(o, true);
                        return;
                }
                case MPZ: {
                        mpz_t bigint;
                        mpz_init_set(bigint, right.v._bigint);
                        o = py_funcs.py_integer_from_mpz(bigint);
                        mpz_clear(bigint);
                        new_right = numeric(o, true);
                        return;
                }
                case MPQ: {
                        mpq_t bigrat;
                        mpq_init(bigrat);
                        mpq_set(bigrat, right.v._bigrat);
                        o = py_funcs.py_rational_from_mpq(bigrat);
                        mpq_clear(bigrat);
                        new_right = numeric(o, true);
                        return;
                }
                default:
                        std::cerr << "type = " << right.t << "\n";
                        stub("** invalid coercion -- left PYOBJECT**");
                }
        default:
                std::cerr << "type = " << left.t << "\n";
                stub("** invalid coercion **");
        }
}

/** Imaginary unit.  This is not a constant but a numeric since we are
 *  natively handing complex numbers anyways, so in each expression containing
 *  an I it is automatically eval'ed away anyhow. */

/** Exponential function.
 *
 *  @return  arbitrary precision numerical exp(x). */
const numeric exp(const numeric &x, PyObject* parent) {
        return x.exp(parent);
}

/** Natural logarithm.
 *
 *  @param x complex number
 *  @return  arbitrary precision numerical log(x).
 *  @exception pole_error("log(): logarithmic pole",0) */
const numeric log(const numeric &x, PyObject* parent) {
        return x.log(parent);
}

// General log
const numeric log(const numeric &x, const numeric &b, PyObject* parent) {
        return x.log(b, parent);
}

/** Numeric sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical sin(x). */
const numeric sin(const numeric &x) {
        return x.sin();
}

/** Numeric cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical cos(x). */
const numeric cos(const numeric &x) {
        return x.cos();
}

/** Numeric tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical tan(x). */
const numeric tan(const numeric &x) {
        return x.tan();
}

/** Numeric inverse sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical asin(x). */
const numeric asin(const numeric &x, PyObject* parent) {
        return x.asin(parent);
}

/** Numeric inverse cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical acos(x). */
const numeric acos(const numeric &x, PyObject* parent) {
        return x.acos(parent);
}

/** Numeric arcustangent.
 *
 *  @param x complex number
 *  @return atan(x)
 *  @exception pole_error("atan(): logarithmic pole",0) if x==I or x==-I. */
const numeric atan(const numeric &x, PyObject* parent) {
        if (!x.is_real() &&
                x.real().is_zero() &&
                abs(x.imag()).is_one())
                throw pole_error("atan(): logarithmic pole", 0);
        return x.atan(parent);
}

/** Numeric arcustangent of two arguments, analytically continued in a suitable way.
 *
 *  @param y complex number
 *  @param x complex number
 *  @return -I*log((x+I*y)/sqrt(x^2+y^2)), which is equal to atan(y/x) if y and
 *    x are both real.
 *  @exception pole_error("atan(): logarithmic pole",0) if y/x==+I or y/x==-I. */
const numeric atan(const numeric &y, const numeric &x, PyObject* parent) {
        return x.atan(y, parent);
}

/** Numeric hyperbolic sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical sinh(x). */
const numeric sinh(const numeric &x, PyObject* parent) {
        return x.sinh(parent);
}

/** Numeric hyperbolic cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical cosh(x). */
const numeric cosh(const numeric &x, PyObject* parent) {
        return x.cosh(parent);
}

/** Numeric hyperbolic tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical tanh(x). */
const numeric tanh(const numeric &x, PyObject* parent) {
        return x.tanh(parent);
}

/** Numeric inverse hyperbolic sine (trigonometric function).
 *
 *  @return  arbitrary precision numerical asinh(x). */
const numeric asinh(const numeric &x, PyObject* parent) {
        return x.asinh(parent);
}

/** Numeric inverse hyperbolic cosine (trigonometric function).
 *
 *  @return  arbitrary precision numerical acosh(x). */
const numeric acosh(const numeric &x, PyObject* parent) {
        return x.acosh(parent);
}

/** Numeric inverse hyperbolic tangent (trigonometric function).
 *
 *  @return  arbitrary precision numerical atanh(x). */
const numeric atanh(const numeric &x, PyObject* parent) {
        return x.atanh(parent);
}


const numeric Li2(const numeric &x, PyObject* parent) {
        try {
                return x.try_py_method("dilog");
        }
        catch (std::logic_error) {}
        try {
                return x.try_py_method("polylog", *_num2_p);
        }
        catch (std::logic_error) {}

        return x.Li2(*_num2_p, parent);
}

const numeric Li2(const numeric &n, const numeric &x, PyObject* parent) {
        return x.Li2(n, parent);
}

/** Evaluation of Riemann's Zeta function.  */
const numeric zeta(const numeric &x) {
        return x.zeta();
}

const numeric stieltjes(const numeric &x) {
        return x.stieltjes();
}

const numeric lgamma(const numeric &x, PyObject* parent) {
        return x.lgamma(parent);
}

/** The Gamma function. */
const numeric gamma(const numeric &x, PyObject* parent) {
        return x.gamma(parent);
}

/** The psi function (aka polygamma function). */
const numeric psi(const numeric &x, PyObject* parent) {
        return x.psi(parent);
}

/** The psi functions (aka polygamma functions). */
const numeric psi(const numeric &n, const numeric &x) {
        return n.psi(x);
}
const numeric beta(const numeric &x, const numeric &y, PyObject* parent)
{
        PyObject *cparent = common_parent(x, y);
        if (parent == nullptr)
               parent = cparent;
        numeric ret = (x+y).rgamma(parent) * x.gamma(parent) * y.gamma(parent);
        Py_DECREF(cparent);
        return ret;
}

/** Factorial combinatorial function.
 *
 *  @param n  integer argument >= 0
 *  @exception range_error (argument must be integer >= 0) */
const numeric factorial(const numeric &n) {
        return n.factorial();
}

/** The double factorial combinatorial function.  (Scarcely used, but still
 *  useful in cases, like for exact results of gamma(n+1/2) for instance.)
 *
 *  @param n  integer argument >= -1
 *  @return n!! == n * (n-2) * (n-4) * ... * ({1|2}) with 0!! == (-1)!! == 1
 *  @exception range_error (argument must be integer >= -1) */
const numeric doublefactorial(const numeric &n) {
        return n.doublefactorial();
}

/** The Binomial coefficients.  It computes the binomial coefficients.  For
 *  integer n and k and positive n this is the number of ways of choosing k
 *  objects from n distinct objects.  If n is negative, the formula
 *  binomial(n,k) == (-1)^k*binomial(k-n-1,k) is used to compute the result. */
const numeric binomial(const numeric &n, const numeric &k) {
        return n.binomial(k);
}

/** Bernoulli number.  The nth Bernoulli number is the coefficient of x^n/n!
 *  in the expansion of the function x/(e^x-1).
 *
 *  @return the nth Bernoulli number (a rational number).
 *  @exception range_error (argument must be integer >= 0) */
const numeric bernoulli(const numeric &n) {
        return n.bernoulli();
}

const numeric hypergeometric_pFq(const std::vector<numeric>& a, const std::vector<numeric>& b, const numeric &z, PyObject* parent) {
        return z.hypergeometric_pFq(a, b, parent);
}

/** Fibonacci number.  The nth Fibonacci number F(n) is defined by the
 *  recurrence formula F(n)==F(n-1)+F(n-2) with F(0)==0 and F(1)==1.
 *
 *  @param n an integer
 *  @return the nth Fibonacci number F(n) (an integer number)
 *  @exception range_error (argument must be an integer) */
const numeric fibonacci(const numeric &n) {
        return n.fibonacci();
}

/** Absolute value. */
const numeric abs(const numeric& x) {
        return x.abs();
}

/** Modulus (in positive representation).
 *  In general, mod(a,b) has the sign of b or is zero, and rem(a,b) has the
 *  sign of a or is zero. This is different from Maple's modp, where the sign
 *  of b is ignored. It is in agreement with Mathematica's Mod.
 *
 *  @return a mod b in the range [0,abs(b)-1] with sign of b if both are
 *  integer, 0 otherwise. */
const numeric mod(const numeric &a, const numeric &b) {
        return a.mod(b);
}

///** Modulus (in symmetric representation).
// *  Equivalent to Maple's mods.
// *
// *  @return a mod b in the range [-iquo(abs(b)-1,2), iquo(abs(b),2)]. */
const numeric smod(const numeric &a, const numeric &b) {
        return a._smod(b);
}

/** Numeric integer remainder.
 *  Equivalent to Maple's irem(a,b) as far as sign conventions are concerned.
 *  In general, mod(a,b) has the sign of b or is zero, and irem(a,b) has the
 *  sign of a or is zero.
 *
 *  @return remainder of a/b if both are integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric irem(const numeric &a, const numeric &b) {
        return a.irem(b);
}

/** Numeric integer quotient.
 *  Equivalent to Maple's iquo as far as sign conventions are concerned.
 *  
 *  @return truncated quotient of a/b if both are integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric iquo(const numeric &a, const numeric &b) {
        return a.iquo(b);
}

/** Numeric integer quotient.
 *  Equivalent to Maple's iquo(a,b,'r') it obeyes the relation
 *  r == a - iquo(a,b,r)*b.
 *
 *  @return truncated quotient of a/b and remainder stored in r if both are
 *  integer, 0 otherwise.
 *  @exception overflow_error (division by zero) if b is zero. */
const numeric iquo(const numeric &a, const numeric &b, numeric &r) {
        return a.iquo(b, r);
}

/** Greatest Common Divisor.
 *   
 *  @return  The GCD of two numbers if both are integer, a numerical 1
 *  if they are not. */
const numeric gcd(const numeric &a, const numeric &b) {
        return a.gcd(b);
}

/** Least Common Multiple.
 *   
 *  @return  The LCM of two numbers if both are integer, the product of those
 *  two numbers if they are not. */
const numeric lcm(const numeric &a, const numeric &b) {
        return a.lcm(b);
}

/** Numeric square root.
 *  If possible, sqrt(x) should respect squares of exact numbers, i.e. sqrt(4)
 *  should return integer 2.
 *
 *  @param x numeric argument
 *  @return square root of x. Branch cut along negative real axis, the negative
 *  real axis itself where imag(x)==0 and real(x)<0 belongs to the upper part
 *  where imag(x)>0. */
const numeric sqrt(const numeric &x) {
        return x.sqrt();
}

/** Integer numeric square root. */
const numeric isqrt(const numeric &x) {
        return x.isqrt();
}

/** Floating point evaluation of Sage's constants. */
ex ConstantEvalf(unsigned serial, PyObject* dict) {
        if (dict == nullptr) {
                dict = PyDict_New();
                PyDict_SetItemString(dict, "parent", CC_get());
        }
        PyObject* x = py_funcs.py_eval_constant(serial, dict);
        if (x == nullptr) py_error("error getting digits of constant");
        return x;
}

ex UnsignedInfinityEvalf(unsigned serial, PyObject* parent) {
        PyObject* x = py_funcs.py_eval_unsigned_infinity();
        return x;
}

ex InfinityEvalf(unsigned serial, PyObject* parent) {
        PyObject* x = py_funcs.py_eval_infinity();
        return x;
}

ex NegInfinityEvalf(unsigned serial, PyObject* parent) {
        PyObject* x = py_funcs.py_eval_neg_infinity();
        return x;
}


} // namespace GiNaC
