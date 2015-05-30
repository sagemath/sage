/*******************************************************************************
    SAGE: Open Source Mathematical Software
        Copyright (C) 2008 William Stein <wstein@gmail.com>
        Copyright (C) 2008-2009 Burcin Erocal <burcin@erocal.org>
   Distributed under the terms of the GNU General Public License (GPL),
   version 2 or any later version.  The full text of the GPL is available at:
                   http://www.gnu.org/licenses/
*******************************************************************************/

#include "ccobject.h"
#include <pynac/ginac.h>
#include <pynac/extern_templates.h>
#include <string>

using namespace GiNaC;

void list_symbols(const ex& e, std::set<ex, ex_is_less> &s)
{
    if (is_a<symbol>(e)) {
        s.insert(e);
    } else {
        for (size_t i=0; i<e.nops(); i++)
            list_symbols(e.op(i), s);
    }
}


ex g_function_evalv(unsigned serial, exvector& vec, bool hold)
{
    if (hold)
        return function(serial, vec).hold();
    return function(serial, vec);
}

ex g_function_eval0(unsigned serial, bool hold)
{
    if (hold)
        return function(serial).hold();
    return function(serial);
}

ex g_function_eval1(unsigned serial, const ex& arg1, bool hold)
{
    if (hold)
        return function(serial, arg1).hold();
    return function(serial, arg1);
}

ex g_function_eval2(unsigned serial, const ex& arg1, const ex& arg2, bool hold)
{
    if (hold)
        return function(serial, arg1, arg2).hold();
    return function(serial, arg1, arg2);
}

ex g_function_eval3(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, bool hold)
{
    if (hold)
        return function(serial, arg1, arg2, arg3).hold();
    return function(serial, arg1, arg2, arg3);
}


bool relational_to_bool(const ex& e) {
    if (ex_to<relational>(e))
        return 1;
    else
        return 0;
}

bool g_is_a_terminating_series(const ex& e) {
    if (is_a<pseries>(e)) {
        return (ex_to<pseries>(e)).is_terminating();
    }
    return false;
}

relational::operators relational_operator(const ex& e) {
    // unsafe cast -- be damn sure the input is a relational.
    return (ex_to<relational>(e)).the_operator();
}

relational::operators switch_operator(relational::operators o) {
    switch(o) {
    case relational::less:
        return relational::greater;
    case relational::less_or_equal:
        return relational::greater_or_equal;
    case relational::greater:
        return relational::less;
    case relational::greater_or_equal:
        return relational::less_or_equal;
    default:
        return o;
    }
}

bool is_negative(ex x) {
    if (is_a<numeric>(x)) {
        return (ex_to<numeric>(x)).is_negative();
    }
    return false;
}

PyObject* py_object_from_numeric(ex x) {
    return (ex_to<numeric>(x)).to_pyobject();
}

template <class T>
PyObject* _to_PyString_latex(const T *x)
{
  std::ostringstream instore;
  instore << latex << (*x);
  return PyString_FromString(instore.str().data());
}

constant* GConstant_construct(void* mem, char* name, char* texname, unsigned domain)
{
  return new(mem) constant(name, ConstantEvalf, texname, domain);
}

#define ASSIGN_WRAP(x,y) x = y

#define ADD_WRAP(x,y) (x)+(y)
#define SUB_WRAP(x,y) (x)-(y)
#define MUL_WRAP(x,y) (x)*(y)
#define DIV_WRAP(x,y) (x)/(y)

#define LT_WRAP(x,y)  (x)<(y)
#define EQ_WRAP(x,y)  (x)==(y)
#define GT_WRAP(x,y)  (x)>(y)
#define LE_WRAP(x,y)  (x)<=(y)
#define NE_WRAP(x,y)  (x)!=(y)
#define GE_WRAP(x,y)  (x)>=(y)

#define HOLD(fn, x, hold_flag) hold_flag ? ex(fn(x).hold()) : fn(x)

#define HOLD2(fn, x, y, hold_flag) hold_flag ? ex(fn(x,y).hold()) : fn(x,y)

// declare the constant 1/2 so we can use it for a custom square root function
namespace GiNaC {

extern const ex _ex1_2;

}
