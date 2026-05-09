/** @file inifcns_hyperg.cpp
 *
 *  Implementation of hypergeometric functions.
 *
 *  (C) 2016 Ralf Stephan <ralf@ark.in-berlin.de>
 */

/*
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

#include <Python.h>
#include "py_funcs.h"
#include "inifcns.h"
#include "ex.h"
#include "constant.h"
#include "infinity.h"
#include "numeric.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "symbol.h"
#include "pseries.h"
#include "utils.h"

#include <utility>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <memory>

#if PY_MAJOR_VERSION > 2
#define PyString_FromString PyUnicode_FromString
#endif

namespace GiNaC {

inline void py_error(const char* errmsg) {
        throw std::runtime_error((PyErr_Occurred() != nullptr) ? errmsg:
                        "pyerror() called but no error occurred!");
}

// Creates the hypergeometric Python BuiltinFunction object
ex _2F1(const ex& a, const ex& b, const ex& c, ex x)
{
        exvector avec, bvec;
        avec.push_back(a);
        avec.push_back(b);
        bvec.push_back(c);
        PyObject *lista = py_funcs.exvector_to_PyTuple(avec);
        PyObject *listb = py_funcs.exvector_to_PyTuple(bvec);
        PyObject *z = py_funcs.ex_to_pyExpression(std::move(x));

        PyObject* m = PyImport_ImportModule("sage.functions.hypergeometric");
        if (m == nullptr)
                py_error("Error importing hypergeometric");
        PyObject* hypfunc = PyObject_GetAttrString(m, "hypergeometric");
        if (hypfunc == nullptr)
                py_error("Error getting hypergeometric attribute");

        PyObject* name = PyString_FromString(const_cast<char*>("__call__"));
        PyObject* pyresult = PyObject_CallMethodObjArgs(hypfunc, name, lista, listb, z, NULL);
        Py_DECREF(m);
        Py_DECREF(name);
        Py_DECREF(hypfunc);
        if (pyresult == nullptr) {
                throw(std::runtime_error("numeric::hypergeometric_pFq(): python function hypergeometric::__call__ raised exception"));
        }
        if ( pyresult == Py_None ) {
                throw(std::runtime_error("numeric::hypergeometric_pFq(): python function hypergeometric::__call__ returned None"));
        }
        // convert output Expression to an ex
        ex eval_result = py_funcs.pyExpression_to_ex(pyresult);
        Py_DECREF(pyresult);
        if (PyErr_Occurred() != nullptr) {
                throw(std::runtime_error("numeric::hypergeometric_pFq(): python function (Expression_to_ex) raised exception"));
        }
        return eval_result;
}


//////////
// Appell F1 function
//////////

static ex appellf1_evalf(const ex& a, const ex& b1, const ex& b2, const ex& c, const ex& x, const ex& y, PyObject* parent)
{
	/*if (is_exactly_a<numeric>(a) and is_exactly_a<numeric>(b1)
            and is_exactly_a<numeric>(b2) and is_exactly_a<numeric>(c)
            and is_exactly_a<numeric>(x) and is_exactly_a<numeric>(y)) {
                return appellf1(ex_to<numeric>(a), ex_to<numeric>(b1),
                                ex_to<numeric>(b2), ex_to<numeric>(c),
                                ex_to<numeric>(x), ex_to<numeric>(y));
	}*/
	return appell_F1(a, b1, b2, c, x, y).hold();
}

static ex appellf1_eval(const ex& a, const ex& b1, const ex& b2, const ex& c, const ex& x, const ex& y)
{
	if (is_exactly_a<numeric>(a) and is_exactly_a<numeric>(b1)
            and is_exactly_a<numeric>(b2) and is_exactly_a<numeric>(c)
            and is_exactly_a<numeric>(x) and is_exactly_a<numeric>(y)) {
                return appellf1_evalf(a, b1, b2, c, x, y, nullptr);
	}

        if (x.is_zero())
                return _2F1(a, b2, c, y);
        if (y.is_zero())
                return _2F1(a, b1, c, x);
        if (x.is_equal(y))
                return _2F1(a, b1+b2, c, x);
        if (c.is_equal(b1+b2))
                return power(1-y, -a) * _2F1(a, b1, b1+b2, (x-y)/(_ex1-y));

	return appell_F1(a, b1, b2, c, x, y).hold();
}

static ex appellf1_deriv(const ex& a, const ex& b1, const ex& b2, const ex& c, const ex& x, const ex& y, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==4 || deriv_param==5);
	
	// d/dx F1(a,b1,b2,c,x,y) --> a*b1/c*F1(a+1,b1+1,b2,c+1,x,y)
	// d/dy F1(a,b1,b2,c,x,y) --> a*b2/c*F1(a+1,b1,b2+1,c+1,x,y)
        if (deriv_param == 4)
        	return mul(mul(mul(a, b1), pow(c, -1)), appell_F1(a+1, b1+1, b2, c+1, x, y));
        return mul(mul(mul(a, b2), pow(c, -1)), appell_F1(a+1, b1, b2+1, c+1, x, y));
}

REGISTER_FUNCTION(appell_F1, eval_func(appellf1_eval).
                        evalf_func(appellf1_evalf).
                        derivative_func(appellf1_deriv).
		        latex_name("\\operatorname{F_1}"));

} // namespace GiNaC
