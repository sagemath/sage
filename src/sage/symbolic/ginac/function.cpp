/** @file function.cpp
 *
 *  Implementation of class of symbolic functions. */

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

#define register
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "py_funcs.h"
#include "function.h"
#include "operators.h"
#include "fderivative.h"
#include "ex.h"
#include "lst.h"
#include "print.h"
#include "power.h"
#include "relational.h"
#include "archive.h"
#include "inifcns.h"
#include "tostring.h"
#include "utils.h"
#include "remember.h"
#include "symbol.h"
#include "cmatcher.h"
#include "wildcard.h"
#include "expairseq.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <list>
#include <limits>
#ifdef DO_GINAC_ASSERT
#  include <typeinfo>
#endif

namespace GiNaC {

//////////
// helper class function_options
//////////

function_options::function_options()
{
	initialize();
}

function_options::function_options(std::string const & n, std::string const & tn)
{
	initialize();
	set_name(n, tn);
}

function_options::function_options(std::string const & n, unsigned np)
{
        static std::string empty;
	initialize();
        set_name(n, empty);
	nparams = np;
}

function_options::~function_options()
{
	// nothing to clean up at the moment
}

void function_options::initialize()
{
        static std::string s1("unnamed_function"), s2("\\mbox{unnamed}");
	set_name(s1, s2);
	nparams = 0;
	eval_f = real_part_f = imag_part_f = conjugate_f = derivative_f
            = pynac_eval_f = expl_derivative_f = power_f = series_f
            = subs_f = nullptr;
	evalf_f = nullptr;
	evalf_params_first = true;
	apply_chain_rule = true;
	use_return_type = false;
	eval_use_exvector_args = false;
	evalf_use_exvector_args = false;
	conjugate_use_exvector_args = false;
	real_part_use_exvector_args = false;
	imag_part_use_exvector_args = false;
	derivative_use_exvector_args = false;
        expl_derivative_use_exvector_args = false;
	power_use_exvector_args = false;
	series_use_exvector_args = false;
	print_use_exvector_args = false;
	use_remember = false;
	python_func = 0;
	functions_with_same_name = 1;
	symtree = 0;
}

function_options & function_options::set_name(std::string const & n,
                                              std::string const & tn)
{
        name.assign(n);
	if (tn.empty()) {
		TeX_name.assign("{\\rm ");
                TeX_name += n;
                TeX_name.append("}");
        }
	else
	        TeX_name.assign(tn);
	return *this;
}

function_options & function_options::latex_name(std::string const & tn)
{
	TeX_name = tn;
	return *this;
}

// the following lines have been generated for max. 14 parameters
function_options & function_options::eval_func(eval_funcp_1 e)
{
	test_and_set_nparams(1);
	eval_f = eval_funcp(e);
        pynac_eval_f = eval_f;
	return *this;
}
function_options & function_options::eval_func(eval_funcp_2 e)
{
	test_and_set_nparams(2);
	eval_f = eval_funcp(e);
        pynac_eval_f = eval_f;
	return *this;
}
function_options & function_options::eval_func(eval_funcp_3 e)
{
	test_and_set_nparams(3);
	eval_f = eval_funcp(e);
        pynac_eval_f = eval_f;
	return *this;
}
function_options & function_options::eval_func(eval_funcp_6 e)
{
	test_and_set_nparams(6);
	eval_f = eval_funcp(e);
        pynac_eval_f = eval_f;
	return *this;
}

function_options & function_options::evalf_func(evalf_funcp_1 ef)
{
	test_and_set_nparams(1);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_2 ef)
{
	test_and_set_nparams(2);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_3 ef)
{
	test_and_set_nparams(3);
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_6 ef)
{
	test_and_set_nparams(6);
	evalf_f = evalf_funcp(ef);
	return *this;
}

function_options & function_options::conjugate_func(conjugate_funcp_1 c)
{
	test_and_set_nparams(1);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_2 c)
{
	test_and_set_nparams(2);
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_3 c)
{
	test_and_set_nparams(3);
	conjugate_f = conjugate_funcp(c);
	return *this;
}

function_options & function_options::real_part_func(real_part_funcp_1 c)
{
	test_and_set_nparams(1);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_2 c)
{
	test_and_set_nparams(2);
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_3 c)
{
	test_and_set_nparams(3);
	real_part_f = real_part_funcp(c);
	return *this;
}

function_options & function_options::imag_part_func(imag_part_funcp_1 c)
{
	test_and_set_nparams(1);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_2 c)
{
	test_and_set_nparams(2);
	imag_part_f = imag_part_funcp(c);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_3 c)
{
	test_and_set_nparams(3);
	imag_part_f = imag_part_funcp(c);
	return *this;
}

function_options & function_options::derivative_func(derivative_funcp_1 d)
{
	test_and_set_nparams(1);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_2 d)
{
	test_and_set_nparams(2);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_3 d)
{
	test_and_set_nparams(3);
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_6 d)
{
	test_and_set_nparams(6);
	derivative_f = derivative_funcp(d);
	return *this;
}

function_options & function_options::expl_derivative_func(expl_derivative_funcp_1 d)
{
	test_and_set_nparams(1);
	expl_derivative_f = expl_derivative_funcp(d);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_2 d)
{
	test_and_set_nparams(2);
	expl_derivative_f = expl_derivative_funcp(d);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_3 d)
{
	test_and_set_nparams(3);
	expl_derivative_f = expl_derivative_funcp(d);
	return *this;
}

function_options & function_options::power_func(power_funcp_1 d)
{
	test_and_set_nparams(1);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_2 d)
{
	test_and_set_nparams(2);
	power_f = power_funcp(d);
	return *this;
}
function_options & function_options::power_func(power_funcp_3 d)
{
	test_and_set_nparams(3);
	power_f = power_funcp(d);
	return *this;
}

function_options & function_options::series_func(series_funcp_1 s)
{
	test_and_set_nparams(1);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_2 s)
{
	test_and_set_nparams(2);
	series_f = series_funcp(s);
	return *this;
}
function_options & function_options::series_func(series_funcp_3 s)
{
	test_and_set_nparams(3);
	series_f = series_funcp(s);
	return *this;
}

function_options & function_options::subs_func(subs_funcp_1 s)
{
	test_and_set_nparams(1);
	subs_f = subs_funcp(s);
	return *this;
}
function_options & function_options::subs_func(subs_funcp_2 s)
{
	test_and_set_nparams(2);
	subs_f = subs_funcp(s);
	return *this;
}
function_options & function_options::subs_func(subs_funcp_3 s)
{
	test_and_set_nparams(3);
	subs_f = subs_funcp(s);
	return *this;
}

// end of generated lines

function_options& function_options::eval_func(eval_funcp_exvector e)
{
	eval_use_exvector_args = true;
	eval_f = eval_funcp(e);
        pynac_eval_f = eval_f;
	return *this;
}
function_options& function_options::evalf_func(evalf_funcp_exvector ef)
{
	evalf_use_exvector_args = true;
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options& function_options::conjugate_func(conjugate_funcp_exvector c)
{
	conjugate_use_exvector_args = true;
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options& function_options::real_part_func(real_part_funcp_exvector c)
{
	real_part_use_exvector_args = true;
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options& function_options::imag_part_func(imag_part_funcp_exvector c)
{
	imag_part_use_exvector_args = true;
	imag_part_f = imag_part_funcp(c);
	return *this;
}

function_options& function_options::derivative_func(derivative_funcp_exvector d)
{
	derivative_use_exvector_args = true;
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options& function_options::power_func(power_funcp_exvector d)
{
	power_use_exvector_args = true;
	power_f = power_funcp(d);
	return *this;
}
function_options& function_options::series_func(series_funcp_exvector s)
{
	series_use_exvector_args = true;
	series_f = series_funcp(s);
	return *this;
}

function_options& function_options::derivative_func(
		derivative_funcp_exvector_symbol d)
{
	derivative_use_exvector_args = true;
	derivative_f = derivative_funcp(d);
	return *this;
}

function_options& function_options::eval_func(PyObject* e)
{
	python_func |= eval_python_f;
	eval_f = eval_funcp(e);
	return *this;
}
function_options& function_options::evalf_func(PyObject* ef)
{
	python_func |= evalf_python_f;
	evalf_f = evalf_funcp(ef);
	return *this;
}
function_options& function_options::conjugate_func(PyObject* c)
{
	python_func |= conjugate_python_f;
	conjugate_f = conjugate_funcp(c);
	return *this;
}
function_options& function_options::real_part_func(PyObject* c)
{
	python_func |= real_part_python_f;
	real_part_f = real_part_funcp(c);
	return *this;
}
function_options& function_options::imag_part_func(PyObject* c)
{
	python_func |= imag_part_python_f;
	imag_part_f = imag_part_funcp(c);
	return *this;
}

function_options& function_options::derivative_func(PyObject* d)
{
	python_func |= derivative_python_f;
	derivative_f = derivative_funcp(d);
	return *this;
}
function_options& function_options::power_func(PyObject* d)
{
	python_func |= power_python_f;
	power_f = power_funcp(d);
	return *this;
}
function_options& function_options::series_func(PyObject* s)
{
	python_func |= series_python_f;
	series_f = series_funcp(s);
	return *this;
}
function_options& function_options::subs_func(PyObject* e)
{
	python_func |= subs_python_f;
	subs_f = subs_funcp(e);
	return *this;
}

function_options & function_options::set_return_type(unsigned rt, tinfo_t rtt)
{
	use_return_type = true;
	return_type = rt;
	return_type_tinfo = rtt;
	return *this;
}

function_options & function_options::do_not_evalf_params()
{
	evalf_params_first = false;
	return *this;
}

function_options & function_options::do_not_apply_chain_rule()
{
	apply_chain_rule = false;
	return *this;
}

function_options & function_options::remember(unsigned size,
                                              unsigned assoc_size,
                                              unsigned strategy)
{
	use_remember = true;
	remember_size = size;
	remember_assoc_size = assoc_size;
	remember_strategy = strategy;
	return *this;
}

function_options & function_options::overloaded(unsigned o)
{
	functions_with_same_name = o;
	return *this;
}

void function_options::test_and_set_nparams(unsigned n)
{
	if (nparams==0) {
		nparams = n;
	} else if (nparams!=n) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// be caught anyhow
		std::cerr << "WARNING: " << name << "(): number of parameters ("
		          << n << ") differs from number set before (" 
		          << nparams << ")" << std::endl;
	}
}

void function_options::set_print_func(unsigned id, print_funcp f)
{
	if (id >= print_dispatch_table.size())
		print_dispatch_table.resize(id + 1);
	print_dispatch_table[id] = f;
}

void function_options::set_print_latex_func(PyObject* f)
{
	unsigned id = print_latex::get_class_info_static().options.get_id();
	if (id >= print_dispatch_table.size())
		print_dispatch_table.resize(id + 1);
	print_dispatch_table[id] = print_funcp(f);
}

void function_options::set_print_dflt_func(PyObject* f)
{
	unsigned id = print_dflt::get_class_info_static().options.get_id();
	if (id >= print_dispatch_table.size())
		print_dispatch_table.resize(id + 1);
	print_dispatch_table[id] = print_funcp(f);
}

/** This can be used as a hook for external applications. */
unsigned function::current_serial = 0;


registered_class_info function::reg_info = \
        registered_class_info(registered_class_options("function",
                                "exprseq",
                                &function::tinfo_static));

const tinfo_static_t function::tinfo_static = {};

//////////
// default constructor
//////////

// public

function::function() : serial(0)
{
	tinfo_key = &function::tinfo_static;
}

//////////
// other constructors
//////////

// public

function::function(unsigned ser) : serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

// the following lines have been generated for max. 14 parameters
function::function(unsigned ser, const ex & param1)
	: exprseq(param1), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2)
	: exprseq(param1, param2), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3)
	: exprseq(param1, param2, param3), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4)
	: exprseq(param1, param2, param3, param4), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5)
	: exprseq(param1, param2, param3, param4, param5), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6)
	: exprseq(param1, param2, param3, param4, param5, param6), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

// end of generated lines

function::function(unsigned ser, exprseq  es) : exprseq(std::move(es)), serial(ser)
{
	tinfo_key = &function::tinfo_static;

	// Force re-evaluation even if the exprseq was already evaluated
	// (the exprseq copy constructor copies the flags)
	clearflag(status_flags::evaluated);
}

function::function(unsigned ser, const exvector & v, bool discardable) 
  : exprseq(v,discardable), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

function::function(unsigned ser, std::unique_ptr<exvector> vp) 
  : exprseq(std::move(vp)), serial(ser)
{
	tinfo_key = &function::tinfo_static;
}

//////////
// archiving
//////////

/** Construct object from archive_node. */
function::function(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	// get python_func flag
	// since python_func used to be a bool flag in the old days,
	// in order to unarchive old format archives, we first check if the
	// archive node contains a bool property
	// if it doesn't we look for the unsigned property indicating which
	// custom functions are defined in python
	unsigned python_func;
	bool old_python_func;
	if (n.find_bool("python", old_python_func))
		python_func = old_python_func ?  0xFFFF : 0;
	else if(!n.find_unsigned("python", python_func))
		throw std::runtime_error("function::function archive error: cannot read python_func flag");
	std::string s;
	if (python_func != 0u) {
		// read the pickle from the archive
		if (!n.find_string("pickle", s))
			throw std::runtime_error("function::function archive error: cannot read pickled function");
		// unpickle
		PyObject* arg = Py_BuildValue("s#",s.c_str(), s.size());
		PyObject* sfunc = py_funcs.py_loads(arg);
		Py_DECREF(arg);
		if (PyErr_Occurred() != nullptr) {
		    throw(std::runtime_error("function::function archive error: caught exception in py_loads"));
		}
		// get the serial of the new SFunction
		unsigned int ser = py_funcs.py_get_serial_from_sfunction(sfunc);
		if (PyErr_Occurred() != nullptr) {
		    throw(std::runtime_error("function::function archive error: cannot get serial from SFunction"));
		}
		// set serial 
		serial = ser;
	} else { // otherwise
	// Find serial number by function name
	if (n.find_string("name", s)) {
		unsigned int ser = 0;
		unsigned int nargs = seq.size();
                for (const auto & elem : registered_functions()) {
			if (s == elem.name && nargs == elem.nparams) {
				serial = ser;
				return;
			}
			++ser;
		}
		// if the name is not already in the registry, we are
		// unarchiving a SymbolicFunction without any custom methods
		// Call Python to create a new symbolic function with name s
		// and get the serial of this new SymbolicFunction
		ser = py_funcs.py_get_serial_for_new_sfunction(s, nargs);
		if (PyErr_Occurred() != nullptr) {
		    throw(std::runtime_error("function::function archive error: cannot create new symbolic function " + s));
		}
		serial = ser;
		//throw (std::runtime_error("unknown function '" + s + "' in archive"));
	} else
		throw (std::runtime_error("unnamed function in archive"));
	}
}

/** Unarchive the object. */
ex function::unarchive(const archive_node &n, lst &sym_lst)
{
	return (new function(n, sym_lst))->setflag(status_flags::dynallocated);
}

/** Archive the object. */
void function::archive(archive_node &n) const
{
	inherited::archive(n);
	GINAC_ASSERT(serial < registered_functions().size());
	// we use Python's pickling mechanism to archive symbolic functions
	// with customized methods defined in Python. Symbolic functions
	// defined from c++ or those without custom methods are archived
	// directly, without calling Python. The python_func flag indicates if
	// we should use the python unpickling mechanism, or the regular
	// unarchiving for c++ functions.
	unsigned python_func = registered_functions()[serial].python_func;
	if (python_func != 0u) {
		n.add_unsigned("python", python_func);
		// find the corresponding SFunction object
		PyObject* sfunc = py_funcs.py_get_sfunction_from_serial(serial);
		if (PyErr_Occurred() != nullptr) {
		    throw(std::runtime_error("function::archive cannot get serial from SFunction"));
		}
		// call python to pickle it
		std::string* pickled = py_funcs.py_dumps(sfunc);
		if (PyErr_Occurred() != nullptr) {
		    throw(std::runtime_error("function::archive py_dumps raised exception"));
		}
		// store the pickle in the archive
		n.add_string("pickle", *pickled);
		delete pickled;
	} else {
		n.add_unsigned("python", 0);
		n.add_string("name", registered_functions()[serial].name);
	}
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void function::print(const print_context & c, unsigned level) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	// Dynamically dispatch on print_context type
	const print_context_class_info *pc_info = &c.get_class_info();
	if (serial >= static_cast<unsigned>(py_funcs.py_get_ginac_serial())) {
		//convert arguments to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);

		std::string* sout;
		if (is_a<print_latex>(c)) {
			sout = py_funcs.py_latex_function(serial, args);
                        if (PyErr_Occurred() != nullptr) {
                                throw(std::runtime_error("function::print(): python print function raised exception"));
                        }
                        c.s << *sout;
                        c.s.flush();
		}
		else if (is_a<print_tree>(c)) {
			sout = py_funcs.py_print_function(serial, args);
                        if (PyErr_Occurred() != nullptr) {
                                throw(std::runtime_error("function::print(): python print function raised exception"));
                        }
                        std::string fname = sout->substr(0, sout->find_first_of('('));
			c.s << std::string(level, ' ') << class_name() << " "
			    << fname << " @" << this << ", serial=" << serial
			    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
			    << ", nops=" << nops()
			    << std::endl;
			unsigned delta_indent = dynamic_cast<const print_tree &>(c).delta_indent;
			for (const auto& term : seq)
				term.print(c, level + delta_indent);
			c.s << std::string(level + delta_indent, ' ') << "=====" << std::endl;
		}
                else {
			sout = py_funcs.py_print_function(serial, args);
                        if (PyErr_Occurred() != nullptr) {
                                throw(std::runtime_error("function::print(): python print function raised exception"));
                        }
                        c.s << *sout;
                        c.s.flush();
		}

		delete sout;
		Py_DECREF(args);
	} else {

		if (is_a<print_latex>(c)) {
		        PyObject* sfunc = py_funcs.py_get_sfunction_from_serial(serial);
                        if (PyObject_HasAttrString(sfunc, "_print_latex_")) {
                                PyObject* args = py_funcs.exvector_to_PyTuple(seq);
                                std::string* sout;
                                sout = py_funcs.py_latex_function(serial, args);
                                if (PyErr_Occurred() != nullptr) {
                                        throw(std::runtime_error("function::print(): python print function raised exception"));
                                }
                                c.s << *sout;
                                c.s.flush();
                                delete sout;
                                Py_DECREF(args);
                                return;
                        }
		}

                const function_options &opt = registered_functions()[serial];
		const std::vector<print_funcp> &pdt = opt.print_dispatch_table;


next_context:
	unsigned id = pc_info->options.get_id();
	if (id >= pdt.size() || pdt[id] == nullptr) {

		// Method not found, try parent print_context class
		const print_context_class_info *parent_pc_info = pc_info->get_parent();
		if (parent_pc_info != nullptr) {
			pc_info = parent_pc_info;
			goto next_context;
		}

		// Method still not found, use default output
		if (is_a<print_tree>(c)) {

			c.s << std::string(level, ' ') << class_name() << " "
			    << opt.name << " @" << this << ", serial=" << serial
			    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
			    << ", nops=" << nops()
			    << std::endl;
			unsigned delta_indent = dynamic_cast<const print_tree &>(c).delta_indent;
			for (auto & elem : seq)
				elem.print(c, level + delta_indent);
			c.s << std::string(level + delta_indent, ' ') << "=====" << std::endl;

		} else if (is_a<print_latex>(c)) {
			c.s << opt.TeX_name;
			printseq(c, "\\left(", ',', "\\right)", exprseq::precedence(), function::precedence());
		} else {
			c.s << opt.name;
			printseq(c, "(", ',', ")", exprseq::precedence(), function::precedence());
		}

	} else {

		// Method found, call it
		current_serial = serial;
                if (opt.print_use_exvector_args)
                        (reinterpret_cast<print_funcp_exvector>(pdt[id]))(seq, c);
                else
                        switch (opt.nparams) {
                        // the following lines have been generated for max. 14 parameters
                        case 1:
                                (reinterpret_cast<print_funcp_1>(pdt[id]))(seq[1 - 1], c);
                                break;
                        case 2:
                                (reinterpret_cast<print_funcp_2>(pdt[id]))(seq[1 - 1], seq[2 - 1], c);
                                break;
                        case 3:
                                (reinterpret_cast<print_funcp_3>(pdt[id]))(seq[1 - 1], seq[2 - 1], seq[3 - 1], c);
                                break;

                        // end of generated lines
                        default:
                                throw(std::logic_error("function::print(): invalid nparams"));
                        }
        }
	}
}

ex function::expand(unsigned options) const
{
	return inherited::expand(options);
}

ex function::eval(int level) const
{
	if (level>1) {
		// first evaluate children, then we will end up here again
		return function(serial,evalchildren(level));
	}

	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	bool use_remember = opt.use_remember;
	ex eval_result;
	if (use_remember && lookup_remember_table(eval_result)) {
		return eval_result;
	}
	current_serial = serial;

	if (opt.eval_f==nullptr)
		return this->hold();
        eval_funcp eval_f;
        if (opt.pynac_eval_f == nullptr)
                eval_f = opt.eval_f;
        else
                eval_f = opt.pynac_eval_f;

	if (opt.pynac_eval_f == nullptr
            and (opt.python_func & function_options::eval_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);
		// call opt.eval_f with this list
		PyObject* pyresult = PyObject_CallMethod(reinterpret_cast<PyObject*>(eval_f),
				const_cast<char*>("_eval_"), const_cast<char*>("O"), args);
		Py_DECREF(args);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::eval(): python function raised exception"));
		}
		if ( pyresult == Py_None ) {
			return this->hold();
		}
		// convert output Expression to an ex
		eval_result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::eval(): python function (Expression_to_ex) raised exception"));
		}
	}
	else if (opt.eval_use_exvector_args)
		eval_result = (reinterpret_cast<eval_funcp_exvector>(eval_f))(seq);
	else
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		eval_result = (reinterpret_cast<eval_funcp_1>(eval_f))(seq[1-1]);
		break;
	case 2:
		eval_result = (reinterpret_cast<eval_funcp_2>(eval_f))(seq[1-1], seq[2-1]);
		break;
	case 3:
		eval_result = (reinterpret_cast<eval_funcp_3>(eval_f))(seq[1-1], seq[2-1], seq[3-1]);
		break;
	case 6:
		eval_result = (reinterpret_cast<eval_funcp_6>(eval_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1]);
		break;

		// end of generated lines
	default:
		throw(std::logic_error("function::eval(): invalid nparams"));
	}
	if (use_remember) {
		store_remember_table(eval_result);
	}
	return eval_result;
}

ex function::evalf(int level, PyObject* kwds) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	// Evaluate children first?
	exvector eseq;
	if (level == 1 || !(opt.evalf_params_first))
		eseq = seq;
	else if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	else {
		eseq.reserve(seq.size());
		--level;
                for (const auto & elem : seq)
			eseq.push_back(elem.evalf(level, kwds));
	}

	if (opt.evalf_f == nullptr) {
                if (opt.nparams == 1 and is_exactly_a<numeric>(eseq[1-1])) {
                        const numeric& n = ex_to<numeric>(eseq[1-1]);
                        try {
                                return n.try_py_method(get_name());
                        }
                        catch (std::logic_error) {
                                try {
                                        const numeric& nn = ex_to<numeric>(n.evalf()).try_py_method(get_name());
                                        return nn.to_dict_parent(kwds);
                                }
                                catch (std::logic_error) {}
                        }
                }
		return function(serial,eseq).hold();
	}
	current_serial = serial;
	if ((opt.python_func & function_options::evalf_python_f) != 0u) { 
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(eseq);
		// call opt.evalf_f with this list
		PyObject* pyresult = PyObject_Call(
			PyObject_GetAttrString(reinterpret_cast<PyObject*>(opt.evalf_f),
				"_evalf_"), args, kwds);
		Py_DECREF(args);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::evalf(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::evalf(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.evalf_use_exvector_args)
		return (reinterpret_cast<evalf_funcp_exvector>(opt.evalf_f))(seq, kwds);
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return (reinterpret_cast<evalf_funcp_1>(opt.evalf_f))(eseq[1-1], kwds);
	case 2:
		return (reinterpret_cast<evalf_funcp_2>(opt.evalf_f))(eseq[1-1], eseq[2-1], kwds);
	case 3:
		return (reinterpret_cast<evalf_funcp_3>(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], kwds);
	case 6:
		return (reinterpret_cast<evalf_funcp_6>(opt.evalf_f))(eseq[1-1], eseq[2-1], eseq[3-1], eseq[4-1], eseq[5-1], eseq[6-1], kwds);

		// end of generated lines
	}
	throw(std::logic_error("function::evalf(): invalid nparams"));
}

long function::calchash() const
{
	long v = golden_ratio_hash(golden_ratio_hash((intptr_t)tinfo()) ^ serial);
	for (size_t i=0; i<nops(); i++) {
		v = rotate_left(v);
		v ^= this->op(i).gethash();
	}

	if (is_evaluated()) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}
	return v;
}

ex function::thiscontainer(const exvector & v) const
{
	return function(serial, v);
}

ex function::thiscontainer(std::unique_ptr<exvector> vp) const
{
	return function(serial, std::move(vp));
}

/** Implementation of ex::series for functions.
 *  @see ex::series */
ex function::series(const relational & r, int order, unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.series_f==nullptr) {
		return basic::series(r, order);
	}
	ex res;
	current_serial = serial;
	if ((opt.python_func & function_options::series_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);
		// create a dictionary {'order': order, 'options':options}
		PyObject* kwds = Py_BuildValue("{s:i,s:I}","order",order,"options",options);
		// add variable to expand for as a keyword argument
		PyDict_SetItemString(kwds, "var", py_funcs.ex_to_pyExpression(r.lhs()));
		// add the point of expansion as a keyword argument
		PyDict_SetItemString(kwds, "at", py_funcs.ex_to_pyExpression(r.rhs()));
		// call opt.series_f with this list
		PyObject* pyresult = PyObject_Call(
			PyObject_GetAttrString(reinterpret_cast<PyObject*>(opt.series_f),
				"_series_"), args, kwds);
		Py_DECREF(args);
		Py_DECREF(kwds);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::series(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::series(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.series_use_exvector_args) {
		try {
			res = (reinterpret_cast<series_funcp_exvector>(opt.series_f))(seq, r, order, options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	}
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		try {
			res = (reinterpret_cast<series_funcp_1>(opt.series_f))(seq[1-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 2:
		try {
			res = (reinterpret_cast<series_funcp_2>(opt.series_f))(seq[1-1], seq[2-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	case 3:
		try {
			res = (reinterpret_cast<series_funcp_3>(opt.series_f))(seq[1-1], seq[2-1], seq[3-1],r,order,options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;

		// end of generated lines
	}
	throw(std::logic_error("function::series(): invalid nparams"));
}


/** Implementation of ex::subs for functions. */
ex function::subs(const exmap & m, unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if ((opt.python_func & function_options::subs_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.subs_args_to_PyTuple(m, options, seq);
		// call opt.subs_f with this list
		PyObject* pyresult = PyObject_CallMethod(
				reinterpret_cast<PyObject*>(opt.subs_f),
				const_cast<char*>("_subs_"), const_cast<char*>("O"), args);
		Py_DECREF(args);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::subs(): python method (_subs_) raised exception"));
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::subs(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	} 
	if (opt.subs_f==nullptr)
        	return exprseq::subs(m, options);

	switch (opt.nparams) {
	case 1:
		return (reinterpret_cast<subs_funcp_1>(opt.subs_f))(m, seq[1-1]);
	case 2:
		return (reinterpret_cast<subs_funcp_2>(opt.subs_f))(m, seq[1-1], seq[2-1]);
	case 3:
		return (reinterpret_cast<subs_funcp_3>(opt.subs_f))(m, seq[1-1], seq[2-1], seq[3-1]);

		// end of generated lines
	}
	throw(std::logic_error("function::subs(): invalid nparams"));
}

/** Implementation of ex::conjugate for functions. */
ex function::conjugate() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.conjugate_f==nullptr) {
		return conjugate_function(*this).hold();
	}

	if ((opt.python_func & function_options::conjugate_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);
		// call opt.conjugate_f with this list
		PyObject* pyresult = PyObject_CallMethod(
				reinterpret_cast<PyObject*>(opt.conjugate_f),
				const_cast<char*>("_conjugate_"), const_cast<char*>("O"), args);
		Py_DECREF(args);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::conjugate(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::conjugate(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.conjugate_use_exvector_args) {
		return (reinterpret_cast<conjugate_funcp_exvector>(opt.conjugate_f))(seq);
	}

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return (reinterpret_cast<conjugate_funcp_1>(opt.conjugate_f))(seq[1-1]);
	case 2:
		return (reinterpret_cast<conjugate_funcp_2>(opt.conjugate_f))(seq[1-1], seq[2-1]);
	case 3:
		return (reinterpret_cast<conjugate_funcp_3>(opt.conjugate_f))(seq[1-1], seq[2-1], seq[3-1]);

		// end of generated lines
	}
	throw(std::logic_error("function::conjugate(): invalid nparams"));
}

/** Implementation of ex::real_part for functions. */
ex function::real_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.real_part_f==nullptr)
		return basic::real_part();

	if ((opt.python_func & function_options::real_part_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);
		// call opt.real_part_f with this list
		PyObject* pyresult = PyObject_CallMethod(reinterpret_cast<PyObject*>(opt.real_part_f),
				const_cast<char*>("_real_part_"), const_cast<char*>("O"), args);
		Py_DECREF(args);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::real_part(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::real_part(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.real_part_use_exvector_args)
		return (reinterpret_cast<real_part_funcp_exvector>(opt.real_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return (reinterpret_cast<real_part_funcp_1>(opt.real_part_f))(seq[1-1]);
	case 2:
		return (reinterpret_cast<real_part_funcp_2>(opt.real_part_f))(seq[1-1], seq[2-1]);
	case 3:
		return (reinterpret_cast<real_part_funcp_3>(opt.real_part_f))(seq[1-1], seq[2-1], seq[3-1]);

		// end of generated lines
	}
	throw(std::logic_error("function::real_part(): invalid nparams"));
}

/** Implementation of ex::imag_part for functions. */
ex function::imag_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.imag_part_f==nullptr)
		return basic::imag_part();

	if ((opt.python_func & function_options::imag_part_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);
		// call opt.imag_part_f with this list
		PyObject* pyresult = PyObject_CallMethod(reinterpret_cast<PyObject*>(opt.imag_part_f),
				const_cast<char*>("_imag_part_"), const_cast<char*>("O"), args);
		Py_DECREF(args);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::imag_part(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::imag_part(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.imag_part_use_exvector_args)
		return (reinterpret_cast<imag_part_funcp_exvector>(opt.imag_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return (reinterpret_cast<imag_part_funcp_1>(opt.imag_part_f))(seq[1-1]);
	case 2:
		return (reinterpret_cast<imag_part_funcp_2>(opt.imag_part_f))(seq[1-1], seq[2-1]);
	case 3:
		return (reinterpret_cast<imag_part_funcp_3>(opt.imag_part_f))(seq[1-1], seq[2-1], seq[3-1]);

		// end of generated lines
	}
	throw(std::logic_error("function::imag_part(): invalid nparams"));
}

// protected

/** Implementation of ex::diff() for functions. It applies the chain rule,
 *  except for the Order term function.
 *  @see ex::diff */
ex function::derivative(const symbol & s) const
{
        ex result;
        
	/*
	if (serial == Order_SERIAL::serial) {
		// Order Term function only differentiates the argument
		return Order(seq[0].diff(s));
		*/
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

        try {
                // Explicit derivation
                return expl_derivative(s);
        } catch (...) {}

	// Check if we need to apply chain rule
	if (!opt.apply_chain_rule) {
		if (opt.derivative_f == nullptr)
			throw(std::runtime_error("function::derivative(): custom derivative function must be defined"));

		if ((opt.python_func & function_options::derivative_python_f) != 0u) {
			// convert seq to a PyTuple of Expressions
			PyObject* args = py_funcs.exvector_to_PyTuple(seq);
			// create a dictionary {'diff_param': s}
			PyObject* symb = py_funcs.ex_to_pyExpression(s);
			PyObject* kwds = Py_BuildValue("{s:O}","diff_param",
					symb);
			// call opt.derivative_f with this list
			PyObject* pyresult = PyObject_Call(
				PyObject_GetAttrString(
					reinterpret_cast<PyObject*>(opt.derivative_f),
					"_tderivative_"), args, kwds);
			Py_DECREF(symb);
			Py_DECREF(args);
			Py_DECREF(kwds);
			if (pyresult == nullptr) { 
				throw(std::runtime_error("function::derivative(): python function raised exception"));
			}
			// convert output Expression to an ex
			result = py_funcs.pyExpression_to_ex(pyresult);
			Py_DECREF(pyresult);
			if (PyErr_Occurred() != nullptr) { 
				throw(std::runtime_error("function::derivative(): python function (pyExpression_to_ex) raised exception"));
			}
			return result;
		}
		// C++ function
		if (!opt.derivative_use_exvector_args)
			throw(std::runtime_error("function::derivative(): cannot call C++ function without exvector args"));
		
		return (reinterpret_cast<derivative_funcp_exvector_symbol>(opt.derivative_f))(seq, s);

	} 
        // Chain rule
        ex arg_diff;
        size_t num = seq.size();
        for (size_t i=0; i<num; i++) {
                arg_diff = seq[i].diff(s);
                // We apply the chain rule only when it makes sense.  This is not
                // just for performance reasons but also to allow functions to
                // throw when differentiated with respect to one of its arguments
                // without running into trouble with our automatic full
                // differentiation:
                if (!arg_diff.is_zero())
                        result += pderivative(i)*arg_diff;
        }

	return result;
}

int function::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return serial < o.serial ? -1 : 1;

        return exprseq::compare_same_type(o);
}


bool function::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return false;

	return exprseq::is_equal_same_type(o);
}

bool function::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	return serial == o.serial;
}

bool function::match(const ex & pattern, exmap& map) const
{
	if (is_exactly_a<wildcard>(pattern)) {
                const auto& it = map.find(pattern);
                if (it != map.end())
		        return is_equal(ex_to<basic>(it->second));
		map[pattern] = *this;
		return true;
	} 
        if (not is_exactly_a<function>(pattern))
                return false;
        CMatcher cm(*this, pattern, map);
        const opt_exmap& m = cm.get();
        if (not m)
                return false;
        map = m.value();
        return true;
}


unsigned function::return_type() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type;
	} 
        // Default behavior is to use the return type of the first
        // argument. Thus, exp() of a matrix behaves like a matrix, etc.
        if (seq.empty())
                return return_types::commutative;

        return seq.begin()->return_type();
	
}

tinfo_t function::return_type_tinfo() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type_tinfo;
	} 
		// Default behavior is to use the return type of the first
		// argument. Thus, exp() of a matrix behaves like a matrix, etc.
		if (seq.empty())
			return this;
		
			return seq.begin()->return_type_tinfo();
	
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// protected

ex function::pderivative(unsigned diff_param) const // partial differentiation
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	// No derivative defined? Then return abstract derivative object
	if (opt.derivative_f == nullptr)
		return fderivative(serial, diff_param, seq);

	current_serial = serial;
	if ((opt.python_func & function_options::derivative_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);
		// create a dictionary {'diff_param': diff_param}
		PyObject* kwds = Py_BuildValue("{s:I}","diff_param",diff_param);
		// call opt.derivative_f with this list
		PyObject* pyresult = PyObject_Call(
			PyObject_GetAttrString(reinterpret_cast<PyObject*>(opt.derivative_f),
				"_derivative_"), args, kwds);
		Py_DECREF(args);
		Py_DECREF(kwds);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::pderivative(): python function raised exception"));
		}
		if ( pyresult == Py_None ) {
			return fderivative(serial, diff_param, seq);
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::pderivative(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.derivative_use_exvector_args)
		return (reinterpret_cast<derivative_funcp_exvector>(opt.derivative_f))(seq, diff_param);
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return (reinterpret_cast<derivative_funcp_1>(opt.derivative_f))(seq[1-1],diff_param);
	case 2:
		return (reinterpret_cast<derivative_funcp_2>(opt.derivative_f))(seq[1-1], seq[2-1],diff_param);
	case 3:
		return (reinterpret_cast<derivative_funcp_3>(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1],diff_param);
	case 6:
		return (reinterpret_cast<derivative_funcp_6>(opt.derivative_f))(seq[1-1], seq[2-1], seq[3-1], seq[4-1], seq[5-1], seq[6-1], diff_param);

		// end of generated lines
	}
	throw(std::logic_error("function::pderivative(): no diff function defined"));
}

ex function::expl_derivative(const symbol & s) const // explicit differentiation
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

        if (opt.expl_derivative_f) {
		// Invoke the defined explicit derivative function.
		current_serial = serial;
		if (opt.expl_derivative_use_exvector_args)
			return (reinterpret_cast<expl_derivative_funcp_exvector>(opt.expl_derivative_f))(seq, s);
		switch (opt.nparams) {
			// the following lines have been generated for max. 14 parameters
			case 1:
				return (reinterpret_cast<expl_derivative_funcp_1>(opt.expl_derivative_f))(seq[0], s);
			case 2:
				return (reinterpret_cast<expl_derivative_funcp_2>(opt.expl_derivative_f))(seq[0], seq[1], s);
			case 3:
				return (reinterpret_cast<expl_derivative_funcp_3>(opt.expl_derivative_f))(seq[0], seq[1], seq[2], s);
		}
	}
	// There is no fallback for explicit derivative.
	throw(std::logic_error("function::expl_derivative(): explicit derivation is called, but no such function defined"));
}

ex function::power(const ex & power_param) const // power of function
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	// No derivative defined? Then return abstract derivative object
	if (opt.power_f == nullptr)
		return (new GiNaC::power(*this, power_param))->setflag(status_flags::dynallocated |
	                                               status_flags::evaluated);

	current_serial = serial;
	if ((opt.python_func & function_options::power_python_f) != 0u) {
		// convert seq to a PyTuple of Expressions
		PyObject* args = py_funcs.exvector_to_PyTuple(seq);
		// create a dictionary {'power_param': power_param}
		PyObject* kwds = PyDict_New();
		PyDict_SetItemString(kwds, "power_param", py_funcs.ex_to_pyExpression(power_param));
		// call opt.power_f with this list
		PyObject* pyresult = PyObject_Call(
			PyObject_GetAttrString(reinterpret_cast<PyObject*>(opt.power_f),
				"_power_"), args, kwds);
		Py_DECREF(args);
		Py_DECREF(kwds);
		if (pyresult == nullptr) { 
			throw(std::runtime_error("function::power(): python function raised exception"));
		}
		// convert output Expression to an ex
		ex result = py_funcs.pyExpression_to_ex(pyresult);
		Py_DECREF(pyresult);
		if (PyErr_Occurred() != nullptr) { 
			throw(std::runtime_error("function::power(): python function (pyExpression_to_ex) raised exception"));
		}
		return result;
	}
	if (opt.power_use_exvector_args)
		return (reinterpret_cast<power_funcp_exvector>(opt.power_f))(seq,  power_param);
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
	case 1:
		return (reinterpret_cast<power_funcp_1>(opt.power_f))(seq[1-1],power_param);
	case 2:
		return (reinterpret_cast<power_funcp_2>(opt.power_f))(seq[1-1], seq[2-1],power_param);
	case 3:
		return (reinterpret_cast<power_funcp_3>(opt.power_f))(seq[1-1], seq[2-1], seq[3-1],power_param);

		// end of generated lines
	}
	throw(std::logic_error("function::power(): no power function defined"));
}

std::vector<function_options> & function::registered_functions()
{
	static auto  rf = new std::vector<function_options>;
	return *rf;
}

bool function::lookup_remember_table(ex & result) const
{
	return remember_table::remember_tables()[this->serial].lookup_entry(*this,result);
}

void function::store_remember_table(ex const & result) const
{
	remember_table::remember_tables()[this->serial].add_entry(*this,result);
}

// public

unsigned function::register_new(function_options const & opt)
{
	size_t same_name = 0;
	for (auto & elem : registered_functions()) {
		if (elem.name==opt.name) {
			++same_name;
		}
	}
	if (same_name>=opt.functions_with_same_name) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// caught anyhow
		//
		// SAGE note: 
		// We suppress this warning since we allow a user to create 
		// functions with same name, but different parameters
		// Sage SFunction class checks existence of a function before
		// allocating a new one.
		//
		//std::cerr << "WARNING: function name " << opt.name
		//          << " already in use!" << std::endl;
	}
	registered_functions().push_back(opt);
	if (opt.use_remember) {
		remember_table::remember_tables().
			emplace_back(opt.remember_size,
			                         opt.remember_assoc_size,
			                         opt.remember_strategy);
	} else {
		remember_table::remember_tables().emplace_back();
	}
	return registered_functions().size()-1;
}

/** Find serial number of function by name and number of parameters.
 *  Throws exception if function was not found. */
unsigned function::find_function(const std::string &name, unsigned nparams)
{
	unsigned serial = 0;
        for (const auto & elem : registered_functions()) {
		if (elem.get_name() == name && elem.get_nparams() == nparams)
			return serial;
		++serial;
	}
	throw (std::runtime_error("no function '" + name + "' with " + ToString(nparams) + " parameters defined"));
}

/** Return the print name of the function. */
std::string function::get_name() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	return registered_functions()[serial].name;
}

void function::set_domain(unsigned d)
{
        domain = d;
        iflags.clear();
        switch (d) {
        case domain::complex:
                break;
        case domain::real:
                iflags.set(info_flags::real, true);
                break;
        case domain::positive:
                iflags.set(info_flags::real, true);
                iflags.set(info_flags::positive, true);
                break;
        case domain::integer:
                iflags.set(info_flags::real, true);
                iflags.set(info_flags::integer, true);
                break;
        }
}

bool has_function(const ex & x)
{
	if (is_exactly_a<function>(x))
		return true;
	for (size_t i=0; i<x.nops(); ++i)
		if (has_function(x.op(i)))
			return true;

	return false;
}

bool has_symbol_or_function(const ex & x)
{
	if (is_exactly_a<symbol>(x) or is_exactly_a<function>(x))
		return true;
	for (size_t i=0; i<x.nops(); ++i)
		if (has_symbol_or_function(x.op(i)))
			return true;

	return false;
}

static bool has_oneof_function_helper(const ex& x,
                const std::map<unsigned,int>& m)
{
	if (is_exactly_a<function>(x)
            and m.find(ex_to<function>(x).get_serial()) != m.end())
		return true;
	for (size_t i=0; i<x.nops(); ++i)
		if (has_oneof_function_helper(x.op(i), m))
			return true;

	return false;
}

static void has_allof_function_helper(const ex& x,
                std::map<unsigned,int>& m)
{
	if (is_exactly_a<function>(x)) {
                unsigned ser = ex_to<function>(x).get_serial();
                if (m.find(ser) != m.end())
        		m[ser] = 1;
        }
	for (size_t i=0; i<x.nops(); ++i)
		has_allof_function_helper(x.op(i), m);
}

bool has_function(const ex& x,
                const std::string& s)
{
        std::map<unsigned,int> m;
        unsigned ser = 0;
        for (const auto & elem : function::registered_functions()) {
                if (s == elem.name)
                        m[ser] = 0;
                ++ser;
        }
        if (m.empty())
                return false;
        return has_oneof_function_helper(x, m);
}

bool has_function(const ex& x,
                const std::vector<std::string>& v,
                bool all)
{
        std::map<unsigned,int> m;
        for (const auto & s : v) {
                unsigned ser = 0;
                for (const auto & elem : function::registered_functions()) {
                        if (s == elem.name)
                                m[ser] = 0;
                        ++ser;
                }
        }
        if (m.empty())
                return false;
        if (all) {
                has_allof_function_helper(x, m);
                for (const auto & p : m)
                        // TODO: false negative if >1 func with same name
                        if (p.second == 0)
                                return false;
                return true;
        }
        return has_oneof_function_helper(x, m);
}

} // namespace GiNaC

