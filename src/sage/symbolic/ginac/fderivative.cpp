/** @file fderivative.cpp
 *
 *  Implementation of abstract derivatives of functions. */

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
#include <Python.h>
#include "py_funcs.h"
#include "fderivative.h"
#include "operators.h"
#include "archive.h"
#include "utils.h"

#include <iostream>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(fderivative, function,
  print_func<print_context>(&fderivative::do_print).
  print_func<print_tree>(&fderivative::do_print_tree))

//////////
// default constructor
//////////

fderivative::fderivative()
{
	tinfo_key = &fderivative::tinfo_static;
}

//////////
// other constructors
//////////

fderivative::fderivative(unsigned ser, unsigned param, const exvector & args) : function(ser, args)
{
	parameter_set.insert(param);
	tinfo_key = &fderivative::tinfo_static;
}

fderivative::fderivative(unsigned ser, paramset  params, const exvector & args) : function(ser, args), parameter_set(std::move(params))
{
	tinfo_key = &fderivative::tinfo_static;
}

fderivative::fderivative(unsigned ser, paramset  params, std::unique_ptr<exvector> vp) : function(ser, std::move(vp)), parameter_set(std::move(params))
{
	tinfo_key = &fderivative::tinfo_static;
}

//////////
// archiving
//////////

fderivative::fderivative(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	unsigned i = 0;
	while (true) {
		unsigned u;
		if (n.find_unsigned("param", u, i))
			parameter_set.insert(u);
		else
			break;
		++i;
	}
}

void fderivative::archive(archive_node &n) const
{
	inherited::archive(n);
        for (const auto & elem : parameter_set)
		n.add_unsigned("param", elem);
}

//////////
// functions overriding virtual functions from base classes
//////////

void fderivative::print(const print_context & c, unsigned level) const
{
	// class function overrides print(), but we don't want that
	basic::print(c, level);
}

void fderivative::do_print(const print_context & c, unsigned /*unused*/) const
{
	//convert paramset to a python list
	PyObject* params = py_funcs.paramset_to_PyTuple(parameter_set);
	//convert arguments to a PyTuple of Expressions
	PyObject* args = py_funcs.exvector_to_PyTuple(seq);
	//check if latex mode
	//call python function
	std::string *sout;
	if (is_a<print_latex>(c)) {
		sout = py_funcs.py_latex_fderivative(serial, params, args);
	} else {
		sout = py_funcs.py_print_fderivative(serial, params, args);
	}
	if (sout == nullptr) { 
		throw(std::runtime_error("fderivative::do_print(): python print function raised exception"));
	}
	c.s<<*sout;
	delete sout;
	Py_DECREF(params);
	Py_DECREF(args);
	
	/*
	c.s << "D[";
	paramset::const_iterator i = parameter_set.begin(), end = parameter_set.end();
	--end;
	while (i != end) {
		c.s << *i++ << ",";
	}
	c.s << *i << "](" << registered_functions()[serial].name << ")";
	printseq(c, "(", ',', ")", exprseq::precedence(), function::precedence());
	*/
}

void fderivative::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " "
	    << registered_functions()[serial].name << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", nops=" << nops()
	    << ", params=";
	auto i = parameter_set.begin(), iend = parameter_set.end();
	--iend;
	while (i != iend)
		c.s << *i++ << ",";
	c.s << *i << std::endl;
	for (auto & elem : seq)
		elem.print(c, level + c.delta_indent);
	c.s << std::string(level + c.delta_indent, ' ') << "=====" << std::endl;
}

ex fderivative::eval(int level) const
{
	if (level > 1) {
		// first evaluate children, then we will end up here again
		return fderivative(serial, parameter_set, evalchildren(level));
	}

	// No parameters specified? Then return the function itself
	if (parameter_set.empty())
		return function(serial, seq);

	// If the function in question actually has a derivative, return it
	if (registered_functions()[serial].has_derivative() && parameter_set.size() == 1)
		return pderivative(*(parameter_set.begin()));

	return this->hold();
}

/** Numeric evaluation falls back to evaluation of arguments.
 *  @see basic::evalf */
ex fderivative::evalf(int level, PyObject* parent) const
{
	return basic::evalf(level, parent);
}

/** The series expansion of derivatives falls back to Taylor expansion.
 *  @see basic::series */
ex fderivative::series(const relational & r, int order, unsigned options) const
{
	return basic::series(r, order, options);
}

ex fderivative::thiscontainer(const exvector & v) const
{
	return fderivative(serial, parameter_set, v);
}

ex fderivative::thiscontainer(std::unique_ptr<exvector> vp) const
{
	return fderivative(serial, parameter_set, std::move(vp));
}

/** Implementation of ex::diff() for derivatives. It applies the chain rule.
 *  @see ex::diff */
ex fderivative::derivative(const symbol & s) const
{
	ex result;
	for (size_t i=0; i<seq.size(); i++) {
		ex arg_diff = seq[i].diff(s);
		if (!arg_diff.is_zero()) {
			paramset ps = parameter_set;
			ps.insert(i);
			result += arg_diff * fderivative(serial, ps, seq);
		}
	}
	return result;
}

int fderivative::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<fderivative>(other));
	const fderivative & o = static_cast<const fderivative &>(other);

	if (parameter_set != o.parameter_set)
		return parameter_set < o.parameter_set ? -1 : 1;

	return inherited::compare_same_type(o);
}

bool fderivative::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<fderivative>(other));
	const fderivative & o = static_cast<const fderivative &>(other);

	if (parameter_set != o.parameter_set)
		return false;

        return inherited::is_equal_same_type(o);
}

bool fderivative::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<fderivative>(other));
	const fderivative & o = static_cast<const fderivative &>(other);

	return parameter_set == o.parameter_set && inherited::match_same_type(other);
}

long fderivative::calchash() const
{
	long res = inherited::calchash();
	long prev=0;

	// FNV hash with custom prime and initial value
	long h = 2166136285L;
	for (const auto & elem : parameter_set) {
		h = ( h * 2097287 ) ^ (elem-prev);
		prev = elem;
	}

	res=h^res;
	if (is_evaluated()) {
		setflag(status_flags::hash_calculated);
		hashvalue = res;
	}
	return res;
}

} // namespace GiNaC
