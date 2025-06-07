/** @file infinity.cpp
 *
 *  Implementation of PyNaC's "infinity". */

/*
 *  Copyright (C) 2011 Volker Braun <vbraun@stp.dias.ie>
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

#include "py_funcs.h"
#include "infinity.h"
#include "numeric.h"
#include "ex.h"
#include "archive.h"
#include "utils.h"
#include "add.h"
#include "mul.h"
#include "inifcns.h"

#include <climits>
#include <string>
#include <stdexcept>
#include <iostream>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(infinity, basic,
  print_func<print_context>(&infinity::do_print).
  print_func<print_latex>(&infinity::do_print_latex).
  print_func<print_tree>(&infinity::do_print_tree).
  print_func<print_python_repr>(&infinity::do_print_python_repr))

static long hash_from_dir(const ex& direction)
{
        if (direction.is_one())
                return LONG_MAX;
        if (direction.is_zero())
                return LONG_MAX-1;
        if (direction.is_minus_one())
                return LONG_MIN;
        return 0L;
}

//////////
// default constructor
//////////

// public

infinity::infinity() 
: basic(&infinity::tinfo_static), direction(+1)
{
        hashvalue = hash_from_dir(direction);
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

// public

infinity::infinity(const numeric & _direction)
: basic(&infinity::tinfo_static)
{
	// Note: we cannot accept an arbitrary ex as argument 
	// or we would take precedence over the copy constructor.
	set_direction(_direction);
        hashvalue = hash_from_dir(direction);
	setflag(status_flags::evaluated | status_flags::expanded);
}

infinity infinity::from_direction(const ex & _direction)
{
	GINAC_ASSERT(!is_a<infinity>(_direction));
	infinity result;
	result.set_direction(_direction);
	result.hashvalue = hash_from_dir(result.direction);
        return result;
}

infinity infinity::from_sign(int sgn)
{
	GINAC_ASSERT(sgn>=-1 and sgn<=1);
	infinity result;
	result.direction = sgn;
	result.hashvalue = hash_from_dir(result.direction);
	return result;
}

//////////
// archiving
//////////

infinity::infinity(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst) {}

ex infinity::unarchive(const archive_node &n, lst &sym_lst)
{
	ex value;
	if (n.find_ex("direction", value, sym_lst))
		return infinity::from_direction(value);

        throw(std::runtime_error("infinity without direction in archive"));
}

void infinity::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_ex("direction", direction);
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void infinity::do_print(const print_context & c, unsigned level) const
{
	if (is_unsigned_infinity())
		c.s << "Infinity";
	else if (is_plus_infinity())
		c.s << "+Infinity";
	else if (is_minus_infinity())
		c.s << "-Infinity";
	else {
		c.s << "(";
		direction.print(c, level);
		c.s << ")*Infinity";
	}
}

void infinity::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ');
	if (is_unsigned_infinity())
		c.s << "unsigned_infinity";
	else
		c.s << "infinity";
	c.s << " (" << class_name() << ")" << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << std::endl;
	if (!is_unsigned_infinity()) {
                c.s << "with direction: ";
		direction.print(c, level+4);
        }
}

void infinity::do_print_latex(const print_latex & c, unsigned level) const
{
	if (is_unsigned_infinity())
		c.s << "\\infty";
	else if (is_plus_infinity())
		c.s << "+\\infty";
	else if (is_minus_infinity())
		c.s << "-\\infty";
	else {
		c.s << "(";
		direction.print(c, level);
		c.s << ") \\infty";
	}
}

void infinity::do_print_python_repr(const print_python_repr & c,
                unsigned level) const
{
	c.s << class_name() << "('" << "Infinity" << "'";
        direction.print(c, level);
	c.s << ')';
}

bool infinity::info(unsigned inf) const
{
        switch (inf) {
        case info_flags::infinity:
                return true;
        case info_flags::real:
                return not direction.is_zero()
                and direction.is_real();
        case info_flags::positive:
        case info_flags::negative:
                return direction.info(inf);
        case info_flags::nonnegative:
                return direction.is_positive();
        default:
                return inherited::info(inf);
        }
}

ex infinity::evalf(int /*level*/, PyObject* /*parent*/) const
{
	if (is_unsigned_infinity())
		return py_funcs.py_eval_unsigned_infinity();
	if (is_plus_infinity())
		return py_funcs.py_eval_infinity();
	if (is_minus_infinity())
		return py_funcs.py_eval_neg_infinity();
	return *this;
}

ex infinity::conjugate() const
{
	return infinity::from_direction(direction.conjugate());
}

ex infinity::real_part() const
{
	if (is_unsigned_infinity()) 
		throw(std::runtime_error("indeterminate expression: "
					 "real part of unsigned_infinity."));
	ex re_dir = direction.real_part();
	if (re_dir.is_zero()) return _ex0;
	return infinity::from_direction(re_dir);
}

ex infinity::imag_part() const
{
	if (is_unsigned_infinity()) 
		throw(std::runtime_error("indeterminate expression: "
					 "imaginary part of unsigned_infinity."));
	ex im_dir = direction.imag_part();
	if (im_dir.is_zero()) return _ex0;
	return infinity::from_direction(im_dir);
}

// protected

ex infinity::derivative(const symbol & /*s*/) const
{
	return _ex0;
}

int infinity::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<infinity>(other));
	const infinity &o = static_cast<const infinity &>(other);

	return direction.compare(o.direction);
}

bool infinity::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<infinity>(other));
	const infinity &o = static_cast<const infinity &>(other);
	return direction.is_equal(o.direction);
}

bool infinity::compare_other_type(const ex & other,
        relational::operators o) const
{
        const ex& e = other.evalf();
        if (not is_exactly_a<numeric>(e))
                return false;
        const numeric& num = ex_to<numeric>(e);
        if (num.imag() > 0)
                return false;
        switch (o) {
        case relational::not_equal:
                return true;
        case relational::equal:
                return false;
        case relational::less_or_equal:
        case relational::less:
                return is_minus_infinity();
        default:
                return is_plus_infinity();
        }
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

bool infinity::is_unsigned_infinity() const 
{ 
	return direction.is_zero(); 
}


bool infinity::is_plus_infinity() const 
{ 
	return direction.is_one(); 
}


bool infinity::is_minus_infinity() const 
{ 
	return direction.is_minus_one(); 
}


void infinity::set_direction(const ex & new_direction)
{
	if (new_direction.is_zero())
		direction = _ex0;
	else {
		ex normalization = GiNaC::pow(GiNaC::abs(new_direction),-1);
		direction = mul(new_direction, normalization);
	}
        hashvalue = hash_from_dir(direction);
}

const infinity & infinity::operator *= (const ex & rhs)
{
	if (is_exactly_a<infinity>(rhs)) {
		const ex & rhs_direction = ex_to<infinity>(rhs).direction;
		set_direction(mul(direction, rhs_direction));
		return *this;
	}
	if (rhs.is_zero())
		throw(std::runtime_error("indeterminate expression: "
					 "0 * infinity encountered."));
	else if (rhs.is_positive()) {
		return *this;
	}
        if (rhs.info(info_flags::negative)) {
		set_direction(mul(-1, direction));
		return *this;
	}
        if (rhs.nsymbols()==0) {
		set_direction(mul(direction, rhs));
		return *this;
	}
	throw(std::runtime_error("indeterminate expression: "
				 "infinity * f(x) encountered."));
}


const infinity & infinity::operator += (const ex & rhs)
{
	if (not is_exactly_a<infinity>(rhs))
		return *this;
	const ex & rhs_direction = ex_to<infinity>(rhs).direction;
	if (not direction.is_equal(rhs_direction) ) {
		if (ex_to<infinity>(rhs).is_unsigned_infinity() or is_unsigned_infinity())
			throw(std::runtime_error("indeterminate expression: "
						 "unsigned_infinity +- infinity encountered."));
		else
			throw(std::runtime_error("indeterminate expression: "
						 "infinity - infinity encountered."));
        }
	return *this;
}


//////////
// global "infinity" constants
//////////

/** Infinity, i.e., positive infinity.
 *  Calls python function py_eval_infinity for evalf(). */
const infinity Infinity = infinity::from_sign(+1);


/** NegInfinity, i.e., negative infinity.
 *  Calls python function py_eval_neg_infinity for evalf(). */
const infinity NegInfinity = infinity::from_sign(-1);


/** UnsignedInfinity. E.g., return value of gamma(-1).
 *  Calls python function py_eval_unsigned_infinity for evalf(). */
const infinity UnsignedInfinity = infinity::from_sign(0);


} // namespace GiNaC
