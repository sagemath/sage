/** @file operators.cpp
 *
 *  Implementation of GiNaC's overloaded operators. */

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

#include "numeric.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "relational.h"
#include "print.h"
#include "utils.h"

#include "operators.h"

#include <iostream>
#include <iomanip>

namespace GiNaC {

/** Used internally by operator+() to add two ex objects together. */
static inline const ex exadd(const ex & lh, const ex & rh)
{
        if (is_exactly_a<numeric>(lh)
            and is_exactly_a<numeric>(rh))
                return ex_to<numeric>(lh).add(ex_to<numeric>(rh));
	return (new add(lh,rh))->setflag(status_flags::dynallocated);
}

/** Used internally by operator*() to multiply two ex objects together. */
static inline const ex exmul(const ex & lh, const ex & rh)
{
        if (is_exactly_a<numeric>(lh)
            and is_exactly_a<numeric>(rh))
                return ex_to<numeric>(lh).mul(ex_to<numeric>(rh));
        return (new mul(lh,rh))->setflag(status_flags::dynallocated);
}

/** Used internally by operator-() and friends to change the sign of an argument. */
static inline const ex exminus(const ex & lh)
{
        if (is_exactly_a<numeric>(lh))
                return ex_to<numeric>(lh).negative();
	return (new mul(lh,_ex_1))->setflag(status_flags::dynallocated);
}

// binary arithmetic operators ex with ex

const ex operator+(const ex & lh, const ex & rh)
{
	return exadd(lh, rh);
}

const ex operator-(const ex & lh, const ex & rh)
{
	return exadd(lh, exminus(rh));
}

const ex operator*(const ex & lh, const ex & rh)
{
	return exmul(lh, rh);
}

const ex operator/(const ex & lh, const ex & rh)
{
	return exmul(lh, power(rh,_ex_1));
}


// binary arithmetic operators numeric with numeric

const numeric operator+(const numeric & lh, const numeric & rh)
{
	return lh.add(rh);
}

const numeric operator-(const numeric & lh, const numeric & rh)
{
	return lh.sub(rh);
}

const numeric operator*(const numeric & lh, const numeric & rh)
{
	return lh.mul(rh);
}

const numeric operator/(const numeric & lh, const numeric & rh)
{
	return lh.div(rh);
}


// binary arithmetic assignment operators with ex

ex & operator+=(ex & lh, const ex & rh)
{
	return lh = exadd(lh, rh);
}

ex & operator-=(ex & lh, const ex & rh)
{
	return lh = exadd(lh, exminus(rh));
}

ex & operator*=(ex & lh, const ex & rh)
{
	return lh = exmul(lh, rh);
}

ex & operator/=(ex & lh, const ex & rh)
{
	return lh = exmul(lh, power(rh,_ex_1));
}


// unary operators

const ex operator+(const ex & lh)
{
	return lh;
}

const ex operator-(const ex & lh)
{
	return exminus(lh);
}

const numeric operator+(const numeric & lh)
{
	return lh;
}

const numeric operator-(const numeric & lh)
{
	return lh.negative(); // better than _num_1_p->mul(lh)
}


// increment / decrement operators

/** Expression prefix increment.  Adds 1 and returns incremented ex. */
ex & operator++(ex & rh)
{
        if (is_exactly_a<numeric>(rh)) {
                rh = numeric(ex_to<numeric>(rh) + *_num1_p);
                return rh;
        }
	return rh = exadd(rh, _ex1);
}

/** Expression prefix decrement.  Subtracts 1 and returns decremented ex. */
ex & operator--(ex & rh)
{
        if (is_exactly_a<numeric>(rh)) {
                rh = numeric(ex_to<numeric>(rh) + *_num_1_p);
                return rh;
        }
	return rh = exadd(rh, _ex_1);
}

/** Expression postfix increment.  Returns the ex and leaves the original
 *  incremented by 1. */
const ex operator++(ex & lh, int)
{
	ex tmp(lh);
	lh = exadd(lh, _ex1);
	return tmp;
}

/** Expression postfix decrement.  Returns the ex and leaves the original
 *  decremented by 1. */
const ex operator--(ex & lh, int)
{
	ex tmp(lh);
	lh = exadd(lh, _ex_1);
	return tmp;
}

/** Numeric prefix increment.  Adds 1 and returns incremented number. */
numeric& operator++(numeric & rh)
{
	rh = rh.add(*_num1_p);
	return rh;
}

/** Numeric prefix decrement.  Subtracts 1 and returns decremented number. */
numeric& operator--(numeric & rh)
{
	rh = rh.add(*_num_1_p);
	return rh;
}

/** Numeric postfix increment.  Returns the number and leaves the original
 *  incremented by 1. */
const numeric operator++(numeric & lh, int)
{
	numeric tmp(lh);
	lh = lh.add(*_num1_p);
	return tmp;
}

/** Numeric postfix decrement.  Returns the number and leaves the original
 *  decremented by 1. */
const numeric operator--(numeric & lh, int)
{
	numeric tmp(lh);
	lh = lh.add(*_num_1_p);
	return tmp;
}

// binary relational operators ex with ex

const relational operator==(const ex & lh, const ex & rh)
{
	return relational(lh,rh,relational::equal);
}

const relational operator!=(const ex & lh, const ex & rh)
{
	return relational(lh,rh,relational::not_equal);
}

const relational operator<(const ex & lh, const ex & rh)
{
	return relational(lh,rh,relational::less);
}

const relational operator<=(const ex & lh, const ex & rh)
{
	return relational(lh,rh,relational::less_or_equal);
}

const relational operator>(const ex & lh, const ex & rh)
{
	return relational(lh,rh,relational::greater);
}

const relational operator>=(const ex & lh, const ex & rh)
{
	return relational(lh,rh,relational::greater_or_equal);
}

// input/output stream operators and manipulators

static int my_ios_index()
{
	static int i = std::ios_base::xalloc();
	return i;
}

// Stream format gets copied or destroyed
static void my_ios_callback(std::ios_base::event ev, std::ios_base & s, int i)
{
	std::unique_ptr<print_context> p(static_cast<print_context *>(s.pword(i)));
	if (ev == std::ios_base::erase_event) {
		s.pword(i) = nullptr;
	} else if (ev == std::ios_base::copyfmt_event && p != nullptr)
		s.pword(i) = p->duplicate();
}

enum {
	callback_registered = 1
};

// Get print_context associated with stream, may return 0 if no context has
// been associated yet
static inline print_context *get_print_context(std::ios_base & s)
{
	return static_cast<print_context *>(s.pword(my_ios_index()));
}

// Set print_context associated with stream, retain options
static void set_print_context(std::ios_base & s, const print_context & c)
{
	int i = my_ios_index();
	long flags = s.iword(i);
	if ((flags & callback_registered) == 0) {
		s.register_callback(my_ios_callback, i);
		s.iword(i) = flags | callback_registered;
	}
	print_context *p = static_cast<print_context *>(s.pword(i));
	unsigned options = p != nullptr ? p->options : c.options;
	delete p;
	p = c.duplicate();
	p->options = options;
	s.pword(i) = p;
}

// Get options for print_context associated with stream
static inline unsigned get_print_options(std::ios_base & s)
{
	print_context *p = get_print_context(s);
	return p != nullptr ? p->options : 0;
}

// Set options for print_context associated with stream
static void set_print_options(std::ostream & s, unsigned options)
{
	print_context *p = get_print_context(s);
	if (p == nullptr)
		set_print_context(s, print_dflt(s, options));
	else
		p->options = options;
}

std::ostream & operator<<(std::ostream & os, const ex & e)
{
	print_context *p = get_print_context(os);
	if (p == nullptr)
		e.print(print_dflt(os));
	else
		e.print(*p);
	return os;
}

std::ostream & operator<<(std::ostream & os, const exvector & e)
{
	print_context *p = get_print_context(os);
	auto i = e.begin();
	auto vend = e.end();

	if (i==vend) {
		os << "[]";
		return os;
	}

	os << "[";
	while (true) {
		if (p == nullptr)
			i -> print(print_dflt(os));
		else
			i -> print(*p);
		++i;
		if (i==vend)
			break;
		os << ",";
	}
	os << "]";

	return os;
}

std::ostream & operator<<(std::ostream & os, const exset & e)
{
	print_context *p = get_print_context(os);
	auto i = e.begin();
	auto send = e.end();

	if (i==send) {
		os << "<>";
		return os;
	}

	os << "<";
	while (true) {
		if (p == nullptr)
			i->print(print_dflt(os));
		else
			i->print(*p);
		++i;
		if (i == send)
			break;
		os << ",";
	}
	os << ">";

	return os;
}

std::ostream & operator<<(std::ostream & os, const exmap & e)
{
	print_context *p = get_print_context(os);
	auto i = e.begin();
	auto mend = e.end();

	if (i==mend) {
		os << "{}";
		return os;
	}

	os << "{";
	while (true) {
		if (p == nullptr)
			i->first.print(print_dflt(os));
		else
			i->first.print(*p);
		os << "==";
		if (p == nullptr)
			i->second.print(print_dflt(os));
		else
			i->second.print(*p);
		++i;
		if( i==mend )
			break;
		os << ",";
	}
	os << "}";

	return os;
}

std::istream & operator>>(std::istream & is, ex & e)
{
	throw (std::logic_error("expression input from streams not implemented"));
}

std::ostream & dflt(std::ostream & os)
{
	set_print_context(os, print_dflt(os));
	set_print_options(os, 0);
	return os;
}

std::ostream & latex(std::ostream & os)
{
	set_print_context(os, print_latex(os));
	return os;
}

std::ostream & python(std::ostream & os)
{
	set_print_context(os, print_python(os));
	return os;
}

std::ostream & python_repr(std::ostream & os)
{
	set_print_context(os, print_python_repr(os));
	return os;
}

std::ostream & tree(std::ostream & os)
{
	set_print_context(os, print_tree(os));
	return os;
}

std::ostream & index_dimensions(std::ostream & os)
{
	set_print_options(os, get_print_options(os) | print_options::print_index_dimensions);
	return os;
}

std::ostream & no_index_dimensions(std::ostream & os)
{
	set_print_options(os, get_print_options(os) & ~print_options::print_index_dimensions);
	return os;
}

} // namespace GiNaC
