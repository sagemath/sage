/** @file inifcns_comb.cpp
 *
 *  Combinaorial functions. */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *  (C) 2016 Ralf Stephan <ralf@ark.in-berlin.de>
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

#include "inifcns.h"
#include "ex.h"
#include "constant.h"
#include "lst.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "pseries.h"
#include "symbol.h"
#include "utils.h"

#include <vector>
#include <stdexcept>
#include <sstream>

namespace GiNaC {

//////////
// factorial
//////////

static ex factorial_evalf(const ex & x, PyObject* parent)
{
	return factorial(x).hold();
}

static void factorial_print_dflt_latex(const ex & x, const print_context & c)
{
	if (is_exactly_a<symbol>(x) ||
	    is_exactly_a<constant>(x) ||
		is_exactly_a<function>(x)) {
		x.print(c); c.s << "!";
	} else {
		std::stringstream tstream;
		print_latex tcontext(tstream, c.options);
		x.print(tcontext);
		std::string argstr = tstream.str();
		bool parenthesis = ((argstr.find(' ') != std::string::npos)||
				(argstr.find('+') != std::string::npos) ||
				(argstr.find('-') != std::string::npos) ||
				(argstr.find('/') != std::string::npos) ||
				(argstr.find('*') != std::string::npos) ||
				(argstr.find('^') != std::string::npos));
		if (parenthesis)
			c.s << "\\left(";

		c.s << argstr;
		if (parenthesis)
			c.s << "\\right)";
		c.s<<"!";
	}
}

static ex factorial_conjugate(const ex & x)
{
	return factorial(x).hold();
}

static ex factorial_real_part(const ex & x)
{
	return factorial(x).hold();
}

static ex factorial_imag_part(const ex & x)
{
	return 0;
}

REGISTER_FUNCTION(factorial, evalf_func(factorial_evalf).
                             print_func<print_latex>(factorial_print_dflt_latex).
                             conjugate_func(factorial_conjugate).
                             real_part_func(factorial_real_part).
                             imag_part_func(factorial_imag_part));

//////////
// binomial
//////////

static ex binomial_evalf(const ex & x, const ex & y, PyObject* parent)
{
	return binomial(x, y).hold();
}

// At the moment the numeric evaluation of a binomail function always
// gives a real number, but if this would be implemented using the gamma
// function, also complex conjugation should be changed (or rather, deleted).
static ex binomial_conjugate(const ex & x, const ex & y)
{
	return binomial(x,y).hold();
}

static ex binomial_real_part(const ex & x, const ex & y)
{
	return binomial(x,y).hold();
}

static ex binomial_imag_part(const ex & x, const ex & y)
{
	return 0;
}

static void binomial_print_latex(const ex & x, const ex & y,
		const print_context & c)
{
	c.s<<"\\binom{";
	x.print(c);
	c.s<<"}{";
	y.print(c);
	c.s<<"}";
}


REGISTER_FUNCTION(binomial, evalf_func(binomial_evalf).
                            conjugate_func(binomial_conjugate).
                            real_part_func(binomial_real_part).
                            print_func<print_latex>(binomial_print_latex).
                            imag_part_func(binomial_imag_part));

//////////
// rising_factorial
//////////

static ex rising_factorial_evalf(const ex & x, const ex & y, PyObject* parent)
{
	return rising_factorial(x, y).hold();
}

static ex rising_factorial_eval(const ex & x, const ex &y)
{
	return rising_factorial(x, y).hold();
}

static void rising_factorial_print_latex(const ex & x, const ex & y,
		const print_context & c)
{
	c.s<<"{";
	x.print(c);
	c.s<<"}^{\\left({";
	y.print(c);
        c.s<<"}\\right)}";
}


REGISTER_FUNCTION(rising_factorial, eval_func(rising_factorial_eval).
                            evalf_func(rising_factorial_evalf).
                            print_func<print_latex>(rising_factorial_print_latex));

//////////
// falling_factorial
//////////

static ex falling_factorial_evalf(const ex & x, const ex & y, PyObject* parent)
{
	return falling_factorial(x, y).hold();
}

static ex falling_factorial_eval(const ex & x, const ex &y)
{
	return falling_factorial(x, y).hold();
}
static void falling_factorial_print_latex(const ex & x, const ex & y,
		const print_context & c)
{
	c.s<<"\\left({";
	x.print(c);
	c.s<<"}\\left)}_{";
	y.print(c);
	c.s<<"}";
}


REGISTER_FUNCTION(falling_factorial, eval_func(falling_factorial_eval).
                            evalf_func(falling_factorial_evalf).
                            print_func<print_latex>(falling_factorial_print_latex));

} // namespace GiNaC
