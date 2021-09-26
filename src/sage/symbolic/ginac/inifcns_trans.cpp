/** @file inifcns_trans.cpp
 *
 *  Implementation of transcendental (and hyperbolic) functions. */

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

#include "inifcns.h"
#include "ex.h"
#include "ex_utils.h"
#include "lst.h"
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
#include "add.h"

#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <memory>

namespace GiNaC {


/* In Sage all the arc trig functions are printed with "arc" instead
   of "a" at the beginning.   This is for consistency with other
   computer algebra systems.   These print methods are registered
   below with each of the corresponding inverse trig function. */


//////////
// exponential function
//////////

static ex exp_eval(const ex & x)
{
	// exp(0) -> 1
	if (x.is_zero()) {
		return _ex1;
	}

	// handle infinity
	// This needs to be before the other tests below, since multiplying
	// infinity with other values throws runtime_errors
	// exp(oo) -> oo
	// exp(-oo) -> 0
	// exp(UnsignedInfinity) -> error
	if (is_exactly_a<infinity>(x)) {
	        const infinity& xinf = ex_to<infinity>(x);
		if (xinf.is_plus_infinity())
			return Infinity;
		if (xinf.is_minus_infinity())
			return _ex0;
		// x is UnsignedInfinity
		throw (std::runtime_error("exp_eval(): exp^(unsigned_infinity) encountered"));
	}

        bool has_pi = false, has_py = false;
        std::function<bool (const ex&)> has_pi_and_py =
        [&](const ex & the_ex) -> bool
        {
                if (is_exactly_a<constant>(the_ex)
                        and ex_to<constant>(the_ex) == Pi)
                        has_pi = true;
                if (is_a_python_object(the_ex))
                        has_py = true;
                if (has_pi and has_py)
                        return true;
                for (size_t i=0; i<the_ex.nops(); ++i)
                        if (has_pi_and_py(the_ex.op(i)))
                                return true;
                return false;
        };

	ex x_red;
        if (has_pi_and_py(x)) { // true if x contains Pi and Python objects
                                // like I. To make this check should be faster
                                // than the following

		ex coef_pi = x.coeff(Pi,_ex1).expand();
		ex rem = _ex0;
		if (is_exactly_a<add>(coef_pi)) {
			for (size_t i=0; i < coef_pi.nops(); i++) {
				if ((coef_pi.op(i) / (_ex2 * I)).is_integer())
					rem += Pi * coef_pi.op(i);
			}
		}
		else if ((coef_pi / (_ex2 * I)).is_integer())
			rem = Pi * coef_pi;
		x_red = (x - rem).expand();


		ex res1 = sin(x_red/I);
		ex res2 = cos(x_red/I);
		if (((not is_exactly_a<function>(res1))
		or (not is_ex_the_function(res1, sin)))
		and ((not is_exactly_a<function>(res2))
		or (not is_ex_the_function(res2, cos))))
			return I*res1 + res2;
        }
	else
		x_red = x;

        if (is_exactly_a<function>(x_red)) {
                const ex& arg = x_red.op(0);
                // exp(log(x)) -> x
                if (is_ex_the_function(x_red, log))
                        return arg;
                // exp(asinh(num)) etc.
                // see https://en.wikipedia.org/wiki/Hyperbolic_functions#Inverse_functions_as_logarithms
                if (is_exactly_a<numeric>(arg)) {
                        // exp(asinh(x)) -> x + sqrt(x^2 + 1)
                        if (is_ex_the_function(x_red, asinh))
                                return arg + sqrt(power(arg, _ex2) + _ex1);
                        // exp(acosh(x)) -> x + sqrt(x^2 - 1)
                        if (is_ex_the_function(x_red, acosh))
                                return arg + sqrt(power(arg, _ex2) - _ex1);
                        // exp(atanh(x)) -> sqrt((1 + x)/(1 - x))
                        if (is_ex_the_function(x_red, atanh))
                                return sqrt((_ex1 + arg) / (_ex1 - x));
                        // exp(acoth(x)) -> sqrt((x + 1)/(x - 1))
                        if (is_ex_the_function(x_red, acoth))
                                return sqrt((arg + _ex1) / (arg - _ex1));
                        // exp(asech(x)) -> 1/x + sqrt(1 - x^2)/x
                        if (is_ex_the_function(x_red, asech))
                                return (_ex1/arg +
                                        sqrt(_ex1 - power(arg, _ex2))/arg);
                        // exp(acsch(x)) -> 1/x + sqrt(1 + 1/x^2)
                        if (is_ex_the_function(x_red, acsch))
                                return (_ex1/arg +
                                        sqrt(_ex1 + _ex1/power(arg, _ex2)));
                }
        }

        static std::unordered_set<unsigned int> funcs = { log_SERIAL::serial,
                asinh_SERIAL::serial, acosh_SERIAL::serial,
                atanh_SERIAL::serial, acoth_SERIAL::serial,
                asech_SERIAL::serial, acsch_SERIAL::serial };

        if (is_exactly_a<mul>(x_red)) {
                // Loop through factors finding and marking the function.
                // Don't proceed if more than one function found.
                bool function_seen = false, applicable = false;
                size_t func_idx;
                for (size_t i=0; i<x_red.nops(); ++i) {
                        const ex& fac = x_red.op(i);
                        if (is_exactly_a<function>(fac)
                            and funcs.find(ex_to<function>(fac).get_serial()) != funcs.end()) {
                                if (function_seen)
                                        { applicable = false; break; }
                                function_seen = applicable = true;
                                func_idx = i;
                        }
                }
                if (applicable) {
                        ex c = _ex1;
                        for (size_t i=0; i<x_red.nops(); ++i)
                                if (i != func_idx)
                                        c *= x_red.op(i);
                        const ex& fac = x_red.op(func_idx);
                        const ex& arg = fac.op(0);
                        // exp(c*log(ex)) -> ex^c
                        if (is_ex_the_function(fac, log))
                                return power(arg, c);
                        // exp(c*asinh(ex)) etc.
                        // see https://en.wikipedia.org/wiki/Hyperbolic_functions#Inverse_functions_as_logarithms
                        // exp(asinh(x)) -> x + sqrt(x^2 + 1)
                        if (is_ex_the_function(fac, asinh))
                                return power(arg + sqrt(power(arg, _ex2) + 1),
                                             c);
                        // exp(acosh(x)) -> x + sqrt(x^2 - 1)
                        if (is_ex_the_function(fac, acosh))
                                return power(arg + sqrt(power(arg, _ex2) - 1),
                                             c);
                        // exp(atanh(x)) -> sqrt((1 + x)/(1 - x))
                        if (is_ex_the_function(fac, atanh))
                                return power((_ex1 + arg) / (_ex1 - arg), c/2);
                        // exp(acoth(x)) -> sqrt((x + 1)/(x - 1))
                        if (is_ex_the_function(fac, acoth))
                                return power((arg + _ex1) / (arg - _ex1), c/2);
                        // exp(asech(x)) -> 1/x + sqrt(1 - x^2)/x
                        if (is_ex_the_function(fac, asech))
                                return power(_ex1/arg +
                                             sqrt(_ex1 - power(arg, _ex2))/arg,
                                             c);
                        // exp(acsch(x)) -> 1/x + sqrt(1 + 1/x^2)
                        if (is_ex_the_function(fac, acsch))
                                return power(_ex1/arg +
                                             sqrt(_ex1 + _ex1/power(arg, _ex2)),
                                             c);
                }
        }

	// exp(float) -> float
	if (is_exactly_a<numeric>(x_red)
            and x_red.info(info_flags::inexact))
		return exp(ex_to<numeric>(x_red));

	return exp(x_red).hold();
}

static ex exp_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx exp(x) -> exp(x)
	return exp(x);
}

static ex exp_real_part(const ex & x)
{
	return exp(GiNaC::real_part(x))*cos(GiNaC::imag_part(x));
}

static ex exp_imag_part(const ex & x)
{
	return exp(GiNaC::real_part(x))*sin(GiNaC::imag_part(x));
}

static ex exp_power(const ex & arg, const ex & p)
{
	// If you change this function make sure to adjust
	// the simplification code in mul::eval accordingly
	if (is_exactly_a<numeric>(p) && ex_to<numeric>(p).is_integer())
		return exp(p*arg);

	return power(exp(arg), p).hold();
}

static void exp_print(const ex & arg, const print_context & c,
		bool latex=false)
{
	c.s << "e";
	if (!arg.is_equal(*_num1_p)) {
		c.s<<"^";
		std::stringstream tstream;
		std::unique_ptr<print_context> tcontext_p;
		if (latex) {
			tcontext_p.reset(new print_latex(tstream, c.options));
		} else {
			tcontext_p.reset(new print_context(tstream, c.options));
		}
		arg.print(*tcontext_p);
		std::string argstr = tstream.str();
		bool parenthesis = ((argstr.find(' ') != std::string::npos)||
				(argstr.find('+') != std::string::npos) ||
				(argstr.find('-') != std::string::npos) ||
				(argstr.find('/') != std::string::npos) ||
				(argstr.find('*') != std::string::npos) ||
				(argstr.find('^') != std::string::npos));
		if (latex) {
			c.s << '{';
			if (parenthesis)
				c.s << "\\left(";
		} else if (parenthesis)
			c.s << '(';

		c.s << argstr;
		if (latex) {
			if (parenthesis)
				c.s << "\\right)";
			c.s << '}';
		} else if (parenthesis)
			c.s << ')';
	}
}

static void exp_print_norm(const ex & arg, const print_context & c)
{
	exp_print(arg, c, false);
}

static void exp_print_latex(const ex & arg, const print_context & c)
{
	exp_print(arg, c, true);
}

static ex exp_conjugate(const ex & x)
{
	// conjugate(exp(x))==exp(conjugate(x))
	return exp(x.conjugate());
}

REGISTER_FUNCTION(exp, eval_func(exp_eval).
                       derivative_func(exp_deriv).
                       real_part_func(exp_real_part).
                       imag_part_func(exp_imag_part).
                       power_func(exp_power).
                       conjugate_func(exp_conjugate).
                       print_func<print_context>(exp_print_norm).
                       print_func<print_latex>(exp_print_latex));

//////////
// natural logarithm
//////////

static ex log_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {
		// log(float) -> float
                const numeric& n = ex_to<numeric>(x);
		if (x.info(info_flags::crational)) {
                        if (n.is_zero())         // log(0) -> infinity
                                //throw(pole_error("log_eval(): log(0)",0));
                                return NegInfinity;
                        if (not x.info(info_flags::inexact) and x.info(info_flags::negative))
                                return (log(-x)+I*Pi);
                        if (n.is_one())  // log(1) -> 0
                                return _ex0;
                        if (x.is_equal(I))       // log(I) -> Pi*I/2
                                return (Pi*I*_ex1_2);
                        if (x.is_equal(-I))      // log(-I) -> -Pi*I/2
                                return (Pi*I*_ex_1_2);
                        std::pair<int,int> p;
                        if (n.is_real() and n.is_integer()
                            and n.is_small_power(p))
                                return mul(p.second, log(p.first).hold());
                }
                else if (not x.info(info_flags::inexact))
                        return log(x).hold();
                else
			return log(n);
	}

	// log(exp(t)) -> t (if -Pi < t.imag() <= Pi):
	if (is_ex_the_function(x, exp)) {
		const ex &t = x.op(0);
		if (t.is_real())
			return t;
	}
	
	// log(oo) -> oo
	if (x.info(info_flags::infinity)) {
		return Infinity;
	}

	return log(x).hold();
}

static ex log_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx log(x) -> 1/x
	return power(x, _ex_1);
}

static ex log_series(const ex &arg,
                     const relational &rel,
                     int order,
                     unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	ex arg_pt;
	bool must_expand_arg = false;
	// maybe substitution of rel into arg fails because of a pole
	try {
		arg_pt = arg.subs(rel, subs_options::no_pattern);
	} catch (pole_error) {
		must_expand_arg = true;
	}
	// or we are at the branch point anyways
	if (arg_pt.is_zero())
		must_expand_arg = true;
	
	if (arg.diff(ex_to<symbol>(rel.lhs())).is_zero()) {
		throw do_taylor();
	}

	if (must_expand_arg) {
		// method:
		// This is the branch point: Series expand the argument first, then
		// trivially factorize it to isolate that part which has constant
		// leading coefficient in this fashion:
		//   x^n + x^(n+1) +...+ Order(x^(n+m))  ->  x^n * (1 + x +...+ Order(x^m)).
		// Return a plain n*log(x) for the x^n part and series expand the
		// other part.  Add them together and reexpand again in order to have
		// one unnested pseries object.  All this also works for negative n.
		pseries argser;          // series expansion of log's argument
		unsigned extra_ord = 0;  // extra expansion order
		do {
			// oops, the argument expanded to a pure Order(x^something)...
			argser = ex_to<pseries>(arg.series(rel, order+extra_ord, options));
			++extra_ord;
		} while (!argser.is_terminating() && argser.nops()==1);

		const symbol &s = ex_to<symbol>(rel.lhs());
		const ex &point = rel.rhs();
		const numeric &num = argser.ldegree(s);
                long n = num.to_long();
		epvector seq;
		// construct what we carelessly called the n*log(x) term above
		const ex coeff = argser.coeff(s, n);
		// expand the log, but only if coeff is real and > 0, since otherwise
		// it would make the branch cut run into the wrong direction
		if (coeff.is_positive())
			seq.emplace_back(n*log(s-point)+log(coeff), _ex0);
		else
			seq.emplace_back(log(coeff*pow(s-point, n)), _ex0);

		if (!argser.is_terminating() || argser.nops()!=1) {
			// in this case n more (or less) terms are needed
			// (sadly, to generate them, we have to start from the beginning)
			if (n == 0 and coeff.is_one()) {
				epvector epv;
				ex acc = (new pseries(rel, epv))->setflag(status_flags::dynallocated);
				epv.reserve(2);
				epv.emplace_back(-1, _ex0);
				epv.emplace_back(Order(_ex1), order);
				ex rest = pseries(rel, epv).add_series(argser);
				for (int i = order-1; i>0; --i) {
					epvector cterm;
					cterm.reserve(1);
					cterm.emplace_back((i%2) != 0 ? _ex1/i : _ex_1/i, _ex0);
					acc = pseries(rel, cterm).add_series(ex_to<pseries>(acc));
					acc = (ex_to<pseries>(rest)).mul_series(ex_to<pseries>(acc));
				}
				return acc;
			}
			const ex newarg = ex_to<pseries>((arg/coeff).series(rel, order+n, options)).shift_exponents(-n).convert_to_poly(true);
			return pseries(rel, seq).add_series(ex_to<pseries>(log(newarg).series(rel, order, options)));
		}  // it was a monomial
		return pseries(rel, seq);
	}
	if (((options & series_options::suppress_branchcut) == 0u) &&
	     arg_pt.info(info_flags::negative)) {
		// method:
		// This is the branch cut: assemble the primitive series manually and
		// then add the corresponding complex step function.
		const symbol &s = ex_to<symbol>(rel.lhs());
		const ex &point = rel.rhs();
		const symbol foo;
		const ex replarg = series(log(arg), s==foo, order).subs(foo==point, subs_options::no_pattern);
		epvector seq;
		seq.emplace_back(-I*csgn(arg*I)*Pi, _ex0);
		seq.emplace_back(Order(_ex1), order);
		return series(replarg - I*Pi + pseries(rel, seq), rel, order);
	}
	throw do_taylor();  // caught by function::series()
}

static ex log_real_part(const ex & x)
{
	if (x.is_positive())
		return log(x).hold();
	return log(abs(x));
}

static ex log_imag_part(const ex & x)
{
	if (x.is_positive())
		return _ex0;
	return atan2(GiNaC::imag_part(x), GiNaC::real_part(x));
}

static ex log_conjugate(const ex & x)
{
	// conjugate(log(x))==log(conjugate(x)) unless on the branch cut which
	// runs along the negative real axis.
	if (x.is_positive()) {
		return log(x);
	}
	if (is_exactly_a<numeric>(x) &&
	    !x.imag_part().is_zero()) {
		return log(x.conjugate());
	}
	return conjugate_function(log(x)).hold();
}

REGISTER_FUNCTION(log, eval_func(log_eval).
                       derivative_func(log_deriv).
                       series_func(log_series).
                       real_part_func(log_real_part).
                       imag_part_func(log_imag_part).
                       conjugate_func(log_conjugate).
                       latex_name("\\log"));

//////////
// General logarithm
// This only has shortcuts and delegates everything else to log(arg)
//////////

static ex logb_evalf(const ex & x, const ex & base, PyObject* parent)
{
        if ((base - exp(_ex1).hold()).is_zero()) {
                if (is_exactly_a<numeric>(x))
                        return log(ex_to<numeric>(x), parent);
                return log(x);
        }
	if (is_exactly_a<numeric>(x) and is_exactly_a<numeric>(base))
		return log(ex_to<numeric>(x), ex_to<numeric>(base));

	return mul(log(x), pow(log(base), _ex_1));
}

static ex logb_eval(const ex & x, const ex & base)
{
	if (is_exactly_a<numeric>(x) and not x.info(info_flags::inexact)
	    and is_exactly_a<numeric>(base) and not base.info(info_flags::inexact)) {
                const numeric& a = ex_to<numeric>(x);
                const numeric& b = ex_to<numeric>(base);
                if (b.is_real() and a.is_real()) {
                        bool israt;
                        numeric ret = a.ratlog(b, israt);
                        if (israt)
                                return ret;
                }
                return mul(log(x), pow(log(base), _ex_1));
        }

        if ((base - exp(_ex1).hold()).is_zero())
                return log_eval(x);

	// log(base^t, base) -> t
	if (is_exactly_a<power>(x)) {
                if (x.op(0).is_equal(base) and x.op(1).is_real())
			return x.op(1);
	}

	// log(oo) -> oo
	if (x.info(info_flags::infinity)) {
		return Infinity;
	}

        // log(x)/log(base)
	return mul(log(x), pow(log(base), _ex_1));
}

REGISTER_FUNCTION(logb, eval_func(logb_eval).
                       evalf_func(logb_evalf).
                       latex_name("\\log"));

//////////
// dilogarithm
//////////

static ex Li2_evalf(const ex & x, PyObject* parent)
{
	if (not is_exactly_a<numeric>(x))
	        return Li2(x).hold();

        return Li2(ex_to<numeric>(x), parent);
}

static ex Li2_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x) and not ex_to<numeric>(x).is_exact())
	        return Li2_evalf(x, nullptr);

	if (x.info(info_flags::numeric)) {
		// Li2(0) -> 0
		if (x.is_zero())
			return _ex0;
		// Li2(1) -> Pi^2/6
		if (x.is_one())
			return power(Pi,_ex2)/_ex6;
		// Li2(1/2) -> Pi^2/12 - log(2)^2/2
		if (x.is_equal(_ex1_2))
			return power(Pi,_ex2)/_ex12 + power(log(_ex2),_ex2)*_ex_1_2;
		// Li2(-1) -> -Pi^2/12
		if (x.is_minus_one())
			return -power(Pi,_ex2)/_ex12;
		// Li2(I) -> -Pi^2/48+Catalan*I
		if (x.is_equal(I))
			return power(Pi,_ex2)/_ex_48 + Catalan*I;
		// Li2(-I) -> -Pi^2/48-Catalan*I
		if (x.is_equal(-I))
			return power(Pi,_ex2)/_ex_48 - Catalan*I;
		// Li2(float)
		if (x.info(info_flags::inexact))
			return Li2(ex_to<numeric>(x));
	}
	
	return Li2(x).hold();
}

static ex Li2_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx Li2(x) -> -log(1-x)/x
	return -log(_ex1-x)/x;
}

static ex Li2_series(const ex &x, const relational &rel, int order, unsigned options)
{
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (x_pt.info(info_flags::numeric)) {
		// First special case: x==0 (derivatives have poles)
		if (x_pt.is_zero()) {
			// method:
			// The problem is that in d/dx Li2(x==0) == -log(1-x)/x we cannot 
			// simply substitute x==0.  The limit, however, exists: it is 1.
			// We also know all higher derivatives' limits:
			// (d/dx)^n Li2(x) == n!/n^2.
			// So the primitive series expansion is
			// Li2(x==0) == x + x^2/4 + x^3/9 + ...
			// and so on.
			// We first construct such a primitive series expansion manually in
			// a dummy symbol s and then insert the argument's series expansion
			// for s.  Reexpanding the resulting series returns the desired
			// result.
			const symbol s;
			ex ser;
			// manually construct the primitive expansion
			for (int i=1; i<order; ++i)
				ser += pow(s,i) / pow(numeric(i), *_num2_p);
			// substitute the argument's series expansion
			ser = ser.subs(s==x.series(rel, order), subs_options::no_pattern);
			// maybe that was terminating, so add a proper order term
			epvector nseq;
			nseq.emplace_back(Order(_ex1), order);
			ser += pseries(rel, nseq);
			// reexpanding it will collapse the series again
			return ser.series(rel, order);
			// NB: Of course, this still does not allow us to compute anything
			// like sin(Li2(x)).series(x==0,2), since then this code here is
			// not reached and the derivative of sin(Li2(x)) doesn't allow the
			// substitution x==0.  Probably limits *are* needed for the general
			// cases.  In case L'Hospital's rule is implemented for limits and
			// basic::series() takes care of this, this whole block is probably
			// obsolete!
		}
		// second special case: x==1 (branch point)
		if (x_pt.is_equal(_ex1)) {
			// method:
			// construct series manually in a dummy symbol s
			const symbol s;
			ex ser = zeta(_ex2);
			// manually construct the primitive expansion
			for (int i=1; i<order; ++i)
				ser += pow(1-s,i) * (numeric(1,i)*(I*Pi+log(s-1)) - numeric(1,i*i));
			// substitute the argument's series expansion
			ser = ser.subs(s==x.series(rel, order), subs_options::no_pattern);
			// maybe that was terminating, so add a proper order term
			epvector nseq;
			nseq.emplace_back(Order(_ex1), order);
			ser += pseries(rel, nseq);
			// reexpanding it will collapse the series again
			return ser.series(rel, order);
		}
		// third special case: x real, >=1 (branch cut)
		if (((options & series_options::suppress_branchcut) == 0u) &&
			ex_to<numeric>(x_pt).is_real() && ex_to<numeric>(x_pt)>1) {
			// method:
			// This is the branch cut: assemble the primitive series manually
			// and then add the corresponding complex step function.
			const symbol &s = ex_to<symbol>(rel.lhs());
			const ex point = rel.rhs();
			const symbol foo;
			epvector seq;
			// zeroth order term:
			seq.emplace_back(Li2(x_pt), _ex0);
			// compute the intermediate terms:
			ex replarg = series(Li2(x), s==foo, order);
			for (unsigned i=1; i < replarg.nops()-1; ++i) {
				ex term = replarg.op(i) / power(s-foo, i);
                                term = term.series(foo==point,1,options).op(0);
                                term.subs(foo==s, subs_options::no_pattern);
				seq.emplace_back(term, numeric(i));
                        }
			// append an order term:
			seq.emplace_back(Order(_ex1),
                                                long(replarg.nops()-1));
			return pseries(rel, seq);
		}
	}
	// all other cases should be safe, by now:
	throw do_taylor();  // caught by function::series()
}

static ex Li2_conjugate(const ex & x)
{
	// conjugate(Li2(x))==Li2(conjugate(x)) unless on the branch cuts which
	// run along the positive real axis beginning at 1.
	if (x.info(info_flags::negative)) {
		return Li2(x).hold();
	}
	if (is_exactly_a<numeric>(x) &&
	    (!x.imag_part().is_zero() || x < *_num1_p)) {
		return Li2(x.conjugate());
	}
	return conjugate_function(Li2(x)).hold();
}


unsigned Li2_SERIAL::serial = function::register_new(function_options("dilog", 1).
                       eval_func(Li2_eval).
                       evalf_func(Li2_evalf).
                       derivative_func(Li2_deriv).
                       series_func(Li2_series).
                       conjugate_func(Li2_conjugate).
                       latex_name("{\\rm Li}_2"));

//////////////////////////////////////
// classical polylog function (see inifcns_nstdsum.cpp for multiple)
//////////////////////////////////////

static ex Li_evalf(const ex& m_, const ex& x_, PyObject* parent)
{
	if (not is_exactly_a<numeric>(m_)
            or not is_exactly_a<numeric>(x_))
                return Li(m_,x_).hold();

        const numeric& num_m = ex_to<numeric>(m_);
        const numeric& num_x = ex_to<numeric>(x_);
        
        return Li2(num_m, num_x, parent);
}


static ex Li_eval(const ex& m_, const ex& x_)
{
        if ((is_exactly_a<numeric>(x_) and not ex_to<numeric>(x_).is_exact())
            or (is_exactly_a<numeric>(m_) and not ex_to<numeric>(m_).is_exact())) {
	        return Li_evalf(m_, x_, nullptr);
	}

	// classical polylogs
        if ((is_exactly_a<numeric>(x_) and not ex_to<numeric>(x_).is_exact())
            or (is_exactly_a<numeric>(m_) and not ex_to<numeric>(m_).is_exact())) {
	        return Li_evalf(m_, x_, nullptr);
	}

	if (x_.is_zero()) {
		return _ex0;
	}
	if (x_.is_one()) {
		return zeta(m_);
	}
	if (x_.is_minus_one()) {
		return (pow(2,1-m_)-1) * zeta(m_);
	}
	if (m_.is_one()) {
		return -log(1-x_);
	}
	if (m_.is_equal(_ex2)) 
                return Li2_eval(x_);

	return Li(m_, x_).hold();
}


static ex Li_series(const ex& m, const ex& x, const relational& rel, int order, unsigned options)
{
	// classical polylog
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (m.info(info_flags::numeric) && x_pt.info(info_flags::numeric)) {
		// First special case: x==0 (derivatives have poles)
		if (x_pt.is_zero()) {
			const symbol s;
			ex ser;
			// manually construct the primitive expansion
			for (int i=1; i<order; ++i)
				ser += pow(s,i) / pow(numeric(i), m);
			// substitute the argument's series expansion
			ser = ser.subs(s==x.series(rel, order), subs_options::no_pattern);
			// maybe that was terminating, so add a proper order term
			epvector nseq;
			nseq.emplace_back(Order(_ex1), order);
			ser += pseries(rel, nseq);
			// reexpanding it will collapse the series again
			return ser.series(rel, order);
		}
		// TODO special cases: x==1 (branch point) and x real, >=1 (branch cut)
		throw std::runtime_error("Li_series: don't know how to do the series expansion at this point!");
	}
	// all other cases should be safe, by now:
	throw do_taylor();  // caught by function::series()
}


static ex Li_deriv(const ex& m_, const ex& x_, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param < 2);
	if (deriv_param == 0) {
		return _ex0;
	}
	const ex& m = m_;
	const ex& x = x_;
	if (m > 0) {
		return Li(m-1, x) / x;
	} 
	return 1/(1-x);
}


static void Li_print_latex(const ex& m_, const ex& x_, const print_context& c)
{
	lst m;
	if (is_a<lst>(m_)) {
		m = ex_to<lst>(m_);
	} else {
		m = lst(m_);
	}
	lst x;
	if (is_a<lst>(x_)) {
		x = ex_to<lst>(x_);
	} else {
		x = lst(x_);
	}
	c.s << "{\\rm Li}_{";
	auto itm = m.begin();
	(*itm).print(c);
	itm++;
	for (; itm != m.end(); itm++) {
		c.s << ",";
		(*itm).print(c);
	}
	c.s << "}(";
	auto itx = x.begin();
	(*itx).print(c);
	itx++;
	for (; itx != x.end(); itx++) {
		c.s << ",";
		(*itx).print(c);
	}
	c.s << ")";
}


unsigned Li_SERIAL::serial = function::register_new(function_options("polylog", 2).
                  evalf_func(Li_evalf).
                  eval_func(Li_eval).
                  series_func(Li_series).
                  derivative_func(Li_deriv).
                  print_func<print_latex>(Li_print_latex).
                  do_not_evalf_params());

} // namespace GiNaC
