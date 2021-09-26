/** @file inifcns.cpp
 *
 *  Implementation of GiNaC's initially known functions. */

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
#include "constant.h"
#include "lst.h"
#include "fderivative.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "pseries.h"
#include "symbol.h"
#include "utils.h"
#include "wildcard.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <sstream>

namespace GiNaC {

//////////
// complex conjugate
//////////

static ex conjugate_evalf(const ex & arg, PyObject* parent)
{
	if (is_exactly_a<numeric>(arg)) {
		return ex_to<numeric>(arg).conjugate();
	}
	return conjugate_function(arg).hold();
}

static ex conjugate_eval(const ex & arg)
{
	return arg.conjugate();
}

static void conjugate_print_latex(const ex & arg, const print_context & c)
{
	c.s << "\\overline{"; arg.print(c); c.s << "}";
}

static ex conjugate_conjugate(const ex & arg)
{
	return arg;
}

// If x is real then U.diff(x)-I*V.diff(x) represents both conjugate(U+I*V).diff(x) 
// and conjugate((U+I*V).diff(x))
static ex conjugate_expl_derivative(const ex & arg, const symbol & s)
{
	if (s.is_real())
		return conjugate(arg.diff(s));

        exvector vec_arg;
        vec_arg.push_back(arg);
        return fderivative(ex_to<function>(conjugate(arg)).get_serial(),
                        0, vec_arg).hold() * arg.diff(s);
}

static ex conjugate_real_part(const ex & arg)
{
	return arg.real_part();
}

static ex conjugate_imag_part(const ex & arg)
{
	return -arg.imag_part();
}

REGISTER_FUNCTION(conjugate_function, eval_func(conjugate_eval).
                                      evalf_func(conjugate_evalf).
                                      expl_derivative_func(conjugate_expl_derivative).
                                      print_func<print_latex>(conjugate_print_latex).
                                      conjugate_func(conjugate_conjugate).
                                      real_part_func(conjugate_real_part).
                                      imag_part_func(conjugate_imag_part).
                                      set_name("conjugate","conjugate"));

//////////
// real part
//////////

static ex real_part_evalf(const ex & arg, PyObject* parent)
{
	if (is_exactly_a<numeric>(arg)) {
		return ex_to<numeric>(arg).real();
	}
	return real_part_function(arg).hold();
}

static ex real_part_eval(const ex & arg)
{
        if (arg.is_real())
                return arg;
        ex pat = pow(wild(1), wild(2));
        lst l;
        if (arg.find(pat, l))
                // real parts of powers can be expensive
                return arg.expand().real_part();

        return arg.real_part();
}

static void real_part_print_latex(const ex & arg, const print_context & c)
{
	c.s << "\\Re \\left( "; arg.print(c); c.s << " \\right)";
}

static ex real_part_conjugate(const ex & arg)
{
	return real_part_function(arg).hold();
}

static ex real_part_real_part(const ex & arg)
{
	return real_part_function(arg).hold();
}

static ex real_part_imag_part(const ex & arg)
{
	return 0;
}

// If x is real then Re(e).diff(x) is equal to Re(e.diff(x)) 
static ex real_part_expl_derivative(const ex & arg, const symbol & s)
{
	if (s.is_real())
		return real_part_function(arg.diff(s));
	
        exvector vec_arg;
        vec_arg.push_back(arg);
        return fderivative(ex_to<function>(real_part(arg)).get_serial(),
                        0, vec_arg).hold() * arg.diff(s);
}

REGISTER_FUNCTION(real_part_function, eval_func(real_part_eval).
                                      evalf_func(real_part_evalf).
                                      expl_derivative_func(real_part_expl_derivative).
                                      print_func<print_latex>(real_part_print_latex).
                                      conjugate_func(real_part_conjugate).
                                      real_part_func(real_part_real_part).
                                      imag_part_func(real_part_imag_part).
                                      set_name("real_part","real_part"));

//////////
// imag part
//////////

static ex imag_part_evalf(const ex & arg, PyObject* parent)
{
	if (is_exactly_a<numeric>(arg)) {
		return ex_to<numeric>(arg).imag();
	}
	return imag_part_function(arg).hold();
}

static ex imag_part_eval(const ex & arg)
{
        if (arg.is_real())
                return _ex0;
        ex pat = pow(wild(1), wild(2));
        lst l;
        if (arg.find(pat, l))
                // imag parts of powers can be expensive
                return arg.expand().imag_part();

        return arg.imag_part();
}

static void imag_part_print_latex(const ex & arg, const print_context & c)
{
	c.s << "\\Im \\left( "; arg.print(c); c.s << " \\right)";
}

static ex imag_part_conjugate(const ex & arg)
{
	return imag_part_function(arg).hold();
}

static ex imag_part_real_part(const ex & arg)
{
	return imag_part_function(arg).hold();
}

static ex imag_part_imag_part(const ex & arg)
{
	return 0;
}

// If x is real then Im(e).diff(x) is equal to Im(e.diff(x)) 
static ex imag_part_expl_derivative(const ex & arg, const symbol & s)
{
	if (s.is_real())
		return imag_part_function(arg.diff(s));

        exvector vec_arg;
        vec_arg.push_back(arg);
        return fderivative(ex_to<function>(imag_part(arg)).get_serial(),
                        0, vec_arg).hold() * arg.diff(s);
}

REGISTER_FUNCTION(imag_part_function, eval_func(imag_part_eval).
                                      evalf_func(imag_part_evalf).
                                      expl_derivative_func(imag_part_expl_derivative).
                                      print_func<print_latex>(imag_part_print_latex).
                                      conjugate_func(imag_part_conjugate).
                                      real_part_func(imag_part_real_part).
                                      imag_part_func(imag_part_imag_part).
                                      set_name("imag_part","imag_part"));

//////////
// absolute value
//////////

static ex abs_evalf(const ex & arg, PyObject* parent)
{
	if (is_exactly_a<numeric>(arg))
		return abs(ex_to<numeric>(arg));
	
	return abs(arg).hold();
}

static ex abs_eval(const ex & arg)
{
	if (is_exactly_a<numeric>(arg))
		return abs(ex_to<numeric>(arg));

	if (arg.info(info_flags::nonnegative)
            or arg.is_positive())
		return arg;

	if (arg.info(info_flags::negative) or (-arg).info(info_flags::nonnegative))
		return -arg;

        if (is_exactly_a<function>(arg)) {                
                if (is_ex_the_function(arg, abs))
                        return arg;
                if (is_ex_the_function(arg, exp))
                        return exp(arg.op(0).real_part());
                if (is_ex_the_function(arg, conjugate_function))
                        return abs(arg.op(0));
                if (is_ex_the_function(arg, unit_step)
                    or is_ex_the_function(arg, heaviside))
                        return arg;
                return abs(arg).hold();
        }

        if (is_exactly_a<mul>(arg)) {
                ex prod = _ex1;
                ex prod_sy = _ex1;
                bool is_prod_neg = false;
                for (size_t i=0; i<arg.nops(); ++i) {
                        const ex& factor = arg.op(i);
                        if (has_symbol(factor))
                                prod_sy *= factor;
                        else if (factor.is_real()) {
                                if (factor.info(info_flags::negative)) {                                        
                                        is_prod_neg = not is_prod_neg;
                                        prod *= factor;
                                }
                                else if (factor.is_positive())                                        
                                        prod *= factor;
                                else
                                        prod *= abs(factor).hold();
                        }
                        else
                                prod *= abs(factor);
                }
                if (is_prod_neg)
                        prod *= _ex_1;
                if (not prod_sy.is_one())
                        prod *= abs(prod_sy).hold();
                return prod;
	}

	if (is_exactly_a<power>(arg)) {
		const ex& base = arg.op(0);
		const ex& exponent = arg.op(1);
		if (base.is_positive()
                    or exponent.is_real())
			return pow(abs(base), exponent.real_part());
	}

	return abs(arg).hold();
}

static ex abs_expl_derivative(const ex & arg, const symbol & s)
{
	ex diff_arg = arg.diff(s);
	return (diff_arg*arg.conjugate()+arg*diff_arg.conjugate())/2/abs(arg);
}

static void abs_print_latex(const ex & arg, const print_context & c)
{
	c.s << "{\\left| "; arg.print(c); c.s << " \\right|}";
}

static ex abs_conjugate(const ex & arg)
{
	return abs(arg).hold();
}

static ex abs_real_part(const ex & arg)
{
	return abs(arg).hold();
}

static ex abs_imag_part(const ex& arg)
{
	return 0;
}

static ex abs_power(const ex & arg, const ex & exp)
{
	if (((is_exactly_a<numeric>(exp) and ex_to<numeric>(exp).is_even()) 
                or exp.info(info_flags::even))
                and (arg.is_real() or arg.is_equal(arg.conjugate())))
	        return power(arg, exp);

        return power(abs(arg), exp).hold();
}

static ex abs_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx abs(x) -> x/abs(x)
	return x/abs(x);
}


REGISTER_FUNCTION(abs, eval_func(abs_eval).
                       evalf_func(abs_evalf).
                       expl_derivative_func(abs_expl_derivative).
                       print_func<print_latex>(abs_print_latex).
                       derivative_func(abs_deriv).
                       conjugate_func(abs_conjugate).
                       real_part_func(abs_real_part).
                       imag_part_func(abs_imag_part).
                       power_func(abs_power));

//////////
// Step function
//////////

static ex unit_step_evalf(const ex & arg, PyObject* parent)
{
	if (not is_exactly_a<numeric>(arg))
	        return unit_step(arg).hold();

        if (arg.is_zero())
                return _ex1;

        if (arg.is_real())
		return step(ex_to<numeric>(arg));
	
	return unit_step(arg).hold();
}

static ex unit_step_eval(const ex & arg)
{
        if (arg.is_positive())
                return _ex1;
        if (arg.info(info_flags::negative))
                return _ex0;
	if (is_exactly_a<numeric>(arg))
		return unit_step_evalf(arg, nullptr);
	if (is_exactly_a<mul>(arg)) {
		const numeric& oc = ex_to<mul>(arg).get_overall_coeff();
		if (oc.is_real()) {
			if (oc > 0)
				// step(42*x) -> step(x)
				return unit_step(arg/oc).hold();

                        // step(-42*x) -> step(-x)
                        return unit_step(-arg/oc).hold();
		}
		if (oc.real().is_zero()) {
			if (oc.imag() > 0)
				// step(42*I*x) -> step(I*x)
				return unit_step(I*arg/oc).hold();

                        // step(-42*I*x) -> step(-I*x)
                        return unit_step(-I*arg/oc).hold();
		}
	}
	
	return unit_step(arg).hold();
}

static ex unit_step_series(const ex & arg,
                           const relational & rel,
                           int order,
                           unsigned options)
{
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (is_exactly_a<numeric>(arg_pt)
	    && ex_to<numeric>(arg_pt).real().is_zero()
	    && ((options & series_options::suppress_branchcut) == 0u))
		throw (std::domain_error("unit_step_series(): on imaginary axis"));
	
	epvector seq;
	seq.emplace_back(unit_step(arg_pt), _ex0);
	return pseries(rel,seq);
}

static ex unit_step_conjugate(const ex& arg)
{
	return unit_step(arg).hold();
}

static ex unit_step_real_part(const ex& arg)
{
	return unit_step(arg).hold();
}

static ex unit_step_imag_part(const ex& arg)
{
	return 0;
}

REGISTER_FUNCTION(unit_step, eval_func(unit_step_eval).
                        evalf_func(unit_step_evalf).
                        series_func(unit_step_series).
                        conjugate_func(unit_step_conjugate).
                        real_part_func(unit_step_real_part).
                        imag_part_func(unit_step_imag_part));

//////////
// Heaviside function
//////////

static ex heaviside_evalf(const ex & arg, PyObject* parent)
{
	if (is_exactly_a<numeric>(arg)
            and arg.is_real()
            and not arg.is_zero())
		return step(ex_to<numeric>(arg));
	
	return heaviside(arg).hold();
}

static ex heaviside_eval(const ex & arg)
{
        if (arg.is_positive())
                return _ex1;
        if (arg.info(info_flags::negative))
                return _ex0;
	if (is_exactly_a<numeric>(arg))
		return heaviside_evalf(arg, nullptr);
	
	if (is_exactly_a<mul>(arg)) {
		const numeric& oc = ex_to<mul>(arg).get_overall_coeff();
		if (oc.is_real()) {
			if (oc > 0)
				// step(42*x) -> step(x)
				return heaviside(arg/oc).hold();

                        // step(-42*x) -> step(-x)
                        return heaviside(-arg/oc).hold();
		}
		if (oc.real().is_zero()) {
			if (oc.imag() > 0)
				// step(42*I*x) -> step(I*x)
				return heaviside(I*arg/oc).hold();

                        // step(-42*I*x) -> step(-I*x)
                        return heaviside(-I*arg/oc).hold();
		}
	}
	
	return heaviside(arg).hold();
}

static ex heaviside_series(const ex & arg,
                           const relational & rel,
                           int order,
                           unsigned options)
{
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (is_exactly_a<numeric>(arg_pt)
	    && ex_to<numeric>(arg_pt).real().is_zero()
	    && ((options & series_options::suppress_branchcut) == 0u))
		throw (std::domain_error("heaviside_series(): on imaginary axis"));
	
	epvector seq;
	seq.emplace_back(heaviside(arg_pt), _ex0);
	return pseries(rel,seq);
}

static ex heaviside_conjugate(const ex& arg)
{
	return heaviside(arg).hold();
}

static ex heaviside_real_part(const ex& arg)
{
	return heaviside(arg).hold();
}

static ex heaviside_imag_part(const ex& arg)
{
	return 0;
}

REGISTER_FUNCTION(heaviside, eval_func(heaviside_eval).
                        evalf_func(heaviside_evalf).
                        series_func(heaviside_series).
                        conjugate_func(heaviside_conjugate).
                        real_part_func(heaviside_real_part).
                        imag_part_func(heaviside_imag_part));

//////////
// Complex sign
//////////

static ex csgn_evalf(const ex & arg, PyObject* parent)
{
	if (is_exactly_a<numeric>(arg))
		return csgn(ex_to<numeric>(arg));
	
	return csgn(arg).hold();
}

static ex csgn_eval(const ex & arg)
{
	if (is_exactly_a<numeric>(arg))
		return csgn(ex_to<numeric>(arg));
	
	if (is_exactly_a<mul>(arg)) {
		const numeric& oc = ex_to<mul>(arg).get_overall_coeff();
		if (oc.is_real()) {
			if (oc > 0)
				// csgn(42*x) -> csgn(x)
				return csgn(arg/oc).hold();

                        // csgn(-42*x) -> -csgn(x)
                        return -csgn(arg/oc).hold();
		}
		if (oc.real().is_zero()) {
			if (oc.imag() > 0)
				// csgn(42*I*x) -> csgn(I*x)
				return csgn(I*arg/oc).hold();

                        // csgn(-42*I*x) -> -csgn(I*x)
                        return -csgn(I*arg/oc).hold();
		}
	}
	
	return csgn(arg).hold();
}

static ex csgn_series(const ex & arg,
                      const relational & rel,
                      int order,
                      unsigned options)
{
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (arg_pt.info(info_flags::numeric)
	    && ex_to<numeric>(arg_pt).real().is_zero()
	    && ((options & series_options::suppress_branchcut) == 0u))
		throw (std::domain_error("csgn_series(): on imaginary axis"));
	
	epvector seq;
	seq.emplace_back(csgn(arg_pt), _ex0);
	return pseries(rel,seq);
}

static ex csgn_conjugate(const ex& arg)
{
	return csgn(arg).hold();
}

static ex csgn_real_part(const ex& arg)
{
	return csgn(arg).hold();
}

static ex csgn_imag_part(const ex& arg)
{
	return 0;
}

static ex csgn_power(const ex & arg, const ex & exp)
{
	if (is_a<numeric>(exp)
            and exp.is_positive()
            and exp.is_integer()) {
		if (ex_to<numeric>(exp).is_odd())
			return csgn(arg).hold();

                return power(csgn(arg), _ex2).hold();
        }
        return power(csgn(arg), exp).hold();
}


REGISTER_FUNCTION(csgn, eval_func(csgn_eval).
                        evalf_func(csgn_evalf).
                        series_func(csgn_series).
                        conjugate_func(csgn_conjugate).
                        real_part_func(csgn_real_part).
                        imag_part_func(csgn_imag_part).
                        power_func(csgn_power));


//////////
// Eta function: eta(x,y) == log(x*y) - log(x) - log(y).
// This function is closely related to the unwinding number K, sometimes found
// in modern literature: K(z) == (z-log(exp(z)))/(2*Pi*I).
//////////

static ex eta_evalf(const ex &x, const ex &y, PyObject* parent)
{
	// It seems like we basically have to replicate the eval function here,
	// since the expression might not be fully evaluated yet.
	if (x.is_positive() or y.is_positive())
		return _ex0;

	if (x.info(info_flags::numeric) &&	y.info(info_flags::numeric)) {
		const numeric& nx = ex_to<numeric>(x);
		const numeric& ny = ex_to<numeric>(y);
		const numeric nxy = ex_to<numeric>(x*y);
		int cut = 0;
		if (nx.is_real() && nx.is_negative())
			cut -= 4;
		if (ny.is_real() && ny.is_negative())
			cut -= 4;
		if (nxy.is_real() && nxy.is_negative())
			cut += 4;
		return evalf(I/4*Pi, 0, parent)*((csgn(-imag(nx))+1)*(csgn(-imag(ny))+1)*(csgn(imag(nxy))+1)-
		                      (csgn(imag(nx))+1)*(csgn(imag(ny))+1)*(csgn(-imag(nxy))+1)+cut);
	}

	return eta(x,y).hold();
}

static ex eta_eval(const ex &x, const ex &y)
{
	// trivial:  eta(x,c) -> 0  if c is real and positive
	if (x.is_positive() or y.is_positive())
		return _ex0;

	if (x.info(info_flags::numeric) &&	y.info(info_flags::numeric)) {
		// don't call eta_evalf here because it would call Pi.evalf()!
		const numeric& nx = ex_to<numeric>(x);
		const numeric& ny = ex_to<numeric>(y);
		const numeric nxy = ex_to<numeric>(x*y);
		int cut = 0;
		if (nx.is_real() && nx.is_negative())
			cut -= 4;
		if (ny.is_real() && ny.is_negative())
			cut -= 4;
		if (nxy.is_real() && nxy.is_negative())
			cut += 4;
		return (I/4)*Pi*((csgn(-imag(nx))+1)*(csgn(-imag(ny))+1)*(csgn(imag(nxy))+1)-
		                 (csgn(imag(nx))+1)*(csgn(imag(ny))+1)*(csgn(-imag(nxy))+1)+cut);
	}
	
	return eta(x,y).hold();
}

static ex eta_series(const ex & x, const ex & y,
                     const relational & rel,
                     int order,
                     unsigned options)
{
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	const ex y_pt = y.subs(rel, subs_options::no_pattern);
	if ((x_pt.info(info_flags::numeric) && x_pt.info(info_flags::negative)) ||
	    (y_pt.info(info_flags::numeric) && y_pt.info(info_flags::negative)) ||
	    ((x_pt*y_pt).info(info_flags::numeric) && (x_pt*y_pt).info(info_flags::negative)))
			throw (std::domain_error("eta_series(): on discontinuity"));
	epvector seq;
	seq.emplace_back(eta(x_pt,y_pt), _ex0);
	return pseries(rel,seq);
}

static ex eta_conjugate(const ex & x, const ex & y)
{
	return -eta(x, y).hold();
}

static ex eta_real_part(const ex & x, const ex & y)
{
	return 0;
}

static ex eta_imag_part(const ex & x, const ex & y)
{
	return -I*eta(x, y).hold();
}

REGISTER_FUNCTION(eta, eval_func(eta_eval).
                       evalf_func(eta_evalf).
                       series_func(eta_series).
                       latex_name("\\eta").
                       conjugate_func(eta_conjugate).
                       real_part_func(eta_real_part).
                       imag_part_func(eta_imag_part));


//////////
// Order term function (for truncated power series)
//////////

static ex Order_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {
		// O(c) -> O(1) or 0
		if (!x.is_zero())
			return Order(_ex1).hold();
                return _ex0;
	}
        if (is_exactly_a<mul>(x)) {
		const mul &m = ex_to<mul>(x);
		// O(c*expr) -> O(expr)
		return Order(x / m.get_overall_coeff()).hold();
	}
	return Order(x).hold();
}

static ex Order_series(const ex & x, const relational & r, int order, unsigned options)
{
	// Just wrap the function into a pseries object
	epvector new_seq;
	GINAC_ASSERT(is_a<symbol>(r.lhs()));
	const symbol &s = ex_to<symbol>(r.lhs());
        int ldeg = x.ldegree(s).to_int();
	new_seq.emplace_back(Order(_ex1), numeric(std::min(ldeg, order)));
	return pseries(r, new_seq);
}

static ex Order_conjugate(const ex & x)
{
	return Order(x).hold();
}

static ex Order_real_part(const ex & x)
{
	return Order(x).hold();
}

static ex Order_imag_part(const ex & x)
{
	if(x.is_real())
		return 0;
	return Order(x).hold();
}

static ex Order_expl_derivative(const ex & arg, const symbol & s)
{
	return Order(arg.diff(s));
}

REGISTER_FUNCTION(Order, eval_func(Order_eval).
                         series_func(Order_series).
                         latex_name("\\mathcal{O}").
                         expl_derivative_func(Order_expl_derivative).
                         conjugate_func(Order_conjugate).
                         real_part_func(Order_real_part).
                         imag_part_func(Order_imag_part));

static ex cases_eval(const ex & arg)
{
        if (not is_exactly_a<exprseq>(arg))
                throw std::runtime_error("cases argument not a sequence");
        bool default_seen = false;
        ex deflt;
        for (auto case_ : arg) {
                if (not is_exactly_a<exprseq>(case_)) {
                        if (not default_seen) {
                                deflt = case_;
                                default_seen = true;
                        }
                        continue;
                }
                const exprseq& seq = ex_to<exprseq>(case_);
                if (seq.nops() == 1) {
                        if (not default_seen) {
                                deflt = seq[0];
                                default_seen = true;
                        }
                        continue;
                }
                if (seq.nops() != 2)
                        throw std::runtime_error("cases case not a pair");
                const ex& cond = seq[0];
                if (is_exactly_a<numeric>(cond)) {
                        if (not cond.is_zero()
                            and not default_seen) {
                                deflt = seq[1];
                                default_seen = true;
                        }
                        continue;
                }
                if (is_exactly_a<relational>(cond)) {
                        auto res = ex_to<relational>(cond).decide();
                        if (res == relational::result::True)
                                return seq[1];
                        if (res == relational::result::False)
                                continue;
                        // undecidable relation
                        return cases(arg).hold();
                }
                throw std::runtime_error("cases with meaningless condition");
        }
        if (default_seen)
                return deflt;
        return cases(arg).hold();
}

static ex cases_evalf(const ex & arg, PyObject* parent)
{
        return cases(arg).hold();
}

static ex cases_conjugate(const ex & arg)
{
        return cases(arg).hold();
}

static ex cases_real_part(const ex & arg)
{
        return cases(arg).hold();
}

static ex cases_imag_part(const ex & arg)
{
        return cases(arg).hold();
}

static ex cases_subs(const exmap& m, const ex & arg) {
        bool default_seen = false;
        ex deflt;
        for (auto case_ : arg) {
                if (not is_exactly_a<exprseq>(case_)) {
                        if (not default_seen) {
                                deflt = case_;
                                default_seen = true;
                        }
                        continue;
                }
                const exprseq& seq = ex_to<exprseq>(case_);
                if (seq.nops() == 1) {
                        if (not default_seen) {
                                deflt = seq[0];
                                default_seen = true;
                        }
                        continue;
                }
                if (seq.nops() != 2)
                        throw std::runtime_error("cases case not a pair");
                const ex& cond = seq[0].subs(m);
                if (is_exactly_a<numeric>(cond)) {
                        if (not cond.is_zero()
                            and not default_seen) {
                                deflt = seq[1];
                                default_seen = true;
                        }
                        continue;
                }
                if (is_exactly_a<relational>(cond)) {
                        auto res = ex_to<relational>(cond).decide();
                        if (res == relational::result::True)
                                return seq[1].subs(m);
                        if (res == relational::result::False)
                                continue;
                        // undecidable relation
                        return cases(arg).hold();
                }
                throw std::runtime_error("cases with meaningless condition");
        }
        if (default_seen)
                return deflt.subs(m);
        return cases(arg).hold();
}

REGISTER_FUNCTION(cases, eval_func(cases_eval).
                       evalf_func(cases_evalf).
                       conjugate_func(cases_conjugate).
                       real_part_func(cases_real_part).
                       imag_part_func(cases_imag_part).
                       subs_func(cases_subs));

static std::string ex_to_str(const ex& e) {
        std::stringstream tstream;
        print_context c(tstream);
        e.print(c);
        return tstream.str();
}

static ex set_of_all_eval(const ex& arg1, const ex& arg2)
{
        if (is_exactly_a<function>(arg2)) {
                std::string name = ex_to<function>(arg2).get_name();
                if (not (name == "element_of"))
                        throw std::runtime_error(std::string("not a set membership statement: ") + ex_to_str(arg2));
                return set_of_all(arg1, arg2).hold();
        }
        if (not is_exactly_a<exprseq>(arg2))
                throw std::runtime_error(std::string("not a set membership statement: ") + ex_to_str(arg2));
        for (auto item : arg2) {
                if (not is_exactly_a<function>(item)
                    or not (ex_to<function>(item).get_name() == "element_of"))
                        throw std::runtime_error(std::string("not a set membership statement: ") + ex_to_str(item));
        }
        return set_of_all(arg1, arg2).hold();
}

REGISTER_FUNCTION(set_of_all, eval_func(set_of_all_eval));

} // namespace GiNaC
