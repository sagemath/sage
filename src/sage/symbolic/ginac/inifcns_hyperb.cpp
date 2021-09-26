/** @file inifcns_hyperb.cpp
 *
 *  Implementation of hyperbolic functions.
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

#include "inifcns.h"
#include "ex.h"
#include "ex_utils.h"
#include "constant.h"
#include "infinity.h"
#include "numeric.h"
#include "mul.h"
#include "add.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "symbol.h"
#include "pseries.h"
#include "utils.h"

#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <memory>

namespace GiNaC {


/* In Sage all the area hyperbolic functions are printed with "arc" instead
   of "a" at the beginning.   This is for consistency with other
   computer algebra systems.   These print methods are registered
   below with each of the corresponding inverse trig function. */

// helper function: returns whether the expression is a multiple of I
static bool is_multiple_of_I(const ex & the_ex)
{

	if (is_exactly_a<numeric>(the_ex)
	    and the_ex.real_part().is_zero())
		return true;

	if (is_exactly_a<mul>(the_ex)) {
		for (size_t i=0; i < the_ex.nops(); ++i)
			if (is_multiple_of_I(the_ex.op(i)))
				return true;
	}

	if (is_exactly_a<add>(the_ex)) {
		for (size_t i=0; i < the_ex.nops(); ++i)
			if (!is_multiple_of_I(the_ex.op(i)))
				return false;
		return true;
	}
	return false;
};


//////////
// hyperbolic sine (trigonometric function)
//////////

static ex sinh_eval(const ex & x)
{
        // sinh() is odd
        if (x.info(info_flags::negative))
                return -sinh(-x);

	if (is_exactly_a<numeric>(x)) {

		// sinh(0) -> 0
		if (x.is_zero())
			return _ex0;        

		// sinh(float) -> float
		if (x.info(info_flags::inexact))
			return sinh(ex_to<numeric>(x));
	}

	// sinh(oo) -> oo
	// sinh(-oo) -> -oo
	// sinh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(UnsignedInfinity))
			throw (std::runtime_error("sinh_eval(): sinh(unsigned_infinity) encountered"));
		return x;
	}

	// sinh(I*x) --> I*sin(x/I)
	if (is_multiple_of_I(x.expand()))
	    return I*sin(x/I);
	
	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

                // sinh(log(x)) -> (x^2 - 1)/(2x)
		if (is_ex_the_function(x, log))
                        return (power(t, _ex2) - _ex1)/(_ex2*t);

		// sinh(asinh(x)) -> x
		if (is_ex_the_function(x, asinh))
			return t;

		// sinh(acosh(x)) -> sqrt(x-1) * sqrt(x+1)
		if (is_ex_the_function(x, acosh))
			return sqrt(t-_ex1)*sqrt(t+_ex1);

		// sinh(atanh(x)) -> x/sqrt(1-x^2)
		if (is_ex_the_function(x, atanh))
			return t*power(_ex1-power(t,_ex2),_ex_1_2);
	}
	
	return sinh(x).hold();
}

static ex sinh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx sinh(x) -> cosh(x)
	return cosh(x);
}

static ex sinh_real_part(const ex & x)
{
	return sinh(GiNaC::real_part(x))*cos(GiNaC::imag_part(x));
}

static ex sinh_imag_part(const ex & x)
{
	return cosh(GiNaC::real_part(x))*sin(GiNaC::imag_part(x));
}

static ex sinh_conjugate(const ex & x)
{
	// conjugate(sinh(x))==sinh(conjugate(x))
	return sinh(x.conjugate());
}

REGISTER_FUNCTION(sinh, eval_func(sinh_eval).
                        derivative_func(sinh_deriv).
                        real_part_func(sinh_real_part).
                        imag_part_func(sinh_imag_part).
			conjugate_func(sinh_conjugate).
		        latex_name("\\sinh"));

//////////
// hyperbolic cosine (trigonometric function)
//////////

static ex cosh_evalf(const ex & x, PyObject* parent)
{
	if (is_exactly_a<numeric>(x))
		return cosh(ex_to<numeric>(x));
	
	return cosh(x).hold();
}

static ex cosh_eval(const ex & x)
{
        // cosh() is even
        if (x.info(info_flags::negative))
                return cosh(-x);
	
        if (is_exactly_a<numeric>(x)) {
		// cosh(0) -> 1
		if (x.is_zero())
			return _ex1;

		// cosh(float) -> float
		if (x.info(info_flags::inexact))
			return cosh(ex_to<numeric>(x));
	}
	
	// cosh(oo) -> oo
	// cosh(-oo) -> oo
	// cosh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(UnsignedInfinity))
			throw (std::runtime_error("cosh_eval(): cosh(unsigned_infinity) encountered"));
		return Infinity;
	}

	// cosh(I*x) --> cos(x)
	if (is_multiple_of_I(x.expand()))
		return cos(x/I);
	
	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

                // cosh(log(x)) -> (x + 1/x)/2
		if (is_ex_the_function(x, log))
                        return (power(t, _ex2) + _ex1)/(_ex2*t);

		// cosh(acosh(x)) -> x
		if (is_ex_the_function(x, acosh))
			return t;

		// cosh(asinh(x)) -> sqrt(1+x^2)
		if (is_ex_the_function(x, asinh))
			return sqrt(_ex1+power(t,_ex2));

		// cosh(atanh(x)) -> 1/sqrt(1-x^2)
		if (is_ex_the_function(x, atanh))
			return power(_ex1-power(t,_ex2),_ex_1_2);
	}
	
	return cosh(x).hold();
}

static ex cosh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx cosh(x) -> sinh(x)
	return sinh(x);
}

static ex cosh_real_part(const ex & x)
{
	return cosh(GiNaC::real_part(x))*cos(GiNaC::imag_part(x));
}

static ex cosh_imag_part(const ex & x)
{
	return sinh(GiNaC::real_part(x))*sin(GiNaC::imag_part(x));
}

static ex cosh_conjugate(const ex & x)
{
	// conjugate(cosh(x))==cosh(conjugate(x))
	return cosh(x.conjugate());
}

REGISTER_FUNCTION(cosh, eval_func(cosh_eval).
                        evalf_func(cosh_evalf).
                        derivative_func(cosh_deriv).
                        real_part_func(cosh_real_part).
                        imag_part_func(cosh_imag_part).
                        conjugate_func(cosh_conjugate).
                        latex_name("\\cosh"));

//////////
// hyperbolic tangent (trigonometric function)
//////////

static ex tanh_eval(const ex & x)
{
        // tanh() is odd
        if (x.info(info_flags::negative))
                return -tanh(-x);

	if (is_exactly_a<numeric>(x)) {
		// tanh(0) -> 0
		if (x.is_zero())
			return _ex0;

		// tanh(float) -> float
		if (x.info(info_flags::inexact))
			return tanh(ex_to<numeric>(x));
	}
	
	// tanh(oo) -> 1
	// tanh(-oo) -> -1
	// tanh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity))
			return _ex1;
		if (x.is_equal(NegInfinity))
			return _ex_1;
		// x is UnsignedInfinity
		throw (std::runtime_error("tanh_eval(): tanh(unsigned_infinity) encountered"));
	}
		
	// tanh(I*x) --> I*tan(x)
	if (is_multiple_of_I(x.expand()))
		return I*tan(x/I);
	
	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

                // tanh(log(x)) -> (x^2 - 1)/(x^2 + 1)
		if (is_ex_the_function(x, log))
                        return (power(t, _ex2) - _ex1)/(power(t, _ex2) + _ex1);

		// tanh(atanh(x)) -> x
		if (is_ex_the_function(x, atanh))
			return t;

		// tanh(asinh(x)) -> x/sqrt(1+x^2)
		if (is_ex_the_function(x, asinh))
			return t*power(_ex1+power(t,_ex2),_ex_1_2);

		// tanh(acosh(x)) -> sqrt(x-1)*sqrt(x+1)/x
		if (is_ex_the_function(x, acosh))
			return sqrt(t-_ex1)*sqrt(t+_ex1)*power(t,_ex_1);
	}
	
	return tanh(x).hold();
}

static ex tanh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx tanh(x) -> 1-tanh(x)^2
	return _ex1-power(tanh(x),_ex2);
}

static ex tanh_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole falls back to tanh_deriv.
	// On a pole simply expand sinh(x)/cosh(x).
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (!(2*I*x_pt/Pi).info(info_flags::odd))
		throw do_taylor();  // caught by function::series()
	// if we got here we have to care for a simple pole
	return (sinh(x)/cosh(x)).series(rel, order, options);
}

// See http://dlmf.nist.gov/4.35.E36
static ex tanh_real_part(const ex & x)
{
	ex a = mul(real_part(x), _ex2);
	ex b = mul(imag_part(x), _ex2);
	return sinh(a) / (cosh(a) + cos(b));
}

static ex tanh_imag_part(const ex & x)
{
	ex a = mul(real_part(x), _ex2);
	ex b = mul(imag_part(x), _ex2);
	return sin(b) / (cosh(a) + cos(b));
}

static ex tanh_conjugate(const ex & x)
{
	// conjugate(tanh(x))==tanh(conjugate(x))
	return tanh(x.conjugate());
}

REGISTER_FUNCTION(tanh, eval_func(tanh_eval).
                        derivative_func(tanh_deriv).
                        series_func(tanh_series).
                        real_part_func(tanh_real_part).
                        imag_part_func(tanh_imag_part).
                        conjugate_func(tanh_conjugate).
                        latex_name("\\tanh"));

//////////
// hyperbolic cotangent (trigonometric function)
//////////

static ex coth_eval(const ex & x)
{
        // coth() is odd
        if (x.info(info_flags::negative))
                return -coth(-x);

	if (is_exactly_a<numeric>(x)) {
		// coth(0) -> zoo
		if (x.is_zero())
			return UnsignedInfinity;

		// coth(float) -> float
		if (x.info(info_flags::inexact))
			return tanh(ex_to<numeric>(x)).inverse();
	}

	// coth(oo) -> 1
	// coth(-oo) -> -1
	// coth(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity))
			return _ex1;
		if (x.is_equal(NegInfinity))
			return _ex_1;
		// x is UnsignedInfinity
		throw (std::runtime_error("coth_eval(): tanh(unsigned_infinity) encountered"));
	}

	// coth(I*x) --> -I*cot(x)
	if (is_multiple_of_I(x.expand()))
		return -I*cot(x/I);

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

                // coth(log(x)) -> (x^2 + 1)/(x^2 - 1)
		if (is_ex_the_function(x, log))
                        return (power(t, _ex2) + _ex1)/(power(t, _ex2) - _ex1);

		// coth(acoth(x)) -> x
		if (is_ex_the_function(x, acoth))
			return t;

		// coth(asinh(x)) -> sqrt(1+x^2)/x
		if (is_ex_the_function(x, asinh))
			return power(_ex1+power(t,_ex2),_ex1_2)/t;

		// coth(acosh(x)) -> x/(sqrt(x-1)*sqrt(x+1))
		if (is_ex_the_function(x, acosh))
			return t/sqrt(t-_ex1)/sqrt(t+_ex1);
	}

	return coth(x).hold();
}

static ex coth_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx tanh(x) -> 1-tanh(x)^2
	return -power(sinh(x),_ex_2);
}

static ex coth_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole falls back to tanh_deriv.
	// On a pole simply expand sinh(x)/cosh(x).
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (!(2*I*x_pt/Pi).info(info_flags::even))
		throw do_taylor();  // caught by function::series()
	// if we got here we have to care for a simple pole
	return (cosh(x)/sinh(x)).series(rel, order, options);
}

static ex coth_real_part(const ex & x)
{
	ex a = real_part(x);
	ex b = imag_part(x);
	return mul(sinh(a), cosh(a)) / (power(sin(b), _ex2) + power(sinh(a), _ex2));
}

static ex coth_imag_part(const ex & x)
{
	ex a = real_part(x);
	ex b = imag_part(x);
	return -mul(sin(b), cos(b)) / (power(sin(b), _ex2) + power(sinh(a), _ex2));
}

static ex coth_conjugate(const ex & x)
{
	return coth(x.conjugate());
}

REGISTER_FUNCTION(coth, eval_func(coth_eval).
                        derivative_func(coth_deriv).
                        series_func(coth_series).
                        real_part_func(coth_real_part).
                        imag_part_func(coth_imag_part).
                        conjugate_func(coth_conjugate).
                        latex_name("\\coth"));

//////////
// hyperbolic secant (trigonometric function)
//////////

static ex sech_eval(const ex & x)
{
        // sech() is even
        if (x.info(info_flags::negative))
                return sech(-x);

	if (is_exactly_a<numeric>(x)) {
		// sech(0) -> 1
		if (x.is_zero())
			return _ex1;

		// sech(float) -> float
		if (x.info(info_flags::inexact))
			return cosh(ex_to<numeric>(x)).inverse();
	}

	// sech(x*I) --> sec(x)
	if (is_multiple_of_I(x.expand()))
	    return sec(x/I);

	// sech(oo) -> 0
	// sech(-oo) -> 0
	// sech(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity) or x.is_equal(NegInfinity))
			return _ex0;
		// x is UnsignedInfinity
		throw (std::runtime_error("sech_eval(): sech(unsigned_infinity) encountered"));
	}

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

                // sech(log(x)) -> 2/(x + 1/x)
		if (is_ex_the_function(x, log))
                        return (_ex2*t) / (power(t, _ex2) + _ex1);

		// sech(asech(x)) -> x
		if (is_ex_the_function(x, asech))
			return t;

		// sech(asinh(x)) -> 1/sqrt(1+x^2)
		if (is_ex_the_function(x, asinh))
			return power(_ex1+power(t,_ex2),_ex_1_2);

		// sech(acosh(x)) -> 1/x
		if (is_ex_the_function(x, acosh))
			return power(t, _ex_1);
	}

	return sech(x).hold();
}

static ex sech_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx sech(x) -> -sech(x)*tanh(x)
	return -mul(sech(x), tanh(x));
}

static ex sech_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole falls back to tanh_deriv.
	// On a pole simply expand sinh(x)/cosh(x).
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (!(2*I*x_pt/Pi).info(info_flags::odd))
		throw do_taylor();  // caught by function::series()
	// if we got here we have to care for a simple pole
	return (_ex1/cosh(x)).series(rel, order, options);
}

static ex sech_real_part(const ex & x)
{
	ex a = real_part(x);
	ex b = imag_part(x);
	return mul(cos(b), cosh(a)) / (power(mul(sin(b), sinh(a)), _ex2) + power(mul(cos(b), cosh(a)), _ex2));
}

static ex sech_imag_part(const ex & x)
{
	ex a = real_part(x);
	ex b = imag_part(x);
	return -mul(sin(b), sinh(a)) / (power(mul(sin(b), sinh(a)), _ex2) + power(mul(cos(b), cosh(a)), _ex2));
}

static ex sech_conjugate(const ex & x)
{
	return sech(x.conjugate());
}

REGISTER_FUNCTION(sech, eval_func(sech_eval).
                        derivative_func(sech_deriv).
                        series_func(sech_series).
                        real_part_func(sech_real_part).
                        imag_part_func(sech_imag_part).
                        conjugate_func(sech_conjugate).
                        latex_name("\\operatorname{sech}"));

//////////
// hyperbolic secant (trigonometric function)
//////////

static ex csch_eval(const ex & x)
{
        // csch() is odd
        if (x.info(info_flags::negative))
                return -csch(-x);

	if (is_exactly_a<numeric>(x)) {
		// csch(0) -> zoo
		if (x.is_zero())
			return UnsignedInfinity;

		// csch(float) -> float
		if (x.info(info_flags::inexact))
			return sinh(ex_to<numeric>(x)).inverse();
	}

	// csch(I*x) --> -I*csc(x)
	if (is_multiple_of_I(x.expand()))
		return -I*csc(x/I);

	// csch(oo) -> 0
	// csch(-oo) -> 0
	// csch(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity) or x.is_equal(NegInfinity))
			return _ex0;
		// x is UnsignedInfinity
		throw (std::runtime_error("csch_eval(): csch(unsigned_infinity) encountered"));
	}

	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

                // csch(log(x)) -> 2/(x - 1/x)
		if (is_ex_the_function(x, log))
                        return (_ex2*t) / (power(t, _ex2) - _ex1);

		// csch(acsch(x)) -> x
		if (is_ex_the_function(x, acsch))
			return t;

		// csch(acosh(x)) -> 1/sqrt(x^2-1)
		if (is_ex_the_function(x, asinh))
			return power(power(t,_ex2)-_ex1,_ex_1_2);

		// csch(asinh(x)) -> 1/x
		if (is_ex_the_function(x, asinh))
			return power(t, _ex_1);
	}

	return csch(x).hold();
}

static ex csch_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);

	// d/dx sech(x) -> -sech(x)*tanh(x)
	return -mul(csch(x), coth(x));
}

static ex csch_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole falls back to tanh_deriv.
	// On a pole simply expand sinh(x)/cosh(x).
	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (!(2*I*x_pt/Pi).info(info_flags::even))
		throw do_taylor();  // caught by function::series()
	// if we got here we have to care for a simple pole
	return (_ex1/sinh(x)).series(rel, order, options);
}

static ex csch_real_part(const ex & x)
{
	ex a = real_part(x);
	ex b = imag_part(x);
	return mul(cos(b), sinh(a)) / (power(mul(sin(b), cosh(a)), _ex2) + power(mul(cos(b), sinh(a)), _ex2));
}

static ex csch_imag_part(const ex & x)
{
	ex a = real_part(x);
	ex b = imag_part(x);
	return -mul(sin(b), cosh(a)) / (power(mul(sin(b), cosh(a)), _ex2) + power(mul(cos(b), sinh(a)), _ex2));
}

static ex csch_conjugate(const ex & x)
{
	return csch(x.conjugate());
}

REGISTER_FUNCTION(csch, eval_func(csch_eval).
                        derivative_func(csch_deriv).
                        series_func(csch_series).
                        real_part_func(csch_real_part).
                        imag_part_func(csch_imag_part).
                        conjugate_func(csch_conjugate).
                        latex_name("\\operatorname{csch}"));

//////////
// inverse hyperbolic sine (trigonometric function)
//////////

static ex asinh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// asinh(0) -> 0
		if (x.is_zero())
			return _ex0;

		// asinh(float) -> float
		if (x.info(info_flags::inexact))
			return asinh(ex_to<numeric>(x));

		// asinh() is odd
		if (x.info(info_flags::negative))
			return -asinh(-x);
	}
	
	// asinh(oo) -> oo
	// asinh(-oo) -> -oo
	// asinh(UnsignedInfinity) -> UnsignedInfinity
	if (x.info(info_flags::infinity)) {
		return x;
	}

	return asinh(x).hold();
}

static ex asinh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx asinh(x) -> 1/sqrt(1+x^2)
	return power(_ex1+power(x,_ex2),_ex_1_2);
}

static ex asinh_conjugate(const ex & x)
{
	// conjugate(asinh(x))==asinh(conjugate(x)) unless on the branch cuts which
	// run along the imaginary axis outside the interval [-I, +I].
	if (x.is_real())
		return asinh(x);
	if (is_exactly_a<numeric>(x)) {
		const numeric x_re = ex_to<numeric>(x.real_part());
		const numeric x_im = ex_to<numeric>(x.imag_part());
		if (!x_re.is_zero() ||
		    (x_im > *_num_1_p && x_im < *_num1_p))
			return asinh(x.conjugate());
	}
	return conjugate_function(asinh(x)).hold();
}

REGISTER_FUNCTION(asinh, eval_func(asinh_eval).
                         derivative_func(asinh_deriv).
                         conjugate_func(asinh_conjugate).
			 set_name("arcsinh"));

//////////
// inverse hyperbolic cosine (trigonometric function)
//////////

static ex acosh_eval(const ex & x)
{
	if (is_exactly_a<numeric>(x)) {

		// acosh(0) -> Pi*I/2
		if (x.is_zero())
			return Pi*I*numeric(1,2);

		// acosh(1) -> 0
		if (x.is_one())
			return _ex0;

		// acosh(-1) -> Pi*I
		if (x.is_minus_one())
			return Pi*I;

		// acosh(float) -> float
		if (x.info(info_flags::inexact))
			return acosh(ex_to<numeric>(x));

		// acosh(-x) -> Pi*I-acosh(x)
		if (x.info(info_flags::negative))
			return Pi*I-acosh(-x);
	}
	
	// acosh(oo) -> oo
	// acosh(-oo) -> oo
	// acosh(UnsignedInfinity) -> oo
	if (x.info(info_flags::infinity)) {
		return Infinity;
	}

	return acosh(x).hold();
}

static ex acosh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx acosh(x) -> 1/(sqrt(x-1)*sqrt(x+1))
	return power(x+_ex_1,_ex_1_2)*power(x+_ex1,_ex_1_2);
}

static ex acosh_conjugate(const ex & x)
{
	// conjugate(acosh(x))==acosh(conjugate(x)) unless on the branch cut
	// which runs along the real axis from +1 to -inf.
	if (is_exactly_a<numeric>(x) &&
	    (!x.imag_part().is_zero() || x > *_num1_p)) {
		return acosh(x.conjugate());
	}
	return conjugate_function(acosh(x)).hold();
}

REGISTER_FUNCTION(acosh, eval_func(acosh_eval).
                         derivative_func(acosh_deriv).
                         conjugate_func(acosh_conjugate).
			 set_name("arccosh"));

//////////
// inverse hyperbolic tangent (trigonometric function)
//////////

static ex atanh_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {

                const numeric& num = ex_to<numeric>(x);

		// atanh(0) -> 0
		if (num.is_zero())
			return _ex0;

		/*
		// atanh({+|-}1) -> throw
		if (x.is_equal(_ex1) || x.is_equal(_ex_1))
			throw (pole_error("atanh_eval(): logarithmic pole",0));
		*/

		// atanh(1) -> oo
		if (num.is_one())
			return Infinity;
		// atahn(-1) -> -oo
		if (num.is_minus_one())
			return NegInfinity;

		// atanh() is odd
		if (num.is_negative())
			return -atanh(-x);

		// atanh(float) -> float
		if (not num.is_exact())
			return atanh(num);
                
                if (num.is_integer() or num.is_rational())
                        return _ex1_2 * log(ex((*_num1_p + num)/(*_num1_p - num)));

                return atanh(x).hold();
	}
	
	// atanh(oo) -> -i*pi/2
	// atanh(-oo) -> i*pi/2
	// atanh(UnsignedInfinity) -> error
	if (x.info(info_flags::infinity)) {
		if (x.is_equal(Infinity))
			return _ex_1_2*Pi*I;
		if (x.is_equal(NegInfinity))
			return _ex1_2*Pi*I;
		// x is UnsignedInfinity
		throw (std::runtime_error("arctanh_eval(): arctanh(unsigned_infinity) encountered"));
	}
		
	return atanh(x).hold();
}

static ex atanh_deriv(const ex & x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param==0);
	
	// d/dx atanh(x) -> 1/(1-x^2)
	return power(_ex1-power(x,_ex2),_ex_1);
}

static ex atanh_series(const ex &arg,
                       const relational &rel,
                       int order,
                       unsigned options)
{
	GINAC_ASSERT(is_a<symbol>(rel.lhs()));
	// method:
	// Taylor series where there is no pole or cut falls back to atanh_deriv.
	// There are two branch cuts, one running from 1 up the real axis and one
	// one running from -1 down the real axis.  The points 1 and -1 are poles
	// On the branch cuts and the poles series expand
	//     (log(1+x)-log(1-x))/2
	// instead.
	const ex arg_pt = arg.subs(rel, subs_options::no_pattern);
	if (!(arg_pt).is_real())
		throw do_taylor();     // Im(x) != 0
	if ((arg_pt).is_real() && abs(arg_pt)<_ex1)
		throw do_taylor();     // Im(x) == 0, but abs(x)<1
	// care for the poles, using the defining formula for atanh()...
	if (arg_pt.is_equal(_ex1) || arg_pt.is_equal(_ex_1))
		return ((log(_ex1+arg)-log(_ex1-arg))*_ex1_2).series(rel, order, options);
	// ...and the branch cuts (the discontinuity at the cut being just I*Pi)
	if ((options & series_options::suppress_branchcut) == 0u) {
 		// method:
 		// This is the branch cut: assemble the primitive series manually and
 		// then add the corresponding complex step function.
 		const symbol &s = ex_to<symbol>(rel.lhs());
 		const ex &point = rel.rhs();
 		const symbol foo;
 		const ex replarg = series(atanh(arg), s==foo, order).subs(foo==point, subs_options::no_pattern);
		ex Order0correction = replarg.op(0)+csgn(I*arg)*Pi*I*_ex1_2;
		if (arg_pt<_ex0)
			Order0correction += log((arg_pt+_ex_1)/(arg_pt+_ex1))*_ex1_2;
		else
			Order0correction += log((arg_pt+_ex1)/(arg_pt+_ex_1))*_ex_1_2;
 		epvector seq;
		seq.emplace_back(Order0correction, _ex0);
 		seq.emplace_back(Order(_ex1), order);
 		return series(replarg - pseries(rel, seq), rel, order);
	}
	throw do_taylor();
}

static ex atanh_conjugate(const ex & x)
{
	// conjugate(atanh(x))==atanh(conjugate(x)) unless on the branch cuts which
	// run along the real axis outside the interval [-1, +1].
	if (is_exactly_a<numeric>(x) &&
	    (!x.imag_part().is_zero() || (x > *_num_1_p && x < *_num1_p))) {
		return atanh(x.conjugate());
	}
	return conjugate_function(atanh(x)).hold();
}

REGISTER_FUNCTION(atanh, eval_func(atanh_eval).
                         derivative_func(atanh_deriv).
                         series_func(atanh_series).
                         conjugate_func(atanh_conjugate).
			 set_name("arctanh"));

//////////
// inverse hyperbolic cotangent (trigonometric function)
//////////

static ex acoth_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                const numeric& num = ex_to<numeric>(x);
                // acoth(1) -> oo
                if (num.is_one())
                        return Infinity;
                // acoth(-1) -> -oo
                if (num.is_minus_one())
                        return NegInfinity;
                //acoth(float) -> float 
                if (not num.is_exact())
                        return atanh(num.inverse());
                // acoth() is odd
                if (num.is_negative())
                        return -acoth(num.negative());
                if (num.is_integer() or num.is_rational())
                        return _ex1_2 * log(ex((num + *_num1_p)/(num - *_num1_p)));
        }
       
	if (is_exactly_a<function>(x)) {
		const ex &t = x.op(0);

		// acoth(coth(x)) -> x
		if (is_ex_the_function(x, coth))
			return t;
	}

        // acoth(oo) -> 0
        // acoth(-oo) -> 0
        // acoth(UnsignedInfinity) -> 0
        if (x.info(info_flags::infinity)) {
                return _ex0;
                // x is UnsignedInfinity
                //throw (std::runtime_error("arccoth_eval(): arccoth(unsigned_infinity) encountered"));
        }
        
        return acoth(x).hold();
}

static ex acoth_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // acoth(x) -> (1/2)*(ln(1 + 1/x) - ln(1 - 1/x))
        // d/dx acoth(x) -> 1/(1-x^2)
        return power(_ex1-power(x, _ex2), _ex_1);
}

static ex acoth_conjugate(const ex & x)
{
        // conjugate(acoth(x))==acoth(conjugate(x)) unless on the branch cuts which
        // run along the real axis inside the interval [-1, +1].
        if (is_exactly_a<numeric>(x) &&
            (!x.imag_part().is_zero() || (x < *_num_1_p && x > *_num1_p))) {
                return acoth(x.conjugate());
        }
        return conjugate_function(acoth(x)).hold();
}

REGISTER_FUNCTION(acoth, eval_func(acoth_eval).
                         derivative_func(acoth_deriv).
                         conjugate_func(acoth_conjugate).
                         set_name("arccoth"));

//////////
// inverse hyperbolic Cosecant (trigonometric function)
//////////

static ex acsch_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                // acsch(0) -> oo
                if (x.is_zero())
                        return Infinity;
                //acsch(float) -> float 
                if (x.info(info_flags::inexact))
                        return asinh(ex_to<numeric>(x).inverse());
                // acsch(-x) -> acsch(-x)
                if (x.info(info_flags::negative))
                        return -acsch(-x);
        }
       
        // acsch(oo) -> 0
        // acsch(-oo) -> 0
        // acsch(UnsignedInfinity) -> 0
        if (x.info(info_flags::infinity)) {
                return _ex0;
        }
        
        return acsch(x).hold();
}

static ex acsch_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // acsch(x) -> ln(1/x + sqrt(1/x^2 + 1))
        // d/dx acsch(x) ->  -1 / [x * sqrt(1 + x^2)];
        return (_ex_1/x)*power(_ex1+power(x, _ex2), _ex_1_2);
}

static ex acsch_conjugate(const ex & x)
{
        // conjugate(acsch(x))==acsch(conjugate(x)) unless on the branch cuts which
        // run along the imaginary axis inside the interval [-I, +I].
        if (x.is_real())
		return acsch(x);
	if (is_exactly_a<numeric>(x)) {
		const numeric x_re = ex_to<numeric>(x.real_part());
		const numeric x_im = ex_to<numeric>(x.imag_part());
		if (!x_re.is_zero() ||
		    (x_im < *_num_1_p && x_im > *_num1_p))
			return acsch(x.conjugate());
	}
	return conjugate_function(acsch(x)).hold();
}

REGISTER_FUNCTION(acsch, eval_func(acsch_eval).
                         derivative_func(acsch_deriv).
                         conjugate_func(acsch_conjugate).
                         set_name("arccsch"));
//////////
// inverse hyperbolic Secant (trigonometric function)
//////////

static ex asech_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                // asech(0) -> oo
                if (x.is_zero())
                        return Infinity;
                // asech(1) -> 0
                if (x.is_one())
                        return _ex0;
                //asech(-1) -> I*Pi
                if (x.is_minus_one())
                        return Pi*I;
                //asech(float) -> float 
                if (x.info(info_flags::inexact))
                        return acosh(ex_to<numeric>(x).inverse());
                // asech(-x) -> Pi*I-asech(-x)
                if (x.info(info_flags::negative))
                        return Pi*I-asech(-x);
        }
       
        // asech(oo) -> Pi*I/2
        // asech(-oo) -> Pi*I/2
        // asech(UnsignedInfinity) -> error
        if (x.info(info_flags::infinity)) {
                if (x.is_equal(Infinity) || x.is_equal(NegInfinity))
			return Pi*I*numeric(1,2);
                // x is UnsignedInfinity
                throw (std::runtime_error("arcsech_eval(): arcsech(unsigned_infinity) encountered"));
        }
        
        return asech(x).hold();
}

static ex asech_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // asech(x) -> ln(1/x + sqrt(1/x^2 - 1))
        // d/dx asech(x) ->  -1 / [x * sqrt(1 - x^2)];
        return (_ex_1/x)*power(_ex1-power(x, _ex2), _ex_1_2);
}

static ex asech_conjugate(const ex & x)
{
        // conjugate(asech(x))==asech(conjugate(x)) unless on the branch cuts which
        // run along the real axis from 0 to -oo and 1 to oo.
        if (is_exactly_a<numeric>(x) &&
            (!x.imag_part().is_zero() || (x < *_num1_p && x > *_num0_p))) {
                return asech(x.conjugate());
        }
        return conjugate_function(asech(x)).hold();
}

REGISTER_FUNCTION(asech, eval_func(asech_eval).
                         derivative_func(asech_deriv).
                         conjugate_func(asech_conjugate).
                         set_name("arcsech"));

} // namespace GiNaC
