/** @file mpoly-ginac.cpp
 *
 *  This file implements several functions that work multivariate polynomials.
*/
/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *                  (C) 2016 Ralf Stephan
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

#include "ex.h"
#include "normal.h"
#include "mpoly.h"
#include "upoly.h"
#include "numeric.h"
#include "operators.h"
#include "utils.h"
#include "symbol.h"
#include "power.h"
#include "add.h"
#include "mul.h"

namespace GiNaC {


/*
 *  Computation of LCM of denominators of coefficients of a polynomial
 */

// Compute LCM of denominators of coefficients by going through the
// expression recursively (used internally by lcm_of_coefficients_denominators())
static numeric lcmcoeff(const ex &e, const numeric &l)
{
	if (is_exactly_a<numeric>(e)
            and e.info(info_flags::rational))
		return lcm(ex_to<numeric>(e).denom(), l);
	if (is_exactly_a<add>(e)) {
		numeric c = *_num1_p;
		for (size_t i=0; i<e.nops(); i++)
			c = lcmcoeff(e.op(i), c);
		return lcm(c, l);
	}
        if (is_exactly_a<mul>(e)) {
		numeric c = *_num1_p;
		for (size_t i=0; i<e.nops(); i++)
			c *= lcmcoeff(e.op(i), *_num1_p);
		return lcm(c, l);
	}
        if (is_exactly_a<power>(e)) {
		if (is_exactly_a<symbol>(e.op(0)))
			return l;

                ex t = pow(lcmcoeff(e.op(0), l), ex_to<numeric>(e.op(1)));
                if (is_exactly_a<numeric>(t))
                        return ex_to<numeric>(t);
	}
	return l;
}

/** Compute LCM of denominators of coefficients of a polynomial.
 *  Given a polynomial with rational coefficients, this function computes
 *  the LCM of the denominators of all coefficients. This can be used
 *  to bring a polynomial from Q[X] to Z[X].
 *
 *  @param e  multivariate polynomial (need not be expanded)
 *  @return LCM of denominators of coefficients */
numeric lcm_of_coefficients_denominators(const ex &e)
{
	return lcmcoeff(e, *_num1_p);
}

/** Bring polynomial from Q[X] to Z[X] by multiplying in the previously
 *  determined LCM of the coefficient's denominators.
 *
 *  @param e  multivariate polynomial (need not be expanded)
 *  @param lcm  LCM to multiply in */
ex multiply_lcm(const ex &e, const numeric &lcm)
{
	if (is_exactly_a<mul>(e)) {
		size_t num = e.nops();
		exvector v; v.reserve(num + 1);
		numeric lcm_accum = *_num1_p;
		for (size_t i=0; i<num; i++) {
			numeric op_lcm = lcmcoeff(e.op(i), *_num1_p);
			v.push_back(multiply_lcm(e.op(i), op_lcm));
			lcm_accum *= op_lcm;
		}
		v.emplace_back(lcm / lcm_accum);
		return (new mul(v))->setflag(status_flags::dynallocated);
	}
        if (is_exactly_a<add>(e)) {
		size_t num = e.nops();
		exvector v; v.reserve(num);
		for (size_t i=0; i<num; i++)
			v.push_back(multiply_lcm(e.op(i), lcm));
		return (new add(v))->setflag(status_flags::dynallocated);
	}
        if (is_exactly_a<power>(e)) {
		if (is_exactly_a<symbol>(e.op(0)))
			return e * lcm;
                if (not is_exactly_a<numeric>(e.op(1)))
                        return e * lcm;
                ex t = lcm.power(ex_to<numeric>(e.op(1)).inverse());
                if (not is_exactly_a<numeric>(t))
                        return e * lcm;
                const numeric& root_of_lcm = ex_to<numeric>(t);
                if (root_of_lcm.is_rational())
                        return pow(multiply_lcm(e.op(0), root_of_lcm), e.op(1));
	}
        return e * lcm;
}

/*
 *  Separation of unit part, content part and primitive part of polynomials
 */

/** Compute unit part (= sign of leading coefficient) of a multivariate
 *  polynomial in Q[x]. The product of unit part, content part, and primitive
 *  part is the polynomial itself.
 *
 *  @param x  main variable
 *  @return unit part
 *  @see ex::content, ex::primpart, ex::unitcontprim */
ex ex::unit(const ex &x) const
{
	ex c = expand().lcoeff(x);
	if (is_exactly_a<numeric>(c))
		return c.info(info_flags::negative) ?_ex_1 : _ex1;

        ex y;
        if (c.get_first_symbol(y))
                return c.unit(y);

        throw(std::invalid_argument("invalid expression in unit()"));
}


/** Compute content part (= unit normal GCD of all coefficients) of a
 *  multivariate polynomial in Q[x]. The product of unit part, content part,
 *  and primitive part is the polynomial itself.
 *
 *  @param x  main variable
 *  @return content part
 *  @see ex::unit, ex::primpart, ex::unitcontprim */
ex ex::content(const ex &x) const
{
	if (is_exactly_a<numeric>(*this))
		return info(info_flags::negative) ? -*this : *this;

	if (this->is_zero())
		return _ex0;

	ex u, c, p;
	unitcontprim(x, u, c, p);
	return c;
}


/** Compute primitive part of a multivariate polynomial in Q[x]. The result
 *  will be a unit-normal polynomial with a content part of 1. The product
 *  of unit part, content part, and primitive part is the polynomial itself.
 *
 *  @param x  main variable
 *  @return primitive part
 *  @see ex::unit, ex::content, ex::unitcontprim */
ex ex::primpart(const ex &x) const
{
	// We need to compute the unit and content anyway, so call unitcontprim()
	ex u, c, p;
	unitcontprim(x, u, c, p);
	return p;
}


/** Compute primitive part of a multivariate polynomial in Q[x] when the
 *  content part is already known. This function is faster in computing the
 *  primitive part than the previous function.
 *
 *  @param x  main variable
 *  @param c  previously computed content part
 *  @return primitive part */
ex ex::primpart(const ex &x, const ex &c) const
{
	if (is_zero() || c.is_zero())
		return _ex0;
	if (is_exactly_a<numeric>(*this))
		return _ex1;

	// Divide by unit and content to get primitive part
	ex u = unit(x);
	if (is_exactly_a<numeric>(c))
		return *this / (c * u);
	return quo(*this, c * u, x, false);
}


/** Compute unit part, content part, and primitive part of a multivariate
 *  polynomial in Q[x]. The product of the three parts is the polynomial
 *  itself.
 *
 *  @param x  main variable
 *  @param u  unit part (returned)
 *  @param c  content part (returned)
 *  @param p  primitive part (returned)
 *  @see ex::unit, ex::content, ex::primpart */
void ex::unitcontprim(const ex &x, ex &u, ex &c, ex &p) const
{
	// Quick check for zero (avoid expanding)
	if (is_zero()) {
		u = _ex1;
		c = p = _ex0;
		return;
	}

	// Special case: input is a number
	if (is_exactly_a<numeric>(*this)) {
		if (info(info_flags::negative)) {
			u = _ex_1;
			c = abs(ex_to<numeric>(*this));
		} else {
			u = _ex1;
			c = *this;
		}
		p = _ex1;
		return;
	}

	// Compute unit and content
	u = unit(x);

        expairvec vec;
        coefficients(x, vec);
        c = vec[0].first;
        for (const auto& pair : range(vec.begin()+1, vec.end())) {
                c = gcdpoly(c, pair.first, nullptr, nullptr, false);
        }

        p = _ex0;
        if (is_exactly_a<numeric>(c))
                for (const auto& pair : vec)
                        p += (pair.first / (c * u)) * pow(x, pair.second);
        else
                for (const auto& pair : vec)
                        p += quo(pair.first, c * u, x, false) * pow(x, pair.second);
}

ex resultant(const ex & e1, const ex & e2, const ex & s)
{
	const ex ee1 = e1.expand();
	const ex ee2 = e2.expand();
	if (!ee1.info(info_flags::polynomial) ||
	    !ee2.info(info_flags::polynomial)) {
                ex res, f1, f2;
                bool changed = factor(ee1, res);
                if (changed)
                        f1 = res;
                else
                        f1 = ee1;
                changed = factor(ee2, res);
                if (changed)
                        f2 = res;
                else
                        f2 = ee1;
                ex den1 = f1.denom();
                ex den2 = f2.denom();
                if (not den1.is_one() and den1.is_equal(den2))
                        return resultant(f1.numer(), f2.numer(), s);
		throw(std::runtime_error("resultant(): arguments must be polynomials"));
        }

        return resultantpoly(ee1, ee2, s);
}

} // namespace GiNaC

