/** @file upoly-ginac.cpp
 *
 *  This file implements several functions that work on univariate and
 *  rational functions.
 *  These functions include polynomial quotient and remainder, GCD and LCM
 *  computation, square-free factorization and rational function normalization. */

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

#include "pynac-config.h"

#include "ex.h"
#include "ex_utils.h"
#include "numeric.h"
#include "normal.h"
#include "upoly.h"
#include "symbol.h"
#include "exprseq.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "matrix.h"
#include "operators.h"
#include "utils.h"

#include <map>
#include <sstream>

namespace GiNaC {

// If comparing expressions (ex::compare()) is fast, you can set this to 1.
// Some routines like quo(), rem() and gcd() will then return a quick answer
// when they are called with two identical arguments.
#define FAST_COMPARE 1

/*
 *  Polynomial quotients and remainders
 */

/** Quotient q(x) of polynomials a(x) and b(x) in Q[x].
 *  It satisfies a(x)=b(x)*q(x)+r(x).
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return quotient of a and b in Q[x] */
ex quo(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("quo: division by zero"));
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b))
		return a / b;
#if FAST_COMPARE
	if (a.is_equal(b))
		return _ex1;
#endif
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("quo: arguments must be polynomials over the rationals"));
	// Polynomial long division
	ex r = a.expand();
	if (r.is_zero())
		return r;
	numeric bdeg = b.degree(x);
	numeric rdeg = r.degree(x);
	ex blcoeff = b.expand().coeff(x, bdeg);
	bool blcoeff_is_numeric = is_exactly_a<numeric>(blcoeff);
	exvector v; //v.reserve(std::max(rdeg - bdeg + 1, 0));
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (blcoeff_is_numeric)
			term = rcoeff / blcoeff;
		else {
			if (!divide(rcoeff, blcoeff, term, false))
                                throw std::logic_error("");
		}
		term *= power(x, rdeg - bdeg);
		v.push_back(term);
		r -= (term * b).expand();
		if (r.is_zero())
			break;
		rdeg = r.degree(x);
	}
        return (new add(v))->setflag(status_flags::dynallocated);
}

/** Remainder r(x) of polynomials a(x) and b(x) in Q[x].
 *  It satisfies a(x)=b(x)*q(x)+r(x).
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return remainder of a(x) and b(x) in Q[x] */
ex rem(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("rem: division by zero"));
	if (is_exactly_a<numeric>(a)) {
		if  (is_exactly_a<numeric>(b))
			return _ex0;
		return a;
	}
#if FAST_COMPARE
	if (a.is_equal(b))
		return _ex0;
#endif
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("rem: arguments must be polynomials over the rationals"));

	// Polynomial long division
	ex r = a.expand();
	if (r.is_zero())
		return r;
	numeric bdeg = b.degree(x);
        numeric rdeg = r.degree(x);
	ex blcoeff = b.expand().coeff(x, bdeg);
	bool blcoeff_is_numeric = is_exactly_a<numeric>(blcoeff);
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (blcoeff_is_numeric)
			term = rcoeff / blcoeff;
		else {
			if (!divide(rcoeff, blcoeff, term, false))
                                throw std::logic_error("");
		}
		term *= power(x, rdeg - bdeg);
		r -= (term * b).expand();
		if (r.is_zero())
			break;
		rdeg = r.degree(x);
	}
	return r;
}


std::pair<ex,ex> quo_rem(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("quo: division by zero"));
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b))
		return std::make_pair(a / b, _ex0);
#if FAST_COMPARE
	if (a.is_equal(b))
		return std::make_pair(_ex1, _ex0);
#endif
        expairvec vec1, vec2;
        using pair_t = std::pair<ex,ex>;
        a.expand().coefficients(x, vec1);
        b.expand().coefficients(x, vec2);
        if (check_args) {
                if (std::any_of(vec1.begin(), vec1.end(),
                             [](const pair_t& p) {
                                return not is_exactly_a<numeric>(p.second)
                                    or (not p.second.is_zero()
                                        and not p.second.info(info_flags::posint)); } 
                                    ))
                        throw std::logic_error("nonposint exponent in quo()");
                if (std::any_of(vec2.begin(), vec2.end(),
                             [](const pair_t& p) {
                                return not is_exactly_a<numeric>(p.second)
                                    or (not p.second.is_zero()
                                        and not p.second.info(info_flags::posint)); } 
                                    ))
                        throw std::logic_error("nonposint exponent in quo()");
        }
        const auto& mit1 = std::max_element(vec1.begin(), vec1.end(),
                        [](const pair_t& p1, const pair_t& p2) {
                           return ex_to<numeric>(p1.second) <
                                  ex_to<numeric>(p2.second); } );
        const auto& mit2 = std::max_element(vec2.begin(), vec2.end(),
                        [](const pair_t& p1, const pair_t& p2) {
                           return ex_to<numeric>(p1.second) <
                                  ex_to<numeric>(p2.second); } );
        size_t adeg = ex_to<numeric>(mit1->second).to_long();
        size_t bdeg = ex_to<numeric>(mit2->second).to_long();
        if (adeg < bdeg)
                return std::make_pair(_ex0, a);
        // make flat coeff vectors
        std::vector<ex> avec, bvec;
        avec.assign(adeg+1, _ex0);
        bvec.assign(bdeg+1, _ex0);
        for (const pair_t& p : vec1)
                *(avec.begin() + ex_to<numeric>(p.second).to_long()) = p.first;
        for (const pair_t& p : vec2)
                *(bvec.begin() + ex_to<numeric>(p.second).to_long()) = p.first;

        // poly division
        ex lc2 = bvec[bdeg];
        exvector qvec;
        for (int k=adeg-bdeg; k>=0; --k) {
                ex q = avec[bdeg+k] / lc2;
                avec[bdeg+k] = _ex0;
                for (int j=bdeg+k-1; j>=k; --j)
                        avec[j] = avec[j] - q * bvec[j-k];
                qvec.insert(qvec.begin(), q);
        }
        for (size_t i=0; i<qvec.size(); ++i)
                if (not qvec[i].is_zero())
                        qvec[i] = qvec[i] * power(x, numeric(static_cast<int>(i)));
        if (avec.size() > bdeg)
                avec.resize(bdeg);
        for (size_t i=0; i<bdeg; ++i)
                if (not avec[i].is_zero())
                        avec[i] = avec[i] * power(x, numeric(static_cast<int>(i)));
        return std::make_pair(add(qvec), add(avec));
}


/** Decompose rational function a(x)=N(x)/D(x) into P(x)+n(x)/D(x)
 *  with degree(n, x) < degree(D, x).
 *
 *  @param a rational function in x
 *  @param x a is a function of x
 *  @return decomposed function. */
ex decomp_rational(const ex &a, const ex &x)
{
	ex nd = numer_denom(a);
	ex numer = nd.op(0), denom = nd.op(1);
        ex q;
        try {
        	q = quo(numer, denom, x);
        }
        catch (std::logic_error) {
		return a;
        }
	return q + rem(numer, denom, x) / denom;
}


/** Pseudo-remainder of polynomials a(x) and b(x) in Q[x].
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return pseudo-remainder of a(x) and b(x) in Q[x] */
ex prem(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("prem: division by zero"));
	if (is_exactly_a<numeric>(a)) {
		if (is_exactly_a<numeric>(b))
			return _ex0;
		return b;
	}
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("prem: arguments must be polynomials over the rationals"));

	// Polynomial long division
	ex r = a.expand();
	ex eb = b.expand();
	numeric rdeg = r.degree(x);
	numeric bdeg = eb.degree(x);
	ex blcoeff;
	if (bdeg <= rdeg) {
		blcoeff = eb.coeff(x, bdeg);
		if (bdeg == 0)
			eb = _ex0;
		else
			eb -= blcoeff * power(x, bdeg);
	} else
		blcoeff = _ex1;

	numeric delta = rdeg - bdeg + 1, i = 0;
	while (rdeg >= bdeg && !r.is_zero()) {
		ex rlcoeff = r.coeff(x, rdeg);
		ex term = (power(x, rdeg - bdeg) * eb * rlcoeff).expand();
		if (rdeg == 0)
			r = _ex0;
		else
			r -= rlcoeff * power(x, rdeg);
		r = (blcoeff * r).expand() - term;
		rdeg = r.degree(x);
		i++;
	}
	return power(blcoeff, delta - i) * r;
}


/** Sparse pseudo-remainder of polynomials a(x) and b(x) in Q[x].
 *
 *  @param a  first polynomial in x (dividend)
 *  @param b  second polynomial in x (divisor)
 *  @param x  a and b are polynomials in x
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return sparse pseudo-remainder of a(x) and b(x) in Q[x] */
ex sprem(const ex &a, const ex &b, const ex &x, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("prem: division by zero"));
	if (is_exactly_a<numeric>(a)) {
		if (is_exactly_a<numeric>(b))
			return _ex0;
		return b;
	}
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("prem: arguments must be polynomials over the rationals"));

	// Polynomial long division
	ex r = a.expand();
	ex eb = b.expand();
	numeric rdeg = r.degree(x);
	numeric bdeg = eb.degree(x);
	ex blcoeff;
	if (bdeg <= rdeg) {
		blcoeff = eb.coeff(x, bdeg);
		if (bdeg == 0)
			eb = _ex0;
		else
			eb -= blcoeff * power(x, bdeg);
	} else
		blcoeff = _ex1;

	while (rdeg >= bdeg && !r.is_zero()) {
		ex rlcoeff = r.coeff(x, rdeg);
		ex term = (power(x, rdeg - bdeg) * eb * rlcoeff).expand();
		if (rdeg == 0)
			r = _ex0;
		else
			r -= rlcoeff * power(x, rdeg);
		r = (blcoeff * r).expand() - term;
		rdeg = r.degree(x);
	}
	return r;
}


/** Exact polynomial division of a(X) by b(X) in Q[X].
 *  
 *  @param a  first multivariate polynomial (dividend)
 *  @param b  second multivariate polynomial (divisor)
 *  @param q  quotient (returned)
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return "true" when exact division succeeds (quotient returned in q),
 *          "false" otherwise (q left untouched) */
bool divide(const ex &a, const ex &b, ex &q, bool check_args)
{
	if (b.is_zero())
		throw(std::overflow_error("divide: division by zero"));
	if (a.is_zero()) {
		q = _ex0;
		return true;
	}
	if (is_exactly_a<numeric>(b)) {
		q = a / b;
		return true;
	}
        if (is_exactly_a<numeric>(a))
		return false;
#if FAST_COMPARE
	if (a.is_equal(b)) {
		q = _ex1;
		return true;
	}
#endif
	if (check_args && (!a.info(info_flags::rational_polynomial) ||
	                   !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("divide: arguments must be polynomials over the rationals"));

	// Find first symbol
	ex x;
	if (!a.get_first_symbol(x) && !b.get_first_symbol(x)) {
                std::ostringstream os;
                os << "invalid expression in divide(): " << a << " / " << b;
		throw(std::invalid_argument(os.str()));
        }

	// Try to avoid expanding partially factored expressions.
	if (is_exactly_a<mul>(b)) {
	// Divide sequentially by each term
		ex rem_new, rem_old = a;
		for (size_t i=0; i < b.nops(); i++) {
			if (! divide(rem_old, b.op(i), rem_new, false))
				return false;
			rem_old = rem_new;
		}
		q = rem_new;
		return true;
	}
        if (is_exactly_a<power>(b)) {
		const ex& bb(b.op(0));
		int exp_b = ex_to<numeric>(b.op(1)).to_int();
		ex rem_new, rem_old = a;
		for (int i=exp_b; i>0; i--) {
			if (! divide(rem_old, bb, rem_new, false))
				return false;
			rem_old = rem_new;
		}
		q = rem_new;
		return true;
	} 
	
	if (is_exactly_a<mul>(a)) {
		// Divide sequentially each term. If some term in a is divisible 
		// by b we are done... and if not, we can't really say anything.
		size_t i;
		ex rem_i;
		bool divisible_p = false;
		for (i=0; i < a.nops(); ++i) {
			if (divide(a.op(i), b, rem_i, false)) {
				divisible_p = true;
				break;
			}
		}
		if (divisible_p) {
			exvector resv;
			resv.reserve(a.nops());
			for (size_t j=0; j < a.nops(); j++) {
				if (j==i)
					resv.push_back(rem_i);
				else
					resv.push_back(a.op(j));
			}
			q = (new mul(resv))->setflag(status_flags::dynallocated);
			return true;
		}
	} else if (is_exactly_a<power>(a)) {
		// The base itself might be divisible by b, in that case we don't
		// need to expand a
		const ex& ab(a.op(0));
		int a_exp = ex_to<numeric>(a.op(1)).to_int();
		ex rem_i;
		if (divide(ab, b, rem_i, false)) {
			q = rem_i*power(ab, a_exp - 1);
			return true;
		}
// code below is commented-out because it leads to a significant slowdown
//		for (int i=2; i < a_exp; i++) {
//			if (divide(power(ab, i), b, rem_i, false)) {
//				q = rem_i*power(ab, a_exp - i);
//				return true;
//			}
//		} // ... so we *really* need to expand expression.
	}
	
	// Polynomial long division (recursive)
	ex r = a.expand();
	if (r.is_zero()) {
		q = _ex0;
		return true;
	}
	numeric bdeg = b.degree(x);
	numeric rdeg = r.degree(x);
	ex blcoeff = b.expand().coeff(x, bdeg);
	bool blcoeff_is_numeric = is_exactly_a<numeric>(blcoeff);
	exvector v; //v.reserve(std::max(rdeg - bdeg + 1, 0));
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (blcoeff_is_numeric)
			term = rcoeff / blcoeff;
		else
			if (!divide(rcoeff, blcoeff, term, false))
				return false;
		term *= power(x, rdeg - bdeg);
		v.push_back(term);
		r -= (term * b).expand();
		if (r.is_zero()) {
			q = (new add(v))->setflag(status_flags::dynallocated);
			return true;
		}
		rdeg = r.degree(x);
	}
	return false;
}

// Distribute a factor over all sum terms
static ex dist(const ex& f, const ex& p)
{
        if (not is_exactly_a<add>(p))
                return mul(f, p);
        ex s = _ex0;
        for (size_t i=0; i<p.nops(); ++i)
                s += f*p.op(i);
        return s;
}

/** Compute partial fraction decomposition of expression
 *
 *  @param a expression
 *  @param x variable to factor in
 *  @return decomposed rational function */
ex parfrac(const ex & a, const ex & x)
{
        if (is_exactly_a<add>(a)) {
                ex s = _ex0;
                for (size_t i=0; i<a.nops(); ++i)
                        s += parfrac(a.op(i), x);
                return s;
        }
        if (not is_exactly_a<mul>(a))
                return a;
	// Find numerator and denominator
	ex nd = numer_denom(a);
	ex numer = nd.op(0), denom = nd.op(1);
        // Expand sum powers
        if (is_exactly_a<add>(numer))
                return parfrac(dist(_ex1/denom, numer), x);
        if (is_exactly_a<power>(numer)
            and is_exactly_a<add>(numer.op(0))
            and numer.op(1).info(info_flags::posint))
                return parfrac(dist(_ex1/denom, numer.expand()), x);

        if (is_exactly_a<mul>(numer))
                for (size_t i=0; i<numer.nops(); ++i) {
                        const ex& t = numer.op(i);
                        if (is_exactly_a<power>(t)
                            and is_exactly_a<add>(t.op(0))
                            and t.op(1).info(info_flags::posint))
                                return parfrac(dist(_ex1/denom, numer.expand()), x);
                }
        std::pair<ex,ex> qr;
        try {
        	// Convert N(x)/D(x) -> Q(x) + R(x)/D(x), so degree(R) < degree(D)
                qr = quo_rem(numer, denom, x, true);
        }
        catch (std::logic_error) {
                return a;
        }
	// Factorize denominator and compute cofactors
        ex facmul;
        bool r = factor(denom, facmul);
        if (not r)
                facmul = denom;
//                return qr.first + qr.second/denom;
        ex oc = _ex1;
	exvector factor, cofac;
        if (is_exactly_a<mul>(facmul)) {
                size_t fsize = facmul.nops();
                for (size_t i=0; i<fsize; i++) {
                        const ex& e = facmul.op(i);
                        if (is_exactly_a<numeric>(e)) {
                                oc = e;
                                continue;
                        }
                        if (is_exactly_a<power>(e)) {
                                if (not has(e.op(0), x))
                                        continue;
                                const ex& ee = e.op(1);
                                if (not is_exactly_a<numeric>(ee)
                                    or not ee.is_integer()
                                    or e.op(0).is_equal(x))
                                        return a;
                                size_t expo = ex_to<numeric>(ee).to_int();
                                for (size_t j=1; j<=expo; ++j) {
                                        ex eee = power(e.op(0), numeric(static_cast<int>(j)));
                                        factor.push_back(eee);
                                        cofac.push_back((facmul/eee).expand());
                                }
                        }
                        else {
                                factor.push_back(e);
                                cofac.push_back((facmul/e).expand());
                        }
                }
        }
        else if (is_exactly_a<power>(facmul)
                 and is_exactly_a<numeric>(facmul.op(1))) {
                if (facmul.op(0).is_equal(x))
                        return a;
                size_t expo = ex_to<numeric>(facmul.op(1)).to_int();
                for (size_t j=1; j<=expo; ++j) {
                        ex ee = power(facmul.op(0), numeric(static_cast<int>(j)));
                        factor.push_back(ee);
                        cofac.push_back((facmul/ee).expand());
                }
        }
        else {
                return qr.first + qr.second/denom;
        }
	size_t num_factors = factor.size();
        //std::cerr << "factors  : " << exprseq(factor) << std::endl;
        //std::cerr << "cofactors: " << exprseq(cofac) << std::endl;

	// Construct coefficient matrix for decomposition
	int max_denom_deg = ex_to<numeric>(denom.degree(x)).to_int();
	matrix sys(max_denom_deg + 1, num_factors);
	matrix rhs(max_denom_deg + 1, 1);
	for (int i=0; i<=max_denom_deg; i++) {
		for (size_t j=0; j<num_factors; j++)
			sys(i, j) = cofac[j].coeff(x, i);
		rhs(i, 0) = qr.second.coeff(x, i);
	}
//clog << "coeffs: " << sys << endl;
//clog << "rhs   : " << rhs << endl;

	// Solve resulting linear system
	matrix vars(num_factors, 1);
	for (size_t i=0; i<num_factors; i++)
		vars(i, 0) = symbol();
	matrix sol = sys.solve(vars, rhs);

	// Sum up decomposed fractions
	ex sum = _ex0;
	for (size_t i=0; i<num_factors; i++)
		sum += sol(i, 0) / oc / factor[i];

	return qr.first + sum;
}


} // namespace GiNaC 
