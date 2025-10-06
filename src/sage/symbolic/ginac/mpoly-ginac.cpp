/** @file mpoly-ginac.cpp
 *
 *  This file implements several functions that work multivariate polynomials.
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

#include "basic.h"
#include "ex.h"
#include "mpoly.h"
#include "upoly.h"
#include "add.h"
#include "constant.h"
#include "expairseq.h"
#include "inifcns.h"
#include "lst.h"
#include "mul.h"
#include "numeric.h"
#include "power.h"
#include "relational.h"
#include "operators.h"
#include "matrix.h"
#include "pseries.h"
#include "symbol.h"
#include "utils.h"
#include "upoly.h"
#include "normal.h"

#include <algorithm>
#include <map>

namespace GiNaC {

// If comparing expressions (ex::compare()) is fast, you can set this to 1.
// Some routines like quo(), rem() and gcd() will then return a quick answer
// when they are called with two identical arguments.
#define FAST_COMPARE 1

// Set this if you want divide_in_z() to use remembering
#define USE_REMEMBER 0

// Set this if you want divide_in_z() to use trial division followed by
// polynomial interpolation (always slower except for completely dense
// polynomials)
#define USE_TRIAL_DIVISION 0

// Set this to enable some statistical output for the GCD routines
#define STATISTICS 0


#if STATISTICS
// Statistics variables
static int gcd_called = 0;
static int sr_gcd_called = 0;

// Print statistics at end of program
static struct _stat_print {
	_stat_print() {}
	~_stat_print() {
		std::cout << "gcd() called " << gcd_called << " times\n";
		std::cout << "sr_gcd() called " << sr_gcd_called << " times\n";
	}
} stat_print;
#endif

#if 0
/*
 *  Statistical information about symbols in polynomials
 */

/** This structure holds information about the highest and lowest degrees
 *  in which a symbol appears in two multivariate polynomials "a" and "b".
 *  A vector of these structures with information about all symbols in
 *  two polynomials can be created with the function get_symbol_stats().
 *
 *  @see get_symbol_stats */
struct sym_desc {
	/** Reference to symbol */
	ex sym;

	/** Highest degree of symbol in polynomial "a" */
	int deg_a;

	/** Highest degree of symbol in polynomial "b" */
	int deg_b;

	/** Lowest degree of symbol in polynomial "a" */
	int ldeg_a;

	/** Lowest degree of symbol in polynomial "b" */
	int ldeg_b;

	/** Maximum of deg_a and deg_b (Used for sorting) */
	int max_deg;

	/** Maximum number of terms of leading coefficient of symbol in both polynomials */
	size_t max_lcnops;

	/** Comparison operator for sorting */
	bool operator<(const sym_desc &x) const
	{
		if (max_deg == x.max_deg)
			return max_lcnops < x.max_lcnops;
		else
			return max_deg < x.max_deg;
	}
};

// Vector of sym_desc structures
typedef std::vector<sym_desc> sym_desc_vec;

// Add symbol the sym_desc_vec (used internally by get_symbol_stats())
static void add_symbol(const ex &s, sym_desc_vec &v)
{
	for (const auto& elem : v)
		if (elem.sym.is_equal(s))  // If it's already in there, don't add it a second time
			return;

	sym_desc d;
	d.sym = s;
	v.push_back(d);
}

// Collect all symbols of an expression (used internally by get_symbol_stats())
static void collect_symbols(const ex &e, sym_desc_vec &v)
{
	if (is_exactly_a<symbol>(e)) {
		add_symbol(e, v);
	} else if (is_exactly_a<add>(e) || is_exactly_a<mul>(e)) {
		for (size_t i=0; i<e.nops(); i++)
			collect_symbols(e.sorted_op(i), v);
	} else if (is_exactly_a<power>(e)) {
		collect_symbols(e.op(0), v);
	}
}

/** Collect statistical information about symbols in polynomials.
 *  This function fills in a vector of "sym_desc" structs which contain
 *  information about the highest and lowest degrees of all symbols that
 *  appear in two polynomials. The vector is then sorted by minimum
 *  degree (lowest to highest). The information gathered by this
 *  function is used by the GCD routines to identify trivial factors
 *  and to determine which variable to choose as the main variable
 *  for GCD computation.
 *
 *  @param a  first multivariate polynomial
 *  @param b  second multivariate polynomial
 *  @param v  vector of sym_desc structs (filled in) */
static void get_symbol_stats(const ex &a, const ex &b, sym_desc_vec &v)
{
	collect_symbols(a.eval(), v);   // eval() to expand assigned symbols
	collect_symbols(b.eval(), v);
	auto it = v.begin(), itend = v.end();
	while (it != itend) {
		int deg_a = a.degree(it->sym);
		int deg_b = b.degree(it->sym);
		it->deg_a = deg_a;
		it->deg_b = deg_b;
		it->max_deg = std::max(deg_a, deg_b);
		it->max_lcnops = std::max(a.lcoeff(it->sym).nops(), b.lcoeff(it->sym).nops());
		it->ldeg_a = a.ldegree(it->sym);
		it->ldeg_b = b.ldegree(it->sym);
		++it;
	}
	std::sort(v.begin(), v.end());

#if 0
	std::clog << "Symbols:\n";
	it = v.begin(); itend = v.end();
	while (it != itend) {
		std::clog << " " << it->sym << ": deg_a=" << it->deg_a << ", deg_b=" << it->deg_b << ", ldeg_a=" << it->ldeg_a << ", ldeg_b=" << it->ldeg_b << ", max_deg=" << it->max_deg << ", max_lcnops=" << it->max_lcnops << endl;
		std::clog << "  lcoeff_a=" << a.lcoeff(it->sym) << ", lcoeff_b=" << b.lcoeff(it->sym) << endl;
		++it;
	}
#endif
}
#endif

#if 0
/** Exact polynomial division of a(X) by b(X) in Z[X].
 *  This functions works like divide() but the input and output polynomials are
 *  in Z[X] instead of Q[X] (i.e. they have integer coefficients). Unlike
 *  divide(), it doesn't check whether the input polynomials really are integer
 *  polynomials, so be careful of what you pass in. Also, you have to run
 *  get_symbol_stats() over the input polynomials before calling this function
 *  and pass an iterator to the first element of the sym_desc vector. This
 *  function was used internally by the heur_gcd().
 *  
 *  @param a  first multivariate polynomial (dividend)
 *  @param b  second multivariate polynomial (divisor)
 *  @param q  quotient (returned)
 *  @param var  iterator to first element of vector of sym_desc structs
 *  @return "true" when exact division succeeds (the quotient is returned in
 *          q), "false" otherwise.
 *  @see get_symbol_stats */
static bool divide_in_z(const ex &a, const ex &b, ex &q, sym_desc_vec::const_iterator var)
{
	q = _ex0;
	if (b.is_zero())
		throw(std::overflow_error("divide_in_z: division by zero"));
	if (b.is_equal(_ex1)) {
		q = a;
		return true;
	}
	if (is_exactly_a<numeric>(a)) {
		if (is_exactly_a<numeric>(b)) {
			q = a / b;
			return q.info(info_flags::integer);
		} else
			return false;
	}
#if FAST_COMPARE
	if (a.is_equal(b)) {
		q = _ex1;
		return true;
	}
#endif

#if USE_REMEMBER
	// Remembering
	static ex2_exbool_remember dr_remember;
	ex2_exbool_remember::const_iterator remembered = dr_remember.find(ex2(a, b));
	if (remembered != dr_remember.end()) {
		q = remembered->second.first;
		return remembered->second.second;
	}
#endif

	if (is_exactly_a<power>(b)) {
		const ex& bb(b.op(0));
		ex qbar = a;
		int exp_b = ex_to<numeric>(b.op(1)).to_int();
		for (int i=exp_b; i>0; i--) {
			if (!divide_in_z(qbar, bb, q, var))
				return false;
			qbar = q;
		}
		return true;
	}

	if (is_exactly_a<mul>(b)) {
		ex qbar = a;
		for (const auto& elem : b) {
			sym_desc_vec sym_stats;
			get_symbol_stats(a, elem, sym_stats);
			if (!divide_in_z(qbar, elem, q, sym_stats.begin()))
				return false;

			qbar = q;
		}
		return true;
	}

	// Main symbol
	const ex &x = var->sym;

	// Compare degrees
	int adeg = a.degree(x), bdeg = b.degree(x);
	if (bdeg > adeg)
		return false;

#if USE_TRIAL_DIVISION

	// Trial division with polynomial interpolation
	int i, k;

	// Compute values at evaluation points 0..adeg
	vector<numeric> alpha; alpha.reserve(adeg + 1);
	exvector u; u.reserve(adeg + 1);
	numeric point = *_num0_p;
	ex c;
	for (i=0; i<=adeg; i++) {
		ex bs = b.subs(x == point, subs_options::no_pattern);
		while (bs.is_zero()) {
			point += *_num1_p;
			bs = b.subs(x == point, subs_options::no_pattern);
		}
		if (!divide_in_z(a.subs(x == point, subs_options::no_pattern), bs, c, var+1))
			return false;
		alpha.push_back(point);
		u.push_back(c);
		point += *_num1_p;
	}

	// Compute inverses
	vector<numeric> rcp; rcp.reserve(adeg + 1);
	rcp.push_back(*_num0_p);
	for (k=1; k<=adeg; k++) {
		numeric product = alpha[k] - alpha[0];
		for (i=1; i<k; i++)
			product *= alpha[k] - alpha[i];
		rcp.push_back(product.inverse());
	}

	// Compute Newton coefficients
	exvector v; v.reserve(adeg + 1);
	v.push_back(u[0]);
	for (k=1; k<=adeg; k++) {
		ex temp = v[k - 1];
		for (i=k-2; i>=0; i--)
			temp = temp * (alpha[k] - alpha[i]) + v[i];
		v.push_back((u[k] - temp) * rcp[k]);
	}

	// Convert from Newton form to standard form
	c = v[adeg];
	for (k=adeg-1; k>=0; k--)
		c = c * (x - alpha[k]) + v[k];

	if (c.degree(x) == (adeg - bdeg)) {
		q = c.expand();
		return true;
	} else
		return false;

#else

	// Polynomial long division (recursive)
	ex r = a.expand();
	if (r.is_zero())
		return true;
	int rdeg = adeg;
	ex eb = b.expand();
	ex blcoeff = eb.coeff(x, bdeg);
	exvector v; v.reserve(std::max(rdeg - bdeg + 1, 0));
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (!divide_in_z(rcoeff, blcoeff, term, var+1))
			break;
		term = (term * power(x, rdeg - bdeg)).expand();
		v.push_back(term);
		r -= (term * eb).expand();
		if (r.is_zero()) {
			q = (new add(v))->setflag(status_flags::dynallocated);
#if USE_REMEMBER
			dr_remember[ex2(a, b)] = exbool(q, true);
#endif
			return true;
		}
		rdeg = r.degree(x);
	}
#if USE_REMEMBER
	dr_remember[ex2(a, b)] = exbool(q, false);
#endif
	return false;

#endif
}

/** xi-adic polynomial interpolation */
static ex interpolate(const ex &gamma, const numeric &xi, const ex &x, int degree_hint = 1)
{
	exvector g; g.reserve(degree_hint);
	ex e = gamma;
	numeric rxi = xi.inverse();
	for (int i=0; !e.is_zero(); i++) {
		ex gi = e.smod(xi);
		g.push_back(gi * power(x, i));
		e = (e - gi) * rxi;
	}
	return (new add(g))->setflag(status_flags::dynallocated);
}
#endif


/** Compute LCM (Least Common Multiple) of multivariate polynomials in Z[X].
 *
 *  @param a  first multivariate polynomial
 *  @param b  second multivariate polynomial
 *  @param check_args  check whether a and b are polynomials with rational
 *         coefficients (defaults to "true")
 *  @return the LCM as a new expression */
ex lcm(const ex &a, const ex &b, bool check_args)
{
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b))
		return lcm(ex_to<numeric>(a), ex_to<numeric>(b));
	if (check_args && (!a.info(info_flags::rational_polynomial) || !b.info(info_flags::rational_polynomial)))
		throw(std::invalid_argument("lcm: arguments must be polynomials over the rationals"));
	
	ex ca, cb;
	ex g = gcdpoly(a, b, &ca, &cb, false);
	return ca * cb * g;
}

#if 0
/*
 *  Square-free factorization
 */

/** Compute square-free factorization of multivariate polynomial a(x) using
 *  Yun's algorithm.  Used internally by sqrfree().
 *
 *  @param a  multivariate polynomial over Z[X], treated here as univariate
 *            polynomial in x (needs not be expanded).
 *  @param x  variable to factor in
 *  @return   vector of factors sorted in ascending degree */
static exvector sqrfree_yun(const ex &a, const symbol &x)
{
	exvector res;
	ex w = a;
	ex z = w.diff(x);
	ex g = gcdpoly(w, z);
	if (g.is_zero()) {
		return res;
	}
	if (g.is_equal(_ex1)) {
		res.push_back(a);
		return res;
	}
	ex y;
	do {
		w = quo(w, g, x);
		if (w.is_zero()) {
			return res;
		}
		y = quo(z, g, x);
		z = y - w.diff(x);
		g = gcdpoly(w, z);
		res.push_back(g);
	} while (!z.is_zero());
	return res;
}



/** Compute a square-free factorization of a multivariate polynomial in Q[X].
 *
 *  @param a  multivariate polynomial over Q[X] (needs not be expanded)
 *  @param l  lst of variables to factor in, may be left empty for autodetection
 *  @return   a square-free factorization of \p a.
 *
 * \note
 * A polynomial \f$p(X) \in C[X]\f$ is said <EM>square-free</EM>
 * if, whenever any two polynomials \f$q(X)\f$ and \f$r(X)\f$
 * are such that
 * \f[
 *     p(X) = q(X)^2 r(X),
 * \f]
 * we have \f$q(X) \in C\f$.
 * This means that \f$p(X)\f$ has no repeated factors, apart
 * eventually from constants.
 * Given a polynomial \f$p(X) \in C[X]\f$, we say that the
 * decomposition
 * \f[
 *   p(X) = b \cdot p_1(X)^{a_1} \cdot p_2(X)^{a_2} \cdots p_r(X)^{a_r}
 * \f]
 * is a <EM>square-free factorization</EM> of \f$p(X)\f$ if the
 * following conditions hold:
 * -#  \f$b \in C\f$ and \f$b \neq 0\f$;
 * -#  \f$a_i\f$ is a positive integer for \f$i = 1, \ldots, r\f$;
 * -#  the degree of the polynomial \f$p_i\f$ is strictly positive
 *     for \f$i = 1, \ldots, r\f$;
 * -#  the polynomial \f$\Pi_{i=1}^r p_i(X)\f$ is square-free.
 *
 * Square-free factorizations need not be unique.  For example, if
 * \f$a_i\f$ is even, we could change the polynomial \f$p_i(X)\f$
 * into \f$-p_i(X)\f$.
 * Observe also that the factors \f$p_i(X)\f$ need not be irreducible
 * polynomials.
 */
ex sqrfree(const ex &a, const lst &l)
{
	if (is_exactly_a<numeric>(a) ||     // algorithm does not trap a==0
	    is_exactly_a<symbol>(a))        // shortcut
		return a;

	// If no lst of variables to factorize in was specified we have to
	// invent one now.  Maybe one can optimize here by reversing the order
	// or so, I don't know.
	lst args;
	if (l.nops()==0) {
		sym_desc_vec sdv;
		get_symbol_stats(a, _ex0, sdv);
		for (const auto& elem : sdv)
			args.append(elem.sym);
	} else {
		args = l;
	}

	// Find the symbol to factor in at this stage
	if (!is_exactly_a<symbol>(args.op(0)))
		throw (std::runtime_error("sqrfree(): invalid factorization variable"));
	const symbol &x = ex_to<symbol>(args.op(0));

	// convert the argument from something in Q[X] to something in Z[X]
	const numeric lcm = lcm_of_coefficients_denominators(a);
	const ex tmp = multiply_lcm(a,lcm);

	// find the factors
	exvector factors = sqrfree_yun(tmp, x);

	// construct the next list of symbols with the first element popped
	lst newargs = args;
	newargs.remove_first();

	// recurse down the factors in remaining variables
	if (newargs.nops()>0) {
		auto i = factors.begin();
		while (i != factors.end()) {
			*i = sqrfree(*i, newargs);
			++i;
		}
	}

	// Done with recursion, now construct the final result
	ex result = _ex1;
        {
        int p = 1;
	for (const auto& elem : factors)
		result *= power(elem, p++);
        }

	// Yun's algorithm does not account for constant factors.  (For univariate
	// polynomials it works only in the monic case.)  We can correct this by
	// inserting what has been lost back into the result.  For completeness
	// we'll also have to recurse down that factor in the remaining variables.
	if (newargs.nops()>0)
		result *= sqrfree(quo(tmp, result, x), newargs);
	else
		result *= quo(tmp, result, x);

	// Put in the reational overall factor again and return
	return result *	lcm.inverse();
}
#endif


} // namespace GiNaC

