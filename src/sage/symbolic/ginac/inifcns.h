/** @file inifcns.h
 *
 *  Interface to GiNaC's initially known functions. */

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

#ifndef __GINAC_INIFCNS_H__
#define __GINAC_INIFCNS_H__

#include "numeric.h"
#include "function.h"
#include "ex.h"

namespace GiNaC {

/** Complex conjugate. */
DECLARE_FUNCTION_1P(conjugate_function)

/** Real part. */
DECLARE_FUNCTION_1P(real_part_function)

/** Imaginary part. */
DECLARE_FUNCTION_1P(imag_part_function)
	
/** Absolute value. */
DECLARE_FUNCTION_1P(abs)

/** Step function. */
DECLARE_FUNCTION_1P(unit_step)
	
/** Heaviside function. */
DECLARE_FUNCTION_1P(heaviside)
	
/** Complex sign. */
DECLARE_FUNCTION_1P(csgn)

/** Eta function: log(a*b) == log(a) + log(b) + eta(a, b). */
DECLARE_FUNCTION_2P(eta)

/** Sine. */
DECLARE_FUNCTION_1P(sin)

/** Cosine. */
DECLARE_FUNCTION_1P(cos)

/** Tangent. */
DECLARE_FUNCTION_1P(tan)

/** Secant. */
DECLARE_FUNCTION_1P(sec)

/** Cosecant. */
DECLARE_FUNCTION_1P(csc)

/** Cotangent. */
DECLARE_FUNCTION_1P(cot)

/** Exponential function. */
DECLARE_FUNCTION_1P(exp)

/** Natural logarithm. */
DECLARE_FUNCTION_1P(log)

/** General logarithm. */
DECLARE_FUNCTION_2P(logb)

/** Inverse sine (arc sine). */
DECLARE_FUNCTION_1P(asin)

/** Inverse cosine (arc cosine). */
DECLARE_FUNCTION_1P(acos)

/** Inverse tangent (arc tangent). */
DECLARE_FUNCTION_1P(atan)

/** Inverse cotangent (arc cotangent). */
DECLARE_FUNCTION_1P(acot)

/** Inverse secant (arc secant). */
DECLARE_FUNCTION_1P(asec)

/** Inverse cosecant (arc cosecant). */
DECLARE_FUNCTION_1P(acsc)

/** Inverse tangent with two arguments. */
DECLARE_FUNCTION_2P(atan2)

/** Hyperbolic Sine. */
DECLARE_FUNCTION_1P(sinh)

/** Hyperbolic Cosine. */
DECLARE_FUNCTION_1P(cosh)

/** Hyperbolic Tangent. */
DECLARE_FUNCTION_1P(tanh)

/** Hyperbolic Cotangent. */
DECLARE_FUNCTION_1P(coth)

/** Hyperbolic Secant. */
DECLARE_FUNCTION_1P(sech)

/** Hyperbolic Cosecant. */
DECLARE_FUNCTION_1P(csch)

/** Inverse hyperbolic Sine (area hyperbolic sine). */
DECLARE_FUNCTION_1P(asinh)

/** Inverse hyperbolic Cosine (area hyperbolic cosine). */
DECLARE_FUNCTION_1P(acosh)

/** Inverse hyperbolic Tangent (area hyperbolic tangent). */
DECLARE_FUNCTION_1P(atanh)

/** Inverse hyperbolic Cotangent (area hyperbolic cotangent). */
DECLARE_FUNCTION_1P(acoth)


/** Inverse hyperbolic Cosecant (area hyperbolic cosecant). */
DECLARE_FUNCTION_1P(acsch)

/** Inverse hyperbolic Secant (area hyperbolic secant). */
DECLARE_FUNCTION_1P(asech)


/** Dilogarithm. */
DECLARE_FUNCTION_1P(Li2)

/** Derivatives of Riemann's Zeta-function. */
DECLARE_FUNCTION_2P(zetaderiv)

// overloading at work: we cannot use the macros here
/** Multiple zeta value including Riemann's zeta-function. */
class zeta1_SERIAL { public: static unsigned serial; };
template<typename T1>
inline function zeta(const T1& p1) {
	return function(zeta1_SERIAL::serial, ex(p1));
}
/** Alternating Euler sum or colored MZV. */
class zeta2_SERIAL { public: static unsigned serial; };
template<typename T1, typename T2>
inline function zeta(const T1& p1, const T2& p2) {
	return function(zeta2_SERIAL::serial, ex(p1), ex(p2));
}
class zeta_SERIAL;
template<> inline bool is_the_function<zeta_SERIAL>(const ex& x)
{
	return is_the_function<zeta1_SERIAL>(x) || is_the_function<zeta2_SERIAL>(x);
}

class stieltjes1_SERIAL { public: static unsigned serial; };
template<typename T1>
inline function stieltjes(const T1& p1) {
	return function(stieltjes1_SERIAL::serial, ex(p1));
}

// overloading at work: we cannot use the macros here
/** Generalized multiple polylogarithm. */
class G2_SERIAL { public: static unsigned serial; };
template<typename T1, typename T2>
inline function G(const T1& x, const T2& y) {
	return function(G2_SERIAL::serial, ex(x), ex(y));
}
/** Generalized multiple polylogarithm with explicit imaginary parts. */
class G3_SERIAL { public: static unsigned serial; };
template<typename T1, typename T2, typename T3>
inline function G(const T1& x, const T2& s, const T3& y) {
	return function(G3_SERIAL::serial, ex(x), ex(s), ex(y));
}
class G_SERIAL;
template<> inline bool is_the_function<G_SERIAL>(const ex& x)
{
	return is_the_function<G2_SERIAL>(x) || is_the_function<G3_SERIAL>(x);
}

/** Polylogarithm and multiple polylogarithm. */
DECLARE_FUNCTION_2P(Li)

/** Nielsen's generalized polylogarithm. */
DECLARE_FUNCTION_3P(S)

/** Harmonic polylogarithm. */
DECLARE_FUNCTION_2P(H)

/** Gamma-function. */
DECLARE_FUNCTION_1P(lgamma)
DECLARE_FUNCTION_1P(gamma)

/** Beta-function. */
DECLARE_FUNCTION_2P(beta)

// overloading at work: we cannot use the macros here
/** Psi-function (aka digamma-function). */
class psi1_SERIAL { public: static unsigned serial; };
template<typename T1>
inline function psi(const T1 & p1) {
	return function(psi1_SERIAL::serial, ex(p1));
}
/** Derivatives of Psi-function (aka polygamma-functions). */
class psi2_SERIAL { public: static unsigned serial; };
template<typename T1, typename T2>
inline function psi(const T1 & p1, const T2 & p2) {
	return function(psi2_SERIAL::serial, ex(p1), ex(p2));
}
class psi_SERIAL;
template<> inline bool is_the_function<psi_SERIAL>(const ex & x)
{
	return is_the_function<psi1_SERIAL>(x) || is_the_function<psi2_SERIAL>(x);
}
	
/** Factorial function. */
DECLARE_FUNCTION_1P(factorial)

/** Binomial function. */
DECLARE_FUNCTION_2P(binomial)

/** Rising factorial function. */
DECLARE_FUNCTION_2P(rising_factorial)

/** Falling factorial function. */
DECLARE_FUNCTION_2P(falling_factorial)

/** Chebyshev T polynomial. */
DECLARE_FUNCTION_2P(chebyshev_T)

/** Chebyshev U polynomial. */
DECLARE_FUNCTION_2P(chebyshev_U)

/** Legendre P polynomial. */
DECLARE_FUNCTION_2P(legendre_P)

/** Hermite polynomial. */
DECLARE_FUNCTION_2P(hermite)

/** Gegenbauer (ultraspherical) polynomial. */
DECLARE_FUNCTION_3P(gegenbauer)

/** Appell F1 function */
DECLARE_FUNCTION_6P(appell_F1)

/** Order term function (for truncated power series). */
DECLARE_FUNCTION_1P(Order)

/** Formal piecewise function */
DECLARE_FUNCTION_1P(cases)

/** Formal set-of-all function */
DECLARE_FUNCTION_2P(set_of_all)

ex lsolve(const ex &eqns, const ex &symbols, unsigned options = solve_algo::automatic);

/** Find a real root of real-valued function f(x) numerically within a given
 *  interval. The function must change sign across interval. Uses Newton-
 *  Raphson method combined with bisection in order to guarantee convergence.
 *
 *  @param f  Function f(x)
 *  @param x  Symbol f(x)
 *  @param x1  lower interval limit
 *  @param x2  upper interval limit
 *  @exception runtime_error (if interval is invalid). */
const numeric fsolve(const ex& f, const symbol& x, const numeric& x1, const numeric& x2, PyObject* parent);

/** Check whether a function is the Order (O(n)) function. */
inline bool is_order_function(const ex & e)
{
	return is_ex_the_function(e, Order);
}

/** Converts a given list containing parameters for H in Remiddi/Vermaseren notation into
 *  the corresponding GiNaC functions.
 */
ex convert_H_to_Li(const ex& parameterlst, const ex& arg);

} // namespace GiNaC

#endif // ndef __GINAC_INIFCNS_H__
