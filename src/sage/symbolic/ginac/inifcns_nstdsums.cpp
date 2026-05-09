/** @file inifcns_nstdsums.cpp
 *
 *  Implementation of some special functions that have a representation as nested sums.
 *
 *  The functions are:
 *    classical polylogarithm              Li(n,x)
 *    multiple polylogarithm               Li(lst(m_1,...,m_k),lst(x_1,...,x_k))
 *                                         G(lst(a_1,...,a_k),y) or G(lst(a_1,...,a_k),lst(s_1,...,s_k),y)
 *    Nielsen's generalized polylogarithm  S(n,p,x)
 *    harmonic polylogarithm               H(m,x) or H(lst(m_1,...,m_k),x)
 *    multiple zeta value                  zeta(m) or zeta(lst(m_1,...,m_k))
 *    alternating Euler sum                zeta(m,s) or zeta(lst(m_1,...,m_k),lst(s_1,...,s_k))
 *
 *  Some remarks:
 *
 *    - All formulae used can be looked up in the following publications:
 *      [Kol] Nielsen's Generalized Polylogarithms, K.S.Kolbig, SIAM J.Math.Anal. 17 (1986), pp. 1232-1258.
 *      [Cra] Fast Evaluation of Multiple Zeta Sums, R.E.Crandall, Math.Comp. 67 (1998), pp. 1163-1172.
 *      [ReV] Harmonic Polylogarithms, E.Remiddi, J.A.M.Vermaseren, Int.J.Mod.Phys. A15 (2000), pp. 725-754
 *      [BBB] Special Values of Multiple Polylogarithms, J.Borwein, D.Bradley, D.Broadhurst, P.Lisonek, Trans.Amer.Math.Soc. 353/3 (2001), pp. 907-941
 *      [VSW] Numerical evaluation of multiple polylogarithms, J.Vollinga, S.Weinzierl, hep-ph/0410259
 *
 *    - The order of parameters and arguments of Li and zeta is defined according to the nested sums
 *      representation. The parameters for H are understood as in [ReV]. They can be in expanded --- only
 *      0, 1 and -1 --- or in compactified --- a string with zeros in front of 1 or -1 is written as a single
 *      number --- notation.
 *
 *    - All functions can be numerically evaluated with arguments in the whole complex plane. The parameters
 *      for Li, zeta and S must be positive integers. If you want to have an alternating Euler sum, you have
 *      to give the signs of the parameters as a second argument s to zeta(m,s) containing 1 and -1.
 *
 *    - The calculation of classical polylogarithms is speeded up by using Bernoulli numbers and 
 *      look-up tables. S uses look-up tables as well. The zeta function applies the algorithms in
 *      [Cra] and [BBB] for speed up. Multiple polylogarithms use Hoelder convolution [BBB].
 *
 *    - The functions have no means to do a series expansion into nested sums. To do this, you have to convert
 *      these functions into the appropriate objects from the nestedsums library, do the expansion and convert
 *      the result back.
 *
 *    - Numerical testing of this implementation has been performed by doing a comparison of results
 *      between this software and the commercial M.......... 4.1. Multiple zeta values have been checked
 *      by means of evaluations into simple zeta values. Harmonic polylogarithms have been checked by
 *      comparison to S(n,p,x) for corresponding parameter combinations and by continuity checks
 *      around |x|=1 along with comparisons to corresponding zeta functions. Multiple polylogarithms were
 *      checked against H and zeta and by means of shuffle and quasi-shuffle relations.
 *
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

#include "add.h"
#include "constant.h"
#include "infinity.h"
#include "lst.h"
#include "mul.h"
#include "numeric.h"
#include "operators.h"
#include "power.h"
#include "pseries.h"
#include "relational.h"
#include "symbol.h"
#include "utils.h"
#include "wildcard.h"

#include <sstream>
#include <stdexcept>
#include <vector>
//#include <cln/cln.h>


namespace GiNaC {


//////////////////////////////////////////////////////////////////////
//
// Classical polylogarithm  Li(n,x)
//
// helper functions
//
//////////////////////////////////////////////////////////////////////


// anonymous namespace for helper functions
namespace {

/*
 * 
// lookup table for factors built from Bernoulli numbers
// see fill_Xn()
//std::vector<std::vector<numeric> > Xn;
// initial size of Xn that should suffice for 32bit machines (must be even)
// const int xninitsizestep = 26;
// int xninitsize = xninitsizestep;
// int xnsize = 0;

  // Waiting to be implemented
 // void binomial(Number_T x, Number_T y) {throw dunno();}
 // void factorial (Number_T x) {throw dunno();}

// This function calculates the X_n. The X_n are needed for speed up of classical polylogarithms.
// With these numbers the polylogs can be calculated as follows:
//   Li_p (x)  =  \sum_{n=0}^\infty X_{p-2}(n) u^{n+1}/(n+1)! with  u = -log(1-x)
//   X_0(n) = B_n (Bernoulli numbers)
//   X_p(n) = \sum_{k=0}^n binomial(n,k) B_{n-k} / (k+1) * X_{p-1}(k)
// The calculation of Xn depends on X0 and X{n-1}.
// X_0 is special, it holds only the non-zero Bernoulli numbers with index 2 or greater.
// This results in a slightly more complicated algorithm for the X_n.
// The first index in Xn corresponds to the index of the polylog minus 2.
// The second index in Xn corresponds to the index from the actual sum.
void fill_Xn(int n)
{
	if (n>1) {
		// calculate X_2 and higher (corresponding to Li_4 and higher)
		std::vector<Number_T> buf(xninitsize);
		std::vector<Number_T>::iterator it = buf.begin();
		Number_T result;
		*it = -(pow(Integer_T(2),n+1) - 1) / pow(Integer_T(2),n+1); // i == 1
		it++;
		for (int i=2; i<=xninitsize; i++) {
			if (i&1) {
				result = 0; // k == 0
			} else {
				result = Xn[0][i/2-1]; // k == 0
			}
			for (int k=1; k<i-1; k++) {
				if ( !(((i-k) & 1) && ((i-k) > 1)) ) {
					result = result + binomial(i,k) * Xn[0][(i-k)/2-1] * Xn[n-1][k-1] / (k+1);
				}
			}
			result = result - binomial(i,i-1) * Xn[n-1][i-2] / 2 / i; // k == i-1
			result = result + Xn[n-1][i-1] / (i+1); // k == i
			
			*it = result;
			it++;
		}
		Xn.push_back(buf);
	} else if (n==1) {
		// special case to handle the X_0 correct
		std::vector<Number_T> buf(xninitsize);
		std::vector<Number_T>::iterator it = buf.begin();
		Number_T result;
		*it = Integer_T(-3)/Integer_T(4); // i == 1
		it++;
		*it = Integer_T(17)/Integer_T(36); // i == 2
		it++;
		for (int i=3; i<=xninitsize; i++) {
			if (i & 1) {
				result = -Xn[0][(i-3)/2]/2;
				*it = (binomial(i,1)/Integer_T(2) + binomial(i,i-1)/Integer_T(i))*result;
				it++;
			} else {
				result = Xn[0][i/2-1] + Xn[0][i/2-1]/(i+1);
				for (int k=1; k<i/2; k++) {
					result = result + binomial(i,k*2) * Xn[0][k-1] * Xn[0][i/2-k-1] / (k*2+1);
				}
				*it = result;
				it++;
			}
		}
		Xn.push_back(buf);
	} else {
		// calculate X_0
		std::vector<Number_T> buf(xninitsize/2);
		std::vector<Number_T>::iterator it = buf.begin();
		for (int i=1; i<=xninitsize/2; i++) {
			*it = bernoulli(i*2);
			it++;
		}
		Xn.push_back(buf);
	}

	xnsize++;
}

// doubles the number of entries in each Xn[]
void double_Xn()
{
	const int pos0 = xninitsize / 2;
	// X_0
	for (int i=1; i<=xninitsizestep/2; ++i) {
		Xn[0].push_back(bernoulli((i+pos0)*2);
	}
	if (Xn.size() > 1) {
		int xend = xninitsize + xninitsizestep;
		Number_T result;
		// X_1
		for (int i=xninitsize+1; i<=xend; ++i) {
			if (i & 1) {
				result = -Xn[0][(i-3)/2]/2;
				Xn[1].push_back((binomial(i,1)/Integer_T(2) + binomial(i,i-1)/Integer_T(i))*result);
			} else {
				result = Xn[0][i/2-1] + Xn[0][i/2-1]/(i+1);
				for (int k=1; k<i/2; k++) {
					result = result + binomial(i,k*2) * Xn[0][k-1] * Xn[0][i/2-k-1] / (k*2+1);
				}
				Xn[1].push_back(result);
			}
		}
		// X_n
		for (int n=2; n<Xn.size(); ++n) {
			for (int i=xninitsize+1; i<=xend; ++i) {
				if (i & 1) {
					result = 0; // k == 0
				} else {
					result = Xn[0][i/2-1]; // k == 0
				}
				for (int k=1; k<i-1; ++k) {
					if ( !(((i-k) & 1) && ((i-k) > 1)) ) {
						result = result + binomial(i,k) * Xn[0][(i-k)/2-1] * Xn[n-1][k-1] / (k+1);
					}
				}
				result = result - binomial(i,i-1) * Xn[n-1][i-2] / 2 / i; // k == i-1
				result = result + Xn[n-1][i-1] / (i+1); // k == i
				Xn[n].push_back(result);
			}
		}
	}
	xninitsize += xninitsizestep;
}


// calculates Li(2,x) without Xn
numeric Li2_do_sum(const numeric& x)
{
	Number_T res = x;
	Number_T resbuf;
	Number_T num = x ;
	Integer_T den = 1; // n^2 = 1
	unsigned i = 3;
	do {
		resbuf = res;
		num = num * x;
		den = den + i;  // n^2 = 4, 9, 16, ...
		i += 2;
		res = res + num / den;
	} while (res != resbuf);
	return res;
}

// calculates Li(2,x) with Xn
numeric Li2_do_sum_Xn(const numeric& x)
{
	std::vector<Number_T>::const_iterator it = Xn[0].begin();
	std::vector<Number_T>::const_iterator xend = Xn[0].end();
	Number_T u = -log(1-x);
	Number_T factor = u * ;
	Number_T uu = square(u);
	Number_T res = u - uu/4;
	Number_T resbuf;
	unsigned i = 1;
	do {
		resbuf = res;
		factor = factor * uu / (2*i * (2*i+1));
		res = res + (*it) * factor;
		i++;
		if (++it == xend) {
			double_Xn();
			it = Xn[0].begin() + (i-1);
			xend = Xn[0].end();
		}
	} while (res != resbuf);
	return res;
}


// calculates Li(n,x), n>2 without Xn
numeric Lin_do_sum(int n, const numeric& x)
{
  	Number_T factor = x;
	Number_T res = x;
	Number_T resbuf;
	int i=2;
	do {
		resbuf = res;
		factor = factor * x;
		res = res + factor / pow(Integer_T(i),n);
		i++;
	} while (res != resbuf);
	return res;
}

// calculates Li(n,x), n>2 with Xn
numeric Lin_do_sum_Xn(int n, const numeric& x)
{
// 	std::vector<Number_T>::const_iterator it = Xn[n-2].begin();
// 	std::vector<Number_T>::const_iterator xend = Xn[n-2].end();
// 	Number_T u = -cln::log(1-x);
// 	Number_T factor = u * cln::cl_float(1, cln::float_format(Digits));
// 	Number_T res = u;
// 	Number_T resbuf;
// 	unsigned i=2;
// 	do {
// 		resbuf = res;
// 		factor = factor * u / i;
// 		res = res + (*it) * factor;
// 		i++;
// 		if (++it == xend) {
// 			double_Xn();
// 			it = Xn[n-2].begin() + (i-2);
// 			xend = Xn[n-2].end();
// 		}
// 	} while (res != resbuf);
// 	return res;
}

// forward declaration needed by function Li_projection and C below
numeric S_num(int n, int p, const numeric& x);


// helper function for classical polylog Li
// Number_T Li_projection(int n, const Number_T& x, const cln::float_format_t& prec)
// {
// 	// treat n=2 as special case
// 	if (n == 2) {
// 		// check if precalculated X0 exists
// 		if (xnsize == 0) {
// 			fill_Xn(0);
// 		}

// 		if (cln::realpart(x) < 0.5) {
// 			// choose the faster algorithm
// 			// the switching point was empirically determined. the optimal point
// 			// depends on hardware, Digits, ... so an approx value is okay.
// 			// it solves also the problem with precision due to the u=-log(1-x) transformation
// 			if (cln::abs(cln::realpart(x)) < 0.25) {
				
// 				return Li2_do_sum(x);
// 			} else {
// 				return Li2_do_sum_Xn(x);
// 			}
// 		} else {
// 			// choose the faster algorithm
// 			if (cln::abs(cln::realpart(x)) > 0.75) {
// 				return -Li2_do_sum(1-x) - cln::log(x) * cln::log(1-x) + cln::zeta(2);
// 			} else {
// 				return -Li2_do_sum_Xn(1-x) - cln::log(x) * cln::log(1-x) + cln::zeta(2);
// 			}
// 		}
// 	} else {
// 		// check if precalculated Xn exist
// 		if (n > xnsize+1) {
// 			for (int i=xnsize; i<n-1; i++) {
// 				fill_Xn(i);
// 			}
// 		}

// 		if (cln::realpart(x) < 0.5) {
// 			// choose the faster algorithm
// 			// with n>=12 the "normal" summation always wins against the method with Xn
// 			if ((cln::abs(cln::realpart(x)) < 0.3) || (n >= 12)) {
// 				return Lin_do_sum(n, x);
// 			} else {
// 				return Lin_do_sum_Xn(n, x);
// 			}
// 		} else {
// 			Number_T result = -pow(log(x), n-1) * log(1-x) / factorial(n-1);
// 			for (int j=0; j<n-1; j++) {
// 				result = result + (S_num(n-j-1, 1, 1).to_cl_N() - S_num(1, n-j-1, 1-x).to_cl_N())
// 				                  * pow(log(x), j) / factorial(j);
// 			}
// 			return result;
// 		}
// 	}
// }
*/

// helper function for classical polylog Li
//numeric Lin_numeric(const numeric& n, const numeric& x, PyObject* parent)
//{
//  return Li2(x, n, parent);

// 	if (n == 1) {
// 		// just a log
// 		return -cln::log(1-x.to_cl_N());
// 	}
// 	if (x.is_zero()) {
// 		return 0;
// 	}
// 	if (x == 1) {
// 		// [Kol] (2.22)
// 		return cln::zeta(n);
// 	}
// 	else if (x == -1) {
// 		// [Kol] (2.22)
// 		return -(1-pow(Integer_T(2),1-n)) * cln::zeta(n);
// 	}
// 	if (abs(x.real()) < 0.4 && abs(abs(x)-1) < 0.01) {
// 		Number_T x_ = ex_to<numeric>(x).to_cl_N();
// 		Number_T result = -cln::expt(cln::log(x_), n-1) * cln::log(1-x_) / cln::factorial(n-1);
// 		for (int j=0; j<n-1; j++) {
// 			result = result + (S_num(n-j-1, 1, 1).to_cl_N() - S_num(1, n-j-1, 1-x_).to_cl_N())
// 				* cln::expt(cln::log(x_), j) / cln::factorial(j);
// 		}
// 		return result;
// 	}

// 	// what is the desired float format?
// 	// first guess: default format
// 	cln::float_format_t prec = cln::default_float_format;
// 	const Number_T value = x.to_cl_N();
// 	// second guess: the argument's format
// 	if (!x.real().is_rational())
// 		prec = cln::float_format(cln::the<cln::cl_F>(cln::realpart(value)));
// 	else if (!x.imag().is_rational())
// 		prec = cln::float_format(cln::the<cln::cl_F>(cln::imagpart(value)));
	
// 	// [Kol] (5.15)
// 	if (cln::abs(value) > 1) {
// 		Number_T result = -cln::expt(cln::log(-value),n) / cln::factorial(n);
// 		// check if argument is complex. if it is real, the new polylog has to be conjugated.
// 		if (cln::zerop(cln::imagpart(value))) {
// 			if (n & 1) {
// 				result = result + conjugate(Li_projection(n, cln::recip(value), prec));
// 			}
// 			else {
// 				result = result - conjugate(Li_projection(n, cln::recip(value), prec));
// 			}
// 		}
// 		else {
// 			if (n & 1) {
// 				result = result + Li_projection(n, cln::recip(value), prec);
// 			}
// 			else {
// 				result = result - Li_projection(n, cln::recip(value), prec);
// 			}
// 		}
// 		Number_T add;
// 		for (int j=0; j<n-1; j++) {
// 			add = add + (1+cln::expt(Integer_T(-1),n-j)) * (1-cln::expt(Integer_T(2),1-n+j))
// 			            * Lin_numeric(n-j,1).to_cl_N() * cln::expt(cln::log(-value),j) / cln::factorial(j);
// 		}
// 		result = result - add;
// 		return result;
// 	}
// 	else {
// 		return Li_projection(n, value, prec);
// 	}
//}


} // end of anonymous namespace


//////////////////////////////////////////////////////////////////////
//
// Multiple polylogarithm  Li(n,x)
//
// helper function
//
//////////////////////////////////////////////////////////////////////

/*

// anonymous namespace for helper function
namespace {


// performs the actual series summation for multiple polylogarithms
numeric multipleLi_do_sum(const std::vector<int>& s, const std::vector<numeric>& x)
{
// 	// ensure all x <> 0.
// 	for (std::vector<Number_T>::const_iterator it = x.begin(); it != x.end(); ++it) {
// 		if ( *it == 0 ) return cln::cl_float(0, cln::float_format(Digits));
// 	}

// 	const int j = s.size();
// 	bool flag_accidental_zero = false;

// 	std::vector<Number_T> t(j);
// 	cln::cl_F one = cln::cl_float(1, cln::float_format(Digits));

// 	Number_T t0buf;
// 	int q = 0;
// 	do {
// 		t0buf = t[0];
// 		q++;
// 		t[j-1] = t[j-1] + cln::expt(x[j-1], q) / cln::expt(Integer_T(q),s[j-1]) * one;
// 		for (int k=j-2; k>=0; k--) {
// 			t[k] = t[k] + t[k+1] * cln::expt(x[k], q+j-1-k) / cln::expt(Integer_T(q+j-1-k), s[k]);
// 		}
// 		q++;
// 		t[j-1] = t[j-1] + cln::expt(x[j-1], q) / cln::expt(Integer_T(q),s[j-1]) * one;
// 		for (int k=j-2; k>=0; k--) {
// 			flag_accidental_zero = cln::zerop(t[k+1]);
// 			t[k] = t[k] + t[k+1] * cln::expt(x[k], q+j-1-k) / cln::expt(Integer_T(q+j-1-k), s[k]);
// 		}
// 	} while ( (t[0] != t0buf) || cln::zerop(t[0]) || flag_accidental_zero );

// 	return t[0];
}

// converts parameter types and calls multipleLi_do_sum (convenience function for G_numeric)
numeric mLi_do_summation(const lst& m, const lst& x)
{
	std::vector<int> m_int;
	std::vector<Number_T> x_cln;
	for (lst::const_iterator itm = m.begin(), itx = x.begin(); itm != m.end(); ++itm, ++itx) {
		m_int.push_back(ex_to<numeric>(*itm).to_int());
		x_cln.push_back(ex_to<numeric>(*itx).to_cl_N());
	}
	return multipleLi_do_sum(m_int, x_cln);
}

// forward declaration for Li_eval()
lst convert_parameter_Li_to_H(const lst& m, const lst& x, ex& pf);


// holding dummy-symbols for the G/Li transformations
std::vector<ex> gsyms;


// type used by the transformation functions for G
typedef std::vector<int> Gparameter;


// G_eval1-function for G transformations
ex G_eval1(int a, int scale)
{
	if (a != 0) {
		const ex& scs = gsyms[std::abs(scale)];
		const ex& as = gsyms[std::abs(a)];
		if (as != scs) {
			return -log(1 - scs/as);
		} else {
			return -zeta(1);
		}
	} else {
		return log(gsyms[std::abs(scale)]);
	}
}

// G_eval-function for G transformations
ex G_eval(const Gparameter& a, int scale)
{
 	// check for properties of G
	ex sc = gsyms[std::abs(scale)];
	lst newa;
	bool all_zero = true;
	bool all_ones = true;
	int count_ones = 0;
	for (Gparameter::const_iterator it = a.begin(); it != a.end(); ++it) {
		if (*it != 0) {
			const ex sym = gsyms[std::abs(*it)];
			newa.append(sym);
			all_zero = false;
			if (sym != sc) {
				all_ones = false;
			}
			if (all_ones) {
				++count_ones;
			}
		} else {
			all_ones = false;
		}
	}

	// care about divergent G: shuffle to separate divergencies that will be canceled
	// later on in the transformation
	if (newa.nops() > 1 && newa.op(0) == sc && !all_ones && a.front()!=0) {
		// do shuffle
		Gparameter short_a;
		Gparameter::const_iterator it = a.begin();
		++it;
		for (; it != a.end(); ++it) {
			short_a.push_back(*it);
		}
		ex result = G_eval1(a.front(), scale) * G_eval(short_a, scale);
		it = short_a.begin();
		for (int i=1; i<count_ones; ++i) {
			++it;
		}
		for (; it != short_a.end(); ++it) {

			Gparameter newa;
			Gparameter::const_iterator it2 = short_a.begin();
			for (; it2 != it; ++it2) {
				newa.push_back(*it2);
			}
			newa.push_back(*it);
			newa.push_back(a[0]);
			it2 = it;
			++it2;
			for (; it2 != short_a.end(); ++it2) {
				newa.push_back(*it2);	
			}
			result -= G_eval(newa, scale);
		}
		return result / count_ones;
	}

	// G({1,...,1};y) -> G({1};y)^k / k!
	if (all_ones && a.size() > 1) {
		return pow(G_eval1(a.front(),scale), count_ones) / factorial(count_ones);
	}

	// G({0,...,0};y) -> log(y)^k / k!
	if (all_zero) {
		return pow(log(gsyms[std::abs(scale)]), a.size()) / factorial(a.size());
	}

	// no special cases anymore -> convert it into Li
	lst m;
	lst x;
	ex argbuf = gsyms[std::abs(scale)];
	ex mval = _ex1;
	for (Gparameter::const_iterator it=a.begin(); it!=a.end(); ++it) {
		if (*it != 0) {
			const ex& sym = gsyms[std::abs(*it)];
			x.append(argbuf / sym);
			m.append(mval);
			mval = _ex1;
			argbuf = sym;
		} else {
			++mval;
		}
	}
	return pow(-1, x.nops()) * Li(m, x);
}

// converts data for G: pending_integrals -> a
Gparameter convert_pending_integrals_G(const Gparameter& pending_integrals)
{
	GINAC_ASSERT(pending_integrals.size() != 1);

	if (pending_integrals.size() > 0) {
		// get rid of the first element, which would stand for the new upper limit
		Gparameter new_a(pending_integrals.begin()+1, pending_integrals.end());
		return new_a;
	} else {
		// just return empty parameter list
		Gparameter new_a;
		return new_a;
	}
}

// check the parameters a and scale for G and return information about convergence, depth, etc.
// convergent     : true if G(a,scale) is convergent
// depth          : depth of G(a,scale)
// trailing_zeros : number of trailing zeros of a
// min_it         : iterator of a pointing on the smallest element in a
Gparameter::const_iterator check_parameter_G(const Gparameter& a, int scale,
		bool& convergent, int& depth, int& trailing_zeros, Gparameter::const_iterator& min_it)
{
	convergent = true;
	depth = 0;
	trailing_zeros = 0;
	min_it = a.end();
	Gparameter::const_iterator lastnonzero = a.end();
	for (Gparameter::const_iterator it = a.begin(); it != a.end(); ++it) {
		if (std::abs(*it) > 0) {
			++depth;
			trailing_zeros = 0;
			lastnonzero = it;
			if (std::abs(*it) < scale) {
				convergent = false;
				if ((min_it == a.end()) || (std::abs(*it) < std::abs(*min_it))) {
					min_it = it;
				}
			}
		} else {
			++trailing_zeros;
		}
	}
	if (lastnonzero == a.end())
		return a.end();
	return ++lastnonzero;
}

// add scale to pending_integrals if pending_integrals is empty
Gparameter prepare_pending_integrals(const Gparameter& pending_integrals, int scale)
{
	GINAC_ASSERT(pending_integrals.size() != 1);

	if (pending_integrals.size() > 0) {
		return pending_integrals;
	} else {
		Gparameter new_pending_integrals;
		new_pending_integrals.push_back(scale);
		return new_pending_integrals;
	}
}

// handles trailing zeroes for an otherwise convergent integral
ex trailing_zeros_G(const Gparameter& a, int scale)
{
	bool convergent;
	int depth, trailing_zeros;
	Gparameter::const_iterator last, dummyit;
	last = check_parameter_G(a, scale, convergent, depth, trailing_zeros, dummyit);

	GINAC_ASSERT(convergent);

	if ((trailing_zeros > 0) && (depth > 0)) {
		ex result;
		Gparameter new_a(a.begin(), a.end()-1);
		result += G_eval1(0, scale) * trailing_zeros_G(new_a, scale);
		for (Gparameter::const_iterator it = a.begin(); it != last; ++it) {
			Gparameter new_a(a.begin(), it);
			new_a.push_back(0);
			new_a.insert(new_a.end(), it, a.end()-1);
			result -= trailing_zeros_G(new_a, scale);
		}

		return result / trailing_zeros;
	} else {
		return G_eval(a, scale);
	}
}

// G transformation [VSW] (57),(58)
ex depth_one_trafo_G(const Gparameter& pending_integrals, const Gparameter& a, int scale)
{
	// pendint = ( y1, b1, ..., br )
	//       a = ( 0, ..., 0, amin )
	//   scale = y2
	//
	// int_0^y1 ds1/(s1-b1) ... int dsr/(sr-br) G(0, ..., 0, sr; y2)
	// where sr replaces amin

	GINAC_ASSERT(a.back() != 0);
	GINAC_ASSERT(a.size() > 0);

	ex result;
	Gparameter new_pending_integrals = prepare_pending_integrals(pending_integrals, std::abs(a.back()));
	const int psize = pending_integrals.size();

	// length == 1
	// G(sr_{+-}; y2 ) = G(y2_{-+}; sr) - G(0; sr) + ln(-y2_{-+})

	if (a.size() == 1) {

	  // ln(-y2_{-+})
	  result += log(gsyms[ex_to<numeric>(scale).to_int()]);
		if (a.back() > 0) {
			new_pending_integrals.push_back(-scale);
			result += I*Pi;
		} else {
			new_pending_integrals.push_back(scale);
			result -= I*Pi;
		}
		if (psize) {
			result *= trailing_zeros_G(convert_pending_integrals_G(pending_integrals), pending_integrals.front());
		}
		
		// G(y2_{-+}; sr)
		result += trailing_zeros_G(convert_pending_integrals_G(new_pending_integrals), new_pending_integrals.front());
		
		// G(0; sr)
		new_pending_integrals.back() = 0;
		result -= trailing_zeros_G(convert_pending_integrals_G(new_pending_integrals), new_pending_integrals.front());

		return result;
	}

	// length > 1
	// G_m(sr_{+-}; y2) = -zeta_m + int_0^y2 dt/t G_{m-1}( (1/y2)_{+-}; 1/t )
	//                            - int_0^sr dt/t G_{m-1}( (1/y2)_{+-}; 1/t )

	//term zeta_m
	result -= zeta(a.size());
	if (psize) {
		result *= trailing_zeros_G(convert_pending_integrals_G(pending_integrals), pending_integrals.front());
	}
	
	// term int_0^sr dt/t G_{m-1}( (1/y2)_{+-}; 1/t )
	//    = int_0^sr dt/t G_{m-1}( t_{+-}; y2 )
	Gparameter new_a(a.begin()+1, a.end());
	new_pending_integrals.push_back(0);
	result -= depth_one_trafo_G(new_pending_integrals, new_a, scale);
	
	// term int_0^y2 dt/t G_{m-1}( (1/y2)_{+-}; 1/t )
	//    = int_0^y2 dt/t G_{m-1}( t_{+-}; y2 )
	Gparameter new_pending_integrals_2;
	new_pending_integrals_2.push_back(scale);
	new_pending_integrals_2.push_back(0);
	if (psize) {
		result += trailing_zeros_G(convert_pending_integrals_G(pending_integrals), pending_integrals.front())
		          * depth_one_trafo_G(new_pending_integrals_2, new_a, scale);
	} else {
		result += depth_one_trafo_G(new_pending_integrals_2, new_a, scale);
	}

	return result;
}

// forward declaration
ex shuffle_G(const Gparameter & a0, const Gparameter & a1, const Gparameter & a2,
	     const Gparameter& pendint, const Gparameter& a_old, int scale);


// G transformation [VSW]
ex G_transform(const Gparameter& pendint, const Gparameter& a, int scale)
{
	// main recursion routine
	//
	// pendint = ( y1, b1, ..., br )
	//       a = ( a1, ..., amin, ..., aw )
	//   scale = y2
	//
	// int_0^y1 ds1/(s1-b1) ... int dsr/(sr-br) G(a1,...,sr,...,aw,y2)
	// where sr replaces amin

	// find smallest alpha, determine depth and trailing zeros, and check for convergence
	bool convergent;
	int depth, trailing_zeros;
	Gparameter::const_iterator min_it;
	Gparameter::const_iterator firstzero = 
		check_parameter_G(a, scale, convergent, depth, trailing_zeros, min_it);
	int min_it_pos = min_it - a.begin();

	// special case: all a's are zero
	if (depth == 0) {
		ex result;

		if (a.size() == 0) {
		  result = 1;
		} else {
		  result = G_eval(a, scale);
		}
		if (pendint.size() > 0) {
		  result *= trailing_zeros_G(convert_pending_integrals_G(pendint), pendint.front());
		} 
		return result;
	}

	// handle trailing zeros
	if (trailing_zeros > 0) {
		ex result;
		Gparameter new_a(a.begin(), a.end()-1);
		result += G_eval1(0, scale) * G_transform(pendint, new_a, scale);
		for (Gparameter::const_iterator it = a.begin(); it != firstzero; ++it) {
			Gparameter new_a(a.begin(), it);
			new_a.push_back(0);
			new_a.insert(new_a.end(), it, a.end()-1);
			result -= G_transform(pendint, new_a, scale);
		}
		return result / trailing_zeros;
	}

	// convergence case
	if (convergent) {
		if (pendint.size() > 0) {
			return G_eval(convert_pending_integrals_G(pendint), pendint.front()) * G_eval(a, scale);
		} else {
			return G_eval(a, scale);
		}
	}

	// call basic transformation for depth equal one
	if (depth == 1) {
		return depth_one_trafo_G(pendint, a, scale);
	}

	// do recursion
	// int_0^y1 ds1/(s1-b1) ... int dsr/(sr-br) G(a1,...,sr,...,aw,y2)
	//  =  int_0^y1 ds1/(s1-b1) ... int dsr/(sr-br) G(a1,...,0,...,aw,y2)
	//   + int_0^y1 ds1/(s1-b1) ... int dsr/(sr-br) int_0^{sr} ds_{r+1} d/ds_{r+1} G(a1,...,s_{r+1},...,aw,y2)

	// smallest element in last place
	if (min_it + 1 == a.end()) {
		do { --min_it; } while (*min_it == 0);
		Gparameter empty;
		Gparameter a1(a.begin(),min_it+1);
		Gparameter a2(min_it+1,a.end());

		ex result = G_transform(pendint,a2,scale)*G_transform(empty,a1,scale);

		result -= shuffle_G(empty,a1,a2,pendint,a,scale);
		return result;
	}

	Gparameter empty;
	Gparameter::iterator changeit;

	// first term G(a_1,..,0,...,a_w;a_0)
	Gparameter new_pendint = prepare_pending_integrals(pendint, a[min_it_pos]);
	Gparameter new_a = a;
	new_a[min_it_pos] = 0;
	ex result = G_transform(empty, new_a, scale);
	if (pendint.size() > 0) {
		result *= trailing_zeros_G(convert_pending_integrals_G(pendint), pendint.front());
	}

	// other terms
	changeit = new_a.begin() + min_it_pos;
	changeit = new_a.erase(changeit);
	if (changeit != new_a.begin()) {
		// smallest in the middle
		new_pendint.push_back(*changeit);
		result -= trailing_zeros_G(convert_pending_integrals_G(new_pendint), new_pendint.front())
			* G_transform(empty, new_a, scale);
		int buffer = *changeit;
		*changeit = *min_it;
		result += G_transform(new_pendint, new_a, scale);
		*changeit = buffer;
		new_pendint.pop_back();
		--changeit;
		new_pendint.push_back(*changeit);
		result += trailing_zeros_G(convert_pending_integrals_G(new_pendint), new_pendint.front())
			* G_transform(empty, new_a, scale);
		*changeit = *min_it;
		result -= G_transform(new_pendint, new_a, scale);
	} else {
		// smallest at the front
		new_pendint.push_back(scale);
		result += trailing_zeros_G(convert_pending_integrals_G(new_pendint), new_pendint.front())
			* G_transform(empty, new_a, scale);
		new_pendint.back() =  *changeit;
		result -= trailing_zeros_G(convert_pending_integrals_G(new_pendint), new_pendint.front())
			* G_transform(empty, new_a, scale);
		*changeit = *min_it;
		result += G_transform(new_pendint, new_a, scale);
	}
	return result;
}

// shuffles the two parameter list a1 and a2 and calls G_transform for every term except
// for the one that is equal to a_old
ex shuffle_G(const Gparameter & a0, const Gparameter & a1, const Gparameter & a2,
	     const Gparameter& pendint, const Gparameter& a_old, int scale) 
{
	if (a1.size()==0 && a2.size()==0) {
		// veto the one configuration we don't want
		if ( a0 == a_old ) return 0;

		return G_transform(pendint,a0,scale);
	}

	if (a2.size()==0) {
		Gparameter empty;
		Gparameter aa0 = a0;
		aa0.insert(aa0.end(),a1.begin(),a1.end());
		return shuffle_G(aa0,empty,empty,pendint,a_old,scale);
	}

	if (a1.size()==0) {
		Gparameter empty;
		Gparameter aa0 = a0;
		aa0.insert(aa0.end(),a2.begin(),a2.end());
		return shuffle_G(aa0,empty,empty,pendint,a_old,scale);
	}

	Gparameter a1_removed(a1.begin()+1,a1.end());
	Gparameter a2_removed(a2.begin()+1,a2.end());

	Gparameter a01 = a0;
	Gparameter a02 = a0;

	a01.push_back( a1[0] );
	a02.push_back( a2[0] );

	return shuffle_G(a01,a1_removed,a2,pendint,a_old,scale)
	     + shuffle_G(a02,a1,a2_removed,pendint,a_old,scale);
}

// handles the transformations and the numerical evaluation of G
// the parameter x, s and y must only contain numerics
ex G_numeric(const lst& x, const lst& s, const ex& y)
{
	// check for convergence and necessary accelerations
	bool need_trafo = false;
	bool need_hoelder = false;
	int depth = 0;
	for (lst::const_iterator it = x.begin(); it != x.end(); ++it) {
		if (!(*it).is_zero()) {
			++depth;
			if (abs(*it) - y < -pow(10,-Digits+1)) {
				need_trafo = true;
			}
			if (abs((abs(*it) - y)/y) < 0.01) {
				need_hoelder = true;
			}
		}
	}
	if (x.op(x.nops()-1).is_zero()) {
		need_trafo = true;
	}
	if (depth == 1 && x.nops() == 2 && !need_trafo) {
		return -Li(x.nops(), y / x.op(x.nops()-1)).evalf(0, prec);
	}
	
	// do acceleration transformation (hoelder convolution [BBB])
	if (need_hoelder) {
		
		ex result;
		const int size = x.nops();
		lst newx;
		for (lst::const_iterator it = x.begin(); it != x.end(); ++it) {
			newx.append(*it / y);
		}
		
		for (int r=0; r<=size; ++r) {
			ex buffer = pow(-1, r);
			ex p = 2;
			bool adjustp;
			do {
				adjustp = false;
				for (lst::const_iterator it = newx.begin(); it != newx.end(); ++it) {
					if (*it == 1/p) {
						p += (3-p)/2; 
						adjustp = true;
						continue;
					}
				}
			} while (adjustp);
			ex q = p / (p-1);
			lst qlstx;
			lst qlsts;
			for (int j=r; j>=1; --j) {
				qlstx.append(1-newx.op(j-1));
				if (newx.op(j-1).info(info_flags::real) && newx.op(j-1) > 1 && newx.op(j-1) <= 2) {
					qlsts.append( s.op(j-1));
				} else {
					qlsts.append( -s.op(j-1));
				}
			}
			if (qlstx.nops() > 0) {
				buffer *= G_numeric(qlstx, qlsts, 1/q);
			}
			lst plstx;
			lst plsts;
			for (int j=r+1; j<=size; ++j) {
				plstx.append(newx.op(j-1));
				plsts.append(s.op(j-1));
			}
			if (plstx.nops() > 0) {
				buffer *= G_numeric(plstx, plsts, 1/p);
			}
			result += buffer;
		}
		return result;
	}
	
	// convergence transformation
	if (need_trafo) {

		// sort (|x|<->position) to determine indices
		std::multimap<ex,int> sortmap;
		int size = 0;
		for (int i=0; i<x.nops(); ++i) {
			if (!x[i].is_zero()) {
				sortmap.insert(std::pair<ex,int>(abs(x[i]), i));
				++size;
			}
		}
		// include upper limit (scale)
		sortmap.insert(std::pair<ex,int>(abs(y), x.nops()));

		// generate missing dummy-symbols
		int i = 1;
		gsyms.clear();
		gsyms.push_back(symbol("GSYMS_ERROR"));
		ex lastentry;
		for (std::multimap<ex,int>::const_iterator it = sortmap.begin(); it != sortmap.end(); ++it) {
			if (it != sortmap.begin()) {
				if (it->second < x.nops()) {
					if (x[it->second] == lastentry) {
						gsyms.push_back(gsyms.back());
						continue;
					}
				} else {
					if (y == lastentry) {
						gsyms.push_back(gsyms.back());
						continue;
					}
				}
			}
			std::ostringstream os;
			os << "a" << i;
			gsyms.push_back(symbol(os.str()));
			++i;
			if (it->second < x.nops()) {
				lastentry = x[it->second];
			} else {
				lastentry = y;
			}
		}

		// fill position data according to sorted indices and prepare substitution list
		Gparameter a(x.nops());
		lst subslst;
		int pos = 1;
		int scale;
		for (std::multimap<ex,int>::const_iterator it = sortmap.begin(); it != sortmap.end(); ++it) {
			if (it->second < x.nops()) {
				if (s[it->second] > 0) {
					a[it->second] = pos;
				} else {
					a[it->second] = -pos;
				}
				subslst.append(gsyms[pos] == x[it->second]);
			} else {
				scale = pos;
				subslst.append(gsyms[pos] == y);
			}
			++pos;
		}

		// do transformation
		Gparameter pendint;
		ex result = G_transform(pendint, a, scale);
		// replace dummy symbols with their values
		result = result.eval().expand();
		result = result.subs(subslst).evalf(0, prec);
		
		return result;
	}

	// do summation
	lst newx;
	lst m;
	int mcount = 1;
	ex sign = 1;
	ex factor = y;
	for (lst::const_iterator it = x.begin(); it != x.end(); ++it) {
		if ((*it).is_zero()) {
			++mcount;
		} else {
			newx.append(factor / (*it));
			factor = *it;
			m.append(mcount);
			mcount = 1;
			sign = -sign;
		}
	}

	return sign * numeric(mLi_do_summation(m, newx));
}

ex mLi_numeric(const lst& m, const lst& x)
{
	// let G_numeric do the transformation
	lst newx;
	lst s;
	ex factor = 1;
	for (lst::const_iterator itm = m.begin(), itx = x.begin(); itm != m.end(); ++itm, ++itx) {
		for (int i = 1; i < *itm; ++i) {
			newx.append(0);
			s.append(1);
		}
		newx.append(factor / *itx);
		factor /= *itx;
		s.append(1);
	}
	return pow(-1, m.nops()) * G_numeric(newx, s, _ex1);
}
        
} // end of anonymous namespace
*/


//////////////////////////////////////////////////////////////////////
//
// Generalized multiple polylogarithm  G(x, y) and G(x, s, y)
//
// GiNaC function
//
//////////////////////////////////////////////////////////////////////

/*

static ex G2_evalf(const ex& x_, const ex& y, PyObject* parent)
{
	if (!y.info(info_flags::positive)) {
		return G(x_, y).hold();
	}
	lst x = is_a<lst>(x_) ? ex_to<lst>(x_) : lst(x_);
	if (x.nops() == 0) {
		return _ex1;
	}
	if (x.op(0) == y) {
		return G(x_, y).hold();
	}
	lst s;
	bool all_zero = true;
	for (lst::const_iterator it = x.begin(); it != x.end(); ++it) {
		if (!(*it).info(info_flags::numeric)) {
			return G(x_, y).hold();
		}
		if (*it != _ex0) {
			all_zero = false;
		}
		s.append(+1);
	}
	if (all_zero) {
		return pow(log(y), x.nops()) / factorial(x.nops());
	}
	return G_numeric(x, s, y);
}

static ex G2_eval(const ex& x_, const ex& y)
{
	//TODO eval to MZV or H or S or Lin
	if (!y.info(info_flags::positive)) {
		return G(x_, y).hold();
	}
	lst x = is_a<lst>(x_) ? ex_to<lst>(x_) : lst(x_);
	if (x.nops() == 0) {
		return _ex1;
	}
	if (x.op(0) == y) {
		return G(x_, y).hold();
	}
	lst s;
	bool all_zero = true;
	bool crational = true;
	for (lst::const_iterator it = x.begin(); it != x.end(); ++it) {
		if (!(*it).info(info_flags::numeric)) {
			return G(x_, y).hold();
		}
		if (!(*it).info(info_flags::crational)) {
			crational = false;
		}
		if (*it != _ex0) {
			all_zero = false;
		}
		s.append(+1);
	}
	if (all_zero) {
		return pow(log(y), x.nops()) / factorial(x.nops());
	}
	if (!y.info(info_flags::crational)) {
		crational = false;
	}
	if (crational) {
		return G(x_, y).hold();
	}
	return G_numeric(x, s, y);
}


unsigned G2_SERIAL::serial = function::register_new(function_options("G", 2).
                                evalf_func(G2_evalf).
                                eval_func(G2_eval).
                                do_not_evalf_params().
                                overloaded(2));
//TODO
//                                derivative_func(G2_deriv).
//                                print_func<print_latex>(G2_print_latex).


static ex G3_evalf(const ex& x_, const ex& s_, const ex& y, PyObject* parent)
{
	if (!y.info(info_flags::positive)) {
		return G(x_, s_, y).hold();
	}
	lst x = is_a<lst>(x_) ? ex_to<lst>(x_) : lst(x_);
	lst s = is_a<lst>(s_) ? ex_to<lst>(s_) : lst(s_);
	if (x.nops() != s.nops()) {
		return G(x_, s_, y).hold();
	}
	if (x.nops() == 0) {
		return _ex1;
	}
	if (x.op(0) == y) {
		return G(x_, s_, y).hold();
	}
	lst sn;
	bool all_zero = true;
	for (lst::const_iterator itx = x.begin(), its = s.begin(); itx != x.end(); ++itx, ++its) {
		if (!(*itx).info(info_flags::numeric)) {
			return G(x_, y).hold();
		}
		if (!(*its).info(info_flags::real)) {
			return G(x_, y).hold();
		}
		if (*itx != _ex0) {
			all_zero = false;
		}
		if (*its >= 0) {
			sn.append(+1);
		} else {
			sn.append(-1);
		}
	}
	if (all_zero) {
		return pow(log(y), x.nops()) / factorial(x.nops());
	}
	return G_numeric(x, sn, y);
}


static ex G3_eval(const ex& x_, const ex& s_, const ex& y)
{
	//TODO eval to MZV or H or S or Lin
	if (!y.info(info_flags::positive)) {
		return G(x_, s_, y).hold();
	}
	lst x = is_a<lst>(x_) ? ex_to<lst>(x_) : lst(x_);
	lst s = is_a<lst>(s_) ? ex_to<lst>(s_) : lst(s_);
	if (x.nops() != s.nops()) {
		return G(x_, s_, y).hold();
	}
	if (x.nops() == 0) {
		return _ex1;
	}
	if (x.op(0) == y) {
		return G(x_, s_, y).hold();
	}
	lst sn;
	bool all_zero = true;
	bool crational = true;
	for (lst::const_iterator itx = x.begin(), its = s.begin(); itx != x.end(); ++itx, ++its) {
		if (!(*itx).info(info_flags::numeric)) {
			return G(x_, s_, y).hold();
		}
		if (!(*its).info(info_flags::real)) {
			return G(x_, s_, y).hold();
		}
		if (!(*itx).info(info_flags::crational)) {
			crational = false;
		}
		if (*itx != _ex0) {
			all_zero = false;
		}
		if (*its >= 0) {
			sn.append(+1);
		} else {
			sn.append(-1);
		}
	}
	if (all_zero) {
		return pow(log(y), x.nops()) / factorial(x.nops());
	}
	if (!y.info(info_flags::crational)) {
		crational = false;
	}
	if (crational) {
		return G(x_, s_, y).hold();
	}
	return G_numeric(x, sn, y);
}


unsigned G3_SERIAL::serial = function::register_new(function_options("G", 3).
                                evalf_func(G3_evalf).
                                eval_func(G3_eval).
                                do_not_evalf_params().
                                overloaded(2));
//TODO
//                                derivative_func(G3_deriv).
//                                print_func<print_latex>(G3_print_latex).

 */

//////////////////////////////////////////////////////////////////////
//
// Classical polylogarithm and multiple polylogarithm  Li(m,x)
//
// GiNaC function
//
//////////////////////////////////////////////////////////////////////

#if 0
static ex Li_evalf(const ex& m_, const ex& x_, PyObject* parent)
{
	// classical polylogs
	if (is_exactly_a<numeric>(m_) and m_.info(info_flags::posint)) {
		if (is_exactly_a<numeric>(x_)) {
			return Lin_numeric(ex_to<numeric>(m_), ex_to<numeric>(x_), parent);
		} else {
			// try to numerically evaluate second argument
			ex x_val = x_.evalf(0, parent);
			if (is_exactly_a<numeric>(x_val)) {
				return Lin_numeric(ex_to<numeric>(m_), ex_to<numeric>(x_val), parent);
			}
		}
	}
	// multiple polylogs
        /*
	if (is_a<lst>(m_) && is_a<lst>(x_)) {

		const lst& m = ex_to<lst>(m_);
		const lst& x = ex_to<lst>(x_);
		if (m.nops() != x.nops()) {
			return Li(m_,x_).hold();
		}
		if (x.nops() == 0) {
			return _ex1;
		}
		if ((m.op(0) == _ex1) && (x.op(0) == _ex1)) {
			return Li(m_,x_).hold();
		}

		for (auto itm = m.begin(), itx = x.begin(); itm != m.end(); ++itm, ++itx) {
			if (!(*itm).info(info_flags::posint)) {
				return Li(m_, x_).hold();
			}
			if (!(*itx).info(info_flags::numeric)) {
				return Li(m_, x_).hold();
			}
			if (*itx == _ex0) {
				return _ex0;
			}
		}

		return mLi_numeric(m, x);
	}*/

	return Li(m_,x_).hold();
}


static ex Li_eval(const ex& m_, const ex& x_)
{
        /*
	if (is_a<lst>(m_)) {
		if (is_a<lst>(x_)) {
			// multiple polylogs
			const lst& m = ex_to<lst>(m_);
			const lst& x = ex_to<lst>(x_);
			if (m.nops() != x.nops()) {
				return Li(m_,x_).hold();
			}
			if (x.nops() == 0) {
				return _ex1;
			}
			bool is_H = true;
			bool is_zeta = true;
			bool do_evalf = true;
			bool crational = true;
			for (auto itm = m.begin(), itx = x.begin(); itm != m.end(); ++itm, ++itx) {
				if (!(*itm).info(info_flags::posint)) {
					return Li(m_,x_).hold();
				}
				if ((*itx != _ex1) && (*itx != _ex_1)) {
					if (itx != x.begin()) {
						is_H = false;
					}
					is_zeta = false;
				}
				if (*itx == _ex0) {
					return _ex0;
				}
				if (!(*itx).info(info_flags::numeric)) {
					do_evalf = false;
				}
				if (!(*itx).info(info_flags::crational)) {
					crational = false;
				}
			}
			if (is_zeta) {
				return zeta(m_,x_);
			}
			if (is_H) {
				ex prefactor;
				lst newm = convert_parameter_Li_to_H(m, x, prefactor);
				return prefactor * H(newm, x[0]);
			}
			if (do_evalf && !crational) {
				return mLi_numeric(m,x);
			}
		}
		return Li(m_, x_).hold();
	} else if (is_a<lst>(x_)) {
		return Li(m_, x_).hold();
	}*/

	// classical polylogs
	if (x_ == _ex0) {
		return _ex0;
	}
	if (x_ == _ex1) {
		return zeta(m_);
	}
	if (x_ == _ex_1) {
		return (pow(2,1-m_)-1) * zeta(m_);
	}
	if (m_ == _ex1) {
		return -log(1-x_);
	}
	if (m_ == _ex2) {
		if (x_.is_equal(I)) {
			return power(Pi,_ex2)/_ex_48 + Catalan*I;
		}
		if (x_.is_equal(-I)) {
			return power(Pi,_ex2)/_ex_48 - Catalan*I;
		}
	}
	/* if (m_.info(info_flags::posint) && x_.info(info_flags::numeric) && !x_.info(info_flags::crational)) { */
	/* 	return Lin_numeric(ex_to<numeric>(m_), ex_to<numeric>(x_)); */
	/* } */

	return Li(m_, x_).hold();
}


static ex Li_series(const ex& m, const ex& x, const relational& rel, int order, unsigned options)
{
	if (is_a<lst>(m) || is_a<lst>(x)) {
		// multiple polylog
		epvector seq;
		seq.push_back(expair(Li(m, x), 0));
		return pseries(rel, seq);
	}
	
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
			nseq.push_back(expair(Order(_ex1), order));
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
	if (m_.nops() > 1) {
		throw std::runtime_error("don't know how to derivate multiple polylogarithm!");
	}
	ex m;
	if (is_a<lst>(m_)) {
		m = m_.op(0);
	} else {
		m = m_;
	}
	ex x;
	if (is_a<lst>(x_)) {
		x = x_.op(0);
	} else {
		x = x_;
	}
	if (m > 0) {
		return Li(m-1, x) / x;
	} else {
		return 1/(1-x);
	}
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
#endif

//////////////////////////////////////////////////////////////////////
//
// Nielsen's generalized polylogarithm  S(n,p,x)
//
// helper functions
//
//////////////////////////////////////////////////////////////////////

/*

// anonymous namespace for helper functions
namespace {


// lookup table for special Euler-Zagier-Sums (used for S_n,p(x))
// see fill_Yn()
//std::vector<std::vector<Number_T> > Yn;
//int ynsize = 0; // number of Yn[]
//int ynlength = 100; // initial length of all Yn[i]


// This function calculates the Y_n. The Y_n are needed for the evaluation of S_{n,p}(x).
// The Y_n are basically Euler-Zagier sums with all m_i=1. They are subsums in the Z-sum
// representing S_{n,p}(x).
// The first index in Y_n corresponds to the parameter p minus one, i.e. the depth of the
// equivalent Z-sum.
// The second index in Y_n corresponds to the running index of the outermost sum in the full Z-sum
// representing S_{n,p}(x).
// The calculation of Y_n uses the values from Y_{n-1}.
// void fill_Yn(int n, const cln::float_format_t& prec)
// {
// 	const int initsize = ynlength;
// 	//const int initsize = initsize_Yn;
// 	Number_T one = cln::cl_float(1, prec);

// 	if (n) {
// 		std::vector<Number_T> buf(initsize);
// 		std::vector<Number_T>::iterator it = buf.begin();
// 		std::vector<Number_T>::iterator itprev = Yn[n-1].begin();
// 		*it = (*itprev) / Number_T(n+1) * one;
// 		it++;
// 		itprev++;
// 		// sums with an index smaller than the depth are zero and need not to be calculated.
// 		// calculation starts with depth, which is n+2)
// 		for (int i=n+2; i<=initsize+n; i++) {
// 			*it = *(it-1) + (*itprev) / Number_T(i) * one;
// 			it++;
// 			itprev++;
// 		}
// 		Yn.push_back(buf);
// 	} else {
// 		std::vector<Number_T> buf(initsize);
// 		std::vector<Number_T>::iterator it = buf.begin();
// 		*it = 1 * one;
// 		it++;
// 		for (int i=2; i<=initsize; i++) {
// 			*it = *(it-1) + 1 / Number_T(i) * one;
// 			it++;
// 		}
// 		Yn.push_back(buf);
// 	}
// 	ynsize++;
// }


// make Yn longer ... 
// void make_Yn_longer(int newsize, const cln::float_format_t& prec)
// {

// 	Number_T one = cln::cl_float(1, prec);

// 	Yn[0].resize(newsize);
// 	std::vector<Number_T>::iterator it = Yn[0].begin();
// 	it += ynlength;
// 	for (int i=ynlength+1; i<=newsize; i++) {
// 		*it = *(it-1) + 1 / Number_T(i) * one;
// 		it++;
// 	}

// 	for (int n=1; n<ynsize; n++) {
// 		Yn[n].resize(newsize);
// 		std::vector<Number_T>::iterator it = Yn[n].begin();
// 		std::vector<Number_T>::iterator itprev = Yn[n-1].begin();
// 		it += ynlength;
// 		itprev += ynlength;
// 		for (int i=ynlength+n+1; i<=newsize+n; i++) {
// 			*it = *(it-1) + (*itprev) / Number_T(i) * one;
// 			it++;
// 			itprev++;
// 		}
// 	}
	
// 	ynlength = newsize;
// }


// helper function for S(n,p,x)
// [Kol] (7.2)
// Number_T C(int n, int p)
// {
// 	Number_T result;

// 	for (int k=0; k<p; k++) {
// 		for (int j=0; j<=(n+k-1)/2; j++) {
// 			if (k == 0) {
// 				if (n & 1) {
// 					if (j & 1) {
// 						result = result - 2 * cln::expt(cln::pi(),2*j) * S_num(n-2*j,p,1).to_cl_N() / cln::factorial(2*j);
// 					}
// 					else {
// 						result = result + 2 * cln::expt(cln::pi(),2*j) * S_num(n-2*j,p,1).to_cl_N() / cln::factorial(2*j);
// 					}
// 				}
// 			}
// 			else {
// 				if (k & 1) {
// 					if (j & 1) {
// 						result = result + cln::factorial(n+k-1)
// 						                  * cln::expt(cln::pi(),2*j) * S_num(n+k-2*j,p-k,1).to_cl_N()
// 						                  / (cln::factorial(k) * cln::factorial(n-1) * cln::factorial(2*j));
// 					}
// 					else {
// 						result = result - cln::factorial(n+k-1)
// 						                  * cln::expt(cln::pi(),2*j) * S_num(n+k-2*j,p-k,1).to_cl_N()
// 						                  / (cln::factorial(k) * cln::factorial(n-1) * cln::factorial(2*j));
// 					}
// 				}
// 				else {
// 					if (j & 1) {
// 						result = result - cln::factorial(n+k-1) * cln::expt(cln::pi(),2*j) * S_num(n+k-2*j,p-k,1).to_cl_N()
// 						                  / (cln::factorial(k) * cln::factorial(n-1) * cln::factorial(2*j));
// 					}
// 					else {
// 						result = result + cln::factorial(n+k-1)
// 						                  * cln::expt(cln::pi(),2*j) * S_num(n+k-2*j,p-k,1).to_cl_N()
// 						                  / (cln::factorial(k) * cln::factorial(n-1) * cln::factorial(2*j));
// 					}
// 				}
// 			}
// 		}
// 	}
// 	int np = n+p;
// 	if ((np-1) & 1) {
// 		if (((np)/2+n) & 1) {
// 			result = -result - cln::expt(cln::pi(),np) / (np * cln::factorial(n-1) * cln::factorial(p));
// 		}
// 		else {
// 			result = -result + cln::expt(cln::pi(),np) / (np * cln::factorial(n-1) * cln::factorial(p));
// 		}
// 	}

// 	return result;
// }


// helper function for S(n,p,x)
// [Kol] remark to (9.1)
numeric a_k(int k)
{
// 	Number_T result;

// 	if (k == 0) {
// 		return 1;
// 	}

// 	result = result;
// 	for (int m=2; m<=k; m++) {
// 		result = result + cln::expt(Number_T(-1),m) * cln::zeta(m) * a_k(k-m);
// 	}

// 	return -result / k;
}


// helper function for S(n,p,x)
// [Kol] remark to (9.1)
numeric b_k(int k)
{
// 	Number_T result;

// 	if (k == 0) {
// 		return 1;
// 	}

// 	result = result;
// 	for (int m=2; m<=k; m++) {
// 		result = result + cln::expt(Number_T(-1),m) * cln::zeta(m) * b_k(k-m);
// 	}

// 	return result / k;
}


// helper function for S(n,p,x)
//Number_T S_do_sum(int n, int p, const Number_T& x, const cln::float_format_t& prec)
//{
// 	static cln::float_format_t oldprec = cln::default_float_format;

// 	if (p==1) {
// 		return Li_projection(n+1, x, prec);
// 	}

// 	// precision has changed, we need to clear lookup table Yn
// 	if ( oldprec != prec ) {
// 		Yn.clear();
// 		ynsize = 0;
// 		ynlength = 100;
// 		oldprec = prec;
// 	}
		
// 	// check if precalculated values are sufficient
// 	if (p > ynsize+1) {
// 		for (int i=ynsize; i<p-1; i++) {
// 			fill_Yn(i, prec);
// 		}
// 	}

// 	// should be done otherwise
// 	cln::cl_F one = cln::cl_float(1, cln::float_format(Digits));
// 	Number_T xf = x * one;
// 	//Number_T xf = x * cln::cl_float(1, prec);

// 	Number_T res;
// 	Number_T resbuf;
// 	Number_T factor = cln::expt(xf, p);
// 	int i = p;
// 	do {
// 		resbuf = res;
// 		if (i-p >= ynlength) {
// 			// make Yn longer
// 			make_Yn_longer(ynlength*2, prec);
// 		}
// 		res = res + factor / cln::expt(Integer_T(i),n+1) * Yn[p-2][i-p]; // should we check it? or rely on magic number? ...
// 		//res = res + factor / cln::expt(Integer_T(i),n+1) * (*it); // should we check it? or rely on magic number? ...
// 		factor = factor * xf;
// 		i++;
// 	} while (res != resbuf);
	
// 	return res;
//}


// helper function for S(n,p,x)
// Number_T S_projection(int n, int p, const Number_T& x, const cln::float_format_t& prec)
// {
// 	// [Kol] (5.3)
// 	if (cln::abs(cln::realpart(x)) > cln::cl_F("0.5")) {

// 		Number_T result = cln::expt(Integer_T(-1),p) * cln::expt(cln::log(x),n)
// 		                   * cln::expt(cln::log(1-x),p) / cln::factorial(n) / cln::factorial(p);

// 		for (int s=0; s<n; s++) {
// 			Number_T res2;
// 			for (int r=0; r<p; r++) {
// 				res2 = res2 + cln::expt(Integer_T(-1),r) * cln::expt(cln::log(1-x),r)
// 				              * S_do_sum(p-r,n-s,1-x,prec) / cln::factorial(r);
// 			}
// 			result = result + cln::expt(cln::log(x),s) * (S_num(n-s,p,1).to_cl_N() - res2) / cln::factorial(s);
// 		}

// 		return result;
// 	}
	
// 	return S_do_sum(n, p, x, prec);
// }


// helper function for S(n,p,x)
numeric S_num(int n, int p, const numeric& x)
{
// 	if (x == 1) {
// 		if (n == 1) {
// 		    // [Kol] (2.22) with (2.21)
// 			return cln::zeta(p+1);
// 		}

// 		if (p == 1) {
// 		    // [Kol] (2.22)
// 			return cln::zeta(n+1);
// 		}

// 		// [Kol] (9.1)
// 		Number_T result;
// 		for (int nu=0; nu<n; nu++) {
// 			for (int rho=0; rho<=p; rho++) {
// 				result = result + b_k(n-nu-1) * b_k(p-rho) * a_k(nu+rho+1)
// 				                  * cln::factorial(nu+rho+1) / cln::factorial(rho) / cln::factorial(nu+1);
// 			}
// 		}
// 		result = result * cln::expt(Integer_T(-1),n+p-1);

// 		return result;
// 	}
// 	else if (x == -1) {
// 		// [Kol] (2.22)
// 		if (p == 1) {
// 			return -(1-cln::expt(Integer_T(2),-n)) * cln::zeta(n+1);
// 		}
// //		throw std::runtime_error("don't know how to evaluate this function!");
// 	}

// 	// what is the desired float format?
// 	// first guess: default format
// 	cln::float_format_t prec = cln::default_float_format;
// 	const Number_T value = x.to_cl_N();
// 	// second guess: the argument's format
// 	if (!x.real().is_rational())
// 		prec = cln::float_format(cln::the<cln::cl_F>(cln::realpart(value)));
// 	else if (!x.imag().is_rational())
// 		prec = cln::float_format(cln::the<cln::cl_F>(cln::imagpart(value)));

// 	// [Kol] (5.3)
// 	if ((cln::realpart(value) < -0.5) || (n == 0) || ((cln::abs(value) <= 1) && (cln::abs(value) > 0.95))) {

// 		Number_T result = cln::expt(Integer_T(-1),p) * cln::expt(cln::log(value),n)
// 		                   * cln::expt(cln::log(1-value),p) / cln::factorial(n) / cln::factorial(p);

// 		for (int s=0; s<n; s++) {
// 			Number_T res2;
// 			for (int r=0; r<p; r++) {
// 				res2 = res2 + cln::expt(Integer_T(-1),r) * cln::expt(cln::log(1-value),r)
// 				              * S_num(p-r,n-s,1-value).to_cl_N() / cln::factorial(r);
// 			}
// 			result = result + cln::expt(cln::log(value),s) * (S_num(n-s,p,1).to_cl_N() - res2) / cln::factorial(s);
// 		}

// 		return result;
		
// 	}
// 	// [Kol] (5.12)
// 	if (cln::abs(value) > 1) {
		
// 		Number_T result;

// 		for (int s=0; s<p; s++) {
// 			for (int r=0; r<=s; r++) {
// 				result = result + cln::expt(Integer_T(-1),s) * cln::expt(cln::log(-value),r) * cln::factorial(n+s-r-1)
// 				                  / cln::factorial(r) / cln::factorial(s-r) / cln::factorial(n-1)
// 				                  * S_num(n+s-r,p-s,cln::recip(value)).to_cl_N();
// 			}
// 		}
// 		result = result * cln::expt(Integer_T(-1),n);

// 		Number_T res2;
// 		for (int r=0; r<n; r++) {
// 			res2 = res2 + cln::expt(cln::log(-value),r) * C(n-r,p) / cln::factorial(r);
// 		}
// 		res2 = res2 + cln::expt(cln::log(-value),n+p) / cln::factorial(n+p);

// 		result = result + cln::expt(Integer_T(-1),p) * res2;

// 		return result;
// 	}
// 	else {
// 		return S_projection(n, p, value, prec);
// 	}
}


} // end of anonymous namespace

*/

//////////////////////////////////////////////////////////////////////
//
// Nielsen's generalized polylogarithm  S(n,p,x)
//
// GiNaC function
//
//////////////////////////////////////////////////////////////////////

/*

static ex S_evalf(const ex& n, const ex& p, const ex& x, PyObject* parent)
{
	if (n.info(info_flags::posint) && p.info(info_flags::posint)) {
		if (is_a<numeric>(x)) {
			return S_num(ex_to<numeric>(n).to_int(), ex_to<numeric>(p).to_int(), ex_to<numeric>(x));
		} else {
			ex x_val = x.evalf(0, parent);
			if (is_a<numeric>(x_val)) {
				return S_num(ex_to<numeric>(n).to_int(), ex_to<numeric>(p).to_int(), ex_to<numeric>(x_val));
			}
		}
	}
	return S(n, p, x).hold();
}


static ex S_eval(const ex& n, const ex& p, const ex& x)
{
	if (n.info(info_flags::posint) && p.info(info_flags::posint)) {
		if (x == 0) {
			return _ex0;
		}
		if (x == 1) {
			lst m(n+1);
			for (int i=ex_to<numeric>(p).to_int()-1; i>0; i--) {
				m.append(1);
			}
			return zeta(m);
		}
		if (p == 1) {
			return Li(n+1, x);
		}
		if (x.info(info_flags::numeric) && (!x.info(info_flags::crational))) {
			return S_num(ex_to<numeric>(n).to_int(), ex_to<numeric>(p).to_int(), ex_to<numeric>(x));
		}
	}
	if (n.is_zero()) {
		// [Kol] (5.3)
		return pow(-log(1-x), p) / factorial(p);
	}
	return S(n, p, x).hold();
}


static ex S_series(const ex& n, const ex& p, const ex& x, const relational& rel, int order, unsigned options)
{
	if (p == _ex1) {
		return Li(n+1, x).series(rel, order, options);
	}

	const ex x_pt = x.subs(rel, subs_options::no_pattern);
	if (n.info(info_flags::posint) && p.info(info_flags::posint) && x_pt.info(info_flags::numeric)) {
		// First special case: x==0 (derivatives have poles)
		if (x_pt.is_zero()) {
			const symbol s;
			ex ser;
			// manually construct the primitive expansion
			// subsum = Euler-Zagier-Sum is needed
			// dirty hack (slow ...) calculation of subsum:
			std::vector<ex> presubsum, subsum;
			subsum.push_back(0);
			for (int i=1; i<order-1; ++i) {
				subsum.push_back(subsum[i-1] + numeric(1, i));
			}
			for (int depth=2; depth<p; ++depth) {
				presubsum = subsum;
				for (int i=1; i<order-1; ++i) {
					subsum[i] = subsum[i-1] + numeric(1, i) * presubsum[i-1];
				}
			}
				
			for (int i=1; i<order; ++i) {
				ser += pow(s,i) / pow(numeric(i), n+1) * subsum[i-1];
			}
			// substitute the argument's series expansion
			ser = ser.subs(s==x.series(rel, order), subs_options::no_pattern);
			// maybe that was terminating, so add a proper order term
			epvector nseq;
			nseq.push_back(expair(Order(_ex1), order));
			ser += pseries(rel, nseq);
			// reexpanding it will collapse the series again
			return ser.series(rel, order);
		}
		// TODO special cases: x==1 (branch point) and x real, >=1 (branch cut)
		throw std::runtime_error("S_series: don't know how to do the series expansion at this point!");
	}
	// all other cases should be safe, by now:
	throw do_taylor();  // caught by function::series()
}


static ex S_deriv(const ex& n, const ex& p, const ex& x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param < 3);
	if (deriv_param < 2) {
		return _ex0;
	}
	if (n > 0) {
		return S(n-1, p, x) / x;
	} else {
		return S(n, p-1, x) / (1-x);
	}
}


static void S_print_latex(const ex& n, const ex& p, const ex& x, const print_context& c)
{
	c.s << "{\\rm S}_{";
	n.print(c);
	c.s << ",";
	p.print(c);
	c.s << "}(";
	x.print(c);
	c.s << ")";
}


REGISTER_FUNCTION(S,
                  evalf_func(S_evalf).
                  eval_func(S_eval).
                  series_func(S_series).
                  derivative_func(S_deriv).
                  print_func<print_latex>(S_print_latex).
                  do_not_evalf_params());

 */

//////////////////////////////////////////////////////////////////////
//
// Harmonic polylogarithm  H(m,x)
//
// helper functions
//
//////////////////////////////////////////////////////////////////////

/*
 * 
// anonymous namespace for helper functions
namespace {

	
// regulates the pole (used by 1/x-transformation)
symbol H_polesign("IMSIGN");


// convert parameters from H to Li representation
// parameters are expected to be in expanded form, i.e. only 0, 1 and -1
// returns true if some parameters are negative
bool convert_parameter_H_to_Li(const lst& l, lst& m, lst& s, ex& pf)
{
	// expand parameter list
	lst mexp;
	for (lst::const_iterator it = l.begin(); it != l.end(); it++) {
		if (*it > 1) {
			for (ex count=*it-1; count > 0; count--) {
				mexp.append(0);
			}
			mexp.append(1);
		} else if (*it < -1) {
			for (ex count=*it+1; count < 0; count++) {
				mexp.append(0);
			}
			mexp.append(-1);
		} else {
			mexp.append(*it);
		}
	}
	
	ex signum = 1;
	pf = 1;
	bool has_negative_parameters = false;
	ex acc = 1;
	for (lst::const_iterator it = mexp.begin(); it != mexp.end(); it++) {
		if (*it == 0) {
			acc++;
			continue;
		}
		if (*it > 0) {
			m.append((*it+acc-1) * signum);
		} else {
			m.append((*it-acc+1) * signum);
		}
		acc = 1;
		signum = *it;
		pf *= *it;
		if (pf < 0) {
			has_negative_parameters = true;
		}
	}
	if (has_negative_parameters) {
		for (int i=0; i<m.nops(); i++) {
			if (m.op(i) < 0) {
				m.let_op(i) = -m.op(i);
				s.append(-1);
			} else {
				s.append(1);
			}
		}
	}
	
	return has_negative_parameters;
}


// recursively transforms H to corresponding multiple polylogarithms
struct map_trafo_H_convert_to_Li : public map_function
{
	ex operator()(const ex& e)
	{
		if (is_a<add>(e) || is_a<mul>(e)) {
			return e.map(*this);
		}
		if (is_a<function>(e)) {
			std::string name = ex_to<function>(e).get_name();
			if (name == "H") {
				lst parameter;
				if (is_a<lst>(e.op(0))) {
						parameter = ex_to<lst>(e.op(0));
				} else {
					parameter = lst(e.op(0));
				}
				ex arg = e.op(1);

				lst m;
				lst s;
				ex pf;
				if (convert_parameter_H_to_Li(parameter, m, s, pf)) {
					s.let_op(0) = s.op(0) * arg;
					return pf * Li(m, s).hold();
				} else {
					for (int i=0; i<m.nops(); i++) {
						s.append(1);
					}
					s.let_op(0) = s.op(0) * arg;
					return Li(m, s).hold();
				}
			}
		}
		return e;
	}
};


// recursively transforms H to corresponding zetas
struct map_trafo_H_convert_to_zeta : public map_function
{
	ex operator()(const ex& e)
	{
		if (is_a<add>(e) || is_a<mul>(e)) {
			return e.map(*this);
		}
		if (is_a<function>(e)) {
			std::string name = ex_to<function>(e).get_name();
			if (name == "H") {
				lst parameter;
				if (is_a<lst>(e.op(0))) {
						parameter = ex_to<lst>(e.op(0));
				} else {
					parameter = lst(e.op(0));
				}

				lst m;
				lst s;
				ex pf;
				if (convert_parameter_H_to_Li(parameter, m, s, pf)) {
					return pf * zeta(m, s);
				} else {
					return zeta(m);
				}
			}
		}
		return e;
	}
};


// remove trailing zeros from H-parameters
struct map_trafo_H_reduce_trailing_zeros : public map_function
{
	ex operator()(const ex& e)
	{
		if (is_a<add>(e) || is_a<mul>(e)) {
			return e.map(*this);
		}
		if (is_a<function>(e)) {
			std::string name = ex_to<function>(e).get_name();
			if (name == "H") {
				lst parameter;
				if (is_a<lst>(e.op(0))) {
					parameter = ex_to<lst>(e.op(0));
				} else {
					parameter = lst(e.op(0));
				}
				ex arg = e.op(1);
				if (parameter.op(parameter.nops()-1) == 0) {
					
					//
					if (parameter.nops() == 1) {
						return log(arg);
					}
					
					//
					lst::const_iterator it = parameter.begin();
					while ((it != parameter.end()) && (*it == 0)) {
						it++;
					}
					if (it == parameter.end()) {
						return pow(log(arg),parameter.nops()) / factorial(parameter.nops());
					}
					
					//
					parameter.remove_last();
					int lastentry = parameter.nops();
					while ((lastentry > 0) && (parameter[lastentry-1] == 0)) {
						lastentry--;
					}
					
					//
					ex result = log(arg) * H(parameter,arg).hold();
					ex acc = 0;
					for (ex i=0; i<lastentry; i++) {
						if (parameter[i] > 0) {
							parameter[i]++;
							result -= (acc + parameter[i]-1) * H(parameter, arg).hold();
							parameter[i]--;
							acc = 0;
						} else if (parameter[i] < 0) {
							parameter[i]--;
							result -= (acc + abs(parameter[i]+1)) * H(parameter, arg).hold();
							parameter[i]++;
							acc = 0;
						} else {
							acc++;
						}
					}
					
					if (lastentry < parameter.nops()) {
						result = result / (parameter.nops()-lastentry+1);
						return result.map(*this);
					} else {
						return result;
					}
				}
			}
		}
		return e;
	}
};


// returns an expression with zeta functions corresponding to the parameter list for H
ex convert_H_to_zeta(const lst& m)
{
	symbol xtemp("xtemp");
	map_trafo_H_reduce_trailing_zeros filter;
	map_trafo_H_convert_to_zeta filter2;
	return filter2(filter(H(m, xtemp).hold())).subs(xtemp == 1);
}


// convert signs form Li to H representation
lst convert_parameter_Li_to_H(const lst& m, const lst& x, ex& pf)
{
 	lst res;
	lst::const_iterator itm = m.begin();
	lst::const_iterator itx = ++x.begin();
	int signum = 1;
	pf = _ex1;
	res.append(*itm);
	itm++;
	while (itx != x.end()) {
		signum *= (*itx > 0) ? 1 : -1;
		pf *= signum;
		res.append((*itm) * signum);
		itm++;
		itx++;
	}
	return res;
}


// multiplies an one-dimensional H with another H
// [ReV] (18)
ex trafo_H_mult(const ex& h1, const ex& h2)
{
	ex res;
	ex hshort;
	lst hlong;
	ex h1nops = h1.op(0).nops();
	ex h2nops = h2.op(0).nops();
	if (h1nops > 1) {
		hshort = h2.op(0).op(0);
		hlong = ex_to<lst>(h1.op(0));
	} else {
		hshort = h1.op(0).op(0);
		if (h2nops > 1) {
			hlong = ex_to<lst>(h2.op(0));
		} else {
			hlong = h2.op(0).op(0);
		}
	}
	for (int i=0; i<=hlong.nops(); i++) {
		lst newparameter;
		int j=0;
		for (; j<i; j++) {
			newparameter.append(hlong[j]);
		}
		newparameter.append(hshort);
		for (; j<hlong.nops(); j++) {
			newparameter.append(hlong[j]);
		}
		res += H(newparameter, h1.op(1)).hold();
	}
	return res;
}


// applies trafo_H_mult recursively on expressions
struct map_trafo_H_mult : public map_function
{
	ex operator()(const ex& e)
	{
		if (is_a<add>(e)) {
			return e.map(*this);
		}

		if (is_a<mul>(e)) {

			ex result = 1;
			ex firstH;
			lst Hlst;
			for (int pos=0; pos<e.nops(); pos++) {
				if (is_a<power>(e.op(pos)) && is_a<function>(e.op(pos).op(0))) {
					std::string name = ex_to<function>(e.op(pos).op(0)).get_name();
					if (name == "H") {
						for (ex i=0; i<e.op(pos).op(1); i++) {
							Hlst.append(e.op(pos).op(0));
						}
						continue;
					}
				} else if (is_a<function>(e.op(pos))) {
					std::string name = ex_to<function>(e.op(pos)).get_name();
					if (name == "H") {
						if (e.op(pos).op(0).nops() > 1) {
							firstH = e.op(pos);
						} else {
							Hlst.append(e.op(pos));
						}
						continue;
					}
				}
				result *= e.op(pos);
			}
			if (firstH == 0) {
				if (Hlst.nops() > 0) {
					firstH = Hlst[Hlst.nops()-1];
					Hlst.remove_last();
				} else {
					return e;
				}
			}

			if (Hlst.nops() > 0) {
				ex buffer = trafo_H_mult(firstH, Hlst.op(0));
				result *= buffer;
				for (int i=1; i<Hlst.nops(); i++) {
					result *= Hlst.op(i);
				}
				result = result.expand();
				map_trafo_H_mult recursion;
				return recursion(result);
			} else {
				return e;
			}

		}
		return e;
	}
};


// do integration [ReV] (55)
// put parameter 0 in front of existing parameters
ex trafo_H_1tx_prepend_zero(const ex& e, const ex& arg)
{
	ex h;
	std::string name;
	if (is_a<function>(e)) {
		name = ex_to<function>(e).get_name();
	}
	if (name == "H") {
		h = e;
	} else {
		for (int i=0; i<e.nops(); i++) {
			if (is_a<function>(e.op(i))) {
				std::string name = ex_to<function>(e.op(i)).get_name();
				if (name == "H") {
					h = e.op(i);
				}
			}
		}
	}
	if (h != 0) {
		lst newparameter = ex_to<lst>(h.op(0));
		newparameter.prepend(0);
		ex addzeta = convert_H_to_zeta(newparameter);
		return e.subs(h == (addzeta-H(newparameter, h.op(1)).hold())).expand();
	} else {
		return e * (-H(lst(0),1/arg).hold());
	}
}


// do integration [ReV] (49)
// put parameter 1 in front of existing parameters
ex trafo_H_prepend_one(const ex& e, const ex& arg)
{
	ex h;
	std::string name;
	if (is_a<function>(e)) {
		name = ex_to<function>(e).get_name();
	}
	if (name == "H") {
		h = e;
	} else {
		for (int i=0; i<e.nops(); i++) {
			if (is_a<function>(e.op(i))) {
				std::string name = ex_to<function>(e.op(i)).get_name();
				if (name == "H") {
					h = e.op(i);
				}
			}
		}
	}
	if (h != 0) {
		lst newparameter = ex_to<lst>(h.op(0));
		newparameter.prepend(1);
		return e.subs(h == H(newparameter, h.op(1)).hold());
	} else {
		return e * H(lst(1),1-arg).hold();
	}
}


// do integration [ReV] (55)
// put parameter -1 in front of existing parameters
ex trafo_H_1tx_prepend_minusone(const ex& e, const ex& arg)
{
	ex h;
	std::string name;
	if (is_a<function>(e)) {
		name = ex_to<function>(e).get_name();
	}
	if (name == "H") {
		h = e;
	} else {
		for (int i=0; i<e.nops(); i++) {
			if (is_a<function>(e.op(i))) {
				std::string name = ex_to<function>(e.op(i)).get_name();
				if (name == "H") {
					h = e.op(i);
				}
			}
		}
	}
	if (h != 0) {
		lst newparameter = ex_to<lst>(h.op(0));
		newparameter.prepend(-1);
		ex addzeta = convert_H_to_zeta(newparameter);
		return e.subs(h == (addzeta-H(newparameter, h.op(1)).hold())).expand();
	} else {
		ex addzeta = convert_H_to_zeta(lst(-1));
		return (e * (addzeta - H(lst(-1),1/arg).hold())).expand();
	}
}


// do integration [ReV] (55)
// put parameter -1 in front of existing parameters
ex trafo_H_1mxt1px_prepend_minusone(const ex& e, const ex& arg)
{
	ex h;
	std::string name;
	if (is_a<function>(e)) {
		name = ex_to<function>(e).get_name();
	}
	if (name == "H") {
		h = e;
	} else {
		for (int i=0; i<e.nops(); i++) {
			if (is_a<function>(e.op(i))) {
				std::string name = ex_to<function>(e.op(i)).get_name();
				if (name == "H") {
					h = e.op(i);
				}
			}
		}
	}
	if (h != 0) {
		lst newparameter = ex_to<lst>(h.op(0));
		newparameter.prepend(-1);
		return e.subs(h == H(newparameter, h.op(1)).hold()).expand();
	} else {
		return (e * H(lst(-1),(1-arg)/(1+arg)).hold()).expand();
	}
}


// do integration [ReV] (55)
// put parameter 1 in front of existing parameters
ex trafo_H_1mxt1px_prepend_one(const ex& e, const ex& arg)
{
	ex h;
	std::string name;
	if (is_a<function>(e)) {
		name = ex_to<function>(e).get_name();
	}
	if (name == "H") {
		h = e;
	} else {
		for (int i=0; i<e.nops(); i++) {
			if (is_a<function>(e.op(i))) {
				std::string name = ex_to<function>(e.op(i)).get_name();
				if (name == "H") {
					h = e.op(i);
				}
			}
		}
	}
	if (h != 0) {
		lst newparameter = ex_to<lst>(h.op(0));
		newparameter.prepend(1);
		return e.subs(h == H(newparameter, h.op(1)).hold()).expand();
	} else {
		return (e * H(lst(1),(1-arg)/(1+arg)).hold()).expand();
	}
}


// do x -> 1-x transformation
struct map_trafo_H_1mx : public map_function
{
	ex operator()(const ex& e)
	{
		if (is_a<add>(e) || is_a<mul>(e)) {
			return e.map(*this);
		}
		
		if (is_a<function>(e)) {
			std::string name = ex_to<function>(e).get_name();
			if (name == "H") {

				lst parameter = ex_to<lst>(e.op(0));
				ex arg = e.op(1);

				// special cases if all parameters are either 0, 1 or -1
				bool allthesame = true;
				if (parameter.op(0) == 0) {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != 0) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						lst newparameter;
						for (int i=parameter.nops(); i>0; i--) {
							newparameter.append(1);
						}
						return pow(-1, parameter.nops()) * H(newparameter, 1-arg).hold();
					}
				} else if (parameter.op(0) == -1) {
					throw std::runtime_error("map_trafo_H_1mx: cannot handle weights equal -1!");
				} else {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != 1) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						lst newparameter;
						for (int i=parameter.nops(); i>0; i--) {
							newparameter.append(0);
						}
						return pow(-1, parameter.nops()) * H(newparameter, 1-arg).hold();
					}
				}

				lst newparameter = parameter;
				newparameter.remove_first();

				if (parameter.op(0) == 0) {

					// leading zero
					ex res = convert_H_to_zeta(parameter);
					//ex res = convert_from_RV(parameter, 1).subs(H(wild(1),wild(2))==zeta(wild(1)));
					map_trafo_H_1mx recursion;
					ex buffer = recursion(H(newparameter, arg).hold());
					if (is_a<add>(buffer)) {
						for (int i=0; i<buffer.nops(); i++) {
							res -= trafo_H_prepend_one(buffer.op(i), arg);
						}
					} else {
						res -= trafo_H_prepend_one(buffer, arg);
					}
					return res;

				} else {

					// leading one
					map_trafo_H_1mx recursion;
					map_trafo_H_mult unify;
					ex res = H(lst(1), arg).hold() * H(newparameter, arg).hold();
					int firstzero = 0;
					while (parameter.op(firstzero) == 1) {
						firstzero++;
					}
					for (int i=firstzero-1; i<parameter.nops()-1; i++) {
						lst newparameter;
						int j=0;
						for (; j<=i; j++) {
							newparameter.append(parameter[j+1]);
						}
						newparameter.append(1);
						for (; j<parameter.nops()-1; j++) {
							newparameter.append(parameter[j+1]);
						}
						res -= H(newparameter, arg).hold();
					}
					res = recursion(res).expand() / firstzero;
					return unify(res);
				}
			}
		}
		return e;
	}
};


// do x -> 1/x transformation
struct map_trafo_H_1overx : public map_function
{
	ex operator()(const ex& e)
	{
		if (is_a<add>(e) || is_a<mul>(e)) {
			return e.map(*this);
		}

		if (is_a<function>(e)) {
			std::string name = ex_to<function>(e).get_name();
			if (name == "H") {

				lst parameter = ex_to<lst>(e.op(0));
				ex arg = e.op(1);

				// special cases if all parameters are either 0, 1 or -1
				bool allthesame = true;
				if (parameter.op(0) == 0) {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != 0) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						return pow(-1, parameter.nops()) * H(parameter, 1/arg).hold();
					}
				} else if (parameter.op(0) == -1) {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != -1) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						map_trafo_H_mult unify;
						return unify((pow(H(lst(-1),1/arg).hold() - H(lst(0),1/arg).hold(), parameter.nops())
						       / factorial(parameter.nops())).expand());
					}
				} else {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != 1) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						map_trafo_H_mult unify;
						return unify((pow(H(lst(1),1/arg).hold() + H(lst(0),1/arg).hold() + H_polesign, parameter.nops())
						       / factorial(parameter.nops())).expand());
					}
				}

				lst newparameter = parameter;
				newparameter.remove_first();

				if (parameter.op(0) == 0) {
					
					// leading zero
					ex res = convert_H_to_zeta(parameter);
					map_trafo_H_1overx recursion;
					ex buffer = recursion(H(newparameter, arg).hold());
					if (is_a<add>(buffer)) {
						for (int i=0; i<buffer.nops(); i++) {
							res += trafo_H_1tx_prepend_zero(buffer.op(i), arg);
						}
					} else {
						res += trafo_H_1tx_prepend_zero(buffer, arg);
					}
					return res;

				} else if (parameter.op(0) == -1) {

					// leading negative one
					ex res = convert_H_to_zeta(parameter);
					map_trafo_H_1overx recursion;
					ex buffer = recursion(H(newparameter, arg).hold());
					if (is_a<add>(buffer)) {
						for (int i=0; i<buffer.nops(); i++) {
							res += trafo_H_1tx_prepend_zero(buffer.op(i), arg) - trafo_H_1tx_prepend_minusone(buffer.op(i), arg);
						}
					} else {
						res += trafo_H_1tx_prepend_zero(buffer, arg) - trafo_H_1tx_prepend_minusone(buffer, arg);
					}
					return res;

				} else {

					// leading one
					map_trafo_H_1overx recursion;
					map_trafo_H_mult unify;
					ex res = H(lst(1), arg).hold() * H(newparameter, arg).hold();
					int firstzero = 0;
					while (parameter.op(firstzero) == 1) {
						firstzero++;
					}
					for (int i=firstzero-1; i<parameter.nops()-1; i++) {
						lst newparameter;
						int j=0;
						for (; j<=i; j++) {
							newparameter.append(parameter[j+1]);
						}
						newparameter.append(1);
						for (; j<parameter.nops()-1; j++) {
							newparameter.append(parameter[j+1]);
						}
						res -= H(newparameter, arg).hold();
					}
					res = recursion(res).expand() / firstzero;
					return unify(res);

				}

			}
		}
		return e;
	}
};


// do x -> (1-x)/(1+x) transformation
struct map_trafo_H_1mxt1px : public map_function
{
	ex operator()(const ex& e)
	{
		if (is_a<add>(e) || is_a<mul>(e)) {
			return e.map(*this);
		}

		if (is_a<function>(e)) {
			std::string name = ex_to<function>(e).get_name();
			if (name == "H") {

				lst parameter = ex_to<lst>(e.op(0));
				ex arg = e.op(1);

				// special cases if all parameters are either 0, 1 or -1
				bool allthesame = true;
				if (parameter.op(0) == 0) {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != 0) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						map_trafo_H_mult unify;
						return unify((pow(-H(lst(1),(1-arg)/(1+arg)).hold() - H(lst(-1),(1-arg)/(1+arg)).hold(), parameter.nops())
						       / factorial(parameter.nops())).expand());
					}
				} else if (parameter.op(0) == -1) {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != -1) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						map_trafo_H_mult unify;
						return unify((pow(log(2) - H(lst(-1),(1-arg)/(1+arg)).hold(), parameter.nops())
						       / factorial(parameter.nops())).expand());
					}
				} else {
					for (int i=1; i<parameter.nops(); i++) {
						if (parameter.op(i) != 1) {
							allthesame = false;
							break;
						}
					}
					if (allthesame) {
						map_trafo_H_mult unify;
						return unify((pow(-log(2) - H(lst(0),(1-arg)/(1+arg)).hold() + H(lst(-1),(1-arg)/(1+arg)).hold(), parameter.nops())
						       / factorial(parameter.nops())).expand());
					}
				}

				lst newparameter = parameter;
				newparameter.remove_first();

				if (parameter.op(0) == 0) {

					// leading zero
					ex res = convert_H_to_zeta(parameter);
					map_trafo_H_1mxt1px recursion;
					ex buffer = recursion(H(newparameter, arg).hold());
					if (is_a<add>(buffer)) {
						for (int i=0; i<buffer.nops(); i++) {
							res -= trafo_H_1mxt1px_prepend_one(buffer.op(i), arg) + trafo_H_1mxt1px_prepend_minusone(buffer.op(i), arg);
						}
					} else {
						res -= trafo_H_1mxt1px_prepend_one(buffer, arg) + trafo_H_1mxt1px_prepend_minusone(buffer, arg);
					}
					return res;

				} else if (parameter.op(0) == -1) {

					// leading negative one
					ex res = convert_H_to_zeta(parameter);
					map_trafo_H_1mxt1px recursion;
					ex buffer = recursion(H(newparameter, arg).hold());
					if (is_a<add>(buffer)) {
						for (int i=0; i<buffer.nops(); i++) {
							res -= trafo_H_1mxt1px_prepend_minusone(buffer.op(i), arg);
						}
					} else {
						res -= trafo_H_1mxt1px_prepend_minusone(buffer, arg);
					}
					return res;

				} else {

					// leading one
					map_trafo_H_1mxt1px recursion;
					map_trafo_H_mult unify;
					ex res = H(lst(1), arg).hold() * H(newparameter, arg).hold();
					int firstzero = 0;
					while (parameter.op(firstzero) == 1) {
						firstzero++;
					}
					for (int i=firstzero-1; i<parameter.nops()-1; i++) {
						lst newparameter;
						int j=0;
						for (; j<=i; j++) {
							newparameter.append(parameter[j+1]);
						}
						newparameter.append(1);
						for (; j<parameter.nops()-1; j++) {
							newparameter.append(parameter[j+1]);
						}
						res -= H(newparameter, arg).hold();
					}
					res = recursion(res).expand() / firstzero;
					return unify(res);

				}

			}
		}
		return e;
	}
};


// do the actual summation.
numeric H_do_sum(const std::vector<int>& m, const numeric& x)
{
// 	const int j = m.size();

// 	std::vector<Number_T> t(j);

// 	cln::cl_F one = cln::cl_float(1, cln::float_format(Digits));
// 	Number_T factor = cln::expt(x, j) * one;
// 	Number_T t0buf;
// 	int q = 0;
// 	do {
// 		t0buf = t[0];
// 		q++;
// 		t[j-1] = t[j-1] + 1 / cln::expt(Integer_T(q),m[j-1]);
// 		for (int k=j-2; k>=1; k--) {
// 			t[k] = t[k] + t[k+1] / cln::expt(Integer_T(q+j-1-k), m[k]);
// 		}
// 		t[0] = t[0] + t[1] * factor / cln::expt(Integer_T(q+j-1), m[0]);
// 		factor = factor * x;
// 	} while (t[0] != t0buf);

  return 0;
}


} // end of anonymous namespace

 */

//////////////////////////////////////////////////////////////////////
//
// Harmonic polylogarithm  H(m,x)
//
// GiNaC function
//
//////////////////////////////////////////////////////////////////////

/*

static ex H_evalf(const ex& x1, const ex& x2, PyObject* parent)
{
// 	if (is_a<lst>(x1)) {
		
// 		Number_T x;
// 		if (is_a<numeric>(x2)) {
// 			x = ex_to<numeric>(x2).to_cl_N();
// 		} else {
// 			ex x2_val = x2.evalf(0, parent);
// 			if (is_a<numeric>(x2_val)) {
// 				x = ex_to<numeric>(x2_val).to_cl_N();
// 			}
// 		}

// 		for (int i=0; i<x1.nops(); i++) {
// 			if (!x1.op(i).info(info_flags::integer)) {
// 				return H(x1, x2).hold();
// 			}
// 		}
// 		if (x1.nops() < 1) {
// 			return H(x1, x2).hold();
// 		}

// 		const lst& morg = ex_to<lst>(x1);
// 		// remove trailing zeros ...
// 		if (*(--morg.end()) == 0) {
// 			symbol xtemp("xtemp");
// 			map_trafo_H_reduce_trailing_zeros filter;
// 			return filter(H(x1, xtemp).hold()).subs(xtemp==x2).evalf(0, parent);
// 		}
// 		// ... and expand parameter notation
// 		bool has_minus_one = false;
// 		lst m;
// 		for (lst::const_iterator it = morg.begin(); it != morg.end(); it++) {
// 			if (*it > 1) {
// 				for (ex count=*it-1; count > 0; count--) {
// 					m.append(0);
// 				}
// 				m.append(1);
// 			} else if (*it <= -1) {
// 				for (ex count=*it+1; count < 0; count++) {
// 					m.append(0);
// 				}
// 				m.append(-1);
// 				has_minus_one = true;
// 			} else {
// 				m.append(*it);
// 			}
// 		}

// 		// do summation
// 		if (cln::abs(x) < 0.95) {
// 			lst m_lst;
// 			lst s_lst;
// 			ex pf;
// 			if (convert_parameter_H_to_Li(m, m_lst, s_lst, pf)) {
// 				// negative parameters -> s_lst is filled
// 				std::vector<int> m_int;
// 				std::vector<Number_T> x_cln;
// 				for (lst::const_iterator it_int = m_lst.begin(), it_cln = s_lst.begin(); 
// 				     it_int != m_lst.end(); it_int++, it_cln++) {
// 					m_int.push_back(ex_to<numeric>(*it_int).to_int());
// 					x_cln.push_back(ex_to<numeric>(*it_cln).to_cl_N());
// 				}
// 				x_cln.front() = x_cln.front() * x;
// 				return pf * numeric(multipleLi_do_sum(m_int, x_cln));
// 			} else {
// 				// only positive parameters
// 				//TODO
// 				if (m_lst.nops() == 1) {
// 					return Li(m_lst.op(0), x2).evalf(0, parent);
// 				}
// 				std::vector<int> m_int;
// 				for (lst::const_iterator it = m_lst.begin(); it != m_lst.end(); it++) {
// 					m_int.push_back(ex_to<numeric>(*it).to_int());
// 				}
// 				return numeric(H_do_sum(m_int, x));
// 			}
// 		}

// 		symbol xtemp("xtemp");
// 		ex res = 1;	
		
// 		// ensure that the realpart of the argument is positive
// 		if (cln::realpart(x) < 0) {
// 			x = -x;
// 			for (int i=0; i<m.nops(); i++) {
// 				if (m.op(i) != 0) {
// 					m.let_op(i) = -m.op(i);
// 					res *= -1;
// 				}
// 			}
// 		}

// 		// x -> 1/x
// 		if (cln::abs(x) >= 2.0) {
// 			map_trafo_H_1overx trafo;
// 			res *= trafo(H(m, xtemp));
// 			if (cln::imagpart(x) <= 0) {
// 				res = res.subs(H_polesign == -I*Pi);
// 			} else {
// 				res = res.subs(H_polesign == I*Pi);
// 			}
// 			return res.subs(xtemp == numeric(x)).evalf(0, parent);
// 		}
		
// 		// check transformations for 0.95 <= |x| < 2.0
		
// 		// |(1-x)/(1+x)| < 0.9 -> circular area with center=9.53+0i and radius=9.47
// 		if (cln::abs(x-9.53) <= 9.47) {
// 			// x -> (1-x)/(1+x)
// 			map_trafo_H_1mxt1px trafo;
// 			res *= trafo(H(m, xtemp));
// 		} else {
// 			// x -> 1-x
// 			if (has_minus_one) {
// 				map_trafo_H_convert_to_Li filter;
// 				return filter(H(m, numeric(x)).hold()).evalf(0, parent);
// 			}
// 			map_trafo_H_1mx trafo;
// 			res *= trafo(H(m, xtemp));
// 		}

// 		return res.subs(xtemp == numeric(x)).evalf(0, parent);
// 	}

// 	return H(x1,x2).hold();
  return 0;
}


static ex H_eval(const ex& m_, const ex& x)
{
	lst m;
	if (is_a<lst>(m_)) {
		m = ex_to<lst>(m_);
	} else {
		m = lst(m_);
	}
	if (m.nops() == 0) {
		return _ex1;
	}
	ex pos1;
	ex pos2;
	ex n;
	ex p;
	int step = 0;
	if (*m.begin() > _ex1) {
		step++;
		pos1 = _ex0;
		pos2 = _ex1;
		n = *m.begin()-1;
		p = _ex1;
	} else if (*m.begin() < _ex_1) {
		step++;
		pos1 = _ex0;
		pos2 = _ex_1;
		n = -*m.begin()-1;
		p = _ex1;
	} else if (*m.begin() == _ex0) {
		pos1 = _ex0;
		n = _ex1;
	} else {
		pos1 = *m.begin();
		p = _ex1;
	}
	for (auto it = ++m.begin(); it != m.end(); it++) {
		if ((*it).info(info_flags::integer)) {
			if (step == 0) {
				if (*it > _ex1) {
					if (pos1 == _ex0) {
						step = 1;
						pos2 = _ex1;
						n += *it-1;
						p = _ex1;
					} else {
						step = 2;
					}
				} else if (*it < _ex_1) {
					if (pos1 == _ex0) {
						step = 1;
						pos2 = _ex_1;
						n += -*it-1;
						p = _ex1;
					} else {
						step = 2;
					}
				} else {
					if (*it != pos1) {
						step = 1;
						pos2 = *it;
					}
					if (*it == _ex0) {
						n++;
					} else {
						p++;
					}
				}
			} else if (step == 1) {
				if (*it != pos2) {
					step = 2;
				} else {
					if (*it == _ex0) {
						n++;
					} else {
						p++;
					}
				}
			}
		} else {
			// if some m_i is not an integer
			return H(m_, x).hold();
		}
	}
	if ((x == _ex1) && (*(--m.end()) != _ex0)) {
		return convert_H_to_zeta(m);
	}
	if (step == 0) {
		if (pos1 == _ex0) {
			// all zero
			if (x == _ex0) {
				return H(m_, x).hold();
			}
			return pow(log(x), m.nops()) / factorial(m.nops());
		} else {
			// all (minus) one
			return pow(-pos1*log(1-pos1*x), m.nops()) / factorial(m.nops());
		}
	} else if ((step == 1) && (pos1 == _ex0)){
		// convertible to S
		if (pos2 == _ex1) {
			return S(n, p, x);
		} else {
			return pow(-1, p) * S(n, p, -x);
		}
	}
	if (x == _ex0) {
		return _ex0;
	}
	if (x.info(info_flags::numeric) && (!x.info(info_flags::crational))) {
		return H(m_, x).evalf();
	}
	return H(m_, x).hold();
}


static ex H_series(const ex& m, const ex& x, const relational& rel, int order, unsigned options)
{
	epvector seq;
	seq.push_back(expair(H(m, x), 0));
	return pseries(rel, seq);
}


static ex H_deriv(const ex& m_, const ex& x, unsigned deriv_param)
{
	GINAC_ASSERT(deriv_param < 2);
	if (deriv_param == 0) {
		return _ex0;
	}
	lst m;
	if (is_a<lst>(m_)) {
		m = ex_to<lst>(m_);
	} else {
		m = lst(m_);
	}
	ex mb = *m.begin();
	if (mb > _ex1) {
		m[0]--;
		return H(m, x) / x;
	}
	if (mb < _ex_1) {
		m[0]++;
		return H(m, x) / x;
	}
	m.remove_first();
	if (mb == _ex1) {
		return 1/(1-x) * H(m, x);
	} else if (mb == _ex_1) {
		return 1/(1+x) * H(m, x);
	} else {
		return H(m, x) / x;
	}
}


static void H_print_latex(const ex& m_, const ex& x, const print_context& c)
{
	lst m;
	if (is_a<lst>(m_)) {
		m = ex_to<lst>(m_);
	} else {
		m = lst(m_);
	}
	c.s << "{\\rm H}_{";
	auto itm = m.begin();
	(*itm).print(c);
	itm++;
	for (; itm != m.end(); itm++) {
		c.s << ",";
		(*itm).print(c);
	}
	c.s << "}(";
	x.print(c);
	c.s << ")";
}


REGISTER_FUNCTION(H,
                  evalf_func(H_evalf).
                  eval_func(H_eval).
                  series_func(H_series).
                  derivative_func(H_deriv).
                  print_func<print_latex>(H_print_latex).
                  do_not_evalf_params());


// takes a parameter list for H and returns an expression with corresponding multiple polylogarithms
ex convert_H_to_Li(const ex& m, const ex& x)
{
	map_trafo_H_reduce_trailing_zeros filter;
	map_trafo_H_convert_to_Li filter2;
	if (is_a<lst>(m)) {
		return filter2(filter(H(m, x).hold()));
	} else {
		return filter2(filter(H(lst(m), x).hold()));
	}
}
*/

//////////////////////////////////////////////////////////////////////
//
// Multiple zeta values  zeta(x) and zeta(x,s)
//
// helper functions
//
//////////////////////////////////////////////////////////////////////

/*

// anonymous namespace for helper functions
namespace {


// parameters and data for [Cra] algorithm
//const Number_T lambda = Number_T("319/320");
//const numeric lambda = numeric(319)/numeric(320);
//int L1;
//int L2;
//std::vector<std::vector<numeric> > f_kj;
//std::vector<numeric> crB;
//std::vector<std::vector<numeric> > crG;
//std::vector<numeric> crX;


void halfcyclic_convolute(const std::vector<numeric>& a, const std::vector<numeric>& b, std::vector<numeric>& c)
{
	const int size = a.size();
	for (int n=0; n<size; n++) {
		c[n] = 0;
		for (int m=0; m<=n; m++) {
			c[n] = c[n] + a[m]*b[n-m];
		}
	}
}


// [Cra] section 4
void initcX(const std::vector<int>& s)
{
// 	const int k = s.size();

// 	crX.clear();
// 	crG.clear();
// 	crB.clear();

// 	for (int i=0; i<=L2; i++) {
// 		crB.push_back(bernoulli(i).to_cl_N() / cln::factorial(i));
// 	}

// 	int Sm = 0;
// 	int Smp1 = 0;
// 	for (int m=0; m<k-1; m++) {
// 		std::vector<Number_T> crGbuf;
// 		Sm = Sm + s[m];
// 		Smp1 = Sm + s[m+1];
// 		for (int i=0; i<=L2; i++) {
// 			crGbuf.push_back(cln::factorial(i + Sm - m - 2) / cln::factorial(i + Smp1 - m - 2));
// 		}
// 		crG.push_back(crGbuf);
// 	}

// 	crX = crB;

// 	for (int m=0; m<k-1; m++) {
// 		std::vector<Number_T> Xbuf;
// 		for (int i=0; i<=L2; i++) {
// 			Xbuf.push_back(crX[i] * crG[m][i]);
// 		}
// 		halfcyclic_convolute(Xbuf, crB, crX);
// 	}
}

// [Cra] section 4
numeric crandall_Y_loop(const numeric& Sqk)
{
// 	cln::cl_F one = cln::cl_float(1, cln::float_format(Digits));
// 	Number_T factor = cln::expt(lambda, Sqk);
// 	Number_T res = factor / Sqk * crX[0] * one;
// 	Number_T resbuf;
// 	int N = 0;
// 	do {
// 		resbuf = res;
// 		factor = factor * lambda;
// 		N++;
// 		res = res + crX[N] * factor / (N+Sqk);
// 	} while ((res != resbuf) || cln::zerop(crX[N]));
// 	return res;
}


// [Cra] section 4
void calc_f(int maxr)
{
// 	f_kj.clear();
// 	f_kj.resize(L1);
	
// 	Number_T t0, t1, t2, t3, t4;
// 	int i, j, k;
// 	std::vector<std::vector<Number_T> >::iterator it = f_kj.begin();
// 	cln::cl_F one = cln::cl_float(1, cln::float_format(Digits));
	
// 	t0 = cln::exp(-lambda);
// 	t2 = 1;
// 	for (k=1; k<=L1; k++) {
// 		t1 = k * lambda;
// 		t2 = t0 * t2;
// 		for (j=1; j<=maxr; j++) {
// 			t3 = 1;
// 			t4 = 1;
// 			for (i=2; i<=j; i++) {
// 				t4 = t4 * (j-i+1);
// 				t3 = t1 * t3 + t4;
// 			}
// 			(*it).push_back(t2 * t3 * cln::expt(Integer_T(k),-j) * one);
// 		}
// 		it++;
// 	}
}


// [Cra] (3.1)
numeric crandall_Z(const std::vector<int>& s)
{
// 	const int j = s.size();

// 	if (j == 1) {	
// 		Number_T t0;
// 		Number_T t0buf;
// 		int q = 0;
// 		do {
// 			t0buf = t0;
// 			q++;
// 			t0 = t0 + f_kj[q+j-2][s[0]-1];
// 		} while (t0 != t0buf);
		
// 		return t0 / cln::factorial(s[0]-1);
// 	}

// 	std::vector<Number_T> t(j);

// 	Number_T t0buf;
// 	int q = 0;
// 	do {
// 		t0buf = t[0];
// 		q++;
// 		t[j-1] = t[j-1] + 1 / cln::expt(Integer_T(q),s[j-1]);
// 		for (int k=j-2; k>=1; k--) {
// 			t[k] = t[k] + t[k+1] / cln::expt(Integer_T(q+j-1-k), s[k]);
// 		}
// 		t[0] = t[0] + t[1] * f_kj[q+j-2][s[0]-1];
// 	} while (t[0] != t0buf);
	
// 	return t[0] / cln::factorial(s[0]-1);
  return 0;
}


// [Cra] (2.4)
numeric zeta_do_sum_Crandall(const std::vector<int>& s)
{
// 	std::vector<int> r = s;
// 	const int j = r.size();

// 	// decide on maximal size of f_kj for crandall_Z
// 	if (Digits < 50) {
// 		L1 = 150;
// 	} else {
// 		L1 = Digits * 3 + j*2;
// 	}

// 	// decide on maximal size of crX for crandall_Y
// 	if (Digits < 38) {
// 		L2 = 63;
// 	} else if (Digits < 86) {
// 		L2 = 127;
// 	} else if (Digits < 192) {
// 		L2 = 255;
// 	} else if (Digits < 394) {
// 		L2 = 511;
// 	} else if (Digits < 808) {
// 		L2 = 1023;
// 	} else {
// 		L2 = 2047;
// 	}

// 	Number_T res;

// 	int maxr = 0;
// 	int S = 0;
// 	for (int i=0; i<j; i++) {
// 		S += r[i];
// 		if (r[i] > maxr) {
// 			maxr = r[i];
// 		}
// 	}

// 	calc_f(maxr);

// 	const Number_T r0factorial = cln::factorial(r[0]-1);

// 	std::vector<int> rz;
// 	int skp1buf;
// 	int Srun = S;
// 	for (int k=r.size()-1; k>0; k--) {

// 		rz.insert(rz.begin(), r.back());
// 		skp1buf = rz.front();
// 		Srun -= skp1buf;
// 		r.pop_back();

// 		initcX(r);
		
// 		for (int q=0; q<skp1buf; q++) {
			
// 			Number_T pp1 = crandall_Y_loop(Srun+q-k);
// 			Number_T pp2 = crandall_Z(rz);

// 			rz.front()--;
			
// 			if (q & 1) {
// 				res = res - pp1 * pp2 / cln::factorial(q);
// 			} else {
// 				res = res + pp1 * pp2 / cln::factorial(q);
// 			}
// 		}
// 		rz.front() = skp1buf;
// 	}
// 	rz.insert(rz.begin(), r.back());

// 	initcX(rz);

// 	res = (res + crandall_Y_loop(S-j)) / r0factorial + crandall_Z(rz);

// 	return res;
  return 0;
}


numeric zeta_do_sum_simple(const std::vector<int>& r)
{
// 	const int j = r.size();

// 	// buffer for subsums
// 	std::vector<Number_T> t(j);
// 	cln::cl_F one = cln::cl_float(1, cln::float_format(Digits));

// 	Number_T t0buf;
// 	int q = 0;
// 	do {
// 		t0buf = t[0];
// 		q++;
// 		t[j-1] = t[j-1] + one / cln::expt(Integer_T(q),r[j-1]);
// 		for (int k=j-2; k>=0; k--) {
// 			t[k] = t[k] + one * t[k+1] / cln::expt(Integer_T(q+j-1-k), r[k]);
// 		}
// 	} while (t[0] != t0buf);

// 	return t[0];
  return 0;
}


// does Hoelder convolution. see [BBB] (7.0)
numeric zeta_do_Hoelder_convolution(const std::vector<int>& m_, const std::vector<int>& s_)
{
	// prepare parameters
	// holds Li arguments in [BBB] notation
	std::vector<int> s = s_;
	std::vector<int> m_p = m_;
	std::vector<int> m_q;
	// holds Li arguments in nested sums notation
	std::vector<Number_T> s_p(s.size(), Number_T(1));
	s_p[0] = s_p[0] * Number_T("1/2");
	// convert notations
	int sig = 1;
	for (int i=0; i<s_.size(); i++) {
		if (s_[i] < 0) {
			sig = -sig;
			s_p[i] = -s_p[i];
		}
		s[i] = sig * std::abs(s[i]);
	}
	std::vector<Number_T> s_q;
	Number_T signum = 1;

	// first term
	Number_T res = multipleLi_do_sum(m_p, s_p);

	// middle terms
	do {

		// change parameters
		if (s.front() > 0) {
			if (m_p.front() == 1) {
				m_p.erase(m_p.begin());
				s_p.erase(s_p.begin());
				if (s_p.size() > 0) {
					s_p.front() = s_p.front() * Number_T("1/2");
				}
				s.erase(s.begin());
				m_q.front()++;
			} else {
				m_p.front()--;
				m_q.insert(m_q.begin(), 1);
				if (s_q.size() > 0) {
					s_q.front() = s_q.front() * 2;
				}
				s_q.insert(s_q.begin(), Number_T("1/2"));
			}
		} else {
			if (m_p.front() == 1) {
				m_p.erase(m_p.begin());
				Number_T spbuf = s_p.front();
				s_p.erase(s_p.begin());
				if (s_p.size() > 0) {
					s_p.front() = s_p.front() * spbuf;
				}
				s.erase(s.begin());
				m_q.insert(m_q.begin(), 1);
				if (s_q.size() > 0) {
					s_q.front() = s_q.front() * 4;
				}
				s_q.insert(s_q.begin(), Number_T("1/4"));
				signum = -signum;
			} else {
				m_p.front()--;
				m_q.insert(m_q.begin(), 1);
				if (s_q.size() > 0) {
					s_q.front() = s_q.front() * 2;
				}
				s_q.insert(s_q.begin(), Number_T("1/2"));
			}
		}

		// exiting the loop
		if (m_p.size() == 0) break;

		res = res + signum * multipleLi_do_sum(m_p, s_p) * multipleLi_do_sum(m_q, s_q);

	} while (true);

	// last term
	res = res + signum * multipleLi_do_sum(m_q, s_q);

	return res;
}


} // end of anonymous namespace

 */

} // namespace GiNaC

