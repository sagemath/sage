/** @file mpoly.h
 *
 *  This file defines several functions that work on multivariate polynomials.

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

#ifndef __GINAC_MPOLY_H__
#define __GINAC_MPOLY_H__

#include "lst.h"

namespace GiNaC {

class ex;
class symbol;

// Polynomial GCD in Z[X], cofactors are returned in ca and cb, if desired
extern ex gcdpoly(const ex &a, const ex &b, ex *ca = nullptr, ex *cb = nullptr, bool check_args = true);
extern bool factorpoly(const ex& p, ex& res);
extern ex poly_mul_expand(const ex &a, const ex &b); 

// Polynomial LCM in Z[X]
extern ex lcm(const ex &a, const ex &b, bool check_args = true);

extern numeric lcm_of_coefficients_denominators(const ex &e);
extern ex multiply_lcm(const ex &e, const numeric &lcm);

// Square-free factorization of a polynomial a(x)
extern ex sqrfree(const ex &a, const lst &l = lst());

// Square-free partial fraction decomposition of a rational function a(x)
extern ex sqrfree_parfrac(const ex & a, const symbol & x);

// Collect common factors in sums.
extern ex collect_common_factors(const ex & e);

// Resultant of two polynomials e1,e2 with respect to symbol s.
extern ex resultant(const ex & e1, const ex & e2, const ex & s);
extern ex resultantpoly(const ex & ee1, const ex & ee2, const ex & s);

} // namespace GiNaC

#endif // ndef __GINAC_MPOLY_H__
