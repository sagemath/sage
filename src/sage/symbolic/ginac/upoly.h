/** @file upoly.h
 *
 *  This file declares several functions that work on univariate and
 *  rational functions.
 *  These functions include polynomial quotient and remainder, GCD and LCM
 *  computation, square-free factorization and rational function normalization. */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *        Copyright (C) 2016 Ralf Stephan
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

#ifndef __PYNAC_UPOLY_H__
#define __PYNAC_UPOLY_H__

#include <utility>

namespace GiNaC {

class ex;
class symbol;

// Quotient q(x) of polynomials a(x) and b(x) in Q[x], so that a(x)=b(x)*q(x)+r(x)
extern ex quo(const ex &a, const ex &b, const ex &x, bool check_args = true);

// Remainder r(x) of polynomials a(x) and b(x) in Q[x], so that a(x)=b(x)*q(x)+r(x)
extern ex rem(const ex &a, const ex &b, const ex &x, bool check_args = true);
extern std::pair<ex,ex> quo_rem(const ex &a, const ex &b, const ex &x,
                bool check_args);

// Decompose rational function a(x)=N(x)/D(x) into Q(x)+R(x)/D(x) with degree(R, x) < degree(D, x)
extern ex decomp_rational(const ex &a, const ex &x);

// Pseudo-remainder of polynomials a(x) and b(x) in Q[x]
extern ex prem(const ex &a, const ex &b, const ex &x, bool check_args = true);

// Pseudo-remainder of polynomials a(x) and b(x) in Q[x]
extern ex sprem(const ex &a, const ex &b, const ex &x, bool check_args = true);

// Exact polynomial division of a(X) by b(X) in Q[X] (quotient returned in q), returns false when exact division fails
extern bool divide(const ex &a, const ex &b, ex &q, bool check_args = true);

extern ex parfrac(const ex & a, const ex & x);

} // namespace GiNaC

#endif // __PYNAC_UPOLY_H__
