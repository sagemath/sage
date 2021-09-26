/** @file normal.h
 *
 *  This file defines several functions that work on univariate and
 *  multivariate polynomials and rational functions.
 *  These functions include polynomial quotient and remainder, GCD and LCM
 *  computation, square-free factorization and rational function normalization. */

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

#ifndef __GINAC_NORMAL_H__
#define __GINAC_NORMAL_H__

namespace GiNaC {
class ex;

// Collect common factors in sums.
extern ex collect_common_factors(const ex & e);
extern ex gcd(const ex & e1, const ex & e2);
extern bool factor(const ex & the_ex, ex & res_ex);

} // namespace GiNaC

#endif // ndef __GINAC_NORMAL_H__
