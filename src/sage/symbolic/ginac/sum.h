/** @file sum.h
 *
 *  Special declarations for summation functions. */

/*
 *  Copyright (C) 2016  Ralf Stephan <ralf@ark.in-berlin.de>
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


namespace GiNaC {

extern ex hypersimp(ex e, ex symbol);
extern ex gosper_term(ex the_ex, ex n);
extern ex gosper_sum_definite(ex the_ex, ex n, ex a, ex b, int* p);
extern ex gosper_sum_indefinite(ex the_ex, ex n, int* p);
extern ex gamma_normalize(ex the_ex);
extern ex to_gamma(const ex& the_ex);

} // namespace GiNaC

