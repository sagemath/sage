/** @file expair.cpp
 *
 *  Implementation of expression pairs (building blocks of expairseq). */

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

#include "expair.h"
#include "operators.h"

#include <iostream>

namespace GiNaC {

void expair::print(std::ostream & os) const
{
	os << "expair:";
	print_tree c(os);
	rest.print(c, c.delta_indent);
	coeff.print(c, c.delta_indent);
}

const expair expair::conjugate() const
{
	ex newrest = rest.conjugate();
	ex newcoeff = coeff.conjugate();
	if (are_ex_trivially_equal(newrest,rest) && are_ex_trivially_equal(newcoeff,coeff)) {
		return *this;
	}
	return expair(newrest, newcoeff);
}

} // namespace GiNaC
