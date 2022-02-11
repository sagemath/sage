/** @file tostring.h
 *
 *  Convert object to its string representation (output form). This is an
 *  internal header file. */

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

#ifndef __GINAC_TOSTRING_H__
#define __GINAC_TOSTRING_H__

#include <sstream>

namespace GiNaC {

template<class T>
std::string ToString(const T & t)
{
	std::ostringstream buf;
	buf << t;
	return buf.str();
}

} // namespace GiNaC

#endif // ndef __GINAC_TOSTRING_H__
