/** @file lst.h
 *
 *  Definition of GiNaC's lst. */

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

#ifndef __GINAC_LST_H__
#define __GINAC_LST_H__

#include "container.h"

#include <list>

namespace GiNaC {

typedef container<std::list> lst;

/** Specialization of container::get_tinfo() for lst. */
template<> tinfo_t lst::get_tinfo();

/** Specialization of container::get_default_flags() for lst. */
template<> inline unsigned lst::get_default_flags() { return status_flags::not_shareable; }

/** Specialization of container::get_open_delim() for lst. */
template<> inline const char* lst::get_open_delim() { return "{"; }

/** Specialization of container::get_close_delim() for lst. */
template<> inline const char* lst::get_close_delim() { return "}"; }

// defined in lst.cpp
template<> bool lst::info(unsigned inf) const;
template<> ex lst::unarchive(const archive_node &n, lst &sym_lst);

} // namespace GiNaC

#endif // ndef __GINAC_LST_H__
