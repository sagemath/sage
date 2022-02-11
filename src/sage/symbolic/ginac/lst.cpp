/** @file lst.cpp
 *
 *  Implementation of GiNaC's lst. */

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

#include "lst.h"
#include "wildcard.h"
#include "cmatcher.h"

namespace GiNaC {

        // 
template<> const tinfo_static_t lst::tinfo_static = {};
template <> registered_class_info lst::reg_info = registered_class_info(registered_class_options("lst", "basic", &lst::tinfo_static).print_func<print_context>(&lst::do_print).print_func<print_tree>(&lst::do_print_tree));

/** Specialization of container::get_tinfo() for lst. */
template<> tinfo_t lst::get_tinfo() { return &lst::tinfo_static; }

/** Specialization of container::info() for lst. */
template <> bool lst::info(unsigned inf) const
{
	if (inf == info_flags::list)
		return true;

	return inherited::info(inf);
}

template <>
bool lst::match(const ex & pattern, exmap& map) const
{
	if (is_exactly_a<wildcard>(pattern)) {
                const auto& it = map.find(pattern);
                if (it != map.end())
		        return is_equal(ex_to<basic>(it->second));
		map[pattern] = *this;
		return true;
	} 
        if (not is_exactly_a<lst>(pattern))
                return false;
        CMatcher cm(*this, pattern, map);
        const opt_exmap& m = cm.get();
        if (not m)
                return false;
        map = m.value();
        return true;
}

template <> ex lst::unarchive(const archive_node &n, lst &sym_lst)
{
        return (new lst(n, sym_lst))->
                setflag(status_flags::dynallocated);
}

} // namespace GiNaC
