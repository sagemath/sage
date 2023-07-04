/** @file exprseq.h
 *
 *  Definition of GiNaC's exprseq. */

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

#ifndef __GINAC_EXPRSEQ_H__
#define __GINAC_EXPRSEQ_H__

#include "container.h"

#include <vector>

namespace GiNaC {

typedef container<std::vector> exprseq;

/** Specialization of container::get_tinfo() for exprseq. */
template<> tinfo_t exprseq::get_tinfo();

// defined in exprseq.cpp
template<> bool exprseq::info(unsigned inf) const;
template<> ex exprseq::unarchive(const archive_node &n, lst &sym_lst);

} // namespace GiNaC

#endif // ndef __GINAC_EXPRSEQ_H__
