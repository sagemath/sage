/** @file useries.h
 *
 *  Interface to class for extended truncated univariate power series. */

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

#ifndef __PYNAC_USERIES_H__
#define __PYNAC_USERIES_H__

#include "pseries.h"
#include "expairseq.h"


namespace GiNaC {

bool useries_can_handle(const ex& the_ex, const symbol& s);
ex useries(const ex& the_ex, const symbol& s, int order, unsigned options = 0);

} // namespace GiNaC

#endif // ndef __PYNAC_USERIES_H__
