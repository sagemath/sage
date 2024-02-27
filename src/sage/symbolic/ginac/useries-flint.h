/** @file useries.h
 *
 *  Special declarations for flint series. */

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

#ifndef __PYNAC_USERIES_FLINT_H__
#define __PYNAC_USERIES_FLINT_H__

#include "gmp.h"
#include "flint/fmpq_poly.h"
#include "flint/fmpq.h"

#include <stdexcept>


namespace GiNaC {

class flint_error : public std::runtime_error {
    public:
        flint_error() : std::runtime_error("") {}
};

class flint_series_t {
    public:
        flint_series_t() : offset(0) { fmpq_poly_init(ft); }
        ~flint_series_t() { fmpq_poly_clear(ft); }
        int offset;
        fmpq_poly_t ft;
};

} // namespace GiNaC

#endif // ndef __PYNAC_USERIES_FLINT_H__
