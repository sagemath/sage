/** @file order.h
 *
 *  Definitions of order used for printing. */

/*
 *   Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
 *   Copyright (C) 2011 Jean-Pierre Flori <flori@enst.fr>
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

#ifndef __GINAC_ORDER_H__
#define __GINAC_ORDER_H__

#include "ex.h"
#include "basic.h"
#include "mul.h"
#include "power.h"
#include "add.h"
#include "symbol.h"
#include "function.h"
#include "fderivative.h"

namespace GiNaC {

class print_order {
private:
	const tinfo_t& function_id() const;
	const tinfo_t& fderivative_id() const;
	const tinfo_t& power_id() const;
	const tinfo_t& symbol_id() const;
	const tinfo_t& mul_id() const;
	const tinfo_t& add_id() const;
	const tinfo_t& numeric_id() const;
	const tinfo_t& constant_id() const;
	const tinfo_t& wildcard_id() const;
	const tinfo_t& pseries_id() const;

public:
	virtual ~print_order() {}
	bool operator() (const ex &lh, const ex &rh) const;
	int compare(const ex &lh, const ex &rh) const;

protected:
	friend class print_order_pair;

	//	int compare(const ptr<basic> &lh, const ptr<basic> &rh) const;

	// all other compare() methods dispatch to this one
	int compare(const basic &lh, const basic &rh) const;

	// mul objects
	int compare_mul_symbol(const mul &lh, const symbol &rh) const;
	int compare_mul_power(const mul &lh, const power &rh) const;
	int compare_same_type_mul(const mul &lh, const mul &rh) const;
	// add objects
	int compare_add_symbol(const add &lh, const symbol &rh) const;
	int compare_add_power(const add &lh, const power &rh) const;
	int compare_add_mul(const add &lh, const mul &rh) const;
	int compare_same_type_add(const add &lh, const add &rh) const;
	// power objects
	virtual int compare_power_symbol(const power &lh, const symbol &rh) const;
	virtual int compare_same_type_power(const power &lh, const power &rh) const;
	// symbol objects
	int compare_same_type_symbol(const symbol &lh, const symbol &rh) const;
	// container objects
	template <template <class T, class = std::allocator<T> > class C>
	int compare_same_type_container(const container<C> &lh, const container<C> &rh) const;
	// function objects
	int compare_same_type_function(const function &lh, const function &rh) const;
	// fderivative objects
	int compare_same_type_fderivative(const fderivative &lh, const fderivative &rh) const;

	int generic_compare(const tinfo_t, const tinfo_t) const
		{ return 1; }
};


class print_order_mul : public print_order {
	int compare_power_symbol(const power &lh, const symbol &rh) const override;
	int compare_same_type_power(const power &lh, const power &rh) const override;
};


// We have to define the following class to sort held expressions
// E.g. 3*x+2*x which does not get simplified to 5*x.
class print_order_pair {
public:
	bool operator() (const expair &lh, const expair &rh) const;
	bool compare_degrees(const expair &lhex, const expair &rhex) const;
};


class print_order_pair_mul : public print_order_pair
{
public:
	bool operator() (const expair &lh, const expair &rh) const;
};

} // namespace GiNaC
#endif // ndef __GINAC_ORDER_H__
