/** @file infinity.h
 *
 *  The value "Infinity".  */

/*
 *  Copyright (C) 2011 Volker Braun <vbraun@stp.dias.ie>
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

#ifndef __GINAC_INFINITY_H__
#define __GINAC_INFINITY_H__

#include "basic.h"
#include "numeric.h"
#include "ex.h"
#include "relational.h"

#include <string>

namespace GiNaC {

	
/** This class holds "infinity"
 *  It includes a direction (like -infinity).
 **/
class infinity : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(infinity, basic)
public:
       	infinity(const numeric & _direction);
	static infinity from_direction(const ex & _direction);
	static infinity from_sign(int sgn);
        static ex unarchive(const archive_node &n, lst &sym_lst);
	
	// functions overriding virtual functions from base classes
public:
	bool info(unsigned inf) const override;
	ex evalf(int level = 0, PyObject* parent=nullptr) const override;
	ex conjugate() const override;
	ex real_part() const override;
	ex imag_part() const override;

	bool is_unsigned_infinity() const;
	bool is_plus_infinity() const;
	bool is_minus_infinity() const;
	const ex & get_direction() const { return direction; };
	const infinity & operator *= (const ex & rhs);
	const infinity & operator += (const ex & rhs);

	bool compare_other_type(const ex & other,
		relational::operators op) const;

protected:
	ex derivative(const symbol & s) const override;
	bool is_equal_same_type(const basic & other) const override;
	long calchash() const override { return hashvalue; }
	
	// non-virtual functions in this class
protected:
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_tree(const print_tree & c, unsigned level) const override;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;

	void set_direction(const ex & new_direction);
	
// member variables
private:
	ex direction;
};

extern const infinity Infinity;
extern const infinity NegInfinity;
extern const infinity UnsignedInfinity;

} // namespace GiNaC

#endif // ndef __GINAC_INFINITY_H__
