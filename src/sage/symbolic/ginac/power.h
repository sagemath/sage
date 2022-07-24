/** @file power.h
 *
 *  Interface to GiNaC's symbolic exponentiation (basis^exponent). */

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

#ifndef __GINAC_POWER_H__
#define __GINAC_POWER_H__

#include "py_funcs.h"
#include "basic.h"
#include "ex.h"
#include "numeric.h"

namespace GiNaC {

class numeric;
class add;
class mul;

/** This class holds a two-component object, a basis and exponent
 *  representing exponentiation. */
class power : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(power, basic)
	
	friend class print_order;
	friend class print_order_mul;
	friend class ex;
	friend class mul;
	
// member functions
	
	// other constructors
public:
	power(ex  lh, ex  rh) : inherited(&power::tinfo_static), basis(std::move(lh)), exponent(std::move(rh)) {}
	template<typename T> power(ex  lh, const T & rh) : inherited(&power::tinfo_static), basis(std::move(lh)), exponent(rh) {}
        static ex unarchive(const archive_node &n, lst &sym_lst)
        {
                return (new power(n, sym_lst))->
                        setflag(status_flags::dynallocated);
        }
	
	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const override {return 60;}
	bool info(unsigned inf) const override;
	size_t nops() const override;
	const ex op(size_t i) const override;
	ex & let_op(size_t i) override;
	ex map(map_function & f) const override;
	bool is_polynomial(const ex & var) const override;
	numeric degree(const ex & s) const override;
	numeric ldegree(const ex & s) const override;
	ex coeff(const ex & s, const ex & n) const override;
	ex eval(int level=0) const override;
	ex evalf(int level=0, PyObject* parent=nullptr) const override;
	ex series(const relational & s, int order, unsigned options = 0) const override;
        void useries(flint_series_t& fp, int order) const override;
        bool match(const ex& pattern, exmap& map) const override;
	ex subs(const exmap & m, unsigned options = 0) const override;
	bool has(const ex & other, unsigned options = 0) const override;
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0, unsigned options = 0) const override;
	ex to_rational(exmap & repl) const override;
	ex to_polynomial(exmap & repl) const override;
	ex conjugate() const override;
	ex real_part() const override;
	ex imag_part() const override;
	//int compare(const basic& other) const;
	//int compare_symbol(const symbol& other) const;
protected:
	ex derivative(const symbol & s) const override;
	unsigned return_type() const override;
	tinfo_t return_type_tinfo() const override;
	ex expand(unsigned options = 0) const override;
	long calchash() const override;

	// new virtual functions which can be overridden by derived classes
	// none
	
	// non-virtual functions in this class
protected:
	void print_power(const print_context & c, const char *powersymbol, const char *openbrace, const char *closebrace, unsigned level) const;
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_python(const print_python & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;

	ex expand_add(const add & a, long int n, unsigned options) const;
	ex expand_add_2(const add & a, unsigned options) const;
	ex expand_mul(const mul & m, const numeric & n, unsigned options, bool from_expand = false) const;
	
// member variables
	
protected:
	ex basis;
	ex exponent;
};

// wrapper functions

/** Symbolic exponentiation.  Returns a power-object as a new expression.
 *
 *  @param b the basis expression
 *  @param e the exponent expression */
inline ex pow(const ex & b, const ex & e)
{
	return power(b, e);
}
template<typename T1, typename T2>
inline ex pow(const T1 & b, const T2 & e)
{
	return power(ex(b), ex(e));
}

/** Square root expression.  Returns a power-object with exponent 1/2. */
inline ex sqrt(const ex & a)
{
	extern const ex _ex1_2;
	return power(a,_ex1_2);
}

} // namespace GiNaC

#endif // ndef __GINAC_POWER_H__
