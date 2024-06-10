/** @file symbol.h
 *
 *  Interface to GiNaC's symbolic objects. */

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

#ifndef __GINAC_SYMBOL_H__
#define __GINAC_SYMBOL_H__

#include "basic.h"
#include "ex.h"
#include "infoflagbase.h"

#include <string>

namespace GiNaC {

/** Basic CAS symbol.  It has a name because it must know how to output itself.
 *  It may be assigned an expression, but this feature is only intended for
 *  programs like 'ginsh' that want to associate symbols with expressions.
 *  If you want to replace symbols by expressions in your code, you should
 *  use ex::subs() or use objects of class ex instead of class symbol in the
 *  first place. */
class symbol : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(symbol, basic)

	friend class print_order;
        friend struct symbolhasher;

// types
	
// member functions
	
	// other constructors
public:
	explicit symbol(std::string  initname, unsigned domain = domain::complex);
	symbol(std::string  initname, std::string  texname, unsigned domain = domain::complex);
	symbol(const std::string & initname, const std::string & texname, unsigned rt, unsigned domain);
	symbol(std::string  initname, unsigned rt, tinfo_t rtt, unsigned domain = domain::complex);
	symbol(std::string  initname, std::string  texname, unsigned rt, tinfo_t rtt, unsigned domain = domain::complex);
	
	// functions overriding virtual functions from base classes
public:
	//int compare_same_type(const basic & other) const;
	//int compare(const basic &other) const;
	bool info(unsigned inf) const override;
	ex eval(int level = 0) const override;
	ex evalf(int level = 0, PyObject* parent = nullptr) const override { return *this; } // overwrites basic::evalf() for performance reasons
	ex series(const relational & s, int order, unsigned options = 0) const override;
        void useries(flint_series_t& fp, int order) const override;
	ex subs(const exmap & m, unsigned options = 0) const override { return subs_one_level(m, options); } // overwrites basic::subs() for performance reasons
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0, unsigned options = 0) const override;
	ex to_rational(exmap & repl) const override;
	ex to_polynomial(exmap & repl) const override;
	unsigned return_type() const override { return ret_type; }
	tinfo_t return_type_tinfo() const override { return ret_type_tinfo; }
	ex conjugate() const override;
	ex real_part() const override;
	ex imag_part() const override;
	bool is_polynomial(const ex & var) const override;
	ex derivative(const symbol & s) const override;
	bool is_equal_same_type(const basic & other) const override;
	long calchash() const override;
	static ex unarchive(const archive_node &n, lst &sym_lst);
        bool operator==(const symbol& other) const { return serial == other.serial; }
	void set_name(const std::string & n) { name = n; }
	std::string get_name() const { return name; }
	unsigned get_domain() const { return domain; }
	void set_domain(unsigned d);
        void set_domain_from_ex(const ex& expr);
	void set_info(unsigned flag, bool value=true) { iflags.set(flag, value); }
	std::string get_texname() const { return TeX_name; }
	void set_texname(const std::string & t) { TeX_name = t; }
protected:
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_tree(const print_tree & c, unsigned level) const override;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;
private:
	std::string & autoname_prefix();
	std::string default_TeX_name() const;

// member variables

protected:
	unsigned serial;                 ///< unique serial number for comparison
	std::string name;                ///< printname of this symbol
	std::string TeX_name;            ///< LaTeX name of this symbol
	unsigned domain;                 ///< domain of symbol
	infoflagbase iflags;              ///< bitset with info_flags
	unsigned ret_type;               ///< value returned by return_type()
	tinfo_t ret_type_tinfo;         ///< value returned by return_type_tinfo()
private:
	static unsigned next_serial;
};

struct symbolhasher {
        std::size_t operator()(const symbol& sym) const { return sym.serial; }
};
using symbolset = std::unordered_set<symbol,symbolhasher>;

// utility functions

// keep symbols unique
const symbol & get_symbol(const std::string & s);

bool has_symbol(const ex & x);
bool has_symbol(const ex & x, const symbol& s);
bool has_free_symbol(const ex & x, const symbol& s);

} // namespace GiNaC

#endif // ndef __GINAC_SYMBOL_H__
