/** @file symbol.cpp
 *
 *  Implementation of GiNaC's symbolic objects. */

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

#include "symbol.h"
#include "power.h"
#include "mul.h"
#include "lst.h"
#include "archive.h"
#include "tostring.h"
#include "utils.h"
#include "inifcns.h"

#include <string>
#include <stdexcept>
#ifdef DO_GINAC_ASSERT
#  include <typeinfo>
#endif

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(symbol, basic,
  print_func<print_context>(&symbol::do_print).
  print_func<print_latex>(&symbol::do_print_latex).
  print_func<print_tree>(&symbol::do_print_tree).
  print_func<print_python_repr>(&symbol::do_print_python_repr))

//////////
// default constructor
//////////

// symbol

symbol::symbol()
 : inherited(&symbol::tinfo_static), serial(next_serial++), name(autoname_prefix() + ToString(serial)), TeX_name(name), domain(domain::complex), ret_type(return_types::commutative), ret_type_tinfo(&symbol::tinfo_static)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

// public

// symbol

symbol::symbol(std::string  initname, unsigned a_domain)
 : inherited(&symbol::tinfo_static), serial(next_serial++), name(std::move(initname)), TeX_name(default_TeX_name()), domain(a_domain), ret_type(return_types::commutative), ret_type_tinfo(&symbol::tinfo_static)
{
        set_domain(a_domain);
	setflag(status_flags::evaluated | status_flags::expanded);
}

symbol::symbol(std::string  initname, unsigned rt, tinfo_t rtt, unsigned a_domain)
 : inherited(&symbol::tinfo_static), serial(next_serial++), name(std::move(initname)), TeX_name(default_TeX_name()), domain(a_domain), ret_type(rt), ret_type_tinfo(rtt)
{
        set_domain(a_domain);
	setflag(status_flags::evaluated | status_flags::expanded);
}

symbol::symbol(std::string  initname, std::string  texname, unsigned a_domain)
 : inherited(&symbol::tinfo_static), serial(next_serial++), name(std::move(initname)), TeX_name(std::move(texname)), domain(a_domain), ret_type(return_types::commutative), ret_type_tinfo(&symbol::tinfo_static)
{
        set_domain(a_domain);
	setflag(status_flags::evaluated | status_flags::expanded);
}

symbol::symbol(std::string  initname, std::string  texname, unsigned rt, tinfo_t rtt, unsigned a_domain)
 : inherited(&symbol::tinfo_static), serial(next_serial++), name(std::move(initname)), TeX_name(std::move(texname)), domain(a_domain), ret_type(rt), ret_type_tinfo(rtt)
{
        set_domain(a_domain);
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// archiving
//////////

/** Construct object from archive_node. */
symbol::symbol(const archive_node &n, lst &sym_lst)
 : inherited(n, sym_lst), serial(next_serial++)
{
	if (!n.find_string("name", name))
		name = autoname_prefix() + ToString(serial);
	if (!n.find_string("TeX_name", TeX_name))
		TeX_name = default_TeX_name();
	if (!n.find_unsigned("domain", domain))
		domain = domain::complex;
	if (!n.find_unsigned("return_type", ret_type))
		ret_type = return_types::commutative;
        set_domain(domain);
	setflag(status_flags::evaluated | status_flags::expanded);
}

/** Unarchive the object. */
ex symbol::unarchive(const archive_node &n, lst &sym_lst)
{
	ex s = (new symbol(n, sym_lst))->setflag(status_flags::dynallocated);

	// If symbol is in sym_lst, return the existing symbol
	for (const auto & elem : sym_lst) {
		if (is_exactly_a<symbol>(elem) && (ex_to<symbol>(elem).name == ex_to<symbol>(s).name))
			return elem;
	}

	// Otherwise add new symbol to list and return it
	sym_lst.append(s);
	return s;
}

/** Archive the object. */
void symbol::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_string("name", name);
	if (TeX_name != default_TeX_name())
		n.add_string("TeX_name", TeX_name);
	if (domain != domain::complex)
		n.add_unsigned("domain", domain);
	if (ret_type != return_types::commutative)
		n.add_unsigned("return_type", ret_type);
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void symbol::do_print(const print_context & c, unsigned level) const
{
	c.s << name;
}

void symbol::do_print_latex(const print_latex & c, unsigned level) const
{
	c.s << TeX_name;
}

void symbol::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << name << " (" << class_name() << ")" << " @" << this
	    << ", serial=" << serial
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", domain=" << domain << ", iflags=" << iflags.to_string()
	    << std::endl;
}

void symbol::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << "('" << name;
	if (TeX_name != default_TeX_name())
		c.s << "','" << TeX_name;
	c.s << "')";
}

void symbol::set_domain(unsigned d)
{
        infoflagbase nb;
        if (not (d == domain::complex))
                nb = iflags;
        switch (d) {
        case domain::complex:
                break;
        case domain::real:
                nb.set(info_flags::real, true);
                break;
        case domain::positive:
                nb.set(info_flags::real, true);
                nb.set(info_flags::positive, true);
                nb.set(info_flags::nonzero, true);
                break;
        case domain::negative:
                nb.set(info_flags::real, true);
                nb.set(info_flags::negative, true);
                nb.set(info_flags::nonzero, true);
                break;
        case domain::rational:
                nb.set(info_flags::real, true);
                nb.set(info_flags::rational, true);
                break;
        case domain::integer:
                nb.set(info_flags::real, true);
                nb.set(info_flags::rational, true);
                nb.set(info_flags::integer, true);
                break;
        case domain::even:
                nb.set(info_flags::real, true);
                nb.set(info_flags::rational, true);
                nb.set(info_flags::integer, true);
                nb.set(info_flags::even, true);
                break;
        }
        iflags = nb;
}

void symbol::set_domain_from_ex(const ex& expr)
{
        iflags.clear();
        if (expr.info(info_flags::even))
                set_domain(domain::even);
        else if (expr.is_integer())
                set_domain(domain::integer);
        else if (expr.info(info_flags::rational))
                set_domain(domain::rational);
        else if (expr.is_real())
                set_domain(domain::real);

        if (expr.is_positive())
                set_domain(domain::positive);
        else if (expr.info(info_flags::negative))
                set_domain(domain::negative);
}

bool symbol::info(unsigned inf) const
{
        switch (inf) {
        case info_flags::symbol:
        case info_flags::polynomial:
        case info_flags::integer_polynomial:
        case info_flags::cinteger_polynomial:
        case info_flags::rational_polynomial:
        case info_flags::crational_polynomial:
        case info_flags::rational_function:
        case info_flags::expanded:
                return true;
        case info_flags::nonzero:
                return iflags.get(info_flags::nonzero)
                or iflags.get(info_flags::positive)
                or iflags.get(info_flags::negative);
        case info_flags::infinity:
                return false;
        default:
                return iflags.get(inf);
        }
        return false;
}

ex symbol::eval(int level) const
{
	if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	return this->hold();
}

ex symbol::conjugate() const
{
	if (iflags.get(info_flags::real))
		return *this;
        return conjugate_function(*this).hold();
}

ex symbol::real_part() const
{
	if (iflags.get(info_flags::real))
		return *this;
	return real_part_function(*this).hold();
}

ex symbol::imag_part() const
{
	if (iflags.get(info_flags::real))
		return _ex0;
	return imag_part_function(*this).hold();
}

bool symbol::is_polynomial(const ex & var) const
{
	return true;
}

// protected

/** Implementation of ex::diff() for single differentiation of a symbol.
 *  It returns 1 or 0.
 *
 *  @see ex::diff */
ex symbol::derivative(const symbol & s) const
{
	if (compare_same_type(s) != 0)
		return _ex0;
	return _ex1;
}

int symbol::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<symbol>(other));
	const symbol *o = static_cast<const symbol *>(&other);
	if (serial==o->serial) return 0;
        return serial < o->serial ? -1 : 1;
}

bool symbol::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<symbol>(other));
	const symbol *o = static_cast<const symbol *>(&other);
	return serial==o->serial;
}

long symbol::calchash() const
{
  // Original code
    hashvalue = golden_ratio_hash((intptr_t)tinfo() ^ serial);
    setflag(status_flags::hash_calculated);
    return hashvalue;
}

//////////
// virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// public

// private

/** Symbols not constructed with a string get one assigned using this
 *  prefix and a number. */
std::string & symbol::autoname_prefix()
{
	static std::string s("symbol");
	return s;
}

/** Return default TeX name for symbol. This recognizes some greek letters. */
std::string symbol::default_TeX_name() const
{
	if (name=="alpha"        || name=="beta"         || name=="gamma"
	 || name=="delta"        || name=="epsilon"      || name=="varepsilon"
	 || name=="zeta"         || name=="eta"          || name=="theta"
	 || name=="vartheta"     || name=="iota"         || name=="kappa"
	 || name=="lambda"       || name=="mu"           || name=="nu"
	 || name=="xi"           || name=="omicron"      || name=="pi"
	 || name=="varpi"        || name=="rho"          || name=="varrho"
	 || name=="sigma"        || name=="varsigma"     || name=="tau"
	 || name=="upsilon"      || name=="phi"          || name=="varphi"
	 || name=="chi"          || name=="psi"          || name=="omega"
	 || name=="Gamma"        || name=="Delta"        || name=="Theta"
	 || name=="Lambda"       || name=="Xi"           || name=="Pi"
	 || name=="Sigma"        || name=="Upsilon"      || name=="Phi"
	 || name=="Psi"          || name=="Omega")
		return "\\" + name;
	
	return name;
}

//////////
// static member variables
//////////

// private

unsigned symbol::next_serial = 0;

// utility function to keep only one instance of a symbol with a given name
const symbol & get_symbol(const std::string & s)
{
	static std::map<std::string, symbol> directory;
	auto i = directory.find(s);
	if (i != directory.end()) {
		return i->second;
	}
	return directory.insert(make_pair(s, symbol(s))).first->second;
}

bool has_symbol(const ex & x)
{
	if (is_exactly_a<symbol>(x))
		return true;
	for (size_t i=0; i<x.nops(); ++i)
		if (has_symbol(x.op(i)))
			return true;
                
	return false;
}

bool has_symbol(const ex & x, const symbol& s)
{
	if (is_exactly_a<symbol>(x) and ex_to<symbol>(x) == s)
		return true;
	for (size_t i=0; i<x.nops(); ++i)
		if (has_symbol(x.op(i), s))
			return true;

	return false;
}

bool has_free_symbol(const ex & x, const symbol& s)
{
        if (is_exactly_a<function>(x))
                return false;
	if (is_exactly_a<symbol>(x) and ex_to<symbol>(x) == s)
		return true;
	for (size_t i=0; i<x.nops(); ++i)
		if (has_free_symbol(x.op(i), s))
			return true;

	return false;
}

} // namespace GiNaC
