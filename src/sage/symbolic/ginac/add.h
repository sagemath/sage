/** @file add.h
 *
 *  Interface to GiNaC's sums of expressions. */

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

#ifndef __GINAC_ADD_H__
#define __GINAC_ADD_H__

#include "expairseq.h"

namespace GiNaC {

/** Sum of expressions. */
class add : public expairseq
{
	GINAC_DECLARE_REGISTERED_CLASS(add, expairseq)
	
        friend class ex;
	friend class mul;
	friend class power;
        friend ex Power(const ex& base_, const ex& expo_);
	
	// other constructors
public:
	add(const ex & lh, const ex & rh);
	add(const exvector & v, bool hold=false);
	add(const epvector & v);
	add(const epvector & v, const numeric & oc);
//	add(std::unique_ptr<epvector> vp, const ex & oc);
        static ex unarchive(const archive_node &n, lst &sym_lst)
        {
                return (new add(n, sym_lst))->
                        setflag(status_flags::dynallocated);
        }
	
	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const override {return 40;}
	bool info(unsigned inf) const override;
	bool is_polynomial(const ex & var) const override;
	numeric degree(const ex & s) const override;
	numeric ldegree(const ex & s) const override;
	ex coeff(const ex & s, const ex & n) const override;
	ex eval(int level=0) const override;
	ex series(const relational & r, int order, unsigned options = 0) const override;
        void useries(flint_series_t& fp, int order) const override;
	ex normal(exmap & repl, exmap & rev_lookup, int level=0, unsigned options = 0) const override;
	numeric integer_content() const override;
	ex smod(const numeric &xi) const override;
	numeric max_coefficient() const override;
	ex conjugate() const override;
	ex real_part() const override;
	ex imag_part() const override;
	const epvector & get_sorted_seq() const override;
	ex lead_coeff() const;
protected:
	ex derivative(const symbol & s) const override;
	unsigned return_type() const override;
	tinfo_t return_type_tinfo() const override;
	ex thisexpairseq(const epvector & v, const numeric & oc, bool do_index_renaming = false) const override;
	ex thisexpairseq(std::unique_ptr<epvector> vp, const numeric & oc, bool do_index_renaming = false) const override;
	expair split_ex_to_pair(const ex & e) const override;
	expair combine_ex_with_coeff_to_pair(const ex & e,
	                                     const numeric & c) const override;
	expair combine_pair_with_coeff_to_pair(const expair & p,
	                                     const numeric & c) const override;
	ex recombine_pair_to_ex(const expair & p) const override;
	ex power(const numeric& expo) const;
	ex expand(unsigned options=0) const override;
	ex eval_infinity(epvector::const_iterator infinity_iter) const;

	// non-virtual functions in this class
protected:
	void print_add(const print_context & c, unsigned level, bool latex) const;
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;

public:
        ex combine_fractions() const;
};

ex Power(const ex& base_, const ex& expo_);

} // namespace GiNaC

#endif // ndef __GINAC_ADD_H__
