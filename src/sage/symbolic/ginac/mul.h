/** @file mul.h
 *
 *  Interface to GiNaC's products of expressions. */

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

#ifndef __GINAC_MUL_H__
#define __GINAC_MUL_H__

#include "expairseq.h"
#include "power.h"

namespace GiNaC {

/** Product of expressions. */
class mul : public expairseq
{
	GINAC_DECLARE_REGISTERED_CLASS(mul, expairseq)
	
	friend class print_order;
	friend class ex;
	friend class add;
	friend class power;
	
	// other constructors
public:
	mul(const ex & lh, const ex & rh);
	mul(const exvector & v, bool hold=false);
	mul(const epvector & v);
	mul(const epvector & v, const numeric & oc,
                        bool do_index_renaming = false);
//	mul(std::unique_ptr<epvector> vp, const ex & oc, bool do_index_renaming = false);
	mul(const ex & lh, const ex & mh, const ex & rh);
        static ex unarchive(const archive_node &n, lst &sym_lst)
        {
                return (new mul(n, sym_lst))->
                        setflag(status_flags::dynallocated);
        }
	
	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const override {return 50;}
	bool info(unsigned inf) const override;
	bool is_polynomial(const ex & var) const override;
	numeric degree(const ex & s) const override;
	numeric ldegree(const ex & s) const override;
	ex coeff(const ex & s, const ex & n) const override;
	bool has(const ex & other, unsigned options = 0) const override;
	ex eval(int level=0) const override;
	ex evalf(int level=0, PyObject* parent=nullptr) const override;
	ex real_part() const override;
	ex imag_part() const override;
	ex series(const relational & s, int order, unsigned options = 0) const override;
        void useries(flint_series_t& fp, int order) const override;
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0, unsigned options = 0) const override;
	numeric integer_content() const override;
	ex smod(const numeric &xi) const override;
	numeric max_coefficient() const override;
	ex conjugate() const override;
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
	bool expair_needs_further_processing(epp it) override;
	const numeric& default_overall_coeff() const override;
	void combine_overall_coeff(const numeric & c) override;
	void combine_overall_coeff(const numeric & c1, const numeric & c2) override;
	bool can_make_flat(const expair & p) const override;
	ex expand(unsigned options=0) const override;
	void find_real_imag(ex&, ex&) const;
	//int compare(const basic& other) const;
	
	ex eval_infinity(epvector::const_iterator infinity_iter) const;
	ex eval_exponentials() const;

	// new virtual functions which can be overridden by derived classes
	// none
	
	// non-virtual functions in this class
public:
	ex algebraic_subs_mul(const exmap & m, unsigned options) const;
	double total_degree() const;
	const epvector & get_sorted_seq() const override;
	ex recombine_pair_to_ex(const expair & p) const override;
	//int compare_symbol(const symbol &other) const;
	//int compare_pow(const power &other) const;
        ex without_known_factor(const ex& f) const;
protected:
	void print_overall_coeff(const ex& coeff_ex, const print_context & c,
			const char *mul_sym, bool latex=false) const;
	void print_exvector(const exvector & v, const print_context & c,
		const char* sep) const;
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_rat_func(const print_context & c, unsigned level, 
			bool latex_tags) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;
	static bool can_be_further_expanded(const ex & e);
	std::unique_ptr<epvector> expandchildren(unsigned options) const;

	mutable double tdegree;
};

} // namespace GiNaC

#endif // ndef __GINAC_MUL_H__
