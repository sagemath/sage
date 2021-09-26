
/** @file expairseq.h
 *
 *  Interface to sequences of expression pairs. */

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

#ifndef __GINAC_EXPAIRSEQ_H__
#define __GINAC_EXPAIRSEQ_H__

#include "expair.h"

#include <vector>
#include <list>
#include <memory>
// CINT needs <algorithm> to work properly with <vector> and <list>
#include <algorithm>

namespace GiNaC {

/** Using hash tables can potentially enhance the asymptotic behaviour of
 *  combining n terms into one large sum (or n terms into one large product)
 *  from O(n*log(n)) to about O(n).  There are, however, several drawbacks.
 *  The constant in front of O(n) is quite large, when copying such an object
 *  one also has to copy the hash table, comparison is quite expensive because
 *  there is no ordering any more, it doesn't help at all when combining two
 *  expairseqs because due to the presorted nature the behaviour would be
 *  O(n) anyways, the code is quite messy, etc, etc.  The code is here as
 *  an example for following generations to tinker with. */
#define EXPAIRSEQ_USE_HASHTAB 0

typedef std::vector<expair> epvector;       ///< expair-vector
typedef epvector::iterator epp;             ///< expair-vector pointer
typedef std::list<epp> epplist;             ///< list of expair-vector pointers
typedef std::vector<epplist> epplistvector; ///< vector of epplist

/** Complex conjugate every element of an epvector. Returns zero if this
 *  does not change anything. */
epvector* conjugateepvector(const epvector&);

/** A sequence of class expair.
 *  This is used for time-critical classes like sums and products of terms
 *  since handling a list of coeff and rest is much faster than handling a
 *  list of products or powers, respectively. (Not incidentally, Maple does it
 *  the same way, maybe others too.)  The semantics is (at least) twofold:
 *  one for addition and one for multiplication and several methods have to
 *  be overridden by derived classes to reflect the change in semantics.
 *  However, most functionality turns out to be shared between addition and
 *  multiplication, which is the reason why there is this base class. */
class expairseq : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(expairseq, basic)

	friend class print_order;
        friend class ex;
	// other constructors
public:
	expairseq(const ex & lh, const ex & rh);
	expairseq(const exvector & v);
	expairseq(const epvector & v, const numeric& oc, bool do_index_renaming = false);
	expairseq(std::unique_ptr<epvector>, const numeric& oc, bool do_index_renaming = false);
	
	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const override {return 10;}
	bool info(unsigned inf) const override;
	size_t nops() const override;
	const ex op(size_t i) const override;
	virtual ex stable_op(size_t i) const;
        void set_pair_from(size_t i, ex e) { seq[i] = split_ex_to_pair(e); }
	const numeric & get_overall_coeff() const { return overall_coeff; }
	ex map(map_function & f) const override;
	ex eval(int level=0) const override;
	ex to_rational(exmap & repl) const override;
	ex to_polynomial(exmap & repl) const override;
	bool match(const ex & pattern, exmap& map) const override;
	ex subs(const exmap & m, unsigned options = 0) const override;
	ex conjugate() const override;
	numeric calc_total_degree() const;
	virtual const epvector & get_sorted_seq() const;
protected:
	bool is_equal_same_type(const basic & other) const override;
	unsigned return_type() const override;
	long calchash() const override;
	ex expand(unsigned options=0) const override;
	
	// new virtual functions which can be overridden by derived classes
protected:
	virtual ex thisexpairseq(const epvector & v, const numeric& oc, bool do_index_renaming = false) const;
	virtual ex thisexpairseq(std::unique_ptr<epvector> vp, const numeric& oc, bool do_index_renaming = false) const;
	virtual void printseq(const print_context & c, char delim,
	                      unsigned this_precedence,
	                      unsigned upper_precedence) const;
	virtual void printpair(const print_context & c, const expair & p,
	                       unsigned upper_precedence) const;
	virtual expair split_ex_to_pair(const ex & e) const;
	virtual expair combine_ex_with_coeff_to_pair(const ex & e,
	                                             const numeric & c) const;
	virtual expair combine_pair_with_coeff_to_pair(const expair & p,
	                                               const numeric & c) const;
	virtual ex recombine_pair_to_ex(const expair & p) const;
	virtual bool expair_needs_further_processing(epp it);
	virtual numeric default_overall_coeff() const;
	virtual void combine_overall_coeff(const numeric & c);
	virtual void combine_overall_coeff(const numeric & c1, const numeric & c2);
	virtual bool can_make_flat(const expair & p) const;
	
	// non-virtual functions in this class
protected:
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_tree(const print_tree & c, unsigned level) const override;
	void construct_from_2_ex_via_exvector(const ex & lh, const ex & rh);
	void construct_from_2_ex(const ex & lh, const ex & rh);
	void construct_from_2_expairseq(const expairseq & s1,
	                                const expairseq & s2);
	void construct_from_expairseq_ex(const expairseq & s,
	                                 const ex & e);
	void construct_from_exvector(const exvector & v, bool hold=false);
	void construct_from_epvector(const epvector & v, bool do_index_renaming = false);
	void make_flat(const exvector & v, bool hold=false);
	void make_flat(const epvector & v, bool do_index_renaming = false);
	void canonicalize();
	void combine_same_terms_sorted_seq();
        bool overall_coeff_equals_default() const;
#if EXPAIRSEQ_USE_HASHTAB
	void combine_same_terms();
	unsigned calc_hashtabsize(unsigned sz) const;
	unsigned calc_hashindex(const ex & e) const;
	void shrink_hashtab();
	void remove_hashtab_entry(epvector::const_iterator element);
	void move_hashtab_entry(epvector::const_iterator oldpos,
	                        epvector::iterator newpos);
	void sorted_insert(epplist & eppl, epvector::const_iterator elem);
	void build_hashtab_and_combine(epvector::iterator & first_numeric,
	                               epvector::iterator & last_non_zero,
	                               vector<bool> & touched,
	                               unsigned & number_of_zeroes);
	void drop_coeff_0_terms(epvector::iterator & first_numeric,
	                        epvector::iterator & last_non_zero,
	                        vector<bool> & touched,
	                        unsigned & number_of_zeroes);
	bool has_coeff_0() const;
	void add_numerics_to_hashtab(epvector::iterator first_numeric,
	                             epvector::const_iterator last_non_zero);
#endif // EXPAIRSEQ_USE_HASHTAB
	bool is_canonical() const;
	std::unique_ptr<epvector> expandchildren(unsigned options) const;
	std::unique_ptr<epvector> evalchildren(int level) const;
	std::unique_ptr<epvector> subschildren(const exmap & m, unsigned options = 0) const;
	
// member variables
	
protected:
	epvector seq;
	mutable epvector seq_sorted;
	numeric overall_coeff;
#if EXPAIRSEQ_USE_HASHTAB
	epplistvector hashtab;
	unsigned hashtabsize;
	unsigned hashmask;
	static unsigned maxhashtabsize;
	static unsigned minhashtabsize;
	static unsigned hashtabfactor;
#endif // EXPAIRSEQ_USE_HASHTAB
};

/** Class to handle the renaming of dummy indices. It holds a vector of
 *  indices that are being used in the expression so-far. If the same
 *  index occurs again as a dummy index in a factor, it is to be renamed.
 *  Unless dummy index renaming was swichted of, of course ;-) . */
class make_flat_inserter
{
	public:
		make_flat_inserter(const epvector &epv, bool b)
		{}
		make_flat_inserter(const exvector &v, bool b)
                {}
		ex handle_factor(const ex &x, const ex &coeff)
		{
			if (is_exactly_a<numeric>(coeff) and coeff.is_zero())
				return coeff;
		        return x;
		}
};

} // namespace GiNaC

#endif // ndef __GINAC_EXPAIRSEQ_H__
