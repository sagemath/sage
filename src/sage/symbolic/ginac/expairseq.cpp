/** @file expairseq.cpp
 *
 *  Implementation of sequences of expression pairs. */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *  (C) 2015-2018 Ralf Stephan <ralf@ark.in-berlin.de>
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

#include "expairseq.h"
#include "lst.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "relational.h"
#include "wildcard.h"
#include "archive.h"
#include "operators.h"
#include "utils.h"
#include "infinity.h"
#include "compiler.h"
#include "cmatcher.h"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <stdexcept>

#if EXPAIRSEQ_USE_HASHTAB
#include <cmath>
#endif // EXPAIRSEQ_USE_HASHTAB

namespace GiNaC {

	
GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(expairseq, basic,
  print_func<print_context>(&expairseq::do_print).
  print_func<print_tree>(&expairseq::do_print_tree))


//////////
// helper classes
//////////

class epp_is_less
{
public:
	bool operator()(const epp &lh, const epp &rh) const
	{
		return (*lh).is_less(*rh);
	}
};

//////////
// default constructor
//////////

// public

expairseq::expairseq() : inherited(&expairseq::tinfo_static)
#if EXPAIRSEQ_USE_HASHTAB
                                                   , hashtabsize(0)
#endif // EXPAIRSEQ_USE_HASHTAB
{}

// protected

#if 0
/** For use by copy ctor and assignment operator. */
void expairseq::copy(const expairseq &other)
{
	seq = other.seq;
	overall_coeff = other.overall_coeff;
#if EXPAIRSEQ_USE_HASHTAB
	// copy hashtab
	hashtabsize = other.hashtabsize;
	if (hashtabsize!=0) {
		hashmask = other.hashmask;
		hashtab.resize(hashtabsize);
		epvector::const_iterator osb = other.seq.begin();
		for (unsigned i=0; i<hashtabsize; ++i) {
			hashtab[i].clear();
			for (epplist::const_iterator cit=other.hashtab[i].begin();
			     cit!=other.hashtab[i].end(); ++cit) {
				hashtab[i].push_back(seq.begin()+((*cit)-osb));
			}
		}
	} else {
		hashtab.clear();
	}
#endif // EXPAIRSEQ_USE_HASHTAB
}
#endif

//////////
// other constructors
//////////

expairseq::expairseq(const ex &lh, const ex &rh) : inherited(&expairseq::tinfo_static)
{
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(const exvector &v) : inherited(&expairseq::tinfo_static)
{
	construct_from_exvector(v);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(const epvector &v,
                const numeric& oc, bool do_index_renaming)
  : inherited(&expairseq::tinfo_static), overall_coeff(oc)
{
	construct_from_epvector(v, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

expairseq::expairseq(std::unique_ptr<epvector> vp,
                const numeric& oc, bool do_index_renaming)
  : inherited(&expairseq::tinfo_static), overall_coeff(oc)
{
	GINAC_ASSERT(vp.get()!=0);
	construct_from_epvector(*vp, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

//////////
// archiving
//////////

expairseq::expairseq(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
#if EXPAIRSEQ_USE_HASHTAB
	, hashtabsize(0)
#endif
{
	auto first = n.find_first("rest");
	auto last = n.find_last("coeff");
	++last;
	seq.reserve((last-first)/2);

	for (auto loc = first; loc < last;) {
		ex lrest;
		ex lcoeff;
		n.find_ex_by_loc(loc++, lrest, sym_lst);
		n.find_ex_by_loc(loc++, lcoeff, sym_lst);
		seq.emplace_back(lrest, lcoeff);
	}

        ex oc;
	n.find_ex("overall_coeff", oc, sym_lst);
        overall_coeff = ex_to<numeric>(oc);

	canonicalize();
	GINAC_ASSERT(is_canonical());
}

void expairseq::archive(archive_node &n) const
{
	inherited::archive(n);
	for (const auto & elem : seq) {
		n.add_ex("rest", elem.rest);
		n.add_ex("coeff", elem.coeff);
	}
	n.add_ex("overall_coeff", ex(overall_coeff));
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

void expairseq::do_print(const print_context & c, unsigned level) const
{
	c.s << "[[";
	printseq(c, ',', precedence(), level);
	c.s << "]]";
}

void expairseq::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", nops=" << nops()
	    << std::endl;
	size_t num = seq.size();
	for (size_t i=0; i<num; ++i) {
		seq[i].rest.print(c, level + c.delta_indent);
		seq[i].coeff.print(c, level + c.delta_indent);
		if (i != num - 1)
			c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl;
	}
	if (!overall_coeff_equals_default()) {
		c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl
		    << std::string(level + c.delta_indent, ' ') << "overall_coeff" << std::endl;
		overall_coeff.print(c, level + c.delta_indent);
	}
	c.s << std::string(level + c.delta_indent,' ') << "=====" << std::endl;
#if EXPAIRSEQ_USE_HASHTAB
	c.s << std::string(level + c.delta_indent,' ')
	    << "hashtab size " << hashtabsize << std::endl;
	if (hashtabsize == 0) return;
#define MAXCOUNT 5
	unsigned count[MAXCOUNT+1];
	for (int i=0; i<MAXCOUNT+1; ++i)
		count[i] = 0;
	unsigned this_bin_fill;
	unsigned cum_fill_sq = 0;
	unsigned cum_fill = 0;
	for (unsigned i=0; i<hashtabsize; ++i) {
		this_bin_fill = 0;
		if (hashtab[i].size() > 0) {
			c.s << std::string(level + c.delta_indent, ' ')
			    << "bin " << i << " with entries ";
			for (epplist::const_iterator it=hashtab[i].begin();
			     it!=hashtab[i].end(); ++it) {
				c.s << *it-seq.begin() << " ";
				++this_bin_fill;
			}
			c.s << std::endl;
			cum_fill += this_bin_fill;
			cum_fill_sq += this_bin_fill*this_bin_fill;
		}
		if (this_bin_fill<MAXCOUNT)
			++count[this_bin_fill];
		else
			++count[MAXCOUNT];
	}
	unsigned fact = 1;
	double cum_prob = 0;
	double lambda = (1.0*seq.size()) / hashtabsize;
	for (int k=0; k<MAXCOUNT; ++k) {
		if (k>0)
			fact *= k;
		double prob = std::pow(lambda,k)/fact * std::exp(-lambda);
		cum_prob += prob;
		c.s << std::string(level + c.delta_indent, ' ') << "bins with " << k << " entries: "
		    << int(1000.0*count[k]/hashtabsize)/10.0 << "% (expected: "
		    << int(prob*1000)/10.0 << ")" << std::endl;
	}
	c.s << std::string(level + c.delta_indent, ' ') << "bins with more entries: "
	    << int(1000.0*count[MAXCOUNT]/hashtabsize)/10.0 << "% (expected: "
	    << int((1-cum_prob)*1000)/10.0 << ")" << std::endl;

	c.s << std::string(level + c.delta_indent, ' ') << "variance: "
	    << 1.0/hashtabsize*cum_fill_sq-(1.0/hashtabsize*cum_fill)*(1.0/hashtabsize*cum_fill)
	    << std::endl;
	c.s << std::string(level + c.delta_indent, ' ') << "average fill: "
	    << (1.0*cum_fill)/hashtabsize
	    << " (should be equal to " << (1.0*seq.size())/hashtabsize << ")" << std::endl;
#endif // EXPAIRSEQ_USE_HASHTAB
}

bool expairseq::info(unsigned inf) const
{
	switch(inf) {
		case info_flags::expanded:
			return (flags & status_flags::expanded) != 0u;
	}
	return inherited::info(inf);
}

size_t expairseq::nops() const
{
	if (overall_coeff_equals_default())
		return seq.size();
        return seq.size()+1;
}

const ex expairseq::op(size_t i) const
{
	if (i < seq.size())
		return recombine_pair_to_ex(seq[i]);
	GINAC_ASSERT(!overall_coeff_equals_default());
	return overall_coeff;
}

ex expairseq::stable_op(size_t i) const
{
	if (i < seq.size()) {
		return recombine_pair_to_ex(get_sorted_seq()[i]);
	}
	GINAC_ASSERT(!overall_coeff_equals_default());
	return overall_coeff;
}

ex expairseq::map(map_function &f) const
{
	std::unique_ptr<epvector> v(new epvector);
	v->reserve(seq.size()+1);

	for (const auto & elem : seq)
		v->push_back(split_ex_to_pair(f(recombine_pair_to_ex(elem))));

	if (overall_coeff_equals_default())
		return thisexpairseq(std::move(v), default_overall_coeff(), true);
        ex newcoeff = f(overall_coeff);
        if(is_exactly_a<numeric>(newcoeff))
                return thisexpairseq(std::move(v),
                                ex_to<numeric>(newcoeff), true);

        v->push_back(split_ex_to_pair(newcoeff));
        return thisexpairseq(std::move(v), default_overall_coeff(), true);
}

/** Perform coefficient-wise automatic term rewriting rules in this class. */
ex expairseq::eval(int level) const
{
	if ((level==1) and is_evaluated())
		return *this;
	
	std::unique_ptr<epvector> vp = evalchildren(level);
	if (vp == nullptr)
		return this->hold();
	
	return (new expairseq(std::move(vp), overall_coeff))->setflag(status_flags::dynallocated | status_flags::evaluated);
}

epvector* conjugateepvector(const epvector& epv)
{
	epvector *newepv = nullptr;
	for (const auto & elem : epv) {
		if (newepv != nullptr) {
			newepv->push_back(elem.conjugate());
			continue;
		}
		expair x = elem.conjugate();
		if (x.is_equal(elem)) {
			continue;
		}
		newepv = new epvector;
		newepv->reserve(epv.size());
                for (const auto & elem2 : epv)
                        if (&elem2 != &elem)
                                newepv->push_back(elem2);
                        else
                                break;
		newepv->push_back(x);
	}
	return newepv;
}

ex expairseq::conjugate() const
{
	std::unique_ptr<epvector> newepv(conjugateepvector(seq));
        // this is known to be numeric
	const numeric& x = ex_to<numeric>(overall_coeff.conjugate());
	if (!newepv and overall_coeff.is_equal(x)) {
		return *this;
	}
	ex result = thisexpairseq(newepv ? *newepv : seq, x);
	return result;
}

// returns total degree of this sequence
numeric expairseq::calc_total_degree() const
{
	numeric deg = 0;
	for (const auto & elem : seq)
		deg = deg.add(ex_to<numeric>(elem.coeff));
	return deg;
}

// protected
int expairseq::compare_same_type(const basic &other) const
{
	GINAC_ASSERT(is_a<expairseq>(other));
	const expairseq &o = static_cast<const expairseq &>(other);

	int cmpval;

	// compare number of elements
	if (seq.size() != o.seq.size())
		return (seq.size()<o.seq.size()) ? -1 : 1;

	// compare overall_coeff
	cmpval = overall_coeff.compare_same_type(o.overall_coeff);
	if (cmpval!=0)
		return cmpval;

#if EXPAIRSEQ_USE_HASHTAB
	GINAC_ASSERT(hashtabsize==o.hashtabsize);
	if (hashtabsize==0) {
#endif // EXPAIRSEQ_USE_HASHTAB

		auto cit1 = seq.begin();
		auto cit2 = o.seq.begin();
		auto last1 = seq.end();
		auto last2 = o.seq.end();

		for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
			cmpval = (*cit1).compare(*cit2);
			if (cmpval!=0) return cmpval;
		}

		GINAC_ASSERT(cit1==last1);
		GINAC_ASSERT(cit2==last2);

		return 0;

#if EXPAIRSEQ_USE_HASHTAB
	}

	// compare number of elements in each hashtab entry
	for (unsigned i=0; i<hashtabsize; ++i) {
		unsigned cursize=hashtab[i].size();
		if (cursize != o.hashtab[i].size())
			return (cursize < o.hashtab[i].size()) ? -1 : 1;
	}

	// compare individual (sorted) hashtab entries
	for (unsigned i=0; i<hashtabsize; ++i) {
		unsigned sz = hashtab[i].size();
		if (sz>0) {
			const epplist &eppl1 = hashtab[i];
			const epplist &eppl2 = o.hashtab[i];
			epplist::const_iterator it1 = eppl1.begin();
			epplist::const_iterator it2 = eppl2.begin();
			while (it1!=eppl1.end()) {
				cmpval = (*(*it1)).compare(*(*it2));
				if (cmpval!=0)
					return cmpval;
				++it1;
				++it2;
			}
		}
	}
	return 0; // equal
#endif // EXPAIRSEQ_USE_HASHTAB
}

bool expairseq::is_equal_same_type(const basic &other) const
{
	const expairseq &o = static_cast<const expairseq &>(other);
	
	// compare number of elements
	if (seq.size()!=o.seq.size())
		return false;
	
	// compare overall_coeff
	if (!overall_coeff.is_equal(o.overall_coeff))
		return false;
	
#if EXPAIRSEQ_USE_HASHTAB
	// compare number of elements in each hashtab entry
	if (hashtabsize!=o.hashtabsize) {
		print(print_tree(std::cout));
		other.print(print_tree(std::cout));
	}
		
	GINAC_ASSERT(hashtabsize==o.hashtabsize);
	
	if (hashtabsize==0) {
#endif // EXPAIRSEQ_USE_HASHTAB
		auto cit1 = seq.begin();
		auto cit2 = o.seq.begin();
		auto last1 = seq.end();
		
		while (cit1!=last1) {
			if (!(*cit1).is_equal(*cit2)) return false;
			++cit1;
			++cit2;
		}
                
                return true;

#if EXPAIRSEQ_USE_HASHTAB
	}
	
	for (unsigned i=0; i<hashtabsize; ++i) {
		if (hashtab[i].size() != o.hashtab[i].size())
			return false;
	}

	// compare individual sorted hashtab entries
	for (unsigned i=0; i<hashtabsize; ++i) {
		unsigned sz = hashtab[i].size();
		if (sz>0) {
			const epplist &eppl1 = hashtab[i];
			const epplist &eppl2 = o.hashtab[i];
			epplist::const_iterator it1 = eppl1.begin();
			epplist::const_iterator it2 = eppl2.begin();
			while (it1!=eppl1.end()) {
				if (!(*(*it1)).is_equal(*(*it2))) return false;
				++it1;
				++it2;
			}
		}
	}
	
	return true;
#endif // EXPAIRSEQ_USE_HASHTAB
}

unsigned expairseq::return_type() const
{
	return return_types::noncommutative_composite;
}

long expairseq::calchash() const
{
	long v = golden_ratio_hash((intptr_t)tinfo());
        for (const auto & elem : seq) {
		v ^= elem.rest.gethash();
#if !EXPAIRSEQ_USE_HASHTAB
		// rotation spoils commutativity!
		v = rotate_left(v);
		v ^= elem.coeff.gethash();
#endif // !EXPAIRSEQ_USE_HASHTAB
	}

	v ^= overall_coeff.gethash();

	// store calculated hash value only if object is already evaluated
	if (is_evaluated()) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}
	
	return v;
}

ex expairseq::expand(unsigned options) const
{
	std::unique_ptr<epvector> vp = expandchildren(options);
	if (vp != nullptr)
		return thisexpairseq(std::move(vp), overall_coeff);

        // The terms have not changed, so it is safe to declare this expanded
        return (options == 0) ? setflag(status_flags::expanded) : *this;
}

ex expairseq::subs(const exmap & m, unsigned options) const
{
	std::unique_ptr<epvector> vp = subschildren(m, options);
	if (vp != nullptr) {
                const ex& ocs = overall_coeff.subs(m, options);
                if (is_exactly_a<numeric>(ocs))
                        return ex_to<basic>(thisexpairseq(std::move(vp),
                                            ex_to<numeric>(ocs),
                                            (options & subs_options::no_index_renaming) == 0));
                else
                        return ex_to<basic>(add(ocs, thisexpairseq(std::move(vp),
                                            *_num0_p,
                                            (options & subs_options::no_index_renaming) == 0)));
        }
	if (((options & subs_options::algebraic) != 0u) && is_exactly_a<mul>(*this))
		return static_cast<const mul *>(this)->algebraic_subs_mul(m, options);
	else
		return subs_one_level(m, options);
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// protected

/** Create an object of this type.
 *  This method works similar to a constructor.  It is useful because expairseq
 *  has (at least) two possible different semantics but we want to inherit
 *  methods thus avoiding code duplication.  Sometimes a method in expairseq
 *  has to create a new one of the same semantics, which cannot be done by a
 *  ctor because the name (add, mul,...) is unknown on the expaiseq level.  In
 *  order for this trick to work a derived class must of course override this
 *  definition. */
ex expairseq::thisexpairseq(const epvector &v, const numeric &oc, bool do_index_renaming) const
{
	return expairseq(v, oc, do_index_renaming);
}

ex expairseq::thisexpairseq(std::unique_ptr<epvector> vp, const numeric &oc, bool do_index_renaming) const
{
	return expairseq(std::move(vp), oc, do_index_renaming);
}

void expairseq::printpair(const print_context & c, const expair & p, unsigned /*unused*/) const
{
	c.s << "[[";
	p.rest.print(c, precedence());
	c.s << ",";
	p.coeff.print(c, precedence());
	c.s << "]]";
}

void expairseq::printseq(const print_context & c, char delim,
                         unsigned this_precedence,
                         unsigned upper_precedence) const
{
	if (this_precedence <= upper_precedence)
		c.s << "(";
        auto it = seq.begin();
        auto it_last = seq.end() - 1;
        for (; it != it_last; ++it) {
		printpair(c, *it, this_precedence);
		c.s << delim;
	}
	printpair(c, *it, this_precedence);
	if (!overall_coeff_equals_default()) {
		c.s << delim;
		overall_coeff.print(c, this_precedence);
	}
	
	if (this_precedence <= upper_precedence)
		c.s << ")";
}


/** Form an expair from an ex, using the corresponding semantics.
 *  @see expairseq::recombine_pair_to_ex() */
expair expairseq::split_ex_to_pair(const ex &e) const
{
	return expair(e,_ex1);
}


expair expairseq::combine_ex_with_coeff_to_pair(const ex &e,
                                                const numeric &c) const
{
	return expair(e,c);
}


expair expairseq::combine_pair_with_coeff_to_pair(const expair &p,
                                                  const numeric &c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));
	
	return expair(p.rest,ex_to<numeric>(p.coeff).mul_dyn(c));
}


/** Form an ex out of an expair, using the corresponding semantics.
 *  @see expairseq::split_ex_to_pair() */
ex expairseq::recombine_pair_to_ex(const expair &p) const
{
	return lst(p.rest,p.coeff);
}

bool expairseq::expair_needs_further_processing(epp /*unused*/)
{
#if EXPAIRSEQ_USE_HASHTAB
	//#  error "FIXME: expair_needs_further_processing not yet implemented for hashtabs, sorry. A.F."
#endif // EXPAIRSEQ_USE_HASHTAB
	return false;
}

numeric expairseq::default_overall_coeff() const
{
	return *_num0_p;
}

bool expairseq::overall_coeff_equals_default() const
{
    return (overall_coeff.is_exact()
         and overall_coeff.is_equal(default_overall_coeff()));
}

void expairseq::combine_overall_coeff(const numeric &c)
{
         overall_coeff += c;
}

void expairseq::combine_overall_coeff(const numeric &c1, const numeric &c2)
{
	overall_coeff += c1.mul(c2);
}

bool expairseq::can_make_flat(const expair & /*unused*/) const
{
	return true;
}


//////////
// non-virtual functions in this class
//////////

void expairseq::construct_from_2_ex_via_exvector(const ex &lh, const ex &rh)
{
	exvector v;
	v.reserve(2);
	v.push_back(lh);
	v.push_back(rh);
	construct_from_exvector(v);
#if EXPAIRSEQ_USE_HASHTAB
	GINAC_ASSERT((hashtabsize==0)||(hashtabsize>=minhashtabsize));
	GINAC_ASSERT(hashtabsize==calc_hashtabsize(seq.size()));
#endif // EXPAIRSEQ_USE_HASHTAB
}

void expairseq::construct_from_2_ex(const ex &lh, const ex &rh)
{
	if (ex_to<basic>(lh).tinfo()==this->tinfo()) {
		if (ex_to<basic>(rh).tinfo()==this->tinfo()) {
#if EXPAIRSEQ_USE_HASHTAB
			unsigned totalsize = ex_to<expairseq>(lh).seq.size() +
			                     ex_to<expairseq>(rh).seq.size();
                        construct_from_2_ex_via_exvector(lh,rh);
                        {
                } else {
#endif // EXPAIRSEQ_USE_HASHTAB
                        construct_from_2_expairseq(ex_to<expairseq>(lh),
                                                   ex_to<expairseq>(rh));
#if EXPAIRSEQ_USE_HASHTAB
			}
#endif // EXPAIRSEQ_USE_HASHTAB
			return;
		} 
#if EXPAIRSEQ_USE_HASHTAB
                unsigned totalsize = ex_to<expairseq>(lh).seq.size()+1;
                if (calc_hashtabsize(totalsize)!=0) {
                        construct_from_2_ex_via_exvector(lh, rh);
                } else {
#endif // EXPAIRSEQ_USE_HASHTAB
                        construct_from_expairseq_ex(ex_to<expairseq>(lh), rh);
#if EXPAIRSEQ_USE_HASHTAB
			}
#endif // EXPAIRSEQ_USE_HASHTAB
			return;
		
	} else if (ex_to<basic>(rh).tinfo()==this->tinfo()) {
#if EXPAIRSEQ_USE_HASHTAB
		unsigned totalsize=ex_to<expairseq>(rh).seq.size()+1;
		if (calc_hashtabsize(totalsize)!=0) {
			construct_from_2_ex_via_exvector(lh,rh);
		} else {
#endif // EXPAIRSEQ_USE_HASHTAB
			construct_from_expairseq_ex(ex_to<expairseq>(rh),lh);
#if EXPAIRSEQ_USE_HASHTAB
		}
#endif // EXPAIRSEQ_USE_HASHTAB
		return;
	}
	
#if EXPAIRSEQ_USE_HASHTAB
	if (calc_hashtabsize(2)!=0) {
		construct_from_2_ex_via_exvector(lh,rh);
		return;
	}
	hashtabsize = 0;
#endif // EXPAIRSEQ_USE_HASHTAB
	
	if (is_exactly_a<numeric>(lh)) {
                const numeric& lhn = ex_to<numeric>(lh);
		if (is_exactly_a<numeric>(rh)) {
			combine_overall_coeff(lhn);
			combine_overall_coeff(ex_to<numeric>(rh));
		} else {
			combine_overall_coeff(lhn);
			seq.push_back(split_ex_to_pair(rh));
		}
	} else {
		if (is_exactly_a<numeric>(rh)) {
			combine_overall_coeff(ex_to<numeric>(rh));
			seq.push_back(split_ex_to_pair(lh));
		} else {
			expair p1 = split_ex_to_pair(lh);
			expair p2 = split_ex_to_pair(rh);

			int cmpval = p1.rest.compare(p2.rest);
			if (cmpval==0
                            and likely(not is_exactly_a<infinity>(p1.rest))) {
				p1.coeff = ex_to<numeric>(p1.coeff).add_dyn(ex_to<numeric>(p2.coeff));
				if (!ex_to<numeric>(p1.coeff).is_zero()) {

					// no further processing is necessary, since this
					// one element will usually be recombined in eval()
					seq.push_back(p1);
				}
			} else {
				seq.reserve(2);
				if (cmpval<0) {
					seq.push_back(p1);
					seq.push_back(p2);
				} else {
					seq.push_back(p2);
					seq.push_back(p1);
				}
			}
		}
	}
}

void expairseq::construct_from_2_expairseq(const expairseq &s1,
		const expairseq &s2)
{
	combine_overall_coeff(s1.overall_coeff);
	combine_overall_coeff(s2.overall_coeff);

	auto first1 = s1.seq.begin();
	auto last1 = s1.seq.end();
	auto first2 = s2.seq.begin();
	auto last2 = s2.seq.end();

	seq.reserve(s1.seq.size()+s2.seq.size());

	bool needs_further_processing=false;
	
	while (first1!=last1 && first2!=last2) {
		int cmpval = (*first1).rest.compare((*first2).rest);

		if (cmpval==0) {
			// infinity evaluation is handled in the eval() method
			// do not let infinities cancel each other here
			if (unlikely(is_exactly_a<infinity>(first1->rest))) {
				seq.push_back(*first1);
				seq.push_back(*first2);
			} else {
				// combine terms
				const numeric &newcoeff = ex_to<numeric>(first1->coeff).
					add(ex_to<numeric>(first2->coeff));
				if (!newcoeff.is_zero()) {
					seq.emplace_back(first1->rest,newcoeff);
					if (expair_needs_further_processing(seq.end()-1)) {
						needs_further_processing = true;
					}
				}
			}
			++first1;
			++first2;
		} else if (cmpval<0) {
			seq.push_back(*first1);
			++first1;
		} else {
			seq.push_back(*first2);
			++first2;
		}
	}
	
	while (first1!=last1) {
		seq.push_back(*first1);
		++first1;
	}
	while (first2!=last2) {
		seq.push_back(*first2);
		++first2;
	}
	
	if (needs_further_processing) {
		epvector v = seq;
		seq.clear();
		construct_from_epvector(v);
	}
}

void expairseq::construct_from_expairseq_ex(const expairseq &s,
			const ex &e)
{
	combine_overall_coeff(s.overall_coeff);
	if (is_exactly_a<numeric>(e)) {
		combine_overall_coeff(ex_to<numeric>(e));
		seq = s.seq;
		return;
	}
	
	auto first = s.seq.begin();
	auto last = s.seq.end();
	expair p = split_ex_to_pair(e);
	// infinity evaluation is handled in eval()
	// do not let infinities cancel each other here
	if (unlikely(is_exactly_a<infinity>(p.rest))) {
		seq.push_back(p);
		seq.insert(seq.end(), first, last);
		return;
	}
	
	seq.reserve(s.seq.size()+1);
	bool p_pushed = false;
	
	bool needs_further_processing=false;
	
	// merge p into s.seq
	while (first!=last) {
		int cmpval = (*first).rest.compare(p.rest);
		if (cmpval==0) {
			// combine terms
			const numeric &newcoeff = ex_to<numeric>(first->coeff).
			                           add(ex_to<numeric>(p.coeff));
			if (!newcoeff.is_zero()) {
				seq.emplace_back(first->rest,newcoeff);
				if (expair_needs_further_processing(seq.end()-1))
					needs_further_processing = true;
			}
			++first;
			p_pushed = true;
			break;
		}
                if (cmpval<0) {
			seq.push_back(*first);
			++first;
		} else {
			seq.push_back(p);
			p_pushed = true;
			break;
		}
	}
	
	if (p_pushed) {
		// while loop exited because p was pushed, now push rest of s.seq
		while (first!=last) {
			seq.push_back(*first);
			++first;
		}
	} else {
		// while loop exited because s.seq was pushed, now push p
		seq.push_back(p);
	}

	if (needs_further_processing) {
		epvector v = seq;
		seq.clear();
		construct_from_epvector(v);
	}
}

void expairseq::construct_from_exvector(const exvector &v, bool do_hold)
{
	// simplifications: +(a,+(b,c),d) -> +(a,b,c,d) (associativity)
	//                  +(d,b,c,a) -> +(a,b,c,d) (canonicalization)
	//                  +(...,x,*(x,c1),*(x,c2)) -> +(...,*(x,1+c1+c2)) (c1, c2 numeric)
	//                  (same for (+,*) -> (*,^)

	make_flat(v, do_hold);
#if EXPAIRSEQ_USE_HASHTAB
	combine_same_terms();
#else
	if (!do_hold) {
		canonicalize();
		combine_same_terms_sorted_seq();
	}
#endif // EXPAIRSEQ_USE_HASHTAB
}

void expairseq::construct_from_epvector(const epvector &v, bool do_index_renaming)
{
	// simplifications: +(a,+(b,c),d) -> +(a,b,c,d) (associativity)
	//                  +(d,b,c,a) -> +(a,b,c,d) (canonicalization)
	//                  +(...,x,*(x,c1),*(x,c2)) -> +(...,*(x,1+c1+c2)) (c1, c2 numeric)
	//                  same for (+,*) -> (*,^)

	make_flat(v, do_index_renaming);
#if EXPAIRSEQ_USE_HASHTAB
	combine_same_terms();
#else
	canonicalize();
	combine_same_terms_sorted_seq();
#endif // EXPAIRSEQ_USE_HASHTAB
}

/** Combine this expairseq with argument exvector.
 *  It cares for associativity as well as for special handling of numerics. */
void expairseq::make_flat(const exvector &v, bool do_hold)
{
	// count number of operands which are of same expairseq derived type
	// and their cumulative number of operands
	int nexpairseqs = 0;
	int noperands = 0;
	bool do_idx_rename = false;
	
	if (!do_hold) {
                for (const auto & elem : v) {
			if (ex_to<basic>(elem).tinfo()==this->tinfo()) {
				++nexpairseqs;
				noperands += ex_to<expairseq>(elem).seq.size();
			}
		}
	} else
		this->hold();
	
	// reserve seq and coeffseq which will hold all operands
	seq.reserve(v.size()+noperands-nexpairseqs);
	
	// copy elements and split off numerical part
	make_flat_inserter mf(v, do_idx_rename);
        for (const auto & elem : v) {
		if (ex_to<basic>(elem).tinfo()==this->tinfo() && !do_hold) {
			ex newfactor = mf.handle_factor(elem, _ex1);
			const expairseq &subseqref = ex_to<expairseq>(newfactor);
			combine_overall_coeff(subseqref.overall_coeff);
                        for (const auto & elem2 : subseqref.seq)
				seq.push_back(elem2);
		} else {
			if (is_exactly_a<numeric>(elem))
				combine_overall_coeff(ex_to<numeric>(elem));
			else {
				ex newfactor = mf.handle_factor(elem, _ex1);
				seq.push_back(split_ex_to_pair(newfactor));
			}
		}
	}
}

/** Combine this expairseq with argument epvector.
 *  It cares for associativity as well as for special handling of numerics. */
void expairseq::make_flat(const epvector &v, bool do_index_renaming)
{
	// count number of operands which are of same expairseq derived type
	// and their cumulative number of operands
	int nexpairseqs = 0;
	int noperands = 0;
	
        for (const auto & elem : v) {
		if (ex_to<basic>(elem.rest).tinfo()==this->tinfo()) {
			++nexpairseqs;
			noperands += ex_to<expairseq>(elem.rest).seq.size();
		}
	}
	
	// reserve seq and coeffseq which will hold all operands
	seq.reserve(v.size()+noperands-nexpairseqs);
	make_flat_inserter mf(v, false);
	
	// copy elements and split off numerical part
        for (const auto & elem : v) {

		if (ex_to<basic>(elem.rest).tinfo()==this->tinfo() &&
		    this->can_make_flat(elem)) {
		        ex newrest;
		        if (is_exactly_a<numeric>(elem.coeff) and elem.coeff.is_zero()) {
		            newrest = default_overall_coeff();
		        } else {
		            newrest = elem.rest;
		        }
			const expairseq &subseqref = ex_to<expairseq>(newrest);
			combine_overall_coeff(subseqref.overall_coeff,
			                      ex_to<numeric>(elem.coeff));
                        for (const auto & elem2 : subseqref.seq)
				seq.emplace_back(elem2.rest,
                                        ex_to<numeric>(elem2.coeff).mul_dyn(ex_to<numeric>(elem.coeff)));
		} else {
			if (elem.is_canonical_numeric())
				combine_overall_coeff(ex_to<numeric>(mf.handle_factor(elem.rest, _ex1)));
			else {
			        if ((is_exactly_a<numeric>(elem.coeff) and elem.coeff.is_zero())
			                or (is_exactly_a<numeric>(elem.rest)
			                    and (ex_to<numeric>(elem.rest).is_equal(default_overall_coeff())))) {
			            combine_overall_coeff(default_overall_coeff());
			        }
			        else {
			                seq.push_back(elem);
			        }
			    }
		}
	}
}

/** Brings this expairseq into a sorted (canonical) form. */
void expairseq::canonicalize()
{
	std::sort(seq.begin(), seq.end(), expair_rest_is_less());
}


/** Compact a presorted expairseq by combining all matching expairs to one
 *  each.  On an add object, this is responsible for 2*x+3*x+y -> 5*x+y, for
 *  instance. */
void expairseq::combine_same_terms_sorted_seq()
{
	if (seq.size()<2)
		return;

	bool needs_further_processing = false;

	auto itin1 = seq.begin();
	auto itin2 = itin1+1;
	auto itout = itin1;
	auto last = seq.end();
	// must_copy will be set to true the first time some combination is 
	// possible from then on the sequence has changed and must be compacted
	bool must_copy = false;
	while (itin2!=last) {
		if (itin1->rest.compare(itin2->rest)==0
                    and likely(not is_exactly_a<infinity>(itin1->rest))) {
			itin1->coeff = ex_to<numeric>(itin1->coeff).
                                add_dyn(ex_to<numeric>(itin2->coeff));
			if (expair_needs_further_processing(itin1))
				needs_further_processing = true;
			must_copy = true;
		} else {
			if (not ex_to<numeric>(itin1->coeff).is_zero()) {
				if (must_copy)
					*itout = *itin1;
				++itout;
			}
			itin1 = itin2;
		}
		++itin2;
	}
	if (not ex_to<numeric>(itin1->coeff).is_zero()) {
		if (must_copy)
			*itout = *itin1;
		++itout;
	}
	if (itout!=last)
		seq.erase(itout,last);

	if (needs_further_processing) {
		epvector v = seq;
		seq.clear();
		construct_from_epvector(v);
	}
}

#if EXPAIRSEQ_USE_HASHTAB

unsigned expairseq::calc_hashtabsize(unsigned sz) const
{
	unsigned size;
	unsigned nearest_power_of_2 = 1 << log2(sz);
	// if (nearest_power_of_2 < maxhashtabsize/hashtabfactor) {
	//  size = nearest_power_of_2*hashtabfactor;
	size = nearest_power_of_2/hashtabfactor;
	if (size<minhashtabsize)
		return 0;

	// hashtabsize must be a power of 2
	GINAC_ASSERT((1U << log2(size))==size);
	return size;
}

unsigned expairseq::calc_hashindex(const ex &e) const
{
	// calculate hashindex
	unsigned hashindex;
	if (is_a<numeric>(e)) {
		hashindex = hashmask;
	} else {
		hashindex = e.gethash() & hashmask;
		// last hashtab entry is reserved for numerics
		if (hashindex==hashmask) hashindex = 0;
	}
	GINAC_ASSERT((hashindex<hashtabsize)||(hashtabsize==0));
	return hashindex;
}

void expairseq::shrink_hashtab()
{
	unsigned new_hashtabsize;
	while (hashtabsize!=(new_hashtabsize=calc_hashtabsize(seq.size()))) {
		GINAC_ASSERT(new_hashtabsize<hashtabsize);
		if (new_hashtabsize==0) {
			hashtab.clear();
			hashtabsize = 0;
			canonicalize();
			return;
		}
		
		// shrink by a factor of 2
		unsigned half_hashtabsize = hashtabsize/2;
		for (unsigned i=0; i<half_hashtabsize-1; ++i)
			hashtab[i].merge(hashtab[i+half_hashtabsize],epp_is_less());
		// special treatment for numeric hashes
		hashtab[0].merge(hashtab[half_hashtabsize-1],epp_is_less());
		hashtab[half_hashtabsize-1] = hashtab[hashtabsize-1];
		hashtab.resize(half_hashtabsize);
		hashtabsize = half_hashtabsize;
		hashmask = hashtabsize-1;
	}
}

void expairseq::remove_hashtab_entry(epvector::const_iterator element)
{
	if (hashtabsize==0)
		return; // nothing to do
	
	// calculate hashindex of element to be deleted
	unsigned hashindex = calc_hashindex((*element).rest);

	// find it in hashtab and remove it
	epplist &eppl = hashtab[hashindex];
	epplist::iterator epplit = eppl.begin();
	bool erased = false;
	while (epplit!=eppl.end()) {
		if (*epplit == element) {
			eppl.erase(epplit);
			erased = true;
			break;
		}
		++epplit;
	}
	if (!erased) {
		unsigned hashindex = calc_hashindex(element->rest);
		epplist &eppl = hashtab[hashindex];
		epplist::iterator epplit = eppl.begin();
		bool erased = false;
		while (epplit!=eppl.end()) {
			if (*epplit == element) {
				eppl.erase(epplit);
				erased = true;
				break;
			}
			++epplit;
		}
		GINAC_ASSERT(erased);
	}
	GINAC_ASSERT(erased);
}

void expairseq::move_hashtab_entry(epvector::const_iterator oldpos,
                                   epvector::iterator newpos)
{
	GINAC_ASSERT(hashtabsize!=0);
	
	// calculate hashindex of element which was moved
	unsigned hashindex=calc_hashindex((*newpos).rest);

	// find it in hashtab and modify it
	epplist &eppl = hashtab[hashindex];
	epplist::iterator epplit = eppl.begin();
	while (epplit!=eppl.end()) {
		if (*epplit == oldpos) {
			*epplit = newpos;
			break;
		}
		++epplit;
	}
	GINAC_ASSERT(epplit!=eppl.end());
}

void expairseq::sorted_insert(epplist &eppl, epvector::const_iterator elem)
{
	epplist::const_iterator current = eppl.begin();
	while ((current!=eppl.end()) && ((*current)->is_less(*elem))) {
		++current;
	}
	eppl.insert(current,elem);
}    

void expairseq::build_hashtab_and_combine(epvector::iterator &first_numeric,
                                          epvector::iterator &last_non_zero,
                                          std::vector<bool> &touched,
                                          unsigned &number_of_zeroes)
{
	epp current = seq.begin();

	while (current!=first_numeric) {
		if (is_exactly_a<numeric>(current->rest)) {
			--first_numeric;
			iter_swap(current,first_numeric);
		} else {
			// calculate hashindex
			unsigned currenthashindex = calc_hashindex(current->rest);

			// test if there is already a matching expair in the hashtab-list
			epplist &eppl=hashtab[currenthashindex];
			epplist::iterator epplit = eppl.begin();
			while (epplit!=eppl.end()) {
				if (current->rest.is_equal((*epplit)->rest))
					break;
				++epplit;
			}
			if (epplit==eppl.end()) {
				// no matching expair found, append this to end of list
				sorted_insert(eppl,current);
				++current;
			} else {
				// epplit points to a matching expair, combine it with current
				(*epplit)->coeff = ex_to<numeric>((*epplit)->coeff).
				                   add_dyn(ex_to<numeric>(current->coeff));
				
				// move obsolete current expair to end by swapping with last_non_zero element
				// if this was a numeric, it is swapped with the expair before first_numeric 
				iter_swap(current,last_non_zero);
				--first_numeric;
				if (first_numeric!=last_non_zero) iter_swap(first_numeric,current);
				--last_non_zero;
				++number_of_zeroes;
				// test if combined term has coeff 0 and can be removed is done later
				touched[(*epplit)-seq.begin()] = true;
			}
		}
	}
}

void expairseq::drop_coeff_0_terms(epvector::iterator &first_numeric,
                                   epvector::iterator &last_non_zero,
                                   std::vector<bool> &touched,
                                   unsigned &number_of_zeroes)
{
	// move terms with coeff 0 to end and remove them from hashtab
	// check only those elements which have been touched
	epp current = seq.begin();
	size_t i = 0;
	while (current!=first_numeric) {
		if (!touched[i]) {
			++current;
			++i;
		} else if (!ex_to<numeric>((*current).coeff).is_zero()) {
			++current;
			++i;
		} else {
			remove_hashtab_entry(current);
			
			// move element to the end, unless it is already at the end
			if (current!=last_non_zero) {
				iter_swap(current,last_non_zero);
				--first_numeric;
				bool numeric_swapped = first_numeric!=last_non_zero;
				if (numeric_swapped)
					iter_swap(first_numeric,current);
				epvector::iterator changed_entry;

				if (numeric_swapped)
					changed_entry = first_numeric;
				else
					changed_entry = last_non_zero;
				
				--last_non_zero;
				++number_of_zeroes;

				if (first_numeric!=current) {

					// change entry in hashtab which referred to first_numeric or last_non_zero to current
					move_hashtab_entry(changed_entry,current);
					touched[current-seq.begin()] = touched[changed_entry-seq.begin()];
				}
			} else {
				--first_numeric;
				--last_non_zero;
				++number_of_zeroes;
			}
		}
	}
	GINAC_ASSERT(i==current-seq.begin());
}

/** True if one of the coeffs vanishes, otherwise false.
 *  This would be an invariant violation, so this should only be used for
 *  debugging purposes. */
bool expairseq::has_coeff_0() const
{
	epvector::const_iterator i = seq.begin(), end = seq.end();
	while (i != end) {
		if (i->coeff.is_zero())
			return true;
		++i;
	}
	return false;
}

void expairseq::add_numerics_to_hashtab(epvector::iterator first_numeric,
                                        epvector::const_iterator last_non_zero)
{
	if (first_numeric == seq.end()) return; // no numerics
	
	epvector::const_iterator current = first_numeric, last = last_non_zero + 1;
	while (current != last) {
		sorted_insert(hashtab[hashmask], current);
		++current;
	}
}

void expairseq::combine_same_terms()
{
	// combine same terms, drop term with coeff 0, move numerics to end
	
	// calculate size of hashtab
	hashtabsize = calc_hashtabsize(seq.size());
	
	// hashtabsize is a power of 2
	hashmask = hashtabsize-1;
	
	// allocate hashtab
	hashtab.clear();
	hashtab.resize(hashtabsize);
	
	if (hashtabsize==0) {
		canonicalize();
		combine_same_terms_sorted_seq();
		GINAC_ASSERT(!has_coeff_0());
		return;
	}
	
	// iterate through seq, move numerics to end,
	// fill hashtab and combine same terms
	epvector::iterator first_numeric = seq.end();
	epvector::iterator last_non_zero = seq.end()-1;
	
	size_t num = seq.size();
	std::vector<bool> touched(num);
	
	unsigned number_of_zeroes = 0;
	
	GINAC_ASSERT(!has_coeff_0());
	build_hashtab_and_combine(first_numeric,last_non_zero,touched,number_of_zeroes);
	
	// there should not be any terms with coeff 0 from the beginning,
	// so it should be safe to skip this step
	if (number_of_zeroes!=0) {
		drop_coeff_0_terms(first_numeric,last_non_zero,touched,number_of_zeroes);
	}
	
	add_numerics_to_hashtab(first_numeric,last_non_zero);
	
	// pop zero elements
	for (unsigned i=0; i<number_of_zeroes; ++i) {
		seq.pop_back();
	}
	
	// shrink hashtabsize to calculated value
	GINAC_ASSERT(!has_coeff_0());
	
	shrink_hashtab();
	
	GINAC_ASSERT(!has_coeff_0());
}

#endif // EXPAIRSEQ_USE_HASHTAB

/** Check if this expairseq is in sorted (canonical) form.  Useful mainly for
 *  debugging or in assertions since being sorted is an invariance. */
bool expairseq::is_canonical() const
{
	if (seq.size() <= 1)
		return true;
	
#if EXPAIRSEQ_USE_HASHTAB
	if (hashtabsize > 0) return 1; // not canoncalized
#endif // EXPAIRSEQ_USE_HASHTAB
	
	auto it = seq.begin(), itend = seq.end();
	auto it_last = it;
	for (++it; it!=itend; it_last=it, ++it) {
		if (!(it_last->is_less(*it) || it_last->is_equal(*it))) {
			if (!is_exactly_a<numeric>(it_last->rest) ||
				!is_exactly_a<numeric>(it->rest)) {
				// double test makes it easier to set a breakpoint...
				if (!is_exactly_a<numeric>(it_last->rest) ||
					!is_exactly_a<numeric>(it->rest)) {
					printpair(std::clog, *it_last, 0);
					std::clog << ">";
					printpair(std::clog, *it, 0);
					std::clog << "\n";
					std::clog << "pair1:" << std::endl;
					it_last->rest.print(print_tree(std::clog));
					it_last->coeff.print(print_tree(std::clog));
					std::clog << "pair2:" << std::endl;
					it->rest.print(print_tree(std::clog));
					it->coeff.print(print_tree(std::clog));
					return false;
				}
			}
		}
	}
	return true;
}


/** Member-wise expand the expairs in this sequence.
 *
 *  @see expairseq::expand()
 *  @return pointer to epvector containing expanded pairs or zero pointer,
 *  if no members were changed. */
std::unique_ptr<epvector> expairseq::expandchildren(unsigned options) const
{
	auto cit = seq.begin();
	auto last = seq.end();
	while (cit!=last) {
		const ex expanded_ex = cit->rest.expand(options);
		if (!are_ex_trivially_equal(cit->rest,expanded_ex)) {
			
			// something changed, copy seq, eval and return it
                        std::unique_ptr<epvector> s(new epvector);
			s->reserve(seq.size());
			
			// copy parts of seq which are known not to have changed
                        s->insert(s->begin(), seq.begin(), cit);

			// copy first changed element
			s->push_back(expair(expanded_ex, cit->coeff));
			++cit;

			// copy rest
			while (cit != last) {
				s->push_back(expair(cit->rest.expand(options),
				                         cit->coeff));
				++cit;
			}
			return s;
		}

		++cit;
	}
	
	return std::unique_ptr<epvector>(nullptr); // signalling nothing has changed
}


/** Member-wise evaluate the expairs in this sequence.
 *
 *  @see expairseq::eval()
 *  @return pointer to epvector containing evaluated pairs or zero pointer,
 *  if no members were changed. */
std::unique_ptr<epvector> expairseq::evalchildren(int level) const
{
	if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	auto last = seq.end();
	auto cit = seq.begin();
	while (cit!=last) {
                const ex evaled_rest = level==1 ? cit->rest : cit->rest.eval(level-1);
                const expair evaled_pair = combine_ex_with_coeff_to_pair(evaled_rest, ex_to<numeric>(cit->coeff));
                if (unlikely(!evaled_pair.is_equal(*cit))) {

			// something changed: copy seq, eval, and return it
                        std::unique_ptr<epvector> s(new epvector);
			s->reserve(seq.size());

			// copy parts of seq which are known not to have changed
                        s->insert(s->begin(), seq.begin(), cit);

			// copy first changed element
			s->push_back(evaled_pair);
			++cit;

			// copy rest
			while (cit != last) {
                                const ex evaled_rest1 = level==1 ? cit->rest : cit->rest.eval(level-1);
				s->push_back(combine_ex_with_coeff_to_pair(evaled_rest1,
				             ex_to<numeric>(cit->coeff)));
				++cit;
			}
			return s;
		}
		++cit;
	}
	
	return std::unique_ptr<epvector>(nullptr); // signalling nothing has changed
}

/** Member-wise substitute in this sequence.
 *
 *  @see expairseq::subs()
 *  @return pointer to epvector containing pairs after application of subs,
 *    or NULL pointer if no members were changed. */
std::unique_ptr<epvector> expairseq::subschildren(const exmap & m, unsigned options) const
{
	// When any of the objects to be substituted is a product or power
	// we have to recombine the pairs because the numeric coefficients may
	// be part of the search pattern.
	if ((options & (subs_options::pattern_is_product | subs_options::pattern_is_not_product)) == 0u) {

		// Search the list of substitutions and cache our findings
		for (const auto & elem : m) {
			if (is_exactly_a<mul>(elem.first) || is_exactly_a<power>(elem.first)) {
				options |= subs_options::pattern_is_product;
				break;
			}
		}
		if ((options & subs_options::pattern_is_product) == 0u)
			options |= subs_options::pattern_is_not_product;
	}

	if ((options & subs_options::pattern_is_product) != 0u) {

		// Substitute in the recombined pairs
		auto cit = seq.begin(), last = seq.end();
                while (cit != last) {

                        const ex &orig_ex = recombine_pair_to_ex(*cit);
                        const ex &subsed_ex = orig_ex.subs(m, options);
                        if (!are_ex_trivially_equal(orig_ex, subsed_ex)) {

				// Something changed: copy seq, subs, and return it
				std::unique_ptr<epvector> s(new epvector);
				s->reserve(seq.size());

				// Copy parts of seq which are known not to have changed
				s->insert(s->begin(), seq.begin(), cit);

				// Copy first changed element
				s->push_back(split_ex_to_pair(subsed_ex));
				++cit;

				// Copy rest
				while (cit != last) {
					s->push_back(split_ex_to_pair(recombine_pair_to_ex(*cit).subs(m, options)));
					++cit;
				}
				return s;
			}

			++cit;
		}

	} else {

		auto cit = seq.begin(), last = seq.end();
		while (cit != last) {
                        const ex &subsed_exr = cit->rest.subs(m, options);
                        const ex &subsed_exc = ex_to<numeric>(cit->coeff).subs(m, options);
                        if (not are_ex_trivially_equal(cit->rest, subsed_exr)
                            or not are_ex_trivially_equal(cit->coeff,
                                    subsed_exc)) {

                                // Something changed, copy seq, subs and return it
				std::unique_ptr<epvector> s(new epvector);
				s->reserve(seq.size());

				// Copy parts of seq which are known not to have changed
				s->insert(s->begin(), seq.begin(), cit);
			
				// Copy first changed element
                                if (is_exactly_a<numeric>(subsed_exc))
        				s->push_back(combine_ex_with_coeff_to_pair(subsed_exr, ex_to<numeric>(subsed_exc)));
                                else
        				s->push_back(split_ex_to_pair(mul(subsed_exr, subsed_exc)));
				++cit;

				// Copy rest
				while (cit != last) {
                                        const ex &sr = cit->rest.subs(m, options);
                                        const ex &sc = ex_to<numeric>(cit->coeff).subs(m, options);
					if (is_exactly_a<numeric>(sc))
                                                s->push_back(combine_ex_with_coeff_to_pair(sr, ex_to<numeric>(sc)));
                                        else
                                                s->push_back(split_ex_to_pair(mul(sr, sc)));
                                        ++cit;
				}
				return s;
			}

			++cit;
		}
	}
	
	// Nothing has changed
	return std::unique_ptr<epvector>(nullptr);
}


const epvector & expairseq::get_sorted_seq() const
{
	if (seq_sorted.empty())
		return seq;

        return seq_sorted;
}

bool expairseq::match(const ex & pattern, exmap& map) const
{
        CMatcher::level=0;
        CMatcher cm(*this, pattern, map);
        opt_exmap m = cm.get();
        if (not m)
                return false;
        map = m.value();
        return true;
}

//////////
// static member variables
//////////

#if EXPAIRSEQ_USE_HASHTAB
unsigned expairseq::maxhashtabsize = 0x4000000U;
unsigned expairseq::minhashtabsize = 0x1000U;
unsigned expairseq::hashtabfactor = 1;
#endif // EXPAIRSEQ_USE_HASHTAB

} // namespace GiNaC
