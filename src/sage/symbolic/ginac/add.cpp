/** @file add.cpp
 *
 *  Implementation of GiNaC's sums of expressions. */

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

#include "add.h"
#include "mul.h"
#include "archive.h"
#include "operators.h"
#include "utils.h"
#include "constant.h"
#include "infinity.h"
#include "compiler.h"
#include "order.h"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <string>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(add, expairseq,
  print_func<print_context>(&add::do_print).
  print_func<print_latex>(&add::do_print_latex).
  print_func<print_tree>(&add::do_print_tree).
  print_func<print_python_repr>(&add::do_print_python_repr))

//////////
// default constructor
//////////

add::add()
{
	tinfo_key = &add::tinfo_static;
}

//////////
// other constructors
//////////

// public

add::add(const ex & lh, const ex & rh)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = numeric(0);
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

add::add(const exvector & v, bool do_hold)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = numeric(0);
	construct_from_exvector(v, do_hold);
	GINAC_ASSERT(is_canonical());
}

add::add(const epvector & v)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = numeric(0);
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

add::add(const epvector & v, const numeric & oc)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = oc;
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

//add::add(std::unique_ptr<epvector> vp, const ex & oc)
//{
//	tinfo_key = &add::tinfo_static;
//	GINAC_ASSERT(vp.get()!=0);
//	overall_coeff = oc;
//	construct_from_epvector(*vp);
//	GINAC_ASSERT(is_canonical());
//}

//////////
// archiving
//////////

DEFAULT_ARCHIVING(add)

//////////
// functions overriding virtual functions from base classes
//////////

// public

void add::print_add(const print_context & c, unsigned level, bool latex) const
{
	if (precedence() <= level){
		if (latex)
			c.s << "{\\left(";
		else
			c.s << '(';
	}

	bool first = true;

	const epvector & sorted_seq = get_sorted_seq();
	// Then proceed with the remaining factors
	for (const auto & elem : sorted_seq) {
		std::stringstream tstream;
		std::unique_ptr<print_context> tcontext_p;
		if (latex) {
			tcontext_p.reset(new print_latex(tstream, c.options));
		} else {
			tcontext_p.reset(new print_dflt(tstream, c.options));
		}
		mul m = mul(elem.rest, elem.coeff);
		m.print(*tcontext_p, precedence());

		if (!first) {
			if (tstream.peek() == '-') {
				tstream.ignore();
				c.s << " - ";
			} else 
				c.s << " + ";
		} else {
			first = false;
		}
		tstream.get(*(c.s.rdbuf()));
	}

	// Finally print the "overall" numeric coefficient, if present.
	// This is just the constant coefficient. 
	if (not overall_coeff.is_zero()) {
		std::stringstream tstream;
		std::unique_ptr<print_context> tcontext_p;
		if (latex) {
			tcontext_p.reset(new print_latex(tstream, c.options));
		} else {
			tcontext_p.reset(new print_dflt(tstream, c.options));
		}
		overall_coeff.print(*tcontext_p, 0);
		if (!first) {
			if (tstream.peek() == '-') {
				c.s << " - ";
				tstream.ignore();
			} else
				c.s << " + ";
		}

		tstream.get(*(c.s.rdbuf()));
	}

	if (precedence() <= level) {
		if (latex)
			c.s << "\\right)}";
		else
			c.s << ')';
	}
}

void add::do_print(const print_context & c, unsigned level) const
{
	print_add(c, level, false);
}

void add::do_print_latex(const print_latex & c, unsigned level) const
{
	print_add(c, level, true);
}

void add::do_print_python_repr(const print_python_repr & c, unsigned /*level*/) const
{
	c.s << class_name() << '(';
	op(0).print(c);
	for (size_t i=1; i<nops(); ++i) {
		c.s << ',';
		op(i).print(c);
	}
	c.s << ')';
}

bool add::info(unsigned inf) const
{
        switch (inf) {
        case info_flags::nonzero:
                return is_positive()
                or info(info_flags::negative);
        case info_flags::polynomial:
        case info_flags::integer_polynomial:
        case info_flags::cinteger_polynomial:
        case info_flags::rational_polynomial:
        case info_flags::real:
        case info_flags::rational:
        case info_flags::integer:
        case info_flags::crational:
        case info_flags::cinteger:
        case info_flags::nonnegative:
        case info_flags::posint:
        case info_flags::nonnegint:
        case info_flags::even:
        case info_flags::crational_polynomial:
        case info_flags::rational_function: {
                for (const auto &elem : seq) {
                        if (!(recombine_pair_to_ex(elem).info(inf)))
                                return false;
                }
                if (overall_coeff.is_zero()
                    and (inf == info_flags::positive
                         or inf == info_flags::posint))
                        return true;
                return overall_coeff.info(inf);
        }
        case info_flags::inexact:
        case info_flags::algebraic: {
                if (overall_coeff.info(inf))
                        return true;
                for (const auto &elem : seq) {
                        if ((recombine_pair_to_ex(elem).info(inf)))
                                return true;
                }
                return false;
        }
        case info_flags::positive: {
                bool positive_seen = overall_coeff.is_positive();
                if (not positive_seen and not overall_coeff.is_zero())
                        return false;
                for (const auto &elem : seq) {
                        ex t = recombine_pair_to_ex(elem);
                        bool is_pos = t.is_positive();
                        if (not is_pos
                        and not t.info(info_flags::nonnegative))
                                return false;
                        if (not positive_seen and is_pos)
                                positive_seen = true;
                }
                return positive_seen;
        }
        case info_flags::negative: {
                bool negative_seen = overall_coeff.is_negative();
                if (not negative_seen and not overall_coeff.is_zero())
                        return false;
                for (const auto &elem : seq) {
                        ex t = recombine_pair_to_ex(elem);
                        bool is_neg = t.info(info_flags::negative);
                        if (not is_neg
                        and not t.is_zero())
                                return false;
                        if (not negative_seen and is_neg)
                                negative_seen = true;
                }
                return negative_seen;
        }
        }
        return inherited::info(inf);
}

bool add::is_polynomial(const ex & var) const
{
	for (const auto & elem : seq) {
		if (!(elem.rest).is_polynomial(var)) {
			return false;
		}
	}
	return true;
}

numeric add::degree(const ex & s) const
{
	numeric deg(seq[0].rest.degree(s));
	// Find maximum of degrees of individual terms
        for (const auto & elem : range(seq.begin()+1, seq.end())) {
		const numeric& t = elem.rest.degree(s);
                if (t > deg)
                        deg = t;
	}
        if (deg.is_negative() and not overall_coeff.is_zero())
                return *_num0_p;
	return deg;
}

numeric add::ldegree(const ex & s) const
{
	numeric deg(seq[0].rest.ldegree(s));
	// Find minimum of degrees of individual terms
        for (const auto & elem : range(seq.begin()+1, seq.end())) {
		const numeric& t = elem.rest.ldegree(s);
                if (t < deg)
                        deg = t;
	}
        if (deg.is_positive() and not overall_coeff.is_zero())
                return *_num0_p;
	return deg;
}

ex add::coeff(const ex & s, const ex & n) const
{
	epvector coeffseq;
        for (const auto & elem : seq) {
		ex restcoeff = elem.rest.coeff(s, n);
 		if (!restcoeff.is_zero()) {
			coeffseq.emplace_back(restcoeff, elem.coeff);
		}
	}

	return (new add(coeffseq,
	                n==0 ? overall_coeff : *_num0_p))->setflag(status_flags::dynallocated);
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x stands for a symbolic variables of type ex and c stands for such
 *  an expression that contain a plain number.
 *  - +(;c) -> c
 *  - +(x;0) -> x
 *
 *  @param level cut-off in recursive evaluation */
ex add::eval(int level) const
{
        if ((level == 1) and is_evaluated()) {
                GINAC_ASSERT(seq.size()>0);
                GINAC_ASSERT(seq.size()>1 || !overall_coeff.is_zero());
                return *this;
        }

	std::unique_ptr<epvector> evaled_seqp = evalchildren(level);
	if (unlikely(evaled_seqp != nullptr)) {
		// start over evaluating a new object
		return (new add(*evaled_seqp, overall_coeff))->
		       setflag(status_flags::dynallocated);
	}
	
#ifdef DO_GINAC_ASSERT
        for (const auto & elem : seq) {
		GINAC_ASSERT(!is_exactly_a<add>(elem.rest));
		if (is_exactly_a<numeric>(elem.rest))
			dbgprint();
		GINAC_ASSERT(!is_exactly_a<numeric>(elem.rest));
	}
#endif // def DO_GINAC_ASSERT
	
	// handle infinity
        for (auto i = seq.begin(); i != seq.end(); i++)
                if (unlikely(is_exactly_a<infinity>(i->rest)))
                        return eval_infinity(i);

	/** Perform automatic term rewriting rules */
	int seq_size = seq.size();
	if (seq_size == 0) {
		// +(;c) -> c
		return overall_coeff;
	}
        if (seq_size == 1 && overall_coeff.is_zero()) {
		// +(x;0) -> x
		return recombine_pair_to_ex(*(seq.begin()));
	} else if (!overall_coeff.is_zero() && seq[0].rest.return_type() != return_types::commutative) {
		throw (std::logic_error("add::eval(): sum of non-commutative objects has non-zero numeric term"));
	}

	// if any terms in the sum still are purely numeric, then they are more
	// appropriately collected into the overall coefficient
	int terms_to_collect = 0;
	for (const auto & elem : seq)
		if (unlikely(is_exactly_a<numeric>(elem.rest)))
			++terms_to_collect;
	if (terms_to_collect != 0) {
                epvector s;
		s.reserve(seq_size - terms_to_collect);
		numeric oc = *_num0_p;
                for (const auto & elem : seq)
			if (unlikely(is_exactly_a<numeric>(elem.rest)))
				oc = oc.add((ex_to<numeric>(elem.rest)).mul(ex_to<numeric>(elem.coeff)));
			else
				s.push_back(elem);
		return (new add(s, overall_coeff + oc))
		        ->setflag(status_flags::dynallocated);
	}

	return this->hold();
}



namespace { // anonymous namespace
	infinity infinity_from_iter(epvector::const_iterator i)
	{
		GINAC_ASSERT(is_exactly_a<infinity>(i->rest));
		GINAC_ASSERT(is_a<numeric>(i->coeff));
		infinity result = ex_to<infinity>(i->rest);
		result *= i->coeff;
		return result;
	}
} // end anonymous namespace


ex add::eval_infinity(epvector::const_iterator infinity_iter) const
{
	GINAC_ASSERT(is_exactly_a<infinity>(infinity_iter->rest));
	infinity result = infinity_from_iter(infinity_iter);

        for (auto i = seq.begin(); i != seq.end(); i++) {
                if (not is_exactly_a<infinity>(i->rest)) continue;
                if (i == infinity_iter) continue;
		infinity i_infty = infinity_from_iter(i);
		result += i_infty;
        }
	return result;
}


ex add::conjugate() const
{
	epvector v;
	v.reserve(seq.size());
	for (const auto & elem : seq)
		if ((elem.coeff).is_real()
		    and (elem.rest).is_real()) {
			v.push_back(elem);
		} else {
			ex cj=recombine_pair_to_ex(elem).conjugate();
                        v.push_back(split_ex_to_pair(cj));
		}
	// we know the conjugate is a numeric
        return (new add(v, overall_coeff.conj()))
		-> setflag(status_flags::dynallocated);
}

ex add::real_part() const
{
	epvector v;
	v.reserve(seq.size());
	for (const auto & elem : seq)
		if ((elem.coeff).is_real()) {
			ex rp = (elem.rest).real_part();
			if (!rp.is_zero())
				v.emplace_back(rp, elem.coeff);
		} else {
			ex rp=recombine_pair_to_ex(elem).real_part();
			if (!rp.is_zero())
				v.push_back(split_ex_to_pair(rp));
		}
	return (new add(v, ex_to<numeric>(overall_coeff.real_part())))
		-> setflag(status_flags::dynallocated);
}

ex add::imag_part() const
{
	epvector v;
	v.reserve(seq.size());
	for (const auto & elem : seq)
		if ((elem.coeff).is_real()) {
			ex ip = (elem.rest).imag_part();
			if (!ip.is_zero())
				v.emplace_back(ip, elem.coeff);
		} else {
			ex ip=recombine_pair_to_ex(elem).imag_part();
			if (!ip.is_zero())
				v.push_back(split_ex_to_pair(ip));
		}
	return (new add(v, ex_to<numeric>(overall_coeff.imag_part())))
		-> setflag(status_flags::dynallocated);
}

// protected

/** Implementation of ex::diff() for a sum. It differentiates each term.
 *  @see ex::diff */
ex add::derivative(const symbol & y) const
{
	epvector s;
	s.reserve(seq.size());
	
	// Only differentiate the "rest" parts of the expairs. This is faster
	// than the default implementation in basic::derivative() although
	// if performs the same function (differentiate each term).
	for (const auto & elem : seq)
		s.emplace_back(elem.rest.diff(y), elem.coeff);
	return (new add(s, *_num0_p))->setflag(status_flags::dynallocated);
}

int add::compare_same_type(const basic & other) const
{
	return inherited::compare_same_type(other);
}

unsigned add::return_type() const
{
	if (seq.empty())
		return return_types::commutative;
	
        return seq.begin()->rest.return_type();
}

tinfo_t add::return_type_tinfo() const
{
	if (seq.empty())
		return this;
	
        return seq.begin()->rest.return_type_tinfo();
}

// Note: do_index_renaming is ignored because it makes no sense for an add.
ex add::thisexpairseq(const epvector & v, const numeric & oc, bool /*do_index_renaming*/) const
{
	return (new add(v,oc))->setflag(status_flags::dynallocated);
}

//// Note: do_index_renaming is ignored because it makes no sense for an add.
ex add::thisexpairseq(std::unique_ptr<epvector> vp, const numeric & oc, bool /*do_index_renaming*/) const
{
	return (new add(*vp,oc))->setflag(status_flags::dynallocated);
}

expair add::split_ex_to_pair(const ex & e) const
{
	if (is_exactly_a<mul>(e)) {
		const mul &mulref(ex_to<mul>(e));
		const numeric &numfactor = mulref.overall_coeff;
		if (numfactor.is_one())
                        return expair(e, _ex1);
                auto mulcopyp = new mul(mulref);
		mulcopyp->overall_coeff = *_num1_p;
		mulcopyp->clearflag(status_flags::evaluated);
		mulcopyp->clearflag(status_flags::hash_calculated);
		mulcopyp->setflag(status_flags::dynallocated);
		return expair(*mulcopyp,numfactor);
	}
	return expair(e,_ex1);
}

expair add::combine_ex_with_coeff_to_pair(const ex & e,
                const numeric & c) const
{
	if (is_exactly_a<mul>(e)) {
		const mul &mulref(ex_to<mul>(e));
		const numeric &numfactor = mulref.overall_coeff;
		if (likely(numfactor.is_one()))
                        return expair(e, c);
		auto mulcopyp = new mul(mulref);
		mulcopyp->overall_coeff = *_num1_p;
		mulcopyp->clearflag(status_flags::evaluated);
		mulcopyp->clearflag(status_flags::hash_calculated);
		mulcopyp->setflag(status_flags::dynallocated);
		if (c.is_one())
			return expair(*mulcopyp, numfactor);
                return expair(*mulcopyp, numfactor*c);
	} else if (is_exactly_a<numeric>(e)) {
		if (c.is_one())
			return expair(e, _ex1);
		if (ex_to<numeric>(e).is_one())
			return expair(c, _ex1);
		return expair(ex_to<numeric>(e)*c, _ex1);
	}
	return expair(e, c);
}

expair add::combine_pair_with_coeff_to_pair(const expair & p,
                const numeric & c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));

	if (is_exactly_a<numeric>(p.rest)) {
		GINAC_ASSERT(ex_to<numeric>(p.coeff).is_one()); // should be normalized
                return expair(ex_to<numeric>(p.rest).mul_dyn(c),_ex1);
	}

	return expair(p.rest,ex_to<numeric>(p.coeff).mul_dyn(c));
}

ex add::recombine_pair_to_ex(const expair & p) const
{
	if (ex_to<numeric>(p.coeff).is_one())
		return p.rest;
	
        return (new mul(p.rest,p.coeff))->setflag(status_flags::dynallocated);
}

ex add::power(const numeric& expo) const
{
        using POW = class power;
        numeric icont = integer_content();
        const numeric lcoeff = ex_to<numeric>(lead_coeff()).div(icont);
        if (not lcoeff.is_positive())
                icont = icont.negative();
        if (icont.is_one()
            or not lcoeff.is_integer())
                return (new POW(*this, expo))->
                        setflag(status_flags::dynallocated
                                | status_flags::evaluated);

        ex c;
        if (expo.is_integer()) {
                c = icont.pow_intexp(expo.to_long());
        }
        else {
                bool ppower_equals_one;
                numeric newbasis, ppower;
                rational_power_parts(icont,
                                     expo,
                                     ppower,
                                     newbasis,
                                     ppower_equals_one);
                if (not newbasis.is_one()
                    and (ppower_equals_one or icont.is_negative()))
                        return (new POW(*this, expo))->
                                setflag(status_flags::dynallocated
                                      | status_flags::evaluated);
                if (ppower_equals_one)
                        c = *_num1_p;
                else
                        c = ppower * POW(newbasis, expo).hold();
        }

        auto addp = new add(*this);
        addp->setflag(status_flags::dynallocated);
        addp->clearflag(status_flags::hash_calculated);
        addp->overall_coeff /= icont;
        addp->seq_sorted.resize(0);
        for (auto &elem : addp->seq)
                elem.coeff = ex_to<numeric>(elem.coeff).div_dyn(icont);
        if (likely(not c.is_one()))
                return (new mul(POW(*addp, expo), c))->
                        setflag(status_flags::dynallocated
                                | status_flags::evaluated);
        return (new POW(*addp, expo))->
                setflag(status_flags::dynallocated
                        | status_flags::evaluated);
}

// This is an alternative version of add::power that is necessary to
// canonicalize powers of sums (always draws out the positive content).
// It would be difficult to introduce this as general behaviour in Sage.
ex Power(const ex& base_, const ex& expo_)
{
        if (not is_exactly_a<add>(base_)
            or not is_exactly_a<numeric>(expo_))
                return power(base_, expo_);
        using POW = class power;
        const add& base = ex_to<add>(base_);
        const numeric& expo = ex_to<numeric>(expo_);
        numeric icont = base.integer_content();
        const numeric lcoeff = ex_to<numeric>(base.lead_coeff()).div(icont);
        if (icont.is_one()
            or icont.is_minus_one()
            or not lcoeff.is_integer())
                return (new POW(base, expo))->
                        setflag(status_flags::dynallocated
                                | status_flags::evaluated);

        ex c;
        if (expo.is_integer()) {
                c = icont.pow_intexp(expo.to_long());
        }
        else {
                bool ppower_equals_one;
                numeric newbasis, ppower;
                rational_power_parts(icont,
                                     expo,
                                     ppower,
                                     newbasis,
                                     ppower_equals_one);
                if (ppower_equals_one) {
                        c = POW(icont.abs(), expo).hold();
                }
                else
                        if (icont.is_negative())
                                c = ppower * POW(-newbasis, expo).hold();
                        else
                                c = ppower * POW(newbasis, expo).hold();
        }

        auto addp = new add(base);
        addp->setflag(status_flags::dynallocated);
        addp->clearflag(status_flags::hash_calculated);
        addp->overall_coeff /= icont.abs();
        addp->seq_sorted.resize(0);
        for (auto &elem : addp->seq)
                elem.coeff = ex_to<numeric>(elem.coeff).div_dyn(icont.abs());
        if (likely(not c.is_one()))
                return (new mul(POW(*addp, expo), c))->
                        setflag(status_flags::dynallocated
                                | status_flags::evaluated);
        return (new POW(*addp, expo))->
                setflag(status_flags::dynallocated
                        | status_flags::evaluated);
}

ex add::expand(unsigned options) const
{
	std::unique_ptr<epvector> vp = expandchildren(options);
	if (vp == nullptr) {
		// the terms have not changed, so it is safe to declare this expanded
		return (options == 0) ? setflag(status_flags::expanded) : *this;
	}

	return (new add(*vp, overall_coeff))->setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
}

const epvector & add::get_sorted_seq() const
{
	if (seq_sorted.empty() and not seq.empty()) {
		seq_sorted = epvector(seq.size());
		partial_sort_copy(seq.begin(), seq.end(),
				  seq_sorted.begin(), seq_sorted.end(),
				  print_order_pair());
	}
	return expairseq::get_sorted_seq();
}

ex add::lead_coeff() const {
	// if sorted_seq was computed before, we don't have to search
	if (seq_sorted.empty() and not seq.empty()) {
		return min_element(seq.begin(), seq.end(),
				print_order_pair())->coeff;
	} 
		return seq_sorted.begin()->coeff;
	
}

ex add::combine_fractions() const
{
        epvector rseq;
        exmap fmap;
        ex oc = overall_coeff;
        for (const auto& pair : seq) {
                using POW = class power;
                const ex &e = recombine_pair_to_ex(pair);
                if (is_exactly_a<POW>(e)) {
                        const POW& p = ex_to<POW>(e);
                        const ex& eexp = p.op(1);
                        if (eexp.info(info_flags::negative)) {
                                auto it = fmap.find(p);
                                if (it != fmap.end())
                                        it->second += _ex1;
                                else
                                        fmap[p] = _ex1;
                        }
                        else
                                rseq.push_back(pair);
                }
                else if (is_exactly_a<mul>(e)) {
                        const mul& m = ex_to<mul>(e);
                        epvector denseq, coseq;
                        for (const auto& mpair : m.seq) {
                                if (mpair.coeff.info(info_flags::negative))
                                        denseq.push_back(mpair);
                                else
                                        coseq.push_back(mpair);
                        }
                        if (not denseq.empty()) {
                                mul den(denseq);
                                auto it = fmap.find(den);
                                mul tcoeff(coseq, ex_to<numeric>(pair.coeff));
                                if (it != fmap.end()) {
                                        it->second += tcoeff;
                                }
                                else
                                        fmap[den] = tcoeff;
                        }
                        else
                                rseq.push_back(pair);
                }
                else
                        rseq.push_back(pair);
        }
        for (const auto& t : fmap)
                rseq.push_back(split_ex_to_pair(mul(t.first, t.second)));
        add rex(rseq);
        return add(rex, oc);
}

} // namespace GiNaC
