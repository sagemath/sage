/** @file mul.cpp
 *
 *  Implementation of GiNaC's products of expressions. */

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

#include "mul.h"
#include "add.h"
#include "operators.h"
#include "lst.h"
#include "archive.h"
#include "utils.h"
#include "symbol.h"
#include "compiler.h"
#include "constant.h"
#include "infinity.h"
#include "function.h"
#include "inifcns.h"
#include "order.h"
#include "mpoly.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <sstream>
#ifdef DO_GINAC_ASSERT
#  include <typeinfo>
#endif

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(mul, expairseq,
  print_func<print_context>(&mul::do_print).
  print_func<print_latex>(&mul::do_print_latex).
  print_func<print_tree>(&mul::do_print_tree).
  print_func<print_python_repr>(&mul::do_print_python_repr))


//////////
// default constructor
//////////

mul::mul()
{
	tinfo_key = &mul::tinfo_static;
}

//////////
// other constructors
//////////

// public

mul::mul(const ex & lh, const ex & rh)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = *_num1_p;
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

mul::mul(const exvector & v, bool do_hold)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = *_num1_p;
	construct_from_exvector(v, do_hold);
	GINAC_ASSERT(is_canonical());
}

mul::mul(const epvector & v)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = *_num1_p;
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

mul::mul(const epvector & v, const numeric & oc, bool do_index_renaming)
{
	tinfo_key = &mul::tinfo_static;
	overall_coeff = oc;
	construct_from_epvector(v, do_index_renaming);
	GINAC_ASSERT(is_canonical());
}

//mul::mul(std::unique_ptr<epvector> vp, const ex & oc, bool do_index_renaming)
//{
//	tinfo_key = &mul::tinfo_static;
//	GINAC_ASSERT(vp.get()!=0);
//	overall_coeff = oc;
//	construct_from_epvector(*vp, do_index_renaming);
//	GINAC_ASSERT(is_canonical());
//}

mul::mul(const ex & lh, const ex & mh, const ex & rh)
{
	tinfo_key = &mul::tinfo_static;
	exvector factors;
	factors.reserve(3);
	factors.push_back(lh);
	factors.push_back(mh);
	factors.push_back(rh);
	overall_coeff = *_num1_p;
	construct_from_exvector(factors);
	GINAC_ASSERT(is_canonical());
}

//////////
// archiving
//////////

DEFAULT_ARCHIVING(mul)

//////////
// functions overriding virtual functions from base classes
//////////

void mul::print_overall_coeff(const ex& coeff_ex, const print_context & c,
		const char *mul_sym, bool latex) const
{
        if (not is_exactly_a<numeric>(coeff_ex))
                throw std::runtime_error("mul::print_overall_coeff: can't happen");
	const numeric & num_coeff = ex_to<numeric>(coeff_ex);
	std::stringstream tstream;
	std::unique_ptr<print_context> tcontext_p;
	if (latex)
		tcontext_p.reset(new print_latex(tstream, c.options));
	else
		tcontext_p.reset(new print_dflt(tstream, c.options));
	//print_context tcontext(tstream, c.options);
	num_coeff.print(*tcontext_p, 0);
	std::string coeffstr = tstream.str();

	bool parenthesis =((coeffstr.find(' ') != std::string::npos && !latex)||
		(coeffstr.find('+') != std::string::npos) ||
		(coeffstr.find('-',1) != std::string::npos));// ||
		//(coeffstr.find('/') != std::string::npos) ||
		//(coeffstr.find('*') != std::string::npos) ||
		//(coeffstr.find('^') != std::string::npos));
	if (num_coeff.is_minus_one())
		c.s<<"-";
	else if (parenthesis && coeffstr[0] == '-') {
		// We want to move the '-' out of the parenthesis if it is
		// the first character. This allows the printing function in
		// add objects to detect the '-' and change the sign.
		// We get x - (2*I - 1)*y instead of x + (-2*I + 1)*y.
		c.s<<"-";
		if (latex)
			c.s<<"\\left(";
		else
			c.s<<"(";
		tstream.str("");
		(-num_coeff).print(*tcontext_p, 0);
		c.s<<tstream.str();
		if (latex)
			c.s<<"\\right)";
		else
			c.s<<")";
		c.s << mul_sym;
	} else if (!num_coeff.is_integer() || !num_coeff.is_one()) {
		if (parenthesis) {
			if (latex)
				c.s << "\\left(";
			else
				c.s << '(';
		}
		c.s<<coeffstr;
		if (parenthesis) {
			if (latex)
				c.s << "\\right)";
			else
				c.s << ')';
		}
		c.s << mul_sym;
	}
}

/*
void mul::do_print(const print_context & c, unsigned level) const
{
	if (precedence() <= level)
		c.s << '(';

	print_overall_coeff(c, "*");

	epvector::const_iterator it = seq.begin(), itend = seq.end();
	bool first = true;
	while (it != itend) {
		if (!first)
			c.s << '*';
		else
			first = false;
		recombine_pair_to_ex(*it).print(c, precedence());
		++it;
	}

	if (precedence() <= level)
		c.s << ')';
}
*/

void mul::do_print(const print_context & c, unsigned level) const
{
	do_print_rat_func(c, level, false);
}

void mul::do_print_latex(const print_latex & c, unsigned level) const
{
	do_print_rat_func(c, level, true);
}

void mul::print_exvector(const exvector & v, const print_context & c,
		const char* sep) const
{
	bool first = true;
        for (const auto & elem : v) {
		if (!first)
			c.s << sep;
		else
			first = false;
		elem.print(c, precedence());
	}
}

void mul::do_print_rat_func(const print_context & c, unsigned level,
		bool latex_tags) const
{
	if (precedence() <= level){
		if (latex_tags)
			c.s << "\\left(";
		else
			c.s << '(';
	}

	const char *sep;
	if (latex_tags) {
		sep = (const char*)" ";
	} else {
		sep = (const char*)"*";
	}

	// Separate factors into those with negative numeric exponent
	// and all others
	const epvector & sorted_seq = get_sorted_seq();
        exvector neg_powers, others;
        for (const auto & elem : sorted_seq) {
		GINAC_ASSERT(is_exactly_a<numeric>(elem.coeff));
		if (ex_to<numeric>(elem.coeff).is_real() &&
				ex_to<numeric>(elem.coeff).is_negative())
			neg_powers.push_back(recombine_pair_to_ex(expair(elem.rest, -(elem.coeff))));
		else
			others.push_back(recombine_pair_to_ex(elem));
	}

	if (!neg_powers.empty()) {
		// Factors with negative exponent are printed as a fraction
		if (latex_tags) {
			const numeric& numer = overall_coeff.numer();
			bool negate = false;
			if (numer.is_minus_one()) {
				c.s<<"-";
				negate = true;
				//print_numer = coeff_numer.mul(*_num_1_p);
			} else {
				std::stringstream tstream;
				std::unique_ptr<print_context> tcontext_p(new print_latex(tstream, c.options));
				numer.print(*tcontext_p, 0);
				if (tstream.peek() == '-') {
					c.s<<"-";
					negate = true;
					//const numeric &print_numer = coeff_numer.mul(*_num_1_p);
				} /*else {
					const numeric &print_numer = coeff_numer;
				}*/
			}

			numeric print_numer = negate ? -numer : numer;
			c.s << "\\frac{";
			if (others.empty()) {
				if (print_numer.is_integer() && print_numer.is_one()){
					c.s<<'1';
				} else {
					print_numer.print(c);
				}
			} else {
				if (print_numer.is_integer() && print_numer.is_one()){
					mul(others).eval().print(c);
				} else {
					mul(print_numer,mul(others).eval()).hold().print(c);
				}
			}
			c.s << "}{";
			numeric denom = overall_coeff.denom();
			if (denom.is_one()) {
				mul(neg_powers).eval().print(c);
			} else {
				mul(denom, mul(neg_powers).eval()).hold().print(c);
			}
			c.s << "}";
		} else {
			print_overall_coeff(overall_coeff, c,
					others.empty() ? "" : sep,
					latex_tags);
			if (others.empty()
                            and (overall_coeff.is_one()
                                 or overall_coeff.is_minus_one())) {
				c.s<<"1";
			} else {
				print_exvector(others, c, sep);
			}
			c.s << "/";

			if (neg_powers.size() > 1) {
				c.s<<"(";
			}
			print_exvector(neg_powers, c, sep);
			if (neg_powers.size() > 1) {
				c.s << ")";
			}
		}

	} else {
		// negative powers are empty
		if (!latex_tags) {
			print_overall_coeff(overall_coeff, c, sep, latex_tags);

			print_exvector(others, c, sep);
		} else {
			// Decide which separator to use after the overall
			// coefficient.
			// overall coefficient is printed first in the
			// output, if the first character after the
			// coefficient is a number, then it is hard to
			// distinguish these in the latex output:
			//
			// sage: e = 2 * 2^(1/3)
			// sage: print latex(e)
			// 2 \, 2^{\left(\frac{1}{3}\right)}
			//
			// In this case instead of the usual separator \,
			// we use \cdot

			std::stringstream tstream;
			print_latex tcontext(tstream, c.options);
			print_exvector(others, tcontext, sep);
			print_overall_coeff(overall_coeff, c,
					std::isdigit(tstream.peek()) != 0 ?
						" \\cdot " : " \\, ",
					latex_tags);
			c.s<<tstream.str();
		}
	}
	if (precedence() <= level){
		if (latex_tags)
			c.s << "\\right)";
		else
			c.s << ')';
	}

}

void mul::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(';
	op(0).print(c);
	for (size_t i=1; i<nops(); ++i) {
		c.s << ',';
		op(i).print(c);
	}
	c.s << ')';
}

bool mul::info(unsigned inf) const
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
        case info_flags::crational_polynomial:
        case info_flags::rational_function: {
                if (not overall_coeff.is_real())
                        return false;
                for (const auto &elem : seq)
                        if (!(recombine_pair_to_ex(elem).info(inf)))
                                return false;
                if (overall_coeff.is_one() && inf == info_flags::even)
                        return true;
                return overall_coeff.info(inf);
        }
        case info_flags::even: {
                if (not overall_coeff.is_integer())
                        return false;
                bool even_seen = false;
                for (const auto &elem : seq) {
                        const ex &e = recombine_pair_to_ex(elem);
                        if (not e.is_integer())
                                return false;
                        if (e.info(info_flags::even))
                                even_seen = true;
                }
                return even_seen or overall_coeff.is_even();
        }
        case info_flags::inexact:
        case info_flags::algebraic: {
                if (overall_coeff.info(inf))
                        return true;
                for (const auto &elem : seq)
                        if ((recombine_pair_to_ex(elem).info(inf)))
                                return true;
                return false;
        }
        case info_flags::positive:
        case info_flags::negative: {
                if (not overall_coeff.is_real())
                        return false;
                bool pos = true;
                for (const auto &elem : seq) {
                        const ex &factor = recombine_pair_to_ex(elem);
                        if (factor.is_positive())
                                continue;
                        if (factor.info(info_flags::negative))
                                pos = !pos;
                        else
                                return false;
                }
                if (overall_coeff.info(info_flags::negative))
                        pos = !pos;
                return (inf == info_flags::positive ? pos : !pos);
        }
        case info_flags::nonnegative: {
                if (not overall_coeff.is_real())
                        return false;
                bool pos = true;
                for (const auto &elem : seq) {
                        const ex &factor = recombine_pair_to_ex(elem);
                        if (factor.info(info_flags::nonnegative)
                            or factor.is_positive())
                                continue;
                        if (factor.info(info_flags::negative))
                                pos = !pos;
                        else
                                return false;
                }
                return (overall_coeff.info(info_flags::negative) ? !pos : pos);
        }
        case info_flags::posint:
        case info_flags::negint: {
                if (not overall_coeff.is_real())
                        return false;
                bool pos = true;
                for (const auto &elem : seq) {
                        const ex &factor = recombine_pair_to_ex(elem);
                        if (factor.info(info_flags::posint))
                                continue;
                        if (factor.info(info_flags::negint))
                                pos = !pos;
                        else
                                return false;
                }
                if (overall_coeff.info(info_flags::negint))
                        pos = !pos;
                else if (!overall_coeff.info(info_flags::posint))
                        return false;
                return (inf == info_flags::posint ? pos : !pos);
        }
        case info_flags::nonnegint: {
                if (not overall_coeff.is_real())
                        return false;
                bool pos = true;
                for (const auto &elem : seq) {
                        const ex &factor = recombine_pair_to_ex(elem);
                        if (factor.info(info_flags::nonnegint)
                            or factor.info(info_flags::posint))
                                continue;
                        if (factor.info(info_flags::negint))
                                pos = !pos;
                        else
                                return false;
                }
                if (overall_coeff.info(info_flags::negint))
                        pos = !pos;
                else if (!overall_coeff.info(info_flags::posint))
                        return false;
                return pos;
        }
        }
        return inherited::info(inf);
}

bool mul::is_polynomial(const ex & var) const
{
	for (const auto & elem : seq) {
		if (!elem.rest.is_polynomial(var) ||
		    (elem.rest.has(var) && !elem.coeff.info(info_flags::nonnegint))) {
			return false;
		}
	}
	return true;
}

numeric mul::degree(const ex & s) const
{
	// Sum up degrees of factors
	numeric deg_sum(0);
        for (const auto & elem : seq) {
		if (is_exactly_a<numeric>(elem.coeff)) {
                        const numeric& n = ex_to<numeric>(elem.coeff);
                        if (n.is_real()) {
			        deg_sum += n * elem.rest.degree(s);
                        }
                        else {
                                if (elem.rest.has(s))
                                        throw std::runtime_error("mul::degree() undefined degree because of complex exponent");
		        }
                }
                else
		        throw std::runtime_error("mul::degree() undefined degree because of non-integer exponent");
        }
	return deg_sum;
}

numeric mul::ldegree(const ex & s) const
{
	// Sum up degrees of factors
	numeric deg_sum(0);
        for (const auto & elem : seq) {
		if (is_exactly_a<numeric>(elem.coeff)) {
                        const numeric& n = ex_to<numeric>(elem.coeff);
                        if (n.is_real()) {
			        deg_sum += n * elem.rest.ldegree(s);
                        }
                        else {
                                if (elem.rest.has(s))
                                        throw std::runtime_error("mul::ldegree() undefined degree because of complex exponent");
		        }
                }
                else
		        throw std::runtime_error("mul::ldegree() undefined degree because of non-integer exponent");
        }
	return deg_sum;
}

ex mul::coeff(const ex & s, const ex & n) const
{
	exvector coeffseq;
	coeffseq.reserve(seq.size()+1);
	
	if (n.is_zero()) {
		// product of individual coeffs
		// if a non-zero power of s is found, the resulting product will be 0
                for (const auto & elem : seq)
			coeffseq.push_back(recombine_pair_to_ex(elem).coeff(s,n));
		coeffseq.emplace_back(overall_coeff);
		return (new mul(coeffseq))->setflag(status_flags::dynallocated);
	}
	
	bool coeff_found = false;
        for (const auto & elem : seq) {
		const ex& t = recombine_pair_to_ex(elem);
		const ex& c = t.coeff(s, n);
		if (!c.is_zero()) {
			coeffseq.push_back(c);
			coeff_found = true;
		} else {
			coeffseq.push_back(t);
		}
	}
	if (coeff_found) {
		coeffseq.emplace_back(overall_coeff);
		return (new mul(coeffseq))->setflag(status_flags::dynallocated);
	}
	
	return _ex0;
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x, x1, x2,... stand for a symbolic variables of type ex and c, c1, c2...
 *  stand for such expressions that contain a plain number.
 *  - *(...,x;0) -> 0
 *  - *(+(x1,x2,...);c) -> *(+(*(x1,c),*(x2,c),...))
 *  - *(x;1) -> x
 *  - *(;c) -> c
 *
 *  @param level cut-off in recursive evaluation */
ex mul::eval(int level) const
{
#ifdef DO_GINAC_ASSERT
        for (const auto & elem : seq) {
		GINAC_ASSERT((!is_exactly_a<mul>(elem.rest)) ||
		             (!(ex_to<numeric>(elem.coeff).is_integer())));
		GINAC_ASSERT(!(elem.is_canonical_numeric()));
		if (is_exactly_a<numeric>(recombine_pair_to_ex(elem)))
		    print(print_tree(std::cerr));
		GINAC_ASSERT(!is_exactly_a<numeric>(recombine_pair_to_ex(elem)));
		/* for paranoia */
		//   The following test will fail on sage: exp(x)*exp(x)
		//   That's probably not an issue, but should be investigated.
		// expair p = split_ex_to_pair(recombine_pair_to_ex(*i));
		// GINAC_ASSERT(p.rest.is_equal(i->rest));
		// GINAC_ASSERT(p.coeff.is_equal(i->coeff));
		/* end paranoia */
	}
#endif // def DO_GINAC_ASSERT
	
	if (level == 1 and is_evaluated()) {
		GINAC_ASSERT(seq.size()>0);
		GINAC_ASSERT(seq.size()>1 || !overall_coeff.is_one());
		return *this;
	}

	std::unique_ptr<epvector> evaled_seqp = evalchildren(level);
	if (unlikely(evaled_seqp != nullptr)) {
		// do more evaluation later
		return (new mul(*evaled_seqp, overall_coeff))->
		           setflag(status_flags::dynallocated);
	}
	
	// handle infinity, simple way
        if (std::any_of(seq.cbegin(), seq.cend(),
                        [](expair p)
                                { return is_exactly_a<infinity>(p.rest); })) {
                infinity result = infinity::from_sign(1);
                result *= overall_coeff;
                for (const auto& p : seq)
		        result *= recombine_pair_to_ex(p);
                return result;
        }

        // handle exp(a)*exp(b) -> exp(a+b) and
	unsigned exp_count = 0;
	for (auto i = seq.begin(); i != seq.end(); i++) {
		const numeric& num_coeff = ex_to<numeric>(i->coeff);
		if (unlikely(is_ex_the_function(i->rest, exp) and
			     num_coeff.is_integer())) {
			exp_count++;
			if (exp_count>1 or not num_coeff.is_one() or
                                not num_coeff.is_integer())
				return eval_exponentials();
		}
	}

	// perform the remaining automatic rewrites
	size_t seq_size = seq.size();
	if (overall_coeff.is_zero()) {
		// *(...,x;0) -> 0
		return overall_coeff;
	}
        if (seq_size == 0) {
		// *(;c) -> c
		return overall_coeff;
	}
        if (seq_size == 1
            and overall_coeff.is_one()) {
		// *(x;1) -> x
		// except in positive characteristic: 1*(x+2) = x in F_2
		return recombine_pair_to_ex(*(seq.begin()));
	}
        if (seq_size == 1
	     and is_exactly_a<add>((*seq.begin()).rest)
	     and ex_to<numeric>((*seq.begin()).coeff).is_one()) {
		// *(+(x,y,...);c) -> +(*(x,c),*(y,c),...) (c numeric(), no powers of +())
		const add & addref = ex_to<add>((*seq.begin()).rest);
		epvector distrseq;
		distrseq.reserve(addref.seq.size());
                for (const auto & elem : addref.seq)
			distrseq.push_back(addref.combine_pair_with_coeff_to_pair(elem, overall_coeff));
		const ex& x = (new add(distrseq,
		                addref.overall_coeff.mul(overall_coeff))
		       )->setflag(status_flags::dynallocated | status_flags::evaluated);
                const add & result = ex_to<add>(x);
                if (result.seq.empty())
                        return result.overall_coeff;
                return x;
	}
        if (seq_size >= 2
            and ((flags & status_flags::expanded) == 0u)) {
		// Strip the content and the unit part from each term. Thus
		// things like (-x+a)*(3*x-3*a) automagically turn into - 3*(x-a)^2
		auto last = seq.end();
		auto i = seq.begin();
		auto j = seq.begin();
		epvector s;
		numeric oc = *_num1_p;
		bool something_changed = false;
		while (i!=last) {
			if (likely(! (is_exactly_a<add>(i->rest) && i->coeff.is_one()))) {
				// power::eval has such a rule, no need to handle powers here
				++i;
				continue;
			}

			// XXX: What is the best way to check if the polynomial is a primitive? 
			numeric c = i->rest.integer_content();
#ifdef DO_GINAC_ASSERT
                        if (c.is_zero())
                                std::cerr << i->rest << "\n" << std::flush;
#endif
			const numeric lead_coeff =
				ex_to<numeric>(ex_to<add>(i->rest).\
						lead_coeff()).div(c);
			const bool canonicalizable = (
                                        lead_coeff.is_exact()
                                        and lead_coeff.is_integer());

			// The following comment is no longer true for pynac.
			// We use the print order to determine the main variable
			// This order is not random.
			// XXX: The main variable is chosen in a random way, so this code 
			// does NOT transform the term into the canonical form (thus, in some
			// very unlucky event it can even loop forever). Hopefully the main
			// variable will be the same for all terms in *this
			const bool unit_normal = lead_coeff.is_pos_integer();
			if (likely(((not canonicalizable) or unit_normal)
                                   and c.is_inexact_one() )) {
				++i;
				continue;
			}

			if (! something_changed) {
				s.reserve(seq_size);
				something_changed = true;
			}

			while ((j!=i) && (j!=last)) {
				s.push_back(*j);
				++j;
			}

                        if (!unit_normal)
                                c = c.mul(*_num_1_p);

                        oc = oc.mul(c);

			// divide add by the number in place to save at least 2 .eval() calls
			const add& addref = ex_to<add>(i->rest);
			auto primitive = new add(addref);
			primitive->setflag(status_flags::dynallocated);
			primitive->clearflag(status_flags::hash_calculated);
			primitive->overall_coeff /= c;
			primitive->seq_sorted.resize(0);
			for (auto & elem : primitive->seq)
				elem.coeff = ex_to<numeric>(elem.coeff).div_dyn(c);

			s.emplace_back(*primitive, _ex1);

			++i;
			++j;
		}
		if (something_changed) {
			while (j!=last) {
				s.push_back(*j);
				++j;
			}
			if (s.empty()) {
				return overall_coeff.mul(oc);
			}
			return (new mul(s, overall_coeff.mul(oc))
			       )->setflag(status_flags::dynallocated);
		}
	}
	
	return this->hold();
}
	


ex mul::eval_exponentials() const
{
	ex exp_arg = _ex0;
	numeric oc = *_num1_p;
	epvector s;
	s.reserve(seq.size());

	for (const auto & elem : seq) {
		const numeric & num_coeff = ex_to<numeric>(elem.coeff);
		const bool simplifyable_exp = is_ex_the_function(elem.rest, exp) and num_coeff.is_integer();
		if (likely(not simplifyable_exp))
			s.push_back(elem);
		else
			exp_arg += elem.rest.op(0) * num_coeff;
	}

	ex new_exp = exp(exp_arg);
	if (is_exactly_a<numeric>(new_exp))
		oc = oc.mul(ex_to<numeric>(new_exp));
	else
		s.emplace_back(new_exp, _ex1);

	mul * result = new mul(s, overall_coeff.mul(oc));
	return result->setflag(status_flags::dynallocated);
}

ex mul::evalf(int level, PyObject* parent) const
{
	if (level==1)
		return mul(seq,overall_coeff);
	
	if (level==-max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	epvector s;
	s.reserve(seq.size());

	--level;
        for (const auto & elem : seq)
               s.push_back(combine_ex_with_coeff_to_pair(
                        elem.rest.evalf(level, parent),
                        ex_to<numeric>(elem.coeff)));
        return mul(s, ex_to<numeric>(overall_coeff.evalf(level, parent)));
}

void mul::find_real_imag(ex & rp, ex & ip) const
{
	rp = overall_coeff.real_part();
	ip = overall_coeff.imag_part();
	for (const auto & elem : seq) {
		ex factor = recombine_pair_to_ex(elem);
		ex new_rp = factor.real_part();
		ex new_ip = factor.imag_part();
		if(new_ip.is_zero()) {
			rp *= new_rp;
			ip *= new_rp;
		} else {
			ex temp = rp*new_rp - ip*new_ip;
			ip = ip*new_rp + rp*new_ip;
			rp = temp;
		}
	}
	rp = rp.expand();
	ip = ip.expand();
}

ex mul::real_part() const
{
	ex rp, ip;
	find_real_imag(rp, ip);
	return rp;
}

ex mul::imag_part() const
{
	ex rp, ip;
	find_real_imag(rp, ip);
	return ip;
}

bool tryfactsubs(const ex & origfactor, const ex & patternfactor, int & nummatches, lst & repls)
{	
	ex origbase;
	int origexponent;
	int origexpsign;

	if (is_exactly_a<power>(origfactor)
            and origfactor.op(1).is_integer()) {
		origbase = origfactor.op(0);
		int expon = ex_to<numeric>(origfactor.op(1)).to_int();
		origexponent = expon > 0 ? expon : -expon;
		origexpsign = expon > 0 ? 1 : -1;
	} else {
		origbase = origfactor;
		origexponent = 1;
		origexpsign = 1;
	}

	ex patternbase;
	int patternexponent;
	int patternexpsign;

	if (is_exactly_a<power>(patternfactor)
            and patternfactor.op(1).is_integer()) {
		patternbase = patternfactor.op(0);
		int expon = ex_to<numeric>(patternfactor.op(1)).to_int();
		patternexponent = expon > 0 ? expon : -expon;
		patternexpsign = expon > 0 ? 1 : -1;
	} else {
		patternbase = patternfactor;
		patternexponent = 1;
		patternexpsign = 1;
	}

	lst saverepls = repls;
	if (origexponent < patternexponent || origexpsign != patternexpsign || !origbase.match(patternbase,saverepls))
		return false;
	repls = saverepls;

	int newnummatches = origexponent / patternexponent;
	if (newnummatches < nummatches)
		nummatches = newnummatches;
	return true;
}

/** Checks whether e matches to the pattern pat and the (possibly to be updated)
  * list of replacements repls. This matching is in the sense of algebraic
  * substitutions. Matching starts with pat.op(factor) of the pattern because
  * the factors before this one have already been matched. The (possibly
  * updated) number of matches is in nummatches. subsed[i] is true for factors
  * that already have been replaced by previous substitutions and matched[i]
  * is true for factors that have been matched by the current match.
  */
bool algebraic_match_mul_with_mul(const mul &e, const ex &pat, lst &repls,
		unsigned factor, int &nummatches, const std::vector<bool> &subsed,
		std::vector<bool> &matched)
{
	GINAC_ASSERT(subsed.size() == e.nops());
	GINAC_ASSERT(matched.size() == e.nops());

	if (factor == pat.nops())
		return true;

	for (size_t i=0; i<e.nops(); ++i) {
		if(subsed[i] || matched[i])
			continue;
		lst newrepls = repls;
		int newnummatches = nummatches;
		if (tryfactsubs(e.op(i), pat.op(factor),
					newnummatches, newrepls)) {
			matched[i] = true;
			if (algebraic_match_mul_with_mul(e, pat, newrepls, factor+1,
					newnummatches, subsed, matched)) {
				repls = newrepls;
				nummatches = newnummatches;
				return true;
			}
			matched[i] = false;
		}
	}

	return false;
}

bool mul::has(const ex & pattern, unsigned options) const
{
	if((options&has_options::algebraic) == 0u)
		return basic::has(pattern,options);
	if(is_exactly_a<mul>(pattern)) {
		lst repls;
		int nummatches = std::numeric_limits<int>::max();
		std::vector<bool> subsed(nops(), false);
		std::vector<bool> matched(nops(), false);
		if(algebraic_match_mul_with_mul(*this, pattern, repls, 0, nummatches,
				subsed, matched))
			return true;
	}
	return basic::has(pattern, options);
}

ex mul::algebraic_subs_mul(const exmap & m, unsigned options) const
{	
	std::vector<bool> subsed(nops(), false);
	ex divide_by = 1;
	ex multiply_by = 1;

	for (const auto & elem : m) {

		if (is_exactly_a<mul>(elem.first)) {
retry1:
			int nummatches = std::numeric_limits<int>::max();
			std::vector<bool> currsubsed(nops(), false);
			lst repls;
			
			if(!algebraic_match_mul_with_mul(*this, elem.first, repls, 0, nummatches, subsed, currsubsed))
				continue;

			for (size_t j=0; j<subsed.size(); j++)
				if (currsubsed[j])
					subsed[j] = true;
			ex subsed_pattern
				= elem.first.subs(ex(repls), subs_options::no_pattern);
			divide_by *= power(subsed_pattern, nummatches);
			ex subsed_result
				= elem.second.subs(ex(repls), subs_options::no_pattern);
			multiply_by *= power(subsed_result, nummatches);
			goto retry1;

		} else {

			for (size_t j=0; j<this->nops(); j++) {
				int nummatches = std::numeric_limits<int>::max();
				lst repls;
				if (!subsed[j] && tryfactsubs(op(j), elem.first, nummatches, repls)){
					subsed[j] = true;
					ex subsed_pattern
						= elem.first.subs(ex(repls), subs_options::no_pattern);
					divide_by *= power(subsed_pattern, nummatches);
					ex subsed_result
						= elem.second.subs(ex(repls), subs_options::no_pattern);
					multiply_by *= power(subsed_result, nummatches);
				}
			}
		}
	}

        if (std::any_of(subsed.cbegin(), subsed.cend(),
                                [](bool b) { return b; } ))
		return subs_one_level(m, options | subs_options::algebraic);
	return ((*this)/divide_by)*multiply_by;
}

ex mul::conjugate() const
{
	// The base class' method is wrong here because we have to be careful at
	// branch cuts. power::conjugate takes care of that already, so use it.
	std::unique_ptr<epvector> newepv(nullptr);
        for (auto i=seq.begin(); i!=seq.end(); ++i) {
		if (newepv) {
			newepv->push_back(split_ex_to_pair(recombine_pair_to_ex(*i).conjugate()));
			continue;
		}
		ex x = recombine_pair_to_ex(*i);
		ex c = x.conjugate();
		if (c.is_equal(x)) {
			continue;
		}
		newepv.reset(new epvector);
		newepv->reserve(seq.size());
		for (auto j=seq.begin(); j!=i; ++j) {
			newepv->push_back(*j);
		}
		newepv->push_back(split_ex_to_pair(c));
	}
        // known to be numeric
	const numeric& x = overall_coeff.conj();
	if (not newepv and x.is_equal(overall_coeff)) {
		return *this;
	}
	ex result = thisexpairseq(newepv ? *newepv : seq, x);
	return result;
}


// protected

/** Implementation of ex::diff() for a product.  It applies the product rule.
 *  @see ex::diff */
ex mul::derivative(const symbol & s) const
{
	exvector addseq;
	addseq.reserve(seq.size());

	// D(a*b*c) = D(a)*b*c + a*D(b)*c + a*b*D(c)
	epvector mulseq = seq;
	auto i = seq.begin(), end = seq.end();
	auto i2 = mulseq.begin();
	while (i != end) {
		expair ep = split_ex_to_pair(power(i->rest, i->coeff - _ex1) *
		                             i->rest.diff(s));
		ep.swap(*i2);
		addseq.emplace_back((new mul(mulseq,
                     overall_coeff.mul(ex_to<numeric>(i->coeff))))->setflag(status_flags::dynallocated));
		ep.swap(*i2);
		++i; ++i2;
	}
	return (new add(addseq))->setflag(status_flags::dynallocated);
}

/** Total degree of a mul object.
 Beware that symbolic exponents are wrapped inside pow objects with 1 as coeff. */
double mul::total_degree() const
{
	if ((flags & status_flags::tdegree_calculated) != 0u) {
		return tdegree;
	}
	numeric tdeg = calc_total_degree();
	if (tdeg.is_real())
		tdegree = tdeg.to_double();
	else
		tdegree = std::sqrt(std::pow(tdeg.real().to_double(), 2) + 
				std::pow(tdeg.imag().to_double(), 2));
	setflag(status_flags::tdegree_calculated);
	return tdegree;
}

int mul::compare_same_type(const basic & other) const
{
	return inherited::compare_same_type(other);
}

unsigned mul::return_type() const
{
	if (seq.empty()) {
		// mul without factors: should not happen, but commutates
		return return_types::commutative;
	}
	
	bool all_commutative = true;
        expair noncommutative_element;
	
        for (const auto & elem : seq) {
		unsigned rt = elem.rest.return_type();
		if (rt == return_types::noncommutative_composite)
			return rt; // one ncc -> mul also ncc
		if (rt == return_types::noncommutative) {
                        if (all_commutative) {
                                // first nc element found, remember position
                                noncommutative_element = elem;
                                all_commutative = false;
			}
                        else {
                                // another nc element found, compare type_infos
                                if (noncommutative_element.rest.return_type_tinfo() != elem.rest.return_type_tinfo()) {
                                                // different types -> mul is ncc
                                        return return_types::noncommutative_composite;
                                }
        		}
                }
	}
	// all factors checked
	return all_commutative ? return_types::commutative : return_types::noncommutative;
}
   
tinfo_t mul::return_type_tinfo() const
{
	if (seq.empty())
		return this;  // mul without factors: should not happen
	
	// return type_info of first noncommutative element
        for (const auto & elem : seq)
		if (elem.rest.return_type() == return_types::noncommutative)
			return elem.rest.return_type_tinfo();

        // no noncommutative element found, should not happen
	return this;
}

ex mul::thisexpairseq(const epvector & v, const numeric & oc, bool do_index_renaming) const
{
	return (new mul(v, oc, do_index_renaming))->setflag(status_flags::dynallocated);
}

ex mul::thisexpairseq(std::unique_ptr<epvector> vp, const numeric & oc, bool do_index_renaming) const
{
	return (new mul(*vp, oc, do_index_renaming))->setflag(status_flags::dynallocated);
}

expair mul::split_ex_to_pair(const ex & e) const
{
	if (is_exactly_a<power>(e)) {
		const power & powerref = ex_to<power>(e);
		if (is_exactly_a<numeric>(powerref.exponent))
			return expair(powerref.basis,powerref.exponent);
	}
	return expair(e,_ex1);
}
	
expair mul::combine_ex_with_coeff_to_pair(const ex & e,
                                          const numeric & c) const
{
        // First, try a common shortcut:
        if (is_exactly_a<symbol>(e))
                return expair(e, c);

        // to avoid duplication of power simplification rules,
	// we create a temporary power object
	// otherwise it would be hard to correctly evaluate
	// expression like (4^(1/3))^(3/2)
	if (c.is_one())
		return split_ex_to_pair(e);

	return split_ex_to_pair(power(e,c));
}
	
expair mul::combine_pair_with_coeff_to_pair(const expair & p,
                                            const numeric & c) const
{
        GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));

        if (is_exactly_a<symbol>(p.rest))
                return expair(p.rest, p.coeff*c);
	if (c.is_one())
		return p;
        if (p.coeff.is_one())
                return expair(p.rest, c);

	// to avoid duplication of power simplification rules,
	// we create a temporary power object
	// otherwise it would be hard to correctly evaluate
	// expression like (4^(1/3))^(3/2)
	return split_ex_to_pair(power(recombine_pair_to_ex(p),c));
}
	
ex mul::recombine_pair_to_ex(const expair & p) const
{
        if (unlikely(is_exactly_a<infinity>(p.rest))) {
                infinity res = ex_to<infinity>(p.rest);
                res *= p.coeff;
                return res;
        }
	if (p.coeff.is_one()) 
		return p.rest;
	return (new power(p.rest,p.coeff))->setflag(status_flags::dynallocated);
}

bool mul::expair_needs_further_processing(epp it)
{
	if (is_exactly_a<mul>(it->rest) &&
		ex_to<numeric>(it->coeff).is_integer()) {
		// combined pair is product with integer power -> expand it
		*it = split_ex_to_pair(recombine_pair_to_ex(*it));
		return true;
	}
	if (is_exactly_a<numeric>(it->rest)) {
		expair ep = split_ex_to_pair(recombine_pair_to_ex(*it));
		if (!ep.is_equal(*it)) {
			// combined pair is a numeric power which can be simplified
			*it = ep;
			return true;
		}
		if (it->coeff.is_one()) {
			// combined pair has coeff 1 and must be moved to the end
			return true;
		}
	}
	return false;
}       

numeric mul::default_overall_coeff() const
{
	return *_num1_p;
}

void mul::combine_overall_coeff(const numeric & c)
{
	overall_coeff *= c;
}

void mul::combine_overall_coeff(const numeric & c1, const numeric & c2)
{
	ex e = c1.power(c2);
	if (not is_exactly_a<numeric>(e))
                throw std::runtime_error("mul::combine_overall_coeff: can't happen");
	overall_coeff *= ex_to<numeric>(e);
}

bool mul::can_make_flat(const expair & p) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));

        // (x*y)^c == x^c*y^c  if c ∈ ℤ
        return p.coeff.is_integer();
}

bool mul::can_be_further_expanded(const ex & e)
{
	if (is_exactly_a<mul>(e)) {
		for (const auto & elem : ex_to<mul>(e).seq) {
			if (is_exactly_a<add>(elem.rest) && elem.coeff.info(info_flags::posint))
				return true;
		}
	} else if (is_exactly_a<power>(e)) {
		if (is_exactly_a<add>(e.op(0)) && e.op(1).info(info_flags::posint))
			return true;
	}
	return false;
}

ex mul::expand(unsigned options) const
{
	// trivial case: expanding the monomial (~ 30% of all calls)
        bool all_intsym = true;
        for (const auto & elem : seq) 
                if (not is_exactly_a<symbol>(elem.rest)
                    or not elem.coeff.is_integer())
                {
                        all_intsym = false;
                        break;
                }
        if (all_intsym) {
                setflag(status_flags::expanded);
                return *this;
        }

	// First, expand the children
	std::unique_ptr<epvector> expanded_seqp = expandchildren(options);
	const epvector & expanded_seq = (expanded_seqp != nullptr ? *expanded_seqp : seq);

	// Now, look for all the factors that are sums and multiply each one out
	// with the next one that is found while collecting the factors which are
	// not sums
	ex last_expanded = _ex1;

	epvector non_adds;
	non_adds.reserve(expanded_seq.size());

	for (const auto & elem : expanded_seq) {
		if (is_exactly_a<add>(elem.rest) &&
			(elem.coeff.is_one())) {
			if (is_exactly_a<add>(last_expanded)) {
				// Expand a product of two sums, aggressive version.
				// Caring for the overall coefficients in separate loops can
				// sometimes give a performance gain of up to 15%!

				const int sizedifference = ex_to<add>(last_expanded).seq.size()-ex_to<add>(elem.rest).seq.size();
				// add2 is for the inner loop and should be the bigger of the two sums
				// in the presence of asymptotically good sorting:
				const add& add1 = (sizedifference<0 ? ex_to<add>(last_expanded) : ex_to<add>(elem.rest));
				const add& add2 = (sizedifference<0 ? ex_to<add>(elem.rest) : ex_to<add>(last_expanded));
				const auto& add1begin = add1.seq.begin();
				const auto& add1end   = add1.seq.end();
				const auto& add2begin = add2.seq.begin();
				const auto& add2end   = add2.seq.end();
				epvector distrseq;
				auto s = add1.seq.size()+add2.seq.size();
//                          // the poly_mul_expand function is buggy and
//                          // therefore cannot be used (see sage :trac:`31478`),
//                          // so we comment out this section of code until the
//                          // function has been fixed
//                                if (s > 400) {
////                                // the condition is probably too simple
//                                        try {
//                                        // can it be converted/expanded via Singular?
//                                        last_expanded = poly_mul_expand(last_expanded,
//                                                                        elem.rest);
//                                                continue;
//                                        }
//                                        catch (std::runtime_error) {
//                                                std::cerr << "can't happen while calling poly_mul_expand\n";
//                                        }
//                                }

				distrseq.reserve(s);
				// Multiply add2 with the overall coefficient of add1 and append it to distrseq:
				if (!add1.overall_coeff.is_zero()) {
					if (add1.overall_coeff.is_one())
						distrseq.insert(distrseq.end(),add2begin,add2end);
					else
                                                for (const auto & elem2 : add2.seq)
							distrseq.emplace_back(elem2.rest,
                                                                ex_to<numeric>(elem2.coeff).mul_dyn(add1.overall_coeff));
				}

				// Multiply add1 with the overall coefficient of add2 and append it to distrseq:
				if (!add2.overall_coeff.is_zero()) {
					if (add2.overall_coeff.is_one())
						distrseq.insert(distrseq.end(),add1begin,add1end);
					else
                                                for (const auto & elem1 : add1.seq)
							distrseq.emplace_back(elem1.rest, ex_to<numeric>(elem1.coeff).mul_dyn(add2.overall_coeff));
				}

				// Compute the new overall coefficient and put it together:
				ex tmp_accu = (new add(distrseq, add1.overall_coeff.mul(add2.overall_coeff)))->setflag(status_flags::dynallocated);

				lst dummy_subs;
				// Multiply explicitly all non-numeric terms of add1 and add2:
                                for (const auto & elem2 : add2.seq) {
					// We really have to combine terms here in order to compactify
					// the result.  Otherwise it would become *way too big*.
					numeric oc(*_num0_p);
					epvector distrseq2;
					distrseq2.reserve(add1.seq.size());
					const ex& i2_new = elem2.rest;
                                        for (const auto & elem1 : add1.seq) {
						// Don't push_back expairs which might have a rest that evaluates to a numeric,
						// since that would violate an invariant of expairseq:
						const ex rest = (new mul(elem1.rest, i2_new))->setflag(status_flags::dynallocated);
						if (is_exactly_a<numeric>(rest)) {
							oc += ex_to<numeric>(rest).mul(ex_to<numeric>(elem1.coeff).mul(ex_to<numeric>(elem2.coeff)));
						} else {
							distrseq2.emplace_back(rest, ex_to<numeric>(elem1.coeff).mul_dyn(ex_to<numeric>(elem2.coeff)));
						}
					}
					tmp_accu += (new add(distrseq2, oc))->setflag(status_flags::dynallocated);
				} 
				last_expanded = tmp_accu;
			} else {
				if (!last_expanded.is_one())
					non_adds.push_back(split_ex_to_pair(last_expanded));
				last_expanded = elem.rest;
			}

		} else {
			non_adds.push_back(elem);
		}
	}

	// Now the only remaining thing to do is to multiply the factors which
	// were not sums into the "last_expanded" sum
	if (is_exactly_a<add>(last_expanded)) {
		exvector distrseq;
		distrseq.reserve(last_expanded.nops());
		exvector va;

		for (const auto term_in_expanded : last_expanded) {
			epvector factors = non_adds;
		        factors.push_back(split_ex_to_pair(term_in_expanded));
			ex term = (new mul(factors, overall_coeff))->setflag(status_flags::dynallocated);
			if (can_be_further_expanded(term)) {
				distrseq.push_back(term.expand());
			} else {
				if (options == 0)
					ex_to<basic>(term).setflag(status_flags::expanded);
				distrseq.push_back(term);
			}
		}

		return ((new add(distrseq))->
		        setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0)));
	}

	non_adds.push_back(split_ex_to_pair(last_expanded));
	ex result = (new mul(non_adds, overall_coeff))->
                setflag(status_flags::dynallocated);
	if (can_be_further_expanded(result)) {
		return result.expand();
	} 
        if (options == 0)
                ex_to<basic>(result).setflag(status_flags::expanded);
        return result;
}

const epvector & mul::get_sorted_seq() const
{
	if (seq_sorted.empty() and not seq.empty()) {
		seq_sorted = epvector(seq.size());
		partial_sort_copy(seq.begin(), seq.end(),
				  seq_sorted.begin(), seq_sorted.end(),
				  print_order_pair_mul());
	}
	return expairseq::get_sorted_seq();
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////


/** Member-wise expand the expairs representing this sequence.  This must be
 *  overridden from expairseq::expandchildren() and done iteratively in order
 *  to allow for early cancellations and thus save memory.
 *
 *  @see mul::expand()
 *  @return pointer to epvector containing expanded representation or zero
 *  pointer, if sequence is unchanged. */
std::unique_ptr<epvector> mul::expandchildren(unsigned options) const
{
	const auto& last = seq.end();
	auto cit = seq.begin();
	while (cit!=last) {
		const ex & factor = recombine_pair_to_ex(*cit);
		const ex & expanded_factor = factor.expand(options);
		if (!are_ex_trivially_equal(factor,expanded_factor)) {
			
			// something changed, copy seq, eval and return it
			std::unique_ptr<epvector> s(new epvector);
			s->reserve(seq.size());
			
			// copy parts of seq which are known not to have changed
			auto cit2 = seq.begin();
			while (cit2!=cit) {
				s->push_back(*cit2);
				++cit2;
			}

			// copy first changed element
			s->push_back(split_ex_to_pair(expanded_factor));
			++cit2;

			// copy rest
			while (cit2!=last) {
				s->push_back(split_ex_to_pair(recombine_pair_to_ex(*cit2).expand(options)));
				++cit2;
			}
			return s;
		}
		++cit;
	}
	
	return std::unique_ptr<epvector>(nullptr); // nothing has changed
}

ex mul::without_known_factor(const ex& f) const
{
	epvector s;
	s.reserve(seq.size()-1);

        bool found = false;
	for (const auto & elem : seq) {
                ex t = recombine_pair_to_ex(elem);
                if (not found and t.is_equal(f)) {
                        found = true;
                }
                else
			s.push_back(elem);
	}

	mul * result = new mul(s, overall_coeff);
	return result->setflag(status_flags::dynallocated);
}

} // namespace GiNaC
