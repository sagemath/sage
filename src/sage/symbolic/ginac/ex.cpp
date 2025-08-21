/** @file ex.cpp
 *
 *  Implementation of GiNaC's light-weight expression handles. */

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

#include <Python.h>
#include "ex.h"
#include "ex_utils.h"
#include "symbol.h"
#include "add.h"
#include "mul.h"
#include "numeric.h"
#include "power.h"
#include "lst.h"
#include "function.h"
#include "fderivative.h"
#include "symbol.h"
#include "relational.h"
#include "utils.h"
#include "operators.h"
#include "wildcard.h"

#include <iostream>
#include <stdexcept>
#include <unordered_map>

namespace GiNaC {

//////////
// other constructors
//////////

// none (all inlined)

//////////
// non-virtual functions in this class
//////////

// public
	
/** Print expression to stream. The formatting of the output is determined
 *  by the kind of print_context object that is passed. Possible formattings
 *  include ginsh-parsable output (the default), tree-like output for
 *  debugging, and C++ source.
 *  @see print_context */
void ex::print(const print_context & c, unsigned level) const
{
	bp->print(c, level);
}

/** Little wrapper around print to be called within a debugger. */
void ex::dbgprint() const
{
	bp->dbgprint();
}

/** Little wrapper around printtree to be called within a debugger. */
void ex::dbgprinttree() const
{
	bp->dbgprinttree();
}

ex ex::expand(unsigned options) const
{
	if (options == 0 && ((bp->flags & status_flags::expanded) != 0u)) // The "expanded" flag only covers the standard options; someone might want to re-expand with different options
		return *this;

        return bp->expand(options);
}

/** Compute partial derivative of an expression.
 *
 *  @param s  symbol by which the expression is derived
 *  @param nth  order of derivative (default 1)
 *  @return partial derivative as a new expression */
ex ex::diff(const symbol & s, unsigned nth) const
{
	if (nth == 0u)
		return *this;

	return bp->diff(s, nth);
}

/** Check whether expression matches a specified pattern. */
bool ex::match(const ex & pattern) const
{
	exmap map;
	return bp->match(pattern, map);
}

bool ex::match(const ex & pattern, lst & repl_lst) const
{
        exmap map;
        bool ret = bp->match(pattern, map);
        for (const auto& pair : map)
                repl_lst.append(pair.first == pair.second);
        return ret;
}

bool ex::match(const ex & pattern, exvector& vec) const
{
        exmap map;
        bool ret = this->match(pattern, map);
        if (not ret)
                return ret;
        unsigned maxl = 0;
        for (const auto& pair : map) {
                if (not is_exactly_a<wildcard>(pair.first))
                        throw std::runtime_error("no wildcard");
                unsigned l = ex_to<wildcard>(pair.first).get_label();
                if (maxl < l)
                        maxl = l;
        }
        ex nan = NaN;
        exvector tvec(maxl+1);
        tvec.assign(maxl+1, nan);
        for (const auto& pair : map) {
                tvec[ex_to<wildcard>(pair.first).get_label()] = pair.second;
        }
        vec = tvec;
        return ret;
}

/** Find all occurrences of a pattern. The found matches are appended to
 *  the "found" list. If the expression itself matches the pattern, the
 *  children are not further examined. This function returns true when any
 *  matches were found. */
bool ex::find(const ex & pattern, lst & found) const
{
	if (match(pattern)) {
		found.append(*this);
		found.sort();
		found.unique();
		return true;
	}
	bool any_found = false;
	for (size_t i=0; i<nops(); i++)
		if (op(i).find(pattern, found))
			any_found = true;
	return any_found;
}

/** Substitute objects in an expression (syntactic substitution) and return
 *  the result as a new expression. */
ex ex::subs(const lst & ls, const lst & lr, unsigned options) const
{
	GINAC_ASSERT(ls.nops() == lr.nops());

	// Convert the lists to a map
	exmap m;
	for (auto its = ls.begin(), itr = lr.begin(); its != ls.end(); ++its, ++itr) {
		m.insert(std::make_pair(*its, *itr));

		// Search for products and powers in the expressions to be substituted
		// (for an optimization in expairseq::subs())
		if (is_exactly_a<mul>(*its) || is_exactly_a<power>(*its))
			options |= subs_options::pattern_is_product;
	}
	if ((options & subs_options::pattern_is_product) == 0u)
		options |= subs_options::pattern_is_not_product;

	return bp->subs(m, options);
}

/** Substitute objects in an expression (syntactic substitution) and return
 *  the result as a new expression.  There are two valid types of
 *  replacement arguments: 1) a relational like object==ex and 2) a list of
 *  relationals lst(object1==ex1,object2==ex2,...). */
ex ex::subs(const ex & e, unsigned options) const
{
	if (e.info(info_flags::relation_equal)) {

		// Argument is a relation: convert it to a map
		exmap m;
		const ex & s = e.op(0);
		m.insert(std::make_pair(s, e.op(1)));

		if (is_exactly_a<mul>(s) || is_exactly_a<power>(s))
			options |= subs_options::pattern_is_product;
		else
			options |= subs_options::pattern_is_not_product;

		return bp->subs(m, options);

	}
        
        if (e.info(info_flags::list)) {

		// Argument is a list: convert it to a map
		exmap m;
		GINAC_ASSERT(is_a<lst>(e));
		for (const auto& r : ex_to<lst>(e)) {
			
			if (!r.info(info_flags::relation_equal))
				throw(std::invalid_argument("basic::subs(ex): argument must be a list of equations"));
			const ex & s = r.op(0);
			m.insert(std::make_pair(s, r.op(1)));

			// Search for products and powers in the expressions to be substituted
			// (for an optimization in expairseq::subs())
			if (is_exactly_a<mul>(s) || is_exactly_a<power>(s))
				options |= subs_options::pattern_is_product;
		}
		if ((options & subs_options::pattern_is_product) == 0u)
			options |= subs_options::pattern_is_not_product;

		return bp->subs(m, options);

	} else
		throw(std::invalid_argument("ex::subs(ex): argument must be a relation_equal or a list"));
}

ex ex::subs(const exmap& m, unsigned options) const
{
        if ((options & subs_options::no_pattern) != 0)
                return bp->subs(m, options);
        for (const auto& p : m)
                if (haswild(p.first))
                        return bp->subs(m, options);

        return bp->subs(m, options | subs_options::no_pattern);
}

/** Traverse expression tree with given visitor, preorder traversal. */
void ex::traverse_preorder(visitor & v) const
{
	accept(v);

	size_t n = nops();
	for (size_t i = 0; i < n; ++i)
		op(i).traverse_preorder(v);
}

/** Traverse expression tree with given visitor, postorder traversal. */
void ex::traverse_postorder(visitor & v) const
{
	size_t n = nops();
	for (size_t i = 0; i < n; ++i)
		op(i).traverse_postorder(v);

	accept(v);
}

/** Return modifiable operand/member at position i. */
ex & ex::let_op(size_t i)
{
	makewriteable();
	return bp->let_op(i);
}

ex & ex::operator[](const ex & index)
{
	makewriteable();
	return (*bp)[index];
}

ex & ex::operator[](size_t i)
{
	makewriteable();
	return (*bp)[i];
}

bool ex::is_equal(const ex & other) const
{
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.total_is_equals++;
#endif
	if (bp == other.bp)  // trivial case: both expressions point to same basic
		return true;
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.nontrivial_is_equals++;
#endif
	if (is_exactly_a<numeric>(*this) and is_exactly_a<numeric>(other))
		return ex_to<numeric>(*this).is_equal(ex_to<numeric>(other));
	const bool equal = bp->is_equal(*other.bp);
#if 0
	if (equal) {
		// Expressions point to different, but equal, trees: conserve
		// memory and make subsequent compare() operations faster by
		// making both expressions point to the same tree.
		share(other);
	}
#endif
	return equal;
}

int ex::compare(const ex & other) const
{
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.total_compares++;
#endif
	if (bp == other.bp)  // trivial case: both expressions point to same basic
		return 0;
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.nontrivial_compares++;
#endif
	if (is_exactly_a<numeric>(*this) and is_exactly_a<numeric>(other))
		return ex_to<numeric>(*this).compare_same_type(ex_to<numeric>(other));
	const int cmpval = bp->compare(*other.bp);
#if 1
	if (cmpval == 0) {
		// Expressions point to different, but equal, trees: conserve
		// memory and make subsequent compare() operations faster by
		// making both expressions point to the same tree.
		share(other);
	}
#endif
	return cmpval;
}

/** Left hand side of relational expression. */
ex ex::lhs() const
{
	if (!is_exactly_a<relational>(*this))
		throw std::runtime_error("ex::lhs(): not a relation");
	return bp->op(0);
}

/** Right hand side of relational expression. */
ex ex::rhs() const
{
	if (!is_exactly_a<relational>(*this))
		throw std::runtime_error("ex::rhs(): not a relation");
	return bp->op(1);
}

/** Check whether expression is a polynomial. */
bool ex::is_polynomial(const ex & vars) const
{
	if (is_exactly_a<lst>(vars)) {
		const lst & varlst = ex_to<lst>(vars);
		for (const auto & elem : varlst)
			if (!bp->is_polynomial(elem))
				return false;
		return true;
	}

	return bp->is_polynomial(vars);
}

bool ex::is_zero() const {
#ifndef _MSC_VER
	  extern const ex _ex0;
#endif
        if (!is_exactly_a<numeric>(*this))
                return is_equal(_ex0);
        return ex_to<numeric>(*this).is_zero();
}

bool ex::is_one() const
{
        if (!is_exactly_a<numeric>(*this))
                return false;
        const numeric& num = ex_to<numeric>(*this);
        return num.is_one();
}
  
bool ex::is_minus_one() const
{
        if (!is_exactly_a<numeric>(*this))
                return false;
        const numeric& num = ex_to<numeric>(*this);
        return num.is_minus_one();
}

bool ex::is_negative_or_minus() const
{
        if (is_exactly_a<mul>(*this)
           and ex_to<mul>(*this).get_overall_coeff().is_negative())
                return true;
        return is_negative();
}

bool ex::is_num_integer() const
{
        if (!is_exactly_a<numeric>(*this))
                return false;
        const numeric& num = ex_to<numeric>(*this);
        return num.is_integer();
}

bool ex::is_num_fraction() const
{
        if (!is_exactly_a<numeric>(*this))
                return false;
        const numeric& num = ex_to<numeric>(*this);
        return num.is_mpq();
}


void ex::set_domain(unsigned d)
{
        if (is_exactly_a<symbol>(*this)) {
                symbol &s = dynamic_cast<symbol&>(*bp);
                s.set_domain(d);
        }
        else if (is_exactly_a<function>(*this)) {
                function &f = dynamic_cast<function&>(*bp);
                f.set_domain(d);
        }
}

static void _treesize(const ex& e, size_t& n)
{
        ++n;
        for (size_t i=0; i < e.nops(); i++)
                _treesize(e.op(i), n);
}

size_t ex::treesize() const
{
        size_t n = 0;
        _treesize(*this, n);
        return n;
}

size_t ex::nsymbols() const
// DON'T USE if you want the number of symbols, instead use symbols().size()
{
	int res = 0;
	if (is_exactly_a<symbol>(*this)) {
		res=1;
	} else {
		for (size_t i=0; i < nops(); i++)
			res += op(i).nsymbols();
	}
	return res;
}

/** Return pointer to first symbol found in expression.  Due to GiNaC's
 *  internal ordering of terms, it may not be obvious which symbol this
 *  function returns for a given expression.
 *
 *  @param e  expression to search
 *  @param x  first symbol found (returned)
 *  @return "false" if no symbol was found, "true" otherwise */
bool ex::get_first_symbol(ex &x) const
{
	if (is_exactly_a<symbol>(*this)) {
		x = *this;
		return true;
	} 
        if (is_exactly_a<add>(*this) || is_exactly_a<mul>(*this)) {
		for (size_t i=0; i<nops(); i++)
			if (sorted_op(i).get_first_symbol(x))
				return true;
	} else if (is_exactly_a<power>(*this)) {
		if (op(0).get_first_symbol(x))
			return true;
	}
	return false;
}

static void collect_symbols(const ex& e, symbolset& syms)
{
	if (is_exactly_a<symbol>(e))
		syms.insert(ex_to<symbol>(e));
	else
		for (size_t i=0; i < e.nops(); i++)
                        collect_symbols(e.op(i), syms);
}

symbolset ex::symbols() const
{
        symbolset the_set;
        collect_symbols(*this, the_set);
	return the_set;
}

static void collect_bound_symbols(const ex& e, symbolset& syms)
{
        static unsigned int sum_serial = function::find_function("sum", 4);
        static unsigned int integral_serial = function::find_function("integrate", 4);
        static unsigned int limit_serial = function::find_function("limit", 0);
	if (is_exactly_a<function>(e)) {
                const function& f = ex_to<function>(e);
                if (f.get_serial() == sum_serial
                    and is_exactly_a<symbol>(f.op(1))) {
                        syms.insert(ex_to<symbol>(f.op(1)));
                        return collect_bound_symbols(f.op(0), syms);
                }
                if (f.get_serial() == integral_serial
                    and is_exactly_a<symbol>(f.op(1))) {
                        syms.insert(ex_to<symbol>(f.op(1)));
                        return collect_bound_symbols(f.op(0), syms);
                }
                if (f.get_serial() == limit_serial
                    and is_exactly_a<symbol>(f.op(1))) {
                        syms.insert(ex_to<symbol>(f.op(1)));
                        return collect_bound_symbols(f.op(0), syms);
                }
        }
        else if (is_exactly_a<fderivative>(e)) {
                const fderivative& d = ex_to<fderivative>(e);
                for (size_t i=0; i<d.nops(); i++)
                        if (is_exactly_a<symbol>(d.op(i)))
                                syms.insert(ex_to<symbol>(d.op(i)));
                return;
        }
	else
		for (size_t i=0; i < e.nops(); i++)
                        collect_bound_symbols(e.op(i), syms);
}

symbolset ex::free_symbols() const
{
        symbolset the_set, bound_set;
        collect_symbols(*this, the_set);
        collect_bound_symbols(*this, bound_set);
        for (auto it = the_set.begin(); it != the_set.end(); )
                if (bound_set.find(*it) != bound_set.end())
                        it = the_set.erase(it);
                else
                        ++it;
	return the_set;
}

static void collect_functions(const ex& e, std::unordered_set<unsigned>& funs)
{
	if (is_exactly_a<function>(e)) {
                const function& f = ex_to<function>(e);
                funs.insert(f.get_serial());
        }
        for (size_t i=0; i < e.nops(); i++)
                collect_functions(e.op(i), funs);
}

std::unordered_set<unsigned> ex::functions() const
{
        std::unordered_set<unsigned> the_set;
        collect_functions(*this, the_set);
        return the_set;
}

static void collect_wilds(const ex& e, wildset& wilds)
{
	if (is_exactly_a<wildcard>(e))
		wilds.insert(ex_to<wildcard>(e));
	else
		for (size_t i=0; i < e.nops(); i++)
                        collect_wilds(e.op(i), wilds);
}

wildset ex::wilds() const
{
        wildset the_set;
        collect_wilds(*this, the_set);
	return the_set;
}

ex ex::sorted_op(size_t i) const
{
	if (is_a<expairseq>(*this))
		return dynamic_cast<const expairseq&>(*bp).stable_op(i);

        return bp->op(i);
}

static bool has_nonposint_power(const ex& x, const symbol& symb)
{
        if (is_exactly_a<power>(x)) {
                const power& p = ex_to<power>(x);
                if (is_exactly_a<add>(p.op(0))
                    and has_symbol(p.op(0), symb)
                    and (not is_exactly_a<numeric>(p.op(1))
                    or not ex_to<numeric>(p.op(1)).is_pos_integer()))
                return true;
        }
        for (size_t i=0; i<x.nops(); ++i)
        if (has_nonposint_power(x.op(i), symb))
                return true;

	return false;
}

// Helper function: Return True if term is a monomial in symb.
// If so, fill vec (applying map at the same time) with (coeff, expo) pairs.
static bool match_monom(const ex& term, const symbol& symb,
                expairvec& vec, const exmap& map)
{
        if (not has_free_symbol(term, symb)) {
                return false;
        }
        if (term.is_equal(symb)) {
                vec.push_back(std::make_pair(_ex1, _ex1));
                return true;
        }
        if (is_exactly_a<power>(term)) {
                const power& p = ex_to<power>(term);
                const ex& expo = p.op(1);
                if (p.op(0).is_equal(symb)) {
                        vec.push_back(std::make_pair(_ex1, expo.subs(map)));
                        return true;
                }

                // of form (...+...)^expo
                // we expand only those with integer exponent
                if (is_exactly_a<numeric>(expo)
                                and has_free_symbol(p.op(0), symb)) {
                        const numeric& ee = ex_to<numeric>(expo);
                        if (ee.is_integer() and ee.to_int() > 1) {
                                expairvec tmpvec;
                                expand(term).coefficients(symb, tmpvec);
                                for (const auto& pair : tmpvec)
                                        vec.push_back(std::make_pair(pair.first.subs(map),
                                                                pair.second.subs(map)));
                                return true;
                        }
                }
        }
        if (is_exactly_a<mul>(term)) {
                // Handle cases C*x, C*x^c
                for (const auto& mterm : term) {
                        if (mterm.is_equal(symb)) {
                                ex oth = ex_to<mul>(term).without_known_factor(symb);
                                if (not has_free_symbol(oth, symb)) {
                                        vec.push_back(std::make_pair(oth.subs(map), _ex1));
                                        return true;
                                }
                        }
                        if (is_exactly_a<power>(mterm)) {
                                const power& p = ex_to<power>(mterm);
                                if (p.op(0).is_equal(symb)) {
                                        ex oth = ex_to<mul>(term).without_known_factor(mterm);
                                        if (not has_free_symbol(oth, symb)) {
                                                vec.push_back(std::make_pair(oth.subs(map), p.op(1).subs(map)));
                                                return true;
                                        }
                                }
                        }
                }
                // Handle C*(...+...)^c
                for (const auto& mterm : term) {
                        if (is_exactly_a<power>(mterm)
                                        and is_exactly_a<numeric>(mterm.op(1))
                                        and has_free_symbol(mterm.op(0), symb)) {
                                numeric ee = ex_to<numeric>(mterm.op(1));
                                if (ee.is_integer() and ee.to_int() > 1) {
                                        expairvec tmpvec;
                                        ex e = expand(term);
                                        if (not is_exactly_a<add>(e))
                                                return false;
                                        e.coefficients(symb, tmpvec);
                                        for (const auto& pair : tmpvec)
                                                vec.push_back(std::make_pair(pair.first.subs(map),
                                                        pair.second.subs(map)));
                                        return true;
                                }
                        }
                        if (is_exactly_a<add>(mterm)
                                        and has_free_symbol(mterm, symb)) {
                                expairvec tmpvec;
                                expand(term).coefficients(symb, tmpvec);
                                for (const auto& pair : tmpvec)
                                        vec.push_back(std::make_pair(pair.first.subs(map),
                                                pair.second.subs(map)));
                                return true;
                        }
                }
        }
        return false;
}

numeric ex::degree(const ex & s) const
{
        return bp->degree(s);
}

numeric ex::ldegree(const ex & s) const
{
        return bp->ldegree(s);
}

static bool _collect_powers(ex& e, ex& repl, bool& collected)
{
        bool is_changed = false;
        if (is_exactly_a<symbol>(e)
            or is_exactly_a<numeric>(e)
            or is_exactly_a<constant>(e))
                return is_changed;
        if (is_a<expairseq>(e)) {
                for (size_t i=0; i<e.nops(); ++i) {
                        ex ee, et = e.op(i);
                        bool c = false;
                        bool changed = _collect_powers(et, ee, c);
                        if (changed or c)
                                is_changed = true;
                        if (c)
                                e.set_epseq_from(i, ee);
                }
        }
        else {
                for (size_t i=0; i<e.nops(); ++i) {
                        ex ee;
                        bool c = false;
                        bool changed = _collect_powers(e.let_op(i), ee, c);
                        if (changed or c)
                                is_changed = true;
                        if (c)
                                e.let_op(i) = ee;
                }
        }

        if (not is_exactly_a<mul>(e))
                return is_changed;
        exmap map;
        for (auto term : e) {
                if (is_exactly_a<power>(term)) {
                        auto it = map.find(term.op(0));
                        if (it == map.end())
                                map.insert(std::make_pair(term.op(0),
                                                          term.op(1)));
                        else {
                                it->second += term.op(1);
                                collected = true;
                        }
                }
                else {
                        auto it = map.find(term);
                        if (it == map.end())
                                map.insert(std::make_pair(term, _ex1));
                        else {
                                it->second += _ex1;
                                collected = true;
                        }
                }
        }
        if (not collected)
                return is_changed;
        exvector vec;
        for (auto p : map)
                vec.push_back(power(p.first, p.second));
        repl = mul(vec);
        return is_changed;
}

// Repeatedly: for any mul in the ex, add exponents of powers
// with the same base
ex ex::collect_powers() const
{
        ex the_ex = *this;
        bool b;
        ex r;
        (void)_collect_powers(the_ex, r, b);
        return b? r:the_ex;
}

/**
 * Return in vec a list of pairs with (coefficient, exponent) of ex
 * when interpreted as polynomial in s (with s any expression).
 *
 * We expand only parts of the expression known to contain free s, in
 * order to keep coefficients as unexpanded as possible.
 */
void ex::coefficients(const ex & s, expairvec & vec) const
{
        vec.clear();
        if (is_exactly_a<add>(s)) {
                vec.push_back(std::make_pair(*this, _ex0));
                return;
        }

        // Replace any s with a new anonymous symbol
        symbol symb;
        exmap submap {{s, symb}}, revmap {{symb, s}};
        ex sub = subs(submap);

        if (is_exactly_a<add>(sub)) {
                const add& addref = ex_to<add>(sub);
                const numeric& oc = addref.get_overall_coeff();
                if (not oc.is_zero())
                        vec.emplace_back(std::make_pair(oc, _ex0));
                for (const auto& term : addref.seq) {
                        ex tmp = addref.recombine_pair_to_ex(term);
                        if (has_nonposint_power(tmp, symb)
                            or not match_monom(tmp, symb, vec, revmap))
                                vec.push_back(std::make_pair(tmp.subs(revmap), _ex0));
                }
        }
        else {
                if (has_nonposint_power(sub, symb)
                    or not match_monom(sub, symb, vec, revmap)) {
                        vec.clear();
                        vec.push_back(std::make_pair(*this, _ex0));
                }
                return;
        }

        // Add: Combine coefficients with the same exponent
        std::sort(vec.begin(), vec.end(), [](const std::pair<ex,ex>& x, const std::pair<ex,ex>& y)
        {
                return x.second < y.second;
        });

        auto tmp_it = vec.end();
        for (auto it = vec.end(); it != vec.begin(); ) {
                --it;
                if ((it->first).is_zero()) {
                        vec.erase(it);
                        tmp_it = vec.end();
                        continue;
                }
                if (tmp_it != vec.end() and (tmp_it->second).is_equal(it->second)) {
                        it->first += tmp_it->first;
                        vec.erase(tmp_it);
                        if ((it->first).is_zero()) {
                                vec.erase(it);
                                tmp_it = vec.end();
                                continue;
                        }
                }
                tmp_it = it;
        }

}

ex ex::lcoeff(const ex & s) const
{
        return coeff(s, degree(s));
}

ex ex::tcoeff(const ex & s) const
{
        return coeff(s, ldegree(s));
}

ex ex::deep_combine_fractions(ex e)
{
        if (is_a<expairseq>(e)) {
                combine_map_function func;
                e = e.map(func);
        }
        else {
                if (is_exactly_a<symbol>(e) or
                                is_exactly_a<constant>(e) or
                                is_exactly_a<numeric>(e))
                        return e;
                
                for (unsigned int i=0; i<e.nops(); ++i) {
                        e.let_op(i) = deep_combine_fractions(e.op(i));
                }
        }

        if (is_exactly_a<add>(e)) {
                ex t = ex_to<add>(e).combine_fractions();
                return t;
        }

        return e;
}

ex ex::combine_fractions(bool deep) const
{
        if (deep)
                return deep_combine_fractions(*this);
        if (is_exactly_a<add>(*this))
                return ex_to<add>(*this).combine_fractions();
        return *this;
}

// private

/** Make this ex writable (if more than one ex handle the same basic) by 
 *  unlinking the object and creating an unshared copy of it. */
void ex::makewriteable()
{
	GINAC_ASSERT(bp->flags & status_flags::dynallocated);
	bp.makewritable();
	GINAC_ASSERT(bp->get_refcount() == 1);
}

/** Share equal objects between expressions.
 *  @see ex::compare(const ex &) */
void ex::share(const ex & other) const
{
	if (((bp->flags | other.bp->flags) & status_flags::not_shareable) != 0u)
		return;

	if (bp->get_refcount() <= other.bp->get_refcount())
		bp = other.bp;
	else
		other.bp = bp;
}

/** Helper function for the ex-from-basic constructor. This is where GiNaC's
 *  automatic evaluator and memory management are implemented.
 *  @see ex::ex(const basic &) */
ptr<basic> ex::construct_from_basic(const basic & other)
{
	if (not other.is_evaluated()) {

		// The object is not yet evaluated, so call eval() to evaluate
		// the top level. This will return either
		//  a) the original object with status_flags::evaluated set (when the
		//     eval() implementation calls hold())
		// or
		//  b) a different expression.
		//
		// eval() returns an ex, not a basic&, so this will go through
		// construct_from_basic() a second time. In case a) we end up in
		// the "else" branch below. In case b) we end up here again and
		// apply eval() once more. The recursion stops when eval() calls
		// hold() or returns an object that already has its "evaluated"
		// flag set, such as a symbol or a numeric.
		const ex & tmpex = other.eval(1);

		// Eventually, the eval() recursion goes through the "else" branch
		// below, which assures that the object pointed to by tmpex.bp is
		// allocated on the heap (either it was already on the heap or it
		// is a heap-allocated duplicate of another object).
		GINAC_ASSERT(tmpex.bp->flags & status_flags::dynallocated); 

		// If the original object is not referenced but heap-allocated,
		// it means that eval() hit case b) above. The original object is
		// no longer needed (it evaluated into something different), so we
		// delete it (because nobody else will).
		if ((other.get_refcount() == 0) && ((other.flags & status_flags::dynallocated) != 0u))
			delete &other; // yes, you can apply delete to a const pointer

		// We can't return a basic& here because the tmpex is destroyed as
		// soon as we leave the function, which would deallocate the
		// evaluated object.
		return tmpex.bp;
	} 

        // The easy case: making an "ex" out of an evaluated object.
        if ((other.flags & status_flags::dynallocated) != 0u) {

                // The object is already heap-allocated, so we can just make
                // another reference to it.
                return ptr<basic>(const_cast<basic &>(other));
        } 

        // The object is not heap-allocated, so we create a duplicate
        // on the heap.
        basic *bp = other.duplicate();
        bp->setflag(status_flags::dynallocated);
        GINAC_ASSERT(bp->get_refcount() == 0);
        return bp;
}

basic & ex::construct_from_int(int i)
{
	switch (i) {  // prefer flyweights over new objects
	case -12:
		return *const_cast<numeric *>(_num_12_p);
	case -11:
		return *const_cast<numeric *>(_num_11_p);
	case -10:
		return *const_cast<numeric *>(_num_10_p);
	case -9:
		return *const_cast<numeric *>(_num_9_p);
	case -8:
		return *const_cast<numeric *>(_num_8_p);
	case -7:
		return *const_cast<numeric *>(_num_7_p);
	case -6:
		return *const_cast<numeric *>(_num_6_p);
	case -5:
		return *const_cast<numeric *>(_num_5_p);
	case -4:
		return *const_cast<numeric *>(_num_4_p);
	case -3:
		return *const_cast<numeric *>(_num_3_p);
	case -2:
		return *const_cast<numeric *>(_num_2_p);
	case -1:
		return *const_cast<numeric *>(_num_1_p);
	case 0:
		return *const_cast<numeric *>(_num0_p);
	case 1:
		return *const_cast<numeric *>(_num1_p);
	case 2:
		return *const_cast<numeric *>(_num2_p);
	case 3:
		return *const_cast<numeric *>(_num3_p);
	case 4:
		return *const_cast<numeric *>(_num4_p);
	case 5:
		return *const_cast<numeric *>(_num5_p);
	case 6:
		return *const_cast<numeric *>(_num6_p);
	case 7:
		return *const_cast<numeric *>(_num7_p);
	case 8:
		return *const_cast<numeric *>(_num8_p);
	case 9:
		return *const_cast<numeric *>(_num9_p);
	case 10:
		return *const_cast<numeric *>(_num10_p);
	case 11:
		return *const_cast<numeric *>(_num11_p);
	case 12:
		return *const_cast<numeric *>(_num12_p);
	default:
		basic *bp = new numeric(i);
		bp->setflag(status_flags::dynallocated);
		GINAC_ASSERT(bp->get_refcount() == 0);
		return *bp;
	}
}
	
basic & ex::construct_from_uint(unsigned int i)
{
	switch (i) {  // prefer flyweights over new objects
	case 0:
		return *const_cast<numeric *>(_num0_p);
	case 1:
		return *const_cast<numeric *>(_num1_p);
	case 2:
		return *const_cast<numeric *>(_num2_p);
	case 3:
		return *const_cast<numeric *>(_num3_p);
	case 4:
		return *const_cast<numeric *>(_num4_p);
	case 5:
		return *const_cast<numeric *>(_num5_p);
	case 6:
		return *const_cast<numeric *>(_num6_p);
	case 7:
		return *const_cast<numeric *>(_num7_p);
	case 8:
		return *const_cast<numeric *>(_num8_p);
	case 9:
		return *const_cast<numeric *>(_num9_p);
	case 10:
		return *const_cast<numeric *>(_num10_p);
	case 11:
		return *const_cast<numeric *>(_num11_p);
	case 12:
		return *const_cast<numeric *>(_num12_p);
	default:
		basic *bp = new numeric(i);
		bp->setflag(status_flags::dynallocated);
		GINAC_ASSERT(bp->get_refcount() == 0);
		return *bp;
	}
}
	
basic & ex::construct_from_long(long i)
{
	switch (i) {  // prefer flyweights over new objects
	case -12:
		return *const_cast<numeric *>(_num_12_p);
	case -11:
		return *const_cast<numeric *>(_num_11_p);
	case -10:
		return *const_cast<numeric *>(_num_10_p);
	case -9:
		return *const_cast<numeric *>(_num_9_p);
	case -8:
		return *const_cast<numeric *>(_num_8_p);
	case -7:
		return *const_cast<numeric *>(_num_7_p);
	case -6:
		return *const_cast<numeric *>(_num_6_p);
	case -5:
		return *const_cast<numeric *>(_num_5_p);
	case -4:
		return *const_cast<numeric *>(_num_4_p);
	case -3:
		return *const_cast<numeric *>(_num_3_p);
	case -2:
		return *const_cast<numeric *>(_num_2_p);
	case -1:
		return *const_cast<numeric *>(_num_1_p);
	case 0:
		return *const_cast<numeric *>(_num0_p);
	case 1:
		return *const_cast<numeric *>(_num1_p);
	case 2:
		return *const_cast<numeric *>(_num2_p);
	case 3:
		return *const_cast<numeric *>(_num3_p);
	case 4:
		return *const_cast<numeric *>(_num4_p);
	case 5:
		return *const_cast<numeric *>(_num5_p);
	case 6:
		return *const_cast<numeric *>(_num6_p);
	case 7:
		return *const_cast<numeric *>(_num7_p);
	case 8:
		return *const_cast<numeric *>(_num8_p);
	case 9:
		return *const_cast<numeric *>(_num9_p);
	case 10:
		return *const_cast<numeric *>(_num10_p);
	case 11:
		return *const_cast<numeric *>(_num11_p);
	case 12:
		return *const_cast<numeric *>(_num12_p);
	default:
		basic *bp = new numeric(i);
		bp->setflag(status_flags::dynallocated);
		GINAC_ASSERT(bp->get_refcount() == 0);
		return *bp;
	}
}
	
basic & ex::construct_from_double(double d)
{
	basic *bp = new numeric(d);
	bp->setflag(status_flags::dynallocated);
	GINAC_ASSERT(bp->get_refcount() == 0);
	return *bp;
}

basic & ex::construct_from_pyobject(PyObject* o)
{
  Py_INCREF(o);   // since numeric steals a reference
                  // note: unlike `numeric(o)`, `ex(o)` does not steal a
                  // reference to `o`. You can explicitly move the reference
                  // to the newly-constructed `ex` object by `ex(numeric(o))`.
  basic *bp = new numeric(o);
  bp->setflag(status_flags::dynallocated);
  GINAC_ASSERT(bp->get_refcount() == 0);
  return *bp;
}


 
	
//////////
// static member variables
//////////

// none

//////////
// functions which are not member functions
//////////

// none

//////////
// global functions
//////////

// none

//
} // namespace GiNaC
