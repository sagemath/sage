/** @file pseries.cpp
 *
 *  Implementation of class for extended truncated power series and
 *  methods for series expansion. */

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

#include "pseries.h"
#include "useries.h"
#include "useries-flint.h"
#include "add.h"
#include "inifcns.h" // for Order function
#include "lst.h"
#include "mul.h"
#include "power.h"
#include "relational.h"
#include "operators.h"
#include "symbol.h"
#include "archive.h"
#include "utils.h"
#include "infinity.h"

#include <numeric>
#include <stdexcept>
#include <limits>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(pseries, basic,
  print_func<print_context>(&pseries::do_print).
  print_func<print_latex>(&pseries::do_print_latex).
  print_func<print_tree>(&pseries::do_print_tree).
  print_func<print_python>(&pseries::do_print_python).
  print_func<print_python_repr>(&pseries::do_print_python_repr))


/*
 *  Default constructor
 */

pseries::pseries() : inherited(&pseries::tinfo_static) { }


/*
 *  Other ctors
 */

/** Construct pseries from a vector of coefficients and powers.
 *  expair.rest holds the coefficient, expair.coeff holds the power.
 *  The powers must be integers (positive or negative) and in ascending order;
 *  the last coefficient can be Order(_ex1) to represent a truncated,
 *  non-terminating series.
 *
 *  @param rel_  expansion variable and point (must hold a relational)
 *  @param ops_  vector of {coefficient, power} pairs (coefficient must not be zero)
 *  @return newly constructed pseries */
pseries::pseries(const ex &rel_, epvector ops_) : basic(&pseries::tinfo_static), seq(std::move(ops_))
{
	GINAC_ASSERT(is_exactly_a<relational>(rel_));
	GINAC_ASSERT(is_exactly_a<symbol>(rel_.lhs()));
	point = rel_.rhs();
	var = rel_.lhs();
}


/*
 *  Archiving
 */

pseries::pseries(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	auto first = n.find_first("coeff");
	auto last = n.find_last("power");
	++last;
	seq.reserve((last-first)/2);

	for (auto loc = first; loc < last;) {
		ex rest;
		ex coef;
		n.find_ex_by_loc(loc++, rest, sym_lst);
		n.find_ex_by_loc(loc++, coef, sym_lst);
		seq.emplace_back(rest, coef);
	}

	n.find_ex("var", var, sym_lst);
	n.find_ex("point", point, sym_lst);
}

void pseries::archive(archive_node &n) const
{
	inherited::archive(n);
        for (const auto & elem : seq) {
		n.add_ex("coeff", elem.rest);
		n.add_ex("power", elem.coeff);
	}
	n.add_ex("var", var);
	n.add_ex("point", point);
}

//////////
// functions overriding virtual functions from base classes
//////////

void pseries::print_series(const print_context & c, const char *openbrace, const char *closebrace, const char *mul_sym, const char *pow_sym, unsigned level) const
{
	if (precedence() <= level)
		c.s << '(';
		
	// objects of type pseries must not have any zero entries, so the
	// trivial (zero) pseries needs a special treatment here:
	if (seq.empty())
		c.s << '0';

	auto i = seq.begin(), end = seq.end();
	while (i != end) {

		// print a sign, if needed
		if (i != seq.begin())
			c.s << " + ";

		if (!is_order_function(i->rest)) {

			// print 'rest', i.e. the expansion coefficient
			if (i->rest.info(info_flags::numeric) &&
				i->rest.is_positive()) {
				i->rest.print(c);
			} else {
				c.s << openbrace << '(';
				i->rest.print(c);
				c.s << ')' << closebrace;
			}

			// print 'coeff', something like (x-1)^42
			if (!i->coeff.is_zero()) {
				c.s << mul_sym;
				if (!point.is_zero()) {
					c.s << openbrace << '(';
					(var-point).print(c);
					c.s << ')' << closebrace;
				} else
					var.print(c);
				if (i->coeff.compare(_ex1) != 0) {
					c.s << pow_sym;
					c.s << openbrace;
					if (i->coeff.info(info_flags::negative)) {
						c.s << '(';
						i->coeff.print(c);
						c.s << ')';
					} else
						i->coeff.print(c);
					c.s << closebrace;
				}
			}
		} else
			Order(power(var-point,i->coeff)).print(c);
		++i;
	}

	if (precedence() <= level)
		c.s << ')';
}

void pseries::do_print(const print_context & c, unsigned level) const
{
	print_series(c, "", "", "*", "^", level);
}

void pseries::do_print_latex(const print_latex & c, unsigned level) const
{
	print_series(c, "{", "}", " ", "^", level);
}

void pseries::do_print_python(const print_python & c, unsigned level) const
{
	print_series(c, "", "", "*", "**", level);
}

void pseries::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << std::endl;
	size_t num = seq.size();
	for (size_t i=0; i<num; ++i) {
		seq[i].rest.print(c, level + c.delta_indent);
		seq[i].coeff.print(c, level + c.delta_indent);
		c.s << std::string(level + c.delta_indent, ' ') << "-----" << std::endl;
	}
	var.print(c, level + c.delta_indent);
	point.print(c, level + c.delta_indent);
}

void pseries::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << "(relational(";
	var.print(c);
	c.s << ',';
	point.print(c);
	c.s << "),[";
	size_t num = seq.size();
	for (size_t i=0; i<num; ++i) {
		if (i != 0u)
			c.s << ',';
		c.s << '(';
		seq[i].rest.print(c);
		c.s << ',';
		seq[i].coeff.print(c);
		c.s << ')';
	}
	c.s << "])";
}

int pseries::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<pseries>(other));
	const pseries &o = static_cast<const pseries &>(other);
	
	// first compare the lengths of the series...
	if (seq.size()>o.seq.size())
		return 1;
	if (seq.size()<o.seq.size())
		return -1;
	
	// ...then the expansion point...
	int cmpval = var.compare(o.var);
	if (cmpval != 0)
		return cmpval;
	cmpval = point.compare(o.point);
	if (cmpval != 0)
		return cmpval;
	
	// ...and if that failed the individual elements
	auto it = seq.begin(), o_it = o.seq.begin();
	while (it!=seq.end() && o_it!=o.seq.end()) {
		cmpval = it->compare(*o_it);
		if (cmpval != 0)
			return cmpval;
		++it;
		++o_it;
	}

	// so they are equal.
	return 0;
}

/** Return the number of operands including a possible order term. */
size_t pseries::nops() const
{
	return seq.size();
}

/** Return the ith term in the series when represented as a sum. */
const ex pseries::op(size_t i) const
{
	if (i >= seq.size())
		throw (std::out_of_range("op() out of range"));

	if (is_order_function(seq[i].rest))
		return Order(power(var-point, seq[i].coeff));
	return seq[i].rest * power(var - point, seq[i].coeff);
}

/** Return degree of highest power of the series.  This is usually the exponent
 *  of the Order term.  If s is not the expansion variable of the series, the
 *  series is examined termwise. */
numeric pseries::degree(const ex &s) const
{
	if (var.is_equal(s)) {
		// Return last exponent
		if (!seq.empty())
			return ex_to<numeric>((seq.end()-1)->coeff).to_int();
		
			return 0;
	} else {
                
		if (seq.empty())
			return 0;
		int max_pow = std::numeric_limits<int>::min();
                for (const auto & elem : seq) {
			int pow = elem.rest.degree(s).to_int();
			if (pow > max_pow)
				max_pow = pow;
		}
		return max_pow;
	}
}

/** Return degree of lowest power of the series.  This is usually the exponent
 *  of the leading term.  If s is not the expansion variable of the series, the
 *  series is examined termwise.  If s is the expansion variable but the
 *  expansion point is not zero the series is not expanded to find the degree.
 *  I.e.: (1-x) + (1-x)^2 + Order((1-x)^3) has ldegree(x) 1, not 0. */
numeric pseries::ldegree(const ex &s) const
{
	if (var.is_equal(s)) {
		// Return first exponent
		if (!seq.empty())
			return ex_to<numeric>((seq.begin())->coeff);
		
			return 0;
	} else {
		if (seq.empty())
			return 0;
		int min_pow = std::numeric_limits<int>::max();
                for (const auto & elem : seq) {
			int pow = elem.rest.ldegree(s).to_int();
			if (pow < min_pow)
				min_pow = pow;
		}
		return min_pow;
	}
}

/** Return coefficient of degree n in power series if s is the expansion
 *  variable.  If the expansion point is nonzero, by definition the n=1
 *  coefficient in s of a+b*(s-z)+c*(s-z)^2+Order((s-z)^3) is b (assuming
 *  the expansion took place in the s in the first place).
 *  If s is not the expansion variable, an attempt is made to convert the
 *  series to a polynomial and return the corresponding coefficient from
 *  there. */
ex pseries::coeff(const ex &s, const ex & n) const
{
	if (var.is_equal(s)) {
		if (seq.empty())
			return _ex0;
		
		// Binary search in sequence for given power
                if (not is_exactly_a<numeric>(n))
                        throw std::runtime_error("can't happen in pseries::coeff");
		const numeric& looking_for = ex_to<numeric>(n);
		int lo = 0, hi = seq.size() - 1;
		while (lo <= hi) {
			int mid = (lo + hi) / 2;
			GINAC_ASSERT(is_exactly_a<numeric>(seq[mid].coeff));
			int cmp = ex_to<numeric>(seq[mid].coeff).compare(looking_for);
                        switch (cmp) {
                        case -1:
                                lo = mid + 1;
                                break;
                        case 0:
                                return seq[mid].rest;
                        case 1:
                                hi = mid - 1;
                                break;
                        default:
                                throw(std::logic_error("pseries::coeff: compare() didn't return -1, 0 or 1"));
                        }
                }
		return _ex0;
	} 
		return convert_to_poly().coeff(s, n);
}

/** Does nothing. */
ex pseries::collect(const ex &s, bool distributed) const
{
	return *this;
}

/** Perform coefficient-wise automatic term rewriting rules in this class. */
ex pseries::eval(int level) const
{
	if (level == 1)
		return this->hold();
	
	if (level == -max_recursion_level)
		throw (std::runtime_error("pseries::eval(): recursion limit exceeded"));
	
	// Construct a new series with evaluated coefficients
	epvector new_seq;
	new_seq.reserve(seq.size());
        for (const auto & elem : seq)
		new_seq.emplace_back(elem.rest.eval(level-1), elem.coeff);
	return (new pseries(relational(var,point), new_seq))->setflag(status_flags::dynallocated | status_flags::evaluated);
}

/** Evaluate coefficients numerically. */
ex pseries::evalf(int level, PyObject* parent) const
{
	if (level == 1)
		return *this;
	
	if (level == -max_recursion_level)
		throw (std::runtime_error("pseries::evalf(): recursion limit exceeded"));
	
	// Construct a new series with evaluated coefficients
	epvector new_seq;
	new_seq.reserve(seq.size());
        for (const auto & elem : seq)
		new_seq.emplace_back(elem.rest.evalf(level-1, parent),
					elem.coeff);
	return (new pseries(relational(var,point), new_seq))->setflag(status_flags::dynallocated | status_flags::evaluated);
}

ex pseries::conjugate() const
{
	if (not var.is_real())
		return conjugate_function(*this).hold();

	std::unique_ptr<epvector> newseq(conjugateepvector(seq));
	ex newpoint = point.conjugate();

	if (!newseq && are_ex_trivially_equal(point, newpoint)) {
		return *this;
	}

	ex result = (new pseries(var==newpoint, newseq ? *newseq : seq))->setflag(status_flags::dynallocated);
	return result;
}

ex pseries::real_part() const
{
	if (not var.is_real())
		return real_part_function(*this).hold();
	ex newpoint = point.real_part();
	if(newpoint != point)
		return real_part_function(*this).hold();

	epvector v;
	v.reserve(seq.size());
	for(const auto & elem : seq)
		v.emplace_back((elem.rest).real_part(), elem.coeff);
	return (new pseries(var==point, v))->setflag(status_flags::dynallocated);
}

ex pseries::imag_part() const
{
	if (not var.is_real())
		return imag_part_function(*this).hold();
	ex newpoint = point.real_part();
	if(newpoint != point)
		return imag_part_function(*this).hold();

	epvector v;
	v.reserve(seq.size());
	for(const auto & elem : seq)
		v.emplace_back((elem.rest).imag_part(), elem.coeff);
	return (new pseries(var==point, v))->setflag(status_flags::dynallocated);
}

ex pseries::subs(const exmap & m, unsigned options) const
{
	// If expansion variable is being substituted, convert the series to a
	// polynomial and do the substitution there because the result might
	// no longer be a power series
	if (m.find(var) != m.end())
		return convert_to_poly(true).subs(m, options);
	
	// Otherwise construct a new series with substituted coefficients and
	// expansion point
	epvector newseq;
	newseq.reserve(seq.size());
        for (const auto & elem : seq)
		newseq.emplace_back(elem.rest.subs(m, options), elem.coeff);
	return (new pseries(relational(var,point.subs(m, options)), newseq))->setflag(status_flags::dynallocated);
}

/** Implementation of ex::expand() for a power series.  It expands all the
 *  terms individually and returns the resulting series as a new pseries. */
ex pseries::expand(unsigned options) const
{
	epvector newseq;
        for (const auto & elem : seq) {
		ex restexp = elem.rest.expand();
		if (!restexp.is_zero())
			newseq.emplace_back(restexp, elem.coeff);
	}
	return (new pseries(relational(var,point), newseq))
	        ->setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
}

/** Implementation of ex::diff() for a power series.
 *  @see ex::diff */
ex pseries::derivative(const symbol & s) const
{
	epvector new_seq;

	if (s == var) {
		
		// FIXME: coeff might depend on var
                for (const auto & elem : seq) {
			if (is_order_function(elem.rest)) {
				new_seq.emplace_back(elem.rest, elem.coeff - 1);
			} else {
				const ex& c = elem.rest * elem.coeff;
				if (!c.is_zero())
					new_seq.emplace_back(c, elem.coeff - 1);
			}
		}

	} else {

                for (const auto & elem : seq) {
			if (is_order_function(elem.rest)) {
				new_seq.push_back(elem);
			} else {
				const ex& c = elem.rest.diff(s);
				if (!c.is_zero())
					new_seq.emplace_back(c, elem.coeff);
			}
		}
	}

	return pseries(relational(var,point), new_seq);
}

ex pseries::convert_to_poly(bool no_order) const
{
	ex e;
        for (const auto & elem : seq) {
		if (is_order_function(elem.rest)) {
			if (!no_order)
				e += Order(power(var - point, elem.coeff));
		}
                else
			e += elem.rest * power(var - point, elem.coeff);
	}
	return e;
}

bool pseries::is_terminating() const
{
	return seq.empty() || !is_order_function((seq.end()-1)->rest);
}

ex pseries::coeffop(size_t i) const
{
	if (i >=nops())
		throw (std::out_of_range("coeffop() out of range"));
	return seq[i].rest;
}

ex pseries::exponop(size_t i) const
{
	if (i >= nops())
		throw (std::out_of_range("exponop() out of range"));
	return seq[i].coeff;
}


/*
 *  Implementations of series expansion
 */

/** Default implementation of ex::series(). This performs Taylor expansion.
 *  @see ex::series */
ex basic::series(const relational & r, int order, unsigned options) const
{
	epvector seq;
	const symbol &s = ex_to<symbol>(r.lhs());

	// default for order-values that make no sense for Taylor expansion
	if ((order <= 0) && this->has(s)) {
		seq.emplace_back(Order(_ex1), order);
		return pseries(r, seq);
	}

	// do Taylor expansion
	numeric fac = 1;
	ex deriv = *this;
	const ex co = deriv.subs(r, subs_options::no_pattern);
	if (not co.is_zero()) {
		seq.emplace_back(co, _ex0);
	}

	int n;
	for (n=1; n<order; ++n) {
		fac = fac.mul(n);
		// We need to test for zero in order to see if the series terminates.
		// The problem is that there is no such thing as a perfect test for
		// zero.  Expanding the term occasionally helps a little...
		deriv = deriv.diff(s).expand(expand_options::expand_only_numerators);
		if (deriv.is_zero())  // Series terminates
			return pseries(r, seq);

		const ex coef = deriv.subs(r, subs_options::no_pattern);
		if (!coef.is_zero())
			seq.emplace_back(fac.inverse() * coef, n);
	}
	
	// Higher-order terms, if present
	deriv = deriv.diff(s);
        if (!deriv.expand(expand_options::expand_only_numerators).is_zero())
		seq.emplace_back(Order(_ex1), n);
	return pseries(r, seq);
}


ex numeric::series(const relational & r, int order, unsigned options) const
{
	epvector seq;
        if (not is_zero())
                seq.emplace_back(*this, _ex0);
        seq.emplace_back(Order(_ex1), numeric(order));
	return pseries(r, seq);
}

/** Implementation of ex::series() for symbols.
 *  @see ex::series */
ex symbol::series(const relational & r, int order, unsigned options) const
{
	epvector seq;
	const ex point = r.rhs();
	GINAC_ASSERT(is_exactly_a<symbol>(r.lhs()));

	if (this->is_equal_same_type(ex_to<symbol>(r.lhs()))) {
		if (order > 0 && !point.is_zero())
			seq.emplace_back(point, _ex0);
		if (order > 1)
			seq.emplace_back(_ex1, _ex1);
		else
			seq.emplace_back(Order(_ex1), numeric(order));
	} else
		seq.emplace_back(*this, _ex0);
	return pseries(r, seq);
}


/** Add one series object to another, producing a pseries object that
 *  represents the sum.
 *
 *  @param other  pseries object to add with
 *  @return the sum as a pseries */
ex pseries::add_series(const pseries &other) const
{
	// Adding two series with different variables or expansion points
	// results in an empty (constant) series 
	if (!is_compatible_to(other)) {
		epvector nul;
		nul.emplace_back(Order(_ex1), _ex0);
		return pseries(relational(var,point), nul);
	}
	
	// Series addition
	epvector new_seq;
	auto a = seq.begin();
	auto b = other.seq.begin();
	auto a_end = seq.end();
	auto b_end = other.seq.end();
	int pow_a = std::numeric_limits<int>::max(), pow_b = std::numeric_limits<int>::max();
	for (;;) {
		// If a is empty, fill up with elements from b and stop
		if (a == a_end) {
			while (b != b_end) {
				new_seq.push_back(*b);
				++b;
			}
			break;
		} 
			pow_a = ex_to<numeric>((*a).coeff).to_int();
		
		// If b is empty, fill up with elements from a and stop
		if (b == b_end) {
			while (a != a_end) {
				new_seq.push_back(*a);
				++a;
			}
			break;
		} 
			pow_b = ex_to<numeric>((*b).coeff).to_int();
		
		// a and b are non-empty, compare powers
		if (pow_a < pow_b) {
			// a has lesser power, get coefficient from a
			new_seq.push_back(*a);
			if (is_order_function((*a).rest))
				break;
			++a;
		} else if (pow_b < pow_a) {
			// b has lesser power, get coefficient from b
			new_seq.push_back(*b);
			if (is_order_function((*b).rest))
				break;
			++b;
		} else {
			// Add coefficient of a and b
			if (is_order_function((*a).rest) || is_order_function((*b).rest)) {
				new_seq.emplace_back(Order(_ex1), (*a).coeff);
				break;  // Order term ends the sequence
			} 
				ex sum = (*a).rest + (*b).rest;
				if (!(sum.is_zero()))
					new_seq.emplace_back(sum, numeric(pow_a));
				++a;
				++b;
			
		}
	}
	return pseries(relational(var,point), new_seq);
}


/** Implementation of ex::series() for sums. This performs series addition when
 *  adding pseries objects.
 *  @see ex::series */
ex add::series(const relational & r, int order, unsigned options) const
{
	ex acc; // Series accumulator
	
	// Get first term from overall_coeff
	acc = overall_coeff.series(r, order, options);
	
	// Add remaining terms
	for (const auto & elem : seq) {
		ex term;
		if (is_exactly_a<pseries>(elem.rest))
			term = elem.rest;
		else
			term = elem.rest.series(r, order, options);
		if (!elem.coeff.is_equal(_ex1))
			term = ex_to<pseries>(term).mul_const(ex_to<numeric>(elem.coeff));
		
		// Series addition
		acc = ex_to<pseries>(acc).add_series(ex_to<pseries>(term));
	}
	return acc;
}


/** Multiply a pseries object with a numeric constant, producing a pseries
 *  object that represents the product.
 *
 *  @param other  constant to multiply with
 *  @return the product as a pseries */
ex pseries::mul_const(const numeric &other) const
{
	epvector new_seq;
	new_seq.reserve(seq.size());
	
        for (const auto & elem : seq) {
		if (!is_order_function(elem.rest))
			new_seq.emplace_back(elem.rest * other, elem.coeff);
		else
			new_seq.push_back(elem);
	}
	return pseries(relational(var,point), new_seq);
}


/** Multiply one pseries object to another, producing a pseries object that
 *  represents the product.
 *
 *  @param other  pseries object to multiply with
 *  @return the product as a pseries */
ex pseries::mul_series(const pseries &other) const
{
	// Multiplying two series with different variables or expansion points
	// results in an empty (constant) series 
	if (!is_compatible_to(other)) {
		epvector nul;
		nul.emplace_back(Order(_ex1), _ex0);
		return pseries(relational(var,point), nul);
	}

	if (seq.empty() || other.seq.empty()) {
		return (new pseries(var==point, epvector()))
		       ->setflag(status_flags::dynallocated);
	}
	
	// Series multiplication
	epvector new_seq;
	const int a_max = degree(var).to_int();
	const int b_max = other.degree(var).to_int();
	const int a_min = ldegree(var).to_int();
	const int b_min = other.ldegree(var).to_int();
	const int cdeg_min = a_min + b_min;
	int cdeg_max = a_max + b_max;
	
	int higher_order_a = std::numeric_limits<int>::max();
	int higher_order_b = std::numeric_limits<int>::max();
	if (is_order_function(coeff(var, a_max)))
		higher_order_a = a_max + b_min;
	if (is_order_function(other.coeff(var, b_max)))
		higher_order_b = b_max + a_min;
	const int higher_order_c = std::min(higher_order_a, higher_order_b);
	if (cdeg_max >= higher_order_c)
		cdeg_max = higher_order_c - 1;

        std::map<int, ex> map_int_ex1, map_int_ex2;
        for (const auto& elem : seq)
                map_int_ex1[ex_to<numeric>(elem.coeff).to_int()] = elem.rest;

        if (other.var.is_equal(var))
                for (const auto& elem : other.seq)
                        map_int_ex2[ex_to<numeric>(elem.coeff).to_int()] = elem.rest;

	for (int cdeg=cdeg_min; cdeg<=cdeg_max; ++cdeg) {
		ex co = _ex0;
		// c(i)=a(0)b(i)+...+a(i)b(0)
		for (int i=a_min; cdeg-i>=b_min; ++i) {
                        const auto& it1 = map_int_ex1.find(i);
                        if (it1 == map_int_ex1.end())
                                continue;
                        const auto& it2 = map_int_ex2.find(cdeg-i);
                        if (it2 == map_int_ex2.end())
                                continue;
			if (!is_order_function(it1->second) && !is_order_function(it2->second))
				co += it1->second * it2->second;
		}
		if (!co.is_zero())
			new_seq.emplace_back(co, numeric(cdeg));
	}
	if (higher_order_c < std::numeric_limits<int>::max())
		new_seq.emplace_back(Order(_ex1), numeric(higher_order_c));
	return pseries(relational(var, point), new_seq);
}


/** Implementation of ex::series() for product. This performs series
 *  multiplication when multiplying series.
 *  @see ex::series */
ex mul::series(const relational & r, int order, unsigned options) const
{
	pseries acc; // Series accumulator

	GINAC_ASSERT(is_exactly_a<symbol>(r.lhs()));
	const ex& sym = r.lhs();
		
	// holds ldegrees of the series of individual factors
	std::vector<int> ldegrees;
	std::vector<bool> ldegree_redo;

	// find minimal degrees
	// first round: obtain a bound up to which minimal degrees have to be
	// considered
        for (const auto & elem : seq) {

		ex expon = elem.coeff;
		int factor = 1;
		ex buf;
		if (expon.is_integer()) {
			buf = elem.rest;
			factor = ex_to<numeric>(expon).to_int();
		} else {
			buf = recombine_pair_to_ex(elem);
		}

		int real_ldegree = 0;
		bool flag_redo = false;
		try {
			real_ldegree = buf.expand().ldegree(sym-r.rhs()).to_int();
		} catch (std::runtime_error) {}

		if (real_ldegree == 0) {
			if ( factor < 0 ) {
				// This case must terminate, otherwise we would have division by
				// zero.
				int orderloop = 0;
				do {
					orderloop++;
					real_ldegree = buf.series(r, orderloop, options).ldegree(sym).to_int();
				} while (real_ldegree == orderloop);
			} else {
				// Here it is possible that buf does not have a ldegree, therefore
				// check only if ldegree is negative, otherwise reconsider the case
				// in the second round.
				real_ldegree = buf.series(r, 0, options).ldegree(sym).to_int();
				if (real_ldegree == 0)
					flag_redo = true;
			}
		}

		ldegrees.push_back(factor * real_ldegree);
		ldegree_redo.push_back(flag_redo);
	}

	int degbound = order-std::accumulate(ldegrees.begin(), ldegrees.end(), 0);
	// Second round: determine the remaining positive ldegrees by the series
	// method.
	// here we can ignore ldegrees larger than degbound
	size_t j = 0;
        for (const auto & elem : seq) {
		if ( ldegree_redo[j] ) {
			ex expon = elem.coeff;
			int factor = 1;
			ex buf;
			if (expon.is_integer()) {
				buf = elem.rest;
				factor = ex_to<numeric>(expon).to_int();
			} else {
				buf = recombine_pair_to_ex(elem);
			}
			int real_ldegree = 0;
			int orderloop = 0;
			do {
				orderloop++;
				real_ldegree = buf.series(r, orderloop, options).ldegree(sym).to_int();
			} while ((real_ldegree == orderloop)
					&& ( factor*real_ldegree < degbound));
			ldegrees[j] = factor * real_ldegree;
			degbound -= factor * real_ldegree;
		}
		j++;
	}

	int degsum = std::accumulate(ldegrees.begin(), ldegrees.end(), 0);

	if (degsum >= order) {
		epvector epv;
		epv.emplace_back(Order(_ex1), order);
		return (new pseries(r, epv))->setflag(status_flags::dynallocated);
	}

	// Multiply with remaining terms
	auto itd = ldegrees.begin();
	for (auto it = seq.begin(); it != seq.end(); ++it, ++itd) {

		// do series expansion with adjusted order
		ex term = recombine_pair_to_ex(*it).series(r, order-degsum+(*itd), options);

		// Series multiplication
		if (it == seq.begin())
			acc = ex_to<pseries>(term);
		else
			acc = ex_to<pseries>(acc.mul_series(ex_to<pseries>(term)));
	}

	return acc.mul_const(overall_coeff);
}


/** Compute the p-th power of a series.
 *
 *  @param p  power to compute
 *  @param deg  truncation order of series calculation */
ex pseries::power_const(const numeric &p, int deg) const
{
	// method:
	// (due to Leonhard Euler)
	// let A(x) be this series and for the time being let it start with a
	// constant (later we'll generalize):
	//     A(x) = a_0 + a_1*x + a_2*x^2 + ...
	// We want to compute
	//     C(x) = A(x)^p
	//     C(x) = c_0 + c_1*x + c_2*x^2 + ...
	// Taking the derivative on both sides and multiplying with A(x) one
	// immediately arrives at
	//     C'(x)*A(x) = p*C(x)*A'(x)
	// Multiplying this out and comparing coefficients we get the recurrence
	// formula
	//     c_i = (i*p*a_i*c_0 + ((i-1)*p-1)*a_{i-1}*c_1 + ...
	//                    ... + (p-(i-1))*a_1*c_{i-1})/(a_0*i)
	// which can easily be solved given the starting value c_0 = (a_0)^p.
	// For the more general case where the leading coefficient of A(x) is not
	// a constant, just consider A2(x) = A(x)*x^m, with some integer m and
	// repeat the above derivation.  The leading power of C2(x) = A2(x)^2 is
	// then of course x^(p*m) but the recurrence formula still holds.
	
	if (seq.empty()) {
		// as a special case, handle the empty (zero) series honoring the
		// usual power laws such as implemented in power::eval()
		if (p.real().is_zero())
			throw std::domain_error("pseries::power_const(): pow(0,I) is undefined");
		else if (p.real().is_negative())
			throw pole_error("pseries::power_const(): division by zero",1);
		else
			return *this;
	}
	
	const int ldeg = ldegree(var).to_int();
	if (!(p*ldeg).is_integer())
		throw std::runtime_error("pseries::power_const(): trying to assemble a Puiseux series");

	// adjust number of coefficients
	int numcoeff = deg - (p*ldeg).to_int();
	if (numcoeff <= 0) {
		epvector epv;
		epv.reserve(1);
		epv.emplace_back(Order(_ex1), deg);
		return (new pseries(relational(var,point), epv))
		       ->setflag(status_flags::dynallocated);
	}
	
	// O(x^n)^(-m) is undefined
	if (seq.size() == 1 && is_order_function(seq[0].rest) && p.real().is_negative())
		throw pole_error("pseries::power_const(): division by zero",1);
	
	// Compute coefficients of the powered series
	exvector co;
	co.reserve(numcoeff);
	co.emplace_back(power(coeff(var, ldeg), p));
	for (int i=1; i<numcoeff; ++i) {
		ex sum = _ex0;
		for (int j=1; j<=i; ++j) {
			ex c = coeff(var, j + ldeg);
			if (is_order_function(c)) {
				co.emplace_back(Order(_ex1));
				break;
			} 
				sum += (p * j - (i - j)) * co[i - j] * c;
		}
		co.push_back(sum / coeff(var, ldeg) / i);
	}
	
	// Construct new series (of non-zero coefficients)
	epvector new_seq;
	bool higher_order = false;
	for (int i=0; i<numcoeff; ++i) {
		if (!co[i].is_zero())
			new_seq.emplace_back(co[i], p * ldeg + i);
		if (is_order_function(co[i])) {
			higher_order = true;
			break;
		}
	}
	if (!higher_order)
		new_seq.emplace_back(Order(_ex1), p * ldeg + numcoeff);

	return pseries(relational(var,point), new_seq);
}


/** Return a new pseries object with the powers shifted by deg. */
pseries pseries::shift_exponents(int deg) const
{
	epvector newseq = seq;
        for (auto & elem : newseq)
		elem.coeff += deg;
	return pseries(relational(var, point), newseq);
}


/** Implementation of ex::series() for powers. This performs Laurent expansion
 *  of reciprocals of series at singularities.
 *  @see ex::series */
ex power::series(const relational & r, int order, unsigned options) const
{
	// If basis is already a series, just power it
	if (is_exactly_a<pseries>(basis)
            and is_exactly_a<numeric>(exponent))
		return ex_to<pseries>(basis).power_const(ex_to<numeric>(exponent), order);

	// Basis is not a series, may there be a singularity?
	bool must_expand_basis = false;
	try {
		ex basis_subs = basis.subs(r, subs_options::no_pattern);
		if (is_exactly_a<infinity>(basis_subs)) {
			must_expand_basis = true;
		}
	} catch (pole_error) {
		must_expand_basis = true;
	}

	bool exponent_is_regular = true;
	try {
		ex exponent_subs = exponent.subs(r, subs_options::no_pattern);
		if (is_exactly_a<infinity>(exponent_subs)) {
			exponent_is_regular = false;
		}
	} catch (pole_error) {
		exponent_is_regular = false;
	}

	if (!exponent_is_regular) {
		ex l = exponent*log(basis);
		// this == exp(l);
		ex le = l.series(r, order, options);
		// Note: expanding exp(l) won't help, since that will attempt
		// Taylor expansion, and fail (because exponent is "singular")
		// Still l itself might be expanded in Taylor series.
		// Examples:
		// sin(x)/x*log(cos(x))
		// 1/x*log(1 + x)
		return exp(le).series(r, order, options);
		// Note: if l happens to have a Laurent expansion (with
		// negative powers of (var - point)), expanding exp(le)
		// will barf (which is The Right Thing).
	}

	// Is the expression of type something^(-int)?
	if (!must_expand_basis && !exponent.info(info_flags::negint)
	 && (!is_exactly_a<add>(basis) || !is_exactly_a<numeric>(exponent)))
		return basic::series(r, order, options);

	// Is the expression of type 0^something?
	if (!must_expand_basis && !basis.subs(r, subs_options::no_pattern).is_zero()
	 && (!is_exactly_a<add>(basis) || !is_exactly_a<numeric>(exponent)))
		return basic::series(r, order, options);

	// Singularity encountered, is the basis equal to (var - point)?
	if (basis.is_equal(r.lhs() - r.rhs())) {
		epvector new_seq;
		if (is_exactly_a<numeric>(exponent)
                    and ex_to<numeric>(exponent).to_int() < order)
			new_seq.emplace_back(_ex1, exponent);
		else
			new_seq.emplace_back(Order(_ex1), exponent);
		return pseries(r, new_seq);
	}

	// No, expand basis into series

	numeric numexp;
	if (is_exactly_a<numeric>(exponent)) {
		numexp = ex_to<numeric>(exponent);
	} else {
		numexp = 0;
	}
	const ex& sym = r.lhs();
	// find existing minimal degree
	ex eb = basis.expand();
	int real_ldegree = 0;
	if (eb.info(info_flags::rational_function))
		real_ldegree = eb.ldegree(sym-r.rhs()).to_int();
	if (real_ldegree == 0) {
		int orderloop = 0;
		do {
			orderloop++;
			real_ldegree = basis.series(r, orderloop, options).ldegree(sym).to_int();
		} while (real_ldegree == orderloop);
	}

	if (!(real_ldegree*numexp).is_integer())
		throw std::runtime_error("pseries::power_const(): trying to assemble a Puiseux series");
	ex e = basis.series(r, (order + real_ldegree*(1-numexp)).to_int(), options);
	
	ex result;
	try {
		result = ex_to<pseries>(e).power_const(numexp, order);
	} catch (pole_error) {
		epvector ser;
		ser.emplace_back(Order(_ex1), order);
		result = pseries(r, ser);
	}

	return result;
}


/** Re-expansion of a pseries object. */
ex pseries::series(const relational & r, int order, unsigned options) const
{
	const ex p = r.rhs();
	GINAC_ASSERT(is_exactly_a<symbol>(r.lhs()));
	const symbol &s = ex_to<symbol>(r.lhs());
	
	if (var.is_equal(s) && point.is_equal(p)) {
		if (order > degree(s))
			return *this;
		
			epvector new_seq;
                        for (const auto & elem : seq) {
				int o = ex_to<numeric>(elem.coeff).to_int();
				if (o >= order) {
					new_seq.emplace_back(Order(_ex1), o);
					break;
				}
				new_seq.push_back(elem);
			}
			return pseries(r, new_seq);
		
	} else
		return convert_to_poly().series(r, order, options);
}


/** Compute the truncated series expansion of an expression.
 *  This function returns an expression containing an object of class pseries 
 *  to represent the series. If the series does not terminate within the given
 *  truncation order, the last term of the series will be an order term.
 *
 *  @param r  expansion relation, lhs holds variable and rhs holds point
 *  @param order  truncation order of series calculations
 *  @param options  of class series_options
 *  @return an expression holding a pseries object */
ex ex::series(const ex & r, int order, unsigned options) const
{
	ex e;
	relational rel_;
	
	if (is_exactly_a<relational>(r))
		rel_ = ex_to<relational>(r);
	else if (is_exactly_a<symbol>(r))
		rel_ = relational(r,_ex0);
	else
		throw (std::logic_error("ex::series(): expansion point has unknown type"));
	
        if ((options & series_options::try_univariate_flint) != 0u
                        and rel_.rhs().is_zero()) {
                options &= ~series_options::try_univariate_flint;
                symbolset syms = rel_.lhs().symbols();
                if (syms.size() == 1
                    and useries_can_handle(*this, *(syms.begin()))) {
                        try {
                                return GiNaC::useries(*this,
                                                *(syms.begin()),
                                                order,
                                                options);
                        }
                        catch(flint_error) {
                                ;
                        }
                }
        }
        e = bp->series(rel_, order, options);
        if ((options & series_options::truncate) != 0u) {
                epvector v = ex_to<pseries>(e).seq;
                if (is_order_function((v.end()-1)->rest)) {
                        v.erase(v.end()-1);
                        return pseries(rel_, v);
                }
        }
        return e;
}

} // namespace GiNaC
