/** @file basic.cpp
 *
 *  Implementation of GiNaC's ABC. */

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

#include "basic.h"
#include "ex.h"
#include "expairseq.h"
#include "ex_utils.h"
#include "numeric.h"
#include "power.h"
#include "add.h"
#include "symbol.h"
#include "lst.h"
#include "relational.h"
#include "operators.h"
#include "wildcard.h"
#include "archive.h"
#include "utils.h"
#include "inifcns.h"

#include <iostream>
#include <stdexcept>
#ifdef DO_GINAC_ASSERT
#  include <typeinfo>
#endif

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(basic, void,
  print_func<print_context>(&basic::do_print).
  print_func<print_tree>(&basic::do_print_tree).
  print_func<print_python_repr>(&basic::do_print_python_repr))

//////////
// default constructor, destructor, copy constructor and assignment operator
//////////

// public

/** basic copy constructor: implicitly assumes that the other class is of
 *  the exact same type (as it's used by duplicate()), so it can copy the
 *  tinfo_key and the hash value. */
basic::basic(const basic & other) : tinfo_key(other.tinfo_key), flags(other.flags & ~status_flags::dynallocated), hashvalue(other.hashvalue)
{
}

/** basic assignment operator: the other object might be of a derived class. */
const basic & basic::operator=(const basic & other)
{
	unsigned fl = other.flags & ~status_flags::dynallocated;
	if (tinfo_key != other.tinfo_key) {
		// The other object is of a derived class, so clear the flags as they
		// might no longer apply (especially hash_calculated). Oh, and don't
		// copy the tinfo_key: it is already set correctly for this object.
		fl &= ~(status_flags::evaluated | status_flags::expanded | status_flags::hash_calculated);
	} else {
		// The objects are of the exact same class, so copy the hash value.
		hashvalue = other.hashvalue;
	}
	flags = fl;
	set_refcount(0);
	return *this;
}

// protected

// none (all inlined)

//////////
// other constructors
//////////

// none (all inlined)

//////////
// archiving
//////////

/** Construct object from archive_node. */
basic::basic(const archive_node &n, lst& /*unused*/) : flags(0)
{
	// Reconstruct tinfo_key from class name
	std::string found_class_name;
	if (n.find_string("class", found_class_name))
		tinfo_key = find_tinfo_key(found_class_name);
	else
		throw (std::runtime_error("archive node contains no class name"));
}

/** Archive the object. */
void basic::archive(archive_node &n) const
{
	n.add_string("class", class_name());
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// public

/** Output to stream. This performs double dispatch on the dynamic type of
 *  *this and the dynamic type of the supplied print context.
 *  @param c print context object that describes the output formatting
 *  @param level value that is used to identify the precedence or indentation
 *               level for placing parentheses and formatting */
void basic::print(const print_context & c, unsigned level) const
{
	print_dispatch(get_class_info(), c, level);
}

/** Like print(), but dispatch to the specified class. Can be used by
 *  implementations of print methods to dispatch to the method of the
 *  superclass.
 *
 *  @see basic::print */
void basic::print_dispatch(const registered_class_info & ri, const print_context & c, unsigned level) const
{
	// Double dispatch on object type and print_context type
	const registered_class_info * r_info = &ri;
	const print_context_class_info * pc_info = &c.get_class_info();

	const std::vector<print_functor> * pdt = &r_info->options.get_print_dispatch_table();
	unsigned id = pc_info->options.get_id();

        while(id >= pdt->size() or not ((*pdt)[id].is_valid())) {

                // Method not found, try parent print_context class
                const print_context_class_info * parent_pc_info = pc_info->get_parent();
                if (parent_pc_info != nullptr) {
                        pc_info = parent_pc_info;
                        id = pc_info->options.get_id();
                        continue;
                }

                // Method still not found, try parent class
                const registered_class_info * parent_reg_info = r_info->get_parent();
                if (parent_reg_info != nullptr) {
                        r_info = parent_reg_info;
                        pc_info = &c.get_class_info();
                        pdt = &r_info->options.get_print_dispatch_table();
                        id = pc_info->options.get_id();
                        continue;
                }

                // Method still not found. This shouldn't happen because basic (the
                // base class of the algebraic hierarchy) registers a method for
                // print_context (the base class of the print context hierarchy),
                // so if we end up here, there's something wrong with the class
                // registry.
                throw (std::runtime_error(std::string("basic::print(): method for ") + class_name() + "/" + c.class_name() + " not found"));

        }

        // Call method
        (*pdt)[id](*this, c, level);
}

/** Default output to stream. */
void basic::do_print(const print_context & c, unsigned /*unused*/) const
{
	c.s << "[" << class_name() << " object]";
}

/** Tree output to stream. */
void basic::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec;
	if (nops() != 0u)
		c.s << ", nops=" << nops();
	c.s << std::endl;
	for (size_t i=0; i<nops(); ++i)
		op(i).print(c, level + c.delta_indent);
}

/** Python parsable output to stream. */
void basic::do_print_python_repr(const print_python_repr & c, unsigned /*unused*/) const
{
	c.s << class_name() << "()";
}

/** Little wrapper around print to be called within a debugger.
 *  This is needed because you cannot call foo.print(cout) from within the
 *  debugger because it might not know what cout is.  This method can be
 *  invoked with no argument and it will simply print to stdout.
 *
 *  @see basic::print
 *  @see basic::dbgprinttree */
void basic::dbgprint() const
{
	this->print(print_dflt(std::cerr));
	std::cerr << std::endl;
}

/** Little wrapper around printtree to be called within a debugger.
 *
 *  @see basic::dbgprint */
void basic::dbgprinttree() const
{
	this->print(print_tree(std::cerr));
}

/** Return relative operator precedence (for parenthezing output). */
unsigned basic::precedence() const
{
	return 70;
}

/** Information about the object.
 *
 *  @see class info_flags */
bool basic::info(unsigned /*unused*/) const
{
	// all possible properties are false for basic objects
	return false;
}

/** Number of operands/members. */
size_t basic::nops() const
{
	// iterating from 0 to nops() on atomic objects should be an empty loop,
	// and accessing their elements is a range error.  Container objects should
	// override this.
	return 0;
}

/** Return operand/member at position i. */
const ex basic::op(size_t /*unused*/) const
{
	throw(std::range_error(std::string("basic::op(): ") + class_name() + std::string(" has no operands")));
}

/** Return modifiable operand/member at position i. */
ex & basic::let_op(size_t /*unused*/)
{
	ensure_if_modifiable();
	throw(std::range_error(std::string("basic::let_op(): ") + class_name() + std::string(" has no operands")));
}

ex basic::operator[](const ex & index) const
{
	if (is_exactly_a<numeric>(index))
		return op(static_cast<size_t>(ex_to<numeric>(index).to_int()));

	throw(std::invalid_argument(std::string("non-numeric indices not supported by ") + class_name()));
}

ex basic::operator[](size_t i) const
{
	return op(i);
}

ex & basic::operator[](const ex & index)
{
	if (is_exactly_a<numeric>(index))
		return let_op(ex_to<numeric>(index).to_int());

	throw(std::invalid_argument(std::string("non-numeric indices not supported by ") + class_name()));
}

ex & basic::operator[](size_t i)
{
	return let_op(i);
}

void basic::set_epseq_from(size_t i, ex e)
{
        dynamic_cast<expairseq&>(*this).set_pair_from(i, e);
}

/** Test for occurrence of a pattern.  An object 'has' a pattern if it matches
 *  the pattern itself or one of the children 'has' it.  As a consequence
 *  (according to the definition of children) given e=x+y+z, e.has(x) is true
 *  but e.has(x+y) is false. */
bool basic::has(const ex & pattern, unsigned options) const
{
	exmap map;
	if (match(pattern, map))
		return true;
	for (size_t i=0; i<nops(); i++)
		if (op(i).has(pattern, options))
			return true;
	
	return false;
}

/** Construct new expression by applying the specified function to all
 *  sub-expressions (one level only, not recursively). */
ex basic::map(map_function & f) const
{
	size_t num = nops();
	if (num == 0)
		return *this;

	basic *copy = nullptr;
	for (size_t i=0; i<num; i++) {
		const ex & o = op(i);
		const ex & n = f(o);
		if (!are_ex_trivially_equal(o, n)) {
			if (copy == nullptr)
				copy = duplicate();
			copy->let_op(i) = n;
		}
	}

	if (copy == nullptr)
	        return *this;

        copy->setflag(status_flags::dynallocated);
        copy->clearflag(status_flags::hash_calculated | status_flags::expanded);
	return *copy;
}

/** Check whether this is a polynomial in the given variables. */
bool basic::is_polynomial(const ex & var) const
{
	return !has(var) || is_equal(ex_to<basic>(var));
}

/** Return degree of highest power in object s. */
numeric basic::degree(const ex & s) const
{
	return is_equal(ex_to<basic>(s)) ? *_num1_p : *_num0_p;
}

/** Return degree of lowest power in object s. */
numeric basic::ldegree(const ex & s) const
{
	return is_equal(ex_to<basic>(s)) ? *_num1_p : *_num0_p;
}

/** Return coefficient of degree n in object s. */
ex basic::coeff(const ex & s, const ex & n) const
{
	if (is_equal(ex_to<basic>(s)))
		return n.is_one() ? _ex1 : _ex0;

        return n.is_zero() ? *this : _ex0;
}

/** Sort expanded expression in terms of powers of some object(s).
 *  @param s object(s) to sort in
 *  @param distributed recursive or distributed form (only used when s is a list) */
ex basic::collect(const ex & s, bool distributed) const
{
	ex x;
	if (is_exactly_a<lst>(s)) {

		// List of objects specified
		if (s.nops() == 0)
			return *this;
		if (s.nops() == 1)
			return collect(s.op(0));
		if (distributed) {
			x = this->expand();
			if (! is_exactly_a<add>(x))
				return x; 
			const lst& l(ex_to<lst>(s));

			exmap cmap;
			cmap[_ex1] = _ex0;
                        for (const auto & xelem : x) {
				ex key = _ex1;
                                ex pre_coeff = xelem;
				for (const auto & lelem : l) {
					ex cexp = pre_coeff.degree(lelem);
					pre_coeff = pre_coeff.coeff(lelem, cexp);
					key *= pow(lelem, cexp);
				}
				auto ci = cmap.find(key);
				if (ci != cmap.end())
					ci->second += pre_coeff;
				else
					cmap.insert(exmap::value_type(key, pre_coeff));
			}

			exvector resv;
                        for (const auto & elem : cmap)
				resv.push_back((elem.first) * (elem.second));
			return (new add(resv))->setflag(status_flags::dynallocated);

		} else {

			// Recursive form
			x = *this;
			size_t n = s.nops() - 1;
			while (true) {
				x = x.collect(s[n]);
				if (n == 0)
					break;
				n--;
			}
		}

	} else {
		// Only one object specified
                expairvec vec;
                ex(*this).coefficients(s, vec);
		for (const auto& term : vec)
			x += term.first * power(s, term.second);
                return x;
	}
	
	// correct for lost fractional arguments and return
	return x + (*this - x).expand();
}

/** Perform automatic non-interruptive term rewriting rules. */
ex basic::eval(int /*unused*/) const
{
	// There is nothing to do for basic objects:
	return hold();
}

/** Function object to be applied by basic::evalf(). */
struct evalf_map_function : public map_function {
	int level;
	PyObject* parent;
	evalf_map_function(int l, PyObject* p) : level(l), parent(p) {}
	ex operator()(const ex & e) override { return evalf(e, level, parent); }
};

/** Evaluate object numerically. */
ex basic::evalf(int level, PyObject* parent) const
{
	if (nops() == 0)
		return *this;
        if (level == 1)
                return *this;
        if (level == -max_recursion_level)
                throw(std::runtime_error("max recursion level reached"));

        evalf_map_function map_evalf(level - 1, parent);
        return map(map_evalf);
}

/** Evaluate sums, products and integer powers of matrices. */
ex basic::evalm() const
{
	return *this;
}

#if 0
/** Perform automatic symbolic evaluations on indexed expression that
 *  contains this object as the base expression. */
ex basic::eval_indexed(const basic & i) const
 // this function can't take a "const ex & i" because that would result
 // in an infinite eval() loop
{
	// There is nothing to do for basic objects
	return i.hold();
}

/** Add two indexed expressions. They are guaranteed to be of class indexed
 *  (or a subclass) and their indices are compatible. This function is used
 *  internally by simplify_indexed().
 *
 *  @param self First indexed expression; its base object is *this
 *  @param other Second indexed expression
 *  @return sum of self and other 
 *  @see ex::simplify_indexed() */
ex basic::add_indexed(const ex & self, const ex & other) const
{
	return self + other;
}

/** Multiply an indexed expression with a scalar. This function is used
 *  internally by simplify_indexed().
 *
 *  @param self Indexed expression; its base object is *this
 *  @param other Numeric value
 *  @return product of self and other
 *  @see ex::simplify_indexed() */
ex basic::scalar_mul_indexed(const ex & self, const numeric & other) const
{
	return self * other;
}

/** Try to contract two indexed expressions that appear in the same product. 
 *  If a contraction exists, the function overwrites one or both of the
 *  expressions and returns true. Otherwise it returns false. It is
 *  guaranteed that both expressions are of class indexed (or a subclass)
 *  and that at least one dummy index has been found. This functions is
 *  used internally by simplify_indexed().
 *
 *  @param self Pointer to first indexed expression; its base object is *this
 *  @param other Pointer to second indexed expression
 *  @param v The complete vector of factors
 *  @return true if the contraction was successful, false otherwise
 *  @see ex::simplify_indexed() */
bool basic::contract_with(exvector::iterator /*unused*/, exvector::iterator /*unused*/, exvector& /*unused*/) const
{
	// Do nothing
	return false;

}
#endif

/** Check whether the expression matches a given pattern. For every wildcard
 *  object in the pattern, an expression of the form "wildcard == matching_expression"
 *  is added to repl_lst. */
bool basic::match(const ex & pattern, exmap& map) const
{
/*
	Sweet sweet shapes, sweet sweet shapes,
	That's the key thing, right right.
	Feed feed face, feed feed shapes,
	But who is the king tonight?
	Who is the king tonight?
	Pattern is the thing, the key thing-a-ling,
	But who is the king of Pattern?
	But who is the king, the king thing-a-ling,
	Who is the king of Pattern?
	Bog is the king, the king thing-a-ling,
	Bog is the king of Pattern.
	Ba bu-bu-bu-bu bu-bu-bu-bu-bu-bu bu-bu
	Bog is the king of Pattern.
*/

	if (is_exactly_a<wildcard>(pattern)) {

		// Wildcard matches anything, but check whether we already have found
		// a match for that wildcard first (if so, the earlier match must be
		// the same expression)
                const auto& it = map.find(pattern);
                if (it != map.end())
		        return is_equal(ex_to<basic>(it->second));
		map[pattern] = *this;
		return true;
	} 

        // Expression must be of the same type as the pattern
        if (tinfo() != ex_to<basic>(pattern).tinfo())
                return false;

        // Number of subexpressions must match
        if (nops() != pattern.nops())
                return false;

        // No subexpressions? Then just compare the objects (there can't be
        // wildcards in the pattern)
        if (nops() == 0)
                return is_equal_same_type(ex_to<basic>(pattern));

        // Check whether attributes that are not subexpressions match
        if (!match_same_type(ex_to<basic>(pattern)))
                return false;

        // Otherwise the subexpressions must match one-to-one
        for (size_t i=0; i<nops(); i++)
                if (!op(i).match(pattern.sorted_op(i), map))
                        return false;

        // Looks similar enough, match found
        return true;
}

/** Helper function for subs(). Does not recurse into subexpressions. */
ex basic::subs_one_level(const exmap & m, unsigned options) const
{
	exmap::const_iterator it;

        if (options & subs_options::no_pattern) {
                ex thisex = *this;
                if (is_exactly_a<numeric>(thisex))
                        return ex_to<numeric>(thisex).subs(m, options);
                for (const auto & pair : m)
                        if (thisex.is_equal(pair.first))
                                return pair.second;
	} else {
                for (const auto & elem : m) {
                        exmap map;
			if (match(ex_to<basic>(elem.first), map)) {
			        lst repl_lst;
                                for (const auto& pair : map)
                                        repl_lst.append(pair.first == pair.second);
                                // avoid infinite recursion when
                                // re-substituting the wildcards
				return elem.second.subs(repl_lst,
                                          options | subs_options::no_pattern);
                        }
		}
	}

	return *this;
}

/** Substitute a set of objects by arbitrary expressions. The ex returned
 *  will already be evaluated. */
ex basic::subs(const exmap & m, unsigned options) const
{
	size_t num = nops();
	if (num != 0u) {

		// Substitute in subexpressions
		for (size_t i=0; i<num; i++) {
			const ex & orig_op = op(i);
			const ex & subsed_op = orig_op.subs(m, options);
			if (!are_ex_trivially_equal(orig_op, subsed_op)) {

				// Something changed, clone the object
				basic *copy = duplicate();
				copy->setflag(status_flags::dynallocated);
				copy->clearflag(status_flags::hash_calculated | status_flags::expanded);

				// Substitute the changed operand
				copy->let_op(i++) = subsed_op;

				// Substitute the other operands
				for (; i<num; i++)
					copy->let_op(i) = op(i).subs(m, options);

				// Perform substitutions on the new object as a whole
				return copy->subs_one_level(m, options);
			}
		}
	}

	// Nothing changed or no subexpressions
	return subs_one_level(m, options);
}

/** Default interface of nth derivative ex::diff(s, n).  It should be called
 *  instead of ::derivative(s) for first derivatives and for nth derivatives it
 *  just recurses down.
 *
 *  @param s symbol to differentiate in
 *  @param nth order of differentiation
 *  @see ex::diff */
ex basic::diff(const symbol & s, unsigned nth) const
{
	// trivial: zeroth derivative
	if (nth==0)
		return ex(*this);
	
	// evaluate unevaluated *this before differentiating
	if (not is_evaluated())
		return ex(*this).diff(s, nth);
	
	ex ndiff = this->derivative(s);
	while (!ndiff.is_zero() &&    // stop differentiating zeros
	       nth>1) {
		ndiff = ndiff.diff(s);
		--nth;
	}
	return ndiff;
}

ex basic::conjugate() const
{
	return *this;
}

ex basic::real_part() const
{
	return real_part_function(*this).hold();
}

ex basic::imag_part() const
{
	return imag_part_function(*this).hold();
}

// protected

/** Function object to be applied by basic::derivative(). */
struct derivative_map_function : public map_function {
	const symbol &s;
	derivative_map_function(const symbol &sym) : s(sym) {}
	ex operator()(const ex & e) override { return diff(e, s); }
};

/** Default implementation of ex::diff(). It maps the operation on the
 *  operands (or returns 0 when the object has no operands).
 *
 *  @see ex::diff */
ex basic::derivative(const symbol & s) const
{
	if (nops() == 0)
		return _ex0;

        derivative_map_function map_derivative(s);
        return map(map_derivative);
}

/** Returns order relation between two objects of same type.  This needs to be
 *  implemented by each class. It may never return anything else than 0,
 *  signalling equality, or +1 and -1 signalling inequality and determining
 *  the canonical ordering.  (Perl hackers will wonder why C++ doesn't feature
 *  the spaceship operator <=> for denoting just this.) */
int basic::compare_same_type(const basic & other) const
{
	return compare_pointers(this, &other);
}

/** Returns true if two objects of same type are equal.  Normally needs
 *  not be reimplemented as long as it wasn't overwritten by some parent
 *  class, since it just calls compare_same_type().  The reason why this
 *  function exists is that sometimes it is easier to determine equality
 *  than an order relation and then it can be overridden. */
bool basic::is_equal_same_type(const basic & other) const
{
	return compare_same_type(other)==0;
}

/** Returns true if the attributes of two objects are similar enough for
 *  a match. This function must not match subexpressions (this is already
 *  done by basic::match()). Only attributes not accessible by op() should
 *  be compared. This is also the reason why this function doesn't take the
 *  wildcard replacement list from match() as an argument: only subexpressions
 *  are subject to wildcard matches. Also, this function only needs to be
 *  implemented for container classes because is_equal_same_type() is
 *  automatically used instead of match_same_type() if nops() == 0.
 *
 *  @see basic::match */
bool basic::match_same_type(const basic & /*unused*/) const
{
	// The default is to only consider subexpressions, but not any other
	// attributes
	return true;
}

unsigned basic::return_type() const
{
	return return_types::commutative;
}

tinfo_t basic::return_type_tinfo() const
{
	return tinfo_key;
}

/** Compute the hash value of an object and if it makes sense to store it in
 *  the objects status_flags, do so.  The method inherited from class basic
 *  computes a hash value based on the type and hash values of possible
 *  members.  For this reason it is well suited for container classes but
 *  atomic classes should override this implementation because otherwise they
 *  would all end up with the same hashvalue. */
long basic::calchash() const
{
	long v = golden_ratio_hash((intptr_t)tinfo());
	for (size_t i=0; i<nops(); i++) {
		v = rotate_left(v);
		v ^= this->op(i).gethash();
	}

	// store calculated hash value only if object is already evaluated
	if (is_evaluated()) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}

	return v;
}

/** Function object to be applied by basic::expand(). */
struct expand_map_function : public map_function {
	unsigned options;
	expand_map_function(unsigned o) : options(o) {}
	ex operator()(const ex & e) override { return e.expand(options); }
};

/** Expand expression, i.e. multiply it out and return the result as a new
 *  expression. */
ex basic::expand(unsigned options) const
{
	if (nops() == 0)
		return (options == 0) ? setflag(status_flags::expanded) : *this;

        expand_map_function map_expand(options);
        return ex_to<basic>(map(map_expand)).setflag(options == 0 ? status_flags::expanded : 0);
	
}


//////////
// non-virtual functions in this class
//////////

// public

/** Compare objects syntactically to establish canonical ordering.
 *  All compare functions return: -1 for *this less than other, 0 equal,
 *  1 greater. */
int basic::compare(const basic & other) const
{
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.total_basic_compares++;
#endif
	const long hash_this = gethash();
	const long hash_other = other.gethash();
	if (hash_this<hash_other) return -1;
	if (hash_this>hash_other) return 1;
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.compare_same_hashvalue++;
#endif

	const tinfo_t& typeid_this = tinfo();
	const tinfo_t& typeid_other = other.tinfo();
	if (typeid_this == typeid_other) {
#ifdef GINAC_COMPARE_STATISTICS
		compare_statistics.compare_same_type++;
#endif
		return compare_same_type(other);
	} 
        return (typeid_this<typeid_other ? -1 : 1);
}

/** Test for syntactic equality.
 *  This is only a quick test, meaning objects should be in the same domain.
 *  You might have to .expand(), .normal() objects first, depending on the
 *  domain of your computation, to get a more reliable answer.
 *
 *  @see is_equal_same_type */
bool basic::is_equal(const basic & other) const
{
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.total_basic_is_equals++;
#endif
	if (this->gethash()!=other.gethash())
		return false;
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.is_equal_same_hashvalue++;
#endif
	if (this->tinfo()!=other.tinfo())
		return false;
	
	GINAC_ASSERT(typeid(*this)==typeid(other));
	
#ifdef GINAC_COMPARE_STATISTICS
	compare_statistics.is_equal_same_type++;
#endif
	return is_equal_same_type(other);
}

// protected

/** Stop further evaluation.
 *
 *  @see basic::eval */
const basic & basic::hold() const
{
	return setflag(status_flags::evaluated);
}

/** Ensure the object may be modified without hurting others, throws if this
 *  is not the case. */
void basic::ensure_if_modifiable() const
{
	if (get_refcount() > 1)
		throw(std::runtime_error("cannot modify multiply referenced object"));
	clearflag(status_flags::hash_calculated | status_flags::evaluated);
}

//////////
// global variables
//////////

int max_recursion_level = 1024;


#ifdef GINAC_COMPARE_STATISTICS
compare_statistics_t::~compare_statistics_t()
{
	std::clog << "ex::compare() called " << total_compares << " times" << std::endl;
	std::clog << "nontrivial compares: " << nontrivial_compares << " times" << std::endl;
	std::clog << "basic::compare() called " << total_basic_compares << " times" << std::endl;
	std::clog << "same hashvalue in compare(): " << compare_same_hashvalue << " times" << std::endl;
	std::clog << "compare_same_type() called " << compare_same_type << " times" << std::endl;
	std::clog << std::endl;
	std::clog << "ex::is_equal() called " << total_is_equals << " times" << std::endl;
	std::clog << "nontrivial is_equals: " << nontrivial_is_equals << " times" << std::endl;
	std::clog << "basic::is_equal() called " << total_basic_is_equals << " times" << std::endl;
	std::clog << "same hashvalue in is_equal(): " << is_equal_same_hashvalue << " times" << std::endl;
	std::clog << "is_equal_same_type() called " << is_equal_same_type << " times" << std::endl;
	std::clog << std::endl;
	std::clog << "basic::gethash() called " << total_gethash << " times" << std::endl;
	std::clog << "used cached hashvalue " << gethash_cached << " times" << std::endl;
}

compare_statistics_t compare_statistics;
#endif

} // namespace GiNaC
